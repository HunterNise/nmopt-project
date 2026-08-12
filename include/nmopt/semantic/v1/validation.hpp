#pragma once

#include "nmopt/semantic/v1/types.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace nmopt::semantic::v1
{
  class ValidationReport
  {
  public:
    bool
    valid() const
    {
      return diagnostics_.empty();
    }

    const std::vector<Diagnostic> &
    diagnostics() const
    {
      return diagnostics_;
    }

    bool
    has_category(const DiagnosticCategory category) const
    {
      return std::any_of(diagnostics_.begin(),
                         diagnostics_.end(),
                         [category](const Diagnostic &diagnostic) {
                           return diagnostic.category == category;
                         });
    }

    void
    add(DiagnosticCategory category,
        std::string        component_id,
        std::string        capability,
        std::string        remedy)
    {
      diagnostics_.push_back(
        {category,
         std::move(component_id),
         std::move(capability),
         std::move(remedy)});
    }

  private:
    std::vector<Diagnostic> diagnostics_;
  };

  // Validation here is deliberately limited to semantic structure and stated
  // policies. Backend lowerability and formulation capabilities belong to the
  // compiler and are appended to the same diagnostic report there.
  class SemanticValidator final
  {
  public:
    ValidationReport
    validate(const ProblemSpec &specification) const
    {
      ValidationReport report;
      if (specification.id.empty())
        report.add(DiagnosticCategory::structural,
                   "problem",
                   "stable_problem_identity",
                   "Give ProblemSpec a non-empty stable identifier.");
      if (specification.label.empty())
        report.add(DiagnosticCategory::structural,
                   specification.id.empty() ? "problem" : specification.id,
                   "human_readable_label",
                   "Give ProblemSpec a non-empty human-readable label.");

      const auto regions = index(specification.regions, report, "region");
      const auto spaces = index(specification.spaces, report, "space");
      const auto pairings = index(specification.pairings, report, "pairing");
      const auto variables = index(specification.variables, report, "variable");
      const auto data = index(specification.data, report, "data");
      const auto transformations = index(specification.transformations,
                                         report,
                                         "transformation");
      const auto terms = index(specification.residual_terms, report,
                               "residual term");
      const auto equations = index(specification.equations, report,
                                   "equation");
      const auto observations = index(specification.observations, report,
                                      "observation");
      index(specification.losses, report, "loss");
      const auto metrics = index(specification.metrics, report, "metric");
      const auto constraints = index(specification.constraints, report,
                                     "constraint");
      index(specification.requirement_policies, report, "requirement policy");

      validate_labels(specification.regions, report, "region");
      validate_labels(specification.spaces, report, "space");
      validate_labels(specification.pairings, report, "pairing");
      validate_labels(specification.variables, report, "variable");
      validate_labels(specification.data, report, "data");
      validate_labels(specification.transformations, report, "transformation");
      validate_labels(specification.residual_terms, report, "residual term");
      validate_labels(specification.equations, report, "equation");
      validate_labels(specification.observations, report, "observation");
      validate_labels(specification.losses, report, "loss");
      validate_labels(specification.metrics, report, "metric");
      validate_labels(specification.constraints, report, "constraint");
      validate_required_enum_fields(specification, report);

      validate_regions(specification, report);
      validate_spaces(specification, regions, report);
      validate_pairings(specification, spaces, report);
      validate_variables(specification, spaces, transformations, report);
      validate_data(specification, spaces, report);
      validate_transformations(specification,
                               variables,
                               data,
                               spaces,
                               report);
      validate_equations(specification, spaces, pairings, terms, report);
      validate_terms(specification,
                     variables,
                     data,
                     equations,
                     spaces,
                     regions,
                     report);
      validate_observations(specification,
                            variables,
                            data,
                            regions,
                            spaces,
                            pairings,
                            report);
      validate_losses(specification, observations, spaces, pairings, data, report);
      validate_metrics(specification, variables, pairings, report);
      validate_constraints(specification, variables, data, spaces, regions, report);
      validate_policy_regions(specification, regions, report);
      validate_formulation(specification.formulation,
                           variables,
                           equations,
                           metrics,
                           constraints,
                           report);
      validate_policies(specification, report);
      return report;
    }

  private:
    template <typename Component>
    using Index = std::unordered_map<std::string, const Component *>;

    template <typename Component>
    static Index<Component>
    index(const std::vector<Component> &components,
          ValidationReport &           report,
          const char *                  component_name)
    {
      Index<Component> result;
      for (const auto &component : components)
        {
          if (component.id.empty())
            {
              report.add(DiagnosticCategory::structural,
                         component_name,
                         "stable_component_identity",
                         "Give every semantic component a non-empty identifier.");
              continue;
            }
          if (!result.emplace(component.id, &component).second)
            report.add(DiagnosticCategory::structural,
                       component.id,
                       "unique_component_identity",
                       "Use a unique identifier for this component kind.");
        }
      return result;
    }

    template <typename Component>
    static bool
    contains(const Index<Component> &index, const std::string &id)
    {
      return index.find(id) != index.end();
    }

    template <typename Component>
    static void
    validate_labels(const std::vector<Component> &components,
                    ValidationReport &           report,
                    const char *                 component_name)
    {
      for (const auto &component : components)
        if (component.label.empty())
          report.add(DiagnosticCategory::structural,
                     component.id.empty() ? component_name : component.id,
                     "human_readable_label",
                     "Give every semantic component a non-empty human-readable label.");
    }

    static void
    require_specified(const bool         specified,
                      const std::string &component_id,
                      const char *       fallback_id,
                      const char *       capability,
                      ValidationReport & report)
    {
      if (!specified)
        report.add(DiagnosticCategory::structural,
                   component_id.empty() ? fallback_id : component_id,
                   capability,
                   "Select an explicit semantic kind, role, status, or scope.");
    }

    static void
    validate_required_enum_fields(const ProblemSpec &specification,
                                  ValidationReport & report)
    {
      for (const auto &region : specification.regions)
        require_specified(region.kind != RegionKind::unspecified,
                          region.id,
                          "region",
                          "region_kind",
                          report);
      for (const auto &space : specification.spaces)
        {
          require_specified(space.topology != SpaceTopology::unspecified,
                            space.id,
                            "space",
                            "space_topology",
                            report);
          require_specified(space.role != SpaceRole::unspecified,
                            space.id,
                            "space",
                            "space_role",
                            report);
        }
      for (const auto &variable : specification.variables)
        require_specified(variable.role != VariableRole::unspecified,
                          variable.id,
                          "variable",
                          "variable_role",
                          report);
      for (const auto &datum : specification.data)
        {
          require_specified(datum.kind != DataKind::unspecified,
                            datum.id,
                            "data",
                            "data_kind",
                            report);
          require_specified(datum.role != DataRole::unspecified,
                            datum.id,
                            "data",
                            "data_role",
                            report);
        }
      for (const auto &transformation : specification.transformations)
        require_specified(
          transformation.kind != TransformationKind::unspecified,
          transformation.id,
          "transformation",
          "transformation_kind",
          report);
      for (const auto &term : specification.residual_terms)
        require_specified(term.kind != ResidualTermKind::unspecified,
                          term.id,
                          "residual term",
                          "residual_term_kind",
                          report);
      for (const auto &observation : specification.observations)
        require_specified(observation.kind != ObservationKind::unspecified,
                          observation.id,
                          "observation",
                          "observation_kind",
                          report);
      for (const auto &loss : specification.losses)
        require_specified(loss.kind != LossKind::unspecified,
                          loss.id,
                          "loss",
                          "loss_kind",
                          report);
      for (const auto &metric : specification.metrics)
        require_specified(metric.kind != MetricKind::unspecified,
                          metric.id,
                          "metric",
                          "metric_kind",
                          report);
      for (const auto &constraint : specification.constraints)
        require_specified(constraint.kind != ConstraintKind::unspecified,
                          constraint.id,
                          "constraint",
                          "constraint_kind",
                          report);
      for (const auto &policy : specification.requirement_policies)
        {
          require_specified(policy.kind != RequirementKind::unspecified,
                            policy.id,
                            "requirement policy",
                            "requirement_kind",
                            report);
          require_specified(policy.status != RequirementStatus::unspecified,
                            policy.id,
                            "requirement policy",
                            "requirement_status",
                            report);
          require_specified(policy.scope != RequirementScope::unspecified,
                            policy.id,
                            "requirement policy",
                            "requirement_scope",
                            report);
        }
      require_specified(
        specification.formulation.kind != FormulationKind::unspecified,
        specification.formulation.id,
        "formulation",
        "formulation_kind",
        report);
    }

    static bool
    pairing_matches_space(const PairingSpec &pairing,
                          const std::string &space_id)
    {
      return pairing.primal_space_id == space_id &&
             pairing.covector_space_id == space_id;
    }

    static void
    validate_regions(const ProblemSpec &specification, ValidationReport &report)
    {
      for (const auto &region : specification.regions)
        {
          if (region.kind != RegionKind::point_set &&
              !region.point_coordinates.empty())
            report.add(DiagnosticCategory::structural,
                       region.id,
                       "non_point_region_has_no_coordinates",
                       "Declare physical coordinates only on a point-set region.");
          if (region.kind == RegionKind::volume && !region.boundary_ids.empty())
            report.add(DiagnosticCategory::structural,
                       region.id,
                       "volume_region_has_no_boundary_ids",
                       "Declare boundary ids only on a boundary region.");
          if (region.kind == RegionKind::volume && region.is_full_domain &&
              !region.material_ids.empty())
            report.add(DiagnosticCategory::structural,
                       region.id,
                       "full_volume_region_has_no_material_ids",
                       "Use a non-full volume region for a material-id restriction.");
          if (region.kind == RegionKind::boundary && region.is_full_domain)
            report.add(DiagnosticCategory::structural,
                       region.id,
                       "boundary_region_is_not_full_domain",
                       "Use a volume region for a full-domain declaration.");
          if (region.kind == RegionKind::boundary && region.boundary_ids.empty())
            report.add(DiagnosticCategory::structural,
                       region.id,
                       "boundary_region_ids",
                       "Declare at least one boundary id for a boundary region.");
          if (region.kind == RegionKind::boundary && !region.material_ids.empty())
            report.add(DiagnosticCategory::structural,
                       region.id,
                       "boundary_region_has_no_material_ids",
                       "Declare material ids only on a volume region.");
          if (region.kind == RegionKind::point_set)
            {
              if (region.is_full_domain || !region.boundary_ids.empty() ||
                  !region.material_ids.empty())
                report.add(DiagnosticCategory::structural,
                           region.id,
                           "point_set_region_geometry",
                           "A point-set region has coordinates only; it is neither a volume nor a boundary region.");
              if (region.point_coordinates.empty())
                report.add(DiagnosticCategory::structural,
                           region.id,
                           "point_set_coordinates",
                           "Declare at least one immutable sensor coordinate.");
              std::size_t coordinate_dimension = 0;
              for (const auto &coordinate : region.point_coordinates)
                {
                  if (coordinate.empty())
                    report.add(DiagnosticCategory::structural,
                               region.id,
                               "point_coordinate_dimension",
                               "Give every sensor coordinate the same positive dimension.");
                  else if (coordinate_dimension == 0)
                    coordinate_dimension = coordinate.size();
                  else if (coordinate.size() != coordinate_dimension)
                    report.add(DiagnosticCategory::structural,
                               region.id,
                               "point_coordinate_dimension",
                               "Give every sensor coordinate the same positive dimension.");
                  if (!std::all_of(coordinate.begin(),
                                   coordinate.end(),
                                   [](const double value) {
                                     return std::isfinite(value);
                                   }))
                    report.add(DiagnosticCategory::structural,
                               region.id,
                               "finite_point_coordinates",
                               "Bind finite physical sensor coordinates.");
                }
              for (std::size_t first = 0;
                   first < region.point_coordinates.size();
                   ++first)
                for (std::size_t second = first + 1;
                     second < region.point_coordinates.size();
                     ++second)
                  if (region.point_coordinates[first] ==
                      region.point_coordinates[second])
                    report.add(DiagnosticCategory::structural,
                               region.id,
                               "unique_point_coordinates",
                               "Declare each sensor coordinate once in the point set.");
            }
        }
    }

    static void
    validate_spaces(const ProblemSpec &      specification,
                    const Index<RegionSpec> &regions,
                    ValidationReport &       report)
    {
      for (const auto &space : specification.spaces)
        {
          if (!contains(regions, space.region_id))
            report.add(DiagnosticCategory::structural,
                       space.id,
                       "space_region_port",
                       "Reference a declared base region.");
          if (!space.is_scalar)
            report.add(DiagnosticCategory::structural,
                       space.id,
                       "scalar_field_shape",
                       "The v1 semantic slice currently declares scalar spaces only.");
        }
    }

    static void
    validate_pairings(const ProblemSpec &     specification,
                      const Index<SpaceSpec> &spaces,
                      ValidationReport &      report)
    {
      for (const auto &pairing : specification.pairings)
        {
          if (!contains(spaces, pairing.primal_space_id) ||
              !contains(spaces, pairing.covector_space_id))
            report.add(DiagnosticCategory::structural,
                       pairing.id,
                       "pairing_space_ports",
                       "Reference declared primal and covector spaces.");
          else if (pairing.primal_space_id != pairing.covector_space_id)
            report.add(
              DiagnosticCategory::structural,
              pairing.id,
              "pairing_primal_covector_space",
              "In the v1 dual-coefficient slice, name the same semantic space on both pairing ports.");
        }
    }

    static void
    validate_variables(const ProblemSpec &                  specification,
                       const Index<SpaceSpec> &             spaces,
                       const Index<TransformationSpec> &    transformations,
                       ValidationReport &                   report)
    {
      for (const auto &variable : specification.variables)
        {
          const auto space = spaces.find(variable.space_id);
          if (space == spaces.end())
            report.add(DiagnosticCategory::structural,
                       variable.id,
                       "variable_space_port",
                       "Reference a declared semantic space.");
          else if ((variable.role == VariableRole::state &&
                    space->second->role != SpaceRole::state) ||
                   (variable.role == VariableRole::control &&
                    space->second->role != SpaceRole::control) ||
                   (variable.role == VariableRole::parameter &&
                    space->second->role != SpaceRole::parameter))
            report.add(DiagnosticCategory::structural,
                       variable.id,
                       "variable_space_role",
                       "Connect each variable to a matching state, control, or parameter space.");
          if (!variable.physical_field_transform_id.empty())
            {
              const auto transformation =
                transformations.find(variable.physical_field_transform_id);
              if (transformation == transformations.end())
                report.add(DiagnosticCategory::structural,
                           variable.id,
                           "physical_field_transformation",
                           "Reference a declared physical-field transformation.");
              else if (transformation->second->input_variable_id != variable.id ||
                       transformation->second->output_space_id != variable.space_id)
                report.add(DiagnosticCategory::structural,
                           variable.id,
                           "physical_field_transformation_ports",
                           "Use a transformation from this variable to its declared space.");
            }
        }
    }

    static void
    validate_data(const ProblemSpec &     specification,
                  const Index<SpaceSpec> &spaces,
                  ValidationReport &      report)
    {
      for (const auto &datum : specification.data)
        {
          if (!datum.space_id.empty() && !contains(spaces, datum.space_id))
            report.add(DiagnosticCategory::structural,
                       datum.id,
                       "data_space_port",
                       "Reference a declared data space or leave scalar constants unspaced.");
          const bool selected_shape =
            (datum.role == DataRole::diffusion &&
             (datum.kind == DataKind::scalar_constant ||
              datum.kind == DataKind::tensor_function)) ||
            (datum.role == DataRole::conservative_transport &&
             datum.kind == DataKind::vector_function) ||
            (datum.role == DataRole::advective_transport &&
             datum.kind == DataKind::vector_function) ||
            (datum.role == DataRole::reaction &&
             (datum.kind == DataKind::scalar_constant ||
              datum.kind == DataKind::function)) ||
            (datum.role == DataRole::robin_coefficient &&
             datum.kind == DataKind::function) ||
            (datum.role == DataRole::robin_source &&
             datum.kind == DataKind::function);
          const bool coefficient_role =
            datum.role == DataRole::diffusion ||
            datum.role == DataRole::conservative_transport ||
            datum.role == DataRole::advective_transport ||
            datum.role == DataRole::reaction ||
            datum.role == DataRole::robin_coefficient ||
            datum.role == DataRole::robin_source;
          if (coefficient_role && !selected_shape)
            report.add(
              DiagnosticCategory::structural,
              datum.id,
              "coefficient_data_shape",
              "Use scalar, vector, or tensor Function data matching the declared coefficient role.");
          if (datum.role == DataRole::observation_weight &&
              datum.kind != DataKind::function)
            report.add(
              DiagnosticCategory::structural,
              datum.id,
              "observation_weight_data_shape",
              "Represent the fixed scalar observation weight as Function data.");
        }
    }

    static void
    validate_transformations(const ProblemSpec &                 specification,
                             const Index<VariableSpec> &         variables,
                             const Index<DataSpec> &             data,
                             const Index<SpaceSpec> &            spaces,
                             ValidationReport &                  report)
    {
      for (const auto &transformation : specification.transformations)
        {
          const auto input = variables.find(transformation.input_variable_id);
          const auto output = spaces.find(transformation.output_space_id);
          const auto fixed_data = data.find(transformation.fixed_data_id);
          const auto control = variables.find(transformation.control_variable_id);
          const bool has_state_ports =
            input != variables.end() &&
            input->second->role == VariableRole::state &&
            output != spaces.end() && output->second->role == SpaceRole::state;
          if (transformation.kind ==
              TransformationKind::fixed_dirichlet_reconstruction)
            {
              if (!has_state_ports || fixed_data == data.end() ||
                  fixed_data->second->kind != DataKind::function ||
                  fixed_data->second->role != DataRole::fixed_dirichlet_lifting ||
                  !transformation.control_variable_id.empty())
                report.add(
                  DiagnosticCategory::structural,
                  transformation.id,
                  "fixed_dirichlet_reconstruction_ports",
                  "Connect the fixed reconstruction to one state variable, its state space, and fixed Function data only.");
              if (fixed_data != data.end() &&
                  fixed_data->second->space_id != transformation.output_space_id)
                report.add(DiagnosticCategory::structural,
                           transformation.id,
                           "fixed_dirichlet_lifting_space",
                           "Declare fixed lifting data in the reconstructed state space.");
            }
          else if (transformation.kind ==
                   TransformationKind::dirichlet_control_lifting)
            {
              if (!has_state_ports || control == variables.end() ||
                  control->second->role != VariableRole::control ||
                  (!transformation.fixed_data_id.empty() &&
                   (fixed_data == data.end() ||
                    fixed_data->second->kind != DataKind::function ||
                    fixed_data->second->role !=
                      DataRole::fixed_dirichlet_lifting ||
                    fixed_data->second->space_id !=
                      transformation.output_space_id)))
                report.add(
                  DiagnosticCategory::structural,
                  transformation.id,
                  "dirichlet_control_lifting_ports",
                  "Connect the Dirichlet lifting to one state variable, one control variable, the reconstructed state space, and optional fixed Function data in that state space.");
            }
          if (input != variables.end() && output != spaces.end() &&
              input->second->space_id != transformation.output_space_id)
            report.add(DiagnosticCategory::structural,
                       transformation.id,
                       "physical_field_reconstruction_space",
                       "Reconstruct the physical field in the state variable's declared space.");
          const auto output_uses = std::count_if(
            specification.variables.begin(),
            specification.variables.end(),
            [&transformation](const VariableSpec &variable) {
              return variable.physical_field_transform_id == transformation.id;
            });
          if (output_uses != 1)
            report.add(DiagnosticCategory::structural,
                       transformation.id,
                       "physical_field_transformation_output",
                       "Connect this transformation to exactly one physical state field.");
        }
    }

    static void
    validate_equations(const ProblemSpec &             specification,
                       const Index<SpaceSpec> &        spaces,
                       const Index<PairingSpec> &      pairings,
                       const Index<ResidualTermSpec> & terms,
                       ValidationReport &              report)
    {
      for (const auto &equation : specification.equations)
        {
          const auto test_space = spaces.find(equation.test_space_id);
          if (test_space == spaces.end())
            report.add(DiagnosticCategory::structural,
                       equation.id,
                       "equation_test_space",
                       "Reference a declared scalar test space.");
          else if (test_space->second->role != SpaceRole::test)
            report.add(DiagnosticCategory::structural,
                       equation.id,
                       "equation_test_space_role",
                       "Reference a space declared with the test role.");
          const auto pairing = pairings.find(equation.test_pairing_id);
          if (pairing == pairings.end() ||
              !pairing_matches_space(*pairing->second,
                                     equation.test_space_id))
            report.add(DiagnosticCategory::structural,
                       equation.id,
                       "equation_test_pairing",
                       "Reference the two-sided pairing for the equation test space.");
          std::unordered_set<std::string> unique_term_ids;
          std::unordered_set<std::string> reported_duplicate_term_ids;
          for (const auto &term_id : equation.residual_term_ids)
            {
              if (!unique_term_ids.insert(term_id).second &&
                  reported_duplicate_term_ids.insert(term_id).second)
                report.add(DiagnosticCategory::structural,
                           equation.id,
                           "unique_equation_residual_term_edges",
                           "List every residual term exactly once in its owning equation.");
              const auto term = terms.find(term_id);
              if (term == terms.end())
                report.add(DiagnosticCategory::structural,
                           equation.id,
                           "equation_residual_term_port",
                           "Reference a declared residual term.");
              else if (term->second->equation_id != equation.id)
                report.add(DiagnosticCategory::structural,
                           term_id,
                           "residual_term_target_equation",
                           "Connect the term to the equation that owns it.");
            }
        }
    }

    static void
    validate_terms(const ProblemSpec &             specification,
                   const Index<VariableSpec> &     variables,
                   const Index<DataSpec> &         data,
                   const Index<EquationBlockSpec> &equations,
                   const Index<SpaceSpec> &        spaces,
                   const Index<RegionSpec> &       regions,
                   ValidationReport &               report)
    {
      for (const auto &term : specification.residual_terms)
        {
          const auto equation = equations.find(term.equation_id);
          if (equation == equations.end())
            report.add(DiagnosticCategory::structural,
                       term.id,
                       "residual_term_target_equation",
                       "Reference a declared target equation block.");
          else if (std::count(equation->second->residual_term_ids.begin(),
                              equation->second->residual_term_ids.end(),
                              term.id) == 0)
            report.add(DiagnosticCategory::structural,
                       term.id,
                       "residual_term_equation_membership",
                       "List this residual term exactly once in its target equation.");
          for (const auto &variable_id : term.variable_ids)
            if (!contains(variables, variable_id))
              report.add(DiagnosticCategory::structural,
                         term.id,
                         "residual_term_variable_port",
                         "Reference a declared variable input.");
          for (const auto &data_id : term.data_ids)
            if (!contains(data, data_id))
              report.add(DiagnosticCategory::structural,
                         term.id,
                         "residual_term_data_port",
                         "Reference declared immutable data.");
          validate_term_signature(term, variables, data, report);
          const bool boundary_term =
            term.kind == ResidualTermKind::neumann_control ||
            term.kind == ResidualTermKind::dirichlet_transposition_control ||
            term.kind == ResidualTermKind::robin_bilinear ||
            term.kind == ResidualTermKind::robin_source;
          if (!boundary_term && !term.region_id.empty())
            report.add(DiagnosticCategory::structural,
                       term.id,
                       "volume_term_has_no_boundary_region",
                       "Leave the region port empty for the registered volume term.");
          if (!boundary_term)
            continue;

          const auto region = regions.find(term.region_id);
          if (region == regions.end() ||
              region->second->kind != RegionKind::boundary)
            report.add(DiagnosticCategory::structural,
                       term.id,
                       term.kind == ResidualTermKind::neumann_control
                         ? "neumann_control_boundary_region"
                       : term.kind ==
                           ResidualTermKind::dirichlet_transposition_control
                         ? "dirichlet_transposition_boundary_region"
                         : "robin_boundary_region",
                       "Declare this natural residual contribution on a boundary region.");
          if (term.kind != ResidualTermKind::neumann_control &&
              term.kind !=
                ResidualTermKind::dirichlet_transposition_control)
            continue;
          const auto control = term.variable_ids.size() == 1
                                 ? variables.find(term.variable_ids.front())
                                 : variables.end();
          const auto control_space = control == variables.end()
                                       ? spaces.end()
                                       : spaces.find(control->second->space_id);
          if (control_space == spaces.end() ||
              control_space->second->region_id != term.region_id)
            report.add(DiagnosticCategory::structural,
                       term.id,
                       term.kind == ResidualTermKind::neumann_control
                         ? "neumann_control_space_region"
                         : "dirichlet_transposition_control_space_region",
                       "Place the boundary control space on the residual term's declared boundary region.");
        }
    }

    static void
    validate_observations(const ProblemSpec &          specification,
                          const Index<VariableSpec> &  variables,
                          const Index<DataSpec> &      data,
                          const Index<RegionSpec> &    regions,
                          const Index<SpaceSpec> &     spaces,
                          const Index<PairingSpec> &   pairings,
                          ValidationReport &            report)
    {
      for (const auto &observation : specification.observations)
        {
          if (!contains(variables, observation.input_variable_id))
            report.add(DiagnosticCategory::structural,
                       observation.id,
                       "observation_input_port",
                       "Reference a declared variable input.");
          const bool boundary_observation =
            observation.kind == ObservationKind::boundary_trace ||
            observation.kind == ObservationKind::weighted_boundary_trace ||
            observation.kind == ObservationKind::boundary_restriction ||
            observation.kind == ObservationKind::normal_flux;
          const bool point_observation =
            observation.kind == ObservationKind::point_sensor;
          const RegionKind expected_region =
            boundary_observation ? RegionKind::boundary
            : point_observation ? RegionKind::point_set
                                : RegionKind::volume;
          if (!contains(regions, observation.region_id))
            report.add(DiagnosticCategory::structural,
                       observation.id,
                       "observation_region_port",
                       "Reference the declared observation region.");
          else if (regions.at(observation.region_id)->kind != expected_region)
            report.add(DiagnosticCategory::structural,
                       observation.id,
                       boundary_observation ? "observation_boundary_region" :
                       point_observation ? "observation_point_set_region" :
                                              "observation_volume_region",
                       boundary_observation
                         ? "The registered boundary observation needs a boundary region."
                       : point_observation
                         ? "The point-sensor observation needs a point-set region."
                         : "The registered volume restriction needs a volume region.");
          const auto input = variables.find(observation.input_variable_id);
          if (observation.kind == ObservationKind::h1_state_restriction)
            {
              const auto input_space =
                input == variables.end() ? spaces.end() :
                                           spaces.find(input->second->space_id);
              if (input == variables.end() ||
                  input->second->role != VariableRole::state)
                report.add(DiagnosticCategory::structural,
                           observation.id,
                           "h1_state_restriction_state_input",
                           "Use an H1 state variable as the source of the energy observation.");
              else if (input_space == spaces.end() ||
                       input_space->second->topology != SpaceTopology::h1)
                report.add(DiagnosticCategory::structural,
                           observation.id,
                           "h1_state_restriction_input_topology",
                           "Declare the observed state in an H1 space.");
            }
          if ((observation.kind == ObservationKind::boundary_trace ||
               observation.kind == ObservationKind::weighted_boundary_trace ||
               observation.kind == ObservationKind::normal_flux) &&
              (input == variables.end() ||
               input->second->role != VariableRole::state))
            report.add(DiagnosticCategory::structural,
                       observation.id,
                       "boundary_trace_state_input",
                       "Use a state variable as the source of a boundary trace.");
          if (observation.kind == ObservationKind::boundary_restriction &&
              (input == variables.end() ||
               input->second->role != VariableRole::control))
            report.add(DiagnosticCategory::structural,
                       observation.id,
                       "boundary_restriction_control_input",
                       "Use the boundary control as the source of its restriction.");
          if (point_observation &&
              (input == variables.end() ||
               input->second->role != VariableRole::state))
            report.add(DiagnosticCategory::structural,
                       observation.id,
                       "point_sensor_state_input",
                       "Use a state variable as the source of a point-sensor observation.");
          const auto output_space = spaces.find(observation.output_space_id);
          if (output_space == spaces.end())
            report.add(DiagnosticCategory::structural,
                       observation.id,
                       "observation_output_space",
                       "Reference a declared observation space.");
          else if (output_space->second->role != SpaceRole::observation)
            report.add(DiagnosticCategory::structural,
                       observation.id,
                       "observation_output_space_role",
                       "Reference a space declared with the observation role.");
          else if (output_space->second->region_id != observation.region_id)
            report.add(DiagnosticCategory::structural,
                       observation.id,
                       "observation_output_region",
                       "Declare the observation output space on the observation region.");
          else if (observation.kind == ObservationKind::h1_state_restriction &&
                   output_space->second->topology != SpaceTopology::h1)
            report.add(DiagnosticCategory::structural,
                       observation.id,
                       "h1_state_restriction_output_topology",
                       "Declare the energy-observation output with H1 topology.");
          if (point_observation && output_space != spaces.end())
            {
              const auto region = regions.find(observation.region_id);
              const std::size_t point_count =
                region == regions.end() ? 0 : region->second->point_coordinates.size();
              if (output_space->second->topology != SpaceTopology::l2 ||
                  output_space->second->dimension != point_count)
                report.add(DiagnosticCategory::structural,
                           observation.id,
                           "point_sensor_output_dimension",
                           "Declare an L2 observation space whose finite dimension equals the point-set cardinality.");
            }
          if (observation.kind == ObservationKind::weighted_boundary_trace)
            {
              if (observation.data_ids.size() != 1)
                report.add(
                  DiagnosticCategory::structural,
                  observation.id,
                  "weighted_boundary_trace_data_port",
                  "Connect the weighted boundary trace to exactly one fixed observation-weight datum.");
              else
                {
                  const auto weight = data.find(observation.data_ids.front());
                  if (weight == data.end())
                    report.add(DiagnosticCategory::structural,
                               observation.id,
                               "observation_data_port",
                               "Reference declared immutable observation data.");
                  else if (weight->second->role !=
                             DataRole::observation_weight ||
                           weight->second->kind != DataKind::function)
                    report.add(
                      DiagnosticCategory::structural,
                      observation.id,
                      "weighted_boundary_trace_data_signature",
                      "Use one scalar Function datum with the observation-weight role.");
                  else if (weight->second->space_id !=
                           observation.output_space_id)
                    report.add(
                      DiagnosticCategory::structural,
                      observation.id,
                      "weighted_boundary_trace_weight_space",
                      "Place the observation weight in the declared boundary observation space.");
                }
            }
          else if (!observation.data_ids.empty())
            report.add(
              DiagnosticCategory::structural,
              observation.id,
              "observation_data_signature",
              "Leave immutable data ports empty for this observation kind.");
          const auto pairing = pairings.find(observation.output_pairing_id);
          if (pairing == pairings.end() ||
              !pairing_matches_space(*pairing->second,
                                     observation.output_space_id))
            report.add(DiagnosticCategory::structural,
                       observation.id,
                       "observation_output_pairing",
                       "Reference the two-sided pairing for the observation output.");
        }
    }

    static void
    validate_losses(const ProblemSpec &             specification,
                    const Index<ObservationSpec> &  observations,
                    const Index<SpaceSpec> &        spaces,
                    const Index<PairingSpec> &      pairings,
                    const Index<DataSpec> &         data,
                    ValidationReport &               report)
    {
      for (const auto &loss : specification.losses)
        {
          if (!contains(observations, loss.source_observation_id))
            report.add(DiagnosticCategory::structural,
                       loss.id,
                       "loss_observation_port",
                       "Connect the loss to a declared observation.");
          if (!contains(data, loss.data_id))
            report.add(DiagnosticCategory::structural,
                       loss.id,
                       "loss_data_port",
                       "Reference declared target or regularisation data.");
          const auto observation = observations.find(loss.source_observation_id);
          const auto pairing = pairings.find(loss.pairing_id);
          if (observation == observations.end() || pairing == pairings.end() ||
              !pairing_matches_space(
                *pairing->second, observation->second->output_space_id))
            report.add(DiagnosticCategory::structural,
                       loss.id,
                       "loss_pairing",
                       "Reference the pairing for the loss observation output.");
          validate_loss_signature(loss, data, report);
          if (loss.kind == LossKind::quadratic_hhalf_control_regularisation &&
              observation != observations.end())
            {
              const auto output_space =
                spaces.find(observation->second->output_space_id);
              if (output_space == spaces.end() ||
                  output_space->second->topology != SpaceTopology::hhalf)
                report.add(
                  DiagnosticCategory::structural,
                  loss.id,
                  "hhalf_control_regularisation_topology",
                  "Declare the fractional control loss on an H1/2 boundary observation space.");
            }
          if (loss.kind == LossKind::quadratic_tracking &&
              observation != observations.end())
            {
              const auto datum = data.find(loss.data_id);
              if (datum != data.end() &&
                  datum->second->space_id !=
                    observation->second->output_space_id)
                report.add(DiagnosticCategory::structural,
                           loss.id,
                           "tracking_target_observation_space",
                           "Bind the tracking target in the selected observation space.");
            }
        }
    }

    static void
    validate_metrics(const ProblemSpec &            specification,
                     const Index<VariableSpec> &    variables,
                     const Index<PairingSpec> &     pairings,
                     ValidationReport &              report)
    {
      for (const auto &metric : specification.metrics)
        {
          const auto variable = variables.find(metric.variable_id);
          if (variable == variables.end())
            report.add(DiagnosticCategory::structural,
                       metric.id,
                       "metric_variable_port",
                       "Reference the primal variable identified by this metric.");
          const auto pairing = pairings.find(metric.pairing_id);
          if (variable == variables.end() || pairing == pairings.end() ||
              !pairing_matches_space(*pairing->second,
                                     variable->second->space_id))
            report.add(DiagnosticCategory::structural,
                       metric.id,
                       "metric_pairing",
                       "Reference the declared primal-dual pairing for the metric variable.");
          else if (variable->second->role != VariableRole::control &&
                   variable->second->role != VariableRole::parameter)
            report.add(DiagnosticCategory::structural,
                       metric.id,
                       "metric_decision_variable",
                       "The v1 metric identifies a control or parameter derivative only.");
        }
    }

    static void
    validate_constraints(const ProblemSpec &          specification,
                         const Index<VariableSpec> &  variables,
                         const Index<DataSpec> &      data,
                         const Index<SpaceSpec> &     spaces,
                         const Index<RegionSpec> &    regions,
                         ValidationReport &            report)
    {
      for (const auto &constraint : specification.constraints)
        {
          const auto variable = variables.find(constraint.variable_id);
          if (variable == variables.end())
            report.add(DiagnosticCategory::structural,
                       constraint.id,
                       "constraint_variable_port",
                       "Reference the constrained control variable.");
          else if (variable->second->role != VariableRole::control &&
                   variable->second->role != VariableRole::parameter)
            report.add(DiagnosticCategory::structural,
                       constraint.id,
                       "constraint_decision_variable",
                       "The registered v1 boxes can constrain a control or parameter variable only.");
          else
            {
              const auto space = spaces.find(variable->second->space_id);
              const auto region = space == spaces.end()
                                    ? regions.end()
                                    : regions.find(space->second->region_id);
              const RegionKind expected_region =
                constraint.kind == ConstraintKind::facewise_box
                  ? RegionKind::boundary
                  : RegionKind::volume;
              if (region == regions.end() || region->second->kind != expected_region)
                report.add(DiagnosticCategory::structural,
                           constraint.id,
                           "constraint_control_region",
                           "Use a cellwise box for a volume control and a facewise box for a boundary control.");
            }
          const auto lower = data.find(constraint.lower_bound_data_id);
          const auto upper = data.find(constraint.upper_bound_data_id);
          if (lower == data.end() || upper == data.end())
            report.add(DiagnosticCategory::structural,
                       constraint.id,
                       "constraint_bound_data_ports",
                       "Reference declared lower and upper bound data.");
          else if (lower->second->role != DataRole::lower_bound ||
                   upper->second->role != DataRole::upper_bound ||
                   lower->second->kind !=
                     (constraint.kind == ConstraintKind::facewise_box
                        ? DataKind::facewise_bound
                        : DataKind::cellwise_bound) ||
                   upper->second->kind !=
                     (constraint.kind == ConstraintKind::facewise_box
                        ? DataKind::facewise_bound
                        : DataKind::cellwise_bound))
            report.add(DiagnosticCategory::structural,
                       constraint.id,
                       constraint.kind == ConstraintKind::facewise_box
                         ? "constraint_facewise_bound_data"
                         : "constraint_cellwise_bound_data",
                       "Use the declared lower and upper bound data for this control layout.");
          if (variable != variables.end() && lower != data.end() &&
              upper != data.end() &&
              (lower->second->space_id != variable->second->space_id ||
               upper->second->space_id != variable->second->space_id))
            report.add(DiagnosticCategory::structural,
                       constraint.id,
                       "constraint_bound_data_space",
                       "Declare both bounds in the constrained variable's semantic space.");
        }
    }

    static void
    validate_policy_regions(const ProblemSpec &      specification,
                            const Index<RegionSpec> &regions,
                            ValidationReport &       report)
    {
      for (const auto &policy : specification.requirement_policies)
        if (!policy.region_id.empty() && !contains(regions, policy.region_id))
          report.add(DiagnosticCategory::structural,
                     policy.id,
                     "requirement_policy_region",
                     "Reference a declared region when a policy is region-specific.");
    }

    static bool
    has_role(const std::vector<std::string> &ids,
             const Index<DataSpec> &         data,
             const DataRole                  role)
    {
      return std::any_of(ids.begin(), ids.end(), [&data, role](const auto &id) {
        const auto entry = data.find(id);
        return entry != data.end() && entry->second->role == role;
      });
    }

    static bool
    has_variable_role(const std::vector<std::string> &ids,
                      const Index<VariableSpec> &     variables,
                      const VariableRole               role)
    {
      return std::any_of(ids.begin(), ids.end(), [&variables, role](const auto &id) {
        const auto entry = variables.find(id);
        return entry != variables.end() && entry->second->role == role;
      });
    }

    static void
    validate_term_signature(const ResidualTermSpec &    term,
                            const Index<VariableSpec> & variables,
                            const Index<DataSpec> &     data,
                            ValidationReport &          report)
    {
      bool valid = false;
      switch (term.kind)
        {
          case ResidualTermKind::unspecified:
            break;
          case ResidualTermKind::diffusion_reaction:
            valid = term.variable_ids.size() == 1 && term.data_ids.size() == 2 &&
                    has_variable_role(term.variable_ids, variables,
                                      VariableRole::state) &&
                    has_role(term.data_ids, data, DataRole::diffusion) &&
                    has_role(term.data_ids, data, DataRole::reaction);
            break;
          case ResidualTermKind::tensor_diffusion:
            valid = term.variable_ids.size() == 1 && term.data_ids.size() == 1 &&
                    has_variable_role(term.variable_ids, variables,
                                      VariableRole::state) &&
                    has_role(term.data_ids, data, DataRole::diffusion) &&
                    term.region_id.empty();
            break;
          case ResidualTermKind::conservative_transport:
            valid = term.variable_ids.size() == 1 && term.data_ids.size() == 1 &&
                    has_variable_role(term.variable_ids, variables,
                                      VariableRole::state) &&
                    has_role(term.data_ids,
                             data,
                             DataRole::conservative_transport) &&
                    term.region_id.empty();
            break;
          case ResidualTermKind::advective_transport:
            valid = term.variable_ids.size() == 1 && term.data_ids.size() == 1 &&
                    has_variable_role(term.variable_ids, variables,
                                      VariableRole::state) &&
                    has_role(term.data_ids,
                             data,
                             DataRole::advective_transport) &&
                    term.region_id.empty();
            break;
          case ResidualTermKind::reaction:
            valid = term.variable_ids.size() == 1 && term.data_ids.size() == 1 &&
                    has_variable_role(term.variable_ids, variables,
                                      VariableRole::state) &&
                    has_role(term.data_ids, data, DataRole::reaction) &&
                    term.region_id.empty();
            break;
          case ResidualTermKind::parameter_diffusion_reaction:
            valid = term.variable_ids.size() == 2 && term.data_ids.size() == 1 &&
                    has_variable_role(term.variable_ids, variables,
                                      VariableRole::state) &&
                    has_variable_role(term.variable_ids, variables,
                                      VariableRole::parameter) &&
                    has_role(term.data_ids, data, DataRole::reaction);
            break;
          case ResidualTermKind::laplacian:
            valid = term.variable_ids.size() == 1 && term.data_ids.empty() &&
                    has_variable_role(term.variable_ids, variables,
                                      VariableRole::state) &&
                    term.region_id.empty();
            break;
          case ResidualTermKind::transposition_laplacian:
            valid = term.variable_ids.size() == 1 && term.data_ids.empty() &&
                    has_variable_role(term.variable_ids, variables,
                                      VariableRole::state) &&
                    term.region_id.empty();
            break;
          case ResidualTermKind::dirichlet_transposition_control:
            valid = term.variable_ids.size() == 1 && term.data_ids.empty() &&
                    has_variable_role(term.variable_ids, variables,
                                      VariableRole::control) &&
                    !term.region_id.empty();
            break;
          case ResidualTermKind::volume_source:
            valid = term.variable_ids.empty() && term.data_ids.size() == 1 &&
                    has_role(term.data_ids, data, DataRole::forcing);
            break;
          case ResidualTermKind::volume_control:
            valid = term.variable_ids.size() == 1 &&
                    has_variable_role(term.variable_ids, variables,
                                      VariableRole::control) &&
                    term.data_ids.empty() && term.region_id.empty();
            break;
          case ResidualTermKind::neumann_control:
            valid = term.variable_ids.size() == 1 &&
                    has_variable_role(term.variable_ids, variables,
                                      VariableRole::control) &&
                    term.data_ids.empty() && !term.region_id.empty();
            break;
          case ResidualTermKind::robin_bilinear:
            valid = term.variable_ids.size() == 1 && term.data_ids.size() == 1 &&
                    has_variable_role(term.variable_ids, variables,
                                      VariableRole::state) &&
                    has_role(term.data_ids,
                             data,
                             DataRole::robin_coefficient) &&
                    !term.region_id.empty();
            break;
          case ResidualTermKind::robin_source:
            valid = term.variable_ids.empty() && term.data_ids.size() == 1 &&
                    has_role(term.data_ids, data, DataRole::robin_source) &&
                    !term.region_id.empty();
            break;
        }
      if (!valid)
        report.add(DiagnosticCategory::structural,
                   term.id,
                   "residual_term_signature",
                   "Supply the declared v1 term inputs and no undeclared ports.");
    }

    static void
    validate_loss_signature(const LossSpec &        loss,
                            const Index<DataSpec> & data,
                            ValidationReport &      report)
    {
      const auto datum = data.find(loss.data_id);
      if (datum == data.end())
        return;

      const bool valid =
        (loss.kind == LossKind::quadratic_tracking &&
         datum->second->role == DataRole::desired_state) ||
        (loss.kind == LossKind::quadratic_control_regularisation &&
         datum->second->role == DataRole::regularisation_weight) ||
        (loss.kind == LossKind::quadratic_hhalf_control_regularisation &&
         datum->second->role == DataRole::regularisation_weight) ||
        (loss.kind == LossKind::quadratic_h1_control_regularisation &&
         datum->second->role == DataRole::regularisation_weight) ||
        (loss.kind == LossKind::quadratic_parameter_regularisation &&
         datum->second->role == DataRole::regularisation_weight);
      if (!valid)
        report.add(DiagnosticCategory::structural,
                   loss.id,
                   "loss_data_role",
                   "Use desired-state data for tracking and a scalar regularisation weight for a decision-variable loss.");
    }

    static void
    validate_formulation(const ReducedFormulationSpec & formulation,
                         const Index<VariableSpec> &    variables,
                         const Index<EquationBlockSpec> & equations,
                         const Index<MetricSpec> &       metrics,
                         const Index<ConstraintSpec> &   constraints,
                         ValidationReport &              report)
    {
      const auto state = variables.find(formulation.state_variable_id);
      const auto control = variables.find(formulation.control_variable_id);
      if (formulation.id.empty())
        report.add(DiagnosticCategory::structural,
                   "formulation",
                   "stable_component_identity",
                   "Give the reduced formulation a non-empty identifier.");
      if (state == variables.end() || state->second->role != VariableRole::state)
        report.add(DiagnosticCategory::structural,
                   formulation.id,
                   "formulation_state_variable",
                   "Select one declared state variable for elimination.");
      if (control == variables.end() ||
          (control->second->role != VariableRole::control &&
           control->second->role != VariableRole::parameter))
        report.add(DiagnosticCategory::structural,
                   formulation.id,
                   "formulation_decision_variable",
                   "Select one declared control or parameter variable for optimisation.");
      if (!contains(equations, formulation.equation_id))
        report.add(DiagnosticCategory::structural,
                   formulation.id,
                   "formulation_equation",
                   "Reference the equation that defines the state.");
      const auto metric = metrics.find(formulation.metric_id);
      if (metric == metrics.end())
        report.add(DiagnosticCategory::structural,
                   formulation.id,
                   "formulation_metric",
                   "Reference a declared search metric.");
      else if (control != variables.end() &&
               metric->second->variable_id != formulation.control_variable_id)
        report.add(DiagnosticCategory::structural,
                   formulation.id,
                   "formulation_metric_variable",
                   "Select a metric acting on the formulation decision variable.");
      if (!formulation.constraint_id.empty())
        {
          const auto constraint = constraints.find(formulation.constraint_id);
          if (constraint == constraints.end())
            report.add(DiagnosticCategory::structural,
                       formulation.id,
                       "formulation_constraint",
                       "Reference a declared constraint or leave this port empty.");
          else if (control != variables.end() &&
                   constraint->second->variable_id !=
                     formulation.control_variable_id)
            report.add(
              DiagnosticCategory::structural,
              formulation.id,
              "formulation_constraint_variable",
              "Select a constraint acting on the formulation decision variable.");
        }
    }

    static void
    validate_policies(const ProblemSpec &specification,
                      ValidationReport & report)
    {
      const auto selected_policy = [&specification](
                                     const std::string &subject,
                                     const RequirementKind kind) {
        return std::find_if(
          specification.requirement_policies.begin(),
          specification.requirement_policies.end(),
          [&subject, kind](const RequirementPolicySpec &policy) {
            return policy.subject_id == subject && policy.kind == kind &&
                   policy.status ==
                     RequirementStatus::selected_discrete_realisation &&
                   !policy.selected_policy.empty();
          });
      };
      const auto has_policy = [&specification, &selected_policy](
                                const std::string &  subject,
                                const RequirementKind kind) {
        return selected_policy(subject, kind) !=
               specification.requirement_policies.end();
      };

      for (const auto &variable : specification.variables)
        if (variable.role == VariableRole::state)
          {
            const bool has_fixed_dirichlet =
              has_policy(variable.id, RequirementKind::fixed_dirichlet);
            const bool has_controlled_dirichlet =
              has_policy(variable.id, RequirementKind::controlled_dirichlet);
            const bool has_mean_zero_multiplier =
              has_policy(variable.id, RequirementKind::mean_zero_multiplier);
            const unsigned int uniqueness_policies =
              static_cast<unsigned int>(has_fixed_dirichlet) +
              static_cast<unsigned int>(has_controlled_dirichlet) +
              static_cast<unsigned int>(has_mean_zero_multiplier);
            if (uniqueness_policies == 0)
              report.add(
                DiagnosticCategory::analytical_policy,
                variable.id,
                "state_uniqueness_realisation",
                "Declare one selected fixed-Dirichlet, controlled-Dirichlet, or mean-zero multiplier policy.");
            if (uniqueness_policies > 1)
              report.add(
                DiagnosticCategory::structural,
                variable.id,
                "state_uniqueness_policy_conflict",
                "Select exactly one state uniqueness policy for the declared residual.");
          }

      const bool has_transposition_state_action = std::any_of(
        specification.residual_terms.begin(),
        specification.residual_terms.end(),
        [](const ResidualTermSpec &term) {
          return term.kind == ResidualTermKind::transposition_laplacian;
        });
      const bool has_transposition_boundary_action = std::any_of(
        specification.residual_terms.begin(),
        specification.residual_terms.end(),
        [](const ResidualTermSpec &term) {
          return term.kind ==
                 ResidualTermKind::dirichlet_transposition_control;
        });
      if (has_transposition_state_action || has_transposition_boundary_action)
        {
          const auto state = std::find_if(
            specification.variables.begin(),
            specification.variables.end(),
            [&specification](const VariableSpec &candidate) {
              return candidate.id == specification.formulation.state_variable_id;
            });
          const auto control = std::find_if(
            specification.variables.begin(),
            specification.variables.end(),
            [&specification](const VariableSpec &candidate) {
              return candidate.id == specification.formulation.control_variable_id;
            });
          const auto equation = std::find_if(
            specification.equations.begin(),
            specification.equations.end(),
            [&specification](const EquationBlockSpec &candidate) {
              return candidate.id == specification.formulation.equation_id;
            });
          const auto state_space = state == specification.variables.end()
                                     ? specification.spaces.end()
                                     : std::find_if(
                                         specification.spaces.begin(),
                                         specification.spaces.end(),
                                         [&state](const SpaceSpec &candidate) {
                                           return candidate.id == state->space_id;
                                         });
          const auto control_space = control == specification.variables.end()
                                       ? specification.spaces.end()
                                       : std::find_if(
                                           specification.spaces.begin(),
                                           specification.spaces.end(),
                                           [&control](const SpaceSpec &candidate) {
                                             return candidate.id == control->space_id;
                                           });
          const auto test_space = equation == specification.equations.end()
                                    ? specification.spaces.end()
                                    : std::find_if(
                                        specification.spaces.begin(),
                                        specification.spaces.end(),
                                        [&equation](const SpaceSpec &candidate) {
                                          return candidate.id ==
                                                 equation->test_space_id;
                                        });
          if (!has_transposition_state_action ||
              !has_transposition_boundary_action ||
              state_space == specification.spaces.end() ||
              state_space->topology != SpaceTopology::l2 ||
              control_space == specification.spaces.end() ||
              control_space->topology != SpaceTopology::l2 ||
              test_space == specification.spaces.end() ||
              test_space->topology != SpaceTopology::h2)
            report.add(
              DiagnosticCategory::structural,
              specification.formulation.equation_id,
              "transposition_space_topologies",
              "Declare one L2 state action and one L2 Dirichlet-control action tested in H2 cap H1_0.");
          if (state != specification.variables.end() &&
              !state->physical_field_transform_id.empty())
            report.add(
              DiagnosticCategory::structural,
              state->id,
              "transposition_physical_state",
              "Declare the L2 very-weak state directly; the conforming lifting is a selected discrete lowerer, not the continuous residual.");

          const auto has_exact_policy = [&specification](
                                          const std::string &subject,
                                          const RequirementKind kind,
                                          const RequirementStatus status,
                                          const RequirementScope scope) {
            return std::any_of(
              specification.requirement_policies.begin(),
              specification.requirement_policies.end(),
              [&subject, kind, status, scope](
                const RequirementPolicySpec &policy) {
                return policy.subject_id == subject && policy.kind == kind &&
                       policy.status == status && policy.scope == scope &&
                       !policy.selected_policy.empty();
              });
          };
          if (!has_exact_policy(specification.formulation.equation_id,
                                RequirementKind::transposition_formulation,
                                RequirementStatus::provided,
                                RequirementScope::continuous_semantics))
            report.add(
              DiagnosticCategory::analytical_policy,
              specification.formulation.equation_id,
              "transposition_formulation_policy",
              "Declare the provided L2-state/H2-test transposition residual and its boundary pairing.");
          if (!has_exact_policy(specification.formulation.equation_id,
                                RequirementKind::domain_regularity,
                                RequirementStatus::user_assumed,
                                RequirementScope::continuous_semantics))
            report.add(
              DiagnosticCategory::analytical_policy,
              specification.formulation.equation_id,
              "transposition_domain_regularity",
              "Declare the model author's domain regularity assumption for the Dirichlet Laplacian isomorphism.");
          if (control == specification.variables.end() ||
              !has_exact_policy(control->id,
                                RequirementKind::conforming_trace_subspace,
                                RequirementStatus::selected_discrete_realisation,
                                RequirementScope::discrete_compilation))
            report.add(
              DiagnosticCategory::analytical_policy,
              specification.formulation.control_variable_id,
              "transposition_conforming_trace_subspace",
              "Select the conforming trace-control subspace before using the variational lifting lowerer.");
          if (!has_exact_policy(specification.formulation.equation_id,
                                RequirementKind::conormal_flux,
                                RequirementStatus::selected_discrete_realisation,
                                RequirementScope::discrete_compilation))
            report.add(
              DiagnosticCategory::analytical_policy,
              specification.formulation.equation_id,
              "transposition_conormal_policy",
              "Select the outward discrete conormal as the lifting pullback of the adjoint residual.");
        }

      for (const auto &metric : specification.metrics)
        if (metric.kind == MetricKind::hhalf)
          {
            const auto variable = std::find_if(
              specification.variables.begin(),
              specification.variables.end(),
              [&metric](const VariableSpec &candidate) {
                return candidate.id == metric.variable_id;
              });
            const auto space = std::find_if(
              specification.spaces.begin(),
              specification.spaces.end(),
              [&variable, &specification](const SpaceSpec &candidate) {
                return variable != specification.variables.end() &&
                       candidate.id == variable->space_id;
              });
            if (space == specification.spaces.end() ||
                space->topology != SpaceTopology::hhalf)
              report.add(
                DiagnosticCategory::structural,
                metric.id,
                "hhalf_metric_search_space",
                "Declare the selected fractional metric on an H1/2 boundary control space.");
            if (selected_policy(metric.id,
                                RequirementKind::fractional_trace_realisation) ==
                specification.requirement_policies.end())
              report.add(
                DiagnosticCategory::analytical_policy,
                metric.id,
                "hhalf_metric_realisation_policy",
                "Select a discrete spectral, extension, or auxiliary realization for the H1/2 Riesz map.");
          }

      for (const auto &metric : specification.metrics)
        if (metric.kind == MetricKind::h1)
          {
            const auto variable = std::find_if(
              specification.variables.begin(),
              specification.variables.end(),
              [&metric](const VariableSpec &candidate) {
                return candidate.id == metric.variable_id;
              });
            const auto space = std::find_if(
              specification.spaces.begin(),
              specification.spaces.end(),
              [&variable, &specification](const SpaceSpec &candidate) {
                return variable != specification.variables.end() &&
                       candidate.id == variable->space_id;
              });
            const auto region =
              space == specification.spaces.end()
                ? specification.regions.end()
                : std::find_if(specification.regions.begin(),
                               specification.regions.end(),
                               [&space](const RegionSpec &candidate) {
                                 return candidate.id == space->region_id;
                               });
            if (region != specification.regions.end() &&
                region->kind == RegionKind::boundary &&
                selected_policy(
                  metric.id,
                  RequirementKind::tangential_gradient_realisation) ==
                  specification.requirement_policies.end())
              report.add(
                DiagnosticCategory::analytical_policy,
                metric.id,
                "boundary_h1_metric_tangential_policy",
                "Select the boundary mass-plus-tangential-stiffness realization for the H1 trace Riesz map.");
          }

      for (const auto &metric : specification.metrics)
        if (metric.kind == MetricKind::hminus1)
          {
            const auto variable = std::find_if(
              specification.variables.begin(),
              specification.variables.end(),
              [&metric](const VariableSpec &candidate) {
                return candidate.id == metric.variable_id;
              });
            const auto space = std::find_if(
              specification.spaces.begin(),
              specification.spaces.end(),
              [&specification, variable](const SpaceSpec &candidate) {
                return variable != specification.variables.end() &&
                       candidate.id == variable->space_id;
              });
            if (space == specification.spaces.end() ||
                space->topology != SpaceTopology::h1)
              report.add(
                DiagnosticCategory::structural,
                metric.id,
                "hminus1_metric_search_space",
                "Declare the selected H-1 metric on its continuous control search space.");

            const auto boundary_policy = selected_policy(
              metric.variable_id, RequirementKind::fixed_dirichlet);
            if (boundary_policy == specification.requirement_policies.end())
              report.add(
                DiagnosticCategory::analytical_policy,
                metric.id,
                "hminus1_metric_boundary_policy",
                "Declare the selected Dirichlet or mean policy for the H-1 Riesz operator.");
            else
              {
                const auto region = std::find_if(
                  specification.regions.begin(),
                  specification.regions.end(),
                  [&boundary_policy](const RegionSpec &candidate) {
                    return candidate.id == boundary_policy->region_id;
                  });
                if (region == specification.regions.end() ||
                    region->kind != RegionKind::boundary ||
                    region->boundary_ids.empty())
                  report.add(
                    DiagnosticCategory::structural,
                    metric.id,
                    "hminus1_metric_dirichlet_region",
                    "Select a non-empty boundary region for the Dirichlet H-1 Riesz operator.");
              }
          }

      for (const auto &datum : specification.data)
        if ((datum.role == DataRole::desired_state ||
             datum.role == DataRole::observation_weight) &&
            !has_policy(datum.id,
                        RequirementKind::analytic_quadrature_evaluation))
          report.add(DiagnosticCategory::analytical_policy,
                     datum.id,
                     datum.role == DataRole::desired_state
                       ? "desired_state_data_rule"
                       : "observation_weight_data_rule",
                     datum.role == DataRole::desired_state
                       ? "Declare the analytic quadrature target-data realization."
                       : "Declare the analytic quadrature observation-weight realization.");

      for (const auto &datum : specification.data)
        if (datum.role == DataRole::observation_weight)
          {
            const auto boundedness = std::find_if(
              specification.requirement_policies.begin(),
              specification.requirement_policies.end(),
              [&datum](const RequirementPolicySpec &policy) {
                return policy.subject_id == datum.id &&
                       policy.kind ==
                         RequirementKind::coefficient_regularity &&
                       (policy.status == RequirementStatus::provided ||
                        policy.status == RequirementStatus::user_assumed) &&
                       !policy.selected_policy.empty();
              });
            if (boundedness == specification.requirement_policies.end())
              report.add(
                DiagnosticCategory::analytical_policy,
                datum.id,
                "observation_weight_boundedness",
                "Declare the model author's L-infinity assumption for the boundary observation weight.");
          }

      for (const auto &loss : specification.losses)
        if (loss.kind == LossKind::quadratic_tracking)
          {
            const auto observation = std::find_if(
              specification.observations.begin(),
              specification.observations.end(),
              [&loss](const ObservationSpec &candidate) {
                return candidate.id == loss.source_observation_id;
              });
            const auto datum = std::find_if(
              specification.data.begin(), specification.data.end(),
              [&loss](const DataSpec &candidate) {
                return candidate.id == loss.data_id;
              });
            if (observation == specification.observations.end() ||
                datum == specification.data.end() ||
                datum->role != DataRole::desired_state)
              continue;
            const auto policy = selected_policy(
              datum->id, RequirementKind::analytic_quadrature_evaluation);
            if (policy != specification.requirement_policies.end() &&
                policy->region_id != observation->region_id)
              report.add(DiagnosticCategory::structural,
                         loss.id,
                         "tracking_target_data_region",
                         "Declare the selected target-data realization on the tracking observation region.");
          }

      for (const auto &observation : specification.observations)
        if (observation.kind == ObservationKind::point_sensor)
          {
            if (!has_policy(observation.id,
                            RequirementKind::analytic_quadrature_evaluation))
              report.add(
                DiagnosticCategory::analytical_policy,
                observation.id,
                "point_sensor_evaluation_policy",
                "Declare the physical point-evaluation and finite-dimensional transpose rule for the sensor map.");
            if (std::find_if(
                  specification.requirement_policies.begin(),
                  specification.requirement_policies.end(),
                  [&specification](const RequirementPolicySpec &policy) {
                    return policy.subject_id == specification.formulation.equation_id &&
                           policy.kind == RequirementKind::transposition_formulation &&
                           policy.status == RequirementStatus::provided &&
                           policy.scope == RequirementScope::continuous_semantics &&
                           !policy.selected_policy.empty();
                  }) == specification.requirement_policies.end())
              report.add(
                DiagnosticCategory::analytical_policy,
                specification.formulation.equation_id,
                "point_sensor_transposition_policy",
                "Declare the strong state space, transposition map, residual codomain, and very-weak adjoint policy for point sensors.");
            if (std::find_if(
                  specification.requirement_policies.begin(),
                  specification.requirement_policies.end(),
                  [&specification](const RequirementPolicySpec &policy) {
                    return policy.subject_id == specification.formulation.equation_id &&
                           policy.kind == RequirementKind::domain_regularity &&
                           policy.status == RequirementStatus::user_assumed &&
                           !policy.selected_policy.empty();
                  }) == specification.requirement_policies.end())
              report.add(
                DiagnosticCategory::analytical_policy,
                specification.formulation.equation_id,
                "point_sensor_domain_regularity",
                "Declare the model author's domain regularity assumption for point evaluation and the very-weak adjoint.");
          }

      for (const auto &observation : specification.observations)
        if (observation.kind == ObservationKind::weighted_boundary_trace &&
            observation.data_ids.size() == 1)
          {
            const auto policy = selected_policy(
              observation.data_ids.front(),
              RequirementKind::analytic_quadrature_evaluation);
            if (policy != specification.requirement_policies.end() &&
                policy->region_id != observation.region_id)
              report.add(
                DiagnosticCategory::structural,
                observation.id,
                "weighted_boundary_trace_data_region",
                "Evaluate the fixed observation weight on the weighted trace region.");
          }

      for (const auto &constraint : specification.constraints)
        if (!has_policy(constraint.id,
                        constraint.kind == ConstraintKind::facewise_box
                          ? RequirementKind::discrete_facewise_bounds
                          : RequirementKind::discrete_cellwise_bounds))
          report.add(DiagnosticCategory::analytical_policy,
                     constraint.id,
                     constraint.kind == ConstraintKind::facewise_box
                       ? "facewise_bound_realisation"
                       : "cellwise_bound_realisation",
                     constraint.kind == ConstraintKind::facewise_box
                       ? "Declare the facewise-constant coefficientwise bound policy."
                       : "Declare the FE_DGQ(0) coefficientwise bound policy.");

      const auto tensor_diffusion = std::find_if(
        specification.residual_terms.begin(),
        specification.residual_terms.end(),
        [](const ResidualTermSpec &term) {
          return term.kind == ResidualTermKind::tensor_diffusion;
        });
      if (tensor_diffusion != specification.residual_terms.end())
        {
          const auto declared_assumption = [&specification](
                                             const std::string &subject,
                                             const RequirementKind kind) {
            return std::find_if(
              specification.requirement_policies.begin(),
              specification.requirement_policies.end(),
              [&subject, kind](const RequirementPolicySpec &policy) {
                return policy.subject_id == subject && policy.kind == kind &&
                       (policy.status == RequirementStatus::provided ||
                        policy.status == RequirementStatus::user_assumed) &&
                       !policy.selected_policy.empty();
              });
          };
          const std::string diffusion_data_id =
            tensor_diffusion->data_ids.empty() ? std::string{} :
                                                 tensor_diffusion->data_ids.front();
          if (declared_assumption(diffusion_data_id,
                                  RequirementKind::uniform_ellipticity) ==
              specification.requirement_policies.end())
            report.add(DiagnosticCategory::analytical_policy,
                       diffusion_data_id,
                       "uniform_ellipticity_assumption",
                       "Declare the model author's uniform-ellipticity assumption for the tensor diffusion data.");
          if (declared_assumption(tensor_diffusion->equation_id,
                                  RequirementKind::coefficient_regularity) ==
              specification.requirement_policies.end())
            report.add(DiagnosticCategory::analytical_policy,
                       tensor_diffusion->equation_id,
                       "scalar_coefficient_regularity_assumption",
                       "Declare the regularity assumptions for all scalar operator coefficients.");
          if (declared_assumption(tensor_diffusion->equation_id,
                                  RequirementKind::coercivity) ==
              specification.requirement_policies.end())
            report.add(DiagnosticCategory::analytical_policy,
                       tensor_diffusion->equation_id,
                       "scalar_coercivity_assumption",
                       "Declare the model author's coercivity assumption for the composed scalar form.");
          const auto state = std::find_if(
            specification.variables.begin(),
            specification.variables.end(),
            [](const VariableSpec &variable) {
              return variable.role == VariableRole::state;
            });
          if (state != specification.variables.end() &&
              !has_policy(state->id, RequirementKind::boundary_partition))
            report.add(DiagnosticCategory::analytical_policy,
                       state->id,
                       "scalar_boundary_partition_policy",
                       "Declare the selected fixed, natural, and transport boundary partition.");

          const auto robin = std::find_if(
            specification.residual_terms.begin(),
            specification.residual_terms.end(),
            [](const ResidualTermSpec &term) {
              return term.kind == ResidualTermKind::robin_bilinear;
            });
          if (robin != specification.residual_terms.end())
            {
              const auto conormal = selected_policy(
                robin->id, RequirementKind::conormal_flux);
              if (conormal == specification.requirement_policies.end())
                report.add(DiagnosticCategory::analytical_policy,
                           robin->id,
                           "conormal_flux_convention",
                           "Declare the selected conormal-flux sign and normal orientation.");
              else if (conormal->region_id != robin->region_id)
                report.add(DiagnosticCategory::structural,
                           robin->id,
                           "conormal_flux_region",
                           "Declare the conormal-flux convention on the Robin region.");
              const auto trace = selected_policy(robin->id,
                                                 RequirementKind::boundary_trace);
              if (trace == specification.requirement_policies.end())
                report.add(DiagnosticCategory::analytical_policy,
                           robin->id,
                           "robin_trace_realisation",
                           "Declare the selected Robin trace and face-quadrature realization.");
              else if (trace->region_id != robin->region_id)
                report.add(DiagnosticCategory::structural,
                           robin->id,
                           "robin_trace_region",
                           "Declare the Robin trace policy on its residual boundary region.");

              const auto conservative = std::find_if(
                specification.residual_terms.begin(),
                specification.residual_terms.end(),
                [](const ResidualTermSpec &term) {
                  return term.kind == ResidualTermKind::conservative_transport;
                });
              if (conservative != specification.residual_terms.end())
                {
                  const auto transport = selected_policy(
                    conservative->id,
                    RequirementKind::transport_boundary_trace);
                  if (transport == specification.requirement_policies.end())
                    report.add(DiagnosticCategory::analytical_policy,
                               conservative->id,
                               "transport_boundary_trace_policy",
                               "Declare the selected transport inflow/outflow trace interpretation.");
                  else if (transport->region_id != robin->region_id)
                    report.add(DiagnosticCategory::structural,
                               conservative->id,
                               "transport_boundary_trace_region",
                               "Declare the natural transport trace on the selected Robin region.");
                }
            }
        }

      for (const auto &term : specification.residual_terms)
        if (term.kind == ResidualTermKind::neumann_control)
          {
            const auto policy = selected_policy(term.id,
                                                RequirementKind::boundary_trace);
            if (policy == specification.requirement_policies.end())
              report.add(DiagnosticCategory::analytical_policy,
                         term.id,
                         "neumann_control_trace_realisation",
                         "Declare the selected Neumann-control trace pairing realization.");
            else if (policy->region_id != term.region_id)
              report.add(DiagnosticCategory::structural,
                         term.id,
                         "neumann_control_trace_region",
                         "Declare the Neumann trace policy on its residual boundary region.");
          }

      for (const auto &observation : specification.observations)
        if (observation.kind == ObservationKind::boundary_trace ||
            observation.kind == ObservationKind::weighted_boundary_trace)
          {
            const auto policy = selected_policy(observation.id,
                                                RequirementKind::boundary_trace);
            if (policy == specification.requirement_policies.end())
              report.add(DiagnosticCategory::analytical_policy,
                         observation.id,
                         "boundary_trace_realisation",
                         "Declare the selected boundary-trace observation realization.");
            else if (policy->region_id != observation.region_id)
              report.add(DiagnosticCategory::structural,
                         observation.id,
                         "boundary_trace_region",
                         "Declare the trace policy on the observation boundary region.");
          }
        else if (observation.kind == ObservationKind::normal_flux)
          {
            const auto flux_policy = selected_policy(
              observation.id, RequirementKind::conormal_flux);
            if (flux_policy == specification.requirement_policies.end())
              report.add(
                DiagnosticCategory::analytical_policy,
                observation.id,
                "normal_flux_orientation_policy",
                "Declare the outward-normal normal-flux convention on the observation boundary.");
            else if (flux_policy->region_id != observation.region_id)
              report.add(
                DiagnosticCategory::structural,
                observation.id,
                "normal_flux_orientation_region",
                "Declare the normal-flux orientation policy on the observation boundary region.");

            const auto evaluation_policy = selected_policy(
              observation.id, RequirementKind::analytic_quadrature_evaluation);
            if (evaluation_policy == specification.requirement_policies.end())
              report.add(
                DiagnosticCategory::analytical_policy,
                observation.id,
                "normal_flux_evaluation_policy",
                "Declare the selected face-quadrature normal-flux evaluation rule.");
            else if (evaluation_policy->region_id != observation.region_id)
              report.add(
                DiagnosticCategory::structural,
                observation.id,
                "normal_flux_evaluation_region",
                "Declare the normal-flux evaluation policy on the observation boundary region.");

            const auto transposition_policy = std::find_if(
              specification.requirement_policies.begin(),
              specification.requirement_policies.end(),
              [&specification](const RequirementPolicySpec &candidate_policy) {
                return candidate_policy.subject_id ==
                         specification.formulation.equation_id &&
                       candidate_policy.kind ==
                         RequirementKind::transposition_formulation &&
                       candidate_policy.status == RequirementStatus::provided &&
                       candidate_policy.scope ==
                         RequirementScope::continuous_semantics &&
                       !candidate_policy.selected_policy.empty();
              });
            if (transposition_policy == specification.requirement_policies.end())
              report.add(
                DiagnosticCategory::analytical_policy,
                specification.formulation.equation_id,
                "normal_flux_transposition_policy",
                "Declare the strong-state normal-flux adjoint transposition and very-weak formulation.");

            const auto regularity_policy = std::find_if(
              specification.requirement_policies.begin(),
              specification.requirement_policies.end(),
              [&specification](const RequirementPolicySpec &candidate_policy) {
                return candidate_policy.subject_id ==
                         specification.formulation.equation_id &&
                       candidate_policy.kind == RequirementKind::domain_regularity &&
                       candidate_policy.status == RequirementStatus::user_assumed &&
                       !candidate_policy.selected_policy.empty();
              });
            if (regularity_policy == specification.requirement_policies.end())
              report.add(
                DiagnosticCategory::analytical_policy,
                specification.formulation.equation_id,
                "normal_flux_domain_regularity",
                "Declare the domain regularity assumption required by the strong normal-flux state and very-weak adjoint.");
          }
    }
  };
} // namespace nmopt::semantic::v1
