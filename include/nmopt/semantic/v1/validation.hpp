#pragma once

#include "nmopt/semantic/v1/types.hpp"

#include <algorithm>
#include <string>
#include <unordered_map>
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
      const auto metrics = index(specification.metrics, report, "metric");
      const auto constraints = index(specification.constraints, report,
                                     "constraint");
      index(specification.requirement_policies, report, "requirement policy");

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
      validate_observations(specification, variables, regions, spaces, pairings,
                            report);
      validate_losses(specification, observations, pairings, data, report);
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

    static void
    validate_regions(const ProblemSpec &specification, ValidationReport &report)
    {
      for (const auto &region : specification.regions)
        {
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
        if (!contains(spaces, pairing.primal_space_id) ||
            !contains(spaces, pairing.covector_space_id))
          report.add(DiagnosticCategory::structural,
                     pairing.id,
                     "pairing_space_ports",
                     "Reference declared primal and covector spaces.");
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
        if (!datum.space_id.empty() && !contains(spaces, datum.space_id))
          report.add(DiagnosticCategory::structural,
                     datum.id,
                     "data_space_port",
                     "Reference a declared data space or leave scalar constants unspaced.");
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
          if (input == variables.end() ||
              input->second->role != VariableRole::state ||
              output == spaces.end() ||
              output->second->role != SpaceRole::state ||
              fixed_data == data.end() ||
              fixed_data->second->kind != DataKind::function ||
              fixed_data->second->role != DataRole::fixed_dirichlet_lifting)
            report.add(
              DiagnosticCategory::structural,
              transformation.id,
              "fixed_dirichlet_reconstruction_ports",
              "Connect the fixed reconstruction to one state variable, its state space, and fixed Function data.");
          if (input != variables.end() && output != spaces.end() &&
              input->second->space_id != transformation.output_space_id)
            report.add(DiagnosticCategory::structural,
                       transformation.id,
                       "fixed_dirichlet_reconstruction_space",
                       "Reconstruct the physical field in the state variable's declared space.");
          if (fixed_data != data.end() &&
              fixed_data->second->space_id != transformation.output_space_id)
            report.add(DiagnosticCategory::structural,
                       transformation.id,
                       "fixed_dirichlet_lifting_space",
                       "Declare fixed lifting data in the reconstructed state space.");
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
              pairing->second->primal_space_id != equation.test_space_id)
            report.add(DiagnosticCategory::structural,
                       equation.id,
                       "equation_test_pairing",
                       "Reference a pairing whose primal port is the equation test space.");
          for (const auto &term_id : equation.residual_term_ids)
            {
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
          if (!contains(equations, term.equation_id))
            report.add(DiagnosticCategory::structural,
                       term.id,
                       "residual_term_target_equation",
                       "Reference a declared target equation block.");
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
          if (term.kind != ResidualTermKind::neumann_control &&
              !term.region_id.empty())
            report.add(DiagnosticCategory::structural,
                       term.id,
                       "volume_term_has_no_boundary_region",
                       "Leave the region port empty for the registered volume term.");
          if (term.kind != ResidualTermKind::neumann_control)
            continue;

          const auto region = regions.find(term.region_id);
          if (region == regions.end() ||
              region->second->kind != RegionKind::boundary)
            report.add(DiagnosticCategory::structural,
                       term.id,
                       "neumann_control_boundary_region",
                       "Declare the Neumann trace pairing on a boundary region.");
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
                       "neumann_control_space_region",
                       "Place the Neumann control space on the declared boundary region.");
        }
    }

    static void
    validate_observations(const ProblemSpec &          specification,
                          const Index<VariableSpec> &  variables,
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
            observation.kind == ObservationKind::boundary_restriction;
          const RegionKind expected_region = boundary_observation
                                               ? RegionKind::boundary
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
                                              "observation_volume_region",
                       boundary_observation
                         ? "The registered boundary observation needs a boundary region."
                         : "The registered volume restriction needs a volume region.");
          const auto input = variables.find(observation.input_variable_id);
          if (observation.kind == ObservationKind::boundary_trace &&
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
          const auto pairing = pairings.find(observation.output_pairing_id);
          if (pairing == pairings.end() ||
              pairing->second->primal_space_id != observation.output_space_id)
            report.add(DiagnosticCategory::structural,
                       observation.id,
                       "observation_output_pairing",
                       "Reference a pairing whose primal port is the observation output.");
        }
    }

    static void
    validate_losses(const ProblemSpec &             specification,
                    const Index<ObservationSpec> &  observations,
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
              pairing->second->primal_space_id !=
                observation->second->output_space_id)
            report.add(DiagnosticCategory::structural,
                       loss.id,
                       "loss_pairing",
                       "Reference the pairing for the loss observation output.");
          validate_loss_signature(loss, data, report);
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
              pairing->second->primal_space_id != variable->second->space_id)
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
          case ResidualTermKind::diffusion_reaction:
            valid = term.variable_ids.size() == 1 && term.data_ids.size() == 2 &&
                    has_variable_role(term.variable_ids, variables,
                                      VariableRole::state) &&
                    has_role(term.data_ids, data, DataRole::diffusion) &&
                    has_role(term.data_ids, data, DataRole::reaction);
            break;
          case ResidualTermKind::parameter_diffusion_reaction:
            valid = term.variable_ids.size() == 2 && term.data_ids.size() == 1 &&
                    has_variable_role(term.variable_ids, variables,
                                      VariableRole::state) &&
                    has_variable_role(term.variable_ids, variables,
                                      VariableRole::parameter) &&
                    has_role(term.data_ids, data, DataRole::reaction);
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
      if (!contains(metrics, formulation.metric_id))
        report.add(DiagnosticCategory::structural,
                   formulation.id,
                   "formulation_metric",
                   "Reference a declared search metric.");
      if (!formulation.constraint_id.empty() &&
          !contains(constraints, formulation.constraint_id))
        report.add(DiagnosticCategory::structural,
                   formulation.id,
                   "formulation_constraint",
                   "Reference a declared constraint or leave this port empty.");
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
            const bool has_mean_zero_multiplier =
              has_policy(variable.id, RequirementKind::mean_zero_multiplier);
            if (!has_fixed_dirichlet && !has_mean_zero_multiplier)
              report.add(
                DiagnosticCategory::analytical_policy,
                variable.id,
                "state_uniqueness_realisation",
                "Declare either the selected fixed-Dirichlet policy or the mean-zero multiplier policy.");
            if (has_fixed_dirichlet && has_mean_zero_multiplier)
              report.add(
                DiagnosticCategory::structural,
                variable.id,
                "state_uniqueness_policy_conflict",
                "Select exactly one state uniqueness policy for the declared residual.");
          }

      for (const auto &datum : specification.data)
        if (datum.role == DataRole::desired_state &&
            !has_policy(datum.id,
                        RequirementKind::analytic_quadrature_evaluation))
          report.add(DiagnosticCategory::analytical_policy,
                     datum.id,
                     "desired_state_data_rule",
                     "Declare the analytic quadrature target-data realization.");

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
        if (observation.kind == ObservationKind::boundary_trace)
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
    }
  };
} // namespace nmopt::semantic::v1
