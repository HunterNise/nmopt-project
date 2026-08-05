#pragma once

#include "nmopt/compiler/v1/compiled_problem.hpp"
#include "nmopt/compiler/v1/dealii_capabilities.hpp"
#include "nmopt/compiler/v1/dealii_fixed_dirichlet.hpp"
#include "nmopt/compiler/v1/dealii_types.hpp"
#include "nmopt/dealii/scalar_diffusion_reaction.hpp"
#include "nmopt/semantic/v1/validation.hpp"

#include <deal.II/grid/tria.h>

#include <algorithm>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <variant>
#include <vector>

namespace nmopt::compiler::v1
{
  class DealiiCompiler final
  {
  public:
    explicit DealiiCompiler(DealiiLowererRegistryV1 registry = {})
      : registry_(std::move(registry))
    {}

    semantic::v1::ValidationReport
    validate(const semantic::v1::ProblemSpec &  specification,
             const DealiiDiscretisationPolicy & policy) const
    {
      semantic::v1::ValidationReport report =
        semantic::v1::SemanticValidator().validate(specification);
      validate_lowerability(specification, policy, report);
      validate_formulation_capability(specification, report);
      return report;
    }

    template <int dim>
    CompilationResultT<dealii_backend::SerialBackend>
    compile(const semantic::v1::ProblemSpec &  specification,
            dealii::Triangulation<dim> &        triangulation,
            const DealiiDataBindings<dim> &     data,
            const DealiiDiscretisationPolicy &  policy = {},
            std::optional<CellwiseBoxDataBindings> bounds = std::nullopt) const
    {
      using Backend = dealii_backend::SerialBackend;
      CompilationResultT<Backend> result;
      result.diagnostics = validate(specification, policy);
      const bool uses_fixed_reconstruction =
        uses_fixed_dirichlet_reconstruction(specification);
      const auto *tracking_region = selected_tracking_region(specification);
      const bool uses_subdomain_observation =
        tracking_region != nullptr && !tracking_region->is_full_domain;
      const bool uses_assembled_v1_target =
        uses_fixed_reconstruction || uses_subdomain_observation;
      if (uses_fixed_reconstruction && !data.fixed_dirichlet_data)
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.state_variable_id,
          "fixed_dirichlet_data_binding",
          "Bind fixed Dirichlet Function data for the declared reconstruction.");
      if (!uses_fixed_reconstruction && data.fixed_dirichlet_data)
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.state_variable_id,
          "selected_fixed_dirichlet_reconstruction",
          "Declare the fixed-Dirichlet reconstruction before binding lifting data.");
      const bool has_constraint = !specification.formulation.constraint_id.empty();
      if (has_constraint && !bounds)
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.constraint_id,
          "bound_data_binding",
          "Bind scalar constants or FE_DGQ(0) coefficient vectors for both bounds.");
      if (bounds && !valid_bound_representation(*bounds))
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.constraint_id,
          "bound_data_representation",
          "Bind both cellwise bounds as scalars or both as FE_DGQ(0) vectors.");
      if (!result.diagnostics.valid())
        return result;

      const auto dirichlet_boundary_ids =
        selected_dirichlet_boundary_ids(specification);
      contract::require(tracking_region != nullptr,
                        "Validated v1 problem has no tracking observation region");
      std::shared_ptr<const contract::MetricT<Backend>> metric;
      std::shared_ptr<const contract::ConstraintT<Backend>> constraint;
      std::shared_ptr<const contract::ExecutableModelT<Backend>> executable;
      contract::StateAdjointSolversT<Backend> solvers;
      if (uses_assembled_v1_target)
        {
          using AssembledModel = detail::AssembledScalarDiffusionReactionModel<dim>;
          const auto assembled = std::make_shared<AssembledModel>(
            triangulation,
            data.forcing,
            data.desired_state,
            data.fixed_dirichlet_data,
            data.diffusion,
            data.reaction,
            data.regularisation_weight,
            policy.state_degree,
            dirichlet_boundary_ids,
            selected_tracking_material_ids(*tracking_region));
          metric = std::make_shared<dealii_backend::MassMetric>(
            assembled->control_l2_metric(policy.control_metric_solve));
          if (has_constraint)
            constraint = std::make_shared<dealii_backend::CellwiseBoxConstraint>(
              make_constraint(*assembled, *bounds));
          solvers = {
            [assembled](const contract::PrimalBlockT<Backend> &control) {
              return assembled->solve_state(control);
            },
            [assembled](const contract::PrimalBlockT<Backend> &full_point,
                        const contract::CovectorBlockT<Backend> &state_rhs) {
              return assembled->solve_adjoint(full_point, state_rhs);
            }};
          executable = assembled;
        }
      else
        {
          using DirectModel = dealii_backend::ScalarDiffusionReactionModel<dim>;
          const auto direct = std::make_shared<DirectModel>(
            triangulation,
            data.forcing,
            data.desired_state,
            data.diffusion,
            data.reaction,
            data.regularisation_weight,
            policy.state_degree,
            dirichlet_boundary_ids);
          metric = std::make_shared<dealii_backend::MassMetric>(
            direct->control_l2_metric(policy.control_metric_solve));
          if (has_constraint)
            constraint = std::make_shared<dealii_backend::CellwiseBoxConstraint>(
              make_constraint(*direct, *bounds));
          solvers = {
            [direct](const contract::PrimalBlockT<Backend> &control) {
              return direct->solve_state(control);
            },
            [direct](const contract::PrimalBlockT<Backend> &full_point,
                     const contract::CovectorBlockT<Backend> &state_rhs) {
              return direct->solve_adjoint(full_point, state_rhs);
            }};
          executable = direct;
        }
      result.problem = std::make_shared<const CompiledProblemT<Backend>>(
        executable,
        metric,
        constraint,
        solvers,
        make_manifest(specification,
                      policy,
                      has_constraint,
                      uses_fixed_reconstruction,
                      uses_assembled_v1_target,
                      *tracking_region));
      return result;
    }

  private:
    static const semantic::v1::RegionSpec *
    find_region(const semantic::v1::ProblemSpec &specification,
                const std::string &              id)
    {
      const auto region = std::find_if(
        specification.regions.begin(),
        specification.regions.end(),
        [&id](const semantic::v1::RegionSpec &candidate) {
          return candidate.id == id;
        });
      return region == specification.regions.end() ? nullptr : &*region;
    }

    static const semantic::v1::VariableSpec *
    find_variable(const semantic::v1::ProblemSpec &specification,
                  const std::string &              id)
    {
      const auto variable = std::find_if(
        specification.variables.begin(),
        specification.variables.end(),
        [&id](const semantic::v1::VariableSpec &candidate) {
          return candidate.id == id;
        });
      return variable == specification.variables.end() ? nullptr : &*variable;
    }

    static const semantic::v1::TransformationSpec *
    find_transformation(const semantic::v1::ProblemSpec &specification,
                        const std::string &              id)
    {
      const auto transformation = std::find_if(
        specification.transformations.begin(),
        specification.transformations.end(),
        [&id](const semantic::v1::TransformationSpec &candidate) {
          return candidate.id == id;
        });
      return transformation == specification.transformations.end()
               ? nullptr
               : &*transformation;
    }

    static const semantic::v1::ObservationSpec *
    find_observation(const semantic::v1::ProblemSpec &specification,
                     const std::string &              id)
    {
      const auto observation = std::find_if(
        specification.observations.begin(),
        specification.observations.end(),
        [&id](const semantic::v1::ObservationSpec &candidate) {
          return candidate.id == id;
        });
      return observation == specification.observations.end() ? nullptr :
                                                               &*observation;
    }

    static const semantic::v1::RegionSpec *
    selected_tracking_region(const semantic::v1::ProblemSpec &specification)
    {
      const auto loss = std::find_if(
        specification.losses.begin(),
        specification.losses.end(),
        [](const semantic::v1::LossSpec &candidate) {
          return candidate.kind == semantic::v1::LossKind::quadratic_tracking;
        });
      if (loss == specification.losses.end())
        return nullptr;
      const auto observation = find_observation(specification,
                                                loss->source_observation_id);
      return observation == nullptr ? nullptr :
                                    find_region(specification,
                                                observation->region_id);
    }

    static std::set<dealii::types::material_id>
    selected_tracking_material_ids(const semantic::v1::RegionSpec &region)
    {
      std::set<dealii::types::material_id> ids;
      for (const auto id : region.material_ids)
        ids.insert(static_cast<dealii::types::material_id>(id));
      return ids;
    }

    static bool
    uses_fixed_dirichlet_reconstruction(
      const semantic::v1::ProblemSpec &specification)
    {
      const auto state = find_variable(specification,
                                       specification.formulation.state_variable_id);
      if (state == nullptr || state->physical_field_transform_id.empty())
        return false;
      const auto transformation = find_transformation(
        specification, state->physical_field_transform_id);
      return transformation != nullptr &&
             transformation->kind ==
               semantic::v1::TransformationKind::fixed_dirichlet_reconstruction;
    }

    static std::set<dealii::types::boundary_id>
    selected_dirichlet_boundary_ids(
      const semantic::v1::ProblemSpec &specification)
    {
      const auto policy = std::find_if(
        specification.requirement_policies.begin(),
        specification.requirement_policies.end(),
        [&specification](const semantic::v1::RequirementPolicySpec &candidate) {
          return candidate.subject_id ==
                   specification.formulation.state_variable_id &&
                 candidate.kind == semantic::v1::RequirementKind::fixed_dirichlet;
        });
      contract::require(policy != specification.requirement_policies.end(),
                        "Validated v1 problem has no fixed Dirichlet policy");
      const auto region = find_region(specification, policy->region_id);
      contract::require(region != nullptr,
                        "Validated v1 fixed Dirichlet policy has no region");
      std::set<dealii::types::boundary_id> ids;
      for (const auto id : region->boundary_ids)
        ids.insert(static_cast<dealii::types::boundary_id>(id));
      return ids;
    }

    static bool
    valid_bound_representation(const CellwiseBoxDataBindings &bounds)
    {
      return (std::holds_alternative<double>(bounds.lower) &&
              std::holds_alternative<double>(bounds.upper)) ||
             (std::holds_alternative<dealii::Vector<double>>(bounds.lower) &&
              std::holds_alternative<dealii::Vector<double>>(bounds.upper));
    }

    template <typename Model>
    static dealii_backend::CellwiseBoxConstraint
    make_constraint(
      const Model &                    executable,
      const CellwiseBoxDataBindings &  bounds)
    {
      contract::require(valid_bound_representation(bounds),
                        "The v1 cellwise box needs compatible bound data");
      if (std::holds_alternative<double>(bounds.lower))
        return executable.control_l2_box_constraint(
          std::get<double>(bounds.lower), std::get<double>(bounds.upper));
      return executable.control_l2_box_constraint(
        std::get<dealii::Vector<double>>(bounds.lower),
        std::get<dealii::Vector<double>>(bounds.upper));
    }

    void
    validate_lowerability(const semantic::v1::ProblemSpec & specification,
                          const DealiiDiscretisationPolicy &policy,
                          semantic::v1::ValidationReport & report) const
    {
      using semantic::v1::DiagnosticCategory;
      if (policy.execution != DealiiDiscretisationPolicy::Execution::assembled)
        report.add(DiagnosticCategory::lowerability,
                   specification.id,
                   "assembled_execution",
                   "Select assembled execution; matrix-free lowering is not registered in v1.");
      if (policy.state_degree == 0)
        report.add(DiagnosticCategory::lowerability,
                   specification.id,
                   "FE_Q_state_degree",
                   "Select a scalar FE_Q state degree of at least one.");

      std::size_t full_volume_regions = 0;
      for (const auto &region : specification.regions)
        full_volume_regions += region.kind == semantic::v1::RegionKind::volume &&
                               region.is_full_domain;
      if (full_volume_regions != 1)
        report.add(DiagnosticCategory::lowerability,
                   specification.id,
                   "single_full_volume_region",
                   "The first v1 deal.II lowerer supports exactly one full volume region.");

      validate_registered_graph(specification, report);

      for (const auto &term : specification.residual_terms)
        if (!registry_.has_residual_term_lowerer(term.kind))
          report.add(DiagnosticCategory::lowerability,
                     term.id,
                     "registered_residual_term_lowerer",
                     "Register a lowerer for this residual term and its derivatives.");
      for (const auto &observation : specification.observations)
        if (!registry_.has_observation_lowerer(observation.kind))
          report.add(DiagnosticCategory::lowerability,
                     observation.id,
                     "registered_observation_lowerer",
                     "Register an observation value, JVP, and VJP lowerer.");
      for (const auto &loss : specification.losses)
        if (!registry_.has_loss_lowerer(loss.kind))
          report.add(DiagnosticCategory::lowerability,
                     loss.id,
                     "registered_loss_lowerer",
                     "Register a matching loss value and derivative lowerer.");
      for (const auto &metric : specification.metrics)
        if (!registry_.has_metric_lowerer(metric.kind))
          report.add(DiagnosticCategory::lowerability,
                     metric.id,
                     "registered_metric_lowerer",
                     "Register a metric realization with inverse apply.");
      for (const auto &constraint : specification.constraints)
        if (!registry_.has_constraint_lowerer(constraint.kind))
          report.add(DiagnosticCategory::lowerability,
                     constraint.id,
                     "registered_constraint_lowerer",
                     "Register the selected constraint projection realization.");
      for (const auto &transformation : specification.transformations)
        if (!registry_.has_transformation_lowerer(transformation.kind))
          report.add(DiagnosticCategory::lowerability,
                     transformation.id,
                     "registered_transformation_lowerer",
                     "Register value, JVP, and VJP lowering for this transformation.");

      const auto fixed_policy = std::find_if(
        specification.requirement_policies.begin(),
        specification.requirement_policies.end(),
        [&specification](const semantic::v1::RequirementPolicySpec &candidate) {
          return candidate.subject_id ==
                   specification.formulation.state_variable_id &&
                 candidate.kind == semantic::v1::RequirementKind::fixed_dirichlet;
        });
      const auto boundary = fixed_policy == specification.requirement_policies.end()
                              ? nullptr
                              : find_region(specification, fixed_policy->region_id);
      if (boundary == nullptr ||
          boundary->kind != semantic::v1::RegionKind::boundary ||
          boundary->boundary_ids.empty())
        report.add(DiagnosticCategory::lowerability,
                   specification.formulation.state_variable_id,
                   "fixed_dirichlet_boundary_ids",
                   "Select a boundary region with at least one fixed Dirichlet id.");
    }

    static void
    validate_registered_graph(const semantic::v1::ProblemSpec &specification,
                              semantic::v1::ValidationReport & report)
    {
      using semantic::v1::DataRole;
      using semantic::v1::DiagnosticCategory;
      using semantic::v1::LossKind;
      using semantic::v1::ResidualTermKind;

      const auto count_terms = [&specification](const ResidualTermKind kind) {
        return std::count_if(
          specification.residual_terms.begin(),
          specification.residual_terms.end(),
          [kind](const semantic::v1::ResidualTermSpec &term) {
            return term.kind == kind;
          });
      };
      if (count_terms(ResidualTermKind::diffusion_reaction) != 1 ||
          count_terms(ResidualTermKind::volume_source) != 1 ||
          count_terms(ResidualTermKind::volume_control) != 1)
        report.add(DiagnosticCategory::lowerability,
                   specification.id,
                   "complete_volume_residual_term_set",
                   "Declare exactly one diffusion-reaction, volume-source, and volume-control term.");

      const auto count_data = [&specification](const DataRole role) {
        return std::count_if(
          specification.data.begin(),
          specification.data.end(),
          [role](const semantic::v1::DataSpec &datum) {
            return datum.role == role;
          });
      };
      if (count_data(DataRole::forcing) != 1 ||
          count_data(DataRole::desired_state) != 1 ||
          count_data(DataRole::diffusion) != 1 ||
          count_data(DataRole::reaction) != 1 ||
          count_data(DataRole::regularisation_weight) != 1)
        report.add(DiagnosticCategory::lowerability,
                   specification.id,
                   "complete_volume_data_set",
                   "Declare one forcing, target, diffusion, reaction, and regularisation datum.");

      const auto count_losses = [&specification](const LossKind kind) {
        return std::count_if(
          specification.losses.begin(),
          specification.losses.end(),
          [kind](const semantic::v1::LossSpec &loss) {
            return loss.kind == kind;
          });
      };
      if (count_losses(LossKind::quadratic_tracking) != 1 ||
          count_losses(LossKind::quadratic_control_regularisation) != 1)
        report.add(DiagnosticCategory::lowerability,
                   specification.id,
                   "complete_quadratic_loss_set",
                   "Declare exactly one tracking and one control-regularisation loss.");

      const auto state = find_variable(
        specification, specification.formulation.state_variable_id);
      for (const auto &observation : specification.observations)
        {
          const auto region = find_region(specification, observation.region_id);
          if (region == nullptr || region->kind != semantic::v1::RegionKind::volume)
            report.add(DiagnosticCategory::lowerability,
                       observation.id,
                       "volume_observation_region",
                       "Select a registered volume observation region.");
          else if (state != nullptr &&
                   observation.input_variable_id == state->id &&
                   !region->is_full_domain && region->material_ids.empty())
            report.add(DiagnosticCategory::lowerability,
                       observation.id,
                       "material_subdomain_observation",
                       "Declare one or more material ids for the subdomain observation.");
          else if ((state == nullptr ||
                    observation.input_variable_id != state->id) &&
                   !region->is_full_domain)
            report.add(DiagnosticCategory::lowerability,
                       observation.id,
                       "full_domain_nonstate_observation",
                       "Only the state tracking observation supports material subdomains in v1.");
        }

      const bool has_constraint = !specification.formulation.constraint_id.empty();
      if (!has_constraint && !specification.constraints.empty())
        report.add(DiagnosticCategory::lowerability,
                   specification.id,
                   "selected_constraint_port",
                   "Select the declared constraint in the reduced formulation or remove it.");
      if (has_constraint &&
          (specification.constraints.size() != 1 ||
           count_data(DataRole::lower_bound) != 1 ||
           count_data(DataRole::upper_bound) != 1))
        report.add(DiagnosticCategory::lowerability,
                   specification.formulation.constraint_id,
                   "complete_cellwise_box_data_set",
                   "Declare one selected box plus one lower and one upper bound datum.");
    }

    static void
    validate_formulation_capability(
      const semantic::v1::ProblemSpec &specification,
      semantic::v1::ValidationReport & report)
    {
      using semantic::v1::DiagnosticCategory;
      if (specification.formulation.kind !=
          semantic::v1::FormulationKind::reduced_dto)
        report.add(DiagnosticCategory::formulation_capability,
                   specification.formulation.id,
                   "reduced_dto_formulation",
                   "Select the v1 reduced DTO formulation; all-at-once is not available.");
      if (specification.variables.size() != 2 ||
          specification.equations.size() != 1)
        report.add(DiagnosticCategory::formulation_capability,
                   specification.formulation.id,
                   "one_state_one_control_one_equation",
                   "The executable DTO contract currently supports one state, control, and equation block.");
      if (specification.formulation.constraint_id.empty())
        return;
      const auto constraint = std::find_if(
        specification.constraints.begin(),
        specification.constraints.end(),
        [&specification](const semantic::v1::ConstraintSpec &candidate) {
          return candidate.id == specification.formulation.constraint_id;
        });
      if (constraint == specification.constraints.end() ||
          constraint->kind != semantic::v1::ConstraintKind::cellwise_box)
        report.add(DiagnosticCategory::formulation_capability,
                   specification.formulation.id,
                   "l2_cellwise_projected_gradient",
                   "Use the registered cellwise L2 box constraint or omit the constraint.");
    }

    template <typename Component>
    static std::vector<std::string>
    identifiers(const std::vector<Component> &components)
    {
      std::vector<std::string> ids;
      ids.reserve(components.size());
      for (const auto &component : components)
        ids.push_back(component.id);
      return ids;
    }

    static CompilationManifest
    make_manifest(const semantic::v1::ProblemSpec &  specification,
                  const DealiiDiscretisationPolicy & policy,
                  const bool                         has_constraint,
                  const bool                         uses_fixed_reconstruction,
                  const bool                         uses_assembled_v1_target,
                  const semantic::v1::RegionSpec &   tracking_region)
    {
      CompilationManifest manifest;
      manifest.semantic_problem_id = specification.id;
      manifest.compiler_id =
        "nmopt.compiler.v1.dealii.scalar_diffusion_reaction";
      manifest.backend = "deal.II serial Vector<double>";
      manifest.execution = "assembled";
      manifest.state_space = "scalar FE_Q(" +
                             std::to_string(policy.state_degree) + ")";
      manifest.control_space = "FE_DGQ(0) on the state active-cell mesh";
      manifest.quadrature = "QGauss(" +
                            std::to_string(policy.state_degree + 2) + ")";
      manifest.dual_representation = "tested dual coefficients with dot pairing";
      manifest.data_rule =
        "analytic desired-state Function at selected QGauss(" +
        std::to_string(policy.state_degree + 2) + ") volume quadrature" +
        (uses_fixed_reconstruction
           ? "; fixed Dirichlet Function interpolated at boundary DoFs"
           : "; scalar coefficients and forcing Function at volume quadrature");
      manifest.observation_realisation =
        observation_realisation(tracking_region);
      manifest.metric_solve_policy =
        "serial CG: maximum iterations=" +
        std::to_string(policy.control_metric_solve.maximum_iterations) +
        ", relative tolerance=" +
        std::to_string(policy.control_metric_solve.relative_tolerance) +
        ", absolute tolerance=" +
        std::to_string(policy.control_metric_solve.absolute_tolerance);
      manifest.constraint_realisation = has_constraint
                                          ? "FE_DGQ(0) coefficientwise l2_cellwise clipping"
                                          : "none";
      manifest.lifting_realisation = uses_fixed_reconstruction
                                       ? "y_phys = P_h y_hat + ell_0,h; independent FE_Q coordinates, AffineConstraints reconstruction, and P_h^* pullbacks"
                                       : uses_assembled_v1_target
                                           ? "y_phys = P_h y_hat; independent FE_Q coordinates and AffineConstraints reconstruction"
                                           : "homogeneous full-vector Dirichlet rows; no inhomogeneous lifting";
      manifest.nullspace_policy =
        "not applicable: non-empty fixed Dirichlet boundary";
      manifest.state_adjoint_solve_policy =
        "serial CG with identity preconditioner for symmetric positive-definite operator";
      manifest.provenance = "DTO";
      manifest.region_ids = identifiers(specification.regions);
      manifest.space_ids = identifiers(specification.spaces);
      manifest.pairing_ids = identifiers(specification.pairings);
      manifest.variable_ids = identifiers(specification.variables);
      manifest.data_ids = identifiers(specification.data);
      manifest.transformation_ids = identifiers(specification.transformations);
      manifest.residual_term_ids = identifiers(specification.residual_terms);
      manifest.observation_ids = identifiers(specification.observations);
      manifest.loss_ids = identifiers(specification.losses);
      manifest.metric_ids = identifiers(specification.metrics);
      manifest.constraint_ids = identifiers(specification.constraints);
      for (const auto &requirement : specification.requirement_policies)
        manifest.declared_assumptions.push_back(
          requirement.id + ": " + requirement.selected_policy);
      return manifest;
    }

    static std::string
    observation_realisation(const semantic::v1::RegionSpec &region)
    {
      if (region.is_full_domain)
        return "full-domain volume restriction";
      std::string result = "material-id volume restriction: ";
      for (std::size_t index = 0; index < region.material_ids.size(); ++index)
        {
          if (index != 0)
            result += ",";
          result += std::to_string(region.material_ids[index]);
        }
      return result;
    }

    DealiiLowererRegistryV1 registry_;
  };
} // namespace nmopt::compiler::v1
