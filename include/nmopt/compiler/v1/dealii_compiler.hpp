#pragma once

#include "nmopt/compiler/v1/compiled_problem.hpp"
#include "nmopt/compiler/v1/dealii_capabilities.hpp"
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
      using DirectModel = dealii_backend::ScalarDiffusionReactionModel<dim>;

      CompilationResultT<Backend> result;
      result.diagnostics = validate(specification, policy);
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
      const auto direct = std::make_shared<DirectModel>(
        triangulation,
        data.forcing,
        data.desired_state,
        data.diffusion,
        data.reaction,
        data.regularisation_weight,
        policy.state_degree,
        dirichlet_boundary_ids);

      const std::shared_ptr<const contract::MetricT<Backend>> metric =
        std::make_shared<dealii_backend::MassMetric>(
          direct->control_l2_metric(policy.control_metric_solve));
      std::shared_ptr<const contract::ConstraintT<Backend>> constraint;
      if (has_constraint)
        constraint = std::make_shared<dealii_backend::CellwiseBoxConstraint>(
          make_constraint(*direct, *bounds));

      const contract::StateAdjointSolversT<Backend> solvers{
        [direct](const contract::PrimalBlockT<Backend> &control) {
          return direct->solve_state(control);
        },
        [direct](const contract::PrimalBlockT<Backend> &full_point,
                 const contract::CovectorBlockT<Backend> &state_rhs) {
          return direct->solve_adjoint(full_point, state_rhs);
        }};
      const std::shared_ptr<const contract::ExecutableModelT<Backend>> executable =
        direct;
      result.problem = std::make_shared<const CompiledProblemT<Backend>>(
        executable,
        metric,
        constraint,
        solvers,
        make_manifest(specification, policy, has_constraint));
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

    template <int dim>
    static dealii_backend::CellwiseBoxConstraint
    make_constraint(
      const dealii_backend::ScalarDiffusionReactionModel<dim> &executable,
      const CellwiseBoxDataBindings &                          bounds)
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
                   "homogeneous_dirichlet_boundary_ids",
                   "Select a boundary region with at least one homogeneous Dirichlet id.");
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

      for (const auto &observation : specification.observations)
        {
          const auto region = find_region(specification, observation.region_id);
          if (region == nullptr || !region->is_full_domain)
            report.add(DiagnosticCategory::lowerability,
                       observation.id,
                       "full_domain_volume_observation",
                       "The first v1 lowerer supports full-domain volume restriction only.");
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
                  const bool                         has_constraint)
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
      manifest.data_rule = "deal.II Function at volume quadrature; scalar coefficients";
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
      manifest.lifting_realisation =
        "homogeneous full-vector Dirichlet rows; no inhomogeneous lifting";
      manifest.nullspace_policy =
        "not applicable: non-empty homogeneous Dirichlet boundary";
      manifest.state_adjoint_solve_policy =
        "serial CG with identity preconditioner for symmetric positive-definite operator";
      manifest.provenance = "DTO";
      manifest.region_ids = identifiers(specification.regions);
      manifest.space_ids = identifiers(specification.spaces);
      manifest.pairing_ids = identifiers(specification.pairings);
      manifest.variable_ids = identifiers(specification.variables);
      manifest.data_ids = identifiers(specification.data);
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

    DealiiLowererRegistryV1 registry_;
  };
} // namespace nmopt::compiler::v1
