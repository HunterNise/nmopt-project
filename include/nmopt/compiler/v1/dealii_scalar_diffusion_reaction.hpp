#pragma once

#include "nmopt/dealii/scalar_diffusion_reaction.hpp"
#include "nmopt/semantic/v1/problem_spec.hpp"

#include <deal.II/base/function.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/vector.h>

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
  // Value bindings are deliberately separate from semantic data declarations:
  // the semantic graph names immutable data and their ports, while this
  // deal.II compiler binds the concrete Function and scalar realizations.
  template <int dim>
  struct DealiiDataBindings
  {
    const dealii::Function<dim> &forcing;
    const dealii::Function<dim> &desired_state;
    double                       diffusion;
    double                       reaction;
    double                       regularisation_weight;
  };

  using CellwiseBoundValue = std::variant<double, dealii::Vector<double>>;

  struct CellwiseBoxDataBindings
  {
    CellwiseBoundValue lower;
    CellwiseBoundValue upper;
  };

  struct DealiiDiscretisationPolicy
  {
    enum class Execution
    {
      assembled,
      matrix_free
    };

    unsigned int state_degree = 1;
    Execution    execution = Execution::assembled;
  };

  // This small built-in registry is intentionally explicit. A semantic kind
  // is executable only when the v1 compiler has registered the corresponding
  // direct lowerer capability below; it never falls back to a nearby term.
  class DealiiLowererRegistryV1 final
  {
  public:
    bool
    has_residual_term_lowerer(const semantic::v1::ResidualTermKind kind) const
    {
      switch (kind)
        {
          case semantic::v1::ResidualTermKind::diffusion_reaction:
          case semantic::v1::ResidualTermKind::volume_source:
          case semantic::v1::ResidualTermKind::volume_control:
            return true;
        }
      return false;
    }

    bool
    has_observation_lowerer(const semantic::v1::ObservationKind kind) const
    {
      switch (kind)
        {
          case semantic::v1::ObservationKind::volume_restriction:
            return true;
        }
      return false;
    }

    bool
    has_loss_lowerer(const semantic::v1::LossKind kind) const
    {
      switch (kind)
        {
          case semantic::v1::LossKind::quadratic_tracking:
          case semantic::v1::LossKind::quadratic_control_regularisation:
            return true;
        }
      return false;
    }

    bool
    has_metric_lowerer(const semantic::v1::MetricKind kind) const
    {
      return kind == semantic::v1::MetricKind::l2;
    }

    bool
    has_constraint_lowerer(const semantic::v1::ConstraintKind kind) const
    {
      return kind == semantic::v1::ConstraintKind::cellwise_box;
    }
  };

  // A compiled v1 problem is a distinct semantic/compiler product. Its first
  // registered target is a newly constructed v0 direct executable instance;
  // retaining the direct object makes v0-v1 value, derivative, and DTO
  // comparisons possible without changing the v0 reference implementation.
  template <int dim>
  class CompiledProblem final
  {
  public:
    using DirectModel = dealii_backend::ScalarDiffusionReactionModel<dim>;

    CompiledProblem(std::unique_ptr<DirectModel> executable,
                    std::optional<dealii_backend::CellwiseBoxConstraint>
                      constraint)
      : executable_(std::move(executable))
      , constraint_(std::move(constraint))
    {
      contract::require(static_cast<bool>(executable_),
                        "A compiled v1 problem needs an executable model");
    }

    const DirectModel &
    executable_model() const
    {
      return *executable_;
    }

    dealii_backend::MassMetric
    control_l2_metric(
      dealii_backend::MassMetricSolveParameters parameters = {}) const
    {
      return executable_->control_l2_metric(parameters);
    }

    const dealii_backend::CellwiseBoxConstraint *
    control_constraint() const
    {
      return constraint_ ? &*constraint_ : nullptr;
    }

  private:
    std::unique_ptr<DirectModel> executable_;
    std::optional<dealii_backend::CellwiseBoxConstraint> constraint_;
  };

  template <int dim>
  struct CompilationResult
  {
    semantic::v1::ValidationReport diagnostics;
    std::unique_ptr<CompiledProblem<dim>> problem;

    bool
    succeeded() const
    {
      return diagnostics.valid() && static_cast<bool>(problem);
    }
  };

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
    CompilationResult<dim>
    compile(const semantic::v1::ProblemSpec &  specification,
            dealii::Triangulation<dim> &        triangulation,
            const DealiiDataBindings<dim> &     data,
            const DealiiDiscretisationPolicy &  policy = {},
            std::optional<CellwiseBoxDataBindings> bounds = std::nullopt) const
    {
      CompilationResult<dim> result;
      result.diagnostics = validate(specification, policy);
      const bool has_constraint = !specification.formulation.constraint_id.empty();
      if (has_constraint && !bounds)
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.constraint_id,
          "bound_data_binding",
          "Bind scalar constants or FE_DGQ(0) coefficient vectors for both bounds.");
      if (!result.diagnostics.valid())
        return result;

      const auto dirichlet_boundary_ids =
        selected_dirichlet_boundary_ids(specification);
      auto executable = std::make_unique<
        dealii_backend::ScalarDiffusionReactionModel<dim>>(
        triangulation,
        data.forcing,
        data.desired_state,
        data.diffusion,
        data.reaction,
        data.regularisation_weight,
        policy.state_degree,
        dirichlet_boundary_ids);

      std::optional<dealii_backend::CellwiseBoxConstraint> constraint;
      if (has_constraint)
        constraint = make_constraint(*executable, *bounds);

      result.problem = std::make_unique<CompiledProblem<dim>>(
        std::move(executable), std::move(constraint));
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

    template <int dim>
    static dealii_backend::CellwiseBoxConstraint
    make_constraint(
      const dealii_backend::ScalarDiffusionReactionModel<dim> &executable,
      const CellwiseBoxDataBindings &                          bounds)
    {
      const bool scalar_bounds =
        std::holds_alternative<double>(bounds.lower) &&
        std::holds_alternative<double>(bounds.upper);
      contract::require(
        scalar_bounds ||
          (std::holds_alternative<dealii::Vector<double>>(bounds.lower) &&
           std::holds_alternative<dealii::Vector<double>>(bounds.upper)),
        "The v1 cellwise box needs both bounds as scalars or both as FE_DGQ(0) vectors");
      if (scalar_bounds)
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

    DealiiLowererRegistryV1 registry_;
  };
} // namespace nmopt::compiler::v1
