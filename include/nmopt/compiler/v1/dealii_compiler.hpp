#pragma once

#include "nmopt/compiler/v1/compiled_problem.hpp"
#include "nmopt/compiler/v1/dealii_capabilities.hpp"
#include "nmopt/compiler/v1/dealii_fixed_dirichlet.hpp"
#include "nmopt/compiler/v1/dealii_h1_control.hpp"
#include "nmopt/compiler/v1/dealii_neumann_boundary.hpp"
#include "nmopt/compiler/v1/dealii_types.hpp"
#include "nmopt/dealii/facewise_box_constraint.hpp"
#include "nmopt/dealii/scalar_diffusion_reaction.hpp"
#include "nmopt/semantic/v1/validation.hpp"

#include <deal.II/grid/tria.h>

#include <algorithm>
#include <cmath>
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
            std::optional<CellwiseBoxDataBindings> bounds = std::nullopt,
            std::optional<FacewiseBoxDataBindings> facewise_bounds = std::nullopt) const
    {
      using Backend = dealii_backend::SerialBackend;
      CompilationResultT<Backend> result;
      result.diagnostics = validate(specification, policy);
      const bool uses_fixed_reconstruction =
        uses_fixed_dirichlet_reconstruction(specification);
      const bool uses_neumann_boundary_control =
        uses_neumann_control(specification);
      const bool uses_mean_zero_gauge =
        uses_mean_zero_multiplier(specification);
      const bool uses_h1_control_regularisation =
        uses_h1_control_regularisation_loss(specification);
      const bool uses_h1_control_metric =
        selects_h1_control_metric(specification);
      const auto *tracking_region = selected_tracking_region(specification);
      const auto *control_boundary_region =
        selected_neumann_control_region(specification);
      const bool uses_subdomain_observation =
        tracking_region != nullptr && !tracking_region->is_full_domain;
      const bool uses_assembled_v1_target =
        !uses_neumann_boundary_control &&
        !uses_h1_control_regularisation &&
        (uses_fixed_reconstruction || uses_subdomain_observation);
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
      if (uses_mean_zero_gauge && std::abs(data.reaction) > 1e-14)
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.state_variable_id,
          "pure_neumann_zero_reaction",
          "Bind zero reaction for the selected pure-Neumann constant-nullspace policy.");
      if (uses_mean_zero_gauge && uses_fixed_reconstruction)
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.state_variable_id,
          "pure_neumann_without_fixed_reconstruction",
          "Remove fixed-Dirichlet reconstruction when selecting the pure-Neumann mean constraint.");
      const bool has_constraint = !specification.formulation.constraint_id.empty();
      if (has_constraint &&
          ((uses_neumann_boundary_control && !facewise_bounds) ||
           (!uses_neumann_boundary_control && !bounds)))
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.constraint_id,
          "bound_data_binding",
          uses_neumann_boundary_control
            ? "Bind scalar constants or exact facewise boundary-control vectors for both bounds."
            : "Bind scalar constants or FE_DGQ(0) coefficient vectors for both bounds.");
      if (bounds && !valid_bound_representation(*bounds))
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.constraint_id,
          "bound_data_representation",
          "Bind both cellwise bounds as scalars or both as FE_DGQ(0) vectors.");
      if (facewise_bounds && !valid_facewise_bound_representation(*facewise_bounds))
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.constraint_id,
          "facewise_bound_data_representation",
          "Bind both facewise bounds as scalars or both as exact boundary-control vectors.");
      if (uses_neumann_boundary_control && bounds)
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.constraint_id,
          "cellwise_bounds_for_boundary_control",
          "Bind the declared facewise box data for the boundary control.");
      if (!uses_neumann_boundary_control && facewise_bounds)
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.constraint_id,
          "facewise_bounds_for_volume_control",
          "Bind the declared cellwise box data for the volume control.");
      if (!result.diagnostics.valid())
        return result;

      const auto dirichlet_boundary_ids = uses_mean_zero_gauge
                                            ? std::set<dealii::types::boundary_id>{}
                                            : selected_dirichlet_boundary_ids(
                                                specification);
      contract::require(tracking_region != nullptr,
                        "Validated v1 problem has no tracking observation region");
      std::shared_ptr<const contract::MetricT<Backend>> metric;
      std::shared_ptr<const contract::ConstraintT<Backend>> constraint;
      std::shared_ptr<const contract::ExecutableModelT<Backend>> executable;
      contract::StateAdjointSolversT<Backend> solvers;
      if (uses_neumann_boundary_control)
        {
          contract::require(control_boundary_region != nullptr,
                            "Validated v1 problem has no Neumann control region");
          using BoundaryModel = detail::NeumannBoundaryControlModel<dim>;
          const auto boundary = std::make_shared<BoundaryModel>(
            triangulation,
            data.forcing,
            data.desired_state,
            data.diffusion,
            data.reaction,
            data.regularisation_weight,
            policy.state_degree,
            dirichlet_boundary_ids,
            boundary_ids(*control_boundary_region),
            boundary_ids(*tracking_region),
            uses_mean_zero_gauge
              ? BoundaryModel::StateGauge::mean_zero_multiplier
              : BoundaryModel::StateGauge::fixed_dirichlet);
          if (uses_mean_zero_gauge && !boundary->forcing_is_compatible())
            {
              result.diagnostics.add(
                semantic::v1::DiagnosticCategory::lowerability,
                "forcing",
                "pure_neumann_forcing_compatibility",
                "Bind forcing with zero discrete pairing against the constant null mode.");
              return result;
            }
          metric = std::make_shared<dealii_backend::MassMetric>(
            boundary->control_l2_metric(policy.control_metric_solve));
          if (has_constraint)
            constraint = std::make_shared<dealii_backend::FacewiseBoxConstraint>(
              make_facewise_constraint(*boundary, *facewise_bounds));
          solvers = {
            [boundary](const contract::PrimalBlockT<Backend> &control) {
              return boundary->solve_state(control);
            },
            [boundary](const contract::PrimalBlockT<Backend> &full_point,
                       const contract::CovectorBlockT<Backend> &state_rhs) {
              return boundary->solve_adjoint(full_point, state_rhs);
            }};
          executable = boundary;
        }
      else if (uses_h1_control_regularisation)
        {
          using H1Model = detail::H1ControlRegularisedModel<dim>;
          const auto h1_control = std::make_shared<H1Model>(
            triangulation,
            data.forcing,
            data.desired_state,
            data.diffusion,
            data.reaction,
            data.regularisation_weight,
            policy.state_degree,
            dirichlet_boundary_ids);
          metric = std::make_shared<dealii_backend::MassMetric>(
            uses_h1_control_metric
              ? h1_control->control_h1_metric(policy.control_metric_solve)
              : h1_control->control_l2_metric(policy.control_metric_solve));
          solvers = {
            [h1_control](const contract::PrimalBlockT<Backend> &control) {
              return h1_control->solve_state(control);
            },
            [h1_control](const contract::PrimalBlockT<Backend> &full_point,
                         const contract::CovectorBlockT<Backend> &state_rhs) {
              return h1_control->solve_adjoint(full_point, state_rhs);
            }};
          executable = h1_control;
        }
      else if (uses_assembled_v1_target)
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
                      uses_neumann_boundary_control,
                      uses_mean_zero_gauge,
                      uses_h1_control_regularisation,
                      uses_h1_control_metric,
                      *tracking_region,
                      control_boundary_region));
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
    uses_neumann_control(const semantic::v1::ProblemSpec &specification)
    {
      return std::any_of(
        specification.residual_terms.begin(),
        specification.residual_terms.end(),
        [](const semantic::v1::ResidualTermSpec &term) {
          return term.kind == semantic::v1::ResidualTermKind::neumann_control;
        });
    }

    static bool
    uses_mean_zero_multiplier(
      const semantic::v1::ProblemSpec &specification)
    {
      return std::any_of(
        specification.requirement_policies.begin(),
        specification.requirement_policies.end(),
        [&specification](const semantic::v1::RequirementPolicySpec &policy) {
          return policy.subject_id == specification.formulation.state_variable_id &&
                 policy.kind ==
                   semantic::v1::RequirementKind::mean_zero_multiplier &&
                 policy.status ==
                   semantic::v1::RequirementStatus::selected_discrete_realisation;
        });
    }

    static bool
    uses_h1_control_regularisation_loss(
      const semantic::v1::ProblemSpec &specification)
    {
      return std::any_of(
        specification.losses.begin(),
        specification.losses.end(),
        [](const semantic::v1::LossSpec &loss) {
          return loss.kind ==
                 semantic::v1::LossKind::quadratic_h1_control_regularisation;
        });
    }

    static bool
    selects_h1_control_metric(
      const semantic::v1::ProblemSpec &specification)
    {
      const auto metric = std::find_if(
        specification.metrics.begin(),
        specification.metrics.end(),
        [&specification](const semantic::v1::MetricSpec &candidate) {
          return candidate.id == specification.formulation.metric_id;
        });
      return metric != specification.metrics.end() &&
             metric->kind == semantic::v1::MetricKind::h1;
    }

    static const semantic::v1::RegionSpec *
    selected_neumann_control_region(
      const semantic::v1::ProblemSpec &specification)
    {
      const auto term = std::find_if(
        specification.residual_terms.begin(),
        specification.residual_terms.end(),
        [](const semantic::v1::ResidualTermSpec &candidate) {
          return candidate.kind == semantic::v1::ResidualTermKind::neumann_control;
        });
      return term == specification.residual_terms.end()
               ? nullptr
               : find_region(specification, term->region_id);
    }

    static std::set<dealii::types::boundary_id>
    boundary_ids(const semantic::v1::RegionSpec &region)
    {
      std::set<dealii::types::boundary_id> ids;
      for (const auto id : region.boundary_ids)
        ids.insert(static_cast<dealii::types::boundary_id>(id));
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
      return boundary_ids(*region);
    }

    static bool
    valid_bound_representation(const CellwiseBoxDataBindings &bounds)
    {
      return (std::holds_alternative<double>(bounds.lower) &&
              std::holds_alternative<double>(bounds.upper)) ||
             (std::holds_alternative<dealii::Vector<double>>(bounds.lower) &&
              std::holds_alternative<dealii::Vector<double>>(bounds.upper));
    }

    static bool
    valid_facewise_bound_representation(const FacewiseBoxDataBindings &bounds)
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

    template <typename Model>
    static dealii_backend::FacewiseBoxConstraint
    make_facewise_constraint(const Model &                  executable,
                             const FacewiseBoxDataBindings &bounds)
    {
      contract::require(valid_facewise_bound_representation(bounds),
                        "The v1 facewise box needs compatible bound data");
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

      const bool mean_zero_gauge = uses_mean_zero_multiplier(specification);
      const bool h1_control_regularisation =
        uses_h1_control_regularisation_loss(specification);
      const bool h1_control_metric = selects_h1_control_metric(specification);
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
      if (!mean_zero_gauge &&
          (boundary == nullptr ||
           boundary->kind != semantic::v1::RegionKind::boundary ||
           boundary->boundary_ids.empty()))
        report.add(DiagnosticCategory::lowerability,
                   specification.formulation.state_variable_id,
                   "fixed_dirichlet_boundary_ids",
                   "Select a boundary region with at least one fixed Dirichlet id.");

      if (mean_zero_gauge)
        {
          if (!uses_neumann_control(specification))
            report.add(
              DiagnosticCategory::lowerability,
              specification.formulation.state_variable_id,
              "pure_neumann_registered_residual",
              "Select the mean-zero multiplier only with the registered pure-Neumann boundary-control residual.");
          const auto mean_policy = std::find_if(
            specification.requirement_policies.begin(),
            specification.requirement_policies.end(),
            [&specification](const semantic::v1::RequirementPolicySpec &candidate) {
              return candidate.subject_id ==
                       specification.formulation.state_variable_id &&
                     candidate.kind ==
                       semantic::v1::RequirementKind::mean_zero_multiplier;
            });
          const auto mean_region =
            mean_policy == specification.requirement_policies.end()
              ? nullptr
              : find_region(specification, mean_policy->region_id);
          if (mean_region == nullptr ||
              mean_region->kind != semantic::v1::RegionKind::volume ||
              !mean_region->is_full_domain)
            report.add(
              DiagnosticCategory::lowerability,
              specification.formulation.state_variable_id,
              "pure_neumann_mean_constraint_region",
              "Place the mean-zero multiplier policy on the single full volume region.");
          if (boundary != nullptr)
            report.add(
              DiagnosticCategory::lowerability,
              specification.formulation.state_variable_id,
              "pure_neumann_without_fixed_dirichlet",
              "Do not declare a fixed Dirichlet policy with the pure-Neumann mean constraint.");
        }

      if (h1_control_regularisation)
        {
          const auto control = find_variable(
            specification, specification.formulation.control_variable_id);
          const auto space = std::find_if(
            specification.spaces.begin(),
            specification.spaces.end(),
            [control](const semantic::v1::SpaceSpec &candidate) {
              return control != nullptr && candidate.id == control->space_id;
            });
          if (control == nullptr ||
              space == specification.spaces.end() ||
              space->topology != semantic::v1::SpaceTopology::h1)
            report.add(
              DiagnosticCategory::lowerability,
              specification.formulation.control_variable_id,
              "h1_continuous_control_space",
              "Select the registered continuous H1 control space for H1 regularisation.");
          if (!specification.formulation.constraint_id.empty())
            report.add(
              DiagnosticCategory::lowerability,
              specification.formulation.constraint_id,
              "continuous_control_box_constraint",
              "Do not select the cellwise or facewise box with the continuous H1 control realization.");
          if (uses_fixed_dirichlet_reconstruction(specification) ||
              uses_neumann_control(specification))
            report.add(
              DiagnosticCategory::lowerability,
              specification.id,
              "h1_regularisation_registered_combination",
              "The first H1-control regularisation target supports the homogeneous volume-control graph only.");
          const auto h1_tracking_region = selected_tracking_region(specification);
          if (h1_tracking_region == nullptr ||
              !h1_tracking_region->is_full_domain)
            report.add(
              DiagnosticCategory::lowerability,
              specification.id,
              "h1_regularisation_full_domain_tracking",
              "The first H1-control regularisation target supports full-domain tracking only.");
        }

      if (h1_control_metric && !h1_control_regularisation)
        report.add(
          DiagnosticCategory::lowerability,
          specification.formulation.metric_id,
          "h1_metric_registered_control_space",
          "Select the registered continuous H1-control target before requesting the H1 metric.");

      if (!uses_neumann_control(specification))
        return;
      const auto control_region = selected_neumann_control_region(specification);
      const auto tracking_region = selected_tracking_region(specification);
      if (control_region == nullptr ||
          control_region->kind != semantic::v1::RegionKind::boundary ||
          control_region->boundary_ids.empty())
        report.add(DiagnosticCategory::lowerability,
                   specification.id,
                   "neumann_control_boundary_ids",
                   "Select a marked boundary region with at least one Neumann control id.");
      if (tracking_region == nullptr ||
          tracking_region->kind != semantic::v1::RegionKind::boundary ||
          tracking_region->boundary_ids.empty())
        report.add(DiagnosticCategory::lowerability,
                   specification.id,
                   "boundary_tracking_region",
                   "Select a marked boundary region for the state trace observation.");
      if (control_region != nullptr && boundary != nullptr)
        for (const auto control_id : control_region->boundary_ids)
          if (std::find(boundary->boundary_ids.begin(),
                        boundary->boundary_ids.end(),
                        control_id) != boundary->boundary_ids.end())
            report.add(DiagnosticCategory::lowerability,
                       control_region->id,
                       "neumann_control_dirichlet_overlap",
                       "Use boundary ids not fixed by the homogeneous Dirichlet realization.");
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
      const bool boundary_control = uses_neumann_control(specification);
      const bool complete_residual = boundary_control
        ? count_terms(ResidualTermKind::diffusion_reaction) == 1 &&
            count_terms(ResidualTermKind::volume_source) == 1 &&
            count_terms(ResidualTermKind::neumann_control) == 1 &&
            count_terms(ResidualTermKind::volume_control) == 0
        : count_terms(ResidualTermKind::diffusion_reaction) == 1 &&
            count_terms(ResidualTermKind::volume_source) == 1 &&
            count_terms(ResidualTermKind::volume_control) == 1 &&
            count_terms(ResidualTermKind::neumann_control) == 0;
      if (!complete_residual)
        report.add(DiagnosticCategory::lowerability,
                   specification.id,
                   boundary_control ? "complete_neumann_boundary_residual_term_set"
                                    : "complete_volume_residual_term_set",
                   boundary_control
                     ? "Declare exactly one diffusion-reaction, volume-source, and Neumann-control term."
                     : "Declare exactly one diffusion-reaction, volume-source, and volume-control term.");

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
      const bool h1_control_regularisation =
        uses_h1_control_regularisation_loss(specification);
      const bool complete_control_loss = h1_control_regularisation
        ? count_losses(LossKind::quadratic_h1_control_regularisation) == 1 &&
            count_losses(LossKind::quadratic_control_regularisation) == 0
        : count_losses(LossKind::quadratic_control_regularisation) == 1 &&
            count_losses(LossKind::quadratic_h1_control_regularisation) == 0;
      if (count_losses(LossKind::quadratic_tracking) != 1 ||
          !complete_control_loss || specification.losses.size() != 2)
        report.add(DiagnosticCategory::lowerability,
                   specification.id,
                   "complete_registered_loss_set",
                   h1_control_regularisation
                     ? "Declare exactly one tracking and one H1 control-regularisation loss."
                     : "Declare exactly one tracking and one L2 control-regularisation loss.");

      const auto selected_metric = std::find_if(
        specification.metrics.begin(),
        specification.metrics.end(),
        [&specification](const semantic::v1::MetricSpec &metric) {
          return metric.id == specification.formulation.metric_id;
        });
      if (specification.metrics.size() != 1 ||
          selected_metric == specification.metrics.end() ||
          (selected_metric->kind != semantic::v1::MetricKind::l2 &&
           selected_metric->kind != semantic::v1::MetricKind::h1))
        report.add(DiagnosticCategory::lowerability,
                   specification.formulation.metric_id,
                   "selected_registered_metric",
                   "Select exactly one registered L2 or H1 control metric.");

      const auto state = find_variable(
        specification, specification.formulation.state_variable_id);
      if (boundary_control)
        {
          const auto state_trace = std::find_if(
            specification.observations.begin(),
            specification.observations.end(),
            [](const semantic::v1::ObservationSpec &observation) {
              return observation.kind ==
                     semantic::v1::ObservationKind::boundary_trace;
            });
          const auto control_restriction = std::find_if(
            specification.observations.begin(),
            specification.observations.end(),
            [](const semantic::v1::ObservationSpec &observation) {
              return observation.kind ==
                     semantic::v1::ObservationKind::boundary_restriction;
            });
          if (specification.observations.size() != 2 ||
              state_trace == specification.observations.end() ||
              control_restriction == specification.observations.end())
            report.add(DiagnosticCategory::lowerability,
                       specification.id,
                       "complete_boundary_observation_set",
                       "Declare one boundary state trace and one boundary control restriction.");
          if (state_trace != specification.observations.end())
            {
              const auto region = find_region(specification, state_trace->region_id);
              if (region == nullptr ||
                  region->kind != semantic::v1::RegionKind::boundary ||
                  region->boundary_ids.empty())
                report.add(DiagnosticCategory::lowerability,
                           state_trace->id,
                           "boundary_trace_observation_region",
                           "Select marked boundary ids for the state trace observation.");
            }
          if (control_restriction != specification.observations.end())
            {
              const auto region = find_region(specification,
                                              control_restriction->region_id);
              const auto control_region = selected_neumann_control_region(specification);
              if (region == nullptr || control_region == nullptr ||
                  region->id != control_region->id)
                report.add(DiagnosticCategory::lowerability,
                           control_restriction->id,
                           "boundary_control_observation_region",
                           "Restrict the facewise control on its Neumann control boundary.");
            }
        }
      else
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
                   boundary_control ? "complete_facewise_box_data_set"
                                    : "complete_cellwise_box_data_set",
                   boundary_control
                     ? "Declare one selected facewise box plus one lower and one upper bound datum."
                     : "Declare one selected box plus one lower and one upper bound datum.");
      if (has_constraint)
        {
          const auto constraint = std::find_if(
            specification.constraints.begin(),
            specification.constraints.end(),
            [&specification](const semantic::v1::ConstraintSpec &candidate) {
              return candidate.id == specification.formulation.constraint_id;
            });
          if (constraint != specification.constraints.end() &&
              constraint->kind != (boundary_control
                                    ? semantic::v1::ConstraintKind::facewise_box
                                    : semantic::v1::ConstraintKind::cellwise_box))
            report.add(DiagnosticCategory::lowerability,
                       constraint->id,
                       "control_layout_constraint_realisation",
                       "Use a facewise box for boundary control and a cellwise box for volume control.");
        }
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
          (constraint->kind != semantic::v1::ConstraintKind::cellwise_box &&
           constraint->kind != semantic::v1::ConstraintKind::facewise_box))
        report.add(DiagnosticCategory::formulation_capability,
                   specification.formulation.id,
                   "l2_coefficientwise_projected_gradient",
                   "Use the registered cellwise or facewise L2 box constraint, or omit the constraint.");
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
                  const bool                         uses_neumann_boundary_control,
                  const bool                         uses_mean_zero_gauge,
                  const bool                         uses_h1_control_regularisation,
                  const bool                         uses_h1_control_metric,
                  const semantic::v1::RegionSpec &   tracking_region,
                  const semantic::v1::RegionSpec *   control_boundary_region)
    {
      CompilationManifest manifest;
      manifest.semantic_problem_id = specification.id;
      manifest.compiler_id =
        "nmopt.compiler.v1.dealii.scalar_diffusion_reaction";
      manifest.backend = "deal.II serial Vector<double>";
      manifest.execution = "assembled";
      manifest.state_space = "scalar FE_Q(" +
                             std::to_string(policy.state_degree) + ")";
      manifest.control_space = uses_neumann_boundary_control
                                 ? "one facewise-constant coefficient per marked state boundary face"
                                 : uses_h1_control_regularisation
                                     ? "continuous scalar FE_Q(" +
                                         std::to_string(policy.state_degree) +
                                         ") on the state mesh"
                                 : "FE_DGQ(0) on the state active-cell mesh";
      manifest.quadrature = "QGauss(" +
                            std::to_string(policy.state_degree + 2) + ")";
      manifest.dual_representation = "tested dual coefficients with dot pairing";
      manifest.data_rule = uses_neumann_boundary_control
        ? "analytic desired-state Function at selected QGauss(" +
            std::to_string(policy.state_degree + 2) +
            ") boundary face quadrature; scalar coefficients and forcing Function at volume quadrature"
        : "analytic desired-state Function at selected QGauss(" +
            std::to_string(policy.state_degree + 2) + ") volume quadrature" +
            (uses_fixed_reconstruction
              ? "; fixed Dirichlet Function interpolated at boundary DoFs"
              : "; scalar coefficients and forcing Function at volume quadrature");
      manifest.observation_realisation = uses_neumann_boundary_control
        ? boundary_observation_realisation(tracking_region)
        : observation_realisation(tracking_region);
      manifest.metric_solve_policy =
        std::string("serial CG for ") +
        (uses_h1_control_metric ? "h1_continuous Riesz map" : "L2 Riesz map") +
        ": maximum iterations=" +
        std::to_string(policy.control_metric_solve.maximum_iterations) +
        ", relative tolerance=" +
        std::to_string(policy.control_metric_solve.relative_tolerance) +
        ", absolute tolerance=" +
        std::to_string(policy.control_metric_solve.absolute_tolerance);
      manifest.constraint_realisation = !has_constraint
                                          ? "none"
                                          : uses_neumann_boundary_control
                                              ? "facewise-constant coefficientwise l2_facewise clipping"
                                              : "FE_DGQ(0) coefficientwise l2_cellwise clipping";
      manifest.lifting_realisation = uses_mean_zero_gauge
                                       ? "none; pure-Neumann state uses an explicit mean-zero gauge"
                                       : uses_fixed_reconstruction
                                       ? "y_phys = P_h y_hat + ell_0,h; independent FE_Q coordinates, AffineConstraints reconstruction, and P_h^* pullbacks"
                                       : uses_assembled_v1_target
                                           ? "y_phys = P_h y_hat; independent FE_Q coordinates and AffineConstraints reconstruction"
                                           : "homogeneous full-vector Dirichlet rows; no inhomogeneous lifting";
      manifest.nullspace_policy = uses_mean_zero_gauge
        ? "one mean-zero Lagrange multiplier; discrete constant-mode compatibility is enforced for forcing and every control state load"
        : "not applicable: non-empty fixed Dirichlet boundary";
      manifest.state_adjoint_solve_policy = uses_mean_zero_gauge
        ? "serial SparseDirectUMFPACK on the augmented symmetric state and adjoint saddle systems"
        : "serial CG with identity preconditioner for symmetric positive-definite operator";
      manifest.provenance = "DTO";
      if (uses_neumann_boundary_control && control_boundary_region != nullptr)
        manifest.declared_assumptions.push_back(
          "neumann_control_realisation: facewise-constant FEFaceValues pairing on boundary ids " +
          boundary_id_list(*control_boundary_region));
      if (uses_h1_control_regularisation)
        manifest.declared_assumptions.push_back(
          "h1_control_regularisation: alpha/2 u^T (M_u + K_u) u; search metric=" +
          std::string(uses_h1_control_metric ? "h1_continuous" : "l2_continuous"));
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

    static std::string
    boundary_id_list(const semantic::v1::RegionSpec &region)
    {
      std::string result;
      for (std::size_t index = 0; index < region.boundary_ids.size(); ++index)
        {
          if (index != 0)
            result += ",";
          result += std::to_string(region.boundary_ids[index]);
        }
      return result;
    }

    static std::string
    boundary_observation_realisation(const semantic::v1::RegionSpec &region)
    {
      return "boundary trace restriction on boundary ids " +
             boundary_id_list(region);
    }

    DealiiLowererRegistryV1 registry_;
  };
} // namespace nmopt::compiler::v1
