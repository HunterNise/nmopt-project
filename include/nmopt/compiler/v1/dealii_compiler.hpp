#pragma once

#include "nmopt/compiler/v1/compiled_problem.hpp"
#include "nmopt/compiler/v1/dealii_capabilities.hpp"
#include "nmopt/compiler/v1/dealii_coefficient_identification.hpp"
#include "nmopt/compiler/v1/dealii_dirichlet_control.hpp"
#include "nmopt/compiler/v1/dealii_fixed_dirichlet.hpp"
#include "nmopt/compiler/v1/dealii_h1_control.hpp"
#include "nmopt/compiler/v1/dealii_neumann_boundary.hpp"
#include "nmopt/compiler/v1/dealii_scalar_plan.hpp"
#include "nmopt/compiler/v1/dealii_types.hpp"
#include "nmopt/dealii/facewise_box_constraint.hpp"
#include "nmopt/dealii/scalar_diffusion_reaction.hpp"
#include "nmopt/semantic/v1/validation.hpp"

#include <deal.II/grid/tria.h>

#include <algorithm>
#include <cmath>
#include <limits>
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
    explicit DealiiCompiler(
      DealiiCapabilityRegistryV1 capabilities = {},
      DealiiScalarLowererRegistryV1 scalar_registry = {})
      : capabilities_(std::move(capabilities))
      , scalar_planner_(std::move(scalar_registry))
    {}

    semantic::v1::ValidationReport
    validate(const semantic::v1::ProblemSpec &  specification,
             const DealiiDiscretisationPolicy & policy) const
    {
      auto resolution = semantic::v1::SemanticResolver().resolve(specification);
      semantic::v1::ValidationReport report = std::move(resolution.diagnostics);
      if (!report.valid())
        return report;
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
      return compile_impl(specification,
                          triangulation,
                          data,
                          policy,
                          std::move(bounds),
                          std::move(facewise_bounds),
                          {},
                          "caller-owned triangulation",
                          false);
    }

    template <int dim>
    CompilationResultT<dealii_backend::SerialBackend>
    compile(
      const semantic::v1::ProblemSpec &specification,
      const std::shared_ptr<DealiiCompilationSession<dim>> &session,
      const DealiiDataBindings<dim> &data,
      const DealiiDiscretisationPolicy &policy = {},
      std::optional<CellwiseBoxDataBindings> bounds = std::nullopt,
      std::optional<FacewiseBoxDataBindings> facewise_bounds = std::nullopt) const
    {
      contract::require(static_cast<bool>(session),
                        "The deal.II compiler needs a compilation session");
      return compile_impl(specification,
                          session->mutable_triangulation(),
                          data,
                          policy,
                          std::move(bounds),
                          std::move(facewise_bounds),
                          session,
                          session->mesh_provenance(),
                          true);
    }

  private:
    template <int dim>
    CompilationResultT<dealii_backend::SerialBackend>
    compile_impl(
            const semantic::v1::ProblemSpec &  specification,
            dealii::Triangulation<dim> &        triangulation,
            const DealiiDataBindings<dim> &     data,
            const DealiiDiscretisationPolicy &  policy,
            std::optional<CellwiseBoxDataBindings> bounds,
            std::optional<FacewiseBoxDataBindings> facewise_bounds,
            std::shared_ptr<const void>             lifetime_owner,
            std::string                             mesh_provenance,
            const bool                              owns_mesh) const
    {
      using Backend = dealii_backend::SerialBackend;
      CompilationResultT<Backend> result;
      auto resolution = semantic::v1::SemanticResolver().resolve(specification);
      result.diagnostics = std::move(resolution.diagnostics);
      if (!result.diagnostics.valid())
        return result;
      validate_lowerability(specification, policy, result.diagnostics);
      validate_formulation_capability(specification, result.diagnostics);
      if (!result.diagnostics.valid())
        return result;
      const bool uses_fixed_reconstruction =
        uses_fixed_dirichlet_reconstruction(specification);
      const bool uses_dirichlet_control =
        uses_dirichlet_control_lifting(specification);
      const bool uses_neumann_boundary_control =
        uses_neumann_control(specification);
      const bool uses_mean_zero_gauge =
        uses_mean_zero_multiplier(specification);
      const bool uses_h1_control_regularisation =
        uses_h1_control_regularisation_loss(specification);
      const bool uses_h1_control_metric =
        selects_h1_control_metric(specification);
      const bool uses_coefficient_identification =
        uses_parameter_diffusion_residual(specification);
      const auto *tracking_region = selected_tracking_region(specification);
      const auto *control_boundary_region =
        uses_dirichlet_control
          ? selected_dirichlet_control_region(specification)
          : selected_neumann_control_region(specification);
      const bool uses_subdomain_observation =
        tracking_region != nullptr && !tracking_region->is_full_domain;
      const bool uses_assembled_v1_target =
        !uses_neumann_boundary_control &&
        !uses_h1_control_regularisation &&
        !uses_coefficient_identification &&
        !uses_dirichlet_control &&
        (uses_fixed_reconstruction || uses_subdomain_observation);
      std::optional<ScalarLoweringPlan> scalar_plan;
      if (uses_assembled_v1_target)
        {
          auto planned = scalar_planner_.plan(*resolution.problem);
          for (const auto &diagnostic : planned.diagnostics.diagnostics())
            result.diagnostics.add(diagnostic.category,
                                   diagnostic.component_id,
                                   diagnostic.capability,
                                   diagnostic.remedy);
          scalar_plan = std::move(planned.plan);
        }
      if (triangulation.n_active_cells() == 0)
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.id,
          "nonempty_triangulation",
          "Compile on a triangulation with at least one active cell.");
      if (data.provenance.forcing.empty())
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          "forcing",
          "forcing_binding_provenance",
          "Supply a stable provenance label for the forcing Function binding.");
      if (data.provenance.desired_state.empty())
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          "desired_state",
          "desired_state_binding_provenance",
          "Supply a stable provenance label for the desired-state Function binding.");
      if (uses_fixed_reconstruction &&
          data.provenance.fixed_dirichlet_data.empty())
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          "fixed_dirichlet_data",
          "fixed_dirichlet_binding_provenance",
          "Supply a stable provenance label for the fixed-Dirichlet Function binding.");
      if (!uses_coefficient_identification && !data.diffusion)
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          "diffusion",
          "diffusion_data_binding",
          "Bind the constant diffusion coefficient selected by this graph.");
      if (!uses_coefficient_identification && data.diffusion &&
          (!std::isfinite(*data.diffusion) || *data.diffusion <= 0.0))
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          "diffusion",
          "positive_finite_diffusion_binding",
          "Bind a positive finite constant diffusion coefficient.");
      if (!std::isfinite(data.reaction) || data.reaction < 0.0)
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          "reaction",
          "nonnegative_finite_reaction_binding",
          "Bind a nonnegative finite reaction coefficient.");
      if (!std::isfinite(data.regularisation_weight) ||
          data.regularisation_weight <= 0.0)
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          "regularisation",
          "positive_finite_regularisation_binding",
          "Bind a positive finite regularisation weight.");
      if (!valid_metric_solve_policy(policy.control_metric_solve))
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.metric_id,
          "valid_metric_solve_policy",
          "Select positive finite metric-solve tolerances and a positive iteration limit.");
      if (!dealii_backend::valid(policy.state_solve))
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.id,
          "valid_state_solve_policy",
          "Select positive finite state-solve tolerances.");
      if (!dealii_backend::valid(policy.adjoint_solve))
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.id,
          "valid_adjoint_solve_policy",
          "Select positive finite adjoint-solve tolerances.");
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
      if (uses_coefficient_identification && bounds &&
          !has_strictly_positive_lower_bound(*bounds))
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.constraint_id,
          "positive_parameter_lower_bound",
          "Bind a strictly positive scalar or every strictly positive cellwise lower parameter bound.");
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
      if (!has_constraint && (bounds || facewise_bounds))
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.id,
          "unselected_bound_data",
          "Remove bound data when the reduced formulation selects no constraint.");
      if (bounds && valid_bound_representation(*bounds))
        validate_cellwise_bound_values(*bounds,
                                       triangulation.n_active_cells(),
                                       specification.formulation.constraint_id,
                                       result.diagnostics);
      if (facewise_bounds && valid_facewise_bound_representation(*facewise_bounds) &&
          control_boundary_region != nullptr)
        validate_facewise_bound_values(
          *facewise_bounds,
          count_boundary_faces(triangulation,
                               boundary_ids(*control_boundary_region)),
          specification.formulation.constraint_id,
          result.diagnostics);
      if (!result.diagnostics.valid())
        return result;

      const auto dirichlet_boundary_ids = uses_mean_zero_gauge
                                            ? std::set<dealii::types::boundary_id>{}
                                            : uses_dirichlet_control
                                              ? selected_dirichlet_control_boundary_ids(
                                                  specification)
                                            : selected_dirichlet_boundary_ids(
                                                specification);
      if (triangulation.has_hanging_nodes())
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.id,
          "conforming_mesh_without_hanging_nodes",
          "Compile the registered serial targets on a mesh without hanging-node relations.");
      if (!uses_mean_zero_gauge && !uses_dirichlet_control &&
          !contains_all_boundary_ids(triangulation, dirichlet_boundary_ids))
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.state_variable_id,
          "fixed_dirichlet_boundary_presence",
          "Select fixed-Dirichlet boundary ids present on the compiled mesh.");
      if (uses_neumann_boundary_control && control_boundary_region != nullptr &&
          !contains_all_boundary_ids(triangulation,
                                     boundary_ids(*control_boundary_region)))
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          control_boundary_region->id,
          "neumann_control_boundary_presence",
          "Select Neumann-control boundary ids present on the compiled mesh.");
      if (uses_neumann_boundary_control && tracking_region != nullptr &&
          !contains_all_boundary_ids(triangulation,
                                     boundary_ids(*tracking_region)))
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          tracking_region->id,
          "boundary_observation_presence",
          "Select boundary-observation ids present on the compiled mesh.");
      if (!result.diagnostics.valid())
        return result;
      if (uses_dirichlet_control &&
          !controls_complete_exterior_boundary(triangulation,
                                               dirichlet_boundary_ids))
        {
          result.diagnostics.add(
            semantic::v1::DiagnosticCategory::lowerability,
            specification.formulation.state_variable_id,
            "complete_dirichlet_control_boundary",
            "Select every exterior boundary id for the registered nodal Dirichlet lifting; partial boundaries, interfaces, and undeclared corner policies are not supported.");
          return result;
        }
      contract::require(tracking_region != nullptr,
                        "Validated v1 problem has no tracking observation region");
      std::shared_ptr<const dealii_backend::MassMetric> metric;
      std::shared_ptr<const contract::ConstraintT<Backend>> constraint;
      std::shared_ptr<const contract::ExecutableModelT<Backend>> executable;
      contract::StateAdjointSolversT<Backend> solvers;
      ConstraintRealisation constraint_realisation = ConstraintRealisation::none;
      if (uses_neumann_boundary_control)
        {
          contract::require(control_boundary_region != nullptr,
                            "Validated v1 problem has no Neumann control region");
          using BoundaryModel = detail::NeumannBoundaryControlModel<dim>;
          const auto boundary = std::make_shared<BoundaryModel>(
            triangulation,
            data.forcing,
            data.desired_state,
            *data.diffusion,
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
            {
              constraint =
                std::make_shared<dealii_backend::FacewiseBoxConstraint>(
                  make_facewise_constraint(*boundary, *facewise_bounds, *metric));
              constraint_realisation = ConstraintRealisation::facewise_l2;
            }
          solvers = {
            [boundary, solve_policy = policy.state_solve](
              const contract::PrimalBlockT<Backend> &control) {
              return boundary->solve_state_with_report(control, solve_policy);
            },
            [boundary, solve_policy = policy.adjoint_solve](
              const contract::PrimalBlockT<Backend> &full_point,
              const contract::CovectorBlockT<Backend> &state_rhs) {
              return boundary->solve_adjoint_with_report(full_point,
                                                         state_rhs,
                                                         solve_policy);
            }};
          executable = boundary;
        }
      else if (uses_dirichlet_control)
        {
          using DirichletModel = detail::DirichletControlLiftingModel<dim>;
          const auto dirichlet = std::make_shared<DirichletModel>(
            triangulation,
            data.forcing,
            data.desired_state,
            *data.diffusion,
            data.reaction,
            data.regularisation_weight,
            policy.state_degree,
            dirichlet_boundary_ids);
          metric = std::make_shared<dealii_backend::MassMetric>(
            dirichlet->control_l2_metric(policy.control_metric_solve));
          solvers = {
            [dirichlet, solve_policy = policy.state_solve](
              const contract::PrimalBlockT<Backend> &control) {
              return dirichlet->solve_state_with_report(control, solve_policy);
            },
            [dirichlet, solve_policy = policy.adjoint_solve](
              const contract::PrimalBlockT<Backend> &full_point,
              const contract::CovectorBlockT<Backend> &state_rhs) {
              return dirichlet->solve_adjoint_with_report(full_point,
                                                          state_rhs,
                                                          solve_policy);
            }};
          executable = dirichlet;
        }
      else if (uses_h1_control_regularisation)
        {
          using H1Model = detail::H1ControlRegularisedModel<dim>;
          const auto h1_control = std::make_shared<H1Model>(
            triangulation,
            data.forcing,
            data.desired_state,
            *data.diffusion,
            data.reaction,
            data.regularisation_weight,
            policy.state_degree,
            dirichlet_boundary_ids);
          metric = std::make_shared<dealii_backend::MassMetric>(
            uses_h1_control_metric
              ? h1_control->control_h1_metric(policy.control_metric_solve)
              : h1_control->control_l2_metric(policy.control_metric_solve));
          solvers = {
            [h1_control, solve_policy = policy.state_solve](
              const contract::PrimalBlockT<Backend> &control) {
              return h1_control->solve_state_with_report(control, solve_policy);
            },
            [h1_control, solve_policy = policy.adjoint_solve](
              const contract::PrimalBlockT<Backend> &full_point,
              const contract::CovectorBlockT<Backend> &state_rhs) {
              return h1_control->solve_adjoint_with_report(full_point,
                                                           state_rhs,
                                                           solve_policy);
            }};
          executable = h1_control;
        }
      else if (uses_coefficient_identification)
        {
          using CoefficientModel = detail::CoefficientIdentificationModel<dim>;
          const auto coefficient = std::make_shared<CoefficientModel>(
            triangulation,
            data.forcing,
            data.desired_state,
            data.reaction,
            data.regularisation_weight,
            policy.state_degree,
            dirichlet_boundary_ids);
          metric = std::make_shared<dealii_backend::MassMetric>(
            coefficient->parameter_l2_metric(policy.control_metric_solve));
          if (has_constraint)
            {
              constraint =
                std::make_shared<dealii_backend::CellwiseBoxConstraint>(
                  make_parameter_constraint(*coefficient, *bounds, *metric));
              constraint_realisation =
                ConstraintRealisation::cellwise_parameter_l2;
            }
          solvers = {
            [coefficient, solve_policy = policy.state_solve](
              const contract::PrimalBlockT<Backend> &parameter) {
              return coefficient->solve_state_with_report(parameter,
                                                          solve_policy);
            },
            [coefficient, solve_policy = policy.adjoint_solve](
              const contract::PrimalBlockT<Backend> &full_point,
              const contract::CovectorBlockT<Backend> &state_rhs) {
              return coefficient->solve_adjoint_with_report(full_point,
                                                            state_rhs,
                                                            solve_policy);
            }};
          executable = coefficient;
        }
      else if (uses_assembled_v1_target)
        {
          using AssembledModel = detail::ScalarComponentModel<dim>;
          const auto assembled = std::make_shared<AssembledModel>(
            triangulation,
            data.forcing,
            data.desired_state,
            data.fixed_dirichlet_data,
            *data.diffusion,
            data.reaction,
            data.regularisation_weight,
            policy.state_degree,
            *scalar_plan);
          metric = std::make_shared<dealii_backend::MassMetric>(
            assembled->control_l2_metric(policy.control_metric_solve));
          if (has_constraint)
            {
              constraint =
                std::make_shared<dealii_backend::CellwiseBoxConstraint>(
                  make_constraint(*assembled, *bounds, *metric));
              constraint_realisation = ConstraintRealisation::cellwise_l2;
            }
          solvers = {
            [assembled, solve_policy = policy.state_solve](
              const contract::PrimalBlockT<Backend> &control) {
              return assembled->solve_state_with_report(control, solve_policy);
            },
            [assembled, solve_policy = policy.adjoint_solve](
              const contract::PrimalBlockT<Backend> &full_point,
              const contract::CovectorBlockT<Backend> &state_rhs) {
              return assembled->solve_adjoint_with_report(full_point,
                                                          state_rhs,
                                                          solve_policy);
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
            *data.diffusion,
            data.reaction,
            data.regularisation_weight,
            policy.state_degree,
            dirichlet_boundary_ids);
          metric = std::make_shared<dealii_backend::MassMetric>(
            direct->control_l2_metric(policy.control_metric_solve));
          if (has_constraint)
            {
              constraint =
                std::make_shared<dealii_backend::CellwiseBoxConstraint>(
                  make_constraint(*direct, *bounds, *metric));
              constraint_realisation = ConstraintRealisation::cellwise_l2;
            }
          solvers = {
            [direct, solve_policy = policy.state_solve](
              const contract::PrimalBlockT<Backend> &control) {
              return direct->solve_state_with_report(control, solve_policy);
            },
            [direct, solve_policy = policy.adjoint_solve](
              const contract::PrimalBlockT<Backend> &full_point,
              const contract::CovectorBlockT<Backend> &state_rhs) {
              return direct->solve_adjoint_with_report(full_point,
                                                       state_rhs,
                                                       solve_policy);
            }};
          executable = direct;
        }
      const CompiledTargetKind target_kind = uses_mean_zero_gauge
                                               ? CompiledTargetKind::pure_neumann
                                             : uses_neumann_boundary_control
                                               ? CompiledTargetKind::neumann_boundary
                                             : uses_dirichlet_control
                                               ? CompiledTargetKind::dirichlet_control
                                             : uses_coefficient_identification
                                               ? CompiledTargetKind::coefficient_identification
                                             : uses_h1_control_regularisation
                                               ? (uses_h1_control_metric
                                                    ? CompiledTargetKind::h1_control_h1_metric
                                                    : CompiledTargetKind::h1_control_l2_metric)
                                             : uses_assembled_v1_target
                                               ? CompiledTargetKind::assembled_volume
                                               : CompiledTargetKind::direct_volume;
      result.problem = std::make_shared<const CompiledProblemT<Backend>>(
        executable,
        metric,
        constraint,
        solvers,
        make_manifest(specification,
                      policy,
                      constraint_realisation,
                      target_kind,
                      *tracking_region,
                      control_boundary_region,
                      data,
                      bounds,
                      facewise_bounds,
                      triangulation,
                      mesh_provenance,
                      owns_mesh,
                      *executable,
                      *metric,
                      scalar_plan ? &*scalar_plan : nullptr),
        std::move(lifetime_owner));
      return result;
    }

  private:
    enum class ConstraintRealisation
    {
      none,
      cellwise_l2,
      cellwise_parameter_l2,
      facewise_l2
    };

    enum class CompiledTargetKind
    {
      direct_volume,
      assembled_volume,
      neumann_boundary,
      pure_neumann,
      dirichlet_control,
      h1_control_l2_metric,
      h1_control_h1_metric,
      coefficient_identification
    };

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
    uses_parameter_diffusion_residual(
      const semantic::v1::ProblemSpec &specification)
    {
      return std::any_of(
        specification.residual_terms.begin(),
        specification.residual_terms.end(),
        [](const semantic::v1::ResidualTermSpec &term) {
          return term.kind ==
                 semantic::v1::ResidualTermKind::parameter_diffusion_reaction;
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

    static const semantic::v1::RegionSpec *
    selected_dirichlet_control_region(
      const semantic::v1::ProblemSpec &specification)
    {
      const auto policy = std::find_if(
        specification.requirement_policies.begin(),
        specification.requirement_policies.end(),
        [&specification](const semantic::v1::RequirementPolicySpec &candidate) {
          return candidate.subject_id ==
                   specification.formulation.state_variable_id &&
                 candidate.kind ==
                   semantic::v1::RequirementKind::controlled_dirichlet;
        });
      return policy == specification.requirement_policies.end()
               ? nullptr
               : find_region(specification, policy->region_id);
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

    static bool
    uses_dirichlet_control_lifting(
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
               semantic::v1::TransformationKind::dirichlet_control_lifting;
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

    static std::set<dealii::types::boundary_id>
    selected_dirichlet_control_boundary_ids(
      const semantic::v1::ProblemSpec &specification)
    {
      const auto region = selected_dirichlet_control_region(specification);
      contract::require(region != nullptr,
                        "Validated v1 problem has no controlled Dirichlet policy");
      return boundary_ids(*region);
    }

    template <int dim>
    static bool
    controls_complete_exterior_boundary(
      const dealii::Triangulation<dim> &              triangulation,
      const std::set<dealii::types::boundary_id> &controlled_ids)
    {
      bool has_boundary_face = false;
      for (auto cell = triangulation.begin_active();
           cell != triangulation.end();
           ++cell)
        for (unsigned int face = 0;
             face < dealii::GeometryInfo<dim>::faces_per_cell;
             ++face)
          if (cell->face(face)->at_boundary())
            {
              has_boundary_face = true;
              if (controlled_ids.count(cell->face(face)->boundary_id()) == 0)
                return false;
            }
      return has_boundary_face;
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
    valid_metric_solve_policy(
      const dealii_backend::MassMetricSolveParameters &policy)
    {
      return policy.maximum_iterations > 0 &&
             std::isfinite(policy.relative_tolerance) &&
             policy.relative_tolerance > 0.0 &&
             std::isfinite(policy.absolute_tolerance) &&
             policy.absolute_tolerance > 0.0;
    }

    template <int dim>
    static std::size_t
    count_boundary_faces(
      const dealii::Triangulation<dim> &              triangulation,
      const std::set<dealii::types::boundary_id> &boundary_ids)
    {
      std::size_t count = 0;
      for (auto cell = triangulation.begin_active();
           cell != triangulation.end();
           ++cell)
        for (unsigned int face = 0;
             face < dealii::GeometryInfo<dim>::faces_per_cell;
             ++face)
          if (cell->face(face)->at_boundary() &&
              boundary_ids.count(cell->face(face)->boundary_id()) != 0)
            ++count;
      return count;
    }

    template <int dim>
    static bool
    contains_all_boundary_ids(
      const dealii::Triangulation<dim> &              triangulation,
      const std::set<dealii::types::boundary_id> &requested_ids)
    {
      std::set<dealii::types::boundary_id> found_ids;
      for (auto cell = triangulation.begin_active();
           cell != triangulation.end();
           ++cell)
        for (unsigned int face = 0;
             face < dealii::GeometryInfo<dim>::faces_per_cell;
             ++face)
          if (cell->face(face)->at_boundary() &&
              requested_ids.count(cell->face(face)->boundary_id()) != 0)
            found_ids.insert(cell->face(face)->boundary_id());
      return found_ids == requested_ids;
    }

    static void
    validate_cellwise_bound_values(
      const CellwiseBoxDataBindings &bounds,
      const std::size_t              expected_size,
      const std::string &            component_id,
      semantic::v1::ValidationReport &report)
    {
      validate_bound_values(bounds.lower,
                            bounds.upper,
                            expected_size,
                            component_id,
                            "cellwise_bound_layout",
                            report);
    }

    static void
    validate_facewise_bound_values(
      const FacewiseBoxDataBindings &bounds,
      const std::size_t              expected_size,
      const std::string &            component_id,
      semantic::v1::ValidationReport &report)
    {
      validate_bound_values(bounds.lower,
                            bounds.upper,
                            expected_size,
                            component_id,
                            "facewise_bound_layout",
                            report);
    }

    template <typename BoundValue>
    static void
    validate_bound_values(const BoundValue &lower,
                          const BoundValue &upper,
                          const std::size_t expected_size,
                          const std::string &component_id,
                          const std::string &layout_capability,
                          semantic::v1::ValidationReport &report)
    {
      using semantic::v1::DiagnosticCategory;
      if (std::holds_alternative<double>(lower))
        {
          const double lower_value = std::get<double>(lower);
          const double upper_value = std::get<double>(upper);
          if (!std::isfinite(lower_value) || !std::isfinite(upper_value))
            report.add(DiagnosticCategory::lowerability,
                       component_id,
                       "finite_bound_values",
                       "Bind finite lower and upper values.");
          else if (lower_value > upper_value)
            report.add(DiagnosticCategory::lowerability,
                       component_id,
                       "ordered_bound_values",
                       "Bind lower values that do not exceed upper values.");
          return;
        }

      const auto &lower_values = std::get<dealii::Vector<double>>(lower);
      const auto &upper_values = std::get<dealii::Vector<double>>(upper);
      if (static_cast<std::size_t>(lower_values.size()) != expected_size ||
          static_cast<std::size_t>(upper_values.size()) != expected_size)
        {
          report.add(DiagnosticCategory::lowerability,
                     component_id,
                     layout_capability,
                     "Bind lower and upper vectors with the exact compiled decision layout.");
          return;
        }
      for (dealii::Vector<double>::size_type index = 0;
           index < lower_values.size();
           ++index)
        if (!std::isfinite(lower_values[index]) ||
            !std::isfinite(upper_values[index]))
          {
            report.add(DiagnosticCategory::lowerability,
                       component_id,
                       "finite_bound_values",
                       "Bind finite lower and upper values.");
            return;
          }
        else if (lower_values[index] > upper_values[index])
          {
            report.add(DiagnosticCategory::lowerability,
                       component_id,
                       "ordered_bound_values",
                       "Bind lower values that do not exceed upper values.");
            return;
          }
    }

    static bool
    has_strictly_positive_lower_bound(const CellwiseBoxDataBindings &bounds)
    {
      if (!valid_bound_representation(bounds))
        return false;
      if (std::holds_alternative<double>(bounds.lower))
        return std::isfinite(std::get<double>(bounds.lower)) &&
               std::get<double>(bounds.lower) > 0.0;
      const auto &lower = std::get<dealii::Vector<double>>(bounds.lower);
      return std::all_of(lower.begin(), lower.end(), [](const double value) {
        return std::isfinite(value) && value > 0.0;
      });
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
      const CellwiseBoxDataBindings &  bounds,
      const dealii_backend::MassMetric &projection_metric)
    {
      contract::require(valid_bound_representation(bounds),
                        "The v1 cellwise box needs compatible bound data");
      if (std::holds_alternative<double>(bounds.lower))
        return executable.control_l2_box_constraint(
          std::get<double>(bounds.lower),
          std::get<double>(bounds.upper),
          projection_metric);
      return executable.control_l2_box_constraint(
        std::get<dealii::Vector<double>>(bounds.lower),
        std::get<dealii::Vector<double>>(bounds.upper),
        projection_metric);
    }

    template <typename Model>
    static dealii_backend::CellwiseBoxConstraint
    make_parameter_constraint(
      const Model &                    executable,
      const CellwiseBoxDataBindings &  bounds,
      const dealii_backend::MassMetric &projection_metric)
    {
      contract::require(valid_bound_representation(bounds),
                        "The v1 parameter box needs compatible bound data");
      contract::require(has_strictly_positive_lower_bound(bounds),
                        "The v1 parameter box needs a strictly positive lower bound");
      if (std::holds_alternative<double>(bounds.lower))
        return executable.parameter_l2_box_constraint(
          std::get<double>(bounds.lower),
          std::get<double>(bounds.upper),
          projection_metric);
      return executable.parameter_l2_box_constraint(
        std::get<dealii::Vector<double>>(bounds.lower),
        std::get<dealii::Vector<double>>(bounds.upper),
        projection_metric);
    }

    template <typename Model>
    static dealii_backend::FacewiseBoxConstraint
    make_facewise_constraint(const Model &                  executable,
                             const FacewiseBoxDataBindings &bounds,
                             const dealii_backend::MassMetric &projection_metric)
    {
      contract::require(valid_facewise_bound_representation(bounds),
                        "The v1 facewise box needs compatible bound data");
      if (std::holds_alternative<double>(bounds.lower))
        return executable.control_l2_box_constraint(
          std::get<double>(bounds.lower),
          std::get<double>(bounds.upper),
          projection_metric);
      return executable.control_l2_box_constraint(
        std::get<dealii::Vector<double>>(bounds.lower),
        std::get<dealii::Vector<double>>(bounds.upper),
        projection_metric);
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
        if (!capabilities_.has_residual_term_lowerer(term.kind))
          report.add(DiagnosticCategory::lowerability,
                     term.id,
                     "registered_residual_term_lowerer",
                     "Register a lowerer for this residual term and its derivatives.");
      for (const auto &observation : specification.observations)
        if (!capabilities_.has_observation_lowerer(observation.kind))
          report.add(DiagnosticCategory::lowerability,
                     observation.id,
                     "registered_observation_lowerer",
                     "Register an observation value, JVP, and VJP lowerer.");
      for (const auto &loss : specification.losses)
        if (!capabilities_.has_loss_lowerer(loss.kind))
          report.add(DiagnosticCategory::lowerability,
                     loss.id,
                     "registered_loss_lowerer",
                     "Register a matching loss value and derivative lowerer.");
      for (const auto &metric : specification.metrics)
        if (!capabilities_.has_metric_lowerer(metric.kind))
          report.add(DiagnosticCategory::lowerability,
                     metric.id,
                     "registered_metric_lowerer",
                     "Register a metric realization with inverse apply.");
      for (const auto &constraint : specification.constraints)
        if (!capabilities_.has_constraint_lowerer(constraint.kind))
          report.add(DiagnosticCategory::lowerability,
                     constraint.id,
                     "registered_constraint_lowerer",
                     "Register the selected constraint projection realization.");
      for (const auto &transformation : specification.transformations)
        if (!capabilities_.has_transformation_lowerer(transformation.kind))
          report.add(DiagnosticCategory::lowerability,
                     transformation.id,
                     "registered_transformation_lowerer",
                     "Register value, JVP, and VJP lowering for this transformation.");

      const bool mean_zero_gauge = uses_mean_zero_multiplier(specification);
      const bool h1_control_regularisation =
        uses_h1_control_regularisation_loss(specification);
      const bool h1_control_metric = selects_h1_control_metric(specification);
      const bool coefficient_identification =
        uses_parameter_diffusion_residual(specification);
      const bool dirichlet_control_lifting =
        uses_dirichlet_control_lifting(specification);
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
      const auto controlled_boundary =
        selected_dirichlet_control_region(specification);
      if (!mean_zero_gauge && !dirichlet_control_lifting &&
          (boundary == nullptr ||
           boundary->kind != semantic::v1::RegionKind::boundary ||
           boundary->boundary_ids.empty()))
        report.add(DiagnosticCategory::lowerability,
                   specification.formulation.state_variable_id,
                   "fixed_dirichlet_boundary_ids",
                   "Select a boundary region with at least one fixed Dirichlet id.");
      if (dirichlet_control_lifting &&
          (controlled_boundary == nullptr ||
           controlled_boundary->kind != semantic::v1::RegionKind::boundary ||
           controlled_boundary->boundary_ids.empty()))
        report.add(
          DiagnosticCategory::lowerability,
          specification.formulation.state_variable_id,
          "controlled_dirichlet_boundary_ids",
          "Select a non-empty exterior boundary region for the registered Dirichlet-control lifting.");

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
              uses_neumann_control(specification) || dirichlet_control_lifting)
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

      if (coefficient_identification)
        {
          const auto parameter = find_variable(
            specification, specification.formulation.control_variable_id);
          const auto space = std::find_if(
            specification.spaces.begin(),
            specification.spaces.end(),
            [parameter](const semantic::v1::SpaceSpec &candidate) {
              return parameter != nullptr && candidate.id == parameter->space_id;
            });
          if (parameter == nullptr ||
              parameter->role != semantic::v1::VariableRole::parameter ||
              space == specification.spaces.end() ||
              space->topology != semantic::v1::SpaceTopology::l2)
            report.add(
              DiagnosticCategory::lowerability,
              specification.formulation.control_variable_id,
              "cellwise_parameter_space",
              "Select the registered cellwise L2 parameter space for coefficient identification.");
          if (specification.formulation.constraint_id.empty())
            report.add(
              DiagnosticCategory::lowerability,
              specification.formulation.control_variable_id,
              "positive_parameter_constraint",
              "Select the registered positive cellwise parameter box.");
          if (uses_fixed_dirichlet_reconstruction(specification) ||
              uses_neumann_control(specification) || h1_control_regularisation ||
              dirichlet_control_lifting)
            report.add(
              DiagnosticCategory::lowerability,
              specification.id,
              "coefficient_identification_registered_combination",
              "The first coefficient-identification target supports the homogeneous full-domain volume graph only.");
          const auto coefficient_tracking_region =
            selected_tracking_region(specification);
          if (coefficient_tracking_region == nullptr ||
              !coefficient_tracking_region->is_full_domain)
            report.add(
              DiagnosticCategory::lowerability,
              specification.id,
              "coefficient_identification_full_domain_tracking",
              "The first coefficient-identification target supports full-domain tracking only.");
        }

      if (dirichlet_control_lifting)
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
              control->role != semantic::v1::VariableRole::control ||
              space == specification.spaces.end() ||
              space->topology != semantic::v1::SpaceTopology::h1 ||
              controlled_boundary == nullptr ||
              space->region_id != controlled_boundary->id)
            report.add(
              DiagnosticCategory::lowerability,
              specification.formulation.control_variable_id,
              "dirichlet_nodal_trace_control_space",
              "Select the registered continuous nodal trace control space on the controlled Dirichlet boundary.");
          if (!specification.formulation.constraint_id.empty())
            report.add(
              DiagnosticCategory::lowerability,
              specification.formulation.constraint_id,
              "dirichlet_control_box_constraint",
              "The first nodal Dirichlet lifting has no box-constraint realization.");
          if (uses_fixed_dirichlet_reconstruction(specification) ||
              uses_neumann_control(specification) || h1_control_regularisation ||
              coefficient_identification)
            report.add(
              DiagnosticCategory::lowerability,
              specification.id,
              "dirichlet_control_registered_combination",
              "The first Dirichlet-control target supports only diffusion-reaction, volume forcing, full-volume tracking, and L2 trace regularisation.");
          const auto tracking_region = selected_tracking_region(specification);
          if (tracking_region == nullptr || !tracking_region->is_full_domain)
            report.add(
              DiagnosticCategory::lowerability,
              specification.id,
              "dirichlet_control_full_domain_tracking",
              "The first Dirichlet-control target supports full-domain state tracking only.");
        }

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
      const bool dirichlet_control =
        uses_dirichlet_control_lifting(specification);
      const bool coefficient_identification =
        uses_parameter_diffusion_residual(specification);
      const bool complete_residual = coefficient_identification
        ? count_terms(ResidualTermKind::parameter_diffusion_reaction) == 1 &&
            count_terms(ResidualTermKind::volume_source) == 1 &&
            count_terms(ResidualTermKind::diffusion_reaction) == 0 &&
            count_terms(ResidualTermKind::volume_control) == 0 &&
            count_terms(ResidualTermKind::neumann_control) == 0
        : boundary_control
        ? count_terms(ResidualTermKind::diffusion_reaction) == 1 &&
            count_terms(ResidualTermKind::volume_source) == 1 &&
            count_terms(ResidualTermKind::neumann_control) == 1 &&
            count_terms(ResidualTermKind::volume_control) == 0
        : dirichlet_control
        ? count_terms(ResidualTermKind::diffusion_reaction) == 1 &&
            count_terms(ResidualTermKind::volume_source) == 1 &&
            count_terms(ResidualTermKind::volume_control) == 0 &&
            count_terms(ResidualTermKind::neumann_control) == 0
        : count_terms(ResidualTermKind::diffusion_reaction) == 1 &&
            count_terms(ResidualTermKind::volume_source) == 1 &&
            count_terms(ResidualTermKind::volume_control) == 1 &&
            count_terms(ResidualTermKind::neumann_control) == 0;
      if (!complete_residual)
        report.add(DiagnosticCategory::lowerability,
                   specification.id,
                   coefficient_identification
                     ? "complete_parameter_diffusion_residual_term_set"
                     : boundary_control
                         ? "complete_neumann_boundary_residual_term_set"
                         : dirichlet_control
                             ? "complete_dirichlet_control_residual_term_set"
                         : "complete_volume_residual_term_set",
                   coefficient_identification
                     ? "Declare exactly one parameter diffusion-reaction and one volume-source term."
                     : boundary_control
                     ? "Declare exactly one diffusion-reaction, volume-source, and Neumann-control term."
                     : dirichlet_control
                     ? "Declare exactly one diffusion-reaction and one volume-source term; the control enters through the declared lifting."
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
          count_data(DataRole::reaction) != 1 ||
          count_data(DataRole::regularisation_weight) != 1 ||
          (coefficient_identification ? count_data(DataRole::diffusion) != 0
                                      : count_data(DataRole::diffusion) != 1))
        report.add(DiagnosticCategory::lowerability,
                   specification.id,
                   coefficient_identification
                     ? "complete_parameter_data_set"
                     : "complete_volume_data_set",
                   coefficient_identification
                     ? "Declare one forcing, target, reaction, and parameter-regularisation datum, with no constant diffusion datum."
                     : "Declare one forcing, target, diffusion, reaction, and regularisation datum.");

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
      const bool complete_control_loss = coefficient_identification
        ? count_losses(LossKind::quadratic_parameter_regularisation) == 1 &&
            count_losses(LossKind::quadratic_control_regularisation) == 0 &&
            count_losses(LossKind::quadratic_h1_control_regularisation) == 0
        : h1_control_regularisation
        ? count_losses(LossKind::quadratic_h1_control_regularisation) == 1 &&
            count_losses(LossKind::quadratic_control_regularisation) == 0
        : count_losses(LossKind::quadratic_control_regularisation) == 1 &&
            count_losses(LossKind::quadratic_h1_control_regularisation) == 0 &&
            count_losses(LossKind::quadratic_parameter_regularisation) == 0;
      if (count_losses(LossKind::quadratic_tracking) != 1 ||
          !complete_control_loss || specification.losses.size() != 2)
        report.add(DiagnosticCategory::lowerability,
                   specification.id,
                   "complete_registered_loss_set",
                   coefficient_identification
                     ? "Declare exactly one tracking and one parameter-regularisation loss."
                     : h1_control_regularisation
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
      if (coefficient_identification &&
          selected_metric != specification.metrics.end() &&
          selected_metric->kind != semantic::v1::MetricKind::l2)
        report.add(DiagnosticCategory::lowerability,
                   selected_metric->id,
                   "parameter_l2_metric",
                   "Select the registered cellwise L2 metric for the coefficient parameter.");
      if (coefficient_identification)
        {
          const auto parameter = find_variable(
            specification, specification.formulation.control_variable_id);
          const auto parameter_observation = std::find_if(
            specification.observations.begin(),
            specification.observations.end(),
            [parameter](const semantic::v1::ObservationSpec &observation) {
              return parameter != nullptr &&
                     observation.kind ==
                       semantic::v1::ObservationKind::volume_restriction &&
                     observation.input_variable_id == parameter->id;
            });
          const auto parameter_loss = std::find_if(
            specification.losses.begin(),
            specification.losses.end(),
            [](const semantic::v1::LossSpec &loss) {
              return loss.kind ==
                     semantic::v1::LossKind::quadratic_parameter_regularisation;
            });
          if (parameter == nullptr ||
              selected_metric == specification.metrics.end() ||
              selected_metric->variable_id != parameter->id ||
              parameter_observation == specification.observations.end() ||
              parameter_loss == specification.losses.end() ||
              parameter_loss->source_observation_id != parameter_observation->id)
            report.add(
              DiagnosticCategory::lowerability,
              specification.id,
              "registered_parameter_observation_metric_loss",
              "Connect the registered parameter observation, L2 metric, and parameter regularisation to the parameter decision variable.");
        }

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
      else if (dirichlet_control)
        {
          const auto state_observation = std::find_if(
            specification.observations.begin(),
            specification.observations.end(),
            [state](const semantic::v1::ObservationSpec &observation) {
              return state != nullptr &&
                     observation.kind ==
                       semantic::v1::ObservationKind::volume_restriction &&
                     observation.input_variable_id == state->id;
            });
          const auto control_restriction = std::find_if(
            specification.observations.begin(),
            specification.observations.end(),
            [](const semantic::v1::ObservationSpec &observation) {
              return observation.kind ==
                     semantic::v1::ObservationKind::boundary_restriction;
            });
          const auto control_region =
            selected_dirichlet_control_region(specification);
          if (specification.observations.size() != 2 ||
              state_observation == specification.observations.end() ||
              control_restriction == specification.observations.end())
            report.add(
              DiagnosticCategory::lowerability,
              specification.id,
              "complete_dirichlet_control_observation_set",
              "Declare one full-volume physical-state restriction and one controlled-boundary control restriction.");
          if (state_observation != specification.observations.end())
            {
              const auto region = find_region(specification,
                                              state_observation->region_id);
              if (region == nullptr || !region->is_full_domain)
                report.add(
                  DiagnosticCategory::lowerability,
                  state_observation->id,
                  "dirichlet_control_full_volume_observation",
                  "Track the physical state on the full volume in the first Dirichlet-control target.");
            }
          if (control_restriction != specification.observations.end())
            {
              const auto region = find_region(specification,
                                              control_restriction->region_id);
              if (region == nullptr || control_region == nullptr ||
                  region->id != control_region->id)
                report.add(
                  DiagnosticCategory::lowerability,
                  control_restriction->id,
                  "dirichlet_control_observation_region",
                  "Restrict the nodal trace control on its declared controlled boundary.");
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
                   "one_state_one_decision_one_equation",
                   "The executable DTO contract currently supports one state, one control-or-parameter decision, and one equation block.");
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

    static std::string
    describe(const ConstraintRealisation realisation)
    {
      switch (realisation)
        {
          case ConstraintRealisation::none:
            return "none";
          case ConstraintRealisation::cellwise_l2:
            return "FE_DGQ(0) coefficientwise l2_cellwise clipping";
          case ConstraintRealisation::cellwise_parameter_l2:
            return "FE_DGQ(0) coefficientwise l2_cellwise_parameter clipping";
          case ConstraintRealisation::facewise_l2:
            return "facewise-constant coefficientwise l2_facewise clipping";
        }
      contract::require(false, "Unknown compiled constraint realization");
      return {};
    }

    static std::string
    constraint_realisation_id(const ConstraintRealisation realisation)
    {
      switch (realisation)
        {
          case ConstraintRealisation::none:
            return "none";
          case ConstraintRealisation::cellwise_l2:
            return "l2_cellwise";
          case ConstraintRealisation::cellwise_parameter_l2:
            return "l2_cellwise_parameter";
          case ConstraintRealisation::facewise_l2:
            return "l2_facewise";
        }
      contract::require(false, "Unknown compiled constraint realization");
      return {};
    }

    static bool
    uses_neumann_target(const CompiledTargetKind target)
    {
      return target == CompiledTargetKind::neumann_boundary ||
             target == CompiledTargetKind::pure_neumann;
    }

    static bool
    uses_h1_control_target(const CompiledTargetKind target)
    {
      return target == CompiledTargetKind::h1_control_l2_metric ||
             target == CompiledTargetKind::h1_control_h1_metric;
    }

    static std::string
    control_space_description(const CompiledTargetKind       target,
                              const DealiiDiscretisationPolicy &policy)
    {
      switch (target)
        {
          case CompiledTargetKind::dirichlet_control:
            return "one shared nodal trace coefficient per state DoF on the complete controlled exterior boundary";
          case CompiledTargetKind::neumann_boundary:
          case CompiledTargetKind::pure_neumann:
            return "one facewise-constant coefficient per marked state boundary face";
          case CompiledTargetKind::coefficient_identification:
            return "cellwise-constant positive diffusion parameter FE_DGQ(0) on the state mesh";
          case CompiledTargetKind::h1_control_l2_metric:
          case CompiledTargetKind::h1_control_h1_metric:
            return "continuous scalar FE_Q(" +
                   std::to_string(policy.state_degree) + ") on the state mesh";
          case CompiledTargetKind::direct_volume:
          case CompiledTargetKind::assembled_volume:
            return "FE_DGQ(0) on the state active-cell mesh";
        }
      contract::require(false, "Unknown compiled target kind");
      return {};
    }

    static unsigned int
    resolved_maximum_iterations(
      const dealii_backend::SPDLinearSolvePolicy &policy,
      const std::size_t                           dimension)
    {
      if (policy.maximum_iterations != 0)
        return policy.maximum_iterations;
      contract::require(
        dimension <= std::numeric_limits<unsigned int>::max() / 10U,
        "Compiled solve dimension exceeds the iteration-policy range");
      return std::max(100U, 10U * static_cast<unsigned int>(dimension));
    }

    static CompiledSolvePolicyRecord
    spd_solve_record(const dealii_backend::SPDLinearSolvePolicy &policy,
                     const std::size_t                           dimension,
                     std::string                                 nullspace_policy)
    {
      return {LinearSolveAlgorithm::serial_cg,
              "identity",
              resolved_maximum_iterations(policy, dimension),
              policy.relative_tolerance,
              policy.absolute_tolerance,
              std::move(nullspace_policy)};
    }

    static std::size_t
    compiled_space_dimension(
      const semantic::v1::ProblemSpec &specification,
      const semantic::v1::SpaceSpec &  space,
      const contract::ExecutableModelT<dealii_backend::SerialBackend> &executable)
    {
      const auto state = find_variable(
        specification, specification.formulation.state_variable_id);
      const auto decision = find_variable(
        specification, specification.formulation.control_variable_id);
      if (state != nullptr && state->space_id == space.id)
        return executable.variable_layout()->dimension(0);
      if (decision != nullptr && decision->space_id == space.id)
        return executable.variable_layout()->dimension(1);
      if (!specification.equations.empty() &&
          specification.equations.front().test_space_id == space.id)
        return executable.test_layout()->dimension(0);
      const auto observation = std::find_if(
        specification.observations.begin(),
        specification.observations.end(),
        [&space](const semantic::v1::ObservationSpec &candidate) {
          return candidate.output_space_id == space.id;
        });
      if (observation != specification.observations.end())
        return decision != nullptr &&
                   observation->input_variable_id == decision->id
                 ? executable.variable_layout()->dimension(1)
                 : executable.variable_layout()->dimension(0);
      return 0;
    }

    static std::string
    compiled_space_runtime_role(const semantic::v1::SpaceRole role)
    {
      switch (role)
        {
          case semantic::v1::SpaceRole::state:
            return "state";
          case semantic::v1::SpaceRole::test:
            return "test_and_adjoint";
          case semantic::v1::SpaceRole::control:
            return "decision_control";
          case semantic::v1::SpaceRole::parameter:
            return "decision_parameter";
          case semantic::v1::SpaceRole::observation:
            return "observation";
          case semantic::v1::SpaceRole::data:
            return "data";
          case semantic::v1::SpaceRole::unspecified:
            return "unspecified";
        }
      return "unspecified";
    }

    static std::string
    compiled_space_finite_element(
      const semantic::v1::SpaceSpec &   space,
      const CompiledTargetKind          target,
      const DealiiDiscretisationPolicy &policy)
    {
      switch (space.role)
        {
          case semantic::v1::SpaceRole::state:
          case semantic::v1::SpaceRole::test:
            return "scalar FE_Q(" + std::to_string(policy.state_degree) + ")";
          case semantic::v1::SpaceRole::control:
          case semantic::v1::SpaceRole::parameter:
            return control_space_description(target, policy);
          case semantic::v1::SpaceRole::observation:
            return "lowered observation coefficients";
          case semantic::v1::SpaceRole::data:
            return "external binding";
          case semantic::v1::SpaceRole::unspecified:
            return "unspecified";
        }
      return "unspecified";
    }

    static std::string
    bound_binding_description(
      const std::optional<CellwiseBoxDataBindings> &bounds,
      const std::optional<FacewiseBoxDataBindings> &facewise_bounds)
    {
      if (bounds)
        return std::holds_alternative<double>(bounds->lower)
                 ? "scalar bound"
                 : "exact-layout cellwise coefficient vector";
      if (facewise_bounds)
        return std::holds_alternative<double>(facewise_bounds->lower)
                 ? "scalar bound"
                 : "exact-layout facewise coefficient vector";
      return "unbound";
    }

    template <int dim>
    static CompilationManifest
    make_manifest(
      const semantic::v1::ProblemSpec &specification,
      const DealiiDiscretisationPolicy &policy,
      const ConstraintRealisation       constraint_realisation,
      const CompiledTargetKind          target,
      const semantic::v1::RegionSpec &  tracking_region,
      const semantic::v1::RegionSpec *  control_boundary_region,
      const DealiiDataBindings<dim> &    data,
      const std::optional<CellwiseBoxDataBindings> &bounds,
      const std::optional<FacewiseBoxDataBindings> &facewise_bounds,
      const dealii::Triangulation<dim> & triangulation,
      const std::string &                mesh_provenance,
      const bool                         owns_mesh,
      const contract::ExecutableModelT<dealii_backend::SerialBackend> &executable,
      const dealii_backend::MassMetric & metric,
      const ScalarLoweringPlan *          scalar_plan)
    {
      CompilationManifest manifest;
      const bool uses_fixed_reconstruction =
        uses_fixed_dirichlet_reconstruction(specification);
      const bool uses_dirichlet_control_lifting =
        target == CompiledTargetKind::dirichlet_control;
      const bool uses_assembled_v1_target =
        target == CompiledTargetKind::assembled_volume;
      const bool uses_neumann_boundary_control = uses_neumann_target(target);
      const bool uses_mean_zero_gauge =
        target == CompiledTargetKind::pure_neumann;
      const bool uses_h1_control_regularisation = uses_h1_control_target(target);
      const bool uses_h1_control_metric =
        target == CompiledTargetKind::h1_control_h1_metric;
      const bool uses_coefficient_identification =
        target == CompiledTargetKind::coefficient_identification;

      manifest.formulation_record = {
        specification.formulation.id,
        specification.formulation.kind,
        ExecutionRealisation::assembled,
        "tested dual coefficients with dot pairing"};
      manifest.mesh_record = {
        dim,
        triangulation.n_active_cells(),
        mesh_provenance,
        owns_mesh ? MeshLifetimePolicy::owned_session
                  : MeshLifetimePolicy::borrowed_immutable};
      for (const auto &space : specification.spaces)
        manifest.spaces.push_back(
          {space.id,
           space.role,
           compiled_space_runtime_role(space.role),
           space.region_id,
           compiled_space_finite_element(space, target, policy),
           compiled_space_dimension(specification, space, executable)});
      for (const auto &binding : specification.data)
        {
          CompiledBindingRecord record;
          record.semantic_id = binding.id;
          record.role = binding.role;
          switch (binding.role)
            {
              case semantic::v1::DataRole::forcing:
                record.representation = "analytic Function at quadrature";
                record.provenance = data.provenance.forcing;
                break;
              case semantic::v1::DataRole::desired_state:
                record.representation = "analytic Function at quadrature";
                record.provenance = data.provenance.desired_state;
                break;
              case semantic::v1::DataRole::fixed_dirichlet_lifting:
                record.representation = "Function interpolated at boundary DoFs";
                record.provenance = data.provenance.fixed_dirichlet_data;
                break;
              case semantic::v1::DataRole::diffusion:
                record.representation = uses_coefficient_identification
                                          ? "decision parameter"
                                          : "scalar constant";
                record.provenance = uses_coefficient_identification
                                      ? specification.formulation.control_variable_id
                                      : std::to_string(*data.diffusion);
                break;
              case semantic::v1::DataRole::reaction:
                record.representation = "scalar constant";
                record.provenance = std::to_string(data.reaction);
                break;
              case semantic::v1::DataRole::conservative_transport:
              case semantic::v1::DataRole::advective_transport:
              case semantic::v1::DataRole::robin_coefficient:
              case semantic::v1::DataRole::robin_source:
                record.representation = "not lowered by this target";
                record.provenance = "not bound";
                break;
              case semantic::v1::DataRole::regularisation_weight:
                record.representation = "scalar constant";
                record.provenance = std::to_string(data.regularisation_weight);
                break;
              case semantic::v1::DataRole::lower_bound:
              case semantic::v1::DataRole::upper_bound:
                record.representation =
                  bound_binding_description(bounds, facewise_bounds);
                record.provenance = "caller-supplied compiled bound data";
                break;
              case semantic::v1::DataRole::unspecified:
                record.representation = "unspecified";
                record.provenance = "unspecified";
                break;
            }
          manifest.bindings.push_back(std::move(record));
        }

      const std::size_t state_dimension =
        executable.test_layout()->dimension(0);
      if (uses_mean_zero_gauge)
        {
          const CompiledSolvePolicyRecord direct_record{
            LinearSolveAlgorithm::serial_sparse_direct_umfpack,
            "not applicable",
            1,
            0.0,
            0.0,
            "one mean-zero Lagrange multiplier"};
          manifest.state_solve_record = direct_record;
          manifest.adjoint_solve_record = direct_record;
        }
      else
        {
          manifest.state_solve_record = spd_solve_record(
            policy.state_solve, state_dimension, "fixed Dirichlet");
          manifest.adjoint_solve_record = spd_solve_record(
            policy.adjoint_solve, state_dimension, "fixed Dirichlet");
        }
      manifest.metric_record = {
        specification.formulation.metric_id,
        metric.id(),
        uses_h1_control_metric ? "mass plus stiffness Riesz map"
                               : "mass Riesz map",
        {LinearSolveAlgorithm::serial_cg,
         "identity",
         policy.control_metric_solve.maximum_iterations,
         policy.control_metric_solve.relative_tolerance,
         policy.control_metric_solve.absolute_tolerance,
         "not applicable"}};
      manifest.constraint_record = {
        constraint_realisation != ConstraintRealisation::none,
        specification.formulation.constraint_id,
        constraint_realisation_id(constraint_realisation),
        constraint_realisation == ConstraintRealisation::none ? "none"
                                                               : metric.id()};
      if (scalar_plan != nullptr)
        manifest.lowering_handler_records = scalar_plan->provenance;
      manifest.semantic_problem_id = specification.id;
      manifest.compiler_id =
        "nmopt.compiler.v1.dealii.scalar_diffusion_reaction";
      manifest.backend = "deal.II serial Vector<double>";
      manifest.execution = "assembled";
      manifest.state_space = "scalar FE_Q(" +
                             std::to_string(policy.state_degree) + ")";
      manifest.control_space = control_space_description(target, policy);
      manifest.quadrature = "QGauss(" +
                            std::to_string(policy.state_degree + 2) + ")";
      manifest.dual_representation = "tested dual coefficients with dot pairing";
      manifest.data_rule = uses_neumann_boundary_control
        ? "analytic desired-state Function at selected QGauss(" +
            std::to_string(policy.state_degree + 2) +
            ") boundary face quadrature; scalar coefficients and forcing Function at volume quadrature"
        : "analytic desired-state Function at selected QGauss(" +
            std::to_string(policy.state_degree + 2) + ") volume quadrature" +
            (uses_coefficient_identification
              ? "; forcing Function, reaction and regularisation scalars; diffusion is the parameter decision block"
              : uses_fixed_reconstruction
              ? "; fixed Dirichlet Function interpolated at boundary DoFs"
              : uses_dirichlet_control_lifting
              ? "; scalar coefficients and forcing Function at volume quadrature; Dirichlet trace is the decision block"
              : "; scalar coefficients and forcing Function at volume quadrature");
      manifest.observation_realisation = uses_neumann_boundary_control
        ? boundary_observation_realisation(tracking_region)
        : observation_realisation(tracking_region);
      manifest.metric_solve_policy =
        std::string("serial CG for ") +
        (uses_h1_control_metric ? "h1_continuous Riesz map"
                                : uses_dirichlet_control_lifting
                                    ? "l2_dirichlet_trace Riesz map"
                                : uses_coefficient_identification
                                    ? "l2_cellwise_parameter Riesz map"
                                    : "L2 Riesz map") +
        ": maximum iterations=" +
        std::to_string(policy.control_metric_solve.maximum_iterations) +
        ", relative tolerance=" +
        std::to_string(policy.control_metric_solve.relative_tolerance) +
        ", absolute tolerance=" +
        std::to_string(policy.control_metric_solve.absolute_tolerance);
      manifest.constraint_realisation =
        describe(constraint_realisation);
      manifest.lifting_realisation = uses_mean_zero_gauge
                                       ? "none; pure-Neumann state uses an explicit mean-zero gauge"
                                       : uses_fixed_reconstruction
                                       ? "y_phys = P_h y_hat + ell_0,h; independent FE_Q coordinates, AffineConstraints reconstruction, and P_h^* pullbacks"
                                       : uses_dirichlet_control_lifting
                                           ? "y_phys = P_h y_hat + L_D,h u_h; complete-boundary shared nodal trace lifting, independent FE_Q coordinates, and P_h^*/L_D,h^* pullbacks"
                                       : uses_assembled_v1_target
                                           ? "y_phys = P_h y_hat; independent FE_Q coordinates and AffineConstraints reconstruction"
                                           : "homogeneous full-vector Dirichlet rows; no inhomogeneous lifting";
      manifest.nullspace_policy = uses_mean_zero_gauge
        ? "one mean-zero Lagrange multiplier; discrete constant-mode compatibility is enforced for forcing and every control state load"
        : "not applicable: non-empty fixed Dirichlet boundary";
      manifest.state_adjoint_solve_policy = uses_mean_zero_gauge
        ? "serial SparseDirectUMFPACK on the augmented symmetric state and adjoint saddle systems"
        : uses_coefficient_identification
          ? "serial CG with identity preconditioner; parameter-dependent SPD state matrix is reassembled for each state and adjoint solve"
        : "serial CG with identity preconditioner for symmetric positive-definite operator";
      manifest.provenance = "DTO";
      if (uses_neumann_boundary_control && control_boundary_region != nullptr)
        manifest.declared_assumptions.push_back(
          "neumann_control_realisation: facewise-constant FEFaceValues pairing on boundary ids " +
          boundary_id_list(*control_boundary_region));
      if (uses_dirichlet_control_lifting && control_boundary_region != nullptr)
        manifest.declared_assumptions.push_back(
          "dirichlet_control_lifting: complete-exterior-boundary shared nodal trace map on boundary ids " +
          boundary_id_list(*control_boundary_region) +
          "; no corner/interface averaging, hanging-node relation, or box policy is registered");
      if (uses_h1_control_regularisation)
        manifest.declared_assumptions.push_back(
          "h1_control_regularisation: alpha/2 u^T (M_u + K_u) u; search metric=" +
          std::string(uses_h1_control_metric ? "h1_continuous" : "l2_continuous"));
      if (uses_coefficient_identification)
        manifest.declared_assumptions.push_back(
          "coefficient_identification: positive cellwise physical diffusion parameter; A(m) is reassembled for every state and adjoint solve");
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

    DealiiCapabilityRegistryV1 capabilities_;
    DealiiScalarLoweringPlanner scalar_planner_;
  };
} // namespace nmopt::compiler::v1
