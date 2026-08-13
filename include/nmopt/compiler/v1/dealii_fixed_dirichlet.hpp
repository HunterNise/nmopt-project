#pragma once

#include "nmopt/compiler/v1/dealii_scalar_plan.hpp"
#include "nmopt/compiler/v1/dealii_types.hpp"
#include "nmopt/contract/executable_model.hpp"
#include "nmopt/dealii/cellwise_box_constraint.hpp"
#include "nmopt/dealii/mass_metric.hpp"
#include "nmopt/dealii/serial_backend.hpp"
#include "nmopt/dealii/serial_spd_solver.hpp"

#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/numerics/vector_tools.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace nmopt::compiler::v1::detail
{
  // This bounded scalar component target represents the independent state coefficients
  // y_hat. It evaluates the compiled equation and observation/loss on
  //
  //   y_phys = P_h y_hat + ell_0,h,
  //
  // and pulls covectors back with P_h^*. The direct v0 lowerer is not used or
  // modified here: it remains the homogeneous comparison implementation.
  template <int dim>
  class ScalarComponentModel final
    : public contract::ExecutableModelT<dealii_backend::SerialBackend>
  {
  public:
    using Backend = dealii_backend::SerialBackend;
    using Vector = dealii::Vector<double>;
    using Primal = contract::PrimalBlockT<Backend>;
    using Covector = contract::CovectorBlockT<Backend>;
    using SolveResult = contract::FormulationSolveResultT<Backend>;

    ScalarComponentModel(
      dealii::Triangulation<dim> &triangulation,
      const dealii::Function<dim> &forcing,
      const dealii::Function<dim> &desired_state,
      std::optional<std::reference_wrapper<const dealii::Function<dim>>>
        fixed_dirichlet_data,
      const double              diffusion,
      const double              reaction,
      const double              regularisation_weight,
      const unsigned int        state_degree,
      const ScalarLoweringPlan &plan)
      : ScalarComponentModel(
          triangulation,
          forcing,
          desired_state,
          fixed_dirichlet_data,
          diffusion,
          reaction,
          regularisation_weight,
          state_degree,
          boundary_ids_from_plan(plan),
          material_ids_from_plan(plan),
          has_observation_operator(
            plan, ScalarObservationOperatorKind::h1_state_restriction),
          normal_flux_boundary_ids_from_plan(plan),
          point_sensor_coordinates_from_plan(plan),
          {},
          nullptr)
    {
      contract::require(
        has_residual_operator(plan,
                              ScalarResidualOperatorKind::diffusion_reaction) &&
          has_residual_operator(plan,
                                ScalarResidualOperatorKind::volume_source) &&
          has_residual_operator(plan,
                                ScalarResidualOperatorKind::volume_control),
        "The composed scalar target needs diffusion-reaction, source, and volume-control contributions");
      contract::require(
        has_loss_operator(plan, ScalarLossOperatorKind::quadratic_tracking) &&
          has_loss_operator(
            plan, ScalarLossOperatorKind::quadratic_control_regularisation),
        "The composed scalar target needs tracking and control-regularisation losses");
      contract::require(
        (plan.transformation ==
           ScalarTransformationOperatorKind::fixed_dirichlet_reconstruction) ==
          fixed_dirichlet_data.has_value(),
        "The scalar reconstruction plan and fixed-data binding disagree");
    }

    ScalarComponentModel(
      dealii::Triangulation<dim> &             triangulation,
      const dealii::Function<dim> &            forcing,
      const dealii::Function<dim> &            desired_state,
      const DealiiGeneralScalarDataBindings<dim> &general_scalar_data,
      const double                              regularisation_weight,
      const unsigned int                        state_degree,
      const ScalarLoweringPlan &                plan)
      : ScalarComponentModel(triangulation,
                             forcing,
                             desired_state,
                             std::nullopt,
                             0.0,
                             0.0,
                             regularisation_weight,
                             state_degree,
                             boundary_ids_from_plan(plan),
                             material_ids_from_plan(plan),
                             has_observation_operator(
                               plan,
                               ScalarObservationOperatorKind::h1_state_restriction),
                             {},
                             {},
                             robin_boundary_ids_from_plan(plan),
                             &general_scalar_data,
                             require_general_scalar_data_placements(plan))
    {
      contract::require(
        has_residual_operator(plan,
                              ScalarResidualOperatorKind::tensor_diffusion) &&
          has_residual_operator(
            plan, ScalarResidualOperatorKind::conservative_transport) &&
          has_residual_operator(plan,
                                ScalarResidualOperatorKind::advective_transport) &&
          has_residual_operator(plan, ScalarResidualOperatorKind::reaction) &&
          has_residual_operator(plan,
                                ScalarResidualOperatorKind::volume_source) &&
          has_residual_operator(plan,
                                ScalarResidualOperatorKind::volume_control) &&
          has_residual_operator(plan,
                                ScalarResidualOperatorKind::robin_bilinear) &&
          has_residual_operator(plan, ScalarResidualOperatorKind::robin_source),
        "The general scalar target needs every selected volume and Robin contribution");
      contract::require(
        has_loss_operator(plan, ScalarLossOperatorKind::quadratic_tracking) &&
          has_loss_operator(
            plan, ScalarLossOperatorKind::quadratic_control_regularisation),
        "The general scalar target needs tracking and control-regularisation losses");
      contract::require(
        plan.transformation == ScalarTransformationOperatorKind::none,
        "The first general scalar target supports homogeneous fixed Dirichlet data only");
    }

    ScalarComponentModel(
      dealii::Triangulation<dim> &                triangulation,
      const dealii::Function<dim> &               forcing,
      const dealii::Function<dim> &               desired_state,
      std::optional<std::reference_wrapper<const dealii::Function<dim>>>
                                                    fixed_dirichlet_data,
      const double                                  diffusion,
      const double                                  reaction,
      const double                                  regularisation_weight,
      const unsigned int                            state_degree,
      std::set<dealii::types::boundary_id>          dirichlet_boundary_ids,
      std::set<dealii::types::material_id>          observation_material_ids,
      const bool                                    uses_h1_state_observation,
      std::set<dealii::types::boundary_id>          normal_flux_boundary_ids,
      std::vector<std::vector<double>>              point_sensor_coordinates,
      std::set<dealii::types::boundary_id>          robin_boundary_ids,
      const DealiiGeneralScalarDataBindings<dim> *  general_scalar_data,
      const bool                                    general_scalar_plan_validated = true)
      : state_fe_(state_degree)
      , control_fe_(0)
      , state_dof_handler_(triangulation)
      , control_dof_handler_(triangulation)
      , diffusion_(diffusion)
      , reaction_(reaction)
      , regularisation_weight_(regularisation_weight)
      , uses_general_scalar_(general_scalar_data != nullptr)
      , dirichlet_boundary_ids_(std::move(dirichlet_boundary_ids))
      , observation_material_ids_(std::move(observation_material_ids))
      , uses_h1_state_observation_(uses_h1_state_observation)
      , normal_flux_boundary_ids_(std::move(normal_flux_boundary_ids))
      , uses_normal_flux_(!normal_flux_boundary_ids_.empty())
      , point_sensor_coordinates_(std::move(point_sensor_coordinates))
      , uses_point_sensor_(!point_sensor_coordinates_.empty())
      , robin_boundary_ids_(std::move(robin_boundary_ids))
    {
      contract::require(uses_general_scalar_ || diffusion_ > 0.0,
                        "Diffusion coefficient must be strictly positive");
      contract::require(uses_general_scalar_ || reaction_ >= 0.0,
                        "Reaction coefficient must be non-negative");
      contract::require(regularisation_weight_ > 0.0,
                        "Control regularisation weight must be strictly positive");
      contract::require(state_degree > 0,
                        "State FE degree must be at least one");
      contract::require(!dirichlet_boundary_ids_.empty(),
                        "The assembled v1 target needs a fixed Dirichlet boundary");
      contract::require(!uses_general_scalar_ || !robin_boundary_ids_.empty(),
                        "The general scalar target needs a Robin boundary");
      contract::require(general_scalar_plan_validated,
                        "The general scalar target needs resolved coefficient data placements");
      contract::require(!uses_h1_state_observation_ ||
                          observation_material_ids_.empty(),
                        "The registered H1 state observation is full-domain only");
      contract::require(!uses_point_sensor_ ||
                          !point_sensor_coordinates_.empty(),
                        "The point-sensor target needs immutable sensor coordinates");

      state_dof_handler_.distribute_dofs(state_fe_);
      control_dof_handler_.distribute_dofs(control_fe_);
      build_constraints(fixed_dirichlet_data);
      build_reconstruction();
      initialise_storage();
      assemble_physical_operators(forcing,
                                  desired_state,
                                  general_scalar_data);
      assemble_reduced_solve_operators();
      if (uses_general_scalar_)
        {
          nonsymmetric_solver_ =
            std::make_unique<dealii::SparseDirectUMFPACK>();
          nonsymmetric_solver_->initialize(reduced_system_matrix_);
        }
    }

    const contract::LayoutPtr &
    variable_layout() const override
    {
      return variable_layout_;
    }

    const contract::LayoutPtr &
    test_layout() const override
    {
      return test_layout_;
    }

    const contract::LayoutPtr &
    state_layout() const
    {
      return state_layout_;
    }

    const contract::LayoutPtr &
    control_layout() const
    {
      return control_layout_;
    }

    dealii_backend::MassMetric
    control_l2_metric(
      dealii_backend::MassMetricSolveParameters solve_parameters = {}) const
    {
      return dealii_backend::MassMetric("l2_cellwise",
                                        control_layout_,
                                        control_mass_,
                                        solve_parameters);
    }

    dealii_backend::CellwiseBoxConstraint
    control_l2_box_constraint(
      Vector                              lower,
      Vector                              upper,
      const dealii_backend::MassMetric & projection_metric) const
    {
      return dealii_backend::CellwiseBoxConstraint(control_layout_,
                                                    std::move(lower),
                                                    std::move(upper),
                                                    projection_metric);
    }

    dealii_backend::CellwiseBoxConstraint
    control_l2_box_constraint(
      const double                        lower,
      const double                        upper,
      const dealii_backend::MassMetric & projection_metric) const
    {
      return dealii_backend::CellwiseBoxConstraint(control_layout_,
                                                    lower,
                                                    upper,
                                                    projection_metric);
    }

    Vector
    reconstruct_physical_state(const Primal &independent_state) const
    {
      contract::require(independent_state.layout()->compatible_with(*state_layout_),
                        "Physical-state reconstruction needs the state layout");
      return reconstruct(independent_state.block(0));
    }

    Covector
    residual(const Primal &variables) const override
    {
      require_variables(variables, "Residual");
      Vector physical_value(state_dof_handler_.n_dofs());
      physical_system_matrix_.vmult(physical_value,
                                    reconstruct(variables.block(0)));
      physical_value.add(-1.0, forcing_load_);

      Vector control_contribution(state_dof_handler_.n_dofs());
      physical_control_coupling_.vmult(control_contribution,
                                        variables.block(1));
      physical_value.add(-1.0, control_contribution);
      return Covector(test_layout_, {pullback(physical_value)});
    }

    Covector
    residual_jvp(const Primal &variables,
                 const Primal &variable_tangent) const override
    {
      require_variables(variables, "Residual JVP");
      require_variables(variable_tangent, "Residual JVP tangent");

      Vector physical_value(state_dof_handler_.n_dofs());
      physical_system_matrix_.vmult(physical_value,
                                    embed_tangent(variable_tangent.block(0)));
      Vector control_contribution(state_dof_handler_.n_dofs());
      physical_control_coupling_.vmult(control_contribution,
                                        variable_tangent.block(1));
      physical_value.add(-1.0, control_contribution);
      return Covector(test_layout_, {pullback(physical_value)});
    }

    Covector
    residual_vjp(const Primal &variables,
                 const Primal &test_seed) const override
    {
      require_variables(variables, "Residual VJP");
      contract::require(test_seed.layout()->compatible_with(*test_layout_),
                        "Residual VJP seed has an incompatible test layout");

      const Vector physical_seed = embed_tangent(test_seed.block(0));
      Vector physical_state(state_dof_handler_.n_dofs());
      physical_system_matrix_.Tvmult(physical_state, physical_seed);
      Vector control(control_dof_handler_.n_dofs());
      physical_control_coupling_.Tvmult(control, physical_seed);
      control *= -1.0;
      return Covector(variable_layout_,
                      {pullback(physical_state), std::move(control)});
    }

    double
    objective(const Primal &variables) const override
    {
      require_variables(variables, "Objective");
      const Vector physical_state = reconstruct(variables.block(0));
      Vector state_tracking_times_state(state_dof_handler_.n_dofs());
      physical_state_tracking_operator_.vmult(state_tracking_times_state,
                                              physical_state);
      const double state_value =
        0.5 * (physical_state * state_tracking_times_state) -
        (desired_state_load_ * physical_state) + 0.5 * desired_state_norm_;

      Vector control_mass_times_control(control_dof_handler_.n_dofs());
      control_mass_->vmult(control_mass_times_control, variables.block(1));
      const double control_value = 0.5 * regularisation_weight_ *
                                   (variables.block(1) *
                                    control_mass_times_control);
      return state_value + control_value;
    }

    Covector
    objective_derivative(const Primal &variables) const override
    {
      require_variables(variables, "Objective derivative");
      Vector physical_state(state_dof_handler_.n_dofs());
      physical_state_tracking_operator_.vmult(
        physical_state, reconstruct(variables.block(0)));
      physical_state.add(-1.0, desired_state_load_);

      Vector control(control_dof_handler_.n_dofs());
      control_mass_->vmult(control, variables.block(1));
      control *= regularisation_weight_;
      return Covector(variable_layout_,
                      {pullback(physical_state), std::move(control)});
    }

    // Target-specific diagnostic ports used by the point-sensor contract.
    // They expose the selected finite-dimensional map without making sensor
    // coordinates part of the backend-neutral executable interface.
    std::vector<double>
    point_sensor_values(const Primal &variables) const
    {
      require_variables(variables, "Point-sensor value");
      contract::require(uses_point_sensor_,
                        "Point-sensor values need the point-sensor target");
      const Vector physical_state = reconstruct(variables.block(0));
      std::vector<double> values;
      values.reserve(point_sensor_evaluations_.size());
      for (const auto &evaluation : point_sensor_evaluations_)
        values.push_back(evaluation * physical_state);
      return values;
    }

    Covector
    point_sensor_jvp(const Primal &variable_tangent) const
    {
      require_variables(variable_tangent, "Point-sensor JVP tangent");
      contract::require(uses_point_sensor_,
                        "Point-sensor JVP needs the point-sensor target");
      const Vector physical_tangent = embed_tangent(
        variable_tangent.block(0));
      Vector values(point_sensor_evaluations_.size());
      for (std::size_t index = 0; index < point_sensor_evaluations_.size(); ++index)
        values[index] = point_sensor_evaluations_[index] * physical_tangent;
      const std::size_t observation_dimension = values.size();
      return Covector(std::make_shared<const contract::BlockLayout>(
                        "point_sensor_observation",
                        std::vector<contract::SpaceId>{{"point_sensor"}},
                        std::vector<std::size_t>{observation_dimension}),
                      {std::move(values)});
    }

    Covector
    point_sensor_vjp(const std::vector<double> &seed) const
    {
      contract::require(uses_point_sensor_,
                        "Point-sensor VJP needs the point-sensor target");
      contract::require(seed.size() == point_sensor_evaluations_.size(),
                        "Point-sensor VJP seed has the wrong dimension");
      Vector physical_covector(state_dof_handler_.n_dofs());
      for (std::size_t index = 0; index < seed.size(); ++index)
        physical_covector.add(seed[index], point_sensor_evaluations_[index]);
      return Covector(state_layout_, {pullback(physical_covector)});
    }

    // The normal-flux ports expose the face-quadrature map selected by the
    // semantic contract. Values and JVPs are pointwise normal derivatives;
    // the VJP includes the corresponding face-quadrature weights.
    std::vector<double>
    normal_flux_values(const Primal &variables) const
    {
      require_variables(variables, "Normal-flux value");
      contract::require(uses_normal_flux_,
                        "Normal-flux values need the normal-flux target");
      const Vector physical_state = reconstruct(variables.block(0));
      std::vector<double> values;
      values.reserve(normal_flux_evaluations_.size());
      for (const auto &evaluation : normal_flux_evaluations_)
        values.push_back(evaluation * physical_state);
      return values;
    }

    Covector
    normal_flux_jvp(const Primal &variable_tangent) const
    {
      require_variables(variable_tangent, "Normal-flux JVP tangent");
      contract::require(uses_normal_flux_,
                        "Normal-flux JVP needs the normal-flux target");
      const Vector physical_tangent = embed_tangent(
        variable_tangent.block(0));
      Vector values(normal_flux_evaluations_.size());
      for (std::size_t index = 0; index < normal_flux_evaluations_.size();
           ++index)
        values[index] = normal_flux_evaluations_[index] * physical_tangent;
      const std::size_t observation_dimension = values.size();
      return Covector(std::make_shared<const contract::BlockLayout>(
                        "normal_flux_observation",
                        std::vector<contract::SpaceId>{{"normal_flux_boundary"}},
                        std::vector<std::size_t>{observation_dimension}),
                      {std::move(values)});
    }

    Covector
    normal_flux_vjp(const std::vector<double> &seed) const
    {
      contract::require(uses_normal_flux_,
                        "Normal-flux VJP needs the normal-flux target");
      contract::require(seed.size() == normal_flux_evaluations_.size(),
                        "Normal-flux VJP seed has the wrong dimension");
      Vector physical_covector(state_dof_handler_.n_dofs());
      for (std::size_t index = 0; index < seed.size(); ++index)
        physical_covector.add(seed[index] * normal_flux_quadrature_weights_[index],
                              normal_flux_evaluations_[index]);
      return Covector(state_layout_, {pullback(physical_covector)});
    }

    const std::vector<double> &
    normal_flux_quadrature_weights() const
    {
      contract::require(uses_normal_flux_,
                        "Normal-flux quadrature weights need the normal-flux target");
      return normal_flux_quadrature_weights_;
    }

    Primal
    solve_state(const Primal &control) const
    {
      auto result = solve_state_with_report(control, {});
      contract::require(result.report.converged(),
                        "State solve did not converge under its declared policy");
      return std::move(result.solution);
    }

    SolveResult
    solve_state_with_report(
      const Primal &                              control,
      const dealii_backend::SPDLinearSolvePolicy &policy) const
    {
      contract::require(control.layout()->compatible_with(*control_layout_),
                        "State solve control has an incompatible layout");
      Vector right_hand_side = reduced_forcing_load_;
      Vector control_contribution(independent_state_dofs_.size());
      reduced_control_coupling_.vmult(control_contribution, control.block(0));
      right_hand_side.add(1.0, control_contribution);

      Vector state(independent_state_dofs_.size());
      if (uses_general_scalar_)
        {
          (void)policy;
          nonsymmetric_solver_->vmult(state, right_hand_side);
          return {Primal(state_layout_, {std::move(state)}),
                  dealii_backend::direct_solve_report(
                    "serial_sparse_direct_umfpack")};
        }
      auto report = solve_symmetric_system(reduced_system_matrix_,
                                           state,
                                           right_hand_side,
                                           policy);
      return {Primal(state_layout_, {std::move(state)}), std::move(report)};
    }

    Primal
    solve_adjoint(const Primal &full_point,
                  const Covector &state_objective_derivative) const
    {
      auto result = solve_adjoint_with_report(full_point,
                                              state_objective_derivative,
                                              {});
      contract::require(result.report.converged(),
                        "Adjoint solve did not converge under its declared policy");
      return std::move(result.solution);
    }

    SolveResult
    solve_adjoint_with_report(
      const Primal &                              full_point,
      const Covector &                            state_objective_derivative,
      const dealii_backend::SPDLinearSolvePolicy &policy) const
    {
      require_variables(full_point, "Adjoint solve point");
      contract::require(
        state_objective_derivative.layout()->compatible_with(*state_layout_),
        "Adjoint solve right-hand side has an incompatible state layout");

      Vector adjoint(independent_state_dofs_.size());
      if (uses_general_scalar_)
        {
          (void)policy;
          nonsymmetric_solver_->Tvmult(adjoint,
                                       state_objective_derivative.block(0));
          return {Primal(test_layout_, {std::move(adjoint)}),
                  dealii_backend::direct_solve_report(
                    "serial_sparse_direct_umfpack_transpose")};
        }
      auto report = solve_symmetric_system(reduced_system_matrix_,
                                           adjoint,
                                           state_objective_derivative.block(0),
                                           policy);
      return {Primal(test_layout_, {std::move(adjoint)}), std::move(report)};
    }

  private:
    static bool
    require_general_scalar_data_placements(const ScalarLoweringPlan &plan)
    {
      const auto placement_for_role = [&plan](const semantic::v1::DataRole role) {
        return std::find_if(
          plan.data_placements.begin(),
          plan.data_placements.end(),
          [role](const ScalarDataPlacement &placement) {
            return placement.role == role;
          });
      };
      const auto valid = [&plan, &placement_for_role](
                           const semantic::v1::DataRole role,
                           const semantic::v1::DataKind kind,
                           const ScalarDataEvaluationKind evaluation) {
        const auto placement = placement_for_role(role);
        if (placement == plan.data_placements.end() ||
            placement->kind != kind || placement->space_id.empty() ||
            placement->region_id.empty() ||
            placement->evaluation != evaluation ||
            placement->handler_id.empty())
          return false;
        return true;
      };
      if (!valid(semantic::v1::DataRole::diffusion,
                 semantic::v1::DataKind::tensor_function,
                 ScalarDataEvaluationKind::volume_quadrature) ||
          !valid(semantic::v1::DataRole::conservative_transport,
                 semantic::v1::DataKind::vector_function,
                 ScalarDataEvaluationKind::volume_quadrature) ||
          !valid(semantic::v1::DataRole::advective_transport,
                 semantic::v1::DataKind::vector_function,
                 ScalarDataEvaluationKind::volume_quadrature) ||
          !valid(semantic::v1::DataRole::reaction,
                 semantic::v1::DataKind::function,
                 ScalarDataEvaluationKind::volume_quadrature) ||
          !valid(semantic::v1::DataRole::robin_coefficient,
                 semantic::v1::DataKind::function,
                 ScalarDataEvaluationKind::boundary_face_quadrature) ||
          !valid(semantic::v1::DataRole::robin_source,
                 semantic::v1::DataKind::function,
                 ScalarDataEvaluationKind::boundary_face_quadrature))
        return false;

      for (const auto role : {semantic::v1::DataRole::robin_coefficient,
                              semantic::v1::DataRole::robin_source})
        {
          const auto placement = placement_for_role(role);
          const auto term = std::find_if(
            plan.residual_terms.begin(),
            plan.residual_terms.end(),
            [role, placement](const ScalarResidualContribution &contribution) {
              return std::any_of(
                contribution.data_ids.begin(),
                contribution.data_ids.end(),
                [placement](const std::string &data_id) {
                  return data_id == placement->semantic_id;
                }) &&
                ((role == semantic::v1::DataRole::robin_coefficient &&
                  contribution.operator_kind ==
                    ScalarResidualOperatorKind::robin_bilinear) ||
                 (role == semantic::v1::DataRole::robin_source &&
                  contribution.operator_kind ==
                    ScalarResidualOperatorKind::robin_source));
            });
          if (term == plan.residual_terms.end() ||
              placement->region_id != term->region_id)
            return false;
        }
      return true;
    }

    static bool
    has_residual_operator(const ScalarLoweringPlan &        plan,
                          const ScalarResidualOperatorKind kind)
    {
      return std::any_of(
        plan.residual_terms.begin(),
        plan.residual_terms.end(),
        [kind](const ScalarResidualContribution &contribution) {
          return contribution.operator_kind == kind;
        });
    }

    static bool
    has_observation_operator(const ScalarLoweringPlan &             plan,
                             const ScalarObservationOperatorKind kind)
    {
      return std::any_of(
        plan.observations.begin(),
        plan.observations.end(),
        [kind](const ScalarObservationContribution &contribution) {
          return contribution.operator_kind == kind;
        });
    }

    static std::vector<std::vector<double>>
    point_sensor_coordinates_from_plan(const ScalarLoweringPlan &plan)
    {
      return std::any_of(
               plan.observations.begin(),
               plan.observations.end(),
               [](const ScalarObservationContribution &contribution) {
                 return contribution.operator_kind ==
                        ScalarObservationOperatorKind::point_sensor;
               })
               ? plan.point_sensor_coordinates
               : std::vector<std::vector<double>>{};
    }

    static std::set<dealii::types::boundary_id>
    normal_flux_boundary_ids_from_plan(const ScalarLoweringPlan &plan)
    {
      return std::any_of(
               plan.observations.begin(),
               plan.observations.end(),
               [](const ScalarObservationContribution &contribution) {
                 return contribution.operator_kind ==
                        ScalarObservationOperatorKind::normal_flux;
               })
               ? std::set<dealii::types::boundary_id>(
                   plan.normal_flux_boundary_ids.begin(),
                   plan.normal_flux_boundary_ids.end())
               : std::set<dealii::types::boundary_id>{};
    }

    static bool
    has_loss_operator(const ScalarLoweringPlan &    plan,
                      const ScalarLossOperatorKind kind)
    {
      return std::any_of(
        plan.losses.begin(),
        plan.losses.end(),
        [kind](const ScalarLossContribution &contribution) {
          return contribution.operator_kind == kind;
        });
    }

    static std::set<dealii::types::boundary_id>
    boundary_ids_from_plan(const ScalarLoweringPlan &plan)
    {
      std::set<dealii::types::boundary_id> result;
      for (const auto id : plan.dirichlet_boundary_ids)
        result.insert(static_cast<dealii::types::boundary_id>(id));
      return result;
    }

    static std::set<dealii::types::material_id>
    material_ids_from_plan(const ScalarLoweringPlan &plan)
    {
      std::set<dealii::types::material_id> result;
      for (const auto id : plan.tracking_material_ids)
        result.insert(static_cast<dealii::types::material_id>(id));
      return result;
    }

    static std::set<dealii::types::boundary_id>
    robin_boundary_ids_from_plan(const ScalarLoweringPlan &plan)
    {
      std::set<dealii::types::boundary_id> result;
      for (const auto id : plan.robin_boundary_ids)
        result.insert(static_cast<dealii::types::boundary_id>(id));
      return result;
    }

    void
    require_variables(const Primal &variables, const char *operation) const
    {
      contract::require(
        variables.layout()->compatible_with(*variable_layout_),
        std::string(operation) + " received an incompatible variable layout");
    }

    void
    build_constraints(
      const std::optional<std::reference_wrapper<const dealii::Function<dim>>>
        fixed_dirichlet_data)
    {
      homogeneous_constraints_.clear();
      physical_constraints_.clear();
      dealii::DoFTools::make_hanging_node_constraints(state_dof_handler_,
                                                       homogeneous_constraints_);
      dealii::DoFTools::make_hanging_node_constraints(state_dof_handler_,
                                                       physical_constraints_);

      dealii::Functions::ZeroFunction<dim> zero;
      const dealii::Function<dim> &physical_dirichlet_data =
        fixed_dirichlet_data ? fixed_dirichlet_data->get() : zero;
      for (const auto boundary_id : dirichlet_boundary_ids_)
        {
          dealii::VectorTools::interpolate_boundary_values(
            state_dof_handler_, boundary_id, zero, homogeneous_constraints_);
          dealii::VectorTools::interpolate_boundary_values(state_dof_handler_,
                                                            boundary_id,
                                                            physical_dirichlet_data,
                                                            physical_constraints_);
        }
      homogeneous_constraints_.close();
      physical_constraints_.close();
    }

    void
    build_reconstruction()
    {
      const auto physical_size = state_dof_handler_.n_dofs();
      for (dealii::types::global_dof_index index = 0; index < physical_size;
           ++index)
        if (!homogeneous_constraints_.is_constrained(index))
          independent_state_dofs_.push_back(index);
      contract::require(!independent_state_dofs_.empty(),
                        "State reconstruction needs an independent state DoF");

      dealii::DynamicSparsityPattern reconstruction_dsp(
        physical_size, independent_state_dofs_.size());
      for (std::size_t column = 0; column < independent_state_dofs_.size();
           ++column)
        {
          Vector basis(physical_size);
          basis[independent_state_dofs_[column]] = 1.0;
          homogeneous_constraints_.distribute(basis);
          for (dealii::types::global_dof_index row = 0; row < physical_size;
               ++row)
            if (basis[row] != 0.0)
              reconstruction_dsp.add(row, column);
        }
      reconstruction_sparsity_.copy_from(reconstruction_dsp);
      reconstruction_.reinit(reconstruction_sparsity_);
      for (std::size_t column = 0; column < independent_state_dofs_.size();
           ++column)
        {
          Vector basis(physical_size);
          basis[independent_state_dofs_[column]] = 1.0;
          homogeneous_constraints_.distribute(basis);
          for (dealii::types::global_dof_index row = 0; row < physical_size;
               ++row)
            if (basis[row] != 0.0)
              reconstruction_.set(row, column, basis[row]);
        }

      lifting_.reinit(physical_size);
      physical_constraints_.distribute(lifting_);
    }

    void
    initialise_storage()
    {
      const auto physical_state_size = state_dof_handler_.n_dofs();
      const auto control_size = control_dof_handler_.n_dofs();
      const auto independent_size = independent_state_dofs_.size();

      variable_layout_ = std::make_shared<const contract::BlockLayout>(
        "reconstructed_variables",
        std::vector<contract::SpaceId>{{"state"}, {"control"}},
        std::vector<std::size_t>{independent_size, control_size});
      test_layout_ = std::make_shared<const contract::BlockLayout>(
        "reconstructed_state_test",
        std::vector<contract::SpaceId>{{"state_test"}},
        std::vector<std::size_t>{independent_size});
      state_layout_ = variable_layout_->single_block(0, "state");
      control_layout_ = variable_layout_->single_block(1, "control");

      dealii::DynamicSparsityPattern physical_state_dsp(physical_state_size,
                                                        physical_state_size);
      dealii::DoFTools::make_sparsity_pattern(state_dof_handler_,
                                               physical_state_dsp);
      physical_state_sparsity_.copy_from(physical_state_dsp);
      physical_system_matrix_.reinit(physical_state_sparsity_);
      physical_state_tracking_operator_.reinit(physical_state_sparsity_);

      dealii::DynamicSparsityPattern physical_control_dsp(physical_state_size,
                                                          control_size);
      std::vector<dealii::types::global_dof_index> state_indices(
        state_fe_.dofs_per_cell);
      std::vector<dealii::types::global_dof_index> control_indices(
        control_fe_.dofs_per_cell);
      auto state_cell = state_dof_handler_.begin_active();
      auto control_cell = control_dof_handler_.begin_active();
      for (; state_cell != state_dof_handler_.end(); ++state_cell, ++control_cell)
        {
          contract::require(control_cell != control_dof_handler_.end(),
                            "State and control DoF handlers do not share cells");
          state_cell->get_dof_indices(state_indices);
          control_cell->get_dof_indices(control_indices);
          for (const auto state_index : state_indices)
            for (const auto control_index : control_indices)
              physical_control_dsp.add(state_index, control_index);
        }
      contract::require(control_cell == control_dof_handler_.end(),
                        "State and control DoF handlers have different cells");
      physical_control_sparsity_.copy_from(physical_control_dsp);
      physical_control_coupling_.reinit(physical_control_sparsity_);

      dealii::DynamicSparsityPattern control_mass_dsp(control_size, control_size);
      dealii::DoFTools::make_sparsity_pattern(control_dof_handler_, control_mass_dsp);
      control_mass_sparsity_.copy_from(control_mass_dsp);
      control_mass_ = std::make_shared<dealii::SparseMatrix<double>>();
      control_mass_->reinit(control_mass_sparsity_);

      forcing_load_.reinit(physical_state_size);
      desired_state_load_.reinit(physical_state_size);
      reduced_forcing_load_.reinit(independent_size);
    }

    void
    assemble_physical_operators(
      const dealii::Function<dim> &forcing,
      const dealii::Function<dim> &desired_state,
      const DealiiGeneralScalarDataBindings<dim> *general_scalar_data)
    {
      const unsigned int quadrature_order =
        std::max(state_fe_.degree, control_fe_.degree) + 2;
      const dealii::QGauss<dim> quadrature(quadrature_order);
      dealii::FEValues<dim> state_values(
        state_fe_,
        quadrature,
        dealii::update_values | dealii::update_gradients |
          dealii::update_quadrature_points | dealii::update_JxW_values);
      dealii::FEValues<dim> control_values(control_fe_, quadrature,
                                            dealii::update_values);
      const dealii::QGauss<dim - 1> face_quadrature(quadrature_order);
      dealii::FEFaceValues<dim> state_face_values(
        state_fe_,
        face_quadrature,
        dealii::update_values | dealii::update_gradients |
          dealii::update_normal_vectors | dealii::update_quadrature_points |
          dealii::update_JxW_values);

      dealii::FullMatrix<double> local_system(state_fe_.dofs_per_cell,
                                              state_fe_.dofs_per_cell);
      dealii::FullMatrix<double> local_state_tracking(
        state_fe_.dofs_per_cell, state_fe_.dofs_per_cell);
      dealii::FullMatrix<double> local_control_coupling(
        state_fe_.dofs_per_cell, control_fe_.dofs_per_cell);
      dealii::FullMatrix<double> local_control_mass(control_fe_.dofs_per_cell,
                                                    control_fe_.dofs_per_cell);
      dealii::Vector<double> local_forcing(state_fe_.dofs_per_cell);
      dealii::Vector<double> local_desired_state(state_fe_.dofs_per_cell);
      std::vector<dealii::types::global_dof_index> state_indices(
        state_fe_.dofs_per_cell);
      std::vector<dealii::types::global_dof_index> control_indices(
        control_fe_.dofs_per_cell);

      auto state_cell = state_dof_handler_.begin_active();
      auto control_cell = control_dof_handler_.begin_active();
      for (; state_cell != state_dof_handler_.end(); ++state_cell, ++control_cell)
        {
          contract::require(control_cell != control_dof_handler_.end(),
                            "State and control DoF handlers do not share cells");
          state_values.reinit(state_cell);
          control_values.reinit(control_cell);
          local_system = 0.0;
          local_state_tracking = 0.0;
          local_control_coupling = 0.0;
          local_control_mass = 0.0;
          local_forcing = 0.0;
          local_desired_state = 0.0;
          const bool observe_cell =
            !uses_point_sensor_ && !uses_normal_flux_ &&
            (observation_material_ids_.empty() ||
             observation_material_ids_.find(state_cell->material_id()) !=
               observation_material_ids_.end());

          for (unsigned int q = 0; q < quadrature.size(); ++q)
            {
              const double weight = state_values.JxW(q);
              const double forcing_value = forcing.value(state_values.quadrature_point(q));
              const double desired_value = observe_cell
                                             ? desired_state.value(
                                                 state_values.quadrature_point(q))
                                             : 0.0;
              const dealii::Tensor<1, dim> desired_gradient =
                observe_cell && uses_h1_state_observation_
                  ? desired_state.gradient(state_values.quadrature_point(q))
                  : dealii::Tensor<1, dim>();
              if (observe_cell)
                desired_state_norm_ +=
                  (desired_value * desired_value +
                   desired_gradient * desired_gradient) *
                  weight;
              dealii::Tensor<2, dim> diffusion_tensor;
              dealii::Tensor<1, dim> conservative_transport;
              dealii::Tensor<1, dim> advective_transport;
              double                reaction_value = reaction_;
              if (general_scalar_data != nullptr)
                {
                  const auto &point = state_values.quadrature_point(q);
                  diffusion_tensor =
                    general_scalar_data->diffusion_tensor.value(point);
                  conservative_transport =
                    general_scalar_data->conservative_transport.value(point);
                  advective_transport =
                    general_scalar_data->advective_transport.value(point);
                  reaction_value = general_scalar_data->reaction.value(point);
                }

              for (unsigned int i = 0; i < state_fe_.dofs_per_cell; ++i)
                {
                  const double phi_i = state_values.shape_value(i, q);
                  const auto   grad_i = state_values.shape_grad(i, q);
                  local_forcing(i) += forcing_value * phi_i * weight;
                  if (observe_cell)
                    local_desired_state(i) +=
                      (desired_value * phi_i +
                       (uses_h1_state_observation_ ? desired_gradient * grad_i :
                                                    0.0)) *
                      weight;
                  for (unsigned int j = 0; j < state_fe_.dofs_per_cell; ++j)
                    {
                      const double phi_j = state_values.shape_value(j, q);
                      const auto   grad_j = state_values.shape_grad(j, q);
                      if (general_scalar_data != nullptr)
                        {
                          local_system(i, j) +=
                            ((diffusion_tensor * grad_j) * grad_i -
                             phi_j * (conservative_transport * grad_i) +
                             (advective_transport * grad_j) * phi_i +
                             reaction_value * phi_i * phi_j) *
                            weight;
                        }
                      else
                        local_system(i, j) +=
                          (diffusion_ * (state_values.shape_grad(i, q) *
                                         state_values.shape_grad(j, q)) +
                           reaction_ * phi_i * phi_j) *
                          weight;
                      if (observe_cell)
                        local_state_tracking(i, j) +=
                          (phi_i * phi_j +
                           (uses_h1_state_observation_ ? grad_i * grad_j : 0.0)) *
                          weight;
                    }
                  for (unsigned int j = 0; j < control_fe_.dofs_per_cell; ++j)
                    local_control_coupling(i, j) +=
                      phi_i * control_values.shape_value(j, q) * weight;
                }
              for (unsigned int i = 0; i < control_fe_.dofs_per_cell; ++i)
                for (unsigned int j = 0; j < control_fe_.dofs_per_cell; ++j)
                  local_control_mass(i, j) +=
                    control_values.shape_value(i, q) *
                    control_values.shape_value(j, q) * weight;
            }

          if (general_scalar_data != nullptr)
            for (unsigned int face = 0;
                 face < dealii::GeometryInfo<dim>::faces_per_cell;
                 ++face)
              if (state_cell->face(face)->at_boundary() &&
                  robin_boundary_ids_.count(
                    state_cell->face(face)->boundary_id()) != 0)
                {
                  state_face_values.reinit(state_cell, face);
                  for (unsigned int q = 0; q < face_quadrature.size(); ++q)
                    {
                      const auto point = state_face_values.quadrature_point(q);
                      const double robin_coefficient =
                        general_scalar_data->robin_coefficient.value(point);
                      const double robin_source =
                        general_scalar_data->robin_source.value(point);
                      const double weight = state_face_values.JxW(q);
                      for (unsigned int i = 0; i < state_fe_.dofs_per_cell; ++i)
                        {
                          const double phi_i =
                            state_face_values.shape_value(i, q);
                          local_forcing(i) += robin_source * phi_i * weight;
                          for (unsigned int j = 0;
                               j < state_fe_.dofs_per_cell;
                               ++j)
                            local_system(i, j) +=
                              robin_coefficient * phi_i *
                              state_face_values.shape_value(j, q) * weight;
                        }
                    }
                }

          state_cell->get_dof_indices(state_indices);
          control_cell->get_dof_indices(control_indices);
          for (unsigned int i = 0; i < state_fe_.dofs_per_cell; ++i)
            {
              forcing_load_[state_indices[i]] += local_forcing(i);
              desired_state_load_[state_indices[i]] += local_desired_state(i);
              for (unsigned int j = 0; j < state_fe_.dofs_per_cell; ++j)
                {
                  physical_system_matrix_.add(state_indices[i],
                                              state_indices[j],
                                              local_system(i, j));
                  physical_state_tracking_operator_.add(
                    state_indices[i],
                    state_indices[j],
                    local_state_tracking(i, j));
                }
              for (unsigned int j = 0; j < control_fe_.dofs_per_cell; ++j)
                physical_control_coupling_.add(state_indices[i],
                                               control_indices[j],
                                               local_control_coupling(i, j));
            }
          for (unsigned int i = 0; i < control_fe_.dofs_per_cell; ++i)
            for (unsigned int j = 0; j < control_fe_.dofs_per_cell; ++j)
              control_mass_->add(control_indices[i],
                                 control_indices[j],
                                 local_control_mass(i, j));
        }
      contract::require(control_cell == control_dof_handler_.end(),
                        "State and control DoF handlers have different cells");
      if (uses_point_sensor_)
        assemble_point_sensor_operator(desired_state);
      if (uses_normal_flux_)
        assemble_normal_flux_operator(desired_state);
    }

    void
    assemble_point_sensor_operator(const dealii::Function<dim> &desired_state)
    {
      dealii::MappingQ1<dim> mapping;
      std::vector<dealii::types::global_dof_index> state_indices(
        state_fe_.dofs_per_cell);
      for (const auto &coordinate : point_sensor_coordinates_)
        {
          contract::require(coordinate.size() == dim,
                            "A point sensor coordinate has the wrong mesh dimension");
          dealii::Point<dim> point;
          for (unsigned int component = 0; component < dim; ++component)
            point[component] = coordinate[component];
          const auto cell_and_reference_point =
            dealii::GridTools::find_active_cell_around_point(
              mapping, state_dof_handler_, point);
          contract::require(cell_and_reference_point.first !=
                              state_dof_handler_.end(),
                            "Every point sensor must lie in the compiled mesh");
          const auto &cell = cell_and_reference_point.first;
          const auto &reference_point = cell_and_reference_point.second;
          cell->get_dof_indices(state_indices);
          const double target = desired_state.value(point);
          desired_state_norm_ += target * target;
          Vector point_evaluation(state_dof_handler_.n_dofs());
          for (unsigned int i = 0; i < state_fe_.dofs_per_cell; ++i)
            {
              const double phi_i =
                state_fe_.shape_value(i, reference_point);
              point_evaluation[state_indices[i]] = phi_i;
              desired_state_load_[state_indices[i]] += target * phi_i;
              for (unsigned int j = 0; j < state_fe_.dofs_per_cell; ++j)
                physical_state_tracking_operator_.add(
                  state_indices[i],
                  state_indices[j],
                  phi_i * state_fe_.shape_value(j, reference_point));
            }
          point_sensor_evaluations_.push_back(std::move(point_evaluation));
        }
    }

    void
    assemble_normal_flux_operator(const dealii::Function<dim> &desired_state)
    {
      const dealii::QGauss<dim - 1> face_quadrature(
        std::max(state_fe_.degree, control_fe_.degree) + 2);
      dealii::FEFaceValues<dim> face_values(
        state_fe_,
        face_quadrature,
        dealii::update_gradients | dealii::update_normal_vectors |
          dealii::update_quadrature_points | dealii::update_JxW_values);
      std::vector<dealii::types::global_dof_index> state_indices(
        state_fe_.dofs_per_cell);
      for (auto cell = state_dof_handler_.begin_active();
           cell != state_dof_handler_.end();
           ++cell)
        for (unsigned int face = 0;
             face < dealii::GeometryInfo<dim>::faces_per_cell;
             ++face)
          if (cell->face(face)->at_boundary() &&
              normal_flux_boundary_ids_.count(
                cell->face(face)->boundary_id()) != 0)
            {
              face_values.reinit(cell, face);
              cell->get_dof_indices(state_indices);
              for (unsigned int q = 0; q < face_quadrature.size(); ++q)
                {
                  const double weight = face_values.JxW(q);
                  const double target =
                    desired_state.value(face_values.quadrature_point(q));
                  desired_state_norm_ += target * target * weight;
                  Vector normal_flux_evaluation(state_dof_handler_.n_dofs());
                  for (unsigned int i = 0; i < state_fe_.dofs_per_cell; ++i)
                    {
                      const double flux_i =
                        face_values.shape_grad(i, q) *
                        face_values.normal_vector(q);
                      const auto global_i = state_indices[i];
                      if (homogeneous_constraints_.is_constrained(global_i))
                        continue;
                      normal_flux_evaluation[global_i] = flux_i;
                      desired_state_load_[global_i] += target * flux_i * weight;
                      for (unsigned int j = 0;
                           j < state_fe_.dofs_per_cell;
                           ++j)
                        {
                          const auto global_j = state_indices[j];
                          if (!homogeneous_constraints_.is_constrained(global_j))
                            physical_state_tracking_operator_.add(
                              global_i,
                              global_j,
                              flux_i *
                                (face_values.shape_grad(j, q) *
                                 face_values.normal_vector(q)) *
                                 weight);
                        }
                    }
                  normal_flux_evaluations_.push_back(
                    std::move(normal_flux_evaluation));
                  normal_flux_quadrature_weights_.push_back(weight);
                }
            }
    }

    void
    assemble_reduced_solve_operators()
    {
      build_reduced_state_matrix();
      build_reduced_control_coupling();
      Vector lifting_contribution(state_dof_handler_.n_dofs());
      physical_system_matrix_.vmult(lifting_contribution, lifting_);
      Vector right_hand_side = forcing_load_;
      right_hand_side.add(-1.0, lifting_contribution);
      reduced_forcing_load_ = pullback(right_hand_side);
    }

    void
    build_reduced_state_matrix()
    {
      const auto size = independent_state_dofs_.size();
      dealii::DynamicSparsityPattern dsp(size, size);
      for (std::size_t column = 0; column < size; ++column)
        {
          Vector basis(size);
          basis[column] = 1.0;
          Vector physical_column(state_dof_handler_.n_dofs());
          reconstruction_.vmult(physical_column, basis);
          Vector physical_result(state_dof_handler_.n_dofs());
          physical_system_matrix_.vmult(physical_result, physical_column);
          const Vector reduced_result = pullback(physical_result);
          for (std::size_t row = 0; row < size; ++row)
            if (reduced_result[row] != 0.0)
              dsp.add(row, column);
        }
      reduced_state_sparsity_.copy_from(dsp);
      reduced_system_matrix_.reinit(reduced_state_sparsity_);
      for (std::size_t column = 0; column < size; ++column)
        {
          Vector basis(size);
          basis[column] = 1.0;
          Vector physical_column(state_dof_handler_.n_dofs());
          reconstruction_.vmult(physical_column, basis);
          Vector physical_result(state_dof_handler_.n_dofs());
          physical_system_matrix_.vmult(physical_result, physical_column);
          const Vector reduced_result = pullback(physical_result);
          for (std::size_t row = 0; row < size; ++row)
            if (reduced_result[row] != 0.0)
              reduced_system_matrix_.set(row, column, reduced_result[row]);
        }
    }

    void
    build_reduced_control_coupling()
    {
      const auto state_size = independent_state_dofs_.size();
      const auto control_size = control_dof_handler_.n_dofs();
      dealii::DynamicSparsityPattern dsp(state_size, control_size);
      for (std::size_t column = 0; column < control_size; ++column)
        {
          Vector basis(control_size);
          basis[column] = 1.0;
          Vector physical_result(state_dof_handler_.n_dofs());
          physical_control_coupling_.vmult(physical_result, basis);
          const Vector reduced_result = pullback(physical_result);
          for (std::size_t row = 0; row < state_size; ++row)
            if (reduced_result[row] != 0.0)
              dsp.add(row, column);
        }
      reduced_control_sparsity_.copy_from(dsp);
      reduced_control_coupling_.reinit(reduced_control_sparsity_);
      for (std::size_t column = 0; column < control_size; ++column)
        {
          Vector basis(control_size);
          basis[column] = 1.0;
          Vector physical_result(state_dof_handler_.n_dofs());
          physical_control_coupling_.vmult(physical_result, basis);
          const Vector reduced_result = pullback(physical_result);
          for (std::size_t row = 0; row < state_size; ++row)
            if (reduced_result[row] != 0.0)
              reduced_control_coupling_.set(row, column, reduced_result[row]);
        }
    }

    Vector
    reconstruct(const Vector &independent_state) const
    {
      contract::require(independent_state.size() == independent_state_dofs_.size(),
                        "State reconstruction received incompatible coordinates");
      Vector physical(state_dof_handler_.n_dofs());
      reconstruction_.vmult(physical, independent_state);
      physical.add(1.0, lifting_);
      return physical;
    }

    Vector
    embed_tangent(const Vector &independent_state) const
    {
      contract::require(independent_state.size() == independent_state_dofs_.size(),
                        "State tangent has incompatible independent coordinates");
      Vector physical(state_dof_handler_.n_dofs());
      reconstruction_.vmult(physical, independent_state);
      return physical;
    }

    Vector
    pullback(const Vector &physical_covector) const
    {
      contract::require(physical_covector.size() == state_dof_handler_.n_dofs(),
                        "State pullback received an incompatible physical covector");
      Vector independent(independent_state_dofs_.size());
      reconstruction_.Tvmult(independent, physical_covector);
      return independent;
    }

    static contract::LinearSolveReport
    solve_symmetric_system(const dealii::SparseMatrix<double> &matrix,
                           Vector &                              solution,
                           const Vector &                        right_hand_side,
                           const dealii_backend::SPDLinearSolvePolicy &policy)
    {
      return dealii_backend::solve_serial_spd(matrix,
                                              solution,
                                              right_hand_side,
                                              policy);
    }

    dealii::FE_Q<dim> state_fe_;
    dealii::FE_DGQ<dim> control_fe_;
    dealii::DoFHandler<dim> state_dof_handler_;
    dealii::DoFHandler<dim> control_dof_handler_;
    dealii::AffineConstraints<double> homogeneous_constraints_;
    dealii::AffineConstraints<double> physical_constraints_;
    std::vector<dealii::types::global_dof_index> independent_state_dofs_;

    const double diffusion_;
    const double reaction_;
    const double regularisation_weight_;
    const bool uses_general_scalar_;
    const std::set<dealii::types::boundary_id> dirichlet_boundary_ids_;
    const std::set<dealii::types::material_id> observation_material_ids_;
    const bool uses_h1_state_observation_;
    const std::set<dealii::types::boundary_id> normal_flux_boundary_ids_;
    const bool uses_normal_flux_;
    const std::vector<std::vector<double>> point_sensor_coordinates_;
    const bool uses_point_sensor_;
    const std::set<dealii::types::boundary_id> robin_boundary_ids_;

    dealii::SparsityPattern reconstruction_sparsity_;
    dealii::SparseMatrix<double> reconstruction_;
    Vector lifting_;

    dealii::SparsityPattern physical_state_sparsity_;
    dealii::SparsityPattern physical_control_sparsity_;
    dealii::SparsityPattern control_mass_sparsity_;
    dealii::SparseMatrix<double> physical_system_matrix_;
    dealii::SparseMatrix<double> physical_state_tracking_operator_;
    dealii::SparseMatrix<double> physical_control_coupling_;
    std::shared_ptr<dealii::SparseMatrix<double>> control_mass_;
    Vector forcing_load_;
    Vector desired_state_load_;
    std::vector<Vector> point_sensor_evaluations_;
    std::vector<Vector> normal_flux_evaluations_;
    std::vector<double> normal_flux_quadrature_weights_;
    double desired_state_norm_ = 0.0;

    dealii::SparsityPattern reduced_state_sparsity_;
    dealii::SparsityPattern reduced_control_sparsity_;
    dealii::SparseMatrix<double> reduced_system_matrix_;
    dealii::SparseMatrix<double> reduced_control_coupling_;
    Vector reduced_forcing_load_;
    std::unique_ptr<dealii::SparseDirectUMFPACK> nonsymmetric_solver_;

    contract::LayoutPtr variable_layout_;
    contract::LayoutPtr test_layout_;
    contract::LayoutPtr state_layout_;
    contract::LayoutPtr control_layout_;
  };
} // namespace nmopt::compiler::v1::detail
