#pragma once

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
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/numerics/vector_tools.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace nmopt::compiler::v1::detail
{
  // P3.1's v1-only coefficient-identification target. The physical parameter
  // m is cellwise constant and positive through the separately compiled box:
  //
  //   r(y,m) = A(m) y - f_h,
  //   <A(m)y,v> = (m grad y, grad v) + (c y,v),
  //   J(y,m) = J_tracking(y) + alpha/2 m^T M_m m.
  //
  // The state matrix is deliberately reassembled for every parameter point.
  // This realizes the nonlinear first-order actions without changing the v0
  // constant-diffusion reference model or the generic reduced DTO contract.
  template <int dim>
  class CoefficientIdentificationModel final
    : public contract::ExecutableModelT<dealii_backend::SerialBackend>
  {
  public:
    using Backend = dealii_backend::SerialBackend;
    using Matrix = dealii::SparseMatrix<double>;
    using Vector = dealii::Vector<double>;
    using Primal = contract::PrimalBlockT<Backend>;
    using Covector = contract::CovectorBlockT<Backend>;
    using SolveResult = contract::FormulationSolveResultT<Backend>;

    CoefficientIdentificationModel(
      dealii::Triangulation<dim> &                triangulation,
      const dealii::Function<dim> &               forcing,
      const dealii::Function<dim> &               desired_state,
      const double                                  reaction,
      const double                                  regularisation_weight,
      const unsigned int                            state_degree,
      std::set<dealii::types::boundary_id>         dirichlet_boundary_ids)
      : state_fe_(state_degree)
      , parameter_fe_(0)
      , state_dof_handler_(triangulation)
      , parameter_dof_handler_(triangulation)
      , reaction_(reaction)
      , regularisation_weight_(regularisation_weight)
      , dirichlet_boundary_ids_(std::move(dirichlet_boundary_ids))
    {
      contract::require(reaction_ >= 0.0,
                        "Reaction coefficient must be non-negative");
      contract::require(regularisation_weight_ > 0.0,
                        "Parameter regularisation weight must be strictly positive");
      contract::require(state_degree > 0,
                        "State FE degree must be at least one");
      contract::require(!dirichlet_boundary_ids_.empty(),
                        "Coefficient-identification target needs a fixed Dirichlet boundary");

      state_dof_handler_.distribute_dofs(state_fe_);
      parameter_dof_handler_.distribute_dofs(parameter_fe_);
      build_state_constraints();
      initialise_storage();
      assemble_static(forcing, desired_state);
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

    dealii_backend::MassMetric
    parameter_l2_metric(
      dealii_backend::MassMetricSolveParameters solve_parameters = {}) const
    {
      return dealii_backend::MassMetric("l2_cellwise_parameter",
                                        parameter_layout_,
                                        parameter_mass_,
                                        solve_parameters);
    }

    dealii_backend::CellwiseBoxConstraint
    parameter_l2_box_constraint(
      Vector                              lower,
      Vector                              upper,
      const dealii_backend::MassMetric & projection_metric) const
    {
      return dealii_backend::CellwiseBoxConstraint(parameter_layout_,
                                                    std::move(lower),
                                                    std::move(upper),
                                                    projection_metric);
    }

    dealii_backend::CellwiseBoxConstraint
    parameter_l2_box_constraint(
      const double                        lower,
      const double                        upper,
      const dealii_backend::MassMetric & projection_metric) const
    {
      return dealii_backend::CellwiseBoxConstraint(parameter_layout_,
                                                    lower,
                                                    upper,
                                                    projection_metric);
    }

    Covector
    residual(const Primal &variables) const override
    {
      require_variables(variables, "Residual");
      Matrix state_matrix;
      assemble_state_matrix(variables.block(1), state_matrix);
      Vector value(state_dof_handler_.n_dofs());
      state_matrix.vmult(value, variables.block(0));
      value.add(-1.0, forcing_load_);
      return Covector(test_layout_, {std::move(value)});
    }

    Covector
    residual_jvp(const Primal &variables,
                 const Primal &variable_tangent) const override
    {
      require_variables(variables, "Residual JVP");
      require_variable_tangent(variable_tangent, "Residual JVP tangent");
      Matrix state_matrix;
      assemble_state_matrix(variables.block(1), state_matrix);
      Vector value(state_dof_handler_.n_dofs());
      state_matrix.vmult(value, variable_tangent.block(0));
      value.add(1.0,
                parameter_derivative_apply(variable_tangent.block(1),
                                           variables.block(0)));
      return Covector(test_layout_, {std::move(value)});
    }

    Covector
    residual_vjp(const Primal &variables, const Primal &test_seed) const override
    {
      require_variables(variables, "Residual VJP");
      contract::require(test_seed.layout()->compatible_with(*test_layout_),
                        "Residual VJP seed has an incompatible test layout");
      Matrix state_matrix;
      assemble_state_matrix(variables.block(1), state_matrix);
      Vector state(state_dof_handler_.n_dofs());
      state_matrix.Tvmult(state, test_seed.block(0));
      Vector parameter = parameter_derivative_transpose(variables.block(0),
                                                         test_seed.block(0));
      return Covector(variable_layout_, {std::move(state), std::move(parameter)});
    }

    double
    objective(const Primal &variables) const override
    {
      require_variables(variables, "Objective");
      Vector state_mass_times_state(state_dof_handler_.n_dofs());
      state_mass_.vmult(state_mass_times_state, variables.block(0));
      const double state_value =
        0.5 * (variables.block(0) * state_mass_times_state) -
        (desired_state_load_ * variables.block(0)) + 0.5 * desired_state_norm_;

      Vector parameter_mass_times_parameter(parameter_dof_handler_.n_dofs());
      parameter_mass_->vmult(parameter_mass_times_parameter, variables.block(1));
      const double parameter_value = 0.5 * regularisation_weight_ *
                                     (variables.block(1) *
                                      parameter_mass_times_parameter);
      return state_value + parameter_value;
    }

    Covector
    objective_derivative(const Primal &variables) const override
    {
      require_variables(variables, "Objective derivative");
      Vector state(state_dof_handler_.n_dofs());
      state_mass_.vmult(state, variables.block(0));
      state.add(-1.0, desired_state_load_);

      Vector parameter(parameter_dof_handler_.n_dofs());
      parameter_mass_->vmult(parameter, variables.block(1));
      parameter *= regularisation_weight_;
      return Covector(variable_layout_, {std::move(state), std::move(parameter)});
    }

    Primal
    solve_state(const Primal &parameter) const
    {
      auto result = solve_state_with_report(parameter, {});
      contract::require(result.report.converged(),
                        "State solve did not converge under its declared policy");
      return std::move(result.solution);
    }

    SolveResult
    solve_state_with_report(
      const Primal &                              parameter,
      const dealii_backend::SPDLinearSolvePolicy &policy) const
    {
      contract::require(parameter.layout()->compatible_with(*parameter_layout_),
                        "State solve parameter has an incompatible layout");
      require_positive_parameter(parameter.block(0));
      Matrix state_matrix;
      assemble_state_matrix(parameter.block(0), state_matrix);
      Vector state(state_dof_handler_.n_dofs());
      auto report = solve_symmetric_system(state_matrix,
                                           state,
                                           forcing_load_,
                                           policy);
      state_constraints_.distribute(state);
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
      Matrix state_matrix;
      assemble_state_matrix(full_point.block(1), state_matrix);
      Vector adjoint(test_layout_->dimension(0));
      auto report = solve_symmetric_system(state_matrix,
                                           adjoint,
                                           state_objective_derivative.block(0),
                                           policy);
      state_constraints_.distribute(adjoint);
      return {Primal(test_layout_, {std::move(adjoint)}), std::move(report)};
    }

  private:
    void
    require_variables(const Primal &variables, const char *operation) const
    {
      contract::require(
        variables.layout()->compatible_with(*variable_layout_),
        std::string(operation) + " received an incompatible variable layout");
      require_positive_parameter(variables.block(1));
    }

    void
    require_variable_tangent(const Primal &tangent, const char *operation) const
    {
      contract::require(
        tangent.layout()->compatible_with(*variable_layout_),
        std::string(operation) + " received an incompatible variable layout");
    }

    void
    require_positive_parameter(const Vector &parameter) const
    {
      for (dealii::types::global_dof_index index = 0; index < parameter.size();
           ++index)
        contract::require(std::isfinite(parameter[index]) && parameter[index] > 0.0,
                          "Diffusion parameter coefficients must be strictly positive");
    }

    void
    build_state_constraints()
    {
      state_constraints_.clear();
      dealii::DoFTools::make_hanging_node_constraints(state_dof_handler_,
                                                       state_constraints_);
      dealii::Functions::ZeroFunction<dim> zero;
      for (const auto boundary_id : dirichlet_boundary_ids_)
        dealii::VectorTools::interpolate_boundary_values(state_dof_handler_,
                                                          boundary_id,
                                                          zero,
                                                          state_constraints_);
      state_constraints_.close();

      std::map<dealii::types::global_dof_index, double> dirichlet_values;
      for (const auto boundary_id : dirichlet_boundary_ids_)
        dealii::VectorTools::interpolate_boundary_values(state_dof_handler_,
                                                          boundary_id,
                                                          zero,
                                                          dirichlet_values);
      constrained_state_dofs_.assign(state_dof_handler_.n_dofs(), false);
      for (const auto &entry : dirichlet_values)
        constrained_state_dofs_.at(entry.first) = true;
      for (dealii::types::global_dof_index index = 0;
           index < state_dof_handler_.n_dofs();
           ++index)
        if (state_constraints_.is_constrained(index) &&
            !constrained_state_dofs_.at(index))
          throw contract::ContractError(
            "The coefficient-identification target does not support hanging, periodic, "
            "or other non-Dirichlet affine state constraints");
    }

    void
    initialise_storage()
    {
      const auto state_size = state_dof_handler_.n_dofs();
      const auto parameter_size = parameter_dof_handler_.n_dofs();
      variable_layout_ = std::make_shared<const contract::BlockLayout>(
        "coefficient_identification_variables",
        std::vector<contract::SpaceId>{{"state"}, {"parameter"}},
        std::vector<std::size_t>{state_size, parameter_size});
      test_layout_ = std::make_shared<const contract::BlockLayout>(
        "coefficient_identification_state_test",
        std::vector<contract::SpaceId>{{"state_test"}},
        std::vector<std::size_t>{state_size});
      state_layout_ = variable_layout_->single_block(0, "state");
      parameter_layout_ = variable_layout_->single_block(1, "parameter");

      dealii::DynamicSparsityPattern state_dsp(state_size, state_size);
      dealii::DoFTools::make_sparsity_pattern(state_dof_handler_, state_dsp);
      state_sparsity_.copy_from(state_dsp);
      state_mass_.reinit(state_sparsity_);

      dealii::DynamicSparsityPattern parameter_dsp(parameter_size,
                                                   parameter_size);
      dealii::DoFTools::make_sparsity_pattern(parameter_dof_handler_,
                                               parameter_dsp);
      parameter_sparsity_.copy_from(parameter_dsp);
      parameter_mass_ = std::make_shared<Matrix>();
      parameter_mass_->reinit(parameter_sparsity_);

      forcing_load_.reinit(state_size);
      desired_state_load_.reinit(state_size);
    }

    void
    assemble_static(const dealii::Function<dim> &forcing,
                    const dealii::Function<dim> &desired_state)
    {
      const dealii::QGauss<dim> quadrature(state_fe_.degree + 2);
      dealii::FEValues<dim> state_values(
        state_fe_,
        quadrature,
        dealii::update_values | dealii::update_quadrature_points |
          dealii::update_JxW_values);
      dealii::FEValues<dim> parameter_values(parameter_fe_,
                                              quadrature,
                                              dealii::update_values);
      dealii::FullMatrix<double> local_state_mass(state_fe_.dofs_per_cell,
                                                  state_fe_.dofs_per_cell);
      dealii::FullMatrix<double> local_parameter_mass(
        parameter_fe_.dofs_per_cell, parameter_fe_.dofs_per_cell);
      dealii::Vector<double> local_forcing(state_fe_.dofs_per_cell);
      dealii::Vector<double> local_desired_state(state_fe_.dofs_per_cell);
      std::vector<dealii::types::global_dof_index> state_indices(
        state_fe_.dofs_per_cell);
      std::vector<dealii::types::global_dof_index> parameter_indices(
        parameter_fe_.dofs_per_cell);

      auto state_cell = state_dof_handler_.begin_active();
      auto parameter_cell = parameter_dof_handler_.begin_active();
      for (; state_cell != state_dof_handler_.end();
           ++state_cell, ++parameter_cell)
        {
          contract::require(parameter_cell != parameter_dof_handler_.end(),
                            "State and parameter DoF handlers do not share cells");
          state_values.reinit(state_cell);
          parameter_values.reinit(parameter_cell);
          local_state_mass = 0.0;
          local_parameter_mass = 0.0;
          local_forcing = 0.0;
          local_desired_state = 0.0;
          for (unsigned int q = 0; q < quadrature.size(); ++q)
            {
              const double weight = state_values.JxW(q);
              const double forcing_value =
                forcing.value(state_values.quadrature_point(q));
              const double desired_value =
                desired_state.value(state_values.quadrature_point(q));
              desired_state_norm_ += desired_value * desired_value * weight;
              for (unsigned int i = 0; i < state_fe_.dofs_per_cell; ++i)
                {
                  const double phi_i = state_values.shape_value(i, q);
                  local_forcing(i) += forcing_value * phi_i * weight;
                  local_desired_state(i) += desired_value * phi_i * weight;
                  for (unsigned int j = 0; j < state_fe_.dofs_per_cell; ++j)
                    local_state_mass(i, j) +=
                      phi_i * state_values.shape_value(j, q) * weight;
                }
              for (unsigned int i = 0; i < parameter_fe_.dofs_per_cell; ++i)
                for (unsigned int j = 0; j < parameter_fe_.dofs_per_cell; ++j)
                  local_parameter_mass(i, j) +=
                    parameter_values.shape_value(i, q) *
                    parameter_values.shape_value(j, q) * weight;
            }

          state_cell->get_dof_indices(state_indices);
          parameter_cell->get_dof_indices(parameter_indices);
          for (unsigned int i = 0; i < state_fe_.dofs_per_cell; ++i)
            {
              const auto global_i = state_indices[i];
              if (constrained_state_dofs_.at(global_i))
                continue;
              forcing_load_[global_i] += local_forcing(i);
              desired_state_load_[global_i] += local_desired_state(i);
              for (unsigned int j = 0; j < state_fe_.dofs_per_cell; ++j)
                {
                  const auto global_j = state_indices[j];
                  if (!constrained_state_dofs_.at(global_j))
                    state_mass_.add(global_i, global_j, local_state_mass(i, j));
                }
            }
          for (unsigned int i = 0; i < parameter_fe_.dofs_per_cell; ++i)
            for (unsigned int j = 0; j < parameter_fe_.dofs_per_cell; ++j)
              parameter_mass_->add(parameter_indices[i],
                                   parameter_indices[j],
                                   local_parameter_mass(i, j));
        }
      contract::require(parameter_cell == parameter_dof_handler_.end(),
                        "State and parameter DoF handlers have different cells");
    }

    void
    assemble_state_matrix(const Vector &parameter, Matrix &matrix) const
    {
      require_positive_parameter(parameter);
      matrix.reinit(state_sparsity_);
      const dealii::QGauss<dim> quadrature(state_fe_.degree + 2);
      dealii::FEValues<dim> state_values(state_fe_,
                                          quadrature,
                                          dealii::update_values |
                                            dealii::update_gradients |
                                            dealii::update_JxW_values);
      dealii::FullMatrix<double> local_matrix(state_fe_.dofs_per_cell,
                                              state_fe_.dofs_per_cell);
      std::vector<dealii::types::global_dof_index> state_indices(
        state_fe_.dofs_per_cell);
      std::vector<dealii::types::global_dof_index> parameter_indices(
        parameter_fe_.dofs_per_cell);
      auto state_cell = state_dof_handler_.begin_active();
      auto parameter_cell = parameter_dof_handler_.begin_active();
      for (; state_cell != state_dof_handler_.end();
           ++state_cell, ++parameter_cell)
        {
          contract::require(parameter_cell != parameter_dof_handler_.end(),
                            "State and parameter DoF handlers do not share cells");
          state_values.reinit(state_cell);
          state_cell->get_dof_indices(state_indices);
          parameter_cell->get_dof_indices(parameter_indices);
          const double coefficient = parameter[parameter_indices.at(0)];
          local_matrix = 0.0;
          for (unsigned int q = 0; q < quadrature.size(); ++q)
            for (unsigned int i = 0; i < state_fe_.dofs_per_cell; ++i)
              for (unsigned int j = 0; j < state_fe_.dofs_per_cell; ++j)
                local_matrix(i, j) +=
                  (coefficient *
                     (state_values.shape_grad(i, q) *
                      state_values.shape_grad(j, q)) +
                   reaction_ * state_values.shape_value(i, q) *
                     state_values.shape_value(j, q)) *
                  state_values.JxW(q);
          for (unsigned int i = 0; i < state_fe_.dofs_per_cell; ++i)
            {
              const auto global_i = state_indices[i];
              if (constrained_state_dofs_.at(global_i))
                continue;
              for (unsigned int j = 0; j < state_fe_.dofs_per_cell; ++j)
                {
                  const auto global_j = state_indices[j];
                  if (!constrained_state_dofs_.at(global_j))
                    matrix.add(global_i, global_j, local_matrix(i, j));
                }
            }
        }
      contract::require(parameter_cell == parameter_dof_handler_.end(),
                        "State and parameter DoF handlers have different cells");
      for (dealii::types::global_dof_index index = 0;
           index < state_dof_handler_.n_dofs();
           ++index)
        if (constrained_state_dofs_.at(index))
          matrix.set(index, index, 1.0);
    }

    Vector
    parameter_derivative_apply(const Vector &parameter_tangent,
                               const Vector &state) const
    {
      Vector value(state_dof_handler_.n_dofs());
      const dealii::QGauss<dim> quadrature(state_fe_.degree + 2);
      dealii::FEValues<dim> state_values(state_fe_,
                                          quadrature,
                                          dealii::update_gradients |
                                            dealii::update_JxW_values);
      std::vector<dealii::Tensor<1, dim>> state_gradients(quadrature.size());
      std::vector<dealii::types::global_dof_index> state_indices(
        state_fe_.dofs_per_cell);
      std::vector<dealii::types::global_dof_index> parameter_indices(
        parameter_fe_.dofs_per_cell);
      auto state_cell = state_dof_handler_.begin_active();
      auto parameter_cell = parameter_dof_handler_.begin_active();
      for (; state_cell != state_dof_handler_.end();
           ++state_cell, ++parameter_cell)
        {
          contract::require(parameter_cell != parameter_dof_handler_.end(),
                            "State and parameter DoF handlers do not share cells");
          state_values.reinit(state_cell);
          state_values.get_function_gradients(state, state_gradients);
          state_cell->get_dof_indices(state_indices);
          parameter_cell->get_dof_indices(parameter_indices);
          const double coefficient_tangent =
            parameter_tangent[parameter_indices.at(0)];
          for (unsigned int q = 0; q < quadrature.size(); ++q)
            for (unsigned int i = 0; i < state_fe_.dofs_per_cell; ++i)
              if (!constrained_state_dofs_.at(state_indices[i]))
                value[state_indices[i]] += coefficient_tangent *
                  (state_gradients[q] * state_values.shape_grad(i, q)) *
                  state_values.JxW(q);
        }
      contract::require(parameter_cell == parameter_dof_handler_.end(),
                        "State and parameter DoF handlers have different cells");
      return value;
    }

    Vector
    parameter_derivative_transpose(const Vector &state,
                                   const Vector &test_seed) const
    {
      Vector unconstrained_test_seed = test_seed;
      for (dealii::types::global_dof_index index = 0;
           index < unconstrained_test_seed.size();
           ++index)
        if (constrained_state_dofs_.at(index))
          unconstrained_test_seed[index] = 0.0;

      Vector value(parameter_dof_handler_.n_dofs());
      const dealii::QGauss<dim> quadrature(state_fe_.degree + 2);
      dealii::FEValues<dim> state_values(state_fe_,
                                          quadrature,
                                          dealii::update_gradients |
                                            dealii::update_JxW_values);
      std::vector<dealii::Tensor<1, dim>> state_gradients(quadrature.size());
      std::vector<dealii::Tensor<1, dim>> test_gradients(quadrature.size());
      std::vector<dealii::types::global_dof_index> parameter_indices(
        parameter_fe_.dofs_per_cell);
      auto state_cell = state_dof_handler_.begin_active();
      auto parameter_cell = parameter_dof_handler_.begin_active();
      for (; state_cell != state_dof_handler_.end();
           ++state_cell, ++parameter_cell)
        {
          contract::require(parameter_cell != parameter_dof_handler_.end(),
                            "State and parameter DoF handlers do not share cells");
          state_values.reinit(state_cell);
          state_values.get_function_gradients(state, state_gradients);
          state_values.get_function_gradients(unconstrained_test_seed,
                                              test_gradients);
          parameter_cell->get_dof_indices(parameter_indices);
          for (unsigned int q = 0; q < quadrature.size(); ++q)
            value[parameter_indices.at(0)] +=
              (state_gradients[q] * test_gradients[q]) * state_values.JxW(q);
        }
      contract::require(parameter_cell == parameter_dof_handler_.end(),
                        "State and parameter DoF handlers have different cells");
      return value;
    }

    static contract::LinearSolveReport
    solve_symmetric_system(const Matrix &matrix,
                           Vector &      solution,
                           const Vector &right_hand_side,
                           const dealii_backend::SPDLinearSolvePolicy &policy)
    {
      return dealii_backend::solve_serial_spd(matrix,
                                              solution,
                                              right_hand_side,
                                              policy);
    }

    dealii::FE_Q<dim> state_fe_;
    dealii::FE_DGQ<dim> parameter_fe_;
    dealii::DoFHandler<dim> state_dof_handler_;
    dealii::DoFHandler<dim> parameter_dof_handler_;
    dealii::AffineConstraints<double> state_constraints_;
    std::vector<bool> constrained_state_dofs_;

    const double reaction_;
    const double regularisation_weight_;
    const std::set<dealii::types::boundary_id> dirichlet_boundary_ids_;

    dealii::SparsityPattern state_sparsity_;
    dealii::SparsityPattern parameter_sparsity_;
    dealii::SparseMatrix<double> state_mass_;
    std::shared_ptr<Matrix> parameter_mass_;
    Vector forcing_load_;
    Vector desired_state_load_;
    double desired_state_norm_ = 0.0;

    contract::LayoutPtr variable_layout_;
    contract::LayoutPtr test_layout_;
    contract::LayoutPtr state_layout_;
    contract::LayoutPtr parameter_layout_;
  };
} // namespace nmopt::compiler::v1::detail
