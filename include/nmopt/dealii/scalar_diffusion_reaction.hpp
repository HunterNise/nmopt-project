#pragma once

#include "nmopt/contract/executable_model.hpp"
#include "nmopt/dealii/serial_backend.hpp"

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
#include <memory>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace nmopt::dealii_backend
{
  // V0 deal.II lowerer:
  //
  //   -div(k grad y) + c y = f + u,
  //
  // with homogeneous Dirichlet data on the selected boundary ids, an FE_Q
  // state, FE_DGQ(0) volume control, distributed state tracking, and
  // L2 control regularisation. It produces the DTO contract
  //
  //   r(y,u) = A y - f_h - B u,
  //   J(y,u) = 1/2 y^T M_y y - q_y^T y + const
  //            + alpha/2 u^T M_u u.
  //
  // The implementation deliberately rejects hanging-node/periodic
  // constraints in this first lowerer. Homogeneous essential conditions are
  // represented by the recorded AffineConstraints and identity residual rows,
  // so residual, state solve, adjoint solve, and pullback share one discrete
  // coordinate convention.
  template <int dim>
  class ScalarDiffusionReactionModel final
    : public contract::ExecutableModelT<SerialBackend>
  {
  public:
    using Vector = dealii::Vector<double>;
    using Primal = contract::PrimalBlockT<SerialBackend>;
    using Covector = contract::CovectorBlockT<SerialBackend>;

    ScalarDiffusionReactionModel(
      dealii::Triangulation<dim> &                triangulation,
      const dealii::Function<dim> &               forcing,
      const dealii::Function<dim> &               desired_state,
      const double                                  diffusion,
      const double                                  reaction,
      const double                                  regularisation_weight,
      const unsigned int                            state_degree,
      std::set<dealii::types::boundary_id> dirichlet_boundary_ids = {0})
      : state_fe_(state_degree)
      , control_fe_(0)
      , state_dof_handler_(triangulation)
      , control_dof_handler_(triangulation)
      , diffusion_(diffusion)
      , reaction_(reaction)
      , regularisation_weight_(regularisation_weight)
      , dirichlet_boundary_ids_(std::move(dirichlet_boundary_ids))
    {
      contract::require(diffusion_ > 0.0,
                        "Diffusion coefficient must be strictly positive");
      contract::require(reaction_ >= 0.0,
                        "Reaction coefficient must be non-negative");
      contract::require(regularisation_weight_ > 0.0,
                        "Control regularisation weight must be strictly positive");
      contract::require(state_degree > 0,
                        "State FE degree must be at least one");
      contract::require(!dirichlet_boundary_ids_.empty(),
                        "The v0 lowerer requires a non-empty Dirichlet boundary");

      state_dof_handler_.distribute_dofs(state_fe_);
      control_dof_handler_.distribute_dofs(control_fe_);
      build_constraints();
      initialise_storage();
      assemble(forcing, desired_state);
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

    const dealii::AffineConstraints<double> &
    state_constraints() const
    {
      return state_constraints_;
    }

    Covector
    residual(const Primal &variables) const override
    {
      require_variables(variables, "Residual");

      Vector value(state_dof_handler_.n_dofs());
      system_matrix_.vmult(value, variables.block(0));
      value.add(-1.0, forcing_load_);

      Vector control_contribution(state_dof_handler_.n_dofs());
      control_coupling_.vmult(control_contribution, variables.block(1));
      value.add(-1.0, control_contribution);
      return Covector(test_layout_, {std::move(value)});
    }

    Covector
    residual_jvp(const Primal &variables,
                 const Primal &variable_tangent) const override
    {
      require_variables(variables, "Residual JVP");
      require_variables(variable_tangent, "Residual JVP tangent");

      Vector value(state_dof_handler_.n_dofs());
      system_matrix_.vmult(value, variable_tangent.block(0));
      Vector control_contribution(state_dof_handler_.n_dofs());
      control_coupling_.vmult(control_contribution,
                               variable_tangent.block(1));
      value.add(-1.0, control_contribution);
      return Covector(test_layout_, {std::move(value)});
    }

    Covector
    residual_vjp(const Primal &variables,
                 const Primal &test_seed) const override
    {
      require_variables(variables, "Residual VJP");
      contract::require(test_seed.layout()->compatible_with(*test_layout_),
                        "Residual VJP seed has an incompatible test layout");

      Vector state(state_dof_handler_.n_dofs());
      system_matrix_.Tvmult(state, test_seed.block(0));
      Vector control(control_dof_handler_.n_dofs());
      control_coupling_.Tvmult(control, test_seed.block(0));
      control *= -1.0;
      return Covector(variable_layout_, {std::move(state), std::move(control)});
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

      Vector control_mass_times_control(control_dof_handler_.n_dofs());
      control_mass_.vmult(control_mass_times_control, variables.block(1));
      const double control_value = 0.5 * regularisation_weight_ *
                                   (variables.block(1) *
                                    control_mass_times_control);
      return state_value + control_value;
    }

    Covector
    objective_derivative(const Primal &variables) const override
    {
      require_variables(variables, "Objective derivative");

      Vector state(state_dof_handler_.n_dofs());
      state_mass_.vmult(state, variables.block(0));
      state.add(-1.0, desired_state_load_);

      Vector control(control_dof_handler_.n_dofs());
      control_mass_.vmult(control, variables.block(1));
      control *= regularisation_weight_;
      return Covector(variable_layout_, {std::move(state), std::move(control)});
    }

    Primal
    solve_state(const Primal &control) const
    {
      contract::require(control.layout()->compatible_with(*control_layout_),
                        "State solve control has an incompatible layout");

      Vector right_hand_side = forcing_load_;
      Vector control_contribution(state_dof_handler_.n_dofs());
      control_coupling_.vmult(control_contribution, control.block(0));
      right_hand_side.add(1.0, control_contribution);

      Vector state(state_dof_handler_.n_dofs());
      solve_symmetric_system(state, right_hand_side);
      state_constraints_.distribute(state);
      return Primal(state_layout_, {std::move(state)});
    }

    Primal
    solve_adjoint(const Primal &full_point,
                  const Covector &state_objective_derivative) const
    {
      require_variables(full_point, "Adjoint solve point");
      contract::require(
        state_objective_derivative.layout()->compatible_with(*state_layout_),
        "Adjoint solve right-hand side has an incompatible state layout");

      // The v0 diffusion-reaction operator is symmetric. The VJP nevertheless
      // uses Tvmult above, so a non-symmetric extension cannot silently reuse
      // this solve path.
      Vector adjoint(test_layout_->dimension(0));
      solve_symmetric_system(adjoint, state_objective_derivative.block(0));
      state_constraints_.distribute(adjoint);
      return Primal(test_layout_, {std::move(adjoint)});
    }

  private:
    void
    require_variables(const Primal &variables, const char *operation) const
    {
      contract::require(
        variables.layout()->compatible_with(*variable_layout_),
        std::string(operation) + " received an incompatible variable layout");
    }

    void
    build_constraints()
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
            "The v0 deal.II lowerer does not support hanging, periodic, or "
            "other non-Dirichlet affine constraints");
    }

    void
    initialise_storage()
    {
      const auto state_size = state_dof_handler_.n_dofs();
      const auto control_size = control_dof_handler_.n_dofs();

      variable_layout_ = std::make_shared<const contract::BlockLayout>(
        "variables",
        std::vector<contract::SpaceId>{{"state"}, {"control"}},
        std::vector<std::size_t>{state_size, control_size});
      test_layout_ = std::make_shared<const contract::BlockLayout>(
        "state_test",
        std::vector<contract::SpaceId>{{"state_test"}},
        std::vector<std::size_t>{state_size});
      state_layout_ = variable_layout_->single_block(0, "state");
      control_layout_ = variable_layout_->single_block(1, "control");

      dealii::DynamicSparsityPattern state_dsp(state_size, state_size);
      dealii::DoFTools::make_sparsity_pattern(state_dof_handler_, state_dsp);
      state_sparsity_.copy_from(state_dsp);
      system_matrix_.reinit(state_sparsity_);
      state_mass_.reinit(state_sparsity_);

      dealii::DynamicSparsityPattern control_dsp(state_size, control_size);
      std::vector<dealii::types::global_dof_index> state_indices(
        state_fe_.dofs_per_cell);
      std::vector<dealii::types::global_dof_index> control_indices(
        control_fe_.dofs_per_cell);
      auto state_cell = state_dof_handler_.begin_active();
      auto control_cell = control_dof_handler_.begin_active();
      for (; state_cell != state_dof_handler_.end();
           ++state_cell, ++control_cell)
        {
          contract::require(control_cell != control_dof_handler_.end(),
                            "State and control DoF handlers do not share cells");
          state_cell->get_dof_indices(state_indices);
          control_cell->get_dof_indices(control_indices);
          for (unsigned int i = 0; i < state_fe_.dofs_per_cell; ++i)
            if (!constrained_state_dofs_.at(state_indices[i]))
              for (unsigned int j = 0; j < control_fe_.dofs_per_cell; ++j)
                control_dsp.add(state_indices[i], control_indices[j]);
        }
      contract::require(control_cell == control_dof_handler_.end(),
                        "State and control DoF handlers have different cells");
      control_sparsity_.copy_from(control_dsp);
      control_coupling_.reinit(control_sparsity_);

      dealii::DynamicSparsityPattern control_mass_dsp(control_size, control_size);
      dealii::DoFTools::make_sparsity_pattern(control_dof_handler_,
                                               control_mass_dsp);
      control_mass_sparsity_.copy_from(control_mass_dsp);
      control_mass_.reinit(control_mass_sparsity_);

      forcing_load_.reinit(state_size);
      desired_state_load_.reinit(state_size);
    }

    void
    assemble(const dealii::Function<dim> &forcing,
             const dealii::Function<dim> &desired_state)
    {
      const unsigned int quadrature_order =
        std::max(state_fe_.degree, control_fe_.degree) + 2;
      const dealii::QGauss<dim> quadrature(quadrature_order);
      dealii::FEValues<dim> state_values(
        state_fe_,
        quadrature,
        dealii::update_values | dealii::update_gradients |
          dealii::update_quadrature_points | dealii::update_JxW_values);
      dealii::FEValues<dim> control_values(control_fe_,
                                            quadrature,
                                            dealii::update_values);

      dealii::FullMatrix<double> local_system(state_fe_.dofs_per_cell,
                                              state_fe_.dofs_per_cell);
      dealii::FullMatrix<double> local_state_mass(state_fe_.dofs_per_cell,
                                                  state_fe_.dofs_per_cell);
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
      for (; state_cell != state_dof_handler_.end();
           ++state_cell, ++control_cell)
        {
          contract::require(control_cell != control_dof_handler_.end(),
                            "State and control DoF handlers do not share cells");
          state_values.reinit(state_cell);
          control_values.reinit(control_cell);
          local_system = 0.0;
          local_state_mass = 0.0;
          local_control_coupling = 0.0;
          local_control_mass = 0.0;
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
                    {
                      local_system(i, j) +=
                        (diffusion_ *
                           (state_values.shape_grad(i, q) *
                            state_values.shape_grad(j, q)) +
                         reaction_ * phi_i * state_values.shape_value(j, q)) *
                        weight;
                      local_state_mass(i, j) +=
                        phi_i * state_values.shape_value(j, q) * weight;
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

          state_cell->get_dof_indices(state_indices);
          control_cell->get_dof_indices(control_indices);

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
                  if (constrained_state_dofs_.at(global_j))
                    continue;
                  system_matrix_.add(global_i, global_j, local_system(i, j));
                  state_mass_.add(global_i, global_j, local_state_mass(i, j));
                }
              for (unsigned int j = 0; j < control_fe_.dofs_per_cell; ++j)
                control_coupling_.add(global_i,
                                      control_indices[j],
                                      local_control_coupling(i, j));
            }

          for (unsigned int i = 0; i < control_fe_.dofs_per_cell; ++i)
            for (unsigned int j = 0; j < control_fe_.dofs_per_cell; ++j)
              control_mass_.add(control_indices[i],
                                control_indices[j],
                                local_control_mass(i, j));
        }
      contract::require(control_cell == control_dof_handler_.end(),
                        "State and control DoF handlers have different cells");

      for (dealii::types::global_dof_index index = 0;
           index < state_dof_handler_.n_dofs();
           ++index)
        if (constrained_state_dofs_.at(index))
          system_matrix_.set(index, index, 1.0);
    }

    void
    solve_symmetric_system(Vector &solution, const Vector &right_hand_side) const
    {
      const double tolerance =
        std::max(1e-14, 1e-12 * right_hand_side.l2_norm());
      dealii::SolverControl control(
        std::max<unsigned int>(100, 10 * system_matrix_.m()), tolerance);
      dealii::SolverCG<Vector> solver(control);
      solver.solve(system_matrix_,
                   solution,
                   right_hand_side,
                   dealii::PreconditionIdentity());
    }

    dealii::FE_Q<dim> state_fe_;
    dealii::FE_DGQ<dim> control_fe_;
    dealii::DoFHandler<dim> state_dof_handler_;
    dealii::DoFHandler<dim> control_dof_handler_;
    dealii::AffineConstraints<double> state_constraints_;
    std::vector<bool> constrained_state_dofs_;

    const double diffusion_;
    const double reaction_;
    const double regularisation_weight_;
    const std::set<dealii::types::boundary_id> dirichlet_boundary_ids_;

    dealii::SparsityPattern state_sparsity_;
    dealii::SparsityPattern control_sparsity_;
    dealii::SparsityPattern control_mass_sparsity_;
    dealii::SparseMatrix<double> system_matrix_;
    dealii::SparseMatrix<double> state_mass_;
    dealii::SparseMatrix<double> control_coupling_;
    dealii::SparseMatrix<double> control_mass_;
    Vector forcing_load_;
    Vector desired_state_load_;
    double desired_state_norm_ = 0.0;

    contract::LayoutPtr variable_layout_;
    contract::LayoutPtr test_layout_;
    contract::LayoutPtr state_layout_;
    contract::LayoutPtr control_layout_;
  };
} // namespace nmopt::dealii_backend
