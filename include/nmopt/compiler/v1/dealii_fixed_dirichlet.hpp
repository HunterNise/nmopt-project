#pragma once

#include "nmopt/contract/executable_model.hpp"
#include "nmopt/dealii/cellwise_box_constraint.hpp"
#include "nmopt/dealii/mass_metric.hpp"
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
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparse_matrix.h>
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
  // This v1-only target represents the independent state coefficients
  // y_hat. It evaluates the compiled equation and observation/loss on
  //
  //   y_phys = P_h y_hat + ell_0,h,
  //
  // and pulls covectors back with P_h^*. The direct v0 lowerer is not used or
  // modified here: it remains the homogeneous comparison implementation.
  template <int dim>
  class AssembledScalarDiffusionReactionModel final
    : public contract::ExecutableModelT<dealii_backend::SerialBackend>
  {
  public:
    using Backend = dealii_backend::SerialBackend;
    using Vector = dealii::Vector<double>;
    using Primal = contract::PrimalBlockT<Backend>;
    using Covector = contract::CovectorBlockT<Backend>;

    AssembledScalarDiffusionReactionModel(
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
      std::set<dealii::types::material_id>          observation_material_ids = {})
      : state_fe_(state_degree)
      , control_fe_(0)
      , state_dof_handler_(triangulation)
      , control_dof_handler_(triangulation)
      , diffusion_(diffusion)
      , reaction_(reaction)
      , regularisation_weight_(regularisation_weight)
      , dirichlet_boundary_ids_(std::move(dirichlet_boundary_ids))
      , observation_material_ids_(std::move(observation_material_ids))
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
                        "The assembled v1 target needs a fixed Dirichlet boundary");

      state_dof_handler_.distribute_dofs(state_fe_);
      control_dof_handler_.distribute_dofs(control_fe_);
      build_constraints(fixed_dirichlet_data);
      build_reconstruction();
      initialise_storage();
      assemble_physical_operators(forcing, desired_state);
      assemble_reduced_solve_operators();
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
    control_l2_box_constraint(Vector lower, Vector upper) const
    {
      return dealii_backend::CellwiseBoxConstraint(control_layout_,
                                                    std::move(lower),
                                                    std::move(upper));
    }

    dealii_backend::CellwiseBoxConstraint
    control_l2_box_constraint(const double lower, const double upper) const
    {
      return dealii_backend::CellwiseBoxConstraint(control_layout_, lower, upper);
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
      Vector state_mass_times_state(state_dof_handler_.n_dofs());
      physical_state_mass_.vmult(state_mass_times_state, physical_state);
      const double state_value =
        0.5 * (physical_state * state_mass_times_state) -
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
      physical_state_mass_.vmult(physical_state, reconstruct(variables.block(0)));
      physical_state.add(-1.0, desired_state_load_);

      Vector control(control_dof_handler_.n_dofs());
      control_mass_->vmult(control, variables.block(1));
      control *= regularisation_weight_;
      return Covector(variable_layout_,
                      {pullback(physical_state), std::move(control)});
    }

    Primal
    solve_state(const Primal &control) const
    {
      contract::require(control.layout()->compatible_with(*control_layout_),
                        "State solve control has an incompatible layout");
      Vector right_hand_side = reduced_forcing_load_;
      Vector control_contribution(independent_state_dofs_.size());
      reduced_control_coupling_.vmult(control_contribution, control.block(0));
      right_hand_side.add(1.0, control_contribution);

      Vector state(independent_state_dofs_.size());
      solve_symmetric_system(reduced_system_matrix_, state, right_hand_side);
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

      Vector adjoint(independent_state_dofs_.size());
      solve_symmetric_system(reduced_system_matrix_,
                             adjoint,
                             state_objective_derivative.block(0));
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
      physical_state_mass_.reinit(physical_state_sparsity_);

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
    assemble_physical_operators(const dealii::Function<dim> &forcing,
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
      dealii::FEValues<dim> control_values(control_fe_, quadrature,
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
      for (; state_cell != state_dof_handler_.end(); ++state_cell, ++control_cell)
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
          const bool observe_cell =
            observation_material_ids_.empty() ||
            observation_material_ids_.find(state_cell->material_id()) !=
              observation_material_ids_.end();

          for (unsigned int q = 0; q < quadrature.size(); ++q)
            {
              const double weight = state_values.JxW(q);
              const double forcing_value = forcing.value(state_values.quadrature_point(q));
              const double desired_value = observe_cell
                                             ? desired_state.value(
                                                 state_values.quadrature_point(q))
                                             : 0.0;
              if (observe_cell)
                desired_state_norm_ += desired_value * desired_value * weight;

              for (unsigned int i = 0; i < state_fe_.dofs_per_cell; ++i)
                {
                  const double phi_i = state_values.shape_value(i, q);
                  local_forcing(i) += forcing_value * phi_i * weight;
                  if (observe_cell)
                    local_desired_state(i) += desired_value * phi_i * weight;
                  for (unsigned int j = 0; j < state_fe_.dofs_per_cell; ++j)
                    {
                      local_system(i, j) +=
                        (diffusion_ * (state_values.shape_grad(i, q) *
                                       state_values.shape_grad(j, q)) +
                         reaction_ * phi_i * state_values.shape_value(j, q)) *
                        weight;
                      if (observe_cell)
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
              forcing_load_[state_indices[i]] += local_forcing(i);
              desired_state_load_[state_indices[i]] += local_desired_state(i);
              for (unsigned int j = 0; j < state_fe_.dofs_per_cell; ++j)
                {
                  physical_system_matrix_.add(state_indices[i],
                                              state_indices[j],
                                              local_system(i, j));
                  physical_state_mass_.add(state_indices[i],
                                           state_indices[j],
                                           local_state_mass(i, j));
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

    static void
    solve_symmetric_system(const dealii::SparseMatrix<double> &matrix,
                           Vector &                              solution,
                           const Vector &                        right_hand_side)
    {
      const double tolerance =
        std::max(1e-14, 1e-12 * right_hand_side.l2_norm());
      dealii::SolverControl control(
        std::max<unsigned int>(100, 10 * matrix.m()), tolerance);
      dealii::SolverCG<Vector> solver(control);
      solver.solve(matrix,
                   solution,
                   right_hand_side,
                   dealii::PreconditionIdentity());
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
    const std::set<dealii::types::boundary_id> dirichlet_boundary_ids_;
    const std::set<dealii::types::material_id> observation_material_ids_;

    dealii::SparsityPattern reconstruction_sparsity_;
    dealii::SparseMatrix<double> reconstruction_;
    Vector lifting_;

    dealii::SparsityPattern physical_state_sparsity_;
    dealii::SparsityPattern physical_control_sparsity_;
    dealii::SparsityPattern control_mass_sparsity_;
    dealii::SparseMatrix<double> physical_system_matrix_;
    dealii::SparseMatrix<double> physical_state_mass_;
    dealii::SparseMatrix<double> physical_control_coupling_;
    std::shared_ptr<dealii::SparseMatrix<double>> control_mass_;
    Vector forcing_load_;
    Vector desired_state_load_;
    double desired_state_norm_ = 0.0;

    dealii::SparsityPattern reduced_state_sparsity_;
    dealii::SparsityPattern reduced_control_sparsity_;
    dealii::SparseMatrix<double> reduced_system_matrix_;
    dealii::SparseMatrix<double> reduced_control_coupling_;
    Vector reduced_forcing_load_;

    contract::LayoutPtr variable_layout_;
    contract::LayoutPtr test_layout_;
    contract::LayoutPtr state_layout_;
    contract::LayoutPtr control_layout_;
  };
} // namespace nmopt::compiler::v1::detail
