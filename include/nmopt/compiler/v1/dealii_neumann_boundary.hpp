#pragma once

#include "nmopt/contract/executable_model.hpp"
#include "nmopt/dealii/facewise_box_constraint.hpp"
#include "nmopt/dealii/mass_metric.hpp"
#include "nmopt/dealii/serial_backend.hpp"

#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparse_direct.h>
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
  // V1-only Neumann target. It owns the facewise control layout and lowers
  // exactly the declared natural-boundary contribution
  //
  //   r(y,u) = A y - f_h - C_Gamma u,
  //
  // where each column of C_Gamma is the FEFaceValues realization of one
  // marked boundary face. It intentionally does not reuse the volume-control
  // model or present this term as a Dirichlet/control lifting.
  template <int dim>
  class NeumannBoundaryControlModel final
    : public contract::ExecutableModelT<dealii_backend::SerialBackend>
  {
  public:
    using Backend = dealii_backend::SerialBackend;
    using Vector = dealii::Vector<double>;
    using Primal = contract::PrimalBlockT<Backend>;
    using Covector = contract::CovectorBlockT<Backend>;

    enum class StateGauge
    {
      fixed_dirichlet,
      mean_zero_multiplier
    };

    NeumannBoundaryControlModel(
      dealii::Triangulation<dim> &                triangulation,
      const dealii::Function<dim> &               forcing,
      const dealii::Function<dim> &               desired_state,
      const double                                  diffusion,
      const double                                  reaction,
      const double                                  regularisation_weight,
      const unsigned int                            state_degree,
      std::set<dealii::types::boundary_id>         dirichlet_boundary_ids,
      std::set<dealii::types::boundary_id>         control_boundary_ids,
      std::set<dealii::types::boundary_id>         observation_boundary_ids,
      const StateGauge                              state_gauge =
        StateGauge::fixed_dirichlet)
      : state_fe_(state_degree)
      , state_dof_handler_(triangulation)
      , diffusion_(diffusion)
      , reaction_(reaction)
      , regularisation_weight_(regularisation_weight)
      , state_gauge_(state_gauge)
      , dirichlet_boundary_ids_(std::move(dirichlet_boundary_ids))
      , control_boundary_ids_(std::move(control_boundary_ids))
      , observation_boundary_ids_(std::move(observation_boundary_ids))
    {
      contract::require(diffusion_ > 0.0,
                        "Diffusion coefficient must be strictly positive");
      contract::require(reaction_ >= 0.0,
                        "Reaction coefficient must be non-negative");
      contract::require(regularisation_weight_ > 0.0,
                        "Control regularisation weight must be strictly positive");
      contract::require(state_degree > 0,
                        "State FE degree must be at least one");
      contract::require(
        state_gauge_ != StateGauge::fixed_dirichlet ||
          !dirichlet_boundary_ids_.empty(),
        "The Neumann v1 target needs a fixed Dirichlet boundary for this gauge");
      contract::require(
        state_gauge_ != StateGauge::mean_zero_multiplier ||
          dirichlet_boundary_ids_.empty(),
        "The pure-Neumann mean-zero gauge cannot also fix boundary DoFs");
      contract::require(
        state_gauge_ != StateGauge::mean_zero_multiplier ||
          std::abs(reaction_) <= 1e-14,
        "The pure-Neumann mean-zero gauge requires zero reaction");
      contract::require(!control_boundary_ids_.empty(),
                        "The Neumann v1 target needs a marked control boundary");
      contract::require(!observation_boundary_ids_.empty(),
                        "The Neumann v1 target needs a marked observation boundary");

      state_dof_handler_.distribute_dofs(state_fe_);
      build_constraints();
      control_face_count_ = count_control_faces();
      contract::require(control_face_count_ > 0,
                        "The selected Neumann boundary has no active boundary faces");
      initialise_storage();
      assemble(forcing, desired_state);
      if (uses_mean_zero_gauge())
        build_mean_zero_system();
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
    control_l2_metric(
      dealii_backend::MassMetricSolveParameters solve_parameters = {}) const
    {
      return dealii_backend::MassMetric("l2_facewise",
                                        control_layout_,
                                        control_mass_,
                                        solve_parameters);
    }

    dealii_backend::FacewiseBoxConstraint
    control_l2_box_constraint(
      Vector                              lower,
      Vector                              upper,
      const dealii_backend::MassMetric & projection_metric) const
    {
      return dealii_backend::FacewiseBoxConstraint(control_layout_,
                                                    std::move(lower),
                                                    std::move(upper),
                                                    projection_metric);
    }

    dealii_backend::FacewiseBoxConstraint
    control_l2_box_constraint(
      const double                        lower,
      const double                        upper,
      const dealii_backend::MassMetric & projection_metric) const
    {
      return dealii_backend::FacewiseBoxConstraint(control_layout_,
                                                    lower,
                                                    upper,
                                                    projection_metric);
    }

    bool
    uses_mean_zero_gauge() const
    {
      return state_gauge_ == StateGauge::mean_zero_multiplier;
    }

    bool
    forcing_is_compatible() const
    {
      return is_compatible(forcing_load_);
    }

    double
    state_mean(const Primal &state) const
    {
      contract::require(
        state.layout()->compatible_with(*state_layout_) ||
          state.layout()->compatible_with(*test_layout_),
        "State mean received an incompatible state or adjoint layout");
      return mean_constraint_ * state.block(0);
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
      control_coupling_.vmult(control_contribution, variable_tangent.block(1));
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
      Vector control(control_face_count_);
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

      Vector control_mass_times_control(control_face_count_);
      control_mass_->vmult(control_mass_times_control, variables.block(1));
      const double control_value = 0.5 * regularisation_weight_ *
                                   (variables.block(1) * control_mass_times_control);
      return state_value + control_value;
    }

    Covector
    objective_derivative(const Primal &variables) const override
    {
      require_variables(variables, "Objective derivative");
      Vector state(state_dof_handler_.n_dofs());
      state_mass_.vmult(state, variables.block(0));
      state.add(-1.0, desired_state_load_);

      Vector control(control_face_count_);
      control_mass_->vmult(control, variables.block(1));
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
      if (uses_mean_zero_gauge())
        {
          contract::require(
            is_compatible(right_hand_side),
            "Pure-Neumann state load violates the discrete constant-mode compatibility condition");
          solve_mean_zero_system(state, right_hand_side);
        }
      else
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

      Vector adjoint(test_layout_->dimension(0));
      if (uses_mean_zero_gauge())
        solve_mean_zero_system(adjoint, state_objective_derivative.block(0));
      else
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

    bool
    is_control_face(const typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
                    const unsigned int face) const
    {
      return cell->face(face)->at_boundary() &&
             control_boundary_ids_.count(cell->face(face)->boundary_id()) != 0;
    }

    bool
    is_observation_face(
      const typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
      const unsigned int face) const
    {
      return cell->face(face)->at_boundary() &&
             observation_boundary_ids_.count(cell->face(face)->boundary_id()) != 0;
    }

    std::size_t
    count_control_faces() const
    {
      std::size_t result = 0;
      for (auto cell = state_dof_handler_.begin_active();
           cell != state_dof_handler_.end();
           ++cell)
        for (unsigned int face = 0; face < dealii::GeometryInfo<dim>::faces_per_cell;
             ++face)
          result += is_control_face(cell, face);
      return result;
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
            "The Neumann v1 target does not support hanging, periodic, or "
            "other non-Dirichlet affine constraints");
    }

    void
    initialise_storage()
    {
      const auto state_size = state_dof_handler_.n_dofs();
      variable_layout_ = std::make_shared<const contract::BlockLayout>(
        "neumann_boundary_variables",
        std::vector<contract::SpaceId>{{"state"}, {"control"}},
        std::vector<std::size_t>{state_size, control_face_count_});
      test_layout_ = std::make_shared<const contract::BlockLayout>(
        "neumann_boundary_state_test",
        std::vector<contract::SpaceId>{{"state_test"}},
        std::vector<std::size_t>{state_size});
      state_layout_ = variable_layout_->single_block(0, "state");
      control_layout_ = variable_layout_->single_block(1, "control");

      dealii::DynamicSparsityPattern state_dsp(state_size, state_size);
      dealii::DoFTools::make_sparsity_pattern(state_dof_handler_, state_dsp);
      state_sparsity_.copy_from(state_dsp);
      system_matrix_.reinit(state_sparsity_);
      state_mass_.reinit(state_sparsity_);

      dealii::DynamicSparsityPattern control_dsp(state_size, control_face_count_);
      std::vector<dealii::types::global_dof_index> state_indices(
        state_fe_.dofs_per_cell);
      std::size_t control_index = 0;
      for (auto cell = state_dof_handler_.begin_active();
           cell != state_dof_handler_.end();
           ++cell)
        {
          cell->get_dof_indices(state_indices);
          for (unsigned int face = 0; face < dealii::GeometryInfo<dim>::faces_per_cell;
               ++face)
            if (is_control_face(cell, face))
              {
                for (const auto state_index : state_indices)
                  if (!constrained_state_dofs_.at(state_index))
                    control_dsp.add(state_index, control_index);
                ++control_index;
              }
        }
      contract::require(control_index == control_face_count_,
                        "Neumann control face enumeration changed during setup");
      control_sparsity_.copy_from(control_dsp);
      control_coupling_.reinit(control_sparsity_);

      dealii::DynamicSparsityPattern control_mass_dsp(control_face_count_,
                                                       control_face_count_);
      for (std::size_t index = 0; index < control_face_count_; ++index)
        control_mass_dsp.add(index, index);
      control_mass_sparsity_.copy_from(control_mass_dsp);
      control_mass_ = std::make_shared<dealii::SparseMatrix<double>>();
      control_mass_->reinit(control_mass_sparsity_);

      forcing_load_.reinit(state_size);
      desired_state_load_.reinit(state_size);
      mean_constraint_.reinit(state_size);
      constant_mode_.reinit(state_size);
      for (dealii::types::global_dof_index index = 0; index < state_size; ++index)
        constant_mode_[index] = 1.0;
    }

    void
    assemble(const dealii::Function<dim> &forcing,
             const dealii::Function<dim> &desired_state)
    {
      const unsigned int quadrature_order = state_fe_.degree + 2;
      const dealii::QGauss<dim> volume_quadrature(quadrature_order);
      dealii::FEValues<dim> state_values(
        state_fe_,
        volume_quadrature,
        dealii::update_values | dealii::update_gradients |
          dealii::update_quadrature_points | dealii::update_JxW_values);
      dealii::FullMatrix<double> local_system(state_fe_.dofs_per_cell,
                                              state_fe_.dofs_per_cell);
      dealii::Vector<double> local_forcing(state_fe_.dofs_per_cell);
      dealii::Vector<double> local_mean_constraint(state_fe_.dofs_per_cell);
      std::vector<dealii::types::global_dof_index> state_indices(
        state_fe_.dofs_per_cell);

      for (auto cell = state_dof_handler_.begin_active();
           cell != state_dof_handler_.end();
           ++cell)
        {
          state_values.reinit(cell);
          local_system = 0.0;
          local_forcing = 0.0;
          local_mean_constraint = 0.0;
          for (unsigned int q = 0; q < volume_quadrature.size(); ++q)
            {
              const double weight = state_values.JxW(q);
              const double forcing_value =
                forcing.value(state_values.quadrature_point(q));
              for (unsigned int i = 0; i < state_fe_.dofs_per_cell; ++i)
                {
                  const double phi_i = state_values.shape_value(i, q);
                  local_forcing(i) += forcing_value * phi_i * weight;
                  local_mean_constraint(i) += phi_i * weight;
                  for (unsigned int j = 0; j < state_fe_.dofs_per_cell; ++j)
                    local_system(i, j) +=
                      (diffusion_ * (state_values.shape_grad(i, q) *
                                     state_values.shape_grad(j, q)) +
                       reaction_ * phi_i * state_values.shape_value(j, q)) *
                      weight;
                }
            }
          cell->get_dof_indices(state_indices);
          for (unsigned int i = 0; i < state_fe_.dofs_per_cell; ++i)
            {
              const auto global_i = state_indices[i];
              if (constrained_state_dofs_.at(global_i))
                continue;
              forcing_load_[global_i] += local_forcing(i);
              mean_constraint_[global_i] += local_mean_constraint(i);
              for (unsigned int j = 0; j < state_fe_.dofs_per_cell; ++j)
                {
                  const auto global_j = state_indices[j];
                  if (!constrained_state_dofs_.at(global_j))
                    system_matrix_.add(global_i, global_j, local_system(i, j));
                }
            }
        }

      if (!uses_mean_zero_gauge())
        for (dealii::types::global_dof_index index = 0;
             index < state_dof_handler_.n_dofs();
             ++index)
          if (constrained_state_dofs_.at(index))
            system_matrix_.set(index, index, 1.0);

      assemble_boundary_operators(desired_state, quadrature_order);
    }

    void
    assemble_boundary_operators(const dealii::Function<dim> &desired_state,
                                const unsigned int quadrature_order)
    {
      const dealii::QGauss<dim - 1> face_quadrature(quadrature_order);
      dealii::FEFaceValues<dim> face_values(
        state_fe_,
        face_quadrature,
        dealii::update_values | dealii::update_quadrature_points |
          dealii::update_JxW_values);
      std::vector<dealii::types::global_dof_index> state_indices(
        state_fe_.dofs_per_cell);
      std::size_t control_index = 0;
      for (auto cell = state_dof_handler_.begin_active();
           cell != state_dof_handler_.end();
           ++cell)
        {
          cell->get_dof_indices(state_indices);
          for (unsigned int face = 0; face < dealii::GeometryInfo<dim>::faces_per_cell;
               ++face)
            {
              const bool control_face = is_control_face(cell, face);
              const bool observation_face = is_observation_face(cell, face);
              if (!control_face && !observation_face)
                continue;
              face_values.reinit(cell, face);
              double control_measure = 0.0;
              for (unsigned int q = 0; q < face_quadrature.size(); ++q)
                {
                  const double weight = face_values.JxW(q);
                  if (control_face)
                    control_measure += weight;
                  const double desired_value = observation_face
                    ? desired_state.value(face_values.quadrature_point(q))
                    : 0.0;
                  if (observation_face)
                    desired_state_norm_ += desired_value * desired_value * weight;
                  for (unsigned int i = 0; i < state_fe_.dofs_per_cell; ++i)
                    {
                      const auto global_i = state_indices[i];
                      if (constrained_state_dofs_.at(global_i))
                        continue;
                      const double phi_i = face_values.shape_value(i, q);
                      if (control_face)
                        control_coupling_.add(global_i,
                                              control_index,
                                              phi_i * weight);
                      if (observation_face)
                        {
                          desired_state_load_[global_i] +=
                            desired_value * phi_i * weight;
                          for (unsigned int j = 0;
                               j < state_fe_.dofs_per_cell;
                               ++j)
                            {
                              const auto global_j = state_indices[j];
                              if (!constrained_state_dofs_.at(global_j))
                                state_mass_.add(
                                  global_i,
                                  global_j,
                                  phi_i * face_values.shape_value(j, q) * weight);
                            }
                        }
                    }
                }
              if (control_face)
                {
                  contract::require(control_measure > 0.0,
                                    "A Neumann control face has zero measure");
                  control_mass_->add(control_index, control_index, control_measure);
                  ++control_index;
                }
            }
        }
      contract::require(control_index == control_face_count_,
                        "Neumann control face enumeration changed during assembly");
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

    bool
    is_compatible(const Vector &load) const
    {
      const double tolerance = std::max(1e-12, 1e-11 * load.l2_norm());
      return std::abs(constant_mode_ * load) <= tolerance;
    }

    void
    build_mean_zero_system()
    {
      const auto state_size = state_dof_handler_.n_dofs();
      dealii::DynamicSparsityPattern augmented_dsp(state_size + 1,
                                                    state_size + 1);
      for (dealii::types::global_dof_index row = 0; row < state_size; ++row)
        {
          for (auto entry = system_matrix_.begin(row);
               entry != system_matrix_.end(row);
               ++entry)
            augmented_dsp.add(row, entry->column());
          if (mean_constraint_[row] != 0.0)
            {
              augmented_dsp.add(row, state_size);
              augmented_dsp.add(state_size, row);
            }
        }
      augmented_sparsity_.copy_from(augmented_dsp);
      augmented_system_.reinit(augmented_sparsity_);
      for (dealii::types::global_dof_index row = 0; row < state_size; ++row)
        {
          for (auto entry = system_matrix_.begin(row);
               entry != system_matrix_.end(row);
               ++entry)
            augmented_system_.set(row, entry->column(), entry->value());
          augmented_system_.set(row, state_size, mean_constraint_[row]);
          augmented_system_.set(state_size, row, mean_constraint_[row]);
        }
      augmented_solver_.initialize(augmented_system_);
    }

    void
    solve_mean_zero_system(Vector &solution, const Vector &right_hand_side) const
    {
      const auto state_size = state_dof_handler_.n_dofs();
      dealii::Vector<double> augmented_right_hand_side(state_size + 1);
      for (dealii::types::global_dof_index index = 0; index < state_size; ++index)
        augmented_right_hand_side[index] = right_hand_side[index];
      dealii::Vector<double> augmented_solution(state_size + 1);
      augmented_solver_.vmult(augmented_solution, augmented_right_hand_side);
      for (dealii::types::global_dof_index index = 0; index < state_size; ++index)
        solution[index] = augmented_solution[index];
    }

    dealii::FE_Q<dim> state_fe_;
    dealii::DoFHandler<dim> state_dof_handler_;
    dealii::AffineConstraints<double> state_constraints_;
    std::vector<bool> constrained_state_dofs_;

    const double diffusion_;
    const double reaction_;
    const double regularisation_weight_;
    const StateGauge state_gauge_;
    const std::set<dealii::types::boundary_id> dirichlet_boundary_ids_;
    const std::set<dealii::types::boundary_id> control_boundary_ids_;
    const std::set<dealii::types::boundary_id> observation_boundary_ids_;
    std::size_t control_face_count_ = 0;

    dealii::SparsityPattern state_sparsity_;
    dealii::SparsityPattern control_sparsity_;
    dealii::SparsityPattern control_mass_sparsity_;
    dealii::SparsityPattern augmented_sparsity_;
    dealii::SparseMatrix<double> system_matrix_;
    dealii::SparseMatrix<double> state_mass_;
    dealii::SparseMatrix<double> control_coupling_;
    dealii::SparseMatrix<double> augmented_system_;
    dealii::SparseDirectUMFPACK augmented_solver_;
    std::shared_ptr<dealii::SparseMatrix<double>> control_mass_;
    Vector forcing_load_;
    Vector desired_state_load_;
    Vector mean_constraint_;
    Vector constant_mode_;
    double desired_state_norm_ = 0.0;

    contract::LayoutPtr variable_layout_;
    contract::LayoutPtr test_layout_;
    contract::LayoutPtr state_layout_;
    contract::LayoutPtr control_layout_;
  };
} // namespace nmopt::compiler::v1::detail
