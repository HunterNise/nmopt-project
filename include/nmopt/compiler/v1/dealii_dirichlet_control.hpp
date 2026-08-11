#pragma once

#include "nmopt/contract/executable_model.hpp"
#include "nmopt/dealii/mass_metric.hpp"
#include "nmopt/dealii/serial_backend.hpp"
#include "nmopt/dealii/serial_spd_solver.hpp"

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
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
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace nmopt::compiler::v1::detail
{
  // This v1-only target keeps the independent state coordinates y_hat and
  // makes the control-to-trace map explicit:
  //
  //   y_phys = P_h y_hat + L_D,h u_h.
  //
  // The first policy gives one shared nodal coefficient to every state DoF on
  // the complete exterior boundary. It is intentionally not a boundary-load
  // residual and does not alter the direct v0 lowerer.
  template <int dim>
  class DirichletControlLiftingModel final
    : public contract::ExecutableModelT<dealii_backend::SerialBackend>
  {
  public:
    using Backend = dealii_backend::SerialBackend;
    using Vector = dealii::Vector<double>;
    using Primal = contract::PrimalBlockT<Backend>;
    using Covector = contract::CovectorBlockT<Backend>;
    using SolveResult = contract::FormulationSolveResultT<Backend>;

    DirichletControlLiftingModel(
      dealii::Triangulation<dim> &       triangulation,
      const dealii::Function<dim> &      forcing,
      const dealii::Function<dim> &      desired_state,
      const double                         diffusion,
      const double                         reaction,
      const double                         regularisation_weight,
      const unsigned int                   state_degree,
      std::set<dealii::types::boundary_id> controlled_boundary_ids,
      std::set<dealii::types::boundary_id> fixed_boundary_ids = {},
      std::optional<std::reference_wrapper<const dealii::Function<dim>>>
        fixed_dirichlet_data = std::nullopt)
      : state_fe_(state_degree)
      , state_dof_handler_(triangulation)
      , diffusion_(diffusion)
      , reaction_(reaction)
      , regularisation_weight_(regularisation_weight)
      , controlled_boundary_ids_(std::move(controlled_boundary_ids))
      , fixed_boundary_ids_(std::move(fixed_boundary_ids))
    {
      contract::require(diffusion_ > 0.0,
                        "Diffusion coefficient must be strictly positive");
      contract::require(reaction_ >= 0.0,
                        "Reaction coefficient must be non-negative");
      contract::require(regularisation_weight_ > 0.0,
                        "Control regularisation weight must be strictly positive");
      contract::require(state_degree > 0,
                        "State FE degree must be at least one");
      contract::require(!controlled_boundary_ids_.empty(),
                        "The Dirichlet-control target needs a controlled boundary");
      contract::require(
        fixed_boundary_ids_.empty() == !fixed_dirichlet_data.has_value(),
        "Partial Dirichlet control must bind fixed data exactly when it declares a fixed boundary");

      state_dof_handler_.distribute_dofs(state_fe_);
      build_control_dof_map();
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
      return dealii_backend::MassMetric("l2_dirichlet_trace",
                                        control_layout_,
                                        control_boundary_mass_,
                                        solve_parameters);
    }

    Vector
    reconstruct_physical_state(const Primal &variables) const
    {
      require_variables(variables, "Physical-state reconstruction");
      return reconstruct(variables.block(0), variables.block(1));
    }

    Covector
    discrete_conormal_covector(const Primal &variables,
                               const Primal &adjoint) const
    {
      require_variables(variables, "Discrete conormal");
      contract::require(adjoint.layout()->compatible_with(*test_layout_),
                        "Discrete conormal received an incompatible adjoint layout");

      Vector physical_adjoint_action(state_dof_handler_.n_dofs());
      physical_system_matrix_.Tvmult(physical_adjoint_action,
                                     embed_state(adjoint.block(0)));
      Vector tracking_covector(state_dof_handler_.n_dofs());
      physical_state_mass_.vmult(
        tracking_covector,
        reconstruct(variables.block(0), variables.block(1)));
      tracking_covector.add(-1.0, desired_state_load_);
      physical_adjoint_action.add(-1.0, tracking_covector);

      return Covector(control_layout_,
                      {pullback_control(physical_adjoint_action)});
    }

    Covector
    residual(const Primal &variables) const override
    {
      require_variables(variables, "Residual");
      Vector physical_value(state_dof_handler_.n_dofs());
      physical_system_matrix_.vmult(
        physical_value, reconstruct(variables.block(0), variables.block(1)));
      physical_value.add(-1.0, forcing_load_);
      return Covector(test_layout_, {pullback_state(physical_value)});
    }

    Covector
    residual_jvp(const Primal &variables,
                 const Primal &variable_tangent) const override
    {
      require_variables(variables, "Residual JVP");
      require_variables(variable_tangent, "Residual JVP tangent");
      Vector physical_value(state_dof_handler_.n_dofs());
      physical_system_matrix_.vmult(
        physical_value,
        physical_tangent(variable_tangent.block(0), variable_tangent.block(1)));
      return Covector(test_layout_, {pullback_state(physical_value)});
    }

    Covector
    residual_vjp(const Primal &variables,
                 const Primal &test_seed) const override
    {
      require_variables(variables, "Residual VJP");
      contract::require(test_seed.layout()->compatible_with(*test_layout_),
                        "Residual VJP seed has an incompatible test layout");

      Vector physical_covector(state_dof_handler_.n_dofs());
      physical_system_matrix_.Tvmult(
        physical_covector, embed_state(test_seed.block(0)));
      return Covector(variable_layout_,
                      {pullback_state(physical_covector),
                       pullback_control(physical_covector)});
    }

    double
    objective(const Primal &variables) const override
    {
      require_variables(variables, "Objective");
      const Vector physical_state =
        reconstruct(variables.block(0), variables.block(1));
      Vector state_mass_times_state(state_dof_handler_.n_dofs());
      physical_state_mass_.vmult(state_mass_times_state, physical_state);
      const double state_value =
        0.5 * (physical_state * state_mass_times_state) -
        (desired_state_load_ * physical_state) + 0.5 * desired_state_norm_;

      Vector control_mass_times_control(controlled_state_dofs_.size());
      control_boundary_mass_->vmult(control_mass_times_control,
                                     variables.block(1));
      const double control_value = 0.5 * regularisation_weight_ *
                                   (variables.block(1) *
                                    control_mass_times_control);
      return state_value + control_value;
    }

    Covector
    objective_derivative(const Primal &variables) const override
    {
      require_variables(variables, "Objective derivative");
      Vector physical_covector(state_dof_handler_.n_dofs());
      physical_state_mass_.vmult(
        physical_covector,
        reconstruct(variables.block(0), variables.block(1)));
      physical_covector.add(-1.0, desired_state_load_);

      Vector control = pullback_control(physical_covector);
      Vector regularisation(controlled_state_dofs_.size());
      control_boundary_mass_->vmult(regularisation, variables.block(1));
      control.add(regularisation_weight_, regularisation);
      return Covector(variable_layout_,
                      {pullback_state(physical_covector), std::move(control)});
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
      Vector lifting_contribution(independent_state_dofs_.size());
      reduced_lifting_coupling_.vmult(lifting_contribution, control.block(0));
      right_hand_side.add(-1.0, lifting_contribution);

      Vector state(independent_state_dofs_.size());
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
      auto report = solve_symmetric_system(reduced_system_matrix_,
                                           adjoint,
                                           state_objective_derivative.block(0),
                                           policy);
      return {Primal(test_layout_, {std::move(adjoint)}), std::move(report)};
    }

  private:
    static constexpr std::size_t invalid_control_index =
      std::numeric_limits<std::size_t>::max();

    void
    require_variables(const Primal &variables, const char *operation) const
    {
      contract::require(
        variables.layout()->compatible_with(*variable_layout_),
        std::string(operation) + " received an incompatible variable layout");
    }

    void
    build_control_dof_map()
    {
      const auto controlled_boundary_dofs = dealii::DoFTools::extract_boundary_dofs(
        state_dof_handler_, dealii::ComponentMask(), controlled_boundary_ids_);
      contract::require(!controlled_boundary_dofs.is_empty(),
                        "The controlled Dirichlet boundary has no state DoFs");
      fixed_state_dofs_.assign(state_dof_handler_.n_dofs(), false);
      // deal.II interprets an empty boundary-id set as "all boundary ids".
      // The complete-control target deliberately has no fixed boundary, so
      // do not query fixed DoFs until the partial target actually declares
      // such a region.
      if (!fixed_boundary_ids_.empty())
        {
          const auto fixed_boundary_dofs =
            dealii::DoFTools::extract_boundary_dofs(
              state_dof_handler_, dealii::ComponentMask(), fixed_boundary_ids_);
          for (auto iterator = fixed_boundary_dofs.begin();
               iterator != fixed_boundary_dofs.end(); ++iterator)
            fixed_state_dofs_[*iterator] = true;
        }
      control_index_for_state_dof_.assign(state_dof_handler_.n_dofs(),
                                           invalid_control_index);
      for (auto iterator = controlled_boundary_dofs.begin();
           iterator != controlled_boundary_dofs.end(); ++iterator)
        {
          const auto state_dof = *iterator;
          // The selected P5.4 interface policy gives fixed data precedence
          // at every fixed/controlled corner or interface DoF. The control
          // is therefore the relative-interior nodal trace and its endpoint
          // values are supplied by ell_0,h rather than implicit averaging.
          if (fixed_state_dofs_[state_dof])
            continue;
          control_index_for_state_dof_[state_dof] =
            controlled_state_dofs_.size();
          controlled_state_dofs_.push_back(state_dof);
        }
      contract::require(!controlled_state_dofs_.empty(),
                        "The selected controlled boundary has no independent trace DoFs after the fixed-interface policy");
    }

    void
    build_constraints(
      const std::optional<std::reference_wrapper<const dealii::Function<dim>>>
        fixed_dirichlet_data)
    {
      homogeneous_constraints_.clear();
      dealii::DoFTools::make_hanging_node_constraints(state_dof_handler_,
                                                       homogeneous_constraints_);
      for (const auto state_dof : controlled_state_dofs_)
        homogeneous_constraints_.add_line(state_dof);
      for (dealii::types::global_dof_index state_dof = 0;
           state_dof < fixed_state_dofs_.size(); ++state_dof)
        if (fixed_state_dofs_[state_dof])
          homogeneous_constraints_.add_line(state_dof);
      homogeneous_constraints_.close();

      // The selected initial policy excludes hanging/periodic/interface
      // relations. A controlled trace DoF must be one exact constrained row,
      // so L_D,h has no implicit averaging or corner choice.
      contract::require(
        homogeneous_constraints_.n_constraints() ==
          controlled_state_dofs_.size() +
            std::count(fixed_state_dofs_.begin(), fixed_state_dofs_.end(), true),
        "The Dirichlet-control lifting supports only hanging-free nodal traces");
      for (const auto state_dof : controlled_state_dofs_)
        contract::require(homogeneous_constraints_.is_constrained(state_dof),
                        "Each controlled trace DoF must be constrained");

      fixed_lifting_.reinit(state_dof_handler_.n_dofs());
      if (fixed_dirichlet_data)
        {
          std::map<dealii::types::global_dof_index, double> fixed_values;
          for (const auto boundary_id : fixed_boundary_ids_)
            dealii::VectorTools::interpolate_boundary_values(
              state_dof_handler_, boundary_id, fixed_dirichlet_data->get(),
              fixed_values);
          for (dealii::types::global_dof_index state_dof = 0;
               state_dof < fixed_state_dofs_.size(); ++state_dof)
            if (fixed_state_dofs_[state_dof])
              {
                const auto value = fixed_values.find(state_dof);
                contract::require(value != fixed_values.end(),
                                  "Every fixed Dirichlet trace DoF needs interpolated lifting data");
                fixed_lifting_[state_dof] = value->second;
              }
        }
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
    }

    void
    initialise_storage()
    {
      const auto physical_size = state_dof_handler_.n_dofs();
      const auto independent_size = independent_state_dofs_.size();
      const auto control_size = controlled_state_dofs_.size();
      variable_layout_ = std::make_shared<const contract::BlockLayout>(
        "dirichlet_control_lifted_variables",
        std::vector<contract::SpaceId>{{"state"}, {"control"}},
        std::vector<std::size_t>{independent_size, control_size});
      test_layout_ = std::make_shared<const contract::BlockLayout>(
        "dirichlet_control_lifted_state_test",
        std::vector<contract::SpaceId>{{"state_test"}},
        std::vector<std::size_t>{independent_size});
      state_layout_ = variable_layout_->single_block(0, "state");
      control_layout_ = variable_layout_->single_block(1, "control");

      dealii::DynamicSparsityPattern physical_dsp(physical_size, physical_size);
      dealii::DoFTools::make_sparsity_pattern(state_dof_handler_, physical_dsp);
      physical_sparsity_.copy_from(physical_dsp);
      physical_system_matrix_.reinit(physical_sparsity_);
      physical_state_mass_.reinit(physical_sparsity_);

      dealii::DynamicSparsityPattern boundary_mass_dsp(control_size, control_size);
      std::vector<dealii::types::global_dof_index> state_indices(
        state_fe_.dofs_per_cell);
      for (auto cell = state_dof_handler_.begin_active();
           cell != state_dof_handler_.end();
           ++cell)
        {
          cell->get_dof_indices(state_indices);
          for (unsigned int face = 0;
               face < dealii::GeometryInfo<dim>::faces_per_cell;
               ++face)
            if (cell->face(face)->at_boundary() &&
                controlled_boundary_ids_.count(cell->face(face)->boundary_id()) != 0)
              for (const auto state_i : state_indices)
                {
                  const auto control_i = control_index_for_state_dof_[state_i];
                  if (control_i == invalid_control_index)
                    continue;
                  for (const auto state_j : state_indices)
                    {
                      const auto control_j =
                        control_index_for_state_dof_[state_j];
                      if (control_j != invalid_control_index)
                        boundary_mass_dsp.add(control_i, control_j);
                    }
                }
        }
      control_boundary_mass_sparsity_.copy_from(boundary_mass_dsp);
      control_boundary_mass_ = std::make_shared<dealii::SparseMatrix<double>>();
      control_boundary_mass_->reinit(control_boundary_mass_sparsity_);

      forcing_load_.reinit(physical_size);
      desired_state_load_.reinit(physical_size);
      reduced_forcing_load_.reinit(independent_size);
    }

    void
    assemble_physical_operators(const dealii::Function<dim> &forcing,
                                const dealii::Function<dim> &desired_state)
    {
      const dealii::QGauss<dim> quadrature(state_fe_.degree + 2);
      const dealii::QGauss<dim - 1> face_quadrature(state_fe_.degree + 2);
      dealii::FEValues<dim> state_values(
        state_fe_, quadrature,
        dealii::update_values | dealii::update_gradients |
          dealii::update_quadrature_points | dealii::update_JxW_values);
      dealii::FEFaceValues<dim> face_values(
        state_fe_, face_quadrature,
        dealii::update_values | dealii::update_JxW_values);

      dealii::FullMatrix<double> local_system(state_fe_.dofs_per_cell,
                                              state_fe_.dofs_per_cell);
      dealii::FullMatrix<double> local_state_mass(state_fe_.dofs_per_cell,
                                                  state_fe_.dofs_per_cell);
      dealii::FullMatrix<double> local_boundary_mass(state_fe_.dofs_per_cell,
                                                     state_fe_.dofs_per_cell);
      dealii::Vector<double> local_forcing(state_fe_.dofs_per_cell);
      dealii::Vector<double> local_desired_state(state_fe_.dofs_per_cell);
      std::vector<dealii::types::global_dof_index> state_indices(
        state_fe_.dofs_per_cell);

      for (auto cell = state_dof_handler_.begin_active();
           cell != state_dof_handler_.end();
           ++cell)
        {
          state_values.reinit(cell);
          local_system = 0.0;
          local_state_mass = 0.0;
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
                        (diffusion_ * (state_values.shape_grad(i, q) *
                                       state_values.shape_grad(j, q)) +
                         reaction_ * phi_i * state_values.shape_value(j, q)) *
                        weight;
                      local_state_mass(i, j) +=
                        phi_i * state_values.shape_value(j, q) * weight;
                    }
                }
            }

          cell->get_dof_indices(state_indices);
          for (unsigned int face = 0;
               face < dealii::GeometryInfo<dim>::faces_per_cell;
               ++face)
            if (cell->face(face)->at_boundary() &&
                controlled_boundary_ids_.count(cell->face(face)->boundary_id()) != 0)
              {
                face_values.reinit(cell, face);
                local_boundary_mass = 0.0;
                for (unsigned int q = 0; q < face_quadrature.size(); ++q)
                  for (unsigned int i = 0; i < state_fe_.dofs_per_cell; ++i)
                    for (unsigned int j = 0; j < state_fe_.dofs_per_cell; ++j)
                      local_boundary_mass(i, j) +=
                        face_values.shape_value(i, q) *
                        face_values.shape_value(j, q) * face_values.JxW(q);
                for (unsigned int i = 0; i < state_fe_.dofs_per_cell; ++i)
                  {
                    const auto control_i =
                      control_index_for_state_dof_[state_indices[i]];
                    if (control_i == invalid_control_index)
                      continue;
                    for (unsigned int j = 0; j < state_fe_.dofs_per_cell; ++j)
                      {
                        const auto control_j =
                          control_index_for_state_dof_[state_indices[j]];
                        if (control_j != invalid_control_index &&
                            local_boundary_mass(i, j) != 0.0)
                          control_boundary_mass_->add(control_i,
                                                      control_j,
                                                      local_boundary_mass(i, j));
                      }
                  }
              }

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
            }
        }
    }

    void
    assemble_reduced_solve_operators()
    {
      build_reduced_state_matrix();
      build_reduced_lifting_coupling();
      Vector fixed_lifting_contribution(state_dof_handler_.n_dofs());
      physical_system_matrix_.vmult(fixed_lifting_contribution, fixed_lifting_);
      fixed_lifting_contribution *= -1.0;
      fixed_lifting_contribution += forcing_load_;
      reduced_forcing_load_ = pullback_state(fixed_lifting_contribution);
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
          Vector physical_column = embed_state(basis);
          Vector physical_result(state_dof_handler_.n_dofs());
          physical_system_matrix_.vmult(physical_result, physical_column);
          const Vector reduced_result = pullback_state(physical_result);
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
          Vector physical_column = embed_state(basis);
          Vector physical_result(state_dof_handler_.n_dofs());
          physical_system_matrix_.vmult(physical_result, physical_column);
          const Vector reduced_result = pullback_state(physical_result);
          for (std::size_t row = 0; row < size; ++row)
            if (reduced_result[row] != 0.0)
              reduced_system_matrix_.set(row, column, reduced_result[row]);
        }
    }

    void
    build_reduced_lifting_coupling()
    {
      const auto state_size = independent_state_dofs_.size();
      const auto control_size = controlled_state_dofs_.size();
      dealii::DynamicSparsityPattern dsp(state_size, control_size);
      for (std::size_t column = 0; column < control_size; ++column)
        {
          Vector basis(control_size);
          basis[column] = 1.0;
          Vector physical_result(state_dof_handler_.n_dofs());
          physical_system_matrix_.vmult(physical_result, lift_control(basis));
          const Vector reduced_result = pullback_state(physical_result);
          for (std::size_t row = 0; row < state_size; ++row)
            if (reduced_result[row] != 0.0)
              dsp.add(row, column);
        }
      reduced_lifting_sparsity_.copy_from(dsp);
      reduced_lifting_coupling_.reinit(reduced_lifting_sparsity_);
      for (std::size_t column = 0; column < control_size; ++column)
        {
          Vector basis(control_size);
          basis[column] = 1.0;
          Vector physical_result(state_dof_handler_.n_dofs());
          physical_system_matrix_.vmult(physical_result, lift_control(basis));
          const Vector reduced_result = pullback_state(physical_result);
          for (std::size_t row = 0; row < state_size; ++row)
            if (reduced_result[row] != 0.0)
              reduced_lifting_coupling_.set(row, column, reduced_result[row]);
        }
    }

    Vector
    reconstruct(const Vector &independent_state, const Vector &control) const
    {
      Vector physical = embed_state(independent_state);
      physical.add(1.0, lift_control(control));
      physical.add(1.0, fixed_lifting_);
      return physical;
    }

    Vector
    physical_tangent(const Vector &state_tangent,
                     const Vector &control_tangent) const
    {
      Vector physical = embed_state(state_tangent);
      physical.add(1.0, lift_control(control_tangent));
      return physical;
    }

    Vector
    embed_state(const Vector &independent_state) const
    {
      contract::require(independent_state.size() == independent_state_dofs_.size(),
                        "State coordinates have an incompatible layout");
      Vector physical(state_dof_handler_.n_dofs());
      reconstruction_.vmult(physical, independent_state);
      return physical;
    }

    Vector
    lift_control(const Vector &control) const
    {
      contract::require(control.size() == controlled_state_dofs_.size(),
                        "Dirichlet control has an incompatible trace layout");
      Vector physical(state_dof_handler_.n_dofs());
      for (std::size_t index = 0; index < controlled_state_dofs_.size(); ++index)
        physical[controlled_state_dofs_[index]] = control[index];
      return physical;
    }

    Vector
    pullback_state(const Vector &physical_covector) const
    {
      contract::require(physical_covector.size() == state_dof_handler_.n_dofs(),
                        "State pullback received an incompatible physical covector");
      Vector independent(independent_state_dofs_.size());
      reconstruction_.Tvmult(independent, physical_covector);
      return independent;
    }

    Vector
    pullback_control(const Vector &physical_covector) const
    {
      contract::require(physical_covector.size() == state_dof_handler_.n_dofs(),
                        "Control pullback received an incompatible physical covector");
      Vector control(controlled_state_dofs_.size());
      for (std::size_t index = 0; index < controlled_state_dofs_.size(); ++index)
        control[index] = physical_covector[controlled_state_dofs_[index]];
      return control;
    }

    static contract::LinearSolveReport
    solve_symmetric_system(const dealii::SparseMatrix<double> &matrix,
                           Vector &                             solution,
                           const Vector &                       right_hand_side,
                           const dealii_backend::SPDLinearSolvePolicy &policy)
    {
      return dealii_backend::solve_serial_spd(matrix,
                                              solution,
                                              right_hand_side,
                                              policy);
    }

    dealii::FE_Q<dim> state_fe_;
    dealii::DoFHandler<dim> state_dof_handler_;
    dealii::AffineConstraints<double> homogeneous_constraints_;
    std::vector<dealii::types::global_dof_index> independent_state_dofs_;
    std::vector<dealii::types::global_dof_index> controlled_state_dofs_;
    std::vector<std::size_t> control_index_for_state_dof_;

    const double diffusion_;
    const double reaction_;
    const double regularisation_weight_;
    const std::set<dealii::types::boundary_id> controlled_boundary_ids_;
    const std::set<dealii::types::boundary_id> fixed_boundary_ids_;
    std::vector<bool> fixed_state_dofs_;

    dealii::SparsityPattern reconstruction_sparsity_;
    dealii::SparseMatrix<double> reconstruction_;
    Vector fixed_lifting_;
    dealii::SparsityPattern physical_sparsity_;
    dealii::SparseMatrix<double> physical_system_matrix_;
    dealii::SparseMatrix<double> physical_state_mass_;
    dealii::SparsityPattern control_boundary_mass_sparsity_;
    std::shared_ptr<dealii::SparseMatrix<double>> control_boundary_mass_;
    Vector forcing_load_;
    Vector desired_state_load_;
    double desired_state_norm_ = 0.0;

    dealii::SparsityPattern reduced_state_sparsity_;
    dealii::SparseMatrix<double> reduced_system_matrix_;
    dealii::SparsityPattern reduced_lifting_sparsity_;
    dealii::SparseMatrix<double> reduced_lifting_coupling_;
    Vector reduced_forcing_load_;

    contract::LayoutPtr variable_layout_;
    contract::LayoutPtr test_layout_;
    contract::LayoutPtr state_layout_;
    contract::LayoutPtr control_layout_;
  };
} // namespace nmopt::compiler::v1::detail
