#pragma once

#include "nmopt/compiler/v1/dealii_reference_cell.hpp"
#include "nmopt/contract/executable_model.hpp"
#include "nmopt/contract/reduced_hessian.hpp"
#include "nmopt/dealii/hminus1_metric.hpp"
#include "nmopt/dealii/mass_metric.hpp"
#include "nmopt/dealii/serial_backend.hpp"
#include "nmopt/dealii/serial_spd_solver.hpp"

#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <limits>
#include <locale>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nmopt::compiler::v1::detail
{
  // The continuous-control target used by P2.3 and P5.2. It owns one
  // conforming scalar Lagrange control realization while keeping state
  // observation, control loss, and search metric as independent compiler
  // selections. Hypercube meshes use FE_Q and simplex meshes use
  // FE_SimplexP:
  //
  //   r(y,u) = A y - f_h - B u,
  //   J(y,u) = J_tracking(y) + alpha/2 u^T R_u u.
  template <int dim>
  class ContinuousControlModel final
    : public contract::ExecutableModelT<dealii_backend::SerialBackend>
    , public contract::ReducedHessianT<dealii_backend::SerialBackend>
  {
  public:
    using Backend = dealii_backend::SerialBackend;
    using Vector = dealii::Vector<double>;
    using Primal = contract::PrimalBlockT<Backend>;
    using Covector = contract::CovectorBlockT<Backend>;
    using SolveResult = contract::FormulationSolveResultT<Backend>;

    ContinuousControlModel(
      dealii::Triangulation<dim> &                triangulation,
      const dealii::Function<dim> &               forcing,
      const dealii::Function<dim> &               desired_state,
      const double                                  diffusion,
      const double                                  reaction,
      const double                                  regularisation_weight,
      const unsigned int                            state_degree,
      std::set<dealii::types::boundary_id>         dirichlet_boundary_ids,
      const bool                                    use_h1_state_observation = false,
      const bool                                    use_h1_control_regularisation = true,
      const bool homogeneous_dirichlet_control = false,
      std::set<dealii::types::boundary_id>         control_boundary_ids = {})
      : state_fe_(make_scalar_lagrange_element(
          triangulation, state_degree, "The continuous-control target"))
      , control_fe_(make_scalar_lagrange_element(
          triangulation, state_degree, "The continuous-control target"))
      , state_dof_handler_(triangulation)
      , control_dof_handler_(triangulation)
      , diffusion_(diffusion)
      , reaction_(reaction)
      , regularisation_weight_(regularisation_weight)
      , dirichlet_boundary_ids_(std::move(dirichlet_boundary_ids))
      , control_boundary_ids_(control_boundary_ids.empty()
                               ? dirichlet_boundary_ids_
                               : std::move(control_boundary_ids))
      , use_h1_state_observation_(use_h1_state_observation)
      , use_h1_control_regularisation_(use_h1_control_regularisation)
      , homogeneous_dirichlet_control_(homogeneous_dirichlet_control)
    {
      contract::require(diffusion_ > 0.0,
                        "Diffusion coefficient must be strictly positive");
      contract::require(reaction_ >= 0.0,
                        "Reaction coefficient must be non-negative");
      contract::require(regularisation_weight_ > 0.0,
                        "Control regularisation weight must be strictly positive");
      contract::require(state_degree > 0,
                        "State and continuous-control FE degrees must be positive");
      contract::require(!dirichlet_boundary_ids_.empty(),
                        "The continuous-control v1 target needs a fixed Dirichlet boundary");
      contract::require(!homogeneous_dirichlet_control ||
                          !control_boundary_ids_.empty(),
                        "The homogeneous continuous-control target needs a fixed control boundary");

      state_dof_handler_.distribute_dofs(*state_fe_);
      control_dof_handler_.distribute_dofs(*control_fe_);
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

    std::size_t
    physical_state_dimension() const
    {
      return state_dof_handler_.n_dofs();
    }

    std::size_t
    independent_state_dimension() const
    {
      return physical_state_dimension() -
             static_cast<std::size_t>(std::count(constrained_state_dofs_.begin(),
                                                 constrained_state_dofs_.end(),
                                                 true));
    }

    std::size_t
    physical_control_dimension() const
    {
      return control_dof_handler_.n_dofs();
    }

    std::size_t
    independent_control_dimension() const
    {
      return independent_control_dofs_.size();
    }

    const contract::LayoutPtr &
    layout() const override
    {
      return control_layout_;
    }

    dealii_backend::MassMetric
    control_l2_metric(
      dealii_backend::MassMetricSolveParameters solve_parameters = {}) const
    {
      return dealii_backend::MassMetric("l2_continuous",
                                        control_layout_,
                                        control_mass_,
                                        solve_parameters);
    }

    dealii_backend::MassMetric
    control_h1_metric(
      dealii_backend::MassMetricSolveParameters solve_parameters = {}) const
    {
      return dealii_backend::MassMetric("h1_continuous",
                                        control_layout_,
                                        control_h1_matrix_,
                                        solve_parameters);
    }

    dealii_backend::Hminus1Metric
    control_hminus1_metric(
      dealii_backend::MetricSolveParameters solve_parameters = {},
      const dealii_backend::Hminus1OperatorRealisation operator_realisation =
        dealii_backend::Hminus1OperatorRealisation::mass_laplacian_inverse_mass,
      const dealii_backend::Hminus1InverseRealisation inverse_realisation =
        dealii_backend::Hminus1InverseRealisation::mass_inverse_laplacian_mass_inverse) const
    {
      contract::require(
        homogeneous_dirichlet_control_,
        "The H-1 metric requires independent homogeneous-Dirichlet control coordinates");
      return dealii_backend::Hminus1Metric("hminus1_continuous",
                                           control_layout_,
                                           control_mass_,
                                           control_stiffness_,
                                           solve_parameters,
                                           operator_realisation,
                                           inverse_realisation);
    }

    void
    write_native_output(const std::filesystem::path &directory,
                        const Primal &                 state,
                        const Primal &                 control,
                        const Primal &                 adjoint,
                        const dealii::Function<dim> *  forcing = nullptr,
                        const dealii::Function<dim> *  desired_state = nullptr) const
    {
      contract::require(state.layout()->compatible_with(*state_layout_),
                        "Continuous-control output state has an incompatible layout");
      contract::require(control.layout()->compatible_with(*control_layout_),
                        "Continuous-control output control has an incompatible layout");
      contract::require(adjoint.layout()->compatible_with(*test_layout_),
                        "Continuous-control output adjoint has an incompatible layout");
      contract::require((forcing == nullptr) == (desired_state == nullptr),
                        "Continuous-control output needs both forcing and target functions");

      std::filesystem::create_directories(directory);

      dealii::DataOut<dim> mesh_out;
      mesh_out.attach_triangulation(state_dof_handler_.get_triangulation());
      mesh_out.build_patches();
      std::ofstream mesh_output(directory / "mesh-volume.vtu");
      if (!mesh_output)
        throw std::runtime_error(
          "could not open continuous-control volume mesh output");
      mesh_output.imbue(std::locale::classic());
      mesh_out.write_vtu(mesh_output);
      if (!mesh_output)
        throw std::runtime_error(
          "could not write continuous-control volume mesh output");

      if constexpr (dim == 2)
        {
          dealii::GridOut grid_out;
          std::ofstream  svg_output(directory / "mesh-volume.svg");
          if (!svg_output)
            throw std::runtime_error(
              "could not open continuous-control volume mesh SVG output");
          svg_output.imbue(std::locale::classic());
          grid_out.write_svg(state_dof_handler_.get_triangulation(),
                             svg_output);
          if (!svg_output)
            throw std::runtime_error(
              "could not write continuous-control volume mesh SVG output");
        }

      Vector full_control(control_dof_handler_.n_dofs());
      for (std::size_t index = 0; index < independent_control_dofs_.size();
           ++index)
        full_control[independent_control_dofs_[index]] = control.block(0)[index];
      control_constraints_.distribute(full_control);

      dealii::DataOut<dim> data_out;
      data_out.attach_dof_handler(state_dof_handler_);
      data_out.add_data_vector(state.block(0), "state");
      data_out.add_data_vector(adjoint.block(0), "adjoint");
      Vector negative_adjoint = adjoint.block(0);
      negative_adjoint *= -1.0;
      data_out.add_data_vector(negative_adjoint, "negative_adjoint");
      data_out.add_data_vector(control_dof_handler_, full_control, "control");
      if (forcing != nullptr)
        {
          Vector forcing_values(state_dof_handler_.n_dofs());
          Vector desired_state_values(state_dof_handler_.n_dofs());
          dealii::VectorTools::interpolate(state_dof_handler_,
                                            *forcing,
                                            forcing_values);
          dealii::VectorTools::interpolate(state_dof_handler_,
                                            *desired_state,
                                            desired_state_values);
          data_out.add_data_vector(forcing_values, "forcing");
          data_out.add_data_vector(desired_state_values, "target");
        }
      data_out.build_patches();

      std::ofstream output(directory / "fields-volume.vtu");
      if (!output)
        throw std::runtime_error(
          "could not open continuous-control volume field output");
      output.imbue(std::locale::classic());
      data_out.write_vtu(output);
      if (!output)
        throw std::runtime_error(
          "could not write continuous-control volume field output");
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
    residual_vjp(const Primal &variables, const Primal &test_seed) const override
    {
      require_variables(variables, "Residual VJP");
      contract::require(test_seed.layout()->compatible_with(*test_layout_),
                        "Residual VJP seed has an incompatible test layout");
      Vector state(state_dof_handler_.n_dofs());
      system_matrix_.Tvmult(state, test_seed.block(0));
      Vector control(control_layout_->dimension(0));
      control_coupling_.Tvmult(control, test_seed.block(0));
      control *= -1.0;
      return Covector(variable_layout_, {std::move(state), std::move(control)});
    }

    double
    objective(const Primal &variables) const override
    {
      require_variables(variables, "Objective");
      Vector state_tracking_times_state(state_dof_handler_.n_dofs());
      state_tracking_matrix_.vmult(state_tracking_times_state,
                                   variables.block(0));
      const double state_value =
        0.5 * (variables.block(0) * state_tracking_times_state) -
        (desired_state_load_ * variables.block(0)) + 0.5 * desired_state_norm_;

      Vector control_regularisation_times_control(control_layout_->dimension(0));
      control_regularisation_matrix().vmult(
        control_regularisation_times_control, variables.block(1));
      const double control_value = 0.5 * regularisation_weight_ *
                                   (variables.block(1) *
                                    control_regularisation_times_control);
      return state_value + control_value;
    }

    Covector
    objective_derivative(const Primal &variables) const override
    {
      require_variables(variables, "Objective derivative");
      Vector state(state_dof_handler_.n_dofs());
      state_tracking_matrix_.vmult(state, variables.block(0));
      state.add(-1.0, desired_state_load_);

      Vector control(control_layout_->dimension(0));
      control_regularisation_matrix().vmult(control, variables.block(1));
      control *= regularisation_weight_;
      return Covector(variable_layout_, {std::move(state), std::move(control)});
    }

    Covector
    apply(const Primal &control, const Primal &direction) const override
    {
      contract::require(control.layout()->compatible_with(*control_layout_),
                        "Continuous-control Hessian control has an incompatible layout");
      contract::require(direction.layout()->compatible_with(*control_layout_),
                        "Continuous-control Hessian direction has an incompatible layout");

      Vector tangent_rhs(state_dof_handler_.n_dofs());
      control_coupling_.vmult(tangent_rhs, direction.block(0));
      Vector tangent_state(state_dof_handler_.n_dofs());
      const auto tangent_report =
        solve_symmetric_system(tangent_state, tangent_rhs, {});
      contract::require(tangent_report.converged(),
                        "Continuous-control Hessian tangent solve did not converge");
      state_constraints_.distribute(tangent_state);

      Vector incremental_adjoint_rhs(state_dof_handler_.n_dofs());
      state_tracking_matrix_.vmult(incremental_adjoint_rhs, tangent_state);
      Vector incremental_adjoint(state_dof_handler_.n_dofs());
      const auto incremental_adjoint_report = solve_symmetric_system(
        incremental_adjoint, incremental_adjoint_rhs, {});
      contract::require(
        incremental_adjoint_report.converged(),
        "Continuous-control Hessian incremental-adjoint solve did not converge");
      state_constraints_.distribute(incremental_adjoint);

      Vector action(control_layout_->dimension(0));
      control_coupling_.Tvmult(action, incremental_adjoint);
      Vector regularisation_action(control_layout_->dimension(0));
      control_regularisation_matrix().vmult(regularisation_action,
                                             direction.block(0));
      action.add(regularisation_weight_, regularisation_action);
      return Covector(control_layout_, {std::move(action)});
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
      Vector right_hand_side = forcing_load_;
      Vector control_contribution(state_dof_handler_.n_dofs());
      control_coupling_.vmult(control_contribution, control.block(0));
      right_hand_side.add(1.0, control_contribution);
      Vector state(state_dof_handler_.n_dofs());
      auto report = solve_symmetric_system(state, right_hand_side, policy);
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
      Vector adjoint(test_layout_->dimension(0));
      auto report = solve_symmetric_system(adjoint,
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
    }

    const dealii::SparseMatrix<double> &
    control_regularisation_matrix() const
    {
      return use_h1_control_regularisation_ ? *control_h1_matrix_
                                            : *control_mass_;
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
            "The continuous-control v1 target does not support hanging, periodic, or "
            "other non-Dirichlet affine state constraints");

      control_constraints_.clear();
      dealii::DoFTools::make_hanging_node_constraints(control_dof_handler_,
                                                       control_constraints_);
      if (homogeneous_dirichlet_control_)
        for (const auto boundary_id : control_boundary_ids_)
          dealii::VectorTools::interpolate_boundary_values(control_dof_handler_,
                                                            boundary_id,
                                                            zero,
                                                            control_constraints_);
      control_constraints_.close();

      std::map<dealii::types::global_dof_index, double>
        control_dirichlet_values;
      if (homogeneous_dirichlet_control_)
        for (const auto boundary_id : control_boundary_ids_)
          dealii::VectorTools::interpolate_boundary_values(
            control_dof_handler_, boundary_id, zero, control_dirichlet_values);
      constrained_control_dofs_.assign(control_dof_handler_.n_dofs(), false);
      for (const auto &entry : control_dirichlet_values)
        constrained_control_dofs_.at(entry.first) = true;
      for (dealii::types::global_dof_index index = 0;
           index < control_dof_handler_.n_dofs();
           ++index)
        if (control_constraints_.is_constrained(index) &&
            !constrained_control_dofs_.at(index))
          throw contract::ContractError(
            "The continuous-control v1 target does not support hanging, periodic, or other non-Dirichlet affine control constraints");

      control_to_independent_.assign(
        control_dof_handler_.n_dofs(), std::numeric_limits<std::size_t>::max());
      for (dealii::types::global_dof_index index = 0;
           index < control_dof_handler_.n_dofs();
           ++index)
        if (!constrained_control_dofs_.at(index))
          {
            control_to_independent_.at(index) = independent_control_dofs_.size();
            independent_control_dofs_.push_back(index);
          }
      contract::require(!independent_control_dofs_.empty(),
                        "The continuous-control target needs an independent control DoF");
    }

    void
    initialise_storage()
    {
      const auto state_size = state_dof_handler_.n_dofs();
      const auto control_size = independent_control_dofs_.size();
      variable_layout_ = std::make_shared<const contract::BlockLayout>(
        "h1_control_variables",
        std::vector<contract::SpaceId>{{"state"}, {"control"}},
        std::vector<std::size_t>{state_size, control_size});
      test_layout_ = std::make_shared<const contract::BlockLayout>(
        "h1_control_state_test",
        std::vector<contract::SpaceId>{{"state_test"}},
        std::vector<std::size_t>{state_size});
      state_layout_ = variable_layout_->single_block(0, "state");
      control_layout_ = variable_layout_->single_block(1, "control");

      dealii::DynamicSparsityPattern state_dsp(state_size, state_size);
      dealii::DoFTools::make_sparsity_pattern(state_dof_handler_, state_dsp);
      state_sparsity_.copy_from(state_dsp);
      system_matrix_.reinit(state_sparsity_);
      state_tracking_matrix_.reinit(state_sparsity_);

      dealii::DynamicSparsityPattern control_dsp(state_size, control_size);
      std::vector<dealii::types::global_dof_index> state_indices(
        state_fe_->dofs_per_cell);
      std::vector<dealii::types::global_dof_index> control_indices(
        control_fe_->dofs_per_cell);
      auto state_cell = state_dof_handler_.begin_active();
      auto control_cell = control_dof_handler_.begin_active();
      for (; state_cell != state_dof_handler_.end(); ++state_cell, ++control_cell)
        {
          contract::require(control_cell != control_dof_handler_.end(),
                            "State and control DoF handlers do not share cells");
          state_cell->get_dof_indices(state_indices);
          control_cell->get_dof_indices(control_indices);
          for (const auto state_index : state_indices)
            if (!constrained_state_dofs_.at(state_index))
              for (const auto control_index : control_indices)
                if (!constrained_control_dofs_.at(control_index))
                  control_dsp.add(state_index,
                                  control_to_independent_.at(control_index));
        }
      contract::require(control_cell == control_dof_handler_.end(),
                        "State and control DoF handlers have different cells");
      control_sparsity_.copy_from(control_dsp);
      control_coupling_.reinit(control_sparsity_);

      dealii::DynamicSparsityPattern control_dsp_square(control_size, control_size);
      control_cell = control_dof_handler_.begin_active();
      for (; control_cell != control_dof_handler_.end(); ++control_cell)
        {
          control_cell->get_dof_indices(control_indices);
          for (const auto row : control_indices)
            if (!constrained_control_dofs_.at(row))
              for (const auto column : control_indices)
                if (!constrained_control_dofs_.at(column))
                  control_dsp_square.add(control_to_independent_.at(row),
                                         control_to_independent_.at(column));
        }
      control_sparsity_square_.copy_from(control_dsp_square);
      control_mass_ = std::make_shared<dealii::SparseMatrix<double>>();
      control_mass_->reinit(control_sparsity_square_);
      control_stiffness_ = std::make_shared<dealii::SparseMatrix<double>>();
      control_stiffness_->reinit(control_sparsity_square_);
      control_h1_matrix_ = std::make_shared<dealii::SparseMatrix<double>>();
      control_h1_matrix_->reinit(control_sparsity_square_);

      forcing_load_.reinit(state_size);
      desired_state_load_.reinit(state_size);
    }

    void
    assemble(const dealii::Function<dim> &forcing,
             const dealii::Function<dim> &desired_state)
    {
      const unsigned int quadrature_order =
        std::max(state_fe_->degree, control_fe_->degree) + 2;
      const auto quadrature = make_gauss_volume_quadrature(
        state_dof_handler_.get_triangulation(),
        quadrature_order,
        "The continuous-control target");
      dealii::FEValues<dim> state_values(
        *state_fe_,
        *quadrature,
        dealii::update_values | dealii::update_gradients |
          dealii::update_quadrature_points | dealii::update_JxW_values);
      dealii::FEValues<dim> control_values(
        *control_fe_,
        *quadrature,
        dealii::update_values | dealii::update_gradients);

      dealii::FullMatrix<double> local_system(state_fe_->dofs_per_cell,
                                              state_fe_->dofs_per_cell);
      dealii::FullMatrix<double> local_state_tracking(state_fe_->dofs_per_cell,
                                                      state_fe_->dofs_per_cell);
      dealii::FullMatrix<double> local_control_coupling(
        state_fe_->dofs_per_cell, control_fe_->dofs_per_cell);
      dealii::FullMatrix<double> local_control_mass(control_fe_->dofs_per_cell,
                                                    control_fe_->dofs_per_cell);
      dealii::FullMatrix<double> local_control_stiffness(
        control_fe_->dofs_per_cell, control_fe_->dofs_per_cell);
      dealii::Vector<double> local_forcing(state_fe_->dofs_per_cell);
      dealii::Vector<double> local_desired_state(state_fe_->dofs_per_cell);
      std::vector<dealii::types::global_dof_index> state_indices(
        state_fe_->dofs_per_cell);
      std::vector<dealii::types::global_dof_index> control_indices(
        control_fe_->dofs_per_cell);

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
          local_control_stiffness = 0.0;
          local_forcing = 0.0;
          local_desired_state = 0.0;

          for (unsigned int q = 0; q < quadrature->size(); ++q)
            {
              const double weight = state_values.JxW(q);
              const double forcing_value =
                forcing.value(state_values.quadrature_point(q));
              const double desired_value =
                desired_state.value(state_values.quadrature_point(q));
              const dealii::Tensor<1, dim> desired_gradient =
                use_h1_state_observation_
                  ? desired_state.gradient(state_values.quadrature_point(q))
                  : dealii::Tensor<1, dim>();
              desired_state_norm_ +=
                (desired_value * desired_value +
                 (use_h1_state_observation_ ? desired_gradient * desired_gradient
                                            : 0.0)) *
                weight;
              for (unsigned int i = 0; i < state_fe_->dofs_per_cell; ++i)
                {
                  const double phi_i = state_values.shape_value(i, q);
                  local_forcing(i) += forcing_value * phi_i * weight;
                  local_desired_state(i) +=
                    (desired_value * phi_i +
                     (use_h1_state_observation_
                        ? desired_gradient * state_values.shape_grad(i, q)
                        : 0.0)) *
                    weight;
                  for (unsigned int j = 0; j < state_fe_->dofs_per_cell; ++j)
                    {
                      local_system(i, j) +=
                        (diffusion_ * (state_values.shape_grad(i, q) *
                                       state_values.shape_grad(j, q)) +
                         reaction_ * phi_i * state_values.shape_value(j, q)) *
                        weight;
                      local_state_tracking(i, j) +=
                        (phi_i * state_values.shape_value(j, q) +
                         (use_h1_state_observation_
                            ? state_values.shape_grad(i, q) *
                                state_values.shape_grad(j, q)
                            : 0.0)) *
                        weight;
                    }
                  for (unsigned int j = 0; j < control_fe_->dofs_per_cell; ++j)
                    local_control_coupling(i, j) +=
                      phi_i * control_values.shape_value(j, q) * weight;
                }
              for (unsigned int i = 0; i < control_fe_->dofs_per_cell; ++i)
                for (unsigned int j = 0; j < control_fe_->dofs_per_cell; ++j)
                  {
                    local_control_mass(i, j) +=
                      control_values.shape_value(i, q) *
                      control_values.shape_value(j, q) * weight;
                    local_control_stiffness(i, j) +=
                      control_values.shape_grad(i, q) *
                      control_values.shape_grad(j, q) * weight;
                  }
            }

          state_cell->get_dof_indices(state_indices);
          control_cell->get_dof_indices(control_indices);
          for (unsigned int i = 0; i < state_fe_->dofs_per_cell; ++i)
            {
              const auto global_i = state_indices[i];
              if (constrained_state_dofs_.at(global_i))
                continue;
              forcing_load_[global_i] += local_forcing(i);
              desired_state_load_[global_i] += local_desired_state(i);
              for (unsigned int j = 0; j < state_fe_->dofs_per_cell; ++j)
                {
                  const auto global_j = state_indices[j];
                  if (!constrained_state_dofs_.at(global_j))
                    {
                      system_matrix_.add(global_i, global_j, local_system(i, j));
                      state_tracking_matrix_.add(global_i,
                                                 global_j,
                                                 local_state_tracking(i, j));
                    }
                }
              for (unsigned int j = 0; j < control_fe_->dofs_per_cell; ++j)
                if (!constrained_control_dofs_.at(control_indices[j]))
                  control_coupling_.add(
                    global_i,
                    control_to_independent_.at(control_indices[j]),
                    local_control_coupling(i, j));
            }
          for (unsigned int i = 0; i < control_fe_->dofs_per_cell; ++i)
            if (!constrained_control_dofs_.at(control_indices[i]))
              for (unsigned int j = 0; j < control_fe_->dofs_per_cell; ++j)
                if (!constrained_control_dofs_.at(control_indices[j]))
                  {
                    const auto row =
                      control_to_independent_.at(control_indices[i]);
                    const auto column =
                      control_to_independent_.at(control_indices[j]);
                    control_mass_->add(row, column, local_control_mass(i, j));
                    control_stiffness_->add(row,
                                            column,
                                            local_control_stiffness(i, j));
                    control_h1_matrix_->add(
                      row,
                      column,
                      local_control_mass(i, j) + local_control_stiffness(i, j));
                  }
        }
      contract::require(control_cell == control_dof_handler_.end(),
                        "State and control DoF handlers have different cells");
      for (dealii::types::global_dof_index index = 0;
           index < state_dof_handler_.n_dofs();
           ++index)
        if (constrained_state_dofs_.at(index))
          system_matrix_.set(index, index, 1.0);
    }

    contract::LinearSolveReport
    solve_symmetric_system(
      Vector &                                      solution,
      const Vector &                                right_hand_side,
      const dealii_backend::SPDLinearSolvePolicy &policy) const
    {
      return dealii_backend::solve_serial_spd(system_matrix_,
                                              solution,
                                              right_hand_side,
                                              policy);
    }

    std::unique_ptr<dealii::FiniteElement<dim>> state_fe_;
    std::unique_ptr<dealii::FiniteElement<dim>> control_fe_;
    dealii::DoFHandler<dim> state_dof_handler_;
    dealii::DoFHandler<dim> control_dof_handler_;
    dealii::AffineConstraints<double> state_constraints_;
    dealii::AffineConstraints<double> control_constraints_;
    std::vector<bool> constrained_state_dofs_;
    std::vector<bool> constrained_control_dofs_;
    std::vector<dealii::types::global_dof_index> independent_control_dofs_;
    std::vector<std::size_t> control_to_independent_;

    const double diffusion_;
    const double reaction_;
    const double regularisation_weight_;
    const std::set<dealii::types::boundary_id> dirichlet_boundary_ids_;
    const std::set<dealii::types::boundary_id> control_boundary_ids_;
    const bool use_h1_state_observation_;
    const bool use_h1_control_regularisation_;
    const bool homogeneous_dirichlet_control_;

    dealii::SparsityPattern state_sparsity_;
    dealii::SparsityPattern control_sparsity_;
    dealii::SparsityPattern control_sparsity_square_;
    dealii::SparseMatrix<double> system_matrix_;
    dealii::SparseMatrix<double> state_tracking_matrix_;
    dealii::SparseMatrix<double> control_coupling_;
    std::shared_ptr<dealii::SparseMatrix<double>> control_mass_;
    std::shared_ptr<dealii::SparseMatrix<double>> control_stiffness_;
    std::shared_ptr<dealii::SparseMatrix<double>> control_h1_matrix_;
    Vector forcing_load_;
    Vector desired_state_load_;
    double desired_state_norm_ = 0.0;

    contract::LayoutPtr variable_layout_;
    contract::LayoutPtr test_layout_;
    contract::LayoutPtr state_layout_;
    contract::LayoutPtr control_layout_;
  };
} // namespace nmopt::compiler::v1::detail
