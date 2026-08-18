#pragma once

#include "nmopt/contract/executable_model.hpp"
#include "nmopt/contract/reduced_hessian.hpp"
#include "nmopt/contract/supplied_otd.hpp"
#include "nmopt/dealii/cellwise_box_constraint.hpp"
#include "nmopt/dealii/mass_metric.hpp"
#include "nmopt/dealii/serial_backend.hpp"
#include "nmopt/dealii/serial_spd_solver.hpp"
#include "nmopt/semantic/v1/types.hpp"

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
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <locale>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <stdexcept>
#include <utility>
#include <vector>

namespace nmopt::dealii_backend
{
  template <int dim>
  class ScalarDiffusionReactionKKT;

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
    , public contract::ReducedHessianT<SerialBackend>
  {
  public:
    using Vector = dealii::Vector<double>;
    using Primal = contract::PrimalBlockT<SerialBackend>;
    using Covector = contract::CovectorBlockT<SerialBackend>;
    using SolveResult = contract::FormulationSolveResultT<SerialBackend>;
    using SuppliedSystem = contract::SuppliedOTDSystemT<SerialBackend>;

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

    const contract::LayoutPtr &
    layout() const override
    {
      return control_layout_;
    }

    MassMetric
    control_l2_metric(
      MassMetricSolveParameters solve_parameters = {}) const
    {
      return MassMetric("l2_cellwise",
                        control_layout_,
                        control_mass_,
                        solve_parameters);
    }

    // Build the supplied-OTD product from the explicitly assembled weak
    // blocks owned by this scalar lowerer.  The returned callbacks do not
    // call the DTO residual or objective derivative methods.
    static SuppliedSystem
    make_supplied_otd_system(
      std::shared_ptr<const ScalarDiffusionReactionModel> model,
      std::shared_ptr<const void>                         lifetime_owner = {})
    {
      semantic::v1::SuppliedOTDDeclaration legacy_declaration;
      legacy_declaration.state_block.variable_space_id = "state";
      legacy_declaration.state_block.residual_space_id = "state_equation";
      legacy_declaration.state_block.runtime_variable_space_id = "state";
      legacy_declaration.state_block.runtime_residual_space_id =
        "state_equation";
      legacy_declaration.adjoint_block.variable_space_id = "state_test";
      legacy_declaration.adjoint_block.residual_space_id = "adjoint_equation";
      legacy_declaration.adjoint_block.runtime_variable_space_id = "state_test";
      legacy_declaration.adjoint_block.runtime_residual_space_id =
        "adjoint_equation";
      legacy_declaration.control_stationarity_block.variable_space_id =
        "control";
      legacy_declaration.control_stationarity_block.residual_space_id =
        "control_stationarity";
      legacy_declaration.control_stationarity_block.runtime_variable_space_id =
        "control";
      legacy_declaration.control_stationarity_block.runtime_residual_space_id =
        "control_stationarity";
      return make_supplied_otd_system(model,
                                      legacy_declaration,
                                      std::move(lifetime_owner));
    }

    static SuppliedSystem
    make_supplied_otd_system(
      std::shared_ptr<const ScalarDiffusionReactionModel> model,
      const semantic::v1::SuppliedOTDDeclaration &         declaration,
      std::shared_ptr<const void>                         lifetime_owner = {})
    {
      contract::require(static_cast<bool>(model),
                        "Supplied OTD lowerer needs a scalar model");
      const std::size_t state_dimension =
        static_cast<std::size_t>(model->system_matrix_.m());
      const std::size_t control_dimension =
        static_cast<std::size_t>(model->control_mass_->m());
      const auto variable_layout = std::make_shared<const contract::BlockLayout>(
        "supplied_otd_variables",
        std::vector<contract::SpaceId>{{declaration.state_block
                                          .runtime_variable_space_id},
                                       {declaration.adjoint_block
                                          .runtime_variable_space_id},
                                       {declaration.control_stationarity_block
                                          .runtime_variable_space_id}},
        std::vector<std::size_t>{state_dimension,
                                 state_dimension,
                                 control_dimension});
      const auto residual_layout = std::make_shared<const contract::BlockLayout>(
        "supplied_otd_residuals",
        std::vector<contract::SpaceId>{{declaration.state_block
                                          .runtime_residual_space_id},
                                       {declaration.adjoint_block
                                          .runtime_residual_space_id},
                                       {declaration.control_stationarity_block
                                          .runtime_residual_space_id}},
        std::vector<std::size_t>{state_dimension,
                                 state_dimension,
                                 control_dimension});
      const contract::SuppliedOTDLayout layout(variable_layout,
                                               residual_layout);
      const auto quadratic_kkt_validity =
        contract::make_canonical_supplied_otd_quadratic_kkt_validity();

      const auto residual = [model, residual_layout](const Primal &point) {
        Vector state(model->system_matrix_.m());
        model->system_matrix_.vmult(state, point.block(0));
        state.add(-1.0, model->forcing_load_);
        Vector control_contribution(model->system_matrix_.m());
        model->control_coupling_.vmult(control_contribution, point.block(2));
        state.add(-1.0, control_contribution);

        Vector adjoint(model->system_matrix_.m());
        model->system_matrix_.Tvmult(adjoint, point.block(1));
        Vector state_objective(model->system_matrix_.m());
        model->state_mass_.vmult(state_objective, point.block(0));
        adjoint.add(-1.0, state_objective);
        adjoint.add(1.0, model->desired_state_load_);

        Vector stationarity(model->control_mass_->m());
        model->control_coupling_.Tvmult(stationarity, point.block(1));
        Vector regularisation(model->control_mass_->m());
        model->control_mass_->vmult(regularisation, point.block(2));
        stationarity.add(model->regularisation_weight_, regularisation);

        return Covector(residual_layout,
                        {std::move(state),
                         std::move(adjoint),
                         std::move(stationarity)});
      };

      const auto residual_jvp = [model, residual_layout](const Primal &,
                                                          const Primal &tangent) {
        Vector state(model->system_matrix_.m());
        model->system_matrix_.vmult(state, tangent.block(0));
        Vector control_contribution(model->system_matrix_.m());
        model->control_coupling_.vmult(control_contribution, tangent.block(2));
        state.add(-1.0, control_contribution);

        Vector adjoint(model->system_matrix_.m());
        model->system_matrix_.Tvmult(adjoint, tangent.block(1));
        Vector state_objective(model->system_matrix_.m());
        model->state_mass_.vmult(state_objective, tangent.block(0));
        adjoint.add(-1.0, state_objective);

        Vector stationarity(model->control_mass_->m());
        model->control_coupling_.Tvmult(stationarity, tangent.block(1));
        Vector regularisation(model->control_mass_->m());
        model->control_mass_->vmult(regularisation, tangent.block(2));
        stationarity.add(model->regularisation_weight_, regularisation);

        return Covector(residual_layout,
                        {std::move(state),
                         std::move(adjoint),
                         std::move(stationarity)});
      };

      const auto residual_vjp = [model, variable_layout](const Primal &,
                                                          const Primal &seed) {
        Vector state(model->system_matrix_.m());
        model->system_matrix_.Tvmult(state, seed.block(0));
        Vector state_objective(model->system_matrix_.m());
        model->state_mass_.Tvmult(state_objective, seed.block(1));
        state.add(-1.0, state_objective);

        Vector adjoint(model->system_matrix_.m());
        model->system_matrix_.vmult(adjoint, seed.block(1));
        Vector stationarity_adjoint(model->system_matrix_.m());
        model->control_coupling_.vmult(stationarity_adjoint, seed.block(2));
        adjoint.add(1.0, stationarity_adjoint);

        Vector control(model->control_mass_->m());
        model->control_coupling_.Tvmult(control, seed.block(0));
        control *= -1.0;
        Vector regularisation(model->control_mass_->m());
        model->control_mass_->vmult(regularisation, seed.block(2));
        control.add(model->regularisation_weight_, regularisation);

        return Covector(variable_layout,
                        {std::move(state),
                         std::move(adjoint),
                         std::move(control)});
      };

      const auto solve = [model, variable_layout](const Primal &) {
        const auto state_dimension = model->system_matrix_.m();
        const auto control_dimension = model->control_mass_->m();
        const auto total_dimension = 2 * state_dimension + control_dimension;
        dealii::DynamicSparsityPattern dynamic_sparsity(total_dimension,
                                                         total_dimension);
        add_sparsity_entries(dynamic_sparsity,
                             model->system_matrix_,
                             0,
                             0,
                             false);
        add_sparsity_entries(dynamic_sparsity,
                             model->control_coupling_,
                             0,
                             2 * state_dimension,
                             false);
        add_sparsity_entries(dynamic_sparsity,
                             model->state_mass_,
                             state_dimension,
                             0,
                             false);
        add_sparsity_entries(dynamic_sparsity,
                             model->system_matrix_,
                             state_dimension,
                             state_dimension,
                             true);
        add_sparsity_entries(dynamic_sparsity,
                             model->control_coupling_,
                             2 * state_dimension,
                             state_dimension,
                             true);
        add_sparsity_entries(dynamic_sparsity,
                             *model->control_mass_,
                             2 * state_dimension,
                             2 * state_dimension,
                             false);

        dealii::SparsityPattern sparsity;
        sparsity.copy_from(dynamic_sparsity);
        dealii::SparseMatrix<double> optimality_matrix(sparsity);
        add_matrix_entries(optimality_matrix,
                           model->system_matrix_,
                           0,
                           0,
                           1.0,
                           false);
        add_matrix_entries(optimality_matrix,
                           model->control_coupling_,
                           0,
                           2 * state_dimension,
                           -1.0,
                           false);
        add_matrix_entries(optimality_matrix,
                           model->state_mass_,
                           state_dimension,
                           0,
                           -1.0,
                           false);
        add_matrix_entries(optimality_matrix,
                           model->system_matrix_,
                           state_dimension,
                           state_dimension,
                           1.0,
                           true);
        add_matrix_entries(optimality_matrix,
                           model->control_coupling_,
                           2 * state_dimension,
                           state_dimension,
                           1.0,
                           true);
        add_matrix_entries(optimality_matrix,
                           *model->control_mass_,
                           2 * state_dimension,
                           2 * state_dimension,
                           model->regularisation_weight_,
                           false);

        Vector right_hand_side(total_dimension);
        for (dealii::types::global_dof_index index = 0;
             index < state_dimension;
             ++index)
          {
            right_hand_side[index] = model->forcing_load_[index];
            right_hand_side[state_dimension + index] =
              -model->desired_state_load_[index];
          }
        Vector solution(total_dimension);
        dealii::SparseDirectUMFPACK solver;
        solver.initialize(optimality_matrix);
        solver.vmult(solution, right_hand_side);

        Vector state(state_dimension);
        Vector adjoint(state_dimension);
        Vector control(control_dimension);
        for (dealii::types::global_dof_index index = 0;
             index < state_dimension;
             ++index)
          {
            state[index] = solution[index];
            adjoint[index] = solution[state_dimension + index];
          }
        for (dealii::types::global_dof_index index = 0;
             index < control_dimension;
             ++index)
          control[index] = solution[2 * state_dimension + index];

        return SolveResult(
          Primal(variable_layout,
                 {std::move(state), std::move(adjoint), std::move(control)}),
          contract::LinearSolveReport{"serial_sparse_direct_umfpack",
                                       "not applicable",
                                       1,
                                       1,
                                       0.0,
                                       0.0,
                                       0.0,
                                       0.0,
                                       contract::LinearSolveTermination::converged});
      };

      return SuppliedSystem(
        layout,
        residual,
        residual_jvp,
        residual_vjp,
        solve,
        quadratic_kkt_validity,
        std::move(lifetime_owner));
    }

    CellwiseBoxConstraint
    control_l2_box_constraint(Vector            lower,
                              Vector            upper,
                              const MassMetric &projection_metric) const
    {
      return CellwiseBoxConstraint(control_layout_,
                                   std::move(lower),
                                   std::move(upper),
                                   projection_metric);
    }

    CellwiseBoxConstraint
    control_l2_box_constraint(const double      lower,
                              const double      upper,
                              const MassMetric &projection_metric) const
    {
      return CellwiseBoxConstraint(control_layout_,
                                   lower,
                                   upper,
                                   projection_metric);
    }

    const dealii::AffineConstraints<double> &
    state_constraints() const
    {
      return state_constraints_;
    }

    void
    write_field_output(const std::filesystem::path &directory,
                       const Primal &                 state,
                       const Primal &                 control,
                       const Primal &                 adjoint) const
    {
      contract::require(state.layout()->compatible_with(*state_layout_),
                        "Field output state has an incompatible layout");
      contract::require(control.layout()->compatible_with(*control_layout_),
                        "Field output control has an incompatible layout");
      contract::require(adjoint.layout()->compatible_with(*test_layout_),
                        "Field output adjoint has an incompatible layout");

      std::filesystem::create_directories(directory);
      dealii::DataOut<dim> data_out;
      data_out.attach_dof_handler(state_dof_handler_);
      data_out.add_data_vector(state.block(0), "state");
      data_out.add_data_vector(adjoint.block(0), "adjoint");
      data_out.add_data_vector(control_dof_handler_,
                               control.block(0),
                               "control");
      data_out.build_patches();

      std::ofstream output(directory / "fields-volume.vtu");
      if (!output)
        throw std::runtime_error("could not open scalar volume field output");
      output.imbue(std::locale::classic());
      data_out.write_vtu(output);
      if (!output)
        throw std::runtime_error("could not write scalar volume field output");
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

      Vector state(state_dof_handler_.n_dofs());
      state_mass_.vmult(state, variables.block(0));
      state.add(-1.0, desired_state_load_);

      Vector control(control_dof_handler_.n_dofs());
      control_mass_->vmult(control, variables.block(1));
      control *= regularisation_weight_;
      return Covector(variable_layout_, {std::move(state), std::move(control)});
    }

    Covector
    apply(const Primal &control, const Primal &direction) const override
    {
      contract::require(control.layout()->compatible_with(*control_layout_),
                        "Reduced Hessian control has an incompatible layout");
      contract::require(direction.layout()->compatible_with(*control_layout_),
                        "Reduced Hessian direction has an incompatible layout");

      Vector tangent_rhs(state_dof_handler_.n_dofs());
      control_coupling_.vmult(tangent_rhs, direction.block(0));
      Vector tangent_state(state_dof_handler_.n_dofs());
      const auto tangent_report =
        solve_symmetric_system(tangent_state, tangent_rhs, {});
      contract::require(tangent_report.converged(),
                        "Reduced Hessian tangent solve did not converge");
      state_constraints_.distribute(tangent_state);

      Vector incremental_adjoint_rhs(state_dof_handler_.n_dofs());
      state_mass_.vmult(incremental_adjoint_rhs, tangent_state);
      Vector incremental_adjoint(state_dof_handler_.n_dofs());
      const auto incremental_adjoint_report = solve_symmetric_system(
        incremental_adjoint, incremental_adjoint_rhs, {});
      contract::require(incremental_adjoint_report.converged(),
                        "Reduced Hessian incremental-adjoint solve did not converge");
      state_constraints_.distribute(incremental_adjoint);

      Vector action(control_dof_handler_.n_dofs());
      control_coupling_.Tvmult(action, incremental_adjoint);
      Vector regularisation_action(control_dof_handler_.n_dofs());
      control_mass_->vmult(regularisation_action, direction.block(0));
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
    solve_state_with_report(const Primal &             control,
                            const SPDLinearSolvePolicy &policy) const
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
      const Primal &              full_point,
      const Covector &            state_objective_derivative,
      const SPDLinearSolvePolicy &policy) const
    {
      require_variables(full_point, "Adjoint solve point");
      contract::require(
        state_objective_derivative.layout()->compatible_with(*state_layout_),
        "Adjoint solve right-hand side has an incompatible state layout");

      // The v0 diffusion-reaction operator is symmetric. The VJP nevertheless
      // uses Tvmult above, so a non-symmetric extension cannot silently reuse
      // this solve path.
      Vector adjoint(test_layout_->dimension(0));
      auto report = solve_symmetric_system(adjoint,
                                           state_objective_derivative.block(0),
                                           policy);
      state_constraints_.distribute(adjoint);
      return {Primal(test_layout_, {std::move(adjoint)}), std::move(report)};
    }

  private:
    friend class ScalarDiffusionReactionKKT<dim>;

    static void
    add_sparsity_entries(
      dealii::DynamicSparsityPattern &       target,
      const dealii::SparseMatrix<double> &   source,
      const dealii::types::global_dof_index row_offset,
      const dealii::types::global_dof_index column_offset,
      const bool                             transpose)
    {
      for (dealii::types::global_dof_index row = 0; row < source.m(); ++row)
        for (auto entry = source.begin(row); entry != source.end(row); ++entry)
          target.add(row_offset + (transpose ? entry->column() : row),
                     column_offset + (transpose ? row : entry->column()));
    }

    static void
    add_matrix_entries(
      dealii::SparseMatrix<double> &         target,
      const dealii::SparseMatrix<double> &   source,
      const dealii::types::global_dof_index row_offset,
      const dealii::types::global_dof_index column_offset,
      const double                           factor,
      const bool                             transpose)
    {
      for (dealii::types::global_dof_index row = 0; row < source.m(); ++row)
        for (auto entry = source.begin(row); entry != source.end(row); ++entry)
          target.add(row_offset + (transpose ? entry->column() : row),
                     column_offset + (transpose ? row : entry->column()),
                     factor * entry->value());
    }

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
      control_mass_ = std::make_shared<dealii::SparseMatrix<double>>();
      control_mass_->reinit(control_mass_sparsity_);

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
              control_mass_->add(control_indices[i],
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

    contract::LinearSolveReport
    solve_symmetric_system(Vector &                     solution,
                           const Vector &               right_hand_side,
                           const SPDLinearSolvePolicy &policy) const
    {
      return solve_serial_spd(system_matrix_, solution, right_hand_side, policy);
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
    std::shared_ptr<dealii::SparseMatrix<double>> control_mass_;
    Vector forcing_load_;
    Vector desired_state_load_;
    double desired_state_norm_ = 0.0;

    contract::LayoutPtr variable_layout_;
    contract::LayoutPtr test_layout_;
    contract::LayoutPtr state_layout_;
    contract::LayoutPtr control_layout_;
  };
} // namespace nmopt::dealii_backend
