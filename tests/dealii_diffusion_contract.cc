#include "nmopt/contract/reduced_dto.hpp"
#include "nmopt/compiler/v1/dealii_scalar_diffusion_reaction.hpp"
#include "nmopt/dealii/scalar_diffusion_reaction.hpp"
#include "nmopt/solvers/reduced_gradient.hpp"

#include <deal.II/base/function_lib.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <cmath>
#include <exception>
#include <iostream>
#include <string>

namespace
{
  using namespace nmopt;
  using Backend = dealii_backend::SerialBackend;
  using Primal = contract::PrimalBlockT<Backend>;
  using Covector = contract::CovectorBlockT<Backend>;

  void
  require_close(const double actual,
                const double expected,
                const double tolerance,
                const std::string &description)
  {
    if (std::abs(actual - expected) > tolerance)
      throw contract::ContractError(description + ": expected " +
                                    std::to_string(expected) + ", got " +
                                    std::to_string(actual));
  }

  void
  require_covector_close(const Covector &   actual,
                         const Covector &   expected,
                         const double       tolerance,
                         const std::string &description)
  {
    contract::require_compatible(actual, expected,
                                 description + " has incompatible layouts");
    for (std::size_t block = 0; block < actual.n_blocks(); ++block)
      {
        dealii::Vector<double> difference = actual.block(block);
        difference.add(-1.0, expected.block(block));
        require_close(difference.l2_norm(), 0.0, tolerance, description);
      }
  }

  Primal
  shifted(Primal value, const Primal &direction, const double step)
  {
    contract::require_compatible(value,
                                 direction,
                                 "deal.II shift has incompatible layouts");
    for (std::size_t block = 0; block < value.n_blocks(); ++block)
      Backend::add_scaled(value.block(block), step, direction.block(block));
    return value;
  }

  template <int dim>
  void
  run_contract_test()
  {
    dealii::Triangulation<dim> triangulation;
    dealii::GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(2);

    const dealii::Functions::ConstantFunction<dim> forcing(1.0);
    const dealii::Functions::ConstantFunction<dim> desired_state(0.25);
    const dealii_backend::ScalarDiffusionReactionModel<dim> model(
      triangulation,
      forcing,
      desired_state,
      1.0,
      0.5,
      0.1,
      1);

    contract::StateControlPartitionT<Backend> partition(model, 0, 1);
    const contract::StateAdjointSolversT<Backend> solvers{
      [&model](const Primal &control) { return model.solve_state(control); },
      [&model](const Primal &full_point, const Covector &state_rhs) {
        return model.solve_adjoint(full_point, state_rhs);
      }};
    const contract::ReducedDTOT<Backend> reduced(model, partition, solvers);

    dealii::Vector<double> control_values(
      partition.control_layout()->dimension(0));
    for (dealii::types::global_dof_index i = 0; i < control_values.size(); ++i)
      control_values[i] = 0.1 + 0.01 * static_cast<double>(i);
    const Primal control(partition.control_layout(), {std::move(control_values)});
    const auto evaluation = reduced.evaluate(control);

    const Covector state_residual = model.residual(evaluation.full_point);
    require_close(state_residual.block(0).l2_norm(),
                  0.0,
                  1e-11,
                  "deal.II state residual");

    dealii_backend::MassMetricSolveParameters metric_parameters;
    metric_parameters.maximum_iterations = 1000;
    metric_parameters.relative_tolerance = 1e-12;
    metric_parameters.absolute_tolerance = 1e-14;
    const auto metric = model.control_l2_metric(metric_parameters);
    dealii::Vector<double> random_covector_values(
      partition.control_layout()->dimension(0));
    for (dealii::types::global_dof_index i = 0;
         i < random_covector_values.size();
         ++i)
      random_covector_values[i] =
        0.17 + 0.03 * static_cast<double>((7 * i) % 11);
    const Covector random_covector(partition.control_layout(),
                                   {std::move(random_covector_values)});
    const Primal recovered_primal = metric.inverse_apply(random_covector);
    const Covector recovered_covector = metric.apply(recovered_primal);
    dealii::Vector<double> inverse_apply_error = recovered_covector.block(0);
    inverse_apply_error.add(-1.0, random_covector.block(0));
    require_close(inverse_apply_error.l2_norm(),
                  0.0,
                  1e-11,
                  "deal.II mass metric inverse/apply consistency");
    dealii::Vector<double> nonidentity_mass_action = recovered_primal.block(0);
    nonidentity_mass_action.add(-1.0, random_covector.block(0));
    contract::require(nonidentity_mass_action.l2_norm() > 1e-3,
                      "deal.II control mass matrix unexpectedly acts as identity");

    dealii::Vector<double> state_tangent(
      model.variable_layout()->dimension(0));
    dealii::Vector<double> control_tangent(
      model.variable_layout()->dimension(1));
    for (dealii::types::global_dof_index i = 0; i < state_tangent.size(); ++i)
      state_tangent[i] = 0.02 * static_cast<double>(i + 1);
    for (dealii::types::global_dof_index i = 0; i < control_tangent.size();
         ++i)
      control_tangent[i] = -0.03 * static_cast<double>(i + 1);
    const Primal tangent(model.variable_layout(),
                         {std::move(state_tangent), std::move(control_tangent)});

    dealii::Vector<double> seed_values(model.test_layout()->dimension(0));
    for (dealii::types::global_dof_index i = 0; i < seed_values.size(); ++i)
      seed_values[i] = 0.04 * static_cast<double>(i + 1);
    const Primal test_seed(model.test_layout(), {std::move(seed_values)});

    const Covector jvp = model.residual_jvp(evaluation.full_point, tangent);
    const Covector vjp = model.residual_vjp(evaluation.full_point, test_seed);
    require_close(contract::pair(jvp, test_seed),
                  contract::pair(vjp, tangent),
                  1e-11,
                  "deal.II residual JVP/VJP pairing");

    constexpr double epsilon = 1e-7;
    const Covector residual_difference =
      model.residual(shifted(evaluation.full_point, tangent, epsilon));
    const Covector residual_at_point = model.residual(evaluation.full_point);
    for (std::size_t block = 0; block < residual_difference.n_blocks(); ++block)
      {
        dealii::Vector<double> finite_difference =
          residual_difference.block(block);
        Backend::add_scaled(finite_difference,
                            -1.0,
                            residual_at_point.block(block));
        Backend::scale(finite_difference, 1.0 / epsilon);
        finite_difference.add(-1.0, jvp.block(block));
        require_close(finite_difference.l2_norm(),
                      0.0,
                      1e-8,
                      "deal.II residual finite-difference JVP");
      }

    dealii::Vector<double> control_direction_values(
      partition.control_layout()->dimension(0));
    for (dealii::types::global_dof_index i = 0;
         i < control_direction_values.size();
         ++i)
      control_direction_values[i] =
        (i % 2 == 0 ? 0.05 : -0.04) * static_cast<double>(i + 1);
    const Primal control_direction(partition.control_layout(),
                                   {std::move(control_direction_values)});

    const Primal metric_direction =
      reduced.gradient_direction(evaluation.reduced_derivative, metric);
    const Covector metric_covector = metric.apply(metric_direction);
    require_close(contract::pair(metric_covector, control_direction),
                  contract::pair(evaluation.reduced_derivative,
                                 control_direction),
                  1e-11,
                  "deal.II mass metric pairing");

    const double reduced_difference =
      reduced.evaluate(shifted(control, control_direction, epsilon)).objective_value -
      evaluation.objective_value;
    require_close(reduced_difference / epsilon,
                  contract::pair(evaluation.reduced_derivative,
                                 control_direction),
                  2e-7,
                  "deal.II reduced DTO derivative");

    nmopt::solvers::ReducedGradientParameters solver_parameters;
    solver_parameters.maximum_iterations = 100;
    solver_parameters.maximum_line_search_trials = 30;
    solver_parameters.gradient_tolerance = 1e-6;
    solver_parameters.initial_step_length = 20.0;
    solver_parameters.armijo_fraction = 1e-4;
    solver_parameters.backtracking_factor = 0.5;
    const nmopt::solvers::ReducedGradientSolverT<Backend> solver(
      reduced, metric, solver_parameters);
    const auto solver_result = solver.solve(control);

    contract::require(
      solver_result.stopping_reason ==
        nmopt::solvers::ReducedGradientStoppingReason::gradient_tolerance,
      "deal.II reduced gradient solver did not reach its tolerance");
    contract::require(solver_result.objective_history.size() > 1,
                      "deal.II reduced gradient solver did not accept an iteration");
    for (std::size_t index = 1;
         index < solver_result.objective_history.size();
         ++index)
      contract::require(solver_result.objective_history[index] <=
                          solver_result.objective_history[index - 1],
                        "deal.II reduced gradient objective history is not monotonic");
    contract::require(solver_result.gradient_norm_history.back() <=
                        solver_parameters.gradient_tolerance,
                      "deal.II reduced gradient final norm exceeds tolerance");
    contract::require(solver_result.state_solve_count ==
                        solver_result.adjoint_solve_count,
                      "deal.II reduced gradient solve counts do not match");
    contract::require(solver_result.line_search_trial_count + 1 ==
                        solver_result.state_solve_count,
                      "deal.II reduced gradient solve count misses a trial evaluation");

    const auto bounds = model.control_l2_box_constraint(-1.0, 0.05);
    dealii::Vector<double> bounded_control_values(
      partition.control_layout()->dimension(0));
    const Primal bounded_control(partition.control_layout(),
                                 {std::move(bounded_control_values)});
    const nmopt::solvers::ReducedGradientSolverT<Backend> projected_solver(
      reduced, metric, bounds, solver_parameters);
    const auto projected_result = projected_solver.solve(bounded_control);

    contract::require(
      projected_result.stopping_reason ==
        nmopt::solvers::ReducedGradientStoppingReason::gradient_tolerance,
      "deal.II projected reduced gradient solver did not reach stationarity");
    contract::require(bounds.is_feasible(projected_result.control),
                      "deal.II projected reduced gradient returned an infeasible control");
    contract::require(projected_result.gradient_norm_history.back() <=
                        solver_parameters.gradient_tolerance,
                      "deal.II projected reduced gradient final norm exceeds tolerance");
    bool upper_bound_is_active = false;
    for (dealii::types::global_dof_index index = 0;
         index < projected_result.control.block(0).size();
         ++index)
      {
        contract::require(projected_result.control.block(0)[index] >= -1.0 &&
                            projected_result.control.block(0)[index] <= 0.05,
                          "deal.II projected reduced gradient left the box");
        upper_bound_is_active = upper_bound_is_active ||
                                projected_result.control.block(0)[index] >=
                                  0.05 - 1e-12;
      }
    contract::require(upper_bound_is_active,
                      "deal.II projected reduced gradient did not reach a bound");
    for (std::size_t index = 1;
         index < projected_result.objective_history.size();
         ++index)
      contract::require(projected_result.objective_history[index] <=
                          projected_result.objective_history[index - 1],
                        "deal.II projected reduced objective is not monotonic");

    const auto specification =
      semantic::v1::make_scalar_diffusion_reaction_problem(true);
    compiler::v1::DealiiDiscretisationPolicy compilation_policy;
    compilation_policy.state_degree = 1;
    const compiler::v1::DealiiCompiler v1_compiler;
    const auto validation = v1_compiler.validate(specification,
                                                 compilation_policy);
    contract::require(validation.valid(),
                      "the canonical v1 problem did not validate for deal.II");

    compiler::v1::DealiiDiscretisationPolicy unsupported_execution =
      compilation_policy;
    unsupported_execution.execution =
      compiler::v1::DealiiDiscretisationPolicy::Execution::matrix_free;
    const auto lowerability_diagnostic =
      v1_compiler.validate(specification, unsupported_execution);
    contract::require(lowerability_diagnostic.has_category(
                        semantic::v1::DiagnosticCategory::lowerability),
                      "v1 compiler did not report an unsupported lowerability");

    auto unsupported_formulation = specification;
    unsupported_formulation.formulation.kind =
      semantic::v1::FormulationKind::all_at_once;
    const auto formulation_diagnostic =
      v1_compiler.validate(unsupported_formulation, compilation_policy);
    contract::require(
      formulation_diagnostic.has_category(
        semantic::v1::DiagnosticCategory::formulation_capability),
      "v1 compiler did not report an unsupported formulation capability");

    const compiler::v1::DealiiDataBindings<dim> data_bindings{
      forcing, desired_state, 1.0, 0.5, 0.1};
    const compiler::v1::CellwiseBoxDataBindings bound_bindings{
      compiler::v1::CellwiseBoundValue{-1.0},
      compiler::v1::CellwiseBoundValue{0.05}};
    const auto compilation = v1_compiler.compile(specification,
                                                  triangulation,
                                                  data_bindings,
                                                  compilation_policy,
                                                  bound_bindings);
    contract::require(compilation.succeeded(),
                      "v1 compiler failed to produce an executable problem");
    const auto &compiled_model = compilation.problem->executable_model();

    const Primal comparison_point =
      shifted(evaluation.full_point, tangent, 0.37);
    const Covector compiled_residual =
      compiled_model.residual(comparison_point);
    require_covector_close(compiled_residual,
                           model.residual(comparison_point),
                           1e-12,
                           "v1 and v0 residual differ");
    require_close(compiled_model.objective(evaluation.full_point),
                  model.objective(evaluation.full_point),
                  1e-12,
                  "v1 and v0 objective differ");
    require_covector_close(compiled_model.objective_derivative(
                             evaluation.full_point),
                           model.objective_derivative(evaluation.full_point),
                           1e-12,
                           "v1 and v0 objective derivative differ");

    contract::StateControlPartitionT<Backend> compiled_partition(compiled_model,
                                                                   0,
                                                                   1);
    const contract::StateAdjointSolversT<Backend> compiled_solvers{
      [&compiled_model](const Primal &compiled_control) {
        return compiled_model.solve_state(compiled_control);
      },
      [&compiled_model](const Primal &full_point, const Covector &state_rhs) {
        return compiled_model.solve_adjoint(full_point, state_rhs);
      }};
    const contract::ReducedDTOT<Backend> compiled_reduced(
      compiled_model, compiled_partition, compiled_solvers);
    const auto compiled_evaluation = compiled_reduced.evaluate(control);
    require_close(compiled_evaluation.objective_value,
                  evaluation.objective_value,
                  1e-12,
                  "v1 and v0 reduced objective differ");
    require_covector_close(compiled_evaluation.reduced_derivative,
                           evaluation.reduced_derivative,
                           1e-12,
                           "v1 and v0 reduced derivative differ");
    const auto *compiled_constraint = compilation.problem->control_constraint();
    contract::require(compiled_constraint != nullptr &&
                        compiled_constraint->is_feasible(bounded_control),
                      "v1 compiler did not preserve the declared box constraint");
  }
} // namespace

int
main()
{
  try
    {
      run_contract_test<2>();
      std::cout << "deal.II diffusion DTO contract test passed\n";
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "deal.II diffusion DTO contract test failed: "
                << exception.what() << '\n';
      return 1;
    }
}
