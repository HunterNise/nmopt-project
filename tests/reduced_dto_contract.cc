#include "nmopt/contract/metric_constraint.hpp"
#include "nmopt/contract/reduced_dto.hpp"
#include "nmopt/reference/linear_quadratic_model.hpp"
#include "nmopt/solvers/reduced_gradient.hpp"
#include "nmopt/solvers/reduced_line_search.hpp"
#include "test_support/contract_errors.hpp"
#include "test_support/scenario_dispatch.hpp"

#include <cmath>
#include <exception>
#include <iostream>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace
{
  using namespace nmopt::contract;
  using nmopt::reference::LinearQuadraticModel;

  struct AlternateDenseBackend
  {
    using Vector = DenseVector;

    static Vector
    zeros(const std::size_t size)
    {
      return Vector(size);
    }

    static std::size_t
    size(const Vector &vector)
    {
      return vector.size();
    }

    static double
    dot(const Vector &left, const Vector &right)
    {
      return nmopt::contract::dot(left, right);
    }

    static void
    add_scaled(Vector &target, const double factor, const Vector &source)
    {
      target.add_scaled(factor, source);
    }

    static void
    scale(Vector &target, const double factor)
    {
      target.scale(factor);
    }
  };

  void
  require_close(const double actual,
                const double expected,
                const double tolerance,
                const std::string &description)
  {
    if (std::abs(actual - expected) > tolerance)
      throw ContractError(description + ": expected " +
                          std::to_string(expected) + ", got " +
                          std::to_string(actual));
  }

  PrimalBlock
  shifted(PrimalBlock value, const PrimalBlock &direction, const double step)
  {
    require_compatible(value, direction, "Shift has incompatible primal layouts");
    for (std::size_t block = 0; block < value.n_blocks(); ++block)
      value.add_scaled_block(block, step, direction.block(block));
    return value;
  }

  class RecordingConstraint final : public Constraint
  {
  public:
    explicit RecordingConstraint(const CellwiseBoxConstraint &constraint)
      : constraint_(constraint)
    {}

    const LayoutPtr &
    layout() const override
    {
      return constraint_.layout();
    }

    bool
    is_feasible(const PrimalBlock &primal) const override
    {
      return constraint_.is_feasible(primal);
    }

    bool
    supports_projection_in(const Metric &metric) const override
    {
      return constraint_.supports_projection_in(metric);
    }

    PrimalBlock
    project_in(const PrimalBlock &primal, const Metric &metric) const override
    {
      PrimalBlock projected = constraint_.project_in(primal, metric);
      projected_controls_.push_back(projected);
      return projected;
    }

    const std::vector<PrimalBlock> &
    projected_controls() const
    {
      return projected_controls_;
    }

  private:
    const CellwiseBoxConstraint &  constraint_;
    mutable std::vector<PrimalBlock> projected_controls_;
  };

  class NonDiagonalMetric final : public Metric
  {
  public:
    explicit NonDiagonalMetric(LayoutPtr layout)
      : layout_(std::move(layout))
    {}

    const std::string &
    id() const override
    {
      static const std::string id = "l2_cellwise";
      return id;
    }

    const LayoutPtr &
    layout() const override
    {
      return layout_;
    }

    CovectorBlock
    apply(const PrimalBlock &primal) const override
    {
      require(primal.layout()->compatible_with(*layout_),
              "Non-diagonal metric primal has an incompatible layout");
      const auto &values = primal.block(0);
      return CovectorBlock(
        layout_, {DenseVector{2.0 * values[0] + values[1],
                              values[0] + 2.0 * values[1]}});
    }

    PrimalBlock
    inverse_apply(const CovectorBlock &covector) const override
    {
      require(covector.layout()->compatible_with(*layout_),
              "Non-diagonal metric covector has an incompatible layout");
      const auto &values = covector.block(0);
      return PrimalBlock(
        layout_, {DenseVector{(2.0 * values[0] - values[1]) / 3.0,
                              (-values[0] + 2.0 * values[1]) / 3.0}});
    }

  private:
    LayoutPtr layout_;
  };

  void
  test_block_layout_invariant()
  {
    using BlockAccess = decltype(std::declval<PrimalBlock &>().block(0));
    static_assert(std::is_same_v<BlockAccess, const DenseVector &>,
                  "BlockValues public block access must be read-only");

    const auto layout = std::make_shared<const BlockLayout>(
      "layout_invariant",
      std::vector<SpaceId>{{"state"}},
      std::vector<std::size_t>{2});
    PrimalBlock   primal(layout, {DenseVector{1.0, 2.0}});
    CovectorBlock covector(layout, {DenseVector{3.0, 4.0}});

    nmopt::test_support::require_contract_error(
      [&primal]() {
        primal.add_scaled_block(0, 1.0, DenseVector{1.0, 2.0, 3.0});
      },
      "BlockValues update vector dimension does not match layout",
      "dimension-changing primal block update");
    nmopt::test_support::require_contract_error(
      [&covector]() {
        covector.add_scaled_block(0, 1.0, DenseVector{1.0, 2.0, 3.0});
      },
      "BlockValues update vector dimension does not match layout",
      "dimension-changing covector block update");
    require_close(pair(covector, primal),
                  11.0,
                  1e-15,
                  "Rejected block updates preserve pairing");

    primal.add_scaled_block(0, 1.0, DenseVector{-1.0, 0.5});
    covector.scale_block(0, 0.5);
    require_close(pair(covector, primal),
                  5.0,
                  1e-15,
                  "Checked block algebra preserves pairing");
  }

  void
  test_v0_contract()
  {
    const LinearQuadraticModel model(
      DenseMatrix(2, 2, {4.0, -1.0, -1.0, 3.0}),
      DenseMatrix(2, 2, {1.0, 0.5, -0.25, 2.0}),
      DenseVector{1.0, -0.5},
      DenseMatrix(2, 2, {1.0, 0.0, 0.5, 1.0}),
      DenseVector{0.25, -1.0},
      DenseVector{1.5, 0.75},
      DenseVector{2.0, 3.0},
      0.4);

    const PrimalBlock point(
      model.variable_layout(),
      {DenseVector{0.2, -0.3}, DenseVector{0.5, -0.4}});
    const PrimalBlock tangent(
      model.variable_layout(),
      {DenseVector{-0.7, 0.25}, DenseVector{0.3, 0.8}});
    const PrimalBlock test_seed(model.test_layout(), {DenseVector{0.6, -1.1}});

    const CovectorBlock residual_direction = model.residual_jvp(point, tangent);
    const CovectorBlock residual_pullback = model.residual_vjp(point, test_seed);
    require_close(pair(residual_direction, test_seed),
                  pair(residual_pullback, tangent),
                  1e-13,
                  "Residual JVP/VJP pairing");

    constexpr double epsilon = 1e-7;
    CovectorBlock finite_difference =
      model.residual(shifted(point, tangent, epsilon));
    const CovectorBlock residual_at_point = model.residual(point);
    for (std::size_t block = 0; block < finite_difference.n_blocks(); ++block)
      {
        finite_difference.add_scaled_block(
          block, -1.0, residual_at_point.block(block));
        finite_difference.scale_block(block, 1.0 / epsilon);
        for (std::size_t entry = 0;
             entry < finite_difference.block(block).size();
             ++entry)
          require_close(finite_difference.block(block)[entry],
                        residual_direction.block(block)[entry],
                        1e-9,
                        "Residual finite-difference JVP");
      }

    const CovectorBlock objective_derivative =
      model.objective_derivative(point);
    const double objective_difference =
      model.objective(shifted(point, tangent, epsilon)) - model.objective(point);
    require_close(objective_difference / epsilon,
                  pair(objective_derivative, tangent),
                  2e-7,
                  "Objective directional derivative");

    nmopt::test_support::require_contract_error(
      [&model]() { (void)StateControlPartition(model, 0, 0); },
      "State and control blocks must be distinct",
      "invalid reduced partition");

    const StateControlPartition partition(model, 0, 1);
    const StateAdjointSolvers solvers{
      [&model](const PrimalBlock &control) { return model.solve_state(control); },
      [&model](const PrimalBlock &full_point, const CovectorBlock &state_rhs) {
        return model.solve_adjoint(full_point, state_rhs);
      }};
    nmopt::test_support::require_contract_error(
      [&model, &partition, &solvers]() {
        (void)ReducedDTO(
          model, partition, StateAdjointSolvers{solvers.solve_state, {}});
      },
      "Reduced DTO requires an adjoint solve operation",
      "missing reduced adjoint callback");
    const ReducedDTO reduced(model, partition, solvers);

    const PrimalBlock control(partition.control_layout(), {DenseVector{0.4, -0.3}});
    const ReducedEvaluation evaluation = reduced.evaluate(control);

    const CovectorBlock state_residual = model.residual(evaluation.full_point);
    require_close(dot(state_residual.block(0), state_residual.block(0)),
                  0.0,
                  1e-24,
                  "State solve residual");

    const PrimalBlock control_direction(
      partition.control_layout(), {DenseVector{-0.2, 0.7}});
    const double reduced_difference =
      reduced.evaluate(shifted(control, control_direction, epsilon)).objective_value -
      evaluation.objective_value;
    require_close(reduced_difference / epsilon,
                  pair(evaluation.reduced_derivative, control_direction),
                  2e-7,
                  "Reduced DTO derivative");

    nmopt::test_support::require_contract_error(
      [&partition]() {
        (void)DiagonalMetric(
          "invalid_metric",
          partition.control_layout(),
          {DenseVector{2.0, 0.0}});
      },
      "Metric diagonal must be strictly positive",
      "nonpositive metric diagonal");
    const DiagonalMetric metric(
      "l2_cellwise", partition.control_layout(), {DenseVector{2.0, 5.0}});
    const PrimalBlock direction =
      reduced.gradient_direction(evaluation.reduced_derivative, metric);
    const CovectorBlock reconstructed = metric.apply(direction);
    require_close(pair(reconstructed, control_direction),
                  pair(evaluation.reduced_derivative, control_direction),
                  1e-13,
                  "Metric inverse/apply consistency");
    const auto search_direction =
      nmopt::solvers::make_steepest_descent_direction(
        evaluation.reduced_derivative, metric);
    require(search_direction.directional_derivative < 0.0,
            "Steepest-descent search direction is not descending");
    require_close(search_direction.gradient_norm,
                  std::sqrt(pair(reconstructed, direction)),
                  1e-15,
                  "Steepest-descent direction metric norm");
    require_close(
      pair(evaluation.reduced_derivative, search_direction.direction),
      -pair(evaluation.reduced_derivative, direction),
      1e-15,
      "Steepest-descent direction sign");

    const ReducedHessian &hessian = model;
    const CovectorBlock hessian_action =
      hessian.apply(control, control_direction);
    CovectorBlock reduced_derivative_difference =
      reduced.evaluate(shifted(control, control_direction, epsilon))
        .reduced_derivative;
    for (std::size_t block = 0;
         block < reduced_derivative_difference.n_blocks();
         ++block)
      {
        reduced_derivative_difference.add_scaled_block(
          block, -1.0, evaluation.reduced_derivative.block(block));
        reduced_derivative_difference.scale_block(block, 1.0 / epsilon);
        for (std::size_t entry = 0;
             entry < reduced_derivative_difference.block(block).size();
             ++entry)
          require_close(reduced_derivative_difference.block(block)[entry],
                        hessian_action.block(block)[entry],
                        2e-7,
                        "Reduced Hessian finite-difference action");
      }
    const PrimalBlock second_hessian_direction(
      partition.control_layout(), {DenseVector{0.35, -0.2}});
    const CovectorBlock second_hessian_action =
      hessian.apply(control, second_hessian_direction);
    require_close(pair(hessian_action, second_hessian_direction),
                  pair(second_hessian_action, control_direction),
                  1e-13,
                  "Reduced Hessian symmetry");

    const auto build_trial =
      [&control, &search_direction](const double step) {
        return shifted(control, search_direction.direction, step);
      };
    const auto evaluate_trial =
      [&reduced](const PrimalBlock &trial_control) {
        return reduced.evaluate(trial_control);
      };

    nmopt::solvers::ArmijoLineSearchParameters armijo_parameters;
    armijo_parameters.maximum_trials = 20;
    armijo_parameters.initial_step_length = 10.0;
    const nmopt::solvers::ArmijoLineSearchPolicy armijo_policy(
      armijo_parameters);
    const auto armijo_result = armijo_policy.search(control,
                                                    evaluation,
                                                    search_direction,
                                                    build_trial,
                                                    evaluate_trial);
    require(armijo_result.accepted() && armijo_result.trial_count > 0,
            "Armijo line search did not accept a trial");
    require(armijo_result.evaluation.objective_value <=
              evaluation.objective_value,
            "Armijo line search accepted an objective increase");

    const nmopt::solvers::ExactQuadraticLineSearchPolicy exact_policy(hessian);
    const auto exact_result = exact_policy.search(control,
                                                  evaluation,
                                                  search_direction,
                                                  build_trial,
                                                  evaluate_trial);
    require(exact_result.accepted() && exact_result.trial_count == 1 &&
              exact_result.hessian_action_count == 1,
            "Exact quadratic line search did not accept its Hessian step");
    require(exact_result.evaluation.objective_value <=
              evaluation.objective_value,
            "Exact quadratic line search accepted an objective increase");

    nmopt::solvers::WolfeLineSearchParameters wolfe_parameters;
    wolfe_parameters.maximum_trials = 30;
    wolfe_parameters.initial_step_length = 10.0;
    const nmopt::solvers::WolfeLineSearchPolicy wolfe_policy(
      wolfe_parameters);
    const auto wolfe_result = wolfe_policy.search(control,
                                                  evaluation,
                                                  search_direction,
                                                  build_trial,
                                                  evaluate_trial);
    require(wolfe_result.accepted() && wolfe_result.trial_count > 0,
            "Wolfe line search did not accept a trial");

    nmopt::solvers::ArmijoLineSearchParameters actual_displacement_parameters;
    actual_displacement_parameters.maximum_trials = 1;
    const nmopt::solvers::ArmijoLineSearchPolicy actual_displacement_policy(
      actual_displacement_parameters);
    const auto scaled_trial_builder =
      [&control, &search_direction](const double step) {
        return shifted(control, search_direction.direction, 0.25 * step);
      };
    const auto actual_displacement_result =
      actual_displacement_policy.search(control,
                                        evaluation,
                                        search_direction,
                                        scaled_trial_builder,
                                        evaluate_trial);
    require(actual_displacement_result.accepted() &&
              actual_displacement_result.trial_count == 1,
            "Armijo line search did not use the actual trial displacement");

    const nmopt::solvers::ExactQuadraticLineSearchPolicy
      missing_exact_hessian_policy;
    nmopt::test_support::require_contract_error(
      [&missing_exact_hessian_policy, &control, &evaluation, &search_direction,
       &build_trial, &evaluate_trial]() {
        (void)missing_exact_hessian_policy.search(control,
                                                  evaluation,
                                                  search_direction,
                                                  build_trial,
                                                  evaluate_trial);
      },
      "Exact quadratic line search requires a reduced Hessian capability",
      "exact line search missing Hessian capability");

    const DiagonalMetric cg_metric(
      "cg_metric", partition.control_layout(), {DenseVector{1.0, 1.0}});
    using CGPolicy =
      nmopt::solvers::NonlinearConjugateGradientDirectionPolicyDense;
    CGPolicy cg_policy({2, 1e-14});
    const CovectorBlock first_cg_derivative(
      partition.control_layout(), {DenseVector{1.0, 0.0}});
    const CovectorBlock second_cg_derivative(
      partition.control_layout(), {DenseVector{1.0, 1.0}});
    const CovectorBlock third_cg_derivative(
      partition.control_layout(), {DenseVector{0.5, 0.5}});
    const auto first_cg_direction =
      cg_policy.next(control, first_cg_derivative, cg_metric);
    const auto second_cg_direction =
      cg_policy.next(control, second_cg_derivative, cg_metric);
    require_close(first_cg_direction.direction.block(0)[0],
                  -1.0,
                  1e-15,
                  "Nonlinear CG initial direction");
    require_close(second_cg_direction.direction.block(0)[0],
                  -2.0,
                  1e-15,
                  "Nonlinear CG Polak-Ribiere direction");
    require_close(second_cg_direction.direction.block(0)[1],
                  -1.0,
                  1e-15,
                  "Nonlinear CG Polak-Ribiere update");
    require(second_cg_direction.directional_derivative < 0.0,
            "Nonlinear CG Polak-Ribiere direction is not descending");
    (void)cg_policy.next(control, third_cg_derivative, cg_metric);
    require(cg_policy.restart_count() == 1,
            "Nonlinear CG did not perform its configured periodic restart");
    CGPolicy negative_beta_policy({10, 1e-14});
    (void)negative_beta_policy.next(control, first_cg_derivative, cg_metric);
    (void)negative_beta_policy.next(
      control,
      CovectorBlock(partition.control_layout(), {DenseVector{0.5, 0.0}}),
      cg_metric);
    require(negative_beta_policy.restart_count() == 1,
            "Nonlinear CG did not restart after a nonpositive coefficient");
    nmopt::test_support::require_contract_error(
      [&cg_policy, &control]() {
        const auto incompatible_layout = std::make_shared<const BlockLayout>(
          "incompatible_cg_layout",
          std::vector<SpaceId>{{"control"}},
          std::vector<std::size_t>{3});
        const DiagonalMetric incompatible_metric(
          "incompatible_cg_metric",
          incompatible_layout,
          {DenseVector{1.0, 1.0, 1.0}});
        const CovectorBlock incompatible_derivative(
          incompatible_layout,
          {DenseVector{1.0, 1.0, 1.0}});
        (void)cg_policy.next(control, incompatible_derivative, incompatible_metric);
      },
      "Nonlinear CG derivative history has an incompatible layout",
      "nonlinear CG history layout mismatch");

    using BfgsPolicy =
      nmopt::solvers::LimitedMemoryBfgsDirectionPolicyDense;
    BfgsPolicy bfgs_policy({2, 1e-14});
    const PrimalBlock bfgs_control_0(
      partition.control_layout(), {DenseVector{0.0, 0.0}});
    const PrimalBlock bfgs_control_1(
      partition.control_layout(), {DenseVector{1.0, 0.0}});
    const CovectorBlock bfgs_derivative_0(
      partition.control_layout(), {DenseVector{1.0, 0.0}});
    const CovectorBlock bfgs_derivative_1(
      partition.control_layout(), {DenseVector{2.0, 0.0}});
    (void)bfgs_policy.next(bfgs_control_0, bfgs_derivative_0, cg_metric);
    const auto bfgs_direction =
      bfgs_policy.next(bfgs_control_1, bfgs_derivative_1, cg_metric);
    require(bfgs_policy.history_size() == 1,
            "L-BFGS did not retain a valid secant pair");
    require(bfgs_policy.last_update_status() ==
              nmopt::solvers::LimitedMemoryBfgsUpdateStatus::accepted_pair,
            "L-BFGS did not report an accepted secant pair");
    require(bfgs_direction.directional_derivative < 0.0,
            "L-BFGS direction is not descending");
    (void)bfgs_policy.next(bfgs_control_1, bfgs_derivative_1, cg_metric);
    require(bfgs_policy.history_size() == 0 &&
              bfgs_policy.reset_count() == 1,
            "L-BFGS did not reset after failed curvature");
    require(bfgs_policy.last_update_status() ==
              nmopt::solvers::LimitedMemoryBfgsUpdateStatus::curvature_reset,
            "L-BFGS curvature reset was not reported");

    nmopt::solvers::ReducedNewtonParameters newton_parameters;
    newton_parameters.maximum_inner_iterations = 10;
    newton_parameters.relative_tolerance = 1e-12;
    newton_parameters.absolute_tolerance = 1e-14;
    newton_parameters.curvature_tolerance = 1e-14;
    nmopt::solvers::ReducedSolverParameters newton_solver_parameters;
    newton_solver_parameters.gradient_tolerance = 1e-8;
    newton_solver_parameters.initial_step_length = 10.0;
    const nmopt::solvers::NewtonDirectionPolicyDense newton_policy(
      hessian, newton_parameters);
    const nmopt::solvers::ReducedNewtonSolver newton_solver(
      reduced, metric, newton_solver_parameters, newton_policy);
    const auto newton_result = newton_solver.solve(control);
    require(newton_result.stopping_reason ==
              nmopt::solvers::ReducedGradientStoppingReason::gradient_tolerance,
            "Reduced Newton solver did not reach its tolerance");
    require(newton_result.hessian_action_count > 0,
            "Reduced Newton solver did not report Hessian actions");
    require(newton_result.step_length_history.size() ==
              newton_result.accepted_iterations,
            "Reduced Newton step history does not match accepted iterations");
    require(newton_result.metric_solve_count >
              newton_result.hessian_action_count,
            "Reduced Newton did not report metric preconditioning actions");

    const nmopt::solvers::NewtonDirectionPolicyDense missing_hessian_policy;
    const nmopt::solvers::ReducedNewtonSolver missing_hessian_solver(
      reduced, metric, newton_solver_parameters, missing_hessian_policy);
    nmopt::test_support::require_contract_error(
      [&missing_hessian_solver, &control]() {
        (void)missing_hessian_solver.solve(control);
      },
      "Newton direction requires a reduced Hessian capability",
      "Newton missing Hessian capability");

    nmopt::test_support::require_contract_error(
      [&partition, &metric]() {
        (void)CellwiseBoxConstraint(
          partition.control_layout(),
          {DenseVector{0.0, 0.6}},
          {DenseVector{1.0, 0.5}},
          metric);
      },
      "Cellwise box lower bound exceeds upper bound",
      "reversed cellwise bounds");
    const CellwiseBoxConstraint bounds(
      partition.control_layout(), {DenseVector{-0.1, -0.2}},
      {DenseVector{0.6, 0.5}}, metric);
    const PrimalBlock projected = bounds.project_in(
      PrimalBlock(partition.control_layout(), {DenseVector{-0.4, 0.8}}), metric);
    require(bounds.is_feasible(projected),
            "Cellwise L2 projection must return a feasible control");
    require_close(projected.block(0)[0], -0.1, 1e-15,
                  "Cellwise box lower projection");
    require_close(projected.block(0)[1], 0.5, 1e-15,
                  "Cellwise box upper projection");

    nmopt::solvers::ReducedGradientParameters solver_parameters;
    solver_parameters.maximum_iterations = 100;
    solver_parameters.maximum_line_search_trials = 20;
    solver_parameters.gradient_tolerance = 1e-8;
    solver_parameters.initial_step_length = 10.0;
    solver_parameters.armijo_fraction = 1e-4;
    solver_parameters.backtracking_factor = 0.5;
    const nmopt::solvers::ReducedGradientSolver solver(
      reduced, metric, solver_parameters);
    const auto solver_result = solver.solve(
      PrimalBlock(partition.control_layout(), {DenseVector{1.0, -1.0}}));

    require(solver_result.stopping_reason ==
              nmopt::solvers::ReducedGradientStoppingReason::gradient_tolerance,
            "Dense reduced gradient solver did not reach its tolerance");
    require(solver_result.objective_history.size() > 1,
            "Dense reduced gradient solver did not accept an iteration");
    for (std::size_t index = 1;
         index < solver_result.objective_history.size();
         ++index)
      require(solver_result.objective_history[index] <=
                solver_result.objective_history[index - 1],
              "Dense reduced gradient objective history is not monotonic");
    require(solver_result.gradient_norm_history.back() <=
              solver_parameters.gradient_tolerance,
            "Dense reduced gradient final norm exceeds the configured tolerance");
    require(solver_result.state_solve_count == solver_result.adjoint_solve_count,
            "Dense reduced gradient solve counts do not match");
    require(solver_result.line_search_trial_count + 1 ==
              solver_result.state_solve_count,
            "Dense reduced gradient solve count misses a trial evaluation");
    require(solver_result.state_solve_count >
              solver_result.objective_history.size(),
            "Dense reduced gradient test did not exercise Armijo backtracking");
    require(solver_result.metric_solve_count ==
              solver_result.gradient_norm_history.size(),
            "Dense reduced gradient metric solve count does not match direction evaluations");
    require(solver_result.step_length_history.size() ==
              solver_result.accepted_iterations,
            "Dense reduced gradient step history does not match accepted iterations");
    require(solver_result.objective_change_history.size() ==
              solver_result.accepted_iterations,
            "Dense reduced gradient objective-change history does not match accepted iterations");
    require(solver_result.relative_gradient_norm_history.size() ==
              solver_result.gradient_norm_history.size(),
            "Dense reduced gradient relative-norm history does not match direction evaluations");
    require(solver_result.step_norm_history.size() ==
              solver_result.accepted_iterations,
            "Dense reduced gradient step-norm history does not match accepted iterations");
    require(solver_result.final_evaluation.state_solve.converged() &&
              solver_result.final_evaluation.adjoint_solve.converged(),
            "Dense reduced gradient result does not retain the final evaluation");
    require_close(solver_result.final_evaluation.objective_value,
                  solver_result.objective_history.back(),
                  1e-14,
                  "Dense reduced gradient final objective");
    require(solver_result.final_evaluation.reduced_derivative.layout()->compatible_with(
              *solver_result.control.layout()),
            "Dense reduced gradient final covector has the wrong layout");
    for (std::size_t index = 0;
         index < solver_result.accepted_iterations;
         ++index)
      {
        require(solver_result.step_length_history[index] > 0.0,
                "Dense reduced gradient accepted a nonpositive step");
        require(solver_result.step_norm_history[index] > 0.0,
                "Dense reduced gradient accepted a zero metric step");
        require_close(solver_result.objective_change_history[index],
                      solver_result.objective_history[index + 1] -
                        solver_result.objective_history[index],
                        1e-14,
                        "Dense reduced gradient objective-change history");
      }

    nmopt::solvers::ReducedGradientParameters relative_stopping_parameters =
      solver_parameters;
    relative_stopping_parameters.gradient_tolerance = 1e-30;
    relative_stopping_parameters.relative_gradient_tolerance = 0.5;
    const nmopt::solvers::ReducedGradientSolver relative_stopping_solver(
      reduced, metric, relative_stopping_parameters);
    const auto relative_stopping_result = relative_stopping_solver.solve(
      PrimalBlock(partition.control_layout(), {DenseVector{1.0, -1.0}}));
    require(relative_stopping_result.stopping_reason ==
              nmopt::solvers::ReducedGradientStoppingReason::relative_gradient_tolerance,
            "Reduced solver did not stop on the relative gradient tolerance");
    require(relative_stopping_result.relative_gradient_norm_history.back() <=
              relative_stopping_parameters.relative_gradient_tolerance,
            "Reduced solver relative gradient tolerance was not enforced");

    nmopt::solvers::ReducedGradientParameters objective_stopping_parameters =
      solver_parameters;
    objective_stopping_parameters.gradient_tolerance = 1e-30;
    objective_stopping_parameters.objective_change_tolerance = 1e6;
    const nmopt::solvers::ReducedGradientSolver objective_stopping_solver(
      reduced, metric, objective_stopping_parameters);
    const auto objective_stopping_result = objective_stopping_solver.solve(
      PrimalBlock(partition.control_layout(), {DenseVector{1.0, -1.0}}));
    require(objective_stopping_result.stopping_reason ==
              nmopt::solvers::ReducedGradientStoppingReason::objective_change_tolerance,
            "Reduced solver did not stop on the objective-change tolerance");
    require(objective_stopping_result.accepted_iterations > 0,
            "Objective-change stopping did not retain an accepted trial");

    nmopt::solvers::ReducedGradientParameters step_stopping_parameters =
      solver_parameters;
    step_stopping_parameters.gradient_tolerance = 1e-30;
    step_stopping_parameters.step_tolerance = 1e6;
    const nmopt::solvers::ReducedGradientSolver step_stopping_solver(
      reduced, metric, step_stopping_parameters);
    const auto step_stopping_result = step_stopping_solver.solve(
      PrimalBlock(partition.control_layout(), {DenseVector{1.0, -1.0}}));
    require(step_stopping_result.stopping_reason ==
              nmopt::solvers::ReducedGradientStoppingReason::step_tolerance,
            "Reduced solver did not stop on the step tolerance");
    require(step_stopping_result.step_norm_history.back() <=
              step_stopping_parameters.step_tolerance,
            "Reduced solver step tolerance was not enforced");

    const nmopt::solvers::ReducedConjugateGradientSolver cg_solver(
      reduced, metric, solver_parameters);
    const auto cg_result = cg_solver.solve(
      PrimalBlock(partition.control_layout(), {DenseVector{1.0, -1.0}}));
    require(cg_result.stopping_reason ==
              nmopt::solvers::ReducedGradientStoppingReason::gradient_tolerance,
            "Dense nonlinear CG solver did not reach its tolerance");
    require(cg_result.metric_solve_count ==
              cg_result.gradient_norm_history.size(),
            "Dense nonlinear CG metric solve count does not match direction evaluations");
    require(cg_result.step_length_history.size() ==
              cg_result.accepted_iterations,
            "Dense nonlinear CG step history does not match accepted iterations");
    for (std::size_t index = 1;
         index < cg_result.objective_history.size();
         ++index)
      require(cg_result.objective_history[index] <=
                cg_result.objective_history[index - 1],
              "Dense nonlinear CG objective history is not monotonic");

    const nmopt::solvers::ReducedLimitedMemoryBfgsSolver lbfgs_solver(
      reduced, metric, solver_parameters);
    const auto lbfgs_result = lbfgs_solver.solve(
      PrimalBlock(partition.control_layout(), {DenseVector{1.0, -1.0}}));
    require(lbfgs_result.stopping_reason ==
              nmopt::solvers::ReducedGradientStoppingReason::gradient_tolerance,
            "Dense L-BFGS solver did not reach its tolerance");
    require(lbfgs_result.step_length_history.size() ==
              lbfgs_result.accepted_iterations,
            "Dense L-BFGS step history does not match accepted iterations");
    require(lbfgs_result.direction_reset_count <=
              lbfgs_result.accepted_iterations,
            "Dense L-BFGS direction reset count exceeds accepted iterations");
    for (std::size_t index = 1;
         index < lbfgs_result.objective_history.size();
         ++index)
      require(lbfgs_result.objective_history[index] <=
                lbfgs_result.objective_history[index - 1],
              "Dense L-BFGS objective history is not monotonic");

    const nmopt::solvers::NewtonDirectionPolicyDense exact_newton_direction(
      hessian, newton_parameters);
    const nmopt::solvers::ExactQuadraticLineSearchPolicy exact_newton_line_search(
      hessian);
    const nmopt::solvers::ReducedExactNewtonSolver exact_newton_solver(
      reduced,
      metric,
      solver_parameters,
      exact_newton_direction,
      exact_newton_line_search);
    const auto exact_newton_result = exact_newton_solver.solve(
      PrimalBlock(partition.control_layout(), {DenseVector{1.0, -1.0}}));
    require(exact_newton_result.stopping_reason ==
              nmopt::solvers::ReducedGradientStoppingReason::gradient_tolerance,
            "Dense exact Newton solver did not reach its tolerance");
    require(exact_newton_result.line_search_trial_count ==
              exact_newton_result.accepted_iterations,
            "Dense exact Newton did not use one exact trial per accepted iteration");
    require(exact_newton_result.state_solve_count ==
              exact_newton_result.line_search_trial_count + 1,
            "Dense exact Newton solve count misses a trial evaluation");
    require(exact_newton_result.hessian_action_count >
              exact_newton_result.accepted_iterations,
            "Dense exact Newton did not report line-search Hessian actions");

    nmopt::solvers::WolfeLineSearchParameters wolfe_solver_line_parameters;
    wolfe_solver_line_parameters.maximum_trials = 30;
    wolfe_solver_line_parameters.initial_step_length = 10.0;
    const nmopt::solvers::SteepestDescentDirectionPolicy wolfe_direction_policy;
    const nmopt::solvers::WolfeLineSearchPolicy wolfe_solver_line_search(
      wolfe_solver_line_parameters);
    const nmopt::solvers::ReducedWolfeGradientSolver wolfe_solver(
      reduced,
      metric,
      solver_parameters,
      wolfe_direction_policy,
      wolfe_solver_line_search);
    const auto wolfe_solver_result = wolfe_solver.solve(
      PrimalBlock(partition.control_layout(), {DenseVector{1.0, -1.0}}));
    require(wolfe_solver_result.stopping_reason ==
              nmopt::solvers::ReducedGradientStoppingReason::gradient_tolerance,
            "Dense Wolfe solver did not reach its tolerance");
    require(wolfe_solver_result.hessian_action_count == 0,
            "Dense Wolfe solver unexpectedly reported Hessian actions");
    require(wolfe_solver_result.state_solve_count ==
              wolfe_solver_result.line_search_trial_count + 1,
            "Dense Wolfe solve count misses a trial evaluation");
    for (std::size_t index = 1;
         index < wolfe_solver_result.objective_history.size();
         ++index)
      require(wolfe_solver_result.objective_history[index] <=
                wolfe_solver_result.objective_history[index - 1],
              "Dense Wolfe objective history is not monotonic");

    const CellwiseBoxConstraint projected_bounds(
      partition.control_layout(), {DenseVector{-1.0, -0.2}},
      {DenseVector{1.0, 0.5}}, metric);
    const RecordingConstraint recording_bounds(projected_bounds);
    const nmopt::solvers::ReducedGradientSolver projected_solver(
      reduced, metric, recording_bounds, solver_parameters);
    nmopt::test_support::require_contract_error(
      [&projected_solver, &partition]() {
        (void)projected_solver.solve(
          PrimalBlock(partition.control_layout(),
                      {DenseVector{0.0, 0.6}}));
      },
      "Projected reduced gradient requires a feasible initial control",
      "infeasible projected-solver initial control");
    const auto projected_result = projected_solver.solve(
      PrimalBlock(partition.control_layout(), {DenseVector{0.4, 0.4}}));

    require(projected_result.stopping_reason ==
              nmopt::solvers::ReducedGradientStoppingReason::gradient_tolerance,
            "Projected reduced gradient solver did not reach stationarity");
    require(projected_bounds.is_feasible(projected_result.control),
            "Projected reduced gradient solver returned an infeasible control");
    require_close(projected_result.control.block(0)[1],
                  -0.2,
                  1e-12,
                  "Projected reduced gradient did not reach the lower bound");
    require(projected_result.gradient_norm_history.back() <=
              solver_parameters.gradient_tolerance,
            "Projected reduced gradient final norm exceeds the tolerance");
    require(projected_result.metric_solve_count ==
              projected_result.gradient_norm_history.size(),
            "Projected reduced gradient metric solve count does not match direction evaluations");
    require(projected_result.step_length_history.size() ==
              projected_result.accepted_iterations,
            "Projected reduced gradient step history does not match accepted iterations");
    require(!recording_bounds.projected_controls().empty(),
            "Projected reduced gradient did not request a projection");
    for (const PrimalBlock &projected_control :
         recording_bounds.projected_controls())
      require(projected_bounds.is_feasible(projected_control),
              "Projected reduced gradient generated an infeasible iterate");
    for (std::size_t index = 1;
         index < projected_result.objective_history.size();
         ++index)
      require(projected_result.objective_history[index] <=
                projected_result.objective_history[index - 1],
              "Projected reduced gradient objective history is not monotonic");
  }

  void
  test_owned_reduced_service_lifetime()
  {
    auto detached = []() {
      auto model = std::make_shared<LinearQuadraticModel>(
        DenseMatrix(2, 2, {4.0, -1.0, -1.0, 3.0}),
        DenseMatrix(2, 2, {1.0, 0.5, -0.25, 2.0}),
        DenseVector{1.0, -0.5},
        DenseMatrix(2, 2, {1.0, 0.0, 0.5, 1.0}),
        DenseVector{0.25, -1.0},
        DenseVector{1.5, 0.75},
        DenseVector{2.0, 3.0},
        0.4);
      const auto *model_view = model.get();
      StateControlPartition partition(*model_view, 0, 1);
      StateAdjointSolvers solvers{
        [model_view](const PrimalBlock &control) {
          return model_view->solve_state(control);
        },
        [model_view](const PrimalBlock &full_point,
                     const CovectorBlock &state_rhs) {
          return model_view->solve_adjoint(full_point, state_rhs);
        }};
      PrimalBlock control(partition.control_layout(),
                          {DenseVector{0.4, -0.3}});
      struct DetachedService
      {
        ReducedDTO  reduced;
        PrimalBlock control;
      };
      return DetachedService{
        ReducedDTO(std::move(model),
                   std::move(partition),
                   std::move(solvers),
                   std::make_shared<const int>(7)),
        std::move(control)};
    }();

    const auto evaluation = detached.reduced.evaluate(detached.control);
    require(evaluation.state_solve.converged() &&
              evaluation.adjoint_solve.converged(),
            "Detached owned reduced service lost its solve callbacks");
  }

  void
  test_backend_parameterisation()
  {
    const auto layout = std::make_shared<const BlockLayout>(
      "alternate",
      std::vector<SpaceId>{{"alternate_control"}},
      std::vector<std::size_t>{2});

    const PrimalBlockT<AlternateDenseBackend> primal(
      layout, {DenseVector{1.0, -2.0}});
    const CovectorBlockT<AlternateDenseBackend> covector(
      layout, {DenseVector{3.0, 4.0}});
    require_close(pair(covector, primal),
                  -5.0,
                  1e-15,
                  "Backend-parametric pairing");

    CovectorBlockT<AlternateDenseBackend> difference =
      subtract(covector,
               CovectorBlockT<AlternateDenseBackend>(
                 layout, {DenseVector{1.0, -1.0}}));
    require_close(difference.block(0)[0],
                  2.0,
                  1e-15,
                  "Backend-parametric covector subtraction");
    require_close(difference.block(0)[1],
                  5.0,
                  1e-15,
                  "Backend-parametric covector subtraction");
  }

  void
  test_projection_compatibility()
  {
    const auto layout = std::make_shared<const BlockLayout>(
      "projection_compatibility",
      std::vector<SpaceId>{{"control"}},
      std::vector<std::size_t>{2});
    const DiagonalMetric trusted_metric(
      "l2_cellwise", layout, {DenseVector{2.0, 2.0}});
    const CellwiseBoxConstraint bounds(
      layout,
      {DenseVector{0.0, 0.0}},
      {DenseVector{1.0, 1.0}},
      trusted_metric);
    const NonDiagonalMetric spoofed_metric(layout);

    require(bounds.supports_projection_in(trusted_metric),
            "Cellwise box rejected its diagonal L2 metric");
    require(!bounds.supports_projection_in(spoofed_metric),
            "A non-diagonal metric obtained clipping projection by reusing the l2_cellwise display identifier");
  }
} // namespace

int
main(const int argc, char **argv)
{
  try
    {
      const std::vector<nmopt::test_support::Scenario> scenarios{
        {"layout_invariant",
         "nmopt.contract.layout_invariant",
         {"backend-neutral", "contract"},
         30,
         test_block_layout_invariant},
        {"v0_contract",
         "nmopt.contract.v0",
         {"backend-neutral", "contract"},
         30,
         test_v0_contract},
        {"backend_parameterisation",
         "nmopt.contract.backend_parameterisation",
         {"backend-neutral", "contract"},
         30,
         test_backend_parameterisation},
        {"owned_reduced_service_lifetime",
         "nmopt.contract.owned_reduced_service_lifetime",
         {"backend-neutral", "contract", "ownership"},
         30,
         test_owned_reduced_service_lifetime},
        {"projection_compatibility",
         "nmopt.contract.projection_compatibility",
         {"backend-neutral", "contract", "constraint"},
         30,
         test_projection_compatibility}};
      const auto result =
        nmopt::test_support::run_requested_scenarios(
          argc, argv, scenarios, std::cout);
      if (!result.listed)
        std::cout << "nmopt executable contract scenario passed: "
                  << result.executed << '\n';
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "nmopt v0 executable contract test failed: "
                << exception.what() << '\n';
      return 1;
    }
}
