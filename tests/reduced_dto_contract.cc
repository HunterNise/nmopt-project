#include "nmopt/contract/metric_constraint.hpp"
#include "nmopt/contract/reduced_dto.hpp"
#include "nmopt/reference/linear_quadratic_model.hpp"
#include "nmopt/solvers/reduced_gradient.hpp"
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
    for (std::size_t index = 0;
         index < solver_result.accepted_iterations;
         ++index)
      {
        require(solver_result.step_length_history[index] > 0.0,
                "Dense reduced gradient accepted a nonpositive step");
        require_close(solver_result.objective_change_history[index],
                      solver_result.objective_history[index + 1] -
                        solver_result.objective_history[index],
                      1e-14,
                      "Dense reduced gradient objective-change history");
      }

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
