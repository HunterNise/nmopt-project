#include "nmopt/contract/metric_constraint.hpp"
#include "nmopt/contract/reduced_dto.hpp"
#include "nmopt/reference/linear_quadratic_model.hpp"

#include <cmath>
#include <exception>
#include <iostream>
#include <memory>
#include <string>
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
      value.block(block).add_scaled(step, direction.block(block));
    return value;
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
        finite_difference.block(block).add_scaled(
          -1.0, residual_at_point.block(block));
        finite_difference.block(block).scale(1.0 / epsilon);
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

    const StateControlPartition partition(model, 0, 1);
    const StateAdjointSolvers solvers{
      [&model](const PrimalBlock &control) { return model.solve_state(control); },
      [&model](const PrimalBlock &full_point, const CovectorBlock &state_rhs) {
        return model.solve_adjoint(full_point, state_rhs);
      }};
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

    const DiagonalMetric metric(
      "l2_cellwise", partition.control_layout(), {DenseVector{2.0, 5.0}});
    const PrimalBlock direction =
      reduced.gradient_direction(evaluation.reduced_derivative, metric);
    const CovectorBlock reconstructed = metric.apply(direction);
    require_close(pair(reconstructed, control_direction),
                  pair(evaluation.reduced_derivative, control_direction),
                  1e-13,
                  "Metric inverse/apply consistency");

    const CellwiseBoxConstraint bounds(
      partition.control_layout(), {DenseVector{-0.1, -0.2}},
      {DenseVector{0.6, 0.5}});
    const PrimalBlock projected = bounds.project_in(
      PrimalBlock(partition.control_layout(), {DenseVector{-0.4, 0.8}}), metric);
    require(bounds.is_feasible(projected),
            "Cellwise L2 projection must return a feasible control");
    require_close(projected.block(0)[0], -0.1, 1e-15,
                  "Cellwise box lower projection");
    require_close(projected.block(0)[1], 0.5, 1e-15,
                  "Cellwise box upper projection");
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
} // namespace

int
main()
{
  try
    {
      test_v0_contract();
      test_backend_parameterisation();
      std::cout << "nmopt v0 executable contract tests passed\n";
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "nmopt v0 executable contract test failed: "
                << exception.what() << '\n';
      return 1;
    }
}
