#include "nmopt/contract/reduced_dto.hpp"
#include "nmopt/reference/linear_quadratic_model.hpp"
#include "nmopt/reference/supplied_linear_quadratic_system.hpp"
#include "test_support/scenario_dispatch.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

namespace
{
  using namespace nmopt::contract;
  using nmopt::reference::LinearQuadraticModel;
  using nmopt::reference::SuppliedLinearQuadraticSystem;

  void
  require_close(const double value,
                const double expected,
                const double tolerance,
                const char *message)
  {
    require(std::abs(value - expected) <= tolerance, message);
  }

  void
  require_vector_close(const DenseVector &vector,
                       const DenseVector &expected,
                       const double       tolerance,
                       const char *       message)
  {
    require(vector.size() == expected.size(), message);
    for (std::size_t index = 0; index < vector.size(); ++index)
      require_close(vector[index], expected[index], tolerance, message);
  }

  PrimalBlock
  shifted(const PrimalBlock &point,
          const PrimalBlock &tangent,
          const double       factor)
  {
    PrimalBlock result = point;
    for (std::size_t block = 0; block < result.n_blocks(); ++block)
      result.add_scaled_block(block, factor, tangent.block(block));
    return result;
  }

  SuppliedLinearQuadraticSystem
  make_supplied_target()
  {
    return SuppliedLinearQuadraticSystem(
      DenseMatrix(2, 2, {4.0, -1.0, -1.0, 3.0}),
      DenseMatrix(2, 2, {1.0, 0.5, -0.25, 2.0}),
      DenseVector{1.0, -0.5},
      DenseMatrix(2, 2, {1.0, 0.0, 0.5, 1.0}),
      DenseVector{0.25, -1.0},
      DenseVector{1.5, 0.75},
      DenseVector{2.0, 3.0},
      0.4,
      nmopt::reference::make_reference_supplied_otd_declaration());
  }

  LinearQuadraticModel
  make_dto_target()
  {
    return LinearQuadraticModel(
      DenseMatrix(2, 2, {4.0, -1.0, -1.0, 3.0}),
      DenseMatrix(2, 2, {1.0, 0.5, -0.25, 2.0}),
      DenseVector{1.0, -0.5},
      DenseMatrix(2, 2, {1.0, 0.0, 0.5, 1.0}),
      DenseVector{0.25, -1.0},
      DenseVector{1.5, 0.75},
      DenseVector{2.0, 3.0},
      0.4);
  }

  void
  test_explicit_supplied_blocks()
  {
    const SuppliedLinearQuadraticSystem target = make_supplied_target();
    const auto &declaration = target.declaration();
    require(declaration.id == "reference_linear_quadratic_supplied_otd" &&
              declaration.state_block.role ==
                nmopt::semantic::v1::SuppliedOTDBlockRole::state &&
              declaration.adjoint_block.role ==
                nmopt::semantic::v1::SuppliedOTDBlockRole::adjoint &&
              declaration.control_stationarity_block.role ==
                nmopt::semantic::v1::SuppliedOTDBlockRole::control_stationarity &&
              declaration.state_block.test_pairing_id ==
                "state_test_pairing" &&
              declaration.adjoint_block.test_pairing_id ==
                "state_test_pairing" &&
              declaration.control_stationarity_block.test_pairing_id ==
                "control_pairing" &&
              declaration.multiplier_conversion ==
                nmopt::semantic::v1::SuppliedOTDMultiplierConversion::identity &&
              declaration.comparison_status ==
                nmopt::semantic::v1::SuppliedOTDComparisonStatus::equivalent_under_declared_conversion &&
              !declaration.comparison_evidence.empty(),
            "reference supplied OTD product did not retain its typed declaration");
    const SuppliedOTDSystem &system = target.system();
    const PrimalBlock point(
      system.variable_layout(),
      {DenseVector{0.2, -0.3}, DenseVector{0.5, -0.4}, DenseVector{0.4, -0.3}});
    const PrimalBlock tangent(
      system.variable_layout(),
      {DenseVector{-0.7, 0.25}, DenseVector{0.3, 0.8}, DenseVector{0.3, 0.8}});
    const CovectorBlock value = system.residual(point);
    const CovectorBlock direction = system.residual_jvp(point, tangent);

    constexpr double epsilon = 1e-7;
    const CovectorBlock finite_difference =
      system.residual_jvp(point, tangent);
    const CovectorBlock shifted_value =
      system.residual(shifted(point, tangent, epsilon));
    for (std::size_t block = 0; block < value.n_blocks(); ++block)
      for (std::size_t entry = 0; entry < value.block(block).size(); ++entry)
        require_close((shifted_value.block(block)[entry] -
                       value.block(block)[entry]) /
                        epsilon,
                      finite_difference.block(block)[entry],
                      1e-8,
                      "supplied OTD block JVP finite difference");

    const PrimalBlock seed(
      system.residual_layout(),
      {DenseVector{0.6, -1.1}, DenseVector{-0.2, 0.7}, DenseVector{0.8, -0.4}});
    require_close(pair(direction, seed),
                  pair(system.residual_vjp(point, seed), tangent),
                  1e-12,
                  "supplied OTD block VJP pairing");

    require_vector_close(value.block(0), DenseVector{-0.15, 0.1}, 1e-12,
                         "supplied OTD state block value");
    require_vector_close(value.block(2), DenseVector{0.92, -0.91}, 1e-12,
                         "supplied OTD stationarity block value");
  }

  void
  test_supplied_solution_matches_reduced_dto()
  {
    const SuppliedLinearQuadraticSystem supplied = make_supplied_target();
    const SuppliedOTDSystem &system = supplied.system();
    const PrimalBlock initial = PrimalBlock::zeros(system.variable_layout());
    const auto supplied_result = system.solve(initial);
    require(supplied_result.report.converged(),
            "supplied OTD reference solve did not converge");
    require(supplied_result.report.algorithm == "dense_gaussian",
            "supplied OTD reference solve reported the wrong algorithm");

    const LinearQuadraticModel dto = make_dto_target();
    const StateControlPartition partition(dto, 0, 1);
    const StateAdjointSolvers solvers{
      [&dto](const PrimalBlock &control) {
        return FormulationSolveResultT<DenseBackend>(dto.solve_state(control));
      },
      [&dto](const PrimalBlock &full_point,
             const CovectorBlock &state_rhs) {
        return FormulationSolveResultT<DenseBackend>(
          dto.solve_adjoint(full_point, state_rhs));
      }};
    const ReducedDTO reduced(dto, partition, solvers);

    const PrimalBlock supplied_control(
      partition.control_layout(),
      {supplied_result.solution.block(
        system.block_selection().control_variable)});
    const ReducedEvaluation evaluation = reduced.evaluate(supplied_control);

    require_vector_close(
      supplied_result.solution.block(system.block_selection().state_variable),
      evaluation.state.block(0),
      1e-12,
      "supplied OTD state does not match reduced DTO state");
    require_vector_close(
      supplied_result.solution.block(system.block_selection().adjoint_variable),
      evaluation.adjoint.block(0),
      1e-12,
      "supplied OTD adjoint does not match reduced DTO adjoint");
    require_vector_close(
      system.control_stationarity(supplied_result.solution).block(0),
      evaluation.reduced_derivative.block(0),
      1e-12,
      "supplied OTD stationarity does not match reduced DTO derivative");

    const std::array<CovectorBlock, 3> residuals{
      system.state_residual(supplied_result.solution),
      system.adjoint_residual(supplied_result.solution),
      system.control_stationarity(supplied_result.solution)};
    for (const CovectorBlock &residual : residuals)
      require_close(dot(residual.block(0), residual.block(0)),
                    0.0,
                    1e-24,
                    "supplied OTD direct solution leaves a residual");
  }
} // namespace

int
main(const int argc, char **argv)
{
  try
    {
      const std::vector<nmopt::test_support::Scenario> scenarios{
        {"explicit_blocks",
         "nmopt.supplied_otd_reference.explicit_blocks",
         {"backend-neutral", "formulation", "supplied-otd", "reference"},
         30,
         test_explicit_supplied_blocks},
        {"matches_reduced_dto",
         "nmopt.supplied_otd_reference.matches_reduced_dto",
         {"backend-neutral", "formulation", "supplied-otd", "reference"},
         30,
         test_supplied_solution_matches_reduced_dto}};
      const auto result = nmopt::test_support::run_requested_scenarios(
        argc, argv, scenarios, std::cout);
      if (!result.listed)
        std::cout << "supplied OTD reference scenario passed: "
                  << result.executed << '\n';
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "supplied OTD reference test failed: " << exception.what()
                << '\n';
      return 1;
    }
}
