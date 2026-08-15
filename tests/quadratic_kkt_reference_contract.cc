#include "nmopt/contract/reduced_dto.hpp"
#include "nmopt/reference/linear_quadratic_model.hpp"
#include "nmopt/reference/quadratic_kkt.hpp"
#include "nmopt/reference/supplied_linear_quadratic_system.hpp"
#include "test_support/scenario_dispatch.hpp"

#include <cmath>
#include <iostream>
#include <vector>

namespace
{
  using namespace nmopt::contract;
  using namespace nmopt::reference;
  using Product = EqualityConstrainedQuadraticKKTProduct;

  void
  require_close(const double actual,
                const double expected,
                const char *message,
                const double tolerance = 1e-12)
  {
    require(std::abs(actual - expected) <= tolerance, message);
  }

  void
  require_vector_close(const DenseVector &actual,
                       const DenseVector &expected,
                       const char *message,
                       const double tolerance = 1e-12)
  {
    require(actual.size() == expected.size(), message);
    for (std::size_t index = 0; index < actual.size(); ++index)
      require_close(actual[index], expected[index], message, tolerance);
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
      0.4);
  }

  Product::Point
  make_point(const Product &product,
             const DenseVector &state,
             const DenseVector &control,
             const DenseVector &multiplier)
  {
    return Product::Point{
      Product::Primal(product.layout().primal, {state, control}),
      Product::Primal(product.layout().multiplier, {multiplier})};
  }

  void
  require_same_residual(const Product::Residual &left,
                        const Product::Residual &right,
                        const char *             message)
  {
    for (std::size_t block = 0; block < left.stationarity.n_blocks(); ++block)
      require_vector_close(left.stationarity.block(block),
                           right.stationarity.block(block),
                           message);
    for (std::size_t block = 0; block < left.equality.n_blocks(); ++block)
      require_vector_close(left.equality.block(block),
                           right.equality.block(block),
                           message);
  }

  void
  test_dto_and_supplied_adapters_agree()
  {
    const LinearQuadraticModel dto = make_dto_target();
    const SuppliedLinearQuadraticSystem supplied = make_supplied_target();
    const Product dto_product = make_dto_kkt_product(dto.kkt_data());
    const Product supplied_product =
      make_canonical_supplied_otd_kkt_product(supplied.system());

    const DenseVector state{0.2, -0.3};
    const DenseVector control{0.4, -0.1};
    const DenseVector multiplier{0.5, -0.25};
    const Product::Point dto_point =
      make_point(dto_product, state, control, multiplier);
    const Product::Point supplied_point =
      make_point(supplied_product, state, control, multiplier);

    const auto dto_action = dto_product.apply_kkt(dto_point);
    const auto supplied_action = supplied_product.apply_kkt(supplied_point);
    require_same_residual(dto_action,
                          supplied_action,
                          "DTO and supplied-OTD KKT actions disagree");
    require_same_residual(dto_product.residual(dto_point),
                          supplied_product.residual(supplied_point),
                          "DTO and supplied-OTD KKT residuals disagree");

    const auto dto_adjoint = dto_product.multiplier_to_adjoint(
      dto_point.multiplier);
    const auto supplied_adjoint = supplied_product.multiplier_to_adjoint(
      supplied_point.multiplier);
    require_vector_close(dto_adjoint.block(0),
                         supplied_adjoint.block(0),
                         "DTO and supplied-OTD multiplier conversions disagree");
  }

  void
  test_supplied_adapter_matches_reduced_dto_solution()
  {
    const LinearQuadraticModel dto = make_dto_target();
    const SuppliedLinearQuadraticSystem supplied = make_supplied_target();
    const Product product =
      make_canonical_supplied_otd_kkt_product(supplied.system());

    const auto supplied_result = supplied.system().solve(
      PrimalBlock::zeros(supplied.system().variable_layout()));
    require(supplied_result.report.converged(),
            "reference supplied-OTD solve did not converge");
    const auto &selection = supplied.system().block_selection();
    const DenseVector state =
      supplied_result.solution.block(selection.state_variable);
    const DenseVector adjoint =
      supplied_result.solution.block(selection.adjoint_variable);
    const DenseVector control =
      supplied_result.solution.block(selection.control_variable);
    DenseVector multiplier = adjoint;
    multiplier.scale(-1.0);

    const Product::Point point =
      make_point(product, state, control, multiplier);
    const auto residual = product.residual(point);
    require_close(dot(residual.stationarity.block(0),
                       residual.stationarity.block(0)),
                  0.0,
                  "supplied adapter solution has a state stationarity residual",
                  1e-22);
    require_close(dot(residual.stationarity.block(1),
                       residual.stationarity.block(1)),
                  0.0,
                  "supplied adapter solution has a control stationarity residual",
                  1e-22);
    require_close(dot(residual.equality.block(0), residual.equality.block(0)),
                  0.0,
                  "supplied adapter solution has an equality residual",
                  1e-22);

    const StateControlPartition partition(dto, 0, 1);
    const StateAdjointSolvers solvers{
      [&dto](const PrimalBlock &control_value) {
        return FormulationSolveResultT<DenseBackend>(
          dto.solve_state(control_value));
      },
      [&dto](const PrimalBlock &full_point,
             const CovectorBlock &state_rhs) {
        return FormulationSolveResultT<DenseBackend>(
          dto.solve_adjoint(full_point, state_rhs));
      }};
    const ReducedDTO reduced(dto, partition, solvers);
    const PrimalBlock reduced_control(partition.control_layout(), {control});
    const ReducedEvaluation evaluation = reduced.evaluate(reduced_control);

    require_vector_close(state,
                         evaluation.state.block(0),
                         "KKT supplied-OTD state disagrees with reduced DTO");
    require_vector_close(adjoint,
                         evaluation.adjoint.block(0),
                         "KKT supplied-OTD adjoint disagrees with reduced DTO");
    require_vector_close(
      product.residual(point).stationarity.block(1),
      evaluation.reduced_derivative.block(0),
      "KKT stationarity disagrees with reduced DTO derivative");
  }
} // namespace

int
main(const int argc, char **argv)
{
  try
    {
      const std::vector<nmopt::test_support::Scenario> scenarios{
        {"dto_and_supplied_adapters_agree",
         "nmopt.quadratic_kkt_reference.dto_and_supplied_adapters_agree",
         {"backend-neutral", "formulation", "kkt", "reference"},
         30,
         test_dto_and_supplied_adapters_agree},
        {"supplied_adapter_matches_reduced_dto_solution",
         "nmopt.quadratic_kkt_reference.supplied_adapter_matches_reduced_dto_solution",
         {"backend-neutral", "formulation", "kkt", "reference"},
         30,
         test_supplied_adapter_matches_reduced_dto_solution}};
      const auto result = nmopt::test_support::run_requested_scenarios(
        argc, argv, scenarios, std::cout);
      if (!result.listed)
        std::cout << "quadratic KKT reference scenario passed: "
                  << result.executed << '\n';
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "quadratic KKT reference test failed: " << exception.what()
                << '\n';
      return 1;
    }
}
