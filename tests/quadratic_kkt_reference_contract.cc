#include "nmopt/contract/reduced_dto.hpp"
#include "nmopt/reference/linear_quadratic_model.hpp"
#include "nmopt/reference/quadratic_kkt.hpp"
#include "nmopt/reference/supplied_linear_quadratic_system.hpp"
#include "test_support/contract_errors.hpp"
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
      0.4,
      nmopt::reference::make_reference_supplied_otd_declaration());
  }

  SuppliedOTDSystem
  make_validity_test_system(
    SuppliedOTDQuadraticKKTValidity validity =
      make_canonical_supplied_otd_quadratic_kkt_validity(),
    const bool point_dependent_jvp = false)
  {
    const auto variable_layout = std::make_shared<const BlockLayout>(
      "validity_test_variables",
      std::vector<SpaceId>{{"validity_state"},
                           {"validity_adjoint"},
                           {"validity_control"}},
      std::vector<std::size_t>{1, 1, 1});
    const auto residual_layout = std::make_shared<const BlockLayout>(
      "validity_test_residuals",
      std::vector<SpaceId>{{"validity_state_equation"},
                           {"validity_adjoint_equation"},
                           {"validity_control_stationarity"}},
      std::vector<std::size_t>{1, 1, 1});
    const SuppliedOTDLayout layout(variable_layout, residual_layout);

    const auto residual = [residual_layout](const PrimalBlock &point) {
      return CovectorBlock(
        residual_layout,
        {DenseVector{point.block(0)[0]},
         DenseVector{point.block(1)[0]},
         DenseVector{point.block(2)[0]}});
    };
    const auto residual_jvp =
      [residual_layout, point_dependent_jvp](const PrimalBlock &point,
                                             const PrimalBlock &tangent) {
        const double factor =
          point_dependent_jvp ? 1.0 + point.block(0)[0] : 1.0;
        return CovectorBlock(
          residual_layout,
          {DenseVector{factor * tangent.block(0)[0]},
           DenseVector{factor * tangent.block(1)[0]},
           DenseVector{factor * tangent.block(2)[0]}});
      };
    const auto residual_vjp = [variable_layout](const PrimalBlock &,
                                                 const PrimalBlock &seed) {
      return CovectorBlock(
        variable_layout,
        {seed.block(0), seed.block(1), seed.block(2)});
    };
    const auto solve = [variable_layout](const PrimalBlock &) {
      return SuppliedOTDSystem::SolveResult(
        PrimalBlock::zeros(variable_layout));
    };

    return SuppliedOTDSystem(layout,
                             residual,
                             residual_jvp,
                             residual_vjp,
                             solve,
                             std::move(validity));
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
    const Product supplied_product = nmopt::contract::
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
    const Product product = nmopt::contract::
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

  void
  test_supplied_otd_validity_diagnostics()
  {
    const auto require_adapter_error =
      [](const SuppliedOTDSystem &system,
         const std::string &        message,
         const std::string &        description) {
        nmopt::test_support::require_contract_error(
          [&] {
            (void)nmopt::contract::make_canonical_supplied_otd_kkt_product(
              system);
          },
          message,
          description);
      };

    require_adapter_error(
      make_validity_test_system(SuppliedOTDQuadraticKKTValidity{}),
      "Canonical supplied OTD KKT adapter needs a declared block selection",
      "supplied adapter accepted an undeclared validity block");

    auto point_dependent = make_canonical_supplied_otd_quadratic_kkt_validity();
    point_dependent.constant_jvp_declared = false;
    require_adapter_error(
      make_validity_test_system(std::move(point_dependent), true),
      "Canonical supplied OTD KKT adapter needs a constant-JVP declaration",
      "supplied adapter accepted a point-dependent JVP without a declaration");

    auto bad_signs = make_canonical_supplied_otd_quadratic_kkt_validity();
    bad_signs.block_signs =
      SuppliedOTDQuadraticKKTBlockSigns::incompatible;
    require_adapter_error(
      make_validity_test_system(std::move(bad_signs)),
      "Canonical supplied OTD KKT adapter needs canonical block-sign declarations",
      "supplied adapter accepted incompatible canonical block signs");

    auto missing_symmetry =
      make_canonical_supplied_otd_quadratic_kkt_validity();
    missing_symmetry.symmetry_declared = false;
    require_adapter_error(
      make_validity_test_system(std::move(missing_symmetry)),
      "Canonical supplied OTD KKT adapter needs symmetry evidence",
      "supplied adapter accepted undeclared symmetry evidence");

    auto missing_rank = make_canonical_supplied_otd_quadratic_kkt_validity();
    missing_rank.rank_condition_declared = false;
    require_adapter_error(
      make_validity_test_system(std::move(missing_rank)),
      "Canonical supplied OTD KKT adapter needs a rank declaration",
      "supplied adapter accepted an undeclared rank condition");

    auto missing_kernel =
      make_canonical_supplied_otd_quadratic_kkt_validity();
    missing_kernel.kernel_positivity_declared = false;
    require_adapter_error(
      make_validity_test_system(std::move(missing_kernel)),
      "Canonical supplied OTD KKT adapter needs a kernel-positivity declaration",
      "supplied adapter accepted an undeclared kernel condition");

    auto missing_conversion =
      make_canonical_supplied_otd_quadratic_kkt_validity();
    missing_conversion.multiplier_conversion_kind =
      SuppliedOTDQuadraticKKTMultiplierConversion::incompatible;
    require_adapter_error(
      make_validity_test_system(std::move(missing_conversion)),
      "Canonical supplied OTD KKT adapter needs a compatible multiplier-conversion declaration",
      "supplied adapter accepted an incompatible multiplier conversion");

    auto nonsymmetric =
      make_canonical_supplied_otd_quadratic_kkt_validity();
    nonsymmetric.symmetry = QuadraticKKTSymmetry::nonsymmetric;
    nonsymmetric.symmetry_policy =
      "declared nonsymmetric supplied KKT actions require GMRES";
    const Product nonsymmetric_product =
      nmopt::contract::make_canonical_supplied_otd_kkt_product(
        make_validity_test_system(std::move(nonsymmetric)));
    require(nonsymmetric_product.symmetry() ==
              QuadraticKKTSymmetry::nonsymmetric,
            "declared nonsymmetric supplied adapter changed its symmetry");
    require(!nonsymmetric_product.supports_minres(),
            "declared nonsymmetric supplied adapter accepted MINRES");
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
         test_supplied_adapter_matches_reduced_dto_solution},
        {"supplied_otd_validity_diagnostics",
         "nmopt.quadratic_kkt_reference.supplied_otd_validity_diagnostics",
         {"backend-neutral", "formulation", "kkt", "reference"},
         30,
         test_supplied_otd_validity_diagnostics}};
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
