#include "nmopt/contract/quadratic_kkt.hpp"
#include "nmopt/contract/quadratic_kkt_solver.hpp"
#include "../support/contract_errors.hpp"
#include "../support/scenario_dispatch.hpp"

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace
{
  using namespace nmopt::contract;
  using Product = EqualityConstrainedQuadraticKKTProduct;
  using Point = Product::Point;
  using Seed = Product::Seed;

  struct TestTarget
  {
    LayoutPtr primal = std::make_shared<const BlockLayout>(
      "kkt_primal", std::vector<SpaceId>{{"primal"}},
      std::vector<std::size_t>{3});
    LayoutPtr multiplier = std::make_shared<const BlockLayout>(
      "kkt_multiplier", std::vector<SpaceId>{{"multiplier"}},
      std::vector<std::size_t>{2});
    LayoutPtr adjoint = std::make_shared<const BlockLayout>(
      "kkt_adjoint", std::vector<SpaceId>{{"adjoint"}},
      std::vector<std::size_t>{2});
    LayoutPtr stationarity = std::make_shared<const BlockLayout>(
      "kkt_stationarity", std::vector<SpaceId>{{"stationarity"}},
      std::vector<std::size_t>{3});
    LayoutPtr equality = std::make_shared<const BlockLayout>(
      "kkt_equality", std::vector<SpaceId>{{"equality"}},
      std::vector<std::size_t>{2});

    DenseMatrix q{3, 3, {4.0, 1.0, 0.0, 1.0, 3.0, 0.5, 0.0, 0.5, 2.0}};
    DenseMatrix d{2, 3, {1.0, -1.0, 0.0, 0.0, 1.0, 2.0}};
  };

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
                       const char *message)
  {
    require(actual.size() == expected.size(), message);
    for (std::size_t index = 0; index < actual.size(); ++index)
      require_close(actual[index], expected[index], message);
  }

  Product
  make_product(const QuadraticKKTSymmetry symmetry =
                 QuadraticKKTSymmetry::symmetric_indefinite,
               const bool rank_declared = true,
               const bool kernel_declared = true,
               const bool d_transpose_declared = true,
               const bool kkt_transpose_declared = true)
  {
    const auto target = std::make_shared<const TestTarget>();
    const Product::Layout layout(target->primal,
                                 target->multiplier,
                                 target->adjoint,
                                 target->stationarity,
                                 target->equality,
                                 {"test_primal_stationarity",
                                  {0},
                                  {0},
                                  {"test_primal_pairing"}},
                                 {"test_multiplier_equality",
                                  {0},
                                  {0},
                                  {"test_multiplier_pairing"}});

    const auto quadratic_action = [target](const Product::Primal &primal) {
      return Product::Covector(
        target->stationarity,
        {target->q.vmult(primal.block(0))});
    };
    const auto equality_action = [target](const Product::Primal &primal) {
      return Product::Covector(target->equality,
                               {target->d.vmult(primal.block(0))});
    };
    const auto multiplier_action = [target](const Product::Primal &multiplier) {
      return Product::Covector(
        target->stationarity,
        {target->d.transpose_vmult(multiplier.block(0))});
    };
    const auto transpose_action = [target](const Seed &seed) {
      DenseVector primal = target->q.transpose_vmult(seed.stationarity.block(0));
      primal.add_scaled(1.0,
                        target->d.transpose_vmult(seed.equality.block(0)));
      DenseVector multiplier = target->d.vmult(seed.stationarity.block(0));
      return Product::TransposeResult{
        Product::Covector(target->primal, {std::move(primal)}),
        Product::Covector(target->multiplier, {std::move(multiplier)})};
    };

    const QuadraticKKTMultiplierConversion conversion{
      "lambda equals negative framework adjoint",
      [target](const Product::Primal &multiplier) {
        DenseVector adjoint = multiplier.block(0);
        adjoint.scale(-1.0);
        return Product::Primal(target->adjoint, {std::move(adjoint)});
      },
      [target](const Product::Primal &adjoint) {
        DenseVector multiplier = adjoint.block(0);
        multiplier.scale(-1.0);
        return Product::Primal(target->multiplier, {std::move(multiplier)});
      }};
    const QuadraticKKTAssumptions assumptions{
      rank_declared,
      kernel_declared,
      "declared full row rank",
      "declared positive on ker(D)",
      d_transpose_declared,
      kkt_transpose_declared,
      "test D-transpose and KKT-transpose actions are declared exact"};

    return Product(layout,
                   quadratic_action,
                   equality_action,
                   multiplier_action,
                   transpose_action,
                   Product::Covector(target->stationarity,
                                     {DenseVector{1.0, 2.0, 3.0}}),
                   Product::Covector(target->equality,
                                     {DenseVector{4.0, 5.0}}),
                   conversion,
                   assumptions,
                   symmetry);
  }

  void
  test_kkt_actions_residual_and_transpose()
  {
    const Product product = make_product();
    const Point point{
      Product::Primal(product.layout().primal,
                      {DenseVector{0.5, -0.2, 1.0}}),
      Product::Primal(product.layout().multiplier,
                      {DenseVector{0.3, -0.4}})};

    const auto action = product.apply_kkt(point);
    require_vector_close(action.stationarity.block(0),
                         DenseVector{2.1, -0.3, 1.1},
                         "KKT block action has the wrong stationarity value");
    require_vector_close(action.equality.block(0),
                         DenseVector{0.7, 1.8},
                         "KKT block action has the wrong equality value");

    const auto residual = product.residual(point);
    require_vector_close(residual.stationarity.block(0),
                         DenseVector{1.1, -2.3, -1.9},
                         "KKT residual has the wrong stationarity value");
    require_vector_close(residual.equality.block(0),
                         DenseVector{-3.3, -3.2},
                         "KKT residual has the wrong equality value");

    const Seed seed{
      Product::Primal(product.layout().stationarity,
                      {DenseVector{0.7, -0.1, 0.4}}),
      Product::Primal(product.layout().equality,
                      {DenseVector{-0.6, 0.8}})};
    const auto transpose = product.apply_kkt_transpose(seed);
    const double left = dot(action.stationarity.block(0),
                            seed.stationarity.block(0)) +
                        dot(action.equality.block(0), seed.equality.block(0));
    const double right =
      dot(transpose.primal.block(0), point.primal.block(0)) +
      dot(transpose.multiplier.block(0), point.multiplier.block(0));
    require_close(left, right, "KKT action and transpose violate pairing");
  }

  void
  test_kkt_diagnostics_and_multiplier_conversion()
  {
    const Product symmetric = make_product();
    require(symmetric.supports_minres(),
            "symmetric-indefinite KKT product rejected MINRES compatibility");
    const Product nonsymmetric =
      make_product(QuadraticKKTSymmetry::nonsymmetric);
    require(!nonsymmetric.supports_minres(),
            "nonsymmetric KKT product accepted MINRES compatibility");

    const auto multiplier = Product::Primal(
      symmetric.layout().multiplier, {DenseVector{0.25, -0.75}});
    const auto adjoint = symmetric.multiplier_to_adjoint(multiplier);
    require_vector_close(adjoint.block(0),
                         DenseVector{-0.25, 0.75},
                         "multiplier-to-adjoint conversion has the wrong sign");
    const auto round_trip = symmetric.adjoint_to_multiplier(adjoint);
    require_vector_close(round_trip.block(0),
                         multiplier.block(0),
                         "adjoint-to-multiplier conversion is not inverse");

    nmopt::test_support::require_contract_error(
      [] { (void)make_product(QuadraticKKTSymmetry::symmetric_indefinite,
                              false,
                              true); },
      "Quadratic KKT product needs a declared rank condition",
      "KKT product accepted an undeclared rank condition");
    nmopt::test_support::require_contract_error(
      [] { (void)make_product(QuadraticKKTSymmetry::symmetric_indefinite,
                              true,
                              false); },
      "Quadratic KKT product needs a declared kernel-positivity condition",
      "KKT product accepted an undeclared kernel condition");

    const auto wrong_layout = std::make_shared<const BlockLayout>(
      "wrong", std::vector<SpaceId>{{"wrong"}}, std::vector<std::size_t>{3});
    const Product::Primal wrong_primal(wrong_layout, {DenseVector(3)});
    nmopt::test_support::require_contract_error(
      [&] { (void)symmetric.apply_q(wrong_primal); },
      "Quadratic KKT Q input has an incompatible layout",
      "KKT Q action accepted an incompatible primal layout");

    const auto target = std::make_shared<const TestTarget>();
    nmopt::test_support::require_contract_error(
      [&] {
        (void)Product::Layout(target->primal,
                              target->multiplier,
                              target->adjoint,
                              target->stationarity,
                              target->equality);
      },
      "Quadratic KKT primal/stationarity pairing needs an identifier",
      "KKT layout accepted a missing primal/stationarity pairing");

    const auto incompatible_stationarity =
      std::make_shared<const BlockLayout>(
        "incompatible_stationarity",
        std::vector<SpaceId>{{"stationarity"}},
        std::vector<std::size_t>{4});
    nmopt::test_support::require_contract_error(
      [&] {
        (void)Product::Layout(
          target->primal,
          target->multiplier,
          target->adjoint,
          incompatible_stationarity,
          target->equality,
          {"bad_primal_stationarity", {0}, {0}, {"bad_pairing"}},
          {"test_multiplier_equality", {0}, {0}, {"test_pairing"}});
      },
      "Quadratic KKT primal/stationarity pairing has incompatible block dimensions",
      "KKT layout accepted incompatible primal/stationarity dimensions");

    const auto incompatible_equality = std::make_shared<const BlockLayout>(
      "incompatible_equality",
      std::vector<SpaceId>{{"equality"}},
      std::vector<std::size_t>{3});
    nmopt::test_support::require_contract_error(
      [&] {
        (void)Product::Layout(
          target->primal,
          target->multiplier,
          target->adjoint,
          target->stationarity,
          incompatible_equality,
          {"test_primal_stationarity", {0}, {0}, {"test_pairing"}},
          {"bad_multiplier_equality", {0}, {0}, {"bad_pairing"}});
      },
      "Quadratic KKT multiplier/equality pairing has incompatible block dimensions",
      "KKT layout accepted incompatible multiplier/equality dimensions");

    nmopt::test_support::require_contract_error(
      [] {
        (void)make_product(QuadraticKKTSymmetry::symmetric_indefinite,
                           true,
                           true,
                           false,
                           true);
      },
      "Symmetric quadratic KKT product needs declared D-transpose consistency",
      "KKT product accepted undeclared D-transpose consistency");
    nmopt::test_support::require_contract_error(
      [] {
        (void)make_product(QuadraticKKTSymmetry::symmetric_indefinite,
                           true,
                           true,
                           true,
                           false);
      },
      "Symmetric quadratic KKT product needs declared KKT-transpose consistency",
      "KKT product accepted undeclared KKT-transpose consistency");
  }

  void
  test_kkt_solver_policy_compatibility()
  {
    const Product symmetric = make_product();
    QuadraticKKTSolverPolicy minres_policy;
    validate(symmetric, minres_policy);
    minres_policy.gmres_maximum_basis = 2;
    require(valid(minres_policy),
            "MINRES validity depends on its unused GMRES basis");
    validate(symmetric, minres_policy);

    QuadraticKKTSolverPolicy gmres_policy;
    gmres_policy.method = QuadraticKKTSolverMethod::gmres;
    require(valid(gmres_policy),
            "valid GMRES policy was rejected by policy validity");
    validate(symmetric, gmres_policy);

    const Product nonsymmetric =
      make_product(QuadraticKKTSymmetry::nonsymmetric);
    validate(nonsymmetric, gmres_policy);
    nmopt::test_support::require_contract_error(
      [&] { validate(nonsymmetric, minres_policy); },
      "MINRES requires a symmetric-indefinite KKT product",
      "KKT policy accepted MINRES for a nonsymmetric product");

    QuadraticKKTSolverPolicy invalid_policy;
    invalid_policy.relative_tolerance = 0.0;
    nmopt::test_support::require_contract_error(
      [&] { validate(symmetric, invalid_policy); },
      "Quadratic KKT solver policy needs positive finite tolerances",
      "KKT policy accepted a non-positive tolerance");

    invalid_policy = {};
    invalid_policy.method = QuadraticKKTSolverMethod::gmres;
    invalid_policy.gmres_maximum_basis = 2;
    nmopt::test_support::require_contract_error(
      [&] { validate(symmetric, invalid_policy); },
      "Quadratic KKT solver policy needs a GMRES basis of at least three vectors",
      "KKT policy accepted an undersized GMRES basis");
  }

  void
  test_kkt_solver_report_conjunction()
  {
    QuadraticKKTSolveReport linear_converged_residual_failed;
    linear_converged_residual_failed.linear_solve.termination =
      LinearSolveTermination::converged;
    linear_converged_residual_failed.residuals_converged = false;
    require(!linear_converged_residual_failed.converged(),
            "KKT report accepted a failed residual after linear convergence");

    QuadraticKKTSolveReport linear_failed_residual_converged;
    linear_failed_residual_converged.linear_solve.termination =
      LinearSolveTermination::failed;
    linear_failed_residual_converged.residuals_converged = true;
    require(!linear_failed_residual_converged.converged(),
            "KKT report accepted a failed linear solve");
  }
} // namespace

int
main(const int argc, char **argv)
{
  try
    {
      const std::vector<nmopt::test_support::Scenario> scenarios{
        {"actions_residual_and_transpose",
         "nmopt.quadratic_kkt.actions_residual_and_transpose",
         {"backend-neutral", "formulation", "kkt"},
         30,
         test_kkt_actions_residual_and_transpose},
        {"diagnostics_and_multiplier_conversion",
         "nmopt.quadratic_kkt.diagnostics_and_multiplier_conversion",
         {"backend-neutral", "formulation", "kkt"},
         30,
         test_kkt_diagnostics_and_multiplier_conversion},
        {"solver_policy_compatibility",
         "nmopt.quadratic_kkt.solver_policy_compatibility",
         {"backend-neutral", "formulation", "kkt", "solver"},
         30,
         test_kkt_solver_policy_compatibility},
        {"solver_report_conjunction",
         "nmopt.quadratic_kkt.solver_report_conjunction",
         {"backend-neutral", "formulation", "kkt", "solver"},
         30,
         test_kkt_solver_report_conjunction}};
      const auto result = nmopt::test_support::run_requested_scenarios(
        argc, argv, scenarios, std::cout);
      if (!result.listed)
        std::cout << "quadratic KKT contract scenario passed: "
                  << result.executed << '\n';
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "quadratic KKT contract test failed: " << exception.what()
                << '\n';
      return 1;
    }
}
