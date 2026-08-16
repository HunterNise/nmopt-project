#include "nmopt/contract/pdas.hpp"
#include "nmopt/contract/quadratic_kkt_solver.hpp"
#include "test_support/contract_errors.hpp"
#include "test_support/scenario_dispatch.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

namespace
{
  using namespace nmopt::contract;
  using Product = EqualityConstrainedQuadraticKKTProduct;
  using Point = Product::Point;
  using Covector = Product::Covector;
  using DenseVector = nmopt::contract::DenseVector;

  struct Target
  {
    LayoutPtr primal = std::make_shared<const BlockLayout>(
      "primal", std::vector<SpaceId>{{"state"}, {"control"}},
      std::vector<std::size_t>{2, 2});
    LayoutPtr multiplier = std::make_shared<const BlockLayout>(
      "multiplier", std::vector<SpaceId>{{"equality_multiplier"}},
      std::vector<std::size_t>{2});
    LayoutPtr adjoint = std::make_shared<const BlockLayout>(
      "adjoint", std::vector<SpaceId>{{"adjoint"}},
      std::vector<std::size_t>{2});
    LayoutPtr stationarity = std::make_shared<const BlockLayout>(
      "stationarity", std::vector<SpaceId>{{"state"}, {"control"}},
      std::vector<std::size_t>{2, 2});
    LayoutPtr equality = std::make_shared<const BlockLayout>(
      "equality", std::vector<SpaceId>{{"state_equation"}},
      std::vector<std::size_t>{2});
  };

  void
  require_close(const double actual,
                const double expected,
                const char *message,
                const double tolerance = 1e-11)
  {
    require(std::abs(actual - expected) <= tolerance, message);
  }

  Product
  make_product()
  {
    const auto target = std::make_shared<const Target>();
    const Product::Layout layout(target->primal,
                                 target->multiplier,
                                 target->adjoint,
                                 target->stationarity,
                                 target->equality);

    const auto quadratic_action = [target](const Product::Primal &primal) {
      return Product::Covector(
        target->stationarity,
        {primal.block(0), primal.block(1)});
    };
    const auto equality_action = [target](const Product::Primal &primal) {
      DenseVector value = primal.block(0);
      value.add_scaled(-1.0, primal.block(1));
      return Product::Covector(target->equality, {std::move(value)});
    };
    const auto multiplier_action = [target](const Product::Primal &multiplier) {
      DenseVector state = multiplier.block(0);
      DenseVector control = multiplier.block(0);
      control.scale(-1.0);
      return Product::Covector(target->stationarity,
                               {std::move(state), std::move(control)});
    };
    const auto transpose_action = [target](const Product::Seed &seed) {
      DenseVector state = seed.stationarity.block(0);
      state.add_scaled(1.0, seed.equality.block(0));
      DenseVector control = seed.stationarity.block(1);
      control.add_scaled(-1.0, seed.equality.block(0));
      DenseVector multiplier = seed.stationarity.block(0);
      multiplier.add_scaled(-1.0, seed.stationarity.block(1));
      return Product::TransposeResult{
        Product::Covector(target->primal,
                          {std::move(state), std::move(control)}),
        Product::Covector(target->multiplier, {std::move(multiplier)})};
    };
    const Product::MultiplierConversion conversion{
      "lambda equals negative framework adjoint",
      [target](const Product::Primal &multiplier) {
        DenseVector value = multiplier.block(0);
        value.scale(-1.0);
        return Product::Primal(target->adjoint, {std::move(value)});
      },
      [target](const Product::Primal &adjoint) {
        DenseVector value = adjoint.block(0);
        value.scale(-1.0);
        return Product::Primal(target->multiplier, {std::move(value)});
      }};
    const QuadraticKKTAssumptions assumptions{
      true, true, "full row rank", "positive on the equality kernel"};

    return Product(
      layout,
      quadratic_action,
      equality_action,
      multiplier_action,
      transpose_action,
      Product::Covector(target->stationarity,
                        {DenseVector{0.0, 0.0}, DenseVector{2.0, 0.25}}),
      Product::Covector(target->equality,
                        {DenseVector{0.0, 0.0}}),
      conversion,
      assumptions,
      QuadraticKKTSymmetry::symmetric_indefinite);
  }

  std::size_t
  layout_dimension(const LayoutPtr &layout)
  {
    std::size_t dimension = 0;
    for (std::size_t block = 0; block < layout->n_blocks(); ++block)
      dimension += layout->dimension(block);
    return dimension;
  }

  std::vector<DenseVector>
  split_vector(const DenseVector &value,
               const LayoutPtr &layout,
               std::size_t &offset)
  {
    std::vector<DenseVector> blocks;
    blocks.reserve(layout->n_blocks());
    for (std::size_t block = 0; block < layout->n_blocks(); ++block)
      {
        DenseVector current(layout->dimension(block));
        for (std::size_t index = 0; index < layout->dimension(block); ++index)
          current[index] = value[offset + index];
        offset += layout->dimension(block);
        blocks.push_back(std::move(current));
      }
    return blocks;
  }

  DenseVector
  flatten(const Product::Residual &residual)
  {
    DenseVector result(layout_dimension(residual.stationarity.layout()) +
                       layout_dimension(residual.equality.layout()));
    std::size_t offset = 0;
    for (std::size_t block = 0;
         block < residual.stationarity.n_blocks();
         ++block)
      for (std::size_t index = 0;
           index < residual.stationarity.layout()->dimension(block);
           ++index)
        result[offset++] = residual.stationarity.block(block)[index];
    for (std::size_t block = 0;
         block < residual.equality.n_blocks();
         ++block)
      for (std::size_t index = 0;
           index < residual.equality.layout()->dimension(block);
           ++index)
        result[offset++] = residual.equality.block(block)[index];
    return result;
  }

  DenseVector
  flatten(const Product::Point &point,
          const Product &       product)
  {
    const auto action = product.apply_kkt(point);
    return flatten(action);
  }

  template <typename Values>
  double
  block_norm(const Values &values)
  {
    double squared_norm = 0.0;
    for (std::size_t block = 0; block < values.n_blocks(); ++block)
      squared_norm += DenseBackend::dot(values.block(block), values.block(block));
    return std::sqrt(squared_norm);
  }

  Point
  unflatten(const DenseVector &value, const Product &product)
  {
    const std::size_t primal_dimension =
      layout_dimension(product.layout().primal);
    std::size_t offset = 0;
    auto primal_blocks = split_vector(value, product.layout().primal, offset);
    auto multiplier_blocks =
      split_vector(value, product.layout().multiplier, offset);
    require(offset == primal_dimension +
              layout_dimension(product.layout().multiplier),
            "dense KKT solution unpacking consumed the wrong dimension");
    return Point{
      Product::Primal(product.layout().primal, std::move(primal_blocks)),
      Product::Primal(product.layout().multiplier,
                      std::move(multiplier_blocks))};
  }

  QuadraticKKTSolveResult
  solve_dense(const Product &product)
  {
    const std::size_t primal_dimension =
      layout_dimension(product.layout().primal);
    const std::size_t multiplier_dimension =
      layout_dimension(product.layout().multiplier);
    const std::size_t dimension = primal_dimension + multiplier_dimension;
    std::vector<double> entries(dimension * dimension, 0.0);
    for (std::size_t column = 0; column < dimension; ++column)
      {
        Point basis{
          Product::Primal::zeros(product.layout().primal),
          Product::Primal::zeros(product.layout().multiplier)};
        std::size_t offset = 0;
        if (column < primal_dimension)
          {
            auto blocks = split_vector(
              DenseVector(primal_dimension, 0.0),
              product.layout().primal,
              offset);
            offset = 0;
            std::size_t remaining = column;
            for (auto &block : blocks)
              if (remaining < block.size())
                {
                  block[remaining] = 1.0;
                  remaining = std::numeric_limits<std::size_t>::max();
                }
              else
                remaining -= block.size();
            basis.primal = Product::Primal(product.layout().primal,
                                           std::move(blocks));
          }
        else
          {
            const std::size_t multiplier_column = column - primal_dimension;
            auto blocks = split_vector(
              DenseVector(multiplier_dimension, 0.0),
              product.layout().multiplier,
              offset);
            std::size_t remaining = multiplier_column;
            for (auto &block : blocks)
              if (remaining < block.size())
                {
                  block[remaining] = 1.0;
                  remaining = std::numeric_limits<std::size_t>::max();
                }
              else
                remaining -= block.size();
            basis.multiplier = Product::Primal(
              product.layout().multiplier, std::move(blocks));
          }
        const DenseVector action = flatten(basis, product);
        for (std::size_t row = 0; row < dimension; ++row)
          entries[row * dimension + column] = action[row];
      }

    const Point zero{
      Product::Primal::zeros(product.layout().primal),
      Product::Primal::zeros(product.layout().multiplier)};
    DenseVector rhs = flatten(product.residual(zero));
    rhs.scale(-1.0);
    const DenseVector solution =
      DenseMatrix(dimension, dimension, std::move(entries)).solve(rhs);
    const Point point = unflatten(solution, product);
    const auto residual = product.residual(point);
    const double stationarity_residual = block_norm(residual.stationarity);
    const double equality_residual = block_norm(residual.equality);
    const LinearSolveReport linear_report{
      "dense_direct",
      "none",
      1,
      1,
      0.0,
      0.0,
      0.0,
      std::max(stationarity_residual, equality_residual),
      LinearSolveTermination::converged};
    return {point,
            {linear_report,
             stationarity_residual,
             equality_residual,
             stationarity_residual <= 1e-10 && equality_residual <= 1e-10}};
  }

  struct Problem
  {
    Product product = make_product();
    LayoutPtr control_layout = std::make_shared<const BlockLayout>(
      "control", std::vector<SpaceId>{{"control"}},
      std::vector<std::size_t>{2});
    DiagonalMetric metric{
      "l2_cellwise", control_layout, {DenseVector{2.0, 4.0}}};
  };

  BoxComplementarity
  make_complementarity(const Problem &problem,
                       const double upper_first,
                       const double upper_second)
  {
    return BoxComplementarity(
      BoxBounds(
        problem.control_layout,
        PrimalBlock(problem.control_layout,
                    {DenseVector{-2.0, -2.0}}),
        PrimalBlock(problem.control_layout,
                    {DenseVector{upper_first, upper_second}})),
      make_metric_multiplier_representation(problem.metric));
  }

  QuadraticKKTAssumptions
  active_set_assumptions()
  {
    return {true,
            true,
            "base equality plus selected active control rows have declared rank",
            "objective is positive on the augmented equality kernel"};
  }

  void
  test_active_subproblem_adds_rows_and_preserves_solution()
  {
    const Problem problem;
    const auto complementarity = make_complementarity(problem, 0.5, 1.0);
    const ActiveSetSelection selection(
      problem.control_layout,
      {BoxActivity::upper, BoxActivity::inactive});
    const ActiveSetKKTSubproblem subproblem(problem.product,
                                            complementarity,
                                            selection,
                                            1,
                                            active_set_assumptions());
    const auto &active_product = subproblem.product();

    require(subproblem.has_active_constraints(),
            "active KKT subproblem did not retain active constraints");
    require(active_product.layout().equality->n_blocks() == 2,
            "active KKT subproblem did not append an equality block");
    require(active_product.layout().equality->dimension(1) == 1,
            "active KKT subproblem appended the wrong equality dimension");
    require(active_product.layout().multiplier->n_blocks() == 2,
            "active KKT subproblem did not append a multiplier block");
    const Point active_zero{
      Product::Primal::zeros(active_product.layout().primal),
      Product::Primal::zeros(active_product.layout().multiplier)};
    const auto active_zero_residual = active_product.residual(active_zero);
    require_close(active_zero_residual.equality.block(1)[0], -0.5,
                  "active KKT subproblem used the wrong bound value");

    const auto solved = solve_dense(active_product);
    require(solved.report.converged(),
            "active KKT subproblem dense solve did not converge");
    const Point base_solution = subproblem.to_base_point(solved.solution);
    require_close(base_solution.primal.block(0)[0], 0.5,
                  "active KKT subproblem returned the wrong active state");
    require_close(base_solution.primal.block(0)[1], 0.125,
                  "active KKT subproblem returned the wrong inactive state");
    require_close(base_solution.primal.block(1)[0], 0.5,
                  "active KKT subproblem returned the wrong active control");
    require_close(base_solution.primal.block(1)[1], 0.125,
                  "active KKT subproblem returned the wrong inactive control");

    const auto equality_action = active_product.apply_d(base_solution.primal);
    require_close(equality_action.block(0)[0], 0.0,
                  "active KKT subproblem changed the base equality action");
    require_close(equality_action.block(0)[1], 0.0,
                  "active KKT subproblem changed the base equality action");
    require_close(equality_action.block(1)[0], 0.5,
                  "active KKT subproblem selected the wrong control row");

    auto transpose_seed = Product::Primal::zeros(active_product.layout().stationarity);
    auto equality_seed = Product::Primal::zeros(active_product.layout().equality);
    equality_seed.add_scaled_block(1, 1.0, DenseVector{1.0});
    const auto transpose = active_product.apply_kkt_transpose(
      Product::Seed{std::move(transpose_seed), std::move(equality_seed)});
    require_close(transpose.primal.block(1)[0], 1.0,
                  "active KKT transpose omitted the active equality row");

    auto active_multiplier =
      Product::Primal::zeros(active_product.layout().multiplier);
    active_multiplier.add_scaled_block(1, 1.0, DenseVector{1.0});
    const auto multiplier_action =
      active_product.apply_d_transpose(active_multiplier);
    require_close(multiplier_action.block(1)[0], 1.0,
                  "active KKT multiplier action omitted the active row");
  }

  void
  test_inactive_subproblem_reuses_base_product()
  {
    const Problem problem;
    const auto complementarity = make_complementarity(problem, 2.0, 2.0);
    const ActiveSetSelection selection(
      problem.control_layout,
      {BoxActivity::inactive, BoxActivity::inactive});
    const ActiveSetKKTSubproblem subproblem(problem.product,
                                            complementarity,
                                            selection,
                                            1,
                                            {});
    require(!subproblem.has_active_constraints(),
            "inactive KKT subproblem reported active constraints");
    require(&subproblem.product() == &problem.product,
            "inactive KKT subproblem unnecessarily replaced the base product");
    require(subproblem.product().layout().equality->n_blocks() == 1,
            "inactive KKT subproblem changed the base equality layout");
  }

  void
  test_invalid_active_subproblem_inputs()
  {
    const Problem problem;
    const auto complementarity = make_complementarity(problem, 0.5, 1.0);
    const ActiveSetSelection selection(
      problem.control_layout,
      {BoxActivity::upper, BoxActivity::inactive});
    nmopt::test_support::require_contract_error(
      [&] {
        (void)ActiveSetKKTSubproblem(problem.product,
                                     complementarity,
                                     selection,
                                     1,
                                     QuadraticKKTAssumptions{});
      },
      "Quadratic KKT product needs a declared rank condition",
      "active KKT subproblem accepted undeclared assumptions");
  }
} // namespace

int
main(const int argc, char **argv)
{
  try
    {
      const std::vector<nmopt::test_support::Scenario> scenarios{
        {"active_subproblem_adds_rows_and_preserves_solution",
         "nmopt.active_set_kkt.active_subproblem_adds_rows_and_preserves_solution",
         {"backend-neutral", "formulation", "kkt", "active-set"},
         30,
         test_active_subproblem_adds_rows_and_preserves_solution},
        {"inactive_subproblem_reuses_base_product",
         "nmopt.active_set_kkt.inactive_subproblem_reuses_base_product",
         {"backend-neutral", "formulation", "kkt", "active-set"},
         30,
         test_inactive_subproblem_reuses_base_product},
        {"invalid_active_subproblem_inputs",
         "nmopt.active_set_kkt.invalid_active_subproblem_inputs",
         {"backend-neutral", "formulation", "kkt", "active-set", "negative"},
         30,
         test_invalid_active_subproblem_inputs}};
      const auto result = nmopt::test_support::run_requested_scenarios(
        argc, argv, scenarios, std::cout);
      if (!result.listed)
        std::cout << "backend-neutral active-set KKT contract test passed: "
                  << result.executed << '\n';
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "backend-neutral active-set KKT contract test failed: "
                << exception.what() << '\n';
      return 1;
    }
}
