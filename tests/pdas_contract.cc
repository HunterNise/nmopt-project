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
                                 target->equality,
                                 {"pdas_primal_stationarity",
                                  {0, 1},
                                  {0, 1},
                                  {"state_stationarity",
                                   "control_stationarity"}},
                                 {"pdas_multiplier_equality",
                                  {0},
                                  {0},
                                  {"state_equation"}});

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
      true,
      true,
      "full row rank",
      "positive on the equality kernel",
      true,
      true,
      "PDAS base D-transpose and KKT-transpose actions are declared exact"};

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
      make_metric_multiplier_representation(
        std::shared_ptr<const Metric>(
          std::make_shared<const DiagonalMetric>(problem.metric))));
  }

  QuadraticKKTAssumptions
  active_set_assumptions()
  {
    return {true,
            true,
            "base equality plus selected active control rows have declared rank",
            "objective is positive on the restricted equality kernel",
            true,
            true,
            "PDAS D-transpose and KKT-transpose actions are declared exact"};
  }

  void
  test_active_subproblem_restricts_free_coordinates_and_reconstructs_solution()
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
    require(active_product.layout().primal->n_blocks() == 2,
            "active KKT subproblem did not retain the state and free control blocks");
    require(active_product.layout().primal->dimension(1) == 1,
            "active KKT subproblem retained the wrong free control dimension");
    require(active_product.layout().stationarity->dimension(1) == 1,
            "active KKT subproblem retained the wrong free stationarity dimension");
    require(active_product.layout().equality->n_blocks() == 1,
            "active KKT subproblem changed the base equality layout");
    require(active_product.layout().multiplier->n_blocks() == 1,
            "active KKT subproblem changed the base multiplier layout");
    const Point active_zero{
      Product::Primal::zeros(active_product.layout().primal),
      Product::Primal::zeros(active_product.layout().multiplier)};
    const auto active_zero_residual = active_product.residual(active_zero);
    require_close(active_zero_residual.equality.block(0)[0], -0.5,
                  "active KKT subproblem used the wrong affine equality shift");
    require_close(active_zero_residual.stationarity.block(1)[0], -0.25,
                  "active KKT subproblem used the wrong affine stationarity shift");

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

    const auto equality_action = problem.product.apply_d(base_solution.primal);
    require_close(equality_action.block(0)[0], 0.0,
                  "active KKT subproblem changed the base equality action");
    require_close(equality_action.block(0)[1], 0.0,
                  "active KKT subproblem changed the base equality action");

    const Point mixed_point{
      Product::Primal(active_product.layout().primal,
                      {DenseVector{0.2, -0.1}, DenseVector{0.7}}),
      Product::Primal(active_product.layout().multiplier,
                      {DenseVector{0.3, -0.2}})};
    const Product::Seed mixed_seed{
      Product::Primal(active_product.layout().stationarity,
                      {DenseVector{0.4, -0.5}, DenseVector{0.6}}),
      Product::Primal(active_product.layout().equality,
                      {DenseVector{0.7, -0.8}})};
    const auto mixed_action = active_product.apply_kkt(mixed_point);
    const auto mixed_transpose = active_product.apply_kkt_transpose(mixed_seed);
    const double transpose_left =
      pair(mixed_action.stationarity, mixed_seed.stationarity) +
      pair(mixed_action.equality, mixed_seed.equality);
    const double transpose_right =
      pair(mixed_transpose.primal, mixed_point.primal) +
      pair(mixed_transpose.multiplier, mixed_point.multiplier);
    require_close(transpose_left,
                  transpose_right,
                  "active KKT restricted action/transpose pairing");
    require_close(mixed_transpose.primal.block(1)[0],
                  1.4,
                  "active KKT restricted transpose used the wrong free control action");
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
    require(&subproblem.product() != &problem.product,
            "inactive KKT subproblem did not retain a value-semantic product");
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

  Point
  zero_point(const Product &product)
  {
    return Point{Product::Primal::zeros(product.layout().primal),
                 Product::Primal::zeros(product.layout().multiplier)};
  }

  Covector
  zero_box_multiplier(const LayoutPtr &layout)
  {
    return Covector(layout, {DenseVector(layout->dimension(0), 0.0)});
  }

  PDASPolicy
  valid_policy()
  {
    PDASPolicy policy;
    policy.active_set_assumptions = active_set_assumptions();
    return policy;
  }

  QuadraticKKTSolveResult
  nonconverged_solve(const Product &product)
  {
    auto result = solve_dense(product);
    result.report.residuals_converged = false;
    return result;
  }

  QuadraticKKTSolveResult
  nonfinite_solution(const Product &product, const double value)
  {
    std::vector<DenseVector> primal_blocks;
    primal_blocks.reserve(product.layout().primal->n_blocks());
    for (std::size_t block = 0;
         block < product.layout().primal->n_blocks();
         ++block)
      primal_blocks.emplace_back(product.layout().primal->dimension(block),
                                 value);
    std::vector<DenseVector> multiplier_blocks;
    multiplier_blocks.reserve(product.layout().multiplier->n_blocks());
    for (std::size_t block = 0;
         block < product.layout().multiplier->n_blocks();
         ++block)
      multiplier_blocks.emplace_back(
        product.layout().multiplier->dimension(block), value);
    return {Point{Product::Primal(product.layout().primal,
                                  std::move(primal_blocks)),
                  Product::Primal(product.layout().multiplier,
                                  std::move(multiplier_blocks))},
            {}};
  }

  Product
  make_nonfinite_residual_product(const Product &prototype,
                                  const bool     stationarity,
                                  const double   value)
  {
    const Product::Layout layout = prototype.layout();
    const auto quadratic_action =
      [layout, stationarity, value](const Product::Primal &) {
        const double q_value = stationarity ? value : 0.0;
        return Product::Covector(
          layout.stationarity,
          {DenseVector(layout.stationarity->dimension(0), q_value),
           DenseVector(layout.stationarity->dimension(1), q_value)});
      };
    const auto equality_action =
      [layout, stationarity, value](const Product::Primal &) {
        const double d_value = stationarity ? 0.0 : value;
        return Product::Covector(
          layout.equality,
          {DenseVector(layout.equality->dimension(0), d_value)});
      };
    const auto multiplier_action = [layout](const Product::Primal &) {
      return Product::Covector::zeros(layout.stationarity);
    };
    const auto transpose_action = [layout](const Product::Seed &) {
      return Product::TransposeResult{
        Product::Covector::zeros(layout.primal),
        Product::Covector::zeros(layout.multiplier)};
    };
    const Product::MultiplierConversion conversion{
      "non-finite residual test conversion",
      [layout](const Product::Primal &multiplier) {
        return Product::Primal(layout.adjoint, {multiplier.block(0)});
      },
      [layout](const Product::Primal &adjoint) {
        return Product::Primal(layout.multiplier, {adjoint.block(0)});
      }};
    return Product(layout,
                   quadratic_action,
                   equality_action,
                   multiplier_action,
                   transpose_action,
                   Product::Covector::zeros(layout.stationarity),
                   Product::Covector::zeros(layout.equality),
                   conversion,
                   prototype.assumptions(),
                   prototype.symmetry());
  }

  QuadraticKKTSolveResult
  zero_solution(const Product &product)
  {
    return {zero_point(product), {}};
  }

  void
  test_inactive_box_agrees_with_unconstrained_kkt()
  {
    const Problem problem;
    const auto complementarity = make_complementarity(problem, 2.0, 2.0);
    const PDASSolver solver(problem.product,
                            complementarity,
                            1,
                            solve_dense);
    const auto result = solver.solve(zero_point(problem.product),
                                     zero_box_multiplier(
                                       problem.control_layout),
                                     valid_policy());

    require(result.converged(),
            "PDAS did not converge for the inactive-box case");
    require(result.stopping_reason == PDASStoppingReason::converged,
            "PDAS returned the wrong inactive-box stopping reason");
    require(result.iterations.size() == 1,
            "inactive-box PDAS needed more than one KKT solve");
    require(result.iterations.front().selection.active_size() == 0,
            "inactive-box PDAS selected an active set");
    require(result.iterations.front().active_set_stable,
            "inactive-box PDAS did not report a stable set");
    require(result.iterations.front().kkt_residuals_converged,
            "inactive-box PDAS did not report full KKT convergence");
    require_close(result.solution.primal.block(1)[0], 1.0,
                  "inactive-box PDAS returned the wrong first control");
    require_close(result.solution.primal.block(1)[1], 0.125,
                  "inactive-box PDAS returned the wrong second control");
    require_close(result.box_multiplier.block(0)[0], 0.0,
                  "inactive-box PDAS returned a nonzero first multiplier");
    require_close(result.box_multiplier.block(0)[1], 0.0,
                  "inactive-box PDAS returned a nonzero second multiplier");
  }

  void
  test_active_box_restricted_kkt_reports_diagnostics()
  {
    const Problem problem;
    const auto complementarity = make_complementarity(problem, 0.5, 1.0);
    const PDASSolver solver(problem.product,
                            complementarity,
                            1,
                            solve_dense);
    const auto result = solver.solve(zero_point(problem.product),
                                     zero_box_multiplier(
                                       problem.control_layout),
                                     valid_policy());

    require(result.converged(),
            "PDAS did not converge for the active-box case");
    require(result.iterations.size() == 2,
            "active-box PDAS did not re-solve after changing its set");
    require(result.iterations[0].active_set_changes == 1,
            "active-box PDAS reported the wrong first set change count");
    require(result.iterations[1].active_set_stable,
            "active-box PDAS did not stabilize its active set");
    require(result.iterations[1].selection.activities() ==
              std::vector<BoxActivity>{BoxActivity::upper,
                                       BoxActivity::inactive},
            "active-box PDAS selected the wrong final set");
    require(result.iterations[1].primal_feasible &&
              result.iterations[1].dual_feasible &&
              result.iterations[1].complementarity_converged,
            "active-box PDAS did not report complementarity diagnostics");
    require(result.iterations[1].kkt_residuals_converged,
            "active-box PDAS did not report full KKT convergence");
    require_close(result.solution.primal.block(0)[0], 0.5,
                  "active-box PDAS returned the wrong active state");
    require_close(result.solution.primal.block(0)[1], 0.125,
                  "active-box PDAS returned the wrong inactive state");
    require_close(result.solution.primal.block(1)[0], 0.5,
                  "active-box PDAS returned the wrong active control");
    require_close(result.solution.primal.block(1)[1], 0.125,
                  "active-box PDAS returned the wrong inactive control");
    require_close(result.box_multiplier.block(0)[0], 1.0,
                  "active-box PDAS returned the wrong active multiplier");
    require_close(result.box_multiplier.block(0)[1], 0.0,
                  "active-box PDAS returned a nonzero inactive multiplier");
  }

  void
  test_invalid_pdas_inputs()
  {
    const Problem problem;
    const auto complementarity = make_complementarity(problem, 0.5, 1.0);
    const PDASSolver solver(problem.product,
                            complementarity,
                            1,
                            solve_dense);
    PDASPolicy invalid;
    invalid.maximum_iterations = 0;
    nmopt::test_support::require_contract_error(
      [&] {
        (void)solver.solve(zero_point(problem.product),
                           zero_box_multiplier(problem.control_layout),
                           invalid);
      },
      "PDAS policy is invalid",
      "PDAS accepted a zero iteration limit");

    const auto infeasible = Point{
      Product::Primal(problem.product.layout().primal,
                      {DenseVector{0.0, 0.0}, DenseVector{0.0, 2.0}}),
      Product::Primal::zeros(problem.product.layout().multiplier)};
    nmopt::test_support::require_contract_error(
      [&] {
        (void)solver.solve(infeasible,
                           zero_box_multiplier(problem.control_layout),
                           valid_policy());
      },
      "PDAS initial control must be feasible",
      "PDAS accepted an infeasible initial control");

    PDASPolicy undeclared_assumptions;
    nmopt::test_support::require_contract_error(
      [&] {
        (void)solver.solve(zero_point(problem.product),
                           zero_box_multiplier(problem.control_layout),
                           undeclared_assumptions);
      },
      "Quadratic KKT product needs a declared rank condition",
      "PDAS constructed an active KKT product without assumptions");
  }

  void
  test_detached_active_product_owns_builder_state()
  {
    const Product detached = [] {
      Problem problem;
      const auto complementarity = make_complementarity(problem, 0.5, 1.0);
      const ActiveSetSelection selection(
        problem.control_layout,
        {BoxActivity::upper, BoxActivity::inactive});
      const ActiveSetKKTSubproblem subproblem(problem.product,
                                              complementarity,
                                              selection,
                                              1,
                                              active_set_assumptions());
      return subproblem.product();
    }();

    const auto residual = detached.residual(zero_point(detached));
    require(residual.stationarity.n_blocks() == 2,
            "detached active product lost its stationarity layout");
    const auto zero_primal = Product::Primal::zeros(detached.layout().primal);
    const auto zero_multiplier =
      Product::Primal::zeros(detached.layout().multiplier);
    const auto quadratic = detached.apply_q(zero_primal);
    const auto equality = detached.apply_d(zero_primal);
    const auto multiplier = detached.apply_d_transpose(zero_multiplier);
    require(quadratic.layout()->compatible_with(*detached.layout().stationarity) &&
              equality.layout()->compatible_with(*detached.layout().equality) &&
              multiplier.layout()->compatible_with(
                *detached.layout().stationarity),
            "detached active product lost a restricted KKT action");
    const auto full_action = detached.apply_kkt(zero_point(detached));
    require(full_action.stationarity.layout()->compatible_with(
              *detached.layout().stationarity),
            "detached active product lost its full KKT action");
    const Product::Seed seed{
      Product::Primal::zeros(detached.layout().stationarity),
      Product::Primal::zeros(detached.layout().equality)};
    const auto transpose = detached.apply_kkt_transpose(seed);
    require(transpose.primal.layout()->compatible_with(
              *detached.layout().primal),
            "detached active product lost its transpose action");
    const auto adjoint = detached.multiplier_to_adjoint(
      Product::Primal::zeros(detached.layout().multiplier));
    require(adjoint.layout()->compatible_with(*detached.layout().adjoint),
            "detached active product lost its multiplier conversion");
    const auto solved = solve_dense(detached);
    require(solved.report.converged(),
            "detached active product dense solve did not converge");
  }

  void
  test_detached_pdas_solver_owns_inputs()
  {
    struct DetachedSolver
    {
      Product    product;
      PDASSolver solver;
      Covector   box_multiplier;
    };

    const auto detached = [] {
      Problem problem;
      const auto complementarity = make_complementarity(problem, 0.5, 1.0);
      return DetachedSolver{problem.product,
                            PDASSolver(problem.product,
                                       complementarity,
                                       1,
                                       solve_dense),
                            zero_box_multiplier(problem.control_layout)};
    }();

    const auto result = detached.solver.solve(
      zero_point(detached.product), detached.box_multiplier, valid_policy());
    require(result.converged(),
            "detached PDAS solver did not converge after its inputs were released");
    require_close(result.solution.primal.block(1)[0], 0.5,
                  "detached PDAS solver returned the wrong active control");
  }

  void
  test_nonfinite_inputs_and_pdas_outputs_are_rejected()
  {
    const std::vector<double> nonfinite_values{
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity()};
    const Problem problem;
    const auto complementarity = make_complementarity(problem, 2.0, 2.0);
    const PDASSolver solver(problem.product,
                            complementarity,
                            1,
                            solve_dense);

    for (const double value : nonfinite_values)
      {
        const Point nonfinite_initial{
          Product::Primal(problem.product.layout().primal,
                          {DenseVector{0.0, 0.0}, DenseVector{value, 0.0}}),
          Product::Primal::zeros(problem.product.layout().multiplier)};
        nmopt::test_support::require_contract_error(
          [&] {
            (void)solver.solve(nonfinite_initial,
                               zero_box_multiplier(problem.control_layout),
                               valid_policy());
          },
          "PDAS initial primal point contains a non-finite value",
          "PDAS accepted a non-finite initial control");

        const Covector nonfinite_box_multiplier(
          problem.control_layout, {DenseVector{value, 0.0}});
        nmopt::test_support::require_contract_error(
          [&] {
            (void)solver.solve(zero_point(problem.product),
                               nonfinite_box_multiplier,
                               valid_policy());
          },
          "PDAS initial box multiplier contains a non-finite value",
          "PDAS accepted a non-finite initial box multiplier");

        const PDASSolver nonfinite_solver(
          problem.product,
          complementarity,
          1,
          [value](const Product &product) {
            return nonfinite_solution(product, value);
          });
        nmopt::test_support::require_contract_error(
          [&] {
            (void)nonfinite_solver.solve(
              zero_point(problem.product),
              zero_box_multiplier(problem.control_layout),
              valid_policy());
          },
          "PDAS KKT solution primal contains a non-finite value",
          "PDAS accepted a non-finite KKT solution");
      }

    for (const bool stationarity : {true, false})
      {
        const Product nonfinite_product =
          make_nonfinite_residual_product(problem.product,
                                          stationarity,
                                          std::numeric_limits<double>::quiet_NaN());
        const PDASSolver nonfinite_solver(nonfinite_product,
                                          complementarity,
                                          1,
                                          zero_solution);
        nmopt::test_support::require_contract_error(
          [&] {
            (void)nonfinite_solver.solve(
              zero_point(nonfinite_product),
              zero_box_multiplier(problem.control_layout),
              valid_policy());
          },
          stationarity ? "PDAS stationarity residual contains a non-finite value"
                       : "PDAS equality residual contains a non-finite value",
          "PDAS reported a non-finite residual block");
      }

    const PDASSolver nonfinite_report_solver(
      problem.product,
      complementarity,
      1,
      [](const Product &product) {
        auto result = zero_solution(product);
        result.report.stationarity_residual =
          std::numeric_limits<double>::quiet_NaN();
        return result;
      });
    nmopt::test_support::require_contract_error(
      [&] {
        (void)nonfinite_report_solver.solve(
          zero_point(problem.product),
          zero_box_multiplier(problem.control_layout),
          valid_policy());
      },
      "PDAS KKT stationarity residual is non-finite",
      "PDAS accepted a non-finite KKT solve report");
  }

  void
  test_pdas_stopping_reasons()
  {
    const Problem problem;
    const auto complementarity = make_complementarity(problem, 0.5, 1.0);
    PDASPolicy one_iteration = valid_policy();
    one_iteration.maximum_iterations = 1;
    const PDASSolver solver(problem.product,
                            complementarity,
                            1,
                            solve_dense);
    const auto maximum_result = solver.solve(
      zero_point(problem.product),
      zero_box_multiplier(problem.control_layout),
      one_iteration);
    require(maximum_result.stopping_reason ==
              PDASStoppingReason::maximum_iterations,
            "PDAS returned the wrong maximum-iteration stopping reason");
    require(maximum_result.iterations.size() == 1,
            "PDAS returned the wrong maximum-iteration report count");

    const PDASSolver failing_solver(problem.product,
                                    complementarity,
                                    1,
                                    nonconverged_solve);
    const auto failed_result = failing_solver.solve(
      zero_point(problem.product),
      zero_box_multiplier(problem.control_layout),
      valid_policy());
    require(failed_result.stopping_reason ==
              PDASStoppingReason::kkt_solve_failed,
            "PDAS returned the wrong KKT-failure stopping reason");
    require(!failed_result.iterations.front().kkt_solve.converged(),
            "PDAS hid the failed KKT solve in its iteration report");
  }
} // namespace

int
main(const int argc, char **argv)
{
  try
    {
      const std::vector<nmopt::test_support::Scenario> scenarios{
        {"active_subproblem_restricts_free_coordinates_and_reconstructs_solution",
         "nmopt.active_set_kkt.active_subproblem_restricts_free_coordinates_and_reconstructs_solution",
         {"backend-neutral", "formulation", "kkt", "active-set"},
         30,
         test_active_subproblem_restricts_free_coordinates_and_reconstructs_solution},
        {"inactive_subproblem_reuses_base_product",
         "nmopt.active_set_kkt.inactive_subproblem_reuses_base_product",
         {"backend-neutral", "formulation", "kkt", "active-set"},
         30,
         test_inactive_subproblem_reuses_base_product},
        {"invalid_active_subproblem_inputs",
         "nmopt.active_set_kkt.invalid_active_subproblem_inputs",
         {"backend-neutral", "formulation", "kkt", "active-set", "negative"},
         30,
         test_invalid_active_subproblem_inputs},
        {"inactive_box_agrees_with_unconstrained_kkt",
         "nmopt.pdas.inactive_box_agrees_with_unconstrained_kkt",
         {"backend-neutral", "formulation", "kkt", "pdas"},
         30,
         test_inactive_box_agrees_with_unconstrained_kkt},
        {"active_box_restricted_kkt_reports_diagnostics",
         "nmopt.pdas.active_box_restricted_kkt_reports_diagnostics",
         {"backend-neutral", "formulation", "kkt", "pdas", "active-set"},
         30,
         test_active_box_restricted_kkt_reports_diagnostics},
        {"invalid_pdas_inputs",
         "nmopt.pdas.invalid_pdas_inputs",
         {"backend-neutral", "formulation", "kkt", "pdas", "negative"},
         30,
         test_invalid_pdas_inputs},
        {"detached_active_product_owns_builder_state",
         "nmopt.active_set_kkt.detached_active_product_owns_builder_state",
         {"backend-neutral", "formulation", "kkt", "active-set", "ownership"},
         30,
         test_detached_active_product_owns_builder_state},
        {"detached_pdas_solver_owns_inputs",
         "nmopt.pdas.detached_pdas_solver_owns_inputs",
         {"backend-neutral", "formulation", "kkt", "pdas", "ownership"},
         30,
         test_detached_pdas_solver_owns_inputs},
        {"nonfinite_inputs_and_pdas_outputs_are_rejected",
         "nmopt.pdas.nonfinite_inputs_and_pdas_outputs_are_rejected",
         {"backend-neutral", "formulation", "kkt", "pdas", "negative"},
         30,
         test_nonfinite_inputs_and_pdas_outputs_are_rejected},
        {"pdas_stopping_reasons",
         "nmopt.pdas.pdas_stopping_reasons",
         {"backend-neutral", "formulation", "kkt", "pdas", "solver"},
         30,
         test_pdas_stopping_reasons}};
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
