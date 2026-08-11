#include "nmopt/dealii/trace_hhalf_metric.hpp"

#include "test_support/contract_errors.hpp"
#include "test_support/scenario_dispatch.hpp"

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>

#include <cmath>
#include <exception>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace
{
  using Backend = nmopt::dealii_backend::SerialBackend;
  using Primal = nmopt::contract::PrimalBlockT<Backend>;
  using Covector = nmopt::contract::CovectorBlockT<Backend>;

  void
  require_close(const double       actual,
                const double       expected,
                const double       tolerance,
                const std::string &description)
  {
    nmopt::contract::require(std::abs(actual - expected) <= tolerance,
                             description + " differs from its oracle");
  }

  void
  require_primal_close(const Primal &     actual,
                       const Primal &     expected,
                       const double       tolerance,
                       const std::string &description)
  {
    nmopt::contract::require_compatible(
      actual, expected, description + " has incompatible layouts");
    dealii::Vector<double> difference = actual.block(0);
    difference.add(-1.0, expected.block(0));
    require_close(difference.l2_norm(), 0.0, tolerance, description);
  }

  struct MatrixStorage
  {
    std::shared_ptr<dealii::SparsityPattern>     sparsity;
    std::shared_ptr<dealii::SparseMatrix<double>> matrix;
  };

  MatrixStorage
  volume_h1_matrix()
  {
    dealii::DynamicSparsityPattern dynamic(3, 3);
    for (std::size_t row = 0; row < 3; ++row)
      for (std::size_t column = 0; column < 3; ++column)
        dynamic.add(row, column);
    auto sparsity = std::make_shared<dealii::SparsityPattern>();
    sparsity->copy_from(dynamic);
    auto matrix = std::make_shared<dealii::SparseMatrix<double>>(*sparsity);
    matrix->set(0, 0, 4.0);
    matrix->set(0, 1, 1.0);
    matrix->set(0, 2, 1.0);
    matrix->set(1, 0, 1.0);
    matrix->set(1, 1, 3.0);
    matrix->set(1, 2, 1.0);
    matrix->set(2, 0, 1.0);
    matrix->set(2, 1, 1.0);
    matrix->set(2, 2, 2.0);
    return {std::move(sparsity), std::move(matrix)};
  }

  void
  run_trace_hhalf_metric_contract()
  {
    const auto layout = std::make_shared<const nmopt::contract::BlockLayout>(
      "trace_hhalf", std::vector<nmopt::contract::SpaceId>{{"control"}},
      std::vector<std::size_t>{2});
    nmopt::dealii_backend::MetricSolveParameters solve_parameters;
    solve_parameters.relative_tolerance = 1e-13;
    solve_parameters.absolute_tolerance = 1e-15;
    const auto volume_h1 = volume_h1_matrix();
    const nmopt::dealii_backend::TraceHhalfMetric metric(
      "hhalf_trace", layout, volume_h1.matrix, {0, 2}, solve_parameters);

    dealii::Vector<double> values(2);
    values[0] = 2.0;
    values[1] = -1.0;
    const Primal primal(layout, {values});
    const Covector action = metric.apply(primal);

    // Eliminating volume DoF 1 gives the exact Schur complement
    // [[11/3, 2/3], [2/3, 5/3]].
    require_close(action.block(0)[0],
                  20.0 / 3.0,
                  1e-12,
                  "H1/2 Schur-complement first row");
    require_close(action.block(0)[1],
                  -1.0 / 3.0,
                  1e-12,
                  "H1/2 Schur-complement second row");
    require_primal_close(metric.inverse_apply(action),
                         primal,
                         1e-12,
                         "H1/2 metric round trip");

    dealii::Vector<double> other_values(2);
    other_values[0] = -0.5;
    other_values[1] = 1.25;
    const Primal other(layout, {other_values});
    const Covector other_action = metric.apply(other);
    require_close(nmopt::contract::pair(action, other),
                  nmopt::contract::pair(other_action, primal),
                  1e-12,
                  "H1/2 metric symmetry");
    nmopt::contract::require(nmopt::contract::pair(action, primal) > 0.0,
                             "H1/2 metric is not positive on a nonzero trace");

    auto duplicate_layout =
      std::make_shared<const nmopt::contract::BlockLayout>(
        "duplicate_trace",
        std::vector<nmopt::contract::SpaceId>{{"control"}},
        std::vector<std::size_t>{2});
    nmopt::test_support::require_contract_error(
      [&duplicate_layout, &volume_h1]() {
        (void)nmopt::dealii_backend::TraceHhalfMetric(
          "duplicate", duplicate_layout, volume_h1.matrix, {0, 0});
      },
      "H1/2 trace DoF map contains a duplicate",
      "duplicate H1/2 trace map");
  }
} // namespace

int
main(const int argc, char **argv)
{
  try
    {
      const std::vector<nmopt::test_support::Scenario> scenarios{
        {"trace_hhalf_metric",
         "nmopt.dealii.trace_hhalf_metric",
         {"dealii", "contract", "metric"},
         30,
         run_trace_hhalf_metric_contract}};
      const auto result = nmopt::test_support::run_requested_scenarios(
        argc, argv, scenarios, std::cout);
      if (!result.listed)
        std::cout << "deal.II H1/2 trace metric contract scenario passed: "
                  << result.executed << '\n';
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "deal.II H1/2 trace metric contract test failed: "
                << exception.what() << '\n';
      return 1;
    }
}
