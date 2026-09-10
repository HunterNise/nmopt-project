#define STEP4_NO_MAIN
#include "../../apps/external-dealii/step-4/step-4.cc"
#undef STEP4_NO_MAIN

#include "../support/scenario_dispatch.hpp"

#include <deal.II/lac/vector.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
  using Matrix = dealii::SparseMatrix<double>;
  using Vector = dealii::Vector<double>;

  void
  require(const bool condition, const std::string &message)
  {
    if (!condition)
      throw std::runtime_error(message);
  }

  void
  require_close(const double       actual,
                const double       expected,
                const double       tolerance,
                const std::string &message)
  {
    if (std::abs(actual - expected) > tolerance)
      {
        std::ostringstream detail;
        detail << message << ": actual=" << actual
               << ", expected=" << expected
               << ", tolerance=" << tolerance;
        throw std::runtime_error(detail.str());
      }
  }

  double
  matrix_difference(const std::vector<double> &left, const Matrix &right)
  {
    require(left.size() ==
              static_cast<std::size_t>(right.m()) * right.n(),
            "Step4 matrix dimensions changed");

    double squared_difference = 0.0;
    for (unsigned int row = 0; row < right.m(); ++row)
      for (unsigned int column = 0; column < right.n(); ++column)
        {
          const double difference =
            left[static_cast<std::size_t>(row) * right.n() + column] -
            right.el(row, column);
          squared_difference += difference * difference;
        }
    return std::sqrt(squared_difference);
  }

  std::vector<double>
  matrix_values(const Matrix &matrix)
  {
    std::vector<double> values(static_cast<std::size_t>(matrix.m()) *
                               matrix.n());
    for (unsigned int row = 0; row < matrix.m(); ++row)
      for (unsigned int column = 0; column < matrix.n(); ++column)
        values[static_cast<std::size_t>(row) * matrix.n() + column] =
          matrix.el(row, column);
    return values;
  }

  double
  matrix_symmetry_error(const Matrix &matrix)
  {
    double squared_difference = 0.0;
    double squared_norm       = 0.0;
    for (unsigned int row = 0; row < matrix.m(); ++row)
      for (unsigned int column = 0; column < matrix.n(); ++column)
        {
          const double value = matrix.el(row, column);
          squared_norm += value * value;
          const double difference = value - matrix.el(column, row);
          squared_difference += difference * difference;
        }
    return std::sqrt(squared_difference) /
           std::max(1.0, std::sqrt(squared_norm));
  }

  double
  vector_difference(const Vector &left, const Vector &right)
  {
    require(left.size() == right.size(), "Step4 vector dimensions changed");
    Vector difference = left;
    difference.add(-1.0, right);
    return difference.l2_norm();
  }

  void
  require_original_assembly(const std::vector<double> &matrix_before,
                            const Vector &             rhs_before,
                            const Step4<2> &           tutorial)
  {
    require_close(matrix_difference(matrix_before,
                                    tutorial.system_matrix_view()),
                  0.0,
                  0.0,
                  "Step4 assembled matrix changed during solve");
    require_close(vector_difference(rhs_before, tutorial.system_rhs_view()),
                  0.0,
                  0.0,
                  "Step4 assembled RHS changed during solve");
  }

  void
  run_matrix_rhs_solve_contract()
  {
    Step4<2> tutorial;
    tutorial.prepare_for_external_use();

    const auto &matrix = tutorial.system_matrix_view();
    const auto &rhs    = tutorial.system_rhs_view();
    require(matrix.m() == 289 && matrix.n() == 289,
            "Step4 2D assembled matrix has the wrong dimension");
    require(rhs.size() == 289, "Step4 2D assembled RHS has the wrong dimension");
    require(matrix_symmetry_error(matrix) <= 1e-14,
            "Step4 2D assembled matrix is not symmetric");

    const auto matrix_before = matrix_values(matrix);
    const Vector rhs_before  = rhs;

    Vector state(rhs.size());
    state = 0.0;
    const auto evidence = tutorial.solve(rhs, state);
    require(evidence.converged, "Step4 supplied-RHS solve did not converge");
    require(evidence.iterations > 0,
            "Step4 supplied-RHS solve forced no work for a nonzero RHS");
    require(evidence.final_residual <= 1e-10,
            "Step4 supplied-RHS solve has a large residual");
    require_original_assembly(matrix_before, rhs_before, tutorial);

    Vector second_state(rhs.size());
    second_state = 0.5;
    const auto second_evidence = tutorial.solve(rhs, second_state);
    require(second_evidence.converged,
            "Step4 nonzero-initial supplied-RHS solve did not converge");
    require(second_evidence.final_residual <= 1e-10,
            "Step4 nonzero-initial solve has a large residual");
    require_original_assembly(matrix_before, rhs_before, tutorial);

    Vector zero_rhs(rhs.size());
    zero_rhs = 0.0;
    Vector zero_state(rhs.size());
    zero_state = 0.0;
    const auto zero_evidence = tutorial.solve(zero_rhs, zero_state);
    require(zero_evidence.converged,
            "Step4 zero-RHS zero-initial solve did not converge");
    require(zero_evidence.iterations == 0,
            "Step4 zero-RHS zero-initial solve forced iterations");
    require(zero_evidence.final_residual == 0.0,
            "Step4 zero-RHS zero-initial solve changed the solution");
    require_original_assembly(matrix_before, rhs_before, tutorial);
  }

  void
  run_supplied_state_output_contract()
  {
    Step4<2> tutorial;
    tutorial.prepare_for_external_use();

    Vector state(tutorial.system_rhs_view().size());
    state = 7.0;

    const auto output_directory =
      std::filesystem::temp_directory_path() / "nmopt-external-step4-output";
    std::filesystem::remove_all(output_directory);
    std::filesystem::create_directories(output_directory);
    const auto output_file = output_directory / "supplied-state.vtk";
    tutorial.output_results(state, output_file);

    std::ifstream input(output_file);
    require(static_cast<bool>(input),
            "Step4 supplied-state output could not be opened");
    const std::string contents((std::istreambuf_iterator<char>(input)), {});
    const auto lookup_table = contents.find("LOOKUP_TABLE default");
    require(contents.find("# vtk DataFile Version") != std::string::npos,
            "Step4 supplied-state output is not legacy VTK");
    require(contents.find("SCALARS solution double 1") != std::string::npos,
            "Step4 supplied-state output changed the field identity");
    require(lookup_table != std::string::npos &&
              contents.find("7", lookup_table) != std::string::npos,
            "Step4 supplied-state output did not use the supplied vector");

    std::filesystem::remove_all(output_directory);
  }
} // namespace

int
main(const int argc, char **argv)
{
  try
    {
      const std::vector<nmopt::test_support::Scenario> scenarios{
        {"matrix_rhs_solve",
         "nmopt.external_tutorial_step_4.native_matrix_rhs_solve",
         {"dealii", "application", "external", "tutorial", "reuse"},
         60,
         run_matrix_rhs_solve_contract},
        {"supplied_state_output",
         "nmopt.external_tutorial_step_4.native_supplied_state_output",
         {"dealii", "application", "external", "tutorial", "reuse"},
         30,
         run_supplied_state_output_contract}};
      const auto result = nmopt::test_support::run_requested_scenarios(
        argc, argv, scenarios, std::cout);
      if (!result.listed)
        std::cout << "external Step4 native reuse scenarios passed: "
                  << result.executed << '\n';
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "external Step4 native reuse test failed: "
                << exception.what() << '\n';
      return 1;
    }
}
