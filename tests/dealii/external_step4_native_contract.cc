#include "../../apps/external-dealii/step-4/evaluation/native_reduced.hpp"
#include "../../apps/external-dealii/step-4/evaluation/scenario.hpp"

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
  using Matrix          = external_dealii_step4::ProblemA::Matrix;
  using ProblemA        = external_dealii_step4::ProblemA;
  using NativeReduced   = external_dealii_step4::NativeReduced;
  using Instrumentation = external_dealii_step4::Instrumentation;
  using Vector          = ProblemA::Vector;

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

  void
  run_problem_a_operations_contract()
  {
    Instrumentation instrumentation;
    ProblemA        problem(instrumentation);
    require(problem.state_dimension() ==
              external_dealii_step4::scenario::dimension,
            "Problem A state dimension does not match the frozen scenario");
    require(problem.control_dimension() ==
              external_dealii_step4::scenario::dimension,
            "Problem A control dimension does not match the frozen scenario");
    require(instrumentation.assembly_calls == 1,
            "Problem A did not assemble exactly once during construction");

    const auto controls = external_dealii_step4::scenario::reduced_controls();
    require(controls.size() == 4, "Problem A scenario controls are incomplete");
    for (const auto &control : controls)
      require(control.size() == problem.control_dimension(),
              "Problem A scenario control has the wrong dimension");

    const auto state_result = problem.solve_state(controls.front());
    require(state_result.evidence.converged,
            "Problem A zero-control state solve did not converge");
    require(state_result.evidence.final_residual <= 1e-10,
            "Problem A zero-control state solve has a large residual");

    const auto objective =
      problem.objective(state_result.solution, controls.front());
    require(objective >= 0.0, "Problem A objective is negative");
    const auto objective_derivative =
      problem.objective_derivative(state_result.solution, controls.front());
    require(vector_difference(objective_derivative.state, state_result.solution) ==
              0.0,
            "Problem A state objective derivative has the wrong value");
    require(vector_difference(objective_derivative.control, controls.front()) ==
              0.0,
            "Problem A control objective derivative has the wrong value");

    const auto adjoint_result =
      problem.solve_adjoint(objective_derivative.state);
    require(adjoint_result.evidence.converged,
            "Problem A adjoint solve did not converge");
    require(adjoint_result.evidence.final_residual <= 1e-10,
            "Problem A adjoint solve has a large residual");

    Vector seed(problem.state_dimension());
    seed = 0.25;
    const auto control_pullback = problem.control_vjp(seed);
    seed *= -1.0;
    require(vector_difference(control_pullback, seed) == 0.0,
            "Problem A control pullback has the wrong sign");
  }

  void
  run_native_reduced_contract()
  {
    Instrumentation instrumentation;
    ProblemA        problem(instrumentation);
    NativeReduced  reduced(problem, instrumentation);
    const auto     controls = external_dealii_step4::scenario::reduced_controls();

    const auto value = reduced.evaluate_value(controls.front());
    require(instrumentation.value_evaluations == 1,
            "native reduced value stage was not counted");
    require(instrumentation.derivative_augmentations == 0,
            "native reduced derivative ran during value evaluation");
    require(instrumentation.state_solve_calls == 1,
            "native reduced value stage did not perform one state solve");
    require(instrumentation.adjoint_solve_calls == 0,
            "native reduced value stage performed an adjoint solve");
    require(instrumentation.objective_calls == 1,
            "native reduced value stage did not evaluate the objective once");

    const auto derivative = reduced.augment_derivative(value);
    require(instrumentation.derivative_augmentations == 1,
            "native reduced derivative stage was not counted");
    require(instrumentation.state_solve_calls == 1,
            "native reduced derivative stage repeated the state solve");
    require(instrumentation.adjoint_solve_calls == 1,
            "native reduced derivative stage did not perform one adjoint solve");
    require(instrumentation.objective_derivative_calls == 1,
            "native reduced derivative stage did not evaluate objective partials once");
    require(instrumentation.control_vjp_calls == 1,
            "native reduced derivative stage did not perform one control pullback");

    Vector expected_gradient = derivative.control_derivative;
    expected_gradient += derivative.adjoint;
    require(vector_difference(derivative.reduced_derivative, expected_gradient) ==
              0.0,
            "native reduced derivative has the wrong control sign");
    require(vector_difference(value.state, derivative.state_derivative) == 0.0,
            "native reduced derivative changed the retained state");

    const auto repeated_before = reduced.evaluate_value(controls[2]);
    const auto intervening      = reduced.evaluate_value(controls[3]);
    const auto repeated_after   = reduced.evaluate_value(controls[2]);
    (void)intervening;
    require(vector_difference(repeated_before.state, repeated_after.state) == 0.0,
            "native reduced repeated control changed the state");
    require_close(repeated_before.objective,
                  repeated_after.objective,
                  0.0,
                  "native reduced repeated control changed the objective");

    const auto output_directory =
      std::filesystem::temp_directory_path() /
      "nmopt-external-step4-native-reduced-output";
    std::filesystem::remove_all(output_directory);
    std::filesystem::create_directories(output_directory);
    const auto output_file = output_directory / "retained-state.vtk";
    problem.output_results(repeated_after.state, output_file);
    require(std::filesystem::exists(output_file),
            "native reduced retained-state output is missing");
    require(instrumentation.output_calls == 1,
            "native reduced output was not counted");
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
         run_supplied_state_output_contract},
        {"problem_a_operations",
         "nmopt.external_tutorial_step_4.native_problem_a_operations",
         {"dealii", "application", "external", "tutorial", "native", "control"},
         60,
         run_problem_a_operations_contract},
        {"native_reduced",
         "nmopt.external_tutorial_step_4.native_reduced",
         {"dealii", "application", "external", "tutorial", "native", "control"},
         60,
         run_native_reduced_contract}};
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
