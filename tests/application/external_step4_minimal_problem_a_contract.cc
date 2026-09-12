#include "../../apps/external-dealii/step-4/evaluation/native_optimization.hpp"
#include "../../apps/external-dealii/step-4/evaluation/optimization_policy.hpp"
#include "../../apps/external-dealii/step-4/minimal/problem_a_binding.hpp"
#include "../../apps/external-dealii/step-4/verification/verification.hpp"

#include "../dealii/external_step4_evidence.hpp"
#include "../support/scenario_dispatch.hpp"

#include <deal.II/lac/vector.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "nmopt/solvers/reduced_gradient.hpp"

namespace
{
  using Backend       = nmopt::dealii_backend::SerialBackend;
  using Binding       = external_dealii_step4::minimal::ProblemABinding;
  using Instrumentation = external_dealii_step4::Instrumentation;
  using NativeReduced = external_dealii_step4::NativeReduced;
  using NativeSolver  = external_dealii_step4::NativeArmijoSolver;
  using OptimizationPolicy = external_dealii_step4::OptimizationPolicy;
  using Parameters    = nmopt::solvers::ReducedSolverParameters;
  using ProblemA      = external_dealii_step4::ProblemA;
  using Matrix        = ProblemA::Matrix;
  using Solver        = nmopt::solvers::ReducedSearchSolverT<Backend>;
  using Vector        = ProblemA::Vector;
  namespace verification = external_dealii_step4::verification;

  struct FreshAudit
  {
    Vector state;
    Vector adjoint;
    Vector gradient;
  };

  struct ConsumerSummary
  {
    bool        completed = false;
    bool        has_stopping_reason = false;
    std::string stopping_reason;
    bool        has_accepted_iterations = false;
    std::size_t accepted_iterations = 0;
    bool        has_line_search_trials = false;
    std::size_t line_search_trials = 0;
    bool        has_objective = false;
    double      objective = 0.0;
    bool        has_gradient = false;
    double      gradient = 0.0;
  };

  void
  require(const bool condition, const std::string &message)
  {
    if (!condition)
      {
        external_dealii_step4_test::note_failure(message);
        throw std::runtime_error(message);
      }
  }

  Parameters
  nmopt_parameters(const OptimizationPolicy &policy)
  {
    Parameters parameters;
    parameters.maximum_iterations = policy.maximum_iterations;
    parameters.maximum_line_search_trials = policy.maximum_line_search_trials;
    parameters.gradient_tolerance = policy.gradient_tolerance;
    parameters.stopping_criterion =
      nmopt::solvers::ReducedStoppingCriterion::gradient_norm;
    parameters.relative_gradient_tolerance =
      policy.relative_gradient_tolerance;
    parameters.objective_change_tolerance = policy.objective_change_tolerance;
    parameters.step_tolerance = policy.step_tolerance;
    parameters.objective_target = policy.objective_target;
    parameters.initial_step_length = policy.initial_step_length;
    parameters.minimum_step_length = policy.minimum_step_length;
    parameters.armijo_fraction = policy.armijo_fraction;
    parameters.backtracking_factor = policy.backtracking_factor;
    return parameters;
  }

  double
  vector_difference(const Vector &left, const Vector &right)
  {
    require(left.size() == right.size(),
            "minimal Problem A vectors have incompatible dimensions");
    Vector difference = left;
    difference.add(-1.0, right);
    return difference.l2_norm();
  }

  double
  matrix_difference(const Matrix &left, const Matrix &right)
  {
    require(left.m() == right.m() && left.n() == right.n(),
            "minimal Problem A matrices have incompatible shapes");
    double squared_difference = 0.0;
    for (unsigned int row = 0; row < left.m(); ++row)
      for (unsigned int column = 0; column < left.n(); ++column)
        {
          const double difference =
            left.el(row, column) - right.el(row, column);
          squared_difference += difference * difference;
        }
    return std::sqrt(squared_difference);
  }

  bool
  scalar_matches(const double left, const double right)
  {
    return std::isfinite(left) && std::isfinite(right) &&
           std::abs(left - right) <=
             1.0e-12 + 1.0e-10 * std::max(std::abs(left), std::abs(right));
  }

  FreshAudit
  fresh_audit(const Vector &control)
  {
    ProblemA fresh_problem;
    auto state_result = fresh_problem.solve_state(control);
    require(state_result.evidence.converged,
            "minimal Problem A fresh state solve did not converge");
    const auto objective_derivative =
      fresh_problem.objective_derivative(state_result.solution, control);
    auto adjoint = fresh_problem.solve_adjoint(objective_derivative.state);
    require(adjoint.evidence.converged,
            "minimal Problem A fresh adjoint solve did not converge");
    const auto control_pullback = fresh_problem.control_vjp(adjoint.solution);
    Vector gradient = objective_derivative.control;
    gradient.add(-1.0, control_pullback);
    return {std::move(state_result.solution),
            std::move(adjoint.solution),
            std::move(gradient)};
  }

  std::string
  shell_quote(const std::filesystem::path &path)
  {
    std::string quoted = "'";
    for (const char character : path.string())
      if (character == '\'')
        quoted += "'\\''";
      else
        quoted += character;
    quoted += "'";
    return quoted;
  }

  ConsumerSummary
  read_consumer_summary(const std::filesystem::path &logfile)
  {
    std::ifstream input(logfile);
    require(static_cast<bool>(input),
            "minimal Problem A executable output could not be opened");

    ConsumerSummary summary;
    std::string line;
    while (std::getline(input, line))
      {
        if (line == "minimal Problem A consumer completed")
          summary.completed = true;

        std::istringstream stream(line);
        std::string key;
        stream >> key;
        if (key == "stopping_reason")
          {
            require(static_cast<bool>(stream >> summary.stopping_reason),
                    "minimal Problem A executable stopping reason is malformed");
            summary.has_stopping_reason = true;
          }
        else if (key == "accepted_iterations")
          {
            require(static_cast<bool>(stream >> summary.accepted_iterations),
                    "minimal Problem A executable iteration count is malformed");
            summary.has_accepted_iterations = true;
          }
        else if (key == "line_search_trials")
          {
            require(static_cast<bool>(stream >> summary.line_search_trials),
                    "minimal Problem A executable trial count is malformed");
            summary.has_line_search_trials = true;
          }
        else if (key == "final_objective")
          {
            require(static_cast<bool>(stream >> summary.objective),
                    "minimal Problem A executable objective is malformed");
            summary.has_objective = true;
          }
        else if (key == "final_gradient_norm")
          {
            require(static_cast<bool>(stream >> summary.gradient),
                    "minimal Problem A executable gradient is malformed");
            summary.has_gradient = true;
          }
      }

    require(summary.completed,
            "minimal Problem A executable did not report completion");
    require(summary.has_stopping_reason && summary.has_accepted_iterations &&
              summary.has_line_search_trials && summary.has_objective &&
              summary.has_gradient,
            "minimal Problem A executable report is incomplete");
    require(std::isfinite(summary.objective) && std::isfinite(summary.gradient),
            "minimal Problem A executable report is not finite");
    return summary;
  }

  std::string
  output_payload(const std::filesystem::path &filename)
  {
    std::ifstream input(filename, std::ios::binary);
    require(static_cast<bool>(input),
            "minimal Problem A output could not be opened");
    std::string contents{std::istreambuf_iterator<char>(input),
                         std::istreambuf_iterator<char>()};
    const std::string generated_header =
      "#This file was generated by the deal.II library on ";
    const auto header_start = contents.find(generated_header);
    if (header_start != std::string::npos)
      {
        const auto line_end = contents.find('\n', header_start);
        contents.erase(header_start,
                       line_end == std::string::npos ?
                         contents.size() - header_start :
                         line_end + 1 - header_start);
      }
    return contents;
  }

  void
  check_final_solution(ProblemA &             problem,
                       const Vector &         state,
                       const Vector &         control,
                       const Vector &         returned_adjoint,
                       const Vector &         returned_gradient,
                       const verification::OracleResult &oracle,
                       const std::string &    label)
  {
    verification::require_finite(state, label + " state");
    verification::require_finite(control, label + " control");
    verification::require_finite(returned_adjoint, label + " adjoint");
    verification::require_finite(returned_gradient, label + " gradient");

    const auto fresh = fresh_audit(control);
    const auto state_audit =
      verification::state_equation_residual(problem, fresh.state, control);
    const auto adjoint_audit =
      verification::adjoint_equation_residual(problem, fresh.adjoint, fresh.state);
    require(state_audit.normalized <= 1.0e-10,
            label + " state equation audit failed");
    require(adjoint_audit.normalized <= 1.0e-10,
            label + " adjoint equation audit failed");
    require(vector_difference(state, fresh.state) <= 1.0e-10,
            label + " returned state differs from fresh audit");
    require(vector_difference(returned_adjoint, fresh.adjoint) <= 1.0e-10,
            label + " returned adjoint differs from fresh audit");
    require(vector_difference(returned_gradient, fresh.gradient) <= 1.0e-10,
            label + " returned gradient differs from fresh audit");
    require(fresh.gradient.l2_norm() <= 1.1e-6,
            label + " fresh gradient audit exceeded the acceptance bound");

    const double oracle_distance = vector_difference(control, oracle.control);
    require(oracle_distance <= 2.0e-6,
            label + " control missed the dense oracle");
  }

  void
  run_minimal_problem_a_contract()
  {
    Instrumentation native_instrumentation;
    const auto artifact = external_dealii_step4_test::create_unique_artifact_root(
      std::filesystem::current_path() /
        "runs/external-dealii/step-4/minimal/problem-a",
      "consumer");
    external_dealii_step4_test::EvidenceGuard evidence(
      artifact,
      "minimal_problem_a_binding",
      {{"native", &native_instrumentation}});
    try
      {
        ProblemA native_problem(native_instrumentation);
        NativeReduced native_reduced(native_problem, native_instrumentation);
        const auto policy = external_dealii_step4::frozen_optimization_policy();

        ProblemA minimal_problem;
        Binding binding(minimal_problem);

        require(native_problem.state_dimension() ==
                  minimal_problem.state_dimension(),
                "minimal Problem A state dimensions differ from native reference");
        require(native_problem.control_dimension() ==
                  minimal_problem.control_dimension(),
                "minimal Problem A control dimensions differ from native reference");
        require(matrix_difference(native_problem.system_matrix(),
                                  minimal_problem.system_matrix()) == 0.0,
                "minimal Problem A matrix differs from native reference");
        require(vector_difference(native_problem.system_rhs(),
                                 minimal_problem.system_rhs()) == 0.0,
                "minimal Problem A RHS differs from native reference");

        Vector initial_control(native_problem.control_dimension());
        initial_control = 0.0;
        const auto native_result = NativeSolver(native_reduced, policy).solve(
          initial_control);

        Solver solver(binding.reduced(),
                      binding.metric(),
                      nmopt_parameters(policy));
        const auto minimal_result = solver.solve(
          nmopt::contract::PrimalBlockT<Backend>(binding.control_layout(),
                                                 {initial_control}));

        require(native_result.stopping_reason ==
                  external_dealii_step4::NativeOptimizationStoppingReason::
                    gradient_tolerance,
                "native Problem A reference did not stop by gradient tolerance");
        require(minimal_result.stopping_reason ==
                  nmopt::solvers::ReducedStoppingReason::gradient_tolerance,
                "minimal Problem A consumer did not stop by gradient tolerance");
        require(native_result.accepted_iteration_count ==
                  minimal_result.accepted_iterations,
                "minimal Problem A accepted-iteration count differs from native");
        require(native_result.line_search_trial_count ==
                  minimal_result.line_search_trial_count,
                "minimal Problem A trial count differs from native");

        const auto oracle = verification::optimum_oracle(native_problem);
        require(oracle.system_residual <= 1.0e-10 &&
                  oracle.stationarity_residual <= 1.0e-10,
                "minimal Problem A dense oracle audit failed");

        const auto &native_state = native_result.value.state;
        const auto &native_control = native_result.value.control;
        const auto &minimal_state = minimal_result.final_evaluation.state.block(0);
        const auto &minimal_control = minimal_result.control.block(0);
        const auto &minimal_adjoint =
          minimal_result.final_evaluation.adjoint.block(0);
        const auto &minimal_gradient =
          minimal_result.final_evaluation.reduced_derivative.block(0);

        require(vector_difference(native_state, minimal_state) <= 1.0e-10,
                "minimal Problem A state differs from native reference");
        require(vector_difference(native_control, minimal_control) <= 1.0e-10,
                "minimal Problem A control differs from native reference");
        require(scalar_matches(native_result.value.objective,
                               minimal_result.final_evaluation.objective_value),
                "minimal Problem A objective differs from native reference");
        require(native_result.derivative.reduced_derivative.l2_norm() <=
                  policy.gradient_tolerance &&
                  minimal_result.gradient_norm_history.back() <=
                    policy.gradient_tolerance,
                "minimal Problem A final gradient did not meet tolerance");

        check_final_solution(native_problem,
                             native_state,
                             native_control,
                             native_result.derivative.adjoint,
                             native_result.derivative.reduced_derivative,
                             oracle,
                             "native Problem A");
        check_final_solution(minimal_problem,
                             minimal_state,
                             minimal_control,
                             minimal_adjoint,
                             minimal_gradient,
                             oracle,
                             "minimal Problem A");

        const auto native_output = artifact / "native-solution.vtk";
        const auto minimal_output = artifact / "minimal-solution.vtk";
        native_problem.output_results(native_state, native_output);
        minimal_problem.output_results(minimal_state, minimal_output);
        require(output_payload(native_output) == output_payload(minimal_output),
                "minimal Problem A output differs from native reference");

        evidence.complete();
      }
    catch (...)
      {
        evidence.fail_current_exception();
        throw;
      }
  }

  void
  run_minimal_problem_a_executable(
    const std::filesystem::path &binary_directory)
  {
    const auto artifact = external_dealii_step4_test::create_unique_artifact_root(
      std::filesystem::current_path() /
        "runs/external-dealii/step-4/minimal/problem-a/executable",
      "consumer");
    const auto output = artifact / "solution.vtk";
    const auto logfile = artifact / "stdout.txt";
    const auto native_output = artifact / "native-solution.vtk";

    Instrumentation native_instrumentation;
    ProblemA native_problem(native_instrumentation);
    NativeReduced native_reduced(native_problem, native_instrumentation);
    const auto policy = external_dealii_step4::frozen_optimization_policy();
    Vector initial_control(native_problem.control_dimension());
    initial_control = 0.0;
    const auto native_result = NativeSolver(native_reduced, policy).solve(
      initial_control);
    const auto oracle = verification::optimum_oracle(native_problem);
    require(oracle.system_residual <= 1.0e-10 &&
              oracle.stationarity_residual <= 1.0e-10,
            "minimal Problem A executable dense oracle audit failed");
    check_final_solution(native_problem,
                         native_result.value.state,
                         native_result.value.control,
                         native_result.derivative.adjoint,
                         native_result.derivative.reduced_derivative,
                         oracle,
                         "executable native Problem A");
    native_problem.output_results(native_result.value.state, native_output);

    const auto application = binary_directory /
                             "nmopt_external_step4_minimal_problem_a";
    const auto command = shell_quote(application) + " " + shell_quote(output) +
                         " > " + shell_quote(logfile) + " 2>&1";
    require(std::system(command.c_str()) == 0,
            "minimal Problem A executable returned failure");

    const auto summary = read_consumer_summary(logfile);
    require(summary.stopping_reason == "gradient_tolerance",
            "minimal Problem A executable stopping reason differs from contract");
    require(summary.accepted_iterations ==
              native_result.accepted_iteration_count,
            "minimal Problem A executable iteration count differs from reference");
    require(summary.line_search_trials == native_result.line_search_trial_count,
            "minimal Problem A executable trial count differs from reference");
    require(scalar_matches(summary.objective, native_result.value.objective),
            "minimal Problem A executable objective differs from reference");
    require(scalar_matches(summary.gradient,
                           native_result.derivative.reduced_derivative.l2_norm()),
            "minimal Problem A executable gradient differs from reference");
    require(std::filesystem::exists(output),
            "minimal Problem A executable did not produce its output");
    require(output_payload(native_output) == output_payload(output),
            "minimal Problem A executable output differs from reference");
  }
} // namespace

int
main(const int argc, char **argv)
{
  try
    {
      const auto binary_directory =
        std::filesystem::absolute(argv[0]).parent_path();
      const std::vector<nmopt::test_support::Scenario> scenarios{
        {"minimal_problem_a_binding",
         "nmopt.external.tutorial_step_4.minimal_problem_a_binding_contract",
         {"dealii", "application", "external", "minimal", "problem_a"},
         360,
         run_minimal_problem_a_contract},
        {"minimal_problem_a_executable",
         "nmopt.external.tutorial_step_4.minimal_problem_a",
         {"dealii", "application", "external", "minimal", "problem_a"},
         360,
         [binary_directory] {
           run_minimal_problem_a_executable(binary_directory);
         }}};
      const auto result = nmopt::test_support::run_requested_scenarios(
        argc, argv, scenarios, std::cout);
      if (!result.listed)
        std::cout << "minimal Problem A binding contract passed: "
                  << result.executed << '\n';
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "minimal Problem A binding contract failed: "
                << exception.what() << '\n';
      return 1;
    }
}
