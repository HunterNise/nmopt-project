#include "../../apps/external-dealii/step-4/evaluation/native_problem_b_optimization.hpp"
#include "../../apps/external-dealii/step-4/minimal/problem_b_binding.hpp"
#include "../../apps/external-dealii/step-4/verification/problem_b_verification.hpp"

#include "../dealii/external_step4_evidence.hpp"
#include "../support/scenario_dispatch.hpp"

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
  using Backend        = nmopt::dealii_backend::SerialBackend;
  using Binding        = external_dealii_step4::minimal::ProblemBBinding;
  using Application    = Binding::Application;
  using Instrumentation = external_dealii_step4::Instrumentation;
  using NativeMetric   = external_dealii_step4::ProblemBMetric<2>;
  using NativeReduced  = external_dealii_step4::NativeProblemBReduced<
    2,
    Application>;
  using NativeSolver = external_dealii_step4::NativeProblemBArmijoSolver<
    2,
    Application>;
  using Parameters = nmopt::solvers::ReducedSolverParameters;
  using Problem    = Binding::Problem;
  using Solver     = nmopt::solvers::ReducedSearchSolverT<Backend>;
  using Vector     = Problem::Vector;
  using DenseOperators =
    external_dealii_step4::problem_b_verification::DenseOperators;
  namespace verification = external_dealii_step4::problem_b_verification;

  struct ConsumerSummary
  {
    bool        completed = false;
    bool        has_stopping_reason = false;
    std::string stopping_reason;
    bool        has_free_state_dimension = false;
    std::size_t free_state_dimension = 0;
    bool        has_control_dimension = false;
    std::size_t control_dimension = 0;
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
  nmopt_parameters(const external_dealii_step4::OptimizationPolicy &policy)
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

  bool
  scalar_matches(const double left, const double right)
  {
    return std::isfinite(left) && std::isfinite(right) &&
           std::abs(left - right) <=
             1.0e-12 + 1.0e-10 * std::max(std::abs(left), std::abs(right));
  }

  template <typename Matrix>
  double
  matrix_difference(const Matrix &left, const Matrix &right)
  {
    require(left.m() == right.m() && left.n() == right.n(),
            "minimal Problem B matrices have incompatible shapes");
    double squared_difference = 0.0;
    for (unsigned int row = 0; row < left.m(); ++row)
      for (unsigned int column = 0; column < left.n(); ++column)
        {
          const double difference =
            verification::matrix_value(left, row, column) -
            verification::matrix_value(right, row, column);
          squared_difference += difference * difference;
        }
    return std::sqrt(squared_difference);
  }

  double
  dense_primal_mass_norm(const dealii::FullMatrix<double> &mass,
                         const Vector &                     vector)
  {
    require(mass.m() == mass.n() &&
              vector.size() == static_cast<unsigned int>(mass.m()),
            "minimal Problem B dense primal mass norm dimensions are invalid");
    Vector mass_vector(vector.size());
    mass_vector = 0.0;
    for (unsigned int row = 0; row < mass.m(); ++row)
      for (unsigned int column = 0; column < mass.n(); ++column)
        mass_vector[row] += mass(row, column) * vector[column];
    const double squared_norm = vector * mass_vector;
    require(std::isfinite(squared_norm) && squared_norm >= 0.0,
            "minimal Problem B dense primal mass norm is invalid");
    return std::sqrt(squared_norm);
  }

  void
  compare_history(const std::vector<double> &left,
                  const std::vector<double> &right,
                  const char *const          name)
  {
    require(left.size() == right.size(),
            std::string("minimal Problem B ") + name +
              " histories have different sizes");
    for (std::size_t index = 0; index < left.size(); ++index)
      require(scalar_matches(left[index], right[index]),
              std::string("minimal Problem B ") + name +
                " histories differ at entry " + std::to_string(index));
  }

  std::string
  output_payload(const std::filesystem::path &filename)
  {
    std::ifstream input(filename, std::ios::binary);
    require(static_cast<bool>(input),
            "minimal Problem B output could not be opened");
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
  require_solution_audit(const verification::StateSolutionAudit &state_audit,
                         const verification::AdjointSolutionAudit &adjoint_audit,
                         const char *const path)
  {
    require(state_audit.full_equation.normalized <= 1.0e-10 &&
              state_audit.reduced_equation.normalized <= 1.0e-10 &&
              state_audit.reconstruction_error <= 1.0e-11 &&
              state_audit.boundary_error <= 1.0e-11,
            std::string(path) + " Problem B state audit failed");
    require(adjoint_audit.full_equation.normalized <= 1.0e-10 &&
              adjoint_audit.homogeneous_reconstruction_error <= 1.0e-11 &&
              adjoint_audit.boundary_error <= 1.0e-11,
            std::string(path) + " Problem B adjoint audit failed");
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
            "minimal Problem B executable output could not be opened");

    ConsumerSummary summary;
    std::string line;
    while (std::getline(input, line))
      {
        if (line == "minimal Problem B consumer completed")
          summary.completed = true;

        std::istringstream stream(line);
        std::string key;
        stream >> key;
        if (key == "stopping_reason")
          {
            require(static_cast<bool>(stream >> summary.stopping_reason),
                    "minimal Problem B executable stopping reason is malformed");
            summary.has_stopping_reason = true;
          }
        else if (key == "free_state_dimension")
          {
            require(static_cast<bool>(stream >> summary.free_state_dimension),
                    "minimal Problem B executable state dimension is malformed");
            summary.has_free_state_dimension = true;
          }
        else if (key == "control_dimension")
          {
            require(static_cast<bool>(stream >> summary.control_dimension),
                    "minimal Problem B executable control dimension is malformed");
            summary.has_control_dimension = true;
          }
        else if (key == "accepted_iterations")
          {
            require(static_cast<bool>(stream >> summary.accepted_iterations),
                    "minimal Problem B executable iteration count is malformed");
            summary.has_accepted_iterations = true;
          }
        else if (key == "line_search_trials")
          {
            require(static_cast<bool>(stream >> summary.line_search_trials),
                    "minimal Problem B executable trial count is malformed");
            summary.has_line_search_trials = true;
          }
        else if (key == "final_objective")
          {
            require(static_cast<bool>(stream >> summary.objective),
                    "minimal Problem B executable objective is malformed");
            summary.has_objective = true;
          }
        else if (key == "final_gradient_norm")
          {
            require(static_cast<bool>(stream >> summary.gradient),
                    "minimal Problem B executable gradient is malformed");
            summary.has_gradient = true;
          }
      }

    require(summary.completed,
            "minimal Problem B executable did not report completion");
    require(summary.has_stopping_reason && summary.has_free_state_dimension &&
              summary.has_control_dimension && summary.has_accepted_iterations &&
              summary.has_line_search_trials && summary.has_objective &&
              summary.has_gradient,
            "minimal Problem B executable report is incomplete");
    require(std::isfinite(summary.objective) && std::isfinite(summary.gradient),
            "minimal Problem B executable report is not finite");
    return summary;
  }

  void
  compare_dense_operators(const DenseOperators &left,
                          const DenseOperators &right)
  {
    require(matrix_difference(left.K, right.K) == 0.0,
            "minimal Problem B free stiffness differs from native reference");
    require(matrix_difference(left.B, right.B) == 0.0,
            "minimal Problem B coupling differs from native reference");
    require(matrix_difference(left.M, right.M) == 0.0,
            "minimal Problem B mass differs from native reference");
    require(matrix_difference(left.Q, right.Q) == 0.0,
            "minimal Problem B state mass differs from native reference");
    verification::require_vector_close(left.lifting,
                                       right.lifting,
                                       1.0e-11,
                                       1.0e-10,
                                       "minimal Problem B liftings differ");
    verification::require_vector_close(left.free_rhs,
                                       right.free_rhs,
                                       1.0e-11,
                                       1.0e-10,
                                       "minimal Problem B free RHS differs");
  }

  void
  run_minimal_problem_b_executable(
    const std::filesystem::path &binary_directory)
  {
    const auto artifact = external_dealii_step4_test::create_unique_artifact_root(
      std::filesystem::current_path() /
        "runs/external-dealii/step-4/minimal/problem-b/executable",
      "consumer");
    const auto output = artifact / "solution.vtk";
    const auto logfile = artifact / "stdout.txt";
    const auto native_output = artifact / "native-solution.vtk";

    Instrumentation native_instrumentation;
    Application native_application;
    native_application.prepare_for_external_use();
    Problem native_problem(native_application, &native_instrumentation);
    NativeMetric native_metric(native_problem.mass());
    NativeReduced native_reduced(native_problem, native_instrumentation);
    Vector initial_control(native_problem.control_dimension());
    initial_control = 0.0;
    const auto policy = external_dealii_step4::frozen_optimization_policy();
    const auto native_result = NativeSolver(native_reduced,
                                             native_metric,
                                             native_instrumentation,
                                             policy)
                                 .solve(initial_control);

    Application fresh_application;
    fresh_application.prepare_for_external_use();
    Problem fresh_problem(fresh_application);
    NativeReduced fresh_reduced(fresh_problem);
    const auto fresh_value =
      fresh_reduced.evaluate_value(native_result.value.control);
    const auto fresh_derivative = fresh_reduced.augment_derivative(fresh_value);
    const auto state_audit = verification::audit_state_solution(
      fresh_problem,
      fresh_application,
      fresh_value.state,
      fresh_value.full_state,
      fresh_value.control);
    const auto adjoint_audit = verification::audit_adjoint_solution(
      fresh_problem,
      fresh_application,
      fresh_derivative.adjoint,
      fresh_derivative.full_adjoint,
      fresh_derivative.state_derivative);
    require_solution_audit(state_audit, adjoint_audit, "executable native");
    const auto operators = verification::make_dense_operators(
      fresh_problem, fresh_application);
    const double stationarity = verification::dense_mass_norm(
      operators.M, fresh_derivative.reduced_derivative);
    require(stationarity <= 1.1e-6,
            "minimal Problem B executable reference failed stationarity audit");
    fresh_application.output_results(fresh_value.full_state, native_output);

    const auto application = binary_directory /
                             "nmopt_external_step4_minimal_problem_b";
    const auto command = shell_quote(application) + " " + shell_quote(output) +
                         " > " + shell_quote(logfile) + " 2>&1";
    require(std::system(command.c_str()) == 0,
            "minimal Problem B executable returned failure");

    const auto summary = read_consumer_summary(logfile);
    require(summary.stopping_reason == "gradient_tolerance",
            "minimal Problem B executable stopping reason differs from contract");
    require(summary.free_state_dimension == native_problem.state_dimension() &&
              summary.control_dimension == native_problem.control_dimension(),
            "minimal Problem B executable dimensions differ from reference");
    require(summary.accepted_iterations ==
              native_result.accepted_iteration_count,
            "minimal Problem B executable iteration count differs from reference");
    require(summary.line_search_trials == native_result.line_search_trial_count,
            "minimal Problem B executable trial count differs from reference");
    require(scalar_matches(summary.objective, fresh_value.objective),
            "minimal Problem B executable objective differs from reference");
    require(scalar_matches(summary.gradient, stationarity),
            "minimal Problem B executable gradient differs from dense audit");
    require(std::filesystem::exists(output),
            "minimal Problem B executable did not produce its output");
    require(output_payload(native_output) == output_payload(output),
            "minimal Problem B executable output differs from reference");
  }

  void
  run_minimal_problem_b_contract()
  {
    Instrumentation native_instrumentation;
    const auto artifact = external_dealii_step4_test::create_unique_artifact_root(
      std::filesystem::current_path() /
        "runs/external-dealii/step-4/minimal/problem-b",
      "consumer");
    external_dealii_step4_test::EvidenceGuard evidence(
      artifact,
      "minimal_problem_b_binding",
      {{"native", &native_instrumentation}});
    try
      {
        Application native_application;
        native_application.prepare_for_external_use();
        Problem native_problem(native_application, &native_instrumentation);
        NativeMetric native_metric(native_problem.mass());
        NativeReduced native_reduced(native_problem, native_instrumentation);

        Application minimal_application;
        minimal_application.prepare_for_external_use();
        Problem minimal_problem(minimal_application);
        Binding binding(minimal_problem);

        require(native_problem.state_dimension() == 225 &&
                  minimal_problem.state_dimension() == 225,
                "minimal Problem B free state dimension is wrong");
        require(native_problem.control_dimension() == 289 &&
                  minimal_problem.control_dimension() == 289,
                "minimal Problem B control dimension is wrong");
        require(native_problem.coordinates().boundary_dimension() == 64 &&
                  minimal_problem.coordinates().boundary_dimension() == 64,
                "minimal Problem B boundary dimension is wrong");
        require(native_problem.coordinates().free_indices() ==
                  minimal_problem.coordinates().free_indices(),
                "minimal Problem B free index maps differ");
        require(matrix_difference(native_application.system_matrix_view(),
                                  minimal_application.system_matrix_view()) == 0.0,
                "minimal Problem B full stiffness differs from native reference");
        verification::require_vector_close(
          native_application.system_rhs_view(),
          minimal_application.system_rhs_view(),
          1.0e-11,
          1.0e-10,
          "minimal Problem B full RHS differs");

        const auto native_operators =
          verification::make_dense_operators(native_problem,
                                              native_application);
        const auto minimal_operators =
          verification::make_dense_operators(minimal_problem,
                                              minimal_application);
        compare_dense_operators(native_operators, minimal_operators);

        Vector initial_control(native_problem.control_dimension());
        initial_control = 0.0;
        const auto policy = external_dealii_step4::frozen_optimization_policy();
        const auto native_result = NativeSolver(native_reduced,
                                                native_metric,
                                                native_instrumentation,
                                                policy)
                                      .solve(initial_control);

        Solver solver(binding.reduced(),
                      binding.metric(),
                      nmopt_parameters(policy));
        const auto minimal_result = solver.solve(
          nmopt::contract::PrimalBlockT<Backend>(
            binding.control_layout(), std::vector<Backend::Vector>{
                                      initial_control}));

        require(native_result.stopping_reason ==
                  external_dealii_step4::NativeProblemBOptimizationStoppingReason::
                    gradient_tolerance,
                "native Problem B reference did not stop by gradient tolerance");
        require(minimal_result.stopping_reason ==
                  nmopt::solvers::ReducedStoppingReason::gradient_tolerance,
                "minimal Problem B consumer did not stop by gradient tolerance");
        require(native_result.accepted_iteration_count ==
                  minimal_result.accepted_iterations,
                "minimal Problem B accepted-iteration count differs from native");
        require(native_result.line_search_trial_count ==
                  minimal_result.line_search_trial_count,
                "minimal Problem B trial count differs from native");
        compare_history(native_result.objective_history,
                        minimal_result.objective_history,
                        "objective");
        compare_history(native_result.gradient_norm_history,
                        minimal_result.gradient_norm_history,
                        "gradient");

        const auto &native_control = native_result.value.control;
        const auto &native_state = native_result.value.state;
        const auto &minimal_control = minimal_result.control.block(0);
        const auto &minimal_state =
          minimal_result.final_evaluation.state.block(0);
        const auto &minimal_adjoint =
          minimal_result.final_evaluation.adjoint.block(0);
        const auto &minimal_gradient =
          minimal_result.final_evaluation.reduced_derivative.block(0);
        const Vector minimal_full_state =
          minimal_problem.coordinates().reconstruct(minimal_state);
        const Vector minimal_full_adjoint =
          minimal_problem.coordinates().embed_free(minimal_adjoint);

        verification::require_vector_close(native_control,
                                           minimal_control,
                                           1.0e-11,
                                           1.0e-10,
                                           "minimal Problem B controls differ");
        verification::require_vector_close(native_state,
                                           minimal_state,
                                           1.0e-11,
                                           1.0e-10,
                                           "minimal Problem B states differ");
        verification::require_vector_close(native_result.value.full_state,
                                           minimal_full_state,
                                           1.0e-11,
                                           1.0e-10,
                                           "minimal Problem B full states differ");
        verification::require_vector_close(
          native_result.derivative.adjoint,
          minimal_adjoint,
          1.0e-11,
          1.0e-10,
          "minimal Problem B adjoints differ");
        verification::require_vector_close(
          native_result.derivative.full_adjoint,
          minimal_full_adjoint,
          1.0e-11,
          1.0e-10,
          "minimal Problem B full adjoints differ");
        verification::require_vector_close(
          native_result.derivative.reduced_derivative,
          minimal_gradient,
          1.0e-11,
          1.0e-10,
          "minimal Problem B reduced gradients differ");
        require(scalar_matches(native_result.value.objective,
                               minimal_result.final_evaluation.objective_value),
                "minimal Problem B objectives differ");

        std::filesystem::create_directories(artifact / "native");
        std::filesystem::create_directories(artifact / "minimal");
        native_application.output_results(native_result.value.full_state,
                                          artifact / "native" / "solution.vtk");
        minimal_application.output_results(minimal_full_state,
                                           artifact / "minimal" / "solution.vtk");
        require(output_payload(artifact / "native" / "solution.vtk") ==
                  output_payload(artifact / "minimal" / "solution.vtk"),
                "minimal Problem B output differs from native reference");

        Application native_verification_application;
        native_verification_application.prepare_for_external_use();
        Problem native_verification_problem(native_verification_application);
        NativeReduced native_verification_reduced(native_verification_problem);
        const auto fresh_native_value =
          native_verification_reduced.evaluate_value(native_control);
        const auto fresh_native_derivative =
          native_verification_reduced.augment_derivative(fresh_native_value);

        Application public_verification_application;
        public_verification_application.prepare_for_external_use();
        Problem public_verification_problem(public_verification_application);
        Binding public_verification_binding(public_verification_problem);
        const auto fresh_public_value =
          public_verification_binding.reduced().evaluate_value(
            nmopt::contract::PrimalBlockT<Backend>(
              public_verification_binding.control_layout(),
              std::vector<Backend::Vector>{minimal_control}));
        const auto fresh_public_derivative =
          public_verification_binding.reduced().augment_derivative(
            fresh_public_value);
        const auto &fresh_public_state = fresh_public_value.state.block(0);
        const auto &fresh_public_adjoint =
          fresh_public_derivative.adjoint.block(0);
        const auto &fresh_public_gradient =
          fresh_public_derivative.reduced_derivative.block(0);
        const Vector fresh_public_full_state =
          public_verification_problem.coordinates().reconstruct(
            fresh_public_state);
        const Vector fresh_public_full_adjoint =
          public_verification_problem.coordinates().embed_free(
            fresh_public_adjoint);
        const auto public_objective_derivative =
          public_verification_problem.objective_derivative(
            fresh_public_state, minimal_control);

        verification::require_vector_close(
          fresh_native_value.state,
          native_state,
          1.0e-11,
          1.0e-10,
          "fresh native Problem B state differs from returned state");
        verification::require_vector_close(
          fresh_native_value.full_state,
          native_result.value.full_state,
          1.0e-11,
          1.0e-10,
          "fresh native Problem B full state differs from returned state");
        verification::require_vector_close(
          fresh_native_derivative.adjoint,
          native_result.derivative.adjoint,
          1.0e-11,
          1.0e-10,
          "fresh native Problem B adjoint differs from returned adjoint");
        verification::require_vector_close(
          fresh_native_derivative.full_adjoint,
          native_result.derivative.full_adjoint,
          1.0e-11,
          1.0e-10,
          "fresh native Problem B full adjoint differs from returned adjoint");
        verification::require_vector_close(
          fresh_native_derivative.reduced_derivative,
          native_result.derivative.reduced_derivative,
          1.0e-11,
          1.0e-10,
          "fresh native Problem B gradient differs from returned gradient");
        require(scalar_matches(fresh_native_value.objective,
                               native_result.value.objective),
                "fresh native Problem B objective differs from returned value");

        verification::require_vector_close(
          fresh_public_state,
          minimal_state,
          1.0e-11,
          1.0e-10,
          "fresh public Problem B state differs from returned state");
        verification::require_vector_close(
          fresh_public_full_state,
          minimal_full_state,
          1.0e-11,
          1.0e-10,
          "fresh public Problem B full state differs from returned state");
        verification::require_vector_close(
          fresh_public_adjoint,
          minimal_adjoint,
          1.0e-11,
          1.0e-10,
          "fresh public Problem B adjoint differs from returned adjoint");
        verification::require_vector_close(
          fresh_public_full_adjoint,
          minimal_full_adjoint,
          1.0e-11,
          1.0e-10,
          "fresh public Problem B full adjoint differs from returned adjoint");
        verification::require_vector_close(
          fresh_public_gradient,
          minimal_gradient,
          1.0e-11,
          1.0e-10,
          "fresh public Problem B gradient differs from returned gradient");
        require(scalar_matches(fresh_public_value.objective_value,
                               minimal_result.final_evaluation.objective_value),
                "fresh public Problem B objective differs from returned value");

        const auto native_state_audit = verification::audit_state_solution(
          native_verification_problem,
          native_verification_application,
          fresh_native_value.state,
          fresh_native_value.full_state,
          fresh_native_value.control);
        const auto native_adjoint_audit = verification::audit_adjoint_solution(
          native_verification_problem,
          native_verification_application,
          fresh_native_derivative.adjoint,
          fresh_native_derivative.full_adjoint,
          fresh_native_derivative.state_derivative);
        const auto public_state_audit = verification::audit_state_solution(
          public_verification_problem,
          public_verification_application,
          fresh_public_state,
          fresh_public_full_state,
          minimal_control);
        const auto public_adjoint_audit = verification::audit_adjoint_solution(
          public_verification_problem,
          public_verification_application,
          fresh_public_adjoint,
          fresh_public_full_adjoint,
          public_objective_derivative.state);
        require_solution_audit(native_state_audit,
                               native_adjoint_audit,
                               "native");
        require_solution_audit(public_state_audit,
                               public_adjoint_audit,
                               "minimal");

        const auto operators = verification::make_dense_operators(
          native_verification_problem, native_verification_application);
        const auto public_operators = verification::make_dense_operators(
          public_verification_problem, public_verification_application);
        compare_dense_operators(operators, public_operators);

        const double native_stationarity = verification::dense_mass_norm(
          operators.M, fresh_native_derivative.reduced_derivative);
        const double public_stationarity = verification::dense_mass_norm(
          operators.M, fresh_public_gradient);
        require(native_stationarity <= 1.1e-6 &&
                  public_stationarity <= 1.1e-6,
                "minimal Problem B final stationarity audit failed");

        const auto oracle = verification::dense_kkt_oracle(operators);
        require(oracle.residuals.state_stationarity <= 1.0e-10 &&
                  oracle.residuals.control_stationarity <= 1.0e-10 &&
                  oracle.residuals.feasibility <= 1.0e-10,
                "minimal Problem B dense KKT oracle failed");
        const auto oracle_value =
          native_verification_reduced.evaluate_value(oracle.control);
        const auto oracle_derivative =
          native_verification_reduced.augment_derivative(oracle_value);
        require(verification::dense_mass_norm(
                  operators.M, oracle_derivative.reduced_derivative) <= 1.0e-8,
                "minimal Problem B oracle reduced gradient is not zero");

        Vector native_control_difference = native_control;
        native_control_difference.add(-1.0, oracle.control);
        Vector public_control_difference = minimal_control;
        public_control_difference.add(-1.0, oracle.control);
        require(dense_primal_mass_norm(operators.M, native_control_difference) <=
                  2.0e-6 &&
                  dense_primal_mass_norm(operators.M, public_control_difference) <=
                    2.0e-6,
                "minimal Problem B control missed the dense KKT oracle");

        evidence.complete();
      }
    catch (...)
      {
        evidence.fail_current_exception();
        throw;
      }
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
        {"minimal_problem_b_binding",
         "nmopt.external.tutorial_step_4.minimal_problem_b_binding_contract",
         {"dealii", "application", "external", "minimal", "problem_b"},
         600,
         run_minimal_problem_b_contract},
        {"minimal_problem_b_executable",
         "nmopt.external.tutorial_step_4.minimal_problem_b",
         {"dealii", "application", "external", "minimal", "problem_b"},
         600,
         [binary_directory] {
           run_minimal_problem_b_executable(binary_directory);
         }}};
      const auto result = nmopt::test_support::run_requested_scenarios(
        argc, argv, scenarios, std::cout);
      if (!result.listed)
        std::cout << "minimal Problem B binding contract passed: "
                  << result.executed << '\n';
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "minimal Problem B binding contract failed: "
                << exception.what() << '\n';
      return 1;
    }
}
