#define STEP4_NO_MAIN
#include "../../apps/external-dealii/step-4/source/adapted/step-4.cc"
#undef STEP4_NO_MAIN

#include "../../apps/external-dealii/step-4/evaluation/native_problem_b_optimization.hpp"
#include "../../apps/external-dealii/step-4/integration/nmopt_problem_b_binding.hpp"
#include "../../apps/external-dealii/step-4/verification/problem_b_verification.hpp"

#include "../dealii/external_step4_evidence.hpp"
#include "../support/scenario_dispatch.hpp"

#include "nmopt/solvers/reduced_gradient.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
  using Backend = nmopt::dealii_backend::SerialBackend;
  using Binding =
    external_dealii_step4::ProblemBNmoptBinding<2, Step4<2>>;
  using Instrumentation = external_dealii_step4::Instrumentation;
  using NativeProblem = external_dealii_step4::ProblemB<2, Step4<2>>;
  using NativeMetric = external_dealii_step4::ProblemBMetric<2>;
  using NativeReduced =
    external_dealii_step4::NativeProblemBReduced<2, Step4<2>>;
  using NativeSolver =
    external_dealii_step4::NativeProblemBArmijoSolver<2, Step4<2>>;
  using NativeResult =
    external_dealii_step4::NativeProblemBOptimizationResult<2, Step4<2>>;
  using OptimizationPolicy = external_dealii_step4::OptimizationPolicy;
  using Vector = NativeProblem::Vector;
  using NmoptSolver = nmopt::solvers::ReducedSearchSolverT<Backend>;
  using NmoptResult = nmopt::solvers::ReducedSolverResultT<Backend>;

  struct TrialRecord
  {
    std::size_t trial = 0;
    std::size_t iteration = 0;
    double      step_length = 0.0;
    double      objective_value = 0.0;
    double      actual_slope = 0.0;
    double      sufficient_decrease_bound = 0.0;
    bool        objective_finite = false;
    bool        slope_negative = false;
    bool        accepted = false;
  };

  struct AcceptedRecord
  {
    std::size_t iteration = 0;
    double      objective_before = 0.0;
    double      objective_after = 0.0;
    double      objective_change = 0.0;
    double      requested_step_length = 0.0;
    double      actual_step_norm = 0.0;
    double      actual_slope = 0.0;
    double      gradient_norm = 0.0;
    std::size_t trial_count = 0;
  };

  struct ScalarTrace
  {
    std::vector<double>      objective_history;
    std::vector<double>      gradient_norms;
    std::vector<TrialRecord> trials;
    std::vector<AcceptedRecord> accepted;
    std::string              stopping_reason;
    double                   final_objective = 0.0;
    double                   final_gradient_norm = 0.0;
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

  bool
  scalar_matches(const double left, const double right)
  {
    return std::isfinite(left) && std::isfinite(right) &&
           std::abs(left - right) <=
             1.0e-12 + 1.0e-11 * std::max(std::abs(left), std::abs(right));
  }

  std::filesystem::path
  find_repository_root()
  {
    auto directory = std::filesystem::current_path();
    while (true)
      {
        if (std::filesystem::exists(
              directory / "apps/external-dealii/step-4/source/upstream/step-4.cc"))
          return directory;

        const auto parent = directory.parent_path();
        if (parent == directory)
          break;
        directory = parent;
      }

    throw std::runtime_error("could not locate the repository root");
  }

  std::filesystem::path
  create_artifact_root()
  {
    const auto root = external_dealii_step4_test::create_unique_artifact_root(
      find_repository_root() / "runs/external-dealii/step-4/problem-b/optimization",
      "paired");
    std::filesystem::create_directories(root / "native");
    std::filesystem::create_directories(root / "nmopt");
    std::filesystem::create_directories(root / "comparison");
    return root;
  }

  nmopt::solvers::ReducedSolverParameters
  nmopt_parameters(const OptimizationPolicy &policy)
  {
    nmopt::solvers::ReducedSolverParameters parameters;
    parameters.maximum_iterations = policy.maximum_iterations;
    parameters.maximum_line_search_trials =
      policy.maximum_line_search_trials;
    parameters.gradient_tolerance = policy.gradient_tolerance;
    parameters.stopping_criterion =
      nmopt::solvers::ReducedStoppingCriterion::gradient_norm;
    parameters.relative_gradient_tolerance =
      policy.relative_gradient_tolerance;
    parameters.objective_change_tolerance =
      policy.objective_change_tolerance;
    parameters.step_tolerance = policy.step_tolerance;
    parameters.objective_target = policy.objective_target;
    parameters.initial_step_length = policy.initial_step_length;
    parameters.minimum_step_length = policy.minimum_step_length;
    parameters.armijo_fraction = policy.armijo_fraction;
    parameters.backtracking_factor = policy.backtracking_factor;
    return parameters;
  }

  ScalarTrace
  native_trace(const NativeResult &result)
  {
    ScalarTrace trace;
    trace.objective_history = result.objective_history;
    trace.gradient_norms = result.gradient_norm_history;
    trace.stopping_reason =
      external_dealii_step4::native_problem_b_optimization_stopping_reason_name(
        result.stopping_reason);
    trace.final_objective = result.value.objective;
    trace.final_gradient_norm = result.gradient_norm_history.back();
    trace.trials.reserve(result.trial_records.size());
    for (const auto &trial : result.trial_records)
      trace.trials.push_back({trial.trial,
                              trial.iteration,
                              trial.step_length,
                              trial.objective_value,
                              trial.actual_slope,
                              trial.sufficient_decrease_bound,
                              trial.objective_finite,
                              trial.slope_negative,
                              trial.accepted});
    trace.accepted.reserve(result.accepted_iterations.size());
    for (const auto &iteration : result.accepted_iterations)
      trace.accepted.push_back({iteration.iteration,
                                iteration.objective_before,
                                iteration.objective_after,
                                iteration.objective_change,
                                iteration.requested_step_length,
                                iteration.actual_step_norm,
                                iteration.actual_slope,
                                iteration.gradient_norm,
                                iteration.trial_count});
    return trace;
  }

  ScalarTrace
  nmopt_trace(const NmoptResult &result)
  {
    ScalarTrace trace;
    trace.objective_history = result.objective_history;
    trace.gradient_norms = result.gradient_norm_history;
    trace.stopping_reason =
      nmopt::solvers::reduced_stopping_reason_name(result.stopping_reason);
    trace.final_objective = result.final_evaluation.objective_value;
    trace.final_gradient_norm = result.gradient_norm_history.back();
    trace.trials.reserve(result.line_search_trials.size());
    for (const auto &trial : result.line_search_trials)
      trace.trials.push_back({trial.trial,
                              trial.iteration,
                              trial.step_length,
                              trial.objective_value,
                              trial.actual_slope,
                              trial.sufficient_decrease_bound,
                              trial.objective_finite,
                              trial.slope_negative,
                              trial.accepted});
    trace.accepted.reserve(result.iteration_records.size());
    for (const auto &record : result.iteration_records)
      {
        const auto &common = record.common;
        trace.accepted.push_back({common.iteration,
                                  common.objective_before,
                                  common.objective_after,
                                  common.objective_change,
                                  common.requested_step_parameter,
                                  common.actual_step_norm,
                                  common.actual_descent_pairing,
                                  common.absolute_stationarity,
                                  common.trial_count});
      }
    return trace;
  }

  void
  write_trace(const std::filesystem::path &root,
              const char *const            path_name,
              const ScalarTrace &          trace)
  {
    std::ofstream output(root / "trace.csv");
    require(static_cast<bool>(output),
            std::string("could not open the ") + path_name +
              " Problem B optimization trace");
    output << std::setprecision(std::numeric_limits<double>::max_digits10)
           << "path " << path_name << '\n'
           << "stopping_reason " << trace.stopping_reason << '\n'
           << "accepted_iterations " << trace.accepted.size() << '\n'
           << "line_search_trials " << trace.trials.size() << '\n'
           << "final_objective " << trace.final_objective << '\n'
           << "final_gradient_norm " << trace.final_gradient_norm << '\n'
           << "objective_history";
    for (const double value : trace.objective_history)
      output << ' ' << value;
    output << '\n'
           << "record,iteration,trial,step_length,objective,actual_slope,"
              "armijo_bound,objective_finite,slope_negative,accepted,"
              "objective_before,objective_after,objective_change,"
              "actual_step_norm,gradient_norm\n";
    for (const auto &trial : trace.trials)
      output << "trial," << trial.iteration << ',' << trial.trial << ','
             << trial.step_length << ',' << trial.objective_value << ','
             << trial.actual_slope << ',' << trial.sufficient_decrease_bound
             << ',' << trial.objective_finite << ',' << trial.slope_negative
             << ',' << trial.accepted << ",,,,,\n";
    for (const auto &iteration : trace.accepted)
      output << "accepted," << iteration.iteration << ",,"
             << iteration.requested_step_length << ','
             << iteration.objective_after << ',' << iteration.actual_slope
             << ",,,1," << iteration.objective_before << ','
             << iteration.objective_after << ',' << iteration.objective_change
             << ',' << iteration.actual_step_norm << ','
             << iteration.gradient_norm << '\n';
  }

  void
  compare_traces(const ScalarTrace &left,
                 const ScalarTrace &right,
                 std::string &       first_difference)
  {
    const auto note = [&](const std::string &message) {
      if (first_difference.empty())
        first_difference = message;
    };
    if (left.stopping_reason != right.stopping_reason)
      note("stopping reason differs");

    if (left.objective_history.size() != right.objective_history.size())
      note("objective-history size differs");
    else
      for (std::size_t index = 0; index < left.objective_history.size(); ++index)
        if (!scalar_matches(left.objective_history[index],
                            right.objective_history[index]))
          note("objective first differs at history entry " +
               std::to_string(index));

    if (left.gradient_norms.size() != right.gradient_norms.size())
      note("gradient-history size differs");
    else
      for (std::size_t index = 0; index < left.gradient_norms.size(); ++index)
        if (!scalar_matches(left.gradient_norms[index],
                            right.gradient_norms[index]))
          note("gradient norm first differs at check " +
               std::to_string(index));

    if (left.trials.size() != right.trials.size())
      note("trial-record size differs");
    const auto trial_count = std::min(left.trials.size(), right.trials.size());
    for (std::size_t index = 0; index < trial_count; ++index)
      {
        const auto &a = left.trials[index];
        const auto &b = right.trials[index];
        if (a.trial != b.trial || a.iteration != b.iteration)
          note("trial index first differs at record " + std::to_string(index));
        if (!scalar_matches(a.step_length, b.step_length) ||
            !scalar_matches(a.objective_value, b.objective_value) ||
            !scalar_matches(a.actual_slope, b.actual_slope) ||
            !scalar_matches(a.sufficient_decrease_bound,
                            b.sufficient_decrease_bound))
          note("trial scalar first differs at record " +
               std::to_string(index));
        if (a.objective_finite != b.objective_finite ||
            a.slope_negative != b.slope_negative || a.accepted != b.accepted)
          note("trial acceptance first differs at record " +
               std::to_string(index));
      }

    if (left.accepted.size() != right.accepted.size())
      note("accepted-record size differs");
    const auto accepted_count =
      std::min(left.accepted.size(), right.accepted.size());
    for (std::size_t index = 0; index < accepted_count; ++index)
      {
        const auto &a = left.accepted[index];
        const auto &b = right.accepted[index];
        if (a.iteration != b.iteration || a.trial_count != b.trial_count)
          note("accepted index first differs at record " +
               std::to_string(index));
        if (!scalar_matches(a.objective_before, b.objective_before) ||
            !scalar_matches(a.objective_after, b.objective_after) ||
            !scalar_matches(a.objective_change, b.objective_change) ||
            !scalar_matches(a.requested_step_length,
                            b.requested_step_length) ||
            !scalar_matches(a.actual_step_norm, b.actual_step_norm) ||
            !scalar_matches(a.actual_slope, b.actual_slope) ||
            !scalar_matches(a.gradient_norm, b.gradient_norm))
          note("accepted scalar first differs at record " +
               std::to_string(index));
      }

    if (!scalar_matches(left.final_objective, right.final_objective))
      note("final objective differs");
    if (!scalar_matches(left.final_gradient_norm, right.final_gradient_norm))
      note("final gradient norm differs");
  }

  void
  require_native_result_shape(const NativeResult &result)
  {
    require(result.stopping_reason ==
              external_dealii_step4::NativeProblemBOptimizationStoppingReason::
                gradient_tolerance,
            "native Problem B optimization did not stop by gradient tolerance");
    require(result.accepted_iteration_count > 0,
            "native Problem B optimization accepted no iterations");
    require(result.objective_history.size() ==
                result.accepted_iteration_count + 1,
            "native objective history has the wrong size");
    require(result.gradient_norm_history.size() ==
                result.accepted_iteration_count + 1,
            "native gradient history has the wrong size");
  }

  void
  require_nmopt_result_shape(const NmoptResult &result)
  {
    require(result.stopping_reason ==
              nmopt::solvers::ReducedStoppingReason::gradient_tolerance,
            "nmopt Problem B optimization did not stop by gradient tolerance");
    require(result.accepted_iterations > 0,
            "nmopt Problem B optimization accepted no iterations");
    require(result.objective_history.size() == result.accepted_iterations + 1,
            "nmopt objective history has the wrong size");
    require(result.gradient_norm_history.size() ==
                result.accepted_iterations + 1,
            "nmopt gradient history has the wrong size");
  }

  void
  require_native_schedule(const NativeResult &   result,
                          const Instrumentation &instrumentation)
  {
    const auto accepted = result.accepted_iteration_count;
    const auto trials = result.line_search_trial_count;
    require(result.trial_records.size() == trials &&
              result.accepted_iterations.size() == accepted &&
              result.metric_inverse_evidence.size() == accepted + 1,
            "native Problem B optimization records have the wrong sizes");
    require(instrumentation.state_solve_calls == 1 + trials &&
              instrumentation.adjoint_solve_calls == 1 + accepted &&
              instrumentation.objective_calls == 1 + trials &&
              instrumentation.objective_derivative_calls == 1 + accepted &&
              instrumentation.value_evaluations == 1 + trials &&
              instrumentation.derivative_augmentations == 1 + accepted &&
              instrumentation.control_vjp_calls == 1 + accepted,
            "native Problem B optimization staged schedule is inconsistent");
    require(instrumentation.residual_calls == 0 &&
              instrumentation.residual_jvp_calls == 0 &&
              instrumentation.residual_vjp_calls == 0 &&
              instrumentation.explicit_matrix_vmult_calls == 0 &&
              instrumentation.explicit_matrix_tvmult_calls == 0,
            "native Problem B optimization used an unexpected callback");
    require(instrumentation.metric_inverse_apply_calls == 1 + accepted &&
              instrumentation.metric_apply_calls == 1 + 2 * accepted,
            "native Problem B optimization metric schedule is inconsistent");
    require(instrumentation.solve_failures == 0 &&
              instrumentation.solve_records.size() == 2 + trials + accepted,
            "native Problem B optimization solve evidence is inconsistent");
  }

  void
  require_nmopt_schedule(const NmoptResult &   result,
                         const Instrumentation &instrumentation)
  {
    const auto accepted = result.accepted_iterations;
    const auto trials = result.line_search_trial_count;
    require(result.line_search_trials.size() == trials &&
              result.iteration_records.size() == accepted &&
              result.state_solve_count == 1 + trials &&
              result.adjoint_solve_count == 1 + accepted,
            "nmopt Problem B optimization records have the wrong sizes");
    require(result.metric_solve_count == 1 + accepted &&
              result.hessian_action_count == 0,
            "nmopt Problem B optimization solver work is inconsistent");
    require(instrumentation.state_solve_calls == 1 + trials &&
              instrumentation.adjoint_solve_calls == 1 + accepted &&
              instrumentation.objective_calls == 1 + trials &&
              instrumentation.objective_derivative_calls == 1 + accepted &&
              instrumentation.value_evaluations == 0 &&
              instrumentation.derivative_augmentations == 0 &&
              instrumentation.control_vjp_calls == 0 &&
              instrumentation.residual_vjp_calls == 1 + accepted,
            "nmopt Problem B optimization staged schedule is inconsistent");
    require(instrumentation.residual_calls == 0 &&
              instrumentation.residual_jvp_calls == 0 &&
              instrumentation.explicit_matrix_vmult_calls == 0 &&
              instrumentation.explicit_matrix_tvmult_calls == 0,
            "nmopt Problem B optimization used an unexpected callback");
    require(instrumentation.metric_inverse_apply_calls == 1 + accepted &&
              instrumentation.metric_apply_calls == 1 + 2 * accepted,
            "nmopt Problem B optimization metric schedule is inconsistent");
    require(instrumentation.solve_failures == 0 &&
              instrumentation.solve_records.size() == 2 + trials + accepted,
            "nmopt Problem B optimization solve evidence is inconsistent");
  }

  void
  write_comparison_summary(const std::filesystem::path &root,
                           const ScalarTrace &          native,
                           const ScalarTrace &          nmopt,
                           const std::string &          first_difference,
                           const Instrumentation &     native_instrumentation,
                           const Instrumentation &     nmopt_instrumentation)
  {
    std::ofstream summary(root / "summary.txt");
    require(static_cast<bool>(summary),
            "could not open the Problem B optimization comparison summary");
    summary << std::setprecision(std::numeric_limits<double>::max_digits10)
            << "first_divergence "
            << (first_difference.empty() ? "none" : first_difference) << '\n'
            << "native_stopping_reason " << native.stopping_reason << '\n'
            << "nmopt_stopping_reason " << nmopt.stopping_reason << '\n'
            << "native_accepted_iterations " << native.accepted.size() << '\n'
            << "nmopt_accepted_iterations " << nmopt.accepted.size() << '\n'
            << "native_line_search_trials " << native.trials.size() << '\n'
            << "nmopt_line_search_trials " << nmopt.trials.size() << '\n'
            << "native_state_solve_calls "
            << native_instrumentation.state_solve_calls << '\n'
            << "nmopt_state_solve_calls "
            << nmopt_instrumentation.state_solve_calls << '\n'
            << "native_adjoint_solve_calls "
            << native_instrumentation.adjoint_solve_calls << '\n'
            << "nmopt_adjoint_solve_calls "
            << nmopt_instrumentation.adjoint_solve_calls << '\n'
            << "native_control_vjp "
            << native_instrumentation.control_vjp_calls << '\n'
            << "nmopt_residual_vjp "
            << nmopt_instrumentation.residual_vjp_calls << '\n'
            << "native_metric_inverse_apply "
            << native_instrumentation.metric_inverse_apply_calls << '\n'
            << "nmopt_metric_inverse_apply "
            << nmopt_instrumentation.metric_inverse_apply_calls << '\n'
            << "native_metric_apply "
            << native_instrumentation.metric_apply_calls << '\n'
            << "nmopt_metric_apply "
            << nmopt_instrumentation.metric_apply_calls << '\n';
  }

  std::string
  read_vtk_payload(const std::filesystem::path &filename)
  {
    std::ifstream input(filename, std::ios::binary);
    require(static_cast<bool>(input), "could not open VTK comparison file");
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
  run_problem_b_nmopt_matched_optimization()
  {
    Instrumentation native_instrumentation;
    Instrumentation nmopt_instrumentation;
    Instrumentation native_verification_instrumentation;
    Instrumentation nmopt_verification_instrumentation;
    const auto artifact_root = create_artifact_root();
    external_dealii_step4_test::EvidenceGuard evidence(
      artifact_root,
      "problem_b_nmopt_matched_optimization",
      {{"native", &native_instrumentation},
       {"nmopt", &nmopt_instrumentation},
       {"native-verification", &native_verification_instrumentation},
       {"nmopt-verification", &nmopt_verification_instrumentation}});
    try
      {
        Step4<2> native_tutorial;
        native_tutorial.prepare_for_external_use();
        NativeProblem native_problem(native_tutorial);
        NativeMetric native_metric(native_problem.mass());
        NativeReduced native_reduced(native_problem, native_instrumentation);
        Binding nmopt_binding(nmopt_instrumentation);

        require(native_problem.state_dimension() ==
                  nmopt_binding.problem().state_dimension() &&
                  native_problem.control_dimension() ==
                    nmopt_binding.problem().control_dimension(),
                "paired Problem B optimization dimensions differ");

        Vector initial_control(native_problem.control_dimension());
        initial_control = 0.0;
        const auto policy = external_dealii_step4::frozen_optimization_policy();

        NativeSolver native_solver(
          native_reduced, native_metric, native_instrumentation, policy);
        const auto native_result = native_solver.solve(initial_control);
        const auto native_trace_value = native_trace(native_result);
        write_trace(artifact_root / "native", "native", native_trace_value);

        NmoptSolver nmopt_solver(nmopt_binding.reduced(),
                                 nmopt_binding.metric(),
                                 nmopt_parameters(policy));
        const auto nmopt_result = nmopt_solver.solve(
          nmopt::contract::PrimalBlockT<Backend>(
            nmopt_binding.control_layout(), {initial_control}));
        const auto nmopt_trace_value = nmopt_trace(nmopt_result);
        write_trace(artifact_root / "nmopt", "nmopt", nmopt_trace_value);

        require_native_result_shape(native_result);
        require_nmopt_result_shape(nmopt_result);
        require_native_schedule(native_result, native_instrumentation);
        require_nmopt_schedule(nmopt_result, nmopt_instrumentation);

        std::string first_difference;
        compare_traces(native_trace_value, nmopt_trace_value, first_difference);
        write_comparison_summary(artifact_root / "comparison",
                                 native_trace_value,
                                 nmopt_trace_value,
                                 first_difference,
                                 native_instrumentation,
                                 nmopt_instrumentation);
        require(first_difference.empty(), first_difference);

        const auto &nmopt_control = nmopt_result.control.block(0);
        const auto &nmopt_state = nmopt_result.final_evaluation.state.block(0);
        const auto &nmopt_adjoint =
          nmopt_result.final_evaluation.adjoint.block(0);
        const auto &nmopt_gradient =
          nmopt_result.final_evaluation.reduced_derivative.block(0);
        const Vector nmopt_full_state =
          nmopt_binding.problem().coordinates().reconstruct(nmopt_state);
        const Vector nmopt_full_adjoint =
          nmopt_binding.problem().coordinates().embed_free(nmopt_adjoint);

        external_dealii_step4::problem_b_verification::require_finite(
          native_result.value.control, "native final control");
        external_dealii_step4::problem_b_verification::require_finite(
          native_result.value.state, "native final state");
        external_dealii_step4::problem_b_verification::require_finite(
          native_result.derivative.adjoint, "native final adjoint");
        external_dealii_step4::problem_b_verification::require_finite(
          native_result.derivative.reduced_derivative,
          "native final reduced gradient");
        external_dealii_step4::problem_b_verification::require_finite(
          native_result.value.objective, "native final objective");
        external_dealii_step4::problem_b_verification::require_finite(
          nmopt_control, "nmopt final control");
        external_dealii_step4::problem_b_verification::require_finite(
          nmopt_state, "nmopt final state");
        external_dealii_step4::problem_b_verification::require_finite(
          nmopt_adjoint, "nmopt final adjoint");
        external_dealii_step4::problem_b_verification::require_finite(
          nmopt_gradient, "nmopt final reduced gradient");
        external_dealii_step4::problem_b_verification::require_finite(
          nmopt_result.final_evaluation.objective_value,
          "nmopt final objective");

        external_dealii_step4::problem_b_verification::require_vector_close(
          native_result.value.control,
          nmopt_control,
          1.0e-11,
          1.0e-10,
          "native and nmopt final controls differ");
        external_dealii_step4::problem_b_verification::require_vector_close(
          native_result.value.state,
          nmopt_state,
          1.0e-11,
          1.0e-10,
          "native and nmopt final states differ");
        external_dealii_step4::problem_b_verification::require_vector_close(
          native_result.value.full_state,
          nmopt_full_state,
          1.0e-11,
          1.0e-10,
          "native and nmopt full final states differ");
        external_dealii_step4::problem_b_verification::require_vector_close(
          native_result.derivative.adjoint,
          nmopt_adjoint,
          1.0e-11,
          1.0e-10,
          "native and nmopt final adjoints differ");
        external_dealii_step4::problem_b_verification::require_vector_close(
          native_result.derivative.full_adjoint,
          nmopt_full_adjoint,
          1.0e-11,
          1.0e-10,
          "native and nmopt full final adjoints differ");
        external_dealii_step4::problem_b_verification::require_vector_close(
          native_result.derivative.reduced_derivative,
          nmopt_gradient,
          1.0e-11,
          1.0e-10,
          "native and nmopt final reduced gradients differ");
        require(scalar_matches(native_result.value.objective,
                               nmopt_result.final_evaluation.objective_value),
                "native and nmopt final objectives differ");

        const auto native_runtime_counts =
          external_dealii_step4::runtime_counter_snapshot(
            native_instrumentation);
        const auto nmopt_runtime_counts =
          external_dealii_step4::runtime_counter_snapshot(
            nmopt_instrumentation);

        native_tutorial.output_results(native_result.value.full_state,
                                       artifact_root / "native" / "solution.vtk");
        Step4<2> nmopt_output_tutorial;
        nmopt_output_tutorial.prepare_for_external_use();
        nmopt_output_tutorial.output_results(
          nmopt_full_state, artifact_root / "nmopt" / "solution.vtk");
        require(std::filesystem::exists(artifact_root / "native" /
                                        "solution.vtk") &&
                  std::filesystem::exists(artifact_root / "nmopt" /
                                          "solution.vtk"),
                "Problem B optimization did not retain both state outputs");

        const auto native_vtk =
          read_vtk_payload(artifact_root / "native" / "solution.vtk");
        const auto nmopt_vtk =
          read_vtk_payload(artifact_root / "nmopt" / "solution.vtk");
        require(native_vtk == nmopt_vtk,
                "native and nmopt retained-state VTK payloads differ");
        std::ofstream output_comparison(artifact_root / "comparison" /
                                        "output-comparison.txt");
        require(static_cast<bool>(output_comparison),
                "could not open the Problem B output comparison artifact");
        output_comparison << "timestamp_header_ignored 1\n"
                          << "native_payload_bytes " << native_vtk.size() << '\n'
                          << "nmopt_payload_bytes " << nmopt_vtk.size() << '\n'
                          << "payload_equal 1\n";
        output_comparison.flush();
        require(external_dealii_step4::runtime_counter_snapshot(
                  native_instrumentation) == native_runtime_counts &&
                  external_dealii_step4::runtime_counter_snapshot(
                    nmopt_instrumentation) == nmopt_runtime_counts,
                "retained-state output changed optimization runtime counters");

        Step4<2> verification_tutorial;
        verification_tutorial.prepare_for_external_use();
        NativeProblem verification_problem(verification_tutorial);
        NativeMetric verification_metric(verification_problem.mass());
        NativeReduced verification_reduced(
          verification_problem, native_verification_instrumentation);
        const auto fresh_native_value =
          verification_reduced.evaluate_value(native_result.value.control);
        const auto fresh_native_derivative =
          verification_reduced.augment_derivative(fresh_native_value);
        const auto fresh_nmopt_value =
          verification_reduced.evaluate_value(nmopt_control);
        const auto fresh_nmopt_derivative =
          verification_reduced.augment_derivative(fresh_nmopt_value);

        Binding nmopt_verification_binding(nmopt_verification_instrumentation);
        const auto fresh_public_nmopt_value =
          nmopt_verification_binding.reduced().evaluate_value(
            nmopt::contract::PrimalBlockT<Backend>(
              nmopt_verification_binding.control_layout(), {nmopt_control}));
        const auto fresh_public_nmopt_derivative =
          nmopt_verification_binding.reduced().augment_derivative(
            fresh_public_nmopt_value);
        const auto &fresh_public_nmopt_state =
          fresh_public_nmopt_value.state.block(0);
        const auto &fresh_public_nmopt_adjoint =
          fresh_public_nmopt_derivative.adjoint.block(0);
        const auto &fresh_public_nmopt_gradient =
          fresh_public_nmopt_derivative.reduced_derivative.block(0);

        external_dealii_step4::problem_b_verification::require_vector_close(
          fresh_native_value.state,
          native_result.value.state,
          1.0e-11,
          1.0e-10,
          "fresh native final state differs from the returned state");
        external_dealii_step4::problem_b_verification::require_vector_close(
          fresh_native_value.full_state,
          native_result.value.full_state,
          1.0e-11,
          1.0e-10,
          "fresh native full state differs from the returned state");
        external_dealii_step4::problem_b_verification::require_vector_close(
          fresh_native_derivative.adjoint,
          native_result.derivative.adjoint,
          1.0e-11,
          1.0e-10,
          "fresh native final adjoint differs from the returned adjoint");
        external_dealii_step4::problem_b_verification::require_vector_close(
          fresh_native_derivative.full_adjoint,
          native_result.derivative.full_adjoint,
          1.0e-11,
          1.0e-10,
          "fresh native full adjoint differs from the returned adjoint");
        external_dealii_step4::problem_b_verification::require_vector_close(
          fresh_native_derivative.reduced_derivative,
          native_result.derivative.reduced_derivative,
          1.0e-11,
          1.0e-10,
          "fresh native final gradient differs from the returned gradient");
        require(scalar_matches(fresh_native_value.objective,
                               native_result.value.objective),
                "fresh native final objective differs from the returned value");

        const Vector fresh_public_nmopt_full_state =
          nmopt_binding.problem().coordinates().reconstruct(
            fresh_public_nmopt_state);
        const Vector fresh_public_nmopt_full_adjoint =
          nmopt_binding.problem().coordinates().embed_free(
            fresh_public_nmopt_adjoint);
        external_dealii_step4::problem_b_verification::require_vector_close(
          fresh_public_nmopt_state,
          nmopt_state,
          1.0e-11,
          1.0e-10,
          "fresh nmopt final state differs from the returned state");
        external_dealii_step4::problem_b_verification::require_vector_close(
          fresh_public_nmopt_full_state,
          nmopt_full_state,
          1.0e-11,
          1.0e-10,
          "fresh nmopt full state differs from the returned state");
        external_dealii_step4::problem_b_verification::require_vector_close(
          fresh_public_nmopt_adjoint,
          nmopt_adjoint,
          1.0e-11,
          1.0e-10,
          "fresh nmopt final adjoint differs from the returned adjoint");
        external_dealii_step4::problem_b_verification::require_vector_close(
          fresh_public_nmopt_gradient,
          nmopt_gradient,
          1.0e-11,
          1.0e-10,
          "fresh nmopt final gradient differs from the returned gradient");
        require(scalar_matches(fresh_public_nmopt_value.objective_value,
                               nmopt_result.final_evaluation.objective_value),
                "fresh nmopt final objective differs from the returned value");

        const auto native_state_audit =
          external_dealii_step4::problem_b_verification::audit_state_solution(
            verification_problem,
            verification_tutorial,
            fresh_native_value.state,
            fresh_native_value.full_state,
            fresh_native_value.control);
        const auto native_adjoint_audit =
          external_dealii_step4::problem_b_verification::audit_adjoint_solution(
            verification_problem,
            verification_tutorial,
            fresh_native_derivative.adjoint,
            fresh_native_derivative.full_adjoint,
            fresh_native_derivative.state_derivative);
        const auto nmopt_state_audit =
          external_dealii_step4::problem_b_verification::audit_state_solution(
            verification_problem,
            verification_tutorial,
            fresh_nmopt_value.state,
            fresh_nmopt_value.full_state,
            fresh_nmopt_value.control);
        const auto nmopt_adjoint_audit =
          external_dealii_step4::problem_b_verification::audit_adjoint_solution(
            verification_problem,
            verification_tutorial,
            fresh_nmopt_derivative.adjoint,
            fresh_nmopt_derivative.full_adjoint,
            fresh_nmopt_derivative.state_derivative);
        const auto require_audit = [](const auto &state_audit,
                                      const auto &adjoint_audit,
                                      const char *const path) {
          require(state_audit.full_equation.normalized <= 1.0e-10 &&
                    state_audit.reduced_equation.normalized <= 1.0e-10 &&
                    state_audit.reconstruction_error <= 1.0e-11 &&
                    state_audit.boundary_error <= 1.0e-11,
                  std::string(path) + " final state audit failed");
          require(adjoint_audit.full_equation.normalized <= 1.0e-10 &&
                    adjoint_audit.homogeneous_reconstruction_error <= 1.0e-11 &&
                    adjoint_audit.boundary_error <= 1.0e-11,
                  std::string(path) + " final adjoint audit failed");
        };
        require_audit(native_state_audit, native_adjoint_audit, "native");
        require_audit(nmopt_state_audit, nmopt_adjoint_audit, "nmopt");

        const auto operators =
          external_dealii_step4::problem_b_verification::make_dense_operators(
            verification_problem, verification_tutorial);
        const double native_final_gradient_norm =
          external_dealii_step4::problem_b_verification::dense_mass_norm(
            operators.M, fresh_native_derivative.reduced_derivative);
        const double nmopt_final_gradient_norm =
          external_dealii_step4::problem_b_verification::dense_mass_norm(
            operators.M,
            fresh_public_nmopt_derivative.reduced_derivative.block(0));
        require(native_final_gradient_norm <= 1.1e-6 &&
                  nmopt_final_gradient_norm <= 1.1e-6,
                "final metric gradient exceeds the audit bound");
        require(scalar_matches(native_result.gradient_norm_history.back(),
                               native_final_gradient_norm) &&
                  scalar_matches(nmopt_result.gradient_norm_history.back(),
                                 nmopt_final_gradient_norm),
                "returned and fresh final gradient norms differ");

        const auto oracle =
          external_dealii_step4::problem_b_verification::dense_kkt_oracle(
            operators);
        require(oracle.residuals.state_stationarity <= 1.0e-10 &&
                  oracle.residuals.control_stationarity <= 1.0e-10 &&
                  oracle.residuals.feasibility <= 1.0e-10,
                "Problem B dense KKT oracle failed");
        const auto oracle_value =
          verification_reduced.evaluate_value(oracle.control);
        const auto oracle_derivative =
          verification_reduced.augment_derivative(oracle_value);
        const double oracle_gradient_norm =
          external_dealii_step4::problem_b_verification::dense_mass_norm(
            operators.M, oracle_derivative.reduced_derivative);
        require(oracle_gradient_norm <= 1.0e-8,
                "Problem B oracle reduced gradient is not zero");

        Vector native_control_difference = native_result.value.control;
        native_control_difference.add(-1.0, oracle.control);
        Vector nmopt_control_difference = nmopt_control;
        nmopt_control_difference.add(-1.0, oracle.control);
        const double native_control_distance =
          external_dealii_step4::problem_b_verification::mass_norm(
            verification_metric, native_control_difference);
        const double nmopt_control_distance =
          external_dealii_step4::problem_b_verification::mass_norm(
            verification_metric, nmopt_control_difference);
        require(native_control_distance <= 2.0e-6 &&
                  nmopt_control_distance <= 2.0e-6,
                "final control differs from the dense KKT oracle");

        const double oracle_objective =
          verification_problem.objective(oracle.state, oracle.control);
        const double native_objective_gap =
          native_result.value.objective - oracle_objective;
        const double nmopt_objective_gap =
          nmopt_result.final_evaluation.objective_value - oracle_objective;
        require(std::isfinite(oracle_objective) &&
                  std::isfinite(native_objective_gap) &&
                  std::isfinite(nmopt_objective_gap),
                "Problem B objective gaps are not finite");

        std::ofstream final_audit(artifact_root / "comparison" /
                                  "final-audit.csv");
        require(static_cast<bool>(final_audit),
                "could not open the Problem B final audit artifact");
        final_audit
          << "path,gradient_norm,oracle_control_mass_distance,"
             "objective_gap,state_residual,adjoint_residual\n"
          << std::setprecision(std::numeric_limits<double>::max_digits10)
          << "native," << native_final_gradient_norm << ','
          << native_control_distance << ',' << native_objective_gap << ','
          << native_state_audit.full_equation.normalized << ','
          << native_adjoint_audit.full_equation.normalized << '\n'
          << "nmopt," << nmopt_final_gradient_norm << ','
          << nmopt_control_distance << ',' << nmopt_objective_gap << ','
          << nmopt_state_audit.full_equation.normalized << ','
          << nmopt_adjoint_audit.full_equation.normalized << '\n';
        final_audit.flush();

        std::ofstream oracle_audit(artifact_root / "comparison" /
                                   "oracle-audit.txt");
        require(static_cast<bool>(oracle_audit),
                "could not open the Problem B oracle audit artifact");
        oracle_audit << std::setprecision(std::numeric_limits<double>::max_digits10)
                     << "state_stationarity "
                     << oracle.residuals.state_stationarity << '\n'
                     << "control_stationarity "
                     << oracle.residuals.control_stationarity << '\n'
                     << "feasibility " << oracle.residuals.feasibility << '\n'
                     << "oracle_gradient_norm " << oracle_gradient_norm << '\n'
                     << "oracle_objective " << oracle_objective << '\n';
        oracle_audit.flush();

        require(external_dealii_step4::runtime_counter_snapshot(
                  native_instrumentation) == native_runtime_counts &&
                  external_dealii_step4::runtime_counter_snapshot(
                    nmopt_instrumentation) == nmopt_runtime_counts,
                "final acceptance checks changed optimization runtime counters");
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
      const std::vector<nmopt::test_support::Scenario> scenarios{
        {"problem_b_nmopt_matched_optimization",
         "nmopt.external_tutorial_step_4.problem_b_nmopt_matched_optimization",
         {"dealii", "application", "external", "tutorial", "optimization",
          "problem_b", "verification"},
         360,
         run_problem_b_nmopt_matched_optimization}};
      const auto result = nmopt::test_support::run_requested_scenarios(
        argc, argv, scenarios, std::cout);
      if (!result.listed)
        std::cout << "Step-4 Problem B nmopt optimization scenario passed: "
                  << result.executed << '\n';
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "Step-4 Problem B nmopt optimization test failed: "
                << exception.what() << '\n';
      return 1;
    }
}
