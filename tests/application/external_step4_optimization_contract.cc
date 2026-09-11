#include "../../apps/external-dealii/step-4/evaluation/native_optimization.hpp"
#include "../../apps/external-dealii/step-4/integration/nmopt_binding.hpp"
#include "../../apps/external-dealii/step-4/verification/verification.hpp"

#include "../dealii/external_step4_evidence.hpp"
#include "../support/scenario_dispatch.hpp"

#include <deal.II/lac/vector.h>

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

#include "nmopt/solvers/reduced_gradient.hpp"

namespace
{
  using Backend = nmopt::dealii_backend::SerialBackend;
  using Binding = external_dealii_step4::NmoptBinding;
  using Instrumentation = external_dealii_step4::Instrumentation;
  using Matrix = external_dealii_step4::ProblemA::Matrix;
  using NativeResult = external_dealii_step4::NativeOptimizationResult;
  using NativeReduced = external_dealii_step4::NativeReduced;
  using OptimizationPolicy = external_dealii_step4::OptimizationPolicy;
  using ProblemA = external_dealii_step4::ProblemA;
  using SolveRecord = external_dealii_step4::SolveRecord;
  using TrialRecord = external_dealii_step4::OptimizationTrialRecord;
  using AcceptedRecord =
    external_dealii_step4::OptimizationAcceptedIterationRecord;
  using Vector = ProblemA::Vector;
  using NmoptSolver = nmopt::solvers::ReducedSearchSolverT<Backend>;
  using NmoptResult = nmopt::solvers::ReducedSolverResultT<Backend>;
  namespace verification = external_dealii_step4::verification;

  struct ScalarTrace
  {
    std::vector<double> gradient_norms;
    std::vector<TrialRecord> trials;
    std::vector<AcceptedRecord> accepted;
    std::string stopping_reason;
    double final_objective = 0.0;
    double final_gradient_norm = 0.0;
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

  double
  vector_difference(const Vector &left, const Vector &right)
  {
    require(left.size() == right.size(),
            "paired optimization vectors have incompatible dimensions");
    Vector difference = left;
    difference.add(-1.0, right);
    return difference.l2_norm();
  }

  bool
  scalar_matches(const double left, const double right)
  {
    return std::isfinite(left) && std::isfinite(right) &&
           std::abs(left - right) <=
             1.0e-12 + 1.0e-11 * std::max(std::abs(left), std::abs(right));
  }

  bool
  vector_matches(const Vector &left, const Vector &right)
  {
    if (!verification::vector_is_finite(left) ||
        !verification::vector_is_finite(right))
      return false;
    const double error = vector_difference(left, right);
    const double bound =
      1.0e-11 + 1.0e-10 * std::max(left.l2_norm(), right.l2_norm());
    return std::isfinite(error) && std::isfinite(bound) && error <= bound;
  }

  double
  matrix_difference(const Matrix &left, const Matrix &right)
  {
    require(left.m() == right.m() && left.n() == right.n(),
            "paired optimization matrices have incompatible shapes");
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
    return {result.gradient_norm_history,
            result.trial_records,
            result.accepted_iterations,
            external_dealii_step4::native_optimization_stopping_reason_name(
              result.stopping_reason),
            result.value.objective,
            result.derivative.reduced_derivative.l2_norm()};
  }

  ScalarTrace
  nmopt_trace(const NmoptResult &result)
  {
    ScalarTrace trace;
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

  std::filesystem::path
  create_artifact_root()
  {
    auto directory = std::filesystem::current_path();
    while (true)
      {
        if (std::filesystem::exists(
              directory / "apps/external-dealii/step-4/source/upstream/step-4.cc"))
          {
            const auto root =
              external_dealii_step4_test::create_unique_artifact_root(
                directory / "runs/external-dealii/step-4/optimization",
                "paired");
            std::filesystem::create_directories(root / "native");
            std::filesystem::create_directories(root / "nmopt");
            std::filesystem::create_directories(root / "comparison");
            return root;
          }
        const auto parent = directory.parent_path();
        if (parent == directory)
          break;
        directory = parent;
      }
    throw std::runtime_error("could not locate the optimization artifact root");
  }

  void
  write_trace(const std::filesystem::path &root,
              const char *const            path_name,
              const ScalarTrace &          trace)
  {
    std::ofstream output(root / "trace.csv");
    require(static_cast<bool>(output), "could not open an optimization trace");
    output << "path " << path_name << '\n'
           << "stopping_reason " << trace.stopping_reason << '\n'
           << "accepted_iterations " << trace.accepted.size() << '\n'
           << "line_search_trials " << trace.trials.size() << '\n'
           << "final_objective " << std::setprecision(
             std::numeric_limits<double>::max_digits10)
           << trace.final_objective << '\n'
           << "final_gradient_norm " << trace.final_gradient_norm << '\n'
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
             << ",,,,1," << iteration.objective_before << ','
             << iteration.objective_after << ',' << iteration.objective_change
             << ',' << iteration.actual_step_norm << ','
             << iteration.gradient_norm << '\n';
  }

  void
  compare_traces(const ScalarTrace &left,
                 const ScalarTrace &right,
                 std::string &       first_divergence)
  {
    const auto note = [&](const std::string &message) {
      if (first_divergence.empty())
        first_divergence = message;
    };
    if (left.stopping_reason != right.stopping_reason)
      note("stopping reason differs");
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
    else
      for (std::size_t index = 0; index < left.trials.size(); ++index)
        {
          const auto &a = left.trials[index];
          const auto &b = right.trials[index];
          if (a.iteration != b.iteration || a.trial != b.trial)
            note("trial index first differs at record " +
                 std::to_string(index));
          if (!scalar_matches(a.step_length, b.step_length) ||
              !scalar_matches(a.objective_value, b.objective_value) ||
              !scalar_matches(a.actual_slope, b.actual_slope) ||
              !scalar_matches(a.sufficient_decrease_bound,
                              b.sufficient_decrease_bound))
            note("trial scalar first differs at record " +
                 std::to_string(index));
          if (a.objective_finite != b.objective_finite ||
              a.slope_negative != b.slope_negative ||
              a.accepted != b.accepted)
            note("trial acceptance first differs at record " +
                 std::to_string(index));
        }

    if (left.accepted.size() != right.accepted.size())
      note("accepted-record size differs");
    else
      for (std::size_t index = 0; index < left.accepted.size(); ++index)
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
  }

  void
  compare_solve_records(const std::vector<SolveRecord> &left,
                        const std::vector<SolveRecord> &right,
                        std::string &first_divergence)
  {
    const auto note = [&](const std::string &message) {
      if (first_divergence.empty())
        first_divergence = message;
    };
    if (left.size() != right.size())
      note("solve-record count differs");
    const auto count = std::min(left.size(), right.size());
    for (std::size_t index = 0; index < count; ++index)
      {
        if (left[index].role != right[index].role ||
            left[index].iterations != right[index].iterations)
          note("CG work first differs at record " + std::to_string(index));
        if (!scalar_matches(left[index].initial_residual,
                            right[index].initial_residual) ||
            !scalar_matches(left[index].final_residual,
                            right[index].final_residual))
          note("CG residual first differs at record " +
               std::to_string(index));
      }
  }

  std::string
  read_output_payload(const std::filesystem::path &filename)
  {
    std::ifstream input(filename, std::ios::binary);
    require(static_cast<bool>(input), "could not open output comparison file");
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

  std::string
  read_file(const std::filesystem::path &filename)
  {
    std::ifstream input(filename, std::ios::binary);
    require(static_cast<bool>(input), "could not open Step-4 evidence file");
    return {std::istreambuf_iterator<char>(input),
            std::istreambuf_iterator<char>()};
  }

  void
  run_matched_optimization()
  {
    Instrumentation native_instrumentation;
    Instrumentation nmopt_instrumentation;
    Instrumentation verification_instrumentation;
    const auto root = create_artifact_root();
    external_dealii_step4_test::EvidenceGuard evidence(
      root,
      "matched_optimization",
      {{"native", &native_instrumentation},
       {"nmopt", &nmopt_instrumentation},
       {"verification", &verification_instrumentation}});
    try
      {
    ProblemA      native_problem(native_instrumentation);
    NativeReduced native_reduced(native_problem, native_instrumentation);
    Binding       nmopt_binding(nmopt_instrumentation);
    require(matrix_difference(native_problem.system_matrix(),
                              nmopt_binding.problem().system_matrix()) == 0.0,
            "paired optimization matrices differ");
    require(vector_difference(native_problem.system_rhs(),
                              nmopt_binding.problem().system_rhs()) == 0.0,
            "paired optimization RHS vectors differ");

    const auto policy = external_dealii_step4::frozen_optimization_policy();
    Vector initial_control(native_problem.control_dimension());
    initial_control = 0.0;
    const auto native_result =
      external_dealii_step4::NativeArmijoSolver(native_reduced, policy)
        .solve(initial_control);
    const auto native_trace_value = native_trace(native_result);
    write_trace(root / "native", "native", native_trace_value);
    const auto parameters = nmopt_parameters(policy);
    NmoptSolver nmopt_solver(nmopt_binding.reduced(),
                             nmopt_binding.metric(),
                             parameters);
    const auto nmopt_result = nmopt_solver.solve(
      nmopt::contract::PrimalBlockT<Backend>(nmopt_binding.control_layout(),
                                             {initial_control}));
    const auto nmopt_trace_value = nmopt_trace(nmopt_result);
    write_trace(root / "nmopt", "nmopt", nmopt_trace_value);

    verification::require_finite(native_result.value.control,
                                  "native optimization final control");
    verification::require_finite(native_result.value.state,
                                  "native optimization final state");
    verification::require_finite(native_result.value.objective,
                                  "native optimization final objective");
    verification::require_finite(
      native_result.derivative.reduced_derivative,
      "native optimization final gradient");

    require(native_result.stopping_reason ==
              external_dealii_step4::NativeOptimizationStoppingReason::
                gradient_tolerance,
            "native optimizer did not stop by gradient tolerance");
    require(nmopt_result.stopping_reason ==
              nmopt::solvers::ReducedStoppingReason::gradient_tolerance,
            "nmopt optimizer did not stop by gradient tolerance");

    const auto &nmopt_control =
      nmopt_result.final_evaluation.full_point.block(1);
    const auto &nmopt_state = nmopt_result.final_evaluation.state.block(0);
    const auto &nmopt_gradient =
      nmopt_result.final_evaluation.reduced_derivative.block(0);
    verification::require_finite(nmopt_control,
                                 "nmopt optimization final control");
    verification::require_finite(nmopt_state,
                                 "nmopt optimization final state");
    verification::require_finite(nmopt_result.final_evaluation.objective_value,
                                 "nmopt optimization final objective");
    verification::require_finite(nmopt_gradient,
                                 "nmopt optimization final gradient");
    const auto oracle = verification::optimum_oracle(native_problem);
    native_problem.output_results(native_result.value.state,
                                  root / "native" / "solution.vtk");
    nmopt_binding.problem().output_results(nmopt_state,
                                           root / "nmopt" / "solution.vtk");

    std::string first_divergence;
    const auto note = [&](const std::string &message) {
      if (first_divergence.empty())
        first_divergence = message;
    };
    compare_traces(native_trace_value, nmopt_trace_value, first_divergence);
    compare_solve_records(native_instrumentation.solve_records,
                          nmopt_instrumentation.solve_records,
                          first_divergence);
    if (!vector_matches(native_result.value.control, nmopt_control))
      note("final control differs");
    if (!vector_matches(native_result.value.state, nmopt_state))
      note("final state differs");
    if (!vector_matches(native_result.derivative.reduced_derivative,
                        nmopt_gradient))
      note("final reduced gradient differs");
    if (!scalar_matches(native_result.value.objective,
                        nmopt_result.final_evaluation.objective_value))
      note("final objective differs");
    if (oracle.system_residual > 1.0e-10 ||
        oracle.stationarity_residual > 1.0e-10)
      note("independent oracle audit failed");
    const double native_control_oracle_error =
      vector_difference(native_result.value.control, oracle.control);
    const double nmopt_control_oracle_error =
      vector_difference(nmopt_control, oracle.control);
    if (native_control_oracle_error > 2.0e-6 ||
        nmopt_control_oracle_error > 2.0e-6)
      note("final control failed the oracle audit");
    if (read_output_payload(root / "native" / "solution.vtk") !=
        read_output_payload(root / "nmopt" / "solution.vtk"))
      note("retained-state output differs");

    const auto accepted = native_result.accepted_iteration_count;
    const auto trials = native_result.line_search_trial_count;
    if (native_instrumentation.assembly_calls != 1 ||
        nmopt_instrumentation.assembly_calls != 1)
      note("assembly schedule is inconsistent");
    if (native_instrumentation.state_solve_calls != 1 + trials ||
        nmopt_instrumentation.state_solve_calls != 1 + trials)
      note("state schedule is inconsistent");
    if (native_instrumentation.adjoint_solve_calls != 1 + accepted ||
        nmopt_instrumentation.adjoint_solve_calls != 1 + accepted)
      note("adjoint schedule is inconsistent");
    if (native_instrumentation.objective_calls != 1 + trials ||
        nmopt_instrumentation.objective_calls != 1 + trials ||
        native_instrumentation.objective_derivative_calls != 1 + accepted ||
        nmopt_instrumentation.objective_derivative_calls != 1 + accepted)
      note("objective schedule is inconsistent");
    if (native_instrumentation.value_evaluations != 1 + trials ||
        nmopt_instrumentation.value_evaluations != 0 ||
        native_instrumentation.derivative_augmentations != 1 + accepted ||
        nmopt_instrumentation.derivative_augmentations != 0)
      note("staged evaluation counters are inconsistent");
    if (native_instrumentation.solve_failures != 0 ||
        nmopt_instrumentation.solve_failures != 0)
      note("optimization recorded a solve failure");
    if (native_instrumentation.residual_calls != 0 ||
        native_instrumentation.residual_jvp_calls != 0 ||
        nmopt_instrumentation.residual_calls != 0 ||
        nmopt_instrumentation.residual_jvp_calls != 0)
      note("optimization unexpectedly used residual value or JVP callbacks");
    if (native_instrumentation.control_vjp_calls != 1 + accepted ||
        native_instrumentation.residual_vjp_calls != 0 ||
        nmopt_instrumentation.control_vjp_calls != 0 ||
        nmopt_instrumentation.residual_vjp_calls != 1 + accepted ||
        native_instrumentation.explicit_matrix_vmult_calls != 0 ||
        nmopt_instrumentation.explicit_matrix_vmult_calls != 0 ||
        nmopt_instrumentation.explicit_matrix_tvmult_calls != 1 + accepted)
      note("pullback schedule is inconsistent");
    if (nmopt_result.metric_solve_count != accepted + 1 ||
        nmopt_result.hessian_action_count != 0 ||
        native_instrumentation.metric_apply_calls != 0 ||
        native_instrumentation.metric_inverse_apply_calls != 0 ||
        nmopt_instrumentation.metric_inverse_apply_calls !=
          nmopt_result.metric_solve_count ||
        nmopt_instrumentation.metric_apply_calls !=
          nmopt_result.metric_solve_count + accepted)
      note("nmopt metric callback counts are inconsistent");
    if (native_instrumentation.output_calls != 1 ||
        nmopt_instrumentation.output_calls != 1)
      note("output schedule is inconsistent");

    const auto native_runtime_counts =
      external_dealii_step4::runtime_counter_snapshot(
        native_instrumentation);
    const auto nmopt_runtime_counts =
      external_dealii_step4::runtime_counter_snapshot(nmopt_instrumentation);

    const auto &nmopt_adjoint = nmopt_result.final_evaluation.adjoint.block(0);
    const auto native_state_audit = verification::state_equation_residual(
      native_problem, native_result.value.state, native_result.value.control);
    const auto native_adjoint_audit =
      verification::adjoint_equation_residual(native_problem,
                                              native_result.derivative.adjoint,
                                              native_result.value.state);
    const auto nmopt_state_audit = verification::state_equation_residual(
      nmopt_binding.problem(), nmopt_state, nmopt_control);
    const auto nmopt_adjoint_audit = verification::adjoint_equation_residual(
      nmopt_binding.problem(), nmopt_adjoint, nmopt_state);
    require(native_state_audit.normalized <= 1.0e-10,
            "native final state residual audit failed");
    require(native_adjoint_audit.normalized <= 1.0e-10,
            "native final adjoint residual audit failed");
    require(nmopt_state_audit.normalized <= 1.0e-10,
            "nmopt final state residual audit failed");
    require(nmopt_adjoint_audit.normalized <= 1.0e-10,
            "nmopt final adjoint residual audit failed");
    require(external_dealii_step4::runtime_counter_snapshot(
              native_instrumentation) == native_runtime_counts,
            "native final residual audits changed runtime counters");
    require(external_dealii_step4::runtime_counter_snapshot(
              nmopt_instrumentation) == nmopt_runtime_counts,
            "nmopt final residual audits changed runtime counters");

    ProblemA      verification_problem(verification_instrumentation);
    NativeReduced verification_reduced(verification_problem,
                                       verification_instrumentation);
    const auto fresh_native_value =
      verification_reduced.evaluate_value(native_result.value.control);
    const auto fresh_native_derivative =
      verification_reduced.augment_derivative(fresh_native_value);
    const auto fresh_nmopt_value =
      verification_reduced.evaluate_value(nmopt_control);
    const auto fresh_nmopt_derivative =
      verification_reduced.augment_derivative(fresh_nmopt_value);
    const double fresh_native_gradient_norm =
      fresh_native_derivative.reduced_derivative.l2_norm();
    const double fresh_nmopt_gradient_norm =
      fresh_nmopt_derivative.reduced_derivative.l2_norm();
    const double native_gradient_difference = vector_difference(
      fresh_native_derivative.reduced_derivative,
      native_result.derivative.reduced_derivative);
    const double nmopt_gradient_difference = vector_difference(
      fresh_nmopt_derivative.reduced_derivative, nmopt_gradient);
    require(fresh_native_gradient_norm <= 1.1e-6,
            "fresh native final gradient exceeds the audit bound");
    require(fresh_nmopt_gradient_norm <= 1.1e-6,
            "fresh nmopt final gradient exceeds the audit bound");
    verification::require_vector_close(
      fresh_native_derivative.reduced_derivative,
      native_result.derivative.reduced_derivative,
      1.0e-11,
      1.0e-10,
      "fresh native final gradient differs from the returned gradient");
    verification::require_vector_close(
      fresh_nmopt_derivative.reduced_derivative,
      nmopt_gradient,
      1.0e-11,
      1.0e-10,
      "fresh nmopt final gradient differs from the returned gradient");
    require(verification_instrumentation.assembly_calls == 1 &&
              verification_instrumentation.state_solve_calls == 2 &&
              verification_instrumentation.adjoint_solve_calls == 2 &&
              verification_instrumentation.objective_calls == 2 &&
              verification_instrumentation.objective_derivative_calls == 2 &&
              verification_instrumentation.control_vjp_calls == 2 &&
              verification_instrumentation.residual_vjp_calls == 0,
            "fresh gradient verification did not recompute both native evaluations");
    require(external_dealii_step4::runtime_counter_snapshot(
              native_instrumentation) == native_runtime_counts,
            "fresh gradient verification changed native runtime counters");
    require(external_dealii_step4::runtime_counter_snapshot(
              nmopt_instrumentation) == nmopt_runtime_counts,
            "fresh gradient verification changed nmopt runtime counters");

    std::ofstream gradient_audit(root / "comparison" / "gradient-audit.csv");
    require(static_cast<bool>(gradient_audit),
            "could not open the optimization gradient audit artifact");
    gradient_audit
      << "path,returned_gradient_norm,fresh_gradient_norm,"
         "gradient_difference,fresh_state_monitored_residual,"
         "fresh_adjoint_monitored_residual\n"
      << std::setprecision(std::numeric_limits<double>::max_digits10)
      << "native," << native_result.derivative.reduced_derivative.l2_norm()
      << ',' << fresh_native_gradient_norm << ','
      << native_gradient_difference << ','
      << fresh_native_value.state_solve.final_residual << ','
      << fresh_native_derivative.adjoint_solve.final_residual << '\n'
      << "nmopt," << nmopt_gradient.l2_norm() << ','
      << fresh_nmopt_gradient_norm << ',' << nmopt_gradient_difference << ','
      << fresh_nmopt_value.state_solve.final_residual << ','
      << fresh_nmopt_derivative.adjoint_solve.final_residual << '\n';
    gradient_audit.flush();

    std::ofstream residual_audits(root / "comparison" /
                                  "residual-audits.csv");
    require(static_cast<bool>(residual_audits),
            "could not open optimization residual audit artifact");
    residual_audits
      << "path,state_monitored_residual,state_recomputed_absolute_norm,"
         "state_recomputed_normalized_residual,state_residual_scale,"
         "adjoint_monitored_residual,adjoint_recomputed_absolute_norm,"
         "adjoint_recomputed_normalized_residual,adjoint_residual_scale\n"
      << std::setprecision(std::numeric_limits<double>::max_digits10)
      << "native," << native_result.value.state_solve.final_residual << ','
      << native_state_audit.absolute_norm << ','
      << native_state_audit.normalized << ',' << native_state_audit.scale
      << ',' << native_result.derivative.adjoint_solve.final_residual << ','
      << native_adjoint_audit.absolute_norm << ','
      << native_adjoint_audit.normalized << ',' << native_adjoint_audit.scale
      << '\n'
      << "nmopt," << nmopt_result.final_evaluation.state_solve.achieved_residual
      << ',' << nmopt_state_audit.absolute_norm << ','
      << nmopt_state_audit.normalized << ',' << nmopt_state_audit.scale
      << ',' << nmopt_result.final_evaluation.adjoint_solve.achieved_residual
      << ',' << nmopt_adjoint_audit.absolute_norm << ','
      << nmopt_adjoint_audit.normalized << ','
      << nmopt_adjoint_audit.scale << '\n';
    residual_audits.flush();

    std::ofstream summary(root / "comparison" / "summary.txt");
    require(static_cast<bool>(summary),
            "could not open the optimization comparison summary");
    summary << std::setprecision(std::numeric_limits<double>::max_digits10)
            << "first_divergence "
            << (first_divergence.empty() ? "none" : first_divergence) << '\n'
            << "accepted_iterations " << accepted << '\n'
            << "line_search_trials " << trials << '\n'
            << "oracle_system_residual " << oracle.system_residual << '\n'
            << "oracle_stationarity_residual "
            << oracle.stationarity_residual << '\n'
            << "oracle_control_distance_bound " << 2.0e-6 << '\n'
            << "native_control_oracle_error " << native_control_oracle_error
            << '\n'
            << "nmopt_control_oracle_error " << nmopt_control_oracle_error
            << '\n'
            << "oracle_objective "
            << (0.5 * (oracle.state * oracle.state) +
                0.5 * (oracle.control * oracle.control)) << '\n'
            << "native_objective_gap "
            << (native_result.value.objective -
                (0.5 * (oracle.state * oracle.state) +
                 0.5 * (oracle.control * oracle.control))) << '\n'
            << "nmopt_objective_gap "
            << (nmopt_result.final_evaluation.objective_value -
                (0.5 * (oracle.state * oracle.state) +
                 0.5 * (oracle.control * oracle.control))) << '\n'
            << "native_final_state_monitored_residual "
            << native_result.value.state_solve.final_residual << '\n'
            << "native_final_state_recomputed_normalized_residual "
            << native_state_audit.normalized << '\n'
            << "native_final_adjoint_monitored_residual "
            << native_result.derivative.adjoint_solve.final_residual << '\n'
            << "native_final_adjoint_recomputed_normalized_residual "
            << native_adjoint_audit.normalized << '\n'
            << "nmopt_final_state_monitored_residual "
            << nmopt_result.final_evaluation.state_solve.achieved_residual
            << '\n'
            << "nmopt_final_state_recomputed_normalized_residual "
            << nmopt_state_audit.normalized << '\n'
            << "nmopt_final_adjoint_monitored_residual "
            << nmopt_result.final_evaluation.adjoint_solve.achieved_residual
            << '\n'
            << "nmopt_final_adjoint_recomputed_normalized_residual "
            << nmopt_adjoint_audit.normalized << '\n'
            << "native_final_gradient_norm "
            << native_trace_value.final_gradient_norm << '\n'
            << "nmopt_final_gradient_norm "
            << nmopt_trace_value.final_gradient_norm << '\n'
            << "final_control_error "
            << vector_difference(native_result.value.control, nmopt_control)
            << '\n'
            << "final_state_error "
            << vector_difference(native_result.value.state, nmopt_state) << '\n'
            << "native_control_vjp "
            << native_instrumentation.control_vjp_calls << '\n'
            << "nmopt_residual_vjp "
            << nmopt_instrumentation.residual_vjp_calls << '\n'
            << "nmopt_explicit_matrix_tvmult "
            << nmopt_instrumentation.explicit_matrix_tvmult_calls << '\n'
            << "native_metric_apply "
            << native_instrumentation.metric_apply_calls << '\n'
            << "nmopt_metric_apply "
            << nmopt_instrumentation.metric_apply_calls << '\n'
            << "nmopt_metric_inverse_apply "
            << nmopt_instrumentation.metric_inverse_apply_calls << '\n'
            << "nmopt_reported_metric_solves "
            << nmopt_result.metric_solve_count << '\n';
    summary.flush();
    require(first_divergence.empty(), first_divergence);
    evidence.complete();
      }
    catch (...)
      {
        evidence.fail_current_exception();
        throw;
      }
  }

  void
  run_completed_trace_failure_probe(std::filesystem::path &artifact)
  {
    Instrumentation native_instrumentation;
    Instrumentation nmopt_instrumentation;
    artifact = create_artifact_root();
    external_dealii_step4_test::EvidenceGuard evidence(
      artifact,
      "completed_trace_failure",
      {{"native", &native_instrumentation},
       {"nmopt", &nmopt_instrumentation}});
    try
      {
        ProblemA      native_problem(native_instrumentation);
        NativeReduced native_reduced(native_problem, native_instrumentation);
        Binding       nmopt_binding(nmopt_instrumentation);
        Vector        initial_control(native_problem.control_dimension());
        initial_control = 0.0;
        const auto policy = external_dealii_step4::frozen_optimization_policy();
        const auto native_result =
          external_dealii_step4::NativeArmijoSolver(native_reduced, policy)
            .solve(initial_control);
        write_trace(artifact / "native", "native", native_trace(native_result));

        NmoptSolver nmopt_solver(nmopt_binding.reduced(),
                                 nmopt_binding.metric(),
                                 nmopt_parameters(policy));
        const auto nmopt_result = nmopt_solver.solve(
          nmopt::contract::PrimalBlockT<Backend>(nmopt_binding.control_layout(),
                                                 {initial_control}));
        write_trace(artifact / "nmopt", "nmopt", nmopt_trace(nmopt_result));

        throw std::runtime_error("failure after completed optimization traces");
      }
    catch (...)
      {
        evidence.fail_current_exception();
        throw;
      }
  }

  void
  require_completed_trace_failure_artifact(
    const std::filesystem::path &artifact)
  {
    const auto status = read_file(artifact / "status.txt");
    const auto failure = read_file(artifact / "failure.txt");
    const auto counters = read_file(artifact / "counters.csv");
    const auto solves = read_file(artifact / "solve-records.csv");
    const auto native_trace_contents = read_file(artifact / "native" / "trace.csv");
    const auto nmopt_trace_contents = read_file(artifact / "nmopt" / "trace.csv");
    require(status.find("status failed") != std::string::npos,
            "completed trace failure did not retain failed status");
    require(failure.find("failure after completed optimization traces") !=
              std::string::npos,
            "completed trace failure lost its original diagnostic");
    require(counters.find("native,state_solve_calls,") != std::string::npos &&
              counters.find("nmopt,state_solve_calls,") != std::string::npos,
            "completed trace failure lost optimization counters");
    require(solves.find("native,success,state") != std::string::npos &&
              solves.find("nmopt,success,state") != std::string::npos,
            "completed trace failure lost solve records");
    require(native_trace_contents.find("accepted_iterations ") !=
              std::string::npos &&
              native_trace_contents.find("record,iteration") !=
                std::string::npos,
            "completed native optimization trace was not retained");
    require(nmopt_trace_contents.find("accepted_iterations ") !=
              std::string::npos &&
              nmopt_trace_contents.find("record,iteration") !=
                std::string::npos,
            "completed nmopt optimization trace was not retained");
  }

  void
  run_completed_trace_failure_contract()
  {
    std::filesystem::path first;
    std::filesystem::path retry;
    bool                  first_thrown = false;
    try
      {
        run_completed_trace_failure_probe(first);
      }
    catch (const std::exception &exception)
      {
        first_thrown = true;
        require(exception.what() ==
                  std::string("failure after completed optimization traces"),
                "completed trace failure did not propagate its diagnostic");
      }
    require(first_thrown, "completed trace failure probe did not throw");

    bool retry_thrown = false;
    try
      {
        run_completed_trace_failure_probe(retry);
      }
    catch (const std::exception &exception)
      {
        retry_thrown = true;
        require(exception.what() ==
                  std::string("failure after completed optimization traces"),
                "completed trace failure retry did not propagate its diagnostic");
      }
    require(retry_thrown, "completed trace failure retry did not throw");
    require(first != retry,
            "completed trace failure retry reused the run directory");
    require_completed_trace_failure_artifact(first);
    require_completed_trace_failure_artifact(retry);
  }

  void
  run_matched_optimization_limit()
  {
    Instrumentation native_instrumentation;
    Instrumentation nmopt_instrumentation;
    ProblemA        native_problem(native_instrumentation);
    NativeReduced   native_reduced(native_problem, native_instrumentation);
    Binding         nmopt_binding(nmopt_instrumentation);
    Vector          initial_control(native_problem.control_dimension());
    initial_control = 0.0;

    auto policy = external_dealii_step4::frozen_optimization_policy();
    policy.maximum_iterations = 1;
    const auto native_result =
      external_dealii_step4::NativeArmijoSolver(native_reduced, policy)
        .solve(initial_control);
    NmoptSolver nmopt_solver(nmopt_binding.reduced(),
                             nmopt_binding.metric(),
                             nmopt_parameters(policy));
    const auto nmopt_result = nmopt_solver.solve(
      nmopt::contract::PrimalBlockT<Backend>(nmopt_binding.control_layout(),
                                             {initial_control}));

    require(native_result.stopping_reason ==
              external_dealii_step4::NativeOptimizationStoppingReason::
                maximum_iterations,
            "native iteration-limit probe returned a successful stopping reason");
    require(nmopt_result.stopping_reason ==
              nmopt::solvers::ReducedStoppingReason::maximum_iterations,
            "nmopt iteration-limit probe returned a successful stopping reason");
    require(native_result.accepted_iteration_count == 1 &&
              nmopt_result.accepted_iterations == 1,
            "paired iteration-limit probe did not accept exactly one step");
    require(native_result.gradient_norm_history.size() == 2 &&
              nmopt_result.gradient_norm_history.size() == 2,
            "paired iteration-limit probe did not recheck the gradient");

    const double native_gradient_norm =
      native_result.derivative.reduced_derivative.l2_norm();
    const double nmopt_gradient_norm =
      nmopt_result.gradient_norm_history.back();
    require(std::isfinite(native_gradient_norm) &&
              native_gradient_norm > policy.gradient_tolerance &&
              std::isfinite(nmopt_gradient_norm) &&
              nmopt_gradient_norm > policy.gradient_tolerance,
            "paired iteration-limit probe was falsely accepted as converged");
  }
} // namespace

int
main(const int argc, char **argv)
{
  try
    {
      const std::vector<nmopt::test_support::Scenario> scenarios{
        {"matched_optimization",
         "nmopt.external.tutorial_step_4.matched_optimization",
         {"dealii", "application", "external", "tutorial", "optimization"},
         360,
         run_matched_optimization},
        {"completed_trace_failure",
         "nmopt.external.tutorial_step_4.completed_trace_failure",
         {"dealii", "application", "external", "tutorial", "optimization", "diagnostics"},
         180,
         run_completed_trace_failure_contract},
        {"matched_optimization_limit",
         "nmopt.external.tutorial_step_4.matched_optimization_limit",
         {"dealii", "application", "external", "tutorial", "optimization"},
         60,
         run_matched_optimization_limit}};
      const auto result = nmopt::test_support::run_requested_scenarios(
        argc, argv, scenarios, std::cout);
      if (!result.listed)
        std::cout << "external Step4 matched optimization scenario passed: "
                  << result.executed << '\n';
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "external Step4 matched optimization test failed: "
                << exception.what() << '\n';
      return 1;
    }
}
