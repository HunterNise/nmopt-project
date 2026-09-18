#include "../../apps/external-dealii/step-4/integration/nmopt_binding.hpp"
#include "../../apps/external-dealii/step-4/evaluation/native_reduced.hpp"
#include "../../apps/external-dealii/step-4/verification/verification.hpp"

#include "../dealii/external_step4_evidence.hpp"
#include "../support/scenario_dispatch.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace
{
  using Backend      = nmopt::dealii_backend::SerialBackend;
  using Binding      = external_dealii_step4::NmoptBinding;
  using Instrumentation = external_dealii_step4::Instrumentation;
  using Model        = Binding::Model;
  using Primal       = Binding::Primal;
  using Covector     = Binding::Covector;
  using LayoutPtr    = external_dealii_step4::LayoutPtr;
  using ProblemA     = external_dealii_step4::ProblemA;
  using Vector       = ProblemA::Vector;

  static_assert(!std::is_copy_constructible_v<Binding>);
  static_assert(!std::is_copy_assignable_v<Binding>);
  static_assert(!std::is_move_constructible_v<Binding>);
  static_assert(!std::is_move_assignable_v<Binding>);
  using Matrix       = ProblemA::Matrix;
  using NativeReduced = external_dealii_step4::NativeReduced;
  namespace verification = external_dealii_step4::verification;

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
            "Step-4 nmopt binding vectors have incompatible dimensions");
    Vector difference = left;
    difference.add(-1.0, right);
    return difference.l2_norm();
  }

  void
  require_vector_equal(const Vector &      actual,
                       const Vector &      expected,
                       const std::string &message)
  {
    require(vector_difference(actual, expected) == 0.0, message);
  }

  void
  require_report(const nmopt::contract::LinearSolveReport &report,
                 const std::string                         &name)
  {
    require(report.converged(), name + " did not report convergence");
    require(report.algorithm == "CG", name + " used the wrong algorithm");
    require(report.preconditioner == "identity",
            name + " used the wrong preconditioner");
    require(report.maximum_iterations == 1000,
            name + " used the wrong iteration limit");
    require(report.iterations > 0, name + " reported no CG work");
    require(report.relative_tolerance == 0.0,
            name + " reported a relative tolerance");
    require(report.absolute_tolerance == 1.0e-12,
            name + " reported the wrong absolute tolerance");
    require(report.requested_tolerance == 1.0e-12,
            name + " reported the wrong requested tolerance");
    require(std::isfinite(report.achieved_residual) &&
              report.achieved_residual >= 0.0,
            name + " reported an invalid achieved residual");
  }

  struct OffSolutionPoint
  {
    Vector state;
    Vector control;
    Vector state_tangent;
    Vector control_tangent;
    Vector test_seed;
  };

  OffSolutionPoint
  make_off_solution_point(const std::size_t dimension)
  {
    const auto constant    = external_dealii_step4::scenario::constant_vector();
    const auto ramp        = external_dealii_step4::scenario::ramp_vector();
    const auto alternating =
      external_dealii_step4::scenario::alternating_vector();

    Vector state = external_dealii_step4::scenario::scaled(constant, 0.2);
    state.add(0.3, ramp);
    Vector control =
      external_dealii_step4::scenario::scaled(alternating, 0.2);
    control.add(0.1, constant);
    Vector seed = external_dealii_step4::scenario::scaled(constant, 0.4);
    seed.add(0.2, alternating);

    require(state.size() == dimension && control.size() == dimension &&
              seed.size() == dimension,
            "Step-4 nmopt binding scenario dimension changed");
    return {std::move(state),
            std::move(control),
            external_dealii_step4::scenario::normalized_ramp_direction(),
            external_dealii_step4::scenario::normalized_alternating_direction(),
            std::move(seed)};
  }

  Primal
  make_full_point(const LayoutPtr &layout,
                  const Vector &   state,
                  const Vector &   control)
  {
    return Primal(layout, {state, control});
  }

  Primal
  make_single_block(const LayoutPtr &layout, const Vector &value)
  {
    return Primal(layout, {value});
  }

  template <typename Callable>
  void
  require_rejected(Callable &&callable, const std::string &message)
  {
    bool rejected = false;
    try
      {
        callable();
      }
    catch (const std::exception &)
      {
        rejected = true;
      }
    require(rejected, message);
  }

  LayoutPtr
  compatible_variable_layout(const Binding &binding)
  {
    return std::make_shared<const nmopt::contract::BlockLayout>(
      "alternate_step4_variables",
      std::vector<nmopt::contract::SpaceId>{{"state"}, {"control"}},
      std::vector<std::size_t>{binding.problem().state_dimension(),
                               binding.problem().control_dimension()});
  }

  LayoutPtr
  incompatible_variable_layout(const Binding &binding)
  {
    return std::make_shared<const nmopt::contract::BlockLayout>(
      "incompatible_step4_variables",
      std::vector<nmopt::contract::SpaceId>{{"not_state"}, {"control"}},
      std::vector<std::size_t>{binding.problem().state_dimension(),
                               binding.problem().control_dimension()});
  }

  LayoutPtr
  compatible_test_layout(const Binding &binding)
  {
    return std::make_shared<const nmopt::contract::BlockLayout>(
      "alternate_step4_test",
      std::vector<nmopt::contract::SpaceId>{{"state_test"}},
      std::vector<std::size_t>{binding.problem().state_dimension()});
  }

  LayoutPtr
  incompatible_test_layout(const Binding &binding)
  {
    return std::make_shared<const nmopt::contract::BlockLayout>(
      "incompatible_step4_test",
      std::vector<nmopt::contract::SpaceId>{{"not_state_test"}},
      std::vector<std::size_t>{binding.problem().state_dimension()});
  }

  LayoutPtr
  compatible_control_layout(const Binding &binding)
  {
    return std::make_shared<const nmopt::contract::BlockLayout>(
      "alternate_step4_control",
      std::vector<nmopt::contract::SpaceId>{{"control"}},
      std::vector<std::size_t>{binding.problem().control_dimension()});
  }

  std::filesystem::path
  create_binding_artifact();

  void
  run_nmopt_binding_construction()
  {
    Instrumentation instrumentation;
    const auto artifact_root = create_binding_artifact();
    external_dealii_step4_test::EvidenceGuard evidence(
      artifact_root, "nmopt_binding_construction", {{"nmopt", &instrumentation}});
    try
      {
    Binding minimal_binding;
    Vector  minimal_control(minimal_binding.problem().control_dimension());
    minimal_control = 0.0;
    const auto minimal_value = minimal_binding.reduced().evaluate_value(
      make_single_block(minimal_binding.control_layout(), minimal_control));
    require(minimal_value.state_solve.converged(),
            "nmopt binding without diagnostics did not solve");

    Binding         binding(instrumentation);
    const auto point =
      make_off_solution_point(binding.problem().state_dimension());

    require(instrumentation.assembly_calls == 1,
            "nmopt binding did not assemble exactly once");
    require(binding.variable_layout()->n_blocks() == 2,
            "nmopt binding variable layout has the wrong block count");
    require(binding.variable_layout()->space(0).value == "state" &&
              binding.variable_layout()->space(1).value == "control",
            "nmopt binding variable space identities are wrong");
    require(binding.test_layout()->n_blocks() == 1 &&
              binding.test_layout()->space(0).value == "state_test",
            "nmopt binding test layout is wrong");
    require(binding.variable_layout()->dimension(0) ==
              binding.problem().state_dimension(),
            "nmopt binding state dimension is not application-owned");
    require(binding.variable_layout()->dimension(1) ==
              binding.problem().control_dimension(),
            "nmopt binding control dimension is not application-owned");

    const auto alternate_variables = compatible_variable_layout(binding);
    const auto alternate_test      = compatible_test_layout(binding);
    const auto variables =
      make_full_point(alternate_variables, point.state, point.control);
    const auto tangent = make_full_point(alternate_variables,
                                         point.state_tangent,
                                         point.control_tangent);
    const auto seed = make_single_block(alternate_test, point.test_seed);

    const auto native_residual =
      binding.problem().residual(point.state, point.control);
    const auto model_residual = binding.model().residual(variables);
    require(model_residual.layout()->compatible_with(*binding.test_layout()),
            "nmopt residual callback returned the wrong layout");
    require_vector_equal(model_residual.block(0),
                         native_residual,
                         "nmopt residual callback changed the native value");

    const auto native_jvp = binding.problem().residual_jvp(
      point.state_tangent, point.control_tangent);
    const auto model_jvp = binding.model().residual_jvp(variables, tangent);
    require_vector_equal(model_jvp.block(0),
                         native_jvp,
                         "nmopt residual JVP callback changed the native value");

    const auto native_vjp = binding.problem().residual_vjp(point.test_seed);
    const auto model_vjp = binding.model().residual_vjp(variables, seed);
    require(model_vjp.layout()->compatible_with(*binding.variable_layout()),
            "nmopt residual VJP callback returned the wrong layout");
    require_vector_equal(model_vjp.block(0),
                         native_vjp.state,
                         "nmopt residual VJP state block changed the native value");
    require_vector_equal(model_vjp.block(1),
                         native_vjp.control,
                         "nmopt residual VJP control block changed the native value");

    require(binding.model().objective(variables) ==
              binding.problem().objective(point.state, point.control),
            "nmopt objective callback changed the native value");
    const auto native_objective_derivative =
      binding.problem().objective_derivative(point.state, point.control);
    const auto model_objective_derivative =
      binding.model().objective_derivative(variables);
    require_vector_equal(model_objective_derivative.block(0),
                         native_objective_derivative.state,
                         "nmopt objective state derivative changed the native value");
    require_vector_equal(model_objective_derivative.block(1),
                         native_objective_derivative.control,
                         "nmopt objective control derivative changed the native value");

    const auto compatible_control = compatible_control_layout(binding);
    const auto control = make_single_block(compatible_control, point.control);
    const auto metric_covector = binding.metric().apply(control);
    const auto metric_primal = binding.metric().inverse_apply(metric_covector);
    require_vector_equal(metric_primal.block(0),
                         point.control,
                         "Step-4 identity metric changed the control");
    require(metric_covector.layout()->compatible_with(*binding.control_layout()),
            "Step-4 identity metric returned the wrong layout");

    Vector zero_control(binding.problem().control_dimension());
    zero_control = 0.0;
    const auto value = binding.reduced().evaluate_value(
      make_single_block(binding.control_layout(), zero_control));
    require_report(value.state_solve, "nmopt state solve");
    require(value.state.layout()->compatible_with(*binding.state_layout()),
            "nmopt reduced value returned the wrong state layout");

    const auto derivative = binding.reduced().augment_derivative(value);
    require_report(derivative.adjoint_solve, "nmopt adjoint solve");
    require(instrumentation.solve_records.size() == 2,
            "nmopt binding did not record both native solves");
    require(value.state_solve.achieved_residual ==
              instrumentation.solve_records.front().final_residual,
            "nmopt state report changed the native monitored residual");
    require(derivative.adjoint_solve.achieved_residual ==
              instrumentation.solve_records.back().final_residual,
            "nmopt adjoint report changed the native monitored residual");
    const auto metric_direction =
      binding.metric().inverse_apply(derivative.reduced_derivative);
    Vector negative_direction = metric_direction.block(0);
    negative_direction *= -1.0;
    require(negative_direction * derivative.reduced_derivative.block(0) < 0.0,
            "identity metric did not produce a negative gradient direction");

    const auto bad_variables = make_full_point(
      incompatible_variable_layout(binding), point.state, point.control);
    require_rejected([&] { binding.model().objective(bad_variables); },
                     "nmopt model accepted an incompatible variable layout");
    const auto bad_seed = make_single_block(incompatible_test_layout(binding),
                                            point.test_seed);
    require_rejected(
      [&] { binding.model().residual_vjp(variables, bad_seed); },
      "nmopt model accepted an incompatible test layout");

    require(instrumentation.state_solve_calls == 1 &&
              instrumentation.adjoint_solve_calls == 1,
            "nmopt binding solve schedule was not staged");
    require(instrumentation.residual_calls == 2 &&
              instrumentation.residual_jvp_calls == 2 &&
              instrumentation.residual_vjp_calls == 3 &&
              instrumentation.control_vjp_calls == 0,
            "nmopt binding callback counts are inconsistent");
    require(instrumentation.objective_calls == 3 &&
              instrumentation.objective_derivative_calls == 3,
            "nmopt binding objective callback counts are inconsistent");
    require(instrumentation.explicit_matrix_vmult_calls == 4 &&
              instrumentation.explicit_matrix_tvmult_calls == 3,
            "nmopt binding explicit matrix-action counts are inconsistent");
    require(instrumentation.metric_apply_calls == 1 &&
              instrumentation.metric_inverse_apply_calls == 2,
            "nmopt binding metric callback counts are inconsistent");
    require(instrumentation.residual_vjp_calls == 3,
            "nmopt binding did not exercise the full VJP callback");
    evidence.complete();
      }
    catch (...)
      {
        evidence.fail_current_exception();
        throw;
      }
  }
  double
  paired_vector_error(const Vector &left, const Vector &right)
  {
    verification::require_finite(left, "paired vector left value");
    verification::require_finite(right, "paired vector right value");
    return vector_difference(left, right);
  }

  double
  paired_vector_bound(const Vector &left, const Vector &right)
  {
    return 1.0e-11 +
           1.0e-10 * std::max(left.l2_norm(), right.l2_norm());
  }

  double
  paired_scalar_error(const double left, const double right)
  {
    verification::require_finite(left, "paired scalar left value");
    verification::require_finite(right, "paired scalar right value");
    return std::abs(left - right);
  }

  double
  paired_scalar_bound(const double left, const double right)
  {
    return 1.0e-12 +
           1.0e-11 * std::max(std::abs(left), std::abs(right));
  }

  void
  require_paired_vector(const Vector &      actual,
                        const Vector &      expected,
                        const std::string &message)
  {
    const double error = paired_vector_error(actual, expected);
    const double bound = paired_vector_bound(actual, expected);
    require(std::isfinite(error) && std::isfinite(bound) && error <= bound,
            message + ": error=" + std::to_string(error) +
              ", bound=" + std::to_string(bound));
  }

  void
  require_paired_scalar(const double        actual,
                        const double        expected,
                        const std::string &message)
  {
    const double error = paired_scalar_error(actual, expected);
    const double bound = paired_scalar_bound(actual, expected);
    require(std::isfinite(error) && std::isfinite(bound) && error <= bound,
            message + ": error=" + std::to_string(error) +
              ", bound=" + std::to_string(bound));
  }

  double
  matrix_difference(const Matrix &left, const Matrix &right)
  {
    require(left.m() == right.m() && left.n() == right.n(),
            "paired Problem A matrices have incompatible shapes");

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
  create_comparison_artifact(const std::string &prefix = "reduced")
  {
    return external_dealii_step4_test::create_unique_artifact_root(
      find_repository_root() / "runs/external-dealii/step-4/reduced-evaluation",
      prefix);
  }

  std::filesystem::path
  create_binding_artifact()
  {
    return external_dealii_step4_test::create_unique_artifact_root(
      find_repository_root() / "runs/external-dealii/step-4/working/binding",
      "construction");
  }

  void
  require_solve_report_match(const ProblemA::SolveEvidence &native,
                             const nmopt::contract::LinearSolveReport &nmopt,
                             const std::string                       &name)
  {
    require(native.converged, name + " native solve did not converge");
    require_report(nmopt, name + " nmopt solve");
    require(native.iterations == nmopt.iterations,
            name + " changed the CG iteration count");
    require_paired_scalar(native.final_residual,
                          nmopt.achieved_residual,
                          name + " changed the monitored residual");
  }

  struct PairedEvaluation
  {
    std::string label;
    Vector control;
    Vector native_state;
    Vector nmopt_state;
    Vector native_adjoint;
    Vector nmopt_adjoint;
    Vector native_gradient;
    Vector nmopt_gradient;
    double native_objective;
    double nmopt_objective;
    double native_state_monitored_residual;
    double nmopt_state_monitored_residual;
    double native_adjoint_monitored_residual;
    double nmopt_adjoint_monitored_residual;
    double state_error;
    double objective_error;
    double adjoint_error;
    double gradient_error;
  };

  PairedEvaluation
  compare_reduced_evaluation(const std::string &label,
                             const Vector       &control,
                             NativeReduced      &native_reduced,
                             Binding            &nmopt_binding,
                             std::ostream       &evidence,
                             std::ostream       &progress)
  {
    progress << "evaluation " << label << '\n';
    progress.flush();
    const auto native_value = native_reduced.evaluate_value(control);
    const auto native_derivative =
      native_reduced.augment_derivative(native_value);

    const auto nmopt_value = nmopt_binding.reduced().evaluate_value(
      make_single_block(nmopt_binding.control_layout(), control));
    const auto nmopt_derivative =
      nmopt_binding.reduced().augment_derivative(nmopt_value);

    const auto &nmopt_state = nmopt_value.state.block(0);
    const auto &nmopt_adjoint = nmopt_derivative.adjoint.block(0);
    const auto &nmopt_gradient = nmopt_derivative.reduced_derivative.block(0);

    const double state_error =
      paired_vector_error(native_value.state, nmopt_state);
    const double objective_error =
      paired_scalar_error(native_value.objective,
                          nmopt_value.objective_value);
    const double adjoint_error =
      paired_vector_error(native_derivative.adjoint, nmopt_adjoint);
    const double gradient_error =
      paired_vector_error(native_derivative.reduced_derivative,
                          nmopt_gradient);
    evidence << std::setprecision(std::numeric_limits<double>::max_digits10)
             << label << ',' << native_value.objective << ','
             << nmopt_value.objective_value << ',' << state_error << ','
             << objective_error << ',' << adjoint_error << ','
             << gradient_error << ',' << native_value.state_solve.iterations
             << ',' << nmopt_value.state_solve.iterations << ','
             << native_value.state_solve.final_residual << ','
             << nmopt_value.state_solve.achieved_residual << ','
             << native_derivative.adjoint_solve.iterations << ','
             << nmopt_derivative.adjoint_solve.iterations << ','
             << native_derivative.adjoint_solve.final_residual << ','
             << nmopt_derivative.adjoint_solve.achieved_residual << '\n';
    evidence.flush();

    require_paired_vector(native_value.state,
                          nmopt_state,
                          label + " state");
    require_paired_scalar(native_value.objective,
                          nmopt_value.objective_value,
                          label + " objective");
    require_paired_vector(native_derivative.adjoint,
                          nmopt_adjoint,
                          label + " adjoint");
    require_paired_vector(native_derivative.reduced_derivative,
                          nmopt_gradient,
                          label + " reduced gradient");
    require_solve_report_match(native_value.state_solve,
                               nmopt_value.state_solve,
                               label + " state solve");
    require_solve_report_match(native_derivative.adjoint_solve,
                               nmopt_derivative.adjoint_solve,
                               label + " adjoint solve");
    progress << "completed " << label << '\n';
    progress.flush();

    return {label,
            control,
            native_value.state,
            nmopt_state,
            native_derivative.adjoint,
            nmopt_adjoint,
            native_derivative.reduced_derivative,
            nmopt_gradient,
            native_value.objective,
            nmopt_value.objective_value,
            native_value.state_solve.final_residual,
            nmopt_value.state_solve.achieved_residual,
            native_derivative.adjoint_solve.final_residual,
            nmopt_derivative.adjoint_solve.achieved_residual,
            state_error,
            objective_error,
            adjoint_error,
            gradient_error};
  }

  void
  run_nmopt_reduced_comparison()
  {
    Instrumentation native_instrumentation;
    Instrumentation nmopt_instrumentation;
    const auto artifact = create_comparison_artifact();
    external_dealii_step4_test::EvidenceGuard evidence(
      artifact,
      "nmopt_reduced_comparison",
      {{"native", &native_instrumentation},
       {"nmopt", &nmopt_instrumentation}},
      external_dealii_step4_test::EvidenceRetention::retain_on_success);
    try
      {
    std::ofstream output(artifact / "comparison.csv");
    require(static_cast<bool>(output),
            "could not open reduced comparison artifact");
    output << "label,native_objective,nmopt_objective,state_error,"
              "objective_error,adjoint_error,gradient_error,"
              "native_state_iterations,nmopt_state_iterations,"
              "native_state_residual,nmopt_state_residual,"
              "native_adjoint_iterations,nmopt_adjoint_iterations,"
              "native_adjoint_residual,nmopt_adjoint_residual\n";
    std::ofstream progress(artifact / "comparison-progress.txt");
    require(static_cast<bool>(progress),
            "could not open reduced comparison progress artifact");
    progress << "status running\n";

    ProblemA      native_problem(native_instrumentation);
    NativeReduced native_reduced(native_problem, native_instrumentation);
    Binding       nmopt_binding(nmopt_instrumentation);

    require(native_problem.state_dimension() ==
              nmopt_binding.problem().state_dimension(),
            "paired Problem A state dimensions differ");
    require(matrix_difference(native_problem.system_matrix(),
                              nmopt_binding.problem().system_matrix()) == 0.0,
            "paired Problem A matrices differ");
    require(vector_difference(native_problem.system_rhs(),
                              nmopt_binding.problem().system_rhs()) == 0.0,
            "paired Problem A RHS vectors differ");

    const auto controls = external_dealii_step4::scenario::reduced_controls();
    require(controls.size() == 4, "paired Problem A controls are incomplete");

    const std::vector<std::string> labels{
      "zero", "constant", "ramp", "alternating"};
    std::vector<PairedEvaluation> evaluations;
    evaluations.reserve(7);
    for (std::size_t index = 0; index < controls.size(); ++index)
      evaluations.push_back(compare_reduced_evaluation(labels[index],
                                                       controls[index],
                                                       native_reduced,
                                                       nmopt_binding,
                                                       output,
                                                       progress));

    const auto repeated_first = compare_reduced_evaluation(
      "ramp repeat first",
      controls[2],
      native_reduced,
      nmopt_binding,
      output,
      progress);
    const auto intervening = compare_reduced_evaluation(
      "alternating intervening",
      controls[3],
      native_reduced,
      nmopt_binding,
      output,
      progress);
    const auto repeated_last = compare_reduced_evaluation(
      "ramp repeat last",
      controls[2],
      native_reduced,
      nmopt_binding,
      output,
      progress);
    evaluations.push_back(repeated_first);
    evaluations.push_back(intervening);
    evaluations.push_back(repeated_last);

    require_paired_vector(repeated_first.native_state,
                          repeated_last.native_state,
                          "native repeated state");
    require_paired_vector(repeated_first.nmopt_state,
                          repeated_last.nmopt_state,
                          "nmopt repeated state");
    require_paired_scalar(repeated_first.native_objective,
                          repeated_last.native_objective,
                          "native repeated objective");
    require_paired_scalar(repeated_first.nmopt_objective,
                          repeated_last.nmopt_objective,
                          "nmopt repeated objective");
    require_paired_vector(repeated_first.native_gradient,
                          repeated_last.native_gradient,
                          "native repeated gradient");
    require_paired_vector(repeated_first.nmopt_gradient,
                          repeated_last.nmopt_gradient,
                          "nmopt repeated gradient");

    const std::size_t evaluation_count = evaluations.size();
    require(native_instrumentation.assembly_calls == 1 &&
              nmopt_instrumentation.assembly_calls == 1,
            "paired paths did not assemble exactly once each");
    require(native_instrumentation.state_solve_calls == evaluation_count &&
              native_instrumentation.adjoint_solve_calls == evaluation_count &&
              nmopt_instrumentation.state_solve_calls == evaluation_count &&
              nmopt_instrumentation.adjoint_solve_calls == evaluation_count,
            "paired solve counts do not match the staged evaluations");
    require(native_instrumentation.objective_calls == evaluation_count &&
              native_instrumentation.objective_derivative_calls ==
                evaluation_count &&
              nmopt_instrumentation.objective_calls == evaluation_count &&
              nmopt_instrumentation.objective_derivative_calls ==
                evaluation_count,
            "paired objective counts do not match the staged evaluations");
    require(native_instrumentation.control_vjp_calls == evaluation_count &&
              native_instrumentation.residual_vjp_calls == 0,
            "native path did not use only the direct control pullback");
    require(nmopt_instrumentation.control_vjp_calls == 0 &&
              nmopt_instrumentation.residual_vjp_calls == evaluation_count,
            "nmopt path did not use the required full residual VJP");
    require(native_instrumentation.explicit_matrix_tvmult_calls == 0 &&
              nmopt_instrumentation.explicit_matrix_tvmult_calls ==
                evaluation_count,
            "full-VJP transpose work was not isolated to nmopt");
    require(native_instrumentation.metric_apply_calls == 0 &&
              native_instrumentation.metric_inverse_apply_calls == 0 &&
              nmopt_instrumentation.metric_apply_calls == 0 &&
              nmopt_instrumentation.metric_inverse_apply_calls == 0,
            "reduced evaluation unexpectedly applied a metric");

    const auto native_runtime_counts =
      external_dealii_step4::runtime_counter_snapshot(
        native_instrumentation);
    const auto nmopt_runtime_counts =
      external_dealii_step4::runtime_counter_snapshot(nmopt_instrumentation);

    std::ofstream residual_audits(artifact / "residual-audits.csv");
    require(static_cast<bool>(residual_audits),
            "could not open reduced residual audit artifact");
    residual_audits
      << "label,path,state_monitored_residual,state_recomputed_absolute_norm,"
         "state_recomputed_normalized_residual,state_residual_scale,"
         "adjoint_monitored_residual,adjoint_recomputed_absolute_norm,"
         "adjoint_recomputed_normalized_residual,adjoint_residual_scale\n";
    for (const auto &evaluation : evaluations)
      {
        const auto native_state_audit = verification::state_equation_residual(
          native_problem, evaluation.native_state, evaluation.control);
        const auto native_adjoint_audit =
          verification::adjoint_equation_residual(native_problem,
                                                  evaluation.native_adjoint,
                                                  evaluation.native_state);
        const auto nmopt_state_audit = verification::state_equation_residual(
          nmopt_binding.problem(), evaluation.nmopt_state, evaluation.control);
        const auto nmopt_adjoint_audit =
          verification::adjoint_equation_residual(nmopt_binding.problem(),
                                                  evaluation.nmopt_adjoint,
                                                  evaluation.nmopt_state);

        require(native_state_audit.normalized <= 1.0e-10,
                evaluation.label + " native state residual audit failed");
        require(native_adjoint_audit.normalized <= 1.0e-10,
                evaluation.label + " native adjoint residual audit failed");
        require(nmopt_state_audit.normalized <= 1.0e-10,
                evaluation.label + " nmopt state residual audit failed");
        require(nmopt_adjoint_audit.normalized <= 1.0e-10,
                evaluation.label + " nmopt adjoint residual audit failed");

        residual_audits << std::setprecision(
                             std::numeric_limits<double>::max_digits10)
                        << evaluation.label << ",native," <<
          evaluation.native_state_monitored_residual << ','
                        << native_state_audit.absolute_norm << ','
                        << native_state_audit.normalized << ','
                        << native_state_audit.scale << ','
                        << evaluation.native_adjoint_monitored_residual << ','
                        << native_adjoint_audit.absolute_norm << ','
                        << native_adjoint_audit.normalized << ','
                        << native_adjoint_audit.scale << '\n'
                        << evaluation.label << ",nmopt," <<
          evaluation.nmopt_state_monitored_residual << ','
                        << nmopt_state_audit.absolute_norm << ','
                        << nmopt_state_audit.normalized << ','
                        << nmopt_state_audit.scale << ','
                        << evaluation.nmopt_adjoint_monitored_residual << ','
                        << nmopt_adjoint_audit.absolute_norm << ','
                        << nmopt_adjoint_audit.normalized << ','
                        << nmopt_adjoint_audit.scale << '\n';
      }
    residual_audits.flush();
    require(external_dealii_step4::runtime_counter_snapshot(
              native_instrumentation) == native_runtime_counts,
            "native residual audits changed runtime counters");
    require(external_dealii_step4::runtime_counter_snapshot(
              nmopt_instrumentation) == nmopt_runtime_counts,
            "nmopt residual audits changed runtime counters");

    output << "\ncount,native,nmopt\n"
           << "state_solves," << native_instrumentation.state_solve_calls
           << ',' << nmopt_instrumentation.state_solve_calls << '\n'
           << "adjoint_solves," << native_instrumentation.adjoint_solve_calls
           << ',' << nmopt_instrumentation.adjoint_solve_calls << '\n'
           << "control_vjp," << native_instrumentation.control_vjp_calls << ','
           << nmopt_instrumentation.control_vjp_calls << '\n'
           << "residual_vjp," << native_instrumentation.residual_vjp_calls << ','
           << nmopt_instrumentation.residual_vjp_calls << '\n'
           << "explicit_matrix_tvmult,"
           << native_instrumentation.explicit_matrix_tvmult_calls << ','
           << nmopt_instrumentation.explicit_matrix_tvmult_calls << '\n'
           << "metric_apply," << native_instrumentation.metric_apply_calls
           << ',' << nmopt_instrumentation.metric_apply_calls << '\n'
           << "metric_inverse_apply,"
           << native_instrumentation.metric_inverse_apply_calls << ','
           << nmopt_instrumentation.metric_inverse_apply_calls << '\n';
    output.flush();
    progress << "status complete\n";
    progress.flush();
    evidence.complete();

    std::cout << "Step-4 native/nmopt reduced comparison passed: "
              << artifact.lexically_relative(find_repository_root()).generic_string()
              << '\n';
      }
    catch (...)
      {
        evidence.fail_current_exception();
        throw;
      }
  }

  std::filesystem::path
  run_deliberate_comparison_failure(const std::string &prefix)
  {
    Instrumentation instrumentation;
    const auto artifact = create_comparison_artifact(prefix);
    external_dealii_step4_test::EvidenceGuard evidence(
      artifact,
      "failure_evidence",
      {{"native", &instrumentation}},
      external_dealii_step4_test::EvidenceRetention::retain_on_success);
    try
      {
    std::ofstream output(artifact / "comparison.csv");
    require(static_cast<bool>(output),
            "could not open failure comparison artifact");
    output << "label,actual,expected\n";

    ProblemA problem(instrumentation);
    const auto control =
      external_dealii_step4::scenario::reduced_controls().front();
    const auto value = problem.solve_state(control);
    const double actual = problem.objective(value.solution, control);
    const double expected = actual + 1.0;
    output << std::setprecision(std::numeric_limits<double>::max_digits10)
           << "deliberate," << actual << ',' << expected << '\n';
    output.flush();

    bool comparison_failed = false;
    try
      {
        require_paired_scalar(actual, expected, "deliberate comparison failure");
      }
    catch (const std::exception &)
      {
        comparison_failed = true;
      }
    require(comparison_failed,
            "deliberate comparison did not produce a failure");
    return artifact;
      }
    catch (...)
      {
        evidence.fail_current_exception();
        throw;
      }
  }

  void
  run_unhandled_exception_failure(const std::string &       prefix,
                                  std::filesystem::path &artifact)
  {
    Instrumentation instrumentation;
    artifact = create_comparison_artifact(prefix);
    external_dealii_step4_test::EvidenceGuard evidence(
      artifact, "unhandled_exception", {{"native", &instrumentation}});
    try
      {
        ProblemA problem(instrumentation);
        const auto control =
          external_dealii_step4::scenario::reduced_controls().front();
        (void)problem.solve_state(control);
        throw std::runtime_error("unhandled Step-4 exception diagnostic");
      }
    catch (...)
      {
        evidence.fail_current_exception();
        throw;
      }
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
  run_failure_evidence_contract()
  {
    const auto first = run_deliberate_comparison_failure("failure-first");
    const auto retry = run_deliberate_comparison_failure("failure-retry");
    require(first != retry, "failure evidence retry reused the run directory");
    for (const auto &artifact : {first, retry})
      {
        const auto status = read_file(artifact / "status.txt");
        const auto failure = read_file(artifact / "failure.txt");
        const auto comparison = read_file(artifact / "comparison.csv");
        const auto solves = read_file(artifact / "solve-records.csv");
        require(status.find("status failed") != std::string::npos,
                "failure evidence did not retain failed status");
        require(failure.find("deliberate comparison failure") !=
                  std::string::npos,
                "failure evidence did not retain the original diagnostic");
        require(comparison.find("deliberate,") != std::string::npos,
                "failure evidence did not retain compared values");
        require(solves.find("success,state") != std::string::npos,
                "failure evidence did not retain solve evidence");
      }

    std::filesystem::path exception_first;
    std::filesystem::path exception_retry;
    bool                  exception_first_thrown = false;
    try
      {
        run_unhandled_exception_failure("exception-first", exception_first);
      }
    catch (const std::exception &exception)
      {
        exception_first_thrown = true;
        require(exception.what() ==
                  std::string("unhandled Step-4 exception diagnostic"),
                "unhandled exception did not propagate its original message");
      }
    require(exception_first_thrown,
            "unhandled exception probe did not throw");

    bool exception_retry_thrown = false;
    try
      {
        run_unhandled_exception_failure("exception-retry", exception_retry);
      }
    catch (const std::exception &exception)
      {
        exception_retry_thrown = true;
        require(exception.what() ==
                  std::string("unhandled Step-4 exception diagnostic"),
                "exception retry did not propagate its original message");
      }
    require(exception_retry_thrown,
            "unhandled exception retry did not throw");
    require(exception_first != exception_retry,
            "unhandled exception retry reused the run directory");
    for (const auto &artifact : {exception_first, exception_retry})
      {
        const auto status = read_file(artifact / "status.txt");
        const auto failure = read_file(artifact / "failure.txt");
        const auto counters = read_file(artifact / "counters.csv");
        const auto solves = read_file(artifact / "solve-records.csv");
        require(status.find("status failed") != std::string::npos,
                "unhandled exception evidence did not retain failed status");
        require(failure.find("unhandled Step-4 exception diagnostic") !=
                  std::string::npos,
                "unhandled exception evidence lost the original diagnostic");
        require(counters.find("state_solve_calls,1") != std::string::npos,
                "unhandled exception evidence lost counters");
        require(solves.find("success,state") != std::string::npos,
                "unhandled exception evidence lost solve records");
      }
  }
} // namespace

int
main(const int argc, char **argv)
{
  try
    {
      const std::vector<nmopt::test_support::Scenario> scenarios{
        {"nmopt_binding_construction",
         "nmopt.external_tutorial_step_4.nmopt_binding_construction",
         {"dealii", "application", "external", "tutorial", "integration"},
         60,
         run_nmopt_binding_construction},
        {"nmopt_reduced_comparison",
         "nmopt.external.tutorial_step_4.nmopt_reduced_comparison",
         {"dealii", "application", "external", "tutorial", "integration", "extended"},
         180,
         run_nmopt_reduced_comparison},
        {"failure_evidence",
         "nmopt.external.tutorial_step_4.failure_evidence",
         {"dealii", "application", "external", "tutorial", "integration", "diagnostics", "reproduction"},
         60,
         run_failure_evidence_contract}};
      const auto result = nmopt::test_support::run_requested_scenarios(
        argc, argv, scenarios, std::cout);
      if (!result.listed)
        std::cout << "Step-4 nmopt binding scenario passed: " << result.executed
                  << '\n';
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "Step-4 nmopt binding test failed: " << exception.what()
                << '\n';
      return 1;
    }
}
