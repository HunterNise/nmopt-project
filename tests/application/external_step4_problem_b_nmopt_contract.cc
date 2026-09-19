#define STEP4_NO_MAIN
#include "../../apps/external-dealii/step-4/source/adapted/step-4.cc"
#undef STEP4_NO_MAIN

#include "../../apps/external-dealii/step-4/evaluation/nmopt_problem_b_binding.hpp"
#include "../../apps/external-dealii/step-4/evaluation/native_problem_b_reduced.hpp"
#include "../../apps/external-dealii/step-4/verification/scenario.hpp"

#include "../dealii/external_step4_evidence.hpp"
#include "../support/scenario_dispatch.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace
{
  using Binding =
    external_dealii_step4::ProblemBNmoptBinding<2, Step4<2>>;
  using Instrumentation = external_dealii_step4::Instrumentation;
  using Primal = Binding::Primal;
  using Covector = Binding::Covector;
  using Vector = external_dealii_step4::ProblemB<2, Step4<2>>::Vector;
  using LayoutPtr = external_dealii_step4::ProblemBNmoptLayoutPtr;
  using NativeProblem = external_dealii_step4::ProblemB<2, Step4<2>>;
  using NativeReduced =
    external_dealii_step4::NativeProblemBReduced<2, Step4<2>>;
  using NativeMetric = external_dealii_step4::ProblemBMetric<2>;
  using Matrix = NativeProblem::Mass::Matrix;

  static_assert(!std::is_copy_constructible_v<Binding>);
  static_assert(!std::is_copy_assignable_v<Binding>);
  static_assert(!std::is_move_constructible_v<Binding>);
  static_assert(!std::is_move_assignable_v<Binding>);

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
            "Problem B public binding vectors have incompatible dimensions");
    Vector difference = left;
    difference.add(-1.0, right);
    return difference.l2_norm();
  }

  void
  require_close(const double       actual,
                const double       expected,
                const double       tolerance,
                const std::string &message)
  {
    require(std::isfinite(actual) && std::isfinite(expected) &&
              std::abs(actual - expected) <= tolerance,
            message);
  }

  void
  require_vector_equal(const Vector &      actual,
                       const Vector &      expected,
                       const std::string &message)
  {
    require(vector_difference(actual, expected) == 0.0, message);
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

  Primal
  make_single_block(const LayoutPtr &layout, const Vector &value)
  {
    return Primal(layout, {value});
  }

  Primal
  make_full_point(const LayoutPtr &layout,
                  const Vector &   state,
                  const Vector &   control)
  {
    return Primal(layout, {state, control});
  }

  Vector
  make_vector(const std::size_t size, const double first, const double step)
  {
    Vector value(size);
    for (unsigned int index = 0; index < value.size(); ++index)
      value[index] = first + step * static_cast<double>(index);
    return value;
  }

  LayoutPtr
  compatible_variable_layout(const Binding &binding)
  {
    return std::make_shared<const nmopt::contract::BlockLayout>(
      "alternate_problem_b_variables",
      std::vector<nmopt::contract::SpaceId>{{"state"}, {"control"}},
      std::vector<std::size_t>{binding.problem().state_dimension(),
                               binding.problem().control_dimension()});
  }

  LayoutPtr
  incompatible_variable_layout(const Binding &binding)
  {
    return std::make_shared<const nmopt::contract::BlockLayout>(
      "incompatible_problem_b_variables",
      std::vector<nmopt::contract::SpaceId>{{"not_state"}, {"control"}},
      std::vector<std::size_t>{binding.problem().state_dimension(),
                               binding.problem().control_dimension()});
  }

  LayoutPtr
  incompatible_control_layout(const Binding &binding)
  {
    return std::make_shared<const nmopt::contract::BlockLayout>(
      "incompatible_problem_b_control",
      std::vector<nmopt::contract::SpaceId>{{"not_control"}},
      std::vector<std::size_t>{binding.problem().control_dimension()});
  }

  LayoutPtr
  incompatible_test_layout(const Binding &binding)
  {
    return std::make_shared<const nmopt::contract::BlockLayout>(
      "incompatible_problem_b_test",
      std::vector<nmopt::contract::SpaceId>{{"not_state_test"}},
      std::vector<std::size_t>{binding.problem().state_dimension()});
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
    require(report.absolute_tolerance == 1.0e-12 &&
              report.requested_tolerance == 1.0e-12,
            name + " reported the wrong tolerance");
    require(std::isfinite(report.achieved_residual) &&
              report.achieved_residual >= 0.0,
              name + " reported an invalid achieved residual");
  }

  double
  paired_vector_error(const Vector &left, const Vector &right)
  {
    require(std::isfinite(left.l2_norm()) && std::isfinite(right.l2_norm()),
            "paired vector contains a non-finite norm");
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
    require(std::isfinite(left) && std::isfinite(right),
            "paired scalar contains a non-finite value");
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
    require(error <= bound,
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
    require(error <= bound,
            message + ": error=" + std::to_string(error) +
              ", bound=" + std::to_string(bound));
  }

  double
  matrix_difference(const Matrix &left, const Matrix &right)
  {
    require(left.m() == right.m() && left.n() == right.n(),
            "paired Problem B matrices have incompatible shapes");

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

  void
  require_solve_match(const NativeProblem::SolveEvidence &native,
                       const nmopt::contract::LinearSolveReport &nmopt,
                       const std::string &name)
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
    Vector      native_state;
    Vector      nmopt_state;
    Vector      native_adjoint;
    Vector      nmopt_adjoint;
    Vector      native_gradient;
    Vector      nmopt_gradient;
    double      native_objective;
    double      nmopt_objective;
  };

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
  create_binding_artifact()
  {
    return external_dealii_step4_test::create_unique_artifact_root(
      find_repository_root() /
        "runs/external-dealii/step-4/problem-b/reduced-evaluation",
      "binding");
  }

  std::filesystem::path
  create_comparison_artifact()
  {
    return external_dealii_step4_test::create_unique_artifact_root(
      find_repository_root() /
        "runs/external-dealii/step-4/problem-b/reduced-evaluation",
      "comparison");
  }

  void
  run_problem_b_nmopt_binding_contract()
  {
    Instrumentation instrumentation;
    const auto artifact_root = create_binding_artifact();
    external_dealii_step4_test::EvidenceGuard evidence(
      artifact_root,
      "problem_b_nmopt_binding",
      {{"nmopt", &instrumentation}});
    try
      {
        Binding uninstrumented;
        Vector  uninstrumented_control(
          uninstrumented.problem().control_dimension());
        uninstrumented_control = 0.0;
        const auto uninstrumented_value = uninstrumented.reduced().evaluate_value(
          make_single_block(uninstrumented.control_layout(),
                            uninstrumented_control));
        require_report(uninstrumented_value.state_solve,
                       "Problem B uninstrumented state solve");

        Binding binding(instrumentation);
        const auto &problem = binding.problem();
        require(problem.state_dimension() == 225,
                "Problem B binding state dimension changed");
        require(problem.control_dimension() == 289,
                "Problem B binding control dimension changed");
        require(binding.variable_layout()->n_blocks() == 2 &&
                  binding.variable_layout()->space(0).value == "state" &&
                  binding.variable_layout()->space(1).value == "control",
                "Problem B binding variable layout is wrong");
        require(binding.test_layout()->n_blocks() == 1 &&
                  binding.test_layout()->space(0).value == "state_test",
                "Problem B binding test layout is wrong");

        const auto alternate_variables = compatible_variable_layout(binding);
        const auto alternate_test =
          std::make_shared<const nmopt::contract::BlockLayout>(
            "alternate_problem_b_test",
            std::vector<nmopt::contract::SpaceId>{{"state_test"}},
            std::vector<std::size_t>{problem.state_dimension()});
        const auto state = make_vector(problem.state_dimension(), 0.1, 0.0003);
        const auto control =
          external_dealii_step4::scenario::ramp_vector();
        const auto state_tangent =
          make_vector(problem.state_dimension(), -0.2, 0.0005);
        const auto control_tangent =
          external_dealii_step4::scenario::alternating_vector();
        const auto test_seed = make_vector(problem.state_dimension(), 0.3, -0.0004);
        const auto variables =
          make_full_point(alternate_variables, state, control);
        const auto tangent =
          make_full_point(alternate_variables, state_tangent, control_tangent);
        const auto seed = make_single_block(alternate_test, test_seed);

        const auto native_residual = problem.residual(state, control);
        const auto model_residual = binding.model().residual(variables);
        require(model_residual.layout()->compatible_with(*binding.test_layout()),
                "Problem B residual callback returned the wrong layout");
        require_vector_equal(model_residual.block(0),
                             native_residual,
                             "Problem B residual callback changed the native value");

        const auto native_jvp =
          problem.residual_jvp(state_tangent, control_tangent);
        const auto model_jvp = binding.model().residual_jvp(variables, tangent);
        require_vector_equal(model_jvp.block(0),
                             native_jvp,
                             "Problem B residual JVP callback changed the native value");

        const auto native_vjp = problem.residual_vjp(test_seed);
        const auto model_vjp = binding.model().residual_vjp(variables, seed);
        require_vector_equal(model_vjp.block(0),
                             native_vjp.state,
                             "Problem B VJP state block changed the native value");
        require_vector_equal(model_vjp.block(1),
                             native_vjp.control,
                             "Problem B VJP control block changed the native value");

        require_close(binding.model().objective(variables),
                      problem.objective(state, control),
                      1.0e-12,
                      "Problem B objective callback changed the native value");
        const auto native_derivative = problem.objective_derivative(state, control);
        const auto model_derivative = binding.model().objective_derivative(variables);
        require_vector_equal(model_derivative.block(0),
                             native_derivative.state,
                             "Problem B objective state derivative changed");
        require_vector_equal(model_derivative.block(1),
                             native_derivative.control,
                             "Problem B objective control derivative changed");

        const auto control_point = make_single_block(binding.control_layout(), control);
        const auto metric_covector = binding.metric().apply(control_point);
        const auto metric_primal = binding.metric().inverse_apply(metric_covector);
        require(vector_difference(metric_primal.block(0), control) <= 1.0e-10,
                "Problem B public metric changed the control");
        require(instrumentation.metric_apply_calls == 1 &&
                  instrumentation.metric_inverse_apply_calls == 1,
                "Problem B public metric calls were not counted");

        Vector zero_control(problem.control_dimension());
        zero_control = 0.0;
        const auto value = binding.reduced().evaluate_value(
          make_single_block(binding.control_layout(), zero_control));
        require_report(value.state_solve, "Problem B public state solve");
        require(value.state.layout()->compatible_with(*binding.state_layout()),
                "Problem B reduced value returned the wrong state layout");
        const auto derivative = binding.reduced().augment_derivative(value);
        require_report(derivative.adjoint_solve,
                       "Problem B public adjoint solve");
        require(derivative.reduced_derivative.layout()->compatible_with(
                  *binding.control_layout()),
                "Problem B reduced derivative returned the wrong layout");

        const auto metric_gradient = binding.reduced().gradient_direction(
          derivative.reduced_derivative, binding.metric());
        require(metric_gradient.layout()->compatible_with(*binding.control_layout()),
                "Problem B public metric gradient returned the wrong layout");
        Vector negative_direction = metric_gradient.block(0);
        negative_direction *= -1.0;
        require(negative_direction * derivative.reduced_derivative.block(0) <
                  0.0,
                "Problem B public metric did not produce a descent direction");

        const auto bad_variables =
          make_full_point(incompatible_variable_layout(binding), state, control);
        require_rejected([&] { binding.model().objective(bad_variables); },
                         "Problem B model accepted an incompatible variable layout");
        const auto bad_seed =
          make_single_block(incompatible_test_layout(binding), test_seed);
        require_rejected(
          [&] { binding.model().residual_vjp(variables, bad_seed); },
          "Problem B model accepted an incompatible test layout");
        const auto bad_control =
          make_single_block(incompatible_control_layout(binding), control);
        require_rejected([&] { binding.metric().apply(bad_control); },
                         "Problem B metric accepted an incompatible layout");

        require(instrumentation.residual_calls == 1 &&
                  instrumentation.residual_jvp_calls == 1 &&
                  instrumentation.residual_vjp_calls == 2 &&
                  instrumentation.objective_calls == 2 &&
                  instrumentation.objective_derivative_calls == 2,
                "Problem B public callback counts are inconsistent");
        require(instrumentation.state_solve_calls == 1 &&
                  instrumentation.adjoint_solve_calls == 1 &&
                  instrumentation.solve_records.size() == 2 &&
                  instrumentation.solve_failures == 0,
                "Problem B public solve schedule is inconsistent");
        evidence.complete();
      }
    catch (...)
      {
        evidence.fail_current_exception();
        throw;
      }
  }

  PairedEvaluation
  compare_reduced_evaluation(const std::string &label,
                             const Vector       &control,
                             NativeReduced      &native_reduced,
                             NativeMetric       &native_metric,
                             Binding            &nmopt_binding,
                             Instrumentation   &native_instrumentation,
                             std::ostream       &output)
  {
    const auto native_value = native_reduced.evaluate_value(control);
    const auto native_derivative =
      native_reduced.augment_derivative(native_value);

    const auto nmopt_value = nmopt_binding.reduced().evaluate_value(
      make_single_block(nmopt_binding.control_layout(), control));
    const auto nmopt_derivative =
      nmopt_binding.reduced().augment_derivative(nmopt_value);

    ++native_instrumentation.metric_inverse_apply_calls;
    const auto native_gradient_result =
      native_metric.inverse_apply(native_derivative.reduced_derivative);
    require(native_gradient_result.evidence.converged,
            label + " native metric inverse did not converge");

    const auto nmopt_gradient = nmopt_binding.reduced().gradient_direction(
      nmopt_derivative.reduced_derivative, nmopt_binding.metric());

    const auto &nmopt_state = nmopt_value.state.block(0);
    const auto &nmopt_adjoint = nmopt_derivative.adjoint.block(0);
    const auto &nmopt_gradient_block = nmopt_gradient.block(0);
    const auto nmopt_full_state =
      nmopt_binding.problem().coordinates().reconstruct(nmopt_state);
    const auto nmopt_full_adjoint =
      nmopt_binding.problem().coordinates().embed_free(nmopt_adjoint);

    require_paired_vector(native_value.state,
                          nmopt_state,
                          label + " state");
    require_paired_vector(native_value.full_state,
                          nmopt_full_state,
                          label + " full state");
    require_paired_vector(native_derivative.adjoint,
                          nmopt_adjoint,
                          label + " adjoint");
    require_paired_vector(native_derivative.full_adjoint,
                          nmopt_full_adjoint,
                          label + " full adjoint");
    require_paired_vector(native_derivative.reduced_derivative,
                          nmopt_derivative.reduced_derivative.block(0),
                          label + " reduced derivative");
    require_paired_vector(native_gradient_result.solution,
                          nmopt_gradient_block,
                          label + " metric gradient");
    require_paired_scalar(native_value.objective,
                          nmopt_value.objective_value,
                          label + " objective");
    require_solve_match(native_value.state_solve,
                        nmopt_value.state_solve,
                        label + " state solve");
    require_solve_match(native_derivative.adjoint_solve,
                        nmopt_derivative.adjoint_solve,
                        label + " adjoint solve");

    output << std::setprecision(std::numeric_limits<double>::max_digits10)
           << label << ',' << native_value.objective << ','
           << nmopt_value.objective_value << ','
           << paired_vector_error(native_value.state, nmopt_state) << ','
           << paired_vector_error(native_value.full_state, nmopt_full_state)
           << ','
           << paired_vector_error(native_derivative.adjoint, nmopt_adjoint)
           << ','
           << paired_vector_error(native_derivative.full_adjoint,
                                  nmopt_full_adjoint)
           << ','
           << paired_vector_error(native_derivative.reduced_derivative,
                                  nmopt_derivative.reduced_derivative.block(0))
           << ','
           << paired_vector_error(native_gradient_result.solution,
                                  nmopt_gradient_block)
           << ',' << native_value.state_solve.iterations << ','
           << nmopt_value.state_solve.iterations << ','
           << native_value.state_solve.final_residual << ','
           << nmopt_value.state_solve.achieved_residual << ','
           << native_derivative.adjoint_solve.iterations << ','
           << nmopt_derivative.adjoint_solve.iterations << ','
           << native_derivative.adjoint_solve.final_residual << ','
           << nmopt_derivative.adjoint_solve.achieved_residual << '\n';

    return {label,
            native_value.state,
            nmopt_state,
            native_derivative.adjoint,
            nmopt_adjoint,
            native_derivative.reduced_derivative,
            nmopt_derivative.reduced_derivative.block(0),
            native_value.objective,
            nmopt_value.objective_value};
  }

  void
  run_problem_b_nmopt_reduced_comparison()
  {
    Instrumentation native_instrumentation;
    Instrumentation nmopt_instrumentation;
    const auto artifact_root = create_comparison_artifact();
    external_dealii_step4_test::EvidenceGuard evidence(
      artifact_root,
      "problem_b_nmopt_reduced_comparison",
      { {"native", &native_instrumentation},
        {"nmopt", &nmopt_instrumentation} },
      external_dealii_step4_test::EvidenceRetention::retain_on_success);
    try
      {
        std::ofstream output(artifact_root / "comparison.csv");
        require(static_cast<bool>(output),
                "could not open Problem B comparison artifact");
        output
          << "label,native_objective,nmopt_objective,state_error,"
             "full_state_error,adjoint_error,full_adjoint_error,"
             "reduced_gradient_error,metric_gradient_error,"
             "native_state_iterations,nmopt_state_iterations,"
             "native_state_residual,nmopt_state_residual,"
             "native_adjoint_iterations,nmopt_adjoint_iterations,"
             "native_adjoint_residual,nmopt_adjoint_residual\n";

        Step4<2> native_tutorial;
        native_tutorial.prepare_for_external_use();
        NativeProblem native_problem(native_tutorial);
        NativeReduced native_reduced(native_problem, native_instrumentation);
        NativeMetric native_metric(native_problem.mass());
        Binding nmopt_binding(nmopt_instrumentation);
        const auto &public_problem = nmopt_binding.problem();

        require(native_problem.state_dimension() ==
                  public_problem.state_dimension(),
                "paired Problem B state dimensions differ");
        require(native_problem.control_dimension() ==
                  public_problem.control_dimension(),
                "paired Problem B control dimensions differ");
        require(matrix_difference(native_problem.mass().mass_matrix(),
                                  public_problem.mass().mass_matrix()) == 0.0,
                "paired Problem B mass matrices differ");
        require(matrix_difference(native_problem.mass().coupling_matrix(),
                                  public_problem.mass().coupling_matrix()) == 0.0,
                "paired Problem B coupling matrices differ");
        require_vector_equal(native_problem.free_system_rhs(),
                             public_problem.free_system_rhs(),
                             "paired Problem B free RHS vectors differ");
        require(native_problem.coordinates().free_indices() ==
                  public_problem.coordinates().free_indices(),
                "paired Problem B free coordinate maps differ");
        require_vector_equal(native_problem.coordinates().lifting(),
                             public_problem.coordinates().lifting(),
                             "paired Problem B coordinate liftings differ");

        const auto controls = external_dealii_step4::scenario::reduced_controls();
        require(controls.size() == 4,
                "Problem B comparison controls are incomplete");
        const std::vector<std::string> labels{
          "zero", "constant", "ramp", "alternating", "ramp_repeat"};
        std::vector<PairedEvaluation> evaluations;
        evaluations.reserve(labels.size());
        for (std::size_t index = 0; index < controls.size(); ++index)
          evaluations.push_back(compare_reduced_evaluation(
            labels[index],
            controls[index],
            native_reduced,
            native_metric,
            nmopt_binding,
            native_instrumentation,
            output));
        evaluations.push_back(compare_reduced_evaluation(
          labels.back(),
          controls[2],
          native_reduced,
          native_metric,
          nmopt_binding,
          native_instrumentation,
          output));

        const auto &first_ramp = evaluations[2];
        const auto &repeated_ramp = evaluations.back();
        require_paired_vector(first_ramp.native_state,
                              repeated_ramp.native_state,
                              "native repeated ramp state");
        require_paired_vector(first_ramp.nmopt_state,
                              repeated_ramp.nmopt_state,
                              "nmopt repeated ramp state");
        require_paired_vector(first_ramp.native_adjoint,
                              repeated_ramp.native_adjoint,
                              "native repeated ramp adjoint");
        require_paired_vector(first_ramp.nmopt_adjoint,
                              repeated_ramp.nmopt_adjoint,
                              "nmopt repeated ramp adjoint");
        require_paired_vector(first_ramp.native_gradient,
                              repeated_ramp.native_gradient,
                              "native repeated ramp gradient");
        require_paired_vector(first_ramp.nmopt_gradient,
                              repeated_ramp.nmopt_gradient,
                              "nmopt repeated ramp gradient");
        require_paired_scalar(first_ramp.native_objective,
                              repeated_ramp.native_objective,
                              "native repeated ramp objective");
        require_paired_scalar(first_ramp.nmopt_objective,
                              repeated_ramp.nmopt_objective,
                              "nmopt repeated ramp objective");

        constexpr std::size_t evaluation_count = 5;
        require(native_instrumentation.state_solve_calls == evaluation_count &&
                  native_instrumentation.adjoint_solve_calls ==
                    evaluation_count &&
                  native_instrumentation.objective_calls == evaluation_count &&
                  native_instrumentation.objective_derivative_calls ==
                    evaluation_count &&
                  native_instrumentation.control_vjp_calls == evaluation_count &&
                  native_instrumentation.residual_calls == 0 &&
                  native_instrumentation.residual_jvp_calls == 0 &&
                  native_instrumentation.residual_vjp_calls == 0 &&
                  native_instrumentation.metric_inverse_apply_calls ==
                    evaluation_count &&
                  native_instrumentation.metric_apply_calls == 0 &&
                  native_instrumentation.solve_failures == 0 &&
                  native_instrumentation.solve_records.size() ==
                    2 * evaluation_count,
                "native Problem B comparison schedule is inconsistent");
        require(nmopt_instrumentation.state_solve_calls == evaluation_count &&
                  nmopt_instrumentation.adjoint_solve_calls ==
                    evaluation_count &&
                  nmopt_instrumentation.objective_calls == evaluation_count &&
                  nmopt_instrumentation.objective_derivative_calls ==
                    evaluation_count &&
                  nmopt_instrumentation.control_vjp_calls == 0 &&
                  nmopt_instrumentation.residual_calls == 0 &&
                  nmopt_instrumentation.residual_jvp_calls == 0 &&
                  nmopt_instrumentation.residual_vjp_calls ==
                    evaluation_count &&
                  nmopt_instrumentation.metric_inverse_apply_calls ==
                    evaluation_count &&
                  nmopt_instrumentation.metric_apply_calls == 0 &&
                  nmopt_instrumentation.solve_failures == 0 &&
                  nmopt_instrumentation.solve_records.size() ==
                    2 * evaluation_count,
                "nmopt Problem B comparison schedule is inconsistent");
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
        {"problem_b_nmopt_binding",
         "nmopt.external_tutorial_step_4.problem_b_nmopt_binding",
         {"dealii", "application", "external", "tutorial", "integration",
          "problem_b"},
         180,
         run_problem_b_nmopt_binding_contract},
        {"problem_b_nmopt_reduced_comparison",
         "nmopt.external_tutorial_step_4.problem_b_nmopt_reduced_comparison",
         {"dealii", "application", "external", "tutorial", "integration",
          "problem_b", "verification", "extended"},
         300,
         run_problem_b_nmopt_reduced_comparison}};
      const auto result = nmopt::test_support::run_requested_scenarios(
        argc, argv, scenarios, std::cout);
      if (!result.listed)
        std::cout << "Step-4 Problem B nmopt contract scenarios passed: "
                  << result.executed << '\n';
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "Step-4 Problem B nmopt binding test failed: "
                << exception.what() << '\n';
      return 1;
    }
}
