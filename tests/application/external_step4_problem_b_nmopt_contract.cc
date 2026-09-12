#define STEP4_NO_MAIN
#include "../../apps/external-dealii/step-4/source/adapted/step-4.cc"
#undef STEP4_NO_MAIN

#include "../../apps/external-dealii/step-4/integration/nmopt_problem_b_binding.hpp"
#include "../../apps/external-dealii/step-4/verification/scenario.hpp"

#include "../dealii/external_step4_evidence.hpp"
#include "../support/scenario_dispatch.hpp"

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <iostream>
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

  std::filesystem::path
  create_binding_artifact()
  {
    auto directory = std::filesystem::current_path();
    while (true)
      {
        if (std::filesystem::exists(
              directory / "apps/external-dealii/step-4/source/upstream/step-4.cc"))
          return external_dealii_step4_test::create_unique_artifact_root(
            directory /
              "runs/external-dealii/step-4/problem-b/reduced-evaluation",
            "binding");

        const auto parent = directory.parent_path();
        if (parent == directory)
          break;
        directory = parent;
      }

    throw std::runtime_error(
      "could not locate the Problem B binding artifact root");
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
         run_problem_b_nmopt_binding_contract}};
      const auto result = nmopt::test_support::run_requested_scenarios(
        argc, argv, scenarios, std::cout);
      if (!result.listed)
        std::cout << "Step-4 Problem B nmopt binding scenario passed: "
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
