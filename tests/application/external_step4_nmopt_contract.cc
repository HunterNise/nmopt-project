#include "../../apps/external-dealii/step-4/evaluation/nmopt_binding.hpp"
#include "../../apps/external-dealii/step-4/evaluation/scenario.hpp"

#include "../support/scenario_dispatch.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
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

  void
  require(const bool condition, const std::string &message)
  {
    if (!condition)
      throw std::runtime_error(message);
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

  void
  run_nmopt_binding_construction()
  {
    Instrumentation instrumentation;
    Binding         binding(instrumentation);
    const auto point = make_off_solution_point(binding.problem().state_dimension());

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
    require(value.state.layout()->compatible_with(
              *binding.state_layout()),
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
    require(instrumentation.residual_vjp_calls == 3,
            "nmopt binding did not exercise the full VJP callback");
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
         run_nmopt_binding_construction}};
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
