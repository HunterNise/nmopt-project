#include "../dealii/external_poisson_fixture.hpp"
#include "nmopt/contract/callback_executable_model.hpp"
#include "nmopt/contract/reduced_dto.hpp"
#include "nmopt/dealii/mass_metric.hpp"
#include "nmopt/dealii/serial_backend.hpp"
#include "nmopt/solvers/reduced_gradient.hpp"
#include "../support/scenario_dispatch.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{
  using Application = external_poisson_application::PoissonControlApplication;
  using Backend = nmopt::dealii_backend::SerialBackend;
  using Model = nmopt::contract::CallbackExecutableModelT<Backend>;
  using Primal = nmopt::contract::PrimalBlockT<Backend>;
  using Covector = nmopt::contract::CovectorBlockT<Backend>;
  using LayoutPtr = nmopt::contract::LayoutPtr;
  using Partition = nmopt::contract::StateControlPartitionT<Backend>;
  using Solvers = nmopt::contract::StateAdjointSolversT<Backend>;
  using Reduced = nmopt::contract::ReducedDTOT<Backend>;
  using Vector = Application::Vector;

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
    require(std::isfinite(actual) && std::isfinite(expected) &&
              std::abs(actual - expected) <= tolerance,
            message);
  }

  void
  require_vector_close(const Vector &      actual,
                       const Vector &      expected,
                       const double       tolerance,
                       const std::string &message)
  {
    require(actual.size() == expected.size(),
            message + " has incompatible dimensions");
    Vector difference = actual;
    difference.add(-1.0, expected);
    require(difference.l2_norm() <= tolerance, message);
  }

  struct Adapter final
  {
    Application application;
    LayoutPtr   variable_layout;
    LayoutPtr   test_layout;
    Model       model;
    Partition   partition;
    Solvers     solvers;
    nmopt::dealii_backend::MassMetric metric;
    Reduced     reduced;

    Adapter()
      : application()
      , variable_layout(std::make_shared<const nmopt::contract::BlockLayout>(
          "external_poisson_variables",
          std::vector<nmopt::contract::SpaceId>{{"state"}, {"control"}},
          std::vector<std::size_t>{application.state_dimension(),
                                   application.control_dimension()}))
      , test_layout(std::make_shared<const nmopt::contract::BlockLayout>(
          "external_poisson_test",
          std::vector<nmopt::contract::SpaceId>{{"state_test"}},
          std::vector<std::size_t>{application.state_dimension()}))
      , model(make_model(application, variable_layout, test_layout, false))
      , partition(model, 0, 1)
      , solvers(make_solvers(application, partition.state_layout(), test_layout))
      , metric("external_poisson_control_l2",
               partition.control_layout(),
               application.control_mass_matrix())
      , reduced(model, partition, solvers)
    {}

    static Model
    make_model(Application &application,
               const LayoutPtr &variable_layout,
               const LayoutPtr &test_layout,
               const bool       alternate_objective)
    {
      Application *const application_ptr = &application;

      return Model(
        variable_layout,
        test_layout,
        [application_ptr, test_layout](const Primal &variables) {
          Vector value = application_ptr->residual(variables.block(0),
                                                   variables.block(1));
          return Covector(test_layout, {std::move(value)});
        },
        [application_ptr, test_layout](const Primal &,
                                       const Primal &variable_tangent) {
          Vector value = application_ptr->residual_jvp(
            variable_tangent.block(0), variable_tangent.block(1));
          return Covector(test_layout, {std::move(value)});
        },
        [application_ptr, variable_layout](const Primal &,
                                           const Primal &test_seed) {
          Vector packed = application_ptr->residual_vjp(test_seed.block(0));
          Vector state(application_ptr->state_dimension());
          Vector control(application_ptr->control_dimension());
          for (std::size_t index = 0; index < state.size(); ++index)
            {
              state[index] = packed[index];
              control[index] =
                packed[application_ptr->state_dimension() + index];
            }
          return Covector(variable_layout,
                          {std::move(state), std::move(control)});
        },
        [application_ptr, alternate_objective](const Primal &variables) {
          if (!alternate_objective)
            return application_ptr->objective(variables.block(0),
                                              variables.block(1));
          return 0.5 * (variables.block(1) * variables.block(1));
        },
        [application_ptr, variable_layout, alternate_objective](
          const Primal &variables) {
          if (!alternate_objective)
            {
              const auto derivative = application_ptr->objective_derivative(
                variables.block(0), variables.block(1));
              return Covector(variable_layout,
                              {std::move(derivative.state),
                               std::move(derivative.control)});
            }

          Vector state(application_ptr->state_dimension());
          state = 0.0;
          Vector control = variables.block(1);
          return Covector(variable_layout,
                          {std::move(state), std::move(control)});
        });
    }

    static Solvers
    make_solvers(Application &application,
                 const LayoutPtr &state_layout,
                 const LayoutPtr &test_layout)
    {
      Application *const application_ptr = &application;
      Solvers             solvers;
      solvers.solve_state = [application_ptr, state_layout](
                              const Primal &control) {
        Vector state = application_ptr->solve_state(control.block(0));
        return nmopt::contract::FormulationSolveResultT<Backend>(
          Primal(state_layout, {std::move(state)}));
      };
      solvers.solve_adjoint = [application_ptr, test_layout](
                                const Primal &,
                                const Covector &state_rhs) {
        Vector adjoint = application_ptr->solve_adjoint(state_rhs.block(0));
        return nmopt::contract::FormulationSolveResultT<Backend>(
          Primal(test_layout, {std::move(adjoint)}));
      };
      return solvers;
    }
  };

  Primal
  zero_control(const Adapter &adapter)
  {
    Vector control(adapter.application.control_dimension());
    control = 0.0;
    return Primal(adapter.partition.control_layout(), {std::move(control)});
  }

  Primal
  make_point(const Adapter &adapter)
  {
    Vector state(adapter.application.state_dimension());
    Vector control(adapter.application.control_dimension());
    for (std::size_t index = 0; index < state.size(); ++index)
      {
        state[index] = 0.1 + 0.013 * static_cast<double>(index);
        control[index] = -0.04 + 0.009 * static_cast<double>(index);
      }
    return Primal(adapter.variable_layout,
                  {std::move(state), std::move(control)});
  }

  Primal
  make_tangent(const Adapter &adapter)
  {
    Vector state(adapter.application.state_dimension());
    Vector control(adapter.application.control_dimension());
    for (std::size_t index = 0; index < state.size(); ++index)
      {
        state[index] = -0.03 + 0.004 * static_cast<double>(index);
        control[index] = 0.02 - 0.003 * static_cast<double>(index);
      }
    return Primal(adapter.variable_layout,
                  {std::move(state), std::move(control)});
  }

  Primal
  make_test_seed(const Adapter &adapter)
  {
    Vector seed(adapter.application.state_dimension());
    for (std::size_t index = 0; index < seed.size(); ++index)
      seed[index] = -0.15 + 0.017 * static_cast<double>(index);
    return Primal(adapter.test_layout, {std::move(seed)});
  }

  Primal
  add_scaled(const Primal &value, const double factor, const Primal &direction)
  {
    Primal result = value;
    nmopt::solvers::add_scaled_primal(result, factor, direction);
    return result;
  }

  void
  test_callback_residual_derivatives()
  {
    Adapter adapter;
    const Primal control = zero_control(adapter);
    const auto state_evaluation = adapter.reduced.evaluate(control);
    require(adapter.application.residual_norm(
               state_evaluation.state.block(0), control.block(0)) < 1e-10,
            "callback state solve left a large native residual");

    const Primal point = make_point(adapter);
    const Primal tangent = make_tangent(adapter);
    const Primal test_seed = make_test_seed(adapter);
    const double finite_difference_step = 1e-6;
    const Covector jvp = adapter.model.residual_jvp(point, tangent);
    const auto residual_plus = adapter.model.residual(
      add_scaled(point, finite_difference_step, tangent));
    const auto residual_minus = adapter.model.residual(
      add_scaled(point, -finite_difference_step, tangent));
    Covector finite_difference = residual_plus;
    finite_difference.add_scaled_block(0, -1.0, residual_minus.block(0));
    finite_difference.scale_block(0, 1.0 / (2.0 * finite_difference_step));
    require_vector_close(finite_difference.block(0),
                         jvp.block(0),
                         1e-8,
                         "callback residual JVP finite difference");

    const Covector vjp = adapter.model.residual_vjp(point, test_seed);
    const double jvp_pairing = nmopt::contract::pair(jvp, test_seed);
    const double vjp_pairing = nmopt::contract::pair(vjp, tangent);
    require_close(jvp_pairing,
                  vjp_pairing,
                  1e-10 * std::max({1.0, std::abs(jvp_pairing),
                                    std::abs(vjp_pairing)}),
                  "callback residual JVP/VJP pairing");
  }

  void
  test_reduced_taylor_check()
  {
    Adapter adapter;
    const Primal control = zero_control(adapter);
    const Primal direction = [&] {
      Vector values(adapter.application.control_dimension());
      for (std::size_t index = 0; index < values.size(); ++index)
        values[index] = 0.015 - 0.002 * static_cast<double>(index);
      return Primal(adapter.partition.control_layout(), {std::move(values)});
    }();

    const auto base = adapter.reduced.evaluate(control);
    const double derivative =
      nmopt::contract::pair(base.reduced_derivative, direction);
    const double step = 1e-4;
    const auto plus = adapter.reduced.evaluate(
      add_scaled(control, step, direction));
    const auto minus = adapter.reduced.evaluate(
      add_scaled(control, -step, direction));
    const double centered_derivative =
      (plus.objective_value - minus.objective_value) / (2.0 * step);
    require_close(centered_derivative,
                  derivative,
                  1e-8,
                  "reduced derivative centered finite difference");

    const double taylor_error =
      std::abs(plus.objective_value - base.objective_value - step * derivative);
    require(taylor_error < 1e-7,
            "reduced state-recomputed Taylor remainder is too large");
    require(adapter.application.residual_norm(
               plus.state.block(0), plus.full_point.block(1)) < 1e-10,
            "perturbed reduced state solve left a large native residual");
  }

  void
  test_optimization_converges()
  {
    Adapter adapter;
    nmopt::solvers::ReducedSolverParameters parameters;
    parameters.maximum_iterations = 100;
    parameters.maximum_line_search_trials = 25;
    parameters.gradient_tolerance = 1e-8;
    nmopt::solvers::ReducedGradientSolverT<Backend> solver(
      adapter.reduced, adapter.metric, parameters);
    const auto result = solver.solve(zero_control(adapter));

    require(result.stopping_reason ==
                  nmopt::solvers::ReducedStoppingReason::gradient_tolerance,
            "external application optimizer did not reach gradient tolerance");
    require(result.accepted_iterations > 0,
            "external application optimizer accepted no iterations");
    require(result.objective_history.back() < result.objective_history.front(),
            "external application optimizer did not reduce the objective");
    require(!result.gradient_norm_history.empty() &&
              result.gradient_norm_history.back() <= 1e-8,
            "external application optimizer final gradient is too large");
  }

  void
  test_objective_replacement_and_output()
  {
    Adapter adapter;
    const Primal control = zero_control(adapter);
    const auto default_value = adapter.reduced.evaluate(control);

    Model alternate_model = Adapter::make_model(
      adapter.application, adapter.variable_layout, adapter.test_layout, true);
    Partition alternate_partition(alternate_model, 0, 1);
    Solvers alternate_solvers = Adapter::make_solvers(
      adapter.application, alternate_partition.state_layout(), adapter.test_layout);
    Reduced alternate_reduced(alternate_model,
                              std::move(alternate_partition),
                              alternate_solvers);
    const auto alternate_value = alternate_reduced.evaluate(control);

    require(std::abs(default_value.objective_value -
                     alternate_value.objective_value) > 1e-8,
            "replacing objective callbacks did not change objective behavior");

    const Primal point = make_point(adapter);
    const Covector default_residual = adapter.model.residual(point);
    const Covector alternate_residual = alternate_model.residual(point);
    require_vector_close(default_residual.block(0),
                         alternate_residual.block(0),
                         1e-12,
                         "objective replacement changed PDE residual callbacks");

    const std::filesystem::path output_directory =
      std::filesystem::temp_directory_path() /
      "external-poisson-adapter-output";
    std::filesystem::remove_all(output_directory);
    adapter.application.write_native_output(output_directory,
                                            default_value.state.block(0));
    const auto output_file = output_directory / "fields-volume.vtu";
    require(std::filesystem::exists(output_file) &&
              std::filesystem::file_size(output_file) > 0,
            "external application native output was not produced");
    std::filesystem::remove_all(output_directory);
  }
} // namespace

int
main(const int argc, char **argv)
{
  try
    {
      const std::vector<nmopt::test_support::Scenario> scenarios{
        {"callback_residual_derivatives",
         "nmopt.external_application.callback_residual_derivatives",
         {"dealii", "application", "external", "derivatives"},
         60,
         test_callback_residual_derivatives},
        {"reduced_taylor_check",
         "nmopt.external_application.reduced_taylor_check",
         {"dealii", "application", "external", "reduced"},
         60,
         test_reduced_taylor_check},
        {"optimization_converges",
         "nmopt.external_application.optimization_converges",
         {"dealii", "application", "external", "optimizer"},
         60,
         test_optimization_converges},
        {"objective_replacement_and_output",
         "nmopt.external_application.objective_replacement_and_output",
         {"dealii", "application", "external", "boundary"},
         60,
         test_objective_replacement_and_output}};
      const auto result = nmopt::test_support::run_requested_scenarios(
        argc, argv, scenarios, std::cout);
      if (!result.listed)
        std::cout << "external application deal.II adapter scenarios passed: "
                  << result.executed << '\n';
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "external application deal.II adapter test failed: "
                << exception.what() << '\n';
      return 1;
    }
}
