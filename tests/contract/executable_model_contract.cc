#include "nmopt/contract/callback_executable_model.hpp"
#include "nmopt/contract/reduced_dto.hpp"
#include "../support/contract_errors.hpp"
#include "../support/scenario_dispatch.hpp"

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace
{
  using namespace nmopt::contract;
  using Model = CallbackExecutableModel;
  using Primal = Model::Primal;
  using Covector = Model::Covector;

  struct Fixture
  {
    LayoutPtr variable_layout = std::make_shared<const BlockLayout>(
      "callback_variables",
      std::vector<SpaceId>{{"state"}, {"control"}},
      std::vector<std::size_t>{2, 1});
    LayoutPtr test_layout = std::make_shared<const BlockLayout>(
      "callback_test",
      std::vector<SpaceId>{{"state_test"}},
      std::vector<std::size_t>{2});
  };

  struct Actions
  {
    Model::ResidualAction           residual;
    Model::LinearizedAction         residual_jvp;
    Model::TransposeAction          residual_vjp;
    Model::ObjectiveAction          objective;
    Model::ObjectiveDerivativeAction objective_derivative;
  };

  Actions
  make_actions(const Fixture &fixture)
  {
    const LayoutPtr variable_layout = fixture.variable_layout;
    const LayoutPtr test_layout = fixture.test_layout;

    Actions actions;
    actions.residual = [test_layout](const Primal &variables) {
      DenseVector value = variables.block(0);
      value[0] -= variables.block(1)[0];
      value[1] -= 2.0 * variables.block(1)[0];
      return Covector(test_layout, {std::move(value)});
    };
    actions.residual_jvp = [test_layout](const Primal &, const Primal &tangent) {
      DenseVector value = tangent.block(0);
      value[0] -= tangent.block(1)[0];
      value[1] -= 2.0 * tangent.block(1)[0];
      return Covector(test_layout, {std::move(value)});
    };
    actions.residual_vjp = [variable_layout](const Primal &,
                                             const Primal &test_seed) {
      DenseVector control{
        -test_seed.block(0)[0] - 2.0 * test_seed.block(0)[1]};
      return Covector(variable_layout,
                      {test_seed.block(0), std::move(control)});
    };
    actions.objective = [](const Primal &variables) {
      const DenseVector &state = variables.block(0);
      const double        control = variables.block(1)[0];
      return 0.5 * (state[0] * state[0] + 2.0 * state[1] * state[1]) +
             1.5 * control * control;
    };
    actions.objective_derivative = [variable_layout](const Primal &variables) {
      const DenseVector &state = variables.block(0);
      return Covector(variable_layout,
                      {DenseVector{state[0], 2.0 * state[1]},
                       DenseVector{3.0 * variables.block(1)[0]}});
    };
    return actions;
  }

  Model
  make_model(const Fixture &fixture, Actions actions = {})
  {
    if (!actions.residual)
      actions = make_actions(fixture);
    return Model(fixture.variable_layout,
                 fixture.test_layout,
                 std::move(actions.residual),
                 std::move(actions.residual_jvp),
                 std::move(actions.residual_vjp),
                 std::move(actions.objective),
                 std::move(actions.objective_derivative));
  }

  void
  require_close(const double actual,
                const double expected,
                const char *message,
                const double tolerance = 1e-12)
  {
    require(std::abs(actual - expected) <= tolerance, message);
  }

  void
  require_vector_close(const DenseVector &actual,
                       const DenseVector &expected,
                       const char *       message,
                       const double       tolerance = 1e-12)
  {
    require(actual.size() == expected.size(), message);
    for (std::size_t index = 0; index < actual.size(); ++index)
      require_close(actual[index], expected[index], message, tolerance);
  }

  void
  test_construction_validation()
  {
    const Fixture fixture;

    auto expect_missing = [&](const std::string &expected_message,
                              const auto &       remove_callback,
                              const std::string &description) {
      Actions actions = make_actions(fixture);
      remove_callback(actions);
      nmopt::test_support::require_contract_error(
        [&]() {
          (void)Model(fixture.variable_layout,
                      fixture.test_layout,
                      std::move(actions.residual),
                      std::move(actions.residual_jvp),
                      std::move(actions.residual_vjp),
                      std::move(actions.objective),
                      std::move(actions.objective_derivative));
        },
        expected_message,
        description);
    };

    Actions actions = make_actions(fixture);
    nmopt::test_support::require_contract_error(
      [&]() {
        (void)Model({},
                    fixture.test_layout,
                    std::move(actions.residual),
                    std::move(actions.residual_jvp),
                    std::move(actions.residual_vjp),
                    std::move(actions.objective),
                    std::move(actions.objective_derivative));
      },
      "Callback executable model needs a variable layout",
      "missing variable layout");

    actions = make_actions(fixture);
    nmopt::test_support::require_contract_error(
      [&]() {
        (void)Model(fixture.variable_layout,
                    {},
                    std::move(actions.residual),
                    std::move(actions.residual_jvp),
                    std::move(actions.residual_vjp),
                    std::move(actions.objective),
                    std::move(actions.objective_derivative));
      },
      "Callback executable model needs a test layout",
      "missing test layout");

    expect_missing("Callback executable model needs a residual action",
                   [](Actions &value) { value.residual = {}; },
                   "missing residual callback");
    expect_missing("Callback executable model needs a residual JVP action",
                   [](Actions &value) { value.residual_jvp = {}; },
                   "missing residual JVP callback");
    expect_missing("Callback executable model needs a residual VJP action",
                   [](Actions &value) { value.residual_vjp = {}; },
                   "missing residual VJP callback");
    expect_missing("Callback executable model needs an objective action",
                   [](Actions &value) { value.objective = {}; },
                   "missing objective callback");
    expect_missing(
      "Callback executable model needs an objective derivative action",
      [](Actions &value) { value.objective_derivative = {}; },
      "missing objective derivative callback");
  }

  void
  test_callback_forwarding()
  {
    const Fixture fixture;
    const Model   model = make_model(fixture);
    const Primal point(
      fixture.variable_layout, {DenseVector{0.5, -1.0}, DenseVector{0.25}});
    const Primal tangent(
      fixture.variable_layout, {DenseVector{0.2, 0.3}, DenseVector{-0.4}});
    const Primal test_seed(fixture.test_layout, {DenseVector{1.5, -0.5}});

    const Covector residual = model.residual(point);
    require_vector_close(residual.block(0), DenseVector{0.25, -1.5},
                          "callback residual forwarding");
    require(residual.layout()->compatible_with(*fixture.test_layout),
            "callback residual returned the wrong layout");

    const Covector jvp = model.residual_jvp(point, tangent);
    require_vector_close(jvp.block(0), DenseVector{0.6, 1.1},
                          "callback residual JVP forwarding");

    const Covector vjp = model.residual_vjp(point, test_seed);
    require_vector_close(vjp.block(0), DenseVector{1.5, -0.5},
                          "callback residual VJP state forwarding");
    require_vector_close(vjp.block(1), DenseVector{-0.5},
                          "callback residual VJP control forwarding");

    require_close(model.objective(point), 1.21875,
                  "callback objective forwarding");
    const Covector derivative = model.objective_derivative(point);
    require_vector_close(derivative.block(0), DenseVector{0.5, -2.0},
                          "callback objective derivative state forwarding");
    require_vector_close(derivative.block(1), DenseVector{0.75},
                          "callback objective derivative control forwarding");
  }

  void
  test_layout_validation()
  {
    const Fixture fixture;
    const Model   model = make_model(fixture);
    const LayoutPtr wrong_variables = std::make_shared<const BlockLayout>(
      "wrong_variables",
      std::vector<SpaceId>{{"wrong"}},
      std::vector<std::size_t>{3});
    const LayoutPtr wrong_test = std::make_shared<const BlockLayout>(
      "wrong_test",
      std::vector<SpaceId>{{"wrong_test"}},
      std::vector<std::size_t>{2});

    const Primal wrong_point(wrong_variables, {DenseVector{0.0, 0.0, 0.0}});
    nmopt::test_support::require_contract_error(
      [&]() { (void)model.residual(wrong_point); },
      "Callback executable model residual has an incompatible variable layout",
      "residual input layout validation");

    const Primal point = Primal::zeros(fixture.variable_layout);
    const Primal wrong_tangent(wrong_variables,
                               {DenseVector{0.0, 0.0, 0.0}});
    nmopt::test_support::require_contract_error(
      [&]() { (void)model.residual_jvp(point, wrong_tangent); },
      "Callback executable model residual JVP tangent has an incompatible variable layout",
      "residual JVP tangent layout validation");

    const Primal wrong_seed(wrong_test, {DenseVector{0.0, 0.0}});
    nmopt::test_support::require_contract_error(
      [&]() { (void)model.residual_vjp(point, wrong_seed); },
      "Callback executable model residual VJP seed has an incompatible test layout",
      "residual VJP seed layout validation");

    nmopt::test_support::require_contract_error(
      [&]() { (void)model.objective(wrong_point); },
      "Callback executable model objective has an incompatible variable layout",
      "objective input layout validation");

    Actions bad_residual = make_actions(fixture);
    bad_residual.residual = [wrong_test](const Primal &) {
      return Covector(wrong_test, {DenseVector{0.0, 0.0}});
    };
    const Model bad_residual_model = make_model(fixture, std::move(bad_residual));
    nmopt::test_support::require_contract_error(
      [&]() { (void)bad_residual_model.residual(point); },
      "Callback executable model residual returned an incompatible test layout",
      "residual result layout validation");

    Actions bad_derivative = make_actions(fixture);
    bad_derivative.objective_derivative = [wrong_test](const Primal &) {
      return Covector(wrong_test, {DenseVector{0.0, 0.0}});
    };
    const Model bad_derivative_model =
      make_model(fixture, std::move(bad_derivative));
    nmopt::test_support::require_contract_error(
      [&]() { (void)bad_derivative_model.objective_derivative(point); },
      "Callback executable model objective derivative returned an incompatible variable layout",
      "objective derivative result layout validation");
  }

  void
  test_jvp_vjp_pairing()
  {
    const Fixture fixture;
    const Model   model = make_model(fixture);
    const Primal point = Primal::zeros(fixture.variable_layout);
    const Primal tangent(
      fixture.variable_layout, {DenseVector{0.3, -0.8}, DenseVector{0.7}});
    const Primal test_seed(fixture.test_layout, {DenseVector{1.2, -0.5}});

    require_close(pair(model.residual_jvp(point, tangent), test_seed),
                  pair(model.residual_vjp(point, test_seed), tangent),
                  "callback residual JVP/VJP pairing");
  }

  void
  test_objective_replacement()
  {
    const Fixture fixture;
    Actions       alternate = make_actions(fixture);
    alternate.objective = [](const Primal &variables) {
      return variables.block(0)[0] + 4.0 * variables.block(1)[0];
    };
    alternate.objective_derivative = [layout = fixture.variable_layout](
                                       const Primal &) {
      return Covector(layout, {DenseVector{1.0, 0.0}, DenseVector{4.0}});
    };

    const Model base = make_model(fixture);
    const Model replaced = make_model(fixture, std::move(alternate));
    const Primal point(
      fixture.variable_layout, {DenseVector{0.5, -1.0}, DenseVector{0.25}});

    require_vector_close(base.residual(point).block(0),
                         replaced.residual(point).block(0),
                         "objective replacement changed the residual");
    require_close(base.objective(point), 1.21875,
                  "base objective before replacement");
    require_close(replaced.objective(point), 1.5,
                  "replaced objective forwarding");
    require_vector_close(replaced.objective_derivative(point).block(0),
                         DenseVector{1.0, 0.0},
                         "replaced objective derivative state");
    require_vector_close(replaced.objective_derivative(point).block(1),
                         DenseVector{4.0},
                         "replaced objective derivative control");
  }

  void
  test_reduced_dto_consumption()
  {
    const Fixture fixture;
    const Model   model = make_model(fixture);
    const StateControlPartition partition(model, 0, 1);
    const LayoutPtr              state_layout = partition.state_layout();
    const LayoutPtr              test_layout = model.test_layout();
    using SolveResult = FormulationSolveResultT<DenseBackend>;
    const StateAdjointSolvers    solvers{
      [state_layout](const Primal &control) {
        const double value = control.block(0)[0];
        return SolveResult(
          Primal(state_layout, {DenseVector{value, 2.0 * value}}));
      },
      [test_layout](const Primal &, const Covector &state_rhs) {
        return SolveResult(
          Primal(test_layout, {state_rhs.block(0)}));
      }};
    const ReducedDTO reduced(model, partition, solvers);
    const Primal control(partition.control_layout(), {DenseVector{0.25}});

    const ReducedEvaluation evaluation = reduced.evaluate(control);
    require_vector_close(evaluation.state.block(0),
                         DenseVector{0.25, 0.5},
                         "reduced DTO state solve through callback model");
    require_vector_close(evaluation.full_point.block(0),
                         DenseVector{0.25, 0.5},
                         "reduced DTO full point state through callback model");
    require_close(evaluation.objective_value,
                  0.375,
                  "reduced DTO objective through callback model");
    require_close(evaluation.reduced_derivative.block(0)[0],
                  3.0,
                  "reduced DTO derivative through callback model");
    require(evaluation.state_solve.converged() &&
              evaluation.adjoint_solve.converged(),
            "reduced DTO callback solves did not report convergence");
  }
} // namespace

int
main(const int argc, char **argv)
{
  try
    {
      const std::vector<nmopt::test_support::Scenario> scenarios{
        {"construction_validation",
         "nmopt.executable_model.construction_validation",
         {"backend-neutral", "contract", "executable-model"},
         30,
         test_construction_validation},
        {"callback_forwarding",
         "nmopt.executable_model.callback_forwarding",
         {"backend-neutral", "contract", "executable-model"},
         30,
         test_callback_forwarding},
        {"layout_validation",
         "nmopt.executable_model.layout_validation",
         {"backend-neutral", "contract", "executable-model"},
         30,
         test_layout_validation},
        {"jvp_vjp_pairing",
         "nmopt.executable_model.jvp_vjp_pairing",
         {"backend-neutral", "contract", "executable-model"},
         30,
         test_jvp_vjp_pairing},
        {"objective_replacement",
         "nmopt.executable_model.objective_replacement",
         {"backend-neutral", "contract", "executable-model"},
         30,
         test_objective_replacement},
        {"reduced_dto_consumption",
         "nmopt.executable_model.reduced_dto_consumption",
         {"backend-neutral", "contract", "executable-model", "reduced"},
         30,
         test_reduced_dto_consumption}};
      const auto result = nmopt::test_support::run_requested_scenarios(
        argc, argv, scenarios, std::cout);
      if (!result.listed)
        std::cout << "executable model contract scenario passed: "
                  << result.executed << '\n';
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "executable model contract test failed: "
                << exception.what() << '\n';
      return 1;
    }
}
