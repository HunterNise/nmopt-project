#include "nmopt/dealii/trace_hhalf_metric.hpp"
#include "nmopt/compiler/v1/dealii_dirichlet_control.hpp"
#include "nmopt/compiler/v1/dealii_scalar_diffusion_reaction.hpp"
#include "nmopt/semantic/v1/reference_specs.hpp"

#include "test_support/contract_errors.hpp"
#include "test_support/scenario_dispatch.hpp"

#include <deal.II/base/function_lib.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>

#include <cmath>
#include <exception>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace
{
  using Backend = nmopt::dealii_backend::SerialBackend;
  using Primal = nmopt::contract::PrimalBlockT<Backend>;
  using Covector = nmopt::contract::CovectorBlockT<Backend>;

  template <int dim>
  class LinearDesiredState final : public dealii::Function<dim>
  {
  public:
    double
    value(const dealii::Point<dim> &point,
          const unsigned int        component = 0) const override
    {
      (void)component;
      return point[0] + 2.0 * point[1];
    }

    dealii::Tensor<1, dim>
    gradient(const dealii::Point<dim> &point,
             const unsigned int        component = 0) const override
    {
      (void)point;
      (void)component;
      dealii::Tensor<1, dim> result;
      result[0] = 1.0;
      result[1] = 2.0;
      return result;
    }
  };

  void
  require_close(const double       actual,
                const double       expected,
                const double       tolerance,
                const std::string &description)
  {
    nmopt::contract::require(std::abs(actual - expected) <= tolerance,
                             description + " differs from its oracle");
  }

  void
  require_primal_close(const Primal &     actual,
                       const Primal &     expected,
                       const double       tolerance,
                       const std::string &description)
  {
    nmopt::contract::require_compatible(
      actual, expected, description + " has incompatible layouts");
    dealii::Vector<double> difference = actual.block(0);
    difference.add(-1.0, expected.block(0));
    require_close(difference.l2_norm(), 0.0, tolerance, description);
  }

  void
  require_covector_close(const Covector &   actual,
                         const Covector &   expected,
                         const double       tolerance,
                         const std::string &description)
  {
    nmopt::contract::require_compatible(
      actual, expected, description + " has incompatible layouts");
    for (std::size_t block = 0; block < actual.n_blocks(); ++block)
      {
        dealii::Vector<double> difference = actual.block(block);
        difference.add(-1.0, expected.block(block));
        require_close(difference.l2_norm(), 0.0, tolerance, description);
      }
  }

  void
  require_valid(const nmopt::semantic::v1::ValidationReport &report,
                const std::string &                           description)
  {
    if (report.valid())
      return;
    std::string message = description;
    for (const auto &diagnostic : report.diagnostics())
      message += " {" + diagnostic.component_id + ", " +
                 diagnostic.capability + "}";
    throw nmopt::contract::ContractError(message);
  }

  Primal
  shifted(Primal value, const Primal &direction, const double step)
  {
    nmopt::contract::require_compatible(
      value, direction, "H1/2 test shift has incompatible layouts");
    for (std::size_t block = 0; block < value.n_blocks(); ++block)
      value.add_scaled_block(block, step, direction.block(block));
    return value;
  }

  struct MatrixStorage
  {
    std::shared_ptr<dealii::SparsityPattern>     sparsity;
    std::shared_ptr<dealii::SparseMatrix<double>> matrix;
  };

  MatrixStorage
  volume_h1_matrix()
  {
    dealii::DynamicSparsityPattern dynamic(3, 3);
    for (std::size_t row = 0; row < 3; ++row)
      for (std::size_t column = 0; column < 3; ++column)
        dynamic.add(row, column);
    auto sparsity = std::make_shared<dealii::SparsityPattern>();
    sparsity->copy_from(dynamic);
    auto matrix = std::make_shared<dealii::SparseMatrix<double>>(*sparsity);
    matrix->set(0, 0, 4.0);
    matrix->set(0, 1, 1.0);
    matrix->set(0, 2, 1.0);
    matrix->set(1, 0, 1.0);
    matrix->set(1, 1, 3.0);
    matrix->set(1, 2, 1.0);
    matrix->set(2, 0, 1.0);
    matrix->set(2, 1, 1.0);
    matrix->set(2, 2, 2.0);
    return {std::move(sparsity), std::move(matrix)};
  }

  void
  run_trace_hhalf_metric_contract()
  {
    const auto layout = std::make_shared<const nmopt::contract::BlockLayout>(
      "trace_hhalf", std::vector<nmopt::contract::SpaceId>{{"control"}},
      std::vector<std::size_t>{2});
    nmopt::dealii_backend::MetricSolveParameters solve_parameters;
    solve_parameters.relative_tolerance = 1e-13;
    solve_parameters.absolute_tolerance = 1e-15;
    const auto volume_h1 = volume_h1_matrix();
    const nmopt::dealii_backend::TraceHhalfMetric metric(
      "hhalf_trace", layout, volume_h1.matrix, {0, 2}, solve_parameters);

    dealii::Vector<double> values(2);
    values[0] = 2.0;
    values[1] = -1.0;
    const Primal primal(layout, {values});
    const Covector action = metric.apply(primal);

    // Eliminating volume DoF 1 gives the exact Schur complement
    // [[11/3, 2/3], [2/3, 5/3]].
    require_close(action.block(0)[0],
                  20.0 / 3.0,
                  1e-12,
                  "H1/2 Schur-complement first row");
    require_close(action.block(0)[1],
                  -1.0 / 3.0,
                  1e-12,
                  "H1/2 Schur-complement second row");
    require_primal_close(metric.inverse_apply(action),
                         primal,
                         1e-12,
                         "H1/2 metric round trip");

    dealii::Vector<double> other_values(2);
    other_values[0] = -0.5;
    other_values[1] = 1.25;
    const Primal other(layout, {other_values});
    const Covector other_action = metric.apply(other);
    require_close(nmopt::contract::pair(action, other),
                  nmopt::contract::pair(other_action, primal),
                  1e-12,
                  "H1/2 metric symmetry");
    nmopt::contract::require(nmopt::contract::pair(action, primal) > 0.0,
                             "H1/2 metric is not positive on a nonzero trace");

    auto duplicate_layout =
      std::make_shared<const nmopt::contract::BlockLayout>(
        "duplicate_trace",
        std::vector<nmopt::contract::SpaceId>{{"control"}},
        std::vector<std::size_t>{2});
    nmopt::test_support::require_contract_error(
      [&duplicate_layout, &volume_h1]() {
        (void)nmopt::dealii_backend::TraceHhalfMetric(
          "duplicate", duplicate_layout, volume_h1.matrix, {0, 0});
      },
      "H1/2 trace DoF map contains a duplicate",
      "duplicate H1/2 trace map");
  }

  void
  run_dirichlet_hhalf_model_contract()
  {
    constexpr int dim = 2;
    constexpr double regularisation_weight = 0.3;
    dealii::Triangulation<dim> triangulation;
    dealii::GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(2);
    const dealii::Functions::ZeroFunction<dim> forcing;
    const dealii::Functions::ZeroFunction<dim> desired_state;
    nmopt::dealii_backend::MetricSolveParameters solve_parameters;
    solve_parameters.relative_tolerance = 1e-13;
    solve_parameters.absolute_tolerance = 1e-15;

    using Model =
      nmopt::compiler::v1::detail::DirichletControlLiftingModel<dim>;
    nmopt::compiler::v1::detail::DirichletObjectivePolicy hhalf_policy;
    hhalf_policy.control_norm =
      nmopt::compiler::v1::detail::DirichletControlNormKind::hhalf;
    hhalf_policy.trace_metric_solve = solve_parameters;
    const Model hhalf_model(
      triangulation,
      forcing,
      desired_state,
      1.0,
      0.0,
      regularisation_weight,
      1,
      {0},
      {},
      std::nullopt,
      hhalf_policy);
    const Model l2_model(triangulation,
                         forcing,
                         desired_state,
                         1.0,
                         0.0,
                         regularisation_weight,
                         1,
                         {0});
    nmopt::compiler::v1::detail::DirichletObjectivePolicy mixed_policy;
    mixed_policy.search_metric =
      nmopt::compiler::v1::detail::DirichletControlNormKind::hhalf;
    mixed_policy.trace_metric_solve = solve_parameters;
    const Model l2_loss_hhalf_metric_model(triangulation,
                                            forcing,
                                            desired_state,
                                            1.0,
                                            0.0,
                                            regularisation_weight,
                                            1,
                                            {0},
                                            {},
                                            std::nullopt,
                                            mixed_policy);

    dealii::Vector<double> state(hhalf_model.variable_layout()->dimension(0));
    dealii::Vector<double> control(
      hhalf_model.variable_layout()->dimension(1));
    for (std::size_t index = 0; index < control.size(); ++index)
      control[index] = 0.05 * static_cast<double>(index + 1);
    const Primal point(hhalf_model.variable_layout(), {state, control});
    const Primal control_primal(hhalf_model.control_layout(), {control});

    const auto &hhalf_metric = hhalf_model.control_hhalf_metric();
    const Covector hhalf_action = hhalf_metric.apply(control_primal);
    require_primal_close(hhalf_metric.inverse_apply(hhalf_action),
                         control_primal,
                         1e-11,
                         "assembled Dirichlet H1/2 metric round trip");
    const Covector l2_action =
      l2_model.control_l2_metric(solve_parameters).apply(control_primal);

    const Covector hhalf_derivative = hhalf_model.objective_derivative(point);
    const Covector l2_derivative = l2_model.objective_derivative(point);
    dealii::Vector<double> observed_difference = hhalf_derivative.block(1);
    observed_difference.add(-1.0, l2_derivative.block(1));
    dealii::Vector<double> expected_difference = hhalf_action.block(0);
    expected_difference.add(-1.0, l2_action.block(0));
    expected_difference *= regularisation_weight;
    observed_difference.add(-1.0, expected_difference);
    require_close(observed_difference.l2_norm(),
                  0.0,
                  1e-11,
                  "Dirichlet objective H1/2 regularisation action");

    const double expected_objective_difference =
      0.5 * regularisation_weight *
      (nmopt::contract::pair(hhalf_action, control_primal) -
       nmopt::contract::pair(l2_action, control_primal));
    require_close(hhalf_model.objective(point) - l2_model.objective(point),
                  expected_objective_difference,
                  1e-12,
                  "Dirichlet objective H1/2 regularisation value");
    require_close(l2_loss_hhalf_metric_model.objective(point),
                  l2_model.objective(point),
                  1e-12,
                  "L2 loss remained independent of the H1/2 search metric");
    const Covector mixed_derivative =
      l2_loss_hhalf_metric_model.objective_derivative(point);
    dealii::Vector<double> mixed_derivative_difference =
      mixed_derivative.block(1);
    mixed_derivative_difference.add(-1.0, l2_derivative.block(1));
    require_close(mixed_derivative_difference.l2_norm(),
                  0.0,
                  1e-12,
                  "L2 loss derivative remained independent of the H1/2 search metric");
    const auto &mixed_metric =
      l2_loss_hhalf_metric_model.control_hhalf_metric();
    require_primal_close(mixed_metric.inverse_apply(mixed_metric.apply(
                           control_primal)),
                         control_primal,
                         1e-11,
                         "independent H1/2 search metric round trip");
  }

  void
  run_dirichlet_h1_model_contract()
  {
    constexpr int dim = 2;
    constexpr double regularisation_weight = 0.3;
    dealii::Triangulation<dim> triangulation;
    dealii::GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(2);
    const dealii::Functions::ZeroFunction<dim> forcing;
    const dealii::Functions::ZeroFunction<dim> desired_state;
    nmopt::dealii_backend::MetricSolveParameters solve_parameters;
    solve_parameters.relative_tolerance = 1e-13;
    solve_parameters.absolute_tolerance = 1e-15;

    using Model =
      nmopt::compiler::v1::detail::DirichletControlLiftingModel<dim>;
    nmopt::compiler::v1::detail::DirichletObjectivePolicy h1_policy;
    h1_policy.control_norm =
      nmopt::compiler::v1::detail::DirichletControlNormKind::h1;
    const Model h1_model(
      triangulation,
      forcing,
      desired_state,
      1.0,
      0.0,
      regularisation_weight,
      1,
      {0},
      {},
      std::nullopt,
      h1_policy);
    const Model l2_model(triangulation,
                         forcing,
                         desired_state,
                         1.0,
                         0.0,
                         regularisation_weight,
                         1,
                         {0});

    dealii::Vector<double> state(h1_model.variable_layout()->dimension(0));
    dealii::Vector<double> control(h1_model.variable_layout()->dimension(1));
    for (std::size_t index = 0; index < control.size(); ++index)
      control[index] = 0.05 * static_cast<double>(index + 1);
    const Primal point(h1_model.variable_layout(), {state, control});
    const Primal control_primal(h1_model.control_layout(), {control});

    const auto h1_metric = h1_model.control_h1_metric(solve_parameters);
    const Covector h1_action = h1_metric.apply(control_primal);
    require_primal_close(h1_metric.inverse_apply(h1_action),
                         control_primal,
                         1e-11,
                         "assembled Dirichlet H1 metric round trip");
    const Covector l2_action =
      l2_model.control_l2_metric(solve_parameters).apply(control_primal);
    dealii::Vector<double> stiffness_action = h1_action.block(0);
    stiffness_action.add(-1.0, l2_action.block(0));
    nmopt::contract::require(
      stiffness_action.l2_norm() > 1e-4,
      "Tangential H1 metric did not differ from boundary mass on a nonconstant trace");

    const Covector h1_derivative = h1_model.objective_derivative(point);
    const Covector l2_derivative = l2_model.objective_derivative(point);
    dealii::Vector<double> observed_difference = h1_derivative.block(1);
    observed_difference.add(-1.0, l2_derivative.block(1));
    stiffness_action *= regularisation_weight;
    observed_difference.add(-1.0, stiffness_action);
    require_close(observed_difference.l2_norm(),
                  0.0,
                  1e-11,
                  "Dirichlet objective tangential H1 regularisation action");

    const double expected_objective_difference =
      0.5 * regularisation_weight *
      (nmopt::contract::pair(h1_action, control_primal) -
       nmopt::contract::pair(l2_action, control_primal));
    require_close(h1_model.objective(point) - l2_model.objective(point),
                  expected_objective_difference,
                  1e-12,
                  "Dirichlet objective tangential H1 regularisation value");
  }

  void
  run_dirichlet_h1_tracking_contract()
  {
    constexpr int dim = 2;
    dealii::Triangulation<dim> triangulation;
    dealii::GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(2);
    const dealii::Functions::ZeroFunction<dim> forcing;
    const LinearDesiredState<dim> desired_state;
    using Model =
      nmopt::compiler::v1::detail::DirichletControlLiftingModel<dim>;
    nmopt::compiler::v1::detail::DirichletObjectivePolicy h1_policy;
    h1_policy.state_tracking =
      nmopt::compiler::v1::detail::DirichletStateTrackingNormKind::h1;
    const Model h1_model(triangulation,
                         forcing,
                         desired_state,
                         1.0,
                         0.0,
                         0.3,
                         1,
                         {0},
                         {},
                         std::nullopt,
                         h1_policy);
    const Model l2_model(triangulation,
                         forcing,
                         desired_state,
                         1.0,
                         0.0,
                         0.3,
                         1,
                         {0});

    dealii::Vector<double> state(h1_model.variable_layout()->dimension(0));
    dealii::Vector<double> control(h1_model.variable_layout()->dimension(1));
    dealii::Vector<double> state_direction(state.size());
    dealii::Vector<double> control_direction(control.size());
    for (std::size_t index = 0; index < state.size(); ++index)
      {
        state[index] = 0.02 * static_cast<double>(index + 1);
        state_direction[index] = -0.03 * static_cast<double>(index + 1);
      }
    for (std::size_t index = 0; index < control.size(); ++index)
      {
        control[index] = 0.01 * static_cast<double>(index + 1);
        control_direction[index] =
          (index % 2 == 0 ? 0.02 : -0.015) *
          static_cast<double>(index + 1);
      }
    const Primal point(h1_model.variable_layout(), {state, control});
    const Primal direction(h1_model.variable_layout(),
                           {state_direction, control_direction});

    const Covector derivative = h1_model.objective_derivative(point);
    constexpr double step = 1e-6;
    const double finite_difference =
      (h1_model.objective(shifted(point, direction, step)) -
       h1_model.objective(shifted(point, direction, -step))) /
      (2.0 * step);
    require_close(finite_difference,
                  nmopt::contract::pair(derivative, direction),
                  2e-8,
                  "H1 state-tracking objective derivative");

    const Covector l2_derivative = l2_model.objective_derivative(point);
    dealii::Vector<double> state_difference = derivative.block(0);
    state_difference.add(-1.0, l2_derivative.block(0));
    dealii::Vector<double> control_difference = derivative.block(1);
    control_difference.add(-1.0, l2_derivative.block(1));
    nmopt::contract::require(
      state_difference.l2_norm() + control_difference.l2_norm() > 1e-4 &&
        std::abs(h1_model.objective(point) - l2_model.objective(point)) > 1e-4,
      "H1 state tracking collapsed to the L2 observation geometry");
  }

  void
  run_section_5_11_compilation_contract()
  {
    constexpr int dim = 2;
    dealii::Triangulation<dim> triangulation;
    dealii::GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(2);
    const dealii::Functions::ZeroFunction<dim> forcing;
    const LinearDesiredState<dim> desired_state;
    const nmopt::compiler::v1::DealiiDataBindings<dim> bindings{
      forcing,
      desired_state,
      std::nullopt,
      7.0,
      0.1,
      {"test.section_5_11.forcing", "test.section_5_11.target", ""}};
    nmopt::compiler::v1::DealiiDiscretisationPolicy policy;
    policy.control_metric_solve.relative_tolerance = 1e-13;
    policy.control_metric_solve.absolute_tolerance = 1e-15;
    const nmopt::compiler::v1::DealiiCompiler compiler;

    const auto verify = [&](const nmopt::semantic::v1::ProblemSpec &specification,
                            const std::string &expected_metric_id) {
      require_valid(compiler.validate(specification, policy),
                    specification.id + " did not validate");
      const auto compilation =
        compiler.compile(specification, triangulation, bindings, policy);
      require_valid(compilation.diagnostics,
                    specification.id + " did not compile");
      const auto &model = compilation.problem->executable_model();
      const auto *dirichlet = dynamic_cast<const nmopt::compiler::v1::detail::
        DirichletControlLiftingModel<dim> *>(&model);
      nmopt::contract::require(
        dirichlet != nullptr,
        specification.id + " did not select the Dirichlet lifting target");

      dealii::Vector<double> control_values(
        model.variable_layout()->dimension(1));
      dealii::Vector<double> direction_values(control_values.size());
      for (std::size_t index = 0; index < control_values.size(); ++index)
        {
          control_values[index] = 0.02 * static_cast<double>(index + 1);
          direction_values[index] =
            (index % 2 == 0 ? 0.015 : -0.01) *
            static_cast<double>(index + 1);
        }
      const Primal control(model.variable_layout()->single_block(1, "control"),
                           {control_values});
      const Primal direction(control.layout(), {direction_values});
      const auto reduced = compilation.problem->make_reduced_dto();
      const auto evaluation = reduced.evaluate(control);
      require_close(model.residual(evaluation.full_point).block(0).l2_norm(),
                    0.0,
                    1e-10,
                    specification.id + " state residual");

      const auto &metric = compilation.problem->metric();
      nmopt::contract::require(metric.id() == expected_metric_id,
                               specification.id + " selected the wrong metric");
      require_covector_close(metric.apply(
                               metric.inverse_apply(evaluation.reduced_derivative)),
                             evaluation.reduced_derivative,
                             1e-9,
                             specification.id + " metric round trip");

      const double derivative =
        nmopt::contract::pair(evaluation.reduced_derivative, direction);
      const auto remainder = [&](const double step) {
        return std::abs(
          reduced.evaluate(shifted(control, direction, step)).objective_value -
          evaluation.objective_value - step * derivative);
      };
      const double coarse = remainder(1e-3);
      const double fine = remainder(5e-4);
      nmopt::contract::require(
        coarse > 1e-12 && fine <= 0.27 * coarse + 1e-12,
        specification.id + " reduced Taylor remainder is not quadratic");

      const auto control_loss = std::find_if(
        specification.losses.begin(),
        specification.losses.end(),
        [](const nmopt::semantic::v1::LossSpec &loss) {
          return loss.kind !=
                 nmopt::semantic::v1::LossKind::quadratic_tracking;
        });
      nmopt::contract::require(control_loss != specification.losses.end(),
                               "Section 5.11 graph has no control loss");
      Covector regularisation =
        control_loss->kind == nmopt::semantic::v1::LossKind::
                                quadratic_hhalf_control_regularisation
          ? dirichlet->control_hhalf_metric().apply(control)
        : control_loss->kind == nmopt::semantic::v1::LossKind::
                                quadratic_h1_control_regularisation
          ? dirichlet->control_h1_metric(policy.control_metric_solve)
              .apply(control)
          : dirichlet->control_l2_metric(policy.control_metric_solve)
              .apply(control);
      regularisation.scale_block(0, bindings.regularisation_weight);
      const Covector conormal = dirichlet->discrete_conormal_covector(
        evaluation.full_point, evaluation.adjoint);
      regularisation.add_scaled_block(0, -1.0, conormal.block(0));
      require_covector_close(evaluation.reduced_derivative,
                             regularisation,
                             1e-9,
                             specification.id + " stationarity composition");
    };

    verify(nmopt::semantic::v1::make_hhalf_dirichlet_laplace_control_problem(),
           "hhalf_dirichlet_trace");
    verify(nmopt::semantic::v1::
             make_h1_tracking_hhalf_dirichlet_laplace_control_problem(),
           "hhalf_dirichlet_trace");
    verify(nmopt::semantic::v1::make_h1_dirichlet_laplace_control_problem(),
           "h1_dirichlet_trace");
  }
} // namespace

int
main(const int argc, char **argv)
{
  try
    {
      const std::vector<nmopt::test_support::Scenario> scenarios{
        {"trace_hhalf_metric",
         "nmopt.dealii.trace_hhalf_metric",
         {"dealii", "contract", "metric"},
         30,
         run_trace_hhalf_metric_contract},
        {"dirichlet_hhalf_model",
         "nmopt.dealii.dirichlet_hhalf_model",
         {"dealii", "compiler", "metric"},
         60,
         run_dirichlet_hhalf_model_contract},
        {"dirichlet_h1_model",
         "nmopt.dealii.dirichlet_h1_model",
         {"dealii", "compiler", "metric"},
         60,
         run_dirichlet_h1_model_contract},
        {"dirichlet_h1_tracking",
         "nmopt.dealii.dirichlet_h1_tracking",
         {"dealii", "compiler", "objective"},
         60,
         run_dirichlet_h1_tracking_contract},
        {"section_5_11_compilation",
         "nmopt.dealii.section_5_11_compilation",
         {"dealii", "compiler", "boundary-control"},
         90,
         run_section_5_11_compilation_contract}};
      const auto result = nmopt::test_support::run_requested_scenarios(
        argc, argv, scenarios, std::cout);
      if (!result.listed)
        std::cout << "deal.II H1/2 trace metric contract scenario passed: "
                  << result.executed << '\n';
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "deal.II H1/2 trace metric contract test failed: "
                << exception.what() << '\n';
      return 1;
    }
}
#include <deal.II/base/function_lib.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
