#include "nmopt/contract/pdas.hpp"
#include "nmopt/dealii/scalar_diffusion_reaction_kkt.hpp"
#include "nmopt/dealii/serial_kkt_solver.hpp"
#include "nmopt/dealii/serial_pdas.hpp"
#include "test_support/scenario_dispatch.hpp"

#include <deal.II/base/function_lib.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <cmath>
#include <exception>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace
{
  using Backend = nmopt::dealii_backend::SerialBackend;
  using KKT = nmopt::dealii_backend::ScalarDiffusionReactionKKT<2>;
  using Product = KKT::Product;
  using Vector = dealii::Vector<double>;
  using Primal = nmopt::contract::PrimalBlockT<Backend>;
  using Covector = nmopt::contract::CovectorBlockT<Backend>;
  using Box = nmopt::dealii_backend::SerialCellwiseBoxComplementarity;

  void
  require_close(const double actual,
                const double expected,
                const double tolerance,
                const std::string &description)
  {
    nmopt::contract::require(std::abs(actual - expected) <= tolerance,
                             description);
  }

  std::shared_ptr<const KKT::Model>
  make_model(dealii::Triangulation<2> &triangulation)
  {
    dealii::Functions::ConstantFunction<2> forcing(1.0);
    dealii::Functions::ConstantFunction<2> desired_state(0.25);
    return std::make_shared<const KKT::Model>(triangulation,
                                              forcing,
                                              desired_state,
                                              1.3,
                                              0.2,
                                              0.7,
                                              1);
  }

  nmopt::contract::QuadraticKKTSolverPolicy
  serial_solver_policy()
  {
    nmopt::contract::QuadraticKKTSolverPolicy policy;
    policy.method = nmopt::contract::QuadraticKKTSolverMethod::minres;
    policy.maximum_iterations = 1000;
    policy.relative_tolerance = 1e-10;
    policy.absolute_tolerance = 1e-12;
    return policy;
  }

  nmopt::contract::PDASPolicy
  pdas_policy()
  {
    nmopt::contract::PDASPolicy policy;
    policy.active_set_assumptions = {
      true,
      true,
      "serial scalar equality plus active control rows have declared rank",
      "serial scalar objective is positive on the restricted equality kernel",
      true,
      true,
      "serial PDAS D-transpose and KKT-transpose actions are declared exact"};
    return policy;
  }

  Product::Point
  zero_point(const Product &product)
  {
    return Product::Point{Product::Primal::zeros(product.layout().primal),
                           Product::Primal::zeros(product.layout().multiplier)};
  }

  Covector
  zero_box_multiplier(const Box &box)
  {
    const auto &layout = box.contract().layout();
    return Covector(layout, {Backend::zeros(layout->dimension(0))});
  }

  void
  test_inactive_serial_pdas_and_metric_conversion()
  {
    dealii::Triangulation<2> triangulation;
    dealii::GridGenerator::hyper_cube(triangulation, 0.0, 1.0);
    triangulation.refine_global(1);
    const auto model = make_model(triangulation);
    const KKT kkt(model);
    const auto metric = std::make_shared<const nmopt::dealii_backend::MassMetric>(
      model->control_l2_metric());
    const Box box(model->control_layout(), -100.0, 100.0, metric);

    const auto &control_layout = model->control_layout();
    Vector control(Backend::checked_native_size(control_layout->dimension(0)));
    for (std::size_t index = 0; index < control.size(); ++index)
      control[Backend::checked_native_size(index)] =
        0.1 * static_cast<double>(index + 1);
    const Primal primal(control_layout, {control});
    const auto multiplier = box.contract().primal_to_multiplier(primal);
    const auto recovered = box.contract().multiplier_to_primal(multiplier);
    Vector conversion_error = recovered.block(0);
    conversion_error.add(-1.0, primal.block(0));
    require_close(conversion_error.l2_norm(),
                  0.0,
                  1e-10,
                  "serial cellwise L2 multiplier conversion");

    const auto initial_multiplier = zero_box_multiplier(box);
    const auto initial_selection = box.contract().classify(
      Primal(control_layout,
             {Backend::zeros(control_layout->dimension(0))}),
      initial_multiplier,
      1.0);
    nmopt::contract::require(initial_selection.active_size() == 0,
                             "broad serial box unexpectedly selected an active set");

    const auto solver_policy = serial_solver_policy();
    const auto solve_action = [&solver_policy](const Product &product) {
      return nmopt::dealii_backend::solve_serial_quadratic_kkt(
        product, solver_policy);
    };
    const nmopt::contract::PDASSolverT<Backend> solver(
      kkt.product(),
      box.contract(),
      1,
      solve_action);
    const auto result = solver.solve(zero_point(kkt.product()),
                                     initial_multiplier,
                                     pdas_policy());

    nmopt::contract::require(result.converged(),
                             "serial inactive-box PDAS did not converge");
    nmopt::contract::require(result.iterations.size() == 1,
                             "serial inactive-box PDAS used too many iterations");
    nmopt::contract::require(
      result.iterations.front().selection.active_size() == 0 &&
        result.iterations.front().active_set_stable &&
        result.iterations.front().kkt_residuals_converged,
      "serial inactive-box PDAS reported incomplete convergence diagnostics");
  }

  void
  test_active_serial_pdas_stabilises_at_bound()
  {
    dealii::Triangulation<2> triangulation;
    dealii::GridGenerator::hyper_cube(triangulation, 0.0, 1.0);
    triangulation.refine_global(1);
    const auto model = make_model(triangulation);
    const KKT kkt(model);
    const auto metric = std::make_shared<const nmopt::dealii_backend::MassMetric>(
      model->control_l2_metric());
    const auto solver_policy = serial_solver_policy();
    const auto solve_action = [&solver_policy](const Product &product) {
      return nmopt::dealii_backend::solve_serial_quadratic_kkt(
        product, solver_policy);
    };

    const Box broad_box(model->control_layout(), -100.0, 100.0, metric);
    const nmopt::contract::PDASSolverT<Backend> unconstrained_solver(
      kkt.product(),
      broad_box.contract(),
      1,
      solve_action);
    const auto unconstrained = unconstrained_solver.solve(
      zero_point(kkt.product()),
      zero_box_multiplier(broad_box),
      pdas_policy());
    nmopt::contract::require(unconstrained.converged(),
                             "serial reference KKT solve did not converge");

    const std::size_t control_dimension = model->control_layout()->dimension(0);
    Vector lower(Backend::checked_native_size(control_dimension));
    Vector upper(Backend::checked_native_size(control_dimension));
    lower = -100.0;
    upper = 100.0;
    const double reference_control = unconstrained.solution.primal.block(1)[0];
    nmopt::contract::require(std::abs(reference_control) > 1e-8,
                             "serial manufactured control did not expose a bound direction");
    if (reference_control < 0.0)
      lower[0] = 0.5 * reference_control;
    else
      upper[0] = 0.5 * reference_control;

    const Box active_box(model->control_layout(),
                         std::move(lower),
                         std::move(upper),
                         metric);
    const nmopt::contract::PDASSolverT<Backend> solver(
      kkt.product(),
      active_box.contract(),
      1,
      solve_action);
    const auto result = solver.solve(zero_point(kkt.product()),
                                     zero_box_multiplier(active_box),
                                     pdas_policy());

    nmopt::contract::require(result.converged(),
                             "serial active-box PDAS did not converge");
    nmopt::contract::require(result.iterations.size() >= 2,
                             "serial active-box PDAS did not re-solve after set change");
    const auto &final_report = result.iterations.back();
    nmopt::contract::require(final_report.active_set_stable &&
                               final_report.selection.active_size() == 1 &&
                               final_report.kkt_residuals_converged,
                             "serial active-box PDAS did not report stable full convergence");
    const double final_control = result.solution.primal.block(1)[0];
    const double bound = reference_control < 0.0
                           ? active_box.contract().bounds().lower().block(0)[0]
                           : active_box.contract().bounds().upper().block(0)[0];
    require_close(final_control,
                  bound,
                  1e-8,
                  "serial active-box PDAS returned the wrong bound value");
  }
} // namespace

int
main(const int argc, char **argv)
{
  try
    {
      const std::vector<nmopt::test_support::Scenario> scenarios{
        {"inactive_serial_pdas_and_metric_conversion",
         "nmopt.dealii.pdas.inactive_serial_pdas_and_metric_conversion",
         {"dealii", "active-set", "complementarity", "pdas"},
         120,
         test_inactive_serial_pdas_and_metric_conversion},
        {"active_serial_pdas_stabilises_at_bound",
         "nmopt.dealii.pdas.active_serial_pdas_stabilises_at_bound",
         {"dealii", "active-set", "kkt", "pdas"},
         120,
         test_active_serial_pdas_stabilises_at_bound}};
      const auto result = nmopt::test_support::run_requested_scenarios(
        argc, argv, scenarios, std::cout);
      if (!result.listed)
        std::cout << "serial deal.II PDAS contract test passed: "
                  << result.executed << '\n';
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "serial deal.II PDAS contract test failed: "
                << exception.what() << '\n';
      return 1;
    }
}
