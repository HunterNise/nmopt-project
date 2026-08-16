#include "nmopt/contract/supplied_otd.hpp"
#include "nmopt/compiler/v1/compiled_problem.hpp"
#include "test_support/contract_errors.hpp"
#include "test_support/scenario_dispatch.hpp"

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace
{
  using namespace nmopt::contract;
  using System = SuppliedOTDSystem;
  using Primal = System::Primal;
  using Covector = System::Covector;

  struct TestOperators
  {
    LayoutPtr variable_layout = std::make_shared<const BlockLayout>(
      "otd_variables",
      std::vector<SpaceId>{{"state"}, {"adjoint"}, {"control"}},
      std::vector<std::size_t>{2, 2, 1});
    LayoutPtr residual_layout = std::make_shared<const BlockLayout>(
      "otd_residuals",
      std::vector<SpaceId>{{"state_equation"},
                           {"adjoint_equation"},
                           {"control_stationarity"}},
      std::vector<std::size_t>{2, 2, 1});

    DenseMatrix state_state{2, 2, {2.0, 1.0, 0.0, 3.0}};
    DenseMatrix state_adjoint{2, 2, {1.0, 0.0, 2.0, -1.0}};
    DenseMatrix state_control{2, 1, {1.0, 2.0}};
    DenseMatrix adjoint_state{2, 2, {1.0, -1.0, 0.0, 2.0}};
    DenseMatrix adjoint_adjoint{2, 2, {3.0, 0.0, 0.0, 4.0}};
    DenseMatrix adjoint_control{2, 1, {-1.0, 1.0}};
    DenseMatrix stationarity_state{1, 2, {1.0, 2.0}};
    DenseMatrix stationarity_adjoint{1, 2, {0.5, -1.0}};
    DenseMatrix stationarity_control{1, 1, {5.0}};
  };

  void
  require_close(const double value, const double expected, const char *message)
  {
    require(std::abs(value - expected) < 1e-12, message);
  }

  void
  require_vector_close(const DenseVector &vector,
                       const DenseVector &expected,
                       const char *       message)
  {
    require(vector.size() == expected.size(), message);
    for (std::size_t index = 0; index < vector.size(); ++index)
      require_close(vector[index], expected[index], message);
  }

  struct LifetimeState
  {
    bool owner_alive = true;
    bool callback_destroyed_after_owner = false;
  };

  struct CallbackLifetimeProbe
  {
    explicit CallbackLifetimeProbe(std::shared_ptr<LifetimeState> lifetime_state)
      : state(std::move(lifetime_state))
    {}

    ~CallbackLifetimeProbe()
    {
      if (state && !state->owner_alive)
        state->callback_destroyed_after_owner = true;
    }

    std::shared_ptr<LifetimeState> state;
  };

  struct SuppliedOTDLifetimeOwner
  {
    explicit SuppliedOTDLifetimeOwner(std::shared_ptr<LifetimeState> lifetime_state)
      : state(std::move(lifetime_state))
    {}

    ~SuppliedOTDLifetimeOwner()
    {
      state->owner_alive = false;
    }

    std::shared_ptr<LifetimeState> state;
  };

  Primal
  shifted(const Primal &point, const Primal &tangent)
  {
    Primal result = point;
    for (std::size_t block = 0; block < result.n_blocks(); ++block)
      result.add_scaled_block(block, 1.0, tangent.block(block));
    return result;
  }

  System
  make_system(
    std::shared_ptr<const CallbackLifetimeProbe> lifetime_probe = {})
  {
    const auto data = std::make_shared<const TestOperators>();
    const SuppliedOTDLayout layout(data->variable_layout,
                                   data->residual_layout);

    const auto residual = [data, lifetime_probe](const Primal &point) {
      (void)lifetime_probe;
      DenseVector state = data->state_state.vmult(point.block(0));
      state.add_scaled(1.0, data->state_adjoint.vmult(point.block(1)));
      state.add_scaled(1.0, data->state_control.vmult(point.block(2)));

      DenseVector adjoint = data->adjoint_state.vmult(point.block(0));
      adjoint.add_scaled(1.0, data->adjoint_adjoint.vmult(point.block(1)));
      adjoint.add_scaled(1.0, data->adjoint_control.vmult(point.block(2)));

      DenseVector stationarity =
        data->stationarity_state.vmult(point.block(0));
      stationarity.add_scaled(
        1.0, data->stationarity_adjoint.vmult(point.block(1)));
      stationarity.add_scaled(
        1.0, data->stationarity_control.vmult(point.block(2)));

      return Covector(data->residual_layout,
                      {std::move(state),
                       std::move(adjoint),
                       std::move(stationarity)});
    };

    const auto residual_jvp = [residual](const Primal &,
                                         const Primal &tangent) {
      return residual(tangent);
    };

    const auto residual_vjp = [data, lifetime_probe](const Primal &,
                                                     const Primal &seed) {
      (void)lifetime_probe;
      DenseVector state =
        data->state_state.transpose_vmult(seed.block(0));
      state.add_scaled(1.0,
                       data->adjoint_state.transpose_vmult(seed.block(1)));
      state.add_scaled(1.0,
                       data->stationarity_state.transpose_vmult(seed.block(2)));

      DenseVector adjoint =
        data->state_adjoint.transpose_vmult(seed.block(0));
      adjoint.add_scaled(
        1.0, data->adjoint_adjoint.transpose_vmult(seed.block(1)));
      adjoint.add_scaled(
        1.0, data->stationarity_adjoint.transpose_vmult(seed.block(2)));

      DenseVector control =
        data->state_control.transpose_vmult(seed.block(0));
      control.add_scaled(
        1.0, data->adjoint_control.transpose_vmult(seed.block(1)));
      control.add_scaled(
        1.0, data->stationarity_control.transpose_vmult(seed.block(2)));

      return Covector(data->variable_layout,
                      {std::move(state),
                       std::move(adjoint),
                       std::move(control)});
    };

    const auto solve = [data, lifetime_probe](const Primal &) {
      (void)lifetime_probe;
      const LinearSolveReport report{"test_direct",
                                     "not applicable",
                                     1,
                                     1,
                                     0.0,
                                     0.0,
                                     0.0,
                                     0.0,
                                     LinearSolveTermination::converged};
      return System::SolveResult(
        Primal(data->variable_layout,
               {DenseVector{0.0, 0.0}, DenseVector{0.0, 0.0}, DenseVector{0.0}}),
        report);
    };

    return System(layout, residual, residual_jvp, residual_vjp, solve);
  }

  void
  test_compiled_supplied_otd_owned_product_lifetime()
  {
    const auto lifetime_state = std::make_shared<LifetimeState>();
    auto product = [&]() {
      const auto probe = std::make_shared<const CallbackLifetimeProbe>(
        lifetime_state);
      auto system = std::make_shared<const System>(make_system(probe));

      nmopt::compiler::v1::CompilationManifest manifest;
      manifest.formulation_record.kind =
        nmopt::semantic::v1::FormulationKind::all_at_once;
      manifest.formulation_record.provenance =
        nmopt::semantic::v1::FormulationProvenance::supplied_otd;
      manifest.supplied_otd_record.present = true;
      manifest.supplied_otd_record.declaration =
        nmopt::semantic::v1::SuppliedOTDDeclaration{};

      return std::make_shared<const nmopt::compiler::v1::CompiledSuppliedOTDProblemT<
        DenseBackend>>(
        std::move(system),
        std::move(manifest),
        std::make_shared<const SuppliedOTDLifetimeOwner>(lifetime_state));
    }();

    const auto &system = product->system();
    const Primal point(
      system.variable_layout(),
      {DenseVector{1.0, -2.0}, DenseVector{0.5, 1.0}, DenseVector{2.0}});
    const Primal tangent(
      system.variable_layout(),
      {DenseVector{0.2, -0.4}, DenseVector{-0.3, 0.5}, DenseVector{0.7}});
    const Primal residual_seed(
      system.residual_layout(),
      {DenseVector{0.4, -0.8}, DenseVector{1.3, -0.2}, DenseVector{0.6}});

    (void)system.residual(point);
    (void)system.residual_jvp(point, tangent);
    (void)system.residual_vjp(point, residual_seed);
    const auto solve = system.solve(point);
    require(solve.report.converged(),
            "owned supplied OTD product did not retain its solve action");

    product.reset();
    require(!lifetime_state->owner_alive,
            "owned supplied OTD product did not release its lifetime owner");
    require(!lifetime_state->callback_destroyed_after_owner,
            "owned supplied OTD callbacks outlived their lifetime owner");
  }

  void
  test_supplied_otd_block_actions_and_pairing()
  {
    const System system = make_system();
    const Primal point(
      system.variable_layout(),
      {DenseVector{1.0, -2.0}, DenseVector{0.5, 1.0}, DenseVector{2.0}});
    const Covector value = system.residual(point);
    require_vector_close(value.block(0), DenseVector{2.5, -2.0},
                         "supplied OTD state block has the wrong value");
    require_vector_close(value.block(1), DenseVector{2.5, 2.0},
                         "supplied OTD adjoint block has the wrong value");
    require_vector_close(value.block(2), DenseVector{6.25},
                         "supplied OTD stationarity block has the wrong value");

    const Primal tangent(
      system.variable_layout(),
      {DenseVector{0.2, -0.4}, DenseVector{-0.3, 0.5}, DenseVector{0.7}});
    const Covector jvp = system.residual_jvp(point, tangent);
    require_vector_close(jvp.block(0), DenseVector{0.4, -0.9},
                         "supplied OTD state JVP has the wrong value");
    require_vector_close(jvp.block(1), DenseVector{-1.0, 1.9},
                         "supplied OTD adjoint JVP has the wrong value");
    require_vector_close(jvp.block(2), DenseVector{2.25},
                         "supplied OTD stationarity JVP has the wrong value");

    const Covector shifted_value = system.residual(shifted(point, tangent));
    for (std::size_t block = 0; block < value.n_blocks(); ++block)
      {
        DenseVector expected = value.block(block);
        expected.add_scaled(1.0, jvp.block(block));
        require_vector_close(shifted_value.block(block),
                             expected,
                             "supplied OTD JVP is not the finite difference action");
      }

    const Primal seed(system.residual_layout(),
                      {DenseVector{0.4, -0.8},
                       DenseVector{1.3, -0.2},
                       DenseVector{0.6}});
    require_close(pair(jvp, seed),
                  pair(system.residual_vjp(point, seed), tangent),
                  "supplied OTD VJP does not satisfy the pairing identity");
  }

  void
  test_supplied_otd_solve_and_block_descriptors()
  {
    const System system = make_system();
    const Primal initial(
      system.variable_layout(),
      {DenseVector{1.0, -2.0}, DenseVector{0.5, 1.0}, DenseVector{2.0}});
    const auto result = system.solve(initial);
    require(result.report.converged(),
            "supplied OTD solve did not preserve its convergence report");
    require(result.report.algorithm == "test_direct",
            "supplied OTD solve did not preserve its algorithm report");
    require(result.solution.layout()->compatible_with(*system.variable_layout()),
            "supplied OTD solve returned an incompatible layout");
    require(system.state_residual(initial).layout()->space(0).value ==
              "state_equation",
            "supplied OTD state descriptor selected the wrong block");
    require(system.adjoint_residual(initial).layout()->space(0).value ==
              "adjoint_equation",
            "supplied OTD adjoint descriptor selected the wrong block");
    require(system.control_stationarity(initial).layout()->space(0).value ==
              "control_stationarity",
            "supplied OTD stationarity descriptor selected the wrong block");
  }

  void
  test_supplied_otd_rejects_invalid_layouts()
  {
    const auto variable_layout = std::make_shared<const BlockLayout>(
      "variables",
      std::vector<SpaceId>{{"state"}, {"adjoint"}, {"control"}},
      std::vector<std::size_t>{1, 1, 1});
    const auto residual_layout = std::make_shared<const BlockLayout>(
      "residuals",
      std::vector<SpaceId>{{"state_equation"},
                           {"adjoint_equation"},
                           {"control_stationarity"}},
      std::vector<std::size_t>{1, 1, 1});
    SuppliedOTDBlockSelection duplicate_selection;
    duplicate_selection.adjoint_variable = duplicate_selection.state_variable;
    nmopt::test_support::require_contract_error(
      [&] {
        (void)SuppliedOTDLayout(variable_layout,
                                residual_layout,
                                duplicate_selection);
      },
      "Supplied OTD variable block selection must be distinct",
      "supplied OTD layout accepted duplicate variable blocks");

    const System system = make_system();
    const auto incompatible_layout = std::make_shared<const BlockLayout>(
      "incompatible",
      std::vector<SpaceId>{{"wrong"}, {"adjoint"}, {"control"}},
      std::vector<std::size_t>{2, 2, 1});
    const Primal incompatible = Primal::zeros(incompatible_layout);
    nmopt::test_support::require_contract_error(
      [&] { (void)system.residual(incompatible); },
      "Supplied OTD residual has an incompatible variable layout",
      "supplied OTD residual accepted an incompatible variable layout");
  }
} // namespace

int
main(const int argc, char **argv)
{
  try
    {
      const std::vector<nmopt::test_support::Scenario> scenarios{
        {"block_actions_and_pairing",
         "nmopt.supplied_otd.block_actions_and_pairing",
         {"backend-neutral", "formulation", "supplied-otd"},
         30,
         test_supplied_otd_block_actions_and_pairing},
        {"solve_and_block_descriptors",
         "nmopt.supplied_otd.solve_and_block_descriptors",
         {"backend-neutral", "formulation", "supplied-otd"},
         30,
         test_supplied_otd_solve_and_block_descriptors},
        {"invalid_layouts",
         "nmopt.supplied_otd.invalid_layouts",
         {"backend-neutral", "formulation", "supplied-otd"},
         30,
         test_supplied_otd_rejects_invalid_layouts},
        {"compiled_owned_product_lifetime",
         "nmopt.supplied_otd.compiled_owned_product_lifetime",
         {"backend-neutral", "formulation", "supplied-otd", "ownership"},
         30,
         test_compiled_supplied_otd_owned_product_lifetime}};
      const auto result = nmopt::test_support::run_requested_scenarios(
        argc, argv, scenarios, std::cout);
      if (!result.listed)
        std::cout << "supplied OTD contract scenario passed: "
                  << result.executed << '\n';
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "supplied OTD contract test failed: " << exception.what()
                << '\n';
      return 1;
    }
}
