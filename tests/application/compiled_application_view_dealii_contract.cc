#include "nmopt/application/dealii/chapter6_b1.hpp"
#include "nmopt/application/dealii/chapter6_b2.hpp"
#include "../support/scenario_dispatch.hpp"
#include "../support/scoped_temporary_directory.hpp"

#include <deal.II/base/point.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
  using namespace nmopt::application;
  using Backend = nmopt::dealii_backend::SerialBackend;

  void
  require(const bool condition, const char *message)
  {
    if (!condition)
      throw std::runtime_error(message);
  }

  void
  test_direct_volume_view()
  {
    auto scenario = chapter6::make_b1_scenario(
      chapter6::ReducedMethod::steepest_descent,
      chapter6::b1_manufactured_zero_forcing());
    scenario.compile.mesh.refinement = 1;
    scenario.problem.regularisation_sweep = {1.0e-2};

    chapter6::dealii::B1SelectedDataT<2> selected_data(
      scenario.problem.forcing, scenario.problem.desired_state);
    const auto runtime =
      chapter6::dealii::make_b1_runtime_data(scenario, selected_data);
    const auto session =
      chapter6::dealii::make_b1_compilation_session<2>(scenario);
    const auto bindings = chapter6::dealii::make_b1_data_bindings(
      scenario.problem, 1.0e-2, runtime);
    const auto specification = chapter6::make_b1_problem_spec(scenario);
    const auto policy =
      chapter6::dealii::make_b1_discretisation_policy(scenario.compile);
    const auto product =
      chapter6::dealii::make_b1_compilation_product(scenario.compile);

    nmopt::compiler::v1::DealiiCompiler compiler;
    const auto compilation = compiler.compile(specification,
                                              session,
                                              bindings,
                                              policy,
                                              std::nullopt,
                                              std::nullopt,
                                              product);
    require(compilation.succeeded(),
            "direct volume compiled application view compilation failed");
    require(compilation.problem != nullptr,
            "direct volume compiled application view has no compiled problem");

    const auto *const view = compilation.problem->compiled_application_view();
    require(view != nullptr,
            "direct volume compilation did not retain a compiled application view");
    require(!view->has_objective_components(),
            "direct volume view unexpectedly exposed objective components");
    const auto &dimensions = view->dimensions();
    require(dimensions.physical_state == 9 &&
              dimensions.independent_state == 1 &&
              dimensions.physical_control == 4 &&
              dimensions.independent_control == 4 &&
              dimensions.realized_observation == 9,
            "direct volume compiled application view retained wrong dimensions");

    const auto reduced = compilation.problem->make_reduced_dto();
    const nmopt::contract::StateControlPartitionT<Backend> partition(
      compilation.problem->executable_model(), 0, 1);
    const auto control = chapter6::dealii::detail::make_b1_uniform_control<Backend>(
      partition.control_layout(), 0.0);
    const auto evaluation = reduced.evaluate(control);
    const nmopt::test_support::ScopedTemporaryDirectory temporary_directory(
      "nmopt-native-application-view-contract");
    const auto &output_directory = temporary_directory.path();
    view->write_native_output(output_directory,
                              evaluation.state,
                              control,
                              evaluation.adjoint);

    const auto fields_path = output_directory / "fields-volume.vtu";
    require(std::filesystem::exists(fields_path) &&
              std::filesystem::file_size(fields_path) > 0,
            "compiled application view did not write field output");
    require(std::filesystem::exists(output_directory / "mesh-volume.vtu") &&
              std::filesystem::exists(output_directory / "mesh-volume.svg"),
            "compiled application view did not write mesh output");
    std::ifstream fields(fields_path);
    const std::string document((std::istreambuf_iterator<char>(fields)),
                               std::istreambuf_iterator<char>());
    require(document.find("Name=\"state\"") != std::string::npos &&
              document.find("Name=\"control\"") != std::string::npos &&
              document.find("Name=\"adjoint\"") != std::string::npos &&
              document.find("Name=\"forcing\"") != std::string::npos &&
              document.find("Name=\"target\"") != std::string::npos,
            "compiled application view omitted retained field identity");
  }

  void
  test_neumann_view()
  {
    auto scenario = chapter6::make_b2_scenario(
      chapter6::B2ObservationRegion::full, "constant");
    scenario.compile.mesh.refinement = 1;
    chapter6::dealii::B2ManufacturedDataT<2> manufactured_data(
      nmopt::application::selected_scalar_function_definition(
        scenario.problem.observation_region_catalog),
      scenario.problem.fixed_dirichlet_data,
      scenario.problem.forcing,
      chapter6::b2_target_definition(scenario.problem.target_catalog),
      scenario.problem.conservative_transport);
    const auto runtime =
      chapter6::dealii::make_b2_manufactured_runtime_data(
        scenario, manufactured_data);
    const auto session =
      chapter6::dealii::make_b2_compilation_session<2>(scenario);
    const auto bindings =
      chapter6::dealii::make_b2_data_bindings(scenario.problem, runtime);
    nmopt::compiler::v1::DealiiCompiler compiler;
    const auto compilation = compiler.compile(
      chapter6::make_b2_problem_spec(scenario),
      session,
      bindings,
      chapter6::dealii::make_b2_discretisation_policy(scenario.compile),
      std::nullopt,
      std::nullopt,
      chapter6::dealii::make_b2_compilation_product(scenario.compile));
    require(compilation.succeeded() && compilation.problem,
            "Neumann compiled application view compilation failed");

    const auto *const view = compilation.problem->compiled_application_view();
    require(view != nullptr,
            "Neumann compilation did not retain a compiled application view");
    require(view->has_objective_components(),
            "Neumann view omitted objective components");
    const nmopt::contract::StateControlPartitionT<Backend> partition(
      compilation.problem->executable_model(), 0, 1);
    const auto &dimensions = view->dimensions();
    require(dimensions.physical_state == partition.state_layout()->dimension(0) &&
              dimensions.physical_control == partition.control_layout()->dimension(0) &&
              dimensions.independent_control == partition.control_layout()->dimension(0) &&
              dimensions.realized_observation > 0,
            "Neumann compiled application view retained wrong dimensions");

    const auto reduced = compilation.problem->make_reduced_dto();
    const auto control = chapter6::dealii::detail::make_b2_uniform_control<Backend>(
      partition.control_layout(), 0.0);
    const auto evaluation = reduced.evaluate(control);
    const auto components = view->objective_components(evaluation.full_point);
    require(std::abs(components.state_tracking +
                       components.control_regularisation -
                       evaluation.objective_value) <= 1.0e-12,
            "Neumann view objective components do not reproduce the objective");

    const nmopt::test_support::ScopedTemporaryDirectory temporary_directory(
      "nmopt-neumann-application-view-contract");
    const auto &output_directory = temporary_directory.path();
    view->write_native_output(output_directory,
                              evaluation.state,
                              control,
                              evaluation.adjoint);
    require(std::filesystem::exists(output_directory / "fields-volume.vtu") &&
              std::filesystem::exists(output_directory / "control-boundary.vtu") &&
              std::filesystem::exists(output_directory / "mesh-volume.vtu"),
            "Neumann compiled application view did not write native output");
  }
} // namespace

int
main(const int argc, char **argv)
{
  try
    {
      const std::vector<nmopt::test_support::Scenario> scenarios{
        // Keep the historical scenario IDs and native-view label for inventory
        // compatibility while the physical test seam uses compiled terminology.
        {"direct_volume_view",
         "nmopt.compiler.native_application_view.direct_volume",
         {"dealii", "compiler", "application", "native-view"},
         60,
         test_direct_volume_view},
        {"neumann_view",
         "nmopt.compiler.native_application_view.neumann",
         {"dealii", "compiler", "application", "native-view"},
         60,
         test_neumann_view}};
      const auto result = nmopt::test_support::run_requested_scenarios(
        argc, argv, scenarios, std::cout);
      if (!result.listed)
        std::cout << "compiled application view deal.II scenarios passed: "
                  << result.executed << '\n';
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "compiled application view deal.II test failed: "
                << exception.what() << '\n';
      return 1;
    }
}
