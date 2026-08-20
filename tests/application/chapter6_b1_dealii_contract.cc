#include "nmopt/application/dealii/chapter6_b1.hpp"
#include "../support/scenario_dispatch.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>

namespace
{
  using namespace nmopt::application;

  void
  require(const bool condition, const char *message)
  {
    if (!condition)
      throw std::runtime_error(message);
  }

  void
  run_b1_case(
    const chapter6::ReducedMethod method,
    const chapter5::DistributedControlDiscretisation discretisation =
      chapter5::DistributedControlDiscretisation::cellwise_constant,
    const chapter6::B1ForcingSelection forcing =
      chapter6::B1ForcingSelection::manufactured_zero)
  {
    auto scenario = chapter6::make_b1_scenario(method, forcing);
    scenario.problem.recipe.discretisation = discretisation;
    scenario.compile.mesh.refinement = 1;
    scenario.problem.regularisation_sweep = {1.0e-2};
    scenario.solver.parameters.maximum_iterations = 20;
    scenario.solver.parameters.gradient_tolerance = 1.0e-3;
    scenario.compile.state_solve = {123, 2.0e-12, 3.0e-14};
    scenario.compile.adjoint_solve = {124, 4.0e-12, 5.0e-14};
    scenario.compile.control_metric_solve = {321, 6.0e-12, 7.0e-14};

    chapter6::dealii::B1SelectedDataT<2> selected_data(forcing);
    const auto runtime =
      chapter6::dealii::make_b1_runtime_data(scenario, selected_data);
    const double expected_forcing =
      forcing == chapter6::B1ForcingSelection::figure_inferred_constant_one
        ? 1.0
        : 0.0;
    require(std::abs(runtime.forcing.value(::dealii::Point<2>(0.5, 0.5)) -
                     expected_forcing) < 1.0e-15,
            "B1 runtime data realized the wrong forcing value");
    const auto session =
      chapter6::dealii::make_b1_compilation_session<2>(scenario);
    const auto native_output_directory =
      std::filesystem::temp_directory_path() /
      ("nmopt-b1-native-contract-" +
       std::to_string(static_cast<int>(method)) + "-" +
       chapter5::distributed_control_discretisation_name(discretisation) + "-" +
       chapter6::b1_forcing_selection_name(forcing));
    std::filesystem::remove_all(native_output_directory);
    chapter6::dealii::B1ReducedExecutionAdapterT<2> execute{
      1.0e-2,
      runtime,
      session,
      {"test.chapter6.b1",
       "debug-dealii",
       "test-compiler",
       "test-version",
       "test-standard-library",
       "test-os",
       "test-architecture",
       "test-hardware"},
      native_output_directory};

    using HeadlessRunner =
      nmopt::application::benchmark::HeadlessBenchmarkRunnerT<
        decltype(scenario)>;
    const auto result = HeadlessRunner(scenario).run(
      [](const auto &parameters) {
        return chapter6::make_b1_problem_spec(parameters);
      },
      execute);

    require(result.artifact.identity().scenario_id ==
              "chapter-6.b1.distributed-laplace",
            "B1 dealii adapter lost scenario identity");
    require(result.artifact.envelope().report().objective_history.size() > 1,
            "B1 dealii adapter did not execute a solver iteration");
    require(result.artifact.envelope().report().state_solve_count > 0 &&
              result.artifact.envelope().report().adjoint_solve_count > 0,
            "B1 dealii adapter did not retain solve counts");
    const auto &manifest = result.artifact.envelope().compilation_manifest();
    require(manifest.state_solve_record.maximum_iterations == 123 &&
              manifest.state_solve_record.relative_tolerance == 2.0e-12 &&
              manifest.state_solve_record.absolute_tolerance == 3.0e-14 &&
              manifest.adjoint_solve_record.maximum_iterations == 124 &&
              manifest.adjoint_solve_record.relative_tolerance == 4.0e-12 &&
              manifest.adjoint_solve_record.absolute_tolerance == 5.0e-14,
            "B1 dealii adapter did not map the state/adjoint solve policies");
    require(manifest.metric_record.solve_policy.maximum_iterations == 321 &&
              manifest.metric_record.solve_policy.relative_tolerance ==
                6.0e-12 &&
              manifest.metric_record.solve_policy.absolute_tolerance ==
                7.0e-14,
            "B1 dealii adapter did not map the control-metric solve policy");
    if (discretisation ==
        chapter5::DistributedControlDiscretisation::
          homogeneous_dirichlet_continuous)
      require(manifest.control_space.find(
                "homogeneous-Dirichlet scalar FE_Q(1)") != std::string::npos &&
                manifest.metric_solve_policy.find("l2_continuous") !=
                  std::string::npos,
              "B1 continuous candidate selected the wrong compiled control");
    else
      require(manifest.control_space.find("FE_DGQ(0)") != std::string::npos,
              "B1 cellwise candidate selected the wrong compiled control");
    require(result.document.find("b1.regularisation_weight=0.01\n") !=
              std::string::npos,
            "B1 dealii adapter omitted regularisation evidence");
    require(result.document.find(
              std::string("b1.control_discretisation=") +
              chapter5::distributed_control_discretisation_name(
                discretisation) +
              "\n") != std::string::npos,
            "B1 dealii adapter omitted control-discretisation evidence");
    require(result.document.find(
              std::string("b1.forcing_selection=") +
              chapter6::b1_forcing_selection_name(forcing) + "\n") !=
              std::string::npos,
            "B1 dealii adapter omitted forcing-selection evidence");
    require(result.document.find("benchmark.state_dimension=9\n") !=
              std::string::npos &&
              result.document.find(
                "benchmark.state_physical_dimension=9\n") !=
                std::string::npos &&
              result.document.find(
                "benchmark.state_independent_dimension=1\n") !=
                std::string::npos &&
              result.document.find("benchmark.adjoint_dimension=9\n") !=
                std::string::npos &&
              result.document.find(
                "benchmark.adjoint_physical_dimension=9\n") !=
                std::string::npos &&
              result.document.find(
                "benchmark.adjoint_independent_dimension=1\n") !=
                std::string::npos,
            "B1 dealii adapter omitted state/adjoint dimension conventions");
    const auto expected_control_dimension =
      discretisation == chapter5::DistributedControlDiscretisation::
                          homogeneous_dirichlet_continuous
        ? "1"
        : "4";
    const auto expected_physical_control_dimension =
      discretisation == chapter5::DistributedControlDiscretisation::
                          homogeneous_dirichlet_continuous
        ? "9"
        : "4";
    require(result.document.find(
              std::string("benchmark.control_dimension=") +
              expected_control_dimension + "\n") != std::string::npos &&
              result.document.find(
                std::string("benchmark.control_physical_dimension=") +
                expected_physical_control_dimension + "\n") !=
                std::string::npos &&
              result.document.find(
                std::string("benchmark.control_independent_dimension=") +
                expected_control_dimension + "\n") != std::string::npos,
            "B1 dealii adapter omitted control dimension conventions");
    require(result.document.find("solver.method=") != std::string::npos,
            "B1 dealii adapter omitted solver-method evidence");
    require(result.document.find(
              "b1.hessian_evidence=centered_finite_difference\n") !=
              std::string::npos,
            "B1 dealii adapter omitted Hessian evidence method");
    require(result.document.find(
              "b1.hessian_finite_difference_passed=true\n") !=
              std::string::npos,
            "B1 dealii adapter omitted Hessian finite-difference evidence");
    require(
      std::filesystem::exists(native_output_directory / "fields-volume.vtu"),
            "B1 dealii adapter did not write field output");
    require(
      std::filesystem::exists(native_output_directory / "mesh-volume.vtu"),
      "B1 dealii adapter did not write the volume mesh");
    require(
      std::filesystem::exists(native_output_directory / "mesh-volume.svg"),
      "B1 dealii adapter did not write the volume mesh SVG");
    std::ifstream fields(native_output_directory / "fields-volume.vtu");
    std::ifstream mesh(native_output_directory / "mesh-volume.vtu");
    std::ifstream mesh_svg(native_output_directory / "mesh-volume.svg");
    const std::string field_document((std::istreambuf_iterator<char>(fields)),
                                     std::istreambuf_iterator<char>());
    const std::string mesh_document((std::istreambuf_iterator<char>(mesh)),
                                    std::istreambuf_iterator<char>());
    const std::string mesh_svg_document(
      (std::istreambuf_iterator<char>(mesh_svg)),
      std::istreambuf_iterator<char>());
    require(mesh_document.find("<UnstructuredGrid>") != std::string::npos,
            "B1 volume mesh output is not a VTU unstructured grid");
    require(mesh_svg_document.find("<svg") != std::string::npos,
            "B1 volume mesh SVG output is not SVG");
    require(field_document.find("Name=\"state\"") != std::string::npos &&
              field_document.find("Name=\"control\"") != std::string::npos &&
              field_document.find("Name=\"adjoint\"") != std::string::npos &&
              field_document.find("Name=\"negative_adjoint\"") !=
                std::string::npos &&
              field_document.find("Name=\"target\"") != std::string::npos &&
              field_document.find("Name=\"forcing\"") != std::string::npos,
            "B1 field output omitted a retained field");
    std::filesystem::remove_all(native_output_directory);
  }

  void
  test_b1_simplex_mesh_generation()
  {
    auto scenario = chapter6::make_b1_scenario(
      chapter6::ReducedMethod::steepest_descent);
    scenario.problem.recipe.discretisation = chapter5::
      DistributedControlDiscretisation::homogeneous_dirichlet_continuous;
    scenario.problem.recipe.with_cellwise_box = false;
    scenario.compile.mesh.generation =
      chapter6::MeshGeneration::structured_simplex;
    scenario.compile.mesh.refinement = 0;
    scenario.compile.mesh.subdivisions = 3;
    scenario.compile.mesh.mesh_provenance =
      "test.chapter6.b1.structured-simplex-n3";

    const auto structured_session =
      chapter6::dealii::make_b1_compilation_session<2>(scenario);
    const auto &structured_mesh = structured_session->triangulation();
    require(structured_mesh.all_reference_cells_are_simplex() &&
              structured_mesh.n_vertices() == 16 &&
              structured_mesh.n_active_cells() == 18,
            "B1 structured simplex generator produced the wrong topology");

    scenario.compile.mesh.generation =
      chapter6::MeshGeneration::centroid_split_simplex;
    scenario.compile.mesh.subdivisions = 100;
    scenario.compile.mesh.centroid_splits = 7160;
    scenario.compile.mesh.selection_seed = 0;
    scenario.compile.mesh.mesh_provenance =
      "test.chapter6.b1.count-matched-simplex-n100-s7160-seed0";

    const auto count_matched_session =
      chapter6::dealii::make_b1_compilation_session<2>(scenario);
    const auto &count_matched_mesh = count_matched_session->triangulation();
    require(count_matched_mesh.all_reference_cells_are_simplex() &&
              count_matched_mesh.n_vertices() == 17361 &&
              count_matched_mesh.n_active_cells() == 34320,
            "B1 count-matched simplex generator missed the source cell and vertex counts");

    std::vector<bool> boundary_vertices(count_matched_mesh.n_vertices(), false);
    for (const auto &cell : count_matched_mesh.active_cell_iterators())
      for (const auto face : cell->face_indices())
        if (cell->face(face)->at_boundary())
          for (const auto vertex : cell->face(face)->vertex_indices())
            boundary_vertices[cell->face(face)->vertex_index(vertex)] = true;
    const auto boundary_vertex_count =
      std::count(boundary_vertices.begin(), boundary_vertices.end(), true);
    require(boundary_vertex_count == 400 &&
              count_matched_mesh.n_vertices() - boundary_vertex_count == 16961,
            "B1 count-matched simplex generator missed the inferred boundary and independent P1 counts");
  }
} // namespace

int
main(const int argc, char **argv)
{
  try
    {
      const std::vector<nmopt::test_support::Scenario> scenarios{
        {"b1_manufactured_steepest_descent",
         "nmopt.application.dealii.b1_manufactured_steepest_descent",
         {"dealii", "application", "benchmark", "b1", "contract"},
         120,
         []() {
           run_b1_case(
             nmopt::application::chapter6::ReducedMethod::steepest_descent);
         }},
        {"b1_manufactured_limited_memory_bfgs",
         "nmopt.application.dealii.b1_manufactured_limited_memory_bfgs",
         {"dealii", "application", "benchmark", "b1", "contract"},
         120,
         []() {
           run_b1_case(
             nmopt::application::chapter6::ReducedMethod::limited_memory_bfgs);
         }},
        {"b1_continuous_control_steepest_descent",
         "nmopt.application.dealii.b1_continuous_control_steepest_descent",
         {"dealii", "application", "benchmark", "b1", "contract"},
         120,
         []() {
           run_b1_case(
             nmopt::application::chapter6::ReducedMethod::steepest_descent,
             nmopt::application::chapter5::DistributedControlDiscretisation::
               homogeneous_dirichlet_continuous);
         }},
        {"b1_inferred_constant_one_steepest_descent",
         "nmopt.application.dealii.b1_inferred_constant_one_steepest_descent",
         {"dealii", "application", "benchmark", "b1", "contract"},
         120,
         []() {
           run_b1_case(
             nmopt::application::chapter6::ReducedMethod::steepest_descent,
             nmopt::application::chapter5::DistributedControlDiscretisation::
               homogeneous_dirichlet_continuous,
             nmopt::application::chapter6::B1ForcingSelection::
               figure_inferred_constant_one);
         }},
        {"b1_simplex_mesh_generation",
         "nmopt.application.dealii.b1_simplex_mesh_generation",
         {"dealii", "application", "benchmark", "b1", "contract"},
         60,
         test_b1_simplex_mesh_generation}};
      const auto result = nmopt::test_support::run_requested_scenarios(
        argc, argv, scenarios, std::cout);
      if (!result.listed)
        std::cout << "chapter 6 B1 deal.II contract test passed: "
                  << result.executed << '\n';
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "chapter 6 B1 deal.II contract test failed: "
                << exception.what() << '\n';
      return 1;
    }
}
