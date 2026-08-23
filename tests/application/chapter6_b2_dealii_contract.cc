#include "nmopt/application/dealii/chapter6_b2.hpp"
#include "../support/scenario_dispatch.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
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

  double
  artifact_number(const std::string &document, const std::string &key)
  {
    const auto prefix = key + "=";
    const auto begin = document.find(prefix);
    if (begin == std::string::npos)
      throw std::runtime_error("B2 artifact omitted numeric field '" + key +
                               "'");
    const auto value_begin = begin + prefix.size();
    const auto value_end = document.find('\n', value_begin);
    return std::stod(document.substr(value_begin, value_end - value_begin));
  }

  void
  run_b2_manufactured_case(const chapter6::GraetzCase graetz_case)
  {
    auto scenario = chapter6::make_b2_scenario(graetz_case);
    scenario.compile.mesh.refinement = 1;
    scenario.solver.parameters.maximum_iterations = 5;
    scenario.solver.parameters.gradient_tolerance = 1.0e-4;
    scenario.compile.state_solve = {125, 2.0e-12, 3.0e-14};
    scenario.compile.adjoint_solve = {126, 4.0e-12, 5.0e-14};
    scenario.compile.control_metric_solve = {322, 6.0e-12, 7.0e-14};

    chapter6::dealii::B2ManufacturedDataT<2> manufactured_data{
      graetz_case, scenario.problem.fixed_temperature};
    const auto runtime =
      chapter6::dealii::make_b2_manufactured_runtime_data(
        scenario, manufactured_data);
    const auto session =
      chapter6::dealii::make_b2_compilation_session<2>(scenario);
    require(session->triangulation().all_reference_cells_are_hyper_cube(),
            "B2 default mesh is no longer the native quadrilateral rectangle");
    const auto native_output_directory =
      std::filesystem::temp_directory_path() /
      ("nmopt-b2-native-contract-" +
       std::to_string(static_cast<int>(graetz_case)));
    std::filesystem::remove_all(native_output_directory);

    const auto specification = chapter6::make_b2_problem_spec(scenario);
    const auto fixed_objective = [&] {
      const auto bindings =
        chapter6::dealii::make_b2_data_bindings(scenario.problem, runtime);
      nmopt::compiler::v1::DealiiCompiler compiler;
      const auto compilation = compiler.compile(specification,
                                                session,
                                                bindings,
                                                {},
                                                std::nullopt,
                                                std::nullopt,
                                                nmopt::compiler::v1::CompilationProduct::reduced_dto);
      require(compilation.succeeded() && compilation.problem,
              "B2 fixed-temperature compilation did not succeed");
      const auto reduced = compilation.problem->make_reduced_dto();
      const nmopt::contract::StateControlPartitionT<chapter6::dealii::Backend> partition(
        compilation.problem->executable_model(), 0, 1);
      const auto control =
        nmopt::contract::PrimalBlockT<chapter6::dealii::Backend>::zeros(
          partition.control_layout());
      return reduced.evaluate(control).objective_value;
    }();

    const double zero_fixed_objective = [&] {
      ::dealii::Functions::ZeroFunction<2> zero_fixed_temperature;
      const chapter6::dealii::B2RuntimeDataT<2> zero_runtime{
        runtime.forcing,
        runtime.desired_state,
        zero_fixed_temperature,
        runtime.conservative_transport};
      const auto zero_session =
        chapter6::dealii::make_b2_compilation_session<2>(scenario);
      const auto bindings = chapter6::dealii::make_b2_data_bindings(
        scenario.problem, zero_runtime);
      nmopt::compiler::v1::DealiiCompiler compiler;
      const auto compilation = compiler.compile(specification,
                                                zero_session,
                                                bindings,
                                                {},
                                                std::nullopt,
                                                std::nullopt,
                                                nmopt::compiler::v1::CompilationProduct::reduced_dto);
      require(compilation.succeeded() && compilation.problem,
              "B2 zero-temperature comparison compilation did not succeed");
      const auto reduced = compilation.problem->make_reduced_dto();
      const nmopt::contract::StateControlPartitionT<chapter6::dealii::Backend> partition(
        compilation.problem->executable_model(), 0, 1);
      const auto control =
        nmopt::contract::PrimalBlockT<chapter6::dealii::Backend>::zeros(
          partition.control_layout());
      return reduced.evaluate(control).objective_value;
    }();
    require(std::abs(fixed_objective - zero_fixed_objective) > 1.0e-8,
            "B2 fixed-temperature binding did not affect the objective");

    chapter6::dealii::B2ReducedExecutionAdapterT<2> execute{
      runtime,
      session,
      {"test.chapter6.b2",
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
        return chapter6::make_b2_problem_spec(parameters);
      },
      execute);

    require(result.artifact.identity().scenario_id == scenario.metadata.id,
            "B2 deal.II adapter lost scenario identity");
    require(result.artifact.envelope().report().objective_history.size() > 1,
            "B2 deal.II adapter did not execute a solver iteration");
    require(result.artifact.envelope().report().state_solve_count > 0 &&
              result.artifact.envelope().report().adjoint_solve_count > 0,
            "B2 deal.II adapter did not retain solve counts");
    const auto &manifest = result.artifact.envelope().compilation_manifest();
    require(manifest.state_solve_record.algorithm ==
                nmopt::compiler::v1::LinearSolveAlgorithm::
                  serial_sparse_direct_umfpack &&
              manifest.state_solve_record.maximum_iterations == 1 &&
              manifest.state_solve_record.relative_tolerance == 0.0 &&
              manifest.adjoint_solve_record.algorithm ==
                nmopt::compiler::v1::LinearSolveAlgorithm::
                  serial_sparse_direct_umfpack &&
              manifest.adjoint_solve_record.maximum_iterations == 1 &&
              manifest.adjoint_solve_record.relative_tolerance == 0.0,
            "B2 manifest did not distinguish its direct state/adjoint solves");
    require(manifest.metric_record.solve_policy.maximum_iterations == 322 &&
              manifest.metric_record.solve_policy.relative_tolerance ==
                6.0e-12 &&
              manifest.metric_record.solve_policy.absolute_tolerance ==
                7.0e-14,
            "B2 dealii adapter did not map the control-metric solve policy");
    require(manifest.boundary_realisation.has_value() &&
              manifest.boundary_realisation->transport_boundary_form ==
                nmopt::semantic::v1::TransportBoundaryForm::ordinary_normal_minus_transport,
            "B2 manifest did not retain the ordinary-normal boundary form");
    require(std::any_of(manifest.declared_assumptions.begin(),
                        manifest.declared_assumptions.end(),
                        [](const std::string &assumption) {
                          return assumption.find(
                                   "boundary_form=ordinary partial_n(y) - (b dot n)y") !=
                                 std::string::npos;
                        }),
            "B2 manifest did not describe the ordinary-normal boundary form");
    const auto has_boundary_region = [&manifest](
                                       const std::string &semantic_id,
                                       const unsigned int boundary_id) {
      return std::any_of(
        manifest.resolved_decision.regions.begin(),
        manifest.resolved_decision.regions.end(),
        [&semantic_id, boundary_id](const auto &region) {
          return region.semantic_id == semantic_id &&
                 region.boundary_ids == std::vector<unsigned int>{boundary_id};
        });
    };
    require(has_boundary_region(
              chapter6::b2_fixed_boundary_region_id,
              chapter6::b2_fixed_boundary_id) &&
              has_boundary_region(
                chapter6::b2_control_boundary_region_id,
                chapter6::b2_control_boundary_id) &&
              has_boundary_region(
                chapter6::b2_outflow_boundary_region_id,
                chapter6::b2_outflow_boundary_id),
            "B2 compiled manifest did not retain the fixed/control/outflow partition");
    require(manifest.compiler_id ==
              "nmopt.compiler.v1.dealii.neumann_convection_subdomain",
            "B2 deal.II adapter selected the wrong compiler target");
    require(manifest.lifting_realisation.find("ell_0") != std::string::npos,
            "B2 deal.II adapter did not retain fixed lifting evidence");
    require(std::any_of(manifest.bindings.begin(),
                        manifest.bindings.end(),
                        [](const auto &binding) {
                          return binding.semantic_id == "fixed_dirichlet_data";
                        }),
            "B2 deal.II adapter did not retain the fixed-data binding");
    require(std::any_of(manifest.bindings.begin(),
                        manifest.bindings.end(),
                        [](const auto &binding) {
                          return binding.semantic_id == "fixed_dirichlet_data" &&
                                 binding.provenance ==
                                   "chapter-6.e6.5.2.fixed-temperature";
                        }),
            "B2 deal.II adapter did not retain fixed-data provenance");
    require(result.document.find("b2.graetz_case=") != std::string::npos &&
              result.document.find(
                "b2.control_discretisation=facewise-constant\n") !=
                std::string::npos &&
              result.document.find(
                "b2.volume_observation_quadrature_order=3\n") !=
                std::string::npos &&
              result.document.find(
                "b2.volume_observation_target_realisation=analytic-quadrature\n") !=
                std::string::npos,
            "B2 deal.II adapter omitted case, control, or observation evidence");
    require(result.document.find("b2.fixed_temperature=1") !=
              std::string::npos,
            "B2 deal.II adapter omitted fixed-temperature evidence");
    const auto expected_observation_region =
      graetz_case == chapter6::GraetzCase::observation_wings_constant_target ||
          graetz_case == chapter6::GraetzCase::observation_wings_parabolic_target
        ? "wings"
        : "full";
    const auto expected_target_profile =
      graetz_case == chapter6::GraetzCase::observation_wings_constant_target ||
          graetz_case == chapter6::GraetzCase::observation_full_constant_target
        ? "constant"
        : "parabolic";
    require(result.document.find(
              std::string("b2.observation_region=") +
              expected_observation_region + "\n") != std::string::npos,
            "B2 deal.II adapter omitted observation-region evidence");
    require(result.document.find(std::string("b2.target_profile=") +
                                 expected_target_profile + "\n") !=
              std::string::npos,
            "B2 deal.II adapter omitted target-profile evidence");
    require(result.document.find("b2.residual_jvp_error=") !=
              std::string::npos &&
              result.document.find("b2.residual_vjp_error=") !=
                std::string::npos &&
              result.document.find(
                "b2.reduced_gradient_finite_difference_error=") !=
                std::string::npos &&
              result.document.find("b2.reduced_taylor_order=") !=
                std::string::npos,
            "B2 deal.II adapter omitted derivative evidence");
    require(result.document.find("b2.derivative_evidence_passed=true\n") !=
              std::string::npos &&
              result.document.find("b2.initial_objective=") !=
                std::string::npos &&
              result.document.find("b2.final_objective=") !=
                std::string::npos &&
              result.document.find("b2.relative_gradient_reduction=") !=
                std::string::npos,
            "B2 deal.II adapter omitted reduction evidence");
    require(result.document.find("b2.state_l2_norm=") !=
              std::string::npos &&
              result.document.find("b2.control_l2_norm=") !=
                std::string::npos,
            "B2 deal.II adapter omitted state/control evidence");
    require(result.document.find("benchmark.mesh_vertices=9\n") !=
                std::string::npos &&
              result.document.find("benchmark.mesh_active_cells=4\n") !=
                std::string::npos &&
              result.document.find("benchmark.boundary_face_count=8\n") !=
                std::string::npos &&
              result.document.find("benchmark.fixed_boundary_face_count=2\n") !=
                std::string::npos &&
              result.document.find("benchmark.control_boundary_face_count=4\n") !=
                std::string::npos &&
              result.document.find("benchmark.outflow_boundary_face_count=2\n") !=
                std::string::npos,
            "B2 deal.II adapter omitted exact mesh and boundary counts");
    require(result.document.find("benchmark.state_dimension=9\n") !=
                std::string::npos &&
              result.document.find("benchmark.state_physical_dimension=9\n") !=
                std::string::npos &&
              result.document.find("benchmark.state_independent_dimension=6\n") !=
                std::string::npos &&
              result.document.find("benchmark.control_dimension=4\n") !=
                std::string::npos &&
              result.document.find("benchmark.control_physical_dimension=4\n") !=
                std::string::npos &&
              result.document.find("benchmark.control_independent_dimension=4\n") !=
                std::string::npos &&
              result.document.find("benchmark.adjoint_dimension=9\n") !=
                std::string::npos &&
              result.document.find("benchmark.adjoint_independent_dimension=6\n") !=
                std::string::npos,
            "B2 deal.II adapter omitted physical and independent dimensions");
    require(std::abs(artifact_number(result.document,
                                     "benchmark.observation_measure") -
                       2.0) < 1.0e-12,
            "B2 deal.II adapter recorded the wrong observation measure");
    const double initial_tracking = artifact_number(
      result.document, "b2.initial_tracking_objective");
    const double initial_regularisation = artifact_number(
      result.document, "b2.initial_control_regularisation_objective");
    const double final_tracking = artifact_number(
      result.document, "b2.final_tracking_objective");
    const double final_regularisation = artifact_number(
      result.document, "b2.final_control_regularisation_objective");
    require(std::abs(initial_tracking + initial_regularisation -
                     artifact_number(result.document, "b2.initial_objective")) <
                1.0e-12 &&
              std::abs(final_tracking + final_regularisation -
                       artifact_number(result.document, "b2.final_objective")) <
                1.0e-12,
            "B2 artifact objective components do not reproduce the totals");
    require(result.document.find("solver.metric_gradient_norm_history=") !=
                std::string::npos &&
              result.document.find(
                "solver.coefficient_derivative_norm_history=") !=
                std::string::npos &&
              result.document.find(
                "b2.initial_coefficient_derivative_norm=") !=
                std::string::npos &&
              result.document.find(
                "b2.final_coefficient_derivative_norm=") !=
                std::string::npos,
            "B2 deal.II adapter omitted the two derivative norm conventions");
    for (const auto *key : {"b2.uncontrolled_state_min",
                            "b2.uncontrolled_state_max",
                            "b2.optimized_state_min",
                            "b2.optimized_state_max",
                            "b2.control_min",
                            "b2.control_max",
                            "b2.adjoint_min",
                            "b2.adjoint_max"})
      require(std::isfinite(artifact_number(result.document, key)),
              "B2 deal.II adapter emitted a non-finite field extremum");
    require(result.document.find("solver.method=bfgs\n") !=
            std::string::npos,
            "B2 deal.II adapter omitted solver-method evidence");
    require(std::filesystem::exists(
              native_output_directory / "fields-volume.vtu") &&
              std::filesystem::exists(
                native_output_directory / "control-boundary.vtu"),
            "B2 deal.II adapter did not write field output");
    require(
      std::filesystem::exists(native_output_directory / "mesh-volume.vtu"),
      "B2 deal.II adapter did not write the volume mesh");
    require(
      std::filesystem::exists(native_output_directory / "mesh-volume.svg"),
      "B2 dealii adapter did not write the volume mesh SVG");
    std::ifstream fields(native_output_directory / "fields-volume.vtu");
    std::ifstream mesh(native_output_directory / "mesh-volume.vtu");
    std::ifstream mesh_svg(native_output_directory / "mesh-volume.svg");
    std::ifstream control(native_output_directory / "control-boundary.vtu");
    const std::string field_document((std::istreambuf_iterator<char>(fields)),
                                     std::istreambuf_iterator<char>());
    const std::string mesh_document((std::istreambuf_iterator<char>(mesh)),
                                    std::istreambuf_iterator<char>());
    const std::string mesh_svg_document(
      (std::istreambuf_iterator<char>(mesh_svg)),
      std::istreambuf_iterator<char>());
    const std::string control_document(
      (std::istreambuf_iterator<char>(control)),
      std::istreambuf_iterator<char>());
    require(mesh_document.find("<UnstructuredGrid>") != std::string::npos,
            "B2 volume mesh output is not a VTU unstructured grid");
    require(mesh_svg_document.find("<svg") != std::string::npos,
            "B2 volume mesh SVG output is not SVG");
    require(field_document.find("Name=\"state\"") != std::string::npos &&
              field_document.find("Name=\"state_uncontrolled\"") !=
                std::string::npos &&
              field_document.find("Name=\"target\"") != std::string::npos &&
              field_document.find("Name=\"forcing\"") != std::string::npos &&
              field_document.find("Name=\"observation_region\"") !=
                std::string::npos &&
              field_document.find("Name=\"adjoint\"") != std::string::npos &&
              control_document.find("Name=\"control\"") != std::string::npos,
            "B2 field output omitted a retained state or input field");
    require(
      control_document.find(
        "<Piece NumberOfPoints=\"8\" NumberOfCells=\"4\">") !=
          std::string::npos &&
        control_document.find("<CellData Scalars=\"control\">") !=
          std::string::npos &&
        control_document.find(
          "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n"
          "3 3 3 3 ") != std::string::npos,
      "B2 facewise control output changed its cell topology or data association");
    std::filesystem::remove_all(native_output_directory);
  }

  void
  run_b2_structured_simplex_mesh()
  {
    auto scenario = chapter6::make_b2_scenario(
      chapter6::GraetzCase::observation_wings_constant_target);
    scenario.compile.mesh.generation =
      chapter6::MeshGeneration::structured_simplex;
    scenario.compile.mesh.refinement = 0;
    scenario.compile.mesh.axis_subdivisions = {4, 10};
    scenario.compile.mesh.mesh_provenance =
      "test.chapter6.b2.structured-simplex-4x10";
    scenario.solver.parameters.maximum_iterations = 2;
    scenario.solver.parameters.gradient_tolerance = 1.0e-4;
    scenario.experiment.retain_fields = false;

    const auto session =
      chapter6::dealii::make_b2_compilation_session<2>(scenario);
    const auto &mesh = session->triangulation();
    require(mesh.all_reference_cells_are_simplex() &&
              mesh.n_vertices() == 55 && mesh.n_active_cells() == 80,
            "B2 structured simplex generator produced the wrong topology");

    std::size_t observed_cell_count = 0;
    std::size_t boundary_face_count = 0;
    std::size_t fixed_boundary_face_count = 0;
    std::size_t control_boundary_face_count = 0;
    std::size_t outflow_boundary_face_count = 0;
    double      observation_measure = 0.0;
    for (const auto &cell : mesh.active_cell_iterators())
      {
        require(cell->material_id() == 0 ||
                  cell->material_id() ==
                    scenario.problem.recipe.observed_material_id,
                "B2 structured simplex mesh has an unknown material id");
        if (cell->material_id() ==
            scenario.problem.recipe.observed_material_id)
          {
            ++observed_cell_count;
            observation_measure += cell->measure();
          }
        for (unsigned int face = 0; face < cell->n_faces(); ++face)
          if (cell->face(face)->at_boundary())
            {
              ++boundary_face_count;
              switch (cell->face(face)->boundary_id())
                {
                  case chapter6::b2_fixed_boundary_id:
                    ++fixed_boundary_face_count;
                    break;
                  case chapter6::b2_control_boundary_id:
                    ++control_boundary_face_count;
                    break;
                  case chapter6::b2_outflow_boundary_id:
                    ++outflow_boundary_face_count;
                    break;
                  default:
                    throw std::runtime_error(
                      "B2 structured simplex mesh has an unclassified boundary face");
                }
            }
      }
    require(observed_cell_count == 36 && boundary_face_count == 28 &&
              fixed_boundary_face_count == 12 &&
              control_boundary_face_count == 6 &&
              outflow_boundary_face_count == 10,
            "B2 structured simplex classification produced the wrong counts");
    require(std::abs(observation_measure - 1.8) < 1.0e-12,
            "B2 structured simplex mesh did not align the wings observation region");

    chapter6::dealii::B2ManufacturedDataT<2> manufactured_data{
      scenario.problem.graetz_case, scenario.problem.fixed_temperature};
    const auto runtime =
      chapter6::dealii::make_b2_manufactured_runtime_data(
        scenario, manufactured_data);
    chapter6::dealii::B2ReducedExecutionAdapterT<2> execute{
      runtime,
      session,
      {"test.chapter6.b2.structured-simplex",
       "debug-dealii",
       "test-compiler",
       "test-version",
       "test-standard-library",
       "test-os",
       "test-architecture",
       "test-hardware"}};
    using HeadlessRunner =
      nmopt::application::benchmark::HeadlessBenchmarkRunnerT<
        decltype(scenario)>;
    const auto result = HeadlessRunner(scenario).run(
      [](const auto &parameters) {
        return chapter6::make_b2_problem_spec(parameters);
      },
      execute);
    require(result.artifact.envelope().report().state_solve_count > 0 &&
              result.artifact.envelope().report().adjoint_solve_count > 0,
            "B2 structured simplex execution did not solve state and adjoint systems");
    require(result.document.find("b2.derivative_evidence_passed=true\n") !=
                std::string::npos &&
              result.document.find("benchmark.mesh_vertices=55\n") !=
                std::string::npos &&
              result.document.find("benchmark.mesh_active_cells=80\n") !=
                std::string::npos &&
              std::abs(artifact_number(result.document,
                                       "benchmark.observation_measure") -
                       1.8) < 1.0e-12,
            "B2 structured simplex execution lost numerical or structural evidence");
  }

  void
  run_b2_centroid_split_simplex_mesh()
  {
    constexpr unsigned int base_vertex_count = 55;
    constexpr unsigned int split_count = 7;

    auto scenario = chapter6::make_b2_scenario(
      chapter6::GraetzCase::observation_wings_constant_target);
    scenario.compile.mesh.generation =
      chapter6::MeshGeneration::centroid_split_simplex;
    scenario.compile.mesh.refinement = 0;
    scenario.compile.mesh.axis_subdivisions = {4, 10};
    scenario.compile.mesh.centroid_splits = split_count;
    scenario.compile.mesh.selection_seed = 19;
    scenario.compile.mesh.mesh_provenance =
      "test.chapter6.b2.centroid-split-4x10-s7-seed19";
    scenario.solver.parameters.maximum_iterations = 2;
    scenario.solver.parameters.gradient_tolerance = 1.0e-4;
    scenario.experiment.retain_fields = false;

    const auto session =
      chapter6::dealii::make_b2_compilation_session<2>(scenario);
    const auto repeated_session =
      chapter6::dealii::make_b2_compilation_session<2>(scenario);
    auto alternate_scenario = scenario;
    alternate_scenario.compile.mesh.selection_seed = 20;
    alternate_scenario.compile.mesh.mesh_provenance =
      "test.chapter6.b2.centroid-split-4x10-s7-seed20";
    const auto alternate_session =
      chapter6::dealii::make_b2_compilation_session<2>(alternate_scenario);

    const auto added_vertices = [base_vertex_count](const auto &mesh) {
      std::vector<std::array<double, 2>> result;
      for (unsigned int index = base_vertex_count;
           index < mesh.n_vertices();
           ++index)
        result.push_back(
          {{mesh.get_vertices()[index][0], mesh.get_vertices()[index][1]}});
      return result;
    };
    const auto same_vertices = [](const auto &first, const auto &second) {
      if (first.size() != second.size())
        return false;
      for (std::size_t index = 0; index < first.size(); ++index)
        for (unsigned int coordinate = 0; coordinate < 2; ++coordinate)
          if (std::abs(first[index][coordinate] -
                       second[index][coordinate]) > 1.0e-15)
            return false;
      return true;
    };

    const auto &mesh = session->triangulation();
    const auto signature = added_vertices(mesh);
    require(mesh.all_reference_cells_are_simplex() &&
              mesh.n_vertices() == base_vertex_count + split_count &&
              mesh.n_active_cells() == 80 + 2 * split_count,
            "B2 centroid-split generator violated its count identities");
    require(same_vertices(signature,
                          added_vertices(repeated_session->triangulation())),
            "B2 centroid-split generator is not deterministic for one seed");
    require(!same_vertices(signature,
                           added_vertices(alternate_session->triangulation())),
            "B2 centroid-split generator ignored its selection seed");

    std::size_t boundary_face_count = 0;
    std::size_t fixed_boundary_face_count = 0;
    std::size_t control_boundary_face_count = 0;
    std::size_t outflow_boundary_face_count = 0;
    double      observation_measure = 0.0;
    for (const auto &cell : mesh.active_cell_iterators())
      {
        std::array<unsigned int, 3> vertices{{cell->vertex_index(0),
                                              cell->vertex_index(1),
                                              cell->vertex_index(2)}};
        std::sort(vertices.begin(), vertices.end());
        require(vertices[0] != vertices[1] && vertices[1] != vertices[2],
                "B2 centroid-split mesh contains degenerate connectivity");
        const auto &a = cell->vertex(0);
        const auto &b = cell->vertex(1);
        const auto &c = cell->vertex(2);
        const double signed_area =
          (b[0] - a[0]) * (c[1] - a[1]) -
          (b[1] - a[1]) * (c[0] - a[0]);
        require(signed_area > 0.0,
                "B2 centroid-split mesh contains a misoriented triangle");

        require(cell->material_id() == 0 ||
                  cell->material_id() ==
                    scenario.problem.recipe.observed_material_id,
                "B2 centroid-split mesh has an unknown material id");
        if (cell->material_id() ==
            scenario.problem.recipe.observed_material_id)
          observation_measure += cell->measure();
        for (unsigned int face = 0; face < cell->n_faces(); ++face)
          if (cell->face(face)->at_boundary())
            {
              ++boundary_face_count;
              switch (cell->face(face)->boundary_id())
                {
                  case chapter6::b2_fixed_boundary_id:
                    ++fixed_boundary_face_count;
                    break;
                  case chapter6::b2_control_boundary_id:
                    ++control_boundary_face_count;
                    break;
                  case chapter6::b2_outflow_boundary_id:
                    ++outflow_boundary_face_count;
                    break;
                  default:
                    throw std::runtime_error(
                      "B2 centroid-split mesh has an unclassified boundary face");
                }
            }
      }
    require(boundary_face_count == 28 && fixed_boundary_face_count == 12 &&
              control_boundary_face_count == 6 &&
              outflow_boundary_face_count == 10,
            "B2 centroid splitting changed the exterior partition");
    require(std::abs(observation_measure - 1.8) < 1.0e-12,
            "B2 centroid splitting changed the wings observation measure");

    chapter6::dealii::B2ManufacturedDataT<2> manufactured_data{
      scenario.problem.graetz_case, scenario.problem.fixed_temperature};
    const auto runtime =
      chapter6::dealii::make_b2_manufactured_runtime_data(
        scenario, manufactured_data);
    chapter6::dealii::B2ReducedExecutionAdapterT<2> execute{
      runtime,
      session,
      {"test.chapter6.b2.centroid-split",
       "debug-dealii",
       "test-compiler",
       "test-version",
       "test-standard-library",
       "test-os",
       "test-architecture",
       "test-hardware"}};
    using HeadlessRunner =
      nmopt::application::benchmark::HeadlessBenchmarkRunnerT<
        decltype(scenario)>;
    const auto result = HeadlessRunner(scenario).run(
      [](const auto &parameters) {
        return chapter6::make_b2_problem_spec(parameters);
      },
      execute);
    require(result.artifact.envelope().report().state_solve_count > 0 &&
              result.artifact.envelope().report().adjoint_solve_count > 0,
            "B2 centroid-split execution did not solve state and adjoint systems");
    require(result.document.find("b2.derivative_evidence_passed=true\n") !=
                std::string::npos &&
              result.document.find("benchmark.mesh_vertices=62\n") !=
                std::string::npos &&
              result.document.find("benchmark.mesh_active_cells=94\n") !=
                std::string::npos &&
              std::abs(artifact_number(result.document,
                                       "benchmark.observation_measure") -
                       1.8) < 1.0e-12,
            "B2 centroid-split execution lost numerical or structural evidence");
  }

  void
  run_b2_control_discretisation_comparison()
  {
    using Discretisation =
      nmopt::semantic::v1::NeumannControlDiscretisation;
    using Primal =
      nmopt::contract::PrimalBlockT<chapter6::dealii::Backend>;

    const auto execute = [&](const Discretisation discretisation,
                             const std::size_t    expected_dimension,
                             const std::string &  expected_metric,
                             const std::string &  expected_space) {
      const auto graetz_case =
        chapter6::GraetzCase::observation_full_constant_target;
      auto scenario = chapter6::make_b2_scenario(
        graetz_case,
        nmopt::semantic::v1::TransportBoundaryForm::
          ordinary_normal_minus_transport,
        discretisation);
      scenario.compile.mesh.refinement = 1;
      scenario.experiment.retain_fields = false;
      chapter6::dealii::B2ManufacturedDataT<2> manufactured_data{
        graetz_case, scenario.problem.fixed_temperature};
      const auto runtime =
        chapter6::dealii::make_b2_manufactured_runtime_data(
          scenario, manufactured_data);
      const auto specification = chapter6::make_b2_problem_spec(scenario);
      const auto session =
        chapter6::dealii::make_b2_compilation_session<2>(scenario);
      const auto policy = chapter6::dealii::make_b2_discretisation_policy(
        scenario.compile);
      const auto bindings = chapter6::dealii::make_b2_data_bindings(
        scenario.problem, runtime);
      nmopt::compiler::v1::DealiiCompiler compiler;
      const auto compilation = compiler.compile(
        specification,
        session,
        bindings,
        policy,
        std::nullopt,
        std::nullopt,
        nmopt::compiler::v1::CompilationProduct::reduced_dto);
      require(compilation.succeeded() && compilation.problem,
              "B2 control realization comparison did not compile");

      using Model = nmopt::compiler::v1::detail::
        NeumannBoundaryControlModel<2>;
      const auto *model = dynamic_cast<const Model *>(
        &compilation.problem->executable_model());
      require(model != nullptr &&
                model->physical_control_dimension() == expected_dimension &&
                model->independent_control_dimension() == expected_dimension,
              "B2 control realization produced the wrong dimensions");

      const auto &metric = compilation.problem->metric();
      require(metric.id() == expected_metric &&
                metric.layout()->dimension(0) == expected_dimension,
              "B2 control realization selected the wrong metric");
      ::dealii::Vector<double> control_values(expected_dimension);
      for (std::size_t index = 0; index < expected_dimension; ++index)
        control_values[index] =
          (index % 2 == 0 ? 0.025 : -0.02) *
          static_cast<double>(index + 1);
      const Primal control(metric.layout(), {std::move(control_values)});
      const auto metric_action = metric.apply(control);
      const Primal recovered = metric.inverse_apply(metric_action);
      ::dealii::Vector<double> metric_error = recovered.block(0);
      metric_error.add(-1.0, control.block(0));
      require(metric_error.l2_norm() < 1.0e-10,
              "B2 control metric apply/inverse identity failed");

      const auto reduced = compilation.problem->make_reduced_dto();
      const auto evaluation = reduced.evaluate(control);
      require(evaluation.state_solve.converged() &&
                evaluation.adjoint_solve.converged() &&
                std::isfinite(evaluation.objective_value),
              "B2 control realization did not complete state and adjoint solves");

      const auto &manifest = compilation.problem->manifest();
      const auto control_space = std::find_if(
        manifest.spaces.begin(),
        manifest.spaces.end(),
        [](const auto &space) {
          return space.role == nmopt::semantic::v1::SpaceRole::control;
        });
      require(control_space != manifest.spaces.end() &&
                control_space->dimension == expected_dimension &&
                control_space->finite_element.find(expected_space) !=
                  std::string::npos &&
                manifest.metric_record.realisation_id == expected_metric,
              "B2 control realization manifest is incomplete");
    };

    execute(Discretisation::facewise_constant,
            4,
            "l2_facewise",
            "facewise-constant");
    execute(Discretisation::continuous_nodal_trace,
            6,
            "l2_neumann_trace",
            "continuous scalar degree-one nodal trace");
  }

  void
  run_b2_volume_observation_realisation_comparison()
  {
    using TargetRealisation =
      chapter6::VolumeObservationTargetRealisation;

    const auto execute = [](const chapter6::VolumeObservationOptions options) {
      auto scenario = chapter6::make_b2_scenario(
        chapter6::GraetzCase::observation_full_parabolic_target,
        nmopt::semantic::v1::TransportBoundaryForm::
          ordinary_normal_minus_transport,
        nmopt::semantic::v1::NeumannControlDiscretisation::facewise_constant,
        chapter6::ReducedGlobalization::armijo,
        options);
      scenario.compile.mesh.refinement = 1;
      scenario.solver.parameters.maximum_iterations = 1;
      scenario.solver.parameters.gradient_tolerance = 1.0e-16;
      scenario.experiment.retain_fields = false;

      const auto policy = chapter6::dealii::make_b2_discretisation_policy(
        scenario.compile);
      const auto expected_compiler_target =
        options.target_realisation == TargetRealisation::analytic_quadrature ?
          nmopt::compiler::v1::VolumeObservationTargetRealisation::
            analytic_quadrature :
          nmopt::compiler::v1::VolumeObservationTargetRealisation::
            state_fe_interpolation;
      require(
        policy.volume_observation.has_value() &&
          policy.volume_observation->quadrature_order ==
            options.quadrature_order &&
          policy.volume_observation->target_realisation ==
            expected_compiler_target,
        "B2 adapter did not map its volume-observation policy");

      chapter6::dealii::B2ManufacturedDataT<2> manufactured_data{
        scenario.problem.graetz_case,
        scenario.problem.fixed_temperature};
      const auto runtime =
        chapter6::dealii::make_b2_manufactured_runtime_data(
          scenario, manufactured_data);
      const auto session =
        chapter6::dealii::make_b2_compilation_session<2>(scenario);
      chapter6::dealii::B2ReducedExecutionAdapterT<2> adapter{
        runtime,
        session,
        {"test.chapter6.b2.volume-observation",
         "debug-dealii",
         "test-compiler",
         "test-version",
         "test-standard-library",
         "test-os",
         "test-architecture",
         "test-hardware"}};
      using HeadlessRunner =
        nmopt::application::benchmark::HeadlessBenchmarkRunnerT<
          decltype(scenario)>;
      const auto result = HeadlessRunner(scenario).run(
        [](const auto &parameters) {
          return chapter6::make_b2_problem_spec(parameters);
        },
        adapter);

      const std::string target_name =
        chapter6::volume_observation_target_realisation_name(
          options.target_realisation);
      const std::string order = std::to_string(options.quadrature_order);
      const auto &manifest =
        result.artifact.envelope().compilation_manifest();
      require(
        manifest.observation_realisation.find("target=" + target_name) !=
            std::string::npos &&
          manifest.observation_realisation.find("(" + order + ")") !=
            std::string::npos &&
          result.document.find(
            "b2.volume_observation_quadrature_order=" + order + "\n") !=
            std::string::npos &&
          result.document.find(
            "b2.volume_observation_target_realisation=" + target_name +
            "\n") != std::string::npos &&
          result.document.find("b2.derivative_evidence_passed=true\n") !=
            std::string::npos,
        "B2 execution omitted its selected volume-observation evidence");
      return artifact_number(result.document, "b2.initial_objective");
    };

    const double analytic_objective = execute(
      {3, TargetRealisation::analytic_quadrature});
    const double interpolated_objective = execute(
      {2, TargetRealisation::state_fe_interpolation});
    require(std::abs(analytic_objective - interpolated_objective) > 1.0e-6,
            "B2 observation selector did not affect a parabolic target");
  }

  void
  run_b2_globalization_comparison()
  {
    for (const auto globalization :
         {chapter6::ReducedGlobalization::armijo,
          chapter6::ReducedGlobalization::fixed_step})
      {
        auto scenario = chapter6::make_b2_scenario(
          chapter6::GraetzCase::observation_wings_constant_target,
          nmopt::semantic::v1::TransportBoundaryForm::
            ordinary_normal_minus_transport,
          nmopt::semantic::v1::NeumannControlDiscretisation::facewise_constant,
          globalization);
        scenario.compile.mesh.refinement = 1;
        scenario.solver.parameters.maximum_iterations = 2;
        scenario.solver.parameters.gradient_tolerance = 1.0e-16;
        scenario.solver.parameters.initial_step_length = 0.25;
        scenario.experiment.retain_fields = false;

        chapter6::dealii::B2ManufacturedDataT<2> manufactured_data{
          scenario.problem.graetz_case,
          scenario.problem.fixed_temperature};
        const auto runtime =
          chapter6::dealii::make_b2_manufactured_runtime_data(
            scenario, manufactured_data);
        const auto session =
          chapter6::dealii::make_b2_compilation_session<2>(scenario);
        chapter6::dealii::B2ReducedExecutionAdapterT<2> execute{
          runtime,
          session,
          {"test.chapter6.b2.globalization",
           "debug-dealii",
           "test-compiler",
           "test-version",
           "test-standard-library",
           "test-os",
           "test-architecture",
           "test-hardware"}};
        using HeadlessRunner =
          nmopt::application::benchmark::HeadlessBenchmarkRunnerT<
            decltype(scenario)>;
        const auto result = HeadlessRunner(scenario).run(
          [](const auto &parameters) {
            return chapter6::make_b2_problem_spec(parameters);
          },
          execute);
        const auto &report = result.artifact.envelope().report();
        const std::string expected_policy =
          globalization == chapter6::ReducedGlobalization::armijo ?
            "armijo" :
            "fixed_step";
        const std::string declared_globalization =
          chapter6::reduced_globalization_name(globalization);

        require(report.policy_name == expected_policy &&
                  report.policy_parameters.policy_name == expected_policy &&
                  result.artifact.envelope().solver_policy().policy_name ==
                    expected_policy &&
                  report.parameters.stopping_criterion ==
                    scenario.solver.parameters.stopping_criterion &&
                  report.parameters.maximum_iterations ==
                    scenario.solver.parameters.maximum_iterations,
                "B2 execution did not retain the selected globalization and stopping policy");
        require(result.document.find("solver.globalization=" +
                                     declared_globalization + "\n") !=
                    std::string::npos &&
                  result.document.find("solver.policy=" + expected_policy +
                                       "\n") != std::string::npos,
                "B2 artifact did not retain its declared and effective globalization");
        require(report.accepted_iterations > 0 &&
                  report.step_length_history.size() ==
                    report.accepted_iterations,
                "B2 globalization comparison did not retain accepted steps");
        if (globalization == chapter6::ReducedGlobalization::fixed_step)
          {
            require(report.policy_parameters.maximum_trials == 1 &&
                      report.policy_parameters.initial_step_length == 0.25 &&
                      report.policy_parameters.maximum_backtracking_reductions ==
                        0 &&
                      std::isnan(report.policy_parameters.armijo_fraction) &&
                      report.line_search_trial_count ==
                      report.accepted_iterations &&
                      report.line_search_trials.size() ==
                        report.accepted_iterations,
                    "B2 fixed-step execution did not use one trial per accepted iteration");
            for (std::size_t index = 0;
                 index < report.accepted_iterations;
                 ++index)
              require(report.step_length_history[index] == 0.25 &&
                        report.line_search_trials[index].accepted &&
                        std::isnan(report.iteration_records[index]
                                     .acceptance_evidence
                                     .sufficient_decrease_bound),
                      "B2 fixed-step execution imposed or omitted acceptance evidence");
          }
      }
  }

  void
  run_b2_transport_boundary_realisation_comparison()
  {
    const auto graetz_case =
      chapter6::GraetzCase::observation_full_constant_target;

    const auto zero_control_state_linf =
      [&](const nmopt::semantic::v1::TransportBoundaryForm boundary_form) {
        auto scenario = chapter6::make_b2_scenario(graetz_case, boundary_form);
        scenario.compile.mesh.refinement = 2;
        chapter6::dealii::B2ManufacturedDataT<2> manufactured_data{
          graetz_case, scenario.problem.fixed_temperature};
        const auto runtime =
          chapter6::dealii::make_b2_manufactured_runtime_data(
            scenario, manufactured_data);
        const auto specification = chapter6::make_b2_problem_spec(scenario);
        const auto session =
          chapter6::dealii::make_b2_compilation_session<2>(scenario);
        const auto policy = chapter6::dealii::make_b2_discretisation_policy(
          scenario.compile);
        const auto bindings =
          chapter6::dealii::make_b2_data_bindings(scenario.problem, runtime);
        nmopt::compiler::v1::DealiiCompiler compiler;
        const auto compilation = compiler.compile(specification,
                                                  session,
                                                  bindings,
                                                  policy,
                                                  std::nullopt,
                                                  std::nullopt,
                                                  nmopt::compiler::v1::CompilationProduct::reduced_dto);
        require(compilation.succeeded() && compilation.problem,
                "B2 transport realization comparison did not compile");
        require(
          compilation.problem->manifest().boundary_realisation.has_value() &&
            compilation.problem->manifest()
                .boundary_realisation->transport_boundary_form == boundary_form,
          "B2 transport realization did not propagate through the application boundary");
        const auto reduced = compilation.problem->make_reduced_dto();
        const nmopt::contract::StateControlPartitionT<
          chapter6::dealii::Backend>
          partition(compilation.problem->executable_model(), 0, 1);
        const auto control =
          nmopt::contract::PrimalBlockT<chapter6::dealii::Backend>::zeros(
            partition.control_layout());
        const auto evaluation = reduced.evaluate(control);
        double state_linf = 0.0;
        for (unsigned int index = 0;
             index < evaluation.state.block(0).size();
             ++index)
          state_linf = std::max(state_linf,
                                std::abs(evaluation.state.block(0)[index]));
        require(std::isfinite(state_linf) &&
                  std::isfinite(evaluation.objective_value),
                "B2 transport realization comparison produced non-finite output");
        return state_linf;
      };

    const double total_conormal_linf = zero_control_state_linf(
      nmopt::semantic::v1::TransportBoundaryForm::total_conormal);
    const double ordinary_normal_linf = zero_control_state_linf(
      nmopt::semantic::v1::TransportBoundaryForm::
        ordinary_normal_minus_transport);
    require(total_conormal_linf > 10.0 * ordinary_normal_linf,
            "B2 transport realizations did not separate the forward state scale");
  }
} // namespace

int
main(const int argc, char **argv)
{
  try
    {
      const std::vector<nmopt::test_support::Scenario> scenarios{
        {"b2_manufactured_wings_constant",
         "nmopt.application.dealii.b2_manufactured_wings_constant",
         {"dealii", "application", "benchmark", "b2", "contract"},
         120,
         []() {
           run_b2_manufactured_case(
             nmopt::application::chapter6::GraetzCase::
               observation_wings_constant_target);
         }},
        {"b2_manufactured_wings_parabolic",
         "nmopt.application.dealii.b2_manufactured_wings_parabolic",
         {"dealii", "application", "benchmark", "b2", "contract"},
         120,
         []() {
           run_b2_manufactured_case(
             nmopt::application::chapter6::GraetzCase::
               observation_wings_parabolic_target);
         }},
        {"b2_manufactured_full_constant",
         "nmopt.application.dealii.b2_manufactured_full_constant",
         {"dealii", "application", "benchmark", "b2", "contract"},
         120,
         []() {
           run_b2_manufactured_case(
             nmopt::application::chapter6::GraetzCase::
               observation_full_constant_target);
         }},
        {"b2_manufactured_full_parabolic",
         "nmopt.application.dealii.b2_manufactured_full_parabolic",
         {"dealii", "application", "benchmark", "b2", "contract"},
         120,
         []() {
           run_b2_manufactured_case(
             nmopt::application::chapter6::GraetzCase::
               observation_full_parabolic_target);
         }},
        {"b2_transport_boundary_realisation_comparison",
         "nmopt.application.dealii.b2_transport_boundary_realisation_comparison",
         {"dealii", "application", "benchmark", "b2", "contract", "diagnostic"},
         120,
         []() { run_b2_transport_boundary_realisation_comparison(); }},
        {"b2_globalization_comparison",
         "nmopt.application.dealii.b2_globalization_comparison",
         {"dealii", "application", "benchmark", "b2", "contract", "solver"},
         120,
         []() { run_b2_globalization_comparison(); }},
        {"b2_control_discretisation_comparison",
         "nmopt.application.dealii.b2_control_discretisation_comparison",
         {"dealii", "application", "benchmark", "b2", "contract", "metric"},
         120,
         []() { run_b2_control_discretisation_comparison(); }},
        {"b2_volume_observation_realisation_comparison",
         "nmopt.application.dealii.b2_volume_observation_realisation_comparison",
         {"dealii", "application", "benchmark", "b2", "contract", "diagnostic"},
         120,
         []() { run_b2_volume_observation_realisation_comparison(); }},
        {"b2_structured_simplex_mesh",
         "nmopt.application.dealii.b2_structured_simplex_mesh",
         {"dealii", "application", "benchmark", "b2", "contract"},
         120,
         []() { run_b2_structured_simplex_mesh(); }},
        {"b2_centroid_split_simplex_mesh",
         "nmopt.application.dealii.b2_centroid_split_simplex_mesh",
         {"dealii", "application", "benchmark", "b2", "contract"},
         120,
         []() { run_b2_centroid_split_simplex_mesh(); }}};
      const auto result = nmopt::test_support::run_requested_scenarios(
        argc, argv, scenarios, std::cout);
      if (!result.listed)
        std::cout << "chapter 6 B2 deal.II contract test passed: "
                  << result.executed << '\n';
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "chapter 6 B2 deal.II contract test failed: "
                << exception.what() << '\n';
      return 1;
    }
}
