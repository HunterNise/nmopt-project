#include "nmopt/application/dealii/chapter6_b2.hpp"
#include "../support/scenario_dispatch.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
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
    require(result.document.find("b2.graetz_case=") != std::string::npos,
            "B2 deal.II adapter omitted case evidence");
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
    std::filesystem::remove_all(native_output_directory);
  }

  void
  run_b2_transport_boundary_realisation_comparison()
  {
    const auto graetz_case =
      chapter6::GraetzCase::observation_full_constant_target;
    auto scenario = chapter6::make_b2_scenario(graetz_case);
    scenario.compile.mesh.refinement = 2;
    chapter6::dealii::B2ManufacturedDataT<2> manufactured_data{
      graetz_case, scenario.problem.fixed_temperature};
    const auto runtime =
      chapter6::dealii::make_b2_manufactured_runtime_data(
        scenario, manufactured_data);
    const auto specification = chapter6::make_b2_problem_spec(scenario);

    const auto zero_control_state_linf =
      [&](const nmopt::compiler::v1::TransportBoundaryRealisation realisation) {
        const auto session =
          chapter6::dealii::make_b2_compilation_session<2>(scenario);
        auto policy = chapter6::dealii::make_b2_discretisation_policy(
          scenario.compile);
        policy.transport_boundary_realisation = realisation;
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
      nmopt::compiler::v1::TransportBoundaryRealisation::total_conormal);
    const double ordinary_normal_linf = zero_control_state_linf(
      nmopt::compiler::v1::TransportBoundaryRealisation::ordinary_normal_transport);
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
         []() { run_b2_transport_boundary_realisation_comparison(); }}};
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
