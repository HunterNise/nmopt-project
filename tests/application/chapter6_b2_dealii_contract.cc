#include "nmopt/application/dealii/chapter6_b2.hpp"
#include "../support/scenario_dispatch.hpp"

#include <algorithm>
#include <cmath>
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

    chapter6::dealii::B2ManufacturedDataT<2> manufactured_data{
      graetz_case, scenario.problem.fixed_temperature};
    const auto runtime =
      chapter6::dealii::make_b2_manufactured_runtime_data(
        scenario, manufactured_data);
    const auto session =
      chapter6::dealii::make_b2_compilation_session<2>(scenario);

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
       "test-hardware"}};

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
    require(result.document.find("solver.method=bfgs\n") !=
              std::string::npos,
            "B2 deal.II adapter omitted solver-method evidence");
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
         }}};
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
