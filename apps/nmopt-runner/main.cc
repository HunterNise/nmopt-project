#include "nmopt/application/application.hpp"
#include "nmopt/application/dealii/chapter6_b1.hpp"
#include "nmopt/application/dealii/chapter6_b2.hpp"
#include "runner.hpp"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef NMOPT_COMPILED_BUILD_PROFILE
#  define NMOPT_COMPILED_BUILD_PROFILE "unknown"
#endif

namespace
{
  const char *
  entry_kind_name(const nmopt::application::CatalogEntryKind kind)
  {
    switch (kind)
      {
        case nmopt::application::CatalogEntryKind::recipe:
          return "recipe";
        case nmopt::application::CatalogEntryKind::scenario:
          return "scenario";
      }
    return "unknown";
  }

  std::string
  join_requirements(const std::vector<std::string> &requirements)
  {
    std::string result;
    for (std::size_t index = 0; index < requirements.size(); ++index)
      {
        if (index != 0)
          result += ',';
        result += requirements[index];
      }
    return result;
  }

  void
  print_usage(std::ostream &output)
  {
    output << "Usage: nmopt_runner --list [--output DIRECTORY]\n"
           << "       nmopt_runner --benchmark b1 --framework-revision REV"
              " [--output DIRECTORY] [--refinement N]\n"
           << "       nmopt_runner --benchmark b2 --framework-revision REV"
              " [--output DIRECTORY] [--refinement N]\n"
           << "       nmopt_runner --help\n"
           << "\n"
           << "--list             list registered Chapter 5/6 application entries\n"
           << "--benchmark b1     run the frozen B1 matrix (six artifacts)\n"
           << "--benchmark b2     run the frozen B2 case batch (four artifacts)\n"
           << "--output DIRECTORY set the runner-owned artifact root (default: runs)\n"
           << "--framework-revision REV\n"
           << "                   record the framework revision in each artifact\n"
           << "--refinement N     override the mesh refinement (default: 7)\n"
           << "--help             show this message\n";
  }

  void
  print_catalog(const nmopt::application::ApplicationCatalog &catalog,
                const std::filesystem::path &output_directory)
  {
    std::cout << "catalog.schema=nmopt-application-v1\n"
              << "catalog.output_directory=" << output_directory.string()
              << '\n';
    for (const auto &entry : catalog.entries())
      std::cout << entry_kind_name(entry.kind) << '\t' << entry.id << '\t'
                << entry.recipe_id << '\t' << entry.chapter << '\t'
                << entry.label << '\t'
                << join_requirements(entry.requirements) << '\n';
  }

  const char *
  b1_method_slug(const nmopt::application::chapter6::ReducedMethod method)
  {
    switch (method)
      {
        case nmopt::application::chapter6::ReducedMethod::steepest_descent:
          return "steepest-descent";
        case nmopt::application::chapter6::ReducedMethod::limited_memory_bfgs:
          return "l-bfgs";
        case nmopt::application::chapter6::ReducedMethod::bfgs:
          break;
      }
    throw std::invalid_argument("B1 has no artifact slug for this method");
  }

  const char *
  b1_beta_slug(const double beta)
  {
    if (std::abs(beta - 1.0e-1) < 1.0e-15)
      return "1e-1";
    if (std::abs(beta - 1.0e-2) < 1.0e-15)
      return "1e-2";
    if (std::abs(beta - 1.0e-3) < 1.0e-15)
      return "1e-3";
    throw std::invalid_argument(
      "B1 has no frozen artifact slug for the regularisation value");
  }

  std::string
  b1_mesh_provenance(const unsigned int refinement)
  {
    return "chapter-6.e6.5.1.framework-native-hypercube-r" +
           std::to_string(refinement);
  }

  std::string
  b2_mesh_provenance(const unsigned int refinement)
  {
    return "chapter-6.e6.5.2.framework-native-rectangle-r" +
           std::to_string(refinement);
  }

  std::string
  b1_number(const double value)
  {
    std::ostringstream output;
    output.precision(17);
    output << value;
    return output.str();
  }

  template <typename Scenario>
  nmopt::experiment::RunEnvironmentRecord
  make_environment(const Scenario &scenario)
  {
#if defined(__clang__)
    const std::string compiler = "clang";
#elif defined(__GNUC__)
    const std::string compiler = "gcc";
#elif defined(_MSC_VER)
    const std::string compiler = "msvc";
#else
    const std::string compiler = "unknown";
#endif

#if defined(__x86_64__) || defined(_M_X64)
    const std::string architecture = "x86_64";
#elif defined(__aarch64__) || defined(_M_ARM64)
    const std::string architecture = "aarch64";
#else
    const std::string architecture = "unknown";
#endif

    const char *hardware = std::getenv("NMOPT_HARDWARE");
    return {scenario.experiment.source_revision,
            NMOPT_COMPILED_BUILD_PROFILE,
            compiler,
#if defined(__VERSION__)
            __VERSION__,
#else
            "unknown",
#endif
            "libstdc++",
            "Linux",
            architecture,
            hardware == nullptr ? "unspecified" : hardware};
  }

  using Chapter6Evidence =
    nmopt::application::benchmark::BenchmarkExecutionEvidenceT<
      nmopt::application::chapter6::dealii::Envelope>;

  void
  add_common_artifact_fields(Chapter6Evidence &evidence,
                             const std::string &framework_revision)
  {
    const auto &manifest = evidence.envelope.compilation_manifest();
    const auto &environment = evidence.envelope.environment();
    evidence.fields.push_back({"provenance.framework_revision",
                               framework_revision});
    evidence.fields.push_back({
      "provenance.recipe_revision",
      "application-chapter6@" + framework_revision});
    evidence.fields.push_back({"provenance.mesh_provenance",
                               manifest.mesh_record.provenance});
    evidence.fields.push_back({"provenance.mesh_structural_identity",
                               manifest.mesh_record.structural_identity});
    evidence.fields.push_back({"provenance.compiler", environment.compiler});
    evidence.fields.push_back({"provenance.compiler_version",
                               environment.compiler_version});
    evidence.fields.push_back({"provenance.standard_library",
                               environment.standard_library});
    evidence.fields.push_back({"provenance.operating_system",
                               environment.operating_system});
    evidence.fields.push_back({"provenance.architecture",
                               environment.architecture});
    evidence.fields.push_back({"provenance.hardware", environment.hardware});
    evidence.fields.push_back({"manifest.schema_version",
                               std::to_string(manifest.schema_version)});
    evidence.fields.push_back({"manifest.semantic_problem_id",
                               manifest.semantic_problem_id});
    evidence.fields.push_back({"manifest.compiler_id", manifest.compiler_id});
    evidence.fields.push_back({"manifest.backend", manifest.backend});
    evidence.fields.push_back({"manifest.execution", manifest.execution});
    evidence.fields.push_back({"manifest.state_space", manifest.state_space});
    evidence.fields.push_back({"manifest.control_space",
                               manifest.control_space});
    evidence.fields.push_back({"manifest.mesh_dimension",
                               std::to_string(manifest.mesh_record.dimension)});
    evidence.fields.push_back({"manifest.active_cells",
                               std::to_string(manifest.mesh_record.active_cells)});
    evidence.fields.push_back({"manifest.space_count",
                               std::to_string(manifest.spaces.size())});
    evidence.fields.push_back({"manifest.binding_count",
                               std::to_string(manifest.bindings.size())});
    evidence.fields.push_back({"manifest.realized_space_count",
                               std::to_string(manifest.realized_spaces.size())});
    evidence.fields.push_back({"manifest.realized_map_count",
                               std::to_string(manifest.realized_maps.size())});
    evidence.fields.push_back({
      "solver.final_objective",
      b1_number(evidence.envelope.report().final_evaluation.objective_value)});
    const auto &parameters = evidence.envelope.report().parameters;
    evidence.fields.push_back({"solver.maximum_iterations",
                               std::to_string(parameters.maximum_iterations)});
    evidence.fields.push_back({"solver.maximum_line_search_trials",
                               std::to_string(parameters.maximum_line_search_trials)});
    evidence.fields.push_back({"solver.gradient_tolerance",
                               b1_number(parameters.gradient_tolerance)});
    evidence.fields.push_back({"solver.initial_step_length",
                               b1_number(parameters.initial_step_length)});
    evidence.fields.push_back({"solver.armijo_fraction",
                               b1_number(parameters.armijo_fraction)});
    evidence.fields.push_back({"solver.backtracking_factor",
                               b1_number(parameters.backtracking_factor)});
  }

  void
  add_b1_artifact_fields(
    nmopt::application::benchmark::BenchmarkExecutionEvidenceT<
      nmopt::application::chapter6::dealii::Envelope> &evidence,
    const nmopt::application::chapter6::B1Scenario &scenario,
    const nmopt::application::chapter6::ReducedMethod method,
    const double beta,
    const std::string &framework_revision)
  {
    add_common_artifact_fields(evidence, framework_revision);
    evidence.fields.push_back({"benchmark.method", b1_method_slug(method)});
    evidence.fields.push_back({"benchmark.regularisation", b1_number(beta)});
    evidence.fields.push_back(
      {"benchmark.mesh_refinement", std::to_string(scenario.compile.mesh.refinement)});
    evidence.fields.push_back({
      "benchmark.mesh_mode",
      scenario.compile.mesh.refinement == 7 ? "source-scale-native"
                                            : "development"});
    evidence.fields.push_back({"provenance.forcing",
                               scenario.problem.data.forcing_provenance});
    evidence.fields.push_back({"provenance.desired_state",
                               scenario.problem.data.desired_state_provenance});
  }

  void
  add_b2_artifact_fields(
    Chapter6Evidence &evidence,
    const nmopt::application::chapter6::B2Scenario &scenario,
    const std::string &framework_revision)
  {
    const auto  case_slug = nmopt::application::chapter6::graetz_case_name(
      scenario.problem.graetz_case);
    evidence.fields.push_back({"benchmark.graetz_case", case_slug});
    evidence.fields.push_back({"benchmark.fixed_temperature",
                               b1_number(scenario.problem.fixed_temperature)});
    evidence.fields.push_back({"benchmark.regularisation",
                               b1_number(
                                 scenario.problem.data.regularisation_weight)});
    evidence.fields.push_back(
      {"benchmark.mesh_refinement", std::to_string(scenario.compile.mesh.refinement)});
    evidence.fields.push_back({
      "benchmark.mesh_mode",
      scenario.compile.mesh.refinement == 7 ? "source-scale-native"
                                            : "development"});
    add_common_artifact_fields(evidence, framework_revision);
    evidence.fields.push_back({"provenance.forcing",
                               scenario.problem.data.forcing_provenance});
    evidence.fields.push_back({"provenance.desired_state",
                               scenario.problem.data.desired_state_provenance});
    evidence.fields.push_back({"provenance.fixed_temperature",
                               scenario.problem.data.fixed_dirichlet_data_provenance});
    evidence.fields.push_back({"provenance.conservative_transport",
                               scenario.problem.data.conservative_transport_provenance});
    evidence.fields.push_back({"provenance.observation_case", case_slug});
  }

  void
  write_artifact(const std::filesystem::path &path, const std::string &document)
  {
    nmopt::application::runner::prepare_artifact_path(path);
    std::ofstream output(path);
    if (!output)
      throw std::runtime_error("could not open benchmark artifact '" +
                               path.string() + "'");
    output << document;
    if (!output)
      throw std::runtime_error("could not write benchmark artifact '" +
                               path.string() + "'");
  }

  template <typename Report>
  void
  write_solver_trace(const std::filesystem::path &path, const Report &report)
  {
    nmopt::application::runner::prepare_artifact_path(path);
    std::ofstream output(path);
    if (!output)
      throw std::runtime_error("could not open solver trace '" +
                               path.string() + "'");
    output.imbue(std::locale::classic());
    output << "iteration,trial,step_length,objective_value,actual_slope,"
               "sufficient_decrease_bound,objective_finite,slope_negative,"
               "accepted\n";
    output.precision(17);
    for (const auto &trial : report.line_search_trials)
      output << trial.iteration << ',' << trial.trial << ','
             << trial.step_length << ',' << trial.objective_value << ','
             << trial.actual_slope << ','
             << trial.sufficient_decrease_bound << ','
             << (trial.objective_finite ? "true" : "false") << ','
             << (trial.slope_negative ? "true" : "false") << ','
             << (trial.accepted ? "true" : "false") << '\n';
    if (!output)
      throw std::runtime_error("could not write solver trace '" +
                               path.string() + "'");
  }

  void
  run_b1(const nmopt::application::runner::CommandLineOptions &options)
  {
    using namespace nmopt::application;
    using namespace chapter6;
    using Runner = benchmark::HeadlessBenchmarkRunnerT<B1Scenario>;
    using Adapter =
      nmopt::application::chapter6::dealii::B1ReducedExecutionAdapterT<2>;

    const std::vector<ReducedMethod> methods{
      ReducedMethod::steepest_descent,
      ReducedMethod::limited_memory_bfgs};
    for (const auto method : methods)
      {
        const auto method_slug = b1_method_slug(method);
        const auto base_scenario = make_b1_scenario(method);
        for (const double beta : base_scenario.problem.regularisation_sweep)
          {
            auto scenario = base_scenario;
            scenario.experiment.build_profile = NMOPT_COMPILED_BUILD_PROFILE;
            const auto beta_slug = b1_beta_slug(beta);
            scenario.compile.mesh.refinement = options.refinement;
            scenario.compile.mesh.mesh_provenance =
              b1_mesh_provenance(options.refinement);
            scenario.experiment.harness.artifact_directory =
              options.output_directory.string();
            scenario.experiment.scenario_output_id =
              "chapter-6.b1.distributed-laplace." +
              std::string(method_slug) + ".beta-" + beta_slug;

            nmopt::application::chapter6::dealii::B1ManufacturedDataT<2> data;
            const auto runtime = data.runtime_data();
            const auto session =
              nmopt::application::chapter6::dealii::
                make_b1_compilation_session<2>(scenario);
            const auto environment = make_environment(scenario);
            const auto path = runner::artifact_path(
              options.output_directory,
              {"chapter-6.b1.distributed-laplace",
               method_slug,
               "beta-" + std::string(beta_slug)});
            Adapter execute{beta,
                            runtime,
                            session,
                            environment,
                            path.parent_path()};
            Runner runner(scenario);
            const auto result = runner.run(
              [](const auto &parameters) {
                return chapter6::make_b1_problem_spec(parameters);
              },
              [&](const auto &specification, const auto &run_scenario) {
                auto evidence = execute(specification, run_scenario);
                add_b1_artifact_fields(evidence,
                                       run_scenario,
                                       method,
                                       beta,
                                       options.framework_revision);
                return evidence;
              });

            write_artifact(path, result.document);
            write_solver_trace(path.parent_path() / "solver-trace.csv",
                               result.artifact.envelope().report());
            std::cout << "B1 wrote " << path.string() << '\n';
          }
      }
  }

  void
  run_b2(const nmopt::application::runner::CommandLineOptions &options)
  {
    using namespace nmopt::application;
    using namespace chapter6;
    using Runner = benchmark::HeadlessBenchmarkRunnerT<B2Scenario>;
    using Adapter =
      nmopt::application::chapter6::dealii::B2ReducedExecutionAdapterT<2>;

    for (const auto graetz_case : b2_case_order)
      {
        auto scenario = make_b2_scenario(graetz_case);
        scenario.experiment.build_profile = NMOPT_COMPILED_BUILD_PROFILE;
        const auto case_slug = graetz_case_name(graetz_case);
        scenario.compile.mesh.refinement = options.refinement;
        scenario.compile.mesh.mesh_provenance =
          b2_mesh_provenance(options.refinement);
        scenario.experiment.harness.artifact_directory =
          options.output_directory.string();
        scenario.experiment.scenario_output_id =
          "chapter-6.b2.graetz-flow." + std::string(case_slug);

        nmopt::application::chapter6::dealii::B2ManufacturedDataT<2> data(
          graetz_case,
          scenario.problem.fixed_temperature);
        const auto runtime =
          nmopt::application::chapter6::dealii::
            make_b2_manufactured_runtime_data<2>(scenario, data);
        const auto session =
          nmopt::application::chapter6::dealii::
            make_b2_compilation_session<2>(scenario);
        const auto environment = make_environment(scenario);
        const auto path = runner::artifact_path(
          options.output_directory,
          {"chapter-6.b2.graetz-flow", case_slug});
        Adapter execute{runtime, session, environment, path.parent_path()};
        Runner runner(scenario);
        const auto result = runner.run(
          [](const auto &parameters) {
            return chapter6::make_b2_problem_spec(parameters);
          },
          [&](const auto &specification, const auto &run_scenario) {
            auto evidence = execute(specification, run_scenario);
            add_b2_artifact_fields(evidence,
                                   run_scenario,
                                   options.framework_revision);
            return evidence;
          });

        write_artifact(path, result.document);
        write_solver_trace(path.parent_path() / "solver-trace.csv",
                           result.artifact.envelope().report());
        std::cout << "B2 wrote " << path.string() << '\n';
      }
  }
} // namespace

int
main(const int argc, char **argv)
{
  try
    {
      const auto options =
        nmopt::application::runner::parse_command_line(argc, argv);
      if (options.help)
        {
          print_usage(std::cout);
          return 0;
        }

      if (options.run_b1)
        {
          run_b1(options);
          return 0;
        }

      if (options.run_b2)
        {
          run_b2(options);
          return 0;
        }

      nmopt::application::ApplicationCatalog catalog;
      const auto chapter5 =
        nmopt::application::chapter5::make_catalog();
      for (const auto &entry : chapter5.entries())
        catalog.add(entry);
      const auto chapter6 =
        nmopt::application::chapter6::make_catalog();
      for (const auto &entry : chapter6.entries())
        catalog.add(entry);

      print_catalog(catalog, options.output_directory);
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "nmopt_runner: " << exception.what() << '\n';
      print_usage(std::cerr);
      return 2;
    }
}
