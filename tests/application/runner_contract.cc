#include "../../apps/nmopt-runner/runner.hpp"
#include "../../apps/nmopt-runner/benchmark_registry.hpp"
#include "../../apps/nmopt-runner/capability_registry.hpp"
#include "nmopt/application/chapter6.hpp"
#include "../support/scenario_dispatch.hpp"

#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
  using nmopt::application::runner::CommandLineOptions;
  using nmopt::application::runner::RunKind;

  void
  require(const bool condition, const char *message)
  {
    if (!condition)
      throw std::runtime_error(message);
  }

  CommandLineOptions
  parse(const std::vector<std::string> &arguments)
  {
    std::vector<std::string> mutable_arguments = arguments;
    std::vector<char *>      argv;
    argv.reserve(mutable_arguments.size());
    for (auto &argument : mutable_arguments)
      argv.push_back(argument.data());

    return nmopt::application::runner::parse_command_line(
      static_cast<int>(argv.size()), argv.data());
  }

  void
  require_invalid_argument(const std::function<void()> &operation,
                           const char                    *message)
  {
    try
      {
        operation();
      }
    catch (const std::invalid_argument &)
      {
        return;
      }
    throw std::runtime_error(message);
  }

  void
  test_reproduction_policy()
  {
    const auto *b1 =
      nmopt::application::runner::find_benchmark_registration("b1");
    const auto *b2 =
      nmopt::application::runner::find_benchmark_registration("b2");
    require(b1 != nullptr && b2 != nullptr,
            "B1 and B2 should be present in the runner-local registry");
    require(b1->parameter_benchmark_id ==
              "chapter-6.b1.distributed-laplace" &&
              b2->parameter_benchmark_id == "chapter-6.b2.graetz-flow",
            "runner registrations should bind CLI IDs to parameter benchmark IDs");
    require(nmopt::application::runner::find_benchmark_registration("b3") ==
              nullptr,
            "unregistered benchmark IDs should not resolve");
    require(nmopt::application::runner::registered_benchmark_ids() == "b1, b2",
            "runner registry should expose its available benchmark IDs");

    auto options = parse({"nmopt_runner",
                          "--benchmark",
                          "b1",
                          "--framework-revision",
                          "source-revision"});
    require(options.benchmark == "b1",
            "benchmark selection should retain one benchmark identifier");
    require(options.run_kind == RunKind::reproduction,
            "benchmark runs should default to reproduction kind");
    require(!options.refinement_override.has_value(),
            "reproduction should use the benchmark mesh default when unset");
    require(nmopt::application::chapter6::make_b1_scenario()
              .compile.mesh.refinement ==
              nmopt::application::chapter6::b1_default_mesh_refinement,
            "B1 should supply its own mesh default");
    require(nmopt::application::chapter6::make_b2_scenario()
              .compile.mesh.refinement ==
              nmopt::application::chapter6::b2_default_mesh_refinement,
            "B2 should supply its own mesh default");
    nmopt::application::runner::validate_run_policy(options, "release-dealii");
    const auto reproduction =
      nmopt::application::runner::resolve_run_configuration(
        options, "release-dealii");
    require(
        reproduction.run_directory ==
        std::filesystem::path(
          "runs/chapter-6/b1/authoritative"),
      "reproduction runs should use the organized run-set layout");
    require(
      nmopt::application::runner::artifact_path(
        reproduction.run_directory, {"sample"}) ==
        reproduction.run_directory / "artifacts" / "sample" / "artifact.kv",
      "artifact paths should be nested below the run-set artifacts directory");

    const auto alternate_mesh = parse({"nmopt_runner",
                                       "--benchmark",
                                       "b1",
                                       "--framework-revision",
                                       "source-revision",
                                       "--refinement",
                                       "1"});
    require(alternate_mesh.refinement_override == 1,
            "an explicit mesh override was not retained");
    nmopt::application::runner::validate_run_policy(alternate_mesh,
                                                    "release-dealii");

    require_invalid_argument(
      [&] {
        nmopt::application::runner::validate_run_policy(options,
                                                        "debug-dealii");
      },
      "a non-release build should not pass reproduction policy");

    auto development = parse({"nmopt_runner",
                               "--benchmark",
                               "b2",
                               "--framework-revision",
                               "source-revision",
                               "--run-kind",
                               "development",
                               "--refinement",
                               "1"});
    require(development.benchmark == "b2",
            "development benchmark selection should retain its identifier");
    require(development.run_kind == RunKind::development,
            "the development run kind was not parsed");
    nmopt::application::runner::validate_run_policy(development,
                                                    "debug-dealii");
    const auto temporary_root =
      std::filesystem::temp_directory_path() /
      "nmopt-runner-contract-layout";
    std::filesystem::remove_all(temporary_root);
    std::filesystem::create_directories(
      temporary_root / "chapter-6" / "b2" / "development" / "001");
    development.output_directory = temporary_root;
    const auto development_configuration =
      nmopt::application::runner::resolve_run_configuration(
        development, "debug-dealii");
    require(
      development_configuration.run_directory ==
        temporary_root / "chapter-6" / "b2" / "development" / "002",
      "development runs should use the organized run-set layout");
    std::filesystem::remove_all(temporary_root);

    auto versioned = parse({"nmopt_runner",
                            "--benchmark",
                            "b2",
                            "--framework-revision",
                            "source-revision",
                            "--run-kind",
                            "development",
                            "--run-slot",
                            "004-v2"});
    require(versioned.run_slot.has_value() && *versioned.run_slot == "004-v2",
            "a named development slot was not parsed");
    versioned.output_directory = temporary_root;
    const auto versioned_configuration =
      nmopt::application::runner::resolve_run_configuration(
        versioned, "debug-dealii");
    require(
      versioned_configuration.run_directory ==
        temporary_root / "chapter-6" / "b2" / "development" / "004-v2",
      "named development runs should use the requested versioned slot");

    require_invalid_argument(
      [] {
        parse({"nmopt_runner",
               "--benchmark",
               "b2",
               "--framework-revision",
               "source-revision",
               "--run-slot",
               "004/invalid"});
      },
      "development slots should reject path separators");

    const auto revision_options = parse({"nmopt_runner",
                                         "--benchmark",
                                         "b1",
                                         "--framework-revision",
                                         "feature/branch"});
    const auto revision_configuration =
      nmopt::application::runner::resolve_run_configuration(
        revision_options, "release-dealii");
    require(revision_configuration.framework_revision == "feature/branch",
            "framework revisions should remain metadata rather than paths");

    const auto future_benchmark = parse({"nmopt_runner",
                                         "--benchmark",
                                         "b3",
                                         "--framework-revision",
                                         "source-revision"});
    require(future_benchmark.benchmark == "b3",
            "CLI parsing should not own a closed B1/B2 benchmark list");
  }

  void
  test_typed_capability_registries()
  {
    using nmopt::application::chapter6::ExecutionSelection;
    using nmopt::application::chapter6::ProductSelection;
    using nmopt::application::chapter6::ReducedMethod;
    using nmopt::application::runner::execution_capability_registry;
    using nmopt::application::runner::product_capability_registry;
    using nmopt::application::runner::reduced_method_capability_registry;

    require(product_capability_registry().resolve("reduced-dto", "product") ==
              ProductSelection::reduced_dto,
            "product capability lookup should return the typed selection");
    require(execution_capability_registry().resolve("assembled", "execution") ==
              ExecutionSelection::assembled,
            "execution capability lookup should return the typed selection");
    require(reduced_method_capability_registry().resolve("bfgs", "method") ==
              ReducedMethod::bfgs,
            "method capability lookup should return the typed selection");
    require(reduced_method_capability_registry().find("not-a-method") == nullptr,
            "unknown capability IDs should not resolve through find");
    require_invalid_argument(
      [] {
        reduced_method_capability_registry().resolve("not-a-method", "method");
      },
      "unknown capability IDs should fail typed resolution");
  }

  void
  test_run_kind_parser_rejects_unknown_values()
  {
    require_invalid_argument(
      [] {
        parse({"nmopt_runner",
               "--benchmark",
               "b1",
               "--framework-revision",
               "source-revision",
               "--run-kind",
               "investigation"});
      },
      "an unknown run kind should be rejected");
  }

  void
  test_matched_reference_stopping_criteria()
  {
    using Criterion = nmopt::solvers::ReducedStoppingCriterion;
    using Reason = nmopt::solvers::ReducedStoppingReason;
    using nmopt::application::runner::reference_reached_stopping_tolerance;

    require(reference_reached_stopping_tolerance(
              Criterion::automatic, Reason::gradient_tolerance) &&
              reference_reached_stopping_tolerance(
                Criterion::automatic, Reason::relative_gradient_tolerance) &&
              reference_reached_stopping_tolerance(
                Criterion::automatic, Reason::objective_change_tolerance) &&
              reference_reached_stopping_tolerance(
                Criterion::automatic, Reason::step_tolerance),
            "automatic matched references should accept every tolerance stop");
    require(reference_reached_stopping_tolerance(
              Criterion::gradient_norm, Reason::gradient_tolerance) &&
              reference_reached_stopping_tolerance(
                Criterion::relative_gradient_norm,
                Reason::relative_gradient_tolerance) &&
              reference_reached_stopping_tolerance(
                Criterion::objective_change,
                Reason::objective_change_tolerance) &&
              reference_reached_stopping_tolerance(
                Criterion::step_norm, Reason::step_tolerance) &&
              reference_reached_stopping_tolerance(
                Criterion::objective_change, Reason::stationary) &&
              reference_reached_stopping_tolerance(
                Criterion::step_norm, Reason::stationary),
            "matched references should accept their selected successful stop");
    require(!reference_reached_stopping_tolerance(
              Criterion::gradient_norm,
              Reason::relative_gradient_tolerance) &&
              !reference_reached_stopping_tolerance(
                Criterion::relative_gradient_norm,
                Reason::gradient_tolerance) &&
              !reference_reached_stopping_tolerance(
                Criterion::automatic, Reason::maximum_iterations),
            "matched references should reject unrelated or failed stops");
  }

  std::string
  read_file(const std::filesystem::path &path)
  {
    std::ifstream input(path);
    if (!input)
      throw std::runtime_error("could not read test file");
    return {std::istreambuf_iterator<char>(input),
            std::istreambuf_iterator<char>()};
  }

  void
  test_run_manifest_records_state_and_failures()
  {
    using nmopt::application::runner::ResolvedRunConfiguration;
    using nmopt::application::runner::RunKind;
    using nmopt::application::runner::RunSetManifest;

    const auto temporary_root =
      std::filesystem::temp_directory_path() / "nmopt-run-manifest-contract";
    std::filesystem::remove_all(temporary_root);
    auto configuration = ResolvedRunConfiguration{
      temporary_root,
      temporary_root / "chapter-6" / "b1" / "development" / "001",
      "b1",
      "debug-dealii",
      "test-revision",
      RunKind::development,
      1};
    configuration.parameter_file = "parameters/example.prm";
    configuration.parameter_hash = "fnv1a64:test";
    configuration.plotting_profile_file = "parameters/plotting/example.json";
    configuration.plotting_profile_hash = "fnv1a64:plot";
    configuration.parameter_selection = "method=l-bfgs";
    configuration.declared_matrix = "method=steepest-descent,l-bfgs";
    configuration.excluded_combinations =
      "[method=steepest-descent,regularisation=1e-2]";
    configuration.resolved_combinations = "[method=l-bfgs]";
    configuration.comparison_rows = "method";
    configuration.comparison_columns = "regularisation";
    configuration.comparison_group_by = "none";
    const std::vector<std::string> expected{
      "artifacts/steepest-descent/beta-1e-1/artifact.kv",
      "artifacts/steepest-descent/beta-1e-2/artifact.kv"};
    const std::vector<std::string> command{
      "nmopt_runner", "--benchmark", "b1"};

    RunSetManifest manifest(configuration, command, expected);
    const auto running = read_file(manifest.path());
    require(running.find("\"status\": \"running\"") != std::string::npos,
            "new run manifests should start in running state");
    require(running.find("\"pending_count\": 2") != std::string::npos,
            "new run manifests should inventory pending artifacts");
    require(running.find("\"content_hash\": \"fnv1a64:test\"") !=
              std::string::npos,
            "run manifests should retain parameter provenance");
    require(running.find("\"resolved_comparison\": {\"rows\": \"method\"") !=
              std::string::npos,
            "run manifests should retain the comparison plan");
    require(
      running.find(
        "\"excluded_combinations\": \"[method=steepest-descent,regularisation=1e-2]\"") !=
        std::string::npos,
      "run manifests should retain excluded matrix coordinates");

    manifest.record_success(
      configuration.run_directory / expected.front());
    manifest.record_failure(
      configuration.run_directory / expected.back(), "solver failed");
    require(!manifest.finalize(),
            "a run manifest with a failed artifact should fail finalization");

    const auto finished = read_file(manifest.path());
    require(finished.find("\"status\": \"failed\"") != std::string::npos,
            "failed run manifests should retain failed state");
    require(finished.find("\"success_count\": 1") != std::string::npos,
            "run manifests should count successful artifacts");
    require(finished.find("\"failure_count\": 1") != std::string::npos,
            "run manifests should count failed artifacts");
    require(finished.find(
              "\"path\": \"artifacts/steepest-descent/beta-1e-1/artifact.kv\"") !=
              std::string::npos,
            "run manifests should retain artifact paths below artifacts");
    require(finished.find("\"error\": \"solver failed\"") !=
              std::string::npos,
            "run manifests should retain artifact failure diagnostics");
    std::filesystem::remove_all(temporary_root);
  }
} // namespace

int
main(const int argc, char **argv)
{
  try
    {
      const std::vector<nmopt::test_support::Scenario> scenarios{
        {"typed_capability_registries",
         "nmopt.runner.typed_capability_registries",
         {"backend-neutral", "application", "runner", "contract"},
         30,
         test_typed_capability_registries},
        {"reproduction_policy",
         "nmopt.runner.reproduction_policy",
         {"backend-neutral", "application", "runner", "contract"},
         30,
         test_reproduction_policy},
        {"run_kind_parser_rejects_unknown_values",
         "nmopt.runner.run_kind_parser_rejects_unknown_values",
         {"backend-neutral", "application", "runner", "contract", "negative"},
         30,
         test_run_kind_parser_rejects_unknown_values},
        {"matched_reference_stopping_criteria",
         "nmopt.runner.matched_reference_stopping_criteria",
         {"backend-neutral", "application", "runner", "contract"},
         30,
         test_matched_reference_stopping_criteria},
        {"run_manifest_records_state_and_failures",
         "nmopt.runner.run_manifest_records_state_and_failures",
         {"backend-neutral", "application", "runner", "contract"},
         30,
         test_run_manifest_records_state_and_failures}};
      const auto result = nmopt::test_support::run_requested_scenarios(
        argc, argv, scenarios, std::cout);
      if (!result.listed)
        std::cout << "runner contract tests passed: " << result.executed
                  << '\n';
      return 0;
    }
  catch (const std::exception &error)
    {
      std::cerr << "runner contract test failed: " << error.what() << '\n';
      return 1;
    }
}
