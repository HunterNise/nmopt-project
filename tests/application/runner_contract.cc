#include "../../apps/nmopt-runner/runner.hpp"
#include "../../apps/nmopt-runner/benchmark_registry.hpp"
#include "../../apps/nmopt-runner/capability_registry.hpp"
#include "../../apps/nmopt-runner/extension_contracts.hpp"
#include "nmopt/application/chapter6.hpp"
#include "../support/scenario_dispatch.hpp"

#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <utility>
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
    using nmopt::application::runner::PreconditionerSelection;
    using nmopt::application::runner::execution_capability_registry;
    using nmopt::application::runner::preconditioner_capability_registry;
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
    require(preconditioner_capability_registry().resolve("identity",
                                                          "preconditioner") ==
              PreconditionerSelection::identity_baseline,
            "preconditioner capability lookup should return the typed selection");

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
  test_future_extension_selection_contracts()
  {
    using nmopt::application::ScalarFunctionKind;
    using nmopt::application::chapter6::ExecutionSelection;
    using nmopt::application::chapter6::ProductSelection;
    using nmopt::application::runner::BoxBoundDataSelection;
    using nmopt::application::runner::CapabilityRegistry;
    using nmopt::application::runner::PreconditionerSelection;
    using nmopt::application::runner::resolve_quadratic_kkt_solver_selection;

    const BoxBoundDataSelection bounds{
      {"lower", ScalarFunctionKind::constant, -1.0, "", "test.lower"},
      {"upper", ScalarFunctionKind::expression, 0.0, "1.0 + x0", "test.upper"}};
    nmopt::application::runner::validate_box_bound_data_selection(bounds);

    const auto kkt = resolve_quadratic_kkt_solver_selection(
      "quadratic-kkt", "minres", "identity");
    require(kkt.product == ProductSelection::quadratic_kkt &&
              kkt.solver.method ==
                nmopt::contract::QuadraticKKTSolverMethod::minres &&
              kkt.preconditioner == PreconditionerSelection::identity_baseline,
            "future KKT selections should remain independently typed");

    const CapabilityRegistry<ExecutionSelection> assembled_only{
      {"assembled", ExecutionSelection::assembled}};
    require_invalid_argument(
      [&] { assembled_only.resolve("matrix-free", "execution"); },
      "an unregistered execution capability should fail resolution");

    require_invalid_argument(
      [] {
        (void)resolve_quadratic_kkt_solver_selection(
          "reduced-dto", "minres", "identity");
      },
      "a KKT extension should reject a non-KKT product");
  }

  nmopt::application::runner::RunSetPlan
  make_extension_fixture_plan(
    const std::string &benchmark_id,
    std::vector<nmopt::application::runner::ParameterAxis> axes)
  {
    using namespace nmopt::application::runner;
    RunSetPlan plan;
    plan.benchmark_id = benchmark_id;
    plan.matrix_axes = std::move(axes);

    ParameterCombination combination;
    for (const auto &axis : plan.matrix_axes)
      combination.values.emplace(axis.id, axis.values.front());
    RunSetCombination resolved;
    resolved.values = std::move(combination);
    resolved.artifact_coordinates = default_artifact_coordinate_components(
      plan.matrix_axes, resolved.values);
    plan.resolved_combinations.push_back(std::move(resolved));
    validate_run_set_plan(plan);
    return plan;
  }

  struct KktExtensionFixture
  {
    nmopt::application::runner::RunSetPlan plan;
    nmopt::application::runner::QuadraticKktSolverSelection solver;
    double reaction = 0.0;
  };

  KktExtensionFixture
  make_kkt_extension_fixture(const std::string &benchmark_id,
                             const double        reaction)
  {
    using nmopt::application::runner::resolve_quadratic_kkt_solver_selection;
    return {make_extension_fixture_plan(
              benchmark_id,
              {{"beta", {"1e-2"}},
               {"product", {"quadratic-kkt"}},
               {"kkt-method", {"minres"}},
               {"preconditioner", {"identity"}}}),
            resolve_quadratic_kkt_solver_selection(
              "quadratic-kkt", "minres", "identity"),
            reaction};
  }

  void
  test_future_b3_to_b6_registration_shape()
  {
    using nmopt::application::ScalarFunctionKind;
    using nmopt::application::runner::BoxBoundDataSelection;
    using nmopt::application::runner::find_benchmark_registration;
    using nmopt::application::runner::reduced_method_capability_registry;
    using nmopt::application::runner::validate_box_bound_data_selection;

    const auto b3 = make_extension_fixture_plan(
      "test-only.b3",
      {{"beta", {"1e-1", "1e-2"}},
       {"lower-bound", {"symmetric-constant"}},
       {"upper-bound", {"symmetric-constant"}},
       {"method", {"steepest-descent", "l-bfgs"}}});
    require(b3.matrix_axes.size() == 4 &&
              b3.matrix_axes[0].id == "beta" &&
              b3.matrix_axes[1].id == "lower-bound" &&
              b3.matrix_axes[2].id == "upper-bound" &&
              b3.matrix_axes[3].id == "method",
            "B3-shaped axes should remain generic run-set data");
    require(reduced_method_capability_registry().resolve(
              "l-bfgs", "solver method") ==
              nmopt::application::chapter6::ReducedMethod::limited_memory_bfgs,
            "B3-shaped method selections should use the typed registry");

    const BoxBoundDataSelection b4_bounds{
      {"lower-expression",
       ScalarFunctionKind::expression,
       0.0,
       "-0.5 - x0",
       "test.b4.lower"},
      {"upper-expression",
       ScalarFunctionKind::expression,
       0.0,
       "0.5 + x0",
       "test.b4.upper"}};
    validate_box_bound_data_selection(b4_bounds);
    const auto b4 = make_extension_fixture_plan(
      "test-only.b4",
      {{"lower-definition", {"lower-expression"}},
       {"upper-definition", {"upper-expression"}},
       {"beta", {"1e-2"}}});
    require(b4.matrix_axes[0].id == "lower-definition" &&
              b4.matrix_axes[1].id == "upper-definition",
            "B4-shaped lower and upper definitions should be independent axes");

    const auto b5 = make_kkt_extension_fixture("test-only.b5", 0.0);
    const auto b6 = make_kkt_extension_fixture("test-only.b6", 2.0);
    require(b5.solver.product == b6.solver.product &&
              b5.solver.solver.method == b6.solver.solver.method &&
              b5.solver.preconditioner == b6.solver.preconditioner &&
              b5.reaction != b6.reaction &&
              b5.plan.matrix_axes.size() == b6.plan.matrix_axes.size(),
            "B6 should reuse the B5 selection shape with changed reaction data");
    for (std::size_t index = 0; index < b5.plan.matrix_axes.size(); ++index)
      require(b5.plan.matrix_axes[index].id == b6.plan.matrix_axes[index].id &&
                b5.plan.matrix_axes[index].values ==
                  b6.plan.matrix_axes[index].values,
              "B6 should not need a distinct runner axis shape");

    for (const auto *id : {"b3", "b4", "b5", "b6"})
      require(find_benchmark_registration(id) == nullptr,
              "future benchmark fixtures must not be advertised as runnable");
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
        {"future_extension_selection_contracts",
         "nmopt.runner.future_extension_selection_contracts",
         {"backend-neutral", "application", "runner", "contract", "extension"},
         30,
         test_future_extension_selection_contracts},
        {"future_b3_to_b6_registration_shape",
         "nmopt.runner.future_b3_to_b6_registration_shape",
         {"backend-neutral", "application", "runner", "contract", "extension"},
         30,
         test_future_b3_to_b6_registration_shape},
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
