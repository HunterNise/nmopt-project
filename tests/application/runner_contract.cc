#include "../../apps/nmopt-runner/runner.hpp"
#include "nmopt/application/chapter6.hpp"
#include "../support/scenario_dispatch.hpp"

#include <filesystem>
#include <functional>
#include <iostream>
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
    auto options = parse({"nmopt_runner",
                          "--benchmark",
                          "b1",
                          "--framework-revision",
                          "source-revision"});
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
} // namespace

int
main(const int argc, char **argv)
{
  try
    {
      const std::vector<nmopt::test_support::Scenario> scenarios{
        {"reproduction_policy",
         "nmopt.runner.reproduction_policy",
         {"backend-neutral", "application", "runner", "contract"},
         30,
         test_reproduction_policy},
        {"run_kind_parser_rejects_unknown_values",
         "nmopt.runner.run_kind_parser_rejects_unknown_values",
         {"backend-neutral", "application", "runner", "contract", "negative"},
         30,
         test_run_kind_parser_rejects_unknown_values}};
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
