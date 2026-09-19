#include "nmopt/application/application.hpp"
#include "chapter6_execution.hpp"
#include "benchmark_registry.hpp"
#include "parameter_files.hpp"
#include "run_lifecycle.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
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
           << "       nmopt_runner --parameter-file FILE --framework-revision REV"
              " [--output DIRECTORY] [--run-slot SLOT] [--select AXIS=VALUE]"
              " [--refinement N]\n"
           << "       nmopt_runner --benchmark ID --framework-revision REV"
              " [--output DIRECTORY] [--run-kind KIND] [--run-slot SLOT]"
              " [--refinement N]\n"
           << "       nmopt_runner --help\n"
           << "\n"
           << "--list             list registered Chapter 5/6 application entries\n"
           << "--benchmark ID     run a registered benchmark by identifier\n"
           << "--parameter-file FILE\n"
           << "                   load a versioned experiment family\n"
           << "--select AXIS=VALUE\n"
           << "                   filter one declared matrix axis (repeatable)\n"
           << "--output DIRECTORY set the generated run-set root (default: runs)\n"
           << "--run-kind KIND    use reproduction or development policy\n"
           << "--run-slot SLOT    use a named development slot instead of auto-allocation\n"
           << "--framework-revision REV\n"
           << "                   record the framework revision in each artifact\n"
           << "--refinement N     optional framework-native mesh override; otherwise\n"
           << "                   use the selected benchmark configuration\n"
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

  std::string
  join_matrix(const nmopt::application::runner::ParameterFile &file)
  {
    std::ostringstream output;
    for (std::size_t index = 0; index < file.matrix.size(); ++index)
      {
        if (index != 0)
          output << ';';
        output << file.matrix[index].id << '=';
        for (std::size_t value = 0; value < file.matrix[index].values.size(); ++value)
          {
            if (value != 0)
              output << ',';
            output << file.matrix[index].values[value];
          }
      }
    return output.str();
  }

  std::string
  join_combinations(const std::vector<nmopt::application::runner::ParameterCombination> &combinations)
  {
    std::ostringstream output;
    for (std::size_t index = 0; index < combinations.size(); ++index)
      {
        if (index != 0)
          output << ';';
        output << '[';
        std::size_t value_index = 0;
        for (const auto &[axis, value] : combinations[index].values)
          {
            if (value_index++ != 0)
              output << ',';
            output << axis << '=' << value;
          }
        output << ']';
      }
    return output.str();
  }

  std::string
  join_combinations(const nmopt::application::runner::RunSetPlan &plan)
  {
    std::ostringstream output;
    for (std::size_t index = 0; index < plan.resolved_combinations.size();
         ++index)
      {
        if (index != 0)
          output << ';';
        output << '[';
        std::size_t value_index = 0;
        for (const auto &[axis, value] :
             plan.resolved_combinations[index].values.values)
          {
            if (value_index++ != 0)
              output << ',';
            output << axis << '=' << value;
          }
        output << ']';
      }
    return output.str();
  }

  std::string
  join_selection(const nmopt::application::runner::ParameterFile &file,
                 const std::vector<std::pair<std::string, std::string>> &cli)
  {
    std::ostringstream output;
    bool                first = true;
    for (const auto &[axis, value] : file.selection)
      {
        if (!first)
          output << ',';
        first = false;
        output << axis << '=' << value;
      }
    for (const auto &[axis, value] : cli)
      {
        if (!first)
          output << ',';
        first = false;
        output << axis << '=' << value;
      }
    return output.str();
  }

  std::vector<std::string>
  command_line_arguments(const int argc, char **argv)
  {
    std::vector<std::string> command;
    command.reserve(static_cast<std::size_t>(argc));
    for (int index = 0; index < argc; ++index)
      command.emplace_back(argv[index]);
    return command;
  }

  nmopt::application::runner::ParameterFile
  load_parameter_file(const nmopt::application::runner::CommandLineOptions &options)
  {
    std::filesystem::path path;
    if (options.parameter_file.has_value())
      path = nmopt::application::runner::find_file_from_current_or_parent(
        *options.parameter_file);
    else
      {
        const auto *registration =
          nmopt::application::runner::find_benchmark_registration(
            options.benchmark.value_or(""));
        if (registration == nullptr)
          throw std::invalid_argument(
            "unsupported benchmark '" + options.benchmark.value_or("") +
            "'; available benchmark IDs: " +
            nmopt::application::runner::registered_benchmark_ids());
        path = nmopt::application::runner::find_file_from_current_or_parent(
          registration->default_parameter_file);
      }
    return nmopt::application::runner::read_parameter_file(path);
  }

  void
  copy_configuration_file(const std::filesystem::path &source,
                          const std::filesystem::path &destination)
  {
    nmopt::application::runner::prepare_artifact_path(destination);
    std::ifstream input(source);
    std::ofstream output(destination);
    if (!input || !output)
      throw std::runtime_error("could not snapshot configuration file '" +
                               source.string() + "'");
    output << input.rdbuf();
    if (!output)
      throw std::runtime_error("could not write configuration snapshot '" +
                               destination.string() + "'");
  }

  void
  snapshot_configuration(
    const nmopt::application::runner::ResolvedRunConfiguration &configuration,
    const nmopt::application::runner::ParameterFile &file,
    const nmopt::application::runner::RunSetPlan &plan)
  {
    std::filesystem::create_directories(configuration.run_directory);
    copy_configuration_file(file.path,
                            configuration.run_directory / "parameters.prm");
    copy_configuration_file(configuration.plotting_profile_file,
                            configuration.run_directory / "plotting-profile.json");
    std::ofstream resolved(configuration.run_directory / "resolved-combinations.txt");
    if (!resolved)
      throw std::runtime_error("could not write resolved parameter combinations");
    resolved << join_combinations(plan) << '\n';
  }

  struct PreparedRun
  {
    nmopt::application::runner::CommandLineOptions options;
    nmopt::application::runner::ParameterFile       file;
    nmopt::application::runner::ResolvedRunConfiguration configuration;
    nmopt::application::runner::RunSetPlan              plan;
  };

  PreparedRun
  prepare_run(const nmopt::application::runner::CommandLineOptions &input_options)
  {
    using namespace nmopt::application::runner;
    PreparedRun prepared{
      input_options,
      load_parameter_file(input_options),
      ResolvedRunConfiguration({},
                               {},
                               "",
                               "",
                               "",
                               RunKind::reproduction,
                               std::nullopt),
      {}};
    const auto benchmark_id = prepared.file.value("Benchmark/id");
    const auto *registration =
      nmopt::application::runner::find_benchmark_registration_for_parameter_id(
        benchmark_id);
    if (registration == nullptr)
      throw std::invalid_argument("parameter file declares unsupported benchmark '" +
                                  benchmark_id + "'; available benchmark IDs: " +
                                  nmopt::application::runner::registered_benchmark_ids());
    if (input_options.parameter_file.has_value())
      {
        prepared.options.benchmark = std::string(registration->id);
        prepared.options.run_kind =
          parse_run_kind(prepared.file.value("Run/kind"));
        if (!input_options.output_directory_explicit)
          prepared.options.output_directory = prepared.file.value("Run/output root");
        if (prepared.options.run_kind == RunKind::reproduction &&
            prepared.file.value("Run/build profile") != NMOPT_COMPILED_BUILD_PROFILE)
          throw std::invalid_argument(
            "parameter Run/build profile does not match the compiled runner profile");
      }
    else if (!input_options.benchmark.has_value() ||
             *input_options.benchmark != registration->id)
      throw std::invalid_argument("--benchmark does not match its authoritative parameter file");

    prepared.plan = make_run_set_plan(prepared.file,
                                      prepared.options.selection_filters);
    const auto *execution_registration =
      nmopt::application::runner::chapter6_execution::
        find_benchmark_execution_registration(registration->id);
    if (execution_registration == nullptr)
      throw std::invalid_argument(
        "benchmark registration has no execution callback: " +
        std::string(registration->id));
    // This also validates the required matrix axes before creating output.
    (void)execution_registration->artifact_planner(prepared.plan);

    prepared.configuration = resolve_run_configuration(
      prepared.options, NMOPT_COMPILED_BUILD_PROFILE);
    prepared.configuration.parameter_file =
      prepared.plan.parameter_provenance.file;
    prepared.configuration.parameter_hash =
      prepared.plan.parameter_provenance.content_hash;
    const auto profile_path = find_file_from_current_or_parent(
      prepared.file.value("Postprocessing/style profile"));
    prepared.configuration.plotting_profile_file = profile_path;
    prepared.configuration.plotting_profile_hash = parameter_file_hash(profile_path);
    prepared.configuration.parameter_selection = join_selection(
      prepared.file, prepared.options.selection_filters);
    prepared.configuration.declared_matrix = join_matrix(prepared.file);
    prepared.configuration.excluded_combinations =
      join_combinations(prepared.plan.excluded_combinations);
    prepared.configuration.resolved_combinations = join_combinations(prepared.plan);
    prepared.configuration.comparison_rows = prepared.plan.comparison.rows;
    prepared.configuration.comparison_columns = prepared.plan.comparison.columns;
    prepared.configuration.comparison_group_by = prepared.plan.comparison.group_by;
    return prepared;
  }
} // namespace

int
main(const int argc, char **argv)
{
  try
    {
      const auto command = command_line_arguments(argc, argv);
      const auto options =
        nmopt::application::runner::parse_command_line(argc, argv);
      if (options.help)
        {
          print_usage(std::cout);
          return 0;
        }

      if (options.benchmark.has_value() || options.parameter_file.has_value())
        {
          auto prepared = prepare_run(options);
          snapshot_configuration(prepared.configuration,
                                 prepared.file,
                                 prepared.plan);
          const auto *registration =
            nmopt::application::runner::chapter6_execution::
              find_benchmark_execution_registration(
                prepared.configuration.benchmark);
          if (registration == nullptr)
            throw std::invalid_argument(
              "benchmark has no execution registration: " +
              prepared.configuration.benchmark);
          nmopt::application::runner::RunSetManifest run_manifest(
            prepared.configuration,
            command,
            registration->artifact_planner(prepared.plan));
          const bool execution_succeeded = registration->execute(
            prepared.configuration,
            command,
            prepared.plan,
            prepared.file,
            run_manifest);
          return run_manifest.finalize() && execution_succeeded ? 0 : 1;
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
