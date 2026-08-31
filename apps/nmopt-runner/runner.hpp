#pragma once

#include "nmopt/solvers/reduced_search.hpp"

#include <algorithm>
#include <charconv>
#include <filesystem>
#include <fstream>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace nmopt::application::runner
{
  inline bool
  reference_reached_stopping_tolerance(
    const solvers::ReducedStoppingCriterion criterion,
    const solvers::ReducedStoppingReason    reason)
  {
    using Criterion = solvers::ReducedStoppingCriterion;
    using Reason = solvers::ReducedStoppingReason;

    switch (criterion)
      {
        case Criterion::automatic:
          return reason == Reason::gradient_tolerance ||
                 reason == Reason::relative_gradient_tolerance ||
                 reason == Reason::objective_change_tolerance ||
                 reason == Reason::step_tolerance;
        case Criterion::gradient_norm:
          return reason == Reason::gradient_tolerance;
        case Criterion::relative_gradient_norm:
          return reason == Reason::relative_gradient_tolerance;
        case Criterion::objective_change:
          return reason == Reason::objective_change_tolerance ||
                 reason == Reason::stationary;
        case Criterion::step_norm:
          return reason == Reason::step_tolerance ||
                 reason == Reason::stationary;
      }
    return false;
  }

  enum class RunKind
  {
    reproduction,
    development
  };

  inline std::string_view
  run_kind_name(const RunKind kind)
  {
    switch (kind)
      {
        case RunKind::reproduction:
          return "reproduction";
        case RunKind::development:
          return "development";
      }
    return "unknown";
  }

  inline RunKind
  parse_run_kind(const std::string_view value)
  {
    if (value == "reproduction")
      return RunKind::reproduction;
    if (value == "development")
      return RunKind::development;
    throw std::invalid_argument(
      "--run-kind needs 'reproduction' or 'development'");
  }

  struct CommandLineOptions
  {
    bool                  list = false;
    bool                  run_b1 = false;
    bool                  run_b2 = false;
    bool                  help = false;
    bool                  output_directory_explicit = false;
    bool                  run_kind_explicit = false;
    RunKind               run_kind = RunKind::reproduction;
    std::filesystem::path output_directory = "runs";
    std::string            framework_revision;
    std::optional<std::filesystem::path> parameter_file;
    std::optional<std::string> run_slot;
    std::vector<std::pair<std::string, std::string>> selection_filters;
    std::optional<unsigned int> refinement_override;
  };

  struct ResolvedRunConfiguration
  {
    ResolvedRunConfiguration(
      const std::filesystem::path       output_root_value,
      const std::filesystem::path       run_directory_value,
      const std::string                 benchmark_value,
      const std::string                 build_profile_value,
      const std::string                 framework_revision_value,
      const RunKind                      run_kind_value,
      const std::optional<unsigned int> refinement_override_value)
      : output_root(output_root_value)
      , run_directory(run_directory_value)
      , benchmark(benchmark_value)
      , build_profile(build_profile_value)
      , framework_revision(framework_revision_value)
      , run_kind(run_kind_value)
      , refinement_override(refinement_override_value)
    {}

    std::filesystem::path        output_root;
    std::filesystem::path        run_directory;
    std::string                  benchmark;
    std::string                  build_profile;
    std::string                  framework_revision;
    RunKind                      run_kind = RunKind::reproduction;
    std::optional<unsigned int>  refinement_override;
    std::filesystem::path        parameter_file;
    std::string                  parameter_hash;
    std::filesystem::path        plotting_profile_file;
    std::string                  plotting_profile_hash;
    std::string                  parameter_selection;
    std::string                  declared_matrix;
    std::string                  excluded_combinations;
    std::string                  resolved_combinations;
    std::string                  comparison_rows;
    std::string                  comparison_columns;
    std::string                  comparison_group_by;
  };

  inline CommandLineOptions
  parse_command_line(const int argc, char **argv)
  {
    CommandLineOptions options;
    for (int index = 1; index < argc; ++index)
      {
        const std::string_view argument(argv[index]);
        if (argument == "--help" || argument == "-h")
          {
            options.help = true;
            continue;
          }
        if (argument == "--list")
          {
            options.list = true;
            continue;
          }
        if (argument == "--benchmark")
          {
            if (index + 1 >= argc)
              throw std::invalid_argument(
                "--benchmark needs a benchmark identifier");
            const std::string_view benchmark(argv[++index]);
            if (benchmark != "b1" && benchmark != "b2")
              throw std::invalid_argument(
                "only benchmarks 'b1' and 'b2' are registered in this runner unit");
            options.run_b1 = benchmark == "b1";
            options.run_b2 = benchmark == "b2";
            continue;
          }
        if (argument == "--parameter-file")
          {
            if (index + 1 >= argc)
              throw std::invalid_argument(
                "--parameter-file needs a file argument");
            const std::filesystem::path path(argv[++index]);
            if (path.empty())
              throw std::invalid_argument(
                "--parameter-file needs a nonempty file argument");
            options.parameter_file = path;
            continue;
          }
        if (argument == "--select")
          {
            if (index + 1 >= argc)
              throw std::invalid_argument(
                "--select needs an axis=value argument");
            const std::string expression(argv[++index]);
            const auto        separator = expression.find('=');
            if (separator == std::string::npos || separator == 0 ||
                separator + 1 == expression.size())
              throw std::invalid_argument(
                "--select needs an axis=value argument");
            options.selection_filters.emplace_back(
              expression.substr(0, separator), expression.substr(separator + 1));
            continue;
          }
        if (argument == "--output")
          {
            if (index + 1 >= argc)
              throw std::invalid_argument(
                "--output needs a directory argument");
            options.output_directory = argv[++index];
            if (options.output_directory.empty())
              throw std::invalid_argument(
                "--output needs a nonempty directory argument");
            options.output_directory_explicit = true;
            continue;
          }
        if (argument == "--run-slot")
          {
            if (index + 1 >= argc)
              throw std::invalid_argument(
                "--run-slot needs a development slot name");
            const std::string_view slot(argv[++index]);
            if (slot.empty() || slot == "." || slot == ".." ||
                slot.find_first_of("/\\") != std::string_view::npos)
              throw std::invalid_argument(
                "--run-slot needs a nonempty single development slot name");
            options.run_slot = std::string(slot);
            continue;
          }
        if (argument == "--run-kind")
          {
            if (index + 1 >= argc)
              throw std::invalid_argument(
                "--run-kind needs a kind argument");
            options.run_kind = parse_run_kind(argv[++index]);
            options.run_kind_explicit = true;
            continue;
          }
        if (argument == "--framework-revision")
          {
            if (index + 1 >= argc)
              throw std::invalid_argument(
                "--framework-revision needs a revision argument");
            options.framework_revision = argv[++index];
            if (options.framework_revision.empty())
              throw std::invalid_argument(
                "--framework-revision needs a nonempty revision argument");
            continue;
          }
        if (argument == "--refinement")
          {
            if (index + 1 >= argc)
              throw std::invalid_argument(
                "--refinement needs a nonnegative integer argument");
            const std::string_view value(argv[++index]);
            const auto              first = value.data();
            const auto              last = first + value.size();
            unsigned int             refinement = 0;
            const auto conversion =
              std::from_chars(first, last, refinement);
            if (conversion.ec != std::errc{} || conversion.ptr != last)
              throw std::invalid_argument(
                "--refinement needs a nonnegative integer argument");
            options.refinement_override = refinement;
            continue;
          }
        throw std::invalid_argument("unknown runner argument '" +
                                    std::string(argument) + "'");
      }

    if (options.list &&
        (options.run_b1 || options.run_b2 || options.parameter_file.has_value()))
      throw std::invalid_argument(
        "--list and a run selection cannot be selected together");
    if (options.parameter_file.has_value() && (options.run_b1 || options.run_b2))
      throw std::invalid_argument(
        "use either --parameter-file or --benchmark, not both");
    if (options.run_b1 && options.run_b2)
      throw std::invalid_argument(
        "select only one benchmark per runner invocation");
    if (!options.help && !options.list && !options.run_b1 && !options.run_b2 &&
        !options.parameter_file.has_value())
      throw std::invalid_argument(
        "select --list, --parameter-file, --benchmark b1/b2, or --help");
    if ((options.run_b1 || options.run_b2 || options.parameter_file.has_value()) &&
        options.framework_revision.empty())
      throw std::invalid_argument(
        "benchmark runs require --framework-revision for provenance");
    if (options.parameter_file.has_value() && options.run_kind_explicit)
      throw std::invalid_argument(
        "the parameter file owns run kind; use a file-specific development family");
    if (options.run_slot.has_value() && !options.parameter_file.has_value() &&
        options.run_kind != RunKind::development)
      throw std::invalid_argument(
        "--run-slot is only valid for development runs");
    return options;
  }

  inline void
  validate_run_policy(const CommandLineOptions &options,
                      const std::string_view  compiled_build_profile)
  {
    if ((!options.run_b1 && !options.run_b2) ||
        options.run_kind == RunKind::development)
      return;

    if (compiled_build_profile != "release-dealii")
      throw std::invalid_argument(
        "reproduction runs require a release-dealii build; compiled profile is '" +
        std::string(compiled_build_profile) +
        "'. Use --run-kind development for a development run");
  }

  inline std::string
  benchmark_name(const CommandLineOptions &options)
  {
    if (options.run_b1)
      return "b1";
    if (options.run_b2)
      return "b2";
    throw std::invalid_argument(
      "a run configuration needs a selected benchmark");
  }

  inline std::string
  development_run_id(const unsigned int number)
  {
    auto value = std::to_string(number);
    if (value.size() < 3)
      value.insert(0, 3 - value.size(), '0');
    return value;
  }

  inline std::filesystem::path
  next_development_directory(const std::filesystem::path &benchmark_directory)
  {
    const auto development_directory = benchmark_directory / "development";
    for (unsigned int number = 1;; ++number)
      {
        const auto candidate =
          development_directory / development_run_id(number);
        std::error_code error;
        const bool exists = std::filesystem::exists(candidate, error);
        if (error)
          throw std::filesystem::filesystem_error(
            "could not inspect development run directory", candidate, error);
        if (!exists)
          return candidate;
      }
  }

  inline ResolvedRunConfiguration
  resolve_run_configuration(const CommandLineOptions &options,
                            const std::string_view  compiled_build_profile)
  {
    if (!options.run_b1 && !options.run_b2)
      throw std::invalid_argument(
        "a run configuration needs a selected benchmark");
    validate_run_policy(options, compiled_build_profile);

    const auto benchmark = benchmark_name(options);
    const auto build_profile = std::string(compiled_build_profile);
    const auto framework_revision = options.framework_revision;
    const auto benchmark_directory =
      options.output_directory / "chapter-6" / benchmark;
    const auto development_directory = benchmark_directory / "development";
    const auto run_directory =
      options.run_kind == RunKind::reproduction
        ? benchmark_directory / "authoritative"
        : options.run_slot.has_value()
            ? development_directory / *options.run_slot
            : next_development_directory(benchmark_directory);

    if (options.run_kind == RunKind::development && options.run_slot.has_value())
      {
        std::error_code error;
        if (std::filesystem::exists(run_directory, error))
          throw std::invalid_argument(
            "the requested development run slot already exists: " +
            run_directory.string());
        if (error)
          throw std::filesystem::filesystem_error(
            "could not inspect requested development run slot",
            run_directory,
            error);
      }

    return ResolvedRunConfiguration(options.output_directory,
                                    run_directory,
                                    benchmark,
                                    build_profile,
                                    framework_revision,
                                    options.run_kind,
                                    options.refinement_override);
  }

  inline std::filesystem::path
  artifact_path(const std::filesystem::path &output_directory,
                const std::vector<std::string> &relative_components)
  {
    if (relative_components.empty())
      throw std::invalid_argument(
        "an artifact path needs at least one relative component");

    std::filesystem::path path = output_directory / "artifacts";
    for (const auto &component : relative_components)
      {
        if (component.empty() || component == "." || component == ".." ||
            component.find_first_of("/\\") != std::string::npos)
          throw std::invalid_argument(
            "artifact path components must be nonempty single names");
        path /= component;
      }
    return path / "artifact.kv";
  }

  inline void
  prepare_artifact_path(const std::filesystem::path &path)
  {
    const auto parent = path.parent_path();
    if (parent.empty())
      throw std::invalid_argument("an artifact path needs a parent directory");
    std::filesystem::create_directories(parent);
  }

  inline std::string
  json_escape(const std::string_view value)
  {
    const char hex[] = "0123456789abcdef";
    std::string result;
    result.reserve(value.size());
    for (const unsigned char character : value)
      switch (character)
        {
          case '"':
            result += "\\\"";
            break;
          case '\\':
            result += "\\\\";
            break;
          case '\b':
            result += "\\b";
            break;
          case '\f':
            result += "\\f";
            break;
          case '\n':
            result += "\\n";
            break;
          case '\r':
            result += "\\r";
            break;
          case '\t':
            result += "\\t";
            break;
          default:
            if (character < 0x20)
              {
                result += "\\u00";
                result += hex[(character >> 4) & 0x0f];
                result += hex[character & 0x0f];
              }
            else
              result += static_cast<char>(character);
            break;
        }
    return result;
  }

  inline std::string
  json_string(const std::string_view value)
  {
    return '"' + json_escape(value) + '"';
  }

  struct RunArtifactRecord
  {
    std::string relative_path;
    std::string status = "pending";
    std::string error;
  };

  class RunSetManifest final
  {
  public:
    RunSetManifest(const ResolvedRunConfiguration &configuration,
                   const std::vector<std::string> &command,
                   const std::vector<std::string> &expected_artifacts)
      : run_directory_(configuration.run_directory)
      , manifest_path_(run_directory_ / "run-manifest.json")
      , benchmark_(configuration.benchmark)
      , build_profile_(configuration.build_profile)
      , framework_revision_(configuration.framework_revision)
      , run_kind_(run_kind_name(configuration.run_kind))
      , refinement_override_(configuration.refinement_override)
      , parameter_file_(configuration.parameter_file)
      , parameter_hash_(configuration.parameter_hash)
      , plotting_profile_file_(configuration.plotting_profile_file)
      , plotting_profile_hash_(configuration.plotting_profile_hash)
      , parameter_selection_(configuration.parameter_selection)
      , declared_matrix_(configuration.declared_matrix)
      , excluded_combinations_(configuration.excluded_combinations)
      , resolved_combinations_(configuration.resolved_combinations)
      , comparison_rows_(configuration.comparison_rows)
      , comparison_columns_(configuration.comparison_columns)
      , comparison_group_by_(configuration.comparison_group_by)
      , command_(command)
    {
      std::filesystem::create_directories(run_directory_);
      records_.reserve(expected_artifacts.size());
      for (const auto &relative_path : expected_artifacts)
        {
          if (relative_path.empty())
            throw std::invalid_argument(
              "run manifest artifact paths need nonempty values");
          const auto duplicate = std::find_if(
            records_.begin(), records_.end(), [&](const auto &record) {
              return record.relative_path == relative_path;
            });
          if (duplicate != records_.end())
            throw std::invalid_argument(
              "run manifest artifact paths must be unique");
          records_.push_back({relative_path, "pending", {}});
        }
      write();
    }

    void
    record_success(const std::filesystem::path &artifact_path)
    {
      update(artifact_path, "ok", {});
      write();
    }

    void
    record_failure(const std::filesystem::path &artifact_path,
                   const std::string_view    error)
    {
      update(artifact_path, "error", error);
      write();
    }

    bool
    finalize()
    {
      finalized_ = true;
      for (auto &record : records_)
        if (record.status == "pending")
          {
            record.status = "error";
            record.error = "artifact was not executed";
          }
      write();
      return failure_count() == 0;
    }

    const std::filesystem::path &
    path() const
    {
      return manifest_path_;
    }

  private:
    void
    update(const std::filesystem::path &artifact_path,
           const std::string_view    status,
           const std::string_view    error)
    {
      const auto relative_path =
        artifact_path.lexically_relative(run_directory_).generic_string();
      const auto record = std::find_if(
        records_.begin(), records_.end(), [&](const auto &candidate) {
          return candidate.relative_path == relative_path;
        });
      if (record == records_.end())
        throw std::invalid_argument(
          "artifact path is not part of the run manifest inventory");
      if (record->status != "pending")
        throw std::invalid_argument(
          "run manifest cannot record an artifact twice");
      record->status = status;
      record->error = error;
    }

    std::size_t
    success_count() const
    {
      return static_cast<std::size_t>(std::count_if(
        records_.begin(), records_.end(), [](const auto &record) {
          return record.status == "ok";
        }));
    }

    std::size_t
    failure_count() const
    {
      return static_cast<std::size_t>(std::count_if(
        records_.begin(), records_.end(), [](const auto &record) {
          return record.status == "error";
        }));
    }

    std::size_t
    pending_count() const
    {
      return records_.size() - success_count() - failure_count();
    }

    void
    write() const
    {
      const auto status =
        !finalized_ ? "running" : (failure_count() == 0 ? "complete" : "failed");
      std::ofstream output(manifest_path_);
      if (!output)
        throw std::runtime_error("could not open run manifest '" +
                                 manifest_path_.string() + "'");

      output << "{\n"
             << "  \"schema\": \"nmopt-run-set-v1\",\n"
             << "  \"status\": " << json_string(status) << ",\n"
             << "  \"benchmark\": " << json_string(benchmark_) << ",\n"
             << "  \"run_kind\": " << json_string(run_kind_) << ",\n"
             << "  \"build_profile\": " << json_string(build_profile_)
             << ",\n"
             << "  \"framework_revision\": "
             << json_string(framework_revision_) << ",\n"
             << "  \"run_directory\": "
             << json_string(run_directory_.string()) << ",\n"
             << "  \"refinement_override\": ";
      if (refinement_override_.has_value())
        output << *refinement_override_;
      else
        output << "null";
      output << ",\n"
             << "  \"parameters\": {\"file\": "
             << json_string(parameter_file_.string())
             << ", \"content_hash\": " << json_string(parameter_hash_)
             << ", \"selection\": "
             << json_string(parameter_selection_)
             << ", \"declared_matrix\": "
             << json_string(declared_matrix_)
             << ", \"excluded_combinations\": "
             << json_string(excluded_combinations_)
             << ", \"resolved_combinations\": "
             << json_string(resolved_combinations_) << "},\n"
             << "  \"plotting\": {\"profile_file\": "
             << json_string(plotting_profile_file_.string())
             << ", \"content_hash\": "
             << json_string(plotting_profile_hash_)
             << ", \"resolved_comparison\": {\"rows\": "
             << json_string(comparison_rows_)
             << ", \"columns\": " << json_string(comparison_columns_)
             << ", \"group_by\": "
             << json_string(comparison_group_by_) << "}},\n"
             << "  \"command\": [";
      for (std::size_t index = 0; index < command_.size(); ++index)
        {
          if (index != 0)
            output << ", ";
          output << json_string(command_[index]);
        }
      output << "],\n"
             << "  \"expected_artifact_count\": " << records_.size()
             << ",\n"
             << "  \"success_count\": " << success_count() << ",\n"
             << "  \"failure_count\": " << failure_count() << ",\n"
             << "  \"pending_count\": " << pending_count() << ",\n"
             << "  \"artifacts\": [\n";
      for (std::size_t index = 0; index < records_.size(); ++index)
        {
          const auto &record = records_[index];
          output << "    {\"path\": "
                 << json_string(record.relative_path)
                 << ", \"status\": " << json_string(record.status);
          if (!record.error.empty())
            output << ", \"error\": " << json_string(record.error);
          output << "}" << (index + 1 == records_.size() ? "\n" : ",\n");
        }
      output << "  ]\n}\n";
      if (!output)
        throw std::runtime_error("could not write run manifest '" +
                                 manifest_path_.string() + "'");
    }

    std::filesystem::path           run_directory_;
    std::filesystem::path           manifest_path_;
    std::string                     benchmark_;
    std::string                     build_profile_;
    std::string                     framework_revision_;
    std::string                     run_kind_;
    std::optional<unsigned int>     refinement_override_;
    std::filesystem::path           parameter_file_;
    std::string                     parameter_hash_;
    std::filesystem::path           plotting_profile_file_;
    std::string                     plotting_profile_hash_;
    std::string                     parameter_selection_;
    std::string                     declared_matrix_;
    std::string                     excluded_combinations_;
    std::string                     resolved_combinations_;
    std::string                     comparison_rows_;
    std::string                     comparison_columns_;
    std::string                     comparison_group_by_;
    std::vector<std::string>        command_;
    std::vector<RunArtifactRecord>  records_;
    bool                            finalized_ = false;
  };
} // namespace nmopt::application::runner
