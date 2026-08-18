#pragma once

#include <charconv>
#include <filesystem>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace nmopt::application::runner
{
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
    RunKind               run_kind = RunKind::reproduction;
    std::filesystem::path output_directory = "runs";
    std::string            framework_revision;
    std::optional<unsigned int> refinement_override;
  };

  struct ResolvedRunConfiguration
  {
    std::filesystem::path        output_root;
    std::filesystem::path        run_directory;
    std::string                  benchmark;
    std::string                  build_profile;
    std::string                  framework_revision;
    RunKind                      run_kind = RunKind::reproduction;
    std::optional<unsigned int>  refinement_override;
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
        if (argument == "--output")
          {
            if (index + 1 >= argc)
              throw std::invalid_argument(
                "--output needs a directory argument");
            options.output_directory = argv[++index];
            if (options.output_directory.empty())
              throw std::invalid_argument(
                "--output needs a nonempty directory argument");
            continue;
          }
        if (argument == "--run-kind")
          {
            if (index + 1 >= argc)
              throw std::invalid_argument(
                "--run-kind needs a kind argument");
            options.run_kind = parse_run_kind(argv[++index]);
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

    if (options.list && (options.run_b1 || options.run_b2))
      throw std::invalid_argument(
        "--list and --benchmark cannot be selected together");
    if (options.run_b1 && options.run_b2)
      throw std::invalid_argument(
        "select only one benchmark per runner invocation");
    if (!options.help && !options.list && !options.run_b1 && !options.run_b2)
      throw std::invalid_argument(
        "select --list, --benchmark b1/b2, or --help");
    if ((options.run_b1 || options.run_b2) &&
        options.framework_revision.empty())
      throw std::invalid_argument(
        "benchmark runs require --framework-revision for provenance");
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
    const auto run_directory =
      options.run_kind == RunKind::reproduction
        ? benchmark_directory / "authoritative"
        : next_development_directory(benchmark_directory);

    return {options.output_directory,
            run_directory,
            benchmark,
            build_profile,
            framework_revision,
            options.run_kind,
            options.refinement_override};
  }

  inline std::filesystem::path
  artifact_path(const std::filesystem::path &output_directory,
                const std::vector<std::string> &relative_components)
  {
    if (relative_components.empty())
      throw std::invalid_argument(
        "an artifact path needs at least one relative component");

    std::filesystem::path path = output_directory;
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
} // namespace nmopt::application::runner
