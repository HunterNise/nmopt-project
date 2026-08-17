#pragma once

#include <charconv>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace nmopt::application::runner
{
  struct CommandLineOptions
  {
    bool                  list = false;
    bool                  run_b1 = false;
    bool                  help = false;
    std::filesystem::path output_directory = "runs";
    std::string            framework_revision;
    unsigned int           refinement = 7;
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
            if (benchmark != "b1")
              throw std::invalid_argument(
                "only benchmark 'b1' is registered in this runner unit");
            options.run_b1 = true;
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
            options.refinement = refinement;
            continue;
          }
        throw std::invalid_argument("unknown runner argument '" +
                                    std::string(argument) + "'");
      }

    if (options.list && options.run_b1)
      throw std::invalid_argument(
        "--list and --benchmark cannot be selected together");
    if (!options.help && !options.list && !options.run_b1)
      throw std::invalid_argument(
        "select --list, --benchmark b1, or --help");
    if (options.run_b1 && options.framework_revision.empty())
      throw std::invalid_argument(
        "B1 runs require --framework-revision for provenance");
    return options;
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
