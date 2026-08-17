#pragma once

#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace nmopt::application::runner
{
  struct CommandLineOptions
  {
    bool                 list = false;
    bool                 help = false;
    std::filesystem::path output_directory = "runs";
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
        throw std::invalid_argument("unknown runner argument '" +
                                    std::string(argument) + "'");
      }

    if (!options.help && !options.list)
      throw std::invalid_argument(
        "select --list or --help; benchmark execution is not registered yet");
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
