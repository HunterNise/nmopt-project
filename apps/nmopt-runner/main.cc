#include "nmopt/application/application.hpp"
#include "runner.hpp"

#include <iostream>
#include <stdexcept>
#include <string>

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
           << "       nmopt_runner --help\n"
           << "\n"
           << "--list             list registered Chapter 5/6 application entries\n"
           << "--output DIRECTORY set the runner-owned artifact root (default: runs)\n"
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
