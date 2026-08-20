#include "../../apps/nmopt-runner/parameter_files.hpp"
#include "../support/scenario_dispatch.hpp"

#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
  using nmopt::application::runner::find_file_from_current_or_parent;
  using nmopt::application::runner::read_parameter_file;

  void
  require(const bool condition, const char *message)
  {
    if (!condition)
      throw std::runtime_error(message);
  }

  void
  test_checked_in_families_expand_and_filter()
  {
    const auto b1 = read_parameter_file(find_file_from_current_or_parent(
      "parameters/chapter-6/b1/authoritative.prm"));
    require(b1.matrix.size() == 2, "B1 should declare two matrix axes");
    require(b1.combinations().size() == 8,
            "B1 should expand to eight authoritative combinations");
    require(b1.combinations({{"method", "l-bfgs"},
                             {"regularisation", "1e-6"}})
              .size() == 1,
            "B1 selection filters should resolve one combination");

    const auto b2 = read_parameter_file(find_file_from_current_or_parent(
      "parameters/chapter-6/b2/authoritative.prm"));
    require(b2.matrix.size() == 2, "B2 should declare independent axes");
    require(b2.combinations().size() == 4,
            "B2 should expand to four authoritative combinations");

    const auto development = read_parameter_file(find_file_from_current_or_parent(
      "parameters/chapter-6/b2/development/forcing-sweep.prm"));
    require(development.combinations().size() == 3,
            "B2 forcing development family should expand to three combinations");
    require(development.content_hash.rfind("fnv1a64:", 0) == 0,
            "parameter provenance should carry a labelled deterministic hash");
  }

  void
  test_unknown_selection_is_rejected()
  {
    const auto b1 = read_parameter_file(find_file_from_current_or_parent(
      "parameters/chapter-6/b1/authoritative.prm"));
    try
      {
        (void)b1.combinations({{"regularisation", "1e-4"}});
      }
    catch (const std::invalid_argument &)
      {
        return;
      }
    throw std::runtime_error("unknown matrix selection values should be rejected");
  }
} // namespace

int
main(const int argc, char **argv)
{
  try
    {
      const std::vector<nmopt::test_support::Scenario> scenarios{
        {"checked_in_families_expand_and_filter",
         "nmopt.parameter_files.checked_in_families_expand_and_filter",
         {"backend", "dealii", "application", "runner", "contract"},
         30,
         test_checked_in_families_expand_and_filter},
        {"unknown_selection_is_rejected",
         "nmopt.parameter_files.unknown_selection_is_rejected",
         {"backend", "dealii", "application", "runner", "contract", "negative"},
         30,
         test_unknown_selection_is_rejected}};
      const auto result = nmopt::test_support::run_requested_scenarios(
        argc, argv, scenarios, std::cout);
      if (!result.listed)
        std::cout << "parameter file contract tests passed: " << result.executed
                  << '\n';
      return 0;
    }
  catch (const std::exception &error)
    {
      std::cerr << "parameter file contract test failed: " << error.what()
                << '\n';
      return 1;
    }
}
