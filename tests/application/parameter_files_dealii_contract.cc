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
    require(b1.value("Compile/state solve maximum iterations") == "0" &&
              b1.value("Compile/adjoint solve relative tolerance") ==
                "1e-12" &&
              b1.value("Compile/control metric solve maximum iterations") ==
                "1000",
            "B1 parameter family lost its linear-solve policies");
    require(b1.value("Mesh/generator") == "framework-native" &&
              b1.value("Mesh/subdivisions") == "0" &&
              b1.value("Mesh/centroid splits") == "0" &&
              b1.value("Mesh/selection seed") == "0",
            "B1 parameter family lost the backward-compatible mesh defaults");

    const auto figure_6_3 = read_parameter_file(find_file_from_current_or_parent(
      "parameters/chapter-6/b1/development/figure-6.3-book-policy.prm"));
    require(figure_6_3.combinations().size() == 6,
            "Figure 6.3 family should expand to six combinations");
    require(figure_6_3.value("Solver/maximum backtracking reductions") == "5" &&
              figure_6_3.value("Solver/objective target policy") ==
                "match-reference-method",
            "Figure 6.3 family lost the recovered book solver policy");

    const auto continuous_control = read_parameter_file(
      find_file_from_current_or_parent(
        "parameters/chapter-6/b1/development/continuous-control.prm"));
    require(continuous_control.combinations().size() == 6 &&
              continuous_control.value("Problem/control representation") ==
                "continuous-volume-homogeneous-dirichlet" &&
              continuous_control.value("Problem/cellwise box constraint") ==
                "false",
            "B1 continuous-control family lost its candidate discretisation");

    const auto constant_one_forcing = read_parameter_file(
      find_file_from_current_or_parent(
        "parameters/chapter-6/b1/development/continuous-control-constant-one.prm"));
    require(constant_one_forcing.combinations().size() == 6 &&
              constant_one_forcing.value("Functions/forcing") ==
                "figure-inferred-constant-one" &&
              constant_one_forcing.value("Functions/forcing/kind") ==
                "constant" &&
              constant_one_forcing.value("Functions/forcing/value") == "1.0",
            "B1 constant-one family lost its inferred forcing candidate");

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
