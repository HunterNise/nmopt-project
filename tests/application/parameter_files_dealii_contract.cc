#include "../../apps/nmopt-runner/parameter_files.hpp"
#include "../support/scenario_dispatch.hpp"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
  using nmopt::application::runner::find_file_from_current_or_parent;
  using nmopt::application::runner::read_parameter_file;
  using nmopt::application::runner::resolve_method_parameter;

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
    const auto steepest_stopping = resolve_method_parameter(
      figure_6_3, "steepest-descent", "stopping criterion");
    const auto steepest_relative = resolve_method_parameter(
      figure_6_3, "steepest-descent", "relative gradient tolerance");
    const auto l_bfgs_stopping =
      resolve_method_parameter(figure_6_3, "l-bfgs", "stopping criterion");
    require(steepest_stopping.value == "relative-gradient-norm" &&
              steepest_relative.value == "1e-3" &&
              l_bfgs_stopping.value == "gradient-norm",
            "Figure 6.3 methods did not resolve distinct stopping policies");
    auto method_override = figure_6_3;
    method_override.values[
      "Solver/method policy steepest-descent/initial step length"] = "10.0";
    require(resolve_method_parameter(method_override,
                                     "steepest-descent",
                                     "initial step length")
                .value == "10.0" &&
              resolve_method_parameter(method_override,
                                       "l-bfgs",
                                       "initial step length")
                .value == "1.0",
            "method policy precedence did not preserve the global fallback");

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

    const auto structured_simplex = read_parameter_file(
      find_file_from_current_or_parent(
        "parameters/chapter-6/b1/development/continuous-control-structured-simplex.prm"));
    require(structured_simplex.combinations().size() == 6 &&
              structured_simplex.value("Mesh/generator") ==
                "structured-simplex" &&
              structured_simplex.value("Mesh/refinement") == "0" &&
              structured_simplex.value("Mesh/subdivisions") == "131" &&
              structured_simplex.value("Mesh/centroid splits") == "0" &&
              structured_simplex.value("Problem/control representation") ==
                "continuous-volume-homogeneous-dirichlet" &&
              structured_simplex.value("Functions/forcing") ==
                "manufactured-zero" &&
              structured_simplex.value("Solver/objective target policy") ==
                "match-reference-method",
            "B1 structured-simplex family lost its mesh candidate");

    const auto count_matched_simplex = read_parameter_file(
      find_file_from_current_or_parent(
        "parameters/chapter-6/b1/development/continuous-control-count-matched-simplex.prm"));
    require(count_matched_simplex.combinations().size() == 6 &&
              count_matched_simplex.value("Mesh/generator") ==
                "centroid-split-simplex" &&
              count_matched_simplex.value("Mesh/subdivisions") == "100" &&
              count_matched_simplex.value("Mesh/centroid splits") == "7160" &&
              count_matched_simplex.value("Mesh/selection seed") == "0" &&
              count_matched_simplex.value("Mesh/provenance") ==
                "chapter-6.e6.5.1.count-matched-simplex-n100-s7160-seed0" &&
              count_matched_simplex.value("Functions/forcing") ==
                "manufactured-zero" &&
              count_matched_simplex.value("Solver/objective target policy") ==
                "match-reference-method",
            "B1 count-matched-simplex family lost its topology hypothesis");

    const auto figure_6_2_constant_half = read_parameter_file(
      find_file_from_current_or_parent(
        "parameters/chapter-6/b1/development/"
        "figure-6.2-early-stop-constant-half.prm"));
    require(
      figure_6_2_constant_half.combinations().size() == 2 &&
        figure_6_2_constant_half.value("Functions/forcing/kind") == "constant" &&
        figure_6_2_constant_half.value("Functions/forcing/value") == "0.5" &&
        figure_6_2_constant_half.value("Mesh/generator") ==
          "structured-simplex" &&
        figure_6_2_constant_half.value("Mesh/subdivisions") == "131" &&
        figure_6_2_constant_half.value("Solver/objective target policy") ==
          "none" &&
        resolve_method_parameter(figure_6_2_constant_half,
                                 "l-bfgs",
                                 "relative gradient tolerance")
            .value == "1e-3",
      "B1 Figure 6.2 constant-half family lost its early-stop hypothesis");

    const auto figure_6_2_objective_matched = read_parameter_file(
      find_file_from_current_or_parent(
        "parameters/chapter-6/b1/development/"
        "figure-6.2-early-stop-objective-matched.prm"));
    require(
      figure_6_2_objective_matched.combinations().size() == 2 &&
        figure_6_2_objective_matched.value("Functions/forcing") ==
          "objective-matched-constant" &&
        figure_6_2_objective_matched.value("Functions/forcing/kind") ==
          "constant" &&
        figure_6_2_objective_matched.value("Functions/forcing/value") ==
          "0.41506741762176758" &&
        figure_6_2_objective_matched.value("Mesh/generator") ==
          "structured-simplex" &&
        figure_6_2_objective_matched.value("Mesh/subdivisions") == "131" &&
        figure_6_2_objective_matched.value(
          "Solver/objective target policy") == "none" &&
        resolve_method_parameter(figure_6_2_objective_matched,
                                 "l-bfgs",
                                 "relative gradient tolerance")
            .value == "1e-3",
      "B1 Figure 6.2 objective-matched family lost its forcing hypothesis");

    const auto require_figure_6_3_candidate = [](const auto &candidate,
                                                  const char *forcing_value,
                                                  const char *message) {
      require(
        candidate.combinations().size() == 6 &&
          candidate.value("Functions/forcing/kind") == "constant" &&
          candidate.value("Functions/forcing/value") == forcing_value &&
          candidate.value("Mesh/generator") == "structured-simplex" &&
          candidate.value("Mesh/subdivisions") == "131" &&
          candidate.value("Solver/objective target policy") ==
            "match-reference-method" &&
          resolve_method_parameter(candidate,
                                   "steepest-descent",
                                   "relative gradient tolerance")
              .value == "1e-3" &&
          resolve_method_parameter(candidate, "l-bfgs", "stopping criterion")
              .value == "gradient-norm",
        message);
    };

    const auto figure_6_3_constant_half = read_parameter_file(
      find_file_from_current_or_parent(
        "parameters/chapter-6/b1/development/figure-6.3-constant-half.prm"));
    require_figure_6_3_candidate(
      figure_6_3_constant_half,
      "0.5",
      "B1 Figure 6.3 constant-half family lost its method comparison");

    const auto figure_6_3_objective_matched = read_parameter_file(
      find_file_from_current_or_parent(
        "parameters/chapter-6/b1/development/figure-6.3-objective-matched.prm"));
    require_figure_6_3_candidate(
      figure_6_3_objective_matched,
      "0.41506741762176758",
      "B1 Figure 6.3 objective-matched family lost its method comparison");

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

    const auto chapter_6_parameters =
      find_file_from_current_or_parent("parameters/chapter-6");
    for (const auto &entry :
         std::filesystem::recursive_directory_iterator(chapter_6_parameters))
      if (entry.is_regular_file() && entry.path().extension() == ".prm")
        {
          const auto stable = read_parameter_file(entry.path());
          if (stable.value("Run/output root") != "runs")
            throw std::runtime_error(
              "tracked parameter file must use the stable runs root: " +
              entry.path().string());
          if (entry.path().parent_path().filename() == "development" &&
              !stable.optional_value("Solver/declared minimum step length")
                 .empty())
            throw std::runtime_error(
              "development parameter file uses a legacy minimum-step entry: " +
              entry.path().string());
        }
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

  void
  test_scalar_function_definitions_are_data_driven()
  {
    nmopt::application::runner::ParameterFile file;
    file.values = {{"Functions/forcing", "candidate-045"},
                   {"Functions/forcing/kind", "constant"},
                   {"Functions/forcing/value", "0.45"},
                   {"Functions/forcing/expression", ""},
                   {"Functions/forcing/provenance",
                    "development.b1.candidate-045"}};
    const auto constant =
      nmopt::application::runner::parameter_scalar_function_definition(
        file, "Functions/forcing");
    require(constant.id == "candidate-045" &&
              constant.kind ==
                nmopt::application::ScalarFunctionKind::constant &&
              std::abs(constant.value - 0.45) < 1.0e-15,
            "parameter parsing hardcoded the constant forcing value");

    file.values["Functions/forcing"] = "spatial-candidate";
    file.values["Functions/forcing/kind"] = "expression";
    file.values["Functions/forcing/value"] = "";
    file.values["Functions/forcing/expression"] =
      "0.4 + x0*(1-x0)*x1*(1-x1)";
    file.values["Functions/forcing/provenance"] =
      "development.b1.spatial-candidate";
    const auto expression =
      nmopt::application::runner::parameter_scalar_function_definition(
        file, "Functions/forcing");
    require(expression.id == "spatial-candidate" &&
              expression.kind ==
                nmopt::application::ScalarFunctionKind::expression &&
              expression.expression ==
                "0.4 + x0*(1-x0)*x1*(1-x1)",
            "parameter parsing did not retain a forcing expression");

    file.values["Functions/forcing/value"] = "0.45";
    bool ambiguous_definition_rejected = false;
    try
      {
        (void)nmopt::application::runner::
          parameter_scalar_function_definition(file, "Functions/forcing");
      }
    catch (const std::invalid_argument &)
      {
        ambiguous_definition_rejected = true;
      }
    require(ambiguous_definition_rejected,
            "parameter parsing accepted both forcing value and expression");

    file.values["Functions/forcing/value"] = "";
    file.values["Functions/forcing/kind"] = "registered-only";
    bool unknown_kind_rejected = false;
    try
      {
        (void)nmopt::application::runner::
          parameter_scalar_function_definition(file, "Functions/forcing");
      }
    catch (const std::invalid_argument &)
      {
        unknown_kind_rejected = true;
      }
    require(unknown_kind_rejected,
            "parameter parsing accepted an unknown scalar function kind");
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
         test_unknown_selection_is_rejected},
        {"scalar_function_definitions_are_data_driven",
         "nmopt.parameter_files.scalar_function_definitions_are_data_driven",
         {"backend", "dealii", "application", "runner", "contract"},
         30,
         test_scalar_function_definitions_are_data_driven}};
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
