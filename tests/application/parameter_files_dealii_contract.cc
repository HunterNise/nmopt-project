#include "../../apps/nmopt-runner/parameter_files.hpp"
#include "../support/scenario_dispatch.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
  using nmopt::application::runner::find_file_from_current_or_parent;
  using nmopt::application::runner::parameter_positive_unsigned_list;
  using nmopt::application::runner::read_parameter_file;
  using nmopt::application::runner::resolve_method_parameter;

  void
  require(const bool condition, const char *message)
  {
    if (!condition)
      throw std::runtime_error(message);
  }

  template <typename Operation>
  void
  require_invalid_argument(Operation &&operation, const char *message)
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

  nmopt::application::runner::ParameterFile
  read_exclusion_parameter_file(const std::string &exclusions)
  {
    const auto path = std::filesystem::temp_directory_path() /
                      "nmopt-parameter-exclusion-contract.prm";
    std::ofstream output(path);
    output << "subsection Benchmark\n"
           << "  set id = b1\n"
           << "  set recipe = chapter-6.b1.distributed-laplace\n"
           << "end\n"
           << "subsection Matrix\n"
           << "  set method = steepest-descent, l-bfgs\n"
           << "  set regularisation = 1e-1, 1e-2\n"
           << "end\n"
           << "subsection Selection\n"
           << "  set exclude combinations = " << exclusions << '\n'
           << "end\n";
    output.close();
    if (!output)
      throw std::runtime_error("could not write exclusion parameter fixture");
    try
      {
        auto result = read_parameter_file(path);
        std::filesystem::remove(path);
        return result;
      }
    catch (...)
      {
        std::filesystem::remove(path);
        throw;
      }
  }

  void
  test_checked_in_families_expand_and_filter()
  {
    const auto b1 = read_parameter_file(find_file_from_current_or_parent(
      "parameters/chapter-6/b1/authoritative.prm"));
    require(b1.matrix.size() == 2, "B1 should declare two matrix axes");
    require(b1.combinations().size() == 7 &&
              b1.excluded_combinations.size() == 1 &&
              b1.excluded_combinations[0].values.at("method") ==
                "steepest-descent" &&
              b1.excluded_combinations[0].values.at("regularisation") ==
                "1e-6",
            "B1 should retain exactly the seven source figure combinations");
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
    require(b1.value("Problem/control representation") ==
                "continuous-volume-homogeneous-dirichlet" &&
              b1.value("Functions/forcing") ==
                "source-oriented-constant-half" &&
              b1.value("Functions/forcing/kind") == "constant" &&
              b1.value("Functions/forcing/value") == "0.5" &&
              b1.value("Mesh/generator") == "structured-simplex" &&
              b1.value("Mesh/refinement") == "0" &&
              b1.value("Mesh/subdivisions") == "131" &&
              b1.value("Mesh/axis subdivisions").empty() &&
              b1.value("Mesh/centroid splits") == "0" &&
              b1.value("Mesh/selection seed") == "0",
            "B1 parameter family lost its source-oriented replacements");
    require(
      b1.value("Solver/maximum backtracking reductions") == "5" &&
        b1.value("Solver/objective target policy") == "none" &&
        resolve_method_parameter(b1,
                                 "steepest-descent",
                                 "stopping criterion")
            .value == "relative-gradient-norm" &&
        resolve_method_parameter(b1,
                                 "steepest-descent",
                                 "relative gradient tolerance")
            .value == "1e-3" &&
        resolve_method_parameter(b1, "l-bfgs", "stopping criterion").value ==
          "relative-gradient-norm" &&
        resolve_method_parameter(b1,
                                 "l-bfgs",
                                 "relative gradient tolerance")
            .value == "1e-3" &&
        resolve_method_parameter(b1,
                                 "l-bfgs",
                                 "initial inverse Hessian scaling")
            .value == "metric-inverse" &&
        resolve_method_parameter(b1, "l-bfgs", "memory").value == "5",
      "B1 parameter family lost its unified source-oriented solver policy");

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
    require(b2.value("Problem/control representation") ==
                "facewise-constant" &&
              b2.value("Solver/globalization") == "armijo",
            "B2 authoritative parameters lost a frozen control or solver choice");

    const auto development = read_parameter_file(find_file_from_current_or_parent(
      "parameters/chapter-6/b2/development/forcing-sweep.prm"));
    require(development.combinations().size() == 3,
            "B2 forcing development family should expand to three combinations");
    require(development.value("Problem/control representation") ==
                "facewise-constant" &&
              development.value("Solver/globalization") == "armijo",
            "B2 forcing family lost a frozen control or solver choice");
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
  test_axis_mesh_subdivisions_are_parsed()
  {
    auto file = read_parameter_file(find_file_from_current_or_parent(
      "parameters/chapter-6/b2/authoritative.prm"));
    require(parameter_positive_unsigned_list(file,
                                             "Mesh/axis subdivisions")
              .empty(),
            "tracked parameter files did not retain the empty axis default");

    file.values["Mesh/axis subdivisions"] = "40, 10";
    require(parameter_positive_unsigned_list(file,
                                             "Mesh/axis subdivisions") ==
              std::vector<unsigned int>({40, 10}),
            "runner parsing lost per-axis subdivision counts");

    file.values["Mesh/axis subdivisions"] = "40, 0";
    require_invalid_argument(
      [&] {
        (void)parameter_positive_unsigned_list(file,
                                               "Mesh/axis subdivisions");
      },
      "runner parsing accepted a zero axis subdivision count");

    file.values["Mesh/axis subdivisions"] = "40,,10";
    require_invalid_argument(
      [&] {
        (void)parameter_positive_unsigned_list(file,
                                               "Mesh/axis subdivisions");
      },
      "runner parsing accepted an empty axis subdivision count");
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
  test_b2_transport_boundary_form_is_selected_coherently()
  {
    const auto ordinary = read_parameter_file(find_file_from_current_or_parent(
      "parameters/chapter-6/b2/authoritative.prm"));
    require(
      nmopt::application::runner::b2_transport_boundary_form(ordinary) ==
        nmopt::semantic::v1::TransportBoundaryForm::
          ordinary_normal_minus_transport,
      "B2 parameter parsing lost the ordinary-normal boundary selection");

    auto total_conormal = ordinary;
    total_conormal.values["Boundary/transport boundary form"] =
      "total-conormal";
    total_conormal.values["Boundary/conormal form"] =
      "diffusion-minus-transport";
    require(
      nmopt::application::runner::b2_transport_boundary_form(total_conormal) ==
        nmopt::semantic::v1::TransportBoundaryForm::total_conormal,
      "B2 parameter parsing lost the total-conormal boundary selection");

    auto inconsistent_total = total_conormal;
    inconsistent_total.values["Boundary/conormal form"] = "unspecified";
    require_invalid_argument(
      [&] {
        (void)nmopt::application::runner::b2_transport_boundary_form(
          inconsistent_total);
      },
      "B2 accepted total conormal without diffusion-minus-transport");

    auto inconsistent_ordinary = ordinary;
    inconsistent_ordinary.values["Boundary/conormal form"] =
      "diffusion-minus-transport";
    require_invalid_argument(
      [&] {
        (void)nmopt::application::runner::b2_transport_boundary_form(
          inconsistent_ordinary);
      },
      "B2 accepted an independently selected ordinary conormal form");

    auto unknown = ordinary;
    unknown.values["Boundary/transport boundary form"] = "book-notation";
    require_invalid_argument(
      [&] {
        (void)nmopt::application::runner::b2_transport_boundary_form(unknown);
      },
      "B2 accepted an unknown transport-boundary form");
  }

  void
  test_reduced_globalization_is_selected()
  {
    const auto armijo = read_parameter_file(find_file_from_current_or_parent(
      "parameters/chapter-6/b2/authoritative.prm"));
    require(
      nmopt::application::runner::reduced_globalization(armijo) ==
        nmopt::application::chapter6::ReducedGlobalization::armijo,
      "runner parsing lost the default Armijo globalization");

    auto legacy_default = armijo;
    legacy_default.values.erase("Solver/globalization");
    require(
      nmopt::application::runner::reduced_globalization(legacy_default) ==
        nmopt::application::chapter6::ReducedGlobalization::armijo,
      "runner parsing did not preserve the legacy Armijo default");

    auto fixed_step = armijo;
    fixed_step.values["Solver/globalization"] = "fixed-step";
    require(
      nmopt::application::runner::reduced_globalization(fixed_step) ==
        nmopt::application::chapter6::ReducedGlobalization::fixed_step,
      "runner parsing lost the fixed-step globalization");

    auto unknown = armijo;
    unknown.values["Solver/globalization"] = "constant";
    require_invalid_argument(
      [&] {
        (void)nmopt::application::runner::reduced_globalization(unknown);
      },
      "runner parsing accepted an unknown globalization");
  }

  void
  test_b2_control_discretisation_is_selected()
  {
    const auto facewise = read_parameter_file(find_file_from_current_or_parent(
      "parameters/chapter-6/b2/authoritative.prm"));
    require(
      nmopt::application::runner::b2_neumann_control_discretisation(facewise) ==
        nmopt::semantic::v1::NeumannControlDiscretisation::facewise_constant,
      "B2 parameter parsing lost the facewise control selection");

    auto continuous = facewise;
    continuous.values["Problem/control representation"] =
      "continuous-nodal-trace";
    require(
      nmopt::application::runner::b2_neumann_control_discretisation(
        continuous) ==
        nmopt::semantic::v1::NeumannControlDiscretisation::
          continuous_nodal_trace,
      "B2 parameter parsing lost the continuous control selection");

    auto unknown = facewise;
    unknown.values["Problem/control representation"] = "cellwise-boundary";
    require_invalid_argument(
      [&] {
        (void)nmopt::application::runner::
          b2_neumann_control_discretisation(unknown);
      },
      "B2 accepted an unknown control representation");
  }

  void
  test_sparse_matrix_exclusions_are_validated_and_applied()
  {
    const std::string excluded =
      "[method=steepest-descent,regularisation=1e-2]";
    const auto file = read_exclusion_parameter_file(excluded);
    const auto combinations = file.combinations();
    require(file.excluded_combinations.size() == 1,
            "parameter file did not retain its excluded coordinate");
    require(combinations.size() == 3,
            "excluded matrix coordinate was not filtered");
    require(combinations[0].values.at("method") == "steepest-descent" &&
              combinations[0].values.at("regularisation") == "1e-1" &&
              combinations[1].values.at("method") == "l-bfgs" &&
              combinations[1].values.at("regularisation") == "1e-1" &&
              combinations[2].values.at("method") == "l-bfgs" &&
              combinations[2].values.at("regularisation") == "1e-2",
            "sparse matrix expansion lost its declared ordering");

    require_invalid_argument(
      [&] {
        (void)file.combinations({{"method", "steepest-descent"},
                                 {"regularisation", "1e-2"}});
      },
      "a CLI selection resolving only to an exclusion was accepted");

    for (const auto &invalid :
         {"[method=steepest-descent]",
          "[method=steepest-descent,regularisation=1e-3]",
          "[method=steepest-descent,unknown=1e-2]",
          "[method=steepest-descent,method=l-bfgs,regularisation=1e-2]",
          "method=steepest-descent,regularisation=1e-2",
          "[method=steepest-descent,regularisation=1e-2];"
          "[regularisation=1e-2,method=steepest-descent]"})
      require_invalid_argument(
        [&] { (void)read_exclusion_parameter_file(invalid); },
        "an invalid excluded matrix coordinate was accepted");
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
        {"axis_mesh_subdivisions_are_parsed",
         "nmopt.parameter_files.axis_mesh_subdivisions_are_parsed",
         {"backend", "dealii", "application", "runner", "contract", "negative"},
         30,
         test_axis_mesh_subdivisions_are_parsed},
        {"b2_transport_boundary_form_is_selected_coherently",
         "nmopt.parameter_files.b2_transport_boundary_form_is_selected_coherently",
         {"backend", "dealii", "application", "runner", "contract"},
         30,
         test_b2_transport_boundary_form_is_selected_coherently},
        {"b2_control_discretisation_is_selected",
         "nmopt.parameter_files.b2_control_discretisation_is_selected",
         {"backend", "dealii", "application", "runner", "contract"},
         30,
         test_b2_control_discretisation_is_selected},
        {"reduced_globalization_is_selected",
         "nmopt.parameter_files.reduced_globalization_is_selected",
         {"backend", "dealii", "application", "runner", "contract", "negative"},
         30,
         test_reduced_globalization_is_selected},
        {"sparse_matrix_exclusions_are_validated_and_applied",
         "nmopt.parameter_files.sparse_matrix_exclusions_are_validated_and_applied",
         {"backend", "dealii", "application", "runner", "contract"},
         30,
         test_sparse_matrix_exclusions_are_validated_and_applied},
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
