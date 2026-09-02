#include "../../apps/nmopt-runner/parameter_files.hpp"
#include "../support/scenario_dispatch.hpp"

// These tests characterize the current production resolution path through the
// runner entry point while preserving the behavior being migrated.
#define main nmopt_runner_characterization_entrypoint
#include "../../apps/nmopt-runner/main.cc"
#undef main

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{
  using nmopt::application::runner::find_file_from_current_or_parent;
  using nmopt::application::runner::parameter_finite_list;
  using nmopt::application::runner::parameter_positive_unsigned_list;
  using nmopt::application::runner::read_parameter_file;
  using nmopt::application::runner::resolve_method_parameter;

  void
  require(const bool condition, const char *message)
  {
    if (!condition)
      throw std::runtime_error(message);
  }

  std::string
  forcing_definition_path()
  {
    return "Functions/forcing";
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

  nmopt::application::chapter6::B1Scenario
  resolve_b1_scenario_for_characterization(
    const nmopt::application::runner::ParameterFile &file,
    const nmopt::application::runner::ParameterCombination &combination)
  {
    const auto method_id = combination_value(combination, "method");
    const auto method = parse_method(method_id);
    const auto beta = parse_number_text(
      combination_value(combination, "regularisation"),
      "Matrix/regularisation");
    auto scenario = nmopt::application::chapter6::make_b1_scenario(method);
    bind_b1_scenario(
      scenario,
      file,
      combination,
      method_id,
      beta,
      std::string(file.value("Benchmark/id")) + "." + method_id + ".beta-" +
        combination_value(combination, "regularisation"));
    return scenario;
  }

  nmopt::application::chapter6::B2Scenario
  resolve_b2_scenario_for_characterization(
    const nmopt::application::runner::ParameterFile &file,
    const nmopt::application::runner::ParameterCombination &combination)
  {
    const auto &observation_region =
      combination_value(combination, "observation-region");
    const auto target_profile = combination_value(combination, "target-profile");
    auto scenario = nmopt::application::chapter6::make_b2_scenario(
      observation_region, target_profile);
    bind_b2_scenario(
      scenario,
      file,
      combination,
      std::string("chapter-6.b2.graetz-flow.") +
        nmopt::application::chapter6::b2_case_name(observation_region,
                                                   target_profile));
    return scenario;
  }

  void
  require_binding_provenance(
    const nmopt::compiler::v1::CompilationManifest &manifest,
    const std::string                             &semantic_id,
    const std::string                             &provenance,
    const char                                    *message)
  {
    require(
      std::any_of(manifest.bindings.begin(),
                  manifest.bindings.end(),
                  [&](const auto &binding) {
                    return binding.semantic_id == semantic_id &&
                           binding.provenance == provenance;
                  }),
      message);
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

  nmopt::application::runner::ParameterFile
  read_b2_scalar_discovery_parameter_file(
    const std::string_view observation_expression =
      "x0 > 1.0 ? (x1 < 0.3 ? 1 : (x1 > 0.7 ? 1 : 0)) : 0",
    const unsigned int fixed_boundary_id = 0,
    const unsigned int control_boundary_id = 1,
    const unsigned int outflow_boundary_id = 2)
  {
    const auto path = std::filesystem::temp_directory_path() /
                      "nmopt-b2-scalar-discovery-contract.prm";
    std::ofstream output(path);
    output << R"prm(
subsection Benchmark
  set id = chapter-6.b2.graetz-flow
  set recipe = chapter-5.scalar-neumann-convection-subdomain
  set source reference = test scalar discovery
  set source revision = test
end

subsection Matrix
  set forcing = zero, spatial-candidate
  set regularisation = 1e-3
  set observation-region = wings
  set target-profile = candidate-target
end

subsection Problem
  set control representation = facewise-constant
  set facewise box constraint = false
end

subsection Observation
  set material id = 1
  set realization = cell-center-indicator

  subsection region definitions
    subsection wings
      set kind = expression
      set expression = )prm"
           << observation_expression << R"prm(
      set provenance = test.wings-observation
    end
  end
end

subsection Functions
  set fixed Dirichlet data = fixed-temperature
  set conservative transport = graetz

  subsection forcing definitions
    subsection zero
      set kind = zero
      set provenance = test.zero-forcing
    end

    subsection spatial-candidate
      set kind = expression
      set expression = 0.25 + sin(pi*x0)*sin(pi*x1)
      set provenance = test.spatial-candidate-forcing
    end
  end

  subsection fixed Dirichlet data
    set kind = constant
    set value = 1.0
    set provenance = test.fixed-temperature
  end

  subsection conservative transport
    set kind = conservative-transport
    set expression = 1.5*x1*(1-x1); 0.0
    set provenance = test.graetz-transport
  end

  subsection target definitions
    subsection candidate-target
      set kind = expression
      set expression = x0 + 2.0*x1
      set provenance = test.target-candidate
    end
  end
end

subsection Runtime
  set diffusion = 0.1
  set reaction = 0.0
  set regularisation = 1e-3
end

subsection Boundary
  set fixed id = )prm"
           << fixed_boundary_id << R"prm(
  set control id = )prm"
           << control_boundary_id << R"prm(
  set outflow id = )prm"
           << outflow_boundary_id << R"prm(
  set upstream transition = 1.0
  set transport boundary form = ordinary-normal-minus-transport
  set conormal form = unspecified
  set normal orientation = outward
  set trace evaluation = fe-q-state-trace
  set face quadrature = qgauss-face
end

subsection Mesh
  set dimension = 2
  set geometry = rectangle
  set lower = (0.0, 0.0)
  set upper = (4.0, 1.0)
  set refinement = 1
  set provenance = test.rectangle
end

subsection Compile
  set state degree = 1
  set volume observation quadrature order = 3
  set volume observation target realisation = analytic-quadrature
  set execution = assembled
  set product = reduced-dto
  set owned session = true
  set state solve maximum iterations = 0
  set state solve relative tolerance = 1e-12
  set state solve absolute tolerance = 1e-14
  set adjoint solve maximum iterations = 0
  set adjoint solve relative tolerance = 1e-12
  set adjoint solve absolute tolerance = 1e-14
  set control metric solve maximum iterations = 1000
  set control metric solve relative tolerance = 1e-12
  set control metric solve absolute tolerance = 1e-14
  set stabilization = galerkin
end

subsection Solver
  set method = bfgs
  set globalization = armijo
  set initial independent control value = 0.0
  set maximum iterations = 2
  set maximum line search trials = 20
  set gradient tolerance = 1e-8
  set stopping criterion = automatic
  set relative gradient tolerance = 0.0
  set objective change tolerance = 0.0
  set step tolerance = 0.0
  set initial step length = 1.0
  set Armijo fraction = 1e-4
  set backtracking factor = 0.5

  subsection method policy bfgs
    set curvature tolerance = 1e-14
  end
end

subsection Run
  set kind = development
  set build profile = debug-dealii
  set output root = runs
  set deterministic = true
  set serialize artifacts = true
  set measure timings = false
  set measure memory = false
end

subsection Output
  set retain fields = true
  set selected fields = state, state-uncontrolled, control, adjoint, negative-adjoint, target, forcing, observation-region
end
)prm";
    output.close();
    if (!output)
      throw std::runtime_error("could not write scalar discovery fixture");
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
  test_parameter_schema_registry_accounts_for_all_entries()
  {
    const auto &common_registry =
      nmopt::application::runner::detail::common_parameter_schema_registry();
    require(std::none_of(common_registry.begin(),
                         common_registry.end(),
                         [](const auto &entry) {
                           return entry.path.rfind("Matrix/", 0) == 0;
                         }),
            "common parameter schema must not own benchmark matrix axes");

    const auto b1_adapter =
      nmopt::application::runner::detail::parameter_schema_adapter_for(
        "b1", "chapter-6.b1.distributed-laplace");
    const auto registry =
      nmopt::application::runner::detail::parameter_schema_registry(b1_adapter);
    require(!registry.empty(), "parameter schema registry must not be empty");

    const auto exclusion_file = read_exclusion_parameter_file("");
    std::size_t required_entries = 0;
    for (const auto &entry : registry)
      {
        require(!entry.path.empty() && entry.pattern,
                "every schema entry needs a path and deal.II pattern");
        require(exclusion_file.values.count(entry.path) == 1,
                "schema extraction must account for every declared entry");
        required_entries += entry.presence ==
                            nmopt::application::runner::ParameterPresence::required;
      }
    require(exclusion_file.values.size() > registry.size() && required_entries == 2,
            "schema registry and extraction accounting is incomplete");
    require(
      exclusion_file.values.count(
        "Solver/method policy steepest-descent/maximum iterations") == 1 &&
        exclusion_file.values.count("Solver/method policy l-bfgs/memory") == 1 &&
        exclusion_file.values.count(
          "Solver/method policy l-bfgs/initial inverse Hessian scaling") == 1 &&
        exclusion_file.values.count(
          "Solver/method policy bfgs/curvature tolerance") == 0 &&
        exclusion_file.values.count(
          "Solver/method policy steepest-descent/memory") == 0,
      "method-specific schema discovery exposed the wrong solver leaves");

    const auto find_entry = [&](const std::string &path) -> const auto & {
      const auto found = std::find_if(
        registry.begin(), registry.end(), [&](const auto &entry) {
          return entry.path == path;
        });
      require(found != registry.end(), "expected typed schema entry is missing");
      return *found;
    };
    require(find_entry("Run/deterministic").pattern->match("true") &&
              !find_entry("Run/deterministic").pattern->match("maybe"),
            "boolean schema pattern should reject non-boolean values");
    require(find_entry("Mesh/refinement").pattern->match("7") &&
              !find_entry("Mesh/refinement").pattern->match("seven"),
            "integer schema pattern should reject non-integer values");
    require(find_entry("Runtime/diffusion").pattern->match("1.25") &&
              !find_entry("Runtime/diffusion").pattern->match("not-a-number"),
            "double schema pattern should reject non-numeric values");
    require(find_entry("Solver/globalization").pattern->match("armijo") &&
              !find_entry("Solver/globalization").pattern->match("unknown"),
            "selection schema pattern should reject unknown values");

    const auto b2_adapter =
      nmopt::application::runner::detail::parameter_schema_adapter_for(
        "chapter-6.b2.graetz-flow",
        "chapter-5.scalar-neumann-convection-subdomain");
    const auto b2_registry =
      nmopt::application::runner::detail::parameter_schema_registry(b2_adapter);
    require(std::any_of(b2_registry.begin(), b2_registry.end(), [](const auto &entry) {
              return entry.path == "Matrix/forcing";
            }) &&
              std::none_of(b1_adapter.entries.begin(),
                           b1_adapter.entries.end(),
                           [](const auto &entry) {
                             return entry.path == "Matrix/forcing";
                           }),
            "benchmark schema adapters must own their matrix axes");

    const auto b2_file = read_parameter_file(find_file_from_current_or_parent(
      "parameters/chapter-6/b2/authoritative.prm"));
    require(
      b2_file.values.count(
        "Solver/method policy bfgs/curvature tolerance") == 1 &&
        b2_file.values.count("Solver/method policy l-bfgs/memory") == 0 &&
        b2_file.values.count("Solver/method policy l-bfgs/initial inverse Hessian scaling") ==
          0,
      "B2 schema discovery exposed L-BFGS policy leaves for full BFGS");
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
              b1.value(forcing_definition_path() + "/kind") == "constant" &&
              b1.value(forcing_definition_path() + "/value") == "0.5" &&
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
              constant_one_forcing.value(
                forcing_definition_path() + "/kind") ==
                "constant" &&
              constant_one_forcing.value(
                forcing_definition_path() + "/value") ==
                "1.0",
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
        figure_6_2_constant_half.value(
          forcing_definition_path() + "/kind") ==
          "constant" &&
        figure_6_2_constant_half.value(
          forcing_definition_path() + "/value") ==
          "0.5" &&
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
        figure_6_2_objective_matched.value(
          forcing_definition_path() + "/kind") ==
          "constant" &&
        figure_6_2_objective_matched.value(
          forcing_definition_path() + "/value") ==
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
      const auto forcing_path = forcing_definition_path();
      require(
        candidate.combinations().size() == 6 &&
          candidate.value(forcing_path + "/kind") == "constant" &&
          candidate.value(forcing_path + "/value") == forcing_value &&
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
              b2.value("Solver/globalization") == "armijo" &&
              b2.value("Compile/volume observation quadrature order") == "3" &&
              b2.value("Compile/volume observation target realisation") ==
                "analytic-quadrature",
            "B2 authoritative parameters lost a frozen discrete choice");

    const auto development = read_parameter_file(find_file_from_current_or_parent(
      "parameters/chapter-6/b2/development/forcing-sweep.prm"));
    require(development.combinations().size() == 3,
            "B2 forcing development family should expand to three combinations");
    require(development.value("Problem/control representation") ==
                "facewise-constant" &&
              development.value("Solver/globalization") == "armijo" &&
              development.value(
                "Compile/volume observation quadrature order") == "3" &&
              development.value(
                "Compile/volume observation target realisation") ==
                "analytic-quadrature",
            "B2 forcing family lost a frozen discrete choice");
    require(development.content_hash.rfind("fnv1a64:", 0) == 0,
            "parameter provenance should carry a labelled deterministic hash");

    const auto require_b2_replication_candidate =
      [](const auto &candidate,
         const char *forcing_value,
         const char *message) {
        const auto forcing_path = forcing_definition_path();
        require(
          candidate.combinations().size() == 4 &&
            candidate.value("Problem/control representation") ==
              "continuous-nodal-trace" &&
            candidate.value(forcing_path + "/kind") == "constant" &&
            candidate.value(forcing_path + "/value") == forcing_value &&
            candidate.value("Runtime/regularisation") == "1e-2" &&
            candidate.value("Mesh/generator") == "structured-simplex" &&
            candidate.value("Mesh/axis subdivisions") == "160, 40" &&
            candidate.value("Solver/globalization") == "fixed-step" &&
            candidate.value("Solver/maximum iterations") == "100" &&
            candidate.value("Solver/initial step length") == "0.05" &&
            candidate.value("Run/kind") == "development" &&
            candidate.value("Run/build profile") == "debug-dealii" &&
            candidate.value("Output/retain fields") == "true",
          message);
      };

    const auto figure_6_5_state_fit = read_parameter_file(
      find_file_from_current_or_parent(
        "parameters/chapter-6/b2/development/figure-6.5-state-fit.prm"));
    require_b2_replication_candidate(
      figure_6_5_state_fit,
      "0.47009",
      "B2 Figure 6.5 candidate lost its fitted-state hypothesis");

    const auto table_6_2_order_fit = read_parameter_file(
      find_file_from_current_or_parent(
        "parameters/chapter-6/b2/development/table-6.2-order-fit.prm"));
    require_b2_replication_candidate(
      table_6_2_order_fit,
      "1.0",
      "B2 Table 6.2 candidate lost its magnitude hypothesis");

    const auto figure_table_parabolic_fit = read_parameter_file(
      find_file_from_current_or_parent(
        "parameters/chapter-6/b2/development/"
        "figure-6.5-table-6.2-parabolic-fit.prm"));
    require_b2_replication_candidate(
      figure_table_parabolic_fit,
      "0.65",
      "B2 combined candidate lost its parabolic-fit hypothesis");

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
  test_mesh_bounds_are_data_driven_and_validated()
  {
    auto file = read_parameter_file(find_file_from_current_or_parent(
      "parameters/chapter-6/b2/authoritative.prm"));
    file.values["Mesh/lower"] = "(-1.5, 0.25)";
    file.values["Mesh/upper"] = "(3.5, 1.75)";
    file.values["Mesh/refinement"] = "0";
    const auto combination = file.combinations().front();
    const auto scenario = resolve_b2_scenario_for_characterization(
      file, combination);
    require(scenario.compile.mesh.lower == std::vector<double>({-1.5, 0.25}) &&
              scenario.compile.mesh.upper == std::vector<double>({3.5, 1.75}),
            "B2 binding did not retain the configured mesh bounds");

    const auto session =
      nmopt::application::chapter6::dealii::make_b2_compilation_session<2>(
        scenario);
    const auto &mesh = session->triangulation();
    double min_x = std::numeric_limits<double>::max();
    double min_y = std::numeric_limits<double>::max();
    double max_x = std::numeric_limits<double>::lowest();
    double max_y = std::numeric_limits<double>::lowest();
    for (unsigned int vertex = 0; vertex < mesh.n_vertices(); ++vertex)
      {
        min_x = std::min(min_x, mesh.get_vertices()[vertex][0]);
        min_y = std::min(min_y, mesh.get_vertices()[vertex][1]);
        max_x = std::max(max_x, mesh.get_vertices()[vertex][0]);
        max_y = std::max(max_y, mesh.get_vertices()[vertex][1]);
      }
    require(std::abs(min_x + 1.5) < 1.0e-15 &&
              std::abs(min_y - 0.25) < 1.0e-15 &&
              std::abs(max_x - 3.5) < 1.0e-15 &&
              std::abs(max_y - 1.75) < 1.0e-15,
            "B2 mesh construction ignored the configured bounds");

    auto invalid_size = file;
    invalid_size.values["Mesh/lower"] = "(-1.5)";
    const auto invalid_size_scenario =
      resolve_b2_scenario_for_characterization(
        invalid_size, invalid_size.combinations().front());
    require_invalid_argument(
      [&] {
        nmopt::application::chapter6::validate_b2(invalid_size_scenario);
      },
      "B2 accepted mesh bounds with the wrong dimension");

    auto invalid_order = file;
    invalid_order.values["Mesh/lower"] = "(2.0, 0.25)";
    invalid_order.values["Mesh/upper"] = "(1.0, 1.75)";
    const auto invalid_order_scenario =
      resolve_b2_scenario_for_characterization(
        invalid_order, invalid_order.combinations().front());
    require_invalid_argument(
      [&] {
        nmopt::application::chapter6::validate_b2(invalid_order_scenario);
      },
      "B2 accepted mesh bounds with reversed coordinates");

    auto invalid_finite = file;
    invalid_finite.values["Mesh/lower"] = "(nan, 0.25)";
    require_invalid_argument(
      [&] {
        (void)parameter_finite_list(invalid_finite, "Mesh/lower");
      },
      "runner parsing accepted a non-finite mesh bound");
  }

  void
  test_b1_beta_coordinates_are_text_stable()
  {
    for (const auto *value : {"1e-1", "1e-2", "1e-3", "1e-6", "5e-4"})
      require(b1_beta_coordinate(value) == value,
              "B1 beta artifact coordinates changed their matrix spelling");

    nmopt::application::runner::RunSetCombination combination;
    combination.values.values = {
      {"method", "l-bfgs"}, {"regularisation", "5e-4"}};
    require(b1_artifact_coordinate_components({}, combination) ==
              std::vector<std::string>({"l-bfgs", "beta-5e-4"}),
            "B1 did not derive an arbitrary beta coordinate from matrix text");
  }

  void
  test_b1_authoritative_resolution_is_characterized()
  {
    const auto file = read_parameter_file(find_file_from_current_or_parent(
      "parameters/chapter-6/b1/authoritative.prm"));
    const auto combinations = file.combinations();
    require(combinations.size() == 7 &&
              file.excluded_combinations.size() == 1,
            "B1 characterization lost its seven-case sparse matrix");

    const std::vector<std::string> expected_paths{
      "artifacts/steepest-descent/beta-1e-1/artifact.kv",
      "artifacts/steepest-descent/beta-1e-2/artifact.kv",
      "artifacts/steepest-descent/beta-1e-3/artifact.kv",
      "artifacts/l-bfgs/beta-1e-1/artifact.kv",
      "artifacts/l-bfgs/beta-1e-2/artifact.kv",
      "artifacts/l-bfgs/beta-1e-3/artifact.kv",
      "artifacts/l-bfgs/beta-1e-6/artifact.kv"};
    const auto b1_plan = nmopt::application::runner::make_run_set_plan(file);
    require(nmopt::application::runner::run_set_artifact_paths(
              b1_plan, b1_artifact_coordinate_components) ==
              expected_paths,
            "B1 characterization changed artifact coordinate paths");

    for (const auto &combination : combinations)
      {
        const auto method_id = combination_value(combination, "method");
        const auto beta = parse_number_text(
          combination_value(combination, "regularisation"),
          "Matrix/regularisation");
        const auto scenario =
          resolve_b1_scenario_for_characterization(file, combination);

        require(scenario.metadata.id ==
                  "chapter-6.b1.distributed-laplace" &&
                  scenario.metadata.recipe_id ==
                    nmopt::application::chapter6::b1_recipe_id,
                "B1 resolution changed the benchmark or recipe identity");
        require(scenario.experiment.scenario_output_id ==
                  std::string("chapter-6.b1.distributed-laplace.") +
                    method_id + ".beta-" + combination_value(
                      combination, "regularisation"),
                "B1 resolution changed the scenario output identity");
        require(scenario.experiment.source_reference ==
                  "E6.5.1 equations (6.64), Figures 6.2-6.3" &&
                  scenario.experiment.source_revision ==
                    "1aaefbe473f9941a89d1df36192251511c052933" &&
                  scenario.experiment.harness.artifact_directory == "runs" &&
                  scenario.experiment.harness.deterministic &&
                  scenario.experiment.harness.serialize_artifacts &&
                  scenario.experiment.harness.measure_timings &&
                  !scenario.experiment.harness.measure_memory &&
                  scenario.experiment.retain_fields,
                "B1 resolution changed run provenance or output policy");

        require(
          scenario.problem.recipe.discretisation ==
              nmopt::application::chapter5::DistributedControlDiscretisation::
                homogeneous_dirichlet_continuous &&
            !scenario.problem.recipe.with_cellwise_box &&
            scenario.problem.data.diffusion == 1.0 &&
            scenario.problem.data.reaction == 0.0 &&
            scenario.problem.data.regularisation_weight == beta &&
            scenario.problem.regularisation_sweep.size() == 1 &&
            scenario.problem.regularisation_sweep.front() == beta &&
            scenario.problem.data.desired_state_provenance ==
              "chapter-6.e6.5.1.desired-state" &&
            scenario.problem.forcing.id == "source-oriented-constant-half" &&
            scenario.problem.forcing.kind ==
              nmopt::application::ScalarFunctionKind::constant &&
            scenario.problem.forcing.value == 0.5 &&
            scenario.problem.forcing.expression.empty() &&
            scenario.problem.data.forcing_provenance ==
              "chapter-6.e6.5.1.source-oriented-constant-half-forcing",
          "B1 resolution changed effective runtime data or control choice");

        require(
          scenario.compile.mesh.dimension == 2 &&
            scenario.compile.mesh.generation ==
              nmopt::application::chapter6::MeshGeneration::structured_simplex &&
            scenario.compile.mesh.refinement == 0 &&
            scenario.compile.mesh.subdivisions == 131 &&
            scenario.compile.mesh.axis_subdivisions.empty() &&
            scenario.compile.mesh.lower == std::vector<double>({0.0, 0.0}) &&
            scenario.compile.mesh.upper == std::vector<double>({1.0, 1.0}) &&
            scenario.compile.mesh.centroid_splits == 0 &&
            scenario.compile.mesh.selection_seed == 0 &&
            scenario.compile.mesh.mesh_provenance ==
              "chapter-6.e6.5.1.structured-simplex-n131" &&
            scenario.compile.state_degree == 1 &&
            scenario.compile.execution ==
              nmopt::application::chapter6::ExecutionSelection::assembled &&
            scenario.compile.product ==
              nmopt::application::chapter6::ProductSelection::reduced_dto &&
            scenario.compile.owned_session &&
            scenario.compile.state_solve.maximum_iterations == 0 &&
            scenario.compile.adjoint_solve.maximum_iterations == 0 &&
            scenario.compile.control_metric_solve.maximum_iterations == 1000,
          "B1 resolution changed effective mesh or compile policy");

        require(
          scenario.solver.initial_control_value == 0.0 &&
            scenario.solver.objective_target_policy ==
              nmopt::application::chapter6::ObjectiveTargetPolicy::none &&
            scenario.solver.parameters.maximum_iterations == 5000 &&
            scenario.solver.parameters.maximum_line_search_trials == 6 &&
            scenario.solver.parameters.gradient_tolerance == 1.0e-30 &&
            scenario.solver.parameters.relative_gradient_tolerance == 1.0e-3 &&
            scenario.solver.parameters.objective_change_tolerance == 0.0 &&
            scenario.solver.parameters.step_tolerance == 0.0 &&
            scenario.solver.parameters.initial_step_length == 1.0 &&
            scenario.solver.parameters.armijo_fraction == 1.0e-5 &&
            scenario.solver.parameters.backtracking_factor == 0.7,
          "B1 resolution changed the common solver policy");

        if (method_id == "steepest-descent")
          require(
            scenario.solver.method ==
                nmopt::application::chapter6::ReducedMethod::steepest_descent &&
              scenario.solver.parameters.minimum_step_length == 0.2,
            "B1 steepest-descent resolution changed its policy");
        else
          require(
            scenario.solver.method ==
                nmopt::application::chapter6::ReducedMethod::
                  limited_memory_bfgs &&
              scenario.solver.parameters.minimum_step_length == 0.01 &&
              scenario.solver.limited_memory_bfgs.memory_size == 5 &&
              scenario.solver.limited_memory_bfgs.curvature_tolerance ==
                1.0e-14 &&
              scenario.solver.limited_memory_bfgs.
                  initial_inverse_hessian_scaling ==
                nmopt::solvers::LimitedMemoryBfgsInitialScaling::metric_inverse,
            "B1 L-BFGS resolution changed its policy");
      }
  }

  void
  test_b2_authoritative_and_forcing_sweep_resolution_is_characterized()
  {
    const auto file = read_parameter_file(find_file_from_current_or_parent(
      "parameters/chapter-6/b2/authoritative.prm"));
    const auto combinations = file.combinations();
    require(combinations.size() == 4,
            "B2 characterization lost its four-case matrix");

    const std::vector<std::string> expected_paths{
      "artifacts/wings-constant/artifact.kv",
      "artifacts/wings-parabolic/artifact.kv",
      "artifacts/full-constant/artifact.kv",
      "artifacts/full-parabolic/artifact.kv"};
    const auto b2_plan = nmopt::application::runner::make_run_set_plan(file);
    require(nmopt::application::runner::run_set_artifact_paths(
              b2_plan, b2_artifact_coordinate_components) ==
              expected_paths,
            "B2 characterization changed artifact coordinate paths");

    for (const auto &combination : combinations)
      {
        const auto region = combination_value(combination, "observation-region");
        const auto target_profile =
          combination_value(combination, "target-profile");
        const auto scenario =
          resolve_b2_scenario_for_characterization(file, combination);
        const auto &observation_region =
          combination_value(combination, "observation-region");
        const auto &target = nmopt::application::chapter6::b2_target_definition(
          scenario.problem.target_catalog);
        const auto expected_metadata_id =
          observation_region == "wings" &&
          target_profile == "constant"
            ? std::string("chapter-6.b2.graetz-flow")
            : std::string("chapter-6.b2.graetz-flow.") +
                  nmopt::application::chapter6::b2_case_name(
                    observation_region, target_profile);

        require(scenario.metadata.id == expected_metadata_id &&
                  scenario.metadata.recipe_id ==
                    nmopt::application::chapter6::b2_recipe_id &&
                  scenario.experiment.scenario_output_id ==
                    std::string("chapter-6.b2.graetz-flow.") +
                      nmopt::application::chapter6::b2_case_name(
                        observation_region, target_profile),
                "B2 resolution changed the benchmark identity");
        require(scenario.problem.observation_region == observation_region &&
                  nmopt::application::chapter6::b2_observation_region_name(
                    observation_region) == region &&
                  scenario.problem.target_profile == target_profile &&
                  ((target_profile == "constant" &&
                    target.id == "constant" &&
                    target.kind ==
                      nmopt::application::ScalarFunctionKind::constant &&
                    target.value == 2.0 && target.expression.empty()) ||
                   (target_profile == "parabolic" &&
                    target.id == "parabolic" &&
                    target.kind ==
                      nmopt::application::ScalarFunctionKind::expression &&
                    target.value == 0.0 &&
                    target.expression == "4.0*x1*(1.0-x1)")),
                "B2 resolution changed the selected target definition");
        require(scenario.problem.target_catalog.definitions.size() == 2 &&
                  std::all_of(
                    scenario.problem.target_catalog.definitions.begin(),
                    scenario.problem.target_catalog.definitions.end(),
                    [](const auto &definition) {
                      return definition.provenance ==
                             "chapter-6.e6.5.2.target";
                    }) &&
                  scenario.problem.target_catalog.selected_id == target.id &&
                  scenario.problem.data.desired_state_provenance ==
                    "chapter-6.e6.5.2.target",
                "B2 resolution changed target provenance");
        require(
          scenario.problem.recipe.control_discretisation ==
              nmopt::semantic::v1::NeumannControlDiscretisation::facewise_constant &&
            !scenario.problem.recipe.with_facewise_box &&
            scenario.problem.recipe.observed_material_id == 1 &&
            scenario.problem.fixed_dirichlet_data.id == "fixed-temperature" &&
            scenario.problem.fixed_dirichlet_data.kind ==
              nmopt::application::ScalarFunctionKind::constant &&
            scenario.problem.fixed_dirichlet_data.value == 1.0 &&
            scenario.problem.fixed_dirichlet_data.expression.empty() &&
            scenario.problem.fixed_dirichlet_data.provenance ==
              "chapter-6.e6.5.2.fixed-temperature" &&
            scenario.problem.conservative_transport.id == "graetz" &&
            scenario.problem.conservative_transport.expression ==
              "1.5*x1*(1-x1); 0.0" &&
            scenario.problem.conservative_transport.provenance ==
              "chapter-6.e6.5.2.graetz-transport" &&
            scenario.problem.forcing.id == "zero" &&
            scenario.problem.forcing.kind ==
              nmopt::application::ScalarFunctionKind::zero &&
            scenario.problem.forcing.value == 0.0 &&
            scenario.problem.forcing.expression.empty() &&
            scenario.problem.forcing.provenance ==
              "chapter-6.e6.5.2.zero-forcing" &&
            scenario.problem.data.forcing_provenance ==
              "chapter-6.e6.5.2.zero-forcing" &&
            scenario.problem.transport_boundary_form ==
              nmopt::semantic::v1::TransportBoundaryForm::
                ordinary_normal_minus_transport,
          "B2 resolution changed effective forcing, control, or boundary data");
        require(
          scenario.compile.mesh.dimension == 2 &&
            scenario.compile.mesh.generation ==
              nmopt::application::chapter6::MeshGeneration::framework_native &&
            scenario.compile.mesh.refinement == 7 &&
            scenario.compile.mesh.lower == std::vector<double>({0.0, 0.0}) &&
            scenario.compile.mesh.upper == std::vector<double>({4.0, 1.0}) &&
            scenario.compile.mesh.mesh_provenance ==
              "chapter-6.e6.5.2.framework-native-rectangle-r7" &&
            scenario.compile.volume_observation.has_value() &&
            scenario.compile.volume_observation->quadrature_order == 3 &&
            scenario.compile.volume_observation->target_realisation ==
              nmopt::application::chapter6::VolumeObservationTargetRealisation::
                analytic_quadrature &&
            scenario.compile.execution ==
              nmopt::application::chapter6::ExecutionSelection::assembled &&
            scenario.compile.product ==
              nmopt::application::chapter6::ProductSelection::reduced_dto,
          "B2 resolution changed effective mesh or compile policy");
        require(
          scenario.solver.method ==
              nmopt::application::chapter6::ReducedMethod::bfgs &&
            scenario.solver.globalization ==
              nmopt::application::chapter6::ReducedGlobalization::armijo &&
            scenario.solver.initial_control_value == 0.0 &&
            scenario.solver.parameters.maximum_iterations == 100 &&
            scenario.solver.parameters.maximum_line_search_trials == 20 &&
            scenario.solver.parameters.gradient_tolerance == 1.0e-8 &&
            scenario.solver.parameters.initial_step_length == 1.0 &&
            scenario.solver.parameters.armijo_fraction == 1.0e-4 &&
            scenario.solver.parameters.backtracking_factor == 0.5 &&
            scenario.solver.parameters.minimum_step_length == 0.0 &&
            scenario.solver.full_bfgs.curvature_tolerance == 1.0e-14,
          "B2 resolution changed the frozen solver policy");
      }

    const auto forcing_sweep = read_parameter_file(find_file_from_current_or_parent(
      "parameters/chapter-6/b2/development/forcing-sweep.prm"));
    const auto forcing_combinations = forcing_sweep.combinations();
    const std::vector<std::string> forcing_paths{
      "artifacts/regularisation-1e-3/forcing-zero/observation-region-wings/"
      "target-profile-constant/artifact.kv",
      "artifacts/regularisation-1e-3/forcing-constant-one/"
      "observation-region-wings/target-profile-constant/artifact.kv",
      "artifacts/regularisation-1e-3/forcing-constant-two/"
      "observation-region-wings/target-profile-constant/artifact.kv"};
    require(forcing_combinations.size() == 3 &&
              nmopt::application::runner::run_set_artifact_paths(
                nmopt::application::runner::make_run_set_plan(forcing_sweep),
                b2_artifact_coordinate_components) ==
                forcing_paths,
            "B2 forcing-sweep characterization changed its matrix paths");

    for (const auto &combination : forcing_combinations)
      {
        const auto forcing_id = combination_value(combination, "forcing");
        const auto scenario =
          resolve_b2_scenario_for_characterization(forcing_sweep, combination);
        const auto expected_kind = forcing_id == "zero"
                                     ? nmopt::application::ScalarFunctionKind::zero
                                     : nmopt::application::ScalarFunctionKind::constant;
        const auto expected_value = forcing_id == "zero"
                                      ? 0.0
                                      : forcing_id == "constant-one" ? 1.0 : 2.0;
        const auto expected_provenance =
          std::string("development.b2.") + forcing_id + "-forcing";
        require(scenario.problem.forcing.id == forcing_id &&
                  scenario.problem.forcing.kind == expected_kind &&
                  scenario.problem.forcing.value == expected_value &&
                  scenario.problem.forcing.expression.empty() &&
                  scenario.problem.forcing.provenance == expected_provenance &&
                  scenario.problem.data.forcing_provenance ==
                    expected_provenance &&
                  scenario.problem.observation_region == "wings" &&
                  scenario.problem.target_profile == "constant" &&
                  scenario.problem.data.regularisation_weight == 1.0e-3 &&
                  scenario.compile.mesh.refinement == 6 &&
                  scenario.compile.mesh.mesh_provenance ==
                    "chapter-6.e6.5.2.framework-native-rectangle-r6",
                "B2 forcing-sweep resolution changed effective forcing data");
      }
  }

  void
  test_b2_output_selection_is_explicit()
  {
    const auto source = find_file_from_current_or_parent(
      "parameters/chapter-6/b2/authoritative.prm");
    const auto combination = read_parameter_file(source).combinations().front();

    auto unsupported_output_selection = read_parameter_file(source);
    unsupported_output_selection.values["Output/selected fields"] =
      "state, control";
    require_invalid_argument(
      [&] {
        (void)resolve_b2_scenario_for_characterization(
          unsupported_output_selection, combination);
      },
      "B2 accepted an output selection that its adapter does not execute");
  }

  void
  test_runner_execution_registrations_bind_callbacks()
  {
    const auto *b1 = find_benchmark_execution_registration("b1");
    const auto *b2 = find_benchmark_execution_registration("b2");
    require(b1 != nullptr && b2 != nullptr,
            "B1 and B2 should have executable runner registrations");
    require(b1->metadata != nullptr && b2->metadata != nullptr &&
              b1->artifact_planner && b1->execute && b2->artifact_planner &&
              b2->execute,
            "runner registrations should bind artifact and execution callbacks");

    const auto b1_file = read_parameter_file(find_file_from_current_or_parent(
      "parameters/chapter-6/b1/authoritative.prm"));
    const auto b1_plan = nmopt::application::runner::make_run_set_plan(b1_file);
    const auto b1_artifacts = b1->artifact_planner(b1_plan);
    require(b1_artifacts == nmopt::application::runner::run_set_artifact_paths(
                              b1_plan, b1_artifact_coordinate_components),
            "B1 registration should use its registered artifact planner");

    const auto b2_file = read_parameter_file(find_file_from_current_or_parent(
      "parameters/chapter-6/b2/authoritative.prm"));
    const auto b2_plan = nmopt::application::runner::make_run_set_plan(b2_file);
    const auto b2_artifacts = b2->artifact_planner(b2_plan);
    require(b2_artifacts == nmopt::application::runner::run_set_artifact_paths(
                              b2_plan, b2_artifact_coordinate_components),
            "B2 registration should use its registered artifact planner");
  }

  nmopt::application::runner::ResolvedRunConfiguration
  characterization_run_configuration(
    const nmopt::application::runner::ParameterFile &file,
    const std::filesystem::path                    &root,
    const std::string                               &benchmark)
  {
    using nmopt::application::runner::ResolvedRunConfiguration;
    using nmopt::application::runner::RunKind;
    ResolvedRunConfiguration configuration{root,
                                           root / "artifact",
                                           benchmark,
                                           "debug-dealii",
                                           "characterization-revision",
                                           RunKind::development,
                                           0};
    configuration.parameter_file = file.path;
    configuration.parameter_hash = file.content_hash;
    configuration.plotting_profile_file = find_file_from_current_or_parent(
      benchmark == "b1" ? "parameters/plotting/chapter-6-b1.json"
                         : "parameters/plotting/chapter-6-b2.json");
    configuration.plotting_profile_hash =
      nmopt::application::runner::parameter_file_hash(
      configuration.plotting_profile_file);
    return configuration;
  }

  void
  test_parameter_resolution_retains_artifact_fields_and_manifest_provenance()
  {
    const auto b1_file = read_parameter_file(find_file_from_current_or_parent(
      "parameters/chapter-6/b1/authoritative.prm"));
    const auto b1_combination = b1_file.combinations().front();
    const auto b1_beta = parse_number_text(
      combination_value(b1_combination, "regularisation"),
      "Matrix/regularisation");
    auto b1_scenario = resolve_b1_scenario_for_characterization(
      b1_file, b1_combination);
    b1_scenario.compile.mesh.generation =
      nmopt::application::chapter6::MeshGeneration::framework_native;
    b1_scenario.compile.mesh.subdivisions = 0;
    b1_scenario.compile.mesh.refinement = 1;
    b1_scenario.compile.mesh.mesh_provenance =
      "test.parameter-resolution.b1-small-mesh";
    b1_scenario.solver.parameters.maximum_iterations = 2;
    b1_scenario.solver.parameters.gradient_tolerance = 1.0e-3;
    b1_scenario.experiment.harness.measure_timings = false;

    const auto b1_native_output = std::filesystem::temp_directory_path() /
                                  "nmopt-parameter-resolution-b1";
    std::filesystem::remove_all(b1_native_output);
    nmopt::application::chapter6::dealii::B1SelectedDataT<2> b1_data(
      b1_scenario.problem.forcing);
    const auto b1_runtime =
      nmopt::application::chapter6::dealii::make_b1_runtime_data(
        b1_scenario, b1_data);
    const auto b1_session =
      nmopt::application::chapter6::dealii::make_b1_compilation_session<2>(
        b1_scenario);
    nmopt::application::chapter6::dealii::B1ReducedExecutionAdapterT<2>
      b1_execute{
        b1_beta,
        b1_runtime,
        b1_session,
        {"characterization-revision",
         "debug-dealii",
         "test-compiler",
         "test-version",
         "test-standard-library",
         "test-os",
         "test-architecture",
         "test-hardware"},
        b1_native_output};
    using B1Runner = nmopt::application::benchmark::HeadlessBenchmarkRunnerT<
      nmopt::application::chapter6::B1Scenario>;
    const auto b1_configuration = characterization_run_configuration(
      b1_file, b1_native_output, "b1");
    const auto b1_result = B1Runner(b1_scenario).run(
      [](const auto &parameters) {
        return nmopt::application::chapter6::make_b1_problem_spec(parameters);
      },
      [&](const auto &specification, const auto &scenario) {
        auto evidence = b1_execute(specification, scenario);
        add_parameter_artifact_fields(
          evidence, b1_configuration, b1_combination);
        return evidence;
      });
    const std::vector<std::string> expected_b1_fields{
      "objective_history",
      "gradient_norm_history",
      "step_length_history",
      "objective_change_history",
      "solve_counts",
      "hessian_evidence",
      "state",
      "control",
      "adjoint",
      "target",
      "forcing",
      "negative_adjoint"};
    require(b1_result.artifact.selected_fields() == expected_b1_fields,
            "B1 parameter resolution changed selected artifact fields");
    require(b1_result.document.find(
              "provenance.parameters_file=" + b1_file.path.string() + "\n") !=
              std::string::npos &&
              b1_result.document.find(
                "provenance.parameters_hash=" + b1_file.content_hash + "\n") !=
                std::string::npos &&
              b1_result.document.find("parameters.method=steepest-descent\n") !=
                std::string::npos &&
              b1_result.document.find("parameters.regularisation=1e-1\n") !=
                std::string::npos,
            "B1 artifact lost parameter-file provenance");
    require_binding_provenance(
      b1_result.artifact.envelope().compilation_manifest(),
      "forcing",
      "chapter-6.e6.5.1.source-oriented-constant-half-forcing",
      "B1 manifest lost forcing provenance");
    require_binding_provenance(
      b1_result.artifact.envelope().compilation_manifest(),
      "desired_state",
      "chapter-6.e6.5.1.desired-state",
      "B1 manifest lost desired-state provenance");

    const auto b2_file = read_parameter_file(find_file_from_current_or_parent(
      "parameters/chapter-6/b2/authoritative.prm"));
    auto b2_customized_file = b2_file;
    b2_customized_file.values["Solver/initial independent control value"] =
      "0.25";
    b2_customized_file.values[
      "Solver/method policy bfgs/curvature tolerance"] =
      "1e-8";
    const auto b2_combination = b2_file.combinations().front();
    auto b2_scenario = resolve_b2_scenario_for_characterization(
      b2_customized_file, b2_combination);
    require(b2_scenario.solver.initial_control_value == 0.25 &&
              b2_scenario.solver.full_bfgs.curvature_tolerance == 1.0e-8,
            "changed B2 solver numerics were not bound from the parameter file");
    b2_scenario.compile.mesh.refinement = 1;
    b2_scenario.compile.mesh.mesh_provenance =
      "test.parameter-resolution.b2-small-mesh";
    b2_scenario.solver.parameters.maximum_iterations = 2;
    b2_scenario.solver.parameters.gradient_tolerance = 1.0e-3;
    b2_scenario.experiment.harness.measure_timings = false;

    const auto b2_native_output = std::filesystem::temp_directory_path() /
                                  "nmopt-parameter-resolution-b2";
    std::filesystem::remove_all(b2_native_output);
    nmopt::application::chapter6::dealii::B2ManufacturedDataT<2> b2_data{
      nmopt::application::selected_scalar_function_definition(
        b2_scenario.problem.observation_region_catalog),
      b2_scenario.problem.fixed_dirichlet_data,
      b2_scenario.problem.forcing,
      nmopt::application::chapter6::b2_target_definition(
        b2_scenario.problem.target_catalog),
      b2_scenario.problem.conservative_transport};
    const auto b2_runtime =
      nmopt::application::chapter6::dealii::make_b2_manufactured_runtime_data<2>(
        b2_scenario, b2_data);
    const auto b2_session =
      nmopt::application::chapter6::dealii::make_b2_compilation_session<2>(
        b2_scenario);
    nmopt::application::chapter6::dealii::B2ReducedExecutionAdapterT<2>
      b2_execute{
        b2_runtime,
        b2_session,
        {"characterization-revision",
         "debug-dealii",
         "test-compiler",
         "test-version",
         "test-standard-library",
         "test-os",
         "test-architecture",
         "test-hardware"},
        b2_native_output};
    using B2Runner = nmopt::application::benchmark::HeadlessBenchmarkRunnerT<
      nmopt::application::chapter6::B2Scenario>;
    const auto b2_configuration = characterization_run_configuration(
      b2_file, b2_native_output, "b2");
    const auto b2_result = B2Runner(b2_scenario).run(
      [](const auto &parameters) {
        return nmopt::application::chapter6::make_b2_problem_spec(parameters);
      },
      [&](const auto &specification, const auto &scenario) {
        auto evidence = b2_execute(specification, scenario);
        add_b2_artifact_fields(evidence,
                               scenario,
                               "characterization-revision",
                               nmopt::application::runner::RunKind::development);
        add_parameter_artifact_fields(
          evidence, b2_configuration, b2_combination);
        return evidence;
      });
    const std::vector<std::string> expected_b2_fields{
      "objective_history",
      "gradient_norm_history",
      "step_length_history",
      "objective_change_history",
      "solve_counts",
      "state",
      "state_uncontrolled",
      "control",
      "adjoint",
      "negative_adjoint",
      "target",
      "forcing",
      "observation_region"};
    require(b2_result.artifact.selected_fields() == expected_b2_fields,
            "B2 parameter resolution changed selected artifact fields");
    require(b2_result.document.find(
              "provenance.parameters_file=" + b2_file.path.string() + "\n") !=
              std::string::npos &&
              b2_result.document.find(
                "provenance.parameters_hash=" + b2_file.content_hash + "\n") !=
                std::string::npos &&
              b2_result.document.find(
                "parameters.observation-region=wings\n") != std::string::npos &&
              b2_result.document.find("parameters.target-profile=constant\n") !=
                std::string::npos,
            "B2 artifact lost parameter-file provenance");
    require(
      b2_result.document.find("solver.initial_control_value=0.25\n") !=
          std::string::npos &&
        b2_result.document.find("solver.full_bfgs_curvature_tolerance=1e-08\n") !=
          std::string::npos &&
        b2_result.document.find(
          "b2.initial_control_regularisation_objective=0\n") ==
          std::string::npos,
      "changed B2 solver numerics did not affect execution evidence");
    require_binding_provenance(
      b2_result.artifact.envelope().compilation_manifest(),
      "forcing",
      "chapter-6.e6.5.2.zero-forcing",
      "B2 manifest lost forcing provenance");
    require_binding_provenance(
      b2_result.artifact.envelope().compilation_manifest(),
      "fixed_dirichlet_data",
      "chapter-6.e6.5.2.fixed-temperature",
      "B2 manifest lost fixed-temperature provenance");
    std::filesystem::remove_all(b1_native_output);
    std::filesystem::remove_all(b2_native_output);
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
  test_b2_volume_observation_is_selected()
  {
    const auto analytic = read_parameter_file(find_file_from_current_or_parent(
      "parameters/chapter-6/b2/authoritative.prm"));
    const auto analytic_options =
      nmopt::application::runner::b2_volume_observation_options(analytic);
    require(
      analytic_options.quadrature_order == 3 &&
        analytic_options.target_realisation ==
          nmopt::application::chapter6::
            VolumeObservationTargetRealisation::analytic_quadrature,
      "B2 parameter parsing lost the frozen volume-observation policy");

    auto interpolated = analytic;
    interpolated.values["Compile/volume observation quadrature order"] = "2";
    interpolated.values["Compile/volume observation target realisation"] =
      "state-fe-interpolation";
    const auto interpolated_options =
      nmopt::application::runner::b2_volume_observation_options(interpolated);
    require(
      interpolated_options.quadrature_order == 2 &&
        interpolated_options.target_realisation ==
          nmopt::application::chapter6::
            VolumeObservationTargetRealisation::state_fe_interpolation,
      "B2 parameter parsing lost the interpolated observation policy");

    auto zero_order = analytic;
    zero_order.values["Compile/volume observation quadrature order"] = "0";
    require_invalid_argument(
      [&] {
        (void)nmopt::application::runner::b2_volume_observation_options(
          zero_order);
      },
      "B2 parameter parsing accepted a zero observation quadrature order");

    auto unknown = analytic;
    unknown.values["Compile/volume observation target realisation"] =
      "nodal-sampling";
    require_invalid_argument(
      [&] {
        (void)nmopt::application::runner::b2_volume_observation_options(
          unknown);
      },
      "B2 parameter parsing accepted an unknown observation target realisation");
  }

  void
  test_b2_target_catalog_is_data_driven()
  {
    auto file = read_parameter_file(find_file_from_current_or_parent(
      "parameters/chapter-6/b2/authoritative.prm"));
    const auto defaults =
      nmopt::application::runner::b2_target_catalog(file);
    require(defaults.selected_id == "constant" &&
              defaults.definitions.size() == 2 &&
              defaults.definitions[0].id == "constant" &&
              defaults.definitions[0].kind ==
                nmopt::application::ScalarFunctionKind::constant &&
              std::abs(defaults.definitions[0].value - 2.0) < 1.0e-15 &&
              defaults.definitions[0].expression.empty() &&
              defaults.definitions[1].id == "parabolic" &&
              defaults.definitions[1].kind ==
                nmopt::application::ScalarFunctionKind::expression &&
              defaults.definitions[1].expression == "4.0*x1*(1.0-x1)" &&
              defaults.definitions[1].value == 0.0,
            "B2 parameter parsing lost its target catalog");

    file.values["Functions/target definitions/constant/value"] = "20.0";
    file.values["Functions/target definitions/parabolic/expression"] =
      "6.0*x1*(1.0-x1)";
    const auto candidate =
      nmopt::application::runner::b2_target_catalog(file, "parabolic");
    require(candidate.selected_id == "parabolic" &&
              candidate.definitions[0].id == "constant" &&
              std::abs(candidate.definitions[0].value - 20.0) < 1.0e-15 &&
              candidate.definitions[1].id == "parabolic" &&
              candidate.definitions[1].expression == "6.0*x1*(1.0-x1)" &&
              nmopt::application::chapter6::b2_target_definition(candidate).id ==
                candidate.selected_id,
            "B2 target catalog was hardcoded instead of parsed");

    auto duplicate = candidate;
    duplicate.definitions[1].id = duplicate.definitions[0].id;
    require_invalid_argument(
      [&] {
        nmopt::application::validate_scalar_function_catalog(
          duplicate, "test B2 target catalog");
      },
      "B2 target catalog accepted duplicate definition IDs");

    auto missing_selection = candidate;
    missing_selection.selected_id = "not-registered";
    require_invalid_argument(
      [&] {
        nmopt::application::validate_scalar_function_catalog(
          missing_selection, "test B2 target catalog");
      },
      "B2 target catalog accepted an unregistered selected definition");

    file.values["Functions/target definitions/constant/value"] = "nan";
    require_invalid_argument(
      [&] {
        (void)nmopt::application::runner::b2_target_catalog(file);
      },
      "B2 parameter parsing accepted a non-finite target value");
  }

  void
  test_b2_scalar_definitions_are_discovered_from_prm_ids()
  {
    const auto file = read_b2_scalar_discovery_parameter_file();
    require(
      file.value("Functions/forcing definitions/spatial-candidate/kind") ==
          "expression" &&
        file.value(
          "Functions/forcing definitions/spatial-candidate/expression") ==
          "0.25 + sin(pi*x0)*sin(pi*x1)" &&
        file.value("Functions/target definitions/candidate-target/kind") ==
          "expression",
      "native scalar subsections were not declared from matrix IDs");

    const auto combinations = file.combinations();
    require(combinations.size() == 2,
            "scalar discovery fixture lost its forcing matrix");
    const auto combination_it = std::find_if(
      combinations.begin(), combinations.end(), [](const auto &candidate) {
        return combination_value(candidate, "forcing") ==
               "spatial-candidate";
      });
    require(combination_it != combinations.end(),
            "scalar discovery fixture lost its selected forcing combination");
    const auto &combination = *combination_it;

    auto missing_selected_definition = file;
    missing_selected_definition.values.erase(
      "Functions/forcing definitions/spatial-candidate/kind");
    require_invalid_argument(
      [&] {
        (void)nmopt::application::runner::
          parameter_scalar_function_definition_from_catalog(
            missing_selected_definition,
            "Functions/forcing definitions/",
            "spatial-candidate");
      },
      "a selected matrix scalar definition may not be missing");

    const auto target_catalog =
      nmopt::application::runner::b2_target_catalog(file, "candidate-target");
    auto scenario =
      nmopt::application::chapter6::make_b2_scenario_with_target_catalog(
        nmopt::application::chapter6::B2ObservationRegion::wings,
        "candidate-target",
        target_catalog);
    bind_b2_scenario(
      scenario,
      file,
      combination,
      "chapter-6.b2.graetz-flow.wings-candidate-target");

    const auto &target =
      nmopt::application::chapter6::b2_target_definition(
        scenario.problem.target_catalog);
    require(scenario.problem.forcing.id == "spatial-candidate" &&
              scenario.problem.forcing.kind ==
                nmopt::application::ScalarFunctionKind::expression &&
              scenario.problem.forcing.expression ==
                "0.25 + sin(pi*x0)*sin(pi*x1)" &&
              scenario.problem.forcing.provenance ==
                "test.spatial-candidate-forcing" &&
              target.id == "candidate-target" &&
              target.kind ==
                nmopt::application::ScalarFunctionKind::expression &&
              target.expression == "x0 + 2.0*x1" &&
              target.provenance == "test.target-candidate" &&
              scenario.metadata.id ==
                "chapter-6.b2.graetz-flow.wings-candidate-target",
            "B2 scalar definitions were not resolved from native subsections");
  }

  void
  test_b2_observation_and_boundary_numerics_are_prm_owned()
  {
    const auto default_file = read_b2_scalar_discovery_parameter_file();
    const auto customized_file = read_b2_scalar_discovery_parameter_file(
      "x0 > 1.0 ? (x1 < 0.6 ? 1 : 0) : 0", 7, 8, 9);

    auto default_scenario =
      nmopt::application::chapter6::make_b2_scenario_with_target_catalog(
        "wings",
        "candidate-target",
        nmopt::application::runner::b2_target_catalog(
          default_file, "candidate-target"));
    bind_b2_scenario(default_scenario,
                     default_file,
                     default_file.combinations().front(),
                     "chapter-6.b2.graetz-flow.wings-candidate-target");
    auto customized_scenario =
      nmopt::application::chapter6::make_b2_scenario_with_target_catalog(
        "wings",
        "candidate-target",
        nmopt::application::runner::b2_target_catalog(
          customized_file, "candidate-target"));
    bind_b2_scenario(customized_scenario,
                     customized_file,
                     customized_file.combinations().front(),
                     "chapter-6.b2.graetz-flow.wings-candidate-target");
    const auto default_session =
      nmopt::application::chapter6::dealii::make_b2_compilation_session<2>(
        default_scenario);
    const auto customized_session =
      nmopt::application::chapter6::dealii::make_b2_compilation_session<2>(
        customized_scenario);

    const auto observation_summary = [](const auto &session, const auto &scenario) {
      std::size_t observed_cells = 0;
      double      observed_measure = 0.0;
      for (const auto &cell : session->triangulation().active_cell_iterators())
        if (cell->material_id() ==
            scenario.problem.recipe.observed_material_id)
          {
            ++observed_cells;
            observed_measure += cell->measure();
          }
      return std::pair<std::size_t, double>{observed_cells, observed_measure};
    };
    const auto default_observation =
      observation_summary(default_session, default_scenario);
    const auto customized_observation =
      observation_summary(customized_session, customized_scenario);
    require(default_observation.first == 2 &&
              std::abs(default_observation.second - 2.0) < 1.0e-12 &&
              customized_observation.first == 1 &&
              std::abs(customized_observation.second - 1.0) < 1.0e-12,
            "B2 observation classification did not follow the edited .prm expression");

    require(customized_scenario.problem.boundary.fixed_id == 7 &&
              customized_scenario.problem.boundary.control_id == 8 &&
              customized_scenario.problem.boundary.outflow_id == 9,
            "B2 boundary IDs were not retained from the edited .prm file");
    std::array<std::size_t, 3> boundary_face_counts{};
    for (const auto &cell :
         customized_session->triangulation().active_cell_iterators())
      for (unsigned int face = 0; face < cell->n_faces(); ++face)
        if (cell->face(face)->at_boundary())
          {
            const auto boundary_id = cell->face(face)->boundary_id();
            if (boundary_id == customized_scenario.problem.boundary.fixed_id)
              ++boundary_face_counts[0];
            else if (boundary_id ==
                     customized_scenario.problem.boundary.control_id)
              ++boundary_face_counts[1];
            else if (boundary_id ==
                     customized_scenario.problem.boundary.outflow_id)
              ++boundary_face_counts[2];
            else
              throw std::runtime_error(
                "B2 customized mesh has an unclassified boundary face");
          }
    require(boundary_face_counts == std::array<std::size_t, 3>{{2, 4, 2}},
            "B2 customized boundary geometry did not retain the .prm IDs");
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
  test_run_set_plan_is_generic_and_structured()
  {
    const auto b1 = read_parameter_file(find_file_from_current_or_parent(
      "parameters/chapter-6/b1/authoritative.prm"));
    const auto b1_plan =
      nmopt::application::runner::make_run_set_plan(b1);
    require(b1_plan.benchmark_id == "chapter-6.b1.distributed-laplace" &&
              b1_plan.matrix_axes.size() == 2 &&
              b1_plan.matrix_axes[0].id == "method" &&
              b1_plan.matrix_axes[1].id == "regularisation" &&
              b1_plan.resolved_combinations.size() == 7 &&
              b1_plan.excluded_combinations.size() == 1,
            "run-set planning changed the authoritative B1 matrix resolution");
    require(b1_plan.parameter_provenance.file == b1.path &&
              b1_plan.parameter_provenance.content_hash == b1.content_hash,
            "run-set planning did not retain parameter-file provenance");
    require(b1_plan.comparison.rows == "method" &&
              b1_plan.comparison.columns == "regularisation" &&
              b1_plan.comparison.group_by == "none",
            "run-set planning did not retain comparison coordinates");
    require(
      b1_plan.resolved_combinations.front().artifact_coordinates.size() == 2 &&
        b1_plan.resolved_combinations.front().artifact_coordinates[0].axis ==
          "method" &&
        b1_plan.resolved_combinations.front().artifact_coordinates[0].value ==
          "steepest-descent" &&
        b1_plan.resolved_combinations.front().artifact_coordinates[1].axis ==
          "regularisation",
      "run-set planning did not derive ordered artifact coordinates");

    nmopt::application::runner::ParameterFile future;
    future.path = "parameters/test/future.prm";
    future.content_hash = "fnv1a64:future";
    future.values = {{"Benchmark/id", "test.future"},
                     {"Postprocessing/comparison rows", "beta"},
                     {"Postprocessing/comparison columns", "profile"},
                     {"Postprocessing/comparison group by", "none"}};
    future.matrix = {{"beta", {"1", "2"}},
                     {"profile", {"left", "right"}}};
    future.selection = {{"beta", "1,2"}};
    nmopt::application::runner::ParameterCombination excluded;
    excluded.values = {{"beta", "2"}, {"profile", "left"}};
    future.excluded_combinations.push_back(excluded);

    const auto future_plan = nmopt::application::runner::make_run_set_plan(
      future, {{"profile", "right"}});
    require(future_plan.benchmark_id == "test.future" &&
              future_plan.matrix_axes[1].id == "profile" &&
              future_plan.selection.at("profile") == "right" &&
              future_plan.resolved_combinations.size() == 2,
            "run-set planning should accept a novel registered axis");
    require(
      future_plan.resolved_combinations[1].artifact_coordinates[0].value ==
          "2" &&
        future_plan.resolved_combinations[1].artifact_coordinates[1].value ==
          "right",
      "run-set planning should keep novel-axis coordinates independent of order");
    require(
      nmopt::application::runner::run_set_artifact_paths(future_plan) ==
        std::vector<std::string>{
          "artifacts/beta-1/profile-right/artifact.kv",
          "artifacts/beta-2/profile-right/artifact.kv"},
      "generic run-set artifact planning should consume novel axes");
  }

  void
  test_scalar_function_definitions_are_data_driven()
  {
    nmopt::application::runner::ParameterFile file;
    file.values = {{"Functions/forcing", "candidate-045"},
                   {"Functions/forcing/kind",
                    "constant"},
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
    file.values["Functions/forcing/kind"] =
      "expression";
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
          parameter_scalar_function_definition_from_selector(
            file, "Functions/forcing");
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
          parameter_scalar_function_definition_from_selector(
            file, "Functions/forcing");
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
        {"parameter_schema_registry_accounts_for_all_entries",
         "nmopt.parameter_files.parameter_schema_registry_accounts_for_all_entries",
         {"backend", "dealii", "application", "runner", "contract"},
         30,
         test_parameter_schema_registry_accounts_for_all_entries},
        {"checked_in_families_expand_and_filter",
         "nmopt.parameter_files.checked_in_families_expand_and_filter",
         {"backend", "dealii", "application", "runner", "contract"},
         30,
         test_checked_in_families_expand_and_filter},
        {"b1_authoritative_resolution_is_characterized",
         "nmopt.parameter_files.b1_authoritative_resolution_is_characterized",
         {"backend", "dealii", "application", "benchmark", "runner", "contract"},
         30,
         test_b1_authoritative_resolution_is_characterized},
        {"b2_authoritative_and_forcing_sweep_resolution_is_characterized",
         "nmopt.parameter_files.b2_authoritative_and_forcing_sweep_resolution_is_characterized",
         {"backend", "dealii", "application", "benchmark", "runner", "contract"},
         30,
         test_b2_authoritative_and_forcing_sweep_resolution_is_characterized},
        {"b2_output_selection_is_explicit",
         "nmopt.parameter_files.b2_output_selection_is_explicit",
         {"backend", "dealii", "application", "benchmark", "runner", "contract", "negative"},
         30,
         test_b2_output_selection_is_explicit},
        {"runner_execution_registrations_bind_callbacks",
         "nmopt.parameter_files.runner_execution_registrations_bind_callbacks",
         {"backend", "dealii", "application", "benchmark", "runner", "contract"},
         30,
         test_runner_execution_registrations_bind_callbacks},
        {"parameter_resolution_retains_artifact_fields_and_manifest_provenance",
         "nmopt.parameter_files.parameter_resolution_retains_artifact_fields_and_manifest_provenance",
         {"backend", "dealii", "application", "benchmark", "runner", "contract"},
         120,
         test_parameter_resolution_retains_artifact_fields_and_manifest_provenance},
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
        {"mesh_bounds_are_data_driven_and_validated",
         "nmopt.parameter_files.mesh_bounds_are_data_driven_and_validated",
         {"backend", "dealii", "application", "runner", "contract", "negative"},
         30,
         test_mesh_bounds_are_data_driven_and_validated},
        {"b1_beta_coordinates_are_text_stable",
         "nmopt.parameter_files.b1_beta_coordinates_are_text_stable",
         {"backend", "dealii", "application", "runner", "contract"},
         30,
         test_b1_beta_coordinates_are_text_stable},
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
        {"b2_volume_observation_is_selected",
         "nmopt.parameter_files.b2_volume_observation_is_selected",
         {"backend", "dealii", "application", "runner", "contract", "negative"},
         30,
         test_b2_volume_observation_is_selected},
        {"b2_target_catalog_is_data_driven",
         "nmopt.parameter_files.b2_target_catalog_is_data_driven",
         {"backend", "dealii", "application", "runner", "contract", "negative"},
         30,
         test_b2_target_catalog_is_data_driven},
        {"b2_scalar_definitions_are_discovered_from_prm_ids",
         "nmopt.parameter_files.b2_scalar_definitions_are_discovered_from_prm_ids",
         {"backend", "dealii", "application", "benchmark", "runner", "contract"},
         30,
         test_b2_scalar_definitions_are_discovered_from_prm_ids},
        {"b2_observation_and_boundary_numerics_are_prm_owned",
         "nmopt.parameter_files.b2_observation_and_boundary_numerics_are_prm_owned",
         {"backend", "dealii", "application", "benchmark", "runner", "contract"},
         30,
         test_b2_observation_and_boundary_numerics_are_prm_owned},
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
        {"run_set_plan_is_generic_and_structured",
         "nmopt.parameter_files.run_set_plan_is_generic_and_structured",
         {"backend", "dealii", "application", "runner", "contract"},
         30,
         test_run_set_plan_is_generic_and_structured},
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
