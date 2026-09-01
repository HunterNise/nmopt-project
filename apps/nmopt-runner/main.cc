#include "nmopt/application/application.hpp"
#include "nmopt/application/dealii/chapter6_b1.hpp"
#include "nmopt/application/dealii/chapter6_b2.hpp"
#include "parameter_files.hpp"
#include "runner.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <locale>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef NMOPT_COMPILED_BUILD_PROFILE
#  define NMOPT_COMPILED_BUILD_PROFILE "unknown"
#endif

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

  std::string
  join_unsigned_ids(const std::vector<unsigned int> &ids)
  {
    std::ostringstream output;
    for (std::size_t index = 0; index < ids.size(); ++index)
      {
        if (index != 0)
          output << ',';
        output << ids[index];
      }
    return output.str();
  }

  void
  print_usage(std::ostream &output)
  {
    output << "Usage: nmopt_runner --list [--output DIRECTORY]\n"
           << "       nmopt_runner --parameter-file FILE --framework-revision REV"
              " [--output DIRECTORY] [--run-slot SLOT] [--select AXIS=VALUE]"
              " [--refinement N]\n"
           << "       nmopt_runner --benchmark b1 --framework-revision REV"
              " [--output DIRECTORY] [--run-kind KIND] [--run-slot SLOT]"
              " [--refinement N]\n"
           << "       nmopt_runner --benchmark b2 --framework-revision REV"
              " [--output DIRECTORY] [--run-kind KIND] [--run-slot SLOT]"
              " [--refinement N]\n"
           << "       nmopt_runner --help\n"
           << "\n"
           << "--list             list registered Chapter 5/6 application entries\n"
           << "--benchmark b1     run the frozen B1 matrix (eight artifacts)\n"
           << "--benchmark b2     run the frozen B2 case batch (four artifacts)\n"
           << "--parameter-file FILE\n"
           << "                   load a versioned experiment family\n"
           << "--select AXIS=VALUE\n"
           << "                   filter one declared matrix axis (repeatable)\n"
           << "--output DIRECTORY set the generated run-set root (default: runs)\n"
           << "--run-kind KIND    use reproduction or development policy\n"
           << "--run-slot SLOT    use a named development slot instead of auto-allocation\n"
           << "--framework-revision REV\n"
           << "                   record the framework revision in each artifact\n"
           << "--refinement N     optional framework-native mesh override; otherwise\n"
           << "                   use the selected benchmark configuration\n"
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

  const char *
  b1_method_slug(const nmopt::application::chapter6::ReducedMethod method)
  {
    switch (method)
      {
        case nmopt::application::chapter6::ReducedMethod::steepest_descent:
          return "steepest-descent";
        case nmopt::application::chapter6::ReducedMethod::limited_memory_bfgs:
          return "l-bfgs";
        case nmopt::application::chapter6::ReducedMethod::bfgs:
          break;
      }
    throw std::invalid_argument("B1 has no artifact slug for this method");
  }

  const char *
  b1_beta_slug(const double beta)
  {
    if (std::abs(beta - 1.0e-1) < 1.0e-15)
      return "1e-1";
    if (std::abs(beta - 1.0e-2) < 1.0e-15)
      return "1e-2";
    if (std::abs(beta - 1.0e-3) < 1.0e-15)
      return "1e-3";
    if (std::abs(beta - 1.0e-6) < 1.0e-18)
      return "1e-6";
    throw std::invalid_argument(
      "B1 has no frozen artifact slug for the regularisation value");
  }

  std::string
  b1_mesh_provenance(const unsigned int refinement)
  {
    return "chapter-6.e6.5.1.framework-native-hypercube-r" +
           std::to_string(refinement);
  }

  std::string
  b2_mesh_provenance(const unsigned int refinement)
  {
    return "chapter-6.e6.5.2.framework-native-rectangle-r" +
           std::to_string(refinement);
  }

  std::string
  b1_number(const double value)
  {
    std::ostringstream output;
    output.precision(17);
    output << value;
    return output.str();
  }

  const std::string &
  combination_value(const nmopt::application::runner::ParameterCombination &combination,
                    const std::string_view                                      axis)
  {
    const auto found = combination.values.find(std::string(axis));
    if (found == combination.values.end())
      throw std::invalid_argument("resolved combination has no matrix axis '" +
                                  std::string(axis) + "'");
    return found->second;
  }

  nmopt::application::chapter6::ReducedMethod
  parse_method(const std::string &value)
  {
    using nmopt::application::chapter6::ReducedMethod;
    if (value == "steepest-descent")
      return ReducedMethod::steepest_descent;
    if (value == "l-bfgs")
      return ReducedMethod::limited_memory_bfgs;
    if (value == "bfgs")
      return ReducedMethod::bfgs;
    throw std::invalid_argument("unknown solver method '" + value + "'");
  }

  nmopt::application::chapter6::MeshGeneration
  parse_mesh_generation(const std::string &value)
  {
    using nmopt::application::chapter6::MeshGeneration;
    if (value == "framework-native")
      return MeshGeneration::framework_native;
    if (value == "structured-simplex")
      return MeshGeneration::structured_simplex;
    if (value == "centroid-split-simplex")
      return MeshGeneration::centroid_split_simplex;
    throw std::invalid_argument("unknown mesh generator '" + value + "'");
  }

  double
  parse_number_text(const std::string &text, const std::string &key)
  {
    if (text.empty() || text == "none")
      throw std::invalid_argument("parameter '" + key + "' needs a number");
    std::size_t consumed = 0;
    double      value = 0.0;
    try
      {
        value = std::stod(text, &consumed);
      }
    catch (const std::exception &)
      {
        throw std::invalid_argument("parameter '" + key + "' needs a number");
      }
    if (consumed != text.size() || !std::isfinite(value))
      throw std::invalid_argument("parameter '" + key + "' needs a finite number");
    return value;
  }

  unsigned int
  parse_unsigned_text(const std::string &text, const std::string &key)
  {
    const auto value = parse_number_text(text, key);
    if (value < 0.0 || value > std::numeric_limits<unsigned int>::max() ||
        value != std::floor(value))
      throw std::invalid_argument("parameter '" + key +
                                  "' needs a nonnegative integer");
    return static_cast<unsigned int>(value);
  }

  nmopt::solvers::ReducedStoppingCriterion
  parse_stopping_criterion(const std::string &value)
  {
    using Criterion = nmopt::solvers::ReducedStoppingCriterion;
    if (value == "automatic")
      return Criterion::automatic;
    if (value == "gradient-norm")
      return Criterion::gradient_norm;
    if (value == "relative-gradient-norm")
      return Criterion::relative_gradient_norm;
    if (value == "objective-change")
      return Criterion::objective_change;
    if (value == "step-norm")
      return Criterion::step_norm;
    throw std::invalid_argument("unknown reduced stopping criterion '" +
                                value + "'");
  }

  void
  require_parameter(const nmopt::application::runner::ParameterFile &file,
                    const std::string_view                         key,
                    const std::string_view                         expected)
  {
    if (file.value(key) != expected)
      throw std::invalid_argument("parameter '" + std::string(key) +
                                  "' is '" + file.value(key) +
                                  "', but this adapter requires '" +
                                  std::string(expected) + "'");
  }

  std::string
  join_matrix(const nmopt::application::runner::ParameterFile &file)
  {
    std::ostringstream output;
    for (std::size_t index = 0; index < file.matrix.size(); ++index)
      {
        if (index != 0)
          output << ';';
        output << file.matrix[index].id << '=';
        for (std::size_t value = 0; value < file.matrix[index].values.size(); ++value)
          {
            if (value != 0)
              output << ',';
            output << file.matrix[index].values[value];
          }
      }
    return output.str();
  }

  std::string
  join_combinations(const std::vector<nmopt::application::runner::ParameterCombination> &combinations)
  {
    std::ostringstream output;
    for (std::size_t index = 0; index < combinations.size(); ++index)
      {
        if (index != 0)
          output << ';';
        output << '[';
        std::size_t value_index = 0;
        for (const auto &[axis, value] : combinations[index].values)
          {
            if (value_index++ != 0)
              output << ',';
            output << axis << '=' << value;
          }
        output << ']';
      }
    return output.str();
  }

  std::string
  join_selection(const nmopt::application::runner::ParameterFile &file,
                 const std::vector<std::pair<std::string, std::string>> &cli)
  {
    std::ostringstream output;
    bool                first = true;
    for (const auto &[axis, value] : file.selection)
      {
        if (!first)
          output << ',';
        first = false;
        output << axis << '=' << value;
      }
    for (const auto &[axis, value] : cli)
      {
        if (!first)
          output << ',';
        first = false;
        output << axis << '=' << value;
      }
    return output.str();
  }

  std::vector<std::string>
  command_line_arguments(const int argc, char **argv)
  {
    std::vector<std::string> command;
    command.reserve(static_cast<std::size_t>(argc));
    for (int index = 0; index < argc; ++index)
      command.emplace_back(argv[index]);
    return command;
  }

  std::vector<std::string>
  b1_expected_artifacts(
    const nmopt::application::runner::ParameterFile &file,
    const std::vector<std::pair<std::string, std::string>> &filters)
  {
    std::vector<std::string> paths;
    for (const auto &combination : file.combinations(filters))
      {
        const auto method = parse_method(combination_value(combination, "method"));
        const auto beta = parse_number_text(
          combination_value(combination, "regularisation"), "Matrix/regularisation");
        paths.push_back("artifacts/" + std::string(b1_method_slug(method)) +
                        "/beta-" + b1_beta_slug(beta) + "/artifact.kv");
      }
    return paths;
  }

  std::vector<std::string>
  b2_expected_artifacts(
    const nmopt::application::runner::ParameterFile &file,
    const std::vector<std::pair<std::string, std::string>> &filters)
  {
    std::vector<std::string> paths;
    const auto combinations = file.combinations(filters);
    const bool standard_case_matrix =
      file.matrix.size() == 2 &&
      file.matrix[0].id == "observation-region" &&
      file.matrix[1].id == "target-profile";
    for (const auto &combination : combinations)
      {
        std::vector<std::string> components;
        if (standard_case_matrix)
          {
            const auto region = combination_value(combination, "observation-region");
            const auto target = combination_value(combination, "target-profile");
            components.push_back(region + "-" + target);
          }
        else
          for (const auto &axis : file.matrix)
            components.push_back(axis.id + "-" + combination_value(combination, axis.id));
        std::string relative = "artifacts";
        for (const auto &component : components)
          relative += "/" + component;
        paths.push_back(relative + "/artifact.kv");
      }
    return paths;
  }

  template <typename Scenario>
  void
  apply_common_parameter_options(
    Scenario &scenario,
    const nmopt::application::runner::ParameterFile &file,
    const std::string &scenario_id)
  {
    using namespace nmopt::application::runner;
    scenario.experiment.scenario_output_id = scenario_id;
    scenario.experiment.source_reference = file.value("Benchmark/source reference");
    scenario.experiment.source_revision = file.value("Benchmark/source revision");
    scenario.experiment.retain_fields = parameter_bool(file, "Output/retain fields");
    scenario.experiment.harness.deterministic =
      parameter_bool(file, "Run/deterministic");
    scenario.experiment.harness.serialize_artifacts =
      parameter_bool(file, "Run/serialize artifacts");
    scenario.experiment.harness.measure_timings =
      parameter_bool(file, "Run/measure timings");
    scenario.experiment.harness.measure_memory =
      parameter_bool(file, "Run/measure memory");
    scenario.experiment.harness.artifact_directory =
      file.value("Run/output root");

    scenario.compile.mesh.dimension = parameter_unsigned(file, "Mesh/dimension");
    scenario.compile.mesh.generation =
      parse_mesh_generation(file.value("Mesh/generator"));
    scenario.compile.mesh.refinement = parameter_unsigned(file, "Mesh/refinement");
    scenario.compile.mesh.subdivisions =
      parameter_unsigned(file, "Mesh/subdivisions");
    scenario.compile.mesh.axis_subdivisions =
      parameter_positive_unsigned_list(file, "Mesh/axis subdivisions");
    scenario.compile.mesh.centroid_splits =
      parameter_unsigned(file, "Mesh/centroid splits");
    scenario.compile.mesh.selection_seed =
      parameter_unsigned(file, "Mesh/selection seed");
    scenario.compile.mesh.mesh_provenance = file.value("Mesh/provenance");
    scenario.compile.state_degree = parameter_unsigned(file, "Compile/state degree");
    scenario.compile.owned_session = parameter_bool(file, "Compile/owned session");
    scenario.compile.state_solve = {
      parameter_unsigned(file, "Compile/state solve maximum iterations"),
      parameter_double(file, "Compile/state solve relative tolerance"),
      parameter_double(file, "Compile/state solve absolute tolerance")};
    scenario.compile.adjoint_solve = {
      parameter_unsigned(file, "Compile/adjoint solve maximum iterations"),
      parameter_double(file, "Compile/adjoint solve relative tolerance"),
      parameter_double(file, "Compile/adjoint solve absolute tolerance")};
    scenario.compile.control_metric_solve = {
      parameter_unsigned(file,
                         "Compile/control metric solve maximum iterations"),
      parameter_double(file,
                       "Compile/control metric solve relative tolerance"),
      parameter_double(file,
                       "Compile/control metric solve absolute tolerance")};
    require_parameter(file, "Compile/execution", "assembled");
    require_parameter(file, "Compile/product", "reduced-dto");
  }

  void
  apply_solver_options(
    nmopt::application::chapter6::SolverOptions &solver,
    const nmopt::application::runner::ParameterFile &file,
    const std::string &method_id)
  {
    using namespace nmopt::application::runner;
    solver.globalization = reduced_globalization(file);
    const auto method_prefix = "Solver/method policy " + method_id + "/";
    const auto resolve = [&](const std::string_view entry) {
      return resolve_method_parameter(file, method_id, entry);
    };

    solver.initial_control = file.value("Solver/initial control");
    const auto maximum_iterations = resolve("maximum iterations");
    solver.parameters.maximum_iterations = parse_unsigned_text(
      maximum_iterations.value, maximum_iterations.key);

    ResolvedParameterValue trials;
    ResolvedParameterValue reductions;
    const auto method_trials =
      file.optional_value(method_prefix + "maximum line search trials");
    const auto method_reductions =
      file.optional_value(method_prefix + "maximum backtracking reductions");
    if (!method_trials.empty() || !method_reductions.empty())
      {
        trials = {method_prefix + "maximum line search trials", method_trials};
        reductions = {method_prefix + "maximum backtracking reductions",
                      method_reductions};
      }
    else
      {
        trials = {"Solver/maximum line search trials",
                  file.optional_value("Solver/maximum line search trials")};
        reductions = {
          "Solver/maximum backtracking reductions",
          file.optional_value("Solver/maximum backtracking reductions")};
      }
    if (!trials.value.empty() && !reductions.value.empty())
      throw std::invalid_argument(
        "select either maximum line-search trials or maximum backtracking reductions");
    if (!reductions.value.empty())
      {
        const auto count = parse_unsigned_text(reductions.value, reductions.key);
        if (count == std::numeric_limits<unsigned int>::max())
          throw std::invalid_argument("parameter '" + reductions.key +
                                      "' is too large");
        solver.parameters.maximum_line_search_trials = count + 1;
      }
    else if (!trials.value.empty())
      solver.parameters.maximum_line_search_trials =
        parse_unsigned_text(trials.value, trials.key);

    const auto stopping = resolve("stopping criterion");
    solver.parameters.stopping_criterion =
      parse_stopping_criterion(stopping.value);
    const auto relative = resolve("relative gradient tolerance");
    solver.parameters.relative_gradient_tolerance =
      parse_number_text(relative.value, relative.key);
    const auto objective_change = resolve("objective change tolerance");
    solver.parameters.objective_change_tolerance =
      parse_number_text(objective_change.value, objective_change.key);
    const auto step = resolve("step tolerance");
    solver.parameters.step_tolerance = parse_number_text(step.value, step.key);
    const auto initial_step = resolve("initial step length");
    solver.parameters.initial_step_length =
      parse_number_text(initial_step.value, initial_step.key);
    const auto armijo = resolve("Armijo fraction");
    solver.parameters.armijo_fraction =
      parse_number_text(armijo.value, armijo.key);
    const auto backtracking = resolve("backtracking factor");
    solver.parameters.backtracking_factor =
      parse_number_text(backtracking.value, backtracking.key);

    const auto objective_target_policy =
      file.optional_value("Solver/objective target policy", "none");
    const auto objective_target = file.optional_value("Solver/objective target");
    if (objective_target_policy.empty() || objective_target_policy == "none")
      {
        solver.objective_target_policy =
          nmopt::application::chapter6::ObjectiveTargetPolicy::none;
        solver.parameters.objective_target = std::nullopt;
      }
    else if (objective_target_policy == "explicit")
      {
        if (objective_target.empty() || objective_target == "none")
          throw std::invalid_argument(
            "an explicit objective-target policy needs a numeric target");
        solver.objective_target_policy =
          nmopt::application::chapter6::ObjectiveTargetPolicy::explicit_value;
        solver.parameters.objective_target = parse_number_text(
          objective_target, "Solver/objective target");
      }
    else if (objective_target_policy == "match-reference-method")
      {
        solver.objective_target_policy = nmopt::application::chapter6::
          ObjectiveTargetPolicy::matched_reference_method;
        solver.parameters.objective_target = std::nullopt;
        solver.objective_target_reference_method =
          file.value("Solver/objective target reference method");
        if (solver.objective_target_reference_method.empty())
          throw std::invalid_argument(
            "a matched objective target needs a reference method");
      }
    else
      throw std::invalid_argument("unknown objective-target policy '" +
                                  objective_target_policy + "'");

    const auto gradient = resolve("gradient tolerance");
    solver.parameters.gradient_tolerance =
      parse_number_text(gradient.value, gradient.key);
    const auto minimum = resolve("minimum step length");
    if (!minimum.value.empty() && minimum.value != "none")
      solver.parameters.minimum_step_length =
        parse_number_text(minimum.value, minimum.key);
    else
      solver.parameters.minimum_step_length = 0.0;

    const auto declared_minimum = resolve("declared minimum step length");
    if (!declared_minimum.value.empty() && declared_minimum.value != "none")
      solver.declared_minimum_step_length =
        parse_number_text(declared_minimum.value, declared_minimum.key);
    else
      solver.declared_minimum_step_length = std::nullopt;

    const auto memory = file.optional_value(method_prefix + "memory");
    if (!memory.empty())
      solver.limited_memory_bfgs.memory_size =
        parse_unsigned_text(memory, method_prefix + "memory");
    const auto curvature =
      file.optional_value(method_prefix + "curvature tolerance");
    if (!curvature.empty())
      solver.limited_memory_bfgs.curvature_tolerance = parse_number_text(
        curvature, method_prefix + "curvature tolerance");
    const auto scaling =
      file.optional_value(method_prefix + "initial inverse Hessian scaling");
    if (!scaling.empty() && scaling != "metric-inverse" &&
        scaling != "scalar-secant")
      throw std::invalid_argument("unknown L-BFGS initial scaling '" + scaling +
                                  "'");
    if (scaling == "metric-inverse")
      solver.limited_memory_bfgs.initial_inverse_hessian_scaling =
        nmopt::solvers::LimitedMemoryBfgsInitialScaling::metric_inverse;
    else if (scaling == "scalar-secant")
      solver.limited_memory_bfgs.initial_inverse_hessian_scaling =
        nmopt::solvers::LimitedMemoryBfgsInitialScaling::scalar_secant;
  }

  void
  configure_b1_scenario(
    nmopt::application::chapter6::B1Scenario &scenario,
    const nmopt::application::runner::ParameterFile &file,
    const nmopt::application::runner::ParameterCombination &combination,
    const std::string &method_id,
    const double beta)
  {
    using namespace nmopt::application;
    using namespace runner;
    require_parameter(file,
                      "Benchmark/id",
                      "chapter-6.b1.distributed-laplace");
    require_parameter(file,
                      "Benchmark/recipe",
                      chapter6::b1_recipe_id);
    require_parameter(file, "Problem/observation", "full-domain");
    require_parameter(file, "Functions/desired state", "b1-polynomial");
    require_parameter(file, "Functions/desired state/kind", "polynomial");
    require_parameter(file,
                      "Functions/desired state/expression",
                      "10*x0*(1-x0)*x1*(1-x1)");
    require_parameter(file,
                      "Mesh/geometry",
                      "unit-hypercube");
    require_parameter(file, "Solver/method", "from-matrix");
    require_parameter(file, "Output/selected fields", "state, control, adjoint, negative-adjoint, target, forcing");
    if (combination.values.size() < 2)
      throw std::invalid_argument("B1 needs method and regularisation matrix axes");

    const auto control_representation =
      file.value("Problem/control representation");
    if (control_representation == "cellwise-volume")
      scenario.problem.recipe.discretisation =
        chapter5::DistributedControlDiscretisation::cellwise_constant;
    else if (control_representation ==
             "continuous-volume-homogeneous-dirichlet")
      scenario.problem.recipe.discretisation = chapter5::
        DistributedControlDiscretisation::homogeneous_dirichlet_continuous;
    else
      throw std::invalid_argument(
        "B1 has an unknown control representation '" +
        control_representation + "'");
    scenario.problem.recipe.with_cellwise_box =
      parameter_bool(file, "Problem/cellwise box constraint");
    scenario.problem.data.diffusion = parameter_double(file, "Runtime/diffusion");
    scenario.problem.data.reaction = parameter_double(file, "Runtime/reaction");
    scenario.problem.data.regularisation_weight = beta;
    scenario.problem.data.desired_state_provenance =
      file.value("Functions/desired state/provenance");
    scenario.problem.regularisation_sweep = {beta};
    scenario.problem.forcing =
      parameter_scalar_function_definition(file, "Functions/forcing");
    scenario.problem.data.forcing_provenance =
      scenario.problem.forcing.provenance;
    scenario.solver.method = parse_method(method_id);
    apply_common_parameter_options(
      scenario,
      file,
      std::string(file.value("Benchmark/id")) + "." + method_id + ".beta-" +
        b1_beta_slug(beta));
    apply_solver_options(scenario.solver, file, method_id);
    if (scenario.solver.method == chapter6::ReducedMethod::bfgs)
      throw std::invalid_argument("B1 does not support full BFGS");
  }

  nmopt::application::chapter6::GraetzCase
  graetz_case_from_combination(
    const nmopt::application::runner::ParameterCombination &combination)
  {
    const auto region = combination.values.find("observation-region");
    const auto target = combination.values.find("target-profile");
    if (region != combination.values.end() && target != combination.values.end())
      {
        using namespace nmopt::application::chapter6;
        if (region->second == "wings" && target->second == "constant")
          return GraetzCase::observation_wings_constant_target;
        if (region->second == "full" && target->second == "constant")
          return GraetzCase::observation_full_constant_target;
        if (region->second == "wings" && target->second == "parabolic")
          return GraetzCase::observation_wings_parabolic_target;
        if (region->second == "full" && target->second == "parabolic")
          return GraetzCase::observation_full_parabolic_target;
      }
    const auto case_value = combination.values.find("case");
    if (case_value != combination.values.end())
      {
        using namespace nmopt::application::chapter6;
        for (const auto candidate : b2_case_order)
          if (case_value->second == graetz_case_name(candidate))
            return candidate;
      }
    throw std::invalid_argument(
      "B2 combinations need observation-region/target-profile or case axes");
  }

  void
  configure_b2_scenario(
    nmopt::application::chapter6::B2Scenario &scenario,
    const nmopt::application::runner::ParameterFile &file,
    const nmopt::application::runner::ParameterCombination &combination,
    const std::string &scenario_id)
  {
    using namespace nmopt::application;
    using namespace runner;
    require_parameter(file, "Benchmark/id", "chapter-6.b2.graetz-flow");
    require_parameter(file, "Benchmark/recipe", chapter6::b2_recipe_id);
    require_parameter(file, "Problem/initial control", "zero");
    require_parameter(file, "Functions/fixed Dirichlet data", "fixed-temperature");
    require_parameter(file, "Functions/conservative transport", "graetz");
    require_parameter(file, "Functions/fixed-temperature/kind", "constant");
    require_parameter(file, "Functions/graetz/kind", "conservative-transport");
    require_parameter(file,
                      "Functions/graetz/expression",
                      "(1.5*x1*(1-x1), 0.0)");
    scenario.problem.transport_boundary_form =
      b2_transport_boundary_form(file);
    scenario.problem.recipe.control_discretisation =
      b2_neumann_control_discretisation(file);
    scenario.compile.volume_observation =
      b2_volume_observation_options(file);
    require_parameter(file, "Boundary/normal orientation", "outward");
    require_parameter(file, "Boundary/trace evaluation", "fe-q-state-trace");
    require_parameter(file, "Boundary/face quadrature", "qgauss-face");
    require_parameter(file, "Boundary/fixed region", "dirichlet-boundary");
    require_parameter(file, "Boundary/control region", "control-boundary");
    require_parameter(file, "Boundary/outflow region", "outflow-boundary");
    require_parameter(file, "Boundary/fixed boundary id", "0");
    require_parameter(file, "Boundary/control boundary id", "1");
    require_parameter(file, "Boundary/outflow boundary id", "2");
    require_parameter(file, "Boundary/upstream transition x", "1.0");
    require_parameter(file, "Boundary/outflow x", "4.0");
    require_parameter(file, "Mesh/geometry", "rectangle");
    require_parameter(file, "Mesh/lower", "(0.0, 0.0)");
    require_parameter(file, "Mesh/upper", "(4.0, 1.0)");
    require_parameter(file, "Compile/stabilization", "galerkin");
    require_parameter(file, "Solver/method", "bfgs");

    scenario.problem.recipe.observed_material_id =
      parameter_unsigned(file, "Problem/observed material id");
    scenario.problem.recipe.with_facewise_box =
      parameter_bool(file, "Problem/facewise box constraint");
    scenario.problem.data.diffusion = parameter_double(file, "Runtime/diffusion");
    scenario.problem.data.reaction = parameter_double(file, "Runtime/reaction");
    scenario.problem.data.regularisation_weight =
      combination.values.count("regularisation")
        ? parse_number_text(combination_value(combination, "regularisation"),
                            "Matrix/regularisation")
        : parameter_double(file, "Runtime/regularisation");
    scenario.problem.fixed_temperature =
      parse_number_text(file.value("Functions/fixed-temperature/value"),
                        "Functions/fixed-temperature/value");
    scenario.problem.data.fixed_dirichlet_data_provenance =
      file.value("Functions/fixed-temperature/provenance");
    scenario.problem.data.conservative_transport_provenance =
      file.value("Functions/graetz/provenance");
    scenario.problem.graetz_case = graetz_case_from_combination(combination);
    const auto region = chapter6::dealii::b2_observation_region(
      scenario.problem.graetz_case);
    require_parameter(file,
                      "Observation/active region",
                      "from-matrix");
    require_parameter(file,
                      "Observation/region " + std::string(region) + "/geometry",
                      region == std::string("wings")
                        ? "x0 > 1.0 and (x1 < 0.3 or x1 > 0.7)"
                        : "x0 > 1.0");
    require_parameter(file, "Functions/desired state", "from-matrix");
    scenario.problem.target_parameters = b2_target_parameters(file);

    const auto forcing_axis = combination.values.find("forcing");
    const std::string forcing_id =
      forcing_axis == combination.values.end()
        ? file.value("Functions/forcing")
        : forcing_axis->second;
    const std::string forcing_prefix =
      forcing_axis == combination.values.end()
        ? "Functions/forcing"
        : "Functions/forcing definition " + forcing_id;
    scenario.problem.forcing =
      forcing_axis == combination.values.end()
        ? parameter_scalar_function_definition(file, forcing_prefix)
        : parameter_scalar_function_definition_with_id(
            file, forcing_prefix, forcing_id);
    scenario.problem.data.forcing_provenance =
      scenario.problem.forcing.provenance;
    scenario.problem.data.desired_state_provenance =
      chapter6::b2_target_definition(scenario.problem.target_parameters,
                                     scenario.problem.graetz_case)
        .provenance;

    apply_common_parameter_options(scenario, file, scenario_id);
    apply_solver_options(scenario.solver, file, "bfgs");
  }

  template <typename Scenario>
  nmopt::experiment::RunEnvironmentRecord
  make_environment(const Scenario &scenario)
  {
#if defined(__clang__)
    const std::string compiler = "clang";
#elif defined(__GNUC__)
    const std::string compiler = "gcc";
#elif defined(_MSC_VER)
    const std::string compiler = "msvc";
#else
    const std::string compiler = "unknown";
#endif

#if defined(__x86_64__) || defined(_M_X64)
    const std::string architecture = "x86_64";
#elif defined(__aarch64__) || defined(_M_ARM64)
    const std::string architecture = "aarch64";
#else
    const std::string architecture = "unknown";
#endif

    const char *hardware = std::getenv("NMOPT_HARDWARE");
    return {scenario.experiment.source_revision,
            NMOPT_COMPILED_BUILD_PROFILE,
            compiler,
#if defined(__VERSION__)
            __VERSION__,
#else
            "unknown",
#endif
            "libstdc++",
            "Linux",
            architecture,
            hardware == nullptr ? "unspecified" : hardware};
  }

  using Chapter6Evidence =
    nmopt::application::benchmark::BenchmarkExecutionEvidenceT<
      nmopt::application::chapter6::dealii::Envelope>;

  void
  add_common_artifact_fields(Chapter6Evidence &evidence,
                             const std::string &framework_revision)
  {
    const auto &manifest = evidence.envelope.compilation_manifest();
    const auto &environment = evidence.envelope.environment();
    const auto algorithm_name = [](const auto algorithm) {
      switch (algorithm)
        {
          case nmopt::compiler::v1::LinearSolveAlgorithm::serial_cg:
            return "serial_cg";
          case nmopt::compiler::v1::LinearSolveAlgorithm::
            serial_sparse_direct_umfpack:
            return "serial_sparse_direct_umfpack";
        }
      return "unknown";
    };
    const auto add_solve_policy = [&](const std::string &prefix,
                                      const auto &record) {
      evidence.fields.push_back(
        {prefix + ".algorithm", algorithm_name(record.algorithm)});
      evidence.fields.push_back(
        {prefix + ".preconditioner", record.preconditioner});
      evidence.fields.push_back(
        {prefix + ".maximum_iterations",
         std::to_string(record.maximum_iterations)});
      evidence.fields.push_back(
        {prefix + ".relative_tolerance",
         b1_number(record.relative_tolerance)});
      evidence.fields.push_back(
        {prefix + ".absolute_tolerance",
         b1_number(record.absolute_tolerance)});
    };
    evidence.fields.push_back({"provenance.framework_revision",
                               framework_revision});
    evidence.fields.push_back({
      "provenance.recipe_revision",
      "application-chapter6@" + framework_revision});
    evidence.fields.push_back({"provenance.mesh_provenance",
                               manifest.mesh_record.provenance});
    evidence.fields.push_back({"provenance.mesh_structural_identity",
                               manifest.mesh_record.structural_identity});
    evidence.fields.push_back({"provenance.compiler", environment.compiler});
    evidence.fields.push_back({"provenance.compiler_version",
                               environment.compiler_version});
    evidence.fields.push_back({"provenance.standard_library",
                               environment.standard_library});
    evidence.fields.push_back({"provenance.operating_system",
                               environment.operating_system});
    evidence.fields.push_back({"provenance.architecture",
                               environment.architecture});
    evidence.fields.push_back({"provenance.hardware", environment.hardware});
    evidence.fields.push_back({"manifest.schema_version",
                               std::to_string(manifest.schema_version)});
    evidence.fields.push_back({"manifest.semantic_problem_id",
                               manifest.semantic_problem_id});
    evidence.fields.push_back({"manifest.compiler_id", manifest.compiler_id});
    evidence.fields.push_back({"manifest.backend", manifest.backend});
    evidence.fields.push_back({"manifest.execution", manifest.execution});
    evidence.fields.push_back({"manifest.state_space", manifest.state_space});
    evidence.fields.push_back({"manifest.control_space",
                               manifest.control_space});
    evidence.fields.push_back({"manifest.mesh_dimension",
                               std::to_string(manifest.mesh_record.dimension)});
    evidence.fields.push_back({"manifest.active_cells",
                               std::to_string(manifest.mesh_record.active_cells)});
    evidence.fields.push_back({"manifest.space_count",
                               std::to_string(manifest.spaces.size())});
    evidence.fields.push_back({"manifest.binding_count",
                               std::to_string(manifest.bindings.size())});
    evidence.fields.push_back({"manifest.realized_space_count",
                               std::to_string(manifest.realized_spaces.size())});
    evidence.fields.push_back({"manifest.realized_map_count",
                               std::to_string(manifest.realized_maps.size())});
    for (const auto &space : manifest.spaces)
      {
        const auto prefix = "manifest.space." + space.semantic_id;
        evidence.fields.push_back(
          {prefix + ".dimension", std::to_string(space.dimension)});
        evidence.fields.push_back(
          {prefix + ".runtime_role", space.runtime_role});
        evidence.fields.push_back(
          {prefix + ".finite_element", space.finite_element});
        evidence.fields.push_back({prefix + ".region", space.region_id});
      }
    add_solve_policy("manifest.state_solve", manifest.state_solve_record);
    add_solve_policy("manifest.adjoint_solve", manifest.adjoint_solve_record);
    add_solve_policy("manifest.control_metric_solve",
                     manifest.metric_record.solve_policy);
    evidence.fields.push_back({
      "solver.final_objective",
      b1_number(evidence.envelope.report().final_evaluation.objective_value)});
    const auto &parameters = evidence.envelope.report().parameters;
    const auto &line_search = evidence.envelope.report().policy_parameters;
    evidence.fields.push_back({"solver.maximum_iterations",
                               std::to_string(parameters.maximum_iterations)});
    evidence.fields.push_back({"solver.maximum_line_search_trials",
                               std::to_string(line_search.maximum_trials)});
    evidence.fields.push_back({"solver.maximum_backtracking_reductions",
                               std::to_string(
                                 line_search.maximum_backtracking_reductions)});
    evidence.fields.push_back({"solver.gradient_tolerance",
                               b1_number(parameters.gradient_tolerance)});
    evidence.fields.push_back(
      {"solver.stopping_criterion",
       nmopt::solvers::reduced_stopping_criterion_name(
         parameters.stopping_criterion)});
    evidence.fields.push_back(
      {"solver.relative_gradient_tolerance",
       b1_number(parameters.relative_gradient_tolerance)});
    evidence.fields.push_back(
      {"solver.objective_change_tolerance",
       b1_number(parameters.objective_change_tolerance)});
    evidence.fields.push_back({"solver.step_tolerance",
                               b1_number(parameters.step_tolerance)});
    evidence.fields.push_back({"solver.initial_step_length",
                               b1_number(line_search.initial_step_length)});
    evidence.fields.push_back({"solver.minimum_step_length",
                               b1_number(line_search.minimum_step_length)});
    evidence.fields.push_back({"solver.armijo_fraction",
                               b1_number(line_search.armijo_fraction)});
    evidence.fields.push_back({"solver.backtracking_factor",
                               b1_number(line_search.backtracking_factor)});
    if (parameters.objective_target.has_value())
      evidence.fields.push_back(
        {"solver.objective_target",
         b1_number(*parameters.objective_target)});
  }

  void
  add_parameter_artifact_fields(
    Chapter6Evidence &evidence,
    const nmopt::application::runner::ResolvedRunConfiguration &configuration,
    const nmopt::application::runner::ParameterCombination &combination)
  {
    evidence.fields.push_back({"provenance.parameters_file",
                               configuration.parameter_file.string()});
    evidence.fields.push_back({"provenance.parameters_hash",
                               configuration.parameter_hash});
    evidence.fields.push_back({"provenance.plotting_profile",
                               configuration.plotting_profile_file.string()});
    evidence.fields.push_back({"provenance.plotting_profile_hash",
                               configuration.plotting_profile_hash});
    for (const auto &[axis, value] : combination.values)
      evidence.fields.push_back({"parameters." + axis, value});
  }

  void
  add_b1_artifact_fields(
    nmopt::application::benchmark::BenchmarkExecutionEvidenceT<
      nmopt::application::chapter6::dealii::Envelope> &evidence,
    const nmopt::application::chapter6::B1Scenario &scenario,
    const nmopt::application::chapter6::ReducedMethod method,
    const double beta,
    const std::string &framework_revision,
    const nmopt::application::runner::RunKind run_kind,
    const std::string &objective_target_reference_artifact = {})
  {
    add_common_artifact_fields(evidence, framework_revision);
    evidence.fields.push_back({"benchmark.method", b1_method_slug(method)});
    evidence.fields.push_back({"benchmark.regularisation", b1_number(beta)});
    evidence.fields.push_back(
      {"benchmark.control_discretisation",
       nmopt::application::chapter5::distributed_control_discretisation_name(
         scenario.problem.recipe.discretisation)});
    evidence.fields.push_back(
      {"benchmark.mesh_refinement", std::to_string(scenario.compile.mesh.refinement)});
    evidence.fields.push_back(
      {"benchmark.mesh_generator",
       nmopt::application::chapter6::mesh_generation_name(
         scenario.compile.mesh.generation)});
    evidence.fields.push_back(
      {"benchmark.mesh_subdivisions",
       std::to_string(scenario.compile.mesh.subdivisions)});
    evidence.fields.push_back(
      {"benchmark.mesh_centroid_splits",
       std::to_string(scenario.compile.mesh.centroid_splits)});
    evidence.fields.push_back(
      {"benchmark.mesh_selection_seed",
       std::to_string(scenario.compile.mesh.selection_seed)});
    evidence.fields.push_back(
      {"benchmark.run_kind",
       std::string(nmopt::application::runner::run_kind_name(run_kind))});
    evidence.fields.push_back({"provenance.forcing",
                               scenario.problem.data.forcing_provenance});
    evidence.fields.push_back({"provenance.desired_state",
                               scenario.problem.data.desired_state_provenance});
    switch (scenario.solver.objective_target_policy)
      {
        case nmopt::application::chapter6::ObjectiveTargetPolicy::none:
          evidence.fields.push_back({"solver.objective_target_policy", "none"});
          break;
        case nmopt::application::chapter6::ObjectiveTargetPolicy::explicit_value:
          evidence.fields.push_back(
            {"solver.objective_target_policy", "explicit"});
          break;
        case nmopt::application::chapter6::ObjectiveTargetPolicy::
          matched_reference_method:
          evidence.fields.push_back(
            {"solver.objective_target_policy", "match-reference-method"});
          evidence.fields.push_back(
            {"solver.objective_target_reference_method",
             scenario.solver.objective_target_reference_method});
          if (!objective_target_reference_artifact.empty())
            evidence.fields.push_back(
              {"solver.objective_target_reference_artifact",
               objective_target_reference_artifact});
          break;
      }
  }

  void
  add_b2_artifact_fields(
    Chapter6Evidence &evidence,
    const nmopt::application::chapter6::B2Scenario &scenario,
    const std::string &framework_revision,
    const nmopt::application::runner::RunKind run_kind)
  {
    const auto  case_slug = nmopt::application::chapter6::graetz_case_name(
      scenario.problem.graetz_case);
    const auto &manifest = evidence.envelope.compilation_manifest();
    const auto &volume_observation = *scenario.compile.volume_observation;
    const std::string volume_observation_target =
      nmopt::application::chapter6::
        volume_observation_target_realisation_name(
          volume_observation.target_realisation);
    const std::string quadrature_order =
      std::to_string(volume_observation.quadrature_order);
    nmopt::contract::require(
      manifest.observation_realisation.find(
        "target=" + volume_observation_target) != std::string::npos &&
        manifest.observation_realisation.find("(" + quadrature_order + ")") !=
          std::string::npos,
      "B2 artifact manifest does not match the declared volume observation");
    evidence.fields.push_back({"benchmark.graetz_case", case_slug});
    evidence.fields.push_back(
      {"benchmark.control_discretisation",
       nmopt::application::chapter5::neumann_control_discretisation_name(
         scenario.problem.recipe.control_discretisation)});
    evidence.fields.push_back(
      {"benchmark.volume_observation_quadrature_order", quadrature_order});
    evidence.fields.push_back(
      {"benchmark.volume_observation_target_realisation",
       volume_observation_target});
    evidence.fields.push_back({"benchmark.fixed_temperature",
                               b1_number(scenario.problem.fixed_temperature)});
    evidence.fields.push_back({"benchmark.regularisation",
                               b1_number(
                                 scenario.problem.data.regularisation_weight)});
    evidence.fields.push_back(
      {"benchmark.mesh_refinement", std::to_string(scenario.compile.mesh.refinement)});
    evidence.fields.push_back(
      {"benchmark.mesh_generator",
       nmopt::application::chapter6::mesh_generation_name(
         scenario.compile.mesh.generation)});
    evidence.fields.push_back(
      {"benchmark.mesh_subdivisions",
       std::to_string(scenario.compile.mesh.subdivisions)});
    evidence.fields.push_back(
      {"benchmark.mesh_axis_subdivisions",
       join_unsigned_ids(scenario.compile.mesh.axis_subdivisions)});
    evidence.fields.push_back(
      {"benchmark.mesh_centroid_splits",
       std::to_string(scenario.compile.mesh.centroid_splits)});
    evidence.fields.push_back(
      {"benchmark.mesh_selection_seed",
       std::to_string(scenario.compile.mesh.selection_seed)});
    evidence.fields.push_back(
      {"benchmark.run_kind",
       std::string(nmopt::application::runner::run_kind_name(run_kind))});
    add_common_artifact_fields(evidence, framework_revision);
    evidence.fields.push_back({"provenance.forcing",
                               scenario.problem.data.forcing_provenance});
    evidence.fields.push_back({"provenance.desired_state",
                               scenario.problem.data.desired_state_provenance});
    evidence.fields.push_back({"provenance.fixed_temperature",
                               scenario.problem.data.fixed_dirichlet_data_provenance});
    evidence.fields.push_back({"provenance.conservative_transport",
                               scenario.problem.data.conservative_transport_provenance});
    evidence.fields.push_back({"provenance.observation_case", case_slug});
    evidence.fields.push_back({"manifest.control_metric_realisation",
                               manifest.metric_record.realisation_id});
    evidence.fields.push_back(
      {"manifest.volume_observation_quadrature_order", quadrature_order});
    evidence.fields.push_back(
      {"manifest.volume_observation_target_realisation",
       volume_observation_target});
    evidence.fields.push_back(
      {"manifest.volume_observation_realisation",
       manifest.observation_realisation});
    if (manifest.boundary_realisation)
      {
        const auto boundary_form =
          manifest.boundary_realisation->transport_boundary_form ==
              nmopt::semantic::v1::TransportBoundaryForm::ordinary_normal_minus_transport
            ? "ordinary-normal-minus-transport"
            : manifest.boundary_realisation->conormal_form ==
                nmopt::semantic::v1::ConormalForm::transport_minus_diffusion
              ? "transport-minus-diffusion-conormal"
              : "total-conormal";
        evidence.fields.push_back({"benchmark.transport_boundary_form",
                                   boundary_form});
        evidence.fields.push_back({"manifest.transport_boundary_form",
                                   boundary_form});
      }
    for (const auto &region : manifest.resolved_decision.regions)
      if (!region.boundary_ids.empty())
        evidence.fields.push_back(
          {"manifest.region." + region.semantic_id + ".boundary_ids",
           join_unsigned_ids(region.boundary_ids)});
  }

  void
  write_artifact(const std::filesystem::path &path, const std::string &document)
  {
    nmopt::application::runner::prepare_artifact_path(path);
    std::ofstream output(path);
    if (!output)
      throw std::runtime_error("could not open benchmark artifact '" +
                               path.string() + "'");
    output << document;
    if (!output)
      throw std::runtime_error("could not write benchmark artifact '" +
                               path.string() + "'");
  }

  template <typename Report>
  void
  write_solver_trace(const std::filesystem::path &path, const Report &report)
  {
    nmopt::application::runner::prepare_artifact_path(path);
    std::ofstream output(path);
    if (!output)
      throw std::runtime_error("could not open solver trace '" +
                               path.string() + "'");
    output.imbue(std::locale::classic());
    output << "iteration,trial,step_length,objective_value,actual_slope,"
               "sufficient_decrease_bound,objective_finite,slope_negative,"
               "accepted\n";
    output.precision(17);
    for (const auto &trial : report.line_search_trials)
      output << trial.iteration << ',' << trial.trial << ','
             << trial.step_length << ',' << trial.objective_value << ','
             << trial.actual_slope << ','
             << trial.sufficient_decrease_bound << ','
             << (trial.objective_finite ? "true" : "false") << ','
             << (trial.slope_negative ? "true" : "false") << ','
             << (trial.accepted ? "true" : "false") << '\n';
    if (!output)
      throw std::runtime_error("could not write solver trace '" +
                               path.string() + "'");
  }

  bool
  run_b1(const nmopt::application::runner::ResolvedRunConfiguration &configuration,
         const std::vector<std::string> &command,
         const nmopt::application::runner::ParameterFile &file,
         const std::vector<std::pair<std::string, std::string>> &filters)
  {
    using namespace nmopt::application;
    using namespace chapter6;
    using Runner = benchmark::HeadlessBenchmarkRunnerT<B1Scenario>;
    using Adapter =
      nmopt::application::chapter6::dealii::B1ReducedExecutionAdapterT<2>;
    const auto &output_directory = configuration.run_directory;
    runner::RunSetManifest run_manifest(
      configuration, command, b1_expected_artifacts(file, filters));
    bool all_artifacts_succeeded = true;

    struct ObjectiveReference
    {
      double      value;
      std::string artifact;
    };
    std::map<std::string, ObjectiveReference> objective_references;
    auto combinations = file.combinations(filters);
    const bool matched_objective_target =
      file.optional_value("Solver/objective target policy", "none") ==
      "match-reference-method";
    ReducedMethod reference_method = ReducedMethod::steepest_descent;
    if (matched_objective_target)
      {
        reference_method = parse_method(
          file.value("Solver/objective target reference method"));
        std::stable_sort(
          combinations.begin(),
          combinations.end(),
          [reference_method](const auto &left, const auto &right) {
            const bool left_is_reference =
              parse_method(combination_value(left, "method")) ==
              reference_method;
            const bool right_is_reference =
              parse_method(combination_value(right, "method")) ==
              reference_method;
            return left_is_reference && !right_is_reference;
          });
      }

    for (const auto &combination : combinations)
      {
        const auto method_id = combination_value(combination, "method");
        const auto method = parse_method(method_id);
        const auto method_slug = b1_method_slug(method);
        const auto beta = parse_number_text(
          combination_value(combination, "regularisation"), "Matrix/regularisation");
        const auto beta_slug = b1_beta_slug(beta);
        const auto path = runner::artifact_path(
          output_directory,
          {method_slug, "beta-" + std::string(beta_slug)});
        try
          {
            auto scenario = make_b1_scenario(method);
            configure_b1_scenario(scenario, file, combination, method_id, beta);
            scenario.experiment.build_profile = NMOPT_COMPILED_BUILD_PROFILE;
            if (configuration.refinement_override.has_value())
              {
                if (scenario.compile.mesh.generation !=
                    nmopt::application::chapter6::MeshGeneration::
                      framework_native)
                  throw std::invalid_argument(
                    "--refinement cannot override a B1 simplex mesh; change its parameter file instead");
                scenario.compile.mesh.refinement =
                  *configuration.refinement_override;
                scenario.compile.mesh.mesh_provenance =
                  b1_mesh_provenance(scenario.compile.mesh.refinement);
              }
            scenario.experiment.harness.artifact_directory = output_directory.string();

            std::string objective_target_reference_artifact;
            if (matched_objective_target && method != reference_method)
              {
                const auto reference = objective_references.find(beta_slug);
                if (reference == objective_references.end())
                  throw std::invalid_argument(
                    "objective-target reference method did not reach its stopping tolerance for beta " +
                    std::string(beta_slug));
                scenario.solver.parameters.objective_target =
                  reference->second.value;
                objective_target_reference_artifact =
                  reference->second.artifact;
              }

            nmopt::application::chapter6::dealii::B1SelectedDataT<2> data(
              scenario.problem.forcing);
            const auto runtime =
              nmopt::application::chapter6::dealii::make_b1_runtime_data(
                scenario, data);
            const auto session =
              nmopt::application::chapter6::dealii::make_b1_compilation_session<2>(
                scenario);
            const auto environment = make_environment(scenario);
            Adapter execute{beta,
                            runtime,
                            session,
                            environment,
                            path.parent_path() / "native"};
            Runner scenario_runner(scenario);
            const auto result = scenario_runner.run(
              [](const auto &parameters) {
                return chapter6::make_b1_problem_spec(parameters);
              },
              [&](const auto &specification, const auto &run_scenario) {
                auto evidence = execute(specification, run_scenario);
                add_b1_artifact_fields(evidence,
                                       run_scenario,
                                       method,
                                       beta,
                                       configuration.framework_revision,
                                       configuration.run_kind,
                                       objective_target_reference_artifact);
                add_parameter_artifact_fields(evidence, configuration, combination);
                return evidence;
              });

            write_artifact(path, result.document);
            write_solver_trace(path.parent_path() / "solver-trace.csv",
                               result.artifact.envelope().report());
            run_manifest.record_success(path);
            if (matched_objective_target && method == reference_method &&
                runner::reference_reached_stopping_tolerance(
                  scenario.solver.parameters.stopping_criterion,
                  result.artifact.envelope().report().stopping_reason))
              objective_references[beta_slug] = {
                result.artifact.envelope().report().final_evaluation.
                  objective_value,
                std::filesystem::relative(path, output_directory).string()};
            std::cout << "B1 wrote " << path.string() << '\n';
          }
        catch (const std::exception &error)
          {
            all_artifacts_succeeded = false;
            run_manifest.record_failure(path, error.what());
            std::cerr << "B1 failed " << path.string() << ": "
                      << error.what() << '\n';
          }
      }
    return run_manifest.finalize() && all_artifacts_succeeded;
  }

  bool
  run_b2(const nmopt::application::runner::ResolvedRunConfiguration &configuration,
         const std::vector<std::string> &command,
         const nmopt::application::runner::ParameterFile &file,
         const std::vector<std::pair<std::string, std::string>> &filters)
  {
    using namespace nmopt::application;
    using namespace chapter6;
    using Runner = benchmark::HeadlessBenchmarkRunnerT<B2Scenario>;
    using Adapter =
      nmopt::application::chapter6::dealii::B2ReducedExecutionAdapterT<2>;
    const auto &output_directory = configuration.run_directory;
    runner::RunSetManifest run_manifest(
      configuration, command, b2_expected_artifacts(file, filters));
    bool all_artifacts_succeeded = true;

    for (const auto &combination : file.combinations(filters))
      {
        const auto graetz_case = graetz_case_from_combination(combination);
        const auto case_slug = graetz_case_name(graetz_case);
        std::vector<std::string> components;
        const bool standard_case_matrix =
          file.matrix.size() == 2 &&
          file.matrix[0].id == "observation-region" &&
          file.matrix[1].id == "target-profile";
        if (standard_case_matrix)
          components.push_back(case_slug);
        else
          for (const auto &axis : file.matrix)
            components.push_back(axis.id + "-" + combination_value(combination, axis.id));
        const auto path = runner::artifact_path(
          output_directory, components);
        try
          {
            auto scenario = make_b2_scenario(graetz_case);
            configure_b2_scenario(
              scenario,
              file,
              combination,
              std::string("chapter-6.b2.graetz-flow.") + case_slug);
            scenario.experiment.build_profile = NMOPT_COMPILED_BUILD_PROFILE;
            if (configuration.refinement_override.has_value())
              {
                if (scenario.compile.mesh.generation !=
                    nmopt::application::chapter6::MeshGeneration::
                      framework_native)
                  throw std::invalid_argument(
                    "--refinement cannot override a B2 simplex mesh; change its parameter file instead");
                scenario.compile.mesh.refinement =
                  *configuration.refinement_override;
                scenario.compile.mesh.mesh_provenance =
                  b2_mesh_provenance(*configuration.refinement_override);
              }
            scenario.experiment.harness.artifact_directory =
              output_directory.string();
            scenario.experiment.scenario_output_id =
              "chapter-6.b2.graetz-flow." + std::string(case_slug);

            nmopt::application::chapter6::dealii::B2ManufacturedDataT<2> data(
              graetz_case,
              scenario.problem.fixed_temperature,
              scenario.problem.forcing,
              scenario.problem.target_parameters);
            const auto runtime =
              nmopt::application::chapter6::dealii::
                make_b2_manufactured_runtime_data<2>(scenario, data);
            const auto session =
              nmopt::application::chapter6::dealii::
                make_b2_compilation_session<2>(scenario);
            const auto environment = make_environment(scenario);
            Adapter execute{runtime,
                            session,
                            environment,
                            path.parent_path() / "native"};
            Runner runner(scenario);
            const auto result = runner.run(
              [](const auto &parameters) {
                return chapter6::make_b2_problem_spec(parameters);
              },
              [&](const auto &specification, const auto &run_scenario) {
                auto evidence = execute(specification, run_scenario);
                add_b2_artifact_fields(evidence,
                                       run_scenario,
                                       configuration.framework_revision,
                                       configuration.run_kind);
                add_parameter_artifact_fields(evidence, configuration, combination);
                return evidence;
              });

            write_artifact(path, result.document);
            write_solver_trace(path.parent_path() / "solver-trace.csv",
                               result.artifact.envelope().report());
            run_manifest.record_success(path);
            std::cout << "B2 wrote " << path.string() << '\n';
          }
        catch (const std::exception &error)
          {
            all_artifacts_succeeded = false;
            run_manifest.record_failure(path, error.what());
            std::cerr << "B2 failed " << path.string() << ": "
                      << error.what() << '\n';
          }
      }
    return run_manifest.finalize() && all_artifacts_succeeded;
  }

  nmopt::application::runner::ParameterFile
  load_parameter_file(const nmopt::application::runner::CommandLineOptions &options)
  {
    std::filesystem::path path;
    if (options.parameter_file.has_value())
      path = nmopt::application::runner::find_file_from_current_or_parent(
        *options.parameter_file);
    else if (options.run_b1)
      path = nmopt::application::runner::find_file_from_current_or_parent(
        "parameters/chapter-6/b1/authoritative.prm");
    else
      path = nmopt::application::runner::find_file_from_current_or_parent(
        "parameters/chapter-6/b2/authoritative.prm");
    return nmopt::application::runner::read_parameter_file(path);
  }

  void
  copy_configuration_file(const std::filesystem::path &source,
                          const std::filesystem::path &destination)
  {
    nmopt::application::runner::prepare_artifact_path(destination);
    std::ifstream input(source);
    std::ofstream output(destination);
    if (!input || !output)
      throw std::runtime_error("could not snapshot configuration file '" +
                               source.string() + "'");
    output << input.rdbuf();
    if (!output)
      throw std::runtime_error("could not write configuration snapshot '" +
                               destination.string() + "'");
  }

  void
  snapshot_configuration(
    const nmopt::application::runner::ResolvedRunConfiguration &configuration,
    const nmopt::application::runner::ParameterFile &file,
    const std::vector<nmopt::application::runner::ParameterCombination> &combinations)
  {
    std::filesystem::create_directories(configuration.run_directory);
    copy_configuration_file(file.path,
                            configuration.run_directory / "parameters.prm");
    copy_configuration_file(configuration.plotting_profile_file,
                            configuration.run_directory / "plotting-profile.json");
    std::ofstream resolved(configuration.run_directory / "resolved-combinations.txt");
    if (!resolved)
      throw std::runtime_error("could not write resolved parameter combinations");
    resolved << join_combinations(combinations) << '\n';
  }

  struct PreparedRun
  {
    nmopt::application::runner::CommandLineOptions options;
    nmopt::application::runner::ParameterFile       file;
    nmopt::application::runner::ResolvedRunConfiguration configuration;
    std::vector<nmopt::application::runner::ParameterCombination> combinations;
  };

  PreparedRun
  prepare_run(const nmopt::application::runner::CommandLineOptions &input_options)
  {
    using namespace nmopt::application::runner;
    PreparedRun prepared{
      input_options,
      load_parameter_file(input_options),
      ResolvedRunConfiguration({},
                               {},
                               "",
                               "",
                               "",
                               RunKind::reproduction,
                               std::nullopt),
      {}};
    const auto benchmark_id = prepared.file.value("Benchmark/id");
    const bool is_b1 = benchmark_id == "chapter-6.b1.distributed-laplace";
    const bool is_b2 = benchmark_id == "chapter-6.b2.graetz-flow";
    if (!is_b1 && !is_b2)
      throw std::invalid_argument("parameter file declares unsupported benchmark '" +
                                  benchmark_id + "'");
    if (input_options.parameter_file.has_value())
      {
        prepared.options.run_b1 = is_b1;
        prepared.options.run_b2 = is_b2;
        prepared.options.run_kind =
          parse_run_kind(prepared.file.value("Run/kind"));
        if (!input_options.output_directory_explicit)
          prepared.options.output_directory = prepared.file.value("Run/output root");
        if (prepared.options.run_kind == RunKind::reproduction &&
            prepared.file.value("Run/build profile") != NMOPT_COMPILED_BUILD_PROFILE)
          throw std::invalid_argument(
            "parameter Run/build profile does not match the compiled runner profile");
      }
    else if ((input_options.run_b1 != is_b1) || (input_options.run_b2 != is_b2))
      throw std::invalid_argument("--benchmark does not match its authoritative parameter file");

    prepared.combinations =
      prepared.file.combinations(prepared.options.selection_filters);
    if (is_b1)
      {
        // This also validates the required B1 matrix axes before creating output.
        (void)b1_expected_artifacts(prepared.file,
                                    prepared.options.selection_filters);
      }
    else
      (void)b2_expected_artifacts(prepared.file,
                                  prepared.options.selection_filters);

    prepared.configuration = resolve_run_configuration(
      prepared.options, NMOPT_COMPILED_BUILD_PROFILE);
    prepared.configuration.parameter_file = prepared.file.path;
    prepared.configuration.parameter_hash = prepared.file.content_hash;
    const auto profile_path = find_file_from_current_or_parent(
      prepared.file.value("Postprocessing/style profile"));
    prepared.configuration.plotting_profile_file = profile_path;
    prepared.configuration.plotting_profile_hash = parameter_file_hash(profile_path);
    prepared.configuration.parameter_selection = join_selection(
      prepared.file, prepared.options.selection_filters);
    prepared.configuration.declared_matrix = join_matrix(prepared.file);
    prepared.configuration.excluded_combinations =
      join_combinations(prepared.file.excluded_combinations);
    prepared.configuration.resolved_combinations = join_combinations(prepared.combinations);
    prepared.configuration.comparison_rows =
      prepared.file.value("Postprocessing/comparison rows");
    prepared.configuration.comparison_columns =
      prepared.file.value("Postprocessing/comparison columns");
    prepared.configuration.comparison_group_by =
      prepared.file.value("Postprocessing/comparison group by");
    return prepared;
  }
} // namespace

int
main(const int argc, char **argv)
{
  try
    {
      const auto command = command_line_arguments(argc, argv);
      const auto options =
        nmopt::application::runner::parse_command_line(argc, argv);
      if (options.help)
        {
          print_usage(std::cout);
          return 0;
        }

      if (options.run_b1 || options.run_b2 || options.parameter_file.has_value())
        {
          auto prepared = prepare_run(options);
          snapshot_configuration(prepared.configuration,
                                 prepared.file,
                                 prepared.combinations);
          if (prepared.configuration.benchmark == "b1")
            return run_b1(prepared.configuration,
                          command,
                          prepared.file,
                          prepared.options.selection_filters)
                     ? 0
                     : 1;
          return run_b2(prepared.configuration,
                        command,
                        prepared.file,
                        prepared.options.selection_filters)
                   ? 0
                   : 1;
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
