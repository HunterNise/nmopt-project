#pragma once

#include "capability_registry.hpp"
#include "nmopt/application/chapter6.hpp"
#include "parameter_files.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>

namespace nmopt::application::runner::binding
{
  inline const std::string &
  combination_value(const ParameterCombination &combination,
                    const std::string_view   axis)
  {
    const auto found = combination.values.find(std::string(axis));
    if (found == combination.values.end())
      throw std::invalid_argument("resolved combination has no matrix axis '" +
                                  std::string(axis) + "'");
    return found->second;
  }

  inline nmopt::application::chapter6::ReducedMethod
  parse_method(const std::string &value)
  {
    return reduced_method_capability_registry().resolve(value,
                                                         "solver method");
  }

  inline nmopt::application::chapter6::ExecutionSelection
  parse_execution(const std::string &value)
  {
    return execution_capability_registry().resolve(value,
                                                    "compile execution");
  }

  inline nmopt::application::chapter6::ProductSelection
  parse_product(const std::string &value)
  {
    return product_capability_registry().resolve(value, "compile product");
  }

  inline nmopt::application::chapter6::MeshGeneration
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

  inline unsigned int
  parse_unsigned_text(const std::string &text, const std::string &key)
  {
    const auto value = parse_number_text(text, key);
    if (value < 0.0 || value > std::numeric_limits<unsigned int>::max() ||
        value != std::floor(value))
      throw std::invalid_argument("parameter '" + key +
                                  "' needs a nonnegative integer");
    return static_cast<unsigned int>(value);
  }

  inline nmopt::solvers::ReducedStoppingCriterion
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

  inline void
  require_parameter(const ParameterFile    &file,
                    const std::string_view key,
                    const std::string_view expected)
  {
    if (file.value(key) != expected)
      throw std::invalid_argument("parameter '" + std::string(key) +
                                  "' is '" + file.value(key) +
                                  "', but this adapter requires '" +
                                  std::string(expected) + "'");
  }

  template <typename Scenario>
  inline void
  apply_common_parameter_options(Scenario              &scenario,
                                 const ParameterFile   &file,
                                 const std::string     &scenario_id)
  {
    scenario.experiment.scenario_output_id = scenario_id;
    scenario.experiment.source_reference = file.value("Benchmark/source reference");
    scenario.experiment.source_revision = file.value("Benchmark/source revision");
    scenario.experiment.retain_fields = parameter_bool(file, "Output/retain fields");
    scenario.experiment.harness.measure_timings =
      parameter_bool(file, "Run/measure timings");
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
    scenario.compile.mesh.lower = parameter_finite_list(file, "Mesh/lower");
    scenario.compile.mesh.upper = parameter_finite_list(file, "Mesh/upper");
    scenario.compile.mesh.centroid_splits =
      parameter_unsigned(file, "Mesh/centroid splits");
    scenario.compile.mesh.selection_seed =
      parameter_unsigned(file, "Mesh/selection seed");
    scenario.compile.mesh.mesh_provenance = file.value("Mesh/provenance");
    scenario.compile.state_degree = parameter_unsigned(file, "Compile/state degree");
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
    scenario.compile.execution = parse_execution(file.value("Compile/execution"));
    scenario.compile.product = parse_product(file.value("Compile/product"));
  }

  inline void
  apply_solver_options(nmopt::application::chapter6::SolverOptions &solver,
                       const ParameterFile                     &file,
                       const std::string                        &method_id)
  {
    solver.globalization = reduced_globalization(file);
    const auto method_prefix = "Solver/method policy " + method_id + "/";
    const auto resolve = [&](const std::string_view entry) {
      return resolve_method_parameter(file, method_id, entry);
    };

    solver.initial_control_value =
      parameter_double(file, "Solver/initial independent control value");
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

    if (method_id == "l-bfgs")
      {
        const auto memory = file.optional_value(method_prefix + "memory");
        if (!memory.empty())
          solver.limited_memory_bfgs.memory_size =
            parse_unsigned_text(memory, method_prefix + "memory");
        const auto curvature =
          file.optional_value(method_prefix + "curvature tolerance");
        if (!curvature.empty())
          solver.limited_memory_bfgs.curvature_tolerance = parse_number_text(
            curvature, method_prefix + "curvature tolerance");
        const auto scaling = file.optional_value(
          method_prefix + "initial inverse Hessian scaling");
        if (!scaling.empty() && scaling != "metric-inverse" &&
            scaling != "scalar-secant")
          throw std::invalid_argument(
            "unknown L-BFGS initial scaling '" + scaling + "'");
        if (scaling == "metric-inverse")
          solver.limited_memory_bfgs.initial_inverse_hessian_scaling =
            nmopt::solvers::LimitedMemoryBfgsInitialScaling::metric_inverse;
        else if (scaling == "scalar-secant")
          solver.limited_memory_bfgs.initial_inverse_hessian_scaling =
            nmopt::solvers::LimitedMemoryBfgsInitialScaling::scalar_secant;
      }
    else if (method_id == "bfgs")
      {
        const auto key = method_prefix + "curvature tolerance";
        solver.full_bfgs.curvature_tolerance = parameter_double(file, key);
      }
  }
} // namespace nmopt::application::runner::binding
