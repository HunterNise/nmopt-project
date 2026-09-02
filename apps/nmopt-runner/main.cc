#include "nmopt/application/application.hpp"
#include "nmopt/application/dealii/chapter6_b1.hpp"
#include "nmopt/application/dealii/chapter6_b2.hpp"
#include "benchmark_registry.hpp"
#include "benchmark_binders.hpp"
#include "parameter_binding.hpp"
#include "parameter_files.hpp"
#include "runner.hpp"

#include <array>
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
  using nmopt::application::runner::binding::combination_value;
  using nmopt::application::runner::binding::apply_common_parameter_options;
  using nmopt::application::runner::binding::apply_solver_options;
  using nmopt::application::runner::binding::bind_b1_scenario;
  using nmopt::application::runner::binding::bind_b2_scenario;
  using nmopt::application::runner::binding::parse_mesh_generation;
  using nmopt::application::runner::binding::parse_method;
  using nmopt::application::runner::parse_number_text;
  using nmopt::application::runner::binding::parse_stopping_criterion;
  using nmopt::application::runner::binding::parse_unsigned_text;
  using nmopt::application::runner::binding::require_parameter;

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
           << "       nmopt_runner --benchmark ID --framework-revision REV"
              " [--output DIRECTORY] [--run-kind KIND] [--run-slot SLOT]"
              " [--refinement N]\n"
           << "       nmopt_runner --help\n"
           << "\n"
           << "--list             list registered Chapter 5/6 application entries\n"
           << "--benchmark ID     run a registered benchmark by identifier\n"
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

  std::string
  b1_beta_coordinate(const std::string &value_text)
  {
    parse_number_text(value_text, "Matrix/regularisation");
    return value_text;
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

  std::string
  join_coordinates(const std::vector<double> &coordinates)
  {
    std::ostringstream output;
    for (std::size_t index = 0; index < coordinates.size(); ++index)
      {
        if (index != 0)
          output << ',';
        output << b1_number(coordinates[index]);
      }
    return output.str();
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
  join_combinations(const nmopt::application::runner::RunSetPlan &plan)
  {
    std::ostringstream output;
    for (std::size_t index = 0; index < plan.resolved_combinations.size();
         ++index)
      {
        if (index != 0)
          output << ';';
        output << '[';
        std::size_t value_index = 0;
        for (const auto &[axis, value] :
             plan.resolved_combinations[index].values.values)
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
  b1_artifact_coordinate_components(
    const nmopt::application::runner::RunSetPlan &,
    const nmopt::application::runner::RunSetCombination &combination)
  {
    const auto method =
      parse_method(combination_value(combination.values, "method"));
    const auto &beta_text = combination_value(combination.values,
                                              "regularisation");
    return {b1_method_slug(method), "beta-" + b1_beta_coordinate(beta_text)};
  }

  std::vector<std::string>
  b2_artifact_coordinate_components(
    const nmopt::application::runner::RunSetPlan &plan,
    const nmopt::application::runner::RunSetCombination &combination)
  {
    std::vector<std::string> axis_ids;
    axis_ids.reserve(plan.matrix_axes.size());
    for (const auto &axis : plan.matrix_axes)
      axis_ids.push_back(axis.id);
    std::sort(axis_ids.begin(), axis_ids.end());
    const bool standard_case_axes =
      axis_ids == std::vector<std::string>{"observation-region", "target-profile"};
    if (standard_case_axes)
      return {combination_value(combination.values, "observation-region") + "-" +
              combination_value(combination.values, "target-profile")};

    std::vector<std::string> components;
    components.reserve(combination.artifact_coordinates.size());
    for (const auto &coordinate : combination.artifact_coordinates)
      components.push_back(coordinate.axis + "-" + coordinate.value);
    return components;
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
      {"benchmark.mesh_lower", join_coordinates(scenario.compile.mesh.lower)});
    evidence.fields.push_back(
      {"benchmark.mesh_upper", join_coordinates(scenario.compile.mesh.upper)});
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
    const auto case_slug = nmopt::application::chapter6::b2_case_name(
      scenario.problem.observation_region, scenario.problem.target_profile);
    const auto &manifest = evidence.envelope.compilation_manifest();
    const auto &volume_observation = *scenario.compile.volume_observation;
    const std::string volume_observation_target =
      nmopt::application::chapter6::
        volume_observation_target_realisation_name(
          volume_observation.target_realisation);
    const std::string quadrature_order =
      std::to_string(volume_observation.quadrature_order);
    const auto &observation_definition =
      nmopt::application::selected_scalar_function_definition(
        scenario.problem.observation_region_catalog, "B2 observation catalog");
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
    evidence.fields.push_back({"benchmark.observation_definition",
                               observation_definition.id});
    evidence.fields.push_back(
      {"benchmark.observation_kind",
       nmopt::application::scalar_function_kind_name(
         observation_definition.kind)});
    evidence.fields.push_back(
      {"benchmark.observation_value", b1_number(observation_definition.value)});
    evidence.fields.push_back({"benchmark.observation_expression",
                               observation_definition.expression});
    evidence.fields.push_back({"benchmark.observation_provenance",
                               observation_definition.provenance});
    evidence.fields.push_back({"benchmark.observation_realisation",
                               "cell-center-indicator"});
    evidence.fields.push_back(
      {"benchmark.observed_material_id",
       std::to_string(scenario.problem.recipe.observed_material_id)});
    evidence.fields.push_back(
      {"benchmark.fixed_boundary_id",
       std::to_string(scenario.problem.boundary.fixed_id)});
    evidence.fields.push_back(
      {"benchmark.control_boundary_id",
       std::to_string(scenario.problem.boundary.control_id)});
    evidence.fields.push_back(
      {"benchmark.outflow_boundary_id",
       std::to_string(scenario.problem.boundary.outflow_id)});
    evidence.fields.push_back(
      {"benchmark.upstream_transition",
       b1_number(scenario.problem.boundary.upstream_transition)});
    evidence.fields.push_back({"benchmark.fixed_dirichlet_data",
                               scenario.problem.fixed_dirichlet_data.id});
    evidence.fields.push_back(
      {"benchmark.fixed_dirichlet_data_kind",
       nmopt::application::scalar_function_kind_name(
         scenario.problem.fixed_dirichlet_data.kind)});
    evidence.fields.push_back(
      {"benchmark.fixed_dirichlet_data_value",
       b1_number(scenario.problem.fixed_dirichlet_data.value)});
    evidence.fields.push_back(
      {"benchmark.fixed_dirichlet_data_expression",
       scenario.problem.fixed_dirichlet_data.expression});
    evidence.fields.push_back({"benchmark.conservative_transport",
                               scenario.problem.conservative_transport.id});
    evidence.fields.push_back(
      {"benchmark.conservative_transport_expression",
       scenario.problem.conservative_transport.expression});
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
    evidence.fields.push_back({"provenance.fixed_dirichlet_data",
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
         const std::vector<std::string> &,
         const nmopt::application::runner::RunSetPlan &plan,
         const nmopt::application::runner::ParameterFile &file,
         nmopt::application::runner::RunSetManifest &run_manifest)
  {
    using namespace nmopt::application;
    using namespace chapter6;
    using Runner = benchmark::HeadlessBenchmarkRunnerT<B1Scenario>;
    using Adapter =
      nmopt::application::chapter6::dealii::B1ReducedExecutionAdapterT<2>;
    const auto &output_directory = configuration.run_directory;
    bool all_artifacts_succeeded = true;

    struct ObjectiveReference
    {
      double      value;
      std::string artifact;
    };
    std::map<std::string, ObjectiveReference> objective_references;
    auto combinations = plan.resolved_combinations;
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
              parse_method(combination_value(left.values, "method")) ==
              reference_method;
            const bool right_is_reference =
              parse_method(combination_value(right.values, "method")) ==
              reference_method;
            return left_is_reference && !right_is_reference;
          });
      }

    for (const auto &planned_combination : combinations)
      {
        const auto &combination = planned_combination.values;
        const auto method_id = combination_value(combination, "method");
        const auto method = parse_method(method_id);
        const auto &beta_text = combination_value(combination,
                                                  "regularisation");
        const auto beta = parse_number_text(beta_text, "Matrix/regularisation");
        const auto beta_coordinate = b1_beta_coordinate(beta_text);
        const auto path = runner::artifact_path(
          output_directory,
          runner::run_set_artifact_components(
            plan, planned_combination, b1_artifact_coordinate_components));
        try
          {
            auto scenario = make_b1_scenario(method);
            bind_b1_scenario(
              scenario,
              file,
              combination,
              method_id,
              beta,
              std::string(file.value("Benchmark/id")) + "." + method_id +
                ".beta-" + beta_coordinate);
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
                const auto reference = objective_references.find(beta_coordinate);
                if (reference == objective_references.end())
                  throw std::invalid_argument(
                    "objective-target reference method did not reach its stopping tolerance for beta " +
                    beta_coordinate);
                scenario.solver.parameters.objective_target =
                  reference->second.value;
                objective_target_reference_artifact =
                  reference->second.artifact;
              }

            nmopt::application::chapter6::dealii::B1SelectedDataT<2> data(
              scenario.problem.forcing, scenario.problem.desired_state);
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
              objective_references[beta_coordinate] = {
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
    return all_artifacts_succeeded;
  }

  bool
  run_b2(const nmopt::application::runner::ResolvedRunConfiguration &configuration,
         const std::vector<std::string> &,
         const nmopt::application::runner::RunSetPlan &plan,
         const nmopt::application::runner::ParameterFile &file,
         nmopt::application::runner::RunSetManifest &run_manifest)
  {
    using namespace nmopt::application;
    using namespace chapter6;
    using Runner = benchmark::HeadlessBenchmarkRunnerT<B2Scenario>;
    using Adapter =
      nmopt::application::chapter6::dealii::B2ReducedExecutionAdapterT<2>;
    const auto &output_directory = configuration.run_directory;
    bool all_artifacts_succeeded = true;

    for (const auto &planned_combination : plan.resolved_combinations)
      {
        const auto &combination = planned_combination.values;
        const auto &observation_region =
          combination_value(combination, "observation-region");
        const auto target_profile = combination_value(combination,
                                                       "target-profile");
        const auto case_slug = b2_case_name(observation_region, target_profile);
        const auto path = runner::artifact_path(
          output_directory,
          runner::run_set_artifact_components(
            plan, planned_combination, b2_artifact_coordinate_components));
        try
          {
            const auto target_catalog = runner::b2_target_catalog(
              file, target_profile);
            auto scenario = make_b2_scenario_with_target_catalog(
              observation_region, target_profile, target_catalog);
            bind_b2_scenario(
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
              selected_scalar_function_definition(
                scenario.problem.observation_region_catalog),
              scenario.problem.fixed_dirichlet_data,
              scenario.problem.forcing,
              chapter6::b2_target_definition(scenario.problem.target_catalog),
              scenario.problem.conservative_transport);
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
    return all_artifacts_succeeded;
  }

  const std::array<nmopt::application::runner::BenchmarkExecutionRegistration, 2> &
  benchmark_execution_registrations()
  {
    using Registration =
      nmopt::application::runner::BenchmarkExecutionRegistration;
    static const std::array<Registration, 2> registrations{{
      {nmopt::application::runner::find_benchmark_registration("b1"),
       [](const auto &plan) {
         return nmopt::application::runner::run_set_artifact_paths(
           plan, b1_artifact_coordinate_components);
       },
       run_b1},
      {nmopt::application::runner::find_benchmark_registration("b2"),
       [](const auto &plan) {
         return nmopt::application::runner::run_set_artifact_paths(
           plan, b2_artifact_coordinate_components);
       },
       run_b2}}};
    return registrations;
  }

  const nmopt::application::runner::BenchmarkExecutionRegistration *
  find_benchmark_execution_registration(const std::string_view id)
  {
    for (const auto &registration : benchmark_execution_registrations())
      if (registration.metadata != nullptr && registration.metadata->id == id)
        return &registration;
    return nullptr;
  }

  nmopt::application::runner::ParameterFile
  load_parameter_file(const nmopt::application::runner::CommandLineOptions &options)
  {
    std::filesystem::path path;
    if (options.parameter_file.has_value())
      path = nmopt::application::runner::find_file_from_current_or_parent(
        *options.parameter_file);
    else
      {
        const auto *registration =
          nmopt::application::runner::find_benchmark_registration(
            options.benchmark.value_or(""));
        if (registration == nullptr)
          throw std::invalid_argument(
            "unsupported benchmark '" + options.benchmark.value_or("") +
            "'; available benchmark IDs: " +
            nmopt::application::runner::registered_benchmark_ids());
        path = nmopt::application::runner::find_file_from_current_or_parent(
          registration->default_parameter_file);
      }
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
    const nmopt::application::runner::RunSetPlan &plan)
  {
    std::filesystem::create_directories(configuration.run_directory);
    copy_configuration_file(file.path,
                            configuration.run_directory / "parameters.prm");
    copy_configuration_file(configuration.plotting_profile_file,
                            configuration.run_directory / "plotting-profile.json");
    std::ofstream resolved(configuration.run_directory / "resolved-combinations.txt");
    if (!resolved)
      throw std::runtime_error("could not write resolved parameter combinations");
    resolved << join_combinations(plan) << '\n';
  }

  struct PreparedRun
  {
    nmopt::application::runner::CommandLineOptions options;
    nmopt::application::runner::ParameterFile       file;
    nmopt::application::runner::ResolvedRunConfiguration configuration;
    nmopt::application::runner::RunSetPlan              plan;
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
    const auto *registration =
      nmopt::application::runner::find_benchmark_registration_for_parameter_id(
        benchmark_id);
    if (registration == nullptr)
      throw std::invalid_argument("parameter file declares unsupported benchmark '" +
                                  benchmark_id + "'; available benchmark IDs: " +
                                  nmopt::application::runner::registered_benchmark_ids());
    if (input_options.parameter_file.has_value())
      {
        prepared.options.benchmark = std::string(registration->id);
        prepared.options.run_kind =
          parse_run_kind(prepared.file.value("Run/kind"));
        if (!input_options.output_directory_explicit)
          prepared.options.output_directory = prepared.file.value("Run/output root");
        if (prepared.options.run_kind == RunKind::reproduction &&
            prepared.file.value("Run/build profile") != NMOPT_COMPILED_BUILD_PROFILE)
          throw std::invalid_argument(
            "parameter Run/build profile does not match the compiled runner profile");
      }
    else if (!input_options.benchmark.has_value() ||
             *input_options.benchmark != registration->id)
      throw std::invalid_argument("--benchmark does not match its authoritative parameter file");

    prepared.plan = make_run_set_plan(prepared.file,
                                      prepared.options.selection_filters);
    const auto *execution_registration =
      find_benchmark_execution_registration(registration->id);
    if (execution_registration == nullptr)
      throw std::invalid_argument(
        "benchmark registration has no execution callback: " +
        std::string(registration->id));
    // This also validates the required matrix axes before creating output.
    (void)execution_registration->artifact_planner(prepared.plan);

    prepared.configuration = resolve_run_configuration(
      prepared.options, NMOPT_COMPILED_BUILD_PROFILE);
    prepared.configuration.parameter_file =
      prepared.plan.parameter_provenance.file;
    prepared.configuration.parameter_hash =
      prepared.plan.parameter_provenance.content_hash;
    const auto profile_path = find_file_from_current_or_parent(
      prepared.file.value("Postprocessing/style profile"));
    prepared.configuration.plotting_profile_file = profile_path;
    prepared.configuration.plotting_profile_hash = parameter_file_hash(profile_path);
    prepared.configuration.parameter_selection = join_selection(
      prepared.file, prepared.options.selection_filters);
    prepared.configuration.declared_matrix = join_matrix(prepared.file);
    prepared.configuration.excluded_combinations =
      join_combinations(prepared.plan.excluded_combinations);
    prepared.configuration.resolved_combinations = join_combinations(prepared.plan);
    prepared.configuration.comparison_rows = prepared.plan.comparison.rows;
    prepared.configuration.comparison_columns = prepared.plan.comparison.columns;
    prepared.configuration.comparison_group_by = prepared.plan.comparison.group_by;
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

      if (options.benchmark.has_value() || options.parameter_file.has_value())
        {
          auto prepared = prepare_run(options);
          snapshot_configuration(prepared.configuration,
                                 prepared.file,
                                 prepared.plan);
          const auto *registration =
            find_benchmark_execution_registration(
              prepared.configuration.benchmark);
          if (registration == nullptr)
            throw std::invalid_argument(
              "benchmark has no execution registration: " +
              prepared.configuration.benchmark);
          nmopt::application::runner::RunSetManifest run_manifest(
            prepared.configuration,
            command,
            registration->artifact_planner(prepared.plan));
          const bool execution_succeeded = registration->execute(
            prepared.configuration,
            command,
            prepared.plan,
            prepared.file,
            run_manifest);
          return run_manifest.finalize() && execution_succeeded ? 0 : 1;
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
