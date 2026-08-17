#pragma once

#include "nmopt/application/chapter5.hpp"
#include "nmopt/application/scenario.hpp"
#include "nmopt/solvers/reduced_search.hpp"

#include <cmath>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nmopt::application::chapter6
{
  inline constexpr const char *b1_recipe_id =
    chapter5::scalar_distributed_recipe_id;
  inline constexpr const char *b2_recipe_id =
    chapter5::neumann_convection_recipe_id;

  enum class ProductSelection
  {
    reduced_dto,
    quadratic_kkt,
    pdas
  };

  enum class ExecutionSelection
  {
    assembled,
    matrix_free
  };

  // The adapter that executes a scenario maps these values to the concrete
  // compiler::v1 and solver policy types. Keeping the records backend-neutral
  // lets discovery and scenario validation run without deal.II.
  struct MeshOptions
  {
    unsigned int dimension = 2;
    unsigned int refinement = 0;
    std::string  mesh_provenance = "scenario-owned-mesh";
  };

  struct CompileOptions
  {
    MeshOptions       mesh;
    unsigned int     state_degree = 1;
    ExecutionSelection execution = ExecutionSelection::assembled;
    ProductSelection  product = ProductSelection::reduced_dto;
    bool              owned_session = true;
  };

  enum class ReducedMethod
  {
    steepest_descent,
    bfgs,
    limited_memory_bfgs
  };

  struct SolverOptions
  {
    ReducedMethod                  method = ReducedMethod::steepest_descent;
    solvers::ReducedSolverParameters parameters;
    // The current generic solver policies do not all expose a minimum-step
    // field. Preserve a source-declared value here so an execution adapter
    // cannot silently discard it when reproducing a benchmark.
    std::optional<double> declared_minimum_step_length = std::nullopt;
    std::string           initial_control = "zero";
  };

  struct ExperimentOptions
  {
    std::string scenario_output_id;
    std::string source_reference;
    std::string source_revision = "to-be-recorded";
    std::string build_profile = "release-dealii";
    bool        retain_fields = true;
    struct Harness
    {
      bool        deterministic = true;
      bool        serialize_artifacts = true;
      bool        measure_timings = true;
      bool        measure_memory = false;
      std::string artifact_directory = "runs";
    } harness;
  };

  struct ScalarRuntimeDataOptions
  {
    double      diffusion = 1.0;
    double      reaction = 0.0;
    double      regularisation_weight = 1.0e-2;
    std::string forcing_provenance;
    std::string desired_state_provenance;
    std::string fixed_dirichlet_data_provenance;
    std::string conservative_transport_provenance;
  };

  enum class B1ForcingSelection
  {
    recovered_source,
    manufactured_zero
  };

  struct B1ProblemParameters
  {
    chapter5::ScalarDistributedControlParameters recipe;
    ScalarRuntimeDataOptions                     data{
      1.0,
      0.0,
      1.0e-2,
      "chapter-6.e6.5.1.forcing",
      "chapter-6.e6.5.1.desired-state",
      "",
      ""};
    std::vector<double> regularisation_sweep = {1.0e-1, 1.0e-2, 1.0e-3};
    B1ForcingSelection forcing_selection = B1ForcingSelection::manufactured_zero;
  };

  enum class GraetzCase
  {
    observation_wings_constant_target,
    observation_full_constant_target,
    observation_wings_parabolic_target,
    observation_full_parabolic_target
  };

  struct B2ProblemParameters
  {
    chapter5::NeumannConvectionParameters recipe;
    ScalarRuntimeDataOptions              data{
      0.1,
      0.0,
      1.0e-3,
      "chapter-6.e6.5.2.zero-forcing",
      "chapter-6.e6.5.2.target",
      "chapter-6.e6.5.2.fixed-temperature",
      "chapter-6.e6.5.2.graetz-transport"};
    GraetzCase graetz_case = GraetzCase::observation_wings_constant_target;
    double     fixed_temperature = 1.0;
    // The current semantic graph seed is homogeneous. The execution adapter
    // must add and validate the explicit fixed-data graph port before binding
    // this nonzero source boundary value.
  };

  using B1Scenario =
    ScenarioT<B1ProblemParameters,
              CompileOptions,
              SolverOptions,
              ExperimentOptions>;
  using B2Scenario =
    ScenarioT<B2ProblemParameters,
              CompileOptions,
              SolverOptions,
              ExperimentOptions>;

  inline void
  validate_common_compile_options(const CompileOptions &options)
  {
    if (options.mesh.dimension == 0)
      throw std::invalid_argument("benchmark scenarios need a mesh dimension");
    if (options.state_degree == 0)
      throw std::invalid_argument("benchmark scenarios need a positive state degree");
    if (options.mesh.mesh_provenance.empty())
      throw std::invalid_argument(
        "benchmark scenarios need mesh provenance");
    if (!options.owned_session)
      throw std::invalid_argument(
        "selected benchmark scenarios require an owned compilation session");
  }

  inline void
  validate_runtime_data(const ScalarRuntimeDataOptions &data)
  {
    if (!std::isfinite(data.diffusion) || data.diffusion <= 0.0)
      throw std::invalid_argument(
        "benchmark scenarios need positive finite diffusion");
    if (!std::isfinite(data.reaction) || data.reaction < 0.0)
      throw std::invalid_argument(
        "benchmark scenarios need nonnegative finite reaction");
    if (!std::isfinite(data.regularisation_weight) ||
        data.regularisation_weight <= 0.0)
      throw std::invalid_argument(
        "benchmark scenarios need positive finite regularisation");
    if (data.forcing_provenance.empty() ||
        data.desired_state_provenance.empty())
      throw std::invalid_argument(
        "benchmark scenarios need forcing and target provenance");
  }

  inline void
  validate_b1(const B1Scenario &scenario)
  {
    scenario.validate();
    validate_common_compile_options(scenario.compile);
    validate_runtime_data(scenario.problem.data);
    if (scenario.problem.regularisation_sweep.empty())
      throw std::invalid_argument(
        "B1 needs at least one regularisation value");
    for (const double value : scenario.problem.regularisation_sweep)
      if (!std::isfinite(value) || value <= 0.0)
        throw std::invalid_argument(
          "B1 regularisation sweep values must be positive and finite");
    if (scenario.problem.data.forcing_provenance.empty())
      throw std::invalid_argument("B1 needs forcing provenance");
    switch (scenario.problem.forcing_selection)
      {
        case B1ForcingSelection::recovered_source:
        case B1ForcingSelection::manufactured_zero:
          break;
        default:
          throw std::invalid_argument("B1 has an unknown forcing selection");
      }
    if (scenario.solver.method != ReducedMethod::steepest_descent &&
        scenario.solver.method != ReducedMethod::limited_memory_bfgs)
      throw std::invalid_argument(
        "B1 selects steepest descent or limited-memory BFGS");
    if (scenario.compile.product != ProductSelection::reduced_dto ||
        scenario.compile.execution != ExecutionSelection::assembled)
      throw std::invalid_argument(
        "B1 freezes the assembled reduced-DTO product");
  }

  inline semantic::v1::ProblemSpec
  make_b1_problem_spec(const B1ProblemParameters &parameters)
  {
    return chapter5::make_scalar_distributed_recipe()(parameters.recipe);
  }

  inline semantic::v1::ProblemSpec
  make_b1_problem_spec(const B1Scenario &scenario)
  {
    validate_b1(scenario);
    return make_b1_problem_spec(scenario.problem);
  }

  inline void
  validate_b2(const B2Scenario &scenario)
  {
    scenario.validate();
    validate_common_compile_options(scenario.compile);
    validate_runtime_data(scenario.problem.data);
    if (!std::isfinite(scenario.problem.fixed_temperature))
      throw std::invalid_argument("B2 fixed temperature must be finite");
    if (scenario.problem.data.fixed_dirichlet_data_provenance.empty() ||
        scenario.problem.data.conservative_transport_provenance.empty())
      throw std::invalid_argument(
        "B2 needs fixed-data and transport provenance");
    if (scenario.solver.method != ReducedMethod::bfgs)
      throw std::invalid_argument("B2 selects the BFGS method");
    if (scenario.compile.product != ProductSelection::reduced_dto ||
        scenario.compile.execution != ExecutionSelection::assembled)
      throw std::invalid_argument(
        "B2 freezes the assembled reduced-DTO product");
  }

  inline B1Scenario
  make_b1_scenario(
    const ReducedMethod       method = ReducedMethod::steepest_descent,
    const B1ForcingSelection forcing = B1ForcingSelection::manufactured_zero)
  {
    B1ProblemParameters problem;
    problem.forcing_selection = forcing;
    problem.data.forcing_provenance =
      forcing == B1ForcingSelection::manufactured_zero
        ? "chapter-6.e6.5.1.manufactured-zero-forcing"
        : "chapter-6.e6.5.1.recovered-forcing";

    SolverOptions solver;
    solver.method = method;
    solver.parameters.maximum_line_search_trials = 5;
    solver.parameters.armijo_fraction = 1.0e-5;
    solver.parameters.backtracking_factor = 0.7;
    solver.declared_minimum_step_length =
      method == ReducedMethod::steepest_descent ? 0.2 : 0.01;

    B1Scenario scenario{
      {"chapter-6.b1.distributed-laplace",
       "E6.5.1 distributed Laplace control",
       "Reduced-space comparison of steepest descent and L-BFGS",
       "chapter-6",
       b1_recipe_id,
       {"B0 harness", "manufactured or recovered forcing is recorded"}},
      std::move(problem),
      {},
      std::move(solver),
      {"chapter-6.b1.distributed-laplace",
       "E6.5.1 equations (6.64), Figures 6.2-6.3",
       "to-be-recorded",
       "release-dealii",
       true,
       {true, true, true, false, "runs"}}};
    validate_b1(scenario);
    return scenario;
  }

  inline B2Scenario
  make_b2_scenario(
    const GraetzCase graetz_case =
      GraetzCase::observation_wings_constant_target)
  {
    B2ProblemParameters problem;
    problem.graetz_case = graetz_case;

    B2Scenario scenario{
      {"chapter-6.b2.graetz-flow",
       "E6.5.2 Graetz-flow boundary control",
       "Facewise Neumann control with downstream subdomain tracking",
       "chapter-6",
       b2_recipe_id,
       {"B0 harness", "fixed temperature port", "declared Galerkin policy"}},
      std::move(problem),
      {},
      {ReducedMethod::bfgs, {}, std::nullopt, "zero"},
      {"chapter-6.b2.graetz-flow",
       "E6.5.2 equation (6.65), Table 6.2, Figures 6.4-6.5",
       "to-be-recorded",
       "release-dealii",
       true,
       {true, true, true, false, "runs"}}};
    validate_b2(scenario);
    return scenario;
  }

  inline ApplicationCatalog
  make_catalog()
  {
    ApplicationCatalog catalog;
    const auto b1 = make_b1_scenario();
    const auto b2 = make_b2_scenario();
    catalog.add(b1.metadata);
    catalog.add(b2.metadata);
    return catalog;
  }
} // namespace nmopt::application::chapter6
