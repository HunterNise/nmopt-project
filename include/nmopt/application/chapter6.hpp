#pragma once

#include "nmopt/application/chapter5.hpp"
#include "nmopt/application/scenario.hpp"
#include "nmopt/solvers/reduced_search.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nmopt::application::chapter6
{
  inline constexpr const char *chapter6_numerical_examples_source_revision =
    "1aaefbe473f9941a89d1df36192251511c052933";

  inline constexpr const char *b1_recipe_id =
    chapter5::scalar_distributed_recipe_id;
  inline constexpr const char *b2_recipe_id =
    chapter5::neumann_convection_recipe_id;

  inline constexpr unsigned int b2_fixed_boundary_id = 0;
  inline constexpr unsigned int b2_control_boundary_id = 1;
  inline constexpr unsigned int b2_outflow_boundary_id = 2;

  inline constexpr const char *b2_fixed_boundary_region_id =
    "dirichlet_boundary";
  inline constexpr const char *b2_control_boundary_region_id =
    "control_boundary";
  inline constexpr const char *b2_outflow_boundary_region_id =
    "outflow_boundary";

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

  inline constexpr std::array<GraetzCase, 4> b2_case_order = {
    {GraetzCase::observation_wings_constant_target,
     GraetzCase::observation_full_constant_target,
     GraetzCase::observation_wings_parabolic_target,
     GraetzCase::observation_full_parabolic_target}};

  inline const char *
  graetz_case_name(const GraetzCase graetz_case)
  {
    switch (graetz_case)
      {
        case GraetzCase::observation_wings_constant_target:
          return "wings-constant";
        case GraetzCase::observation_full_constant_target:
          return "full-constant";
        case GraetzCase::observation_wings_parabolic_target:
          return "wings-parabolic";
        case GraetzCase::observation_full_parabolic_target:
          return "full-parabolic";
        default:
          throw std::invalid_argument("B2 has an unknown Graetz case");
      }
  }

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

  inline void
  add_b2_fixed_dirichlet_reconstruction(
    semantic::v1::ProblemSpec &specification)
  {
    specification.data.push_back(
      {"fixed_dirichlet_data",
       "Fixed Dirichlet data",
       semantic::v1::DataKind::function,
       semantic::v1::DataRole::fixed_dirichlet_lifting,
       "state_space"});
    specification.transformations = {
      {"fixed_dirichlet_reconstruction",
       "Fixed Dirichlet reconstruction",
       semantic::v1::TransformationKind::fixed_dirichlet_reconstruction,
       "state",
       "state_space",
       "fixed_dirichlet_data",
       ""}};
    std::find_if(specification.variables.begin(),
                 specification.variables.end(),
                 [](const auto &variable) { return variable.id == "state"; })
      ->physical_field_transform_id = "fixed_dirichlet_reconstruction";
    std::find_if(specification.requirement_policies.begin(),
                 specification.requirement_policies.end(),
                 [](const auto &policy) {
                   return policy.id == "state_fixed_dirichlet";
                 })
      ->selected_policy =
      "P_h independent FE_Q coordinates plus nodal fixed Dirichlet lifting";
  }

  inline semantic::v1::ProblemSpec
  make_b2_problem_spec(const B2ProblemParameters &parameters)
  {
    auto specification =
      chapter5::make_neumann_convection_recipe()(parameters.recipe);

    const auto partition_policy = std::find_if(
      specification.requirement_policies.begin(),
      specification.requirement_policies.end(),
      [](const semantic::v1::RequirementPolicySpec &policy) {
        return policy.id == "neumann_convection_partition";
      });
    if (partition_policy == specification.requirement_policies.end())
      throw std::invalid_argument(
        "B2 requires the registered Neumann-convection boundary partition");

    specification.regions.push_back(
      {b2_outflow_boundary_region_id,
       "Natural transport outflow boundary",
       semantic::v1::RegionKind::boundary,
       false,
       {b2_outflow_boundary_id},
       {},
       {}});
    partition_policy->selected_policy =
      "disjoint complete fixed-Dirichlet, Neumann-control, and natural outflow boundary regions; the Neumann datum is the conormal flux of the conservative transport form";
    partition_policy->typed_selection = semantic::v1::BoundaryRealisationSelection{
      partition_policy->id,
      "state",
      b2_fixed_boundary_region_id,
      b2_outflow_boundary_region_id,
      {},
      {},
      b2_outflow_boundary_region_id,
      semantic::v1::ConormalForm::diffusion_minus_transport,
      semantic::v1::NormalOrientation::outward,
      semantic::v1::TraceEvaluationRealisation::fe_q_state_trace,
      semantic::v1::FaceQuadratureRealisation::qgauss_face};
    add_b2_fixed_dirichlet_reconstruction(specification);
    return specification;
  }

  inline semantic::v1::ProblemSpec
  make_b2_problem_spec(const B2Scenario &scenario)
  {
    validate_b2(scenario);
    return make_b2_problem_spec(scenario.problem);
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
    solver.parameters.gradient_tolerance =
      method == ReducedMethod::steepest_descent ? 1.0e-3 : 1.0e-8;
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
       chapter6_numerical_examples_source_revision,
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

    std::string scenario_id = "chapter-6.b2.graetz-flow";
    if (graetz_case != GraetzCase::observation_wings_constant_target)
      scenario_id += "." + std::string(graetz_case_name(graetz_case));

    B2Scenario scenario{
      {scenario_id,
       "E6.5.2 Graetz-flow boundary control",
       "Facewise Neumann control with downstream subdomain tracking",
       "chapter-6",
       b2_recipe_id,
       {"B0 harness", "fixed temperature port", "declared Galerkin policy"}},
      std::move(problem),
      {},
      {ReducedMethod::bfgs, {}, std::nullopt, "zero"},
      {scenario_id,
       "E6.5.2 equation (6.65), Table 6.2, Figures 6.4-6.5",
       chapter6_numerical_examples_source_revision,
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
    catalog.add(b1.metadata);
    for (const auto graetz_case : b2_case_order)
      catalog.add(make_b2_scenario(graetz_case).metadata);
    return catalog;
  }
} // namespace nmopt::application::chapter6
