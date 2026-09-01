#pragma once

#include "nmopt/application/chapter5.hpp"
#include "nmopt/application/scalar_function.hpp"
#include "nmopt/application/scenario.hpp"
#include "nmopt/solvers/reduced_search.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
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

  inline constexpr unsigned int b1_default_mesh_refinement = 7;
  inline constexpr unsigned int b2_default_mesh_refinement = 7;

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

  enum class MeshGeneration
  {
    framework_native,
    structured_simplex,
    centroid_split_simplex
  };

  inline const char *
  mesh_generation_name(const MeshGeneration generation)
  {
    switch (generation)
      {
        case MeshGeneration::framework_native:
          return "framework-native";
        case MeshGeneration::structured_simplex:
          return "structured-simplex";
        case MeshGeneration::centroid_split_simplex:
          return "centroid-split-simplex";
      }
    throw std::invalid_argument("unknown mesh generation selection");
  }

  // The adapter that executes a scenario maps these values to the concrete
  // compiler::v1 and solver policy types. Keeping the records backend-neutral
  // lets discovery and scenario validation run without deal.II.
  struct MeshOptions
  {
    MeshGeneration            generation = MeshGeneration::framework_native;
    unsigned int              dimension = 2;
    unsigned int              refinement = 0;
    unsigned int              subdivisions = 0;
    std::vector<unsigned int> axis_subdivisions;
    unsigned int              centroid_splits = 0;
    unsigned int              selection_seed = 0;
    std::string               mesh_provenance = "scenario-owned-mesh";
  };

  struct IterativeSolveOptions
  {
    unsigned int maximum_iterations = 0;
    double       relative_tolerance = 1.0e-12;
    double       absolute_tolerance = 1.0e-14;
  };

  enum class VolumeObservationTargetRealisation
  {
    analytic_quadrature,
    state_fe_interpolation
  };

  inline const char *
  volume_observation_target_realisation_name(
    const VolumeObservationTargetRealisation realisation)
  {
    switch (realisation)
      {
        case VolumeObservationTargetRealisation::analytic_quadrature:
          return "analytic-quadrature";
        case VolumeObservationTargetRealisation::state_fe_interpolation:
          return "state-fe-interpolation";
      }
    throw std::invalid_argument(
      "unknown volume-observation target realisation");
  }

  struct VolumeObservationOptions
  {
    unsigned int quadrature_order = 3;
    VolumeObservationTargetRealisation target_realisation =
      VolumeObservationTargetRealisation::analytic_quadrature;
  };

  struct CompileOptions
  {
    MeshOptions           mesh;
    unsigned int          state_degree = 1;
    ExecutionSelection    execution = ExecutionSelection::assembled;
    ProductSelection      product = ProductSelection::reduced_dto;
    bool                  owned_session = true;
    std::optional<VolumeObservationOptions> volume_observation;
    IterativeSolveOptions state_solve;
    IterativeSolveOptions adjoint_solve;
    IterativeSolveOptions control_metric_solve{1000, 1.0e-12, 1.0e-14};
  };

  enum class ReducedMethod
  {
    steepest_descent,
    bfgs,
    limited_memory_bfgs
  };

  enum class ReducedGlobalization
  {
    armijo,
    fixed_step
  };

  inline const char *
  reduced_globalization_name(const ReducedGlobalization globalization)
  {
    switch (globalization)
      {
        case ReducedGlobalization::armijo:
          return "armijo";
        case ReducedGlobalization::fixed_step:
          return "fixed-step";
      }
    throw std::invalid_argument("unknown reduced globalization selection");
  }

  enum class ObjectiveTargetPolicy
  {
    none,
    explicit_value,
    matched_reference_method
  };

  struct SolverOptions
  {
    ReducedMethod                  method = ReducedMethod::steepest_descent;
    ReducedGlobalization           globalization = ReducedGlobalization::armijo;
    solvers::ReducedSolverParameters parameters;
    solvers::LimitedMemoryBfgsParameters limited_memory_bfgs;
    ObjectiveTargetPolicy              objective_target_policy =
      ObjectiveTargetPolicy::none;
    std::string objective_target_reference_method;
    // Retained for old benchmark files that recorded a source value without
    // executing it. New experiments use parameters.minimum_step_length.
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

  inline ScalarFunctionDefinition
  b1_manufactured_zero_forcing()
  {
    return {"manufactured-zero",
            ScalarFunctionKind::zero,
            0.0,
            "",
            "chapter-6.e6.5.1.manufactured-zero-forcing"};
  }

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
    std::vector<double> regularisation_sweep = {
      1.0e-1, 1.0e-2, 1.0e-3, 1.0e-6};
    ScalarFunctionDefinition forcing = b1_manufactured_zero_forcing();
  };

  enum class B2ObservationRegion
  {
    wings,
    full
  };

  inline constexpr std::array<B2ObservationRegion, 2>
    b2_observation_region_order = {{B2ObservationRegion::wings,
                                    B2ObservationRegion::full}};

  inline const char *
  b2_observation_region_name(const B2ObservationRegion region)
  {
    switch (region)
      {
        case B2ObservationRegion::wings:
          return "wings";
        case B2ObservationRegion::full:
          return "full";
      }
    throw std::invalid_argument("B2 has an unknown observation region");
  }

  inline ScalarFunctionDefinition
  b2_manufactured_constant_target()
  {
    return {"constant-2",
            ScalarFunctionKind::constant,
            2.0,
            "",
            "chapter-6.e6.5.2.target"};
  }

  inline ScalarFunctionDefinition
  b2_manufactured_parabolic_target()
  {
    return {"parabolic-4*x1*(1-x1)",
            ScalarFunctionKind::expression,
            0.0,
            "4.0*x1*(1.0-x1)",
            "chapter-6.e6.5.2.target"};
  }

  inline ScalarFunctionCatalog
  b2_manufactured_target_catalog(
    const std::string_view selected_id = "constant-2")
  {
    ScalarFunctionCatalog catalog{{b2_manufactured_constant_target(),
                                   b2_manufactured_parabolic_target()},
                                  std::string(selected_id)};
    validate_scalar_function_catalog(catalog, "B2 target catalog");
    return catalog;
  }

  inline const ScalarFunctionDefinition &
  b2_target_definition(const ScalarFunctionCatalog &catalog)
  {
    return selected_scalar_function_definition(catalog, "B2 target catalog");
  }

  inline const char *
  b2_manufactured_target_id(const std::string_view target_profile)
  {
    if (target_profile == "constant")
      return "constant-2";
    if (target_profile == "parabolic")
      return "parabolic-4*x1*(1-x1)";
    throw std::invalid_argument("B2 has an unknown target profile '" +
                                std::string(target_profile) + "'");
  }

  inline ScalarFunctionDefinition
  b2_manufactured_zero_forcing()
  {
    return {"zero",
            ScalarFunctionKind::zero,
            0.0,
            "",
            "chapter-6.e6.5.2.zero-forcing"};
  }

  inline std::string
  b2_case_name(const B2ObservationRegion region,
               const std::string_view  target_profile)
  {
    return std::string(b2_observation_region_name(region)) + "-" +
           std::string(target_profile);
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
    B2ObservationRegion observation_region = B2ObservationRegion::wings;
    std::string          target_profile = "constant";
    ScalarFunctionCatalog target_catalog = b2_manufactured_target_catalog();
    double                   fixed_temperature = 1.0;
    ScalarFunctionDefinition forcing = b2_manufactured_zero_forcing();
    semantic::v1::TransportBoundaryForm transport_boundary_form =
      semantic::v1::TransportBoundaryForm::ordinary_normal_minus_transport;
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
    const auto validate_solve = [](const IterativeSolveOptions &solve,
                                   const bool allow_automatic_iterations,
                                   const char *description) {
      if (!allow_automatic_iterations && solve.maximum_iterations == 0)
        throw std::invalid_argument(
          std::string(description) + " needs a positive iteration limit");
      if (!std::isfinite(solve.relative_tolerance) ||
          solve.relative_tolerance <= 0.0)
        throw std::invalid_argument(
          std::string(description) + " needs a positive finite relative tolerance");
      if (!std::isfinite(solve.absolute_tolerance) ||
          solve.absolute_tolerance <= 0.0)
        throw std::invalid_argument(
          std::string(description) + " needs a positive finite absolute tolerance");
    };

    if (options.mesh.dimension == 0)
      throw std::invalid_argument("benchmark scenarios need a mesh dimension");
    const bool has_axis_subdivisions =
      !options.mesh.axis_subdivisions.empty();
    if (has_axis_subdivisions &&
        options.mesh.axis_subdivisions.size() != options.mesh.dimension)
      throw std::invalid_argument(
        "axis mesh subdivisions must match the mesh dimension");
    if (std::find(options.mesh.axis_subdivisions.begin(),
                  options.mesh.axis_subdivisions.end(),
                  0U) != options.mesh.axis_subdivisions.end())
      throw std::invalid_argument("axis mesh subdivisions must be positive");
    if (options.mesh.subdivisions != 0 && has_axis_subdivisions)
      throw std::invalid_argument(
        "mesh subdivisions must select either isotropic or per-axis counts");
    switch (options.mesh.generation)
      {
        case MeshGeneration::framework_native:
          if (options.mesh.subdivisions != 0 ||
              has_axis_subdivisions ||
              options.mesh.centroid_splits != 0 ||
              options.mesh.selection_seed != 0)
            throw std::invalid_argument(
              "framework-native meshes do not use simplex generation parameters");
          break;
        case MeshGeneration::structured_simplex:
          if (options.mesh.dimension != 2 ||
              (options.mesh.subdivisions == 0 && !has_axis_subdivisions) ||
              options.mesh.refinement != 0 ||
              options.mesh.centroid_splits != 0 ||
              options.mesh.selection_seed != 0)
            throw std::invalid_argument(
              "structured-simplex meshes need two-dimensional positive isotropic or per-axis subdivisions and no refinement or split parameters");
          break;
        case MeshGeneration::centroid_split_simplex:
          if (options.mesh.dimension != 2 ||
              (options.mesh.subdivisions == 0 && !has_axis_subdivisions) ||
              options.mesh.refinement != 0 ||
              options.mesh.centroid_splits == 0)
            throw std::invalid_argument(
              "centroid-split-simplex meshes need two dimensions, positive isotropic or per-axis subdivisions and splits, and zero refinement");
          {
            const auto base_cell_pair_count = has_axis_subdivisions ?
              static_cast<std::uint64_t>(options.mesh.axis_subdivisions[0]) *
                options.mesh.axis_subdivisions[1] :
              static_cast<std::uint64_t>(options.mesh.subdivisions) *
                options.mesh.subdivisions;
            const auto requested_split_count =
              static_cast<std::uint64_t>(options.mesh.centroid_splits);
            if (requested_split_count > base_cell_pair_count &&
                requested_split_count - base_cell_pair_count >
                  base_cell_pair_count)
              throw std::invalid_argument(
                "centroid-split-simplex requests more splits than base triangles");
          }
          break;
        default:
          throw std::invalid_argument("unknown mesh generation selection");
      }
    if (options.state_degree == 0)
      throw std::invalid_argument("benchmark scenarios need a positive state degree");
    if (options.volume_observation)
      {
        if (options.volume_observation->quadrature_order == 0)
          throw std::invalid_argument(
            "volume observation needs a positive quadrature order");
        switch (options.volume_observation->target_realisation)
          {
            case VolumeObservationTargetRealisation::analytic_quadrature:
            case VolumeObservationTargetRealisation::state_fe_interpolation:
              break;
            default:
              throw std::invalid_argument(
                "unknown volume-observation target realisation");
          }
      }
    if (options.mesh.mesh_provenance.empty())
      throw std::invalid_argument(
        "benchmark scenarios need mesh provenance");
    if (!options.owned_session)
      throw std::invalid_argument(
        "selected benchmark scenarios require an owned compilation session");
    validate_solve(options.state_solve, true, "state solve");
    validate_solve(options.adjoint_solve, true, "adjoint solve");
    validate_solve(options.control_metric_solve,
                   false,
                   "control metric solve");
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
  validate_reduced_globalization(const ReducedGlobalization globalization)
  {
    switch (globalization)
      {
        case ReducedGlobalization::armijo:
        case ReducedGlobalization::fixed_step:
          return;
      }
    throw std::invalid_argument("unknown reduced globalization selection");
  }

  inline void
  validate_b1(const B1Scenario &scenario)
  {
    scenario.validate();
    validate_common_compile_options(scenario.compile);
    if (scenario.compile.volume_observation)
      throw std::invalid_argument(
        "B1 does not select a material-subdomain volume observation policy");
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
    switch (scenario.problem.recipe.discretisation)
      {
        case chapter5::DistributedControlDiscretisation::cellwise_constant:
          break;
        case chapter5::DistributedControlDiscretisation::
          homogeneous_dirichlet_continuous:
          if (scenario.problem.recipe.with_cellwise_box)
            throw std::invalid_argument(
              "B1 continuous control cannot use the cellwise box");
          break;
        default:
          throw std::invalid_argument(
            "B1 has an unknown control discretisation");
      }
    if (scenario.compile.mesh.generation != MeshGeneration::framework_native &&
        scenario.problem.recipe.discretisation != chapter5::
          DistributedControlDiscretisation::homogeneous_dirichlet_continuous)
      throw std::invalid_argument(
        "B1 simplex meshes require the continuous homogeneous-Dirichlet control target");
    if (!scenario.compile.mesh.axis_subdivisions.empty())
      throw std::invalid_argument(
        "B1 mesh construction supports only isotropic subdivisions");
    validate_scalar_function_definition(scenario.problem.forcing,
                                        "B1 forcing");
    if (scenario.problem.data.forcing_provenance !=
        scenario.problem.forcing.provenance)
      throw std::invalid_argument(
        "B1 forcing definition and runtime provenance must agree");
    if (scenario.solver.method != ReducedMethod::steepest_descent &&
        scenario.solver.method != ReducedMethod::limited_memory_bfgs)
      throw std::invalid_argument(
        "B1 selects steepest descent or limited-memory BFGS");
    validate_reduced_globalization(scenario.solver.globalization);
    if (scenario.solver.globalization != ReducedGlobalization::armijo)
      throw std::invalid_argument("B1 currently selects Armijo globalization");
    switch (scenario.solver.objective_target_policy)
      {
        case ObjectiveTargetPolicy::none:
          if (scenario.solver.parameters.objective_target.has_value())
            throw std::invalid_argument(
              "B1 objective target needs an explicit target policy");
          break;
        case ObjectiveTargetPolicy::explicit_value:
          if (!scenario.solver.parameters.objective_target.has_value())
            throw std::invalid_argument(
              "B1 explicit objective-target policy needs a target value");
          break;
        case ObjectiveTargetPolicy::matched_reference_method:
          if (scenario.solver.objective_target_reference_method.empty())
            throw std::invalid_argument(
              "B1 matched objective target needs a reference method");
          break;
      }
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
  validate_b2_transport_boundary_form(
    const semantic::v1::TransportBoundaryForm form)
  {
    switch (form)
      {
        case semantic::v1::TransportBoundaryForm::total_conormal:
        case semantic::v1::TransportBoundaryForm::
          ordinary_normal_minus_transport:
          return;
        default:
          throw std::invalid_argument(
            "B2 needs a total-conormal or ordinary-normal transport boundary form");
      }
  }

  inline void
  validate_b2_control_discretisation(
    const chapter5::NeumannConvectionParameters &parameters)
  {
    switch (parameters.control_discretisation)
      {
        case semantic::v1::NeumannControlDiscretisation::facewise_constant:
          return;
        case semantic::v1::NeumannControlDiscretisation::continuous_nodal_trace:
          if (parameters.with_facewise_box)
            throw std::invalid_argument(
              "B2 continuous Neumann control does not support a facewise box");
          return;
        case semantic::v1::NeumannControlDiscretisation::unspecified:
          throw std::invalid_argument(
            "B2 needs a facewise-constant or continuous nodal-trace control");
      }
    throw std::invalid_argument("B2 has an unknown control discretisation");
  }

  inline void
  validate_b2(const B2Scenario &scenario)
  {
    scenario.validate();
    validate_common_compile_options(scenario.compile);
    validate_runtime_data(scenario.problem.data);
    validate_b2_transport_boundary_form(
      scenario.problem.transport_boundary_form);
    validate_b2_control_discretisation(scenario.problem.recipe);
    if (!scenario.compile.volume_observation)
      throw std::invalid_argument(
        "B2 needs a volume-observation discretisation policy");
    switch (scenario.compile.mesh.generation)
      {
        case MeshGeneration::framework_native:
          break;
        case MeshGeneration::structured_simplex:
          if (scenario.compile.mesh.axis_subdivisions.empty())
            throw std::invalid_argument(
              "B2 structured-simplex meshes require per-axis subdivisions");
          break;
        case MeshGeneration::centroid_split_simplex:
          if (scenario.compile.mesh.axis_subdivisions.empty())
            throw std::invalid_argument(
              "B2 centroid-split-simplex meshes require per-axis subdivisions");
          break;
      }
    if (!std::isfinite(scenario.problem.fixed_temperature))
      throw std::invalid_argument("B2 fixed temperature must be finite");
    validate_scalar_function_catalog(scenario.problem.target_catalog,
                                     "B2 target catalog");
    if (scenario.problem.data.desired_state_provenance !=
        b2_target_definition(scenario.problem.target_catalog).provenance)
      throw std::invalid_argument(
        "B2 desired-state provenance does not match the selected target");
    validate_scalar_function_definition(scenario.problem.forcing, "B2 forcing");
    if (scenario.problem.data.forcing_provenance !=
        scenario.problem.forcing.provenance)
      throw std::invalid_argument(
        "B2 forcing definition and runtime provenance must agree");
    if (scenario.problem.data.fixed_dirichlet_data_provenance.empty() ||
        scenario.problem.data.conservative_transport_provenance.empty())
      throw std::invalid_argument(
        "B2 needs fixed-data and transport provenance");
    if (scenario.solver.method != ReducedMethod::bfgs)
      throw std::invalid_argument("B2 selects the BFGS method");
    validate_reduced_globalization(scenario.solver.globalization);
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
    validate_b2_transport_boundary_form(parameters.transport_boundary_form);
    validate_b2_control_discretisation(parameters.recipe);
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
    const bool ordinary_normal =
      parameters.transport_boundary_form ==
      semantic::v1::TransportBoundaryForm::ordinary_normal_minus_transport;
    partition_policy->selected_policy = ordinary_normal
      ? "disjoint complete fixed-Dirichlet, Neumann-control, and natural outflow boundary regions; the source Neumann datum is ordinary normal derivative minus transport, with no diffusion scaling"
      : "disjoint complete fixed-Dirichlet, Neumann-control, and natural outflow boundary regions; the diagnostic datum is the outward diffusion-minus-transport total conormal";
    partition_policy->typed_selection = semantic::v1::BoundaryRealisationSelection{
      partition_policy->id,
      "state",
      b2_fixed_boundary_region_id,
      b2_outflow_boundary_region_id,
      {},
      {},
      b2_outflow_boundary_region_id,
      ordinary_normal ?
        semantic::v1::ConormalForm::unspecified :
        semantic::v1::ConormalForm::diffusion_minus_transport,
      semantic::v1::NormalOrientation::outward,
      semantic::v1::TraceEvaluationRealisation::fe_q_state_trace,
      semantic::v1::FaceQuadratureRealisation::qgauss_face,
      parameters.transport_boundary_form};
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
    const ReducedMethod method = ReducedMethod::steepest_descent,
    ScalarFunctionDefinition forcing = b1_manufactured_zero_forcing())
  {
    validate_scalar_function_definition(forcing, "B1 forcing");
    B1ProblemParameters problem;
    problem.forcing = std::move(forcing);
    problem.data.forcing_provenance = problem.forcing.provenance;

    SolverOptions solver;
    solver.method = method;
    solver.parameters.maximum_line_search_trials = 5;
    solver.parameters.gradient_tolerance =
      method == ReducedMethod::steepest_descent ? 1.0e-3 : 1.0e-8;
    solver.parameters.armijo_fraction = 1.0e-5;
    solver.parameters.backtracking_factor = 0.7;
    solver.parameters.minimum_step_length =
      method == ReducedMethod::steepest_descent ? 0.2 : 0.01;

    CompileOptions compile;
    compile.mesh.refinement = b1_default_mesh_refinement;

    B1Scenario scenario{
      {"chapter-6.b1.distributed-laplace",
       "E6.5.1 distributed Laplace control",
       "Reduced-space comparison of steepest descent and L-BFGS",
       "chapter-6",
       b1_recipe_id,
       {"B0 harness",
        "declarative forcing definition and provenance are recorded"}},
      std::move(problem),
      std::move(compile),
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
  make_b2_scenario_with_target_catalog(
    const B2ObservationRegion observation_region,
    const std::string_view target_profile,
    ScalarFunctionCatalog target_catalog,
    const semantic::v1::TransportBoundaryForm transport_boundary_form =
      semantic::v1::TransportBoundaryForm::ordinary_normal_minus_transport,
    const semantic::v1::NeumannControlDiscretisation control_discretisation =
      semantic::v1::NeumannControlDiscretisation::facewise_constant,
    const ReducedGlobalization globalization = ReducedGlobalization::armijo,
    const VolumeObservationOptions volume_observation = {})
  {
    B2ProblemParameters problem;
    problem.observation_region = observation_region;
    problem.target_profile = target_profile;
    problem.target_catalog = std::move(target_catalog);
    problem.data.desired_state_provenance =
      b2_target_definition(problem.target_catalog).provenance;
    problem.transport_boundary_form = transport_boundary_form;
    problem.recipe.control_discretisation =
      control_discretisation;

    CompileOptions compile;
    compile.mesh.refinement = b2_default_mesh_refinement;
    compile.volume_observation = volume_observation;

    std::string scenario_id = "chapter-6.b2.graetz-flow";
    if (observation_region != B2ObservationRegion::wings ||
        target_profile != "constant")
      scenario_id += "." + b2_case_name(observation_region, target_profile);

    B2Scenario scenario{
      {scenario_id,
       "E6.5.2 Graetz-flow boundary control",
       "Neumann control with downstream subdomain tracking",
       "chapter-6",
       b2_recipe_id,
       {"B0 harness", "fixed temperature port", "declared Galerkin policy"}},
      std::move(problem),
      std::move(compile),
      {ReducedMethod::bfgs,
       globalization,
       {},
       {},
       ObjectiveTargetPolicy::none,
       "",
       std::nullopt,
       "zero"},
      {scenario_id,
       "E6.5.2 equation (6.65), Table 6.2, Figures 6.4-6.5",
       chapter6_numerical_examples_source_revision,
       "release-dealii",
       true,
       {true, true, true, false, "runs"}}};
    validate_b2(scenario);
    return scenario;
  }

  inline B2Scenario
  make_b2_scenario(
    const B2ObservationRegion observation_region = B2ObservationRegion::wings,
    const std::string_view target_profile = "constant",
    const semantic::v1::TransportBoundaryForm transport_boundary_form =
      semantic::v1::TransportBoundaryForm::ordinary_normal_minus_transport,
    const semantic::v1::NeumannControlDiscretisation control_discretisation =
      semantic::v1::NeumannControlDiscretisation::facewise_constant,
    const ReducedGlobalization globalization = ReducedGlobalization::armijo,
    const VolumeObservationOptions volume_observation = {})
  {
    return make_b2_scenario_with_target_catalog(
      observation_region,
      target_profile,
      b2_manufactured_target_catalog(b2_manufactured_target_id(target_profile)),
      transport_boundary_form,
      control_discretisation,
      globalization,
      volume_observation);
  }

  inline ApplicationCatalog
  make_catalog()
  {
    ApplicationCatalog catalog;
    const auto b1 = make_b1_scenario();
    catalog.add(b1.metadata);
    for (const auto observation_region : b2_observation_region_order)
      for (const auto target_profile : {"constant", "parabolic"})
        catalog.add(make_b2_scenario(observation_region, target_profile).metadata);
    return catalog;
  }
} // namespace nmopt::application::chapter6
