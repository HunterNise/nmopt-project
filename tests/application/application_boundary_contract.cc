#include "nmopt/application/application.hpp"
#include "../support/scenario_dispatch.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
  using namespace nmopt::application;

  struct ProblemParameters
  {
    double regularisation = 0.0;
  };

  struct CompileOptions
  {
    unsigned int state_degree = 1;
  };

  struct SolverOptions
  {
    unsigned int maximum_iterations = 10;
  };

  struct ExperimentOptions
  {
    std::string output_id;
  };

  using Scenario =
    ScenarioT<ProblemParameters,
              CompileOptions,
              SolverOptions,
              ExperimentOptions>;

  void
  require(const bool condition, const char *message)
  {
    if (!condition)
      throw std::runtime_error(message);
  }

  void
  test_recipe_builds_problem_spec()
  {
    const ProblemRecipeT<ProblemParameters> recipe(
      RecipeMetadata{"chapter-5.scalar",
                     "Scalar volume control",
                     "A typed scalar Chapter 5 recipe",
                     "chapter-5",
                     {"fixed Dirichlet", "L2 control metric"}},
      [](const ProblemParameters &parameters) {
        nmopt::semantic::v1::ProblemSpec specification;
        specification.id = "scalar_volume_control";
        specification.label = "Scalar volume control";
        specification.data.push_back(
          {"regularisation_weight",
           "Control regularisation",
           nmopt::semantic::v1::DataKind::scalar_constant,
           nmopt::semantic::v1::DataRole::regularisation_weight,
           ""});
        require(parameters.regularisation > 0.0,
                "recipe parameters were not passed to the builder");
        return specification;
      });

    const ProblemParameters parameters{0.1};
    const auto specification = recipe(parameters);
    require(specification.id == "scalar_volume_control",
            "recipe did not return its semantic ProblemSpec");
    require(recipe.metadata().id == "chapter-5.scalar",
            "recipe metadata has the wrong stable id");
  }

  void
  test_scenario_binds_typed_options()
  {
    const Scenario scenario{
      ScenarioMetadata{"chapter-6.b1.default",
                       "Distributed Laplace benchmark",
                       "A frozen reduced-space scenario",
                       "chapter-6",
                       "chapter-5.scalar",
                       {"owned compilation session"}},
      ProblemParameters{0.01},
      CompileOptions{1},
      SolverOptions{50},
      ExperimentOptions{"b1-default"}};

    scenario.validate();
    require(scenario.metadata.recipe_id == "chapter-5.scalar",
            "scenario is not linked to its recipe");
    require(scenario.solver.maximum_iterations == 50,
            "scenario solver options were not retained");
  }

  void
  test_catalog_discovers_metadata_and_rejects_duplicates()
  {
    const RecipeMetadata recipe{
      "chapter-5.scalar", "Scalar volume control", "", "chapter-5", {}};
    const ScenarioMetadata scenario{"chapter-6.b1.default",
                                   "Distributed Laplace benchmark",
                                   "",
                                   "chapter-6",
                                   recipe.id,
                                   {}};

    ApplicationCatalog catalog;
    catalog.add(recipe);
    catalog.add(scenario);

    require(catalog.entries().size() == 2,
            "catalog did not retain both metadata entries");
    const auto *scenario_entry = catalog.find(scenario.id);
    require(scenario_entry != nullptr,
            "catalog could not find a registered scenario");
    require(scenario_entry->kind == CatalogEntryKind::scenario,
            "catalog assigned the wrong entry kind");
    require(scenario_entry->recipe_id == recipe.id,
            "catalog lost the scenario recipe reference");

    bool duplicate_rejected = false;
    try
      {
        catalog.add(recipe);
      }
    catch (const std::invalid_argument &)
      {
        duplicate_rejected = true;
      }
    require(duplicate_rejected, "catalog accepted a duplicate stable id");
  }

  void
  test_chapter5_recipe_records_build_registered_graphs()
  {
    const auto distributed =
      nmopt::application::chapter5::make_scalar_distributed_recipe();
    const auto distributed_spec =
      distributed(nmopt::application::chapter5::ScalarDistributedControlParameters{
        true});
    require(distributed.metadata().id ==
              nmopt::application::chapter5::scalar_distributed_recipe_id,
            "Chapter 5 distributed recipe has the wrong stable ID");
    require(distributed_spec.id ==
              "scalar_diffusion_reaction_volume_control",
            "distributed recipe did not select the scalar reference graph");
    require(!distributed_spec.constraints.empty(),
            "distributed recipe did not retain its cellwise box choice");

    const auto general = nmopt::application::chapter5::make_general_scalar_recipe();
    const auto general_spec = general(
      nmopt::application::chapter5::GeneralScalarParameters{{0}, {1}, false});
    require(general_spec.id == "general_scalar_elliptic_robin_volume_control",
            "general recipe did not select the Robin graph");
    require(std::any_of(general_spec.data.begin(),
                        general_spec.data.end(),
                        [](const auto &data) {
                          return data.id == "robin_coefficient";
                        }),
            "general recipe lost the Robin data port");

    const auto neumann =
      nmopt::application::chapter5::make_neumann_convection_recipe();
    const auto neumann_spec = neumann(
      nmopt::application::chapter5::NeumannConvectionParameters{1, false});
    require(neumann_spec.id == "scalar_convection_neumann_subdomain_control",
            "Neumann recipe did not select the C5.6 graph");
    require(std::any_of(neumann_spec.data.begin(),
                        neumann_spec.data.end(),
                        [](const auto &data) {
                          return data.role ==
                                 nmopt::semantic::v1::DataRole::conservative_transport;
                        }),
            "Neumann recipe lost its transport data port");
  }

  void
  test_chapter6_scenario_records_freeze_b1_and_b2_choices()
  {
    const auto b1 = nmopt::application::chapter6::make_b1_scenario(
      nmopt::application::chapter6::ReducedMethod::limited_memory_bfgs);
    require(b1.metadata.recipe_id ==
              nmopt::application::chapter6::b1_recipe_id,
            "B1 scenario is not linked to the distributed recipe");
    require(b1.solver.parameters.maximum_line_search_trials == 5,
            "B1 did not retain its source line-search trial limit");
    require(std::abs(b1.solver.parameters.gradient_tolerance - 1.0e-8) <
              1.0e-15,
            "B1 did not retain its declared L-BFGS stopping tolerance");
    require(std::abs(b1.solver.parameters.armijo_fraction - 1.0e-5) < 1.0e-15,
            "B1 did not retain its source Armijo fraction");
    require(b1.solver.declared_minimum_step_length == 0.01,
            "B1 did not retain the L-BFGS minimum-step declaration");
    require(b1.experiment.harness.deterministic,
            "B1 did not retain the deterministic B0 harness policy");
    require(b1.experiment.source_revision ==
              nmopt::application::chapter6::
                chapter6_numerical_examples_source_revision,
            "B1 did not retain the frozen source catalogue revision");
    require(b1.problem.forcing_selection ==
              nmopt::application::chapter6::B1ForcingSelection::manufactured_zero,
            "B1 did not make its manufactured forcing choice explicit");
    require(b1.problem.data.forcing_provenance ==
              "chapter-6.e6.5.1.manufactured-zero-forcing",
            "B1 did not retain manufactured forcing provenance");

    const auto recovered = nmopt::application::chapter6::make_b1_scenario(
      nmopt::application::chapter6::ReducedMethod::steepest_descent,
      nmopt::application::chapter6::B1ForcingSelection::recovered_source);
    require(recovered.problem.forcing_selection ==
              nmopt::application::chapter6::B1ForcingSelection::recovered_source,
            "B1 did not retain the recovered forcing choice");
    require(recovered.problem.data.forcing_provenance ==
              "chapter-6.e6.5.1.recovered-forcing",
            "B1 did not retain recovered forcing provenance");

    const auto steepest = nmopt::application::chapter6::make_b1_scenario(
      nmopt::application::chapter6::ReducedMethod::steepest_descent);
    require(std::abs(steepest.solver.parameters.gradient_tolerance - 1.0e-3) <
              1.0e-15,
            "B1 did not retain the source steepest-descent tolerance");

    const auto b2 = nmopt::application::chapter6::make_b2_scenario(
      nmopt::application::chapter6::GraetzCase::observation_full_parabolic_target);
    require(b2.metadata.recipe_id ==
              nmopt::application::chapter6::b2_recipe_id,
            "B2 scenario is not linked to the Neumann recipe");
    require(b2.problem.graetz_case ==
              nmopt::application::chapter6::GraetzCase::observation_full_parabolic_target,
            "B2 did not retain its observation/target case");
    require(b2.problem.fixed_temperature == 1.0,
            "B2 did not retain the fixed temperature");
    require(b2.problem.data.conservative_transport_provenance ==
              "chapter-6.e6.5.2.graetz-transport",
            "B2 did not retain transport provenance");
    require(b2.experiment.source_revision ==
              nmopt::application::chapter6::
                chapter6_numerical_examples_source_revision,
            "B2 did not retain the frozen source catalogue revision");
    require(b2.metadata.id == "chapter-6.b2.graetz-flow.full-parabolic",
            "B2 did not assign a stable case-specific scenario ID");
  }

  void
  test_b1_scenario_assembles_distributed_problem_spec()
  {
    const auto scenario = nmopt::application::chapter6::make_b1_scenario(
      nmopt::application::chapter6::ReducedMethod::limited_memory_bfgs);
    const auto specification =
      nmopt::application::chapter6::make_b1_problem_spec(scenario);
    require(specification.id ==
              "scalar_diffusion_reaction_volume_control",
            "B1 spec assembly selected the wrong semantic graph");
    require(specification.formulation.id == "reduced_dto" &&
              specification.formulation.kind ==
                nmopt::semantic::v1::FormulationKind::reduced_dto,
            "B1 spec assembly did not preserve the reduced-DTO product");
    require(std::any_of(specification.data.begin(),
                        specification.data.end(),
                        [](const auto &data) {
                          return data.id == "forcing";
                        }),
            "B1 spec assembly lost the forcing port");
    require(std::any_of(specification.data.begin(),
                        specification.data.end(),
                        [](const auto &data) {
                          return data.id == "desired_state";
                        }),
            "B1 spec assembly lost the desired-state port");

    const auto parameter_spec =
      nmopt::application::chapter6::make_b1_problem_spec(scenario.problem);
    require(parameter_spec.id == specification.id &&
              parameter_spec.variables.size() == specification.variables.size(),
            "B1 parameter and scenario spec assembly disagree");
  }

  void
  test_b2_scenarios_assemble_fixed_temperature_problem_specs()
  {
    const auto validator = nmopt::semantic::v1::SemanticValidator{};
    for (const auto graetz_case :
         nmopt::application::chapter6::b2_case_order)
      {
        const auto scenario =
          nmopt::application::chapter6::make_b2_scenario(graetz_case);
        const auto specification =
          nmopt::application::chapter6::make_b2_problem_spec(scenario);
        require(specification.id ==
                  "scalar_convection_neumann_subdomain_control",
                "B2 spec assembly selected the wrong semantic graph");
        const auto fixed_data = std::find_if(
          specification.data.begin(),
          specification.data.end(),
          [](const auto &data) { return data.id == "fixed_dirichlet_data"; });
        require(fixed_data != specification.data.end() &&
                  fixed_data->role ==
                    nmopt::semantic::v1::DataRole::fixed_dirichlet_lifting,
                "B2 spec assembly lost the fixed-data lifting port");
        const auto reconstruction = std::find_if(
          specification.transformations.begin(),
          specification.transformations.end(),
          [](const auto &transformation) {
            return transformation.id == "fixed_dirichlet_reconstruction";
          });
        require(reconstruction != specification.transformations.end() &&
                  reconstruction->fixed_data_id == "fixed_dirichlet_data" &&
                  reconstruction->kind ==
                    nmopt::semantic::v1::TransformationKind::fixed_dirichlet_reconstruction,
                "B2 spec assembly lost the fixed-data reconstruction");
        const auto state = std::find_if(
          specification.variables.begin(),
          specification.variables.end(),
          [](const auto &variable) { return variable.id == "state"; });
        require(state != specification.variables.end() &&
                  state->physical_field_transform_id ==
                    "fixed_dirichlet_reconstruction",
                "B2 state variable is not connected to fixed reconstruction");
        require(validator.validate(specification).valid(),
                "B2 fixed-temperature semantic graph is invalid");
      }
  }

  void
  test_chapter6_catalog_discovers_standard_scenarios()
  {
    const auto catalog = nmopt::application::chapter6::make_catalog();
    require(catalog.entries().size() == 5,
            "Chapter 6 catalog did not retain B1 and all B2 cases");
    require(catalog.find("chapter-6.b1.distributed-laplace") != nullptr,
            "Chapter 6 catalog cannot discover B1");
    for (const auto graetz_case :
         nmopt::application::chapter6::b2_case_order)
      {
        const auto id = std::string("chapter-6.b2.graetz-flow") +
                        (graetz_case ==
                           nmopt::application::chapter6::GraetzCase::observation_wings_constant_target
                           ? ""
                           : "." + std::string(
                                     nmopt::application::chapter6::graetz_case_name(
                                       graetz_case)));
        require(catalog.find(id) != nullptr,
                "Chapter 6 catalog cannot discover a B2 case");
      }
  }

  void
  test_b0_harness_captures_deterministic_artifact_boundary()
  {
    const auto scenario = nmopt::application::chapter6::make_b1_scenario();
    using Harness =
      nmopt::application::benchmark::BenchmarkHarnessT<decltype(scenario)>;
    struct DetachedEnvelope
    {
      int marker = 0;
    };

    Harness harness(scenario);
    nmopt::application::benchmark::BenchmarkMeasurements measurements;
    measurements.timing_collected = true;
    measurements.wall_seconds = 1.25;
    measurements.cpu_seconds = 0.75;
    measurements.memory_collected = true;
    measurements.peak_memory_bytes = 4096;

    const auto artifact = harness.finalize(
      DetachedEnvelope{7},
      nmopt::semantic::v1::ValidationReport{},
      measurements,
      {"state", "control", "objective"});

    require(artifact.identity().scenario_id ==
              "chapter-6.b1.distributed-laplace",
            "B0 artifact lost the scenario identity");
    require(artifact.identity().recipe_id ==
              nmopt::application::chapter6::b1_recipe_id,
            "B0 artifact lost the recipe identity");
    require(artifact.identity().deterministic,
            "B0 artifact lost deterministic harness metadata");
    require(artifact.measurements().wall_seconds == 1.25 &&
              artifact.measurements().peak_memory_bytes == 4096,
            "B0 artifact lost measurement metadata");
    require(artifact.envelope().marker == 7,
            "B0 artifact lost the detached experiment envelope");
    require(artifact.selected_fields().size() == 3,
            "B0 artifact lost selected output fields");

    nmopt::application::benchmark::BenchmarkMeasurements invalid_measurements;
    invalid_measurements.timing_collected = true;
    invalid_measurements.wall_seconds = -1.0;
    bool rejected = false;
    try
      {
        (void)harness.finalize(DetachedEnvelope{},
                               nmopt::semantic::v1::ValidationReport{},
                               invalid_measurements);
      }
    catch (const std::invalid_argument &)
      {
        rejected = true;
      }
    require(rejected, "B0 harness accepted invalid timing evidence");
  }

  void
  test_b0_artifact_writer_is_deterministic_and_escapes_values()
  {
    const auto scenario = nmopt::application::chapter6::make_b1_scenario();
    const nmopt::application::benchmark::BenchmarkHarnessT<decltype(scenario)>
      harness(scenario);
    struct DetachedEnvelope
    {};
    const auto artifact = harness.finalize(DetachedEnvelope{});

    const nmopt::application::benchmark::BenchmarkArtifactWriter writer;
    const std::vector<nmopt::application::benchmark::ArtifactField> fields{
      {"result.objective", "1=2\nnext"},
      {"result.status", "ok"}};
    const auto first = writer.render(artifact, fields);
    const auto second = writer.render(artifact, fields);
    require(first == second,
            "B0 artifact writer did not produce deterministic output");
    require(first.find("result.objective=1\\=2\\nnext\n") !=
              std::string::npos,
            "B0 artifact writer did not escape field values");
    require(first.find("identity.scenario_id=chapter-6.b1.distributed-laplace\n") !=
              std::string::npos,
            "B0 artifact writer omitted the scenario identity");

    std::ostringstream stream;
    writer.write(stream, artifact, fields);
    require(stream.str() == first,
            "B0 artifact stream writer disagrees with render");

    bool duplicate_rejected = false;
    try
      {
        (void)writer.render(
          artifact,
          {{"duplicate", "a"}, {"duplicate", "b"}});
      }
    catch (const std::invalid_argument &)
      {
        duplicate_rejected = true;
      }
    require(duplicate_rejected,
            "B0 artifact writer accepted duplicate field keys");
  }

  void
  test_b0_runner_delegates_build_and_execution_callbacks()
  {
    const auto scenario = nmopt::application::chapter6::make_b1_scenario();
    using Runner =
      nmopt::application::benchmark::HeadlessBenchmarkRunnerT<decltype(scenario)>;
    struct DetachedEnvelope
    {
      int marker = 0;
    };

    bool built = false;
    bool executed = false;
    Runner runner(scenario);
    const auto result = runner.run(
      [&built](const auto &parameters) {
        built = true;
        require(!parameters.regularisation_sweep.empty(),
                "B0 runner did not pass typed problem parameters");
        return nmopt::application::chapter6::make_b1_problem_spec(parameters);
      },
      [&executed](const auto &specification, const auto &received_scenario) {
        executed = true;
        require(specification.id ==
                  "scalar_diffusion_reaction_volume_control",
                "B0 runner passed the wrong ProblemSpec to execution");
        require(received_scenario.metadata.id ==
                  "chapter-6.b1.distributed-laplace",
                "B0 runner passed the wrong scenario to execution");
        using Evidence =
          nmopt::application::benchmark::BenchmarkExecutionEvidenceT<
            DetachedEnvelope>;
        return Evidence{DetachedEnvelope{11},
                         nmopt::semantic::v1::ValidationReport{},
                         {},
                         {"objective", "control"},
                         {{"compile.product", "reduced_dto"},
                          {"compile.execution", "assembled"}}};
      });

    require(built && executed,
            "B0 runner did not invoke both orchestration callbacks");
    require(result.artifact.envelope().marker == 11,
            "B0 runner lost the execution envelope");
    require(result.artifact.measurements().timing_collected &&
              result.artifact.measurements().wall_seconds >= 0.0,
            "B0 runner did not capture orchestration timing");
    require(result.document.find("compile.product=reduced_dto\n") !=
              std::string::npos,
            "B0 runner did not render execution evidence");
  }
} // namespace

int
main(const int argc, char **argv)
{
  try
    {
      const std::vector<nmopt::test_support::Scenario> scenarios{
        {"recipe_builds_problem_spec",
         "nmopt.application.recipe_builds_problem_spec",
         {"backend-neutral", "application", "contract"},
         30,
         test_recipe_builds_problem_spec},
        {"scenario_binds_typed_options",
         "nmopt.application.scenario_binds_typed_options",
         {"backend-neutral", "application", "contract"},
         30,
         test_scenario_binds_typed_options},
        {"catalog_discovers_metadata_and_rejects_duplicates",
         "nmopt.application.catalog_discovers_metadata_and_rejects_duplicates",
         {"backend-neutral", "application", "contract", "negative"},
         30,
         test_catalog_discovers_metadata_and_rejects_duplicates},
        {"chapter5_recipe_records_build_registered_graphs",
         "nmopt.application.chapter5_recipe_records_build_registered_graphs",
         {"backend-neutral", "application", "semantic", "contract"},
         30,
         test_chapter5_recipe_records_build_registered_graphs},
        {"chapter6_scenario_records_freeze_b1_and_b2_choices",
         "nmopt.application.chapter6_scenario_records_freeze_b1_and_b2_choices",
         {"backend-neutral", "application", "benchmark", "contract"},
         30,
         test_chapter6_scenario_records_freeze_b1_and_b2_choices},
        {"b1_scenario_assembles_distributed_problem_spec",
         "nmopt.application.b1_scenario_assembles_distributed_problem_spec",
         {"backend-neutral", "application", "benchmark", "semantic", "contract"},
         30,
         test_b1_scenario_assembles_distributed_problem_spec},
        {"b2_scenarios_assemble_fixed_temperature_problem_specs",
         "nmopt.application.b2_scenarios_assemble_fixed_temperature_problem_specs",
         {"backend-neutral", "application", "benchmark", "semantic", "contract"},
         30,
         test_b2_scenarios_assemble_fixed_temperature_problem_specs},
        {"chapter6_catalog_discovers_standard_scenarios",
         "nmopt.application.chapter6_catalog_discovers_standard_scenarios",
         {"backend-neutral", "application", "benchmark", "contract"},
         30,
         test_chapter6_catalog_discovers_standard_scenarios},
        {"b0_harness_captures_deterministic_artifact_boundary",
         "nmopt.application.b0_harness_captures_deterministic_artifact_boundary",
         {"backend-neutral", "application", "benchmark", "contract", "negative"},
         30,
         test_b0_harness_captures_deterministic_artifact_boundary},
        {"b0_artifact_writer_is_deterministic_and_escapes_values",
         "nmopt.application.b0_artifact_writer_is_deterministic_and_escapes_values",
         {"backend-neutral", "application", "benchmark", "contract", "negative"},
         30,
         test_b0_artifact_writer_is_deterministic_and_escapes_values},
        {"b0_runner_delegates_build_and_execution_callbacks",
         "nmopt.application.b0_runner_delegates_build_and_execution_callbacks",
         {"backend-neutral", "application", "benchmark", "contract"},
         30,
         test_b0_runner_delegates_build_and_execution_callbacks}};
      const auto result = nmopt::test_support::run_requested_scenarios(
        argc, argv, scenarios, std::cout);
      if (!result.listed)
        std::cout << "application boundary contract test passed: "
                  << result.executed << '\n';
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "application boundary contract test failed: "
                << exception.what() << '\n';
      return 1;
    }
}
