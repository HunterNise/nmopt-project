#include "nmopt/application/application.hpp"
#include "../support/scenario_dispatch.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
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
    require(std::abs(b1.solver.parameters.armijo_fraction - 1.0e-5) < 1.0e-15,
            "B1 did not retain its source Armijo fraction");
    require(b1.solver.declared_minimum_step_length == 0.01,
            "B1 did not retain the L-BFGS minimum-step declaration");
    require(b1.experiment.harness.deterministic,
            "B1 did not retain the deterministic B0 harness policy");

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
  }

  void
  test_chapter6_catalog_discovers_standard_scenarios()
  {
    const auto catalog = nmopt::application::chapter6::make_catalog();
    require(catalog.entries().size() == 2,
            "Chapter 6 catalog did not retain B1 and B2");
    require(catalog.find("chapter-6.b1.distributed-laplace") != nullptr,
            "Chapter 6 catalog cannot discover B1");
    require(catalog.find("chapter-6.b2.graetz-flow") != nullptr,
            "Chapter 6 catalog cannot discover B2");
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
        {"chapter6_catalog_discovers_standard_scenarios",
         "nmopt.application.chapter6_catalog_discovers_standard_scenarios",
         {"backend-neutral", "application", "benchmark", "contract"},
         30,
         test_chapter6_catalog_discovers_standard_scenarios}};
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
