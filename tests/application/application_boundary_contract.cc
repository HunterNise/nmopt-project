#include "nmopt/application/application.hpp"
#include "../support/scenario_dispatch.hpp"

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
         test_catalog_discovers_metadata_and_rejects_duplicates}};
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
