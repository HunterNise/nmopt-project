#include "nmopt/contract/complementarity.hpp"
#include "test_support/contract_errors.hpp"
#include "test_support/scenario_dispatch.hpp"

#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

namespace
{
  using namespace nmopt::contract;

  struct Fixture
  {
    LayoutPtr layout = std::make_shared<const BlockLayout>(
      "control", std::vector<SpaceId>{{"control"}},
      std::vector<std::size_t>{4});
    DiagonalMetric metric{
      "l2_cellwise",
      layout,
      {DenseVector{2.0, 4.0, 1.0, 3.0}}};
    BoxBounds bounds{
      layout,
      PrimalBlock(layout, {DenseVector{0.0, -1.0, 0.0, 0.0}}),
      PrimalBlock(layout, {DenseVector{1.0, 2.0, 3.0, 4.0}})};
    BoxComplementarity complementarity{
      bounds,
      make_metric_multiplier_representation(metric)};
  };

  void
  require_equal(const double actual,
                const double expected,
                const char *message)
  {
    require(std::abs(actual - expected) <= 1e-12, message);
  }

  void
  test_bounds_and_metric_representation()
  {
    const Fixture fixture;
    const PrimalBlock primal(
      fixture.layout, {DenseVector{1.0, -0.5, 1.5, 2.0}});
    const auto multiplier = fixture.complementarity.primal_to_multiplier(
      PrimalBlock(fixture.layout, {DenseVector{1.0, -1.0, 0.0, 1.0}}));
    const auto representative =
      fixture.complementarity.multiplier_to_primal(multiplier);

    require(fixture.bounds.is_feasible(primal),
            "valid cellwise control was rejected by box bounds");
    require_equal(representative.block(0)[0], 1.0,
                  "metric dual-to-primal conversion changed entry 0");
    require_equal(representative.block(0)[1], -1.0,
                  "metric dual-to-primal conversion changed entry 1");
    require_equal(representative.block(0)[2], 0.0,
                  "metric dual-to-primal conversion changed entry 2");
    require_equal(representative.block(0)[3], 1.0,
                  "metric dual-to-primal conversion changed entry 3");

    const PrimalBlock infeasible(
      fixture.layout, {DenseVector{-0.1, -0.5, 1.5, 4.1}});
    require(!fixture.bounds.is_feasible(infeasible),
            "infeasible cellwise control was accepted by box bounds");
  }

  void
  test_classification_and_selection()
  {
    const Fixture fixture;
    const PrimalBlock primal(
      fixture.layout, {DenseVector{1.0, -0.5, 1.5, 2.0}});
    const auto multiplier = fixture.complementarity.primal_to_multiplier(
      PrimalBlock(fixture.layout, {DenseVector{1.0, -1.0, 0.0, 1.0}}));
    const auto selection = fixture.complementarity.classify(primal,
                                                            multiplier,
                                                            1.0);

    require(selection.activities() ==
              std::vector<BoxActivity>{BoxActivity::upper,
                                       BoxActivity::lower,
                                       BoxActivity::inactive,
                                       BoxActivity::inactive},
            "cellwise box classification selected the wrong active sets");
    require(selection.active_size() == 2 && selection.free_size() == 2,
            "cellwise box classification reported wrong set sizes");
    require(selection.active_indices() == std::vector<std::size_t>{0, 1},
            "active-set selection returned the wrong active indices");
    require(selection.free_indices() == std::vector<std::size_t>{2, 3},
            "active-set selection returned the wrong free indices");

    const DenseVector free = selection.restrict_free(primal.block(0));
    const DenseVector active = selection.restrict_active(primal.block(0));
    require(free.values() == std::vector<double>{1.5, 2.0},
            "active-set free restriction returned the wrong values");
    require(active.values() == std::vector<double>{1.0, -0.5},
            "active-set active restriction returned the wrong values");

    const DenseVector prolonged = selection.prolong(
      DenseVector{1.25, 2.5}, DenseVector{0.75, -0.25});
    require(prolonged.values() ==
              std::vector<double>{0.75, -0.25, 1.25, 2.5},
            "active-set prolongation returned the wrong values");

    const ActiveSetSelection same_selection(
      fixture.layout,
      {BoxActivity::upper,
       BoxActivity::lower,
       BoxActivity::inactive,
       BoxActivity::inactive});
    require(selection == same_selection,
            "equal active-set selections did not compare equal");
  }

  void
  test_invalid_contract_inputs()
  {
    const auto layout = std::make_shared<const BlockLayout>(
      "control", std::vector<SpaceId>{{"control"}},
      std::vector<std::size_t>{2});
    const auto lower = PrimalBlock(layout, {DenseVector{0.0, 1.0}});
    const auto upper = PrimalBlock(layout, {DenseVector{1.0, 0.0}});
    nmopt::test_support::require_contract_error(
      [&] { (void)BoxBounds(layout, lower, upper); },
      "Box lower bound exceeds upper bound",
      "box bounds accepted reversed limits");

    const DiagonalMetric metric(
      "l2_cellwise", layout, {DenseVector{1.0, 1.0}});
    const BoxBounds bounds(
      layout,
      PrimalBlock(layout, {DenseVector{0.0, 0.0}}),
      PrimalBlock(layout, {DenseVector{1.0, 1.0}}));
    const BoxComplementarity complementarity(
      bounds, make_metric_multiplier_representation(metric));
    const PrimalBlock primal(layout, {DenseVector{0.5, 0.5}});
    const auto multiplier = complementarity.primal_to_multiplier(
      PrimalBlock(layout, {DenseVector{0.0, 0.0}}));
    nmopt::test_support::require_contract_error(
      [&] { (void)complementarity.classify(primal, multiplier, 0.0); },
      "Box classification parameter must be positive and finite",
      "box classification accepted a zero parameter");

    const auto wrong_layout = std::make_shared<const BlockLayout>(
      "wrong", std::vector<SpaceId>{{"wrong"}},
      std::vector<std::size_t>{2});
    const auto wrong_multiplier = CovectorBlock(
      wrong_layout, {DenseVector{0.0, 0.0}});
    nmopt::test_support::require_contract_error(
      [&] { (void)complementarity.classify(primal, wrong_multiplier, 1.0); },
      "Box multiplier dual-to-primal input has an incompatible layout",
      "box classification accepted an incompatible multiplier layout");

    const ActiveSetSelection selection(
      layout, {BoxActivity::inactive, BoxActivity::upper});
    nmopt::test_support::require_contract_error(
      [&] { (void)selection.restrict_free(DenseVector{0.0}); },
      "Active-set free restriction input has the wrong size",
      "active-set restriction accepted the wrong vector size");
  }
} // namespace

int
main(const int argc, char **argv)
{
  try
    {
      const std::vector<nmopt::test_support::Scenario> scenarios{
        {"bounds_and_metric_representation",
         "nmopt.complementarity.bounds_and_metric_representation",
         {"backend-neutral", "contract", "complementarity"},
         30,
         test_bounds_and_metric_representation},
        {"classification_and_selection",
         "nmopt.complementarity.classification_and_selection",
         {"backend-neutral", "contract", "complementarity", "active-set"},
         30,
         test_classification_and_selection},
        {"invalid_contract_inputs",
         "nmopt.complementarity.invalid_contract_inputs",
         {"backend-neutral", "contract", "complementarity", "negative"},
         30,
         test_invalid_contract_inputs}};
      const auto result = nmopt::test_support::run_requested_scenarios(
        argc, argv, scenarios, std::cout);
      if (!result.listed)
        std::cout << "backend-neutral complementarity contract test passed: "
                  << result.executed << '\n';
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "backend-neutral complementarity contract test failed: "
                << exception.what() << '\n';
      return 1;
    }
}
