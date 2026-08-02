#include "nmopt/semantic/v1/problem_spec.hpp"

#include <exception>
#include <iostream>

namespace
{
  void
  require(const bool condition, const char *message)
  {
    if (!condition)
      throw nmopt::contract::ContractError(message);
  }

  void
  test_semantic_v1_validation()
  {
    const auto specification =
      nmopt::semantic::v1::make_scalar_diffusion_reaction_problem(true);
    const nmopt::semantic::v1::SemanticValidator validator;
    const auto valid_report = validator.validate(specification);
    require(valid_report.valid(),
            "the canonical v1 scalar diffusion-reaction graph is invalid");

    auto missing_policy = specification;
    missing_policy.requirement_policies.clear();
    const auto policy_report = validator.validate(missing_policy);
    require(policy_report.has_category(
              nmopt::semantic::v1::DiagnosticCategory::analytical_policy),
            "v1 semantic validation did not classify a missing policy");

    auto missing_test_space = specification;
    missing_test_space.equations.front().test_space_id = "missing_test_space";
    const auto structural_report = validator.validate(missing_test_space);
    require(structural_report.has_category(
              nmopt::semantic::v1::DiagnosticCategory::structural),
            "v1 semantic validation did not classify a broken structural port");
  }
} // namespace

int
main()
{
  try
    {
      test_semantic_v1_validation();
      std::cout << "semantic v1 contract test passed\n";
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "semantic v1 contract test failed: " << exception.what()
                << '\n';
      return 1;
    }
}
