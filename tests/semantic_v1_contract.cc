#include "nmopt/contract/linalg.hpp"
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

    const auto fixed_specification =
      nmopt::semantic::v1::make_fixed_dirichlet_scalar_diffusion_reaction_problem();
    const auto fixed_report = validator.validate(fixed_specification);
    require(fixed_report.valid(),
            "the fixed-Dirichlet v1 reconstruction graph is invalid");

    auto missing_lifting_port = fixed_specification;
    missing_lifting_port.transformations.front().fixed_data_id = "missing_data";
    const auto lifting_port_report = validator.validate(missing_lifting_port);
    require(lifting_port_report.has_category(
              nmopt::semantic::v1::DiagnosticCategory::structural),
            "v1 semantic validation did not classify a broken lifting port");

    auto unused_reconstruction = fixed_specification;
    unused_reconstruction.variables.front().physical_field_transform_id.clear();
    const auto unused_reconstruction_report =
      validator.validate(unused_reconstruction);
    require(unused_reconstruction_report.has_category(
              nmopt::semantic::v1::DiagnosticCategory::structural),
            "v1 semantic validation did not classify an unused reconstruction");

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
