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

    const auto subdomain_specification =
      nmopt::semantic::v1::make_subdomain_tracking_scalar_diffusion_reaction_problem(
        1);
    const auto subdomain_report = validator.validate(subdomain_specification);
    require(subdomain_report.valid(),
            "the material-subdomain v1 tracking graph is invalid");

    const auto boundary_specification =
      nmopt::semantic::v1::make_neumann_boundary_control_problem(true);
    const auto boundary_report = validator.validate(boundary_specification);
    require(boundary_report.valid(),
            "the Neumann boundary-control v1 graph is invalid");

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

    auto missing_target_rule = subdomain_specification;
    for (auto &policy : missing_target_rule.requirement_policies)
      if (policy.subject_id == "desired_state")
        policy.status = nmopt::semantic::v1::RequirementStatus::provided;
    const auto target_rule_report = validator.validate(missing_target_rule);
    require(target_rule_report.has_category(
              nmopt::semantic::v1::DiagnosticCategory::analytical_policy),
            "v1 semantic validation did not require an explicit target-data rule");

    auto mismatched_target_region = subdomain_specification;
    for (auto &policy : mismatched_target_region.requirement_policies)
      if (policy.subject_id == "desired_state")
        policy.region_id = "domain";
    const auto target_region_report =
      validator.validate(mismatched_target_region);
    require(target_region_report.has_category(
              nmopt::semantic::v1::DiagnosticCategory::structural),
            "v1 semantic validation did not match target data to its observation region");

    auto missing_neumann_trace_policy = boundary_specification;
    for (auto &policy : missing_neumann_trace_policy.requirement_policies)
      if (policy.subject_id == "neumann_control")
        policy.status = nmopt::semantic::v1::RequirementStatus::provided;
    const auto missing_neumann_trace_report =
      validator.validate(missing_neumann_trace_policy);
    require(missing_neumann_trace_report.has_category(
              nmopt::semantic::v1::DiagnosticCategory::analytical_policy),
            "v1 semantic validation did not require the Neumann trace policy");

    auto mismatched_neumann_region = boundary_specification;
    for (auto &term : mismatched_neumann_region.residual_terms)
      if (term.id == "neumann_control")
        term.region_id = "observation_boundary";
    const auto mismatched_neumann_region_report =
      validator.validate(mismatched_neumann_region);
    require(mismatched_neumann_region_report.has_category(
              nmopt::semantic::v1::DiagnosticCategory::structural),
            "v1 semantic validation did not match the Neumann control to its boundary space");

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
