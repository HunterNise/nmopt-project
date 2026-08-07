#include "nmopt/contract/linalg.hpp"
#include "nmopt/semantic/v1/problem_spec.hpp"
#include "test_support/diagnostics.hpp"
#include "test_support/scenario_dispatch.hpp"

#include <exception>
#include <iostream>
#include <string>
#include <vector>

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

    const auto dirichlet_control_specification =
      nmopt::semantic::v1::make_dirichlet_control_scalar_diffusion_reaction_problem();
    const auto dirichlet_control_report =
      validator.validate(dirichlet_control_specification);
    require(dirichlet_control_report.valid(),
            "the Dirichlet-control lifting v1 graph is invalid");

    const auto subdomain_specification =
      nmopt::semantic::v1::make_subdomain_tracking_scalar_diffusion_reaction_problem(
        1);
    const auto subdomain_report = validator.validate(subdomain_specification);
    require(subdomain_report.valid(),
            "the material-subdomain v1 tracking graph is invalid");

    const auto h1_control_specification =
      nmopt::semantic::v1::make_h1_regularised_scalar_diffusion_reaction_problem();
    const auto h1_control_report = validator.validate(h1_control_specification);
    require(h1_control_report.valid(),
            "the H1-control regularisation v1 graph is invalid");

    const auto h1_metric_specification =
      nmopt::semantic::v1::make_h1_metric_scalar_diffusion_reaction_problem();
    const auto h1_metric_report = validator.validate(h1_metric_specification);
    require(h1_metric_report.valid(),
            "the H1-control metric v1 graph is invalid");

    const auto coefficient_specification =
      nmopt::semantic::v1::make_coefficient_identification_problem();
    const auto coefficient_report = validator.validate(coefficient_specification);
    require(coefficient_report.valid(),
            "the coefficient-identification v1 graph is invalid");

    const auto boundary_specification =
      nmopt::semantic::v1::make_neumann_boundary_control_problem(true);
    const auto boundary_report = validator.validate(boundary_specification);
    require(boundary_report.valid(),
            "the Neumann boundary-control v1 graph is invalid");

    const auto pure_neumann_specification =
      nmopt::semantic::v1::make_pure_neumann_boundary_control_problem();
    const auto pure_neumann_report = validator.validate(pure_neumann_specification);
    require(pure_neumann_report.valid(),
            "the pure-Neumann mean-constraint v1 graph is invalid");

    auto missing_lifting_port = fixed_specification;
    missing_lifting_port.transformations.front().fixed_data_id = "missing_data";
    const auto lifting_port_report = validator.validate(missing_lifting_port);
    nmopt::test_support::require_exact_diagnostic(
      lifting_port_report,
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "fixed_dirichlet_reconstruction",
      "fixed_dirichlet_reconstruction_ports",
      "v1 semantic validation did not classify a broken lifting port");

    auto missing_dirichlet_control_port = dirichlet_control_specification;
    missing_dirichlet_control_port.transformations.front().control_variable_id =
      "missing_control";
    const auto dirichlet_control_port_report =
      validator.validate(missing_dirichlet_control_port);
    nmopt::test_support::require_exact_diagnostic(
      dirichlet_control_port_report,
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "dirichlet_control_lifting",
      "dirichlet_control_lifting_ports",
      "v1 semantic validation did not classify a broken Dirichlet-control lifting port");

    auto unused_reconstruction = fixed_specification;
    unused_reconstruction.variables.front().physical_field_transform_id.clear();
    const auto unused_reconstruction_report =
      validator.validate(unused_reconstruction);
    nmopt::test_support::require_exact_diagnostic(
      unused_reconstruction_report,
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "fixed_dirichlet_reconstruction",
      "physical_field_transformation_output",
      "v1 semantic validation did not classify an unused reconstruction");

    auto missing_policy = specification;
    missing_policy.requirement_policies.clear();
    const auto policy_report = validator.validate(missing_policy);
    nmopt::test_support::require_exact_diagnostic(
      policy_report,
      nmopt::semantic::v1::DiagnosticCategory::analytical_policy,
      "state",
      "state_uniqueness_realisation",
      "v1 semantic validation did not classify a missing policy");

    auto missing_target_rule = subdomain_specification;
    for (auto &policy : missing_target_rule.requirement_policies)
      if (policy.subject_id == "desired_state")
        policy.status = nmopt::semantic::v1::RequirementStatus::provided;
    const auto target_rule_report = validator.validate(missing_target_rule);
    nmopt::test_support::require_exact_diagnostic(
      target_rule_report,
      nmopt::semantic::v1::DiagnosticCategory::analytical_policy,
      "desired_state",
      "desired_state_data_rule",
      "v1 semantic validation did not require an explicit target-data rule");

    auto mismatched_target_region = subdomain_specification;
    for (auto &policy : mismatched_target_region.requirement_policies)
      if (policy.subject_id == "desired_state")
        policy.region_id = "domain";
    const auto target_region_report =
      validator.validate(mismatched_target_region);
    nmopt::test_support::require_exact_diagnostic(
      target_region_report,
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "state_tracking",
      "tracking_target_data_region",
      "v1 semantic validation did not match target data to its observation region");

    auto missing_neumann_trace_policy = boundary_specification;
    for (auto &policy : missing_neumann_trace_policy.requirement_policies)
      if (policy.subject_id == "neumann_control")
        policy.status = nmopt::semantic::v1::RequirementStatus::provided;
    const auto missing_neumann_trace_report =
      validator.validate(missing_neumann_trace_policy);
    nmopt::test_support::require_exact_diagnostic(
      missing_neumann_trace_report,
      nmopt::semantic::v1::DiagnosticCategory::analytical_policy,
      "neumann_control",
      "neumann_control_trace_realisation",
      "v1 semantic validation did not require the Neumann trace policy");

    auto missing_mean_constraint = pure_neumann_specification;
    for (auto &policy : missing_mean_constraint.requirement_policies)
      if (policy.kind == nmopt::semantic::v1::RequirementKind::mean_zero_multiplier)
        policy.status = nmopt::semantic::v1::RequirementStatus::provided;
    const auto missing_mean_constraint_report =
      validator.validate(missing_mean_constraint);
    nmopt::test_support::require_exact_diagnostic(
      missing_mean_constraint_report,
      nmopt::semantic::v1::DiagnosticCategory::analytical_policy,
      "state",
      "state_uniqueness_realisation",
      "v1 semantic validation did not require the pure-Neumann mean constraint");

    auto h1_regularisation_data_mismatch = h1_control_specification;
    h1_regularisation_data_mismatch.losses.at(1).data_id = "desired_state";
    const auto h1_regularisation_data_mismatch_report =
      validator.validate(h1_regularisation_data_mismatch);
    nmopt::test_support::require_exact_diagnostic(
      h1_regularisation_data_mismatch_report,
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "control_h1_regularisation",
      "loss_data_role",
      "v1 semantic validation did not validate H1 regularisation data");

    auto conflicting_state_gauges = pure_neumann_specification;
    conflicting_state_gauges.requirement_policies.push_back(
      {"state_fixed_dirichlet", "state",
       nmopt::semantic::v1::RequirementKind::fixed_dirichlet,
       nmopt::semantic::v1::RequirementStatus::selected_discrete_realisation,
       nmopt::semantic::v1::RequirementScope::discrete_compilation,
       "homogeneous full-vector Dirichlet rows", "control_boundary"});
    const auto conflicting_state_gauges_report =
      validator.validate(conflicting_state_gauges);
    nmopt::test_support::require_exact_diagnostic(
      conflicting_state_gauges_report,
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "state",
      "state_uniqueness_policy_conflict",
      "v1 semantic validation did not reject conflicting state gauges");

    auto mismatched_neumann_region = boundary_specification;
    for (auto &term : mismatched_neumann_region.residual_terms)
      if (term.id == "neumann_control")
        term.region_id = "observation_boundary";
    const auto mismatched_neumann_region_report =
      validator.validate(mismatched_neumann_region);
    nmopt::test_support::require_exact_diagnostic(
      mismatched_neumann_region_report,
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "neumann_control",
      "neumann_control_space_region",
      "v1 semantic validation did not match the Neumann control to its boundary space");

    auto missing_test_space = specification;
    missing_test_space.equations.front().test_space_id = "missing_test_space";
    const auto structural_report = validator.validate(missing_test_space);
    nmopt::test_support::require_exact_diagnostic(
      structural_report,
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "state_equation",
      "equation_test_space",
      "v1 semantic validation did not classify a broken structural port");
  }
} // namespace

int
main(const int argc, char **argv)
{
  try
    {
      const std::string executed =
        nmopt::test_support::run_requested_scenarios(
          argc,
          argv,
          {{"validation", test_semantic_v1_validation}});
      std::cout << "semantic v1 contract scenario passed: " << executed
                << '\n';
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "semantic v1 contract test failed: " << exception.what()
                << '\n';
      return 1;
    }
}
