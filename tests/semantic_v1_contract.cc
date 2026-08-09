#include "nmopt/contract/linalg.hpp"
#include "nmopt/compiler/v1/dealii_scalar_plan.hpp"
#include "nmopt/semantic/v1/problem_spec.hpp"
#include "test_support/diagnostics.hpp"
#include "test_support/scenario_dispatch.hpp"

#include <algorithm>
#include <exception>
#include <functional>
#include <iostream>
#include <set>
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

  template <typename Component>
  Component &
  component_by_id(std::vector<Component> &components, const std::string &id)
  {
    const auto component = std::find_if(
      components.begin(), components.end(), [&id](const Component &candidate) {
        return candidate.id == id;
      });
    require(component != components.end(), "semantic test component is missing");
    require(std::count_if(components.begin(),
                          components.end(),
                          [&id](const Component &candidate) {
                            return candidate.id == id;
                          }) == 1,
            "semantic test component is not unique");
    return *component;
  }

  template <typename Component>
  std::vector<std::string>
  sorted_component_ids(const std::vector<Component> &components)
  {
    std::vector<std::string> ids;
    ids.reserve(components.size());
    for (const auto &component : components)
      ids.push_back(component.id);
    std::sort(ids.begin(), ids.end());
    return ids;
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

    const auto h1_state_specification =
      nmopt::semantic::v1::make_h1_state_tracking_scalar_diffusion_reaction_problem();
    const auto h1_state_report = validator.validate(h1_state_specification);
    require(h1_state_report.valid(),
            "the H1-state tracking v1 graph is invalid");

    auto h1_state_l2_output = h1_state_specification;
    component_by_id(h1_state_l2_output.spaces, "state_observation_space")
      .topology = nmopt::semantic::v1::SpaceTopology::l2;
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(h1_state_l2_output),
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "state_observation",
      "h1_state_restriction_output_topology",
      "v1 semantic validation accepted an L2 output for the H1 state observation");

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

    auto weighted_boundary_specification =
      nmopt::semantic::v1::
        make_weighted_boundary_trace_neumann_control_problem(true);
    const auto weighted_boundary_report =
      validator.validate(weighted_boundary_specification);
    require(weighted_boundary_report.valid(),
            "the weighted boundary-trace v1 graph is invalid");
    require(component_by_id(weighted_boundary_specification.observations,
                            "weighted_state_boundary_trace")
                .data_ids == std::vector<std::string>{"boundary_weight"},
            "the weighted boundary trace omitted its immutable data port");

    auto missing_weight_port = weighted_boundary_specification;
    component_by_id(missing_weight_port.observations,
                    "weighted_state_boundary_trace")
      .data_ids.clear();
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(missing_weight_port),
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "weighted_state_boundary_trace",
      "weighted_boundary_trace_data_port",
      "v1 semantic validation accepted a weighted trace without weight data");

    auto missing_weight_policy = weighted_boundary_specification;
    for (auto &policy : missing_weight_policy.requirement_policies)
      if (policy.subject_id == "boundary_weight")
        policy.status = nmopt::semantic::v1::RequirementStatus::provided;
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(missing_weight_policy),
      nmopt::semantic::v1::DiagnosticCategory::analytical_policy,
      "boundary_weight",
      "observation_weight_data_rule",
      "v1 semantic validation did not require weight quadrature provenance");

    auto missing_weight_boundedness = weighted_boundary_specification;
    component_by_id(missing_weight_boundedness.requirement_policies,
                    "boundary_weight_boundedness")
      .status = nmopt::semantic::v1::RequirementStatus::unspecified;
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(missing_weight_boundedness),
      nmopt::semantic::v1::DiagnosticCategory::analytical_policy,
      "boundary_weight",
      "observation_weight_boundedness",
      "v1 semantic validation did not require bounded observation weight data");

    auto mismatched_weight_region = weighted_boundary_specification;
    for (auto &policy : mismatched_weight_region.requirement_policies)
      if (policy.subject_id == "boundary_weight")
        policy.region_id = "control_boundary";
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(mismatched_weight_region),
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "weighted_state_boundary_trace",
      "weighted_boundary_trace_data_region",
      "v1 semantic validation did not match the weight to the trace region");

    const auto pure_neumann_specification =
      nmopt::semantic::v1::make_pure_neumann_boundary_control_problem();
    const auto pure_neumann_report = validator.validate(pure_neumann_specification);
    require(pure_neumann_report.valid(),
            "the pure-Neumann mean-constraint v1 graph is invalid");

    const auto general_scalar_specification =
      nmopt::semantic::v1::make_general_scalar_elliptic_robin_problem(
        {0, 2, 3}, {1});
    const auto general_scalar_report =
      validator.validate(general_scalar_specification);
    require(general_scalar_report.valid(),
            "the P5.1 general scalar Robin graph is invalid");

    auto wrong_tensor_shape = general_scalar_specification;
    component_by_id(wrong_tensor_shape.data, "diffusion_tensor").kind =
      nmopt::semantic::v1::DataKind::vector_function;
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(wrong_tensor_shape),
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "diffusion_tensor",
      "coefficient_data_shape",
      "v1 semantic validation accepted a vector-valued tensor diffusion binding");

    auto missing_conormal_policy = general_scalar_specification;
    for (auto &policy : missing_conormal_policy.requirement_policies)
      if (policy.kind ==
          nmopt::semantic::v1::RequirementKind::conormal_flux)
        policy.status = nmopt::semantic::v1::RequirementStatus::provided;
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(missing_conormal_policy),
      nmopt::semantic::v1::DiagnosticCategory::analytical_policy,
      "robin_bilinear",
      "conormal_flux_convention",
      "v1 semantic validation did not require a conormal-flux convention");

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

  void
  test_semantic_v1_graph_closure()
  {
    using namespace nmopt::semantic::v1;
    const SemanticValidator validator;

    auto orphan_term = make_scalar_diffusion_reaction_problem(true);
    orphan_term.residual_terms.push_back(
      {"orphan_diffusion_reaction", "Orphan diffusion and reaction",
       ResidualTermKind::diffusion_reaction, "state_equation", {"state"},
       {"diffusion", "reaction"}, ""});
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(orphan_term),
      DiagnosticCategory::structural,
      "orphan_diffusion_reaction",
      "residual_term_equation_membership",
      "v1 semantic validation accepted an orphan residual term");

    auto duplicate_edge = make_scalar_diffusion_reaction_problem(true);
    component_by_id(duplicate_edge.equations, "state_equation")
      .residual_term_ids.push_back("volume_source");
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(duplicate_edge),
      DiagnosticCategory::structural,
      "state_equation",
      "unique_equation_residual_term_edges",
      "v1 semantic validation accepted a duplicate equation edge");

    auto wrong_bound_space = make_scalar_diffusion_reaction_problem(true);
    component_by_id(wrong_bound_space.data, "lower_bound").space_id =
      "state_space";
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(wrong_bound_space),
      DiagnosticCategory::structural,
      "control_box",
      "constraint_bound_data_space",
      "v1 semantic validation accepted bound data in the wrong space");

    auto missing_label = make_scalar_diffusion_reaction_problem(true);
    component_by_id(missing_label.variables, "state").label.clear();
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(missing_label),
      DiagnosticCategory::structural,
      "state",
      "human_readable_label",
      "v1 semantic validation accepted a missing component label");

    auto mismatched_metric = make_scalar_diffusion_reaction_problem(true);
    mismatched_metric.metrics.push_back(
      {"state_l2_metric", "State L2 metric", MetricKind::l2, "state",
       "state_pairing"});
    mismatched_metric.formulation.metric_id = "state_l2_metric";
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(mismatched_metric),
      DiagnosticCategory::structural,
      "reduced_dto",
      "formulation_metric_variable",
      "v1 semantic validation accepted a metric for another variable");

    auto mismatched_constraint = make_scalar_diffusion_reaction_problem(true);
    mismatched_constraint.variables.push_back(
      {"other_control", "Other control", VariableRole::control,
       "control_space", ""});
    mismatched_constraint.constraints.push_back(
      {"other_control_box", "Other control box", ConstraintKind::cellwise_box,
       "other_control", "lower_bound", "upper_bound"});
    mismatched_constraint.formulation.constraint_id = "other_control_box";
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(mismatched_constraint),
      DiagnosticCategory::structural,
      "reduced_dto",
      "formulation_constraint_variable",
      "v1 semantic validation accepted a constraint for another variable");
  }

  void
  test_semantic_v1_pairing_compatibility()
  {
    using namespace nmopt::semantic::v1;
    const SemanticValidator validator;

    auto equation_pairing = make_scalar_diffusion_reaction_problem();
    component_by_id(equation_pairing.pairings, "state_test_pairing")
      .covector_space_id = "control_space";
    const auto equation_report = validator.validate(equation_pairing);
    nmopt::test_support::require_exact_diagnostic(
      equation_report,
      DiagnosticCategory::structural,
      "state_test_pairing",
      "pairing_primal_covector_space",
      "v1 semantic validation accepted incompatible pairing ports");
    nmopt::test_support::require_exact_diagnostic(
      equation_report,
      DiagnosticCategory::structural,
      "state_equation",
      "equation_test_pairing",
      "v1 equation validation ignored the pairing covector port");

    auto observation_pairing = make_scalar_diffusion_reaction_problem();
    component_by_id(observation_pairing.pairings,
                    "state_observation_pairing")
      .covector_space_id = "control_space";
    const auto observation_report = validator.validate(observation_pairing);
    nmopt::test_support::require_exact_diagnostic(
      observation_report,
      DiagnosticCategory::structural,
      "state_observation",
      "observation_output_pairing",
      "v1 observation validation ignored the pairing covector port");
    nmopt::test_support::require_exact_diagnostic(
      observation_report,
      DiagnosticCategory::structural,
      "state_tracking",
      "loss_pairing",
      "v1 loss validation ignored the pairing covector port");

    auto metric_pairing = make_scalar_diffusion_reaction_problem();
    component_by_id(metric_pairing.pairings, "control_pairing")
      .covector_space_id = "state_space";
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(metric_pairing),
      DiagnosticCategory::structural,
      "control_l2_metric",
      "metric_pairing",
      "v1 metric validation ignored the pairing covector port");
  }

  void
  test_semantic_v1_incomplete_components()
  {
    using namespace nmopt::semantic::v1;
    require(RegionSpec{}.kind == RegionKind::unspecified,
            "default region kind is not safe");
    require(SpaceSpec{}.topology == SpaceTopology::unspecified &&
              SpaceSpec{}.role == SpaceRole::unspecified,
            "default space enums are not safe");
    require(VariableSpec{}.role == VariableRole::unspecified,
            "default variable role is not safe");
    require(DataSpec{}.kind == DataKind::unspecified &&
              DataSpec{}.role == DataRole::unspecified,
            "default data enums are not safe");
    require(TransformationSpec{}.kind == TransformationKind::unspecified,
            "default transformation kind is not safe");
    require(ResidualTermSpec{}.kind == ResidualTermKind::unspecified,
            "default residual-term kind is not safe");
    require(ObservationSpec{}.kind == ObservationKind::unspecified,
            "default observation kind is not safe");
    require(LossSpec{}.kind == LossKind::unspecified,
            "default loss kind is not safe");
    require(MetricSpec{}.kind == MetricKind::unspecified,
            "default metric kind is not safe");
    require(ConstraintSpec{}.kind == ConstraintKind::unspecified,
            "default constraint kind is not safe");
    require(RequirementPolicySpec{}.kind == RequirementKind::unspecified &&
              RequirementPolicySpec{}.status == RequirementStatus::unspecified &&
              RequirementPolicySpec{}.scope == RequirementScope::unspecified,
            "default requirement-policy enums are not safe");
    require(ReducedFormulationSpec{}.kind == FormulationKind::unspecified,
            "default formulation kind is not safe");

    const SemanticValidator validator;
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(ProblemSpec{}),
      DiagnosticCategory::structural,
      "formulation",
      "formulation_kind",
      "default semantic problem did not diagnose an incomplete formulation");

    auto default_pairing = make_scalar_diffusion_reaction_problem();
    default_pairing.pairings.push_back(PairingSpec{});
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(default_pairing),
      DiagnosticCategory::structural,
      "pairing",
      "stable_component_identity",
      "default semantic pairing was not diagnosed safely");

    auto default_equation = make_scalar_diffusion_reaction_problem();
    default_equation.equations.push_back(EquationBlockSpec{});
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(default_equation),
      DiagnosticCategory::structural,
      "equation",
      "stable_component_identity",
      "default semantic equation was not diagnosed safely");

    struct IncompleteCase
    {
      std::string                        component_id;
      std::string                        capability;
      std::function<void(ProblemSpec &)> make_incomplete;
    };

    const std::vector<IncompleteCase> cases{
      {"domain", "region_kind", [](ProblemSpec &specification) {
         component_by_id(specification.regions, "domain").kind =
           RegionKind::unspecified;
       }},
      {"state_space", "space_topology", [](ProblemSpec &specification) {
         component_by_id(specification.spaces, "state_space").topology =
           SpaceTopology::unspecified;
       }},
      {"state_space", "space_role", [](ProblemSpec &specification) {
         component_by_id(specification.spaces, "state_space").role =
           SpaceRole::unspecified;
       }},
      {"state", "variable_role", [](ProblemSpec &specification) {
         component_by_id(specification.variables, "state").role =
           VariableRole::unspecified;
       }},
      {"forcing", "data_kind", [](ProblemSpec &specification) {
         component_by_id(specification.data, "forcing").kind =
           DataKind::unspecified;
       }},
      {"forcing", "data_role", [](ProblemSpec &specification) {
         component_by_id(specification.data, "forcing").role =
           DataRole::unspecified;
       }},
      {"fixed_dirichlet_reconstruction",
       "transformation_kind",
       [](ProblemSpec &specification) {
         specification =
           make_fixed_dirichlet_scalar_diffusion_reaction_problem();
         component_by_id(specification.transformations,
                         "fixed_dirichlet_reconstruction")
           .kind = TransformationKind::unspecified;
       }},
      {"diffusion_reaction", "residual_term_kind", [](ProblemSpec &specification) {
         component_by_id(specification.residual_terms, "diffusion_reaction").kind =
           ResidualTermKind::unspecified;
       }},
      {"state_observation", "observation_kind", [](ProblemSpec &specification) {
         component_by_id(specification.observations, "state_observation").kind =
           ObservationKind::unspecified;
       }},
      {"state_tracking", "loss_kind", [](ProblemSpec &specification) {
         component_by_id(specification.losses, "state_tracking").kind =
           LossKind::unspecified;
       }},
      {"control_l2_metric", "metric_kind", [](ProblemSpec &specification) {
         component_by_id(specification.metrics, "control_l2_metric").kind =
           MetricKind::unspecified;
       }},
      {"control_box", "constraint_kind", [](ProblemSpec &specification) {
         component_by_id(specification.constraints, "control_box").kind =
           ConstraintKind::unspecified;
       }},
      {"state_fixed_dirichlet", "requirement_kind", [](ProblemSpec &specification) {
         component_by_id(specification.requirement_policies,
                         "state_fixed_dirichlet")
           .kind = RequirementKind::unspecified;
       }},
      {"state_fixed_dirichlet", "requirement_status", [](ProblemSpec &specification) {
         component_by_id(specification.requirement_policies,
                         "state_fixed_dirichlet")
           .status = RequirementStatus::unspecified;
       }},
      {"state_fixed_dirichlet", "requirement_scope", [](ProblemSpec &specification) {
         component_by_id(specification.requirement_policies,
                         "state_fixed_dirichlet")
           .scope = RequirementScope::unspecified;
       }},
      {"reduced_dto", "formulation_kind", [](ProblemSpec &specification) {
         specification.formulation.kind = FormulationKind::unspecified;
       }}};

    for (const auto &test_case : cases)
      {
        ProblemSpec specification = make_scalar_diffusion_reaction_problem(true);
        test_case.make_incomplete(specification);
        nmopt::test_support::require_exact_diagnostic(
          validator.validate(specification),
          DiagnosticCategory::structural,
          test_case.component_id,
          test_case.capability,
          "partially populated semantic component was not diagnosed");
      }
  }

  void
  test_semantic_v1_reference_delta_stability()
  {
    using namespace nmopt::semantic::v1;
    ProblemSpec reordered = make_scalar_diffusion_reaction_problem();
    std::reverse(reordered.regions.begin(), reordered.regions.end());
    std::reverse(reordered.spaces.begin(), reordered.spaces.end());
    std::reverse(reordered.pairings.begin(), reordered.pairings.end());
    std::reverse(reordered.variables.begin(), reordered.variables.end());
    std::reverse(reordered.data.begin(), reordered.data.end());
    std::reverse(reordered.residual_terms.begin(), reordered.residual_terms.end());
    std::reverse(reordered.equations.begin(), reordered.equations.end());
    std::reverse(reordered.observations.begin(), reordered.observations.end());
    std::reverse(reordered.losses.begin(), reordered.losses.end());
    std::reverse(reordered.metrics.begin(), reordered.metrics.end());
    std::reverse(reordered.requirement_policies.begin(),
                 reordered.requirement_policies.end());

    reference_detail::apply_coefficient_identification_delta(reordered);
    const SemanticValidator validator;
    require(validator.validate(reordered).valid(),
            "an ID-based feature delta depends on declaration order");

    const ProblemSpec expected = make_coefficient_identification_problem();
    require(sorted_component_ids(reordered.regions) ==
              sorted_component_ids(expected.regions) &&
              sorted_component_ids(reordered.spaces) ==
                sorted_component_ids(expected.spaces) &&
              sorted_component_ids(reordered.pairings) ==
                sorted_component_ids(expected.pairings) &&
              sorted_component_ids(reordered.variables) ==
                sorted_component_ids(expected.variables) &&
              sorted_component_ids(reordered.data) ==
                sorted_component_ids(expected.data) &&
              sorted_component_ids(reordered.residual_terms) ==
                sorted_component_ids(expected.residual_terms) &&
              sorted_component_ids(reordered.equations) ==
                sorted_component_ids(expected.equations) &&
              sorted_component_ids(reordered.observations) ==
                sorted_component_ids(expected.observations) &&
              sorted_component_ids(reordered.losses) ==
                sorted_component_ids(expected.losses) &&
              sorted_component_ids(reordered.metrics) ==
                sorted_component_ids(expected.metrics) &&
              sorted_component_ids(reordered.constraints) ==
                sorted_component_ids(expected.constraints) &&
              sorted_component_ids(reordered.requirement_policies) ==
                sorted_component_ids(expected.requirement_policies),
            "a reordered feature delta changed semantic component identities");
    require(reordered.formulation.state_variable_id ==
                expected.formulation.state_variable_id &&
              reordered.formulation.control_variable_id ==
                expected.formulation.control_variable_id &&
              reordered.formulation.equation_id ==
                expected.formulation.equation_id &&
              reordered.formulation.metric_id ==
                expected.formulation.metric_id &&
              reordered.formulation.constraint_id ==
                expected.formulation.constraint_id,
            "a reordered feature delta changed formulation ports");

    ProblemSpec reordered_boundary = make_neumann_boundary_control_problem();
    std::reverse(reordered_boundary.data.begin(), reordered_boundary.data.end());
    std::reverse(reordered_boundary.observations.begin(),
                 reordered_boundary.observations.end());
    std::reverse(reordered_boundary.losses.begin(),
                 reordered_boundary.losses.end());
    std::reverse(reordered_boundary.requirement_policies.begin(),
                 reordered_boundary.requirement_policies.end());
    reference_detail::apply_weighted_boundary_trace_delta(reordered_boundary);
    const ProblemSpec expected_weighted =
      make_weighted_boundary_trace_neumann_control_problem();
    require(validator.validate(reordered_boundary).valid() &&
              sorted_component_ids(reordered_boundary.data) ==
                sorted_component_ids(expected_weighted.data) &&
              sorted_component_ids(reordered_boundary.observations) ==
                sorted_component_ids(expected_weighted.observations) &&
              sorted_component_ids(reordered_boundary.losses) ==
                sorted_component_ids(expected_weighted.losses) &&
              sorted_component_ids(reordered_boundary.requirement_policies) ==
                sorted_component_ids(expected_weighted.requirement_policies) &&
              component_by_id(reordered_boundary.observations,
                              "weighted_state_boundary_trace")
                  .data_ids == std::vector<std::string>{"boundary_weight"},
            "the weighted-trace delta depends on declaration order");
  }

  void
  test_semantic_v1_resolution()
  {
    auto specification =
      nmopt::semantic::v1::make_scalar_diffusion_reaction_problem(true);
    const nmopt::semantic::v1::SemanticResolver resolver;
    const auto resolution = resolver.resolve(specification);
    require(resolution.succeeded(),
            "the canonical semantic graph did not produce a resolved view");
    const auto &view = *resolution.problem;
    require(view.specification().id == specification.id &&
              view.variable("state").space_id == "state_space" &&
              view.variable("control").space_id == "control_space" &&
              view.equation("state_equation").residual_term_ids.size() == 3 &&
              view.residual_term("volume_control").kind ==
                nmopt::semantic::v1::ResidualTermKind::volume_control &&
              view.observation("state_observation").region_id == "domain" &&
              view.metric("control_l2_metric").variable_id == "control" &&
              view.constraint("control_box").variable_id == "control",
            "the resolved semantic view did not preserve stable-ID edges");

    std::reverse(specification.regions.begin(), specification.regions.end());
    std::reverse(specification.spaces.begin(), specification.spaces.end());
    std::reverse(specification.residual_terms.begin(),
                 specification.residual_terms.end());
    std::reverse(specification.observations.begin(),
                 specification.observations.end());
    const auto reordered = resolver.resolve(specification);
    require(reordered.succeeded() &&
              reordered.problem->region("domain").is_full_domain &&
              reordered.problem->residual_term("diffusion_reaction").data_ids ==
                std::vector<std::string>({"diffusion", "reaction"}),
            "semantic resolution depended on declaration order");

    specification.spaces.push_back(specification.spaces.front());
    const auto duplicate = resolver.resolve(specification);
    require(!duplicate.succeeded() && !duplicate.problem.has_value(),
            "an invalid graph produced a resolved semantic view");
    nmopt::test_support::require_exact_diagnostic(
      duplicate.diagnostics,
      nmopt::semantic::v1::DiagnosticCategory::structural,
      specification.spaces.front().id,
      "unique_component_identity",
      "semantic resolver did not retain duplicate-ID diagnostics");
  }

  void
  test_dealii_scalar_lowering_plan()
  {
    const nmopt::semantic::v1::SemanticResolver resolver;
    const nmopt::compiler::v1::DealiiScalarLoweringPlanner planner;
    const auto specification =
      nmopt::semantic::v1::make_scalar_diffusion_reaction_problem(true);
    const auto resolution = resolver.resolve(specification);
    require(resolution.succeeded(),
            "scalar lowering-plan setup did not resolve");
    const auto planned = planner.plan(*resolution.problem);
    require(planned.succeeded(),
            "canonical graph did not produce a scalar lowering plan");
    require(planned.plan->residual_terms.size() == 3 &&
              planned.plan->observations.size() == 2 &&
              planned.plan->losses.size() == 2 &&
              planned.plan->constraint ==
                nmopt::compiler::v1::ScalarConstraintOperatorKind::cellwise_box &&
              planned.plan->transformation ==
                nmopt::compiler::v1::ScalarTransformationOperatorKind::none &&
              planned.plan->dirichlet_boundary_ids ==
                std::set<unsigned int>{0} &&
              planned.plan->provenance.size() == 9,
            "canonical scalar lowering plan omitted component contributions");

    const auto fixed_specification =
      nmopt::semantic::v1::make_fixed_dirichlet_scalar_diffusion_reaction_problem();
    const auto fixed_resolution = resolver.resolve(fixed_specification);
    require(fixed_resolution.succeeded(),
            "fixed reconstruction lowering-plan setup did not resolve");
    const auto fixed_plan = planner.plan(*fixed_resolution.problem);
    require(fixed_plan.succeeded() &&
              fixed_plan.plan->transformation ==
                nmopt::compiler::v1::ScalarTransformationOperatorKind::fixed_dirichlet_reconstruction &&
              fixed_plan.plan->fixed_data_id == "fixed_dirichlet_data",
            "fixed reconstruction did not contribute its scalar strategy");

    const auto boundary_specification =
      nmopt::semantic::v1::make_neumann_boundary_control_problem();
    const auto boundary_resolution = resolver.resolve(boundary_specification);
    require(boundary_resolution.succeeded(),
            "specialized-boundary lowering-plan setup did not resolve");
    const auto boundary_plan = planner.plan(*boundary_resolution.problem);
    require(!boundary_plan.succeeded() && !boundary_plan.plan.has_value(),
            "specialized Neumann graph entered the bounded scalar plan");
    nmopt::test_support::require_exact_diagnostic(
      boundary_plan.diagnostics,
      nmopt::semantic::v1::DiagnosticCategory::lowerability,
      "neumann_control",
      "scalar_residual_component_lowerer",
      "scalar planner did not identify its specialized Neumann boundary");

    const auto general_specification =
      nmopt::semantic::v1::make_general_scalar_elliptic_robin_problem(
        {0, 2, 3}, {1});
    const auto general_resolution = resolver.resolve(general_specification);
    require(general_resolution.succeeded(),
            "general scalar lowering-plan setup did not resolve");
    const auto general_plan = planner.plan(*general_resolution.problem);
    require(general_plan.succeeded() &&
              general_plan.plan->residual_terms.size() == 8 &&
              general_plan.plan->robin_boundary_ids ==
                std::set<unsigned int>{1} &&
              general_plan.plan->provenance.size() == 13 &&
              std::any_of(
                general_plan.plan->residual_terms.begin(),
                general_plan.plan->residual_terms.end(),
                [](const nmopt::compiler::v1::ScalarResidualContribution &term) {
                  return term.operator_kind ==
                         nmopt::compiler::v1::ScalarResidualOperatorKind::conservative_transport;
                }),
            "P5.1 general scalar plan omitted a term or Robin boundary contribution");

    const auto h1_state_specification =
      nmopt::semantic::v1::make_h1_state_tracking_scalar_diffusion_reaction_problem();
    const auto h1_state_resolution = resolver.resolve(h1_state_specification);
    require(h1_state_resolution.succeeded(),
            "H1-state observation lowering-plan setup did not resolve");
    const auto h1_state_plan = planner.plan(*h1_state_resolution.problem);
    require(
      h1_state_plan.succeeded() &&
        std::any_of(
          h1_state_plan.plan->observations.begin(),
          h1_state_plan.plan->observations.end(),
          [](const nmopt::compiler::v1::ScalarObservationContribution &observation) {
            return observation.operator_kind ==
                   nmopt::compiler::v1::ScalarObservationOperatorKind::h1_state_restriction;
          }) &&
        std::find(h1_state_plan.plan->provenance.begin(),
                  h1_state_plan.plan->provenance.end(),
                  "state_observation <- "
                  "dealii.scalar.observation.h1_state_restriction") !=
          h1_state_plan.plan->provenance.end(),
      "P5.2 H1-state observation did not contribute its scalar lowering handler");
  }

} // namespace

int
main(const int argc, char **argv)
{
  try
    {
      const std::vector<nmopt::test_support::Scenario> scenarios{
        {"validation",
         "nmopt.semantic.v1_validation",
         {"backend-neutral", "semantic"},
         30,
         test_semantic_v1_validation},
        {"graph_closure",
         "nmopt.semantic.v1_graph_closure",
         {"backend-neutral", "semantic"},
         30,
         test_semantic_v1_graph_closure},
        {"pairing_compatibility",
         "nmopt.semantic.v1_pairing_compatibility",
         {"backend-neutral", "semantic"},
         30,
         test_semantic_v1_pairing_compatibility},
        {"incomplete_components",
         "nmopt.semantic.v1_incomplete_components",
         {"backend-neutral", "semantic"},
         30,
         test_semantic_v1_incomplete_components},
        {"reference_delta_stability",
         "nmopt.semantic.v1_reference_delta_stability",
         {"backend-neutral", "semantic"},
         30,
         test_semantic_v1_reference_delta_stability},
        {"resolution",
         "nmopt.semantic.v1_resolution",
         {"backend-neutral", "semantic", "compiler"},
         30,
         test_semantic_v1_resolution},
        {"scalar_lowering_plan",
         "nmopt.semantic.v1_scalar_lowering_plan",
         {"backend-neutral", "semantic", "compiler"},
         30,
         test_dealii_scalar_lowering_plan}};
      const auto result =
        nmopt::test_support::run_requested_scenarios(
          argc, argv, scenarios, std::cout);
      if (!result.listed)
        std::cout << "semantic v1 contract scenario passed: "
                  << result.executed << '\n';
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "semantic v1 contract test failed: " << exception.what()
                << '\n';
      return 1;
    }
}
