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
  const Component &
  component_by_id(const std::vector<Component> &components,
                  const std::string            &id)
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
    using namespace nmopt::semantic::v1;
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

    auto l2_dirichlet_specification =
      nmopt::semantic::v1::make_l2_dirichlet_laplace_control_problem();
    require(validator.validate(l2_dirichlet_specification).valid(),
            "the L2 Dirichlet transposition graph is invalid");
    require(component_by_id(l2_dirichlet_specification.spaces, "state_space")
                .topology == nmopt::semantic::v1::SpaceTopology::l2 &&
              component_by_id(l2_dirichlet_specification.spaces,
                              "state_test_space")
                  .topology == nmopt::semantic::v1::SpaceTopology::h2 &&
              component_by_id(l2_dirichlet_specification.spaces,
                              "control_space")
                  .topology == nmopt::semantic::v1::SpaceTopology::l2 &&
              l2_dirichlet_specification.transformations.empty() &&
              component_by_id(l2_dirichlet_specification.residual_terms,
                              "transposition_state_action")
                  .kind == nmopt::semantic::v1::ResidualTermKind::
                             transposition_laplacian &&
              component_by_id(l2_dirichlet_specification.residual_terms,
                              "dirichlet_transposition_control")
                  .kind == nmopt::semantic::v1::ResidualTermKind::
                             dirichlet_transposition_control,
            "the L2 Dirichlet graph did not declare its exact continuous residual");

    const auto remove_policy = [](auto &graph,
                                  const std::string &policy_id) {
      graph.requirement_policies.erase(
        std::remove_if(graph.requirement_policies.begin(),
                       graph.requirement_policies.end(),
                       [&policy_id](const auto &policy) {
                         return policy.id == policy_id;
                       }),
        graph.requirement_policies.end());
    };

    auto point_sensor_specification =
      nmopt::semantic::v1::make_point_sensor_scalar_diffusion_reaction_problem(
        {{0.25, 0.35}, {0.7, 0.6}});
    require(validator.validate(point_sensor_specification).valid(),
            "the C5.10 point-sensor v1 graph is invalid");
    const auto point_sensor_region =
      component_by_id(point_sensor_specification.regions,
                      "point_sensor_region");
    const auto point_sensor_space =
      component_by_id(point_sensor_specification.spaces,
                      "state_observation_space");
    const auto point_sensor_observation =
      component_by_id(point_sensor_specification.observations,
                      "state_observation");
    require(point_sensor_region.kind ==
              nmopt::semantic::v1::RegionKind::point_set &&
              point_sensor_region.point_coordinates.size() == 2 &&
              point_sensor_space.dimension == 2 &&
              point_sensor_observation.kind ==
              nmopt::semantic::v1::ObservationKind::point_sensor,
            "the point-sensor graph did not pair physical points with a finite output space");
    const auto custom_point_sensor_specification =
      nmopt::semantic::v1::make_point_sensor_scalar_diffusion_reaction_problem(
        {{0.25, 0.35}}, {7});
    require(validator.validate(custom_point_sensor_specification).valid() &&
              component_by_id(custom_point_sensor_specification.regions,
                              "dirichlet_boundary")
                  .boundary_ids == std::vector<unsigned int>{7},
            "the point-sensor factory did not preserve custom fixed boundary ids");
    auto duplicate_point_fixed_boundary = point_sensor_specification;
    component_by_id(duplicate_point_fixed_boundary.regions,
                    "dirichlet_boundary")
      .boundary_ids = {0, 0};
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(duplicate_point_fixed_boundary),
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "dirichlet_boundary",
      "p53_fixed_dirichlet_boundary_ids",
      "the point-sensor graph accepted duplicate fixed boundary ids");
    const auto point_transposition_policy = component_by_id(
      point_sensor_specification.requirement_policies,
      "point_sensor_transposition_policy");
    require(point_transposition_policy.typed_transposition_selection.has_value() &&
              point_transposition_policy.typed_transposition_selection
                  ->strong_space_id == "transposition_strong_space" &&
              point_transposition_policy.typed_transposition_selection
                  ->observation_id == "state_observation" &&
              point_transposition_policy.typed_transposition_selection
                  ->diffusion_data_id == "diffusion" &&
              point_transposition_policy.typed_transposition_selection
                  ->reaction_data_id == "reaction" &&
              point_transposition_policy.typed_transposition_selection
                  ->discrete_realisation ==
                nmopt::semantic::v1::TranspositionDiscreteRealisation::
                  fe_q_point_sensor_very_weak,
            "the point-sensor graph omitted its typed transposition selection");
    auto missing_point_diffusion_data = point_sensor_specification;
    component_by_id(missing_point_diffusion_data.requirement_policies,
                    "point_sensor_transposition_policy")
      .typed_transposition_selection->diffusion_data_id.clear();
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(missing_point_diffusion_data),
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "state_equation",
      "transposition_diffusion_data",
      "the point-sensor graph accepted a missing transposition diffusion port");
    auto missing_point_evaluation = point_sensor_specification;
    remove_policy(missing_point_evaluation,
                  "point_sensor_evaluation_policy");
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(missing_point_evaluation),
      nmopt::semantic::v1::DiagnosticCategory::analytical_policy,
      "state_observation",
      "point_sensor_evaluation_policy",
      "the point-sensor graph accepted a missing physical evaluation policy");
    auto missing_point_transposition = point_sensor_specification;
    remove_policy(missing_point_transposition,
                  "point_sensor_transposition_policy");
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(missing_point_transposition),
      nmopt::semantic::v1::DiagnosticCategory::analytical_policy,
      "state_equation",
      "point_sensor_transposition_policy",
      "the point-sensor graph accepted a missing very-weak transposition policy");
    auto missing_point_typed_transposition = point_sensor_specification;
    component_by_id(missing_point_typed_transposition.requirement_policies,
                    "point_sensor_transposition_policy")
      .typed_transposition_selection.reset();
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(missing_point_typed_transposition),
      nmopt::semantic::v1::DiagnosticCategory::analytical_policy,
      "state_equation",
      "transposition_strong_space",
      "the point-sensor graph accepted an untyped transposition selection");
    auto alternate_point_isomorphism = point_sensor_specification;
    component_by_id(alternate_point_isomorphism.requirement_policies,
                    "point_sensor_transposition_policy")
      .typed_transposition_selection->isomorphism_id = "alternate";
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(alternate_point_isomorphism),
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "state_equation",
      "transposition_isomorphism",
      "the point-sensor graph accepted an alternate transposition isomorphism");
    auto display_only_point_transposition = point_sensor_specification;
    for (auto &policy : display_only_point_transposition.requirement_policies)
      if (policy.id == "point_sensor_evaluation_policy" ||
          policy.id == "point_sensor_transposition_policy")
        policy.selected_policy.clear();
    require(validator.validate(display_only_point_transposition).valid(),
            "point-sensor transposition resolution depended on display text");

    auto normal_flux_specification =
      nmopt::semantic::v1::make_normal_flux_scalar_diffusion_reaction_problem();
    require(validator.validate(normal_flux_specification).valid(),
            "the C5.8 normal-flux v1 graph is invalid");
    const auto normal_flux_region =
      component_by_id(normal_flux_specification.regions,
                      "normal_flux_boundary");
    const auto normal_flux_space =
      component_by_id(normal_flux_specification.spaces,
                      "state_observation_space");
    const auto normal_flux_observation =
      component_by_id(normal_flux_specification.observations,
                      "state_observation");
    require(normal_flux_region.kind ==
              nmopt::semantic::v1::RegionKind::boundary &&
              normal_flux_region.boundary_ids == std::vector<unsigned int>{1} &&
              normal_flux_space.topology ==
                nmopt::semantic::v1::SpaceTopology::l2 &&
              normal_flux_space.region_id == "normal_flux_boundary" &&
              normal_flux_observation.kind ==
              nmopt::semantic::v1::ObservationKind::normal_flux,
            "the normal-flux graph did not declare its boundary output pairing");
    const auto custom_normal_flux_specification =
      nmopt::semantic::v1::make_normal_flux_scalar_diffusion_reaction_problem(
        {7}, {0, 7});
    require(validator.validate(custom_normal_flux_specification).valid() &&
              component_by_id(custom_normal_flux_specification.regions,
                              "dirichlet_boundary")
                  .boundary_ids == std::vector<unsigned int>{0, 7},
            "the normal-flux factory did not preserve custom fixed boundary ids");
    auto normal_flux_not_fixed = normal_flux_specification;
    component_by_id(normal_flux_not_fixed.regions, "normal_flux_boundary")
      .boundary_ids = {2};
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(normal_flux_not_fixed),
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "state_observation",
      "normal_flux_fixed_boundary_subset",
      "the normal-flux graph accepted an observed boundary outside the fixed region");
    const auto normal_flux_transposition_policy = component_by_id(
      normal_flux_specification.requirement_policies,
      "normal_flux_transposition_policy");
    require(normal_flux_transposition_policy.typed_transposition_selection
                  .has_value() &&
              normal_flux_transposition_policy.typed_transposition_selection
                  ->diffusion_data_id == "diffusion" &&
              normal_flux_transposition_policy.typed_transposition_selection
                  ->reaction_data_id == "reaction" &&
              normal_flux_transposition_policy.typed_transposition_selection
                  ->discrete_realisation ==
                nmopt::semantic::v1::TranspositionDiscreteRealisation::
                  fe_q_normal_flux_very_weak,
            "the normal-flux graph omitted its typed transposition selection");
    auto missing_normal_flux_orientation = normal_flux_specification;
    remove_policy(missing_normal_flux_orientation,
                  "normal_flux_orientation_policy");
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(missing_normal_flux_orientation),
      nmopt::semantic::v1::DiagnosticCategory::analytical_policy,
      "state_observation",
      "normal_flux_orientation_policy",
      "the normal-flux graph accepted a missing outward-normal policy");
    auto missing_normal_flux_transposition = normal_flux_specification;
    remove_policy(missing_normal_flux_transposition,
                  "normal_flux_transposition_policy");
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(missing_normal_flux_transposition),
      nmopt::semantic::v1::DiagnosticCategory::analytical_policy,
      "state_equation",
      "normal_flux_transposition_policy",
      "the normal-flux graph accepted a missing very-weak transposition policy");
    auto missing_normal_flux_evaluation = normal_flux_specification;
    remove_policy(missing_normal_flux_evaluation,
                  "normal_flux_evaluation_policy");
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(missing_normal_flux_evaluation),
      nmopt::semantic::v1::DiagnosticCategory::analytical_policy,
      "state_observation",
      "normal_flux_evaluation_policy",
      "the normal-flux graph accepted a missing face evaluation policy");
    auto display_only_normal_flux = normal_flux_specification;
    for (auto &policy : display_only_normal_flux.requirement_policies)
      if (policy.id == "normal_flux_orientation_policy" ||
          policy.id == "normal_flux_evaluation_policy" ||
          policy.id == "normal_flux_transposition_policy")
        policy.selected_policy.clear();
    require(validator.validate(display_only_normal_flux).valid(),
            "normal-flux policy resolution depended on display prose");

    auto missing_transposition = l2_dirichlet_specification;
    remove_policy(missing_transposition, "transposition_formulation");
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(missing_transposition),
      nmopt::semantic::v1::DiagnosticCategory::analytical_policy,
      "state_equation",
      "transposition_formulation_policy",
      "the transposition graph accepted a missing continuous formulation policy");
    auto missing_domain_regularity = l2_dirichlet_specification;
    remove_policy(missing_domain_regularity,
                  "transposition_domain_regularity");
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(missing_domain_regularity),
      nmopt::semantic::v1::DiagnosticCategory::analytical_policy,
      "state_equation",
      "transposition_domain_regularity",
      "the transposition graph accepted a missing domain assumption");
    auto missing_trace_subspace = l2_dirichlet_specification;
    remove_policy(missing_trace_subspace, "conforming_trace_subspace");
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(missing_trace_subspace),
      nmopt::semantic::v1::DiagnosticCategory::analytical_policy,
      "control",
      "transposition_conforming_trace_subspace",
      "the transposition graph accepted a missing conforming trace policy");
    auto missing_conormal = l2_dirichlet_specification;
    remove_policy(missing_conormal, "discrete_conormal_policy");
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(missing_conormal),
      nmopt::semantic::v1::DiagnosticCategory::analytical_policy,
      "state_equation",
      "transposition_conormal_policy",
      "the transposition graph accepted a missing conormal policy");
    auto wrong_transposition_state = l2_dirichlet_specification;
    component_by_id(wrong_transposition_state.spaces, "state_space").topology =
      nmopt::semantic::v1::SpaceTopology::h1;
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(wrong_transposition_state),
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "state_equation",
      "transposition_space_topologies",
      "the transposition graph accepted an H1 continuous state declaration");
    const auto l2_transposition_policy = component_by_id(
      l2_dirichlet_specification.requirement_policies,
      "transposition_formulation");
    require(l2_transposition_policy.typed_transposition_selection.has_value() &&
              l2_transposition_policy.typed_transposition_selection
                  ->continuous_parent_space_id == "control_space" &&
              l2_transposition_policy.typed_transposition_selection
                  ->diffusion_data_id.empty() &&
              l2_transposition_policy.typed_transposition_selection
                  ->reaction_data_id.empty() &&
              l2_transposition_policy.typed_transposition_selection
                  ->equivalence_realisation ==
                nmopt::semantic::v1::TranspositionEquivalenceRealisation::
                  conforming_lifting_variational_equivalence,
            "the transposition graph omitted its typed equivalence selection");
    auto missing_l2_typed_transposition = l2_dirichlet_specification;
    component_by_id(missing_l2_typed_transposition.requirement_policies,
                    "transposition_formulation")
      .typed_transposition_selection.reset();
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(missing_l2_typed_transposition),
      nmopt::semantic::v1::DiagnosticCategory::analytical_policy,
      "state_equation",
      "transposition_strong_space",
      "the transposition graph accepted an untyped formulation selection");
    auto alternate_l2_equivalence = l2_dirichlet_specification;
    component_by_id(alternate_l2_equivalence.requirement_policies,
                    "transposition_formulation")
      .typed_transposition_selection->equivalence_realisation =
      nmopt::semantic::v1::TranspositionEquivalenceRealisation::none;
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(alternate_l2_equivalence),
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "state_equation",
      "transposition_equivalence",
      "the transposition graph accepted an alternate lifting equivalence");

    auto hhalf_dirichlet_specification =
      nmopt::semantic::v1::make_hhalf_dirichlet_laplace_control_problem();
    require(validator.validate(hhalf_dirichlet_specification).valid(),
            "the H1/2 Dirichlet-control graph is invalid");
    require(
      component_by_id(hhalf_dirichlet_specification.spaces, "control_space")
            .topology == nmopt::semantic::v1::SpaceTopology::hhalf &&
        component_by_id(hhalf_dirichlet_specification.residual_terms,
                        "laplacian")
            .kind == nmopt::semantic::v1::ResidualTermKind::laplacian &&
        component_by_id(hhalf_dirichlet_specification.losses,
                        "control_regularisation")
            .kind == nmopt::semantic::v1::LossKind::
                       quadratic_hhalf_control_regularisation &&
        component_by_id(hhalf_dirichlet_specification.metrics,
                        "control_hhalf_metric")
            .kind == nmopt::semantic::v1::MetricKind::hhalf &&
        std::none_of(hhalf_dirichlet_specification.data.begin(),
                     hhalf_dirichlet_specification.data.end(),
                     [](const auto &datum) {
                       return datum.role ==
                                nmopt::semantic::v1::DataRole::diffusion ||
                              datum.role ==
                                nmopt::semantic::v1::DataRole::reaction;
                     }),
      "the H1/2 graph omitted its normalized Laplacian or fractional geometry");
    const auto hhalf_metric_policy = component_by_id(
      hhalf_dirichlet_specification.requirement_policies,
      "control_hhalf_metric_realisation");
    require(hhalf_metric_policy.typed_fractional_metric_selection.has_value() &&
              hhalf_metric_policy.typed_fractional_metric_selection
                  ->control_space_id == "control_space" &&
              hhalf_metric_policy.typed_fractional_metric_selection
                  ->volume_operator_id == "volume_mass_plus_stiffness",
            "the H1/2 graph omitted its typed fractional metric selection");

    auto h1_tracking_hhalf_specification = nmopt::semantic::v1::
      make_h1_tracking_hhalf_dirichlet_laplace_control_problem();
    require(validator.validate(h1_tracking_hhalf_specification).valid() &&
              component_by_id(h1_tracking_hhalf_specification.observations,
                              "state_observation")
                  .kind == nmopt::semantic::v1::ObservationKind::
                             h1_state_restriction &&
              component_by_id(h1_tracking_hhalf_specification.spaces,
                              "control_observation_space")
                  .topology == nmopt::semantic::v1::SpaceTopology::l2 &&
              component_by_id(h1_tracking_hhalf_specification.losses,
                              "control_regularisation")
                  .kind == nmopt::semantic::v1::LossKind::
                             quadratic_control_regularisation &&
              component_by_id(h1_tracking_hhalf_specification.metrics,
                              "control_hhalf_metric")
                  .kind == nmopt::semantic::v1::MetricKind::hhalf,
            "Section 5.11.1 option 2 did not separate loss, observation, and metric");

    auto h1_dirichlet_specification =
      nmopt::semantic::v1::make_h1_dirichlet_laplace_control_problem();
    require(validator.validate(h1_dirichlet_specification).valid() &&
              component_by_id(h1_dirichlet_specification.spaces,
                              "control_space")
                  .topology == nmopt::semantic::v1::SpaceTopology::h1 &&
              component_by_id(h1_dirichlet_specification.losses,
                              "control_regularisation")
                  .kind == nmopt::semantic::v1::LossKind::
                             quadratic_h1_control_regularisation &&
              component_by_id(h1_dirichlet_specification.metrics,
                              "control_h1_metric")
                  .kind == nmopt::semantic::v1::MetricKind::h1,
            "the tangential H1 Dirichlet-control graph is invalid");
    const auto h1_metric_policy = component_by_id(
      h1_dirichlet_specification.requirement_policies,
      "control_h1_metric_realisation");
    require(h1_metric_policy.typed_boundary_h1_metric_selection.has_value() &&
              h1_metric_policy.typed_boundary_h1_metric_selection
                  ->boundary_region_id == "control_boundary" &&
              h1_metric_policy.typed_boundary_h1_metric_selection
                  ->tangential_gradient_realisation ==
                nmopt::semantic::v1::BoundaryH1TangentialGradientRealisation::
                  projected_ambient_gradient,
            "the tangential H1 graph omitted its typed metric selection");

    auto missing_hhalf_realisation = hhalf_dirichlet_specification;
    remove_policy(missing_hhalf_realisation,
                  "control_hhalf_metric_realisation");
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(missing_hhalf_realisation),
      nmopt::semantic::v1::DiagnosticCategory::analytical_policy,
      "control_hhalf_metric",
      "hhalf_metric_realisation_policy",
      "the fractional metric accepted a missing discrete realization");
    auto wrong_hhalf_space = hhalf_dirichlet_specification;
    component_by_id(wrong_hhalf_space.spaces, "control_space").topology =
      nmopt::semantic::v1::SpaceTopology::l2;
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(wrong_hhalf_space),
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "control_hhalf_metric",
      "hhalf_metric_search_space",
      "the fractional metric accepted a non-H1/2 control space");
    auto missing_tangential_realisation = h1_dirichlet_specification;
    remove_policy(missing_tangential_realisation,
                  "control_h1_metric_realisation");
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(missing_tangential_realisation),
      nmopt::semantic::v1::DiagnosticCategory::analytical_policy,
      "control_h1_metric",
      "boundary_h1_metric_tangential_policy",
      "the boundary H1 metric accepted a missing tangential realization");
    auto missing_hhalf_typed_selection = hhalf_dirichlet_specification;
    component_by_id(missing_hhalf_typed_selection.requirement_policies,
                    "control_hhalf_metric_realisation")
      .typed_fractional_metric_selection.reset();
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(missing_hhalf_typed_selection),
      nmopt::semantic::v1::DiagnosticCategory::analytical_policy,
      "control_hhalf_metric",
      "hhalf_metric_realisation_selection",
      "the fractional metric accepted an untyped realization");
    auto alternate_hhalf_operator = hhalf_dirichlet_specification;
    component_by_id(alternate_hhalf_operator.requirement_policies,
                    "control_hhalf_metric_realisation")
      .typed_fractional_metric_selection->operator_realisation =
      nmopt::semantic::v1::FractionalTraceOperatorRealisation::unspecified;
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(alternate_hhalf_operator),
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "control_hhalf_metric",
      "hhalf_metric_realisation_selection",
      "the fractional metric accepted an alternate operator realization");
    auto alternate_h1_tangential = h1_dirichlet_specification;
    component_by_id(alternate_h1_tangential.requirement_policies,
                    "control_h1_metric_realisation")
      .typed_boundary_h1_metric_selection->tangential_gradient_realisation =
      nmopt::semantic::v1::BoundaryH1TangentialGradientRealisation::
        unspecified;
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(alternate_h1_tangential),
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "control_h1_metric",
      "boundary_h1_metric_realisation_selection",
      "the boundary H1 metric accepted an alternate tangential realization");

    auto partial_dirichlet_control_specification =
      nmopt::semantic::v1::
        make_partial_dirichlet_control_scalar_diffusion_reaction_problem();
    require(validator.validate(partial_dirichlet_control_specification).valid(),
            "the partial Dirichlet-control lifting v1 graph is invalid");
    require(component_by_id(partial_dirichlet_control_specification.transformations,
                            "dirichlet_control_lifting")
                .fixed_data_id == "fixed_dirichlet_data",
            "the partial Dirichlet-control graph omitted its fixed lifting port");
    const auto partial_partition_policy = component_by_id(
      partial_dirichlet_control_specification.requirement_policies,
      "partial_dirichlet_boundary_partition");
    require(partial_partition_policy.typed_partial_boundary_selection.has_value() &&
              partial_partition_policy.typed_partial_boundary_selection
                  ->fixed_boundary_region_id == "fixed_dirichlet_boundary" &&
              partial_partition_policy.typed_partial_boundary_selection
                  ->interface_realisation ==
                nmopt::semantic::v1::PartialDirichletInterfaceRealisation::
                  fixed_data_precedence,
            "the partial Dirichlet graph omitted its typed boundary partition");
    auto missing_partial_partition = partial_dirichlet_control_specification;
    remove_policy(missing_partial_partition,
                  "partial_dirichlet_boundary_partition");
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(missing_partial_partition),
      nmopt::semantic::v1::DiagnosticCategory::analytical_policy,
      "state",
      "partial_dirichlet_partition_policy",
      "the partial Dirichlet graph accepted a missing typed boundary partition");
    auto alternate_partial_interface = partial_dirichlet_control_specification;
    component_by_id(alternate_partial_interface.requirement_policies,
                    "partial_dirichlet_boundary_partition")
      .typed_partial_boundary_selection->interface_realisation =
      nmopt::semantic::v1::PartialDirichletInterfaceRealisation::unspecified;
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(alternate_partial_interface),
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "state",
      "partial_dirichlet_interface_selection",
      "the partial Dirichlet graph accepted an alternate interface owner");

    auto neumann_convection_specification =
      nmopt::semantic::v1::make_neumann_convection_subdomain_tracking_problem(
        1);
    require(validator.validate(neumann_convection_specification).valid(),
            "the C5.6 Neumann convection composition graph is invalid");
    require(component_by_id(neumann_convection_specification.observations,
                            "state_subdomain_restriction")
                .kind == nmopt::semantic::v1::ObservationKind::volume_restriction,
            "the C5.6 composition did not select the volume observation");
    const auto neumann_desired_state_policy = component_by_id(
      neumann_convection_specification.requirement_policies,
      "desired_state_quadrature_policy");
    require(neumann_desired_state_policy.region_id == "observation_subdomain" &&
              neumann_desired_state_policy.selected_policy.find(
                "volume quadrature") != std::string::npos &&
              neumann_desired_state_policy.selected_policy.find(
                "boundary face") == std::string::npos,
            "the C5.6 desired-state policy did not select volume quadrature");

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
    const auto require_h1_target_selection =
      [](const auto &candidate, const std::string &expected_boundary) {
        const auto policy = component_by_id(
          candidate.requirement_policies, "h1_target_data_membership");
        return policy.status ==
                   nmopt::semantic::v1::RequirementStatus::user_assumed &&
               policy.scope ==
                 nmopt::semantic::v1::RequirementScope::continuous_semantics &&
               policy.typed_h1_target_data_membership_selection.has_value() &&
               policy.typed_h1_target_data_membership_selection->data_id ==
                 "desired_state" &&
               policy.typed_h1_target_data_membership_selection
                   ->observation_space_id == "state_observation_space" &&
               policy.typed_h1_target_data_membership_selection
                   ->fixed_boundary_region_id == expected_boundary &&
               policy.typed_h1_target_data_membership_selection
                   ->regularity_realisation ==
                 nmopt::semantic::v1::H1TargetDataRegularityRealisation::
                   h1_value_and_weak_gradient &&
               policy.typed_h1_target_data_membership_selection
                   ->trace_realisation ==
                 nmopt::semantic::v1::H1TargetDataTraceRealisation::
                   zero_trace_on_fixed_boundary;
      };
    require(
      require_h1_target_selection(h1_state_specification, "dirichlet_boundary"),
      "the H1-state graph omitted its typed target-data membership assumption");
    require(
      require_h1_target_selection(h1_tracking_hhalf_specification,
                                  "control_boundary"),
      "the Section 5.11.1 H1-state graph omitted its target-data membership assumption");

    auto hminus1_metric_specification =
      nmopt::semantic::v1::
        make_hminus1_metric_h1_state_tracking_scalar_diffusion_reaction_problem();
    auto hminus1_l2_companion = nmopt::semantic::v1::
      make_l2_metric_h1_state_tracking_continuous_control_problem();
    require(validator.validate(hminus1_metric_specification).valid(),
            "the H-1 metric v1 graph is invalid");
    require(validator.validate(hminus1_l2_companion).valid() &&
              component_by_id(hminus1_l2_companion.metrics,
                              "control_l2_metric")
                  .kind == nmopt::semantic::v1::MetricKind::l2,
            "the H-1 metric L2 comparison graph is invalid");
    require(require_h1_target_selection(hminus1_l2_companion,
                                        "dirichlet_boundary") &&
              require_h1_target_selection(hminus1_metric_specification,
                                           "dirichlet_boundary"),
            "the H1-state comparison factories did not inherit the target-data assumption");
    auto missing_h1_target_membership = h1_state_specification;
    remove_policy(missing_h1_target_membership,
                  "h1_target_data_membership");
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(missing_h1_target_membership),
      nmopt::semantic::v1::DiagnosticCategory::analytical_policy,
      "desired_state",
      "h1_target_space_membership",
      "the H1-state graph accepted a missing target-data membership assumption");
    auto wrong_h1_target_data = h1_state_specification;
    component_by_id(wrong_h1_target_data.requirement_policies,
                    "h1_target_data_membership")
      .typed_h1_target_data_membership_selection->data_id = "forcing";
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(wrong_h1_target_data),
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "desired_state",
      "h1_target_space_membership",
      "the H1-state graph accepted the wrong target datum");
    auto wrong_h1_target_space = h1_state_specification;
    component_by_id(wrong_h1_target_space.requirement_policies,
                    "h1_target_data_membership")
      .typed_h1_target_data_membership_selection
      ->observation_space_id = "state_space";
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(wrong_h1_target_space),
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "desired_state",
      "h1_target_space_membership",
      "the H1-state graph accepted the wrong target observation space");
    auto missing_h1_target_zero_trace = h1_state_specification;
    component_by_id(missing_h1_target_zero_trace.requirement_policies,
                    "h1_target_data_membership")
      .typed_h1_target_data_membership_selection->trace_realisation =
      nmopt::semantic::v1::H1TargetDataTraceRealisation::unspecified;
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(missing_h1_target_zero_trace),
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "desired_state",
      "h1_target_zero_trace",
      "the H1-state graph accepted a missing zero-trace selection");
    auto wrong_h1_target_boundary = h1_state_specification;
    component_by_id(wrong_h1_target_boundary.requirement_policies,
                    "h1_target_data_membership")
      .typed_h1_target_data_membership_selection->fixed_boundary_region_id =
      "domain";
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(wrong_h1_target_boundary),
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "desired_state",
      "h1_target_zero_trace",
      "the H1-state graph accepted a wrong zero-trace boundary");
    auto wrong_h1_target_scope = h1_state_specification;
    component_by_id(wrong_h1_target_scope.requirement_policies,
                    "h1_target_data_membership")
      .scope = nmopt::semantic::v1::RequirementScope::both;
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(wrong_h1_target_scope),
      nmopt::semantic::v1::DiagnosticCategory::analytical_policy,
      "desired_state",
      "h1_target_space_membership",
      "the H1-state graph accepted a target assumption with the wrong scope");
    auto provided_h1_target_assumption = h1_state_specification;
    component_by_id(provided_h1_target_assumption.requirement_policies,
                    "h1_target_data_membership")
      .status = nmopt::semantic::v1::RequirementStatus::provided;
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(provided_h1_target_assumption),
      nmopt::semantic::v1::DiagnosticCategory::analytical_policy,
      "desired_state",
      "h1_target_space_membership",
      "the H1-state graph accepted a target assumption presented as proven");
    const auto hminus1_control_space =
      component_by_id(hminus1_metric_specification.spaces, "control_space");
    const auto hminus1_metric = component_by_id(
      hminus1_metric_specification.metrics, "control_hminus1_metric");
    require(hminus1_control_space.topology ==
              nmopt::semantic::v1::SpaceTopology::h1 &&
              hminus1_metric.kind ==
                nmopt::semantic::v1::MetricKind::hminus1 &&
              hminus1_metric_specification.formulation.metric_id ==
                hminus1_metric.id,
            "the H-1 metric factory did not select its control space and Riesz map");

    const auto hminus1_realisation_policy = component_by_id(
      hminus1_metric_specification.requirement_policies,
      "hminus1_metric_realisation");
    require(hminus1_realisation_policy.typed_metric_selection.has_value() &&
              hminus1_realisation_policy.typed_metric_selection->metric_id ==
                "control_hminus1_metric" &&
              hminus1_realisation_policy.typed_metric_selection
                  ->fixed_boundary_region_id == "dirichlet_boundary",
            "the H-1 metric factory omitted its typed operator and boundary selection");

    auto missing_hminus1_boundary_policy = hminus1_metric_specification;
    remove_policy(missing_hminus1_boundary_policy,
                  "control_homogeneous_dirichlet_policy");
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(missing_hminus1_boundary_policy),
      nmopt::semantic::v1::DiagnosticCategory::analytical_policy,
      "control_hminus1_metric",
      "hminus1_metric_boundary_policy",
      "v1 semantic validation accepted an H-1 metric without its boundary policy");

    auto missing_hminus1_realisation = hminus1_metric_specification;
    remove_policy(missing_hminus1_realisation, "hminus1_metric_realisation");
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(missing_hminus1_realisation),
      nmopt::semantic::v1::DiagnosticCategory::analytical_policy,
      "control_hminus1_metric",
      "hminus1_metric_realisation",
      "v1 semantic validation accepted an H-1 metric without its typed realization");

    auto alternate_hminus1_operator = hminus1_metric_specification;
    component_by_id(alternate_hminus1_operator.requirement_policies,
                    "hminus1_metric_realisation")
      .typed_metric_selection->operator_realisation =
      nmopt::semantic::v1::Hminus1MetricOperatorRealisation::unspecified;
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(alternate_hminus1_operator),
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "control_hminus1_metric",
      "hminus1_metric_operator",
      "v1 semantic validation accepted an alternate H-1 metric operator");

    auto alternate_hminus1_inverse = hminus1_metric_specification;
    component_by_id(alternate_hminus1_inverse.requirement_policies,
                    "hminus1_metric_realisation")
      .typed_metric_selection->inverse_realisation =
      nmopt::semantic::v1::Hminus1MetricInverseRealisation::unspecified;
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(alternate_hminus1_inverse),
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "control_hminus1_metric",
      "hminus1_metric_inverse",
      "v1 semantic validation accepted an alternate H-1 inverse sequence");

    auto mismatched_hminus1_boundary = hminus1_metric_specification;
    component_by_id(mismatched_hminus1_boundary.requirement_policies,
                    "hminus1_metric_realisation")
      .typed_metric_selection->fixed_boundary_region_id = "domain";
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(mismatched_hminus1_boundary),
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "control_hminus1_metric",
      "hminus1_control_boundary",
      "v1 semantic validation accepted a mismatched H-1 control boundary");

    auto discontinuous_hminus1_search_space = hminus1_metric_specification;
    component_by_id(discontinuous_hminus1_search_space.spaces, "control_space")
      .topology = nmopt::semantic::v1::SpaceTopology::l2;
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(discontinuous_hminus1_search_space),
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "control_hminus1_metric",
      "hminus1_metric_search_space",
      "v1 semantic validation accepted the selected H-1 metric on the wrong search space");

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
    const auto weighted_trace_policy = component_by_id(
      weighted_boundary_specification.requirement_policies,
      "state_boundary_trace_policy");
    require(weighted_trace_policy.typed_trace_selection.has_value() &&
              weighted_trace_policy.typed_trace_selection->weight_data_id ==
                "boundary_weight" &&
              weighted_trace_policy.typed_trace_selection->pairing_id ==
                "state_observation_pairing",
            "the weighted boundary trace omitted its typed map selection");

    auto missing_weight_realisation = weighted_boundary_specification;
    component_by_id(missing_weight_realisation.requirement_policies,
                    "state_boundary_trace_policy")
      .typed_trace_selection.reset();
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(missing_weight_realisation),
      nmopt::semantic::v1::DiagnosticCategory::analytical_policy,
      "weighted_state_boundary_trace",
      "weighted_trace_realisation",
      "v1 semantic validation accepted a weighted trace without its typed realization");

    auto alternate_weight_rule = weighted_boundary_specification;
    component_by_id(alternate_weight_rule.requirement_policies,
                    "state_boundary_trace_policy")
      .typed_trace_selection->weight_realisation =
      nmopt::semantic::v1::TraceWeightRealisation::unspecified;
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(alternate_weight_rule),
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "weighted_state_boundary_trace",
      "weighted_trace_weight_rule",
      "v1 semantic validation accepted an alternate weighted trace rule");

    auto display_only_weight_change = weighted_boundary_specification;
    component_by_id(display_only_weight_change.requirement_policies,
                    "state_boundary_trace_policy")
      .selected_policy = "arbitrary display text";
    require(validator.validate(display_only_weight_change).valid(),
            "weighted trace lowering policy depended on display text");

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

    const auto require_data_diagnostic =
      [&validator, &general_scalar_specification](
        const DataRole role,
        const char *  capability,
        const char *  description,
        const auto &  mutate) {
        auto broken = general_scalar_specification;
        mutate(broken);
        const auto &datum = *std::find_if(
          broken.data.begin(), broken.data.end(), [role](const DataSpec &candidate) {
            return candidate.role == role;
          });
        nmopt::test_support::require_exact_diagnostic(
          validator.validate(broken),
          DiagnosticCategory::structural,
          datum.id,
          capability,
          description);
      };
    require_data_diagnostic(
      DataRole::diffusion,
      "coefficient_data_space",
      "v1 semantic validation accepted diffusion data without a space",
      [](ProblemSpec &broken) {
        component_by_id(broken.data, "diffusion_tensor").space_id.clear();
      });
    require_data_diagnostic(
      DataRole::conservative_transport,
      "coefficient_data_space",
      "v1 semantic validation accepted conservative transport data without a space",
      [](ProblemSpec &broken) {
        component_by_id(broken.data, "conservative_transport").space_id.clear();
      });
    require_data_diagnostic(
      DataRole::advective_transport,
      "coefficient_data_space",
      "v1 semantic validation accepted advective transport data without a space",
      [](ProblemSpec &broken) {
        component_by_id(broken.data, "advective_transport").space_id.clear();
      });
    require_data_diagnostic(
      DataRole::reaction,
      "coefficient_data_space",
      "v1 semantic validation accepted reaction data without a space",
      [](ProblemSpec &broken) {
        component_by_id(broken.data, "reaction").space_id.clear();
      });
    require_data_diagnostic(
      DataRole::robin_coefficient,
      "coefficient_data_space",
      "v1 semantic validation accepted Robin coefficient data without a space",
      [](ProblemSpec &broken) {
        component_by_id(broken.data, "robin_coefficient").space_id.clear();
      });
    require_data_diagnostic(
      DataRole::robin_source,
      "coefficient_data_space",
      "v1 semantic validation accepted Robin source data without a space",
      [](ProblemSpec &broken) {
        component_by_id(broken.data, "robin_source").space_id.clear();
      });
    require_data_diagnostic(
      DataRole::reaction,
      "coefficient_data_space_shape",
      "v1 semantic validation accepted reaction data in a non-data space",
      [](ProblemSpec &broken) {
        component_by_id(broken.data, "reaction").space_id = "control_space";
      });
    require_data_diagnostic(
      DataRole::diffusion,
      "volume_coefficient_region",
      "v1 semantic validation accepted volume coefficient data on a boundary",
      [](ProblemSpec &broken) {
        component_by_id(broken.spaces, "diffusion_data_space").region_id =
          "robin_boundary";
      });
    require_data_diagnostic(
      DataRole::robin_coefficient,
      "robin_coefficient_region",
      "v1 semantic validation accepted Robin coefficient data on the wrong region",
      [](ProblemSpec &broken) {
        component_by_id(broken.spaces, "robin_coefficient_data_space")
          .region_id = "domain";
      });
    require_data_diagnostic(
      DataRole::robin_source,
      "coefficient_data_space_shape",
      "v1 semantic validation accepted Robin source data in state_test_space",
      [](ProblemSpec &broken) {
        component_by_id(broken.data, "robin_source").space_id = "state_test_space";
      });
    require_data_diagnostic(
      DataRole::robin_source,
      "robin_source_region",
      "v1 semantic validation accepted Robin source data on a different region",
      [](ProblemSpec &broken) {
        component_by_id(broken.spaces, "robin_source_data_space").region_id =
          "domain";
      });
    require_data_diagnostic(
      DataRole::robin_source,
      "robin_source_trace_pairing",
      "v1 semantic validation accepted Robin source data with the wrong trace space",
      [](ProblemSpec &broken) {
        component_by_id(broken.spaces, "robin_source_data_space").topology =
          SpaceTopology::h1;
      });

    auto wrong_tensor_shape = general_scalar_specification;
    component_by_id(wrong_tensor_shape.data, "diffusion_tensor").kind =
      nmopt::semantic::v1::DataKind::vector_function;
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(wrong_tensor_shape),
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "diffusion_tensor",
      "coefficient_data_shape",
      "v1 semantic validation accepted a vector-valued tensor diffusion binding");

    auto missing_boundary_selection = general_scalar_specification;
    component_by_id(missing_boundary_selection.requirement_policies,
                    "scalar_boundary_partition")
      .typed_selection.reset();
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(missing_boundary_selection),
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "scalar_boundary_partition",
      "boundary_partition_selection",
      "v1 semantic validation accepted a missing typed boundary selection");

    auto wrong_fixed_region = general_scalar_specification;
    component_by_id(wrong_fixed_region.requirement_policies,
                    "scalar_boundary_partition")
      .typed_selection->fixed_dirichlet_region_id = "robin_boundary";
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(wrong_fixed_region),
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "scalar_boundary_partition",
      "boundary_partition_fixed_region",
      "v1 semantic validation accepted a mismatched typed fixed boundary");

    auto wrong_robin_region = general_scalar_specification;
    component_by_id(wrong_robin_region.requirement_policies,
                    "scalar_boundary_partition")
      .typed_selection->robin_region_id = "dirichlet_boundary";
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(wrong_robin_region),
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "scalar_boundary_partition",
      "boundary_partition_robin_region",
      "v1 semantic validation accepted a mismatched typed Robin boundary");

    auto nonempty_neumann = general_scalar_specification;
    component_by_id(nonempty_neumann.requirement_policies,
                    "scalar_boundary_partition")
      .typed_selection->neumann_region_ids = {"dirichlet_boundary"};
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(nonempty_neumann),
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "scalar_boundary_partition",
      "boundary_partition_neumann_role",
      "v1 semantic validation accepted a Neumann role in the first P5.1 target");

    auto separate_outflow = general_scalar_specification;
    component_by_id(separate_outflow.requirement_policies,
                    "scalar_boundary_partition")
      .typed_selection->transport_outflow_region_id = "dirichlet_boundary";
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(separate_outflow),
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "scalar_boundary_partition",
      "transport_outflow_region",
      "v1 semantic validation accepted a separate transport outflow region");

    auto nonempty_inflow = general_scalar_specification;
    component_by_id(nonempty_inflow.requirement_policies,
                    "scalar_boundary_partition")
      .typed_selection->transport_inflow_region_ids = {"dirichlet_boundary"};
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(nonempty_inflow),
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "scalar_boundary_partition",
      "transport_inflow_region",
      "v1 semantic validation accepted a transport inflow role in the first P5.1 target");

    auto inward_normal = general_scalar_specification;
    component_by_id(inward_normal.requirement_policies,
                    "scalar_boundary_partition")
      .typed_selection->normal_orientation =
      nmopt::semantic::v1::NormalOrientation::inward;
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(inward_normal),
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "scalar_boundary_partition",
      "conormal_flux_orientation",
      "v1 semantic validation accepted an inward conormal orientation");

    auto opposite_conormal = general_scalar_specification;
    component_by_id(opposite_conormal.requirement_policies,
                    "scalar_boundary_partition")
      .typed_selection->conormal_form =
      nmopt::semantic::v1::ConormalForm::transport_minus_diffusion;
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(opposite_conormal),
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "scalar_boundary_partition",
      "conormal_flux_form",
      "v1 semantic validation accepted the opposite conormal form");

    auto unsupported_trace = general_scalar_specification;
    component_by_id(unsupported_trace.requirement_policies,
                    "scalar_boundary_partition")
      .typed_selection->trace_realisation =
      nmopt::semantic::v1::TraceEvaluationRealisation::unspecified;
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(unsupported_trace),
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "scalar_boundary_partition",
      "robin_trace_realisation",
      "v1 semantic validation accepted an unsupported Robin trace realization");

    auto unsupported_quadrature = general_scalar_specification;
    component_by_id(unsupported_quadrature.requirement_policies,
                    "scalar_boundary_partition")
      .typed_selection->face_quadrature_realisation =
      nmopt::semantic::v1::FaceQuadratureRealisation::unspecified;
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(unsupported_quadrature),
      nmopt::semantic::v1::DiagnosticCategory::structural,
      "scalar_boundary_partition",
      "robin_face_quadrature",
      "v1 semantic validation accepted an unsupported face quadrature realization");

    auto wrong_boundary_scope = general_scalar_specification;
    component_by_id(wrong_boundary_scope.requirement_policies,
                    "scalar_boundary_partition")
      .scope = nmopt::semantic::v1::RequirementScope::continuous_semantics;
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(wrong_boundary_scope),
      nmopt::semantic::v1::DiagnosticCategory::analytical_policy,
      "scalar_boundary_partition",
      "boundary_partition_policy_scope",
      "v1 semantic validation accepted the wrong boundary policy scope");

    auto wrong_boundary_status = general_scalar_specification;
    component_by_id(wrong_boundary_status.requirement_policies,
                    "scalar_boundary_partition")
      .status = nmopt::semantic::v1::RequirementStatus::provided;
    nmopt::test_support::require_exact_diagnostic(
      validator.validate(wrong_boundary_status),
      nmopt::semantic::v1::DiagnosticCategory::analytical_policy,
      "scalar_boundary_partition",
      "boundary_partition_policy_status",
      "v1 semantic validation accepted the wrong boundary policy status");

    auto display_only_boundary_change = general_scalar_specification;
    component_by_id(display_only_boundary_change.requirement_policies,
                    "scalar_boundary_partition")
      .selected_policy = "arbitrary display text";
    require(validator.validate(display_only_boundary_change).valid(),
            "typed boundary lowering changed when its display text changed");

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
    require(ReducedFormulationSpec{}.provenance ==
              FormulationProvenance::unspecified,
            "default formulation provenance is not safe");

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
       }},
      {"reduced_dto", "formulation_provenance", [](ProblemSpec &specification) {
         specification.formulation.provenance =
           FormulationProvenance::unspecified;
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

    auto supplied_otd =
      nmopt::semantic::v1::make_scalar_diffusion_reaction_problem(true);
    supplied_otd.formulation.kind = FormulationKind::all_at_once;
    supplied_otd.formulation.provenance = FormulationProvenance::supplied_otd;
    require(validator.validate(supplied_otd).valid(),
            "a typed supplied OTD formulation was rejected structurally");
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
    using namespace nmopt::semantic::v1;
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
    const auto fixed_services =
      nmopt::compiler::v1::service_plan(*fixed_plan.plan);
    require(fixed_services.transformation ==
              nmopt::compiler::v1::ScalarTransformationOperatorKind::fixed_dirichlet_reconstruction &&
              fixed_services.transformation_handler_id ==
                "dealii.scalar.transformation.fixed_dirichlet" &&
              fixed_services.fixed_data_id == "fixed_dirichlet_data",
            "fixed reconstruction service selection lost its transformation factory");

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
    const auto placement_for_role =
      [](const auto &placements, const DataRole role) {
        return std::find_if(
          placements.begin(),
          placements.end(),
          [role](const nmopt::compiler::v1::ScalarDataPlacement &placement) {
            return placement.role == role;
          });
      };
    const auto diffusion_placement = placement_for_role(
      general_plan.plan->data_placements, DataRole::diffusion);
    const auto conservative_placement = placement_for_role(
      general_plan.plan->data_placements, DataRole::conservative_transport);
    const auto advective_placement = placement_for_role(
      general_plan.plan->data_placements, DataRole::advective_transport);
    const auto reaction_placement = placement_for_role(
      general_plan.plan->data_placements, DataRole::reaction);
    const auto robin_coefficient_placement = placement_for_role(
      general_plan.plan->data_placements, DataRole::robin_coefficient);
    const auto robin_source_placement = placement_for_role(
      general_plan.plan->data_placements, DataRole::robin_source);
    const auto &boundary_selection = general_plan.plan->boundary_selection;
    require(general_plan.succeeded() &&
              general_plan.plan->residual_terms.size() == 8 &&
              // The forcing Function is also a residual data port; the six
              // P5.1 coefficient/Robin ports are checked individually below.
              general_plan.plan->data_placements.size() == 7 &&
              boundary_selection.has_value() &&
              boundary_selection->id == "scalar_boundary_partition" &&
              boundary_selection->subject_id == "state" &&
              boundary_selection->fixed_dirichlet_region_id ==
                "dirichlet_boundary" &&
              boundary_selection->robin_region_id == "robin_boundary" &&
              boundary_selection->neumann_region_ids.empty() &&
              boundary_selection->transport_inflow_region_ids.empty() &&
              boundary_selection->transport_outflow_region_id ==
                "robin_boundary" &&
              boundary_selection->conormal_form ==
                nmopt::semantic::v1::ConormalForm::diffusion_minus_transport &&
              boundary_selection->normal_orientation ==
                nmopt::semantic::v1::NormalOrientation::outward &&
              boundary_selection->trace_realisation ==
                nmopt::semantic::v1::TraceEvaluationRealisation::fe_q_state_trace &&
              boundary_selection->face_quadrature_realisation ==
                nmopt::semantic::v1::FaceQuadratureRealisation::qgauss_face &&
              general_plan.plan->robin_boundary_ids ==
                std::set<unsigned int>{1} &&
              general_plan.plan->provenance.size() == 13 &&
              diffusion_placement != general_plan.plan->data_placements.end() &&
              diffusion_placement->semantic_id == "diffusion_tensor" &&
              diffusion_placement->space_id == "diffusion_data_space" &&
              diffusion_placement->region_id == "domain" &&
              diffusion_placement->kind == DataKind::tensor_function &&
              diffusion_placement->evaluation ==
                nmopt::compiler::v1::ScalarDataEvaluationKind::volume_quadrature &&
              conservative_placement != general_plan.plan->data_placements.end() &&
              conservative_placement->space_id ==
                "conservative_transport_data_space" &&
              conservative_placement->kind == DataKind::vector_function &&
              advective_placement != general_plan.plan->data_placements.end() &&
              advective_placement->space_id == "advective_transport_data_space" &&
              reaction_placement != general_plan.plan->data_placements.end() &&
              reaction_placement->space_id == "reaction_data_space" &&
              reaction_placement->kind == DataKind::function &&
              robin_coefficient_placement !=
                general_plan.plan->data_placements.end() &&
              robin_coefficient_placement->space_id ==
                "robin_coefficient_data_space" &&
              robin_coefficient_placement->region_id == "robin_boundary" &&
              robin_coefficient_placement->evaluation ==
                nmopt::compiler::v1::ScalarDataEvaluationKind::boundary_face_quadrature &&
              robin_source_placement != general_plan.plan->data_placements.end() &&
              robin_source_placement->space_id == "robin_source_data_space" &&
              robin_source_placement->region_id == "robin_boundary" &&
              robin_source_placement->evaluation ==
                nmopt::compiler::v1::ScalarDataEvaluationKind::boundary_face_quadrature &&
              std::any_of(
                general_plan.plan->residual_terms.begin(),
                general_plan.plan->residual_terms.end(),
                [](const nmopt::compiler::v1::ScalarResidualContribution &term) {
                  return term.operator_kind ==
                         nmopt::compiler::v1::ScalarResidualOperatorKind::conservative_transport;
                }),
            "P5.1 general scalar plan omitted a term or Robin boundary contribution");

    const auto general_residual_assembly =
      nmopt::compiler::v1::residual_assembly_plan(*general_plan.plan);
    const auto canonical_residual_assembly =
      nmopt::compiler::v1::residual_assembly_plan(*planned.plan);
    auto incomplete_general_residual_assembly = general_residual_assembly;
    incomplete_general_residual_assembly.residual_terms.pop_back();
    require(general_residual_assembly.has(
              nmopt::compiler::v1::ScalarResidualOperatorKind::tensor_diffusion) &&
              general_residual_assembly.has(
                nmopt::compiler::v1::ScalarResidualOperatorKind::robin_bilinear) &&
              general_residual_assembly.has(
                nmopt::compiler::v1::ScalarResidualOperatorKind::robin_source) &&
              !general_residual_assembly.has(
                nmopt::compiler::v1::ScalarResidualOperatorKind::diffusion_reaction) &&
              canonical_residual_assembly.has(
                nmopt::compiler::v1::ScalarResidualOperatorKind::diffusion_reaction) &&
              !canonical_residual_assembly.has(
                nmopt::compiler::v1::ScalarResidualOperatorKind::tensor_diffusion) &&
              general_residual_assembly.placement("diffusion_tensor") != nullptr &&
              general_residual_assembly.placement("robin_source") != nullptr &&
              general_residual_assembly.robin_boundary_ids ==
                std::set<unsigned int>{1} &&
              general_residual_assembly.registration().has_value() &&
              *general_residual_assembly.registration() ==
                nmopt::compiler::v1::ScalarResidualAssemblyPlan::Registration::
                  general_tensor_transport_robin &&
              canonical_residual_assembly.registration().has_value() &&
              *canonical_residual_assembly.registration() ==
                nmopt::compiler::v1::ScalarResidualAssemblyPlan::Registration::
                  diffusion_reaction &&
              !incomplete_general_residual_assembly.registration().has_value(),
            "scalar residual assembly did not preserve closed typed registrations");

    const auto canonical_services =
      nmopt::compiler::v1::service_plan(*planned.plan);
    const auto general_services =
      nmopt::compiler::v1::service_plan(*general_plan.plan);
    require(canonical_services.has_observation(
              nmopt::compiler::v1::ScalarObservationOperatorKind::volume_restriction) &&
              canonical_services.has_loss(
                nmopt::compiler::v1::ScalarLossOperatorKind::quadratic_tracking) &&
              canonical_services.has_loss(
                nmopt::compiler::v1::ScalarLossOperatorKind::quadratic_control_regularisation) &&
              canonical_services.metric ==
                nmopt::compiler::v1::ScalarMetricOperatorKind::cellwise_l2 &&
              canonical_services.constraint ==
                nmopt::compiler::v1::ScalarConstraintOperatorKind::cellwise_box &&
              canonical_services.metric_handler_id ==
                "dealii.scalar.metric.cellwise_l2" &&
              canonical_services.constraint_handler_id ==
                "dealii.scalar.constraint.cellwise_box" &&
              general_services.observations.size() == 2 &&
              general_services.losses.size() == 2,
            "scalar service plan did not preserve objective and factory selections");

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

    const auto normal_flux_specification =
      nmopt::semantic::v1::make_normal_flux_scalar_diffusion_reaction_problem();
    const auto normal_flux_resolution = resolver.resolve(normal_flux_specification);
    require(normal_flux_resolution.succeeded(),
            "normal-flux lowering-plan setup did not resolve");
    const auto normal_flux_plan = planner.plan(*normal_flux_resolution.problem);
    require(
      normal_flux_plan.succeeded() &&
        normal_flux_plan.plan->normal_flux_boundary_ids ==
          std::set<unsigned int>{1} &&
        normal_flux_plan.plan->normal_flux_evaluation_policy.find(
          "outward normal derivative") != std::string::npos &&
        std::find(normal_flux_plan.plan->provenance.begin(),
                  normal_flux_plan.plan->provenance.end(),
                  "state_observation <- "
                  "dealii.scalar.observation.normal_flux") !=
          normal_flux_plan.plan->provenance.end(),
      "C5.8 normal-flux observation did not contribute its scalar lowering handler");
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
