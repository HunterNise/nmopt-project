#pragma once

// This compiler-planning scenario is included into the existing neutral
// semantic test translation unit to improve responsibility locality without
// multiplying test targets or compilation cost.

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
    auto natural_source_plan_specification = general_specification;
    component_by_id(natural_source_plan_specification.data, "robin_source").role =
      DataRole::natural_boundary_source;
    component_by_id(natural_source_plan_specification.residual_terms, "robin_source")
      .kind = ResidualTermKind::natural_boundary_source;
    const auto natural_source_resolution =
      resolver.resolve(natural_source_plan_specification);
    require(natural_source_resolution.succeeded(),
            "natural-boundary source lowering-plan setup did not resolve");
    const auto natural_source_plan =
      planner.plan(*natural_source_resolution.problem);
    require(natural_source_plan.succeeded(),
            "natural-boundary source lowering plan did not succeed");
    const auto natural_source_residual_assembly =
      nmopt::compiler::v1::residual_assembly_plan(*natural_source_plan.plan);
    const auto natural_source_placement =
      std::find_if(natural_source_plan.plan->data_placements.begin(),
                   natural_source_plan.plan->data_placements.end(),
                   [](const nmopt::compiler::v1::ScalarDataPlacement &placement) {
                     return placement.role == DataRole::natural_boundary_source;
                   });
    require(!natural_source_residual_assembly.registration().has_value() &&
              natural_source_placement !=
                natural_source_plan.plan->data_placements.end() &&
              natural_source_placement->region_id == "robin_boundary" &&
              natural_source_placement->evaluation ==
                nmopt::compiler::v1::ScalarDataEvaluationKind::boundary_face_quadrature &&
              std::any_of(natural_source_plan.plan->residual_terms.begin(),
                          natural_source_plan.plan->residual_terms.end(),
                          [](const nmopt::compiler::v1::ScalarResidualContribution &term) {
                            return term.operator_kind ==
                                   nmopt::compiler::v1::ScalarResidualOperatorKind::natural_boundary_source;
                          }),
            "the scalar lowering plan did not preserve the natural-boundary source contribution");
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
