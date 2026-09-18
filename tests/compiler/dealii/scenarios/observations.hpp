#pragma once

// These scenario implementations are deliberately included into one heavy
// deal.II/compiler translation unit to improve source navigation without
// multiplying compilation cost.

  template <int dim>
  void
  run_volume_observation_assembly_contract_test()
  {
    dealii::Triangulation<dim> triangulation;
    dealii::GridGenerator::hyper_cube(triangulation);
    for (const auto &cell : triangulation.active_cell_iterators())
      {
        cell->set_material_id(7);
        for (unsigned int face = 0; face < cell->n_faces(); ++face)
          if (cell->face(face)->at_boundary())
            cell->face(face)->set_boundary_id(
              cell->face(face)->center()[0] < 1.0e-12 ? 3 : 4);
      }

    dealii::FE_Q<dim> state_fe(1);
    dealii::DoFHandler<dim> state_dof_handler(triangulation);
    state_dof_handler.distribute_dofs(state_fe);
    const dealii::Functions::ConstantFunction<dim> fixed_state(0.7);
    const FirstCoordinateFunction<dim> desired_state;
    std::map<dealii::types::global_dof_index, double> fixed_values;
    dealii::VectorTools::interpolate_boundary_values(state_dof_handler,
                                                      3,
                                                      fixed_state,
                                                      fixed_values);
    std::vector<bool> constrained_state_dofs(state_dof_handler.n_dofs(), false);
    dealii::Vector<double> fixed_state_values(state_dof_handler.n_dofs());
    for (const auto &entry : fixed_values)
      {
        constrained_state_dofs.at(entry.first) = true;
        fixed_state_values[entry.first] = entry.second;
      }

    const compiler::v1::detail::VolumeObservationAssembly<dim> observation(
      state_dof_handler,
      state_fe,
      constrained_state_dofs,
      fixed_state_values,
      {7},
      desired_state,
      3);
    const compiler::v1::detail::VolumeObservationAssembly<dim>
      interpolated_observation(
        state_dof_handler,
        state_fe,
        constrained_state_dofs,
        fixed_state_values,
        {7},
        desired_state,
        2,
        compiler::v1::VolumeObservationTargetRealisation::
          state_fe_interpolation);
    contract::require(
      observation.quadrature_order() == 3 &&
        observation.target_realisation() ==
          compiler::v1::VolumeObservationTargetRealisation::
            analytic_quadrature &&
        interpolated_observation.quadrature_order() == 2 &&
        interpolated_observation.target_realisation() ==
          compiler::v1::VolumeObservationTargetRealisation::
            state_fe_interpolation,
      "Volume-observation assembly did not retain its selected policy");
    dealii::Vector<double> direction(state_dof_handler.n_dofs());
    dealii::Vector<double> state = fixed_state_values;
    for (dealii::types::global_dof_index index = 0;
         index < state_dof_handler.n_dofs();
         ++index)
      if (!constrained_state_dofs.at(index))
        {
          direction[index] = index % 2 == 0 ? 0.2 : -0.15;
          state[index] = 0.1 * static_cast<double>(index + 1);
        }

    const dealii::QGauss<dim> reference_quadrature(3);
    dealii::FEValues<dim> reference_values(
      state_fe,
      reference_quadrature,
      dealii::update_values | dealii::update_quadrature_points |
        dealii::update_JxW_values);
    struct TrackingReference
    {
      double objective = 0.0;
      double derivative = 0.0;
    };
    const auto reference_tracking =
      [&](const dealii::Vector<double> &values,
          const dealii::Vector<double> &tangent) {
        TrackingReference result;
        std::vector<double> state_at_quadrature(reference_quadrature.size());
        std::vector<double> tangent_at_quadrature(reference_quadrature.size());
        for (const auto &cell : state_dof_handler.active_cell_iterators())
          {
            reference_values.reinit(cell);
            reference_values.get_function_values(values, state_at_quadrature);
            reference_values.get_function_values(tangent,
                                                 tangent_at_quadrature);
            for (unsigned int q = 0; q < reference_quadrature.size(); ++q)
              {
                const double error =
                  state_at_quadrature[q] -
                  desired_state.value(reference_values.quadrature_point(q));
                result.objective +=
                  0.5 * error * error * reference_values.JxW(q);
                result.derivative +=
                  error * tangent_at_quadrature[q] * reference_values.JxW(q);
              }
          }
        return result;
      };

    const auto fixed_reference =
      reference_tracking(fixed_state_values, direction);
    const auto state_reference = reference_tracking(state, direction);
    dealii::Vector<double> shifted_state = state;
    shifted_state.add(1.0, direction);
    const auto shifted_reference = reference_tracking(shifted_state, direction);
    dealii::Vector<double> mass_direction(state_dof_handler.n_dofs());
    observation.state_tracking_matrix().vmult(mass_direction, direction);
    require_close(direction * mass_direction,
                  shifted_reference.derivative - state_reference.derivative,
                  1.0e-14,
                  "Volume-observation tracking matrix");
    require_close(0.5 * observation.desired_state_norm(),
                  fixed_reference.objective,
                  1.0e-14,
                  "Volume-observation affine target norm");
    require_close(-(observation.desired_state_load() * direction),
                  fixed_reference.derivative,
                  1.0e-14,
                  "Volume-observation affine target load");

    dealii::Vector<double> mass_state(state_dof_handler.n_dofs());
    observation.state_tracking_matrix().vmult(mass_state, state);
    const double assembled_objective =
      0.5 * (state * mass_state) -
      observation.desired_state_load() * state +
      0.5 * observation.desired_state_norm();
    require_close(assembled_objective,
                  state_reference.objective,
                  1.0e-14,
                  "Volume-observation objective reconstruction");
    mass_state.add(-1.0, observation.desired_state_load());
    require_close(mass_state * direction,
                  state_reference.derivative,
                  1.0e-14,
                  "Volume-observation state derivative");

    dealii::Vector<double> interpolated_mass_direction(
      state_dof_handler.n_dofs());
    interpolated_observation.state_tracking_matrix().vmult(
      interpolated_mass_direction,
      direction);
    interpolated_mass_direction.add(-1.0, mass_direction);
    dealii::Vector<double> interpolated_load_difference =
      interpolated_observation.desired_state_load();
    interpolated_load_difference.add(-1.0,
                                     observation.desired_state_load());
    require_close(interpolated_mass_direction.l2_norm(),
                  0.0,
                  1.0e-14,
                  "Volume-observation polynomial-exact mass matrix");
    require_close(interpolated_load_difference.l2_norm(),
                  0.0,
                  1.0e-14,
                  "Volume-observation polynomial-exact target load");
    require_close(interpolated_observation.desired_state_norm(),
                  observation.desired_state_norm(),
                  1.0e-14,
                  "Volume-observation polynomial-exact target norm");

    const EnergyPolynomial<dim> nonrepresentable_target;
    const compiler::v1::detail::VolumeObservationAssembly<dim>
      analytic_nonrepresentable_observation(
        state_dof_handler,
        state_fe,
        constrained_state_dofs,
        fixed_state_values,
        {7},
        nonrepresentable_target,
        3,
        compiler::v1::VolumeObservationTargetRealisation::
          analytic_quadrature);
    const compiler::v1::detail::VolumeObservationAssembly<dim>
      interpolated_nonrepresentable_observation(
        state_dof_handler,
        state_fe,
        constrained_state_dofs,
        fixed_state_values,
        {7},
        nonrepresentable_target,
        3,
        compiler::v1::VolumeObservationTargetRealisation::
          state_fe_interpolation);
    contract::require(
      std::abs(analytic_nonrepresentable_observation.desired_state_norm() -
               interpolated_nonrepresentable_observation.desired_state_norm()) >
        1.0e-4,
      "Volume-observation target realization did not change a "
      "nonrepresentable target");
  }

  template <int dim>
  void
  run_subdomain_observation_contract_test()
  {
    dealii::Triangulation<dim> triangulation;
    dealii::GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(2);
    bool first_cell = true;
    for (auto cell = triangulation.begin_active();
         cell != triangulation.end();
         ++cell)
      {
        cell->set_material_id(first_cell ? 1 : 0);
        first_cell = false;
      }

    const dealii::Functions::ConstantFunction<dim> forcing(1.0);
    const dealii::Functions::ConstantFunction<dim> desired_state(0.0);
    const auto material_one_specification =
      semantic::v1::make_subdomain_tracking_scalar_diffusion_reaction_problem(1);
    const auto material_zero_specification =
      semantic::v1::make_subdomain_tracking_scalar_diffusion_reaction_problem(0);
    compiler::v1::DealiiDiscretisationPolicy policy;
    policy.state_degree = 1;
    const compiler::v1::DealiiCompiler compiler;
    contract::require(
      compiler.validate(material_one_specification, policy).valid() &&
        compiler.validate(material_zero_specification, policy).valid(),
      "material-subdomain v1 graphs did not validate for deal.II");

    const compiler::v1::DealiiDataBindings<dim> bindings{
      forcing,
      desired_state,
      1.0,
      0.5,
      0.1,
      test_binding_provenance("subdomain_observation")};
    const auto material_one = compiler.compile(material_one_specification,
                                               triangulation,
                                               bindings,
                                               policy);
    const auto material_zero = compiler.compile(material_zero_specification,
                                                triangulation,
                                                bindings,
                                                policy);
    contract::require(material_one.succeeded() && material_zero.succeeded(),
                      "material-subdomain v1 compilation failed");
    require_compiled_hessian_evidence(*material_one.problem,
                                      "subdomain observation");

    const auto &one_model = material_one.problem->executable_model();
    const auto &zero_model = material_zero.problem->executable_model();
    dealii::Vector<double> one_control_values(
      one_model.variable_layout()->dimension(1));
    dealii::Vector<double> zero_control_values(
      zero_model.variable_layout()->dimension(1));
    const Primal one_control(one_model.variable_layout()->single_block(1,
                                                                       "control"),
                             {std::move(one_control_values)});
    const Primal zero_control(zero_model.variable_layout()->single_block(1,
                                                                         "control"),
                              {std::move(zero_control_values)});
    const auto one_evaluation =
      material_one.problem->make_reduced_dto().evaluate(one_control);
    const auto zero_evaluation =
      material_zero.problem->make_reduced_dto().evaluate(zero_control);

    // Changing the observation mask must not change the state solve or its
    // residual operators. Only the tracking mass/load are selected below.
    require_primal_close(one_evaluation.state,
                         zero_evaluation.state,
                         1e-12,
                         "material observation changed the state solution");
    dealii::Vector<double> state_tangent(
      one_model.variable_layout()->dimension(0));
    dealii::Vector<double> control_tangent(
      one_model.variable_layout()->dimension(1));
    for (dealii::types::global_dof_index index = 0;
         index < state_tangent.size();
         ++index)
      state_tangent[index] = 0.02 * static_cast<double>(index + 1);
    for (dealii::types::global_dof_index index = 0;
         index < control_tangent.size();
         ++index)
      control_tangent[index] = -0.03 * static_cast<double>(index + 1);
    const Primal tangent(one_model.variable_layout(),
                         {std::move(state_tangent), std::move(control_tangent)});
    require_covector_close(
      one_model.residual_jvp(one_evaluation.full_point, tangent),
      zero_model.residual_jvp(zero_evaluation.full_point, tangent),
      1e-12,
      "material observation changed the residual operators");

    const Covector one_objective_derivative =
      one_model.objective_derivative(one_evaluation.full_point);
    const Covector zero_objective_derivative =
      zero_model.objective_derivative(zero_evaluation.full_point);
    dealii::Vector<double> state_rhs_difference =
      one_objective_derivative.block(0);
    state_rhs_difference.add(-1.0, zero_objective_derivative.block(0));
    contract::require(
      std::abs(one_evaluation.objective_value -
               zero_evaluation.objective_value) > 1e-6 &&
        state_rhs_difference.l2_norm() > 1e-6,
      "material observation did not change the tracking objective and adjoint RHS");

    const auto &one_manifest = material_one.problem->manifest();
    require_constraint_realisation(one_manifest, "none", "subdomain-tracking");
    require_constraint_realisation(material_zero.problem->manifest(),
                                   "none",
                                   "subdomain-tracking alternate region");
    contract::require(
      one_manifest.compatibility.observation_realisation ==
        "material-id volume restriction: 1" &&
        one_manifest.compatibility.data_rule.find("analytic desired-state Function") !=
          std::string::npos &&
        one_manifest.compatibility.lowering_handler_records.size() == 8 &&
        std::find(one_manifest.compatibility.lowering_handler_records.begin(),
                  one_manifest.compatibility.lowering_handler_records.end(),
                  "diffusion_reaction <- "
                  "dealii.scalar.residual.diffusion_reaction") !=
          one_manifest.compatibility.lowering_handler_records.end(),
      "v1 subdomain observation manifest omitted its restriction or data rule");

    // The same residual and metric are recombined with a fixed-data
    // transformation while only the tracking observation changes.
    const auto fixed_specification =
      semantic::v1::make_fixed_dirichlet_scalar_diffusion_reaction_problem();
    auto fixed_subdomain_specification = fixed_specification;
    fixed_subdomain_specification.id =
      "fixed_dirichlet_scalar_diffusion_reaction_subdomain_tracking";
    fixed_subdomain_specification.regions.push_back(
      {"observation_subdomain", "Material subdomain observation region",
       semantic::v1::RegionKind::volume, false, {}, {1}, {}});
    component_by_id(fixed_subdomain_specification.spaces,
                    "state_observation_space")
      .region_id = "observation_subdomain";
    component_by_id(fixed_subdomain_specification.observations,
                    "state_observation")
      .region_id = "observation_subdomain";
    component_by_id(fixed_subdomain_specification.requirement_policies,
                    "desired_state_quadrature_policy")
      .region_id = "observation_subdomain";
    contract::require(compiler.validate(fixed_specification, policy).valid() &&
                        compiler.validate(fixed_subdomain_specification, policy)
                          .valid(),
                      "fixed reconstruction/subdomain recombination did not validate");

    const dealii::Functions::ConstantFunction<dim> fixed_forcing(0.5);
    const dealii::Functions::ConstantFunction<dim> fixed_data(1.0);
    auto fixed_bindings = compiler::v1::DealiiDataBindings<dim>{
      fixed_forcing,
      desired_state,
      1.0,
      0.5,
      0.1,
      test_binding_provenance("fixed_subdomain", true)};
    fixed_bindings.fixed_dirichlet_data = std::cref(fixed_data);
    const auto fixed_full = compiler.compile(fixed_specification,
                                              triangulation,
                                              fixed_bindings,
                                              policy);
    const auto fixed_subdomain = compiler.compile(fixed_subdomain_specification,
                                                  triangulation,
                                                  fixed_bindings,
                                                  policy);
    contract::require(fixed_full.succeeded() && fixed_subdomain.succeeded(),
                      "fixed reconstruction/subdomain recombination did not compile");

    const auto &fixed_full_model = fixed_full.problem->executable_model();
    const auto &fixed_subdomain_model =
      fixed_subdomain.problem->executable_model();
    dealii::Vector<double> fixed_full_control_values(
      fixed_full_model.variable_layout()->dimension(1));
    dealii::Vector<double> fixed_subdomain_control_values(
      fixed_subdomain_model.variable_layout()->dimension(1));
    const Primal fixed_full_control(
      fixed_full_model.variable_layout()->single_block(1, "control"),
      {std::move(fixed_full_control_values)});
    const Primal fixed_subdomain_control(
      fixed_subdomain_model.variable_layout()->single_block(1, "control"),
      {std::move(fixed_subdomain_control_values)});
    const auto fixed_full_evaluation =
      fixed_full.problem->make_reduced_dto().evaluate(fixed_full_control);
    const auto fixed_subdomain_evaluation =
      fixed_subdomain.problem->make_reduced_dto().evaluate(
        fixed_subdomain_control);
    require_primal_close(fixed_full_evaluation.state,
                         fixed_subdomain_evaluation.state,
                         1e-12,
                         "fixed observation recombination changed the state solve");
    require_close(fixed_full_evaluation.objective_value,
                  0.5,
                  1e-12,
                  "fixed full-volume observation objective accounting");
    contract::require(
      fixed_subdomain_evaluation.objective_value <
        fixed_full_evaluation.objective_value - 1e-6,
      "fixed subdomain observation did not change only the tracking objective");
    contract::require(
      fixed_full.problem->manifest().compatibility.lowering_handler_records ==
        fixed_subdomain.problem->manifest().compatibility.lowering_handler_records &&
        fixed_full.problem->manifest().resolved_decision.metric_record.realisation_id ==
          fixed_subdomain.problem->manifest().resolved_decision.metric_record.realisation_id &&
        fixed_subdomain.problem->manifest().compatibility.observation_realisation.find(
          "material-id volume restriction: 1") != std::string::npos,
      "fixed observation recombination did not preserve unchanged service records");
  }

  template <int dim>
  void
  run_h1_state_observation_contract_test()
  {
    static_assert(dim == 2, "The energy-observation oracle is two-dimensional");
    dealii::Triangulation<dim> triangulation;
    dealii::GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(1);

    constexpr double reaction = 0.5;
    const EnergyPolynomialForcing<dim> forcing(reaction);
    const EnergyPolynomial<dim>        desired_state(0.5);
    const auto specification =
      semantic::v1::make_h1_state_tracking_scalar_diffusion_reaction_problem();
    compiler::v1::DealiiDiscretisationPolicy policy;
    policy.state_degree = 2;
    const compiler::v1::DealiiCompiler compiler;
    contract::require(compiler.validate(specification, policy).valid(),
                      "H1-state tracking v1 graph did not validate for deal.II");

    auto subdomain_specification = specification;
    subdomain_specification.regions.push_back(
      {"energy_subdomain", "Unsupported H1 observation subdomain",
       semantic::v1::RegionKind::volume, false, {}, {0}, {}});
    component_by_id(subdomain_specification.spaces,
                    "state_observation_space")
      .region_id = "energy_subdomain";
    component_by_id(subdomain_specification.observations,
                    "state_observation")
      .region_id = "energy_subdomain";
    component_by_id(subdomain_specification.requirement_policies,
                    "desired_state_quadrature_policy")
      .region_id = "energy_subdomain";
    test_support::require_exact_diagnostic(
      compiler.validate(subdomain_specification, policy),
      semantic::v1::DiagnosticCategory::lowerability,
      specification.id,
      "h1_state_observation_full_domain",
      "H1-state observation compiler accepted an unregistered subdomain target");

    const compiler::v1::DealiiDataBindings<dim> bindings{
      forcing,
      desired_state,
      1.0,
      reaction,
      0.2,
      test_binding_provenance("h1_state_observation")};
    const auto compilation = compiler.compile(specification,
                                              triangulation,
                                              bindings,
                                              policy);
    contract::require(compilation.succeeded(),
                      "H1-state observation v1 compilation failed");
    require_compiled_hessian_evidence(*compilation.problem,
                                      "H1 state observation");
    const auto &model = compilation.problem->executable_model();
    const auto reduced = compilation.problem->make_reduced_dto();
    dealii::Vector<double> control_values(model.variable_layout()->dimension(1));
    const Primal control(model.variable_layout()->single_block(1, "control"),
                         {std::move(control_values)});
    const auto evaluation = reduced.evaluate(control);
    require_close(model.residual(evaluation.full_point).block(0).l2_norm(),
                  0.0,
                  1e-10,
                  "H1-state observation manufactured residual");

    // For y=product_i x_i(1-x_i) and z_d=y/2 on the unit square,
    // ||y||_H1^2=7/300. The state loss is therefore 7/2400. This oracle
    // separately proves that the stiffness contribution is present.
    require_close(evaluation.objective_value,
                  7.0 / 2400.0,
                  2e-10,
                  "H1-state observation energy value");

    dealii::Vector<double> state_tangent = evaluation.full_point.block(0);
    dealii::Vector<double> control_tangent(model.variable_layout()->dimension(1));
    const Primal tangent(model.variable_layout(),
                         {std::move(state_tangent), std::move(control_tangent)});
    const double derivative_pairing =
      contract::pair(model.objective_derivative(evaluation.full_point), tangent);
    require_close(derivative_pairing,
                  7.0 / 600.0,
                  5e-10,
                  "H1-state observation VJP pairing");
    constexpr double derivative_step = 1e-6;
    const double centered_difference =
      (model.objective(
         shifted(evaluation.full_point, tangent, derivative_step)) -
       model.objective(
         shifted(evaluation.full_point, tangent, -derivative_step))) /
      (2.0 * derivative_step);
    require_close(centered_difference,
                  derivative_pairing,
                  2e-10,
                  "H1-state observation JVP/VJP derivative");

    dealii::Vector<double> direction_values(control.layout()->dimension(0));
    for (dealii::types::global_dof_index index = 0;
         index < direction_values.size();
         ++index)
      direction_values[index] =
        (index % 2 == 0 ? 0.04 : -0.03) * static_cast<double>(index + 1);
    const Primal direction(control.layout(), {std::move(direction_values)});
    const double reduced_derivative =
      contract::pair(evaluation.reduced_derivative, direction);
    const auto remainder = [&](const double step) {
      return std::abs(reduced.evaluate(shifted(control, direction, step))
                        .objective_value -
                      evaluation.objective_value - step * reduced_derivative);
    };
    const double coarse_remainder = remainder(1e-3);
    const double fine_remainder = remainder(5e-4);
    contract::require(coarse_remainder > 1e-13 &&
                        fine_remainder <= 0.26 * coarse_remainder + 1e-13,
                      "H1-state observation reduced Taylor remainder is not quadratic");

    const auto &manifest = compilation.problem->manifest();
    require_constraint_realisation(manifest, "none", "H1-state observation");
    contract::require(
      manifest.resolved_decision.h1_target_data_membership_selection.has_value() &&
        manifest.resolved_decision.h1_target_data_membership_selection->data_id ==
          "desired_state" &&
        manifest.resolved_decision.h1_target_data_membership_selection->observation_space_id ==
          "state_observation_space" &&
        manifest.resolved_decision.h1_target_data_membership_selection
            ->fixed_boundary_region_id == "dirichlet_boundary" &&
        std::any_of(manifest.compatibility.declared_assumptions.begin(),
                    manifest.compatibility.declared_assumptions.end(),
                    [](const std::string &assumption) {
                      return assumption.find(
                               "h1_target_data_membership: status=user_assumed") ==
                             0;
                    }) &&
      manifest.compatibility.compiler_id == "nmopt.compiler.v1.dealii.h1_state_tracking" &&
        manifest.compatibility.observation_realisation.find("H1_0") != std::string::npos &&
        manifest.compatibility.observation_realisation.find("mass-plus-stiffness") !=
          std::string::npos &&
        manifest.compatibility.data_rule.find("value and gradient") != std::string::npos &&
        manifest.resolved_decision.metric_record.realisation_id == "l2_cellwise" &&
        std::find(manifest.compatibility.lowering_handler_records.begin(),
                  manifest.compatibility.lowering_handler_records.end(),
                  "state_observation <- "
                  "dealii.scalar.observation.h1_state_restriction") !=
          manifest.compatibility.lowering_handler_records.end(),
      "H1-state observation manifest omitted its observation, data, or metric provenance");
  }

  template <int dim>
  void
  run_point_sensor_contract_test()
  {
    static_assert(dim == 2, "The point-sensor oracle is two-dimensional");
    dealii::Triangulation<dim> triangulation;
    dealii::GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(2);

    const dealii::Functions::ConstantFunction<dim> forcing(0.0);
    const FirstCoordinateFunction<dim>             desired_state;
    const std::vector<std::vector<double>> sensor_coordinates{
      {0.23, 0.37}, {0.71, 0.62}};
    const auto specification =
      semantic::v1::make_point_sensor_scalar_diffusion_reaction_problem(
        sensor_coordinates);
    compiler::v1::DealiiDiscretisationPolicy policy;
    policy.state_degree = 2;
    const compiler::v1::DealiiCompiler compiler;
    contract::require(compiler.validate(specification, policy).valid(),
                      "point-sensor v1 graph did not validate for deal.II");

    const compiler::v1::DealiiDataBindings<dim> bindings{
      forcing,
      desired_state,
      1.0,
      0.5,
      0.2,
      test_binding_provenance("point_sensor")};
    auto display_only_specification = specification;
    for (auto &requirement : display_only_specification.requirement_policies)
      if (requirement.id == "point_sensor_evaluation_policy" ||
          requirement.id == "point_sensor_transposition_policy")
        requirement.selected_policy.clear();
    contract::require(compiler.validate(display_only_specification, policy).valid(),
                      "point-sensor compiler validation depended on policy prose");
    const auto display_only_compilation = compiler.compile(
      display_only_specification, triangulation, bindings, policy);
    contract::require(display_only_compilation.succeeded(),
                      "point-sensor lowering depended on policy prose");
    auto outside_mesh_specification = specification;
    component_by_id(outside_mesh_specification.regions, "point_sensor_region")
      .point_coordinates.front()[0] = 1.25;
    const auto outside_mesh = compiler.compile(outside_mesh_specification,
                                               triangulation,
                                               bindings,
                                               policy);
    contract::require(!outside_mesh.succeeded(),
                      "point-sensor compiler accepted an out-of-mesh coordinate");
    test_support::require_exact_diagnostic(
      outside_mesh.diagnostics,
      semantic::v1::DiagnosticCategory::lowerability,
      "point_sensor_region",
      "point_sensor_inside_mesh",
      "point-sensor compiler did not reject an out-of-mesh coordinate");

    dealii::Triangulation<dim> extra_boundary_triangulation;
    dealii::GridGenerator::hyper_cube(extra_boundary_triangulation);
    extra_boundary_triangulation.refine_global(2);
    for (auto cell = extra_boundary_triangulation.begin_active();
         cell != extra_boundary_triangulation.end();
         ++cell)
      for (unsigned int face = 0;
           face < dealii::GeometryInfo<dim>::faces_per_cell;
           ++face)
        if (cell->face(face)->at_boundary())
          cell->face(face)->set_boundary_id(
            cell->face(face)->center()[0] < 1e-12 ? 0 : 1);
    const auto extra_boundary = compiler.compile(
      specification,
      extra_boundary_triangulation,
      bindings,
      policy);
    test_support::require_exact_diagnostic(
      extra_boundary.diagnostics,
      semantic::v1::DiagnosticCategory::lowerability,
      "state",
      "p53_complete_fixed_dirichlet_boundary",
      "point-sensor compiler did not reject an exterior id outside the fixed region");

    const auto absent_fixed_boundary_specification =
      semantic::v1::make_point_sensor_scalar_diffusion_reaction_problem(
        sensor_coordinates, {99});
    const auto absent_fixed_boundary = compiler.compile(
      absent_fixed_boundary_specification,
      triangulation,
      bindings,
      policy);
    test_support::require_exact_diagnostic(
      absent_fixed_boundary.diagnostics,
      semantic::v1::DiagnosticCategory::lowerability,
      "state",
      "fixed_dirichlet_boundary_presence",
      "point-sensor compiler did not reject an absent fixed boundary id");

    const auto compilation = compiler.compile(specification,
                                              triangulation,
                                              bindings,
                                              policy);
    contract::require(compilation.succeeded(),
                      "point-sensor v1 compilation failed");
    require_compiled_hessian_evidence(*compilation.problem, "point sensor");
    test_support::require_manifest_compatibility_equal(
      compilation.problem->manifest(),
      display_only_compilation.problem->manifest(),
      "point-sensor display-policy edit");
    const auto &model = compilation.problem->executable_model();
    const auto reduced = compilation.problem->make_reduced_dto();
    const auto *point_model =
      dynamic_cast<const compiler::v1::detail::ScalarComponentModel<dim> *>(
        &model);
    contract::require(point_model != nullptr,
                      "point-sensor compilation did not produce its scalar target");

    dealii::Vector<double> control_values(model.variable_layout()->dimension(1));
    const Primal control(model.variable_layout()->single_block(1, "control"),
                         {std::move(control_values)});
    const auto evaluation = compilation.problem->make_reduced_dto().evaluate(control);
    require_close(model.residual(evaluation.full_point).block(0).l2_norm(),
                  0.0,
                  1e-11,
                  "point-sensor manufactured state residual");
    require_close(evaluation.objective_value,
                  0.5 * (0.23 * 0.23 + 0.71 * 0.71),
                  1e-11,
                  "point-sensor tracking value");

    const auto values = point_model->point_sensor_values(evaluation.full_point);

    dealii::Vector<double> state_tangent(model.variable_layout()->dimension(0));
    for (dealii::types::global_dof_index index = 0;
         index < state_tangent.size();
         ++index)
      state_tangent[index] = 0.03 * static_cast<double>(index + 1);
    dealii::Vector<double> control_tangent(model.variable_layout()->dimension(1));
    const Primal tangent(model.variable_layout(),
                         {std::move(state_tangent), std::move(control_tangent)});
    const Covector sensor_jvp = point_model->point_sensor_jvp(tangent);
    dealii::Vector<double> sensor_seed_values(2);
    sensor_seed_values[0] = 1.3;
    sensor_seed_values[1] = -0.7;
    const Primal sensor_seed(sensor_jvp.layout(),
                             {std::move(sensor_seed_values)});
    const Covector sensor_vjp = point_model->point_sensor_vjp({1.3, -0.7});
    const Primal state_tangent_block =
      contract::extract_primal_block(tangent, 0, "state");
    require_close(contract::pair(sensor_jvp, sensor_seed),
                  contract::pair(sensor_vjp, state_tangent_block),
                  1e-11,
                  "point-sensor evaluation JVP/VJP pairing");

    const Covector expected_state_derivative =
      point_model->point_sensor_vjp({-0.23, -0.71});
    const Covector actual_objective_derivative =
      model.objective_derivative(evaluation.full_point);
    require_covector_close(
      contract::extract_covector_block(actual_objective_derivative, 0, "state"),
      expected_state_derivative,
      1e-11,
      "point-sensor very-weak objective transpose");

    const Covector adjoint_rhs = actual_objective_derivative;
    const Covector adjoint_pullback =
      model.residual_vjp(evaluation.full_point, evaluation.adjoint);
    dealii::Vector<double> adjoint_dual_residual = adjoint_pullback.block(0);
    adjoint_dual_residual.add(-1.0, adjoint_rhs.block(0));
    require_close(adjoint_dual_residual.l2_norm(),
                  0.0,
                  1e-10,
                  "point-sensor very-weak adjoint dual residual");

    dealii::Vector<double> direction_values(control.layout()->dimension(0));
    for (dealii::types::global_dof_index index = 0;
         index < direction_values.size();
         ++index)
      direction_values[index] =
        (index % 2 == 0 ? 0.04 : -0.03) * static_cast<double>(index + 1);
    const Primal direction(control.layout(), {std::move(direction_values)});
    const double reduced_directional_derivative =
      contract::pair(evaluation.reduced_derivative, direction);
    constexpr double derivative_step = 1e-6;
    const double centered_reduced_derivative =
      (reduced.evaluate(shifted(control, direction, derivative_step))
           .objective_value -
       reduced.evaluate(shifted(control, direction, -derivative_step))
           .objective_value) /
      (2.0 * derivative_step);
    require_close(centered_reduced_derivative,
                  reduced_directional_derivative,
                  2e-7,
                  "point-sensor reduced objective directional derivative");
    const auto remainder = [&](const double step) {
      return std::abs(reduced.evaluate(shifted(control, direction, step))
                        .objective_value -
                      evaluation.objective_value -
                      step * reduced_directional_derivative);
    };
    const double coarse_remainder = remainder(1e-3);
    const double fine_remainder = remainder(5e-4);
    contract::require(coarse_remainder > 1e-12 &&
                        fine_remainder <= 0.26 * coarse_remainder + 1e-13,
                      "point-sensor reduced Taylor remainder is not quadratic");

    const auto &manifest = compilation.problem->manifest();
    const auto point_map = std::find_if(
      manifest.resolved_decision.realized_maps.begin(),
      manifest.resolved_decision.realized_maps.end(),
      [](const compiler::v1::CompiledRealizedMapRecord &map) {
        return map.semantic_id == "state_observation" &&
               map.realization_id == "ordered_point_sensor_values";
      });
    const auto point_space = std::find_if(
      manifest.resolved_decision.spaces.begin(),
      manifest.resolved_decision.spaces.end(),
      [](const compiler::v1::CompiledSpaceRecord &space) {
        return space.semantic_id == "state_observation_space" &&
               space.role == semantic::v1::SpaceRole::observation;
      });
    contract::require(
      point_space != manifest.resolved_decision.spaces.end() &&
        point_space->dimension == values.size() &&
        point_space->dimension == sensor_jvp.block(0).size(),
      "point-sensor manifest recorded a dimension different from its realized output");
    contract::require(
      manifest.compatibility.compiler_id == "nmopt.compiler.v1.dealii.point_sensor" &&
        point_map != manifest.resolved_decision.realized_maps.end() &&
        point_map->output_dimension == values.size() &&
        point_map->output_layout.find("sensor values") != std::string::npos &&
        manifest.compatibility.observation_realisation.find("immutable physical coordinates") !=
          std::string::npos &&
        manifest.resolved_decision.transposition_realisation.has_value() &&
        manifest.resolved_decision.transposition_realisation->diffusion_data_id == "diffusion" &&
        manifest.resolved_decision.transposition_realisation->reaction_data_id == "reaction" &&
        manifest.compatibility.data_rule.find("assembled C_h^T point-load transpose") !=
          std::string::npos &&
        std::any_of(manifest.compatibility.declared_assumptions.begin(),
                    manifest.compatibility.declared_assumptions.end(),
                    [](const std::string &assumption) {
                      return assumption.find("very-weak adjoint source") !=
                             std::string::npos;
                    }) &&
        std::find(manifest.compatibility.lowering_handler_records.begin(),
                  manifest.compatibility.lowering_handler_records.end(),
                  "state_observation <- dealii.scalar.observation.point_sensor") !=
          manifest.compatibility.lowering_handler_records.end(),
      "point-sensor compilation manifest omitted its finite transpose policy");
  }

  template <int dim>
  void
  run_normal_flux_contract_test()
  {
    static_assert(dim == 2, "The normal-flux oracle is two-dimensional");
    dealii::Triangulation<dim> triangulation;
    dealii::GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(2);
    for (auto cell = triangulation.begin_active();
         cell != triangulation.end();
         ++cell)
      for (unsigned int face = 0;
           face < dealii::GeometryInfo<dim>::faces_per_cell;
           ++face)
        if (cell->face(face)->at_boundary())
          cell->face(face)->set_boundary_id(
            cell->face(face)->center()[0] > 1.0 - 1e-12 ? 7 : 0);

    constexpr double reaction = 0.5;
    const EnergyPolynomialForcing<dim> forcing(reaction);
    const RightBoundaryNormalFluxFunction<dim> desired_state(0.5);
    const auto specification =
      semantic::v1::make_normal_flux_scalar_diffusion_reaction_problem(
        {7}, {0, 7});
    compiler::v1::DealiiDiscretisationPolicy policy;
    policy.state_degree = 2;
    const compiler::v1::DealiiCompiler compiler;
    contract::require(compiler.validate(specification, policy).valid(),
                      "normal-flux v1 graph did not validate for deal.II");

    const compiler::v1::DealiiDataBindings<dim> bindings{
      forcing,
      desired_state,
      1.0,
      reaction,
      0.2,
      test_binding_provenance("normal_flux")};
    const EnergyPolynomialForcing<dim> nonunit_forcing(reaction, 2.0);
    const compiler::v1::DealiiDataBindings<dim> nonunit_bindings{
      nonunit_forcing,
      desired_state,
      2.0,
      reaction,
      0.2,
      test_binding_provenance("normal_flux_nonunit")};
    const auto nonunit_compilation = compiler.compile(
      specification, triangulation, nonunit_bindings, policy);
    contract::require(
      nonunit_compilation.succeeded() &&
        nonunit_compilation.problem->manifest().resolved_decision.transposition_realisation
          .has_value() &&
        nonunit_compilation.problem->manifest().resolved_decision.transposition_realisation
            ->diffusion_data_id == "diffusion" &&
        nonunit_compilation.problem->manifest().resolved_decision.transposition_realisation
            ->reaction_data_id == "reaction" &&
        nonunit_compilation.problem->manifest().compatibility.data_rule.find(
          "T=-kappa Delta+rI with kappa <- diffusion and r <- reaction") !=
          std::string::npos,
      "normal-flux compiler did not preserve coefficient provenance for non-unit diffusion");
    auto display_only_specification = specification;
    for (auto &requirement : display_only_specification.requirement_policies)
      if (requirement.id == "normal_flux_orientation_policy" ||
          requirement.id == "normal_flux_evaluation_policy" ||
          requirement.id == "normal_flux_transposition_policy")
        requirement.selected_policy.clear();
    contract::require(compiler.validate(display_only_specification, policy).valid(),
                      "normal-flux compiler validation depended on policy prose");
    const auto display_only_compilation = compiler.compile(
      display_only_specification, triangulation, bindings, policy);
    contract::require(display_only_compilation.succeeded(),
                      "normal-flux lowering depended on policy prose");

    const auto omitted_fixed_boundary_specification =
      semantic::v1::make_normal_flux_scalar_diffusion_reaction_problem(
        {7}, {0});
    const auto omitted_fixed_boundary = compiler.compile(
      omitted_fixed_boundary_specification,
      triangulation,
      bindings,
      policy);
    test_support::require_exact_diagnostic(
      omitted_fixed_boundary.diagnostics,
      semantic::v1::DiagnosticCategory::structural,
      "state_observation",
      "normal_flux_fixed_boundary_subset",
      "normal-flux compiler accepted an observed boundary outside the fixed region");

    const auto absent_normal_flux_specification =
      semantic::v1::make_normal_flux_scalar_diffusion_reaction_problem(
        {9}, {0, 9});
    const auto absent_normal_flux = compiler.compile(
      absent_normal_flux_specification,
      triangulation,
      bindings,
      policy);
    test_support::require_exact_diagnostic(
      absent_normal_flux.diagnostics,
      semantic::v1::DiagnosticCategory::lowerability,
      "normal_flux_boundary",
      "normal_flux_boundary_presence",
      "normal-flux compiler did not reject an absent observed boundary id");

    dealii::Triangulation<dim> extra_boundary_triangulation;
    dealii::GridGenerator::hyper_cube(extra_boundary_triangulation);
    extra_boundary_triangulation.refine_global(2);
    for (auto cell = extra_boundary_triangulation.begin_active();
         cell != extra_boundary_triangulation.end();
         ++cell)
      for (unsigned int face = 0;
           face < dealii::GeometryInfo<dim>::faces_per_cell;
           ++face)
        if (cell->face(face)->at_boundary())
          {
            const auto center = cell->face(face)->center();
            cell->face(face)->set_boundary_id(
              center[0] > 1.0 - 1e-12 ? 7 :
              center[1] > 1.0 - 1e-12 ? 2 : 0);
          }
    const auto extra_boundary = compiler.compile(
      specification,
      extra_boundary_triangulation,
      bindings,
      policy);
    test_support::require_exact_diagnostic(
      extra_boundary.diagnostics,
      semantic::v1::DiagnosticCategory::lowerability,
      "state",
      "p53_complete_fixed_dirichlet_boundary",
      "normal-flux compiler did not reject an exterior id outside the fixed region");
    const auto compilation = compiler.compile(specification,
                                              triangulation,
                                              bindings,
                                              policy);
    contract::require(compilation.succeeded(),
                      "normal-flux v1 compilation failed");
    require_compiled_hessian_evidence(*compilation.problem, "normal flux");
    test_support::require_manifest_compatibility_equal(
      compilation.problem->manifest(),
      display_only_compilation.problem->manifest(),
      "normal-flux display-policy edit");
    const auto &model = compilation.problem->executable_model();
    const auto reduced = compilation.problem->make_reduced_dto();
    const auto *normal_flux_model =
      dynamic_cast<const compiler::v1::detail::ScalarComponentModel<dim> *>(
        &model);
    contract::require(normal_flux_model != nullptr,
                      "normal-flux compilation did not produce its scalar target");

    dealii::Vector<double> control_values(model.variable_layout()->dimension(1));
    const Primal control(model.variable_layout()->single_block(1, "control"),
                         {std::move(control_values)});
    const auto evaluation = compilation.problem->make_reduced_dto().evaluate(control);
    require_close(model.residual(evaluation.full_point).block(0).l2_norm(),
                  0.0,
                  1e-10,
                  "normal-flux manufactured state residual");
    const auto &nonunit_model = nonunit_compilation.problem->executable_model();
    dealii::Vector<double> nonunit_control_values(
      nonunit_model.variable_layout()->dimension(1));
    const Primal nonunit_control(
      nonunit_model.variable_layout()->single_block(1, "control"),
      {std::move(nonunit_control_values)});
    const auto nonunit_evaluation =
      nonunit_compilation.problem->make_reduced_dto().evaluate(nonunit_control);
    require_close(
      nonunit_model.residual(nonunit_evaluation.full_point).block(0).l2_norm(),
      0.0,
      1e-10,
      "non-unit normal-flux manufactured state residual");
    require_primal_close(nonunit_evaluation.full_point,
                         evaluation.full_point,
                         1e-10,
                         "non-unit diffusion manufactured state");

    const std::vector<double> normal_flux_values =
      normal_flux_model->normal_flux_values(evaluation.full_point);
    const auto &normal_flux_weights =
      normal_flux_model->normal_flux_quadrature_weights();
    contract::require(!normal_flux_values.empty() &&
                        normal_flux_values.size() == normal_flux_weights.size(),
                      "normal-flux target did not assemble face quadrature values");
    dealii::Vector<double> normal_flux_value_vector(normal_flux_values.size());
    for (std::size_t index = 0; index < normal_flux_values.size(); ++index)
      normal_flux_value_vector[index] = normal_flux_values[index];
    contract::require(normal_flux_value_vector.l2_norm() > 1e-6,
                      "normal-flux value map was unexpectedly zero");
    require_close(evaluation.objective_value,
                  1.0 / 240.0,
                  2e-10,
                  "normal-flux nonzero mismatch tracking value");

    dealii::Vector<double> state_tangent(model.variable_layout()->dimension(0));
    for (dealii::types::global_dof_index index = 0;
         index < state_tangent.size();
         ++index)
      state_tangent[index] = 0.03 * static_cast<double>(index + 1);
    dealii::Vector<double> control_tangent(model.variable_layout()->dimension(1));
    const Primal tangent(model.variable_layout(),
                         {std::move(state_tangent), std::move(control_tangent)});
    const Covector normal_flux_jvp = normal_flux_model->normal_flux_jvp(tangent);
    std::vector<double> normal_flux_seed(normal_flux_values.size());
    for (std::size_t index = 0; index < normal_flux_seed.size(); ++index)
      normal_flux_seed[index] = 0.7 - 0.01 * static_cast<double>(index);
    const Covector normal_flux_vjp =
      normal_flux_model->normal_flux_vjp(normal_flux_seed);
    double weighted_jvp_pairing = 0.0;
    for (std::size_t index = 0; index < normal_flux_seed.size(); ++index)
      weighted_jvp_pairing += normal_flux_weights[index] *
                              normal_flux_jvp.block(0)[index] *
                              normal_flux_seed[index];
    const Primal state_tangent_block =
      contract::extract_primal_block(tangent, 0, "state");
    require_close(weighted_jvp_pairing,
                  contract::pair(normal_flux_vjp, state_tangent_block),
                  1e-11,
                  "normal-flux evaluation JVP/VJP quadrature pairing");

    dealii::Vector<double> zero_state(model.variable_layout()->dimension(0));
    dealii::Vector<double> zero_control(model.variable_layout()->dimension(1));
    const Primal zero_point(model.variable_layout(),
                            {std::move(zero_state), std::move(zero_control)});
    std::vector<double> objective_seed(normal_flux_values.size());
    for (std::size_t index = 0; index < objective_seed.size(); ++index)
      objective_seed[index] = -0.5 * normal_flux_values[index];
    const Covector expected_objective_state =
      normal_flux_model->normal_flux_vjp(objective_seed);
    const Covector actual_objective = model.objective_derivative(zero_point);
    require_covector_close(
      contract::extract_covector_block(actual_objective, 0, "state"),
      expected_objective_state,
      1e-11,
      "normal-flux very-weak objective transpose");

    const Covector state_objective =
      contract::extract_covector_block(actual_objective, 0, "state");
    const auto adjoint =
      normal_flux_model->solve_adjoint(zero_point, state_objective);
    const Covector adjoint_pullback =
      model.residual_vjp(zero_point, adjoint);
    dealii::Vector<double> adjoint_dual_residual = adjoint_pullback.block(0);
    adjoint_dual_residual.add(-1.0, state_objective.block(0));
    require_close(adjoint_dual_residual.l2_norm(),
                  0.0,
                  1e-10,
                  "normal-flux very-weak adjoint dual residual");

    dealii::Vector<double> direction_values(control.layout()->dimension(0));
    for (dealii::types::global_dof_index index = 0;
         index < direction_values.size();
         ++index)
      direction_values[index] =
        (index % 2 == 0 ? 0.04 : -0.03) * static_cast<double>(index + 1);
    const Primal direction(control.layout(), {std::move(direction_values)});
    const double reduced_directional_derivative =
      contract::pair(evaluation.reduced_derivative, direction);
    constexpr double derivative_step = 1e-6;
    const double centered_reduced_derivative =
      (reduced.evaluate(shifted(control, direction, derivative_step))
           .objective_value -
       reduced.evaluate(shifted(control, direction, -derivative_step))
           .objective_value) /
      (2.0 * derivative_step);
    require_close(centered_reduced_derivative,
                  reduced_directional_derivative,
                  2e-7,
                  "normal-flux reduced objective directional derivative");
    const auto remainder = [&](const double step) {
      return std::abs(reduced.evaluate(shifted(control, direction, step))
                        .objective_value -
                      evaluation.objective_value -
                      step * reduced_directional_derivative);
    };
    const double coarse_remainder = remainder(1e-3);
    const double fine_remainder = remainder(5e-4);
    contract::require(coarse_remainder > 1e-12 &&
                        fine_remainder <= 0.26 * coarse_remainder + 1e-13,
                      "normal-flux reduced Taylor remainder is not quadratic");

    const auto &manifest = compilation.problem->manifest();
    const auto normal_flux_map = std::find_if(
      manifest.resolved_decision.realized_maps.begin(),
      manifest.resolved_decision.realized_maps.end(),
      [](const compiler::v1::CompiledRealizedMapRecord &map) {
        return map.semantic_id == "state_observation" &&
               map.realization_id == "ordered_normal_flux_face_quadrature";
      });
    const auto normal_flux_space = std::find_if(
      manifest.resolved_decision.spaces.begin(),
      manifest.resolved_decision.spaces.end(),
      [](const compiler::v1::CompiledSpaceRecord &space) {
        return space.semantic_id == "state_observation_space" &&
               space.role == semantic::v1::SpaceRole::observation;
      });
    contract::require(
      normal_flux_space != manifest.resolved_decision.spaces.end() &&
        normal_flux_space->dimension == normal_flux_values.size() &&
        normal_flux_space->dimension == normal_flux_jvp.block(0).size(),
      "normal-flux manifest recorded the state dimension instead of its realized face output");
    contract::require(
      manifest.compatibility.compiler_id == "nmopt.compiler.v1.dealii.normal_flux" &&
        normal_flux_map != manifest.resolved_decision.realized_maps.end() &&
        normal_flux_map->output_dimension == normal_flux_values.size() &&
        normal_flux_map->pairing_realization.find("declared pairing") !=
          std::string::npos &&
        manifest.compatibility.observation_realisation.find("outward normal-flux") !=
          std::string::npos &&
        manifest.resolved_decision.transposition_realisation.has_value() &&
        manifest.resolved_decision.transposition_realisation->diffusion_data_id == "diffusion" &&
        manifest.resolved_decision.transposition_realisation->reaction_data_id == "reaction" &&
        manifest.compatibility.data_rule.find(
          "T=-kappa Delta+rI with kappa <- diffusion and r <- reaction") !=
          std::string::npos &&
        manifest.compatibility.data_rule.find("boundary face quadrature") !=
          std::string::npos &&
        std::any_of(manifest.compatibility.declared_assumptions.begin(),
                    manifest.compatibility.declared_assumptions.end(),
                    [](const std::string &assumption) {
                      return assumption.find("very-weak adjoint boundary source") !=
                             std::string::npos;
                    }) &&
        std::find(manifest.compatibility.lowering_handler_records.begin(),
                  manifest.compatibility.lowering_handler_records.end(),
                  "state_observation <- dealii.scalar.observation.normal_flux") !=
          manifest.compatibility.lowering_handler_records.end(),
      "normal-flux compilation manifest omitted its face transpose policy");
  }

  template <int dim>
  void
  run_weighted_boundary_trace_contract_test()
  {
    dealii::Triangulation<dim> triangulation;
    dealii::GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(2);
    for (auto cell = triangulation.begin_active();
         cell != triangulation.end();
         ++cell)
      for (unsigned int face = 0;
           face < dealii::GeometryInfo<dim>::faces_per_cell;
           ++face)
        if (cell->face(face)->at_boundary())
          {
            const double x = cell->face(face)->center()[0];
            cell->face(face)->set_boundary_id(x < 1e-12 ? 0 :
                                              x > 1.0 - 1e-12 ? 1 : 2);
          }

    const dealii::Functions::ConstantFunction<dim> forcing(0.5);
    const dealii::Functions::ConstantFunction<dim> weighted_target(0.2);
    const dealii::Functions::ConstantFunction<dim> comparison_target(0.1);
    const dealii::Functions::ConstantFunction<dim> boundary_weight(2.0);
    const auto weighted_specification =
      semantic::v1::make_weighted_boundary_trace_neumann_control_problem();
    const auto comparison_specification =
      semantic::v1::make_neumann_boundary_control_problem();
    compiler::v1::DealiiDiscretisationPolicy policy;
    policy.state_degree = 1;
    const compiler::v1::DealiiCompiler compiler;
    contract::require(compiler.validate(weighted_specification, policy).valid(),
                      "weighted boundary-trace graph did not validate for deal.II");

    const compiler::v1::DealiiDataBindings<dim> missing_weight_binding{
      forcing,
      weighted_target,
      1.0,
      0.5,
      0.1,
      test_binding_provenance("weighted_boundary_missing")};
    const auto missing_weight = compiler.compile(weighted_specification,
                                                 triangulation,
                                                 missing_weight_binding,
                                                 policy);
    test_support::require_exact_diagnostic(
      missing_weight.diagnostics,
      semantic::v1::DiagnosticCategory::lowerability,
      "boundary_weight",
      "weighted_boundary_trace_data_binding",
      "weighted trace compiler accepted a missing boundary-weight binding");

    const compiler::v1::DealiiDataBindings<dim> missing_weight_provenance{
      forcing,
      weighted_target,
      1.0,
      0.5,
      0.1,
      test_binding_provenance("weighted_boundary_missing_provenance"),
      std::nullopt,
      std::nullopt,
      compiler::v1::DealiiWeightedTraceDataBindings<dim>{boundary_weight, ""}};
    const auto missing_provenance = compiler.compile(
      weighted_specification,
      triangulation,
      missing_weight_provenance,
      policy);
    test_support::require_exact_diagnostic(
      missing_provenance.diagnostics,
      semantic::v1::DiagnosticCategory::lowerability,
      "boundary_weight",
      "boundary_weight_binding_provenance",
      "weighted trace compiler accepted missing weight provenance");

    const dealii::Functions::ConstantFunction<dim> vector_weight(2.0, 2);
    const compiler::v1::DealiiDataBindings<dim> vector_weight_binding{
      forcing,
      weighted_target,
      1.0,
      0.5,
      0.1,
      test_binding_provenance("weighted_boundary_vector_weight"),
      std::nullopt,
      std::nullopt,
      compiler::v1::DealiiWeightedTraceDataBindings<dim>{
        vector_weight, "test.weighted_boundary.vector_weight"}};
    const auto wrong_shape = compiler.compile(weighted_specification,
                                              triangulation,
                                              vector_weight_binding,
                                              policy);
    test_support::require_exact_diagnostic(
      wrong_shape.diagnostics,
      semantic::v1::DiagnosticCategory::lowerability,
      "boundary_weight",
      "scalar_boundary_weight_binding",
      "weighted trace compiler accepted a multi-component weight Function");

    const compiler::v1::DealiiDataBindings<dim> weighted_bindings{
      forcing,
      weighted_target,
      1.0,
      0.5,
      0.1,
      test_binding_provenance("weighted_boundary"),
      std::nullopt,
      std::nullopt,
      compiler::v1::DealiiWeightedTraceDataBindings<dim>{
        boundary_weight, "test.weighted_boundary.weight"}};
    const compiler::v1::DealiiDataBindings<dim> comparison_bindings{
      forcing,
      comparison_target,
      1.0,
      0.5,
      0.1,
      test_binding_provenance("weighted_boundary_comparison")};
    const auto weighted = compiler.compile(weighted_specification,
                                           triangulation,
                                           weighted_bindings,
                                           policy);
    const auto comparison = compiler.compile(comparison_specification,
                                             triangulation,
                                             comparison_bindings,
                                             policy);
    contract::require(weighted.succeeded() && comparison.succeeded(),
                      "weighted boundary-trace comparison compilation failed");

    const auto &weighted_model = weighted.problem->executable_model();
    const auto &comparison_model = comparison.problem->executable_model();
    const auto *weighted_neumann_model =
      dynamic_cast<const compiler::v1::detail::NeumannBoundaryControlModel<dim> *>(
        &weighted_model);
    const auto *comparison_neumann_model =
      dynamic_cast<const compiler::v1::detail::NeumannBoundaryControlModel<dim> *>(
        &comparison_model);
    contract::require(weighted_neumann_model != nullptr,
                      "weighted trace compilation did not produce its Neumann target");
    contract::require(comparison_neumann_model != nullptr,
                      "comparison compilation did not produce its Neumann target");
    const dealii::QGauss<dim - 1> weighted_boundary_quadrature(
      policy.state_degree + 2);
    std::size_t expected_weighted_trace_samples = 0;
    for (auto cell = triangulation.begin_active();
         cell != triangulation.end();
         ++cell)
      for (unsigned int face = 0;
           face < dealii::GeometryInfo<dim>::faces_per_cell;
           ++face)
        if (cell->face(face)->at_boundary() &&
            cell->face(face)->boundary_id() == 2)
          expected_weighted_trace_samples += weighted_boundary_quadrature.size();
    contract::require(expected_weighted_trace_samples > 0,
                      "weighted trace test found no observation boundary samples");
    dealii::Vector<double> zero_control_values(
      weighted_model.variable_layout()->dimension(1));
    const Primal zero_control(
      weighted_model.variable_layout()->single_block(1, "control"),
      {std::move(zero_control_values)});
    const auto weighted_reduced = weighted.problem->make_reduced_dto();
    const auto comparison_reduced = comparison.problem->make_reduced_dto();
    const auto weighted_evaluation = weighted_reduced.evaluate(zero_control);
    const auto comparison_evaluation = comparison_reduced.evaluate(zero_control);
    require_primal_close(weighted_evaluation.full_point,
                         comparison_evaluation.full_point,
                         1e-12,
                         "weighted trace changed the Neumann state equation");
    require_covector_close(
      weighted_model.residual(weighted_evaluation.full_point),
      comparison_model.residual(weighted_evaluation.full_point),
      1e-12,
      "weighted trace changed the Neumann residual value");

    dealii::Vector<double> state_tangent(
      weighted_model.variable_layout()->dimension(0));
    dealii::Vector<double> control_tangent(
      weighted_model.variable_layout()->dimension(1));
    for (dealii::types::global_dof_index index = 0;
         index < state_tangent.size();
         ++index)
      state_tangent[index] = 0.01 * static_cast<double>(index + 1);
    for (dealii::types::global_dof_index index = 0;
         index < control_tangent.size();
         ++index)
      control_tangent[index] =
        (index % 2 == 0 ? 0.03 : -0.02) * static_cast<double>(index + 1);
    const Primal tangent(weighted_model.variable_layout(),
                         {std::move(state_tangent),
                          std::move(control_tangent)});
    dealii::Vector<double> seed_values(weighted_model.test_layout()->dimension(0));
    for (dealii::types::global_dof_index index = 0;
         index < seed_values.size();
         ++index)
      seed_values[index] = 0.04 * static_cast<double>(index + 1);
    const Primal seed(weighted_model.test_layout(), {std::move(seed_values)});
    require_covector_close(
      weighted_model.residual_jvp(weighted_evaluation.full_point, tangent),
      comparison_model.residual_jvp(weighted_evaluation.full_point, tangent),
      1e-12,
      "weighted trace changed the Neumann residual JVP");
    require_covector_close(
      weighted_model.residual_vjp(weighted_evaluation.full_point, seed),
      comparison_model.residual_vjp(weighted_evaluation.full_point, seed),
      1e-12,
      "weighted trace changed the Neumann residual VJP");

    const auto weighted_trace_values =
      weighted_neumann_model->boundary_trace_values(
        weighted_evaluation.full_point);
    const auto comparison_trace_values =
      comparison_neumann_model->boundary_trace_values(
        comparison_evaluation.full_point);
    contract::require(
      expected_weighted_trace_samples > 0 &&
        weighted_trace_values.size() == expected_weighted_trace_samples &&
        comparison_trace_values.size() == expected_weighted_trace_samples,
      "boundary-trace actions did not expose the ordered face samples (weighted=" +
        std::to_string(weighted_trace_values.size()) +
        ", comparison=" + std::to_string(comparison_trace_values.size()) +
        ", expected=" + std::to_string(expected_weighted_trace_samples) + ")");
    for (std::size_t index = 0; index < weighted_trace_values.size(); ++index)
      require_close(weighted_trace_values[index],
                    2.0 * comparison_trace_values[index],
                    1e-12,
                    "weighted boundary-trace value omitted its weight");
    const auto weighted_trace_jvp =
      weighted_neumann_model->boundary_trace_jvp(tangent);
    const auto &weighted_trace_weights =
      weighted_neumann_model->boundary_trace_quadrature_weights();
    contract::require(weighted_trace_weights.size() == weighted_trace_values.size(),
                      "weighted boundary-trace map omitted its pairing weights");
    std::vector<double> boundary_trace_seed(weighted_trace_values.size());
    for (std::size_t index = 0; index < boundary_trace_seed.size(); ++index)
      boundary_trace_seed[index] = 0.02 * static_cast<double>(index + 1);
    const auto weighted_trace_vjp =
      weighted_neumann_model->boundary_trace_vjp(boundary_trace_seed);
    double weighted_trace_jvp_pairing = 0.0;
    for (std::size_t index = 0; index < boundary_trace_seed.size(); ++index)
      weighted_trace_jvp_pairing += weighted_trace_weights[index] *
                                    weighted_trace_jvp.block(0)[index] *
                                    boundary_trace_seed[index];
    require_close(
      weighted_trace_jvp_pairing,
      contract::pair(
        weighted_trace_vjp,
        contract::extract_primal_block(tangent, 0, "state")),
      1e-11,
      "weighted boundary-trace value/JVP/VJP map is inconsistent");

    require_close(weighted_evaluation.objective_value,
                  4.0 * comparison_evaluation.objective_value,
                  1e-11,
                  "weighted trace value did not realize h times the trace");
    Covector expected_full_derivative = comparison_model.objective_derivative(
      comparison_evaluation.full_point);
    expected_full_derivative.scale_block(0, 4.0);
    require_covector_close(
      weighted_model.objective_derivative(weighted_evaluation.full_point),
      expected_full_derivative,
      1e-11,
      "weighted trace transpose pullback did not include both weight factors");
    Covector expected_reduced_derivative =
      comparison_evaluation.reduced_derivative;
    expected_reduced_derivative.scale_block(0, 4.0);
    require_covector_close(weighted_evaluation.reduced_derivative,
                           expected_reduced_derivative,
                           1e-9,
                           "weighted trace reduced pullback is inconsistent");

    constexpr double derivative_step = 1e-7;
    const double centered_jvp =
      (weighted_model.objective(
         shifted(weighted_evaluation.full_point, tangent, derivative_step)) -
       weighted_model.objective(
         shifted(weighted_evaluation.full_point, tangent, -derivative_step))) /
      (2.0 * derivative_step);
    require_close(
      centered_jvp,
      contract::pair(
        weighted_model.objective_derivative(weighted_evaluation.full_point),
        tangent),
      2e-7,
      "weighted trace value/JVP/VJP chain rule");

    dealii::Vector<double> control_direction_values(zero_control.layout()->dimension(0));
    for (dealii::types::global_dof_index index = 0;
         index < control_direction_values.size();
         ++index)
      control_direction_values[index] =
        (index % 2 == 0 ? 0.04 : -0.03) * static_cast<double>(index + 1);
    const Primal control_direction(zero_control.layout(),
                                   {std::move(control_direction_values)});
    const double directional_derivative =
      contract::pair(weighted_evaluation.reduced_derivative,
                     control_direction);
    const auto remainder = [&](const double step) {
      return std::abs(
        weighted_reduced.evaluate(shifted(zero_control,
                                          control_direction,
                                          step))
            .objective_value -
        weighted_evaluation.objective_value - step * directional_derivative);
    };
    const double coarse_remainder = remainder(1e-3);
    const double fine_remainder = remainder(5e-4);
    contract::require(coarse_remainder > 1e-12 &&
                        fine_remainder <= 0.26 * coarse_remainder + 1e-13,
                      "weighted boundary-trace reduced Taylor remainder is not quadratic");

    const auto &weighted_metric = weighted.problem->metric();
    const auto &comparison_metric = comparison.problem->metric();
    require_covector_close(weighted_metric.apply(control_direction),
                           comparison_metric.apply(control_direction),
                           1e-12,
                           "weighted trace changed the facewise L2 control metric");

    const auto &manifest = weighted.problem->manifest();
    const auto weighted_map = std::find_if(
      manifest.resolved_decision.realized_maps.begin(),
      manifest.resolved_decision.realized_maps.end(),
      [](const compiler::v1::CompiledRealizedMapRecord &map) {
        return map.semantic_id == "weighted_state_boundary_trace" &&
               map.realization_id == "ordered_boundary_face_quadrature_trace";
      });
    const auto weight_record = std::find_if(
      manifest.resolved_decision.bindings.begin(),
      manifest.resolved_decision.bindings.end(),
      [](const compiler::v1::CompiledBindingRecord &record) {
        return record.semantic_id == "boundary_weight";
      });
    contract::require(
      weight_record != manifest.resolved_decision.bindings.end() &&
        weight_record->provenance == "test.weighted_boundary.weight" &&
        weight_record->representation.find("boundary face quadrature") !=
          std::string::npos &&
        manifest.compatibility.compiler_id ==
          "nmopt.compiler.v1.dealii.weighted_boundary_trace" &&
        manifest.compatibility.observation_realisation.find("weighted boundary trace") !=
          std::string::npos &&
        manifest.compatibility.data_rule.find("desired-state and fixed boundary-weight") !=
          std::string::npos &&
        manifest.compatibility.data_rule.find("boundary face quadrature") !=
          std::string::npos &&
        manifest.resolved_decision.metric_record.realisation_id == "l2_facewise" &&
        std::find(manifest.compatibility.lowering_handler_records.begin(),
                  manifest.compatibility.lowering_handler_records.end(),
                  "weighted_state_boundary_trace <- "
                  "dealii.neumann.observation.weighted_boundary_trace") !=
          manifest.compatibility.lowering_handler_records.end(),
      "weighted trace manifest omitted weight, target, quadrature, or metric provenance");
    contract::require(
      weighted_map != manifest.resolved_decision.realized_maps.end() &&
        weighted_map->output_dimension ==
          weighted_trace_values.size() &&
        weighted_map->output_dimension == expected_weighted_trace_samples &&
        weighted_map->output_dimension == weighted_trace_jvp.block(0).size() &&
        weighted_map->output_layout.find("boundary face") != std::string::npos,
      "weighted trace realized map is missing its face-quadrature output");
    const auto weighted_observation_record = std::find_if(
      manifest.resolved_decision.observations.begin(),
      manifest.resolved_decision.observations.end(),
      [](const compiler::v1::CompiledRealisationRecord &record) {
        return record.semantic_id == "weighted_state_boundary_trace";
      });
    contract::require(
      weighted_observation_record != manifest.resolved_decision.observations.end() &&
        weighted_observation_record->realisation_id ==
          "weighted_state_boundary_trace_fe_qgauss" &&
        std::find(weighted_observation_record->input_ids.begin(),
                  weighted_observation_record->input_ids.end(),
                  "boundary_weight") != weighted_observation_record->input_ids.end(),
      "weighted trace resolved provenance omitted its typed data selection");
    require_constraint_realisation(
      manifest, "none", "weighted boundary-trace manifest projection");
  }

