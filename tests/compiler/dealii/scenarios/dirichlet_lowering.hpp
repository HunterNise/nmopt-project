#pragma once

// These scenario implementations are deliberately included into one heavy
// deal.II/compiler translation unit to improve source navigation without
// multiplying compilation cost.

  template <int dim>
  void
  run_fixed_dirichlet_contract_test()
  {
    dealii::Triangulation<dim> triangulation;
    dealii::GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(2);

    // y_phys = 1 is the manufactured state for -Delta y + 0.5 y = 0.5
    // with fixed Dirichlet data y = 1. The state coordinates are the
    // independent y_hat values, while the executable actions see y_phys.
    const dealii::Functions::ConstantFunction<dim> forcing(0.5);
    const dealii::Functions::ConstantFunction<dim> desired_state(0.25);
    const dealii::Functions::ConstantFunction<dim> fixed_dirichlet_data(1.0);
    const dealii::Functions::ConstantFunction<dim> changed_dirichlet_data(2.0);
    const auto specification =
      semantic::v1::make_fixed_dirichlet_scalar_diffusion_reaction_problem();
    compiler::v1::DealiiDiscretisationPolicy policy;
    policy.state_degree = 1;
    const compiler::v1::DealiiCompiler compiler;
    contract::require(compiler.validate(specification, policy).valid(),
                      "fixed-Dirichlet v1 graph did not validate for deal.II");

    const compiler::v1::DealiiDataBindings<dim> missing_lifting_binding{
      forcing,
      desired_state,
      1.0,
      0.5,
      0.1,
      test_binding_provenance("fixed_dirichlet_missing")};
    const auto missing_lifting = compiler.compile(specification,
                                                  triangulation,
                                                  missing_lifting_binding,
                                                  policy);
    contract::require(
      !missing_lifting.succeeded(),
      "v1 compiler did not diagnose missing fixed-Dirichlet data");
    test_support::require_exact_diagnostic(
      missing_lifting.diagnostics,
      semantic::v1::DiagnosticCategory::lowerability,
      "state",
      "fixed_dirichlet_data_binding",
      "v1 compiler did not identify missing fixed-Dirichlet data");

    const dealii::Functions::ConstantFunction<dim> vector_fixed_data(1.0, 2);
    auto vector_fixed_bindings = compiler::v1::DealiiDataBindings<dim>{
      forcing,
      desired_state,
      1.0,
      0.5,
      0.1,
      test_binding_provenance("fixed_dirichlet_vector")};
    vector_fixed_bindings.fixed_dirichlet_data = std::cref(vector_fixed_data);
    const auto rejected_fixed_shape = compiler.compile(specification,
                                                       triangulation,
                                                       vector_fixed_bindings,
                                                       policy);
    test_support::require_exact_diagnostic(
      rejected_fixed_shape.diagnostics,
      semantic::v1::DiagnosticCategory::lowerability,
      "fixed_dirichlet_data",
      "scalar_function_binding_shape",
      "v1 compiler did not route fixed-Dirichlet Function shape through the resolved binding request");

    auto bindings = compiler::v1::DealiiDataBindings<dim>{
      forcing,
      desired_state,
      1.0,
      0.5,
      0.1,
      test_binding_provenance("fixed_dirichlet", true)};
    bindings.fixed_dirichlet_data = std::cref(fixed_dirichlet_data);
    const auto compilation = compiler.compile(specification,
                                              triangulation,
                                              bindings,
                                              policy);
    contract::require(compilation.succeeded(),
                      "v1 fixed-Dirichlet compilation failed");

    const auto &model = compilation.problem->executable_model();
    const auto reduced = compilation.problem->make_reduced_dto();
    dealii::Vector<double> control_values(model.variable_layout()->dimension(1));
    const Primal control(model.variable_layout()->single_block(1, "control"),
                         {std::move(control_values)});
    const auto evaluation = reduced.evaluate(control);
    const Covector state_residual = model.residual(evaluation.full_point);
    require_close(state_residual.block(0).l2_norm(),
                  0.0,
                  1e-11,
                  "fixed-Dirichlet reconstructed state residual");
    require_close(evaluation.objective_value,
                  0.28125,
                  1e-11,
                  "fixed-Dirichlet physical state tracking value");

    dealii::Vector<double> state_tangent(model.variable_layout()->dimension(0));
    dealii::Vector<double> control_tangent(model.variable_layout()->dimension(1));
    for (dealii::types::global_dof_index index = 0;
         index < state_tangent.size();
         ++index)
      state_tangent[index] = 0.02 * static_cast<double>(index + 1);
    for (dealii::types::global_dof_index index = 0;
         index < control_tangent.size();
         ++index)
      control_tangent[index] = -0.03 * static_cast<double>(index + 1);
    const Primal tangent(model.variable_layout(),
                         {std::move(state_tangent), std::move(control_tangent)});
    dealii::Vector<double> seed_values(model.test_layout()->dimension(0));
    for (dealii::types::global_dof_index index = 0;
         index < seed_values.size();
         ++index)
      seed_values[index] = 0.04 * static_cast<double>(index + 1);
    const Primal test_seed(model.test_layout(), {std::move(seed_values)});

    const Covector jvp = model.residual_jvp(evaluation.full_point, tangent);
    const Covector vjp = model.residual_vjp(evaluation.full_point, test_seed);
    require_close(contract::pair(jvp, test_seed),
                  contract::pair(vjp, tangent),
                  1e-11,
                  "fixed-Dirichlet reconstruction JVP/VJP pairing");

    constexpr double derivative_step = 1e-7;
    const Covector residual_at_step = model.residual(
      shifted(evaluation.full_point, tangent, derivative_step));
    for (std::size_t block = 0; block < residual_at_step.n_blocks(); ++block)
      {
        dealii::Vector<double> finite_difference = residual_at_step.block(block);
        finite_difference.add(-1.0, state_residual.block(block));
        finite_difference *= 1.0 / derivative_step;
        finite_difference.add(-1.0, jvp.block(block));
        require_close(finite_difference.l2_norm(),
                      0.0,
                      1e-7,
                      "fixed-Dirichlet reconstruction residual JVP");
      }
    const double objective_difference =
      model.objective(shifted(evaluation.full_point, tangent, derivative_step)) -
      evaluation.objective_value;
    require_close(objective_difference / derivative_step,
                  contract::pair(model.objective_derivative(evaluation.full_point),
                                 tangent),
                  2e-7,
                  "fixed-Dirichlet physical objective derivative");

    dealii::Vector<double> control_direction_values(
      control.layout()->dimension(0));
    for (dealii::types::global_dof_index index = 0;
         index < control_direction_values.size();
         ++index)
      control_direction_values[index] =
        (index % 2 == 0 ? 0.05 : -0.04) * static_cast<double>(index + 1);
    const Primal control_direction(control.layout(),
                                   {std::move(control_direction_values)});
    const double directional_derivative =
      contract::pair(evaluation.reduced_derivative, control_direction);
    const auto remainder = [&](const double step) {
      return std::abs(reduced.evaluate(shifted(control, control_direction, step))
                        .objective_value -
                      evaluation.objective_value - step * directional_derivative);
    };
    const double coarse_remainder = remainder(1e-3);
    const double fine_remainder = remainder(5e-4);
    contract::require(coarse_remainder > 1e-12 &&
                        fine_remainder <= 0.26 * coarse_remainder + 1e-13,
                      "fixed-Dirichlet reduced Taylor remainder is not quadratic");

    const auto *hessian = compilation.problem->reduced_hessian();
    contract::require(hessian != nullptr,
                      "fixed-Dirichlet compiled target omitted its Hessian capability");
    const Covector hessian_action =
      hessian->apply(control, control_direction);
    dealii::Vector<double> second_direction_values(
      control.layout()->dimension(0));
    for (dealii::types::global_dof_index index = 0;
         index < second_direction_values.size();
         ++index)
      second_direction_values[index] =
        (index % 3 == 0 ? -0.02 : 0.03) * static_cast<double>(index + 1);
    const Primal second_control_direction(
      control.layout(), {std::move(second_direction_values)});
    const Covector second_hessian_action =
      hessian->apply(control, second_control_direction);
    require_close(contract::pair(hessian_action, second_control_direction),
                  contract::pair(second_hessian_action, control_direction),
                  1e-10,
                  "fixed-Dirichlet compiled Hessian symmetry");

    constexpr double hessian_step = 1e-5;
    const Covector reduced_derivative_plus =
      reduced.evaluate(shifted(control, control_direction, hessian_step))
        .reduced_derivative;
    const Covector reduced_derivative_minus =
      reduced.evaluate(shifted(control, control_direction, -hessian_step))
        .reduced_derivative;
    Covector hessian_finite_difference = reduced_derivative_plus;
    hessian_finite_difference.add_scaled_block(
      0, -1.0, reduced_derivative_minus.block(0));
    hessian_finite_difference.scale_block(0, 1.0 / (2.0 * hessian_step));
    hessian_finite_difference.add_scaled_block(
      0, -1.0, hessian_action.block(0));
    require_close(hessian_finite_difference.block(0).l2_norm(),
                  0.0,
                  2e-8,
                  "fixed-Dirichlet compiled Hessian finite-difference action");

    auto changed_bindings = compiler::v1::DealiiDataBindings<dim>{
      forcing,
      desired_state,
      1.0,
      0.5,
      0.1,
      test_binding_provenance("changed_fixed_dirichlet", true)};
    changed_bindings.fixed_dirichlet_data = std::cref(changed_dirichlet_data);
    const auto changed_compilation = compiler.compile(specification,
                                                      triangulation,
                                                      changed_bindings,
                                                      policy);
    contract::require(changed_compilation.succeeded(),
                      "v1 fixed-Dirichlet recompilation failed for changed data");
    const auto changed_evaluation =
      changed_compilation.problem->make_reduced_dto().evaluate(control);
    contract::require(
      std::abs(changed_evaluation.objective_value - evaluation.objective_value) >
        1e-4,
      "changed fixed-Dirichlet data reused stale compiled values");

    const auto &manifest = compilation.problem->manifest();
    const auto fixed_map = std::find_if(
      manifest.resolved_decision.realized_maps.begin(),
      manifest.resolved_decision.realized_maps.end(),
      [](const compiler::v1::CompiledRealizedMapRecord &map) {
        return map.semantic_id == "fixed_dirichlet_reconstruction";
      });
    require_constraint_realisation(manifest, "none", "fixed-Dirichlet");
    contract::require(fixed_map != manifest.resolved_decision.realized_maps.end(),
                      "fixed-Dirichlet realized transformation is missing");
    contract::require(
      fixed_map->input_dimensions.size() == 1 &&
        fixed_map->input_dimensions.front() != fixed_map->output_dimension,
      "fixed-Dirichlet realized transformation dimensions are not distinct");
    contract::require(
        manifest.compatibility.lifting_realisation.find("y_phys = P_h y_hat + ell_0,h") !=
        std::string::npos &&
        manifest.compatibility.data_rule.find("boundary DoFs") != std::string::npos &&
        manifest.compatibility.transformation_ids.size() == 1 &&
        manifest.compatibility.lowering_handler_records.size() == 9 &&
        std::find(manifest.compatibility.lowering_handler_records.begin(),
                  manifest.compatibility.lowering_handler_records.end(),
                  "fixed_dirichlet_reconstruction <- "
                  "dealii.scalar.transformation.fixed_dirichlet") !=
          manifest.compatibility.lowering_handler_records.end(),
      "v1 fixed-Dirichlet compilation manifest is incomplete");
  }

  template <int dim>
  void
  run_dirichlet_control_contract_test()
  {
    dealii::Triangulation<dim> triangulation;
    dealii::GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(2);

    // y_phys = 1 is the manufactured state for -Delta y + 0.5 y = 0.5,
    // now with the entire Dirichlet trace supplied by the decision block.
    const dealii::Functions::ConstantFunction<dim> forcing(0.5);
    const dealii::Functions::ConstantFunction<dim> desired_state(0.25);
    const auto specification =
      semantic::v1::make_dirichlet_control_scalar_diffusion_reaction_problem();
    compiler::v1::DealiiDiscretisationPolicy policy;
    policy.state_degree = 1;
    const compiler::v1::DealiiCompiler compiler;
    contract::require(
      compiler.validate(specification, policy).valid(),
      "Dirichlet-control lifting v1 graph did not validate for deal.II");

    auto partial_boundary_specification = specification;
    partial_boundary_specification.regions.at(1).boundary_ids = {1};
    const compiler::v1::DealiiDataBindings<dim> bindings{
      forcing,
      desired_state,
      1.0,
      0.5,
      0.1,
      test_binding_provenance("dirichlet_control")};
    const auto partial_boundary = compiler.compile(partial_boundary_specification,
                                                   triangulation,
                                                   bindings,
                                                   policy);
    contract::require(
      !partial_boundary.succeeded(),
      "Dirichlet-control compiler did not reject an incomplete exterior boundary");
    test_support::require_exact_diagnostic(
      partial_boundary.diagnostics,
      semantic::v1::DiagnosticCategory::lowerability,
      "state",
      "complete_dirichlet_control_boundary",
      "Dirichlet-control compiler did not identify the incomplete exterior boundary");

    const auto compilation = compiler.compile(specification,
                                              triangulation,
                                              bindings,
                                              policy);
    contract::require(compilation.succeeded(),
                      "Dirichlet-control lifting v1 compilation failed");
    const auto &model = compilation.problem->executable_model();
    const auto *dirichlet_model =
      dynamic_cast<const compiler::v1::detail::DirichletControlLiftingModel<dim> *>(
        &model);
    contract::require(dirichlet_model != nullptr,
                      "Dirichlet-control compiler did not select its lifting target");
    const auto reduced = compilation.problem->make_reduced_dto();

    dealii::Vector<double> control_values(model.variable_layout()->dimension(1));
    for (dealii::types::global_dof_index index = 0;
         index < control_values.size();
         ++index)
      control_values[index] = 1.0;
    const Primal control(model.variable_layout()->single_block(1, "control"),
                         {std::move(control_values)});
    const auto evaluation = reduced.evaluate(control);
    require_close(model.residual(evaluation.full_point).block(0).l2_norm(),
                  0.0,
                  1e-11,
                  "Dirichlet-control lifted state residual");
    const dealii::Vector<double> physical_state =
      dirichlet_model->reconstruct_physical_state(evaluation.full_point);
    for (dealii::types::global_dof_index index = 0;
         index < physical_state.size();
         ++index)
      require_close(physical_state[index],
                    1.0,
                    1e-11,
                    "Dirichlet-control physical-state reconstruction");

    dealii::Vector<double> state_tangent(model.variable_layout()->dimension(0));
    dealii::Vector<double> control_tangent(model.variable_layout()->dimension(1));
    for (dealii::types::global_dof_index index = 0;
         index < state_tangent.size();
         ++index)
      state_tangent[index] = 0.02 * static_cast<double>(index + 1);
    for (dealii::types::global_dof_index index = 0;
         index < control_tangent.size();
         ++index)
      control_tangent[index] =
        (index % 2 == 0 ? 0.03 : -0.02) * static_cast<double>(index + 1);
    const Primal tangent(model.variable_layout(),
                         {std::move(state_tangent), std::move(control_tangent)});
    dealii::Vector<double> seed_values(model.test_layout()->dimension(0));
    for (dealii::types::global_dof_index index = 0;
         index < seed_values.size();
         ++index)
      seed_values[index] = 0.04 * static_cast<double>(index + 1);
    const Primal test_seed(model.test_layout(), {std::move(seed_values)});
    const Covector jvp = model.residual_jvp(evaluation.full_point, tangent);
    const Covector vjp = model.residual_vjp(evaluation.full_point, test_seed);
    require_close(contract::pair(jvp, test_seed),
                  contract::pair(vjp, tangent),
                  1e-11,
                  "Dirichlet-control composed lifting JVP/VJP pairing");

    dealii::Vector<double> control_direction_values(
      control.layout()->dimension(0));
    for (dealii::types::global_dof_index index = 0;
         index < control_direction_values.size();
         ++index)
      control_direction_values[index] =
        (index % 2 == 0 ? 0.05 : -0.04) * static_cast<double>(index + 1);
    const Primal control_direction(control.layout(),
                                   {std::move(control_direction_values)});
    const double directional_derivative =
      contract::pair(evaluation.reduced_derivative, control_direction);
    const auto remainder = [&](const double step) {
      return std::abs(
        reduced.evaluate(shifted(control, control_direction, step)).objective_value -
        evaluation.objective_value - step * directional_derivative);
    };
    const double coarse_remainder = remainder(1e-3);
    const double fine_remainder = remainder(5e-4);
    contract::require(coarse_remainder > 1e-12 &&
                        fine_remainder <= 0.26 * coarse_remainder + 1e-13,
                      "Dirichlet-control reduced Taylor remainder is not quadratic");

    const auto &metric = compilation.problem->metric();
    const Primal metric_direction =
      metric.inverse_apply(evaluation.reduced_derivative);
    require_covector_close(metric.apply(metric_direction),
                           evaluation.reduced_derivative,
                           1e-10,
                           "Dirichlet-control trace metric inverse/apply");
    const auto &manifest = compilation.problem->manifest();
    const auto dirichlet_map = std::find_if(
      manifest.resolved_decision.realized_maps.begin(),
      manifest.resolved_decision.realized_maps.end(),
      [](const compiler::v1::CompiledRealizedMapRecord &map) {
        return map.semantic_id == "dirichlet_control_lifting";
      });
    require_constraint_realisation(manifest, "none", "Dirichlet-control");
    contract::require(
      manifest.compatibility.control_space.find("nodal trace") != std::string::npos &&
        manifest.compatibility.lifting_realisation.find("L_D,h") != std::string::npos &&
        dirichlet_map != manifest.resolved_decision.realized_maps.end() &&
        dirichlet_map->input_space_ids.size() == 2 &&
        dirichlet_map->output_dimension > 0 &&
        manifest.compatibility.metric_solve_policy.find("l2_dirichlet_trace") !=
          std::string::npos &&
        manifest.compatibility.declared_assumptions.front().find("dirichlet_control_lifting") !=
          std::string::npos,
      "Dirichlet-control compilation manifest is incomplete");
  }

  template <int dim>
  void
  run_l2_dirichlet_transposition_lowering_test()
  {
    dealii::Triangulation<dim> triangulation;
    dealii::GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(2);

    const dealii::Functions::ZeroFunction<dim> forcing;
    const dealii::Functions::ConstantFunction<dim> desired_state(0.25);
    const auto specification =
      semantic::v1::make_l2_dirichlet_laplace_control_problem();
    compiler::v1::DealiiDiscretisationPolicy policy;
    policy.state_degree = 1;
    const compiler::v1::DealiiCompiler compiler;
    contract::require(
      compiler.validate(specification, policy).valid(),
      "L2 Dirichlet transposition graph did not validate for deal.II");

    const compiler::v1::DealiiDataBindings<dim> bindings{
      forcing,
      desired_state,
      std::nullopt,
      // This generic aggregate field is intentionally unused because the
      // semantic graph has no reaction-data port.
      7.0,
      0.1,
      test_binding_provenance("l2_dirichlet_transposition")};
    const auto compilation =
      compiler.compile(specification, triangulation, bindings, policy);
    contract::require(
      compilation.succeeded(),
      "L2 Dirichlet transposition graph did not lower through the conforming trace equivalence");

    const auto &model = compilation.problem->executable_model();
    const auto *dirichlet_model = dynamic_cast<
      const compiler::v1::detail::DirichletControlLiftingModel<dim> *>(&model);
    contract::require(
      dirichlet_model != nullptr,
      "L2 Dirichlet transposition compiler did not reuse the lifting target");
    dealii::Vector<double> control_values(model.variable_layout()->dimension(1));
    control_values = 1.0;
    const Primal control(model.variable_layout()->single_block(1, "control"),
                         {std::move(control_values)});
    const auto reduced = compilation.problem->make_reduced_dto();
    const auto evaluation = reduced.evaluate(control);
    require_close(model.residual(evaluation.full_point).block(0).l2_norm(),
                  0.0,
                  1e-11,
                  "L2 Dirichlet transposition-equivalent state residual");
    const auto physical_state =
      dirichlet_model->reconstruct_physical_state(evaluation.full_point);
    for (dealii::types::global_dof_index index = 0;
         index < physical_state.size();
         ++index)
      require_close(physical_state[index],
                    1.0,
                    1e-11,
                    "L2 Dirichlet conforming trace state");

    const Covector conormal = dirichlet_model->discrete_conormal_covector(
      evaluation.full_point, evaluation.adjoint);
    Covector regularisation = compilation.problem->metric().apply(control);
    regularisation.scale_block(0, 0.1);
    const Covector residual_pullback = model.residual_vjp(
      evaluation.full_point, evaluation.adjoint);
    const Covector objective_derivative =
      model.objective_derivative(evaluation.full_point);
    dealii::Vector<double> composed_conormal_values =
      residual_pullback.block(1);
    composed_conormal_values.add(-1.0, objective_derivative.block(1));
    composed_conormal_values.add(1.0, regularisation.block(0));
    const Covector composed_conormal(
      control.layout(), {std::move(composed_conormal_values)});
    require_covector_close(conormal,
                           composed_conormal,
                           1e-11,
                           "L2 Dirichlet discrete conormal pullback");

    Covector expected_stationarity = regularisation;
    expected_stationarity.add_scaled_block(0, -1.0, conormal.block(0));
    require_covector_close(evaluation.reduced_derivative,
                           expected_stationarity,
                           1e-11,
                           "L2 Dirichlet beta M_Gamma u minus conormal sign");
    Covector wrong_plus_stationarity = regularisation;
    wrong_plus_stationarity.add_scaled_block(0, 1.0, conormal.block(0));
    dealii::Vector<double> sign_difference =
      wrong_plus_stationarity.block(0);
    sign_difference.add(-1.0, evaluation.reduced_derivative.block(0));
    contract::require(
      sign_difference.l2_norm() > 1e-6,
      "L2 Dirichlet stationarity test did not distinguish the rejected plus sign");

    dealii::Vector<double> direction_values(control.layout()->dimension(0));
    for (dealii::types::global_dof_index index = 0;
         index < direction_values.size();
         ++index)
      direction_values[index] =
        (index % 2 == 0 ? 0.05 : -0.04) * static_cast<double>(index + 1);
    const Primal direction(control.layout(), {std::move(direction_values)});
    const double derivative =
      contract::pair(evaluation.reduced_derivative, direction);
    const auto remainder = [&](const double step) {
      return std::abs(reduced.evaluate(shifted(control, direction, step))
                        .objective_value -
                      evaluation.objective_value - step * derivative);
    };
    const double coarse_remainder = remainder(1e-3);
    const double fine_remainder = remainder(5e-4);
    contract::require(
      coarse_remainder > 1e-12 &&
        fine_remainder <= 0.26 * coarse_remainder + 1e-13,
      "L2 Dirichlet transposition reduced Taylor remainder is not quadratic");

    const auto &manifest = compilation.problem->manifest();
    test_support::require_dirichlet_manifest_dimensions(
      manifest,
      model.variable_layout()->dimension(0),
      model.test_layout()->dimension(0),
      model.variable_layout()->dimension(1),
      dirichlet_model->physical_state_dimension(),
      "L2 Dirichlet transposition");
    const auto has_assumption = [&manifest](const std::string &prefix) {
      return std::any_of(
        manifest.compatibility.declared_assumptions.begin(),
        manifest.compatibility.declared_assumptions.end(),
        [&prefix](const std::string &assumption) {
          return assumption.find(prefix) == 0;
        });
    };
    const auto has_binding_role = [&manifest](const auto role) {
      return std::any_of(
        manifest.resolved_decision.bindings.begin(),
        manifest.resolved_decision.bindings.end(),
        [role](const compiler::v1::CompiledBindingRecord &binding) {
          return binding.role == role;
        });
    };
    contract::require(
      manifest.compatibility.compiler_id ==
          "nmopt.compiler.v1.dealii.l2_dirichlet_transposition" &&
        manifest.resolved_decision.transposition_realisation.has_value() &&
        manifest.resolved_decision.transposition_realisation->id ==
          "transposition_formulation" &&
        manifest.resolved_decision.transposition_realisation->continuous_parent_space_id ==
          "control_space" &&
        manifest.resolved_decision.transposition_realisation->equivalence_policy_id ==
          "conforming_trace_subspace" &&
        manifest.resolved_decision.transposition_realisation->discrete_realisation ==
          semantic::v1::TranspositionDiscreteRealisation::
            conforming_nodal_lifting_equivalence &&
        manifest.compatibility.state_space.find("continuous L2(Omega) parent") !=
          std::string::npos &&
        manifest.compatibility.control_space.find("U_h=trace(V_h)") != std::string::npos &&
        manifest.compatibility.lifting_realisation.find("E_tr(y,u;f)") !=
          std::string::npos &&
        manifest.compatibility.transformation_ids.empty() &&
        std::find(manifest.compatibility.lowering_handler_records.begin(),
                  manifest.compatibility.lowering_handler_records.end(),
                  "l2_dirichlet_transposition <- "
                  "dealii.dirichlet_control.conforming_trace_equivalence") !=
          manifest.compatibility.lowering_handler_records.end() &&
        has_assumption("transposition_formulation:") &&
        has_assumption("transposition_domain_regularity:") &&
        has_assumption("conforming_trace_subspace:") &&
        has_assumption("discrete_conormal_policy:") &&
        !has_binding_role(semantic::v1::DataRole::diffusion) &&
        !has_binding_role(semantic::v1::DataRole::reaction),
      "L2 Dirichlet manifest omitted its continuous parent, equivalence, or policy provenance");
  }

  template <int dim>
  void
  run_partial_dirichlet_control_contract_test()
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
          cell->face(face)->set_boundary_id(
            cell->face(face)->center()[0] < 1e-12 ? 0 : 1);

    const dealii::Functions::ConstantFunction<dim> forcing(0.5);
    const dealii::Functions::ConstantFunction<dim> desired_state(0.25);
    const dealii::Functions::ConstantFunction<dim> fixed_data(2.0);
    const dealii::Functions::ConstantFunction<dim> changed_fixed_data(3.0);
    const auto specification = semantic::v1::
      make_partial_dirichlet_control_scalar_diffusion_reaction_problem();
    compiler::v1::DealiiDiscretisationPolicy policy;
    policy.state_degree = 1;
    const compiler::v1::DealiiCompiler compiler;
    contract::require(compiler.validate(specification, policy).valid(),
                      "partial Dirichlet-control graph did not validate for deal.II");

    const compiler::v1::DealiiDataBindings<dim> missing_fixed_data{
      forcing, desired_state, 1.0, 0.5, 0.1,
      test_binding_provenance("partial_dirichlet_missing")};
    const auto missing = compiler.compile(specification,
                                          triangulation,
                                          missing_fixed_data,
                                          policy);
    test_support::require_exact_diagnostic(
      missing.diagnostics,
      semantic::v1::DiagnosticCategory::lowerability,
      "state",
      "fixed_dirichlet_data_binding",
      "partial Dirichlet control did not require its fixed lifting data");

    auto bindings = compiler::v1::DealiiDataBindings<dim>{
      forcing, desired_state, 1.0, 0.5, 0.1,
      test_binding_provenance("partial_dirichlet", true)};
    bindings.fixed_dirichlet_data = std::cref(fixed_data);

    auto overlapping_specification = specification;
    component_by_id(overlapping_specification.regions,
                    "fixed_dirichlet_boundary")
      .boundary_ids = {0, 1};
    test_support::require_exact_diagnostic(
      compiler.validate(overlapping_specification, policy),
      semantic::v1::DiagnosticCategory::lowerability,
      "state",
      "partial_dirichlet_boundary_partition",
      "partial Dirichlet control did not reject overlapping boundary regions");

    dealii::Triangulation<dim> missing_control_boundary_triangulation;
    dealii::GridGenerator::hyper_cube(missing_control_boundary_triangulation);
    missing_control_boundary_triangulation.refine_global(1);
    const auto missing_control_boundary = compiler.compile(
      specification,
      missing_control_boundary_triangulation,
      bindings,
      policy);
    test_support::require_exact_diagnostic(
      missing_control_boundary.diagnostics,
      semantic::v1::DiagnosticCategory::lowerability,
      "state",
      "complete_partial_dirichlet_boundary_partition",
      "partial Dirichlet control did not diagnose an absent control boundary");

    const auto compilation = compiler.compile(specification,
                                              triangulation,
                                              bindings,
                                              policy);
    contract::require(compilation.succeeded(),
                      "partial Dirichlet-control compilation failed");
    const auto &model = compilation.problem->executable_model();
    const auto *dirichlet_model =
      dynamic_cast<const compiler::v1::detail::DirichletControlLiftingModel<dim> *>(
        &model);
    contract::require(dirichlet_model != nullptr,
                      "partial Dirichlet compiler did not select its lifting target");
    const auto reduced = compilation.problem->make_reduced_dto();

    dealii::Vector<double> control_values(model.variable_layout()->dimension(1));
    for (dealii::types::global_dof_index index = 0;
         index < control_values.size(); ++index)
      control_values[index] = 1.0;
    const Primal control(model.variable_layout()->single_block(1, "control"),
                         {std::move(control_values)});

    dealii::Vector<double> zero_state(model.variable_layout()->dimension(0));
    const Primal lifting_point(model.variable_layout(),
                               {std::move(zero_state), control.block(0)});
    const auto lifted_trace =
      dirichlet_model->reconstruct_physical_state(lifting_point);
    const auto count_value = [&lifted_trace](const double expected) {
      return std::count_if(
        lifted_trace.begin(), lifted_trace.end(), [expected](const double value) {
          return std::abs(value - expected) <= 1e-12;
        });
    };
    const auto integer_power = [](std::size_t base, unsigned int exponent) {
      std::size_t value = 1;
      for (unsigned int factor = 0; factor < exponent; ++factor)
        value *= base;
      return value;
    };
    const std::size_t nodes_per_axis = 5;
    const std::size_t expected_fixed_dofs =
      integer_power(nodes_per_axis, dim - 1);
    const std::size_t expected_interior_dofs = integer_power(3, dim);
    const std::size_t expected_control_dofs =
      integer_power(nodes_per_axis, dim) - expected_interior_dofs -
      expected_fixed_dofs;
    contract::require(
      count_value(2.0) == expected_fixed_dofs &&
        count_value(1.0) == expected_control_dofs &&
        count_value(0.0) == expected_interior_dofs,
      "partial Dirichlet lifting did not give fixed data precedence at interface DoFs");

    const auto evaluation = reduced.evaluate(control);
    require_close(model.residual(evaluation.full_point).block(0).l2_norm(),
                  0.0,
                  1e-11,
                  "partial Dirichlet-control reconstructed state residual");
    dealii::Vector<double> state_tangent(model.variable_layout()->dimension(0));
    dealii::Vector<double> control_tangent(model.variable_layout()->dimension(1));
    dealii::Vector<double> seed_values(model.test_layout()->dimension(0));
    for (dealii::types::global_dof_index index = 0;
         index < state_tangent.size(); ++index)
      {
        state_tangent[index] = 0.02 * static_cast<double>(index + 1);
        seed_values[index] = -0.03 * static_cast<double>(index + 1);
      }
    for (dealii::types::global_dof_index index = 0;
         index < control_tangent.size(); ++index)
      control_tangent[index] =
        (index % 2 == 0 ? 0.04 : -0.025) * static_cast<double>(index + 1);
    dealii::Vector<double> zero_control_tangent(
      model.variable_layout()->dimension(1));
    dealii::Vector<double> zero_state_tangent(
      model.variable_layout()->dimension(0));
    const Primal state_only_tangent(
      model.variable_layout(), {state_tangent, std::move(zero_control_tangent)});
    const Primal control_only_tangent(
      model.variable_layout(), {std::move(zero_state_tangent), control_tangent});
    const Primal test_seed(model.test_layout(), {std::move(seed_values)});
    const Covector residual_pullback =
      model.residual_vjp(evaluation.full_point, test_seed);
    require_close(contract::pair(model.residual_jvp(evaluation.full_point,
                                                    state_only_tangent),
                                 test_seed),
                  contract::pair(residual_pullback, state_only_tangent),
                  1e-11,
                  "partial Dirichlet-control state reconstruction pullback");
    require_close(contract::pair(model.residual_jvp(evaluation.full_point,
                                                    control_only_tangent),
                                 test_seed),
                  contract::pair(residual_pullback, control_only_tangent),
                  1e-11,
                  "partial Dirichlet-control trace lifting pullback");

    constexpr double derivative_step = 1e-7;
    const Covector objective_pullback =
      model.objective_derivative(evaluation.full_point);
    const auto check_objective_pullback = [&](const Primal &tangent,
                                              const std::string &description) {
      const double centered_difference =
        (model.objective(
           shifted(evaluation.full_point, tangent, derivative_step)) -
         model.objective(
           shifted(evaluation.full_point, tangent, -derivative_step))) /
        (2.0 * derivative_step);
      require_close(centered_difference,
                    contract::pair(objective_pullback, tangent),
                    2e-8,
                    description);
    };
    check_objective_pullback(state_only_tangent,
                             "partial Dirichlet-control state objective pullback");
    check_objective_pullback(control_only_tangent,
                             "partial Dirichlet-control trace objective pullback");

    dealii::Vector<double> direction_values(control.layout()->dimension(0));
    for (dealii::types::global_dof_index index = 0;
         index < direction_values.size(); ++index)
      direction_values[index] =
        (index % 2 == 0 ? 0.05 : -0.04) * static_cast<double>(index + 1);
    const Primal direction(control.layout(), {std::move(direction_values)});
    const double derivative = contract::pair(evaluation.reduced_derivative,
                                             direction);
    const auto remainder = [&](const double step) {
      return std::abs(reduced.evaluate(shifted(control, direction, step))
                        .objective_value - evaluation.objective_value -
                      step * derivative);
    };
    const double coarse_remainder = remainder(1e-3);
    const double fine_remainder = remainder(5e-4);
    contract::require(coarse_remainder > 1e-12 &&
                        fine_remainder <= 0.26 * coarse_remainder + 1e-13,
                      "partial Dirichlet-control reduced Taylor remainder is not quadratic");

    const auto &metric = compilation.problem->metric();
    const Primal metric_direction =
      metric.inverse_apply(evaluation.reduced_derivative);
    require_covector_close(metric.apply(metric_direction),
                           evaluation.reduced_derivative,
                           1e-10,
                           "partial Dirichlet-control trace metric inverse/apply");

    auto changed_bindings = compiler::v1::DealiiDataBindings<dim>{
      forcing, desired_state, 1.0, 0.5, 0.1,
      test_binding_provenance("partial_dirichlet_changed", true)};
    changed_bindings.fixed_dirichlet_data = std::cref(changed_fixed_data);
    const auto changed = compiler.compile(specification,
                                          triangulation,
                                          changed_bindings,
                                          policy);
    contract::require(changed.succeeded() &&
                        std::abs(changed.problem->make_reduced_dto()
                                   .evaluate(control)
                                   .objective_value -
                                 evaluation.objective_value) > 1e-4,
                      "partial Dirichlet-control recompilation reused fixed data");
    const auto *changed_dirichlet_model = dynamic_cast<
      const compiler::v1::detail::DirichletControlLiftingModel<dim> *>(
        &changed.problem->executable_model());
    contract::require(changed_dirichlet_model != nullptr,
                      "changed partial Dirichlet compilation lost its lifting target");
    const auto changed_trace =
      changed_dirichlet_model->reconstruct_physical_state(lifting_point);
    contract::require(
      std::count_if(changed_trace.begin(),
                    changed_trace.end(),
                    [](const double value) {
                      return std::abs(value - 3.0) <= 1e-12;
                    }) == expected_fixed_dofs,
      "changed partial Dirichlet data did not own the interface DoFs");

    const auto &manifest = compilation.problem->manifest();
    test_support::require_dirichlet_manifest_dimensions(
      manifest,
      model.variable_layout()->dimension(0),
      model.test_layout()->dimension(0),
      model.variable_layout()->dimension(1),
      dirichlet_model->physical_state_dimension(),
      "partial Dirichlet-control");
    contract::require(
      manifest.resolved_decision.partial_boundary_selection.has_value() &&
        manifest.resolved_decision.partial_boundary_selection->fixed_boundary_region_id ==
          "fixed_dirichlet_boundary" &&
        manifest.resolved_decision.partial_boundary_selection->controlled_boundary_region_id ==
          "control_boundary" &&
        manifest.resolved_decision.partial_boundary_selection->interface_realisation ==
          semantic::v1::PartialDirichletInterfaceRealisation::
            fixed_data_precedence &&
        manifest.resolved_decision.partial_boundary_selection->trace_realisation ==
          semantic::v1::PartialDirichletTraceRealisation::
            relative_interior_nodal_zero_endpoint &&
      manifest.compatibility.lifting_realisation.find("ell_0,h + L_D,h") != std::string::npos &&
        manifest.compatibility.data_rule.find("fixed Dirichlet Function") != std::string::npos &&
        std::any_of(manifest.compatibility.declared_assumptions.begin(),
                    manifest.compatibility.declared_assumptions.end(),
                    [](const std::string &assumption) {
                      return assumption.find("fixed-data precedence") !=
                             std::string::npos;
                    }),
      "partial Dirichlet-control manifest omitted the lifting interface policy");
  }

