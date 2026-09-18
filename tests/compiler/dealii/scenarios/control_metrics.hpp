#pragma once

// These scenario implementations are deliberately included into one heavy
// deal.II/compiler translation unit to improve source navigation without
// multiplying compilation cost.

  template <int dim>
  void
  run_continuous_control_component_contract_test()
  {
    dealii::Triangulation<dim> triangulation;
    dealii::GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(1);
    const dealii::Functions::ConstantFunction<dim> forcing(0.0);
    const EnergyPolynomial<dim> desired_state;
    const compiler::v1::detail::ContinuousControlModel<dim> model(
      triangulation,
      forcing,
      desired_state,
      1.0,
      0.5,
      0.2,
      1,
      {0},
      true,
      false,
      true);

    contract::require(model.variable_layout()->dimension(1) == 1,
                      "homogeneous-Dirichlet Q1 control did not expose only its independent DoF");
    dealii::Vector<double> state(model.variable_layout()->dimension(0));
    dealii::Vector<double> control_values(1);
    control_values[0] = 0.4;
    const Primal point(model.variable_layout(), {std::move(state), control_values});
    const Covector derivative = model.objective_derivative(point);
    require_close(derivative.block(1)[0],
                  0.2 * 0.4 / 9.0,
                  1e-13,
                  "continuous-control L2 regularisation on independent coordinates");

    const auto control_layout =
      model.variable_layout()->single_block(1, "control");
    const Primal control(control_layout, {std::move(control_values)});
    const auto metric = model.control_hminus1_metric();
    const Covector metric_covector = metric.apply(control);
    require_primal_close(metric.inverse_apply(metric_covector),
                         control,
                         1e-11,
                         "continuous-control H-1 metric round trip");
  }

  template <int dim>
  void
  run_l2_tracking_continuous_control_contract_test()
  {
    dealii::Triangulation<dim> triangulation;
    dealii::GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(1);
    const dealii::Functions::ConstantFunction<dim> forcing(0.0);
    const EnergyPolynomial<dim> desired_state;
    const auto specification =
      semantic::v1::make_l2_state_tracking_continuous_control_problem();
    compiler::v1::DealiiDiscretisationPolicy policy;
    policy.state_degree = 1;
    const compiler::v1::DealiiCompiler compiler;
    contract::require(
      compiler.validate(specification, policy).valid(),
      "L2-tracking continuous-control graph did not validate for deal.II");

    auto non_full_domain_specification = specification;
    component_by_id(non_full_domain_specification.regions, "domain")
      .is_full_domain = false;
    test_support::require_exact_diagnostic(
      compiler.validate(non_full_domain_specification, policy),
      semantic::v1::DiagnosticCategory::lowerability,
      non_full_domain_specification.id,
      "continuous_control_full_domain_tracking",
      "Continuous-control compilation accepted a non-full-domain observation");

    const compiler::v1::DealiiDataBindings<dim> bindings{
      forcing,
      desired_state,
      1.0,
      0.0,
      0.2,
      test_binding_provenance("l2_tracking_continuous_control")};
    const auto compilation = compiler.compile(specification,
                                              triangulation,
                                              bindings,
                                              policy);
    contract::require(compilation.succeeded() && compilation.problem,
                      "L2-tracking continuous-control compilation failed");
    contract::require(
      compilation.problem->metric().id() == "l2_continuous" &&
        compilation.problem->manifest().compatibility.control_space.find(
          "homogeneous-Dirichlet scalar FE_Q(1)") != std::string::npos,
      "L2-tracking continuous-control compilation selected the wrong control realization");
    require_compiled_hessian_evidence(
      *compilation.problem, "L2-tracking continuous control");

    const auto &model = compilation.problem->executable_model();
    const auto reduced = compilation.problem->make_reduced_dto();
    dealii::Vector<double> control_values(model.variable_layout()->dimension(1));
    const Primal control(model.variable_layout()->single_block(1, "control"),
                         {std::move(control_values)});
    const auto evaluation = reduced.evaluate(control);
    const auto *continuous_model = dynamic_cast<const
      compiler::v1::detail::ContinuousControlModel<dim> *>(&model);
    contract::require(continuous_model != nullptr,
                      "L2-tracking compilation did not retain its continuous model");
    const auto output_directory =
      std::filesystem::temp_directory_path() /
      "nmopt-l2-tracking-continuous-control-contract";
    std::filesystem::remove_all(output_directory);
    continuous_model->write_native_output(output_directory,
                                          evaluation.state,
                                          control,
                                          evaluation.adjoint,
                                          &forcing,
                                          &desired_state);
    std::ifstream output(output_directory / "fields-volume.vtu");
    const std::string document((std::istreambuf_iterator<char>(output)),
                               std::istreambuf_iterator<char>());
    contract::require(document.find("Name=\"control\"") != std::string::npos &&
                        document.find("Name=\"state\"") != std::string::npos &&
                        document.find("Name=\"adjoint\"") != std::string::npos,
                      "Continuous-control native output omitted a primary field");
    std::filesystem::remove_all(output_directory);
  }

  template <int dim>
  void
  run_simplex_continuous_control_contract_test()
  {
    dealii::Triangulation<dim> triangulation;
    dealii::GridGenerator::subdivided_hyper_cube_with_simplices(triangulation,
                                                                 2);
    const dealii::Functions::ConstantFunction<dim> forcing(0.0);
    const EnergyPolynomial<dim> desired_state;
    const auto specification =
      semantic::v1::make_l2_state_tracking_continuous_control_problem();
    compiler::v1::DealiiDiscretisationPolicy policy;
    policy.state_degree = 1;
    const compiler::v1::DealiiCompiler compiler;
    const compiler::v1::DealiiDataBindings<dim> bindings{
      forcing,
      desired_state,
      1.0,
      0.0,
      0.2,
      test_binding_provenance("simplex_continuous_control")};

    const auto compilation = compiler.compile(specification,
                                              triangulation,
                                              bindings,
                                              policy);
    contract::require(compilation.succeeded() && compilation.problem,
                      "Continuous-control compilation rejected a simplex mesh");
    const auto &manifest = compilation.problem->manifest();
    const auto state_observation = std::find_if(
      manifest.resolved_decision.realized_maps.begin(),
      manifest.resolved_decision.realized_maps.end(),
      [](const compiler::v1::CompiledRealizedMapRecord &map) {
        return map.semantic_id == "state_observation";
      });
    contract::require(
      compilation.problem->executable_model().variable_layout()->dimension(0) ==
          9 &&
        compilation.problem->metric().layout()->dimension(0) == 1 &&
        manifest.compatibility.state_space == "scalar FE_SimplexP(1)" &&
        manifest.compatibility.control_space.find("FE_SimplexP(1)") != std::string::npos &&
        manifest.compatibility.quadrature == "QGaussSimplex(3)" &&
        manifest.compatibility.data_rule.find("QGaussSimplex(3) volume quadrature") !=
          std::string::npos &&
        state_observation != manifest.resolved_decision.realized_maps.end() &&
        state_observation->realization_id ==
          "fe_simplex_p_coefficient_restriction" &&
        state_observation->output_layout.find("FE_SimplexP") !=
          std::string::npos,
      "Simplex compilation did not report its P1 element, quadrature, or independent control layout");
    require_compiled_hessian_evidence(*compilation.problem,
                                      "simplex continuous control");

    const auto unsupported = compiler.compile(
      semantic::v1::make_scalar_diffusion_reaction_problem(),
      triangulation,
      bindings,
      policy);
    contract::require(!unsupported.succeeded(),
                      "A hypercube-only compiler target accepted simplex cells");
    test_support::require_exact_diagnostic(
      unsupported.diagnostics,
      semantic::v1::DiagnosticCategory::lowerability,
      "scalar_diffusion_reaction_volume_control",
      "simplex_registered_target",
      "Simplex target validation did not identify the bounded registered-target capability");
  }

  template <int dim>
  void
  run_continuous_neumann_control_lowering_contract_test()
  {
    const auto specification =
      semantic::v1::make_neumann_boundary_control_problem(
        false,
        semantic::v1::NeumannControlDiscretisation::continuous_nodal_trace);
    compiler::v1::DealiiDiscretisationPolicy policy;
    policy.state_degree = 1;
    const compiler::v1::DealiiCompiler compiler;
    contract::require(
      compiler.validate(specification, policy).valid(),
      "Continuous Neumann boundary-control graph did not validate for deal.II");

    for (unsigned int family = 0; family < 2; ++family)
      {
        const bool simplex = family == 1;
        dealii::Triangulation<dim> triangulation;
        if (simplex)
          dealii::GridGenerator::subdivided_hyper_cube_with_simplices(
            triangulation,
            2);
        else
          dealii::GridGenerator::subdivided_hyper_cube(triangulation, 2);
        for (auto cell = triangulation.begin_active();
             cell != triangulation.end();
             ++cell)
          for (unsigned int face = 0; face < cell->n_faces(); ++face)
            if (cell->face(face)->at_boundary())
              {
                const double x = cell->face(face)->center()[0];
                cell->face(face)->set_boundary_id(x < 1e-12 ? 0 :
                                                  x > 1.0 - 1e-12 ? 1 : 2);
              }

        const dealii::Functions::ConstantFunction<dim> forcing(0.5);
        const dealii::Functions::ConstantFunction<dim> desired_state(0.2);
        const compiler::v1::DealiiDataBindings<dim> bindings{
          forcing,
          desired_state,
          1.0,
          0.5,
          0.1,
          test_binding_provenance(
            simplex ? "continuous_neumann_simplex" :
                      "continuous_neumann_hypercube")};
        const auto compilation = compiler.compile(specification,
                                                  triangulation,
                                                  bindings,
                                                  policy);
        contract::require(compilation.succeeded(),
                          "Continuous Neumann boundary-control compilation failed");
        contract::require(
          compilation.problem->constraint() == nullptr,
          "Continuous Neumann boundary control unexpectedly produced a box");

        const auto &model = compilation.problem->executable_model();
        const auto *neumann_model = dynamic_cast<const
          compiler::v1::detail::NeumannBoundaryControlModel<dim> *>(&model);
        contract::require(
          neumann_model != nullptr &&
            neumann_model->physical_control_dimension() == 3 &&
            neumann_model->independent_control_dimension() == 3 &&
            neumann_model->control_coordinates().size() == 3,
          "Continuous Neumann lowering changed the trace dimensions");

        const auto &metric = compilation.problem->metric();
        contract::require(
          metric.id() == "l2_neumann_trace" &&
            metric.layout()->dimension(0) == 3,
          "Continuous Neumann lowering selected the wrong metric realization");
        const auto reduced = compilation.problem->make_reduced_dto();
        dealii::Vector<double> control_values(3);
        control_values[0] = 0.06;
        control_values[1] = -0.04;
        control_values[2] = 0.03;
        const Primal control(model.variable_layout()->single_block(1, "control"),
                             {std::move(control_values)});
        const auto evaluation = reduced.evaluate(control);
        contract::require(
          evaluation.state_solve.converged() &&
            evaluation.adjoint_solve.converged(),
          "Continuous Neumann lowering did not converge in state and adjoint solves");
        require_close(model.residual(evaluation.full_point).block(0).l2_norm(),
                      0.0,
                      2e-11,
                      "Continuous Neumann boundary-control state residual");

        const auto output_directory =
          std::filesystem::temp_directory_path() /
          (simplex ? "nmopt-continuous-neumann-simplex-output" :
                     "nmopt-continuous-neumann-hypercube-output");
        std::filesystem::remove_all(output_directory);
        neumann_model->write_native_output(output_directory,
                                           evaluation.state,
                                           control,
                                           evaluation.adjoint,
                                           nullptr,
                                           &forcing,
                                           &desired_state);
        contract::require(
          std::filesystem::exists(output_directory / "mesh-volume.vtu") &&
            std::filesystem::exists(output_directory / "mesh-volume.svg") &&
            std::filesystem::exists(output_directory / "fields-volume.vtu") &&
            std::filesystem::exists(output_directory /
                                    "control-boundary.vtu"),
          "Continuous Neumann model output omitted a native file");
        std::filesystem::remove_all(output_directory);

        const Covector measured_control = metric.apply(control);
        require_close(
          neumann_model->objective_components(evaluation.full_point)
            .control_regularisation,
          0.5 * bindings.regularisation_weight *
            contract::pair(measured_control, control),
          1e-14,
          "Continuous Neumann control regularisation objective");
        dealii::Vector<double> regularisation_error =
          model.objective_derivative(evaluation.full_point).block(1);
        regularisation_error.add(-bindings.regularisation_weight,
                                 measured_control.block(0));
        require_close(regularisation_error.l2_norm(),
                      0.0,
                      1e-14,
                      "Continuous Neumann control regularisation derivative");

        dealii::Vector<double> state_tangent(
          model.variable_layout()->dimension(0));
        dealii::Vector<double> control_tangent(3);
        control_tangent[0] = 0.03;
        control_tangent[1] = -0.02;
        control_tangent[2] = 0.04;
        const Primal coupling_tangent(model.variable_layout(),
                                      {std::move(state_tangent),
                                       std::move(control_tangent)});
        dealii::Vector<double> seed_values(model.test_layout()->dimension(0));
        for (dealii::types::global_dof_index index = 0;
             index < seed_values.size();
             ++index)
          seed_values[index] =
            (index % 2 == 0 ? 0.025 : -0.015) *
            static_cast<double>(index + 1);
        const Primal test_seed(model.test_layout(), {std::move(seed_values)});
        const Covector coupling_jvp =
          model.residual_jvp(evaluation.full_point, coupling_tangent);
        contract::require(coupling_jvp.block(0).l2_norm() > 1e-12,
                          "Continuous Neumann coupling action vanished");
        require_close(
          contract::pair(coupling_jvp, test_seed),
          contract::pair(model.residual_vjp(evaluation.full_point, test_seed),
                         coupling_tangent),
          2e-11,
          "Continuous Neumann coupling JVP/VJP pairing");

        constexpr double derivative_step = 1e-7;
        Covector residual_difference = model.residual(
          shifted(evaluation.full_point, coupling_tangent, derivative_step));
        const Covector base_residual = model.residual(evaluation.full_point);
        residual_difference.add_scaled_block(0,
                                             -1.0,
                                             base_residual.block(0));
        residual_difference.scale_block(0, 1.0 / derivative_step);
        residual_difference.add_scaled_block(0,
                                             -1.0,
                                             coupling_jvp.block(0));
        require_close(residual_difference.block(0).l2_norm(),
                      0.0,
                      2e-7,
                      "Continuous Neumann residual finite-difference action");

        dealii::Vector<double> direction_values(3);
        direction_values[0] = -0.05;
        direction_values[1] = 0.04;
        direction_values[2] = 0.02;
        const Primal direction(control.layout(),
                               {std::move(direction_values)});
        const double directional_derivative =
          contract::pair(evaluation.reduced_derivative, direction);
        const double central_difference =
          (reduced.evaluate(shifted(control, direction, derivative_step))
             .objective_value -
           reduced.evaluate(shifted(control, direction, -derivative_step))
             .objective_value) /
          (2.0 * derivative_step);
        require_close(central_difference,
                      directional_derivative,
                      2e-8,
                      "Continuous Neumann reduced-gradient finite difference");
        const auto remainder = [&](const double step) {
          return std::abs(
            reduced.evaluate(shifted(control, direction, step)).objective_value -
            evaluation.objective_value - step * directional_derivative);
        };
        const double coarse_remainder = remainder(1e-3);
        const double fine_remainder = remainder(5e-4);
        contract::require(
          coarse_remainder > 1e-12 &&
            fine_remainder <= 0.26 * coarse_remainder + 1e-13,
          "Continuous Neumann reduced Taylor remainder is not quadratic");

        const Primal metric_direction =
          reduced.gradient_direction(evaluation.reduced_derivative, metric);
        require_close(contract::pair(metric.apply(metric_direction), direction),
                      directional_derivative,
                      2e-11,
                      "Continuous Neumann metric-gradient identity");

        const auto &manifest = compilation.problem->manifest();
        const bool has_control_space = std::any_of(
          manifest.resolved_decision.spaces.begin(),
          manifest.resolved_decision.spaces.end(),
          [](const compiler::v1::CompiledSpaceRecord &space) {
            return space.role == semantic::v1::SpaceRole::control &&
                   space.dimension == 3 &&
                   space.finite_element.find(
                     "continuous scalar degree-one nodal trace") !=
                     std::string::npos;
          });
        contract::require(
          has_control_space &&
            manifest.resolved_decision.metric_record.realisation_id == "l2_neumann_trace" &&
            std::find(manifest.compatibility.lowering_handler_records.begin(),
                      manifest.compatibility.lowering_handler_records.end(),
                      "neumann_control <- dealii.neumann.control.continuous_p1_trace") !=
              manifest.compatibility.lowering_handler_records.end() &&
            std::any_of(
              manifest.compatibility.declared_assumptions.begin(),
              manifest.compatibility.declared_assumptions.end(),
              [](const std::string &assumption) {
                return assumption.find("closure endpoints") !=
                       std::string::npos;
              }),
          "Continuous Neumann compilation manifest is incomplete");
      }
  }

  template <int dim>
  void
  run_h1_control_regularisation_contract_test()
  {
    dealii::Triangulation<dim> triangulation;
    dealii::GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(1);

    const dealii::Functions::ConstantFunction<dim> forcing(0.0);
    const dealii::Functions::ConstantFunction<dim> desired_state(0.0);
    const auto specification =
      semantic::v1::make_h1_regularised_scalar_diffusion_reaction_problem();
    const auto h1_metric_specification =
      semantic::v1::make_h1_metric_scalar_diffusion_reaction_problem();
    compiler::v1::DealiiDiscretisationPolicy policy;
    policy.state_degree = 1;
    const compiler::v1::DealiiCompiler compiler;
    contract::require(
      compiler.validate(specification, policy).valid(),
      "H1-control regularisation v1 graph did not validate for deal.II");
    contract::require(
      compiler.validate(h1_metric_specification, policy).valid(),
      "H1-control metric v1 graph did not validate for deal.II");

    auto unsupported_h1_metric =
      semantic::v1::make_scalar_diffusion_reaction_problem();
    unsupported_h1_metric.metrics.at(0) =
      {"control_h1_metric", "Unsupported discontinuous H1 metric",
       semantic::v1::MetricKind::h1, "control", "control_pairing"};
    unsupported_h1_metric.formulation.metric_id = "control_h1_metric";
    const auto unsupported_h1_metric_report =
      compiler.validate(unsupported_h1_metric, policy);
    test_support::require_exact_diagnostic(
      unsupported_h1_metric_report,
      semantic::v1::DiagnosticCategory::lowerability,
      "control_h1_metric",
      "h1_metric_registered_control_space",
      "H1 metric compiler did not require the continuous H1-control target");

    auto discontinuous_control = specification;
    discontinuous_control.spaces.at(2).topology =
      semantic::v1::SpaceTopology::l2;
    const auto discontinuous_control_report =
      compiler.validate(discontinuous_control, policy);
    test_support::require_exact_diagnostic(
      discontinuous_control_report,
      semantic::v1::DiagnosticCategory::lowerability,
      "control",
      "h1_continuous_control_space",
      "H1-control compiler did not require the continuous control realization");

    auto h1_control_box = specification;
    const auto cellwise_box_source =
      semantic::v1::make_scalar_diffusion_reaction_problem(true);
    h1_control_box.data.push_back(cellwise_box_source.data.at(5));
    h1_control_box.data.push_back(cellwise_box_source.data.at(6));
    h1_control_box.constraints = cellwise_box_source.constraints;
    h1_control_box.requirement_policies.push_back(
      cellwise_box_source.requirement_policies.at(2));
    h1_control_box.formulation.constraint_id = "control_box";
    const auto h1_control_box_report = compiler.validate(h1_control_box, policy);
    test_support::require_exact_diagnostic(
      h1_control_box_report,
      semantic::v1::DiagnosticCategory::lowerability,
      "control_box",
      "continuous_control_box_constraint",
      "H1-control compiler did not reject the unsupported cellwise box");

    const compiler::v1::DealiiDataBindings<dim> bindings{
      forcing,
      desired_state,
      1.0,
      0.5,
      0.2,
      test_binding_provenance("h1_control")};
    const auto compilation = compiler.compile(specification,
                                              triangulation,
                                              bindings,
                                              policy);
    contract::require(compilation.succeeded(),
                      "H1-control regularisation v1 compilation failed");
    const auto h1_metric_compilation = compiler.compile(h1_metric_specification,
                                                        triangulation,
                                                        bindings,
                                                        policy);
    contract::require(h1_metric_compilation.succeeded(),
                      "H1-control metric v1 compilation failed");
    require_compiled_hessian_evidence(*compilation.problem,
                                      "H1-control L2 metric");
    require_compiled_hessian_evidence(*h1_metric_compilation.problem,
                                      "H1-control H1 metric");

    const auto &model = compilation.problem->executable_model();
    const auto reduced = compilation.problem->make_reduced_dto();
    dealii::Vector<double> control_values(model.variable_layout()->dimension(1));
    for (dealii::types::global_dof_index index = 0;
         index < control_values.size();
         ++index)
      control_values[index] = index % 2 == 0
                                ? 0.1 * static_cast<double>(index + 1)
                                : -0.05 * static_cast<double>(index + 1);
    const Primal control(model.variable_layout()->single_block(1, "control"),
                         {control_values});

    dealii::Vector<double> zero_state(model.variable_layout()->dimension(0));
    const Primal direct_objective_point(model.variable_layout(),
                                        {std::move(zero_state), control_values});
    const Covector direct_objective_derivative =
      model.objective_derivative(direct_objective_point);
    const Covector direct_control_derivative = contract::extract_covector_block(
      direct_objective_derivative, 1, "control");
    const Primal l2_direction = compilation.problem->metric().inverse_apply(
      direct_control_derivative);
    dealii::Vector<double> stiffness_component = l2_direction.block(0);
    stiffness_component.add(-0.2, control.block(0));
    contract::require(stiffness_component.l2_norm() > 1e-4,
                      "H1 regularisation did not contribute its control stiffness");

    const auto &h1_metric_model =
      h1_metric_compilation.problem->executable_model();
    const auto h1_metric_reduced =
      h1_metric_compilation.problem->make_reduced_dto();
    const Primal h1_metric_control(
      h1_metric_model.variable_layout()->single_block(1, "control"),
      {control_values});
    dealii::Vector<double> h1_metric_zero_state(
      h1_metric_model.variable_layout()->dimension(0));
    const Primal h1_metric_direct_objective_point(
      h1_metric_model.variable_layout(),
      {std::move(h1_metric_zero_state), control_values});
    const Covector h1_metric_direct_objective_derivative =
      h1_metric_model.objective_derivative(h1_metric_direct_objective_point);
    const Covector h1_metric_control_derivative =
      contract::extract_covector_block(h1_metric_direct_objective_derivative,
                                       1,
                                       "control");
    const Primal h1_metric_direction =
      h1_metric_compilation.problem->metric().inverse_apply(
        h1_metric_control_derivative);
    const Primal expected_h1_metric_direction(
      h1_metric_control.layout(), {control_values});
    dealii::Vector<double> h1_metric_difference = h1_metric_direction.block(0);
    h1_metric_difference.add(-0.2, expected_h1_metric_direction.block(0));
    require_close(h1_metric_difference.l2_norm(),
                  0.0,
                  1e-10,
                  "H1 metric did not invert the H1 regularisation Riesz map");
    require_covector_close(
      h1_metric_compilation.problem->metric().apply(h1_metric_direction),
      h1_metric_control_derivative,
      1e-10,
      "H1 metric inverse/apply relation");

    const auto evaluation = reduced.evaluate(control);
    const auto h1_metric_evaluation = h1_metric_reduced.evaluate(h1_metric_control);
    require_close(h1_metric_evaluation.objective_value,
                  evaluation.objective_value,
                  1e-12,
                  "H1 metric changed the reduced objective");
    require_covector_close(h1_metric_evaluation.reduced_derivative,
                           evaluation.reduced_derivative,
                           1e-11,
                           "H1 metric changed the reduced derivative");
    require_close(model.residual(evaluation.full_point).block(0).l2_norm(),
                  0.0,
                  1e-11,
                  "H1-control regularisation state residual");

    dealii::Vector<double> state_tangent(model.variable_layout()->dimension(0));
    dealii::Vector<double> control_tangent(model.variable_layout()->dimension(1));
    for (dealii::types::global_dof_index index = 0;
         index < state_tangent.size();
         ++index)
      state_tangent[index] = -0.01 * static_cast<double>(index + 1);
    for (dealii::types::global_dof_index index = 0;
         index < control_tangent.size();
         ++index)
      control_tangent[index] = index % 2 == 0
                                  ? 0.03 * static_cast<double>(index + 1)
                                  : -0.02 * static_cast<double>(index + 1);
    const Primal tangent(model.variable_layout(),
                         {std::move(state_tangent), std::move(control_tangent)});
    dealii::Vector<double> test_seed_values(model.test_layout()->dimension(0));
    for (dealii::types::global_dof_index index = 0;
         index < test_seed_values.size();
         ++index)
      test_seed_values[index] = 0.04 * static_cast<double>(index + 1);
    const Primal test_seed(model.test_layout(), {std::move(test_seed_values)});
    require_close(
      contract::pair(model.residual_jvp(evaluation.full_point, tangent), test_seed),
      contract::pair(model.residual_vjp(evaluation.full_point, test_seed), tangent),
      1e-11,
      "H1-control regularisation residual JVP/VJP pairing");

    dealii::Vector<double> direction_values(control.layout()->dimension(0));
    for (dealii::types::global_dof_index index = 0;
         index < direction_values.size();
         ++index)
      direction_values[index] = index % 2 == 0 ? 0.02 : -0.03;
    const Primal direction(control.layout(), {std::move(direction_values)});
    const double directional_derivative =
      contract::pair(evaluation.reduced_derivative, direction);
    const auto remainder = [&](const double step) {
      return std::abs(reduced.evaluate(shifted(control, direction, step)).objective_value -
                      evaluation.objective_value - step * directional_derivative);
    };
    const double coarse_remainder = remainder(1e-3);
    const double fine_remainder = remainder(5e-4);
    contract::require(coarse_remainder > 1e-12 &&
                        fine_remainder <= 0.26 * coarse_remainder + 1e-13,
                      "H1-control reduced Taylor remainder is not quadratic");

    const auto &metric = compilation.problem->metric();
    contract::require(metric.id() == "l2_continuous",
                      "H1 regularisation incorrectly selected an H1 search metric");
    contract::require(h1_metric_compilation.problem->metric().id() ==
                        "h1_continuous",
                      "H1 metric compilation did not select the H1 Riesz map");
    const auto &manifest = compilation.problem->manifest();
    require_constraint_realisation(manifest, "none", "H1-control L2 metric");
    require_constraint_realisation(h1_metric_compilation.problem->manifest(),
                                   "none",
                                   "H1-control H1 metric");
    contract::require(
      manifest.compatibility.control_space.find("continuous scalar FE_Q") != std::string::npos &&
        manifest.compatibility.declared_assumptions.front().find("h1_control_regularisation") !=
          std::string::npos,
      "H1-control compilation manifest omitted the loss or control realization");
    contract::require(
      h1_metric_compilation.problem->manifest().compatibility.metric_solve_policy.find(
        "h1_continuous") != std::string::npos,
      "H1 metric compilation manifest omitted the selected Riesz map");
  }

  template <int dim>
  void
  run_hminus1_compilation_contract_test()
  {
    dealii::Triangulation<dim> triangulation;
    dealii::GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(2);

    const dealii::Functions::ConstantFunction<dim> forcing(1.0);
    const EnergyPolynomial<dim> desired_state(0.5);
    auto hminus1_specification = semantic::v1::
      make_hminus1_metric_h1_state_tracking_scalar_diffusion_reaction_problem();
    auto l2_specification = semantic::v1::
      make_l2_metric_h1_state_tracking_continuous_control_problem();
    compiler::v1::DealiiDiscretisationPolicy policy;
    policy.state_degree = 1;
    policy.control_metric_solve.maximum_iterations = 500;
    policy.control_metric_solve.relative_tolerance = 1e-13;
    policy.control_metric_solve.absolute_tolerance = 1e-15;
    const compiler::v1::DealiiCompiler compiler;
    contract::require(compiler.validate(hminus1_specification, policy).valid() &&
                        compiler.validate(l2_specification, policy).valid(),
                      "H-1/L2 continuous-control comparison graphs did not validate");

    auto missing_energy_observation = hminus1_specification;
    auto &state_observation = component_by_id(
      missing_energy_observation.observations, "state_observation");
    state_observation.kind = semantic::v1::ObservationKind::volume_restriction;
    component_by_id(missing_energy_observation.spaces,
                    "state_observation_space")
      .topology = semantic::v1::SpaceTopology::l2;
    test_support::require_exact_diagnostic(
      compiler.validate(missing_energy_observation, policy),
      semantic::v1::DiagnosticCategory::lowerability,
      missing_energy_observation.id,
      "hminus1_metric_energy_observation",
      "H-1 compiler did not require the P5.2 energy observation");

    const compiler::v1::DealiiDataBindings<dim> bindings{
      forcing,
      desired_state,
      1.0,
      0.5,
      0.2,
      test_binding_provenance("hminus1_metric")};
    const auto hminus1_compilation = compiler.compile(hminus1_specification,
                                                      triangulation,
                                                      bindings,
                                                      policy);
    const auto l2_compilation = compiler.compile(l2_specification,
                                                 triangulation,
                                                 bindings,
                                                 policy);
    contract::require(hminus1_compilation.succeeded() &&
                        l2_compilation.succeeded(),
                      "H-1/L2 continuous-control comparison compilation failed");

    dealii::Triangulation<dim> incomplete_boundary_triangulation;
    dealii::GridGenerator::hyper_cube(incomplete_boundary_triangulation);
    incomplete_boundary_triangulation.refine_global(1);
    for (auto cell = incomplete_boundary_triangulation.begin_active();
         cell != incomplete_boundary_triangulation.end();
         ++cell)
      for (unsigned int face = 0;
           face < dealii::GeometryInfo<dim>::faces_per_cell;
           ++face)
        if (cell->face(face)->at_boundary())
          cell->face(face)->set_boundary_id(
            cell->face(face)->center()[0] < 0.5 ? 0 : 1);
    const auto incomplete_boundary = compiler.compile(
      hminus1_specification,
      incomplete_boundary_triangulation,
      bindings,
      policy);
    test_support::require_exact_diagnostic(
      incomplete_boundary.diagnostics,
      semantic::v1::DiagnosticCategory::lowerability,
      "control",
      "continuous_control_complete_boundary",
      "H-1 compilation accepted a continuous-control boundary that did not cover the exterior mesh");

    const auto &hminus1_model =
      hminus1_compilation.problem->executable_model();
    const auto &l2_model = l2_compilation.problem->executable_model();
    contract::require(
      hminus1_model.variable_layout()->dimension(1) ==
          l2_model.variable_layout()->dimension(1) &&
        hminus1_model.variable_layout()->dimension(1) > 1,
      "H-1/L2 comparison did not retain one independent continuous-control layout");
    dealii::Vector<double> control_values(
      hminus1_model.variable_layout()->dimension(1));
    for (dealii::types::global_dof_index index = 0;
         index < control_values.size();
         ++index)
      control_values[index] =
        (index % 2 == 0 ? 0.03 : -0.02) * static_cast<double>(index + 1);
    const Primal hminus1_control(
      hminus1_model.variable_layout()->single_block(1, "control"),
      {control_values});
    const Primal l2_control(
      l2_model.variable_layout()->single_block(1, "control"), {control_values});
    const auto hminus1_reduced =
      hminus1_compilation.problem->make_reduced_dto();
    const auto l2_reduced = l2_compilation.problem->make_reduced_dto();
    const auto hminus1_evaluation = hminus1_reduced.evaluate(hminus1_control);
    const auto l2_evaluation = l2_reduced.evaluate(l2_control);

    require_close(hminus1_evaluation.objective_value,
                  l2_evaluation.objective_value,
                  1e-12,
                  "H-1 metric changed the reduced objective");
    require_covector_close(hminus1_evaluation.reduced_derivative,
                           l2_evaluation.reduced_derivative,
                           1e-11,
                           "H-1 metric changed the reduced covector");
    const Primal hminus1_direction =
      hminus1_compilation.problem->metric().inverse_apply(
        hminus1_evaluation.reduced_derivative);
    const Primal l2_direction = l2_compilation.problem->metric().inverse_apply(
      l2_evaluation.reduced_derivative);
    dealii::Vector<double> direction_difference = hminus1_direction.block(0);
    direction_difference.add(-1.0, l2_direction.block(0));
    contract::require(direction_difference.l2_norm() > 1e-3,
                      "H-1 and L2 metrics produced the same search direction");
    require_covector_close(
      hminus1_compilation.problem->metric().apply(hminus1_direction),
      hminus1_evaluation.reduced_derivative,
      1e-10,
      "compiled H-1 metric apply/inverse pairing");

    const auto verify_reduced_taylor = [](const auto &reduced,
                                          const Primal &control,
                                          const auto &evaluation,
                                          const Primal &direction,
                                          const char *description) {
      const double slope =
        contract::pair(evaluation.reduced_derivative, direction);
      const auto remainder = [&](const double step) {
        return std::abs(
          reduced.evaluate(shifted(control, direction, step)).objective_value -
          evaluation.objective_value - step * slope);
      };
      const double coarse = remainder(2e-4);
      const double fine = remainder(1e-4);
      contract::require(coarse > 1e-13 && fine <= 0.26 * coarse + 1e-13,
                        description);
    };
    verify_reduced_taylor(hminus1_reduced,
                          hminus1_control,
                          hminus1_evaluation,
                          hminus1_direction,
                          "H-1 search-direction Taylor remainder is not quadratic");
    verify_reduced_taylor(l2_reduced,
                          l2_control,
                          l2_evaluation,
                          l2_direction,
                          "L2 comparison-direction Taylor remainder is not quadratic");

    const auto &manifest = hminus1_compilation.problem->manifest();
    const auto &l2_manifest = l2_compilation.problem->manifest();
    contract::require(
      manifest.resolved_decision.h1_target_data_membership_selection.has_value() &&
        l2_manifest.resolved_decision.h1_target_data_membership_selection.has_value() &&
        manifest.resolved_decision.h1_target_data_membership_selection->data_id ==
          l2_manifest.resolved_decision.h1_target_data_membership_selection->data_id &&
        manifest.resolved_decision.h1_target_data_membership_selection
            ->fixed_boundary_region_id == "dirichlet_boundary" &&
        std::any_of(manifest.compatibility.declared_assumptions.begin(),
                    manifest.compatibility.declared_assumptions.end(),
                    [](const std::string &assumption) {
                      return assumption.find(
                               "h1_target_data_membership: status=user_assumed") ==
                             0;
                    }) &&
      hminus1_compilation.problem->metric().id() == "hminus1_continuous" &&
        l2_compilation.problem->metric().id() == "l2_continuous" &&
        manifest.resolved_decision.metric_record.operator_description.find("M_h K_h^{-1} M_h") !=
          std::string::npos &&
        manifest.resolved_decision.metric_record.operator_id ==
          "mass_laplacian_inverse_mass" &&
        manifest.resolved_decision.metric_record.inverse_operator_id ==
          "mass_inverse_laplacian_mass_inverse" &&
        manifest.resolved_decision.metric_record.boundary_region_id == "dirichlet_boundary" &&
        manifest.resolved_decision.metric_record.laplacian_solve_policy_id ==
          "control_metric_solve.laplacian_inverse" &&
        manifest.resolved_decision.metric_record.mass_solve_policy_id ==
          "control_metric_solve.mass_inverse" &&
        manifest.compatibility.metric_solve_policy.find("hminus1_continuous") !=
          std::string::npos &&
        manifest.compatibility.control_space.find("independent homogeneous-Dirichlet") !=
          std::string::npos,
      "H-1 compilation manifest omitted its operator, solve, or control-space policy");
  }

  template <int dim>
  void
  run_coefficient_identification_contract_test()
  {
    dealii::Triangulation<dim> triangulation;
    dealii::GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(1);

    const dealii::Functions::ConstantFunction<dim> forcing(1.0);
    const dealii::Functions::ConstantFunction<dim> desired_state(0.0);
    const auto specification =
      semantic::v1::make_coefficient_identification_problem();
    compiler::v1::DealiiDiscretisationPolicy policy;
    policy.state_degree = 1;
    const compiler::v1::DealiiCompiler compiler;
    contract::require(
      compiler.validate(specification, policy).valid(),
      "coefficient-identification v1 graph did not validate for deal.II");

    auto continuous_parameter = specification;
    continuous_parameter.spaces.at(2).topology =
      semantic::v1::SpaceTopology::h1;
    const auto continuous_parameter_report =
      compiler.validate(continuous_parameter, policy);
    test_support::require_exact_diagnostic(
      continuous_parameter_report,
      semantic::v1::DiagnosticCategory::lowerability,
      "diffusion_parameter",
      "cellwise_parameter_space",
      "coefficient-identification compiler did not require cellwise parameters");

    auto missing_positive_box = specification;
    missing_positive_box.constraints.clear();
    missing_positive_box.formulation.constraint_id.clear();
    const auto missing_positive_box_report =
      compiler.validate(missing_positive_box, policy);
    test_support::require_exact_diagnostic(
      missing_positive_box_report,
      semantic::v1::DiagnosticCategory::lowerability,
      "diffusion_parameter",
      "positive_parameter_constraint",
      "coefficient-identification compiler did not require the positive parameter box");

    const compiler::v1::DealiiDataBindings<dim> bindings{
      forcing,
      desired_state,
      std::nullopt,
      0.5,
      0.2,
      test_binding_provenance("coefficient_identification")};
    const compiler::v1::CellwiseBoxDataBindings bounds{
      compiler::v1::CellwiseBoundValue{0.2},
      compiler::v1::CellwiseBoundValue{2.0}};
    const auto compilation = compiler.compile(specification,
                                              triangulation,
                                              bindings,
                                              policy,
                                              bounds);
    contract::require(compilation.succeeded(),
                      "coefficient-identification v1 compilation failed");

    const compiler::v1::CellwiseBoxDataBindings nonpositive_bounds{
      compiler::v1::CellwiseBoundValue{0.0},
      compiler::v1::CellwiseBoundValue{2.0}};
    const auto rejected_nonpositive_bounds = compiler.compile(specification,
                                                              triangulation,
                                                              bindings,
                                                              policy,
                                                              nonpositive_bounds);
    contract::require(
      !rejected_nonpositive_bounds.succeeded(),
      "coefficient-identification compiler accepted a nonpositive lower bound");
    test_support::require_exact_diagnostic(
      rejected_nonpositive_bounds.diagnostics,
      semantic::v1::DiagnosticCategory::lowerability,
      "parameter_box",
      "positive_parameter_lower_bound",
      "coefficient-identification compiler did not identify a nonpositive lower bound");

    const auto &model = compilation.problem->executable_model();
    const auto reduced = compilation.problem->make_reduced_dto();
    dealii::Vector<double> nonpositive_lower(
      model.variable_layout()->dimension(1));
    dealii::Vector<double> positive_upper(
      model.variable_layout()->dimension(1));
    for (dealii::types::global_dof_index index = 0;
         index < nonpositive_lower.size();
         ++index)
      {
        nonpositive_lower[index] = index == 0 ? 0.0 : 0.2;
        positive_upper[index] = 2.0;
      }
    const compiler::v1::CellwiseBoxDataBindings nonpositive_vector_bounds{
      compiler::v1::CellwiseBoundValue{std::move(nonpositive_lower)},
      compiler::v1::CellwiseBoundValue{std::move(positive_upper)}};
    const auto rejected_nonpositive_vector_bounds =
      compiler.compile(specification,
                       triangulation,
                       bindings,
                       policy,
                       nonpositive_vector_bounds);
    contract::require(
      !rejected_nonpositive_vector_bounds.succeeded(),
      "coefficient-identification compiler accepted a nonpositive vector lower bound");
    test_support::require_exact_diagnostic(
      rejected_nonpositive_vector_bounds.diagnostics,
      semantic::v1::DiagnosticCategory::lowerability,
      "parameter_box",
      "positive_parameter_lower_bound",
      "coefficient-identification compiler did not identify a nonpositive vector lower bound");
    dealii::Vector<double> parameter_values(
      model.variable_layout()->dimension(1));
    for (dealii::types::global_dof_index index = 0;
         index < parameter_values.size();
         ++index)
      parameter_values[index] = index % 2 == 0 ? 0.6 : 1.1;
    const Primal parameter(
      model.variable_layout()->single_block(1, "parameter"),
      {parameter_values});
    const auto evaluation = reduced.evaluate(parameter);
    require_close(model.residual(evaluation.full_point).block(0).l2_norm(),
                  0.0,
                  1e-11,
                  "coefficient-identification state residual");

    dealii::Vector<double> changed_parameter_values = parameter_values;
    for (dealii::types::global_dof_index index = 0;
         index < changed_parameter_values.size();
         ++index)
      changed_parameter_values[index] = 1.5;
    const Primal changed_parameter(parameter.layout(),
                                   {std::move(changed_parameter_values)});
    const auto changed_evaluation = reduced.evaluate(changed_parameter);
    require_close(model.residual(changed_evaluation.full_point).block(0).l2_norm(),
                  0.0,
                  1e-11,
                  "coefficient-identification reassembled state residual");
    dealii::Vector<double> changed_state_difference = changed_evaluation.state.block(0);
    changed_state_difference.add(-1.0, evaluation.state.block(0));
    contract::require(changed_state_difference.l2_norm() > 1e-4,
                      "coefficient-identification state matrix was not reassembled");

    dealii::Vector<double> state_tangent(model.variable_layout()->dimension(0));
    dealii::Vector<double> parameter_tangent(model.variable_layout()->dimension(1));
    for (dealii::types::global_dof_index index = 0;
         index < state_tangent.size();
         ++index)
      state_tangent[index] = -0.01 * static_cast<double>(index + 1);
    for (dealii::types::global_dof_index index = 0;
         index < parameter_tangent.size();
         ++index)
      parameter_tangent[index] = index % 2 == 0 ? 0.03 : -0.02;
    const Primal tangent(model.variable_layout(),
                         {std::move(state_tangent), std::move(parameter_tangent)});
    dealii::Vector<double> test_seed_values(model.test_layout()->dimension(0));
    for (dealii::types::global_dof_index index = 0;
         index < test_seed_values.size();
         ++index)
      test_seed_values[index] = 0.04 * static_cast<double>(index + 1);
    const Primal test_seed(model.test_layout(), {std::move(test_seed_values)});
    const Covector jvp = model.residual_jvp(evaluation.full_point, tangent);
    const Covector vjp = model.residual_vjp(evaluation.full_point, test_seed);
    require_close(contract::pair(jvp, test_seed),
                  contract::pair(vjp, tangent),
                  1e-11,
                  "coefficient-identification residual JVP/VJP pairing");

    constexpr double derivative_step = 1e-7;
    const Covector residual_at_step =
      model.residual(shifted(evaluation.full_point, tangent, derivative_step));
    const Covector residual_at_point = model.residual(evaluation.full_point);
    for (std::size_t block = 0; block < residual_at_step.n_blocks(); ++block)
      {
        dealii::Vector<double> finite_difference = residual_at_step.block(block);
        finite_difference.add(-1.0, residual_at_point.block(block));
        finite_difference *= 1.0 / derivative_step;
        finite_difference.add(-1.0, jvp.block(block));
        require_close(finite_difference.l2_norm(),
                      0.0,
                      1e-7,
                      "coefficient-identification residual finite-difference JVP");
      }

    const double objective_difference =
      model.objective(shifted(evaluation.full_point, tangent, derivative_step)) -
      evaluation.objective_value;
    require_close(
      objective_difference / derivative_step,
      contract::pair(model.objective_derivative(evaluation.full_point), tangent),
      1e-7,
      "coefficient-identification objective directional derivative");

    dealii::Vector<double> direction_values(parameter.layout()->dimension(0));
    for (dealii::types::global_dof_index index = 0;
         index < direction_values.size();
         ++index)
      direction_values[index] = index % 2 == 0 ? 0.02 : -0.03;
    const Primal direction(parameter.layout(), {std::move(direction_values)});
    const double directional_derivative =
      contract::pair(evaluation.reduced_derivative, direction);
    const auto remainder = [&](const double step) {
      return std::abs(reduced.evaluate(shifted(parameter, direction, step))
                        .objective_value -
                      evaluation.objective_value - step * directional_derivative);
    };
    const double coarse_remainder = remainder(1e-3);
    const double fine_remainder = remainder(5e-4);
    contract::require(coarse_remainder > 1e-12 &&
                        fine_remainder <= 0.26 * coarse_remainder + 1e-13,
                      "coefficient-identification reduced Taylor remainder is not quadratic");

    const auto &metric = compilation.problem->metric();
    const Primal metric_direction =
      metric.inverse_apply(evaluation.reduced_derivative);
    require_covector_close(metric.apply(metric_direction),
                           evaluation.reduced_derivative,
                           1e-10,
                           "coefficient-identification metric inverse/apply");
    const auto *constraint = compilation.problem->constraint();
    contract::require(constraint != nullptr && constraint->is_feasible(parameter),
                      "coefficient-identification compilation omitted the parameter box");
    contract::require(metric.id() == "l2_cellwise_parameter",
                      "coefficient-identification compilation selected the wrong metric");
    contract::require(
        compilation.problem->manifest().compatibility.state_adjoint_solve_policy.find(
        "reassembled") != std::string::npos,
      "coefficient-identification manifest omitted state-matrix reassembly");
    require_constraint_realisation(
      compilation.problem->manifest(),
      "FE_DGQ(0) coefficientwise l2_cellwise_parameter clipping",
      "coefficient-identification");
  }

