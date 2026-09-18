#pragma once

// These scenario implementations are deliberately included into one heavy
// deal.II/compiler translation unit to improve source navigation without
// multiplying compilation cost.

  template <int dim>
  void
  run_neumann_boundary_contract_test()
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
    const dealii::Functions::ConstantFunction<dim> desired_state(0.2);
    const auto specification =
      semantic::v1::make_neumann_boundary_control_problem(true);
    compiler::v1::DealiiDiscretisationPolicy policy;
    policy.state_degree = 1;
    const compiler::v1::DealiiCompiler compiler;
    contract::require(compiler.validate(specification, policy).valid(),
                      "Neumann boundary-control v1 graph did not validate for deal.II");
    auto unselected_observation_policy = policy;
    unselected_observation_policy.volume_observation =
      compiler::v1::VolumeObservationDiscretisationPolicy{
        3,
        compiler::v1::VolumeObservationTargetRealisation::
          analytic_quadrature};
    test_support::require_exact_diagnostic(
      compiler.validate(specification, unselected_observation_policy),
      semantic::v1::DiagnosticCategory::lowerability,
      specification.id,
      "unselected_volume_observation_policy",
      "A boundary-observation target accepted a volume-observation policy");
    const compiler::v1::DealiiDataBindings<dim> bindings{
      forcing,
      desired_state,
      1.0,
      0.5,
      0.1,
      test_binding_provenance("neumann_boundary")};
    const compiler::v1::FacewiseBoxDataBindings facewise_bounds{
      compiler::v1::FacewiseBoundValue{-0.25},
      compiler::v1::FacewiseBoundValue{0.4}};
    const auto compilation = compiler.compile(specification,
                                              triangulation,
                                              bindings,
                                              policy,
                                              std::nullopt,
                                              facewise_bounds);
    contract::require(compilation.succeeded(),
                      "Neumann boundary-control v1 compilation failed");

    const auto &model = compilation.problem->executable_model();
    const auto *neumann_model =
      dynamic_cast<const compiler::v1::detail::NeumannBoundaryControlModel<dim> *>(
        &model);
    contract::require(neumann_model != nullptr,
                      "Neumann boundary-control did not produce its Neumann target");
    const dealii::QGauss<dim - 1> boundary_quadrature(policy.state_degree + 2);
    std::vector<dealii::Point<dim>> expected_control_coordinates;
    std::vector<double>             expected_control_measures;
    for (auto cell = triangulation.begin_active();
         cell != triangulation.end();
         ++cell)
      for (unsigned int face = 0;
           face < dealii::GeometryInfo<dim>::faces_per_cell;
           ++face)
        if (cell->face(face)->at_boundary() &&
            cell->face(face)->boundary_id() == 1)
          {
            expected_control_coordinates.push_back(cell->face(face)->center());
            expected_control_measures.push_back(cell->face(face)->measure());
          }
    contract::require(
      neumann_model->physical_control_dimension() ==
          expected_control_coordinates.size() &&
        neumann_model->independent_control_dimension() ==
          expected_control_coordinates.size() &&
        neumann_model->control_coordinates().size() ==
          expected_control_coordinates.size(),
      "Facewise Neumann realization changed the control dimensions");
    for (std::size_t index = 0;
         index < expected_control_coordinates.size();
         ++index)
      require_close(
        neumann_model->control_coordinates()[index].distance(
          expected_control_coordinates[index]),
        0.0,
        1e-15,
        "Facewise Neumann realization changed control-face ordering");

    const auto &metric = compilation.problem->metric();
    dealii::Vector<double> unit_control_values(
      expected_control_coordinates.size());
    unit_control_values = 1.0;
    const Primal unit_control(
      model.variable_layout()->single_block(1, "control"),
      {std::move(unit_control_values)});
    const Covector measured_control = metric.apply(unit_control);
    for (std::size_t index = 0; index < expected_control_measures.size(); ++index)
      require_close(measured_control.block(0)[index],
                    expected_control_measures[index],
                    1e-14,
                    "Facewise Neumann realization changed its mass matrix");

    std::size_t expected_boundary_trace_samples = 0;
    for (auto cell = triangulation.begin_active();
         cell != triangulation.end();
         ++cell)
      for (unsigned int face = 0;
           face < dealii::GeometryInfo<dim>::faces_per_cell;
           ++face)
        if (cell->face(face)->at_boundary() &&
            cell->face(face)->boundary_id() == 2)
          expected_boundary_trace_samples += boundary_quadrature.size();
    const auto boundary_trace_values =
      neumann_model->boundary_trace_values(
        Primal(model.variable_layout(),
               {dealii::Vector<double>(model.variable_layout()->dimension(0)),
                dealii::Vector<double>(model.variable_layout()->dimension(1))}));
    contract::require(
      boundary_trace_values.size() == expected_boundary_trace_samples,
      "Neumann boundary trace dimension did not count ordered face samples");
    const auto reduced = compilation.problem->make_reduced_dto();
    dealii::Vector<double> control_values(model.variable_layout()->dimension(1));
    for (dealii::types::global_dof_index index = 0;
         index < control_values.size();
         ++index)
      control_values[index] = 0.05 * static_cast<double>(index + 1);
    const Primal control(model.variable_layout()->single_block(1, "control"),
                         {std::move(control_values)});
    const auto evaluation = reduced.evaluate(control);
    require_close(model.residual(evaluation.full_point).block(0).l2_norm(),
                  0.0,
                  1e-11,
                  "Neumann boundary-control state residual");
    double expected_control_objective = 0.0;
    for (std::size_t index = 0; index < expected_control_measures.size(); ++index)
      expected_control_objective +=
        0.5 * bindings.regularisation_weight * expected_control_measures[index] *
        control.block(0)[index] * control.block(0)[index];
    require_close(
      neumann_model->objective_components(evaluation.full_point)
        .control_regularisation,
      expected_control_objective,
      1e-14,
      "Facewise Neumann realization changed the control objective");
    const Covector objective_derivative =
      model.objective_derivative(evaluation.full_point);
    for (std::size_t index = 0; index < expected_control_measures.size(); ++index)
      require_close(
        objective_derivative.block(1)[index],
        bindings.regularisation_weight * expected_control_measures[index] *
          control.block(0)[index],
        1e-14,
        "Facewise Neumann realization changed the control derivative");

    dealii::Vector<double> state_tangent(
      model.variable_layout()->dimension(0));
    dealii::Vector<double> control_tangent(
      model.variable_layout()->dimension(1));
    for (dealii::types::global_dof_index index = 0;
         index < control_tangent.size();
         ++index)
      control_tangent[index] =
        (index % 2 == 0 ? 0.03 : -0.02) * static_cast<double>(index + 1);
    const Primal coupling_tangent(model.variable_layout(),
                                  {std::move(state_tangent),
                                   std::move(control_tangent)});
    dealii::Vector<double> seed_values(model.test_layout()->dimension(0));
    for (dealii::types::global_dof_index index = 0;
         index < seed_values.size();
         ++index)
      seed_values[index] = 0.04 * static_cast<double>(index + 1);
    const Primal test_seed(model.test_layout(), {std::move(seed_values)});
    const Covector coupling_jvp =
      model.residual_jvp(evaluation.full_point, coupling_tangent);
    const Covector coupling_vjp =
      model.residual_vjp(evaluation.full_point, test_seed);
    require_close(contract::pair(coupling_jvp, test_seed),
                  contract::pair(coupling_vjp, coupling_tangent),
                  1e-11,
                  "Neumann boundary-control residual JVP/VJP pairing");

    dealii::FE_Q<dim> independent_state_fe(policy.state_degree);
    dealii::DoFHandler<dim> independent_state_dof_handler(triangulation);
    independent_state_dof_handler.distribute_dofs(independent_state_fe);
    contract::require(
      independent_state_dof_handler.n_dofs() == coupling_jvp.block(0).size(),
      "Independent Neumann coupling assembly changed the state dimension");
    dealii::FEFaceValues<dim> independent_face_values(
      independent_state_fe,
      boundary_quadrature,
      dealii::update_values | dealii::update_JxW_values);
    std::vector<dealii::types::global_dof_index> independent_state_indices(
      independent_state_fe.dofs_per_cell);
    dealii::Vector<double> expected_coupling(coupling_jvp.block(0).size());
    std::size_t control_index = 0;
    for (auto cell = independent_state_dof_handler.begin_active();
         cell != independent_state_dof_handler.end();
         ++cell)
      {
        cell->get_dof_indices(independent_state_indices);
        for (unsigned int face = 0;
             face < dealii::GeometryInfo<dim>::faces_per_cell;
             ++face)
          if (cell->face(face)->at_boundary() &&
              cell->face(face)->boundary_id() == 1)
            {
              independent_face_values.reinit(cell, face);
              for (unsigned int q = 0; q < boundary_quadrature.size(); ++q)
                for (unsigned int i = 0;
                     i < independent_state_fe.dofs_per_cell;
                     ++i)
                  expected_coupling[independent_state_indices[i]] -=
                    independent_face_values.shape_value(i, q) *
                    independent_face_values.JxW(q) *
                    coupling_tangent.block(1)[control_index];
              ++control_index;
            }
      }
    contract::require(control_index == expected_control_coordinates.size(),
                      "Independent Neumann coupling changed face ordering");
    dealii::Vector<double> coupling_error = coupling_jvp.block(0);
    coupling_error.add(-1.0, expected_coupling);
    require_close(coupling_error.l2_norm(),
                  0.0,
                  1e-14,
                  "Facewise Neumann realization changed its coupling matrix");

    dealii::Vector<double> trace_state_tangent(
      model.variable_layout()->dimension(0));
    dealii::Vector<double> trace_control_tangent(
      model.variable_layout()->dimension(1));
    for (dealii::types::global_dof_index index = 0;
         index < trace_state_tangent.size();
         ++index)
      trace_state_tangent[index] = 0.01 * static_cast<double>(index + 1);
    const Primal trace_tangent(model.variable_layout(),
                               {std::move(trace_state_tangent),
                                std::move(trace_control_tangent)});
    constexpr double derivative_step = 1e-7;
    const double trace_objective_difference =
      model.objective(shifted(evaluation.full_point, trace_tangent, derivative_step)) -
      evaluation.objective_value;
    require_close(
      trace_objective_difference / derivative_step,
      contract::pair(model.objective_derivative(evaluation.full_point), trace_tangent),
      2e-7,
      "boundary trace observation objective derivative");

    dealii::Vector<double> control_direction_values(control.layout()->dimension(0));
    for (dealii::types::global_dof_index index = 0;
         index < control_direction_values.size();
         ++index)
      control_direction_values[index] =
        (index % 2 == 0 ? 0.04 : -0.03) * static_cast<double>(index + 1);
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
                      "Neumann boundary-control reduced Taylor remainder is not quadratic");

    const Primal metric_direction =
      reduced.gradient_direction(evaluation.reduced_derivative, metric);
    const Covector metric_covector = metric.apply(metric_direction);
    require_close(contract::pair(metric_covector, control_direction),
                  contract::pair(evaluation.reduced_derivative, control_direction),
                  1e-11,
                  "facewise L2 metric pairing");

    const auto *constraint = compilation.problem->constraint();
    contract::require(constraint != nullptr,
                      "Neumann boundary compiler did not produce the facewise box");
    dealii::Vector<double> infeasible_values(control.layout()->dimension(0));
    for (dealii::types::global_dof_index index = 0;
         index < infeasible_values.size();
         ++index)
      infeasible_values[index] = 1.0;
    const Primal infeasible(control.layout(), {std::move(infeasible_values)});
    const Primal projected = constraint->project_in(infeasible, metric);
    contract::require(constraint->is_feasible(projected),
                      "facewise boundary box did not project into its feasible set");
    for (dealii::types::global_dof_index index = 0;
         index < projected.block(0).size();
         ++index)
      require_close(projected.block(0)[index],
                    0.4,
                    1e-12,
                    "facewise boundary box upper clipping");

    const auto &manifest = compilation.problem->manifest();
    require_constraint_realisation(
      manifest,
      "facewise-constant coefficientwise l2_facewise clipping",
      "Neumann-boundary control");
    contract::require(
      manifest.compatibility.control_space.find("facewise-constant") != std::string::npos &&
        manifest.compatibility.observation_realisation.find("boundary trace") !=
          std::string::npos &&
        manifest.compatibility.data_rule.find("boundary face quadrature") != std::string::npos &&
        manifest.compatibility.constraint_realisation.find("l2_facewise") != std::string::npos &&
        std::any_of(
          manifest.resolved_decision.realized_maps.begin(),
          manifest.resolved_decision.realized_maps.end(),
          [expected_boundary_trace_samples](
            const compiler::v1::CompiledRealizedMapRecord &map) {
            return map.semantic_id == "state_boundary_trace" &&
                   map.realization_id == "ordered_boundary_face_quadrature_trace" &&
                   map.output_dimension == expected_boundary_trace_samples;
          }),
      "Neumann boundary compilation manifest is incomplete");
  }

  template <int dim>
  void
  run_neumann_convection_subdomain_contract_test()
  {
    dealii::Triangulation<dim> triangulation;
    dealii::GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(2);
    for (auto cell = triangulation.begin_active();
         cell != triangulation.end();
         ++cell)
      {
        cell->set_material_id(cell->center()[0] < 0.5 ? 1 : 0);
        for (unsigned int face = 0;
             face < dealii::GeometryInfo<dim>::faces_per_cell;
             ++face)
          if (cell->face(face)->at_boundary())
            cell->face(face)->set_boundary_id(
              cell->face(face)->center()[0] < 1e-12 ? 0 : 1);
      }

    dealii::Tensor<1, dim> transport_value;
    transport_value[0] = -0.2;
    const ConstantVectorCoefficient<dim> conservative_transport(transport_value);
    const dealii::Tensor<1, dim> zero_transport_value;
    const ConstantVectorCoefficient<dim> zero_transport(zero_transport_value);
    const dealii::Functions::ConstantFunction<dim> forcing(0.3);
    const dealii::Functions::ConstantFunction<dim> desired_state(-0.15);
    const auto specification = semantic::v1::
      make_neumann_convection_subdomain_tracking_problem(1);
    compiler::v1::DealiiDiscretisationPolicy policy;
    policy.state_degree = 1;
    const compiler::v1::DealiiCompiler compiler;
    contract::require(compiler.validate(specification, policy).valid(),
                      "C5.6 Neumann convection graph did not validate for deal.II");
    auto zero_observation_quadrature = policy;
    zero_observation_quadrature.volume_observation =
      compiler::v1::VolumeObservationDiscretisationPolicy{
        0,
        compiler::v1::VolumeObservationTargetRealisation::
          analytic_quadrature};
    test_support::require_exact_diagnostic(
      compiler.validate(specification, zero_observation_quadrature),
      semantic::v1::DiagnosticCategory::lowerability,
      specification.id,
      "volume_observation_quadrature_order",
      "C5.6 compiler accepted a zero observation quadrature order");
    auto unknown_target_realisation = policy;
    unknown_target_realisation.volume_observation =
      compiler::v1::VolumeObservationDiscretisationPolicy{
        3,
        static_cast<compiler::v1::VolumeObservationTargetRealisation>(99)};
    test_support::require_exact_diagnostic(
      compiler.validate(specification, unknown_target_realisation),
      semantic::v1::DiagnosticCategory::lowerability,
      specification.id,
      "volume_observation_target_realisation",
      "C5.6 compiler accepted an unknown observation target realization");

    const compiler::v1::DealiiDataBindings<dim> missing_transport{
      forcing, desired_state, 1.0, 0.0, 0.2,
      test_binding_provenance("neumann_convection_missing")};
    const auto missing = compiler.compile(specification,
                                          triangulation,
                                          missing_transport,
                                          policy);
    test_support::require_exact_diagnostic(
      missing.diagnostics,
      semantic::v1::DiagnosticCategory::lowerability,
      "conservative_transport",
      "conservative_transport_data_binding",
      "C5.6 compiler did not require conservative transport data");

    const compiler::v1::DealiiConservativeTransportBindingProvenance provenance{
      "test.neumann_convection.conservative_transport"};
    const compiler::v1::DealiiDataBindings<dim> bindings{
      forcing,
      desired_state,
      1.0,
      0.0,
      0.2,
      test_binding_provenance("neumann_convection"),
      std::nullopt,
      std::nullopt,
      std::nullopt,
      compiler::v1::DealiiConservativeTransportDataBindings<dim>{
        conservative_transport, provenance}};
    const auto compilation = compiler.compile(specification,
                                              triangulation,
                                              bindings,
                                              policy);
    contract::require(compilation.succeeded(),
                      "C5.6 Neumann convection compilation failed");

    auto natural_source_specification = specification;
    natural_source_specification.spaces.push_back(
      {"natural_boundary_source_space",
       "Natural boundary source space",
       "control_boundary",
       semantic::v1::SpaceTopology::l2,
       semantic::v1::SpaceRole::data});
    natural_source_specification.data.push_back(
      {"natural_boundary_source",
       "Natural boundary source",
       semantic::v1::DataKind::function,
       semantic::v1::DataRole::natural_boundary_source,
       "natural_boundary_source_space"});
    natural_source_specification.residual_terms.push_back(
      {"natural_boundary_source",
       "Immutable natural boundary source",
       semantic::v1::ResidualTermKind::natural_boundary_source,
       "state_equation",
       {},
       {"natural_boundary_source"},
       "control_boundary"});
    component_by_id(natural_source_specification.equations, "state_equation")
      .residual_term_ids.push_back("natural_boundary_source");
    const auto missing_natural_source = compiler.compile(
      natural_source_specification, triangulation, bindings, policy);
    test_support::require_exact_diagnostic(
      missing_natural_source.diagnostics,
      semantic::v1::DiagnosticCategory::lowerability,
      natural_source_specification.id,
      "natural_boundary_source_data_binding",
      "C5.6 compiler did not require the natural-boundary source binding");

    const dealii::Functions::ConstantFunction<dim> natural_source(0.25);
    const dealii::Functions::ConstantFunction<dim> zero_natural_source(0.0);
    const auto make_natural_source_bindings =
      [&bindings](const dealii::Function<dim> &source,
                  const std::string &           provenance) {
        return compiler::v1::DealiiDataBindings<dim>{
          bindings.forcing,
          bindings.desired_state,
          bindings.diffusion,
          bindings.reaction,
          bindings.regularisation_weight,
          bindings.provenance,
          bindings.fixed_dirichlet_data,
          bindings.general_scalar,
          bindings.weighted_trace,
          bindings.conservative_transport,
          compiler::v1::DealiiNaturalBoundarySourceBinding<dim>{
            source, provenance}};
      };
    const auto natural_source_bindings = make_natural_source_bindings(
      natural_source, "test.neumann_convection.natural_boundary_source");
    const auto zero_natural_source_bindings = make_natural_source_bindings(
      zero_natural_source, "test.neumann_convection.zero_natural_source");
    const auto natural_source_compilation = compiler.compile(
      natural_source_specification,
      triangulation,
      natural_source_bindings,
      policy);
    const auto zero_natural_source_compilation = compiler.compile(
      natural_source_specification,
      triangulation,
      zero_natural_source_bindings,
      policy);
    contract::require(
      natural_source_compilation.succeeded() &&
        zero_natural_source_compilation.succeeded(),
      "C5.6 natural-boundary source compilation failed");
    const auto &natural_source_model =
      natural_source_compilation.problem->executable_model();
    const auto &zero_natural_source_model =
      zero_natural_source_compilation.problem->executable_model();
    dealii::Vector<double> natural_state(
      natural_source_model.variable_layout()->dimension(0));
    dealii::Vector<double> zero_natural_state = natural_state;
    dealii::Vector<double> natural_control(
      natural_source_model.variable_layout()->dimension(1));
    dealii::Vector<double> zero_natural_control = natural_control;
    const Primal natural_point(natural_source_model.variable_layout(),
                               {std::move(natural_state),
                                std::move(natural_control)});
    const Primal zero_natural_point(
      zero_natural_source_model.variable_layout(),
      {std::move(zero_natural_state), std::move(zero_natural_control)});
    auto natural_residual = natural_source_model.residual(natural_point);
    const auto zero_natural_residual =
      zero_natural_source_model.residual(zero_natural_point);
    natural_residual.add_scaled_block(0,
                                      -1.0,
                                      zero_natural_residual.block(0));
    contract::require(natural_residual.block(0).l2_norm() > 1e-8,
                      "C5.6 natural-boundary source did not enter the state load");

    dealii::Vector<double> natural_state_tangent(
      natural_source_model.variable_layout()->dimension(0));
    dealii::Vector<double> natural_control_tangent(
      natural_source_model.variable_layout()->dimension(1));
    dealii::Vector<double> zero_natural_state_tangent = natural_state_tangent;
    dealii::Vector<double> zero_natural_control_tangent = natural_control_tangent;
    for (dealii::types::global_dof_index index = 0;
         index < natural_state_tangent.size();
         ++index)
      natural_state_tangent[index] =
        zero_natural_state_tangent[index] = 0.017 * (index + 1);
    for (dealii::types::global_dof_index index = 0;
         index < natural_control_tangent.size();
         ++index)
      natural_control_tangent[index] =
        zero_natural_control_tangent[index] = -0.013 * (index + 1);
    const Primal natural_tangent(natural_source_model.variable_layout(),
                                 {std::move(natural_state_tangent),
                                  std::move(natural_control_tangent)});
    const Primal zero_natural_tangent(
      zero_natural_source_model.variable_layout(),
      {std::move(zero_natural_state_tangent),
       std::move(zero_natural_control_tangent)});
    auto natural_jvp =
      natural_source_model.residual_jvp(natural_point, natural_tangent);
    const auto zero_natural_jvp = zero_natural_source_model.residual_jvp(
      zero_natural_point, zero_natural_tangent);
    natural_jvp.add_scaled_block(0, -1.0, zero_natural_jvp.block(0));
    require_close(natural_jvp.block(0).l2_norm(),
                  0.0,
                  1e-12,
                  "C5.6 natural-boundary source changed the residual JVP");

    dealii::Vector<double> natural_seed(
      natural_source_model.test_layout()->dimension(0));
    dealii::Vector<double> zero_natural_seed = natural_seed;
    for (dealii::types::global_dof_index index = 0;
         index < natural_seed.size();
         ++index)
      natural_seed[index] = zero_natural_seed[index] = 0.021 * (index + 1);
    const Primal natural_test_seed(natural_source_model.test_layout(),
                                   {std::move(natural_seed)});
    const Primal zero_natural_test_seed(
      zero_natural_source_model.test_layout(),
      {std::move(zero_natural_seed)});
    auto natural_vjp =
      natural_source_model.residual_vjp(natural_point, natural_test_seed);
    const auto zero_natural_vjp = zero_natural_source_model.residual_vjp(
      zero_natural_point, zero_natural_test_seed);
    for (std::size_t block = 0; block < natural_vjp.n_blocks(); ++block)
      {
        natural_vjp.add_scaled_block(block,
                                     -1.0,
                                     zero_natural_vjp.block(block));
        require_close(natural_vjp.block(block).l2_norm(),
                      0.0,
                      1e-12,
                      "C5.6 natural-boundary source changed the residual VJP");
      }
    const auto &natural_manifest =
      natural_source_compilation.problem->manifest();
    contract::require(
      std::any_of(natural_manifest.resolved_decision.bindings.begin(),
                  natural_manifest.resolved_decision.bindings.end(),
                  [](const compiler::v1::CompiledBindingRecord &record) {
                    return record.role ==
                             semantic::v1::DataRole::natural_boundary_source &&
                           record.evaluation_realisation ==
                             "boundary_face_quadrature" &&
                           record.provenance ==
                             "test.neumann_convection.natural_boundary_source";
                  }) &&
        std::find(natural_manifest.compatibility.lowering_handler_records.begin(),
                  natural_manifest.compatibility.lowering_handler_records.end(),
                  "natural_boundary_source <- dealii.neumann.residual.natural_boundary_source") !=
          natural_manifest.compatibility.lowering_handler_records.end() &&
        natural_manifest.compatibility.data_rule.find(
          "natural-boundary source Function at selected boundary face quadrature") !=
          std::string::npos &&
        std::any_of(natural_manifest.compatibility.declared_assumptions.begin(),
                    natural_manifest.compatibility.declared_assumptions.end(),
                    [](const std::string &assumption) {
                      return assumption.find(
                               "natural_boundary_source: immutable scalar Function") !=
                             std::string::npos &&
                             assumption.find("scale=one") != std::string::npos &&
                             assumption.find("state/control derivatives are zero") !=
                               std::string::npos;
                    }),
      "C5.6 manifest omitted natural-boundary source provenance");
    auto explicit_analytic_policy = policy;
    explicit_analytic_policy.volume_observation =
      compiler::v1::VolumeObservationDiscretisationPolicy{
        policy.state_degree + 2,
        compiler::v1::VolumeObservationTargetRealisation::
          analytic_quadrature};
    const auto explicit_analytic_compilation = compiler.compile(
      specification, triangulation, bindings, explicit_analytic_policy);
    auto interpolated_policy = explicit_analytic_policy;
    interpolated_policy.volume_observation->quadrature_order = 2;
    interpolated_policy.volume_observation->target_realisation =
      compiler::v1::VolumeObservationTargetRealisation::
        state_fe_interpolation;
    const auto interpolated_compilation = compiler.compile(
      specification, triangulation, bindings, interpolated_policy);
    contract::require(
      explicit_analytic_compilation.succeeded() &&
        interpolated_compilation.succeeded(),
      "C5.6 compiler rejected a selected observation realization");
    const auto &model = compilation.problem->executable_model();
    const auto reduced = compilation.problem->make_reduced_dto();

    const compiler::v1::DealiiDataBindings<dim> zero_transport_bindings{
      forcing,
      desired_state,
      1.0,
      0.0,
      0.2,
      test_binding_provenance("neumann_zero_convection"),
      std::nullopt,
      std::nullopt,
      std::nullopt,
      compiler::v1::DealiiConservativeTransportDataBindings<dim>{
        zero_transport, {"test.neumann_convection.zero_transport"}}};
    const auto zero_transport_compilation = compiler.compile(
      specification, triangulation, zero_transport_bindings, policy);
    contract::require(zero_transport_compilation.succeeded(),
                      "zero-transport C5.6 comparison compilation failed");

    dealii::FE_Q<dim> oracle_fe(policy.state_degree);
    dealii::DoFHandler<dim> oracle_dof_handler(triangulation);
    oracle_dof_handler.distribute_dofs(oracle_fe);
    const FirstCoordinateFunction<dim> first_coordinate;
    dealii::Vector<double> coordinate_values(oracle_dof_handler.n_dofs());
    dealii::VectorTools::interpolate(oracle_dof_handler,
                                     first_coordinate,
                                     coordinate_values);
    dealii::Vector<double> zero_oracle_control(
      model.variable_layout()->dimension(1));
    const Primal coordinate_point(model.variable_layout(),
                                  {coordinate_values,
                                   std::move(zero_oracle_control)});
    const Primal coordinate_seed(model.test_layout(), {coordinate_values});
    const double transport_value_oracle =
      contract::pair(model.residual(coordinate_point), coordinate_seed) -
      contract::pair(
        zero_transport_compilation.problem->executable_model().residual(
          coordinate_point),
        coordinate_seed);
    require_close(transport_value_oracle,
                  0.1,
                  2e-12,
                  "C5.6 conservative transport weak-form value");
    require_close(model.objective(coordinate_point),
                  217.0 / 4800.0,
                  2e-12,
                  "C5.6 material-subdomain observation value");
    require_close(
      explicit_analytic_compilation.problem->executable_model().objective(
        coordinate_point),
      model.objective(coordinate_point),
      1e-14,
      "C5.6 implicit observation-policy compatibility");
    require_close(
      interpolated_compilation.problem->executable_model().objective(
        coordinate_point),
      model.objective(coordinate_point),
      1e-14,
      "C5.6 exactly represented interpolated target");

    dealii::Vector<double> control_values(model.variable_layout()->dimension(1));
    for (dealii::types::global_dof_index index = 0;
         index < control_values.size(); ++index)
      control_values[index] =
        (index % 2 == 0 ? 0.03 : -0.02) * static_cast<double>(index + 1);
    const Primal control(model.variable_layout()->single_block(1, "control"),
                         {std::move(control_values)});
    const auto evaluation = reduced.evaluate(control);
    require_close(model.residual(evaluation.full_point).block(0).l2_norm(),
                  0.0,
                  2e-11,
                  "C5.6 Neumann convection state residual");
    contract::require(
      evaluation.state_solve.algorithm == "serial_sparse_direct_umfpack" &&
        evaluation.adjoint_solve.algorithm ==
          "serial_sparse_direct_umfpack_transpose",
      "C5.6 Neumann convection did not use the exact nonsymmetric transpose solve");

    dealii::Vector<double> state_tangent(model.variable_layout()->dimension(0));
    dealii::Vector<double> control_tangent(model.variable_layout()->dimension(1));
    dealii::Vector<double> seed_values(model.test_layout()->dimension(0));
    for (dealii::types::global_dof_index index = 0;
         index < state_tangent.size(); ++index)
      {
        state_tangent[index] = 0.025 * static_cast<double>(index + 1);
        seed_values[index] = -0.035 * static_cast<double>(index + 1);
      }
    for (dealii::types::global_dof_index index = 0;
         index < control_tangent.size(); ++index)
      control_tangent[index] =
        (index % 2 == 0 ? -0.04 : 0.03) * static_cast<double>(index + 1);
    const Primal tangent(model.variable_layout(),
                         {std::move(state_tangent), std::move(control_tangent)});
    const Primal test_seed(model.test_layout(), {std::move(seed_values)});
    require_close(contract::pair(model.residual_jvp(evaluation.full_point, tangent),
                                 test_seed),
                  contract::pair(model.residual_vjp(evaluation.full_point,
                                                    test_seed),
                                 tangent),
                  3e-11,
                  "C5.6 conservative Neumann residual JVP/VJP pairing");
    constexpr double derivative_step = 1e-7;
    Covector residual_difference =
      model.residual(shifted(evaluation.full_point, tangent, derivative_step));
    const Covector base_residual = model.residual(evaluation.full_point);
    for (std::size_t block = 0; block < residual_difference.n_blocks(); ++block)
      {
        residual_difference.add_scaled_block(block,
                                             -1.0,
                                             base_residual.block(block));
        residual_difference.scale_block(block, 1.0 / derivative_step);
        residual_difference.add_scaled_block(
          block,
          -1.0,
          model.residual_jvp(evaluation.full_point, tangent).block(block));
        require_close(residual_difference.block(block).l2_norm(),
                      0.0,
                      2e-7,
                      "C5.6 conservative Neumann residual JVP");
      }

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
                      "C5.6 Neumann convection reduced Taylor remainder is not quadratic");

    const auto &manifest = compilation.problem->manifest();
    contract::require(
      manifest.compatibility.compiler_id ==
          "nmopt.compiler.v1.dealii.neumann_convection_subdomain" &&
        manifest.compatibility.observation_realisation ==
          "material-id volume restriction: 1; target=analytic-quadrature; "
          "quadrature=QGauss(3)" &&
        manifest.compatibility.data_rule.find(
          "analytic desired-state Function evaluated at selected QGauss(3) "
          "volume-observation quadrature") !=
          std::string::npos &&
        std::find(
          manifest.compatibility.lowering_handler_records.begin(),
          manifest.compatibility.lowering_handler_records.end(),
          "state_observation <- dealii.volume_observation.analytic-quadrature") !=
          manifest.compatibility.lowering_handler_records.end() &&
        manifest.resolved_decision.boundary_realisation.has_value() &&
        manifest.resolved_decision.boundary_realisation->id == "neumann_convection_partition" &&
        manifest.resolved_decision.boundary_realisation->fixed_dirichlet_region_id ==
          "dirichlet_boundary" &&
        manifest.resolved_decision.boundary_realisation->neumann_region_ids ==
          std::vector<std::string>{"control_boundary"} &&
        manifest.resolved_decision.boundary_realisation->conormal_form ==
          semantic::v1::ConormalForm::diffusion_minus_transport &&
        manifest.resolved_decision.boundary_realisation->normal_orientation ==
          semantic::v1::NormalOrientation::outward &&
        manifest.resolved_decision.boundary_realisation->transport_boundary_form ==
          semantic::v1::TransportBoundaryForm::total_conormal &&
        manifest.resolved_decision.state_solve_record.algorithm ==
          compiler::v1::LinearSolveAlgorithm::serial_sparse_direct_umfpack &&
        std::any_of(manifest.resolved_decision.bindings.begin(), manifest.resolved_decision.bindings.end(),
                    [](const compiler::v1::CompiledBindingRecord &record) {
                      return record.semantic_id == "conservative_transport" &&
                             record.provenance ==
                               "test.neumann_convection.conservative_transport";
                    }),
      "C5.6 manifest omitted transport, subdomain, or solve provenance");
    const auto &interpolated_manifest =
      interpolated_compilation.problem->manifest();
    contract::require(
      interpolated_manifest.compatibility.observation_realisation ==
          "material-id volume restriction: 1; target=state-fe-interpolation; "
          "quadrature=QGauss(2)" &&
        interpolated_manifest.compatibility.data_rule.find(
          "desired-state Function interpolated into scalar FE_Q(1) and "
          "evaluated at selected QGauss(2) volume-observation quadrature") !=
          std::string::npos &&
        std::find(
          interpolated_manifest.compatibility.lowering_handler_records.begin(),
          interpolated_manifest.compatibility.lowering_handler_records.end(),
          "state_observation <- dealii.volume_observation.state-fe-interpolation") !=
          interpolated_manifest.compatibility.lowering_handler_records.end(),
      "C5.6 manifest omitted the interpolated observation policy");
  }

  template <int dim>
  void
  run_pure_neumann_contract_test()
  {
    dealii::Triangulation<dim> triangulation;
    dealii::GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(1);

    const dealii::Functions::ConstantFunction<dim> zero_forcing(0.0);
    const dealii::Functions::ConstantFunction<dim> incompatible_forcing(1.0);
    const dealii::Functions::ConstantFunction<dim> desired_state(1.0);
    const auto specification =
      semantic::v1::make_pure_neumann_boundary_control_problem();
    compiler::v1::DealiiDiscretisationPolicy policy;
    policy.state_degree = 1;
    const compiler::v1::DealiiCompiler compiler;
    contract::require(
      compiler.validate(specification, policy).valid(),
      "pure-Neumann mean-constraint v1 graph did not validate for deal.II");

    const compiler::v1::DealiiDataBindings<dim> nonzero_reaction{
      zero_forcing,
      desired_state,
      1.0,
      0.5,
      0.1,
      test_binding_provenance("pure_neumann_nonzero_reaction")};
    const auto rejected_reaction = compiler.compile(specification,
                                                    triangulation,
                                                    nonzero_reaction,
                                                    policy);
    contract::require(
      !rejected_reaction.succeeded(),
      "pure-Neumann compiler did not reject a nonzero reaction");
    test_support::require_exact_diagnostic(
      rejected_reaction.diagnostics,
      semantic::v1::DiagnosticCategory::lowerability,
      "state",
      "pure_neumann_zero_reaction",
      "pure-Neumann compiler did not identify the nonzero reaction");

    const compiler::v1::DealiiDataBindings<dim> incompatible_binding{
      incompatible_forcing,
      desired_state,
      1.0,
      0.0,
      0.1,
      test_binding_provenance("pure_neumann_incompatible")};
    const auto rejected_forcing = compiler.compile(specification,
                                                   triangulation,
                                                   incompatible_binding,
                                                   policy);
    contract::require(
      !rejected_forcing.succeeded(),
      "pure-Neumann compiler did not reject incompatible forcing");
    test_support::require_exact_diagnostic(
      rejected_forcing.diagnostics,
      semantic::v1::DiagnosticCategory::lowerability,
      "forcing",
      "pure_neumann_forcing_compatibility",
      "pure-Neumann compiler did not identify incompatible forcing");

    const compiler::v1::DealiiDataBindings<dim> compatible_binding{
      zero_forcing,
      desired_state,
      1.0,
      0.0,
      0.1,
      test_binding_provenance("pure_neumann")};
    const auto compilation = compiler.compile(specification,
                                              triangulation,
                                              compatible_binding,
                                              policy);
    contract::require(compilation.succeeded(),
                      "pure-Neumann v1 compilation failed");

    const auto &model = compilation.problem->executable_model();
    const auto *pure_model =
      dynamic_cast<const compiler::v1::detail::NeumannBoundaryControlModel<dim> *>(
        &model);
    contract::require(pure_model != nullptr && pure_model->uses_mean_zero_gauge(),
                      "pure-Neumann compiler did not select the mean-zero target");

    const auto reduced = compilation.problem->make_reduced_dto();
    dealii::Vector<double> zero_values(model.variable_layout()->dimension(1));
    const Primal zero_control(model.variable_layout()->single_block(1, "control"),
                              {std::move(zero_values)});
    const auto evaluation = reduced.evaluate(zero_control);
    const auto repeated_evaluation = reduced.evaluate(zero_control);
    require_primal_close(evaluation.state,
                         repeated_evaluation.state,
                         1e-13,
                         "pure-Neumann state must not depend on an implicit pin");
    require_close(pure_model->state_mean(evaluation.state),
                  0.0,
                  1e-12,
                  "pure-Neumann state mean");
    require_close(pure_model->state_mean(evaluation.adjoint),
                  0.0,
                  1e-12,
                  "pure-Neumann adjoint mean");
    require_close(model.residual(evaluation.full_point).block(0).l2_norm(),
                  0.0,
                  1e-12,
                  "pure-Neumann compatible state residual");

    dealii::Vector<double> incompatible_values(zero_control.block(0).size());
    incompatible_values[0] = 0.1;
    const Primal incompatible_control(zero_control.layout(),
                                      {std::move(incompatible_values)});
    test_support::require_contract_error(
      [&reduced, &incompatible_control]() {
        (void)reduced.evaluate(incompatible_control);
      },
      "Pure-Neumann state load violates the discrete constant-mode compatibility condition",
      "pure-Neumann incompatible boundary control");

    const auto &manifest = compilation.problem->manifest();
    require_constraint_realisation(manifest, "none", "pure-Neumann");
    contract::require(
      manifest.compatibility.nullspace_policy.find("mean-zero Lagrange multiplier") !=
        std::string::npos &&
        manifest.compatibility.state_adjoint_solve_policy.find("SparseDirectUMFPACK") !=
          std::string::npos &&
        manifest.compatibility.lifting_realisation.find("pure-Neumann") != std::string::npos,
      "pure-Neumann compilation manifest omitted the selected gauge");
  }

  template <int dim>
  void
  run_general_scalar_robin_contract_test()
  {
    dealii::Triangulation<dim> triangulation;
    dealii::GridGenerator::hyper_cube(triangulation, 0.0, 1.0, true);
    triangulation.refine_global(2);

    dealii::Tensor<2, dim> identity_tensor;
    dealii::Tensor<2, dim> anisotropic_tensor;
    for (unsigned int direction = 0; direction < dim; ++direction)
      {
        identity_tensor[direction][direction] = 1.0;
        anisotropic_tensor[direction][direction] =
          1.4 + 0.3 * static_cast<double>(direction);
      }
    dealii::Tensor<1, dim> zero_vector;
    dealii::Tensor<1, dim> conservative_vector;
    dealii::Tensor<1, dim> advective_vector;
    conservative_vector[0] = 0.45;
    conservative_vector[1] = -0.15;
    advective_vector[0] = -0.2;
    advective_vector[1] = 0.35;

    const ConstantTensorCoefficient<dim> identity_diffusion(identity_tensor);
    const ConstantTensorCoefficient<dim> anisotropic_diffusion(
      anisotropic_tensor);
    const ConstantVectorCoefficient<dim> zero_transport(zero_vector);
    const ConstantVectorCoefficient<dim> conservative_transport(
      conservative_vector);
    const ConstantVectorCoefficient<dim> advective_transport(advective_vector);
    const dealii::Functions::ConstantFunction<dim> zero(0.0);
    const dealii::Functions::ConstantFunction<dim> reaction(0.6);
    const dealii::Functions::ConstantFunction<dim> robin_coefficient(1.1);
    const dealii::Functions::ConstantFunction<dim> robin_source(0.7);
    const dealii::Functions::ConstantFunction<dim> forcing(0.3);
    const dealii::Functions::ConstantFunction<dim> desired_state(-0.2);

    const auto specification =
      semantic::v1::make_general_scalar_elliptic_robin_problem(
        {0, 2, 3}, {1});
    compiler::v1::DealiiDiscretisationPolicy policy;
    policy.state_degree = 1;
    const compiler::v1::DealiiCompiler compiler;
    contract::require(compiler.validate(specification, policy).valid(),
                      "P5.1 general scalar graph did not validate for deal.II");

    const compiler::v1::DealiiGeneralScalarBindingProvenance
      coefficient_provenance{"test.general.diffusion_tensor",
                             "test.general.conservative_transport",
                             "test.general.advective_transport",
                             "test.general.reaction",
                             "test.general.robin_coefficient",
                             "test.general.robin_source"};
    const auto compile_with = [&](
                                const dealii::TensorFunction<2, dim> &diffusion,
                                const dealii::TensorFunction<1, dim> &conservative,
                                const dealii::TensorFunction<1, dim> &advective,
                                const dealii::Function<dim> &reaction_data,
                                const dealii::Function<dim> &robin_coefficient_data,
                                const dealii::Function<dim> &robin_source_data) {
      const compiler::v1::DealiiDataBindings<dim> bindings{
        forcing,
        desired_state,
        std::nullopt,
        0.0,
        0.2,
        test_binding_provenance("general_scalar_robin"),
        std::nullopt,
        compiler::v1::DealiiGeneralScalarDataBindings<dim>{
          diffusion,
          conservative,
          advective,
          reaction_data,
          robin_coefficient_data,
          robin_source_data,
          coefficient_provenance}};
      return compiler.compile(specification,
                              triangulation,
                              bindings,
                              policy);
    };

    const auto base = compile_with(identity_diffusion,
                                   zero_transport,
                                   zero_transport,
                                   zero,
                                   zero,
                                   zero);
    const auto diffusion_only = compile_with(anisotropic_diffusion,
                                             zero_transport,
                                             zero_transport,
                                             zero,
                                             zero,
                                             zero);
    const auto conservative_only = compile_with(identity_diffusion,
                                                conservative_transport,
                                                zero_transport,
                                                zero,
                                                zero,
                                                zero);
    const auto advective_only = compile_with(identity_diffusion,
                                             zero_transport,
                                             advective_transport,
                                             zero,
                                             zero,
                                             zero);
    const auto reaction_only = compile_with(identity_diffusion,
                                            zero_transport,
                                            zero_transport,
                                            reaction,
                                            zero,
                                            zero);
    const auto robin_bilinear_only = compile_with(identity_diffusion,
                                                  zero_transport,
                                                  zero_transport,
                                                  zero,
                                                  robin_coefficient,
                                                  zero);
    const auto robin_source_only = compile_with(identity_diffusion,
                                                zero_transport,
                                                zero_transport,
                                                zero,
                                                zero,
                                                robin_source);
    const auto combined = compile_with(anisotropic_diffusion,
                                       conservative_transport,
                                       advective_transport,
                                       reaction,
                                       robin_coefficient,
                                       robin_source);
    contract::require(base.succeeded() && diffusion_only.succeeded() &&
                        conservative_only.succeeded() &&
                        advective_only.succeeded() && reaction_only.succeeded() &&
                        robin_bilinear_only.succeeded() &&
                        robin_source_only.succeeded() && combined.succeeded(),
                      "P5.1 coefficient recombination did not compile");
    require_compiled_hessian_evidence(*combined.problem,
                                      "nonsymmetric general scalar");

    const auto &base_model = base.problem->executable_model();
    dealii::Vector<double> state_values(
      base_model.variable_layout()->dimension(0));
    dealii::Vector<double> control_values(
      base_model.variable_layout()->dimension(1));
    for (dealii::types::global_dof_index index = 0;
         index < state_values.size();
         ++index)
      state_values[index] = 0.04 * static_cast<double>(index + 1);
    for (dealii::types::global_dof_index index = 0;
         index < control_values.size();
         ++index)
      control_values[index] = -0.03 * static_cast<double>(index + 1);
    const Primal point(base_model.variable_layout(),
                       {state_values, control_values});
    const Primal tangent(base_model.variable_layout(),
                         {std::move(state_values), std::move(control_values)});
    dealii::Vector<double> seed_values(base_model.test_layout()->dimension(0));
    for (dealii::types::global_dof_index index = 0;
         index < seed_values.size();
         ++index)
      seed_values[index] =
        (index % 2 == 0 ? 0.05 : -0.035) *
        static_cast<double>(index + 1);
    const Primal seed(base_model.test_layout(), {std::move(seed_values)});

    using Model = contract::ExecutableModelT<Backend>;
    const std::vector<std::pair<std::string, const Model *>>
      variable_term_changes{
        {"tensor diffusion", &diffusion_only.problem->executable_model()},
        {"conservative transport",
         &conservative_only.problem->executable_model()},
        {"advective transport", &advective_only.problem->executable_model()},
        {"reaction", &reaction_only.problem->executable_model()},
        {"Robin bilinear", &robin_bilinear_only.problem->executable_model()}};
    const Covector base_residual = base_model.residual(point);
    const Covector base_jvp = base_model.residual_jvp(point, tangent);
    const Covector base_vjp = base_model.residual_vjp(point, seed);
    for (const auto &[name, changed_model] : variable_term_changes)
      {
        const Covector changed_residual = changed_model->residual(point);
        const Covector changed_jvp = changed_model->residual_jvp(point, tangent);
        const Covector changed_vjp = changed_model->residual_vjp(point, seed);
        dealii::Vector<double> value_jvp_difference = changed_residual.block(0);
        value_jvp_difference.add(-1.0, base_residual.block(0));
        value_jvp_difference.add(-1.0, changed_jvp.block(0));
        value_jvp_difference.add(1.0, base_jvp.block(0));
        require_close(value_jvp_difference.l2_norm(),
                      0.0,
                      2e-11,
                      name + " value/JVP action");

        dealii::Vector<double> term_action = changed_jvp.block(0);
        term_action.add(-1.0, base_jvp.block(0));
        contract::require(term_action.l2_norm() > 1e-8,
                          name + " term action vanished in its focused test");
        require_close(contract::pair(changed_jvp, seed) -
                        contract::pair(base_jvp, seed),
                      contract::pair(changed_vjp, tangent) -
                        contract::pair(base_vjp, tangent),
                      2e-11,
                      name + " JVP/VJP pairing");
      }

    const auto &source_model = robin_source_only.problem->executable_model();
    dealii::Vector<double> robin_load_delta =
      source_model.residual(point).block(0);
    robin_load_delta.add(-1.0, base_residual.block(0));
    contract::require(robin_load_delta.l2_norm() > 1e-8,
                      "Robin boundary source did not change the residual value");
    dealii::Vector<double> robin_source_jvp =
      source_model.residual_jvp(point, tangent).block(0);
    robin_source_jvp.add(-1.0, base_jvp.block(0));
    require_close(robin_source_jvp.l2_norm(),
                  0.0,
                  1e-13,
                  "Robin source derivative independence");

    const auto &model = combined.problem->executable_model();
    const Covector combined_jvp = model.residual_jvp(point, tangent);
    const Covector combined_vjp = model.residual_vjp(point, seed);
    require_close(contract::pair(combined_jvp, seed),
                  contract::pair(combined_vjp, tangent),
                  3e-11,
                  "general scalar combined JVP/VJP pairing");

    dealii::Vector<double> state_only_values(
      model.variable_layout()->dimension(0));
    dealii::Vector<double> state_seed_values(model.test_layout()->dimension(0));
    for (dealii::types::global_dof_index index = 0;
         index < state_only_values.size();
         ++index)
      {
        const double value = 0.025 * static_cast<double>(index + 1);
        state_only_values[index] = value;
        state_seed_values[index] = value;
      }
    dealii::Vector<double> zero_control(
      model.variable_layout()->dimension(1));
    const Primal state_only_tangent(model.variable_layout(),
                                    {std::move(state_only_values),
                                     std::move(zero_control)});
    const Primal state_seed(model.test_layout(),
                            {std::move(state_seed_values)});
    const Covector forward_state =
      model.residual_jvp(point, state_only_tangent);
    const Covector transpose_state = model.residual_vjp(point, state_seed);
    dealii::Vector<double> nonsymmetric_difference = forward_state.block(0);
    nonsymmetric_difference.add(-1.0, transpose_state.block(0));
    contract::require(
      nonsymmetric_difference.l2_norm() > 1e-7,
      "P5.1 transport target did not expose a nonsymmetric transpose action");

    dealii::Vector<double> reduced_control_values(
      model.variable_layout()->dimension(1));
    for (dealii::types::global_dof_index index = 0;
         index < reduced_control_values.size();
         ++index)
      reduced_control_values[index] =
        0.02 * std::sin(static_cast<double>(index + 1));
    const Primal reduced_control(
      model.variable_layout()->single_block(1, "control"),
      {std::move(reduced_control_values)});
    const auto reduced = combined.problem->make_reduced_dto();
    const auto evaluation = reduced.evaluate(reduced_control);
    require_close(model.residual(evaluation.full_point).block(0).l2_norm(),
                  0.0,
                  2e-11,
                  "general scalar state solve residual");
    contract::require(
      evaluation.state_solve.algorithm == "serial_sparse_direct_umfpack" &&
        evaluation.adjoint_solve.algorithm ==
          "serial_sparse_direct_umfpack_transpose",
      "P5.1 target did not report its nonsymmetric state/adjoint solves");

    dealii::Vector<double> direction_values(
      reduced_control.layout()->dimension(0));
    for (dealii::types::global_dof_index index = 0;
         index < direction_values.size();
         ++index)
      direction_values[index] =
        (index % 2 == 0 ? 0.04 : -0.03) *
        static_cast<double>(index + 1);
    const Primal direction(reduced_control.layout(),
                           {std::move(direction_values)});
    const double reduced_directional_derivative =
      contract::pair(evaluation.reduced_derivative, direction);
    const auto remainder = [&](const double step) {
      return std::abs(
        reduced.evaluate(shifted(reduced_control, direction, step))
          .objective_value -
        evaluation.objective_value - step * reduced_directional_derivative);
    };
    const double coarse_remainder = remainder(1e-3);
    const double fine_remainder = remainder(5e-4);
    contract::require(coarse_remainder > 1e-12 &&
                        fine_remainder <= 0.26 * coarse_remainder + 1e-13,
                      "P5.1 reduced Taylor remainder is not quadratic");

    const auto &manifest = combined.problem->manifest();
    require_constraint_realisation(manifest, "none", "general scalar Robin");
    const auto &boundary_selection = manifest.resolved_decision.boundary_realisation;
    const auto has_binding = [&manifest](
                               const std::string &semantic_id,
                               const semantic::v1::DataRole role,
                               const semantic::v1::DataKind kind,
                               const std::string &space_id,
                               const std::string &region_id,
                               const std::string &evaluation) {
      return std::any_of(
        manifest.resolved_decision.bindings.begin(),
        manifest.resolved_decision.bindings.end(),
        [&](const compiler::v1::CompiledBindingRecord &binding) {
          return binding.semantic_id == semantic_id && binding.role == role &&
                 binding.kind == kind && binding.space_id == space_id &&
                 binding.region_id == region_id &&
                 binding.evaluation_realisation == evaluation &&
                 !binding.provenance.empty();
        });
    };
    contract::require(
        manifest.resolved_decision.state_solve_record.algorithm ==
          compiler::v1::LinearSolveAlgorithm::serial_sparse_direct_umfpack &&
        manifest.resolved_decision.adjoint_solve_record.algorithm ==
          compiler::v1::LinearSolveAlgorithm::serial_sparse_direct_umfpack &&
        manifest.compatibility.lowering_handler_records.size() == 13 &&
        boundary_selection.has_value() &&
        boundary_selection->id == "scalar_boundary_partition" &&
        boundary_selection->fixed_dirichlet_region_id == "dirichlet_boundary" &&
        boundary_selection->robin_region_id == "robin_boundary" &&
        boundary_selection->neumann_region_ids.empty() &&
        boundary_selection->transport_inflow_region_ids.empty() &&
        boundary_selection->transport_outflow_region_id == "robin_boundary" &&
        boundary_selection->conormal_form ==
          semantic::v1::ConormalForm::diffusion_minus_transport &&
        boundary_selection->normal_orientation ==
          semantic::v1::NormalOrientation::outward &&
        boundary_selection->trace_realisation ==
          semantic::v1::TraceEvaluationRealisation::fe_q_state_trace &&
        boundary_selection->face_quadrature_realisation ==
          semantic::v1::FaceQuadratureRealisation::qgauss_face &&
        manifest.compatibility.data_rule.find("Robin coefficient and source") !=
          std::string::npos &&
        has_binding("diffusion_tensor",
                    semantic::v1::DataRole::diffusion,
                    semantic::v1::DataKind::tensor_function,
                    "diffusion_data_space",
                    "domain",
                    "volume_quadrature") &&
        has_binding("conservative_transport",
                    semantic::v1::DataRole::conservative_transport,
                    semantic::v1::DataKind::vector_function,
                    "conservative_transport_data_space",
                    "domain",
                    "volume_quadrature") &&
        has_binding("advective_transport",
                    semantic::v1::DataRole::advective_transport,
                    semantic::v1::DataKind::vector_function,
                    "advective_transport_data_space",
                    "domain",
                    "volume_quadrature") &&
        has_binding("reaction",
                    semantic::v1::DataRole::reaction,
                    semantic::v1::DataKind::function,
                    "reaction_data_space",
                    "domain",
                    "volume_quadrature") &&
        has_binding("robin_coefficient",
                    semantic::v1::DataRole::robin_coefficient,
                    semantic::v1::DataKind::function,
                    "robin_coefficient_data_space",
                    "robin_boundary",
                    "boundary_face_quadrature") &&
        has_binding("robin_source",
                    semantic::v1::DataRole::robin_source,
                    semantic::v1::DataKind::function,
                    "robin_source_data_space",
                    "robin_boundary",
                    "boundary_face_quadrature") &&
        std::any_of(manifest.compatibility.declared_assumptions.begin(),
                    manifest.compatibility.declared_assumptions.end(),
                    [](const std::string &assumption) {
                      return assumption.find(
                               "general_scalar_robin: boundary selection scalar_boundary_partition") ==
                             0;
                    }) &&
        std::find(manifest.compatibility.lowering_handler_records.begin(),
                  manifest.compatibility.lowering_handler_records.end(),
                  "conservative_transport <- "
                  "dealii.scalar.residual.conservative_transport") !=
          manifest.compatibility.lowering_handler_records.end(),
      "P5.1 manifest omitted coefficient, handler, or solve provenance");

    const auto overlapping =
      semantic::v1::make_general_scalar_elliptic_robin_problem(
        {0, 1, 2, 3}, {1});
    test_support::require_exact_diagnostic(
      compiler.validate(overlapping, policy),
      semantic::v1::DiagnosticCategory::lowerability,
      "robin_boundary",
      "scalar_boundary_partition_overlap",
      "P5.1 compiler accepted overlapping Dirichlet and Robin regions");

    const auto incomplete_partition =
      semantic::v1::make_general_scalar_elliptic_robin_problem({0}, {1});
    const compiler::v1::DealiiDataBindings<dim> incomplete_bindings{
      forcing,
      desired_state,
      std::nullopt,
      0.0,
      0.2,
      test_binding_provenance("general_incomplete_partition"),
      std::nullopt,
      compiler::v1::DealiiGeneralScalarDataBindings<dim>{
        anisotropic_diffusion,
        conservative_transport,
        advective_transport,
        reaction,
        robin_coefficient,
        robin_source,
        coefficient_provenance}};
    const auto incomplete = compiler.compile(incomplete_partition,
                                             triangulation,
                                             incomplete_bindings,
                                             policy);
    test_support::require_exact_diagnostic(
      incomplete.diagnostics,
      semantic::v1::DiagnosticCategory::lowerability,
      "state",
      "complete_scalar_boundary_partition",
      "P5.1 compiler accepted an incomplete boundary partition");

    const dealii::Functions::ConstantFunction<dim> vector_reaction(0.6, 2);
    const compiler::v1::DealiiDataBindings<dim> wrong_shape_bindings{
      forcing,
      desired_state,
      std::nullopt,
      0.0,
      0.2,
      test_binding_provenance("general_wrong_shape"),
      std::nullopt,
      compiler::v1::DealiiGeneralScalarDataBindings<dim>{
        anisotropic_diffusion,
        conservative_transport,
        advective_transport,
        vector_reaction,
        robin_coefficient,
        robin_source,
        coefficient_provenance}};
    const auto wrong_shape = compiler.compile(specification,
                                              triangulation,
                                              wrong_shape_bindings,
                                              policy);
    test_support::require_exact_diagnostic(
      wrong_shape.diagnostics,
      semantic::v1::DiagnosticCategory::lowerability,
      "reaction",
      "scalar_coefficient_function_shape",
      "P5.1 compiler accepted a multi-component reaction Function");
  }

