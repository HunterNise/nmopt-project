#pragma once

// These scenario implementations are deliberately included into one heavy
// deal.II/compiler translation unit to improve source navigation without
// multiplying compilation cost.

  template <int dim>
  void
  run_compiler_diagnostics_contract_test()
  {
    const compiler::v1::DealiiCompiler compiler;
    const auto specification =
      semantic::v1::make_scalar_diffusion_reaction_problem(true);

    dealii::Triangulation<dim> borrowed_mesh;
    dealii::GridGenerator::hyper_cube(borrowed_mesh);
    borrowed_mesh.refine_global(2);
    const dealii::Functions::ConstantFunction<dim> forcing(1.0);
    const dealii::Functions::ConstantFunction<dim> desired_state(0.25);
    const compiler::v1::CellwiseBoxDataBindings valid_bounds{
      compiler::v1::CellwiseBoundValue{-1.0},
      compiler::v1::CellwiseBoundValue{1.0}};
    const compiler::v1::DealiiDataBindings<dim> valid_bindings{
      forcing,
      desired_state,
      1.0,
      0.5,
      0.1,
      test_binding_provenance("compiler_contract")};

    const dealii::Functions::ConstantFunction<dim> vector_forcing(1.0, 2);
    const compiler::v1::DealiiDataBindings<dim> vector_forcing_bindings{
      vector_forcing,
      desired_state,
      1.0,
      0.5,
      0.1,
      test_binding_provenance("vector_forcing")};
    const auto rejected_forcing_shape = compiler.compile(
      specification,
      borrowed_mesh,
      vector_forcing_bindings,
      {},
      valid_bounds);
    test_support::require_exact_diagnostic(
      rejected_forcing_shape.diagnostics,
      semantic::v1::DiagnosticCategory::lowerability,
      "forcing",
      "scalar_function_binding_shape",
      "compiler did not route forcing Function shape through the resolved binding request");

    const dealii::Functions::ConstantFunction<dim> vector_desired_state(0.25, 2);
    const compiler::v1::DealiiDataBindings<dim> vector_desired_state_bindings{
      forcing,
      vector_desired_state,
      1.0,
      0.5,
      0.1,
      test_binding_provenance("vector_desired_state")};
    const auto rejected_desired_state_shape = compiler.compile(
      specification,
      borrowed_mesh,
      vector_desired_state_bindings,
      {},
      valid_bounds);
    test_support::require_exact_diagnostic(
      rejected_desired_state_shape.diagnostics,
      semantic::v1::DiagnosticCategory::lowerability,
      "desired_state",
      "scalar_function_binding_shape",
      "compiler did not route desired-state Function shape through the resolved binding request");

    const compiler::v1::DealiiDataBindings<dim> invalid_diffusion{
      forcing,
      desired_state,
      -1.0,
      0.5,
      0.1,
      test_binding_provenance("invalid_diffusion")};
    const auto rejected_diffusion = compiler.compile(specification,
                                                      borrowed_mesh,
                                                      invalid_diffusion,
                                                      {},
                                                      valid_bounds);
    test_support::require_exact_diagnostic(
      rejected_diffusion.diagnostics,
      semantic::v1::DiagnosticCategory::lowerability,
      "diffusion",
      "positive_finite_diffusion_binding",
      "compiler did not diagnose a negative diffusion binding");

    const compiler::v1::DealiiDataBindings<dim> invalid_regularisation{
      forcing,
      desired_state,
      1.0,
      0.5,
      0.0,
      test_binding_provenance("invalid_regularisation")};
    const auto rejected_regularisation = compiler.compile(
      specification,
      borrowed_mesh,
      invalid_regularisation,
      {},
      valid_bounds);
    test_support::require_exact_diagnostic(
      rejected_regularisation.diagnostics,
      semantic::v1::DiagnosticCategory::lowerability,
      "regularisation",
      "positive_finite_regularisation_binding",
      "compiler did not diagnose a zero regularisation binding");

    const compiler::v1::DealiiDataBindings<dim> missing_provenance{
      forcing,
      desired_state,
      1.0,
      0.5,
      0.1,
      {"", "", ""}};
    const auto rejected_provenance = compiler.compile(specification,
                                                       borrowed_mesh,
                                                       missing_provenance,
                                                       {},
                                                       valid_bounds);
    test_support::require_exact_diagnostic(
      rejected_provenance.diagnostics,
      semantic::v1::DiagnosticCategory::lowerability,
      "forcing",
      "forcing_binding_provenance",
      "compiler did not require forcing binding provenance");
    test_support::require_exact_diagnostic(
      rejected_provenance.diagnostics,
      semantic::v1::DiagnosticCategory::lowerability,
      "desired_state",
      "desired_state_binding_provenance",
      "compiler did not require desired-state binding provenance");

    const compiler::v1::CellwiseBoxDataBindings reversed_bounds{
      compiler::v1::CellwiseBoundValue{1.0},
      compiler::v1::CellwiseBoundValue{-1.0}};
    const auto rejected_order = compiler.compile(specification,
                                                  borrowed_mesh,
                                                  valid_bindings,
                                                  {},
                                                  reversed_bounds);
    test_support::require_exact_diagnostic(
      rejected_order.diagnostics,
      semantic::v1::DiagnosticCategory::lowerability,
      "control_box",
      "ordered_bound_values",
      "compiler did not diagnose reversed scalar bounds");

    dealii::Vector<double> short_lower(1);
    dealii::Vector<double> short_upper(1);
    const compiler::v1::CellwiseBoxDataBindings short_bounds{
      compiler::v1::CellwiseBoundValue{std::move(short_lower)},
      compiler::v1::CellwiseBoundValue{std::move(short_upper)}};
    const auto rejected_layout = compiler.compile(specification,
                                                   borrowed_mesh,
                                                   valid_bindings,
                                                   {},
                                                   short_bounds);
    test_support::require_exact_diagnostic(
      rejected_layout.diagnostics,
      semantic::v1::DiagnosticCategory::lowerability,
      "control_box",
      "cellwise_bound_layout",
      "compiler did not diagnose a bound/layout mismatch");

    compiler::v1::DealiiDiscretisationPolicy invalid_policy;
    invalid_policy.state_solve.relative_tolerance = 0.0;
    const auto rejected_policy = compiler.compile(specification,
                                                   borrowed_mesh,
                                                   valid_bindings,
                                                   invalid_policy,
                                                   valid_bounds);
    test_support::require_exact_diagnostic(
      rejected_policy.diagnostics,
      semantic::v1::DiagnosticCategory::lowerability,
      "reduced_dto",
      "valid_state_solve_policy",
      "compiler did not diagnose an invalid state-solve policy");

    dealii::Triangulation<dim> empty_mesh;
    const auto rejected_empty_mesh = compiler.compile(specification,
                                                       empty_mesh,
                                                       valid_bindings,
                                                       {},
                                                       valid_bounds);
    test_support::require_exact_diagnostic(
      rejected_empty_mesh.diagnostics,
      semantic::v1::DiagnosticCategory::lowerability,
      specification.id,
      "nonempty_triangulation",
      "compiler did not diagnose an empty triangulation");

    auto missing_boundary_specification = specification;
    missing_boundary_specification.regions.at(1).boundary_ids = {99};
    const auto rejected_boundary = compiler.compile(
      missing_boundary_specification,
      borrowed_mesh,
      valid_bindings,
      {},
      valid_bounds);
    test_support::require_exact_diagnostic(
      rejected_boundary.diagnostics,
      semantic::v1::DiagnosticCategory::lowerability,
      "state",
      "fixed_dirichlet_boundary_presence",
      "compiler did not diagnose a boundary id absent from the mesh");
  }

