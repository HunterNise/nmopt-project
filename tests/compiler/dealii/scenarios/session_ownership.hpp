#pragma once

// These scenario implementations are deliberately included into one heavy
// deal.II/compiler translation unit to improve source navigation without
// multiplying compilation cost.

  template <int dim>
  void
  run_compiler_session_contract_test()
  {
    const compiler::v1::DealiiCompiler compiler;
    const auto specification =
      semantic::v1::make_scalar_diffusion_reaction_problem(true);
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
      test_binding_provenance("compiler_session")};

    const std::shared_ptr<compiler::v1::DealiiCompilationSession<dim>>
      null_session;
    const auto rejected_null_session = compiler.compile(
      specification, null_session, valid_bindings);
    test_support::require_exact_diagnostic(
      rejected_null_session.diagnostics,
      semantic::v1::DiagnosticCategory::lowerability,
      specification.id,
      "compilation_session_presence",
      "compiler threw instead of returning a null-session lowerability diagnostic");

    compiler::v1::DealiiDiscretisationPolicy solve_policy;
    solve_policy.state_solve = {317, 2e-11, 3e-14};
    solve_policy.adjoint_solve = {419, 4e-11, 5e-14};

    auto detached = [&]() {
      auto mesh = std::make_unique<dealii::Triangulation<dim>>();
      dealii::GridGenerator::hyper_cube(*mesh);
      mesh->refine_global(2);
      auto session =
        std::make_shared<compiler::v1::DealiiCompilationSession<dim>>(
          std::move(mesh), "test.unit_square.refine_2");
      const auto compilation = compiler.compile(specification,
                                                session,
                                                valid_bindings,
                                                solve_policy,
                                                valid_bounds);
      contract::require(compilation.succeeded(),
                        "owned-session compiler contract setup failed");
      const auto &model = compilation.problem->executable_model();
      dealii::Vector<double> values(model.variable_layout()->dimension(1));
      Primal control(model.variable_layout()->single_block(1, "control"),
                     {std::move(values)});
      struct DetachedService
      {
        contract::ReducedDTOT<Backend>          reduced;
        Primal                                  control;
        compiler::v1::CompilationManifest manifest;
      };
      return DetachedService{compilation.problem->make_reduced_dto(),
                             std::move(control),
                             compilation.problem->manifest()};
    }();

    const auto evaluation = detached.reduced.evaluate(detached.control);
    contract::require(
      evaluation.state_solve.converged() &&
        evaluation.adjoint_solve.converged() &&
        evaluation.state_solve.maximum_iterations == 317 &&
        evaluation.adjoint_solve.maximum_iterations == 419,
      "detached compiled service did not retain its solve policies and reports");
    contract::require(
      detached.manifest.schema_version == 4 &&
        detached.manifest.resolved_decision.mesh_record.dimension ==
          static_cast<unsigned int>(dim) &&
        detached.manifest.resolved_decision.mesh_record.active_cells == 16 &&
        detached.manifest.resolved_decision.mesh_record.provenance ==
          "test.unit_square.refine_2" &&
        detached.manifest.resolved_decision.mesh_record.lifetime ==
          compiler::v1::MeshLifetimePolicy::owned_session &&
        detached.manifest.resolved_decision.formulation_record.kind ==
          semantic::v1::FormulationKind::reduced_dto &&
        detached.manifest.resolved_decision.state_solve_record.maximum_iterations == 317 &&
        detached.manifest.resolved_decision.adjoint_solve_record.maximum_iterations == 419 &&
        detached.manifest.resolved_decision.constraint_record.realisation_id == "l2_cellwise" &&
        detached.manifest.resolved_decision.constraint_record.projection_metric_id ==
          detached.manifest.resolved_decision.metric_record.realisation_id &&
        !detached.manifest.resolved_decision.mesh_record.structural_identity.empty() &&
        detached.manifest.resolved_decision.semantic_problem_id ==
          detached.manifest.resolved_decision.semantic_problem_id &&
        !detached.manifest.resolved_decision.spaces.empty() &&
        !detached.manifest.resolved_decision.bindings.empty(),
      "structured compilation manifest omitted resolved session decisions");
  }

  template <int dim>
  void
  run_supplied_otd_owned_session_lifetime_contract_test()
  {
    const compiler::v1::DealiiCompiler compiler;
    const auto specification =
      semantic::v1::make_scalar_diffusion_reaction_supplied_otd_problem(false);
    const dealii::Functions::ConstantFunction<dim> forcing(1.0);
    const dealii::Functions::ConstantFunction<dim> desired_state(0.25);
    const compiler::v1::DealiiDataBindings<dim> data_bindings{
      forcing,
      desired_state,
      1.0,
      0.5,
      0.1,
      test_binding_provenance("supplied_otd_owned_session")};

    using KKTProduct =
      contract::EqualityConstrainedQuadraticKKTProductT<Backend>;
    struct DetachedKKT
    {
      std::shared_ptr<const KKTProduct> product;
      compiler::v1::CompilationManifest manifest;
    };

    auto detached = [&]() {
      auto mesh = std::make_unique<dealii::Triangulation<dim>>();
      dealii::GridGenerator::hyper_cube(*mesh);
      mesh->refine_global(1);
      auto session =
        std::make_shared<compiler::v1::DealiiCompilationSession<dim>>(
          std::move(mesh), "test.supplied_otd.owned_session");
      const auto compilation =
        compiler.compile(specification, session, data_bindings);
      contract::require(
        compilation.succeeded() && !compilation.problem &&
          compilation.supplied_otd_problem,
        "owned-session compiler did not produce the supplied OTD product");
      const auto &supplied_otd = *compilation.supplied_otd_problem;
      const auto supplied_initial =
        Primal::zeros(supplied_otd.system().variable_layout());
      const auto supplied_solve =
        supplied_otd.system().solve(supplied_initial);
      contract::require(supplied_solve.report.converged(),
                        "owned-session supplied OTD solve did not converge");
      return DetachedKKT{
        std::make_shared<const KKTProduct>(
          contract::make_canonical_supplied_otd_kkt_product(
            supplied_otd.system())),
        supplied_otd.manifest()};
    }();

    const auto &product = *detached.product;
    contract::require(
      detached.manifest.resolved_decision.mesh_record.lifetime ==
        compiler::v1::MeshLifetimePolicy::owned_session,
      "supplied OTD owned-session manifest lost its mesh lifetime policy");

    const KKTProduct::Point point{
      KKTProduct::Primal::zeros(product.layout().primal),
      KKTProduct::Primal::zeros(product.layout().multiplier)};
    const auto residual = product.residual(point);
    const auto transpose = product.apply_kkt_transpose(
      {KKTProduct::Primal::zeros(product.layout().stationarity),
       KKTProduct::Primal::zeros(product.layout().equality)});
    const auto adjoint = product.multiplier_to_adjoint(point.multiplier);
    contract::require(
      residual.stationarity.layout()->compatible_with(
        *product.layout().stationarity) &&
        residual.equality.layout()->compatible_with(*product.layout().equality) &&
        transpose.primal.layout()->compatible_with(*product.layout().primal) &&
        transpose.multiplier.layout()->compatible_with(
          *product.layout().multiplier) &&
        adjoint.layout()->compatible_with(*product.layout().adjoint),
      "detached supplied KKT bridge lost a typed action or conversion");

    contract::QuadraticKKTSolverPolicy solve_policy;
    solve_policy.maximum_iterations = 100;
    const auto kkt_solve =
      dealii_backend::solve_serial_quadratic_kkt(product, solve_policy);
    contract::require(kkt_solve.report.converged(),
                      "detached supplied KKT bridge solve did not converge");
  }

