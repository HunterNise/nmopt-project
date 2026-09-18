#pragma once

// These scenario implementations are deliberately included into one heavy
// deal.II/compiler translation unit to improve source navigation without
// multiplying compilation cost.

  template <int dim>
  void
  run_compiled_dto_kkt_contract_test()
  {
    dealii::Triangulation<dim> triangulation;
    dealii::GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(1);

    const auto specification =
      semantic::v1::make_fixed_dirichlet_scalar_diffusion_reaction_problem();
    const dealii::Functions::ConstantFunction<dim> forcing(1.0);
    const dealii::Functions::ConstantFunction<dim> desired_state(0.25);
    const dealii::Functions::ConstantFunction<dim> fixed_dirichlet_data(0.0);
    auto bindings = compiler::v1::DealiiDataBindings<dim>{
      forcing,
      desired_state,
      1.0,
      0.5,
      0.1,
      test_binding_provenance("compiled_dto_kkt", true)};
    bindings.fixed_dirichlet_data = std::cref(fixed_dirichlet_data);
    compiler::v1::DealiiDiscretisationPolicy policy;
    policy.state_degree = 1;
    const compiler::v1::DealiiCompiler compiler;

    const auto validation = compiler.validate(
      specification,
      policy,
      compiler::v1::CompilationProduct::quadratic_kkt);
    contract::require(validation.valid(),
                      "canonical DTO KKT target did not validate");

    const auto compilation = compiler.compile(
      specification,
      triangulation,
      bindings,
      policy,
      std::nullopt,
      std::nullopt,
      compiler::v1::CompilationProduct::quadratic_kkt);
    contract::require(
      compilation.succeeded() && compilation.kkt_problem &&
        !compilation.problem && !compilation.supplied_otd_problem,
      "compiler did not produce the distinct canonical DTO KKT product");

    const auto &kkt = *compilation.kkt_problem;
    const auto &product = kkt.product();
    using KKTProduct = contract::EqualityConstrainedQuadraticKKTProductT<Backend>;
    const auto &manifest = kkt.manifest();
    const auto &record = manifest.resolved_decision.kkt_record;
    contract::require(
      record.present && record.product_id == "compiled.scalar.dto.kkt" &&
        record.construction_realisation.find("canonical DTO") !=
          std::string::npos &&
        !record.primal_layout.empty() && !record.multiplier_layout.empty() &&
        !record.adjoint_layout.empty() && !record.stationarity_layout.empty() &&
        !record.equality_layout.empty() &&
        !record.primal_stationarity_pairing.empty() &&
        !record.multiplier_equality_pairing.empty() &&
        record.primal_stationarity_pairing_ids.size() == 2 &&
        record.multiplier_equality_pairing_ids.size() == 1 &&
        record.multiplier_conversion.find("negative") != std::string::npos &&
        record.rank_condition_declared && record.kernel_positivity_declared &&
        record.symmetry == "symmetric_indefinite" &&
        record.solver_policy.find("MINRES") != std::string::npos &&
        record.preconditioner == "identity baseline" &&
        record.d_transpose_consistency_declared &&
        record.kkt_transpose_consistency_declared &&
        !record.transpose_consistency_policy.empty() &&
        record.action_provenance.size() == 5 &&
        record.assembled_block_provenance.size() == 6 &&
        manifest.resolved_decision.kkt_record.product_id == record.product_id,
      "compiled DTO KKT manifest omitted its structured product boundary");

    const auto zero_primal = Primal::zeros(product.layout().primal);
    const auto zero_multiplier = Primal::zeros(product.layout().multiplier);
    const KKTProduct::Point point{zero_primal, zero_multiplier};
    const auto residual = product.residual(point);
    contract::require(
      residual.stationarity.layout()->compatible_with(
        *product.layout().stationarity) &&
        residual.equality.layout()->compatible_with(*product.layout().equality),
      "compiled DTO KKT residual did not retain its declared layouts");
    const auto transpose = product.apply_kkt_transpose(
      {Primal::zeros(product.layout().stationarity),
       Primal::zeros(product.layout().equality)});
    contract::require(
      transpose.primal.layout()->compatible_with(*product.layout().primal) &&
        transpose.multiplier.layout()->compatible_with(
          *product.layout().multiplier) &&
        product.supports_minres(),
      "compiled DTO KKT action or solver symmetry declaration is incomplete");

    const auto constrained_specification =
      semantic::v1::make_scalar_diffusion_reaction_problem();
    const auto rejected = compiler.compile(
      constrained_specification,
      triangulation,
      bindings,
      policy,
      std::nullopt,
      std::nullopt,
      compiler::v1::CompilationProduct::quadratic_kkt);
    test_support::require_exact_diagnostic(
      rejected.diagnostics,
      semantic::v1::DiagnosticCategory::formulation_capability,
      "reduced_dto",
      "compiled_quadratic_kkt",
      "compiler inferred a KKT product for a noncanonical DTO target");

    const auto supplied_specification =
      semantic::v1::make_scalar_diffusion_reaction_supplied_otd_problem();
    const auto supplied_bindings = compiler::v1::DealiiDataBindings<dim>{
      forcing,
      desired_state,
      1.0,
      0.5,
      0.1,
      test_binding_provenance("compiled_supplied_otd_kkt")};
    const auto supplied_validation = compiler.validate(
      supplied_specification,
      policy,
      compiler::v1::CompilationProduct::quadratic_kkt);
    contract::require(supplied_validation.valid(),
                      "canonical supplied-OTD KKT target did not validate");
    const auto supplied_compilation = compiler.compile(
      supplied_specification,
      triangulation,
      supplied_bindings,
      policy,
      std::nullopt,
      std::nullopt,
      compiler::v1::CompilationProduct::quadratic_kkt);
    contract::require(
      supplied_compilation.succeeded() && supplied_compilation.kkt_problem &&
        !supplied_compilation.problem &&
        !supplied_compilation.pdas_problem &&
        !supplied_compilation.supplied_otd_problem,
      "compiler did not produce the supplied-OTD KKT product through its KKT boundary");

    const auto &supplied_kkt = *supplied_compilation.kkt_problem;
    const auto &supplied_manifest = supplied_kkt.manifest();
    const auto &supplied_record = supplied_manifest.resolved_decision.kkt_record;
    contract::require(
      supplied_manifest.resolved_decision.formulation_record.provenance ==
          semantic::v1::FormulationProvenance::supplied_otd &&
        supplied_manifest.resolved_decision.supplied_otd_record.present &&
        supplied_manifest.resolved_decision.supplied_otd_record.declaration.has_value() &&
        supplied_record.present &&
        supplied_record.product_id == "compiled.scalar.supplied_otd.kkt" &&
        supplied_record.construction_realisation.find("supplied-OTD") !=
          std::string::npos &&
        supplied_manifest.resolved_decision.kkt_record.product_id ==
          supplied_record.product_id,
      "compiled supplied-OTD KKT manifest lost formulation provenance");
    const auto &supplied_product = supplied_kkt.product();
    const auto supplied_transpose = supplied_product.apply_kkt_transpose(
      {Primal::zeros(supplied_product.layout().stationarity),
       Primal::zeros(supplied_product.layout().equality)});
    contract::require(
      supplied_transpose.primal.layout()->compatible_with(
        *supplied_product.layout().primal) &&
        supplied_transpose.multiplier.layout()->compatible_with(
          *supplied_product.layout().multiplier) &&
        supplied_product.supports_minres(),
      "compiled supplied-OTD KKT product lost its typed executable boundary");

    struct DetachedSuppliedKKT
    {
      KKTProduct                              product;
      contract::QuadraticKKTSolverPolicy      solver_policy;
    };
    const auto detached_supplied = [&] {
      auto mesh = std::make_unique<dealii::Triangulation<dim>>();
      dealii::GridGenerator::hyper_cube(*mesh);
      mesh->refine_global(1);
      const auto session =
        std::make_shared<compiler::v1::DealiiCompilationSession<dim>>(
          std::move(mesh), "test.compiled_supplied_otd_kkt.owned_session");
      const dealii::Functions::ConstantFunction<dim> detached_forcing(1.0);
      const dealii::Functions::ConstantFunction<dim> detached_desired(0.25);
      const compiler::v1::DealiiDataBindings<dim> detached_bindings{
        detached_forcing,
        detached_desired,
        1.0,
        0.5,
        0.1,
        test_binding_provenance("compiled_supplied_otd_kkt_detached")};
      const auto detached_compilation = compiler.compile(
        semantic::v1::make_scalar_diffusion_reaction_supplied_otd_problem(),
        session,
        detached_bindings,
        policy,
        std::nullopt,
        std::nullopt,
        compiler::v1::CompilationProduct::quadratic_kkt);
      contract::require(detached_compilation.succeeded() &&
                          detached_compilation.kkt_problem,
                        "detached supplied-OTD KKT compilation failed");
      return DetachedSuppliedKKT{
        detached_compilation.kkt_problem->product(),
        policy.pdas_kkt_solver};
    }();
    const auto detached_result = dealii_backend::solve_serial_quadratic_kkt(
      detached_supplied.product, detached_supplied.solver_policy);
    contract::require(
      detached_result.report.converged() &&
        detached_result.report.linear_solve.converged(),
      "detached supplied-OTD KKT product lost its owned executable session");
  }

  template <int dim>
  void
  run_compiled_pdas_contract_test()
  {
    dealii::Triangulation<dim> triangulation;
    dealii::GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(1);

    const dealii::Functions::ConstantFunction<dim> forcing(1.0);
    const dealii::Functions::ConstantFunction<dim> desired_state(0.25);
    const compiler::v1::DealiiDataBindings<dim> data_bindings{
      forcing,
      desired_state,
      1.0,
      0.5,
      0.1,
      test_binding_provenance("compiled_pdas")};
    const compiler::v1::CellwiseBoxDataBindings bounds{
      compiler::v1::CellwiseBoundValue{-10.0},
      compiler::v1::CellwiseBoundValue{10.0}};
    compiler::v1::DealiiDiscretisationPolicy policy;
    policy.state_degree = 1;
    const compiler::v1::DealiiCompiler compiler;
    const auto specification =
      semantic::v1::make_scalar_diffusion_reaction_problem(true);

    const auto validation = compiler.validate(
      specification, policy, compiler::v1::CompilationProduct::pdas);
    contract::require(validation.valid(),
                      "canonical DTO PDAS target did not validate");
    const auto compilation = compiler.compile(
      specification,
      triangulation,
      data_bindings,
      policy,
      bounds,
      std::nullopt,
      compiler::v1::CompilationProduct::pdas);
    contract::require(
      compilation.succeeded() && compilation.pdas_problem &&
        !compilation.problem && !compilation.kkt_problem &&
        !compilation.supplied_otd_problem,
      "compiler did not produce the distinct canonical DTO PDAS product");

    const auto &pdas = *compilation.pdas_problem;
    const auto &manifest = pdas.manifest();
    const auto &record = manifest.resolved_decision.pdas_record;
    contract::require(
      record.present && record.product_id == "compiled.scalar.dto.pdas" &&
        record.bound_source.find("cellwise_box") != std::string::npos &&
        record.bound_realisation.find("FE_DGQ(0)") != std::string::npos &&
        record.control_block == 1 &&
        record.multiplier_representation.find("metric Riesz map") !=
          std::string::npos &&
        record.positive_diagonal_metric_declared &&
        record.classification_parameter ==
          policy.pdas.classification_parameter &&
        record.maximum_iterations == policy.pdas.maximum_iterations &&
        record.inner_kkt_solver == "serial MINRES" &&
        record.inner_kkt_maximum_iterations ==
          policy.pdas_kkt_solver.maximum_iterations &&
        record.active_set_rank_condition_declared &&
        record.active_set_kernel_positivity_declared &&
        manifest.resolved_decision.pdas_record.product_id ==
          record.product_id,
      "compiled DTO PDAS manifest omitted its structured product boundary");

    const auto &box_data = pdas.box_data();
    const auto *projection = dynamic_cast<
      const dealii_backend::CellwiseBoxConstraint *>(pdas.constraint());
    contract::require(
      projection != nullptr &&
        projection->supports_projection_in(pdas.metric()) &&
        projection->box_data_token() == box_data.token() &&
        pdas.complementarity().box_data_token() == box_data.token() &&
        record.box_data_token == box_data.token_id() &&
        record.bounds_digest == box_data.bounds_digest() &&
        record.control_ordering == box_data.layout_signature() &&
        record.data_provenance == box_data.data_provenance() &&
        manifest.resolved_decision.constraint_record.box_data_token == box_data.token_id() &&
        manifest.resolved_decision.constraint_record.bounds_digest == box_data.bounds_digest() &&
        manifest.resolved_decision.constraint_record.control_layout == box_data.layout_signature() &&
        manifest.resolved_decision.constraint_record.metric_identity == pdas.metric().id(),
      "compiled DTO projection and complementarity did not share box data");
    for (std::size_t index = 0;
         index < box_data.layout()->dimension(0);
         ++index)
      {
        require_close(
          box_data.lower().block(0)[index],
          pdas.complementarity().bounds().lower().block(0)[index],
          0.0,
          "compiled DTO shared lower box data");
        require_close(
          box_data.upper().block(0)[index],
          pdas.complementarity().bounds().upper().block(0)[index],
          0.0,
          "compiled DTO shared upper box data");
      }
    const compiler::v1::CellwiseBoxDataBindings changed_bounds{
      compiler::v1::CellwiseBoundValue{-9.0},
      compiler::v1::CellwiseBoundValue{10.0}};
    const auto changed_compilation = compiler.compile(
      specification,
      triangulation,
      data_bindings,
      policy,
      changed_bounds,
      std::nullopt,
      compiler::v1::CompilationProduct::pdas);
    contract::require(
      changed_compilation.succeeded() && changed_compilation.pdas_problem &&
        changed_compilation.pdas_problem->box_data().token_id() !=
          box_data.token_id(),
      "compiled DTO changed bounds did not receive a distinct box token");
    test_support::require_contract_error(
      [&box_data, &changed_compilation] {
        contract::require(
          box_data.compatible_with(changed_compilation.pdas_problem->box_data()),
          "compiled box data token mismatch");
      },
      "compiled box data token mismatch",
      "compiled DTO changed bounds were not rejected before comparison");

    using Product = contract::EqualityConstrainedQuadraticKKTProductT<Backend>;
    using Primal = Product::Primal;
    const auto initial = Product::Point{
      Primal::zeros(pdas.product().layout().primal),
      Primal::zeros(pdas.product().layout().multiplier)};
    const auto initial_box_multiplier =
      contract::CovectorBlockT<Backend>::zeros(pdas.complementarity().layout());
    const auto solver = pdas.make_solver(
      [&pdas](const Product &product) {
        return dealii_backend::solve_serial_quadratic_kkt(
          product, pdas.kkt_solver_policy());
      });
    const auto solve = solver.solve(initial,
                                    initial_box_multiplier,
                                    pdas.pdas_policy());
    contract::require(
      solve.converged() && !solve.iterations.empty() &&
        solve.iterations.back().kkt_solve.converged(),
      "compiled DTO PDAS product did not execute its owned KKT service");
    contract::require(
      solve.iterations.back().selection.active_size() == 0,
      "compiled DTO inactive PDAS product selected an active set");

    const auto direct_kkt = dealii_backend::solve_serial_quadratic_kkt(
      pdas.product(), pdas.kkt_solver_policy());
    contract::require(
      direct_kkt.report.converged(),
      "compiled DTO inactive reference KKT solve did not converge");
    require_primal_close(
      solve.solution.primal,
      direct_kkt.solution.primal,
      1e-9,
      "compiled DTO inactive PDAS primal disagrees with direct KKT");
    require_primal_close(
      solve.solution.multiplier,
      direct_kkt.solution.multiplier,
      1e-9,
      "compiled DTO inactive PDAS equality multiplier disagrees with direct KKT");

    const auto run_active_case = [&](const contract::BoxActivity expected_activity) {
      const std::size_t control_dimension =
        pdas.product().layout().primal->dimension(1);
      dealii::Vector<double> lower(control_dimension);
      dealii::Vector<double> upper(control_dimension);
      lower = -100.0;
      upper = 100.0;
      dealii::Vector<double> initial_control(control_dimension);
      initial_control = 0.0;
      const double reference_control = solve.solution.primal.block(1)[0];
      const double margin = std::max(0.25, 0.25 * std::abs(reference_control));
      double active_bound = 0.0;
      if (expected_activity == contract::BoxActivity::lower)
        {
          active_bound = reference_control + margin;
          lower[0] = active_bound;
          initial_control[0] = upper[0];
        }
      else
        {
          active_bound = reference_control - margin;
          upper[0] = active_bound;
          initial_control[0] = lower[0];
        }

      const compiler::v1::CellwiseBoxDataBindings active_bounds{
        compiler::v1::CellwiseBoundValue{std::move(lower)},
        compiler::v1::CellwiseBoundValue{std::move(upper)}};
      const auto active_compilation = compiler.compile(
        specification,
        triangulation,
        data_bindings,
        policy,
        active_bounds,
        std::nullopt,
        compiler::v1::CompilationProduct::pdas);
      contract::require(
        active_compilation.succeeded() && active_compilation.pdas_problem,
        "compiled DTO active PDAS case did not compile");
      const auto &active_pdas = *active_compilation.pdas_problem;
      const auto active_solver = active_pdas.make_solver(
        [&active_pdas](const Product &product) {
          return dealii_backend::solve_serial_quadratic_kkt(
            product, active_pdas.kkt_solver_policy());
        });
      const auto active_initial = Product::Point{
        Primal(active_pdas.product().layout().primal,
               {Backend::zeros(active_pdas.product().layout().primal->dimension(0)),
                std::move(initial_control)}),
        Primal::zeros(active_pdas.product().layout().multiplier)};
      const auto active_result = active_solver.solve(
        active_initial,
        Covector::zeros(active_pdas.complementarity().layout()),
        active_pdas.pdas_policy());

      contract::require(
        active_result.converged() && !active_result.iterations.empty(),
        "compiled DTO active PDAS case did not converge");
      bool saw_active_set_change = false;
      for (const auto &iteration : active_result.iterations)
        saw_active_set_change =
          saw_active_set_change || iteration.active_set_changes > 0;
      contract::require(
        saw_active_set_change,
        "compiled DTO active PDAS case did not report an active-set change");

      const auto &final_report = active_result.iterations.back();
      contract::require(
        final_report.selection.active_size() == 1 &&
          final_report.selection.activities().at(0) == expected_activity,
        "compiled DTO active PDAS case selected the wrong final active set");
      for (std::size_t index = 1; index < control_dimension; ++index)
        contract::require(
          final_report.selection.activities().at(index) ==
            contract::BoxActivity::inactive,
          "compiled DTO active PDAS case activated an unconstrained control");
      require_close(
        active_result.solution.primal.block(1)[0],
        active_bound,
        1e-8,
        "compiled DTO active PDAS case returned the wrong bound value");

      const auto represented_box_multiplier =
        active_pdas.complementarity().multiplier_to_primal(
          active_result.box_multiplier);
      const double active_multiplier = represented_box_multiplier.block(0)[0];
      if (expected_activity == contract::BoxActivity::lower)
        contract::require(
          active_multiplier < -1e-8,
          "compiled DTO lower-active PDAS multiplier has the wrong sign");
      else
        contract::require(
          active_multiplier > 1e-8,
          "compiled DTO upper-active PDAS multiplier has the wrong sign");
      for (std::size_t index = 1; index < control_dimension; ++index)
        require_close(
          represented_box_multiplier.block(0)[index],
          0.0,
          1e-10,
          "compiled DTO PDAS reported a nonzero inactive box multiplier");

      const auto residual = active_pdas.product().residual(
        active_result.solution);
      double stationarity_squared = 0.0;
      for (std::size_t block = 0;
           block < residual.stationarity.n_blocks();
           ++block)
        {
          dealii::Vector<double> constrained = residual.stationarity.block(block);
          if (block == 1)
            constrained.add(1.0, active_result.box_multiplier.block(0));
          stationarity_squared += constrained.l2_norm() * constrained.l2_norm();
        }
      double equality_squared = 0.0;
      for (std::size_t block = 0; block < residual.equality.n_blocks(); ++block)
        equality_squared += residual.equality.block(block).l2_norm() *
                            residual.equality.block(block).l2_norm();
      const double recovered_stationarity = std::sqrt(stationarity_squared);
      const double recovered_equality = std::sqrt(equality_squared);
      require_close(final_report.stationarity_residual,
                    recovered_stationarity,
                    1e-12,
                    "compiled DTO PDAS stationarity report was not recovered");
      require_close(final_report.equality_residual,
                    recovered_equality,
                    1e-12,
                    "compiled DTO PDAS equality report was not recovered");
      const auto &active_policy = active_pdas.pdas_policy();
      contract::require(
        std::isfinite(final_report.primal_violation) &&
          std::isfinite(final_report.dual_violation) &&
          std::isfinite(final_report.complementarity_residual) &&
          std::isfinite(final_report.stationarity_residual) &&
          std::isfinite(final_report.equality_residual) &&
          final_report.primal_violation <=
            active_policy.primal_feasibility_tolerance &&
          final_report.dual_violation <=
            active_policy.dual_feasibility_tolerance &&
          final_report.complementarity_residual <=
            active_policy.complementarity_tolerance &&
          final_report.stationarity_residual <=
            active_policy.stationarity_tolerance &&
          final_report.equality_residual <= active_policy.equality_tolerance &&
          final_report.primal_feasible && final_report.dual_feasible &&
          final_report.complementarity_converged &&
          final_report.kkt_residuals_converged &&
          final_report.active_set_stable && final_report.kkt_solve.converged() &&
          final_report.kkt_solve.linear_solve.converged(),
        "compiled DTO active PDAS diagnostics were incomplete");

      const contract::ActiveSetKKTSubproblemT<Backend> active_subproblem(
        active_pdas.product(),
        active_pdas.complementarity(),
        final_report.selection,
        1,
        active_policy.active_set_assumptions);
      const auto &active_product = active_subproblem.product();
      dealii::Vector<double> point_state(active_product.layout().primal->dimension(0));
      dealii::Vector<double> point_control(active_product.layout().primal->dimension(1));
      dealii::Vector<double> point_multiplier(active_product.layout().multiplier->dimension(0));
      dealii::Vector<double> seed_state(active_product.layout().stationarity->dimension(0));
      dealii::Vector<double> seed_control(active_product.layout().stationarity->dimension(1));
      dealii::Vector<double> seed_equality(active_product.layout().equality->dimension(0));
      for (std::size_t index = 0; index < point_state.size(); ++index)
        point_state[index] = 0.11 + 0.01 * static_cast<double>(index);
      for (std::size_t index = 0; index < point_control.size(); ++index)
        point_control[index] = -0.23 + 0.02 * static_cast<double>(index);
      for (std::size_t index = 0; index < point_multiplier.size(); ++index)
        point_multiplier[index] = 0.31 - 0.01 * static_cast<double>(index);
      for (std::size_t index = 0; index < seed_state.size(); ++index)
        seed_state[index] = -0.17 + 0.015 * static_cast<double>(index);
      for (std::size_t index = 0; index < seed_control.size(); ++index)
        seed_control[index] = 0.29 - 0.02 * static_cast<double>(index);
      for (std::size_t index = 0; index < seed_equality.size(); ++index)
        seed_equality[index] = 0.07 + 0.025 * static_cast<double>(index);
      const Product::Point mixed_point{
        Primal(active_product.layout().primal,
               {std::move(point_state), std::move(point_control)}),
        Primal(active_product.layout().multiplier,
               {std::move(point_multiplier)})};
      const Product::Seed mixed_seed{
        Primal(active_product.layout().stationarity,
               {std::move(seed_state), std::move(seed_control)}),
        Primal(active_product.layout().equality,
               {std::move(seed_equality)})};
      const auto mixed_action = active_product.apply_kkt(mixed_point);
      const auto mixed_transpose = active_product.apply_kkt_transpose(mixed_seed);
      const double transpose_left =
        contract::pair(mixed_action.stationarity, mixed_seed.stationarity) +
        contract::pair(mixed_action.equality, mixed_seed.equality);
      const double transpose_right =
        contract::pair(mixed_transpose.primal, mixed_point.primal) +
        contract::pair(mixed_transpose.multiplier, mixed_point.multiplier);
      require_close(transpose_left,
                    transpose_right,
                    1e-9,
                    "compiled DTO active PDAS action/transpose pairing");
    };
    run_active_case(contract::BoxActivity::lower);
    run_active_case(contract::BoxActivity::upper);

    using Solver = contract::PDASSolverT<Backend>;
    using BoxMultiplier = contract::CovectorBlockT<Backend>;
    struct DetachedCompiledPDAS
    {
      Product       product;
      Solver        solver;
      BoxMultiplier initial_box_multiplier;
      contract::PDASPolicy policy;
    };

    const auto detached = [&] {
      auto detached_mesh = std::make_unique<dealii::Triangulation<dim>>();
      dealii::GridGenerator::hyper_cube(*detached_mesh);
      detached_mesh->refine_global(1);
      const auto detached_session =
        std::make_shared<compiler::v1::DealiiCompilationSession<dim>>(
          std::move(detached_mesh), "test.compiled_pdas.owned_session");
      const dealii::Functions::ConstantFunction<dim> detached_forcing(1.0);
      const dealii::Functions::ConstantFunction<dim> detached_desired(0.25);
      const compiler::v1::DealiiDataBindings<dim> detached_bindings{
        detached_forcing,
        detached_desired,
        1.0,
        0.5,
        0.1,
        test_binding_provenance("compiled_pdas_detached")};
      const auto detached_specification =
        semantic::v1::make_scalar_diffusion_reaction_problem(true);
      const auto detached_compilation = compiler.compile(
        detached_specification,
        detached_session,
        detached_bindings,
        policy,
        bounds,
        std::nullopt,
        compiler::v1::CompilationProduct::pdas);
      contract::require(detached_compilation.succeeded() &&
                          detached_compilation.pdas_problem,
                        "compiled PDAS detachment setup failed");
      const auto compiled = detached_compilation.pdas_problem;
      const auto kkt_policy = compiled->kkt_solver_policy();
      const auto pdas_policy = compiled->pdas_policy();
      const auto product = compiled->product();
      const auto initial_box_multiplier =
        BoxMultiplier::zeros(compiled->complementarity().layout());
      const auto solver = compiled->make_solver(
        [kkt_policy](const Product &detached_product) {
          return dealii_backend::solve_serial_quadratic_kkt(
            detached_product, kkt_policy);
        });
      return DetachedCompiledPDAS{product,
                                  solver,
                                  initial_box_multiplier,
                                  pdas_policy};
    }();

    const auto detached_result = detached.solver.solve(
      Product::Point{Primal::zeros(detached.product.layout().primal),
                     Primal::zeros(detached.product.layout().multiplier)},
      detached.initial_box_multiplier,
      detached.policy);
    contract::require(
      detached_result.converged() &&
        !detached_result.iterations.empty(),
      "detached compiled PDAS solver did not retain its executable ownership");

    const auto supplied_specification =
      semantic::v1::make_scalar_diffusion_reaction_supplied_otd_problem(true);
    const auto supplied_validation = compiler.validate(
      supplied_specification,
      policy,
      compiler::v1::CompilationProduct::pdas);
    contract::require(supplied_validation.valid(),
                      "canonical supplied-OTD PDAS target did not validate");
    dealii::Triangulation<dim> supplied_triangulation;
    dealii::GridGenerator::hyper_cube(supplied_triangulation);
    supplied_triangulation.refine_global(1);
    const auto supplied_compilation = compiler.compile(
      supplied_specification,
      supplied_triangulation,
      data_bindings,
      policy,
      bounds,
      std::nullopt,
      compiler::v1::CompilationProduct::pdas);
    contract::require(
      supplied_compilation.succeeded() && supplied_compilation.pdas_problem &&
        supplied_compilation.pdas_problem->manifest()
              .resolved_decision.formulation_record.provenance ==
          semantic::v1::FormulationProvenance::supplied_otd &&
        supplied_compilation.pdas_problem->manifest().resolved_decision.pdas_record.product_id ==
          "compiled.scalar.supplied_otd.pdas",
      "compiler did not preserve the supplied-OTD route through the PDAS product");

    auto unsupported_policy = policy;
    unsupported_policy.pdas.active_set_assumptions.rank_condition_declared =
      false;
    const auto rejected_assumptions = compiler.validate(
      specification,
      unsupported_policy,
      compiler::v1::CompilationProduct::pdas);
    test_support::require_exact_diagnostic(
      rejected_assumptions,
      semantic::v1::DiagnosticCategory::formulation_capability,
      specification.formulation.id,
      "compiled_pdas_active_row_rank",
      "compiler accepted invalid PDAS active-row assumptions");

    const auto continuous_specification =
      semantic::v1::make_h1_regularised_scalar_diffusion_reaction_problem();
    const auto rejected_continuous = compiler.validate(
      continuous_specification,
      policy,
      compiler::v1::CompilationProduct::pdas);
    test_support::require_exact_diagnostic(
      rejected_continuous,
      semantic::v1::DiagnosticCategory::formulation_capability,
      continuous_specification.formulation.id,
      "compiled_pdas_cellwise_box",
      "compiler inferred a PDAS product for continuous control");

    auto mismatched_metric = specification;
    const auto selected_metric = std::find_if(
      mismatched_metric.metrics.begin(),
      mismatched_metric.metrics.end(),
      [&mismatched_metric](const semantic::v1::MetricSpec &candidate) {
        return candidate.id == mismatched_metric.formulation.metric_id;
      });
    contract::require(selected_metric != mismatched_metric.metrics.end(),
                      "compiled PDAS test could not find its selected metric");
    selected_metric->kind = semantic::v1::MetricKind::h1;
    const auto rejected_metric = compiler.validate(
      mismatched_metric,
      policy,
      compiler::v1::CompilationProduct::pdas);
    test_support::require_exact_diagnostic(
      rejected_metric,
      semantic::v1::DiagnosticCategory::formulation_capability,
      mismatched_metric.formulation.metric_id,
      "compiled_pdas_metric_realisation",
      "compiler accepted a PDAS metric without the selected multiplier realization");
  }

  template <int dim>
  void
  run_canonical_volume_control_contract_test()
  {
    verify_homogeneous_weak_form_oracle<dim>();

    dealii::Triangulation<dim> triangulation;
    dealii::GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(2);

    const dealii::Functions::ConstantFunction<dim> forcing(1.0);
    const dealii::Functions::ConstantFunction<dim> desired_state(0.25);
    const auto specification =
      semantic::v1::make_scalar_diffusion_reaction_problem(true);
    compiler::v1::DealiiDiscretisationPolicy compilation_policy;
    compilation_policy.state_degree = 1;
    compilation_policy.control_metric_solve = {1000, 1e-12, 1e-14};
    const compiler::v1::DealiiCompiler v1_compiler;
    const compiler::v1::DealiiDataBindings<dim> data_bindings{
      forcing,
      desired_state,
      1.0,
      0.5,
      0.1,
      test_binding_provenance("canonical_volume")};
    const compiler::v1::CellwiseBoxDataBindings bound_bindings{
      compiler::v1::CellwiseBoundValue{-1.0},
      compiler::v1::CellwiseBoundValue{0.05}};
    const auto compilation = v1_compiler.compile(specification,
                                                  triangulation,
                                                  data_bindings,
                                                  compilation_policy,
                                                  bound_bindings);
    contract::require(compilation.succeeded() && compilation.problem,
                      "v1 compiler failed to produce the canonical component");
    const auto &model = compilation.problem->executable_model();

    contract::StateControlPartitionT<Backend> partition(model, 0, 1);
    const auto reduced = compilation.problem->make_reduced_dto();

    dealii::Vector<double> control_values(
      partition.control_layout()->dimension(0));
    for (dealii::types::global_dof_index i = 0; i < control_values.size(); ++i)
      control_values[i] = 0.1 + 0.01 * static_cast<double>(i);
    const Primal control(partition.control_layout(), {std::move(control_values)});
    const auto evaluation = reduced.evaluate(control);

    const Covector state_residual = model.residual(evaluation.full_point);
    require_close(state_residual.block(0).l2_norm(),
                  0.0,
                  1e-11,
                  "deal.II state residual");

    const auto &metric = compilation.problem->metric();
    dealii::Vector<double> random_covector_values(
      partition.control_layout()->dimension(0));
    for (dealii::types::global_dof_index i = 0;
         i < random_covector_values.size();
         ++i)
      random_covector_values[i] =
        0.17 + 0.03 * static_cast<double>((7 * i) % 11);
    const Covector random_covector(partition.control_layout(),
                                   {std::move(random_covector_values)});
    const Primal recovered_primal = metric.inverse_apply(random_covector);
    const Covector recovered_covector = metric.apply(recovered_primal);
    dealii::Vector<double> inverse_apply_error = recovered_covector.block(0);
    inverse_apply_error.add(-1.0, random_covector.block(0));
    require_close(inverse_apply_error.l2_norm(),
                  0.0,
                  1e-11,
                  "deal.II mass metric inverse/apply consistency");
    dealii::Vector<double> nonidentity_mass_action = recovered_primal.block(0);
    nonidentity_mass_action.add(-1.0, random_covector.block(0));
    contract::require(nonidentity_mass_action.l2_norm() > 1e-3,
                      "deal.II control mass matrix unexpectedly acts as identity");

    dealii::Vector<double> state_tangent(
      model.variable_layout()->dimension(0));
    dealii::Vector<double> control_tangent(
      model.variable_layout()->dimension(1));
    for (dealii::types::global_dof_index i = 0; i < state_tangent.size(); ++i)
      state_tangent[i] = 0.02 * static_cast<double>(i + 1);
    for (dealii::types::global_dof_index i = 0; i < control_tangent.size();
         ++i)
      control_tangent[i] = -0.03 * static_cast<double>(i + 1);
    const Primal tangent(model.variable_layout(),
                         {std::move(state_tangent), std::move(control_tangent)});

    dealii::Vector<double> seed_values(model.test_layout()->dimension(0));
    for (dealii::types::global_dof_index i = 0; i < seed_values.size(); ++i)
      seed_values[i] = 0.04 * static_cast<double>(i + 1);
    const Primal test_seed(model.test_layout(), {std::move(seed_values)});

    const Covector jvp = model.residual_jvp(evaluation.full_point, tangent);
    const Covector vjp = model.residual_vjp(evaluation.full_point, test_seed);
    require_close(contract::pair(jvp, test_seed),
                  contract::pair(vjp, tangent),
                  1e-11,
                  "deal.II residual JVP/VJP pairing");

    constexpr double epsilon = 1e-7;
    const Covector residual_difference =
      model.residual(shifted(evaluation.full_point, tangent, epsilon));
    const Covector residual_at_point = model.residual(evaluation.full_point);
    for (std::size_t block = 0; block < residual_difference.n_blocks(); ++block)
      {
        dealii::Vector<double> finite_difference =
          residual_difference.block(block);
        Backend::add_scaled(finite_difference,
                            -1.0,
                            residual_at_point.block(block));
        Backend::scale(finite_difference, 1.0 / epsilon);
        finite_difference.add(-1.0, jvp.block(block));
        require_close(finite_difference.l2_norm(),
                      0.0,
                      1e-8,
                      "deal.II residual finite-difference JVP");
      }

    dealii::Vector<double> control_direction_values(
      partition.control_layout()->dimension(0));
    for (dealii::types::global_dof_index i = 0;
         i < control_direction_values.size();
         ++i)
      control_direction_values[i] =
        (i % 2 == 0 ? 0.05 : -0.04) * static_cast<double>(i + 1);
    const Primal control_direction(partition.control_layout(),
                                   {std::move(control_direction_values)});

    const Primal metric_direction =
      reduced.gradient_direction(evaluation.reduced_derivative, metric);
    const Covector metric_covector = metric.apply(metric_direction);
    require_close(contract::pair(metric_covector, control_direction),
                  contract::pair(evaluation.reduced_derivative,
                                 control_direction),
                  1e-11,
                  "deal.II mass metric pairing");

    const double reduced_difference =
      reduced.evaluate(shifted(control, control_direction, epsilon)).objective_value -
      evaluation.objective_value;
    require_close(reduced_difference / epsilon,
                  contract::pair(evaluation.reduced_derivative,
                                 control_direction),
                  2e-7,
                  "deal.II reduced DTO derivative");

    const auto &hessian = *compilation.problem->reduced_hessian();
    const Covector hessian_action =
      hessian.apply(control, control_direction);
    const Primal second_control_direction = [&partition]() {
      dealii::Vector<double> values(partition.control_layout()->dimension(0));
      for (dealii::types::global_dof_index i = 0; i < values.size(); ++i)
        values[i] = (i % 3 == 0 ? -0.02 : 0.03) *
                    static_cast<double>(i + 1);
      return Primal(partition.control_layout(), {std::move(values)});
    }();
    const Covector second_hessian_action =
      hessian.apply(control, second_control_direction);
    require_close(contract::pair(hessian_action, second_control_direction),
                  contract::pair(second_hessian_action, control_direction),
                  1e-10,
                  "deal.II reduced Hessian symmetry");

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
                  "deal.II reduced Hessian finite-difference action");

    nmopt::solvers::ReducedNewtonParameters newton_parameters;
    newton_parameters.maximum_inner_iterations = 100;
    newton_parameters.relative_tolerance = 1e-8;
    newton_parameters.absolute_tolerance = 1e-10;
    const nmopt::solvers::NewtonDirectionPolicyT<Backend> newton_direction(
      hessian, newton_parameters);
    nmopt::solvers::ReducedSolverParameters newton_solver_parameters;
    newton_solver_parameters.maximum_iterations = 20;
    newton_solver_parameters.maximum_line_search_trials = 30;
    newton_solver_parameters.gradient_tolerance = 1e-6;
    newton_solver_parameters.initial_step_length = 1.0;
    const nmopt::solvers::ReducedNewtonSolverT<Backend> newton_solver(
      reduced, metric, newton_solver_parameters, newton_direction);
    const auto newton_result = newton_solver.solve(control);
    contract::require(
      newton_result.stopping_reason ==
        nmopt::solvers::ReducedGradientStoppingReason::gradient_tolerance,
      "deal.II reduced Newton solver did not reach its tolerance");
    contract::require(newton_result.hessian_action_count > 0,
                      "deal.II reduced Newton solver did not use its Hessian");
    contract::require(newton_result.state_solve_count ==
                        newton_result.line_search_trial_count + 1,
                      "deal.II reduced Newton solve count misses a trial evaluation");

    nmopt::solvers::ReducedGradientParameters solver_parameters;
    solver_parameters.maximum_iterations = 100;
    solver_parameters.maximum_line_search_trials = 30;
    solver_parameters.gradient_tolerance = 1e-6;
    solver_parameters.initial_step_length = 20.0;
    solver_parameters.armijo_fraction = 1e-4;
    solver_parameters.backtracking_factor = 0.5;
    const nmopt::solvers::ReducedGradientSolverT<Backend> solver(
      reduced, metric, solver_parameters);
    const auto solver_result = solver.solve(control);

    contract::require(
      solver_result.stopping_reason ==
        nmopt::solvers::ReducedGradientStoppingReason::gradient_tolerance,
      "deal.II reduced gradient solver did not reach its tolerance");
    contract::require(solver_result.objective_history.size() > 1,
                      "deal.II reduced gradient solver did not accept an iteration");
    for (std::size_t index = 1;
         index < solver_result.objective_history.size();
         ++index)
      contract::require(solver_result.objective_history[index] <=
                          solver_result.objective_history[index - 1],
                        "deal.II reduced gradient objective history is not monotonic");
    contract::require(solver_result.gradient_norm_history.back() <=
                        solver_parameters.gradient_tolerance,
                      "deal.II reduced gradient final norm exceeds tolerance");
    contract::require(solver_result.adjoint_solve_count ==
                        solver_result.accepted_iterations + 1,
                      "deal.II reduced gradient adjoint count includes rejected trials");
    contract::require(solver_result.line_search_trial_count + 1 ==
                        solver_result.state_solve_count,
                      "deal.II reduced gradient solve count misses a trial evaluation");
    contract::require(solver_result.metric_solve_count ==
                        solver_result.gradient_norm_history.size(),
                      "deal.II reduced gradient metric solve count does not match direction evaluations");
    contract::require(solver_result.step_length_history.size() ==
                        solver_result.accepted_iterations,
                      "deal.II reduced gradient step history does not match accepted iterations");
    contract::require(solver_result.objective_change_history.size() ==
                        solver_result.accepted_iterations,
                      "deal.II reduced gradient objective-change history does not match accepted iterations");

    const nmopt::solvers::ReducedLimitedMemoryBfgsSolverT<Backend>
      lbfgs_solver(reduced, metric, solver_parameters);
    const auto lbfgs_result = lbfgs_solver.solve(control);
    contract::require(
      lbfgs_result.stopping_reason ==
        nmopt::solvers::ReducedGradientStoppingReason::gradient_tolerance,
      "deal.II L-BFGS solver did not reach its tolerance");
    contract::require(lbfgs_result.step_length_history.size() ==
                        lbfgs_result.accepted_iterations,
                      "deal.II L-BFGS step history does not match accepted iterations");
    contract::require(lbfgs_result.direction_reset_count <=
                        lbfgs_result.accepted_iterations,
                      "deal.II L-BFGS direction reset count exceeds accepted iterations");
    for (std::size_t index = 1;
         index < lbfgs_result.objective_history.size();
         ++index)
      contract::require(lbfgs_result.objective_history[index] <=
                          lbfgs_result.objective_history[index - 1],
                        "deal.II L-BFGS objective history is not monotonic");

    const auto *const bounds = compilation.problem->constraint();
    contract::require(bounds != nullptr,
                      "the canonical component lost its box constraint");
    dealii::Vector<double> bounded_control_values(
      partition.control_layout()->dimension(0));
    const Primal bounded_control(partition.control_layout(),
                                 {std::move(bounded_control_values)});
    const nmopt::solvers::ReducedGradientSolverT<Backend> projected_solver(
      reduced, metric, *bounds, solver_parameters);
    const auto projected_result = projected_solver.solve(bounded_control);

    contract::require(
      projected_result.stopping_reason ==
        nmopt::solvers::ReducedGradientStoppingReason::gradient_tolerance,
      "deal.II projected reduced gradient solver did not reach stationarity");
    contract::require(bounds->is_feasible(projected_result.control),
                      "deal.II projected reduced gradient returned an infeasible control");
    contract::require(projected_result.gradient_norm_history.back() <=
                        solver_parameters.gradient_tolerance,
                      "deal.II projected reduced gradient final norm exceeds tolerance");
    contract::require(projected_result.metric_solve_count ==
                        projected_result.gradient_norm_history.size(),
                      "deal.II projected reduced gradient metric solve count does not match direction evaluations");
    contract::require(projected_result.step_length_history.size() ==
                        projected_result.accepted_iterations,
                      "deal.II projected reduced gradient step history does not match accepted iterations");
    bool upper_bound_is_active = false;
    for (dealii::types::global_dof_index index = 0;
         index < projected_result.control.block(0).size();
         ++index)
      {
        contract::require(projected_result.control.block(0)[index] >= -1.0 &&
                            projected_result.control.block(0)[index] <= 0.05,
                          "deal.II projected reduced gradient left the box");
        upper_bound_is_active = upper_bound_is_active ||
                                projected_result.control.block(0)[index] >=
                                  0.05 - 1e-12;
      }
    contract::require(upper_bound_is_active,
                      "deal.II projected reduced gradient did not reach a bound");
    for (std::size_t index = 1;
         index < projected_result.objective_history.size();
         ++index)
      contract::require(projected_result.objective_history[index] <=
                          projected_result.objective_history[index - 1],
                        "deal.II projected reduced objective is not monotonic");

    const auto validation = v1_compiler.validate(specification,
                                                 compilation_policy);
    contract::require(validation.valid(),
                      "the canonical v1 problem did not validate for deal.II");

    compiler::v1::DealiiDiscretisationPolicy unsupported_execution =
      compilation_policy;
    unsupported_execution.execution =
      compiler::v1::DealiiDiscretisationPolicy::Execution::matrix_free;
    const auto lowerability_diagnostic =
      v1_compiler.validate(specification, unsupported_execution);
    test_support::require_exact_diagnostic(
      lowerability_diagnostic,
      semantic::v1::DiagnosticCategory::lowerability,
      "scalar_diffusion_reaction_volume_control",
      "assembled_execution",
      "v1 compiler did not report an unsupported execution mode");

    auto unsupported_formulation = specification;
    unsupported_formulation.formulation.kind =
      semantic::v1::FormulationKind::all_at_once;
    const auto formulation_diagnostic =
      v1_compiler.validate(unsupported_formulation, compilation_policy);
    test_support::require_exact_diagnostic(
      formulation_diagnostic,
      semantic::v1::DiagnosticCategory::formulation_capability,
      "reduced_dto",
      "reduced_dto_formulation",
      "v1 compiler did not report an unsupported formulation capability");

    auto supplied_otd_specification =
      semantic::v1::make_scalar_diffusion_reaction_supplied_otd_problem(false);
    const auto supplied_otd_validation =
      v1_compiler.validate(supplied_otd_specification, compilation_policy);
    contract::require(supplied_otd_validation.valid(),
                      "v1 compiler rejected the canonical supplied OTD target");
    auto mismatched_supplied_otd = supplied_otd_specification;
    mismatched_supplied_otd.spaces.push_back(
      {"alternate_state_test_space", "Alternate adjoint test", "domain",
       semantic::v1::SpaceTopology::h1, semantic::v1::SpaceRole::test});
    mismatched_supplied_otd.pairings.push_back(
      {"alternate_state_test_pairing", "Alternate adjoint pairing",
       "alternate_state_test_space", "alternate_state_test_space"});
    mismatched_supplied_otd.supplied_otd_declaration->adjoint_block
      .variable_space_id = "alternate_state_test_space";
    mismatched_supplied_otd.supplied_otd_declaration->adjoint_block
      .residual_space_id = "alternate_state_test_space";
    mismatched_supplied_otd.supplied_otd_declaration->adjoint_block
      .trial_pairing_id = "alternate_state_test_pairing";
    mismatched_supplied_otd.supplied_otd_declaration->adjoint_block
      .test_pairing_id = "alternate_state_test_pairing";
    const auto mismatched_otd_diagnostic =
      v1_compiler.validate(mismatched_supplied_otd, compilation_policy);
    test_support::require_exact_diagnostic(
      mismatched_otd_diagnostic,
      semantic::v1::DiagnosticCategory::formulation_capability,
      "adjoint_block",
      "supplied_otd_adjoint_space",
      "v1 compiler reported a mismatched supplied OTD adjoint space as valid");

    const auto supplied_otd_compilation = v1_compiler.compile(
      supplied_otd_specification,
      triangulation,
      data_bindings,
      compilation_policy);
    contract::require(
      supplied_otd_compilation.succeeded() &&
        !supplied_otd_compilation.problem &&
        supplied_otd_compilation.supplied_otd_problem,
      "v1 compiler did not produce the distinct supplied OTD product");
    const auto &supplied_otd = *supplied_otd_compilation.supplied_otd_problem;
    const auto supplied_initial =
      Primal::zeros(supplied_otd.system().variable_layout());
    const auto supplied_solution = supplied_otd.system().solve(supplied_initial);
    contract::require(supplied_solution.report.converged() &&
                        supplied_solution.report.algorithm ==
                          "serial_sparse_direct_umfpack",
                      "serial supplied OTD solve did not report direct convergence");
    const auto supplied_residual =
      supplied_otd.system().residual(supplied_solution.solution);
    for (std::size_t block = 0; block < supplied_residual.n_blocks(); ++block)
      require_close(supplied_residual.block(block).l2_norm(),
                    0.0,
                    1e-11,
                    "serial supplied OTD solution leaves a residual");

    const Primal supplied_control(
      partition.control_layout(),
      {supplied_solution.solution.block(
        supplied_otd.system().block_selection().control_variable)});
    const auto dto_at_supplied_control = reduced.evaluate(supplied_control);
    dealii::Vector<double> state_difference =
      supplied_solution.solution.block(
        supplied_otd.system().block_selection().state_variable);
    state_difference.add(-1.0, dto_at_supplied_control.state.block(0));
    require_close(state_difference.l2_norm(),
                  0.0,
                  1e-11,
                  "serial supplied OTD state differs from reduced DTO");
    dealii::Vector<double> adjoint_difference =
      supplied_solution.solution.block(
        supplied_otd.system().block_selection().adjoint_variable);
    adjoint_difference.add(-1.0, dto_at_supplied_control.adjoint.block(0));
    require_close(adjoint_difference.l2_norm(),
                  0.0,
                  1e-11,
                  "serial supplied OTD adjoint differs from reduced DTO");
    dealii::Vector<double> stationarity_difference =
      supplied_otd.system().control_stationarity(supplied_solution.solution)
        .block(0);
    stationarity_difference.add(
      -1.0, dto_at_supplied_control.reduced_derivative.block(0));
    require_close(stationarity_difference.l2_norm(),
                  0.0,
                  1e-11,
                  "serial supplied OTD stationarity differs from reduced DTO");
    const auto &supplied_system = supplied_otd.system();
    dealii::Vector<double> adjoint_tangent(
      model.variable_layout()->dimension(0));
    for (dealii::types::global_dof_index index = 0;
         index < adjoint_tangent.size();
         ++index)
      adjoint_tangent[index] = -0.017 * static_cast<double>(index + 1);
    const Primal supplied_tangent(
      supplied_system.variable_layout(),
      {tangent.block(0), std::move(adjoint_tangent), tangent.block(1)});
    const Primal supplied_comparison_point =
      shifted(supplied_solution.solution, supplied_tangent, 0.37);
    contract::require(
      supplied_otd.manifest().resolved_decision.supplied_otd_record.declaration.has_value() &&
        supplied_otd.manifest().resolved_decision.supplied_otd_record.declaration
            ->multiplier_conversion ==
          semantic::v1::SuppliedOTDMultiplierConversion::identity,
      "serial supplied OTD DTO comparison used an undeclared multiplier conversion");
    const Covector supplied_comparison_residual =
      supplied_system.residual(supplied_comparison_point);
    for (std::size_t block = 0;
         block < supplied_comparison_residual.n_blocks();
         ++block)
      contract::require(
        supplied_comparison_residual.block(block).l2_norm() > 1e-10,
        "serial supplied OTD manufactured point did not activate every block");

    constexpr double supplied_derivative_step = 1e-7;
    const Covector supplied_jvp =
      supplied_system.residual_jvp(supplied_comparison_point,
                                   supplied_tangent);
    const Covector supplied_comparison_residual_plus = supplied_system.residual(
      shifted(supplied_comparison_point,
              supplied_tangent,
              supplied_derivative_step));
    const Covector supplied_comparison_residual_minus = supplied_system.residual(
      shifted(supplied_comparison_point,
              supplied_tangent,
              -supplied_derivative_step));
    for (std::size_t block = 0;
         block < supplied_comparison_residual.n_blocks();
         ++block)
      {
        dealii::Vector<double> finite_difference =
          supplied_comparison_residual_plus.block(block);
        finite_difference.add(-1.0,
                              supplied_comparison_residual_minus.block(block));
        finite_difference *= 0.5 / supplied_derivative_step;
        finite_difference.add(-1.0, supplied_jvp.block(block));
        require_close(finite_difference.l2_norm(),
                      0.0,
                      1e-7,
                      "serial supplied OTD JVP finite difference");
      }

    dealii::Vector<double> seed_state(model.variable_layout()->dimension(0));
    dealii::Vector<double> seed_adjoint(model.variable_layout()->dimension(0));
    dealii::Vector<double> seed_control(model.variable_layout()->dimension(1));
    for (dealii::types::global_dof_index index = 0;
         index < seed_state.size();
         ++index)
      {
        seed_state[index] = 0.013 * static_cast<double>(index + 1);
        seed_adjoint[index] = -0.021 * static_cast<double>(index + 1);
      }
    for (dealii::types::global_dof_index index = 0;
         index < seed_control.size();
         ++index)
      seed_control[index] = 0.031 * static_cast<double>(index + 1);
    const Primal supplied_residual_seed(
      supplied_system.residual_layout(),
      {std::move(seed_state),
       std::move(seed_adjoint),
       std::move(seed_control)});
    const Covector supplied_vjp = supplied_system.residual_vjp(
      supplied_comparison_point, supplied_residual_seed);
    require_close(contract::pair(supplied_jvp, supplied_residual_seed),
                  contract::pair(supplied_vjp, supplied_tangent),
                  1e-10,
                  "serial supplied OTD JVP/VJP pairing");

    const Primal dto_comparison_point(
      model.variable_layout(),
      {supplied_comparison_point.block(0),
       supplied_comparison_point.block(2)});
    const Primal dto_multiplier(model.test_layout(),
                                {supplied_comparison_point.block(1)});
    const Covector dto_state_residual = model.residual(dto_comparison_point);
    const Covector dto_objective_derivative =
      model.objective_derivative(dto_comparison_point);
    const Covector dto_residual_pullback =
      model.residual_vjp(dto_comparison_point, dto_multiplier);
    dealii::Vector<double> expected_adjoint = dto_residual_pullback.block(0);
    expected_adjoint.add(-1.0, dto_objective_derivative.block(0));
    dealii::Vector<double> expected_stationarity =
      dto_objective_derivative.block(1);
    expected_stationarity.add(-1.0, dto_residual_pullback.block(1));
    dealii::Vector<double> dto_state_difference =
      supplied_comparison_residual.block(0);
    dto_state_difference.add(-1.0, dto_state_residual.block(0));
    dealii::Vector<double> dto_adjoint_difference =
      supplied_comparison_residual.block(1);
    dto_adjoint_difference.add(-1.0, expected_adjoint);
    dealii::Vector<double> dto_stationarity_difference =
      supplied_comparison_residual.block(2);
    dto_stationarity_difference.add(-1.0, expected_stationarity);
    require_close(dto_state_difference.l2_norm(),
                  0.0,
                  1e-11,
                  "serial supplied OTD state action differs from DTO");
    require_close(dto_adjoint_difference.l2_norm(),
                  0.0,
                  1e-11,
                  "serial supplied OTD adjoint action differs from DTO");
    require_close(dto_stationarity_difference.l2_norm(),
                  0.0,
                  1e-11,
                  "serial supplied OTD stationarity action differs from DTO");

    const auto &supplied_manifest = supplied_otd.manifest();
    const std::string expected_comparison_status =
      "equivalence verified under declared conversion: "
      "verified by scenario nmopt.dealii.canonical_volume_control: "
      "manufactured-point block, JVP finite-difference, transpose-pairing, "
      "and DTO-action comparisons";
    contract::require(
      supplied_manifest.resolved_decision.formulation_record.kind ==
          semantic::v1::FormulationKind::all_at_once &&
        supplied_manifest.resolved_decision.formulation_record.provenance ==
          semantic::v1::FormulationProvenance::supplied_otd &&
        supplied_manifest.compatibility.provenance == "supplied OTD" &&
        supplied_manifest.resolved_decision.supplied_otd_record.present &&
        supplied_manifest.resolved_decision.supplied_otd_record.declaration.has_value() &&
        supplied_manifest.resolved_decision.supplied_otd_record.declaration->id ==
          "canonical_scalar_supplied_otd" &&
        supplied_manifest.resolved_decision.supplied_otd_record.declaration->adjoint_block
            .test_pairing_id == "state_test_pairing" &&
        supplied_manifest.resolved_decision.supplied_otd_record.declaration
            ->multiplier_conversion ==
          semantic::v1::SuppliedOTDMultiplierConversion::identity &&
        supplied_manifest.resolved_decision.supplied_otd_record.variable_space_ids ==
          std::vector<std::string>{"state", "state_test", "control"} &&
        supplied_manifest.resolved_decision.supplied_otd_record.variable_dimensions ==
          std::vector<std::size_t>{model.variable_layout()->dimension(0),
                                   model.test_layout()->dimension(0),
                                   model.variable_layout()->dimension(1)} &&
        supplied_manifest.resolved_decision.supplied_otd_record.residual_space_ids ==
          std::vector<std::string>{"state_equation",
                                   "adjoint_equation",
                                   "control_stationarity"} &&
        supplied_manifest.resolved_decision.supplied_otd_record.residual_dimensions ==
          std::vector<std::size_t>{model.test_layout()->dimension(0),
                                   model.test_layout()->dimension(0),
                                   model.variable_layout()->dimension(1)} &&
        supplied_manifest.resolved_decision.supplied_otd_record.comparison_status ==
          expected_comparison_status,
      "serial supplied OTD manifest omitted formulation and comparison provenance");
    const Covector compiled_residual = model.residual(evaluation.full_point);
    require_close(compiled_residual.block(0).l2_norm(),
                  0.0,
                  1e-11,
                  "canonical scalar state residual");
    const auto *compiled_constraint = compilation.problem->constraint();
    contract::require(compiled_constraint != nullptr &&
                        compiled_constraint->is_feasible(bounded_control),
                      "v1 compiler did not preserve the declared box constraint");
    const auto &manifest = compilation.problem->manifest();
    require_resolved_manifest_projection(
      manifest, "canonical volume control manifest");
    dealii::Triangulation<dim> coarser_triangulation;
    dealii::GridGenerator::hyper_cube(coarser_triangulation);
    coarser_triangulation.refine_global(1);
    const auto coarser_compilation = v1_compiler.compile(
      specification,
      coarser_triangulation,
      data_bindings,
      compilation_policy,
      bound_bindings);
    contract::require(
      coarser_compilation.succeeded() &&
        coarser_compilation.problem->manifest().resolved_decision.mesh_record.structural_identity !=
          manifest.resolved_decision.mesh_record.structural_identity,
      "v1 mesh provenance did not distinguish different compiled structures");

    auto relabeled_specification = specification;
    relabeled_specification.label = "same typed scalar graph, different label";
    for (auto &region : relabeled_specification.regions)
      region.label += " (renamed)";
    for (auto &space : relabeled_specification.spaces)
      space.label += " (renamed)";
    const auto relabeled_compilation = v1_compiler.compile(
      relabeled_specification,
      triangulation,
      data_bindings,
      compilation_policy,
      bound_bindings);
    contract::require(relabeled_compilation.succeeded(),
                      "v1 compiler rejected a label-only graph change");
    require_resolved_manifest_projection(
      relabeled_compilation.problem->manifest(),
      "label-only graph change manifest");
    require_compiled_binding_records_equal(
      relabeled_compilation.problem->manifest().resolved_decision.bindings,
      manifest.resolved_decision.bindings,
      "label-only graph change");
    contract::require(
      relabeled_compilation.problem->manifest().resolved_decision.target_id ==
        manifest.resolved_decision.target_id &&
        relabeled_compilation.problem->manifest().resolved_decision.metric_record.realisation_id ==
          manifest.resolved_decision.metric_record.realisation_id,
      "label-only graph change altered typed execution records");

    auto reordered_specification = specification;
    std::reverse(reordered_specification.regions.begin(),
                 reordered_specification.regions.end());
    std::reverse(reordered_specification.spaces.begin(),
                 reordered_specification.spaces.end());
    std::reverse(reordered_specification.pairings.begin(),
                 reordered_specification.pairings.end());
    std::reverse(reordered_specification.variables.begin(),
                 reordered_specification.variables.end());
    std::reverse(reordered_specification.data.begin(),
                 reordered_specification.data.end());
    std::reverse(reordered_specification.transformations.begin(),
                 reordered_specification.transformations.end());
    std::reverse(reordered_specification.residual_terms.begin(),
                 reordered_specification.residual_terms.end());
    std::reverse(reordered_specification.observations.begin(),
                 reordered_specification.observations.end());
    std::reverse(reordered_specification.losses.begin(),
                 reordered_specification.losses.end());
    std::reverse(reordered_specification.metrics.begin(),
                 reordered_specification.metrics.end());
    std::reverse(reordered_specification.constraints.begin(),
                 reordered_specification.constraints.end());
    std::reverse(reordered_specification.requirement_policies.begin(),
                 reordered_specification.requirement_policies.end());
    const auto reordered_compilation = v1_compiler.compile(
      reordered_specification,
      triangulation,
      data_bindings,
      compilation_policy,
      bound_bindings);
    contract::require(reordered_compilation.succeeded(),
                      "v1 compiler rejected declaration-order permutation");
    const auto &reordered_manifest = reordered_compilation.problem->manifest();
    test_support::require_manifest_compatibility_equal(
      manifest,
      reordered_manifest,
      "declaration-order permutation");
    require_compiled_binding_records_equal(
      reordered_manifest.resolved_decision.bindings,
      manifest.resolved_decision.bindings,
      "declaration-order permutation");
    contract::require(
      reordered_manifest.resolved_decision.target_id ==
          manifest.resolved_decision.target_id &&
        reordered_manifest.resolved_decision.formulation_record.semantic_id ==
          manifest.resolved_decision.formulation_record.semantic_id &&
        reordered_manifest.resolved_decision.mesh_record.structural_identity ==
          manifest.resolved_decision.mesh_record.structural_identity &&
        reordered_manifest.resolved_decision.metric_record.realisation_id ==
          manifest.resolved_decision.metric_record.realisation_id &&
        reordered_manifest.resolved_decision.constraint_record.realisation_id ==
          manifest.resolved_decision.constraint_record.realisation_id &&
        reordered_manifest.resolved_decision.state_solve_record.maximum_iterations ==
          manifest.resolved_decision.state_solve_record.maximum_iterations &&
        reordered_manifest.resolved_decision.adjoint_solve_record.maximum_iterations ==
          manifest.resolved_decision.adjoint_solve_record.maximum_iterations,
      "declaration-order permutation changed the closed compiler selection");

    const compiler::v1::CellwiseBoxDataBindings changed_bound_bindings{
      compiler::v1::CellwiseBoundValue{-1.0},
      compiler::v1::CellwiseBoundValue{0.15}};
    const auto changed_bound_compilation = v1_compiler.compile(
      specification,
      triangulation,
      data_bindings,
      compilation_policy,
      changed_bound_bindings);
    contract::require(changed_bound_compilation.succeeded(),
                      "v1 compiler rejected changed bound data");
    const auto find_binding = [](const compiler::v1::CompilationManifest &record,
                                 const semantic::v1::DataRole role) {
      return std::find_if(
        record.resolved_decision.bindings.begin(),
        record.resolved_decision.bindings.end(),
        [role](const compiler::v1::CompiledBindingRecord &binding) {
          return binding.role == role;
        });
    };
    const auto original_upper =
      find_binding(manifest, semantic::v1::DataRole::upper_bound);
    const auto changed_upper = find_binding(
      changed_bound_compilation.problem->manifest(),
      semantic::v1::DataRole::upper_bound);
    contract::require(
      original_upper != manifest.resolved_decision.bindings.end() &&
        changed_upper != changed_bound_compilation.problem->manifest().resolved_decision.bindings.end() &&
        original_upper->scalar_value.has_value() &&
        changed_upper->scalar_value.has_value() &&
        *original_upper->scalar_value == 0.05 &&
        *changed_upper->scalar_value == 0.15 &&
        original_upper->value_digest != changed_upper->value_digest,
      "v1 manifest did not distinguish changed bound values");

    require_constraint_realisation(
      manifest,
      "FE_DGQ(0) coefficientwise l2_cellwise clipping",
      "canonical volume control");
    const auto diffusion_binding = std::find_if(
      manifest.resolved_decision.bindings.begin(),
      manifest.resolved_decision.bindings.end(),
      [](const compiler::v1::CompiledBindingRecord &record) {
        return record.role == semantic::v1::DataRole::diffusion;
      });
    const auto lower_binding = std::find_if(
      manifest.resolved_decision.bindings.begin(),
      manifest.resolved_decision.bindings.end(),
      [](const compiler::v1::CompiledBindingRecord &record) {
        return record.role == semantic::v1::DataRole::lower_bound;
      });
    const auto upper_binding = std::find_if(
      manifest.resolved_decision.bindings.begin(),
      manifest.resolved_decision.bindings.end(),
      [](const compiler::v1::CompiledBindingRecord &record) {
        return record.role == semantic::v1::DataRole::upper_bound;
      });
    contract::require(
      diffusion_binding != manifest.resolved_decision.bindings.end() &&
        diffusion_binding->field_shape ==
          compiler::v1::CompiledFieldShape::scalar_constant &&
        diffusion_binding->scalar_value.has_value() &&
        *diffusion_binding->scalar_value == 1.0 &&
        !diffusion_binding->value_digest.empty() &&
        lower_binding != manifest.resolved_decision.bindings.end() &&
        lower_binding->scalar_value.has_value() &&
        *lower_binding->scalar_value == -1.0 &&
        upper_binding != manifest.resolved_decision.bindings.end() &&
        upper_binding->scalar_value.has_value() &&
        *upper_binding->scalar_value == 0.05 &&
        manifest.resolved_decision.target_id.find("compiled_target:") == 0,
      "v1 compiler did not retain lossless scalar binding provenance");
    contract::require(manifest.resolved_decision.semantic_problem_id == specification.id &&
                        manifest.compatibility.provenance == "DTO" &&
                        manifest.compatibility.execution == "assembled",
                      "v1 compiler did not record its compilation manifest");
  }

