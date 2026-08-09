#include "nmopt/contract/reduced_dto.hpp"
#include "nmopt/compiler/v1/dealii_scalar_diffusion_reaction.hpp"
#include "nmopt/dealii/scalar_diffusion_reaction.hpp"
#include "nmopt/semantic/v1/problem_spec.hpp"
#include "nmopt/solvers/reduced_gradient.hpp"
#include "test_support/contract_errors.hpp"
#include "test_support/diagnostics.hpp"
#include "test_support/scenario_dispatch.hpp"

#include <deal.II/base/function_lib.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace
{
  using namespace nmopt;
  using Backend = dealii_backend::SerialBackend;
  using Primal = contract::PrimalBlockT<Backend>;
  using Covector = contract::CovectorBlockT<Backend>;

  template <int dim>
  class ConstantTensorCoefficient final
    : public dealii::TensorFunction<2, dim>
  {
  public:
    explicit ConstantTensorCoefficient(const dealii::Tensor<2, dim> value)
      : value_(value)
    {}

    dealii::Tensor<2, dim>
    value(const dealii::Point<dim> &) const override
    {
      return value_;
    }

  private:
    dealii::Tensor<2, dim> value_;
  };

  template <int dim>
  class ConstantVectorCoefficient final
    : public dealii::TensorFunction<1, dim>
  {
  public:
    explicit ConstantVectorCoefficient(const dealii::Tensor<1, dim> value)
      : value_(value)
    {}

    dealii::Tensor<1, dim>
    value(const dealii::Point<dim> &) const override
    {
      return value_;
    }

  private:
    dealii::Tensor<1, dim> value_;
  };

  compiler::v1::DealiiBindingProvenance
  test_binding_provenance(const std::string &target,
                          const bool         has_fixed_data = false)
  {
    return {"test." + target + ".forcing",
            "test." + target + ".desired_state",
            has_fixed_data ? "test." + target + ".fixed_dirichlet" : ""};
  }

  void
  require_close(const double actual,
                const double expected,
                const double tolerance,
                const std::string &description)
  {
    if (std::abs(actual - expected) > tolerance)
      throw contract::ContractError(description + ": expected " +
                                    std::to_string(expected) + ", got " +
                                    std::to_string(actual));
  }

  void
  require_covector_close(const Covector &   actual,
                         const Covector &   expected,
                         const double       tolerance,
                         const std::string &description)
  {
    contract::require_compatible(actual, expected,
                                 description + " has incompatible layouts");
    for (std::size_t block = 0; block < actual.n_blocks(); ++block)
      {
        dealii::Vector<double> difference = actual.block(block);
        difference.add(-1.0, expected.block(block));
        require_close(difference.l2_norm(), 0.0, tolerance, description);
      }
  }

  void
  require_primal_close(const Primal &     actual,
                       const Primal &     expected,
                       const double       tolerance,
                       const std::string &description)
  {
    contract::require_compatible(actual, expected,
                                 description + " has incompatible layouts");
    for (std::size_t block = 0; block < actual.n_blocks(); ++block)
      {
        dealii::Vector<double> difference = actual.block(block);
        difference.add(-1.0, expected.block(block));
        require_close(difference.l2_norm(), 0.0, tolerance, description);
      }
  }

  Primal
  shifted(Primal value, const Primal &direction, const double step)
  {
    contract::require_compatible(value,
                                 direction,
                                 "deal.II shift has incompatible layouts");
    for (std::size_t block = 0; block < value.n_blocks(); ++block)
      value.add_scaled_block(block, step, direction.block(block));
    return value;
  }

  void
  require_constraint_realisation(
    const compiler::v1::CompilationManifest &manifest,
    const std::string &                       expected,
    const std::string &                       target)
  {
    std::string structured_expected = "none";
    if (expected.find("l2_cellwise_parameter") != std::string::npos)
      structured_expected = "l2_cellwise_parameter";
    else if (expected.find("l2_cellwise") != std::string::npos)
      structured_expected = "l2_cellwise";
    else if (expected.find("l2_facewise") != std::string::npos)
      structured_expected = "l2_facewise";
    contract::require(
      manifest.constraint_realisation == expected &&
        manifest.constraint_record.realisation_id == structured_expected &&
        manifest.constraint_record.present == (structured_expected != "none") &&
        (structured_expected == "none" ||
         manifest.constraint_record.projection_metric_id == structured_expected),
      target + " manifest constraint realization: expected " + expected +
        ", got " + manifest.constraint_realisation);
    const auto has_runtime_role = [&manifest](const std::string &role) {
      return std::any_of(
        manifest.spaces.begin(),
        manifest.spaces.end(),
        [&role](const compiler::v1::CompiledSpaceRecord &space) {
          return space.runtime_role == role;
        });
    };
    contract::require(
      manifest.schema_version == 1 &&
        manifest.formulation_record.kind ==
          semantic::v1::FormulationKind::reduced_dto &&
        manifest.formulation_record.execution ==
          compiler::v1::ExecutionRealisation::assembled &&
        manifest.mesh_record.dimension > 0 &&
        manifest.mesh_record.active_cells > 0 &&
        !manifest.mesh_record.provenance.empty() &&
        has_runtime_role("state") &&
        has_runtime_role("test_and_adjoint") &&
        (has_runtime_role("decision_control") ||
         has_runtime_role("decision_parameter")) &&
        has_runtime_role("observation") &&
        !manifest.bindings.empty() &&
        manifest.state_solve_record.maximum_iterations > 0 &&
        manifest.adjoint_solve_record.maximum_iterations > 0 &&
        !manifest.metric_record.semantic_id.empty() &&
        !manifest.metric_record.realisation_id.empty(),
      target + " structured manifest is incomplete");
  }

  void
  run_projection_compatibility_contract_test()
  {
    const auto layout = std::make_shared<const contract::BlockLayout>(
      "projection_compatibility",
      std::vector<contract::SpaceId>{{"control"}},
      std::vector<std::size_t>{2});
    dealii::DynamicSparsityPattern dynamic_pattern(2, 2);
    for (std::size_t row = 0; row < 2; ++row)
      for (std::size_t column = 0; column < 2; ++column)
        dynamic_pattern.add(row, column);
    dealii::SparsityPattern pattern;
    pattern.copy_from(dynamic_pattern);

    auto diagonal_matrix =
      std::make_shared<dealii::SparseMatrix<double>>(pattern);
    diagonal_matrix->set(0, 0, 2.0);
    diagonal_matrix->set(1, 1, 2.0);
    const dealii_backend::MassMetric cellwise_metric(
      "l2_cellwise", layout, diagonal_matrix);
    const dealii_backend::CellwiseBoxConstraint cellwise_box(
      layout, 0.0, 1.0, cellwise_metric);
    contract::require(cellwise_box.supports_projection_in(cellwise_metric),
                      "Cellwise box rejected its coupled diagonal metric");

    auto non_diagonal_matrix =
      std::make_shared<dealii::SparseMatrix<double>>(pattern);
    non_diagonal_matrix->set(0, 0, 2.0);
    non_diagonal_matrix->set(0, 1, 1.0);
    non_diagonal_matrix->set(1, 0, 1.0);
    non_diagonal_matrix->set(1, 1, 2.0);
    const dealii_backend::MassMetric spoofed_cellwise_metric(
      "l2_cellwise", layout, non_diagonal_matrix);
    contract::require(
      !cellwise_box.supports_projection_in(spoofed_cellwise_metric),
      "A non-diagonal deal.II metric obtained cellwise clipping by reusing the l2_cellwise display identifier");
    test_support::require_contract_error(
      [&layout, &spoofed_cellwise_metric]() {
        (void)dealii_backend::CellwiseBoxConstraint(
          layout, 0.0, 1.0, spoofed_cellwise_metric);
      },
      "Cellwise box projection needs a positive diagonal metric realization",
      "non-diagonal cellwise projection coupling");

    const dealii_backend::MassMetric facewise_metric(
      "l2_facewise", layout, diagonal_matrix);
    const dealii_backend::FacewiseBoxConstraint facewise_box(
      layout, 0.0, 1.0, facewise_metric);
    contract::require(facewise_box.supports_projection_in(facewise_metric),
                      "Facewise box rejected its coupled diagonal metric");
    const dealii_backend::MassMetric spoofed_facewise_metric(
      "l2_facewise", layout, non_diagonal_matrix);
    contract::require(
      !facewise_box.supports_projection_in(spoofed_facewise_metric),
      "A non-diagonal deal.II metric obtained facewise clipping by reusing the l2_facewise display identifier");
    const dealii_backend::MassMetric h1_metric(
      "h1_continuous", layout, non_diagonal_matrix);
    contract::require(!cellwise_box.supports_projection_in(h1_metric),
                      "Cellwise clipping accepted an H1 metric realization");
  }

  void
  run_backend_size_conversion_contract_test()
  {
    const std::size_t maximum =
      dealii_backend::SerialBackend::maximum_native_size();
    contract::require(
      static_cast<std::size_t>(
        dealii_backend::SerialBackend::checked_native_size(maximum)) == maximum,
      "Serial deal.II size conversion rejected its native maximum");

    if constexpr (dealii_backend::SerialBackend::native_size_is_narrower)
      test_support::require_contract_error(
        [maximum]() {
          (void)dealii_backend::SerialBackend::checked_native_size(maximum + 1);
        },
        "Serial deal.II vector size exceeds its native range",
        "oversized serial deal.II vector dimension");
  }

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
    require_constraint_realisation(manifest, "none", "fixed-Dirichlet");
    contract::require(
        manifest.lifting_realisation.find("y_phys = P_h y_hat + ell_0,h") !=
        std::string::npos &&
        manifest.data_rule.find("boundary DoFs") != std::string::npos &&
        manifest.transformation_ids.size() == 1 &&
        manifest.lowering_handler_records.size() == 9 &&
        std::find(manifest.lowering_handler_records.begin(),
                  manifest.lowering_handler_records.end(),
                  "fixed_dirichlet_reconstruction <- "
                  "dealii.scalar.transformation.fixed_dirichlet") !=
          manifest.lowering_handler_records.end(),
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
    require_constraint_realisation(manifest, "none", "Dirichlet-control");
    contract::require(
      manifest.control_space.find("nodal trace") != std::string::npos &&
        manifest.lifting_realisation.find("L_D,h") != std::string::npos &&
        manifest.metric_solve_policy.find("l2_dirichlet_trace") !=
          std::string::npos &&
        manifest.declared_assumptions.front().find("dirichlet_control_lifting") !=
          std::string::npos,
      "Dirichlet-control compilation manifest is incomplete");
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
      one_manifest.observation_realisation ==
        "material-id volume restriction: 1" &&
        one_manifest.data_rule.find("analytic desired-state Function") !=
          std::string::npos &&
        one_manifest.lowering_handler_records.size() == 8 &&
        std::find(one_manifest.lowering_handler_records.begin(),
                  one_manifest.lowering_handler_records.end(),
                  "diffusion_reaction <- "
                  "dealii.scalar.residual.diffusion_reaction") !=
          one_manifest.lowering_handler_records.end(),
      "v1 subdomain observation manifest omitted its restriction or data rule");
  }

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

    const auto &metric = compilation.problem->metric();
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
      manifest.control_space.find("facewise-constant") != std::string::npos &&
        manifest.observation_realisation.find("boundary trace") !=
          std::string::npos &&
        manifest.data_rule.find("boundary face quadrature") != std::string::npos &&
        manifest.constraint_realisation.find("l2_facewise") != std::string::npos,
      "Neumann boundary compilation manifest is incomplete");
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
      manifest.control_space.find("continuous scalar FE_Q") != std::string::npos &&
        manifest.declared_assumptions.front().find("h1_control_regularisation") !=
          std::string::npos,
      "H1-control compilation manifest omitted the loss or control realization");
    contract::require(
      h1_metric_compilation.problem->manifest().metric_solve_policy.find(
        "h1_continuous") != std::string::npos,
      "H1 metric compilation manifest omitted the selected Riesz map");
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
      compilation.problem->manifest().state_adjoint_solve_policy.find(
        "reassembled") != std::string::npos,
      "coefficient-identification manifest omitted state-matrix reassembly");
    require_constraint_realisation(
      compilation.problem->manifest(),
      "FE_DGQ(0) coefficientwise l2_cellwise_parameter clipping",
      "coefficient-identification");
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
      manifest.nullspace_policy.find("mean-zero Lagrange multiplier") !=
        std::string::npos &&
        manifest.state_adjoint_solve_policy.find("SparseDirectUMFPACK") !=
          std::string::npos &&
        manifest.lifting_realisation.find("pure-Neumann") != std::string::npos,
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
    contract::require(
      manifest.state_solve_record.algorithm ==
          compiler::v1::LinearSolveAlgorithm::serial_sparse_direct_umfpack &&
        manifest.adjoint_solve_record.algorithm ==
          compiler::v1::LinearSolveAlgorithm::serial_sparse_direct_umfpack &&
        manifest.lowering_handler_records.size() == 13 &&
        manifest.data_rule.find("Robin coefficient and source") !=
          std::string::npos &&
        std::find(manifest.lowering_handler_records.begin(),
                  manifest.lowering_handler_records.end(),
                  "conservative_transport <- "
                  "dealii.scalar.residual.conservative_transport") !=
          manifest.lowering_handler_records.end(),
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
      detached.manifest.schema_version == 1 &&
        detached.manifest.mesh_record.dimension ==
          static_cast<unsigned int>(dim) &&
        detached.manifest.mesh_record.active_cells == 16 &&
        detached.manifest.mesh_record.provenance ==
          "test.unit_square.refine_2" &&
        detached.manifest.mesh_record.lifetime ==
          compiler::v1::MeshLifetimePolicy::owned_session &&
        detached.manifest.formulation_record.kind ==
          semantic::v1::FormulationKind::reduced_dto &&
        detached.manifest.state_solve_record.maximum_iterations == 317 &&
        detached.manifest.adjoint_solve_record.maximum_iterations == 419 &&
        detached.manifest.constraint_record.realisation_id == "l2_cellwise" &&
        detached.manifest.constraint_record.projection_metric_id ==
          detached.manifest.metric_record.realisation_id &&
        !detached.manifest.spaces.empty() &&
        !detached.manifest.bindings.empty(),
      "structured compilation manifest omitted resolved session decisions");
  }

  void
  run_serial_spd_reporting_contract_test()
  {
    dealii::DynamicSparsityPattern dynamic_pattern(2, 2);
    for (std::size_t row = 0; row < 2; ++row)
      for (std::size_t column = 0; column < 2; ++column)
        dynamic_pattern.add(row, column);
    dealii::SparsityPattern sparsity;
    sparsity.copy_from(dynamic_pattern);
    dealii::SparseMatrix<double> matrix(sparsity);
    matrix.set(0, 0, 4.0);
    matrix.set(0, 1, 1.0);
    matrix.set(1, 0, 1.0);
    matrix.set(1, 1, 3.0);
    dealii::Vector<double> right_hand_side(2);
    right_hand_side[0] = 1.0;
    right_hand_side[1] = 2.0;
    dealii::Vector<double> approximate_solution(2);
    const auto failed_report = dealii_backend::solve_serial_spd(
      matrix,
      approximate_solution,
      right_hand_side,
      dealii_backend::SPDLinearSolvePolicy{1, 1e-15, 1e-15});
    contract::require(!failed_report.converged() &&
                        failed_report.maximum_iterations == 1,
                      "serial SPD service did not report deliberate nonconvergence");

    dealii::Vector<double> zero_right_hand_side(2);
    dealii::Vector<double> zero_solution(2);
    const auto zero_report = dealii_backend::solve_serial_spd(
      matrix, zero_solution, zero_right_hand_side, {});
    contract::require(zero_report.converged() &&
                        zero_solution.l2_norm() == 0.0,
                      "serial SPD service did not handle a zero right-hand side");
  }

  template <int dim>
  void
  verify_homogeneous_weak_form_oracle()
  {
    static_assert(dim == 2,
                  "The hand-integrated weak-form oracle is two-dimensional");

    dealii::Triangulation<dim> triangulation;
    dealii::GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(1);

    const dealii::Functions::ConstantFunction<dim> forcing(1.0);
    const dealii::Functions::ConstantFunction<dim> desired_state(0.25);
    const dealii_backend::ScalarDiffusionReactionModel<dim> model(
      triangulation,
      forcing,
      desired_state,
      1.0,
      0.5,
      0.1,
      1);

    const std::size_t state_size = model.variable_layout()->dimension(0);
    std::size_t       interior_state_dof = state_size;
    std::size_t       unconstrained_dofs = 0;
    for (std::size_t index = 0; index < state_size; ++index)
      if (!model.state_constraints().is_constrained(index))
        {
          interior_state_dof = index;
          ++unconstrained_dofs;
        }
    contract::require(unconstrained_dofs == 1,
                      "the 2x2 Q1 oracle mesh must have one interior state DoF");

    dealii::Vector<double> state(state_size);
    state[interior_state_dof] = 1.0;
    dealii::Vector<double> control(model.variable_layout()->dimension(1));
    control = 2.0;
    const Primal point(model.variable_layout(),
                       {std::move(state), std::move(control)});

    // For the central Q1 hat on four squares of side 1/2:
    // integral |grad phi|^2 = 8/3, integral phi^2 = 1/9, and
    // integral phi = 1/4. With f = 1 and u = 2, the residual is
    // 8/3 + (1/2)(1/9) - 1/4 - 2(1/4) = 71/36.
    const Covector residual = model.residual(point);
    require_close(residual.block(0)[interior_state_dof],
                  71.0 / 36.0,
                  1e-13,
                  "independent Q1/DGQ0 weak-form residual oracle");

    // The tracking part is 7/288 for desired state 1/4; the constant
    // control contributes (0.1/2) integral 2^2 = 1/5.
    require_close(model.objective(point),
                  7.0 / 288.0 + 1.0 / 5.0,
                  1e-13,
                  "independent Q1/DGQ0 objective oracle");
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
    const dealii_backend::ScalarDiffusionReactionModel<dim> model(
      triangulation,
      forcing,
      desired_state,
      1.0,
      0.5,
      0.1,
      1);

    contract::StateControlPartitionT<Backend> partition(model, 0, 1);
    const contract::StateAdjointSolversT<Backend> solvers{
      [&model](const Primal &control) { return model.solve_state(control); },
      [&model](const Primal &full_point, const Covector &state_rhs) {
        return model.solve_adjoint(full_point, state_rhs);
      }};
    const contract::ReducedDTOT<Backend> reduced(model, partition, solvers);

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

    dealii_backend::MassMetricSolveParameters metric_parameters;
    metric_parameters.maximum_iterations = 1000;
    metric_parameters.relative_tolerance = 1e-12;
    metric_parameters.absolute_tolerance = 1e-14;
    const auto metric = model.control_l2_metric(metric_parameters);
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
    contract::require(solver_result.state_solve_count ==
                        solver_result.adjoint_solve_count,
                      "deal.II reduced gradient solve counts do not match");
    contract::require(solver_result.line_search_trial_count + 1 ==
                        solver_result.state_solve_count,
                      "deal.II reduced gradient solve count misses a trial evaluation");

    const auto bounds = model.control_l2_box_constraint(-1.0, 0.05, metric);
    dealii::Vector<double> bounded_control_values(
      partition.control_layout()->dimension(0));
    const Primal bounded_control(partition.control_layout(),
                                 {std::move(bounded_control_values)});
    const nmopt::solvers::ReducedGradientSolverT<Backend> projected_solver(
      reduced, metric, bounds, solver_parameters);
    const auto projected_result = projected_solver.solve(bounded_control);

    contract::require(
      projected_result.stopping_reason ==
        nmopt::solvers::ReducedGradientStoppingReason::gradient_tolerance,
      "deal.II projected reduced gradient solver did not reach stationarity");
    contract::require(bounds.is_feasible(projected_result.control),
                      "deal.II projected reduced gradient returned an infeasible control");
    contract::require(projected_result.gradient_norm_history.back() <=
                        solver_parameters.gradient_tolerance,
                      "deal.II projected reduced gradient final norm exceeds tolerance");
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

    const auto specification =
      semantic::v1::make_scalar_diffusion_reaction_problem(true);
    compiler::v1::DealiiDiscretisationPolicy compilation_policy;
    compilation_policy.state_degree = 1;
    const compiler::v1::DealiiCompiler v1_compiler;
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
    contract::require(compilation.succeeded(),
                      "v1 compiler failed to produce an executable problem");
    const auto &compiled_model = compilation.problem->executable_model();

    const Primal comparison_point =
      shifted(evaluation.full_point, tangent, 0.37);
    const Covector compiled_residual =
      compiled_model.residual(comparison_point);
    require_covector_close(compiled_residual,
                           model.residual(comparison_point),
                           1e-12,
                           "compiled/direct wiring residual differs");
    require_close(compiled_model.objective(evaluation.full_point),
                  model.objective(evaluation.full_point),
                  1e-12,
                  "compiled/direct wiring objective differs");
    require_covector_close(compiled_model.objective_derivative(
                             evaluation.full_point),
                           model.objective_derivative(evaluation.full_point),
                           1e-12,
                           "compiled/direct wiring objective derivative differs");

    const auto compiled_reduced = compilation.problem->make_reduced_dto();
    const auto compiled_evaluation = compiled_reduced.evaluate(control);
    require_close(compiled_evaluation.objective_value,
                  evaluation.objective_value,
                  1e-12,
                  "compiled/direct wiring reduced objective differs");
    require_covector_close(compiled_evaluation.reduced_derivative,
                           evaluation.reduced_derivative,
                           1e-12,
                           "compiled/direct wiring reduced derivative differs");
    const auto *compiled_constraint = compilation.problem->constraint();
    contract::require(compiled_constraint != nullptr &&
                        compiled_constraint->is_feasible(bounded_control),
                      "v1 compiler did not preserve the declared box constraint");
    const auto &manifest = compilation.problem->manifest();
    require_constraint_realisation(
      manifest,
      "FE_DGQ(0) coefficientwise l2_cellwise clipping",
      "canonical volume control");
    contract::require(manifest.semantic_problem_id == specification.id &&
                        manifest.provenance == "DTO" &&
                        manifest.execution == "assembled",
                      "v1 compiler did not record its compilation manifest");
  }
} // namespace

int
main(const int argc, char **argv)
{
  try
    {
      const std::vector<test_support::Scenario> scenarios{
        {"canonical_volume_control",
         "nmopt.dealii.canonical_volume_control",
         {"dealii", "compiler"},
         180,
         []() { run_canonical_volume_control_contract_test<2>(); }},
        {"fixed_dirichlet",
         "nmopt.dealii.fixed_dirichlet",
         {"dealii", "compiler"},
         60,
         []() { run_fixed_dirichlet_contract_test<2>(); }},
        {"dirichlet_control",
         "nmopt.dealii.dirichlet_control",
         {"dealii", "compiler"},
         60,
         []() { run_dirichlet_control_contract_test<2>(); }},
        {"subdomain_observation",
         "nmopt.dealii.subdomain_observation",
         {"dealii", "compiler"},
         60,
         []() { run_subdomain_observation_contract_test<2>(); }},
        {"neumann_boundary",
         "nmopt.dealii.neumann_boundary",
         {"dealii", "compiler"},
         60,
         []() { run_neumann_boundary_contract_test<2>(); }},
        {"h1_control",
         "nmopt.dealii.h1_control",
         {"dealii", "compiler"},
         60,
         []() { run_h1_control_regularisation_contract_test<2>(); }},
        {"coefficient_identification",
         "nmopt.dealii.coefficient_identification",
         {"dealii", "compiler"},
         60,
         []() { run_coefficient_identification_contract_test<2>(); }},
        {"pure_neumann",
         "nmopt.dealii.pure_neumann",
         {"dealii", "compiler"},
         60,
         []() { run_pure_neumann_contract_test<2>(); }},
        {"general_scalar_robin",
         "nmopt.dealii.general_scalar_robin",
         {"dealii", "compiler"},
         90,
         []() { run_general_scalar_robin_contract_test<2>(); }},
        {"projection_compatibility",
         "nmopt.dealii.projection_compatibility",
         {"dealii", "contract", "constraint"},
         60,
         run_projection_compatibility_contract_test},
        {"compiler_diagnostics",
         "nmopt.dealii.compiler_diagnostics",
         {"dealii", "compiler", "contract"},
         60,
         []() { run_compiler_diagnostics_contract_test<2>(); }},
        {"compiler_session",
         "nmopt.dealii.compiler_session",
         {"dealii", "compiler", "contract"},
         60,
         []() { run_compiler_session_contract_test<2>(); }},
        {"serial_spd_reporting",
         "nmopt.dealii.serial_spd_reporting",
         {"dealii", "contract", "solver"},
         30,
         run_serial_spd_reporting_contract_test},
        {"backend_size_conversion",
         "nmopt.dealii.backend_size_conversion",
         {"dealii", "contract", "adapter"},
         30,
         run_backend_size_conversion_contract_test}};
      const auto result = test_support::run_requested_scenarios(
        argc, argv, scenarios, std::cout);
      if (!result.listed)
        std::cout << "deal.II diffusion DTO contract scenario passed: "
                  << result.executed << '\n';
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "deal.II diffusion DTO contract test failed: "
                << exception.what() << '\n';
      return 1;
    }
}
