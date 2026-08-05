#include "nmopt/contract/reduced_dto.hpp"
#include "nmopt/compiler/v1/dealii_scalar_diffusion_reaction.hpp"
#include "nmopt/dealii/scalar_diffusion_reaction.hpp"
#include "nmopt/semantic/v1/problem_spec.hpp"
#include "nmopt/solvers/reduced_gradient.hpp"

#include <deal.II/base/function_lib.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <cmath>
#include <exception>
#include <functional>
#include <iostream>
#include <string>

namespace
{
  using namespace nmopt;
  using Backend = dealii_backend::SerialBackend;
  using Primal = contract::PrimalBlockT<Backend>;
  using Covector = contract::CovectorBlockT<Backend>;

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
      Backend::add_scaled(value.block(block), step, direction.block(block));
    return value;
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
      forcing, desired_state, 1.0, 0.5, 0.1};
    const auto missing_lifting = compiler.compile(specification,
                                                  triangulation,
                                                  missing_lifting_binding,
                                                  policy);
    contract::require(
      !missing_lifting.succeeded() &&
        missing_lifting.diagnostics.has_category(
          semantic::v1::DiagnosticCategory::lowerability),
      "v1 compiler did not diagnose missing fixed-Dirichlet data");

    auto bindings = compiler::v1::DealiiDataBindings<dim>{
      forcing, desired_state, 1.0, 0.5, 0.1};
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
      forcing, desired_state, 1.0, 0.5, 0.1};
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
    contract::require(
        manifest.lifting_realisation.find("y_phys = P_h y_hat + ell_0,h") !=
        std::string::npos &&
        manifest.data_rule.find("boundary DoFs") != std::string::npos &&
        manifest.transformation_ids.size() == 1,
      "v1 fixed-Dirichlet compilation manifest is incomplete");
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
      forcing, desired_state, 1.0, 0.5, 0.1};
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
    contract::require(
      one_manifest.observation_realisation ==
        "material-id volume restriction: 1" &&
        one_manifest.data_rule.find("analytic desired-state Function") !=
          std::string::npos,
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
      forcing, desired_state, 1.0, 0.5, 0.1};
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
    compiler::v1::DealiiDiscretisationPolicy policy;
    policy.state_degree = 1;
    const compiler::v1::DealiiCompiler compiler;
    contract::require(
      compiler.validate(specification, policy).valid(),
      "H1-control regularisation v1 graph did not validate for deal.II");

    auto discontinuous_control = specification;
    discontinuous_control.spaces.at(2).topology =
      semantic::v1::SpaceTopology::l2;
    const auto discontinuous_control_report =
      compiler.validate(discontinuous_control, policy);
    contract::require(
      discontinuous_control_report.has_category(
        semantic::v1::DiagnosticCategory::lowerability),
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
    contract::require(
      h1_control_box_report.has_category(
        semantic::v1::DiagnosticCategory::lowerability),
      "H1-control compiler did not reject the unsupported cellwise box");

    const compiler::v1::DealiiDataBindings<dim> bindings{
      forcing, desired_state, 1.0, 0.5, 0.2};
    const auto compilation = compiler.compile(specification,
                                              triangulation,
                                              bindings,
                                              policy);
    contract::require(compilation.succeeded(),
                      "H1-control regularisation v1 compilation failed");

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

    const auto evaluation = reduced.evaluate(control);
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
    const auto &manifest = compilation.problem->manifest();
    contract::require(
      manifest.control_space.find("continuous scalar FE_Q") != std::string::npos &&
        manifest.declared_assumptions.front().find("h1_control_regularisation") !=
          std::string::npos,
      "H1-control compilation manifest omitted the loss or control realization");
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
      zero_forcing, desired_state, 1.0, 0.5, 0.1};
    const auto rejected_reaction = compiler.compile(specification,
                                                    triangulation,
                                                    nonzero_reaction,
                                                    policy);
    contract::require(
      !rejected_reaction.succeeded() &&
        rejected_reaction.diagnostics.has_category(
          semantic::v1::DiagnosticCategory::lowerability),
      "pure-Neumann compiler did not reject a nonzero reaction");

    const compiler::v1::DealiiDataBindings<dim> incompatible_binding{
      incompatible_forcing, desired_state, 1.0, 0.0, 0.1};
    const auto rejected_forcing = compiler.compile(specification,
                                                   triangulation,
                                                   incompatible_binding,
                                                   policy);
    contract::require(
      !rejected_forcing.succeeded() &&
        rejected_forcing.diagnostics.has_category(
          semantic::v1::DiagnosticCategory::lowerability),
      "pure-Neumann compiler did not reject incompatible forcing");

    const compiler::v1::DealiiDataBindings<dim> compatible_binding{
      zero_forcing, desired_state, 1.0, 0.0, 0.1};
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
    bool rejected_control = false;
    try
      {
        (void)reduced.evaluate(incompatible_control);
      }
    catch (const contract::ContractError &)
      {
        rejected_control = true;
      }
    contract::require(rejected_control,
                      "pure-Neumann solve did not reject an incompatible boundary control");

    const auto &manifest = compilation.problem->manifest();
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
  run_contract_test()
  {
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

    const auto bounds = model.control_l2_box_constraint(-1.0, 0.05);
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
    contract::require(lowerability_diagnostic.has_category(
                        semantic::v1::DiagnosticCategory::lowerability),
                      "v1 compiler did not report an unsupported lowerability");

    auto unsupported_formulation = specification;
    unsupported_formulation.formulation.kind =
      semantic::v1::FormulationKind::all_at_once;
    const auto formulation_diagnostic =
      v1_compiler.validate(unsupported_formulation, compilation_policy);
    contract::require(
      formulation_diagnostic.has_category(
        semantic::v1::DiagnosticCategory::formulation_capability),
      "v1 compiler did not report an unsupported formulation capability");

    const compiler::v1::DealiiDataBindings<dim> data_bindings{
      forcing, desired_state, 1.0, 0.5, 0.1};
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
                           "v1 and v0 residual differ");
    require_close(compiled_model.objective(evaluation.full_point),
                  model.objective(evaluation.full_point),
                  1e-12,
                  "v1 and v0 objective differ");
    require_covector_close(compiled_model.objective_derivative(
                             evaluation.full_point),
                           model.objective_derivative(evaluation.full_point),
                           1e-12,
                           "v1 and v0 objective derivative differ");

    const auto compiled_reduced = compilation.problem->make_reduced_dto();
    const auto compiled_evaluation = compiled_reduced.evaluate(control);
    require_close(compiled_evaluation.objective_value,
                  evaluation.objective_value,
                  1e-12,
                  "v1 and v0 reduced objective differ");
    require_covector_close(compiled_evaluation.reduced_derivative,
                           evaluation.reduced_derivative,
                           1e-12,
                           "v1 and v0 reduced derivative differ");
    const auto *compiled_constraint = compilation.problem->constraint();
    contract::require(compiled_constraint != nullptr &&
                        compiled_constraint->is_feasible(bounded_control),
                      "v1 compiler did not preserve the declared box constraint");
    const auto &manifest = compilation.problem->manifest();
    contract::require(manifest.semantic_problem_id == specification.id &&
                        manifest.provenance == "DTO" &&
                        manifest.execution == "assembled" &&
                        manifest.constraint_realisation.find("FE_DGQ(0)") !=
                          std::string::npos,
                      "v1 compiler did not record its compilation manifest");
  }
} // namespace

int
main()
{
  try
    {
      run_contract_test<2>();
      run_fixed_dirichlet_contract_test<2>();
      run_subdomain_observation_contract_test<2>();
      run_neumann_boundary_contract_test<2>();
      run_h1_control_regularisation_contract_test<2>();
      run_pure_neumann_contract_test<2>();
      std::cout << "deal.II diffusion DTO contract test passed\n";
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "deal.II diffusion DTO contract test failed: "
                << exception.what() << '\n';
      return 1;
    }
}
