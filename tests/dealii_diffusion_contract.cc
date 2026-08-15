#include "nmopt/contract/reduced_dto.hpp"
#include "nmopt/compiler/v1/dealii_scalar_diffusion_reaction.hpp"
#include "nmopt/dealii/hminus1_metric.hpp"
#include "nmopt/dealii/scalar_diffusion_reaction.hpp"
#include "nmopt/semantic/v1/problem_spec.hpp"
#include "nmopt/solvers/reduced_gradient.hpp"
#include "test_support/contract_errors.hpp"
#include "test_support/diagnostics.hpp"
#include "test_support/manifest_contracts.hpp"
#include "test_support/scenario_dispatch.hpp"

#include <deal.II/base/function_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/numerics/vector_tools.h>

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

  void
  require_close(double              actual,
                double              expected,
                double              tolerance,
                const std::string &description);

  void
  run_hminus1_metric_contract_test()
  {
    const auto layout = std::make_shared<const contract::BlockLayout>(
      "hminus1_metric_test",
      std::vector<contract::SpaceId>{{"control"}},
      std::vector<std::size_t>{2});
    dealii::DynamicSparsityPattern dsp(2, 2);
    dsp.add(0, 0);
    dsp.add(0, 1);
    dsp.add(1, 0);
    dsp.add(1, 1);
    dealii::SparsityPattern sparsity;
    sparsity.copy_from(dsp);

    auto mass = std::make_shared<dealii::SparseMatrix<double>>();
    mass->reinit(sparsity);
    mass->set(0, 0, 2.0);
    mass->set(1, 1, 3.0);
    auto laplace = std::make_shared<dealii::SparseMatrix<double>>();
    laplace->reinit(sparsity);
    laplace->set(0, 0, 4.0);
    laplace->set(0, 1, 1.0);
    laplace->set(1, 0, 1.0);
    laplace->set(1, 1, 2.0);

    dealii_backend::MetricSolveParameters solve_parameters;
    solve_parameters.maximum_iterations = 100;
    solve_parameters.relative_tolerance = 1e-13;
    solve_parameters.absolute_tolerance = 1e-15;
    const dealii_backend::Hminus1Metric metric(
      "hminus1_continuous", layout, mass, laplace, solve_parameters);

    dealii::Vector<double> primal_values(2);
    primal_values[0] = 1.0;
    primal_values[1] = -2.0;
    const Primal primal(layout, {std::move(primal_values)});
    const Covector applied = metric.apply(primal);
    require_close(applied.block(0)[0],
                  20.0 / 7.0,
                  1e-12,
                  "H-1 metric M K^-1 M first component");
    require_close(applied.block(0)[1],
                  -78.0 / 7.0,
                  1e-12,
                  "H-1 metric M K^-1 M second component");

    const Primal recovered = metric.inverse_apply(applied);
    dealii::Vector<double> recovery_error = recovered.block(0);
    recovery_error.add(-1.0, primal.block(0));
    require_close(recovery_error.l2_norm(),
                  0.0,
                  1e-11,
                  "H-1 metric apply/inverse pairing");

    dealii::Vector<double> second_values(2);
    second_values[0] = 0.25;
    second_values[1] = 0.5;
    const Primal second(layout, {std::move(second_values)});
    require_close(contract::pair(metric.apply(primal), second),
                  contract::pair(metric.apply(second), primal),
                  1e-12,
                  "H-1 metric symmetry");
    contract::require(metric.id() == "hminus1_continuous" &&
                        metric.solve_parameters().maximum_iterations == 100,
                      "H-1 metric omitted its identity or solve policy");
  }

  template <typename Component>
  Component &
  component_by_id(std::vector<Component> &components, const std::string &id)
  {
    const auto component = std::find_if(
      components.begin(), components.end(), [&id](const Component &candidate) {
        return candidate.id == id;
      });
    contract::require(component != components.end(),
                      "deal.II test semantic component is missing");
    return *component;
  }

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

  template <int dim>
  class FirstCoordinateFunction final : public dealii::Function<dim>
  {
  public:
    double
    value(const dealii::Point<dim> &point,
          const unsigned int        component = 0) const override
    {
      (void)component;
      return point[0];
    }
  };

  template <int dim>
  class EnergyPolynomial final : public dealii::Function<dim>
  {
  public:
    explicit EnergyPolynomial(const double scale = 1.0)
      : scale_(scale)
    {}

    double
    value(const dealii::Point<dim> &point,
          const unsigned int        component = 0) const override
    {
      (void)component;
      double result = scale_;
      for (unsigned int direction = 0; direction < dim; ++direction)
        result *= point[direction] * (1.0 - point[direction]);
      return result;
    }

    dealii::Tensor<1, dim>
    gradient(const dealii::Point<dim> &point,
             const unsigned int        component = 0) const override
    {
      (void)component;
      dealii::Tensor<1, dim> result;
      for (unsigned int derivative = 0; derivative < dim; ++derivative)
        {
          result[derivative] = scale_ * (1.0 - 2.0 * point[derivative]);
          for (unsigned int direction = 0; direction < dim; ++direction)
            if (direction != derivative)
              result[derivative] *=
                point[direction] * (1.0 - point[direction]);
        }
      return result;
    }

    double
    laplacian(const dealii::Point<dim> &point) const
    {
      double result = 0.0;
      for (unsigned int derivative = 0; derivative < dim; ++derivative)
        {
          double contribution = -2.0 * scale_;
          for (unsigned int direction = 0; direction < dim; ++direction)
            if (direction != derivative)
              contribution *= point[direction] * (1.0 - point[direction]);
          result += contribution;
        }
      return result;
    }

  private:
    double scale_;
  };

  template <int dim>
  class EnergyPolynomialForcing final : public dealii::Function<dim>
  {
  public:
    explicit EnergyPolynomialForcing(const double reaction,
                                     const double diffusion = 1.0)
      : diffusion_(diffusion)
      , reaction_(reaction)
    {}

    double
    value(const dealii::Point<dim> &point,
          const unsigned int        component = 0) const override
    {
      return -diffusion_ * state_.laplacian(point) +
             reaction_ * state_.value(point, component);
    }

  private:
    double                diffusion_;
    double                reaction_;
    EnergyPolynomial<dim> state_;
  };

  template <int dim>
  class RightBoundaryNormalFluxFunction final : public dealii::Function<dim>
  {
  public:
    explicit RightBoundaryNormalFluxFunction(const double scale = 1.0)
      : scale_(scale)
    {}

    double
    value(const dealii::Point<dim> &point,
          const unsigned int        component = 0) const override
    {
      (void)component;
      return -scale_ * point[1] * (1.0 - point[1]);
    }

  private:
    double scale_;
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

  void
  require_constraint_realisation(
    const compiler::v1::CompilationManifest &manifest,
    const std::string &                       expected,
    const std::string &                       target)
  {
    const auto &compatibility = manifest.resolved_decision.compatibility;
    contract::require(
      manifest.compiler_id == compatibility.compiler_id &&
        manifest.backend == compatibility.backend &&
        manifest.execution == compatibility.execution &&
        manifest.state_space == compatibility.state_space &&
        manifest.control_space == compatibility.control_space &&
        manifest.quadrature == compatibility.quadrature &&
        manifest.dual_representation == compatibility.dual_representation &&
        manifest.data_rule == compatibility.data_rule &&
        manifest.observation_realisation ==
          compatibility.observation_realisation &&
        manifest.metric_solve_policy == compatibility.metric_solve_policy &&
        manifest.constraint_realisation ==
          compatibility.constraint_realisation &&
        manifest.lifting_realisation == compatibility.lifting_realisation &&
        manifest.nullspace_policy == compatibility.nullspace_policy &&
        manifest.state_adjoint_solve_policy ==
          compatibility.state_adjoint_solve_policy &&
        manifest.provenance == compatibility.provenance &&
        manifest.lowering_handler_records ==
          compatibility.lowering_handler_records &&
        manifest.region_ids == compatibility.region_ids &&
        manifest.space_ids == compatibility.space_ids &&
        manifest.pairing_ids == compatibility.pairing_ids &&
        manifest.variable_ids == compatibility.variable_ids &&
        manifest.data_ids == compatibility.data_ids &&
        manifest.transformation_ids == compatibility.transformation_ids &&
        manifest.residual_term_ids == compatibility.residual_term_ids &&
        manifest.observation_ids == compatibility.observation_ids &&
        manifest.loss_ids == compatibility.loss_ids &&
        manifest.metric_ids == compatibility.metric_ids &&
        manifest.constraint_ids == compatibility.constraint_ids &&
        manifest.declared_assumptions == compatibility.declared_assumptions,
      target + " manifest compatibility view was not projected from the decision");
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
    const auto has_realized_map = [&manifest](const std::string &id) {
      return std::any_of(
        manifest.realized_maps.begin(),
        manifest.realized_maps.end(),
        [&id](const compiler::v1::CompiledRealizedMapRecord &map) {
          return map.semantic_id == id && !map.source_space_id.empty() &&
                 !map.output_space_id.empty() && map.source_dimension > 0 &&
                 map.output_dimension > 0 && !map.realization_id.empty() &&
                 !map.value_provenance.empty() && !map.jvp_provenance.empty() &&
                 !map.vjp_provenance.empty();
        });
    };
    const auto control_restrictions_are_coefficient_maps =
      std::all_of(
        manifest.realized_maps.begin(),
        manifest.realized_maps.end(),
        [](const compiler::v1::CompiledRealizedMapRecord &map) {
          return map.realization_id != "coefficient_restriction" ||
                 (map.input_dimensions.size() == 1 &&
                  map.input_dimensions.front() == map.output_dimension);
        });
    const auto observation_map_dimensions_match_spaces =
      std::all_of(
        manifest.realized_maps.begin(),
        manifest.realized_maps.end(),
        [&manifest](const compiler::v1::CompiledRealizedMapRecord &map) {
          if (std::find(manifest.observation_ids.begin(),
                        manifest.observation_ids.end(),
                        map.semantic_id) == manifest.observation_ids.end())
            return true;
          const auto space = std::find_if(
            manifest.spaces.begin(),
            manifest.spaces.end(),
            [&map](const compiler::v1::CompiledSpaceRecord &candidate) {
              return candidate.semantic_id == map.output_space_id;
            });
          return space != manifest.spaces.end() &&
                 space->dimension == map.output_dimension;
        });
    contract::require(
      std::all_of(manifest.observation_ids.begin(),
                  manifest.observation_ids.end(),
                  has_realized_map) &&
        std::all_of(manifest.transformation_ids.begin(),
                    manifest.transformation_ids.end(),
                    has_realized_map) &&
        control_restrictions_are_coefficient_maps &&
        observation_map_dimensions_match_spaces,
      target + " manifest omitted a realized observation or transformation map");
    contract::require(
      manifest.schema_version == 3 &&
        manifest.formulation_record.kind ==
          semantic::v1::FormulationKind::reduced_dto &&
        manifest.formulation_record.provenance ==
          semantic::v1::FormulationProvenance::dto &&
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
        !manifest.metric_record.realisation_id.empty() &&
        !manifest.mesh_record.structural_identity.empty() &&
        manifest.resolved_decision.semantic_problem_id ==
          manifest.semantic_problem_id &&
        manifest.resolved_decision.formulation_id ==
          manifest.formulation_record.semantic_id &&
        manifest.resolved_decision.spaces.size() == manifest.spaces.size() &&
        manifest.resolved_decision.bindings.size() == manifest.bindings.size() &&
        manifest.resolved_decision.realized_maps.size() ==
          manifest.realized_maps.size() &&
        manifest.resolved_decision.realized_spaces.size() ==
          manifest.realized_spaces.size() &&
        !manifest.resolved_decision.residuals.empty() &&
        !manifest.resolved_decision.observations.empty(),
      target + " structured manifest is incomplete");
  }

  void
  require_compiled_binding_records_equal(
    const std::vector<compiler::v1::CompiledBindingRecord> &actual,
    const std::vector<compiler::v1::CompiledBindingRecord> &expected,
    const std::string &                                      description)
  {
    contract::require(actual.size() == expected.size(),
                      description + " changed its binding record count");
    for (const auto &record : actual)
      {
        const auto match = std::find_if(
          expected.begin(),
          expected.end(),
          [&record](const compiler::v1::CompiledBindingRecord &candidate) {
            return candidate.semantic_id == record.semantic_id;
          });
        contract::require(match != expected.end(),
                          description + " lost binding " + record.semantic_id);
        contract::require(
          record.role == match->role && record.kind == match->kind &&
            record.space_id == match->space_id &&
            record.region_id == match->region_id &&
            record.representation == match->representation &&
            record.evaluation_realisation == match->evaluation_realisation &&
            record.runtime_representation == match->runtime_representation &&
            record.provenance == match->provenance &&
            record.field_shape == match->field_shape &&
            record.scalar_value == match->scalar_value &&
            record.value_digest == match->value_digest &&
            record.value_status == match->value_status,
          description + " changed binding " + record.semantic_id);
      }
  }

  void
  require_resolved_manifest_projection(
    const compiler::v1::CompilationManifest &manifest,
    const std::string &                       description)
  {
    const auto &decision = manifest.resolved_decision;
    contract::require(
      decision.formulation_record.semantic_id ==
          manifest.formulation_record.semantic_id &&
        decision.formulation_record.kind == manifest.formulation_record.kind &&
        decision.formulation_record.provenance ==
          manifest.formulation_record.provenance &&
        decision.formulation_record.execution ==
          manifest.formulation_record.execution &&
        decision.formulation_record.dual_representation ==
          manifest.formulation_record.dual_representation &&
        decision.mesh_record.dimension == manifest.mesh_record.dimension &&
        decision.mesh_record.active_cells == manifest.mesh_record.active_cells &&
        decision.mesh_record.provenance == manifest.mesh_record.provenance &&
        decision.mesh_record.lifetime == manifest.mesh_record.lifetime &&
        decision.mesh_record.structural_identity ==
          manifest.mesh_record.structural_identity &&
        decision.spaces.size() == manifest.spaces.size() &&
        decision.bindings.size() == manifest.bindings.size() &&
        decision.state_solve_record.algorithm ==
          manifest.state_solve_record.algorithm &&
        decision.state_solve_record.maximum_iterations ==
          manifest.state_solve_record.maximum_iterations &&
        decision.state_solve_record.relative_tolerance ==
          manifest.state_solve_record.relative_tolerance &&
        decision.state_solve_record.absolute_tolerance ==
          manifest.state_solve_record.absolute_tolerance &&
        decision.state_solve_record.nullspace_policy ==
          manifest.state_solve_record.nullspace_policy &&
        decision.state_solve_record.operator_realisation ==
          manifest.state_solve_record.operator_realisation &&
        decision.adjoint_solve_record.algorithm ==
          manifest.adjoint_solve_record.algorithm &&
        decision.adjoint_solve_record.maximum_iterations ==
          manifest.adjoint_solve_record.maximum_iterations &&
        decision.adjoint_solve_record.relative_tolerance ==
          manifest.adjoint_solve_record.relative_tolerance &&
        decision.adjoint_solve_record.absolute_tolerance ==
          manifest.adjoint_solve_record.absolute_tolerance &&
        decision.adjoint_solve_record.nullspace_policy ==
          manifest.adjoint_solve_record.nullspace_policy &&
        decision.adjoint_solve_record.operator_realisation ==
          manifest.adjoint_solve_record.operator_realisation &&
        manifest.dual_representation ==
          decision.formulation_record.dual_representation &&
        manifest.execution == decision.execution_id &&
        manifest.metric_solve_policy.find(
          manifest.metric_record.realisation_id) != std::string::npos &&
        manifest.state_adjoint_solve_policy.find(
          manifest.state_solve_record.operator_realisation) !=
          std::string::npos,
      description + " was not rendered from its resolved typed records");
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
      manifest.realized_maps.begin(),
      manifest.realized_maps.end(),
      [](const compiler::v1::CompiledRealizedMapRecord &map) {
        return map.semantic_id == "fixed_dirichlet_reconstruction";
      });
    require_constraint_realisation(manifest, "none", "fixed-Dirichlet");
    contract::require(fixed_map != manifest.realized_maps.end(),
                      "fixed-Dirichlet realized transformation is missing");
    contract::require(
      fixed_map->input_dimensions.size() == 1 &&
        fixed_map->input_dimensions.front() != fixed_map->output_dimension,
      "fixed-Dirichlet realized transformation dimensions are not distinct");
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
    const auto dirichlet_map = std::find_if(
      manifest.realized_maps.begin(),
      manifest.realized_maps.end(),
      [](const compiler::v1::CompiledRealizedMapRecord &map) {
        return map.semantic_id == "dirichlet_control_lifting";
      });
    require_constraint_realisation(manifest, "none", "Dirichlet-control");
    contract::require(
      manifest.control_space.find("nodal trace") != std::string::npos &&
        manifest.lifting_realisation.find("L_D,h") != std::string::npos &&
        dirichlet_map != manifest.realized_maps.end() &&
        dirichlet_map->input_space_ids.size() == 2 &&
        dirichlet_map->output_dimension > 0 &&
        manifest.metric_solve_policy.find("l2_dirichlet_trace") !=
          std::string::npos &&
        manifest.declared_assumptions.front().find("dirichlet_control_lifting") !=
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
        manifest.declared_assumptions.begin(),
        manifest.declared_assumptions.end(),
        [&prefix](const std::string &assumption) {
          return assumption.find(prefix) == 0;
        });
    };
    const auto has_binding_role = [&manifest](const auto role) {
      return std::any_of(
        manifest.bindings.begin(),
        manifest.bindings.end(),
        [role](const compiler::v1::CompiledBindingRecord &binding) {
          return binding.role == role;
        });
    };
    contract::require(
      manifest.compiler_id ==
          "nmopt.compiler.v1.dealii.l2_dirichlet_transposition" &&
        manifest.transposition_realisation.has_value() &&
        manifest.transposition_realisation->id ==
          "transposition_formulation" &&
        manifest.transposition_realisation->continuous_parent_space_id ==
          "control_space" &&
        manifest.transposition_realisation->equivalence_policy_id ==
          "conforming_trace_subspace" &&
        manifest.transposition_realisation->discrete_realisation ==
          semantic::v1::TranspositionDiscreteRealisation::
            conforming_nodal_lifting_equivalence &&
        manifest.state_space.find("continuous L2(Omega) parent") !=
          std::string::npos &&
        manifest.control_space.find("U_h=trace(V_h)") != std::string::npos &&
        manifest.lifting_realisation.find("E_tr(y,u;f)") !=
          std::string::npos &&
        manifest.transformation_ids.empty() &&
        std::find(manifest.lowering_handler_records.begin(),
                  manifest.lowering_handler_records.end(),
                  "l2_dirichlet_transposition <- "
                  "dealii.dirichlet_control.conforming_trace_equivalence") !=
          manifest.lowering_handler_records.end() &&
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
      manifest.partial_boundary_selection.has_value() &&
        manifest.partial_boundary_selection->fixed_boundary_region_id ==
          "fixed_dirichlet_boundary" &&
        manifest.partial_boundary_selection->controlled_boundary_region_id ==
          "control_boundary" &&
        manifest.partial_boundary_selection->interface_realisation ==
          semantic::v1::PartialDirichletInterfaceRealisation::
            fixed_data_precedence &&
        manifest.partial_boundary_selection->trace_realisation ==
          semantic::v1::PartialDirichletTraceRealisation::
            relative_interior_nodal_zero_endpoint &&
      manifest.lifting_realisation.find("ell_0,h + L_D,h") != std::string::npos &&
        manifest.data_rule.find("fixed Dirichlet Function") != std::string::npos &&
        std::any_of(manifest.declared_assumptions.begin(),
                    manifest.declared_assumptions.end(),
                    [](const std::string &assumption) {
                      return assumption.find("fixed-data precedence") !=
                             std::string::npos;
                    }),
      "partial Dirichlet-control manifest omitted the lifting interface policy");
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
      fixed_full.problem->manifest().lowering_handler_records ==
        fixed_subdomain.problem->manifest().lowering_handler_records &&
        fixed_full.problem->manifest().metric_record.realisation_id ==
          fixed_subdomain.problem->manifest().metric_record.realisation_id &&
        fixed_subdomain.problem->manifest().observation_realisation.find(
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
      manifest.h1_target_data_membership_selection.has_value() &&
        manifest.h1_target_data_membership_selection->data_id ==
          "desired_state" &&
        manifest.h1_target_data_membership_selection->observation_space_id ==
          "state_observation_space" &&
        manifest.h1_target_data_membership_selection
            ->fixed_boundary_region_id == "dirichlet_boundary" &&
        std::any_of(manifest.declared_assumptions.begin(),
                    manifest.declared_assumptions.end(),
                    [](const std::string &assumption) {
                      return assumption.find(
                               "h1_target_data_membership: status=user_assumed") ==
                             0;
                    }) &&
      manifest.compiler_id == "nmopt.compiler.v1.dealii.h1_state_tracking" &&
        manifest.observation_realisation.find("H1_0") != std::string::npos &&
        manifest.observation_realisation.find("mass-plus-stiffness") !=
          std::string::npos &&
        manifest.data_rule.find("value and gradient") != std::string::npos &&
        manifest.metric_record.realisation_id == "l2_cellwise" &&
        std::find(manifest.lowering_handler_records.begin(),
                  manifest.lowering_handler_records.end(),
                  "state_observation <- "
                  "dealii.scalar.observation.h1_state_restriction") !=
          manifest.lowering_handler_records.end(),
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
      manifest.realized_maps.begin(),
      manifest.realized_maps.end(),
      [](const compiler::v1::CompiledRealizedMapRecord &map) {
        return map.semantic_id == "state_observation" &&
               map.realization_id == "ordered_point_sensor_values";
      });
    const auto point_space = std::find_if(
      manifest.spaces.begin(),
      manifest.spaces.end(),
      [](const compiler::v1::CompiledSpaceRecord &space) {
        return space.semantic_id == "state_observation_space" &&
               space.role == semantic::v1::SpaceRole::observation;
      });
    contract::require(
      point_space != manifest.spaces.end() &&
        point_space->dimension == values.size() &&
        point_space->dimension == sensor_jvp.block(0).size(),
      "point-sensor manifest recorded a dimension different from its realized output");
    contract::require(
      manifest.compiler_id == "nmopt.compiler.v1.dealii.point_sensor" &&
        point_map != manifest.realized_maps.end() &&
        point_map->output_dimension == values.size() &&
        point_map->output_layout.find("sensor values") != std::string::npos &&
        manifest.observation_realisation.find("immutable physical coordinates") !=
          std::string::npos &&
        manifest.transposition_realisation.has_value() &&
        manifest.transposition_realisation->diffusion_data_id == "diffusion" &&
        manifest.transposition_realisation->reaction_data_id == "reaction" &&
        manifest.data_rule.find("assembled C_h^T point-load transpose") !=
          std::string::npos &&
        std::any_of(manifest.declared_assumptions.begin(),
                    manifest.declared_assumptions.end(),
                    [](const std::string &assumption) {
                      return assumption.find("very-weak adjoint source") !=
                             std::string::npos;
                    }) &&
        std::find(manifest.lowering_handler_records.begin(),
                  manifest.lowering_handler_records.end(),
                  "state_observation <- dealii.scalar.observation.point_sensor") !=
          manifest.lowering_handler_records.end(),
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
        nonunit_compilation.problem->manifest().transposition_realisation
          .has_value() &&
        nonunit_compilation.problem->manifest().transposition_realisation
            ->diffusion_data_id == "diffusion" &&
        nonunit_compilation.problem->manifest().transposition_realisation
            ->reaction_data_id == "reaction" &&
        nonunit_compilation.problem->manifest().data_rule.find(
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
      manifest.realized_maps.begin(),
      manifest.realized_maps.end(),
      [](const compiler::v1::CompiledRealizedMapRecord &map) {
        return map.semantic_id == "state_observation" &&
               map.realization_id == "ordered_normal_flux_face_quadrature";
      });
    const auto normal_flux_space = std::find_if(
      manifest.spaces.begin(),
      manifest.spaces.end(),
      [](const compiler::v1::CompiledSpaceRecord &space) {
        return space.semantic_id == "state_observation_space" &&
               space.role == semantic::v1::SpaceRole::observation;
      });
    contract::require(
      normal_flux_space != manifest.spaces.end() &&
        normal_flux_space->dimension == normal_flux_values.size() &&
        normal_flux_space->dimension == normal_flux_jvp.block(0).size(),
      "normal-flux manifest recorded the state dimension instead of its realized face output");
    contract::require(
      manifest.compiler_id == "nmopt.compiler.v1.dealii.normal_flux" &&
        normal_flux_map != manifest.realized_maps.end() &&
        normal_flux_map->output_dimension == normal_flux_values.size() &&
        normal_flux_map->pairing_realization.find("declared pairing") !=
          std::string::npos &&
        manifest.observation_realisation.find("outward normal-flux") !=
          std::string::npos &&
        manifest.transposition_realisation.has_value() &&
        manifest.transposition_realisation->diffusion_data_id == "diffusion" &&
        manifest.transposition_realisation->reaction_data_id == "reaction" &&
        manifest.data_rule.find(
          "T=-kappa Delta+rI with kappa <- diffusion and r <- reaction") !=
          std::string::npos &&
        manifest.data_rule.find("boundary face quadrature") !=
          std::string::npos &&
        std::any_of(manifest.declared_assumptions.begin(),
                    manifest.declared_assumptions.end(),
                    [](const std::string &assumption) {
                      return assumption.find("very-weak adjoint boundary source") !=
                             std::string::npos;
                    }) &&
        std::find(manifest.lowering_handler_records.begin(),
                  manifest.lowering_handler_records.end(),
                  "state_observation <- dealii.scalar.observation.normal_flux") !=
          manifest.lowering_handler_records.end(),
      "normal-flux compilation manifest omitted its face transpose policy");
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
    const auto *neumann_model =
      dynamic_cast<const compiler::v1::detail::NeumannBoundaryControlModel<dim> *>(
        &model);
    contract::require(neumann_model != nullptr,
                      "Neumann boundary-control did not produce its Neumann target");
    const dealii::QGauss<dim - 1> boundary_quadrature(policy.state_degree + 2);
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
        manifest.constraint_realisation.find("l2_facewise") != std::string::npos &&
        std::any_of(
          manifest.realized_maps.begin(),
          manifest.realized_maps.end(),
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
      manifest.compiler_id ==
          "nmopt.compiler.v1.dealii.neumann_convection_subdomain" &&
        manifest.observation_realisation == "material-id volume restriction: 1" &&
        manifest.state_solve_record.algorithm ==
          compiler::v1::LinearSolveAlgorithm::serial_sparse_direct_umfpack &&
        std::any_of(manifest.bindings.begin(), manifest.bindings.end(),
                    [](const compiler::v1::CompiledBindingRecord &record) {
                      return record.semantic_id == "conservative_transport" &&
                             record.provenance ==
                               "test.neumann_convection.conservative_transport";
                    }),
      "C5.6 manifest omitted transport, subdomain, or solve provenance");
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
      manifest.realized_maps.begin(),
      manifest.realized_maps.end(),
      [](const compiler::v1::CompiledRealizedMapRecord &map) {
        return map.semantic_id == "weighted_state_boundary_trace" &&
               map.realization_id == "ordered_boundary_face_quadrature_trace";
      });
    const auto weight_record = std::find_if(
      manifest.bindings.begin(),
      manifest.bindings.end(),
      [](const compiler::v1::CompiledBindingRecord &record) {
        return record.semantic_id == "boundary_weight";
      });
    contract::require(
      weight_record != manifest.bindings.end() &&
        weight_record->provenance == "test.weighted_boundary.weight" &&
        weight_record->representation.find("boundary face quadrature") !=
          std::string::npos &&
        manifest.compiler_id ==
          "nmopt.compiler.v1.dealii.weighted_boundary_trace" &&
        manifest.observation_realisation.find("weighted boundary trace") !=
          std::string::npos &&
        manifest.data_rule.find("desired-state and fixed boundary-weight") !=
          std::string::npos &&
        manifest.data_rule.find("boundary face quadrature") !=
          std::string::npos &&
        manifest.metric_record.realisation_id == "l2_facewise" &&
        std::find(manifest.lowering_handler_records.begin(),
                  manifest.lowering_handler_records.end(),
                  "weighted_state_boundary_trace <- "
                  "dealii.neumann.observation.weighted_boundary_trace") !=
          manifest.lowering_handler_records.end(),
      "weighted trace manifest omitted weight, target, quadrature, or metric provenance");
    contract::require(
      weighted_map != manifest.realized_maps.end() &&
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
      manifest.h1_target_data_membership_selection.has_value() &&
        l2_manifest.h1_target_data_membership_selection.has_value() &&
        manifest.h1_target_data_membership_selection->data_id ==
          l2_manifest.h1_target_data_membership_selection->data_id &&
        manifest.h1_target_data_membership_selection
            ->fixed_boundary_region_id == "dirichlet_boundary" &&
        std::any_of(manifest.declared_assumptions.begin(),
                    manifest.declared_assumptions.end(),
                    [](const std::string &assumption) {
                      return assumption.find(
                               "h1_target_data_membership: status=user_assumed") ==
                             0;
                    }) &&
      hminus1_compilation.problem->metric().id() == "hminus1_continuous" &&
        l2_compilation.problem->metric().id() == "l2_continuous" &&
        manifest.metric_record.operator_description.find("M_h K_h^{-1} M_h") !=
          std::string::npos &&
        manifest.metric_record.operator_id ==
          "mass_laplacian_inverse_mass" &&
        manifest.metric_record.inverse_operator_id ==
          "mass_inverse_laplacian_mass_inverse" &&
        manifest.metric_record.boundary_region_id == "dirichlet_boundary" &&
        manifest.metric_record.laplacian_solve_policy_id ==
          "control_metric_solve.laplacian_inverse" &&
        manifest.metric_record.mass_solve_policy_id ==
          "control_metric_solve.mass_inverse" &&
        manifest.metric_solve_policy.find("hminus1_continuous") !=
          std::string::npos &&
        manifest.control_space.find("independent homogeneous-Dirichlet") !=
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
    const auto &boundary_selection = manifest.boundary_realisation;
    const auto has_binding = [&manifest](
                               const std::string &semantic_id,
                               const semantic::v1::DataRole role,
                               const semantic::v1::DataKind kind,
                               const std::string &space_id,
                               const std::string &region_id,
                               const std::string &evaluation) {
      return std::any_of(
        manifest.bindings.begin(),
        manifest.bindings.end(),
        [&](const compiler::v1::CompiledBindingRecord &binding) {
          return binding.semantic_id == semantic_id && binding.role == role &&
                 binding.kind == kind && binding.space_id == space_id &&
                 binding.region_id == region_id &&
                 binding.evaluation_realisation == evaluation &&
                 !binding.provenance.empty();
        });
    };
    contract::require(
        manifest.state_solve_record.algorithm ==
          compiler::v1::LinearSolveAlgorithm::serial_sparse_direct_umfpack &&
        manifest.adjoint_solve_record.algorithm ==
          compiler::v1::LinearSolveAlgorithm::serial_sparse_direct_umfpack &&
        manifest.lowering_handler_records.size() == 13 &&
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
        manifest.data_rule.find("Robin coefficient and source") !=
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
        std::any_of(manifest.declared_assumptions.begin(),
                    manifest.declared_assumptions.end(),
                    [](const std::string &assumption) {
                      return assumption.find(
                               "general_scalar_robin: boundary selection scalar_boundary_partition") ==
                             0;
                    }) &&
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
      detached.manifest.schema_version == 3 &&
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
        !detached.manifest.mesh_record.structural_identity.empty() &&
        detached.manifest.resolved_decision.semantic_problem_id ==
          detached.manifest.semantic_problem_id &&
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

    const contract::ReducedHessianT<Backend> &hessian = model;
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
    contract::require(solver_result.state_solve_count ==
                        solver_result.adjoint_solve_count,
                      "deal.II reduced gradient solve counts do not match");
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

    auto supplied_otd_specification =
      semantic::v1::make_scalar_diffusion_reaction_problem(false);
    supplied_otd_specification.formulation.kind =
      semantic::v1::FormulationKind::all_at_once;
    supplied_otd_specification.formulation.provenance =
      semantic::v1::FormulationProvenance::supplied_otd;
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
    mismatched_supplied_otd.equations.front().test_space_id =
      "alternate_state_test_space";
    mismatched_supplied_otd.equations.front().test_pairing_id =
      "alternate_state_test_pairing";
    const auto mismatched_otd_diagnostic =
      v1_compiler.validate(mismatched_supplied_otd, compilation_policy);
    test_support::require_exact_diagnostic(
      mismatched_otd_diagnostic,
      semantic::v1::DiagnosticCategory::formulation_capability,
      "state_equation",
      "supplied_otd_adjoint_space",
      "v1 compiler reported a mismatched supplied OTD adjoint space as valid");

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
    const auto &supplied_manifest = supplied_otd.manifest();
    contract::require(
      supplied_manifest.formulation_record.kind ==
          semantic::v1::FormulationKind::all_at_once &&
        supplied_manifest.formulation_record.provenance ==
          semantic::v1::FormulationProvenance::supplied_otd &&
        supplied_manifest.provenance == "supplied OTD" &&
        supplied_manifest.supplied_otd_record.present &&
        supplied_manifest.supplied_otd_record.variable_space_ids ==
          std::vector<std::string>{"state", "state_test", "control"} &&
        supplied_manifest.supplied_otd_record.residual_space_ids ==
          std::vector<std::string>{"state_equation",
                                   "adjoint_equation",
                                   "control_stationarity"} &&
        supplied_manifest.supplied_otd_record.comparison_status.find(
          "equivalence verified") != std::string::npos,
      "serial supplied OTD manifest omitted formulation and comparison provenance");
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
        coarser_compilation.problem->manifest().mesh_record.structural_identity !=
          manifest.mesh_record.structural_identity,
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
        relabeled_compilation.problem->manifest().metric_record.realisation_id ==
          manifest.metric_record.realisation_id,
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
        reordered_manifest.formulation_record.semantic_id ==
          manifest.formulation_record.semantic_id &&
        reordered_manifest.mesh_record.structural_identity ==
          manifest.mesh_record.structural_identity &&
        reordered_manifest.metric_record.realisation_id ==
          manifest.metric_record.realisation_id &&
        reordered_manifest.constraint_record.realisation_id ==
          manifest.constraint_record.realisation_id &&
        reordered_manifest.state_solve_record.maximum_iterations ==
          manifest.state_solve_record.maximum_iterations &&
        reordered_manifest.adjoint_solve_record.maximum_iterations ==
          manifest.adjoint_solve_record.maximum_iterations,
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
        record.bindings.begin(),
        record.bindings.end(),
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
      original_upper != manifest.bindings.end() &&
        changed_upper != changed_bound_compilation.problem->manifest().bindings.end() &&
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
      manifest.bindings.begin(),
      manifest.bindings.end(),
      [](const compiler::v1::CompiledBindingRecord &record) {
        return record.role == semantic::v1::DataRole::diffusion;
      });
    const auto lower_binding = std::find_if(
      manifest.bindings.begin(),
      manifest.bindings.end(),
      [](const compiler::v1::CompiledBindingRecord &record) {
        return record.role == semantic::v1::DataRole::lower_bound;
      });
    const auto upper_binding = std::find_if(
      manifest.bindings.begin(),
      manifest.bindings.end(),
      [](const compiler::v1::CompiledBindingRecord &record) {
        return record.role == semantic::v1::DataRole::upper_bound;
      });
    contract::require(
      diffusion_binding != manifest.bindings.end() &&
        diffusion_binding->field_shape ==
          compiler::v1::CompiledFieldShape::scalar_constant &&
        diffusion_binding->scalar_value.has_value() &&
        *diffusion_binding->scalar_value == 1.0 &&
        !diffusion_binding->value_digest.empty() &&
        lower_binding != manifest.bindings.end() &&
        lower_binding->scalar_value.has_value() &&
        *lower_binding->scalar_value == -1.0 &&
        upper_binding != manifest.bindings.end() &&
        upper_binding->scalar_value.has_value() &&
        *upper_binding->scalar_value == 0.05 &&
        manifest.resolved_decision.target_id.find("compiled_target:") == 0,
      "v1 compiler did not retain lossless scalar binding provenance");
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
        {"l2_dirichlet_transposition",
         "nmopt.dealii.l2_dirichlet_transposition",
         {"dealii", "compiler"},
         60,
         []() { run_l2_dirichlet_transposition_lowering_test<2>(); }},
        {"partial_dirichlet_control",
         "nmopt.dealii.partial_dirichlet_control",
         {"dealii", "compiler"},
         60,
         []() { run_partial_dirichlet_control_contract_test<2>(); }},
        {"subdomain_observation",
         "nmopt.dealii.subdomain_observation",
         {"dealii", "compiler"},
         60,
         []() { run_subdomain_observation_contract_test<2>(); }},
        {"point_sensor",
         "nmopt.dealii.point_sensor",
         {"dealii", "compiler", "observation"},
         60,
         []() { run_point_sensor_contract_test<2>(); }},
        {"normal_flux",
         "nmopt.dealii.normal_flux",
         {"dealii", "compiler", "observation"},
         60,
         []() { run_normal_flux_contract_test<2>(); }},
        {"h1_state_observation",
         "nmopt.dealii.h1_state_observation",
         {"dealii", "compiler"},
         60,
         []() { run_h1_state_observation_contract_test<2>(); }},
        {"neumann_boundary",
         "nmopt.dealii.neumann_boundary",
         {"dealii", "compiler"},
         60,
         []() { run_neumann_boundary_contract_test<2>(); }},
        {"neumann_convection_subdomain",
         "nmopt.dealii.neumann_convection_subdomain",
         {"dealii", "compiler"},
         60,
         []() { run_neumann_convection_subdomain_contract_test<2>(); }},
        {"weighted_boundary_trace",
         "nmopt.dealii.weighted_boundary_trace",
         {"dealii", "compiler"},
         60,
         []() { run_weighted_boundary_trace_contract_test<2>(); }},
        {"h1_control",
         "nmopt.dealii.h1_control",
         {"dealii", "compiler"},
         60,
         []() { run_h1_control_regularisation_contract_test<2>(); }},
        {"hminus1_metric",
         "nmopt.dealii.hminus1_metric",
         {"dealii", "contract", "metric"},
         30,
         run_hminus1_metric_contract_test},
        {"continuous_control_components",
         "nmopt.dealii.continuous_control_components",
         {"dealii", "compiler", "metric"},
         30,
         []() { run_continuous_control_component_contract_test<2>(); }},
        {"hminus1_compilation",
         "nmopt.dealii.hminus1_compilation",
         {"dealii", "compiler", "metric"},
         60,
         []() { run_hminus1_compilation_contract_test<2>(); }},
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
