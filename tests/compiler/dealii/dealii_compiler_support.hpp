#pragma once

#include "nmopt/contract/reduced_dto.hpp"
#include "nmopt/contract/supplied_otd_kkt.hpp"
#include "nmopt/compiler/v1/dealii_neumann_control_realisation.hpp"
#include "nmopt/compiler/v1/dealii_compiler.hpp"
#include "nmopt/compiler/v1/dealii_volume_observation.hpp"
#include "nmopt/dealii/cellwise_box_constraint.hpp"
#include "nmopt/dealii/hminus1_metric.hpp"
#include "nmopt/dealii/quadratic_form.hpp"
#include "nmopt/dealii/serial_kkt_solver.hpp"
#include "nmopt/semantic/v1/problem_spec.hpp"
#include "nmopt/solvers/reduced_gradient.hpp"
#include "../../support/contract_errors.hpp"
#include "../../support/diagnostics.hpp"
#include "../../support/manifest_contracts.hpp"
#include "../../support/scenario_dispatch.hpp"
#include "../../support/scoped_temporary_directory.hpp"

#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/numerics/vector_tools.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
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

  void
  require_compiled_hessian_evidence(
    const compiler::v1::CompiledProblemT<Backend> &problem,
    const std::string                            &target)
  {
    const auto &model = problem.executable_model();
    const auto *hessian = problem.reduced_hessian();
    contract::require(hessian != nullptr,
                      target + " compiled target omitted its Hessian capability");

    const auto reduced = problem.make_reduced_dto();
    dealii::Vector<double> control_values(model.variable_layout()->dimension(1));
    const Primal control(model.variable_layout()->single_block(1, "control"),
                         {std::move(control_values)});

    dealii::Vector<double> direction_values(control.layout()->dimension(0));
    dealii::Vector<double> second_direction_values(
      control.layout()->dimension(0));
    for (dealii::types::global_dof_index index = 0;
         index < direction_values.size();
         ++index)
      {
        direction_values[index] =
          (index % 2 == 0 ? 0.04 : -0.03) * static_cast<double>(index + 1);
        second_direction_values[index] =
          (index % 3 == 0 ? -0.02 : 0.03) * static_cast<double>(index + 1);
      }
    const Primal direction(control.layout(), {std::move(direction_values)});
    const Primal second_direction(
      control.layout(), {std::move(second_direction_values)});
    const Covector hessian_action = hessian->apply(control, direction);
    const Covector second_hessian_action =
      hessian->apply(control, second_direction);
    require_close(contract::pair(hessian_action, second_direction),
                  contract::pair(second_hessian_action, direction),
                  1e-10,
                  target + " compiled Hessian symmetry");

    constexpr double hessian_step = 1e-5;
    const Covector reduced_derivative_plus =
      reduced.evaluate(shifted(control, direction, hessian_step))
        .reduced_derivative;
    const Covector reduced_derivative_minus =
      reduced.evaluate(shifted(control, direction, -hessian_step))
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
                  target + " compiled Hessian finite-difference action");
  }




  void
  require_constraint_realisation(
    const compiler::v1::CompilationManifest &manifest,
    const std::string &                       expected,
    const std::string &                       target)
  {
    const auto &compatibility = manifest.compatibility;
    contract::require(
      !compatibility.compiler_id.empty() && !compatibility.backend.empty() &&
        compatibility.execution == manifest.resolved_decision.execution_id &&
        compatibility.dual_representation ==
          manifest.resolved_decision.formulation_record.dual_representation &&
        compatibility.metric_solve_policy.find(
          manifest.resolved_decision.metric_record.realisation_id) !=
          std::string::npos &&
        compatibility.constraint_realisation == expected,
      target + " manifest compatibility view is incomplete");
    std::string structured_expected = "none";
    if (expected.find("l2_cellwise_parameter") != std::string::npos)
      structured_expected = "l2_cellwise_parameter";
    else if (expected.find("l2_cellwise") != std::string::npos)
      structured_expected = "l2_cellwise";
    else if (expected.find("l2_facewise") != std::string::npos)
      structured_expected = "l2_facewise";
    contract::require(
      manifest.compatibility.constraint_realisation == expected &&
        manifest.resolved_decision.constraint_record.realisation_id == structured_expected &&
        manifest.resolved_decision.constraint_record.present == (structured_expected != "none") &&
        (structured_expected == "none" ||
         manifest.resolved_decision.constraint_record.projection_metric_id == structured_expected),
      target + " manifest constraint realization: expected " + expected +
        ", got " + manifest.compatibility.constraint_realisation);
    const auto has_runtime_role = [&manifest](const std::string &role) {
      return std::any_of(
        manifest.resolved_decision.spaces.begin(),
        manifest.resolved_decision.spaces.end(),
        [&role](const compiler::v1::CompiledSpaceRecord &space) {
          return space.runtime_role == role;
        });
    };
    const auto has_realized_map = [&manifest](const std::string &id) {
      return std::any_of(
        manifest.resolved_decision.realized_maps.begin(),
        manifest.resolved_decision.realized_maps.end(),
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
        manifest.resolved_decision.realized_maps.begin(),
        manifest.resolved_decision.realized_maps.end(),
        [](const compiler::v1::CompiledRealizedMapRecord &map) {
          return map.realization_id != "coefficient_restriction" ||
                 (map.input_dimensions.size() == 1 &&
                  map.input_dimensions.front() == map.output_dimension);
        });
    const auto observation_map_dimensions_match_spaces =
      std::all_of(
        manifest.resolved_decision.realized_maps.begin(),
        manifest.resolved_decision.realized_maps.end(),
        [&manifest](const compiler::v1::CompiledRealizedMapRecord &map) {
          if (std::find(manifest.compatibility.observation_ids.begin(),
                        manifest.compatibility.observation_ids.end(),
                        map.semantic_id) == manifest.compatibility.observation_ids.end())
            return true;
          const auto space = std::find_if(
            manifest.resolved_decision.spaces.begin(),
            manifest.resolved_decision.spaces.end(),
            [&map](const compiler::v1::CompiledSpaceRecord &candidate) {
              return candidate.semantic_id == map.output_space_id;
            });
          return space != manifest.resolved_decision.spaces.end() &&
                 space->dimension == map.output_dimension;
        });
    contract::require(
      std::all_of(manifest.compatibility.observation_ids.begin(),
                  manifest.compatibility.observation_ids.end(),
                  has_realized_map) &&
        std::all_of(manifest.compatibility.transformation_ids.begin(),
                    manifest.compatibility.transformation_ids.end(),
                    has_realized_map) &&
        control_restrictions_are_coefficient_maps &&
        observation_map_dimensions_match_spaces,
      target + " manifest omitted a realized observation or transformation map");
    contract::require(
      manifest.schema_version == 4 &&
        manifest.resolved_decision.formulation_record.kind ==
          semantic::v1::FormulationKind::reduced_dto &&
        manifest.resolved_decision.formulation_record.provenance ==
          semantic::v1::FormulationProvenance::dto &&
        manifest.resolved_decision.formulation_record.execution ==
          compiler::v1::ExecutionRealisation::assembled &&
        manifest.resolved_decision.mesh_record.dimension > 0 &&
        manifest.resolved_decision.mesh_record.active_cells > 0 &&
        !manifest.resolved_decision.mesh_record.provenance.empty() &&
        has_runtime_role("state") &&
        has_runtime_role("test_and_adjoint") &&
        (has_runtime_role("decision_control") ||
         has_runtime_role("decision_parameter")) &&
        has_runtime_role("observation") &&
        !manifest.resolved_decision.bindings.empty() &&
        manifest.resolved_decision.state_solve_record.maximum_iterations > 0 &&
        manifest.resolved_decision.adjoint_solve_record.maximum_iterations > 0 &&
        !manifest.resolved_decision.metric_record.semantic_id.empty() &&
        !manifest.resolved_decision.metric_record.realisation_id.empty() &&
        !manifest.resolved_decision.mesh_record.structural_identity.empty() &&
        manifest.resolved_decision.formulation_id ==
          manifest.resolved_decision.formulation_record.semantic_id &&
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
      decision.formulation_id == decision.formulation_record.semantic_id &&
        decision.formulation_record.dual_representation ==
          manifest.compatibility.dual_representation &&
        manifest.compatibility.execution == decision.execution_id &&
        manifest.compatibility.metric_solve_policy.find(
          decision.metric_record.realisation_id) != std::string::npos &&
        manifest.compatibility.state_adjoint_solve_policy.find(
          decision.state_solve_record.operator_realisation) !=
          std::string::npos &&
        !decision.spaces.empty() && !decision.bindings.empty() &&
        !decision.realized_maps.empty() && !decision.realized_spaces.empty(),
      description + " was not rendered from its resolved typed records");
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
    const auto specification =
      semantic::v1::make_scalar_diffusion_reaction_problem(true);
    compiler::v1::DealiiDiscretisationPolicy policy;
    policy.state_degree = 1;
    const compiler::v1::DealiiDataBindings<dim> data_bindings{
      forcing,
      desired_state,
      1.0,
      0.5,
      0.1,
      test_binding_provenance("weak_form_oracle")};
    const compiler::v1::CellwiseBoxDataBindings bounds{
      compiler::v1::CellwiseBoundValue{-1.0},
      compiler::v1::CellwiseBoundValue{2.0}};
    const compiler::v1::DealiiCompiler compiler;
    const auto compilation = compiler.compile(specification,
                                              triangulation,
                                              data_bindings,
                                              policy,
                                              bounds);
    contract::require(compilation.succeeded() && compilation.problem,
                      "the scalar component oracle compilation failed");
    const auto &model = compilation.problem->executable_model();

    const std::size_t state_size = model.variable_layout()->dimension(0);
    contract::require(state_size == 1,
                      "the scalar component oracle needs one independent state DoF");
    constexpr std::size_t interior_state_dof = 0;

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


} // namespace
