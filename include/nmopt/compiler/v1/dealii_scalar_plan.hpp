#pragma once

#include "nmopt/contract/linalg.hpp"
#include "nmopt/semantic/v1/resolved_problem.hpp"

#include <algorithm>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace nmopt::compiler::v1
{
  enum class ScalarResidualOperatorKind
  {
    diffusion_reaction,
    tensor_diffusion,
    conservative_transport,
    advective_transport,
    reaction,
    volume_source,
    volume_control,
    robin_bilinear,
    robin_source
  };

  enum class ScalarObservationOperatorKind
  {
    volume_restriction,
    h1_state_restriction,
    point_sensor
  };

  enum class ScalarLossOperatorKind
  {
    quadratic_tracking,
    quadratic_control_regularisation
  };

  enum class ScalarMetricOperatorKind
  {
    cellwise_l2
  };

  enum class ScalarConstraintOperatorKind
  {
    none,
    cellwise_box
  };

  enum class ScalarTransformationOperatorKind
  {
    none,
    fixed_dirichlet_reconstruction
  };

  struct ScalarResidualContribution
  {
    std::string                component_id;
    ScalarResidualOperatorKind operator_kind;
    std::vector<std::string>   variable_ids;
    std::vector<std::string>   data_ids;
    std::string                region_id;
    std::string                handler_id;
  };

  struct ScalarObservationContribution
  {
    std::string                   component_id;
    ScalarObservationOperatorKind operator_kind;
    std::string                   input_variable_id;
    std::string                   output_space_id;
    std::string                   region_id;
    std::string                   handler_id;
  };

  struct ScalarLossContribution
  {
    std::string            component_id;
    ScalarLossOperatorKind operator_kind;
    std::string            observation_id;
    std::string            data_id;
    std::string            handler_id;
  };

  struct ScalarLoweringPlan
  {
    std::string                              semantic_problem_id;
    std::string                              state_variable_id;
    std::string                              decision_variable_id;
    std::string                              equation_id;
    std::vector<ScalarResidualContribution>  residual_terms;
    std::vector<ScalarObservationContribution> observations;
    std::vector<ScalarLossContribution>       losses;
    ScalarMetricOperatorKind                  metric =
      ScalarMetricOperatorKind::cellwise_l2;
    ScalarConstraintOperatorKind              constraint =
      ScalarConstraintOperatorKind::none;
    ScalarTransformationOperatorKind          transformation =
      ScalarTransformationOperatorKind::none;
    std::string                               metric_id;
    std::string                               constraint_id;
    std::string                               fixed_data_id;
    std::set<unsigned int>                    dirichlet_boundary_ids;
    std::set<unsigned int>                    robin_boundary_ids;
    bool                                      tracking_full_domain = true;
    std::set<unsigned int>                    tracking_material_ids;
    std::string                               point_sensor_region_id;
    std::vector<std::vector<double>>          point_sensor_coordinates;
    std::string                               point_sensor_evaluation_policy;
    std::vector<std::string>                  provenance;
  };

  class ScalarResidualTermHandler final
  {
  public:
    ScalarResidualTermHandler(const semantic::v1::ResidualTermKind semantic_kind,
                              const ScalarResidualOperatorKind operator_kind,
                              std::string handler_id)
      : semantic_kind_(semantic_kind)
      , operator_kind_(operator_kind)
      , handler_id_(std::move(handler_id))
    {}

    semantic::v1::ResidualTermKind
    semantic_kind() const
    {
      return semantic_kind_;
    }

    void
    contribute(const semantic::v1::ResidualTermSpec &term,
               ScalarLoweringPlan &                  plan) const
    {
      contract::require(term.kind == semantic_kind_,
                        "Scalar residual handler received the wrong semantic kind");
      plan.residual_terms.push_back({term.id,
                                     operator_kind_,
                                     term.variable_ids,
                                     term.data_ids,
                                     term.region_id,
                                     handler_id_});
      plan.provenance.push_back(term.id + " <- " + handler_id_);
    }

  private:
    semantic::v1::ResidualTermKind semantic_kind_;
    ScalarResidualOperatorKind     operator_kind_;
    std::string                    handler_id_;
  };

  class ScalarObservationHandler final
  {
  public:
    ScalarObservationHandler(const semantic::v1::ObservationKind semantic_kind,
                             const ScalarObservationOperatorKind operator_kind,
                             std::string handler_id)
      : semantic_kind_(semantic_kind)
      , operator_kind_(operator_kind)
      , handler_id_(std::move(handler_id))
    {}

    semantic::v1::ObservationKind
    semantic_kind() const
    {
      return semantic_kind_;
    }

    void
    contribute(const semantic::v1::ObservationSpec &observation,
               ScalarLoweringPlan &                  plan) const
    {
      contract::require(
        observation.kind == semantic_kind_,
        "Scalar observation handler received the wrong semantic kind");
      plan.observations.push_back(
        {observation.id,
         operator_kind_,
         observation.input_variable_id,
         observation.output_space_id,
         observation.region_id,
         handler_id_});
      plan.provenance.push_back(observation.id + " <- " + handler_id_);
    }

  private:
    semantic::v1::ObservationKind semantic_kind_;
    ScalarObservationOperatorKind operator_kind_;
    std::string                   handler_id_;
  };

  class ScalarLossHandler final
  {
  public:
    ScalarLossHandler(const semantic::v1::LossKind semantic_kind,
                      const ScalarLossOperatorKind operator_kind,
                      std::string handler_id)
      : semantic_kind_(semantic_kind)
      , operator_kind_(operator_kind)
      , handler_id_(std::move(handler_id))
    {}

    semantic::v1::LossKind
    semantic_kind() const
    {
      return semantic_kind_;
    }

    void
    contribute(const semantic::v1::LossSpec &loss,
               ScalarLoweringPlan &          plan) const
    {
      contract::require(loss.kind == semantic_kind_,
                        "Scalar loss handler received the wrong semantic kind");
      plan.losses.push_back({loss.id,
                             operator_kind_,
                             loss.source_observation_id,
                             loss.data_id,
                             handler_id_});
      plan.provenance.push_back(loss.id + " <- " + handler_id_);
    }

  private:
    semantic::v1::LossKind semantic_kind_;
    ScalarLossOperatorKind operator_kind_;
    std::string             handler_id_;
  };

  class ScalarMetricHandler final
  {
  public:
    ScalarMetricHandler(const semantic::v1::MetricKind semantic_kind,
                        const ScalarMetricOperatorKind operator_kind,
                        std::string handler_id)
      : semantic_kind_(semantic_kind)
      , operator_kind_(operator_kind)
      , handler_id_(std::move(handler_id))
    {}

    semantic::v1::MetricKind
    semantic_kind() const
    {
      return semantic_kind_;
    }

    void
    contribute(const semantic::v1::MetricSpec &metric,
               ScalarLoweringPlan &            plan) const
    {
      plan.metric = operator_kind_;
      plan.provenance.push_back(metric.id + " <- " + handler_id_);
    }

  private:
    semantic::v1::MetricKind semantic_kind_;
    ScalarMetricOperatorKind operator_kind_;
    std::string               handler_id_;
  };

  class ScalarConstraintHandler final
  {
  public:
    ScalarConstraintHandler(
      const semantic::v1::ConstraintKind semantic_kind,
      const ScalarConstraintOperatorKind operator_kind,
      std::string handler_id)
      : semantic_kind_(semantic_kind)
      , operator_kind_(operator_kind)
      , handler_id_(std::move(handler_id))
    {}

    semantic::v1::ConstraintKind
    semantic_kind() const
    {
      return semantic_kind_;
    }

    void
    contribute(const semantic::v1::ConstraintSpec &constraint,
               ScalarLoweringPlan &                plan) const
    {
      plan.constraint = operator_kind_;
      plan.provenance.push_back(constraint.id + " <- " + handler_id_);
    }

  private:
    semantic::v1::ConstraintKind semantic_kind_;
    ScalarConstraintOperatorKind operator_kind_;
    std::string                   handler_id_;
  };

  class ScalarTransformationHandler final
  {
  public:
    ScalarTransformationHandler(
      const semantic::v1::TransformationKind semantic_kind,
      const ScalarTransformationOperatorKind operator_kind,
      std::string handler_id)
      : semantic_kind_(semantic_kind)
      , operator_kind_(operator_kind)
      , handler_id_(std::move(handler_id))
    {}

    semantic::v1::TransformationKind
    semantic_kind() const
    {
      return semantic_kind_;
    }

    void
    contribute(const semantic::v1::TransformationSpec &transformation,
               ScalarLoweringPlan &                    plan) const
    {
      plan.transformation = operator_kind_;
      plan.fixed_data_id = transformation.fixed_data_id;
      plan.provenance.push_back(transformation.id + " <- " + handler_id_);
    }

  private:
    semantic::v1::TransformationKind semantic_kind_;
    ScalarTransformationOperatorKind operator_kind_;
    std::string                       handler_id_;
  };

  class DealiiScalarLowererRegistryV1 final
  {
  public:
    DealiiScalarLowererRegistryV1()
      : residual_handlers_{
          {semantic::v1::ResidualTermKind::diffusion_reaction,
           ScalarResidualOperatorKind::diffusion_reaction,
           "dealii.scalar.residual.diffusion_reaction"},
          {semantic::v1::ResidualTermKind::tensor_diffusion,
           ScalarResidualOperatorKind::tensor_diffusion,
           "dealii.scalar.residual.tensor_diffusion"},
          {semantic::v1::ResidualTermKind::conservative_transport,
           ScalarResidualOperatorKind::conservative_transport,
           "dealii.scalar.residual.conservative_transport"},
          {semantic::v1::ResidualTermKind::advective_transport,
           ScalarResidualOperatorKind::advective_transport,
           "dealii.scalar.residual.advective_transport"},
          {semantic::v1::ResidualTermKind::reaction,
           ScalarResidualOperatorKind::reaction,
           "dealii.scalar.residual.reaction"},
          {semantic::v1::ResidualTermKind::volume_source,
           ScalarResidualOperatorKind::volume_source,
           "dealii.scalar.residual.volume_source"},
          {semantic::v1::ResidualTermKind::volume_control,
           ScalarResidualOperatorKind::volume_control,
           "dealii.scalar.residual.volume_control"},
          {semantic::v1::ResidualTermKind::robin_bilinear,
           ScalarResidualOperatorKind::robin_bilinear,
           "dealii.scalar.residual.robin_bilinear"},
          {semantic::v1::ResidualTermKind::robin_source,
           ScalarResidualOperatorKind::robin_source,
           "dealii.scalar.residual.robin_source"}}
      , observation_handlers_{
          {semantic::v1::ObservationKind::volume_restriction,
           ScalarObservationOperatorKind::volume_restriction,
           "dealii.scalar.observation.volume_restriction"},
          {semantic::v1::ObservationKind::h1_state_restriction,
           ScalarObservationOperatorKind::h1_state_restriction,
           "dealii.scalar.observation.h1_state_restriction"},
          {semantic::v1::ObservationKind::point_sensor,
           ScalarObservationOperatorKind::point_sensor,
           "dealii.scalar.observation.point_sensor"}}
      , loss_handlers_{
          {semantic::v1::LossKind::quadratic_tracking,
           ScalarLossOperatorKind::quadratic_tracking,
           "dealii.scalar.loss.quadratic_tracking"},
          {semantic::v1::LossKind::quadratic_control_regularisation,
           ScalarLossOperatorKind::quadratic_control_regularisation,
           "dealii.scalar.loss.quadratic_control_regularisation"}}
      , metric_handlers_{
          {semantic::v1::MetricKind::l2,
           ScalarMetricOperatorKind::cellwise_l2,
           "dealii.scalar.metric.cellwise_l2"}}
      , constraint_handlers_{
          {semantic::v1::ConstraintKind::cellwise_box,
           ScalarConstraintOperatorKind::cellwise_box,
           "dealii.scalar.constraint.cellwise_box"}}
      , transformation_handlers_{
          {semantic::v1::TransformationKind::fixed_dirichlet_reconstruction,
           ScalarTransformationOperatorKind::fixed_dirichlet_reconstruction,
           "dealii.scalar.transformation.fixed_dirichlet"}}
    {}

    const ScalarResidualTermHandler *
    residual_handler(const semantic::v1::ResidualTermKind kind) const
    {
      return find_handler(residual_handlers_, kind);
    }

    const ScalarObservationHandler *
    observation_handler(const semantic::v1::ObservationKind kind) const
    {
      return find_handler(observation_handlers_, kind);
    }

    const ScalarLossHandler *
    loss_handler(const semantic::v1::LossKind kind) const
    {
      return find_handler(loss_handlers_, kind);
    }

    const ScalarMetricHandler *
    metric_handler(const semantic::v1::MetricKind kind) const
    {
      return find_handler(metric_handlers_, kind);
    }

    const ScalarConstraintHandler *
    constraint_handler(const semantic::v1::ConstraintKind kind) const
    {
      return find_handler(constraint_handlers_, kind);
    }

    const ScalarTransformationHandler *
    transformation_handler(const semantic::v1::TransformationKind kind) const
    {
      return find_handler(transformation_handlers_, kind);
    }

  private:
    template <typename Handler, typename Kind>
    static const Handler *
    find_handler(const std::vector<Handler> &handlers, const Kind kind)
    {
      const auto handler = std::find_if(
        handlers.begin(), handlers.end(), [kind](const Handler &candidate) {
          return candidate.semantic_kind() == kind;
        });
      return handler == handlers.end() ? nullptr : &*handler;
    }

    std::vector<ScalarResidualTermHandler> residual_handlers_;
    std::vector<ScalarObservationHandler>  observation_handlers_;
    std::vector<ScalarLossHandler>         loss_handlers_;
    std::vector<ScalarMetricHandler>       metric_handlers_;
    std::vector<ScalarConstraintHandler>   constraint_handlers_;
    std::vector<ScalarTransformationHandler> transformation_handlers_;
  };

  struct ScalarPlanResult
  {
    semantic::v1::ValidationReport diagnostics;
    std::optional<ScalarLoweringPlan> plan;

    bool
    succeeded() const
    {
      return diagnostics.valid() && plan.has_value();
    }
  };

  class DealiiScalarLoweringPlanner final
  {
  public:
    explicit DealiiScalarLoweringPlanner(
      DealiiScalarLowererRegistryV1 registry = {})
      : registry_(std::move(registry))
    {}

    ScalarPlanResult
    plan(const semantic::v1::ResolvedProblemView &problem) const
    {
      using semantic::v1::DiagnosticCategory;
      ScalarPlanResult result;
      ScalarLoweringPlan plan;
      const auto &specification = problem.specification();
      plan.semantic_problem_id = specification.id;
      plan.state_variable_id = specification.formulation.state_variable_id;
      plan.decision_variable_id = specification.formulation.control_variable_id;
      plan.equation_id = specification.formulation.equation_id;
      plan.metric_id = specification.formulation.metric_id;
      plan.constraint_id = specification.formulation.constraint_id;

      const auto &equation = problem.equation(plan.equation_id);
      for (const auto &term_id : equation.residual_term_ids)
        {
          const auto &term = problem.residual_term(term_id);
          const auto *handler = registry_.residual_handler(term.kind);
          if (handler == nullptr)
            result.diagnostics.add(
              DiagnosticCategory::lowerability,
              term.id,
              "scalar_residual_component_lowerer",
              "Register this residual contribution in the bounded scalar lowering plan.");
          else
            handler->contribute(term, plan);
        }
      for (const auto &observation : specification.observations)
        {
          const auto *handler = registry_.observation_handler(observation.kind);
          if (handler == nullptr)
            result.diagnostics.add(
              DiagnosticCategory::lowerability,
              observation.id,
              "scalar_observation_component_lowerer",
              "Register this observation in the bounded scalar lowering plan.");
          else
            {
              handler->contribute(observation, plan);
              if (observation.kind == semantic::v1::ObservationKind::point_sensor)
                {
                  const auto &region = problem.region(observation.region_id);
                  plan.point_sensor_region_id = region.id;
                  plan.point_sensor_coordinates = region.point_coordinates;
                  plan.point_sensor_evaluation_policy =
                    "FE_Q shape evaluation at immutable physical points with assembled transpose C_h^T";
                }
            }
        }
      for (const auto &loss : specification.losses)
        {
          const auto *handler = registry_.loss_handler(loss.kind);
          if (handler == nullptr)
            result.diagnostics.add(
              DiagnosticCategory::lowerability,
              loss.id,
              "scalar_loss_component_lowerer",
              "Register this loss in the bounded scalar lowering plan.");
          else
            handler->contribute(loss, plan);
        }

      const auto &metric = problem.metric(plan.metric_id);
      const auto *metric_handler = registry_.metric_handler(metric.kind);
      if (metric_handler == nullptr)
        result.diagnostics.add(
          DiagnosticCategory::lowerability,
          metric.id,
          "scalar_metric_component_lowerer",
          "Select the registered cellwise L2 scalar metric.");
      else
        metric_handler->contribute(metric, plan);

      if (!plan.constraint_id.empty())
        {
          const auto &constraint = problem.constraint(plan.constraint_id);
          const auto *constraint_handler =
            registry_.constraint_handler(constraint.kind);
          if (constraint_handler == nullptr)
            result.diagnostics.add(
              DiagnosticCategory::lowerability,
              constraint.id,
              "scalar_constraint_component_lowerer",
              "Select the registered cellwise scalar box.");
          else
            constraint_handler->contribute(constraint, plan);
        }

      const auto &state = problem.variable(plan.state_variable_id);
      if (!state.physical_field_transform_id.empty())
        {
          const auto &transformation =
            problem.transformation(state.physical_field_transform_id);
          const auto *transformation_handler =
            registry_.transformation_handler(transformation.kind);
          if (transformation_handler == nullptr)
            result.diagnostics.add(
              DiagnosticCategory::lowerability,
              transformation.id,
              "scalar_transformation_component_lowerer",
              "Select no transformation or the registered fixed-Dirichlet reconstruction.");
          else
            transformation_handler->contribute(transformation, plan);
        }

      for (const auto &requirement : specification.requirement_policies)
        if (requirement.subject_id == plan.state_variable_id &&
            requirement.kind == semantic::v1::RequirementKind::fixed_dirichlet)
          {
            const auto &region = problem.region(requirement.region_id);
            plan.dirichlet_boundary_ids.insert(region.boundary_ids.begin(),
                                               region.boundary_ids.end());
          }

      for (const auto &term : plan.residual_terms)
        if (term.operator_kind == ScalarResidualOperatorKind::robin_bilinear ||
            term.operator_kind == ScalarResidualOperatorKind::robin_source)
          {
            const auto &region = problem.region(term.region_id);
            plan.robin_boundary_ids.insert(region.boundary_ids.begin(),
                                           region.boundary_ids.end());
          }

      const auto tracking_loss = std::find_if(
        specification.losses.begin(),
        specification.losses.end(),
        [](const semantic::v1::LossSpec &loss) {
          return loss.kind == semantic::v1::LossKind::quadratic_tracking;
        });
      if (tracking_loss != specification.losses.end())
        {
          const auto &observation =
            problem.observation(tracking_loss->source_observation_id);
          const auto &region = problem.region(observation.region_id);
          plan.tracking_full_domain = region.is_full_domain;
          plan.tracking_material_ids.insert(region.material_ids.begin(),
                                            region.material_ids.end());
        }

      if (plan.dirichlet_boundary_ids.empty())
        result.diagnostics.add(
          DiagnosticCategory::lowerability,
          plan.state_variable_id,
          "scalar_fixed_dirichlet_plan",
          "Select a non-empty fixed-Dirichlet boundary for this scalar plan.");
      if (result.diagnostics.valid())
        result.plan = std::move(plan);
      return result;
    }

  private:
    DealiiScalarLowererRegistryV1 registry_;
  };
} // namespace nmopt::compiler::v1
