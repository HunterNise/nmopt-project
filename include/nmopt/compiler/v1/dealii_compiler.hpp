#pragma once

#include "nmopt/compiler/v1/compiled_problem.hpp"
#include "nmopt/compiler/v1/dealii_capabilities.hpp"
#include "nmopt/compiler/v1/dealii_coefficient_identification.hpp"
#include "nmopt/compiler/v1/dealii_dirichlet_control.hpp"
#include "nmopt/compiler/v1/dealii_fixed_dirichlet.hpp"
#include "nmopt/compiler/v1/dealii_continuous_control.hpp"
#include "nmopt/compiler/v1/dealii_neumann_boundary.hpp"
#include "nmopt/compiler/v1/dealii_scalar_plan.hpp"
#include "nmopt/compiler/v1/dealii_types.hpp"
#include "nmopt/contract/supplied_otd_kkt.hpp"
#include "nmopt/dealii/facewise_box_constraint.hpp"
#include "nmopt/dealii/scalar_diffusion_reaction.hpp"
#include "nmopt/semantic/v1/validation.hpp"

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/base/geometry_info.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <limits>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace nmopt::compiler::v1
{
  class DealiiCompiler final
  {
  public:
    explicit DealiiCompiler(
      DealiiCapabilityRegistryV1 capabilities = {},
      DealiiScalarLowererRegistryV1 scalar_registry = {})
      : capabilities_(std::move(capabilities))
      , scalar_planner_(std::move(scalar_registry))
    {}

    semantic::v1::ValidationReport
    validate(const semantic::v1::ProblemSpec &  specification,
             const DealiiDiscretisationPolicy & policy,
             const CompilationProduct product = CompilationProduct::reduced_dto) const
    {
      auto resolution = semantic::v1::SemanticResolver().resolve(specification);
      semantic::v1::ValidationReport report = std::move(resolution.diagnostics);
      if (!report.valid())
        return report;
      auto request = resolve_compilation_request(*resolution.problem);
      const auto dirichlet_registration =
        resolve_dirichlet_control_registration(*resolution.problem, request);
      close_compilation_request(*resolution.problem,
                                request,
                                dirichlet_registration);
      validate_lowerability(specification, request, policy, report);
      validate_formulation_capability(specification, report);
      validate_supplied_otd_capability(
        specification, request, product, report);
      validate_dirichlet_control_registration(*resolution.problem,
                                               request,
                                               report);
      validate_product_capability(
        specification, request, policy, product, report);
      return report;
    }

    template <int dim>
    CompilationResultT<dealii_backend::SerialBackend>
    compile(const semantic::v1::ProblemSpec &  specification,
            dealii::Triangulation<dim> &        triangulation,
            const DealiiDataBindings<dim> &     data,
            const DealiiDiscretisationPolicy &  policy = {},
            std::optional<CellwiseBoxDataBindings> bounds = std::nullopt,
            std::optional<FacewiseBoxDataBindings> facewise_bounds = std::nullopt,
            const CompilationProduct product = CompilationProduct::reduced_dto) const
    {
      return compile_impl(specification,
                          triangulation,
                          data,
                          policy,
                          std::move(bounds),
                          std::move(facewise_bounds),
                          {},
                          "caller-owned triangulation",
                          false,
                          product);
    }

    template <int dim>
    CompilationResultT<dealii_backend::SerialBackend>
    compile(
      const semantic::v1::ProblemSpec &specification,
      const std::shared_ptr<DealiiCompilationSession<dim>> &session,
      const DealiiDataBindings<dim> &data,
      const DealiiDiscretisationPolicy &policy = {},
      std::optional<CellwiseBoxDataBindings> bounds = std::nullopt,
      std::optional<FacewiseBoxDataBindings> facewise_bounds = std::nullopt,
      const CompilationProduct product = CompilationProduct::reduced_dto) const
    {
      CompilationResultT<dealii_backend::SerialBackend> result;
      if (!session)
        {
          result.diagnostics.add(
            semantic::v1::DiagnosticCategory::lowerability,
            specification.id,
            "compilation_session_presence",
            "Supply a non-null owned deal.II compilation session.");
          return result;
        }
      return compile_impl(specification,
                          session->mutable_triangulation(),
                          data,
                          policy,
                          std::move(bounds),
                          std::move(facewise_bounds),
                          session,
                          session->mesh_provenance(),
                          true,
                          product);
    }

  private:
    static void
    append_data_binding_request(
      ResolvedCompilationRequest &             request,
      const semantic::v1::ResolvedProblemView &resolved,
      const std::string &                       data_id,
      const ResolvedBindingPort                port,
      const std::string &                       evaluation_realisation,
      const std::string &                       runtime_representation,
      const bool                                required)
    {
      if (data_id.empty())
        return;
      const auto &datum = resolved.datum(data_id);

      std::string region_id;
      if (!datum.space_id.empty())
        region_id = resolved.space(datum.space_id).region_id;
      request.data_bindings.push_back({port,
                                       datum.id,
                                       datum.role,
                                       datum.kind,
                                       datum.space_id,
                                       region_id,
                                       evaluation_realisation,
                                       runtime_representation,
                                       required});
    }

    static std::string
    data_id_for_term_role(const semantic::v1::ResolvedProblemView &resolved,
                          const semantic::v1::ResidualTermKind             kind,
                          const semantic::v1::DataRole                     role)
    {
      for (const auto &term : resolved.specification().residual_terms)
        if (term.kind == kind)
          for (const auto &data_id : term.data_ids)
            if (resolved.datum(data_id).role == role)
              return data_id;
      return {};
    }

    static std::string
    data_id_for_loss_role(const semantic::v1::ResolvedProblemView &resolved,
                          const semantic::v1::DataRole             role)
    {
      for (const auto &loss : resolved.specification().losses)
        if (resolved.datum(loss.data_id).role == role)
          return loss.data_id;
      return {};
    }

    static ResolvedCompilationRequest
    resolve_compilation_request(
      const semantic::v1::ResolvedProblemView &resolved)
    {
      const auto &specification = resolved.specification();
      ResolvedCompilationRequest request;
      request.semantic_problem_id = specification.id;
      request.transposition_diffusion_data_id = data_id_for_term_role(
        resolved,
        semantic::v1::ResidualTermKind::diffusion_reaction,
        semantic::v1::DataRole::diffusion);
      request.transposition_reaction_data_id = data_id_for_term_role(
        resolved,
        semantic::v1::ResidualTermKind::diffusion_reaction,
        semantic::v1::DataRole::reaction);

      append_data_binding_request(request,
                                  resolved,
                                  data_id_for_term_role(
                                    resolved,
                                    semantic::v1::ResidualTermKind::volume_source,
                                    semantic::v1::DataRole::forcing),
                                  ResolvedBindingPort::forcing,
                                  "volume_quadrature",
                                  "dealii::Function<dim>",
                                  true);
      append_data_binding_request(request,
                                  resolved,
                                  data_id_for_loss_role(
                                    resolved, semantic::v1::DataRole::desired_state),
                                  ResolvedBindingPort::desired_state,
                                  "observation_quadrature",
                                  "dealii::Function<dim>",
                                  true);

      request.requires_fixed_dirichlet_data = std::any_of(
        specification.transformations.begin(),
        specification.transformations.end(),
        [](const semantic::v1::TransformationSpec &transformation) {
          return !transformation.fixed_data_id.empty();
        });
      std::string fixed_data_id;
      for (const auto &transformation : specification.transformations)
        if (!transformation.fixed_data_id.empty())
          {
            fixed_data_id = transformation.fixed_data_id;
            break;
          }
      append_data_binding_request(
        request,
        resolved,
        request.requires_fixed_dirichlet_data
          ? fixed_data_id
          : std::string{},
        ResolvedBindingPort::fixed_dirichlet_data,
        "dirichlet_boundary_dof_interpolation",
        "dealii::Function<dim>",
        request.requires_fixed_dirichlet_data);

      request.requires_observation_weight = std::any_of(
        specification.observations.begin(),
        specification.observations.end(),
        [&resolved](const semantic::v1::ObservationSpec &observation) {
          return std::any_of(
            observation.data_ids.begin(),
            observation.data_ids.end(),
            [&resolved](const std::string &data_id) {
              return resolved.datum(data_id).role ==
                     semantic::v1::DataRole::observation_weight;
            });
        });
      for (const auto &observation : specification.observations)
        for (const auto &data_id : observation.data_ids)
          if (resolved.datum(data_id).role ==
              semantic::v1::DataRole::observation_weight)
            append_data_binding_request(request,
                                        resolved,
                                        data_id,
                                        ResolvedBindingPort::observation_weight,
                                        "boundary_face_quadrature",
                                        "dealii::Function<dim>",
                                        request.requires_observation_weight);

      request.requires_general_scalar_data = std::any_of(
        specification.residual_terms.begin(),
        specification.residual_terms.end(),
        [](const semantic::v1::ResidualTermSpec &term) {
          return term.kind == semantic::v1::ResidualTermKind::tensor_diffusion;
        });
      for (const auto &term : specification.residual_terms)
        for (const auto &data_id : term.data_ids)
          if (request.requires_general_scalar_data &&
              (resolved.datum(data_id).role == semantic::v1::DataRole::diffusion ||
               resolved.datum(data_id).role ==
                 semantic::v1::DataRole::conservative_transport ||
               resolved.datum(data_id).role ==
                 semantic::v1::DataRole::advective_transport ||
               resolved.datum(data_id).role == semantic::v1::DataRole::reaction ||
               resolved.datum(data_id).role ==
                 semantic::v1::DataRole::robin_coefficient ||
               resolved.datum(data_id).role ==
                 semantic::v1::DataRole::robin_source))
            append_data_binding_request(
              request,
              resolved,
              data_id,
              ResolvedBindingPort::general_scalar_data,
              resolved.datum(data_id).role ==
                    semantic::v1::DataRole::robin_coefficient ||
                  resolved.datum(data_id).role ==
                    semantic::v1::DataRole::robin_source
                ? "boundary_face_quadrature"
                : "volume_quadrature",
              resolved.datum(data_id).kind == semantic::v1::DataKind::tensor_function
                ? "dealii::TensorFunction<2, dim>"
                : resolved.datum(data_id).kind ==
                    semantic::v1::DataKind::vector_function
                  ? "dealii::TensorFunction<1, dim>"
                  : "dealii::Function<dim>",
              request.requires_general_scalar_data);

      request.requires_conservative_transport_data = std::any_of(
        specification.residual_terms.begin(),
        specification.residual_terms.end(),
        [](const semantic::v1::ResidualTermSpec &term) {
          return term.kind == semantic::v1::ResidualTermKind::
                                 conservative_transport;
        }) && !request.requires_general_scalar_data;
      append_data_binding_request(
        request,
        resolved,
        data_id_for_term_role(resolved,
                              semantic::v1::ResidualTermKind::conservative_transport,
                              semantic::v1::DataRole::conservative_transport),
        ResolvedBindingPort::conservative_transport_data,
        "volume_quadrature",
        "dealii::TensorFunction<1, dim>",
        request.requires_conservative_transport_data);

      const auto weighted_observation = std::find_if(
        specification.observations.begin(), specification.observations.end(),
        [](const semantic::v1::ObservationSpec &observation) {
          return observation.kind ==
                 semantic::v1::ObservationKind::weighted_boundary_trace;
        });
      const auto normal_flux_observation = std::find_if(
        specification.observations.begin(), specification.observations.end(),
        [](const semantic::v1::ObservationSpec &observation) {
          return observation.kind == semantic::v1::ObservationKind::normal_flux;
        });
      const auto point_sensor_observation = std::find_if(
        specification.observations.begin(), specification.observations.end(),
        [](const semantic::v1::ObservationSpec &observation) {
          return observation.kind == semantic::v1::ObservationKind::point_sensor;
        });
      const auto has_policy = [&specification](
                                const std::string &subject_id,
                                const semantic::v1::RequirementKind kind,
                                const semantic::v1::RequirementStatus status,
                                const semantic::v1::RequirementScope scope) {
        return std::any_of(
          specification.requirement_policies.begin(),
          specification.requirement_policies.end(),
          [&subject_id, kind, status, scope](
            const semantic::v1::RequirementPolicySpec &policy) {
            return policy.subject_id == subject_id && policy.kind == kind &&
                   policy.status == status && policy.scope == scope;
          });
      };
      if (normal_flux_observation != specification.observations.end())
        {
          request.has_normal_flux_orientation_policy = has_policy(
            normal_flux_observation->id,
            semantic::v1::RequirementKind::conormal_flux,
            semantic::v1::RequirementStatus::selected_discrete_realisation,
            semantic::v1::RequirementScope::both);
          request.has_normal_flux_evaluation_policy = has_policy(
            normal_flux_observation->id,
            semantic::v1::RequirementKind::analytic_quadrature_evaluation,
            semantic::v1::RequirementStatus::selected_discrete_realisation,
            semantic::v1::RequirementScope::discrete_compilation);
        }
      if (point_sensor_observation != specification.observations.end())
        request.has_point_sensor_evaluation_policy = has_policy(
          point_sensor_observation->id,
          semantic::v1::RequirementKind::analytic_quadrature_evaluation,
          semantic::v1::RequirementStatus::selected_discrete_realisation,
          semantic::v1::RequirementScope::discrete_compilation);
      for (const auto &policy : specification.requirement_policies)
        {
          if (policy.kind == semantic::v1::RequirementKind::boundary_trace &&
              weighted_observation != specification.observations.end() &&
              policy.subject_id == weighted_observation->id &&
              policy.typed_trace_selection &&
              policy.typed_trace_selection->weight_data_id != "")
            request.weighted_trace_selection =
              policy.typed_trace_selection;
          if (policy.kind == semantic::v1::RequirementKind::metric_realisation &&
              policy.subject_id == specification.formulation.metric_id &&
              policy.typed_metric_selection)
            request.hminus1_metric_selection = policy.typed_metric_selection;
          if (policy.kind ==
                semantic::v1::RequirementKind::transposition_formulation &&
              policy.typed_transposition_selection)
            request.transposition_selection =
              policy.typed_transposition_selection;
          if (policy.kind == semantic::v1::RequirementKind::boundary_partition &&
              policy.typed_partial_boundary_selection)
            request.partial_boundary_selection =
              policy.typed_partial_boundary_selection;
          if (policy.kind ==
                semantic::v1::RequirementKind::fractional_trace_realisation &&
              policy.typed_fractional_metric_selection)
            request.fractional_metric_selection =
              policy.typed_fractional_metric_selection;
          if (policy.kind ==
                semantic::v1::RequirementKind::tangential_gradient_realisation &&
              policy.typed_boundary_h1_metric_selection)
            request.boundary_h1_metric_selection =
              policy.typed_boundary_h1_metric_selection;
          if (policy.kind == semantic::v1::RequirementKind::target_data_membership &&
              policy.subject_id == "desired_state" &&
              policy.typed_h1_target_data_membership_selection)
            request.h1_target_data_membership_selection =
              policy.typed_h1_target_data_membership_selection;
          if (policy.subject_id == specification.formulation.control_variable_id &&
              policy.kind == semantic::v1::RequirementKind::fixed_dirichlet &&
              policy.status ==
                semantic::v1::RequirementStatus::selected_discrete_realisation)
            request.continuous_control_boundary_region_id = policy.region_id;
        }
      return request;
    }

    static const ResolvedDataBindingRequest *
    find_data_binding_request(const ResolvedCompilationRequest &request,
                              const semantic::v1::DataRole          role,
                              const ResolvedBindingPort              port)
    {
      const auto match = std::find_if(
        request.data_bindings.begin(),
        request.data_bindings.end(),
        [role, port](const ResolvedDataBindingRequest &candidate) {
          return candidate.role == role && candidate.port == port;
        });
      return match == request.data_bindings.end() ? nullptr : &*match;
    }

    template <int dim>
    static void
    validate_resolved_function_bindings(
      const ResolvedCompilationRequest &request,
      const DealiiDataBindings<dim> &    data,
      semantic::v1::ValidationReport &   report)
    {
      using semantic::v1::DataRole;
      using semantic::v1::DiagnosticCategory;

      const auto require_provenance = [&report](
                                        const ResolvedDataBindingRequest &binding,
                                        const std::string &provenance,
                                        const std::string &capability) {
        if (provenance.empty())
          report.add(DiagnosticCategory::lowerability,
                    binding.semantic_id,
                    capability,
                    "Supply a stable provenance label for the resolved Function binding.");
      };

      const auto require_scalar = [&report](
                                    const ResolvedDataBindingRequest &binding,
                                    const unsigned int                n_components,
                                    const std::string &capability) {
        if (n_components != 1)
          report.add(DiagnosticCategory::lowerability,
                     binding.semantic_id,
                     capability,
                     "Bind a one-component scalar Function for the resolved data port.");
      };

      if (const auto *binding = find_data_binding_request(
            request, DataRole::forcing, ResolvedBindingPort::forcing))
        {
          require_scalar(*binding,
                         data.forcing.n_components,
                         "scalar_function_binding_shape");
          require_provenance(*binding,
                             data.provenance.forcing,
                             "forcing_binding_provenance");
        }
      if (const auto *binding = find_data_binding_request(
            request, DataRole::desired_state, ResolvedBindingPort::desired_state))
        {
          require_scalar(*binding,
                         data.desired_state.n_components,
                         "scalar_function_binding_shape");
          require_provenance(*binding,
                             data.provenance.desired_state,
                             "desired_state_binding_provenance");
        }

      if (const auto *binding = find_data_binding_request(
            request,
            DataRole::fixed_dirichlet_lifting,
            ResolvedBindingPort::fixed_dirichlet_data);
          binding != nullptr && request.requires_fixed_dirichlet_data &&
          data.fixed_dirichlet_data)
        require_scalar(*binding,
                       data.fixed_dirichlet_data->get().n_components,
                       "scalar_function_binding_shape");
      if (const auto *binding = find_data_binding_request(
            request,
            DataRole::fixed_dirichlet_lifting,
            ResolvedBindingPort::fixed_dirichlet_data);
          binding != nullptr && request.requires_fixed_dirichlet_data &&
          data.fixed_dirichlet_data)
        require_provenance(*binding,
                           data.provenance.fixed_dirichlet_data,
                           "fixed_dirichlet_binding_provenance");

      if (const auto *binding = find_data_binding_request(
            request,
            DataRole::observation_weight,
            ResolvedBindingPort::observation_weight);
          binding != nullptr && request.requires_observation_weight &&
          data.weighted_trace)
        require_scalar(*binding,
                       data.weighted_trace->weight.n_components,
                       "scalar_boundary_weight_binding");
      if (const auto *binding = find_data_binding_request(
            request,
            DataRole::observation_weight,
            ResolvedBindingPort::observation_weight);
          binding != nullptr && request.requires_observation_weight &&
          data.weighted_trace)
        require_provenance(*binding,
                           data.weighted_trace->provenance,
                           "boundary_weight_binding_provenance");

      if (!request.requires_general_scalar_data || !data.general_scalar)
        return;
      for (const auto role : {DataRole::reaction,
                              DataRole::robin_coefficient,
                              DataRole::robin_source})
        {
          const auto *binding = find_data_binding_request(
            request, role, ResolvedBindingPort::general_scalar_data);
          if (binding == nullptr)
            continue;
          const auto *function = role == DataRole::reaction
                                   ? &data.general_scalar->reaction
                                   : role == DataRole::robin_coefficient
                                       ? &data.general_scalar->robin_coefficient
                                       : &data.general_scalar->robin_source;
          require_scalar(*binding,
                         function->n_components,
                         "scalar_coefficient_function_shape");
          const auto &provenance = role == DataRole::reaction
                                     ? data.general_scalar->provenance.reaction
                                     : role == DataRole::robin_coefficient
                                         ? data.general_scalar->provenance.robin_coefficient
                                         : data.general_scalar->provenance.robin_source;
          require_provenance(*binding,
                             provenance,
                             "coefficient_binding_provenance");
        }

      if (request.requires_conservative_transport_data &&
          data.conservative_transport)
        if (const auto *binding = find_data_binding_request(
              request,
              DataRole::conservative_transport,
              ResolvedBindingPort::conservative_transport_data))
          require_provenance(
            *binding,
            data.conservative_transport->provenance.conservative_transport,
            "conservative_transport_binding_provenance");
    }

    template <int dim>
    CompilationResultT<dealii_backend::SerialBackend>
    compile_impl(
            const semantic::v1::ProblemSpec &  specification,
            dealii::Triangulation<dim> &        triangulation,
            const DealiiDataBindings<dim> &     data,
            const DealiiDiscretisationPolicy &  policy,
            std::optional<CellwiseBoxDataBindings> bounds,
            std::optional<FacewiseBoxDataBindings> facewise_bounds,
            std::shared_ptr<const void>             lifetime_owner,
            std::string                             mesh_provenance,
            const bool                              owns_mesh,
            const CompilationProduct               product) const
    {
      using Backend = dealii_backend::SerialBackend;
      CompilationResultT<Backend> result;
      auto resolution = semantic::v1::SemanticResolver().resolve(specification);
      result.diagnostics = std::move(resolution.diagnostics);
      if (!result.diagnostics.valid())
        return result;
      auto request = resolve_compilation_request(*resolution.problem);
      const auto dirichlet_registration =
        resolve_dirichlet_control_registration(*resolution.problem, request);
      close_compilation_request(*resolution.problem,
                                request,
                                dirichlet_registration);
      validate_lowerability(specification, request, policy, result.diagnostics);
      validate_formulation_capability(specification, result.diagnostics);
      validate_supplied_otd_capability(specification,
                                       request,
                                       product,
                                       result.diagnostics);
      validate_dirichlet_control_registration(*resolution.problem,
                                              request,
                                              result.diagnostics);
      validate_product_capability(specification,
                                  request,
                                  policy,
                                  product,
                                  result.diagnostics);
      if (!result.diagnostics.valid())
        return result;
      const bool uses_fixed_reconstruction =
        request.uses_fixed_reconstruction;
      const bool uses_dirichlet_control = request.uses_dirichlet_control;
      const bool uses_l2_dirichlet_control =
        request.uses_l2_dirichlet_control;
      const bool uses_normalized_laplacian = request.uses_normalized_laplacian;
      const bool uses_partial_dirichlet_control =
        request.uses_partial_dirichlet_control;
      const bool uses_neumann_boundary_control =
        request.uses_neumann_boundary_control;
      const bool uses_neumann_convection = request.uses_neumann_convection;
      const bool uses_mean_zero_gauge = request.uses_mean_zero_gauge;
      const bool uses_h1_control_regularisation =
        request.uses_h1_control_regularisation;
      const bool uses_h1_control_metric = request.uses_h1_control_metric;
      const bool uses_hhalf_control_metric = request.uses_hhalf_control_metric;
      const bool uses_hminus1_control_metric =
        request.uses_hminus1_control_metric;
      const bool uses_homogeneous_dirichlet_continuous_control =
        request.uses_homogeneous_dirichlet_continuous_control;
      const bool uses_coefficient_identification =
        request.uses_coefficient_identification;
      const bool uses_general_scalar = request.uses_general_scalar;
      const bool uses_h1_state_observation = request.uses_h1_state_observation;
      const bool uses_weighted_boundary_trace =
        request.uses_weighted_boundary_trace;
      const bool uses_point_sensor = request.uses_point_sensor;
      const bool uses_normal_flux = request.uses_normal_flux;
      const bool uses_supplied_otd =
        specification.formulation.kind ==
          semantic::v1::FormulationKind::all_at_once &&
        specification.formulation.provenance ==
          semantic::v1::FormulationProvenance::supplied_otd;
      const bool uses_h1_dirichlet_control = request.uses_h1_dirichlet_control;
      const bool uses_hhalf_dirichlet_registration =
        request.uses_hhalf_dirichlet_registration;
      const bool uses_h1_tracking_hhalf_dirichlet_registration =
        request.uses_h1_tracking_hhalf_dirichlet_registration;
      if (uses_h1_state_observation)
        {
          if (!request.h1_target_data_membership_selection)
            result.diagnostics.add(
              semantic::v1::DiagnosticCategory::lowerability,
              "desired_state",
              "h1_target_space_membership",
              "Resolve the model-author H1 target-space membership assumption before constructing the registered model.");
          else
            {
              const auto &selection =
                *request.h1_target_data_membership_selection;
              const std::string &expected_state_boundary =
                request.uses_dirichlet_control
                  ? request.control_boundary_region_id
                  : request.fixed_boundary_region_id;
              if (selection.data_id != "desired_state" ||
                  selection.observation_space_id != "state_observation_space" ||
                  selection.regularity_realisation !=
                    semantic::v1::H1TargetDataRegularityRealisation::
                      h1_value_and_weak_gradient)
                result.diagnostics.add(
                  semantic::v1::DiagnosticCategory::lowerability,
                  "desired_state",
                  "h1_target_space_membership",
                  "The registered H1-state lowerer requires the declared desired_state H1 value/gradient membership selection.");
              if (selection.fixed_boundary_region_id != expected_state_boundary ||
                  selection.trace_realisation !=
                    semantic::v1::H1TargetDataTraceRealisation::
                      zero_trace_on_fixed_boundary)
                result.diagnostics.add(
                  semantic::v1::DiagnosticCategory::lowerability,
                  "desired_state",
                  "h1_target_zero_trace",
                  "The registered H1-state lowerer requires the declared zero trace on the fixed Dirichlet boundary.");
            }
        }
      const bool uses_transposition_policy =
        uses_point_sensor || uses_normal_flux || uses_l2_dirichlet_control;
      if (uses_transposition_policy)
        {
          if (!request.transposition_selection)
            result.diagnostics.add(
              semantic::v1::DiagnosticCategory::lowerability,
              specification.formulation.equation_id,
              "transposition_realisation_request",
              "Resolve the typed transposition selection before constructing the registered model.");
          else
            {
              const auto expected_realisation =
                uses_point_sensor
                  ? semantic::v1::TranspositionDiscreteRealisation::
                      fe_q_point_sensor_very_weak
                : uses_normal_flux
                  ? semantic::v1::TranspositionDiscreteRealisation::
                      fe_q_normal_flux_very_weak
                  : semantic::v1::TranspositionDiscreteRealisation::
                      conforming_nodal_lifting_equivalence;
              if (request.transposition_selection->operator_realisation !=
                    semantic::v1::TranspositionOperatorRealisation::
                      scalar_diffusion_reaction_dirichlet_laplacian)
                result.diagnostics.add(
                  semantic::v1::DiagnosticCategory::lowerability,
                  specification.formulation.equation_id,
                  "transposition_operator_realisation",
                  "The registered lowerer supports only the scalar diffusion-reaction Dirichlet-Laplacian transposition.");
              if (request.transposition_selection->isomorphism_id !=
                    "dirichlet_laplacian_isomorphism" ||
                  request.transposition_selection->subject_equation_id !=
                    specification.formulation.equation_id)
                result.diagnostics.add(
                  semantic::v1::DiagnosticCategory::lowerability,
                  specification.formulation.equation_id,
                  "transposition_isomorphism",
                  "The typed transposition operator must match the selected state equation.");
              if (request.transposition_selection->diffusion_data_id !=
                  request.transposition_diffusion_data_id)
                result.diagnostics.add(
                  semantic::v1::DiagnosticCategory::lowerability,
                  specification.formulation.equation_id,
                  "transposition_diffusion_data",
                  "The typed transposition operator must bind the diffusion data port selected by the residual.");
              if (request.transposition_selection->reaction_data_id !=
                  request.transposition_reaction_data_id)
                result.diagnostics.add(
                  semantic::v1::DiagnosticCategory::lowerability,
                  specification.formulation.equation_id,
                  "transposition_reaction_data",
                  "The typed transposition operator must bind the reaction data port selected by the residual.");
              if ((uses_point_sensor || uses_normal_flux) &&
                  request.transposition_selection->observation_id !=
                    "state_observation")
                result.diagnostics.add(
                  semantic::v1::DiagnosticCategory::lowerability,
                  specification.formulation.equation_id,
                  "transposition_observation_source",
                  "Bind the P5.3 transposition to the selected state observation.");
              if (request.transposition_selection->discrete_realisation !=
                  expected_realisation)
                result.diagnostics.add(
                  semantic::v1::DiagnosticCategory::lowerability,
                  specification.formulation.equation_id,
                  "transposition_discrete_realisation",
                  "The typed transposition realization does not match the selected registered observation or lifting lowerer.");
            }
        }
      if (uses_partial_dirichlet_control)
        {
          if (!request.partial_boundary_selection)
            result.diagnostics.add(
              semantic::v1::DiagnosticCategory::lowerability,
              specification.formulation.state_variable_id,
              "partial_dirichlet_partition_policy",
              "Resolve the typed partial Dirichlet boundary partition before constructing the lifting model.");
          else if (
            request.partial_boundary_selection->interface_realisation !=
                semantic::v1::PartialDirichletInterfaceRealisation::
                  fixed_data_precedence ||
            request.partial_boundary_selection->trace_realisation !=
                semantic::v1::PartialDirichletTraceRealisation::
                  relative_interior_nodal_zero_endpoint ||
            request.partial_boundary_selection->hanging_realisation !=
                semantic::v1::PartialDirichletHangingRealisation::unsupported)
            result.diagnostics.add(
              semantic::v1::DiagnosticCategory::lowerability,
              specification.formulation.state_variable_id,
              "partial_dirichlet_interface_selection",
              "The registered partial lifting supports only fixed-data precedence with zero-endpoint relative-interior traces and no hanging relation.");
        }
      if (uses_hhalf_control_metric)
        {
          if (!request.fractional_metric_selection)
            result.diagnostics.add(
              semantic::v1::DiagnosticCategory::lowerability,
              specification.formulation.metric_id,
              "hhalf_metric_realisation_selection",
              "Resolve the typed fractional trace metric selection before constructing the metric.");
          else if (
            request.fractional_metric_selection->operator_realisation !=
                semantic::v1::FractionalTraceOperatorRealisation::
                  volume_mass_plus_stiffness_schur ||
            request.fractional_metric_selection->apply_realisation !=
                semantic::v1::FractionalTraceApplyRealisation::
                  minimum_h1_extension ||
            request.fractional_metric_selection->inverse_realisation !=
                semantic::v1::FractionalTraceInverseRealisation::
                  full_volume_operator_inverse)
            result.diagnostics.add(
              semantic::v1::DiagnosticCategory::lowerability,
              specification.formulation.metric_id,
              "hhalf_metric_realisation_selection",
              "The registered lowerer supports only the minimum-extension H1/2 Schur-complement metric.");
        }
      if (uses_h1_dirichlet_control)
        {
          if (!request.boundary_h1_metric_selection)
            result.diagnostics.add(
              semantic::v1::DiagnosticCategory::lowerability,
              specification.formulation.metric_id,
              "boundary_h1_metric_realisation_selection",
              "Resolve the typed boundary H1 metric selection before constructing the metric.");
          else if (
            request.boundary_h1_metric_selection->operator_realisation !=
                semantic::v1::BoundaryH1MetricOperatorRealisation::
                  boundary_mass_plus_tangential_stiffness ||
            request.boundary_h1_metric_selection->tangential_gradient_realisation !=
                semantic::v1::BoundaryH1TangentialGradientRealisation::
                  projected_ambient_gradient ||
            request.boundary_h1_metric_selection->nullspace_realisation !=
                semantic::v1::BoundaryH1MetricNullspaceRealisation::
                  positive_mass_no_nullspace)
            result.diagnostics.add(
              semantic::v1::DiagnosticCategory::lowerability,
              specification.formulation.metric_id,
              "boundary_h1_metric_realisation_selection",
              "The registered lowerer supports only the projected-gradient boundary H1 metric with positive mass.");
        }
      const auto *tracking_region =
        request.tracking_region_id.empty()
          ? nullptr
          : find_region(specification, request.tracking_region_id);
      const auto *robin_boundary_region =
        request.robin_boundary_region_id.empty()
          ? nullptr
          : find_region(specification, request.robin_boundary_region_id);
      const auto *partial_fixed_boundary_region =
        request.partial_fixed_boundary_region_id.empty()
          ? nullptr
          : find_region(specification,
                        request.partial_fixed_boundary_region_id);
      const auto *partial_control_boundary_region =
        request.partial_control_boundary_region_id.empty()
          ? nullptr
          : find_region(specification,
                        request.partial_control_boundary_region_id);
      const auto *control_boundary_region =
        request.control_boundary_region_id.empty()
          ? nullptr
          : find_region(specification, request.control_boundary_region_id);
      const auto *fixed_boundary_region =
        request.fixed_boundary_region_id.empty()
          ? nullptr
          : find_region(specification, request.fixed_boundary_region_id);
      const auto *continuous_control_boundary_region =
        request.continuous_control_boundary_region_id.empty()
          ? nullptr
          : find_region(specification,
                        request.continuous_control_boundary_region_id);
      const bool uses_assembled_v1_target =
        request.uses_assembled_v1_target;
      std::optional<ScalarLoweringPlan> scalar_plan;
      if (uses_assembled_v1_target)
        {
          auto planned = scalar_planner_.plan(*resolution.problem);
          for (const auto &diagnostic : planned.diagnostics.diagnostics())
            result.diagnostics.add(diagnostic.category,
                                   diagnostic.component_id,
                                   diagnostic.capability,
                                   diagnostic.remedy);
          scalar_plan = std::move(planned.plan);
        }
      validate_resolved_function_bindings(request, data, result.diagnostics);
      if (triangulation.n_active_cells() == 0)
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.id,
          "nonempty_triangulation",
          "Compile on a triangulation with at least one active cell.");
      if (uses_point_sensor && tracking_region != nullptr)
        validate_point_sensor_mesh(triangulation, *tracking_region,
                                   result.diagnostics);
      if (request.requires_observation_weight && !data.weighted_trace)
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          "boundary_weight",
          "weighted_boundary_trace_data_binding",
          "Bind the fixed scalar Function consumed by the weighted boundary trace.");
      if (!request.requires_observation_weight && data.weighted_trace)
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.id,
          "selected_weighted_boundary_trace",
          "Declare a weighted boundary trace before binding observation-weight data.");
      if (uses_weighted_boundary_trace &&
          !request.weighted_trace_selection.has_value())
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.id,
          "weighted_trace_realisation_request",
          "Resolve the typed weighted boundary-trace realization before lowering.");
      if (uses_hminus1_control_metric &&
          !request.hminus1_metric_selection.has_value())
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.metric_id,
          "hminus1_metric_realisation_request",
          "Resolve the typed H-1 metric realization before lowering.");
      if (!uses_normalized_laplacian && !uses_coefficient_identification &&
          !uses_general_scalar &&
          !data.diffusion)
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          "diffusion",
          "diffusion_data_binding",
          "Bind the constant diffusion coefficient selected by this graph.");
      if (!uses_normalized_laplacian && !uses_coefficient_identification &&
          !uses_general_scalar &&
          data.diffusion &&
          (!std::isfinite(*data.diffusion) || *data.diffusion <= 0.0))
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          "diffusion",
          "positive_finite_diffusion_binding",
          "Bind a positive finite constant diffusion coefficient.");
      if (!uses_normalized_laplacian && !uses_general_scalar &&
          (!std::isfinite(data.reaction) || data.reaction < 0.0))
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          "reaction",
          "nonnegative_finite_reaction_binding",
          "Bind a nonnegative finite reaction coefficient.");
      if (request.requires_general_scalar_data && !data.general_scalar)
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.id,
          "general_scalar_coefficient_bindings",
          "Bind tensor diffusion, both vector transports, scalar reaction, Robin coefficient, and Robin source Functions.");
      if (!request.requires_general_scalar_data && data.general_scalar)
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.id,
          "selected_general_scalar_target",
          "Remove general scalar coefficient bindings unless the graph selects the P5.1 target.");
      if (request.requires_conservative_transport_data &&
          !data.conservative_transport)
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          "conservative_transport",
          "conservative_transport_data_binding",
          "Bind the conservative transport Function selected by the C5.6 Neumann-control composition.");
      if (!request.requires_conservative_transport_data &&
          data.conservative_transport)
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.id,
          "selected_neumann_convection_target",
          "Remove conservative transport data unless the graph declares the registered C5.6 composition.");
      if (request.requires_conservative_transport_data &&
          data.conservative_transport &&
          data.conservative_transport->provenance.conservative_transport.empty())
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          "conservative_transport",
          "conservative_transport_binding_provenance",
          "Supply a stable provenance label for the conservative transport Function binding.");
      if (!std::isfinite(data.regularisation_weight) ||
          data.regularisation_weight <= 0.0)
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          "regularisation",
          "positive_finite_regularisation_binding",
          "Bind a positive finite regularisation weight.");
      if (!valid_metric_solve_policy(policy.control_metric_solve))
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.metric_id,
          "valid_metric_solve_policy",
          "Select positive finite metric-solve tolerances and a positive iteration limit.");
      if (!dealii_backend::valid(policy.state_solve))
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.id,
          "valid_state_solve_policy",
          "Select positive finite state-solve tolerances.");
      if (!dealii_backend::valid(policy.adjoint_solve))
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.id,
          "valid_adjoint_solve_policy",
          "Select positive finite adjoint-solve tolerances.");
      if (request.requires_fixed_dirichlet_data &&
          !data.fixed_dirichlet_data)
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.state_variable_id,
          "fixed_dirichlet_data_binding",
          "Bind fixed Dirichlet Function data for the declared reconstruction.");
      if (!request.requires_fixed_dirichlet_data &&
          !uses_fixed_reconstruction && !uses_partial_dirichlet_control &&
          data.fixed_dirichlet_data)
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.state_variable_id,
          "selected_fixed_dirichlet_reconstruction",
          "Declare the fixed-Dirichlet reconstruction before binding lifting data.");
      if (uses_mean_zero_gauge && std::abs(data.reaction) > 1e-14)
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.state_variable_id,
          "pure_neumann_zero_reaction",
          "Bind zero reaction for the selected pure-Neumann constant-nullspace policy.");
      if (uses_mean_zero_gauge && uses_fixed_reconstruction)
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.state_variable_id,
          "pure_neumann_without_fixed_reconstruction",
          "Remove fixed-Dirichlet reconstruction when selecting the pure-Neumann mean constraint.");
      const bool has_constraint = !specification.formulation.constraint_id.empty();
      if (has_constraint &&
          ((uses_neumann_boundary_control && !facewise_bounds) ||
           (!uses_neumann_boundary_control && !bounds)))
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.constraint_id,
          "bound_data_binding",
          uses_neumann_boundary_control
            ? "Bind scalar constants or exact facewise boundary-control vectors for both bounds."
            : "Bind scalar constants or FE_DGQ(0) coefficient vectors for both bounds.");
      if (bounds && !valid_bound_representation(*bounds))
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.constraint_id,
          "bound_data_representation",
          "Bind both cellwise bounds as scalars or both as FE_DGQ(0) vectors.");
      if (facewise_bounds && !valid_facewise_bound_representation(*facewise_bounds))
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.constraint_id,
          "facewise_bound_data_representation",
          "Bind both facewise bounds as scalars or both as exact boundary-control vectors.");
      if (uses_coefficient_identification && bounds &&
          !has_strictly_positive_lower_bound(*bounds))
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.constraint_id,
          "positive_parameter_lower_bound",
          "Bind a strictly positive scalar or every strictly positive cellwise lower parameter bound.");
      if (uses_neumann_boundary_control && bounds)
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.constraint_id,
          "cellwise_bounds_for_boundary_control",
          "Bind the declared facewise box data for the boundary control.");
      if (!uses_neumann_boundary_control && facewise_bounds)
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.constraint_id,
          "facewise_bounds_for_volume_control",
          "Bind the declared cellwise box data for the volume control.");
      if (!has_constraint && (bounds || facewise_bounds))
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.id,
          "unselected_bound_data",
          "Remove bound data when the reduced formulation selects no constraint.");
      if (bounds && valid_bound_representation(*bounds))
        validate_cellwise_bound_values(*bounds,
                                       triangulation.n_active_cells(),
                                       specification.formulation.constraint_id,
                                       result.diagnostics);
      if (facewise_bounds && valid_facewise_bound_representation(*facewise_bounds) &&
          control_boundary_region != nullptr)
        validate_facewise_bound_values(
          *facewise_bounds,
          count_boundary_faces(triangulation,
                               boundary_ids(*control_boundary_region)),
          specification.formulation.constraint_id,
          result.diagnostics);
      if (!result.diagnostics.valid())
        return result;

      // Resolve the target and all semantic component choices before any
      // backend model is constructed. The executable consumes the same
      // closed request and scalar plan captured by this decision.
      const CompiledTargetKind target_kind =
        target_kind_from_request(request);
      const ResolvedCompilationDecision resolved_decision =
        make_resolved_decision(specification,
                               policy,
                               target_kind,
                               request,
                               scalar_plan ? &*scalar_plan : nullptr,
                               data,
                               bounds,
                               facewise_bounds,
                               triangulation,
                               mesh_provenance,
                               owns_mesh);

      const auto fixed_dirichlet_boundary_ids =
        uses_general_scalar && scalar_plan
          ? boundary_ids_from_plan(scalar_plan->dirichlet_boundary_ids)
        : uses_partial_dirichlet_control &&
            partial_fixed_boundary_region != nullptr
          ? boundary_ids(*partial_fixed_boundary_region)
        : uses_mean_zero_gauge ||
            (uses_dirichlet_control && !uses_partial_dirichlet_control)
          ? std::set<dealii::types::boundary_id>{}
          : fixed_boundary_region == nullptr
            ? std::set<dealii::types::boundary_id>{}
            : boundary_ids(*fixed_boundary_region);
      if ((uses_point_sensor || uses_normal_flux) &&
          !controls_complete_exterior_boundary(
            triangulation, fixed_dirichlet_boundary_ids))
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.state_variable_id,
          "p53_complete_fixed_dirichlet_boundary",
          "Select a fixed-Dirichlet boundary region covering every exterior face for the registered P5.3 target.");
      const auto dirichlet_boundary_ids = uses_mean_zero_gauge
                                            ? std::set<dealii::types::boundary_id>{}
                                            : uses_dirichlet_control
                                              ? (uses_partial_dirichlet_control &&
                                                 partial_control_boundary_region != nullptr
                                                   ? boundary_ids(
                                                       *partial_control_boundary_region)
                                                   : control_boundary_region == nullptr
                                                     ? std::set<dealii::types::boundary_id>{}
                                                     : boundary_ids(*control_boundary_region))
                                            : fixed_dirichlet_boundary_ids;
      const auto continuous_control_boundary_ids =
        uses_homogeneous_dirichlet_continuous_control &&
            continuous_control_boundary_region != nullptr
          ? boundary_ids(*continuous_control_boundary_region)
          : std::set<dealii::types::boundary_id>{};
      if (triangulation.has_hanging_nodes())
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.id,
          "conforming_mesh_without_hanging_nodes",
          "Compile the registered serial targets on a mesh without hanging-node relations.");
      if (!uses_mean_zero_gauge && !uses_dirichlet_control &&
          !contains_all_boundary_ids(triangulation, dirichlet_boundary_ids))
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.state_variable_id,
          "fixed_dirichlet_boundary_presence",
          "Select fixed-Dirichlet boundary ids present on the compiled mesh.");
      if (uses_general_scalar &&
          (robin_boundary_region == nullptr ||
           !contains_all_boundary_ids(triangulation,
                                      uses_general_scalar && scalar_plan
                                        ? boundary_ids_from_plan(
                                            scalar_plan->robin_boundary_ids)
                                        : boundary_ids(*robin_boundary_region))))
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          robin_boundary_region == nullptr ? specification.id :
                                             robin_boundary_region->id,
          "robin_boundary_presence",
          "Select Robin boundary ids present on the compiled mesh.");
      if (uses_general_scalar && robin_boundary_region != nullptr &&
          !forms_complete_boundary_partition(triangulation,
                                             dirichlet_boundary_ids,
                                             uses_general_scalar && scalar_plan
                                               ? boundary_ids_from_plan(
                                                   scalar_plan->robin_boundary_ids)
                                               : boundary_ids(*robin_boundary_region)))
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.state_variable_id,
          "complete_scalar_boundary_partition",
          "Partition every exterior boundary face into the declared fixed-Dirichlet or Robin region.");
      if (uses_normal_flux && tracking_region != nullptr &&
          (tracking_region->kind != semantic::v1::RegionKind::boundary ||
           !contains_all_boundary_ids(triangulation,
                                      boundary_ids(*tracking_region))))
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          tracking_region == nullptr ? specification.id : tracking_region->id,
          "normal_flux_boundary_presence",
          "Select normal-flux boundary ids present on the compiled mesh.");
      if (uses_neumann_boundary_control && control_boundary_region != nullptr &&
          !contains_all_boundary_ids(triangulation,
                                     boundary_ids(*control_boundary_region)))
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          control_boundary_region->id,
          "neumann_control_boundary_presence",
          "Select Neumann-control boundary ids present on the compiled mesh.");
      if (uses_neumann_convection && control_boundary_region != nullptr &&
          !forms_complete_boundary_partition(
            triangulation,
            dirichlet_boundary_ids,
            boundary_ids(*control_boundary_region)))
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.state_variable_id,
          "complete_neumann_convection_boundary_partition",
          "Partition every exterior boundary face into the declared fixed-Dirichlet or Neumann-control region.");
      if (uses_neumann_boundary_control && !uses_neumann_convection &&
          tracking_region != nullptr &&
          !contains_all_boundary_ids(triangulation,
                                     boundary_ids(*tracking_region)))
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          tracking_region->id,
          "boundary_observation_presence",
          "Select boundary-observation ids present on the compiled mesh.");
      if (uses_neumann_convection && tracking_region != nullptr &&
          !contains_all_material_ids(triangulation,
                                     selected_tracking_material_ids(*tracking_region)))
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          tracking_region->id,
          "subdomain_observation_material_presence",
          "Select material ids present on the compiled mesh for C5.6 subdomain tracking.");
      if (uses_homogeneous_dirichlet_continuous_control &&
          (continuous_control_boundary_region == nullptr ||
           !controls_complete_exterior_boundary(
             triangulation, continuous_control_boundary_ids)))
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.control_variable_id,
          "continuous_control_complete_boundary",
          "Select a homogeneous continuous-control boundary covering every exterior face.");
      if (uses_homogeneous_dirichlet_continuous_control &&
          !controls_complete_exterior_boundary(triangulation,
                                               dirichlet_boundary_ids))
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.state_variable_id,
          "continuous_state_complete_boundary",
          "Select a fixed state boundary covering every exterior face for the H1_0 control comparison.");
      if (uses_homogeneous_dirichlet_continuous_control &&
          continuous_control_boundary_ids != dirichlet_boundary_ids)
        result.diagnostics.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.formulation.control_variable_id,
          "continuous_control_boundary_match",
          "Use the same resolved fixed boundary ids for the state and continuous control.");
      if (!result.diagnostics.valid())
        return result;
      if (uses_dirichlet_control &&
          ((!uses_partial_dirichlet_control &&
            !controls_complete_exterior_boundary(triangulation,
                                                 dirichlet_boundary_ids)) ||
           (uses_partial_dirichlet_control &&
            (!contains_all_boundary_ids(triangulation,
                                        fixed_dirichlet_boundary_ids) ||
             !contains_all_boundary_ids(triangulation,
                                        dirichlet_boundary_ids) ||
             !forms_complete_boundary_partition(
               triangulation,
               fixed_dirichlet_boundary_ids,
               dirichlet_boundary_ids)))))
        {
          result.diagnostics.add(
            semantic::v1::DiagnosticCategory::lowerability,
            specification.formulation.state_variable_id,
            uses_partial_dirichlet_control
              ? "complete_partial_dirichlet_boundary_partition"
              : "complete_dirichlet_control_boundary",
            uses_partial_dirichlet_control
              ? "Partition every exterior boundary face into the declared fixed or controlled Dirichlet region."
              : "Select every exterior boundary id for the registered nodal Dirichlet lifting; partial boundaries, interfaces, and undeclared corner policies are not supported.");
          return result;
        }
      contract::require(tracking_region != nullptr,
                        "Validated v1 problem has no tracking observation region");
      std::shared_ptr<const contract::MetricT<Backend>> metric;
      std::shared_ptr<const CompiledCellwiseBoxDataT<Backend>> box_data;
      std::shared_ptr<const contract::ConstraintT<Backend>> constraint;
      std::shared_ptr<const contract::ExecutableModelT<Backend>> executable;
      std::shared_ptr<const contract::ReducedHessianT<Backend>> reduced_hessian;
      std::shared_ptr<const contract::SuppliedOTDSystemT<Backend>>
        supplied_otd_system;
      contract::StateAdjointSolversT<Backend> solvers;
      ConstraintRealisation constraint_realisation = ConstraintRealisation::none;
      if (uses_neumann_boundary_control)
        {
          contract::require(control_boundary_region != nullptr,
                            "Validated v1 problem has no Neumann control region");
          using BoundaryModel = detail::NeumannBoundaryControlModel<dim>;
          const auto boundary = std::make_shared<BoundaryModel>(
            triangulation,
            data.forcing,
            data.desired_state,
            *data.diffusion,
            data.reaction,
            data.regularisation_weight,
            policy.state_degree,
            dirichlet_boundary_ids,
            boundary_ids(*control_boundary_region),
            uses_neumann_convection
              ? std::set<dealii::types::boundary_id>{}
              : boundary_ids(*tracking_region),
            uses_weighted_boundary_trace ? &data.weighted_trace->weight : nullptr,
            uses_mean_zero_gauge
              ? BoundaryModel::StateGauge::mean_zero_multiplier
              : BoundaryModel::StateGauge::fixed_dirichlet,
            uses_neumann_convection
              ? BoundaryModel::StateObservation::volume_restriction
              : BoundaryModel::StateObservation::boundary_trace,
            uses_neumann_convection
              ? selected_tracking_material_ids(*tracking_region)
              : std::set<dealii::types::material_id>{},
            uses_neumann_convection
              ? &data.conservative_transport->conservative_transport
              : nullptr,
            uses_weighted_boundary_trace
              ? std::optional<typename BoundaryModel::WeightedTraceRealisation>{
                  BoundaryModel::WeightedTraceRealisation::fe_q_face_quadrature}
              : std::nullopt);
          if (uses_mean_zero_gauge && !boundary->forcing_is_compatible())
            {
              result.diagnostics.add(
                semantic::v1::DiagnosticCategory::lowerability,
                "forcing",
                "pure_neumann_forcing_compatibility",
                "Bind forcing with zero discrete pairing against the constant null mode.");
              return result;
            }
          metric = std::make_shared<dealii_backend::MassMetric>(
            boundary->control_l2_metric(policy.control_metric_solve));
          if (has_constraint)
            {
              constraint =
                std::make_shared<dealii_backend::FacewiseBoxConstraint>(
                  make_facewise_constraint(*boundary,
                                           *facewise_bounds,
                                           as_mass_metric(*metric)));
              constraint_realisation = ConstraintRealisation::facewise_l2;
            }
          solvers = contract::StateAdjointSolversT<Backend>{
            [boundary, solve_policy = policy.state_solve](
              const contract::PrimalBlockT<Backend> &control) {
              return boundary->solve_state_with_report(control, solve_policy);
            },
            [boundary, solve_policy = policy.adjoint_solve](
              const contract::PrimalBlockT<Backend> &full_point,
              const contract::CovectorBlockT<Backend> &state_rhs) {
              return boundary->solve_adjoint_with_report(full_point,
                                                         state_rhs,
                                                         solve_policy);
            }};
          executable = boundary;
        }
      else if (uses_dirichlet_control)
        {
          using DirichletModel = detail::DirichletControlLiftingModel<dim>;
          detail::DirichletObjectivePolicy objective_policy;
          objective_policy.state_tracking =
            uses_h1_tracking_hhalf_dirichlet_registration
              ? detail::DirichletStateTrackingNormKind::h1
              : detail::DirichletStateTrackingNormKind::l2;
          objective_policy.control_norm =
            request.dirichlet_registration ==
                ResolvedDirichletRegistration::hhalf_control
              ? detail::DirichletControlNormKind::hhalf
            : uses_h1_dirichlet_control
              ? detail::DirichletControlNormKind::h1
              : detail::DirichletControlNormKind::l2;
          objective_policy.search_metric =
            uses_hhalf_dirichlet_registration
              ? detail::DirichletControlNormKind::hhalf
            : uses_h1_dirichlet_control
              ? detail::DirichletControlNormKind::h1
              : detail::DirichletControlNormKind::l2;
          objective_policy.trace_metric_solve = policy.control_metric_solve;
          const auto dirichlet = std::make_shared<DirichletModel>(
            triangulation,
            data.forcing,
            data.desired_state,
            uses_normalized_laplacian ? 1.0 : *data.diffusion,
            uses_normalized_laplacian ? 0.0 : data.reaction,
            data.regularisation_weight,
            policy.state_degree,
            dirichlet_boundary_ids,
            fixed_dirichlet_boundary_ids,
            data.fixed_dirichlet_data,
            objective_policy);
          if (uses_hhalf_dirichlet_registration)
            metric = std::make_shared<dealii_backend::TraceHhalfMetric>(
              dirichlet->control_hhalf_metric());
          else
            metric = std::make_shared<dealii_backend::MassMetric>(
              uses_h1_dirichlet_control
                ? dirichlet->control_h1_metric(policy.control_metric_solve)
                : dirichlet->control_l2_metric(policy.control_metric_solve));
          solvers = contract::StateAdjointSolversT<Backend>{
            [dirichlet, solve_policy = policy.state_solve](
              const contract::PrimalBlockT<Backend> &control) {
              return dirichlet->solve_state_with_report(control, solve_policy);
            },
            [dirichlet, solve_policy = policy.adjoint_solve](
              const contract::PrimalBlockT<Backend> &full_point,
              const contract::CovectorBlockT<Backend> &state_rhs) {
              return dirichlet->solve_adjoint_with_report(full_point,
                                                          state_rhs,
                                                          solve_policy);
            }};
          executable = dirichlet;
        }
      else if (uses_h1_control_regularisation ||
               uses_homogeneous_dirichlet_continuous_control)
        {
          using H1Model = detail::ContinuousControlModel<dim>;
          const auto h1_control = std::make_shared<H1Model>(
            triangulation,
            data.forcing,
            data.desired_state,
            *data.diffusion,
            data.reaction,
            data.regularisation_weight,
            policy.state_degree,
            dirichlet_boundary_ids,
            uses_h1_state_observation,
            uses_h1_control_regularisation,
            uses_homogeneous_dirichlet_continuous_control,
            continuous_control_boundary_ids);
          if (uses_hminus1_control_metric)
            {
              contract::require(request.hminus1_metric_selection.has_value(),
                                "Validated H-1 target has no resolved metric selection");
              contract::require(
                request.hminus1_metric_selection->operator_realisation ==
                  semantic::v1::Hminus1MetricOperatorRealisation::mass_laplacian_inverse_mass &&
                  request.hminus1_metric_selection->inverse_realisation ==
                    semantic::v1::Hminus1MetricInverseRealisation::mass_inverse_laplacian_mass_inverse,
                "Validated H-1 target selected an unsupported metric realization");
              metric = std::make_shared<dealii_backend::Hminus1Metric>(
                h1_control->control_hminus1_metric(
                  policy.control_metric_solve,
                  dealii_backend::Hminus1OperatorRealisation::mass_laplacian_inverse_mass,
                  dealii_backend::Hminus1InverseRealisation::mass_inverse_laplacian_mass_inverse));
            }
          else
            metric = std::make_shared<dealii_backend::MassMetric>(
              uses_h1_control_metric
                ? h1_control->control_h1_metric(policy.control_metric_solve)
                : h1_control->control_l2_metric(policy.control_metric_solve));
          solvers = contract::StateAdjointSolversT<Backend>{
            [h1_control, solve_policy = policy.state_solve](
              const contract::PrimalBlockT<Backend> &control) {
              return h1_control->solve_state_with_report(control, solve_policy);
            },
            [h1_control, solve_policy = policy.adjoint_solve](
              const contract::PrimalBlockT<Backend> &full_point,
              const contract::CovectorBlockT<Backend> &state_rhs) {
              return h1_control->solve_adjoint_with_report(full_point,
                                                           state_rhs,
                                                           solve_policy);
            }};
          executable = h1_control;
        }
      else if (uses_coefficient_identification)
        {
          using CoefficientModel = detail::CoefficientIdentificationModel<dim>;
          const auto coefficient = std::make_shared<CoefficientModel>(
            triangulation,
            data.forcing,
            data.desired_state,
            data.reaction,
            data.regularisation_weight,
            policy.state_degree,
            dirichlet_boundary_ids);
          metric = std::make_shared<dealii_backend::MassMetric>(
            coefficient->parameter_l2_metric(policy.control_metric_solve));
          if (has_constraint)
            {
              box_data = make_cellwise_box_data(
                metric->layout(),
                *bounds,
                metric,
                specification.formulation.constraint_id,
                "compiler.v1.cellwise_parameter_box");
              constraint =
                std::make_shared<dealii_backend::CellwiseBoxConstraint>(
                  make_constraint(*box_data, as_mass_metric(*metric)));
              constraint_realisation =
                ConstraintRealisation::cellwise_parameter_l2;
            }
          solvers = contract::StateAdjointSolversT<Backend>{
            [coefficient, solve_policy = policy.state_solve](
              const contract::PrimalBlockT<Backend> &parameter) {
              return coefficient->solve_state_with_report(parameter,
                                                          solve_policy);
            },
            [coefficient, solve_policy = policy.adjoint_solve](
              const contract::PrimalBlockT<Backend> &full_point,
              const contract::CovectorBlockT<Backend> &state_rhs) {
              return coefficient->solve_adjoint_with_report(full_point,
                                                            state_rhs,
                                                            solve_policy);
            }};
          executable = coefficient;
        }
      else if (uses_assembled_v1_target)
        {
          using AssembledModel = detail::ScalarComponentModel<dim>;
          std::shared_ptr<AssembledModel> assembled;
          if (uses_general_scalar)
            assembled = std::make_shared<AssembledModel>(
              triangulation,
              data.forcing,
              data.desired_state,
              *data.general_scalar,
              data.regularisation_weight,
              policy.state_degree,
              *scalar_plan);
          else
            assembled = std::make_shared<AssembledModel>(
              triangulation,
              data.forcing,
              data.desired_state,
              data.fixed_dirichlet_data,
              *data.diffusion,
              data.reaction,
              data.regularisation_weight,
              policy.state_degree,
              *scalar_plan);
          reduced_hessian = assembled;
          metric = std::make_shared<dealii_backend::MassMetric>(
            assembled->control_l2_metric(policy.control_metric_solve));
          const bool scalar_plan_has_constraint =
            scalar_plan->constraint == ScalarConstraintOperatorKind::cellwise_box;
          if (scalar_plan_has_constraint)
            {
              box_data = make_cellwise_box_data(
                metric->layout(),
                *bounds,
                metric,
                specification.formulation.constraint_id,
                "compiler.v1.cellwise_box");
              constraint =
                std::make_shared<dealii_backend::CellwiseBoxConstraint>(
                  make_constraint(*box_data, as_mass_metric(*metric)));
              constraint_realisation = ConstraintRealisation::cellwise_l2;
            }
          solvers = contract::StateAdjointSolversT<Backend>{
            [assembled, solve_policy = policy.state_solve](
              const contract::PrimalBlockT<Backend> &control) {
              return assembled->solve_state_with_report(control, solve_policy);
            },
            [assembled, solve_policy = policy.adjoint_solve](
              const contract::PrimalBlockT<Backend> &full_point,
              const contract::CovectorBlockT<Backend> &state_rhs) {
              return assembled->solve_adjoint_with_report(full_point,
                                                          state_rhs,
                                                          solve_policy);
            }};
          executable = assembled;
        }
      else
        {
          using DirectModel = dealii_backend::ScalarDiffusionReactionModel<dim>;
          const auto direct = std::make_shared<DirectModel>(
            triangulation,
            data.forcing,
            data.desired_state,
            *data.diffusion,
            data.reaction,
            data.regularisation_weight,
            policy.state_degree,
            dirichlet_boundary_ids);
          reduced_hessian = direct;
          metric = std::make_shared<dealii_backend::MassMetric>(
            direct->control_l2_metric(policy.control_metric_solve));
          if (has_constraint)
            {
              box_data = make_cellwise_box_data(
                metric->layout(),
                *bounds,
                metric,
                specification.formulation.constraint_id,
                "compiler.v1.cellwise_box");
              constraint =
                std::make_shared<dealii_backend::CellwiseBoxConstraint>(
                  make_constraint(*box_data, as_mass_metric(*metric)));
              constraint_realisation = ConstraintRealisation::cellwise_l2;
            }
          solvers = contract::StateAdjointSolversT<Backend>{
            [direct, solve_policy = policy.state_solve](
              const contract::PrimalBlockT<Backend> &control) {
              return direct->solve_state_with_report(control, solve_policy);
            },
            [direct, solve_policy = policy.adjoint_solve](
              const contract::PrimalBlockT<Backend> &full_point,
              const contract::CovectorBlockT<Backend> &state_rhs) {
              return direct->solve_adjoint_with_report(full_point,
                                                       state_rhs,
                                                       solve_policy);
            }};
          executable = direct;
          if (uses_supplied_otd)
            supplied_otd_system = std::make_shared<
              const contract::SuppliedOTDSystemT<Backend>>(
              DirectModel::make_supplied_otd_system(
                direct,
                *specification.supplied_otd_declaration,
                lifetime_owner));
        }
      auto finalized_decision = finalize_resolved_decision<dim>(
        policy,
        constraint_realisation,
        resolved_decision,
        request,
        *executable,
        *metric,
        scalar_plan ? &*scalar_plan : nullptr,
        supplied_otd_system ? supplied_otd_system.get() : nullptr,
        specification.supplied_otd_declaration
          ? &*specification.supplied_otd_declaration
          : nullptr);
      std::shared_ptr<const contract::EqualityConstrainedQuadraticKKTProductT<Backend>>
        kkt_product;
      std::shared_ptr<const contract::BoxComplementarityT<Backend>>
        pdas_complementarity;
      if (product == CompilationProduct::quadratic_kkt)
        {
          kkt_product = make_compiled_dto_kkt_product<Backend>(executable);
          finalized_decision.kkt_record = make_kkt_record(*kkt_product);
        }
      else if (product == CompilationProduct::pdas)
        {
          if (supplied_otd_system)
            kkt_product = std::make_shared<const
              contract::EqualityConstrainedQuadraticKKTProductT<Backend>>(
              contract::make_canonical_supplied_otd_kkt_product(
                *supplied_otd_system));
          else
            kkt_product = make_compiled_dto_kkt_product<Backend>(executable);
          const auto *mass_metric =
            dynamic_cast<const dealii_backend::MassMetric *>(metric.get());
          if (mass_metric == nullptr ||
              !mass_metric->supports_coefficientwise_box_projection())
            {
              result.diagnostics.add(
                semantic::v1::DiagnosticCategory::formulation_capability,
                specification.formulation.metric_id,
                "compiled_pdas_metric_realisation",
                "Select a concrete positive-diagonal cellwise L2 metric before constructing the PDAS multiplier conversion.");
              return result;
            }
          contract::require(static_cast<bool>(box_data),
                            "Compiled PDAS needs its shared cellwise box data");
          pdas_complementarity = make_pdas_complementarity(
            *kkt_product, *box_data, metric);
          finalized_decision.kkt_record = make_kkt_record(
            *kkt_product, static_cast<bool>(supplied_otd_system));
          finalized_decision.pdas_record = make_pdas_record(
            *kkt_product,
            *pdas_complementarity,
            *box_data,
            *metric,
            policy,
            static_cast<bool>(supplied_otd_system));
        }
      CompilationManifest manifest = make_manifest(finalized_decision);
      if (box_data)
        {
          manifest.constraint_record.box_data_token = box_data->token_id();
          manifest.constraint_record.bounds_digest = box_data->bounds_digest();
          manifest.constraint_record.control_layout =
            box_data->layout_signature();
          manifest.constraint_record.metric_identity = metric->id();
          manifest.constraint_record.data_provenance =
            box_data->data_provenance();
          manifest.resolved_decision.constraint_record =
            manifest.constraint_record;
        }
      if (pdas_complementarity)
        result.pdas_problem = std::make_shared<const
          CompiledPDASProblemT<Backend>>(
          std::move(kkt_product),
          std::move(pdas_complementarity),
          metric,
          1,
          policy.pdas,
          policy.pdas_kkt_solver,
          manifest,
          std::move(lifetime_owner),
          std::move(box_data),
          std::move(constraint));
      else if (kkt_product)
        result.kkt_problem = std::make_shared<
          const CompiledQuadraticKKTProblemT<Backend>>(
          std::move(kkt_product),
          manifest,
          std::move(lifetime_owner));
      else if (supplied_otd_system)
        result.supplied_otd_problem = std::make_shared<
          const CompiledSuppliedOTDProblemT<Backend>>(
          std::move(supplied_otd_system),
          manifest,
          std::move(lifetime_owner));
      else
        result.problem = std::make_shared<const CompiledProblemT<Backend>>(
          executable,
          metric,
          constraint,
          solvers,
          manifest,
          std::move(lifetime_owner),
          std::move(reduced_hessian),
          std::move(box_data));
      return result;
    }

  private:
    enum class DirichletControlRegistration
    {
      complete_nodal_l2,
      partial_nodal_l2,
      l2_transposition,
      hhalf_control,
      h1_tracking_hhalf_control,
      h1_control
    };

    enum class ConstraintRealisation
    {
      none,
      cellwise_l2,
      cellwise_parameter_l2,
      facewise_l2
    };

    enum class CompiledTargetKind
    {
      direct_volume,
      assembled_volume,
      neumann_boundary,
      weighted_boundary_trace,
      pure_neumann,
      dirichlet_control,
      l2_dirichlet_transposition,
      hhalf_dirichlet_control,
      h1_tracking_hhalf_dirichlet_control,
      h1_dirichlet_control,
      h1_control_l2_metric,
      h1_control_h1_metric,
      hminus1_control_metric,
      continuous_control_l2_metric,
      coefficient_identification,
      general_scalar_robin,
      point_sensor,
      normal_flux
    };

    static const semantic::v1::RegionSpec *
    find_region(const semantic::v1::ProblemSpec &specification,
                const std::string &              id)
    {
      const auto region = std::find_if(
        specification.regions.begin(),
        specification.regions.end(),
        [&id](const semantic::v1::RegionSpec &candidate) {
          return candidate.id == id;
        });
      return region == specification.regions.end() ? nullptr : &*region;
    }

    static const semantic::v1::VariableSpec *
    find_variable(const semantic::v1::ProblemSpec &specification,
                  const std::string &              id)
    {
      const auto variable = std::find_if(
        specification.variables.begin(),
        specification.variables.end(),
        [&id](const semantic::v1::VariableSpec &candidate) {
          return candidate.id == id;
        });
      return variable == specification.variables.end() ? nullptr : &*variable;
    }

    static const semantic::v1::TransformationSpec *
    find_transformation(const semantic::v1::ProblemSpec &specification,
                        const std::string &              id)
    {
      const auto transformation = std::find_if(
        specification.transformations.begin(),
        specification.transformations.end(),
        [&id](const semantic::v1::TransformationSpec &candidate) {
          return candidate.id == id;
        });
      return transformation == specification.transformations.end()
               ? nullptr
               : &*transformation;
    }

    static const semantic::v1::ObservationSpec *
    find_observation(const semantic::v1::ProblemSpec &specification,
                     const std::string &              id)
    {
      const auto observation = std::find_if(
        specification.observations.begin(),
        specification.observations.end(),
        [&id](const semantic::v1::ObservationSpec &candidate) {
          return candidate.id == id;
        });
      return observation == specification.observations.end() ? nullptr :
                                                               &*observation;
    }

    static const semantic::v1::RegionSpec *
    selected_tracking_region(const semantic::v1::ProblemSpec &specification)
    {
      const auto loss = std::find_if(
        specification.losses.begin(),
        specification.losses.end(),
        [](const semantic::v1::LossSpec &candidate) {
          return candidate.kind == semantic::v1::LossKind::quadratic_tracking;
        });
      if (loss == specification.losses.end())
        return nullptr;
      const auto observation = find_observation(specification,
                                                loss->source_observation_id);
      return observation == nullptr ? nullptr :
                                    find_region(specification,
                                                observation->region_id);
    }

    static std::set<dealii::types::material_id>
    selected_tracking_material_ids(const semantic::v1::RegionSpec &region)
    {
      std::set<dealii::types::material_id> ids;
      for (const auto id : region.material_ids)
        ids.insert(static_cast<dealii::types::material_id>(id));
      return ids;
    }

    template <int dim>
    static void
    validate_point_sensor_mesh(
      const dealii::Triangulation<dim> &triangulation,
      const semantic::v1::RegionSpec &  region,
      semantic::v1::ValidationReport &   report)
    {
      dealii::MappingQ1<dim> mapping;
      for (std::size_t index = 0; index < region.point_coordinates.size(); ++index)
        {
          const auto &coordinate = region.point_coordinates[index];
          if (coordinate.size() != dim)
            {
              report.add(semantic::v1::DiagnosticCategory::lowerability,
                         region.id,
                         "point_sensor_coordinate_dimension",
                         "Provide sensor coordinates with the compiled mesh dimension.");
              continue;
            }
          dealii::Point<dim> point;
          for (unsigned int component = 0; component < dim; ++component)
            point[component] = coordinate[component];
          const auto cell_and_reference_point =
            dealii::GridTools::find_active_cell_around_point(
              mapping, triangulation, point);
          if (cell_and_reference_point.first == triangulation.end())
            report.add(semantic::v1::DiagnosticCategory::lowerability,
                       region.id,
                       "point_sensor_inside_mesh",
                       "Place every immutable sensor coordinate inside the compiled mesh.");
        }
    }

    static bool
    uses_neumann_control(const semantic::v1::ProblemSpec &specification)
    {
      return std::any_of(
        specification.residual_terms.begin(),
        specification.residual_terms.end(),
        [](const semantic::v1::ResidualTermSpec &term) {
          return term.kind == semantic::v1::ResidualTermKind::neumann_control;
        });
    }

    static bool
    uses_neumann_conservative_transport(
      const semantic::v1::ProblemSpec &specification)
    {
      return uses_neumann_control(specification) &&
             std::any_of(
               specification.residual_terms.begin(),
               specification.residual_terms.end(),
               [](const semantic::v1::ResidualTermSpec &term) {
                 return term.kind ==
                        semantic::v1::ResidualTermKind::conservative_transport;
               });
    }

    static bool
    uses_mean_zero_multiplier(
      const semantic::v1::ProblemSpec &specification)
    {
      return std::any_of(
        specification.requirement_policies.begin(),
        specification.requirement_policies.end(),
        [&specification](const semantic::v1::RequirementPolicySpec &policy) {
          return policy.subject_id == specification.formulation.state_variable_id &&
                 policy.kind ==
                   semantic::v1::RequirementKind::mean_zero_multiplier &&
                 policy.status ==
                   semantic::v1::RequirementStatus::selected_discrete_realisation;
        });
    }

    static bool
    uses_h1_control_regularisation_loss(
      const semantic::v1::ProblemSpec &specification)
    {
      return std::any_of(
        specification.losses.begin(),
        specification.losses.end(),
        [](const semantic::v1::LossSpec &loss) {
          return loss.kind ==
                 semantic::v1::LossKind::quadratic_h1_control_regularisation;
        });
    }

    static bool
    uses_hhalf_control_regularisation_loss(
      const semantic::v1::ProblemSpec &specification)
    {
      return std::any_of(
        specification.losses.begin(),
        specification.losses.end(),
        [](const semantic::v1::LossSpec &loss) {
          return loss.kind == semantic::v1::LossKind::
                                quadratic_hhalf_control_regularisation;
        });
    }

    static bool
    has_homogeneous_dirichlet_continuous_control(
      const semantic::v1::ProblemSpec &specification)
    {
      const auto control = find_variable(
        specification, specification.formulation.control_variable_id);
      if (control == nullptr)
        return false;
      const auto space = std::find_if(
        specification.spaces.begin(),
        specification.spaces.end(),
        [control](const semantic::v1::SpaceSpec &candidate) {
          return candidate.id == control->space_id;
        });
      if (space == specification.spaces.end() ||
          space->topology != semantic::v1::SpaceTopology::h1)
        return false;
      return std::any_of(
        specification.requirement_policies.begin(),
        specification.requirement_policies.end(),
        [control](const semantic::v1::RequirementPolicySpec &policy) {
          return policy.subject_id == control->id &&
                 policy.kind == semantic::v1::RequirementKind::fixed_dirichlet &&
                 policy.status == semantic::v1::RequirementStatus::
                                    selected_discrete_realisation;
        });
    }

    static bool
    uses_parameter_diffusion_residual(
      const semantic::v1::ProblemSpec &specification)
    {
      return std::any_of(
        specification.residual_terms.begin(),
        specification.residual_terms.end(),
        [](const semantic::v1::ResidualTermSpec &term) {
          return term.kind ==
                 semantic::v1::ResidualTermKind::parameter_diffusion_reaction;
        });
    }

    static bool
    has_h1_state_observation(
      const semantic::v1::ProblemSpec &specification)
    {
      return std::any_of(
        specification.observations.begin(),
        specification.observations.end(),
        [](const semantic::v1::ObservationSpec &observation) {
          return observation.kind ==
                 semantic::v1::ObservationKind::h1_state_restriction;
        });
    }

    static bool
    has_weighted_boundary_trace(
      const semantic::v1::ProblemSpec &specification)
    {
      return std::any_of(
        specification.observations.begin(),
        specification.observations.end(),
        [](const semantic::v1::ObservationSpec &observation) {
          return observation.kind ==
                 semantic::v1::ObservationKind::weighted_boundary_trace;
        });
    }

    static bool
    has_point_sensor_observation(
      const semantic::v1::ProblemSpec &specification)
    {
      return std::any_of(
        specification.observations.begin(),
        specification.observations.end(),
        [](const semantic::v1::ObservationSpec &observation) {
          return observation.kind ==
                 semantic::v1::ObservationKind::point_sensor;
        });
    }

    static bool
    has_normal_flux_observation(
      const semantic::v1::ProblemSpec &specification)
    {
      return std::any_of(
        specification.observations.begin(),
        specification.observations.end(),
        [](const semantic::v1::ObservationSpec &observation) {
          return observation.kind ==
                 semantic::v1::ObservationKind::normal_flux;
        });
    }

    static bool
    uses_general_scalar_residual(
      const semantic::v1::ProblemSpec &specification)
    {
      return std::any_of(
        specification.residual_terms.begin(),
        specification.residual_terms.end(),
        [](const semantic::v1::ResidualTermSpec &term) {
          return term.kind ==
                 semantic::v1::ResidualTermKind::tensor_diffusion;
        });
    }

    static const semantic::v1::RegionSpec *
    selected_robin_boundary_region(
      const semantic::v1::ProblemSpec &specification)
    {
      const auto term = std::find_if(
        specification.residual_terms.begin(),
        specification.residual_terms.end(),
        [](const semantic::v1::ResidualTermSpec &candidate) {
          return candidate.kind ==
                 semantic::v1::ResidualTermKind::robin_bilinear;
        });
      return term == specification.residual_terms.end()
               ? nullptr
               : find_region(specification, term->region_id);
    }

    static bool
    selects_h1_control_metric(
      const semantic::v1::ProblemSpec &specification)
    {
      const auto metric = std::find_if(
        specification.metrics.begin(),
        specification.metrics.end(),
        [&specification](const semantic::v1::MetricSpec &candidate) {
          return candidate.id == specification.formulation.metric_id;
        });
      return metric != specification.metrics.end() &&
             metric->kind == semantic::v1::MetricKind::h1;
    }

    static bool
    selects_hhalf_control_metric(
      const semantic::v1::ProblemSpec &specification)
    {
      const auto metric = std::find_if(
        specification.metrics.begin(),
        specification.metrics.end(),
        [&specification](const semantic::v1::MetricSpec &candidate) {
          return candidate.id == specification.formulation.metric_id;
        });
      return metric != specification.metrics.end() &&
             metric->kind == semantic::v1::MetricKind::hhalf;
    }

    static bool
    selects_hminus1_control_metric(
      const semantic::v1::ProblemSpec &specification)
    {
      const auto metric = std::find_if(
        specification.metrics.begin(),
        specification.metrics.end(),
        [&specification](const semantic::v1::MetricSpec &candidate) {
          return candidate.id == specification.formulation.metric_id;
        });
      return metric != specification.metrics.end() &&
             metric->kind == semantic::v1::MetricKind::hminus1;
    }

    static const semantic::v1::RegionSpec *
    selected_neumann_control_region(
      const semantic::v1::ProblemSpec &specification)
    {
      const auto term = std::find_if(
        specification.residual_terms.begin(),
        specification.residual_terms.end(),
        [](const semantic::v1::ResidualTermSpec &candidate) {
          return candidate.kind == semantic::v1::ResidualTermKind::neumann_control;
        });
      return term == specification.residual_terms.end()
               ? nullptr
               : find_region(specification, term->region_id);
    }

    static const semantic::v1::RegionSpec *
    selected_dirichlet_control_region(
      const semantic::v1::ProblemSpec &specification)
    {
      const auto policy = std::find_if(
        specification.requirement_policies.begin(),
        specification.requirement_policies.end(),
        [&specification](const semantic::v1::RequirementPolicySpec &candidate) {
          return candidate.subject_id ==
                   specification.formulation.state_variable_id &&
                 candidate.kind ==
                   semantic::v1::RequirementKind::controlled_dirichlet;
        });
      return policy == specification.requirement_policies.end()
               ? nullptr
               : find_region(specification, policy->region_id);
    }

    static std::set<dealii::types::boundary_id>
    boundary_ids(const semantic::v1::RegionSpec &region)
    {
      std::set<dealii::types::boundary_id> ids;
      for (const auto id : region.boundary_ids)
        ids.insert(static_cast<dealii::types::boundary_id>(id));
      return ids;
    }

    static std::set<dealii::types::boundary_id>
    boundary_ids_from_plan(const std::set<unsigned int> &ids)
    {
      std::set<dealii::types::boundary_id> result;
      for (const auto id : ids)
        result.insert(static_cast<dealii::types::boundary_id>(id));
      return result;
    }

    static bool
    forms_declared_boundary_partition(
      const semantic::v1::RegionSpec &fixed_region,
      const semantic::v1::RegionSpec &controlled_region)
    {
      if (fixed_region.kind != semantic::v1::RegionKind::boundary ||
          controlled_region.kind != semantic::v1::RegionKind::boundary ||
          fixed_region.boundary_ids.empty() ||
          controlled_region.boundary_ids.empty())
        return false;
      return std::none_of(
        fixed_region.boundary_ids.begin(), fixed_region.boundary_ids.end(),
        [&controlled_region](const unsigned int id) {
          return std::find(controlled_region.boundary_ids.begin(),
                           controlled_region.boundary_ids.end(),
                           id) != controlled_region.boundary_ids.end();
        });
    }

    static bool
    uses_fixed_dirichlet_reconstruction(
      const semantic::v1::ProblemSpec &specification)
    {
      const auto state = find_variable(specification,
                                       specification.formulation.state_variable_id);
      if (state == nullptr || state->physical_field_transform_id.empty())
        return false;
      const auto transformation = find_transformation(
        specification, state->physical_field_transform_id);
      return transformation != nullptr &&
             transformation->kind ==
               semantic::v1::TransformationKind::fixed_dirichlet_reconstruction;
    }

    static bool
    uses_dirichlet_control_lifting(
      const semantic::v1::ProblemSpec &specification)
    {
      const auto state = find_variable(specification,
                                       specification.formulation.state_variable_id);
      if (state == nullptr || state->physical_field_transform_id.empty())
        return false;
      const auto transformation = find_transformation(
        specification, state->physical_field_transform_id);
      return transformation != nullptr &&
             transformation->kind ==
               semantic::v1::TransformationKind::dirichlet_control_lifting;
    }

    static bool
    uses_l2_dirichlet_transposition(
      const semantic::v1::ProblemSpec &specification)
    {
      return std::any_of(
        specification.residual_terms.begin(),
        specification.residual_terms.end(),
        [](const semantic::v1::ResidualTermSpec &term) {
          return term.kind ==
                   semantic::v1::ResidualTermKind::transposition_laplacian ||
                 term.kind == semantic::v1::ResidualTermKind::
                                dirichlet_transposition_control;
        });
    }

    static bool
    uses_normalized_dirichlet_laplace_control(
      const semantic::v1::ProblemSpec &specification)
    {
      return uses_dirichlet_control_lifting(specification) &&
             std::any_of(
               specification.residual_terms.begin(),
               specification.residual_terms.end(),
               [](const semantic::v1::ResidualTermSpec &term) {
                 return term.kind ==
                        semantic::v1::ResidualTermKind::laplacian;
               });
    }

    static bool
    uses_dirichlet_control_target(
      const semantic::v1::ProblemSpec &specification)
    {
      return uses_dirichlet_control_lifting(specification) ||
             uses_l2_dirichlet_transposition(specification);
    }

    static bool
    uses_partial_dirichlet_control_lifting(
      const semantic::v1::ProblemSpec &specification)
    {
      if (!uses_dirichlet_control_lifting(specification))
        return false;
      const auto state = find_variable(
        specification, specification.formulation.state_variable_id);
      const auto transformation = state == nullptr ? nullptr : find_transformation(
        specification, state->physical_field_transform_id);
      return transformation != nullptr && !transformation->fixed_data_id.empty();
    }

    static CompiledTargetKind
    target_kind_for_dirichlet_registration(
      const DirichletControlRegistration registration)
    {
      switch (registration)
        {
          case DirichletControlRegistration::complete_nodal_l2:
          case DirichletControlRegistration::partial_nodal_l2:
            return CompiledTargetKind::dirichlet_control;
          case DirichletControlRegistration::l2_transposition:
            return CompiledTargetKind::l2_dirichlet_transposition;
          case DirichletControlRegistration::hhalf_control:
            return CompiledTargetKind::hhalf_dirichlet_control;
          case DirichletControlRegistration::h1_tracking_hhalf_control:
            return CompiledTargetKind::h1_tracking_hhalf_dirichlet_control;
          case DirichletControlRegistration::h1_control:
            return CompiledTargetKind::h1_dirichlet_control;
        }
      contract::require(false, "Unknown Dirichlet control registration");
      return CompiledTargetKind::dirichlet_control;
    }

    static bool
    has_residual_signature(
      const semantic::v1::ResolvedProblemView &resolved,
      const std::vector<semantic::v1::ResidualTermKind> &expected_kinds)
    {
      const auto &specification = resolved.specification();
      if (specification.residual_terms.size() != expected_kinds.size() ||
          specification.equations.size() != 1)
        return false;
      const auto &equation = resolved.equation(specification.formulation.equation_id);
      if (equation.residual_term_ids.size() != expected_kinds.size())
        return false;
      for (const auto &term_id : equation.residual_term_ids)
        {
          const auto &term = resolved.residual_term(term_id);
          if (term.equation_id != equation.id ||
              std::find(expected_kinds.begin(), expected_kinds.end(), term.kind) ==
                expected_kinds.end())
            return false;
        }
      for (const auto kind : expected_kinds)
        if (std::count_if(
              specification.residual_terms.begin(),
              specification.residual_terms.end(),
              [kind](const semantic::v1::ResidualTermSpec &term) {
                return term.kind == kind;
              }) != 1)
          return false;
      return true;
    }

    static bool
    matches_fractional_metric_registration(
      const ResolvedCompilationRequest &request,
      const semantic::v1::MetricSpec &  metric,
      const semantic::v1::VariableSpec &control,
      const semantic::v1::VariableSpec &state)
    {
      if (!request.fractional_metric_selection)
        return false;
      const auto &selection = *request.fractional_metric_selection;
      return selection.metric_id == metric.id &&
             selection.control_space_id == control.space_id &&
             selection.volume_space_id == state.space_id &&
             selection.operator_realisation ==
               semantic::v1::FractionalTraceOperatorRealisation::
                 volume_mass_plus_stiffness_schur &&
             selection.apply_realisation ==
               semantic::v1::FractionalTraceApplyRealisation::
                 minimum_h1_extension &&
             selection.inverse_realisation ==
               semantic::v1::FractionalTraceInverseRealisation::
                 full_volume_operator_inverse;
    }

    static bool
    matches_boundary_h1_metric_registration(
      const ResolvedCompilationRequest &request,
      const semantic::v1::MetricSpec &  metric,
      const semantic::v1::VariableSpec &control,
      const semantic::v1::RegionSpec &  controlled_region)
    {
      if (!request.boundary_h1_metric_selection)
        return false;
      const auto &selection = *request.boundary_h1_metric_selection;
      return selection.metric_id == metric.id &&
             selection.control_space_id == control.space_id &&
             selection.boundary_region_id == controlled_region.id &&
             selection.operator_realisation ==
               semantic::v1::BoundaryH1MetricOperatorRealisation::
                 boundary_mass_plus_tangential_stiffness &&
             selection.tangential_gradient_realisation ==
               semantic::v1::BoundaryH1TangentialGradientRealisation::
                 projected_ambient_gradient &&
             selection.nullspace_realisation ==
               semantic::v1::BoundaryH1MetricNullspaceRealisation::
                 positive_mass_no_nullspace;
    }

    static bool
    matches_transposition_registration(
      const ResolvedCompilationRequest &request,
      const semantic::v1::ResolvedProblemView &resolved,
      const semantic::v1::EquationBlockSpec &equation,
      const semantic::v1::ObservationSpec &control_observation,
      const semantic::v1::VariableSpec &control)
    {
      if (!request.transposition_selection)
        return false;
      const auto &selection = *request.transposition_selection;
      return selection.subject_equation_id == equation.id &&
             selection.strong_space_id == equation.test_space_id &&
             selection.observation_id == control_observation.id &&
             selection.transpose_source_space_id ==
               control_observation.output_space_id &&
             selection.continuous_parent_space_id == control.space_id &&
             !selection.conforming_trace_space_id.empty() &&
             !selection.equivalence_policy_id.empty() &&
             !selection.conormal_policy_id.empty() &&
             selection.operator_realisation ==
               semantic::v1::TranspositionOperatorRealisation::
                 scalar_diffusion_reaction_dirichlet_laplacian &&
             selection.discrete_realisation ==
               semantic::v1::TranspositionDiscreteRealisation::
                 conforming_nodal_lifting_equivalence &&
             selection.equivalence_realisation ==
               semantic::v1::TranspositionEquivalenceRealisation::
                 conforming_lifting_variational_equivalence &&
             resolved.space(selection.conforming_trace_space_id).topology ==
               semantic::v1::SpaceTopology::hhalf;
    }

    static bool
    matches_partial_registration(
      const ResolvedCompilationRequest &request,
      const semantic::v1::ResolvedProblemView &resolved,
      const semantic::v1::TransformationSpec &transformation,
      const semantic::v1::RegionSpec &controlled_region)
    {
      if (!request.partial_boundary_selection)
        return false;
      const auto &selection = *request.partial_boundary_selection;
      return selection.transformation_id == transformation.id &&
             selection.controlled_boundary_region_id == controlled_region.id &&
             !selection.fixed_boundary_region_id.empty() &&
             selection.requires_complete_exterior &&
             selection.requires_disjoint_regions &&
             selection.interface_realisation ==
               semantic::v1::PartialDirichletInterfaceRealisation::
                 fixed_data_precedence &&
             selection.trace_realisation ==
               semantic::v1::PartialDirichletTraceRealisation::
                 relative_interior_nodal_zero_endpoint &&
             selection.hanging_realisation ==
               semantic::v1::PartialDirichletHangingRealisation::unsupported &&
             resolved.region(selection.fixed_boundary_region_id).kind ==
               semantic::v1::RegionKind::boundary;
    }

    static std::optional<DirichletControlRegistration>
    resolve_dirichlet_control_registration(
      const semantic::v1::ResolvedProblemView &resolved,
      const ResolvedCompilationRequest &        request)
    {
      const auto &specification = resolved.specification();
      if (!uses_dirichlet_control_target(specification))
        return std::nullopt;

      const auto state = find_variable(
        specification, specification.formulation.state_variable_id);
      const auto control = find_variable(
        specification, specification.formulation.control_variable_id);
      if (state == nullptr || control == nullptr ||
          state->role != semantic::v1::VariableRole::state ||
          control->role != semantic::v1::VariableRole::control ||
          specification.formulation.constraint_id.size() != 0 ||
          !specification.constraints.empty() ||
          specification.observations.size() != 2 ||
          specification.losses.size() != 2 || specification.metrics.size() != 1)
        return std::nullopt;

      const auto controlled_policy = std::find_if(
        specification.requirement_policies.begin(),
        specification.requirement_policies.end(),
        [&specification](const semantic::v1::RequirementPolicySpec &policy) {
          return policy.subject_id == specification.formulation.state_variable_id &&
                 policy.kind == semantic::v1::RequirementKind::controlled_dirichlet &&
                 policy.status ==
                   semantic::v1::RequirementStatus::selected_discrete_realisation;
        });
      if (controlled_policy == specification.requirement_policies.end())
        return std::nullopt;
      const auto *controlled_region =
        find_region(specification, controlled_policy->region_id);
      if (controlled_region == nullptr ||
          controlled_region->kind != semantic::v1::RegionKind::boundary ||
          controlled_region->boundary_ids.empty())
        return std::nullopt;

      const auto equation = std::find_if(
        specification.equations.begin(),
        specification.equations.end(),
        [&specification](const semantic::v1::EquationBlockSpec &candidate) {
          return candidate.id == specification.formulation.equation_id;
        });
      if (equation == specification.equations.end())
        return std::nullopt;

      const auto state_space = std::find_if(
        specification.spaces.begin(),
        specification.spaces.end(),
        [state](const semantic::v1::SpaceSpec &space) {
          return space.id == state->space_id;
        });
      const auto test_space = std::find_if(
        specification.spaces.begin(),
        specification.spaces.end(),
        [equation](const semantic::v1::SpaceSpec &space) {
          return space.id == equation->test_space_id;
        });
      const auto control_space = std::find_if(
        specification.spaces.begin(),
        specification.spaces.end(),
        [control](const semantic::v1::SpaceSpec &space) {
          return space.id == control->space_id;
        });
      if (state_space == specification.spaces.end() ||
          test_space == specification.spaces.end() ||
          control_space == specification.spaces.end() ||
          test_space->region_id != state_space->region_id ||
          !find_region(specification, state_space->region_id)->is_full_domain)
        return std::nullopt;

      const auto state_observation = std::find_if(
        specification.observations.begin(),
        specification.observations.end(),
        [state](const semantic::v1::ObservationSpec &observation) {
          return observation.input_variable_id == state->id &&
                 (observation.kind ==
                    semantic::v1::ObservationKind::volume_restriction ||
                  observation.kind ==
                    semantic::v1::ObservationKind::h1_state_restriction);
        });
      const auto control_observation = std::find_if(
        specification.observations.begin(),
        specification.observations.end(),
        [control](const semantic::v1::ObservationSpec &observation) {
          return observation.input_variable_id == control->id &&
                 observation.kind ==
                   semantic::v1::ObservationKind::boundary_restriction;
        });
      if (state_observation == specification.observations.end() ||
          control_observation == specification.observations.end() ||
          state_observation->region_id != state_space->region_id ||
          control_observation->region_id != controlled_region->id)
        return std::nullopt;

      const auto state_observation_space = std::find_if(
        specification.spaces.begin(),
        specification.spaces.end(),
        [state_observation](const semantic::v1::SpaceSpec &space) {
          return space.id == state_observation->output_space_id;
        });
      const auto control_observation_space = std::find_if(
        specification.spaces.begin(),
        specification.spaces.end(),
        [control_observation](const semantic::v1::SpaceSpec &space) {
          return space.id == control_observation->output_space_id;
        });
      if (state_observation_space == specification.spaces.end() ||
          control_observation_space == specification.spaces.end() ||
          state_observation_space->region_id != state_observation->region_id ||
          control_observation_space->region_id != controlled_region->id ||
          control_space->region_id != controlled_region->id)
        return std::nullopt;

      const auto tracking_loss = std::find_if(
        specification.losses.begin(), specification.losses.end(),
        [state_observation](const semantic::v1::LossSpec &loss) {
          return loss.kind == semantic::v1::LossKind::quadratic_tracking &&
                 loss.source_observation_id == state_observation->id;
        });
      const auto control_loss = std::find_if(
        specification.losses.begin(), specification.losses.end(),
        [control_observation](const semantic::v1::LossSpec &loss) {
          return loss.kind != semantic::v1::LossKind::quadratic_tracking &&
                 loss.source_observation_id == control_observation->id;
        });
      const auto metric = std::find_if(
        specification.metrics.begin(), specification.metrics.end(),
        [&specification](const semantic::v1::MetricSpec &candidate) {
          return candidate.id == specification.formulation.metric_id;
        });
      if (tracking_loss == specification.losses.end() ||
          control_loss == specification.losses.end() ||
          metric == specification.metrics.end() ||
          tracking_loss->pairing_id != state_observation->output_pairing_id ||
          control_loss->pairing_id != control_observation->output_pairing_id ||
          metric->variable_id != control->id)
        return std::nullopt;

      const auto &state_observation_pairing =
        resolved.pairing(state_observation->output_pairing_id);
      const auto &control_observation_pairing =
        resolved.pairing(control_observation->output_pairing_id);
      const auto &control_pairing = resolved.pairing(metric->pairing_id);
      if (state_observation_pairing.primal_space_id !=
            state_observation->output_space_id ||
          state_observation_pairing.covector_space_id !=
            state_observation->output_space_id ||
          control_observation_pairing.primal_space_id !=
            control_observation->output_space_id ||
          control_observation_pairing.covector_space_id !=
            control_observation->output_space_id ||
          control_pairing.primal_space_id != control->space_id ||
          control_pairing.covector_space_id != control->space_id)
        return std::nullopt;

      if (uses_l2_dirichlet_transposition(specification))
        {
          if (!specification.transformations.empty() ||
              state->physical_field_transform_id.size() != 0 ||
              state_space->topology != semantic::v1::SpaceTopology::l2 ||
              test_space->topology != semantic::v1::SpaceTopology::h2 ||
              control_space->topology != semantic::v1::SpaceTopology::l2 ||
              control_observation_space->topology !=
                semantic::v1::SpaceTopology::l2 ||
              state_observation_space->topology !=
                semantic::v1::SpaceTopology::l2 ||
              metric->kind != semantic::v1::MetricKind::l2 ||
              control_loss->kind !=
                semantic::v1::LossKind::quadratic_control_regularisation ||
              !has_residual_signature(
                resolved,
                {semantic::v1::ResidualTermKind::transposition_laplacian,
                 semantic::v1::ResidualTermKind::volume_source,
                 semantic::v1::ResidualTermKind::
                   dirichlet_transposition_control}) ||
              !matches_transposition_registration(request,
                                                   resolved,
                                                   *equation,
                                                   *control_observation,
                                                   *control))
            return std::nullopt;
          return DirichletControlRegistration::l2_transposition;
        }

      if (specification.transformations.size() != 1 ||
          state->physical_field_transform_id.empty())
        return std::nullopt;
      const auto &transformation =
        resolved.transformation(state->physical_field_transform_id);
      if (transformation.kind !=
            semantic::v1::TransformationKind::dirichlet_control_lifting ||
          transformation.input_variable_id != state->id ||
          transformation.output_space_id != state->space_id ||
          transformation.control_variable_id != control->id)
        return std::nullopt;

      const bool partial = !transformation.fixed_data_id.empty();
      if (partial != uses_partial_dirichlet_control_lifting(specification))
        return std::nullopt;
      if (partial)
        {
          if (!matches_partial_registration(request,
                                            resolved,
                                            transformation,
                                            *controlled_region))
            return std::nullopt;
        }
      else if (request.partial_boundary_selection)
        return std::nullopt;

      if (state_space->topology != semantic::v1::SpaceTopology::h1 ||
          test_space->topology != semantic::v1::SpaceTopology::h1)
        return std::nullopt;

      const bool normalized = has_residual_signature(
        resolved,
        {semantic::v1::ResidualTermKind::laplacian,
         semantic::v1::ResidualTermKind::volume_source});
      const bool diffusion_reaction = has_residual_signature(
        resolved,
        {semantic::v1::ResidualTermKind::diffusion_reaction,
         semantic::v1::ResidualTermKind::volume_source});
      if (!normalized && !diffusion_reaction)
        return std::nullopt;

      if (diffusion_reaction &&
          control_space->topology == semantic::v1::SpaceTopology::h1 &&
          control_observation_space->topology ==
            semantic::v1::SpaceTopology::h1 &&
          state_observation->kind ==
            semantic::v1::ObservationKind::volume_restriction &&
          control_loss->kind ==
            semantic::v1::LossKind::quadratic_control_regularisation &&
          metric->kind == semantic::v1::MetricKind::l2)
        return partial ? DirichletControlRegistration::partial_nodal_l2
                       : DirichletControlRegistration::complete_nodal_l2;

      if (!normalized || partial)
        return std::nullopt;

      if (control_space->topology == semantic::v1::SpaceTopology::hhalf &&
          control_observation_space->topology ==
            semantic::v1::SpaceTopology::hhalf &&
          state_observation->kind ==
            semantic::v1::ObservationKind::volume_restriction &&
          control_loss->kind ==
            semantic::v1::LossKind::quadratic_hhalf_control_regularisation &&
          metric->kind == semantic::v1::MetricKind::hhalf &&
          matches_fractional_metric_registration(request,
                                                 *metric,
                                                 *control,
                                                 *state))
        return DirichletControlRegistration::hhalf_control;

      if (control_space->topology == semantic::v1::SpaceTopology::hhalf &&
          control_observation_space->topology ==
            semantic::v1::SpaceTopology::l2 &&
          state_observation->kind ==
            semantic::v1::ObservationKind::h1_state_restriction &&
          control_loss->kind ==
            semantic::v1::LossKind::quadratic_control_regularisation &&
          metric->kind == semantic::v1::MetricKind::hhalf &&
          matches_fractional_metric_registration(request,
                                                 *metric,
                                                 *control,
                                                 *state) &&
          request.h1_target_data_membership_selection &&
          request.h1_target_data_membership_selection->data_id ==
            tracking_loss->data_id &&
          request.h1_target_data_membership_selection->observation_space_id ==
            state_observation->output_space_id &&
          request.h1_target_data_membership_selection->fixed_boundary_region_id ==
            controlled_region->id &&
          request.h1_target_data_membership_selection->regularity_realisation ==
            semantic::v1::H1TargetDataRegularityRealisation::
              h1_value_and_weak_gradient &&
          request.h1_target_data_membership_selection->trace_realisation ==
            semantic::v1::H1TargetDataTraceRealisation::
              zero_trace_on_fixed_boundary)
        return DirichletControlRegistration::h1_tracking_hhalf_control;

      if (control_space->topology == semantic::v1::SpaceTopology::h1 &&
          control_observation_space->topology ==
            semantic::v1::SpaceTopology::h1 &&
          state_observation->kind ==
            semantic::v1::ObservationKind::volume_restriction &&
          control_loss->kind ==
            semantic::v1::LossKind::quadratic_h1_control_regularisation &&
          metric->kind == semantic::v1::MetricKind::h1 &&
          matches_boundary_h1_metric_registration(request,
                                                  *metric,
                                                  *control,
                                                  *controlled_region))
        return DirichletControlRegistration::h1_control;

      return std::nullopt;
    }

    static ResolvedDirichletRegistration
    resolved_dirichlet_registration(
      const std::optional<DirichletControlRegistration> &registration)
    {
      if (!registration)
        return ResolvedDirichletRegistration::none;
      switch (*registration)
        {
          case DirichletControlRegistration::complete_nodal_l2:
            return ResolvedDirichletRegistration::complete_nodal_l2;
          case DirichletControlRegistration::partial_nodal_l2:
            return ResolvedDirichletRegistration::partial_nodal_l2;
          case DirichletControlRegistration::l2_transposition:
            return ResolvedDirichletRegistration::l2_transposition;
          case DirichletControlRegistration::hhalf_control:
            return ResolvedDirichletRegistration::hhalf_control;
          case DirichletControlRegistration::h1_tracking_hhalf_control:
            return ResolvedDirichletRegistration::h1_tracking_hhalf_control;
          case DirichletControlRegistration::h1_control:
            return ResolvedDirichletRegistration::h1_control;
        }
      contract::require(false, "Unknown resolved Dirichlet registration");
      return ResolvedDirichletRegistration::none;
    }

    static std::optional<DirichletControlRegistration>
    dirichlet_registration_from_request(
      const ResolvedCompilationRequest &request)
    {
      switch (request.dirichlet_registration)
        {
          case ResolvedDirichletRegistration::complete_nodal_l2:
            return DirichletControlRegistration::complete_nodal_l2;
          case ResolvedDirichletRegistration::partial_nodal_l2:
            return DirichletControlRegistration::partial_nodal_l2;
          case ResolvedDirichletRegistration::l2_transposition:
            return DirichletControlRegistration::l2_transposition;
          case ResolvedDirichletRegistration::hhalf_control:
            return DirichletControlRegistration::hhalf_control;
          case ResolvedDirichletRegistration::h1_tracking_hhalf_control:
            return DirichletControlRegistration::h1_tracking_hhalf_control;
          case ResolvedDirichletRegistration::h1_control:
            return DirichletControlRegistration::h1_control;
          case ResolvedDirichletRegistration::none:
            return std::nullopt;
        }
      contract::require(false, "Unknown request Dirichlet registration");
      return std::nullopt;
    }

    static CompiledTargetKind
    target_kind_from_request(const ResolvedCompilationRequest &request)
    {
      switch (request.target_family)
        {
          case ResolvedTargetFamily::direct_volume:
            return CompiledTargetKind::direct_volume;
          case ResolvedTargetFamily::assembled_volume:
            return CompiledTargetKind::assembled_volume;
          case ResolvedTargetFamily::neumann_boundary:
            return CompiledTargetKind::neumann_boundary;
          case ResolvedTargetFamily::weighted_boundary_trace:
            return CompiledTargetKind::weighted_boundary_trace;
          case ResolvedTargetFamily::pure_neumann:
            return CompiledTargetKind::pure_neumann;
          case ResolvedTargetFamily::dirichlet_control:
            return CompiledTargetKind::dirichlet_control;
          case ResolvedTargetFamily::l2_dirichlet_transposition:
            return CompiledTargetKind::l2_dirichlet_transposition;
          case ResolvedTargetFamily::hhalf_dirichlet_control:
            return CompiledTargetKind::hhalf_dirichlet_control;
          case ResolvedTargetFamily::h1_tracking_hhalf_dirichlet_control:
            return CompiledTargetKind::h1_tracking_hhalf_dirichlet_control;
          case ResolvedTargetFamily::h1_dirichlet_control:
            return CompiledTargetKind::h1_dirichlet_control;
          case ResolvedTargetFamily::h1_control_l2_metric:
            return CompiledTargetKind::h1_control_l2_metric;
          case ResolvedTargetFamily::h1_control_h1_metric:
            return CompiledTargetKind::h1_control_h1_metric;
          case ResolvedTargetFamily::hminus1_control_metric:
            return CompiledTargetKind::hminus1_control_metric;
          case ResolvedTargetFamily::continuous_control_l2_metric:
            return CompiledTargetKind::continuous_control_l2_metric;
          case ResolvedTargetFamily::coefficient_identification:
            return CompiledTargetKind::coefficient_identification;
          case ResolvedTargetFamily::general_scalar_robin:
            return CompiledTargetKind::general_scalar_robin;
          case ResolvedTargetFamily::point_sensor:
            return CompiledTargetKind::point_sensor;
          case ResolvedTargetFamily::normal_flux:
            return CompiledTargetKind::normal_flux;
          case ResolvedTargetFamily::unresolved:
            break;
        }
      contract::require(false, "The compiler request has no resolved target family");
      return CompiledTargetKind::direct_volume;
    }

    static void
    close_compilation_request(
      const semantic::v1::ResolvedProblemView &resolved,
      ResolvedCompilationRequest &              request,
      const std::optional<DirichletControlRegistration> &registration)
    {
      const auto &specification = resolved.specification();
      request.dirichlet_registration =
        resolved_dirichlet_registration(registration);
      request.uses_fixed_reconstruction =
        uses_fixed_dirichlet_reconstruction(specification);
      request.uses_dirichlet_control =
        uses_dirichlet_control_target(specification);
      request.uses_l2_dirichlet_control =
        uses_l2_dirichlet_transposition(specification);
      request.uses_normalized_dirichlet_laplace =
        uses_normalized_dirichlet_laplace_control(specification);
      request.uses_normalized_laplacian =
        request.uses_l2_dirichlet_control ||
        request.uses_normalized_dirichlet_laplace;
      request.uses_partial_dirichlet_control =
        uses_partial_dirichlet_control_lifting(specification);
      request.uses_neumann_boundary_control = uses_neumann_control(specification);
      request.uses_neumann_convection =
        uses_neumann_conservative_transport(specification);
      request.uses_mean_zero_gauge = uses_mean_zero_multiplier(specification);
      request.uses_h1_control_regularisation_loss =
        uses_h1_control_regularisation_loss(specification);
      request.uses_hhalf_control_regularisation_loss =
        uses_hhalf_control_regularisation_loss(specification);
      request.uses_h1_control_regularisation =
        request.uses_h1_control_regularisation_loss &&
        !request.uses_dirichlet_control;
      request.uses_h1_control_metric = selects_h1_control_metric(specification);
      request.uses_hhalf_control_metric =
        selects_hhalf_control_metric(specification);
      request.uses_hminus1_control_metric =
        selects_hminus1_control_metric(specification);
      request.uses_homogeneous_dirichlet_continuous_control =
        has_homogeneous_dirichlet_continuous_control(specification);
      request.uses_coefficient_identification =
        uses_parameter_diffusion_residual(specification);
      request.uses_general_scalar = uses_general_scalar_residual(specification);
      request.uses_h1_state_observation = has_h1_state_observation(specification);
      request.uses_weighted_boundary_trace = has_weighted_boundary_trace(specification);
      request.uses_point_sensor = has_point_sensor_observation(specification);
      request.uses_normal_flux = has_normal_flux_observation(specification);
      request.uses_h1_dirichlet_control =
        request.dirichlet_registration == ResolvedDirichletRegistration::h1_control;
      request.uses_hhalf_dirichlet_registration =
        request.dirichlet_registration == ResolvedDirichletRegistration::hhalf_control ||
        request.dirichlet_registration ==
          ResolvedDirichletRegistration::h1_tracking_hhalf_control;
      request.uses_h1_tracking_hhalf_dirichlet_registration =
        request.dirichlet_registration ==
        ResolvedDirichletRegistration::h1_tracking_hhalf_control;

      const auto *tracking_region = selected_tracking_region(specification);
      const auto *robin_region = selected_robin_boundary_region(specification);
      request.tracking_region_id =
        tracking_region == nullptr ? std::string{} : tracking_region->id;
      request.robin_boundary_region_id =
        robin_region == nullptr ? std::string{} : robin_region->id;
      const auto mean_policy = std::find_if(
        specification.requirement_policies.begin(),
        specification.requirement_policies.end(),
        [&specification](const semantic::v1::RequirementPolicySpec &candidate) {
          return candidate.subject_id ==
                   specification.formulation.state_variable_id &&
                 candidate.kind ==
                   semantic::v1::RequirementKind::mean_zero_multiplier &&
                 candidate.status ==
                   semantic::v1::RequirementStatus::selected_discrete_realisation;
        });
      request.mean_zero_region_id =
        mean_policy == specification.requirement_policies.end()
          ? std::string{}
          : mean_policy->region_id;
      if (request.partial_boundary_selection)
        {
          request.partial_fixed_boundary_region_id =
            request.partial_boundary_selection->fixed_boundary_region_id;
          request.partial_control_boundary_region_id =
            request.partial_boundary_selection->controlled_boundary_region_id;
        }
      const auto fixed_policy = std::find_if(
        specification.requirement_policies.begin(),
        specification.requirement_policies.end(),
        [&specification, &request](const semantic::v1::RequirementPolicySpec &candidate) {
          return candidate.kind == semantic::v1::RequirementKind::fixed_dirichlet &&
                 candidate.status ==
                   semantic::v1::RequirementStatus::selected_discrete_realisation &&
                 (candidate.subject_id == specification.formulation.state_variable_id ||
                  (request.uses_partial_dirichlet_control &&
                   candidate.subject_id == "dirichlet_control_lifting"));
        });
      if (fixed_policy != specification.requirement_policies.end())
        request.fixed_boundary_region_id = fixed_policy->region_id;
      const auto controlled_policy = std::find_if(
        specification.requirement_policies.begin(),
        specification.requirement_policies.end(),
        [&specification](const semantic::v1::RequirementPolicySpec &candidate) {
          return candidate.subject_id == specification.formulation.state_variable_id &&
                 candidate.kind == semantic::v1::RequirementKind::controlled_dirichlet &&
                 candidate.status ==
                   semantic::v1::RequirementStatus::selected_discrete_realisation;
        });
      if (request.uses_dirichlet_control &&
          controlled_policy != specification.requirement_policies.end())
        request.control_boundary_region_id = controlled_policy->region_id;
      else
        {
          const auto *neumann_region = selected_neumann_control_region(specification);
          request.control_boundary_region_id =
            neumann_region == nullptr ? std::string{} : neumann_region->id;
        }
      if (request.uses_partial_dirichlet_control &&
          !request.partial_control_boundary_region_id.empty())
        request.control_boundary_region_id =
          request.partial_control_boundary_region_id;
      request.uses_subdomain_observation =
        tracking_region != nullptr && !tracking_region->is_full_domain;
      request.uses_assembled_v1_target =
        !request.uses_neumann_boundary_control &&
        !request.uses_h1_control_regularisation &&
        !request.uses_homogeneous_dirichlet_continuous_control &&
        !request.uses_coefficient_identification &&
        !request.uses_dirichlet_control &&
        (request.uses_fixed_reconstruction ||
         request.uses_subdomain_observation || request.uses_general_scalar ||
         request.uses_h1_state_observation || request.uses_point_sensor ||
         request.uses_normal_flux);

      if (request.uses_mean_zero_gauge)
        request.target_family = ResolvedTargetFamily::pure_neumann;
      else if (request.uses_weighted_boundary_trace)
        request.target_family = ResolvedTargetFamily::weighted_boundary_trace;
      else if (request.uses_neumann_boundary_control)
        request.target_family = ResolvedTargetFamily::neumann_boundary;
      else if (request.uses_dirichlet_control)
        switch (request.dirichlet_registration)
          {
            case ResolvedDirichletRegistration::complete_nodal_l2:
            case ResolvedDirichletRegistration::partial_nodal_l2:
              request.target_family = ResolvedTargetFamily::dirichlet_control;
              break;
            case ResolvedDirichletRegistration::l2_transposition:
              request.target_family =
                ResolvedTargetFamily::l2_dirichlet_transposition;
              break;
            case ResolvedDirichletRegistration::hhalf_control:
              request.target_family =
                ResolvedTargetFamily::hhalf_dirichlet_control;
              break;
            case ResolvedDirichletRegistration::h1_tracking_hhalf_control:
              request.target_family =
                ResolvedTargetFamily::h1_tracking_hhalf_dirichlet_control;
              break;
            case ResolvedDirichletRegistration::h1_control:
              request.target_family = ResolvedTargetFamily::h1_dirichlet_control;
              break;
            case ResolvedDirichletRegistration::none:
              request.target_family = ResolvedTargetFamily::unresolved;
              break;
          }
      else if (request.uses_coefficient_identification)
        request.target_family = ResolvedTargetFamily::coefficient_identification;
      else if (request.uses_hminus1_control_metric)
        request.target_family = ResolvedTargetFamily::hminus1_control_metric;
      else if (request.uses_h1_control_regularisation)
        request.target_family = request.uses_h1_control_metric
                                  ? ResolvedTargetFamily::h1_control_h1_metric
                                  : ResolvedTargetFamily::h1_control_l2_metric;
      else if (request.uses_homogeneous_dirichlet_continuous_control)
        request.target_family =
          ResolvedTargetFamily::continuous_control_l2_metric;
      else if (request.uses_point_sensor)
        request.target_family = ResolvedTargetFamily::point_sensor;
      else if (request.uses_normal_flux)
        request.target_family = ResolvedTargetFamily::normal_flux;
      else if (request.uses_assembled_v1_target)
        request.target_family = request.uses_general_scalar
                                  ? ResolvedTargetFamily::general_scalar_robin
                                  : ResolvedTargetFamily::assembled_volume;
      else
        request.target_family = ResolvedTargetFamily::direct_volume;
    }

    static void
    validate_dirichlet_control_registration(
      const semantic::v1::ResolvedProblemView &resolved,
      const ResolvedCompilationRequest &        request,
      semantic::v1::ValidationReport &          report)
    {
      const auto &specification = resolved.specification();
      if (request.uses_dirichlet_control &&
          request.dirichlet_registration == ResolvedDirichletRegistration::none)
        report.add(
          semantic::v1::DiagnosticCategory::lowerability,
          specification.id,
          "section_5_11_registered_signature",
          "Select one registered Dirichlet composition: complete or partial nodal L2 lifting, conforming L2 transposition, Section 5.11.1 option 1 or 2, or Section 5.11.3 tangential H1 control.");
    }

    template <int dim>
    static bool
    controls_complete_exterior_boundary(
      const dealii::Triangulation<dim> &              triangulation,
      const std::set<dealii::types::boundary_id> &controlled_ids)
    {
      bool has_boundary_face = false;
      for (auto cell = triangulation.begin_active();
           cell != triangulation.end();
           ++cell)
        for (unsigned int face = 0;
             face < dealii::GeometryInfo<dim>::faces_per_cell;
             ++face)
          if (cell->face(face)->at_boundary())
            {
              has_boundary_face = true;
              if (controlled_ids.count(cell->face(face)->boundary_id()) == 0)
                return false;
            }
      return has_boundary_face;
    }

    static bool
    valid_bound_representation(const CellwiseBoxDataBindings &bounds)
    {
      return (std::holds_alternative<double>(bounds.lower) &&
              std::holds_alternative<double>(bounds.upper)) ||
             (std::holds_alternative<dealii::Vector<double>>(bounds.lower) &&
              std::holds_alternative<dealii::Vector<double>>(bounds.upper));
    }

    static bool
    valid_metric_solve_policy(
      const dealii_backend::MassMetricSolveParameters &policy)
    {
      return policy.maximum_iterations > 0 &&
             std::isfinite(policy.relative_tolerance) &&
             policy.relative_tolerance > 0.0 &&
             std::isfinite(policy.absolute_tolerance) &&
             policy.absolute_tolerance > 0.0;
    }

    template <int dim>
    static std::size_t
    count_boundary_faces(
      const dealii::Triangulation<dim> &              triangulation,
      const std::set<dealii::types::boundary_id> &boundary_ids)
    {
      std::size_t count = 0;
      for (auto cell = triangulation.begin_active();
           cell != triangulation.end();
           ++cell)
        for (unsigned int face = 0;
             face < dealii::GeometryInfo<dim>::faces_per_cell;
             ++face)
          if (cell->face(face)->at_boundary() &&
              boundary_ids.count(cell->face(face)->boundary_id()) != 0)
            ++count;
      return count;
    }

    template <int dim>
    static bool
    contains_all_boundary_ids(
      const dealii::Triangulation<dim> &              triangulation,
      const std::set<dealii::types::boundary_id> &requested_ids)
    {
      std::set<dealii::types::boundary_id> found_ids;
      for (auto cell = triangulation.begin_active();
           cell != triangulation.end();
           ++cell)
        for (unsigned int face = 0;
             face < dealii::GeometryInfo<dim>::faces_per_cell;
             ++face)
          if (cell->face(face)->at_boundary() &&
              requested_ids.count(cell->face(face)->boundary_id()) != 0)
            found_ids.insert(cell->face(face)->boundary_id());
      return found_ids == requested_ids;
    }

    template <int dim>
    static bool
    contains_all_material_ids(
      const dealii::Triangulation<dim> &              triangulation,
      const std::set<dealii::types::material_id> &requested_ids)
    {
      std::set<dealii::types::material_id> found_ids;
      for (auto cell = triangulation.begin_active();
           cell != triangulation.end();
           ++cell)
        if (requested_ids.count(cell->material_id()) != 0)
          found_ids.insert(cell->material_id());
      return found_ids == requested_ids;
    }

    template <int dim>
    static bool
    forms_complete_boundary_partition(
      const dealii::Triangulation<dim> &              triangulation,
      const std::set<dealii::types::boundary_id> &fixed_ids,
      const std::set<dealii::types::boundary_id> &robin_ids)
    {
      if (std::any_of(fixed_ids.begin(), fixed_ids.end(), [&robin_ids](auto id) {
            return robin_ids.count(id) != 0;
          }))
        return false;
      bool has_boundary_face = false;
      for (auto cell = triangulation.begin_active();
           cell != triangulation.end();
           ++cell)
        for (unsigned int face = 0;
             face < dealii::GeometryInfo<dim>::faces_per_cell;
             ++face)
          if (cell->face(face)->at_boundary())
            {
              has_boundary_face = true;
              const auto id = cell->face(face)->boundary_id();
              if (fixed_ids.count(id) + robin_ids.count(id) != 1)
                return false;
            }
      return has_boundary_face;
    }

    static void
    validate_cellwise_bound_values(
      const CellwiseBoxDataBindings &bounds,
      const std::size_t              expected_size,
      const std::string &            component_id,
      semantic::v1::ValidationReport &report)
    {
      validate_bound_values(bounds.lower,
                            bounds.upper,
                            expected_size,
                            component_id,
                            "cellwise_bound_layout",
                            report);
    }

    static void
    validate_facewise_bound_values(
      const FacewiseBoxDataBindings &bounds,
      const std::size_t              expected_size,
      const std::string &            component_id,
      semantic::v1::ValidationReport &report)
    {
      validate_bound_values(bounds.lower,
                            bounds.upper,
                            expected_size,
                            component_id,
                            "facewise_bound_layout",
                            report);
    }

    template <typename BoundValue>
    static void
    validate_bound_values(const BoundValue &lower,
                          const BoundValue &upper,
                          const std::size_t expected_size,
                          const std::string &component_id,
                          const std::string &layout_capability,
                          semantic::v1::ValidationReport &report)
    {
      using semantic::v1::DiagnosticCategory;
      if (std::holds_alternative<double>(lower))
        {
          const double lower_value = std::get<double>(lower);
          const double upper_value = std::get<double>(upper);
          if (!std::isfinite(lower_value) || !std::isfinite(upper_value))
            report.add(DiagnosticCategory::lowerability,
                       component_id,
                       "finite_bound_values",
                       "Bind finite lower and upper values.");
          else if (lower_value > upper_value)
            report.add(DiagnosticCategory::lowerability,
                       component_id,
                       "ordered_bound_values",
                       "Bind lower values that do not exceed upper values.");
          return;
        }

      const auto &lower_values = std::get<dealii::Vector<double>>(lower);
      const auto &upper_values = std::get<dealii::Vector<double>>(upper);
      if (static_cast<std::size_t>(lower_values.size()) != expected_size ||
          static_cast<std::size_t>(upper_values.size()) != expected_size)
        {
          report.add(DiagnosticCategory::lowerability,
                     component_id,
                     layout_capability,
                     "Bind lower and upper vectors with the exact compiled decision layout.");
          return;
        }
      for (dealii::Vector<double>::size_type index = 0;
           index < lower_values.size();
           ++index)
        if (!std::isfinite(lower_values[index]) ||
            !std::isfinite(upper_values[index]))
          {
            report.add(DiagnosticCategory::lowerability,
                       component_id,
                       "finite_bound_values",
                       "Bind finite lower and upper values.");
            return;
          }
        else if (lower_values[index] > upper_values[index])
          {
            report.add(DiagnosticCategory::lowerability,
                       component_id,
                       "ordered_bound_values",
                       "Bind lower values that do not exceed upper values.");
            return;
          }
    }

    static bool
    has_strictly_positive_lower_bound(const CellwiseBoxDataBindings &bounds)
    {
      if (!valid_bound_representation(bounds))
        return false;
      if (std::holds_alternative<double>(bounds.lower))
        return std::isfinite(std::get<double>(bounds.lower)) &&
               std::get<double>(bounds.lower) > 0.0;
      const auto &lower = std::get<dealii::Vector<double>>(bounds.lower);
      return std::all_of(lower.begin(), lower.end(), [](const double value) {
        return std::isfinite(value) && value > 0.0;
      });
    }

    static bool
    valid_facewise_bound_representation(const FacewiseBoxDataBindings &bounds)
    {
      return (std::holds_alternative<double>(bounds.lower) &&
              std::holds_alternative<double>(bounds.upper)) ||
             (std::holds_alternative<dealii::Vector<double>>(bounds.lower) &&
              std::holds_alternative<dealii::Vector<double>>(bounds.upper));
    }

    static std::shared_ptr<const CompiledCellwiseBoxDataT<
      dealii_backend::SerialBackend>>
    make_cellwise_box_data(
      const contract::LayoutPtr &layout,
      const CellwiseBoxDataBindings &bounds,
      std::shared_ptr<const contract::MetricT<dealii_backend::SerialBackend>>
        metric,
      std::string semantic_id,
      std::string data_provenance)
    {
      using Backend = dealii_backend::SerialBackend;
      using Primal = contract::PrimalBlockT<Backend>;
      contract::require(static_cast<bool>(layout),
                        "Compiled cellwise box data needs a control layout");
      contract::require(valid_bound_representation(bounds),
                        "The compiled cellwise box needs compatible bound data");
      const std::size_t dimension = layout->dimension(0);
      dealii::Vector<double> lower = make_pdas_bound_vector(
        bounds.lower, dimension, "lower");
      dealii::Vector<double> upper = make_pdas_bound_vector(
        bounds.upper, dimension, "upper");
      return std::make_shared<const CompiledCellwiseBoxDataT<Backend>>(
        layout,
        Primal(layout, {std::move(lower)}),
        Primal(layout, {std::move(upper)}),
        std::move(metric),
        std::move(semantic_id),
        std::move(data_provenance));
    }

    static dealii_backend::CellwiseBoxConstraint
    make_constraint(
      const CompiledCellwiseBoxDataT<dealii_backend::SerialBackend> &box_data,
      const dealii_backend::MassMetric &                             projection_metric)
    {
      return dealii_backend::CellwiseBoxConstraint(
        box_data.layout(),
        box_data.lower().block(0),
        box_data.upper().block(0),
        projection_metric,
        box_data.token());
    }

    template <typename Model>
    static dealii_backend::FacewiseBoxConstraint
    make_facewise_constraint(const Model &                  executable,
                             const FacewiseBoxDataBindings &bounds,
                             const dealii_backend::MassMetric &projection_metric)
    {
      contract::require(valid_facewise_bound_representation(bounds),
                        "The v1 facewise box needs compatible bound data");
      if (std::holds_alternative<double>(bounds.lower))
        return executable.control_l2_box_constraint(
          std::get<double>(bounds.lower),
          std::get<double>(bounds.upper),
          projection_metric);
      return executable.control_l2_box_constraint(
        std::get<dealii::Vector<double>>(bounds.lower),
        std::get<dealii::Vector<double>>(bounds.upper),
        projection_metric);
    }

    void
    validate_lowerability(const semantic::v1::ProblemSpec & specification,
                          const ResolvedCompilationRequest &request,
                          const DealiiDiscretisationPolicy &policy,
                          semantic::v1::ValidationReport & report) const
    {
      using semantic::v1::DiagnosticCategory;
      if (policy.execution != DealiiDiscretisationPolicy::Execution::assembled)
        report.add(DiagnosticCategory::lowerability,
                   specification.id,
                   "assembled_execution",
                   "Select assembled execution; matrix-free lowering is not registered in v1.");
      if (policy.state_degree == 0)
        report.add(DiagnosticCategory::lowerability,
                   specification.id,
                   "FE_Q_state_degree",
                   "Select a scalar FE_Q state degree of at least one.");

      std::size_t full_volume_regions = 0;
      for (const auto &region : specification.regions)
        full_volume_regions += region.kind == semantic::v1::RegionKind::volume &&
                               region.is_full_domain;
      if (full_volume_regions != 1)
        report.add(DiagnosticCategory::lowerability,
                   specification.id,
                   "single_full_volume_region",
                   "The first v1 deal.II lowerer supports exactly one full volume region.");

      validate_registered_graph(specification, request, report);

      for (const auto &term : specification.residual_terms)
        if (!capabilities_.has_residual_term_lowerer(term.kind))
          report.add(DiagnosticCategory::lowerability,
                     term.id,
                     "registered_residual_term_lowerer",
                     "Register a lowerer for this residual term and its derivatives.");
      for (const auto &observation : specification.observations)
        if (!capabilities_.has_observation_lowerer(observation.kind))
          report.add(DiagnosticCategory::lowerability,
                     observation.id,
                     "registered_observation_lowerer",
                     "Register an observation value, JVP, and VJP lowerer.");
      for (const auto &loss : specification.losses)
        if (!capabilities_.has_loss_lowerer(loss.kind))
          report.add(DiagnosticCategory::lowerability,
                     loss.id,
                     "registered_loss_lowerer",
                     "Register a matching loss value and derivative lowerer.");
      for (const auto &metric : specification.metrics)
        if (!capabilities_.has_metric_lowerer(metric.kind))
          report.add(DiagnosticCategory::lowerability,
                     metric.id,
                     "registered_metric_lowerer",
                     "Register a metric realization with inverse apply.");
      for (const auto &constraint : specification.constraints)
        if (!capabilities_.has_constraint_lowerer(constraint.kind))
          report.add(DiagnosticCategory::lowerability,
                     constraint.id,
                     "registered_constraint_lowerer",
                     "Register the selected constraint projection realization.");
      for (const auto &transformation : specification.transformations)
        if (!capabilities_.has_transformation_lowerer(transformation.kind))
          report.add(DiagnosticCategory::lowerability,
                     transformation.id,
                     "registered_transformation_lowerer",
                     "Register value, JVP, and VJP lowering for this transformation.");

      const bool mean_zero_gauge = request.uses_mean_zero_gauge;
      const bool has_h1_control_regularisation_loss =
        request.uses_h1_control_regularisation_loss;
      const bool h1_control_regularisation =
        has_h1_control_regularisation_loss &&
        !request.uses_dirichlet_control;
      const bool hhalf_control_regularisation =
        request.uses_hhalf_control_regularisation_loss;
      const bool h1_control_metric = request.uses_h1_control_metric;
      const bool hhalf_control_metric =
        request.uses_hhalf_control_metric;
      const bool hminus1_control_metric =
        request.uses_hminus1_control_metric;
      const bool homogeneous_dirichlet_continuous_control =
        request.uses_homogeneous_dirichlet_continuous_control;
      const bool coefficient_identification =
        request.uses_coefficient_identification;
      const bool general_scalar =
        request.uses_general_scalar;
      const bool neumann_convection =
        request.uses_neumann_convection;
      const bool weighted_boundary_trace =
        request.uses_weighted_boundary_trace;
      const bool h1_state_observation =
        request.uses_h1_state_observation;
      const bool l2_dirichlet_transposition =
        request.uses_l2_dirichlet_control;
      const bool normalized_dirichlet_laplace =
        request.uses_normalized_dirichlet_laplace;
      const bool dirichlet_control_lifting =
        request.uses_dirichlet_control;
      const bool partial_dirichlet_control =
        request.uses_partial_dirichlet_control;
      const bool registered_h1_dirichlet_tracking =
        request.dirichlet_registration ==
          ResolvedDirichletRegistration::h1_tracking_hhalf_control;
      const auto *boundary = request.fixed_boundary_region_id.empty()
                               ? nullptr
                               : find_region(specification,
                                             request.fixed_boundary_region_id);
      const auto *controlled_boundary = request.control_boundary_region_id.empty()
                                         ? nullptr
                                         : find_region(
                                             specification,
                                             request.control_boundary_region_id);
      const auto *robin_boundary = request.robin_boundary_region_id.empty()
                                    ? nullptr
                                    : find_region(
                                        specification,
                                        request.robin_boundary_region_id);
      const auto *tracking_region = request.tracking_region_id.empty()
                                      ? nullptr
                                      : find_region(specification,
                                                    request.tracking_region_id);
      if (!mean_zero_gauge && !dirichlet_control_lifting &&
          (boundary == nullptr ||
           boundary->kind != semantic::v1::RegionKind::boundary ||
           boundary->boundary_ids.empty()))
        report.add(DiagnosticCategory::lowerability,
                   specification.formulation.state_variable_id,
                   "fixed_dirichlet_boundary_ids",
                   "Select a boundary region with at least one fixed Dirichlet id.");
      if (h1_state_observation)
        {
          if (tracking_region == nullptr || !tracking_region->is_full_domain)
            report.add(
              DiagnosticCategory::lowerability,
              specification.id,
              "h1_state_observation_full_domain",
              "Select the full volume region for the registered H1 state observation.");
          if (request.uses_fixed_reconstruction ||
              (dirichlet_control_lifting &&
               !registered_h1_dirichlet_tracking) ||
              request.uses_neumann_boundary_control ||
              h1_control_regularisation ||
              coefficient_identification ||
              general_scalar)
            report.add(
              DiagnosticCategory::lowerability,
              specification.id,
              "h1_state_observation_registered_combination",
              "Combine the first H1 state observation only with the registered homogeneous diffusion-reaction volume-control target and L2 control regularisation.");
        }
      if (dirichlet_control_lifting &&
          (controlled_boundary == nullptr ||
           controlled_boundary->kind != semantic::v1::RegionKind::boundary ||
           controlled_boundary->boundary_ids.empty()))
        report.add(
          DiagnosticCategory::lowerability,
          specification.formulation.state_variable_id,
          "controlled_dirichlet_boundary_ids",
          "Select a non-empty exterior boundary region for the registered Dirichlet-control lifting.");

      if (general_scalar)
        {
          const auto robin_source = std::find_if(
            specification.residual_terms.begin(),
            specification.residual_terms.end(),
            [](const semantic::v1::ResidualTermSpec &term) {
              return term.kind ==
                     semantic::v1::ResidualTermKind::robin_source;
            });
          if (robin_boundary == nullptr ||
              robin_boundary->kind != semantic::v1::RegionKind::boundary ||
              robin_boundary->boundary_ids.empty())
            report.add(
              DiagnosticCategory::lowerability,
              specification.id,
              "robin_boundary_ids",
              "Select one non-empty Robin boundary region for both Robin contributions.");
          if (robin_source == specification.residual_terms.end() ||
              robin_boundary == nullptr ||
              robin_source->region_id != robin_boundary->id)
            report.add(
              DiagnosticCategory::lowerability,
              specification.id,
              "shared_robin_boundary_region",
              "Place the Robin bilinear and source contributions on the same declared boundary region.");
          if (boundary != nullptr && robin_boundary != nullptr)
            for (const auto robin_id : robin_boundary->boundary_ids)
              if (std::find(boundary->boundary_ids.begin(),
                            boundary->boundary_ids.end(),
                            robin_id) != boundary->boundary_ids.end())
                report.add(DiagnosticCategory::lowerability,
                           robin_boundary->id,
                           "scalar_boundary_partition_overlap",
                           "Use disjoint fixed-Dirichlet and Robin boundary ids.");
          if (request.uses_fixed_reconstruction ||
              request.uses_neumann_boundary_control ||
              request.uses_mean_zero_gauge ||
              h1_control_regularisation || coefficient_identification ||
              dirichlet_control_lifting)
            report.add(
              DiagnosticCategory::lowerability,
              specification.id,
              "general_scalar_registered_combination",
              "Compose the first general scalar target only with homogeneous fixed Dirichlet data, volume control, full-volume tracking, and one Robin region.");
          const auto *general_tracking_region = tracking_region;
          if (general_tracking_region == nullptr ||
              !general_tracking_region->is_full_domain)
            report.add(
              DiagnosticCategory::lowerability,
              specification.id,
              "general_scalar_full_domain_tracking",
              "Select full-domain state tracking for the first general scalar target.");
        }

      if (mean_zero_gauge)
        {
          if (!request.uses_neumann_boundary_control)
            report.add(
              DiagnosticCategory::lowerability,
              specification.formulation.state_variable_id,
              "pure_neumann_registered_residual",
              "Select the mean-zero multiplier only with the registered pure-Neumann boundary-control residual.");
          const auto *mean_region = request.mean_zero_region_id.empty()
                                      ? nullptr
                                      : find_region(specification,
                                                    request.mean_zero_region_id);
          if (mean_region == nullptr ||
              mean_region->kind != semantic::v1::RegionKind::volume ||
              !mean_region->is_full_domain)
            report.add(
              DiagnosticCategory::lowerability,
              specification.formulation.state_variable_id,
              "pure_neumann_mean_constraint_region",
              "Place the mean-zero multiplier policy on the single full volume region.");
          if (boundary != nullptr)
            report.add(
              DiagnosticCategory::lowerability,
              specification.formulation.state_variable_id,
              "pure_neumann_without_fixed_dirichlet",
              "Do not declare a fixed Dirichlet policy with the pure-Neumann mean constraint.");
        }

      if (h1_control_regularisation)
        {
          const auto control = find_variable(
            specification, specification.formulation.control_variable_id);
          const auto space = std::find_if(
            specification.spaces.begin(),
            specification.spaces.end(),
            [control](const semantic::v1::SpaceSpec &candidate) {
              return control != nullptr && candidate.id == control->space_id;
            });
          if (control == nullptr ||
              space == specification.spaces.end() ||
              space->topology != semantic::v1::SpaceTopology::h1)
            report.add(
              DiagnosticCategory::lowerability,
              specification.formulation.control_variable_id,
              "h1_continuous_control_space",
              "Select the registered continuous H1 control space for H1 regularisation.");
          if (!specification.formulation.constraint_id.empty())
            report.add(
              DiagnosticCategory::lowerability,
              specification.formulation.constraint_id,
              "continuous_control_box_constraint",
              "Do not select the cellwise or facewise box with the continuous H1 control realization.");
          if (request.uses_fixed_reconstruction ||
              request.uses_neumann_boundary_control || dirichlet_control_lifting)
            report.add(
              DiagnosticCategory::lowerability,
              specification.id,
              "h1_regularisation_registered_combination",
              "The first H1-control regularisation target supports the homogeneous volume-control graph only.");
          if (tracking_region == nullptr || !tracking_region->is_full_domain)
            report.add(
              DiagnosticCategory::lowerability,
              specification.id,
              "h1_regularisation_full_domain_tracking",
              "The first H1-control regularisation target supports full-domain tracking only.");
        }

      if (h1_control_metric && !h1_control_regularisation &&
          !(normalized_dirichlet_laplace &&
            has_h1_control_regularisation_loss))
        report.add(
          DiagnosticCategory::lowerability,
          specification.formulation.metric_id,
          "h1_metric_registered_control_space",
          "Select the registered continuous H1-control target before requesting the H1 metric.");

      if (hminus1_control_metric)
        {
          if (h1_control_regularisation)
            report.add(
              DiagnosticCategory::lowerability,
              specification.id,
              "hminus1_metric_l2_control_loss",
              "Keep the registered L2 control loss when selecting the P5.2 H-1 metric target.");
          if (!h1_state_observation)
            report.add(
              DiagnosticCategory::lowerability,
              specification.id,
              "hminus1_metric_energy_observation",
              "Select the full-domain H1 state observation for the P5.2 H-1 metric target.");
          if (!specification.formulation.constraint_id.empty())
            report.add(
              DiagnosticCategory::lowerability,
              specification.formulation.constraint_id,
              "hminus1_metric_constraint",
              "Do not select a coefficientwise box for the coupled H-1 metric.");
          if (request.uses_fixed_reconstruction ||
              request.uses_neumann_boundary_control || dirichlet_control_lifting ||
              mean_zero_gauge || coefficient_identification || general_scalar)
            report.add(
              DiagnosticCategory::lowerability,
              specification.id,
              "hminus1_metric_registered_combination",
              "Compose the first H-1 metric target only with homogeneous fixed-Dirichlet state data and volume control.");
          if (tracking_region == nullptr || !tracking_region->is_full_domain)
            report.add(
              DiagnosticCategory::lowerability,
              specification.id,
              "hminus1_metric_full_domain_tracking",
              "Select full-domain energy tracking for the first H-1 metric target.");
        }

      if (homogeneous_dirichlet_continuous_control)
        {
          if (!h1_state_observation)
            report.add(
              DiagnosticCategory::lowerability,
              specification.id,
              "continuous_control_energy_observation",
              "Select the full-domain H1 state observation for the registered homogeneous-Dirichlet continuous-control target.");
          if (!specification.formulation.constraint_id.empty())
            report.add(
              DiagnosticCategory::lowerability,
              specification.formulation.constraint_id,
              "continuous_control_box_constraint",
              "Do not select a coefficientwise box on the continuous control layout.");
          if (request.uses_fixed_reconstruction ||
              request.uses_neumann_boundary_control || dirichlet_control_lifting ||
              mean_zero_gauge || coefficient_identification || general_scalar)
            report.add(
              DiagnosticCategory::lowerability,
              specification.id,
              "continuous_control_registered_combination",
              "Compose the first homogeneous-Dirichlet continuous-control target only with homogeneous fixed-Dirichlet state data and volume control.");
        }

      if (coefficient_identification)
        {
          const auto parameter = find_variable(
            specification, specification.formulation.control_variable_id);
          const auto space = std::find_if(
            specification.spaces.begin(),
            specification.spaces.end(),
            [parameter](const semantic::v1::SpaceSpec &candidate) {
              return parameter != nullptr && candidate.id == parameter->space_id;
            });
          if (parameter == nullptr ||
              parameter->role != semantic::v1::VariableRole::parameter ||
              space == specification.spaces.end() ||
              space->topology != semantic::v1::SpaceTopology::l2)
            report.add(
              DiagnosticCategory::lowerability,
              specification.formulation.control_variable_id,
              "cellwise_parameter_space",
              "Select the registered cellwise L2 parameter space for coefficient identification.");
          if (specification.formulation.constraint_id.empty())
            report.add(
              DiagnosticCategory::lowerability,
              specification.formulation.control_variable_id,
              "positive_parameter_constraint",
              "Select the registered positive cellwise parameter box.");
          if (request.uses_fixed_reconstruction ||
              request.uses_neumann_boundary_control || h1_control_regularisation ||
              dirichlet_control_lifting)
            report.add(
              DiagnosticCategory::lowerability,
              specification.id,
              "coefficient_identification_registered_combination",
              "The first coefficient-identification target supports the homogeneous full-domain volume graph only.");
          if (tracking_region == nullptr || !tracking_region->is_full_domain)
            report.add(
              DiagnosticCategory::lowerability,
              specification.id,
              "coefficient_identification_full_domain_tracking",
              "The first coefficient-identification target supports full-domain tracking only.");
        }

      if (dirichlet_control_lifting)
        {
          const auto control = find_variable(
            specification, specification.formulation.control_variable_id);
          const auto space = std::find_if(
            specification.spaces.begin(),
            specification.spaces.end(),
            [control](const semantic::v1::SpaceSpec &candidate) {
              return control != nullptr && candidate.id == control->space_id;
            });
          if (control == nullptr ||
              control->role != semantic::v1::VariableRole::control ||
              space == specification.spaces.end() ||
              space->topology != (l2_dirichlet_transposition
                                    ? semantic::v1::SpaceTopology::l2
                                  : hhalf_control_metric
                                    ? semantic::v1::SpaceTopology::hhalf
                                    : semantic::v1::SpaceTopology::h1) ||
              controlled_boundary == nullptr ||
              space->region_id != controlled_boundary->id)
            report.add(
              DiagnosticCategory::lowerability,
              specification.formulation.control_variable_id,
              l2_dirichlet_transposition
                ? "l2_dirichlet_control_parent_space"
              : hhalf_control_metric
                ? "hhalf_dirichlet_nodal_trace_control_space"
                : "dirichlet_nodal_trace_control_space",
              l2_dirichlet_transposition
                ? "Declare the continuous parent control in L2 on the controlled boundary; its conforming nodal trace is a selected discrete subspace."
              : hhalf_control_metric
                ? "Select the registered H1/2 nodal trace control space on the controlled Dirichlet boundary."
                : "Select the registered continuous nodal trace control space on the controlled Dirichlet boundary.");
          if (!specification.formulation.constraint_id.empty())
            report.add(
              DiagnosticCategory::lowerability,
              specification.formulation.constraint_id,
              "dirichlet_control_box_constraint",
              "The first nodal Dirichlet lifting has no box-constraint realization.");
          if (request.uses_fixed_reconstruction ||
              request.uses_neumann_boundary_control || h1_control_regularisation ||
              coefficient_identification)
            report.add(
              DiagnosticCategory::lowerability,
              specification.id,
              "dirichlet_control_registered_combination",
              "The first Dirichlet-control target supports only diffusion-reaction, volume forcing, full-volume tracking, and L2 trace regularisation.");
          if (normalized_dirichlet_laplace && partial_dirichlet_control)
            report.add(
              DiagnosticCategory::lowerability,
              specification.id,
              "normalized_dirichlet_complete_boundary",
              "Use the complete controlled exterior boundary for the registered Section 5.11 trace-space realizations.");
          if (hhalf_control_regularisation && !hhalf_control_metric)
            report.add(
              DiagnosticCategory::lowerability,
              specification.formulation.metric_id,
              "hhalf_loss_metric_realisation",
              "Select the minimum-extension H1/2 search metric with the H1/2 control loss.");
          if (partial_dirichlet_control)
            {
              if (boundary == nullptr ||
                  boundary->kind != semantic::v1::RegionKind::boundary ||
                  boundary->boundary_ids.empty())
                report.add(
                  DiagnosticCategory::lowerability,
                  specification.formulation.state_variable_id,
                  "partial_dirichlet_fixed_boundary",
                  "Declare a non-empty fixed Dirichlet boundary for the partial controlled lifting.");
              if (controlled_boundary == nullptr || boundary == nullptr ||
                  !forms_declared_boundary_partition(*boundary,
                                                     *controlled_boundary))
                report.add(
                  DiagnosticCategory::lowerability,
                  specification.formulation.state_variable_id,
                  "partial_dirichlet_boundary_partition",
                  "Declare disjoint fixed and controlled Dirichlet boundary regions for the partial lifting.");
              if (!request.partial_boundary_selection)
                report.add(
                  DiagnosticCategory::lowerability,
                  "dirichlet_control_lifting",
                  "partial_dirichlet_interface_policy",
                  "Declare the fixed/controlled corner and interface ownership policy.");
            }
          if (tracking_region == nullptr || !tracking_region->is_full_domain)
            report.add(
              DiagnosticCategory::lowerability,
              specification.id,
              "dirichlet_control_full_domain_tracking",
              "The first Dirichlet-control target supports full-domain state tracking only.");
        }

      if (weighted_boundary_trace &&
          (!request.uses_neumann_boundary_control || mean_zero_gauge))
        report.add(
          DiagnosticCategory::lowerability,
          specification.id,
          "weighted_boundary_trace_registered_combination",
          "Combine the first weighted boundary trace with the fixed-Dirichlet Neumann-control target.");

      if (!request.uses_neumann_boundary_control)
        return;
      const auto *control_region = controlled_boundary;
      if (control_region == nullptr ||
          control_region->kind != semantic::v1::RegionKind::boundary ||
          control_region->boundary_ids.empty())
        report.add(DiagnosticCategory::lowerability,
                   specification.id,
                   "neumann_control_boundary_ids",
                   "Select a marked boundary region with at least one Neumann control id.");
      if (tracking_region == nullptr ||
          (neumann_convection
             ? tracking_region->kind != semantic::v1::RegionKind::volume ||
                 tracking_region->is_full_domain ||
                 tracking_region->material_ids.empty()
             : tracking_region->kind != semantic::v1::RegionKind::boundary ||
                 tracking_region->boundary_ids.empty()))
        report.add(DiagnosticCategory::lowerability,
                   specification.id,
                   neumann_convection ? "neumann_convection_subdomain_region"
                                      : "boundary_tracking_region",
                   neumann_convection
                     ? "Select one or more material ids for the C5.6 state observation."
                     : "Select a marked boundary region for the state trace observation.");
      if (control_region != nullptr && boundary != nullptr)
        for (const auto control_id : control_region->boundary_ids)
          if (std::find(boundary->boundary_ids.begin(),
                        boundary->boundary_ids.end(),
                        control_id) != boundary->boundary_ids.end())
            report.add(DiagnosticCategory::lowerability,
                       control_region->id,
                       "neumann_control_dirichlet_overlap",
                       "Use boundary ids not fixed by the homogeneous Dirichlet realization.");
      if (neumann_convection && control_region != nullptr && boundary != nullptr &&
          !forms_declared_boundary_partition(*boundary, *control_region))
        report.add(DiagnosticCategory::lowerability,
                   specification.id,
                   "neumann_convection_boundary_partition",
                   "Declare disjoint fixed-Dirichlet and Neumann-control boundary regions for the C5.6 composition.");
    }

    static void
    validate_registered_graph(const semantic::v1::ProblemSpec &specification,
                              const ResolvedCompilationRequest &request,
                              semantic::v1::ValidationReport & report)
    {
      using semantic::v1::DataRole;
      using semantic::v1::DiagnosticCategory;
      using semantic::v1::LossKind;
      using semantic::v1::ResidualTermKind;

      const auto count_terms = [&specification](const ResidualTermKind kind) {
        return std::count_if(
          specification.residual_terms.begin(),
          specification.residual_terms.end(),
          [kind](const semantic::v1::ResidualTermSpec &term) {
            return term.kind == kind;
          });
      };
      const bool boundary_control = request.uses_neumann_boundary_control;
      const bool l2_dirichlet_transposition =
        request.uses_l2_dirichlet_control;
      const bool normalized_dirichlet_laplace =
        request.uses_normalized_dirichlet_laplace;
      const bool normalized_laplacian =
        request.uses_normalized_laplacian;
      const bool dirichlet_control =
        request.uses_dirichlet_control;
      const bool coefficient_identification =
        request.uses_coefficient_identification;
      const bool general_scalar =
        request.uses_general_scalar;
      const bool neumann_convection =
        request.uses_neumann_convection;
      const bool weighted_boundary_trace =
        request.uses_weighted_boundary_trace;
      const bool point_sensor = request.uses_point_sensor;
      const bool normal_flux = request.uses_normal_flux;
      const bool complete_residual = general_scalar
        ? count_terms(ResidualTermKind::tensor_diffusion) == 1 &&
            count_terms(ResidualTermKind::conservative_transport) == 1 &&
            count_terms(ResidualTermKind::advective_transport) == 1 &&
            count_terms(ResidualTermKind::reaction) == 1 &&
            count_terms(ResidualTermKind::volume_source) == 1 &&
            count_terms(ResidualTermKind::volume_control) == 1 &&
            count_terms(ResidualTermKind::robin_bilinear) == 1 &&
            count_terms(ResidualTermKind::robin_source) == 1 &&
            specification.residual_terms.size() == 8
        : coefficient_identification
        ? count_terms(ResidualTermKind::parameter_diffusion_reaction) == 1 &&
            count_terms(ResidualTermKind::volume_source) == 1 &&
            count_terms(ResidualTermKind::diffusion_reaction) == 0 &&
            count_terms(ResidualTermKind::volume_control) == 0 &&
            count_terms(ResidualTermKind::neumann_control) == 0
        : boundary_control
        ? (neumann_convection
             ? count_terms(ResidualTermKind::diffusion_reaction) == 1 &&
                 count_terms(ResidualTermKind::conservative_transport) == 1 &&
                 count_terms(ResidualTermKind::volume_source) == 1 &&
                 count_terms(ResidualTermKind::neumann_control) == 1 &&
                 count_terms(ResidualTermKind::volume_control) == 0 &&
                 specification.residual_terms.size() == 4
             : count_terms(ResidualTermKind::diffusion_reaction) == 1 &&
                 count_terms(ResidualTermKind::volume_source) == 1 &&
                 count_terms(ResidualTermKind::neumann_control) == 1 &&
                 count_terms(ResidualTermKind::volume_control) == 0)
        : l2_dirichlet_transposition
        ? count_terms(ResidualTermKind::transposition_laplacian) == 1 &&
            count_terms(
              ResidualTermKind::dirichlet_transposition_control) == 1 &&
            count_terms(ResidualTermKind::volume_source) == 1 &&
            specification.residual_terms.size() == 3
        : normalized_dirichlet_laplace
        ? count_terms(ResidualTermKind::laplacian) == 1 &&
            count_terms(ResidualTermKind::volume_source) == 1 &&
            specification.residual_terms.size() == 2
        : dirichlet_control
        ? count_terms(ResidualTermKind::diffusion_reaction) == 1 &&
            count_terms(ResidualTermKind::volume_source) == 1 &&
            count_terms(ResidualTermKind::volume_control) == 0 &&
            count_terms(ResidualTermKind::neumann_control) == 0
        : count_terms(ResidualTermKind::diffusion_reaction) == 1 &&
            count_terms(ResidualTermKind::volume_source) == 1 &&
            count_terms(ResidualTermKind::volume_control) == 1 &&
            count_terms(ResidualTermKind::neumann_control) == 0;
      if (!complete_residual)
        report.add(DiagnosticCategory::lowerability,
                   specification.id,
                   general_scalar
                     ? "complete_general_scalar_residual_term_set"
                     : coefficient_identification
                     ? "complete_parameter_diffusion_residual_term_set"
                     : boundary_control
                         ? neumann_convection
                             ? "complete_neumann_convection_residual_term_set"
                             : "complete_neumann_boundary_residual_term_set"
                         : l2_dirichlet_transposition
                             ? "complete_l2_dirichlet_transposition_term_set"
                         : normalized_dirichlet_laplace
                             ? "complete_normalized_dirichlet_laplace_term_set"
                         : dirichlet_control
                             ? "complete_dirichlet_control_residual_term_set"
                         : "complete_volume_residual_term_set",
                   general_scalar
                     ? "Declare exactly one tensor-diffusion, conservative-transport, advective-transport, reaction, volume-source, volume-control, Robin-bilinear, and Robin-source term."
                     : coefficient_identification
                     ? "Declare exactly one parameter diffusion-reaction and one volume-source term."
                     : boundary_control
                     ? neumann_convection
                         ? "Declare exactly one diffusion-reaction, conservative-transport, volume-source, and Neumann-control term."
                         : "Declare exactly one diffusion-reaction, volume-source, and Neumann-control term."
                     : l2_dirichlet_transposition
                     ? "Declare exactly one transposition Laplace state action, volume source, and Dirichlet normal-test-derivative control action."
                     : normalized_dirichlet_laplace
                     ? "Declare exactly one normalized Laplace state action and one volume-source term; the control enters through the declared lifting."
                     : dirichlet_control
                     ? "Declare exactly one diffusion-reaction and one volume-source term; the control enters through the declared lifting."
                     : "Declare exactly one diffusion-reaction, volume-source, and volume-control term.");

      const auto count_data = [&specification](const DataRole role) {
        return std::count_if(
          specification.data.begin(),
          specification.data.end(),
          [role](const semantic::v1::DataSpec &datum) {
            return datum.role == role;
          });
      };
      const bool complete_general_scalar_data =
        count_data(DataRole::diffusion) == 1 &&
        count_data(DataRole::conservative_transport) == 1 &&
        count_data(DataRole::advective_transport) == 1 &&
        count_data(DataRole::reaction) == 1 &&
        count_data(DataRole::robin_coefficient) == 1 &&
        count_data(DataRole::robin_source) == 1;
      if (count_data(DataRole::forcing) != 1 ||
          count_data(DataRole::desired_state) != 1 ||
          count_data(DataRole::observation_weight) !=
            (weighted_boundary_trace ? 1 : 0) ||
          count_data(DataRole::regularisation_weight) != 1 ||
          (normalized_laplacian
             ? count_data(DataRole::diffusion) != 0 ||
                 count_data(DataRole::reaction) != 0
           : general_scalar
             ? !complete_general_scalar_data
             : neumann_convection
               ? count_data(DataRole::diffusion) != 1 ||
                 count_data(DataRole::reaction) != 1 ||
                 count_data(DataRole::conservative_transport) != 1
             : (count_data(DataRole::reaction) != 1 ||
                (coefficient_identification
                   ? count_data(DataRole::diffusion) != 0
                   : count_data(DataRole::diffusion) != 1))))
        report.add(DiagnosticCategory::lowerability,
                   specification.id,
                   normalized_laplacian
                     ? "complete_normalized_laplacian_data_set"
                   : general_scalar
                     ? "complete_general_scalar_data_set"
                   : coefficient_identification
                     ? "complete_parameter_data_set"
                   : weighted_boundary_trace
                     ? "complete_weighted_boundary_data_set"
                     : neumann_convection
                       ? "complete_neumann_convection_data_set"
                     : "complete_volume_data_set",
                   normalized_laplacian
                     ? "Declare exactly one forcing, target, and regularisation datum; the normalized Laplacian has no coefficient binding."
                   : general_scalar
                     ? "Declare forcing, target, tensor diffusion, both transports, reaction, Robin coefficient/source, and regularisation data exactly once."
                   : coefficient_identification
                     ? "Declare one forcing, target, reaction, and parameter-regularisation datum, with no constant diffusion datum."
                   : weighted_boundary_trace
                     ? "Declare one forcing, target, boundary weight, diffusion, reaction, and regularisation datum."
                     : neumann_convection
                       ? "Declare one forcing, target, scalar diffusion/reaction, conservative transport, and regularisation datum."
                     : "Declare one forcing, target, diffusion, reaction, and regularisation datum.");

      const auto count_losses = [&specification](const LossKind kind) {
        return std::count_if(
          specification.losses.begin(),
          specification.losses.end(),
          [kind](const semantic::v1::LossSpec &loss) {
            return loss.kind == kind;
          });
      };
      const bool h1_control_regularisation =
        request.uses_h1_control_regularisation_loss;
      const bool hhalf_control_regularisation =
        request.uses_hhalf_control_regularisation_loss;
      const bool complete_control_loss = coefficient_identification
        ? count_losses(LossKind::quadratic_parameter_regularisation) == 1 &&
            count_losses(LossKind::quadratic_control_regularisation) == 0 &&
            count_losses(
              LossKind::quadratic_hhalf_control_regularisation) == 0 &&
            count_losses(LossKind::quadratic_h1_control_regularisation) == 0
        : hhalf_control_regularisation
        ? count_losses(
            LossKind::quadratic_hhalf_control_regularisation) == 1 &&
            count_losses(LossKind::quadratic_control_regularisation) == 0 &&
            count_losses(LossKind::quadratic_h1_control_regularisation) == 0
        : h1_control_regularisation
        ? count_losses(LossKind::quadratic_h1_control_regularisation) == 1 &&
            count_losses(LossKind::quadratic_control_regularisation) == 0 &&
            count_losses(
              LossKind::quadratic_hhalf_control_regularisation) == 0
        : count_losses(LossKind::quadratic_control_regularisation) == 1 &&
            count_losses(
              LossKind::quadratic_hhalf_control_regularisation) == 0 &&
            count_losses(LossKind::quadratic_h1_control_regularisation) == 0 &&
            count_losses(LossKind::quadratic_parameter_regularisation) == 0;
      if (count_losses(LossKind::quadratic_tracking) != 1 ||
          !complete_control_loss || specification.losses.size() != 2)
        report.add(DiagnosticCategory::lowerability,
                   specification.id,
                   "complete_registered_loss_set",
                   coefficient_identification
                     ? "Declare exactly one tracking and one parameter-regularisation loss."
                     : hhalf_control_regularisation
                     ? "Declare exactly one tracking and one H1/2 control-regularisation loss."
                     : h1_control_regularisation
                     ? "Declare exactly one tracking and one H1 control-regularisation loss."
                     : "Declare exactly one tracking and one L2 control-regularisation loss.");

      const auto selected_metric = std::find_if(
        specification.metrics.begin(),
        specification.metrics.end(),
        [&specification](const semantic::v1::MetricSpec &metric) {
          return metric.id == specification.formulation.metric_id;
        });
      if (specification.metrics.size() != 1 ||
          selected_metric == specification.metrics.end() ||
          (selected_metric->kind != semantic::v1::MetricKind::l2 &&
           selected_metric->kind != semantic::v1::MetricKind::hhalf &&
           selected_metric->kind != semantic::v1::MetricKind::h1 &&
           selected_metric->kind != semantic::v1::MetricKind::hminus1))
        report.add(DiagnosticCategory::lowerability,
                   specification.formulation.metric_id,
                   "selected_registered_metric",
                   "Select exactly one registered L2, H1/2, H1, or H-1 control metric.");
      if (coefficient_identification &&
          selected_metric != specification.metrics.end() &&
          selected_metric->kind != semantic::v1::MetricKind::l2)
        report.add(DiagnosticCategory::lowerability,
                   selected_metric->id,
                   "parameter_l2_metric",
                   "Select the registered cellwise L2 metric for the coefficient parameter.");
      if (coefficient_identification)
        {
          const auto parameter = find_variable(
            specification, specification.formulation.control_variable_id);
          const auto parameter_observation = std::find_if(
            specification.observations.begin(),
            specification.observations.end(),
            [parameter](const semantic::v1::ObservationSpec &observation) {
              return parameter != nullptr &&
                     observation.kind ==
                       semantic::v1::ObservationKind::volume_restriction &&
                     observation.input_variable_id == parameter->id;
            });
          const auto parameter_loss = std::find_if(
            specification.losses.begin(),
            specification.losses.end(),
            [](const semantic::v1::LossSpec &loss) {
              return loss.kind ==
                     semantic::v1::LossKind::quadratic_parameter_regularisation;
            });
          if (parameter == nullptr ||
              selected_metric == specification.metrics.end() ||
              selected_metric->variable_id != parameter->id ||
              parameter_observation == specification.observations.end() ||
              parameter_loss == specification.losses.end() ||
              parameter_loss->source_observation_id != parameter_observation->id)
            report.add(
              DiagnosticCategory::lowerability,
              specification.id,
              "registered_parameter_observation_metric_loss",
              "Connect the registered parameter observation, L2 metric, and parameter regularisation to the parameter decision variable.");
        }

      const auto state = find_variable(
        specification, specification.formulation.state_variable_id);
      if (boundary_control)
        {
          const auto state_observation = std::find_if(
            specification.observations.begin(),
            specification.observations.end(),
            [weighted_boundary_trace, neumann_convection](
              const semantic::v1::ObservationSpec &observation) {
              return observation.kind ==
                     (neumann_convection
                        ? semantic::v1::ObservationKind::volume_restriction
                        : weighted_boundary_trace
                        ? semantic::v1::ObservationKind::weighted_boundary_trace
                        : semantic::v1::ObservationKind::boundary_trace);
            });
          const auto control_restriction = std::find_if(
            specification.observations.begin(),
            specification.observations.end(),
            [](const semantic::v1::ObservationSpec &observation) {
              return observation.kind ==
                     semantic::v1::ObservationKind::boundary_restriction;
            });
          if (specification.observations.size() != 2 ||
              state_observation == specification.observations.end() ||
              control_restriction == specification.observations.end())
            report.add(DiagnosticCategory::lowerability,
                       specification.id,
                       neumann_convection
                         ? "complete_neumann_convection_observation_set"
                         : "complete_boundary_observation_set",
                       neumann_convection
                         ? "Declare one material-subdomain state restriction and one boundary control restriction."
                         : weighted_boundary_trace
                         ? "Declare one weighted boundary state trace and one boundary control restriction."
                         : "Declare one boundary state trace and one boundary control restriction.");
          if (state_observation != specification.observations.end())
            {
              const auto region = find_region(specification,
                                              state_observation->region_id);
              if (region == nullptr ||
                  (neumann_convection
                     ? region->kind != semantic::v1::RegionKind::volume ||
                         region->is_full_domain || region->material_ids.empty()
                     : region->kind != semantic::v1::RegionKind::boundary ||
                         region->boundary_ids.empty()))
                report.add(DiagnosticCategory::lowerability,
                           state_observation->id,
                           neumann_convection
                             ? "neumann_convection_observation_region"
                             : "boundary_trace_observation_region",
                           neumann_convection
                             ? "Select material ids for the C5.6 state observation."
                             : "Select marked boundary ids for the state trace observation.");
            }
          if (control_restriction != specification.observations.end())
            {
              const auto region = find_region(specification,
                                              control_restriction->region_id);
              const auto *control_region = request.control_boundary_region_id.empty()
                                             ? nullptr
                                             : find_region(
                                                 specification,
                                                 request.control_boundary_region_id);
              if (region == nullptr || control_region == nullptr ||
                  region->id != control_region->id)
                report.add(DiagnosticCategory::lowerability,
                           control_restriction->id,
                           "boundary_control_observation_region",
                           "Restrict the facewise control on its Neumann control boundary.");
            }
        }
      else if (dirichlet_control)
        {
          const auto state_observation = std::find_if(
            specification.observations.begin(),
            specification.observations.end(),
            [state](const semantic::v1::ObservationSpec &observation) {
              return state != nullptr &&
                     (observation.kind ==
                        semantic::v1::ObservationKind::volume_restriction ||
                      observation.kind == semantic::v1::ObservationKind::
                                            h1_state_restriction) &&
                     observation.input_variable_id == state->id;
            });
          const auto control_restriction = std::find_if(
            specification.observations.begin(),
            specification.observations.end(),
            [](const semantic::v1::ObservationSpec &observation) {
              return observation.kind ==
                     semantic::v1::ObservationKind::boundary_restriction;
            });
          const auto *control_region = request.control_boundary_region_id.empty()
                                         ? nullptr
                                         : find_region(
                                             specification,
                                             request.control_boundary_region_id);
          if (specification.observations.size() != 2 ||
              state_observation == specification.observations.end() ||
              control_restriction == specification.observations.end())
            report.add(
              DiagnosticCategory::lowerability,
              specification.id,
              "complete_dirichlet_control_observation_set",
              "Declare one full-volume physical-state restriction and one controlled-boundary control restriction.");
          if (state_observation != specification.observations.end())
            {
              const auto region = find_region(specification,
                                              state_observation->region_id);
              if (region == nullptr || !region->is_full_domain)
                report.add(
                  DiagnosticCategory::lowerability,
                  state_observation->id,
                  "dirichlet_control_full_volume_observation",
                  "Track the physical state on the full volume in the first Dirichlet-control target.");
            }
          if (control_restriction != specification.observations.end())
            {
              const auto region = find_region(specification,
                                              control_restriction->region_id);
              if (region == nullptr || control_region == nullptr ||
                  region->id != control_region->id)
                report.add(
                  DiagnosticCategory::lowerability,
                  control_restriction->id,
                  "dirichlet_control_observation_region",
                  "Restrict the nodal trace control on its declared controlled boundary.");
            }
        }
      else if (normal_flux)
        {
          const auto state_observation = std::find_if(
            specification.observations.begin(),
            specification.observations.end(),
            [state](const semantic::v1::ObservationSpec &observation) {
              return state != nullptr &&
                     observation.kind ==
                       semantic::v1::ObservationKind::normal_flux &&
                     observation.input_variable_id == state->id;
            });
          const auto control_observation = std::find_if(
            specification.observations.begin(),
            specification.observations.end(),
            [&specification](const semantic::v1::ObservationSpec &observation) {
              return observation.kind ==
                       semantic::v1::ObservationKind::volume_restriction &&
                     observation.input_variable_id ==
                       specification.formulation.control_variable_id;
            });
          if (specification.observations.size() != 2 ||
              state_observation == specification.observations.end() ||
              control_observation == specification.observations.end())
            report.add(
              DiagnosticCategory::lowerability,
              specification.id,
              "complete_normal_flux_observation_set",
              "Declare exactly one normal-flux state observation and one full-volume control observation for the first C5.8 target.");
          if (state_observation != specification.observations.end())
            {
              const auto region = find_region(specification,
                                              state_observation->region_id);
              if (region == nullptr ||
                  region->kind != semantic::v1::RegionKind::boundary ||
                  region->boundary_ids.empty())
                report.add(
                  DiagnosticCategory::lowerability,
                  state_observation->id,
                  "normal_flux_observation_region",
                  "Select a non-empty declared boundary subset for the C5.8 normal-flux state observation.");
            }
          if (!request.has_normal_flux_orientation_policy)
            report.add(
              DiagnosticCategory::lowerability,
              state_observation == specification.observations.end()
                ? specification.id
                : state_observation->id,
              "normal_flux_orientation_policy",
              "Select the registered outward-unit-normal normal-flux convention.");
          if (!request.has_normal_flux_evaluation_policy)
            report.add(
              DiagnosticCategory::lowerability,
              state_observation == specification.observations.end()
                ? specification.id
                : state_observation->id,
              "normal_flux_evaluation_policy",
              "Select the registered FE_Q face-quadrature normal-flux evaluation and transpose policy.");
          if (!request.transposition_selection ||
              request.transposition_selection->discrete_realisation !=
                semantic::v1::TranspositionDiscreteRealisation::
                  fe_q_normal_flux_very_weak)
            report.add(
              DiagnosticCategory::lowerability,
              specification.formulation.equation_id,
              "normal_flux_transposition_policy",
              "Select the registered strong-state normal-flux transposition and very-weak adjoint policy.");
          if (control_observation != specification.observations.end())
            {
              const auto region = find_region(specification,
                                              control_observation->region_id);
              if (region == nullptr || !region->is_full_domain)
                report.add(
                  DiagnosticCategory::lowerability,
                  control_observation->id,
                  "normal_flux_control_observation_region",
                  "Keep the first normal-flux target's control regularisation on the full volume.");
            }
        }
      else if (point_sensor)
        {
          const auto state_observation = std::find_if(
            specification.observations.begin(),
            specification.observations.end(),
            [state](const semantic::v1::ObservationSpec &observation) {
              return state != nullptr &&
                     observation.kind ==
                       semantic::v1::ObservationKind::point_sensor &&
                     observation.input_variable_id == state->id;
            });
          const auto control_observation = std::find_if(
            specification.observations.begin(),
            specification.observations.end(),
            [&specification](const semantic::v1::ObservationSpec &observation) {
              return observation.kind ==
                     semantic::v1::ObservationKind::volume_restriction &&
                     observation.input_variable_id ==
                       specification.formulation.control_variable_id;
            });
          if (specification.observations.size() != 2 ||
              state_observation == specification.observations.end() ||
              control_observation == specification.observations.end())
            report.add(
              DiagnosticCategory::lowerability,
              specification.id,
              "complete_point_sensor_observation_set",
              "Declare exactly one point-sensor state observation and one full-volume control observation for the first C5.10 target.");
          if (state_observation != specification.observations.end())
            {
              const auto region = find_region(specification,
                                              state_observation->region_id);
              if (region == nullptr ||
                  region->kind != semantic::v1::RegionKind::point_set ||
                  region->point_coordinates.empty())
                report.add(
                  DiagnosticCategory::lowerability,
                  state_observation->id,
                  "point_sensor_observation_region",
                  "Select a non-empty immutable point-set region for the C5.10 state observation.");
            }
          if (!request.has_point_sensor_evaluation_policy)
            report.add(
              DiagnosticCategory::lowerability,
              state_observation == specification.observations.end()
                ? specification.id
                : state_observation->id,
              "point_sensor_evaluation_policy",
              "Select the registered FE_Q physical-point evaluation and assembled transpose policy; nearest-node and quadrature-coincidence rules are not registered.");
          if (!request.transposition_selection ||
              request.transposition_selection->discrete_realisation !=
                semantic::v1::TranspositionDiscreteRealisation::
                  fe_q_point_sensor_very_weak)
            report.add(
              DiagnosticCategory::lowerability,
              specification.formulation.equation_id,
              "point_sensor_transposition_policy",
              "Select the registered very-weak point-sensor transpose formulation.");
          if (control_observation != specification.observations.end())
            {
              const auto region = find_region(specification,
                                              control_observation->region_id);
              if (region == nullptr || !region->is_full_domain)
                report.add(
                  DiagnosticCategory::lowerability,
                  control_observation->id,
                  "point_sensor_control_observation_region",
                  "Keep the first point-sensor target's control regularisation on the full volume.");
            }
        }
      else
        for (const auto &observation : specification.observations)
          {
            const auto region = find_region(specification, observation.region_id);
            if (region == nullptr || region->kind != semantic::v1::RegionKind::volume)
              report.add(DiagnosticCategory::lowerability,
                         observation.id,
                         "volume_observation_region",
                         "Select a registered volume observation region.");
            else if (state != nullptr &&
                     observation.input_variable_id == state->id &&
                     !region->is_full_domain && region->material_ids.empty())
              report.add(DiagnosticCategory::lowerability,
                         observation.id,
                         "material_subdomain_observation",
                         "Declare one or more material ids for the subdomain observation.");
            else if ((state == nullptr ||
                      observation.input_variable_id != state->id) &&
                     !region->is_full_domain)
              report.add(DiagnosticCategory::lowerability,
                         observation.id,
                         "full_domain_nonstate_observation",
                         "Only the state tracking observation supports material subdomains in v1.");
          }

      const bool has_constraint = !specification.formulation.constraint_id.empty();
      if (!has_constraint && !specification.constraints.empty())
        report.add(DiagnosticCategory::lowerability,
                   specification.id,
                   "selected_constraint_port",
                   "Select the declared constraint in the reduced formulation or remove it.");
      if (has_constraint &&
          (specification.constraints.size() != 1 ||
           count_data(DataRole::lower_bound) != 1 ||
           count_data(DataRole::upper_bound) != 1))
        report.add(DiagnosticCategory::lowerability,
                   specification.formulation.constraint_id,
                   boundary_control ? "complete_facewise_box_data_set"
                                    : "complete_cellwise_box_data_set",
                   boundary_control
                     ? "Declare one selected facewise box plus one lower and one upper bound datum."
                     : "Declare one selected box plus one lower and one upper bound datum.");
      if (has_constraint)
        {
          const auto constraint = std::find_if(
            specification.constraints.begin(),
            specification.constraints.end(),
            [&specification](const semantic::v1::ConstraintSpec &candidate) {
              return candidate.id == specification.formulation.constraint_id;
            });
          if (constraint != specification.constraints.end() &&
              constraint->kind != (boundary_control
                                    ? semantic::v1::ConstraintKind::facewise_box
                                    : semantic::v1::ConstraintKind::cellwise_box))
            report.add(DiagnosticCategory::lowerability,
                       constraint->id,
                       "control_layout_constraint_realisation",
                       "Use a facewise box for boundary control and a cellwise box for volume control.");
        }
    }

    static void
    validate_formulation_capability(
      const semantic::v1::ProblemSpec &specification,
      semantic::v1::ValidationReport & report)
    {
      using semantic::v1::DiagnosticCategory;
      const bool supplied_otd =
        specification.formulation.kind ==
          semantic::v1::FormulationKind::all_at_once &&
        specification.formulation.provenance ==
          semantic::v1::FormulationProvenance::supplied_otd;
      if (!supplied_otd &&
          (specification.formulation.kind !=
               semantic::v1::FormulationKind::reduced_dto ||
           specification.formulation.provenance !=
             semantic::v1::FormulationProvenance::dto))
        report.add(DiagnosticCategory::formulation_capability,
                   specification.formulation.id,
                   "reduced_dto_formulation",
                   "Select the registered v1 DTO formulation; the requested provenance or execution shape is not available.");
      if (specification.variables.size() != 2 ||
          specification.equations.size() != 1)
        report.add(DiagnosticCategory::formulation_capability,
                   specification.formulation.id,
                   "one_state_one_decision_one_equation",
                   "The executable DTO contract currently supports one state, one control-or-parameter decision, and one equation block.");
      if (specification.formulation.constraint_id.empty())
        return;
      const auto constraint = std::find_if(
        specification.constraints.begin(),
        specification.constraints.end(),
        [&specification](const semantic::v1::ConstraintSpec &candidate) {
          return candidate.id == specification.formulation.constraint_id;
        });
      if (constraint == specification.constraints.end() ||
          (constraint->kind != semantic::v1::ConstraintKind::cellwise_box &&
           constraint->kind != semantic::v1::ConstraintKind::facewise_box))
        report.add(DiagnosticCategory::formulation_capability,
                   specification.formulation.id,
                   "l2_coefficientwise_projected_gradient",
                   "Use the registered cellwise or facewise L2 box constraint, or omit the constraint.");
    }

    static void
    validate_product_capability(
      const semantic::v1::ProblemSpec &       specification,
      const ResolvedCompilationRequest &     request,
      const DealiiDiscretisationPolicy &      policy,
      const CompilationProduct                product,
      semantic::v1::ValidationReport &        report)
    {
      if (product == CompilationProduct::pdas)
        {
          const bool reduced_dto =
            specification.formulation.kind ==
              semantic::v1::FormulationKind::reduced_dto &&
            specification.formulation.provenance ==
              semantic::v1::FormulationProvenance::dto;
          const bool supplied_otd =
            specification.formulation.kind ==
              semantic::v1::FormulationKind::all_at_once &&
            specification.formulation.provenance ==
              semantic::v1::FormulationProvenance::supplied_otd;
          const auto constraint = std::find_if(
            specification.constraints.begin(),
            specification.constraints.end(),
            [&specification](const semantic::v1::ConstraintSpec &candidate) {
              return candidate.id == specification.formulation.constraint_id;
            });
          const bool cellwise_box =
            constraint != specification.constraints.end() &&
            constraint->kind == semantic::v1::ConstraintKind::cellwise_box;
          if (!(request.target_family == ResolvedTargetFamily::direct_volume &&
                (reduced_dto || supplied_otd) && cellwise_box))
            report.add(
              semantic::v1::DiagnosticCategory::formulation_capability,
              specification.formulation.id,
              "compiled_pdas_cellwise_box",
              "Request the registered scalar distributed-control DTO or supplied-OTD target with a cellwise L2 box constraint.");

          if (request.uses_h1_control_metric ||
              request.uses_hhalf_control_metric ||
              request.uses_hminus1_control_metric)
            report.add(
              semantic::v1::DiagnosticCategory::formulation_capability,
              specification.formulation.metric_id,
              "compiled_pdas_metric_realisation",
              "Select the registered positive-diagonal cellwise L2 metric for the PDAS multiplier conversion.");

          if (!contract::valid(policy.pdas))
            report.add(
              semantic::v1::DiagnosticCategory::formulation_capability,
              specification.formulation.id,
              "compiled_pdas_policy",
              "Supply finite positive PDAS iteration and residual tolerances.");
          const auto &assumptions = policy.pdas.active_set_assumptions;
          if (!assumptions.rank_condition_declared ||
              assumptions.rank_policy.empty())
            report.add(
              semantic::v1::DiagnosticCategory::formulation_capability,
              specification.formulation.id,
              "compiled_pdas_active_row_rank",
              "Declare the active-row rank condition and its policy before constructing the PDAS product.");
          if (!assumptions.kernel_positivity_declared ||
              assumptions.kernel_policy.empty())
            report.add(
              semantic::v1::DiagnosticCategory::formulation_capability,
              specification.formulation.id,
              "compiled_pdas_kernel_positivity",
              "Declare positive definiteness on the restricted active-set kernel before constructing the PDAS product.");
          if (!contract::valid(policy.pdas_kkt_solver) ||
              policy.pdas_kkt_solver.maximum_iterations == 0)
            report.add(
              semantic::v1::DiagnosticCategory::formulation_capability,
              specification.formulation.id,
              "compiled_pdas_inner_kkt_solver",
              "Supply a valid inner quadratic KKT solver policy with a positive iteration limit.");
          return;
        }

      if (product != CompilationProduct::quadratic_kkt)
        return;

      const bool canonical_scalar_dto =
        specification.id == "scalar_diffusion_reaction_fixed_dirichlet" &&
        request.uses_assembled_v1_target &&
        specification.formulation.kind ==
          semantic::v1::FormulationKind::reduced_dto &&
        specification.formulation.provenance ==
          semantic::v1::FormulationProvenance::dto &&
        specification.formulation.constraint_id.empty();
      if (!canonical_scalar_dto)
        report.add(
          semantic::v1::DiagnosticCategory::formulation_capability,
          specification.formulation.id,
          "compiled_quadratic_kkt",
          "Request the canonical unconstrained scalar DTO target before constructing a compiled quadratic KKT product.");
    }

    static void
    validate_supplied_otd_capability(
      const semantic::v1::ProblemSpec & specification,
      const ResolvedCompilationRequest &request,
      const CompilationProduct            product,
      semantic::v1::ValidationReport &   report)
    {
      using semantic::v1::DiagnosticCategory;
      const bool supplied_otd =
        specification.formulation.kind ==
          semantic::v1::FormulationKind::all_at_once &&
        specification.formulation.provenance ==
          semantic::v1::FormulationProvenance::supplied_otd;
      if (!supplied_otd)
        return;

      if (!specification.supplied_otd_declaration)
        {
          report.add(
            DiagnosticCategory::formulation_capability,
            specification.formulation.id,
            "supplied_otd_declaration",
            "Select a complete registered supplied-OTD declaration; changing DTO formulation labels is not sufficient.");
          return;
        }
      const auto &declaration = *specification.supplied_otd_declaration;

      if (request.target_family != ResolvedTargetFamily::direct_volume)
        report.add(
          DiagnosticCategory::formulation_capability,
          specification.formulation.id,
          "supplied_otd_scalar_target",
          "The first supplied-OTD lowerer supports only the canonical scalar diffusion-reaction target.");
      if (!specification.formulation.constraint_id.empty() &&
          product != CompilationProduct::pdas)
        report.add(
          DiagnosticCategory::formulation_capability,
          specification.formulation.constraint_id,
          "supplied_otd_equality_only",
          "The supplied-OTD execution product currently supports equality optimality blocks without a box constraint.");

      if (declaration.state_block.variable_space_id != "state_space" ||
          declaration.state_block.residual_space_id != "state_test_space" ||
          declaration.state_block.trial_pairing_id != "state_pairing" ||
          declaration.state_block.test_pairing_id != "state_test_pairing")
        report.add(
          DiagnosticCategory::formulation_capability,
          declaration.state_block.id,
          "supplied_otd_state_block",
          "The canonical supplied-OTD lowerer requires the declared state space, state-test space, and pairings.");

      if (declaration.adjoint_block.variable_space_id != "state_test_space" ||
          declaration.adjoint_block.residual_space_id != "state_test_space" ||
          declaration.adjoint_block.trial_pairing_id != "state_test_pairing" ||
          declaration.adjoint_block.test_pairing_id != "state_test_pairing")
        report.add(
          DiagnosticCategory::formulation_capability,
          declaration.adjoint_block.id,
          "supplied_otd_adjoint_space",
          "The canonical supplied-OTD lowerer requires the declared state-test adjoint space and pairing.");

      if (declaration.control_stationarity_block.variable_space_id !=
            "control_space" ||
          declaration.control_stationarity_block.residual_space_id !=
            "control_space" ||
          declaration.control_stationarity_block.trial_pairing_id !=
            "control_pairing" ||
          declaration.control_stationarity_block.test_pairing_id !=
            "control_pairing")
        report.add(
          DiagnosticCategory::formulation_capability,
          declaration.control_stationarity_block.id,
          "supplied_otd_stationarity_space",
          "The canonical supplied-OTD lowerer requires the declared control space and control pairing.");

      if (declaration.multiplier_convention !=
            semantic::v1::SuppliedOTDMultiplierConvention::framework_adjoint ||
          declaration.multiplier_conversion !=
            semantic::v1::SuppliedOTDMultiplierConversion::identity)
        report.add(
          DiagnosticCategory::formulation_capability,
          declaration.id,
          "supplied_otd_multiplier_conversion",
          "The canonical supplied-OTD lowerer requires a framework-adjoint multiplier with identity conversion.");
    }

    template <typename Component>
    static std::vector<std::string>
    identifiers(const std::vector<Component> &components)
    {
      std::vector<std::string> ids;
      ids.reserve(components.size());
      for (const auto &component : components)
        ids.push_back(component.id);
      std::sort(ids.begin(), ids.end());
      return ids;
    }

    static std::string
    describe(const ConstraintRealisation realisation)
    {
      switch (realisation)
        {
          case ConstraintRealisation::none:
            return "none";
          case ConstraintRealisation::cellwise_l2:
            return "FE_DGQ(0) coefficientwise l2_cellwise clipping";
          case ConstraintRealisation::cellwise_parameter_l2:
            return "FE_DGQ(0) coefficientwise l2_cellwise_parameter clipping";
          case ConstraintRealisation::facewise_l2:
            return "facewise-constant coefficientwise l2_facewise clipping";
        }
      contract::require(false, "Unknown compiled constraint realization");
      return {};
    }

    static std::string
    constraint_realisation_id(const ConstraintRealisation realisation)
    {
      switch (realisation)
        {
          case ConstraintRealisation::none:
            return "none";
          case ConstraintRealisation::cellwise_l2:
            return "l2_cellwise";
          case ConstraintRealisation::cellwise_parameter_l2:
            return "l2_cellwise_parameter";
          case ConstraintRealisation::facewise_l2:
            return "l2_facewise";
        }
      contract::require(false, "Unknown compiled constraint realization");
      return {};
    }

    static std::string
    describe(const CompiledConstraintRecord &record)
    {
      if (!record.present || record.realisation_id == "none")
        return "none";
      if (record.realisation_id == "l2_cellwise")
        return "FE_DGQ(0) coefficientwise l2_cellwise clipping";
      if (record.realisation_id == "l2_cellwise_parameter")
        return "FE_DGQ(0) coefficientwise l2_cellwise_parameter clipping";
      if (record.realisation_id == "l2_facewise")
        return "facewise-constant coefficientwise l2_facewise clipping";
      contract::require(false, "Unknown compiled constraint record realization");
      return {};
    }

    static std::string
    metric_solve_policy_description(const CompiledMetricRecord &record)
    {
      const auto &policy = record.solve_policy;
      return record.realisation_id + ": " + record.operator_description +
             "; serial CG with " + policy.preconditioner +
             " preconditioner, maximum iterations=" +
             std::to_string(policy.maximum_iterations) +
             ", relative tolerance=" +
             std::to_string(policy.relative_tolerance) +
             ", absolute tolerance=" +
             std::to_string(policy.absolute_tolerance);
    }

    static std::string
    state_adjoint_solve_policy_description(
      const CompiledSolvePolicyRecord &state,
      const CompiledSolvePolicyRecord &adjoint)
    {
      const auto describe_policy = [](const CompiledSolvePolicyRecord &record) {
        const std::string algorithm =
          record.algorithm == LinearSolveAlgorithm::serial_sparse_direct_umfpack
            ? "serial SparseDirectUMFPACK"
            : "serial CG with " + record.preconditioner + " preconditioner";
        return algorithm + " (maximum iterations=" +
               std::to_string(record.maximum_iterations) +
               ", relative tolerance=" +
               std::to_string(record.relative_tolerance) +
               ", absolute tolerance=" +
               std::to_string(record.absolute_tolerance) +
               (record.operator_realisation.empty()
                  ? std::string{}
                  : "; " + record.operator_realisation) + ")";
      };
      return "state=" + describe_policy(state) +
             "; adjoint=" + describe_policy(adjoint);
    }

    static const dealii_backend::MassMetric &
    as_mass_metric(const contract::MetricT<dealii_backend::SerialBackend> &metric)
    {
      const auto *mass_metric =
        dynamic_cast<const dealii_backend::MassMetric *>(&metric);
      contract::require(
        mass_metric != nullptr,
        "A coefficientwise constraint requires its concrete mass metric");
      return *mass_metric;
    }

    static dealii::Vector<double>
    make_pdas_bound_vector(const CellwiseBoundValue &value,
                           const std::size_t           dimension,
                           const char *                description)
    {
      if (std::holds_alternative<double>(value))
        {
          dealii::Vector<double> result(
            dealii_backend::SerialBackend::checked_native_size(dimension));
          result = std::get<double>(value);
          return result;
        }

      const auto &bound = std::get<dealii::Vector<double>>(value);
      contract::require(
        static_cast<std::size_t>(bound.size()) == dimension,
        std::string("PDAS ") + description + " bound has the wrong layout");
      return bound;
    }

    static std::shared_ptr<const contract::BoxComplementarityT<
      dealii_backend::SerialBackend>>
    make_pdas_complementarity(
      const contract::EqualityConstrainedQuadraticKKTProductT<
        dealii_backend::SerialBackend> &product,
      const CompiledCellwiseBoxDataT<dealii_backend::SerialBackend> &box_data,
      std::shared_ptr<const contract::MetricT<dealii_backend::SerialBackend>>
        metric)
    {
      using Backend = dealii_backend::SerialBackend;
      using Primal = contract::PrimalBlockT<Backend>;
      using Complementarity = contract::BoxComplementarityT<Backend>;
      contract::require(static_cast<bool>(metric),
                        "Compiled PDAS needs an owned metric");
      contract::require(product.layout().primal->n_blocks() == 2,
                        "Compiled PDAS needs state and control KKT blocks");
      const auto control_layout =
        product.layout().primal->single_block(1, "compiled_pdas_control");
      contract::require(metric->layout()->compatible_with(*control_layout),
                        "Compiled PDAS metric does not match its control block");
      contract::require(box_data.layout()->compatible_with(*control_layout),
                        "Compiled PDAS shared box does not match its control block");
      contract::require(box_data.metric_owner()->realisation_witness().matches(
                          metric->realisation_witness()),
                        "Compiled PDAS shared box does not match its metric");
      const auto &mass_metric = as_mass_metric(*metric);
      contract::require(
        mass_metric.supports_coefficientwise_box_projection(),
        "Compiled PDAS needs a positive diagonal L2 metric realization");
      dealii::Vector<double> lower = box_data.lower().block(0);
      dealii::Vector<double> upper = box_data.upper().block(0);
      return std::make_shared<const Complementarity>(
        contract::BoxBoundsT<Backend>(
          control_layout,
          Primal(control_layout, {lower}),
          Primal(control_layout, {upper})),
        contract::make_metric_multiplier_representation(
          box_data.metric_owner()),
        box_data.token());
    }

    static bool
    uses_neumann_target(const CompiledTargetKind target)
    {
      return target == CompiledTargetKind::neumann_boundary ||
             target == CompiledTargetKind::weighted_boundary_trace ||
             target == CompiledTargetKind::pure_neumann;
    }

    static std::string
    control_space_description(const CompiledTargetKind       target,
                              const DealiiDiscretisationPolicy &policy,
                              const std::optional<DirichletControlRegistration> &
                                registration = std::nullopt)
    {
      switch (target)
        {
          case CompiledTargetKind::dirichlet_control:
            return registration.has_value() &&
                     *registration ==
                       DirichletControlRegistration::partial_nodal_l2
              ? "relative-interior nodal trace coefficients on the partial controlled boundary with fixed endpoint precedence"
              : "one shared nodal trace coefficient per state DoF on the complete controlled exterior boundary";
          case CompiledTargetKind::l2_dirichlet_transposition:
            return "conforming nodal trace FE subspace U_h=trace(V_h) of the continuous L2 boundary control";
          case CompiledTargetKind::hhalf_dirichlet_control:
          case CompiledTargetKind::h1_tracking_hhalf_dirichlet_control:
            return "conforming nodal trace FE subspace U_h=trace(V_h) with the minimum-extension H1/2 geometry";
          case CompiledTargetKind::h1_dirichlet_control:
            return "conforming nodal trace FE subspace with boundary mass-plus-tangential-stiffness H1 geometry";
          case CompiledTargetKind::neumann_boundary:
          case CompiledTargetKind::weighted_boundary_trace:
          case CompiledTargetKind::pure_neumann:
            return "one facewise-constant coefficient per marked state boundary face";
          case CompiledTargetKind::coefficient_identification:
            return "cellwise-constant positive diffusion parameter FE_DGQ(0) on the state mesh";
          case CompiledTargetKind::h1_control_l2_metric:
          case CompiledTargetKind::h1_control_h1_metric:
            return "continuous scalar FE_Q(" +
                   std::to_string(policy.state_degree) + ") on the state mesh";
          case CompiledTargetKind::hminus1_control_metric:
          case CompiledTargetKind::continuous_control_l2_metric:
            return "independent homogeneous-Dirichlet scalar FE_Q(" +
                   std::to_string(policy.state_degree) +
                   ") coefficients on the state mesh";
          case CompiledTargetKind::direct_volume:
          case CompiledTargetKind::assembled_volume:
          case CompiledTargetKind::general_scalar_robin:
          case CompiledTargetKind::point_sensor:
          case CompiledTargetKind::normal_flux:
            return "FE_DGQ(0) on the state active-cell mesh";
        }
      contract::require(false, "Unknown compiled target kind");
      return {};
    }

    static unsigned int
    resolved_maximum_iterations(
      const dealii_backend::SPDLinearSolvePolicy &policy,
      const std::size_t                           dimension)
    {
      if (policy.maximum_iterations != 0)
        return policy.maximum_iterations;
      contract::require(
        dimension <= std::numeric_limits<unsigned int>::max() / 10U,
        "Compiled solve dimension exceeds the iteration-policy range");
      return std::max(100U, 10U * static_cast<unsigned int>(dimension));
    }

    static CompiledSolvePolicyRecord
    spd_solve_record(const dealii_backend::SPDLinearSolvePolicy &policy,
                     const std::size_t                           dimension,
                     std::string                                 nullspace_policy)
    {
      return {LinearSolveAlgorithm::serial_cg,
              "identity",
              resolved_maximum_iterations(policy, dimension),
              policy.relative_tolerance,
              policy.absolute_tolerance,
              std::move(nullspace_policy),
              {}};
    }

    template <int dim>
    static std::size_t
    physical_state_dimension(
      const contract::ExecutableModelT<dealii_backend::SerialBackend> &executable)
    {
      if (const auto *scalar =
            dynamic_cast<const detail::ScalarComponentModel<dim> *>(&executable))
        return scalar->physical_state_dimension();
      if (const auto *neumann =
            dynamic_cast<const detail::NeumannBoundaryControlModel<dim> *>(
              &executable))
        return neumann->physical_state_dimension();
      if (const auto *dirichlet =
            dynamic_cast<const detail::DirichletControlLiftingModel<dim> *>(
              &executable))
        return dirichlet->physical_state_dimension();
      return executable.variable_layout()->dimension(0);
    }

    template <int dim>
    static std::size_t
    realized_observation_dimension(
      const contract::ExecutableModelT<dealii_backend::SerialBackend> &executable)
    {
      if (const auto *scalar =
            dynamic_cast<const detail::ScalarComponentModel<dim> *>(&executable))
        return scalar->realized_observation_dimension();
      if (const auto *neumann =
            dynamic_cast<const detail::NeumannBoundaryControlModel<dim> *>(
              &executable))
        return neumann->realized_observation_dimension();
      return physical_state_dimension<dim>(executable);
    }

    static std::optional<std::size_t>
    realized_output_dimension(
      const std::vector<CompiledRealizedSpaceRecord> &realized_spaces,
      const std::string &                              semantic_space_id)
    {
      const auto space = std::find_if(
        realized_spaces.begin(),
        realized_spaces.end(),
        [&semantic_space_id](const CompiledRealizedSpaceRecord &candidate) {
          return candidate.semantic_id == semantic_space_id &&
                 candidate.realization_id.find("output:") == 0;
        });
      return space == realized_spaces.end()
               ? std::nullopt
               : std::optional<std::size_t>(space->dimension);
    }

    static std::size_t
    compiled_space_dimension(
      const CompiledSpaceRecord &space,
      const contract::ExecutableModelT<dealii_backend::SerialBackend> &executable,
      const std::vector<CompiledRealizedSpaceRecord> &realized_spaces)
    {
      switch (space.role)
        {
          case semantic::v1::SpaceRole::state:
            return executable.variable_layout()->dimension(0);
          case semantic::v1::SpaceRole::test:
            return executable.test_layout()->dimension(0);
          case semantic::v1::SpaceRole::control:
          case semantic::v1::SpaceRole::parameter:
            return executable.variable_layout()->dimension(1);
          case semantic::v1::SpaceRole::observation:
            if (const auto dimension =
                  realized_output_dimension(realized_spaces,
                                            space.semantic_id))
              return *dimension;
            return 0;
          case semantic::v1::SpaceRole::data:
          case semantic::v1::SpaceRole::auxiliary:
          case semantic::v1::SpaceRole::unspecified:
            return 0;
        }
      return 0;
    }

    static std::string
    compiled_space_runtime_role(const semantic::v1::SpaceRole role)
    {
      switch (role)
        {
          case semantic::v1::SpaceRole::state:
            return "state";
          case semantic::v1::SpaceRole::test:
            return "test_and_adjoint";
          case semantic::v1::SpaceRole::control:
            return "decision_control";
          case semantic::v1::SpaceRole::parameter:
            return "decision_parameter";
          case semantic::v1::SpaceRole::observation:
            return "observation";
          case semantic::v1::SpaceRole::data:
            return "data";
          case semantic::v1::SpaceRole::auxiliary:
            return "auxiliary";
          case semantic::v1::SpaceRole::unspecified:
            return "unspecified";
        }
      return "unspecified";
    }

    static std::string
    compiled_space_finite_element(
      const semantic::v1::SpaceSpec &   space,
      const CompiledTargetKind          target,
      const DealiiDiscretisationPolicy &policy)
    {
      switch (space.role)
        {
          case semantic::v1::SpaceRole::state:
          case semantic::v1::SpaceRole::test:
            return "scalar FE_Q(" + std::to_string(policy.state_degree) + ")";
          case semantic::v1::SpaceRole::control:
          case semantic::v1::SpaceRole::parameter:
            return control_space_description(target, policy);
          case semantic::v1::SpaceRole::observation:
            return "lowered observation coefficients";
          case semantic::v1::SpaceRole::data:
            return "external binding";
          case semantic::v1::SpaceRole::auxiliary:
            return "semantic auxiliary space";
          case semantic::v1::SpaceRole::unspecified:
            return "unspecified";
        }
      return "unspecified";
    }

    static std::string
    bound_binding_description(
      const std::optional<CellwiseBoxDataBindings> &bounds,
      const std::optional<FacewiseBoxDataBindings> &facewise_bounds)
    {
      if (bounds)
        return std::holds_alternative<double>(bounds->lower)
                 ? "scalar bound"
                 : "exact-layout cellwise coefficient vector";
      if (facewise_bounds)
        return std::holds_alternative<double>(facewise_bounds->lower)
                 ? "scalar bound"
                 : "exact-layout facewise coefficient vector";
      return "unbound";
    }

    static std::string
    general_scalar_data_rule(
      const std::vector<CompiledBindingRecord> &bindings,
      const unsigned int                         quadrature_order)
    {
      const auto binding_for_role = [&bindings](
                                      const semantic::v1::DataRole role) {
        const auto binding = std::find_if(
          bindings.begin(),
          bindings.end(),
          [role](const CompiledBindingRecord &candidate) {
            return candidate.role == role;
          });
        contract::require(binding != bindings.end(),
                          "The general scalar manifest is missing a data binding");
        return &*binding;
      };
      const auto describe = [](const CompiledBindingRecord &binding) {
        return binding.representation + " [space=" + binding.space_id +
               ", region=" + binding.region_id + ", evaluation=" +
               binding.evaluation_realisation + "]";
      };
      const auto diffusion = binding_for_role(semantic::v1::DataRole::diffusion);
      const auto forcing = binding_for_role(semantic::v1::DataRole::forcing);
      const auto desired_state =
        binding_for_role(semantic::v1::DataRole::desired_state);
      const auto conservative =
        binding_for_role(semantic::v1::DataRole::conservative_transport);
      const auto advective =
        binding_for_role(semantic::v1::DataRole::advective_transport);
      const auto reaction = binding_for_role(semantic::v1::DataRole::reaction);
      const auto robin_coefficient =
        binding_for_role(semantic::v1::DataRole::robin_coefficient);
      const auto robin_source =
        binding_for_role(semantic::v1::DataRole::robin_source);
      return "general scalar data at selected QGauss(" +
             std::to_string(quadrature_order) + ") quadrature: " +
             "forcing=" + describe(*forcing) + ", desired_state=" +
             describe(*desired_state) + ", " + describe(*diffusion) + ", " +
             describe(*conservative) + ", " +
             describe(*advective) + ", " + describe(*reaction) +
             "; Robin coefficient and source: " + describe(*robin_coefficient) +
             ", " + describe(*robin_source);
    }

    static std::string
    boundary_realisation_description(
      const semantic::v1::BoundaryRealisationSelection &selection)
    {
      return "boundary selection " + selection.id + ": fixed=" +
             selection.fixed_dirichlet_region_id + ", robin=" +
             selection.robin_region_id + ", neumann=" +
             (selection.neumann_region_ids.empty() ? "empty" : "nonempty") +
             ", transport_inflow=" +
             (selection.transport_inflow_region_ids.empty() ? "empty"
                                                             : "nonempty") +
             ", transport_outflow=" +
             selection.transport_outflow_region_id +
             ", conormal=outward(A grad(y) - b y), trace=FE_Q state trace, face=QGauss";
    }

    static std::uint64_t
    hash_word(const std::uint64_t hash, const std::uint64_t word)
    {
      std::uint64_t result = hash;
      for (unsigned int byte = 0; byte < sizeof(word); ++byte)
        {
          result ^= (word >> (8U * byte)) & 0xffU;
          result *= 1099511628211ULL;
        }
      return result;
    }

    static std::uint64_t
    double_bits(const double value)
    {
      std::uint64_t bits = 0;
      static_assert(sizeof(bits) == sizeof(value));
      std::memcpy(&bits, &value, sizeof(bits));
      return bits;
    }

    static std::string
    hash_text(const std::uint64_t hash)
    {
      std::ostringstream stream;
      stream << "fnv1a64:0x" << std::hex << std::setw(16) << std::setfill('0')
             << hash;
      return stream.str();
    }

    static std::string
    vector_identity(const dealii::Vector<double> &values)
    {
      std::uint64_t hash = 1469598103934665603ULL;
      hash = hash_word(hash, values.size());
      for (std::size_t index = 0; index < values.size(); ++index)
        hash = hash_word(hash, double_bits(values[index]));
      return hash_text(hash);
    }

    template <int dim>
    static std::string
    mesh_structural_identity(const dealii::Triangulation<dim> &triangulation)
    {
      std::uint64_t hash = 1469598103934665603ULL;
      hash = hash_word(hash, dim);
      hash = hash_word(hash, triangulation.n_active_cells());
      for (const auto &cell : triangulation.active_cell_iterators())
        {
          hash = hash_word(hash, cell->material_id());
          for (unsigned int vertex = 0;
               vertex < dealii::GeometryInfo<dim>::vertices_per_cell;
               ++vertex)
            for (unsigned int coordinate = 0; coordinate < dim; ++coordinate)
              hash = hash_word(hash, double_bits(cell->vertex(vertex)[coordinate]));
          for (unsigned int face = 0;
               face < dealii::GeometryInfo<dim>::faces_per_cell;
               ++face)
            {
              hash = hash_word(hash, cell->face(face)->at_boundary() ? 1U : 0U);
              hash = hash_word(hash, cell->face(face)->boundary_id());
            }
        }
      return hash_text(hash);
    }

    static CompiledFieldShape
    compiled_field_shape(const semantic::v1::DataKind kind)
    {
      switch (kind)
        {
          case semantic::v1::DataKind::function:
            return CompiledFieldShape::scalar;
          case semantic::v1::DataKind::vector_function:
            return CompiledFieldShape::vector;
          case semantic::v1::DataKind::tensor_function:
            return CompiledFieldShape::tensor;
          case semantic::v1::DataKind::scalar_constant:
            return CompiledFieldShape::scalar_constant;
          case semantic::v1::DataKind::cellwise_bound:
            return CompiledFieldShape::cellwise_scalar;
          case semantic::v1::DataKind::facewise_bound:
            return CompiledFieldShape::facewise_scalar;
          case semantic::v1::DataKind::unspecified:
            return CompiledFieldShape::unspecified;
        }
      return CompiledFieldShape::unspecified;
    }

    template <int dim>
    static std::vector<CompiledBindingRecord>
    make_resolved_binding_records(
      const semantic::v1::ProblemSpec &              specification,
      const CompiledTargetKind                       target,
      const ResolvedCompilationRequest &             request,
      const DealiiDataBindings<dim> &                data,
      const std::optional<CellwiseBoxDataBindings> & bounds,
      const std::optional<FacewiseBoxDataBindings> & facewise_bounds)
    {
      const bool uses_coefficient_identification =
        target == CompiledTargetKind::coefficient_identification;
      const bool uses_general_scalar =
        target == CompiledTargetKind::general_scalar_robin;
      const bool uses_h1_state_observation =
        has_h1_state_observation(specification);
      const bool uses_point_sensor =
        target == CompiledTargetKind::point_sensor ||
        has_point_sensor_observation(specification);
      const bool uses_normal_flux =
        target == CompiledTargetKind::normal_flux ||
        has_normal_flux_observation(specification);

      const auto request_for = [&request](const std::string &semantic_id) {
        const auto match = std::find_if(
          request.data_bindings.begin(),
          request.data_bindings.end(),
          [&semantic_id](const ResolvedDataBindingRequest &candidate) {
            return candidate.semantic_id == semantic_id;
          });
        return match == request.data_bindings.end() ? nullptr : &*match;
      };
      const auto space_region = [&specification](const std::string &space_id) {
        const auto space = std::find_if(
          specification.spaces.begin(),
          specification.spaces.end(),
          [&space_id](const semantic::v1::SpaceSpec &candidate) {
            return candidate.id == space_id;
          });
        return space == specification.spaces.end() ? std::string{} : space->region_id;
      };

      std::vector<CompiledBindingRecord> records;
      records.reserve(specification.data.size());
      for (const auto &binding : specification.data)
        {
          CompiledBindingRecord record;
          record.semantic_id = binding.id;
          record.role = binding.role;
          record.kind = binding.kind;
          record.field_shape = compiled_field_shape(binding.kind);
          record.space_id = binding.space_id;
          record.region_id = space_region(binding.space_id);
          record.evaluation_realisation = "not applicable";
          record.runtime_representation =
            "semantic DataKind " + std::to_string(static_cast<int>(binding.kind));
          if (const auto *resolved_binding = request_for(binding.id))
            {
              record.evaluation_realisation =
                resolved_binding->evaluation_realisation;
              record.runtime_representation =
                resolved_binding->runtime_representation;
            }

          switch (binding.role)
            {
              case semantic::v1::DataRole::forcing:
                record.representation = "analytic Function at quadrature";
                record.provenance = data.provenance.forcing;
                break;
              case semantic::v1::DataRole::desired_state:
                record.representation = uses_h1_state_observation
                                          ? "analytic Function value and gradient at quadrature"
                                        : uses_point_sensor
                                          ? "analytic Function value at immutable sensor coordinates"
                                        : uses_normal_flux
                                          ? "analytic Function value at selected boundary face quadrature"
                                          : "analytic Function at quadrature";
                record.provenance = data.provenance.desired_state;
                break;
              case semantic::v1::DataRole::fixed_dirichlet_lifting:
                record.representation = "Function interpolated at boundary DoFs";
                record.provenance = data.provenance.fixed_dirichlet_data;
                break;
              case semantic::v1::DataRole::diffusion:
                record.representation = uses_coefficient_identification
                                          ? "decision parameter"
                                        : uses_general_scalar
                                          ? "tensor Function at volume quadrature"
                                          : "scalar constant";
                record.provenance = uses_coefficient_identification
                                      ? specification.formulation.control_variable_id
                                    : uses_general_scalar
                                      ? data.general_scalar->provenance.diffusion_tensor
                                    : data.diffusion
                                      ? std::to_string(*data.diffusion)
                                      : "unbound";
                break;
              case semantic::v1::DataRole::conservative_transport:
                record.representation = "vector Function at volume quadrature";
                record.provenance = uses_general_scalar
                  ? data.general_scalar->provenance.conservative_transport
                  : data.conservative_transport->provenance.conservative_transport;
                break;
              case semantic::v1::DataRole::advective_transport:
                record.representation = "vector Function at volume quadrature";
                record.provenance = data.general_scalar->provenance.advective_transport;
                break;
              case semantic::v1::DataRole::reaction:
                record.representation = uses_general_scalar
                                          ? "scalar Function at volume quadrature"
                                          : "scalar constant";
                record.provenance = uses_general_scalar
                                      ? data.general_scalar->provenance.reaction
                                      : std::to_string(data.reaction);
                break;
              case semantic::v1::DataRole::robin_coefficient:
                record.representation = "scalar Function at boundary quadrature";
                record.provenance = data.general_scalar->provenance.robin_coefficient;
                break;
              case semantic::v1::DataRole::robin_source:
                record.representation = "scalar Function at boundary quadrature";
                record.provenance = data.general_scalar->provenance.robin_source;
                break;
              case semantic::v1::DataRole::observation_weight:
                record.representation = "scalar Function at boundary face quadrature";
                record.provenance = data.weighted_trace->provenance;
                break;
              case semantic::v1::DataRole::regularisation_weight:
                record.representation = "scalar constant";
                record.provenance = std::to_string(data.regularisation_weight);
                break;
              case semantic::v1::DataRole::lower_bound:
              case semantic::v1::DataRole::upper_bound:
                record.representation =
                  bound_binding_description(bounds, facewise_bounds);
                record.provenance = "caller-supplied compiled bound data";
                break;
              case semantic::v1::DataRole::unspecified:
                record.representation = "unspecified";
                record.provenance = "unspecified";
                break;
            }

          const auto record_scalar = [&record](const double value) {
            record.scalar_value = value;
            record.value_status = CompiledBindingStatus::checked;
            record.value_digest = hash_text(double_bits(value));
          };
          if (binding.role == semantic::v1::DataRole::diffusion &&
              !uses_coefficient_identification && !uses_general_scalar &&
              data.diffusion)
            record_scalar(*data.diffusion);
          else if (binding.role == semantic::v1::DataRole::reaction &&
                   !uses_general_scalar)
            record_scalar(data.reaction);
          else if (binding.role == semantic::v1::DataRole::regularisation_weight)
            record_scalar(data.regularisation_weight);
          else if (binding.role == semantic::v1::DataRole::lower_bound &&
                   bounds)
            std::visit(
              [&record, &record_scalar](const auto &value) {
                using Value = std::decay_t<decltype(value)>;
                if constexpr (std::is_same_v<Value, double>)
                  record_scalar(value);
                else
                  {
                    record.value_digest = vector_identity(value);
                    record.value_status = CompiledBindingStatus::checked;
                  }
              },
              bounds->lower);
          else if (binding.role == semantic::v1::DataRole::upper_bound &&
                   bounds)
            std::visit(
              [&record, &record_scalar](const auto &value) {
                using Value = std::decay_t<decltype(value)>;
                if constexpr (std::is_same_v<Value, double>)
                  record_scalar(value);
                else
                  {
                    record.value_digest = vector_identity(value);
                    record.value_status = CompiledBindingStatus::checked;
                  }
              },
              bounds->upper);
          else if (binding.role == semantic::v1::DataRole::lower_bound &&
                   facewise_bounds)
            std::visit(
              [&record, &record_scalar](const auto &value) {
                using Value = std::decay_t<decltype(value)>;
                if constexpr (std::is_same_v<Value, double>)
                  record_scalar(value);
                else
                  {
                    record.value_digest = vector_identity(value);
                    record.value_status = CompiledBindingStatus::checked;
                  }
              },
              facewise_bounds->lower);
          else if (binding.role == semantic::v1::DataRole::upper_bound &&
                   facewise_bounds)
            std::visit(
              [&record, &record_scalar](const auto &value) {
                using Value = std::decay_t<decltype(value)>;
                if constexpr (std::is_same_v<Value, double>)
                  record_scalar(value);
                else
                  {
                    record.value_digest = vector_identity(value);
                    record.value_status = CompiledBindingStatus::checked;
                  }
              },
              facewise_bounds->upper);
          records.push_back(std::move(record));
        }
      return records;
    }

    static std::string
    target_id(const CompiledTargetKind target)
    {
      return "compiled_target:" + std::to_string(static_cast<int>(target));
    }

    static std::string
    enum_id(const int value)
    {
      return "semantic_enum:" + std::to_string(value);
    }

    template <typename Backend>
    static CompiledSuppliedOTDRecord
    make_supplied_otd_record(
      const contract::SuppliedOTDSystemT<Backend> &system,
      const semantic::v1::SuppliedOTDDeclaration & declaration,
      const unsigned int                            state_degree)
    {
      CompiledSuppliedOTDRecord record;
      record.present = true;
      record.declaration = declaration;
      record.variable_layout = system.variable_layout()->label();
      record.residual_layout = system.residual_layout()->label();
      for (std::size_t block = 0;
           block < system.variable_layout()->n_blocks();
           ++block)
        {
          record.variable_space_ids.push_back(
            system.variable_layout()->space(block).value);
          record.variable_dimensions.push_back(
            system.variable_layout()->dimension(block));
          record.residual_space_ids.push_back(
            system.residual_layout()->space(block).value);
          record.residual_dimensions.push_back(
            system.residual_layout()->dimension(block));
        }
      const auto &selection = system.block_selection();
      const auto require_declared_layout =
        [&system](const std::size_t variable_block,
                  const std::size_t residual_block,
                  const semantic::v1::SuppliedOTDBlockSpec &block) {
          contract::require(
            system.variable_layout()->space(variable_block).value ==
              block.runtime_variable_space_id,
            "Supplied OTD variable layout disagrees with its declaration");
          contract::require(
            system.residual_layout()->space(residual_block).value ==
              block.runtime_residual_space_id,
            "Supplied OTD residual layout disagrees with its declaration");
          contract::require(
            system.variable_layout()->dimension(variable_block) ==
              system.residual_layout()->dimension(residual_block),
            "Supplied OTD block dimensions disagree between variables and residuals");
        };
      require_declared_layout(selection.state_variable,
                              selection.state_equation,
                              declaration.state_block);
      require_declared_layout(selection.adjoint_variable,
                              selection.adjoint_equation,
                              declaration.adjoint_block);
      require_declared_layout(selection.control_variable,
                              selection.control_stationarity,
                              declaration.control_stationarity_block);
      record.state_variable_block = selection.state_variable;
      record.adjoint_variable_block = selection.adjoint_variable;
      record.control_variable_block = selection.control_variable;
      record.state_equation_block = selection.state_equation;
      record.adjoint_equation_block = selection.adjoint_equation;
      record.control_stationarity_block = selection.control_stationarity;
      const std::string multiplier_convention =
        declaration.multiplier_convention ==
            semantic::v1::SuppliedOTDMultiplierConvention::framework_adjoint
          ? "framework adjoint"
          : "negative framework adjoint";
      const std::string multiplier_conversion =
        declaration.multiplier_conversion ==
            semantic::v1::SuppliedOTDMultiplierConversion::identity
          ? "identity"
          : "negate";
      record.sign_convention =
        "declared multiplier convention: " + multiplier_convention +
        "; framework conversion: " + multiplier_conversion;
      record.discretisation_provenance =
        declaration.state_block.discretisation_provenance + "; " +
        declaration.adjoint_block.discretisation_provenance + "; " +
        declaration.control_stationarity_block.discretisation_provenance +
        "; deal.II state degree " + std::to_string(state_degree);
      record.state_block_provenance = declaration.state_block.action_provenance;
      record.adjoint_block_provenance = declaration.adjoint_block.action_provenance;
      record.stationarity_block_provenance =
        declaration.control_stationarity_block.action_provenance;
      record.value_action_provenance = declaration.value_action_provenance;
      record.jvp_action_provenance = declaration.jvp_action_provenance;
      record.vjp_action_provenance = declaration.vjp_action_provenance;
      record.solve_provenance = declaration.solve_provenance;
      switch (declaration.comparison_status)
        {
          case semantic::v1::SuppliedOTDComparisonStatus::not_compared:
            record.comparison_status = "not compared: ";
            break;
          case semantic::v1::SuppliedOTDComparisonStatus::different:
            record.comparison_status = "different: ";
            break;
          case semantic::v1::SuppliedOTDComparisonStatus::
            equivalent_under_declared_conversion:
            record.comparison_status =
              "equivalence verified under declared conversion: ";
            break;
          case semantic::v1::SuppliedOTDComparisonStatus::unspecified:
            contract::require(false,
                              "Supplied OTD comparison status is unspecified");
        }
      record.comparison_status += declaration.comparison_evidence;
      return record;
    }

    static std::string
    handler_for(const ScalarLoweringPlan *plan, const std::string &component_id)
    {
      if (plan == nullptr)
        return "registered-target";
      for (const auto &term : plan->residual_terms)
        if (term.component_id == component_id)
          return term.handler_id;
      for (const auto &observation : plan->observations)
        if (observation.component_id == component_id)
          return observation.handler_id;
      for (const auto &loss : plan->losses)
        if (loss.component_id == component_id)
          return loss.handler_id;
      return "registered-target";
    }

    template <int dim>
    static ResolvedCompilationDecision
    make_resolved_decision(const semantic::v1::ProblemSpec &specification,
                           const DealiiDiscretisationPolicy &policy,
                           const CompiledTargetKind target,
                           const ResolvedCompilationRequest &request,
                           const ScalarLoweringPlan *scalar_plan,
                           const DealiiDataBindings<dim> &data,
                           const std::optional<CellwiseBoxDataBindings> &bounds,
                           const std::optional<FacewiseBoxDataBindings> &facewise_bounds,
                           const dealii::Triangulation<dim> &triangulation,
                           const std::string &mesh_provenance,
                           const bool owns_mesh)
    {
      ResolvedCompilationDecision decision;
      decision.semantic_problem_id = specification.id;
      decision.formulation_id = specification.formulation.id;
      decision.target_id = target_id(target);
      decision.formulation_record = {
        specification.formulation.id,
        specification.formulation.kind,
        specification.formulation.provenance,
        ExecutionRealisation::assembled,
        "tested dual coefficients with dot pairing"};
      decision.mesh_record = {
        dim,
        triangulation.n_active_cells(),
        mesh_provenance,
        owns_mesh ? MeshLifetimePolicy::owned_session
                  : MeshLifetimePolicy::borrowed_immutable,
        mesh_structural_identity(triangulation)};
      for (const auto &region : specification.regions)
        decision.regions.push_back(
          {region.id,
           region.kind,
           region.is_full_domain,
           region.boundary_ids,
           region.material_ids,
           region.point_coordinates.size()});
      decision.bindings = make_resolved_binding_records<dim>(
        specification, target, request, data, bounds, facewise_bounds);
      for (const auto &space : specification.spaces)
        decision.spaces.push_back(
          {space.id,
           space.role,
           compiled_space_runtime_role(space.role),
           space.region_id,
           compiled_space_finite_element(space, target, policy),
           0});
      for (const auto &pairing : specification.pairings)
        decision.pairings.push_back(
          {pairing.id, pairing.primal_space_id, pairing.covector_space_id});
      for (const auto &term : specification.residual_terms)
        decision.residuals.push_back(
          {term.id,
           enum_id(static_cast<int>(term.kind)),
           "resolved",
           handler_for(scalar_plan, term.id),
           term.variable_ids,
           term.data_ids,
           term.region_id,
           term.kind});
      for (const auto &observation : specification.observations)
        {
          std::vector<std::string> input_ids{observation.input_variable_id};
          input_ids.insert(input_ids.end(),
                           observation.data_ids.begin(),
                           observation.data_ids.end());
          std::string realization_id = "resolved";
          if (observation.kind ==
                semantic::v1::ObservationKind::weighted_boundary_trace &&
              request.weighted_trace_selection)
            realization_id = request.weighted_trace_selection->id;
          decision.observations.push_back(
            {observation.id,
             enum_id(static_cast<int>(observation.kind)),
             realization_id,
             handler_for(scalar_plan, observation.id),
             std::move(input_ids),
             {observation.output_space_id, observation.output_pairing_id},
             observation.region_id,
             semantic::v1::ResidualTermKind::unspecified,
             observation.kind});
        }
      for (const auto &loss : specification.losses)
        decision.losses.push_back(
          {loss.id,
           enum_id(static_cast<int>(loss.kind)),
           "resolved",
           handler_for(scalar_plan, loss.id),
           {loss.source_observation_id, loss.data_id},
           {loss.pairing_id},
           {},
           semantic::v1::ResidualTermKind::unspecified,
           semantic::v1::ObservationKind::unspecified,
           loss.kind});
      for (const auto &transformation : specification.transformations)
        decision.transformations.push_back(
          {transformation.id,
           enum_id(static_cast<int>(transformation.kind)),
           "resolved",
           scalar_plan == nullptr ? "registered-target"
                                   : scalar_plan->transformation_handler_id,
           {transformation.input_variable_id,
            transformation.fixed_data_id,
            transformation.control_variable_id},
           {transformation.output_space_id},
           {},
           semantic::v1::ResidualTermKind::unspecified,
           semantic::v1::ObservationKind::unspecified,
           semantic::v1::LossKind::unspecified,
           transformation.kind});
      if (scalar_plan != nullptr)
        decision.boundary_realisation = scalar_plan->boundary_selection;
      decision.realized_maps = make_resolved_maps(specification);
      decision.transposition_realisation = request.transposition_selection;
      decision.partial_boundary_selection = request.partial_boundary_selection;
      decision.fractional_metric_selection = request.fractional_metric_selection;
      decision.boundary_h1_metric_selection =
        request.boundary_h1_metric_selection;
      decision.h1_target_data_membership_selection =
        request.h1_target_data_membership_selection;
      for (const auto &requirement : specification.requirement_policies)
        decision.assumptions.push_back(
          {requirement.id,
           requirement.subject_id,
           requirement.kind,
           requirement.status,
           requirement.scope,
           requirement.region_id,
           requirement.status == semantic::v1::RequirementStatus::user_assumed
             ? requirement.selected_policy
             : std::string{}});
      decision.metric_record.semantic_id = specification.formulation.metric_id;
      decision.constraint_record.semantic_id =
        specification.formulation.constraint_id;

      auto &inventory = decision.compatibility;
      inventory.region_ids = identifiers(specification.regions);
      inventory.space_ids = identifiers(specification.spaces);
      inventory.pairing_ids = identifiers(specification.pairings);
      inventory.variable_ids = identifiers(specification.variables);
      inventory.data_ids = identifiers(specification.data);
      inventory.transformation_ids = identifiers(specification.transformations);
      inventory.residual_term_ids = identifiers(specification.residual_terms);
      inventory.observation_ids = identifiers(specification.observations);
      inventory.loss_ids = identifiers(specification.losses);
      inventory.metric_ids = identifiers(specification.metrics);
      inventory.constraint_ids = identifiers(specification.constraints);

      const auto by_semantic_id = [](const auto &left, const auto &right) {
        return left.semantic_id < right.semantic_id;
      };
      std::sort(decision.regions.begin(),
                decision.regions.end(),
                by_semantic_id);
      std::sort(decision.spaces.begin(), decision.spaces.end(), by_semantic_id);
      std::sort(decision.bindings.begin(),
                decision.bindings.end(),
                by_semantic_id);
      std::sort(decision.residuals.begin(),
                decision.residuals.end(),
                by_semantic_id);
      std::sort(decision.observations.begin(),
                decision.observations.end(),
                by_semantic_id);
      std::sort(decision.losses.begin(), decision.losses.end(), by_semantic_id);
      std::sort(decision.transformations.begin(),
                decision.transformations.end(),
                by_semantic_id);
      std::sort(decision.realized_maps.begin(),
                decision.realized_maps.end(),
                by_semantic_id);
      std::sort(decision.pairings.begin(),
                decision.pairings.end(),
                [](const CompiledPairingRecord &left,
                   const CompiledPairingRecord &right) {
                  return left.pairing_id < right.pairing_id;
                });
      std::sort(decision.assumptions.begin(),
                decision.assumptions.end(),
                [](const CompiledAssumptionRecord &left,
                   const CompiledAssumptionRecord &right) {
                  return left.id < right.id;
                });
      return decision;
    }

    static std::string
    observation_realization_id(const semantic::v1::ObservationKind kind,
                               const bool                         control_map)
    {
      if (control_map)
        return "coefficient_restriction";
      switch (kind)
        {
          case semantic::v1::ObservationKind::volume_restriction:
            return "fe_q_coefficient_restriction";
          case semantic::v1::ObservationKind::h1_state_restriction:
            return "fe_q_coefficient_h1_restriction";
          case semantic::v1::ObservationKind::boundary_trace:
          case semantic::v1::ObservationKind::boundary_restriction:
          case semantic::v1::ObservationKind::weighted_boundary_trace:
            return "ordered_boundary_face_quadrature_trace";
          case semantic::v1::ObservationKind::point_sensor:
            return "ordered_point_sensor_values";
          case semantic::v1::ObservationKind::normal_flux:
            return "ordered_normal_flux_face_quadrature";
          case semantic::v1::ObservationKind::unspecified:
            return "unspecified_observation_map";
        }
      return "unspecified_observation_map";
    }

    static std::string
    observation_output_layout(const semantic::v1::ObservationKind kind,
                              const bool                         control_map)
    {
      if (control_map)
        return "ordered control coefficient values";
      switch (kind)
        {
          case semantic::v1::ObservationKind::point_sensor:
            return "ordered scalar sensor values";
          case semantic::v1::ObservationKind::normal_flux:
            return "ordered scalar face-quadrature values";
          case semantic::v1::ObservationKind::boundary_trace:
          case semantic::v1::ObservationKind::boundary_restriction:
          case semantic::v1::ObservationKind::weighted_boundary_trace:
            return "ordered scalar boundary face-quadrature values";
          case semantic::v1::ObservationKind::volume_restriction:
          case semantic::v1::ObservationKind::h1_state_restriction:
            return "scalar FE_Q state coefficients";
          case semantic::v1::ObservationKind::unspecified:
            return "unspecified observation values";
        }
      return "unspecified observation values";
    }

    static std::vector<CompiledRealizedMapRecord>
    make_resolved_maps(const semantic::v1::ProblemSpec &specification)
    {
      std::vector<CompiledRealizedMapRecord> maps;
      for (const auto &observation : specification.observations)
        {
          const auto *input =
            find_variable(specification, observation.input_variable_id);
          const bool control_map =
            input != nullptr &&
            (input->role == semantic::v1::VariableRole::control ||
             input->role == semantic::v1::VariableRole::parameter);
          const std::string source_space_id =
            input == nullptr ? std::string{} : input->space_id;
          const std::string transformation_chain =
            input == nullptr || input->physical_field_transform_id.empty()
              ? "identity"
              : input->physical_field_transform_id;
          maps.push_back(
            {observation.id,
             {source_space_id},
             source_space_id,
             observation.output_space_id,
             {0},
             0,
             0,
             observation_realization_id(observation.kind, control_map),
             "independent variable-layout coefficients",
             observation_output_layout(observation.kind, control_map),
             control_map
               ? "declared control coefficient ordering"
               : observation.kind == semantic::v1::ObservationKind::point_sensor ||
                     observation.kind == semantic::v1::ObservationKind::normal_flux ||
                     observation.kind == semantic::v1::ObservationKind::boundary_trace ||
                     observation.kind == semantic::v1::ObservationKind::boundary_restriction ||
                     observation.kind == semantic::v1::ObservationKind::weighted_boundary_trace
                   ? "declared mesh traversal and quadrature order"
                   : "deal.II global state DoF ordering",
             observation.output_pairing_id.empty()
               ? "assembled coefficient pairing"
               : "declared pairing " + observation.output_pairing_id,
             transformation_chain,
             observation.kind == semantic::v1::ObservationKind::point_sensor
               ? "FE_Q shape evaluation at immutable sensor coordinates"
             : observation.kind == semantic::v1::ObservationKind::normal_flux
               ? "FEFaceValues outward normal derivative"
               : "selected observation operator value",
             "same realized observation operator applied to the tangent",
             observation.kind == semantic::v1::ObservationKind::normal_flux
               ? "weighted transpose of ordered face-quadrature flux samples"
               : "transpose of the realized observation operator",
             observation.output_pairing_id});
        }
      for (const auto &transformation : specification.transformations)
        {
          const auto *input =
            find_variable(specification, transformation.input_variable_id);
          const auto *control = transformation.control_variable_id.empty()
                                  ? nullptr
                                  : find_variable(specification,
                                                  transformation.control_variable_id);
          const std::string source_space_id =
            input == nullptr ? std::string{} : input->space_id;
          std::vector<std::string> input_space_ids;
          if (input != nullptr)
            input_space_ids.push_back(input->space_id);
          if (control != nullptr)
            input_space_ids.push_back(control->space_id);
          maps.push_back(
            {transformation.id,
             std::move(input_space_ids),
             source_space_id,
             transformation.output_space_id,
             {},
             0,
             0,
             transformation.kind ==
                 semantic::v1::TransformationKind::fixed_dirichlet_reconstruction
               ? "fixed_dirichlet_reconstruction"
               : "dirichlet_control_lifting",
             "independent variable-layout coefficients",
             "physical scalar FE_Q coefficients",
             "deal.II global state DoF ordering",
             "state-space reconstruction pairing",
             transformation.control_variable_id.empty()
               ? "one state input"
               : "state plus control inputs",
             "physical field reconstruction",
             "reconstruction applied to the tangent",
             "pullback of the physical covector through the reconstruction",
             {}});
        }
      return maps;
    }

    static const CompiledSpaceRecord *
    find_compiled_space(const ResolvedCompilationDecision &decision,
                        const std::string &                 semantic_id)
    {
      const auto space = std::find_if(
        decision.spaces.begin(),
        decision.spaces.end(),
        [&semantic_id](const CompiledSpaceRecord &candidate) {
          return candidate.semantic_id == semantic_id;
        });
      return space == decision.spaces.end() ? nullptr : &*space;
    }

    template <int dim>
    static std::vector<CompiledRealizedMapRecord>
    finalize_realized_maps(
      const ResolvedCompilationDecision &decision,
      const contract::ExecutableModelT<dealii_backend::SerialBackend> &executable)
    {
      auto maps = decision.realized_maps;
      const std::size_t physical_dimension =
        physical_state_dimension<dim>(executable);
      const std::size_t observation_dimension =
        realized_observation_dimension<dim>(executable);
      const auto runtime_dimension =
        [&decision, &executable](const std::string &space_id) {
          const auto *space = find_compiled_space(decision, space_id);
          contract::require(space != nullptr,
                            "A realized map references an unresolved space");
          switch (space->role)
            {
              case semantic::v1::SpaceRole::control:
              case semantic::v1::SpaceRole::parameter:
                return executable.variable_layout()->dimension(1);
              case semantic::v1::SpaceRole::test:
                return executable.test_layout()->dimension(0);
              case semantic::v1::SpaceRole::state:
                return executable.variable_layout()->dimension(0);
              case semantic::v1::SpaceRole::observation:
              case semantic::v1::SpaceRole::data:
              case semantic::v1::SpaceRole::auxiliary:
              case semantic::v1::SpaceRole::unspecified:
                return std::size_t{0};
            }
          return std::size_t{0};
        };
      for (auto &map : maps)
        {
          map.input_dimensions.clear();
          for (const auto &space_id : map.input_space_ids)
            map.input_dimensions.push_back(runtime_dimension(space_id));
          map.source_dimension = map.input_dimensions.empty()
                                   ? 0
                                   : map.input_dimensions.front();
          const auto observation = std::find_if(
            decision.observations.begin(),
            decision.observations.end(),
            [&map](const CompiledRealisationRecord &candidate) {
              return candidate.semantic_id == map.semantic_id;
            });
          if (observation == decision.observations.end())
            {
              map.output_dimension = physical_dimension;
              continue;
            }
          const auto *source_space =
            find_compiled_space(decision, map.source_space_id);
          const bool control_map =
            source_space != nullptr &&
            (source_space->role == semantic::v1::SpaceRole::control ||
             source_space->role == semantic::v1::SpaceRole::parameter);
          const bool finite_sample_map =
            observation->observation_kind ==
                semantic::v1::ObservationKind::point_sensor ||
            observation->observation_kind ==
                semantic::v1::ObservationKind::normal_flux ||
            observation->observation_kind ==
                semantic::v1::ObservationKind::boundary_trace ||
            observation->observation_kind ==
                semantic::v1::ObservationKind::boundary_restriction ||
            observation->observation_kind ==
              semantic::v1::ObservationKind::weighted_boundary_trace;
          map.output_dimension = control_map
                                   ? map.source_dimension
                                   : finite_sample_map ? observation_dimension
                                                       : physical_dimension;
        }
      return maps;
    }

    static std::vector<CompiledRealizedSpaceRecord>
    make_realized_spaces(const std::vector<CompiledRealizedMapRecord> &maps)
    {
      std::vector<CompiledRealizedSpaceRecord> spaces;
      for (const auto &map : maps)
        {
          for (std::size_t index = 0; index < map.input_space_ids.size(); ++index)
            spaces.push_back(
              {map.semantic_id,
               map.input_space_ids[index],
               "input:" + map.realization_id,
               map.input_dimensions.at(index),
               map.source_layout,
               map.ordering,
               map.pairing_id});
          spaces.push_back(
            {map.semantic_id,
             map.output_space_id,
             "output:" + map.realization_id,
             map.output_dimension,
             map.output_layout,
             map.ordering,
             map.pairing_id});
        }
      return spaces;
    }

    template <int dim>
    static ResolvedCompilationDecision
    finalize_resolved_decision(
      const DealiiDiscretisationPolicy &policy,
      const ConstraintRealisation       constraint_realisation,
      const ResolvedCompilationDecision &decision,
      const ResolvedCompilationRequest  &request,
      const contract::ExecutableModelT<dealii_backend::SerialBackend> &executable,
      const contract::MetricT<dealii_backend::SerialBackend> &metric,
      const ScalarLoweringPlan *          scalar_plan,
      const contract::SuppliedOTDSystemT<dealii_backend::SerialBackend> *
        supplied_otd_system,
      const semantic::v1::SuppliedOTDDeclaration * supplied_otd_declaration)
    {
      CompilationManifest manifest;
      manifest.resolved_decision = decision;
      manifest.realized_maps = finalize_realized_maps<dim>(decision, executable);
      manifest.realized_spaces = make_realized_spaces(manifest.realized_maps);
      manifest.resolved_decision.realized_maps = manifest.realized_maps;
      manifest.transposition_realisation = decision.transposition_realisation;
      manifest.partial_boundary_selection = decision.partial_boundary_selection;
      manifest.fractional_metric_selection = decision.fractional_metric_selection;
      manifest.boundary_h1_metric_selection =
        decision.boundary_h1_metric_selection;
      manifest.h1_target_data_membership_selection =
        decision.h1_target_data_membership_selection;
      manifest.resolved_decision.realized_spaces = manifest.realized_spaces;
      if (supplied_otd_system != nullptr)
        {
          contract::require(supplied_otd_declaration != nullptr,
                            "A supplied OTD system needs its formulation declaration");
          manifest.supplied_otd_record = make_supplied_otd_record(
            *supplied_otd_system,
            *supplied_otd_declaration,
            policy.state_degree);
          manifest.resolved_decision.supplied_otd_record =
            manifest.supplied_otd_record;
        }
      const CompiledTargetKind target = target_kind_from_request(request);
      const auto registration = dirichlet_registration_from_request(request);
      const bool uses_fixed_reconstruction = request.uses_fixed_reconstruction;
      const bool uses_hhalf_dirichlet_control =
        request.dirichlet_registration ==
        ResolvedDirichletRegistration::hhalf_control;
      const bool uses_h1_tracking_hhalf_dirichlet_control =
        request.uses_h1_tracking_hhalf_dirichlet_registration;
      const bool uses_h1_dirichlet_control = request.uses_h1_dirichlet_control;
      const bool uses_section_5_11_dirichlet_control =
        uses_hhalf_dirichlet_control ||
        uses_h1_tracking_hhalf_dirichlet_control ||
        uses_h1_dirichlet_control;
      const bool uses_dirichlet_control_lifting =
        request.uses_dirichlet_control;
      const bool uses_l2_dirichlet_control = request.uses_l2_dirichlet_control;
      const bool uses_normalized_dirichlet_control =
        request.uses_normalized_laplacian;
      const bool uses_partial_dirichlet_control =
        request.uses_partial_dirichlet_control;
      const bool uses_assembled_v1_target = request.uses_assembled_v1_target;
      const bool uses_general_scalar = request.uses_general_scalar;
      const bool uses_h1_state_observation = request.uses_h1_state_observation;
      const bool uses_weighted_boundary_trace =
        request.uses_weighted_boundary_trace;
      const bool uses_point_sensor = request.uses_point_sensor;
      const bool uses_normal_flux = request.uses_normal_flux;
      const bool uses_neumann_boundary_control =
        request.uses_neumann_boundary_control;
      const bool uses_neumann_convection = request.uses_neumann_convection;
      const bool uses_mean_zero_gauge = request.uses_mean_zero_gauge;
      const bool uses_h1_control_regularisation =
        request.uses_h1_control_regularisation_loss;
      const bool uses_h1_control_metric = request.uses_h1_control_metric;
      const bool uses_hhalf_control_metric = request.uses_hhalf_control_metric;
      const bool uses_hminus1_control_metric =
        request.uses_hminus1_control_metric;
      const auto *hminus1_selection =
        request.hminus1_metric_selection ? &*request.hminus1_metric_selection
                                         : nullptr;
      const bool uses_coefficient_identification =
        request.uses_coefficient_identification;
      const auto region_by_id = [&decision](const std::string &id) {
        const auto region = std::find_if(
          decision.regions.begin(),
          decision.regions.end(),
          [&id](const CompiledRegionRecord &candidate) {
            return candidate.semantic_id == id;
          });
        return region == decision.regions.end() ? nullptr : &*region;
      };
      const auto *tracking_region = region_by_id(request.tracking_region_id);
      contract::require(tracking_region != nullptr,
                        "The closed request omitted its tracking region");
      const auto *control_boundary_region = region_by_id(
        request.uses_homogeneous_dirichlet_continuous_control
          ? request.continuous_control_boundary_region_id
          : request.control_boundary_region_id);

      manifest.formulation_record = decision.formulation_record;
      manifest.mesh_record = decision.mesh_record;
      manifest.spaces = decision.spaces;
      for (auto &space : manifest.spaces)
        space.dimension = compiled_space_dimension(space,
                                                   executable,
                                                   manifest.realized_spaces);
      manifest.resolved_decision.spaces = manifest.spaces;
      manifest.bindings = decision.bindings;
      manifest.resolved_decision.bindings = manifest.bindings;

      const std::size_t state_dimension =
        executable.test_layout()->dimension(0);
      if (uses_mean_zero_gauge || uses_general_scalar || uses_neumann_convection)
        {
          const CompiledSolvePolicyRecord direct_record{
            LinearSolveAlgorithm::serial_sparse_direct_umfpack,
            "not applicable",
            1,
            0.0,
            0.0,
            uses_mean_zero_gauge ? "one mean-zero Lagrange multiplier"
                                 : "fixed Dirichlet",
            {}};
          manifest.state_solve_record = direct_record;
          manifest.adjoint_solve_record = direct_record;
          manifest.state_solve_record.operator_realisation =
            uses_mean_zero_gauge
              ? "augmented symmetric state and adjoint saddle systems"
            : uses_general_scalar
              ? "nonsymmetric state operator and its exact transpose"
            : "nonsymmetric conservative-transport state operator and its exact transpose";
          manifest.adjoint_solve_record.operator_realisation =
            manifest.state_solve_record.operator_realisation;
        }
      else
        {
          manifest.state_solve_record = spd_solve_record(
            policy.state_solve,
            state_dimension,
            uses_dirichlet_control_lifting
              ? "homogeneous conforming Galerkin state subspace"
              : "fixed Dirichlet");
          manifest.adjoint_solve_record = spd_solve_record(
            policy.adjoint_solve,
            state_dimension,
            uses_dirichlet_control_lifting
              ? "homogeneous conforming Galerkin adjoint subspace"
              : "fixed Dirichlet");
          manifest.state_solve_record.operator_realisation =
            uses_coefficient_identification
              ? "parameter-dependent SPD state matrix is reassembled for each state and adjoint solve"
            : uses_normalized_dirichlet_control
              ? (uses_l2_dirichlet_control
                   ? "conforming variational state and adjoint systems selected by transposition equivalence"
                   : "normalized-Laplacian conforming state and adjoint systems")
              : "symmetric positive-definite operator";
          manifest.adjoint_solve_record.operator_realisation =
            manifest.state_solve_record.operator_realisation;
        }
      manifest.metric_record = {
        decision.metric_record.semantic_id,
        metric.id(),
        uses_hminus1_control_metric
          ? "M_h K_h^{-1} M_h negative-norm Riesz map"
        : uses_hhalf_control_metric
          ? "minimum-volume-H1-extension Schur-complement Riesz map"
        : uses_h1_dirichlet_control
          ? "boundary mass plus tangential stiffness Riesz map"
        : uses_h1_control_metric ? "mass plus stiffness Riesz map"
                                 : "mass Riesz map",
        hminus1_selection == nullptr ? "" : hminus1_selection->primal_space_id,
        hminus1_selection == nullptr ? "" : hminus1_selection->dual_space_id,
        hminus1_selection == nullptr ? "" : "mass_laplacian_inverse_mass",
        hminus1_selection == nullptr
          ? ""
          : "mass_inverse_laplacian_mass_inverse",
        hminus1_selection == nullptr ? "" : hminus1_selection->mass_pairing_id,
        hminus1_selection == nullptr
          ? ""
          : hminus1_selection->laplacian_pairing_id,
        hminus1_selection == nullptr
          ? ""
          : hminus1_selection->fixed_boundary_region_id,
        hminus1_selection == nullptr
          ? ""
          : hminus1_selection->laplacian_solve_policy_id,
        hminus1_selection == nullptr
          ? ""
          : hminus1_selection->mass_solve_policy_id,
        hminus1_selection == nullptr
          ? ""
          : "fixed_dirichlet_no_nullspace",
        {LinearSolveAlgorithm::serial_cg,
         "identity",
         policy.control_metric_solve.maximum_iterations,
         policy.control_metric_solve.relative_tolerance,
         policy.control_metric_solve.absolute_tolerance,
         "not applicable",
         {}}};
      manifest.constraint_record = {
        constraint_realisation != ConstraintRealisation::none,
        decision.constraint_record.semantic_id,
        constraint_realisation_id(constraint_realisation),
        constraint_realisation == ConstraintRealisation::none ? "none"
                                                               : metric.id(),
        {},
        {},
        {},
        {},
        {}};
      manifest.resolved_decision.state_solve_record =
        manifest.state_solve_record;
      manifest.resolved_decision.adjoint_solve_record =
        manifest.adjoint_solve_record;
      manifest.resolved_decision.metric_record = manifest.metric_record;
      manifest.resolved_decision.constraint_record = manifest.constraint_record;
      manifest.execution = manifest.resolved_decision.execution_id;
      manifest.dual_representation =
        manifest.resolved_decision.formulation_record.dual_representation;
      manifest.metric_solve_policy =
        metric_solve_policy_description(manifest.metric_record);
      manifest.constraint_realisation = describe(manifest.constraint_record);
      manifest.nullspace_policy =
        manifest.resolved_decision.state_solve_record.nullspace_policy;
      if (scalar_plan != nullptr)
        {
          manifest.lowering_handler_records = scalar_plan->provenance;
          manifest.boundary_realisation = scalar_plan->boundary_selection;
        }
      if (uses_weighted_boundary_trace)
        manifest.lowering_handler_records.push_back(
          "weighted_state_boundary_trace <- "
          "dealii.neumann.observation.weighted_boundary_trace");
      if (uses_l2_dirichlet_control)
        manifest.lowering_handler_records.push_back(
          "l2_dirichlet_transposition <- "
          "dealii.dirichlet_control.conforming_trace_equivalence");
      if (uses_section_5_11_dirichlet_control)
        manifest.lowering_handler_records.push_back(
          "normalized_dirichlet_laplace <- "
          "dealii.dirichlet_control.explicit_nodal_lifting");
      if (uses_hhalf_control_metric)
        manifest.lowering_handler_records.push_back(
          "hhalf_trace_metric <- "
          "dealii.metric.minimum_volume_h1_extension");
      if (uses_h1_dirichlet_control)
        manifest.lowering_handler_records.push_back(
          "h1_trace_metric <- "
          "dealii.metric.boundary_mass_tangential_stiffness");
      manifest.semantic_problem_id = decision.semantic_problem_id;
      manifest.compiler_id =
        uses_normal_flux
          ? "nmopt.compiler.v1.dealii.normal_flux"
        : uses_point_sensor
          ? "nmopt.compiler.v1.dealii.point_sensor"
        : uses_general_scalar
          ? "nmopt.compiler.v1.dealii.general_scalar_elliptic_robin"
        : uses_neumann_convection
          ? "nmopt.compiler.v1.dealii.neumann_convection_subdomain"
        : uses_weighted_boundary_trace
          ? "nmopt.compiler.v1.dealii.weighted_boundary_trace"
        : uses_hminus1_control_metric
          ? "nmopt.compiler.v1.dealii.hminus1_control_metric"
        : uses_h1_tracking_hhalf_dirichlet_control
          ? "nmopt.compiler.v1.dealii.h1_tracking_hhalf_dirichlet_control"
        : uses_hhalf_dirichlet_control
          ? "nmopt.compiler.v1.dealii.hhalf_dirichlet_control"
        : uses_h1_dirichlet_control
          ? "nmopt.compiler.v1.dealii.h1_dirichlet_control"
        : uses_h1_state_observation
          ? "nmopt.compiler.v1.dealii.h1_state_tracking"
        : uses_l2_dirichlet_control
          ? "nmopt.compiler.v1.dealii.l2_dirichlet_transposition"
          : "nmopt.compiler.v1.dealii.scalar_diffusion_reaction";
      manifest.backend = "deal.II serial Vector<double>";
      manifest.execution = "assembled";
      manifest.state_space = uses_l2_dirichlet_control
                               ? "continuous L2(Omega) parent lowered to conforming scalar FE_Q(" +
                                   std::to_string(policy.state_degree) +
                                   ") variational Galerkin coordinates"
                               : "scalar FE_Q(" +
                                   std::to_string(policy.state_degree) + ")";
      manifest.control_space =
        control_space_description(target, policy, registration);
      manifest.quadrature = "QGauss(" +
                            std::to_string(policy.state_degree + 2) + ")";
      manifest.dual_representation = "tested dual coefficients with dot pairing";
      manifest.data_rule = uses_section_5_11_dirichlet_control
        ? "analytic forcing and desired-state Function " +
            std::string(uses_h1_state_observation ? "value and gradient"
                                                  : "value") +
            " at selected QGauss(" +
            std::to_string(policy.state_degree + 2) +
            ") volume quadrature; normalized unit-diffusion zero-reaction Laplacian; Dirichlet trace is the decision block"
        : uses_l2_dirichlet_control
        ? "analytic forcing and desired-state Functions at selected QGauss(" +
            std::to_string(policy.state_degree + 2) +
            ") volume quadrature; normalized unit-diffusion zero-reaction Laplacian; Dirichlet trace is the decision block"
        : uses_point_sensor
        ? "analytic forcing Function at selected volume quadrature; scalar operator T=-kappa Delta+rI with kappa <- " +
            request.transposition_selection->diffusion_data_id +
            " and r <- " + request.transposition_selection->reaction_data_id +
            "; desired-state Function evaluated at immutable physical sensor coordinates; FE_Q shape evaluation and assembled C_h^T point-load transpose"
        : uses_normal_flux
        ? "analytic forcing Function at selected volume quadrature; scalar operator T=-kappa Delta+rI with kappa <- " +
            request.transposition_selection->diffusion_data_id +
            " and r <- " + request.transposition_selection->reaction_data_id +
            "; desired-state Function evaluated at selected boundary face quadrature; FE_Q outward normal derivative and assembled face-map transpose"
        : uses_h1_state_observation
        ? "analytic desired-state Function value and gradient at selected QGauss(" +
            std::to_string(policy.state_degree + 2) +
            ") volume quadrature; scalar coefficients and forcing Function at volume quadrature"
        : uses_weighted_boundary_trace
        ? "analytic desired-state and fixed boundary-weight Functions at selected QGauss(" +
            std::to_string(policy.state_degree + 2) +
            ") boundary face quadrature; scalar coefficients and forcing Function at volume quadrature"
        : uses_neumann_boundary_control
        ? (uses_neumann_convection
             ? "analytic desired-state, conservative transport, and forcing Functions at selected QGauss(" +
                 std::to_string(policy.state_degree + 2) +
                 ") volume quadrature; scalar diffusion, reaction, and regularisation constants"
             : "analytic desired-state Function at selected QGauss(" +
            std::to_string(policy.state_degree + 2) +
            ") boundary face quadrature; scalar coefficients and forcing Function at volume quadrature")
        : uses_general_scalar
        ? general_scalar_data_rule(manifest.bindings, policy.state_degree + 2)
        : "analytic desired-state Function at selected QGauss(" +
            std::to_string(policy.state_degree + 2) + ") volume quadrature" +
            (uses_coefficient_identification
              ? "; forcing Function, reaction and regularisation scalars; diffusion is the parameter decision block"
              : uses_fixed_reconstruction
              ? "; fixed Dirichlet Function interpolated at boundary DoFs"
              : uses_dirichlet_control_lifting
              ? (uses_partial_dirichlet_control
                   ? "; scalar coefficients and forcing Function at volume quadrature; fixed Dirichlet Function is interpolated at fixed boundary DoFs and the relative-interior trace is the decision block"
                   : "; scalar coefficients and forcing Function at volume quadrature; Dirichlet trace is the decision block")
              : "; scalar coefficients and forcing Function at volume quadrature");
      manifest.observation_realisation = uses_h1_state_observation
        ? (uses_dirichlet_control_lifting
             ? "full-domain physical H1 state restriction with mass-plus-stiffness pairing"
             : "full-domain H1_0 state restriction with mass-plus-stiffness pairing")
        : uses_point_sensor
        ? point_sensor_observation_realisation(*tracking_region)
        : uses_normal_flux
        ? normal_flux_observation_realisation(*tracking_region)
        : uses_weighted_boundary_trace
        ? weighted_boundary_observation_realisation(*tracking_region)
        : uses_neumann_boundary_control
        ? (uses_neumann_convection
             ? observation_realisation(*tracking_region)
             : boundary_observation_realisation(*tracking_region))
        : observation_realisation(*tracking_region);
      manifest.metric_solve_policy =
        metric_solve_policy_description(manifest.metric_record);
      manifest.constraint_realisation =
        describe(manifest.constraint_record);
      manifest.lifting_realisation = uses_mean_zero_gauge
                                       ? "none; pure-Neumann state uses an explicit mean-zero gauge"
                                       : uses_fixed_reconstruction
                                       ? "y_phys = P_h y_hat + ell_0,h; independent FE_Q coordinates, AffineConstraints reconstruction, and P_h^* pullbacks"
                                       : uses_dirichlet_control_lifting
                                           ? (uses_l2_dirichlet_control
                                                ? "continuous E_tr(y,u;f) in (H2 cap H1_0)^*; U_h=trace(V_h) subset H1/2(Gamma) uses the equivalent y_phys = P_h y_hat + L_D,h u_h conforming Galerkin lifting"
                                              : uses_section_5_11_dirichlet_control
                                                ? "y_phys = P_h y_hat + L_D,h u_h; complete-boundary conforming nodal trace lifting with P_h^*/L_D,h^* pullbacks for the normalized Laplacian"
                                              : uses_partial_dirichlet_control
                                                ? "y_phys = P_h y_hat + ell_0,h + L_D,h u_h; partial nodal trace lifting with fixed-data interface precedence, independent FE_Q coordinates, and P_h^*/L_D,h^* pullbacks"
                                                : "y_phys = P_h y_hat + L_D,h u_h; complete-boundary shared nodal trace lifting, independent FE_Q coordinates, and P_h^*/L_D,h^* pullbacks")
                                       : uses_assembled_v1_target
                                           ? "y_phys = P_h y_hat; independent FE_Q coordinates and AffineConstraints reconstruction"
                                           : "homogeneous full-vector Dirichlet rows; no inhomogeneous lifting";
      manifest.nullspace_policy =
        manifest.state_solve_record.nullspace_policy;
      manifest.state_adjoint_solve_policy =
        state_adjoint_solve_policy_description(manifest.state_solve_record,
                                               manifest.adjoint_solve_record);
      manifest.provenance =
        manifest.formulation_record.provenance ==
            semantic::v1::FormulationProvenance::supplied_otd
          ? "supplied OTD"
          : "DTO";
      if (uses_neumann_boundary_control && control_boundary_region != nullptr)
        manifest.declared_assumptions.push_back(
          "neumann_control_realisation: facewise-constant FEFaceValues pairing on boundary ids " +
          boundary_id_list(*control_boundary_region));
      if (uses_neumann_convection)
        manifest.declared_assumptions.push_back(
          "neumann_convection_subdomain: conservative transport is assembled in the scalar residual; the Neumann datum is its declared conormal boundary functional and state tracking is restricted to declared material ids");
      if (uses_dirichlet_control_lifting && control_boundary_region != nullptr)
        manifest.declared_assumptions.push_back(
          std::string(uses_l2_dirichlet_control
                        ? "l2_dirichlet_transposition: continuous L2 boundary control with H2 cap H1_0 tests; complete-boundary conforming trace subspace on boundary ids "
                      : uses_section_5_11_dirichlet_control
                        ? "section_5_11_dirichlet_control: complete-boundary conforming nodal trace on boundary ids "
                      : uses_partial_dirichlet_control
                        ? "dirichlet_control_lifting: partial controlled nodal trace map on boundary ids "
                        : "dirichlet_control_lifting: complete-exterior-boundary shared nodal trace map on boundary ids ") +
          boundary_id_list(*control_boundary_region) +
          (uses_l2_dirichlet_control
             ? "; the variational lifting lowerer is valid only on U_h=trace(V_h) subset H1/2(Gamma); discontinuous or facewise controls are rejected"
           : uses_section_5_11_dirichlet_control
             ? "; normalized Laplacian, explicit lifting pullbacks, and no trace box constraint"
           : uses_partial_dirichlet_control
             ? "; fixed-data precedence owns every fixed/controlled corner or interface DoF; no hanging-node relation or box policy is registered"
             : "; no corner/interface averaging, hanging-node relation, or box policy is registered"));
      if (uses_h1_control_regularisation)
        manifest.declared_assumptions.push_back(uses_h1_dirichlet_control
          ? "h1_trace_control_regularisation: alpha/2 u^T (M_Gamma + K_Gamma^tau) u; search metric=h1_dirichlet_trace"
          : "h1_control_regularisation: alpha/2 u^T (M_u + K_u) u; search metric=" +
              std::string(uses_h1_control_metric ? "h1_continuous"
                                                 : "l2_continuous"));
      if (uses_hhalf_dirichlet_control)
        manifest.declared_assumptions.push_back(
          "hhalf_trace_control_regularisation: alpha/2 u^T G_1/2,h u with G_1/2,h the minimum-volume-H1-extension Schur complement; search metric=hhalf_dirichlet_trace");
      if (uses_h1_tracking_hhalf_dirichlet_control)
        manifest.declared_assumptions.push_back(
          "l2_trace_control_regularisation: alpha/2 u^T M_Gamma u; search metric remains the independent hhalf_dirichlet_trace Riesz map");
      if (uses_hminus1_control_metric)
        {
          contract::require(hminus1_selection != nullptr,
                            "H-1 manifest needs its typed metric realization");
          manifest.declared_assumptions.push_back(
            "hminus1_control_metric: realization=" + hminus1_selection->id +
            "; boundary=" + hminus1_selection->fixed_boundary_region_id +
            "; operator=mass_laplacian_inverse_mass; inverse=" +
            "mass_inverse_laplacian_mass_inverse; laplacian_solve=" +
            hminus1_selection->laplacian_solve_policy_id +
            "; mass_solve=" + hminus1_selection->mass_solve_policy_id +
            "; nullspace=fixed_dirichlet_no_nullspace");
        }
      if (uses_coefficient_identification)
        manifest.declared_assumptions.push_back(
          "coefficient_identification: positive cellwise physical diffusion parameter; A(m) is reassembled for every state and adjoint solve");
      if (uses_general_scalar)
        {
          contract::require(manifest.boundary_realisation.has_value(),
                            "General scalar manifest needs its typed boundary realization");
          manifest.declared_assumptions.push_back(
            "general_scalar_robin: " +
            boundary_realisation_description(*manifest.boundary_realisation));
        }
      if (uses_h1_state_observation)
        {
          contract::require(
            manifest.h1_target_data_membership_selection.has_value(),
            "H1-state manifest needs its typed target-data membership assumption");
          const auto &selection =
            *manifest.h1_target_data_membership_selection;
          manifest.declared_assumptions.push_back(
            "h1_target_data_membership: status=user_assumed; data=" +
            selection.data_id + "; observation_space=" +
            selection.observation_space_id + "; fixed_boundary=" +
            selection.fixed_boundary_region_id +
            "; regularity=h1_value_and_weak_gradient; "
            "trace=zero_trace_on_fixed_boundary");
          manifest.declared_assumptions.push_back(
            std::string(uses_dirichlet_control_lifting
                          ? "h1_state_observation: full-domain physical H1 value/gradient quadrature with mass-plus-stiffness loss pullback; control metric="
                          : "h1_state_observation: full-domain H1_0 value/gradient quadrature with mass-plus-stiffness loss pullback; control metric=") +
            std::string(uses_hminus1_control_metric ? "hminus1_continuous"
                        : uses_hhalf_control_metric ? "hhalf_dirichlet_trace"
                                                   : "L2"));
        }
      if (uses_weighted_boundary_trace)
        {
          contract::require(request.weighted_trace_selection.has_value(),
                            "Weighted trace manifest needs its typed trace realization");
          const auto &selection = *request.weighted_trace_selection;
          manifest.declared_assumptions.push_back(
            "weighted_boundary_trace: realization=" + selection.id +
            "; source=" + selection.source_space_id +
            "; output=" + selection.output_space_id +
            "; boundary=" + selection.region_id +
            "; weight=" + selection.weight_data_id +
            "; pairing=" + selection.pairing_id +
            "; Neumann residual and facewise L2 control metric are unchanged");
        }
      if (request.uses_homogeneous_dirichlet_continuous_control &&
          control_boundary_region != nullptr)
        manifest.declared_assumptions.push_back(
          "continuous_control_boundary: homogeneous FE_Q control boundary " +
          boundary_id_list(*control_boundary_region) +
          " is resolved independently and matches the state fixed boundary");
      if (uses_point_sensor)
        manifest.declared_assumptions.push_back(
          "point_sensor: immutable physical coordinates; FE_Q shape evaluation defines C_h, and the very-weak adjoint source is the assembled finite-dimensional transpose C_h^T(C_h y-z); nearest-node and quadrature-point coincidence policies are rejected");
      if (uses_normal_flux)
        manifest.declared_assumptions.push_back(
          "normal_flux: strong H2 cap H1_0 state with outward unit normal; FE_Q normal derivatives are evaluated at selected boundary face quadrature and their transpose is assembled as the very-weak adjoint boundary source");
      for (const auto &assumption : decision.assumptions)
        manifest.declared_assumptions.push_back(
          assumption.id + ": subject=" + assumption.subject_id +
          "; kind=" + enum_id(static_cast<int>(assumption.kind)) +
          "; status=" + enum_id(static_cast<int>(assumption.status)) +
          "; scope=" + enum_id(static_cast<int>(assumption.scope)) +
          "; region=" + assumption.region_id);
      manifest.region_ids = decision.compatibility.region_ids;
      manifest.space_ids = decision.compatibility.space_ids;
      manifest.pairing_ids = decision.compatibility.pairing_ids;
      manifest.variable_ids = decision.compatibility.variable_ids;
      manifest.data_ids = decision.compatibility.data_ids;
      manifest.transformation_ids = decision.compatibility.transformation_ids;
      manifest.residual_term_ids = decision.compatibility.residual_term_ids;
      manifest.observation_ids = decision.compatibility.observation_ids;
      manifest.loss_ids = decision.compatibility.loss_ids;
      manifest.metric_ids = decision.compatibility.metric_ids;
      manifest.constraint_ids = decision.compatibility.constraint_ids;
      std::sort(manifest.lowering_handler_records.begin(),
                manifest.lowering_handler_records.end());
      auto &compatibility = manifest.resolved_decision.compatibility;
      compatibility.compiler_id = manifest.compiler_id;
      compatibility.backend = manifest.backend;
      compatibility.execution = manifest.execution;
      compatibility.state_space = manifest.state_space;
      compatibility.control_space = manifest.control_space;
      compatibility.quadrature = manifest.quadrature;
      compatibility.dual_representation = manifest.dual_representation;
      compatibility.data_rule = manifest.data_rule;
      compatibility.observation_realisation = manifest.observation_realisation;
      compatibility.metric_solve_policy = manifest.metric_solve_policy;
      compatibility.constraint_realisation = manifest.constraint_realisation;
      compatibility.lifting_realisation = manifest.lifting_realisation;
      compatibility.nullspace_policy = manifest.nullspace_policy;
      compatibility.state_adjoint_solve_policy =
        manifest.state_adjoint_solve_policy;
      compatibility.provenance = manifest.provenance;
      compatibility.lowering_handler_records =
        manifest.lowering_handler_records;
      compatibility.region_ids = manifest.region_ids;
      compatibility.space_ids = manifest.space_ids;
      compatibility.pairing_ids = manifest.pairing_ids;
      compatibility.variable_ids = manifest.variable_ids;
      compatibility.data_ids = manifest.data_ids;
      compatibility.transformation_ids = manifest.transformation_ids;
      compatibility.residual_term_ids = manifest.residual_term_ids;
      compatibility.observation_ids = manifest.observation_ids;
      compatibility.loss_ids = manifest.loss_ids;
      compatibility.metric_ids = manifest.metric_ids;
      compatibility.constraint_ids = manifest.constraint_ids;
      compatibility.declared_assumptions = manifest.declared_assumptions;
      return manifest.resolved_decision;
    }

    static CompiledKKTRecord
    make_kkt_record(
      const contract::EqualityConstrainedQuadraticKKTProductT<
        dealii_backend::SerialBackend> &product,
      const bool supplied_otd = false)
    {
      const auto &layout = product.layout();
      CompiledKKTRecord record;
      record.present = true;
      record.product_id = supplied_otd ? "compiled.scalar.supplied_otd.kkt"
                                       : "compiled.scalar.dto.kkt";
      record.construction_realisation = supplied_otd
        ? "compiler-owned adapter from the canonical supplied-OTD KKT bridge"
        : "compiler-owned adapter from canonical DTO objective/residual action ports";
      record.primal_layout = layout.primal->label();
      record.multiplier_layout = layout.multiplier->label();
      record.adjoint_layout = layout.adjoint->label();
      record.stationarity_layout = layout.stationarity->label();
      record.equality_layout = layout.equality->label();
      record.primal_stationarity_pairing =
        "primal '" + record.primal_layout + "' <-> stationarity '" +
        record.stationarity_layout + "'";
      record.multiplier_equality_pairing =
        "multiplier '" + record.multiplier_layout + "' <-> equality '" +
        record.equality_layout + "'";
      record.primal_stationarity_pairing_ids =
        layout.primal_stationarity_pairing.pairing_ids;
      record.multiplier_equality_pairing_ids =
        layout.multiplier_equality_pairing.pairing_ids;
      record.multiplier_conversion =
        "KKT multiplier equals negative framework adjoint";
      record.rank_condition_declared = product.assumptions().rank_condition_declared;
      record.rank_policy = product.assumptions().rank_policy;
      record.kernel_positivity_declared =
        product.assumptions().kernel_positivity_declared;
      record.kernel_policy = product.assumptions().kernel_policy;
      record.symmetry = product.supports_minres() ? "symmetric_indefinite"
                                                  : "nonsymmetric";
      record.solver_policy =
        "MINRES for the declared symmetric-indefinite product; GMRES requires a later explicit nonsymmetric target";
      record.preconditioner = "identity baseline";
      record.d_transpose_consistency_declared =
        product.assumptions().d_transpose_consistency_declared;
      record.kkt_transpose_consistency_declared =
        product.assumptions().kkt_transpose_consistency_declared;
      record.transpose_consistency_policy =
        product.assumptions().transpose_consistency_policy;
      record.action_provenance = {
        "Q <- objective_derivative(primal) - objective_derivative(zero)",
        "D <- residual_jvp(zero, primal)",
        "D^T <- residual_vjp(zero, multiplier)",
        "KKT transpose <- Q + residual_vjp applied to the equality seed",
        "RHS <- negative objective_derivative(zero) and residual(zero)"};
      record.assembled_block_provenance = {
        "Q[0,0] <- assembled state-tracking Hessian",
        "Q[1,1] <- assembled control-regularisation Hessian",
        "Q[0,1] and Q[1,0] <- zero canonical cross blocks",
        "D[0,0] <- assembled state residual Jacobian",
        "D[0,1] <- assembled volume-control residual coupling",
        "D^T <- exact coefficient transpose of the assembled D blocks"};
      return record;
    }

    static CompiledPDASRecord
    make_pdas_record(
      const contract::EqualityConstrainedQuadraticKKTProductT<
        dealii_backend::SerialBackend> &product,
      const contract::BoxComplementarityT<dealii_backend::SerialBackend>
        &complementarity,
      const CompiledCellwiseBoxDataT<dealii_backend::SerialBackend> &box_data,
      const contract::MetricT<dealii_backend::SerialBackend> &metric,
      const DealiiDiscretisationPolicy &policy,
      const bool supplied_otd)
    {
      const auto &active_set = policy.pdas.active_set_assumptions;
      CompiledPDASRecord record;
      record.present = true;
      record.product_id = supplied_otd
                            ? "compiled.scalar.supplied_otd.pdas"
                            : "compiled.scalar.dto.pdas";
      record.construction_realisation = supplied_otd
        ? "compiler-owned cellwise box product over the canonical supplied-OTD KKT bridge"
        : "compiler-owned cellwise box product over the canonical DTO KKT product";
      record.bound_source = "semantic cellwise_box constraint and caller-supplied lower/upper data";
      record.bound_realisation = "FE_DGQ(0) coefficientwise two-sided bounds";
      record.control_block = 1;
      record.control_layout = product.layout().primal->space(1).value +
                              " (" +
                              std::to_string(product.layout().primal->dimension(1)) +
                              " coefficients)";
      record.box_data_token = box_data.token_id();
      record.bounds_digest = box_data.bounds_digest();
      record.control_ordering = box_data.layout_signature();
      record.data_provenance = box_data.data_provenance();
      record.multiplier_representation =
        complementarity.multiplier_representation().description;
      record.metric_realisation = metric.id();
      const auto *mass_metric =
        dynamic_cast<const dealii_backend::MassMetric *>(&metric);
      record.positive_diagonal_metric_declared =
        mass_metric != nullptr &&
        mass_metric->supports_coefficientwise_box_projection();
      record.classification_parameter = policy.pdas.classification_parameter;
      record.maximum_iterations = policy.pdas.maximum_iterations;
      record.primal_feasibility_tolerance =
        policy.pdas.primal_feasibility_tolerance;
      record.dual_feasibility_tolerance =
        policy.pdas.dual_feasibility_tolerance;
      record.complementarity_tolerance =
        policy.pdas.complementarity_tolerance;
      record.stationarity_tolerance = policy.pdas.stationarity_tolerance;
      record.equality_tolerance = policy.pdas.equality_tolerance;
      record.inner_kkt_solver =
        policy.pdas_kkt_solver.method ==
            contract::QuadraticKKTSolverMethod::minres
          ? "serial MINRES"
          : "serial GMRES";
      record.inner_kkt_maximum_iterations =
        policy.pdas_kkt_solver.maximum_iterations;
      record.inner_kkt_relative_tolerance =
        policy.pdas_kkt_solver.relative_tolerance;
      record.inner_kkt_absolute_tolerance =
        policy.pdas_kkt_solver.absolute_tolerance;
      record.active_set_rank_condition_declared =
        active_set.rank_condition_declared;
      record.active_set_rank_policy = active_set.rank_policy;
      record.active_set_kernel_positivity_declared =
        active_set.kernel_positivity_declared;
      record.active_set_kernel_policy = active_set.kernel_policy;
      record.exclusions = {
        "continuous control layouts",
        "facewise or quadrature-point bounds",
        "mixed or state constraints",
        "non-positive-diagonal or non-L2 multiplier metrics"};
      return record;
    }

    static CompilationManifest
    make_manifest(const ResolvedCompilationDecision &decision)
    {
      CompilationManifest manifest;
      manifest.resolved_decision = decision;
      manifest.formulation_record = decision.formulation_record;
      manifest.supplied_otd_record = decision.supplied_otd_record;
      manifest.kkt_record = decision.kkt_record;
      manifest.pdas_record = decision.pdas_record;
      manifest.mesh_record = decision.mesh_record;
      manifest.spaces = decision.spaces;
      manifest.bindings = decision.bindings;
      manifest.state_solve_record = decision.state_solve_record;
      manifest.adjoint_solve_record = decision.adjoint_solve_record;
      manifest.metric_record = decision.metric_record;
      manifest.constraint_record = decision.constraint_record;
      manifest.realized_spaces = decision.realized_spaces;
      manifest.realized_maps = decision.realized_maps;
      manifest.boundary_realisation = decision.boundary_realisation;
      manifest.transposition_realisation = decision.transposition_realisation;
      manifest.partial_boundary_selection = decision.partial_boundary_selection;
      manifest.fractional_metric_selection = decision.fractional_metric_selection;
      manifest.boundary_h1_metric_selection =
        decision.boundary_h1_metric_selection;
      manifest.h1_target_data_membership_selection =
        decision.h1_target_data_membership_selection;
      const auto &compatibility = decision.compatibility;
      manifest.lowering_handler_records =
        compatibility.lowering_handler_records;
      manifest.semantic_problem_id = decision.semantic_problem_id;
      manifest.compiler_id = compatibility.compiler_id;
      manifest.backend = compatibility.backend;
      manifest.execution = compatibility.execution;
      manifest.state_space = compatibility.state_space;
      manifest.control_space = compatibility.control_space;
      manifest.quadrature = compatibility.quadrature;
      manifest.dual_representation = compatibility.dual_representation;
      manifest.data_rule = compatibility.data_rule;
      manifest.observation_realisation = compatibility.observation_realisation;
      manifest.metric_solve_policy = compatibility.metric_solve_policy;
      manifest.constraint_realisation = compatibility.constraint_realisation;
      manifest.lifting_realisation = compatibility.lifting_realisation;
      manifest.nullspace_policy = compatibility.nullspace_policy;
      manifest.state_adjoint_solve_policy =
        compatibility.state_adjoint_solve_policy;
      manifest.provenance = compatibility.provenance;
      manifest.region_ids = compatibility.region_ids;
      manifest.space_ids = compatibility.space_ids;
      manifest.pairing_ids = compatibility.pairing_ids;
      manifest.variable_ids = compatibility.variable_ids;
      manifest.data_ids = compatibility.data_ids;
      manifest.transformation_ids = compatibility.transformation_ids;
      manifest.residual_term_ids = compatibility.residual_term_ids;
      manifest.observation_ids = compatibility.observation_ids;
      manifest.loss_ids = compatibility.loss_ids;
      manifest.metric_ids = compatibility.metric_ids;
      manifest.constraint_ids = compatibility.constraint_ids;
      manifest.declared_assumptions = compatibility.declared_assumptions;
      return manifest;
    }

    static std::string
    observation_realisation(const CompiledRegionRecord &region)
    {
      if (region.is_full_domain)
        return "full-domain volume restriction";
      std::string result = "material-id volume restriction: ";
      for (std::size_t index = 0; index < region.material_ids.size(); ++index)
        {
          if (index != 0)
            result += ",";
          result += std::to_string(region.material_ids[index]);
        }
      return result;
    }

    static std::string
    boundary_id_list(const CompiledRegionRecord &region)
    {
      std::string result;
      for (std::size_t index = 0; index < region.boundary_ids.size(); ++index)
        {
          if (index != 0)
            result += ",";
          result += std::to_string(region.boundary_ids[index]);
        }
      return result;
    }

    static std::string
    boundary_observation_realisation(const CompiledRegionRecord &region)
    {
      return "boundary trace restriction on boundary ids " +
             boundary_id_list(region);
    }

    static std::string
    weighted_boundary_observation_realisation(
      const CompiledRegionRecord &region)
    {
      return "fixed-data weighted boundary trace on boundary ids " +
             boundary_id_list(region);
    }

    static std::string
    point_sensor_observation_realisation(
      const CompiledRegionRecord &region)
    {
      return "finite point-sensor evaluation at " +
             std::to_string(region.point_count) +
             " immutable physical coordinates with assembled FE_Q transpose";
    }

    static std::string
    normal_flux_observation_realisation(
      const CompiledRegionRecord &region)
    {
      return "outward normal-flux evaluation on boundary ids " +
             boundary_id_list(region) +
             " with assembled FE_Q face-quadrature transpose";
    }

    DealiiCapabilityRegistryV1 capabilities_;
    DealiiScalarLoweringPlanner scalar_planner_;
  };
} // namespace nmopt::compiler::v1
