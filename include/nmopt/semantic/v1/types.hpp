#pragma once

#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <vector>

// The v1 semantic graph deliberately contains no discretisation-backend
// objects. Values such as deal.II Functions and FE choices are bound by a
// compiler after this graph has been validated.
namespace nmopt::semantic::v1
{
  enum class RegionKind
  {
    unspecified = -1,
    volume,
    boundary,
    point_set
  };

  enum class SpaceTopology
  {
    unspecified = -1,
    h1,
    h2,
    hhalf,
    l2,
    // Bounded coefficient/data fields are distinct from L2 data: the
    // selected weak form records the regularity assumption explicitly.
    bounded_function
  };

  enum class SpaceRole
  {
    unspecified = -1,
    state,
    test,
    control,
    parameter,
    observation,
    data,
    auxiliary
  };

  enum class VariableRole
  {
    unspecified = -1,
    state,
    control,
    parameter
  };

  enum class DataKind
  {
    unspecified = -1,
    function,
    vector_function,
    tensor_function,
    scalar_constant,
    cellwise_bound,
    facewise_bound
  };

  enum class DataRole
  {
    unspecified = -1,
    forcing,
    desired_state,
    fixed_dirichlet_lifting,
    diffusion,
    conservative_transport,
    advective_transport,
    reaction,
    robin_coefficient,
    robin_source,
    regularisation_weight,
    lower_bound,
    upper_bound,
    observation_weight
  };

  enum class ResidualTermKind
  {
    unspecified = -1,
    diffusion_reaction,
    tensor_diffusion,
    conservative_transport,
    advective_transport,
    reaction,
    parameter_diffusion_reaction,
    laplacian,
    transposition_laplacian,
    dirichlet_transposition_control,
    volume_source,
    volume_control,
    neumann_control,
    robin_bilinear,
    robin_source
  };

  enum class ObservationKind
  {
    unspecified = -1,
    volume_restriction,
    h1_state_restriction,
    boundary_trace,
    boundary_restriction,
    weighted_boundary_trace,
    point_sensor,
    normal_flux
  };

  enum class LossKind
  {
    unspecified = -1,
    quadratic_tracking,
    quadratic_control_regularisation,
    quadratic_hhalf_control_regularisation,
    quadratic_h1_control_regularisation,
    quadratic_parameter_regularisation
  };

  enum class MetricKind
  {
    unspecified = -1,
    l2,
    hhalf,
    h1,
    hminus1
  };

  enum class ConstraintKind
  {
    unspecified = -1,
    cellwise_box,
    facewise_box
  };

  enum class TransformationKind
  {
    unspecified = -1,
    fixed_dirichlet_reconstruction,
    dirichlet_control_lifting
  };

  enum class RequirementKind
  {
    unspecified = -1,
    fixed_dirichlet,
    controlled_dirichlet,
    mean_zero_multiplier,
    boundary_trace,
    analytic_quadrature_evaluation,
    discrete_cellwise_bounds,
    discrete_facewise_bounds,
    uniform_ellipticity,
    coefficient_regularity,
    coercivity,
    conormal_flux,
    boundary_partition,
    transport_boundary_trace,
    transposition_formulation,
    domain_regularity,
    conforming_trace_subspace,
    fractional_trace_realisation,
    tangential_gradient_realisation,
    target_data_membership,
    metric_realisation
  };

  enum class RequirementStatus
  {
    unspecified = -1,
    provided,
    user_assumed,
    selected_discrete_realisation
  };

  enum class RequirementScope
  {
    unspecified = -1,
    continuous_semantics,
    discrete_compilation,
    both
  };

  enum class NormalOrientation
  {
    unspecified = -1,
    outward,
    inward
  };

  enum class ConormalForm
  {
    unspecified = -1,
    diffusion_minus_transport,
    transport_minus_diffusion
  };

  enum class TraceEvaluationRealisation
  {
    unspecified = -1,
    fe_q_state_trace
  };

  enum class FaceQuadratureRealisation
  {
    unspecified = -1,
    qgauss_face
  };

  enum class TraceWeightRealisation
  {
    unspecified = -1,
    scalar_pointwise_multiplication
  };

  enum class TracePairingRealisation
  {
    unspecified = -1,
    face_quadrature_weights
  };

  enum class TraceTransposeRealisation
  {
    unspecified = -1,
    same_face_quadrature_pullback
  };

  struct TraceRealisationSelection
  {
    std::string                 id;
    std::string                 source_space_id;
    std::string                 output_space_id;
    std::string                 region_id;
    std::string                 weight_data_id;
    std::string                 pairing_id;
    TraceEvaluationRealisation trace_realisation =
      TraceEvaluationRealisation::unspecified;
    FaceQuadratureRealisation  face_quadrature_realisation =
      FaceQuadratureRealisation::unspecified;
    TraceWeightRealisation     weight_realisation =
      TraceWeightRealisation::unspecified;
    TracePairingRealisation    pairing_realisation =
      TracePairingRealisation::unspecified;
    TraceTransposeRealisation  transpose_realisation =
      TraceTransposeRealisation::unspecified;
  };

  enum class Hminus1MetricOperatorRealisation
  {
    unspecified = -1,
    mass_laplacian_inverse_mass
  };

  enum class Hminus1MetricInverseRealisation
  {
    unspecified = -1,
    mass_inverse_laplacian_mass_inverse
  };

  enum class Hminus1MetricNullspaceRealisation
  {
    unspecified = -1,
    fixed_dirichlet_no_nullspace
  };

  struct Hminus1MetricRealisationSelection
  {
    std::string                        id;
    std::string                        metric_id;
    std::string                        primal_space_id;
    std::string                        dual_space_id;
    std::string                        mass_pairing_id;
    std::string                        laplacian_pairing_id;
    std::string                        fixed_boundary_region_id;
    std::string                        laplacian_solve_policy_id;
    std::string                        mass_solve_policy_id;
    Hminus1MetricOperatorRealisation  operator_realisation =
      Hminus1MetricOperatorRealisation::unspecified;
    Hminus1MetricInverseRealisation   inverse_realisation =
      Hminus1MetricInverseRealisation::unspecified;
    Hminus1MetricNullspaceRealisation nullspace_realisation =
      Hminus1MetricNullspaceRealisation::unspecified;
  };

  enum class H1TargetDataRegularityRealisation
  {
    unspecified = -1,
    h1_value_and_weak_gradient
  };

  enum class H1TargetDataTraceRealisation
  {
    unspecified = -1,
    zero_trace_on_fixed_boundary
  };

  struct H1TargetDataMembershipSelection
  {
    std::string                         id;
    std::string                         data_id;
    std::string                         observation_space_id;
    std::string                         fixed_boundary_region_id;
    H1TargetDataRegularityRealisation  regularity_realisation =
      H1TargetDataRegularityRealisation::unspecified;
    H1TargetDataTraceRealisation       trace_realisation =
      H1TargetDataTraceRealisation::unspecified;
  };

  enum class TranspositionOperatorRealisation
  {
    unspecified = -1,
    scalar_diffusion_reaction_dirichlet_laplacian
  };

  enum class TranspositionDiscreteRealisation
  {
    unspecified = -1,
    fe_q_point_sensor_very_weak,
    fe_q_normal_flux_very_weak,
    conforming_nodal_lifting_equivalence
  };

  enum class TranspositionEquivalenceRealisation
  {
    unspecified = -1,
    none,
    conforming_lifting_variational_equivalence
  };

  // Shared strong/very-weak policy for P5.3 observations and the P5.4
  // transposition-control slice. The string ports identify semantic spaces
  // and policies; the enums close the currently registered realizations.
  struct TranspositionRealisationSelection
  {
    std::string                       id;
    std::string                       subject_equation_id;
    std::string                       strong_space_id;
    std::string                       operator_range_space_id;
    std::string                       isomorphism_id;
    std::string                       diffusion_data_id;
    std::string                       reaction_data_id;
    std::string                       residual_codomain_space_id;
    std::string                       multiplier_space_id;
    std::string                       observation_id;
    std::string                       transpose_source_space_id;
    std::string                       domain_regularity_policy_id;
    std::string                       continuous_parent_space_id;
    std::string                       conforming_trace_space_id;
    std::string                       equivalence_policy_id;
    std::string                       conormal_policy_id;
    TranspositionOperatorRealisation operator_realisation =
      TranspositionOperatorRealisation::unspecified;
    TranspositionDiscreteRealisation discrete_realisation =
      TranspositionDiscreteRealisation::unspecified;
    TranspositionEquivalenceRealisation equivalence_realisation =
      TranspositionEquivalenceRealisation::unspecified;
  };

  enum class PartialDirichletInterfaceRealisation
  {
    unspecified = -1,
    fixed_data_precedence
  };

  enum class PartialDirichletTraceRealisation
  {
    unspecified = -1,
    relative_interior_nodal_zero_endpoint
  };

  enum class PartialDirichletHangingRealisation
  {
    unspecified = -1,
    unsupported
  };

  struct PartialDirichletBoundarySelection
  {
    std::string                         id;
    std::string                         subject_id;
    std::string                         transformation_id;
    std::string                         fixed_boundary_region_id;
    std::string                         controlled_boundary_region_id;
    bool                                requires_complete_exterior = true;
    bool                                requires_disjoint_regions = true;
    PartialDirichletInterfaceRealisation interface_realisation =
      PartialDirichletInterfaceRealisation::unspecified;
    PartialDirichletTraceRealisation     trace_realisation =
      PartialDirichletTraceRealisation::unspecified;
    PartialDirichletHangingRealisation   hanging_realisation =
      PartialDirichletHangingRealisation::unspecified;
  };

  enum class FractionalTraceOperatorRealisation
  {
    unspecified = -1,
    volume_mass_plus_stiffness_schur
  };

  enum class FractionalTraceApplyRealisation
  {
    unspecified = -1,
    minimum_h1_extension
  };

  enum class FractionalTraceInverseRealisation
  {
    unspecified = -1,
    full_volume_operator_inverse
  };

  struct FractionalTraceMetricRealisationSelection
  {
    std::string                          id;
    std::string                          metric_id;
    std::string                          control_space_id;
    std::string                          volume_space_id;
    std::string                          trace_inclusion_id;
    std::string                          volume_operator_id;
    std::string                          apply_policy_id;
    std::string                          inverse_policy_id;
    std::string                          solve_policy_id;
    FractionalTraceOperatorRealisation   operator_realisation =
      FractionalTraceOperatorRealisation::unspecified;
    FractionalTraceApplyRealisation      apply_realisation =
      FractionalTraceApplyRealisation::unspecified;
    FractionalTraceInverseRealisation    inverse_realisation =
      FractionalTraceInverseRealisation::unspecified;
  };

  enum class BoundaryH1MetricOperatorRealisation
  {
    unspecified = -1,
    boundary_mass_plus_tangential_stiffness
  };

  enum class BoundaryH1TangentialGradientRealisation
  {
    unspecified = -1,
    projected_ambient_gradient
  };

  enum class BoundaryH1MetricNullspaceRealisation
  {
    unspecified = -1,
    positive_mass_no_nullspace
  };

  struct BoundaryH1MetricRealisationSelection
  {
    std::string                              id;
    std::string                              metric_id;
    std::string                              control_space_id;
    std::string                              boundary_region_id;
    BoundaryH1MetricOperatorRealisation     operator_realisation =
      BoundaryH1MetricOperatorRealisation::unspecified;
    BoundaryH1TangentialGradientRealisation tangential_gradient_realisation =
      BoundaryH1TangentialGradientRealisation::unspecified;
    BoundaryH1MetricNullspaceRealisation    nullspace_realisation =
      BoundaryH1MetricNullspaceRealisation::unspecified;
  };

  // A backend-neutral boundary selection. The first registered general
  // scalar target uses one fixed region, one Robin/outflow region, and
  // explicitly empty Neumann and transport-inflow selections.
  struct BoundaryRealisationSelection
  {
    std::string                 id;
    std::string                 subject_id;
    std::string                 fixed_dirichlet_region_id;
    std::string                 robin_region_id;
    std::vector<std::string>    neumann_region_ids;
    std::vector<std::string>    transport_inflow_region_ids;
    std::string                 transport_outflow_region_id;
    ConormalForm                conormal_form = ConormalForm::unspecified;
    NormalOrientation           normal_orientation =
      NormalOrientation::unspecified;
    TraceEvaluationRealisation trace_realisation =
      TraceEvaluationRealisation::unspecified;
    FaceQuadratureRealisation  face_quadrature_realisation =
      FaceQuadratureRealisation::unspecified;
  };

  enum class FormulationKind
  {
    unspecified = -1,
    reduced_dto,
    all_at_once
  };

  struct RegionSpec
  {
    std::string               id;
    std::string               label;
    RegionKind                kind = RegionKind::unspecified;
    bool                      is_full_domain = false;
    std::vector<unsigned int> boundary_ids;
    std::vector<unsigned int> material_ids;
    // Backend-neutral physical coordinates for a point-set region. The
    // compiler validates their dimension against its concrete mesh dimension
    // and copies them into its immutable sensor operator.
    std::vector<std::vector<double>> point_coordinates = {};
  };

  struct SpaceSpec
  {
    std::string   id;
    std::string   label;
    std::string   region_id;
    SpaceTopology topology = SpaceTopology::unspecified;
    SpaceRole     role = SpaceRole::unspecified;
    bool          is_scalar = true;
    // A positive dimension is required for finite-dimensional point-sensor
    // observation spaces. Other semantic spaces leave this unspecified.
    std::size_t dimension = 0;
  };

  struct PairingSpec
  {
    std::string id;
    std::string label;
    std::string primal_space_id;
    std::string covector_space_id;
  };

  struct VariableSpec
  {
    std::string  id;
    std::string  label;
    VariableRole role = VariableRole::unspecified;
    std::string  space_id;
    // An empty port means that this variable is already its physical field.
    // Otherwise the referenced transformation reconstructs that field before
    // residual and observation evaluation.
    std::string  physical_field_transform_id;
  };

  struct DataSpec
  {
    std::string id;
    std::string label;
    DataKind    kind = DataKind::unspecified;
    DataRole    role = DataRole::unspecified;
    std::string space_id;
  };

  struct TransformationSpec
  {
    std::string        id;
    std::string        label;
    TransformationKind kind = TransformationKind::unspecified;
    std::string        input_variable_id;
    std::string        output_space_id;
    std::string        fixed_data_id;
    // A controlled Dirichlet lifting has a second primal input. Fixed-data
    // reconstructions leave this port empty.
    std::string        control_variable_id;
  };

  struct ResidualTermSpec
  {
    std::string              id;
    std::string              label;
    ResidualTermKind         kind = ResidualTermKind::unspecified;
    std::string              equation_id;
    std::vector<std::string> variable_ids;
    std::vector<std::string> data_ids = {};
    // Natural boundary terms name the boundary on which their declared trace
    // pairing is evaluated. Volume terms leave this port empty.
    std::string              region_id;
  };

  struct EquationBlockSpec
  {
    std::string              id;
    std::string              label;
    std::string              test_space_id;
    std::string              test_pairing_id;
    std::vector<std::string> residual_term_ids;
  };

  struct ObservationSpec
  {
    std::string     id;
    std::string     label;
    ObservationKind kind = ObservationKind::unspecified;
    std::string     input_variable_id;
    std::string     region_id;
    std::string     output_space_id;
    std::string     output_pairing_id;
    // Immutable inputs consumed by the observation map. Plain restrictions
    // leave this port empty; weighted observations name their fixed data.
    std::vector<std::string> data_ids;
  };

  struct LossSpec
  {
    std::string id;
    std::string label;
    LossKind    kind = LossKind::unspecified;
    std::string source_observation_id;
    std::string data_id;
    std::string pairing_id;
  };

  struct MetricSpec
  {
    std::string id;
    std::string label;
    MetricKind  kind = MetricKind::unspecified;
    std::string variable_id;
    std::string pairing_id;
  };

  struct ConstraintSpec
  {
    std::string    id;
    std::string    label;
    ConstraintKind kind = ConstraintKind::unspecified;
    std::string    variable_id;
    std::string    lower_bound_data_id;
    std::string    upper_bound_data_id;
  };

  struct RequirementPolicySpec
  {
    RequirementPolicySpec() = default;

    RequirementPolicySpec(std::string                            policy_id,
                          std::string                            policy_subject_id,
                          const RequirementKind                  policy_kind,
                          const RequirementStatus                policy_status,
                          const RequirementScope                 policy_scope,
                          std::string                            policy_description,
                          std::string                            policy_region_id,
                          std::optional<BoundaryRealisationSelection>
                            policy_typed_selection = std::nullopt)
      : id(std::move(policy_id))
      , subject_id(std::move(policy_subject_id))
      , kind(policy_kind)
      , status(policy_status)
      , scope(policy_scope)
      , selected_policy(std::move(policy_description))
      , region_id(std::move(policy_region_id))
      , typed_selection(std::move(policy_typed_selection))
    {}

    std::string       id;
    std::string       subject_id;
    RequirementKind   kind = RequirementKind::unspecified;
    RequirementStatus status = RequirementStatus::unspecified;
    RequirementScope  scope = RequirementScope::unspecified;
    std::string       selected_policy;
    std::string       region_id;
    std::optional<BoundaryRealisationSelection> typed_selection;
    std::optional<TraceRealisationSelection>     typed_trace_selection;
    std::optional<Hminus1MetricRealisationSelection>
      typed_metric_selection;
    std::optional<TranspositionRealisationSelection>
      typed_transposition_selection;
    std::optional<PartialDirichletBoundarySelection>
      typed_partial_boundary_selection;
    std::optional<FractionalTraceMetricRealisationSelection>
      typed_fractional_metric_selection;
    std::optional<BoundaryH1MetricRealisationSelection>
      typed_boundary_h1_metric_selection;
    std::optional<H1TargetDataMembershipSelection>
      typed_h1_target_data_membership_selection;
  };

  struct ReducedFormulationSpec
  {
    std::string     id;
    FormulationKind kind = FormulationKind::unspecified;
    std::string     state_variable_id;
    // This is the binary reduced DTO decision port. It may identify either a
    // control or a parameter variable; the field name is retained for source
    // compatibility with the first control-only semantic slice.
    std::string     control_variable_id;
    std::string     equation_id;
    std::string     metric_id;
    std::string     constraint_id;
  };

  // This is a composition root, not a PDE-model class. The v1 compiler only
  // accepts the explicitly represented narrow stationary volume-control
  // graph; extensions add semantic node kinds and registered lowerers.
  struct ProblemSpec
  {
    std::string                        id;
    std::string                        label;
    std::vector<RegionSpec>            regions;
    std::vector<SpaceSpec>             spaces;
    std::vector<PairingSpec>           pairings;
    std::vector<VariableSpec>          variables;
    std::vector<DataSpec>              data;
    std::vector<TransformationSpec>    transformations;
    std::vector<ResidualTermSpec>      residual_terms;
    std::vector<EquationBlockSpec>     equations;
    std::vector<ObservationSpec>       observations;
    std::vector<LossSpec>              losses;
    std::vector<MetricSpec>            metrics;
    std::vector<ConstraintSpec>        constraints;
    std::vector<RequirementPolicySpec> requirement_policies;
    ReducedFormulationSpec              formulation;
  };

  enum class DiagnosticCategory
  {
    structural,
    analytical_policy,
    lowerability,
    formulation_capability
  };

  struct Diagnostic
  {
    DiagnosticCategory category;
    std::string        component_id;
    std::string        capability;
    std::string        remedy;
  };
} // namespace nmopt::semantic::v1
