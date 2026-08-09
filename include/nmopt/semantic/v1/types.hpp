#pragma once

#include <string>
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
    boundary
  };

  enum class SpaceTopology
  {
    unspecified = -1,
    h1,
    l2
  };

  enum class SpaceRole
  {
    unspecified = -1,
    state,
    test,
    control,
    parameter,
    observation,
    data
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
    weighted_boundary_trace
  };

  enum class LossKind
  {
    unspecified = -1,
    quadratic_tracking,
    quadratic_control_regularisation,
    quadratic_h1_control_regularisation,
    quadratic_parameter_regularisation
  };

  enum class MetricKind
  {
    unspecified = -1,
    l2,
    h1
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
    transport_boundary_trace
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
  };

  struct SpaceSpec
  {
    std::string   id;
    std::string   label;
    std::string   region_id;
    SpaceTopology topology = SpaceTopology::unspecified;
    SpaceRole     role = SpaceRole::unspecified;
    bool          is_scalar = true;
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
    std::string       id;
    std::string       subject_id;
    RequirementKind   kind = RequirementKind::unspecified;
    RequirementStatus status = RequirementStatus::unspecified;
    RequirementScope  scope = RequirementScope::unspecified;
    std::string       selected_policy;
    std::string       region_id;
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
