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
    volume,
    boundary
  };

  enum class SpaceTopology
  {
    h1,
    l2
  };

  enum class SpaceRole
  {
    state,
    test,
    control,
    observation,
    data
  };

  enum class VariableRole
  {
    state,
    control
  };

  enum class DataKind
  {
    function,
    scalar_constant,
    cellwise_bound,
    facewise_bound
  };

  enum class DataRole
  {
    forcing,
    desired_state,
    fixed_dirichlet_lifting,
    diffusion,
    reaction,
    regularisation_weight,
    lower_bound,
    upper_bound
  };

  enum class ResidualTermKind
  {
    diffusion_reaction,
    volume_source,
    volume_control,
    neumann_control
  };

  enum class ObservationKind
  {
    volume_restriction,
    boundary_trace,
    boundary_restriction
  };

  enum class LossKind
  {
    quadratic_tracking,
    quadratic_control_regularisation,
    quadratic_h1_control_regularisation
  };

  enum class MetricKind
  {
    l2
  };

  enum class ConstraintKind
  {
    cellwise_box,
    facewise_box
  };

  enum class TransformationKind
  {
    fixed_dirichlet_reconstruction
  };

  enum class RequirementKind
  {
    fixed_dirichlet,
    mean_zero_multiplier,
    boundary_trace,
    analytic_quadrature_evaluation,
    discrete_cellwise_bounds,
    discrete_facewise_bounds
  };

  enum class RequirementStatus
  {
    provided,
    user_assumed,
    selected_discrete_realisation
  };

  enum class RequirementScope
  {
    continuous_semantics,
    discrete_compilation,
    both
  };

  enum class FormulationKind
  {
    reduced_dto,
    all_at_once
  };

  struct RegionSpec
  {
    std::string               id;
    std::string               label;
    RegionKind                kind;
    bool                      is_full_domain = false;
    std::vector<unsigned int> boundary_ids;
    std::vector<unsigned int> material_ids;
  };

  struct SpaceSpec
  {
    std::string   id;
    std::string   label;
    std::string   region_id;
    SpaceTopology topology;
    SpaceRole     role;
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
    VariableRole role;
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
    DataKind    kind;
    DataRole    role;
    std::string space_id;
  };

  struct TransformationSpec
  {
    std::string        id;
    std::string        label;
    TransformationKind kind;
    std::string        input_variable_id;
    std::string        output_space_id;
    std::string        fixed_data_id;
  };

  struct ResidualTermSpec
  {
    std::string              id;
    std::string              label;
    ResidualTermKind         kind;
    std::string              equation_id;
    std::vector<std::string> variable_ids;
    std::vector<std::string> data_ids;
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
    ObservationKind kind;
    std::string     input_variable_id;
    std::string     region_id;
    std::string     output_space_id;
    std::string     output_pairing_id;
  };

  struct LossSpec
  {
    std::string id;
    std::string label;
    LossKind    kind;
    std::string source_observation_id;
    std::string data_id;
    std::string pairing_id;
  };

  struct MetricSpec
  {
    std::string id;
    std::string label;
    MetricKind  kind;
    std::string variable_id;
    std::string pairing_id;
  };

  struct ConstraintSpec
  {
    std::string    id;
    std::string    label;
    ConstraintKind kind;
    std::string    variable_id;
    std::string    lower_bound_data_id;
    std::string    upper_bound_data_id;
  };

  struct RequirementPolicySpec
  {
    std::string       id;
    std::string       subject_id;
    RequirementKind   kind;
    RequirementStatus status;
    RequirementScope  scope;
    std::string       selected_policy;
    std::string       region_id;
  };

  struct ReducedFormulationSpec
  {
    std::string     id;
    FormulationKind kind;
    std::string     state_variable_id;
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
