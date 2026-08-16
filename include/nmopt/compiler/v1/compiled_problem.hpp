#pragma once

#include "nmopt/contract/reduced_dto.hpp"
#include "nmopt/contract/reduced_hessian.hpp"
#include "nmopt/contract/supplied_otd.hpp"
#include "nmopt/semantic/v1/validation.hpp"

#include <memory>
#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace nmopt::compiler::v1
{
  enum class MeshLifetimePolicy
  {
    borrowed_immutable,
    owned_session
  };

  enum class ExecutionRealisation
  {
    assembled
  };

  enum class LinearSolveAlgorithm
  {
    serial_cg,
    serial_sparse_direct_umfpack
  };

  enum class CompiledFieldShape
  {
    unspecified,
    scalar,
    vector,
    tensor,
    scalar_constant,
    cellwise_scalar,
    facewise_scalar
  };

  enum class CompiledBindingStatus
  {
    checked,
    provenance_only
  };

  struct CompiledSpaceRecord
  {
    std::string             semantic_id;
    semantic::v1::SpaceRole role = semantic::v1::SpaceRole::unspecified;
    std::string             runtime_role;
    std::string             region_id;
    std::string             finite_element;
    std::size_t             dimension = 0;
  };

  struct CompiledMeshRecord
  {
    unsigned int       dimension = 0;
    std::size_t        active_cells = 0;
    std::string        provenance;
    MeshLifetimePolicy lifetime = MeshLifetimePolicy::borrowed_immutable;
    // Stable identity of the mesh structure consumed by lowering.  This is
    // deliberately separate from caller provenance: two meshes can have the
    // same source label while producing different discrete operators.
    std::string        structural_identity;
  };

  struct CompiledBindingRecord
  {
    std::string            semantic_id;
    semantic::v1::DataRole role = semantic::v1::DataRole::unspecified;
    semantic::v1::DataKind kind = semantic::v1::DataKind::unspecified;
    std::string            space_id;
    std::string            region_id;
    std::string            representation;
    std::string            evaluation_realisation;
    std::string            provenance;
    CompiledFieldShape     field_shape = CompiledFieldShape::unspecified;
    std::string            runtime_representation;
    std::optional<double>  scalar_value;
    std::string            value_digest;
    CompiledBindingStatus  value_status = CompiledBindingStatus::provenance_only;
  };

  struct CompiledSolvePolicyRecord
  {
    LinearSolveAlgorithm algorithm = LinearSolveAlgorithm::serial_cg;
    std::string          preconditioner;
    unsigned int         maximum_iterations = 0;
    double               relative_tolerance = 0.0;
    double               absolute_tolerance = 0.0;
    std::string          nullspace_policy;
    // Typed description of the realized operator invocation.
    std::string          operator_realisation;
  };

  struct CompiledFormulationRecord
  {
    std::string                         semantic_id;
    semantic::v1::FormulationKind       kind =
      semantic::v1::FormulationKind::unspecified;
    semantic::v1::FormulationProvenance provenance =
      semantic::v1::FormulationProvenance::unspecified;
    ExecutionRealisation execution = ExecutionRealisation::assembled;
    std::string          dual_representation;
  };

  // Typed provenance for the separately executable supplied-OTD product.
  // The callbacks themselves live in SuppliedOTDSystem; this record describes
  // the block boundary and the evidence attached to that product.
  struct CompiledSuppliedOTDRecord
  {
    bool                     present = false;
    std::optional<semantic::v1::SuppliedOTDDeclaration> declaration;
    std::string              variable_layout;
    std::string              residual_layout;
    std::vector<std::string> variable_space_ids;
    std::vector<std::size_t> variable_dimensions;
    std::vector<std::string> residual_space_ids;
    std::vector<std::size_t> residual_dimensions;
    std::size_t              state_variable_block = 0;
    std::size_t              adjoint_variable_block = 0;
    std::size_t              control_variable_block = 0;
    std::size_t              state_equation_block = 0;
    std::size_t              adjoint_equation_block = 0;
    std::size_t              control_stationarity_block = 0;
    std::string              sign_convention;
    std::string              discretisation_provenance;
    std::string              state_block_provenance;
    std::string              adjoint_block_provenance;
    std::string              stationarity_block_provenance;
    std::string              value_action_provenance;
    std::string              jvp_action_provenance;
    std::string              vjp_action_provenance;
    std::string              solve_provenance;
    std::string              comparison_status;
  };

  struct CompiledMetricRecord
  {
    std::string                 semantic_id;
    std::string                 realisation_id;
    std::string                 operator_description;
    std::string                 primal_space_id;
    std::string                 dual_space_id;
    std::string                 operator_id;
    std::string                 inverse_operator_id;
    std::string                 mass_pairing_id;
    std::string                 laplacian_pairing_id;
    std::string                 boundary_region_id;
    std::string                 laplacian_solve_policy_id;
    std::string                 mass_solve_policy_id;
    std::string                 nullspace_realisation_id;
    CompiledSolvePolicyRecord solve_policy;
  };

  struct CompiledConstraintRecord
  {
    bool        present = false;
    std::string semantic_id;
    std::string realisation_id;
    std::string projection_metric_id;
  };

  struct CompiledRealisationRecord
  {
    std::string              semantic_id;
    std::string              kind;
    std::string              realisation_id;
    std::string              handler_id;
    std::vector<std::string> input_ids;
    std::vector<std::string> output_ids;
    std::string              region_id;
    semantic::v1::ResidualTermKind residual_kind =
      semantic::v1::ResidualTermKind::unspecified;
    semantic::v1::ObservationKind observation_kind =
      semantic::v1::ObservationKind::unspecified;
    semantic::v1::LossKind loss_kind = semantic::v1::LossKind::unspecified;
    semantic::v1::TransformationKind transformation_kind =
      semantic::v1::TransformationKind::unspecified;
  };

  struct CompiledPairingRecord
  {
    std::string pairing_id;
    std::string primal_space_id;
    std::string covector_space_id;
  };

  struct CompiledRegionRecord
  {
    std::string               semantic_id;
    semantic::v1::RegionKind kind = semantic::v1::RegionKind::unspecified;
    bool                      is_full_domain = false;
    std::vector<unsigned int> boundary_ids;
    std::vector<unsigned int> material_ids;
    std::size_t               point_count = 0;
  };

  struct CompiledAssumptionRecord
  {
    std::string                    id;
    std::string                    subject_id;
    semantic::v1::RequirementKind kind =
      semantic::v1::RequirementKind::unspecified;
    semantic::v1::RequirementStatus status =
      semantic::v1::RequirementStatus::unspecified;
    semantic::v1::RequirementScope scope =
      semantic::v1::RequirementScope::unspecified;
    std::string region_id;
    // Only a model-author assumption retains free-form text. Selected
    // compiler capabilities are represented by the typed records below.
    std::string model_author_declaration;
  };

  // Compatibility text is a rendered view of the typed decision.  It is kept
  // in the decision so manifest construction cannot reconstruct it from the
  // semantic graph or a target enum.
  struct CompiledCompatibilityView
  {
    std::string              compiler_id;
    std::string              backend;
    std::string              execution;
    std::string              state_space;
    std::string              control_space;
    std::string              quadrature;
    std::string              dual_representation;
    std::string              data_rule;
    std::string              observation_realisation;
    std::string              metric_solve_policy;
    std::string              constraint_realisation;
    std::string              lifting_realisation;
    std::string              nullspace_policy;
    std::string              state_adjoint_solve_policy;
    std::string              provenance;
    std::vector<std::string> lowering_handler_records;
    std::vector<std::string> region_ids;
    std::vector<std::string> space_ids;
    std::vector<std::string> pairing_ids;
    std::vector<std::string> variable_ids;
    std::vector<std::string> data_ids;
    std::vector<std::string> transformation_ids;
    std::vector<std::string> residual_term_ids;
    std::vector<std::string> observation_ids;
    std::vector<std::string> loss_ids;
    std::vector<std::string> metric_ids;
    std::vector<std::string> constraint_ids;
    std::vector<std::string> declared_assumptions;
  };

  struct CompiledRealizedSpaceRecord
  {
    std::string map_id;
    std::string semantic_id;
    std::string realization_id;
    std::size_t dimension = 0;
    std::string layout;
    std::string ordering;
    std::string pairing_id;
  };

  struct CompiledRealizedMapRecord
  {
    std::string              semantic_id;
    std::vector<std::string> input_space_ids;
    std::string              source_space_id;
    std::string              output_space_id;
    std::vector<std::size_t> input_dimensions;
    std::size_t              source_dimension = 0;
    std::size_t              output_dimension = 0;
    std::string              realization_id;
    std::string              source_layout;
    std::string              output_layout;
    std::string              ordering;
    std::string              pairing_realization;
    std::string              transformation_chain;
    std::string              value_provenance;
    std::string              jvp_provenance;
    std::string              vjp_provenance;
    std::string              pairing_id;
  };

  // The single typed decision selected by semantic resolution and compiler
  // lowerability checks.  Manifest display strings are rendered from this
  // record and retained only as a compatibility view.
  struct ResolvedCompilationDecision
  {
    std::string                         semantic_problem_id;
    std::string                         formulation_id;
    std::string                         target_id;
    std::string                         execution_id = "assembled";
    CompiledFormulationRecord           formulation_record;
    CompiledSuppliedOTDRecord           supplied_otd_record;
    CompiledMeshRecord                  mesh_record;
    std::vector<CompiledRegionRecord>   regions;
    std::vector<CompiledSpaceRecord>    spaces;
    std::vector<CompiledBindingRecord>  bindings;
    std::vector<CompiledPairingRecord>  pairings;
    std::vector<CompiledRealisationRecord> residuals;
    std::vector<CompiledRealisationRecord> observations;
    std::vector<CompiledRealisationRecord> losses;
    std::vector<CompiledRealisationRecord> transformations;
    std::vector<CompiledRealizedSpaceRecord> realized_spaces;
    std::vector<CompiledRealizedMapRecord>   realized_maps;
    CompiledSolvePolicyRecord             state_solve_record;
    CompiledSolvePolicyRecord             adjoint_solve_record;
    CompiledMetricRecord                  metric_record;
    CompiledConstraintRecord              constraint_record;
    CompiledCompatibilityView             compatibility;
    std::optional<semantic::v1::BoundaryRealisationSelection>
      boundary_realisation;
    std::optional<semantic::v1::TranspositionRealisationSelection>
      transposition_realisation;
    std::optional<semantic::v1::PartialDirichletBoundarySelection>
      partial_boundary_selection;
    std::optional<semantic::v1::FractionalTraceMetricRealisationSelection>
      fractional_metric_selection;
    std::optional<semantic::v1::BoundaryH1MetricRealisationSelection>
      boundary_h1_metric_selection;
    std::optional<semantic::v1::H1TargetDataMembershipSelection>
      h1_target_data_membership_selection;
    std::vector<CompiledAssumptionRecord> assumptions;
  };

  // This is descriptive provenance, not a second executable configuration.
  // It records the selected compilation choices that affect the discrete
  // model, so a later DTO/OTD or assembled/matrix-free result cannot be
  // mistaken for the same computation.
  struct CompilationManifest
  {
    unsigned int                         schema_version = 3;
    ResolvedCompilationDecision           resolved_decision;
    CompiledFormulationRecord            formulation_record;
    CompiledSuppliedOTDRecord            supplied_otd_record;
    CompiledMeshRecord                   mesh_record;
    std::vector<CompiledSpaceRecord>      spaces;
    std::vector<CompiledBindingRecord>    bindings;
    CompiledSolvePolicyRecord             state_solve_record;
    CompiledSolvePolicyRecord             adjoint_solve_record;
    CompiledMetricRecord                  metric_record;
    CompiledConstraintRecord              constraint_record;
    std::vector<CompiledRealizedSpaceRecord> realized_spaces;
    std::vector<CompiledRealizedMapRecord>   realized_maps;
    std::optional<semantic::v1::BoundaryRealisationSelection>
                                          boundary_realisation;
    std::optional<semantic::v1::TranspositionRealisationSelection>
                                          transposition_realisation;
    std::optional<semantic::v1::PartialDirichletBoundarySelection>
                                          partial_boundary_selection;
    std::optional<semantic::v1::FractionalTraceMetricRealisationSelection>
                                          fractional_metric_selection;
    std::optional<semantic::v1::BoundaryH1MetricRealisationSelection>
                                          boundary_h1_metric_selection;
    std::optional<semantic::v1::H1TargetDataMembershipSelection>
                                          h1_target_data_membership_selection;
    std::vector<std::string>              lowering_handler_records;
    // Human-readable rendering retained for logs and source compatibility.
    // Tests and experiment tooling use the structured records above.
    std::string              semantic_problem_id;
    std::string              compiler_id;
    std::string              backend;
    std::string              execution;
    std::string              state_space;
    std::string              control_space;
    std::string              quadrature;
    std::string              dual_representation;
    std::string              data_rule;
    std::string              observation_realisation;
    std::string              metric_solve_policy;
    std::string              constraint_realisation;
    std::string              lifting_realisation;
    std::string              nullspace_policy;
    std::string              state_adjoint_solve_policy;
    std::string              provenance;
    std::vector<std::string> region_ids;
    std::vector<std::string> space_ids;
    std::vector<std::string> pairing_ids;
    std::vector<std::string> variable_ids;
    std::vector<std::string> data_ids;
    std::vector<std::string> transformation_ids;
    std::vector<std::string> residual_term_ids;
    std::vector<std::string> observation_ids;
    std::vector<std::string> loss_ids;
    std::vector<std::string> metric_ids;
    std::vector<std::string> constraint_ids;
    std::vector<std::string> declared_assumptions;
  };

  // The compiler product contains only the backend-neutral executable ports.
  // A concrete lowerer stays private to its compiler and supplies the model,
  // metric, optional constraint, and formulation services below.
  template <typename Backend>
  class CompiledProblemT final
  {
  public:
    using Model = contract::ExecutableModelT<Backend>;
    using Metric = contract::MetricT<Backend>;
    using Constraint = contract::ConstraintT<Backend>;
    using ReducedHessian = contract::ReducedHessianT<Backend>;

    CompiledProblemT(std::shared_ptr<const Model>             executable,
                     std::shared_ptr<const Metric>            metric,
                     std::shared_ptr<const Constraint>        constraint,
                     contract::StateAdjointSolversT<Backend>   solvers,
                     CompilationManifest                       manifest,
                     std::shared_ptr<const void>               lifetime_owner = {},
                     std::shared_ptr<const ReducedHessian>     reduced_hessian = {})
      : executable_(std::move(executable))
      , metric_(std::move(metric))
      , constraint_(std::move(constraint))
      , solvers_(std::move(solvers))
      , manifest_(std::move(manifest))
      , lifetime_owner_(std::move(lifetime_owner))
      , reduced_hessian_(std::move(reduced_hessian))
    {
      contract::require(static_cast<bool>(executable_),
                        "A compiled problem needs an executable model");
      contract::require(static_cast<bool>(metric_),
                        "A compiled problem needs a search metric");
      contract::require(static_cast<bool>(solvers_.solve_state),
                        "A compiled problem needs a state solve service");
      contract::require(static_cast<bool>(solvers_.solve_adjoint),
                        "A compiled problem needs an adjoint solve service");
      contract::require(!manifest_.semantic_problem_id.empty(),
                        "A compiled problem manifest needs a semantic identifier");
      const contract::StateControlPartitionT<Backend> partition(*executable_, 0, 1);
      contract::require(metric_->layout()->compatible_with(
                          *partition.control_layout()),
                        "A compiled problem metric does not match its control block");
      if (constraint_)
        {
          contract::require(constraint_->layout()->compatible_with(*metric_->layout()),
                            "A compiled problem constraint does not match its metric");
          contract::require(constraint_->supports_projection_in(*metric_),
                            "A compiled problem constraint cannot project in its metric");
        }
      if (reduced_hessian_)
        contract::require(
          reduced_hessian_->layout()->compatible_with(*metric_->layout()),
          "A compiled problem Hessian does not match its metric");
    }

    const Model &
    executable_model() const
    {
      return *executable_;
    }

    const Metric &
    metric() const
    {
      return *metric_;
    }

    const Constraint *
    constraint() const
    {
      return constraint_ ? constraint_.get() : nullptr;
    }

    const ReducedHessian *
    reduced_hessian() const
    {
      return reduced_hessian_ ? reduced_hessian_.get() : nullptr;
    }

    const contract::StateAdjointSolversT<Backend> &
    state_adjoint_solvers() const
    {
      return solvers_;
    }

    contract::ReducedDTOT<Backend>
    make_reduced_dto() const
    {
      return contract::ReducedDTOT<Backend>(
        executable_,
        contract::StateControlPartitionT<Backend>(*executable_, 0, 1),
        solvers_,
        lifetime_owner_);
    }

    const CompilationManifest &
    manifest() const
    {
      return manifest_;
    }

  private:
    std::shared_ptr<const Model>           executable_;
    std::shared_ptr<const Metric>          metric_;
    std::shared_ptr<const Constraint>      constraint_;
    contract::StateAdjointSolversT<Backend> solvers_;
    CompilationManifest                     manifest_;
    std::shared_ptr<const void>             lifetime_owner_;
    std::shared_ptr<const ReducedHessian>   reduced_hessian_;
  };

  // A supplied-OTD compilation is a distinct product.  It intentionally does
  // not expose the reduced DTO service bundle used by CompiledProblemT.
  template <typename Backend>
  class CompiledSuppliedOTDProblemT final
  {
  public:
    using System = contract::SuppliedOTDSystemT<Backend>;

    CompiledSuppliedOTDProblemT(std::shared_ptr<const System> system,
                                CompilationManifest           manifest,
                                std::shared_ptr<const void>    lifetime_owner = {})
      : system_(std::move(system))
      , manifest_(std::move(manifest))
      , lifetime_owner_(std::move(lifetime_owner))
    {
      contract::require(static_cast<bool>(system_),
                        "A supplied-OTD problem needs an executable system");
      contract::require(
        manifest_.formulation_record.kind ==
            semantic::v1::FormulationKind::all_at_once &&
          manifest_.formulation_record.provenance ==
            semantic::v1::FormulationProvenance::supplied_otd,
        "A supplied-OTD problem needs an all-at-once supplied-OTD manifest");
      contract::require(manifest_.supplied_otd_record.present,
                        "A supplied-OTD problem needs typed block provenance");
      contract::require(
        manifest_.supplied_otd_record.declaration.has_value(),
        "A supplied-OTD problem needs its typed formulation declaration");
    }

    const System &
    system() const
    {
      return *system_;
    }

    const CompilationManifest &
    manifest() const
    {
      return manifest_;
    }

  private:
    std::shared_ptr<const System> system_;
    CompilationManifest           manifest_;
    std::shared_ptr<const void>    lifetime_owner_;
  };

  template <typename Backend>
  struct CompilationResultT
  {
    semantic::v1::ValidationReport             diagnostics;
    std::shared_ptr<const CompiledProblemT<Backend>> problem;
    std::shared_ptr<const CompiledSuppliedOTDProblemT<Backend>>
      supplied_otd_problem;

    bool
    succeeded() const
    {
      return diagnostics.valid() &&
             (static_cast<bool>(problem) ||
              static_cast<bool>(supplied_otd_problem));
    }
  };
} // namespace nmopt::compiler::v1
