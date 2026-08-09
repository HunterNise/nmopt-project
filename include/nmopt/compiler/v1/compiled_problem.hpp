#pragma once

#include "nmopt/contract/reduced_dto.hpp"
#include "nmopt/semantic/v1/validation.hpp"

#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace nmopt::compiler::v1
{
  // This is descriptive provenance, not a second executable configuration.
  // It records the selected compilation choices that affect the discrete
  // model, so a later DTO/OTD or assembled/matrix-free result cannot be
  // mistaken for the same computation.
  struct CompilationManifest
  {
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

    CompiledProblemT(std::shared_ptr<const Model>             executable,
                     std::shared_ptr<const Metric>            metric,
                     std::shared_ptr<const Constraint>        constraint,
                     contract::StateAdjointSolversT<Backend>   solvers,
                     CompilationManifest                       manifest,
                     std::shared_ptr<const void>               lifetime_owner = {})
      : executable_(std::move(executable))
      , metric_(std::move(metric))
      , constraint_(std::move(constraint))
      , solvers_(std::move(solvers))
      , manifest_(std::move(manifest))
      , lifetime_owner_(std::move(lifetime_owner))
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
  };

  template <typename Backend>
  struct CompilationResultT
  {
    semantic::v1::ValidationReport             diagnostics;
    std::shared_ptr<const CompiledProblemT<Backend>> problem;

    bool
    succeeded() const
    {
      return diagnostics.valid() && static_cast<bool>(problem);
    }
  };
} // namespace nmopt::compiler::v1
