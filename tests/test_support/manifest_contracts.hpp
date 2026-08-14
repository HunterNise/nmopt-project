#pragma once

#include "nmopt/compiler/v1/compiled_problem.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

namespace nmopt::test_support
{
  inline void
  require_manifest_compatibility_equal(
    const compiler::v1::CompilationManifest &left,
    const compiler::v1::CompilationManifest &right,
    const std::string &                        description)
  {
    const auto &lhs = left.resolved_decision.compatibility;
    const auto &rhs = right.resolved_decision.compatibility;
    contract::require(
      lhs.compiler_id == rhs.compiler_id && lhs.backend == rhs.backend &&
        lhs.execution == rhs.execution && lhs.state_space == rhs.state_space &&
        lhs.control_space == rhs.control_space &&
        lhs.quadrature == rhs.quadrature &&
        lhs.dual_representation == rhs.dual_representation &&
        lhs.data_rule == rhs.data_rule &&
        lhs.observation_realisation == rhs.observation_realisation &&
        lhs.metric_solve_policy == rhs.metric_solve_policy &&
        lhs.constraint_realisation == rhs.constraint_realisation &&
        lhs.lifting_realisation == rhs.lifting_realisation &&
        lhs.nullspace_policy == rhs.nullspace_policy &&
        lhs.state_adjoint_solve_policy == rhs.state_adjoint_solve_policy &&
        lhs.provenance == rhs.provenance &&
        lhs.lowering_handler_records == rhs.lowering_handler_records &&
        lhs.region_ids == rhs.region_ids && lhs.space_ids == rhs.space_ids &&
        lhs.pairing_ids == rhs.pairing_ids &&
        lhs.variable_ids == rhs.variable_ids && lhs.data_ids == rhs.data_ids &&
        lhs.transformation_ids == rhs.transformation_ids &&
        lhs.residual_term_ids == rhs.residual_term_ids &&
        lhs.observation_ids == rhs.observation_ids &&
        lhs.loss_ids == rhs.loss_ids && lhs.metric_ids == rhs.metric_ids &&
        lhs.constraint_ids == rhs.constraint_ids &&
        lhs.declared_assumptions == rhs.declared_assumptions,
      description + " changed compatibility provenance");
  }

  inline void
  require_dirichlet_manifest_dimensions(
    const compiler::v1::CompilationManifest &manifest,
    const std::size_t                         independent_state_dimension,
    const std::size_t                         test_dimension,
    const std::size_t                         control_dimension,
    const std::size_t                         physical_state_dimension,
    const std::string &                        description)
  {
    const auto dimension = [&manifest, &description](const std::string &id) {
      const auto space = std::find_if(
        manifest.spaces.begin(),
        manifest.spaces.end(),
        [&id](const compiler::v1::CompiledSpaceRecord &candidate) {
          return candidate.semantic_id == id;
        });
      contract::require(space != manifest.spaces.end(),
                        description + " omitted manifest space " + id);
      return space->dimension;
    };

    contract::require(
      dimension("state_space") == independent_state_dimension,
      description + " reported the wrong independent state dimension");
    contract::require(dimension("state_test_space") == test_dimension,
                      description + " reported the wrong state-test dimension");
    contract::require(dimension("control_space") == control_dimension,
                      description + " reported the wrong control dimension");
    contract::require(
      dimension("state_observation_space") == physical_state_dimension,
      description + " reported the wrong physical state-observation dimension");
    contract::require(
      dimension("control_observation_space") == control_dimension,
      description + " reported the wrong control-observation dimension");
    contract::require(
      physical_state_dimension != independent_state_dimension,
      description + " did not exercise distinct physical and independent state layouts");
  }
} // namespace nmopt::test_support
