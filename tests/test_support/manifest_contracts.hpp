#pragma once

#include "nmopt/compiler/v1/compiled_problem.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

namespace nmopt::test_support
{
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
