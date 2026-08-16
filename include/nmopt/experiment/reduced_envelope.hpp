#pragma once

#include "nmopt/compiler/v1/compiled_problem.hpp"
#include "nmopt/solvers/reduced_gradient.hpp"
#include "nmopt/solvers/reduced_trust_region.hpp"

#include <string>
#include <utility>

namespace nmopt::experiment
{
  // The environment is deliberately supplied by the caller. Collecting host
  // or timing data is an orchestration concern, not a responsibility of the
  // backend-neutral solver or this in-memory association contract.
  struct RunEnvironmentRecord
  {
    std::string source_revision;
    std::string build_profile;
    std::string compiler;
    std::string compiler_version;
    std::string standard_library;
    std::string operating_system;
    std::string architecture;
    std::string hardware;
  };

  struct ReducedSearchPolicySnapshot
  {
    std::string                              solver_name;
    std::string                              policy_name;
    solvers::ReducedSolverParameters         solver_parameters;
    solvers::ReducedLineSearchPolicySnapshot line_search_parameters;
  };

  struct ReducedTrustRegionPolicySnapshot
  {
    std::string                         solver_name;
    std::string                         policy_name;
    solvers::ReducedTrustRegionParameters parameters;
  };

  template <typename Backend>
  inline ReducedSearchPolicySnapshot
  make_reduced_search_policy_snapshot(
    const solvers::ReducedSolverResultT<Backend> &report)
  {
    return {"reduced_search",
            report.policy_name,
            report.parameters,
            report.policy_parameters};
  }

  template <typename Backend>
  inline ReducedTrustRegionPolicySnapshot
  make_reduced_trust_region_policy_snapshot(
    const solvers::ReducedTrustRegionResultT<Backend> &report)
  {
    return {"reduced_trust_region", report.policy_name, report.parameters};
  }

  // This contract owns values only. In particular, it intentionally does not
  // retain a CompiledProblemT or ReducedDTOT service: a detached report can
  // therefore be associated with the manifest that produced it without
  // extending the executable service lifetime.
  template <typename PolicySnapshot, typename Report>
  class ReducedExperimentEnvelopeT final
  {
  public:
    using Policy = PolicySnapshot;
    using ReportType = Report;

    ReducedExperimentEnvelopeT(
      compiler::v1::CompilationManifest compilation_manifest,
      PolicySnapshot                       solver_policy,
      Report                                report,
      RunEnvironmentRecord                 environment)
      : compilation_manifest_(std::move(compilation_manifest))
      , solver_policy_(std::move(solver_policy))
      , report_(std::move(report))
      , environment_(std::move(environment))
    {
      contract::require(!compilation_manifest_.semantic_problem_id.empty(),
                        "An experiment envelope needs a compilation manifest identifier");
    }

    const compiler::v1::CompilationManifest &
    compilation_manifest() const
    {
      return compilation_manifest_;
    }

    const PolicySnapshot &
    solver_policy() const
    {
      return solver_policy_;
    }

    const Report &
    report() const
    {
      return report_;
    }

    const RunEnvironmentRecord &
    environment() const
    {
      return environment_;
    }

  private:
    compiler::v1::CompilationManifest compilation_manifest_;
    PolicySnapshot                       solver_policy_;
    Report                               report_;
    RunEnvironmentRecord                 environment_;
  };

  template <typename Backend>
  using ReducedSearchExperimentEnvelopeT =
    ReducedExperimentEnvelopeT<ReducedSearchPolicySnapshot,
                               solvers::ReducedSolverResultT<Backend>>;

  template <typename Backend>
  using ReducedTrustRegionExperimentEnvelopeT =
    ReducedExperimentEnvelopeT<ReducedTrustRegionPolicySnapshot,
                               solvers::ReducedTrustRegionResultT<Backend>>;
} // namespace nmopt::experiment
