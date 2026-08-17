#pragma once

#include "nmopt/application/metadata.hpp"

namespace nmopt::application
{
  // A scenario binds one recipe to typed problem, compilation, solver, and
  // experiment choices. It does not compile or execute the problem.
  template <typename ProblemParameters,
            typename CompileOptions,
            typename SolverOptions,
            typename ExperimentOptions>
  struct ScenarioT
  {
    using problem_parameters_type = ProblemParameters;
    using compile_options_type = CompileOptions;
    using solver_options_type = SolverOptions;
    using experiment_options_type = ExperimentOptions;

    ScenarioMetadata   metadata;
    ProblemParameters  problem;
    CompileOptions     compile;
    SolverOptions      solver;
    ExperimentOptions  experiment;

    void
    validate() const
    {
      detail::validate_scenario_metadata(metadata);
    }
  };
} // namespace nmopt::application
