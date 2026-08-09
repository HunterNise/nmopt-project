#pragma once

#include "nmopt/contract/layout.hpp"

#include <cstddef>
#include <string>
#include <utility>

namespace nmopt::contract
{
  enum class LinearSolveTermination
  {
    converged,
    failed
  };

  // Backend-neutral evidence returned by a formulation solve service. The
  // policy that produced this report remains owned by the compiler/backend;
  // the reduced formulation only needs deterministic convergence and work
  // information.
  struct LinearSolveReport
  {
    std::string            algorithm;
    std::string            preconditioner;
    std::size_t            maximum_iterations = 0;
    std::size_t            iterations = 0;
    double                 relative_tolerance = 0.0;
    double                 absolute_tolerance = 0.0;
    double                 requested_tolerance = 0.0;
    double                 achieved_residual = 0.0;
    LinearSolveTermination termination = LinearSolveTermination::failed;

    bool
    converged() const
    {
      return termination == LinearSolveTermination::converged;
    }
  };

  template <typename Backend>
  struct FormulationSolveResultT
  {
    PrimalBlockT<Backend> solution;
    LinearSolveReport     report;

    // Retain source compatibility for exact/direct reference callbacks while
    // making compiled backends return an explicit report.
    FormulationSolveResultT(PrimalBlockT<Backend> value)
      : solution(std::move(value))
      , report{"caller-supplied exact solve",
               "not applicable",
               1,
               1,
               0.0,
               0.0,
               0.0,
               0.0,
               LinearSolveTermination::converged}
    {}

    FormulationSolveResultT(PrimalBlockT<Backend> value,
                            LinearSolveReport     solve_report)
      : solution(std::move(value))
      , report(std::move(solve_report))
    {}
  };
} // namespace nmopt::contract
