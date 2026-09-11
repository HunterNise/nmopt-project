#pragma once

#include "../diagnostics/instrumentation.hpp"
#include "../integration/problem_a.hpp"

#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>

namespace external_dealii_step4
{
  class NativeReduced final
  {
  public:
    using Vector         = ProblemA::Vector;
    using SolveEvidence = ProblemA::SolveEvidence;

    struct ValueEvaluation
    {
      Vector         control;
      Vector         state;
      double         objective;
      SolveEvidence state_solve;
    };

    struct DerivativeEvaluation
    {
      Vector         state_derivative;
      Vector         control_derivative;
      Vector         adjoint;
      Vector         reduced_derivative;
      SolveEvidence adjoint_solve;
    };

    NativeReduced(ProblemA &problem, Instrumentation &instrumentation)
      : problem_(problem)
      , instrumentation_(instrumentation)
    {}

    ValueEvaluation
    evaluate_value(const Vector &control)
    {
      require_size(control, problem_.control_dimension(), "control");

      auto solve = problem_.solve_state(control);
      ValueEvaluation result{control,
                             std::move(solve.solution),
                             0.0,
                             solve.evidence};
      result.objective = problem_.objective(result.state, result.control);
      ++instrumentation_.value_evaluations;
      return result;
    }

    DerivativeEvaluation
    augment_derivative(const ValueEvaluation &value)
    {
      require_size(value.control,
                   problem_.control_dimension(),
                   "value control");
      require_size(value.state, problem_.state_dimension(), "value state");

      const auto objective_derivative =
        problem_.objective_derivative(value.state, value.control);
      const auto adjoint =
        problem_.solve_adjoint(objective_derivative.state);
      const auto control_pullback = problem_.control_vjp(adjoint.solution);

      Vector reduced_derivative = objective_derivative.control;
      reduced_derivative.add(-1.0, control_pullback);

      ++instrumentation_.derivative_augmentations;
      return {objective_derivative.state,
              objective_derivative.control,
              std::move(adjoint.solution),
              std::move(reduced_derivative),
              adjoint.evidence};
    }

  private:
    static void
    require_size(const Vector &vector,
                 const std::size_t expected,
                 const char *const name)
    {
      if (static_cast<std::size_t>(vector.size()) != expected)
        throw std::invalid_argument(std::string("native reduced ") + name +
                                    " has the wrong dimension");
    }

    ProblemA &       problem_;
    Instrumentation &instrumentation_;
  };
} // namespace external_dealii_step4
