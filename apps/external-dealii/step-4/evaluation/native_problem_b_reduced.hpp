#pragma once

#include "../diagnostics/instrumentation.hpp"
#include "../integration/problem_b.hpp"

#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>

namespace external_dealii_step4
{
  template <int dim, typename Application>
  class NativeProblemBReduced final
  {
  public:
    using Problem        = ProblemB<dim, Application>;
    using Vector         = typename Problem::Vector;
    using SolveEvidence  = typename Problem::SolveEvidence;
    using SolveResult    = typename Problem::SolveResult;

    struct ValueEvaluation
    {
      Vector         control;
      Vector         state;
      Vector         full_state;
      double         objective;
      SolveEvidence state_solve;
    };

    struct DerivativeEvaluation
    {
      Vector         state_derivative;
      Vector         control_derivative;
      Vector         adjoint;
      Vector         full_adjoint;
      Vector         reduced_derivative;
      SolveEvidence adjoint_solve;
    };

    explicit NativeProblemBReduced(Problem &problem)
      : problem_(problem)
      , instrumentation_(nullptr)
    {}

    NativeProblemBReduced(Problem &problem, Instrumentation &instrumentation)
      : problem_(problem)
      , instrumentation_(&instrumentation)
    {}

    NativeProblemBReduced(const NativeProblemBReduced &) = delete;
    NativeProblemBReduced &operator=(const NativeProblemBReduced &) = delete;

    ValueEvaluation
    evaluate_value(const Vector &control)
    {
      require_size(control, problem_.control_dimension(), "control");
      count(&Instrumentation::state_solve_calls);
      SolveResult solve;
      try
        {
          solve = problem_.solve_state(control);
          record_solve(SolveRole::state, solve.evidence);
        }
      catch (...)
        {
          if (instrumentation_ != nullptr)
            instrumentation_->record_solve_failure();
          throw;
        }

      ValueEvaluation result{control,
                             std::move(solve.solution),
                             std::move(solve.full_solution),
                             0.0,
                             std::move(solve.evidence)};
      count(&Instrumentation::objective_calls);
      result.objective = problem_.objective(result.state, result.control);
      count(&Instrumentation::value_evaluations);
      return result;
    }

    DerivativeEvaluation
    augment_derivative(const ValueEvaluation &value)
    {
      require_size(value.control,
                   problem_.control_dimension(),
                   "value control");
      require_size(value.state, problem_.state_dimension(), "value state");
      require_size(value.full_state,
                   problem_.control_dimension(),
                   "value full state");

      count(&Instrumentation::objective_derivative_calls);
      const auto objective_derivative =
        problem_.objective_derivative(value.state, value.control);

      count(&Instrumentation::adjoint_solve_calls);
      SolveResult adjoint;
      try
        {
          adjoint = problem_.solve_adjoint(objective_derivative.state);
          record_solve(SolveRole::adjoint, adjoint.evidence);
        }
      catch (...)
        {
          if (instrumentation_ != nullptr)
            instrumentation_->record_solve_failure();
          throw;
        }

      count(&Instrumentation::control_vjp_calls);
      const auto control_pullback = problem_.control_vjp(adjoint.solution);

      Vector reduced_derivative = objective_derivative.control;
      reduced_derivative.add(-1.0, control_pullback);
      count(&Instrumentation::derivative_augmentations);
      return {objective_derivative.state,
              objective_derivative.control,
              std::move(adjoint.solution),
              std::move(adjoint.full_solution),
              std::move(reduced_derivative),
              std::move(adjoint.evidence)};
    }

  private:
    using Counter = std::size_t Instrumentation::*;

    void
    count(const Counter counter) const
    {
      if (instrumentation_ != nullptr)
        ++(instrumentation_->*counter);
    }

    void
    record_solve(const SolveRole role, const SolveEvidence &evidence) const
    {
      if (instrumentation_ == nullptr)
        return;
      if (evidence.converged)
        instrumentation_->record_solve_success(role,
                                               evidence.iterations,
                                               evidence.initial_residual,
                                               evidence.final_residual);
      else
        instrumentation_->record_solve_failure(role,
                                               evidence.iterations,
                                               evidence.final_residual);
    }

    static void
    require_size(const Vector &vector,
                 const std::size_t expected,
                 const char *const name)
    {
      if (static_cast<std::size_t>(vector.size()) != expected)
        throw std::invalid_argument(std::string("native Problem B ") + name +
                                    " has the wrong dimension");
    }

    Problem &       problem_;
    Instrumentation *instrumentation_;
  };
} // namespace external_dealii_step4
