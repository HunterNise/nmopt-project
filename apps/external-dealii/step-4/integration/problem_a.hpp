#pragma once

#include "adapted_step4.hpp"
#include "../diagnostics/instrumentation.hpp"

#include <deal.II/lac/solver_control.h>

#include <cstddef>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <utility>

namespace external_dealii_step4
{
  class ProblemA final
  {
  public:
    using Application  = Step4<2>;
    using Matrix       = dealii::SparseMatrix<double>;
    using Vector       = dealii::Vector<double>;
    using SolveEvidence = Application::SolveEvidence;

    struct ResidualDerivative
    {
      Vector state;
      Vector control;
    };

    struct ObjectiveDerivative
    {
      Vector state;
      Vector control;
    };

    struct SolveResult
    {
      Vector         solution;
      SolveEvidence evidence;
    };

    ProblemA()
      : ProblemA(nullptr)
    {}

    explicit ProblemA(Instrumentation &instrumentation)
      : ProblemA(&instrumentation)
    {}

    explicit ProblemA(Instrumentation *const instrumentation)
      : instrumentation_(instrumentation)
    {
      increment(&Instrumentation::assembly_calls);
      application_.prepare_for_external_use();
    }

    std::size_t
    state_dimension() const
    {
      return application_.system_rhs_view().size();
    }

    std::size_t
    control_dimension() const
    {
      return state_dimension();
    }

    const Matrix &
    system_matrix() const
    {
      return application_.system_matrix_view();
    }

    const Vector &
    system_rhs() const
    {
      return application_.system_rhs_view();
    }

    Vector
    residual(const Vector &state, const Vector &control) const
    {
      require_size(state, state_dimension(), "state");
      require_size(control, control_dimension(), "control");

      increment(&Instrumentation::residual_calls);
      Vector value(state_dimension());
      system_matrix().vmult(value, state);
      increment(&Instrumentation::explicit_matrix_vmult_calls);
      value.add(-1.0, system_rhs());
      value.add(-1.0, control);
      return value;
    }

    Vector
    residual_jvp(const Vector &state_tangent,
                 const Vector &control_tangent) const
    {
      require_size(state_tangent, state_dimension(), "state tangent");
      require_size(control_tangent, control_dimension(), "control tangent");

      increment(&Instrumentation::residual_jvp_calls);
      Vector value(state_dimension());
      system_matrix().vmult(value, state_tangent);
      increment(&Instrumentation::explicit_matrix_vmult_calls);
      value.add(-1.0, control_tangent);
      return value;
    }

    ResidualDerivative
    residual_vjp(const Vector &test_seed) const
    {
      require_size(test_seed, state_dimension(), "test seed");

      increment(&Instrumentation::residual_vjp_calls);
      ResidualDerivative result{Vector(state_dimension()),
                                Vector(control_dimension())};
      system_matrix().Tvmult(result.state, test_seed);
      increment(&Instrumentation::explicit_matrix_tvmult_calls);
      result.control = test_seed;
      result.control *= -1.0;
      return result;
    }

    Vector
    control_vjp(const Vector &test_seed) const
    {
      require_size(test_seed, state_dimension(), "control test seed");

      increment(&Instrumentation::control_vjp_calls);
      Vector result = test_seed;
      result *= -1.0;
      return result;
    }

    double
    objective(const Vector &state, const Vector &control) const
    {
      require_size(state, state_dimension(), "state");
      require_size(control, control_dimension(), "control");

      increment(&Instrumentation::objective_calls);
      return 0.5 * (state * state) + 0.5 * (control * control);
    }

    ObjectiveDerivative
    objective_derivative(const Vector &state, const Vector &control) const
    {
      require_size(state, state_dimension(), "state");
      require_size(control, control_dimension(), "control");

      increment(&Instrumentation::objective_derivative_calls);
      return {state, control};
    }

    SolveResult
    solve_state(const Vector &control) const
    {
      require_size(control, control_dimension(), "control");
      Vector rhs = system_rhs();
      rhs.add(1.0, control);
      return solve(rhs, SolveRole::state);
    }

    SolveResult
    solve_adjoint(const Vector &state_objective_derivative) const
    {
      require_size(state_objective_derivative,
                   state_dimension(),
                   "state objective derivative");
      return solve(state_objective_derivative, SolveRole::adjoint);
    }

    void
    output_results(const Vector &                 state,
                   const std::filesystem::path &filename) const
    {
      require_size(state, state_dimension(), "output state");
      increment(&Instrumentation::output_calls);
      application_.output_results(state, filename);
    }

  private:
    void
    increment(std::size_t Instrumentation::*const counter) const
    {
      if (instrumentation_ != nullptr)
        ++(instrumentation_->*counter);
    }

    static void
    require_size(const Vector &vector,
                 const std::size_t expected,
                 const char *const name)
    {
      if (static_cast<std::size_t>(vector.size()) != expected)
        throw std::invalid_argument(std::string("Problem A ") + name +
                                    " has the wrong dimension");
    }

    SolveResult
    solve(const Vector &rhs, const SolveRole role) const
    {
      if (instrumentation_ != nullptr)
        instrumentation_->record_solve_start(role);
      Vector solution(rhs.size());
      solution = 0.0;
      try
        {
          const auto evidence = application_.solve(rhs, solution);
          if (evidence.converged)
            {
              if (instrumentation_ != nullptr)
                instrumentation_->record_solve_success(
                  role,
                  evidence.iterations,
                  evidence.initial_residual,
                  evidence.final_residual);
            }
          else
            {
              if (instrumentation_ != nullptr)
                instrumentation_->record_solve_failure(
                  role, evidence.iterations, evidence.final_residual);
            }
          return {std::move(solution), evidence};
        }
      catch (const dealii::SolverControl::NoConvergence &exception)
        {
          if (instrumentation_ != nullptr)
            instrumentation_->record_solve_failure(
              role, exception.last_step, exception.last_residual);
          throw;
        }
      catch (...)
        {
          if (instrumentation_ != nullptr)
            instrumentation_->record_solve_failure();
          throw;
        }
    }

    Application       application_;
    Instrumentation *instrumentation_;
  };
} // namespace external_dealii_step4
