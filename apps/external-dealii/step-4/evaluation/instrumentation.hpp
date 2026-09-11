#pragma once

#include <cstddef>
#include <vector>

namespace external_dealii_step4
{
  enum class SolveRole
  {
    state,
    adjoint
  };

  struct SolveRecord
  {
    SolveRole    role;
    unsigned int iterations;
    double       initial_residual;
    double       final_residual;
  };

  struct SolveFailureRecord
  {
    SolveRole    role;
    unsigned int iterations;
    double       final_residual;
  };

  struct Instrumentation
  {
    std::size_t assembly_calls                 = 0;
    std::size_t state_solve_calls              = 0;
    std::size_t adjoint_solve_calls            = 0;
    std::size_t solve_failures                 = 0;
    std::size_t objective_calls                = 0;
    std::size_t objective_derivative_calls     = 0;
    std::size_t residual_calls                 = 0;
    std::size_t residual_jvp_calls             = 0;
    std::size_t residual_vjp_calls             = 0;
    std::size_t control_vjp_calls              = 0;
    std::size_t explicit_matrix_vmult_calls    = 0;
    std::size_t explicit_matrix_tvmult_calls   = 0;
    std::size_t value_evaluations              = 0;
    std::size_t derivative_augmentations       = 0;
    std::size_t output_calls                   = 0;
    std::vector<SolveRecord> solve_records;
    std::vector<SolveFailureRecord> solve_failure_records;

    void
    record_solve_start(const SolveRole role)
    {
      if (role == SolveRole::state)
        ++state_solve_calls;
      else
        ++adjoint_solve_calls;
    }

    void
    record_solve_success(const SolveRole                 role,
                         const unsigned int              iterations,
                         const double                    initial_residual,
                         const double                    final_residual)
    {
      solve_records.push_back(
        {role, iterations, initial_residual, final_residual});
    }

    void
    record_solve_failure()
    {
      ++solve_failures;
    }

    void
    record_solve_failure(const SolveRole    role,
                         const unsigned int iterations,
                         const double       final_residual)
    {
      ++solve_failures;
      solve_failure_records.push_back({role, iterations, final_residual});
    }
  };
} // namespace external_dealii_step4
