#pragma once

#include <array>
#include <cstddef>
#include <vector>

namespace external_dealii_step4
{
  enum class SolveRole
  {
    state,
    adjoint
  };

  enum class ProblemBMatrixAction
  {
    stiffness_apply,
    stiffness_transpose_apply,
    mass_apply,
    coupling_apply,
    coupling_transpose_apply
  };

  enum class ProblemBMatrixPurpose
  {
    unspecified,
    residual,
    residual_jvp,
    residual_vjp,
    objective,
    objective_derivative,
    state_solve,
    control_vjp,
    metric_apply
  };

  enum class MetricSolvePurpose
  {
    gradient_norm
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

  struct ProblemBMatrixActionRecord
  {
    ProblemBMatrixAction  action;
    ProblemBMatrixPurpose purpose;
  };

  struct MetricSolveRecord
  {
    MetricSolvePurpose purpose;
    bool               converged;
    unsigned int       iterations;
    double             initial_residual;
    double             final_residual;
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
    std::size_t metric_apply_calls             = 0;
    std::size_t metric_inverse_apply_calls     = 0;
    std::vector<SolveRecord> solve_records;
    std::vector<SolveFailureRecord> solve_failure_records;
    std::vector<ProblemBMatrixActionRecord> problem_b_matrix_actions;
    std::vector<MetricSolveRecord>          metric_solve_records;

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

    void
    record_problem_b_matrix_action(const ProblemBMatrixAction  action,
                                   const ProblemBMatrixPurpose purpose)
    {
      problem_b_matrix_actions.push_back({action, purpose});
    }

    void
    record_metric_solve(const MetricSolvePurpose purpose,
                        const bool              converged,
                        const unsigned int      iterations,
                        const double            initial_residual,
                        const double            final_residual)
    {
      metric_solve_records.push_back(
        {purpose, converged, iterations, initial_residual, final_residual});
    }
  };

  inline const char *
  problem_b_matrix_action_name(const ProblemBMatrixAction action)
  {
    switch (action)
      {
        case ProblemBMatrixAction::stiffness_apply:
          return "stiffness_apply";
        case ProblemBMatrixAction::stiffness_transpose_apply:
          return "stiffness_transpose_apply";
        case ProblemBMatrixAction::mass_apply:
          return "mass_apply";
        case ProblemBMatrixAction::coupling_apply:
          return "coupling_apply";
        case ProblemBMatrixAction::coupling_transpose_apply:
          return "coupling_transpose_apply";
      }
    return "unknown";
  }

  inline const char *
  problem_b_matrix_purpose_name(const ProblemBMatrixPurpose purpose)
  {
    switch (purpose)
      {
        case ProblemBMatrixPurpose::unspecified:
          return "unspecified";
        case ProblemBMatrixPurpose::residual:
          return "residual";
        case ProblemBMatrixPurpose::residual_jvp:
          return "residual_jvp";
        case ProblemBMatrixPurpose::residual_vjp:
          return "residual_vjp";
        case ProblemBMatrixPurpose::objective:
          return "objective";
        case ProblemBMatrixPurpose::objective_derivative:
          return "objective_derivative";
        case ProblemBMatrixPurpose::state_solve:
          return "state_solve";
        case ProblemBMatrixPurpose::control_vjp:
          return "control_vjp";
        case ProblemBMatrixPurpose::metric_apply:
          return "metric_apply";
      }
    return "unknown";
  }

  inline const char *
  metric_solve_purpose_name(const MetricSolvePurpose purpose)
  {
    switch (purpose)
      {
        case MetricSolvePurpose::gradient_norm:
          return "gradient_norm";
      }
    return "unknown";
  }

  inline std::array<std::size_t, 17>
  runtime_counter_snapshot(const Instrumentation &instrumentation)
  {
    return {{instrumentation.assembly_calls,
             instrumentation.state_solve_calls,
             instrumentation.adjoint_solve_calls,
             instrumentation.solve_failures,
             instrumentation.objective_calls,
             instrumentation.objective_derivative_calls,
             instrumentation.residual_calls,
             instrumentation.residual_jvp_calls,
             instrumentation.residual_vjp_calls,
             instrumentation.control_vjp_calls,
             instrumentation.explicit_matrix_vmult_calls,
             instrumentation.explicit_matrix_tvmult_calls,
             instrumentation.value_evaluations,
             instrumentation.derivative_augmentations,
             instrumentation.output_calls,
             instrumentation.metric_apply_calls,
             instrumentation.metric_inverse_apply_calls}};
  }
} // namespace external_dealii_step4
