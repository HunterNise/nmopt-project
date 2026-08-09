#pragma once

#include "nmopt/contract/linear_solve.hpp"

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace nmopt::dealii_backend
{
  struct SPDLinearSolvePolicy
  {
    // Zero selects the established dimension-dependent rule
    // max(100, 10 * operator dimension).
    unsigned int maximum_iterations = 0;
    double       relative_tolerance = 1e-12;
    double       absolute_tolerance = 1e-14;
  };

  inline bool
  valid(const SPDLinearSolvePolicy &policy)
  {
    return std::isfinite(policy.relative_tolerance) &&
           policy.relative_tolerance > 0.0 &&
           std::isfinite(policy.absolute_tolerance) &&
           policy.absolute_tolerance > 0.0;
  }

  template <typename Matrix, typename Vector>
  contract::LinearSolveReport
  solve_serial_spd(const Matrix &               matrix,
                   Vector &                     solution,
                   const Vector &               right_hand_side,
                   const SPDLinearSolvePolicy & policy)
  {
    contract::require(valid(policy),
                      "Serial SPD solve policy needs positive finite tolerances");
    contract::require(matrix.m() == matrix.n(),
                      "Serial SPD solve needs a square operator");
    contract::require(matrix.m() == right_hand_side.size() &&
                        matrix.n() == solution.size(),
                      "Serial SPD solve operator and vectors have incompatible dimensions");

    const auto native_dimension = matrix.m();
    contract::require(
      native_dimension <=
        std::numeric_limits<unsigned int>::max() / 10U,
      "Serial SPD solve dimension exceeds the iteration-policy range");
    const unsigned int maximum_iterations =
      policy.maximum_iterations == 0
        ? std::max(100U, 10U * static_cast<unsigned int>(native_dimension))
        : policy.maximum_iterations;
    const double requested_tolerance =
      std::max(policy.absolute_tolerance,
               policy.relative_tolerance * right_hand_side.l2_norm());

    dealii::SolverControl control(maximum_iterations, requested_tolerance);
    dealii::SolverCG<Vector> solver(control);
    bool converged = true;
    try
      {
        solver.solve(matrix,
                     solution,
                     right_hand_side,
                     dealii::PreconditionIdentity());
      }
    catch (const dealii::SolverControl::NoConvergence &)
      {
        converged = false;
      }

    return {"serial_cg",
            "identity",
            maximum_iterations,
            control.last_step(),
            policy.relative_tolerance,
            policy.absolute_tolerance,
            requested_tolerance,
            control.last_value(),
            converged ? contract::LinearSolveTermination::converged
                      : contract::LinearSolveTermination::failed};
  }

  inline contract::LinearSolveReport
  direct_solve_report(const std::string &algorithm)
  {
    return {algorithm,
            "not applicable",
            1,
            1,
            0.0,
            0.0,
            0.0,
            0.0,
            contract::LinearSolveTermination::converged};
  }
} // namespace nmopt::dealii_backend
