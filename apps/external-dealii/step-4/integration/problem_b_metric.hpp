#pragma once

#include "problem_b_mass.hpp"

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/vector.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>

namespace external_dealii_step4
{
  template <int dim>
  class ProblemBMetric final
  {
  public:
    using Mass   = ProblemBMass<dim>;
    using Vector = dealii::Vector<double>;

    struct SolveEvidence
    {
      bool         converged;
      unsigned int iterations;
      double       initial_residual;
      double       final_residual;
    };

    struct InverseResult
    {
      Vector         solution;
      SolveEvidence evidence;
    };

    explicit ProblemBMetric(const Mass &mass)
      : mass_(mass)
    {}

    Vector
    apply(const Vector &vector) const
    {
      require_size(vector);
      require_finite(vector);
      return mass_.mass_apply(vector);
    }

    InverseResult
    inverse_apply(const Vector &rhs) const
    {
      require_size(rhs);
      require_finite(rhs);

      const double tolerance =
        std::max(1.0e-14, 1.0e-12 * rhs.l2_norm());
      dealii::SolverControl solver_control(1000, tolerance);
      dealii::SolverCG<Vector> solver(solver_control);
      Vector solution(rhs.size());
      solution = 0.0;
      solver.solve(mass_.mass_matrix(),
                   solution,
                   rhs,
                   dealii::PreconditionIdentity());

      return {std::move(solution),
              {solver_control.last_check() == dealii::SolverControl::success,
               solver_control.last_step(),
               solver_control.initial_value(),
               solver_control.last_value()}};
    }

  private:
    void
    require_size(const Vector &vector) const
    {
      if (static_cast<std::size_t>(vector.size()) != mass_.full_dimension())
        throw std::invalid_argument(
          "Problem B metric vector has the wrong dimension");
    }

    static void
    require_finite(const Vector &vector)
    {
      for (unsigned int index = 0; index < vector.size(); ++index)
        if (!std::isfinite(vector[index]))
          throw std::invalid_argument(
            "Problem B metric vector contains a non-finite value");
    }

    const Mass &mass_;
  };
} // namespace external_dealii_step4
