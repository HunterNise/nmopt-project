#pragma once

#include "nmopt/contract/quadratic_kkt.hpp"
#include "nmopt/contract/linear_solve.hpp"

#include <cmath>
#include <cstddef>
#include <string>

namespace nmopt::contract
{
  enum class QuadraticKKTSolverMethod
  {
    minres,
    gmres
  };

  struct QuadraticKKTSolverPolicy
  {
    QuadraticKKTSolverMethod method = QuadraticKKTSolverMethod::minres;
    std::size_t maximum_iterations = 0;
    double      relative_tolerance = 1e-10;
    double      absolute_tolerance = 1e-12;
    unsigned int gmres_maximum_basis = 30;
  };

  inline bool
  valid(const QuadraticKKTSolverPolicy &policy)
  {
    return std::isfinite(policy.relative_tolerance) &&
           policy.relative_tolerance > 0.0 &&
           std::isfinite(policy.absolute_tolerance) &&
           policy.absolute_tolerance > 0.0 &&
           policy.gmres_maximum_basis >= 3;
  }

  template <typename Backend>
  void
  validate(const EqualityConstrainedQuadraticKKTProductT<Backend> &product,
           const QuadraticKKTSolverPolicy &                         policy)
  {
    require(std::isfinite(policy.relative_tolerance) &&
              policy.relative_tolerance > 0.0 &&
              std::isfinite(policy.absolute_tolerance) &&
              policy.absolute_tolerance > 0.0,
            "Quadratic KKT solver policy needs positive finite tolerances");
    require(policy.gmres_maximum_basis >= 3,
            "Quadratic KKT solver policy needs a GMRES basis of at least three vectors");
    if (policy.method == QuadraticKKTSolverMethod::minres)
      require(product.supports_minres(),
              "MINRES requires a symmetric-indefinite KKT product");
  }

  struct QuadraticKKTSolveReport
  {
    LinearSolveReport linear_solve;
    double            stationarity_residual = 0.0;
    double            equality_residual = 0.0;
    bool              residuals_converged = false;

    bool
    converged() const
    {
      return linear_solve.converged() && residuals_converged;
    }
  };

  template <typename Backend>
  struct QuadraticKKTSolveResultT
  {
    QuadraticKKTPointT<Backend> solution;
    QuadraticKKTSolveReport     report;
  };

  using QuadraticKKTSolveResult =
    QuadraticKKTSolveResultT<DenseBackend>;
} // namespace nmopt::contract
