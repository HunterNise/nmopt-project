#pragma once

#include "nmopt/contract/metric_constraint.hpp"
#include "nmopt/dealii/serial_backend.hpp"

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparse_matrix.h>

#include <algorithm>
#include <memory>
#include <string>
#include <utility>

namespace nmopt::dealii_backend
{
  struct MassMetricSolveParameters
  {
    unsigned int maximum_iterations = 1000;
    double       relative_tolerance = 1e-12;
    double       absolute_tolerance = 1e-14;
  };

  // A serial deal.II realization of a one-block Hilbert metric. The supplied
  // matrix must be symmetric positive definite and represents the
  // primal-to-dual Riesz action; its inverse is evaluated by CG with the
  // explicitly declared solve parameters above.
  class MassMetric final : public contract::MetricT<SerialBackend>
  {
  public:
    using Matrix = dealii::SparseMatrix<double>;
    using Primal = contract::PrimalBlockT<SerialBackend>;
    using Covector = contract::CovectorBlockT<SerialBackend>;

    MassMetric(std::string                        id,
               contract::LayoutPtr                layout,
               std::shared_ptr<const Matrix>      matrix,
               MassMetricSolveParameters solve_parameters = {})
      : id_(std::move(id))
      , layout_(std::move(layout))
      , matrix_(std::move(matrix))
      , solve_parameters_(solve_parameters)
    {
      contract::require(!id_.empty(), "Metric identifier must not be empty");
      contract::require(static_cast<bool>(layout_), "Mass metric needs a layout");
      contract::require(layout_->n_blocks() == 1,
                        "Mass metric supports exactly one coefficient block");
      contract::require(static_cast<bool>(matrix_), "Mass metric needs a matrix");
      contract::require(matrix_->m() == matrix_->n(),
                        "Mass metric matrix must be square");
      contract::require(matrix_->m() == layout_->dimension(0),
                        "Mass metric matrix dimension does not match its layout");
      contract::require(solve_parameters_.maximum_iterations > 0,
                        "Mass metric maximum iterations must be positive");
      contract::require(solve_parameters_.relative_tolerance > 0.0,
                        "Mass metric relative tolerance must be positive");
      contract::require(solve_parameters_.absolute_tolerance > 0.0,
                        "Mass metric absolute tolerance must be positive");
    }

    const std::string &
    id() const override
    {
      return id_;
    }

    const contract::LayoutPtr &
    layout() const override
    {
      return layout_;
    }

    const MassMetricSolveParameters &
    solve_parameters() const
    {
      return solve_parameters_;
    }

    Covector
    apply(const Primal &primal) const override
    {
      require_primal(primal, "Mass metric primal has an incompatible layout");

      Vector value(matrix_->m());
      matrix_->vmult(value, primal.block(0));
      return Covector(layout_, {std::move(value)});
    }

    Primal
    inverse_apply(const Covector &covector) const override
    {
      require_covector(covector,
                       "Mass metric covector has an incompatible layout");

      Vector value(matrix_->m());
      const double tolerance =
        std::max(solve_parameters_.absolute_tolerance,
                 solve_parameters_.relative_tolerance *
                   covector.block(0).l2_norm());
      dealii::SolverControl control(solve_parameters_.maximum_iterations,
                                    tolerance);
      dealii::SolverCG<Vector> solver(control);
      solver.solve(*matrix_, value, covector.block(0),
                   dealii::PreconditionIdentity());
      return Primal(layout_, {std::move(value)});
    }

  private:
    using Vector = dealii::Vector<double>;

    void
    require_primal(const Primal &primal, const char *message) const
    {
      contract::require(primal.layout()->compatible_with(*layout_), message);
    }

    void
    require_covector(const Covector &covector, const char *message) const
    {
      contract::require(covector.layout()->compatible_with(*layout_), message);
    }

    std::string                        id_;
    contract::LayoutPtr                layout_;
    std::shared_ptr<const Matrix>      matrix_;
    MassMetricSolveParameters          solve_parameters_;
  };
} // namespace nmopt::dealii_backend
