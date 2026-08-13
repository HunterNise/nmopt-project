#pragma once

#include "nmopt/contract/metric_constraint.hpp"
#include "nmopt/dealii/metric_solve_parameters.hpp"
#include "nmopt/dealii/serial_backend.hpp"
#include "nmopt/dealii/serial_spd_solver.hpp"

#include <deal.II/lac/sparse_matrix.h>

#include <memory>
#include <string>
#include <utility>

namespace nmopt::dealii_backend
{
  enum class Hminus1OperatorRealisation
  {
    mass_laplacian_inverse_mass
  };

  enum class Hminus1InverseRealisation
  {
    mass_inverse_laplacian_mass_inverse
  };

  // Discrete negative-norm metric on an independent H1_0 coefficient space:
  //
  //   G = M K^{-1} M,        G^{-1} = M^{-1} K M^{-1}.
  //
  // K is the homogeneous-Dirichlet Laplacian and M is the control mass
  // matrix. Both must be symmetric positive definite on the supplied layout.
  class Hminus1Metric final : public contract::MetricT<SerialBackend>
  {
  public:
    using Matrix = dealii::SparseMatrix<double>;
    using Primal = contract::PrimalBlockT<SerialBackend>;
    using Covector = contract::CovectorBlockT<SerialBackend>;
    using Vector = dealii::Vector<double>;

    Hminus1Metric(std::string                   id,
                  contract::LayoutPtr           layout,
                  std::shared_ptr<const Matrix> mass_matrix,
                  std::shared_ptr<const Matrix> laplace_matrix,
                  MetricSolveParameters         solve_parameters = {},
                  const Hminus1OperatorRealisation operator_realisation =
                    Hminus1OperatorRealisation::mass_laplacian_inverse_mass,
                  const Hminus1InverseRealisation inverse_realisation =
                    Hminus1InverseRealisation::mass_inverse_laplacian_mass_inverse)
      : id_(std::move(id))
      , layout_(std::move(layout))
      , mass_matrix_(std::move(mass_matrix))
      , laplace_matrix_(std::move(laplace_matrix))
      , solve_parameters_(solve_parameters)
      , operator_realisation_(operator_realisation)
      , inverse_realisation_(inverse_realisation)
    {
      contract::require(!id_.empty(),
                        "H-1 metric identifier must not be empty");
      contract::require(static_cast<bool>(layout_),
                        "H-1 metric needs a layout");
      contract::require(layout_->n_blocks() == 1,
                        "H-1 metric supports exactly one coefficient block");
      require_matrix(mass_matrix_, "mass");
      require_matrix(laplace_matrix_, "Dirichlet Laplace");
      contract::require(solve_parameters_.maximum_iterations > 0,
                        "H-1 metric maximum iterations must be positive");
      contract::require(solve_parameters_.relative_tolerance > 0.0,
                        "H-1 metric relative tolerance must be positive");
      contract::require(solve_parameters_.absolute_tolerance > 0.0,
                        "H-1 metric absolute tolerance must be positive");
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

    const MetricSolveParameters &
    solve_parameters() const
    {
      return solve_parameters_;
    }

    Covector
    apply(const Primal &primal) const override
    {
      require_primal(primal);
      contract::require(
        operator_realisation_ ==
          Hminus1OperatorRealisation::mass_laplacian_inverse_mass,
        "H-1 metric received an unsupported operator realization");
      Vector mass_action(layout_->dimension(0));
      mass_matrix_->vmult(mass_action, primal.block(0));
      const Vector potential = solve(*laplace_matrix_, mass_action, "Laplacian");
      Vector value(layout_->dimension(0));
      mass_matrix_->vmult(value, potential);
      return Covector(layout_, {std::move(value)});
    }

    Primal
    inverse_apply(const Covector &covector) const override
    {
      require_covector(covector);
      contract::require(
        inverse_realisation_ ==
          Hminus1InverseRealisation::mass_inverse_laplacian_mass_inverse,
        "H-1 metric received an unsupported inverse realization");
      const Vector mass_representative =
        solve(*mass_matrix_, covector.block(0), "first mass");
      Vector laplace_action(layout_->dimension(0));
      laplace_matrix_->vmult(laplace_action, mass_representative);
      Vector value = solve(*mass_matrix_, laplace_action, "second mass");
      return Primal(layout_, {std::move(value)});
    }

  private:
    void
    require_matrix(const std::shared_ptr<const Matrix> &matrix,
                   const char *                        name) const
    {
      contract::require(static_cast<bool>(matrix),
                        std::string("H-1 metric needs a ") + name + " matrix");
      contract::require(matrix->m() == matrix->n(),
                        std::string("H-1 metric ") + name +
                          " matrix must be square");
      contract::require(matrix->m() == layout_->dimension(0),
                        std::string("H-1 metric ") + name +
                          " matrix dimension does not match its layout");
    }

    Vector
    solve(const Matrix &matrix,
          const Vector &right_hand_side,
          const char *  operation) const
    {
      Vector solution(right_hand_side.size());
      const SPDLinearSolvePolicy policy{solve_parameters_.maximum_iterations,
                                        solve_parameters_.relative_tolerance,
                                        solve_parameters_.absolute_tolerance};
      const auto report =
        solve_serial_spd(matrix, solution, right_hand_side, policy);
      contract::require(report.converged(),
                        std::string("H-1 metric ") + operation +
                          " solve did not converge");
      return solution;
    }

    void
    require_primal(const Primal &primal) const
    {
      contract::require(primal.layout()->compatible_with(*layout_),
                        "H-1 metric primal has an incompatible layout");
    }

    void
    require_covector(const Covector &covector) const
    {
      contract::require(covector.layout()->compatible_with(*layout_),
                        "H-1 metric covector has an incompatible layout");
    }

    std::string                   id_;
    contract::LayoutPtr           layout_;
    std::shared_ptr<const Matrix> mass_matrix_;
    std::shared_ptr<const Matrix> laplace_matrix_;
    MetricSolveParameters         solve_parameters_;
    Hminus1OperatorRealisation    operator_realisation_;
    Hminus1InverseRealisation     inverse_realisation_;
  };
} // namespace nmopt::dealii_backend
