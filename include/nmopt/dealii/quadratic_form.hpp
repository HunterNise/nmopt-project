#pragma once

#include "nmopt/contract/linalg.hpp"

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <string>

namespace nmopt::dealii_backend
{
  // Non-owning finite-dimensional quadratic objective algebra. The supplied
  // matrix and affine load remain owned by the numerical target; this value
  // only centralizes the common value, gradient, and curvature actions.
  class QuadraticForm final
  {
  public:
    using Matrix = dealii::SparseMatrix<double>;
    using Vector = dealii::Vector<double>;

    explicit QuadraticForm(const Matrix &matrix)
      : matrix_(&matrix)
    {
      require_matrix();
    }

    QuadraticForm(const Matrix & matrix,
                  const Vector & linear,
                  const double   constant)
      : matrix_(&matrix)
      , linear_(&linear)
      , constant_(constant)
    {
      require_matrix();
      contract::require(linear_->size() == matrix_->n(),
                        "Quadratic-form affine load has an incompatible dimension");
    }

    double
    value(const Vector &coordinates) const
    {
      const Vector action = hessian_action(coordinates);
      const double linear_value = linear_ == nullptr ? 0.0 :
                                                       *linear_ * coordinates;
      return 0.5 * (coordinates * action) - linear_value + 0.5 * constant_;
    }

    Vector
    gradient(const Vector &coordinates) const
    {
      Vector result = hessian_action(coordinates);
      if (linear_ != nullptr)
        result.add(-1.0, *linear_);
      return result;
    }

    Vector
    hessian_action(const Vector &direction) const
    {
      require_coordinates(direction, "Quadratic-form action");
      Vector result(matrix_->m());
      matrix_->vmult(result, direction);
      return result;
    }

  private:
    void
    require_matrix() const
    {
      contract::require(matrix_->m() == matrix_->n(),
                        "Quadratic-form matrix must be square");
    }

    void
    require_coordinates(const Vector &coordinates, const char *operation) const
    {
      contract::require(coordinates.size() == matrix_->n(),
                        std::string(operation) +
                          " received an incompatible dimension");
    }

    const Matrix *matrix_ = nullptr;
    const Vector *linear_ = nullptr;
    double        constant_ = 0.0;
  };
} // namespace nmopt::dealii_backend
