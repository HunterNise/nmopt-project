#pragma once

#include "nmopt/contract/linalg.hpp"

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace nmopt::dealii_backend
{
  // Finite-dimensional independent state coordinates for a constrained
  // physical FE vector. The homogeneous constraints define the embedding P;
  // the supplied vector is the fixed value lifting ell_0.
  class IndependentStateCoordinates final
  {
  public:
    using Vector = dealii::Vector<double>;

    IndependentStateCoordinates() = default;

    void
    initialize(const std::size_t                         physical_dimension,
               const dealii::AffineConstraints<double> &homogeneous_constraints,
               Vector                                    fixed_lifting)
    {
      contract::require(physical_dimension > 0,
                        "Independent state coordinates need a physical dimension");
      contract::require(fixed_lifting.size() == physical_dimension,
                        "Independent state fixed lifting has an incompatible dimension");

      independent_state_dofs_.clear();
      for (dealii::types::global_dof_index index = 0;
           index < physical_dimension;
           ++index)
        if (!homogeneous_constraints.is_constrained(index))
          independent_state_dofs_.push_back(index);
      contract::require(!independent_state_dofs_.empty(),
                        "Independent state coordinates need an independent DoF");

      dealii::DynamicSparsityPattern reconstruction_dsp(
        physical_dimension, independent_state_dofs_.size());
      for (std::size_t column = 0; column < independent_state_dofs_.size();
           ++column)
        {
          Vector basis(physical_dimension);
          basis[independent_state_dofs_[column]] = 1.0;
          homogeneous_constraints.distribute(basis);
          for (dealii::types::global_dof_index row = 0;
               row < physical_dimension;
               ++row)
            if (basis[row] != 0.0)
              reconstruction_dsp.add(row, column);
        }

      reconstruction_sparsity_.copy_from(reconstruction_dsp);
      reconstruction_.reinit(reconstruction_sparsity_);
      for (std::size_t column = 0; column < independent_state_dofs_.size();
           ++column)
        {
          Vector basis(physical_dimension);
          basis[independent_state_dofs_[column]] = 1.0;
          homogeneous_constraints.distribute(basis);
          for (dealii::types::global_dof_index row = 0;
               row < physical_dimension;
               ++row)
            if (basis[row] != 0.0)
              reconstruction_.set(row, column, basis[row]);
        }

      fixed_lifting_ = std::move(fixed_lifting);
      physical_dimension_ = physical_dimension;
    }

    std::size_t
    physical_dimension() const
    {
      return physical_dimension_;
    }

    std::size_t
    independent_dimension() const
    {
      return independent_state_dofs_.size();
    }

    const Vector &
    fixed_lifting() const
    {
      require_initialized();
      return fixed_lifting_;
    }

    Vector
    reconstruct(const Vector &independent_state) const
    {
      require_initialized();
      require_independent_dimension(independent_state,
                                    "State reconstruction");
      Vector physical = embed(independent_state);
      physical.add(1.0, fixed_lifting_);
      return physical;
    }

    Vector
    embed(const Vector &independent_state) const
    {
      require_initialized();
      require_independent_dimension(independent_state, "State embedding");
      Vector physical(physical_dimension_);
      reconstruction_.vmult(physical, independent_state);
      return physical;
    }

    Vector
    pullback(const Vector &physical_covector) const
    {
      require_initialized();
      contract::require(physical_covector.size() == physical_dimension_,
                        "State pullback received an incompatible physical covector");
      Vector independent(independent_state_dofs_.size());
      reconstruction_.Tvmult(independent, physical_covector);
      return independent;
    }

  private:
    void
    require_initialized() const
    {
      contract::require(physical_dimension_ > 0,
                        "Independent state coordinates are not initialised");
    }

    void
    require_independent_dimension(const Vector &value,
                                  const char   *operation) const
    {
      contract::require(value.size() == independent_state_dofs_.size(),
                        std::string(operation) +
                          " received incompatible independent coordinates");
    }

    std::size_t                                  physical_dimension_ = 0;
    std::vector<dealii::types::global_dof_index> independent_state_dofs_;
    dealii::SparsityPattern                       reconstruction_sparsity_;
    dealii::SparseMatrix<double>                  reconstruction_;
    Vector                                        fixed_lifting_;
  };
} // namespace nmopt::dealii_backend
