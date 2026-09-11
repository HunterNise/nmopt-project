#pragma once

#include "problem_b_coordinates.hpp"

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector.h>

#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace external_dealii_step4
{
  template <int dim>
  class ProblemBMass final
  {
  public:
    using Coordinates = ProblemBCoordinates<dim>;
    using Index       = typename Coordinates::Index;
    using Matrix      = dealii::SparseMatrix<double>;
    using Vector      = dealii::Vector<double>;

    ProblemBMass(const dealii::DoFHandler<dim> &dof_handler,
                 const Coordinates &            coordinates)
      : full_dimension_(dof_handler.n_dofs())
      , free_dimension_(coordinates.free_dimension())
    {
      if (full_dimension_ == 0 ||
          full_dimension_ != coordinates.full_dimension())
        throw std::invalid_argument(
          "Problem B mass received incompatible DoF coordinates");

      if (full_dimension_ >
          static_cast<std::size_t>(std::numeric_limits<unsigned int>::max()))
        throw std::invalid_argument(
          "Problem B mass exceeds the vector dimension limit");

      assemble_mass(dof_handler);
      assemble_coupling(coordinates);
    }

    std::size_t
    full_dimension() const
    {
      return full_dimension_;
    }

    std::size_t
    free_dimension() const
    {
      return free_dimension_;
    }

    const Matrix &
    mass_matrix() const
    {
      return mass_matrix_;
    }

    const Matrix &
    coupling_matrix() const
    {
      return coupling_matrix_;
    }

    Vector
    mass_apply(const Vector &full_vector) const
    {
      require_size(full_vector, full_dimension_, "full mass input");

      Vector result(full_dimension_);
      mass_matrix_.vmult(result, full_vector);
      return result;
    }

    Vector
    coupling_apply(const Vector &full_vector) const
    {
      require_size(full_vector, full_dimension_, "full coupling input");

      Vector result(free_dimension_);
      coupling_matrix_.vmult(result, full_vector);
      return result;
    }

    Vector
    coupling_transpose_apply(const Vector &free_vector) const
    {
      require_size(free_vector, free_dimension_, "free coupling input");

      Vector result(full_dimension_);
      coupling_matrix_.Tvmult(result, free_vector);
      return result;
    }

  private:
    void
    assemble_mass(const dealii::DoFHandler<dim> &dof_handler)
    {
      dealii::DynamicSparsityPattern sparsity(full_dimension_,
                                               full_dimension_);
      dealii::DoFTools::make_sparsity_pattern(dof_handler, sparsity);

      mass_sparsity_.copy_from(sparsity);
      mass_matrix_.reinit(mass_sparsity_);

      const auto &finite_element = dof_handler.get_fe();
      const dealii::QGauss<dim> quadrature_formula(2);
      dealii::FEValues<dim> fe_values(finite_element,
                                      quadrature_formula,
                                      dealii::update_values |
                                        dealii::update_JxW_values);

      const unsigned int dofs_per_cell = finite_element.n_dofs_per_cell();
      dealii::FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
      std::vector<Index> local_dof_indices(dofs_per_cell);

      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          fe_values.reinit(cell);
          cell_matrix = 0.0;

          for (const unsigned int q_index :
               fe_values.quadrature_point_indices())
            for (const unsigned int i : fe_values.dof_indices())
              for (const unsigned int j : fe_values.dof_indices())
                cell_matrix(i, j) +=
                  fe_values.shape_value(i, q_index) *
                  fe_values.shape_value(j, q_index) * fe_values.JxW(q_index);

          cell->get_dof_indices(local_dof_indices);
          for (const unsigned int i : fe_values.dof_indices())
            for (const unsigned int j : fe_values.dof_indices())
              mass_matrix_.add(local_dof_indices[i],
                               local_dof_indices[j],
                               cell_matrix(i, j));
        }
    }

    void
    assemble_coupling(const Coordinates &coordinates)
    {
      dealii::DynamicSparsityPattern sparsity(free_dimension_,
                                               full_dimension_);
      const auto &free_indices = coordinates.free_indices();

      for (std::size_t row = 0; row < free_indices.size(); ++row)
        for (Index column = 0; column < full_dimension_; ++column)
          if (mass_matrix_.el(free_indices[row], column) != 0.0)
            sparsity.add(row, column);

      coupling_sparsity_.copy_from(sparsity);
      coupling_matrix_.reinit(coupling_sparsity_);

      for (std::size_t row = 0; row < free_indices.size(); ++row)
        for (Index column = 0; column < full_dimension_; ++column)
          {
            const double value = mass_matrix_.el(free_indices[row], column);
            if (value != 0.0)
              coupling_matrix_.set(row, column, value);
          }
    }

    static void
    require_size(const Vector & vector,
                 const std::size_t expected,
                 const char *const name)
    {
      if (static_cast<std::size_t>(vector.size()) != expected)
        throw std::invalid_argument(std::string("Problem B ") + name +
                                    " has the wrong dimension");
    }

    const std::size_t full_dimension_;
    const std::size_t free_dimension_;
    dealii::SparsityPattern mass_sparsity_;
    Matrix             mass_matrix_;
    dealii::SparsityPattern coupling_sparsity_;
    Matrix             coupling_matrix_;
  };
} // namespace external_dealii_step4
