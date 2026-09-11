#pragma once

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/vector.h>

#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace external_dealii_step4
{
  template <int dim>
  class ProblemBCoordinates final
  {
  public:
    using Index          = dealii::types::global_dof_index;
    using Vector         = dealii::Vector<double>;
    using BoundaryValues = std::map<Index, double>;

    ProblemBCoordinates(const dealii::DoFHandler<dim> &dof_handler,
                        const BoundaryValues &         boundary_values)
      : full_dimension_(dof_handler.n_dofs())
      , lifting_(full_dimension_)
    {
      if (full_dimension_ == 0)
        throw std::invalid_argument(
          "Problem B coordinate map needs an initialized DoF handler");

      if (full_dimension_ >
          static_cast<Index>(std::numeric_limits<unsigned int>::max()))
        throw std::invalid_argument(
          "Problem B coordinate map exceeds the vector dimension limit");

      lifting_ = 0.0;
      free_indices_.reserve(full_dimension_);

      for (const auto &[index, value] : boundary_values)
        {
          if (index >= full_dimension_)
            throw std::invalid_argument(
              "Problem B boundary data contains an invalid DoF index");
          if (!std::isfinite(value))
            throw std::invalid_argument(
              "Problem B boundary data contains a non-finite value");
          lifting_[index] = value;
        }

      for (Index index = 0; index < full_dimension_; ++index)
        if (boundary_values.find(index) == boundary_values.end())
          free_indices_.push_back(index);

      if (free_indices_.empty())
        throw std::invalid_argument(
          "Problem B coordinate map has no free DoFs");
    }

    std::size_t
    full_dimension() const
    {
      return full_dimension_;
    }

    std::size_t
    free_dimension() const
    {
      return free_indices_.size();
    }

    std::size_t
    boundary_dimension() const
    {
      return full_dimension_ - free_dimension();
    }

    const std::vector<Index> &
    free_indices() const
    {
      return free_indices_;
    }

    const Vector &
    lifting() const
    {
      return lifting_;
    }

    Vector
    restrict(const Vector &full_vector) const
    {
      require_size(full_vector, full_dimension_, "full vector");

      Vector free_vector(free_dimension());
      for (std::size_t index = 0; index < free_indices_.size(); ++index)
        free_vector[index] = full_vector[free_indices_[index]];
      return free_vector;
    }

    Vector
    reconstruct(const Vector &free_vector) const
    {
      require_size(free_vector, free_dimension(), "free vector");

      Vector full_vector = lifting_;
      for (std::size_t index = 0; index < free_indices_.size(); ++index)
        full_vector[free_indices_[index]] = free_vector[index];
      return full_vector;
    }

    Vector
    embed_free(const Vector &free_vector) const
    {
      require_size(free_vector, free_dimension(), "free vector");

      Vector full_vector(full_dimension_);
      full_vector = 0.0;
      for (std::size_t index = 0; index < free_indices_.size(); ++index)
        full_vector[free_indices_[index]] = free_vector[index];
      return full_vector;
    }

  private:
    static void
    require_size(const Vector & vector,
                 const std::size_t expected,
                 const char *const name)
    {
      if (static_cast<std::size_t>(vector.size()) != expected)
        throw std::invalid_argument(std::string("Problem B ") + name +
                                    " has the wrong dimension");
    }

    const std::size_t    full_dimension_;
    std::vector<Index>   free_indices_;
    Vector               lifting_;
  };
} // namespace external_dealii_step4
