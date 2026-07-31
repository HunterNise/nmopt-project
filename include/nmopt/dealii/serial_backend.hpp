#pragma once

#include "nmopt/contract/linalg.hpp"

#include <deal.II/lac/vector.h>

#include <cstddef>

namespace nmopt::dealii_backend
{
  // Serial deal.II vector policy for the v0 lowerer. Distributed vectors and
  // their ownership/ghost policy are a separate backend extension.
  struct SerialBackend
  {
    using Vector = dealii::Vector<double>;

    static Vector
    zeros(const std::size_t size)
    {
      return Vector(size);
    }

    static std::size_t
    size(const Vector &vector)
    {
      return vector.size();
    }

    static double
    dot(const Vector &left, const Vector &right)
    {
      return left * right;
    }

    static void
    add_scaled(Vector &target, const double factor, const Vector &source)
    {
      target.add(factor, source);
    }

    static void
    scale(Vector &target, const double factor)
    {
      target *= factor;
    }
  };
} // namespace nmopt::dealii_backend
