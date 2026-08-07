#pragma once

#include "nmopt/contract/linalg.hpp"

#include <deal.II/lac/vector.h>

#include <cstddef>
#include <limits>

namespace nmopt::dealii_backend
{
  // Serial deal.II vector policy for the v0 lowerer. Distributed vectors and
  // their ownership/ghost policy are a separate backend extension.
  struct SerialBackend
  {
    using Vector = dealii::Vector<double>;
    using NativeSize = Vector::size_type;

    static constexpr bool native_size_is_narrower =
      std::numeric_limits<NativeSize>::digits <
      std::numeric_limits<std::size_t>::digits;

    static constexpr std::size_t
    maximum_native_size() noexcept
    {
      if constexpr (native_size_is_narrower)
        return static_cast<std::size_t>(
          std::numeric_limits<NativeSize>::max());
      return std::numeric_limits<std::size_t>::max();
    }

    static NativeSize
    checked_native_size(const std::size_t size)
    {
      contract::require(
        size <= maximum_native_size(),
        "Serial deal.II vector size exceeds its native range");
      return static_cast<NativeSize>(size);
    }

    static Vector
    zeros(const std::size_t size)
    {
      return Vector(checked_native_size(size));
    }

    static std::size_t
    size(const Vector &vector)
    {
      return static_cast<std::size_t>(vector.size());
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
