#pragma once

#include <deal.II/lac/vector.h>

#include <cmath>
#include <cstddef>
#include <vector>

namespace external_dealii_step4
{
  namespace scenario
  {
    using Vector = dealii::Vector<double>;

    inline constexpr std::size_t dimension = 289;

    inline Vector
    constant_vector()
    {
      Vector value(dimension);
      value = 1.0 / std::sqrt(static_cast<double>(dimension));
      return value;
    }

    inline Vector
    ramp_vector()
    {
      Vector value(dimension);
      const double scale = std::sqrt(static_cast<double>(dimension));
      for (std::size_t index = 0; index < dimension; ++index)
        value[index] =
          (static_cast<double>(index) / static_cast<double>(dimension - 1) -
           0.5) /
          scale;
      return value;
    }

    inline Vector
    alternating_vector()
    {
      Vector value(dimension);
      const double scale = std::sqrt(static_cast<double>(dimension));
      for (std::size_t index = 0; index < dimension; ++index)
        value[index] = (index % 2 == 0 ? 1.0 : -1.0) / scale;
      return value;
    }

    inline Vector
    scaled(Vector value, const double factor)
    {
      value *= factor;
      return value;
    }

    inline Vector
    normalized(Vector value)
    {
      value /= value.l2_norm();
      return value;
    }

    inline std::vector<Vector>
    reduced_controls()
    {
      const auto constant   = constant_vector();
      const auto ramp       = ramp_vector();
      const auto alternating = alternating_vector();
      Vector zero(dimension);
      zero = 0.0;
      return {zero,
              scaled(constant, 0.1),
              scaled(ramp, 0.1),
              scaled(alternating, 0.1)};
    }

    inline Vector
    normalized_ramp_direction()
    {
      return normalized(ramp_vector());
    }

    inline Vector
    normalized_alternating_direction()
    {
      return normalized(alternating_vector());
    }
  } // namespace scenario
} // namespace external_dealii_step4
