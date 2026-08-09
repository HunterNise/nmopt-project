#pragma once

namespace nmopt::dealii_backend
{
  struct MetricSolveParameters
  {
    unsigned int maximum_iterations = 1000;
    double       relative_tolerance = 1e-12;
    double       absolute_tolerance = 1e-14;
  };
} // namespace nmopt::dealii_backend
