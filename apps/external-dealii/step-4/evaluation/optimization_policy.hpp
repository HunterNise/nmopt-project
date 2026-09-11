#pragma once

#include <optional>

namespace external_dealii_step4
{
  struct OptimizationPolicy
  {
    unsigned int       maximum_iterations          = 5000;
    unsigned int       maximum_line_search_trials   = 30;
    double             gradient_tolerance           = 1.0e-6;
    double             relative_gradient_tolerance = 0.0;
    double             objective_change_tolerance  = 0.0;
    double             step_tolerance               = 0.0;
    std::optional<double> objective_target;
    double             initial_step_length         = 1.0;
    double             minimum_step_length         = 0.0;
    double             armijo_fraction              = 1.0e-4;
    double             backtracking_factor          = 0.5;
  };

  inline OptimizationPolicy
  frozen_optimization_policy()
  {
    return {};
  }
} // namespace external_dealii_step4
