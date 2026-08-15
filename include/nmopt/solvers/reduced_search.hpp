#pragma once

#include "nmopt/contract/metric_constraint.hpp"

#include <cmath>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace nmopt::solvers
{
  enum class ReducedStoppingReason
  {
    gradient_tolerance,
    maximum_iterations,
    line_search_failure
  };

  inline const char *
  reduced_stopping_reason_name(const ReducedStoppingReason reason)
  {
    switch (reason)
      {
        case ReducedStoppingReason::gradient_tolerance:
          return "gradient_tolerance";
        case ReducedStoppingReason::maximum_iterations:
          return "maximum_iterations";
        case ReducedStoppingReason::line_search_failure:
          return "line_search_failure";
      }
    return "unknown";
  }

  inline const char *
  reduced_gradient_stopping_reason_name(const ReducedStoppingReason reason)
  {
    return reduced_stopping_reason_name(reason);
  }

  struct ReducedSolverParameters
  {
    unsigned int maximum_iterations = 100;
    unsigned int maximum_line_search_trials = 20;
    double       gradient_tolerance = 1e-8;
    double       initial_step_length = 1.0;
    double       armijo_fraction = 1e-4;
    double       backtracking_factor = 0.5;
  };

  template <typename Backend>
  struct ReducedSearchDirectionT
  {
    contract::PrimalBlockT<Backend> direction;
    double                          gradient_norm;
    double                          directional_derivative;
  };

  using ReducedSearchDirection =
    ReducedSearchDirectionT<contract::DenseBackend>;

  template <typename Backend>
  struct ReducedSolverResultT
  {
    contract::PrimalBlockT<Backend> control;
    std::vector<double>             objective_history;
    std::vector<double>             gradient_norm_history;
    std::size_t                     accepted_iterations;
    std::size_t                     line_search_trial_count;
    std::size_t                     state_solve_count;
    std::size_t                     adjoint_solve_count;
    ReducedStoppingReason           stopping_reason;
    std::vector<double>             step_length_history;
    std::vector<double>             objective_change_history;
    std::size_t                     metric_solve_count;
  };

  using ReducedSolverResult = ReducedSolverResultT<contract::DenseBackend>;

  // Compatibility aliases retain the v0 gradient-solver names while the
  // result and direction protocol is shared by later reduced methods.
  using ReducedGradientParameters = ReducedSolverParameters;
  using ReducedGradientStoppingReason = ReducedStoppingReason;

  template <typename Backend>
  using ReducedGradientResultT = ReducedSolverResultT<Backend>;

  using ReducedGradientResult =
    ReducedGradientResultT<contract::DenseBackend>;

  template <typename Backend>
  inline void
  add_scaled_primal(contract::PrimalBlockT<Backend> &      target,
                    const double                           factor,
                    const contract::PrimalBlockT<Backend> &source)
  {
    contract::require_compatible(target,
                                  source,
                                  "Reduced search update has incompatible layouts");
    for (std::size_t block = 0; block < target.n_blocks(); ++block)
      target.add_scaled_block(block, factor, source.block(block));
  }

  template <typename Backend>
  inline void
  scale_primal(contract::PrimalBlockT<Backend> &target, const double factor)
  {
    for (std::size_t block = 0; block < target.n_blocks(); ++block)
      target.scale_block(block, factor);
  }

  template <typename Backend>
  ReducedSearchDirectionT<Backend>
  make_steepest_descent_direction(
    const contract::CovectorBlockT<Backend> &reduced_derivative,
    const contract::MetricT<Backend> &       metric)
  {
    contract::require(reduced_derivative.layout()->compatible_with(
                        *metric.layout()),
                      "Reduced search derivative does not match metric");

    contract::PrimalBlockT<Backend> gradient =
      metric.inverse_apply(reduced_derivative);
    const contract::CovectorBlockT<Backend> metric_gradient =
      metric.apply(gradient);
    const double metric_norm_squared =
      contract::pair(metric_gradient, gradient);
    contract::require(std::isfinite(metric_norm_squared) &&
                        metric_norm_squared >= 0.0,
                      "Reduced search metric returned an invalid squared norm");

    contract::PrimalBlockT<Backend> direction = gradient;
    scale_primal(direction, -1.0);
    const double directional_derivative =
      contract::pair(reduced_derivative, direction);
    contract::require(std::isfinite(directional_derivative) &&
                        directional_derivative <= 0.0,
                      "Reduced search metric returned a non-descent direction");

    return {std::move(direction),
            std::sqrt(metric_norm_squared),
            directional_derivative};
  }
} // namespace nmopt::solvers
