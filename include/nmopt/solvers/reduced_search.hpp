#pragma once

#include "nmopt/contract/metric_constraint.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
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
  struct ReducedMetricGradientT
  {
    contract::PrimalBlockT<Backend> gradient;
    double                          norm;
  };

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
  ReducedMetricGradientT<Backend>
  make_metric_gradient(const contract::CovectorBlockT<Backend> &reduced_derivative,
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

    return {std::move(gradient), std::sqrt(metric_norm_squared)};
  }

  template <typename Backend>
  ReducedSearchDirectionT<Backend>
  make_steepest_descent_direction(
    const contract::CovectorBlockT<Backend> &reduced_derivative,
    const contract::MetricT<Backend> &       metric)
  {
    ReducedMetricGradientT<Backend> metric_gradient =
      make_metric_gradient(reduced_derivative, metric);

    contract::PrimalBlockT<Backend> direction = metric_gradient.gradient;
    scale_primal(direction, -1.0);
    const double directional_derivative =
      contract::pair(reduced_derivative, direction);
    contract::require(std::isfinite(directional_derivative) &&
                        directional_derivative <= 0.0,
                      "Reduced search metric returned a non-descent direction");

    return {std::move(direction),
            metric_gradient.norm,
            directional_derivative};
  }

  template <typename Backend>
  class SteepestDescentDirectionPolicyT
  {
  public:
    using Direction = ReducedSearchDirectionT<Backend>;

    void
    reset()
    {}

    Direction
    next(const contract::CovectorBlockT<Backend> &reduced_derivative,
         const contract::MetricT<Backend> &       metric)
    {
      return make_steepest_descent_direction(reduced_derivative, metric);
    }
  };

  using SteepestDescentDirectionPolicy =
    SteepestDescentDirectionPolicyT<contract::DenseBackend>;

  struct NonlinearConjugateGradientParameters
  {
    // Zero selects one restart interval equal to the total coefficient count
    // of the current reduced layout.
    unsigned int restart_interval = 0;
    double       curvature_tolerance = 1e-14;
  };

  // The selected nonlinear-CG policy is Polak–Ribière+ in the declared metric.
  // It stores the previous covector, metric gradient, and primal direction so
  // primal/dual layout compatibility is checked at every update.
  template <typename Backend>
  class NonlinearConjugateGradientDirectionPolicyT
  {
  public:
    using Direction = ReducedSearchDirectionT<Backend>;
    using Primal = contract::PrimalBlockT<Backend>;
    using Covector = contract::CovectorBlockT<Backend>;

    NonlinearConjugateGradientDirectionPolicyT(
      NonlinearConjugateGradientParameters parameters = {})
      : parameters_(parameters)
    {
      contract::require(parameters_.curvature_tolerance > 0.0,
                        "Nonlinear CG curvature tolerance must be positive");
    }

    void
    reset()
    {
      previous_derivative_.reset();
      previous_gradient_.reset();
      previous_direction_.reset();
      directions_since_restart_ = 0;
      restart_count_ = 0;
    }

    Direction
    next(const Covector &reduced_derivative, const contract::MetricT<Backend> &metric)
    {
      const ReducedMetricGradientT<Backend> current_gradient =
        make_metric_gradient(reduced_derivative, metric);
      const bool has_history = previous_derivative_.has_value();
      bool       restart = !has_history;
      double     beta = 0.0;

      if (has_history)
        {
          contract::require(reduced_derivative.layout()->compatible_with(
                              *previous_derivative_->layout()),
                            "Nonlinear CG derivative history has an incompatible layout");
          contract::require(current_gradient.gradient.layout()->compatible_with(
                              *previous_gradient_->layout()),
                            "Nonlinear CG gradient history has an incompatible layout");
          contract::require(previous_direction_->layout()->compatible_with(
                              *current_gradient.gradient.layout()),
                            "Nonlinear CG direction history has an incompatible layout");

          const std::size_t restart_interval = effective_restart_interval(
            current_gradient.gradient.layout());
          restart = directions_since_restart_ >= restart_interval;

          const double denominator =
            contract::pair(*previous_derivative_, *previous_gradient_);
          const double numerator =
            contract::pair(reduced_derivative, current_gradient.gradient) -
            contract::pair(*previous_derivative_, current_gradient.gradient);
          const double scale = std::max(1.0, std::abs(denominator));
          if (!std::isfinite(denominator) ||
              denominator <= parameters_.curvature_tolerance * scale ||
              !std::isfinite(numerator))
            restart = true;
          else
            {
              beta = numerator / denominator;
              if (!std::isfinite(beta) || beta <= 0.0)
                restart = true;
            }
        }

      Primal direction = current_gradient.gradient;
      scale_primal(direction, -1.0);
      if (!restart)
        add_scaled_primal(direction, beta, *previous_direction_);

      double directional_derivative =
        contract::pair(reduced_derivative, direction);
      if (!std::isfinite(directional_derivative) ||
          directional_derivative >= 0.0)
        {
          restart = true;
          direction = current_gradient.gradient;
          scale_primal(direction, -1.0);
          directional_derivative =
            contract::pair(reduced_derivative, direction);
        }

      if (has_history && restart)
        ++restart_count_;
      directions_since_restart_ = restart ? 1 : directions_since_restart_ + 1;
      previous_derivative_ = reduced_derivative;
      previous_gradient_ = current_gradient.gradient;
      previous_direction_ = direction;

      contract::require(std::isfinite(directional_derivative) &&
                          directional_derivative <= 0.0,
                        "Nonlinear CG could not produce a descent direction");
      return {std::move(direction),
              current_gradient.norm,
              directional_derivative};
    }

    std::size_t
    restart_count() const noexcept
    {
      return restart_count_;
    }

  private:
    std::size_t
    effective_restart_interval(const contract::LayoutPtr &layout) const
    {
      if (parameters_.restart_interval > 0)
        return parameters_.restart_interval;

      std::size_t dimension = 0;
      for (std::size_t block = 0; block < layout->n_blocks(); ++block)
        dimension += layout->dimension(block);
      return std::max<std::size_t>(1, dimension);
    }

    NonlinearConjugateGradientParameters parameters_;
    std::optional<Covector>             previous_derivative_;
    std::optional<Primal>               previous_gradient_;
    std::optional<Primal>               previous_direction_;
    std::size_t                         directions_since_restart_ = 0;
    std::size_t                         restart_count_ = 0;
  };

  template <typename Backend>
  using NonlinearConjugateGradientDirectionPolicy =
    NonlinearConjugateGradientDirectionPolicyT<Backend>;

  using NonlinearConjugateGradientDirectionPolicyDense =
    NonlinearConjugateGradientDirectionPolicyT<contract::DenseBackend>;
} // namespace nmopt::solvers
