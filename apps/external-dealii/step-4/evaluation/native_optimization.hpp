#pragma once

#include "native_reduced.hpp"
#include "optimization_policy.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace external_dealii_step4
{
  enum class NativeOptimizationStoppingReason
  {
    gradient_tolerance,
    maximum_iterations,
    line_search_failure
  };

  inline const char *
  native_optimization_stopping_reason_name(
    const NativeOptimizationStoppingReason reason)
  {
    switch (reason)
      {
        case NativeOptimizationStoppingReason::gradient_tolerance:
          return "gradient_tolerance";
        case NativeOptimizationStoppingReason::maximum_iterations:
          return "maximum_iterations";
        case NativeOptimizationStoppingReason::line_search_failure:
          return "line_search_failure";
      }
    return "unknown";
  }

  struct OptimizationTrialRecord
  {
    std::size_t trial = 0;
    std::size_t iteration = 0;
    double      step_length = 0.0;
    double      objective_value = std::numeric_limits<double>::quiet_NaN();
    double      actual_slope = std::numeric_limits<double>::quiet_NaN();
    double      sufficient_decrease_bound =
      std::numeric_limits<double>::quiet_NaN();
    bool objective_finite = false;
    bool slope_negative = false;
    bool accepted = false;
  };

  struct OptimizationAcceptedIterationRecord
  {
    std::size_t iteration = 0;
    double      objective_before = 0.0;
    double      objective_after = 0.0;
    double      objective_change = 0.0;
    double      requested_step_length = 0.0;
    double      actual_step_norm = 0.0;
    double      actual_slope = 0.0;
    double      gradient_norm = 0.0;
    std::size_t trial_count = 0;
  };

  struct NativeOptimizationResult
  {
    NativeReduced::ValueEvaluation      value;
    NativeReduced::DerivativeEvaluation derivative;
    std::vector<double>                  objective_history;
    std::vector<double>                  gradient_norm_history;
    std::vector<OptimizationAcceptedIterationRecord>
      accepted_iterations;
    std::vector<OptimizationTrialRecord> trial_records;
    std::size_t accepted_iteration_count = 0;
    std::size_t line_search_trial_count = 0;
    NativeOptimizationStoppingReason stopping_reason;
  };

  class NativeArmijoSolver final
  {
  public:
    NativeArmijoSolver(NativeReduced &      reduced,
                       OptimizationPolicy  policy = frozen_optimization_policy())
      : reduced_(reduced)
      , policy_(std::move(policy))
    {
      validate_policy();
    }

    NativeOptimizationResult
    solve(const NativeReduced::Vector &initial_control) const
    {
      auto current_value = reduced_.evaluate_value(initial_control);
      auto current_derivative = reduced_.augment_derivative(current_value);

      std::vector<double> objective_history{current_value.objective};
      std::vector<double> gradient_norm_history;
      std::vector<OptimizationAcceptedIterationRecord> accepted_iterations;
      std::vector<OptimizationTrialRecord> trial_records;
      std::size_t accepted_iteration_count = 0;
      std::size_t line_search_trial_count = 0;

      const auto finish = [&](const NativeOptimizationStoppingReason reason) {
        return NativeOptimizationResult{
          std::move(current_value),
          std::move(current_derivative),
          std::move(objective_history),
          std::move(gradient_norm_history),
          std::move(accepted_iterations),
          std::move(trial_records),
          accepted_iteration_count,
          line_search_trial_count,
          reason};
      };

      for (;;)
        {
          require_finite(current_value.objective, "current objective");
          require_finite(current_derivative.reduced_derivative,
                         "current reduced derivative");

          const double gradient_norm =
            current_derivative.reduced_derivative.l2_norm();
          require_finite(gradient_norm, "current gradient norm");
          gradient_norm_history.push_back(gradient_norm);
          if (gradient_norm <= policy_.gradient_tolerance)
            return finish(NativeOptimizationStoppingReason::gradient_tolerance);

          if (accepted_iteration_count == policy_.maximum_iterations)
            return finish(NativeOptimizationStoppingReason::maximum_iterations);

          NativeReduced::Vector direction =
            current_derivative.reduced_derivative;
          direction *= -1.0;
          require_finite(direction, "search direction");
          const double direction_pairing =
            current_derivative.reduced_derivative * direction;
          if (!std::isfinite(direction_pairing) || direction_pairing >= 0.0)
            throw std::runtime_error(
              "native optimization did not produce a finite descent direction");

          double step_length = policy_.initial_step_length;
          bool accepted = false;
          for (unsigned int trial = 0;
               trial < policy_.maximum_line_search_trials;
               ++trial)
            {
              NativeReduced::Vector trial_control = current_value.control;
              trial_control.add(step_length, direction);
              const auto trial_value = reduced_.evaluate_value(trial_control);

              NativeReduced::Vector actual_update = trial_control;
              actual_update.add(-1.0, current_value.control);
              const double actual_slope =
                current_derivative.reduced_derivative * actual_update;
              if (!std::isfinite(actual_slope))
                throw std::runtime_error(
                  "native optimization produced a non-finite actual slope");

              const double sufficient_decrease_bound =
                current_value.objective +
                policy_.armijo_fraction * actual_slope;
              const bool objective_finite =
                std::isfinite(trial_value.objective);
              const bool slope_negative = actual_slope < 0.0;
              accepted = slope_negative && objective_finite &&
                         trial_value.objective <= sufficient_decrease_bound;
              trial_records.push_back({trial,
                                       accepted_iteration_count + 1,
                                       step_length,
                                       trial_value.objective,
                                       actual_slope,
                                       sufficient_decrease_bound,
                                       objective_finite,
                                       slope_negative,
                                       accepted});
              ++line_search_trial_count;

              if (accepted)
                {
                  const double objective_before = current_value.objective;
                  const double actual_step_norm = actual_update.l2_norm();
                  const auto trial_derivative =
                    reduced_.augment_derivative(trial_value);
                  require_finite(trial_derivative.reduced_derivative,
                                 "accepted reduced derivative");

                  current_value = trial_value;
                  current_derivative = trial_derivative;
                  ++accepted_iteration_count;
                  objective_history.push_back(current_value.objective);
                  accepted_iterations.push_back(
                    {accepted_iteration_count,
                     objective_before,
                     current_value.objective,
                     current_value.objective - objective_before,
                     step_length,
                     actual_step_norm,
                     actual_slope,
                     gradient_norm,
                     trial + 1});
                  break;
                }

              if (trial + 1 == policy_.maximum_line_search_trials ||
                  (policy_.minimum_step_length > 0.0 &&
                   step_length <= policy_.minimum_step_length))
                break;

              step_length *= policy_.backtracking_factor;
              if (policy_.minimum_step_length > 0.0)
                step_length =
                  std::max(step_length, policy_.minimum_step_length);
            }

          if (!accepted)
            return finish(NativeOptimizationStoppingReason::line_search_failure);
        }
    }

  private:
    static void
    require_finite(const double value, const char *const name)
    {
      if (!std::isfinite(value))
        throw std::runtime_error(std::string("native optimization received a non-finite ") +
                                 name);
    }

    static void
    require_finite(const NativeReduced::Vector &vector,
                   const char *const           name)
    {
      for (unsigned int index = 0; index < vector.size(); ++index)
        if (!std::isfinite(vector[index]))
          throw std::runtime_error(
            std::string("native optimization received a non-finite ") + name);
    }

    void
    validate_policy() const
    {
      if (policy_.maximum_iterations == 0 ||
          policy_.maximum_line_search_trials == 0)
        throw std::invalid_argument(
          "native optimization iteration and trial limits must be positive");
      if (policy_.gradient_tolerance <= 0.0 ||
          policy_.relative_gradient_tolerance < 0.0 ||
          policy_.objective_change_tolerance < 0.0 ||
          policy_.step_tolerance < 0.0)
        throw std::invalid_argument(
          "native optimization tolerances have invalid values");
      if (policy_.objective_target.has_value() &&
          !std::isfinite(*policy_.objective_target))
        throw std::invalid_argument(
          "native optimization objective target must be finite");
      if (policy_.initial_step_length <= 0.0 ||
          policy_.minimum_step_length < 0.0 ||
          policy_.minimum_step_length > policy_.initial_step_length ||
          policy_.armijo_fraction <= 0.0 ||
          policy_.armijo_fraction >= 1.0 ||
          policy_.backtracking_factor <= 0.0 ||
          policy_.backtracking_factor >= 1.0)
        throw std::invalid_argument(
          "native optimization line-search parameters have invalid values");
    }

    NativeReduced &     reduced_;
    OptimizationPolicy policy_;
  };
} // namespace external_dealii_step4
