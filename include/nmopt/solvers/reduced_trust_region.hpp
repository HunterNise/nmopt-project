#pragma once

#include "nmopt/solvers/reduced_search.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace nmopt::solvers
{
  enum class ReducedTrustRegionStoppingReason
  {
    gradient_tolerance,
    relative_gradient_tolerance,
    objective_change_tolerance,
    step_tolerance,
    maximum_iterations,
    radius_too_small,
    trial_limit
  };

  inline const char *
  reduced_trust_region_stopping_reason_name(
    const ReducedTrustRegionStoppingReason reason)
  {
    switch (reason)
      {
        case ReducedTrustRegionStoppingReason::gradient_tolerance:
          return "gradient_tolerance";
        case ReducedTrustRegionStoppingReason::relative_gradient_tolerance:
          return "relative_gradient_tolerance";
        case ReducedTrustRegionStoppingReason::objective_change_tolerance:
          return "objective_change_tolerance";
        case ReducedTrustRegionStoppingReason::step_tolerance:
          return "step_tolerance";
        case ReducedTrustRegionStoppingReason::maximum_iterations:
          return "maximum_iterations";
        case ReducedTrustRegionStoppingReason::radius_too_small:
          return "radius_too_small";
        case ReducedTrustRegionStoppingReason::trial_limit:
          return "trial_limit";
      }
    return "unknown";
  }

  struct ReducedTrustRegionParameters
  {
    unsigned int maximum_iterations = 100;
    unsigned int maximum_trials_per_iteration = 10;
    double       gradient_tolerance = 1e-8;
    // A zero value disables the optional relative/objective/step criteria.
    double relative_gradient_tolerance = 0.0;
    double objective_change_tolerance = 0.0;
    double step_tolerance = 0.0;
    double initial_radius = 1.0;
    double minimum_radius = 1e-12;
    double maximum_radius = 1e6;
    double acceptance_threshold = 0.1;
    double shrink_threshold = 0.25;
    double expansion_threshold = 0.75;
    double shrink_factor = 0.25;
    double expansion_factor = 2.0;
  };

  template <typename Backend>
  struct ReducedTrustRegionResultT
  {
    contract::PrimalBlockT<Backend>   control;
    contract::ReducedEvaluationT<Backend> final_evaluation;
    std::vector<double>               objective_history;
    std::vector<double>               gradient_norm_history;
    std::vector<double>               relative_gradient_norm_history;
    std::vector<double>               radius_history;
    std::vector<double>               step_norm_history;
    std::vector<double>               predicted_reduction_history;
    std::vector<double>               actual_reduction_history;
    std::vector<double>               reduction_ratio_history;
    std::vector<bool>                 accepted_history;
    std::size_t                       accepted_iterations;
    std::size_t                       trial_count;
    std::size_t                       state_solve_count;
    std::size_t                       adjoint_solve_count;
    std::size_t                       metric_solve_count;
    std::size_t                       hessian_action_count;
    ReducedTrustRegionStoppingReason  stopping_reason;
  };

  using ReducedTrustRegionResult =
    ReducedTrustRegionResultT<contract::DenseBackend>;

  // A backend-parametric Cauchy trust-region method for the unconstrained
  // reduced DTO boundary. The quadratic model is matrix-free: only H(u)g is
  // requested, where g = G^{-1}j'. Projection and Newton-subproblem solves
  // remain separate future boundaries.
  template <typename Backend>
  class ReducedTrustRegionSolverT
  {
  public:
    using Primal = contract::PrimalBlockT<Backend>;
    using Covector = contract::CovectorBlockT<Backend>;
    using Evaluation = contract::ReducedEvaluationT<Backend>;
    using Result = ReducedTrustRegionResultT<Backend>;

    ReducedTrustRegionSolverT(
      const contract::ReducedDTOT<Backend> &   reduced,
      const contract::MetricT<Backend> &       metric,
      const contract::ReducedHessianT<Backend> &hessian,
      ReducedTrustRegionParameters              parameters = {})
      : reduced_(reduced)
      , metric_(metric)
      , hessian_(hessian)
      , parameters_(parameters)
    {
      validate_parameters();
      contract::require(metric_.layout()->compatible_with(*hessian_.layout()),
                        "Trust-region metric does not match reduced Hessian");
    }

    Result
    solve(const Primal &initial_control) const
    {
      contract::require(initial_control.layout()->compatible_with(*metric_.layout()),
                        "Trust-region initial control does not match metric");
      contract::require(initial_control.layout()->compatible_with(*hessian_.layout()),
                        "Trust-region initial control does not match Hessian");

      std::size_t state_solve_count = 0;
      std::size_t adjoint_solve_count = 0;
      std::size_t metric_solve_count = 0;
      std::size_t hessian_action_count = 0;
      std::size_t accepted_iterations = 0;
      std::size_t trial_count = 0;
      double       initial_gradient_norm = 0.0;
      bool         have_initial_gradient_norm = false;
      double       radius = parameters_.initial_radius;

      Primal current_control = initial_control;
      Evaluation current_evaluation = reduced_.evaluate(current_control);
      ++state_solve_count;
      ++adjoint_solve_count;

      std::vector<double> objective_history{current_evaluation.objective_value};
      std::vector<double> gradient_norm_history;
      std::vector<double> relative_gradient_norm_history;
      std::vector<double> radius_history;
      std::vector<double> step_norm_history;
      std::vector<double> predicted_reduction_history;
      std::vector<double> actual_reduction_history;
      std::vector<double> reduction_ratio_history;
      std::vector<bool>   accepted_history;

      const auto finish =
        [&](const ReducedTrustRegionStoppingReason stopping_reason) {
          return result(std::move(current_control),
                        std::move(current_evaluation),
                        std::move(objective_history),
                        std::move(gradient_norm_history),
                        std::move(relative_gradient_norm_history),
                        std::move(radius_history),
                        std::move(step_norm_history),
                        std::move(predicted_reduction_history),
                        std::move(actual_reduction_history),
                        std::move(reduction_ratio_history),
                        std::move(accepted_history),
                        accepted_iterations,
                        trial_count,
                        state_solve_count,
                        adjoint_solve_count,
                        metric_solve_count,
                        hessian_action_count,
                        stopping_reason);
        };

      for (;;)
        {
          const ReducedMetricGradientT<Backend> metric_gradient =
            make_metric_gradient(current_evaluation.reduced_derivative, metric_);
          ++metric_solve_count;
          const double gradient_norm = metric_gradient.norm;
          if (!have_initial_gradient_norm)
            {
              initial_gradient_norm = gradient_norm;
              have_initial_gradient_norm = true;
            }
          const double relative_gradient_norm =
            initial_gradient_norm > 0.0
              ? gradient_norm / initial_gradient_norm
              : 0.0;
          gradient_norm_history.push_back(gradient_norm);
          relative_gradient_norm_history.push_back(relative_gradient_norm);

          if (gradient_norm <= parameters_.gradient_tolerance)
            return finish(ReducedTrustRegionStoppingReason::gradient_tolerance);

          if (parameters_.relative_gradient_tolerance > 0.0 &&
              relative_gradient_norm <=
                parameters_.relative_gradient_tolerance)
            return finish(
              ReducedTrustRegionStoppingReason::relative_gradient_tolerance);

          if (accepted_iterations == parameters_.maximum_iterations)
            return finish(ReducedTrustRegionStoppingReason::maximum_iterations);

          const Covector hessian_gradient =
            hessian_.apply(current_control, metric_gradient.gradient);
          ++hessian_action_count;
          contract::require(hessian_gradient.layout()->compatible_with(
                              *metric_.layout()),
                            "Trust-region Hessian returned an incompatible covector layout");
          const double gradient_pairing =
            contract::pair(current_evaluation.reduced_derivative,
                           metric_gradient.gradient);
          const double curvature =
            contract::pair(hessian_gradient, metric_gradient.gradient);
          contract::require(std::isfinite(gradient_pairing) &&
                              gradient_pairing > 0.0,
                            "Trust-region gradient pairing is not positive");
          contract::require(std::isfinite(curvature) && curvature > 0.0,
                            "Trust-region model has non-positive curvature");

          const double unconstrained_cauchy_step =
            gradient_pairing / curvature;
          bool accepted = false;
          for (unsigned int trial = 0;
               trial < parameters_.maximum_trials_per_iteration;
               ++trial)
            {
              if (radius < parameters_.minimum_radius)
                return finish(ReducedTrustRegionStoppingReason::radius_too_small);

              const double step_scale = std::min(
                unconstrained_cauchy_step,
                radius / gradient_norm);
              const double step_norm = step_scale * gradient_norm;
              const double predicted_reduction =
                step_scale * gradient_pairing -
                0.5 * step_scale * step_scale * curvature;
              contract::require(std::isfinite(step_scale) && step_scale > 0.0 &&
                                  std::isfinite(step_norm) && step_norm > 0.0,
                                "Trust-region Cauchy step is invalid");
              contract::require(std::isfinite(predicted_reduction) &&
                                  predicted_reduction > 0.0,
                                "Trust-region model predicted non-positive reduction");

              Primal trial_control = current_control;
              add_scaled_primal(trial_control,
                                -step_scale,
                                metric_gradient.gradient);
              Evaluation trial_evaluation = reduced_.evaluate(trial_control);
              ++state_solve_count;
              ++adjoint_solve_count;
              const double actual_reduction =
                current_evaluation.objective_value -
                trial_evaluation.objective_value;
              const double reduction_ratio =
                actual_reduction / predicted_reduction;
              const bool trial_accepted =
                std::isfinite(actual_reduction) &&
                std::isfinite(reduction_ratio) &&
                reduction_ratio >= parameters_.acceptance_threshold;

              radius_history.push_back(radius);
              step_norm_history.push_back(step_norm);
              predicted_reduction_history.push_back(predicted_reduction);
              actual_reduction_history.push_back(actual_reduction);
              reduction_ratio_history.push_back(reduction_ratio);
              accepted_history.push_back(trial_accepted);
              ++trial_count;

              if (trial_accepted)
                {
                  current_control = std::move(trial_control);
                  current_evaluation = std::move(trial_evaluation);
                  objective_history.push_back(current_evaluation.objective_value);
                  ++accepted_iterations;
                  accepted = true;

                  if (reduction_ratio > parameters_.expansion_threshold &&
                      step_norm >= 0.9 * radius)
                    radius = std::min(parameters_.maximum_radius,
                                      parameters_.expansion_factor * radius);
                  else if (reduction_ratio < parameters_.shrink_threshold)
                    radius = std::max(parameters_.minimum_radius,
                                      parameters_.shrink_factor * radius);

                  if (parameters_.objective_change_tolerance > 0.0 &&
                      actual_reduction <=
                        parameters_.objective_change_tolerance)
                    return finish(
                      ReducedTrustRegionStoppingReason::objective_change_tolerance);

                  if (parameters_.step_tolerance > 0.0 &&
                      step_norm <= parameters_.step_tolerance)
                    return finish(ReducedTrustRegionStoppingReason::step_tolerance);
                  break;
                }

              radius = std::max(parameters_.minimum_radius,
                                parameters_.shrink_factor * radius);
              if (radius <= parameters_.minimum_radius)
                return finish(ReducedTrustRegionStoppingReason::radius_too_small);
            }

          if (!accepted)
            return finish(ReducedTrustRegionStoppingReason::trial_limit);
        }
    }

  private:
    void
    validate_parameters() const
    {
      contract::require(parameters_.maximum_iterations > 0,
                        "Trust-region maximum iterations must be positive");
      contract::require(parameters_.maximum_trials_per_iteration > 0,
                        "Trust-region trial limit must be positive");
      contract::require(parameters_.gradient_tolerance > 0.0,
                        "Trust-region gradient tolerance must be positive");
      contract::require(parameters_.relative_gradient_tolerance >= 0.0,
                        "Trust-region relative gradient tolerance must be nonnegative");
      contract::require(parameters_.objective_change_tolerance >= 0.0,
                        "Trust-region objective-change tolerance must be nonnegative");
      contract::require(parameters_.step_tolerance >= 0.0,
                        "Trust-region step tolerance must be nonnegative");
      contract::require(parameters_.minimum_radius > 0.0 &&
                          parameters_.initial_radius >= parameters_.minimum_radius &&
                          parameters_.maximum_radius >= parameters_.initial_radius,
                        "Trust-region radii are not ordered and positive");
      contract::require(parameters_.acceptance_threshold >= 0.0 &&
                          parameters_.acceptance_threshold <
                            parameters_.shrink_threshold &&
                          parameters_.shrink_threshold <
                            parameters_.expansion_threshold &&
                          parameters_.expansion_threshold < 1.0,
                        "Trust-region reduction thresholds are not ordered");
      contract::require(parameters_.shrink_factor > 0.0 &&
                          parameters_.shrink_factor < 1.0,
                        "Trust-region shrink factor must lie in (0, 1)");
      contract::require(parameters_.expansion_factor > 1.0,
                        "Trust-region expansion factor must exceed one");
    }

    static Result
    result(Primal                         control,
           Evaluation                     final_evaluation,
           std::vector<double>            objective_history,
           std::vector<double>            gradient_norm_history,
           std::vector<double>            relative_gradient_norm_history,
           std::vector<double>            radius_history,
           std::vector<double>            step_norm_history,
           std::vector<double>            predicted_reduction_history,
           std::vector<double>            actual_reduction_history,
           std::vector<double>            reduction_ratio_history,
           std::vector<bool>              accepted_history,
           const std::size_t              accepted_iterations,
           const std::size_t              trial_count,
           const std::size_t              state_solve_count,
           const std::size_t              adjoint_solve_count,
           const std::size_t              metric_solve_count,
           const std::size_t              hessian_action_count,
           const ReducedTrustRegionStoppingReason stopping_reason)
    {
      return {std::move(control),
              std::move(final_evaluation),
              std::move(objective_history),
              std::move(gradient_norm_history),
              std::move(relative_gradient_norm_history),
              std::move(radius_history),
              std::move(step_norm_history),
              std::move(predicted_reduction_history),
              std::move(actual_reduction_history),
              std::move(reduction_ratio_history),
              std::move(accepted_history),
              accepted_iterations,
              trial_count,
              state_solve_count,
              adjoint_solve_count,
              metric_solve_count,
              hessian_action_count,
              stopping_reason};
    }

    const contract::ReducedDTOT<Backend> &   reduced_;
    const contract::MetricT<Backend> &       metric_;
    const contract::ReducedHessianT<Backend> &hessian_;
    ReducedTrustRegionParameters             parameters_;
  };

  using ReducedTrustRegionSolver =
    ReducedTrustRegionSolverT<contract::DenseBackend>;
} // namespace nmopt::solvers
