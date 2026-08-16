#pragma once

#include "nmopt/solvers/reduced_search.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace nmopt::solvers
{
  enum class ReducedTrustRegionStoppingReason
  {
    gradient_tolerance,
    relative_gradient_tolerance,
    stationary,
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
        case ReducedTrustRegionStoppingReason::stationary:
          return "stationary";
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

  enum class ReducedTrustRegionSubproblemMethod
  {
    cauchy,
    truncated_conjugate_gradient
  };

  inline const char *
  reduced_trust_region_subproblem_method_name(
    const ReducedTrustRegionSubproblemMethod method)
  {
    switch (method)
      {
        case ReducedTrustRegionSubproblemMethod::cauchy:
          return "cauchy";
        case ReducedTrustRegionSubproblemMethod::truncated_conjugate_gradient:
          return "truncated_conjugate_gradient";
      }
    return "unknown";
  }

  enum class ReducedTrustRegionSubproblemStatus
  {
    cauchy,
    converged,
    boundary,
    negative_curvature,
    iteration_limit
  };

  inline const char *
  reduced_trust_region_subproblem_status_name(
    const ReducedTrustRegionSubproblemStatus status)
  {
    switch (status)
      {
        case ReducedTrustRegionSubproblemStatus::cauchy:
          return "cauchy";
        case ReducedTrustRegionSubproblemStatus::converged:
          return "converged";
        case ReducedTrustRegionSubproblemStatus::boundary:
          return "boundary";
        case ReducedTrustRegionSubproblemStatus::negative_curvature:
          return "negative_curvature";
        case ReducedTrustRegionSubproblemStatus::iteration_limit:
          return "iteration_limit";
      }
    return "unknown";
  }

  template <typename Backend>
  struct ReducedTrustRegionIterationRecordT
  {
    ReducedAcceptedIterationCommonT<Backend> common;
    double                         radius = 0.0;
    double                         predicted_reduction = 0.0;
    double                         actual_reduction = 0.0;
    double                         reduction_ratio = 0.0;
    ReducedTrustRegionSubproblemStatus subproblem_status;
    std::size_t                    subproblem_iteration_count = 0;
    double                         subproblem_residual_norm = 0.0;
  };

  struct ReducedTrustRegionParameters
  {
    unsigned int maximum_iterations = 100;
    unsigned int maximum_trials_per_iteration = 10;
    double       gradient_tolerance = 1e-8;
    ReducedStoppingCriterion stopping_criterion =
      ReducedStoppingCriterion::automatic;
    // A zero value disables the optional relative/objective/step criteria.
    double relative_gradient_tolerance = 0.0;
    double objective_change_tolerance = 0.0;
    double step_tolerance = 0.0;
    ReducedTrustRegionSubproblemMethod subproblem_method =
      ReducedTrustRegionSubproblemMethod::cauchy;
    unsigned int maximum_subproblem_iterations = 100;
    double subproblem_relative_tolerance = 1e-10;
    double subproblem_absolute_tolerance = 1e-12;
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
    std::vector<ReducedTrustRegionSubproblemStatus>
      subproblem_status_history;
    std::vector<std::size_t>          subproblem_iteration_history;
    std::vector<double>               subproblem_residual_norm_history;
    std::size_t                       accepted_iterations;
    std::size_t                       trial_count;
    std::size_t                       state_solve_count;
    std::size_t                       adjoint_solve_count;
    std::size_t                       metric_solve_count;
    std::size_t                       hessian_action_count;
    ReducedTrustRegionStoppingReason  stopping_reason;
    ReducedTrustRegionParameters       parameters;
    std::string                        policy_name;
    std::vector<ReducedTrustRegionIterationRecordT<Backend>> iteration_records;
  };

  using ReducedTrustRegionResult =
    ReducedTrustRegionResultT<contract::DenseBackend>;

  // A backend-parametric trust-region method for the unconstrained reduced DTO
  // boundary. The default Cauchy step remains available alongside a
  // matrix-free truncated-CG subproblem solver for explicit linear-quadratic
  // Hessians. Projection remains a separate future boundary.
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
      std::vector<ReducedTrustRegionSubproblemStatus>
        subproblem_status_history;
      std::vector<std::size_t> subproblem_iteration_history;
      std::vector<double>      subproblem_residual_norm_history;
      std::vector<ReducedTrustRegionIterationRecordT<Backend>> iteration_records;

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
                        std::move(subproblem_status_history),
                        std::move(subproblem_iteration_history),
                        std::move(subproblem_residual_norm_history),
                        accepted_iterations,
                        trial_count,
                        state_solve_count,
                        adjoint_solve_count,
                        metric_solve_count,
                        hessian_action_count,
                        stopping_reason,
                        std::move(iteration_records));
        };

      for (;;)
        {
          const ReducedIterationWork work_before{
            state_solve_count,
            adjoint_solve_count,
            metric_solve_count,
            hessian_action_count,
            0};
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

          const bool automatic_stopping =
            parameters_.stopping_criterion ==
            ReducedStoppingCriterion::automatic;
          if ((automatic_stopping ||
               parameters_.stopping_criterion ==
                 ReducedStoppingCriterion::gradient_norm) &&
              gradient_norm <= parameters_.gradient_tolerance)
            return finish(ReducedTrustRegionStoppingReason::gradient_tolerance);

          if ((automatic_stopping
                 ? parameters_.relative_gradient_tolerance > 0.0
                 : parameters_.stopping_criterion ==
                     ReducedStoppingCriterion::relative_gradient_norm) &&
              relative_gradient_norm <=
                parameters_.relative_gradient_tolerance)
            return finish(
              ReducedTrustRegionStoppingReason::relative_gradient_tolerance);

          const bool post_step_stopping =
            parameters_.stopping_criterion ==
                ReducedStoppingCriterion::objective_change ||
            parameters_.stopping_criterion ==
              ReducedStoppingCriterion::step_norm;
          if (post_step_stopping &&
              reduced_is_numerically_stationary(gradient_norm))
            return finish(ReducedTrustRegionStoppingReason::stationary);

          if (accepted_iterations == parameters_.maximum_iterations)
            return finish(ReducedTrustRegionStoppingReason::maximum_iterations);

          std::optional<Covector> cauchy_hessian_gradient;
          if (parameters_.subproblem_method ==
              ReducedTrustRegionSubproblemMethod::cauchy)
            {
              cauchy_hessian_gradient =
                hessian_.apply(current_control, metric_gradient.gradient);
              ++hessian_action_count;
              contract::require(
                cauchy_hessian_gradient->layout()->compatible_with(
                  *metric_.layout()),
                "Trust-region Hessian returned an incompatible covector layout");
            }

          bool accepted = false;
          for (unsigned int trial = 0;
               trial < parameters_.maximum_trials_per_iteration;
               ++trial)
            {
              if (radius < parameters_.minimum_radius)
                return finish(ReducedTrustRegionStoppingReason::radius_too_small);

              const SubproblemResult subproblem =
                parameters_.subproblem_method ==
                    ReducedTrustRegionSubproblemMethod::cauchy
                  ? make_cauchy_subproblem(
                      current_evaluation.reduced_derivative,
                      metric_gradient,
                      *cauchy_hessian_gradient,
                      radius)
                  : make_truncated_cg_subproblem(
                      current_control,
                      current_evaluation.reduced_derivative,
                      metric_gradient,
                      radius);
              if (parameters_.subproblem_method !=
                  ReducedTrustRegionSubproblemMethod::cauchy)
                {
                  metric_solve_count += subproblem.metric_solve_count;
                  hessian_action_count += subproblem.hessian_action_count;
                }

              subproblem_status_history.push_back(subproblem.status);
              subproblem_iteration_history.push_back(
                subproblem.iteration_count);
              subproblem_residual_norm_history.push_back(
                subproblem.residual_norm);

              Primal trial_control = current_control;
              add_scaled_primal(trial_control, 1.0, subproblem.step);
              const auto trial_value = reduced_.evaluate_value(trial_control);
              ++state_solve_count;
              const double actual_reduction =
                current_evaluation.objective_value -
                trial_value.objective_value;
              const double reduction_ratio =
                actual_reduction / subproblem.predicted_reduction;
              const bool trial_accepted =
                std::isfinite(actual_reduction) &&
                std::isfinite(reduction_ratio) &&
                reduction_ratio >= parameters_.acceptance_threshold;

              radius_history.push_back(radius);
              step_norm_history.push_back(subproblem.step_norm);
              predicted_reduction_history.push_back(
                subproblem.predicted_reduction);
              actual_reduction_history.push_back(actual_reduction);
              reduction_ratio_history.push_back(reduction_ratio);
              accepted_history.push_back(trial_accepted);
              ++trial_count;

              if (trial_accepted)
                {
                  Evaluation trial_evaluation =
                    reduced_.augment_derivative(trial_value);
                  ++adjoint_solve_count;
                  Primal actual_update = subproblem.step;
                  const double actual_descent_pairing =
                    contract::pair(current_evaluation.reduced_derivative,
                                   actual_update);
                  contract::require(std::isfinite(actual_descent_pairing),
                                    "Trust-region accepted a non-finite descent pairing");
                  const ReducedIterationWork work_after{
                    state_solve_count,
                    adjoint_solve_count,
                    metric_solve_count,
                    hessian_action_count,
                    0};
                  current_control = std::move(trial_control);
                  current_evaluation = std::move(trial_evaluation);
                  objective_history.push_back(current_evaluation.objective_value);
                  ++accepted_iterations;
                  accepted = true;

                  iteration_records.push_back(
                    ReducedTrustRegionIterationRecordT<Backend>{
                      {accepted_iterations,
                       std::string("trust_region_") +
                         reduced_trust_region_subproblem_method_name(
                           parameters_.subproblem_method),
                       objective_history[objective_history.size() - 2],
                       current_evaluation.objective_value,
                       current_evaluation.objective_value -
                         objective_history[objective_history.size() - 2],
                       radius,
                       subproblem.step_norm,
                       actual_descent_pairing,
                       gradient_norm,
                       relative_gradient_norm,
                       trial + 1,
                       reduced_iteration_work_difference(work_after, work_before),
                       work_after,
                       {},
                       current_evaluation},
                      radius,
                      subproblem.predicted_reduction,
                      actual_reduction,
                      reduction_ratio,
                      subproblem.status,
                      subproblem.iteration_count,
                      subproblem.residual_norm});

                  if (reduction_ratio > parameters_.expansion_threshold &&
                      subproblem.step_norm >= 0.9 * radius)
                    radius = std::min(parameters_.maximum_radius,
                                      parameters_.expansion_factor * radius);
                  else if (reduction_ratio < parameters_.shrink_threshold)
                    radius = std::max(parameters_.minimum_radius,
                                      parameters_.shrink_factor * radius);

                  if ((automatic_stopping
                         ? parameters_.objective_change_tolerance > 0.0
                         : parameters_.stopping_criterion ==
                             ReducedStoppingCriterion::objective_change) &&
                      actual_reduction <=
                        parameters_.objective_change_tolerance)
                    return finish(
                      ReducedTrustRegionStoppingReason::objective_change_tolerance);

                  if ((automatic_stopping
                         ? parameters_.step_tolerance > 0.0
                     : parameters_.stopping_criterion ==
                             ReducedStoppingCriterion::step_norm) &&
                      subproblem.step_norm <= parameters_.step_tolerance)
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
    struct SubproblemResult
    {
      Primal                              step;
      double                              step_norm;
      double                              predicted_reduction;
      double                              residual_norm;
      std::size_t                         iteration_count;
      std::size_t                         metric_solve_count;
      std::size_t                         hessian_action_count;
      ReducedTrustRegionSubproblemStatus  status;
    };

    double
    metric_pairing(const Primal &left, const Primal &right) const
    {
      const Covector metric_left = metric_.apply(left);
      const double pairing = contract::pair(metric_left, right);
      contract::require(std::isfinite(pairing),
                        "Trust-region metric pairing is not finite");
      return pairing;
    }

    SubproblemResult
    make_cauchy_subproblem(const Covector &reduced_derivative,
                           const ReducedMetricGradientT<Backend> &metric_gradient,
                           const Covector &hessian_gradient,
                           const double radius) const
    {
      const double gradient_pairing =
        contract::pair(reduced_derivative, metric_gradient.gradient);
      const double curvature =
        contract::pair(hessian_gradient, metric_gradient.gradient);
      contract::require(std::isfinite(gradient_pairing) &&
                          gradient_pairing > 0.0,
                        "Trust-region gradient pairing is not positive");
      contract::require(std::isfinite(curvature) && curvature > 0.0,
                        "Trust-region model has non-positive curvature");

      const double unconstrained_step = gradient_pairing / curvature;
      const double step_scale = std::min(unconstrained_step,
                                         radius / metric_gradient.norm);
      const double step_norm = step_scale * metric_gradient.norm;
      const double predicted_reduction =
        step_scale * gradient_pairing -
        0.5 * step_scale * step_scale * curvature;
      contract::require(std::isfinite(step_scale) && step_scale > 0.0 &&
                          std::isfinite(step_norm) && step_norm > 0.0,
                        "Trust-region Cauchy step is invalid");
      contract::require(std::isfinite(predicted_reduction) &&
                          predicted_reduction > 0.0,
                        "Trust-region model predicted non-positive reduction");

      Primal step = metric_gradient.gradient;
      scale_primal(step, -step_scale);
      return {std::move(step),
              step_norm,
              predicted_reduction,
              0.0,
              0,
              0,
              0,
              ReducedTrustRegionSubproblemStatus::cauchy};
    }

    SubproblemResult
    make_boundary_subproblem(const Covector &reduced_derivative,
                             const Primal &current_step,
                             const Covector &current_hessian_step,
                             const Primal &search_direction,
                             const Covector &hessian_direction,
                             const double radius,
                             const std::size_t iteration_count,
                             const std::size_t metric_solve_count,
                             const std::size_t hessian_action_count,
                             const double residual_norm,
                             const ReducedTrustRegionSubproblemStatus status) const
    {
      const double current_norm_squared =
        metric_pairing(current_step, current_step);
      const double cross_pairing =
        metric_pairing(current_step, search_direction);
      const double direction_norm_squared =
        metric_pairing(search_direction, search_direction);
      contract::require(std::isfinite(current_norm_squared) &&
                          std::isfinite(cross_pairing) &&
                          std::isfinite(direction_norm_squared) &&
                          direction_norm_squared > 0.0,
                        "Trust-region boundary direction has invalid metric norm");

      const double radius_squared = radius * radius;
      const double raw_discriminant =
        cross_pairing * cross_pairing -
        direction_norm_squared * (current_norm_squared - radius_squared);
      contract::require(std::isfinite(raw_discriminant),
                        "Trust-region boundary intersection is not finite");
      const double discriminant = std::max(0.0, raw_discriminant);
      const double boundary_scale =
        (-cross_pairing + std::sqrt(discriminant)) /
        direction_norm_squared;
      contract::require(std::isfinite(boundary_scale) &&
                          boundary_scale >= 0.0,
                        "Trust-region boundary intersection is invalid");

      Primal step = current_step;
      add_scaled_primal(step, boundary_scale, search_direction);
      Covector hessian_step = current_hessian_step;
      add_scaled_covector(hessian_step,
                          boundary_scale,
                          hessian_direction);
      const double step_norm_squared = metric_pairing(step, step);
      const double step_norm = std::sqrt(std::max(0.0, step_norm_squared));
      const double predicted_reduction =
        -contract::pair(reduced_derivative, step) -
        0.5 * contract::pair(hessian_step, step);
      contract::require(std::isfinite(step_norm) && step_norm > 0.0,
                        "Trust-region boundary step norm is invalid");
      contract::require(std::abs(step_norm - radius) <=
                          1e-10 * std::max(1.0, radius),
                        "Trust-region boundary step does not fill the radius");
      contract::require(std::isfinite(predicted_reduction) &&
                          predicted_reduction > 0.0,
                        "Trust-region boundary step has non-positive predicted reduction");
      return {std::move(step),
              step_norm,
              predicted_reduction,
              residual_norm,
              iteration_count,
              metric_solve_count,
              hessian_action_count,
              status};
    }

    SubproblemResult
    make_truncated_cg_subproblem(
      const Primal &control,
      const Covector &reduced_derivative,
      const ReducedMetricGradientT<Backend> &metric_gradient,
      const double radius) const
    {
      Covector residual = reduced_derivative;
      Primal preconditioned_residual = metric_gradient.gradient;
      const double initial_pairing =
        contract::pair(residual, preconditioned_residual);
      contract::require(std::isfinite(initial_pairing) &&
                          initial_pairing > 0.0,
                        "Trust-region PCG residual pairing is not positive");
      const double initial_norm = std::sqrt(initial_pairing);
      const double target_norm = std::max(
        parameters_.subproblem_absolute_tolerance,
        parameters_.subproblem_relative_tolerance * initial_norm);
      Primal step = Primal::zeros(metric_.layout());
      Covector hessian_step = Covector::zeros(metric_.layout());

      // A nonzero outer residual must receive at least one inner action even
      // when the requested inner tolerance is already satisfied. This
      // forcing step prevents a loose inner absolute tolerance from creating
      // a zero predicted reduction and an unusable trust-region ratio.
      Primal search_direction = preconditioned_residual;
      scale_primal(search_direction, -1.0);
      double residual_pairing = initial_pairing;
      std::size_t metric_solve_count = 0;
      std::size_t hessian_action_count = 0;
      for (unsigned int iteration = 0;
           iteration < parameters_.maximum_subproblem_iterations;
           ++iteration)
        {
          const Covector hessian_direction =
            hessian_.apply(control, search_direction);
          ++hessian_action_count;
          contract::require(hessian_direction.layout()->compatible_with(
                              *metric_.layout()),
                            "Trust-region Hessian returned an incompatible covector layout");
          const double curvature =
            contract::pair(hessian_direction, search_direction);
          if (!std::isfinite(curvature) || curvature <= 0.0)
            return make_boundary_subproblem(
              reduced_derivative,
              step,
              hessian_step,
              search_direction,
              hessian_direction,
              radius,
              iteration + 1,
              metric_solve_count,
              hessian_action_count,
              std::sqrt(residual_pairing),
              ReducedTrustRegionSubproblemStatus::negative_curvature);

          const double step_length = residual_pairing / curvature;
          contract::require(std::isfinite(step_length) && step_length > 0.0,
                            "Trust-region PCG step is invalid");
          Primal candidate_step = step;
          add_scaled_primal(candidate_step,
                            step_length,
                            search_direction);
          const double candidate_norm_squared =
            metric_pairing(candidate_step, candidate_step);
          if (candidate_norm_squared >= radius * radius)
            return make_boundary_subproblem(
              reduced_derivative,
              step,
              hessian_step,
              search_direction,
              hessian_direction,
              radius,
              iteration + 1,
              metric_solve_count,
              hessian_action_count,
              std::sqrt(residual_pairing),
              ReducedTrustRegionSubproblemStatus::boundary);

          add_scaled_primal(step, step_length, search_direction);
          add_scaled_covector(hessian_step,
                              step_length,
                              hessian_direction);
          add_scaled_covector(residual, step_length, hessian_direction);
          const Primal next_preconditioned_residual =
            metric_.inverse_apply(residual);
          ++metric_solve_count;
          const double next_pairing =
            contract::pair(residual, next_preconditioned_residual);
          contract::require(std::isfinite(next_pairing) &&
                              next_pairing >= 0.0,
                            "Trust-region PCG residual norm is invalid");
          const double next_norm = std::sqrt(next_pairing);
          const std::size_t iteration_count = iteration + 1;
          if (next_norm <= target_norm)
            {
              const double step_norm =
                std::sqrt(std::max(0.0, metric_pairing(step, step)));
              const double predicted_reduction =
                -contract::pair(reduced_derivative, step) -
                0.5 * contract::pair(hessian_step, step);
              contract::require(std::isfinite(step_norm) && step_norm > 0.0 &&
                                  std::isfinite(predicted_reduction) &&
                                  predicted_reduction > 0.0,
                                "Trust-region PCG converged step is invalid");
              return {std::move(step),
                      step_norm,
                      predicted_reduction,
                      next_norm,
                      iteration_count,
                      metric_solve_count,
                      hessian_action_count,
                      ReducedTrustRegionSubproblemStatus::converged};
            }

          if (iteration_count == parameters_.maximum_subproblem_iterations)
            {
              const double step_norm =
                std::sqrt(std::max(0.0, metric_pairing(step, step)));
              const double predicted_reduction =
                -contract::pair(reduced_derivative, step) -
                0.5 * contract::pair(hessian_step, step);
              contract::require(std::isfinite(step_norm) && step_norm > 0.0 &&
                                  std::isfinite(predicted_reduction) &&
                                  predicted_reduction > 0.0,
                                "Trust-region PCG iteration-limit step is invalid");
              return {std::move(step),
                      step_norm,
                      predicted_reduction,
                      next_norm,
                      iteration_count,
                      metric_solve_count,
                      hessian_action_count,
                      ReducedTrustRegionSubproblemStatus::iteration_limit};
            }

          const double beta = next_pairing / residual_pairing;
          contract::require(std::isfinite(beta) && beta >= 0.0,
                            "Trust-region PCG coefficient is invalid");
          scale_primal(search_direction, beta);
          add_scaled_primal(search_direction,
                            -1.0,
                            next_preconditioned_residual);
          residual_pairing = next_pairing;
        }

      contract::require(false, "Trust-region PCG did not return a subproblem step");
      return {Primal::zeros(metric_.layout()),
              0.0,
              0.0,
              0.0,
              0,
              0,
              0,
              ReducedTrustRegionSubproblemStatus::iteration_limit};
    }

    void
    validate_parameters() const
    {
      contract::require(parameters_.maximum_iterations > 0,
                        "Trust-region maximum iterations must be positive");
      contract::require(parameters_.maximum_trials_per_iteration > 0,
                        "Trust-region trial limit must be positive");
      contract::require(parameters_.gradient_tolerance > 0.0,
                        "Trust-region gradient tolerance must be positive");
      if (parameters_.stopping_criterion ==
          ReducedStoppingCriterion::relative_gradient_norm)
        contract::require(parameters_.relative_gradient_tolerance > 0.0,
                          "Selected relative-gradient tolerance must be positive");
      if (parameters_.stopping_criterion ==
          ReducedStoppingCriterion::objective_change)
        contract::require(parameters_.objective_change_tolerance > 0.0,
                          "Selected objective-change tolerance must be positive");
      if (parameters_.stopping_criterion ==
          ReducedStoppingCriterion::step_norm)
        contract::require(parameters_.step_tolerance > 0.0,
                          "Selected step tolerance must be positive");
      contract::require(parameters_.relative_gradient_tolerance >= 0.0,
                        "Trust-region relative gradient tolerance must be nonnegative");
      contract::require(parameters_.objective_change_tolerance >= 0.0,
                        "Trust-region objective-change tolerance must be nonnegative");
      contract::require(parameters_.step_tolerance >= 0.0,
                        "Trust-region step tolerance must be nonnegative");
      contract::require(parameters_.maximum_subproblem_iterations > 0,
                        "Trust-region subproblem iteration limit must be positive");
      contract::require(parameters_.subproblem_relative_tolerance > 0.0,
                        "Trust-region subproblem relative tolerance must be positive");
      contract::require(parameters_.subproblem_absolute_tolerance > 0.0,
                        "Trust-region subproblem absolute tolerance must be positive");
      contract::require(
        parameters_.subproblem_method ==
            ReducedTrustRegionSubproblemMethod::cauchy ||
          parameters_.subproblem_method ==
            ReducedTrustRegionSubproblemMethod::truncated_conjugate_gradient,
        "Trust-region subproblem method is unknown");
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

    Result
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
           std::vector<ReducedTrustRegionSubproblemStatus>
             subproblem_status_history,
           std::vector<std::size_t>       subproblem_iteration_history,
           std::vector<double>            subproblem_residual_norm_history,
           const std::size_t              accepted_iterations,
           const std::size_t              trial_count,
           const std::size_t              state_solve_count,
           const std::size_t              adjoint_solve_count,
           const std::size_t              metric_solve_count,
           const std::size_t              hessian_action_count,
           const ReducedTrustRegionStoppingReason stopping_reason,
           std::vector<ReducedTrustRegionIterationRecordT<Backend>>
             iteration_records) const
    {
      contract::require(iteration_records.size() == accepted_iterations,
                        "Trust-region iteration records are incomplete");
      for (std::size_t index = 0; index < accepted_iterations; ++index)
        {
          const auto &audit_record = iteration_records[index];
          const auto &record = audit_record.common;
          contract::require(record.iteration == index + 1 &&
                              record.objective_before ==
                                objective_history[index] &&
                              record.objective_after ==
                                objective_history[index + 1] &&
                              record.objective_change ==
                                record.objective_after -
                                  record.objective_before &&
                              record.absolute_stationarity ==
                                gradient_norm_history[index] &&
                              record.relative_stationarity ==
                                relative_gradient_norm_history[index] &&
                              record.accepted_evaluation.objective_value ==
                                objective_history[index + 1] &&
                              audit_record.predicted_reduction > 0.0 &&
                              audit_record.actual_reduction > 0.0,
                            "Trust-region compatibility histories disagree with audit records");
        }
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
              std::move(subproblem_status_history),
              std::move(subproblem_iteration_history),
              std::move(subproblem_residual_norm_history),
              accepted_iterations,
              trial_count,
              state_solve_count,
              adjoint_solve_count,
              metric_solve_count,
              hessian_action_count,
              stopping_reason,
              parameters_,
              std::string("trust_region_") +
                reduced_trust_region_subproblem_method_name(
                  parameters_.subproblem_method),
              std::move(iteration_records)};
    }

    const contract::ReducedDTOT<Backend> &   reduced_;
    const contract::MetricT<Backend> &       metric_;
    const contract::ReducedHessianT<Backend> &hessian_;
    ReducedTrustRegionParameters             parameters_;
  };

  using ReducedTrustRegionSolver =
    ReducedTrustRegionSolverT<contract::DenseBackend>;
} // namespace nmopt::solvers
