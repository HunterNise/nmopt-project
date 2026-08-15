#pragma once

#include "nmopt/contract/reduced_dto.hpp"
#include "nmopt/solvers/reduced_search.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <string>
#include <utility>

namespace nmopt::solvers
{
  enum class ReducedLineSearchStatus
  {
    accepted,
    failure
  };

  inline const char *
  reduced_line_search_status_name(const ReducedLineSearchStatus status)
  {
    switch (status)
      {
        case ReducedLineSearchStatus::accepted:
          return "accepted";
        case ReducedLineSearchStatus::failure:
          return "failure";
      }
    return "unknown";
  }

  template <typename Backend>
  using ReducedTrialControlBuilderT =
    std::function<contract::PrimalBlockT<Backend>(double)>;

  template <typename Backend>
  using ReducedTrialEvaluatorT = std::function<
    contract::ReducedEvaluationT<Backend>(
      const contract::PrimalBlockT<Backend> &)>;

  template <typename Backend>
  struct ReducedLineSearchResultT
  {
    ReducedLineSearchStatus             status;
    contract::PrimalBlockT<Backend>     control;
    contract::ReducedEvaluationT<Backend> evaluation;
    double                              step_length;
    std::size_t                         trial_count;
    std::size_t                         hessian_action_count;

    bool
    accepted() const noexcept
    {
      return status == ReducedLineSearchStatus::accepted;
    }
  };

  using ReducedLineSearchResult =
    ReducedLineSearchResultT<contract::DenseBackend>;

  namespace detail
  {
    template <typename Backend>
    inline void
    validate_line_search_inputs(
      const contract::PrimalBlockT<Backend> &       current_control,
      const contract::ReducedEvaluationT<Backend> & current_evaluation,
      const ReducedSearchDirectionT<Backend> &      direction,
      const ReducedTrialControlBuilderT<Backend> &  build_trial_control,
      const ReducedTrialEvaluatorT<Backend> &       evaluate_trial,
      const char *                                  policy_name)
    {
      contract::require(current_control.layout()->compatible_with(
                          *direction.direction.layout()),
                        std::string(policy_name) +
                          " direction has an incompatible control layout");
      contract::require(current_evaluation.reduced_derivative.layout()->compatible_with(
                          *current_control.layout()),
                        std::string(policy_name) +
                          " derivative has an incompatible control layout");
      contract::require(static_cast<bool>(build_trial_control),
                        std::string(policy_name) +
                          " requires a trial-control builder");
      contract::require(static_cast<bool>(evaluate_trial),
                        std::string(policy_name) +
                          " requires a trial evaluator");
      contract::require(std::isfinite(direction.directional_derivative) &&
                          direction.directional_derivative < 0.0,
                        std::string(policy_name) +
                          " requires a descent direction");
    }

    template <typename Backend>
    inline ReducedLineSearchResultT<Backend>
    failure(const contract::PrimalBlockT<Backend> &       current_control,
            const contract::ReducedEvaluationT<Backend> & current_evaluation,
            const std::size_t                             trial_count,
            const std::size_t                             hessian_action_count = 0)
    {
      return {ReducedLineSearchStatus::failure,
              current_control,
              current_evaluation,
              0.0,
              trial_count,
              hessian_action_count};
    }

    template <typename Backend>
    inline ReducedLineSearchResultT<Backend>
    accepted(const contract::PrimalBlockT<Backend> &       control,
             contract::ReducedEvaluationT<Backend>        evaluation,
             const double                                 step_length,
             const std::size_t                             trial_count,
             const std::size_t                             hessian_action_count = 0)
    {
      return {ReducedLineSearchStatus::accepted,
              control,
              std::move(evaluation),
              step_length,
              trial_count,
              hessian_action_count};
    }
  } // namespace detail

  struct ArmijoLineSearchParameters
  {
    unsigned int maximum_trials = 20;
    double       initial_step_length = 1.0;
    double       armijo_fraction = 1e-4;
    double       backtracking_factor = 0.5;
  };

  template <typename Backend>
  class ArmijoLineSearchPolicyT
  {
  public:
    using Result = ReducedLineSearchResultT<Backend>;
    using Primal = contract::PrimalBlockT<Backend>;
    using Covector = contract::CovectorBlockT<Backend>;
    using Evaluation = contract::ReducedEvaluationT<Backend>;

    ArmijoLineSearchPolicyT(ArmijoLineSearchParameters parameters = {})
      : parameters_(parameters)
    {
      validate_parameters();
    }

    Result
    search(const Primal &       current_control,
           const Evaluation &   current_evaluation,
           const ReducedSearchDirectionT<Backend> &direction,
           const ReducedTrialControlBuilderT<Backend> &build_trial_control,
           const ReducedTrialEvaluatorT<Backend> &     evaluate_trial) const
    {
      detail::validate_line_search_inputs(current_control,
                                          current_evaluation,
                                          direction,
                                          build_trial_control,
                                          evaluate_trial,
                                          "Armijo line search");

      double step_length = parameters_.initial_step_length;
      for (unsigned int trial = 0; trial < parameters_.maximum_trials; ++trial)
        {
          Primal trial_control = build_trial_control(step_length);
          contract::require(trial_control.layout()->compatible_with(
                              *current_control.layout()),
                            "Armijo trial control has an incompatible layout");
          Evaluation trial_evaluation = evaluate_trial(trial_control);
          contract::require(trial_evaluation.reduced_derivative.layout()->compatible_with(
                              *current_control.layout()),
                            "Armijo trial derivative has an incompatible layout");

          Primal actual_update = trial_control;
          add_scaled_primal(actual_update, -1.0, current_control);
          const double actual_slope =
            contract::pair(current_evaluation.reduced_derivative, actual_update);
          contract::require(std::isfinite(actual_slope),
                            "Armijo trial produced a non-finite actual slope");
          const double armijo_bound =
            current_evaluation.objective_value +
            parameters_.armijo_fraction * actual_slope;
          if (actual_slope < 0.0 &&
              std::isfinite(trial_evaluation.objective_value) &&
              trial_evaluation.objective_value <= armijo_bound)
            return detail::accepted(trial_control,
                                    std::move(trial_evaluation),
                                    step_length,
                                    trial + 1);
          step_length *= parameters_.backtracking_factor;
        }

      return detail::failure(current_control,
                             current_evaluation,
                             parameters_.maximum_trials);
    }

  private:
    void
    validate_parameters() const
    {
      contract::require(parameters_.maximum_trials > 0,
                        "Armijo line-search trial limit must be positive");
      contract::require(parameters_.initial_step_length > 0.0,
                        "Armijo line-search initial step must be positive");
      contract::require(parameters_.armijo_fraction > 0.0 &&
                          parameters_.armijo_fraction < 1.0,
                        "Armijo line-search fraction must lie in (0, 1)");
      contract::require(parameters_.backtracking_factor > 0.0 &&
                          parameters_.backtracking_factor < 1.0,
                        "Armijo line-search backtracking factor must lie in (0, 1)");
    }

    ArmijoLineSearchParameters parameters_;
  };

  using ArmijoLineSearchPolicy =
    ArmijoLineSearchPolicyT<contract::DenseBackend>;

  struct ExactQuadraticLineSearchParameters
  {
    double curvature_tolerance = 1e-14;
    double objective_tolerance = 1e-12;
  };

  template <typename Backend>
  class ExactQuadraticLineSearchPolicyT
  {
  public:
    using Result = ReducedLineSearchResultT<Backend>;
    using Primal = contract::PrimalBlockT<Backend>;
    using Covector = contract::CovectorBlockT<Backend>;
    using Evaluation = contract::ReducedEvaluationT<Backend>;

    ExactQuadraticLineSearchPolicyT() = default;

    ExactQuadraticLineSearchPolicyT(
      const contract::ReducedHessianT<Backend> &hessian,
      ExactQuadraticLineSearchParameters         parameters = {})
      : hessian_(&hessian)
      , parameters_(parameters)
    {
      validate_parameters();
    }

    Result
    search(const Primal &       current_control,
           const Evaluation &   current_evaluation,
           const ReducedSearchDirectionT<Backend> &direction,
           const ReducedTrialControlBuilderT<Backend> &build_trial_control,
           const ReducedTrialEvaluatorT<Backend> &     evaluate_trial) const
    {
      detail::validate_line_search_inputs(current_control,
                                          current_evaluation,
                                          direction,
                                          build_trial_control,
                                          evaluate_trial,
                                          "Exact quadratic line search");
      contract::require(hessian_ != nullptr,
                        "Exact quadratic line search requires a reduced Hessian capability");
      contract::require(current_control.layout()->compatible_with(
                          *hessian_->layout()),
                        "Exact quadratic line search control does not match Hessian");

      const Covector hessian_direction =
        hessian_->apply(current_control, direction.direction);
      const double curvature =
        contract::pair(hessian_direction, direction.direction);
      contract::require(std::isfinite(curvature) &&
                          curvature > parameters_.curvature_tolerance,
                        "Exact quadratic line search requires positive curvature");
      const double step_length =
        -direction.directional_derivative / curvature;
      contract::require(std::isfinite(step_length) && step_length > 0.0,
                        "Exact quadratic line search returned an invalid step");

      Primal trial_control = build_trial_control(step_length);
      contract::require(trial_control.layout()->compatible_with(
                          *current_control.layout()),
                        "Exact quadratic trial control has an incompatible layout");
      Evaluation trial_evaluation = evaluate_trial(trial_control);
      contract::require(trial_evaluation.reduced_derivative.layout()->compatible_with(
                          *current_control.layout()),
                        "Exact quadratic trial derivative has an incompatible layout");
      Primal actual_update = trial_control;
      add_scaled_primal(actual_update, -1.0, current_control);
      const double actual_slope =
        contract::pair(current_evaluation.reduced_derivative, actual_update);
      const double objective_bound =
        current_evaluation.objective_value +
        parameters_.objective_tolerance *
          std::max(1.0, std::abs(current_evaluation.objective_value));
      if (std::isfinite(actual_slope) && actual_slope < 0.0 &&
          std::isfinite(trial_evaluation.objective_value) &&
          trial_evaluation.objective_value <= objective_bound)
        return detail::accepted(trial_control,
                                std::move(trial_evaluation),
                                step_length,
                                1,
                                1);
      return detail::failure(current_control, current_evaluation, 1, 1);
    }

  private:
    void
    validate_parameters() const
    {
      contract::require(parameters_.curvature_tolerance > 0.0,
                        "Exact quadratic curvature tolerance must be positive");
      contract::require(parameters_.objective_tolerance >= 0.0,
                        "Exact quadratic objective tolerance must be nonnegative");
    }

    const contract::ReducedHessianT<Backend> *hessian_ = nullptr;
    ExactQuadraticLineSearchParameters         parameters_;
  };

  using ExactQuadraticLineSearchPolicy =
    ExactQuadraticLineSearchPolicyT<contract::DenseBackend>;

  struct WolfeLineSearchParameters
  {
    unsigned int maximum_trials = 20;
    double       initial_step_length = 1.0;
    double       backtracking_factor = 0.5;
    double       sufficient_decrease_fraction = 1e-4;
    double       curvature_fraction = 0.9;
  };

  template <typename Backend>
  class WolfeLineSearchPolicyT
  {
  public:
    using Result = ReducedLineSearchResultT<Backend>;
    using Primal = contract::PrimalBlockT<Backend>;
    using Evaluation = contract::ReducedEvaluationT<Backend>;

    WolfeLineSearchPolicyT(WolfeLineSearchParameters parameters = {})
      : parameters_(parameters)
    {
      validate_parameters();
    }

    Result
    search(const Primal &       current_control,
           const Evaluation &   current_evaluation,
           const ReducedSearchDirectionT<Backend> &direction,
           const ReducedTrialControlBuilderT<Backend> &build_trial_control,
           const ReducedTrialEvaluatorT<Backend> &     evaluate_trial) const
    {
      detail::validate_line_search_inputs(current_control,
                                          current_evaluation,
                                          direction,
                                          build_trial_control,
                                          evaluate_trial,
                                          "Wolfe line search");

      double step_length = parameters_.initial_step_length;
      for (unsigned int trial = 0; trial < parameters_.maximum_trials; ++trial)
        {
          Primal trial_control = build_trial_control(step_length);
          contract::require(trial_control.layout()->compatible_with(
                              *current_control.layout()),
                            "Wolfe trial control has an incompatible layout");
          Evaluation trial_evaluation = evaluate_trial(trial_control);
          contract::require(trial_evaluation.reduced_derivative.layout()->compatible_with(
                              *current_control.layout()),
                            "Wolfe trial derivative has an incompatible layout");

          Primal actual_update = trial_control;
          add_scaled_primal(actual_update, -1.0, current_control);
          const double actual_initial_slope =
            contract::pair(current_evaluation.reduced_derivative, actual_update);
          const double actual_trial_slope =
            contract::pair(trial_evaluation.reduced_derivative, actual_update);
          contract::require(std::isfinite(actual_initial_slope) &&
                              std::isfinite(actual_trial_slope),
                            "Wolfe trial produced a non-finite slope");
          const double sufficient_decrease_bound =
            current_evaluation.objective_value +
            parameters_.sufficient_decrease_fraction * actual_initial_slope;
          if (actual_initial_slope < 0.0 &&
              std::isfinite(trial_evaluation.objective_value) &&
              trial_evaluation.objective_value <= sufficient_decrease_bound &&
              std::abs(actual_trial_slope) <=
                parameters_.curvature_fraction * std::abs(actual_initial_slope))
            return detail::accepted(trial_control,
                                    std::move(trial_evaluation),
                                    step_length,
                                    trial + 1);
          step_length *= parameters_.backtracking_factor;
        }

      return detail::failure(current_control,
                             current_evaluation,
                             parameters_.maximum_trials);
    }

  private:
    void
    validate_parameters() const
    {
      contract::require(parameters_.maximum_trials > 0,
                        "Wolfe line-search trial limit must be positive");
      contract::require(parameters_.initial_step_length > 0.0,
                        "Wolfe line-search initial step must be positive");
      contract::require(parameters_.backtracking_factor > 0.0 &&
                          parameters_.backtracking_factor < 1.0,
                        "Wolfe line-search backtracking factor must lie in (0, 1)");
      contract::require(parameters_.sufficient_decrease_fraction > 0.0 &&
                          parameters_.sufficient_decrease_fraction < 1.0,
                        "Wolfe sufficient-decrease fraction must lie in (0, 1)");
      contract::require(parameters_.curvature_fraction >
                          parameters_.sufficient_decrease_fraction &&
                          parameters_.curvature_fraction < 1.0,
                        "Wolfe curvature fraction must lie between sufficient decrease and 1");
    }

    WolfeLineSearchParameters parameters_;
  };

  using WolfeLineSearchPolicy =
    WolfeLineSearchPolicyT<contract::DenseBackend>;
} // namespace nmopt::solvers
