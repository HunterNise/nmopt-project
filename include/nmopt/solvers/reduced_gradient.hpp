#pragma once

#include "nmopt/contract/reduced_dto.hpp"
#include "nmopt/solvers/reduced_line_search.hpp"
#include "nmopt/solvers/reduced_search.hpp"

#include <cmath>
#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>

namespace nmopt::solvers
{
  // An unconstrained reduced-space method. The metric identifies the reduced
  // covector j' with g = G^{-1} j'; the direction protocol supplies d = -g.
  // Armijo accepts u + alpha d when j(u + alpha d) <=
  // j(u) + c <j', alpha d>.
  template <typename Backend,
            typename DirectionPolicy =
              SteepestDescentDirectionPolicyT<Backend>,
            typename LineSearchPolicy = ArmijoLineSearchPolicyT<Backend>>
  class ReducedSearchSolverT
  {
  public:
    using Primal = contract::PrimalBlockT<Backend>;
    using Covector = contract::CovectorBlockT<Backend>;
    using Evaluation = contract::ReducedEvaluationT<Backend>;
    using Result = ReducedGradientResultT<Backend>;

    ReducedSearchSolverT(const contract::ReducedDTOT<Backend> &reduced,
                         const contract::MetricT<Backend> &    metric,
                         ReducedSolverParameters               parameters = {},
                         DirectionPolicy                        direction_policy = {})
      : reduced_(reduced)
      , metric_(metric)
      , parameters_(parameters)
      , direction_policy_(std::move(direction_policy))
      , line_search_policy_(default_line_search_policy(parameters))
    {
      validate_parameters();
    }

    ReducedSearchSolverT(const contract::ReducedDTOT<Backend> &reduced,
                         const contract::MetricT<Backend> &    metric,
                         ReducedSolverParameters               parameters,
                         DirectionPolicy                        direction_policy,
                         LineSearchPolicy                        line_search_policy)
      : reduced_(reduced)
      , metric_(metric)
      , parameters_(parameters)
      , direction_policy_(std::move(direction_policy))
      , line_search_policy_(std::move(line_search_policy))
    {
      validate_parameters();
    }

    ReducedSearchSolverT(const contract::ReducedDTOT<Backend> &reduced,
                         const contract::MetricT<Backend> &    metric,
                         const contract::ConstraintT<Backend> & constraint,
                         ReducedSolverParameters               parameters = {},
                         DirectionPolicy                        direction_policy = {})
      : reduced_(reduced)
      , metric_(metric)
      , constraint_(&constraint)
      , parameters_(parameters)
      , direction_policy_(std::move(direction_policy))
      , line_search_policy_(default_line_search_policy(parameters))
    {
      validate_parameters();
      contract::require(constraint_->layout()->compatible_with(*metric_.layout()),
                        "Reduced gradient constraint does not match metric");
      contract::require(constraint_->supports_projection_in(metric_),
                        "Reduced gradient constraint cannot project in this metric");
    }

    ReducedSearchSolverT(const contract::ReducedDTOT<Backend> &reduced,
                         const contract::MetricT<Backend> &    metric,
                         const contract::ConstraintT<Backend> & constraint,
                         ReducedSolverParameters               parameters,
                         DirectionPolicy                        direction_policy,
                         LineSearchPolicy                        line_search_policy)
      : reduced_(reduced)
      , metric_(metric)
      , constraint_(&constraint)
      , parameters_(parameters)
      , direction_policy_(std::move(direction_policy))
      , line_search_policy_(std::move(line_search_policy))
    {
      validate_parameters();
      contract::require(constraint_->layout()->compatible_with(*metric_.layout()),
                        "Reduced gradient constraint does not match metric");
      contract::require(constraint_->supports_projection_in(metric_),
                        "Reduced gradient constraint cannot project in this metric");
    }

    Result
    solve(const Primal &initial_control) const
    {
      contract::require(initial_control.layout()->compatible_with(*metric_.layout()),
                        "Reduced gradient initial control does not match metric");
      if (constraint_ != nullptr)
        contract::require(constraint_->is_feasible(initial_control),
                          "Projected reduced gradient requires a feasible initial control");

      std::size_t state_solve_count = 0;
      std::size_t adjoint_solve_count = 0;
      std::size_t line_search_trial_count = 0;
      std::size_t accepted_iterations = 0;
      std::size_t metric_solve_count = 0;
      std::size_t hessian_action_count = 0;
      double       initial_gradient_norm = 0.0;
      bool         have_initial_gradient_norm = false;

      Primal current_control = initial_control;
      Evaluation current_evaluation = reduced_.evaluate(current_control);
      ++state_solve_count;
      ++adjoint_solve_count;

      std::vector<double> objective_history{current_evaluation.objective_value};
      std::vector<double> gradient_norm_history;
      std::vector<double> relative_gradient_norm_history;
      std::vector<ReducedHessianSolveDiagnostics> hessian_solve_history;
      std::vector<double> step_length_history;
      std::vector<double> step_norm_history;
      std::vector<double> objective_change_history;
      DirectionPolicy direction_policy = direction_policy_;
      direction_policy.reset();
      LineSearchPolicy line_search_policy = line_search_policy_;

      for (;;)
        {
          const ReducedSearchDirectionT<Backend> direction =
            direction_policy.next(current_control,
                                  current_evaluation.reduced_derivative,
                                  metric_);
          metric_solve_count += direction.metric_solve_count;
          hessian_action_count += direction.hessian_action_count;
          hessian_solve_history.push_back(direction.hessian_solve);
          double stopping_norm = direction.gradient_norm;
          double unit_step_descent_measure = direction.directional_derivative;
          if (constraint_ != nullptr)
            {
              const ProjectedGradientData projected_gradient =
                evaluate_projected_gradient(current_control,
                                            current_evaluation.reduced_derivative,
                                            direction.direction);
              stopping_norm = projected_gradient.metric_norm;
              unit_step_descent_measure = projected_gradient.descent_measure;
            }
          if (!have_initial_gradient_norm)
            {
              initial_gradient_norm = stopping_norm;
              have_initial_gradient_norm = true;
            }
          const double relative_gradient_norm =
            initial_gradient_norm > 0.0
              ? stopping_norm / initial_gradient_norm
              : 0.0;
          gradient_norm_history.push_back(stopping_norm);
          relative_gradient_norm_history.push_back(relative_gradient_norm);

          const bool automatic_stopping =
            parameters_.stopping_criterion ==
            ReducedStoppingCriterion::automatic;
          if ((automatic_stopping ||
               parameters_.stopping_criterion ==
                 ReducedStoppingCriterion::gradient_norm) &&
              stopping_norm <= parameters_.gradient_tolerance)
            return result(current_control,
                          std::move(current_evaluation),
                          std::move(objective_history),
                          std::move(gradient_norm_history),
                          std::move(relative_gradient_norm_history),
                          std::move(hessian_solve_history),
                          accepted_iterations,
                          line_search_trial_count,
                          state_solve_count,
                          adjoint_solve_count,
                          ReducedGradientStoppingReason::gradient_tolerance,
                          std::move(step_length_history),
                          std::move(step_norm_history),
                          std::move(objective_change_history),
                          metric_solve_count,
                          hessian_action_count,
                          direction_policy.reset_count());

          if ((automatic_stopping
                 ? parameters_.relative_gradient_tolerance > 0.0
                 : parameters_.stopping_criterion ==
                     ReducedStoppingCriterion::relative_gradient_norm) &&
              relative_gradient_norm <= parameters_.relative_gradient_tolerance)
            return result(current_control,
                          std::move(current_evaluation),
                          std::move(objective_history),
                          std::move(gradient_norm_history),
                          std::move(relative_gradient_norm_history),
                          std::move(hessian_solve_history),
                          accepted_iterations,
                          line_search_trial_count,
                          state_solve_count,
                          adjoint_solve_count,
                          ReducedGradientStoppingReason::relative_gradient_tolerance,
                          std::move(step_length_history),
                          std::move(step_norm_history),
                          std::move(objective_change_history),
                          metric_solve_count,
                          hessian_action_count,
                          direction_policy.reset_count());

          contract::require(unit_step_descent_measure < 0.0,
                            "Reduced gradient did not produce a descent direction");

          if (accepted_iterations == parameters_.maximum_iterations)
            return result(current_control,
                          std::move(current_evaluation),
                          std::move(objective_history),
                          std::move(gradient_norm_history),
                          std::move(relative_gradient_norm_history),
                          std::move(hessian_solve_history),
                          accepted_iterations,
                          line_search_trial_count,
                          state_solve_count,
                          adjoint_solve_count,
                          ReducedGradientStoppingReason::maximum_iterations,
                          std::move(step_length_history),
                          std::move(step_norm_history),
                          std::move(objective_change_history),
                          metric_solve_count,
                          hessian_action_count,
                          direction_policy.reset_count());

          const auto build_trial_control = [this,
                                            &current_control,
                                            &direction](const double step_length) {
            Primal trial_control = current_control;
            add_scaled_primal(trial_control,
                              step_length,
                              direction.direction);
            if (constraint_ != nullptr)
              {
                trial_control = constraint_->project_in(trial_control, metric_);
                contract::require(constraint_->is_feasible(trial_control),
                                  "Reduced gradient projection returned an infeasible control");
              }
            return trial_control;
          };
          const auto evaluate_trial_value = [this,
                                             &state_solve_count](
                                               const Primal &trial_control) {
            auto trial_value = reduced_.evaluate_value(trial_control);
            ++state_solve_count;
            return trial_value;
          };
          const auto augment_trial_derivative = [this,
                                                 &adjoint_solve_count](
                                                   const auto &trial_value) {
            auto trial_evaluation = reduced_.augment_derivative(trial_value);
            ++adjoint_solve_count;
            return trial_evaluation;
          };
          auto line_search_result = line_search_policy.search(
            current_control,
            current_evaluation,
            direction,
            build_trial_control,
            evaluate_trial_value,
            augment_trial_derivative);
          line_search_trial_count += line_search_result.trial_count;
          hessian_action_count += line_search_result.hessian_action_count;

          if (!line_search_result.accepted())
            return result(current_control,
                          std::move(current_evaluation),
                          std::move(objective_history),
                          std::move(gradient_norm_history),
                          std::move(relative_gradient_norm_history),
                          std::move(hessian_solve_history),
                          accepted_iterations,
                          line_search_trial_count,
                          state_solve_count,
                          adjoint_solve_count,
                          ReducedGradientStoppingReason::line_search_failure,
                          std::move(step_length_history),
                          std::move(step_norm_history),
                          std::move(objective_change_history),
                          metric_solve_count,
                          hessian_action_count,
                          direction_policy.reset_count());

          Primal trial_control = std::move(line_search_result.control);
          Evaluation trial_evaluation =
            std::move(line_search_result.evaluation);
          const double objective_change =
            trial_evaluation.objective_value - current_evaluation.objective_value;
          Primal actual_update = trial_control;
          add_scaled_primal(actual_update, -1.0, current_control);
          const double step_norm = metric_norm(actual_update);
          current_control = std::move(trial_control);
          current_evaluation = std::move(trial_evaluation);
          objective_history.push_back(current_evaluation.objective_value);
          step_length_history.push_back(line_search_result.step_length);
          step_norm_history.push_back(step_norm);
          objective_change_history.push_back(objective_change);
          ++accepted_iterations;

          if ((automatic_stopping
                 ? parameters_.objective_change_tolerance > 0.0
                 : parameters_.stopping_criterion ==
                     ReducedStoppingCriterion::objective_change) &&
              std::abs(objective_change) <=
                parameters_.objective_change_tolerance)
            return result(current_control,
                          std::move(current_evaluation),
                          std::move(objective_history),
                          std::move(gradient_norm_history),
                          std::move(relative_gradient_norm_history),
                          std::move(hessian_solve_history),
                          accepted_iterations,
                          line_search_trial_count,
                          state_solve_count,
                          adjoint_solve_count,
                          ReducedGradientStoppingReason::objective_change_tolerance,
                          std::move(step_length_history),
                          std::move(step_norm_history),
                          std::move(objective_change_history),
                          metric_solve_count,
                          hessian_action_count,
                          direction_policy.reset_count());

          if ((automatic_stopping
                 ? parameters_.step_tolerance > 0.0
                 : parameters_.stopping_criterion ==
                     ReducedStoppingCriterion::step_norm) &&
              step_norm <= parameters_.step_tolerance)
            return result(current_control,
                          std::move(current_evaluation),
                          std::move(objective_history),
                          std::move(gradient_norm_history),
                          std::move(relative_gradient_norm_history),
                          std::move(hessian_solve_history),
                          accepted_iterations,
                          line_search_trial_count,
                          state_solve_count,
                          adjoint_solve_count,
                          ReducedGradientStoppingReason::step_tolerance,
                          std::move(step_length_history),
                          std::move(step_norm_history),
                          std::move(objective_change_history),
                          metric_solve_count,
                          hessian_action_count,
                          direction_policy.reset_count());
        }
    }

  private:
    static LineSearchPolicy
    default_line_search_policy(const ReducedSolverParameters &parameters)
    {
      if constexpr (std::is_same_v<LineSearchPolicy,
                                   ArmijoLineSearchPolicyT<Backend>>)
        return LineSearchPolicy(ArmijoLineSearchParameters{
          parameters.maximum_line_search_trials,
          parameters.initial_step_length,
          parameters.armijo_fraction,
          parameters.backtracking_factor});
      else
        return LineSearchPolicy{};
    }

    struct ProjectedGradientData
    {
      double metric_norm;
      double descent_measure;
    };

    double
    metric_norm(const Primal &value) const
    {
      const Covector metric_value = metric_.apply(value);
      const double squared_norm = contract::pair(metric_value, value);
      contract::require(std::isfinite(squared_norm) && squared_norm >= 0.0,
                        "Reduced solver metric returned an invalid step norm");
      return std::sqrt(squared_norm);
    }

    ProjectedGradientData
    evaluate_projected_gradient(const Primal &  control,
                                const Covector &reduced_derivative,
                                const Primal &  direction) const
    {
      Primal projected_control = control;
      add_scaled_primal(projected_control, 1.0, direction);
      projected_control = constraint_->project_in(projected_control, metric_);
      contract::require(constraint_->is_feasible(projected_control),
                        "Reduced gradient projection returned an infeasible control");

      Primal update = projected_control;
      add_scaled_primal(update, -1.0, control);
      const Covector metric_update = metric_.apply(update);
      const double metric_norm_squared = contract::pair(metric_update, update);
      contract::require(metric_norm_squared >= 0.0,
                        "Projected gradient metric returned a negative squared norm");
      return {std::sqrt(metric_norm_squared),
              contract::pair(reduced_derivative, update)};
    }

    void
    validate_parameters() const
    {
      contract::require(parameters_.maximum_iterations > 0,
                        "Reduced gradient maximum iterations must be positive");
      contract::require(parameters_.maximum_line_search_trials > 0,
                        "Reduced gradient line-search trials must be positive");
      contract::require(parameters_.gradient_tolerance > 0.0,
                        "Reduced gradient tolerance must be positive");
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
                        "Reduced relative gradient tolerance must be nonnegative");
      contract::require(parameters_.objective_change_tolerance >= 0.0,
                        "Reduced objective-change tolerance must be nonnegative");
      contract::require(parameters_.step_tolerance >= 0.0,
                        "Reduced step tolerance must be nonnegative");
      contract::require(parameters_.initial_step_length > 0.0,
                        "Reduced gradient initial step must be positive");
      contract::require(parameters_.armijo_fraction > 0.0 &&
                          parameters_.armijo_fraction < 1.0,
                        "Reduced gradient Armijo fraction must lie in (0, 1)");
      contract::require(parameters_.backtracking_factor > 0.0 &&
                          parameters_.backtracking_factor < 1.0,
                        "Reduced gradient backtracking factor must lie in (0, 1)");
    }

    static Result
    result(Primal                         control,
           Evaluation                     final_evaluation,
           std::vector<double>            objective_history,
           std::vector<double>            gradient_norm_history,
           std::vector<double>            relative_gradient_norm_history,
           std::vector<ReducedHessianSolveDiagnostics> hessian_solve_history,
           const std::size_t               accepted_iterations,
           const std::size_t               line_search_trial_count,
           const std::size_t               state_solve_count,
           const std::size_t               adjoint_solve_count,
           const ReducedGradientStoppingReason stopping_reason,
           std::vector<double>            step_length_history,
           std::vector<double>            step_norm_history,
           std::vector<double>            objective_change_history,
           const std::size_t               metric_solve_count,
           const std::size_t               hessian_action_count,
           const std::size_t               direction_reset_count)
    {
      return {std::move(control),
              std::move(final_evaluation),
              std::move(objective_history),
              std::move(gradient_norm_history),
              std::move(relative_gradient_norm_history),
              std::move(hessian_solve_history),
              accepted_iterations,
              line_search_trial_count,
              state_solve_count,
              adjoint_solve_count,
              stopping_reason,
              std::move(step_length_history),
              std::move(step_norm_history),
              std::move(objective_change_history),
              metric_solve_count,
              hessian_action_count,
              direction_reset_count};
    }

    const contract::ReducedDTOT<Backend> &reduced_;
    const contract::MetricT<Backend> &    metric_;
    const contract::ConstraintT<Backend> * constraint_ = nullptr;
    ReducedSolverParameters                parameters_;
    DirectionPolicy                         direction_policy_;
    LineSearchPolicy                        line_search_policy_;
  };

  template <typename Backend>
  using ReducedGradientSolverT =
    ReducedSearchSolverT<Backend,
                         SteepestDescentDirectionPolicyT<Backend>,
                         ArmijoLineSearchPolicyT<Backend>>;

  template <typename Backend>
  using ReducedConjugateGradientSolverT =
    ReducedSearchSolverT<Backend,
                         NonlinearConjugateGradientDirectionPolicyT<Backend>>;

  using ReducedGradientSolver = ReducedGradientSolverT<contract::DenseBackend>;
  using ReducedConjugateGradientSolver =
    ReducedConjugateGradientSolverT<contract::DenseBackend>;

  template <typename Backend>
  using ReducedFletcherReevesSolverT =
    ReducedSearchSolverT<Backend, FletcherReevesDirectionPolicyT<Backend>>;

  using ReducedFletcherReevesSolver =
    ReducedFletcherReevesSolverT<contract::DenseBackend>;

  template <typename Backend>
  using ReducedExactConjugateGradientSolverT =
    ReducedSearchSolverT<Backend,
                         NonlinearConjugateGradientDirectionPolicyT<Backend>,
                         ExactQuadraticLineSearchPolicyT<Backend>>;

  using ReducedExactConjugateGradientSolver =
    ReducedExactConjugateGradientSolverT<contract::DenseBackend>;

  template <typename Backend>
  using ReducedExactFletcherReevesSolverT =
    ReducedSearchSolverT<Backend,
                         FletcherReevesDirectionPolicyT<Backend>,
                         ExactQuadraticLineSearchPolicyT<Backend>>;

  using ReducedExactFletcherReevesSolver =
    ReducedExactFletcherReevesSolverT<contract::DenseBackend>;

  template <typename Backend>
  using ReducedQuadraticConjugateGradientSolverT =
    ReducedSearchSolverT<Backend,
                         QuadraticConjugateGradientDirectionPolicyT<Backend>,
                         ExactQuadraticLineSearchPolicyT<Backend>>;

  using ReducedQuadraticConjugateGradientSolver =
    ReducedQuadraticConjugateGradientSolverT<contract::DenseBackend>;

  template <typename Backend>
  using ReducedLimitedMemoryBfgsSolverT =
    ReducedSearchSolverT<Backend, LimitedMemoryBfgsDirectionPolicyT<Backend>>;

  using ReducedLimitedMemoryBfgsSolver =
    ReducedLimitedMemoryBfgsSolverT<contract::DenseBackend>;

  template <typename Backend>
  using ReducedFullBfgsSolverT =
    ReducedSearchSolverT<Backend, FullBfgsDirectionPolicyT<Backend>>;

  using ReducedFullBfgsSolver =
    ReducedFullBfgsSolverT<contract::DenseBackend>;

  template <typename Backend>
  using ReducedWolfeGradientSolverT =
    ReducedSearchSolverT<Backend,
                         SteepestDescentDirectionPolicyT<Backend>,
                         WolfeLineSearchPolicyT<Backend>>;

  using ReducedWolfeGradientSolver =
    ReducedWolfeGradientSolverT<contract::DenseBackend>;

  template <typename Backend>
  using ReducedWeakWolfeGradientSolverT =
    ReducedSearchSolverT<Backend,
                         SteepestDescentDirectionPolicyT<Backend>,
                         WeakWolfeLineSearchPolicyT<Backend>>;

  using ReducedWeakWolfeGradientSolver =
    ReducedWeakWolfeGradientSolverT<contract::DenseBackend>;

  template <typename Backend>
  using ReducedNewtonSolverT =
    ReducedSearchSolverT<Backend, NewtonDirectionPolicyT<Backend>>;

  using ReducedNewtonSolver = ReducedNewtonSolverT<contract::DenseBackend>;

  template <typename Backend>
  using ReducedExactNewtonSolverT =
    ReducedSearchSolverT<Backend,
                         NewtonDirectionPolicyT<Backend>,
                         ExactQuadraticLineSearchPolicyT<Backend>>;

  using ReducedExactNewtonSolver =
    ReducedExactNewtonSolverT<contract::DenseBackend>;
} // namespace nmopt::solvers
