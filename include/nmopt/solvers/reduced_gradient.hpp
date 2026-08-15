#pragma once

#include "nmopt/contract/reduced_dto.hpp"
#include "nmopt/solvers/reduced_search.hpp"

#include <cmath>
#include <cstddef>
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
              SteepestDescentDirectionPolicyT<Backend>>
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

      Primal current_control = initial_control;
      Evaluation current_evaluation = reduced_.evaluate(current_control);
      ++state_solve_count;
      ++adjoint_solve_count;

      std::vector<double> objective_history{current_evaluation.objective_value};
      std::vector<double> gradient_norm_history;
      std::vector<double> step_length_history;
      std::vector<double> objective_change_history;
      DirectionPolicy direction_policy = direction_policy_;
      direction_policy.reset();

      for (;;)
        {
          const ReducedSearchDirectionT<Backend> direction =
            direction_policy.next(current_control,
                                  current_evaluation.reduced_derivative,
                                  metric_);
          ++metric_solve_count;
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
          gradient_norm_history.push_back(stopping_norm);

          if (stopping_norm <= parameters_.gradient_tolerance)
            return result(current_control,
                          std::move(objective_history),
                          std::move(gradient_norm_history),
                          accepted_iterations,
                          line_search_trial_count,
                          state_solve_count,
                          adjoint_solve_count,
                          ReducedGradientStoppingReason::gradient_tolerance,
                          std::move(step_length_history),
                          std::move(objective_change_history),
                          metric_solve_count,
                          direction_policy.reset_count());

          contract::require(unit_step_descent_measure < 0.0,
                            "Reduced gradient did not produce a descent direction");

          if (accepted_iterations == parameters_.maximum_iterations)
            return result(current_control,
                          std::move(objective_history),
                          std::move(gradient_norm_history),
                          accepted_iterations,
                          line_search_trial_count,
                          state_solve_count,
                          adjoint_solve_count,
                          ReducedGradientStoppingReason::maximum_iterations,
                          std::move(step_length_history),
                          std::move(objective_change_history),
                          metric_solve_count,
                          direction_policy.reset_count());

          double step_length = parameters_.initial_step_length;
          bool   accepted = false;
          for (unsigned int trial = 0;
               trial < parameters_.maximum_line_search_trials;
               ++trial)
            {
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
              Evaluation trial_evaluation = reduced_.evaluate(trial_control);
              ++state_solve_count;
              ++adjoint_solve_count;
              ++line_search_trial_count;

              Primal update = trial_control;
              add_scaled_primal(update, -1.0, current_control);
              const double armijo_decrease =
                contract::pair(current_evaluation.reduced_derivative, update);
              const double armijo_bound = current_evaluation.objective_value +
                                          parameters_.armijo_fraction *
                                            armijo_decrease;
              if (armijo_decrease < 0.0 &&
                  trial_evaluation.objective_value <= armijo_bound)
                {
                  const double objective_change =
                    trial_evaluation.objective_value -
                    current_evaluation.objective_value;
                  current_control = std::move(trial_control);
                  current_evaluation = std::move(trial_evaluation);
                  objective_history.push_back(current_evaluation.objective_value);
                  step_length_history.push_back(step_length);
                  objective_change_history.push_back(objective_change);
                  ++accepted_iterations;
                  accepted = true;
                  break;
                }
              step_length *= parameters_.backtracking_factor;
            }

          if (!accepted)
            return result(current_control,
                          std::move(objective_history),
                          std::move(gradient_norm_history),
                          accepted_iterations,
                          line_search_trial_count,
                          state_solve_count,
                          adjoint_solve_count,
                          ReducedGradientStoppingReason::line_search_failure,
                          std::move(step_length_history),
                          std::move(objective_change_history),
                          metric_solve_count,
                          direction_policy.reset_count());
        }
    }

  private:
    struct ProjectedGradientData
    {
      double metric_norm;
      double descent_measure;
    };

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
           std::vector<double>            objective_history,
           std::vector<double>            gradient_norm_history,
           const std::size_t               accepted_iterations,
           const std::size_t               line_search_trial_count,
           const std::size_t               state_solve_count,
           const std::size_t               adjoint_solve_count,
           const ReducedGradientStoppingReason stopping_reason,
           std::vector<double>            step_length_history,
           std::vector<double>            objective_change_history,
           const std::size_t               metric_solve_count,
           const std::size_t               direction_reset_count)
    {
      return {std::move(control),
              std::move(objective_history),
              std::move(gradient_norm_history),
              accepted_iterations,
              line_search_trial_count,
              state_solve_count,
              adjoint_solve_count,
              stopping_reason,
              std::move(step_length_history),
              std::move(objective_change_history),
              metric_solve_count,
              direction_reset_count};
    }

    const contract::ReducedDTOT<Backend> &reduced_;
    const contract::MetricT<Backend> &    metric_;
    const contract::ConstraintT<Backend> * constraint_ = nullptr;
    ReducedSolverParameters                parameters_;
    DirectionPolicy                         direction_policy_;
  };

  template <typename Backend>
  using ReducedGradientSolverT =
    ReducedSearchSolverT<Backend, SteepestDescentDirectionPolicyT<Backend>>;

  template <typename Backend>
  using ReducedConjugateGradientSolverT =
    ReducedSearchSolverT<Backend,
                         NonlinearConjugateGradientDirectionPolicyT<Backend>>;

  using ReducedGradientSolver = ReducedGradientSolverT<contract::DenseBackend>;
  using ReducedConjugateGradientSolver =
    ReducedConjugateGradientSolverT<contract::DenseBackend>;

  template <typename Backend>
  using ReducedLimitedMemoryBfgsSolverT =
    ReducedSearchSolverT<Backend, LimitedMemoryBfgsDirectionPolicyT<Backend>>;

  using ReducedLimitedMemoryBfgsSolver =
    ReducedLimitedMemoryBfgsSolverT<contract::DenseBackend>;
} // namespace nmopt::solvers
