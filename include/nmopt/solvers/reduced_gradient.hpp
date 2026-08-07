#pragma once

#include "nmopt/contract/reduced_dto.hpp"

#include <cmath>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace nmopt::solvers
{
  enum class ReducedGradientStoppingReason
  {
    gradient_tolerance,
    maximum_iterations,
    line_search_failure
  };

  inline const char *
  reduced_gradient_stopping_reason_name(
    const ReducedGradientStoppingReason reason)
  {
    switch (reason)
      {
        case ReducedGradientStoppingReason::gradient_tolerance:
          return "gradient_tolerance";
        case ReducedGradientStoppingReason::maximum_iterations:
          return "maximum_iterations";
        case ReducedGradientStoppingReason::line_search_failure:
          return "line_search_failure";
      }
    return "unknown";
  }

  struct ReducedGradientParameters
  {
    unsigned int maximum_iterations = 100;
    unsigned int maximum_line_search_trials = 20;
    double       gradient_tolerance = 1e-8;
    double       initial_step_length = 1.0;
    double       armijo_fraction = 1e-4;
    double       backtracking_factor = 0.5;
  };

  template <typename Backend>
  struct ReducedGradientResultT
  {
    contract::PrimalBlockT<Backend> control;
    std::vector<double>             objective_history;
    std::vector<double>             gradient_norm_history;
    std::size_t                     accepted_iterations;
    std::size_t                     line_search_trial_count;
    std::size_t                     state_solve_count;
    std::size_t                     adjoint_solve_count;
    ReducedGradientStoppingReason   stopping_reason;
  };

  using ReducedGradientResult = ReducedGradientResultT<contract::DenseBackend>;

  // An unconstrained reduced-space method. The metric identifies the reduced
  // covector j' with g = G^{-1} j'; Armijo accepts u - alpha g when
  // j(u - alpha g) <= j(u) - c alpha <j', g>.
  template <typename Backend>
  class ReducedGradientSolverT
  {
  public:
    using Primal = contract::PrimalBlockT<Backend>;
    using Covector = contract::CovectorBlockT<Backend>;
    using Evaluation = contract::ReducedEvaluationT<Backend>;
    using Result = ReducedGradientResultT<Backend>;

    ReducedGradientSolverT(const contract::ReducedDTOT<Backend> &reduced,
                           const contract::MetricT<Backend> &    metric,
                           ReducedGradientParameters              parameters = {})
      : reduced_(reduced)
      , metric_(metric)
      , parameters_(parameters)
    {
      validate_parameters();
    }

    ReducedGradientSolverT(const contract::ReducedDTOT<Backend> &reduced,
                           const contract::MetricT<Backend> &    metric,
                           const contract::ConstraintT<Backend> & constraint,
                           ReducedGradientParameters              parameters = {})
      : reduced_(reduced)
      , metric_(metric)
      , constraint_(&constraint)
      , parameters_(parameters)
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

      Primal current_control = initial_control;
      Evaluation current_evaluation = reduced_.evaluate(current_control);
      ++state_solve_count;
      ++adjoint_solve_count;

      std::vector<double> objective_history{current_evaluation.objective_value};
      std::vector<double> gradient_norm_history;

      for (;;)
        {
          const GradientData gradient =
            evaluate_gradient(current_evaluation.reduced_derivative);
          double stopping_norm = gradient.metric_norm;
          double unit_step_descent_measure = -gradient.descent_measure;
          if (constraint_ != nullptr)
            {
              const ProjectedGradientData projected_gradient =
                evaluate_projected_gradient(current_control,
                                            current_evaluation.reduced_derivative,
                                            gradient.direction);
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
                          ReducedGradientStoppingReason::gradient_tolerance);

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
                          ReducedGradientStoppingReason::maximum_iterations);

          double step_length = parameters_.initial_step_length;
          bool   accepted = false;
          for (unsigned int trial = 0;
               trial < parameters_.maximum_line_search_trials;
               ++trial)
            {
              Primal trial_control = current_control;
              add_scaled(trial_control, -step_length, gradient.direction);
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
              add_scaled(update, -1.0, current_control);
              const double armijo_decrease =
                contract::pair(current_evaluation.reduced_derivative, update);
              const double armijo_bound = current_evaluation.objective_value +
                                          parameters_.armijo_fraction *
                                            armijo_decrease;
              if (armijo_decrease < 0.0 &&
                  trial_evaluation.objective_value <= armijo_bound)
                {
                  current_control = std::move(trial_control);
                  current_evaluation = std::move(trial_evaluation);
                  objective_history.push_back(current_evaluation.objective_value);
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
                          ReducedGradientStoppingReason::line_search_failure);
        }
    }

  private:
    struct GradientData
    {
      Primal direction;
      double metric_norm;
      double descent_measure;
    };

    struct ProjectedGradientData
    {
      double metric_norm;
      double descent_measure;
    };

    GradientData
    evaluate_gradient(const Covector &reduced_derivative) const
    {
      Primal gradient = reduced_.gradient_direction(reduced_derivative, metric_);
      const Covector metric_gradient = metric_.apply(gradient);
      const double metric_norm_squared = contract::pair(metric_gradient, gradient);
      contract::require(metric_norm_squared >= 0.0,
                        "Reduced gradient metric returned a negative squared norm");

      const double descent_measure = contract::pair(reduced_derivative, gradient);
      contract::require(descent_measure >= 0.0,
                        "Reduced gradient metric returned a negative descent measure");
      return {std::move(gradient), std::sqrt(metric_norm_squared), descent_measure};
    }

    ProjectedGradientData
    evaluate_projected_gradient(const Primal &  control,
                                const Covector &reduced_derivative,
                                const Primal &  gradient) const
    {
      Primal projected_control = control;
      add_scaled(projected_control, -1.0, gradient);
      projected_control = constraint_->project_in(projected_control, metric_);
      contract::require(constraint_->is_feasible(projected_control),
                        "Reduced gradient projection returned an infeasible control");

      Primal update = projected_control;
      add_scaled(update, -1.0, control);
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

    static void
    add_scaled(Primal &target, const double factor, const Primal &source)
    {
      contract::require_compatible(target,
                                   source,
                                   "Reduced gradient update has incompatible layouts");
      for (std::size_t block = 0; block < target.n_blocks(); ++block)
        target.add_scaled_block(block, factor, source.block(block));
    }

    static Result
    result(Primal                         control,
           std::vector<double>            objective_history,
           std::vector<double>            gradient_norm_history,
           const std::size_t               accepted_iterations,
           const std::size_t               line_search_trial_count,
           const std::size_t               state_solve_count,
           const std::size_t               adjoint_solve_count,
           const ReducedGradientStoppingReason stopping_reason)
    {
      return {std::move(control),
              std::move(objective_history),
              std::move(gradient_norm_history),
              accepted_iterations,
              line_search_trial_count,
              state_solve_count,
              adjoint_solve_count,
              stopping_reason};
    }

    const contract::ReducedDTOT<Backend> &reduced_;
    const contract::MetricT<Backend> &    metric_;
    const contract::ConstraintT<Backend> * constraint_ = nullptr;
    ReducedGradientParameters              parameters_;
  };

  using ReducedGradientSolver = ReducedGradientSolverT<contract::DenseBackend>;
} // namespace nmopt::solvers
