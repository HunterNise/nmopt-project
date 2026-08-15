#pragma once

#include "nmopt/contract/metric_constraint.hpp"
#include "nmopt/contract/reduced_dto.hpp"
#include "nmopt/contract/reduced_hessian.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <deque>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

namespace nmopt::solvers
{
  enum class ReducedStoppingReason
  {
    gradient_tolerance,
    relative_gradient_tolerance,
    objective_change_tolerance,
    step_tolerance,
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
        case ReducedStoppingReason::relative_gradient_tolerance:
          return "relative_gradient_tolerance";
        case ReducedStoppingReason::objective_change_tolerance:
          return "objective_change_tolerance";
        case ReducedStoppingReason::step_tolerance:
          return "step_tolerance";
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
    // A zero value disables the optional relative/objective/step criteria.
    double       relative_gradient_tolerance = 0.0;
    double       objective_change_tolerance = 0.0;
    double       step_tolerance = 0.0;
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
    std::size_t                     metric_solve_count;
    std::size_t                     hessian_action_count;
  };

  using ReducedSearchDirection =
    ReducedSearchDirectionT<contract::DenseBackend>;

  template <typename Backend>
  struct ReducedSolverResultT
  {
    contract::PrimalBlockT<Backend> control;
    contract::ReducedEvaluationT<Backend> final_evaluation;
    std::vector<double>             objective_history;
    std::vector<double>             gradient_norm_history;
    std::vector<double>             relative_gradient_norm_history;
    std::size_t                     accepted_iterations;
    std::size_t                     line_search_trial_count;
    std::size_t                     state_solve_count;
    std::size_t                     adjoint_solve_count;
    ReducedStoppingReason           stopping_reason;
    std::vector<double>             step_length_history;
    std::vector<double>             step_norm_history;
    std::vector<double>             objective_change_history;
    std::size_t                     metric_solve_count;
    std::size_t                     hessian_action_count;
    std::size_t                     direction_reset_count;
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
  add_scaled_covector(contract::CovectorBlockT<Backend> &      target,
                      const double                             factor,
                      const contract::CovectorBlockT<Backend> &source)
  {
    contract::require_compatible(target,
                                  source,
                                  "Reduced search covector update has incompatible layouts");
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
  make_steepest_descent_direction_from_metric_gradient(
    const contract::CovectorBlockT<Backend> &reduced_derivative,
    ReducedMetricGradientT<Backend>          metric_gradient,
    const std::size_t                        metric_solve_count = 1)
  {
    contract::PrimalBlockT<Backend> direction = metric_gradient.gradient;
    scale_primal(direction, -1.0);
    const double directional_derivative =
      contract::pair(reduced_derivative, direction);
    contract::require(std::isfinite(directional_derivative) &&
                        directional_derivative <= 0.0,
                      "Reduced search metric returned a non-descent direction");

    return {std::move(direction),
            metric_gradient.norm,
            directional_derivative,
            metric_solve_count,
            0};
  }

  template <typename Backend>
  ReducedSearchDirectionT<Backend>
  make_steepest_descent_direction(
    const contract::CovectorBlockT<Backend> &reduced_derivative,
    const contract::MetricT<Backend> &       metric)
  {
    return make_steepest_descent_direction_from_metric_gradient(
      reduced_derivative,
      make_metric_gradient(reduced_derivative, metric));
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
    next(const contract::PrimalBlockT<Backend> &,
         const contract::CovectorBlockT<Backend> &reduced_derivative,
         const contract::MetricT<Backend> &       metric)
    {
      return make_steepest_descent_direction(reduced_derivative, metric);
    }

    std::size_t
    reset_count() const noexcept
    {
      return 0;
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

  enum class NonlinearConjugateGradientUpdate
  {
    polak_ribiere_plus,
    fletcher_reeves,
    classical_quadratic
  };

  inline const char *
  nonlinear_conjugate_gradient_update_name(
    const NonlinearConjugateGradientUpdate update)
  {
    switch (update)
      {
        case NonlinearConjugateGradientUpdate::polak_ribiere_plus:
          return "polak_ribiere_plus";
        case NonlinearConjugateGradientUpdate::fletcher_reeves:
          return "fletcher_reeves";
        case NonlinearConjugateGradientUpdate::classical_quadratic:
          return "classical_quadratic";
      }
    return "unknown";
  }

  // The default nonlinear-CG policy is Polak–Ribière+ in the declared metric.
  // The same typed implementation is specialized below for Fletcher–Reeves
  // and classical quadratic CG. All variants store the previous covector,
  // metric gradient, and primal direction so primal/dual layout compatibility
  // is checked at every update.
  template <typename Backend,
            NonlinearConjugateGradientUpdate Update =
              NonlinearConjugateGradientUpdate::polak_ribiere_plus>
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
    next(const Primal &,
         const Covector &reduced_derivative,
         const contract::MetricT<Backend> &metric)
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

          if constexpr (Update ==
                        NonlinearConjugateGradientUpdate::classical_quadratic)
            {
              if (!restart)
                {
                  const double denominator =
                    contract::pair(*previous_derivative_, *previous_gradient_);
                  const double numerator =
                    contract::pair(reduced_derivative,
                                   current_gradient.gradient);
                  const double scale = std::max(1.0, std::abs(denominator));
                  contract::require(
                    std::isfinite(denominator) &&
                      denominator > parameters_.curvature_tolerance * scale &&
                      std::isfinite(numerator) && numerator > 0.0,
                    "Classical quadratic CG requires positive finite gradient curvature");
                  beta = numerator / denominator;
                }
            }
          else
            {
              const double denominator =
                contract::pair(*previous_derivative_, *previous_gradient_);
              double numerator = 0.0;
              if constexpr (Update ==
                            NonlinearConjugateGradientUpdate::fletcher_reeves)
                numerator =
                  contract::pair(reduced_derivative,
                                 current_gradient.gradient);
              else
                numerator =
                  contract::pair(reduced_derivative,
                                 current_gradient.gradient) -
                  contract::pair(*previous_derivative_,
                                 current_gradient.gradient);
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
          if constexpr (Update ==
                        NonlinearConjugateGradientUpdate::classical_quadratic)
            contract::require(
              std::isfinite(directional_derivative) &&
                directional_derivative < 0.0,
              "Classical quadratic CG could not produce a descent direction");
          else
            {
              restart = true;
              direction = current_gradient.gradient;
              scale_primal(direction, -1.0);
              directional_derivative =
                contract::pair(reduced_derivative, direction);
            }
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
              directional_derivative,
              1,
              0};
    }

    std::size_t
    restart_count() const noexcept
    {
      return restart_count_;
    }

    std::size_t
    reset_count() const noexcept
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

  template <typename Backend>
  using FletcherReevesDirectionPolicyT =
    NonlinearConjugateGradientDirectionPolicyT<
      Backend,
      NonlinearConjugateGradientUpdate::fletcher_reeves>;

  using FletcherReevesDirectionPolicyDense =
    FletcherReevesDirectionPolicyT<contract::DenseBackend>;

  template <typename Backend>
  using QuadraticConjugateGradientDirectionPolicyT =
    NonlinearConjugateGradientDirectionPolicyT<
      Backend,
      NonlinearConjugateGradientUpdate::classical_quadratic>;

  using QuadraticConjugateGradientDirectionPolicyDense =
    QuadraticConjugateGradientDirectionPolicyT<contract::DenseBackend>;

  struct LimitedMemoryBfgsParameters
  {
    unsigned int memory_size = 5;
    double       curvature_tolerance = 1e-14;
  };

  enum class LimitedMemoryBfgsUpdateStatus
  {
    initial,
    accepted_pair,
    curvature_reset
  };

  inline const char *
  limited_memory_bfgs_update_status_name(
    const LimitedMemoryBfgsUpdateStatus status)
  {
    switch (status)
      {
        case LimitedMemoryBfgsUpdateStatus::initial:
          return "initial";
        case LimitedMemoryBfgsUpdateStatus::accepted_pair:
          return "accepted_pair";
        case LimitedMemoryBfgsUpdateStatus::curvature_reset:
          return "curvature_reset";
      }
    return "unknown";
  }

  // The initial inverse-Hessian policy is the declared metric inverse. Each
  // stored pair is (s, y) with primal displacement s and covector difference
  // y, and the two-loop recursion uses only the declared primal-dual pairing.
  template <typename Backend>
  class LimitedMemoryBfgsDirectionPolicyT
  {
  public:
    using Direction = ReducedSearchDirectionT<Backend>;
    using Primal = contract::PrimalBlockT<Backend>;
    using Covector = contract::CovectorBlockT<Backend>;

    LimitedMemoryBfgsDirectionPolicyT(
      LimitedMemoryBfgsParameters parameters = {})
      : parameters_(parameters)
    {
      contract::require(parameters_.memory_size > 0,
                        "L-BFGS memory size must be positive");
      contract::require(parameters_.curvature_tolerance > 0.0,
                        "L-BFGS curvature tolerance must be positive");
    }

    void
    reset()
    {
      previous_control_.reset();
      previous_derivative_.reset();
      history_.clear();
      reset_count_ = 0;
      last_status_ = LimitedMemoryBfgsUpdateStatus::initial;
    }

    Direction
    next(const Primal &control,
         const Covector &reduced_derivative,
         const contract::MetricT<Backend> &metric)
    {
      contract::require(control.layout()->compatible_with(*metric.layout()),
                        "L-BFGS control does not match metric");
      const ReducedMetricGradientT<Backend> current_gradient =
        make_metric_gradient(reduced_derivative, metric);

      if (!previous_control_.has_value())
        {
          previous_control_ = control;
          previous_derivative_ = reduced_derivative;
          last_status_ = LimitedMemoryBfgsUpdateStatus::initial;
          return make_steepest_descent_direction_from_metric_gradient(
            reduced_derivative, current_gradient);
        }

      contract::require(control.layout()->compatible_with(
                          *previous_control_->layout()),
                        "L-BFGS control history has an incompatible layout");
      contract::require(reduced_derivative.layout()->compatible_with(
                          *previous_derivative_->layout()),
                        "L-BFGS covector history has an incompatible layout");

      Primal displacement = control;
      add_scaled_primal(displacement, -1.0, *previous_control_);
      Covector covector_difference =
        contract::subtract(reduced_derivative, *previous_derivative_);
      const double curvature =
        contract::pair(covector_difference, displacement);
      if (!std::isfinite(curvature) ||
          curvature <= parameters_.curvature_tolerance)
        {
          clear_history_after_reset();
          update_previous(control, reduced_derivative);
          last_status_ = LimitedMemoryBfgsUpdateStatus::curvature_reset;
          return make_steepest_descent_direction_from_metric_gradient(
            reduced_derivative, current_gradient);
        }

      if (history_.size() == parameters_.memory_size)
        history_.pop_front();
      history_.push_back(
        {std::move(displacement),
         std::move(covector_difference),
         1.0 / curvature});

      Covector q = reduced_derivative;
      std::vector<double> alpha(history_.size());
      for (std::size_t reverse = history_.size(); reverse > 0; --reverse)
        {
          const std::size_t index = reverse - 1;
          const SecantPair &pair = history_[index];
          alpha[index] = pair.inverse_curvature *
                         contract::pair(q, pair.displacement);
          add_scaled_covector(q, -alpha[index], pair.covector_difference);
        }

      Primal inverse_hessian_action = metric.inverse_apply(q);
      for (std::size_t index = 0; index < history_.size(); ++index)
        {
          const SecantPair &pair = history_[index];
          const double beta = pair.inverse_curvature *
                              contract::pair(pair.covector_difference,
                                             inverse_hessian_action);
          add_scaled_primal(inverse_hessian_action,
                            alpha[index] - beta,
                            pair.displacement);
        }

      scale_primal(inverse_hessian_action, -1.0);
      const double directional_derivative =
        contract::pair(reduced_derivative, inverse_hessian_action);
      if (!std::isfinite(directional_derivative) ||
          directional_derivative >= 0.0)
        {
          clear_history_after_reset();
          update_previous(control, reduced_derivative);
          last_status_ = LimitedMemoryBfgsUpdateStatus::curvature_reset;
          return make_steepest_descent_direction_from_metric_gradient(
            reduced_derivative, current_gradient);
        }

      update_previous(control, reduced_derivative);
      last_status_ = LimitedMemoryBfgsUpdateStatus::accepted_pair;
      return {std::move(inverse_hessian_action),
              current_gradient.norm,
              directional_derivative,
              2,
              0};
    }

    std::size_t
    history_size() const noexcept
    {
      return history_.size();
    }

    std::size_t
    reset_count() const noexcept
    {
      return reset_count_;
    }

    LimitedMemoryBfgsUpdateStatus
    last_update_status() const noexcept
    {
      return last_status_;
    }

  private:
    struct SecantPair
    {
      Primal   displacement;
      Covector covector_difference;
      double   inverse_curvature;
    };

    void
    clear_history_after_reset()
    {
      history_.clear();
      ++reset_count_;
    }

    void
    update_previous(const Primal &control, const Covector &reduced_derivative)
    {
      previous_control_ = control;
      previous_derivative_ = reduced_derivative;
    }

    LimitedMemoryBfgsParameters parameters_;
    std::optional<Primal>       previous_control_;
    std::optional<Covector>     previous_derivative_;
    std::deque<SecantPair>      history_;
    std::size_t                 reset_count_ = 0;
    LimitedMemoryBfgsUpdateStatus last_status_ =
      LimitedMemoryBfgsUpdateStatus::initial;
  };

  template <typename Backend>
  using LimitedMemoryBfgsDirectionPolicy =
    LimitedMemoryBfgsDirectionPolicyT<Backend>;

  using LimitedMemoryBfgsDirectionPolicyDense =
    LimitedMemoryBfgsDirectionPolicyT<contract::DenseBackend>;

  struct ReducedNewtonParameters
  {
    unsigned int maximum_inner_iterations = 100;
    double       relative_tolerance = 1e-10;
    double       absolute_tolerance = 1e-12;
    double       curvature_tolerance = 1e-14;
  };

  // Newton is available only when a model supplies the explicit reduced
  // Hessian capability. The inner solve is metric-preconditioned CG for the
  // covector equation H d = -j'. It never manufactures a Hessian from the
  // first-order DTO ports.
  template <typename Backend>
  class NewtonDirectionPolicyT
  {
  public:
    using Direction = ReducedSearchDirectionT<Backend>;
    using Primal = contract::PrimalBlockT<Backend>;
    using Covector = contract::CovectorBlockT<Backend>;

    NewtonDirectionPolicyT() = default;

    NewtonDirectionPolicyT(
      const contract::ReducedHessianT<Backend> &hessian,
      ReducedNewtonParameters                    parameters = {})
      : hessian_(&hessian)
      , parameters_(parameters)
    {
      validate_parameters();
    }

    void
    reset()
    {}

    Direction
    next(const Primal &control,
         const Covector &reduced_derivative,
         const contract::MetricT<Backend> &metric)
    {
      contract::require(hessian_ != nullptr,
                        "Newton direction requires a reduced Hessian capability");
      validate_parameters();
      contract::require(control.layout()->compatible_with(*metric.layout()),
                        "Newton control does not match metric");
      contract::require(control.layout()->compatible_with(*hessian_->layout()),
                        "Newton control does not match reduced Hessian layout");

      const ReducedMetricGradientT<Backend> current_gradient =
        make_metric_gradient(reduced_derivative, metric);

      Covector residual = reduced_derivative;
      for (std::size_t block = 0; block < residual.n_blocks(); ++block)
        residual.scale_block(block, -1.0);

      Primal preconditioned_residual = metric.inverse_apply(residual);
      std::size_t metric_solve_count = 1;
      const double initial_squared_norm =
        contract::pair(residual, preconditioned_residual);
      contract::require(std::isfinite(initial_squared_norm) &&
                          initial_squared_norm >= 0.0,
                        "Newton inner solve returned an invalid initial residual norm");

      const double initial_norm = std::sqrt(initial_squared_norm);
      const double target_norm = std::max(
        parameters_.absolute_tolerance,
        parameters_.relative_tolerance * initial_norm);
      Primal newton_direction = Primal::zeros(metric.layout());
      std::size_t hessian_action_count = 0;

      if (initial_norm > target_norm)
        {
          Primal search_direction = preconditioned_residual;
          double residual_preconditioned_pairing = initial_squared_norm;
          bool converged = false;

          for (unsigned int iteration = 0;
               iteration < parameters_.maximum_inner_iterations;
               ++iteration)
            {
              const Covector hessian_direction =
                hessian_->apply(control, search_direction);
              ++hessian_action_count;
              contract::require(hessian_direction.layout()->compatible_with(
                                  *metric.layout()),
                                "Reduced Hessian returned an incompatible covector layout");

              const double curvature =
                contract::pair(hessian_direction, search_direction);
              const double direction_norm_squared = contract::pair(
                metric.apply(search_direction), search_direction);
              contract::require(std::isfinite(direction_norm_squared) &&
                                  direction_norm_squared > 0.0,
                                "Newton inner CG produced a zero search direction");
              const double curvature_scale =
                std::max(direction_norm_squared,
                         std::numeric_limits<double>::min());
              contract::require(
                std::isfinite(curvature) &&
                  curvature > parameters_.curvature_tolerance * curvature_scale,
                "Newton reduced Hessian is not positive curvature for CG");

              const double step =
                residual_preconditioned_pairing / curvature;
              contract::require(std::isfinite(step) && step > 0.0,
                                "Newton inner CG returned an invalid step");
              add_scaled_primal(newton_direction, step, search_direction);
              add_scaled_covector(residual, -step, hessian_direction);

              Primal next_preconditioned_residual =
                metric.inverse_apply(residual);
              ++metric_solve_count;
              const double next_pairing =
                contract::pair(residual, next_preconditioned_residual);
              contract::require(std::isfinite(next_pairing) &&
                                  next_pairing >= 0.0,
                                "Newton inner solve returned an invalid residual norm");
              const double next_norm = std::sqrt(next_pairing);
              if (next_norm <= target_norm)
                {
                  converged = true;
                  break;
                }

              contract::require(
                iteration + 1 < parameters_.maximum_inner_iterations,
                "Newton inner CG reached its iteration limit");
              const double beta = next_pairing /
                                  residual_preconditioned_pairing;
              contract::require(std::isfinite(beta) && beta >= 0.0,
                                "Newton inner CG returned an invalid coefficient");
              scale_primal(search_direction, beta);
              add_scaled_primal(search_direction,
                                1.0,
                                next_preconditioned_residual);
              residual_preconditioned_pairing = next_pairing;
            }

          contract::require(converged,
                            "Newton inner CG did not converge");
        }

      const double directional_derivative =
        contract::pair(reduced_derivative, newton_direction);
      contract::require(std::isfinite(directional_derivative) &&
                          directional_derivative <= 0.0,
                        "Newton reduced Hessian did not produce a descent direction");
      return {std::move(newton_direction),
              current_gradient.norm,
              directional_derivative,
              metric_solve_count,
              hessian_action_count};
    }

    std::size_t
    reset_count() const noexcept
    {
      return 0;
    }

  private:
    void
    validate_parameters() const
    {
      contract::require(parameters_.maximum_inner_iterations > 0,
                        "Newton inner iteration limit must be positive");
      contract::require(parameters_.relative_tolerance > 0.0 &&
                          parameters_.relative_tolerance < 1.0,
                        "Newton relative tolerance must lie in (0, 1)");
      contract::require(parameters_.absolute_tolerance > 0.0,
                        "Newton absolute tolerance must be positive");
      contract::require(parameters_.curvature_tolerance > 0.0,
                        "Newton curvature tolerance must be positive");
    }

    const contract::ReducedHessianT<Backend> *hessian_ = nullptr;
    ReducedNewtonParameters                    parameters_;
  };

  template <typename Backend>
  using NewtonDirectionPolicy = NewtonDirectionPolicyT<Backend>;

  using NewtonDirectionPolicyDense =
    NewtonDirectionPolicyT<contract::DenseBackend>;
} // namespace nmopt::solvers
