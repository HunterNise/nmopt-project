#pragma once

#include "nmopt/contract/executable_model.hpp"
#include "nmopt/contract/linear_solve.hpp"
#include "nmopt/contract/metric_constraint.hpp"

#include <cstddef>
#include <functional>
#include <memory>
#include <utility>

namespace nmopt::contract
{
  // This is intentionally the first narrow formulation boundary: one
  // eliminated state block, one decision/control block, and one residual test
  // block. Mixed and multi-equation partitions are later extensions.
  template <typename Backend>
  class StateControlPartitionT
  {
  public:
    StateControlPartitionT(const ExecutableModelT<Backend> &model,
                           const std::size_t                state_block,
                           const std::size_t                control_block)
      : variable_layout_(model.variable_layout())
      , state_block_(state_block)
      , control_block_(control_block)
    {
      require(variable_layout_->n_blocks() == 2,
              "The v0 reduced DTO contract requires exactly two variable blocks");
      require(model.test_layout()->n_blocks() == 1,
              "The v0 reduced DTO contract requires exactly one test block");
      require(state_block_ != control_block_,
              "State and control blocks must be distinct");
      require(state_block_ < variable_layout_->n_blocks() &&
                control_block_ < variable_layout_->n_blocks(),
              "State/control block index is outside the variable layout");
      state_layout_ = variable_layout_->single_block(state_block_, "state");
      control_layout_ =
        variable_layout_->single_block(control_block_, "control");
    }

    const LayoutPtr &
    state_layout() const
    {
      return state_layout_;
    }

    const LayoutPtr &
    control_layout() const
    {
      return control_layout_;
    }

    PrimalBlockT<Backend>
    compose(const PrimalBlockT<Backend> &state,
            const PrimalBlockT<Backend> &control) const
    {
      require(state.layout()->compatible_with(*state_layout_),
              "State does not match the selected state block");
      require(control.layout()->compatible_with(*control_layout_),
              "Control does not match the selected control block");

      std::vector<typename Backend::Vector> blocks(variable_layout_->n_blocks());
      blocks.at(state_block_)   = state.block(0);
      blocks.at(control_block_) = control.block(0);
      return PrimalBlockT<Backend>(variable_layout_, std::move(blocks));
    }

    CovectorBlockT<Backend>
    state_component(const CovectorBlockT<Backend> &full_covector) const
    {
      require(full_covector.layout()->compatible_with(*variable_layout_),
              "Covector does not match the variable layout");
      return extract_covector_block(full_covector, state_block_, "state");
    }

    CovectorBlockT<Backend>
    control_component(const CovectorBlockT<Backend> &full_covector) const
    {
      require(full_covector.layout()->compatible_with(*variable_layout_),
              "Covector does not match the variable layout");
      return extract_covector_block(full_covector, control_block_, "control");
    }

  private:
    LayoutPtr   variable_layout_;
    std::size_t state_block_;
    std::size_t control_block_;
    LayoutPtr   state_layout_;
    LayoutPtr   control_layout_;
  };

  using StateControlPartition = StateControlPartitionT<DenseBackend>;

  template <typename Backend>
  struct StateAdjointSolversT
  {
    // Input: control. Output: state satisfying the compiled residual.
    std::function<FormulationSolveResultT<Backend>(
      const PrimalBlockT<Backend> &)>
      solve_state;

    // Input: full point and covector D_y J. Output: test-space primal p
    // satisfying E_y^* p = D_y J.
    std::function<FormulationSolveResultT<Backend>(
      const PrimalBlockT<Backend> &,
      const CovectorBlockT<Backend> &)>
      solve_adjoint;
  };

  using StateAdjointSolvers = StateAdjointSolversT<DenseBackend>;

  template <typename Backend>
  class ReducedDTOT;

  template <typename Backend>
  struct ReducedValueEvaluationT
  {
    PrimalBlockT<Backend> state;
    PrimalBlockT<Backend> full_point;
    PrimalBlockT<Backend> control;
    double                objective_value;
    LinearSolveReport     state_solve;

  private:
    friend class ReducedDTOT<Backend>;

    ReducedValueEvaluationT(
      PrimalBlockT<Backend>       state_value,
      PrimalBlockT<Backend>       full_point_value,
      PrimalBlockT<Backend>       control_value,
      const double                objective,
      LinearSolveReport           state_report,
      std::shared_ptr<const void> evaluation_token,
      std::shared_ptr<const void> lifetime_owner)
      : state(std::move(state_value))
      , full_point(std::move(full_point_value))
      , control(std::move(control_value))
      , objective_value(objective)
      , state_solve(std::move(state_report))
      , evaluation_token_(std::move(evaluation_token))
      , lifetime_owner_(std::move(lifetime_owner))
    {}

    std::shared_ptr<const void> evaluation_token_;
    std::shared_ptr<const void> lifetime_owner_;
  };

  using ReducedValueEvaluation = ReducedValueEvaluationT<DenseBackend>;

  template <typename Backend>
  struct ReducedEvaluationT
  {
    PrimalBlockT<Backend>   state;
    PrimalBlockT<Backend>   adjoint;
    PrimalBlockT<Backend>   full_point;
    CovectorBlockT<Backend> reduced_derivative;
    double                  objective_value;
    LinearSolveReport       state_solve;
    LinearSolveReport       adjoint_solve;
  };

  using ReducedEvaluation = ReducedEvaluationT<DenseBackend>;

  template <typename Backend>
  class ReducedDTOT
  {
  public:
    ReducedDTOT(const ExecutableModelT<Backend> &model,
                StateControlPartitionT<Backend>  partition,
                StateAdjointSolversT<Backend>    solvers)
      : model_(&model)
      , partition_(std::move(partition))
      , solvers_(std::move(solvers))
      , evaluation_token_(std::make_shared<const std::size_t>(0))
    {
      require(static_cast<bool>(solvers_.solve_state),
              "Reduced DTO requires a state solve operation");
      require(static_cast<bool>(solvers_.solve_adjoint),
              "Reduced DTO requires an adjoint solve operation");
    }

    // A compiled reduced service owns its executable and any backend session
    // token needed by that executable. The reference-taking constructor above
    // remains available for deliberately short-lived direct use.
    ReducedDTOT(std::shared_ptr<const ExecutableModelT<Backend>> model,
                StateControlPartitionT<Backend>                  partition,
                StateAdjointSolversT<Backend>                    solvers,
                std::shared_ptr<const void> lifetime_owner = {})
      : owned_model_(std::move(model))
      , model_(owned_model_.get())
      , partition_(std::move(partition))
      , solvers_(std::move(solvers))
      , lifetime_owner_(std::move(lifetime_owner))
      , evaluation_token_(std::make_shared<const std::size_t>(0))
    {
      require(static_cast<bool>(owned_model_),
              "Owned reduced DTO requires an executable model");
      require(static_cast<bool>(solvers_.solve_state),
              "Reduced DTO requires a state solve operation");
      require(static_cast<bool>(solvers_.solve_adjoint),
              "Reduced DTO requires an adjoint solve operation");
    }

    ReducedValueEvaluationT<Backend>
    evaluate_value(const PrimalBlockT<Backend> &control) const
    {
      require(control.layout()->compatible_with(*partition_.control_layout()),
              "Control does not match the reduced DTO control layout");

      FormulationSolveResultT<Backend> state_result =
        solvers_.solve_state(control);
      require(state_result.report.converged(),
              "State solve did not converge under its declared policy");
      PrimalBlockT<Backend> state = std::move(state_result.solution);
      require(state.layout()->compatible_with(*partition_.state_layout()),
              "State solver returned an incompatible state layout");

      PrimalBlockT<Backend> full_point = partition_.compose(state, control);
      const double objective_value = model_->objective(full_point);

      return ReducedValueEvaluationT<Backend>(
        std::move(state),
        std::move(full_point),
        control,
        objective_value,
        std::move(state_result.report),
        evaluation_token_,
        lifetime_owner_);
    }

    ReducedEvaluationT<Backend>
    augment_derivative(const ReducedValueEvaluationT<Backend> &value) const
    {
      require(value.evaluation_token_ == evaluation_token_,
              "Reduced DTO value evaluation belongs to another service");
      require(value.lifetime_owner_ == lifetime_owner_,
              "Reduced DTO value evaluation has an incompatible lifetime owner");
      require(value.control.layout()->compatible_with(
                *partition_.control_layout()),
              "Reduced DTO value evaluation has an incompatible control layout");
      require(value.state.layout()->compatible_with(*partition_.state_layout()),
              "Reduced DTO value evaluation has an incompatible state layout");
      require(value.full_point.layout()->compatible_with(
                *model_->variable_layout()),
              "Reduced DTO value evaluation has an incompatible full-point layout");
      require(same_control(value.state, value.control, value.full_point),
              "Reduced DTO value evaluation control does not match its full point");

      CovectorBlockT<Backend> objective_derivative =
        model_->objective_derivative(value.full_point);
      CovectorBlockT<Backend> state_rhs =
        partition_.state_component(objective_derivative);

      FormulationSolveResultT<Backend> adjoint_result =
        solvers_.solve_adjoint(value.full_point, state_rhs);
      require(adjoint_result.report.converged(),
              "Adjoint solve did not converge under its declared policy");
      PrimalBlockT<Backend> adjoint = std::move(adjoint_result.solution);
      require(adjoint.layout()->compatible_with(*model_->test_layout()),
              "Adjoint solver returned an incompatible test-space layout");

      CovectorBlockT<Backend> residual_pullback =
        model_->residual_vjp(value.full_point, adjoint);
      CovectorBlockT<Backend> reduced_derivative =
        subtract(partition_.control_component(objective_derivative),
                 partition_.control_component(residual_pullback));

      return {value.state,
              std::move(adjoint),
              value.full_point,
              std::move(reduced_derivative),
              value.objective_value,
              value.state_solve,
              std::move(adjoint_result.report)};
    }

    ReducedEvaluationT<Backend>
    evaluate(const PrimalBlockT<Backend> &control) const
    {
      return augment_derivative(evaluate_value(control));
    }

    PrimalBlockT<Backend>
    gradient_direction(const CovectorBlockT<Backend> &reduced_derivative,
                       const MetricT<Backend>         &metric) const
    {
      require(reduced_derivative.layout()->compatible_with(
                *partition_.control_layout()),
              "Reduced derivative does not match the control layout");
      require(metric.layout()->compatible_with(*partition_.control_layout()),
              "Metric does not match the control layout");
      return metric.inverse_apply(reduced_derivative);
    }

  private:
    bool
    same_control(const PrimalBlockT<Backend> &state,
                 const PrimalBlockT<Backend> &control,
                 const PrimalBlockT<Backend> &full_point) const
    {
      require(state.layout()->compatible_with(*partition_.state_layout()),
              "Reduced DTO value state does not match the state layout");
      require(control.layout()->compatible_with(*partition_.control_layout()),
              "Reduced DTO value control does not match the control layout");
      require(full_point.layout()->compatible_with(*model_->variable_layout()),
              "Reduced DTO value full point does not match the variable layout");

      const PrimalBlockT<Backend> expected_full_point =
        partition_.compose(state, control);
      if (!expected_full_point.layout()->compatible_with(*full_point.layout()))
        return false;

      for (std::size_t block = 0; block < full_point.n_blocks(); ++block)
        {
          typename Backend::Vector difference = full_point.block(block);
          Backend::add_scaled(difference,
                              -1.0,
                              expected_full_point.block(block));
          if (Backend::dot(difference, difference) != 0.0)
            return false;
        }
      return true;
    }

    std::shared_ptr<const ExecutableModelT<Backend>> owned_model_;
    const ExecutableModelT<Backend> *                model_ = nullptr;
    StateControlPartitionT<Backend>                  partition_;
    StateAdjointSolversT<Backend>                    solvers_;
    std::shared_ptr<const void>                      lifetime_owner_;
    std::shared_ptr<const void>                      evaluation_token_;
  };

  using ReducedDTO = ReducedDTOT<DenseBackend>;
} // namespace nmopt::contract
