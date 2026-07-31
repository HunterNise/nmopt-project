#pragma once

#include "nmopt/contract/executable_model.hpp"
#include "nmopt/contract/metric_constraint.hpp"

#include <functional>
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
    std::function<PrimalBlockT<Backend>(const PrimalBlockT<Backend> &)>
      solve_state;

    // Input: full point and covector D_y J. Output: test-space primal p
    // satisfying E_y^* p = D_y J.
    std::function<PrimalBlockT<Backend>(const PrimalBlockT<Backend> &,
                                        const CovectorBlockT<Backend> &)>
      solve_adjoint;
  };

  using StateAdjointSolvers = StateAdjointSolversT<DenseBackend>;

  template <typename Backend>
  struct ReducedEvaluationT
  {
    PrimalBlockT<Backend>   state;
    PrimalBlockT<Backend>   adjoint;
    PrimalBlockT<Backend>   full_point;
    CovectorBlockT<Backend> reduced_derivative;
    double                  objective_value;
  };

  using ReducedEvaluation = ReducedEvaluationT<DenseBackend>;

  template <typename Backend>
  class ReducedDTOT
  {
  public:
    ReducedDTOT(const ExecutableModelT<Backend> &model,
                StateControlPartitionT<Backend>  partition,
                StateAdjointSolversT<Backend>    solvers)
      : model_(model)
      , partition_(std::move(partition))
      , solvers_(std::move(solvers))
    {
      require(static_cast<bool>(solvers_.solve_state),
              "Reduced DTO requires a state solve operation");
      require(static_cast<bool>(solvers_.solve_adjoint),
              "Reduced DTO requires an adjoint solve operation");
    }

    ReducedEvaluationT<Backend>
    evaluate(const PrimalBlockT<Backend> &control) const
    {
      require(control.layout()->compatible_with(*partition_.control_layout()),
              "Control does not match the reduced DTO control layout");

      PrimalBlockT<Backend> state = solvers_.solve_state(control);
      require(state.layout()->compatible_with(*partition_.state_layout()),
              "State solver returned an incompatible state layout");

      PrimalBlockT<Backend> full_point = partition_.compose(state, control);
      CovectorBlockT<Backend> objective_derivative =
        model_.objective_derivative(full_point);
      CovectorBlockT<Backend> state_rhs =
        partition_.state_component(objective_derivative);

      PrimalBlockT<Backend> adjoint =
        solvers_.solve_adjoint(full_point, state_rhs);
      require(adjoint.layout()->compatible_with(*model_.test_layout()),
              "Adjoint solver returned an incompatible test-space layout");

      CovectorBlockT<Backend> residual_pullback =
        model_.residual_vjp(full_point, adjoint);
      CovectorBlockT<Backend> reduced_derivative =
        subtract(partition_.control_component(objective_derivative),
                 partition_.control_component(residual_pullback));
      const double objective_value = model_.objective(full_point);

      return {std::move(state),
              std::move(adjoint),
              std::move(full_point),
              std::move(reduced_derivative),
              objective_value};
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
    const ExecutableModelT<Backend> &model_;
    StateControlPartitionT<Backend>  partition_;
    StateAdjointSolversT<Backend>    solvers_;
  };

  using ReducedDTO = ReducedDTOT<DenseBackend>;
} // namespace nmopt::contract
