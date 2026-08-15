#pragma once

#include "nmopt/contract/executable_model.hpp"
#include "nmopt/contract/reduced_hessian.hpp"
#include "nmopt/reference/quadratic_kkt.hpp"

#include <utility>

namespace nmopt::reference
{
  using namespace nmopt::contract;

  // Reference DTO model:
  //
  // r(y,u) = A y - f - B u,
  // J(y,u) = 1/2 ||C y-d||_W^2 + alpha/2 ||u||_R^2.
  //
  // This is an executable algebraic stand-in for the first FE compiler output.
  // It is intentionally not a PDE compiler or an alternative discretisation.
  class LinearQuadraticModel final : public ExecutableModel,
                                     public ReducedHessian
  {
  public:
    LinearQuadraticModel(DenseMatrix A,
                         DenseMatrix B,
                         DenseVector f,
                         DenseMatrix C,
                         DenseVector desired_observation,
                         DenseVector observation_weights,
                         DenseVector regularisation_weights,
                         const double alpha)
      : A_(std::move(A))
      , B_(std::move(B))
      , f_(std::move(f))
      , C_(std::move(C))
      , desired_observation_(std::move(desired_observation))
      , observation_weights_(std::move(observation_weights))
      , regularisation_weights_(std::move(regularisation_weights))
      , alpha_(alpha)
      , variable_layout_(std::make_shared<const BlockLayout>(
          "variables",
          std::vector<SpaceId>{{"state"}, {"control"}},
          std::vector<std::size_t>{A_.columns(), B_.columns()}))
      , test_layout_(std::make_shared<const BlockLayout>(
          "state_test",
          std::vector<SpaceId>{{"state_test"}},
          std::vector<std::size_t>{A_.rows()}))
      , state_layout_(variable_layout_->single_block(0, "state"))
      , control_layout_(variable_layout_->single_block(1, "control"))
    {
      require(A_.rows() == A_.columns(),
              "Reference model state matrix must be square");
      require(B_.rows() == A_.rows(),
              "Reference model control coupling has the wrong row count");
      require(f_.size() == A_.rows(),
              "Reference model load has the wrong size");
      require(C_.columns() == A_.columns(),
              "Reference model observation has the wrong column count");
      require(desired_observation_.size() == C_.rows() &&
                observation_weights_.size() == C_.rows(),
              "Reference model observation data has the wrong size");
      require(regularisation_weights_.size() == B_.columns(),
              "Reference model regularisation has the wrong size");
      require(alpha_ > 0.0, "Reference model alpha must be positive");
      for (std::size_t index = 0; index < observation_weights_.size(); ++index)
        require(observation_weights_[index] > 0.0,
                "Reference model observation weights must be positive");
      for (std::size_t index = 0; index < regularisation_weights_.size();
           ++index)
        require(regularisation_weights_[index] > 0.0,
                "Reference model regularisation weights must be positive");
    }

    const LayoutPtr &
    variable_layout() const override
    {
      return variable_layout_;
    }

    const LayoutPtr &
    test_layout() const override
    {
      return test_layout_;
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

    const LayoutPtr &
    layout() const override
    {
      return control_layout_;
    }

    LinearQuadraticKKTData
    kkt_data() const
    {
      return LinearQuadraticKKTData(A_,
                                    B_,
                                    f_,
                                    C_,
                                    desired_observation_,
                                    observation_weights_,
                                    regularisation_weights_,
                                    alpha_);
    }

    CovectorBlock
    apply(const PrimalBlock &control,
          const PrimalBlock &direction) const override
    {
      require(control.layout()->compatible_with(*control_layout_) &&
                direction.layout()->compatible_with(*control_layout_),
              "Reference reduced Hessian received incompatible controls");

      DenseVector state_tangent =
        A_.solve(B_.vmult(direction.block(0)));
      DenseVector observation = C_.vmult(state_tangent);
      for (std::size_t index = 0; index < observation.size(); ++index)
        observation[index] *= observation_weights_[index];

      DenseVector incremental_adjoint =
        A_.transpose().solve(C_.transpose_vmult(observation));
      DenseVector control_action = B_.transpose_vmult(incremental_adjoint);
      for (std::size_t index = 0; index < control_action.size(); ++index)
        control_action[index] += alpha_ * regularisation_weights_[index] *
                                 direction.block(0)[index];

      return CovectorBlock(control_layout_, {std::move(control_action)});
    }

    CovectorBlock
    residual(const PrimalBlock &variables) const override
    {
      require(variables.layout()->compatible_with(*variable_layout_),
              "Reference residual received incompatible variables");
      DenseVector value = A_.vmult(variables.block(0));
      value.add_scaled(-1.0, f_);
      value.add_scaled(-1.0, B_.vmult(variables.block(1)));
      return CovectorBlock(test_layout_, {std::move(value)});
    }

    CovectorBlock
    residual_jvp(const PrimalBlock &variables,
                 const PrimalBlock &variable_tangent) const override
    {
      require(variables.layout()->compatible_with(*variable_layout_) &&
                variable_tangent.layout()->compatible_with(*variable_layout_),
              "Reference residual JVP received incompatible variables");
      DenseVector value = A_.vmult(variable_tangent.block(0));
      value.add_scaled(-1.0, B_.vmult(variable_tangent.block(1)));
      return CovectorBlock(test_layout_, {std::move(value)});
    }

    CovectorBlock
    residual_vjp(const PrimalBlock &variables,
                 const PrimalBlock &test_seed) const override
    {
      require(variables.layout()->compatible_with(*variable_layout_) &&
                test_seed.layout()->compatible_with(*test_layout_),
              "Reference residual VJP received incompatible inputs");
      DenseVector state = A_.transpose_vmult(test_seed.block(0));
      DenseVector control = B_.transpose_vmult(test_seed.block(0));
      control.scale(-1.0);
      return CovectorBlock(variable_layout_,
                           {std::move(state), std::move(control)});
    }

    double
    objective(const PrimalBlock &variables) const override
    {
      require(variables.layout()->compatible_with(*variable_layout_),
              "Reference objective received incompatible variables");

      DenseVector observation = C_.vmult(variables.block(0));
      observation.add_scaled(-1.0, desired_observation_);

      double value = 0.0;
      for (std::size_t index = 0; index < observation.size(); ++index)
        value += 0.5 * observation_weights_[index] * observation[index] *
                 observation[index];
      for (std::size_t index = 0; index < variables.block(1).size(); ++index)
        value += 0.5 * alpha_ * regularisation_weights_[index] *
                 variables.block(1)[index] * variables.block(1)[index];
      return value;
    }

    CovectorBlock
    objective_derivative(const PrimalBlock &variables) const override
    {
      require(variables.layout()->compatible_with(*variable_layout_),
              "Reference objective derivative received incompatible variables");

      DenseVector observation = C_.vmult(variables.block(0));
      observation.add_scaled(-1.0, desired_observation_);
      for (std::size_t index = 0; index < observation.size(); ++index)
        observation[index] *= observation_weights_[index];

      DenseVector state = C_.transpose_vmult(observation);
      DenseVector control = variables.block(1);
      for (std::size_t index = 0; index < control.size(); ++index)
        control[index] *= alpha_ * regularisation_weights_[index];

      return CovectorBlock(variable_layout_,
                           {std::move(state), std::move(control)});
    }

    PrimalBlock
    solve_state(const PrimalBlock &control) const
    {
      require(control.layout()->compatible_with(*control_layout_),
              "Reference state solve received incompatible control");
      DenseVector right_hand_side = B_.vmult(control.block(0));
      right_hand_side.add_scaled(1.0, f_);
      return PrimalBlock(state_layout_, {A_.solve(std::move(right_hand_side))});
    }

    PrimalBlock
    solve_adjoint(const PrimalBlock &full_point,
                  const CovectorBlock &state_objective_derivative) const
    {
      require(full_point.layout()->compatible_with(*variable_layout_) &&
                state_objective_derivative.layout()->compatible_with(
                  *state_layout_),
              "Reference adjoint solve received incompatible arguments");
      return PrimalBlock(
        test_layout_,
        {A_.transpose().solve(state_objective_derivative.block(0))});
    }

  private:
    DenseMatrix A_;
    DenseMatrix B_;
    DenseVector f_;
    DenseMatrix C_;
    DenseVector desired_observation_;
    DenseVector observation_weights_;
    DenseVector regularisation_weights_;
    double      alpha_;
    LayoutPtr   variable_layout_;
    LayoutPtr   test_layout_;
    LayoutPtr   state_layout_;
    LayoutPtr   control_layout_;
  };
} // namespace nmopt::reference
