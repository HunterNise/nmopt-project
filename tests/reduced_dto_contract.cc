#include "nmopt/contract/metric_constraint.hpp"
#include "nmopt/contract/reduced_dto.hpp"
#include "nmopt/experiment/reduced_envelope.hpp"
#include "nmopt/reference/linear_quadratic_model.hpp"
#include "nmopt/solvers/reduced_gradient.hpp"
#include "nmopt/solvers/reduced_line_search.hpp"
#include "nmopt/solvers/reduced_trust_region.hpp"
#include "test_support/contract_errors.hpp"
#include "test_support/scenario_dispatch.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace
{
  using namespace nmopt::contract;
  using nmopt::reference::LinearQuadraticModel;

  struct AlternateDenseBackend
  {
    using Vector = DenseVector;

    static Vector
    zeros(const std::size_t size)
    {
      return Vector(size);
    }

    static std::size_t
    size(const Vector &vector)
    {
      return vector.size();
    }

    static double
    dot(const Vector &left, const Vector &right)
    {
      return nmopt::contract::dot(left, right);
    }

    static void
    add_scaled(Vector &target, const double factor, const Vector &source)
    {
      target.add_scaled(factor, source);
    }

    static void
    scale(Vector &target, const double factor)
    {
      target.scale(factor);
    }
  };

  class AlternateQuadraticModel final
    : public ExecutableModelT<AlternateDenseBackend>
    , public ReducedHessianT<AlternateDenseBackend>
  {
  public:
    using Backend = AlternateDenseBackend;
    using Primal = PrimalBlockT<Backend>;
    using Covector = CovectorBlockT<Backend>;

    AlternateQuadraticModel()
      : variable_layout_(std::make_shared<const BlockLayout>(
          "alternate_variables",
          std::vector<SpaceId>{{"state"}, {"control"}},
          std::vector<std::size_t>{2, 2}))
      , test_layout_(std::make_shared<const BlockLayout>(
          "alternate_state_test",
          std::vector<SpaceId>{{"state_test"}},
          std::vector<std::size_t>{2}))
      , state_layout_(variable_layout_->single_block(0, "state"))
      , control_layout_(variable_layout_->single_block(1, "control"))
      , target_({1.0, -2.0})
    {}

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
    layout() const override
    {
      return control_layout_;
    }

    Covector
    residual(const Primal &variables) const override
    {
      require(variables.layout()->compatible_with(*variable_layout_),
              "Alternate residual received incompatible variables");
      DenseVector value = variables.block(0);
      value.add_scaled(-1.0, variables.block(1));
      return Covector(test_layout_, {std::move(value)});
    }

    Covector
    residual_jvp(const Primal &variables,
                 const Primal &variable_tangent) const override
    {
      require(variables.layout()->compatible_with(*variable_layout_) &&
                variable_tangent.layout()->compatible_with(*variable_layout_),
              "Alternate residual JVP received incompatible variables");
      DenseVector value = variable_tangent.block(0);
      value.add_scaled(-1.0, variable_tangent.block(1));
      return Covector(test_layout_, {std::move(value)});
    }

    Covector
    residual_vjp(const Primal &variables,
                 const PrimalBlockT<Backend> &test_seed) const override
    {
      require(variables.layout()->compatible_with(*variable_layout_) &&
                test_seed.layout()->compatible_with(*test_layout_),
              "Alternate residual VJP received incompatible variables");
      DenseVector control = test_seed.block(0);
      control.scale(-1.0);
      return Covector(variable_layout_, {test_seed.block(0), std::move(control)});
    }

    double
    objective(const Primal &variables) const override
    {
      require(variables.layout()->compatible_with(*variable_layout_),
              "Alternate objective received incompatible variables");
      DenseVector difference = variables.block(1);
      difference.add_scaled(-1.0, target_);
      return 0.5 * dot(difference, difference);
    }

    Covector
    objective_derivative(const Primal &variables) const override
    {
      require(variables.layout()->compatible_with(*variable_layout_),
              "Alternate objective derivative received incompatible variables");
      DenseVector control = variables.block(1);
      control.add_scaled(-1.0, target_);
      return Covector(variable_layout_, {DenseVector(2), std::move(control)});
    }

    Covector
    apply(const Primal &control, const Primal &direction) const override
    {
      require(control.layout()->compatible_with(*control_layout_) &&
                direction.layout()->compatible_with(*control_layout_),
              "Alternate Hessian received incompatible controls");
      return Covector(control_layout_, {direction.block(0)});
    }

    FormulationSolveResultT<Backend>
    solve_state(const Primal &control) const
    {
      require(control.layout()->compatible_with(*control_layout_),
              "Alternate state solve received incompatible control");
      return FormulationSolveResultT<Backend>(
        Primal(state_layout_, {control.block(0)}));
    }

    FormulationSolveResultT<Backend>
    solve_adjoint(const Primal &full_point,
                  const Covector &state_rhs) const
    {
      require(full_point.layout()->compatible_with(*variable_layout_) &&
                state_rhs.layout()->compatible_with(*state_layout_),
              "Alternate adjoint solve received incompatible arguments");
      return FormulationSolveResultT<Backend>(Primal::zeros(test_layout_));
    }

  private:
    LayoutPtr   variable_layout_;
    LayoutPtr   test_layout_;
    LayoutPtr   state_layout_;
    LayoutPtr   control_layout_;
    DenseVector target_;
  };

  class AlternateIdentityMetric final : public MetricT<AlternateDenseBackend>
  {
  public:
    explicit AlternateIdentityMetric(LayoutPtr layout)
      : layout_(std::move(layout))
    {}

    const std::string &
    id() const override
    {
      static const std::string id = "alternate_identity";
      return id;
    }

    const LayoutPtr &
    layout() const override
    {
      return layout_;
    }

    CovectorBlockT<AlternateDenseBackend>
    apply(const PrimalBlockT<AlternateDenseBackend> &primal) const override
    {
      require(primal.layout()->compatible_with(*layout_),
              "Alternate metric received incompatible primal");
      return CovectorBlockT<AlternateDenseBackend>(layout_, {primal.block(0)});
    }

    PrimalBlockT<AlternateDenseBackend>
    inverse_apply(
      const CovectorBlockT<AlternateDenseBackend> &covector) const override
    {
      require(covector.layout()->compatible_with(*layout_),
              "Alternate metric received incompatible covector");
      return PrimalBlockT<AlternateDenseBackend>(layout_, {covector.block(0)});
    }

  private:
    LayoutPtr layout_;
  };

  class CountingLinearQuadraticModel final : public ExecutableModel
  {
  public:
    explicit CountingLinearQuadraticModel(const LinearQuadraticModel &model)
      : model_(model)
    {}

    const LayoutPtr &
    variable_layout() const override
    {
      return model_.variable_layout();
    }

    const LayoutPtr &
    test_layout() const override
    {
      return model_.test_layout();
    }

    CovectorBlock
    residual(const PrimalBlock &variables) const override
    {
      return model_.residual(variables);
    }

    CovectorBlock
    residual_jvp(const PrimalBlock &variables,
                 const PrimalBlock &variable_tangent) const override
    {
      return model_.residual_jvp(variables, variable_tangent);
    }

    CovectorBlock
    residual_vjp(const PrimalBlock &variables,
                 const PrimalBlock &test_seed) const override
    {
      return model_.residual_vjp(variables, test_seed);
    }

    double
    objective(const PrimalBlock &variables) const override
    {
      ++objective_calls;
      return model_.objective(variables);
    }

    CovectorBlock
    objective_derivative(const PrimalBlock &variables) const override
    {
      ++objective_derivative_calls;
      return model_.objective_derivative(variables);
    }

    mutable std::size_t objective_calls = 0;
    mutable std::size_t objective_derivative_calls = 0;

  private:
    const LinearQuadraticModel &model_;
  };

  void
  require_close(const double actual,
                const double expected,
                const double tolerance,
                const std::string &description)
  {
    if (std::abs(actual - expected) > tolerance)
      throw ContractError(description + ": expected " +
                          std::to_string(expected) + ", got " +
                          std::to_string(actual));
  }

  PrimalBlock
  shifted(PrimalBlock value, const PrimalBlock &direction, const double step)
  {
    require_compatible(value, direction, "Shift has incompatible primal layouts");
    for (std::size_t block = 0; block < value.n_blocks(); ++block)
      value.add_scaled_block(block, step, direction.block(block));
    return value;
  }

  class RecordingConstraint final : public Constraint
  {
  public:
    explicit RecordingConstraint(const CellwiseBoxConstraint &constraint)
      : constraint_(constraint)
    {}

    const LayoutPtr &
    layout() const override
    {
      return constraint_.layout();
    }

    bool
    is_feasible(const PrimalBlock &primal) const override
    {
      return constraint_.is_feasible(primal);
    }

    bool
    supports_projection_in(const Metric &metric) const override
    {
      return constraint_.supports_projection_in(metric);
    }

    PrimalBlock
    project_in(const PrimalBlock &primal, const Metric &metric) const override
    {
      PrimalBlock projected = constraint_.project_in(primal, metric);
      projected_controls_.push_back(projected);
      return projected;
    }

    const std::vector<PrimalBlock> &
    projected_controls() const
    {
      return projected_controls_;
    }

  private:
    const CellwiseBoxConstraint &  constraint_;
    mutable std::vector<PrimalBlock> projected_controls_;
  };

  class NonDiagonalMetric final : public Metric
  {
  public:
    explicit NonDiagonalMetric(LayoutPtr layout)
      : layout_(std::move(layout))
    {}

    const std::string &
    id() const override
    {
      static const std::string id = "l2_cellwise";
      return id;
    }

    const LayoutPtr &
    layout() const override
    {
      return layout_;
    }

    CovectorBlock
    apply(const PrimalBlock &primal) const override
    {
      require(primal.layout()->compatible_with(*layout_),
              "Non-diagonal metric primal has an incompatible layout");
      const auto &values = primal.block(0);
      return CovectorBlock(
        layout_, {DenseVector{2.0 * values[0] + values[1],
                              values[0] + 2.0 * values[1]}});
    }

    PrimalBlock
    inverse_apply(const CovectorBlock &covector) const override
    {
      require(covector.layout()->compatible_with(*layout_),
              "Non-diagonal metric covector has an incompatible layout");
      const auto &values = covector.block(0);
      return PrimalBlock(
        layout_, {DenseVector{(2.0 * values[0] - values[1]) / 3.0,
                              (-values[0] + 2.0 * values[1]) / 3.0}});
    }

  private:
    LayoutPtr layout_;
  };

  class ScaledReducedHessian final : public ReducedHessian
  {
  public:
    ScaledReducedHessian(const ReducedHessian &hessian, const double scale)
      : hessian_(hessian)
      , scale_(scale)
    {
      require(scale_ > 0.0, "Scaled reduced Hessian factor must be positive");
    }

    const LayoutPtr &
    layout() const override
    {
      return hessian_.layout();
    }

    CovectorBlock
    apply(const PrimalBlock &control,
          const PrimalBlock &direction) const override
    {
      CovectorBlock result = hessian_.apply(control, direction);
      for (std::size_t block = 0; block < result.n_blocks(); ++block)
        result.scale_block(block, scale_);
      return result;
    }

  private:
    const ReducedHessian &hessian_;
    double                 scale_;
  };

  void
  test_block_layout_invariant()
  {
    using BlockAccess = decltype(std::declval<PrimalBlock &>().block(0));
    static_assert(std::is_same_v<BlockAccess, const DenseVector &>,
                  "BlockValues public block access must be read-only");

    const auto layout = std::make_shared<const BlockLayout>(
      "layout_invariant",
      std::vector<SpaceId>{{"state"}},
      std::vector<std::size_t>{2});
    PrimalBlock   primal(layout, {DenseVector{1.0, 2.0}});
    CovectorBlock covector(layout, {DenseVector{3.0, 4.0}});

    nmopt::test_support::require_contract_error(
      [&primal]() {
        primal.add_scaled_block(0, 1.0, DenseVector{1.0, 2.0, 3.0});
      },
      "BlockValues update vector dimension does not match layout",
      "dimension-changing primal block update");
    nmopt::test_support::require_contract_error(
      [&covector]() {
        covector.add_scaled_block(0, 1.0, DenseVector{1.0, 2.0, 3.0});
      },
      "BlockValues update vector dimension does not match layout",
      "dimension-changing covector block update");
    require_close(pair(covector, primal),
                  11.0,
                  1e-15,
                  "Rejected block updates preserve pairing");

    primal.add_scaled_block(0, 1.0, DenseVector{-1.0, 0.5});
    covector.scale_block(0, 0.5);
    require_close(pair(covector, primal),
                  5.0,
                  1e-15,
                  "Checked block algebra preserves pairing");
  }

  void
  test_v0_contract()
  {
    const LinearQuadraticModel model(
      DenseMatrix(2, 2, {4.0, -1.0, -1.0, 3.0}),
      DenseMatrix(2, 2, {1.0, 0.5, -0.25, 2.0}),
      DenseVector{1.0, -0.5},
      DenseMatrix(2, 2, {1.0, 0.0, 0.5, 1.0}),
      DenseVector{0.25, -1.0},
      DenseVector{1.5, 0.75},
      DenseVector{2.0, 3.0},
      0.4);

    const PrimalBlock point(
      model.variable_layout(),
      {DenseVector{0.2, -0.3}, DenseVector{0.5, -0.4}});
    const PrimalBlock tangent(
      model.variable_layout(),
      {DenseVector{-0.7, 0.25}, DenseVector{0.3, 0.8}});
    const PrimalBlock test_seed(model.test_layout(), {DenseVector{0.6, -1.1}});

    const CovectorBlock residual_direction = model.residual_jvp(point, tangent);
    const CovectorBlock residual_pullback = model.residual_vjp(point, test_seed);
    require_close(pair(residual_direction, test_seed),
                  pair(residual_pullback, tangent),
                  1e-13,
                  "Residual JVP/VJP pairing");

    constexpr double epsilon = 1e-7;
    CovectorBlock finite_difference =
      model.residual(shifted(point, tangent, epsilon));
    const CovectorBlock residual_at_point = model.residual(point);
    for (std::size_t block = 0; block < finite_difference.n_blocks(); ++block)
      {
        finite_difference.add_scaled_block(
          block, -1.0, residual_at_point.block(block));
        finite_difference.scale_block(block, 1.0 / epsilon);
        for (std::size_t entry = 0;
             entry < finite_difference.block(block).size();
             ++entry)
          require_close(finite_difference.block(block)[entry],
                        residual_direction.block(block)[entry],
                        1e-9,
                        "Residual finite-difference JVP");
      }

    const CovectorBlock objective_derivative =
      model.objective_derivative(point);
    const double objective_difference =
      model.objective(shifted(point, tangent, epsilon)) - model.objective(point);
    require_close(objective_difference / epsilon,
                  pair(objective_derivative, tangent),
                  2e-7,
                  "Objective directional derivative");

    nmopt::test_support::require_contract_error(
      [&model]() { (void)StateControlPartition(model, 0, 0); },
      "State and control blocks must be distinct",
      "invalid reduced partition");

    const StateControlPartition partition(model, 0, 1);
    const StateAdjointSolvers solvers{
      [&model](const PrimalBlock &control) { return model.solve_state(control); },
      [&model](const PrimalBlock &full_point, const CovectorBlock &state_rhs) {
        return model.solve_adjoint(full_point, state_rhs);
      }};
    nmopt::test_support::require_contract_error(
      [&model, &partition, &solvers]() {
        (void)ReducedDTO(
          model, partition, StateAdjointSolvers{solvers.solve_state, {}});
      },
      "Reduced DTO requires an adjoint solve operation",
      "missing reduced adjoint callback");
    const ReducedDTO reduced(model, partition, solvers);

    const PrimalBlock control(partition.control_layout(), {DenseVector{0.4, -0.3}});
    const ReducedEvaluation evaluation = reduced.evaluate(control);

    const CovectorBlock state_residual = model.residual(evaluation.full_point);
    require_close(dot(state_residual.block(0), state_residual.block(0)),
                  0.0,
                  1e-24,
                  "State solve residual");

    const PrimalBlock control_direction(
      partition.control_layout(), {DenseVector{-0.2, 0.7}});
    const double reduced_difference =
      reduced.evaluate(shifted(control, control_direction, epsilon)).objective_value -
      evaluation.objective_value;
    require_close(reduced_difference / epsilon,
                  pair(evaluation.reduced_derivative, control_direction),
                  2e-7,
                  "Reduced DTO derivative");

    nmopt::test_support::require_contract_error(
      [&partition]() {
        (void)DiagonalMetric(
          "invalid_metric",
          partition.control_layout(),
          {DenseVector{2.0, 0.0}});
      },
      "Metric diagonal must be strictly positive",
      "nonpositive metric diagonal");
    const DiagonalMetric metric(
      "l2_cellwise", partition.control_layout(), {DenseVector{2.0, 5.0}});
    const PrimalBlock direction =
      reduced.gradient_direction(evaluation.reduced_derivative, metric);
    const CovectorBlock reconstructed = metric.apply(direction);
    require_close(pair(reconstructed, control_direction),
                  pair(evaluation.reduced_derivative, control_direction),
                  1e-13,
                  "Metric inverse/apply consistency");
    const auto search_direction =
      nmopt::solvers::make_steepest_descent_direction(
        evaluation.reduced_derivative, metric);
    require(search_direction.directional_derivative < 0.0,
            "Steepest-descent search direction is not descending");
    require_close(search_direction.gradient_norm,
                  std::sqrt(pair(reconstructed, direction)),
                  1e-15,
                  "Steepest-descent direction metric norm");
    require_close(
      pair(evaluation.reduced_derivative, search_direction.direction),
      -pair(evaluation.reduced_derivative, direction),
      1e-15,
      "Steepest-descent direction sign");

    const ReducedHessian &hessian = model;
    const CovectorBlock hessian_action =
      hessian.apply(control, control_direction);
    CovectorBlock reduced_derivative_difference =
      reduced.evaluate(shifted(control, control_direction, epsilon))
        .reduced_derivative;
    for (std::size_t block = 0;
         block < reduced_derivative_difference.n_blocks();
         ++block)
      {
        reduced_derivative_difference.add_scaled_block(
          block, -1.0, evaluation.reduced_derivative.block(block));
        reduced_derivative_difference.scale_block(block, 1.0 / epsilon);
        for (std::size_t entry = 0;
             entry < reduced_derivative_difference.block(block).size();
             ++entry)
          require_close(reduced_derivative_difference.block(block)[entry],
                        hessian_action.block(block)[entry],
                        2e-7,
                        "Reduced Hessian finite-difference action");
      }
    const PrimalBlock second_hessian_direction(
      partition.control_layout(), {DenseVector{0.35, -0.2}});
    const CovectorBlock second_hessian_action =
      hessian.apply(control, second_hessian_direction);
    require_close(pair(hessian_action, second_hessian_direction),
                  pair(second_hessian_action, control_direction),
                  1e-13,
                  "Reduced Hessian symmetry");

    const auto build_trial =
      [&control, &search_direction](const double step) {
        return shifted(control, search_direction.direction, step);
      };
    std::size_t trial_state_calls   = 0;
    std::size_t trial_adjoint_calls = 0;
    const auto evaluate_trial_value =
      [&reduced, &trial_state_calls](const PrimalBlock &trial_control) {
        ++trial_state_calls;
        return reduced.evaluate_value(trial_control);
      };
    const auto augment_trial_derivative =
      [&reduced, &trial_adjoint_calls](const ReducedValueEvaluation &value) {
        ++trial_adjoint_calls;
        return reduced.augment_derivative(value);
      };

    nmopt::solvers::ArmijoLineSearchParameters armijo_parameters;
    armijo_parameters.maximum_trials = 20;
    armijo_parameters.initial_step_length = 10.0;
    const nmopt::solvers::ArmijoLineSearchPolicy armijo_policy(
      armijo_parameters);
    const auto armijo_result = armijo_policy.search(control,
                                                    evaluation,
                                                    search_direction,
                                                    build_trial,
                                                    evaluate_trial_value,
                                                    augment_trial_derivative);
    require(armijo_result.accepted() && armijo_result.trial_count > 0,
            "Armijo line search did not accept a trial");
    require(armijo_result.evaluation.objective_value <=
              evaluation.objective_value,
            "Armijo line search accepted an objective increase");
    require(trial_state_calls == armijo_result.trial_count &&
              trial_adjoint_calls == 1,
            "Armijo line search did not defer derivative work until acceptance");
    require(armijo_result.acceptance_evidence.policy_name == "armijo" &&
              std::isfinite(armijo_result.acceptance_evidence.sufficient_decrease_bound) &&
              armijo_result.acceptance_evidence.trial_slope < 0.0,
            "Armijo line search did not report its acceptance evidence");

    const nmopt::solvers::ExactQuadraticLineSearchPolicy exact_policy(hessian);
    const std::size_t exact_state_calls = trial_state_calls;
    const std::size_t exact_adjoint_calls = trial_adjoint_calls;
    const auto exact_result = exact_policy.search(control,
                                                  evaluation,
                                                  search_direction,
                                                  build_trial,
                                                  evaluate_trial_value,
                                                  augment_trial_derivative);
    require(exact_result.accepted() && exact_result.trial_count == 1 &&
              exact_result.hessian_action_count == 1,
            "Exact quadratic line search did not accept its Hessian step");
    require(exact_result.evaluation.objective_value <=
              evaluation.objective_value,
            "Exact quadratic line search accepted an objective increase");
    require(trial_state_calls == exact_state_calls + 1 &&
              trial_adjoint_calls == exact_adjoint_calls + 1,
            "Exact quadratic line search did not augment its accepted value once");
    require(exact_result.acceptance_evidence.policy_name == "exact_quadratic" &&
              std::isfinite(exact_result.acceptance_evidence.curvature_value) &&
              exact_result.acceptance_evidence.curvature_value > 0.0,
            "Exact quadratic line search did not report its curvature evidence");
    PrimalBlock exact_update = exact_result.control;
    nmopt::solvers::add_scaled_primal(exact_update, -1.0, control);
    PrimalBlock expected_exact_update = search_direction.direction;
    nmopt::solvers::scale_primal(expected_exact_update,
                                 exact_result.step_length);
    for (std::size_t block = 0; block < exact_update.n_blocks(); ++block)
      for (std::size_t entry = 0;
           entry < exact_update.block(block).size();
           ++entry)
        require_close(exact_update.block(block)[entry],
                      expected_exact_update.block(block)[entry],
                      1e-12,
                      "Exact quadratic line search reported the wrong displacement");
    require_close(pair(exact_result.evaluation.reduced_derivative,
                       search_direction.direction),
                  0.0,
                  1e-12,
                  "Exact quadratic line search did not reach the one-dimensional stationary point");

    nmopt::solvers::WolfeLineSearchParameters wolfe_parameters;
    wolfe_parameters.maximum_trials = 30;
    wolfe_parameters.initial_step_length = 10.0;
    const nmopt::solvers::WolfeLineSearchPolicy wolfe_policy(
      wolfe_parameters);
    const std::size_t wolfe_state_calls = trial_state_calls;
    const std::size_t wolfe_adjoint_calls = trial_adjoint_calls;
    const auto wolfe_result = wolfe_policy.search(control,
                                                  evaluation,
                                                  search_direction,
                                                  build_trial,
                                                  evaluate_trial_value,
                                                  augment_trial_derivative);
    require(wolfe_result.accepted() && wolfe_result.trial_count > 0,
            "Wolfe line search did not accept a trial");
    require(trial_state_calls == wolfe_state_calls + wolfe_result.trial_count &&
              trial_adjoint_calls ==
                wolfe_adjoint_calls + wolfe_result.trial_count,
            "Wolfe line search did not retain derivative work for every trial");
    require(wolfe_result.acceptance_evidence.policy_name == "strong_wolfe" &&
              std::isfinite(wolfe_result.acceptance_evidence.sufficient_decrease_bound) &&
              std::isfinite(wolfe_result.acceptance_evidence.trial_slope),
            "Wolfe line search did not report its acceptance evidence");

    nmopt::solvers::WeakWolfeLineSearchParameters weak_wolfe_parameters;
    weak_wolfe_parameters.maximum_trials = 30;
    weak_wolfe_parameters.initial_step_length = 10.0;
    const nmopt::solvers::WeakWolfeLineSearchPolicy weak_wolfe_policy(
      weak_wolfe_parameters);
    const std::size_t weak_wolfe_state_calls = trial_state_calls;
    const std::size_t weak_wolfe_adjoint_calls = trial_adjoint_calls;
    const auto weak_wolfe_result = weak_wolfe_policy.search(control,
                                                            evaluation,
                                                            search_direction,
                                                            build_trial,
                                                            evaluate_trial_value,
                                                            augment_trial_derivative);
    require(weak_wolfe_result.accepted() && weak_wolfe_result.trial_count > 0,
            "Weak Wolfe line search did not accept a trial");
    require(weak_wolfe_result.evaluation.objective_value <=
              evaluation.objective_value,
            "Weak Wolfe line search accepted an objective increase");
    require(trial_state_calls ==
              weak_wolfe_state_calls + weak_wolfe_result.trial_count &&
              trial_adjoint_calls ==
                weak_wolfe_adjoint_calls + weak_wolfe_result.trial_count,
            "Weak Wolfe line search did not retain derivative work for every trial");
    require(weak_wolfe_result.acceptance_evidence.policy_name == "weak_wolfe" &&
              std::isfinite(weak_wolfe_result.acceptance_evidence.sufficient_decrease_bound) &&
              std::isfinite(weak_wolfe_result.acceptance_evidence.trial_slope),
            "Weak Wolfe line search did not report its acceptance evidence");
    PrimalBlock weak_wolfe_update = weak_wolfe_result.control;
    nmopt::solvers::add_scaled_primal(weak_wolfe_update, -1.0, control);
    const double weak_wolfe_initial_slope =
      pair(evaluation.reduced_derivative, weak_wolfe_update);
    const double weak_wolfe_trial_slope =
      pair(weak_wolfe_result.evaluation.reduced_derivative,
           weak_wolfe_update);
    require(weak_wolfe_initial_slope < 0.0 &&
              weak_wolfe_trial_slope >=
                weak_wolfe_parameters.curvature_fraction *
                  weak_wolfe_initial_slope,
            "Weak Wolfe line search did not satisfy its one-sided curvature condition");

    nmopt::solvers::ArmijoLineSearchParameters actual_displacement_parameters;
    actual_displacement_parameters.maximum_trials = 1;
    const nmopt::solvers::ArmijoLineSearchPolicy actual_displacement_policy(
      actual_displacement_parameters);
    const auto scaled_trial_builder =
      [&control, &search_direction](const double step) {
        return shifted(control, search_direction.direction, 0.25 * step);
      };
    const auto actual_displacement_result =
      actual_displacement_policy.search(control,
                                        evaluation,
                                        search_direction,
                                        scaled_trial_builder,
                                        evaluate_trial_value,
                                        augment_trial_derivative);
    require(actual_displacement_result.accepted() &&
              actual_displacement_result.trial_count == 1,
            "Armijo line search did not use the actual trial displacement");

    nmopt::test_support::require_contract_error(
      [&exact_policy, &control, &evaluation, &search_direction,
       &scaled_trial_builder, &evaluate_trial_value,
       &augment_trial_derivative]() {
        (void)exact_policy.search(control,
                                  evaluation,
                                  search_direction,
                                  scaled_trial_builder,
                                  evaluate_trial_value,
                                  augment_trial_derivative);
      },
      "Exact quadratic line search requires the trial builder to return the straight-line step",
      "exact line search rejects a transformed trial");

    const nmopt::solvers::ExactQuadraticLineSearchPolicy
      missing_exact_hessian_policy;
    nmopt::test_support::require_contract_error(
      [&missing_exact_hessian_policy, &control, &evaluation, &search_direction,
       &build_trial, &evaluate_trial_value, &augment_trial_derivative]() {
        (void)missing_exact_hessian_policy.search(control,
                                                  evaluation,
                                                  search_direction,
                                                  build_trial,
                                                  evaluate_trial_value,
                                                  augment_trial_derivative);
      },
      "Exact quadratic line search requires a reduced Hessian capability",
      "exact line search missing Hessian capability");

    const DiagonalMetric cg_metric(
      "cg_metric", partition.control_layout(), {DenseVector{1.0, 1.0}});
    using CGPolicy =
      nmopt::solvers::NonlinearConjugateGradientDirectionPolicyDense;
    CGPolicy cg_policy({2, 1e-14});
    const CovectorBlock first_cg_derivative(
      partition.control_layout(), {DenseVector{1.0, 0.0}});
    const CovectorBlock second_cg_derivative(
      partition.control_layout(), {DenseVector{1.0, 1.0}});
    const CovectorBlock third_cg_derivative(
      partition.control_layout(), {DenseVector{0.5, 0.5}});
    const auto first_cg_direction =
      cg_policy.next(control, first_cg_derivative, cg_metric);
    const auto second_cg_direction =
      cg_policy.next(control, second_cg_derivative, cg_metric);
    require_close(first_cg_direction.direction.block(0)[0],
                  -1.0,
                  1e-15,
                  "Nonlinear CG initial direction");
    require_close(second_cg_direction.direction.block(0)[0],
                  -2.0,
                  1e-15,
                  "Nonlinear CG Polak-Ribiere direction");
    require_close(second_cg_direction.direction.block(0)[1],
                  -1.0,
                  1e-15,
                  "Nonlinear CG Polak-Ribiere update");
    require(second_cg_direction.directional_derivative < 0.0,
            "Nonlinear CG Polak-Ribiere direction is not descending");
    (void)cg_policy.next(control, third_cg_derivative, cg_metric);
    require(cg_policy.restart_count() == 1,
            "Nonlinear CG did not perform its configured periodic restart");

    using FRPolicy =
      nmopt::solvers::FletcherReevesDirectionPolicyDense;
    FRPolicy fr_policy({10, 1e-14});
    const auto first_fr_direction =
      fr_policy.next(control, first_cg_derivative, cg_metric);
    const auto second_fr_direction =
      fr_policy.next(control, second_cg_derivative, cg_metric);
    require_close(first_fr_direction.direction.block(0)[0],
                  -1.0,
                  1e-15,
                  "Fletcher-Reeves initial direction");
    require_close(second_fr_direction.direction.block(0)[0],
                  -3.0,
                  1e-15,
                  "Fletcher-Reeves coefficient");
    require_close(second_fr_direction.direction.block(0)[1],
                  -1.0,
                  1e-15,
                  "Fletcher-Reeves direction update");
    require(second_fr_direction.directional_derivative < 0.0,
            "Fletcher-Reeves direction is not descending");

    using QCGPolicy =
      nmopt::solvers::QuadraticConjugateGradientDirectionPolicyDense;
    QCGPolicy qcg_policy({2, 1e-14});
    const auto first_qcg_direction =
      qcg_policy.next(control, first_cg_derivative, cg_metric);
    const auto second_qcg_direction =
      qcg_policy.next(control, second_cg_derivative, cg_metric);
    require_close(first_qcg_direction.direction.block(0)[0],
                  -1.0,
                  1e-15,
                  "Quadratic CG initial direction");
    require_close(second_qcg_direction.direction.block(0)[0],
                  -3.0,
                  1e-15,
                  "Quadratic CG coefficient");
    require_close(second_qcg_direction.direction.block(0)[1],
                  -1.0,
                  1e-15,
                  "Quadratic CG direction update");
    require(second_qcg_direction.directional_derivative < 0.0,
            "Quadratic CG direction is not descending");
    (void)qcg_policy.next(control, third_cg_derivative, cg_metric);
    require(qcg_policy.restart_count() == 1,
            "Quadratic CG did not perform its configured periodic restart");
    QCGPolicy invalid_qcg_policy({10, 1e-14});
    (void)invalid_qcg_policy.next(control, first_cg_derivative, cg_metric);
    nmopt::test_support::require_contract_error(
      [&invalid_qcg_policy, &control, &cg_metric]() {
        (void)invalid_qcg_policy.next(
          control,
          CovectorBlock(control.layout(), {DenseVector{-1.0, 0.0}}),
          cg_metric);
      },
      "Classical quadratic CG could not produce a descent direction",
      "quadratic CG rejects a non-descent recurrence");
    CGPolicy negative_beta_policy({10, 1e-14});
    (void)negative_beta_policy.next(control, first_cg_derivative, cg_metric);
    (void)negative_beta_policy.next(
      control,
      CovectorBlock(partition.control_layout(), {DenseVector{0.5, 0.0}}),
      cg_metric);
    require(negative_beta_policy.restart_count() == 1,
            "Nonlinear CG did not restart after a nonpositive coefficient");
    nmopt::test_support::require_contract_error(
      [&cg_policy, &control]() {
        const auto incompatible_layout = std::make_shared<const BlockLayout>(
          "incompatible_cg_layout",
          std::vector<SpaceId>{{"control"}},
          std::vector<std::size_t>{3});
        const DiagonalMetric incompatible_metric(
          "incompatible_cg_metric",
          incompatible_layout,
          {DenseVector{1.0, 1.0, 1.0}});
        const CovectorBlock incompatible_derivative(
          incompatible_layout,
          {DenseVector{1.0, 1.0, 1.0}});
        (void)cg_policy.next(control, incompatible_derivative, incompatible_metric);
      },
      "Nonlinear CG derivative history has an incompatible layout",
      "nonlinear CG history layout mismatch");

    using BfgsPolicy =
      nmopt::solvers::LimitedMemoryBfgsDirectionPolicyDense;
    BfgsPolicy bfgs_policy({2, 1e-14});
    const PrimalBlock bfgs_control_0(
      partition.control_layout(), {DenseVector{0.0, 0.0}});
    const PrimalBlock bfgs_control_1(
      partition.control_layout(), {DenseVector{1.0, 0.0}});
    const CovectorBlock bfgs_derivative_0(
      partition.control_layout(), {DenseVector{1.0, 0.0}});
    const CovectorBlock bfgs_derivative_1(
      partition.control_layout(), {DenseVector{2.0, 0.0}});
    (void)bfgs_policy.next(bfgs_control_0, bfgs_derivative_0, cg_metric);
    const auto bfgs_direction =
      bfgs_policy.next(bfgs_control_1, bfgs_derivative_1, cg_metric);
    require(bfgs_policy.history_size() == 1,
            "L-BFGS did not retain a valid secant pair");
    require(bfgs_policy.last_update_status() ==
              nmopt::solvers::LimitedMemoryBfgsUpdateStatus::accepted_pair,
            "L-BFGS did not report an accepted secant pair");
    require(bfgs_direction.directional_derivative < 0.0,
            "L-BFGS direction is not descending");
    (void)bfgs_policy.next(bfgs_control_1, bfgs_derivative_1, cg_metric);
    require(bfgs_policy.history_size() == 0 &&
              bfgs_policy.reset_count() == 1,
            "L-BFGS did not reset after failed curvature");
    require(bfgs_policy.last_update_status() ==
              nmopt::solvers::LimitedMemoryBfgsUpdateStatus::curvature_reset,
            "L-BFGS curvature reset was not reported");

    using FullBfgsPolicy = nmopt::solvers::FullBfgsDirectionPolicyDense;
    FullBfgsPolicy full_bfgs_policy({1e-14});
    const PrimalBlock full_bfgs_control_2(
      partition.control_layout(), {DenseVector{1.0, 1.0}});
    const CovectorBlock full_bfgs_derivative_2(
      partition.control_layout(), {DenseVector{2.0, 2.0}});
    (void)full_bfgs_policy.next(
      bfgs_control_0, bfgs_derivative_0, cg_metric);
    (void)full_bfgs_policy.next(
      bfgs_control_1, bfgs_derivative_1, cg_metric);
    const auto full_bfgs_direction = full_bfgs_policy.next(
      full_bfgs_control_2, full_bfgs_derivative_2, cg_metric);
    require(full_bfgs_policy.history_size() == 2,
            "full BFGS did not retain all valid secant pairs");
    require(full_bfgs_policy.last_update_status() ==
              nmopt::solvers::FullBfgsUpdateStatus::accepted_pair,
            "full BFGS did not report an accepted secant pair");
    require(full_bfgs_direction.directional_derivative < 0.0,
            "full BFGS direction is not descending");
    (void)full_bfgs_policy.next(
      full_bfgs_control_2, full_bfgs_derivative_2, cg_metric);
    require(full_bfgs_policy.history_size() == 0 &&
              full_bfgs_policy.reset_count() == 1,
            "full BFGS did not reset after failed curvature");
    require(full_bfgs_policy.last_update_status() ==
              nmopt::solvers::FullBfgsUpdateStatus::curvature_reset,
            "full BFGS curvature reset was not reported");

    nmopt::solvers::ReducedNewtonParameters newton_parameters;
    newton_parameters.maximum_inner_iterations = 10;
    newton_parameters.relative_tolerance = 1e-12;
    newton_parameters.absolute_tolerance = 1e-14;
    newton_parameters.curvature_tolerance = 1e-14;
    nmopt::solvers::ReducedSolverParameters newton_solver_parameters;
    newton_solver_parameters.gradient_tolerance = 1e-8;
    newton_solver_parameters.initial_step_length = 10.0;
    const nmopt::solvers::NewtonDirectionPolicyDense newton_policy(
      hessian, newton_parameters);
    const nmopt::solvers::ReducedNewtonSolver newton_solver(
      reduced, metric, newton_solver_parameters, newton_policy);
    const auto newton_result = newton_solver.solve(control);
    require(newton_result.stopping_reason ==
              nmopt::solvers::ReducedGradientStoppingReason::gradient_tolerance,
            "Reduced Newton solver did not reach its tolerance");
    require(newton_result.hessian_action_count > 0,
            "Reduced Newton solver did not report Hessian actions");
    require(newton_result.step_length_history.size() ==
              newton_result.accepted_iterations,
            "Reduced Newton step history does not match accepted iterations");
    require(newton_result.metric_solve_count >
              newton_result.hessian_action_count,
            "Reduced Newton did not report metric preconditioning actions");
    require(newton_result.hessian_solve_history.size() ==
              newton_result.gradient_norm_history.size(),
            "Reduced Newton PCG history is not aligned with direction evaluations");
    std::size_t reported_hessian_actions = 0;
    bool saw_inner_pcg = false;
    for (const auto &diagnostics : newton_result.hessian_solve_history)
      {
        require(diagnostics.initial_residual_norm >= 0.0 &&
                  std::isfinite(diagnostics.initial_residual_norm) &&
                  diagnostics.final_residual_norm >= 0.0 &&
                  std::isfinite(diagnostics.final_residual_norm),
                "Reduced Newton PCG residual diagnostics are invalid");
        require(diagnostics.iteration_count <=
                  newton_parameters.maximum_inner_iterations,
                "Reduced Newton PCG iteration count exceeds its limit");
        const double target_residual = std::max(
          newton_parameters.absolute_tolerance,
          newton_parameters.relative_tolerance *
            diagnostics.initial_residual_norm);
        require(diagnostics.final_residual_norm <=
                  target_residual + 1e-12,
                "Reduced Newton PCG final residual exceeds its target");
        reported_hessian_actions += diagnostics.iteration_count;
        saw_inner_pcg = saw_inner_pcg || diagnostics.iteration_count > 0;
      }
    require(saw_inner_pcg,
            "Reduced Newton did not expose a nontrivial PCG solve");
    require(reported_hessian_actions == newton_result.hessian_action_count,
            "Reduced Newton PCG diagnostics do not match Hessian actions");
    require(newton_result.metric_solve_count ==
              2 * newton_result.gradient_norm_history.size() +
                newton_result.hessian_action_count,
            "Reduced Newton PCG diagnostics do not match metric actions");

    nmopt::solvers::ReducedNewtonParameters loose_newton_parameters =
      newton_parameters;
    loose_newton_parameters.absolute_tolerance = 1e6;
    loose_newton_parameters.relative_tolerance = 0.9;
    nmopt::solvers::NewtonDirectionPolicyDense loose_newton_policy(
      hessian, loose_newton_parameters);
    const auto loose_newton_direction = loose_newton_policy.next(
      control, evaluation.reduced_derivative, metric);
    const double loose_newton_target = std::max(
      loose_newton_parameters.absolute_tolerance,
      loose_newton_parameters.relative_tolerance *
        loose_newton_direction.hessian_solve.initial_residual_norm);
    require(loose_newton_direction.directional_derivative < 0.0 &&
              loose_newton_direction.hessian_solve.iteration_count > 0 &&
              loose_newton_direction.hessian_solve.final_residual_norm <=
                loose_newton_target + 1e-12,
            "Newton loose inner tolerance returned an unusable zero direction");

    const nmopt::solvers::NewtonDirectionPolicyDense missing_hessian_policy;
    const nmopt::solvers::ReducedNewtonSolver missing_hessian_solver(
      reduced, metric, newton_solver_parameters, missing_hessian_policy);
    nmopt::test_support::require_contract_error(
      [&missing_hessian_solver, &control]() {
        (void)missing_hessian_solver.solve(control);
      },
      "Newton direction requires a reduced Hessian capability",
      "Newton missing Hessian capability");

    nmopt::test_support::require_contract_error(
      [&partition, &metric]() {
        (void)CellwiseBoxConstraint(
          partition.control_layout(),
          {DenseVector{0.0, 0.6}},
          {DenseVector{1.0, 0.5}},
          metric);
      },
      "Cellwise box lower bound exceeds upper bound",
      "reversed cellwise bounds");
    const CellwiseBoxConstraint bounds(
      partition.control_layout(), {DenseVector{-0.1, -0.2}},
      {DenseVector{0.6, 0.5}}, metric);
    const PrimalBlock projected = bounds.project_in(
      PrimalBlock(partition.control_layout(), {DenseVector{-0.4, 0.8}}), metric);
    require(bounds.is_feasible(projected),
            "Cellwise L2 projection must return a feasible control");
    require_close(projected.block(0)[0], -0.1, 1e-15,
                  "Cellwise box lower projection");
    require_close(projected.block(0)[1], 0.5, 1e-15,
                  "Cellwise box upper projection");

    nmopt::solvers::ReducedGradientParameters solver_parameters;
    solver_parameters.maximum_iterations = 100;
    solver_parameters.maximum_line_search_trials = 20;
    solver_parameters.gradient_tolerance = 1e-8;
    solver_parameters.initial_step_length = 10.0;
    solver_parameters.armijo_fraction = 1e-4;
    solver_parameters.backtracking_factor = 0.5;
    const nmopt::solvers::ReducedGradientSolver solver(
      reduced, metric, solver_parameters);
    const auto solver_result = solver.solve(
      PrimalBlock(partition.control_layout(), {DenseVector{1.0, -1.0}}));

    require(solver_result.stopping_reason ==
              nmopt::solvers::ReducedGradientStoppingReason::gradient_tolerance,
            "Dense reduced gradient solver did not reach its tolerance");
    require(solver_result.objective_history.size() > 1,
            "Dense reduced gradient solver did not accept an iteration");
    for (std::size_t index = 1;
         index < solver_result.objective_history.size();
         ++index)
      require(solver_result.objective_history[index] <=
                solver_result.objective_history[index - 1],
              "Dense reduced gradient objective history is not monotonic");
    require(solver_result.gradient_norm_history.back() <=
              solver_parameters.gradient_tolerance,
            "Dense reduced gradient final norm exceeds the configured tolerance");
    require(solver_result.line_search_trial_count + 1 ==
              solver_result.state_solve_count,
            "Dense reduced gradient state count misses a trial evaluation");
    require(solver_result.accepted_iterations + 1 ==
              solver_result.adjoint_solve_count,
            "Dense reduced gradient adjoint count includes rejected trials");
    require(solver_result.state_solve_count >
              solver_result.objective_history.size(),
            "Dense reduced gradient test did not exercise Armijo backtracking");
    require(solver_result.metric_solve_count ==
              solver_result.gradient_norm_history.size(),
            "Dense reduced gradient metric solve count does not match direction evaluations");
    require(solver_result.step_length_history.size() ==
              solver_result.accepted_iterations,
            "Dense reduced gradient step history does not match accepted iterations");
    require(solver_result.objective_change_history.size() ==
              solver_result.accepted_iterations,
            "Dense reduced gradient objective-change history does not match accepted iterations");
    require(solver_result.relative_gradient_norm_history.size() ==
              solver_result.gradient_norm_history.size(),
            "Dense reduced gradient relative-norm history does not match direction evaluations");
    require(solver_result.step_norm_history.size() ==
              solver_result.accepted_iterations,
            "Dense reduced gradient step-norm history does not match accepted iterations");
    require(solver_result.final_evaluation.state_solve.converged() &&
              solver_result.final_evaluation.adjoint_solve.converged(),
            "Dense reduced gradient result does not retain the final evaluation");
    require_close(solver_result.final_evaluation.objective_value,
                  solver_result.objective_history.back(),
                  1e-14,
                  "Dense reduced gradient final objective");
    require(solver_result.parameters.maximum_iterations ==
              solver_parameters.maximum_iterations &&
              solver_result.policy_name == "armijo" &&
              solver_result.policy_parameters.policy_name == "armijo" &&
              solver_result.policy_parameters.maximum_trials ==
                solver_parameters.maximum_line_search_trials &&
              solver_result.policy_parameters.armijo_fraction ==
                solver_parameters.armijo_fraction &&
              solver_result.iteration_records.size() ==
                solver_result.accepted_iterations,
            "Dense reduced gradient result does not retain its audit snapshot");
    require(solver_result.final_evaluation.reduced_derivative.layout()->compatible_with(
              *solver_result.control.layout()),
            "Dense reduced gradient final covector has the wrong layout");
    for (std::size_t index = 0;
         index < solver_result.accepted_iterations;
         ++index)
      {
        require(solver_result.step_length_history[index] > 0.0,
                "Dense reduced gradient accepted a nonpositive step");
        require(solver_result.step_norm_history[index] > 0.0,
                "Dense reduced gradient accepted a zero metric step");
        require_close(solver_result.objective_change_history[index],
                      solver_result.objective_history[index + 1] -
                        solver_result.objective_history[index],
                      1e-14,
                      "Dense reduced gradient objective-change history");
        const auto &record = solver_result.iteration_records[index];
        require(record.common.iteration == index + 1 &&
                  record.common.policy_name == "armijo" &&
                  record.common.trial_count > 0 &&
                  record.common.actual_step_norm == solver_result.step_norm_history[index] &&
                  record.common.per_iteration_work.state_solve_count ==
                    record.common.trial_count &&
                  record.common.per_iteration_work.adjoint_solve_count == 1 &&
                  record.common.per_iteration_work.metric_solve_count == 1 &&
                  record.acceptance_evidence.trial_slope < 0.0 &&
                  std::isfinite(record.acceptance_evidence.sufficient_decrease_bound) &&
                  record.common.accepted_evaluation.objective_value ==
                    solver_result.objective_history[index + 1],
                "Dense reduced gradient accepted-iteration audit is incomplete");
      }

    nmopt::solvers::ReducedGradientParameters relative_stopping_parameters =
      solver_parameters;
    relative_stopping_parameters.gradient_tolerance = 1e-30;
    relative_stopping_parameters.relative_gradient_tolerance = 0.5;
    relative_stopping_parameters.stopping_criterion =
      nmopt::solvers::ReducedStoppingCriterion::relative_gradient_norm;
    relative_stopping_parameters.gradient_tolerance = 1e6;
    const nmopt::solvers::ReducedGradientSolver relative_stopping_solver(
      reduced, metric, relative_stopping_parameters);
    const auto relative_stopping_result = relative_stopping_solver.solve(
      PrimalBlock(partition.control_layout(), {DenseVector{1.0, -1.0}}));
    require(relative_stopping_result.stopping_reason ==
              nmopt::solvers::ReducedGradientStoppingReason::relative_gradient_tolerance,
            "Reduced solver did not stop on the relative gradient tolerance");
    require(relative_stopping_result.relative_gradient_norm_history.back() <=
              relative_stopping_parameters.relative_gradient_tolerance,
            "Reduced solver relative gradient tolerance was not enforced");

    nmopt::solvers::ReducedGradientParameters objective_stopping_parameters =
      solver_parameters;
    objective_stopping_parameters.gradient_tolerance = 1e-30;
    objective_stopping_parameters.objective_change_tolerance = 1e6;
    objective_stopping_parameters.stopping_criterion =
      nmopt::solvers::ReducedStoppingCriterion::objective_change;
    objective_stopping_parameters.gradient_tolerance = 1e6;
    const nmopt::solvers::ReducedGradientSolver objective_stopping_solver(
      reduced, metric, objective_stopping_parameters);
    const auto objective_stopping_result = objective_stopping_solver.solve(
      PrimalBlock(partition.control_layout(), {DenseVector{1.0, -1.0}}));
    require(objective_stopping_result.stopping_reason ==
              nmopt::solvers::ReducedGradientStoppingReason::objective_change_tolerance,
            "Reduced solver did not stop on the objective-change tolerance");
    require(objective_stopping_result.accepted_iterations > 0,
            "Objective-change stopping did not retain an accepted trial");

    nmopt::solvers::ReducedGradientParameters step_stopping_parameters =
      solver_parameters;
    step_stopping_parameters.gradient_tolerance = 1e-30;
    step_stopping_parameters.step_tolerance = 1e6;
    step_stopping_parameters.stopping_criterion =
      nmopt::solvers::ReducedStoppingCriterion::step_norm;
    step_stopping_parameters.gradient_tolerance = 1e6;
    const nmopt::solvers::ReducedGradientSolver step_stopping_solver(
      reduced, metric, step_stopping_parameters);
    const auto step_stopping_result = step_stopping_solver.solve(
      PrimalBlock(partition.control_layout(), {DenseVector{1.0, -1.0}}));
    require(step_stopping_result.stopping_reason ==
              nmopt::solvers::ReducedGradientStoppingReason::step_tolerance,
            "Reduced solver did not stop on the step tolerance");
    require(step_stopping_result.step_norm_history.back() <=
              step_stopping_parameters.step_tolerance,
            "Reduced solver step tolerance was not enforced");

    const nmopt::solvers::ReducedConjugateGradientSolver cg_solver(
      reduced, metric, solver_parameters);
    const auto cg_result = cg_solver.solve(
      PrimalBlock(partition.control_layout(), {DenseVector{1.0, -1.0}}));
    require(cg_result.stopping_reason ==
              nmopt::solvers::ReducedGradientStoppingReason::gradient_tolerance,
            "Dense nonlinear CG solver did not reach its tolerance");
    require(cg_result.metric_solve_count ==
              cg_result.gradient_norm_history.size(),
            "Dense nonlinear CG metric solve count does not match direction evaluations");
    require(cg_result.step_length_history.size() ==
              cg_result.accepted_iterations,
            "Dense nonlinear CG step history does not match accepted iterations");
    for (std::size_t index = 1;
         index < cg_result.objective_history.size();
         ++index)
      require(cg_result.objective_history[index] <=
                cg_result.objective_history[index - 1],
              "Dense nonlinear CG objective history is not monotonic");

    const nmopt::solvers::ReducedFletcherReevesSolver fr_solver(
      reduced, metric, solver_parameters);
    const auto fr_result = fr_solver.solve(
      PrimalBlock(partition.control_layout(), {DenseVector{1.0, -1.0}}));
    require(fr_result.stopping_reason ==
              nmopt::solvers::ReducedGradientStoppingReason::gradient_tolerance,
            "Dense Fletcher-Reeves solver did not reach its tolerance");
    require(fr_result.metric_solve_count ==
              fr_result.gradient_norm_history.size(),
            "Dense Fletcher-Reeves metric solve count does not match direction evaluations");
    for (std::size_t index = 1;
         index < fr_result.objective_history.size();
         ++index)
      require(fr_result.objective_history[index] <=
                fr_result.objective_history[index - 1],
              "Dense Fletcher-Reeves objective history is not monotonic");

    nmopt::solvers::ReducedGradientParameters exact_cg_parameters =
      solver_parameters;
    exact_cg_parameters.gradient_tolerance = 1e-12;
    const nmopt::solvers::ExactQuadraticLineSearchPolicy exact_cg_line_search(
      hessian);
    const nmopt::solvers::ReducedExactConjugateGradientSolver exact_pr_solver(
      reduced,
      metric,
      exact_cg_parameters,
      CGPolicy{},
      exact_cg_line_search);
    const auto exact_pr_result = exact_pr_solver.solve(
      PrimalBlock(partition.control_layout(), {DenseVector{1.0, -1.0}}));
    const nmopt::solvers::ReducedExactFletcherReevesSolver exact_fr_solver(
      reduced,
      metric,
      exact_cg_parameters,
      FRPolicy{},
      exact_cg_line_search);
    const auto exact_fr_result = exact_fr_solver.solve(
      PrimalBlock(partition.control_layout(), {DenseVector{1.0, -1.0}}));
    require(exact_pr_result.stopping_reason ==
              nmopt::solvers::ReducedGradientStoppingReason::gradient_tolerance &&
              exact_fr_result.stopping_reason ==
                nmopt::solvers::ReducedGradientStoppingReason::gradient_tolerance,
            "Exact-search PR+ and Fletcher-Reeves did not reach tolerance");
    require(exact_pr_result.accepted_iterations > 0 &&
              exact_pr_result.accepted_iterations ==
                exact_fr_result.accepted_iterations,
            "Exact-search PR+ and Fletcher-Reeves used different iteration counts");
    require(exact_pr_result.line_search_trial_count ==
              exact_pr_result.accepted_iterations &&
              exact_fr_result.line_search_trial_count ==
                exact_fr_result.accepted_iterations,
            "Exact quadratic CG searches did not use one trial per iteration");
    require(exact_pr_result.objective_history.size() ==
              exact_fr_result.objective_history.size(),
            "Exact-search PR+ and Fletcher-Reeves histories have different sizes");
    for (std::size_t index = 0;
         index < exact_pr_result.objective_history.size();
         ++index)
      require_close(exact_pr_result.objective_history[index],
                    exact_fr_result.objective_history[index],
                    1e-12,
                    "Exact-search PR+ and Fletcher-Reeves objective equivalence");
    for (std::size_t entry = 0;
         entry < exact_pr_result.control.block(0).size();
         ++entry)
      require_close(exact_pr_result.control.block(0)[entry],
                    exact_fr_result.control.block(0)[entry],
                    1e-10,
                    "Exact-search PR+ and Fletcher-Reeves control equivalence");

    const nmopt::solvers::ReducedQuadraticConjugateGradientSolver qcg_solver(
      reduced,
      metric,
      exact_cg_parameters,
      QCGPolicy{},
      exact_cg_line_search);
    const auto qcg_result = qcg_solver.solve(
      PrimalBlock(partition.control_layout(), {DenseVector{1.0, -1.0}}));
    require(qcg_result.stopping_reason ==
              nmopt::solvers::ReducedGradientStoppingReason::gradient_tolerance,
            "Dense quadratic CG solver did not reach tolerance");
    require(qcg_result.accepted_iterations ==
              exact_fr_result.accepted_iterations &&
              qcg_result.line_search_trial_count ==
                qcg_result.accepted_iterations,
            "Dense quadratic CG did not use the exact quadratic recurrence");
    for (std::size_t index = 0;
         index < qcg_result.objective_history.size();
         ++index)
      require_close(qcg_result.objective_history[index],
                    exact_fr_result.objective_history[index],
                    1e-12,
                    "Quadratic CG and Fletcher-Reeves objective equivalence");

    nmopt::solvers::ReducedTrustRegionParameters trust_region_parameters;
    trust_region_parameters.maximum_iterations = 100;
    trust_region_parameters.maximum_trials_per_iteration = 10;
    trust_region_parameters.gradient_tolerance = 1e-8;
    trust_region_parameters.initial_radius = 0.25;
    trust_region_parameters.minimum_radius = 1e-12;
    trust_region_parameters.maximum_radius = 4.0;
    const nmopt::solvers::ReducedTrustRegionSolver trust_region_solver(
      reduced, metric, hessian, trust_region_parameters);
    const auto trust_region_result = trust_region_solver.solve(
      PrimalBlock(partition.control_layout(), {DenseVector{1.0, -1.0}}));
    require(trust_region_result.stopping_reason ==
              nmopt::solvers::ReducedTrustRegionStoppingReason::gradient_tolerance,
            "Trust-region solver did not reach its gradient tolerance");
    require(trust_region_result.accepted_iterations > 1,
            "Trust-region solver did not exercise multiple accepted radii");
    require(trust_region_result.trial_count ==
              trust_region_result.accepted_iterations,
            "Trust-region quadratic target unexpectedly rejected a trial");
    require(trust_region_result.radius_history.size() ==
              trust_region_result.trial_count &&
              trust_region_result.step_norm_history.size() ==
                trust_region_result.trial_count &&
              trust_region_result.predicted_reduction_history.size() ==
                trust_region_result.trial_count &&
              trust_region_result.actual_reduction_history.size() ==
                trust_region_result.trial_count &&
              trust_region_result.reduction_ratio_history.size() ==
                trust_region_result.trial_count &&
              trust_region_result.accepted_history.size() ==
                trust_region_result.trial_count &&
              trust_region_result.subproblem_status_history.size() ==
                trust_region_result.trial_count &&
              trust_region_result.subproblem_iteration_history.size() ==
                trust_region_result.trial_count &&
              trust_region_result.subproblem_residual_norm_history.size() ==
                trust_region_result.trial_count,
            "Trust-region diagnostics do not match trial count");
    for (std::size_t trial = 0;
         trial < trust_region_result.trial_count;
         ++trial)
      {
        require(trust_region_result.accepted_history[trial],
                "Trust-region quadratic target rejected an exact model step");
        require(trust_region_result.subproblem_status_history[trial] ==
                  nmopt::solvers::ReducedTrustRegionSubproblemStatus::cauchy,
                "Default trust-region subproblem is not Cauchy");
        require(trust_region_result.predicted_reduction_history[trial] > 0.0 &&
                  trust_region_result.actual_reduction_history[trial] > 0.0,
                "Trust-region reductions are not positive");
        require(std::isfinite(trust_region_result.reduction_ratio_history[trial]) &&
                  trust_region_result.reduction_ratio_history[trial] > 0.0,
                "Trust-region reduction ratio is not positive and finite");
        if (trial == 0)
          require_close(trust_region_result.reduction_ratio_history[trial],
                        1.0,
                        1e-10,
                        "Trust-region initial actual/predicted reduction ratio");
      }
    require(trust_region_result.metric_solve_count ==
              trust_region_result.gradient_norm_history.size() &&
              trust_region_result.hessian_action_count + 1 ==
                trust_region_result.gradient_norm_history.size(),
            "Trust-region action counts do not match its model evaluations");
    require(trust_region_result.state_solve_count ==
              trust_region_result.trial_count + 1 &&
              trust_region_result.adjoint_solve_count ==
                trust_region_result.accepted_iterations + 1,
            "Trust-region formulation solve counts miss a trial evaluation");
    require(trust_region_result.parameters.subproblem_method ==
              trust_region_parameters.subproblem_method &&
              trust_region_result.policy_name == "trust_region_cauchy" &&
              trust_region_result.iteration_records.size() ==
                trust_region_result.accepted_iterations,
            "Trust-region result does not retain its audit snapshot");
    for (std::size_t index = 0;
         index < trust_region_result.accepted_iterations;
         ++index)
      {
        const auto &record = trust_region_result.iteration_records[index];
        require(record.common.iteration == index + 1 &&
                  record.common.policy_name == "trust_region_cauchy" &&
                  record.common.trial_count > 0 &&
                  record.common.per_iteration_work.state_solve_count == 1 &&
                  record.common.per_iteration_work.adjoint_solve_count == 1 &&
                  record.common.per_iteration_work.metric_solve_count == 1 &&
                  record.radius > 0.0 &&
                  record.predicted_reduction > 0.0 &&
                  record.actual_reduction > 0.0 &&
                  std::isfinite(record.reduction_ratio) &&
                  record.common.accepted_evaluation.objective_value ==
                    trust_region_result.objective_history[index + 1],
                "Trust-region accepted-iteration audit is incomplete");
      }

    nmopt::solvers::ReducedTrustRegionParameters truncated_parameters =
      trust_region_parameters;
    truncated_parameters.subproblem_method =
      nmopt::solvers::ReducedTrustRegionSubproblemMethod::
        truncated_conjugate_gradient;
    truncated_parameters.maximum_subproblem_iterations = 10;
    truncated_parameters.subproblem_relative_tolerance = 1e-12;
    truncated_parameters.subproblem_absolute_tolerance = 1e-14;
    truncated_parameters.initial_radius = 100.0;
    truncated_parameters.maximum_radius = 100.0;
    const nmopt::solvers::ReducedTrustRegionSolver truncated_solver(
      reduced, metric, hessian, truncated_parameters);
    const auto truncated_result = truncated_solver.solve(
      PrimalBlock(partition.control_layout(), {DenseVector{1.0, -1.0}}));
    require(truncated_result.stopping_reason ==
              nmopt::solvers::ReducedTrustRegionStoppingReason::gradient_tolerance,
            "Truncated-CG trust-region solver did not reach its gradient tolerance");
    require(truncated_result.trial_count ==
              truncated_result.accepted_iterations &&
              truncated_result.trial_count > 0,
            "Truncated-CG trust-region solver did not accept its trials");
    require(truncated_result.subproblem_status_history.size() ==
              truncated_result.trial_count &&
              truncated_result.subproblem_iteration_history.size() ==
                truncated_result.trial_count &&
              truncated_result.subproblem_residual_norm_history.size() ==
                truncated_result.trial_count,
            "Truncated-CG subproblem diagnostics do not match trial count");
    std::size_t truncated_hessian_actions = 0;
    bool saw_converged_subproblem = false;
    for (std::size_t trial = 0;
         trial < truncated_result.trial_count;
         ++trial)
      {
        require(truncated_result.subproblem_status_history[trial] ==
                    nmopt::solvers::ReducedTrustRegionSubproblemStatus::
                      converged ||
                  truncated_result.subproblem_status_history[trial] ==
                    nmopt::solvers::ReducedTrustRegionSubproblemStatus::boundary,
                "Truncated-CG subproblem returned an unexpected status");
        require(std::isfinite(
                  truncated_result.subproblem_residual_norm_history[trial]) &&
                  truncated_result.subproblem_residual_norm_history[trial] >= 0.0,
                "Truncated-CG residual diagnostic is invalid");
        truncated_hessian_actions +=
          truncated_result.subproblem_iteration_history[trial];
        saw_converged_subproblem =
          saw_converged_subproblem ||
          truncated_result.subproblem_status_history[trial] ==
            nmopt::solvers::ReducedTrustRegionSubproblemStatus::converged;
        require(truncated_result.accepted_history[trial],
                "Truncated-CG trust-region target rejected an exact model step");
        require_close(truncated_result.reduction_ratio_history[trial],
                      1.0,
                      1e-10,
                      "Truncated-CG actual/predicted reduction ratio");
      }
    require(saw_converged_subproblem,
            "Truncated-CG trust-region solver did not converge its subproblem");
    require(truncated_hessian_actions ==
              truncated_result.hessian_action_count,
            "Truncated-CG subproblem actions do not match Hessian actions");

    nmopt::solvers::ReducedTrustRegionParameters loose_truncated_parameters =
      truncated_parameters;
    loose_truncated_parameters.subproblem_absolute_tolerance = 1e6;
    loose_truncated_parameters.subproblem_relative_tolerance = 0.9;
    loose_truncated_parameters.maximum_iterations = 1;
    const nmopt::solvers::ReducedTrustRegionSolver loose_truncated_solver(
      reduced, metric, hessian, loose_truncated_parameters);
    const auto loose_truncated_result = loose_truncated_solver.solve(
      PrimalBlock(partition.control_layout(), {DenseVector{1.0, -1.0}}));
    require(loose_truncated_result.trial_count == 1 &&
              loose_truncated_result.accepted_history.front() &&
              loose_truncated_result.step_norm_history.front() > 0.0 &&
              loose_truncated_result.predicted_reduction_history.front() > 0.0 &&
              std::isfinite(loose_truncated_result.reduction_ratio_history.front()),
            "Truncated-CG loose inner tolerance returned a zero predicted reduction");

    nmopt::solvers::ReducedTrustRegionParameters boundary_parameters =
      truncated_parameters;
    boundary_parameters.maximum_iterations = 1;
    boundary_parameters.initial_radius = 0.25;
    boundary_parameters.maximum_radius = 4.0;
    const nmopt::solvers::ReducedTrustRegionSolver boundary_solver(
      reduced, metric, hessian, boundary_parameters);
    const auto boundary_result = boundary_solver.solve(
      PrimalBlock(partition.control_layout(), {DenseVector{1.0, -1.0}}));
    require(boundary_result.trial_count == 1 &&
              boundary_result.accepted_history.front(),
            "Truncated-CG boundary trust-region step was not accepted");
    require(boundary_result.subproblem_status_history.front() ==
              nmopt::solvers::ReducedTrustRegionSubproblemStatus::boundary,
            "Truncated-CG trust-region did not report boundary termination");
    require_close(boundary_result.step_norm_history.front(),
                  boundary_parameters.initial_radius,
                  1e-12,
                  "Truncated-CG boundary step does not fill the radius");

    const ScaledReducedHessian under_model_hessian(hessian, 0.1);
    nmopt::solvers::ReducedTrustRegionParameters rejection_parameters =
      trust_region_parameters;
    rejection_parameters.stopping_criterion =
      nmopt::solvers::ReducedStoppingCriterion::step_norm;
    rejection_parameters.gradient_tolerance = 1e6;
    rejection_parameters.step_tolerance = 0.1;
    rejection_parameters.acceptance_threshold = 0.99;
    rejection_parameters.shrink_threshold = 0.995;
    rejection_parameters.expansion_threshold = 0.999;
    const nmopt::solvers::ReducedTrustRegionSolver rejection_solver(
      reduced, metric, under_model_hessian, rejection_parameters);
    const auto rejection_result = rejection_solver.solve(
      PrimalBlock(partition.control_layout(), {DenseVector{1.0, -1.0}}));
    require(rejection_result.stopping_reason ==
              nmopt::solvers::ReducedTrustRegionStoppingReason::step_tolerance,
            "Trust-region solver did not stop after a rejected model step: " +
              std::string(nmopt::solvers::reduced_trust_region_stopping_reason_name(
                rejection_result.stopping_reason)) +
              " after " + std::to_string(rejection_result.accepted_iterations) +
              " accepted iterations");
    require(rejection_result.trial_count >
              rejection_result.accepted_iterations,
            "Trust-region rejection scenario did not reject a trial");
    bool saw_rejected_trial = false;
    for (std::size_t trial = 0;
         trial < rejection_result.accepted_history.size();
         ++trial)
      saw_rejected_trial = saw_rejected_trial ||
                           !rejection_result.accepted_history[trial];
    require(saw_rejected_trial,
            "Trust-region rejection diagnostics omitted the rejected trial");
    require(rejection_result.radius_history.size() > 1 &&
              rejection_result.radius_history[1] <
                rejection_result.radius_history[0],
            "Trust-region radius did not shrink after rejection");

    const nmopt::solvers::ReducedLimitedMemoryBfgsSolver lbfgs_solver(
      reduced, metric, solver_parameters);
    const auto lbfgs_result = lbfgs_solver.solve(
      PrimalBlock(partition.control_layout(), {DenseVector{1.0, -1.0}}));
    require(lbfgs_result.stopping_reason ==
              nmopt::solvers::ReducedGradientStoppingReason::gradient_tolerance,
            "Dense L-BFGS solver did not reach its tolerance");
    require(lbfgs_result.step_length_history.size() ==
              lbfgs_result.accepted_iterations,
            "Dense L-BFGS step history does not match accepted iterations");
    require(lbfgs_result.direction_reset_count <=
              lbfgs_result.accepted_iterations,
            "Dense L-BFGS direction reset count exceeds accepted iterations");
    for (std::size_t index = 1;
         index < lbfgs_result.objective_history.size();
         ++index)
      require(lbfgs_result.objective_history[index] <=
                lbfgs_result.objective_history[index - 1],
              "Dense L-BFGS objective history is not monotonic");

    const nmopt::solvers::ReducedFullBfgsSolver full_bfgs_solver(
      reduced, metric, solver_parameters);
    const auto full_bfgs_result = full_bfgs_solver.solve(
      PrimalBlock(partition.control_layout(), {DenseVector{1.0, -1.0}}));
    require(full_bfgs_result.stopping_reason ==
              nmopt::solvers::ReducedGradientStoppingReason::gradient_tolerance,
            "Dense full BFGS solver did not reach its tolerance");
    require(full_bfgs_result.step_length_history.size() ==
              full_bfgs_result.accepted_iterations,
            "Dense full BFGS step history does not match accepted iterations");
    require(full_bfgs_result.direction_reset_count <=
              full_bfgs_result.accepted_iterations,
            "Dense full BFGS direction reset count exceeds accepted iterations");
    for (std::size_t index = 1;
         index < full_bfgs_result.objective_history.size();
         ++index)
      require(full_bfgs_result.objective_history[index] <=
                full_bfgs_result.objective_history[index - 1],
              "Dense full BFGS objective history is not monotonic");

    const nmopt::solvers::NewtonDirectionPolicyDense exact_newton_direction(
      hessian, newton_parameters);
    const nmopt::solvers::ExactQuadraticLineSearchPolicy exact_newton_line_search(
      hessian);
    const nmopt::solvers::ReducedExactNewtonSolver exact_newton_solver(
      reduced,
      metric,
      solver_parameters,
      exact_newton_direction,
      exact_newton_line_search);
    const auto exact_newton_result = exact_newton_solver.solve(
      PrimalBlock(partition.control_layout(), {DenseVector{1.0, -1.0}}));
    require(exact_newton_result.stopping_reason ==
              nmopt::solvers::ReducedGradientStoppingReason::gradient_tolerance,
            "Dense exact Newton solver did not reach its tolerance");
    require(exact_newton_result.line_search_trial_count ==
              exact_newton_result.accepted_iterations,
            "Dense exact Newton did not use one exact trial per accepted iteration");
    require(exact_newton_result.state_solve_count ==
              exact_newton_result.line_search_trial_count + 1,
            "Dense exact Newton solve count misses a trial evaluation");
    require(exact_newton_result.adjoint_solve_count ==
              exact_newton_result.accepted_iterations + 1,
            "Dense exact Newton adjoint count includes rejected trials");
    require(exact_newton_result.hessian_action_count >
              exact_newton_result.accepted_iterations,
            "Dense exact Newton did not report line-search Hessian actions");

    nmopt::solvers::WolfeLineSearchParameters wolfe_solver_line_parameters;
    wolfe_solver_line_parameters.maximum_trials = 30;
    wolfe_solver_line_parameters.initial_step_length = 10.0;
    const nmopt::solvers::SteepestDescentDirectionPolicy wolfe_direction_policy;
    const nmopt::solvers::WolfeLineSearchPolicy wolfe_solver_line_search(
      wolfe_solver_line_parameters);
    const nmopt::solvers::ReducedWolfeGradientSolver wolfe_solver(
      reduced,
      metric,
      solver_parameters,
      wolfe_direction_policy,
      wolfe_solver_line_search);
    const auto wolfe_solver_result = wolfe_solver.solve(
      PrimalBlock(partition.control_layout(), {DenseVector{1.0, -1.0}}));
    require(wolfe_solver_result.stopping_reason ==
              nmopt::solvers::ReducedGradientStoppingReason::gradient_tolerance,
            "Dense Wolfe solver did not reach its tolerance");
    require(wolfe_solver_result.hessian_action_count == 0,
            "Dense Wolfe solver unexpectedly reported Hessian actions");
    require(wolfe_solver_result.state_solve_count ==
              wolfe_solver_result.line_search_trial_count + 1,
            "Dense Wolfe solve count misses a trial evaluation");
    require(wolfe_solver_result.adjoint_solve_count ==
              wolfe_solver_result.line_search_trial_count + 1,
            "Dense Wolfe adjoint count misses a trial derivative");
    for (std::size_t index = 1;
         index < wolfe_solver_result.objective_history.size();
         ++index)
      require(wolfe_solver_result.objective_history[index] <=
                wolfe_solver_result.objective_history[index - 1],
              "Dense Wolfe objective history is not monotonic");

    const nmopt::solvers::WeakWolfeLineSearchPolicy weak_wolfe_solver_line_search(
      weak_wolfe_parameters);
    const nmopt::solvers::ReducedWeakWolfeGradientSolver weak_wolfe_solver(
      reduced,
      metric,
      solver_parameters,
      wolfe_direction_policy,
      weak_wolfe_solver_line_search);
    const auto weak_wolfe_solver_result = weak_wolfe_solver.solve(
      PrimalBlock(partition.control_layout(), {DenseVector{1.0, -1.0}}));
    require(weak_wolfe_solver_result.stopping_reason ==
              nmopt::solvers::ReducedGradientStoppingReason::gradient_tolerance,
            "Dense weak Wolfe solver did not reach its tolerance");
    require(weak_wolfe_solver_result.state_solve_count ==
              weak_wolfe_solver_result.line_search_trial_count + 1,
            "Dense weak Wolfe solve count misses a trial evaluation");
    require(weak_wolfe_solver_result.adjoint_solve_count ==
              weak_wolfe_solver_result.line_search_trial_count + 1,
            "Dense weak Wolfe adjoint count misses a trial derivative");

    const CellwiseBoxConstraint projected_bounds(
      partition.control_layout(), {DenseVector{-1.0, -0.2}},
      {DenseVector{1.0, 0.5}}, metric);
    const RecordingConstraint recording_bounds(projected_bounds);
    const nmopt::solvers::ReducedGradientSolver projected_solver(
      reduced, metric, recording_bounds, solver_parameters);
    nmopt::test_support::require_contract_error(
      [&projected_solver, &partition]() {
        (void)projected_solver.solve(
          PrimalBlock(partition.control_layout(),
                      {DenseVector{0.0, 0.6}}));
      },
      "Projected reduced gradient requires a feasible initial control",
      "infeasible projected-solver initial control");
    const auto projected_result = projected_solver.solve(
      PrimalBlock(partition.control_layout(), {DenseVector{0.4, 0.4}}));

    require(projected_result.stopping_reason ==
              nmopt::solvers::ReducedGradientStoppingReason::gradient_tolerance,
            "Projected reduced gradient solver did not reach stationarity");
    require(projected_bounds.is_feasible(projected_result.control),
            "Projected reduced gradient solver returned an infeasible control");
    require_close(projected_result.control.block(0)[1],
                  -0.2,
                  1e-12,
                  "Projected reduced gradient did not reach the lower bound");
    require(projected_result.gradient_norm_history.back() <=
              solver_parameters.gradient_tolerance,
            "Projected reduced gradient final norm exceeds the tolerance");
    require(projected_result.metric_solve_count ==
              projected_result.gradient_norm_history.size(),
            "Projected reduced gradient metric solve count does not match direction evaluations");
    require(projected_result.step_length_history.size() ==
              projected_result.accepted_iterations,
            "Projected reduced gradient step history does not match accepted iterations");
    require(!recording_bounds.projected_controls().empty(),
            "Projected reduced gradient did not request a projection");
    for (const PrimalBlock &projected_control :
         recording_bounds.projected_controls())
      require(projected_bounds.is_feasible(projected_control),
              "Projected reduced gradient generated an infeasible iterate");
    for (std::size_t index = 1;
         index < projected_result.objective_history.size();
         ++index)
      require(projected_result.objective_history[index] <=
                projected_result.objective_history[index - 1],
              "Projected reduced gradient objective history is not monotonic");

    const nmopt::solvers::ReducedWolfeGradientSolver projected_wolfe_solver(
      reduced, metric, recording_bounds, solver_parameters);
    const auto projected_wolfe_result = projected_wolfe_solver.solve(
      PrimalBlock(partition.control_layout(), {DenseVector{0.4, 0.4}}));
    require(projected_wolfe_result.stopping_reason ==
              nmopt::solvers::ReducedGradientStoppingReason::gradient_tolerance,
            "Projected strong-Wolfe reduced gradient solver did not reach stationarity");
    require(projected_wolfe_result.control.layout()->compatible_with(
              *partition.control_layout()),
            "Projected strong-Wolfe reduced gradient returned an incompatible control");

    const nmopt::solvers::ReducedWeakWolfeGradientSolver
      projected_weak_wolfe_solver(
        reduced, metric, recording_bounds, solver_parameters);
    const auto projected_weak_wolfe_result = projected_weak_wolfe_solver.solve(
      PrimalBlock(partition.control_layout(), {DenseVector{0.4, 0.4}}));
    require(projected_weak_wolfe_result.stopping_reason ==
              nmopt::solvers::ReducedGradientStoppingReason::gradient_tolerance,
            "Projected weak-Wolfe reduced gradient solver did not reach stationarity");

    const auto expect_projected_policy_rejection =
      [&reduced, &metric, &recording_bounds, &solver_parameters](auto *solver_tag,
                                                                  const char *label) {
        using Solver = std::remove_pointer_t<decltype(solver_tag)>;
        nmopt::test_support::require_contract_error(
          [&reduced, &metric, &recording_bounds, &solver_parameters]() {
            Solver rejected_solver(reduced,
                                   metric,
                                   recording_bounds,
                                   solver_parameters);
            (void)rejected_solver;
          },
          "Projected reduced solver supports only steepest descent with Armijo, weak Wolfe, or strong Wolfe line search",
          label);
      };
    expect_projected_policy_rejection(
      static_cast<nmopt::solvers::ReducedConjugateGradientSolver *>(nullptr),
      "projected nonlinear CG rejection");
    expect_projected_policy_rejection(
      static_cast<nmopt::solvers::ReducedFletcherReevesSolver *>(nullptr),
      "projected Fletcher-Reeves rejection");
    expect_projected_policy_rejection(
      static_cast<nmopt::solvers::ReducedQuadraticConjugateGradientSolver *>(nullptr),
      "projected quadratic CG rejection");
    expect_projected_policy_rejection(
      static_cast<nmopt::solvers::ReducedLimitedMemoryBfgsSolver *>(nullptr),
      "projected L-BFGS rejection");
    expect_projected_policy_rejection(
      static_cast<nmopt::solvers::ReducedFullBfgsSolver *>(nullptr),
      "projected full BFGS rejection");
    expect_projected_policy_rejection(
      static_cast<nmopt::solvers::ReducedNewtonSolver *>(nullptr),
      "projected Newton rejection");
    expect_projected_policy_rejection(
      static_cast<nmopt::solvers::ReducedExactConjugateGradientSolver *>(nullptr),
      "projected exact nonlinear CG rejection");
    expect_projected_policy_rejection(
      static_cast<nmopt::solvers::ReducedExactFletcherReevesSolver *>(nullptr),
      "projected exact Fletcher-Reeves rejection");
    expect_projected_policy_rejection(
      static_cast<nmopt::solvers::ReducedExactNewtonSolver *>(nullptr),
      "projected exact Newton rejection");
  }

  void
  test_staged_reduced_evaluation()
  {
    const LinearQuadraticModel base_model(
      DenseMatrix(2, 2, {4.0, -1.0, -1.0, 3.0}),
      DenseMatrix(2, 2, {1.0, 0.5, -0.25, 2.0}),
      DenseVector{1.0, -0.5},
      DenseMatrix(2, 2, {1.0, 0.0, 0.5, 1.0}),
      DenseVector{0.25, -1.0},
      DenseVector{1.5, 0.75},
      DenseVector{2.0, 3.0},
      0.4);
    const CountingLinearQuadraticModel model(base_model);
    const StateControlPartition partition(model, 0, 1);

    std::size_t state_calls   = 0;
    std::size_t adjoint_calls = 0;
    const StateAdjointSolvers solvers{
      [&base_model, &state_calls](const PrimalBlock &control) {
        ++state_calls;
        return FormulationSolveResultT<DenseBackend>(
          base_model.solve_state(control));
      },
      [&base_model, &adjoint_calls](const PrimalBlock &full_point,
                                    const CovectorBlock &state_rhs) {
        ++adjoint_calls;
        return FormulationSolveResultT<DenseBackend>(
          base_model.solve_adjoint(full_point, state_rhs));
      }};
    const ReducedDTO reduced(model, partition, solvers);
    const PrimalBlock control(partition.control_layout(),
                              {DenseVector{0.4, -0.3}});

    const ReducedValueEvaluation value = reduced.evaluate_value(control);
    require(state_calls == 1 && adjoint_calls == 0,
            "Reduced value evaluation performed unexpected solve work");
    require(model.objective_calls == 1 &&
              model.objective_derivative_calls == 0,
            "Reduced value evaluation performed unexpected objective work");
    require(value.state_solve.converged(),
            "Reduced value evaluation did not retain state solve evidence");

    const ReducedEvaluation staged = reduced.augment_derivative(value);
    require(state_calls == 1 && adjoint_calls == 1,
            "Reduced derivative augmentation repeated the state solve");
    require(model.objective_calls == 1 &&
              model.objective_derivative_calls == 1,
            "Reduced derivative augmentation repeated objective evaluation");
    require(staged.adjoint_solve.converged(),
            "Reduced derivative augmentation did not retain adjoint evidence");
    require_close(staged.objective_value,
                  value.objective_value,
                  0.0,
                  "Staged reduced evaluation objective reuse");

    const ReducedEvaluation compatibility = reduced.evaluate(control);
    require(state_calls == 2 && adjoint_calls == 2,
            "Compatibility reduced evaluation did not compose both stages");
    require(model.objective_calls == 2 &&
              model.objective_derivative_calls == 2,
            "Compatibility reduced evaluation changed objective work");
    require_close(compatibility.objective_value,
                  staged.objective_value,
                  0.0,
                  "Compatibility reduced evaluation objective");
    require_close(pair(compatibility.reduced_derivative,
                       PrimalBlock(partition.control_layout(),
                                   {DenseVector{1.0, -2.0}})),
                  pair(staged.reduced_derivative,
                       PrimalBlock(partition.control_layout(),
                                   {DenseVector{1.0, -2.0}})),
                  1e-15,
                  "Compatibility reduced evaluation derivative");

    const ReducedDTO foreign(model, partition, solvers);
    const ReducedValueEvaluation foreign_value =
      foreign.evaluate_value(control);
    const std::size_t adjoint_calls_before_rejection = adjoint_calls;
    nmopt::test_support::require_contract_error(
      [&reduced, &foreign_value]() {
        (void)reduced.augment_derivative(foreign_value);
      },
      "Reduced DTO value evaluation belongs to another service",
      "foreign reduced value evaluation");
    require(adjoint_calls == adjoint_calls_before_rejection,
            "Rejected foreign value evaluation invoked the adjoint solve");
  }

  void
  test_owned_reduced_service_lifetime()
  {
    auto detached = []() {
      auto model = std::make_shared<LinearQuadraticModel>(
        DenseMatrix(2, 2, {4.0, -1.0, -1.0, 3.0}),
        DenseMatrix(2, 2, {1.0, 0.5, -0.25, 2.0}),
        DenseVector{1.0, -0.5},
        DenseMatrix(2, 2, {1.0, 0.0, 0.5, 1.0}),
        DenseVector{0.25, -1.0},
        DenseVector{1.5, 0.75},
        DenseVector{2.0, 3.0},
        0.4);
      const auto *model_view = model.get();
      StateControlPartition partition(*model_view, 0, 1);
      StateAdjointSolvers solvers{
        [model_view](const PrimalBlock &control) {
          return model_view->solve_state(control);
        },
        [model_view](const PrimalBlock &full_point,
                     const CovectorBlock &state_rhs) {
          return model_view->solve_adjoint(full_point, state_rhs);
        }};
      PrimalBlock control(partition.control_layout(),
                          {DenseVector{0.4, -0.3}});
      struct DetachedService
      {
        ReducedDTO  reduced;
        PrimalBlock control;
      };
      return DetachedService{
        ReducedDTO(std::move(model),
                   std::move(partition),
                   std::move(solvers),
                   std::make_shared<const int>(7)),
        std::move(control)};
    }();

    const auto evaluation = detached.reduced.evaluate(detached.control);
    require(evaluation.state_solve.converged() &&
              evaluation.adjoint_solve.converged(),
            "Detached owned reduced service lost its solve callbacks");
  }

  void
  test_experiment_envelope()
  {
    using CompiledProblem =
      nmopt::compiler::v1::CompiledProblemT<DenseBackend>;
    using Envelope =
      nmopt::experiment::ReducedSearchExperimentEnvelopeT<DenseBackend>;
    using Manifest = nmopt::compiler::v1::CompilationManifest;
    using Environment = nmopt::experiment::RunEnvironmentRecord;

    const auto envelope = []() {
      auto model = std::make_shared<LinearQuadraticModel>(
        DenseMatrix(2, 2, {4.0, -1.0, -1.0, 3.0}),
        DenseMatrix(2, 2, {1.0, 0.5, -0.25, 2.0}),
        DenseVector{1.0, -0.5},
        DenseMatrix(2, 2, {1.0, 0.0, 0.5, 1.0}),
        DenseVector{0.25, -1.0},
        DenseVector{1.5, 0.75},
        DenseVector{2.0, 3.0},
        0.4);
      auto metric = std::make_shared<DiagonalMetric>(
        "l2_cellwise",
        model->control_layout(),
        std::vector<DenseVector>{DenseVector{2.0, 5.0}});
      const StateAdjointSolvers solvers{
        [model](const PrimalBlock &control) {
          return model->solve_state(control);
        },
        [model](const PrimalBlock &full_point,
                const CovectorBlock &state_rhs) {
          return model->solve_adjoint(full_point, state_rhs);
        }};

      Manifest manifest;
      manifest.semantic_problem_id = "reference.scalar.reduced.envelope";
      manifest.compiler_id = "reference";
      manifest.backend = "dense";
      manifest.execution = "assembled";
      manifest.provenance = "DTO";
      manifest.mesh_record.provenance = "manufactured scalar mesh";
      manifest.mesh_record.structural_identity = "mesh-a";
      manifest.formulation_record.semantic_id = "reduced_dto";
      manifest.formulation_record.kind =
        nmopt::semantic::v1::FormulationKind::reduced_dto;
      manifest.formulation_record.provenance =
        nmopt::semantic::v1::FormulationProvenance::dto;

      CompiledProblem compiled_problem(
        model,
        metric,
        std::shared_ptr<const Constraint>{},
        solvers,
        manifest);
      const auto reduced = compiled_problem.make_reduced_dto();
      nmopt::solvers::ReducedSolverParameters parameters;
      parameters.gradient_tolerance = 1e-8;
      parameters.initial_step_length = 2.0;
      const nmopt::solvers::ReducedGradientSolver solver(
        reduced, compiled_problem.metric(), parameters);
      const PrimalBlock control(model->control_layout(),
                                {DenseVector{1.0, -1.0}});
      auto report = solver.solve(control);
      const auto policy =
        nmopt::experiment::make_reduced_search_policy_snapshot(report);
      Environment environment{"revision-a",
                              "debug-neutral",
                              "GNU",
                              "test-version",
                              "libstdc++",
                              "test-os",
                              "x86_64",
                              "test-host"};
      return Envelope(compiled_problem.manifest(),
                      policy,
                      std::move(report),
                      std::move(environment));
    }();

    require(envelope.compilation_manifest().semantic_problem_id ==
              "reference.scalar.reduced.envelope" &&
              envelope.compilation_manifest().mesh_record.structural_identity ==
                "mesh-a",
            "Experiment envelope did not retain the detached compilation manifest");
    require(envelope.solver_policy().solver_name == "reduced_search" &&
              envelope.solver_policy().policy_name == "armijo" &&
              envelope.solver_policy().line_search_parameters.policy_name ==
                "armijo" &&
              envelope.report().policy_name ==
                envelope.solver_policy().policy_name,
            "Experiment envelope did not retain its typed policy snapshot");
    require(envelope.environment().source_revision == "revision-a" &&
              envelope.environment().build_profile == "debug-neutral" &&
              envelope.environment().hardware == "test-host",
            "Experiment envelope did not retain its environment record");

    auto changed_manifest = envelope.compilation_manifest();
    changed_manifest.mesh_record.structural_identity = "mesh-b";
    require(changed_manifest.mesh_record.structural_identity !=
              envelope.compilation_manifest().mesh_record.structural_identity,
            "Experiment envelope manifest identity did not distinguish products");

    auto changed_policy = envelope.solver_policy();
    changed_policy.line_search_parameters.backtracking_factor = 0.25;
    require(changed_policy.line_search_parameters.backtracking_factor !=
              envelope.solver_policy().line_search_parameters.backtracking_factor,
            "Experiment policy snapshot did not retain an independent change");

    auto changed_environment = envelope.environment();
    changed_environment.hardware = "other-host";
    require(changed_environment.hardware != envelope.environment().hardware,
            "Experiment environment record did not retain an independent change");

    nmopt::test_support::require_contract_error(
      [&envelope]() {
        Manifest missing_identifier = envelope.compilation_manifest();
        missing_identifier.semantic_problem_id.clear();
        (void)Envelope(missing_identifier,
                       envelope.solver_policy(),
                       envelope.report(),
                       envelope.environment());
      },
      "An experiment envelope needs a compilation manifest identifier",
      "experiment envelope missing manifest identifier");
  }

  void
  test_backend_parameterisation()
  {
    const auto layout = std::make_shared<const BlockLayout>(
      "alternate",
      std::vector<SpaceId>{{"alternate_control"}},
      std::vector<std::size_t>{2});

    const PrimalBlockT<AlternateDenseBackend> primal(
      layout, {DenseVector{1.0, -2.0}});
    const CovectorBlockT<AlternateDenseBackend> covector(
      layout, {DenseVector{3.0, 4.0}});
    require_close(pair(covector, primal),
                  -5.0,
                  1e-15,
                  "Backend-parametric pairing");

    CovectorBlockT<AlternateDenseBackend> difference =
      subtract(covector,
               CovectorBlockT<AlternateDenseBackend>(
                 layout, {DenseVector{1.0, -1.0}}));
    require_close(difference.block(0)[0],
                  2.0,
                  1e-15,
                  "Backend-parametric covector subtraction");
    require_close(difference.block(0)[1],
                  5.0,
                  1e-15,
                  "Backend-parametric covector subtraction");

    AlternateQuadraticModel model;
    using Backend = AlternateDenseBackend;
    using AlternatePrimal = PrimalBlockT<Backend>;
    using AlternateCovector = CovectorBlockT<Backend>;
    using AlternatePartition = StateControlPartitionT<Backend>;
    using AlternateSolvers = StateAdjointSolversT<Backend>;
    const AlternatePartition partition(model, 0, 1);
    const AlternateSolvers solvers{
      [&model](const AlternatePrimal &control) {
        return model.solve_state(control);
      },
      [&model](const AlternatePrimal &full_point,
               const AlternateCovector &state_rhs) {
        return model.solve_adjoint(full_point, state_rhs);
      }};
    const ReducedDTOT<Backend> reduced(model, partition, solvers);
    const AlternateIdentityMetric metric(partition.control_layout());

    nmopt::solvers::ReducedTrustRegionParameters trust_region_parameters;
    trust_region_parameters.maximum_iterations = 20;
    trust_region_parameters.gradient_tolerance = 1e-10;
    trust_region_parameters.initial_radius = 0.5;
    trust_region_parameters.maximum_radius = 8.0;
    const nmopt::solvers::ReducedTrustRegionSolverT<Backend> solver(
      reduced, metric, model, trust_region_parameters);
    const AlternatePrimal initial_control(
      partition.control_layout(), {DenseVector{0.0, 0.0}});
    const auto result = solver.solve(initial_control);

    require(result.stopping_reason ==
              nmopt::solvers::ReducedTrustRegionStoppingReason::gradient_tolerance,
            "Alternate backend trust-region solver did not converge");
    require(result.accepted_iterations > 1,
            "Alternate backend trust-region solver accepted too few steps");
    require(result.trial_count == result.accepted_iterations,
            "Alternate backend trust-region solver unexpectedly rejected a trial");
    require(result.metric_solve_count == result.gradient_norm_history.size(),
            "Alternate backend trust-region metric count is inconsistent");
    require(result.hessian_action_count + 1 ==
              result.gradient_norm_history.size(),
            "Alternate backend trust-region Hessian count is inconsistent");
    require(result.state_solve_count == result.trial_count + 1 &&
              result.adjoint_solve_count == result.accepted_iterations + 1,
            "Alternate backend trust-region solve counts are inconsistent");
    require_close(result.control.block(0)[0], 1.0, 1e-8,
                  "Alternate backend trust-region first control");
    require_close(result.control.block(0)[1], -2.0, 1e-8,
                  "Alternate backend trust-region second control");
  }

  void
  test_projection_compatibility()
  {
    const auto layout = std::make_shared<const BlockLayout>(
      "projection_compatibility",
      std::vector<SpaceId>{{"control"}},
      std::vector<std::size_t>{2});
    const DiagonalMetric trusted_metric(
      "l2_cellwise", layout, {DenseVector{2.0, 2.0}});
    const CellwiseBoxConstraint bounds(
      layout,
      {DenseVector{0.0, 0.0}},
      {DenseVector{1.0, 1.0}},
      trusted_metric);
    const NonDiagonalMetric spoofed_metric(layout);

    require(bounds.supports_projection_in(trusted_metric),
            "Cellwise box rejected its diagonal L2 metric");
    require(!bounds.supports_projection_in(spoofed_metric),
            "A non-diagonal metric obtained clipping projection by reusing the l2_cellwise display identifier");
  }
} // namespace

int
main(const int argc, char **argv)
{
  try
    {
      const std::vector<nmopt::test_support::Scenario> scenarios{
        {"layout_invariant",
         "nmopt.contract.layout_invariant",
         {"backend-neutral", "contract"},
         30,
         test_block_layout_invariant},
        {"v0_contract",
         "nmopt.contract.v0",
         {"backend-neutral", "contract"},
         30,
         test_v0_contract},
        {"staged_reduced_evaluation",
         "nmopt.contract.staged_reduced_evaluation",
         {"backend-neutral", "contract", "reduced"},
         30,
         test_staged_reduced_evaluation},
        {"backend_parameterisation",
         "nmopt.contract.backend_parameterisation",
         {"backend-neutral", "contract"},
         30,
         test_backend_parameterisation},
        {"owned_reduced_service_lifetime",
         "nmopt.contract.owned_reduced_service_lifetime",
         {"backend-neutral", "contract", "ownership"},
         30,
         test_owned_reduced_service_lifetime},
        {"experiment_envelope",
         "nmopt.contract.experiment_envelope",
         {"backend-neutral", "contract", "ownership"},
         30,
         test_experiment_envelope},
        {"projection_compatibility",
         "nmopt.contract.projection_compatibility",
         {"backend-neutral", "contract", "constraint"},
         30,
         test_projection_compatibility}};
      const auto result =
        nmopt::test_support::run_requested_scenarios(
          argc, argv, scenarios, std::cout);
      if (!result.listed)
        std::cout << "nmopt executable contract scenario passed: "
                  << result.executed << '\n';
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "nmopt v0 executable contract test failed: "
                << exception.what() << '\n';
      return 1;
    }
}
