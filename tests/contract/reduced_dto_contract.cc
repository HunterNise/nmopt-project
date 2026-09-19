#include "nmopt/contract/metric_constraint.hpp"
#include "nmopt/contract/reduced_dto.hpp"
#include "nmopt/compiler/v1/compiled_products.hpp"
#include "nmopt/experiment/reduced_envelope.hpp"
#include "nmopt/reference/linear_quadratic_model.hpp"
#include "nmopt/solvers/reduced_gradient.hpp"
#include "nmopt/solvers/reduced_line_search.hpp"
#include "nmopt/solvers/reduced_trust_region.hpp"
#include "../support/contract_errors.hpp"
#include "../support/scenario_dispatch.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <limits>
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

  class NegativeIdentityReducedHessian final : public ReducedHessian
  {
  public:
    explicit NegativeIdentityReducedHessian(LayoutPtr layout)
      : layout_(std::move(layout))
    {}

    const LayoutPtr &
    layout() const override
    {
      return layout_;
    }

    CovectorBlock
    apply(const PrimalBlock &control,
          const PrimalBlock &direction) const override
    {
      require(control.layout()->compatible_with(*layout_) &&
                direction.layout()->compatible_with(*layout_),
              "Negative identity Hessian received incompatible controls");
      DenseVector result = direction.block(0);
      result.scale(-1.0);
      return CovectorBlock(layout_, {std::move(result)});
    }

  private:
    LayoutPtr layout_;
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


#include "../solvers/scenarios/reduced_search.hpp"

  void
  test_exact_quadratic_line_search_nonpositive_curvature_diagnostic()
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
    const StateControlPartition partition(model, 0, 1);
    const StateAdjointSolvers solvers{
      [&model](const PrimalBlock &control) { return model.solve_state(control); },
      [&model](const PrimalBlock &full_point, const CovectorBlock &state_rhs) {
        return model.solve_adjoint(full_point, state_rhs);
      }};
    const ReducedDTO reduced(model, partition, solvers);
    const PrimalBlock control(partition.control_layout(),
                              {DenseVector{0.4, -0.3}});
    const ReducedEvaluation evaluation = reduced.evaluate(control);
    const DiagonalMetric metric(
      "l2_cellwise", partition.control_layout(), {DenseVector{2.0, 5.0}});
    const auto search_direction = nmopt::solvers::make_steepest_descent_direction(
      evaluation.reduced_derivative, metric);
    const auto build_trial =
      [&control, &search_direction](const double step) {
        return shifted(control, search_direction.direction, step);
      };
    const auto evaluate_trial_value =
      [&reduced](const PrimalBlock &trial_control) {
        return reduced.evaluate_value(trial_control);
      };
    const auto augment_trial_derivative =
      [&reduced](const ReducedValueEvaluation &value) {
        return reduced.augment_derivative(value);
      };
    const NegativeIdentityReducedHessian negative_hessian(
      partition.control_layout());
    const nmopt::solvers::ExactQuadraticLineSearchPolicy exact_policy(
      negative_hessian);

    nmopt::test_support::require_contract_error(
      [&exact_policy, &control, &evaluation, &search_direction,
       &build_trial, &evaluate_trial_value,
       &augment_trial_derivative]() {
        (void)exact_policy.search(control,
                                  evaluation,
                                  search_direction,
                                  build_trial,
                                  evaluate_trial_value,
                                  augment_trial_derivative);
      },
      "Exact quadratic line search requires positive curvature",
      "exact line search non-positive curvature");
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
      manifest.resolved_decision.semantic_problem_id = "reference.scalar.reduced.envelope";
      manifest.compatibility.compiler_id = "reference";
      manifest.compatibility.backend = "dense";
      manifest.compatibility.execution = "assembled";
      manifest.compatibility.provenance = "DTO";
      manifest.resolved_decision.mesh_record.provenance = "manufactured scalar mesh";
      manifest.resolved_decision.mesh_record.structural_identity = "mesh-a";
      manifest.resolved_decision.formulation_record.semantic_id = "reduced_dto";
      manifest.resolved_decision.formulation_record.kind =
        nmopt::semantic::v1::FormulationKind::reduced_dto;
      manifest.resolved_decision.formulation_record.provenance =
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

    require(envelope.compilation_manifest().resolved_decision.semantic_problem_id ==
              "reference.scalar.reduced.envelope" &&
              envelope.compilation_manifest().resolved_decision.mesh_record.structural_identity ==
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
    changed_manifest.resolved_decision.mesh_record.structural_identity = "mesh-b";
    require(changed_manifest.resolved_decision.mesh_record.structural_identity !=
              envelope.compilation_manifest().resolved_decision.mesh_record.structural_identity,
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
        missing_identifier.resolved_decision.semantic_problem_id.clear();
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
        {"exact_quadratic_line_search_nonpositive_curvature_diagnostic",
         "nmopt.reduced.exact_quadratic_line_search_nonpositive_curvature_diagnostic",
         {"backend-neutral", "contract", "reduced", "negative"},
         30,
         test_exact_quadratic_line_search_nonpositive_curvature_diagnostic},
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
