#pragma once

#include "../diagnostics/instrumentation.hpp"
#include "problem_b.hpp"
#include "problem_b_metric.hpp"

#include "nmopt/contract/callback_executable_model.hpp"
#include "nmopt/contract/linear_solve.hpp"
#include "nmopt/contract/metric_constraint.hpp"
#include "nmopt/contract/reduced_dto.hpp"
#include "nmopt/dealii/serial_backend.hpp"

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace external_dealii_step4
{
  using ProblemBNmoptLayoutPtr = nmopt::contract::LayoutPtr;

  template <int dim>
  class ProblemBNmoptMetric final
    : public nmopt::contract::MetricT<nmopt::dealii_backend::SerialBackend>
  {
  public:
    using Backend  = nmopt::dealii_backend::SerialBackend;
    using Primal   = nmopt::contract::PrimalBlockT<Backend>;
    using Covector = nmopt::contract::CovectorBlockT<Backend>;
    using NativeMetric = ProblemBMetric<dim>;

    ProblemBNmoptMetric(ProblemBNmoptLayoutPtr layout,
                        const NativeMetric &  native_metric,
                        Instrumentation *      instrumentation)
      : layout_(std::move(layout))
      , native_metric_(native_metric)
      , instrumentation_(instrumentation)
      , id_("external_step4_problem_b_mass")
    {
      if (!layout_)
        throw std::invalid_argument("Step-4 Problem B metric needs a layout");
    }

    const std::string &
    id() const override
    {
      return id_;
    }

    const ProblemBNmoptLayoutPtr &
    layout() const override
    {
      return layout_;
    }

    Covector
    apply(const Primal &primal) const override
    {
      count(&Instrumentation::metric_apply_calls);
      require_compatible(primal.layout(),
                         "Step-4 Problem B metric primal");
      return Covector(layout_, {native_metric_.apply(primal.block(0))});
    }

    Primal
    inverse_apply(const Covector &covector) const override
    {
      count(&Instrumentation::metric_inverse_apply_calls);
      require_compatible(covector.layout(),
                         "Step-4 Problem B metric covector");
      auto result = native_metric_.inverse_apply(covector.block(0));
      if (!result.evidence.converged)
        throw std::runtime_error(
          "Step-4 Problem B metric inverse did not converge");
      return Primal(layout_, {std::move(result.solution)});
    }

  private:
    using Counter = std::size_t Instrumentation::*;

    void
    count(const Counter counter) const
    {
      if (instrumentation_ != nullptr)
        ++(instrumentation_->*counter);
    }

    void
    require_compatible(const ProblemBNmoptLayoutPtr &actual,
                       const char *operation) const
    {
      if (!actual || !actual->compatible_with(*layout_))
        throw std::invalid_argument(std::string(operation) +
                                    " has an incompatible layout");
    }

    ProblemBNmoptLayoutPtr layout_;
    const NativeMetric &    native_metric_;
    Instrumentation *       instrumentation_;
    std::string             id_;
  };

  template <int dim, typename Application>
  class ProblemBNmoptBinding final
  {
  public:
    using Backend  = nmopt::dealii_backend::SerialBackend;
    using Problem  = ProblemB<dim, Application>;
    using Primal   = nmopt::contract::PrimalBlockT<Backend>;
    using Covector = nmopt::contract::CovectorBlockT<Backend>;
    using Model = nmopt::contract::CallbackExecutableModelT<Backend>;
    using Partition = nmopt::contract::StateControlPartitionT<Backend>;
    using Solvers = nmopt::contract::StateAdjointSolversT<Backend>;
    using Reduced = nmopt::contract::ReducedDTOT<Backend>;
    using SolveResult = nmopt::contract::FormulationSolveResultT<Backend>;
    using Metric = ProblemBNmoptMetric<dim>;

    ProblemBNmoptBinding()
      : ProblemBNmoptBinding(nullptr)
    {}

    explicit ProblemBNmoptBinding(Instrumentation &instrumentation)
      : ProblemBNmoptBinding(&instrumentation)
    {}

    ProblemBNmoptBinding(const ProblemBNmoptBinding &) = delete;
    ProblemBNmoptBinding &operator=(const ProblemBNmoptBinding &) = delete;
    ProblemBNmoptBinding(ProblemBNmoptBinding &&) = delete;
    ProblemBNmoptBinding &operator=(ProblemBNmoptBinding &&) = delete;

    const Problem &
    problem() const
    {
      return problem_;
    }

    const ProblemBNmoptLayoutPtr &
    variable_layout() const
    {
      return variable_layout_;
    }

    const ProblemBNmoptLayoutPtr &
    test_layout() const
    {
      return test_layout_;
    }

    const ProblemBNmoptLayoutPtr &
    control_layout() const
    {
      return partition_.control_layout();
    }

    const ProblemBNmoptLayoutPtr &
    state_layout() const
    {
      return partition_.state_layout();
    }

    const Model &
    model() const
    {
      return model_;
    }

    const Metric &
    metric() const
    {
      return metric_;
    }

    const Reduced &
    reduced() const
    {
      return reduced_;
    }

  private:
    using ApplicationPtr = std::unique_ptr<Application>;
    using Counter = std::size_t Instrumentation::*;

    explicit ProblemBNmoptBinding(Instrumentation *const instrumentation)
      : application_(make_application())
      , problem_(*application_)
      , native_metric_(problem_.mass())
      , variable_layout_(std::make_shared<const nmopt::contract::BlockLayout>(
          "external_step4_problem_b_variables",
          std::vector<nmopt::contract::SpaceId>{{"state"}, {"control"}},
          std::vector<std::size_t>{problem_.state_dimension(),
                                   problem_.control_dimension()}))
      , test_layout_(std::make_shared<const nmopt::contract::BlockLayout>(
          "external_step4_problem_b_test",
          std::vector<nmopt::contract::SpaceId>{{"state_test"}},
          std::vector<std::size_t>{problem_.state_dimension()}))
      , model_(make_model(problem_, variable_layout_, test_layout_,
                          instrumentation))
      , partition_(model_, 0, 1)
      , solvers_(make_solvers(problem_, partition_.state_layout(), test_layout_,
                              instrumentation))
      , metric_(partition_.control_layout(), native_metric_, instrumentation)
      , reduced_(model_, partition_, solvers_)
    {}

    static ApplicationPtr
    make_application()
    {
      auto application = std::make_unique<Application>();
      application->prepare_for_external_use();
      return application;
    }

    static void
    count(Instrumentation *const instrumentation, const Counter counter)
    {
      if (instrumentation != nullptr)
        ++(instrumentation->*counter);
    }

    static void
    record_solve(Instrumentation *const instrumentation,
                 const SolveRole         role,
                 const typename Problem::SolveEvidence &evidence)
    {
      if (instrumentation == nullptr)
        return;
      if (evidence.converged)
        instrumentation->record_solve_success(role,
                                               evidence.iterations,
                                               evidence.initial_residual,
                                               evidence.final_residual);
      else
        instrumentation->record_solve_failure(role,
                                               evidence.iterations,
                                               evidence.final_residual);
    }

    static nmopt::contract::LinearSolveReport
    solve_report(const typename Problem::SolveEvidence &evidence)
    {
      return {"CG",
              "identity",
              1000,
              evidence.iterations,
              0.0,
              1.0e-12,
              1.0e-12,
              evidence.final_residual,
              evidence.converged ?
                nmopt::contract::LinearSolveTermination::converged :
                nmopt::contract::LinearSolveTermination::failed};
    }

    static Model
    make_model(Problem &                     problem,
               const ProblemBNmoptLayoutPtr &variable_layout,
               const ProblemBNmoptLayoutPtr &test_layout,
               Instrumentation *const        instrumentation)
    {
      Problem *const problem_ptr = &problem;

      return Model(
        variable_layout,
        test_layout,
        [problem_ptr, test_layout, instrumentation](const Primal &variables) {
          count(instrumentation, &Instrumentation::residual_calls);
          auto value =
            problem_ptr->residual(variables.block(0), variables.block(1));
          return Covector(test_layout, {std::move(value)});
        },
        [problem_ptr, test_layout, instrumentation](
          const Primal &, const Primal &variable_tangent) {
          count(instrumentation, &Instrumentation::residual_jvp_calls);
          auto value = problem_ptr->residual_jvp(
            variable_tangent.block(0), variable_tangent.block(1));
          return Covector(test_layout, {std::move(value)});
        },
        [problem_ptr, variable_layout, instrumentation](
          const Primal &, const Primal &test_seed) {
          count(instrumentation, &Instrumentation::residual_vjp_calls);
          auto value = problem_ptr->residual_vjp(test_seed.block(0));
          return Covector(variable_layout,
                          {std::move(value.state), std::move(value.control)});
        },
        [problem_ptr, instrumentation](const Primal &variables) {
          count(instrumentation, &Instrumentation::objective_calls);
          return problem_ptr->objective(variables.block(0), variables.block(1));
        },
        [problem_ptr, variable_layout, instrumentation](
          const Primal &variables) {
          count(instrumentation,
                &Instrumentation::objective_derivative_calls);
          auto value = problem_ptr->objective_derivative(variables.block(0),
                                                         variables.block(1));
          return Covector(variable_layout,
                          {std::move(value.state), std::move(value.control)});
        });
    }

    static Solvers
    make_solvers(Problem &                     problem,
                 const ProblemBNmoptLayoutPtr &state_layout,
                 const ProblemBNmoptLayoutPtr &test_layout,
                 Instrumentation *const        instrumentation)
    {
      Problem *const problem_ptr = &problem;
      Solvers       solvers;
      solvers.solve_state = [problem_ptr, state_layout, instrumentation](
                              const Primal &control) {
        if (instrumentation != nullptr)
          instrumentation->record_solve_start(SolveRole::state);
        try
          {
            auto result = problem_ptr->solve_state(control.block(0));
            record_solve(instrumentation, SolveRole::state, result.evidence);
            return SolveResult(
              Primal(state_layout, {std::move(result.solution)}),
              solve_report(result.evidence));
          }
        catch (...)
          {
            if (instrumentation != nullptr)
              instrumentation->record_solve_failure();
            throw;
          }
      };
      solvers.solve_adjoint = [problem_ptr, test_layout, instrumentation](
                                const Primal &,
                                const Covector &state_rhs) {
        if (instrumentation != nullptr)
          instrumentation->record_solve_start(SolveRole::adjoint);
        try
          {
            auto result = problem_ptr->solve_adjoint(state_rhs.block(0));
            record_solve(instrumentation,
                         SolveRole::adjoint,
                         result.evidence);
            return SolveResult(
              Primal(test_layout, {std::move(result.solution)}),
              solve_report(result.evidence));
          }
        catch (...)
          {
            if (instrumentation != nullptr)
              instrumentation->record_solve_failure();
            throw;
          }
      };
      return solvers;
    }

    ApplicationPtr        application_;
    Problem               problem_;
    ProblemBMetric<dim>   native_metric_;
    ProblemBNmoptLayoutPtr variable_layout_;
    ProblemBNmoptLayoutPtr test_layout_;
    Model                 model_;
    Partition             partition_;
    Solvers               solvers_;
    Metric                metric_;
    Reduced               reduced_;
  };
} // namespace external_dealii_step4
