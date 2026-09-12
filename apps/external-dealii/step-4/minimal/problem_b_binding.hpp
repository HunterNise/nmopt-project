#pragma once

#define STEP4_NO_MAIN
#include "../source/adapted/step-4.cc"
#undef STEP4_NO_MAIN

#include "../integration/problem_b.hpp"
#include "../integration/problem_b_metric.hpp"

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
  namespace minimal
  {
    using LayoutPtr = nmopt::contract::LayoutPtr;

    class ProblemBMassMetric final
      : public nmopt::contract::MetricT<nmopt::dealii_backend::SerialBackend>
    {
    public:
      using Backend  = nmopt::dealii_backend::SerialBackend;
      using Primal   = nmopt::contract::PrimalBlockT<Backend>;
      using Covector = nmopt::contract::CovectorBlockT<Backend>;
      using NativeMetric = ProblemBMetric<2>;

      ProblemBMassMetric(LayoutPtr layout, const NativeMetric &native_metric)
        : layout_(std::move(layout))
        , native_metric_(native_metric)
        , id_("external_step4_minimal_problem_b_mass")
      {
        if (!layout_)
          throw std::invalid_argument(
            "minimal Step-4 Problem B metric needs a layout");
      }

      const std::string &
      id() const override
      {
        return id_;
      }

      const LayoutPtr &
      layout() const override
      {
        return layout_;
      }

      Covector
      apply(const Primal &primal) const override
      {
        require_compatible(primal.layout(),
                           "minimal Step-4 Problem B metric primal");
        return Covector(layout_, {native_metric_.apply(primal.block(0))});
      }

      Primal
      inverse_apply(const Covector &covector) const override
      {
        require_compatible(covector.layout(),
                           "minimal Step-4 Problem B metric covector");
        const auto result = native_metric_.inverse_apply(covector.block(0));
        if (!result.evidence.converged)
          throw std::runtime_error(
            "minimal Step-4 Problem B metric inverse did not converge");
        return Primal(layout_, {result.solution});
      }

    private:
      void
      require_compatible(const LayoutPtr &actual,
                         const char *const operation) const
      {
        if (!actual || !actual->compatible_with(*layout_))
          throw std::invalid_argument(std::string(operation) +
                                      " has an incompatible layout");
      }

      LayoutPtr          layout_;
      const NativeMetric &native_metric_;
      std::string        id_;
    };

    class ProblemBBinding final
    {
    public:
      using Backend     = nmopt::dealii_backend::SerialBackend;
      using Application = ::Step4<2>;
      using Problem     = ProblemB<2, Application>;
      using Primal      = nmopt::contract::PrimalBlockT<Backend>;
      using Covector    = nmopt::contract::CovectorBlockT<Backend>;
      using Model       = nmopt::contract::CallbackExecutableModelT<Backend>;
      using Partition   = nmopt::contract::StateControlPartitionT<Backend>;
      using Solvers     = nmopt::contract::StateAdjointSolversT<Backend>;
      using Reduced     = nmopt::contract::ReducedDTOT<Backend>;
      using SolveResult = nmopt::contract::FormulationSolveResultT<Backend>;

      explicit ProblemBBinding(Problem &problem)
        : problem_(problem)
        , native_metric_(problem_.mass())
        , variable_layout_(std::make_shared<const nmopt::contract::BlockLayout>(
            "external_step4_minimal_problem_b_variables",
            std::vector<nmopt::contract::SpaceId>{{"state"}, {"control"}},
            std::vector<std::size_t>{problem_.state_dimension(),
                                     problem_.control_dimension()}))
        , test_layout_(std::make_shared<const nmopt::contract::BlockLayout>(
            "external_step4_minimal_problem_b_test",
            std::vector<nmopt::contract::SpaceId>{{"state_test"}},
            std::vector<std::size_t>{problem_.state_dimension()}))
        , model_(make_model(problem_, variable_layout_, test_layout_))
        , partition_(model_, 0, 1)
        , solvers_(make_solvers(problem_, partition_.state_layout(),
                                test_layout_))
        , metric_(partition_.control_layout(), native_metric_)
        , reduced_(model_, partition_, solvers_)
      {}

      ProblemBBinding(const ProblemBBinding &) = delete;
      ProblemBBinding &operator=(const ProblemBBinding &) = delete;
      ProblemBBinding(ProblemBBinding &&) = delete;
      ProblemBBinding &operator=(ProblemBBinding &&) = delete;

      const Problem &
      problem() const
      {
        return problem_;
      }

      const LayoutPtr &
      control_layout() const
      {
        return partition_.control_layout();
      }

      const Model &
      model() const
      {
        return model_;
      }

      const ProblemBMassMetric &
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
      static nmopt::contract::LinearSolveReport
      solve_report(const Problem::SolveEvidence &evidence)
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
      make_model(Problem &problem,
                 const LayoutPtr &variable_layout,
                 const LayoutPtr &test_layout)
      {
        Problem *const problem_ptr = &problem;

        return Model(
          variable_layout,
          test_layout,
          [problem_ptr, test_layout](const Primal &variables) {
            auto value =
              problem_ptr->residual(variables.block(0), variables.block(1));
            return Covector(test_layout, {std::move(value)});
          },
          [problem_ptr, test_layout](const Primal &,
                                     const Primal &variable_tangent) {
            auto value = problem_ptr->residual_jvp(
              variable_tangent.block(0), variable_tangent.block(1));
            return Covector(test_layout, {std::move(value)});
          },
          [problem_ptr, variable_layout](const Primal &,
                                         const Primal &test_seed) {
            auto value = problem_ptr->residual_vjp(test_seed.block(0));
            return Covector(variable_layout,
                            {std::move(value.state),
                             std::move(value.control)});
          },
          [problem_ptr](const Primal &variables) {
            return problem_ptr->objective(variables.block(0),
                                          variables.block(1));
          },
          [problem_ptr, variable_layout](const Primal &variables) {
            auto value = problem_ptr->objective_derivative(
              variables.block(0), variables.block(1));
            return Covector(variable_layout,
                            {std::move(value.state),
                             std::move(value.control)});
          });
      }

      static Solvers
      make_solvers(Problem &problem,
                   const LayoutPtr &state_layout,
                   const LayoutPtr &test_layout)
      {
        Problem *const problem_ptr = &problem;
        Solvers       solvers;
        solvers.solve_state = [problem_ptr, state_layout](const Primal &control) {
          const auto result = problem_ptr->solve_state(control.block(0));
          return SolveResult(
            Primal(state_layout, {result.solution}),
            solve_report(result.evidence));
        };
        solvers.solve_adjoint = [problem_ptr, test_layout](
                                  const Primal &,
                                  const Covector &state_rhs) {
          const auto result =
            problem_ptr->solve_adjoint(state_rhs.block(0));
          return SolveResult(
            Primal(test_layout, {result.solution}),
            solve_report(result.evidence));
        };
        return solvers;
      }

      ProblemB<2, Application> &problem_;
      ProblemBMetric<2>         native_metric_;
      LayoutPtr                variable_layout_;
      LayoutPtr                test_layout_;
      Model                    model_;
      Partition                partition_;
      Solvers                  solvers_;
      ProblemBMassMetric       metric_;
      Reduced                  reduced_;
    };
  } // namespace minimal
} // namespace external_dealii_step4
