#pragma once

#include "instrumentation.hpp"
#include "problem_a.hpp"

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
  using LayoutPtr = nmopt::contract::LayoutPtr;

  class IdentityMetric final
    : public nmopt::contract::MetricT<nmopt::dealii_backend::SerialBackend>
  {
  public:
    using Backend  = nmopt::dealii_backend::SerialBackend;
    using Primal   = nmopt::contract::PrimalBlockT<Backend>;
    using Covector = nmopt::contract::CovectorBlockT<Backend>;

    explicit IdentityMetric(LayoutPtr layout)
      : id_("external_step4_identity")
      , layout_(std::move(layout))
    {
      if (!layout_)
        throw std::invalid_argument("Step-4 identity metric needs a layout");
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
                         "Step-4 identity metric primal");
      return Covector(layout_, {primal.block(0)});
    }

    Primal
    inverse_apply(const Covector &covector) const override
    {
      require_compatible(covector.layout(),
                         "Step-4 identity metric covector");
      return Primal(layout_, {covector.block(0)});
    }

  private:
    void
    require_compatible(const LayoutPtr &actual, const char *operation) const
    {
      if (!actual || !actual->compatible_with(*layout_))
        throw std::invalid_argument(std::string(operation) +
                                    " has an incompatible layout");
    }

    std::string id_;
    LayoutPtr   layout_;
  };

  class NmoptBinding final
  {
  public:
    using Backend  = nmopt::dealii_backend::SerialBackend;
    using Primal   = nmopt::contract::PrimalBlockT<Backend>;
    using Covector = nmopt::contract::CovectorBlockT<Backend>;
    using Model = nmopt::contract::CallbackExecutableModelT<Backend>;
    using Partition = nmopt::contract::StateControlPartitionT<Backend>;
    using Solvers = nmopt::contract::StateAdjointSolversT<Backend>;
    using Reduced = nmopt::contract::ReducedDTOT<Backend>;
    using SolveResult = nmopt::contract::FormulationSolveResultT<Backend>;

    explicit NmoptBinding(Instrumentation &instrumentation)
      : instrumentation_(instrumentation)
      , problem_(instrumentation)
      , variable_layout_(std::make_shared<const nmopt::contract::BlockLayout>(
          "external_step4_problem_a_variables",
          std::vector<nmopt::contract::SpaceId>{{"state"}, {"control"}},
          std::vector<std::size_t>{problem_.state_dimension(),
                                   problem_.control_dimension()}))
      , test_layout_(std::make_shared<const nmopt::contract::BlockLayout>(
          "external_step4_problem_a_test",
          std::vector<nmopt::contract::SpaceId>{{"state_test"}},
          std::vector<std::size_t>{problem_.state_dimension()}))
      , model_(make_model(problem_, variable_layout_, test_layout_))
      , partition_(model_, 0, 1)
      , solvers_(make_solvers(problem_, partition_.state_layout(), test_layout_))
      , metric_(partition_.control_layout())
      , reduced_(model_, partition_, solvers_)
    {}

    const ProblemA &
    problem() const
    {
      return problem_;
    }

    const LayoutPtr &
    variable_layout() const
    {
      return variable_layout_;
    }

    const LayoutPtr &
    test_layout() const
    {
      return test_layout_;
    }

    const LayoutPtr &
    control_layout() const
    {
      return partition_.control_layout();
    }

    const LayoutPtr &
    state_layout() const
    {
      return partition_.state_layout();
    }

    const Model &
    model() const
    {
      return model_;
    }

    const IdentityMetric &
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
    solve_report(const ProblemA::SolveEvidence &evidence)
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
    make_model(ProblemA &problem,
               const LayoutPtr &variable_layout,
               const LayoutPtr &test_layout)
    {
      ProblemA *const problem_ptr = &problem;

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
                          {std::move(value.state), std::move(value.control)});
        },
        [problem_ptr](const Primal &variables) {
          return problem_ptr->objective(variables.block(0), variables.block(1));
        },
        [problem_ptr, variable_layout](const Primal &variables) {
          auto value = problem_ptr->objective_derivative(variables.block(0),
                                                         variables.block(1));
          return Covector(variable_layout,
                          {std::move(value.state), std::move(value.control)});
        });
    }

    static Solvers
    make_solvers(ProblemA &problem,
                 const LayoutPtr &state_layout,
                 const LayoutPtr &test_layout)
    {
      ProblemA *const problem_ptr = &problem;
      Solvers       solvers;
      solvers.solve_state = [problem_ptr, state_layout](const Primal &control) {
        auto result = problem_ptr->solve_state(control.block(0));
        return SolveResult(
          Primal(state_layout, {std::move(result.solution)}),
          solve_report(result.evidence));
      };
      solvers.solve_adjoint = [problem_ptr, test_layout](
                                const Primal &,
                                const Covector &state_rhs) {
        auto result = problem_ptr->solve_adjoint(state_rhs.block(0));
        return SolveResult(
          Primal(test_layout, {std::move(result.solution)}),
          solve_report(result.evidence));
      };
      return solvers;
    }

    Instrumentation & instrumentation_;
    ProblemA          problem_;
    LayoutPtr         variable_layout_;
    LayoutPtr         test_layout_;
    Model             model_;
    Partition         partition_;
    Solvers           solvers_;
    IdentityMetric    metric_;
    Reduced           reduced_;
  };
} // namespace external_dealii_step4
