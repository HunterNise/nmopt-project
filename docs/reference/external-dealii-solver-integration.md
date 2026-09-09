# External deal.II solver integration reference

This guide is for an application that already owns a deal.II discretization
and wants to use nmopt's formulation and optimization services. The
application does not need to become a semantic compiler recipe and does not
need to transfer mesh, assembly, or native output ownership to nmopt.

The tested reference for the construction below is the
[external deal.II contract test](../../tests/application/external_application_dealii_contract.cc).
The test uses a small Poisson-control application, but the boundary is about
operations and ownership rather than that PDE family.

This reference describes the currently implemented API exercised by that
fixture. Authentic external-application adaptation cost and ergonomic
sufficiency are under evaluation in the
[external Step-4 boundary evaluation](../planning/external-dealii-boundary-evaluation.md).
That controlled native-versus-current-nmopt comparison retains the completed
PDE–solver boundary refactor as a baseline. It has not yet established whether
helpers, a boundary change, or deeper architectural work are warranted.

## Integration boundary

An external application owns the numerical realization:

- mesh, finite elements, DoF handlers, constraints, and assembled operators;
- state and adjoint linear solves, including their policies and diagnostics;
- the residual, its JVP and VJP, the objective, and its derivative;
- native reconstruction, field inspection, and filesystem output; and
- the lifetime of every object captured by a callback.

nmopt consumes those operations through small solver-facing contracts:

```text
existing deal.II application
    |  layouts, callbacks, solve services, metric
    v
ExecutableModelT + StateAdjointSolversT + MetricT
    |
    v
ReducedDTOT -> ReducedGradientSolverT (or another compatible optimizer)
```

`CompiledProblemT` is the packaging for the semantic/compiler path. It is not
required by this integration path. An external application should construct
the solver-facing contracts directly.

## 1. Keep the PDE application standalone

An application can continue to expose ordinary domain-specific methods. The
reference application has methods with this shape:

```cpp
using Vector = dealii::Vector<double>;

Vector residual(const Vector &state, const Vector &control) const;
Vector residual_jvp(const Vector &state_tangent,
                    const Vector &control_tangent) const;
Vector residual_vjp(const Vector &test_seed) const;
double objective(const Vector &state, const Vector &control) const;
ObjectiveDerivative objective_derivative(const Vector &state,
                                         const Vector &control) const;

Vector solve_state(const Vector &control) const;
Vector solve_adjoint(const Vector &state_objective_derivative) const;

void write_native_output(const std::filesystem::path &directory,
                         const Vector &state) const;
```

These methods are illustrative of the tested boundary: the application may
use different names, several assembled operators, iterative solves, or richer
native output. PDE assembly and output stay in the application.

Construction and runtime requirements differ in the current API.
[`CallbackExecutableModelT`](../../include/nmopt/contract/callback_executable_model.hpp)
requires all five executable operations: residual, residual JVP, residual VJP,
objective, and objective derivative. The selected first-order
[`ReducedDTOT` evaluation path](../../include/nmopt/contract/reduced_dto.hpp)
does not call residual or JVP at runtime. It requests the full residual VJP,
then consumes only its control component when forming the reduced derivative.
These are current contract facts; their integration cost is part of the
evaluation.

The snippets below assume `Application` names the application-owned class and
that `application` is a live instance of it.

## 2. Define compatible layouts

The current reduced DTO boundary uses one state block, one control block, and
one residual-test block. The layouts describe the block identity and dimension
that every callback must preserve.

```cpp
#include "nmopt/contract/callback_executable_model.hpp"
#include "nmopt/contract/reduced_dto.hpp"
#include "nmopt/dealii/mass_metric.hpp"
#include "nmopt/dealii/serial_backend.hpp"
#include "nmopt/solvers/reduced_gradient.hpp"

#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

using Backend = nmopt::dealii_backend::SerialBackend;
using Primal = nmopt::contract::PrimalBlockT<Backend>;
using Covector = nmopt::contract::CovectorBlockT<Backend>;
using LayoutPtr = nmopt::contract::LayoutPtr;

LayoutPtr variable_layout =
  std::make_shared<const nmopt::contract::BlockLayout>(
    "external_poisson_variables",
    std::vector<nmopt::contract::SpaceId>{{"state"}, {"control"}},
    std::vector<std::size_t>{application.state_dimension(),
                             application.control_dimension()});

LayoutPtr test_layout =
  std::make_shared<const nmopt::contract::BlockLayout>(
    "external_poisson_test",
    std::vector<nmopt::contract::SpaceId>{{"state_test"}},
    std::vector<std::size_t>{application.state_dimension()});
```

[`BlockLayout::compatible_with()`](../../include/nmopt/contract/layout.hpp)
compares ordered space IDs and dimensions. Distinct layout objects, including
objects with different display labels, can be compatible; pointer identity is
not required. Equal raw vector dimensions alone are insufficient. The callback
model rejects incompatible input or output layouts.

## 3. Expose the executable operations with callbacks

`CallbackExecutableModelT` checks input and output layouts at the solver
boundary. The following is the complete callback pattern used by the tested
external application; replace the application-specific method bodies with the
corresponding operations of the existing PDE code.

```cpp
using Model = nmopt::contract::CallbackExecutableModelT<Backend>;

Model make_model(Application &application,
                 const LayoutPtr &variable_layout,
                 const LayoutPtr &test_layout)
{
  Application *const application_ptr = &application;

  return Model(
    variable_layout,
    test_layout,
    [application_ptr, test_layout](const Primal &variables) {
      auto value = application_ptr->residual(variables.block(0),
                                              variables.block(1));
      return Covector(test_layout, {std::move(value)});
    },
    [application_ptr, test_layout](const Primal &,
                                   const Primal &variable_tangent) {
      auto value = application_ptr->residual_jvp(
        variable_tangent.block(0), variable_tangent.block(1));
      return Covector(test_layout, {std::move(value)});
    },
    [application_ptr, variable_layout](const Primal &,
                                       const Primal &test_seed) {
      auto packed = application_ptr->residual_vjp(test_seed.block(0));
      Application::Vector state(application_ptr->state_dimension());
      Application::Vector control(
        application_ptr->control_dimension());
      for (std::size_t index = 0; index < state.size(); ++index)
        state[index] = packed[index];
      for (std::size_t index = 0; index < control.size(); ++index)
        control[index] =
          packed[application_ptr->state_dimension() + index];
      return Covector(variable_layout,
                      {std::move(state), std::move(control)});
    },
    [application_ptr](const Primal &variables) {
      return application_ptr->objective(variables.block(0),
                                        variables.block(1));
    },
    [application_ptr, variable_layout](const Primal &variables) {
      const auto derivative = application_ptr->objective_derivative(
        variables.block(0), variables.block(1));
      return Covector(variable_layout,
                      {std::move(derivative.state),
                       std::move(derivative.control)});
    });
}
```

The residual VJP returns a covector in the full variable layout, while its
seed is a primal in the test layout. This distinction matters for the reduced
derivative and should be explicit in the application adapter. If the PDE
operator is nonlinear, the JVP callback uses both the point and tangent; the
linear reference problem ignores the point because its derivative is constant.

## 4. Supply state and adjoint solves

The reduced DTO asks the application for a state at a control and for an
adjoint at a full state/control point and state-objective covector. Wrap the
native vectors in the appropriate one-block layouts and return a
`FormulationSolveResultT`. Its
[one-argument constructor](../../include/nmopt/contract/linear_solve.hpp)
creates a converged "caller-supplied exact solve" report. It does not inspect
or validate a native solve. The snippets below use that convenience form as
the direct-solve reference fixture does; successful iterative solves should
use the two-argument form with an actual `LinearSolveReport`.

```cpp
using Partition = nmopt::contract::StateControlPartitionT<Backend>;
using Solvers = nmopt::contract::StateAdjointSolversT<Backend>;

Partition partition(model, 0, 1);
Solvers solvers;
Application *const application_ptr = &application;

solvers.solve_state = [application_ptr, state_layout = partition.state_layout()](
                        const Primal &control) {
  auto state = application_ptr->solve_state(control.block(0));
  return nmopt::contract::FormulationSolveResultT<Backend>(
    Primal(state_layout, {std::move(state)}));
};

solvers.solve_adjoint = [application_ptr,
                         test_layout](const Primal &,
                                      const Covector &state_rhs) {
  auto adjoint = application_ptr->solve_adjoint(state_rhs.block(0));
  return nmopt::contract::FormulationSolveResultT<Backend>(
    Primal(test_layout, {std::move(adjoint)}));
};
```

If this setup is placed in a helper, return the populated `Solvers` value and
keep the same `state_layout` and `test_layout` objects associated with the
callbacks:

```cpp
Solvers make_solvers(Application &application,
                     const LayoutPtr &state_layout,
                     const LayoutPtr &test_layout)
{
  Application *const application_ptr = &application;
  Solvers solvers;
  solvers.solve_state = [application_ptr, state_layout](
                          const Primal &control) {
    auto state = application_ptr->solve_state(control.block(0));
    return nmopt::contract::FormulationSolveResultT<Backend>(
      Primal(state_layout, {std::move(state)}));
  };
  solvers.solve_adjoint = [application_ptr, test_layout](
                            const Primal &,
                            const Covector &state_rhs) {
    auto adjoint = application_ptr->solve_adjoint(state_rhs.block(0));
    return nmopt::contract::FormulationSolveResultT<Backend>(
      Primal(test_layout, {std::move(adjoint)}));
  };
  return solvers;
}
```

The full point is supplied to `solve_adjoint` because a nonlinear or
point-dependent application may need it. The reference linear application
does not use it. A solve callback must not return a state or adjoint with an
incompatible layout or label a failed solve converged. A non-converged report
is rejected by `ReducedDTOT`; a native exception also propagates out of the
evaluation. Under the frozen protocol, the Step-4 evaluation will preserve
native solve exception propagation and supply actual convergence evidence
for successful CG solves. It will not label them exact using the convenience
constructor.

## 5. Add the metric and optional capabilities

The metric is an application-owned numerical object. For a deal.II mass
matrix, the supplied matrix must be square, dimensionally compatible with the
control layout, and positive definite for the selected solve policy.

```cpp
nmopt::dealii_backend::MassMetric metric(
  "external_poisson_control_l2",
  partition.control_layout(),
  application.control_mass_matrix());
```

The metric identifies the reduced derivative with a primal search direction.
It is not the objective regularization merely because the application may use
the same matrix in both places.

An application may also provide a `ConstraintT<Backend>` when the chosen
optimizer needs feasibility or projection, and a `ReducedHessianT<Backend>`
when a second-order method explicitly requests Hessian actions. Neither is
required for the first-order reduced integration described here. Use the
corresponding optimizer constructor only when that capability is actually
available, and keep its ownership in the application.

## 6. Construct the reduced DTO and optimizer

Construct the reduced DTO from the executable model, state/control partition,
and solve callbacks. Then pass it, the metric, and solver parameters to the
reduced optimizer.

```cpp
using Reduced = nmopt::contract::ReducedDTOT<Backend>;

Model model = make_model(application, variable_layout, test_layout);
Partition partition(model, 0, 1);
Solvers solvers = make_solvers(application,
                               partition.state_layout(),
                               test_layout);
nmopt::dealii_backend::MassMetric metric(
  "external_poisson_control_l2",
  partition.control_layout(),
  application.control_mass_matrix());
Reduced reduced(model, partition, solvers);

nmopt::solvers::ReducedSolverParameters parameters;
parameters.maximum_iterations = 100;
parameters.maximum_line_search_trials = 25;
parameters.gradient_tolerance = 1.0e-8;

nmopt::solvers::ReducedGradientSolverT<Backend> solver(
  reduced, metric, parameters);
const Primal initial_control(
  partition.control_layout(),
  {Backend::zeros(partition.control_layout()->dimension(0))});
const auto report = solver.solve(initial_control);
```

The initial control must use the metric/control layout and, when a constraint
is supplied, must already be feasible. The returned report contains the
optimization history and stopping reason. A caller can evaluate the final
control with `reduced.evaluate(report.control)` when it needs the final state,
adjoint, objective, or reduced derivative; use the actual result member name
provided by the selected solver result type in application code.

For a direct, reference-taking `ReducedDTOT`, keep `application`, layouts,
`model`, `partition`, `solvers`, and `metric` alive for the complete lifetime
of the DTO and optimizer. The reference test stores these dependencies in one
owning adapter so destruction order is unambiguous. A compiled product may
instead use the DTO's owning constructor, but that is a compiler-path detail,
not a requirement for external integration.

## 7. Keep reconstruction and output in the application

The reduced evaluation exposes the state in the solver-facing state layout.
The external application remains responsible for converting that state into
physical fields, diagnostics, visualization data, or checkpoint files.

```cpp
const auto evaluation = reduced.evaluate(initial_control);

application.write_native_output(
  output_directory,
  evaluation.state.block(0));
```

For an application with independent state coordinates, perform the physical
reconstruction before output or observation. For Dirichlet data, for example,
the application may need to apply its lifting and constraints before passing a
field to deal.II's `DataOut`. Do not add filesystem callbacks to
`ExecutableModelT` just to make native output accessible to an optimizer.

## 8. Verify derivatives and solve behavior

Before trusting an optimization result, verify the same operations that the
reduced DTO composes.

1. Check a state solve by evaluating the native residual at the returned
   state and control.
2. Check the residual JVP against a centered finite difference of the
   residual.
3. Check the residual JVP/VJP pairing:

   ```cpp
   const auto jvp = model.residual_jvp(point, tangent);
   const auto vjp = model.residual_vjp(point, test_seed);
   const double left = nmopt::contract::pair(jvp, test_seed);
   const double right = nmopt::contract::pair(vjp, tangent);
   ```

4. Check the reduced derivative against a centered finite difference of
   `reduced.evaluate(control).objective_value`.
5. Check a reduced Taylor remainder at a small step and inspect the state and
   adjoint solve reports.
6. Run a small optimization and verify that the stopping reason, objective
   history, and final gradient norm satisfy the application tolerance.

The external contract test exercises all six checks, including replacement of
the objective callbacks without changing the PDE callbacks and application-
owned native output. Use tolerances appropriate to the assembled operators and
linear-solve accuracy; a passing finite-difference check is evidence about the
specific discretization and callback wiring, not a substitute for the
mathematical model's own validation.

## Related contracts

- [PDE, formulation, and solver boundary](../design/pde-solver-boundary.md)
- [V0 executable contract](../implementation/v0/executable-contract.md)
- [Application assembly API](../reference/application-api.md)
- [External application contract test](../../tests/application/external_application_dealii_contract.cc)
- [Standalone external forward application](../../tests/dealii/external_poisson_forward.cc)
- [Current external Step-4 boundary evaluation](../planning/external-dealii-boundary-evaluation.md)
- [Superseded external tutorial roadmap](../planning/external-dealii-tutorial-roadmap.md)
