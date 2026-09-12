# External deal.II solver integration reference

This reference describes the implemented public API for connecting an existing
deal.II application to nmopt's reduced formulation and optimizer. It requires
no `ProblemSpec`, compiler product, recipe, manifest, or project runner.
The application retains its numerical realization and output ownership.

The current reduced construction supports one state block, one control block,
and one residual-test block. The tested deal.II backend is serial.

The runnable examples are the Step-4 [A consumer](../../apps/external-dealii/step-4/minimal/problem_a.cc)
and [B consumer](../../apps/external-dealii/step-4/minimal/problem_b.cc), with
separate [A](../../apps/external-dealii/step-4/minimal/problem_a_binding.hpp)
and [B](../../apps/external-dealii/step-4/minimal/problem_b_binding.hpp)
bindings. Their classes are application-local examples, not public nmopt
factory APIs. Read the [explanatory overview](../../apps/external-dealii/step-4/external-integration-overview.md)
for the mathematical construction and the [closure audit](../planning/review/external-dealii-boundary-evaluation/closure-report.md)
for evidence and limits. This reference records existing contracts, not a new
interface or an acceptance claim for arbitrary external applications.

## 1. Required objects and ownership

```text
application-native mathematics and solve services
    -> CallbackExecutableModelT + StateAdjointSolversT
    -> ReducedDTOT
    -> ReducedSearchSolverT, supplied with MetricT
    -> retained numerical result -> application output
```

| Object | Public header | Role |
| --- | --- | --- |
| `SerialBackend` | [serial_backend.hpp](../../include/nmopt/dealii/serial_backend.hpp) | Native serial deal.II vectors and primitive algebra |
| `BlockLayout`, `PrimalBlockT`, `CovectorBlockT` | [layout.hpp](../../include/nmopt/contract/layout.hpp) | Space identities, dimensions, and typed coefficient values |
| `CallbackExecutableModelT` | [callback_executable_model.hpp](../../include/nmopt/contract/callback_executable_model.hpp) | Five mathematical callbacks with layout checks |
| `StateControlPartitionT`, `StateAdjointSolversT`, `ReducedDTOT` | [reduced_dto.hpp](../../include/nmopt/contract/reduced_dto.hpp) | State/control partition, native solve callbacks, and reduced orchestration |
| `LinearSolveReport`, `FormulationSolveResultT` | [linear_solve.hpp](../../include/nmopt/contract/linear_solve.hpp) | A solved primal value and actual convergence evidence |
| `MetricT`, `ConstraintT` | [metric_constraint.hpp](../../include/nmopt/contract/metric_constraint.hpp) | Primal/dual conversion and optional compatible projection |
| `ReducedSearchSolverT` | [reduced_gradient.hpp](../../include/nmopt/solvers/reduced_gradient.hpp) | Search-direction and line-search orchestration |
| `ReducedSolverParameters`, `ReducedSolverResultT` | [reduced_search.hpp](../../include/nmopt/solvers/reduced_search.hpp) | Policy settings, histories, stopping reason, and retained result |

The application owns mesh, FE, DoF handler, constraints, matrices, control
mathematics, objective, derivatives, native solves, and reconstruction.
The binding maps those operations into contracts. nmopt schedules their use;
it does not infer the OCP from a forward solver or require a matrix view at
its generic formulation boundary.

## 2. Backend, blocks, and layouts

The current serial deal.II path uses:

```cpp
using Backend = nmopt::dealii_backend::SerialBackend;
using Primal = nmopt::contract::PrimalBlockT<Backend>;
using Covector = nmopt::contract::CovectorBlockT<Backend>;
using LayoutPtr = nmopt::contract::LayoutPtr;
```

`Backend::Vector` is `dealii::Vector<double>`. The backend supplies `zeros`,
`size`, `value`, `set_value`, `dot`, `add_scaled`, and `scale`, including native
size-range checks. It does not define FE assembly, boundary conditions, or
PDE inversion. Distributed vectors and their ownership/ghost policy are a
separate extension, not part of this tested backend.

The reduced partition currently requires exactly two variable blocks and one
residual-test block. State and control must select distinct block indices.
For example, with application-supplied dimensions `n_state` and `n_control`:

```cpp
auto variable_layout =
  std::make_shared<const nmopt::contract::BlockLayout>(
    "external_variables",
    std::vector<nmopt::contract::SpaceId>{{"state"}, {"control"}},
    std::vector<std::size_t>{n_state, n_control});
auto test_layout =
  std::make_shared<const nmopt::contract::BlockLayout>(
    "external_test",
    std::vector<nmopt::contract::SpaceId>{{"state_test"}},
    std::vector<std::size_t>{n_state});
```

These are setup fragments; the linked bindings contain complete includes and
construction. A uses 289 state/control entries. B uses 225 free state/test
entries and 289 full control entries. Choose spaces and dimensions from the
new application's actual discretization rather than copying these counts.

`BlockLayout::compatible_with()` compares ordered `SpaceId` values and
block dimensions. Display labels and pointer identity do not determine
compatibility. Equal raw dimensions alone are insufficient. The layout does
not independently prove that two FE bases have the same mathematical meaning;
the application must assign identities consistently.

Block constructors accept a layout and `std::vector<Backend::Vector>`, check
block counts/dimensions, and own the coefficient vectors. `block(i)` returns a
const native-vector reference. Wrapping existing lvalues can copy vectors;
this interface does not promise zero-copy storage or use deal.II's
`BlockVector` as its container.

A primal represents a point, tangent, or adjoint test vector. A covector
represents a derivative or residual functional. `nmopt::contract::pair`
evaluates the coefficient dual pairing for compatible covector/primal values.
It does not silently apply a mass matrix or another metric.

## 3. Executable callbacks

With `Model = nmopt::contract::CallbackExecutableModelT<Backend>`, the
constructor takes `variable_layout`, `test_layout`, and these five callbacks
in order. Let $x=(y,u)$ denote the full point and $q$ a test-space primal.

| Callback | C++ function shape | Mathematical result and output layout |
| --- | --- | --- |
| Residual | `Covector(const Primal &x)` | $E(x)$ in the test layout |
| Residual JVP | `Covector(const Primal &x, const Primal &v)` | $E'(x)v$ in the test layout |
| Residual VJP | `Covector(const Primal &x, const Primal &q)` | $E'(x)^{\ast}q$ in the full variable layout |
| Objective | `double(const Primal &x)` | $J(x)$ |
| Objective derivative | `Covector(const Primal &x)` | $J'(x)$ in the full variable layout |

All five callbacks must be nonempty. Inputs and returned block layouts are
checked. Residual and JVP remain construction requirements even though the
selected successful first-order reduced optimization path calls neither.
Full VJP must return a valid state and control contribution; it cannot use a
zero placeholder for a mathematically nonzero component simply because that
component is discarded by the current reduced derivative.

For example, B's VJP wrapping is:

```cpp
[problem_ptr, variable_layout](const Primal &, const Primal &test_seed) {
  auto value = problem_ptr->residual_vjp(test_seed.block(0));
  return Covector(variable_layout,
                  {std::move(value.state), std::move(value.control)});
}
```

This linear example ignores the point because its derivative is constant.
A nonlinear or point-dependent application must use it. The native method
name and returned struct are B's choices, not requirements on application
classes. Lambdas may call any equivalent native operations.

## 4. State and adjoint solve services

After constructing `model`, form:

```cpp
using Partition = nmopt::contract::StateControlPartitionT<Backend>;
using Solvers = nmopt::contract::StateAdjointSolversT<Backend>;
using SolveResult = nmopt::contract::FormulationSolveResultT<Backend>;
Partition partition(model, 0, 1);
Solvers solvers;
```

| Callback member | Inputs | Required result |
| --- | --- | --- |
| `solve_state` | Control `const Primal &` | `SolveResult` containing the state in `partition.state_layout()` |
| `solve_adjoint` | Full point `const Primal &`, state-objective `const Covector &` | `SolveResult` containing a primal in the model's test layout |

The mathematical convention is:

```math
E(y,u)=0,\qquad E_{y}(y,u)^{\ast}p=J_{y}(y,u),\qquad
j'(u)=J_{u}(y,u)-E_{u}(y,u)^{\ast}p.
```

The adjoint RHS is a state covector, while the returned adjoint is a test
primal. The full point lets the callback update a point-dependent transpose
operator. Reusing the state solve for the adjoint requires actual mathematical
justification; symmetry of Step-4 does not establish it for other PDEs.

`FormulationSolveResultT` has two constructors:

- `(Primal solution, LinearSolveReport report)` preserves supplied evidence.
- `(Primal solution)` declares a converged `caller-supplied exact solve`.
  This convenience constructor does not inspect a native solver and must not
  be used to manufacture iterative convergence.

`LinearSolveReport` contains `algorithm`, `preconditioner`,
`maximum_iterations`, `iterations`, `relative_tolerance`, `absolute_tolerance`,
`requested_tolerance`, `achieved_residual`, and `termination`.
`termination` defaults to `LinearSolveTermination::failed`;
`converged()` tests whether it is `LinearSolveTermination::converged`.

The Step-4 wrappers report CG, identity preconditioning, maximum 1000
iterations, absolute/requested tolerance $10^{-12}$, and the native monitored
iterations/residual. Other applications must translate their actual policies
and evidence. The wrapper does not need an extra residual matrix action just
to populate the report.

A typical state wrapper, with an application-local report translation, is:

```cpp
solvers.solve_state = [problem_ptr, state_layout = partition.state_layout()](
                        const Primal &control) {
  const auto result = problem_ptr->solve_state(control.block(0));
  return SolveResult(Primal(state_layout, {result.solution}),
                     solve_report(result.evidence));
};
```

`solve_report` here denotes the explicit local conversion shown in the
[A](../../apps/external-dealii/step-4/minimal/problem_a_binding.hpp) and
[B](../../apps/external-dealii/step-4/minimal/problem_b_binding.hpp) bindings;
it is not a public nmopt helper. A non-converged report is rejected by the
reduced service. Native exceptions propagate; they are not converted into
rejected Armijo trials. Public convergence evidence is not an independent
proof that the application solved the correct equation.

## 5. Metrics and optional capabilities

A `MetricT<Backend>` implementation supplies:

```cpp
const std::string &id() const;
const LayoutPtr &layout() const;
Covector apply(const Primal &primal) const;
Primal inverse_apply(const Covector &covector) const;
```

These are virtual methods overridden by a concrete metric. The inherited
realization witness identifies a constructed metric for capability checks;
a display ID alone does not establish projection compatibility.

The metric maps a control primal to its dual and inverts that map:

```math
G:U\rightarrow U^{\ast},\qquad g=G^{-1}j'(u).
```

The metric is independent of objective regularization even when the same
matrix appears in both. A defines a local identity metric. B's local adapter
forwards to its native mass application/inverse, checks native CG convergence,
and returns a primal. `MetricT::inverse_apply()` returns no public solve
report. Instrumented metric records in the evaluation are separate support.

An application may instead choose an existing backend realization such as
[MassMetric](../../include/nmopt/dealii/mass_metric.hpp). Its constructor takes
an ID, a one-block layout, a `std::shared_ptr<const dealii::SparseMatrix<double>>`,
and optional solve parameters. The matrix must represent the appropriate
symmetric positive-definite metric, with compatible dimensions and valid
sparsity lifetime. Selecting that implementation selects its numerical policy;
it is not automatically equivalent to an application's existing inverse.
The unqualified `DiagonalMetric` uses `DenseBackend`, so it is not a drop-in
identity metric for `SerialBackend`.

A projected optimizer also needs a compatible `ConstraintT<Backend>` and an
initially feasible control. Metric projection support is checked explicitly;
coefficient clipping is not automatically valid for a general mass metric.
A Newton direction needs a supplied `ReducedHessianT<Backend>` through its
direction policy. Neither constraint nor Hessian capability is needed by the
unconstrained first-order A/B path. See the [solver source](../../include/nmopt/solvers/reduced_gradient.hpp)
and [direction policies](../../include/nmopt/solvers/reduced_search.hpp) for
implemented combinations; their existence does not extend the Step-4 evidence.

## 6. Reduced evaluation and reuse

The direct constructor is:

```cpp
nmopt::contract::ReducedDTOT<Backend> reduced(model, partition, solvers);
```

The model is borrowed; partition and solve-service values are stored by value.

| Operation | Work and result |
| --- | --- |
| `evaluate_value(control)` | State solve and objective; returns state, control, full point, objective, and state-solve report |
| `augment_derivative(value)` | Objective derivative, adjoint solve, and full residual VJP; reuses the supplied state and returns a reduced covector and solve reports |
| `evaluate(control)` | Composes the preceding two operations |
| `gradient_direction(derivative, metric)` | Applies the metric inverse; returns the gradient, without a descent minus sign |

Value evaluations carry a service token. Derivative augmentation rejects a
value from a different reduced service and checks consistency of its
state/control/full point. Treat returned values as coherent snapshots;
retention is not permission to mutate borrowed application mathematics
between stages. The API does not detect arbitrary changes to a mesh or
operator inside a captured application.

The reduced derivative uses the control part of the full objective derivative
minus the control part of the full residual VJP. Residual evaluation and JVP
are not used by these operations. Equation and derivative verification remain
separate tasks.

## 7. Optimizer, results, and output

The selected first-order constructor is:

```cpp
nmopt::solvers::ReducedSearchSolverT<Backend> solver(
  reduced, metric, parameters);
const auto result = solver.solve(initial_control);
```

Default template policies are steepest descent and Armijo.
`ReducedGradientSolverT<Backend>` is the compatibility alias for that
combination. `ReducedSolverParameters` defaults to 100 accepted iterations,
20 trials, and gradient tolerance $10^{-8}$. It also controls the stopping
criterion, optional relative/objective/step criteria and target, initial and
minimum step, Armijo fraction, and backtracking factor. Zero disables the
optional relative/objective/step tolerances.

Those defaults are not the frozen Step-4 policy: A/B explicitly use 5000
accepted iterations, 30 trials, gradient-norm stopping at $10^{-6}$, and the
remaining settings shown in their entry points. Applications choose their
own policy rather than treating the experiment's limits as universal defaults.

`solve()` returns `ReducedSolverResultT<Backend>` (also named
`ReducedGradientResultT<Backend>`). Important fields include:

- `control` and `final_evaluation`;
- `objective_history`, `gradient_norm_history`, and other stopping histories;
- `accepted_iterations`, `line_search_trial_count`, and work counts;
- `stopping_reason`, selected parameters, trial records, and accepted records.

A result is not synonymous with gradient convergence. Inspect
`stopping_reason`; iteration-limit and line-search-failure results are also
possible. `final_evaluation` contains `state`, `adjoint`, `full_point`,
`reduced_derivative`, `objective_value`, `state_solve`, and `adjoint_solve`.
It retains the final current evaluation, including the initial evaluation
when no step was accepted. Do not solve again merely to write that state.

For the B consumer, output is exactly:

```cpp
const auto full_state = problem.coordinates().reconstruct(
  result.final_evaluation.state.block(0));
application.output_results(full_state, output);
```

A can write its retained full state directly through `ProblemA`. The
application performs its own lifting, reconstruction, diagnostics, and file
output; the executable model has no filesystem callback. A fresh state solve
belongs to an independent audit when requested, not ordinary retained output.

## 8. Lifetime rules

For the reference-taking reduced constructor:

- The model must outlive the reduced service. Every object borrowed by a
  callback must outlive all calls to that callback.
- Layouts are shared values. The DTO stores its own partition and solve
  callbacks; original local `partition` and `solvers` variables need not be
  retained solely for the DTO.
- The solver borrows the reduced service and metric; both must outlive it.
  A metric adapter's native service and any matrices/sparsity it borrows must
  also remain alive.

A's `ProblemA` owns its prepared Step-4 instance and outlives its binding.
For B the lifetime chain is application → problem → binding → solver; the
binding's native metric outlives its metric adapter. Both local binding
classes delete copy/move operations because member references must remain
valid. That packaging is an example, not a universal requirement on all
external applications.

An overload accepts `std::shared_ptr<const ExecutableModelT<Backend>>`,
partition, solvers, and an optional `std::shared_ptr<const void>` lifetime
owner. It retains the model and supplied owner, but does not infer ownership
of arbitrary raw pointers captured by callbacks. Owning a model alone does
not automatically own the PDE application.

## 9. Verification and navigation

Verify native state equations, off-solution JVP finite differences, VJP
pairings, objective derivatives, and reduced derivatives before interpreting
optimizer convergence. In particular:

```cpp
const auto jvp = model.residual_jvp(point, tangent);
const auto vjp = model.residual_vjp(point, test_seed);
const double left = nmopt::contract::pair(jvp, test_seed);
const double right = nmopt::contract::pair(vjp, tangent);
```

Compare those pairings with problem-appropriate tolerances. Independent final
state/adjoint/stationarity checks catch errors that a self-consistent pair of
callbacks could share. The Step-4 tests also use dense oracles and compare
actual executable output; those evaluation facilities are not dependencies of
an external consumer.

- [Step-4 walkthrough](../../apps/external-dealii/step-4/external-integration-overview.md)
  explains how the native OCP and binding fit together.
- [Step-4 implementation report](../../apps/external-dealii/step-4/integration-report.md)
  maps mathematical operations to code and accounts for source responsibilities.
- [Minimal consumer commands](../../apps/external-dealii/step-4/minimal/README.md)
  build and run the tested examples from this repository. Package/install
  support is outside the evaluation.
- [A contract test](../../tests/application/external_step4_minimal_problem_a_contract.cc)
  and [B contract test](../../tests/application/external_step4_minimal_problem_b_contract.cc)
  exercise the exact minimal bindings and executables.
- [Earlier external Poisson fixture](../../tests/application/external_application_dealii_contract.cc)
  additionally demonstrates objective replacement and direct-solve wrapping.
- [PDE/formulation/solver design](../design/pde-solver-boundary.md) owns the
  architecture; [system blueprint](../design/system-blueprint.md) gives the
  wider project map, and [application assembly API](application-api.md)
  describes the separate semantic/compiler route.
