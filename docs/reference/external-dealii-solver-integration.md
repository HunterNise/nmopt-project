# External deal.II solver integration reference

This reference explains how to connect an **existing deal.II PDE code** to
nmopt without rebuilding that code as a semantic `ProblemSpec`.

Use this path when the application already owns its mesh, finite-element
spaces, assembly, linear solvers, boundary treatment, and output. The job of
the integration is to expose the numerical operations needed by an optimal
control formulation while leaving those application internals in place.

For the conceptual boundary, see
[Integrating an existing PDE application](../manual/concepts/13-integrating-an-existing-pde-application.md).
For reduced state/adjoint conventions, see
[Reduced state–adjoint formulation](../manual/concepts/05-reduced-state-adjoint-formulation.md).

The canonical code example is Step-4:

- [minimal consumers](../../apps/external-dealii/step-4/minimal/README.md);
- [application-side integration](../../apps/external-dealii/step-4/integration/).

By the end of this page, the intended application structure is:

```text
existing deal.II application
        │
        │ keep mesh / FE / assembly / native output here
        ▼
application-local OCP facade
        │
        │ residual, derivatives, state/adjoint solves
        ▼
nmopt binding
        │
        ├── layouts
        ├── CallbackExecutableModelT
        ├── StateControlPartitionT
        ├── MetricT
        └── ReducedDTOT
        │
        ▼
nmopt optimizer
        │
        ▼
application-owned reconstruction / output
```

## Choose this path instead of the compiler path

nmopt has two peer producer paths:

```text
semantic/compiler path                    existing-application path

ProblemSpec                               existing PDE application
    │                                             │
    ▼                                             ▼
DealiiCompiler                         application-local OCP facade
    │                                             │
    └─────────────────────┬───────────────────────┘
                          ▼
              common numerical/formulation
                       contracts
                          │
                          ▼
                      solvers
```

If nmopt should create the numerical realization from a semantic graph, use
[Problem authoring](problem-authoring.md) and [Compiler](compiler.md).

If the PDE code already exists and should remain authoritative for assembly
and solves, use the direct path here.

Do not create a dummy `ProblemSpec` only to reach the optimizer.
`CompiledProblemT` is a compiler-path product, not a universal wrapper for all
nmopt applications.

## First decide what the existing code must expose

Before writing nmopt adapters, isolate an application-local OCP facade around
the existing PDE code.

For a reduced problem with state $y$ and control $u$, the facade needs the
operations corresponding to

$$
E(y,u)=0,
\qquad
J(y,u),
$$

and their first derivatives.

A useful native shape is:

```cpp
class MyOptimalControlProblem
{
public:
  using Vector = dealii::Vector<double>;

  std::size_t state_dimension() const;
  std::size_t control_dimension() const;

  Vector residual(
    const Vector &state,
    const Vector &control) const;

  Vector residual_jvp(
    const Vector &state,
    const Vector &control,
    const Vector &state_tangent,
    const Vector &control_tangent) const;

  struct VjpResult
  {
    Vector state;
    Vector control;
  };

  VjpResult residual_vjp(
    const Vector &state,
    const Vector &control,
    const Vector &test_seed) const;

  double objective(
    const Vector &state,
    const Vector &control) const;

  VjpResult objective_derivative(
    const Vector &state,
    const Vector &control) const;

  SolveResult solve_state(const Vector &control);

  SolveResult solve_adjoint(
    const Vector &state,
    const Vector &control,
    const Vector &state_rhs);
};
```

Those exact method names are **application-local**, not required nmopt API.
The point is to put the OCP-facing operations in one narrow layer rather than
spread nmopt callbacks across assembly classes and executable code.

The Step-4 `integration/problem_a.hpp` and `integration/problem_b.hpp` are the
concrete examples of this facade.

## What the facade means mathematically

Let

$$
x=(y,u).
$$

The nmopt executable model expects:

```text
residual(x)                 E(y,u)
residual_jvp(x, v)          E'(y,u)v
residual_vjp(x, q)          E'(y,u)^* q
objective(x)                J(y,u)
objective_derivative(x)     J'(y,u)
```

The VJP result contains both state and control covector components:

$$
E'(y,u)^{\ast}q
=
\begin{bmatrix}
E_{y}(y,u)^{\ast}q\\
E_{u}(y,u)^{\ast}q
\end{bmatrix}.
$$

That control contribution is essential for the reduced derivative. It should
not be omitted merely because the native PDE solver was originally written
only to solve for the state.

For a nonlinear PDE, JVP/VJP must use the current point $(y,u)$. Step-4 can
omit that dependence in some callbacks only because its demonstrated
operators are linear in the relevant variables.

## Working path: connect one existing application

The following sequence mirrors Step-4 Problem A and shows where each piece
belongs.

### Step 1 – define nmopt layouts from native coordinate spaces

```cpp
using Backend =
  nmopt::dealii_backend::SerialBackend;

using Primal =
  nmopt::contract::PrimalBlockT<Backend>;

using Covector =
  nmopt::contract::CovectorBlockT<Backend>;

using LayoutPtr =
  nmopt::contract::LayoutPtr;
```

Build the full variable layout:

```cpp
auto variable_layout =
  std::make_shared<const nmopt::contract::BlockLayout>(
    "my_problem_variables",
    std::vector<nmopt::contract::SpaceId>{
      {"state"},
      {"control"}
    },
    std::vector<std::size_t>{
      problem.state_dimension(),
      problem.control_dimension()
    });
```

and the residual-test layout:

```cpp
auto test_layout =
  std::make_shared<const nmopt::contract::BlockLayout>(
    "my_problem_test",
    std::vector<nmopt::contract::SpaceId>{
      {"state_test"}
    },
    std::vector<std::size_t>{
      problem.state_dimension()
    });
```

The current reduced DTO requires exactly two variable blocks and one test
block.

`SpaceId` describes mathematical coordinate identity. Equal vector sizes do
not make two spaces interchangeable.

For an application with constrained/full and independent/free coordinates,
do not hide that difference by assigning the same `SpaceId`; expose the
independent coordinates to nmopt and keep reconstruction in the application
facade. Problem B below shows that case.

### Step 2 – wrap the five executable-model callbacks

```cpp
#include "nmopt/contract/callback_executable_model.hpp"

using Model =
  nmopt::contract::CallbackExecutableModelT<Backend>;

MyOptimalControlProblem *const problem_ptr = &problem;
```

Then construct:

```cpp
Model model(
  variable_layout,
  test_layout,

  // E(y,u)
  [problem_ptr, test_layout](const Primal &variables) {
    auto value = problem_ptr->residual(
      variables.block(0),
      variables.block(1));

    return Covector(
      test_layout,
      {std::move(value)});
  },

  // E'(y,u) v
  [problem_ptr, test_layout](
    const Primal &variables,
    const Primal &tangent) {

    auto value = problem_ptr->residual_jvp(
      variables.block(0),
      variables.block(1),
      tangent.block(0),
      tangent.block(1));

    return Covector(
      test_layout,
      {std::move(value)});
  },

  // E'(y,u)^* q
  [problem_ptr, variable_layout](
    const Primal &variables,
    const Primal &test_seed) {

    auto value = problem_ptr->residual_vjp(
      variables.block(0),
      variables.block(1),
      test_seed.block(0));

    return Covector(
      variable_layout,
      {
        std::move(value.state),
        std::move(value.control)
      });
  },

  // J(y,u)
  [problem_ptr](const Primal &variables) {
    return problem_ptr->objective(
      variables.block(0),
      variables.block(1));
  },

  // J'(y,u)
  [problem_ptr, variable_layout](const Primal &variables) {
    auto value = problem_ptr->objective_derivative(
      variables.block(0),
      variables.block(1));

    return Covector(
      variable_layout,
      {
        std::move(value.state),
        std::move(value.control)
      });
  });
```

The application facade may expose a simpler signature when an operator is
constant. For example, Step-4 Problem A does not need the current point for
its derivative actions. The nmopt callback still has the general signature;
the adapter is where the simplification belongs.

### Step 3 – expose the state/control partition

```cpp
using Partition =
  nmopt::contract::StateControlPartitionT<Backend>;

Partition partition(model, 0, 1);
```

This declares:

```text
variable block 0 -> state
variable block 1 -> control
```

Use the derived layouts downstream:

```cpp
partition.state_layout();
partition.control_layout();
```

rather than creating new nominally equivalent layouts for the solve and metric
services.

### Step 4 – adapt the native state solve

The reduced formulation needs a service solving

$$
E(y,u)=0
$$

for $y$ at a given control.

```cpp
using SolveResult =
  nmopt::contract::FormulationSolveResultT<Backend>;

using Solvers =
  nmopt::contract::StateAdjointSolversT<Backend>;

Solvers solvers;

solvers.solve_state =
  [problem_ptr,
   state_layout = partition.state_layout()]
  (const Primal &control) {

    const auto native =
      problem_ptr->solve_state(
        control.block(0));

    return SolveResult(
      Primal(
        state_layout,
        {native.solution}),
      solve_report(native.evidence));
  };
```

`solve_report(...)` is application code that translates native solver evidence
to `LinearSolveReport`.

A representative adapter is:

```cpp
nmopt::contract::LinearSolveReport
solve_report(const NativeSolveEvidence &evidence)
{
  return {
    "CG",
    "identity",
    1000,
    evidence.iterations,
    0.0,
    1.0e-12,
    1.0e-12,
    evidence.final_residual,
    evidence.converged
      ? nmopt::contract::LinearSolveTermination::converged
      : nmopt::contract::LinearSolveTermination::failed
  };
}
```

Report the algorithm and termination actually used by the native code. nmopt
uses this evidence to reject failed state/adjoint solves rather than silently
continuing optimization.

### Step 5 – adapt the native adjoint solve

nmopt uses the convention

$$
E_{y}(y,u)^{\ast}p=J_{y}(y,u),
$$

which leads to the reduced derivative

$$
j'(u)
=
J_{u}(y,u)-E_{u}(y,u)^{\ast}p.
$$

The service callback receives the current full point and the state-objective
covector:

```cpp
solvers.solve_adjoint =
  [problem_ptr, test_layout](
    const Primal &full_point,
    const Covector &state_rhs) {

    const auto native =
      problem_ptr->solve_adjoint(
        full_point.block(0),
        full_point.block(1),
        state_rhs.block(0));

    return SolveResult(
      Primal(
        test_layout,
        {native.solution}),
      solve_report(native.evidence));
  };
```

If the native application has a symmetric linear state operator, its adjoint
implementation may reuse the state matrix. That reuse belongs inside the
application facade. The nmopt binding should still expose a separate adjoint
service because a different application may need an actual transpose
assembly/solve.

### Step 6 – expose the control metric

The reduced derivative is a covector. Optimization needs a metric $G$ to form
its primal gradient:

$$
Gg=j'(u).
$$

A metric implements:

```cpp
const std::string &id() const;
const LayoutPtr &layout() const;

Covector apply(const Primal &) const;
Primal inverse_apply(const Covector &) const;
```

For a genuinely Euclidean/identity control geometry:

```cpp
class IdentityMetric final
  : public nmopt::contract::MetricT<Backend>
{
public:
  explicit IdentityMetric(LayoutPtr layout)
    : id_("my_identity_metric")
    , layout_(std::move(layout))
  {}

  const std::string &id() const override
  {
    return id_;
  }

  const LayoutPtr &layout() const override
  {
    return layout_;
  }

  Covector apply(const Primal &value) const override
  {
    require_compatible(value.layout());
    return Covector(layout_, {value.block(0)});
  }

  Primal inverse_apply(const Covector &value) const override
  {
    require_compatible(value.layout());
    return Primal(layout_, {value.block(0)});
  }

private:
  void require_compatible(const LayoutPtr &actual) const
  {
    if (!actual || !actual->compatible_with(*layout_))
      throw std::invalid_argument("incompatible metric layout");
  }

  std::string id_;
  LayoutPtr   layout_;
};
```

For finite-element $L^{2}$ geometry, the application normally exposes mass
application and mass inversion instead. Do not choose identity only because it
makes the adapter shorter.

The stock serial deal.II
[`MassMetric`](../../include/nmopt/dealii/mass_metric.hpp) can be used when the
application already has a compatible mass realization.

For the conceptual role of the metric, see
[Metrics, gradients, and constraints](../manual/concepts/04-metrics-gradients-and-constraints.md).

### Step 7 – compose `ReducedDTOT`

```cpp
using Reduced =
  nmopt::contract::ReducedDTOT<Backend>;

IdentityMetric metric(
  partition.control_layout());

Reduced reduced(
  model,
  partition,
  solvers);
```

At this point the existing application has become an nmopt reduced
state/adjoint service.

The reference-taking constructor does not own `model`; the model and every
application object captured by its callbacks must outlive `reduced`.

### Step 8 – solve through nmopt

```cpp
#include "nmopt/solvers/reduced_gradient.hpp"

using Solver =
  nmopt::solvers::ReducedSearchSolverT<Backend>;

MyOptimalControlProblem::Vector initial_values(
  problem.control_dimension());
initial_values = 0.0;

Primal initial_control(
  partition.control_layout(),
  {initial_values});

nmopt::solvers::ReducedSolverParameters parameters;
parameters.maximum_iterations = 5000;
parameters.gradient_tolerance = 1.0e-6;

Solver solver(
  reduced,
  metric,
  parameters);

const auto result =
  solver.solve(initial_control);
```

The iteration/tolerance values above are application choices. The Step-4
examples use frozen values for their own comparison; they are not framework
default recommendations.

For solver families and the meaning of their fields, continue with
[Optimization](optimization.md).

### Step 9 – return to application-owned output

The optimizer's final reduced evaluation already retains the final state:

```cpp
const auto &state =
  result.final_evaluation.state.block(0);

problem.output_results(
  state,
  output_path);
```

This keeps output where it belongs: the deal.II application still owns mesh
and field interpretation.

Do not perform an extra state solve only to obtain output unless that solve is
an intentional independent verification.

## Put the nmopt pieces in one binding object

A useful application layout is:

```cpp
class MyBinding final
{
public:
  using Model =
    nmopt::contract::CallbackExecutableModelT<Backend>;

  using Partition =
    nmopt::contract::StateControlPartitionT<Backend>;

  using Reduced =
    nmopt::contract::ReducedDTOT<Backend>;

  explicit MyBinding(MyOptimalControlProblem &problem)
    : problem_(problem)
    , variable_layout_(make_variable_layout(problem_))
    , test_layout_(make_test_layout(problem_))
    , model_(make_model(
        problem_,
        variable_layout_,
        test_layout_))
    , partition_(model_, 0, 1)
    , metric_(partition_.control_layout())
    , reduced_(
        model_,
        partition_,
        make_solvers(
          problem_,
          partition_.state_layout(),
          test_layout_))
  {}

  const LayoutPtr &control_layout() const
  {
    return partition_.control_layout();
  }

  const auto &metric() const
  {
    return metric_;
  }

  const Reduced &reduced() const
  {
    return reduced_;
  }

private:
  MyOptimalControlProblem &problem_;

  LayoutPtr variable_layout_;
  LayoutPtr test_layout_;

  Model       model_;
  Partition   partition_;
  MyMetric    metric_;
  Reduced     reduced_;
};
```

This is an application pattern, not a required nmopt base class.

The Step-4 bindings delete copy/move operations because the members borrow
from one another and from the native application. That is a simple way to make
the local lifetime graph stable.

## Problem B: when native and optimization coordinates differ

A real deal.II application often does not optimize directly in the full
physical coefficient vector.

Step-4 Problem B uses essential-boundary constraints. Its nmopt state/test
space is the vector of independent/free coefficients, while output uses the
full reconstructed field:

```text
physical/full state
        │
        │ restrict to independent coordinates
        ▼
independent state coordinates
        │
        │ nmopt residual / state / adjoint services
        ▼
independent optimized state
        │
        │ reconstruct / embed into physical coordinates
        ▼
physical/full optimized state
        │
        ▼
application output
```

The application-side facade owns those coordinate maps.

A residual can therefore look conceptually like:

```cpp
Vector residual(
  const Vector &free_state,
  const Vector &full_control) const
{
  const Vector full_state =
    coordinates.embed_free(free_state);

  Vector full_residual =
    apply_native_state_operator(full_state);

  return coordinates.restrict(
    full_residual
    - native_rhs
    - control_coupling(full_control));
}
```

Objective evaluation can reconstruct before applying the physical mass metric:

```cpp
const Vector physical_state =
  coordinates.reconstruct(free_state);

const Vector state_mass =
  mass.mass_apply(physical_state);

const Vector control_mass =
  mass.mass_apply(control);
```

The JVP/VJP must use the same maps and their transpose relationships.

The canonical code is:

- [Problem B application integration](../../apps/external-dealii/step-4/integration/problem_b.hpp);
- [Problem B nmopt binding](../../apps/external-dealii/step-4/minimal/problem_b_binding.hpp);
- [Problem B minimal consumer](../../apps/external-dealii/step-4/minimal/problem_b.cc).

The important rule is not that nmopt requires free coordinates. It is that the
binding must expose one coherent independent-coordinate contract; physical
reconstruction remains application-owned.

## Bring up the adapter before running a long optimization

The shortest debugging path is to test each OCP operation at the boundary
where it is introduced.

### Check the native state solve against the residual

For a representative control $u$, compute the returned state $y(u)$ and check
that

$$
E(y(u),u)
$$

has the residual expected from the native solve tolerance.

If this fails, the problem is still in the native OCP facade; do not debug the
optimizer yet.

### Check JVP against a finite difference

For a point $x$ and tangent $v$:

$$
\frac{E(x+\varepsilon v)-E(x)}{\varepsilon}
\approx
E'(x)v.
$$

This catches missing control/state terms and stale linearization points.

### Check the VJP pairing

For tangent $v$ and test seed $q$:

$$
\langle E'(x)v,q\rangle
\approx
\langle E'(x)^{\ast}q,v\rangle.
$$

In nmopt types:

```cpp
const auto jvp =
  model.residual_jvp(point, tangent);

const auto vjp =
  model.residual_vjp(point, test_seed);

const double left =
  nmopt::contract::pair(jvp, test_seed);

const double right =
  nmopt::contract::pair(vjp, tangent);
```

This is the most direct test that the state **and control** pieces of the VJP
are consistent with the JVP.

### Check the objective derivative before the reduced derivative

Compare $J'(y,u)$ against a finite difference of $J$ while holding the same
coordinate representation.

Only after the full derivative is correct should you test the reduced
quantity

$$
j'(u)[d]
$$

against a reduced finite difference that includes a state solve.

This progression localizes errors. A failed reduced-gradient check otherwise
mixes state solves, objective derivatives, adjoint sign conventions, VJP
implementation, and metric geometry in one number.

The Step-4 `verification/` and `evaluation/` directories contain stronger
repository evidence, but they are not dependencies of the minimal integration
path.

## Evaluate the reduced service directly during bring-up

Before choosing an optimizer, inspect one control through the DTO.

Value only:

```cpp
const auto value =
  reduced.evaluate_value(control);
```

This performs the state solve and objective evaluation.

Add the derivative:

```cpp
const auto evaluation =
  reduced.augment_derivative(value);
```

This evaluates the objective derivative, solves the adjoint, applies the VJP,
and forms the reduced covector.

Or request both:

```cpp
const auto evaluation =
  reduced.evaluate(control);
```

Convert the reduced covector to the metric gradient with:

```cpp
const auto gradient =
  reduced.gradient_direction(
    evaluation.reduced_derivative,
    metric);
```

This is a useful seam for debugging because it lets the application verify one
state/adjoint/gradient evaluation independently of line search or quasi-Newton
history.

## Optional projection constraints

An existing application may also supply a `ConstraintT<Backend>` to projected
reduced search.

That constraint must use a layout compatible with the control metric and must
report projection support for the actual metric realization.

Do not implement coefficient clipping as a generic “box projection” unless it
is really the projection in the declared metric.

The canonical Step-4 A/B minimal consumers are unconstrained; they demonstrate
the direct reduced path, not arbitrary constrained external applications.

See [Optimization](optimization.md) for the currently supported projected
search combinations.

## Lifetime and solve evidence

The binding normally borrows the existing application:

```text
application
    │
    ▼
OCP facade
    │
    ▼
binding
├── layouts
├── model callbacks
├── partition
├── metric
└── reduced DTO
    │
    ▼
solver
```

Keep the application alive while any callback can reach it. Keep the model
alive while the reference-taking `ReducedDTOT` uses it. Keep metric-owned
native matrices/solvers alive while optimization can call them.

There is also an owning DTO constructor taking a
`std::shared_ptr<const ExecutableModelT<Backend>>` plus an optional lifetime
owner. That can retain an owner explicitly, but it cannot discover ownership
of raw pointers captured inside arbitrary callbacks.

`LinearSolveReport` is part of the reduced contract for the same reason:
nmopt cannot inspect an external application's native linear solve after the
fact. The application adapter must say whether the solve converged and retain
its algorithm, iteration, tolerance, residual, and termination evidence.

## Map the repository example to your own code

The post-refactor Step-4 tree separates the responsibilities clearly:

```text
integration/
    application-local OCP facade
    coordinate maps
    native objective/residual/solve operations

minimal/
    nmopt layouts
    CallbackExecutableModelT
    state/control partition
    metric
    ReducedDTOT binding
    small optimizer consumer

evaluation/ + verification/
    stronger comparison / derivative evidence

diagnostics/
    optional instrumentation support
```

For an existing deal.II application, the practical conversion sequence is:

1. **Do not rewrite the PDE solver.** Add a thin OCP facade around its native
   state/operator/output code.
2. Give that facade explicit state/control coordinate dimensions.
3. Implement residual, JVP, VJP, objective, and objective derivative in those
   coordinates.
4. Expose truthful state/adjoint solve evidence.
5. Decide the control metric and implement/apply its Riesz map.
6. Put nmopt layouts/model/partition/metric/DTO in a binding object beside the
   application, not inside the PDE assembly internals.
7. Bring up one reduced evaluation and derivative test.
8. Pass the resulting DTO and metric to the optimizer.
9. Reconstruct/output the final state through the original application.

That is the boundary the Step-4 `minimal/` path is intended to demonstrate.

For the repository evidence and historical comparison around this integration,
see the
[Step-4 integration report](../../apps/external-dealii/step-4/integration-report.md).
