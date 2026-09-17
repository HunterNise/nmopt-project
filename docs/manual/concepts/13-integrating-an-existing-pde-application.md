# Integrating an existing PDE application

## The second execution path starts with numerical objects, not `ProblemSpec`

The compiled path begins with a semantic problem and asks `nmopt` to construct a
registered numerical realization.

An existing PDE application starts at the other end.

It may already own:

```text
mesh
finite-element spaces
DoFHandler
assembled operators
boundary treatment
linear solvers
output
```

Recreating those objects through the semantic/compiler path may add work without
adding value.

The direct integration path instead asks:

> Which mathematical operations does the application already provide, and how can
> they be exposed through the contracts consumed by an `nmopt` formulation?

The architecture is:

```text
semantic/compiler path         native application path
    ProblemSpec                    existing PDE application
        ↓                                 ↓
    compiler                       thin binding
            ↘                      ↙
               numerical contracts
                       ↓
               formulations / solvers
```

The paths converge at the numerical/formulation boundary.

`CompiledProblemT` is therefore one producer of those services, not the universal
application interface.

## 1. Preserve the native application's responsibilities

A direct integration should begin by identifying which responsibilities should remain
inside the existing application.

For the Step-4 case study, the adapted application continues to own:

```text
triangulation
finite-element space
DoFHandler
stiffness matrix
load vector
Dirichlet treatment
native CG solve
native output
```

`nmopt` does not need to become the owner of those objects.

The application only needs reusable seams around the operations that optimization
requires.

This is the central difference from the compiler path:

```text
compiled path
    nmopt constructs numerical services

native path
    application supplies numerical services
```

## 2. A monolithic `run()` method is often the first integration obstacle

A traditional PDE example may expose only:

```text
run()
    setup
    assemble
    solve
    output
```

That is enough for one forward simulation.

Optimization needs repeated calls to individual numerical operations.

A reusable application seam therefore often needs to separate:

```text
prepare / assemble once

read numerical dimensions and operator views

solve for a supplied right-hand side

write a supplied state
```

This does not mean redesigning the whole PDE application around `nmopt`.

It means making the existing numerical work callable at the points where an outer
algorithm needs it.

## 3. First define the optimal-control mathematics around the native PDE

A forward PDE application alone is not yet an optimization problem.

For Step-4 Problem B, the added optimal-control model is

```math
K\mathbf z
=
\mathbf b_{F}
+
B\mathbf u,
```

with

```math
\mathbf y_{\mathrm{phys}}
=
P\mathbf z+\boldsymbol\ell,
```

and objective

```math
J(\mathbf z,\mathbf u)
=
\frac12
(P\mathbf z+\boldsymbol\ell)^{\mathsf T}
M
(P\mathbf z+\boldsymbol\ell)
+
\frac12
\mathbf u^{\mathsf T}M\mathbf u.
```

The native Step-4 application supplies the prepared stiffness matrix and solver.

The optimal-control layer adds:

```text
independent state coordinates

physical-state reconstruction

mass matrix

control coupling

objective

objective derivative

residual JVP / VJP

control metric
```

That division is useful for any external integration.

The PDE application need not already have an "optimization interface"; the adapter
can assemble one from its native numerical pieces.

## 4. Decide the state and control coordinates before writing callbacks

The Step-4 physical state has 289 finite-element coefficients in the studied
configuration.

Fixed boundary values leave 225 independent state coefficients.

Problem B therefore uses:

```text
state block
    225 independent coefficients

control block
    289 physical FE coefficients

test block
    225 independent equations
```

The relationship is

```math
\mathbf y_{\mathrm{phys}}
=
P\mathbf z+\boldsymbol\ell.
```

The control-to-state coupling is

$$
B=P^{\mathsf T}M.
$$

This coordinate decision must be settled before the public layouts are created.

The contract should reflect the application's mathematical coordinates rather than
forcing everything into the physical state vector merely because that is what the
native matrix stores.

## 5. `BlockLayout` records the spaces that native vectors represent

The Step-4 binding constructs a variable layout with two blocks:

```text
state
control
```

and dimensions

```text
state_dimension()
control_dimension()
```

It separately constructs a one-block test layout:

```text
state_test
```

This is important even though all blocks use

```text
dealii::Vector<double>
```

under the serial backend.

A vector's native C++ type does not tell us whether it represents

```text
state
control
test seed
state derivative
control derivative
```

The layouts carry that identity across the public contract boundary.

## 6. `CallbackExecutableModelT` is the main model adapter

For a direct application, the most important generic adapter is
`CallbackExecutableModelT`.

It packages five operations:

```text
residual

residual JVP

residual VJP

objective

objective derivative
```

For Problem B, the binding unwraps the `PrimalBlockT` blocks, calls the native
`ProblemB` methods, and wraps the returned vectors in the appropriate covector
layout.

The pattern is:

```text
nmopt contract value
        ↓ unwrap native vectors
ProblemB method
        ↓
native vector result
        ↓ wrap with layout/role
nmopt contract value
```

The callback is therefore a representation adapter, not a new implementation of the
PDE mathematics.

## 7. The residual callback should call the application's residual meaning directly

Problem B defines

```math
E(\mathbf z,\mathbf u)
=
K\mathbf z
-
\mathbf b_{F}
-
B\mathbf u.
```

Its native method is conceptually:

```cpp
Vector residual(const Vector &state,
                const Vector &control) const;
```

The binding callback is correspondingly small:

```cpp
[problem_ptr, test_layout](const Primal &variables) {
  auto value =
    problem_ptr->residual(variables.block(0),
                          variables.block(1));

  return Covector(test_layout, {std::move(value)});
}
```

The adapter should not reproduce matrix algebra that `ProblemB` already owns.

That keeps one source of truth for the OCP residual.

## 8. JVP and VJP callbacks preserve the native derivative implementation

The same rule applies to derivatives.

Problem B supplies

```text
residual_jvp(state_tangent, control_tangent)

residual_vjp(test_seed)
```

The VJP returns two native covector components:

```text
state
control
```

which the binding wraps in the full variable layout.

This is preferable to asking `nmopt` to recover a transpose numerically from the JVP.

The application can use whatever internal representation makes the transpose action
correct and efficient.

## 9. The objective remains application mathematics

Problem B evaluates

```math
J(\mathbf z,\mathbf u)
=
\frac12
\mathbf y_{\mathrm{phys}}^{\mathsf T}
M
\mathbf y_{\mathrm{phys}}
+
\frac12
\mathbf u^{\mathsf T}
M
\mathbf u.
```

Its native objective derivative computes

```text
state covector
    P^T M y_phys

control covector
    M u
```

The callback binding simply exposes those values as a `CovectorBlockT`.

This separation is important:

```text
application
    defines J and J'

ReducedDTOT
    combines J' with state dependence through an adjoint
```

The formulation does not own the application's objective.

## 10. Residual evaluation and state solving remain separate

An application may already have a good native linear or nonlinear solve.

Do not replace it merely to satisfy the executable-model callback.

The direct path uses two distinct surfaces:

```text
CallbackExecutableModelT
    residual / JVP / VJP / objective operations

StateAdjointSolversT
    solve_state / solve_adjoint
```

For Problem B, the state solve:

```text
build controlled full RHS

call Step-4 native solve

restrict full solution to independent state coordinates

return solve evidence
```

The executable residual and the state solver therefore share native operators without
being the same operation.

## 11. The adjoint solver should expose the transpose solve the formulation needs

Problem B's adjoint solve accepts the state-objective derivative in independent
coordinates.

It embeds that right-hand side into the physical system and calls the native solve.

For the symmetric Step-4 stiffness matrix, the same native solve algorithm realizes
both

```math
K\mathbf z=\mathbf b
```

and

```math
K^{\mathsf T}\mathbf p=\mathbf r.
```

A nonsymmetric application could expose a distinct transpose solve.

The public contract asks for the mathematical service, not for one particular
implementation strategy.

## 12. Native solve evidence should cross the boundary too

The Step-4 native solver returns evidence such as

```text
converged

iteration count

initial residual

final residual
```

The binding converts that information into `LinearSolveReport`.

This is better than returning only the solution vector and pretending every native
solve succeeded.

The formulation can then retain and check the solve report as part of its evaluation.

Direct integration should preserve useful native evidence rather than discard it at
the adapter boundary.

## 13. The metric can wrap an existing native Riesz implementation

Problem B uses the FE mass matrix as the control metric.

The native `ProblemBMetric` already implements

```text
apply

inverse_apply
```

on native control vectors.

The `ProblemBNmoptMetric` adapter implements `MetricT` by:

```text
check layout

call native metric

wrap primal/covector result
```

The inverse adapter also translates the native metric-solve evidence and rejects a
failed solve.

Again, the binding adds the public contract but does not reimplement the mass solve.

## 14. The direct path can construct the reduced formulation explicitly

Once the binding has

```text
CallbackExecutableModelT

StateControlPartitionT

StateAdjointSolversT

MetricT
```

it can construct

```text
ReducedDTOT
```

directly.

Problem B does this in its binding constructor:

```text
model
    ↓
partition
    ↓
state/adjoint solvers
    ↓
metric
    ↓
reduced DTO
```

No `ProblemSpec`, semantic validator, or compiler is needed.

The problem is already numerical and the application has supplied the required
operations.

## 15. The optimizer sees the same reduced interface as the compiled path

After construction, the optimizer receives:

```text
ReducedDTOT

MetricT

optional ConstraintT

initial control
```

This is the same kind of solver-facing interface used by a compiler-produced reduced
problem.

The optimizer therefore does not need a special "external application algorithm".

The producer path has disappeared from view.

This is the main architectural payoff of the contract boundary.

## 16. `ProblemBNmoptBinding` keeps the pieces together

The concrete Step-4 binding owns:

```text
prepared application

ProblemB

native metric

layouts

callback model

state/control partition

state/adjoint solvers

metric adapter

ReducedDTOT
```

in one object.

This arrangement is convenient because several later objects borrow earlier ones.

It also makes their lifetime order visible in the member list.

The class deliberately deletes copy and move operations, avoiding accidental movement
of an object graph containing internal references.

This is one valid integration pattern.

The framework does not require every external application to package its binding in
exactly the same class shape.

## 17. Type erasure does not transfer ownership

The callback model stores callable objects.

Those callbacks may capture pointers or references to native application objects.

For Problem B, callbacks capture a pointer to the `ProblemB` instance owned by the
binding.

Therefore:

> Type erasure is the solver boundary, not the ownership boundary.

Wrapping a method in `std::function` does not make the callback own the matrix,
`DoFHandler`, or application instance that method uses.

The owner must outlive every object that can invoke the callback.

## 18. A safe lifetime chain is easier to reason about than scattered callbacks

For the concrete Problem B binding, the ownership/borrowing chain is:

```text
ProblemBNmoptBinding
    │
    ├── owns prepared Step-4 application
    │
    ├── owns ProblemB
    │       └── borrows application
    │
    ├── owns callback model
    │       └── callbacks borrow ProblemB
    │
    ├── owns solve callbacks
    │       └── borrow ProblemB
    │
    ├── owns metric adapter
    │       └── borrows native metric
    │
    └── owns ReducedDTOT
            └── borrows model / solvers
```

The important rule is that destruction occurs only after all downstream borrowers are
finished.

An application can use shared ownership or an outer owner instead, but the lifetime
must remain equally explicit.

## 19. Do not hide lifetime problems behind global state

A tempting integration shortcut is to store the current application in a global or
static pointer so callbacks can find it.

That makes ownership less visible and complicates:

```text
multiple simultaneous problems

tests

reentrant use

threaded execution

clean teardown
```

Capturing a well-scoped owner or borrowed reference in an explicit binding is easier
to reason about.

The Step-4 binding demonstrates that a small local object graph is sufficient.

## 20. Coordinate reconstruction should remain with the application mathematics

Problem B owns `ProblemBCoordinates`.

It knows:

```text
which state DoFs are independent

how to embed a free state

how to reconstruct the physical state

how to restrict a physical covector
```

Those choices depend on the native Step-4 boundary treatment.

The generic optimizer should not learn them.

The binding exposes only the resulting state layout and model operations.

This keeps the application-specific coordinate convention on the application side of
the boundary.

## 21. The control coupling should remain native too

Problem B owns the mass/coupling implementation that realizes

$$
B=P^{\mathsf T}M.
$$

The model uses it for:

```text
residual

residual JVP

control VJP

state RHS construction
```

All four operations must agree on the same coupling.

Keeping them in the application OCP layer reduces the chance that the binding
duplicates one direction with a different sign, coordinate convention, or assembly.

The adapter should expose a coherent native model, not assemble one from unrelated
callbacks.

## 22. Output should remain where the native finite-element context lives

A generic optimizer does not know how to write a finite-element solution.

Native output may require:

```text
physical state reconstruction

DoFHandler

triangulation

cell/point data associations

application naming conventions
```

The external application already owns that context.

After optimization, use the final retained state/control and pass them back to the
native problem/application output path.

There is no need to make VTK or deal.II output part of `ExecutableModelT`.

## 23. Retained evaluations avoid unnecessary final solves

A reduced solver already evaluates the accepted final control.

That evaluation retains the state and, after derivative augmentation, the adjoint.

A native application should normally reuse those values for output.

Conceptually:

```text
final solver report
    final control
    retained state
    retained adjoint
        ↓
native reconstruction / writer
```

Re-solving the PDE after optimization merely to obtain a state for output can waste
work and may even use a slightly different solve policy.

## 24. The application may retain its native solver policy

Direct integration does not require replacing

```text
CG

GMRES

sparse direct solve

nonlinear iteration

custom preconditioner
```

with an `nmopt` solver.

`nmopt` needs a callable service with the declared mathematical meaning and sufficient
evidence.

The Step-4 case retains its CG solver with identity preconditioning.

Another application can retain a completely different native solve strategy.

This is one of the main reasons to support a direct integration path.

## 25. The backend supplies storage algebra, not PDE ownership

Problem B and the generic solver both use

```text
dealii::Vector<double>
```

through `SerialBackend`.

The backend supplies operations such as

```text
zeros
dot
add_scaled
scale
```

The application still owns

```text
stiffness action

mass action

coordinate maps

PDE solves

metric solve
```

Using the same native vector type makes the adapter inexpensive, but it does not merge
the backend and PDE realization layers.

## 26. A native application can expose only the capabilities its solver needs

For steepest descent with Armijo, a reduced problem needs:

```text
reduced value / derivative

metric inverse

state/adjoint solves underneath the formulation
```

A Newton method may additionally require a reduced Hessian action.

A constrained algorithm needs a compatible constraint service.

Therefore an integration does not need to implement every optional `nmopt`
capability on day one.

Start from the formulation/algorithm actually required.

Add optional services only when a consumer needs them.

## 27. The five executable-model callbacks are still a deliberate minimum

Even if one algorithm does not call every model action on every iteration,
`CallbackExecutableModelT` describes a coherent executable mathematical model:

```text
residual

JVP

VJP

objective

objective derivative
```

Providing all five makes the model useful for:

```text
reduced formulations

derivative verification

other formulation products

diagnostics
```

The direct path should avoid defining a "solver-specific model" whose mathematics
changes with the selected outer algorithm.

## 28. A thin binding is preferable to a framework-specific PDE base class

The Step-4 case does not derive the native application from an `nmopt` PDE base class.

Instead it adapts ordinary application methods to public contracts.

That keeps the dependency direction favorable:

```text
existing PDE code
    remains native

integration layer
    depends on nmopt contracts

optimizer/formulation
    depends only on contracts
```

The application can continue to be built, tested, and understood as a PDE code
outside the optimization framework.

## 29. The direct path is not a second semantic compiler

An external binding should not reconstruct a pseudo-`ProblemSpec` just to label its
callbacks.

The native path already has its mathematical realization.

Its responsibility is to expose correct operations and layouts.

If semantic discovery/provenance is important for a particular application, metadata
can be recorded separately, but the direct integration does not need to pass through
the semantic validator/compiler before reaching the formulation layer.

The two producer paths are peers.

## 30. Decide between compiled and native integration by ownership reality

The compiled path is natural when:

```text
the problem fits the registered semantic/compiler vocabulary

nmopt should construct the FE realization

the application benefits from recipe/scenario reuse
```

The direct path is natural when:

```text
a mature PDE application already owns the realization

replacing its mesh/assembly/solver would be undesirable

the required numerical operations can be exposed cleanly
```

Neither path is categorically more "native" to `nmopt`.

They solve different integration problems.

## 31. A practical native-integration sequence

For a new external application, the sequence is:

```text
1. Identify what the application already owns.

2. Make repeated preparation / solve / output operations callable.

3. Define the optimization variables and coordinate spaces.

4. Implement the OCP mathematics:
       residual
       JVP
       VJP
       objective
       objective derivative

5. Define public BlockLayouts.

6. Wrap the model in CallbackExecutableModelT.

7. Wrap state and adjoint solves in StateAdjointSolversT.

8. Adapt the control metric to MetricT.

9. Construct the required formulation, usually ReducedDTOT.

10. Pass the formulation + metric + optional constraint to a solver.

11. Return retained results to the native application for output.
```

Only steps 5–10 are `nmopt` contract integration.

The mesh and PDE implementation remain application concerns.

## 32. A compact binding skeleton

A direct reduced binding has the conceptual shape:

```cpp
using Backend = nmopt::dealii_backend::SerialBackend;
using Model = nmopt::contract::CallbackExecutableModelT<Backend>;
using Reduced = nmopt::contract::ReducedDTOT<Backend>;

auto variable_layout = /* state + control blocks */;
auto test_layout = /* state test block */;

Model model(
  variable_layout,
  test_layout,
  residual_callback,
  residual_jvp_callback,
  residual_vjp_callback,
  objective_callback,
  objective_derivative_callback);

nmopt::contract::StateControlPartitionT<Backend>
  partition(model, 0, 1);

nmopt::contract::StateAdjointSolversT<Backend> solvers;
solvers.solve_state = /* native state solve adapter */;
solvers.solve_adjoint = /* native adjoint solve adapter */;

Reduced reduced(model, partition, solvers);

// construct MetricT adapter, then pass reduced + metric to a solver
```

The application-specific work is hidden behind the callback bodies and metric/solve
adapters.

## 33. The Step-4 Problem B binding is deliberately close to that skeleton

`ProblemBNmoptBinding` constructs, in order:

```text
application_

problem_

native_metric_

variable_layout_

test_layout_

model_

partition_

solvers_

metric_

reduced_
```

The order reflects dependency.

The model callbacks borrow `problem_`.

The metric adapter borrows `native_metric_`.

The reduced formulation refers to the model/partition/solvers.

A reader integrating another application can use this as an ownership pattern without
copying the Step-4 mathematics.

## 34. Problem B keeps OCP mathematics out of the adapter

`ProblemB` itself owns methods for:

```text
residual

residual_jvp

residual_vjp

objective

objective_derivative

solve_state

solve_adjoint
```

It also owns its coordinate and mass/coupling helpers.

`ProblemBNmoptBinding` mostly translates those operations into generic types.

This separation is worth emulating.

If the callback file contains all PDE/OCP algebra directly, it becomes difficult to
test the native model independently of the framework binding.

## 35. Native and compiled paths should agree at the formulation boundary

For a reduced problem, both paths should eventually supply equivalent categories of
service:

| Service | Compiled path | Native path |
| --- | --- | --- |
| executable model | compiler-created model | `CallbackExecutableModelT` |
| state/control partition | compiled model layout | explicit `StateControlPartitionT` |
| state/adjoint solves | compiler-created services | native solve adapters |
| metric | compiler-created `MetricT` | metric adapter |
| constraint | compiler-created optional service | native/adapted optional service |
| reduced formulation | `make_reduced_dto()` | explicit `ReducedDTOT` |

The producer details differ.

The formulation mathematics does not.

## 36. The optimizer should not branch on the producer path

Avoid code such as:

```text
if compiled_problem:
    run optimizer A
else if external_application:
    run optimizer B
```

when the two problems satisfy the same formulation and metric contracts.

Prefer:

```text
obtain ReducedDTOT
obtain MetricT
obtain optional ConstraintT
        ↓
run one generic solver
```

The whole purpose of the common contract vocabulary is to remove producer-specific
branches above that boundary.

## 37. External integration does not eliminate verification

A direct binding bypasses semantic/compiler validation.

That makes local numerical checks even more important.

The earlier chapters already developed the relevant identities:

```text
residual Taylor/JVP consistency

JVP/VJP transpose pairing

objective derivative Taylor check

metric apply/inverse consistency

state/adjoint solve reports

reduced derivative check
```

The Step-4 integration keeps dedicated verification/evaluation material for these
contracts.

The point here is not to repeat those derivations, but to remember that the native
path assumes responsibility for satisfying them because no compiler built the model
on its behalf.

## 38. Following the native path through the repository

The case-study narrative is:

- [`apps/external-dealii/step-4/external-integration-overview.md`](../../../apps/external-dealii/step-4/external-integration-overview.md)

The native optimal-control problem is:

- [`apps/external-dealii/step-4/integration/problem_b.hpp`](../../../apps/external-dealii/step-4/integration/problem_b.hpp)

The public binding is:

- [`apps/external-dealii/step-4/integration/nmopt_problem_b_binding.hpp`](../../../apps/external-dealii/step-4/integration/nmopt_problem_b_binding.hpp)

Supporting coordinate, mass/coupling, and metric implementations live beside those
files under:

- [`apps/external-dealii/step-4/integration/`](../../../apps/external-dealii/step-4/integration/)

For exact generic contract signatures, use:

- [External deal.II solver integration](../../reference/external-dealii-solver-integration.md)

The common reduced formulation and solver concepts were developed earlier in:

- [Reduced state–adjoint formulation](05-reduced-state-adjoint-formulation.md)
- [Reduced optimization methods](06-reduced-optimization-methods.md)

## 39. The two execution paths close the concept manual

The full architecture can now be read from either side:

```text
compiled path                          native path
    mathematical problem                  existing numerical application
        ↓                                          ↓
    recipe / ProblemSpec                  OCP operations
        ↓                                          ↓
    semantic validation / compiler        callback / solve / metric adapters
                     ↘                        ↙
                         numerical contracts
                                 ↓
                            formulation
                                 ↓
                               solver
```

The current Step-4 case demonstrates this convergence for the reduced formulation; it does not establish a native external path for every compiler product.

Below the convergence point, the optimizer should not care who constructed the
numerical realization.

Above it, ownership and authoring are intentionally different.

That separation is the final organizing principle of the manual: `nmopt` provides
reusable mathematical contracts, formulations, and algorithms without requiring
every PDE application to surrender its native numerical structure.
