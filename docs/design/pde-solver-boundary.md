# PDE, formulation, and solver boundary

## Status and authority

This document defines the long-lived boundary between semantic problem
description, numerical realization, formulations, optimization algorithms, and
applications. It governs both cleanup of the current deal.II implementation
and integration of an existing deal.II application.

The [interface specification](interface-specification.md) remains authoritative
for the mathematical and semantic model. The
[v0 executable contract](../implementation/v0/executable-contract.md) remains
the exact record of the currently implemented solver-facing API until that
implementation record is retired or superseded. This document decides how
implementations reach those contracts and which responsibilities must not cross
the boundary.

This is a design document, not a mutable implementation-status ledger. Current
refactor status belongs in the corresponding planning roadmap.

## Decision

Keep the existing small executable and formulation contracts. Make the semantic
compiler and an existing deal.II application two independent producers of the
same downstream mathematical operations:

```text
nmopt ProblemSpec -> compiler/lowering --------\
                                              +-> executable operations
existing deal.II application -> adapter -------/   + state/adjoint solves
                                                  + metric
                                                  + optional constraint/Hessian
                                                            |
                                                            v
                                                       formulation
                                                            |
                                                            v
                                                       optimizer
```

The two paths converge on solver-facing mathematical operations. They do not
need to share problem-description types, finite-element assembly code, output
code, compiler provenance, or application configuration.

The formulation owns **when** mathematical operations are required. The
application or compiled numerical realization owns **how** those operations are
implemented.

## Stable solver-facing surface

The following existing contracts are the common downstream boundary:

- `ExecutableModelT<Backend>` for variable and test layouts, residual,
  residual JVP, residual VJP, objective value, and objective derivative;
- `StateAdjointSolversT<Backend>` for state and adjoint solution;
- `MetricT<Backend>` for decision-space primal-dual identification;
- optional `ConstraintT<Backend>` for feasibility operations; and
- optional `ReducedHessianT<Backend>` for methods that explicitly request a
  reduced Hessian action.

`ReducedDTOT<Backend>` remains the owner of first-order reduced
state-adjoint-gradient orchestration. Optimization algorithms consume
`ReducedDTOT`, metrics, constraints, and optional Hessian actions without
knowing whether the producer was a semantic compiler or an external
application.

No new public hierarchy of PDE families, residual operators, objectives,
observations, state solvers, linear solvers, or output providers is required to
establish this boundary. Prefer existing contracts, callbacks, values, free
functions, and small concrete numerical components.

`CompiledProblemT` is packaging for the semantic/compiler path. It is not a
prerequisite for using a solver and must not become the universal integration
type for existing applications.

## The missing middle layer

The current implementation often jumps from a validated semantic graph to one
complete concrete `*Model`, then erases that concrete type to
`ExecutableModelT`. Native application code later needs dimensions,
reconstruction, diagnostics, or field output and therefore recovers the
concrete model with `dynamic_cast`.

That sequence erases too much information too early:

```text
semantic graph
    -> complete concrete model
    -> ExecutableModelT
    -> application dynamic_cast back to concrete model
```

The target architecture retains the typed numerical realization while exposing
only the narrow mathematical operations required by formulations and
optimizers:

```text
semantic description
        |
        v
closed lowering decision
        |
        v
typed numerical realization
        |\
        | \
        |  +--> typed/native application view
        |        reconstruction, dimensions, diagnostics, output
        |
        +-----> solver-facing ports
                 E, E'v, E'^*p, J, J'
                 state/adjoint solves
                 metric/optional capabilities
```

Type erasure is therefore a useful solver boundary, not the ownership boundary
for all compiled numerical state.

## Semantic meaning versus numerical realization

The semantic layer is backend-neutral, but it is not limited to purely
continuous mathematics. It may contain:

1. mathematical declarations, such as variables, regions, residual terms,
   observations, losses, and metrics;
2. analytical policy commitments, such as trace, nullspace, target-regularity,
   or conormal requirements; and
3. backend-neutral discrete-realization choices, such as a selected trace,
   quadrature, control representation, or negative-norm realization.

The semantic layer must not contain deal.II objects, concrete sparse matrices,
linear-solver implementations, or application filesystem behavior.

The compiler/lowerer turns those declarations and policies into typed numerical
objects such as finite elements, DoF handlers, affine constraints, sparse
operators, coordinate maps, observation operators, and solve callbacks.

## Typed numerical responsibilities

The target middle layer should preserve independently varying numerical
responsibilities. These are responsibility categories, not a requirement to
create one class for every row.

| Concern | Owns or realizes | Must not own merely by convenience |
| --- | --- | --- |
| State space | FE, DoF handler, mesh association, physical state dimension | Objective, optimizer policy, output orchestration |
| State coordinates | Independent/full coordinates, reconstruction, tangent embedding, dual pullback, fixed lifting | PDE-family selection, line search |
| Physical state operator | Discrete residual/operator actions and operator properties | Objective regularization, search metric |
| Gauge/nullspace policy | Mean-zero or related admissibility needed for state/adjoint solvability | Generic optimizer behavior |
| Decision discretization | Decision coordinates, layout, reusable mass/stiffness/trace operators | One fixed regularization or search metric |
| Decision-to-PDE influence | Affine coupling, boundary lifting, or parameter derivative machinery | A universal assumption that every decision enters through one matrix |
| Observation realization | Physical-field measurements and optional value/JVP/VJP data | Loss semantics or optimizer geometry |
| Tracking/loss realization | Scalar objective contribution and derivative data, often a compiled quadratic form | PDE solve order |
| Regularization realization | Objective penalty and its derivative/curvature | Search metric merely because the same matrix is reused |
| Metric realization | `apply`, `inverse_apply`, metric-specific solve policy and witness | Objective semantics |
| Constraint realization | Feasibility/projection capability supported by the selected decision realization | Residual physics |
| Solve realization | Concrete state/adjoint inversion and convergence report | Optimization iteration/globalization policy |
| Reduced Hessian capability | Optional second-order reduced action | Universal PDE contract |
| Native application view | Reconstruction, dimensions, diagnostics, field output | Generic formulation or optimizer logic |

The same numerical operator may serve several roles. Sharing a matrix does not
merge the semantic ownership of those roles.

## State coordinates and reconstruction

A physical state and the state coordinates exposed to a formulation need not be
the same vector.

For fixed essential data, a useful discrete form is

```math
  y_{\mathrm{phys}} = P_{h}\widehat y + \ell_{h}.
```

The corresponding numerical operations are:

```text
reconstruct:       y_hat -> P y_hat + ell
embed tangent:    dy_hat -> P dy_hat
pull back dual:        q -> P^T q
```

Dirichlet control extends the same pattern:

```math
  y_{\mathrm{phys}} = P_{h}\widehat y + L_{D,h}u + \ell_{0,h}.
```

Then both state and control derivatives receive pullbacks through $P^T$ and
$L_{D}^T$. Observations and native output should act on the physical field;
state-coordinate machinery is responsible for translating to and from
solver-facing coordinates.

Full-coordinate representations with identity rows and independent-coordinate
representations with reconstruction are implementation policies, not different
PDE families.

## Decision realizations

Current applications demonstrate at least three mathematically different ways
a decision variable can affect the PDE.

### Affine forcing controls

Distributed and Neumann controls have the fixed-linear form

```math
  E(y,u) = A y - f - B u,
```

with decision actions $B u$ and $B^T p$. Different control discretizations may
share this internal pattern without creating a universal public control base
class.

### State-reconstruction controls

Dirichlet control enters through the physical-state reconstruction
$L_D u$, so residuals and observations see the reconstructed physical field.
Its natural numerical operations are lifting and dual pullback rather than a
standalone forcing matrix.

### Operator parameters

Coefficient identification has a point-dependent operator, for example

```math
  E(y,m) = A(m)y - f.
```

Its derivative contains both a state action and a parameter derivative of the
operator. It cannot be reduced to a fixed decision-coupling matrix.

The universal solver-facing contract therefore remains $E$, $E'v$, and
$E'^*p$. Any lower-level fixed-linear coupling abstraction is only an internal
specialization justified by demonstrated code reuse.

## Observations and losses

Observations and losses are semantically distinct:

```text
physical state -> observation map O -> loss
```

For quadratic tracking, numerical lowering may compile the composition into

```math
  J_{\mathrm{track}}(y)
  = \frac{1}{2} y^{\mathsf T} Q y - q^{\mathsf T} y + \frac{1}{2} c.
```

That representation is often the efficient optimizer-facing form. Explicit
observation value/JVP/VJP data should be retained when another consumer needs
the measurement itself, as with point sensors, boundary traces, normal flux,
verification, or diagnostics.

Do not force every observation into an explicit sampled-vector interface merely
for uniformity. A compiled quadratic form and an explicit observation map are
two useful numerical representations with different consumers.

## Regularization, metric, and reduced Hessian

Objective regularization, optimization metric, and reduced Hessian are distinct
mathematical roles.

For

```math
  J_{\mathrm{reg}}(u)
  = \frac{\alpha}{2} u^{\mathsf T} R u,
```

$R$ belongs to the objective definition.

A search metric $G$ maps primal and dual decision representations and may be
chosen independently of $R$. Changing $G$ changes optimization geometry but not
the objective.

A reduced Hessian combines PDE sensitivity, observation curvature, and direct
objective curvature. In a fixed linear-quadratic case it has the form

```math
  H_{\mathrm{red}}
  = B^{\mathsf T}A^{-\mathsf T}Q A^{-1}B + \alpha R.
```

It is therefore an optional composed formulation capability, not a primitive
property of every PDE discretization.

Mass, stiffness, trace, or quotient operators should be owned as reusable
numerical ingredients where several consumers need them. A type named for one
role should not silently become the authority for another role merely because
it exposes the same operator.

## State and adjoint solves

`StateAdjointSolversT` is already the inversion-of-control seam required by the
formulation. State/adjoint implementations may use CG, GMRES, direct solves,
multigrid, PETSc, Trilinos, or application-specific code without changing the
formulation or optimizer.

The numerical realization may expose facts such as symmetry, nullspace, or
point dependence. The application/compiler composition decides how those facts
map to a concrete solve policy.

Do not introduce a universal public linear-solver interface unless multiple
real applications demonstrate shared code that the new interface would remove.
Metric-internal solves are different: if applying an $H^{-1}$ or trace
$H^{1/2}$ operator mathematically requires an inversion, that solve is intrinsic
to the metric operation.

## Compiler path

The native compiler path remains:

```text
ProblemSpec
  -> semantic validation and resolution
  -> closed compilation/lowering decision
  -> typed numerical realization
  -> solver-facing operations + typed native view
  -> formulation and optimizer
```

The compiler must not represent the same realization decision repeatedly as
parallel target enums, predicate sets, boolean fields, construction switches,
and manifest compatibility copies.

A future closed plan should describe the selected numerical axes once. Exact
C++ representation is deferred to the implementation unit that can prove which
fields and variants remove existing duplication.

Component-aware lowering does not require arbitrary composition of every
semantic graph. A compiler may reject unsupported combinations explicitly. The
architectural requirement is that supported combinations are not needlessly
re-encoded as complete problem identities after their independently meaningful
choices have already been resolved.

## Existing deal.II application path

An existing application must be usable without adopting `ProblemSpec`, the
semantic resolver, `DealiiCompiler`, `CompiledProblemT`, recipes, scenarios, or
benchmark manifests.

The application's original PDE implementation remains independently usable. A
small adapter supplies:

1. variable and test layouts;
2. residual, JVP, and VJP callbacks;
3. objective value and derivative callbacks;
4. state and adjoint solve callbacks; and
5. a metric, plus a constraint or Hessian only when required.

The adapter may construct a callback-backed `ExecutableModelT`, but the
original PDE class should not be required to inherit an nmopt interface.
Callbacks may borrow or own captured objects; lifetime choice must remain
explicit. Do not copy meshes, DoF handlers, or matrices merely to conceal a
lifetime requirement.

## Application and output boundary

Typed numerical code owns the backend objects and operations needed to expose
physical fields, reconstruction, coordinates, dimensions, and backend-specific
diagnostic quantities.

The application owns the decision to produce output: filesystem paths,
filenames, serialization, visualization, artifact organization, and when
output is written.

`ExecutableModelT`, `ReducedDTOT`, optimizers, and reusable numerical
components must not own filesystem or application-output policy.

The compiler/native application path therefore retains a typed native view
beside the erased solver-facing view. Application code must not recover private
compiler target types through `dynamic_cast` from `ExecutableModelT`.

A native view may expose typed field data, reconstruction operations, geometry,
or side-effect-free diagnostic access needed by the application. This does not
imply a universal field-output interface.

Introduce a shared output abstraction only if multiple real applications
demonstrate one stable requirement and the abstraction removes more code or
duplication than it adds.

## Architectural correspondence

The rules in this document are implementation-independent.

The R0 refactor review contains a descriptive
[architecture map](../planning/review/pde-solver-refactor/architecture-map.md)
showing how the code at the start of the refactor corresponds to these
boundaries. It records the native and external producer paths, compiler and
lowering flow, typed-numerics axes, formulation runtime, ownership and
lifetime, and the current coupling pressure points.

That map is review evidence rather than design authority. This document remains
authoritative for the boundary itself.

The architecture map is intentionally allowed to name current concrete types
and implementation pressure that may disappear during the refactor. Durable
parts of that correspondence can later be promoted into the permanent design
documentation after the implementation stabilizes.

## Design non-goals

This boundary does not require:

- a maximally general PDE object model or public PDE-family hierarchy;
- one public class corresponding to every responsibility named in the
  architecture diagrams;
- a universal residual, observation, output, or linear-solver hierarchy;
- a symbolic weak-form or general PDE expression DSL;
- arbitrary lowering of every semantically valid component combination; or
- replacing demonstrated concrete/value/callback composition with inheritance
  merely to make the architecture visually uniform.

Responsibility names such as “typed numerical realization,” “state
coordinates,” and “decision realization” describe architectural roles. They
do not prescribe matching C++ type names or inheritance structures.
