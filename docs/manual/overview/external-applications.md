# Integrating an existing PDE application

`nmopt` does not require every PDE problem to be reconstructed through the semantic
compiler.

If a numerical application already owns a trustworthy finite-element realization –
its mesh, spaces, assembled operators, state solves, boundary handling, and output –
the direct integration path keeps that ownership intact and adapts only the
mathematical operations needed by the formulation and optimizer.

This path is architecturally first-class. It converges with the semantic/compiler
path at the common numerical contracts rather than passing through `ProblemSpec`.

## 1. When this path is the right choice

Direct integration is usually preferable when:

- the PDE code already exists and is independently useful;
- its assembly and solve policies should remain authoritative;
- the application owns numerical objects with nontrivial lifetimes;
- rewriting the problem through `nmopt` semantics would duplicate mature code;
- only a narrow optimization-facing interface is needed.

The semantic/compiler path is usually preferable when the goal is instead to build a
reusable family of supported problems from explicit semantic components and let the
framework construct the deal.II realization.

Neither path is more fundamental. They solve different integration problems.

## 2. What stays in the application

An existing application normally keeps ownership of:

```text
mesh and topology
finite-element spaces and DoF handlers
assembled matrices/operators
boundary-condition realization
state and adjoint numerical solves
application-specific reconstruction
native field output
```

`nmopt` should not become a second owner of those objects merely because optimization
is added.

The direct binding borrows or references the numerical application as needed and
presents its operations through the common mathematical interfaces.

## 3. What `nmopt` needs from the application

For the common first-order formulation boundary, the relevant operations are
conceptually:

```math
\begin{aligned}
E(x) &\quad &&\text{residual evaluation},\\
E'(x) \delta x &&&\text{Jacobian action},\\
E'(x)^{\ast}p &&&\text{transpose / pullback action},\\
J(x) &&&\text{objective value},\\
J'(x) &&&\text{objective derivative}.
\end{aligned}
```

A reduced formulation additionally needs callable state and adjoint solves.

The optimizer also needs a metric on the decision space so that the reduced covector
can be converted to a primal gradient or search direction.

Depending on the selected algorithm or formulation, extra services may be required,
for example a reduced Hessian action or a constraint projection.

The important point is that these are **operations**, not ownership transfers.

## 4. The direct path in one picture

```text
existing numerical application
        │
        ├─ mesh / FE / assembly
        ├─ state and adjoint solves
        ├─ objective and derivative operations
        ├─ metric / constraints
        └─ native reconstruction and output
        │
        ▼
thin numerical binding
        │
        ▼
common nmopt contracts
        │
   ┌────┼───────────────┐
   ▼    ▼               ▼
reduced  KKT / OTD    PDAS
path     products      products
   │
   ▼
selected algorithm / solve service
```

The current external case study exercises the reduced path, but the architectural
boundary is the common numerical contract layer, not one specific optimizer. Other
formulation products can also be constructed directly from application-owned
numerical actions, but each has its own required actions and assumptions; one
`ExecutableModelT` binding does not automatically provide KKT, OTD, or PDAS.

## 5. A concrete case: deal.II Step-4

The repository uses deal.II's tutorial program `step-4` as a deliberately small
external application.

Step-4 already knows how to solve a Poisson problem. It creates the mesh, distributes
finite-element degrees of freedom, assembles a stiffness matrix and load vector,
applies Dirichlet data, runs conjugate gradients, and writes the field.

The adapted version exposes enough seams to reuse that work repeatedly from an
optimization problem. The application remains the owner of its numerical data and
solver.

Two optimal-control problems are then built around it. The richer **Problem B**
introduces a distributed finite-element control while preserving fixed physical
state boundary values.

The physical state is reconstructed as

$$
y_{\mathrm{phys}} = Pz+\ell,
$$

and the control coupling uses the finite-element mass matrix $M$ through

$$
B=P^{\mathsf T}M.
$$

The direct binding exposes the resulting residual/objective actions, state and
adjoint solves, and control metric to the common formulation.

The detailed coefficient dimensions and assembly formulas are useful in the case
study, but they are not requirements of the general integration architecture.

## 6. Preserve the application's coordinate model

A recurring integration mistake is to assume that the vector optimized by the
solver must be identical to a full finite-element field.

That need not be true.

An application may use:

- free-state coordinates plus a fixed boundary lifting;
- a control space with a different dimension from the state;
- a boundary trace representation;
- a cellwise-discontinuous control;
- application-specific reconstruction for output.

The binding should preserve those meanings explicitly rather than flattening
everything into one interchangeable coefficient vector.

This is why the common contracts distinguish layouts, primal values, covectors, and
metrics.

## 7. Residual evaluation and solving remain separate

An application may already have an efficient state solver without exposing a useful
general residual API, or vice versa.

The common contract distinguishes:

```text
evaluate E(y,u)
```

from

```text
given u, solve E(y,u)=0 for y
```

and similarly distinguishes a transpose action from an adjoint solve.

That separation lets the formulation express the mathematics cleanly while the
application chooses how each operation is implemented.

One current ergonomic limitation is that the broad executable-model contract asks
for residual and Jacobian actions even when a particular reduced first-order path
does not call all of them during optimization. That is a contract-design issue to
evaluate separately; it is not a reason to move PDE ownership into the optimizer.

## 8. Lifetime and borrowing are part of the integration

The callback boundary is generic, but the underlying numerical objects are not
necessarily cheap or freely movable.

deal.II matrices, sparsity structures, DoF handlers, triangulations, and solver
services may have ordering and lifetime constraints.

A safe integration therefore keeps the native application, problem realization, and
binding in a scope where borrowed references remain valid for the lifetime of the
formulation and optimizer.

The useful rule is the same as in the compiled path:

> Type erasure is a solver boundary, not an ownership boundary.

## 9. The metric belongs to the optimization model, even when the application computes it

Suppose a reduced derivative is $r\in U^{\ast}$.

The optimization gradient $g\in U$ is defined by the selected metric $G$:

$$
Gg=r.
$$

For a finite-element $L^{2}$ control, $G$ may be the mass matrix. For a coefficient
Euclidean geometry it may be the identity.

The application can own and solve with $G$ while `nmopt` consumes it through the
metric interface.

This keeps a crucial distinction visible:

- the application owns the numerical operator;
- the mathematical model chooses what that operator means;
- the optimizer uses the declared metric without knowing how it was assembled.

## 10. Output remains native

After optimization, the generic solver may have a final control and a retained
state evaluation.

Writing a meaningful PDE result often still requires:

- reconstructing fixed boundary values;
- mapping reduced coordinates to a physical field;
- attaching cell, point, or boundary data to the mesh;
- using the application's own writer.

That belongs in the application/native layer.

The Step-4 integration therefore writes the final field through Step-4 rather than
teaching the generic optimizer how to understand a `DoFHandler`.

## 11. What should be in the binding, and what should not

A good binding is narrow.

It should contain:

- layout/partition construction;
- wrappers from native vectors to contract values;
- callbacks to objective/residual derivative operations;
- state and adjoint solve adapters;
- metric/constraint adapters;
- lifetime-preserving references to the native application.

It should not contain:

- a second PDE discretization;
- a benchmark-specific optimization algorithm;
- duplicated mesh ownership;
- application-independent solver policy that belongs in `nmopt::solvers`.

If the binding becomes a second application, the ownership boundary has probably
been crossed in the wrong direction.

## 12. Relationship to the semantic/compiler path

The semantic/compiler path and direct-application path are peers.

```text
semantic ProblemSpec ──► compiler ───────────┐
                                             │
existing PDE application ──► thin binding ───┤
                                             ▼
                                  common numerical services
                                             │
                                             ▼
                                  selected formulation
```

The compiler path is valuable when `nmopt` should create the numerical realization.
The direct path is valuable when the realization already exists.

The downstream formulation should not care which path was used.

## 13. Relationship to experiments and benchmarks

The external Step-4 case is currently organized as a case study under
`apps/external-dealii/step-4/`, with its explanatory overview beside the code.

The Chapter 6 B-series work, by contrast, uses the repository's experiment and
replication layer because it needs parameter families, source-oriented run matrices,
manifests, post-processing, and benchmark evidence.

This difference in documentation location reflects **use**, not architectural
importance:

- the Step-4 material demonstrates how an independently owned application connects;
- the B-series material demonstrates how framework-managed problems are executed and
  evaluated as reproducible experiments.

This overview provides the general direct-integration explanation. The Step-4 page
can therefore remain next to the case-study source without making the direct path
look like a special or lesser architecture.

## 14. A practical integration checklist

For a new application:

1. identify the optimization variables and the coordinates used by the native code;
2. expose repeatable state and, when required, adjoint solves;
3. implement or expose objective and derivative actions;
4. decide which residual/JVP/VJP actions the selected formulation requires;
5. choose the optimization-space metric and constraints;
6. construct layouts that make incompatible spaces explicit;
7. bind the native operations without transferring numerical ownership;
8. run the selected formulation and algorithm;
9. reconstruct and write results through the native application.

If several applications later share the same mathematical structure and numerical
realization strategy, that may justify adding a semantic/compiler family. It is not
a prerequisite for the first integration.

## Read later: authoritative sources

For the concrete worked case, read
[How `nmopt` connects to an existing deal.II application](../../../apps/external-dealii/step-4/external-integration-overview.md).

For exact public calling conventions, read
[External deal.II solver integration](../../reference/external-dealii-solver-integration.md).

For the accepted ownership rules, read
[PDE, formulation, and solver boundary](../../design/pde-solver-boundary.md).

For the mathematical role of layouts, covectors, metrics, and adjoints, read
[Theoretical formalism](../../design/mathematical-model.md).

For the peer framework-managed route, continue with
[Describing and compiling a problem](semantic-compiler.md).
