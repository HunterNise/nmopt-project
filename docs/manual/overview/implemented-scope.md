# Implemented scope

This page summarizes what `nmopt` implements and how strongly each part is exercised.

It is an overview, not an API reference or implementation-status ledger. For exact
public interfaces, use [`docs/reference/`](../../reference/). For source ownership
and focused verification, use the
[implementation map](../../internals/implementation-map.md). Current source and
tests remain authoritative for implementation behavior.

## Scope at a glance

`nmopt` is a C++17 framework for discretized PDE-constrained optimization, with
deal.II as its principal finite-element environment.

Its central separation is:

```text
PDE realization
      │
      ▼
optimal-control formulation
      │
      ▼
optimization / KKT / active-set solver
```

Numerical PDE services can be produced in two ways:

```text
semantic/compiler path                 application-owned path

ProblemSpec                            existing PDE application
    │                                          │
    ▼                                          ▼
validation / resolution                numerical adapters
    │                                          │
    ▼                                          │
deal.II compiler                               │
    │                                          │
    └──────────────────────┬───────────────────┘
                           ▼
                shared numerical contracts
                           │
                           ▼
           reduced / OTD / KKT / PDAS products
                           │
                           ▼
                         solvers
```

The two paths meet at numerical and formulation contracts. An existing
application does not need to become a `ProblemSpec` or pass through the deal.II
compiler.

The strongest demonstrated part of the project is **serial, scalar,
steady/elliptic finite-element optimal control**. The semantic and numerical
interfaces extend beyond that description in several directions, but the project
does not provide a general compiler for arbitrary PDE-constrained optimization
problems.

## Numerical contracts

The reusable core distinguishes mathematical roles that are easily blurred after
discretization.

It provides:

- block layouts and space identities;
- separate primal and covector values;
- dual pairings;
- residual, JVP, and VJP actions;
- objective values and derivatives;
- separate state and adjoint solve services;
- metrics and metric inversion;
- projection constraints;
- optional reduced-Hessian actions; and
- structured solve and optimization evidence.

A derivative remains a covector until a metric maps it to a primal gradient. State
and adjoint inversion are also separate from residual and transpose actions.

These contracts are backend-parametric. The repository supplies a reference dense
backend and a serial deal.II backend.

## Formulations and optimization

The formulation layer does not force every problem through one universal
formulation type.

### Reduced state–adjoint formulation

The reduced path is the most extensively used formulation in the repository.

It supports:

- retained value evaluation followed by derivative augmentation;
- steepest descent;
- nonlinear conjugate-gradient methods;
- full BFGS and L-BFGS;
- Newton when an explicit reduced-Hessian capability is available;
- Armijo, fixed-step, weak-Wolfe, strong-Wolfe, and exact-quadratic
  globalization where compatible; and
- unconstrained trust-region optimization with an explicit reduced Hessian.

Projection constraints are supported through a separate `ConstraintT` capability.
The current projected line-search surface is deliberately narrower than the
unconstrained one: projected steepest descent is supported with the registered
line-search policies, while projected nonlinear CG, BFGS, Newton, exact-quadratic
search, and constrained trust region are not general supported combinations.

See the [optimization reference](../../reference/optimization.md) for the exact
policy combinations.

### Supplied optimize-then-discretize systems

A supplied OTD product owns explicit state, adjoint, and control-stationarity
blocks supplied by the application author.

The framework can execute the registered canonical scalar supplied-OTD system and,
where the required evidence is declared, convert the canonical system to a
quadratic KKT product.

The framework does **not** derive a continuous adjoint or OTD system automatically
from an arbitrary strong PDE description.

### Quadratic KKT products

The reusable KKT product represents equality-constrained quadratic systems through
typed primal, multiplier, stationarity, and equality spaces together with the
required operator and transpose actions.

The current serial path includes MINRES and GMRES policies and independent
stationarity/equality residual reporting.

Compiler-produced KKT products are currently registered for the canonical scalar
DTO target and the canonical compatible supplied-OTD target. A general reduced
problem does not automatically acquire a KKT representation.

### Complementarity and PDAS

The framework includes typed box complementarity and a primal-dual active-set
solver built on restricted KKT subproblems.

The compiler-backed PDAS path is currently bounded to the registered
cellwise-discontinuous control box with a compatible positive-diagonal cellwise
$L^{2}$ metric. It records active sets, primal and dual feasibility,
complementarity, stationarity, equality residuals, and inner KKT evidence.

Continuous-control and other bound representations are not automatically converted
to this PDAS path.

## Semantic problem descriptions and deal.II realization

The semantic layer describes mathematical structure independently of deal.II
objects. Validation checks graph structure and declared policies; compiler
lowerability is a separate question.

The deal.II compiler is intentionally **registered rather than universally
compositional**. A semantically valid graph may therefore be unsupported by the
current numerical realization.

The current registered scalar realization surface includes the following families.

### Volume control

Implemented scalar volume-control paths include combinations of:

- diffusion and diffusion–reaction equations;
- tensor diffusion;
- conservative or advective transport;
- Robin terms;
- homogeneous or fixed Dirichlet state data;
- full-domain and material-subdomain tracking;
- $L^{2}$ and $H^{1}$ state tracking;
- cellwise `FE_DGQ(0)` controls;
- continuous conforming volume controls;
- $L^{2}$, $H^{1}$, and selected $H^{-1}$ metric/regularization combinations;
- cellwise boxes where the control and metric realization support their declared
  coefficient semantics; and
- positive cellwise diffusion-coefficient identification.

### Boundary control

Registered boundary-control paths include:

- facewise-constant Neumann control;
- continuous degree-one Neumann traces;
- complete and selected partial Dirichlet-control transformations;
- $L^{2}$, $H^{1/2}$, and $H^{1}$ trace-related control metrics and losses; and
- the explicit boundary partition and transport conventions required by those
  targets.

Optional boxes are tied to specific realizations. A boundary or continuous control
does not gain box semantics merely because its coefficients are scalar.

### Observations

Specialized registrations include:

- volume and subdomain tracking;
- boundary tracking;
- point sensors; and
- normal-flux observations.

Point and flux observations use explicit transposition policies. They are not
generic plugins that can be attached to every compiler target.

### Mesh families

Selected continuous-control and Neumann-control registrations support both
hypercube `FE_Q` and simplex `FE_SimplexP` realizations.

Other registered targets remain hypercube paths. Mixed reference-cell meshes are not
a general fallback.

For the exact supported combinations and diagnostics, see the
[compiler reference](../../reference/compiler.md).

## Application and execution layer

The repository contains a reusable application layer around the numerical core.

Implemented pieces include:

- typed problem recipes;
- typed application scenarios;
- scenario validation and metadata;
- runtime-data definitions;
- backend execution adapters;
- a repository application catalog;
- deal.II-style `.prm` parameter files for the registered application slice;
- matrix expansion and exclusions;
- a headless repository runner;
- deterministic run directories and manifests;
- native deal.II mesh and field output;
- solver traces and numerical evidence records;
- Python post-processing;
- comparison plots; and
- deterministic benchmark reports.

The runner and experiment machinery are optional outer infrastructure. They are not
required to use the core numerical/formulation interfaces.

The current registered repository-runner slice is centered on the Chapter 6 B1 and
B2 applications. The compiler supports more mathematical families than are exposed
as complete application recipes and benchmark scenarios.

See [application authoring](../../reference/application-authoring.md),
[application execution](../../reference/application-execution.md), and
[parameter files](../../reference/parameter-files.md).

## Existing-application integration

An existing deal.II application can bypass the semantic/compiler producer entirely.

The application keeps ownership of:

- its mesh;
- finite-element spaces;
- assembly;
- boundary treatment;
- linear solves; and
- native output.

It exposes the residual/objective actions and state/adjoint solve services required
by the selected formulation.

The Step-4 integration study verifies this boundary for two serial linear problems.
Problem B additionally exercises independent versus physical state coordinates,
Dirichlet lifting, a rectangular finite-element control coupling, a nonidentity mass
metric, and reconstruction back to the native field representation.

This is bounded evidence that the ownership separation works. It is not evidence
that arbitrary nonlinear, distributed, adaptive, or constrained applications can be
integrated with the same effort.

See the
[external integration reference](../../reference/external-dealii-solver-integration.md)
and the
[Step-4 case study](../../../apps/external-dealii/step-4/README.md).

## Application and study evidence

Implementation and end-to-end demonstration are not the same thing.

### B1

The Chapter 6 B1 distributed-control study exercises the semantic/compiler producer,
reduced state–adjoint formulation, reduced optimization, native field output,
post-processing, and reproduction machinery.

Its source-oriented run is reproduction-verified under an explicit set of project
replacement choices for information omitted by the source.

See [B1 replication](../../studies/chapter-6/b1-replication.md).

### B2

The Chapter 6 B2 boundary-control study exercises transport, boundary control,
subdomain observations, derivative checks, native boundary output, and the same
application/reproduction layer.

The declared problem is framework-verified, but the literal source formulation does
not reproduce the published source fields/results. The repository preserves that as
a negative replication result rather than silently modifying the mathematical
problem.

See [B2 replication](../../studies/chapter-6/b2-replication.md).

### Step-4

The Step-4 study exercises the application-owned producer path and verifies that a
pre-existing deal.II code can retain numerical ownership while using the common
reduced formulation and optimizer.

### KKT, supplied OTD, and PDAS

These formulation families have focused algebraic, compiler, deal.II, and numerical
tests, including real serial KKT and active-set solves.

They do not have B1/B2-scale completed reproduction campaigns. Their implementation
evidence is therefore strong at the contract and integration-test levels but
narrower at the application-study level.

## Important limits

The current project does not implement or demonstrate all directions suggested by
its semantic model or retained roadmaps.

In particular:

- the deal.II execution path is serial;
- there is no general MPI/distributed backend;
- there is no general time-dependent/parabolic compiler path;
- mixed multi-equation systems such as Stokes are outside the implemented compiler
  scope;
- general state constraints and measure-valued multipliers are outside scope;
- automatic OTD derivation is not implemented;
- GLS and stabilized-Lagrangian formulation families are not implemented as a
  general path;
- the compiler does not automatically lower every semantically valid graph;
- PDAS box semantics are narrower than the set of available control
  representations;
- Chapter 6 B3–B6 are not registered completed benchmarks;
- the broader Chapter 5 compiler surface is not yet mirrored by a complete
  user-facing recipe/catalog library;
- no desktop GUI is implemented; and
- the CMake project does not currently install/export an `nmopt` package for
  downstream `find_package(nmopt)` use.

Possible continuation directions are collected separately in
[Future extensions](../../planning/future-extensions.md).
