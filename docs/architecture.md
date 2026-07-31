# Architecture record: composable PDE-constrained optimization

## Purpose and status

This document records the high-level direction agreed for `nmopt-project`.
It is intentionally architectural rather than an implementation specification.
It should guide future design discussions and prevent the project from growing
into a collection of tightly coupled, problem-specific classes.

The target implementation library is deal.II. The framework is not limited to
one PDE or one optimal-control formulation, but the first working version must
be deliberately small and well tested.

## The problem family

The general semantic model is a constrained optimization problem

$$
  \min_{x\in X_\mathrm{ad}} J(x)
  \qquad\text{subject to}\qquad
  E(x)=0\quad\text{in }Z^{*}.
$$

The variable block $x$ may contain a state $y$, control $u$, parameter
$m$, and further variables. The residual may map trial variables into the
dual of a different test space. This admits ordinary weak, mixed, and, with a
specialized realization, very weak formulations.

For a state/control model and the Lagrangian convention

$$
  \mathcal L(y,u,p)=J(y,u)-\langle p,E(y,u)\rangle_{Z,Z^{*}},
$$

the reusable first-order pattern is

$$
  E(y,u)=0,
  \qquad E_y(y,u)^{*}p=J_y(y,u),
  \qquad j'(u)=J_u(y,u)-E_u(y,u)^{*}p\in U^{*}.
$$

This is the central abstraction. In the linear-quadratic case it reduces to
constant state, control, observation, and metric operators. In coefficient
identification or nonlinear problems, the same relations hold but derivative
operators depend on the current point.

## Derivative versus gradient

A derivative is naturally a dual object. It must not be treated as a control
vector without declaring a primal-dual identification. In a Hilbert space a
chosen metric provides a Riesz map

$$
 R_U:U\to U^{*},\qquad \nabla_U j=R_U^{-1}j'(u).
$$

Thus changing $L^2$ to $H^1$ regularization/gradient changes the metric
realization, not the state equation or adjoint derivation. An $H^{-1}$ metric may
require an auxiliary elliptic solve. Fractional and trace norms require an
explicit, documented discrete realization; they are not generic strings or
interchangeable mass matrices.

At the semantic level the framework may describe Banach spaces and pairings. At
the executable level, every optimization method must have a concrete discrete
way to map the required dual derivative into a search direction. This can be a
Hilbert metric, a duality map, or another explicit algorithmic geometry.

## Composition model

Do not represent whole PDE/control combinations by subclasses. A user builds a
`ProblemSpec` from independently meaningful components. Components contribute
through a limited set of ports:

| Port | Meaning |
|---|---|
| Residual | Add a contribution to $E(x)\in Z^{*}$. |
| Derivative | Supply Jacobian and transpose-Jacobian actions needed by adjoints. |
| Objective | Add a scalar term and its partial derivatives. |
| Metric | Define a primal-dual pairing and, when needed, an inverse map. |
| Transformation | Parameterize variables via lifting, restriction, trace, or transfer. |
| Constraint | Define feasibility, projection, multiplier, or active-set behavior. |
| Space/geometry | Declare fields, components, regions, boundary tags, and requirements. |

Typical components are:

```text
Residual terms: diffusion, transport, reaction, load, Robin, stabilization,
                mixed coupling, time derivative
Control couplings: volume source, Neumann/Robin boundary, actuator,
                   Dirichlet lifting, parameter/coefficient action
Observations: identity, subdomain restriction, boundary trace, flux, sensor
Objectives: tracking, regularization, penalty
Metrics: L2, H1, boundary L2, auxiliary-solve H-1, declared fractional metric
Constraints: box, affine, mean/gauge, complementarity/active set
```

Adding a valid feature should mean implementing its own local contributions,
not modifying every solver or creating a new complete problem type.

## Boundary-condition semantics

Boundary conditions must not be flattened into one generic load interface.

| Kind | Architectural effect |
|---|---|
| Essential Dirichlet | Restricts or parameterizes the trial space; compile to constraints and/or a lifting. |
| Natural Neumann | Adds a boundary residual/load contribution. |
| Robin | Adds a boundary bilinear contribution and possibly a load. |
| Periodic/hanging-node | Adds algebraic DoF constraints. |
| Pure Neumann | Requires compatibility and an explicit nullspace/gauge policy. |

For Dirichlet control, use $y=\ell(u)+\hat y$ with homogeneous unknown
$\hat y$. The control derivative must include the chain rule through
$\ell$. It is not generally equivalent to an ordinary boundary source term.

## Concrete reference: Laplace source control

The basic model is

$$
  -\Delta y=f+u\ \text{in }\Omega,
  \qquad y=0\ \text{on }\partial\Omega,
$$

with $`y\in H_0^1(\Omega)`$, $u\in L^2(\Omega)$, and a distributed $L^2$
tracking and $L^2$ regularization objective. Its residual is

$$
 \langle E(y,u),v\rangle=(\nabla y,\nabla v)-(f,v)-(u,v).
$$

The recipe consists of a state field, diffusion term, load term, source-control
coupling, Dirichlet condition, tracking term, and regularization metric. The
adjoint and control derivative emerge from the generic residual relations:

$$
 -\Delta p=y-y_d,
 \qquad j'(u)=\alpha u+p.
$$

Replacing source control by Neumann control changes the control-residual term
and its transpose into a boundary trace. Replacing it by Dirichlet control
changes the variable parameterization via a lifting. Changing distributed
tracking into boundary tracking changes only the observation map and its
transpose. Changing the metric changes only the dual-to-primal conversion.

## Semantic-to-executable pipeline

```text
1. Geometry and named regions
2. ProblemSpec: fields, residual terms, boundary conditions, objectives,
   metrics, constraints
3. Semantic validation: resolve fields/regions and check declared requirements
4. Discretization policy: FE choices, quadrature, backend, assembled/matrix-free
5. Compilation: DoF layouts, AffineConstraints, liftings, discrete maps
6. Executable problem: residual/Jacobian/transpose, objective, metrics
7. Formulation: reduced-space, KKT/all-at-once, active-set/complementarity
8. Solver: state/adjoint and optimization/Krylov/preconditioning strategies
9. Output: fields, diagnostics, objective history, active sets, metadata
```

The semantic layer contains no deal.II objects. The discrete compiler lowers a
semantic recipe to deal.II `DoFHandler`/FE layouts, `AffineConstraints`,
assembly or matrix-free operators, distributed vectors, and transfer/lifting
maps. Algorithms consume only the executable operator contract.

## What the generic core should require

An executable model needs to supply, directly or through composition:

```text
E(x)                         residual
E'(x) dx                     linearized residual action
E'(x)* p                     transpose/adjoint action
J(x), J'(x)                  objective and derivative
R v, R^-1 ξ                  metric and inverse metric where required
constraint operations        projection, active set, or KKT relation
```

For linear-quadratic elliptic models, these actions can be represented by
cached matrices. For nonlinear, coefficient-dependent, and time-dependent
models, they can be assembled at each iterate or provided matrix-free.

## Extensibility map

| Requested variation | Local extension | Core layers left unchanged |
|---|---|---|
| diffusion/reaction/transport change | residual term | objective, optimization workflow |
| volume or Neumann control | control coupling and transpose | state/adjoint driver, optimizer |
| boundary/subdomain tracking | observation and metric | PDE and optimizer |
| $H^1$ or $H^{-1}$ search geometry | metric realization | residual and adjoint derivation |
| $H^1$ or $H^{-1}$ regularization | loss and control-space realization | residual and adjoint derivation |
| Dirichlet control | lifting/parameterization | formulation and outer solver |
| coefficient determination | parameter residual derivatives | adjoint orchestration and line search |
| vector/mixed state | field blocks and residual blocks | conceptual model and optimizer |
| time dependence | evolution residual and time adapter | objective/metric contract |
| box constraints | projection or active-set component | state residual |

## Limits and validation

The framework must be explicit about where generic composition ends. It cannot
prove well-posedness for arbitrary combinations or infer a unique FE treatment
for every functional-analytic space. Each component declares requirements, and
the validator rejects invalid combinations or requires an explicit policy.

Examples:

```text
point observation       requires regularity or a declared discrete-only sensor
fractional norm         requires a selected discrete realization
pure Neumann operator   requires compatibility and a nullspace treatment
Dirichlet control       requires a lifting, not a load coupling
coefficient control     is nonlinear despite linearity in the state
very weak formulation   requires trial/test-space and pairing realization
```

This is a feature: it keeps subtle mathematical assumptions visible instead of
hiding them in branches of unrelated classes.

## Initial vertical slice

Implement and test the following before attempting the advanced cases:

```text
scalar stationary diffusion-reaction equation
Dirichlet, Neumann, Robin, and pure-Neumann nullspace policy
volume source and Neumann boundary controls
distributed and boundary-trace observations
L2 and H1 metrics
box-constrained and unconstrained controls
reduced-space state/adjoint/gradient workflow
Armijo gradient method, nonlinear CG or L-BFGS
adjoint-consistency, Taylor-remainder, and manufactured-solution tests
```

Future extensions should preserve the same contracts. In particular, do not
add PDE-specific conditions to optimizers or solver-specific behavior to weak
form terms.
