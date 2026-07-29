# Project guide

## Goal

Build a deal.II-based framework for PDE-constrained optimal control and inverse
problems. It should make valid variations composable: elliptic and later
evolution PDEs; scalar/vector/mixed fields; distributed, boundary, and
parameter controls; observations; norms; boundary conditions; box constraints;
and reduced-space or all-at-once solvers.

The detailed architectural record is [docs/architecture.md](docs/architecture.md).

## Core design decisions

1. Do **not** create subclasses such as
   `NeumannBoundaryControlProblem` or `DiffusionTrackingProblem`.
2. The general mathematical model is
   $\min J(x)$ subject to $E(x)=0$ in a declared test-space dual. Here
   $x$ can contain state, control, parameters, and other variables.
3. Components compose through explicit ports:

   - residual and residual derivatives;
   - objective and objective derivatives;
   - primal/dual pairings and metrics;
   - transformations such as liftings and restrictions;
   - constraints and admissible-set operations.

4. A derivative is a dual element. A gradient is obtained only after choosing
   a metric/Riesz or other primal-dual identification. Do not silently identify
   all derivatives with `L2` functions.
5. Dirichlet conditions are not ordinary load terms. Essential conditions
   restrict/parameterize the state space and may require a lifting. Dirichlet
   boundary control must go through such a lifting.
6. Keep continuous semantics and deal.II realization separate. A semantic
   model describes spaces, residuals, objectives, and requirements; a compiler
   creates FE layouts, `AffineConstraints`, liftings, operators, and vectors.
7. The core coordinates state, adjoint, gradient, KKT, and solver workflows.
   Individual terms own the weak form and adjoint/Jacobian actions they need.

## Target pipeline

```text
ProblemSpec -> semantic validation -> discretization/compilation
            -> executable residual/objective/metric operators
            -> reduced or KKT formulation -> solver -> output
```

## Scope discipline

Do not claim arbitrary combinations are automatically well posed or have a
canonical FE discretization. Components must declare requirements, and the
validator must reject or request an explicit policy for cases such as point
observations, fractional norms, pure Neumann nullspaces, very weak forms, and
Dirichlet control.

Prefer a small, tested first vertical slice before advanced cases. The proposed
initial slice is scalar stationary diffusion-reaction with standard boundary
conditions, volume and Neumann control, distributed/boundary tracking, `L2` and
`H1` metrics, box constraints, and a reduced-space optimizer.

## When implementing

- Preserve the residual-based abstraction; never add solver-specific branches
  to a PDE term or PDE-specific branches to an optimizer.
- Make every required linearized action and its transpose testable.
- Keep signs, dual pairings, and Lagrangian convention documented globally.
- Add adjoint-consistency and Taylor-remainder tests for each new coupling.
- Use deal.II implementation details only in the discrete layer.
