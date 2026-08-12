# Chapter 6 benchmark suite roadmap

## Purpose

This roadmap owns the frozen numerical experiments from Chapter 6. It is
separate from the [Chapter 6 numerical-methods guide](../guides/chapter-6-numerical-methods.md),
which records reusable numerical-method contracts, and the
[Chapter 6 numerical-examples reference](../guides/chapter-6-numerical-examples.md),
which records source equations and data.

The benchmark suite is not a second problem library. Each benchmark selects a
specific configuration from the Chapter 5 problem library, freezes its data
and discretization, and records expected numerical behavior. Its purpose is
system-level validation and source-result reproduction.

Benchmark results should be interpreted in layers:

1. exact residual, objective, JVP, VJP, reduced-Taylor, and KKT tests establish
   local mathematical correctness;
2. the benchmark establishes that the compiled components cooperate on a
   source-sized problem and reproduce objective, field, feasibility, and
   convergence trends;
3. timings and iteration counts are performance observations, not universal
   correctness tolerances.

## Selection

| Priority | Example | Decision | Purpose |
| --- | --- | --- | --- |
| 1 | E6.5.1 distributed Laplace | Selected | Reduced-space/L-BFGS calibration on a cheap linear-quadratic problem. |
| 2 | E6.5.2 Graetz boundary control | Selected | Flagship scalar boundary-control benchmark. |
| 3 | E6.9.1 symmetric box Laplace | Desirable | PDAS and complementarity validation with a manufactured active set. |
| 4 | E6.9.2 asymmetric box Laplace | Desirable | Spatially varying bounds after the symmetric PDAS case. |
| 5 | E6.7.1 all-at-once Laplace | Conditional | KKT and optional preconditioning benchmark. |
| – | E6.7.2 diffusion-reaction all-at-once | Follow-up | Regression variant of E6.7.1. |
| – | E6.5.3 reduced Stokes | Excluded | Deferred mixed-block/Section 5.13 work. |
| – | E6.9.3 box-constrained Stokes | Excluded | Deferred Stokes plus vector active-set work. |

## Common benchmark contract

Every benchmark record must freeze:

- the problem-library recipe and its version;
- domain, mesh generator, mesh size, finite elements, and boundary labels;
- coefficient, forcing, target, observation, and bound data;
- formulation policy and all state/control/observation spaces;
- initial control and all algorithm parameters;
- state, adjoint, metric, KKT, and stopping tolerances;
- output fields, objective components, gradient/KKT histories, and residuals.

Known omissions in the book must be recorded as explicit replacements. In
particular, do not silently choose unspecified forcing, line-search rules,
viscosity, pressure gauges, or objective scaling.

Runs should use `release-dealii`. Development uses smaller meshes; source-
sized meshes are the reproduction gate. The first performance evidence is
assembled-operator reuse, factorization reuse, objective-only trial
evaluations, and separate solve counts. Advanced preconditioning is not a
prerequisite for the selected reduced-space benchmarks.

## Benchmark sequence

### B0 — Common experiment harness

**Dependency:** Chapter 5 remediation, the selected P6.1 policies, and the
problem-library recipe boundary.

Implement a thin runner that consumes a compiled problem and a frozen
scenario. It should provide deterministic manifest/report output, mesh and
parameter selection, history serialization, and timing/memory measurements.
It must not contain a second PDE lowerer or optimizer implementation.

### B1 — E6.5.1 distributed Laplace control

Use the Chapter 5 distributed-control recipe with the source's square-domain
target and $\beta$ sweep. Compare steepest descent and L-BFGS from the same
initial control.

Required evidence:

- objective and tracking reduction;
- iteration, state-solve, adjoint-solve, and line-search counts;
- reduced-gradient histories;
- finite-difference verification of the linear-quadratic Hessian-vector
  action;
- qualitative agreement with the source trend that L-BFGS is much faster
  than steepest descent.

The source text does not fully specify $f$. Recover it from the cited source or
record a manufactured replacement; do not compare absolute objective values
across different choices.

**Dependencies:** P6.1. No KKT or preconditioner is required.

### B2 — E6.5.2 Graetz-flow boundary control

Use the Neumann boundary-control recipe with scalar transport, fixed
temperature, controlled boundary flux, insulated outflow, and downstream
observation regions. Run all four source combinations of observation region
and target.

Required evidence:

- correct volume-state/face-control dimensions and manifest;
- transport and inflow/outflow policy;
- boundary-control residual JVP/VJP and reduced Taylor test;
- objective and relative-gradient reduction;
- state/control field behavior when the observation region changes.

The first reproduction is the stated Galerkin formulation. GLS and other
stabilization are explicitly outside this benchmark suite. If the unstabilized
discretization is numerically inadequate, record the limitation instead of
silently changing the benchmark formulation.

The book does not fully specify the BFGS step policy. Declare a recovered
fixed-step or Armijo policy in the scenario manifest.

**Dependencies:** Chapter 5 boundary-control remediation and P6.1.

### B3 — E6.9.1 symmetric box-constrained Laplace control

Freeze the manufactured target, regularization sweep, and symmetric bounds.
The first framework-native benchmark should use cellwise-discontinuous
control, where coefficientwise $L^{2}$ box semantics are exact. A continuous
`Q1` reproduction is conditional on a separate continuous-control bound
policy.

Required evidence:

- inactive-box agreement with the unconstrained solution;
- active-set stabilization at the expected bounds;
- primal feasibility, dual feasibility, complementarity, stationarity, and
  active-set-change histories;
- comparison with projected reduced optimization where applicable.

**Dependencies:** P6.3 and P6.5.

### B4 — E6.9.2 asymmetric box-constrained Laplace control

Extend B3 to spatially varying lower and upper bounds, beginning in two
dimensions. The three-dimensional source case is optional and must not block
the scalar framework benchmark.

Required evidence is the B3 KKT/PDAS record plus sensitivity to the bound
jump and regularization parameter. The KKT and PDAS services are reused; only
the bound-data recipe and scenario change.

**Dependencies:** B3 and spatially varying bound data in the problem library.

### B5 — E6.7.1 all-at-once Laplace control

This is conditional because it opens a second all-at-once/Krylov direction.
The first target is a small scalar problem, not the full source mesh sweep.

Required evidence:

- assembled/operator-action agreement for the KKT product;
- agreement with the reduced DTO solution;
- independent feasibility, stationarity, multiplier-conversion, and solver
  diagnostics;
- source-consistent objective scaling, including the book's
  $\beta\lVert u\rVert^{2}$ convention;
- outer and inner work records if a preconditioner is activated.

Activate the bounded P6.4 preconditioner only after measuring whether direct
or basic serial solves are inadequate. Begin with one scalar block-diagonal
mass/PDE approximation; do not make mesh-robustness claims without a sweep.

**Dependencies:** P6.3. P6.4 is conditional.

### B6 — E6.7.2 diffusion-reaction follow-up

After B5, reuse the same all-at-once implementation and replace the state
stiffness action by the diffusion-reaction operator. This is a regression
benchmark, not a new application or solver path.

## Exclusions

The benchmark suite excludes:

- E6.5.3 reduced Stokes control;
- E6.9.3 box-constrained Stokes control;
- measure-valued state-constraint examples;
- GLS, stabilized-Lagrangian, and other stabilization comparisons;
- automatic derivation of continuous OtD systems;
- continuous-control box reproduction before its bound policy exists.

The numerical source catalogue remains useful as a reference even where a
benchmark is excluded. Exclusion means that the experiment is not an
acceptance target for the current project scope.
