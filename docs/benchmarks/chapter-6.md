# Chapter 6 benchmark contracts

This document freezes the project choices needed to execute and accept the
Chapter 6 benchmark suite. It does not reproduce the book's equations or
figures, define application construction, or document runner and
post-processing APIs.

The authorities are deliberately separate:

- the [numerical-examples reference](../guides/chapter-6-numerical-examples.md)
  records what the book states and leaves unresolved;
- the [Chapter 6 application scenarios](../applications/chapter-6.md) define
  recipe selection, scenario construction, and runtime ports;
- the [benchmark-suite roadmap](../planning/chapter-6-benchmark-suite-roadmap.md)
  owns selection, order, dependencies, and suite-level gates;
- the [application roadmap](../planning/application-roadmap.md) owns runner,
  artifact, post-processing, and current execution status; and
- the [application execution and artifact reference](../reference/application-execution.md)
  defines the public artifact and execution interfaces.

A benchmark choice that cannot be inferred from the source record is recorded
here when it is specific to a frozen run. Reusable semantic or API choices
belong to the application or implementation contracts instead.

## Benchmark catalogue

The full catalogue is B0–B6. B0–B2 have frozen contracts in this document;
B3–B6 are reserved until their source-sized run matrices and acceptance
evidence are frozen.

| ID | Source/application | Contract state |
| --- | --- | --- |
| B0 | Common benchmark record boundary | Frozen below |
| B1 | [E6.5.1 distributed Laplace](../guides/chapter-6-numerical-examples.md#e651--distributed-laplace-control) / [scenario](../applications/chapter-6.md#b1--distributed-laplace-control) | Frozen below |
| B2 | [E6.5.2 Graetz flow](../guides/chapter-6-numerical-examples.md#e652--graetz-flow-boundary-control) / [scenario](../applications/chapter-6.md#b2--graetz-flow-boundary-control) | Frozen below |
| B3 | [E6.9.1 symmetric box control](../guides/chapter-6-numerical-examples.md#e691--symmetric-box-constrained-laplace-control) | Reserved; see [roadmap](../planning/chapter-6-benchmark-suite-roadmap.md#b3--e691-symmetric-box-constrained-laplace-control) |
| B4 | [E6.9.2 asymmetric box control](../guides/chapter-6-numerical-examples.md#e692--asymmetric-box-constrained-laplace-control) | Reserved; see [roadmap](../planning/chapter-6-benchmark-suite-roadmap.md#b4--e692-asymmetric-box-constrained-laplace-control) |
| B5 | [E6.7.1 all-at-once Laplace](../guides/chapter-6-numerical-examples.md#e671--all-at-once-laplace-control) | Reserved; see [roadmap](../planning/chapter-6-benchmark-suite-roadmap.md#b5--e671-all-at-once-laplace-control) |
| B6 | [E6.7.2 all-at-once diffusion-reaction](../guides/chapter-6-numerical-examples.md#e672--all-at-once-diffusion-reaction-control) | Reserved; see [roadmap](../planning/chapter-6-benchmark-suite-roadmap.md#b6--e672-diffusion-reaction-follow-up) |

## Common benchmark contract

Every activated benchmark must freeze or record in its run manifest:

- the scenario and recipe identity, source reference, and source-catalogue
  revision;
- every replacement for a source omission, with provenance and rationale;
- the run matrix, initial values, solver policy, stopping rules, and mesh
  policy;
- the realized dimensions, finite elements, regions, boundary labels, and
  runtime data provenance; and
- the required fields, objective/gradient/KKT histories, residual checks, and
  comparison quantities.

The benchmark artifact must preserve the detached compilation manifest,
solver report and policy snapshot, environment, benchmark choices, diagnostics,
and selected native fields. The artifact schema, deterministic serialization,
path layout, report generation, and post-processing behavior are defined by the
[application execution reference](../reference/application-execution.md)
and [application roadmap](../planning/application-roadmap.md#application-boundaries),
not repeated here.

Development runs may use a smaller explicitly named mesh or build profile.
They can validate orchestration and artifact shape but cannot close a
source-sized reproduction gate. Timing and iteration counts are observations;
they are not portable correctness tolerances.

## B0 — common benchmark record boundary

B0 is the common association and evidence boundary, not a PDE solver. A
completed run must produce one detached benchmark artifact and its deterministic
projection using the public `nmopt-benchmark-v1` writer. The record must retain:

- deterministic scenario, recipe, output, and source identity;
- framework, mesh, runtime-data, and environment provenance;
- the complete compilation manifest and validation diagnostics;
- solver policy, stopping reason, histories, solve counts, and measurements;
- benchmark-specific replacement choices and acceptance evidence; and
- the selected native field inventory.

B0 is contract-complete when repeated writes are byte-identical, invalid
identity/diagnostic/measurement/field data is rejected, failed matrix entries
remain visible with diagnostics, and path selection remains outside the
artifact writer. Execution status and the concrete run layout belong to the
[application roadmap](../planning/application-roadmap.md#current-handoff-state).

## B1 — E6.5.1 distributed Laplace control

The source problem, target, figures, and source omissions are recorded in the
[E6.5.1 source entry](../guides/chapter-6-numerical-examples.md#e651--distributed-laplace-control).
The development evidence and current source-parity assessment are recorded in
the [B1 replication findings](b1-replication.md).
The application scenario owns the recipe, graph, control representation,
runtime ports, and backend construction.

### Frozen benchmark choices

| Choice | Frozen value |
| --- | --- |
| Forcing replacement | Constant $f=0.5$, selected as the simplest balanced candidate across the source field extrema; this is not a recovered source fact. |
| Forcing provenance | `chapter-6.e6.5.1.source-oriented-constant-half-forcing` |
| Mesh policy | Regular 131 by 131 structured-simplex unit-square mesh, yielding 34,322 triangles and 17,424 vertices; this is the closest natural regular candidate, not the omitted source connectivity. |
| Control discretisation | Homogeneous-Dirichlet continuous `FE_SimplexP<2>(1)` on the selected triangles, with the continuous $L^{2}$ metric. |
| Linear solves | State and adjoint use identity-preconditioned serial CG with dimension-dependent iteration limits and tolerances `1e-12` relative / `1e-14` absolute; the control mass metric uses at most `1000` iterations with the same tolerances. |
| Regularisation and method matrix | Steepest descent at $\beta\in\{10^{-1},10^{-2},10^{-3}\}$ and L-BFGS at those values plus $10^{-6}$, matching the seven unique source figure cases. |
| Methods | Steepest descent and metric-inverse limited-memory BFGS with memory 5, both from the same zero control. |

The benchmark uses the current B1 scenario's declared finite elements, control
layout, and runtime bindings. The source-oriented forcing, mesh, and stopping
choices are explicit evidence-backed replacements for omitted details; they do
not claim parity with the undisclosed source realization.

The development family `continuous-control.prm` changes only the original
framework-native control realization to independent homogeneous-Dirichlet
continuous `FE_Q` while retaining the $L^{2}$ loss and search metric. It
records the earlier current-mesh isolation of the source's continuous `P1`
clue.
The `continuous-control-constant-one.prm` family additionally tests the
$f\equiv1$ hypothesis inferred from the Figure 6.2 extrema. Its provenance
records that inference; it does not promote the omitted forcing to a source
fact.

Two further development families isolate the triangular connectivity
hypothesis while retaining manufactured-zero forcing, continuous P1 control,
and the recovered Figure 6.3 solver policy. The
`continuous-control-structured-simplex.prm` family uses a regular 131 by 131
subdivision, yielding 34,322 triangles, 17,424 vertices, and 16,900 independent
homogeneous-Dirichlet P1 coordinates. This is the closest regular square-grid
candidate to all three published counts. The
`continuous-control-count-matched-simplex.prm` family instead starts from a
100 by 100 subdivision and deterministically splits 7,160 triangles at their
centroids. It exactly yields the published 34,320 triangles, 17,361 vertices,
and 16,961 independent coordinates, but its seed-zero connectivity is an
explicit hypothesis rather than recovered source information.

### Solver policy

| Method | Stopping rule | Relative threshold | Armijo fraction | Backtracking | Maximum reductions | Operative minimum step |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Steepest descent | Relative gradient norm | $10^{-3}$ | $10^{-5}$ | `0.7` | `5` | `0.2` |
| L-BFGS | Relative gradient norm | $10^{-3}$ | $10^{-5}$ | `0.7` | `5` | `0.01` |

The common relative threshold selects the same L-BFGS iterates `2/4/4` as the
recovered Figure 6.3 objective-target policy and selects iteration 4 for the
$\beta=10^{-6}$ Figure 6.2 fields. The source states this threshold only for
steepest descent, so its use for L-BFGS is an explicit unifying inference.

The development family `figure-6.3-book-policy.prm` encodes the recovered
Figure 6.3 conventions separately: five Armijo reductions (six possible
trials), operative floors `0.2` and `0.01`, steepest-descent termination at a
relative gradient norm of $10^{-3}$, and L-BFGS termination at the terminal
steepest-descent objective for the same regularisation value. The L-BFGS
absolute gradient threshold is only a positive solver-contract safeguard. The
runner executes the reference method first and records the exact artifact that
supplied each dependent objective target.

### Acceptance evidence

The sparse matrix contains seven artifacts: the six Figure 6.3 combinations
and the L-BFGS $\beta=10^{-6}$ Figure 6.2 combination. Each artifact must
retain:

- objective and tracking histories and final values;
- gradient, step, line-search, state, adjoint, metric, and direction counts;
- the compilation manifest and executable, physical, and independent
  state/control/adjoint dimensions; and
- a finite-difference check of the linear-quadratic Hessian action.

The expected qualitative comparison is decreasing cost/tracking as
regularisation decreases and substantially less optimization work for L-BFGS.
This trend is evidence, not a portable numerical tolerance; a failure is
reported as a benchmark limitation.

Execution and acceptance status are maintained in the [application roadmap](../planning/application-roadmap.md#b1),
not duplicated in this freeze.

## B2 — E6.5.2 Graetz-flow boundary control

The source equation, four observation/target combinations, boundary geometry,
and source omissions are recorded in the [E6.5.2 source entry](../guides/chapter-6-numerical-examples.md#e652--graetz-flow-boundary-control).
The development evidence, count-based deductions, candidate experiment order,
and current source-parity assessment are recorded in the
[B2 replication findings](b2-replication.md).
The application scenario owns the recipe, case enum, material-region
realization, runtime ports, and backend construction.

### Frozen benchmark choices

| Choice | Frozen value |
| --- | --- |
| Run matrix | All four public `GraetzCase` values from `make_catalog()`. |
| Forcing | Zero forcing, with provenance `chapter-6.e6.5.2.zero-forcing`. |
| Regularisation and initial value | $\beta=10^{-3}$ and zero facewise control. |
| Method | Full BFGS with the solver policy below. |
| Mesh policy | Framework-native rectangle, source-sized `refine_global(7)`; realized labels and dimensions remain manifest data. |
| Stabilization | Stated Galerkin formulation only; GLS and other stabilization policies are excluded. |

The four stable output IDs are:

| Case enum | Scenario ID |
| --- | --- |
| `observation_wings_constant_target` | `chapter-6.b2.graetz-flow` |
| `observation_full_constant_target` | `chapter-6.b2.graetz-flow.full-constant` |
| `observation_wings_parabolic_target` | `chapter-6.b2.graetz-flow.wings-parabolic` |
| `observation_full_parabolic_target` | `chapter-6.b2.graetz-flow.full-parabolic` |

The frozen runtime and region metadata are:

| Record | Frozen value |
| --- | --- |
| Forcing provenance | `chapter-6.e6.5.2.zero-forcing` |
| Desired-state provenance | `chapter-6.e6.5.2.target` |
| Fixed-temperature provenance | `chapter-6.e6.5.2.fixed-temperature` |
| Conservative-transport provenance | `chapter-6.e6.5.2.graetz-transport` |
| Fixed Dirichlet region | `dirichlet_boundary`, boundary ID `0` |
| Control region | `control_boundary`, boundary ID `1` |
| Outflow region | `outflow_boundary`, boundary ID `2` |
| Observed material ID | `1`; unobserved cells use material ID `0` |

The source-sized framework mesh provenance is
`chapter-6.e6.5.2.framework-native-rectangle-r7`; a development refinement
uses the same identifier with its realized refinement in place of `7`.
The fixed temperature is `1`, and the manufactured transport field is
$b(x)=(1.5x_{2}(1-x_{2}),0)$.

The boundary IDs realize the source geometry as follows:

| Region | Geometry | State condition |
| --- | --- | --- |
| `dirichlet_boundary` (`0`) | Left edge and upstream top/bottom walls, $0\leq x_{1}\leq 1$ | Fixed temperature $y=1$ |
| `control_boundary` (`1`) | Downstream top/bottom walls, $1\leq x_{1}\leq 4$ | Facewise control $u$ |
| `outflow_boundary` (`2`) | Right edge, $x_{1}=4$ | Zero natural transport outflow |

The frozen boundary interpretation is the ordinary-normal condition

$$
\partial_n y-(b\mathbin\cdot n)y=u
$$

on the control boundary, with zero outflow condition. The source record does
not resolve ordinary versus diffusion-weighted normal derivatives, so this is
an explicit framework choice rather than a claim of source parity. The
compiler may retain the total conservative-transport conormal as a diagnostic
alternative, but it is not the frozen benchmark interpretation.

### Solver policy

The frozen full-BFGS policy is initial step `1`, Armijo fraction `1e-4`,
backtracking factor `0.5`, at most `20` line-search trials, gradient tolerance
`1e-8`, and no declared minimum step. The source step policy is unspecified;
the artifact must retain that distinction.

### Acceptance evidence

Each of the four artifacts must retain:

- state/control dimensions and the complete compilation manifest;
- transport, fixed/controlled/outflow boundary, and material-region records;
- residual JVP/VJP and reduced-Taylor evidence for facewise control;
- objective and relative-gradient reduction; and
- state/control evidence showing the effect of changing observation region and
  target.

If the unstabilized Galerkin realization is inadequate, the limitation is
reported without changing the frozen scenario.

Execution and acceptance status are maintained in the [application roadmap](../planning/application-roadmap.md#b2),
not duplicated in this freeze.

## Future benchmark contracts

B3–B6 are part of the catalogue but are not executable freezes yet. Their
selection, dependencies, and intended evidence are defined by the
[benchmark-suite roadmap](../planning/chapter-6-benchmark-suite-roadmap.md#benchmark-sequence).
Before activating any one of them, add its source reference, application
scenario, explicit project choices, run matrix, and acceptance evidence here.

The excluded E6.5.3 reduced Stokes and E6.9.3 box-constrained Stokes examples,
stabilization comparisons, automatic OtD derivation, measure-valued state
constraints, and continuous-control boxes remain source references rather than
acceptance targets.
