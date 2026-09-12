# External deal.II boundary evaluation: closure report

Status: closed on 2026-09-12 by user decision. Problem A, Problem B, the
minimal consumers, and their corrective reviews are complete. No further
implementation unit is selected.

Reviewed implementation: `2ba749b`, including the MC4–MC6 corrections
`edf16f2`, `f5f2fe5`, and `2ba749b`.
Baseline refactor: `67104fd376a9b1d25d7e6cdaca664a7fe1aa98ed`.
Phase branch: `codex/evaluate/external-dealii-boundary`.
Authority: [evaluation roadmap](../../external-dealii-boundary-evaluation.md)
and [Problem B protocol](problem-b-protocol.md).

## Question and bounded result

The original concern was whether nmopt's growing implementation could still
be used as a library by an existing deal.II application with modest adaptation
and understandable wiring. The failed initial tutorial attempt mixed together
PDE reuse, new optimal-control mathematics, framework construction, and
verification. The controlled evaluation separated those responsibilities
before considering any shared interface change.

The current public boundary successfully supports the two tested Step-4
problems while preserving application ownership and numerical policies.
The compiler, recipes, manifests, and project runner are unnecessary for the
external path. The minimal consumers expose the remaining construction work
directly. That work is still verbose and requires several public-contract
concepts; effortless integration has not been demonstrated.

Closing this phase records an answered investigation, not a claim that every
original usability aspiration has been achieved.

## Objectives and disposition

| Objective | Disposition | Evidence and practical meaning |
| --- | --- | --- |
| Reuse authentic upstream Step-4 with small native seams | Reached | Exact upstream and stripped fixtures are retained; forward fidelity checks cover both original dimensions. Adaptations expose numerical reuse and B's native discretization views. |
| Separate intrinsic mathematics from framework cost | Reached for A/B | Application reuse, OCP operations, native orchestration, binding, and verification have distinct owners and recorded work. |
| Compare identical mathematics and selected policies | Reached | Native and nmopt paths agree on reduced evaluations and matched optimization, with independent final audits. |
| Exercise actual FE control geometry | Reached for B | Free state coordinates, fixed lifting, rectangular mass coupling, and a full control mass metric remain application-owned. |
| Keep shared nmopt and the refactor baseline intact | Reached | No shared API, compiler, backend/storage, or optimizer change was required. |
| Exhibit consumers without evaluation machinery | Reached | Separate runnable A/B bindings use the existing OCP operations without instantiated instrumentation or native reference optimizers. |
| Establish very little boilerplate and few concepts | Partly reached | The functional path is visible, but repeated layout, block, solve-report, metric, and DTO construction remains. |
| Establish newcomer usability or general project maintainability | Not established | No unfamiliar-user exercise or repository-wide maintainability audit was performed. |

## Evidence and validation

The numerical evidence is owned by the [G1 report](g1-report.md),
[G2 report](problem-b-report.md), and
[minimal consumer assessment](minimal-consumers-report.md). Earlier passing
results remain tied to their historical revisions; later test counts do not
retroactively change those records.

| Observation | Problem A | Problem B |
| --- | ---: | ---: |
| State/control dimensions | 289 / 289 | 225 / 289 |
| Accepted iterations | 828 | 5 |
| Line-search trials | 6,025 | 5 |
| State solves per matched path | 6,026 | 6 |
| Adjoint solves per matched path | 829 | 6 |
| Final optimizer gradient norm | $9.5036543162094535\times10^{-7}$ | $4.8845615376102665\times10^{-8}$ |

Both paths used the same frozen steepest-descent/Armijo schedule within each
problem. These iteration counts compare native with nmopt for the same
problem; they are not an A-versus-B optimizer benchmark. B changes the
objective, coupling, coordinates, and geometry simultaneously.

Verification includes off-solution residual/JVP/VJP checks, reduced derivative
checks, fresh final state/adjoint/covector recomputation, independent dense
oracles, and native field comparisons. B's stationarity audit solves the mass
system independently using dense algebra, rather than reusing the runtime CG
inverse. Its audited norm is $4.8845615376171785\times10^{-8}$; the tiny
difference from the optimizer norm reflects the separate computations.

MC5 additionally launches the actual minimal executables and checks their
stopping reason, counts, objective, gradient norm, and VTK payload against
audited native results. Final implementation validation at `f5f2fe5` passed
195/195 `debug-dealii` and 67/67 `debug-neutral` tests. Review at `2ba749b`
checked the source changes, saved test logs, and generated output; the
documentation-only closure does not claim a new deal.II numerical run.

Evidence remains ignored below `runs/external-dealii/step-4/`, including the
corresponding tree under `build/debug-dealii/` when tests run from that
working directory. Reproduction commands and source-to-evidence links remain
in the three reports above and the [application README](../../../../apps/external-dealii/step-4/README.md).
Generated traces are not tracked acceptance fixtures.

## Causal findings

| Finding | Classification | What follows |
| --- | --- | --- |
| A closed `run()` application needs callable assembly/solve/output seams | Application reuse | This work exists before selecting an optimizer library. |
| A forward PDE has no chosen control, objective, adjoint, or metric | Intrinsic OCP work | The application must supply those mathematical decisions. B's mass/coupling and coordinate work cannot fairly be counted as nmopt boilerplate. |
| Both bindings repeat layouts, block conversions, two solve wrappers, report translation, and DTO construction | Observed construction burden | A bounded construction-convenience experiment could be motivated by this evidence, without first adding another PDE. None is selected here. |
| Residual and JVP are mandatory at construction but unused by this reduced optimization path | Capability obligation | A helper cannot remove the underlying requirement; reconsidering it would require a separate boundary decision. |
| Full VJP computes a state transpose contribution discarded by the reduced derivative | Measured additional work | Native uses a control pullback; nmopt performs 829 extra explicit state transpose actions in A and 6 in B. No timing claim follows. |
| B's metric adapter checks native convergence but returns a primal block | Current metric contract | Public metric inversion supplies no solve report. Evaluation instrumentation must not be presented as an API requirement. |
| PDE ownership and output stay native as B becomes more realistic | Evidence against forced cross-layer coupling in these cases | No boundary redesign or deeper structural repair is required by this experiment. |

The minimal binding counts are 206 nonblank, non-comment lines for A and 222
for B; entry points add 107 and 112. These are implementation observations,
not universal lower bounds. Local classes, error handling, output allocation,
and frozen solver settings are not all framework obligations. The
[construction comparison](minimal-consumers-report.md#a-and-b-construction-comparison)
explains the repeated concepts and B-specific responsibilities.

## Remaining limits and hypotheses

Both problems use one linear, symmetric, serial Step-4 application on a fixed
mesh. B removes A's identity coupling/metric convenience, but remains a
favorable same-space case: its exact mass cancellation does not establish
generality for arbitrary control spaces. Neither experiment establishes cost
for unrelated applications, nonlinear or nonsymmetric PDEs, MPI, adaptivity,
control constraints, package installation, or compiler-path equivalence.

The measured binding repetition is an observed limitation. Its effect on
learning time, the benefit of a helper, and the performance materiality of
extra VJP/copy work remain unmeasured. No general API-minimality theorem or
overall codebase-health verdict follows from a successful external path.

A deeper structural investigation would need evidence such as forced compiler
or PDE knowledge in generic solvers, ownership violations, or a mathematical
capability that cannot be represented without cross-layer changes. That was
not required here. Unused capability requirements remain a local boundary
question even though they did not block correctness.

## Closure decision and documentation

Retain the refactor baseline, evaluated native/nmopt paths, independent
verification, and two explicit consumers. Close the evaluation and consumer
follow-up. No helper, new public port, further PDE, benchmark replication,
profiling, or shared-library change is authorized by this closure.

The [Step-4 overview](../../../../apps/external-dealii/step-4/external-integration-overview.md)
explains the complete application path, and the
[implementation report](../../../../apps/external-dealii/step-4/integration-report.md)
accounts for the concrete source responsibilities. The
[external integration reference](../../../reference/external-dealii-solver-integration.md)
documents the implemented public API. Source-supported behavior belongs in
that reference; broader design hypotheses remain in this review record and
the historical [design investigation](design-investigation.md). This closure
introduces no new authoritative architectural contract.
