# S1 preparation review and remediation handoff

## Status and authority

This document records the current-state review of the S1 preparation gate at
`5a4d83c` on 2026-08-16. S1 has no dedicated implementation range: no commit
message names S1, `RF-014`, or `RF-015`. Parts of the gate arrived indirectly
with P6.1, but the gate was never separately closed.

This is static review evidence and an implementation handoff, not a second
status ledger. The [implementation roadmap](../../implementation-roadmap.md)
remains the sole owner of mutable status and acceptance state. The
[Stage B S1 gate](../pre-ch5-ch6/stage-b-roadmap.md#s1--prepare-selected-p61-reduced-methods),
[refactor assessment](../pre-ch5-ch6/assessment.md#s1--prepare-only-the-selected-p61-reduced-methods),
[Chapter 6 numerical-methods guide](../../../guides/chapter-6-numerical-methods.md#c63--reduced-space-methods-for-unconstrained-ocps),
and [required verification and provenance](../../../implementation/implementation-readiness-review.md#12-required-verification-and-provenance)
remain authoritative.

S1 is partially implemented. The selected direction and line-search policy
concepts exist, terminal enums and aggregate histories exist, and current
results retain a final full reduced evaluation. The two findings that created
the gate remain open:

1. `RF-014` – reduced evaluation is still monolithic, so every objective-only
   trial performs an adjoint solve; and
2. `RF-015` – solver results still lack a typed accepted-iteration audit and
   an outer run envelope that binds the report to compilation and environment
   provenance.

S1 must close before P6.1 defect remediation. It is a separate gate: any S1
evidence that also bears on a finding in the
[P6.1 implementation review](p6.1-implementation-review.md) must be assessed
and recorded there rather than assumed to close it. S1 does not activate a
benchmark.

## Current implementation assessment

| S1 obligation | Current evidence | Status |
| --- | --- | --- |
| Split value/state evaluation from derivative/adjoint augmentation | `ReducedDTOT::evaluate()` always performs both solves and returns one `ReducedEvaluationT` | missing |
| Request only the information used by a trial | the line-search evaluator and trust-region trial path both require a full `ReducedEvaluationT` | missing |
| Add only selected direction and line-search policies | steepest descent, nonlinear CG, BFGS/L-BFGS, Newton, Armijo, exact-quadratic, and Wolfe policies are present within the selected P6.1 scope | complete for S1; retain |
| Produce typed terminal and per-accepted-iteration reports | typed terminal reasons, final evaluation, parallel histories, and aggregate work counts exist; one typed accepted-iteration record does not | partial |
| Pair reports with compiler and run provenance | `CompilationManifest` exists on the compiled product, but no experiment envelope owns it together with a solver policy snapshot, report, and environment | missing |

The decisive exit criterion currently fails. `ReducedSearchSolverT` calls
`ReducedDTOT::evaluate()` for every trial and increments both state and adjoint
counts. The dense backtracking test explicitly expects
`state_solve_count == adjoint_solve_count`. The trust-region acceptance path
also computes a full derivative for rejected objective-ratio trials.

The reporting surface has advanced beyond the original `RF-015` evidence, but
not enough to close it. `ReducedSolverResultT` now contains final accepted
state/adjoint/covector data, objective and stopping histories, step histories,
Hessian diagnostics, and aggregate state, adjoint, metric, and Hessian work.
Those parallel arrays do not retain the descent pairing or policy-specific
acceptance evidence needed to audit one accepted iteration independently, and
they do not identify the solver parameters or compiled product that produced
the report.

No build or test command was run for this static assessment. The assessment
uses the current code, committed tests, and authoritative S1 documents.

## Behavior to retain

The remediation must preserve the following selected behavior:

- reduced derivatives remain covectors and metric gradients and search
  directions remain primal;
- the reduced covector keeps the sign
  $`J_{u}'-E_{u}'^{\ast}p`$;
- state and adjoint solves retain their typed `LinearSolveReport` records and
  owned compiled-service lifetime;
- projected acceptance and stopping continue to use the actual projected
  displacement and the declared metric/constraint coupling;
- Wolfe trials continue to request a trial derivative because their curvature
  condition uses it;
- exact-quadratic and second-order policies continue to require an explicit
  reduced-Hessian capability;
- existing direction restart, metric-action, Hessian-action, and inner-solve
  diagnostics remain visible;
- trust-region trial diagnostics remain distinct from line-search accepted-
  iteration reporting; and
- `ReducedDTOT::evaluate()` may remain as a compatibility composition of the
  new value and derivative operations.

## Required reading for the implementing agent

Before changing code, follow the repository routing in
[`conventions/README.md`](../../../../conventions/README.md). For one bounded
S1 unit, read only the relevant authorities and current interfaces:

- [Stage B S1 actions and exit](../pre-ch5-ch6/stage-b-roadmap.md#s1--prepare-selected-p61-reduced-methods);
- [assessment RF-014](../pre-ch5-ch6/assessment.md#rf-014--the-reduced-evaluation-protocol-performs-an-adjoint-solve-for-every-objective-only-trial);
- [assessment RF-015](../pre-ch5-ch6/assessment.md#rf-015--solver-reporting-cannot-yet-audit-a-chapter-6-run-or-carry-its-compilation-provenance);
- [selected reduced methods](../../../guides/chapter-6-numerical-methods.md#c63--reduced-space-methods-for-unconstrained-ocps);
- [numerical-example provenance requirements](../../../guides/chapter-6-numerical-examples.md#purpose-and-use);
- [P6.1 retained behavior and open findings](p6.1-implementation-review.md);
- [`ReducedDTOT`](../../../../include/nmopt/contract/reduced_dto.hpp);
- [shared reduced result and direction contracts](../../../../include/nmopt/solvers/reduced_search.hpp);
- [line-search policies](../../../../include/nmopt/solvers/reduced_line_search.hpp);
- [line-search solver integration](../../../../include/nmopt/solvers/reduced_gradient.hpp); and
- [trust-region solver integration](../../../../include/nmopt/solvers/reduced_trust_region.hpp).

Do not turn S1 into a generic evaluation query language, universal solver
registry, serialization framework, or benchmark harness. Do not add new
directions, line searches, projected advanced methods, generic nonlinear
Hessians, or PDE-family branches. The backend-neutral solver must not depend
on the deal.II compiler or its manifest type.

## Required end state

### Split reduced evaluation

The reduced DTO boundary must expose two typed stages:

1. a value/state evaluation that validates the control, performs exactly one
   state solve, and returns the state, full point, objective value, state-solve
   report, and the validity/layout information needed for safe reuse; and
2. derivative/adjoint augmentation of a valid value record that performs
   exactly one adjoint solve and returns the adjoint, reduced covector, and
   adjoint-solve report without repeating the state solve or objective value.

The value record must not permit accidental augmentation after incompatible
control, layout, model, or lifetime changes. `evaluate(control)` may compose
the two operations and must remain numerically equivalent to the current full
evaluation during migration.

Invalid layouts, missing callbacks, incompatible returned blocks, and violated
metric or constraint contracts remain exceptions. A declared linear solve
that does not converge is an expected algorithm event only where the solver
can report it meaningfully; it must not be silently converted into a valid
evaluation.

### Request policy-specific trial information

Each selected policy must request only the information used by its acceptance
predicate:

| Policy or phase | Trial information | Derivative work after acceptance |
| --- | --- | --- |
| initial accepted point | full value plus derivative | already available |
| Armijo | value/state only | augment the accepted value record once |
| exact quadratic | value/state for the returned trial; Hessian action remains separate | augment the accepted value record once |
| weak or strong Wolfe | value/state plus derivative for every tested slope | reuse the accepted full record |
| trust-region reduction-ratio trial | value/state only | augment only an accepted value record |

The line-search boundary may use separate value and derivative evaluators, a
small typed information request, or policy-specific overloads. Whichever
shape is selected must make the work requirement visible in the type-level
interface and must not branch on algorithm names inside `ReducedDTOT`.

An accepted value record must be cached and augmented rather than recomputed.
Objective- or step-based termination immediately after acceptance must still
return the accepted state, adjoint, and reduced covector required by the final
report. S1 changes information flow and accounting; it must not conceal or
implicitly repair `P6.1-R5` exact-search semantics or the other open P6.1
correctness findings.

### Produce auditable typed reports

Replace parallel-history-only reporting with one typed record per accepted
iteration. Each record must contain enough typed or scalar evidence to
recompute the declared acceptance and stopping decisions:

- accepted iteration index and policy identity;
- objective before and after the step and the objective change;
- requested step parameter, actual displacement norm, and descent pairing
  formed with that actual displacement;
- policy-specific acceptance evidence, such as the Armijo bound, Wolfe trial
  slope, exact-quadratic curvature, or trust-region reduction ratio;
- absolute and relative stationarity measures used by stopping;
- trial count for the accepted step;
- per-iteration and cumulative state, adjoint, metric, and Hessian work;
- inner Hessian-solve diagnostics and direction-reset information where the
  selected policy produces them; and
- the accepted evaluation or an explicit association with the accepted
  control/state/adjoint/covector record.

The terminal report must retain the final accepted evaluation, a snapshot of
the selected solver and policy parameters, the typed stopping outcome, and
the complete iteration records. Existing history vectors may remain as
compatibility views only if they are derived from or checked against the
typed records rather than maintained as an independent source of truth.

Expected algorithm outcomes—iteration limit, exhausted line search,
non-finite trial, policy-declared curvature/reset outcome, radius exhaustion,
or failed inner solve when recovery or comparison is meaningful—must use
typed report status. Invalid wiring and broken algebraic contracts remain
exceptions. Do not relabel the open P6.1 zero-step, exact-search, boundary, or
stationarity defects as successful terminal outcomes.

### Add the outer experiment envelope

Add one narrow run-level type outside the backend-neutral solver. It must own
or value-copy:

- a schema version;
- the C1 `CompilationManifest` or an explicitly typed manifest parameter;
- the formulation and compiled-product identity already carried by that
  manifest;
- the solver/direction/globalization parameter snapshot;
- the typed solver report; and
- run-environment metadata sufficient for the selected numerical-example
  contract, including build profile and compiler identity, relevant dependency
  versions, process/thread policy, hardware description, and the wall-clock
  definition when timings are recorded.

The envelope proves association; S1 does not require a stable disk format,
database, CLI runner, or execution of B1/B2. The generic solver returns its
backend-neutral report without importing `CompilationManifest`. A compiler or
future experiment layer constructs the envelope after it owns both products.

## Ordered implementation plan

### Work unit 1 – split value and derivative evaluation

**Outcome:** one reusable value/state record and one checked
derivative/adjoint augmentation operation; the existing full evaluation is a
compatibility composition.

**Boundary:** `RF-014` at the reduced DTO contract. Do not change direction,
line-search, stopping, or benchmark behavior in this unit.

**Likely files:**

- `include/nmopt/contract/reduced_dto.hpp`;
- `tests/reduced_dto_contract.cc`;
- `docs/guides/chapter-6-numerical-methods.md`; and
- any concrete reduced-service adapter that must expose the staged operation.

**Focused verification:** one state and zero adjoint callbacks for a value
request; one later adjoint and zero repeated state callbacks for augmentation;
full-evaluation equivalence; state/adjoint failure behavior; and stale,
foreign, or incompatible value-record rejection before an adjoint callback.

**Prospective commit:**
`refactor(dto): split reduced value and derivative evaluation`

### Work unit 2 – request only policy-required trial work

**Outcome:** Armijo, exact-quadratic, and trust-region objective trials use
value records, accepted trials are augmented once, and Wolfe retains the full
trial derivative it needs.

**Boundary:** `RF-014` solver integration. Retain selected policy formulas,
actual-displacement semantics, and every P6.1 exclusion. Do not repair the
separate P6.1 review findings in this unit.

**Likely files:**

- `include/nmopt/solvers/reduced_line_search.hpp`;
- `include/nmopt/solvers/reduced_gradient.hpp`;
- `include/nmopt/solvers/reduced_trust_region.hpp`;
- `include/nmopt/solvers/reduced_search.hpp`; and
- `tests/reduced_dto_contract.cc`.

**Focused verification:** callback-level work oracles for accepted and failed
backtracking, Wolfe curvature trials, exact-quadratic evaluation, and accepted
and rejected trust-region trials. Include a projected Armijo case to prove
that the actual displacement and metric/constraint semantics did not change.

**Prospective commit:**
`refactor(solver): stage reduced trial evaluations`

### Work unit 3 – introduce accepted-iteration audit records

**Outcome:** line-search and trust-region results have typed terminal status,
policy snapshots, exact work records, and independently auditable accepted
iterations instead of unrelated parallel histories.

**Boundary:** the reporting part of `RF-015`. Share common report components
where their semantics are identical, but retain policy-specific line-search
and trust-region diagnostics. Do not force every solver into one lossy result
variant.

**Likely files:**

- `include/nmopt/solvers/reduced_search.hpp`;
- `include/nmopt/solvers/reduced_line_search.hpp`;
- `include/nmopt/solvers/reduced_gradient.hpp`;
- `include/nmopt/solvers/reduced_trust_region.hpp`;
- `tests/reduced_dto_contract.cc`; and
- `docs/guides/chapter-6-numerical-methods.md`.

**Focused verification:** reconstruct every accepted Armijo/Wolfe predicate
from its record; match delta and cumulative work to instrumented callbacks;
check policy snapshots and history projections; and distinguish expected
terminal outcomes from invalid-layout or invalid-metric exceptions.

**Prospective commit:**
`feat(report): add reduced iteration audit records`

### Work unit 4 – bind reports to compilation provenance

**Outcome:** a typed experiment envelope owns one structured compilation
manifest, one solver-policy snapshot, one typed report, and one run-environment
record without adding compiler dependencies to generic solvers.

**Boundary:** the outer-envelope part of `RF-015`. Define only the in-memory
association needed before comparative compiled runs; leave serialization,
scenario execution, timing collection, and B0–B6 orchestration to the
benchmark roadmap.

**Likely files:**

- one new narrow experiment/reporting contract header, with its layer and
  ownership documented in the same unit;
- `include/nmopt/compiler/v1/compiled_problem.hpp` only if a non-deal.II
  manifest-facing seam is required;
- `tests/reduced_dto_contract.cc` and/or
  `tests/dealii_diffusion_contract.cc`;
- `docs/implementation/implementation-readiness-review.md`;
- `docs/guides/chapter-6-numerical-examples.md`; and
- `docs/planning/implementation-roadmap.md` after all S1 checks pass.

**Focused verification:** construct an envelope from a compiled scalar target
and its detached reduced service; prove that the manifest survives service
detachment and distinguishes materially different compiled products; verify
that policy and environment changes alter their own typed records; and prove
that the solver headers remain independent of compiler/deal.II types.

**Prospective commit:**
`feat(experiment): bind solver reports to compilation provenance`

## Verification gate

Run the focused backend-neutral S1 scenarios after each unit. At the completed
gate, run:

```bash
cmake --preset debug-neutral
cmake --build --preset debug-neutral
ctest --preset debug-neutral

cmake --preset debug-dealii
cmake --build --preset debug-dealii --parallel 1
ctest --preset debug-dealii

cmake --preset sanitize-neutral
cmake --build --preset sanitize-neutral
ctest --preset sanitize-neutral
```

S1 does not require a timing claim. Run `release-dealii` when the remediation
sequence reaches B1/B2 reproduction or if a later change reports optimized
numerical behavior.

Before requesting each commit, inspect the staged diff and statistics,
confirm that unrelated P6.1 remediation remains unstaged, and run
`git diff --check`.

## Acceptance checklist

S1 can be marked complete in the implementation roadmap only when:

- [x] the reduced DTO exposes checked value/state and derivative/adjoint
      stages;
- [x] the compatibility full evaluation produces the same mathematical
      result without repeating work;
- [x] an accepted value record cannot be augmented under an incompatible
      control, layout, model, or lifetime;
- [x] rejected Armijo trials perform state solves and no adjoint solves;
- [x] Armijo, exact-quadratic, and trust-region accepted values are augmented
      once without a repeated state solve;
- [x] Wolfe trials retain the derivative work required by their curvature
      predicate;
- [x] state, adjoint, metric, Hessian, inner-solve, and trial counts match
      independent callback-level oracles;
- [x] every accepted step has one typed record sufficient to audit its
      acceptance and stopping decisions;
- [x] the terminal report retains the final accepted state, adjoint, reduced
      covector, objective, policy snapshot, and typed outcome;
- [x] expected algorithm outcomes and invalid contract wiring use their
      distinct documented channels;
- [x] one outer envelope pairs the structured C1 manifest, policy snapshot,
      report, and environment without coupling the generic solver to the
      compiler;
- [x] projected-displacement, covector, metric, constraint, Hessian-capability,
      and selected-policy semantics remain intact;
- [x] all focused and required neutral, deal.II Debug, and sanitizer profiles
      pass; and
- [x] the implementation roadmap alone is updated to mark S1 complete and to
      hand off to P6.1 remediation.

## Expected state after S1

For an Armijo or objective-ratio run starting from one fully evaluated point,
with `T` total trials and `I` accepted iterations, the expected first-order
work identity is:

```text
state_solve_count   = 1 + T
adjoint_solve_count = 1 + I
rejected_trial_count = T - I
```

The same identity applies when the last line search fails: its rejected trials
increase the state count but not the adjoint count. Wolfe is intentionally
different because every curvature trial requests a derivative:

```text
state_solve_count   = 1 + T
adjoint_solve_count = 1 + T
```

Counts for metric inverses, Hessian actions, inner iterations, and resets must
likewise match instrumented callbacks and the selected policy rather than a
hard-coded relationship to outer iterations.

At the end of S1, a caller can determine what information every trial
requested, audit every accepted step from one record, distinguish normal
algorithm termination from invalid wiring, and associate a compiled run with
the exact compilation and execution context that produced it. P6.1 remains
acceptance-pending until its own open findings are remediated, and B1/B2 remain
inactive until the later problem-library and benchmark gates pass.
