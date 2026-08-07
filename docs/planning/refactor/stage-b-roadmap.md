# Stage B refactor roadmap

## Purpose and status boundary

This is the concise execution guide derived from the exhaustive
[assessment](assessment.md). It defines stable batch boundaries and acceptance
criteria. The [implementation roadmap](../implementation-roadmap.md), not this
document, owns current progress and the next-agent handoff.

The final assessment decision is:

> A small cross-cutting refactor is required before further features; focused
> refactors are additionally required only before selected features exercise
> them. No broad rewrite is warranted.

## Global execution rules

- Continue on `codex/refactor-ch5-ch6-readiness` and compare behavior with
  `pre-refactor-ch5-ch6` at `7c2496b`.
- Execute one coherent batch or named sub-batch at a time.
- Add focused characterization before or with behavior-changing repairs.
- Run focused checks first, then the applicable full Debug configurations.
- Keep feature implementation separate from prerequisite refactor commits.
- Stop at conditional C1, C2, or S1 gates until the user selects the affected
  Chapter 5/6 vertical slice.
- Preserve composition of residual, objective, metric, constraint,
  formulation, and discretization components.

The common dependency order is:

```text
R0 documentation truth
  -> R1 characterization
  -> R2a contract safety
  -> R2b semantic safety
  -> R2c current compiler/constraint correctness
  -> R3 build and test feedback
  -> common stabilization checkpoint

checkpoint -> C1 -> C2 -> selected compiler-facing feature
checkpoint -> S1       -> selected P6.1 feature
```

R2 is split for context control. The sub-batches may remain separate commits.

## R0 — Restore one trustworthy implementation narrative

**Findings:** RF-020.

**Actions:** Normalize roadmap completion markers and handoff; make the v1
compiler record the exact capability/exclusion owner; remove mutable release
status from the readiness review; distinguish current dispatch from intended
component lowering in the blueprint.

**Exit:** Status-bearing documents agree with code and tests, and no second
capability ledger or governance layer is introduced.

## R1 — Characterize refactor boundaries

**Findings:** RF-006, RF-009, RF-016.

**Actions:** Add exact diagnostic matching; make each logical scenario
independently selectable and visible to CTest; add desired-failure regressions
with their repairs; add one small assembly oracle before C2; retain and
accurately label the v0/v1 wiring comparison.

**Exit:** Every current logical scenario can run by name, exact diagnostic
identity is tested, and each accepted repair has a focused regression.

## R2a — Restore block/layout invariants

**Findings:** RF-001 and the relevant RF-006 characterization.

**Actions:** Prevent `BlockValuesT` storage dimension from diverging from its
layout and expose only checked mutation required by current callers.

**Exit:** Dimension-changing mutation is rejected, valid current algebra and
pairing identities are unchanged, and no algebra-layer redesign is added.

## R2b — Make semantic construction and validation safe

**Findings:** RF-002 through RF-005 and relevant RF-006 cases.

**Actions:** Add explicit safe incomplete enum states; close whole-graph term,
edge, binding, label, formulation, and pairing invariants; replace positional
reference-graph mutation with small ID-based helpers.

**Exit:** All recorded malformed graphs fail with exact diagnostics, every
current reference graph remains valid, and no graph DSL or staged builder is
introduced.

## R2c — Repair present provenance and projection correctness

**Findings:** the current factual defect in RF-008 and RF-012.

**Actions:** Derive the current manifest constraint realization from the
selected realization; replace caller-controlled metric display strings as
proof of projection compatibility with an opaque witness or coupled service.

**Exit:** Every current target reports the correct realization, and a
non-diagonal metric cannot obtain clipping projection by adopting a known
display label.

## R3 — Make build and test feedback explicit

**Findings:** RF-016 through RF-019.

**Actions:** Register logical CTest cases with labels and timeouts; make a
requested unavailable deal.II backend a clear configuration error; provide
explicit Debug and Release CMake presets under `build/`; integrate project
warnings into Debug profiles; keep backend-neutral sanitizers separate; add
checked size conversions only at touched adapter seams.

The intended persistent build layout is:

```text
build/
  debug-neutral/
  debug-dealii/
  sanitize-neutral/
  release-dealii/
```

`debug-neutral` is the fast loop, `debug-dealii` is the full correctness gate,
`sanitize-neutral` owns instrumented checks, and `release-dealii` owns
optimized tests and Chapter 6 timings. Warnings belong to the Debug profiles;
they do not require another permanent build tree.

**Exit:** Clean configurations intentionally support backend-neutral and
deal.II modes, requested dependencies cannot disappear silently, Debug checks
are reproducible, and Release is explicit before numerical timing.

R0 through R3 form the common stabilization checkpoint.

## C1 — Make compiler inputs, products, and provenance explicit

**Findings:** RF-008 through RF-013. Conditional on a compiler-facing feature
or comparative numerical experiment.

**Actions:** Define the diagnostic/exception boundary for bindings; establish
an explicit discretization/session lifetime; share typed state/adjoint solve
policy and convergence reporting; populate structured manifests from resolved
compiler decisions; carry the R2c realization witness through compilation.

**Exit:** Invalid inputs use the documented channel, detached compiled services
have tested lifetime, solve policy is configurable and recorded, and exact
structured manifests identify every current target.

## C2 — Introduce bounded component lowering

**Findings:** RF-007 and RF-009. Required before selected component-heavy
P5.1/P5.2/P5.4 work or P6.2 advection/stabilization.

**Actions:** Produce one ID-resolved semantic view and a modest scalar lowering
plan; let registered handlers contribute operators, data, objective pieces,
services, and provenance; extract shared FE and solve machinery only where
current targets prove identical policy; retain reconstruction, trace,
coefficient, nullspace, and gauge strategies as explicit specializations.

**Exit:** One selected scalar recombination is expressed without a new complete
problem class, current mathematical identities and independent oracles pass,
and build cost is remeasured before further template machinery is considered.

## S1 — Prepare selected P6.1 reduced methods

**Findings:** RF-014 and RF-015. Independent of C2.

**Actions:** Split value/state evaluation from derivative/adjoint augmentation;
make trial evaluations request only needed information; add only selected
direction and line-search policies; produce typed terminal and per-iteration
reports; pair reports with compiler provenance in an outer experiment envelope.

**Exit:** Rejected Armijo trials perform no adjoint solve, accepted steps are
independently auditable, exact work counts are available, and existing
covector, metric, constraint, and projected-displacement semantics remain.

## Feature gates

| Selected capability | Required gate |
| --- | --- |
| P6.1 reduced directions or line searches | Common checkpoint and S1; C1 before comparative compiled runs |
| P5.1 scalar terms, P5.2 observations/metrics, or P5.4 reconstruction | Common checkpoint, C1, and C2 |
| P5.3 low-regularity targets | Common checkpoint and C1/C2, then a selected formulation-policy design |
| P5.5/P6.5 state constraints and PDAS | Common checkpoint, then explicit KKT and complementarity contracts |
| P5.6 or Stokes variants | Common checkpoint, then an explicit mixed-block/formulation design |
| P6.2 DTO/OTD/stabilization | Common checkpoint, C1/C2, and selected transport terms |
| P6.3/P6.4 scalar KKT/preconditioning | Common checkpoint, then a bounded KKT product design |

## Explicit non-goals

- No rewrite of typed algebra, the reduced DTO, or compiled ports.
- No PDE-family inheritance, runtime plugin system, universal graph DSL, or
  general weak-form engine.
- No speculative mixed, transposition, KKT, complementarity, or PDAS contracts.
- No implementation of all Chapter 5/6 catalogue entries to prove generality.
- No distributed-memory, packaging/ABI, test-framework, documentation-renderer,
  coverage, formatter, or repository-wide native-size initiative without a
  selected requirement.

## Verification at every accepted checkpoint

- The diff contains only the approved batch and no unselected feature.
- Every repaired defect has a focused regression.
- Backend-neutral and deal.II logical cases pass in their intended profiles.
- Affected value, JVP, VJP, reconstruction, nullspace, metric, constraint, and
  projected-displacement identities retain established tolerances.
- Ownership changes receive focused sanitizer coverage.
- Numerical/compiler changes are compared with the tagged baseline and the
  independent oracle where applicable.
- Authoritative documentation and the implementation-roadmap handoff agree.
- `git diff --check` and local Markdown-link validation pass.
