# Chapter 6 remediation findings

## Purpose and status authority

This file records issues found while rechecking the Chapter 6 remediation
against the acceptance checklists. The historical review files retain the
findings and status at their original review heads; they are not rewritten as
if the historical implementation had already passed. The acceptance
checklists are updated only when current code and verification establish an
item as complete.

The [implementation roadmap](../../implementation-roadmap.md) remains the
mutable project status ledger. This file is the separate working record for
remaining issues, verification notes, and follow-up decisions.

## Status conventions

- **Open:** evidence is missing, contradictory, or exposes a remaining defect.
- **Deferred:** intentionally outside the selected scope or conditional on a
  benchmark trigger.
- **Closed:** the corresponding acceptance item is checked in its review
  document and the evidence is recorded here.

## Verified closure

### S1 preparation

The S1 checklist in
[s1-preparation-remediation-review.md](s1-preparation-remediation-review.md)
is complete. The remediation work units are `4765cfb`, `7245c31`, `150f532`,
and `1aaefbe`. Current verification passes:

- `debug-neutral`: 43/43;
- `debug-dealii`: 81/81; and
- `sanitize-neutral`: 43/43.

The focused scenarios cover staged value/derivative evaluation, incompatible
value-record rejection, detached reduced-service lifetime, policy-specific
trial work, accepted-iteration audit records, projection compatibility, and
the compilation/provenance envelope.

## Open and deferred findings

| ID | Scope | Status | Finding and next action |
| --- | --- | --- | --- |
| `CH6-F1` | Benchmark acceptance | Open | The B0 harness and selected B1/B2 executable benchmark scenarios are not yet present. Reconcile this separately from implementation-remediation closure using the [benchmark roadmap](../../chapter-6-benchmark-suite-roadmap.md). |
| `CH6-F2` | P6.4 preconditioning | Deferred | P6.4 remains conditional and is not a Chapter 6 remediation prerequisite. Activate it only if a selected all-at-once benchmark demonstrates that direct or basic serial solves are inadequate. |
| `CH6-F3` | P6.1 diagnostic evidence | Closed | Added the focused `nmopt.reduced.exact_quadratic_line_search_nonpositive_curvature_diagnostic` scenario. It constructs a negative-curvature reduced Hessian and asserts the exact `Exact quadratic line search requires positive curvature` diagnostic; P6.1 checklist item 10 is now checked. |
| `CH6-F4` | Supplemental release verification | Open, non-blocking | The attempted `release-dealii` build was interrupted after substantial progress without reaching CTest. The required P6.1 Debug and sanitizer gates pass; independently rerun the optimized profile before relying on the roadmap's release count. |
| `CH6-F5` | P6.3 compiler boundary | Closed | `CompilationProduct::quadratic_kkt` now validates and returns the canonical supplied-OTD KKT product through the same distinct owner-bearing compiler result boundary as DTO. The KKT manifest preserves supplied-OTD provenance, typed layouts/pairings, conversion, assumptions, and solver policy; unsupported targets retain the stable formulation diagnostic. |

## P6.1 verification state

All P6.1 checklist items are checked in
[p6.1-implementation-review.md](p6.1-implementation-review.md). The P6.1
remediation commits are `8865591`, `9c9de52`, `b80463a`, `5494e49`,
`b92c8d5`, `726edeb`, `8cd4fb4`, and `bdfaa9d`. Current completed profile
evidence is 44/44 `debug-neutral`; the previously recorded
`debug-dealii` and `sanitize-neutral` profiles pass 81/81 and 43/43
scenarios. The focused current execution
`nmopt.reduced.exact_quadratic_line_search_nonpositive_curvature_diagnostic`
passes and asserts the stable positive-curvature diagnostic directly.

The P6.2 remediation checklist was assessed independently of the S1 and P6.1
closures and was not inferred solely from the roadmap's aggregate Chapter 6
closure statement.

## P6.2 verification state

All P6.2 checklist items are checked in
[p6.2-implementation-review.md](p6.2-implementation-review.md). The
remediation work units are `53bb1c6`, `28ef9d7`, `b4ceda5`, and `2437449`.

The semantic contract validates typed state, adjoint, and stationarity
declarations and rejects DTO-label mutation, missing pairings, and missing
multiplier conversion. The compiler's separate supplied-OTD product retains
the declaration and validates its runtime layouts, block dimensions, action
provenance, and comparison status before constructing the manifest. Focused
diagnostics cover mismatched adjoint spaces, pairings, block signs, and
multiplier conversion.

The owned-session scenarios exercise callbacks, solve, and teardown after the
mesh/session scope has ended. The native deal.II scenario checks every
supplied residual block, centered JVP finite differences, full JVP/VJP pairing,
and state/adjoint/stationarity equality with the DTO action under the declared
identity conversion. The dense reference scenarios independently check typed
blocks, JVP finite differences, VJP pairing, and the reduced-DTO solution
comparison.

Focused current executions passed:

- `nmopt.semantic.v1_supplied_otd_declaration`;
- `nmopt.supplied_otd.compiled_owned_product_lifetime`;
- `nmopt.supplied_otd_reference` (both scenarios);
- `nmopt.dealii.canonical_volume_control`; and
- `nmopt.dealii.supplied_otd_owned_session_lifetime`.

The required current profile counts are 43/43 `debug-neutral`, 81/81
`debug-dealii`, and 43/43 `sanitize-neutral`. The P6.3 review is recorded
below; the remaining open/deferred findings `CH6-F1`, `CH6-F2`, `CH6-F3`, and
`CH6-F4` remain independent of this P6.2 closure.

## P6.3 verification state

All P6.3 checklist items are checked in
[p6.3-implementation-review.md](p6.3-implementation-review.md). Items 1–2
are now closed: the remediation registers both the distinct compiled DTO KKT
product and the canonical supplied-OTD KKT product through
`CompilationProduct::quadratic_kkt`. The supplied-OTD KKT bridge remains a
typed adapter owned by the returned compiler product, while the manifest
retains its distinct formulation provenance and KKT construction record.

The P6.3 remediation commits are `06fa551`, `5aa196c`, `3d20acb`, `01dd59a`,
`3968414`, `4415044`, and `559737f`. Focused current executions passed:

- all four backend-neutral quadratic-KKT scenarios;
- all three dense DTO/supplied-OTD reference scenarios;
- the supplied-OTD bridge ownership scenario;
- all four deal.II KKT scenarios; and
- the compiled DTO and supplied-OTD KKT compiler routes and detached
  supplied-session scenarios.

The full `debug-neutral` and `debug-dealii` profiles passed 43/43 and 81/81.
The previously recorded sanitizer gate passed 43/43; a fresh default CTest
invocation is blocked by this runner's LeakSanitizer/ptrace restriction, while
the same profile with leak detection disabled passed. This is an environment
limitation, not a new code finding. The next review unit is P6.5, but P6.3's
compiler-boundary finding should remain open until its scope is resolved.

## P6.5 verification state

The current P6.5 acceptance checklist is recorded in
[p6.5-implementation-review.md](p6.5-implementation-review.md). Its
remediation range is `9579784^..58fc8f6`.

The P6.5-specific R1–R7 outcomes are now evidenced by the current code and
focused scenarios:

- restricted free-coordinate active KKT products and full-primal
  reconstruction, including affine shifts and the empty-set identity;
- mixed base/active multiplier action and transpose pairing;
- compiled DTO and supplied-OTD PDAS routes with structured manifest records
  and stable negative diagnostics;
- owned metric, active-product, solver, and compiled-product lifetimes;
- rejection of non-finite inputs, conversions, solver outputs, and residuals;
- independent inactive-box agreement and active native diagnostics; and
- detached execution through the neutral and compiled deal.II paths.

Focused current executions passed all ten backend-neutral PDAS/active-set
scenarios, both native deal.II PDAS scenarios, and the compiled deal.II
`compiled_pdas` scenario. The current profile evidence is 43/43
`debug-neutral`, 81/81 `debug-dealii`, and 43/43 `sanitize-neutral`; the fresh
default sanitizer invocation has the runner's LeakSanitizer/ptrace limitation
and passes with leak detection disabled.

No new P6.5 implementation finding was identified. `CH6-F4` remains a
separate, non-blocking optimized-profile verification item.
