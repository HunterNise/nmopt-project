# Chapter 5 remediation audit follow-ups

## Purpose and scope

This file records issues found while auditing the Chapter 5 remediation
commits against the acceptance checklists in this directory. The audit is
bounded to deciding whether implementation may proceed to Chapter 6. It does
not reopen deferred architecture or unselected capabilities.

The [implementation roadmap](../../implementation-roadmap.md) remains the
status ledger. Checked acceptance items remain in their owning review files.

## Open issues

### CH5-AUDIT-001 – the closed compiler request is not yet the sole decision source

**Chapter 6 impact:** does not block starting P6.1. Resolve or explicitly
accept this limitation before P6.2, whose advection and stabilization work is
gated by C2.

`ResolvedCompilationRequest` centralizes target and boundary selections, and
later target dispatch consumes that snapshot. However,
`close_compilation_request()` still populates it with raw-`ProblemSpec`
predicate helpers, while `validate_lowerability()` independently consumes and
rescans the raw specification. The C2 acceptance conditions that one closed
request drives all later stages and that raw graph flags no longer duplicate
those decisions therefore remain unchecked.

The late residual-lowering follow-up does make scalar residual execution
plan-owned: `ScalarComponentModel` iterates the selected contribution records
and accepts only a closed bounded residual registration. The remaining issue
is therefore limited to upstream request construction and validation, not
residual assembly.

The current declaration-order regression passes, and this audit found no
failing registered Chapter 5 execution caused by the duplication. Closing the
issue requires only a bounded reconciliation of the remaining validation and
selection paths; it does not justify a new compiler framework.

### CH5-AUDIT-002 – the manifest is not yet a pure resolved-decision projection

**Chapter 6 impact:** does not block starting P6.1, whose first work is in the
reduced solver and Hessian-service layers. Resolve or explicitly accept this
limitation before P6.2, where the manifest must distinguish DTO from a supplied
OtD system without reconstructing formulation choices.

The version-2 manifest now records exact scalar values and bound digests,
caller provenance for opaque `Function` objects, mesh structural identity,
independent solve policies, and realized map/space records. Those changes close
the earlier non-identifying provenance cases.

However, `DealiiCompiler::make_manifest()` still receives the raw
`ProblemSpec`, `CompiledTargetKind`, and target registration. It rebuilds a
large set of `uses_*` decisions and independently derives solve, metric,
observation, lifting, handler, and assumption records or compatibility prose.
Consequently the C1 requirements that the manifest be populated only from
resolved decisions and realized services, and that all compatibility prose be
rendered from typed records, remain unchecked. The current tests compare a
useful subset of projections and differentials, but do not provide the required
exact structured table for every retained target.

This shares the decision-duplication root of `CH5-AUDIT-001`, but is recorded
separately because its closure boundary is the manifest projection rather than
validation or dispatch. A bounded fix should move the remaining typed service
facts into `ResolvedCompilationDecision` and make `make_manifest()` a renderer;
it does not require a new serialization framework.

### CH5-AUDIT-003 – boundary-trace realized-map ownership and dimension remain unproved

**Chapter 6 impact:** does not block starting P6.1. Resolve or explicitly
accept the limitation before treating C1 as closed or relying on structured
observation provenance in P6.2 and the Chapter 6 benchmarks.

The common realized-map schema records source/output spaces, dimensions,
ordering, pairing, transformation chains, and value/JVP/VJP descriptions. The
point-sensor and normal-flux models expose matching executable actions, and
their current dimension and pairing tests pass.

The records themselves are nevertheless synthesized after model construction
by `make_realized_maps()` from the raw semantic graph and target-specific
`dynamic_cast` queries. They are not published by the objects that implement
the actions. In particular, the baseline and weighted Neumann boundary traces
do not expose the recorded face-sample value/JVP/VJP map. Their
`observation_sample_count_` is incremented inside the local shape-function
loop, so it counts unconstrained basis contributions rather than the ordered
face-quadrature samples claimed by the realized-map record. The weighted-trace
test compares the record with that same counter and therefore cannot detect
the mismatch.

The C1 acceptance item requiring every realized map to own its dimensions,
ordering, and pairing remains unchecked. A bounded closure should choose the
actual fused boundary-trace output representation, count that representation
directly, and publish the record alongside equivalent value/JVP/VJP actions;
the point-sensor and normal-flux implementations need not be redesigned.

### CH5-AUDIT-004 – typed transposition capability still depends on policy prose — resolved

**Chapter 6 impact:** does not block starting P6.1. Resolve before P6.2,
whose supplied-OtD boundary depends directly on unambiguous formulation and
transposition provenance.

P5.3 and P5.4 now share typed transposition, partial-boundary, fractional-
metric, and boundary-$H^{1}$ selections. Semantic and compiler checks reject
unsupported enum selections, and the resolved decision and manifest retain
the typed payloads.

The compiler lowerability path now consumes the resolved typed transposition,
boundary, metric, and target selections; it does not search `selected_policy`
for English fragments. Semantic selection presence is likewise determined by
the selected status/scope and typed payloads. The untyped domain-regularity
assumption remains an explicit textual declaration because it has no structured
replacement in the v1 contract.

The point-sensor, normal-flux, and Section 5.11 deal.II regressions clear the
relevant policy descriptions before validation and compilation. The full
`debug-neutral` (12/12) and `debug-dealii` (41/41) baselines pass, so the P5.3
prose-independent-validation and P5.4 prose-independent-factory acceptance
items are closed. Compatibility rendering remains a manifest concern already
addressed by `CH5-AUDIT-002`.

### CH5-AUDIT-005 – P5.4 realized-layout and Chapter 5.11.2 Taylor evidence is incomplete — resolved

**Chapter 6 impact:** does not block starting P6.1. Close this evidence gap
before marking P5.4 complete or using its structured layouts in Chapter 6
benchmark reporting.

The cross-cutting realized-map commit makes
`DirichletControlLiftingModel::physical_state_dimension()` the generic source
for Dirichlet state-observation dimensions. Static inspection shows the
current path keeps the independent state, test, and control dimensions
separate while assigning the physical reconstructed dimension to state
observations and the trace-control dimension to control observations.

The existing P5.4 scenarios now share an assertion over the five canonical
manifest spaces: independent state, state test, control, physical
state-observation, and control observation. It is exercised by the partial
target, Chapter 5.11.2, and all three remaining Section 5.11 registrations;
the refined mesh also asserts that the physical state layout differs from the
independent state layout. Chapter 5.11.2 now performs the same reduced,
state-recomputed directional Taylor check used by the other registrations.
The focused P5.4 scenarios and both full debug baselines pass, closing this
acceptance-evidence gap without adding an executable observation interface.

### CH5-AUDIT-006 – the roadmap declares remediation acceptance beyond the evidence

**Chapter 6 impact:** this prevents treating the Chapter 5 remediation review
as completed, but it does not block starting the bounded P6.1 implementation.
The remaining issues must be resolved or explicitly accepted before the
roadmap is used as the gate for P6.2 and Chapter 6 benchmark provenance.

Commit `6a84b34` marks P5.2--P5.4 acceptance-complete, and `d3cd628` extends
that claim to C1/C2. The neutral, deal.II Debug, and neutral sanitizer profiles
do pass at the audited head. The acceptance claims nevertheless conflict with
unchecked items in the owning remediation reviews and with
`CH5-AUDIT-001`--`CH5-AUDIT-005` and `CH5-AUDIT-007`: the request and
manifest still duplicate raw decisions, realized boundary-trace ownership is
unresolved, transposition capability still reads prose, and required P5.4
layout/Taylor evidence is missing.

The roadmap and v1 capability record should describe those bounded open items
instead of declaring the whole remediation closed. After the fixes land, or
the limitations are explicitly accepted at their stated Chapter 6 boundary,
the authoritative status can be advanced without changing the selected P6.1
scope.

### CH5-AUDIT-007 – P5.3 operator provenance does not close over bound coefficients

**Chapter 6 impact:** does not block starting P6.1. Resolve before P6.2 so a
supplied optimality system cannot inherit ambiguous transposition-operator
provenance.

The P5.3 factories inherit the scalar diffusion-reaction residual and the
deal.II model assembles the bound positive diffusion and reaction values. The
typed transposition enum names a scalar diffusion-reaction Dirichlet-Laplacian,
but its stable isomorphism ID and the factory, guide, and roadmap prose still
describe `T=-Delta+rI`. The compiler neither restricts the diffusion binding
to one nor records coefficient-data ports in the transposition selection, and
the focused tests use unit diffusion only.

The P5.3 acceptance item requiring operator provenance to agree with accepted
diffusion/reaction data therefore remains unproved. A bounded closure can
either restrict these registered factories to unit diffusion or make the
typed operator explicitly reference the diffusion/reaction data and render
`T=-kappa Delta+rI`, with one non-unit binding regression. A general symbolic
operator representation is unnecessary.

## Closed issues

None.
