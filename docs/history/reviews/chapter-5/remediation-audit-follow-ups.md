# Chapter 5 remediation audit follow-ups

## Purpose and scope

This file records issues found while auditing the Chapter 5 remediation
commits against the acceptance checklists in this directory. The audit is
bounded to deciding whether implementation may proceed to Chapter 6. It does
not reopen deferred architecture or unselected capabilities.

The [implementation roadmap](../../implementation-roadmap.md) remains the
status ledger. Checked acceptance items remain in their owning review files.

## Resolved issues

### CH5-AUDIT-001 – the closed compiler request is now the sole decision source

**Chapter 6 impact:** resolved for the accepted C2 scope.

Semantic resolution and lowerability close one request before construction.
That request supplies target and Dirichlet registration, planning, model
construction, mesh checks, and final provenance; scalar execution consumes its
selected plan records. `finalize_resolved_decision()` no longer receives the
raw graph, target enum, second registration, or raw region objects. A
declaration-order permutation now compares the complete compatibility record,
not a selected subset. The regression and full neutral/deal.II baselines pass,
closing the C2 request-source item without expanding the compiler framework.

### CH5-AUDIT-002 – the manifest is now a pure resolved-decision projection

**Chapter 6 impact:** resolved for the accepted C1 scope.

The pre-construction decision now owns typed regions, requirement records,
component inventories, and unresolved map skeletons. Finalization fills their
runtime dimensions and service records using only that decision, the closed
request, and the constructed services. Compiler-selected requirements render
stable typed ID/kind/status/scope/region provenance rather than copying
`selected_policy`; only a `user_assumed` record retains its model-author
declaration. `make_manifest()` is a pure projection of the finalized decision.

Schema 3 and its structured tests cover scalar, binding, mesh, and realized-
map/space records. Point-sensor, normal-flux, and P5.4 display-prose edits now
compare the complete compatibility record, and declaration-order permutation
does the same. The exact projection and full-baseline checks close the C1
manifest item.

### CH5-AUDIT-003 – boundary-trace realized-map ownership and dimension are closed

**Chapter 6 impact:** resolved for the accepted C1 scope.

Boundary-trace models now own their ordered face-sample maps and dimensions,
and expose equivalent value/JVP/VJP actions. Manifest realized maps and spaces
use those records; weighted and unweighted trace tests compare output counts,
maps, and pairings. The focused and full baselines pass, closing the ownership
and dimension item.

### CH5-AUDIT-004 – typed transposition capability is independent of policy prose — resolved

**Chapter 6 impact:** resolved for the accepted P5.3/P5.4 scope; the supplied-
OtD boundary retains the typed formulation and transposition provenance.

P5.3 and P5.4 now share typed transposition, partial-boundary, fractional-
metric, and boundary-$H^{1}$ selections. Semantic and compiler checks reject
unsupported enum selections, and the resolved decision and manifest retain
the typed payloads.

The compiler lowerability path now consumes the resolved typed transposition,
boundary, metric, and target selections; it does not search `selected_policy`
for English fragments. Semantic selection presence is likewise determined by
the selected status/scope and typed payloads. A model-author domain-regularity
declaration remains explicit semantic input, while its compatibility provenance
is rendered from the typed requirement identity and classification.

The point-sensor, normal-flux, and Section 5.11 deal.II regressions clear the
relevant policy descriptions before validation and compilation and compare the
complete compatibility rendering with an unchanged baseline. The full
`debug-neutral` (12/12) and `debug-dealii` (41/41) baselines pass, so the P5.3
prose-independent-validation and P5.4 prose-independent-factory acceptance
items are closed. Compatibility rendering remains a manifest concern already
addressed by `CH5-AUDIT-002`.

### CH5-AUDIT-005 – P5.4 realized-layout and Chapter 5.11.2 Taylor evidence is complete — resolved

**Chapter 6 impact:** resolved for the accepted P5.4 scope; the structured
layouts and Taylor evidence are available for Chapter 6 benchmark reporting.

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

### CH5-AUDIT-006 – the roadmap and review evidence are reconciled

**Chapter 6 impact:** resolved.

The owning C1/C2 and P5.1–P5.4 checklists record the completed bounded
remediation, consistent with the implementation roadmap and v1 capability
ledger. The focused evidence and full `debug-neutral` (12/12),
`debug-dealii` (41/41), and `sanitize-neutral` (12/12) baselines pass. No
Chapter 5 audit item remains open; implementation may proceed to the selected
Chapter 6 work within its stated scope and gates.

### CH5-AUDIT-007 – P5.3 operator provenance closes over bound coefficients

**Chapter 6 impact:** resolved for the accepted P5.3 scope.

Typed transposition selections now carry diffusion and reaction data-port
provenance, and semantic/compiler validation binds those ports to the residual
data. The manifest retains the typed provenance. The non-unit regression now
solves the manufactured zero-control state for both unit and non-unit
diffusion, with forcing scaled consistently, and requires equal states plus a
vanishing residual. It therefore exercises the realized
`T=-kappa Delta+rI` action rather than only inspecting its manifest text.
Focused and full baselines pass, closing the operator-provenance item without
introducing a general symbolic operator representation.

## Deferred issues

None.
