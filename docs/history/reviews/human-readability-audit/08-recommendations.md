# Audit recommendations and proposed refactor sequence

**Audit baseline:** [`f53b7f009e5c`](https://github.com/HunterNise/nmopt-project/commit/f53b7f009e5c418ec4f3855db29f7eb924faf6f1)  
**Status:** final audit handoff for the documentation refactor. The implementation remains frozen; deeper source/API restructuring is outside this audit's evidentiary scope.

## 1. Main conclusion

The codebase does **not** justify a broad architectural rewrite **as part of
this human-readability phase**.

This is a bounded conclusion. The audit inspected architecture and source
readability deeply enough to reject obvious or inexpensive implementation
refactors, but it was not a dedicated whole-repository source-maintainability
study. A later decision to restructure headers, namespaces, APIs, or ownership
boundaries should begin with a focused audit of maintenance cost, dependency
structure, compile cost, edit patterns, and reviewer workflow.

The strongest evidence for keeping implementation frozen in this phase is:

- optimizer/formulation code is backend-neutral;
- the semantic/compiler path and external-application path are genuinely independent producers;
- the external Step-4 experiment reached the shared formulation boundary without modifying core APIs;
- KKT/PDAS/supplied-OTD are distinct typed products rather than PDE-specific branches in the reduced solver;
- the project runner is visibly an outer application/research shell.

The main human-readability problem is therefore:

> **inspectability and authority presentation**.

The recommended refactor is documentation-only. Small source-comment cleanup
candidates are recorded in `05-source-readability.md` but deliberately deferred
from `codex/docs-refactor`.

## 2. Priority 0 — freeze numerical behavior

During the readability refactor:

- no algorithm changes;
- no compiler-capability expansion;
- no new benchmark family;
- no new external PDE experiment;
- no API redesign based only on aesthetics.

Allow within this branch:

- stale-link fixes;
- documentation routing;
- authority/status corrections;
- new explanatory/reference documentation.

If a correctness problem is discovered while documenting the code, record it
and handle it separately rather than silently expanding this branch into an
implementation refactor.

This keeps the audit from becoming another development phase.

## 3. Priority 1 — create the overview layer

Create:

```text
docs/overview/project-architecture.md
docs/overview/semantic-compiler.md
docs/overview/numerical-realization.md
docs/overview/reduced-optimization.md
docs/overview/research-execution.md
```

Write these from implementation evidence, not by copying existing prose wholesale.

### First file to write

`project-architecture.md`

Reason: every later cleanup decision depends on the canonical project story.

It should establish:

```text
two producer paths
 -> shared mathematical/formulation boundary
 -> primary reduced optimization
 -> secondary KKT/OTD/PDAS capabilities
 -> optional research execution shell
```

## 4. Priority 2 — repair top-level navigation

### Root README

Update:

- project identity;
- current implemented scope;
- bounded limitations;
- “Where to start.”

Do not send new readers into a roadmap.

### `docs/README.md`

Reduce from a catalogue to a router.

Keep planning/review links available under a clearly named history/evidence section.

Separate:

```text
human project understanding
human API/use
research reproduction
agent implementation workflow
```

## 5. Priority 3 — establish a current numerical-contract reference

Create:

```text
docs/reference/numerical-contracts.md
```

This should own the exact current cross-producer solver/formulation surface:

- layouts/primal/covector;
- executable model;
- solve services/reports;
- metrics/constraints;
- reduced DTO;
- reduced Hessian;
- reduced search/trust region;
- supplied OTD;
- quadratic KKT;
- complementarity/PDAS;
- backend-parametric convention and tested backends.

### Why

The current exact API is still routed through a file called:

```text
implementation/v0/executable-contract.md
```

even though that document now describes a much larger current surface.

### Preserve history

Do **not** delete or rewrite the v0 record as though it was always current.

Instead:

- add a historical banner;
- link to the new current reference;
- update inbound links that mean “current numerical API.”

## 6. Priority 4 — reconcile concrete stale documentation

### 6.1 Repository organization decision

Fix:

- `parameters/` is implemented, not “reserved future”;
- actual test directories; remove nonexistent `tests/benchmark/`;
- current application/reference organization if needed.

### 6.2 Documentation map

Fix:

- “proposed” parameter-file schema -> implemented B1/B2 scope;
- default routes away from mutable roadmaps;
- historical v0 lowerer links;
- new overview links.

### 6.3 Implementation roadmap status

The current `implementation-roadmap.md` and `application-roadmap.md` disagree on B1/B2 status.

If the project is frozen, preferred treatment:

```text
implementation-roadmap.md
  Historical implementation ledger.
  Current project scope/status is summarized in the overview/current-scope page.
```

Do not spend effort keeping an 85 KB “next task” ledger current when there is no next implementation task.

### 6.4 Root scope language

Replace “first scalar slice / still being extended” with a bounded statement of what is actually implemented and what is deferred.

## 7. Priority 5 — clean design/document role leakage

### `interface-specification.md`

Remove or visibly isolate implementation-generation references such as P5.1/Chapter-specific current implementation notes from normative clauses.

Move those notes to:

- v1 capability record;
- application/benchmark docs;
- historical review.

### `architecture.md`

Add explicit “design intent / target architecture” framing and route current readers to the overview.

### `implementation-readiness-review.md`

Mark it clearly as a historical/default-policy review rather than current implementation authority.

### `system-blueprint.md`

After overview files exist, either:

1. retain it as a deep “implementation blueprint” and remove claims that it is the shortest first-read path; or
2. split material that became redundant and keep only implementation correspondence/source tour.

Option 1 is lower risk.

## 8. Priority 6 — defer source changes from this branch

The audit found a few cheap comment-only candidates in reusable public code,
including development-history shorthand such as:

```text
P5.1
P6.2
C5.6
selected P5.x product
Chapter-specific implementation-unit names
```

Those are genuine readability observations, but they are not important enough
to mix implementation files into `codex/docs-refactor`.

Record them for a possible later source-cleanup pass. Do not change C++ headers,
symbols, namespaces, include structure, or file organization in this phase.

If broader source refactoring is later desired, first perform the dedicated
source/API-maintainability audit described in
`05-source-readability.md`.

## 9. Priority 7 — do not split large headers yet

Large public headers are a real source-navigation cost, particularly:

- semantic validation;
- deal.II compiler;
- fixed-Dirichlet realization;
- large integration tests.

But the audit has not found a correctness or ownership failure caused by their physical size.

Therefore:

### First response

- add overview source maps;
- improve section-level navigation/comments where missing;
- ensure responsibility boundaries are named.

### Only split if one of these becomes true

- maintainers repeatedly edit unrelated responsibilities in the same file;
- include-time/build-cost evidence points to a problem;
- a file contains independently testable units with stable dependency direction;
- review routinely cannot isolate changes.

Avoid a large mechanical file move immediately before presentation.

## 10. Priority 8 — leave naming/API compatibility alone unless evidence demands change

Do not rename immediately:

```text
dealii_backend
application
StateControlPartitionT
CompiledProblemT
```

The docs can clarify their roles.

Potential future ergonomic work, deliberately out of this refactor:

- narrower reduced-only callback construction;
- callback metric adapter;
- helper for layout/solve translation;
- further interface segregation of `ExecutableModelT`.

The Step-4 evaluation already established these as ergonomics questions rather than architectural blockers.

## 11. Doxygen/generated API documentation

Do **not** add generated API documentation in this refactor.

Reasons:

- exact Markdown references already exist;
- much of the important explanation is mathematical/ownership-oriented rather than signature-only;
- another generated surface would create another authority layer;
- the project is header-heavy and would require a substantial annotation pass for good output.

Reconsider generated API docs only if nmopt becomes a separately packaged external library with users who need exhaustive symbol browsing.

## 12. Suggested implementation batches

### Batch 0 — preserve audit evidence

- add this audit under `docs/planning/review/human-readability-audit/`;
- link it from `docs/planning/review/README.md`.

### Batch A — canonical understanding

- add `docs/overview/project-architecture.md`;
- add `docs/overview/reduced-optimization.md`;
- add `docs/overview/semantic-compiler.md`;
- add `docs/overview/numerical-realization.md`;
- add `docs/overview/research-execution.md`.

### Batch B — current reference

- create `docs/reference/numerical-contracts.md`;
- reroute “current executable contract” links away from
  `implementation/v0/executable-contract.md`;
- mark v0 contract/lowerer history explicitly.

### Batch C — routing

- rewrite root README scope/start section;
- rewrite `docs/README.md`;
- update major references to the new overview where appropriate.

### Batch D — stale-authority repair

- repository organization decision;
- implementation-roadmap historical framing;
- design/reference authority banners;
- documentation-map stale wording;
- any broken status cross-links discovered during link audit.

### Batch E — consistency validation

- Markdown link check;
- repository search for stale phrases:
  - `proposed .prm`;
  - `first scalar finite-element slice`;
  - `still being extended`;
  - `reserved future`;
  - current routes to v0 as present API;
  - agent-facing wording outside intended agent docs;
- no C++ build is required if the branch remains documentation-only.

### Batch F — audit closure

- add `docs/planning/review/human-readability-audit/closure-report.md`;
- record which recommendations were implemented, revised, or deferred;
- record final link/consistency validation.

## 13. Proposed stop condition

The humanization refactor is done when all of the following are true:

1. A new reader can understand the system from one primary overview without opening a roadmap.
2. External-application and semantic/compiler users can identify their separate entry paths.
3. The reduced runtime can be followed from state solve through optimizer from one overview.
4. “Backend” is not conflated with FE realization.
5. KKT/OTD/PDAS are visible without dominating onboarding.
6. Chapter 6 runner/artifact complexity is clearly optional outer infrastructure.
7. Every major current API claim routes to current code/reference, not v0 history.
8. Historical roadmaps/reviews remain accessible but are visibly historical.
9. Root/documentation READMEs no longer contain known stale status claims.
10. Any source-level readability issues intentionally left unchanged are recorded as deferred rather than implied to be resolved.

## 14. Recommendation severity matrix

| Area | Recommendation |
| --- | --- |
| Core numerical architecture | **Keep** |
| External integration boundary | **Keep; document prominently** |
| Semantic/compiler architecture | **Keep; explain bounded composition accurately** |
| KKT/PDAS/OTD products | **Keep; secondary onboarding** |
| Application/runner architecture | **Keep; present as optional research shell** |
| Root/doc navigation | **Major documentation refactor** |
| Current-vs-history authority | **Major documentation reconciliation** |
| Public code comments | **Defer from docs branch; cheap future cleanup candidate** |
| Header/file organization | **No move now; dedicated audit required before restructuring** |
| Public API ergonomics | **Defer; evidence exists but not a blocker** |
| Generated API docs | **Do not add now** |
