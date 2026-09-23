# Human-readability refactor closure

**Original audit baseline:** `f53b7f009e5c418ec4f3855db29f7eb924faf6f1`

**Documentation branch:** `codex/docs-refactor`

**Repository-structure interlude:** closed separately by [`10-repository-structure-refactor-closure.md`](10-repository-structure-refactor-closure.md)

**Status:** complete

## 1. Purpose and closure verdict

This report closes the human-readability refactor that began with the
code-derived audit at `f53b7f009e5c`.

The original audit did not find a repository with no documentation or an
obviously broken numerical architecture. It found a repository whose
architecture, implementation history, planning records, exact APIs, experiments,
and agent workflow were all documented, but whose **default human reading path
did not make the authority of those documents sufficiently clear**.

The central problem was therefore **inspectability and authority presentation**.

The refactor addressed that problem in several phases. It first established a
human-facing project explanation and manual, then reorganized documentation by
role. A separate repository-structure refactor was deliberately inserted while
source paths and ownership boundaries were still unstable. After that structural
work closed, the documentation refactor resumed against the settled tree:
current public references and maintainer internals were rewritten, foundational
and historical material was reclassified rather than erased, and a final
repository-wide consistency pass reconciled terminology, Markdown rendering,
diagrams, and links.

The resulting repository preserves the architecture reconstructed by the
original audit:

```text
semantic/compiler producer             application-owned producer
        │                                       │
        ▼                                       ▼
semantic resolution                    native numerical realization
        │                                       │
        ▼                                       │
deal.II compiler                       thin numerical adapters
        │                                       │
        └──────────────────┬────────────────────┘
                           ▼
                shared numerical contracts
             model / solves / metric / constraint
                           │
                           ▼
                  formulation products
          reduced DTO / supplied OTD / KKT / PDAS
                           │
                           ▼
                         solvers
                           │
                           ▼
             optional application/evidence shell
```

The two producer paths remain peers. They do not converge at `ProblemSpec`,
`DealiiCompiler`, or `CompiledProblemT`. Backend algebra remains distinct from
finite-element realization. The reduced state–adjoint route remains the clearest
shared runtime path, while supplied OTD, KKT, complementarity, and PDAS remain
visible as real but secondary formulation products. The Chapter-6 runner and
reproduction machinery remain optional project infrastructure rather than a
required layer of the reusable numerical library.

The final two historical Step-4 section links were pinned to their original
pre-refactor targets, and the resulting documentation diff passed its final
whitespace check. The human-readability refactor is complete. No further broad
documentation or source cleanup is required for this audit.

## 2. Evidence and interpretation

This closure is not a rewrite of the original audit and is not based on the
final commit log alone. Different records answer different historical
questions.

The main tracked records used for this reconstruction are:

- [`00-audit-index.md`](00-audit-index.md), the original audit entry point;
- [`05-source-readability.md`](05-source-readability.md), the bounded source-readability conclusion;
- [`06-documentation-audit.md`](06-documentation-audit.md), the original authority/routing assessment;
- [`07-human-information-architecture.md`](07-human-information-architecture.md), the first proposed reader-facing structure;
- [`08-recommendations.md`](08-recommendations.md), the original implementation handoff and stop criteria;
- [`09-repository-structure-follow-up.md`](09-repository-structure-follow-up.md), the later dedicated source/test/application/tooling assessment;
- [`10-repository-structure-refactor-closure.md`](10-repository-structure-refactor-closure.md), the closure and validation record for that interlude; and
- [`../../documentation-lineage.md`](../../documentation-lineage.md), the tracked map of major documentation transitions.

The evidence hierarchy used here is:

```text
final repository
    authority for what exists now

original human-readability audit
    why the refactor began and what the initial evidence supported

human-readability continuation handoff
    state and intended continuation after the manual phase

repository-structure follow-up and roadmap
    why documentation work paused and what source structure was reconsidered

reconciled pre-closure structure roadmap
    strongest reconstruction of why the final R1–R13 decisions were accepted

repository-structure closure
    what that bounded subproject actually completed and validated

historical conversations and agent handoffs
    decision provenance, alternatives considered, and explicit deferrals
```

The original audit remains baseline-pinned evidence. It has not been rewritten
to make its recommendations look identical to the final implementation.

Likewise, temporary agent plans and handoffs remain source material for this
closure rather than becoming additional tracked authorities. Durable
documentation history is preserved through
[`../../documentation-lineage.md`](../../documentation-lineage.md), named
historical records, and Git itself.

This distinction matters because several good final decisions deliberately
departed from earlier proposals once more evidence existed.

## 3. Architectural invariants preserved throughout

Across both the documentation and repository-structure phases, several
distinctions were treated as architectural invariants rather than cleanup
targets:

- primal values versus covectors;
- derivatives versus metric-dependent gradients;
- operator actions versus inverse/solve operations;
- residual, JVP, and VJP actions;
- state solves versus adjoint solves;
- numerical values versus solve/convergence/evidence records;
- compiler-owned products versus application-owned numerical realization;
- reusable numerical/formulation contracts versus runner/reproduction
  infrastructure;
- distinct reduced DTO, supplied OTD, quadratic KKT, complementarity, and PDAS
  products;
- bounded compiler composition rather than an arbitrary Cartesian product of
  semantic components;
- no universal formulation interface.

These invariants explain several decisions that would otherwise look like
unfinished simplification. The project was made easier to inspect without
flattening distinctions that carry mathematical, ownership, or verification
meaning.

# Phase I — Original audit and target reading model

## 4. The original diagnosis

The audit reconstructed the repository from public headers, representative
runtime paths, tests, build dependencies, application code, compiler behavior,
and existing documentation.

Its central finding was that the project already had substantial documentation,
but new readers encountered several roles before knowing how to distinguish
them:

- current implemented architecture;
- foundational or target design;
- implementation-generation records;
- current public programming interfaces;
- mutable or completed roadmaps;
- review and remediation evidence;
- source/application studies;
- agent workflow.

The original audit therefore recommended a documentation-only first phase. It
explicitly froze numerical behavior and rejected broad API, namespace, or
physical source restructuring without a separate source-maintainability study.

This was a bounded conclusion, not a claim that the physical repository could
never benefit from restructuring. That distinction became important later.

## 5. The canonical reader model

The audit proposed that a new reader should be able to answer, without opening
a roadmap:

1. what `nmopt` is;
2. what reusable numerical boundary it exposes;
3. how mathematical problems reach that boundary;
4. what deal.II owns;
5. what formulations and algorithms exist;
6. what is optional research/reproduction infrastructure;
7. where exact current programming information lives;
8. which records are historical rather than current authority.

That reader model survived the refactor. The final information architecture is
more elaborate than the first proposal, but it still serves these questions.

# Phase II — Manual and onboarding

## 6. Canonical overview and progressive manual

The first substantial implementation phase built the human-facing explanation
of the project.

The initial standalone overview layer was later moved under:

```text
docs/manual/overview/
```

so that first-time architectural orientation and the longer concept progression
formed one reader-facing hierarchy.

The overview now supplies a fast mental model of:

- the project architecture;
- the semantic/compiler producer;
- the application-owned producer;
- numerical realization;
- reduced optimization;
- application/experiment evidence;
- the major architecture maps and reading routes.

The concept manual then developed into a deliberate thirteen-chapter sequence:

```text
Part I   – numerical language and duality
Part II  – formulations and optimization algorithms
Part III – semantic authoring, validation, compilation, and lowering
Part IV  – using compiled problems and integrating an existing application
```

This was substantially broader than the original five-overview recommendation,
but it did not arise by simply adding documentation volume. Several proposed
chapters were removed or reshaped when review showed that they duplicated
existing explanations.

### 6.1 Standalone numerical-realization chapter rejected

A proposed concept chapter between compilation and application integration was
dropped because its useful content already belonged in earlier numerical
chapters and the compilation/lowering chapter.

Independent coordinates, reconstruction and pullback, observation transpose
realizations, metric/Riesz services, solve separation, and the deliberately
narrow backend algebra were integrated where they had conceptual context
instead of being repeated in a separate catalogue.

### 6.2 Verification/evidence ending replaced by actual usage paths

An early idea was to end the manual with a verification/evidence chapter.

That was rejected because Taylor checks, JVP/VJP checks, state/adjoint
consistency, KKT validity, and PDAS residual evidence were already introduced
where those concepts matter.

Part IV instead closes the two actual user stories:

```text
compiled path
    recipe / parameters
        ↓
    ProblemSpec
        ↓
    validation / runtime bindings / compiler
        ↓
    compiled product
        ↓
    formulation + solver
        ↓
    result / output / manifest

application-owned path
    existing PDE application
        ↓
    native numerical mathematics
        ↓
    callback / solve / metric adapters
        ↓
    shared formulation contracts
        ↓
    solver
        ↓
    result returned to the application
```

### 6.3 Naming and writing conventions

The semantic formulation type was renamed from
`ReducedFormulationSpec` to `FormulationSpec`; current explanatory
documentation uses the final name without teaching the migration history.

The manual phase also established writing principles that influenced later
documentation:

- teach the mathematical or software problem before naming the project
  abstraction;
- explain mathematical meaning before implementation names;
- use derivations when comparison of identities matters;
- avoid long abridged source dumps;
- distinguish framework capability from benchmark/application coverage;
- avoid editorial meta-commentary and microheading-heavy prose;
- preserve mathematically meaningful distinctions rather than simplifying
  terminology for visual neatness.

# Phase III — Documentation reorganization by role

## 7. From working categories to authority categories

After the manual had a stable reader path, the broader documentation tree was
reorganized by long-term role.

The important historical boundary is:

```text
640d5a3    last convenient pre-role-reorganization snapshot
997367f    docs(structure): reorganize documentation by role
```

The final top-level roles became:

```text
manual/       teach the current project
design/       long-lived architecture, mathematics, and accepted decisions
reference/    current practical public programming interface
internals/    current implementation mechanics and maintainer maps
studies/      Chapter 5/6 source, application, benchmark, and case-study corpus
planning/     retained development roadmaps and planning records
history/      superseded implementations, reviews, audits, and evidence
```

This was not a mass rewrite of document bodies.

Examples of role reclassification included:

```text
design/theoretical-formalism.md
    → design/mathematical-model.md

design/system-blueprint.md
    → internals/system-blueprint.md

implementation/v1/semantic-compiler.md
    → internals/compiler/semantic-compiler.md

implementation/v0/*
    → history/implementation/*

planning/review/*
    → history/reviews/*
```

Chapter-specific source, benchmark, guide, and case-study material was brought
under `studies/`, while accepted decisions became part of `design/`.

The preservation rule that eventually governed the rest of the refactor began
to solidify here:

> change a document's current role and routing when needed, but do not rewrite
> historical evidence merely so that it reads as though it was written after
> project completion.

# Phase IV — Repository-structure readability interlude

## 8. Why documentation work paused

After the role-based documentation reorganization, a dedicated source,
test, application, tooling, and build-structure review was opened.

This did not contradict the original documentation-only audit. The original
source-readability conclusion was explicitly limited to the evidence needed for
the documentation phase. The later follow-up was a deeper repository-wide
maintainability assessment.

Broad reference and documentation rewriting paused while source paths and
ownership boundaries stabilized.

The structure refactor's objective was:

> make the physical repository reflect the architecture that already exists,
> while preserving behavior unless a separately identified correctness defect
> requires a deliberate correction.

The structure work was not a compiler rewrite, API simplification project, or
numerical redesign.

## 9. Source and test responsibility locality

Several source and test boundaries were changed where the responsibility was
real and reviewable.

Current runtime/compiler diagnostics moved away from roadmap-era `P5.x`,
`P6.x`, `C5.x`, and `v0` wording when a user otherwise needed project-history
knowledge to understand present behavior. Historical IDs remained in provenance,
reviews, studies, and intentionally stable test traceability.

Neutral tests were separated by responsibility without unnecessarily
multiplying translation units or CTest identities.

The former large deal.II capability ledger was reorganized into:

```text
tests/compiler/dealii/
    one heavy compiler translation unit
    capability-grouped scenario implementation headers

tests/dealii/
    lower-level backend/service/session/reporting coverage
```

Keeping one heavy compiler translation unit was deliberate. Recompiling the
large header-only compiler and deal.II surface in many translation units would
increase build cost without improving the logical verification inventory.

## 10. Application and runner ownership

`apps/nmopt-runner/main.cc` had become an implementation dependency of tests.

The chosen correction was not to invent a generic execution framework. One
cohesive private seam was extracted:

```text
apps/nmopt-runner/chapter6_execution.hpp
```

It owns Chapter-6 execution, related artifact/evidence enrichment, and execution
dispatch. The executable entry point retains CLI, configuration, run
preparation, and lifecycle concerns.

Its size is therefore not treated as a remaining defect by itself.

## 11. Tooling and evidence semantics

Two structure units intentionally changed behavior or evidence interpretation
and must not be hidden under the label "file organization."

### 11.1 Plotting-profile correctness

The accepted plotting-profile design said that generic rendering should honor
resolved profile policy rather than retain Chapter-specific presentation
branches.

The refactor corrected the implementation so that declared profile properties
actually control rendering, while preserving explicit override precedence.

### 11.2 Persisted-artifact parsing

Artifact/history parsing was consolidated under a dependency-light shared
package.

The motivation was not only duplicate code. Two parsers could interpret
persisted evidence differently. The stricter existing malformed/non-finite
history interpretation became canonical.

## 12. Semantic/compiler boundaries and naming

Semantic validation acquired a narrow public facade with its order-sensitive
implementation moved under `detail/`.

The implementation was deliberately not fragmented further because validation
order and helper coupling are meaningful and smaller files alone would not
reduce template/include cost.

Compiler manifest/provenance schema was separated from executable compiled
products along a stronger natural boundary. `compiled_problem.hpp` remains the
aggregate compatibility surface.

Several names and placements were also changed after explicit ownership review:

```text
reference_specs.hpp
    → problem_library.hpp

NativeApplicationViewT
    → CompiledApplicationViewT

runner.hpp
    → run_lifecycle.hpp

include/nmopt/reference/*
    → tests/support/reference_models/*
```

These changes did not all have the same compatibility policy.

The algebraic reference models were intentionally demoted from apparent public
framework API to test/oracle support. Recreating public forwarding headers
would undo that ownership decision.

`CompiledApplicationViewT` was renamed because "native application" had acquired
a precise meaning for the peer application-owned producer path.

The `reference_specs.hpp → problem_library.hpp` include migration was also
deliberate and repository-wide. Public problem-factory function names remained
stable, but the old include path was not retained.

## 13. Deliberate structural stopping points

The structure refactor also made several explicit **no-change** decisions.

Only the exceptional external Step-4 registration block was extracted to:

```text
cmake/ExternalStep4.cmake
```

Further CMake modules were rejected because the remaining root file usefully
shows project policy and the ordinary build/test graph.

No new `apps/` taxonomy was introduced because it would mostly have produced
one-item speculative directory levels.

For Step-4:

- canonical minimal bindings were separated from evaluated/instrumented
  comparison bindings;
- repeated adapted-source inclusion knowledge was consolidated;
- preserved source snapshots were not restructured;
- an instrumentation observer/event framework was rejected;
- a shared Problem A/B adapter mini-framework was rejected.

The distinction established by this review is important:

> minimal in dependencies and experiment machinery is not the same as minimum
> possible source line count.

Finally, `dealii_compiler.hpp` remained intentionally large. The review found
plausible physical regions but not stable enough typed phase boundaries to make
a mechanical split an architectural improvement. Coarse source navigation was
added instead.

The structure branch closed separately with neutral, deal.II, and supported
Python/tooling validation passing, then handed the broad documentation work
back to `codex/docs-refactor`.

# Phase V — Resumed current public/reference authority

## 14. Task-oriented reference replaces the earlier single-page proposal

Before the structure interlude, the continuation plan proposed:

```text
docs/reference/numerical-contracts.md
```

because current solver/formulation information was still routed through a
historically named `implementation/v0/executable-contract.md`.

After source paths stabilized, the reference problem was reassessed more
broadly.

The final reference philosophy became:

> Reference is the practical programming interface: guided enough to build
> with, precise enough to look things up in.

The unit of organization is a programming task, not a class or one giant
contract catalogue.

The final reference suite is:

```text
problem-authoring.md
compiler.md
application-authoring.md
external-dealii-solver-integration.md
optimization.md
application-execution.md
parameter-files.md
```

The previous `application-api.md` was retired.

This seven-page suite covers the current numerical/formulation contracts where
users actually need them: authoring, compilation, external integration,
optimization/formulation choice, execution, and configuration.

A standalone `numerical-contracts.md` would now duplicate these task-oriented
references and pull the documentation back toward class/catalogue organization.
Its absence is therefore an intentional revision of the original plan, not
unfinished work.

## 15. Producer-path authority clarified

The current reference and overview material consistently distinguish:

```text
compiler-owned path
    semantic graph → validation/resolution → lowering → numerical realization

application-owned path
    existing numerical application → thin actions/solves/metric adapters
```

They converge downstream at numerical/formulation contracts.

`CompiledApplicationViewT` is a compiler-side application-facing sidecar, not
the universal representation of an application. Likewise
`HeadlessBenchmarkRunnerT` is generic over its builder result even though the
repository's current native applications commonly use the compiler-backed
`ProblemSpec` path.

The external Step-4 evidence is kept intentionally narrow: it demonstrates the
application-owned producer route for the reduced formulation and does not imply
external/native construction parity for every compiler formulation product.

# Phase VI — Current implementation documentation and lineage

## 16. Replacing broad legacy internals

The role-reorganization initially moved two broad implementation documents into
`internals/`:

```text
internals/system-blueprint.md
internals/compiler/semantic-compiler.md
```

Later review showed that both still mixed too many authorities.

`system-blueprint.md` combined teaching, architecture correspondence, source
tour, verification guidance, and current implementation statements.

`compiler/semantic-compiler.md` combined public compiler capability/reference
material with private implementation mechanics.

The final current internals are instead:

```text
internals/implementation-map.md
internals/compiler.md
internals/runner.md
```

Their roles are intentionally narrower:

- `implementation-map.md` maps architectural responsibilities to current source,
  ownership, lifetimes, tests, and deliberate structural stopping points;
- `compiler.md` explains validation/resolution, closed requests, planning versus
  specialized targets, runtime bindings, numerical realization, products,
  manifests, lifetimes, and compiler verification;
- `runner.md` explains the repository-specific application shell, parameter
  resolution, run planning, lifecycle/manifests, execution registration,
  Chapter-6 orchestration, and post-processing boundary.

Public programming questions remain in `reference/`.

The superseded internals were deleted rather than copied into a second archive.
Their evolution and replacement are recorded in
[`../../documentation-lineage.md`](../../documentation-lineage.md), and exact
older versions remain recoverable through Git.

## 17. Retrospective architecture draft deliberately not added

A separate `design/implemented-architecture-draft.md` was considered.

After the three current internals existed, it was rejected as part of this
refactor.

The implementation map already gives the clean retrospective correspondence
between producer paths, numerical/formulation convergence, application/evidence
layers, ownership, source, tests, and deliberate structural decisions.
Another broad architecture synthesis would sit awkwardly between:

```text
manual/
design/
internals/implementation-map.md
```

and recreate the overlapping-authority problem the refactor was removing.

An end-of-project architecture essay may still be useful someday, but it is not
required current documentation authority.

# Phase VII — Reconcile durable design, studies, planning, and history

## 18. Foundational design: preserve intent, frame the time period

The final policy for foundational design records is not to modernize them into
an artificial present tense.

Documents such as:

```text
design/architecture.md
design/composition-boundaries.md
design/interface-specification.md
design/pde-solver-boundary.md
```

remain valuable partly because they record the design model and reasoning that
guided development.

They now carry completion-era framing that distinguishes:

- normative or foundational design intent;
- the exact finished implementation;
- current programming/reference authority.

Old P5/P6 or implementation-generation vocabulary may therefore remain inside
the body when it is part of the historical design context. The completion note
is the mechanism that prevents readers from mistaking those passages for an
exact inventory of current C++ support.

Objectively false current-routing statements were corrected. The historical
voice itself was not scrubbed.

`design/mathematical-model.md` remains the durable mathematical/sign-convention
authority.

## 19. Studies: current evidence without roadmap ownership leakage

`studies/` now owns the Chapter-5/6 source/application/benchmark/case-study
corpus.

The final review kept scientific and reproduction evidence in place while
clarifying which page owns which kind of fact:

- source facts and extracted numerical methods;
- scenario definitions;
- frozen benchmark choices;
- observed B1/B2 replication outcomes;
- framework behavior owned by reference;
- historical intended B3–B6 sequencing owned only as retained planning.

One substantive authority ambiguity was corrected for B1: the authoritative
profile's continuous volume control representation is distinguished from a
different default in a scenario-construction helper.

No broad rewrite of B1/B2 evidence bodies was performed.

## 20. Planning: retained plans, not the current work queue

The project keeps planning documents under `planning/` because they remain plan
artifacts even when inactive.

They now state explicitly that they are retained records of intended sequencing,
handoffs, gates, and extensions considered during development, not a current
queue merely because of their directory.

This avoids two bad alternatives:

1. continuously rewriting old roadmaps until their chronology disappears;
2. moving every completed plan into `history/` and losing the distinction
   between plan artifacts and review/evidence artifacts.

If substantial development resumes, current capability and source state must be
reconstructed before any old remaining task is treated as active.

## 21. History: preservation and provenance

Historical implementation records, reviews, audits, and retained drafts remain
historical.

The final preservation policy is:

- use Git for earlier versions of documents whose current identities remain
  recognizable;
- retain named records under `docs/history/` when they have durable value as
  evidence in their own right;
- do not maintain a copied shadow tree of every superseded document;
- do not replace a historical link with a modern document merely because the
  modern document covers similar subject matter.

For broken historical links, provenance takes precedence over conceptual
similarity.

Where the original historical target still exists at the pre-refactor audit
baseline, the durable repair is a permalink to that exact target at:

```text
f53b7f009e5c418ec4f3855db29f7eb924faf6f1
```

This policy was applied to the historical missing-link set after confirming
that all 131 intended targets existed at that single frozen snapshot.

A final historical-anchor check also identified two Step-4 section links whose
relative file target still exists today but whose section names changed. Those
links are to be pinned to the same frozen snapshot rather than silently
retargeted to the reorganized current report.

# Phase VIII — Final repository-wide readability pass

## 22. Conceptual review by documentation role

The final pass did not reopen the information architecture.

Instead, each current documentation role was reviewed against its own purpose:

- manual: teaching sequence and project mental model;
- design: durable architecture/math and accepted decisions;
- reference: exact practical programming workflows;
- internals: implementation ownership and source maps;
- studies: source/application/benchmark/case-study evidence;
- planning: retained plans;
- history: historical evidence and routing.

The result was mostly bounded wording, authority, notation, and routing cleanup
rather than structural rewriting.

## 23. Markdown, mathematics, naming, and diagrams

The project documentation conventions were tightened and then applied to
current/non-historical documentation.

The final work included:

- consistent current project naming as `nmopt`;
- GitHub-safe handling of Markdown-sensitive inline mathematics;
- consistent braced scripts, `\ast`, norm notation, and other established math
  conventions;
- explicit rules for Markdown-table mathematics;
- deliberate Unicode conventions for `text` diagrams;
- removal of diagrammatic ASCII arrows from current/non-historical text
  diagrams;
- preservation of historical typography unless rendering/navigation was
  actually broken.

These changes were mechanical support for readability, not a new prose style
imposed on historical evidence.

## 24. Link and navigation integrity

The final audit separated different kinds of link findings rather than treating
every heuristic candidate as an error.

The current documentation had no remaining missing relative file targets after
the cleanup.

Twenty-three current anchor candidates produced by the checker were
independently verified against GitHub-style heading IDs and all resolved
correctly.

Historical missing targets were repaired against the frozen pre-refactor
snapshot instead of being redirected to modern successor documents.

The remaining historical anchor review found two stale Step-4 section anchors;
those are the final bounded closure repair.

The temporary audit tooling was intentionally local under ignored `tmp/` and
does not become another maintained repository subsystem.

## 25. Where the plan changed

The refactor intentionally diverged from several earlier plans.

| Earlier proposal or expectation | Final decision | Rationale |
| --- | --- | --- |
| Standalone `docs/overview/` as a peer top-level layer | Move overview under `docs/manual/overview/` | Keep one first-time reader hierarchy: overview first, concepts second |
| Small overview-only documentation phase | Build a thirteen-chapter concept manual as well | The project needed a progressive explanation of its mathematical and software vocabulary, not only routing pages |
| Standalone numerical-realization concept chapter | Do not add it | Its useful material already belonged in numerical concepts and compiler/lowering explanation |
| End the manual with verification/evidence | End with compiled use and existing-application integration | Those two chapters close the actual user workflows; verification is taught where each concept arises |
| Create `reference/numerical-contracts.md` | Use seven task-oriented reference pages instead | Contracts are more useful in authoring, compilation, integration, optimization, execution, and configuration context than in one class-style catalogue |
| Keep one broad application API reference | Retire `application-api.md` | Separate programming tasks and producer paths are clearer than one large semantic/application assembly page |
| Keep broad `system-blueprint` and semantic-compiler internals | Replace with `implementation-map.md`, `compiler.md`, and `runner.md` | Separate public use, teaching, implementation mechanics, and source correspondence |
| Add `design/implemented-architecture-draft.md` | Do not add it | `implementation-map.md` already provides the needed retrospective implemented correspondence; another architecture summary would duplicate authorities |
| Treat foundational design role leakage mainly by rewriting old vocabulary | Preserve bodies and add completion-era framing, correcting only objectively false current claims | Historical design perspective is useful evidence and should not be rewritten as if authored after completion |
| Complete optional prerequisite/background chapters before closure | Defer them as a separate pedagogical project | Their absence is not a current architecture/authority defect; the manual already states its assumptions |
| Leave source organization entirely untouched | Run a separate dedicated repository-structure refactor | The original source conclusion was bounded to the docs phase; the later deeper audit found targeted responsibility-locality work worth doing |
| Split large files because they are hard to browse | Split only where a real responsibility boundary exists; otherwise add navigation | Physical size alone does not justify artificial architecture or repeated template compile cost |
| Modernize broken historical links to current equivalent pages | Pin historical references to the original artifact/revision | Preserve provenance and the meaning of the historical record |

## 26. Deliberate stopping points and deferred work

The following are **not unfinished closure tasks**.

### 26.1 Architecture and public API

Deferred unless future evidence establishes a concrete need:

- broad numerical-contract redesign;
- a universal formulation abstraction;
- narrower reduced-only callback construction;
- a callback-style metric adapter;
- convenience helpers for layout/solve translation;
- further interface segregation of `ExecutableModelT`;
- generated Doxygen/API documentation.

The original external-integration evaluation classified these as ergonomics
questions rather than blockers.

### 26.2 Physical source structure

Deliberately left as-is:

- further decomposition of `dealii_compiler.hpp`;
- further splitting of `validation_detail.hpp`;
- reduced-solver public-header reorganization;
- a new physical semantic problem-library taxonomy;
- relocation of `experiment/reduced_envelope.hpp`;
- relocation of Chapter-specific application headers;
- additional CMake modules;
- a new `apps/` taxonomy.

These changes lack a sufficiently strong responsibility, maintenance, compile
cost, or ownership argument at closure.

### 26.3 Step-4

Deliberately not added:

- an observer/event architecture solely to erase the optional instrumentation
  type from the include graph;
- a shared Problem A/B binding mini-framework.

The canonical minimal consumers remain minimal in required experiment/evidence
machinery, not artificially minimum-line adapters.

### 26.4 Optional prerequisite/background teaching

Earlier continuation planning proposed:

```text
docs/manual/background/
    pde-constrained-optimization.md
    optimization-methods.md
    cpp-and-dealii.md
    [optional finite-elements-and-dealii.md]
```

That proposal remains useful future pedagogical work.

It was consciously moved outside this refactor because prerequisite teaching is
additive material, not a defect in current documentation authority. The
project-specific manual and current routing do not depend on those chapters for
closure.

## 27. Final documentation authority map

The finished documentation model is:

```text
Need first-time understanding
    → docs/manual/

Need durable mathematical or architectural rules
    → docs/design/

Need to program against the supported current surface
    → docs/reference/

Need to modify or inspect current implementation mechanics
    → docs/internals/ + source + focused tests

Need Chapter-5/6 source, application, benchmark, or replication evidence
    → docs/studies/

Need development intent, sequencing, or proposed extensions
    → docs/planning/

Need superseded implementation records, reviews, audits, or historical rationale
    → docs/history/
```

`docs/history/documentation-lineage.md` explains major moves and how to recover
superseded documents from Git.

Current source and focused tests remain the exact authority for implementation
behavior when prose and code ever diverge.

## 28. Original audit stop conditions

The original ten stop conditions can now be assessed directly.

| Stop condition | Closure assessment |
| --- | --- |
| A new reader can understand the system from one primary overview without opening a roadmap | **Satisfied.** `manual/overview/` is the default first-read path. |
| External-application and semantic/compiler users can identify separate entry paths | **Satisfied.** Both are first-class throughout overview, manual, reference, and internals. |
| Reduced runtime can be followed from state solve through optimizer | **Satisfied.** The manual and optimization reference trace the state/adjoint/reduced path end to end. |
| Backend is not conflated with FE realization | **Satisfied.** Backend algebra, deal.II numerical services, compiler realization, and application ownership are explicitly separated. |
| KKT/OTD/PDAS are visible without dominating onboarding | **Satisfied.** They are taught and referenced as real formulation products but remain secondary to the reduced route in first-time orientation. |
| Runner/artifact complexity is visibly optional | **Satisfied.** The runner is described as an outer project/reproduction shell, not the reusable numerical boundary. |
| Current API claims route to current code/reference rather than v0 history | **Satisfied.** Current public use routes through the seven reference pages; v0 implementation records are historical. |
| Historical roadmaps/reviews remain accessible but visibly historical | **Satisfied.** Planning has completion/status framing; completed reviews and implementation generations live under history. |
| Root/docs routing no longer contains the known stale status model | **Satisfied.** `README.md` and `docs/README.md` route by audience/task and current authority. |
| Deferred source/API readability work is recorded rather than implied complete | **Satisfied.** The structure closure, internals, lineage, and this report record deliberate stopping points and future ergonomic possibilities. |

A later continuation handoff temporarily added prerequisite/background chapters
as an additional working-plan stop condition. That sequencing was subsequently
changed explicitly: optional background material was separated from the
authority/readability closure. This report records that plan divergence rather
than pretending the proposal never existed.

## 29. Final validation and residuals

The repository-structure phase already closed with its own implementation-level
validation:

- supported Python/tooling contracts passed;
- backend-neutral Debug pipeline passed;
- normal deal.II Debug pipeline passed;
- the final structure cleanup passed `git diff --check`.

The resumed documentation phase did not intentionally change numerical
behavior.

Its final consistency audit established:

- no severity-A Markdown audit findings;
- no remaining missing relative file targets in current documentation;
- all 23 current GitHub-anchor candidates independently verified as valid;
- the remaining current unbraced-script candidates were deliberate negative
  examples in the documentation-convention guide;
- current Markdown-sensitive inline formulas normalized;
- current/non-historical `text` diagrams contained no remaining diagrammatic
  ASCII-arrow residue;
- 131 broken historical file targets pinned to the original pre-refactor
  snapshot rather than modernized to successor documents.

The final historical-anchor review identified two links in:

```text
docs/history/reviews/external-dealii-boundary-evaluation/
    minimal-consumers-report.md
```

Those links named historical sections of the Step-4 integration report whose
headings were later reorganized. Their original referents were confirmed at the
audit baseline and the links were pinned to that frozen snapshot.

With those links repaired and the final documentation diff passing
`git diff --check`, no known human-readability defect remains that justifies
another broad pass.

## 30. Closure

The human-readability refactor did not attempt to make every document sound as
though it was written at project completion, make every large source file
small, or make every producer path pass through one universal object.

Instead, it made the repository easier to understand by clarifying **where
different kinds of truth live**.

A first-time reader now has a progressive manual rather than a roadmap-first
entry path. A user has task-oriented current references. A maintainer has
implementation maps separated from public programming guidance. Foundational
design and completed planning retain their historical perspective without being
mistaken for exact current implementation state. Source-derived studies and
replication evidence have their own home. Superseded implementation records and
reviews remain recoverable without competing with current authority.

The repository-structure interlude improved responsibility locality where the
evidence justified it and deliberately stopped where further splitting would
have been cosmetic, speculative, or disproportionately costly. Its result was
then reflected back into the documentation rather than allowing the
documentation to describe an obsolete physical tree.

The final state is not a repository with all history flattened into one
present-tense narrative. It is a repository with a clearer authority graph,
stable reading paths, documented provenance, and explicit reasons for the work
that was not done.

That satisfies the original audit's purpose: make the implemented architecture
inspectable to humans without redesigning it merely to make the presentation
simpler.
