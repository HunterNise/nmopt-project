# Documentation audit

**Audit baseline:** [`f53b7f009e5c`](https://github.com/HunterNise/nmopt-project/commit/f53b7f009e5c418ec4f3855db29f7eb924faf6f1)  
**Status:** complete for the documentation-refactor scope: architecture, implementation authority, public API, application execution, project status, and repository-organization documents were reconciled to the audit baseline.

## 1. Overall finding

The repository does **not** primarily suffer from a lack of documentation. It has the opposite problem:

> there is enough documentation to preserve design intent, implementation history, review evidence, exact APIs, benchmark policy, and agent workflow, but the default navigation does not sufficiently separate those roles for a human approaching the finished project.

Many individual documents are strong. The problem is the authority and routing graph among them.

A new reader currently encounters several kinds of statements before knowing which kind they are reading:

- current implemented architecture;
- normative target architecture;
- implementation-generation records;
- future or deferred architecture;
- mutable roadmap status;
- closed review evidence;
- benchmark/source-study policy;
- coding-agent workflow.

The documentation refactor should therefore **add a small inspectability layer and demote historical/planning material from the default path**, rather than rewrite everything.

## 2. Current documentation structure

At the baseline, `docs/` contains:

```text
README.md
applications/
benchmarks/
case-studies/
decisions/
design/
drafts/
guides/
implementation/
planning/
reference/
```

There is no `docs/overview/`.

The role taxonomy itself is mostly sensible. The problem is that `docs/README.md` is a 17 KB catalogue whose default new-contributor route enters `design/system-blueprint.md`, a mutable implementation roadmap, and the normative interface specification.

That is a good **agent implementation** route. It is not the best **finished-project understanding** route.

## 3. Root README

Source: [`README.md`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/README.md)

### What works

- build/test instructions are concrete;
- deal.II is correctly optional for neutral profiles;
- application execution has a runnable entry point;
- repository layout is stated;
- documentation and agent instructions are linked.

### Problems for a finished project

The scope section says:

> “The current implementation establishes the first scalar finite-element slice”

and then says:

> “The application and benchmark layer is still being extended.”

That framing understates the current implementation and presents the repository as an active construction sequence. The current code has:

- multiple reduced-search policies and trust region;
- supplied OTD;
- quadratic KKT;
- complementarity/PDAS;
- many registered deal.II targets;
- B1 reproduction evidence;
- B2 completed negative replication evidence;
- a closed external-integration evaluation.

The application roadmap still lists later work, but that is different from saying the current repository is merely “the first slice.”

The “Where to start” section points humans directly to:

- `system-blueprint.md`;
- `architecture.md`;
- `interface-specification.md`;
- application assembly/execution references.

This skips a current, compact project overview because none exists.

### Recommended disposition

**Rewrite the scope and “Where to start” sections.** Keep build/run instructions. Route first to a new `docs/overview/project-architecture.md`, then to path-specific overview/reference material.

## 4. `docs/README.md`

Source: [`docs/README.md`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/docs/README.md)

### What works

The document has unusually good awareness of authority. It explicitly distinguishes:

- design;
- reference;
- implementation records;
- guides/case studies;
- planning/review.

That taxonomy should be preserved.

### Main problem

It is still primarily an **agent-oriented catalogue**.

Examples:

- “New contributor or agent” routes to `system-blueprint -> implementation roadmap -> interface specification`;
- many task rows are “repair the reviewed P5.x implementation,” “execute Stage B,” or “choose the next implementation task”;
- closed review/planning records occupy substantial visual weight.

For a frozen or presentation-ready project, normal reading should not begin in remediation history.

### Concrete stale text

The Public reference section says the parameter-file reference defines the **proposed** `.prm` and JSON schemas.

But [`reference/parameter-files.md`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/docs/reference/parameter-files.md) explicitly states:

> **Status: implemented for the registered B1/B2 Chapter 6 runner slice.**

### Recommended disposition

Turn `docs/README.md` into a short router with four primary intentions:

```text
Understand the project
Use the library
Reproduce project applications/experiments
Read design/history/evidence
```

Retain a compact “for coding agents” link to `.agents/`, but move implementation/remediation routing out of the default human path.

## 5. Design documents

Current files:

```text
architecture.md
composition-boundaries.md
interface-specification.md
pde-solver-boundary.md
system-blueprint.md
theoretical-formalism.md
```

### 5.1 `theoretical-formalism.md`

Source: [`docs/design/theoretical-formalism.md`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/docs/design/theoretical-formalism.md)

This is a strong long-lived mathematical document. It clearly states:

- the abstract constrained problem;
- variables/data/equation blocks;
- pairings;
- observations/losses;
- derivatives/adjoints;
- the strong-to-variational boundary;
- analysis assumptions left to the user.

**Disposition: keep as design/theory authority.**

It should be linked from overviews, not used as the initial project introduction.

### 5.2 `pde-solver-boundary.md`

Source: [`docs/design/pde-solver-boundary.md`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/docs/design/pde-solver-boundary.md)

This is one of the strongest current architecture documents.

It correctly records:

```text
ProblemSpec -> compiler/lowering --------\
                                          +-> mathematical operations -> formulation -> optimizer
existing application -> adapter --------/
```

It also states the key implemented boundary:

> type erasure is a useful solver boundary, not the ownership boundary for all compiled numerical state.

It distinguishes semantic meaning, numerical realization, state coordinates, metric, regularization, Hessian, and application output.

**Disposition: keep as current design authority.**

It should be the principal deep-design link from the new architecture overview.

### 5.3 `interface-specification.md`

Source: [`docs/design/interface-specification.md`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/docs/design/interface-specification.md)

The document is intentionally normative and implementation-neutral. That role is useful.

However, some sections leak current implementation/history into the normative language, for example references to:

- “the existing P5.1 `robin_source` component”;
- Chapter 6 transport forms;
- specific initial implementation relationships.

This weakens the document's claim to be implementation-neutral.

**Disposition: keep normative, but remove or clearly mark implementation examples/history.**

Current implementation correspondence belongs in the v1 capability record or a separate note.

### 5.4 `architecture.md`

Source: [`docs/design/architecture.md`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/docs/design/architecture.md)

This is a valuable design-direction record, but it is not a reliable current-project introduction.

It contains future-facing language such as:

- “the first working version must be deliberately small”;
- “Initial vertical slice” with an implementation checklist;
- a target composition architecture broader than the current bounded compiler.

**Disposition: keep as architectural intent/design record, but stop calling it the default ‘architecture overview.’**

Add an explicit status note such as:

> This document records architectural intent and target invariants. For the implemented architecture, start with `docs/overview/project-architecture.md`.

### 5.5 `composition-boundaries.md`

Source: [`docs/design/composition-boundaries.md`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/docs/design/composition-boundaries.md)

This is a useful compact design principle document. It explains local component ownership and cross-cutting mathematical cases well.

Its language is intentionally prescriptive (“should be”, “decisions to make now”), so it belongs under design rather than overview.

**Disposition: keep; demote from first-read status.**

### 5.6 `system-blueprint.md`

Source: [`docs/design/system-blueprint.md`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/docs/design/system-blueprint.md)

This document currently tries to serve several roles:

- shortest mental model;
- current-vs-target architecture;
- mathematical primer;
- executable contract explanation;
- DTO derivation;
- v1 realization walkthrough;
- source-code tour;
- authority router.

It is useful, but that breadth is exactly why it is not the ideal canonical overview.

It also explicitly contains:

```text
TARGET COMPONENT-LOWERING ARCHITECTURE (NOT YET IMPLEMENTED)
```

beside current implementation diagrams.

**Disposition: retain as a deep implementation blueprint, or eventually split/trim it after the new overview layer exists.**

At minimum, remove “shortest path to a working mental model” language once a true current overview exists.

## 6. Implementation records

### 6.1 `implementation/v1/semantic-compiler.md`

Source: [`docs/implementation/v1/semantic-compiler.md`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/docs/implementation/v1/semantic-compiler.md)

This is a strong exact registered-capability ledger.

It accurately states the hybrid current compiler structure:

```text
semantic graph
 -> validator/resolver/compiler diagnostics
 -> bounded component plan OR registered target strategy
 -> compiled product
```

It also records explicit exclusions and realization details.

**Disposition: keep as exact implementation/capability authority.**

Do not duplicate its large capability table into overview docs.

### 6.2 `implementation/v0/dealii-lowerer.md`

Source: [`docs/implementation/v0/dealii-lowerer.md`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/docs/implementation/v0/dealii-lowerer.md)

This is already clearly marked as a historical record and routes readers to v1/current integration.

**Disposition: keep historical; remove from normal “where to start” routes.**

### 6.3 `implementation/v0/executable-contract.md`

Source: [`docs/implementation/v0/executable-contract.md`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/docs/implementation/v0/executable-contract.md)

This file has become conceptually awkward.

Despite the `v0` name, current docs route to it as the exact solver-facing API. The file now describes a large current surface:

- backend-parametric executable model;
- metrics/constraints;
- reduced DTO;
- many direction policies;
- several line searches;
- trust region;
- explicit Hessians;
- current deal.II backend notes.

Meanwhile, separate current contracts exist for KKT/PDAS/supplied OTD.

### Recommended disposition

Create a **current** `docs/reference/numerical-contracts.md` that explains the presently shipped contract/formulation/solver surface.

Then reduce `implementation/v0/executable-contract.md` to its true historical role, or leave it intact with a strong banner pointing to the current reference.

This avoids silently rewriting history while removing “v0” from the current user's critical path.

### 6.4 `implementation-readiness-review.md`

Source: [`docs/implementation/implementation-readiness-review.md`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/docs/implementation/implementation-readiness-review.md)

This is a design-review outcome and default-policy record for the early implementation generation.

It now explicitly notes that the original direct scalar realization was retired.

It also contains target/default mechanisms that are not identical to the current bounded implementation, such as a central typed DAG as a “first default.”

**Disposition: preserve as design history/policy rationale, not current implementation reference.**

`docs/README.md` should stop describing it simply as “selected policies” without historical qualification.

## 7. Public references

Current files:

```text
application-api.md
application-execution.md
external-dealii-solver-integration.md
parameter-files.md
```

These are generally strong and should remain exact-detail destinations.

### 7.1 External integration reference

Source: [`docs/reference/external-dealii-solver-integration.md`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/docs/reference/external-dealii-solver-integration.md)

This is exemplary reference documentation:

- exact public types;
- precise ownership;
- explicit tested limitations;
- no requirement to adopt semantic/compiler/application infrastructure;
- clear note that residual/JVP are required even though unused by the selected reduced path;
- clear lifetime rules.

**Disposition: keep.**

### 7.2 Application API reference

Source: [`docs/reference/application-api.md`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/docs/reference/application-api.md)

Technically useful, but it calls itself **agent-facing**.

It also presents one large “application assembly” path that combines recipe, semantic graph, compiler, optimization, and provenance. That is correct for the semantic application path but not universal project usage.

**Disposition: keep as semantic/compiler application reference, but make wording human-neutral and clarify its scope in the title/introduction.**

Suggested conceptual title:

> Semantic/compiler application assembly reference

### 7.3 Application execution reference

Source: [`docs/reference/application-execution.md`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/docs/reference/application-execution.md)

This is an exact repository execution/artifact reference and correctly says the runner does not lower PDEs or solve optimization itself.

**Disposition: keep, but route it under “reproduce project experiments,” not “understand/use core library.”**

### 7.4 Parameter files

Source: [`docs/reference/parameter-files.md`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/docs/reference/parameter-files.md)

This is current and explicit about implemented B1/B2 scope.

**Disposition: keep.**

## 8. Planning/status documents

The current planning directory contains several large mutable or historical records, including:

- 85 KB `implementation-roadmap.md`;
- 81 KB external boundary evaluation;
- 27 KB application roadmap;
- Chapter 5/6 roadmaps;
- many review archives.

These are valuable evidence, but should not dominate navigation once the project is frozen.

### 8.1 Major status contradiction

[`planning/application-roadmap.md`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/docs/planning/application-roadmap.md) says:

- B1 is **reproduction-verified**;
- B2 is **framework-verified; replication attempt closed negatively**;
- parameter files are implemented for B1/B2.

By contrast, [`planning/implementation-roadmap.md`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/docs/planning/implementation-roadmap.md) still says in its “Current handoff state” that:

- B1's selected source-oriented profile still needs its authoritative release run;
- the B2 benchmark activation gate remains open pending refreshed runs/investigation.

These cannot both be the current status.

The root README currently tells readers that application/benchmark status is tracked in **both** roadmaps, so this inconsistency is exposed to ordinary readers.

### Recommended disposition

If development is frozen:

1. stop treating the implementation roadmap as a current handoff authority;
2. add a top banner identifying it as a historical implementation ledger;
3. record one short final/current scope statement elsewhere;
4. route humans away from both roadmaps unless they are studying development history.

If continued development resumes later, first reconcile the two ledgers.

## 9. Repository-organization decision

Source: [`docs/decisions/repository-organization.md`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/docs/decisions/repository-organization.md)

This accepted decision is partly stale.

### Concrete mismatch 1 — `parameters/`

The document later calls `parameters/` a **reserved future location**, but the directory exists and the parameter-file reference declares B1/B2 implementation complete.

### Concrete mismatch 2 — test tree

The decision's top-level layout lists:

```text
tests/
  benchmark/
```

The current tree has:

```text
tests/
  application/
  contract/
  dealii/
  semantic/
  support/
  tools/
```

No top-level `tests/benchmark/` exists.

### Recommended disposition

Update the accepted repository decision to describe the repository that actually exists.

Do not preserve stale “future” language in an accepted current organization decision.

## 10. External Step-4 documentation

The Step-4 local documentation is a good model for the global refactor.

### `external-integration-overview.md`

Source: [`apps/external-dealii/step-4/external-integration-overview.md`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/apps/external-dealii/step-4/external-integration-overview.md)

It does one job well:

> explain the conceptual/runtime path.

It uses diagrams, mathematical operations, ownership, one evaluation, and A/B comparison.

### `integration-report.md`

Owns source accounting, exact changes, measurements and reproduction.

### Closure review

Owns bounded evidence and unsupported claims.

### Lesson

The project-wide docs should copy this **division of document personality**:

```text
overview       understanding
reference      exact API
report/review  evidence
design         decisions/invariants
history        why/how the implementation got here
```

## 11. Agent-facing language in human docs

Several current public docs use phrases such as:

- “agent-facing reference”;
- “an agent selects parameters”;
- “what agents need to consume”;
- remediation routes by RF/P unit.

Those statements are useful for autonomous development, but they should not define the public conceptual surface.

### Recommended boundary

```text
docs/          human-readable project/design/reference documentation
.agents/       agent workflow and repository-changing instructions
```

Human docs may still be usable by agents. They should not require agent context to make sense.

## 12. Recommended authority hierarchy

For a frozen, inspectable project:

```text
CURRENT UNDERSTANDING
docs/overview/*

LONG-LIVED DESIGN / NORMATIVE MODEL
docs/design/*

EXACT IMPLEMENTED API / SCHEMAS
docs/reference/*
docs/implementation/v1/semantic-compiler.md   # exact compiler capability ledger

USAGE / REPRODUCTION
docs/guides/*
docs/applications/*
docs/benchmarks/*

EVIDENCE / HISTORY
apps/.../integration-report.md
docs/planning/review/*
docs/planning/*roadmap*.md
historical implementation records
```

The critical change is that **current understanding must no longer be reconstructed from design + roadmaps + implementation history**.
