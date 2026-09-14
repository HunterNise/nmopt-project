# nmopt Human-Readability and Architecture Audit Plan

## Purpose

Reconstruct the current `nmopt-project` architecture from the implementation itself, then assess how well the repository communicates that architecture to human readers, users, reviewers, and maintainers.

Existing documentation is treated as evidence to reconcile, not as an authoritative description of the current code unless verified against the implementation.

The audit is intended to determine:

- the actual scope of the project;
- the major architectural layers and their responsibilities;
- the dependency and ownership boundaries between those layers;
- the public surfaces intended for users versus project-internal machinery;
- which parts of the repository are essential library code, optional infrastructure, experiments, validation machinery, or historical material;
- how the project should be explained to humans;
- what documentation structure best supports that explanation;
- whether documentation changes are sufficient, or whether targeted source-level cleanup is justified.

## Primary audience for the final repository

The main human reader is assumed to be technically competent in numerical PDEs and C++, possibly familiar with deal.II, but unfamiliar with the history and internal terminology of this repository.

Secondary audiences are:

1. a technical reviewer or evaluator assessing the project;
2. a developer trying to modify or extend one subsystem;
3. a user trying to consume the numerical library or integrate an existing application;
4. the project author preparing the final presentation and future maintenance.

The repository should not require knowledge of the Codex development history or `.agents/` instructions to form a correct mental model.

## Audit principles

1. **Code first.** Public headers, implementations, tests, build dependencies, and runtime paths are the primary source of truth.
2. **Responsibilities before directories.** Architectural layers are inferred from ownership, dependencies, and behavior rather than copied from namespace or directory names.
3. **Current behavior versus intent.** Implemented architecture, accepted design intent, future targets, and historical experiments must remain distinct.
4. **Minimal human model.** The final explanatory architecture should use the fewest conceptual layers that accurately describe the code.
5. **No premature refactor.** Readability problems are recorded before deciding whether they require prose, comments, naming/file cleanup, or code changes.
6. **Preserve evidence.** Historical reviews and experiments may remain valuable, but they should not compete with current documentation for authority.
7. **Traceability.** Important conclusions should identify the code, tests, build targets, or documents that support them.

## Audit outputs

The audit will be recorded as several Markdown files rather than one monolithic report.

### `00-audit-index.md`

Living index and status page.

Contains:

- audit purpose and scope;
- inspected areas;
- files produced;
- unresolved questions;
- current high-level findings;
- progress through the audit.

### `01-repository-scope.md`

Answers: **What is this project now?**

Examines:

- core numerical library;
- semantic/compiler subsystem;
- deal.II numerical realization;
- formulations and optimization;
- application/runner layer;
- Chapter 5/6 recipes and benchmarks;
- external integration work;
- verification/reference infrastructure;
- tooling and generated artifacts.

Classifies areas as appropriate into:

- core library capability;
- optional user-facing subsystem;
- project/application infrastructure;
- verification/evidence machinery;
- historical or research artifact.

### `02-architecture-map.md`

Code-derived architectural synthesis.

For each inferred layer:

- purpose;
- owned state and responsibilities;
- inputs and outputs;
- dependencies;
- public interfaces;
- representative source files;
- representative tests;
- lifetime/ownership rules;
- what the layer deliberately does not own.

Includes dependency and runtime-flow diagrams.

This file is descriptive first: it records what the implementation actually does before recommending a documentation structure.

### `03-runtime-paths.md`

Follows a small number of representative end-to-end executions through the code.

Likely paths:

1. semantic `ProblemSpec` → compiler → reduced formulation → optimizer;
2. external deal.II application → callback binding → reduced formulation → optimizer;
3. Chapter 5/6 scenario → runner → compiler/execution adapter → artifacts;
4. one reduced objective/gradient evaluation;
5. one selected non-reduced or specialized path if it materially changes the architecture.

For each path:

- caller/callee sequence;
- state ownership;
- important type transitions;
- mathematical operation performed;
- optional versus mandatory infrastructure.

### `04-public-api-and-user-surfaces.md`

Examines what a human user actually has to interact with.

Separates:

- external-application integration surface;
- semantic/compiler authoring surface;
- formulation/solver surface;
- application/runner surface;
- lower-level expert interfaces.

Assesses:

- conceptual load;
- naming;
- required boilerplate;
- discoverability;
- lifetime requirements;
- representation ceremony;
- places where public API and project-internal machinery are difficult to distinguish.

This is an ergonomics audit, not automatically a proposal to redesign the API.

### `05-source-readability.md`

Audits the code itself as documentation.

Looks for:

- public types with unclear purpose;
- ownership/lifetime assumptions not visible locally;
- mathematically important conventions missing near the implementation;
- misleading or overly generic names;
- large files mixing responsibilities;
- unnecessary indirection;
- comments that narrate syntax instead of rationale;
- missing comments where design constraints are non-obvious;
- stale comments;
- locations where Doxygen/doc comments would materially improve usability.

Findings are classified as:

- documentation-only;
- source comment/docstring improvement;
- naming/file organization cleanup;
- possible non-behavioral refactor;
- possible architectural/API issue.

### `06-documentation-audit.md`

Inventory and critique of the current documentation system.

Each substantial document or document family is classified by role:

- overview;
- design;
- reference;
- guide;
- case study;
- planning/status;
- review/audit;
- historical/superseded;
- draft;
- agent instruction.

For each, record:

- intended role;
- apparent authority;
- whether it matches current code;
- duplication;
- stale or future-looking claims;
- terminology mismatches;
- recommended disposition.

Possible dispositions:

- keep as current authority;
- edit;
- merge;
- split;
- move/demote;
- mark historical/superseded;
- archive;
- delete if genuinely redundant and valueless.

### `07-human-information-architecture.md`

Designs the human-facing documentation structure after the code and docs audits are complete.

Determines:

- canonical entry points;
- overview layer;
- reference/design split;
- task-oriented guides;
- case studies;
- historical/evidence routing;
- role of root `README.md`;
- role of `docs/README.md`;
- relationship between human docs and `.agents/`;
- navigation conventions;
- terminology policy.

The previously proposed `docs/overview/` structure is treated as a hypothesis and may be revised.

### `08-recommendations.md`

Final actionable proposal.

Separates recommendations into levels:

#### A. Documentation structure
New/edit/move/archive/delete operations.

#### B. Source documentation
Comments, docstrings, public-header explanations.

#### C. Low-risk readability cleanup
Naming, file moves, splitting mixed files, dead-code removal, similar changes that do not alter architecture.

#### D. API or architectural changes
Only if the audit finds that documentation cannot reasonably compensate for a structural problem.

Every recommendation includes:

- motivation;
- evidence;
- expected human benefit;
- risk;
- whether it is necessary before final presentation;
- proposed priority.

## Inspection strategy

The repository will be inspected by architectural concern rather than by simply reading directories alphabetically.

### Phase 1 — Repository topology and dependency skeleton

Inspect:

- root build files;
- CMake targets and interface libraries;
- public include tree;
- applications;
- tests;
- tools.

Goal:

- identify major compiled components and dependency directions;
- distinguish public headers from application/test-only code;
- identify likely architectural seams.

Output begins in `01-repository-scope.md` and `02-architecture-map.md`.

### Phase 2 — Mathematical and solver contracts

Inspect:

- layouts and primal/covector representation;
- executable model contract;
- metric/constraint contracts;
- solve reports/services;
- reduced DTO/formulation;
- optimization solvers;
- reference models and contract tests.

Questions:

- what is the true stable numerical boundary?
- what does the formulation actually require?
- which contracts are general and which reflect current implementation choices?
- where are mathematical conventions encoded?

### Phase 3 — Semantic and compiler path

Inspect:

- `ProblemSpec` and semantic types;
- validation/resolution;
- compiler plans/registrations;
- typed numerical realization;
- `CompiledProblemT`;
- native application views;
- relevant compiler tests.

Questions:

- what problem does the semantic layer solve?
- what does the compiler own?
- how generic is the current implementation?
- what remains target-specific?
- what information survives compilation and why?
- where does this path converge with the external path?

### Phase 4 — deal.II numerical realization

Inspect:

- backend policy;
- coordinate/reconstruction machinery;
- metrics;
- constraints;
- linear/KKT/PDAS services;
- selected compiler realizations;
- deal.II tests.

Questions:

- what does “backend” actually mean in this repository?
- which responsibilities belong to deal.II realization rather than backend algebra?
- which components are reusable numerical services?
- what ownership/lifetime assumptions exist?

### Phase 5 — External application path

Inspect:

- Step-4 source lineage;
- application adaptations;
- Problem A and Problem B mathematics;
- minimal bindings;
- native versus nmopt execution;
- tests/evaluation machinery.

Goal:

- use the external path as an independent check of the inferred solver boundary;
- separate intrinsic OCP work from nmopt-specific representation and wiring;
- identify which parts of the larger repository are truly optional.

### Phase 6 — Application and execution layer

Inspect:

- recipes/scenarios/catalog;
- application harness/runner;
- deal.II execution adapters;
- `nmopt-runner`;
- artifacts/manifests;
- parameter files;
- reporting/post-processing tools.

Questions:

- is this part of the reusable library, a client of the library, or both?
- which pieces are needed only for the research/benchmark project?
- what does a normal external user need to know about them?
- does the directory/API structure communicate that distinction?

### Phase 7 — Representative runtime traces

After individual layers are understood, trace the selected end-to-end paths and validate the proposed dependency diagram against actual calls and object lifetimes.

Output: `03-runtime-paths.md`.

### Phase 8 — Source readability pass

Revisit representative public headers and complex implementation points from a human-maintenance perspective.

Output: `05-source-readability.md`.

### Phase 9 — Documentation reconciliation

Only after the implementation model is stable:

- inspect design documents;
- reference documents;
- guides;
- roadmaps;
- reviews;
- case studies;
- drafts;
- READMEs;
- agent instructions.

Compare them against the code-derived architecture.

Output: `06-documentation-audit.md`.

### Phase 10 — Synthesis

Produce:

- final project scope;
- minimal architectural decomposition;
- human-facing information architecture;
- concrete recommendations;
- proposed implementation sequence.

Outputs:

- `07-human-information-architecture.md`;
- `08-recommendations.md`.

## Evidence standard

A finding should preferably be supported by at least one of:

- public header/API;
- implementation that establishes ownership or behavior;
- focused test that demonstrates the contract;
- build dependency or target boundary;
- actual application call path.

Existing documentation may corroborate the finding, but should not be its sole basis when the implementation can be inspected directly.

Important discrepancies between docs and code will be recorded explicitly rather than silently corrected in the audit narrative.

## Questions to keep open until evidence resolves them

- What is the smallest accurate definition of `nmopt`?
- Are the semantic/compiler facilities part of the core identity or one optional front end?
- Where should the boundary between “library” and “research application” be drawn?
- Is `application/` genuinely a library layer or primarily project orchestration?
- Does “backend” currently name one coherent responsibility?
- Are `contract/`, formulation, and solver best explained as separate human-facing layers?
- Which abstractions are stable public concepts versus artifacts of the current implementation?
- Which historical documents still have enduring explanatory value?
- Is the existing directory structure suitable for humans once documented?
- Are there source-level readability problems severe enough to justify code movement or cleanup?
- Does any public API problem rise above ergonomics and warrant architectural change?

## Stopping condition

The audit is complete when we can answer, from implementation evidence:

1. what `nmopt` is;
2. what its major layers are and why they exist;
3. how the layers communicate;
4. which paths a user may enter through;
5. which project subsystems are optional;
6. which documents should be authoritative for humans;
7. what source-level documentation is missing;
8. whether code changes are required for human maintainability;
9. what exact final refactor should be performed.

The subsequent refactor should not begin until these answers are stable enough that we are no longer designing the documentation structure around guesses.
