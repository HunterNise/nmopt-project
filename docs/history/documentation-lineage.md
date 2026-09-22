# Documentation lineage

## Purpose

This page records the main documentation-structure transitions in the repository and
shows how to recover exact earlier versions from Git.

It is not a second archive of superseded Markdown files. Historical documents that
remain useful as project records live under `docs/history/`; earlier versions of
files that were moved, renamed, or later superseded remain available through Git.

The goal is to preserve the evolution of the documentation without maintaining a
parallel copy of every earlier tree.

This index intentionally describes **tracked repository history only**. Temporary
handoffs, scratch roadmaps, and other ignored working material are not part of the
documentation authority model.

## Main historical anchors

The most useful commits for documentation archaeology are:

```text
f53b7f009e5c    original human-readability audit baseline

bad5429ea437    docs(manual): organize overview and concept chapters
640d5a300319    docs(manual): reconcile links and conventions

997367f37b14    docs(structure): reorganize documentation by role

c4a4089fa3b0    docs(audit): close repository structure refactor
db4f2eab098a    docs(structure): reconcile documentation routing

16b9959524c1    docs(navigation): finalize public interface references
```

These commits mark different kinds of change:

- `bad5429ea437` organized the new explanatory material into the manual;
- `640d5a300319` is the last convenient snapshot before the large role-based
  documentation reorganization;
- `997367f37b14` reclassified existing documents by long-term role;
- `c4a4089fa3b0` closed the later repository-structure refactor that had paused
  the broad documentation rewrite while source paths and ownership stabilized;
- `db4f2eab098a` reconciled documentation routing against that settled source tree;
- `16b9959524c1` is the navigation/reference state after the public-reference
  rewrite and before the final implementation-documentation cleanup.

For exact history, prefer the full commit SHA printed by `git show` in the local
repository.

## Before the role-based documentation reorganization

At `640d5a300319`, the documentation tree still reflected the project's earlier
working categories.

The most important part of that layout was approximately:

```text
docs/
├── design/
│   ├── architecture.md
│   ├── composition-boundaries.md
│   ├── interface-specification.md
│   ├── pde-solver-boundary.md
│   ├── system-blueprint.md
│   └── theoretical-formalism.md
│
├── implementation/
│   ├── implementation-readiness-review.md
│   ├── v0/
│   │   ├── dealii-lowerer.md
│   │   └── executable-contract.md
│   └── v1/
│       └── semantic-compiler.md
│
├── decisions/
├── applications/
├── benchmarks/
├── case-studies/
├── guides/
├── manual/
├── planning/
└── reference/
```

This structure is historically useful because it shows how documents were understood
while the project was still being built:

- `design/` mixed long-lived architectural intent with implementation correspondence;
- `implementation/` contained both generation-specific implementation records and the
  then-current semantic compiler description;
- Chapter 5/6 application and benchmark material was split across several categories;
- reviews were still grouped below planning.

The later reorganization changed those roles without erasing the underlying history.

To inspect the complete tree:

```bash
git ls-tree -r --name-only 640d5a300319 -- docs/
```

## Manual organization before the role split

The explanatory documentation had already started converging on its current shape
before `997367f37b14`.

Commit `bad5429ea437` moved the architecture overview from:

```text
docs/overview/
```

to:

```text
docs/manual/overview/
```

and gave the concept chapters their numbered reading order:

```text
01-discrete-problem-anatomy.md
02-spaces-coordinates-and-duality.md
...
13-integrating-an-existing-pde-application.md
```

This was a pedagogical reorganization rather than an archival one. The manual was
becoming the first-time reading path while design, implementation, reference, and
review material retained different roles.

To inspect that move directly:

```bash
git show --stat bad5429ea437
git diff --name-status bad5429ea437^ bad5429ea437 -- docs/
```

## Role-based reorganization

Commit `997367f37b14` performed the major documentation reclassification.

Several path changes are especially useful when following older links or understanding
the origin of current documents.

### Design and implementation correspondence

```text
docs/design/theoretical-formalism.md
    → docs/design/mathematical-model.md

docs/design/system-blueprint.md
    → docs/internals/system-blueprint.md

docs/implementation/v1/semantic-compiler.md
    → docs/internals/compiler/semantic-compiler.md
```

The first move retained a long-lived mathematical role.

The latter two moves recognized that implementation correspondence and compiler
mechanics were not design authority merely because their original documents had been
written during the design/implementation phase.

### Superseded implementation generations

```text
docs/implementation/v0/executable-contract.md
    → docs/history/implementation/executable-contract.md

docs/implementation/v0/dealii-lowerer.md
    → docs/history/implementation/dealii-lowerer.md

docs/implementation/implementation-readiness-review.md
    → docs/history/implementation/implementation-readiness-review.md
```

These files were preserved as implementation-generation records rather than rewritten
to look current.

### Reviews and audits

```text
docs/planning/review/
    → docs/history/reviews/
```

Completed assessments, remediation reviews, and audit evidence moved out of mutable
planning while retaining their original contents and chronology.

### Accepted repository decisions

```text
docs/decisions/
    → docs/design/decisions/
```

Accepted decisions became part of the design record rather than a separate top-level
documentation category.

### Studies and source-derived application material

Material previously distributed across:

```text
docs/applications/
docs/benchmarks/
docs/case-studies/
docs/guides/
```

was reorganized under:

```text
docs/studies/
```

This grouped source Chapter 5/6 material, recipes, benchmark definitions,
replication records, and case studies without treating them as general framework
design or public API reference.

To inspect the complete migration:

```bash
git show --stat 997367f37b14
git diff --name-status 640d5a300319 997367f37b14 -- docs/
```

## Documentation pause during the repository-structure refactor

The role-based documentation reorganization intentionally preceded a separate
repository-structure readability refactor.

That source/test/application/tooling refactor used `997367f37b14` as its documentation
starting point. Broad reference and documentation rewriting was paused while source
paths, ownership boundaries, application seams, tests, and build organization were
reassessed.

The structure work later closed at:

```text
c4a4089fa3b0    docs(audit): close repository structure refactor
```

Its closure record is:

[`reviews/human-readability-audit/10-repository-structure-refactor-closure.md`](reviews/human-readability-audit/10-repository-structure-refactor-closure.md).

After that closure, the documentation branch resumed against the settled source tree.
Commit `db4f2eab098a` reconciled routing and paths before the task-oriented public
reference rewrite continued.

This distinction is useful when reading history:

```text
documentation information architecture
        ↓
role-based docs reorganization
        ↓
repository/source structure refactor
        ↓
documentation/reference work resumes
```

Changes made during the repository-structure branch should therefore not be mistaken
for documentation architecture proposals, and earlier documentation plans should not
be assumed to describe source paths that had not yet stabilized.

## Current implementation-documentation split

The final implementation-documentation pass replaced the two broad legacy
internals with narrower maintainer-oriented documents:

```text
docs/internals/system-blueprint.md
    → docs/internals/implementation-map.md

docs/internals/compiler/semantic-compiler.md
    → docs/reference/compiler.md
      + docs/internals/compiler.md

new:
docs/internals/runner.md
```

`system-blueprint.md` had mixed conceptual teaching, architecture
correspondence, implementation navigation, and verification guidance. Those
roles are now separated between the manual/design records and the current
implementation map.

`compiler/semantic-compiler.md` had combined public capability documentation
with compiler implementation mechanics. Current public compilation workflows,
supported realization families, products, diagnostics, and rejection
boundaries now belong to the task-oriented compiler reference. Request closure,
lowering/planning mechanics, target ownership, product packaging, provenance,
lifetimes, and verification belong to the compiler internals.

`runner.md` was added because the repository application has a distinct
maintainer-facing implementation structure that the public execution reference
intentionally abstracts: parameter/schema resolution, `RunSetPlan`, lifecycle
and manifest ownership, execution registration, and Chapter-6-specific
orchestration.

The superseded files are not copied under `history/`; Git remains their exact
archive and this page records their replacements.

## Recovering important earlier documents

Git is the canonical archive for files that were renamed or later superseded.

Examples:

```bash
# Foundational design-era system correspondence.
git show 640d5a300319:docs/design/system-blueprint.md

# Semantic compiler implementation document before role reclassification.
git show 640d5a300319:docs/implementation/v1/semantic-compiler.md

# Formal mathematical document before its final name.
git show 640d5a300319:docs/design/theoretical-formalism.md

# Earlier implementation-generation contract.
git show 640d5a300319:docs/implementation/v0/executable-contract.md
```

To follow a file across renames:

```bash
git log --follow -- docs/internals/system-blueprint.md
git log --follow -- docs/internals/compiler/semantic-compiler.md
git log --follow -- docs/design/mathematical-model.md
```

To compare the documentation tree at two points:

```bash
git diff --name-status 640d5a300319 997367f37b14 -- docs/
git diff --name-status 997367f37b14 16b9959524c1 -- docs/
```

## Preservation policy

A previous version does not need a copied file under `docs/history/` merely because
it is no longer current.

Prefer Git history when:

- the earlier file is simply a previous version of a still-recognizable document;
- an exact historical snapshot is more useful than a maintained duplicate;
- preserving the old file in the current tree would create competing documentation
  authority or broken historical links.

Keep a tracked file under `docs/history/` when it has durable value as a named project
record in its own right, such as:

- a completed implementation-generation record;
- an audit or review;
- a closure report;
- a historically meaningful draft intentionally retained by the project.

When a current document is superseded rather than merely edited, record its old path
and replacement in this lineage page as part of the same documentation change.

## Current authority

This lineage page explains where documents came from. It does not make an older
version authoritative.

For current use:

- start with the [manual overview](../manual/overview/README.md) for project
  orientation;
- use [`docs/reference/`](../reference/) for current public interfaces and
  operational workflows;
- use [`docs/internals/`](../internals/) for current implementation mechanics;
- use [`docs/design/`](../design/) for mathematical conventions, foundational
  architecture, and accepted decisions; and
- use [`docs/history/`](./) when the historical record itself is the subject of
  investigation.

If substantial development resumes, compare foundational design records with the
current implementation before treating early future-oriented statements as current
implementation commitments.
