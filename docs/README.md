# Documentation map

The documentation is organized by role and authority. Start with the smallest
document that answers the task, then follow its links to deeper material. Do
not read every document by default.

## Organization and authority

### Design

The documents under `design/` describe long-lived mathematical, semantic, and
architectural boundaries.

- [Interface specification](design/interface-specification.md) is the
  normative, implementation-neutral component and protocol contract.
- [Architecture record](design/architecture.md) explains the long-term design
  rationale.
- [Composition boundaries](design/composition-boundaries.md) summarizes
  component ownership and cross-cutting seams.
- [System blueprint](design/system-blueprint.md) is the shortest end-to-end
  mental model and code correspondence.
- [Theoretical formalism](design/theoretical-formalism.md) records the
  mathematical model and strong-to-variational bridge.

### Public reference

- [Public API reference](reference/README.md) is the entry point for exact
  `ProblemSpec`, compiler, solver, and experiment configuration contracts.

### Applications and benchmarks

- [Application recipes](applications/README.md) describes how concrete Chapter
  5/6 applications are assembled from the public API.
- [Benchmark specifications](benchmarks/README.md) describes frozen Chapter 6
  experiment contracts.

### Decisions

- [Repository organization](decisions/001-repository-organization.md) records
  the directory and authority boundaries for the application work.

### Implemented generations and selected policies

The documents under `implementation/` describe concrete implemented contracts,
realizations, and the policies selected for the first executable generations.

- [Implementation-readiness review](implementation/implementation-readiness-review.md)
  selects required defaults where the normative specification permits several
  policies. It is not a second capability ledger.
- [V0 executable contract](implementation/v0/executable-contract.md) and
  [deal.II v0 lowerer](implementation/v0/dealii-lowerer.md) define the direct
  reference slice and its explicit exclusions.
- [V1 semantic graph and compiler](implementation/v1/semantic-compiler.md)
  owns the exact implemented v1 capability and exclusion record.

### Guides and case studies

The documents under `guides/` explain how to implement or reproduce the
bounded Chapter 5/6 work. Their location does not imply that every catalogue
entry will be implemented.

- [Chapter 5 elliptic-control guide](guides/chapter-5-elliptic-control.md)
- [Chapter 6 numerical-methods guide](guides/chapter-6-numerical-methods.md)
- [Chapter 6 numerical-examples reference](guides/chapter-6-numerical-examples.md)

The documents under `case-studies/` are worked examples used to derive and
stress the interfaces.

- [Laplace composition-growth study](case-studies/laplace-growth.md)
- [Laplace formula and interface deltas](case-studies/laplace-interface-formulas.md)

### Planning and review evidence

- [Implementation roadmap](planning/implementation-roadmap.md) owns mutable
  compiler, solver, lowering, backend, and implementation status.
- [Application roadmap](planning/application-roadmap.md) owns mutable
  application-layer status, runner/artifact work, visualization, and B0–B2
  execution handoffs.
- [Chapter 5 problem library roadmap](planning/chapter-5-problem-library-roadmap.md)
  owns reusable, parameterized standard-problem recipes; feature status
  remains in the implementation roadmap.
- [Chapter 6 benchmark suite roadmap](planning/chapter-6-benchmark-suite-roadmap.md)
  owns frozen numerical examples, reproduction order, benchmark dependencies,
  and system-level acceptance gates.
- [Chapter 5 reviews](planning/review/chapter-5/README.md) index the C1/C2
  preparation and selected Chapter 5 implementation remediation reviews. They
  give the bounded repair sequences; the implementation roadmap owns
  remediation status.
- [Chapter 6 reviews](planning/review/chapter-6/README.md) index the incremental
  P6 implementation review and remediation handoffs; the implementation
  roadmap remains the mutable status ledger.
- [Pre-Chapter 5/6 review routing](planning/review/pre-ch5-ch6/README.md) tells
  an agent what to read for one bounded review batch.
- [Stage B roadmap](planning/review/pre-ch5-ch6/stage-b-roadmap.md) defines the accepted
  batch boundaries and gates without duplicating current project status.
- [Pre-Chapter 5/6 assessment](planning/review/pre-ch5-ch6/assessment.md) is
  the exhaustive evidence archive. Read only the assigned findings unless a
  tradeoff or scope decision requires wider context.
- [Assessment plan](planning/review/pre-ch5-ch6/assessment-plan.md) records the
  completed Stage A audit method and is not normal Stage B reading.

Project-wide working conventions live in the
[project conventions](../conventions/README.md). Read the applicable
convention before inspecting or changing repository content.

## Choose by task

| Audience or task | Start with | Then consult |
| --- | --- | --- |
| New contributor or agent | [System blueprint](design/system-blueprint.md) | Root `AGENTS.md`, the [implementation roadmap](planning/implementation-roadmap.md), and the [interface specification](design/interface-specification.md) |
| Author a Chapter 5/6 application | [Application recipes](applications/README.md) | [Public API reference](reference/README.md), the relevant Chapter guide, and the [v1 compiler](implementation/v1/semantic-compiler.md) |
| Add or reproduce a Chapter 6 benchmark | [Benchmark specifications](benchmarks/README.md) | [Chapter 6 benchmark roadmap](planning/chapter-6-benchmark-suite-roadmap.md), [numerical examples](guides/chapter-6-numerical-examples.md), and [Chapter 6 methods](guides/chapter-6-numerical-methods.md) |
| Understand the whole system and its code correspondence | [System blueprint](design/system-blueprint.md) | [Interface specification](design/interface-specification.md) and [v0 executable contract](implementation/v0/executable-contract.md) |
| Change semantic interfaces or ports | [Interface specification](design/interface-specification.md) | [Architecture](design/architecture.md) and [selected policies](implementation/implementation-readiness-review.md) |
| Change the backend-neutral executable API | [V0 executable contract](implementation/v0/executable-contract.md) | [Selected policies](implementation/implementation-readiness-review.md) and [roadmap](planning/implementation-roadmap.md) |
| Change deal.II code or compiler/lowering | [deal.II v0 lowerer](implementation/v0/dealii-lowerer.md) | [V0 executable contract](implementation/v0/executable-contract.md), [v1 compiler](implementation/v1/semantic-compiler.md), and [roadmap](planning/implementation-roadmap.md) |
| Change the v1 semantic graph/compiler | [V1 compiler](implementation/v1/semantic-compiler.md) | [Interface specification](design/interface-specification.md) and [selected policies](implementation/implementation-readiness-review.md) |
| Decide component ownership | [Composition boundaries](design/composition-boundaries.md) | [Architecture](design/architecture.md) and [interface specification](design/interface-specification.md) |
| Check mathematical signs or formulas | [Theoretical formalism](design/theoretical-formalism.md) | [Laplace formulas](case-studies/laplace-interface-formulas.md) and [growth study](case-studies/laplace-growth.md) |
| Implement a selected Chapter 5 application | [Chapter 5 guide](guides/chapter-5-elliptic-control.md) | [V1 compiler](implementation/v1/semantic-compiler.md) and [roadmap](planning/implementation-roadmap.md) |
| Implement selected Chapter 6 methods | [Chapter 6 methods](guides/chapter-6-numerical-methods.md) | [V0 contract](implementation/v0/executable-contract.md), [v1 compiler](implementation/v1/semantic-compiler.md), and [roadmap](planning/implementation-roadmap.md) |
| Explore standard Chapter 5 problems | [Chapter 5 problem library](planning/chapter-5-problem-library-roadmap.md) | [Chapter 5 guide](guides/chapter-5-elliptic-control.md), [v1 compiler](implementation/v1/semantic-compiler.md), and [implementation roadmap](planning/implementation-roadmap.md) |
| Reproduce Chapter 6 examples | [Chapter 6 benchmark suite](planning/chapter-6-benchmark-suite-roadmap.md) | [Numerical examples](guides/chapter-6-numerical-examples.md) and [Chapter 6 methods](guides/chapter-6-numerical-methods.md) |
| Plan application-layer work | [Application roadmap](planning/application-roadmap.md) | [Application recipes](applications/README.md), [benchmark specifications](benchmarks/README.md), and the [public application API](reference/application-api.md) |
| Repair the reviewed C1/C2 preparation | [Chapter 5 reviews](planning/review/chapter-5/README.md) | [Stage B roadmap](planning/review/pre-ch5-ch6/stage-b-roadmap.md), [pre-Chapter 5/6 assessment](planning/review/pre-ch5-ch6/assessment.md), and [v1 compiler](implementation/v1/semantic-compiler.md) |
| Repair the reviewed P5.1 implementation | [P5.1 remediation review](planning/review/chapter-5/p5.1-remediation-review.md) | [Implementation roadmap](planning/implementation-roadmap.md), [interface specification](design/interface-specification.md), and [v1 compiler](implementation/v1/semantic-compiler.md) |
| Repair the reviewed P5.2 implementation | [P5.2 remediation review](planning/review/chapter-5/p5.2-remediation-review.md) | [P5.1 remediation review](planning/review/chapter-5/p5.1-remediation-review.md), [implementation roadmap](planning/implementation-roadmap.md), and [selected policies](implementation/implementation-readiness-review.md) |
| Repair the reviewed P5.3 implementation | [P5.3 remediation review](planning/review/chapter-5/p5.3-remediation-review.md) | [Implementation roadmap](planning/implementation-roadmap.md), [interface specification](design/interface-specification.md), and [v1 compiler](implementation/v1/semantic-compiler.md) |
| Repair the reviewed P5.4 implementation | [P5.4 remediation review](planning/review/chapter-5/p5.4-remediation-review.md) | [P5.3 remediation review](planning/review/chapter-5/p5.3-remediation-review.md), [implementation roadmap](planning/implementation-roadmap.md), and [v1 compiler](implementation/v1/semantic-compiler.md) |
| Review the Chapter 6 implementation | [Chapter 6 reviews](planning/review/chapter-6/README.md) | [Chapter 6 methods](guides/chapter-6-numerical-methods.md), [implementation roadmap](planning/implementation-roadmap.md), and the relevant executable/compiler contract |
| Execute one Stage B review batch | [Pre-Chapter 5/6 review routing](planning/review/pre-ch5-ch6/README.md) | The current batch in the [Stage B roadmap](planning/review/pre-ch5-ch6/stage-b-roadmap.md), assigned findings, and task-specific authorities above |
| Review the complete pre-Chapter 5/6 evidence | [Pre-Chapter 5/6 assessment](planning/review/pre-ch5-ch6/assessment.md) | [Assessment plan](planning/review/pre-ch5-ch6/assessment-plan.md) only when the audit method matters |
| Choose the next implementation task | [Implementation roadmap](planning/implementation-roadmap.md) | The task-specific compiler, solver, backend, or semantic contract listed above |
| Edit Markdown or LaTeX | [Documentation conventions](../conventions/documentation.md) | The document being changed |
