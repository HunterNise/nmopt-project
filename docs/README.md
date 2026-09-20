# Documentation map

The documentation is organized by role and authority. Start with the smallest
document that answers the task, then follow its links to deeper material. Do
not read every document by default.

## Start here

For a first-time reading of the project, start with the
[manual overview](manual/overview/README.md), then follow the
[concept chapters](manual/concepts/README.md). Use the sections below when you
need an exact design decision, public contract, study record, active plan, or
historical review.

## Organization and authority

```text
docs/
  manual/       teach the current project
  design/       long-lived architecture, mathematics, and accepted decisions
  reference/    exact current public interfaces, configuration, and execution
  internals/    current implementation mechanics and capability internals
  studies/      Chapter 5/6 source, application, benchmark, and case-study corpus
  planning/     genuinely active and mutable roadmaps
  history/      superseded implementations, reviews, audits, and evidence
```

### Manual

- [Manual overview](manual/overview/README.md) is the recommended first-time
  reading path.
- [Concept chapters](manual/concepts/README.md) develop the project’s
  mathematical, formulation, compiler, and integration language.

### Design

- [Architecture](design/architecture.md)
- [Composition boundaries](design/composition-boundaries.md)
- [Interface specification](design/interface-specification.md)
- [Mathematical model](design/mathematical-model.md)
- [PDE–solver boundary](design/pde-solver-boundary.md)
- [Repository organization decision](design/decisions/repository-organization.md)
- [Parameter and plotting profiles decision](design/decisions/parameter-and-plotting-profiles.md)

These documents own long-lived architecture, mathematical conventions,
component contracts, and accepted design decisions.

### Reference

- [Application assembly API](reference/application-api.md)
- [Application execution](reference/application-execution.md)
- [External deal.II solver integration](reference/external-dealii-solver-integration.md)
- [Parameter files](reference/parameter-files.md)

These documents describe exact current public interfaces, schemas,
configuration, execution, and integration contracts.

### Internals

- [System blueprint](internals/system-blueprint.md)
- [Semantic compiler internals](internals/compiler/semantic-compiler.md)

These documents describe current implementation mechanics and capability
internals without replacing the public references or design contracts.

### Studies

- [Chapter 5 source catalogue](studies/chapter-5/source-catalogue.md)
- [Chapter 5 application recipes](studies/chapter-5/recipes.md)
- [Chapter 6 numerical methods](studies/chapter-6/numerical-methods.md)
- [Chapter 6 numerical examples](studies/chapter-6/numerical-examples.md)
- [Chapter 6 scenarios](studies/chapter-6/scenarios.md)
- [Chapter 6 benchmarks](studies/chapter-6/benchmarks.md)
- [B1 replication](studies/chapter-6/b1-replication.md)
- [B2 replication](studies/chapter-6/b2-replication.md)
- [Laplace growth case study](studies/case-studies/laplace-growth.md)
- [Laplace interface formulas](studies/case-studies/laplace-interface-formulas.md)

This is the current Chapter 5/6 source, application, benchmark, and case-study
corpus. It records study content and evidence rather than general framework
architecture.

### Planning

- [Implementation roadmap](planning/implementation-roadmap.md)
- [Application roadmap](planning/application-roadmap.md)
- [Chapter 5 problem-library roadmap](planning/chapter-5-problem-library-roadmap.md)
- [Chapter 6 benchmark-suite roadmap](planning/chapter-6-benchmark-suite-roadmap.md)

Planning contains genuinely active, mutable work order and status. Completed
reviews and audits belong under history, even when they remain useful context.

### History

- [`history/implementation/`](history/implementation/) contains superseded
  implementation records.
- [`history/reviews/`](history/reviews/) contains reviews, audits, closure
  reports, and historical evidence.
- [Review history index](history/reviews/README.md) routes the available review,
  audit, and closure series.
- [Human-readability audit](history/reviews/human-readability-audit/00-audit-index.md)
  is the index for that specific audit series.

The application-local [external Step-4 overview](../apps/external-dealii/step-4/external-integration-overview.md)
and [integration report](../apps/external-dealii/step-4/integration-report.md)
describe the external application study and its evidence.

Agent working instructions live in the
[agent instructions](../.agents/README.md). Read the applicable instruction
before inspecting or changing repository content.

## Choose by task

| Audience or task | Start with | Then consult |
| --- | --- | --- |
| First-time project reader | [Manual overview](manual/overview/README.md) | [Concept chapters](manual/concepts/README.md), then the reference for the task at hand |
| Understand the whole system and code correspondence | [System blueprint](internals/system-blueprint.md) | [Interface specification](design/interface-specification.md) and [semantic compiler internals](internals/compiler/semantic-compiler.md) |
| Author or modify a semantic problem | [Problem authoring](reference/problem-authoring.md) | [Compiler](reference/compiler.md) and the relevant manual concept chapter |
| Author a reusable nmopt-native application | [Application authoring](reference/application-authoring.md) | [Problem authoring](reference/problem-authoring.md), [compiler](reference/compiler.md), and [application execution](reference/application-execution.md) |
| Connect an existing PDE application | [External deal.II solver integration](reference/external-dealii-solver-integration.md) | [Optimization](reference/optimization.md) and the relevant public contract/deal.II headers |
| Add or use an optimization method | [Optimization](reference/optimization.md) | [Reduced optimization methods](manual/concepts/06-reduced-optimization-methods.md) and the public solver/contract headers |
| Add or reproduce a Chapter 6 benchmark | [Chapter 6 benchmarks](studies/chapter-6/benchmarks.md) | [Application execution](reference/application-execution.md), [parameter files](reference/parameter-files.md), and the [benchmark roadmap](planning/chapter-6-benchmark-suite-roadmap.md) |
| Change semantic interfaces or ports | [Interface specification](design/interface-specification.md) | [Problem authoring](reference/problem-authoring.md), [architecture](design/architecture.md), and [composition boundaries](design/composition-boundaries.md) |
| Change deal.II compiler/lowering implementation | [Compiler](reference/compiler.md) | [Semantic compiler internals](internals/compiler/semantic-compiler.md), [interface specification](design/interface-specification.md), and [implementation roadmap](planning/implementation-roadmap.md) |
| Check mathematical signs or formulas | [Mathematical model](design/mathematical-model.md) | [Laplace formulas](studies/case-studies/laplace-interface-formulas.md) and [growth study](studies/case-studies/laplace-growth.md) |
| Implement or reproduce a Chapter 5 application | [Chapter 5 recipes](studies/chapter-5/recipes.md) | [Application authoring](reference/application-authoring.md), [source catalogue](studies/chapter-5/source-catalogue.md), and [problem-library roadmap](planning/chapter-5-problem-library-roadmap.md) |
| Generate or inspect application runs | [Application execution](reference/application-execution.md) | [Parameter files](reference/parameter-files.md), [Chapter 6 benchmarks](studies/chapter-6/benchmarks.md), and [application roadmap](planning/application-roadmap.md) |
| Design parameter files or plotting profiles | [Parameter files](reference/parameter-files.md) | [Parameter and plotting profiles](design/decisions/parameter-and-plotting-profiles.md) and [application execution](reference/application-execution.md) |
| Review historical decisions or evidence | [Review history](history/reviews/README.md) | The relevant closure report, audit, or implementation record under `history/` |
| Edit Markdown or LaTeX | [Documentation instructions](../.agents/documentation.md) | The document being changed |
