# Documentation

`docs/` contains the project manual, current references, implementation maps,
design records, numerical studies, retained plans, and historical evidence.

There is no single required reading order. Start by identifying the kind of
information you need, then follow the corresponding file below.

## Structure

```text
docs/
├── manual/       teaching and orientation
│   ├── overview/     project architecture, scope, and experiment overview
│   ├── concepts/     project-specific mathematical and software concepts
│   └── background/   optional prerequisite material
├── reference/    current public programming, configuration, and execution interfaces
├── internals/    current implementation ownership and mechanics
├── design/       long-lived architecture, mathematics, and accepted decisions
├── studies/      source problems, applications, benchmarks, and numerical evidence
├── planning/     retained roadmaps and possible continuation work
└── history/      superseded work, reviews, audits, and closure evidence
```

Repository-level build and environment documentation sits outside `docs/`:

- [`BUILD.md`](../BUILD.md) — build profiles, testing, focused verification,
  application execution, and generated output.
- [`DEPENDENCIES.md`](../DEPENDENCIES.md) — installation, dependencies, and
  machine/environment setup.

## Where should I start?

| If you want to... | Start with |
| --- | --- |
| Understand the project as a whole | [`manual/overview/README.md`](manual/overview/README.md) |
| Check what is currently implemented | [`manual/overview/implemented-scope.md`](manual/overview/implemented-scope.md) |
| Learn the project-specific concepts in depth | [`manual/concepts/README.md`](manual/concepts/README.md) |
| Fill gaps in PDEs, FEM, optimization, C++, or deal.II | [`manual/background/README.md`](manual/background/README.md) |
| Build, test, or run the repository | [`../BUILD.md`](../BUILD.md) |
| Describe a new supported problem | [`reference/problem-authoring.md`](reference/problem-authoring.md) |
| Understand what the compiler can build | [`reference/compiler.md`](reference/compiler.md) |
| Use reduced optimization, KKT, OTD, or PDAS | [`reference/optimization.md`](reference/optimization.md) |
| Build a reusable `nmopt` application | [`reference/application-authoring.md`](reference/application-authoring.md) |
| Run/configure an existing application | [`reference/application-execution.md`](reference/application-execution.md) |
| Integrate an existing deal.II PDE code | [`reference/external-dealii-solver-integration.md`](reference/external-dealii-solver-integration.md) |
| Find where something is implemented | [`internals/implementation-map.md`](internals/implementation-map.md) |
| Understand a design or mathematical decision | [`design/`](design/) |
| Inspect Chapter 5/6 studies and evidence | [`studies/`](studies/) |
| Trace old plans, reviews, or superseded work | [`planning/`](planning/) and [`history/`](history/) |

## Manual

The manual is the best place to learn the project before reading APIs or source.

### `manual/overview/`

- [`manual/overview/README.md`](manual/overview/README.md) — recommended first
  reading; routes the overview pages and explains what each one contributes.
- [`manual/overview/project-architecture.md`](manual/overview/project-architecture.md)
  — high-level architecture, scope, and the common numerical boundary.
- [`manual/overview/architecture-maps.md`](manual/overview/architecture-maps.md)
  — visual maps of producer paths, formulation families, reduced evaluation,
  and experiment/evidence flow.
- [`manual/overview/semantic-compiler.md`](manual/overview/semantic-compiler.md)
  — the structured `ProblemSpec` → validation/resolution → deal.II compiler path.
- [`manual/overview/numerical-realization.md`](manual/overview/numerical-realization.md)
  — how spaces, operators, boundary conditions, metrics, observations, and
  linear solves become the numerical services consumed by formulations.
- [`manual/overview/reduced-optimization.md`](manual/overview/reduced-optimization.md)
  — the reduced state–adjoint workflow and one optimization iteration.
- [`manual/overview/experiments-and-replication.md`](manual/overview/experiments-and-replication.md)
  — scenarios, parameter files, runs, artifacts, provenance, and
  post-processing.
- [`manual/overview/external-applications.md`](manual/overview/external-applications.md)
  — how an existing PDE application can keep numerical ownership while using
  `nmopt`.
- [`manual/overview/implemented-scope.md`](manual/overview/implemented-scope.md)
  — concise boundary of what is implemented, demonstrated, and intentionally
  unsupported.

### `manual/concepts/`

- [`manual/concepts/README.md`](manual/concepts/README.md) — ordered concept
  chapters covering numerical spaces, duality, operators, metrics,
  formulations, optimization, semantic modeling, compilation, and integration.

Use this when the overview gives you the right mental model but you need the
mathematical/software ideas developed more carefully.

### `manual/background/`

- [`manual/background/README.md`](manual/background/README.md) — optional
  prerequisite routing for weak PDEs, finite elements, numerical linear
  algebra, unconstrained/constrained optimization, PDE-constrained
  optimization, modern C++, scientific-software workflow, and deal.II.

This is support material, not a required preface to the project manual.

## Reference

`reference/` describes **current public behavior**. Use these files when you
already know what you want to do.

- [`reference/problem-authoring.md`](reference/problem-authoring.md) — construct
  and validate semantic problem descriptions.
- [`reference/compiler.md`](reference/compiler.md) — supported deal.II compiler
  targets, capabilities, products, and unsupported-combination diagnostics.
- [`reference/optimization.md`](reference/optimization.md) — reduced
  optimization, supplied OTD, quadratic KKT, PDAS, metrics, constraints, and
  solver compatibility.
- [`reference/application-authoring.md`](reference/application-authoring.md) —
  define reusable recipes, typed scenarios, runtime data, and execution
  adapters.
- [`reference/application-execution.md`](reference/application-execution.md) —
  run already-authored applications; runner CLI, run sets, manifests,
  artifacts, native output, and post-processing.
- [`reference/parameter-files.md`](reference/parameter-files.md) — `.prm`
  configuration, experiment matrices, selections/exclusions, plotting profiles,
  and precedence rules.
- [`reference/external-dealii-solver-integration.md`](reference/external-dealii-solver-integration.md)
  — integrate an application-owned deal.II solver through numerical callbacks,
  solve services, metrics, and lifetime rules.

## Internals

`internals/` maps current architecture to source ownership and implementation
mechanics.

- [`internals/implementation-map.md`](internals/implementation-map.md) — best
  starting point for source-level work; maps the major layers to directories,
  headers, applications, and focused tests.
- [`internals/compiler.md`](internals/compiler.md) — compiler resolution,
  lowering stages, product construction, ownership, and extension points.
- [`internals/runner.md`](internals/runner.md) — parameter parsing, run-set
  planning, scenario resolution, execution, manifests, and artifact lifecycle.

For repository working conventions, also read [`AGENTS.md`](../AGENTS.md) and
the relevant files under [`.agents/`](../.agents/). Despite the name,
`.agents/` also records human-useful conventions for builds, Git,
documentation, verification, and source organization.

## Design

`design/` records long-lived architectural and mathematical decisions. These
documents explain **why the current interfaces and boundaries look the way
they do**.

- [`design/architecture.md`](design/architecture.md) — overall component model
  and architecture.
- [`design/composition-boundaries.md`](design/composition-boundaries.md) —
  ownership and dependency boundaries between numerical layers.
- [`design/interface-specification.md`](design/interface-specification.md) —
  foundational component contracts and interface vocabulary.
- [`design/mathematical-model.md`](design/mathematical-model.md) — conventions
  for residuals, derivatives, adjoints, covectors, metrics, and formulations.
- [`design/pde-solver-boundary.md`](design/pde-solver-boundary.md) — accepted
  separation between PDE realization, formulation, compiler, and optimization.
- [`design/decisions/repository-organization.md`](design/decisions/repository-organization.md)
  — rationale for the current repository/documentation layout.
- [`design/decisions/parameter-and-plotting-profiles.md`](design/decisions/parameter-and-plotting-profiles.md)
  — rationale for separating numerical `.prm` inputs from plotting policy.

## Studies

`studies/` contains the concrete problem families, benchmark definitions,
replication work, and numerical evidence. Use it to understand **what was
actually studied**, rather than the general framework API.

### `studies/chapter-5/`

- [`studies/chapter-5/source-catalogue.md`](studies/chapter-5/source-catalogue.md)
  — catalogue of the Chapter 5 source problem families used to shape the
  semantic/compiler work.
- [`studies/chapter-5/recipes.md`](studies/chapter-5/recipes.md) — mapping from
  those problem families to reusable application/compiler recipes.

### `studies/chapter-6/`

- [`studies/chapter-6/numerical-methods.md`](studies/chapter-6/numerical-methods.md)
  — Chapter 6 optimization/numerical-method material relevant to the project.
- [`studies/chapter-6/numerical-examples.md`](studies/chapter-6/numerical-examples.md)
  — source numerical examples and their project interpretation.
- [`studies/chapter-6/scenarios.md`](studies/chapter-6/scenarios.md) — typed
  application/scenario choices used by the repository.
- [`studies/chapter-6/benchmarks.md`](studies/chapter-6/benchmarks.md) — benchmark
  contracts, evidence expectations, and run families.
- [`studies/chapter-6/b1-replication.md`](studies/chapter-6/b1-replication.md) —
  B1 distributed-control reproduction record and project replacement choices.
- [`studies/chapter-6/b2-replication.md`](studies/chapter-6/b2-replication.md) —
  B2 boundary-control verification/replication record, including the preserved
  negative literal source-replication result.

### `studies/case-studies/`

- [`studies/case-studies/laplace-growth.md`](studies/case-studies/laplace-growth.md)
  — focused study of Laplace problem growth/behavior.
- [`studies/case-studies/laplace-interface-formulas.md`](studies/case-studies/laplace-interface-formulas.md)
  — concrete formulas used to check interface/sign conventions.

The external-application case study lives with the application itself:

- [`apps/external-dealii/step-4/README.md`](../apps/external-dealii/step-4/README.md)
  — overview, result, evidence, and reading paths for the Step-4 integration.

## Planning

`planning/` retains development plans and possible continuation work. These
files are useful for understanding intended sequencing and unfinished ideas,
but they are not automatically current implementation requirements.

- [`planning/future-extensions.md`](planning/future-extensions.md) — current
  summary of substantial capabilities left unimplemented, deferred, or
  conditional.
- [`planning/implementation-roadmap.md`](planning/implementation-roadmap.md) —
  retained core implementation roadmap.
- [`planning/application-roadmap.md`](planning/application-roadmap.md) —
  retained application/runner roadmap.
- [`planning/chapter-5-problem-library-roadmap.md`](planning/chapter-5-problem-library-roadmap.md)
  — retained plan for the Chapter 5 problem/recipe library.
- [`planning/chapter-6-benchmark-suite-roadmap.md`](planning/chapter-6-benchmark-suite-roadmap.md)
  — retained benchmark-suite plan, including later B3–B6 candidates.
- [`planning/repository-structure-refactor.md`](planning/repository-structure-refactor.md)
  — retained plan/record for the repository-structure refactor.

Read each file's own status before treating an item as active or incomplete.

## History

`history/` preserves superseded implementation records, reviews, audits, and
closure evidence. It is primarily for provenance and reconstruction.

Useful entry points are:

- [`history/reviews/README.md`](history/reviews/README.md) — index of review,
  audit, and closure series.
- [`history/implementation/`](history/implementation/) — superseded
  implementation records kept for provenance and comparison with the current
  code.
- [`history/documentation-lineage.md`](history/documentation-lineage.md) —
  major documentation-structure transitions and recovery points.
- [`history/application-and-gui-layer-draft.md`](history/application-and-gui-layer-draft.md)
  — historical application/GUI proposal; useful for provenance, not as a
  statement of current architecture.

## Which documents are authoritative?

For current implementation behavior:

1. current source and focused tests are primary;
2. `reference/` describes the supported public surface;
3. `internals/` describes current implementation ownership/mechanics;
4. `design/` records accepted long-lived architecture and mathematical
   conventions;
5. `studies/` owns numerical/application evidence;
6. `planning/` and `history/` explain intended or historical states and should
   not override current source/reference behavior.

Historical documents may intentionally preserve terminology, paths, and
assumptions from the state they record.
