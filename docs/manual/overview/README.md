# Project overviews

The documents in this directory are the recommended starting point for understanding
`nmopt` as a working project.

They are not API references and they are not implementation roadmaps. Their job is
to give a first-time reader a usable mental model:

- what kind of PDE-constrained optimization problems the project is intended to
  represent;
- what parts of the numerical work `nmopt` owns and what remains in an application;
- the two main ways a problem reaches the common optimization machinery;
- how discretization, state and adjoint solves, metrics, and optimization fit
  together; and
- how the project-specific runner adds reproducible experiments around the numerical
  core.

A reader should be able to finish these pages and then decide whether they want to
integrate an existing PDE code, describe a problem through the semantic/compiler
path, inspect the optimization machinery, or reproduce the repository's Chapter 6
experiments.

## Suggested reading order

Start with [Project architecture](project-architecture.md). It explains the scope of
the project and introduces the common numerical boundary shared by all problem
sources.

If you prefer to orient yourself visually, [Architecture maps](architecture-maps.md)
starts from the same high-level picture and progressively zooms into the main
producer paths, formulation families, reduced-evaluation lifecycle, and
experiment/evidence flow.

Then read the pages that match your interest:

- [Integrating an existing PDE application](external-applications.md) if you already
  have a numerical PDE code and want to retain its mesh, assembly, solves, and
  native output while using `nmopt` formulations and algorithms.
- [Semantic and compiler path](semantic-compiler.md) if you want to describe a
  supported PDE optimal-control problem declaratively and let the framework build
  the deal.II realization.
- [Numerical realization](numerical-realization.md) if you want to understand how
  finite-element coordinates, operators, boundary conditions, observations,
  metrics, and linear solves become the numerical operations used by optimization.
- [Reduced optimization](reduced-optimization.md) if you want to understand the
  state–adjoint reduced method and what one optimization iteration actually does.
- [Experiments and replication](experiments-and-replication.md) if you want to
  understand scenarios, parameter files, benchmark runs, artifacts, provenance,
  and post-processing.

The pages deliberately overlap at their boundaries. A small amount of repetition is
useful here: the overview set should remain understandable without requiring a
reader to chase links for every central concept.

The two problem-entry routes are intentionally given equal architectural weight.
The semantic/compiler path is the framework's structured authoring route; the
external-application path is the route for preserving an independently owned
numerical code. The experiment and replication layer is different: it sits outside
those numerical boundaries and turns selected capabilities into reproducible
project runs and benchmarks.

## Repository-specific names used in these pages

A few names recur because the repository grew around concrete source material and
integration studies:

- **Chapter 5** and **Chapter 6** refer to the corresponding chapters of
  *Optimal Control of Partial Differential Equations* by Manzoni, Quarteroni, and
  Salsa. Chapter 5 supplies much of the elliptic optimal-control problem catalogue;
  Chapter 6 supplies numerical methods and the numerical examples used for the
  repository's reproduction work.
- **B1** and **B2** are repository scenario IDs for two selected Chapter 6 examples:
  B1 is the distributed Laplace-control example and B2 is the Graetz-flow
  boundary-control example.
- **Step-4** refers to the deal.II tutorial program `step-4`, adapted in this
  repository as a small external-application integration case study. **Problem B** is
  the richer distributed-control problem built around that adapted application.

These labels identify examples; they are not architectural layers of `nmopt`.

## What these pages assume

The overviews use the normal language of numerical PDEs and optimization: weak
forms, finite-element degrees of freedom, state and adjoint equations, derivatives,
metrics, line search, and similar concepts. They explain how those ideas are used
inside this repository, but they are not intended to teach all prerequisite theory
from first principles.

The prerequisite theory is intentionally kept separate from the architectural
story. These overviews should remain technically accurate rather than replacing
project terminology with simpler but less precise analogies; dedicated background
recaps can provide the slower mathematical and C++/deal.II refreshers.

## Read later: authoritative documentation

Once the mental model is clear, route by the kind of authority you need:

- [`docs/design/theoretical-formalism.md`](../design/theoretical-formalism.md) is the
  long-lived mathematical convention for residuals, adjoints, covectors, metrics,
  and formulations.
- [`docs/design/pde-solver-boundary.md`](../design/pde-solver-boundary.md) records
  the accepted ownership boundary between PDE realization, formulation, compiler,
  and optimization.
- [`docs/implementation/v1/semantic-compiler.md`](../implementation/v1/semantic-compiler.md)
  is the detailed current capability ledger for the v1 semantic/compiler path.
- `docs/reference/` contains the exact public operational contracts that already have
  dedicated references, including external deal.II integration, application
  execution, and parameter files.
- `docs/applications/`, `docs/guides/`, and `docs/benchmarks/` own the concrete
  Chapter 5/6 problem families, source transcriptions, and reproduction evidence.
- `docs/planning/review/` preserves audits and historical review evidence; it is not
  the default source for current architecture.

The overview layer should tell you where to look next, but it should already provide
enough context to understand why that location is relevant.
