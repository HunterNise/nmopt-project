# nmopt-project

A deal.II-oriented framework for PDE-constrained optimal control and inverse
problems. The aim is to support reusable combinations of PDE operators,
boundary conditions, controls, observations, norms, constraints, and solvers
without creating a new `Problem` class for every combination.

The project has a small backend-neutral executable DTO contract and reference
test. A deal.II compiler for the first finite-element slice remains to be
implemented.

Start here:

- [Architecture overview](docs/architecture.md) — the detailed design record.
- [Theoretical formalism](docs/theoretical-formalism.md) — the abstract
  problem, its relation to textbook strong forms, and worked special cases.
- [Composition boundaries](docs/composition-boundaries.md) — actionable
  component ownership and unavoidable cross-cutting interactions.
- [Laplace growth case study](docs/laplace-growth-case-study.md) — a concrete
  semantic-to-solver walkthrough, extensions, and DTO/OTD comparison.
- [Laplace formula deltas](docs/laplace-interface-formulas.md) — exact
  formula changes, compiled operators, and required interface primitives.
- [Interface specification](docs/interface-specification.md) — normative
  component contracts, compilation protocols, and solver-facing outputs.
- [Implementation-readiness review and default policies](docs/implementation-readiness-review.md)
  — resolved discrete-algebra, compilation, formulation, metric, constraint,
  boundary, and rollout decisions for the first implementation.
- [V0 executable contract](docs/executable-contract-v0.md) — the implemented
  typed value/JVP/VJP, metric, constraint, and reduced-DTO boundary.
- [Project guide for contributors and agents](AGENTS.md) — concise working
  instructions and non-negotiable design decisions.

The governing principle is **composition of residual, objective, metric,
constraint, and discretization components—not inheritance from particular PDE
problem types**.

## Validate the contract

~~~bash
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build
ctest --test-dir build --output-on-failure
~~~
