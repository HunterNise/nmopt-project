# nmopt-project

A future deal.II framework for PDE-constrained optimal control and inverse
problems. The aim is to support reusable combinations of PDE operators,
boundary conditions, controls, observations, norms, constraints, and solvers
without creating a new `Problem` class for every combination.

The project is currently in the architecture stage: no implementation choices
below the public model/compilation boundaries are fixed yet.

Start here:

- [Architecture overview](docs/architecture.md) — the detailed design record.
- [Project guide for contributors and agents](AGENTS.md) — concise working
  instructions and non-negotiable design decisions.

The governing principle is **composition of residual, objective, metric,
constraint, and discretization components—not inheritance from particular PDE
problem types**.
