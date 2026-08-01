# Project guide

## Mission

Build a deal.II-based framework for PDE-constrained optimal control and
inverse problems. The project should make combinations of PDE operators,
controls, observations, norms, constraints, and solvers reusable without
creating a new problem class for every combination.

The governing principle is **composition of residual, objective, metric,
constraint, and discretization components—not inheritance from particular PDE
problem types**.

## Required routing

Before changing anything, inspect `git status`, preserve existing user work,
and follow the mandatory routing rules in
[conventions/README.md](conventions/README.md). That file is the single source
of truth for action-specific convention routing; do not duplicate its
convention table here.

Use [docs/README.md](docs/README.md) to select the relevant design and
implementation documents. Do not read every document by default.

## Project authority

The convention files define how to work in the repository. The documents under
`docs/` define the project’s architecture, contracts, mathematical choices,
implementation scope, and roadmap. When a task requires a new architectural
decision, document the decision instead of silently extending an interface.
