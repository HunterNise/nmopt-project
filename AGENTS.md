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

Before changing anything:

1. Inspect `git status` and preserve existing user work.
2. Read [conventions/README.md](conventions/README.md).
3. Read every convention file that applies to the intended action.
4. Use [docs/README.md](docs/README.md) to select the relevant design and
   implementation documents; do not read every document by default.

The convention files are required pre-reading, not optional background:

- Code, headers, or tests: `conventions/code.md`
- CMake, build, or test configuration: `conventions/build.md`
- Markdown or LaTeX: `conventions/documentation.md`
- Branches, commits, or remotes: `conventions/git.md`

If a change crosses categories, read all applicable convention files first.

## Project authority

The convention files define how to work in the repository. The documents under
`docs/` define the project’s architecture, contracts, mathematical choices,
implementation scope, and roadmap. When a task requires a new architectural
decision, document the decision instead of silently extending an interface.
