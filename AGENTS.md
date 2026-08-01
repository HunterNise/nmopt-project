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

Before any repository action, follow the mandatory routing rules in
[conventions/README.md](conventions/README.md). That file is the single
source of truth for action-specific convention routing; do not duplicate its
convention table here.

Use [docs/README.md](docs/README.md) to select the relevant design and
implementation documents. Do not read every document by default.

For read-only tasks such as inspecting, explaining, reviewing, or diagnosing,
read the relevant authoritative code and documentation first. Do not modify
repository files unless the user requests a change.

## Project authority

The convention files define how to work in the repository. The documents under
`docs/` define the project’s architecture, contracts, mathematical choices,
implementation scope, and roadmap.

When a task requires a new architectural decision,
document the decision instead of silently extending an interface.

## Chat handoffs

Lead chat replies with the concrete outcome. Report what was changed, read,
verified, or decided, along with any meaningful failure, blocker, or next
choice.

Assume the user knows the repository conventions. Omit routine compliance
statements such as preserving existing work, leaving changes uncommitted, or
passing a basic diff check unless they affect the result or were explicitly
requested.

For read-only tasks, report the findings and their authoritative sources
directly. Report negative results when they are explicit test outcomes,
failures, blockers, or requested absence checks.
