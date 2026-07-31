# nmopt-project

A deal.II-oriented framework for PDE-constrained optimal control and inverse
problems. The aim is to support reusable combinations of PDE operators,
boundary conditions, controls, observations, norms, constraints, and solvers
without creating a new `Problem` class for every combination.

The project has a backend-parametric executable DTO contract, dense reference
test, and a serial deal.II lowerer for the first scalar finite-element slice.

Start here:

- [Documentation map](docs/README.md) — which document applies to which
  audience and task.
- [Architecture overview](docs/architecture.md) — the detailed design record.
- [Interface specification](docs/interface-specification.md) — normative
  component contracts and compilation protocols.
- [Implementation roadmap](docs/implementation-roadmap.md) — current state,
  prioritized tasks, and acceptance checks.
- [Project guide for contributors and agents](AGENTS.md) — concise working,
  architecture, and formatting rules.
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
