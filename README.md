# nmopt-project

A deal.II-oriented framework for PDE-constrained optimal control and inverse
problems. The aim is to support reusable combinations of PDE operators,
boundary conditions, controls, observations, norms, constraints, and solvers
without creating a new `Problem` class for every combination.

The project has a backend-parametric executable DTO contract, dense reference
test, and a serial deal.II lowerer for the first scalar finite-element slice.

Start here:

- [System blueprint](docs/system-blueprint.md) — theory, specification, v0
  implementation, and test correspondence in one implementer-oriented guide.
- [Documentation map](docs/README.md) — which document applies to which
  audience and task.
- [Architecture overview](docs/architecture.md) — the detailed design record.
- [Interface specification](docs/interface-specification.md) — normative
  component contracts and compilation protocols.
- [Implementation roadmap](docs/implementation-roadmap.md) — current state,
  prioritized tasks, and acceptance checks.
- [Project conventions](conventions/README.md) — action-specific guidance for
  code, builds, documentation, and Git operations.
- [Project guide for contributors and agents](AGENTS.md) — concise working,
  mission and required routing.
