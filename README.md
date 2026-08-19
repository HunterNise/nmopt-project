# nmopt-project

A deal.II-oriented framework for PDE-constrained optimal control and inverse
problems. The project is designed around reusable combinations of PDE
operators, boundary conditions, controls, observations, norms, constraints,
and solvers rather than a separate problem class for every combination.

## Current scope

The current implementation establishes the first scalar finite-element slice:

- a backend-neutral executable DTO and formulation contract;
- semantic validation and v1 graph/compiler interfaces;
- a dense reference path and focused contract tests;
- a serial deal.II lowerer and headless application runner;
- the initial Chapter 5/6 application and benchmark layer.

The application and benchmark layer is still being extended. Its current
status and acceptance work are tracked in the [application roadmap](docs/planning/application-roadmap.md)
and [implementation roadmap](docs/planning/implementation-roadmap.md).

## Installation

> [!WARNING]
> This project is meant to run in a Linux environment. Windows users should
> install WSL2 with an Ubuntu distribution and work inside the Linux
> environment. See the [WSL installation guide](https://learn.microsoft.com/en-us/windows/wsl/install).

Read [Dependencies and environment](DEPENDENCIES.md) for the complete
toolchain, installation options, version checks, and runtime requirements.

## Build and test

The fastest local verification loop is the backend-neutral CMake preset:

It requires a C++ compiler, CMake, and Ninja. The deal.II package is needed
only for the deal.II-enabled presets.

```bash
cmake --preset debug-neutral
cmake --build --preset debug-neutral
ctest --preset debug-neutral
```

The deal.II backend uses the separate `debug-dealii` preset and requires an
official deal.II CMake package. Build profiles, dependency setup, generated
output, focused tests, and environment-specific notes are documented in the
[build instructions](.agents/build.md).

## Running applications

After building the `debug-dealii` profile, run a Chapter 6 development case
with the repository helper:

```bash
tools/run_chapter6.sh --benchmark b1 --refinement 1 --format png svg
```

Replace `b1` with `b2` for the second benchmark. The helper runs the compiled
application, post-processes its native fields, and writes artifacts and
reports under the ignored `runs/` directory. Python post-processing
dependencies are required for this step.

See the [application execution reference](docs/reference/application-execution.md)
for run schemas, output organization, direct runner commands, and report
details.

## Where to start

The [documentation map](docs/README.md) is the task-oriented index. The most
useful entry points are:

- [System blueprint](docs/design/system-blueprint.md) — the shortest
  implementer-oriented overview of the theory, specification, implementation,
  and test correspondence.
- [Architecture overview](docs/design/architecture.md) — long-lived design
  rationale and component boundaries.
- [Interface specification](docs/design/interface-specification.md) —
  normative component contracts and compilation protocols.
- [Application assembly API](docs/reference/application-api.md) — how recipes,
  scenarios, compilers, solvers, and experiments are assembled.
- [Application execution reference](docs/reference/application-execution.md) —
  run schemas, generated artifacts, reports, and post-processing.

For application and benchmark work, use the [Chapter 5 application
recipes](docs/applications/chapter-5.md), [Chapter 6 application
scenarios](docs/applications/chapter-6.md), [Chapter 6 benchmark
specifications](docs/benchmarks/chapter-6.md), and the [application
roadmap](docs/planning/application-roadmap.md). The [Chapter 6 numerical
examples](docs/guides/chapter-6-numerical-examples.md) guide records what the
book says; implementation choices and reproduction status belong in the
application and benchmark planning documents.

## Repository layout

| Path | Role |
| --- | --- |
| `include/nmopt/` | Public C++ headers, organized by semantic, compiler, backend, solver, and application layer. |
| `apps/` | Executable entry points and application-specific orchestration. |
| `tests/` | Contract, semantic, backend, application, and benchmark tests. |
| `tools/` | Run generation, post-processing, reporting, and related utilities. |
| `cmake/` | CMake helpers for scenario discovery and generated test registration. |
| `docs/` | Design, reference, guides, applications, benchmarks, decisions, and planning records. |
| `.agents/` | Detailed instructions and prompt references for coding-agent work. |
| `build/` | Ignored, profile-specific CMake and build output. |
| `runs/` | Ignored generated application-run artifacts and reports. |

The [repository organization decision](docs/decisions/repository-organization.md)
records the ownership boundaries in more detail.

## Working with coding agents

[`AGENTS.md`](AGENTS.md) is the concise agent entry point and routes to the
detailed files under `.agents/`. The [prompt cookbook](.agents/prompts.md) is a
human-maintained collection of useful prompts; it is not part of the normal
agent instruction routing.
