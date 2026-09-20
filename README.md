# nmopt-project

A deal.II-oriented framework for PDE-constrained optimal control and inverse
problems. The project is designed around reusable combinations of PDE
operators, boundary conditions, controls, observations, norms, constraints,
and solvers rather than a separate problem class for every combination.

## Current scope

The current implementation is centered on serial deal.II finite-element
realizations behind backend-neutral numerical contracts. It includes:

- a backend-neutral executable model and formulation vocabulary;
- a semantic `ProblemSpec` layer with validation, resolution, and a bounded v1
  deal.II compiler;
- reduced state–adjoint, supplied-OTD, quadratic KKT, complementarity, and PDAS
  formulation surfaces with their required numerical services;
- reduced optimization algorithms, metrics, constraints, and solve-report contracts;
- both compiler-owned and independently owned application paths to the common
  formulation boundary; and
- recipes, scenarios, a runner, and reproduction infrastructure for the current
  Chapter 5/6 application work.

Semantic validity is broader than the combinations implemented by the current
compiler, and benchmark coverage is narrower than framework capability. The
[project manual](docs/manual/overview/README.md) explains those boundaries before
the exact design and reference documents.

## Installation

> [!WARNING]
> This project is meant to run in a Linux environment. Windows users should
> install WSL2 with an Ubuntu distribution and work inside the Linux
> environment. See the [WSL installation guide](https://learn.microsoft.com/en-us/windows/wsl/install).

Read [Dependencies and environment](DEPENDENCIES.md) for the complete
toolchain, installation options, version checks, and runtime requirements.

## Build and test

The preferred local build and test entry point is the root-level `build.sh`
helper. It requires a C++ compiler, CMake, and Ninja. The deal.II package is
needed only for the deal.II-enabled profiles.

Create the machine-local configuration once before running a build pipeline:

```bash
./build.sh init-config
```

The default backend-neutral verification loop is:

```bash
./build.sh pipeline debug-neutral
```

The helper also supports the deal.II and sanitizer profiles, multiple profiles
in one invocation, and the complete profile set:

```bash
./build.sh pipeline debug-dealii
./build.sh pipeline debug-neutral sanitize-neutral
./build.sh all
```

Use `configure`, `build`, or `test` when only one phase is needed. These are
thin wrappers around the corresponding CMake or CTest preset commands:

```bash
./build.sh configure debug-neutral
./build.sh build debug-neutral
./build.sh test debug-neutral
```

The configuration file records machine-specific settings such as the
deal.II package directory and the maximum number of build jobs for each
profile. Run `./build.sh show-config` to inspect the active values. The
pipeline additionally times builds and uses compact configure output plus
progress-oriented test output. Pass phase-specific options such as `--target`
to `build` or `--ctest-arg=--output-on-failure` to `test`; `--dry-run` and
`--verbose` are useful with any action when more control is needed.

The helper maps directly to the checked-in CMake presets. The equivalent
manual neutral-profile commands are:

```bash
cmake --preset debug-neutral
cmake --build --preset debug-neutral
ctest --preset debug-neutral
```

Manual commands remain useful for one-off CMake/CTest options or a workflow
that does not fit the helper. Unlike `build.sh`, they do not read
`build.local.conf` or apply its per-profile job limit automatically.

The deal.II backend uses the separate `debug-dealii` preset and requires an
official deal.II CMake package. The helper, dependency setup, generated
output, focused tests, and environment-specific notes are documented in the
[dependencies and environment reference](DEPENDENCIES.md), the [application
execution reference](docs/reference/application-execution.md), and the
[build instructions](.agents/build.md).

## Running applications

After building the `debug-dealii` profile, run a Chapter 6 development case
with the repository helper:

```bash
tools/run_chapter6.sh --benchmark b1 --refinement 1
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

- [Project manual](docs/manual/overview/README.md) — the recommended first-time
  reading path from the project architecture into the concept chapters.
- [Concept chapters](docs/manual/concepts/README.md) — the progressive mathematical,
  formulation, compiler, and integration manual.
- [System blueprint](docs/internals/system-blueprint.md) — the shortest
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
recipes](docs/studies/chapter-5/recipes.md), [Chapter 6 application
scenarios](docs/studies/chapter-6/scenarios.md), [Chapter 6 benchmark
specifications](docs/studies/chapter-6/benchmarks.md), and the [application
roadmap](docs/planning/application-roadmap.md). The [Chapter 6 numerical
examples](docs/studies/chapter-6/numerical-examples.md) guide records what the
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
| `docs/` | Manual, design, reference, internals, studies, planning, and historical records. |
| `.agents/` | Detailed instructions and prompt references for coding-agent work. |
| `build/` | Ignored, profile-specific CMake and build output. |
| `runs/` | Ignored generated application-run artifacts and reports. |

The [repository organization decision](docs/design/decisions/repository-organization.md)
records the ownership boundaries in more detail.

## Working with coding agents

[`AGENTS.md`](AGENTS.md) is the concise agent entry point and routes to the
detailed files under `.agents/`. The [prompt cookbook](.agents/prompts.md) is a
human-maintained collection of useful prompts; it is not part of the normal
agent instruction routing.
