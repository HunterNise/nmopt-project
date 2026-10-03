# nmopt

> [!NOTE]
> **Course project.** This repository was developed for the [*Numerical Methods
> for Optimal Control*](https://luca-heltai.github.io/numerical-methods-for-optimal-control/)
> course. It is a completed educational/research project, not a maintained production library.

> [!WARNING]
> **Agent-authored code.** Most of the implementation and a substantial part of
> the documentation were written by coding agents under human direction. The
> repository includes automated tests, numerical checks, application studies,
> build/run evidence, and later review passes, but the implementation was not
> independently reviewed line by line.

`nmopt` is a C++17/deal.II project for discretized PDE-constrained optimal
control and inverse problems.

## Why this project

Course examples are often implemented as separate programs: each one owns its
PDE discretization, solvers, optimization loop, constraints, and output.

The project instead asks:

**Which parts are genuinely problem-specific, and which numerical interfaces
can be reused across different optimal-control problems?**

The project separates three concerns:

| Concern | Owns |
| --- | --- |
| **PDE realization** | FE spaces, matrix assembly, boundary conditions, state/adjoint solves, native fields |
| **Optimal-control formulation** | reduced systems, supplied OTD, KKT systems, complementarity |
| **Optimization** | search directions, globalization, constraints, stopping criteria |

That separation led to two peer ways of supplying the PDE side:

```text
structured problem ──► nmopt deal.II compiler ──┐
                                                ├──► formulations ──► solvers
existing deal.II app ─► thin numerical adapter ─┘
```

An existing application can therefore keep its own mesh, assembly, solves, and
output; it does not have to be rewritten around the framework's compiler.

## What `nmopt` is today

| Area | Current project |
| --- | --- |
| **Numerical core** | Reusable contracts for PDE actions, derivatives/adjoints, state and adjoint solves, metrics, constraints, and solver evidence |
| **Formulations** | Reduced state–adjoint, supplied optimize-then-discretize (OTD), quadratic KKT, complementarity/PDAS |
| **Optimization** | Steepest descent, nonlinear CG, BFGS/L-BFGS, Newton, line searches, trust region, projected search where supported |
| **deal.II realization** | Registered scalar elliptic/control families, several volume/boundary controls, observations, metrics, and constraints |
| **Applications** | Typed recipes/scenarios, parameter files, headless execution, native FE output, run manifests, post-processing and reports |
| **Evidence** | Chapter 5/6 problem studies, B1/B2 numerical replication work, and an external deal.II Step-4 integration |

The strongest implemented and demonstrated region is **serial, scalar,
steady/elliptic finite-element optimal control**. The compiler is deliberately
bounded: a mathematically valid problem is not automatically a supported
deal.II realization.

See [Implemented scope](docs/manual/overview/implemented-scope.md) for the
detailed capability boundary, and
[Future extensions](docs/planning/future-extensions.md) for directions that
were left unimplemented or deferred.

## Repository guide

| Path | What to expect |
| --- | --- |
| [`include/nmopt/`](include/nmopt/) | Reusable contracts, formulations, solvers, semantic model, compiler, deal.II services, application interfaces |
| [`apps/`](apps/) | Concrete applications, repository runner, external-integration examples |
| [`parameters/`](parameters/) | Checked-in application and experiment configuration |
| [`tests/`](tests/) | Contract, compiler, solver, deal.II, application, and tooling verification |
| [`tools/`](tools/) | Run helpers, post-processing, reporting, documentation utilities |
| [`docs/`](docs/) | Manual, reference, design, internals, studies, planning, and history |
| [`.agents/`](.agents/) | Repository conventions for code, Git, builds, documentation, explanations, and agent-assisted work |
| [`cmake/`](cmake/) | Repository-specific CMake support |

Generated build and run output lives under the ignored `build/` and `runs/`
directories.

For source-level ownership and dependency flow, use the
[implementation map](docs/internals/implementation-map.md).

## Start here

**New to the project**

- [Project overviews](docs/manual/overview/README.md) — architecture and main ideas.
- [Concept chapters](docs/manual/concepts/README.md) — project-specific mathematical and software concepts.
- [Background guide](docs/manual/background/README.md) — optional PDE, FEM, optimization, C++, and deal.II prerequisites.

**Want to use it**

- [Problem authoring](docs/reference/problem-authoring.md) — describe a supported problem.
- [Compiler](docs/reference/compiler.md) — build its deal.II realization.
- [External deal.II integration](docs/reference/external-dealii-solver-integration.md) — keep an existing PDE application in control.
- [Optimization](docs/reference/optimization.md) — consume reduced, KKT, PDAS, and supplied-OTD products.

**Want to build, test, or run the included applications**

- [Build, test, and run](BUILD.md) — build profiles, focused verification, the runner, and generated output.
- [Dependencies and environment](DEPENDENCIES.md) — installation and machine setup.
- [Application execution](docs/reference/application-execution.md) — exact runner, run-set, manifest, and artifact contracts.
- [Parameter files](docs/reference/parameter-files.md) — tracked application/run configuration.
- [Repository tools](tools/README.md) — run helpers, post-processing, reporting, and linting.

**Want to inspect the project itself**

- [Documentation map](docs/README.md) — documentation organization and routing.
- [Implementation map](docs/internals/implementation-map.md) — map architecture and responsibilities to source directories, applications, and focused tests.
- [Numerical studies](docs/studies/) — Chapter 5/6 problem material, benchmarks, replication records, and case studies.
- [`AGENTS.md`](AGENTS.md) and [`.agents/`](.agents/) — repository conventions for code, Git, builds, documentation, verification, and agent-assisted work.

## Build and test

> The supported development environment is Linux; Windows users can use WSL2.

See [DEPENDENCIES.md](DEPENDENCIES.md) for installation and environment setup.
The canonical build, test, and execution workflow is documented in
[BUILD.md](BUILD.md).

A first backend-neutral check is:

```bash
./build.sh init-config
./build.sh pipeline debug-neutral
```

With deal.II available:

```bash
./build.sh pipeline debug-dealii
```

## Run

After building the deal.II profile, a small development run can be started with:

```bash
tools/run_chapter6.sh --benchmark b1 --refinement 1
```

The helper runs the application, post-processes its native output, and writes
the generated run set under `runs/`.

See [BUILD.md](BUILD.md) for choosing a build profile and running the included
applications, and [Application execution](docs/reference/application-execution.md)
for the exact runner, manifest, artifact, and output contracts.

## License

Project code is released under the MIT License.

Third-party material retains its original copyright and licensing terms.
Preserved upstream deal.II tutorial source remains under the license stated in
those files.
