# Repository organization and authority boundaries

## Status

Accepted repository decision, updated to the current documentation and
application layout.

## Decision

Organize the repository by architectural layer and user-facing application
role. Keep long-lived mathematical and semantic contracts separate from
realized implementation records, mutable plans, and generated output.

The current top-level organization is:

```text
AGENTS.md                         # repository-wide agent instructions
.agents/                         # mandatory action-specific agent instructions
build.sh                          # preferred local build/test orchestration
build.local.conf.example         # tracked machine-configuration template
build.local.conf                  # ignored machine-local build configuration
docs/
  manual/                         # current project manual and concept chapters
  design/                         # long-lived architecture and mathematics
  reference/                      # exact public API and execution references
  internals/                      # current implementation mechanics and maintainer maps
  studies/                        # Chapter 5/6 source, application, benchmark, and case-study corpus
  planning/                       # retained development roadmaps and planning records
  history/                        # superseded implementations, reviews, audits, and evidence
include/nmopt/
  application/                    # public recipe, scenario, and catalog API
  contract/                       # backend-neutral executable/formulation ports
  semantic/v1/                    # ProblemSpec and semantic validation
  compiler/v1/                    # lowering and compiled products
  dealii/                         # deal.II backend services
  solvers/                        # optimization and Krylov policies
  experiment/                     # provenance and report envelopes
tests/
  support/                        # shared test utilities
  contract/                       # backend-neutral contract tests
  semantic/                       # semantic graph and validation tests
  compiler/dealii/                # deal.II compiler and lowering capability tests
  dealii/                         # low-level deal.II backend/service/session/reporting tests
  solvers/                        # reduced/KKT/PDAS solver contracts
  application/                    # application, runner, benchmark, and integration tests
  tools/                          # persisted-evidence and post-processing tests
apps/
  nmopt-runner/                   # headless application orchestration
  external-dealii/                # external-application integration study
parameters/                        # versioned experiment and plotting inputs
runs/                              # ignored generated output
```

Directories are created when their first authoritative file or implementation
unit lands. Empty placeholder directories are not tracked.

## Authority boundaries

- `AGENTS.md` owns repository-wide agent instructions.
- `.agents/README.md` is the mandatory routing entry point for workflow,
  Git, code, build, documentation, and explanation instructions.
- `docs/design/` owns long-lived semantic, mathematical, and architectural
  contracts.
- `docs/reference/` owns exact public types, options, schemas, and execution
  interfaces.
- `docs/manual/` teaches the current project through overviews and concept
  chapters.
- `docs/internals/` records current implementation mechanics, source ownership,
  lifetime boundaries, extension points, and focused verification maps; exact
  public contracts remain in `docs/reference/`.
- `docs/studies/` contains the Chapter 5/6 source, application, benchmark, and
  case-study corpus.
- `docs/planning/` contains retained development roadmaps and planning records;
  each document's status banner determines whether it describes active or
  historical/intended work.
- `docs/history/` contains superseded implementations, reviews, audits, and
  historical evidence.

The Chapter 5 and Chapter 6 study documents remain the source and application
catalogue. Scenario, benchmark, and replication records add executable choices
and evidence requirements without silently rewriting the source catalogue.
Retained roadmaps record development intent, sequencing, and completion
history; current capability or work status must not be inferred solely from
their location under `planning/`.

The root-level `build.sh` owns the convenient configure/build/test workflow and
delegates the actual profile definitions to `CMakePresets.json`. The tracked
`build.local.conf.example` documents the available machine-specific settings;
`build.local.conf` is the ignored, user-edited instance created by
`./build.sh init-config`. Build orchestration intentionally remains at the
repository root rather than becoming part of `tools/`, whose role is
application execution, post-processing, and reporting.

## Source, tooling, and generated-state boundaries

The non-documentation directories have these stable roles:

| Directory | Role | Boundary |
| --- | --- | --- |
| `include/nmopt/` | Public and internal C++ headers under the `nmopt` include and namespace prefix. | Reusable library interfaces belong here; executable-specific CLI and run-set code does not. |
| `apps/` | Source files for executable products and application-local integration studies. | `apps/nmopt-runner/main.cc` is the headless entry point; its adjacent `run_lifecycle.hpp` is private CLI/run-set support. `apps/external-dealii/` contains the application-owned integration study. Reusable interfaces remain under `include/nmopt/`. |
| `tests/` | Contract-layer test translation units and shared test support. | CTest inventories are discovered from test executables; generated registrations belong in `build/`, not in source. |
| `tools/` | Repository-local user and analysis tools, including Python post-processing and shell wrappers. | Tools consume public artifacts and interfaces; they must not become a second PDE lowerer or solver implementation. |
| `cmake/` | Source-controlled CMake modules used by the root build. | Build configuration remains expressed through root `CMakeLists.txt` and `CMakePresets.json`; generated CMake state does not belong here. |
| `build/` | Ignored generated build trees, one named directory per configured preset. | Never configure directly in the `build/` container; use `build/debug-neutral/`, `build/debug-dealii/`, `build/sanitize-neutral/`, or `build/release-dealii/`. |
| `runs/` | Ignored generated experiment evidence. | Run manifests, artifacts, native output, reports, and post-processing follow the [execution reference](../../reference/application-execution.md). |

The standard generated shapes are:

```text
build/<preset>/
  bin/                         # executable targets
  lib/                         # library targets
  generated/scenarios/         # generated CTest registrations

runs/chapter-6/<benchmark>/<run-slot>/
  run-manifest.json
  artifacts/
    <method-or-case>/...
  report/
  postprocess/
```

The [build instructions](../../../.agents/build.md) own profile selection,
cache recovery, generated-file rules, manual-command safety rules, and build
commands. The [application
execution reference](../../reference/application-execution.md) owns artifact
schemas, run paths, native output, reports, and post-processing commands. This
decision records the directory ownership without duplicating those procedures.

The `nmopt` include prefix mirrors the C++ namespace and gives external users
stable, collision-resistant include paths such as
`#include "nmopt/application/runner.hpp"`. Versioned semantic paths such as
`semantic/v1/` allow a later public generation to coexist with the current one.
The `apps/nmopt-runner/` grouping similarly keeps the executable entry point
and private CLI helpers together while leaving reusable application interfaces
under `include/nmopt/application/`.

## Public reference layout

The `reference/` directory intentionally contains seven task-oriented public
references:

- `problem-authoring.md` — semantic graph construction, validation, recipes,
  and stable-ID composition;
- `compiler.md` — runtime bindings, deal.II lowering, compiled products,
  manifests, and compiler diagnostics;
- `application-authoring.md` — typed application/scenario construction,
  backend execution adapters, and optional runner registration;
- `external-dealii-solver-integration.md` — application-owned deal.II
  callbacks, solve services, metrics, formulation binding, and native output;
- `optimization.md` — reduced search, trust region, KKT, PDAS, and supplied-OTD
  consumption;
- `application-execution.md` — run organization, artifacts, manifests, native
  output, reports, and post-processing;
- `parameter-files.md` — parameter-file and plotting-profile configuration.

There is no additional directory README because `docs/README.md` is the
repository-wide documentation map and these files are self-describing.

## Implementation internals layout

The `internals/` directory intentionally contains three maintainer-oriented
maps:

- `implementation-map.md` — repository-wide implementation ownership,
  dependency flow, lifetime landmarks, tests, and change routing;
- `compiler.md` — semantic resolution, request closure, bounded planning,
  target realization, compiled products, provenance, and compiler verification;
- `runner.md` — parameter resolution, run-set planning, runner lifecycle,
  Chapter-6 execution, artifacts, and run-manifest mechanics.

They explain how the current implementation works and where to modify it.
They do not replace the task-oriented public references, and the current
source plus focused tests remain authoritative for implementation behavior.

## Include and test layers

The current `include/nmopt/{contract,semantic,compiler,dealii,solvers,experiment}`
layout follows the framework's layer boundaries and is not renamed by this
decision. There is no public `include/nmopt/reference/` reference-model
surface; test-only reference models live under
`tests/support/reference_models/`. New public recipe and scenario records belong under
`include/nmopt/application/`; they may refer to semantic and compiler ports,
but must not move backend details into `ProblemSpec` or add PDE-specific
branches to generic solvers.

Existing test executables retain their names and scenario inventories while
source files are grouped by the contract layer they verify. Compiler and
lowering capability scenarios belong under `tests/compiler/dealii/`; low-level
deal.II backend, service, session, reporting, and lifetime scenarios belong
under `tests/dealii/`. New application and benchmark tests must consume the
same public recipe and runner boundaries as users; they must not construct
private lowerers to make a test pass.

## Reserved future locations

The `parameters/` tree is reserved for versioned experiment inputs. Numerical
experiment families use Deal.II-style `.prm` files; reusable post-processing
styles use versioned JSON profiles beside them. Their current schema and precedence rules are defined by the
[parameter-file reference](../../reference/parameter-files.md) and the
[parameter/plotting-profile decision](parameter-and-plotting-profiles.md);
roadmaps retain the development history behind those interfaces. A parameter file
may describe a matrix of concrete artifacts; generated runs still record the
resolved combination and the hashes of both input documents.

Generated run sets remain under the ignored `runs/` directory. They are
evidence inputs for reports and reproduction investigations, not committed
documentation or source authority.

The headless runner, parameter-file format, and any GUI/client remain
application-layer concerns and must not become core-library dependencies.
