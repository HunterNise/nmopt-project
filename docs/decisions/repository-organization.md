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
  design/                         # long-lived architecture and mathematics
  reference/                      # exact public API and execution references
  guides/                         # implementation and source-example guides
  applications/                   # Chapter 5/6 application contracts
  benchmarks/                     # frozen benchmark contracts
  implementation/                 # realized capabilities and selected policies
  planning/                       # mutable roadmaps, reviews, and handoffs
  case-studies/                   # worked mathematical/interface studies
  decisions/                      # accepted repository and architecture decisions
  drafts/                         # explicitly non-authoritative material
include/nmopt/
  application/                    # public recipe, scenario, and catalog API
  contract/                       # backend-neutral executable/formulation ports
  semantic/v1/                    # ProblemSpec and semantic validation
  compiler/v1/                    # lowering and compiled products
  dealii/                         # deal.II backend services
  solvers/                        # optimization and Krylov policies
  experiment/                     # provenance and report envelopes
  reference/                      # dense/reference models
tests/
  support/                        # shared test utilities
  contract/                       # backend-neutral contract tests
  semantic/                       # semantic graph and validation tests
  dealii/                         # deal.II compiler and lowering tests
  application/                    # recipe and scenario contract tests
  benchmark/                      # harness and benchmark tests
apps/
  nmopt-runner/                   # headless application orchestration
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
- `docs/guides/` explains how to perform implementation or source-example
  tasks; the Chapter 6 numerical-examples guide is authoritative for what the
  book states.
- `docs/applications/` explains how concrete Chapter 5/6 applications are
  assembled from the public API.
- `docs/benchmarks/` owns frozen numerical scenario definitions and required
  evidence.
- `docs/implementation/` records realized capabilities, exclusions, and
  selected implementation policies.
- `docs/planning/` owns mutable status, work order, review evidence, and
  handoffs.
- `docs/case-studies/` contains worked derivations and interface studies. It
  is explanatory material, not a second normative contract.
- `docs/decisions/` contains accepted architectural and repository decisions.
- `docs/drafts/` contains discussion material that is explicitly not
  authoritative.

The mathematical Chapter 5 and Chapter 6 guides remain the source catalogue.
Application and benchmark documents add executable choices and evidence
requirements without silently rewriting the source catalogue. The
implementation and application roadmaps remain the mutable status ledgers for
their respective layers.

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
| `apps/` | Source files for executable products. | `apps/nmopt-runner/main.cc` is the headless entry point; its adjacent `runner.hpp` is private CLI/run-set support. The public backend-neutral runner API is `include/nmopt/application/runner.hpp`. |
| `tests/` | Contract-layer test translation units and shared test support. | CTest inventories are discovered from test executables; generated registrations belong in `build/`, not in source. |
| `tools/` | Repository-local user and analysis tools, including Python post-processing and shell wrappers. | Tools consume public artifacts and interfaces; they must not become a second PDE lowerer or solver implementation. |
| `cmake/` | Source-controlled CMake modules used by the root build. | Build configuration remains expressed through root `CMakeLists.txt` and `CMakePresets.json`; generated CMake state does not belong here. |
| `build/` | Ignored generated build trees, one named directory per configured preset. | Never configure directly in the `build/` container; use `build/debug-neutral/`, `build/debug-dealii/`, `build/sanitize-neutral/`, or `build/release-dealii/`. |
| `runs/` | Ignored generated experiment evidence. | Run manifests, artifacts, native output, reports, and post-processing follow the [execution reference](../reference/application-execution.md). |

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

The [build instructions](../../.agents/build.md) own profile selection,
cache recovery, generated-file rules, manual-command safety rules, and build
commands. The [application
execution reference](../reference/application-execution.md) owns artifact
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

The `reference/` directory intentionally contains the two agent-facing
application references:

- `application-api.md` — assembly, recipe, scenario, compiler, solver, and
  provenance boundaries;
- `application-execution.md` — schemas, run organization, native output,
  reports, post-processing, and verification commands.

There is no additional directory README because `docs/README.md` is the
repository-wide documentation map and these two files are self-describing.

## Include and test layers

The current `include/nmopt/{contract,semantic,compiler,dealii,solvers,experiment,reference}`
layout follows the framework's layer boundaries and is not renamed by this
decision. New public recipe and scenario records belong under
`include/nmopt/application/`; they may refer to semantic and compiler ports,
but must not move backend details into `ProblemSpec` or add PDE-specific
branches to generic solvers.

Existing test executables retain their names and scenario inventories while
source files are grouped by the contract layer they verify. New application
and benchmark tests must consume the same public recipe and runner boundaries
as users; they must not construct private lowerers to make a test pass.

## Reserved future locations

The `parameters/` tree is reserved for versioned experiment inputs. Numerical
experiment families use Deal.II-style `.prm` files; reusable post-processing
styles use versioned JSON profiles beside them. Their schema and precedence
rules are defined by the application roadmap's parameter-file unit and its
[parameter-file reference](../reference/parameter-files.md). A parameter file
may describe a matrix of concrete artifacts; generated runs still record the
resolved combination and the hashes of both input documents.

Generated run sets remain under the ignored `runs/` directory. They are
evidence inputs for reports and reproduction investigations, not committed
documentation or source authority.

The headless runner, parameter-file format, and any GUI/client remain
application-layer concerns and must not become core-library dependencies.
