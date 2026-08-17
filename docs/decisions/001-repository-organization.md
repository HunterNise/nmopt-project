# Repository organization for applications

## Decision

Organize the repository by architectural layer and by user-facing
application role. Keep the existing core include layers stable, add an
explicit application API boundary, and keep reusable Chapter 5 recipes
separate from frozen Chapter 6 benchmark scenarios.

The planned structure is:

```text
docs/
  design/                    # normative architecture and interfaces
  reference/                 # exact public API and option reference
  guides/                    # implementation and authoring procedures
  applications/
    chapter-5/               # concrete recipe contracts
    chapter-6/               # application and scenario contracts
  benchmarks/
    chapter-6/               # frozen benchmark specifications
  implementation/             # realized v0/v1 behavior
  planning/                  # roadmap, reviews, and handoffs
  case-studies/              # mathematical derivations
  decisions/                 # accepted architectural decisions
  drafts/                    # explicitly non-authoritative material

include/nmopt/
  application/                # public recipe, scenario, and catalog API
  contract/                   # backend-neutral executable/formulation ports
  semantic/v1/                # ProblemSpec and semantic validation
  compiler/v1/                # lowering and compiled products
  dealii/                     # deal.II backend services
  solvers/                    # optimization and Krylov policies
  experiment/                 # provenance and report envelopes
  reference/                  # dense/reference models

tests/
  support/                    # shared test utilities
  contract/                   # backend-neutral contract tests
  semantic/                   # semantic graph and validation tests
  dealii/                     # deal.II compiler and lowering tests
  application/                # recipe and scenario contract tests
  benchmark/                  # harness and benchmark tests

apps/
  nmopt-runner/               # headless orchestration only

parameters/
  chapter-5/                  # exploratory/default inputs
  chapter-6/                  # versioned benchmark inputs

runs/                         # generated output; ignored by Git
```

Directories are created when their first authoritative file or implementation
unit lands. Empty placeholder directories are not tracked.

## Authority boundaries

- `docs/design/` owns long-lived semantic and architectural contracts.
- `docs/reference/` owns exact public types, options, parameters, and
  diagnostics.
- `docs/guides/` explains how to perform an implementation or authoring task.
- `docs/applications/` explains how a concrete Chapter 5/6 application is
  assembled from the public API.
- `docs/benchmarks/` owns frozen numerical scenario definitions and expected
  evidence.
- `docs/implementation/` records realized capabilities and exclusions.
- `docs/planning/` remains the sole owner of mutable status and handoff order.
- `docs/drafts/` is discussion material and is not an implementation
  authority.

The mathematical Chapter 5 and Chapter 6 guides remain the source catalogue.
Application documents add executable recipe information without duplicating
the mathematical authority. The implementation roadmap remains the status
ledger; application and benchmark documents must not create a second status
table.

## Include-layer policy

The current `include/nmopt/{contract,semantic,compiler,dealii,solvers,experiment,reference}`
layout already follows the framework's layer boundaries and will not be
renamed as part of this restructuring.

New public recipe and scenario records belong under
`include/nmopt/application/`. They may refer to semantic and compiler ports,
but must not move backend details into `ProblemSpec` or add PDE-specific
branches to generic solvers.

## Test-layer policy

Existing test executables retain their names and scenario inventories while
their source files are grouped by the contract layer they verify. CMake
registration remains the sole build integration change for that relocation.
New application and benchmark tests must consume the same public recipe and
runner boundaries as users; they must not construct private lowerers to make a
test pass.

## Deferred application concerns

The headless runner, parameter-file format, generated run artifacts, and GUI
are separate application-layer work. This decision reserves their locations
but does not make a GUI or a serialization format a core-library dependency.
