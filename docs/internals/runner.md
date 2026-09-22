# Repository runner implementation

## Purpose and boundary

This document describes the current implementation of the repository's configured
benchmark runner.

There are two different things called a runner in the source tree:

1. the reusable `HeadlessBenchmarkRunnerT` under `include/nmopt/application/`;
2. the concrete `nmopt_runner` application under `apps/nmopt-runner/`.

They should not be conflated.

For command-line and parameter-file usage, see
[Application execution](../reference/application-execution.md) and
[Parameter files](../reference/parameter-files.md). This page explains ownership,
resolution, lifecycle, and evidence flow inside the implementation.

## Two runner layers

### Reusable benchmark orchestration

`include/nmopt/application/runner.hpp` defines:

```text
HeadlessBenchmarkRunnerT<Scenario>
```

It owns orchestration only:

```text
Scenario
   │
   ├─────────────────────┐
   ▼                     │
problem builder          │
   │                     │
   ▼                     │
producer-specific object │
   │                     │
   └──────────┬──────────┘
              ▼
execution adapter(object, Scenario)
              │
              ▼
   BenchmarkExecutionEvidenceT
              │
              ▼
      BenchmarkHarnessT
              │
              ▼
          artifact text
```

The builder return type is inferred from the callable passed to `run()`. The generic
runner does **not** require the builder to return `ProblemSpec`.

The current Chapter-6 application passes builders that return semantic
`ProblemSpec`s and deal.II execution adapters that compile them. That is a property of
B1/B2, not a restriction in `HeadlessBenchmarkRunnerT`.

The generic runner also optionally measures wall time according to the scenario's
harness policy, then asks `BenchmarkHarnessT` and `BenchmarkArtifactWriter` to
finalize the evidence.

The execution adapter returns `BenchmarkExecutionEvidenceT<Envelope>`, containing:

```text
envelope
ValidationReport
BenchmarkMeasurements
selected field names
additional ArtifactField values
```

`HeadlessBenchmarkRunnerT::run()` converts that into
`BenchmarkRunResultT<Envelope>{artifact, document}`. The typed `BenchmarkArtifactT`
remains available in memory, while `document` is the deterministic flat rendering.
The reusable runner does not choose filesystem paths or write the document itself.

### Concrete repository application

`apps/nmopt-runner/` owns:

- CLI parsing;
- tracked `.prm` configuration;
- matrix expansion and selection;
- reproduction/development run-directory policy;
- run-set manifest lifecycle;
- B1/B2 benchmark registration;
- Chapter-6 scenario binding and execution;
- artifact coordinates;
- solver-trace and native-output placement.

This layer is intentionally project-specific.

## Source map

The application is physically divided as follows.

`run_set_plan.hpp` owns parameter axes, resolved combinations, comparison coordinates,
artifact-coordinate components, provenance, plan validation, and artifact paths.

`run_lifecycle.hpp` owns command-line options, run-kind policy, resolved run
configuration, output-directory allocation, artifact path helpers, and
`RunSetManifest`.

`parameter_files.hpp` owns schema selection/discovery, deal.II `ParameterHandler`
declaration/parsing, typed text conversion helpers, matrix expansion, exclusion
handling, and `make_run_set_plan()`.

`parameter_binding.hpp` maps parsed configuration values into typed scenario options.

`benchmark_binders.hpp` contains the B1/B2-specific scenario binding functions.

`capability_registry.hpp` maps file-facing IDs such as solver method or compilation
product names to typed application selections. It owns name resolution only, not
construction or execution.

`benchmark_registry.hpp` owns stable CLI/parameter benchmark metadata and the generic
execution-registration record types.

`chapter6_execution.hpp` owns B1/B2 artifact-coordinate policies, environment and
artifact-field enrichment, concrete execution loops, and the execution registrations.

`main.cc` remains the entry point and owns run preparation, configuration snapshotting,
top-level dispatch, and catalog listing.

## End-to-end repository run

```text
CLI
 │
 ▼
parse_command_line()
 │
 ▼
load and validate ParameterFile
 │
 ▼
make_run_set_plan()
 │
 ▼
resolve_run_configuration()
 │
 ├── snapshot parameters / plotting profile / resolved combinations
 │
 ▼
RunSetManifest(status = running)
 │
 ▼
BenchmarkExecutionRegistration::execute
 │
 ├── bind one typed scenario per resolved combination
 ├── construct runtime data + owned compilation session
 ├── invoke HeadlessBenchmarkRunnerT through the Chapter-6 adapter
 ├── write artifact.kv / solver-trace.csv / optional native output
 └── record success or failure in the run manifest
 │
 ▼
RunSetManifest::finalize()
 │
 ▼
process exit status
```

Post-processing consumes persisted run evidence afterward; it is not another execution
branch inside `HeadlessBenchmarkRunnerT`.

## Parameter-file resolution

The runner does not feed arbitrary strings directly into scenario fields.

The parse pipeline is:

```text
raw .prm text
     │
     ▼
benchmark + recipe identity
     │
     ▼
ParameterSchemaAdapter
     │
     ▼
dynamic scalar-definition discovery
     │
     ▼
declare ParameterHandler schema
     │
     ▼
parse and extract values
     │
     ▼
ParameterFile
     │
     ▼
matrix expansion + selection + exclusions
     │
     ▼
RunSetPlan
```

### Schema selection

`ParameterSchemaAdapter` contains:

- accepted benchmark IDs/prefixes;
- accepted recipe IDs;
- ordinary schema entries;
- scalar-definition schemas.

The runner first reads enough raw text to determine benchmark/recipe identity, selects
the matching adapter, discovers dynamic scalar definition entries, then declares the
complete deal.II `ParameterHandler` schema before parsing.

This is why scalar function catalogs can be data-driven without making the main schema
an untyped free-form map.

### `ParameterFile`

After parsing, `ParameterFile` stores:

- source path;
- deterministic content hash;
- resolved key/value strings;
- matrix axes;
- file-level selections;
- excluded combinations.

Typed helpers convert individual values into numbers, booleans, lists, enum-like
selections, and scalar/vector function definitions.

The content hash is provenance drift detection (`fnv1a64:`), not a cryptographic
authentication mechanism.

## Matrix expansion and `RunSetPlan`

`ParameterFile::combinations()` performs a Cartesian expansion of declared matrix
axes, then applies file and CLI selections and removes explicitly excluded
combinations.

It rejects:

- undeclared axes;
- unknown selected values;
- duplicate selected values;
- malformed/incomplete exclusions;
- a selection that resolves to no combinations.

`make_run_set_plan()` then records:

- benchmark ID;
- original matrix axes;
- effective selections;
- excluded combinations;
- resolved combinations;
- comparison row/column/grouping coordinates;
- parameter-file path and content hash.

Each resolved combination also receives artifact-coordinate components.

The default coordinate policy mirrors matrix axes. B1/B2 can supply their own
coordinate policy when their stable artifact organization is more meaningful than the
generic axis names.

## Registration has two layers

`benchmark_registry.hpp` deliberately separates metadata from execution.

`BenchmarkRegistration` records only:

```text
CLI ID
parameter-file benchmark ID
default parameter file
```

The current metadata registry contains B1 and B2.

`BenchmarkExecutionRegistration` adds:

```text
metadata pointer
artifact planner
execution callback
```

The Chapter-6 execution registrations live in `chapter6_execution.hpp`, not in the
metadata registry.

This separation allows configuration and identity checks to occur without pretending
that every registered piece of metadata automatically has executable support.

`capability_registry.hpp` is a different concept again. It translates individual
file-facing capability IDs to typed values such as reduced method, execution
selection, or compilation product.

## Run preparation in `main.cc`

`prepare_run()` is the boundary between CLI/file input and an executable run set.

Its order is important.

### 1. Load the parameter file

The CLI can specify a parameter file directly, or select a benchmark whose registration
points to a default authoritative file.

Paths are located from the current directory or its parents so invocation is not tied
to one working directory.

### 2. Reconcile file identity and CLI identity

The parameter file's `Benchmark/id` is resolved through the metadata registry.

When the user supplied a parameter file, the file owns its benchmark/run kind and,
unless explicitly overridden, output root. Reproduction files also require their
declared build profile to match the compiled runner profile.

When the user selected `--benchmark`, the selected ID must match the benchmark
registration behind the authoritative parameter file.

### 3. Build the run-set plan

`make_run_set_plan()` applies CLI selections on top of file selections and expands the
matrix.

### 4. Require executable registration and validate artifact coordinates

Before output is created, the runner finds the corresponding
`BenchmarkExecutionRegistration` and invokes its artifact planner.

This validates benchmark-specific axis/coordinate requirements early and guarantees
that expected artifact paths are unique.

### 5. Resolve run configuration

`resolve_run_configuration()` applies the run-kind policy.

Reproduction runs require the release deal.II build profile. Development runs use
numbered or explicitly named slots.

The final `ResolvedRunConfiguration` is enriched with:

- parameter-file path/hash;
- plotting-profile path/hash;
- selected/declaration/exclusion/resolution strings;
- comparison coordinates.

## Configuration snapshot

Before benchmark execution begins, `snapshot_configuration()` writes a reproducible
run-set snapshot:

```text
parameters.prm
plotting-profile.json
resolved-combinations.txt
```

These are copied/resolved inputs for the run, distinct from per-artifact evidence.

## Run-set lifecycle and manifest

`RunSetManifest` is created with the complete expected artifact inventory.

Construction:

1. creates the run directory;
2. rejects empty or duplicate expected artifact paths;
3. records every expected artifact as `pending`;
4. writes `run-manifest.json` immediately with overall status `running`.

During execution:

```text
artifact succeeds  → record_success() → rewrite manifest
artifact fails     → record_failure() → rewrite manifest
```

An artifact may be recorded only once and must belong to the expected inventory.

`finalize()` converts any still-pending artifact to an error with
`artifact was not executed`, writes the final manifest, and returns success only when
the failure count is zero.

The overall manifest state is therefore:

```text
not finalized             → running
finalized, no failures    → complete
finalized, any failure    → failed
```

This makes partial execution visible instead of silently treating missing artifacts as
absent work.

## Chapter-6 execution seam

`chapter6_execution.hpp` is separate from `main.cc` so tests can invoke production
execution behavior directly without including the executable entry point.

It is intentionally **not** a generic framework.

It owns Chapter-6 concerns such as:

- B1/B2 artifact coordinate conventions;
- B1 matched-reference objective-target ordering;
- B2 case slugs;
- environment records;
- manifest and solver metadata projected into artifact fields;
- B1/B2 execution loops;
- solver-trace writing;
- B1/B2 execution registrations.

The file remains substantial because these concerns form one repository-application
responsibility. Splitting it merely by line count would not create a stronger reusable
boundary; the generic runner already lives under `include/nmopt/application/`.

## B1 execution

For each resolved B1 combination, `run_b1()` resolves:

- method;
- regularization value;
- artifact path.

It then:

```text
make B1 scenario
      │
      ▼
bind .prm + combination into typed scenario
      │
      ▼
apply optional mesh override
      │
      ▼
construct runtime data
      │
      ▼
construct owned compilation session
      │
      ▼
construct B1ReducedExecutionAdapterT
      │
      ▼
HeadlessBenchmarkRunnerT::run()
      │
      ▼
artifact.kv + solver-trace.csv + optional native output
```

### Matched objective-target dependency

B1 has one deliberate cross-coordinate dependency.

When `Solver/objective target policy` is `match-reference-method`, combinations are
stable-sorted so the configured reference method executes first for each
regularization coordinate.

Only a reference result that stops for the selected successful stopping criterion is
recorded as an objective target. A dependent method fails that artifact if its
reference coordinate has not produced a valid target.

This dependency belongs in the Chapter-6 application shell, not in the generic
benchmark runner.

## B2 execution

B2 iterates the resolved combinations independently.

For each combination it resolves the observation region and target profile, constructs
the target catalog and scenario, binds parameter values, builds manufactured/runtime
data and an owned compilation session, then invokes
`B2ReducedExecutionAdapterT<2>` through `HeadlessBenchmarkRunnerT`.

The adapter requires the compiler-produced application view because B2 evidence uses
native field output, physical/independent dimensions, and objective-component
decomposition in addition to solver-facing reduced services.

The application view is still a compiler-side optional capability; it is not part of
the generic runner interface.

## Execution adapters

The Chapter-6 deal.II adapters live under:

```text
include/nmopt/application/dealii/chapter6_b1.hpp
include/nmopt/application/dealii/chapter6_b2.hpp
```

They are the boundary between typed scenario data and numerical execution.

A typical adapter:

1. validates the typed scenario;
2. creates `DealiiDataBindings`;
3. maps compile options to `DealiiDiscretisationPolicy`;
4. calls `DealiiCompiler`;
5. checks that the expected product was produced;
6. constructs/evaluates the reduced formulation;
7. runs the selected reduced solver;
8. constructs derivative/Hessian/evidence records as appropriate;
9. writes native output through `CompiledApplicationViewT` when requested and
   available;
10. returns `BenchmarkExecutionEvidenceT`.

The generic `HeadlessBenchmarkRunnerT` does not know these steps.

## Experiment envelope

B1/B2 use `ReducedSearchExperimentEnvelopeT<SerialBackend>` from
`include/nmopt/experiment/reduced_envelope.hpp`.

The envelope copies:

- `CompilationManifest`;
- a solver-policy snapshot;
- the solver report;
- caller-supplied `RunEnvironmentRecord`.

It intentionally does not retain the compiled problem or reduced DTO.

The environment is constructed by the application shell because compiler version,
operating system, architecture, hardware label, and source revision are orchestration
evidence, not numerical-solver responsibilities.

## Artifact layers

A completed run set can contain several different evidence layers.

At the run-set root:

```text
run-manifest.json
parameters.prm
plotting-profile.json
resolved-combinations.txt
```

Per artifact coordinate:

```text
artifact.kv
solver-trace.csv
native/...
```

`artifact.kv` is rendered from `BenchmarkArtifactT` plus typed artifact fields.

`solver-trace.csv` is application-specific projection of line-search trial evidence.

`native/` is written by compiler/application-specific output code; it is not an
`ExecutableModelT` contract.

## Post-processing boundary

The runner writes evidence. It does not own plotting/report interpretation.

After execution:

```text
run evidence
    │
    ▼
tools/nmopt_artifacts
    │
    ├── dependency-free persisted-record semantics
    │
    ▼
tools/nmopt_postprocess
    │
    ▼
plots / comparisons / reports
```

Chapter-specific reporting remains allowed above the shared record parser. Generic
parsing/rendering should not reintroduce hard-coded B1/B2 policy where the tracked
plotting profile owns that policy.

## What the runner does not prove

The repository runner is a first-class consumer of the semantic/compiler path and
Chapter-6 application adapters.

It should not be described as a universal experiment shell already demonstrated for
every producer.

In particular:

- `HeadlessBenchmarkRunnerT` itself is generic;
- the surrounding `.prm`, benchmark registry, B1/B2 execution registration, run-set
  coordinates, and artifact enrichment are repository-specific;
- the external/native Step-4 producer is a first-class numerical/optimization
  consumer, but it is not currently demonstrated as a consumer of this full
  scenario/run-set infrastructure.

This distinction keeps the reusable orchestration seam separate from the concrete
reproduction application.

## Tests

Backend-neutral runner/lifecycle behavior is exercised in:

```text
tests/application/runner_contract.cc
```

It covers:

- benchmark metadata registry behavior;
- reproduction/development policy;
- run-directory allocation;
- typed capability registries;
- CLI rejection;
- matched-reference stopping classification;
- run-manifest running/success/failure state.

Parameter/schema/binding and production execution behavior are characterized in:

```text
tests/application/parameter_files_dealii_contract.cc
```

This test calls the production `chapter6_execution.hpp` seam rather than including
`main.cc`. That separation is the reason the execution header exists.

B1/B2 numerical application behavior is further covered by:

```text
tests/application/chapter6_b1_dealii_contract.cc
tests/application/chapter6_b2_dealii_contract.cc
tests/application/compiled_application_view_dealii_contract.cc
```

## Extending the runner

A new repository benchmark needs two separate registrations:

1. metadata in `benchmark_registry.hpp`;
2. artifact planning and execution in the application execution registry.

If the benchmark is configured by `.prm`, add or extend the appropriate
`ParameterSchemaAdapter` and typed scenario binder.

Keep generic and benchmark-specific concerns separated:

```text
generic:
    CLI/run lifecycle
    matrix/run-set model
    benchmark harness
    generic headless runner

benchmark-specific:
    scenario schema meaning
    typed binding
    coordinate convention
    cross-coordinate dependencies
    execution adapter
    evidence fields
```

Do not put a new benchmark ID into `HeadlessBenchmarkRunnerT`, and do not make the
generic persisted-evidence parser depend on Chapter-specific fields.

If a future external/application-owned producer needs the full run-set shell, add an
execution registration and adapter for that producer rather than routing it through
`ProblemSpec` solely to reuse the runner.
