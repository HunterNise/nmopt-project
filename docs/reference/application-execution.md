# Application execution reference

This reference explains how to **run an already-authored nmopt application**,
retain deterministic evidence, organize a run set, and post-process persisted
outputs.

If you are adding a new recipe/scenario/backend adapter, start with
[Application authoring](application-authoring.md). If you are configuring an
existing application through `.prm`, use
[Parameter files](parameter-files.md).

The reusable execution boundary is:

```text
typed ScenarioT
      │
      ▼
HeadlessBenchmarkRunnerT
      │
      ├── ProblemSpec builder
      └── execution adapter
              │
              ▼
BenchmarkExecutionEvidenceT
      │
      ▼
BenchmarkArtifactT
      │
      ▼
BenchmarkArtifactWriter
      │
      ▼
artifact.kv
```

The repository application `nmopt_runner` then adds matrix expansion,
filesystem run sets, manifests, configuration snapshots, and registered
benchmark dispatch.

## What an authored application should already provide

Before the execution layer starts, the application family should expose:

- a typed scenario factory;
- scenario validation;
- a `ProblemSpec` builder from scenario problem parameters;
- a backend execution adapter returning `BenchmarkExecutionEvidenceT`;
- optional catalog metadata;
- optional runner/parameter registration.

For example, the current Chapter 6 path has:

```text
make_b1_scenario(...)
make_b1_problem_spec(...)
B1ReducedExecutionAdapterT<2>
```

and the equivalent B2 pieces.

Execution code should compose those objects rather than duplicate their
problem/compiler/solver logic.

## Public execution types

The reusable application headers are:

```cpp
#include "nmopt/application/harness.hpp"
#include "nmopt/application/artifact_writer.hpp"
#include "nmopt/application/runner.hpp"
```

The main public types are:

```cpp
BenchmarkIdentity
BenchmarkMeasurements
BenchmarkArtifactT<Envelope>
BenchmarkExecutionEvidenceT<Envelope>
BenchmarkRunResultT<Envelope>
BenchmarkHarnessT<Scenario>
HeadlessBenchmarkRunnerT<Scenario>
```

They live under:

```cpp
nmopt::application::benchmark
```

None of these types knows deal.II or a particular optimizer.

## Working path: execute one scenario from C++

Assume the application author already supplied:

```cpp
MyScenario make_my_scenario();

nmopt::semantic::v1::ProblemSpec
make_my_problem_spec(
  const MyScenario::problem_parameters_type &);

MyExecutionAdapter execute;
```

### Step 1 – obtain and adjust the typed scenario

```cpp
MyScenario scenario = make_my_scenario();

scenario.solver.parameters.gradient_tolerance = 1.0e-8;
scenario.experiment.source_revision = source_revision;
scenario.experiment.build_profile = build_profile;

validate_my_scenario(scenario);
```

Application users should adjust declared typed fields. They should not need to
reconstruct the `ProblemSpec`, deal.II data bindings, or optimizer template
composition.

### Step 2 – construct the headless runner

```cpp
using Runner =
  nmopt::application::benchmark::
    HeadlessBenchmarkRunnerT<MyScenario>;

Runner runner(scenario);
```

`runner.identity()` is available before execution and is derived from the
scenario's metadata/experiment record.

### Step 3 – supply the problem builder

```cpp
auto build_problem =
  [](const auto &problem_parameters) {
    return make_my_problem_spec(
      problem_parameters);
  };
```

For the Chapter 6 semantic/compiler path, this callback does one thing:

```text
typed problem parameters
        ↓
ProblemSpec
```

More generally, the runner contract is:

```text
typed problem parameters
        ↓
application-defined build context
```

The runner only forwards that result to the execution adapter. Keep expensive
backend realization, solving, and output in the adapter rather than hiding
them in the builder.

### Step 4 – supply the execution adapter

If the backend adapter is already callable as

```cpp
Evidence operator()(
  const ProblemSpec &,
  const MyScenario &) const;
```

the execution callback can simply forward:

```cpp
auto run = runner.run(
  build_problem,
  [&](const auto &specification,
      const auto &run_scenario) {
    return execute(
      specification,
      run_scenario);
  });
```

The runner measures wall time around problem construction + execution when

```cpp
scenario.experiment.harness.measure_timings
```

is true.

### Step 5 – consume the detached result

The return type is:

```cpp
BenchmarkRunResultT<Envelope>
```

with:

```cpp
run.artifact;
run.document;
```

`run.artifact` owns the typed detached envelope and generic evidence.
`run.document` is the deterministic flat rendering produced by
`BenchmarkArtifactWriter`.

At this point the caller can write the document wherever its application
policy requires:

```cpp
std::ofstream output(path);
output << run.document;
```

Filesystem placement is intentionally outside `HeadlessBenchmarkRunnerT`.

## What execution evidence should retain

The execution adapter returns:

```cpp
template <typename Envelope>
struct BenchmarkExecutionEvidenceT
{
  Envelope                        envelope;
  semantic::v1::ValidationReport  diagnostics;
  BenchmarkMeasurements           measurements;
  std::vector<std::string>        selected_fields;
  std::vector<ArtifactField>      fields;
};
```

The fields have different purposes.

### `envelope`

This is the detached typed numerical record. For a reduced compiler-backed
application it typically owns copies of:

```text
compilation manifest
solver policy snapshot
solver result/report
run environment
```

It should remain usable after mesh/session/solver objects have been destroyed.

### `diagnostics`

Retain semantic/compiler diagnostics rather than reducing failure to a boolean.
The artifact writer serializes the diagnostic category, component, capability,
and remedy.

### `measurements`

The reusable runner can fill wall time. Memory evidence, when collected, must
be supplied by application-specific execution code.

### `selected_fields`

This is the stable logical inventory of fields the application chose to
retain, for example:

```cpp
evidence.selected_fields = {
  "state",
  "control",
  "adjoint"
};
```

It is metadata; the numerical field arrays live in backend-native files.

### `fields`

These are additional deterministic flat artifact fields:

```cpp
evidence.fields.push_back({
  "benchmark.method",
  "l-bfgs"
});
```

Use them for application evidence that belongs in `artifact.kv` but is not
part of the generic benchmark identity/diagnostic schema.

## `artifact.kv`

`BenchmarkArtifactWriter` always emits the generic groups:

```text
artifact.schema
identity.*
diagnostics.*
measurements.*
selected_field[*]
```

and then adds the application-specific `ArtifactField` values in sorted key
order.

Values are escaped deterministically and numeric formatting uses the classic
locale with precision 17.

This file is deliberately flat. It is useful for stable contract tests,
command-line inspection, and lightweight Python loading.

It is **not** the same object as `run-manifest.json`:

```text
artifact.kv
    one executed artifact / one resolved scenario coordinate

run-manifest.json
    lifecycle and inventory of the whole run set
```

## Native numerical output

The application/backend execution adapter remains responsible for authoritative
field output.

The current Chapter 6 artifact shape is:

```text
artifact-directory/
├── artifact.kv
├── solver-trace.csv
└── native/
    ├── mesh-volume.vtu
    ├── mesh-volume.svg
    ├── fields-volume.vtu
    └── control-boundary.vtu   # when boundary topology is separate
```

Not every application needs every file.

The rule is that post-processing consumes persisted numerical outputs. It
should not re-run the PDE or reconstruct an authoritative FE field that was
never written by the execution adapter.

## The repository `nmopt_runner`

`apps/nmopt-runner` builds a concrete run-set application on top of the
reusable execution API.

Its high-level flow is:

```text
CLI
 │
 ▼
benchmark registration
 │
 ▼
ParameterFile
 │
 ▼
RunSetPlan
 │
 ▼
resolved ScenarioT for each coordinate
 │
 ▼
HeadlessBenchmarkRunnerT
 │
 ▼
artifact.kv + native files
 │
 ▼
RunSetManifest
```

The current runnable benchmark registrations are:

```text
b1
b2
```

Listing application metadata is separate from runnable benchmark registration:

```bash
build/debug-dealii/bin/nmopt_runner --list
```

A catalog entry can exist without an execution registration.

## Run through the default registered benchmark

The compatibility CLI path is:

```bash
build/release-dealii/bin/nmopt_runner \
  --benchmark b1 \
  --framework-revision REV
```

A benchmark registration maps the short ID to its default authoritative
parameter file.

The parameter file is still the experiment-family definition; `--benchmark`
only selects that registered default.

For explicit configuration:

```bash
build/release-dealii/bin/nmopt_runner \
  --parameter-file parameters/chapter-6/b1/authoritative.prm \
  --framework-revision REV
```

Do not supply `--benchmark` and `--parameter-file` together.

Every run requires `--framework-revision`; the runner persists it in the run
manifest and artifact fields.

## Run-set planning before execution

A parsed parameter file is expanded into a `RunSetPlan` containing:

```text
benchmark ID
matrix axes
in-file selection
CLI selection
excluded combinations
resolved combinations
comparison coordinates
parameter-file provenance
```

For each resolved combination, the plan also retains ordered artifact
coordinate components.

This plan is validated before the runner starts creating artifact output.

For example, a family with:

```text
method = steepest-descent, l-bfgs
regularisation = 1e-1, 1e-2
```

resolves to four coordinates before selection/exclusion:

```text
[method=steepest-descent, regularisation=1e-1]
[method=steepest-descent, regularisation=1e-2]
[method=l-bfgs,           regularisation=1e-1]
[method=l-bfgs,           regularisation=1e-2]
```

The benchmark's artifact planner maps those logical coordinates to path
components.

The details of matrix syntax and precedence are in
[Parameter files](parameter-files.md).

## Run kinds and output directories

The current runner distinguishes:

```text
reproduction
development
```

### Reproduction

A reproduction run uses:

```text
<output>/chapter-6/<benchmark>/authoritative/
```

and current policy requires a runner compiled as:

```text
release-dealii
```

For parameter-file runs, the declared build profile must match the compiled
runner profile.

A debug or refinement-smoke result should therefore be recorded as development
rather than represented as source-sized reproduction evidence.

### Development

Development run directories are allocated monotonically:

```text
<output>/chapter-6/<benchmark>/development/001/
<output>/chapter-6/<benchmark>/development/002/
...
```

A parameter-file development run can request a new named slot:

```bash
nmopt_runner \
  --parameter-file parameters/chapter-6/b1/development/my-family.prm \
  --framework-revision REV \
  --run-slot 004-v2
```

The slot must be a single unused directory name. It does not overwrite an
existing run.

## CLI overrides

The runner's generic options are:

```text
--list
--benchmark ID
--parameter-file FILE
--select AXIS=VALUE
--output DIRECTORY
--run-kind reproduction|development
--run-slot SLOT
--framework-revision REV
--refinement N
```

The ownership rules are more useful than memorizing the options:

```text
parameter file
    owns numerical/application family and, for --parameter-file, run kind

--select
    narrows declared matrix axes

--output
    changes destination only

--run-slot
    changes development placement only

--refinement
    explicit supported mesh override for smoke/development workflows
```

A parameter-file invocation cannot override its `Run/kind` through
`--run-kind`.

For the exact current CLI surface, use:

```bash
nmopt_runner --help
```

## Run-manifest lifecycle

`RunSetManifest` writes `run-manifest.json` when the run set is created, before
individual artifacts finish.

The initial structure is:

```text
run status      = running
artifact status = pending
```

Each artifact is then recorded once as:

```text
ok
```

or:

```text
error
```

with an error message.

When `finalize()` runs, any still-pending artifact becomes an error with:

```text
artifact was not executed
```

The final run status is:

```text
complete   when every artifact is ok
failed     when one or more artifacts failed
```

This design preserves the expected artifact inventory and partial-failure
evidence instead of only writing a manifest after a completely successful
run.

The manifest retains:

```text
benchmark and run kind
build profile and framework revision
run directory
refinement override
parameter file/hash
selection and declared/resolved matrix
excluded combinations
plotting profile/hash
comparison rows/columns/group_by
original command
artifact inventory/status/error
```

## Configuration snapshots

Before benchmark execution, the current runner copies:

```text
parameters.prm
plotting-profile.json
resolved-combinations.txt
```

into the run root.

The run manifest also keeps hashes/provenance for the parameter and plotting
sources.

Those snapshots describe one historical execution. They are not the editable
source for a future run; tracked configuration remains under `parameters/`.

## Cross-coordinate execution dependencies

Most resolved matrix coordinates can execute independently.

The current B1 runner also supports a matched-objective-target policy in which
one method supplies a target value for another method at the same compatible
coordinate.

The runner explicitly orders the reference method first and only accepts its
objective when its stopping reason satisfies the configured criterion.

This dependency is application/run-set policy. It does not live in
`ReducedSearchSolverT` and should not be reproduced by relying on filesystem
iteration order.

## Post-process one persisted artifact

```bash
python3 tools/postprocess.py \
  --artifact RUN/artifacts/... \
  --output RUN/artifacts/.../postprocess
```

For an entire run root:

```bash
python3 tools/postprocess.py \
  --input RUN \
  --output RUN/postprocess
```

The postprocessor searches the run root/parents for:

```text
parameters.prm
plotting-profile.json
run-manifest.json
```

and reconstructs the effective matrix/comparison/style policy from the
persisted snapshot.

An explicit:

```bash
--profile-file FILE
```

overrides the snapshot style for that invocation, while:

```bash
--format png svg
```

overrides output formats.

The derived-output provenance records both the snapshot sources and explicit
overrides. It does not modify the original artifact or run manifest.

## Reports and plots consume the same persisted run differently

The Chapter 6 report consumes the run manifest:

```bash
python3 tools/chapter6_report.py \
  --run-manifest RUN/run-manifest.json \
  --output RUN/report
```

and writes summary tables such as:

```text
summary.csv
summary.md
```

It includes failed/pending/missing artifact status rather than silently
dropping unsuccessful coordinates.

`tools/postprocess.py` consumes persisted artifact metadata plus native field
and history files to render figures.

So the data flow is:

```text
C++ execution
   │
   ├── run-manifest.json
   │       └── report summaries
   │
   └── artifact.kv + native fields + solver trace
           └── field/history post-processing
```

## End-to-end view

Once an application family is authored, a normal configured run is:

```text
tracked .prm + plotting JSON
          │
          ▼
nmopt_runner
          │
          ├── parse / resolve RunSetPlan
          └── create run root + manifest
          │
          │  for each resolved coordinate
          ▼
typed ScenarioT
          │
          ▼
HeadlessBenchmarkRunnerT
          │
          ├── ProblemSpec builder
          └── backend execution adapter
                  │
                  ▼
          detached artifact evidence
                  │
                  ├── artifact.kv
                  ├── solver trace
                  └── native fields
          │
          ▼
finalized run-manifest.json
          │
          ├── report
          └── postprocess
```

That separation is the intended user-facing boundary: application authors
implement the typed scenario and execution adapter once; application users can
then operate the family through `.prm`, CLI filtering, persisted run evidence,
and post-processing without touching the semantic compiler or optimizer
assembly code.
