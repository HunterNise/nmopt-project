# Draft: application, run, and GUI layer

Status: working draft for later refinement. This document records an
application-layer and user-interface direction; it is not an architectural
contract and is not part of the implementation roadmap.

## Intent

The framework should eventually support two complementary ways of using a
compiled optimal-control problem:

1. a reproducible, headless command-line run suitable for tests, batch jobs,
   clusters, and scripts; and
2. a Slint-based desktop layer that discovers available targets, scenarios,
   and parameter files, lets a user assemble a valid run, streams the command
   line output, and opens the produced artifacts for inspection.

The GUI should be a client of the same headless runner used by scripts. It
should not contain a second solver path, silently construct a different
problem, or make the core compiler depend on Slint.

## Current repository starting point

The repository currently has reusable library and compiler code plus test
executables, but not a user-facing application catalog or a GUI application.
The existing scenario-dispatch pattern is a useful starting point: an
executable can expose named scenarios and list them, while the scenario
builder constructs the semantic problem. This pattern should be turned into a
small, explicit application contract before a GUI relies on it.

The framework should preserve the composition principle. A problem should be
assembled from operators, controls, observations, metrics, constraints, and
discretization components. An application source file should describe a
scenario or a target-family adapter, not become a new monolithic PDE problem
class for every combination.

## Proposed separation of responsibilities

### Core library and compiler

The core owns semantic validation, lowering, compiled problems, solver
interfaces, and structured reports. It must remain usable without a display,
filesystem layout, JSON library, or Slint. Generic solvers should report
events and typed data through an observer or sink; they should not write VTU,
CSV, or GUI-specific files themselves.

### Headless runner

The runner is the stable user-facing executable. It should:

- discover and describe registered targets and scenarios;
- load one selected parameter file;
- validate the selection before allocating an expensive solve;
- compile the semantic problem;
- run the selected solver;
- emit human-readable progress for terminals;
- optionally emit machine-readable events, preferably JSON Lines;
- create an immutable run directory containing inputs, metadata, reports,
  histories, and visualization files; and
- return meaningful exit codes for invalid input, compilation failure,
  numerical failure, cancellation, and successful completion.

The runner should be usable independently of the GUI. A first implementation
could be a target-family executable, but the long-term interface should avoid
requiring one binary for every individual problem combination.

### Slint GUI

The GUI should discover capabilities by asking the runner for a catalog or by
reading a generated catalog with the same schema. It should present only
options that the runner can validate. It should launch the runner as a child
process, capture stdout/stderr and structured events, and display the same
run directory after completion.

The GUI is therefore an orchestration and inspection layer, not a new
numerical layer. This keeps headless execution available on systems where
Slint is not built and makes failures reproducible outside the GUI.

## Repository organization

A possible future layout is:

```text
apps/
  nmopt-runner/                 # stable headless entry point
    main.cc
    catalog.cc
    run_controller.cc
  nmopt-gui/                    # optional Slint application
    main.rs                     # or the selected Slint host language
    runner_client.rs
    model.rs
    ui.slint

scenarios/
  boundary_control/
    catalog.cc                  # target/scenario registration
    builders.cc                 # semantic graph construction
    builders.hpp

parameters/
  boundary_control/
    weighted_neumann_tracking.prm
    partial_dirichlet_tracking.prm
    README.md

runs/                            # local generated output; normally ignored
  boundary_control/
    <scenario>/<timestamp>-<run-id>/
```

This is a direction, not a request to create all of these directories now.
The names can change once the application boundary and build packaging are
settled.

### Source files and problem types

Do not write one source file for every mathematical combination. Prefer:

- one source module per reusable scenario family or target family;
- small registration records that name the target and expose its builder;
- parameter files for values and policies that should change between runs;
- tests for semantic combinations and expected failure boundaries; and
- explicit source code when a new composition is genuinely a new supported
  policy rather than merely a new number.

For example, changing a mesh, coefficient, desired state, regularization
weight, solver tolerance, or output cadence belongs in a parameter file. A
new partial-boundary ownership policy, observation type, or unsupported
composition belongs in the typed semantic/compiler layer and should not be
hidden in a `.prm` file.

### Target, scenario, and parameter file

These should be distinct concepts:

- A **target** is a registered, validated capability with a stable identifier
  and declared requirements.
- A **scenario** selects a concrete composition and supplies the semantic
  choices needed to build that target.
- A **parameter file** supplies run-time data and numerical policy for that
  scenario.

The GUI should display these as separate selections, while the runner should
record the fully resolved combination in the run metadata. A parameter file
should not be able to turn an unsupported target into a supported one by
using an undocumented key.

## Parameter-file and discovery model

The first version can continue to use deal.II `ParameterHandler` files. Each
registered scenario should expose a schema or a validated default parameter
tree. The runner should support at least:

```text
nmopt-runner --list-targets
nmopt-runner --describe-target boundary_control.partial_dirichlet
nmopt-runner --validate --target boundary_control.partial_dirichlet \
  --prm parameters/boundary_control/partial_dirichlet_tracking.prm
nmopt-runner --run --target boundary_control.partial_dirichlet \
  --prm parameters/boundary_control/partial_dirichlet_tracking.prm \
  --output runs/
```

The exact command names are open. The important properties are that listing
and validation are cheap, non-solving operations; `--describe` supplies
human- and machine-readable metadata; and a run records the original and
resolved parameter values.

For GUI use, a JSON description generated by the runner is preferable to
parsing arbitrary help text. A catalog entry could include:

```json
{
  "id": "boundary_control.partial_dirichlet",
  "label": "Partial Dirichlet control",
  "status": "available",
  "parameter_schema": "...",
  "requirements": ["scalar diffusion-reaction", "disjoint boundary roles"],
  "artifacts": ["state", "adjoint", "control", "optimization history"]
}
```

The final catalog schema and whether the runner or the application binary
owns registration should be decided before GUI implementation.

## User-facing run output

The terminal output should explain the work being done without forcing a user
to inspect source code. A typical run should identify:

```text
[run] id=2026-08-09T12:34:56Z-7f31
[input] target=boundary_control.partial_dirichlet
[input] prm=partial_dirichlet_tracking.prm
[compile] validating semantic problem
[compile] assembling mesh and DoF data
[solve] iteration=0 objective=... gradient_norm=...
[solve] iteration=1 objective=... gradient_norm=... accepted=true
[artifact] state=fields/state/step-0001.vtu
[artifact] control=fields/control/step-0001.vtu
[run] status=success iterations=... wall_time=...
```

Human output should remain concise and stable enough for a person or log
collector to follow. Detailed data belongs in structured files and events.
The machine-readable stream should use event records with an event type,
monotonic sequence number, run id, timestamp, and payload. Useful event types
include `run_started`, `input_resolved`, `compile_started`,
`compile_finished`, `iteration_started`, `trial_evaluated`,
`iteration_accepted`, `artifact_written`, `warning`, `error`, and
`run_finished`.

The GUI can show the human stream immediately and use events for progress
bars, current objective values, cancellation state, and artifact links. This
also allows Python scripts to consume a run without scraping terminal text.

## Run directory and artifact history

Every run should receive a unique directory and should never overwrite a
previous run by default. A proposed layout is:

```text
runs/
  boundary_control/
    partial_dirichlet/
      2026-08-09T123456Z-7f31/
        input/
          original.prm
          resolved.prm
          command.json
        metadata/
          run.json
          compilation_manifest.json
          environment.json
        logs/
          stdout.log
          stderr.log
          events.jsonl
        reports/
          summary.json
          optimization.json
          iterations.csv
          trials.jsonl
          linear_solves.csv
        fields/
          state/step-0000.vtu
          state/state.pvd
          adjoint/step-0000.vtu
          control/step-0000.vtu
          control/control.pvd
        artifacts.json
```

The exact layout is open, but the following properties are important:

- original input, resolved input, target id, source revision, and run id are
  retained;
- the compilation manifest is paired with the solver report and environment
  metadata;
- state, adjoint, and control histories have an explicit time/iteration
  index, even for a stationary PDE;
- boundary controls are written in a format that identifies the boundary
  support rather than pretending they are ordinary volume fields;
- `artifacts.json` maps logical names to files, formats, iteration numbers,
  and optional visualization hints; and
- a failed or cancelled run retains diagnostics and partial artifacts.

By default, write the initial state/control, every accepted optimization
iterate, and the final artifacts. Rejected line-search trials should normally
be kept in compact tabular form, not written as full field files unless a
debug output policy requests them. Output frequency, field selection, and
retention should be parameterized.

deal.II `DataOut`-style VTU/PVTU/PVD files are a natural first visualization
format for volume state and adjoint fields. Boundary controls need a separate
face or boundary-cell representation with a documented mapping. CSV is useful
for simple histories; JSON/JSONL is better for nested solver events and
artifact manifests. Python can consume the tabular/JSON data, while ParaView
can consume the visualization collection files.

## Solver reporting changes needed

The current reduced-gradient result contains useful aggregate histories and
final objects, but a GUI/artifact workflow needs more than the final state and
control. The solver layer should gain a typed report or observer boundary
that can expose, as available:

- objective components and total objective per iteration;
- gradient norms and stopping metrics;
- accepted step length, line-search trials, and decrease tests;
- state, adjoint, and linear-solve counts;
- iteration status and stopping reason;
- references or callbacks for accepted state, adjoint, and control snapshots;
- warnings and recoverable numerical failures; and
- cancellation/progress signals.

The generic solver should not know the final filesystem layout. An outer run
controller can subscribe to the report, serialize it, and ask the
discretization-specific application code to write fields. This follows the
existing separation between a compiled problem and execution policy while
leaving room for a future typed optimization report and run envelope.

## GUI sketch

The first Slint UI can be deliberately small:

1. **Select target/scenario** — show available, experimental, and unavailable
   entries with requirements and a short description.
2. **Edit parameters** — load a `.prm`, expose schema-backed fields, show
   defaults and validation errors, and allow saving a derived run input.
3. **Review and launch** — show the resolved target, parameter file, output
   policy, and estimated resource settings before starting.
4. **Run view** — stream stdout/stderr, show structured progress, objective and
   gradient plots, warnings, and a cancel action.
5. **Artifact view** — browse the manifest, open CSV/JSON summaries, render
   basic plots with a Python helper or embedded plotting component, and offer
   “open in ParaView” for VTU/PVD artifacts.

An embedded full mesh viewer should not be a first milestone. Launching
ParaView from an artifact manifest is less coupled and gives useful results
earlier. A native viewer can be added later if the run/artifact contracts are
already stable.

The GUI build should be optional, for example behind a CMake option, and its
dependency should not affect core-library or headless CI builds. The
cross-language boundary should be a process protocol or a small stable
library API, not direct access from Slint widgets into deal.II objects.

## Staged implementation proposal

This layer should be developed after the headless contracts are clear:

1. Define a run id, resolved-input record, exit-code policy, artifact manifest,
   and minimal structured summary.
2. Add a headless runner around one boundary-control scenario and make it
   produce immutable run directories.
3. Add typed solver events/reporting and accepted-iterate field snapshots.
4. Add catalog and parameter-schema discovery plus `--list`, `--describe`,
   and `--validate` operations.
5. Build the first Slint client as a runner subprocess client with target and
   parameter selection, live output, and artifact browsing.
6. Add Python plotting helpers and ParaView handoff; defer an embedded viewer
   and remote/batch execution until the local workflow is reliable.

The first boundary-control vertical slice is a substantial cross-cutting
feature rather than a small front-end task. A rough planning estimate is
approximately 12–20 review-sized implementation units for the runner,
run-envelope/reporting contracts, serialization, one scenario, and a useful
initial Slint UI. An embedded viewer, packaging, resume/cancel semantics,
multiple target families, and polished schema-driven editing would add
several more units. These are deliberately rough estimates until the
serialization and language-boundary choices are fixed.

## Decisions to settle before implementation

- Should the canonical run metadata and artifact manifest be JSON, TOML, or a
  versioned combination?
- Is `ParameterHandler` the long-term input format, or should it be wrapped by
  a generated schema/JSON representation?
- Which fields are retained by default, and how are large histories pruned?
- How are cancellation, timeout, resume, and partial failed runs represented?
- Which source revision and dependency/environment details are required for
  reproducibility?
- Should target registration live in the runner, separate scenario libraries,
  or a generated catalog?
- Which platforms must the Slint application package support?
- Is ParaView an optional external dependency, and what is the fallback when
  it is unavailable?

## Related authoritative material

The design should eventually be reconciled with the existing compiled-problem
manifest, solver result/reporting interfaces, and refactor assessment:

- [`semantic-compiler.md`](../implementation/v1/semantic-compiler.md)
- [`reduced_gradient.hpp`](../../include/nmopt/solvers/reduced_gradient.hpp)
- [`compiled_problem.hpp`](../../include/nmopt/compiler/v1/compiled_problem.hpp)
- [`refactor assessment`](../planning/refactor/assessment.md)

Until those interfaces are intentionally extended, this draft should remain a
discussion artifact and should not be treated as an implementation commitment.
