# Application roadmap and handoff

## Purpose and authority

This roadmap owns the user-facing application layer after a semantic problem
can be compiled and executed. It records the application work sequence, current
runner and artifact status, benchmark-execution handoffs, and application-level
acceptance gates.

It does not own compiler, lowerer, solver, or backend capability status. Those
remain in the [implementation roadmap](implementation-roadmap.md). Frozen
benchmark definitions remain in the [Chapter 6 benchmark specifications](../benchmarks/chapter-6/README.md),
and scenario assembly remains in the [Chapter 6 application contract](../applications/chapter-6/README.md).

The application roadmap is therefore the mutable status owner for:

- recipe, scenario, catalog, and application discovery surfaces;
- run configuration, reproduction policy, and generated-output organization;
- benchmark harnesses, artifact schemas, diagnostics, and provenance;
- native deal.II mesh and field export;
- Python and external-tool post-processing;
- B0, B1, and B2 execution and acceptance handoffs; and
- future parameter-file and GUI boundaries.

## Quick start and build policy

Use the backend-neutral Debug profile for ordinary contract work:

```bash
cmake --preset debug-neutral
cmake --build --preset debug-neutral
ctest --preset debug-neutral --output-on-failure
```

For application development that needs the deal.II runner, use the Debug
deal.II profile and build only the runner target:

```bash
cmake --preset debug-dealii
cmake --build --preset debug-dealii --target nmopt_runner --parallel 1
build/debug-dealii/bin/nmopt_runner --list
```

Development runs must opt into the development policy and may use a smaller
mesh, for example:

```bash
build/debug-dealii/bin/nmopt_runner \
  --benchmark b1 \
  --framework-revision REV \
  --run-kind development \
  --refinement 1 \
  --output runs
```

The runner will place this below
`runs/chapter-6/b1/development/001/`.

Do not configure or compile `release-dealii` unless it is absolutely
necessary for source-scale reproduction or optimized verification, and ask
the user for explicit permission before doing so. A request to implement,
document, inspect, or report on an application does not grant that
permission. Debug deal.II is the default for development and diagnostics;
Debug output must not be presented as reproduction evidence.

Use `python3` explicitly for reports and post-processing. Check the
interpreter and optional plotting libraries with:

```bash
python3 --version
python3 -c "import matplotlib, meshio; print('matplotlib and meshio available')"
```

The current deterministic report tool uses only the standard library:

```bash
python3 tools/chapter6_report.py \
  --run-manifest runs/chapter-6/b1/authoritative/run-manifest.json \
  --output runs/chapter-6/b1/authoritative/report
```

The report should select one authoritative run set through its
`run-manifest.json`; `--input` remains available for legacy or aggregate
inspection. The A7 post-processing tool uses `meshio` and `matplotlib` to read native
deal.II VTU files and produce derived state, control, and adjoint plots. Both
single-artifact and run-root modes are available; aggregate report integration
is now manifest-aware and must not become a prerequisite for the runner.

## Application boundaries

The application layer has five distinct responsibilities:

| Boundary | Responsibility | Does not own |
| --- | --- | --- |
| Recipe/scenario API | Assemble typed application choices and public `ProblemSpec` values. | Meshes, compilers, solver loops, or files. |
| Runner configuration | Resolve the requested benchmark, run kind, build profile, revision, and overrides. | PDE semantics or unsupported capabilities. |
| Run controller | Select paths, create run sets, dispatch typed execution adapters, and retain failures. | A second compiler or optimization algorithm. |
| Native output | Ask deal.II to write meshes and finite-element fields with their actual topology. | Plots, summary statistics, or reconstructed fields. |
| Post-processing | Read native outputs and records with Python or external tools and produce derived views. | Changes to the authoritative run data. |

The core application API remains under `include/nmopt/application/`. The
headless executable remains under `apps/nmopt-runner/`. Generated output stays
under the ignored `runs/` directory.

## Status vocabulary

Application work uses these meanings:

- **planned** — scoped but not implemented;
- **implemented** — source and focused tests exist;
- **development-verified** — exercised on a development mesh or test run;
- **reproduction-verified** — exercised with the benchmark's declared mesh
  policy and `release-dealii` profile;
- **acceptance-complete** — all benchmark evidence and gates are present; and
- **deferred** — intentionally outside the current application sequence.

Development output must not be presented as reproduction evidence. An
investigation may explain a numerical choice, but it does not replace the
frozen benchmark run.

## Current handoff state

The typed application API, Chapter 5 recipes, Chapter 6 B1/B2 scenarios, B0
harness, deterministic artifact writer, headless runner boundary, and selected
deal.II execution adapters are implemented. The current code also has B1
Hessian evidence, Armijo trial traces, and native VTU field export for the
selected applications.

The current acceptance state is:

### B0

B0 execution plumbing, manifest selection, artifact inventory, failure records,
and deterministic report scope are implemented. The benchmark-specific
acceptance evidence is tracked below.

### B1

The authoritative release run set at
`runs/chapter-6/b1/authoritative/` is complete at the benchmark-default
refinement 7. All six artifacts contain solver traces, native mesh/field
exports, provenance, and passed centered finite-difference Hessian evidence.
B1 is acceptance-complete for the current framework-native contract. The
lower-regularisation cases still record iteration-limit or line-search
limitations, which are part of the reported result.

### B2

The authoritative release run set at
`runs/chapter-6/b2/authoritative/` is complete at the benchmark-default
refinement 7. All four cases record the ordinary-normal transport convention,
solver traces, native volume fields, native boundary-control fields, and valid
compilation diagnostics. B2 reproduction is therefore verified, but acceptance
remains open because the artifacts do not yet persist the required residual
JVP/VJP, reduced-Taylor, and explicit state/control comparison evidence.

The frozen B2 Galerkin formulation remains in force. Stabilization is not an
implicit application-layer fix.

Authoritative frozen inputs and required evidence are maintained in the
[benchmark specification](../benchmarks/chapter-6/README.md), not duplicated
here.

## Roadmap sequence

### A0 — Establish application authority

**Status:** implemented by this document.

Create the separate application roadmap, route it from `docs/README.md`, and
remove application acceptance status from the compiler-focused roadmap except
where it is a dependency or cross-reference.

### A1 — Define run sets and reproduction policy

**Status:** runner policy and run-set persistence implemented; adaptive or
mesh-structure-based resolution remains a future extension.

The runner now accepts an explicit `--run-kind reproduction|development`
policy. Reproduction is the default and requires the `release-dealii` build
profile, but it does not require a runner-level refinement number. The
selected benchmark supplies its mesh default; an explicit `--refinement` is
only an optional override. Future mesh resolvers may estimate resolution from
actual vertices, faces, cells, or structural identity, and adaptive meshes
must be identified by their strategy and realized structure rather than by a
single integer. The policy is validated before dispatch.

The generated layout for reproduction, development, investigation, and
derived reports is defined by A2. Each run set now records its benchmark,
framework revision, profile, refinement override, expected matrix size,
command, completion state, and per-artifact status in `run-manifest.json`.

Reproduction means `release-dealii` plus the benchmark's declared mesh
policy. Development runs may use smaller meshes or other explicitly selected
mesh policies. The runner should record the actual compiled profile and
resolved mesh structure, and reject or invalidate a reproduction request made
with a different profile.

### A2 — Separate runner configuration and path layout

**Status:** runner configuration and deterministic path layout implemented.

The runner now resolves CLI options into a run configuration and places each
benchmark below:

```text
<output>/chapter-6/<benchmark>/<run-slot>/
```

For example:

```text
runs/chapter-6/b1/authoritative/
runs/chapter-6/b2/development/001/
```

The `authoritative` slot is the single reproduction run set. Development runs
receive progressive slots such as `001` and `002`. The framework source
commit or tag, actual build profile, command, mesh resolution, and full
scenario IDs remain in artifact/run metadata rather than in the path. The
matrix dimensions (method, regularization, or case) remain as concise
directories below the run-set's `artifacts/` directory. Run-set reports and
aggregate post-processing remain siblings of `artifacts/`.
Keep frozen numerical values in typed scenario factories; do not add a general
parameter-file system as a prerequisite for B1/B2.

### A3 — Define native deal.II output

**Status:** topology-specific native field names, standalone volume mesh
exports, a 2D SVG mesh preview, and B2 boundary-control topology export
implemented. A separate boundary-only mesh without a field remains optional.

Deal.II writes the authoritative native field datasets directly:

```text
native/
  mesh-volume.vtu
  mesh-volume.svg           # 2D quick preview of the volume mesh
  mesh-boundary.vtu          # when boundary-region evidence is needed
  fields-volume.vtu          # state/adjoint, plus B1 volume control
  control-boundary.vtu       # B2 facewise control
```

The volume and boundary files intentionally remain separate because they have
different discrete topologies. A run manifest will associate logical fields
with files and supports such as `volume` or `boundary_faces`.

ParaView, meshio, matplotlib, and other consumers read these native files;
they must not reconstruct the mesh or finite-element fields.

### A4 — Define records, diagnostics, and provenance

**Status:** run-set manifest, artifact inventory, failure records, and
manifest-aware report selection implemented.

The runner writes `run-manifest.json` before execution and updates it after
each artifact. It records the expected matrix, command, benchmark, run kind,
build profile, framework revision, refinement override, completion state, and
per-artifact success or failure. `artifact.kv` remains the deterministic
per-artifact projection containing the accepted-iteration history, measurement
flags, environment, compilation manifest, and benchmark fields. The Chapter 6
report accepts `--run-manifest`, automatically consumes a manifest when `--input`
points directly at a run set, and retains missing or failed inventory records
with their diagnostics in `summary.csv` and `summary.md`.

### A5 — Refresh B1 reproduction evidence

**Status:** reproduction-verified; acceptance-complete for the current
framework-native contract.

The six-case B1 matrix was run with `release-dealii` without a runner-level
refinement override, so the benchmark supplied refinement 7. The authoritative
run set retains native fields, the `GridOut` mesh SVG, solver traces, passed
centered finite-difference Hessian evidence, and a manifest-driven report. The
lower-$\beta$ iteration-limit and line-search outcomes remain documented
execution-policy limitations.

### A6 — Refresh B2 reproduction evidence

**Status:** reproduction-verified; artifact evidence implemented; report and
release refresh pending.

The four-case source-scale matrix was refreshed with the authorized
ordinary-normal convention and the benchmark-default refinement 7. The B2
adapter now computes and persists residual JVP/VJP, reduced-Taylor,
objective/gradient reduction, and case-label/state/control diagnostics; the
focused Debug contract verifies the fields. The next small unit must project
these values into the report and compare the four cases, then request
permission before rebuilding `release-dealii` and refreshing the authoritative
matrix. Record any remaining Galerkin limitation without silently activating
stabilization.

### A7 — Add Python post-processing

**Status:** profile-driven single-artifact renderer, run-root processing, and
manifest-aware report integration implemented; combined visualization/report
orchestration remains planned.

The first supported post-processing backend is Python with `meshio` and
`matplotlib`. It should read only native VTU files and runner records and
produce derived plots for state, adjoint, volume control, and boundary
control. `tools/postprocess.py` is the generic entry point and selects an
application profile with `--profile` and a Matplotlib format with `--format`;
PNG is the default, SVG is available when vector output is preferred, and
`--format png svg` requests both. The current `chapter6` profile reads one
artifact, writes plots and
`postprocess.json` below an artifact-local `postprocess/` directory, and
supports legacy field filenames for historical development runs. `--input`
keeps those per-artifact outputs colocated and writes only the aggregate index
and comparisons below a run-set-level `postprocess/` output root. It writes
`postprocess-index.json` with per-artifact success or failure records.
Run-root mode also creates shared-scale comparison figures grouped by
benchmark family, so B1 and B2 fields are never compared in the same panel.
`tools/chapter6_postprocess.py` remains a compatibility wrapper for the
Chapter 6 profile. Plotting remains optional and must not be a C++ build or
CTest dependency.

### A8 — Add discovery and parameter-file interfaces

**Status:** deferred until A1–A7 stabilize.

Add `--describe`, machine-readable catalog output, and typed parameter-file
loading only after the run and artifact contracts are stable. A parameter file
may select values within a registered scenario; it may not activate an
unsupported semantic or compiler capability.

### A9 — GUI/client layer

**Status:** deferred.

The GUI should consume the same catalog, runner, run manifest, native fields,
and derived reports as scripts. It should not introduce a second execution
path.

## Handoff protocol

Each application unit is one review-sized work unit unless its contract change
requires a smaller split. A handoff must report:

```text
Unit:
Purpose:
Files changed:
Contract changed:
Focused verification:
Generated artifacts:
Known limitations:
Next unit:
```

The next remaining application unit is the B2 report-evidence subunit: project
the new artifact diagnostics into the report and add the four-case comparison.
It must be followed by an explicit permission request before rebuilding
`release-dealii` and rerunning the authoritative matrix.

## Exclusions

The current application roadmap does not activate B3/B4, all-at-once B5/B6,
preconditioning, stabilization comparisons, Stokes applications, automatic
OtD derivation, continuous-control boxes, remote execution, or a GUI. Those
remain separate decisions under the benchmark and implementation roadmaps.
