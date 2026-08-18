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
- **reproduction-verified** — exercised with the frozen source-scale policy and
  `release-dealii` profile;
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

The acceptance state is still open:

### B0

B0 execution plumbing is implemented. The remaining application work is to
make the run-set identity, artifact inventory, failure records, and report
scope explicit rather than relying on recursive discovery below `runs/`.

### B1

The six source-scale artifacts under `runs/` predate the later solver-trace,
native-field, and finite-difference Hessian additions. Refinement-1
development artifacts contain the newer evidence. B1 requires a refreshed
source-scale matrix before it can be marked acceptance-complete.

### B2

The current semantic and deal.II paths record the ordinary-normal transport
boundary convention. The canonical source-scale artifacts predate that
correction. Development evidence includes both successful ordinary-normal
investigations and lower-refinement line-search limitations; the source-scale
four-case matrix must be refreshed and supplemented with the required residual,
Taylor, and state/control comparison evidence.

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

**Status:** planned.

Define the generated layout for reproduction, development, investigation, and
derived reports. A run set must record its benchmark, framework revision,
profile, refinement, expected matrix size, command, and completion state.

Reproduction means `release-dealii` and the frozen source-scale refinement.
Development runs may use smaller meshes. The runner should record the actual
compiled profile and reject or invalidate a reproduction request made with a
different profile.

### A2 — Separate runner configuration and path layout

**Status:** planned.

Split the current executable responsibilities into CLI parsing, resolved run
configuration, run-layout/path selection, and B1/B2 command dispatch. Keep
frozen numerical values in typed scenario factories; do not add a general
parameter-file system as a prerequisite for B1/B2.

### A3 — Define native deal.II output

**Status:** partially implemented.

Deal.II should write the authoritative native datasets directly:

```text
native/
  mesh-volume.vtu
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

**Status:** partially implemented.

Add a per-run manifest, explicit artifact inventory, accepted-iteration solver
history, accurate measurement flags, complete environment provenance, and
failure records. Preserve `artifact.kv` as a deterministic projection if the
machine-readable manifest evolves separately.

### A5 — Refresh B1 reproduction evidence

**Status:** planned.

Run the current six-case B1 matrix with `release-dealii`, retain native fields
and solver traces, verify the Hessian evidence, and generate a report from that
run set only. Mark B1 acceptance-complete only when the frozen evidence list is
present.

### A6 — Refresh B2 reproduction evidence

**Status:** planned.

Add residual JVP/VJP, reduced-Taylor, objective/gradient reduction, and
four-case state/control comparison evidence. Refresh all four source-scale
runs under the authorized ordinary-normal convention. Record any remaining
Galerkin limitation without silently activating stabilization.

### A7 — Add Python post-processing

**Status:** planned.

The first supported post-processing backend is Python with `meshio` and
`matplotlib`. It should read only native VTU files and runner records and
produce derived PNG/SVG plots for state, adjoint, volume control, and boundary
control. Plotting remains optional and must not be a C++ build or CTest
dependency.

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

The next unit after A0 is A1: define the run-set layout and reproduction
policy. It must stop before changing native field names, artifact schemas, or
Python tooling unless those changes are independently split and reviewed.

## Exclusions

The current application roadmap does not activate B3/B4, all-at-once B5/B6,
preconditioning, stabilization comparisons, Stokes applications, automatic
OtD derivation, continuous-control boxes, remote execution, or a GUI. Those
remain separate decisions under the benchmark and implementation roadmaps.
