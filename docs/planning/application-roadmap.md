# Application roadmap and handoff

## Purpose and authority

This roadmap owns the user-facing application layer after a semantic problem
can be compiled and executed. It records the application work sequence, current
runner and artifact status, benchmark-execution handoffs, and application-level
acceptance gates.

It does not own compiler, lowerer, solver, or backend capability status. Those
remain in the [implementation roadmap](implementation-roadmap.md). Frozen
benchmark definitions remain in the [Chapter 6 benchmark specifications](../benchmarks/chapter-6.md),
scenario assembly remains in the [Chapter 6 application contract](../applications/chapter-6.md),
and the current execution contract remains in the [application execution and
artifact reference](../reference/application-execution.md).

The application roadmap is therefore the mutable status owner for:

- recipe, scenario, catalog, and application discovery surfaces;
- run configuration, reproduction policy, and generated-output implementation
  status;
- benchmark harnesses, artifact schemas, diagnostics, and provenance status;
- native deal.II mesh and field export status;
- Python and external-tool post-processing status;
- B0, B1, and B2 execution and acceptance handoffs; and
- future parameter-file and GUI boundaries.

## Execution reference

The current build loop, runner commands, schemas, run-set layout, native output,
report, and post-processing policies are maintained in the [application
execution and artifact reference](../reference/application-execution.md). This
roadmap records their implementation status and handoffs; it does not duplicate
the operational contract.

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
[benchmark specification](../benchmarks/chapter-6.md), not duplicated
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

The runner accepts an explicit `--run-kind reproduction|development` policy.
The selected benchmark supplies its mesh default; an explicit `--refinement`
is only an optional override. The current policy and generated-record
requirements are defined in the [execution reference](../reference/application-execution.md).

### A2 — Separate runner configuration and path layout

**Status:** runner configuration and deterministic path layout implemented.

The runner resolves CLI options into a run configuration. Its current path,
slot, metadata, and derived-output policy is maintained in the [execution
reference](../reference/application-execution.md). Frozen numerical values
remain in typed scenario factories and benchmark contracts.

### A3 — Define native deal.II output

**Status:** topology-specific native field names, standalone volume mesh
exports, a 2D SVG mesh preview, and B2 boundary-control topology export
implemented. A separate boundary-only mesh without a field remains optional.

The execution reference records the current native filenames and topology
rules. This roadmap retains only the implementation status and optional
boundary-only export follow-up.

### A4 — Define records, diagnostics, and provenance

**Status:** run-set manifest, artifact inventory, failure records, and
manifest-aware report selection implemented.

The runner now persists a run manifest, artifact inventory, failure records, and
manifest-aware report selection. The [execution
reference](../reference/application-execution.md) defines the current record
contents and tool behavior.

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

The Python post-processing backend and report integration are implemented. The
supported inputs, outputs, profiles, and compatibility wrapper are documented
in the [execution reference](../reference/application-execution.md). Plotting
remains optional and is not a C++ build or CTest dependency.

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
