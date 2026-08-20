# Application roadmap and handoff

## Purpose and authority

This roadmap owns the mutable application-layer work after a semantic problem
can be compiled and executed. It records what the application units are for,
what has actually been implemented, how a unit becomes complete, and what the
next handoff is.

It does not own compiler, lowerer, solver, or backend capability status. Those
remain in the [implementation roadmap](implementation-roadmap.md). Frozen
benchmark definitions remain in the [Chapter 6 benchmark contracts](../benchmarks/chapter-6.md),
scenario assembly remains in the [Chapter 6 application contract](../applications/chapter-6.md),
and the operational schemas and commands remain in the
[application execution and artifact reference](../reference/application-execution.md).

The application roadmap owns:

- recipe, scenario, catalog, and application discovery surfaces;
- runner configuration, reproduction policy, and generated-run organization;
- benchmark harnesses, artifact schemas, diagnostics, and provenance plumbing;
- native deal.II mesh and field export;
- Python post-processing and report tooling;
- B0–B2 execution and reproduction handoffs; and
- future parameter-file and client boundaries.

## How to use this roadmap

Each unit has four questions:

1. **Purpose:** why the unit exists and which application boundary it owns;
2. **Work:** the implementation and documentation tasks included in it;
3. **Done when:** the evidence required before the unit can be closed; and
4. **Actual status:** what is already present and what remains open.

A unit is not complete merely because its source compiles. It is complete when
its focused contract tests, generated artifacts, and documentation support the
stated use. A development run can verify orchestration and artifact shape, but
cannot close a source-sized reproduction gate.

The roadmap distinguishes these states:

- **planned** — scoped but not started;
- **implemented** — source and focused tests exist;
- **development-verified** — exercised on a development mesh or run;
- **framework-verified** — exercised at the declared framework-native mesh and
  profile, without claiming source-result parity;
- **reproduction-verified** — source-sized output has been compared against the
  source record under the frozen project contract;
- **acceptance-complete** — all benchmark evidence and gates are present; and
- **deferred** — intentionally outside the current sequence.

## Current status

| Area | Status | Actual state |
| --- | --- | --- |
| Typed application API, recipes, scenarios, and catalog | Implemented for the selected B1/B2 slice | Chapter 5 scalar and Neumann/convection recipes, typed Chapter 6 scenarios, metadata, and catalog entries are available. The complete Chapter 5 problem library is still future work. |
| Runner and run sets | Implemented | `nmopt_runner` resolves benchmark, run kind, build profile, revision, and optional refinement; it writes manifests, artifacts, failures, and native outputs. |
| B0 artifact boundary | Implemented and contract-tested | The deterministic `nmopt-benchmark-v1` projection, manifest records, provenance fields, selected fields, and failure visibility are present. |
| Native deal.II output | Implemented for B1/B2 | B1 writes volume mesh/fields; B2 additionally writes facewise control on the boundary topology, together with solver traces and mesh previews. |
| Reports and post-processing | Implemented as tooling | Python reads persisted `artifact.kv`, solver traces, and VTU files; it produces field plots, comparisons, post-process indexes, and deterministic benchmark reports. The renderer now uses `turbo` and explicit field-extrema colorbar endpoints; scientific visual parity is not yet signed off. |
| B1 framework-native execution | Framework-verified; active matrix expanded | The existing six-artifact `release-dealii`, refinement-7 matrix and Hessian evidence remain available; the active B1 matrix now includes the source's $\beta=10^{-6}$ field case, while the authoritative release refresh remains gated by the reproduction audit. |
| B2 framework-native execution | Framework-verified; reproduction audit open | The four-artifact release matrix exists, while the newer derivative-evidence fields are present in a refinement-6 development run rather than the authoritative release set. The source-sized release artifacts need a refresh after the audit. |
| Parameter files | Planned | No stable Deal.II-style `.prm` boundary exists yet. |
| Later Chapter 6 benchmarks | Planned | B3/B4 are the next selected benchmark families after B1/B2 reproduction is resolved; B5/B6 remain later. |
| Complete Chapter 5 recipe library | Planned | The selected recipes used by B1/B2 exist, but the reusable recipe families listed in the problem-library roadmap are not all implemented. |

The current authoritative run manifests record complete six-case B1 and
four-case B2 release matrices at benchmark-default refinement 7. These runs
demonstrate that the framework-native applications execute and produce
structured evidence. They do not by themselves establish that the rendered
fields or numerical conventions reproduce the book's figures.

The B1 report records two gradient-tolerance terminations, one line-search
failure, and three maximum-iteration terminations across the two methods and
three regularisation values. The B2 authoritative report records complete
artifacts and native fields, but all four cases terminate at the configured
maximum iteration count. These outcomes are evidence to preserve, not reasons
to silently change the frozen benchmark policy.

The development snapshot
`runs/chapter-6/b2/development/001/` was produced with refinement 6 and
contains the newer `b2.*` derivative, objective-reduction, and
state/control evidence. It is useful for diagnosis but is not source-sized
acceptance evidence. The existing `postprocessed/<case>/state.png` files
belong to that development snapshot; the current wrapper writes the
`postprocess/` directory name.

## Application boundaries

The application layer has five distinct responsibilities:

| Boundary | Responsibility | Does not own |
| --- | --- | --- |
| Recipe/scenario API | Assemble typed application choices and public `ProblemSpec` values. | Meshes, compilers, solver loops, or files. |
| Runner configuration | Resolve benchmark, run kind, build profile, framework revision, and overrides. | PDE semantics or unsupported capabilities. |
| Run controller | Select paths, create run sets, dispatch typed execution adapters, and retain failures. | A second compiler or optimization algorithm. |
| Native output | Ask deal.II to write meshes and finite-element fields with their actual topology. | Plots, summaries, or reconstructed fields. |
| Post-processing | Read native outputs and records with Python or external tools and produce derived views. | Changes to authoritative run data. |

The core application API remains under `include/nmopt/application/`. The
headless executable remains under `apps/nmopt-runner/`. Generated output stays
under the ignored `runs/` directory.

## Post-processing and reporting

The application output pipeline is:

```text
nmopt_runner
  -> artifact.kv, solver-trace.csv, native/*.vtu, native/*.svg
  -> tools/postprocess.py
  -> artifact plots, postprocess.json, comparison plots, postprocess-index.json
  -> tools/chapter6_report.py
  -> summary.csv, summary.md
```

The current tooling has these responsibilities:

- `tools/run_chapter6.sh` runs an already-built B1/B2 executable, allocates
  the run slot, invokes post-processing, and generates the manifest-aware
  report. It is a convenience wrapper, not a build system.
- `tools/postprocess.py` is the profile-driven entry point. The `chapter6`
  profile reads `fields-volume.vtu` and `control-boundary.vtu`, extracts
  state, target, forcing, native and book-sign adjoints, and control fields,
  and supports PNG output for the current reproduction workflow.
- `tools/nmopt_postprocess/` contains the reusable mesh, field, rendering,
  comparison, and profile code. Run-root comparisons are grouped by scenario
  family so B1 and B2 are not combined accidentally.
- `tools/chapter6_report.py` reads persisted artifacts and sidecars without
  rerunning a solver or inferring missing values. It retains failed, pending,
  and missing manifest entries in the report.

The current renderer uses Matplotlib's `turbo` colormap, point-field Gouraud
interpolation, cell-field flat shading, no volume mesh-edge overlay, and a
shared finite-value normalization for comparison panels. These are current
implementation defaults, not claims about the book's plotting configuration.
Their scientific adequacy is part of the next reproduction unit.

## Work units

### A0 — Establish application authority

**Status:** acceptance-complete.

**Purpose:** Give application-layer status, execution handoffs, and output
policy one mutable owner separate from compiler and implementation capability
records.

**Work completed:**

- Added this roadmap and routed it from `docs/README.md`.
- Separated benchmark source facts, application construction, execution
  interfaces, and mutable status.
- Added the application execution reference for schemas, commands, run layout,
  native output, reporting, and post-processing.

**Done when:** an agent can identify the authority for scenario construction,
benchmark choices, execution commands, implementation capability, and current
application status without consulting duplicated status ledgers. This is met.

### A1 — Define typed application construction surfaces

**Status:** implemented for the selected B1/B2 application slice.

**Purpose:** Provide reusable typed recipes, scenario records, metadata, and
catalog discovery without creating a new problem class for each benchmark
combination.

**Work completed:**

- Chapter 5 scalar distributed-control and Neumann/convection recipes exist.
- Chapter 6 B1 and B2 scenarios freeze recipe IDs, runtime ports, provenance,
  mesh defaults, solver policy, and benchmark identity.
- The catalog exposes the selected Chapter 5/6 entries.
- Backend-neutral scenario validation rejects incomplete provenance and
  unsupported product or execution choices.

**Still needed:** complete the standard Chapter 5 recipe families listed in the
[problem-library roadmap](chapter-5-problem-library-roadmap.md). That work is
tracked as A9 below; it must reuse these boundaries rather than introduce
benchmark-specific application classes.

**Done when:** each selected application can produce a validated
backend-neutral `ProblemSpec`, expose its required runtime bindings, and be
consumed by the benchmark harness. This is met for B1/B2.

### A2 — Define runner configuration and run sets

**Status:** implemented and development-verified.

**Purpose:** Make benchmark execution reproducible at the application boundary:
resolve policy, choose a run slot, retain the command and revision, and keep
failed matrix entries visible.

**Work completed:**

- `nmopt_runner` accepts `--benchmark`, `--run-kind`,
  `--framework-revision`, `--output`, and optional `--refinement`.
- Reproduction requires the compiled `release-dealii` profile; development
  runs may use an explicit smaller mesh.
- Reproduction uses the `authoritative` slot; development uses the next
  numbered slot.
- The runner writes a run manifest before the matrix and records each
  successful or failed artifact.

**Done when:** repeated invocations with the same frozen scenario and build
revision select deterministic artifact paths and metadata, while an exception
cannot make a matrix entry disappear. The focused runner and artifact tests
cover this boundary.

### A3 — Define native deal.II output

**Status:** implemented for the selected applications.

**Purpose:** Preserve authoritative numerical topology and fields before any
Python rendering or format conversion.

**Work completed:**

- B1 writes `native/mesh-volume.vtu`, `native/mesh-volume.svg`, and
  `native/fields-volume.vtu`, including state, target, forcing, native
  adjoint, book-sign adjoint, and control fields.
- B2 writes the same volume outputs plus
  `native/control-boundary.vtu` for its facewise control.
- Solver Armijo trials are retained in `solver-trace.csv`.
- The output writer keeps final-state fields separate from per-iteration
  solver history.

**Still optional:** a separate B2 boundary-only mesh export without a control
field. It is not required by the current benchmark contract.

**Done when:** native files can be read independently of the plotting tools,
their topology matches the realized compilation manifest, and reports never
need to reconstruct missing authoritative values. This is met for the current
B1/B2 output inventory.

### A4 — Define records, diagnostics, and provenance

**Status:** implemented; B2 source-sized refresh remains open.

**Purpose:** Make every run auditable: identify what was selected, how it was
compiled and executed, what data was bound, and why a matrix entry failed.

**Work completed:**

- The deterministic `nmopt-benchmark-v1` artifact projection is implemented.
- Artifacts retain identity, provenance, compilation manifest, diagnostics,
  solver histories, benchmark fields, measurements, and selected fields.
- B1 records forcing and Hessian finite-difference evidence.
- B2 development artifacts record boundary form, case, provenance, derivative
  checks, objective/gradient reduction, and state/control norms.
- The report consumes the run manifest rather than discovering an arbitrary
  directory tree.

**Still needed:** refresh the authoritative B2 matrix so it contains the
current B2 evidence fields, then project those fields into the report's
comparison and acceptance views.

**Done when:** the source-sized artifacts for an activated benchmark contain
all evidence listed in its benchmark contract, and the report preserves
successes, failures, and missing evidence without inventing values. B0 and B1
meet this boundary; B2 does not yet.

### A5 — Build and validate post-processing tooling

**Status:** implementation-complete; scientific visual validation open.

**Purpose:** Turn native fields into inspectable plots and comparisons without
making post-processing a second numerical implementation.

**Work completed:**

- Single-artifact and run-root processing are available through
  `tools/postprocess.py`.
- The Chapter 6 profile handles volume state/adjoint/control and boundary
  control fields; B2 additionally exports and plots the zero-control state,
  target, forcing, and observation-domain mask.
- PNG and SVG output, artifact-local manifests, run-root indexes, and
  grouped comparisons are available.
- `chapter6_report.py` generates deterministic CSV/Markdown summaries, while
  the profile-driven postprocessor renders objective and gradient histories.
- The compatibility `chapter6_postprocess.py` wrapper remains available for
  older Chapter 6 inputs.
- `meshio` and `matplotlib` are runtime dependencies only; they are not
  C++ build or CTest dependencies.

**Still needed:** validate field orientation, value ranges, normalization,
interpolation, colorbar semantics, and source-figure palette before treating
plots as reproduction evidence. A visually plausible plot is not enough.

**Done when:** the tool preserves field identity and geometry, documents its
rendering policy, produces deterministic comparisons from persisted inputs,
and a benchmark-specific visual check confirms that any difference from the
source is understood. The plumbing is met; the visual criterion is not.

### A6 — Correctly reproduce B1 and B2

**Status:** in progress; development-only evidence is being audited.

**Purpose:** Establish whether the current framework-native fields and plots
represent the same mathematical quantities and visual conventions as the
book's B1/B2 figures.

**Work to do:**

1. Build a source-to-artifact comparison sheet for E6.5.1 and E6.5.2 using
   the guide's extracted pages and figures. Record the source field, coordinate
   orientation, sign convention, observation/target case, and whether the
   comparison is qualitative or quantitative.
2. Inspect raw native values before plotting. Compare VTU field names, ranges,
   extrema, mesh orientation, and point/cell association against the state,
   adjoint, and control quantities selected by the benchmark contract.
3. Render the B1 input functions as first-class comparison evidence: the
   desired state and forcing, with their provenance and the same coordinate
   orientation as the native fields. A zero forcing must still be visible as a
   field with an explicit zero range rather than inferred from the PDE name.
4. Reproduce the current PNG from the same VTU input with controlled changes to
   colormap, normalization, interpolation, axis orientation, and colorbar.
   Include the exact field extrema as colorbar endpoints.
5. Compare the source-like blue-cyan-green-yellow-red colormap with the
   selected renderer policy. The current reproduction policy is `turbo`; the
   source plotting configuration remains unavailable, so palette equivalence
   is visual rather than a recovered plotting setting.
6. Keep the native framework adjoint `p` unchanged, and add a comparison-only
   book-convention field `p_book=-p` when comparing against the source figures.
   Record this transformation in the plot metadata; never overwrite the raw
   adjoint field.
7. If the raw field is wrong, trace the issue through native export, boundary
   orientation, target/observation realization, and PDE/adjoint conventions.
   If the raw field is correct and only rendering differs, make the smallest
   post-processing change and add a focused regression check.
8. Regenerate development plots, including the B1 $\beta=10^{-6}$ case, then
   request the explicit release policy
   before refreshing the B1/B2 authoritative matrices. Update the report and
   benchmark handoff with the result.

**Done when:**

- the color inversion and palette difference have a documented cause;
- raw numerical fields and plotted fields have an auditable mapping;
- the framework adjoint and book-convention adjoint are explicitly separated;
- B1 and B2 plots use an explicitly documented normalization, endpoint ticks,
  and colormap policy;
- B1 input functions are plotted with their provenance and coordinate mapping;
- the B1 matrix includes the source-referenced $\beta=10^{-6}$ field case;
- the source-sized matrices are regenerated only after the interpretation is
  settled; and
- the benchmark report separates framework-native evidence, source comparison,
  and unresolved source omissions.

The first diagnostic input is:

```text
runs/chapter-6/b2/development/001/postprocessed/<case>/state.png
```

It is development evidence only. Do not change the frozen Galerkin,
ordinary-normal, case, or mesh policy merely to make a plot look like the
book.

### A7 — Add Deal.II-style parameter files

**Status:** planned; after A6.

**Purpose:** Allow an agent or user to reproduce a declared scenario from a
reviewable `.prm` file instead of relying on hidden defaults or a long CLI
command.

**Work to do:**

- Define which values belong to the scenario parameter file and which remain
  runner policy or source-controlled benchmark contract.
- Use a Deal.II-compatible `ParameterHandler` structure with explicit
  subsections for benchmark identity, mesh, runtime data, solver, output, and
  post-processing.
- Load and validate parameters into the typed scenario records; reject values
  outside the registered capability boundary.
- Record the parameter-file path, content hash or revision, and effective
  values in the run manifest and artifact.
- Preserve CLI overrides only where precedence is explicit and auditable.

**Done when:** a checked-in `.prm` reproduces a development run, the
effective configuration is visible in the artifact, invalid combinations fail
before execution, and a benchmark can declare whether its parameters are
frozen or overridden.

### A8 — Implement the next Chapter 6 benchmarks

**Status:** planned; after A6, with A7 optional for the first slice.

**Purpose:** Extend the application layer beyond the selected B1/B2 vertical
slice while reusing the same recipe, runner, artifact, native-output, and
post-processing boundaries.

**Work sequence:**

- Start with B3 / E6.9.1 symmetric box-constrained Laplace control.
- Follow with B4 / E6.9.2 asymmetric box-constrained Laplace control.
- Add each benchmark's source omissions, frozen scenario choices, run matrix,
  required KKT/PDAS evidence, and native fields before implementation.
- Reuse the Chapter 5 problem-library recipe and the existing benchmark
  harness; do not add a benchmark-specific optimizer or PDE lowerer.
- Defer B5/B6 until the all-at-once/KKT application boundary is explicitly
  activated by the benchmark roadmap.

**Done when:** each selected benchmark has a frozen contract, typed scenario,
runner matrix, derivative/KKT evidence, source-sized report, native fields,
and a documented acceptance decision. A compiled example without those
records is not a completed benchmark unit.

### A9 — Complete the Chapter 5 problem library

**Status:** planned; selected B1/B2 recipes are the initial partial slice.

**Purpose:** Turn the existing composition interfaces into a coherent,
parameterized library of reusable Chapter 5 application recipes.

**Work to do:**

- Complete the scalar distributed-control, general scalar/Robin,
  Neumann-boundary, observation, and registered Dirichlet-control families
  listed in the [Chapter 5 problem-library roadmap](chapter-5-problem-library-roadmap.md).
- Keep mathematical composition in the semantic/compiler services; recipes
  provide typed parameters, defaults, provenance, and validation.
- Give every family a fast manufactured configuration, focused derivative or
  KKT evidence, and at least one direct-composition comparison.
- Consume recipes from exploratory drivers and later benchmarks without
  duplicating assembly or solver logic.
- Update the Chapter 5 application contract and catalog only when a family is
  actually available.

**Done when:** each advertised family validates and compiles through the
selected path, supported parameter changes affect only declared components,
unsupported combinations fail with capability diagnostics, and at least one
configuration is consumed unchanged by the benchmark harness.

### A10 — Client and GUI boundary

**Status:** deferred.

A future client should consume the catalog, parameter files, runner, manifests,
native fields, and reports. It must not introduce a second execution path.
This unit remains deferred until the command-line and parameter-file contracts
are stable.

## Handoff protocol

Each application unit is one review-sized work unit unless its contract change
requires a smaller split. Every handoff must report:

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

The next handoff is A6, correct B1/B2 reproduction. It should begin with
development-only inspection and post-processing experiments. Rebuilding
`release-dealii` or replacing authoritative artifacts requires an explicit
permission request after the cause of the visual discrepancy is understood.

## Exclusions

The application roadmap does not currently activate B5/B6 all-at-once
execution, Stokes applications, stabilization comparisons, automatic OtD
derivation, continuous-control box reproduction, remote execution, or a GUI.
Those remain separate decisions under the benchmark, problem-library, and
implementation roadmaps.
