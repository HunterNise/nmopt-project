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
| Reports and post-processing | Implemented as tooling | Python reads persisted `artifact.kv`, solver traces, and VTU files; it produces field plots, comparisons, post-process indexes, and deterministic benchmark reports. The renderer uses `turbo` and explicit field-extrema colorbar endpoints; B1 visual validation is complete under its frozen contract, and B2 now has a calibrated source-raster comparison while its forward-state interpretation remains open. |
| B1 source-oriented execution | Reproduction-verified | The selected seven-case release profile combines continuous P1 control, a regular triangular mesh, constant-half forcing, and the common early-stop policy. Revision `631537a` completed all seven artifacts and their PNG comparisons without failures; the selected replacements do not claim recovery of omitted source data. |
| B2 framework-native execution | Framework-verified; forward replication audit open | A refreshed four-artifact `release-dealii` matrix at revision `df50946` contains current derivative, objective, field, and manifest evidence. The completed joint campaign documents non-reproducibility for the tested interpretations; a no-control forward-state audit now precedes any further optimizer work. |
| Parameter files | Implemented for registered B1/B2 slice | One ordered schema registry and benchmark adapters drive matrix expansion, exact exclusions, typed capability resolution, solver policy, run layout, and post-processing provenance. B3–B6 extension contracts are tested but not registered as executable benchmarks. |
| Later Chapter 6 benchmarks | Planned | B3/B4 are the next selected benchmark families after B1/B2 reproduction is resolved; B5/B6 remain later. |
| Complete Chapter 5 recipe library | Planned | The selected recipes used by B1/B2 exist, but the reusable recipe families listed in the problem-library roadmap are not all implemented. |

The current B1 authoritative manifest records the seven unique combinations
shown across Figures 6.2–6.3. All cases terminate at the common relative
gradient threshold: steepest descent takes `64/547/2175` iterations and
L-BFGS takes `2/4/4` iterations for
$\beta=10^{-1},10^{-2},10^{-3}$, while the additional L-BFGS
$\beta=10^{-6}$ field case takes four iterations. The source-sized report,
native fields, and PNG comparisons satisfy the frozen B1 contract; the
[replication findings](../benchmarks/b1-replication.md) retain the unresolved
source omissions and replacement rationale.

The refreshed B2 authoritative evidence records a complete four-case release
matrix with current derivative, objective-reduction, state/control, and
manifest fields. A subsequent 37-artifact release campaign completed without
manifest, derivative-check, or post-processing failures, but no tested joint
interpretation reproduced both the published no-control range and the four
optimized panels. These outcomes are evidence to preserve, not reasons to
silently change the frozen benchmark policy.

The historical development snapshot
`runs/chapter-6/b2/development/001/` was produced with refinement 6 and
contains the newer `b2.*` derivative, objective-reduction, and
state/control evidence. It remains useful for historical diagnosis but is not
the current source-sized acceptance evidence. The existing
`postprocessed/<case>/state.png` files
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
Their scientific adequacy is accepted for the frozen B1 comparison and remains
part of the B2 reproduction unit.

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

**Status:** implemented for B1/B2.

**Purpose:** Make every run auditable: identify what was selected, how it was
compiled and executed, what data was bound, and why a matrix entry failed.

**Work completed:**

- The deterministic `nmopt-benchmark-v1` artifact projection is implemented.
- Artifacts retain identity, provenance, compilation manifest, diagnostics,
  solver histories, benchmark fields, measurements, and selected fields.
- B1 records forcing and Hessian finite-difference evidence.
- Refreshed B2 source-sized artifacts record boundary form, case, provenance,
  derivative checks, objective/gradient reduction, and state/control norms.
- The report consumes the run manifest rather than discovering an arbitrary
  directory tree.

**Done when:** the source-sized artifacts for an activated benchmark contain
all evidence listed in its benchmark contract, and the report preserves
successes, failures, and missing evidence without inventing values. B0, B1,
and B2 meet this evidence boundary; B2's remaining question is numerical
interpretation rather than artifact completeness.

### A5 — Build and validate post-processing tooling

**Status:** implementation-complete; B1/B2 source comparison complete.

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

The calibrated B2 source-raster audit now covers orientation, decoded value
ranges, normalization, interpolation, colorbar semantics, and spatial-error
metrics. It demonstrates a numerical mismatch rather than a renderer-only
difference; diagnosing that mismatch belongs to A6.

**Done when:** the tool preserves field identity and geometry, documents its
rendering policy, produces deterministic comparisons from persisted inputs,
and a benchmark-specific visual check confirms whether a difference is
numerical or presentational. B1 meets this boundary under its frozen
replacement contract; B2 meets it through its audited mismatch report without
claiming reproduction.

### A6 — Correctly reproduce B1 and B2

**Status:** in progress; B1 reproduction-verified, B2 F0 provenance repaired,
coupled-scaling forward audit open.

**Purpose:** Establish whether the current framework-native fields and plots
represent the same mathematical quantities and visual conventions as the
book's B1/B2 figures.

**Work completed for B1:**

- source facts and omissions are separated from project replacement choices;
- raw and plotted field identities, extrema, orientation, point association,
  and the comparison-only $`p_{\mathrm{book}}=-p`$ sign are auditable;
- desired-state and forcing plots retain provenance and coordinate mapping;
- the `turbo` rendering, interpolation, shared normalization, and endpoint
  colorbar policy are explicit;
- the seven-case matrix includes the source-referenced
  $\beta=10^{-6}$ L-BFGS field case; and
- the release report and comparisons at revision `631537a` record seven
  successful artifacts and distinguish qualitative reproduction from the
  source details that cannot be recovered.

**Work completed for B2:**

- the source pages and Figure 6.5 raster are calibrated against native fields,
  including orientation, range, spatial correlation, and error metrics;
- the refreshed `release-dealii` evidence at revision `df50946` contains
  current derivative, objective, field, post-processing, and manifest data;
- zero-forcing target-transcription and forcing/target/association factorial
  campaigns separate objective-table clues from image-fit evidence; and
- the completed 37-artifact campaign documents non-reproducibility for the
  tested joint interpretations instead of selecting a fitted replacement;
- Follow-up F0 now has one correctly labelled Debug inlet-only diagnostic and
  six correctly labelled Release transport/ordinary provenance repairs, all
  with complete manifests, derivative evidence, native fields, and clean
  post-processing; and
- the F0 reconstruction script records seven complete exclusions with no
  pre-refinement forward-gate pass.

**Still needed for B2:**

1. Run the coupled interior-scaling ray in F1 at the prescribed Release
   values and apply the complete forward gate.
2. If F1 fails the shape gate, run the conditional F2 normalized load
   fingerprints or report the required weak-form contract.
3. Prepare, but do not implement without a separate decision, the F3
   forward-only serialization and scale-aware derivative-check proposals.
4. Require agreement in value range, decoded spatial shape, sign/trend, peak
   location, provenance, and a refinement pair before accepting a forward
   interpretation.
5. Only after F4 passes, compare initial objectives and the expected
   nested-observation ordering, then design the optimizer campaign.

**Done when:**

- any color inversion and palette difference have a documented cause;
- raw numerical fields and plotted fields have an auditable mapping;
- any framework and source sign conventions are explicitly separated;
- B1 and B2 plots use an explicitly documented normalization, endpoint ticks,
  and colormap policy;
- B1 input functions are plotted with their provenance and coordinate mapping;
- the B1 matrix includes the source-referenced $\beta=10^{-6}$ field case;
- the source-sized matrices are regenerated only after the interpretation is
  settled; and
- the benchmark report separates framework-native evidence, source comparison,
  and unresolved source omissions;
- the cause of the no-control discrepancy is either localized or bounded by
  explicit negative evidence; and
- optimizer conclusions are gated on a forward interpretation that passes the
  no-control checks.

B1 meets these criteria under its frozen replacement contract. B2 remains open
at the forward-state gate; its artifact and rendering evidence are complete.

The current diagnostic basis is the
[B2 replication report](../benchmarks/b2-replication.md). Continue from its
forward-state candidate table and preserve the distinction between source
facts, historically grounded hypotheses, and fitted controls.

### A7 — Add Deal.II-style parameter files

**Status:** implemented for the registered B1/B2 slice; B3–B6 extension
contracts remain unregistered.

**Purpose:** Allow an agent or user to reproduce a declared scenario from a
reviewable `.prm` file instead of relying on hidden defaults or a long CLI
command.

**Work completed:**

- One ordered schema registry drives `ParameterHandler` declaration and
  extraction, with benchmark-specific adapters that discover native scalar
  definition IDs from selectors and matrix axes.
- Typed binders resolve product, execution, reduced-method, and extension
  capability IDs, then reject known-but-unsupported combinations before output
  creation.
- Generic run-set planning expands, filters, excludes, and records resolved
  combinations; manifests retain parameter paths, content hashes, effective
  values, and artifact coordinates.
- CLI selection, output, and refinement overrides have explicit precedence and
  are retained as provenance.
- Future B3–B6 extension records have contract coverage but no parameter
  adapters or executable registrations.

**Still bounded:**

- Only B1 and B2 currently have parameter schema and execution registrations.
- B1/B2 forcing and target data use native named scalar-definition subsections;
  the B1 desired state is parsed through the same scalar definition contract.
- B3–B6 are not runnable through the parameter-file or benchmark-ID boundary.

**Done when:** a checked-in `.prm` reproduces a development run, the
effective configuration is visible in the artifact, invalid combinations fail
before execution, and a benchmark can declare whether its parameters are
frozen or overridden. This criterion is met for the current registered B1/B2
boundary; extending it to B3–B6 remains part of A8.

### A8 — Implement the next Chapter 6 benchmarks

**Status:** planned; after the remaining A6 B2 reproduction audit and benchmark
contracts. A7 is complete for the current registered B1/B2 slice.

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

The next handoff is the A6 B2 coupled-scaling forward screen: use the existing
`release-dealii` build without rebuilding it, test only the prescribed F1
values, and do not resume optimization until a candidate passes the range,
shape, trend, peak-location, provenance, and refinement checks recorded in the
B2 replication report.

## Exclusions

The application roadmap does not currently activate B5/B6 all-at-once
execution, Stokes applications, stabilization comparisons, automatic OtD
derivation, continuous-control box reproduction, remote execution, or a GUI.
Those remain separate decisions under the benchmark, problem-library, and
implementation roadmaps.
