# Application execution and artifact reference

This is the agent-facing reference for generating, inspecting, and verifying
application runs. It records the current repository policies for benchmark
artifacts, run sets, native output, reports, and post-processing. It does not
describe recipe assembly or compiler internals; start with the
[application API reference](application-api.md) for that work.

The current benchmark choices and required evidence remain in the relevant
[benchmark specification](../benchmarks/chapter-6.md). The
[application roadmap](../planning/application-roadmap.md) owns mutable
implementation status and execution handoffs. This document records the
current interfaces and repository policies that agents need to consume.

## Current repository policies

| Concern | Current policy |
| --- | --- |
| Artifact schema | `nmopt-benchmark-v1`, rendered as deterministic escaped `key=value` lines. |
| Run-set root | `<output>/chapter-6/<benchmark>/<run-slot>/`. |
| Reproduction slot | `authoritative`, using the benchmark's declared mesh policy and `release-dealii`. |
| Development slots | `development/001`, `development/002`, and so on; smaller explicitly named meshes are allowed. |
| Run manifest | `run-manifest.json` at the run-set root, retained through success or failure. |
| Per-artifact record | `artifact.kv` below the artifact directory. |
| Derived output | `report/` and aggregate `postprocess/` below the run-set root; artifact-local derived files remain below that artifact. |
| Native output | Deal.II writes authoritative mesh and field files; Python tools only derive views. |

These policies are repository contracts, not mathematical source claims. When
the schema, layout, or tool boundary changes, update this document together
with the implementation and its focused contract tests. The roadmap records
the rollout and status of that change.

## Artifact and manifest contract

The public execution boundary is:

```text
typed scenario
  -> compiled problem and solver report
  -> detached experiment envelope
  -> benchmark artifact
  -> deterministic artifact.kv projection
```

For a Chapter 6 run, `BenchmarkHarnessT<Scenario>` projects scenario metadata
into a deterministic identity and `finalize(...)` creates a
`BenchmarkArtifactT<Envelope>`. The harness does not compile a PDE, execute a
solver, choose a path, or create a directory. The public types and ownership
rules are defined in the [application API reference](application-api.md).

The detached record must retain the following logical groups:

| Group | Required contents |
| --- | --- |
| `artifact.*` | Schema identifier and format version. |
| `identity.*` | Scenario, recipe, output, source reference/revision, build profile, artifact directory, deterministic flag, and requirements. |
| `provenance.*` | Framework revision, recipe revision, mesh provenance and identity, runtime-data provenance, and environment fields. |
| `manifest.*` | Semantic problem, backend, execution, spaces, bindings, formulation, metric, solve policies, regions, and realized maps/spaces. |
| `diagnostics.*` | Validation status and every diagnostic category, component, capability, and remedy. |
| `solver.*` | Solver policy, stopping reason, objective/gradient/step histories, accepted iterations, line-search trials, solve counts, Hessian actions, and direction resets. |
| `benchmark.*` | Frozen benchmark parameters and explicit replacement choices. |
| `measurements.*` | Timing and memory values only when collection was enabled. |
| `selected_field.*` | Caller-selected output inventory in stable order. |

`BenchmarkArtifactWriter` owns stable ordering, escaping, and locale-independent
rendering. The `nmopt_runner` application owns command-line selection, run-set
creation, artifact-directory creation, and file writing. Failed or unfinished
matrix entries remain listed in the run manifest with their diagnostic.

The structured `manifest.state_solve.*`, `manifest.adjoint_solve.*`, and
`manifest.control_metric_solve.*` fields retain each compiled algorithm,
preconditioner, resolved iteration limit, and relative/absolute tolerance.
These are realized compiler records: a direct solve is not mislabeled with an
iterative tolerance merely because the parameter file declared a fallback.
For B1, `benchmark.control_discretisation` records the stable effective
cellwise or continuous-control choice in addition to the compiled
`manifest.control_space` description. `b1.forcing_selection` records the stable
definition ID, `b1.forcing_kind` records `zero`, `constant`, or `expression`,
and the applicable `b1.forcing_value` or `b1.forcing_expression` records the
resolved payload. `provenance.forcing` retains the source/replacement
provenance.
`benchmark.state_dimension`,
`benchmark.control_dimension`, and `benchmark.adjoint_dimension` expose the
three executable-vector sizes. The corresponding
`benchmark.*_physical_dimension` fields count coefficients in the complete
finite-element fields, including homogeneous boundary values, while
`benchmark.*_independent_dimension` excludes constrained coefficients. This
distinction lets B1 comparisons test either interpretation of the source's
equal-dimension clue without changing the solver's coordinate convention.
The same dimension triplets are emitted for B2 state, control, and adjoint
fields. B2 additionally records the mesh generator, subdivision counts,
centroid-split count, and selection seed, along with mesh vertices and active
cells, total and per-region exterior-face counts, and the realized observation
measure.

Every compiled semantic space is projected as
`manifest.space.<semantic-id>.*`, including its realized dimension, runtime
role, finite element, and region. These deterministic fields expose the
structured manifest values directly; the human-readable
`manifest.state_space` and `manifest.control_space` descriptions remain
compatibility renderings.

B2 records the selected control as `benchmark.control_discretisation` and
`b2.control_discretisation`, while `manifest.control_metric_realisation`
records the compiled metric. Its objective evidence separates tracking and
control-regularisation values at the initial and terminal controls. The two
components must sum to the retained total objective. Its optimized and
uncontrolled state, control, and adjoint coefficient extrema are retained
beside the field files. The existing
`solver.gradient_norm_history` is the declared-metric gradient norm; the
explicit `solver.metric_gradient_norm_history` alias and the separate
`solver.coefficient_derivative_norm_history` make that convention auditable
without treating dual coefficient norms as Hilbert gradients.
The `solver.globalization` field records the scenario selection as `armijo` or
`fixed-step`, while `solver.policy` records the effective reusable policy as
`armijo` or `fixed_step`. Effective trial limits, reductions, and line-search
parameters come from the retained policy snapshot, so fixed-step artifacts
report one trial, zero backtracking reductions, and no Armijo fraction.

## Runner API

The backend-neutral public runner types are:

- `BenchmarkExecutionEvidenceT<Envelope>`: the execution callback's result,
  containing the detached envelope, validation diagnostics, optional
  measurements, selected fields, and additional artifact fields.
- `BenchmarkRunResultT<Envelope>`: the finalized
  `BenchmarkArtifactT<Envelope>` together with its rendered document.
- `HeadlessBenchmarkRunnerT<Scenario>`: the orchestration adapter that exposes
  `run(build_problem, execute)` and an `identity()` accessor.

The runner's `run(...)` sequence is fixed:

1. invoke `build_problem` with `scenario.problem` to obtain a
   `semantic::v1::ProblemSpec`;
2. invoke `execute` with that specification and the complete scenario to
   obtain `BenchmarkExecutionEvidenceT<Envelope>`;
3. record runner wall time when timing collection is enabled;
4. finalize the harness and render the artifact writer, returning
   `BenchmarkRunResultT<Envelope>`.

The problem builder and execution adapter own backend compilation, solver
invocation, envelope construction, and execution evidence. The runner owns
orchestration and artifact finalization; it does not lower PDEs, solve the
optimization problem, or select output paths.

## Run-set organization

The runner places every Chapter 6 run below the selected benchmark and run kind:

```text
<output>/chapter-6/<benchmark>/<run-slot>/
  run-manifest.json
  artifacts/
  report/
  postprocess/
```

The current B1 and B2 matrix layouts are:

```text
runs/chapter-6/b1/authoritative/
  artifacts/<method>/beta-<value>/
runs/chapter-6/b2/development/001/
  artifacts/<case>/
```

The run manifest is written before execution and updated after each artifact.
It records the benchmark, run kind, command, build profile, framework revision,
refinement override, expected artifact inventory, and `running`, `complete`, or
`failed` status. The selected build profile, framework revision, and realized
mesh remain metadata rather than path components.

The runner exposes `--list`, the compatibility `--benchmark` selector, and
the parameter-driven `--parameter-file` selector. Parameter files own the run
kind and matrix; repeatable `--select AXIS=VALUE` filters a declared matrix,
while `--output` and the optional `--refinement` remain destination/smoke
overrides. `--framework-revision` records executable provenance. Use
`nmopt_runner --help` for the complete current option surface.

Sparse experiment families use `Selection/exclude combinations` to remove
validated full coordinates from the declared Cartesian product. Exclusions
remain active under `Selection` and `--select` narrowing. The run manifest
records the declared matrix, exclusion coordinates, and final resolved
combinations separately.

For a versioned experiment family:

```bash
build/release-dealii/bin/nmopt_runner \
  --parameter-file parameters/chapter-6/b1/authoritative.prm \
  --framework-revision REV
```

The runner copies the parameter and plotting profile into the run directory,
records their content hashes and resolved combinations in the run manifest,
and retains the comparison-axis override for post-processing.

A B1 family may declare a matched objective target. In that case the runner
orders each regularisation's reference method before dependent methods,
requires the reference to satisfy the selected stopping criterion, and records
the supplying artifact in every dependent result. `automatic` stopping accepts
any enabled tolerance stop; numerical stationarity also qualifies for the
objective-change and step-only criteria. Iteration-limit and line-search
failures never supply a target. A selection containing only the dependent
method is therefore invalid.

## Native and derived outputs

Deal.II writes native outputs directly from the realized mesh and finite-element
fields. Current Chapter 6 artifacts use:

```text
artifacts/<method-or-case>/
  artifact.kv
  solver-trace.csv
  native/
    mesh-volume.vtu
    mesh-volume.svg
    fields-volume.vtu
    control-boundary.vtu       # B2 boundary control
  postprocess/                 # when an artifact is processed separately
```

`fields-volume.vtu` contains the final state, the native framework adjoint,
the comparison-only negative adjoint, and the B1 volume control. B1 also
exports its target and forcing functions sampled on the state DoFs. B2 stores
the optimized state, the zero-control state, both adjoint conventions, its
target and forcing functions, and a cellwise observation-domain mask. Its
control is stored in `control-boundary.vtu` because the boundary topology
differs from the volume topology. Facewise constants are cell data on separate
line cells; a continuous nodal trace is point data on connected line cells
whose unique points match the control coordinates. These are final-state and
input-function exports, not per-iteration output. The mesh SVG is a lightweight
2D preview; the VTU files retain the authoritative numerical topology and field
data.

`solver-trace.csv` records line-search trials for diagnostics. It is not a
field export. Reports and post-processing read persisted artifacts and native
files; they must not reconstruct missing authoritative values.

## Reports and post-processing

The deterministic report consumes one selected run manifest:

```bash
python3 tools/chapter6_report.py \
  --run-manifest runs/chapter-6/b1/authoritative/run-manifest.json \
  --output runs/chapter-6/b1/authoritative/report
```

The report retains pending, failed, and missing artifacts in its summary. The
current report outputs are `summary.csv` and `summary.md`. Field plots are
generated separately by the post-processing command; only the native mesh
preview remains SVG.

To render one artifact's native fields:

```bash
python3 tools/postprocess.py \
  --artifact runs/chapter-6/b1/authoritative/artifacts/steepest-descent/beta-1e-1 \
  --output runs/chapter-6/b1/authoritative/artifacts/steepest-descent/beta-1e-1/postprocess
```

The tool writes derived plots and `postprocess.json`. Run-root processing with
`--input` writes a root `postprocess-index.json` and comparison outputs. When a
run snapshot is present, the copied JSON profile and manifest comparison plan
are loaded automatically. `--profile-file FILE` explicitly selects another
profile; the compatibility `chapter6_postprocess.py` wrapper remains
available for legacy Chapter 6 inputs.

Both post-processing manifests retain a `provenance` object describing the
parameter and plotting-profile sources, their content hashes, explicit CLI
overrides, and the effective selections actually applied. The nested
`effective` record contains selected volume and boundary fields, resolved
matrix axes and combinations, the comparison plan, and output formats. This
keeps a derived plot set auditable without changing the runner's authoritative
artifact or run-manifest records.

The B1 plotting profile writes `comparisons/figure-6.3` for a run root with one
scenario. If an input root contains multiple scenarios, the scenario name is
retained as an intermediate folder; `group_by` values are retained as folders
in either case. Both panels use a logarithmic iteration axis with history
sample zero displayed at coordinate one; the objective axis is linear and the
gradient-norm axis is logarithmic. The profile selects only the three
regularisation values shown in the source figure, so additional development
artifacts do not silently alter the plot.

## Agent verification loop

After changing application, recipe, scenario, runner, or artifact code, use the
smallest relevant loop. The full build policy is in the
[build instructions](../../.agents/build.md).

For backend-neutral changes:

```bash
cmake --preset debug-neutral
cmake --build --preset debug-neutral
ctest --preset debug-neutral --output-on-failure
```

For deal.II application changes:

```bash
cmake --preset debug-dealii
cmake --build --preset debug-dealii --target nmopt_runner --parallel 1
ctest --preset debug-dealii --output-on-failure
build/debug-dealii/bin/nmopt_runner --list
```

For source-sized reproduction, use an existing release runner and omit the
refinement override so that the benchmark's declared mesh policy is used:

```bash
build/release-dealii/bin/nmopt_runner \
  --benchmark b1 \
  --framework-revision REV \
  --run-kind reproduction \
  --output runs
```

The direct runner command does not build the executable. Before rendering
derived output, check the Python environment used by the post-processing
tools:

```bash
python3 --version
python3 -c "import matplotlib, meshio; print('post-processing dependencies available')"
```

For a development smoke run after the runner is built:

```bash
tools/run_chapter6.sh \
  --benchmark b1 \
  --refinement 1
```

Do not present Debug or refinement-1 output as source-sized reproduction
evidence. Do not build `release-dealii` merely to inspect or document an
application; source-scale reproduction requires the explicit release policy
and permission described by the application roadmap.

## Routing by task

| Task | Read next |
| --- | --- |
| Assemble or change a recipe | [Application API](application-api.md), [Chapter 5 recipes](../applications/chapter-5.md) |
| Assemble or change a Chapter 6 scenario | [Chapter 6 scenarios](../applications/chapter-6.md), then the relevant benchmark contract |
| Generate or inspect a run | This document, then the relevant [benchmark contract](../benchmarks/chapter-6.md) |
| Check current execution status | [Application roadmap](../planning/application-roadmap.md) |
| Inspect compiler capability | [V1 semantic compiler](../implementation/v1/semantic-compiler.md) |
