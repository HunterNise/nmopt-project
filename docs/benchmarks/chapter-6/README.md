# Chapter 6 benchmark specifications

This document freezes the selected B0, B1, and B2 benchmark contracts. The
[benchmark-suite roadmap](../../planning/chapter-6-benchmark-suite-roadmap.md)
owns their order and acceptance purpose; the
[numerical-examples reference](../../guides/chapter-6-numerical-examples.md)
owns the source equations and reported values; and the
[application scenario contract](../../applications/chapter-6/README.md) owns
the public construction path.

These specifications select the current framework-native scalar realization.
They do not claim discretisation parity with the source examples where the
source uses triangular meshes and the current adapter generates quadrilateral
meshes. Every run must retain the realized mesh dimensions and provenance in
its manifest.

## Freeze record

| Field | Frozen value |
| --- | --- |
| Artifact schema | `nmopt-benchmark-v1` |
| Source catalogue | `docs/guides/chapter-6-numerical-examples.md` |
| Source catalogue revision | `1aaefbe473f9941a89d1df36192251511c052933` |
| External source | *Optimal Control of Partial Differential Equations*, Manzoni, Quarteroni, Salsa; Chapter 6, pages 187–190 |
| Build profile | `release-dealii` for reproduction; smaller explicitly named profiles are development-only |
| Backend | serial deal.II |
| Execution | assembled reduced DTO |
| Framework revision | recorded per run from the repository revision used to build the executable |

The source revision identifies the repository catalogue that records the
selected equations and recovered data policy. It is not an invented edition
or revision of the external book. A future change to that catalogue must
update this freeze record before new comparison runs are accepted.

## Common B0 contract

The B0 boundary is an association and artifact contract, not another solver.
One frozen run produces one detached `BenchmarkArtifactT<Envelope>` and one
deterministic text projection. The detached envelope must retain:

- the complete `CompilationManifest`;
- the typed solver-policy snapshot and solver report;
- the `RunEnvironmentRecord`; and
- the validation diagnostics returned by compilation.

The serialized projection uses sorted, escaped `key=value` lines with the
classic locale. Its required groups are:

| Group | Required contents |
| --- | --- |
| `artifact.*` | Schema identifier and format version. |
| `identity.*` | Scenario, recipe, output, source reference/revision, build profile, artifact directory, deterministic flag, and requirements. |
| `provenance.*` | Framework revision, recipe revision, mesh provenance, mesh structural identity, runtime-data provenance, and environment fields. |
| `manifest.*` | Manifest schema, semantic problem, backend, execution, spaces, bindings, formulation, metric, solve policies, regions, and realized maps/spaces. |
| `diagnostics.*` | Validation status and every diagnostic category, component, capability, and remedy. |
| `solver.*` | Solver policy, stopping reason, objective/gradient/step histories, accepted iterations, line-search trials, state/adjoint/metric solves, Hessian actions, and direction resets. |
| `benchmark.*` | Frozen benchmark parameters and explicit replacement choices. |
| `measurements.*` | Timing and memory values only when collection was enabled. |
| `selected_field.*` | The caller-selected output inventory, in stable order. |

The existing in-memory writer already fixes the escaping, locale, identity,
diagnostic, measurement, and selected-field conventions. The executable
runner must add the retained manifest, solver, environment, and benchmark
fields to the persisted projection without changing their mathematical
meaning or introducing a second lowerer.

Artifact paths are deterministic and owned by the runner:

```text
<output>/chapter-6/b1/<run-slot>/<method>/beta-<value>/artifact.kv
<output>/chapter-6/b2/<run-slot>/<case>/artifact.kv
<each artifact directory>/solver-trace.csv
<each artifact directory>/native/mesh-volume.vtu
<each artifact directory>/native/mesh-volume.svg
<each artifact directory>/native/fields-volume.vtu
<each B2 artifact directory>/native/control-boundary.vtu
```

The B1 `<method>` values are `steepest-descent` and `l-bfgs`; the B2 `<case>`
values are `wings-constant`, `full-constant`, `wings-parabolic`, and
`full-parabolic`. The per-run `identity.output_id` must include the same
method/regularisation or case suffix used by the path.

`solver-trace.csv` is a structured Armijo trial trace. It records the outer
optimization iteration, trial index, step length, trial objective, actual
descent slope, sufficient-decrease bound, and the finite/slope/acceptance
predicates. It is intended for diagnosing line-search behavior; it is not a
field or visualization export.

`native/mesh-volume.vtu` is a standalone serial deal.II VTU export of the
volume mesh. `native/mesh-volume.svg` is a lightweight 2D `GridOut` preview of
the same mesh. `native/fields-volume.vtu` is a serial deal.II VTU export
readable by ParaView. B1 contains
the final state, adjoint, and cellwise control fields on the shared volume
mesh. B2 contains the final state and adjoint on the volume mesh, while its
facewise control is exported separately as `native/control-boundary.vtu` on the
controlled boundary faces. These are final-state artifacts, not per-iteration
output.

### Report generation

The repository-local `tools/chapter6_report.py` script creates a deterministic
summary from one or more runner output directories. It uses only the Python
standard library and does not rerun the solver. Use the repository Python
tooling convention when invoking it:

```bash
python3 tools/chapter6_report.py \
  --input runs \
  --output runs/chapter-6-report
```

The report directory contains `summary.csv`, `summary.md`,
`objective-history.svg`, and `armijo-trials.svg`. The plots use the persisted
solver histories and `solver-trace.csv` sidecars when present. The summary
also inventories native mesh/field files and identifies runs that predate
those sidecars or the ParaView fields, so missing evidence remains visible
rather than being silently reconstructed.

### Field post-processing

Render one artifact's native fields without modifying its authoritative
`native/` files:

```bash
python3 tools/chapter6_postprocess.py \
  --artifact runs/chapter-6/b1/authoritative/steepest-descent/beta-1e-1 \
  --output runs/chapter-6/b1/authoritative/steepest-descent/beta-1e-1/derived
```

The command writes PNG and SVG plots for every available state, adjoint,
volume-control, and boundary-control field, plus `postprocess.json`. The
default output is `<artifact>/derived`. It reads current files below `native/`
and also accepts historical `fields.vtu` and `control.vtu` files.

To process a complete run root, use `--input`:

```bash
python3 tools/chapter6_postprocess.py \
  --input runs/chapter-6 \
  --output runs/chapter-6-postprocessed
```

The output mirrors the artifact paths below the input root. Each processed
artifact receives its own `postprocess.json`, and the root receives
`postprocess-index.json` with success and failure records. Successful artifacts
also produce shared-scale comparison PNG/SVG files below
`comparisons/<benchmark-family>/`; B1 and B2 are kept in separate comparison
groups.

## B0 acceptance

B0 is ready for executable integration when:

1. the same detached record renders byte-identically on repeated writes;
2. field keys are unique, sorted, escaped, and locale-independent;
3. invalid identities, measurements, diagnostics, and field keys fail with a
   contract error;
4. the persisted record retains the manifest, solver report, policy snapshot,
   environment, and benchmark-specific fields; and
5. path selection and directory creation remain runner responsibilities.

## B1 — E6.5.1 distributed Laplace control

### Identity and source policy

| Field | Frozen value |
| --- | --- |
| Base scenario | `chapter-6.b1.distributed-laplace` |
| Recipe | `chapter-5.scalar-diffusion-reaction-volume` |
| Source reference | E6.5.1, equation (6.64), Figures 6.2–6.3 |
| Forcing | Manufactured zero forcing, explicitly replacing the unspecified source forcing |
| Target | $z_{d}(x)=10x_{1}(1-x_{1})x_{2}(1-x_{2})$ |
| Domain | $(0,1)^{2}$ |
| Equation | $-\Delta y=f+u$, homogeneous Dirichlet boundary |
| Control | Cellwise `FE_DGQ(0)` volume control with positive cellwise $L^{2}$ metric |
| Regularisation | $\beta\in\{10^{-1},10^{-2},10^{-3}\}$ |
| Methods | Steepest descent and limited-memory BFGS from zero control |

The $\beta=10^{-6}$ field illustration is not part of the acceptance matrix.
Absolute objective values are not compared with the source because the
forcing is manufactured. The artifact must record
`chapter-6.e6.5.1.manufactured-zero-forcing` as the forcing provenance.

### Discretisation and solver policy

The reproduction mesh is the current adapter's generated unit-square mesh
with `GridGenerator::hyper_cube` and `refine_global(7)`. It uses the current
assembled `FE_Q` state realization and cellwise control realization. This is
a framework-native, source-scale run, not a claim that the source's triangular
mesh has been reproduced. The manifest must record active cells, dimensions,
finite elements, and `chapter-6.e6.5.1.framework-native-hypercube-r7` as mesh
provenance. Development runs may use a smaller refinement, but their artifact
must say so and cannot be used as reproduction evidence.

The frozen solver policies are:

| Method | Gradient tolerance | Armijo fraction | Backtracking | Maximum trials | Minimum step declaration |
| --- | ---: | ---: | ---: | ---: | ---: |
| Steepest descent | $10^{-3}$ | $10^{-5}$ | `0.7` | `5` | `0.2` |
| L-BFGS | $10^{-8}$ | $10^{-5}$ | `0.7` | `5` | `0.01` |

The minimum-step values are retained as benchmark declarations because the
current generic solver policy does not expose them as an active option. The
artifact must distinguish declared source policy from the solver behavior
actually executed.

### B1 acceptance evidence

The acceptance matrix contains six artifacts: three regularisation values
times two methods. Each artifact must contain:

- objective and tracking histories and final values;
- gradient, step, line-search, state, adjoint, metric, and direction counts;
- the compilation manifest and realized state/control dimensions;
- a finite-difference check of the linear-quadratic Hessian action; and
- the qualitative trend that smaller regularisation reduces cost/tracking and
  L-BFGS requires substantially less optimization work than steepest descent.

The trend is evidence, not a portable numerical tolerance. Failure to achieve
the trend must be reported as a benchmark limitation rather than hidden by
changing the frozen formulation.

### Current B1 status

The six source-scale framework-native runs execute with valid diagnostics and
show the selected qualitative trends: decreasing $\beta$ lowers the final
objective, and L-BFGS is substantially faster than steepest descent for
$\beta=10^{-1}$. The lower-$\beta$ runs reach the generic iteration limit or a
floating-point line-search floor before the configured gradient tolerance; this
is an execution-policy limitation, not a failed PDE compilation.

These runs are not absolute source-value reproductions because the frozen B1
forcing is manufactured zero forcing. The current artifact directories contain
the serialized `artifact.kv` records, but the source-scale records predate the
later `solver-trace.csv` and VTU sidecar exports, and they do not yet contain
the required finite-difference Hessian evidence. B1 is therefore a sound
framework validation path, but its acceptance handoff remains open until the
evidence set is refreshed or the missing evidence is supplied.

## B2 — E6.5.2 Graetz-flow boundary control

### Identity and source policy

| Field | Frozen value |
| --- | --- |
| Recipe | `chapter-5.scalar-neumann-convection-subdomain` |
| Source reference | E6.5.2, equation (6.65), Table 6.2, Figures 6.4–6.5 |
| Domain | $(0,4)\times(0,1)$ |
| Diffusion | $\mu=0.1$ |
| Reaction | `0` |
| Transport | $b(x)=(1.5x_{2}(1-x_{2}),0)$ through `conservative_transport` |
| Fixed data | Temperature `1` through `fixed_dirichlet_data` |
| Forcing | Zero forcing, with provenance `chapter-6.e6.5.2.zero-forcing` |
| Regularisation | $\beta=10^{-3}$ |
| Control | Facewise Neumann flux on `control_boundary`, no box constraint |
| Method | Full BFGS from zero facewise control |

The four frozen cases are:

| Scenario ID | Observation region | Target |
| --- | --- | --- |
| `chapter-6.b2.graetz-flow` | downstream wings | $2$ |
| `chapter-6.b2.graetz-flow.full-constant` | full downstream region | $2$ |
| `chapter-6.b2.graetz-flow.wings-parabolic` | downstream wings | $4x_{2}(1-x_{2})$ |
| `chapter-6.b2.graetz-flow.full-parabolic` | full downstream region | $4x_{2}(1-x_{2})$ |

The source boundary partition, recovered from Figure 6.4, is
$\Gamma_{D}$ as the left edge plus the first unit of the top and bottom walls,
$\Gamma_{c}$ as the remaining top and bottom wall, and
$\Gamma_{\mathrm{out}}$ as the right edge. The current B2 semantic graph and
deal.II adapter realize this as the fixed, control, and natural-outflow
regions with boundary IDs `0`, `1`, and `2`, respectively. The deal.II
compiler checks that these three declared sets are disjoint and cover every
exterior face.

The runtime provenance strings are fixed as follows:

| Binding | Provenance |
| --- | --- |
| Forcing | `chapter-6.e6.5.2.zero-forcing` |
| Desired state | `chapter-6.e6.5.2.target` |
| Fixed temperature | `chapter-6.e6.5.2.fixed-temperature` |
| Conservative transport | `chapter-6.e6.5.2.graetz-transport` |

### Discretisation and solver policy

The reproduction mesh is the current adapter's generated rectangle with
`GridGenerator::hyper_rectangle` and `refine_global(7)`. The mesh uses
material ID `1` for the selected downstream observation region and `0`
elsewhere, with the boundary labels created by the adapter. This is the
framework-native quadrilateral realization; the source's triangular mesh
counts remain comparison context only. The manifest must record actual cells,
state/control dimensions, material IDs, boundary IDs, and
`chapter-6.e6.5.2.framework-native-rectangle-r7` as mesh provenance.

The current B2 factory freezes the generic full-BFGS policy: initial step `1`,
Armijo fraction `1e-4`, backtracking factor `0.5`, maximum line-search trials
`20`, gradient tolerance `1e-8`, and no declared minimum step. These values
are framework policy, not recovered source values; the artifact must record
that the source step policy was unspecified.

The production B2 realization follows the source's ordinary-normal transport
boundary form

$$
\partial_n y-(b\mathbin\cdot n)y.
$$

The deal.II compiler retains the total conservative-transport conormal as an
explicit opt-in diagnostic alternative. The ordinary-normal realization scales
the control coupling by $\mu$ and assembles the corresponding state-dependent
outlet term. The typed B2 semantic policy, default runner, manifest, and
artifacts now record the ordinary-normal form; this choice is not a third
benchmark case.

### B2 acceptance evidence

The acceptance matrix contains four artifacts, one per case. Each artifact
must contain:

- state/control dimensions and the complete compilation manifest;
- transport, fixed/controlled/outflow boundary, and material-region records;
- residual JVP/VJP and reduced-Taylor evidence for facewise control;
- objective and relative-gradient reduction; and
- state/control behavior evidence showing the effect of changing observation
  region and target.

GLS and other stabilization policies are excluded. If the unstabilized
Galerkin realization is inadequate, the artifact and handoff must report that
limitation without changing the frozen scenario.

### Current B2 status

The existing source-scale artifacts compile and serialize valid results, but
they predate the boundary correction and all four stop at zero accepted
iterations with `line_search_failure`. They must be refreshed before serving
as post-fix evidence. The Debug refinement-1 rerun realizes the corrected
partition: both parabolic cases reach `gradient_tolerance` after seven
accepted iterations, while both constant-target cases reduce the objective
and then stop at `line_search_failure` after nine accepted iterations.

The solver traces show finite trial objectives and negative directional slopes;
the first B2 search direction is therefore not rejected as non-descent. The
remaining constant-target failures cannot yet be attributed solely to the
excluded unstabilized Galerkin realization. Corrected refinement-2 and
refinement-3 runs now fail before accepting an optimization step for all four
cases; their zero-control forward states grow markedly with refinement. The
refinement-2 contract diagnostic confirms that the ordinary-normal comparison
has a substantially smaller forward-state scale than the total-conormal
default, so the convention is numerically material rather than a control-unit
rescaling. The source-scale artifacts remain stale and have not been rebuilt.
The next B2 unit is to refresh the artifacts under the now-authorized
ordinary-normal reproduction convention.

## Activation gate

B0–B2 benchmark acceptance requires the contract tests, the complete six-run
B1 matrix, the complete four-run B2 matrix, and source-sized `release-dealii`
artifacts. Development runs can validate orchestration and artifact shape but
cannot close the benchmark gate. B3/B4, all-at-once B5/B6, preconditioning,
Stokes, stabilization, automatic OtD derivation, and continuous-control boxes
remain outside this freeze.
