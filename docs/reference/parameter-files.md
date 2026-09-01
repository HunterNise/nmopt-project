# Parameter files and plotting profiles

**Status:** implemented for the registered B1/B2 Chapter 6 runner slice;
future B3–B6 extension contracts exist but are not registered as executable
benchmarks.

This document defines the schema and resolution boundary for versioned
Chapter 6 experiment inputs. It separates the values that select and execute a
registered benchmark from the reusable visual policy that renders its persisted
fields.

The source catalogue and frozen benchmark contracts remain authoritative for
the mathematical meaning of B1 and B2. A parameter file is an executable
instance of those contracts or an explicitly marked development variation; it
does not silently replace the benchmark documentation.

## File roles and locations

Versioned inputs belong below `parameters/`:

```text
parameters/
  chapter-6/
    b1/
      authoritative.prm
      development/
        figure-6.3-book-policy.prm
        continuous-control.prm
        continuous-control-constant-one.prm
        continuous-control-structured-simplex.prm
        continuous-control-count-matched-simplex.prm
        figure-6.2-early-stop-constant-half.prm
        figure-6.2-early-stop-objective-matched.prm
        figure-6.3-constant-half.prm
        figure-6.3-objective-matched.prm
    b2/
      authoritative.prm
      development/
        figure-6.5-state-fit.prm
        figure-6.5-table-6.2-parabolic-fit.prm
        forcing-sweep.prm
        table-6.2-order-fit.prm
    plotting/
      chapter-6-b1.json
      chapter-6-b2.json
```

The `.prm` file describes a run set. Its matrix can expand to several
concrete artifacts. The JSON file describes reusable field and rendering
policy. The run manifest must retain the selected file paths, content hashes,
matrix axes, and resolved values so that generated evidence does not depend on
future edits to the source files.

Every tracked or experiment `.prm` file must be self-contained. `include` and
`INCLUDE` directives are prohibited. Repeat inherited settings explicitly in
each file, and use comments and provenance fields to describe the relationship
between related experiment families. Self-contained inputs keep runner,
post-processing, and archived run snapshots on the same portable configuration
contract.

Generated output remains below `runs/`; it must not become a second source of
configuration.

## Experiment parameter schema

The supported `.prm` sections are:

| Section | Owns |
| --- | --- |
| `Benchmark` | Scenario identity, recipe, source reference, and source revision. |
| `Matrix` | Explicit axes whose Cartesian product becomes the run-set. |
| `Problem` | Recipe-level semantic choices such as control representation, observation region, and bounds. |
| `Functions` | Registered forcing, target, fixed-data, and transport function selections and coefficients. |
| `Boundary` | Boundary IDs, regions, normal/conormal interpretation, trace, and face quadrature choices. |
| `Mesh` | Dimension, geometry, generator parameters, refinement, and mesh provenance. |
| `Compile` | State degree, observation discretization, execution, product, session ownership, and state/adjoint/control-metric solve policies. |
| `Solver` | Method, initial control, stopping rules, line search, iteration limits, and method-specific policies. |
| `Run` | Authoritative/development policy, build profile, output root, and harness behavior. |
| `Output` | Retained native fields and artifact output policy. |
| `Postprocessing` | Plot-style reference and matrix-axis binding for comparisons. |

### Schema registry and ownership

`ParameterHandler` declaration and value extraction are driven by one ordered
schema registry. The registry combines common run, mesh, compiler, solver,
output, and post-processing entries with the selected benchmark's adapter
entries. Each entry records its path, default, `ParameterHandler` pattern,
presence policy, and ownership class.

The ownership classes are:

- `consumed`: parsed into a typed scenario or run-set record and used by the
  selected execution path;
- `locked_profile`: retained as an explicit benchmark profile choice and
  checked by that benchmark's binder; and
- `provenance_only`: retained for audit but excluded from numerical ownership.

The current adapters register B1's `method` and `regularisation` axes and B2's
`regularisation`, `forcing`, `observation-region`, and `target-profile` axes.
Matrix expansion, filtering, exclusions, artifact coordinates, and manifest
construction are generic over those registered axes. Adding a future axis is a
benchmark-adapter change; it does not require a CLI, run-controller, or
manifest-schema branch. Only B1 and B2 currently have parameter schema and
execution registrations.

Typed resolution occurs after parsing and matrix expansion. Product and
execution IDs resolve to `ProductSelection` and `ExecutionSelection`; reduced
method IDs resolve to `ReducedMethod`; and the quadratic-KKT extension contract
provides typed KKT-method and identity-preconditioner selections. The B1
validator accepts only its registered reduced methods, assembled execution, and
reduced-DTO product. B2 accepts BFGS with the same assembled reduced-DTO
profile. Unknown IDs fail lookup, while known but unsupported selections fail
the benchmark capability validation before an output directory is populated.

The runner also exposes typed scalar lower/upper bound records for future box
constraints. Their scalar definitions are validated, but no B3/B4 bound
lowering or executable adapter is registered.

`Solver/maximum line search trials` counts all attempted steps, including the
initial one. `Solver/maximum backtracking reductions` is the source-facing
alternative and maps to one additional possible trial; a file must not set
both. `minimum step length` is operative, while the legacy
`declared minimum step length` remains provenance-only for stale benchmark
files.

`Solver/globalization` accepts `armijo` or `fixed-step` and defaults to
`armijo` for compatibility. It is global rather than method-specific. Armijo
uses the declared trial limit, fraction, backtracking factor, and optional
minimum step. Fixed-step instead uses `initial step length` as its exact
positive finite step and performs one unconditional finite-objective trial;
the Armijo-only entries do not become hidden acceptance conditions. Both
policies retain the existing iteration and stopping configuration unchanged.

Every shared reduced-solver field may be overridden in a method-policy
subsection (`Solver/method policy <method>`): `maximum iterations`, either
line-search trial/reduction count,
`gradient tolerance`, `stopping criterion`, the relative-gradient,
objective-change, and step tolerances, `initial step length`, `Armijo
fraction`, `backtracking factor`, and `minimum step length`. An empty method
entry inherits the corresponding `Solver` value. For the mutually exclusive
trial/reduction pair, specifying either method entry replaces the global pair.
The globalization, objective-target policy, and initial control remain global
because they define the common policy, relation, and start across the method
comparison.

An `objective target policy` of `explicit` consumes the numeric
`objective target`. The B1-only `match-reference-method` policy instead runs
the named `objective target reference method` first for every regularisation
value and passes its terminal cost to the dependent method after it satisfies
the selected stopping criterion. In `automatic` mode, any enabled tolerance
stop qualifies; numerical stationarity also qualifies for objective-change and
step-only stopping. Iteration-limit and line-search failures do not. Filtering
out the required reference artifact is an error rather than an implicit extra
run. L-BFGS method-policy subsections may also set `memory`, `curvature
tolerance`, and `initial inverse Hessian scaling` to `metric-inverse` or
`scalar-secant`.

B1 accepts `Problem/control representation` values `cellwise-volume` and
`continuous-volume-homogeneous-dirichlet`. The former is the authoritative
`FE_DGQ(0)` choice. The latter uses continuous `FE_Q` or `FE_SimplexP` at the
declared state degree according to the mesh family, shares the state's
homogeneous boundary region, and rejects a cellwise box. The effective
representation is recorded independently of the parameter-file provenance.

For B2, `Problem/control representation` accepts `facewise-constant` and
`continuous-nodal-trace`. The former is the frozen default with
`l2_facewise`; the latter uses a continuous degree-one trace with
`l2_neumann_trace` and rejects a facewise box. The
`Boundary/transport boundary form` entry independently selects the boundary
operator. The ordinary-normal value
`ordinary-normal-minus-transport` requires
`Boundary/conormal form = unspecified`. The diagnostic `total-conormal` value
requires `Boundary/conormal form = diffusion-minus-transport`. The latter
entry is a consistency declaration, not an independently selectable boundary
operator. Other values and incoherent pairs are rejected before execution.

B2 also requires `Compile/volume observation quadrature order` to be a
positive integer and accepts `analytic-quadrature` or
`state-fe-interpolation` for
`Compile/volume observation target realisation`. The frozen tracked files use
order `3` and analytic evaluation. Interpolation first realizes the target in
the scalar state finite-element space; it does not change the stated
$L^{2}$ observation functional or introduce an objective multiplier.

The B2 `Functions/target definitions` subsection declares each target profile
listed in the `Matrix/target-profile` axis. Each nested profile subsection is
a scalar-function definition with `id`, `kind`, `provenance`, and either
`value` or `expression`:

```text
subsection target definitions
  subsection constant
    set id = constant-2
    set kind = constant
    set value = 2.0
    set provenance = chapter-6.e6.5.2.target
  end
  subsection parabolic
    set id = parabolic-4*x1*(1-x1)
    set kind = expression
    set expression = 4.0*x1*(1.0-x1)
    set provenance = chapter-6.e6.5.2.target
  end
end
```

The selected matrix profile chooses which declared definition is evaluated;
all declared definitions are parsed and retained for a complete run-set
description.
The source defaults are `constant-2` with value `2` and
`parabolic-4*x1*(1-x1)` with expression `4.0*x1*(1.0-x1)`. Changing a value or
expression is an explicit development hypothesis, not a change to the frozen
B2 benchmark. The effective selected definition, kind, value, and expression
are written to B2 artifact evidence.

B1 and B2 read `Functions/forcing` as a stable definition ID. The selected
forcing is declared in the named `Functions/forcing definition <id>` subsection.
For B2, a populated `Matrix/forcing` axis owns the forcing IDs; otherwise the
direct `Functions/forcing` selector owns the one ID. Both forms use the same
declarative scalar records:

| `kind` | Required data | Realization |
| --- | --- | --- |
| `zero` | No value or expression | Scalar zero function. |
| `constant` | Any finite `value` | Scalar constant function. |
| `expression` | Nonempty scalar `expression` | deal.II `FunctionParser` using coordinates `x0`, `x1`, … and constants `pi` and `e`. |

Expressions use the deterministic scalar arithmetic and function syntax
supported by deal.II `FunctionParser`; semicolon-separated components and the
`rand`/`rand_seed` functions are rejected. Parser errors, unknown coordinates,
and coordinates beyond `Mesh/dimension` are rejected before assembly. For
example:

```text
subsection Functions
  set forcing = spatial-candidate

  subsection forcing definition spatial-candidate
    set kind = expression
    set expression = 0.4 + sin(pi*x0)*sin(pi*x1)
    set provenance = development.b1.spatial-candidate
  end
end
```

The checked `source-oriented-constant-half` ID is the authoritative replacement
selected after the B1 investigation. `manufactured-zero` and
`figure-inferred-constant-one` remain explicit development choices. None is a
special case in the runner: every selected ID resolves to its named subsection
and is lowered through the shared `ScalarFunctionDefinition` contract. B1's
`Functions/desired state` selector remains a direct scalar definition because
it is not a matrix axis; its `kind`, `expression`, and `provenance` entries are
parsed from the adjacent `desired state` subsection.

`Mesh/generator` defaults to `framework-native`, which consumes
`Mesh/refinement` and leaves the simplex entries unset. A simplex configuration
selects exactly one subdivision representation: positive `Mesh/subdivisions`
for one isotropic count, or `Mesh/axis subdivisions` for a comma-separated
positive count along each mesh axis. The per-axis list must have
`Mesh/dimension` entries. The runner parses both representations into the typed
mesh record. The B1 mesh constructor accepts only the isotropic form. B2
accepts `structured-simplex` and `centroid-split-simplex` on its rectangle
with per-axis subdivisions; its framework-native rectangle remains the frozen
default. The centroid-split choice additionally needs a positive
`Mesh/centroid splits` value and a deterministic `Mesh/selection seed`. For
per-axis counts $n_1,n_2$, the split count cannot exceed the $2n_1n_2$ base
triangles. It is a repeatable topology-sensitivity candidate, not a
reconstruction of omitted source connectivity.

B1 accepts `structured-simplex` with a positive isotropic subdivision count,
and `centroid-split-simplex` with positive subdivision and
`Mesh/centroid splits` counts plus a deterministic `Mesh/selection seed`. Both
simplex generators are two-dimensional, require zero global refinement, and
require B1's continuous homogeneous-Dirichlet control. A centroid split
replaces one base triangle by three, adding one vertex and two cells; the split
count cannot exceed the $2 n^{2}$ base triangles for isotropic count $n$.

The `Compile` section exposes separate maximum iterations, relative tolerance,
and absolute tolerance entries for the state, adjoint, and control-metric
solves. A zero state or adjoint iteration limit selects the compiler's
dimension-dependent rule; the control-metric limit must be positive. The
compiled manifest, rather than the requested values alone, is authoritative:
the current B2 target, for example, records direct UMFPACK state and adjoint
solves even though the shared scenario carries fallback iterative values.

The schema assigns typed patterns where a value has a common shape. Strings
such as function IDs and method names are registry selections, not arbitrary
code. Function expressions in the examples are declarative records consumed
by registered function constructors; they are not an embedded programming
language.

### Matrix expansion

Only entries under `Matrix` are expanded. The selected benchmark schema
adapter declares which matrix entries exist. If there is no `Selection`
section, all combinations of the declared axes are executed. A selection is
an optional filter for focused development or smoke work.

```text
subsection Matrix
  set method = steepest-descent, l-bfgs
  set regularisation = 1e-1, 1e-2, 1e-3, 1e-6
end
```

This declares eight B1 artifacts. A B2 file with two observation-region
values and two target-profile values declares four artifacts. A partial
selection can reduce the set:

```text
subsection Selection
  set method = l-bfgs
  # Omitting regularisation keeps all four beta values.
end
```

`exclude combinations` removes exact coordinates when a family is not a full
Cartesian product. Each bracketed coordinate must name every declared axis
exactly once and select declared values; semicolons separate coordinates:

```text
subsection Selection
  set exclude combinations = [method=steepest-descent,regularisation=1e-6]
end
```

Axis selections and CLI filters may narrow the remaining matrix, but they do
not disable exclusions. A selection that retains only excluded coordinates is
invalid because it resolves to an empty product. Duplicate, incomplete, or
unknown exclusion coordinates are rejected while reading the parameter file.

The resolver must validate every expanded combination before execution. Typed
registry lookups reject unknown product, execution, method, KKT-method, and
preconditioner IDs; benchmark capability validation rejects known but
unsupported combinations. It must also reject empty products, duplicate IDs,
and malformed coordinates before an output directory is populated.

The matrix uses stable IDs for categorical values. Numeric values are stored
in canonical form in the resolved configuration so that plotting and artifact
paths do not depend on local formatting. Independent axes should be declared
independently: B2 declares `observation-region` and `target-profile`, whose
Cartesian product is the four public cases.

### CLI precedence

The effective precedence is:

```text
parameter-file defaults
  → Selection axis filters and exact exclusions
  → explicit CLI selection filter
  → explicit CLI refinement override
```

With `--parameter-file`, the file supplies the benchmark, run kind, build
profile, output root, and matrix. An explicit `--output` changes only the
destination, while `--select` overrides the same matrix axis before exclusions
are applied. The `--refinement` option is a runner mesh override for supported
framework-native runs.

The CLI refinement override is intentionally retained for framework-native
smoke runs. It changes only the realized mesh refinement, updates mesh
provenance, and is recorded as an override in the manifest. It is rejected for
the subdivision-based B1 simplex generators, whose topology must instead be
changed explicitly in a parameter file. The override must not modify the
checked-in parameter file.

The output directory may remain a destination override, and the framework
revision should normally be recorded from the compiled executable or an
explicit provenance argument. Neither is a mathematical experiment choice.
All choices that affect the numerical run or retained evidence belong in the
parameter file.

The checked-in future-extension fixtures are contract tests only. They are not
accepted as runnable benchmark IDs and do not appear in `--list`.

## Plotting profile schema

Plot profiles are declarative JSON documents with schema ID
`nmopt-plot-v1`. They own reusable presentation policy:

```json
{
  "schema": "nmopt-plot-v1",
  "profile_id": "chapter-6.b2",
  "compatible_scenario_prefixes": ["chapter-6.b2.graetz-flow"],
  "defaults": {
    "colormap": "turbo",
    "normalization": "finite-extrema",
    "comparison_normalization": "shared-finite-extrema",
    "volume_interpolation": "gouraud",
    "volume_mesh_overlay": false,
    "colorbar_ticks": "endpoint-inclusive",
    "output_formats": ["png"],
    "dpi": 180
  },
  "fields": {},
  "axes": {},
  "default_comparison": {
    "rows": [],
    "columns": [],
    "group_by": []
  },
  "title_templates": {}
}
```

The profile must not assume that a particular number of artifacts implies a
layout. A comparison plan names matrix axes explicitly:

```json
{
  "rows": ["target_profile"],
  "columns": ["observation_region"],
  "group_by": []
}
```

The `.prm` file can override the default comparison axes. This allows the B2
style to be reused for a development forcing sweep with `columns = forcing`
and a one-row, three-column comparison.

If more than two axes vary and the plan leaves one unbound, post-processing
must either group by that axis or fail with a diagnostic. It must not silently
flatten it into an arbitrary panel order.

## Provenance requirements

The run manifest should retain at least:

```text
parameters.file
parameters.content_hash
parameters.selection
parameters.declared_matrix
parameters.excluded_combinations
parameters.resolved_combinations
plotting.profile_file
plotting.profile_content_hash
plotting.resolved_comparison
cli.refinement_override
```

The self-contained parameter and plotting documents should also be copied below
the run directory. Post-processing uses those snapshots by default; applying a
different profile is an explicit derived-output override.

## Examples

- [B1 authoritative parameter family](../../parameters/chapter-6/b1/authoritative.prm)
- [B1 Figure 6.3 solver-policy family](../../parameters/chapter-6/b1/development/figure-6.3-book-policy.prm)
- [B1 continuous-control constant-one-forcing candidate](../../parameters/chapter-6/b1/development/continuous-control-constant-one.prm)
- [B1 Figure 6.2 constant-half early-stop candidate](../../parameters/chapter-6/b1/development/figure-6.2-early-stop-constant-half.prm)
- [B1 Figure 6.2 objective-matched early-stop candidate](../../parameters/chapter-6/b1/development/figure-6.2-early-stop-objective-matched.prm)
- [B1 Figure 6.3 constant-half candidate](../../parameters/chapter-6/b1/development/figure-6.3-constant-half.prm)
- [B1 Figure 6.3 objective-matched candidate](../../parameters/chapter-6/b1/development/figure-6.3-objective-matched.prm)
- [B2 authoritative parameter family](../../parameters/chapter-6/b2/authoritative.prm)
- [B2 development forcing sweep](../../parameters/chapter-6/b2/development/forcing-sweep.prm)
- [B2 Figure 6.5 state-magnitude diagnostic](../../parameters/chapter-6/b2/development/figure-6.5-state-fit.prm)
- [B2 Figure 6.5 and Table 6.2 parabolic diagnostic](../../parameters/chapter-6/b2/development/figure-6.5-table-6.2-parabolic-fit.prm)
- [B2 Table 6.2 order-of-magnitude diagnostic](../../parameters/chapter-6/b2/development/table-6.2-order-fit.prm)
- [B1 plotting profile](../../parameters/plotting/chapter-6-b1.json)
- [B2 plotting profile](../../parameters/plotting/chapter-6-b2.json)

The runner accepts these files with `--parameter-file`; the copied snapshots
and resolved combinations below the run directory are the post-processing
default. Content hashes currently use a labelled deterministic FNV-1a-64
digest for drift detection; they are provenance values, not authentication.
