# Parameter files and plotting profiles

**Status:** implemented for the Chapter 6 runner and post-processing pipeline.

This document defines the first schema slice for versioned Chapter 6
experiment inputs. It separates the values that select and execute a
benchmark from the reusable visual policy that renders its persisted fields.

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
    b2/
      authoritative.prm
      development/
        forcing-sweep.prm
    plotting/
      chapter-6-b1.json
      chapter-6-b2.json
```

The `.prm` file describes a run set. Its matrix can expand to several
concrete artifacts. The JSON file describes reusable field and rendering
policy. The run manifest must retain the selected file paths, content hashes,
matrix axes, and resolved values so that generated evidence does not depend on
future edits to the source files.

Generated output remains below `runs/`; it must not become a second source of
configuration.

## Experiment parameter schema

The proposed `.prm` sections are:

| Section | Owns |
| --- | --- |
| `Benchmark` | Scenario identity, recipe, source reference, and source revision. |
| `Matrix` | Explicit axes whose Cartesian product becomes the run-set. |
| `Problem` | Recipe-level semantic choices such as control representation, observation region, and bounds. |
| `Functions` | Registered forcing, target, fixed-data, and transport function selections and coefficients. |
| `Boundary` | Boundary IDs, regions, normal/conormal interpretation, trace, and face quadrature choices. |
| `Mesh` | Dimension, geometry, generator parameters, refinement, and mesh provenance. |
| `Compile` | State degree, execution, product, session ownership, and state/adjoint/control-metric solve policies. |
| `Solver` | Method, initial control, stopping rules, line search, iteration limits, and method-specific policies. |
| `Run` | Authoritative/development policy, build profile, output root, and harness behavior. |
| `Output` | Retained native fields and artifact output policy. |
| `Postprocessing` | Plot-style reference and matrix-axis binding for comparisons. |

`Solver/maximum line search trials` counts all attempted steps, including the
initial one. `Solver/maximum backtracking reductions` is the source-facing
alternative and maps to one additional possible trial; a file must not set
both. `minimum step length` is operative, while the legacy
`declared minimum step length` remains provenance-only for stale benchmark
files.

An `objective target policy` of `explicit` consumes the numeric
`objective target`. The B1-only `match-reference-method` policy instead runs
the named `objective target reference method` first for every regularisation
value and passes its gradient-tolerance terminal cost to the dependent
method. Filtering out the required reference artifact is an error rather than
an implicit extra run. L-BFGS method-policy subsections may also set `memory`,
`curvature tolerance`, and `initial inverse Hessian scaling` to
`metric-inverse` or `scalar-secant`.

B1 accepts `Problem/control representation` values `cellwise-volume` and
`continuous-volume-homogeneous-dirichlet`. The former is the authoritative
`FE_DGQ(0)` choice. The latter uses continuous `FE_Q` or `FE_SimplexP` at the
declared state degree according to the mesh family, shares the state's
homogeneous boundary region, and rejects a cellwise box. The effective
representation is recorded independently of the parameter-file provenance.

B1 accepts the registered `Functions/forcing` values `manufactured-zero` and
`figure-inferred-constant-one`. Their respective kinds are `zero` and
`constant`; the latter must declare `value = 1.0`. The constant-one selection
is a Figure 6.2 hypothesis with explicit inference provenance, not recovered
source data.

`Mesh/generator` defaults to `framework-native`, which consumes
`Mesh/refinement` and leaves the simplex entries at zero. B1 additionally
accepts `structured-simplex`, with a positive `Mesh/subdivisions`, and
`centroid-split-simplex`, with positive subdivision and `Mesh/centroid splits`
counts plus a deterministic `Mesh/selection seed`. Both simplex generators are
two-dimensional, require zero global refinement, and require B1's continuous
homogeneous-Dirichlet control. A centroid split replaces one base triangle by
three, adding one vertex and two cells; the split count cannot exceed the
`2 * subdivisions^2` base triangles. B2 retains only its framework-native
rectangle generator.

The `Compile` section exposes separate maximum iterations, relative tolerance,
and absolute tolerance entries for the state, adjoint, and control-metric
solves. A zero state or adjoint iteration limit selects the compiler's
dimension-dependent rule; the control-metric limit must be positive. The
compiled manifest, rather than the requested values alone, is authoritative:
the current B2 target, for example, records direct UMFPACK state and adjoint
solves even though the shared scenario carries fallback iterative values.

The exact `ParameterHandler` declarations will use typed patterns. Strings
such as function IDs and method names are registry selections, not arbitrary
code. Function expressions in the examples are declarative records consumed
by registered function constructors; they are not an embedded programming
language.

### Matrix expansion

Only entries under `Matrix` are expanded. If there is no `Selection` section,
all combinations of the declared axes are executed. A selection is an
optional filter for focused development or smoke work.

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

The resolver must validate every expanded combination before execution. It
must reject unknown values, empty products, duplicate IDs, and combinations
outside the registered benchmark capability.

The matrix uses stable IDs for categorical values. Numeric values are stored
in canonical form in the resolved configuration so that plotting and artifact
paths do not depend on local formatting. Independent axes should be declared
independently: B2 declares `observation-region` and `target-profile`, whose
Cartesian product is the four public cases.

### CLI precedence

The intended precedence is:

```text
parameter-file defaults
  → Selection section
  → explicit CLI selection filter
  → explicit CLI refinement override
```

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
parameters.resolved_combinations
plotting.profile_file
plotting.profile_content_hash
plotting.resolved_comparison
cli.refinement_override
```

The resolved parameter and plotting documents should also be copied below the
run directory. Post-processing uses those snapshots by default; applying a
different profile is an explicit derived-output override.

## Examples

- [B1 authoritative parameter family](../../parameters/chapter-6/b1/authoritative.prm)
- [B1 Figure 6.3 solver-policy family](../../parameters/chapter-6/b1/development/figure-6.3-book-policy.prm)
- [B1 continuous-control constant-one-forcing candidate](../../parameters/chapter-6/b1/development/continuous-control-constant-one.prm)
- [B2 authoritative parameter family](../../parameters/chapter-6/b2/authoritative.prm)
- [B2 development forcing sweep](../../parameters/chapter-6/b2/development/forcing-sweep.prm)
- [B1 plotting profile](../../parameters/plotting/chapter-6-b1.json)
- [B2 plotting profile](../../parameters/plotting/chapter-6-b2.json)

The runner accepts these files with `--parameter-file`; the copied snapshots
and resolved combinations below the run directory are the post-processing
default. Content hashes currently use a labelled deterministic FNV-1a-64
digest for drift detection; they are provenance values, not authentication.
