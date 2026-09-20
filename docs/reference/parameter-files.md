# Parameter files and plotting profiles

This reference explains how to configure an already-authored application
family through the current `nmopt_runner` and post-processing inputs.

There are two separate tracked configuration formats:

```text
Deal.II-style .prm
    numerical/application/run-family choices

JSON plotting profile
    derived presentation policy
```

They are intentionally separate. Changing a colormap should not change a PDE
or optimization run, and changing a solver tolerance should not require
editing plotting code.

For the C++ layer that defines which fields/choices a family exposes, see
[Application authoring](application-authoring.md). For run-set lifecycle and
artifact output, see [Application execution](application-execution.md).

## Source configuration versus run snapshots

The tracked input tree is:

```text
parameters/
├── chapter-6/
│   ├── b1/
│   │   ├── authoritative.prm
│   │   └── development/
│   └── b2/
│       ├── authoritative.prm
│       └── development/
└── plotting/
```

The repository distinguishes three roles:

```text
authoritative.prm
    accepted repository reproduction input for one family

development/*.prm
    explicit investigation / hypothesis family

runs/.../parameters.prm
    snapshot of the exact input used by one execution
```

`authoritative.prm` does not claim that every value was explicitly stated by
the source publication. It records the repository's accepted executable
realization, including explicit choices for omitted details.

Scientific/source interpretation belongs in the matching study documents, for
example [Chapter 6 scenarios](../studies/chapter-6/scenarios.md) and
[Chapter 6 benchmarks](../studies/chapter-6/benchmarks.md).

## How the runner selects a schema

The parser first discovers:

```text
Benchmark/id
Benchmark/recipe
```

and uses them to choose a `ParameterSchemaAdapter`.

The adapter contributes application-specific matrix axes and scalar-definition
catalogues to the common schema before the full file is parsed by deal.II
`ParameterHandler`.

That means `.prm` support is part of the application layer, not a free-form
map from strings to arbitrary C++ fields.

The current adapters expose:

```text
b1
    recipe:
      chapter-5.scalar-diffusion-reaction-volume

    possible matrix axes:
      method
      regularisation

b2
    recipe:
      chapter-5.scalar-neumann-convection-subdomain

    possible matrix axes:
      regularisation
      forcing
      observation-region
      target-profile
```

A concrete file can populate only the axes it needs, but at least one
`Matrix/*` axis must be nonempty.

If you are adding a **new** application family or parameter key, do not start
by editing a `.prm` file. First add the typed scenario field and application
binding described in [Application authoring](application-authoring.md), then
expose it through the schema adapter.

## Working path: define a run family

Start with stable application identity:

```text
subsection Benchmark
  set id = chapter-6.b1.my-development-family
  set recipe = chapter-5.scalar-diffusion-reaction-volume
  set source reference = local development experiment
  set source revision = working-copy
end
```

`Benchmark/id` and `Benchmark/recipe` are required and participate in schema
selection.

### Declare the experiment matrix

A parameter file describes a family of resolved scenarios:

```text
subsection Matrix
  set method = steepest-descent, l-bfgs
  set regularisation = 1e-2, 1e-3
end
```

The declared Cartesian product is:

```text
[method=steepest-descent, regularisation=1e-2]
[method=steepest-descent, regularisation=1e-3]
[method=l-bfgs,           regularisation=1e-2]
[method=l-bfgs,           regularisation=1e-3]
```

Matrix values are strings in `ParameterFile`; typed conversion occurs when the
resolved combination is bound to the application's scenario record.

Within one axis, values must be nonempty and unique.

### Narrow the family in the file

```text
subsection Selection
  set method = l-bfgs
end
```

Selection values must be a subset of the corresponding declared axis.

Multiple selected values are allowed:

```text
set regularisation = 1e-2, 1e-3
```

The declared matrix is still retained in run provenance; selection describes
the active subset.

### Exclude sparse coordinates

Use complete coordinate exclusions when the family is not a full Cartesian
product:

```text
subsection Selection
  set exclude combinations = [method=steepest-descent,regularisation=1e-3]
end
```

Each exclusion must name every declared matrix axis exactly once and use only
declared values.

This is intentionally different from selection:

```text
selection
    narrows one or more axes

exclusion
    removes a complete coordinate from the remaining product
```

After selections and exclusions, the resolved family must contain at least one
combination.

## Map `.prm` sections to the typed scenario

The current Chapter 6 binder resolves a file/combination into the same typed
scenario used by C++ callers:

```text
Problem         -> scenario.problem recipe choices
Functions       -> scenario.problem runtime definitions
Runtime         -> scenario.problem scalar runtime values
Boundary        -> scenario.problem boundary selections
Observation     -> scenario.problem observation selections
Mesh            -> scenario.compile.mesh
Compile         -> scenario.compile
Solver          -> scenario.solver
Run / Output    -> scenario.experiment
```

The parameter layer therefore does not bypass scenario validation. It is a
textual input path into the typed application API.

## Semantic problem choices

For a B1-style distributed control:

```text
subsection Problem
  set control representation = continuous-volume-homogeneous-dirichlet
  set cellwise box constraint = false
end
```

For a B2-style Neumann control:

```text
subsection Problem
  set control representation = facewise-constant
  set facewise box constraint = false
end
```

These strings are mapped to typed recipe enums before the `ProblemSpec` is
built. Unsupported combinations are rejected by the application validator.

For the underlying semantic choices, see
[Problem authoring](problem-authoring.md).

## Scalar and function-valued runtime data

Common scalar values use:

```text
subsection Runtime
  set diffusion = 1.0
  set reaction = 0.0
  set regularisation = 1e-2
end
```

Function-valued inputs use named definitions so numerical payload and
provenance stay together.

### Directly selected function

```text
subsection Functions
  set forcing = my-forcing

  subsection forcing
    set kind = expression
    set expression = sin(pi*x0)*sin(pi*x1)
    set provenance = development.my-forcing
  end
end
```

A constant definition is:

```text
subsection forcing
  set kind = constant
  set value = 0.5
  set provenance = development.constant-half
end
```

The current scalar-definition layer supports the zero/constant/expression
forms registered by the application schema.

`provenance` should identify where the chosen datum came from. Two functions
with the same coefficient values can still have different scientific roles or
source status.

### Catalogue selected by a matrix axis

A matrix axis can select one definition from a catalogue. B2 target profiles
use:

```text
subsection Matrix
  set target-profile = constant, parabolic
end

subsection Functions
  subsection target definitions

    subsection constant
      set kind = constant
      set value = 2.0
      set provenance = my.constant-target
    end

    subsection parabolic
      set kind = expression
      set expression = 4.0*x1*(1.0-x1)
      set provenance = my.parabolic-target
    end

  end
end
```

For each resolved `target-profile` coordinate, the binder selects the matching
definition.

The schema rejects an ambiguous scalar slot that is simultaneously driven by a
direct selector and a matrix catalogue.

## Boundary and observation selections

Application-specific semantic/runtime settings use their own typed sections.
For example:

```text
subsection Boundary
  set fixed id = 0
  set control id = 1
  set outflow id = 2
  set upstream transition = 1.0
  set transport boundary form = ordinary-normal-minus-transport
end
```

and:

```text
subsection Observation
  set material id = 1
end
```

These values are not global framework defaults. The selected application
binder decides what they mean and validates their compatibility with the
scenario/mesh.

## Mesh configuration

The common current mesh fields include:

```text
subsection Mesh
  set dimension = 2
  set lower = 0.0, 0.0
  set upper = 1.0, 1.0

  set generator = framework-native
  set refinement = 2

  set subdivisions = 0
  set axis subdivisions =
  set centroid splits = 0
  set selection seed = 0

  set provenance = development.unit-square-r2
end
```

Current generator names are:

```text
framework-native
structured-simplex
centroid-split-simplex
```

The typed scenario validator determines which fields are meaningful for the
selected generator/application.

`lower` and `upper` must have `dimension` entries and satisfy

$$
\text{lower}_{i} < \text{upper}_{i}
$$

for every coordinate.

The application-side backend adapter turns these typed options into the actual
deal.II triangulation. `.prm` parsing itself does not construct a mesh.

## Compiler configuration

A current Chapter 6 reduced compilation block looks like:

```text
subsection Compile
  set state degree = 1
  set execution = assembled
  set product = reduced-dto

  set state solve maximum iterations = 0
  set state solve relative tolerance = 1e-12
  set state solve absolute tolerance = 1e-14

  set adjoint solve maximum iterations = 0
  set adjoint solve relative tolerance = 1e-12
  set adjoint solve absolute tolerance = 1e-14

  set control metric solve maximum iterations = 1000
  set control metric solve relative tolerance = 1e-12
  set control metric solve absolute tolerance = 1e-14
end
```

These fields are first mapped to `chapter6::CompileOptions`; the backend
execution adapter then maps that record to `DealiiDiscretisationPolicy`.

The runner currently exposes the application-supported assembled/reduced-DTO
surface. The fact that lower compiler headers contain another enum value does
not automatically make it a valid `.prm` selection.

For the lower API, see [Compiler](compiler.md).

### Automatic inner iteration limits

The current state/adjoint scenario options allow a maximum-iteration value of
zero to request the application's automatic limit. It does **not** mean
"perform zero iterations."

The realized compiler manifest records the concrete policy actually used.

This interpretation is field-specific: the outer optimization maximum
iteration count must be positive.

## Solver configuration

A typical section is:

```text
subsection Solver
  set method = l-bfgs
  set globalization = armijo

  set initial independent control value = 0.0

  set maximum iterations = 200
  set maximum line search trials = 25

  set gradient tolerance = 1e-8
  set stopping criterion = gradient-norm

  set relative gradient tolerance = 0.0
  set objective change tolerance = 0.0
  set step tolerance = 0.0

  set objective target = none
  set objective target policy = none

  set initial step length = 1.0
  set Armijo fraction = 1e-4
  set backtracking factor = 0.5
  set minimum step length = 0.0
end
```

The current Chapter 6 application enums expose:

```text
method
    steepest-descent
    bfgs
    l-bfgs

globalization
    armijo
    fixed-step
```

This is narrower than the full solver API in
[Optimization](optimization.md). Adding a low-level solver policy to the
framework does not automatically make it configurable in every application.

### Method-specific overrides

When `method` is a matrix axis, per-method policy can override common fields:

```text
subsection Solver

  set maximum iterations = 5000
  set initial step length = 1.0

  subsection method policy l-bfgs
    set stopping criterion = relative-gradient-norm
    set gradient tolerance = 1e-30
    set relative gradient tolerance = 1e-3

    set memory = 5
    set curvature tolerance = 1e-14
    set initial inverse Hessian scaling = metric-inverse
  end

end
```

The binder resolves the effective field for the current method before
constructing the solver options.

For L-BFGS, the fields map directly to the lower-level semantics documented in
[Optimization](optimization.md): `memory` is the retained secant-history cap,
`curvature tolerance` guards $\langle y,s\rangle$, and the initial scaling
selects the metric inverse or its scalar-secant scaling.

### Backtracking reductions and trial counts

The schema accepts either:

```text
maximum line search trials
```

or:

```text
maximum backtracking reductions
```

for the effective method policy, but not both simultaneously.

If a source says "five backtracking reductions", the binder resolves that to
six total line-search trials: the initial trial plus up to five reduced trials.

This conversion is performed once at the configuration boundary so the lower
solver receives its native `maximum_line_search_trials` convention.

## Objective-target policies

The application solver record can select:

```text
objective target
objective target policy
objective target reference method
```

Current policy strings include:

```text
none
explicit
match-reference-method
```

An explicit policy parses a numeric `objective target` directly into the
solver options.

A matched-reference policy is different: it defines a dependency between
matrix coordinates. The run-set execution layer runs the reference method
first and supplies its converged objective to the dependent method.

That behavior is documented under
[Application execution](application-execution.md); it is not a hidden feature
of `ReducedSolverParameters`.

## Run and output configuration

```text
subsection Run
  set kind = development
  set build profile = debug-dealii
  set output root = runs
  set measure timings = true
end

subsection Output
  set retain fields = true
end
```

`Run/kind` is:

```text
reproduction
development
```

For `--parameter-file` execution, the file owns this value; the CLI cannot
replace it with `--run-kind`.

Current reproduction policy checks the declared build profile against the
compiled runner and requires the release deal.II profile.

`retain fields` controls whether the execution adapter writes authoritative
native fields. It does not choose which fields a plotting profile later
renders.

## Post-processing choices in `.prm`

The `.prm` file connects the run family to a tracked plotting profile and says
how matrix coordinates should form comparison grids:

```text
subsection Postprocessing
  set style profile = parameters/plotting/chapter-6-b1.json

  set comparison rows = method
  set comparison columns = regularisation
  set comparison group by = none

  set output formats = png
end
```

These axes refer to the declared experiment matrix/application metadata. They
are not arbitrary figure labels.

The numerical family remains in `.prm`; rendering details such as colormap,
field title, axis label, and history-figure style remain in JSON.

## CLI selection is an additional narrowing layer

Given:

```text
subsection Matrix
  set method = steepest-descent, l-bfgs
  set regularisation = 1e-1, 1e-2, 1e-3
end

subsection Selection
  set regularisation = 1e-1, 1e-2
end
```

a caller can further select:

```bash
nmopt_runner \
  --parameter-file family.prm \
  --framework-revision REV \
  --select method=l-bfgs \
  --select regularisation=1e-2
```

CLI selection is merged over the in-file selection.

It can only narrow declared axes to declared values. It cannot create a new
axis or add a value absent from the tracked family.

The effective selection and final resolved combinations are persisted in the
run manifest and snapshot.

## `ParameterFile` in C++

Runner code loads:

```cpp
using nmopt::application::runner::ParameterFile;

ParameterFile file =
  nmopt::application::runner::read_parameter_file(
    path);
```

The parsed record contains:

```cpp
file.path;
file.content_hash;
file.values;
file.matrix;
file.selection;
file.excluded_combinations;
```

Use:

```cpp
const std::string &id =
  file.value("Benchmark/id");
```

for a declared schema entry.

`optional_value()` is appropriate only when the application contract really
has a fallback:

```cpp
const auto group_by =
  file.optional_value(
    "Postprocessing/comparison group by",
    "none");
```

Application binders should still explicitly require values their typed scenario
needs rather than using `optional_value()` as an escape from validation.

### Resolve combinations

```cpp
const auto combinations =
  file.combinations();
```

or with additional CLI-style filters:

```cpp
const auto combinations =
  file.combinations({
    {"method", "l-bfgs"}
  });
```

Resolution proceeds in declared matrix-axis order, applies selection, removes
excluded coordinates, and rejects an empty result.

## Prefer `RunSetPlan` after parsing

The runner turns a parsed file into:

```cpp
const auto plan =
  nmopt::application::runner::make_run_set_plan(
    file,
    cli_filters);
```

The plan retains both requested and resolved structure:

```text
matrix_axes
selection
excluded_combinations
resolved_combinations
comparison rows/columns/group_by
parameter provenance
```

Each resolved combination also carries artifact-coordinate components.

Benchmark execution code should consume this plan instead of re-parsing matrix
strings independently.

## Plotting profile JSON

The checked-in plotting schema is:

```text
nmopt-plot-v1
```

A profile can specify:

```text
compatible scenarios
render defaults
field sources
axis order and labels
default comparison
history figures
title templates
```

A small profile is:

```json
{
  "schema": "nmopt-plot-v1",
  "profile_id": "my-profile",

  "defaults": {
    "colormap": "turbo",
    "normalization": "finite-extrema",
    "comparison_normalization": "shared-finite-extrema",
    "volume_interpolation": "gouraud",
    "volume_mesh_overlay": false,
    "colorbar_ticks": "endpoint-inclusive",
    "colorbar_tick_count": 5,
    "output_formats": ["png"],
    "dpi": 180,
    "axis_labels": ["x", "y"]
  },

  "fields": {
    "state": {
      "source": "state",
      "title": "State",
      "colorbar_label": "state"
    }
  },

  "default_comparison": {
    "rows": [],
    "columns": [],
    "group_by": []
  }
}
```

The current loader intentionally accepts a narrow set of rendering policies,
including:

```text
normalization              finite-extrema
comparison_normalization   shared-finite-extrema
volume_interpolation       gouraud
colorbar_ticks              endpoint-inclusive
output_formats              png / svg
```

Unsupported policy strings fail during profile loading rather than being
silently passed to matplotlib.

## Map logical plotting fields to persisted sources

A volume field is:

```json
"state": {
  "source": "state",
  "title": "State",
  "colorbar_label": "state"
}
```

A field on separate boundary topology declares:

```json
"control-boundary": {
  "source": "control",
  "mesh": "boundary",
  "title": "Boundary control",
  "colorbar_label": "control"
}
```

The JSON profile does not cause C++ execution to retain these fields. The
execution adapter must have written the corresponding native source first.

## Control comparison order and labels

```json
"axes": {
  "method": {
    "order": [
      "steepest-descent",
      "l-bfgs"
    ],
    "labels": {
      "steepest-descent": "Steepest descent",
      "l-bfgs": "L-BFGS"
    }
  }
}
```

An axis can also request labels from the matrix definition:

```json
"forcing": {
  "labels_from": "matrix"
}
```

The profile's default comparison layout can be:

```json
"default_comparison": {
  "rows": ["method"],
  "columns": ["regularisation"],
  "group_by": []
}
```

During normal run-root processing, this style metadata is combined with the
persisted `.prm` matrix/selection and run-manifest comparison plan.

## History figures

A profile can render persisted solver histories:

```json
"history_figures": {
  "convergence": {
    "plots": [
      {
        "source": "solver.objective_history",
        "title": "Objective",
        "x_label": "iteration",
        "y_label": "J",
        "x_scale": "linear",
        "y_scale": "linear",
        "iteration_origin": 0
      }
    ],
    "series": [
      {
        "where": {
          "method": "l-bfgs"
        },
        "label": "L-BFGS",
        "linestyle": "-"
      }
    ]
  }
}
```

Every history figure needs at least one plot and one nonempty series selector.
Optional `x_limits`/`y_limits` are increasing numeric pairs.

The `source` names persisted artifact histories. Post-processing does not have
an in-memory optimizer object.

## Profile and override precedence during post-processing

For a persisted run, `tools/postprocess.py` searches upward for:

```text
parameters.prm
plotting-profile.json
run-manifest.json
```

With no explicit override:

```text
plotting-profile.json
    supplies style

parameters.prm
    supplies matrix/comparison choices

run-manifest.json
    supplies persisted resolved provenance
```

An explicit:

```bash
--profile-file new-profile.json
```

replaces the snapshot style for that invocation.

An explicit:

```bash
--format svg
```

replaces output formats.

The post-processing manifest records the parameter source/hash, plotting
source/hash, snapshot provenance, explicit overrides, and the final effective
render/matrix/comparison selections.

## Hashes and historical snapshots

The runner records content hashes for the selected parameter file and plotting
profile in `run-manifest.json` and copies both files into the run root.

Post-processing recomputes snapshot hashes. If a snapshot was edited after the
run, its provenance record can retain the mismatch against the original
manifest hash.

For that reason, treat the run-root copies as historical evidence. To change an
experiment family:

```text
edit tracked parameters/... source
        ↓
commit/review the new input
        ↓
create a new run directory
```

rather than mutating an old run snapshot.

## From a tracked family to one resolved execution

The complete configuration flow is:

```text
tracked family.prm
      │
      ├── Benchmark identity selects schema adapter
      ├── Matrix declares candidate coordinates
      ├── Selection narrows axes
      ├── exclusions remove complete coordinates
      └── typed problem / compile / solver / run values are parsed
      │
      ▼
ParameterFile
      │
      ▼
RunSetPlan
      │
      ├── CLI --select narrows the declared family
      └── resolved combinations are validated
      │
      │  for each resolved combination
      ▼
application parameter binder
      │
      ▼
validated ScenarioT
      │
      ▼
application execution
      │
      ├── artifact.kv
      ├── native output
      └── run-manifest.json
      │
      ▼
run snapshot + plotting profile
      │
      ▼
report / post-processing
```

That is the intended high-level user experience: once an application family
has been authored in C++, a user can change supported problem data, numerical
policies, matrix coordinates, and presentation choices through versioned
configuration without learning the compiler's binding constructors or solver
template composition.

For the accepted Chapter 6 source/reproduction choices, continue with:

- [Chapter 6 scenarios](../studies/chapter-6/scenarios.md);
- [Chapter 6 benchmarks](../studies/chapter-6/benchmarks.md);
- [B1 replication](../studies/chapter-6/b1-replication.md);
- [B2 replication](../studies/chapter-6/b2-replication.md).
