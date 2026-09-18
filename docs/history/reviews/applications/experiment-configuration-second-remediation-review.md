# Application parameter ownership: second remediation review

- Status: implementation handoff
- Reviewed state: `codex/applications` at `361b0cc` (`refactor(applications): remove redundant b2 profile leaves`)
- Scope: B1 and B2 application configuration, runner binding, deal.II lowering, and experiment artifacts

## Outcome

The first remediation made scalar forcing and desired-state definitions more
honest, but the application layer still does not meet the experimentation
goal. Several numerical choices remain compiled into B1/B2, and several
`.prm` entries are accepted only when they equal those compiled values. Other
entries are parsed and recorded but do not affect execution.

The next remediation must establish this rule:

> A numerical value that defines an experiment and can reasonably be varied
> within an already implemented capability is owned by the `.prm` file,
> parsed into a typed value, consumed by the executable path, and reported in
> the run evidence. Changing it must not require rebuilding the framework.

This does **not** require another configuration language, JSON, a general
geometry parser, a general expression AST, compatibility syntax, or a broad
redesign of the compiler. deal.II already supplies the scalar and rank-one
function parsers needed here. The required framework work is limited to
binding typed data, preserving its lifetime, passing it to existing
application/backend construction points, and removing misleading options.

This review refines
[the ownership review](experiment-configuration-ownership-review.md) and
supersedes the narrower interpretation in
[the first remediation review](experiment-configuration-remediation-review.md)
that allowed routine geometry and vector-function values to remain compiled
when a profile name was present. A profile identifier is a useful artifact
coordinate; it is not a substitute for consuming the profile definition.
The earlier remediation remains applicable where it requires dynamic scalar
definitions, removal of duplicate values, complete evidence, and final
cleanup.

## Scope boundary

### The `.prm` file must own

- PDE coefficients and scalar or vector data used by a selected composition;
- domain bounds and experiment-specific geometric thresholds;
- externally meaningful boundary and material identifiers;
- observation-region definitions;
- mesh refinements, polynomial degrees, and selected quadrature orders;
- objective weights and regularisation values;
- solver tolerances and other numerical policies that the selected solver
  actually implements;
- the supported initial-control value;
- matrix axes, function-definition catalogs, output location, and operative
  run/evidence controls.

These values may still have defaults in low-level test fixtures. They must not
remain fallback sources for production runner scenarios after binding.

### Code should continue to own

- which residual, objective, metric, constraint, discretization, and solver
  capabilities exist;
- B1/B2 composition semantics and mathematical invariants;
- internal semantic region identities used to connect graph components;
- numerical safety tolerances used only to implement an algorithm, such as a
  coordinate-comparison tolerance or a derivative-check probe policy;
- capabilities not yet implemented, such as arbitrary spatial initial-control
  projection or a new stabilization scheme;
- source/book facts documented as provenance rather than exposed as choices.

The test for every declaration is therefore not merely “is this numeric?” but
“would an experimenter reasonably tweak this without changing the selected
capability?” If yes, it belongs in `.prm`. If no alternative is implemented
and the entry cannot affect execution, remove the entry rather than pretending
that it is configurable.

This boundary is also the useful preparation for B3–B6. Typed scalar
definitions, one rank-one vector definition, box coordinates, external region
IDs, and method-specific numerical policies can be reused by later
compositions without profile-name switches. Do not implement B3–B6 or invent
their unconfirmed parameter surfaces during this remediation; future-proofing
here means removing B1/B2 numerical closures at the existing composition
seams.

## Findings in the current tree

### 1. The removed B2 observation geometry was duplication, but the capability is still closed

The latest cleanup did not remove a working override. Before the cleanup, the
observation expressions were checked against expected strings while
`chapter6_b2.hpp` independently marked cells using the compiled conditions
`x > 1`, `y < 0.3`, and `y > 0.7`. Deleting the duplicate strings made that
particular state more honest, but left the experiment geometry hardcoded.

The runner also closes `Matrix/observation-region` to `wings` and `full` in
`apps/nmopt-runner/benchmark_binders.hpp`. Consequently, adding a third region
or changing a cutoff still requires source changes.

This can be fixed with the scalar function path that already exists. An
observation region is a scalar indicator evaluated at cell centers during
material tagging. It does not need a general geometry language. The selected
definition should be lowered through the existing scalar deal.II
`FunctionParser`; a positive/nonzero result marks the observed material.

### 2. B2 conservative transport is a checked label over a compiled vector function

`Functions/conservative transport` must currently be `graetz`, its definition
is checked against a fixed tuple-like string, and
`B2ConservativeTransportFunction` still compiles the coefficient `1.5` and the
parabolic profile directly into `chapter6_b2.hpp`.

deal.II already has the required parser:
`dealii::TensorFunctionParser<1, dim>` accepts one scalar expression per vector
component, separated using deal.II's semicolon convention. The native 2D
expression is, for example:

```text
1.5*x1*(1-x1); 0.0
```

The current parenthesized tuple spelling is not a reason to implement a tuple,
JSON, vector, or geometry parser. Migrate the parameter to deal.II expression
syntax and add only the small typed definition and lowering needed to own a
rank-one vector function.

The B2 fixed Dirichlet datum has the same, simpler problem: its value is read,
but its selector/kind is closed and the backend constructs a constant
function directly. It should use an ordinary `ScalarFunctionDefinition` so a
constant or expression can be changed from the file.

### 3. Domain and boundary numbers are declared but independently compiled

B2 declares the rectangle bounds, boundary IDs, upstream transition, and an
outflow coordinate. The binder checks them against the compiled construction
in `chapter6_b2.hpp` instead of passing them to it. B1 declares only a
`unit-hypercube` label while all mesh generators compile the interval
`[0, 1]`.

The minimal correction is a finite lower/upper point value in typed mesh
options. Use `ParameterHandler` list syntax, for example:

```text
subsection Mesh
  set lower = 0.0, 0.0
  set upper = 4.0, 1.0
end
```

Validate the component count, finiteness, and `lower[d] < upper[d]`, then use
the points in every selected B1/B2 mesh generator. This is ordinary structured
data, not an expression language.

For B2, bind and consume the fixed, control, and outflow boundary IDs plus the
upstream transition coordinate. Derive the outflow coordinate from
`Mesh/upper[0]`; `Boundary/outflow x` is redundant and should disappear. Keep
the geometric closeness tolerance internal. Keep internal graph-region IDs in
code, because they are semantic wiring rather than experimental numeric data.

The B2 classification shape may remain B2-specific. This remediation does not
need a general boundary-classification grammar.

### 4. B1 still has a compiled inventory of regularisation values

The parameter binder accepts a general positive regularisation value, but
`b1_beta_slug()` in `apps/nmopt-runner/main.cc` accepts only `1e-1`, `1e-2`,
`1e-3`, and `1e-6`. That table is used in artifact coordinates and scenario
IDs, so a `.prm` edit to `5e-4` fails before execution.

Generate the artifact coordinate from the stable, exact matrix value spelling
instead of mapping values through a compiled table. Existing authoritative
coordinates such as `beta-1e-3` must remain unchanged, while a new valid value
must naturally produce a coordinate such as `beta-5e-4`. Do not add more
cases to the table.

### 5. B1 desired-state data is consumed, but its identifier is still closed

The current direct definition reaches the backend, which is correct. The
binder nevertheless requires the identifier `b1-polynomial`. Remove that
identifier check. The selector value should be the file-owned semantic ID,
and the direct definition should supply its kind, data, and provenance.

The artifact must report the effective desired-state ID, kind, value or
expression, and provenance. Provenance alone is insufficient evidence that
the intended definition was executed.

### 6. Several solver entries are not honest executable choices

`Solver/initial control` is parsed as the string `zero` and recorded, but both
B1 and B2 executors always allocate a zero control vector. Within the current
capability, make this an explicitly numeric, uniform value for all independent
control coefficients. Construct the initial control with that value rather
than relying on zero initialization. A value of `0.0` preserves current runs.
Do not claim support for a spatial initial-control expression: projecting such
an expression is a separate capability and is outside this remediation.

The common solver schema also declares L-BFGS memory, curvature tolerance, and
initial inverse-Hessian scaling for every method. The binder stores those
values in the L-BFGS record even when B2 selects full BFGS. In that case the
entries are ignored, while full BFGS uses its own compiled curvature tolerance.

Declare solver parameters by selected method:

- L-BFGS owns memory, curvature tolerance, and its implemented initial
  inverse-Hessian scaling;
- full BFGS owns its implemented curvature tolerance, and the B2 construction
  must pass it to the solver;
- do not expose L-BFGS-only entries under full BFGS;
- do not expose full-BFGS inverse-Hessian scaling until that capability is
  implemented as a separate, deliberately selected unit.

### 7. The parameter surface contains read-and-ignore and exact-check leaves

The following entries should not survive merely for explicitness. A value
that can only repeat code is misleading, not useful documentation.

| Entry | Current behavior | Minimal disposition |
| --- | --- | --- |
| `Output/selected fields` | exact-checked while the backend writes a fixed complete set | remove; retain the operative all-or-none `Output/retain fields` |
| `Run/serialize artifacts` | parsed, but artifacts are always serialized | remove; serialization is a runner invariant |
| `Run/measure memory` | parsed, but no memory measurement is performed | remove until the capability exists |
| `Run/deterministic` | affects recorded identity, not execution | remove as a pseudo-choice; document determinism as an invariant if needed |
| `Compile/owned session` | validation requires `true` | remove; it is not a choice |
| `Compile/stabilization` | exact-checked as `galerkin`, with no alternative | remove; B2's supported composition remains Galerkin |
| `Boundary/normal orientation` | exact-checked as `outward` | remove; outward normal is part of the selected boundary operator |
| `Boundary/trace evaluation` | exact-checked to the sole implementation | remove |
| `Boundary/face quadrature` | exact-checked to the sole implementation | remove the label; keep the operative quadrature order |
| boundary/region name leaves | exact-checked while semantic names remain compiled | remove; keep internal semantic identity in code |
| `Boundary/conormal form` | duplicates the selected transport boundary form | remove and derive it from that typed selection |
| B1 `Problem/observation = full-domain` | exact-checked benchmark structure | remove |
| B1 `Solver/method = from-matrix` | matrix method is the real selection | remove |
| B2 `Functions/forcing = from-matrix` | matrix forcing is already the real selection | remove |
| `Mesh/geometry` for fixed box generators | exact label over typed bounds/generator | remove once bounds and the selected generator are authoritative |

Keep entries that do have an executable role: benchmark identity, recipe,
execution product, supported compile selections, run kind, build profile,
output root, timing measurement, retained-field toggle, and postprocessing
settings consumed by the manifest/tooling.

Do not turn this table into a generic “reject obsolete keys” facility. Migrate
the tracked parameter files and the declarations together. Normal
`ParameterHandler` declaration behavior is sufficient; no special include
checks, syntax linter, compatibility reader, or legacy-key framework is
needed.

### 8. Direct definitions and sweep catalogs need one predictable convention

Descriptive subsection names are valuable catalog IDs when several candidates
are present, but add an unnecessary second identifier for a single direct
definition. Use these two cardinality forms consistently.

For one definition selected directly by a port, the selector owns the ID and
the subsection uses the standard port name:

```text
subsection Functions
  set forcing = source-oriented-constant-half

  subsection forcing
    set kind = constant
    set value = 0.5
    set provenance = chapter-6-b1
  end
end
```

For a matrix-selected family, the matrix axis owns the selected IDs and a
plural catalog contains the definitions:

```text
subsection Matrix
  set forcing = zero, spatial-candidate
end

subsection Functions
  subsection forcing definitions
    subsection zero
      set kind = zero
      set provenance = experiment-baseline
    end

    subsection spatial-candidate
      set kind = expression
      set expression = 0.4 + sin(pi*x0)*sin(pi*x1)
      set provenance = experiment-candidate
    end
  end
end
```

Do not add `from-matrix`: it duplicates the matrix axis. Do not use progressive
numeric names such as `forcing-1`; stable semantic IDs survive reordered and
extended sweeps and make artifact coordinates meaningful. Apply the same
direct/catalog rule to forcing, desired state, target, observation-region, and
the two B2 boundary data ports where their cardinalities require it.

These are two schema shapes supported by the existing deal.II parameter
syntax, not two input languages and not legacy compatibility modes. The
descriptor needs only a direct definition path or a catalog prefix. Migrate
all tracked files atomically, then remove the displaced paths; do not maintain
both conventions indefinitely.

### 9. Production binding still starts from compiled benchmark defaults

The runner creates a B1 or B2 scenario using a default factory and then patches
it from the parameter file. This is safe only if the binder overwrites every
run-affecting field. At present, the hardcoded geometry and function paths show
that this condition is not met.

Default factories may remain for focused unit fixtures. Production binding
should be a complete typed B1/B2 construction step whose required inputs come
from the parsed file. A small benchmark-specific binder or factory is enough;
do not introduce a universal scenario builder solely for this cleanup.

The B2 production path must also stop converting observation IDs to the
closed `B2ObservationRegion` enum before binding. It should carry the selected
catalog ID and parsed indicator into material tagging. Enum helpers may remain
in local test fixtures if they no longer constrain production.

## Target parameter shapes

The exact section placement may follow the current descriptor layout, but the
following shapes express the required ownership and syntax.

### B2 conservative transport

```text
subsection Functions
  set conservative transport = graetz

  subsection conservative transport
    set kind = expression
    set expression = 1.5*x1*(1-x1); 0.0
    set provenance = chapter-6-b2-graetz-profile
  end
end
```

Implement only the rank-one vector definition required by conservative
transport: stable ID, deal.II expression, and provenance. The deal.II lowerer
must create and initialize `TensorFunctionParser<1, dim>` and retain it for as
long as the compiled application uses it. Keep scalar-expression restrictions
in the scalar lowerer; the semicolon-separated component syntax belongs only
to the vector lowerer.

### B2 observation regions

```text
subsection Matrix
  set observation-region = wings, full
end

subsection Observation
  set material id = 1
  set realization = cell-center-indicator

  subsection region definitions
    subsection wings
      set kind = expression
      set expression = x0 > 1.0 ? (x1 < 0.3 ? 1 : (x1 > 0.7 ? 1 : 0)) : 0
      set provenance = chapter-6-b2-wings
    end

    subsection full
      set kind = expression
      set expression = x0 > 1.0 ? 1 : 0
      set provenance = chapter-6-b2-full
    end
  end
end
```

Use the existing scalar definition and lowerer. Evaluate the selected function
at active-cell centers; reject non-finite results and mark the configured
observed material ID exactly when the result is greater than zero. Record the
selected definition and the realization rule in evidence. Do not invent region
primitives, boolean geometry nodes, or a second parser.

### Mesh and boundary values

```text
subsection Mesh
  set lower = 0.0, 0.0
  set upper = 4.0, 1.0
end

subsection Boundary
  set fixed id = 0
  set control id = 1
  set outflow id = 2
  set upstream transition = 1.0
end
```

The B2 application may retain its current boundary-classification algorithm,
but every numerical operand above must come from this typed data. The outflow
plane is `Mesh/upper[0]`; do not repeat it as another parameter.

### Initial control

```text
subsection Solver
  set initial independent control value = 0.0
end
```

Use a precise name that describes the supported realization. Populate every
independent control coefficient with the parsed finite value. Do not retain a
string selector named `zero`, and do not generalize this into spatial function
projection.

## Implementation plan

Each unit below should be reviewable on its own and should leave the tree in a
coherent state. The suggested commit boundaries are guidance for the
implementing agent; do not combine unrelated framework refactors with them.

### Unit 0 — Freeze current B1/B2 evidence before changing the schema

Outcome: a local comparison baseline for the current authoritative and
development profiles.

1. With the existing Debug deal.II runner, execute the smallest representative
   B1 cases covering both control discretizations and the existing four
   authoritative regularisation coordinates.
2. Execute B2's four observation-region/target-profile combinations and the
   development profiles used for figure/table evidence at the lowest useful
   refinement.
3. Preserve locally the effective parameter snapshots, scenario/combination
   IDs, material counts/measures, solver summaries, objective/gradient data,
   and field inventory required for comparison.
4. Record which hashes are expected to change because the parameter schema
   changes. Do not treat hash equality as the numerical-equivalence test.

This unit makes no tracked changes and needs no commit. Do not claim B1/B2
replication from a smoke run; it is only a migration baseline.

### Unit 1 — Canonicalize direct definitions and matrix catalogs

Outcome: one predictable scalar-definition convention without closed IDs or
sentinel values.

Primary files: `apps/nmopt-runner/parameter_files.hpp`,
`apps/nmopt-runner/parameter_binding.hpp`,
`apps/nmopt-runner/benchmark_binders.hpp`, the tracked B1/B2 `.prm` files, and
the existing parameter-file tests.

1. Represent each single definition under the standard subsection for its
   port. Take its ID from the selector.
2. Represent sweep candidates under a plural definitions catalog; take IDs
   from the corresponding matrix axis.
3. Remove `from-matrix` sentinels and the hardcoded `b1-polynomial` check.
4. Use the same discovery and validation rules for every affected scalar port.
   Require every selected ID to resolve exactly once and reject missing or
   duplicate definitions using the existing declaration/validation mechanism.
5. Migrate all tracked parameter files in the same unit and remove the old
   paths. Do not implement old/new compatibility.
6. Ensure effective evidence contains selected ID, kind, value/expression, and
   provenance.
7. Extend the existing parameter-file tests for a novel direct ID, a
   matrix-selected catalog, and a missing selected definition.

Suggested commit: `refactor(applications): canonicalize function parameter definitions`

### Unit 2 — Make B1/B2 mesh bounds and B1 matrix coordinates data-driven

Outcome: changing a supported box domain or B1 regularisation requires only a
`.prm` edit.

Primary files: typed scenario/mesh option definitions,
`apps/nmopt-runner/main.cc`, the B1/B2 binders,
`include/nmopt/application/dealii/chapter6_b1.hpp`,
`include/nmopt/application/dealii/chapter6_b2.hpp`, parameter files, and
existing application tests.

1. Add finite dimension-aware lower/upper coordinates to the smallest shared
   typed mesh record that both applications already consume.
2. Parse deal.II list syntax, validate size and ordering, and pass the points
   through every B1/B2 structured/simplex mesh path.
3. Remove exact geometry and bound checks once the typed points are consumed.
4. Replace `b1_beta_slug()` and its finite lookup table with a stable artifact
   coordinate derived from the exact matrix value text. Preserve the current
   four coordinate spellings.
5. Make production binding overwrite the mesh fields explicitly; no compiled
   factory default may determine a runner experiment after binding.
6. Add focused cases for nondefault valid bounds, invalid bounds, and a B1
   beta such as `5e-4`.

Do not add arbitrary mesh-generator plugins or a geometry DSL.

Suggested commit: `refactor(applications): bind mesh and b1 matrix numerics from prm`

### Unit 3 — Bind B2 scalar boundary data and conservative transport

Outcome: fixed boundary data and the transport vector are executable `.prm`
definitions.

Primary files: the smallest backend-neutral function-definition records,
deal.II lowering code, B2 binding/construction, B2 parameter files, and
existing B2 tests.

1. Reuse `ScalarFunctionDefinition` for fixed Dirichlet data. Remove the
   fixed selector/kind check and construct the backend function from the
   selected definition.
2. Add a narrowly scoped rank-one vector definition containing ID,
   deal-expression text, and provenance. Do not create a universal tensor
   hierarchy.
3. Lower it with `dealii::TensorFunctionParser<1, dim>` using semicolon-separated
   component syntax, validate dimension/parser errors at binding or lowering,
   and retain the parsed object for the compiled application's lifetime.
4. Remove `B2ConservativeTransportFunction` or reduce it to no production role
   after the parsed definition is authoritative.
5. Report both effective definitions in artifacts.
6. Extend existing B2 tests so changing the fixed datum and the transport
   coefficient changes the assembled/evaluated data without recompilation.

Suggested commit: `refactor(applications): lower b2 function data from prm`

### Unit 4 — Bind B2 observation and boundary geometry

Outcome: observation profiles and all routine B2 boundary numerics are
file-owned while classification semantics remain B2-specific.

Primary files: B2 typed options/binder, B2 mesh construction, parameter files,
artifact evidence, and existing B2 deal.II tests.

1. Replace the production observation enum/closed switch with the selected
   scalar indicator definition from the catalog.
2. Evaluate that indicator at cell centers during material tagging and apply
   the configured observed material ID.
3. Parse and consume fixed/control/outflow boundary IDs and the upstream
   transition coordinate. Validate that the external boundary IDs are
   distinct and suitable for the deal.II ID type.
4. Derive the outflow plane from the upper mesh bound. Remove the duplicate
   `outflow x` leaf and all exact checks of the old constants.
5. Thread the typed boundary IDs through mesh tagging, the B2 problem-spec
   regions/boundary sets, compiler inputs, structural validation, face-count
   evidence, and field metadata. Updating the three B2 regions after invoking
   the existing Chapter 5 recipe is sufficient; do not generalize every
   semantic reference-spec factory for this benchmark.
6. Keep internal semantic region IDs compiled and remove parameter leaves that
   only rename them.
7. Record the selected observation definition, realization rule, material ID,
   boundary IDs, transition coordinate, and mesh bounds in effective evidence.
8. Compare the authoritative `wings` and `full` material counts/measures with
   Unit 0, then test one new cutoff using only a `.prm` edit.

Do not make boundary classification generic and do not parse geometry beyond
the scalar indicator already supported by deal.II.

Suggested commit: `refactor(applications): bind b2 region and boundary numerics from prm`

### Unit 5 — Make supported solver numerical choices effective

Outcome: every declared solver numeric affects the selected solver, and the
initial control has one honest supported realization.

Primary files: solver parameter declarations/binding, typed solver options,
B1/B2 execution adapters, solver construction, parameter files, evidence, and
existing runner/solver tests.

1. Replace the `zero` initial-control string with the finite uniform
   independent-coefficient value and use it to build both B1 and B2 controls.
2. Keep initial-control feasibility validation where the selected constraint
   requires it.
3. Split method-specific parameter declarations and typed records. Do not
   parse L-BFGS values for full BFGS.
4. Add the full-BFGS curvature tolerance to its typed options and pass it to
   the full-BFGS solver construction. Add `1e-14` to the migrated B2 profiles
   to preserve the currently compiled behavior.
5. Confirm every remaining solver leaf is either used by the selected method
   or absent from that method's schema.
6. Extend existing tests with one changed uniform initial value and one changed
   full-BFGS curvature tolerance, asserting effective binding and behavior at
   the nearest stable seam.

Do not implement spatial initial-control projection or a new full-BFGS
inverse-Hessian scaling policy in this unit.

Suggested commit: `refactor(applications): make solver prm numerics effective`

### Unit 6 — Remove pseudo-options and close the default-then-patch gap

Outcome: the remaining parameter surface consists only of identity,
implemented selections, consumed values, and operative run controls.

Primary files: parameter descriptors, binders, scenario records, artifact
serialization, all tracked B1/B2 parameter files, and existing schema tests.

1. Remove the non-operative and exact-check leaves listed in Finding 7.
2. Remove the corresponding scenario fields and serialization only when they
   have no remaining consumer. Preserve evidence for actual selections and
   values.
3. Audit each B1/B2 production field from declaration to execution. The binder
   must either assign it from the file/selected definition, derive it from
   another authoritative field, or omit it because it is a fixed capability.
4. Keep benchmark factories usable for focused tests if helpful, but ensure
   the production runner does not silently inherit experiment numerics from
   them.
5. Update existing schema tests for the reduced surface. Do not add a special
   legacy-key rejection layer or a new CMake-only syntax test.

Suggested commit: `refactor(applications): remove non-operative prm surface`

### Unit 7 — Documentation, no-rebuild proof, and final cleanup

Outcome: the streamlined application layer proves the intended experiment
loop and contains no displaced compatibility or dead B1/B2 code.

1. Update authoritative parameter/profile documentation to describe the
   direct-versus-catalog convention, deal.II scalar/vector expression syntax,
   observation indicator realization, typed mesh/boundary values, and the
   supported uniform initial control.
2. Build the Debug runner once. Then, without rebuilding, create ignored,
   self-contained deal.II `.prm` variants under `runs/parameters` that change:
   - B1 regularisation to a value outside the old table, such as `5e-4`;
   - a B1 direct forcing or desired-state definition and its semantic ID;
   - B2 forcing or target data;
   - the coefficient in the B2 transport vector expression;
   - a B2 observation cutoff;
   - the B2 fixed boundary value;
   - one supported solver numeric, preferably the uniform initial coefficient.
3. Run those variants with the smallest useful mesh. Verify that the effective
   artifacts contain the edited definitions, hashes/coordinates distinguish
   them, and the corresponding numerical or material-tagging behavior changes.
4. Re-run the original B1/B2 profiles and compare them with Unit 0. Existing
   B1 artifact coordinates, B2 material measures/counts, field inventory, and
   numerical summaries must be preserved within the existing tolerances. New
   parameter hashes are expected after the schema migration.
5. Run `./build.sh pipeline debug-neutral` and
   `./build.sh pipeline debug-dealii` after all units are integrated.
6. Delete obsolete exact-check helpers, closed ID tables/enums no longer used
   by production, old parameter paths, redundant scenario fields, compiled
   function classes displaced by parsers, and compatibility branches. Keep
   only code needed by the intended application layer or by focused tests with
   a clear purpose.
7. Inspect the final diff for a smaller and more direct B1/B2 configuration
   path. Cleanup is complete only when there is one authoritative route from
   `.prm` declaration to typed value to backend consumer to evidence.

Do not run or regenerate the Release evidence as part of these units. The user
will rebuild Release after the remediation has no further code changes.

Suggested commit: `docs(applications): document prm-owned experiment workflow`

## Required preservation

The migration must preserve:

- all seven authoritative B1 matrix combinations and their current
  method/regularisation artifact coordinates;
- current B1 numerical behavior for the migrated values and the existing
  desired/forcing profiles;
- the four authoritative B2 observation-region/target-profile combinations;
- B2 development profiles used for the current figure and table comparisons;
- B2 material IDs, marked-cell counts/measures, boundary classification, field
  inventory, and numerical results for unchanged parameter values;
- existing B1/B2 semantic graph and compiler contracts;
- ordinary deal.II `.prm` parsing and self-contained parameter snapshots.

It must additionally enable, without a rebuild:

- an arbitrary valid B1 regularisation value;
- a new scalar definition ID and expression in a supported direct or catalog
  port;
- a changed B1/B2 box bound within the existing mesh capability;
- a changed B2 fixed boundary datum;
- a changed B2 conservative transport expression;
- a changed B2 observation indicator and cutoff;
- changed supported solver numerics, including the uniform initial control and
  full-BFGS curvature tolerance.

## Test strategy

Prefer extending the existing parameter, runner, B1, and B2 test executables.
The behavior belongs at their existing seams and does not justify a new test
framework or a parallel CMake convention.

The minimum focused coverage is:

1. parameter discovery resolves arbitrary direct IDs and every matrix catalog
   member;
2. the exact matrix spelling generates stable B1 artifact coordinates without
   a compiled value table;
3. mesh bounds reach every B1/B2 generator;
4. scalar fixed data and rank-one transport expressions evaluate to the
   parameter values;
5. observation expressions change material tagging while current expressions
   reproduce current counts/measures;
6. B2 boundary IDs and transition coordinate reach classification;
7. uniform initial-control coefficients reach the solver input;
8. full-BFGS curvature tolerance reaches full BFGS, while L-BFGS-only fields
   are not declared for it;
9. effective artifacts report all selected IDs and executable numerical data;
10. a single built runner executes several edited `.prm` variants without any
    framework rebuild.

Tests should verify values at the nearest typed or numerical consumer, not
only that parsing succeeds. A parameter test that accepts an expression while
the backend still uses a compiled class is not sufficient.

## Explicit non-goals

The implementing agent should not:

- introduce JSON or translate `.prm` files through JSON;
- create a new expression syntax, tuple syntax, include system, or preprocessor;
- write a general geometry/vector/tensor parser;
- add checks to reject include directives or other convention violations that
  are absent from the tracked files;
- preserve both old and new parameter layouts indefinitely;
- add profile-name switches for new numerical candidates;
- turn exact literal checks into larger exact literal tables;
- expose internal tolerances or semantic graph IDs merely because they are
  numeric;
- implement arbitrary spatial initial-control projection;
- add unrequested stabilization, mesh, solver, or output-field capabilities;
- introduce a universal application builder when complete B1/B2 binders solve
  the ownership gap;
- add a new CMake test solely to preserve a legacy parameter syntax;
- regenerate Release outputs before the code and parameter surface are final.

## Completion gate

The remediation is complete only when all of the following are true:

- every remaining B1/B2 `.prm` numerical experiment choice is consumed by the
  selected execution path;
- no production B1/B2 binder requires an experiment value to equal a compiled
  numerical literal;
- no supported numerical variation requires adding a C++ enum case, slug
  table entry, function subclass, or profile-name branch;
- direct definitions and matrix catalogs follow the single documented
  convention and use deal.II syntax;
- effective artifacts identify the selected definitions and values actually
  executed;
- existing B1/B2 profiles reproduce the pre-migration baseline;
- multiple new parameter-only variants run from one already-built executable;
- the final cleanup has removed displaced, duplicate, ignored, and
  compatibility-only application code.

That is the intended experimentation loop: edit a self-contained `.prm` file,
run the existing executable, inspect complete evidence, and iterate—without
changing or rebuilding the framework for routine numerical choices.
