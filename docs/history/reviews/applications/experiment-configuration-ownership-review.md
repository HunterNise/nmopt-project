# Application-layer experiment-configuration ownership review

## Status and decision

This document records a static application-layer review at `ec328551704c` on
2026-09-01. It is an implementation handoff, not an implementation batch: no
production code, parameter file, test, or generated run was changed as part of
this review.

The application layer has a useful executable B1/B2 vertical slice, but it
does not yet meet the intended configuration boundary. A `.prm` file can
change many values, but the runner and B2 adapter still own several candidate
names, function kinds, cross-products, geometries, and dispatch choices. The
result is that an experiment such as “try this B2 forcing expression” can
require a framework edit and rebuild even though the backend already supports
deterministic scalar expressions.

The required direction is:

> A `.prm` file owns run values, data definitions, selections, and matrix
> combinations within registered capabilities. The framework owns semantic
> ports, validation rules, compiler products, and concrete lowerers. Changing
> a value or expression within an existing capability must not require a
> rebuild; adding a genuinely new operator, discretization, function family,
> or solver capability still may.

B1 must remain the compatibility oracle. The refactor must fix B2 rather than
generalizing only the B1 path, and it must leave extension points that admit
B3–B6 without another B1/B2-shaped branch in the central parser or runner.
It must not implement B3–B6 in this batch.

## Authority and review scope

The review applies the ownership rules in the
[parameter-file reference](../../../reference/parameter-files.md), the
[application execution reference](../../../reference/application-execution.md),
the [application assembly API](../../../reference/application-api.md), and the
[composition boundaries](../../../design/composition-boundaries.md). B1–B6
requirements come from the
[Chapter 6 benchmark roadmap](../../chapter-6-benchmark-suite-roadmap.md), the
[Chapter 6 application contract](../../../applications/chapter-6.md), and the
[Chapter 6 benchmark contract](../../../benchmarks/chapter-6.md). Mutable
application status remains owned by the
[application roadmap](../../application-roadmap.md).

The inspected implementation surface was deliberately narrower than the full
branch history:

- `apps/nmopt-runner/main.cc`, `parameter_files.hpp`, and `runner.hpp`;
- `include/nmopt/application/{scenario,chapter6,scalar_function,catalog}.hpp`;
- the B1 and B2 deal.II application adapters;
- the B1/B2 application, runner, and parameter-file contract tests;
- the authoritative and representative development `.prm` families; and
- commit subjects, path statistics, and selected application-layer changes
  between `codex/main` and the reviewed head.

At the reviewed head, `codex/applications` is 12 commits ahead of
`origin/codex/applications` and 134 commits ahead of `codex/main`, with merge
base `90361b8eb590dba7cb9976b48c55917b3540f52d`. Those counts explain why the
review treats the accumulated application design rather than only the 12
unpublished commits.

## Executive assessment

The implementation currently follows this path:

```text
.prm
  -> one global ParameterHandler schema
  -> raw map<string, string> plus a hard-coded matrix inventory
  -> B1/B2 branches in main.cc patch factory-created scenarios
  -> B1/B2-specific runtime-data classes
  -> compiler and solver
  -> B1/B2-specific artifact paths and evidence augmentation
```

The target is a composed resolution path:

```text
.prm values and definitions
  -> one schema registry used for declaration and extraction
  -> generic matrix/run-set plan
  -> registered benchmark binder resolves one complete typed scenario
  -> shared data lowerers plus the selected compiler product and method
  -> generic run controller, with benchmark-owned evidence extensions
```

This preserves static typing at the recipe, scenario, compiler, and execution
boundaries. It does not require a universal untyped scenario or a new problem
class for every experiment.

| Finding | Severity | Required disposition |
| --- | --- | --- |
| `APP-R1` – the global schema and matrix inventory are closed over today’s B1/B2 names | P1 | replace with one registry and benchmark-owned schema extensions |
| `APP-R2` – B2 scalar forcing and target selection are benchmark-coded capabilities | P1 | use shared scalar definitions and one shared deal.II lowerer |
| `APP-R3` – several `.prm` entries are asserted against hard-coded adapter data instead of being consumed | P1 | consume supported data or expose it honestly as a locked named profile |
| `APP-R4` – runner dispatch, scenario construction, and artifact coordinates are closed B1/B2 branches | P1 | introduce runner-local benchmark registrations and a generic run-set plan |
| `APP-R5` – default-then-patch construction and duplicate fields obscure the effective source of truth | P2 | resolve complete typed records before scenario construction |
| `APP-R6` – tests prove known profiles but not the no-rebuild experimentation contract | P2 | add novel-value/expression and schema-consumption contract tests |
| `APP-R7` – documentation overstates the implemented parameter boundary | P2 | update the reference and roadmap only after the new gates pass |

## Boundary that the refactor must preserve

“Owned by `.prm`” does not mean that arbitrary text can manufacture an
unsupported mathematical capability. The implementation should use this
boundary consistently:

| Change | Owner | Rebuild expected? |
| --- | --- | --- |
| a new finite scalar constant or deterministic scalar expression for an already registered data port | `.prm` | no |
| a new regularization value, solver tolerance, matrix combination, exclusion, mesh size, or supported method selection | `.prm` | no |
| selection among already registered observation, control, product, or solver realizations | `.prm` | no |
| a new semantic residual/objective/constraint component or an unsupported data shape | framework | yes |
| a new deal.II lowerer, finite-element realization, compiler product, or solver algorithm | framework | yes |
| activation of a new benchmark contract | benchmark registration plus `.prm` | once at implementation time; later value experiments should not rebuild |

Unsupported selections must fail during configuration resolution with the
parameter path, requested value, benchmark ID, and available capability IDs.
They must not be silently ignored, accepted as provenance only, or fail later
inside assembly.

## B1 preservation contract

B1 is already reproduction-verified under its frozen project replacement
contract. The refactor must preserve all of the following:

- keep `parameters/chapter-6/b1/authoritative.prm` byte-for-byte unchanged
  unless a separately reviewed schema migration proves equivalent input and
  intentionally changes its recorded content hash;
- resolve the same seven unique method–regularization combinations and the
  same excluded steepest-descent `1e-6` coordinate;
- preserve the selected constant-half forcing, polynomial desired state,
  131-subdivision structured-simplex mesh, continuous homogeneous-Dirichlet
  control, solver policies, common initial control, and objective-target
  behavior;
- preserve run directories, artifact coordinates, manifest fields, native
  field names, solver-trace schema, comparison axes, and failure visibility;
- preserve the existing scalar-expression restrictions and B1 forcing
  evidence; and
- keep the established numerical and qualitative reproduction outcomes within
  the tolerances already owned by the B1 contracts. Do not use this refactor
  to retune B1 or regenerate a different authoritative profile.

The historical release record at framework revision `631537a`—seven
successful artifacts, parameter hash `fnv1a64:d9f046cb19ac4d2f`, and the
documented method/regularization trends—is the system-level comparison point.
It is not a request for bitwise equality across compilers.

## Selective history evidence

The application history shows repeated work in the same cross-cutting seam:

| Change | What it demonstrates |
| --- | --- |
| `fd0290a` and `c89522c` | parameter-file resolution and run-family expansion established the right external direction |
| `cff2561` | adding one inferred B1 forcing candidate crossed application, runner, parameters, tests, and documentation |
| `d46d1a9` | declarative B1 forcing required a second broad change before a scalar expression could be data |
| `18e005f` | method-specific policy variants again required schema, resolution, tests, and profiles to move together |
| `21982c6` | volume-observation realization is a legitimate framework capability addition, not merely a data change |
| `ec32855` | B2 targets gained scalar expressions, but through B2-specific target slots and adapter code |

This history supports a strict distinction. It is reasonable for a new
observation discretization or compiler product to touch the framework. It is
not reasonable for another scalar forcing or target expression to do so after
the scalar-expression capability exists.

## Findings

### APP-R1 – the parameter schema is a closed B1/B2 inventory

**Severity:** P1.

`apps/nmopt-runner/parameter_files.hpp` declares every entry in
`detail::declare_schema()` and repeats the same inventory in
`detail::read_values()`. Every declaration currently uses
`Patterns::Anything()`. Matrix and selection axes are hard-coded to `method`,
`regularisation`, `case`, `forcing`, `observation-region`, and
`target-profile`; observation profiles, B2 target slots, B2 forcing
definitions, and method-policy subsections are likewise pre-enumerated.

Consequences:

1. adding a field requires coordinated edits in declaration and extraction;
2. a typo survives schema parsing and is diagnosed only by a later ad hoc
   parser, if it is diagnosed at all;
3. a new axis needed by B3–B5—bounds, dimension, product, KKT method,
   preconditioner, or mesh level—requires central framework edits; and
4. candidate names such as `constant`, `parabolic`, `constant-one`, and
   `constant-two` become compiled vocabulary rather than `.prm` data.

The comment that `ParameterHandler` has no enumeration API explains the
current duplication but does not make it an appropriate long-term boundary.
The replacement needs one internal entry registry that both declares and
extracts values. Benchmark-specific entries and valid axes should be supplied
by benchmark configuration descriptors, while reusable data definitions must
not require the central schema to know candidate IDs.

Because arbitrary subsection names and deal.II `ParameterHandler` do not
naturally compose, the implementation must explicitly choose and document a
ParameterHandler-compatible representation for named data catalogs. Acceptable
directions include a single declared, losslessly parsed list/map entry or a
documented two-phase declaration mechanism. Do not solve this with a fixed
number of slots, a guessed list of names, or another B1/B2 global inventory.

### APP-R2 – B2 does not use the scalar-function capability already proved by B1

**Severity:** P1.

The backend-neutral layer already has `ScalarFunctionDefinition` with
`zero`, `constant`, and `expression` kinds. B1 lowers all three. B2 instead
stores:

- `B2ProblemParameters::ForcingSelection {zero, constant}` plus a separate
  `forcing_value`;
- `B2TargetParameters` with exactly two compiled members named `constant` and
  `parabolic`; and
- a four-value `GraetzCase` enum encoding the Cartesian product of observation
  region and target profile.

The deal.II adapters duplicate coordinate-name generation, identifier scans,
`FunctionParser` setup, constants, and random-function rejection. B2 forcing
uses a constant-only `B2ForcingFunction`; B2 desired state uses a separate
function class that switches on `GraetzCase` and owns both hard-coded target
slots. The runner mirrors those restrictions and rejects a B2 forcing kind
other than zero or constant.

This is the immediate B2 defect. A B2 forcing expression is within an existing
scalar data capability, yet it cannot be selected by editing the file alone.
The fix is not another B2 enum value. Replace B2’s selection/value pair with a
`ScalarFunctionDefinition`, resolve the selected target to the same type, and
lower both through one shared deal.II scalar-function factory. Keep semantic
role and provenance outside the generic function factory so artifacts can
still say whether a definition was forcing, target, or bound data.

### APP-R3 – the file sometimes declares values that execution does not own

**Severity:** P1.

`configure_b1_scenario()` and `configure_b2_scenario()` mix three different
operations: parsing values, validating a frozen benchmark, and asserting that
file text matches compiled constants. Examples include:

- B1 desired-state kind and expression are required to match a compiled
  polynomial function rather than lowered from the definition;
- B2 transport kind and expression are required to match a compiled Graetz
  tensor function;
- B2 mesh bounds, upstream transition, outflow coordinate, observation
  geometry, boundary IDs, region names, and selected solver are checked
  against fixed literals, while the mesh and material partition code still
  uses compiled coordinates and predicates;
- `Compile/execution` and `Compile/product` are declared as choices but the
  common resolver requires `assembled` and `reduced-dto`; and
- output selections are exact-string checks in B1 and are not a general field
  selection mechanism.

A frozen benchmark is allowed to constrain values. The defect is presenting
compiled profiles as general numeric/expression inputs. For every declared
entry, the new schema must classify it as one of:

- **consumed:** parsed into the typed scenario and used by execution;
- **locked profile:** a named registered realization whose internal constants
  are not separately presented as editable run choices; or
- **provenance-only:** explicitly documented and excluded from numerical
  ownership.

No run-affecting entry may be read and discarded or compared with a literal
while the adapter independently uses a different literal. In the immediate
B2 repair, forcing and target definitions must be consumed. Geometry and
transport should be either genuinely lowered in a later capability unit or
collapsed to honest named profiles such as the frozen Graetz geometry and
transport realization. Do not broaden this refactor into a generic vector
expression or geometry engine merely to avoid that honest lock.

### APP-R4 – central runner control flow is closed over B1 and B2

**Severity:** P1 for the B3–B6 handoff.

`CommandLineOptions` stores `run_b1` and `run_b2` booleans. CLI parsing accepts
only those two IDs. `main.cc` prepares and dispatches two branches, constructs
two scenario types, builds benchmark-specific artifact paths, and augments
evidence with benchmark-specific functions. It also carries B1/B2 parsing,
typed resolution, execution, and output orchestration in the same translation
unit.

That shape makes B3–B6 additions cumulative central edits. It also encourages
new configuration fields to be added where a branch can conveniently inspect
raw strings rather than where the typed benchmark scenario owns them.

Keep the public `ApplicationCatalog` metadata-only. Add a runner-local
registration layer composed from value records and callables, not an
inheritance hierarchy. A registration should connect:

- accepted CLI benchmark ID and default parameter path;
- benchmark-specific schema extension and typed combination binder;
- generic matrix coordinates and artifact slugs;
- a typed execution callback; and
- benchmark-specific evidence extension.

The run controller should operate on a common run-set plan and should not know
whether a coordinate is called `target-profile`, `bounds`, `product`, or
`reaction`. Adding an activated B3 benchmark may add a B3 registration and
typed binder, but it must not add B3 fields to a global `if` chain or matrix
inventory.

### APP-R5 – default-then-patch resolution creates hidden choices

**Severity:** P2.

Both execution paths call `make_b1_scenario()` or `make_b2_scenario()` and
then patch the returned defaults from raw parameter strings. Some defaults are
later validated, some are overwritten, and some remain effective without a
clear declaration at the call site. B2 also has both `Problem/initial control`
and `Solver/initial control`; only the former is required to equal `zero`,
while execution consumes the latter through solver options.

Production parameter loading should resolve one complete typed configuration
before scenario construction. Factory defaults may remain convenient for
unit fixtures, but the runner must be able to show which source supplied every
effective run choice. Remove duplicate effective fields; retain a compatibility
alias only if an existing tracked file depends on it, diagnose disagreement,
and record the canonical path in evidence.

### APP-R6 – current tests do not prove the experimentation loop

**Severity:** P2.

The existing contracts provide valuable evidence for matrix expansion,
specific B1/B2 profiles, B1 expression forcing, B2 target expressions,
boundary assembly, observation realization, artifacts, and run policy. They
mostly exercise names and combinations already compiled into the schema.

The missing acceptance test is behavioral: build the runner once, create a
self-contained development `.prm` containing a new ID and deterministic
expression absent from C++ source, execute it with the unchanged binary, and
verify that the resolved function and artifact evidence match the file. The
test must cover B2 forcing and should cover a selected target definition. A
source scan for the novel ID is a useful guard against accidentally adding it
to an enum or declaration table.

Tests must also enumerate all declared paths and prove that each one is
consumed, intentionally locked, or provenance-only. Unknown keys, unknown
capability IDs, duplicate data-definition IDs, malformed expressions,
dimension-incompatible coordinates, nondeterministic functions, invalid
matrix coordinates, and unused definitions selected by a combination must
fail before a run directory is populated.

### APP-R7 – parameter documentation is ahead of the actual boundary

**Severity:** P2.

The parameter reference correctly states that all choices affecting a run or
retained evidence belong in the parameter file. It also explicitly records
that named forcing definitions selected by a matrix are not yet generic. The
application roadmap’s status table calls parameter files implemented, while
its A7 work unit is still marked planned. The implementation is therefore
useful but only partially complete under its own ownership rule.

Do not resolve this by weakening the rule. After implementation, update the
reference with the exact generic definition and axis grammar, mark unsupported
capabilities honestly, and reconcile the roadmap status. The roadmap should
move only after the no-rebuild and B1 compatibility gates pass.

## Target application design

### One registry for parameter declaration and extraction

Introduce a small configuration schema value type, for example an entry with
path, default, deal.II pattern, required/optional policy, and ownership class.
The exact name is not prescribed. One ordered registry must drive both
`ParameterHandler::declare_entry()` and extraction into the parsed document.

Compose the registry from:

1. common run, output, post-processing, mesh, compiler, and solver entries;
2. reusable data-port schemas such as scalar function definitions; and
3. benchmark-specific entries and capability selections.

Do not put B3–B6 candidate values in the common schema. Patterns should reject
booleans, unsigned integers, finite numbers, and closed capability IDs as
early as deal.II permits; cross-field validation remains in typed resolvers.

### Typed resolved configuration before scenario construction

Separate these phases:

```text
parse document
  -> expand/filter matrix
  -> resolve one combination to complete typed records
  -> validate benchmark contract and registered capabilities
  -> construct scenario and executable data
```

Raw path lookup should end at the resolver boundary. The compiler adapter,
solver adapter, artifact writer, and benchmark execution callback must not
query `map<string, string>`.

Use common typed loaders for common records and separate B1/B2 binders for
their mathematical choices. Do not create one giant universal scenario with
optional B1–B6 fields. Composition of small records is the intended boundary.

### Shared scalar-data lowering

Move deterministic scalar lowering to a reusable deal.II application utility,
with one implementation of:

- zero and finite constant realization;
- dimension-derived coordinate names `x0`, `x1`, …;
- `pi` and `e` constants;
- deterministic expression parsing;
- rejection of `rand`, `rand_seed`, multiple components, unknown coordinates,
  and empty or inconsistent fields; and
- owner-bearing lifetime suitable for compiler data bindings.

Use that utility for B1 forcing, B2 forcing, and B2 target. Design the typed
input so B3 target and B4 lower/upper bound data can use the same definition
without depending on a B2 header. A later bound lowerer may sample the scalar
function into the registered cellwise bound representation; that is a
constraint capability, not work to hide inside the scalar parser.

### Independent axes and selected definitions

Retain B2’s independent `observation-region` and `target-profile` axes. Remove
the four-way `GraetzCase` as the configuration representation. If a compact
runtime case helper remains temporarily, construct it after the two axes have
been independently validated and do not use it to choose scalar data.

A named scalar-data catalog must be generic over IDs. The central schema may
know that a scalar definition has `id`, `kind`, `value`, `expression`, and
`provenance`; it must not know that the IDs are `constant`, `parabolic`,
`constant-one`, or `constant-two`. The selected combination resolves exactly
one definition per semantic port and carries its ID and provenance into the
artifact. Unselected definitions may remain in the self-contained run-set
snapshot but must not be instantiated as runtime functions.

### Generic run-set planning, typed execution

Create a backend-neutral run-set plan containing benchmark ID, normalized
matrix axes, resolved combinations, exclusions, comparison coordinates,
parameter provenance, and artifact coordinate components. The generic run
controller uses that plan for manifests and iteration.

Each benchmark registration then binds one combination to its concrete typed
scenario and execution adapter. This deliberately keeps execution typed while
removing benchmark names from CLI parsing, generic matrix expansion, manifest
construction, and artifact path generation.

## B3–B6 future-proofing requirements

The refactor is future-proof when the following future inputs fit the new
seams without another central schema redesign. It need not make the benchmark
executable yet.

| Benchmark | Future `.prm` ownership | Framework capability or registration |
| --- | --- | --- |
| B3 | manufactured scalar target, beta sweep, symmetric bound values, initial control, reduced/PDAS comparison coordinates, PDAS tolerances and outputs | B3 binder, registered cellwise box data, compiled PDAS product, evidence adapter |
| B4 | lower and upper scalar definitions, including spatial expressions, beta and optional dimension sweeps | bound-function-to-cellwise-data lowerer and B4 binder; reuse B3 KKT/PDAS path |
| B5 | beta, mesh, product, KKT method, preconditioner, tolerances, and source objective-scaling selection | registered quadratic-KKT product and method/preconditioner capabilities |
| B6 | the B5 run family with a changed reaction coefficient and corresponding provenance | reuse B5 binder/executor; no B6-specific runner or solver branch |

Specific design gates are:

- matrix planning cannot have a closed list of B1/B2 axes;
- scalar definitions cannot be named by B1/B2 enums;
- bound data must have a typed semantic port ready for constant or expression
  definitions, without pretending the bound lowerer already exists;
- product and method selections must resolve through capability registries
  rather than global `require_parameter(..., "reduced-dto")` checks;
- artifact coordinates and comparison plans must consume arbitrary registered
  axes; and
- B6 must be expressible as B5 plus a different registered operator/data
  selection, not as a new execution path.

## Implementation handoff

Implement the work in the following order. Each unit should be independently
reviewable and leave tests green before the next begins.

### Unit 0 – freeze characterization and B1 compatibility

**Purpose:** establish a safety net before changing configuration ownership.

**Actions:**

1. Add characterization tests for all seven authoritative B1 resolved
   combinations, their exclusion, effective forcing/target/mesh/control and
   solver records, artifact coordinate paths, selected fields, and manifest
   parameter provenance.
2. Add characterization tests for the four authoritative B2 combinations and
   the current development forcing sweep.
3. Record the effective typed scenario, not merely raw strings. Do not create
   brittle snapshots of unordered formatting or compiler-dependent floating
   output.
4. Preserve the existing B2 independent boundary residual oracles and
   assembly tests added around `34d8e34` and `b31a6c2`; they are mathematical
   gates, not configuration tests to rewrite.

**Likely files:** `tests/application/parameter_files_dealii_contract.cc`,
`tests/application/runner_contract.cc`, B1/B2 application contract tests, and
small test-only helpers.

**Exit gate:** all existing tests plus the new characterization tests pass,
and no production behavior or tracked `.prm` content changes.

### Unit 1 – replace the duplicate global schema inventory

**Purpose:** make valid paths and their ownership explicit in one place.

**Actions:**

1. Introduce the entry-schema value type and common/benchmark schema
   composition described above.
2. Generate declaration and extraction from the same registry; delete the
   parallel path lists in `declare_schema()` and `read_values()`.
3. Add typed deal.II patterns wherever a field has a local scalar shape.
4. Add ownership metadata (`consumed`, `locked profile`, or
   `provenance-only`) and a test that accounts for every registered entry.
5. Move matrix-axis declarations out of the common fixed array. Preserve the
   existing B1/B2 syntax through benchmark schema adapters during the
   migration.
6. Specify and document the dynamic named-definition representation before
   implementing it. It must be self-contained, deterministic, round-trippable
   in the run snapshot, and able to represent expressions containing normal
   function syntax without ambiguous delimiters.

**Likely files:** `apps/nmopt-runner/parameter_files.hpp`, a new small runner
configuration header if needed, and parameter-file contract tests.

**Exit gate:** one registry controls declaration and extraction; an added
benchmark axis or scalar-definition ID does not require editing a central list;
all existing files still parse to equivalent raw and normalized values.

### Unit 2 – introduce the shared scalar-function realization

**Purpose:** turn scalar data into a reusable application capability.

**Actions:**

1. Keep `ScalarFunctionDefinition` backend-neutral and add only fields needed
   by all scalar roles. Do not add B1/B2 labels to it.
2. Add a reusable deal.II factory under the application/deal.II boundary.
3. Move the duplicated B1/B2 coordinate, constant, expression, and randomness
   logic into the factory.
4. Migrate B1 forcing to the factory without changing its public evidence,
   expression semantics, or runtime ownership.
5. Migrate B2 target realization to the same factory while retaining the
   existing target artifact fields.
6. Add dimension, malformed-expression, nondeterminism, ownership, and
   zero/constant/expression contract tests directly against the shared
   factory.

**Likely files:** `include/nmopt/application/scalar_function.hpp`, a new
`include/nmopt/application/dealii/scalar_function.hpp`, the B1/B2 deal.II
headers, and their contract tests.

**Exit gate:** there is one deal.II scalar lowerer; B1 characterization is
unchanged; B2 targets no longer rely on a B2-specific parser implementation.

### Unit 3 – fix B2 configuration ownership

**Purpose:** make the motivating forcing/target experiments require only a
`.prm` edit.

**Actions:**

1. Replace `B2ProblemParameters::ForcingSelection` and `forcing_value` with a
   selected `ScalarFunctionDefinition`.
2. Resolve `Functions/forcing` through the generic definition loader and use
   the shared scalar lowerer. Remove `B2ForcingFunction` and the runner’s
   zero/constant switch.
3. Replace `B2TargetParameters {constant, parabolic}` with a generic catalog or
   resolved selected definition. Remove compiled target-profile member names.
4. Keep observation region and target profile independent. Bind the selected
   target by catalog ID, not through the four-way `GraetzCase` switch.
5. Ensure B2 artifact evidence records forcing and target selection ID, kind,
   value or expression, and provenance. Evidence must record the effective
   definition, not just the source file hash.
6. Audit every remaining B2 raw parameter check. Consume the field if the
   existing lowerer supports it; otherwise replace editable-looking leaf data
   with an honest registered-profile selection and a precise unsupported-value
   diagnostic.
7. Resolve the duplicate initial-control entries. Prefer
   `Solver/initial control` as the canonical executable field; accept the
   problem-level entry only as a temporary compatibility alias if required by
   existing tracked files.

**Likely files:** `include/nmopt/application/chapter6.hpp`,
`include/nmopt/application/dealii/chapter6_b2.hpp`,
`apps/nmopt-runner/main.cc`, parameter loaders, B2 `.prm` files only where a
documented migration is unavoidable, and B2/parameter contract tests.

**Exit gate:** after building once, a test-owned B2 development file can use a
new ID and expression such as `sin(pi*x0)*x1` for forcing, run with the
unchanged binary, and emit matching evidence. The frozen four-case B2 file
retains its mathematics and all independent boundary oracles pass.

### Unit 4 – separate run-set planning from benchmark execution

**Purpose:** stop central runner growth before B3–B6 activation.

**Actions:**

1. Replace `run_b1`/`run_b2` booleans with one benchmark ID in command-line
   options.
2. Introduce runner-local benchmark registrations for B1 and B2. Keep the
   public application catalog metadata-only.
3. Move B1 and B2 raw-to-typed binding out of `main.cc` into bounded binder
   units. The main translation unit should retain CLI setup, run policy,
   registry lookup, generic orchestration, and top-level diagnostics.
4. Build one generic `RunSetPlan` from registered axes, filters, exclusions,
   and comparison coordinates.
5. Generate artifact coordinate components from normalized axis records plus
   optional benchmark slug policies. Preserve current B1/B2 paths through
   explicit compatibility tests; do not infer special paths from matrix size
   or axis order.
6. Route product, execution, method, and preconditioner values through typed
   capability lookups. B1/B2 may still register only their implemented
   selections, but the common resolver must not hard-require one literal.

**Likely files:** `apps/nmopt-runner/{main,runner,parameter_files}.cc/.hpp`, new
small B1/B2 binder or registration units, `CMakeLists.txt`, and runner tests.

**Exit gate:** adding a test-only benchmark registration with a novel axis
requires no edit to CLI parsing, matrix expansion, manifest writing, or the
generic run loop; B1/B2 paths and evidence remain compatible.

### Unit 5 – add explicit B3–B6 extension contracts, not implementations

**Purpose:** prove the refactor has the necessary shape without expanding the
benchmark scope.

**Actions:**

1. Add compile-only or lightweight resolver fixtures showing that registered
   axes can represent B3 beta/bounds/method, B4 lower/upper definitions, and
   B5 product/preconditioner selections.
2. Add typed configuration records or extension ports for box-bound data,
   product selection, and solver-family parameters only where an implemented
   framework contract already exists. Otherwise add a deliberate
   unsupported-capability diagnostic and leave the executable adapter absent.
3. Demonstrate through a fixture that a B6-like configuration reuses the B5
   registration shape with a changed reaction coefficient rather than a new
   runner branch.
4. Do not add fake PDAS/KKT execution, generic geometry parsing, vector
   expressions, or continuous-control bound semantics.

**Likely files:** runner configuration contracts and tests. Changes to compiler
or solver headers are out of scope unless an existing public capability cannot
be selected without a narrowly documented adapter.

**Exit gate:** configuration resolution has named, typed places for future
bound, product, method, and preconditioner selections; unsupported execution
fails at capability resolution; no B3–B6 benchmark is advertised as runnable.

### Unit 6 – remove compatibility scaffolding and close documentation

**Purpose:** leave one comprehensible path after behavioral equivalence is
proved.

**Actions:**

1. Remove obsolete B2 enums, duplicated scalar parsers, dead exact-string
   checks, and compatibility aliases that no tracked file needs.
2. Update the parameter-file reference with the exact definition/catalog,
   axis, precedence, validation, and unsupported-capability rules.
3. Update the application execution reference with registry dispatch and the
   generic run-set plan.
4. Reconcile A7 and the status table in the application roadmap. Do not mark
   B3–B6 implemented.
5. Add a concise migration note if any development `.prm` syntax changed.

**Exit gate:** documentation describes the code that passed the final test
matrix, and repository search finds no central list of B1/B2 candidate data
IDs or B3–B6 placeholder branches.

## Final verification matrix

The implementing agent should finish with evidence in this order:

1. non-deal.II application and runner contract tests;
2. parameter-file deal.II contract tests, including unknown-key and
   unused-entry failures;
3. shared scalar-function and B1/B2 deal.II application contracts;
4. the novel B2 forcing-expression no-rebuild test using an ID absent from C++
   source;
5. unchanged B1 seven-combination resolution, artifact paths, manifests, and
   selected-field schema;
6. unchanged B2 four-combination resolution and independent transport boundary
   oracles;
7. a small development smoke run from a copied `.prm`, followed by a second
   run after changing only its scalar expression and provenance; and
8. the authoritative B1 release comparison when the environment and runtime
   budget permit, using the existing contract tolerances and without replacing
   the tracked authoritative profile.

The two smoke invocations must use the same executable. The test or handoff
record should show that no source or build step occurred between them and that
the parameter hash and effective scalar-definition evidence changed.

## Non-goals

This refactor must not:

- implement or claim acceptance of B3, B4, B5, or B6;
- move PDE assembly or optimization algorithms into the runner;
- weaken frozen B1 or B2 mathematical validation;
- turn the public metadata catalog into an owner of backend-specific builders;
- accept arbitrary unknown strings and defer errors to assembly;
- introduce `.prm` includes or make run snapshots depend on external fragments;
- make all geometry, tensor data, or finite-element choices expression-driven;
  or
- regenerate authoritative outputs merely because internal configuration code
  moved.

## Handoff completion condition

The work is complete when B1 remains reproducible under its frozen contract,
B2 accepts a previously unseen deterministic forcing expression and target
definition from a self-contained `.prm` with no rebuild, every declared
run-affecting field has an honest execution owner, and B3–B6 can be added by
registering typed benchmark binders and existing capabilities rather than by
expanding a central B1/B2 parser and dispatch chain.
