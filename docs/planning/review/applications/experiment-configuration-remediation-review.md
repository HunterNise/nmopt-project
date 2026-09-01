# Application-layer `.prm` ownership remediation review

## Status and decision

This document records the corrected remediation scope at
`71286c9be0cf52e4e364167d72be14e5f8ad413f` on 2026-09-01. It follows the
[application configuration ownership review](experiment-configuration-ownership-review.md)
but narrows the implementation handoff to the defect that motivated the work:
experiment values and scalar expressions that fit an existing application
capability must come from deal.II-style `.prm` files rather than compiled
candidate lists or scenario defaults.

The current implementation already contains most of the required machinery:

- one `ParameterHandler`-backed schema registry with benchmark adapters;
- generic matrix expansion and run-set planning;
- typed scalar definitions for zero, constant, and expression data;
- one shared deal.II scalar-function lowerer;
- typed B1 and B2 binders; and
- independent B2 observation-region and target-profile axes.

The remaining scalar-configuration problem is smaller than a new
configuration system. B2 schema declaration still enumerates `constant`,
`parabolic`, `zero`, `constant-one`, and `constant-two` in C++, and the B2
target resolver switches over the two compiled profile names. The repository
also reserves an opaque JSON value inside `Functions/scalar definitions`, even
though the accepted configuration decision assigns deal.II `.prm` syntax to
numerical experiments and JSON only to plotting profiles.

The required correction is:

> Continue using ordinary `.prm` `set` entries and `subsection` blocks. Discover
> the definition IDs already selected by the file, declare those subsection
> paths through `ParameterHandler`, and bind their values into the existing
> typed scalar-function records. Do not add another serialization format,
> parser, syntax layer, or general configuration framework.

B1 and B2 should migrate to the resulting native `.prm` convention after
their current behavior has been reproduced and recorded. Preservation means
equivalent resolved scenarios, run matrices, artifacts, and numerical results;
it does not mean permanently preserving the old file text or content hash.

## Authority and scope

The accepted boundary is the
[parameter and plotting profile decision](../../../decisions/parameter-and-plotting-profiles.md):
numerical experiment families use deal.II-style `.prm` files parsed through a
typed `ParameterHandler` boundary, while JSON is the plotting-profile format.
The detailed public behavior remains in the
[parameter-file reference](../../../reference/parameter-files.md) and the
[application execution reference](../../../reference/application-execution.md).
B1/B2 mathematical and evidence requirements come from the
[Chapter 6 application contract](../../../applications/chapter-6.md), the
[Chapter 6 benchmark contract](../../../benchmarks/chapter-6.md), and the
[Chapter 6 benchmark roadmap](../../chapter-6-benchmark-suite-roadmap.md).

This remediation owns only:

1. native `.prm` declaration and selection of scalar definitions already
   supported by `ScalarFunctionDefinition` and the shared deal.II lowerer;
2. removal of compiled B1/B2 scalar candidate inventories and the unused JSON
   reservation;
3. consumption of B1's scalar desired-state expression through the existing
   scalar capability rather than equality with a compiled polynomial;
4. migration of all tracked B1/B2 parameter families to the selected native
   convention after behavioral baselines exist; and
5. focused verification that a new ID, value, or expression needs only a
   `.prm` edit and no rebuild; and
6. a final, evidence-driven deletion pass over the runner and its application
   tests so the resulting layer contains only current B1/B2 execution needs
   and genuinely shared runner mechanisms.

It does not own a general rewrite of parameter parsing, runner orchestration,
scenario construction, validation policy, or future benchmark execution.

## Correct ownership boundary

The `.prm` file owns experiment data and selections within an implemented
capability:

- scalar definition IDs, kinds, finite values, deterministic expressions, and
  provenance;
- which definition a direct selection or matrix coordinate uses;
- regularisation values and other existing matrix coordinates;
- existing solver, product, mesh, observation, and output selections exposed
  by the registered benchmark adapter; and
- all other values already parsed into typed scenario records.

The framework continues to own:

- semantic ports such as forcing, desired state, target, and future bounds;
- the supported scalar kinds `zero`, `constant`, and `expression`;
- expression validation and the deal.II lowerer;
- benchmark schema shape and the relation between a matrix axis and a semantic
  port;
- actual solver, compiler-product, discretization, and backend capabilities;
  and
- frozen benchmark realizations for geometry or vector/tensor data that do not
  yet have a reusable configurable capability.

Adding another value or deterministic expression for an existing scalar port
must not require a rebuild. Adding a new function family, solver algorithm,
compiler product, discretization, or executable benchmark may still require
C++ implementation and registration. This distinction is the scope guard for
the work.

## Findings that require remediation

### APP-M1 – B2 candidate IDs are compiled into schema declaration

**Priority:** P1.

`append_legacy_b2_definition_entries()` declares target subsections only for
`constant` and `parabolic`, and forcing subsections only for `zero`,
`constant-one`, and `constant-two`. Those strings are experiment coordinates,
not scalar-function capabilities. The `.prm` files already contain the same
IDs in `Matrix/target-profile`, `Matrix/forcing`, or `Functions/forcing`, so the
compiled list is redundant and is the direct reason a previously unseen ID
cannot be parsed without rebuilding.

`b2_target_catalog()` repeats the closure by loading exactly two subsection
names and switching on the selected profile string. The forcing binder is
closer to the desired design because it constructs the subsection path from
the resolved forcing ID, but that path still has to have been declared from the
compiled list.

The same profile closure also appears after parsing. `b2_manufactured_target_id()`
maps only the two current profiles, `b2_case_name()` rejects any other profile,
and `make_b2_scenario()` calls that mapping before `bind_b2_scenario()` can
replace the default catalog with file data. Removing only the schema list would
therefore still leave a new target profile blocked by scenario construction.

**Required outcome:** derive scalar-definition subsection declarations from
the IDs present in the selected `.prm` file. The C++ adapter may know which
axis supplies forcing or target IDs and which subsection prefix contains the
definition fields. It must not enumerate the values selected on those axes.
The production runner must also construct and name a B2 scenario from the
parsed profile and catalog without passing through a manufactured-profile
switch. Frozen manufactured helpers may remain for direct unit fixtures and
public catalog metadata if they are not used as the experiment resolver.

### APP-M2 – B1 and B2 use multiple scalar input conventions

**Priority:** P1 for a clean migration.

B1 and most B2 files use a direct `Functions/forcing` subsection, while the B2
forcing sweep uses named `Functions/forcing definition <id>` subsections. B2
targets use the existing nested `Functions/target definitions/<profile>`
form. The runner therefore carries separate direct and matrix binding paths,
and the fixed schema entries preserve those differences.

The remediation does not need a new grammar to unify the executable behavior.
Use the named forcing-definition subsection convention that already exists in
the B2 forcing sweep for both direct and matrix-selected forcing. A direct file
has one selected ID in `Functions/forcing`; a sweep obtains its IDs from
`Matrix/forcing`. In both cases the selected ID identifies
`Functions/forcing definition <id>`.

Keep the existing B2 target subsection form. Its profile keys come from
`Matrix/target-profile`, and the nested `id` remains the effective scalar
definition ID recorded in evidence. This is ordinary deal.II syntax already
present in tracked files; no `definitions` list, embedded document, escaping
scheme, or alternative parser is needed.

All tracked B1 and B2 files should migrate after equivalence has been
established. Do not retain both direct and named forcing forms indefinitely.

### APP-M3 – the B1 desired-state expression is declared but not consumed

**Priority:** P1 under the same scalar-data rule.

The B1 `.prm` files declare the desired-state ID, kind, expression, and
provenance. `bind_b1_scenario()` nevertheless requires the exact polynomial
kind and expression, while `B1SelectedDataT` constructs the corresponding
polynomial function in C++. The file therefore describes a value that the
runtime does not actually obtain from the file.

The shared scalar lowerer already supports this deterministic scalar
expression. This is not a request for a new function capability.

**Required outcome:** represent the B1 desired state as a
`ScalarFunctionDefinition`, parse the existing `Functions/desired state`
subsection with `kind = expression`, and lower it through the same shared
deal.II scalar utility used by forcing and B2 target data. Remove the compiled
`B1DesiredStateFunction` after tests show equivalent values and ownership.

The desired-state subsection does not need to become a catalog or matrix axis
in this batch. It is a single selected scalar definition whose existing ID,
expression, and provenance are file-owned.

### APP-M4 – the JSON scalar-definition reservation is outside the decision

**Priority:** P1 cleanup because it points implementation in the wrong
direction.

`Functions/scalar definitions` is declared as one opaque string with a JSON
default, its current test verifies only string round-tripping, and the
parameter reference describes a canonical JSON document that no production
binder consumes. This representation is unnecessary once native subsection
IDs are declared from the file. Implementing it would require a second parser
and a second syntax inside `.prm` without removing any mathematical
limitation.

**Required outcome:** remove the entry from the common schema, remove the
opaque round-trip test, and delete the scalar-JSON section from the parameter
reference. Keep JSON only for plotting profiles as recorded by the accepted
decision.

### APP-M5 – preservation currently risks freezing legacy text instead of behavior

**Priority:** P1 for migration order.

B1 is reproduction-verified under its current mathematical contract, and B2
has established executable/development baselines. Those results must protect
the migration. The old `.prm` hashes and exact subsection spellings should not
become permanent tests, because the purpose of the migration is to replace
those spellings with one correct native convention.

**Required outcome:** establish pre-migration behavior with the current files,
then compare the migrated files at the typed scenario, matrix, artifact, and
numerical levels. Expected content-hash changes are migration provenance, not
failures. Remove old syntax support after every tracked family has moved and
the comparisons pass.

Do not create a CMake script that copies authoritative files and mutates exact
lines with `string(REPLACE`. Such a test would bind acceptance to incidental
text and could fail for reasons unrelated to the requested configuration
case. Use the existing parameter contract for focused resolution checks and
ordinary application runs for numerical equivalence.

### APP-M6 – the runner contains compatibility and speculative surface with no current execution owner

**Priority:** P2 cleanup, only after APP-M1–APP-M5 are behaviorally verified.

The application layer currently retains several kinds of code that do not
serve the implemented B1/B2 runner:

- `extension_contracts.hpp` is included only by `runner_contract.cc`. Commit
  `134bb1d` added its box-bound and quadratic-KKT records as future B3–B5
  shapes, and `b3b8297` expanded the hand-built B3–B6 fixtures. Those tests
  also assert that none of the four benchmarks is runnable. A
  production-looking header whose only consumer is a speculative test is not
  future-proofing.
- `ParameterOwnership`, `ownership_for()`, and their three-way accounting are
  the unused accounting layer introduced with schema centralization in
  `dfbe123`: they are read by the schema-characterization test and
  documentation but do not drive declaration, binding, validation, execution,
  or evidence. Maintaining this partial label set would amount to the exact
  ownership ledger explicitly excluded from this remediation.
- `Problem/initial control` duplicates the value consumed from
  `Solver/initial control`, while `Problem/observed material id` duplicates
  `Observation/material id`. The B2 binder accepts the duplication only when
  both copies agree. Once all tracked files are migrated together, those
  aliases add failure modes without adding a choice.
- `Solver/declared minimum step length` is a legacy, provenance-only leaf. No
  runtime policy consumes it; absence already represents the source statement
  that no minimum step was declared.
- the runner capability registries contain quadratic-KKT, PDAS, matrix-free,
  KKT-method, and preconditioner selections that the registered B1/B2
  application paths cannot execute. The KKT and preconditioner registries are
  used only by the speculative extension header and its tests. Accepting a
  name at the first application boundary and rejecting it only later makes an
  unavailable application variant look supported.

These observations do **not** authorize deleting compiler, solver, contract,
or public application capabilities outside `apps/nmopt-runner`. Those layers
have their own tests and purposes. They also do not mean that every test with
"future" in its local variable name is disposable: generic CLI parsing,
generic run-set axes, and artifact-coordinate planning are current runner
mechanisms and should remain covered independently of B3–B6 examples.

**Required outcome:** after migrated B1/B2 equivalence is established, remove
the confirmed compatibility aliases, unused ownership labels, speculative
runner extension records, runner-only unavailable registrations, and the
tests and documentation that exist solely to preserve them. Then perform one
bounded reachability audit for other code with the same properties. Prefer
deletion, direct use of an existing current abstraction, or a single canonical
parameter path. Do not replace removed scaffolding with a new framework.

## Native `.prm` target

The target reuses syntax already present in the repository. A direct B1
forcing should have this shape after migration:

```text
subsection Functions
  set forcing = source-oriented-constant-half
  set desired state = b1-polynomial

  subsection forcing definition source-oriented-constant-half
    set kind = constant
    set value = 0.5
    set provenance = chapter-6.e6.5.1.source-oriented-constant-half-forcing
  end

  subsection desired state
    set kind = expression
    set expression = 10*x0*(1-x0)*x1*(1-x1)
    set provenance = chapter-6.e6.5.1.desired-state
  end
end
```

A B2 forcing sweep keeps its current matrix and native named subsections:

```text
subsection Matrix
  set forcing = zero, constant-one, spatial-candidate
  set target-profile = constant, parabolic
end

subsection Functions
  set forcing = from-matrix
  set desired state = from-matrix

  subsection forcing definition zero
    set kind = zero
    set provenance = development.b2.zero-forcing
  end

  subsection forcing definition spatial-candidate
    set kind = expression
    set expression = 0.4 + sin(pi*x0)*sin(pi*x1)
    set provenance = development.b2.spatial-candidate
  end

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
end
```

The example introduces no syntax construct. It changes one file-owned
candidate and expression using the same `set` and `subsection` forms already
accepted by `ParameterHandler`.

The selection rules are deliberately small:

- if `Functions/forcing` is a concrete ID, declare and select that one named
  forcing subsection;
- if it is `from-matrix`, declare every ID listed in `Matrix/forcing` and
  select the current combination's value;
- declare every target profile listed in `Matrix/target-profile` under the
  existing target-definition subsection; and
- parse the B1 direct desired-state subsection as one scalar definition.

Do not add an independent list of definition IDs. The existing direct selector
and matrix axes already own that information.

## Minimal implementation design

### Extend the existing `ParameterHandler` discovery pass

`read_parameter_file()` already performs a partial `ParameterHandler` parse to
discover `Benchmark/id` and `Benchmark/recipe`, selects the benchmark schema
adapter, declares the full schema, and parses again strictly. Reuse that
pattern.

After selecting the adapter:

1. declare the adapter's static entries, including its matrix axes and direct
   scalar selectors, in a discovery handler;
2. parse the same source with `skip_undefined = true`;
3. read only the forcing selector and the registered scalar-related matrix
   axes needed by that adapter;
4. append the standard scalar fields for the discovered subsection IDs to the
   schema registry; and
5. perform the existing strict parse with the completed registry.

Use `ParameterHandler` for both passes. Do not scan the source with regular
expressions, split arbitrary subsection text, parse a second serialization
format, or duplicate the full `.prm` grammar.

The implementation needs only a small reusable helper that appends
`kind`, `value`, `expression`, and `provenance` entries below a supplied
definition path, plus `id` where the existing target form requires it. The
benchmark adapter may carry the minimal mapping between:

- `Functions/forcing` and optional `Matrix/forcing`;
- the forcing-definition path prefix; and
- `Matrix/target-profile` and the target-definition path prefix.

Do not introduce a universal configuration AST, generic JSON value tree,
plugin registry, reflection layer, or type-erased parameter system. If adding
a reusable descriptor is more code than the two existing adapters need, use a
small helper called by each adapter and generalize only when B3 or B4 actually
arrives.

### Bind definitions by discovered ID

Keep `ScalarFunctionDefinition` and `ScalarFunctionCatalog` as the typed
boundary.

- Build a forcing definition path from the selected forcing ID for both direct
  and matrix cases.
- Build the B2 target catalog by iterating the file's declared
  `target-profile` axis values, not a compiled array.
- Set the catalog's selected scalar ID from the selected profile's parsed
  definition.
- Ensure production B2 scenario construction receives that parsed catalog
  without first calling `b2_manufactured_target_id()` or another fixed-profile
  validator.
- Let B2 case/scenario naming use the resolved profile ID as data. Do not make
  a name-formatting helper the registry of valid scalar profiles.
- Keep the existing validation of kind-specific fields and deterministic
  scalar expressions.
- Carry the effective ID, kind, value or expression, and provenance into the
  same artifact evidence fields already used by B1 and B2.

Do not add B1/B2 names to `ScalarFunctionDefinition`. Semantic role remains in
the binder and evidence adapter; scalar realization remains generic.

### Consume the B1 desired state through the scalar lowerer

Add a desired-state scalar definition to the B1 problem/runtime data at the
narrowest application boundary that needs it. Construct both forcing and
desired-state deal.II functions in `B1SelectedDataT` through the shared scalar
lowerer, retain owner-bearing lifetime, and compare both definitions when
constructing runtime data.

The B1 binder should parse the desired-state ID and subsection instead of
requiring one literal kind and expression. Preserve the existing desired-state
provenance and artifact fields. This change should not alter the B1 problem
recipe, objective, control discretization, mesh, solver, or selected fields.

### Migrate tracked files and remove compatibility in the same roadmap

There are ten tracked B1 families and five tracked B2 families. Migrate all of
them after the new reader and binder can parse both the current and target
forms. Do not migrate only the forcing sweep and leave two production
conventions.

For each file:

- keep benchmark identity, matrix values, selections, exclusions, numerical
  scalar values, expressions, provenance, and all non-scalar settings
  unchanged;
- move direct forcing fields below the existing named forcing-definition form;
- change the B1 desired-state kind to the already supported scalar
  `expression` kind without changing its expression;
- keep the current B2 target-definition structure, changing only declaration
  ownership in C++; and
- expect the parameter content hash to change.

After every tracked file has migrated and passed the equivalence gates, remove
the old direct-forcing schema entries and binding branch, the fixed B2
candidate loops, and any temporary compatibility code. Compatibility is a
short migration seam, not a supported second format.

## End-of-remediation cleanup boundary

The general cleanup starts only after the migrated typed scenarios and
numerical comparisons match their Unit 0 baselines. This ordering keeps a
behavioral change attributable to the scalar remediation rather than mixing it
with deletion work. The cleanup may then deliberately remove unwanted
behavior as well as dead code, but only within the application-layer scope
below.

### What counts as unneeded

An application-layer item is a deletion candidate when repository search and
call-site inspection show one or more of the following:

- it has no production consumer reachable from a registered B1 or B2 run;
- its only consumer is a test that constructs and tests that same speculative
  item;
- it preserves pre-migration scalar syntax or a duplicate parameter alias
  after every tracked file has moved to the canonical path;
- it registers or appears to accept an application selection that no current
  benchmark registration can execute;
- it records an ownership classification, default, or provenance-only leaf
  that no runtime or required evidence path uses;
- it is an unreachable branch, obsolete helper, include, declaration, or test
  fixture left by the remediation; or
- it exposes a supposedly editable leaf in `.prm` but only compares it to a
  compiled literal, while an existing named frozen-profile selector already
  owns the choice.

History identifies likely candidates but is not, by itself, a reason to
delete code. Conversely, pre-existing code receives no presumption of need.
The deciding evidence is the current production call graph, the B1/B2
contracts, artifact/evidence requirements, and the post-migration tests.

### What must remain

Retain code that is used by current B1/B2 execution or that is a genuinely
generic runner primitive already exercised by those registrations. In
particular, keep generic matrix expansion, run-set planning, benchmark
registration and dispatch, artifact coordinates, typed scalar lowering, and
the reduced-method registry if their production consumers remain. Keep tests
of those behaviors, even if a test uses an arbitrary or future-looking ID to
prove that the primitive is not closed over B1/B2.

Do not remove public compiler products, solver policies, semantic contracts,
or backend implementations merely because the runner does not yet expose
them. The cleanup may remove their unused *runner registration* or adapter,
not the underlying framework capability. Do not collapse benchmark-required
provenance, numerical evidence, or frozen scientific identity. A frozen
geometry or vector/tensor realization may remain an internal named profile;
this remediation must not make it configurable by inventing a new parser or
function family.

### Canonicalization and deletion rules

Apply these rules in order:

1. Keep `Solver/initial control` as the single initial-control parameter,
   because that is the path consumed by common solver binding. Remove the B2
   `Problem/initial control` entry, equality check, duplicate values in all B2
   files, and tests of mismatch behavior.
2. Keep `Observation/material id` as the single B2 observation selection and
   bind the typed B2 recipe from it. Remove `Problem/observed material id`, its
   equality check, duplicate file entries, and mismatch test.
3. Remove `Solver/declared minimum step length` from the schema, the B2
   authoritative file, characterization, and parameter reference. Preserve
   the benchmark contract's factual statement that the source declares no
   minimum; absence of the parameter is sufficient.
4. Remove `ParameterOwnership`, `ownership_for()`, the member on
   `ParameterSchemaEntry`, ownership-count assertions, and the documentation
   that promises this unused ledger. Keep schema presence and pattern checks
   that actually protect parsing.
5. Delete `extension_contracts.hpp`, its box-bound and quadratic-KKT tests,
   and the hand-built B3–B6 registration-shape fixtures. Remove
   `PreconditionerSelection`, the KKT-method and preconditioner runner
   registries, and unavailable product/execution registry entries when their
   remaining searches confirm that they have no current runner consumer.
   Retain the generic registry template and current reduced-method/product/
   execution selections only where production binding uses them.
6. For exact-literal B2 profile leaves, remove a leaf only if an existing
   named selector already identifies the frozen realization, the leaf is not
   consumed as numerical data or required provenance, and deletion needs no
   replacement capability. Keep the selector and internal frozen
   realization. If any of those conditions is false, leave the field alone
   and report it rather than expanding the scope.
7. Remove test-only access seams, dead helpers, and includes exposed by the
   preceding deletions. In particular, remove the `#define main` inclusion
   seam only if the scalar work naturally leaves a small existing header
   boundary through which the same production resolver can be tested. Do not
   invent a public resolver architecture solely to eliminate that test seam.

This is a local audit, not a new persisted ownership inventory. Do not add a
spreadsheet, schema annotation, linter, source scan, warning framework, or
"unused parameter" runtime check. Repository search and ordinary compiler
diagnostics are enough to execute the cleanup.

### Streamlining success criterion

The cleanup should be net deletion in the runner and its tests. A small direct
edit that selects a canonical source is acceptable; a replacement abstraction
larger than the code it removes is not. After formatting, compare the scoped
diff and line counts for `apps/nmopt-runner` and the touched application tests.
If a proposed cleanup increases concepts, branches, configuration forms, or
lines without being required for current B1/B2 behavior, revert that part and
record it as out of scope.

## B1 and B2 preservation contract

### B1

Before changing tracked B1 inputs, retain evidence for:

- the seven authoritative method–regularisation combinations and the excluded
  steepest-descent `1e-6` coordinate;
- effective forcing and desired-state definitions;
- mesh, control representation, compiler, solver, and output selections;
- artifact coordinate paths, manifest structure, selected fields, and solver
  traces; and
- the established numerical and qualitative reproduction outcomes and
  tolerances.

After migration, those properties must match. The `.prm` path may remain the
same, while its copied snapshot and content hash are expected to change. Do
not hardcode the old hash in a permanent test.

### B2

Before changing tracked B2 inputs, retain evidence for:

- the authoritative four observation-region/target-profile combinations;
- the three forcing-sweep combinations;
- effective forcing and target definitions and their provenance;
- the independent transport-boundary and observation-region contract tests;
- artifact coordinate paths and manifest structure; and
- the current authoritative/development numerical baselines used by the B2
  replication report.

After migration, these properties must match except for parameter hashes and
snapshots. B2 remains executable but not replication-verified; describe the
comparison as migration equivalence, not as new source reproduction evidence.

## Focused verification strategy

Use the existing test and run boundaries. Do not create a new configuration
test framework.

### Contract tests

Extend `tests/application/parameter_files_dealii_contract.cc` only where the
existing coverage does not already establish the required behavior:

1. Parse a temporary native `.prm` fixture whose forcing matrix contains a
   previously unseen ID and whose matching named subsection contains an
   expression.
2. Resolve the B2 combination and verify the resulting forcing definition and
   selected target definition.
3. Resolve B1's desired state and verify its ID, expression, and provenance.
4. Keep the existing checked-in-family expansion, B1/B2 characterization,
   scalar validation, shared-lowerer, artifact, and boundary tests.
5. Remove tests that inspect the opaque `Functions/scalar definitions` JSON
   value or require old direct-forcing paths.

The temporary fixture must itself be normal deal.II `.prm` text. Do not encode
it as JSON, generate it by patching an authoritative file with CMake string
replacement, or add a separate CTest merely to invoke the runner.

### No-rebuild development check

After building the remediated runner once:

1. copy a small B2 development `.prm` to an ignored experiment file;
2. add or rename one forcing definition using native subsection syntax and
   change its deterministic expression and provenance;
3. run a bounded development case with the existing binary;
4. edit only that `.prm` definition and run again without configuring,
   compiling, or relinking; and
5. verify both manifests complete and their effective scalar evidence and
   parameter hashes reflect the two files.

This is an implementation acceptance run, not a new permanent process-test
harness. Use a valid small B2 mesh so the check exercises scalar ownership
rather than the known empty-observation diagnostic.

### Numerical migration comparison

Run the existing neutral and deal.II pipelines before and after the code and
file migration. For system evidence:

- reproduce the current B1 authoritative run before migration and repeat it
  after migration using the profile and tolerances required by the B1
  contract; and
- capture and repeat the current B2 authoritative or documented development
  comparison used by the B2 replication record.

The release-dealii build and source-sized B1 run require the authorization and
resource handling defined by the repository workflow. If that authorization
is not available during implementation, do not claim the migration complete;
report the remaining system comparison explicitly.

## Implementation handoff

### Unit 0 – record pre-migration behavior

**Outcome:** the current B1/B2 behavior is available as a semantic and
numerical comparison point without freezing old file syntax.

**Actions:**

1. Run the existing neutral and deal.II pipelines at the pre-migration head.
2. Confirm that existing parameter characterization covers all tracked B1/B2
   families and add only missing effective scalar assertions.
3. Run or retain the required current B1 and B2 comparison evidence according
   to the benchmark and build workflows.
4. Record matrices, effective scalar definitions, artifact coordinates,
   manifests, and numerical comparison values in the implementation handoff.

**Do not:** add an old-file content-hash gate, copy the authoritative files
into permanent legacy fixtures, or add a new CMake process test.

**Prospective commit:** none if existing tests are sufficient; otherwise one
focused characterization commit containing only missing semantic assertions.

### Unit 1 – declare native scalar definitions from file-owned IDs

**Outcome:** a previously unseen forcing or target-profile ID in a valid
`.prm` file can be declared and parsed without a C++ candidate-list edit.

**Actions:**

1. Add the small scalar-definition schema-entry helper.
2. Extend the current partial `ParameterHandler` discovery flow to obtain the
   direct forcing selector and relevant matrix-axis IDs.
3. Append forcing and target definition paths from those IDs before the strict
   parse.
4. Change B2 target catalog construction to iterate declared profiles.
5. Remove fixed-profile gates from the production B2 construction and naming
   path while retaining frozen helpers only where they serve fixtures or
   metadata.
6. Route direct and matrix-selected forcing through the same ID-derived path.
7. Add the focused unseen-ID case to the existing parameter contract.

**Do not:** change tracked `.prm` files yet, add JSON, parse raw subsection
text, redesign the complete schema registry, change matrix syntax, or
restructure the run loop, artifact writing, or manifest lifecycle. Make only
the narrow construction change needed to pass the parsed B2 catalog into the
scenario.

**Likely files:** `apps/nmopt-runner/parameter_files.hpp`,
`apps/nmopt-runner/benchmark_binders.hpp`, `apps/nmopt-runner/main.cc`,
`include/nmopt/application/chapter6.hpp`, and
`tests/application/parameter_files_dealii_contract.cc`.

**Prospective commit:** `refactor(applications): declare scalar definitions from prm ids`.

### Unit 2 – consume B1 scalar target data and migrate B1/B2 files

**Outcome:** all tracked B1/B2 files use the native named-definition path and
resolve to equivalent typed scenarios.

**Actions:**

1. Represent and lower B1 desired state through `ScalarFunctionDefinition`.
2. Migrate all direct B1/B2 forcing subsections to the existing named forcing
   definition convention.
3. Change B1 desired-state declarations to the supported expression kind.
4. Keep B2 target subsection contents and profile IDs mathematically
   unchanged.
5. Update existing test expectations from old parameter paths to effective
   typed definitions.
6. Compare every tracked family's resolved matrix and scenario with the Unit 0
   baseline.

**Do not:** change numerical values, expressions, solver policy, meshes,
artifact coordinates, benchmark status, or plotting profiles.

**Likely files:** `include/nmopt/application/chapter6.hpp`,
`include/nmopt/application/dealii/chapter6_b1.hpp`, the runner binders, the ten
B1 and five B2 `.prm` files, and their existing contract tests.

**Prospective commit:** `refactor(applications): migrate b1 b2 scalar prm definitions`.

### Unit 3 – remove obsolete scalar configuration paths and reconcile docs

**Outcome:** one native production path remains and documentation describes
it accurately.

**Actions:**

1. Remove `append_legacy_b2_definition_entries()` and all compiled scalar
   candidate values.
2. Remove the old direct-forcing binding/schema compatibility path.
3. Remove `Functions/scalar definitions` and its opaque JSON test.
4. Remove `B1DesiredStateFunction` after shared-lowerer tests pass.
5. Confirm that `b2_manufactured_target_id()` and fixed-profile checks are not
   reachable from production parameter execution; remove obsolete ones.
6. Update the parameter reference with the native ID-discovery and subsection
   convention.
7. Update the execution reference or application roadmap only where their
   current statements change; do not broaden status claims.
8. Search production C++ for the migrated experiment IDs and confirm they
   remain only in `.prm`, evidence expectations, or benchmark-specific
   scientific records where appropriate.

**Do not:** mix the scalar-path deletion with broader `main.cc`, registration,
run-set, capability-registry, scenario-factory, or parameter-field cleanup.
The bounded general deletion pass belongs to Unit 5, after equivalence is
known.

**Prospective commit:** `refactor(applications): remove legacy scalar prm paths`.

### Unit 4 – prove migrated behavior before general cleanup

**Outcome:** the scalar migration is proven at contract, no-rebuild
development, and system-comparison levels before unrelated application-layer
behavior is removed.

**Actions:**

1. Run the neutral and deal.II pipelines.
2. Perform the same-binary two-file B2 development check.
3. Repeat the B1 and B2 numerical comparisons from Unit 0.
4. Compare manifests and artifacts, allowing only the intentional parameter
   snapshot and hash changes.
5. Report any unavailable release comparison as an explicit remaining gate.

**Prospective commit:** none; verification evidence belongs in the handoff or
the authoritative replication document only if it changes that document's
scientific conclusions.

### Unit 5 – prune the application layer to its implemented purpose

**Outcome:** the runner contains the canonical B1/B2 configuration and
execution paths plus shared primitives they actually use; it contains no
compatibility aliases, speculative B3–B6 contracts, or runner registrations
for application paths it cannot execute.

**Actions:**

1. Start from the candidates in APP-M6. Use `rg`, include/call-site
   inspection, compiler diagnostics, and scoped commit diffs to confirm their
   current consumers. Review history only far enough to understand why a
   candidate exists; do not replay or preserve a past abstraction merely
   because it was committed separately.
2. Canonicalize the two B2 duplicate values to `Solver/initial control` and
   `Observation/material id`, migrate all five tracked B2 files together, and
   remove the duplicate schema entries, binder equality checks, and mismatch
   tests.
3. Remove the legacy declared-minimum-step leaf and the unused
   `ParameterOwnership` classification machinery, while retaining effective
   solver policy and schema declaration/pattern tests.
4. Delete the test-only extension header and self-referential B3–B6 selection
   fixtures. Trim only the runner registry entries and helper types with no
   B1/B2 production consumer; leave framework compiler/solver capabilities
   untouched.
5. Audit remaining exact-literal checks and schema entries using the rules in
   the cleanup boundary. Collapse a duplicate locked leaf to an already
   existing named frozen profile only when this is pure deletion and the
   current scientific identity and evidence remain explicit.
6. Remove dead branches, helpers, includes, comments, tests, and documentation
   exposed by these deletions. Preserve generic runner tests by rewriting
   B3–B6-themed names or fixtures to capability-neutral examples only when the
   underlying generic behavior is used in production.
7. Format the touched C++ files and inspect the scoped diff. Record before/after
   line counts for the runner and touched application tests as a guard that
   this unit is a simplification, not another abstraction batch.

**Do not:** add replacement interfaces for B3–B6, move framework code into the
runner, generalize frozen geometry/vector data, introduce ownership metadata
or unused-parameter checks, or refactor a live subsystem merely because its
style could be improved. Do not delete generic behavior with a current B1/B2
consumer.

**Prospective commits:** keep the review surface small and causal. Prefer one
commit for canonical parameter cleanup and one for removal of speculative/dead
runner code if the combined diff is larger than a straightforward review:

- `refactor(applications): remove duplicate b2 parameter aliases`
- `refactor(applications): prune unused runner scaffolding`

### Unit 6 – verify the streamlined result

**Outcome:** cleanup has reduced the application surface without changing the
verified B1/B2 experiment behavior or weakening the no-rebuild loop.

**Actions:**

1. Run the neutral and deal.II pipelines again after Unit 5.
2. Re-run the same-binary two-file B2 development check without rebuilding
   between parameter edits.
3. Compare all B1/B2 typed configurations, matrices, manifests, artifact
   coordinates, and numerical results with the Unit 4 result. Changes from
   removal of duplicate/provenance-only parameter leaves may alter snapshots
   and hashes only; document those expected differences.
4. Search `apps/nmopt-runner` and application tests for the removed aliases,
   legacy scalar paths, ownership classifications, speculative extension
   types, and unavailable runner registrations. There should be no production
   or self-testing remnants.
5. Report any candidate intentionally retained, its current production or
   evidence consumer, and why deleting it would exceed this remediation.

**Prospective commit:** none unless an authoritative documentation statement
must be corrected to describe the final application surface.

## Explicit non-goals

The implementing agent must not add or redesign any of the following for this
remediation:

- JSON or another embedded serialization format for numerical `.prm` data;
- a custom `.prm` parser, syntax preprocessor, or regular-expression scanner;
- a new list of scalar definition IDs separate from existing selectors and
  matrix axes;
- include-directive rejection or other convention-enforcement checks;
- numeric-coordinate canonicalization;
- a generalized diagnostic-context framework;
- a new CMake/process-test harness or string-rewritten authoritative fixtures;
- an exact ownership ledger for every schema path;
- generic geometry, region-expression, vector-function, or tensor-function
  configuration;
- preflight/output-lifecycle restructuring;
- a rewrite of default scenario construction;
- relocation of B1/B2 execution registrations or a plugin system;
- implementation or advertisement of B3, B4, B5, or B6; or
- permanent parsing support for the pre-migration B1/B2 scalar syntax;
- deletion or redesign of framework compiler, solver, semantic, or backend
  capabilities outside the runner merely because B1/B2 do not expose them;
  or
- a repository-wide style, naming, header-layout, or architecture cleanup.

Existing conventions should continue to govern self-contained parameter
files, valid IDs, output placement, and provenance. Do not add code to enforce
a convention unless the scalar ownership change itself makes that enforcement
necessary.

## B3–B6 future-proofing, kept proportional

No future benchmark fixture or fake execution path is needed. The remediation
is sufficiently future-proof if the small native definition-declaration helper
can later be reused as follows:

- B3 registers target and bound scalar ports when its actual adapter is
  implemented;
- B4 uses the same scalar kind/value/expression fields for lower and upper
  bounds;
- B5 can instantiate the existing typed registry pattern for the compiled
  products, KKT methods, and preconditioners its real application adapter can
  execute; unused B5-specific runner registries need not exist before then;
  and
- B6 changes its reaction or other already exposed scalar data in `.prm` while
  reusing B5's implemented application path.

Do not implement those registrations now and do not retain speculative data
contracts as placeholders. The proof required in this batch is that the
scalar-definition helper, generic run-set plan, and registration pattern
contain no B1/B2 candidate values or structural assumption that scalar IDs and
matrix axes are limited to the current files.

## Completion condition

The remediation is complete when:

1. all tracked B1/B2 families use ordinary deal.II `.prm` subsections for
   their scalar definitions;
2. the schema and binders derive forcing and B2 target candidates from file
   selectors and matrix values rather than compiled ID lists;
3. B1 desired-state data is parsed and lowered from its `.prm` definition;
4. the unused scalar JSON reservation and old scalar syntax paths are gone;
5. a new B2 scalar ID and deterministic expression run twice through one
   already-built binary with only `.prm` edits;
6. migrated B1/B2 files reproduce their pre-migration typed configuration,
   artifact structure, and numerical behavior, aside from expected parameter
   hashes and snapshots;
7. duplicate B2 parameter aliases, the legacy declared-minimum-step leaf, the
   unused ownership ledger, speculative extension contracts, and unavailable
   runner registrations are gone;
8. a bounded final audit finds no other dead compatibility branch,
   self-testing production header, or exact-literal pseudo-choice that can be
   deleted within the stated rules; and
9. no unrelated parser, validation, runner, or future-benchmark framework has
   been introduced, and the scoped cleanup is a net simplification.
