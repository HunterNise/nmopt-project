# Repository-structure readability refactor

**Working branch:** `codex/repo-structure-refactor`  
**Starting point:** documentation structure commit `997367f37b14da0010e0e474cf1e75690c4d79fb`  
**Status:** active planning and implementation guide  
**Assessment:** `docs/history/reviews/human-readability-audit/09-repository-structure-follow-up.md`

## 1. Purpose and how to use this roadmap

This roadmap turns the non-documentation findings of the human-readability
audit into bounded, reviewable work.

The objective is:

> make the repository's physical source, test, application, tooling, and build
> structure reflect the architecture that already exists, while preserving
> behavior unless a separately identified correctness defect requires a
> deliberate correction.

This is **not** an authoritative implementation sequence. The numbered units
are a suggested decomposition based on the repository state inspected at
`997367f`. They exist to orient an implementation agent, not to replace fresh
inspection.

Before implementing any unit, the agent must:

1. inspect the current branch head and working tree;
2. re-open the files named by the unit rather than relying on this roadmap's
   remembered line counts or exact contents;
3. search current uses before moving, renaming, or splitting a public file;
4. check whether earlier units or user work have changed the proposed boundary;
5. revise the local plan if current evidence no longer matches the suggestion;
6. stop and ask when the unit crosses an AMBER decision gate.

A later unit may become unnecessary after an earlier cleanup. Conversely, a
unit may need to be split if current code reveals a larger responsibility
boundary than the audit saw.

This is not a compiler rewrite, API simplification project, benchmark project,
or continuation of the detailed documentation rewrite.

The detailed documentation/reference rewrite is paused while this branch
stabilizes source paths and ownership boundaries. After this roadmap reaches
an agreed stopping point, its result should be integrated back into
`codex/docs-refactor` and documentation rewriting should resume against the
settled source tree.

## 2. Repository mental model to preserve

The restructuring must not blur the architecture that motivated the audit.

Two producer paths reach common numerical/formulation contracts:

```text
semantic/compiler path
ProblemSpec
  -> validation/resolution
  -> compiler/lowering
  -> numerical realization
  -> compiled formulation product
                         \
                          +-> common formulation / solver contracts
                         /
native/application path
existing application
  -> thin numerical adapters
  -> layouts, actions, solves, metrics, constraints
  -> formulation product
```

They converge at shared numerical/formulation contracts, **not** at
`ProblemSpec` or `CompiledProblemT`.

Important distinctions to preserve:

- primal versus covector;
- derivative versus gradient;
- action/application versus inverse/solve;
- residual/JVP/VJP;
- state solve versus adjoint solve;
- solve result versus solve evidence/report;
- compiler-owned products versus application-owned/native products;
- reusable numerical core versus optional application/reproduction shell.

There is no universal formulation interface. Reduced DTO, supplied OTD,
quadratic KKT, and complementarity/PDAS remain distinct products.

The repository structure may become clearer; the mathematics and ownership
model must not be flattened to achieve that.

## 3. Fixed decisions and non-goals

### 3.1 `apps/` remains the application/showcase area

`apps/` contains concrete final consumers of the reusable nmopt core.

The current examples demonstrate different consumption modes:

```text
semantic/compiler-driven
configured recipe/scenario/parameter-driven
external/native application integration
```

Do not create a new top-level `integrations/` directory.

Do not introduce new `apps/` subdirectories merely to make the tree look
symmetric. If a stronger taxonomy appears useful after the other cleanup,
present alternatives and ask first.

### 3.2 Keep `.agents/explanation-essay.md`

It is an intentionally unrouted reference. Do not delete, merge, or rewrite it
as duplicate cleanup.

### 3.3 Step-4 physical restructuring is late work

Treat `apps/external-dealii/step-4/` as a coherent snapshot. Do not move its
source/evaluation directories or rewrite its source-accounting layout during
early work.

A later unit may improve the binding seam without modifying the preserved
upstream/baseline/adapted source snapshots.

### 3.4 Preserve behavior by default

Do not combine a file move with algorithmic cleanup, numerical-policy changes,
new defaults, changed tolerances, changed experiment choices, or altered
evidence semantics.

When a move exposes a possible behavioral bug, record it and handle it in a
separate correction unit.

### 3.5 Do not resume the broad docs rewrite here

This branch may:

- maintain this roadmap/audit;
- add small local README/routing material required by the new structure;
- update authoritative documentation when a changed interface would otherwise
  make it false.

It should not perform the documentation-wide link, authority, reference, or
manual rewrite already owned by `codex/docs-refactor`.

## 4. Reading and instruction routing

The repository's agent instructions remain authoritative. This section only
orients the next agent toward the minimum likely reading set.

### 4.1 Before any implementation work

Read:

```text
.agents/README.md
.agents/workflow.md
.agents/git.md
```

Then read the task-specific instruction:

```text
.agents/code.md            C++, tests, interfaces, source structure
.agents/build.md           CMake, targets, test registration, validation
.agents/run.md             runner, parameters, runs, post-processing
.agents/documentation.md   Markdown or documentation changes
```

Do not read `.agents/explanation-essay.md` as normative routing; it is retained
as reference only.

### 4.2 Architecture/reference documents that are likely useful

Read only the ones needed by the current unit.

Core numerical/interface changes:

```text
docs/design/interface-specification.md
docs/design/architecture.md
docs/design/composition-boundaries.md
docs/design/pde-solver-boundary.md
```

Compiler/semantic implementation work:

```text
docs/internals/compiler/semantic-compiler.md
docs/internals/system-blueprint.md
```

Application/runner work:

```text
docs/reference/application-api.md
docs/reference/application-execution.md
docs/reference/parameter-files.md
```

External/native application integration:

```text
docs/reference/external-dealii-solver-integration.md
apps/external-dealii/step-4/README.md
apps/external-dealii/step-4/external-integration-overview.md
```

Parameter/post-processing policy:

```text
docs/design/decisions/parameter-and-plotting-profiles.md
docs/reference/parameter-files.md
docs/reference/application-execution.md
```

The follow-up audit is evidence and orientation, not a replacement for these
current authorities.

## 5. Agent authority levels

Every work unit below has one of these authority levels.

### GREEN — proceed without another architecture decision

An agent may implement the unit after inspecting current branch state and the
named files, provided it stays inside the stated compatibility boundary.

Examples:

- moving implementation detail behind an unchanged public facade;
- splitting a giant test file while preserving scenario behavior;
- extracting production runner logic from `main.cc` without changing CLI or
  output behavior;
- correcting implementation so it obeys an already accepted design document.

If the implementation cannot stay inside those boundaries, GREEN immediately
becomes AMBER.

### AMBER — stop and ask before implementation

Use AMBER whenever more than one materially different repository/API design is
valid.

The agent must present:

1. the concrete problem;
2. current evidence from the code;
3. two or more viable alternatives;
4. affected public paths/types/namespaces;
5. migration/compatibility cost;
6. focused validation implications;
7. its recommendation.

Do not choose silently.

Typical AMBER changes include public header renames, new top-level categories,
application taxonomy, reference-problem placement, namespace changes, shared
Python package boundaries, and major compiler decomposition.

### RED — out of scope unless explicitly reopened

Do not perform these as opportunistic cleanup:

- redesign numerical contracts;
- simplify primal/covector, derivative/gradient, application/solve, or
  state/adjoint distinctions;
- change benchmark mathematics, solver tolerances, or reproduction choices;
- introduce `.cc` implementation solely because a header is large;
- restructure preserved Step-4 source snapshots;
- remove compatibility aliases/headers without a dedicated decision;
- configure or build `release-dealii` without explicit permission.

## 6. Compatibility contract

Unless a unit explicitly says otherwise, preserve all of the following.

### 6.1 C++ and numerical API

Preserve:

- namespaces;
- public type/function names;
- template parameters and semantics;
- include paths;
- ownership/lifetime behavior;
- formulation products and supported capability closure;
- diagnostics except explicitly approved wording-only cleanup.

When physically moving a public header, keep the old include path as a
forwarding compatibility header unless the user explicitly approves a
breaking migration.

### 6.2 Applications

Preserve:

- CLI flags and accepted values;
- return codes;
- default behavior;
- parameter precedence;
- run-directory allocation;
- artifact naming and schemas;
- benchmark IDs;
- run/provenance semantics.

### 6.3 Tests

The externally meaningful compatibility surface is the logical test inventory:

- scenario names;
- CTest names;
- labels;
- timeouts;
- assertions and expected behavior.

Test executable target names may need to change when one giant translation
unit becomes several executables. Do not treat that as automatically harmless:
search for direct target users first. Preserve user-facing application/library
target names unconditionally unless a separate decision says otherwise.

Physical test splitting must not silently become a coverage rewrite.

### 6.4 Build

Preserve:

- user-facing application/library targets;
- build options;
- preset names;
- deal.II optionality;
- scenario discovery behavior;
- normal build outputs.

### 6.5 Python evidence tooling

Preserve persisted source data. Derived output may intentionally change only in
a correction unit that explicitly changes rendering behavior.

Malformed/non-finite persisted evidence must never be silently reinterpreted
merely to keep a report running.

## 7. Validation policy

The repository is large. Validation should be proportional and cumulative,
not repeated after every small edit.

### 7.1 General rule

A **coherent work unit**, not an individual file edit, is the validation
boundary.

Several tightly related edits may be made before running focused verification.
Do not run a complete pipeline after every comment, include change, or file
move.

A recent successful full-profile run may be reused as the baseline for
subsequent focused work if:

- no other actor changed the relevant branch;
- build configuration did not change;
- generated/build state is still trustworthy;
- the intervening edits did not touch that profile's behavior.

Before each unit, inspect repository status and branch head. Never overwrite
unrelated user work.

### 7.2 Validation levels

#### V0 — documentation/text/path-only

Use when no executable source/build behavior changed.

Validate with:

- path/search checks for renamed references;
- Markdown structure/link inspection for files directly touched;
- `git diff --check`;
- diff review.

Do **not** compile merely because a Markdown/README/comment-only unit changed.

#### V1 — Python/tool-only

For Python changes:

1. run `python3 -m py_compile` on changed Python modules;
2. run the focused contract script(s) that exercise the changed behavior;
3. compare generated metadata/file inventory when relevant;
4. if CTest registration changed, run the focused registered tests after
   configuring once.

Do not rebuild C++ unless the Python change also modifies CMake registration or
the runner contract.

#### V2 — backend-neutral C++ focused

For neutral headers/tests:

1. configure only if target/test registration changed;
2. build only affected test/application targets first;
3. run affected CTest names/labels/regexes;
4. inspect discovered scenario inventory if translation units were split.

A full:

```bash
./build.sh pipeline debug-neutral
```

is a **stage gate**, not a per-edit requirement. Run it when:

- leaving a cluster of neutral refactors;
- changing a public neutral facade;
- before an AMBER public-API decision that relies on a clean baseline;
- before final handoff.

#### V3 — deal.II C++ focused

For deal.II compiler/lowering/application/test work:

1. configure only when target/test registration changed;
2. build the smallest affected target set;
3. run focused CTest names/labels;
4. inspect logical scenario inventory after test restructuring.

A full:

```bash
./build.sh pipeline debug-dealii
```

is a stage gate. Run it when:

- leaving a cluster of deal.II-sensitive refactors;
- CMake registration has been substantially reorganized;
- before Step-4 work;
- before final handoff.

The helper's configured job ceiling is authoritative. If direct CMake is
unavoidable, manual deal.II builds use `--parallel 1`.

#### V4 — integration/final gate

Before the source-structure branch is handed back to the docs branch, run:

```bash
./build.sh pipeline debug-neutral
./build.sh pipeline debug-dealii
```

Run `sanitize-neutral` only if the implemented work changed ownership,
lifetime, allocation-sensitive solver logic, or another area where sanitizer
evidence is materially relevant.

Do not build `release-dealii` without explicit permission.

### 7.3 Baseline strategy

Before the first source change, establish or confirm a trustworthy
`debug-neutral` baseline. Do not rerun it merely because R0 changed Markdown.

Establish `debug-dealii` before the first deal.II-sensitive cluster, not before
unrelated neutral or Python-only work.

If `build.local.conf` is required and absent, stop and ask the user to run
`./build.sh init-config`; do not create it implicitly.

## 8. Suggested work units

The order below is a **suggested dependency-aware decomposition**, not a
mandatory sequence.

An agent may reorder units when current evidence makes another order safer or
more economical. It should explain the reason in its working plan.

A unit may be skipped entirely if the user decides its benefit does not justify
the churn.

---

## R0 — Record assessment and active roadmap

**Authority:** GREEN  
**Intervention:** documentation only  
**Suggested commit:** `docs(audit): record repository structure follow-up`

### Orientation

This is the bridge between the completed documentation-only audit and the
separate source-structure branch.

### Read first

```text
.agents/documentation.md
docs/history/reviews/human-readability-audit/00-audit-index.md
docs/history/reviews/human-readability-audit/05-source-readability.md
docs/history/reviews/human-readability-audit/08-recommendations.md
```

### Likely files

```text
docs/history/reviews/human-readability-audit/00-audit-index.md
docs/history/reviews/human-readability-audit/09-repository-structure-follow-up.md
docs/planning/repository-structure-refactor.md
```

### Suggested change

- preserve the original audit's documentation-only conclusion as historical;
- add the later dedicated source-structure assessment;
- add this mutable roadmap under `docs/planning/`.

Do not rewrite old audit files to make them appear to have reached the later
conclusion.

### Alternatives

None material. If the docs tree changed again before applying R0, re-check the
appropriate history/planning destinations.

### Validation

V0 only.

---

## R1 — Cheap path-neutral readability cleanup

**Authority:** GREEN  
**Intervention:** light edits  
**Suggested commit:** `refactor(readability): remove stale implementation-era wording`

### Orientation

The audit found current reusable code and diagnostics that still use roadmap
milestone names (`v0`, `P5.x`, `P6.x`, `C5.x`) where a caller needs the
mathematical capability name instead.

This is wording cleanup, not a source-history purge. Historical labels remain
valid in provenance, audit material, benchmark identities, and genuinely
historical tests.

### Read first

```text
docs/history/reviews/human-readability-audit/09-repository-structure-follow-up.md
docs/design/interface-specification.md
```

Then inspect only the files where current searches find the labels.

### Likely affected files observed at the audit baseline

Representative locations included:

```text
include/nmopt/contract/linalg.hpp
include/nmopt/contract/reduced_dto.hpp
include/nmopt/contract/supplied_otd.hpp
include/nmopt/compiler/v1/dealii_compiler.hpp
include/nmopt/compiler/v1/dealii_types.hpp
tests/contract/reduced_dto_contract.cc
tests/semantic/semantic_v1_contract.cc
tests/dealii/dealii_diffusion_contract.cc
tests/dealii/dealii_trace_hhalf_metric_contract.cc
```

Do not assume every occurrence should change. Search first and classify each
use as current runtime language, test traceability, or historical provenance.

Also add:

```text
parameters/README.md
```

to explain:

- tracked `.prm` files are executable inputs;
- `authoritative.prm` means the accepted repository reproduction input, not
  that every value is directly stated by the source publication;
- `development/` is explicitly hypothesis/investigation material;
- plotting JSON is versioned presentation policy;
- generated/copied configuration belongs under `runs/`.

### Suggested implementation

For runtime/compiler diagnostics, prefer descriptions such as:

```text
continuous Neumann control
coefficient identification
weighted boundary trace
supplied OTD product
reduced DTO contract
```

rather than roadmap IDs.

### Alternatives / decision gate

If a message is simultaneously a stable user-facing diagnostic and an
important historical identifier, do not invent a policy. Present the concrete
case and ask.

### Validation

- README/comment-only edits: V0;
- diagnostics that are asserted by tests: build/run only the affected test
  targets/scenarios;
- no full pipeline unless R1 is being used as the end of a larger neutral
  cleanup stage.

---

## R2 — Split neutral test responsibilities

**Authority:** GREEN for physical splitting; AMBER for changing logical CTest identities  
**Intervention:** static split  
**Suggested commit:** `test(structure): separate neutral contract and solver coverage`

### Orientation

The test harness is already strong: each C++ executable exposes named
`Scenario` records and scenario discovery generates CTest entries.

That means source-file identity and CTest identity are decoupled. The main
problem is that some translation units now contain several architectural
responsibilities.

### Read first

```text
.agents/code.md
.agents/build.md
tests/support/scenario_dispatch.hpp
cmake/ScenarioDiscovery.cmake
cmake/WriteDiscoveredScenarios.cmake
```

Then inspect:

```text
tests/contract/reduced_dto_contract.cc
tests/semantic/semantic_v1_contract.cc
CMakeLists.txt
```

### Current evidence

`tests/contract/reduced_dto_contract.cc` was about 2,944 lines but contained
only eight top-level scenarios. Its `v0_contract` scenario exercises much more
than the DTO contract:

- Armijo/fixed/exact/Wolfe variants;
- actual displacement handling;
- nonlinear CG;
- L-BFGS/full BFGS;
- Newton/PCG;
- reduced search solver;
- stopping/evidence/work accounting.

`tests/semantic/semantic_v1_contract.cc` combines semantic validation and
reference-graph stability with semantic resolution and a
`DealiiScalarLoweringPlanner` scenario.

### Suggested shape

Use responsibility boundaries similar to:

```text
tests/
  contract/
    reduced_dto*.cc

  solvers/
    reduced_common*.cc
    reduced_directions*.cc
    reduced_line_search*.cc
    reduced_search*.cc
    reduced_trust_region*.cc

  semantic/
    validation*.cc
    reference_problems*.cc
    resolution*.cc

  compiler/
    scalar_lowering_plan*.cc
```

These names are suggestions. Prefer fewer coherent files over one tiny file
per algorithm.

### Compatibility target

Preserve existing logical CTest names, labels, timeouts, and assertions first.

The old `nmopt.contract.v0` catch-all is itself a readability problem, but
changing its logical scenario identity is a separate decision.

### Alternatives

1. **Preferred:** split both directory responsibility and translation units.
2. **Lower churn:** split giant files but leave the current top-level test
   directories temporarily.
3. **Not recommended:** rename scenarios and reorganize physical files in one
   commit.

If direct users depend on old test target names, present that evidence before
choosing new executable targets.

### Validation

- capture pre-change discovered CTest names for the affected executable(s);
- reconfigure once because target registration changes;
- build only affected neutral test targets;
- compare post-change logical scenario inventory;
- run all moved scenarios;
- run a full `debug-neutral` pipeline only if this closes the current neutral
  restructuring cluster.

---

## R3 — Split the deal.II compiler capability test ledger

**Authority:** GREEN for responsibility-based physical splitting  
**Intervention:** static split  
**Suggested commit:** `test(dealii): split compiler capability coverage`

### Orientation

`tests/dealii/dealii_diffusion_contract.cc` is not really one diffusion test
anymore. At the audit baseline it was nearly 8,000 lines with 34 scenarios and
functioned as the executable deal.II compiler capability ledger.

This unit should improve discoverability **without changing compiler code**.

### Read first

```text
.agents/code.md
.agents/build.md
docs/design/interface-specification.md
docs/internals/compiler/semantic-compiler.md
tests/support/scenario_dispatch.hpp
tests/support/manifest_contracts.hpp
tests/dealii/dealii_diffusion_contract.cc
CMakeLists.txt
```

### Current scenario families observed

The file covered areas including:

```text
canonical volume control
compiled DTO/KKT/PDAS
fixed Dirichlet
Dirichlet control
partial/L2-transposition Dirichlet control
subdomain and point observations
normal flux / H1 state observation
Neumann and convection/subdomain control
weighted boundary trace
H1 / H-1 metrics
quadratic forms
continuous controls
volume-observation assembly
coefficient identification
pure Neumann / general scalar Robin
projection compatibility
compiler diagnostics
compiler sessions and lifetime
serial SPD reporting
```

Re-enumerate the current scenario list before splitting.

### Suggested shape

A likely home is:

```text
tests/compiler/dealii/
```

with capability-grouped files such as:

```text
canonical_volume*.cc
compiled_formulations*.cc
fixed_dirichlet*.cc
dirichlet_control*.cc
neumann_control*.cc
continuous_control*.cc
observations*.cc
metrics*.cc
coefficient_identification*.cc
general_scalar*.cc
diagnostics*.cc
session*.cc
```

Do not treat those filenames as mandatory. Group by shared setup and
responsibility.

### Alternatives

1. `tests/compiler/dealii/` — best architectural signal.
2. Keep `tests/dealii/` but split the giant file — less path churn.
3. A hybrid where low-level backend services stay in `tests/dealii/` and
   compiler-lowering scenarios move to `tests/compiler/dealii/`.

If the split reveals genuine low-level backend tests mixed into the ledger,
option 3 may be best. Reassess rather than following the suggested tree
blindly.

### Do not absorb

Keep these out of R3:

```text
tests/dealii/external_step4_native_contract.cc
tests/application/external_step4_*.cc
```

Step-4 classification is intentionally late.

### Validation

- record pre-change logical CTest inventory;
- reconfigure `debug-dealii` once;
- build affected test targets only;
- compare CTest names/labels/timeouts;
- run all moved scenarios;
- full `debug-dealii` only when leaving the current deal.II test-restructure
  cluster.

---

## R4 — Extract reusable runner behavior from `main.cc`

**Authority:** GREEN for minimal extraction; AMBER for broad redesign  
**Intervention:** substantial refactor with behavior preservation  
**Suggested commit:** `refactor(runner): extract reusable execution seams`

### Orientation

The runner is a genuine application and should remain under
`apps/nmopt-runner/`.

The current problem is not that every helper file is too large. It is that
tests need production functions that still live inside `main.cc`.

### Read first

```text
.agents/code.md
.agents/build.md
.agents/run.md
docs/reference/application-api.md
docs/reference/application-execution.md
docs/reference/parameter-files.md
```

Then inspect:

```text
apps/nmopt-runner/main.cc
apps/nmopt-runner/runner.hpp
apps/nmopt-runner/parameter_files.hpp
apps/nmopt-runner/run_set_plan.hpp
apps/nmopt-runner/benchmark_registry.hpp
apps/nmopt-runner/capability_registry.hpp
apps/nmopt-runner/benchmark_binders.hpp
apps/nmopt-runner/parameter_binding.hpp
tests/application/parameter_files_dealii_contract.cc
tests/application/runner_contract.cc
```

### Current evidence

The characterization test currently uses the entrypoint as an include:

```cpp
#define main nmopt_runner_characterization_entrypoint
#include "../../apps/nmopt-runner/main.cc"
#undef main
```

The test needs reusable functions such as execution-registration lookup,
artifact-coordinate helpers, and parameter/artifact enrichment.

`main.cc` also owns B1/B2 execution and configuration snapshot/run-preparation
logic.

### Required outcome

Tests must no longer include `main.cc`.

### Suggested implementation

Extract the **smallest coherent production seam** required by tests.

Likely private application modules include some subset of:

```text
execution_registry.hpp
run_preparation.hpp
artifact_fields.hpp
```

The names are suggestions.

Keep `main.cc` responsible for:

```text
parse CLI
handle help/list
prepare/delegate
report failure/success
return process code
```

Do not automatically split `runner.hpp` and `parameter_files.hpp` in the same
unit.

### Alternatives / decision gate

If current code shows two plausible ownership models, for example:

```text
artifact metadata owned by generic runner support
vs
artifact metadata owned by B1/B2 execution adapters
```

stop and present both before extracting broadly.

Likewise, if removing the `main.cc` include requires changing public runner
types or CLI behavior, escalate to AMBER.

### Validation

Focused first:

- build `nmopt_runner`;
- build parameter/runner test target(s);
- run runner/parameter CTest scenarios;
- compare representative parameter resolution and artifact metadata results.

Do not run a benchmark reproduction.

A full `debug-dealii` pipeline is only needed if this closes the current
deal.II-sensitive stage or if the extraction touches shared application/deal.II
headers broadly.

---

## R5 — Correct the plotting-profile implementation

**Authority:** GREEN because the governing design is already accepted  
**Intervention:** explicit behavior correction  
**Suggested commit:** `fix(postprocess): honor plotting profile policy`

### Orientation

This is one of the few units intended to change observable behavior.

The accepted parameter/plotting decision says presentation policy belongs in
the versioned profile and the generic renderer must not contain B1/B2 branches.

The checked-in JSON files already declare policy that the implementation
currently hard-codes or ignores.

### Read first

```text
.agents/run.md
docs/design/decisions/parameter-and-plotting-profiles.md
docs/reference/parameter-files.md
docs/reference/application-execution.md
```

Then inspect:

```text
tools/postprocess.py
tools/nmopt_postprocess/pipeline.py
tools/nmopt_postprocess/render.py
tools/nmopt_postprocess/chapter6.py
tools/nmopt_postprocess/parameters.py
parameters/plotting/chapter-6-b1.json
parameters/plotting/chapter-6-b2.json
tests/tools/postprocess_configuration_contract.py
tests/tools/postprocess_comparison_plan_contract.py
tests/tools/postprocess_history_contract.py
tests/tools/postprocess_boundary_field_contract.py
```

### Current defects observed

Generic `pipeline.py` currently has B2-specific metadata fallback such as:

```text
b2.<axis>
```

Generic `render.py` currently hard-codes Chapter-6 presentation choices such
as:

```text
turbo colormap
finite-extrema normalization
shared comparison normalization
Gouraud interpolation
no volume mesh overlay
five colorbar ticks
180 DPI
x/y axis labels
```

The JSON profiles already declare these under `defaults`, but the loader does
not currently carry all of them into rendering.

### Suggested implementation

Introduce a resolved presentation-policy object as part of the profile model,
or equivalent fields on `PostprocessProfile`.

Flow:

```text
JSON profile
  -> schema validation / parsing
  -> resolved PostprocessProfile presentation policy
  -> generic pipeline
  -> generic renderer
```

Move benchmark-specific metadata-key fallback into the Chapter-6 profile
adapter.

Implement the **currently declared supported policy**, not a generic plotting
framework.

### Alternatives

The precise Python dataclass/file organization is not important.

Valid shapes include:

1. presentation fields directly on `PostprocessProfile`;
2. a nested `RenderPolicy`/`PresentationPolicy` dataclass;
3. a small `model.py` extraction if it materially clarifies imports.

Do not create `model.py` merely because the audit suggested it. Choose based on
the current dependency graph.

### Tests

Add focused evidence that:

- a non-default profile colormap/DPI/axis-label or other supported policy is
  actually resolved and consumed;
- generic axis lookup no longer requires a B2 key convention;
- existing B1/B2 profiles retain their intended output policy.

Prefer testing resolved policy and renderer calls/file inventory over brittle
pixel-perfect raster comparisons.

### Validation

V1 only unless CMake registration is changed.

Run all existing post-processing contracts plus new policy coverage.

---

## R6 — Give Python tool contracts a supported test path

**Authority:** GREEN for straightforward registration; AMBER for a new permanent profile/option policy  
**Intervention:** build/test integration  
**Suggested commit:** `test(tools): register postprocessing contracts`

### Orientation

The audit found four meaningful Python contracts that exist as executable test
scripts but are not part of the supported CTest path.

The root build already registers the external Step-4 forward-comparator Python
contract, so there is precedent for Python tests.

### Read first

```text
.agents/build.md
tools/README.md
DEPENDENCIES.md
CMakeLists.txt
tests/tools/*.py
```

### Existing tests to account for

```text
tests/tools/external_dealii_forward_contract.py
tests/tools/postprocess_configuration_contract.py
tests/tools/postprocess_comparison_plan_contract.py
tests/tools/postprocess_history_contract.py
tests/tools/postprocess_boundary_field_contract.py
```

Add a focused contract for:

```text
tools/chapter6_report.py
```

especially run-manifest safety/status handling and evidence parsing.

### Suggested policy

Separate:

- standard-library-only Python contracts, which can run whenever Python is
  available;
- plotting contracts that require NumPy/Matplotlib/meshio.

Do not make ordinary backend-neutral C++ verification unexpectedly fail because
optional plotting libraries are missing.

### Alternatives / decision gate

Plausible approaches:

1. dependency-detect and register available plotting tests;
2. one explicit `NMOPT_ENABLE_TOOL_TESTS` option;
3. a documented standalone Python test command outside CTest.

The audit prefers a supported CTest path, but the exact option/detection policy
affects normal developer behavior. If it is not obvious from current build
conventions, present alternatives and ask.

Do **not** add a new permanent CMake preset without approval.

### Validation

- configure once after registration changes;
- run only tool/Python tests first;
- deliberately verify behavior when optional plotting packages are absent or
  document why that environment cannot be simulated;
- no deal.II compilation unless the chosen registration is incorrectly tied to
  the deal.II block, in which case reconsider the placement.

---

## R7 — Unify persisted artifact parsing

**Authority:** AMBER  
**Intervention:** substantial refactor / evidence semantics  
**Suggested commit after decision:** `refactor(tools): share artifact parsing`

### Orientation

This is not just code deduplication. Two consumers currently interpret
persisted evidence differently.

### Read first

```text
.agents/run.md
docs/reference/application-execution.md
tools/nmopt_postprocess/records.py
tools/chapter6_report.py
tests/tools/postprocess_history_contract.py
```

Also inspect the runner-side writer before finalizing parsing semantics:

```text
include/nmopt/application/artifact_writer.hpp
apps/nmopt-runner/main.cc
```

or their current successors after R4.

### Current evidence

Both:

```text
nmopt_postprocess/records.py
chapter6_report.py
```

implement artifact unescaping and history parsing.

The post-processing reader rejects malformed/non-finite numeric history;
`chapter6_report.py` currently skips some invalid entries.

That means the same persisted evidence can be interpreted differently.

### Semantic requirement

There must be one evidence interpretation.

Do not silently discard malformed or non-finite history values in one
consumer.

### Alternatives

#### A. Standard-library artifact package

For example:

```text
tools/nmopt_artifacts/
  records.py
  run_manifest.py
```

**Advantages**
- clear ownership of persisted schema semantics;
- dependency-free;
- both plotting and reports depend downward on it.

**Costs**
- new package/category;
- some import migration.

#### B. Keep ownership in `nmopt_postprocess`

Make record/manifest imports dependency-light and avoid eager plotting imports.

**Advantages**
- less tree growth;
- smaller migration.

**Costs**
- `nmopt_postprocess` owns semantics broader than post-processing.

#### C. Smaller shared module at `tools/`

Possible, but a loose top-level utility file may become less discoverable than
A or B.

### Decision gate

Present current imports/usages and recommend one alternative before coding.

### Validation

Contracts should cover:

- escape decoding;
- missing/empty history;
- malformed numeric history;
- non-finite values;
- safe/unsafe run-manifest artifact paths;
- manifest statuses;
- report and post-processing agreement on the same fixture.

V1 only.

---

## R8 — Hide large implementation bodies behind stable C++ facades

**Authority:** GREEN for private/static split; AMBER if public include paths change  
**Intervention:** static split  
**Suggested commits:** separate semantic and compiler units

This unit has two independent sub-units. They do not need to be done together.

### R8a — semantic validation implementation

#### Read first

```text
.agents/code.md
.agents/build.md
docs/design/interface-specification.md
docs/internals/compiler/semantic-compiler.md
include/nmopt/semantic/v1/validation.hpp
tests/semantic/
```

#### Current evidence

`validation.hpp` is roughly 3,400 lines, while the public concept is small:

```text
ValidationReport
SemanticValidator::validate(...)
```

The large body contains recognizable validation categories.

#### Suggested implementation

Keep:

```text
include/nmopt/semantic/v1/validation.hpp
```

as the stable public facade.

Move private implementation into a narrow detail area, for example:

```text
include/nmopt/semantic/v1/detail/
  validation_graph.hpp
  validation_pairings.hpp
  validation_formulation.hpp
  ...
```

The exact split should follow actual helper dependencies; do not mechanically
create one file per section heading.

#### Alternatives

1. several detail headers by validation responsibility;
2. one `validation_detail.hpp` first, then split later if that already solves
   navigation;
3. leave as-is if splitting would require invasive circular-include surgery.

The third outcome is acceptable. Large size alone is not sufficient reason to
destabilize the layer.

#### Validation

- focused semantic scenarios;
- diagnostic ordering/content comparisons where tests depend on them;
- neutral compile/build target(s).

Run a full neutral pipeline only when closing the surrounding neutral-refactor
stage.

### R8b — compiled product declarations

#### Read first

```text
docs/design/interface-specification.md
docs/internals/compiler/semantic-compiler.md
include/nmopt/compiler/v1/compiled_problem.hpp
include/nmopt/compiler/v1/dealii_compiler.hpp
tests exercising CompiledProblemT / supplied OTD / KKT / PDAS
```

#### Current evidence

`compiled_problem.hpp` mixes:

```text
manifest/provenance declarations
compiled box/shared data
compiled formulation product types
compilation results
```

#### Suggested implementation

Split along actual declaration dependencies, while keeping:

```text
include/nmopt/compiler/v1/compiled_problem.hpp
```

as a compatibility aggregate.

Possible internal/public helper names:

```text
manifest.hpp
compiled_box_data.hpp
compiled_products.hpp
```

These are suggestions, not required filenames.

#### Alternatives

If the types are too mutually dependent for a clean split, retain a single
header rather than introducing forwarding cycles.

#### Validation

- backend-neutral compile coverage where possible;
- focused deal.II compiler scenarios instantiating each product family;
- no numerical benchmark runs.

### Explicit exclusion

Do not decompose `dealii_compiler.hpp` as part of R8.

---

## R9 — Public include/file naming decisions

**Authority:** AMBER  
**Intervention:** public organization decision

### Orientation

This is a review/decision unit, not an instruction to perform every rename.

By this point, earlier mechanical cleanup should make it easier to judge
whether the remaining public-path ambiguity is worth migration cost.

### Read first

For each candidate, search all current include/users before deciding.

Common authorities:

```text
docs/design/interface-specification.md
docs/reference/application-api.md
docs/internals/compiler/semantic-compiler.md
```

### D1 — semantic reference problems

Current file:

```text
include/nmopt/semantic/v1/reference_specs.hpp
```

Evidence:

- roughly 2,000 lines;
- concrete canonical/reference `ProblemSpec` factories;
- stable delta/reference tests;
- Chapter-5/application consumers;
- stronger role than a simple schema helper.

Alternatives:

1. dedicated `reference_problems/` library area;
2. `semantic/v1/reference_problems.hpp`;
3. keep current file, but stop broad umbrella includes from dragging it in.

### D2 — `problem_spec.hpp`

Current behavior: tiny compatibility aggregate; actual `ProblemSpec` definition
lives elsewhere.

Alternatives:

1. make `problem_spec.hpp` the focused schema include and add explicit
   `all.hpp` for the umbrella;
2. keep `problem_spec.hpp` as compatibility aggregate but document it;
3. rename only after external/public include usage is better understood.

### D3 — algebraic reference models

Current:

```text
include/nmopt/reference/
```

Alternative:

```text
include/nmopt/reference_models/
```

The existing name is broad but not incorrect. Rename only if the clarification
is worth public include churn.

### D4 — reduced solver headers

Inspect:

```text
include/nmopt/solvers/reduced_search.hpp
include/nmopt/solvers/reduced_gradient.hpp
include/nmopt/solvers/reduced_line_search.hpp
include/nmopt/solvers/reduced_trust_region.hpp
```

Possible structure:

```text
reduced_common.hpp
reduced_directions.hpp
reduced_line_search.hpp
reduced_line_search_solver.hpp
reduced_trust_region.hpp
```

Keep old headers forwarding if adopted.

### D5 — compiler application view name

Current:

```text
NativeApplicationViewT
```

Potential ambiguity: “native application” now denotes the independent
application-owned producer path.

Alternatives:

```text
CompiledApplicationViewT
CompilerApplicationViewT
keep NativeApplicationViewT
```

A compatibility alias is likely if renamed.

### D6 — `experiment/reduced_envelope.hpp`

Current top-level `experiment/` contains one file tied to application/benchmark
evidence.

Alternatives:

- move under application evidence/benchmark support;
- retain `experiment/` anticipating genuinely generic experiment facilities.

### D7 — Chapter-specific application C++ placement

Inspect:

```text
include/nmopt/application/chapter5.hpp
include/nmopt/application/chapter6.hpp
include/nmopt/application/dealii/chapter6_b1.hpp
include/nmopt/application/dealii/chapter6_b2.hpp
```

Alternatives:

1. keep them public as showcase/application library material;
2. move them closer to `apps/`;
3. create a separate project-specific library surface.

Do not move them merely because they are Chapter-specific. The user explicitly
views applications as final consumers/showcases of the core.

### Decision protocol

For each D-item actually proposed:

- show current users;
- show migration impact;
- state whether forwarding compatibility is straightforward;
- recommend keep/rename/move;
- get approval before editing.

Do not bundle all accepted D-items into one giant migration if they are
independently reviewable.

---

## R10 — Simplify root CMake after paths stabilize

**Authority:** GREEN if behavior is unchanged  
**Intervention:** static split  
**Suggested commit:** `build(cmake): group repository target registration`

### Orientation

The root `CMakeLists.txt` is not large enough to justify arbitrary
modularization. The issue is that it currently mixes unrelated registration
responsibilities.

Do this **after** major test/source paths settle so CMake is not rewritten
twice.

### Read first

```text
.agents/build.md
CMakeLists.txt
CMakePresets.json
cmake/ScenarioDiscovery.cmake
cmake/WriteDiscoveredScenarios.cmake
```

Also inspect final test/app paths produced by previous units.

### Current responsibility mix

The root file currently covers:

```text
build options/warnings/sanitizers
interface targets
neutral tests
deal.II discovery
nmopt_runner
external Step-4 upstream/stripped/adapted targets
Step-4 minimal consumers
Step-4 verification tests
ordinary deal.II compiler/application tests
Python checks
configuration summary
```

### Suggested shape

A small number of role-based modules is enough, for example:

```text
cmake/BuildOptions.cmake
cmake/Applications.cmake
cmake/TestsNeutral.cmake
cmake/TestsDealii.cmake
cmake/ExternalStep4.cmake
```

These filenames are suggestions.

Preserve the existing scenario-discovery files as-is unless a concrete defect
is found.

### Alternatives

1. extract only Step-4 and test registration, leaving build policy in root;
2. extract build policy plus test/application registration;
3. leave root CMake as-is if previous path cleanup makes it sufficiently
   readable.

Prefer the smallest decomposition that makes the root read like an
architecture map.

### Compatibility

Preserve:

- normal targets;
- options;
- preset behavior;
- labels/timeouts;
- deal.II discovery;
- default build graph.

Do not add a new Step-4 build option in this mechanical unit.

### Validation

This is a good point for a broader gate because registration itself changes:

- fresh `debug-neutral` configure;
- fresh `debug-dealii` configure;
- inspect CTest inventory;
- build representative runner/test/Step-4 targets;
- one full neutral and one full deal.II pipeline after the migration.

---

## R11 — Optional `apps/` taxonomy decision

**Authority:** AMBER  
**Intervention:** move/rename only if selected

### Orientation

There is no requirement to change `apps/`.

The user prefers to keep the small set of showcase consumers together rather
than scattering them into top-level `examples/` or `integrations/`.

Only revisit this if navigation remains confusing after runner/core cleanup.

### Read first

```text
apps/nmopt-runner/
apps/external-dealii/step-4/README.md
docs/reference/application-api.md
docs/reference/external-dealii-solver-integration.md
```

Inspect any Chapter-specific executable entry points that exist by this stage.

### Plausible classifications

If grouping becomes useful, candidate dimensions are:

```text
direct/compiler-driven
configured/runner-driven
external/native
```

Another valid outcome is to keep the current shallow tree.

### Decision rule

Do not create empty or one-item directories for hypothetical future growth.

Present the actual final `apps/` tree before asking for a taxonomy decision.

---

## R12 — Step-4 binding cleanup

**Authority:** AMBER  
**Intervention:** late substantial refactor  
**Suggested commit:** only after explicit approval

### Orientation

Step-4 is valuable partly because it is a coherent, source-accounted
integration snapshot. Do not optimize away that property.

The audit found appealing binding improvements that do not require moving the
snapshots themselves.

### Read first

```text
apps/external-dealii/step-4/README.md
apps/external-dealii/step-4/external-integration-overview.md
apps/external-dealii/step-4/integration-report.md
docs/reference/external-dealii-solver-integration.md
```

Then inspect:

```text
apps/external-dealii/step-4/source/adapted/step-4.cc
apps/external-dealii/step-4/integration/problem_a.hpp
apps/external-dealii/step-4/integration/problem_b.hpp
apps/external-dealii/step-4/integration/problem_b_mass.hpp
apps/external-dealii/step-4/integration/problem_b_metric.hpp
apps/external-dealii/step-4/diagnostics/instrumentation.hpp
apps/external-dealii/step-4/minimal/problem_a_binding.hpp
apps/external-dealii/step-4/minimal/problem_b_binding.hpp
tests/dealii/external_step4_native_contract.cc
tests/application/external_step4_*.cc
```

Re-enumerate current paths first.

### Candidate A — isolate adapted `.cc` inclusion

Several consumers know the pattern:

```cpp
#define STEP4_NO_MAIN
#include "../source/adapted/step-4.cc"
#undef STEP4_NO_MAIN
```

A named integration wrapper could own that unusual seam while preserving the
adapted source verbatim.

### Candidate B — invert instrumentation dependency

Functional Problem A/B code currently refers directly to diagnostic
instrumentation types.

Desired direction:

```text
functional integration
  <- optional observer/instrumentation
```

rather than functional code depending on evaluation machinery.

Use the smallest observer/callback seam that preserves existing evidence.

### Candidate C — local binding duplication

Only reduce duplicated support when doing so does not make the minimal
consumers less pedagogically clear.

### Alternatives

Leaving all three issues unchanged is valid if the churn outweighs the
readability benefit.

### Hard exclusions

Do not change:

- upstream/baseline/adapted source snapshots;
- source line accounting;
- Problem A/B mathematics;
- accepted native/nmopt comparison semantics.

### Validation

Focused Step-4 native/minimal/optimization tests first.

Run one full `debug-dealii` gate before declaring the Step-4 binding work
complete.

No release reproduction run is required.

---

## R13 — Major `dealii_compiler.hpp` decomposition

**Authority:** AMBER and optional  
**Intervention:** high-cost substantial refactor

### Orientation

This is the largest source-readability hotspot, but it is deliberately last.
Earlier test/facade cleanup may make the compiler sufficiently navigable
without touching it.

Do not start R13 just because the file is large.

### Read first

```text
docs/design/interface-specification.md
docs/design/pde-solver-boundary.md
docs/internals/compiler/semantic-compiler.md
include/nmopt/compiler/v1/dealii_compiler.hpp
include/nmopt/compiler/v1/dealii_types.hpp
include/nmopt/compiler/v1/dealii_scalar_plan.hpp
include/nmopt/compiler/v1/compiled_problem.hpp
```

Then use the final compiler test layout produced by R3.

### Current phase map observed by the audit

The compiler body appeared to combine:

```text
semantic/request closure
binding resolution and validation
mesh-dependent checks
target/capability closure
lowerability/formulation/product checks
target-specific policies
FE/numerical realization
formulation-product construction
manifest/provenance construction
```

Reconstruct the current call graph before proposing files.

### Required proposal before implementation

Present:

- current public `DealiiCompiler` API;
- current private phase/call map;
- candidate helper boundaries;
- expected include/dependency direction;
- test scenarios covering each proposed phase;
- whether any extraction affects compile time or diagnostics;
- an option to defer the work.

### Preferred style

Private/detail extraction behind the existing public compiler header.

Do not introduce new public compiler concepts merely to justify file splits.

### Valid alternatives

1. several private phase headers;
2. one or two large detail headers;
3. no decomposition; instead improve source map/comments after earlier cleanup.

The third outcome is acceptable.

### Validation

If implemented:

- focused compiler/deal.II scenario groups after each coherent extraction;
- compare diagnostics/capability behavior;
- one `debug-dealii` stage gate at the end.

---

## R14 — Final structure gate and documentation handoff

**Authority:** GREEN for verification/handoff; AMBER for unresolved design choices

### Orientation

The refactor should stop at a deliberate stable point. It does not need to
implement every optional idea in this roadmap.

### Inspect before testing

- final repository tree;
- `git diff` / branch history for all accepted units;
- public forwarding headers/aliases promised during migrations;
- current CTest inventory;
- remaining `TODO`/historical terminology introduced or exposed by the work;
- stale paths caused specifically by this branch.

### Required final validation

Run:

```bash
./build.sh pipeline debug-neutral
./build.sh pipeline debug-dealii
```

Run all supported Python/tool contracts.

Run `sanitize-neutral` only if ownership/lifetime-sensitive code was changed.

Do not run `release-dealii` without separate authorization.

### Closure record

Add a closure report to the human-readability audit package recording:

- implemented structural changes;
- explicit behavior corrections;
- public compatibility seams retained;
- proposals declined or deferred;
- verification performed;
- path/name changes the resumed docs refactor must reconcile.

### Handoff to `codex/docs-refactor`

After integrating this branch:

- resume the reference rewrite;
- update root/docs routing for the actual final source tree;
- run the broad documentation link/path consistency sweep;
- do not describe deferred optional refactors as implemented architecture.

## 9. Dependency hints, not mandatory order

A likely low-churn path is:

```text
R0
 |
 +--> R1
 |
 +--> R2 ----+
 |           |
 +--> R3 ----+--> later CMake cleanup (R10)
 |
 +--> R4
 |
 +--> R5 --> R6 --> R7 decision
 |
 +--> R8
 |
 +--> R9 decisions
 |
 +--> R11 decision if still useful
 |
 +--> R12 decision
 |
 +--> R13 optional decision
 |
 +--> R14
```

Important non-dependencies:

- R5 does not need the C++ test restructuring.
- R8a does not need R3.
- R7 does not need R9.
- R11 may be skipped entirely.
- R12 and R13 may be deferred beyond this branch if the earlier work already
  achieves sufficient readability.

The implementation agent should prefer the order that minimizes repeated
movement of the same files and repeated expensive validation.

## 10. Suggested prospective commit boundaries

Actual commits remain subject to the repository Git rules and explicit user
authorization.

Likely review-sized boundaries are:

```text
docs(audit): record repository structure follow-up
refactor(readability): remove stale implementation-era wording
test(structure): separate neutral contract and solver coverage
test(dealii): split compiler capability coverage
refactor(runner): extract reusable execution seams
fix(postprocess): honor plotting profile policy
test(tools): register postprocessing contracts
refactor(semantic): hide validator implementation
refactor(compiler): split compiled product declarations
build(cmake): group repository target registration
```

Decision-gated public renames/moves should receive their own commits after
approval rather than being folded into nearby mechanical work.

## 11. Per-unit handoff template

Every implementation agent should finish a unit with:

```text
Unit:
Authority level:
Current branch/head:
Evidence rechecked:
Outcome:
Files changed:
Behavior intentionally changed:
Behavior explicitly preserved:
Compatibility seams:
Focused verification:
Full profile verification:
Unresolved findings:
Roadmap assumptions that changed:
Decision needed before next unit:
Suggested commit:
```

If `Behavior intentionally changed` is nonempty for a unit that was supposed
to be structural-only, stop and explain why before proceeding.

## 12. Stop conditions

Stop immediately and ask when:

- current user modifications overlap the intended files;
- the current source tree materially differs from the roadmap's baseline
  assumption;
- a public path/type/namespace must change outside an approved decision unit;
- a supposedly static move changes test behavior or numerical results;
- source and accepted design documentation disagree materially;
- the proposed fix requires a new permanent build profile or top-level
  repository category;
- a Step-4 change alters preserved source snapshots or source accounting;
- `release-dealii` appears necessary;
- the work unit grows into adjacent cleanup not listed in its scope.

The purpose of this roadmap is to make a weak implementation agent predictable:
the repository is oriented before edits begin, mechanical work is explicit,
architectural choices are gated, and validation is proportional to actual
risk.
