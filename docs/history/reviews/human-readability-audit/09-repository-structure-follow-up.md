# Repository-structure follow-up to the human-readability audit

**Follow-up baseline:** [`997367f37b14`](https://github.com/HunterNise/nmopt-project/commit/997367f37b14da0010e0e474cf1e75690c4d79fb) (`docs(structure): reorganize documentation by role`)  
**Review scope:** non-documentation repository structure and source inspectability after the documentation information-architecture decision  
**Status:** complete assessment; implementation is owned by the active `docs/planning/repository-structure-refactor.md`

## 1. Why this follow-up exists

The original human-readability audit deliberately stopped short of recommending
source reorganization. Its bounded source pass concluded that the documentation
problem could be addressed without mixing code movement into
`codex/docs-refactor`.

After the documentation tree was reorganized, a dedicated follow-up inspected:

- `include/nmopt/`;
- `tests/`;
- `apps/`;
- `tools/`;
- `parameters/`;
- `.agents/`;
- `cmake/`;
- root build, configuration, dependency, and repository-entry files.

That deeper pass found that the original deferral was correct for the
documentation-only branch, but too conservative as a long-term source
conclusion. The project does not need a wholesale redesign, yet several
physical file boundaries, test groupings, and application/tool seams now lag
behind the architecture they implement.

## 2. Constraints fixed by the follow-up discussion

### 2.1 `apps/` means final consumers of the reusable core

`apps/` is the home for concrete users/showcases of the header-heavy nmopt
library: contracts, semantic description, compiler/lowering, backends,
formulations, and solvers.

B1/B2-style configured applications and external Step-4 integration are all
valid application examples in this sense. They demonstrate different ways a
user can consume the same core, analogous to tutorial applications.

Therefore:

- do **not** create a new top-level `integrations/` tree merely to separate
  Step-4;
- do **not** scatter the small number of current examples across many
  top-level categories;
- if a stronger `apps/` taxonomy later becomes useful, plausible dimensions
  are direct/compiler-driven consumption, configured recipe/scenario/parameter
  consumption, and external/native integration;
- choosing such a taxonomy is a design decision, not a mechanical move.

### 2.2 Keep `.agents/explanation-essay.md`

`explanation.md` remains the routed normative instruction. The essay variant is
an intentionally unrouted reference and should not be deleted as duplicate
material.

### 2.3 Preserve the Step-4 package until late

`apps/external-dealii/step-4/` is a coherent source/evaluation snapshot with
useful source accounting. Its internal source/evaluation layout should not be
rearranged early in the refactor.

The binding layer may still be improved later, especially where functional
integration code depends on evaluation instrumentation or several consumers
know how to include the adapted tutorial source.

### 2.4 Documentation rewriting is paused

The detailed reference rewrite should resume only after the repository
structure settles enough that paths and source boundaries are stable.

The repository-structure branch may update documentation that is directly
required to keep a changed interface truthful, but it should not perform the
broad link/reference rewrite already owned by `codex/docs-refactor`.

## 3. Executive verdict

The repository is **architecturally coherent but physically uneven**.

The strongest existing boundaries should be preserved:

- the numerical contract layer;
- most deal.II service headers;
- semantic/compiler versus native-application producer paths;
- application versus reusable numerical-core ownership;
- the scenario-discovery test infrastructure;
- root parameter inputs versus generated run evidence;
- root build presets and machine-local build configuration.

The main structural debt is concentrated in:

1. very large test translation units that have become executable capability
   ledgers;
2. semantic/compiler headers whose public facade is much smaller than their
   implementation body;
3. runner logic still trapped in `main.cc`;
4. generic post-processing code that contains Chapter-6 policy despite an
   accepted configuration contract saying it must not;
5. some umbrella/header names that reflect historical growth rather than the
   current responsibility;
6. root CMake registration that mirrors all repository history in one file.

The goal is therefore **responsibility locality and inspectability**, not a new
architecture.

## 4. `include/nmopt/`

At the follow-up baseline, `include/nmopt/` contained roughly 68 headers and
1.6 MB of source. The size is highly concentrated.

Representative large files were:

| File | Approximate size / lines | Finding |
| --- | ---: | --- |
| `compiler/v1/dealii_compiler.hpp` | 341 KB / 6,969 lines | major orchestration hotspot |
| `semantic/v1/validation.hpp` | 178 KB / 3,456 lines | small public concept, very large implementation body |
| `semantic/v1/reference_specs.hpp` | 98 KB / 1,981 lines | canonical/reference problem library mixed into semantic schema |
| `compiler/v1/dealii_fixed_dirichlet.hpp` | 96 KB / 2,095 lines | large but substantially one responsibility |
| `compiler/v1/compiled_problem.hpp` | 43 KB / 1,125 lines | manifest, shared data, products, and result types mixed |
| `solvers/reduced_search.hpp` | 49 KB / 1,327 lines | common solver policy plus several direction policies |
| `application/chapter6.hpp` | 45 KB / 1,182 lines | deliberately project-specific application material |

### 4.1 Contract layer

`include/nmopt/contract/` is largely healthy. Its files form a recognizable
numerical vocabulary:

```text
linear algebra and layouts
  -> executable model / solves / metrics
  -> formulation contracts
  -> complementarity / PDAS
```

Mathematical distinctions such as primal versus covector, derivative versus
gradient, application versus inverse/solve, and state/adjoint evidence should
not be simplified for cosmetic reasons.

Cheap cleanup is justified where reusable contracts still describe current
behavior as `v0`, `P6.2`, or another implementation-phase label.

### 4.2 deal.II service layer

`include/nmopt/dealii/` is mostly well factored. Individual files for metrics,
constraints, coordinates, serial solves, KKT, and PDAS have clear
responsibilities.

The namespace/directory is broader than the primitive `SerialBackend` storage
policy, but that alone does not justify a namespace migration.

### 4.3 Semantic layer

`semantic/v1/validation.hpp` has a very small conceptual public surface and a
large collection of internal validation phases. A facade plus internal/detail
headers is justified without changing semantic behavior.

`semantic/v1/reference_specs.hpp` is a stronger conceptual placement issue.
Tests and application recipes treat it as a stable library of canonical
problem graphs and feature deltas, not merely the schema definition itself.

Valid alternatives remain:

1. move canonical/reference problems to a dedicated `reference_problems/` or
   equivalent library area;
2. keep them within `semantic/v1/` but rename the file to state its role
   clearly;
3. retain the current placement temporarily and only stop
   `problem_spec.hpp` from pulling the whole reference library implicitly.

Choosing among these changes public include organization and must not be done
mechanically.

`semantic/v1/problem_spec.hpp` is currently a convenience aggregate rather
than the physical definition of `ProblemSpec`. Long-term, either the filename
should become a focused schema include or the aggregate role should be made
explicit through a separate umbrella such as `all.hpp`.

### 4.4 Compiler layer

`compiler/v1/dealii_compiler.hpp` is the largest source-readability hotspot.
Its private body combines several recognizable compiler phases:

```text
request and semantic closure
bindings and data-port validation
target/capability selection
lowerability and mesh-dependent validation
target-specific policy
numerical realization
formulation-product construction
manifest/provenance production
```

The public compiler surface is much smaller than this implementation body.
Decomposition is plausible, but it is a high-cost intervention and should not
be an early mechanical cleanup.

`compiled_problem.hpp` has clearer static seams: manifest/provenance schema,
compiled box/shared data, concrete compiled products, and compilation results.
It is a better low-risk split candidate.

Current compiler diagnostics still contain roadmap-era names such as `P5.1`,
`P5.3`, `C5.6`, `C5.8`, and `C5.10`. Runtime diagnostics should state the
mathematical capability or realization directly. Historical IDs remain valid
in provenance/audit material.

`NativeApplicationViewT` may now be ambiguous because project documentation
uses “native application” for the independent application-owned producer path.
Renaming it to a compiler-produced/application-view name is plausible, but it
is a public naming decision and requires explicit approval.

### 4.5 Solver headers

The solver implementation has evolved beyond some filenames.

`reduced_gradient.hpp` now contains a generic reduced-search solver;
`reduced_search.hpp` contains common policy plus multiple direction families.
A clearer physical structure could expose common records, directions, line
search, the line-search solver, and trust-region policy separately.

This is a public include-path issue. Preserve existing forwarding headers if
it is adopted, and ask before choosing a new canonical naming scheme.

### 4.6 Application and experiment headers

`include/nmopt/application/` mixes reusable recipe/scenario/catalog and
benchmark-harness concepts with explicitly Chapter-5/6-specific material.

That is not automatically wrong: Chapter applications are intentional users
of the reusable core. Whether their C++ definitions belong in the public
`include/nmopt/application/` surface, under `apps/`, or under a project-specific
library namespace remains a design choice.

`include/nmopt/experiment/reduced_envelope.hpp` is a single top-level module
whose known role is application/benchmark evidence. Folding it into an
application evidence/benchmark area is plausible, but not urgent.

### 4.7 Reference models

`include/nmopt/reference/` contains algebraic executable reference/oracle
models used to verify formulation equivalence. Tests demonstrate that these
are more than incidental fixtures.

`reference_models/` would be a clearer name than the broad `reference/`, but a
rename is optional and changes public include paths.

## 5. `tests/`

At the baseline, `tests/` contained roughly 37 files and 1.5 MB.

The custom scenario-discovery infrastructure is a strength. C++ test
executables expose fine-grained scenario metadata, and CTest sees roughly 188
named scenarios rather than one test per translation unit.

The problem is physical and taxonomic, not coverage.

### 5.1 Mixed classification axes

Current directories mix architecture and backend dimensions:

```text
contract/
semantic/
application/
dealii/
tools/
support/
```

`dealii/` can mean compiler realization, backend behavior, or external native
integration. A clearer long-term classification would expose responsibilities
such as contract, solver, semantic, compiler, application, integration/example,
tooling, and support, with deal.II as a subtree where needed.

### 5.2 Giant executable ledgers

`tests/dealii/dealii_diffusion_contract.cc` is nearly 8,000 lines with 34
scenarios. It covers many compiler target families, formulations, metrics,
observations, constraints, diagnostics, sessions, and solve reporting. It is
effectively the executable deal.II compiler capability ledger.

`tests/contract/reduced_dto_contract.cc` is roughly 2,944 lines and mixes DTO
contract checks with line searches, direction policies, stopping, ownership,
projection, experiment-envelope behavior, and backend parameterization.

`tests/semantic/semantic_v1_contract.cc` mixes semantic validation/reference
graphs with resolution and a scalar lowering-plan compiler scenario.

These files should be split by responsibility while preserving logical CTest
identities initially.

### 5.3 External Step-4 tests

Step-4 verification is currently spread across `tests/application/` and
`tests/dealii/`. Conceptually these scenarios belong to one application
example/integration study.

Because the application snapshot itself is intentionally deferred, test moves
that would force Step-4 source restructuring should also be late.

### 5.4 Runner testability

`tests/application/parameter_files_dealii_contract.cc` includes
`apps/nmopt-runner/main.cc` after redefining `main`. This is characterization
evidence that reusable run preparation/execution-registration functions still
live in the entrypoint translation unit.

The source boundary should be fixed rather than normalizing this test trick.

### 5.5 Python tooling tests

Four post-processing contract scripts exist under `tests/tools/` but are not
registered in the root CTest pipeline:

```text
postprocess_configuration_contract.py
postprocess_comparison_plan_contract.py
postprocess_history_contract.py
postprocess_boundary_field_contract.py
```

They cover meaningful configuration, matrix, history, provenance, and
boundary-field behavior. They need an explicit supported verification path.

## 6. `apps/`

The application layer should be understood as concrete consumers/showcases of
the reusable nmopt core.

### 6.1 `nmopt-runner`

`apps/nmopt-runner/` is a genuine application entry point and should remain an
application.

Its helper files already expose useful seams:

- capability and benchmark registries;
- parameter binding and schema work;
- run-set planning;
- runner configuration and manifests.

The main issue is `main.cc`, which still owns reusable run preparation,
execution registration, artifact enrichment, B1/B2 dispatch, configuration
snapshotting, and output helpers. Extracting the minimum reusable seams would
allow tests to call production logic without including the entrypoint source.

Large support headers such as `parameter_files.hpp` and `runner.hpp` also have
static seams, but splitting them is secondary to fixing the `main.cc`
ownership problem.

### 6.2 Step-4

`apps/external-dealii/step-4/` is internally coherent as a case-study snapshot:

```text
source
integration
minimal consumers
evaluation
verification
diagnostics
```

It should remain intact until later in the refactor.

Two later code-level opportunities remain:

1. functional Problem A/B integration currently depends directly on diagnostic
   instrumentation types;
2. several consumers know the `STEP4_NO_MAIN` mechanism for including the
   adapted tutorial `.cc` source.

These can be improved without altering the preserved source snapshots, but
they are not early prerequisites.

### 6.3 Application taxonomy remains optional

If `apps/` eventually needs stronger navigation, valid classifications include:

- direct/compiler-driven consumers;
- configured recipe/scenario/parameter consumers;
- external/native application integrations.

There are few applications today and no evidence that they will grow rapidly.
Do not introduce this hierarchy until the user explicitly chooses it.

## 7. `tools/`

The tool tree contains three distinguishable responsibilities:

1. reusable persisted-artifact and field post-processing;
2. Chapter-6-specific workflow/reporting;
3. external deal.II fidelity utilities.

Physical grouping can improve this later, but one correctness issue takes
priority over moves.

### 7.1 Plotting-profile contract violation

The accepted parameter/plotting-profile decision says the generic renderer
must not contain B1/B2 branches and that JSON plotting profiles own
presentation policy.

Current code violates that contract:

- generic `pipeline.py` recognizes `b2.*` metadata;
- generic `render.py` hard-codes the Chapter-6 `turbo` colormap and related
  normalization/interpolation/mesh-overlay/tick/DPI/axis behavior;
- the checked-in `nmopt-plot-v1` JSON files declare those defaults, but the
  loader does not consume most of them.

This is an explicit behavior correction, not merely a file move. The renderer
should honor validated resolved profile data.

### 7.2 Persisted-record parsing is duplicated

`nmopt_postprocess/records.py` and `chapter6_report.py` independently decode
the same artifact format. Their numeric-history behavior differs: one rejects
malformed/non-finite histories while the report silently drops invalid
entries.

Persisted evidence should have one parsing contract. Two implementation
approaches are valid:

1. create a tiny standard-library-only artifact/manifest package used by both
   reporting and plotting;
2. make the existing record parser dependency-light enough to be imported by
   the report without loading plotting dependencies.

The package choice requires a decision; the semantic requirement does not.

### 7.3 Chapter-6 tooling

`run_chapter6.sh` is useful workflow orchestration but discovers the run
directory by parsing human runner output. A stable machine-readable runner
result would be more robust if the runner is already being refactored.

`chapter6_postprocess.py` is explicitly a compatibility wrapper and should not
drive new architecture.

`chapter6_report.py` is domain-specific by design and need not be generalized.

### 7.4 External deal.II tools

`tools/external_dealii/` exists for Step-4 source/output fidelity. Moving these
utilities into the application snapshot would be coherent, but because Step-4
reorganization is deferred, this move should also be deferred.

## 8. `parameters/`

`parameters/` is coherent and should remain a first-class root input tree.

It represents tracked executable configuration, distinct from:

```text
docs/studies/    human/source interpretation
parameters/      versioned executable inputs
runs/            generated evidence
```

The `authoritative.prm` name means the authoritative **repository reproduction
input**, not that every value is an explicitly recovered publication fact.
A small `parameters/README.md` should make this distinction and explain the
`authoritative` versus `development` split.

The plotting-profile JSON structure itself is appropriate; the defect is in
the tool consumer.

## 9. `.agents/`

The hidden agent-instruction tree has a clear role and should remain separate
from user/project documentation.

Its flat routing structure is small enough to remain readable.

`explanation-essay.md` is intentionally retained as an unrouted reference.
`explanation.md` remains the routed normative instruction.

No directory restructuring is justified.

## 10. CMake and root build surface

### 10.1 Scenario discovery

`cmake/ScenarioDiscovery.cmake` and
`cmake/WriteDiscoveredScenarios.cmake` are well-factored and should be
preserved.

### 10.2 Root `CMakeLists.txt`

The root file is not unmanageably large, but it owns unrelated registration
responsibilities:

```text
build policy
neutral tests
deal.II setup
runner
Step-4 source/minimal targets
Step-4 tests and Python checks
compiler/backend tests
Chapter-6 application tests
```

Once source/test boundaries stabilize, role-based CMake modules can make the
root build file mirror the repository architecture. This should preserve
target names, options, labels, and build semantics.

### 10.3 `build.sh`, presets, and dependencies

`build.sh` is large but coherent as one build CLI and should remain at root.

`CMakePresets.json`, `build.local.conf.example`, and `DEPENDENCIES.md` also
have clear roles and need no structural redesign.

## 11. Intervention classification

### 11.1 Light edits

Low-risk when bounded to current wording/configuration truth:

- remove roadmap-era labels from current reusable runtime/compiler diagnostics;
- improve comments where current behavior is still described as `v0`, `P5.x`,
  `P6.x`, or `C5.x`;
- add `parameters/README.md`;
- fix local formatting/readability defects encountered in touched CMake;
- add concise local routing documentation where a complex application entry
  point lacks it.

Broad documentation-link repair remains owned by the paused documentation
refactor unless a structure change in this branch directly creates a new
broken link.

### 11.2 Static split/merge/move with behavior preserved

Good candidates:

- split large neutral/deal.II test translation units by responsibility while
  preserving scenario names, labels, and assertions;
- split semantic validation implementation behind its current facade;
- split `compiled_problem.hpp` behind compatibility includes;
- extract runner-owned reusable functions from `main.cc`;
- later group CMake registration into role-based modules.

### 11.3 Substantial refactoring with an accepted direction

The plotting-profile correction is required by an already accepted design:

```text
versioned profile
  -> validated resolved presentation policy
  -> generic pipeline/renderer
```

The generic engine must not contain B1/B2 policy.

Runner extraction also has an accepted boundary: production logic used by
tests must live outside the program entrypoint. The exact private file split
should remain minimal.

### 11.4 Decision-gated changes

Ask before implementing any of these:

- exact `apps/` subdirectory taxonomy;
- relocation/renaming of semantic reference problems;
- canonical meaning of `problem_spec.hpp` versus a new umbrella header;
- `reference/` to `reference_models/`;
- `NativeApplicationViewT` rename;
- solver public-header rename/reorganization;
- moving `experiment/reduced_envelope.hpp`;
- relocating Chapter-5/6 application headers;
- choosing the new shared Python artifact-parser package boundary;
- major physical decomposition of `dealii_compiler.hpp`;
- Step-4 internal restructuring or instrumentation architecture.

### 11.5 Explicitly deferred unless new evidence appears

- redesigning numerical contracts;
- collapsing mathematical distinctions for API simplicity;
- introducing compiled `.cc` implementation merely to reduce header size;
- broad namespace redesign;
- deleting compatibility seams before proving they are unused;
- restructuring preserved Step-4 source snapshots;
- release-build or benchmark reproduction work unrelated to structure.

## 12. Behavior-preservation principle

Unless a roadmap unit explicitly identifies a correctness defect, structure
work must preserve:

- numerical semantics and formulas;
- public types and namespaces;
- existing include paths, using forwarding headers where necessary;
- CLI flags and exit behavior;
- parameter schemas and precedence;
- artifact/run-manifest schemas;
- CMake target names and build options;
- CTest scenario names, labels, and timeouts;
- benchmark/application identities and provenance.

The two already identified behavior corrections are:

1. plotting profiles must actually control the presentation properties they
   declare, with Chapter-specific policy removed from the generic engine;
2. persisted artifact/history parsing should not silently produce a different
   evidentiary interpretation in reporting than in post-processing.

Any additional behavior change discovered during restructuring must be stopped
and presented as a separate decision before implementation.

## 13. Handoff

Implementation order, decision gates, compatibility rules, and validation
cadence are defined in the active roadmap
`docs/planning/repository-structure-refactor.md`.

The purpose of that roadmap is not to implement every optional cleanup listed
here. It is to make the repository easier to inspect while keeping the
architecture and observable behavior stable.
