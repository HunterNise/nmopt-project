# Repository-structure refactor closure

**Branch:** `codex/repo-structure-refactor`  
**Base HEAD:** `0881ee9f768c8e36099939a3261b841830402b0b`
(`docs(step4): refresh external integration study`)  
**Governing roadmap:** [`docs/planning/repository-structure-refactor.md`](../../../planning/repository-structure-refactor.md)  
**Follow-up assessment:** [`09-repository-structure-follow-up.md`](09-repository-structure-follow-up.md)

## 1. Result

The repository-structure refactor is complete.

The work reorganized source, tests, tooling, application integration, and build
registration where the repository's existing responsibilities had become hard
to inspect. It did not redesign the numerical architecture. In particular, the
refactor preserves the distinctions between:

- primal and covector values;
- derivatives and gradients;
- operator actions and inverse/solve operations;
- residual, JVP, and VJP actions;
- state and adjoint solves;
- numerical results and solve/evidence records;
- compiler-owned and application-owned numerical products;
- reusable numerical contracts and runner/reproduction infrastructure;
- reduced DTO, supplied OTD, KKT, complementarity, and PDAS products.

The semantic/compiler producer path and the application-owned/native producer
path remain peers. They converge at shared numerical/formulation contracts,
not at `ProblemSpec` or `CompiledProblemT`. No universal formulation interface
was introduced.

Large-file size alone was not treated as sufficient evidence for a physical
split. Changes were made where a responsibility, ownership boundary, or
supported workflow could be expressed more clearly without flattening those
architectural distinctions.

## 2. Completed work and rationale

### R1 — current terminology

Current runtime/compiler diagnostics were moved away from roadmap-era
`P5.x`/`P6.x`/`C5.x`/`v0` vocabulary where that vocabulary was required to
understand present behavior.

Historical/study identifiers remain where they are intentional provenance or
test traceability.

The final readability pass also changed the remaining explanatory
`problem_library.hpp` comments to capability-first language without changing
semantic requirement strings, persisted identifiers, diagnostic strings, or
test assertions.

### R2 — backend-neutral test responsibility

The neutral test cleanup moved two misplaced scenario bodies to
responsibility-named implementation headers while preserving the existing
translation units, scenario registries, CTest identities, labels, timeouts,
and behavior.

This deliberately avoided multiplying header-heavy compile units merely to
make source files smaller.

### R3 — deal.II compiler and backend tests

The former deal.II capability ledger was separated by responsibility:

- compiler capability scenarios live under `tests/compiler/dealii/`;
- lower-level deal.II service/session/reporting coverage remains under
  `tests/dealii/`.

The compiler coverage intentionally retains one heavy compiler translation
unit with capability-grouped scenario headers. Splitting it into many
translation units would repeatedly compile the large header-only compiler and
deal.II template surface without improving the logical test inventory.

Routine regression tests were also separated from retained
extended/reproduction evidence through existing labels and pipeline policy.

### R4 — runner execution seam

Chapter-6 execution behavior was extracted from the executable entry point into
`apps/nmopt-runner/chapter6_execution.hpp`.

The goal was to give tests a callable production seam without including
`main.cc`, while keeping CLI/configuration/lifecycle ownership in the
application. The resulting private header remains substantial but represents
one coherent application responsibility; no generic execution framework was
introduced.

### R5 — plotting-profile correctness

The post-processing work included one intentional behavior correction:
resolved plotting profiles now control the presentation properties they
declare instead of leaving Chapter-specific presentation policy embedded in
generic rendering code.

Existing explicit override precedence was preserved.

### R6 — supported Python/tooling contracts

Supported Python contract scripts were registered in the test workflow with
dependency-aware gating:

- stdlib-only tooling contracts require Python;
- plotting contracts remain conditional on plotting dependencies;
- CMake does not create environments or install optional packages.

### R7 — shared persisted-artifact parsing

Dependency-free artifact record parsing was centralized under
`tools/nmopt_artifacts/`.

This was both a structural cleanup and an evidence-semantics correction:
reporting and post-processing no longer maintain separate interpretations of
malformed, empty-token, or non-finite persisted history. The existing stricter
history interpretation became canonical.

`nmopt_postprocess` retains its compatibility adapter; unrelated run-manifest
parsing was not folded into the new package.

### R8 — semantic and compiler product boundaries

Semantic validation now has a narrow public facade:

```text
include/nmopt/semantic/v1/validation.hpp
include/nmopt/semantic/v1/detail/validation_detail.hpp
```

The large private implementation remains one inline unit because validation
order and helpers are tightly coupled. The final readability pass added only a
small source map inside that private implementation.

Compiler declarations were separated along a stronger natural boundary:

```text
compilation_manifest.hpp
compiled_products.hpp
detail/compiled_kkt_product.hpp
compiled_problem.hpp
```

The manifest is an independently meaningful typed/provenance record; executable
compiled products and ownership-sensitive product state remain separate.
`compiled_problem.hpp` remains the public compatibility aggregate.

### R9 — naming, ownership, and test hygiene

R9 was an intentional repository-wide naming/ownership migration rather than a
mechanical file shuffle.

The three former public algebraic `include/nmopt/reference/*` models were found
to be test/oracle infrastructure and moved to:

```text
tests/support/reference_models/
```

They were not re-exposed through forwarding public headers because the
architectural decision was that they are not part of the public framework
surface.

The semantic recipe library was renamed:

```text
include/nmopt/semantic/v1/reference_specs.hpp
  -> include/nmopt/semantic/v1/problem_library.hpp
```

Public problem-factory function names were preserved. The old include path was
not retained; this is an intentional include-path migration, not a claim of
path compatibility.

The compiler-owned application view was renamed:

```text
native_application_view.hpp / NativeApplicationViewT
  -> compiled_application_view.hpp / CompiledApplicationViewT
```

This avoids using “native application” for a compiler-produced object now that
the native/application-owned producer path has a precise architectural
meaning. No compatibility alias was retained for this internal
repository-wide migration. Historical CTest IDs/labels remain where they are
intentional inventory compatibility.

The app-local generic `runner.hpp` name became `run_lifecycle.hpp`, separating
application lifecycle support from the reusable numerical runner abstraction.

Ordinary test scratch paths were also changed from fixed shared locations to
scoped unique temporary directories.

### R10 — build organization

Only the exceptional Step-4 build-registration block was extracted:

```text
cmake/ExternalStep4.cmake
```

This was deliberate. Step-4 is a coherent mini-ecosystem with preserved source
targets, minimal consumers, Python comparison requirements, run directories,
direct tests, and evidence/reproduction labels.

The remaining root `CMakeLists.txt` continues to expose useful top-level build
policy and ordinary test/application registration. Additional modules were not
created merely for symmetry. Navigation comments were added after the real
extraction so the root file reads as an architecture map rather than a set of
arbitrary fragments.

A proposed CMake link-signature style cleanup was reverted because the existing
deal.II integration requires the compatible plain signature.

### R11 — `apps/` taxonomy

No new `apps/` taxonomy was introduced.

The current application set is small enough that another directory level would
mostly create one-item categories. `apps/` remains the clear home for concrete
final consumers/showcases, including `nmopt-runner` and the external Step-4
case.

### R12 — external Step-4 integration

The final Step-4 structure separates responsibilities explicitly:

```text
source/        preserved source snapshots
integration/   application-owned OCP/native reuse
minimal/       canonical nmopt external consumers
evaluation/    instrumented nmopt/native comparison paths
verification/  mathematical/derivative/oracle checks
diagnostics/   optional counters/evidence support
```

Instrumented comparison bindings were moved out of the apparent canonical
integration path so that `minimal/` is the unambiguous external-consumer
example.

Repeated `STEP4_NO_MAIN` inclusion knowledge was consolidated behind:

```text
integration/adapted_step4.hpp
```

The preserved source snapshots were not restructured.

Two tempting abstractions were deliberately not introduced:

1. no observer/event framework was added solely to erase the optional
   `Instrumentation` type from the source include graph;
2. Problem A and Problem B bindings were not hidden behind a shared custom
   adapter framework merely to reduce line count.

The separate visible bindings are part of the experiment's evidence: Problem B
is mathematically richer than Problem A, while the nmopt adaptation burden
remains of similar scale.

The canonical minimal consumers do have a transitive source-level dependency
on the optional instrumentation type through application integration headers.
They do **not** instantiate instrumentation at runtime, execute instrumentation
recording branches, or produce instrumentation records. They also do not adopt
the native reference optimizer, dense oracle, comparison runner, or
verification machinery. The final documentation pass makes this distinction
explicit.

The fresh full Step-4 reproduction was completed at `170c9f1`; it was not
repeated during the final comments/documentation-only closure gate.

### R13 — `dealii_compiler.hpp`

`include/nmopt/compiler/v1/dealii_compiler.hpp` remains intentionally large.

The review identified plausible physical regions, including registered-graph
validation and manifest/provenance projection, but those regions remain
coupled to surrounding request state, semantic lookups, lowerability,
realization data, and final product construction.

A mechanical `detail/*.hpp` split would:

- improve physical file size/navigation;
- not materially reduce template parse/rebuild cost;
- add cross-file navigation;
- risk presenting artificial boundaries as architectural ones.

Future extraction should wait for a compiler phase with a stable typed
input/output boundary, narrower dependencies, independent focused coverage,
or measured compile-time benefit.

The low-risk readability recommendation was implemented instead: the final
worktree adds coarse navigation comments for request construction, compiler
orchestration, validation/realization phases, semantic lookup, metadata and
provenance, manifest projection, and final KKT/PDAS products.

## 3. Compatibility and migration ledger

The branch does not claim universal old-path compatibility.

### Compatibility intentionally retained

- `compiled_problem.hpp` remains an aggregate after the manifest/product split;
- logical CTest names, labels, timeouts, and scenario identities remain stable
  across physical test-source reorganization;
- historical compiled-view CTest IDs/labels remain intentionally;
- Python/post-processing compatibility adapters remain where they represent
  supported persisted formats/workflows.

### Intentional path/name migrations

```text
semantic/v1/reference_specs.hpp
  -> semantic/v1/problem_library.hpp

native_application_view.hpp / NativeApplicationViewT
  -> compiled_application_view.hpp / CompiledApplicationViewT
```

These migrations were repository-wide and deliberate; the old public/internal
names were not restored solely to make the migration mechanically
backward-compatible.

### Public surface intentionally removed

```text
include/nmopt/reference/*
  -> tests/support/reference_models/*
```

This is an ownership correction, not a public-header relocation. Reintroducing
forwarders would recreate the public surface that R9 intentionally removed.

## 4. Deliberate deferrals

The following were considered and intentionally left unchanged:

- further physical decomposition of `dealii_compiler.hpp`;
- further splitting of `validation_detail.hpp`;
- a Step-4 instrumentation observer/event abstraction;
- a common Problem A/B Step-4 binding framework;
- additional root-CMake modules;
- a new `apps/` taxonomy;
- reduced-solver public-header reorganization;
- a new physical taxonomy for the semantic problem library;
- relocation of `experiment/reduced_envelope.hpp`;
- relocation of Chapter-specific application headers.

These are not unfinished R14 tasks.

## 5. Final readability pass

After reconstructing the original R1-R13 rationale, one bounded final pass was
applied on top of base HEAD `0881ee9`.

Exactly six tracked files changed:

```text
apps/external-dealii/step-4/external-integration-overview.md
apps/external-dealii/step-4/integration-report.md
include/nmopt/compiler/v1/dealii_compiler.hpp
include/nmopt/semantic/v1/detail/validation_detail.hpp
include/nmopt/semantic/v1/problem_library.hpp
tests/compiler/scenarios/scalar_lowering_plan.hpp
```

The C++ edits are comments only:

- source navigation in `dealii_compiler.hpp`;
- source navigation in `validation_detail.hpp`;
- capability-first explanatory comments in `problem_library.hpp`;
- one capability-first compiler-scenario comment.

No executable C++ statement, declaration, include, string literal, assertion,
or build registration changed.

The Step-4 Markdown edits only clarify the distinction between a transitive
source-level optional instrumentation type dependency and the absence of
active instrumentation/evaluation machinery in canonical minimal execution.

`git diff --check` passes.

## 6. Final validation

Validation was run against the final six-file worktree.

### Python/tooling contracts

All eight direct contract scripts passed:

```text
artifact_records_contract.py
chapter6_report_contract.py
external_dealii_forward_contract.py
postprocess_boundary_field_contract.py
postprocess_comparison_plan_contract.py
postprocess_configuration_contract.py
postprocess_history_contract.py
postprocess_render_policy_contract.py
```

The rendering contracts emitted only Matplotlib temporary-cache notices.

### Backend-neutral pipeline

```bash
./build.sh pipeline debug-neutral
```

Result:

- configure passed;
- build passed;
- 74/74 tests passed;
- no skipped/not-built tests were reported.

### deal.II pipeline

```bash
./build.sh pipeline debug-dealii
```

Result:

- configure passed;
- build passed with the repository-configured single build job;
- configured CTest inventory reached 202;
- the normal pipeline selected 178 tests;
- 178/178 selected tests passed.

The normal policy excluded 10 `extended` and 14 `reproduction` tests.

No `release-dealii`, sanitizer, repository-wide extended/reproduction, or full
Step-4 reproduction run was started during this final gate.

Sanitizers were not rerun because the final cleanup contains no executable,
ownership, or lifetime-sensitive implementation change.

A final `git diff --check` passed, and validation did not modify the six tracked
files.

## 7. Remaining documentation handoff

The final structure review found several stale documentation links inherited
from the earlier documentation reorganization:

- `README.md:130` -> absent `docs/design/system-blueprint.md`;
- `docs/manual/concepts/10-validation-resolution-and-capabilities.md:1395`
  -> absent `docs/implementation/v1/semantic-compiler.md`;
- `docs/manual/concepts/11-compilation-and-lowering.md:2048`
  -> absent `docs/implementation/v1/semantic-compiler.md`;
- `docs/manual/overview/architecture-maps.md:165`
  -> absent `docs/implementation/v1/semantic-compiler.md`;
- `docs/manual/overview/README.md:104`
  -> absent `docs/implementation/v1/semantic-compiler.md`;
- `docs/manual/overview/semantic-compiler.md:358`
  -> absent `docs/implementation/v1/semantic-compiler.md`;
- `docs/planning/implementation-roadmap.md:386`
  -> absent `docs/implementation/v1/semantic-compiler.md`;
- `docs/manual/overview/README.md:111`
  -> absent `docs/planning/review/`.

These are documentation-refactor residuals rather than repository-structure
defects and should be handled when `codex/docs-refactor` resumes.

No current source uses the removed public `include/nmopt/reference/` tree or
the obsolete Step-4 binding paths. Remaining compatibility/legacy terminology
is either intentional adapter/provenance language, mathematical terminology,
known limitation text, or retained test-inventory identity.

## 8. Closure

No technical issue remains that should block declaring the repository-structure
refactor complete.

The final state has:

- clearer responsibility locality without redesigning numerical contracts;
- test organization that reflects compiler/backend/application ownership while
  preserving logical test identity;
- explicit current naming for problem libraries and compiler-owned application
  views;
- a small public semantic-validation facade and separated compiler
  manifest/product declarations;
- a dependency-free shared artifact-record parser;
- a clearer Step-4 consumer/evaluation boundary;
- a deliberately limited CMake extraction;
- documented reasons for the structural changes that were *not* made;
- final neutral, deal.II, and Python/tooling validation passing.

The appropriate closure action is to keep this file as the permanent audit
record, remove the temporary uncommitted
`11-r14-detailed-review-report.md` fresh-reviewer handoff, and commit the final
six-file readability/documentation pass together with this closure record.

Suggested commit message:

```text
docs(audit): close repository structure refactor
```
