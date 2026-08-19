# Pre-Chapter 5/6 refactor assessment plan

## Purpose

This document is the execution handoff for assessing whether `nmopt` needs a
refactor before implementing selected problems and numerical methods from
Chapters 5 and 6. The assessment must be exhaustive across architecture,
code, tests, documentation, build configuration, and contributor or agent
instructions, while remaining inside the scope already established by the
project documentation.

The immediate deliverable is an evidence-backed assessment, not a code
refactor. The assessment may recommend that no broad refactor is needed. It
must distinguish structural work required by the selected Chapter 5/6 targets
from cleanup that is merely attractive in a more general framework.

The governing design principle remains composition of residual, objective,
metric, constraint, formulation, and discretisation components. The
assessment must not use “future extensibility” to broaden the project beyond
the documented Chapter 5/6 endpoint.

## Authority and required routing

This plan does not replace agent instructions or architectural contracts.
Before every repository action, follow the routing rules in the
[agent instructions](../../../../.agents/README.md). Use the
[documentation map](../../../README.md) to select authoritative design documents.

At minimum, the assessment agent must:

1. Read `.agents/README.md` first.
2. Inspect `git status` before making any change.
3. Read the instruction files required for the next action.
4. Read authoritative documentation before reviewing its corresponding code.
5. Work progressively. Do not begin by reading every code and documentation
   file in one pass.

If this plan conflicts with `AGENTS.md`, a convention, or a normative contract,
the more authoritative source wins and the conflict must be reported.

## Branch and baseline policy

All assessment and later refactor work belongs on:

```text
codex/ch5-ch6-development
```

At the time this plan was created, the branch started from commit `7c2496b`,
which is also marked by the annotated tag:

```text
pre-refactor-ch5-ch6
```

The tagged commit is the immutable comparison baseline. The branch is the
moving integration line for the assessment and any user-approved refactors.

Every agent must verify the branch and baseline before changing files:

```bash
git status --short --branch
git log -1 --oneline --decorate
git show --no-patch --oneline pre-refactor-ch5-ch6
```

Required branch discipline:

- Do not perform assessment or refactor edits on `codex/main`.
- If the worktree is not already on `codex/ch5-ch6-development`, stop
  and ask the user or worktree owner to place it there. Do not switch branches
  for convenience.
- Do not move or recreate `pre-refactor-ch5-ch6`.
- Do not push, merge, rebase, reset, or commit unless the user explicitly
  requests an operation permitted by the repository conventions.
- Preserve the tag as the behavioral and structural comparison point until
  the selected Chapter 5/6 work is complete.
- Keep feature implementation out of refactor changes. A commit or review unit
  should be either an assessment, a behavior-preserving refactor, a contract
  change, or a feature — never an opaque mixture.

The branch may eventually contain several coherent refactor batches, but do
not create additional branches or worktrees unless the user requests them.

## Scope boundary

The assessment is bounded by the capabilities already recorded in:

- the [Chapter 5 implementation guide](../../../guides/chapter-5-elliptic-control.md);
- the [Chapter 6 numerical-methods guide](../../../guides/chapter-6-numerical-methods.md);
- the [Chapter 6 numerical-examples reference](../../../guides/chapter-6-numerical-examples.md);
- P5.1–P5.6 and P6.1–P6.5 of the
  [implementation roadmap](../../implementation-roadmap.md).

The project is not required to implement every catalogue entry. Until the
user selects exact targets, the assessment must present conditional findings:
for example, a mixed-block issue may be mandatory for a Stokes target but
irrelevant to a selected scalar Laplace target.

The following are not valid reasons to refactor by themselves:

- supporting arbitrary PDEs, backends, discretisations, or execution modes;
- preparing for chapters or mathematical features not listed in the current
  documentation;
- making every component maximally generic;
- replacing a working design solely because another pattern is fashionable;
- file length without demonstrated loss of cohesion, duplicated policy,
  testing difficulty, or Chapter 5/6 change amplification;
- introducing packaging or deployment machinery with no project use case.

Explicit exclusions already stated by the Chapter 6 guide, including its
error-estimation sections and later-chapter nonlinear or SQP concerns, remain
out of scope.

## Mandatory phase gate

Work is split into two separately authorized stages.

### Stage A — Assessment and documentation

Stage A is read-only with respect to code, tests, CMake, conventions, and
agent instructions. It may create or update only the assessment document and
documentation-map entry requested by the user.

The expected assessment file is:

```text
docs/planning/review/pre-ch5-ch6/assessment.md
```

Do not implement a proposed refactor during Stage A, even if the improvement
appears obvious.

### Stage B — User-approved refactors

Stage B begins only after the user reviews the assessment and selects the
recommendations to implement. Each selected refactor needs its own bounded
plan, affected contract review, tests, and verification. An unresolved
architectural choice must be documented before an interface is silently
extended.

## Assessment questions

The report must answer all of the following:

1. Which parts of the current system are well designed and should be
   preserved?
2. Which problems are correctness defects, which are design debt, which are
   test or documentation debt, and which are simply unimplemented features?
3. Which findings increase the cost or risk of the documented Chapter 5/6
   targets?
4. Which refactors are required before any new feature, which are required
   only before a particular capability, and which can safely be deferred?
5. Is the current v0-reference/v1-compiler split still providing independent
   verification, or has it become duplicated implementation without adequate
   value?
6. Are semantic declarations, executable contracts, compiler policy,
   deal.II realization, formulations, and solvers still separated in code as
   specified in documentation?
7. Do the tests independently protect the mathematical contracts and provide
   useful failure localization?
8. Do build and contributor workflows make correct changes easy to verify?
9. What is the smallest coherent refactor sequence that improves engineering
   quality without expanding project scope?

## Phase 0 — Reproduce the starting state

Record the following in the assessment before drawing architectural
conclusions:

- active branch, baseline tag, current commit, and worktree status;
- compiler, CMake, Ninja, and deal.II versions or availability;
- configured targets and registered CTest tests;
- backend-neutral build and test results;
- deal.II build and test results when deal.II is available;
- optional-backend behavior when deal.II is unavailable or disabled;
- any existing warnings, failures, skips, or nondeterministic results.

Use the baseline commands from the build convention:

```bash
cmake -S . -B build -G Ninja -DBUILD_TESTING=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

Also verify the backend-neutral configuration in a separate ignored build
directory when useful:

```bash
cmake -S . -B build-nodealii -G Ninja \
  -DBUILD_TESTING=ON -DNMOPT_ENABLE_DEAL_II=OFF
cmake --build build-nodealii
ctest --test-dir build-nodealii --output-on-failure
```

Do not change configuration or suppress a failing test to manufacture a clean
baseline. Report optional dependency limitations explicitly.

## Phase 1 — Build the scope and traceability matrix

Before deep code review, create a matrix connecting the documented endpoint
to the current implementation. Use one row per reusable capability rather
than one row per complete textbook problem.

The matrix should include:

| Field | Meaning |
| --- | --- |
| Capability | Residual term, observation, metric, constraint, formulation, solver service, or discrete policy |
| Documentation authority | Contract or guide defining the capability |
| Current implementation | Existing types, lowerers, tests, and explicit exclusions |
| Chapter target | P5/P6 request or example that needs it |
| Necessity | Common prerequisite, target-dependent, deferred, or out of scope |
| Refactor exposure | Existing ownership or composition concern that affects adding it |
| Evidence needed | Code, test, build, or documentation checks required before judgment |

At minimum, cover the capability families in P5.1–P5.6 and P6.1–P6.5. Do not
assume that all of them will be selected for implementation.

This matrix is the filter for every later recommendation. If a proposed
refactor has no correctness benefit, engineering-learning benefit, or
documented Chapter 5/6 consequence, classify it as deferred rather than
expanding the scope.

## Phase 2 — Inventory the repository progressively

Construct a factual inventory before evaluating design quality:

- tracked source, test, build, documentation, convention, and instruction
  files;
- CMake targets, interface dependencies, compile features, and test
  registration;
- public include boundaries and major include dependencies;
- line counts and unusually large units;
- duplicated names, policies, diagnostics, and assembly patterns;
- recent change concentration and files repeatedly changed together;
- generated or local build artifacts and ignore behavior;
- agent-facing files under `.agents/` and `.codex/`, where relevant.

Use `rg` and `rg --files` for repository searches. File size and history churn
are navigation aids, not proof of a design problem.

Do not read all documentation deeply during this phase. Use headings, links,
and the documentation map to select the authoritative material required by
each later review wave.

## Phase 3 — Audit architecture and C++ by layer

Review the code in dependency order so that higher-level judgments are based
on established lower-level contracts.

### 3.1 Typed algebra and executable contracts

Review `include/nmopt/contract/` and the dense reference model. Check:

- primal/covector separation and layout validation;
- whether pairing is the sole primitive dual action;
- backend capability assumptions and accidental dense-specific behavior;
- ownership, copying, lifetime, and error behavior of executable callbacks;
- residual, JVP, VJP, objective, and objective-derivative consistency;
- reduced DTO assumptions about state, decision, and equation blocks;
- separation of objective regularisation, metric, and constraint;
- whether current interfaces can support only the documented next steps
  without speculative generalization;
- whether the dense reference remains an independent oracle.

Record both strong invariants worth preserving and limitations that affect a
documented P5/P6 capability.

### 3.2 Semantic graph and validation

Review `include/nmopt/semantic/v1/` against the interface specification and
semantic-compiler record. Check:

- component identity and port ownership;
- structural, policy, lowerability, and formulation-capability diagnostic
  separation;
- invalid or ambiguous graph states that types permit;
- consistency and duplication among reference specifications;
- whether problem factories describe compositions or encode complete
  application classes indirectly;
- region, space, pairing, transformation, observation, and requirement-policy
  traceability;
- deterministic and actionable diagnostics;
- change amplification expected for selected P5/P6 capabilities.

Do not label an unsupported semantic construct as refactor debt when it is
simply a future feature with a correct capability diagnostic.

### 3.3 Compiler package, registry, and manifest

Review the generic v1 compiler files before specialized deal.II targets.
Check:

- dependency direction between semantic declarations, capability checks,
  registered lowerers, compiled packages, and backend types;
- whether the registry actually composes independently declared components;
- repeated whole-problem dispatch or branching that should belong to a term,
  observation, metric, constraint, or formulation service;
- manifest completeness and provenance consistency;
- data binding, cache invalidation, and recompilation boundaries;
- public versus private compiler responsibilities;
- compile-time coupling and include fan-out;
- diagnostic ownership and duplication;
- v0/v1 comparison boundaries and independence.

Pay particular attention to the responsibilities concentrated in
`dealii_compiler.hpp`, but do not recommend splitting it solely because of its
size.

### 3.4 deal.II realizations

Review the direct v0 lowerer, backend adapter, metrics, constraints, and every
registered v1 assembled target. Compare implementations structurally and
mathematically. Check:

- repeated mesh, DoF, quadrature, assembly, constraint, state-solve, and
  adjoint-solve machinery;
- whether apparent duplication encodes genuinely different policies;
- whether reusable assembly actions can be extracted without erasing
  mathematical distinctions;
- physical versus independent state coordinates and all reconstruction
  pullbacks;
- residual signs and JVP/VJP transpose consistency;
- state-dependent reassembly and cache lifetime;
- nullspace, trace, boundary partition, and compatibility policies;
- solver tolerances and diagnostics;
- ownership of metrics, constraints, bounds, and data bindings;
- assumptions that block a selected Chapter 5/6 target;
- assumptions that should remain deliberately narrow.

Any proposed extraction must identify the exact shared invariant and the
tests that would prove behavior preservation.

### 3.5 Formulations and solvers

Review reduced DTO orchestration and the reduced-gradient solver against the
Chapter 6 guide. Check:

- solver dependence only on declared executable ports;
- line-search, stopping, projection, and diagnostic ownership;
- state/adjoint solve accounting;
- distinction between a covector, a metric direction, and a projected
  stationarity measure;
- hard-coded assumptions that matter to P6.1–P6.5;
- whether Hessian, KKT, preconditioner, and complementarity work needs new
  explicit contracts rather than extensions to reduced gradient;
- deterministic behavior and clear failure modes.

Do not turn the current solver into a universal algorithm dispatcher. Prefer
small reusable policies only when a documented method requires them.

## Phase 4 — Audit tests and numerical verification

Inventory each test case and map it to a contract. The report must cover:

- residual value checks;
- residual JVP finite-difference or Taylor checks;
- residual VJP pairing checks;
- objective and objective-derivative checks;
- physical reconstruction value, JVP, and VJP checks;
- state and adjoint residual checks;
- state-recomputed reduced Taylor checks;
- metric apply/inverse-apply and pairing checks;
- constraint feasibility and projection checks;
- solver descent, line-search, stopping, and diagnostic checks;
- semantic, policy, lowerability, and formulation diagnostic checks;
- compilation manifest and provenance checks;
- v0/v1 equivalence checks where independence is meaningful.

Assess the following engineering properties:

- whether tests derive expected values independently or repeat production
  formulas;
- failure localization and CTest granularity;
- organization of the large deal.II contract executable;
- fixture/helper duplication and hidden shared state;
- tolerance derivation, scale dependence, and mesh dependence;
- deterministic seeds and iteration counts;
- positive, negative, boundary, and degenerate cases;
- backend-neutral versus deal.II-specific coverage;
- separation of correctness tests from performance or parameter sweeps;
- missing regression tests for already fixed defects;
- practical ability to test a refactor without adding a new application
  class.

Do not introduce a coverage service or new test framework during Stage A.
Report whether such tooling would provide proportionate value.

## Phase 5 — Audit build configuration and tooling

Review the root CMake configuration and observed configuration behavior.
Check:

- target boundaries and transitive usage requirements;
- include-directory build/install semantics appropriate to this project;
- C++ standard and compiler-feature declarations;
- deal.II package and fallback discovery behavior;
- explicit reporting of disabled or unavailable backend tests;
- test-target granularity and registration;
- warning visibility and whether warnings are actionable;
- configuration reproducibility and generated-artifact isolation;
- whether header-heavy compilation produces avoidable coupling or cost;
- feasibility of backend-neutral checks without deal.II;
- whether automation is absent, stale, excessive, or proportionate to the
  project’s software-engineering goal.

Do not recommend packaging, installation, CI, formatting, static analysis, or
sanitizer machinery automatically. State the concrete failure mode or
learning benefit each proposed tool would address.

## Phase 6 — Audit documentation and agent guidance

Review documentation by authority and purpose rather than linearly. Check:

- consistency between interface specification, implementation-readiness
  choices, executable contracts, compiler records, and implementation;
- current-status and completed-task claims in the system blueprint and
  roadmap;
- Chapter 5/6 capability classifications and implementation prerequisites;
- stale next-agent handoffs or references to already completed work;
- terminology, signs, spaces, pairings, and metric/regularisation distinction;
- source-code links, document-index coverage, and broken relative links;
- whether `README.md` gives a reliable newcomer path;
- whether `AGENTS.md`, conventions, `.agents/`, and `.codex/` instructions are
  consistent, minimal, and operational;
- duplicated authority or instructions that can drift;
- missing decision records for behavior already embedded in code.

Do not rewrite mathematical background for stylistic consistency. Prioritize
errors that could cause an agent or contributor to implement the wrong
contract, repeat completed work, or bypass a layer.

## Phase 7 — Synthesize findings

The assessment must lead with a concrete verdict and present strengths before
or alongside improvements. Use stable finding identifiers such as `RF-001`.

Each improvement finding must include:

| Field | Required content |
| --- | --- |
| ID and title | Stable identifier and concise outcome-oriented title |
| Classification | Correctness defect, architectural debt, code-quality debt, test debt, build/tooling debt, documentation debt, or missing feature |
| Evidence | Exact files, symbols, tests, commands, and observed behavior |
| Authority | Contract or documented boundary used for judgment |
| Consequence | Current failure or credible Chapter 5/6 change risk |
| Scope relevance | Common prerequisite, target-dependent, deferred, or out of scope |
| Action tier | Before any feature, before affected feature, opportunistic, or do not act |
| Recommendation | Smallest coherent change that addresses the evidence |
| Verification | Tests, comparisons, and documentation needed to accept it |
| Tradeoff | Complexity introduced, alternatives, and reason for the choice |

For positive findings, record the evidence, the invariant worth preserving,
and the kinds of change that could accidentally weaken it.

Do not use a simple “critical/high/medium/low” list without timing and scope.
A target-dependent issue can be mathematically serious while still being
irrelevant to the first selected application.

## Required assessment structure

The final assessment document should contain:

1. Executive verdict.
2. Audited baseline and limitations.
3. Project scope and Chapter 5/6 traceability matrix.
4. Confirmed strengths and invariants to preserve.
5. Architecture and code findings by layer.
6. Test and numerical-verification findings.
7. Build and tooling findings.
8. Documentation and agent-guidance findings.
9. Refactor decision table with action tiers.
10. Recommended dependency-ordered refactor batches.
11. Deferred work and explicit non-goals.
12. Open decisions requiring user selection.
13. Verification checklist for Stage B.

The executive verdict must answer one of these forms, with qualifications:

- no broad refactor is warranted;
- a small cross-cutting refactor is required before further features;
- refactor only before particular selected Chapter 5/6 capabilities; or
- a broader structural refactor is required before the documented endpoint is
  safely achievable.

## Decision rules for recommending a refactor

Recommend refactoring when evidence shows at least one of the following:

- a documented layer boundary is currently bypassed;
- the same mathematical policy is duplicated and already diverging;
- adding a documented capability would require editing multiple complete
  problem implementations rather than registering a component;
- a public contract permits invalid states that tests or diagnostics cannot
  reliably contain;
- correctness cannot be tested independently at the responsible layer;
- build or test structure conceals failures or makes backend-neutral work
  impractical;
- stale or conflicting authority is likely to direct future implementation
  incorrectly;
- a focused change provides a concrete software-engineering lesson without
  speculative abstraction.

Prefer preserving the current design when:

- a narrow implementation accurately records its supported policy;
- duplication represents genuinely different mathematical semantics;
- an unsupported combination fails with a precise capability diagnostic;
- extraction would create an interface broader than the selected targets;
- the only benefit is aesthetic uniformity;
- tests already isolate the behavior and likely future changes do not cross
  the boundary.

## Stage B execution policy

After user approval, convert accepted findings into dependency-ordered batches.
For each batch:

1. Restate the finding IDs and affected layer.
2. Read the authoritative contract and relevant conventions.
3. Define behavior that must remain unchanged and any deliberate contract
   change.
4. Add or improve focused tests before or with the refactor.
5. Make the smallest coherent structural change.
6. Run focused tests, then the full baseline configuration.
7. Compare against `pre-refactor-ch5-ch6` where numerical or API behavior
   could drift.
8. Update authoritative documentation and the assessment status.
9. Inspect the diff and leave the change pending for user review.

Suggested batch boundaries, to be accepted or rejected by evidence, are:

```text
R0  correctness and stale-authority repairs
R1  test organization and reusable verification helpers
R2  compiler responsibility and diagnostic cleanup
R3  shared deal.II assembly/solve services with preserved policies
R4  formulation and solver ports required by selected Chapter 6 methods
R5  build and automation improvements justified by the preceding work
```

These are planning placeholders, not predetermined recommendations. Delete or
reorder them when the assessment evidence does not support them.

## Exit criteria for Stage A

Stage A is complete only when:

- the baseline build/test state is recorded or an explicit blocker is
  documented;
- all tracked C++ headers, tests, and CMake configuration have been reviewed;
- authoritative and status-bearing documentation has been checked against the
  implementation;
- contributor and agent instructions have been reviewed;
- every P5/P6 capability family appears in the traceability matrix;
- strengths and improvement findings both cite concrete evidence;
- missing features are not mislabeled as refactor debt;
- every proposed refactor has scope relevance, action timing, tradeoffs, and
  verification;
- deferred and out-of-scope work is explicit;
- the user can decide what to refactor without another discovery pass;
- no code, test, CMake, or instruction file was changed during the assessment.

For the documentation-only Stage A change, run `git diff --check`, inspect the
rendered-looking Markdown, and verify new relative links before handoff.

## Immediate next-agent handoff

The next agent should perform only Stage A unless the user explicitly approves
implementation. Its first actions are:

1. Verify that the worktree is clean or identify existing user changes and
   confirm branch `codex/ch5-ch6-development`.
2. Follow the mandatory convention routing.
3. Reproduce and record the baseline.
4. Create `docs/planning/review/pre-ch5-ch6/assessment.md`.
5. Build the scope matrix before deep code review.
6. Audit progressively in the layer order given above.
7. Finish the evidence-backed report and request the user’s Stage B choices.

Do not begin Chapter 5/6 feature implementation and do not implement refactors
while writing the assessment.
