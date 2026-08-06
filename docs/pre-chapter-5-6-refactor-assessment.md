# Pre-Chapter 5/6 refactor assessment

## Assessment status

**In progress — evidence wave 2 complete.**

This document records the assessment defined by the
[refactor assessment plan](pre-chapter-5-6-refactor-plan.md). The current
contents establish the reproducible baseline, the Chapter 5/6 scope matrix,
and the backend-neutral contract and semantic-layer findings. Compiler,
deal.II realization, solver, broader test-quality, build-design, and
documentation findings remain pending until their corresponding review waves
are complete.

No code, tests, CMake, conventions, or agent instructions may be changed while
this assessment is in progress.

## Executive verdict

**Preliminary outcome: a small cross-cutting refactor is warranted before new
semantic features are added.** The build and test baseline is healthy, and the
core mathematical distinctions are good. However, the public block wrapper
can be mutated out of agreement with its layout, and semantic validation
currently accepts several malformed graphs. Reference graph factories also
modify inherited graphs by positional vector index. These are scope-neutral
engineering issues whose cost and risk increase as Chapter 5 adds spaces,
pairings, terms, and observations.

It is still too early to decide whether the compiler and deal.II realization
need a broader structural refactor. The final verdict must distinguish the
small common prerequisite above from target-dependent work such as mixed
blocks, KKT products, and complementarity.

The final verdict will select and qualify one of these outcomes:

- no broad refactor is warranted;
- a small cross-cutting refactor is required before further features;
- refactor only before particular selected Chapter 5/6 capabilities; or
- a broader structural refactor is required before the documented endpoint is
  safely achievable.

## Audited baseline

### Repository state

The baseline was recorded on 2026-08-06.

| Item | Observed state |
| --- | --- |
| Assessment branch | `codex/refactor-ch5-ch6-readiness` |
| Assessment start commit | `8616cba` (`docs(refactor): add pre-chapter 5 and 6 assessment plan`) |
| Tagged implementation baseline | `pre-refactor-ch5-ch6` at `7c2496b` |
| Baseline relation | Assessment branch is one documentation-only commit ahead of the tagged implementation baseline |
| Worktree at start | Clean |
| Generated build directories used | `build-audit/` and `build-nodealii/`, both ignored by `/build*/` |

The tag is the immutable comparison point for later behavior-preserving
refactors. The assessment branch is the only integration branch authorized by
the plan.

### Toolchain and dependency state

| Tool or dependency | Observed state |
| --- | --- |
| C++ compiler | GCC 13.3.0 (`/usr/bin/c++`) |
| CMake | 3.28.3 |
| Ninja | 1.11.1 |
| deal.II | 9.5.1, discovered from `/usr/share/cmake/deal.II` |
| deal.II configuration | Installed `DebugRelease`; project configuration was forced to `Debug` when no build type was supplied |

`pkg-config` does not advertise deal.II, but this is not a project failure:
the supported CMake package discovery succeeded.

### Configuration observations

The exact convention baseline command could not reuse `build/`:

```text
CMake Error: generator Ninja does not match the generator used previously:
Unix Makefiles
```

The existing `build/CMakeCache.txt` selects Unix Makefiles, while
`build-ninja/CMakeCache.txt` selects Ninja. No existing build tree was deleted
or rewritten. A fresh ignored `build-audit/` directory was used to obtain
independent compilation evidence.

Fresh deal.II configuration succeeded but emitted this policy warning:

```text
CMAKE_BUILD_TYPE "" unsupported by current installation.
deal.II was configured with "DebugRelease".
CMAKE_BUILD_TYPE was forced to "Debug".
```

These are baseline observations, not yet classified refactor findings. The
build audit must decide whether the repository should prevent, document, or
accept them.

### Fresh build and test results

#### deal.II-enabled configuration

Commands:

```bash
cmake -S . -B build-audit -G Ninja -DBUILD_TESTING=ON
cmake --build build-audit
ctest --test-dir build-audit --output-on-failure
```

Results:

| Test | Result | Wall time |
| --- | --- | --- |
| `nmopt_contract_tests` | Passed | 0.05 s |
| `nmopt_semantic_v1_contract_test` | Passed | 0.10 s |
| `nmopt_dealii_diffusion_contract_test` | Passed | 9.99 s |
| Overall | 3/3 passed | 10.14 s |

The fresh build compiled and linked all three executables. No compiler warning
was printed during this build.

#### Backend-neutral configuration

Commands:

```bash
cmake -S . -B build-nodealii -G Ninja \
  -DBUILD_TESTING=ON -DNMOPT_ENABLE_DEAL_II=OFF
cmake --build build-nodealii
ctest --test-dir build-nodealii --output-on-failure
```

Results:

| Test | Result | Wall time |
| --- | --- | --- |
| `nmopt_contract_tests` | Passed | 0.01 s |
| `nmopt_semantic_v1_contract_test` | Passed | 0.01 s |
| Overall | 2/2 passed | 0.01 s |

This establishes that the typed contracts, dense reference path, semantic
graph, and semantic validation can be built and tested without deal.II.

### CTest registration baseline

The deal.II-enabled configuration registers three coarse executables and the
backend-neutral configuration registers two. CTest assigns no timeout,
labels, fixtures, resource locks, or other properties beyond the working
directory. Whether that granularity is appropriate remains a Phase 4 question.

### Baseline limitations

- Tests were run once in each fresh configuration; nondeterminism and repeated
  execution have not yet been assessed.
- The passing suite establishes current behavior, not completeness or
  independence of mathematical verification.
- No warning-enabling configuration, sanitizer, static analysis, coverage, or
  performance measurement was introduced.
- The optional fallback deal.II discovery path was not exercised because the
  CMake package was available.
- No Chapter 6 reproduction target was run; the examples guide explicitly
  classifies those experiments as deferred benchmarks rather than current
  correctness tests.

## Scope and traceability matrix

### Interpretation

The matrix is derived from P5.1–P5.6 and P6.1–P6.5 of the
[implementation roadmap](implementation-roadmap.md), the
[Chapter 5 implementation guide](chapter-5-implementation-guide.md), and the
[Chapter 6 numerical-methods guide](chapter-6-numerical-methods-guide.md).

“Documented baseline” reports what those sources claim. It is not treated as
verified code evidence until the relevant implementation and tests are
reviewed. “Assessment exposure” states the question to investigate; it is not
a finding.

No capability below is an unconditional implementation commitment. The final
report must preserve target-dependent recommendations until the user selects
which catalogue problems and methods to implement.

### Chapter 5 capability families

| Capability | Authority and target condition | Documented baseline | Necessity | Assessment exposure |
| --- | --- | --- | --- | --- |
| Scalar elliptic volume terms: tensor diffusion, transport, reaction, source, volume control | P5.1; needed by selected C5.1, C5.5.1, C5.6 and the P6.2 advection target | Scalar constant diffusion-reaction and volume control exist; general coefficients and transport are absent | Target-dependent | Determine whether physics can be registered term by term or currently requires another complete assembled target |
| Robin bilinear/source terms and boundary partition | P5.1; Robin applications | Fixed/controlled Dirichlet and Neumann slices exist; Robin and general overlap policy are absent | Target-dependent | Audit whether boundary policies share explicit composition points without conflating essential, natural, Robin, and transport boundaries |
| Energy-state observation | P5.2; C5.5.2 | Volume and boundary tracking exist; state $H^{1}$ observation is absent | Target-dependent | Check whether observations and losses are independently lowerable or embedded in specialized targets |
| Weighted boundary-trace observation | P5.2; C5.7 variants | Unweighted boundary trace is documented | Target-dependent | Check whether fixed weight data can enter an observation without changing a loss or residual implementation |
| $H^{-1}$ control metric | P5.2; alternate C5.5.2 control geometry | $L^{2}$ and $H^{1}$ metrics exist; exact $H^{-1}$ policy is absent | Target-dependent and separate from observations | Verify metric ports and solver diagnostics can express the selected Riesz/inverse policy without changing objective semantics |
| Point-set regions and sensor observation | P5.3; C5.10 | Unsupported | Target-dependent | Determine whether region, observation-space, data, and pairing declarations can be extended without special application branches |
| Normal-flux observation | P5.3; C5.8 | Unsupported | Target-dependent | Audit trace/regularity policy boundaries and avoid treating flux as an ordinary boundary restriction |
| Strong/transposition/very-weak formulation policy | P5.3; low-regularity C5.8, C5.10, and C5.11 variants | Unsupported | Prerequisite only for selected low-regularity targets | Determine whether formulation ownership belongs outside the current reduced DTO and how diagnostics prevent undeclared approximations |
| Partial controlled-Dirichlet reconstruction with corner/interface policy | P5.4; general C5.11 and later boundary Stokes control | Complete-boundary nodal lifting exists; partial boundaries and corner/interface rules are rejected | Target-dependent | Audit whether the existing physical reconstruction has extractable value/JVP/VJP actions or is coupled to one compiled target |
| Trace metrics and surface-gradient maps | P5.4; fractional or tangential boundary geometry | Boundary $L^{2}$ metric exists; fractional and tangential choices are absent | Target-dependent follow-up | Preserve separation of regularisation, observation, and search metric; do not generalize before a selected trace space exists |
| Regularised state-control constraint observation and multipliers | P5.5; regularised C5.12 | Control-only box projection exists; state-derived constraint, multipliers, and complementarity are absent | Target-dependent; coupled to P6.3/P6.5 | Audit constraint interface limits and avoid misrepresenting unregularised measure-multiplier state constraints |
| Multiple state/equation blocks and coupled reduced DTO | P5.6; all Stokes targets | Current executable DTO is one scalar state, one decision, and one test block | Required only if a Stokes target is selected | Determine whether block generalization belongs in contracts/formulations and whether scalar behavior can remain simple |
| Vector mixed residual terms and Stokes policies | P5.6; C5.13 | Unsupported | Required only if a Stokes target is selected | Audit whether semantic/compiler registries can add vector, pressure, gauge, and inf-sup policies compositionally |

### Chapter 6 capability families

| Capability | Authority and target condition | Documented baseline | Necessity | Assessment exposure |
| --- | --- | --- | --- | --- |
| Typed reduced search-direction policies | P6.1; nonlinear CG and L-BFGS targets | Metric-qualified steepest descent exists | Target-dependent; likely first reduced-method extension | Check whether direction choice can be separated from PDE evaluation, metric application, and reporting without a universal dispatcher |
| Line-search policies and uniform solver report | P6.1; exact quadratic, Armijo, or Wolfe choices | Armijo is embedded in the current reduced-gradient solver with existing diagnostics | Target-dependent | Audit ownership of trial evaluation, projected displacement, counters, curvature failures, and stopping reasons |
| Formulation/trial/test provenance | P6.2; every DTO/OTD/stabilized comparison | DTO convention and manifests exist; general formulation policy is absent | Required for selected formulation comparisons | Check whether manifest and compiler contracts can record provenance without claiming equivalence |
| Stabilization components and OTD product | P6.2; Graetz/advection and stabilized targets | Unsupported | Target-dependent | Determine whether exact stabilized derivatives can be owned as components and whether OTD needs a separate executable product |
| Equality-constrained quadratic KKT product | P6.3; all-at-once examples and prerequisite for P6.4/P6.5 | Reduced DTO only; generic KKT action absent | Required for all-at-once or PDAS targets | Audit whether existing JVP/VJP and objective actions expose enough structure without broadening `ExecutableModelT` incorrectly |
| KKT Krylov policies and multiplier provenance | P6.3 | Current state/adjoint solves are target-owned; no general KKT solver | Required with P6.3 | Check layout, sign conversion, rank/kernel policy, and solver diagnostic ownership |
| Composable approximate solves and block preconditioners | P6.4; all-at-once performance methods | Unsupported | Requires P6.3; Stokes variants also require P5.6 | Audit whether approximate inverse services can remain deterministic operator ports rather than backend/problem branches |
| Typed complementarity and active selection | P6.5; control-box and regularised mixed constraints | Projection exists without multipliers or active-set services | Required only for PDAS targets | Audit primal/dual representation and conversion boundaries, especially continuous versus cellwise classification |
| PDAS/semismooth-Newton orchestration | P6.5; Section 6.9 examples and P5.5 regularised constraint | Unsupported | Requires P6.3 and complementarity services | Determine whether repeated KKT subproblems can reuse formulation services and report full KKT convergence independently |

### Dependency summary

```text
P5.1 scalar terms ────────────────> P6.2 advection/stabilisation target

P5.5 regularised constraint ──┐
                              ├──> P6.5 complementarity/PDAS
P6.3 KKT product ─────────────┘
        │
        └─────────────────────────> P6.4 block preconditioners

P5.6 mixed blocks ────────────────> Stokes portions of P6.1/P6.3/P6.4/P6.5

P5.3 transposition policy ────────> low-regularity flux/sensor/Dirichlet targets
```

P6.1 can be developed against the current scalar reduced DTO without waiting
for mixed blocks. P6.3 can likewise begin with the documented scalar
linear-quadratic target. P5.6 is therefore not a prerequisite for all Chapter
6 work, only for its Stokes variants.

### Scope conclusions to carry into the code audit

1. The first code-review question is not “Can the framework implement every
   row?” It is “Does current ownership force unrelated rows to change
   together?”
2. Specialized implementations are acceptable when they enforce a selected
   discrete or functional-analytic policy and reject alternatives clearly.
3. Repeated whole-problem assembly is a refactor candidate only when the same
   policy is copied, diverging, or blocking independent component registration.
4. Mixed blocks, very-weak formulations, KKT products, and complementarity
   are explicit new contracts, not small generalizations to infer during a
   refactor.
5. Metrics, losses, observations, residual terms, and formulations must remain
   separate even when one textbook example selects them together.

## Confirmed strengths and invariants to preserve

The following strengths are established by the baseline and scope documents;
code-level confirmation remains part of later waves.

### S-001 — Backend-neutral contracts have an executable build boundary

**Evidence:** Configuring with `NMOPT_ENABLE_DEAL_II=OFF` builds and passes the
contract and semantic test executables without deal.II.

**Preserve:** Core contract, dense reference, and semantic headers must not
acquire deal.II dependencies during refactoring.

### S-002 — The current checked-in baseline passes a fresh deal.II build

**Evidence:** A clean generated Ninja tree compiled all three test executables
and passed all registered tests with deal.II 9.5.1.

**Preserve:** Refactor batches must compare against this baseline rather than
changing tests or configuration to hide regression.

### S-003 — Scope and mathematical distinctions are explicit

**Evidence:** The P5/P6 roadmap separates residuals, observations, losses,
metrics, constraints, formulations, and solver services; it also records
unsupported and target-dependent cases.

**Preserve:** Refactoring must not erase differences such as Neumann versus
Dirichlet control, reduced covector versus metric direction, DTO versus OTD,
or regularised versus unregularised state constraints.

### S-004 — Primal and dual coefficients are distinct public C++ types

**Evidence:** `PrimalBlockT` and `CovectorBlockT` are separate final wrapper
types in `include/nmopt/contract/layout.hpp:154-183`. The only primitive
pairing has the typed signature `pair(const CovectorBlockT &,
const PrimalBlockT &)` at lines 198-210. `ExecutableModelT` returns covectors
from residual and derivative actions and accepts a primal test-space seed for
the residual VJP.

**Preserve:** Do not replace these wrappers with untyped backend vectors or
add an implicit primal/covector conversion. Any Riesz identification must
continue to pass through a declared metric.

### S-005 — The reduced DTO boundary is narrow and mathematically honest

**Evidence:** `StateControlPartitionT` rejects anything other than two
variable blocks and one test block in
`include/nmopt/contract/reduced_dto.hpp:18-37`. `ReducedDTOT::evaluate()`
centralizes the documented Lagrangian sign by computing
`J_u' - E_u'^* p` at lines 147-163. State and adjoint solves are injected as
formulation services rather than hidden in a residual term.

**Preserve:** Mixed states, multiple equations, all-at-once actions, OTD, and
second derivatives must be new explicit formulation capabilities. They must
not be smuggled into the current binary DTO by weakening its checks.

### S-006 — The semantic graph is backend-free and diagnostics have explicit categories

**Evidence:** `include/nmopt/semantic/v1/types.hpp` includes only standard
strings and vectors and declares regions, spaces, pairings, variables, data,
transformations, terms, equations, observations, losses, metrics,
constraints, policies, and a formulation as data. Backend objects are bound
later. `DiagnosticCategory` distinguishes structural, analytical-policy,
lowerability, and formulation-capability failures.

**Preserve:** The semantic layer must remain usable in the deal.II-disabled
build, and unsupported choices must remain diagnostics rather than silent
compiler substitutions.

### S-007 — The dense reference exercises the core derivative identities

**Evidence:** `tests/reduced_dto_contract.cc:118-220` checks residual
JVP/VJP pairing, a residual finite difference, the objective directional
derivative, the state residual, the reduced derivative, metric
apply/inverse-apply consistency, and box projection. The same executable
instantiates the basic algebra with an alternate backend type at lines
292-321.

**Preserve:** Keep a small, fast, deal.II-free algebraic path. Later compiler
tests should complement it, not make deal.II necessary to verify DTO signs or
typed algebra.

## Architecture and code findings: contracts and semantic layer

### RF-001 — Block storage can be mutated out of agreement with its layout

**Classification:** Correctness defect in a public contract.

**Evidence:** `BlockValuesT` validates block count and vector dimensions only
in its constructor (`include/nmopt/contract/layout.hpp:124-136`) but exposes a
mutable backend vector reference from `block()` at lines 112-116.
`require_compatible()` later compares only the two immutable layouts at lines
189-196. A throwaway C++17 probe against the current headers created primal
and covector blocks with declared dimension two, assigned three-entry dense
vectors through the mutable references, and called `pair()`. The call
succeeded and returned `32` while reporting declared dimension `2` and stored
dimension `3`.

**Authority:** The [executable contract](executable-contract-v0.md) describes
`BlockLayout` as the compatible block dimensions, and the
[implementation-readiness review](implementation-readiness-review.md)
requires typed coefficient objects whose discrete layouts are preserved.

**Consequence:** Layout compatibility is not an invariant after construction.
Two equally corrupted objects can pass `require_compatible()` and produce an
algebraic result under a false space declaration. This becomes more dangerous
with deal.II vectors that expose `reinit()`, mixed layouts, KKT blocks, and
ownership/ghost state.

**Scope relevance:** Common prerequisite. The immediate code uses
well-behaved callers, but P5.6 mixed blocks and P6.3 KKT products would amplify
the invalid-state surface.

**Action tier:** Before any feature. This is a small core-contract refactor,
not a Chapter 5/6 feature.

**Recommendation:** Make dimension-changing storage mutation unavailable from
the public wrapper. Prefer const public block access plus checked algebraic
updates and a narrowly scoped internal mutation path for construction and
backend assembly. If that is disproportionate for the current backend
surface, the minimum acceptable alternative is a `validate_storage()` check
at every public algebra/model boundary and after any controlled mutable
operation. Do not rely on caller discipline alone.

**Verification:** Add negative contract tests that attempt to violate storage
dimensions through every retained mutable path; require a `ContractError`
before pairing or model evaluation. Re-run dense and deal.II derivative,
metric, constraint, DTO, and solver tests without changing numerical results.

**Tradeoff:** Restricting mutation requires a few checked update helpers and
may make backend assembly more explicit. Rechecking dimensions is simpler but
continues to represent invalid states temporarily and adds runtime checks.

### RF-002 — Semantic aggregates have no safe incomplete state

**Classification:** Correctness defect in public semantic construction.

**Evidence:** The enum-bearing aggregate fields in
`include/nmopt/semantic/v1/types.hpp:133-280` have no default member
initializers, and the enums at lines 11-131 have no `unspecified` value.
Default-initializing a component and filling it incrementally can therefore
leave fields such as `RegionSpec::kind`, `SpaceSpec::topology`, or
`ResidualTermSpec::kind` indeterminate. `SemanticValidator` reads these fields
in comparisons and switches, for example at
`include/nmopt/semantic/v1/validation.hpp:166-191` and lines 711-745. The
reference factories avoid the fault by fully aggregate-initializing every
component, but the public API does not enforce that discipline.

**Authority:** The [interface specification](interface-specification.md)
defines semantic resolution and validation as the boundary that rejects
invalid or incomplete declarations. Validation itself must therefore be safe
on incomplete user input.

**Consequence:** A partially assembled public graph can reach validation with
an indeterminate enum value rather than a diagnosable `unspecified` value.
Adding more Chapter 5 node kinds increases the number of such construction
sites and switch statements.

**Scope relevance:** Common prerequisite for semantic feature work.

**Action tier:** Before any feature that adds or constructs semantic nodes.

**Recommendation:** Add an explicit invalid/unspecified enumerator and default
member initializer for every enum-bearing field, then diagnose it
structurally. This retains aggregate ergonomics and is smaller than adding a
constructor hierarchy. Constructors or builders that require every field are
an alternative only if they remain plain composition utilities.

**Verification:** Validate default and partially populated instances of every
semantic component under normal, warning-enabled, and sanitizer test
configurations. Require exact structural diagnostics and no compiler-only
failure.

**Tradeoff:** Sentinel values add switch cases and require compiler registries
to reject them. That explicit work is preferable to making an arbitrary real
kind the default or relying on initialization order.

### RF-003 — Structural validation does not establish whole-graph closure

**Classification:** Correctness defect and architectural debt.

**Evidence:** `validate_equations()` verifies each term named by an equation,
and `validate_terms()` verifies that each term names an existing equation, but
neither checks that every declared term appears exactly once in its target
equation (`include/nmopt/semantic/v1/validation.hpp:353-455`). Constraint
validation checks bound roles and kinds but not that bound data use the
constrained variable's space at lines 606-668. No validation checks required
human-readable labels, even though the common component rules require them.

A throwaway C++17 probe copied the canonical valid graph and applied one
mutation at a time. `SemanticValidator::validate()` reported each of the
following graphs as valid:

| Accepted malformed graph | Mutation |
| --- | --- |
| Orphan residual term | Added a valid diffusion-reaction term that names `state_equation` but is absent from that equation's term list |
| Duplicate equation edge | Added `volume_source` to its owning equation a second time |
| Bound in the wrong space | Changed the lower-bound datum from `control_space` to `state_space` |
| Missing required label | Cleared the state variable's label |

Source inspection also shows that formulation validation checks only that its
metric and optional constraint IDs exist; it does not establish that they act
on the formulation's selected decision variable
(`include/nmopt/semantic/v1/validation.hpp:778-816`).

**Authority:** Sections 2, 5.2, and 5.3 of the
[interface specification](interface-specification.md) require declared ports,
ordered semantic resolution, and structural rejection of mismatched
connections. An equation owns an ordered sum of terms; a constraint and
metric act on declared variable blocks.

**Consequence:** A valid semantic report is currently weaker than a resolved
typed graph. A later compiler may reject some malformed graphs as
unlowerable, but that misclassifies an earlier structural error and lets
backend-specific exact-shape checks substitute for semantic integrity.
Duplicate edges may also change an accumulated residual rather than fail if a
future lowerer composes terms generically.

**Scope relevance:** Common prerequisite for P5.1 term additions, P5.2/P5.3
observations, P5.5 constraints, and any later multi-variable formulation.

**Action tier:** Before any new semantic component.

**Recommendation:** Add one graph-level structural pass after local indexes
are built. It should enforce bidirectional term ownership, reject duplicate
edges, validate label requirements, match constraint data spaces, and match
the selected metric/constraint to the formulation decision variable. Keep
lowerability and formulation capability checks in the compiler. Do not reject
otherwise harmless unused declarations unless the architecture first makes
that a documented rule.

**Verification:** Turn every accepted mutation above into an exact negative
test. Add formulation metric/constraint mismatch cases. Each test must assert
category, component ID, and capability, not only that some structural error
exists.

**Tradeoff:** Stricter validation may reveal malformed graphs that the current
exact compiler registry already rejects for a different reason. That is a
useful compatibility break at the diagnostic boundary; the canonical graphs
and compiled numerical behavior should remain unchanged.

### RF-004 — Pairing validation ignores the declared covector side

**Classification:** Correctness defect in typed semantic validation.

**Evidence:** `validate_pairings()` checks only that both referenced space IDs
exist (`include/nmopt/semantic/v1/validation.hpp:215-227`). Equation,
observation, loss, and metric checks compare only
`PairingSpec::primal_space_id` with the consuming port at lines 373-379,
522-528, 551-559, and 590-596. A throwaway probe changed the canonical state
test pairing's `covector_space_id` to `control_space`; the semantic report
remained valid.

**Authority:** The [interface specification](interface-specification.md)
requires every dual value to be paired with its declared primal space, and the
[implementation-readiness review](implementation-readiness-review.md) makes
the primal/covector distinction the first discrete-algebra default.

**Consequence:** The semantic graph can claim a test-space pairing while its
declared covector side names an unrelated control space. The current compiler
may accept only its exact reference shapes, but new observation, trace,
mixed-field, and Petrov-Galerkin work would inherit a false typed boundary.

**Scope relevance:** Before affected features. This applies to nearly every
Chapter 5 capability that introduces a space or pairing, but it does not
require implementing any of those capabilities now.

**Action tier:** Before the next new pairing or space realization.

**Recommendation:** Document the precise v1 compatibility rule for a
pairing's covector port and enforce both sides at semantic consumers. For the
current dual-coefficient slice this will likely require the covector
declaration to identify the dual of the same semantic space; if distinct
space descriptors are intended, compatibility must compare region, shape,
topology, and declared dual relation rather than merely string existence.

**Verification:** Add negative pairings for wrong region, role, topology, and
field shape as those notions become available. Retain residual VJP pairing
tests after compilation so structural compatibility and numerical transpose
consistency remain separate checks.

**Tradeoff:** An exact-ID v1 rule is simple but must not accidentally forbid a
future explicitly declared nontrivial pairing. Writing the rule down before
coding avoids over-generalizing the current space model.

### RF-005 — Reference graph variants depend on positional mutation

**Classification:** Architectural and code-quality debt.

**Evidence:** `include/nmopt/semantic/v1/reference_specs.hpp` contains 36
numeric `.at(index)` mutations, plus positional erase operations such as
`residual_terms.erase(begin() + 2)`. The H1, coefficient-identification,
fixed-Dirichlet, and controlled-Dirichlet factories derive variants by
rewriting assumed positions in the base scalar graph. For example, the
coefficient factory rewrites spaces, pairings, the second variable, the first
term, the second observation/loss, and the first metric at lines 202-265.

**Authority:** The project mission and
[composition boundaries](composition-boundaries.md) require combinations to
be built from connected components without a new whole-problem
implementation for each combination.

**Consequence:** Inserting or reordering an otherwise independent component
in the base graph can make a variant mutate the wrong component while all
indices remain in range. The named factories are useful examples and are not
an inheritance hierarchy, but their implementation couples independent
features to one canonical vector order. Chapter 5 additions would multiply
that coupling.

**Scope relevance:** Common prerequisite for maintainable semantic examples;
it does not imply a general user-facing problem-description language.

**Action tier:** Before the first Chapter 5 semantic graph addition.

**Recommendation:** Introduce small ID-based lookup/replace/remove helpers
that require exactly one match, and express each existing variant as a named
feature delta over shared declarations. Keep the current named reference
factories as readable test fixtures and public examples. Avoid a broad fluent
builder or PDE-family dispatcher unless a concrete later use requires it.

**Verification:** Rebuild every current reference graph, require it to validate
and compile, compare semantic identities/manifests, and rerun existing v0/v1
and derivative comparisons. Add a unit test showing that harmless base-vector
reordering does not redirect a feature delta.

**Tradeoff:** ID lookup is more verbose and detects failures at runtime during
graph construction rather than at compile time. It is still substantially
safer than positional mutation and much smaller than a new graph DSL.

### RF-006 — Contract and semantic negative tests are too coarse for refactoring

**Classification:** Test debt.

**Evidence:** `tests/semantic_v1_contract.cc` has 12 negative assertions using
only `ValidationReport::has_category()`. It never asserts the available
`Diagnostic::component_id` or `Diagnostic::capability`, so an unrelated error
in the same category can satisfy a test. `tests/reduced_dto_contract.cc`
checks many successful numerical identities but has no expected-failure case
for layout mismatch, invalid partition, solver callback, metric, or
constraint construction. All current tests pass, so this is a containment and
localization weakness rather than a failing baseline.

**Authority:** The [implementation-readiness review](implementation-readiness-review.md)
requires explicit component diagnostics and value/JVP/VJP, constraint, and
formulation verification. The assessment plan requires tests to localize the
responsible layer and cover negative and degenerate cases.

**Consequence:** A refactor can change which component failed, or can leave
new invalid states accepted, while the current semantic test remains green.
The block-layout defect in RF-001 had no regression test capable of exposing
it.

**Scope relevance:** Common Stage B prerequisite.

**Action tier:** Before any implementation refactor. Add characterization
tests first, then change production code.

**Recommendation:** Add a reusable exact-diagnostic matcher over category,
component ID, and capability; make negative semantic cases table-driven where
that improves readability. Add `expect_contract_error` tests for public
contract preconditions and separate the dense test's algebra, DTO, metric,
constraint, and solver cases into focused functions. CTest executable
granularity can be decided in the broader test/build wave.

**Verification:** Demonstrate that each test fails if its expected capability
or component is changed. Run the backend-neutral suite first, then the full
deal.II suite, and preserve all current numerical checks.

**Tradeoff:** Exact assertions require intentional updates when diagnostic
names change. That maintenance cost is the mechanism that prevents silent
diagnostic drift.

## Findings and observations pending cross-layer review

The following baseline observations remain queued for later classification:

| Observation | Evidence needed before classification |
| --- | --- |
| `OBS-BLD-01`: the prescribed `build/` Ninja command conflicts with an existing Makefiles cache | Determine whether this is only local stale state or a documentation/workflow weakness reproducible for contributors |
| `OBS-BLD-02`: an unspecified build type is silently forced to `Debug` by deal.II after a warning | Review CMake policy, expected developer configurations, and whether performance-sensitive examples need an explicit choice |
| `OBS-TST-01`: all deal.II behavior is registered as one CTest executable | Review internal test organization, failure handling, runtime, and whether granularity impedes diagnosis or focused refactoring |
| `OBS-ARC-01`: `ReducedDTOT` and injected solve functions use non-owning references | Trace ownership through compiled products and solvers before deciding whether the lifetime contract is safe, merely undocumented, or defective |
| `OBS-ARC-02`: cellwise projection capability is matched by layout and a metric string identifier | Check whether compiler construction makes spoofing impossible in practice and how additional metrics are expected to advertise projection compatibility |
| `OBS-ARC-03`: the semantic validator is a 955-line public header with 70 diagnostic sites and centralized kind switches | Compare this ownership with compiler registry structure before deciding whether a focused internal split is warranted |

## Pending review waves

1. Compiler ownership, manifests, capability diagnostics, v0/v1 independence,
   and deal.II realizations.
2. Formulations, solvers, tests, numerical verification, CMake, and tooling.
3. Documentation consistency, conventions, and agent guidance.
4. Final synthesis, action tiers, dependency-ordered refactor batches,
   deferred work, and user decisions.
