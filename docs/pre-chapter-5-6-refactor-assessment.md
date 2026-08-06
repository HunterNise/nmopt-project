# Pre-Chapter 5/6 refactor assessment

## Assessment status

**In progress — evidence wave 4 complete.**

This document records the assessment defined by the
[refactor assessment plan](pre-chapter-5-6-refactor-plan.md). The current
contents establish the reproducible baseline, the Chapter 5/6 scope matrix,
and findings for the backend-neutral contracts, semantic layer, v1 compiler,
compiled-product ownership, manifests, current deal.II realizations, reduced
formulation and solver, tests, numerical verification, CMake, and developer
tooling. Documentation, conventions, agent guidance, and the final synthesis
remain pending until their corresponding review waves are complete.

No code, tests, CMake, conventions, or agent instructions may be changed while
this assessment is in progress.

## Executive verdict

**Preliminary outcome: two bounded refactor tiers are warranted, not a
rewrite.** A small cross-cutting safety refactor is required before any new
semantic feature: restore block/layout invariants, make incomplete semantic
objects safe to validate, close structural graph validation, replace
positional reference-graph mutation, and strengthen negative tests.

A focused compiler decomposition is also warranted before adding
component-heavy Chapter 5 features or a Chapter 6 DTO/OTD or stabilization
comparison. The advertised component registry is currently a boolean
whitelist in front of an exact whole-graph matcher and a six-way whole-model
dispatcher. The compiled targets contain valuable, tested realization
policies, but common scalar assembly and solve policy are repeated across
them. The right refactor is therefore to introduce a resolved lowering plan
and extract shared realization services while preserving specialized
reconstruction, trace, nullspace, and parameter-dependent policies as
composed strategies.

This compiler tier is target-dependent: an isolated P6.1 reduced search or
line-search improvement should not be blocked on it. P6.1 does, however,
need a smaller solver-local refactor: separate value-only trial evaluation
from derivative evaluation, make direction and line-search choices modest
policies instead of branches in the current driver, and produce a typed report
whose work counts and accepted steps can be audited. The current binary DTO
boundary should be preserved; OTD and KKT products remain separate future
formulations.

Test and build changes are supporting batches rather than a third architecture
tier. The mathematical test coverage is unusually strong for the repository's
size, but three process-level tests hide eleven logical scenarios, exact
solver failure behavior is largely uncharacterized, and the optional deal.II
configuration can degrade silently. No current evidence justifies a runtime
plugin system, a general graph DSL, a PDE inheritance hierarchy, packaging
machinery, or implementation of unselected Chapter 5/6 capabilities. The
final verdict must still integrate the documentation and guidance wave.

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
| Baseline relation | At the start of evidence wave 4, the assessment branch is four documentation-only commits ahead of the tagged implementation baseline |
| Worktree at start | Clean |
| Generated build directories used | `build-audit/`, `build-nodealii/`, `build-audit-warnings/`, `build-audit-dealii-warnings/`, and `build-audit-sanitizers/`, all ignored by `/build*/` |

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

### Supplemental diagnostics and limitations

- `ctest --repeat until-fail:3` passed all three registered tests on every
  repetition. The first deal.II process took 22.17 s after the rebuild and the
  next two took 0.32 s and 0.30 s, so this establishes repeatability of
  outcomes, not a portable runtime baseline.
- A backend-neutral Debug build with `-Wall -Wextra -Wpedantic -Wconversion
  -Wshadow` compiled without warnings and passed both tests.
- The same warning profile with deal.II compiled and passed all tests, but
  exposed pervasive narrowing warnings at the adapter boundary from
  `std::size_t` to the 32-bit `dealii::Vector<double>::size_type`. The clean
  three-target build took 54.00 s and 1,486,284 KiB peak resident memory on
  this machine. These are local build-cost observations, not benchmark gates.
- A backend-neutral AddressSanitizer/UndefinedBehaviorSanitizer build passed
  both tests. LeakSanitizer itself had to be disabled because it cannot run
  under the workspace's tracing environment; this is an environment
  limitation, not a repository failure. The deal.II target was not sanitized
  in this wave.
- The passing suite still establishes current behavior, not completeness or
  independence of every mathematical oracle. No static analysis or coverage
  service was introduced.
- The optional fallback deal.II discovery path was not exercised because the
  CMake package was available. Its CMake logic is assessed structurally below.
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

### S-008 — Compiled consumers see ports, not concrete problem variants

**Evidence:** `CompiledProblemT<Backend>` owns an executable model, metric,
optional constraint, state/adjoint services, and manifest, but exposes only
the backend-neutral contract interfaces
(`include/nmopt/compiler/v1/compiled_problem.hpp:48-130`). The six concrete
deal.II choices occur inside `DealiiCompiler::compile()`; the reduced-gradient
solver does not branch on diffusion, control placement, observation, lifting,
or nullspace policy. Every current compiler-created solver callback captures
the same `shared_ptr` as the executable port, so a reduced DTO detached from
its `CompiledProblemT` currently retains the concrete target indirectly.

**Preserve:** Keep target selection and deal.II types behind compilation. Make
the ownership guarantee explicit during refactoring rather than relying on
the particular captures in the current compiler.

### S-009 — Unsupported realization combinations are rejected explicitly

**Evidence:** The compiler distinguishes structural/analytical validation,
deal.II lowerability, and formulation capability. It rejects matrix-free
execution, missing fixed lifting data, incomplete controlled boundaries,
continuous-control boxes, nonpositive coefficient bounds, incompatible pure
Neumann forcing, and all-at-once formulation requests instead of silently
substituting a nearby implementation
(`include/nmopt/compiler/v1/dealii_compiler.hpp:28-177,738-1407`). The
deal.II contract test exercises representative rejection paths.

**Preserve:** A component-lowering refactor must retain conservative failure.
Composition means combining registered semantics, not accepting an
unregistered cross-product by inference.

### S-010 — Specialized deal.II policies implement real transpose and reconstruction actions

**Evidence:** The fixed-Dirichlet and controlled-Dirichlet targets construct
physical reconstruction/lifting maps and use their transpose pullbacks. The
Neumann target assembles a separate boundary coupling and a mean-zero saddle
system for its pure-Neumann variant. The coefficient target reassembles its
state operator and implements parameter JVP/VJP actions. The H1 target keeps
objective regularisation distinct from the selected search metric. Targeted
tests check residual JVP/VJP pairings, finite differences or reduced Taylor
remainders, metric identities, physical reconstruction, nullspace behavior,
and constraint semantics in
`tests/dealii_diffusion_contract.cc:79-1181`.

**Preserve:** These are useful vertical-slice oracles and policy
implementations. Extract common FE assembly and solve infrastructure without
flattening reconstruction, trace, gauge, or parameter dependence into a
single branch-heavy model.

### S-011 — The first optimizer preserves the covector/metric distinction and uses the correct projected displacement

**Evidence:** `ReducedGradientSolverT` obtains a primal direction only through
`ReducedDTOT::gradient_direction(reduced_derivative, metric)`, measures it
through `metric.apply()`, and never assumes Euclidean coefficient identity
(`include/nmopt/solvers/reduced_gradient.hpp:93-246`). In the constrained
path it computes
`Pi(u - G^-1 j') - u`, measures that projected residual in the selected
metric, and evaluates Armijo decrease with the actual projected trial
displacement rather than the unprojected direction. It also requires a
feasible initial control and rechecks projection feasibility.

**Preserve:** A P6.1 policy split must keep typed covector/primal pairings and
must give projected line searches the accepted displacement. Do not turn a
search direction into a raw coefficient vector or apply unconstrained Armijo
formulas after projection.

### S-012 — Numerical contracts cover the difficult chain-rule and policy boundaries

**Evidence:** The deal.II test contains separate scenarios for the baseline
volume control, fixed reconstruction, controlled-Dirichlet lifting, subdomain
observation, Neumann boundary control and trace observation, H1
regularisation and metric choice, coefficient identification, and the
pure-Neumann gauge. Across them it checks manufactured values, residual
pairings, finite differences, state-recomputed reduced Taylor ratios,
metric identities, box behavior, changed-data invalidation, nullspace
conditions, and rejected capabilities. All three CTests passed three
consecutive repetitions. The backend-neutral suite also passed strict common
warnings and AddressSanitizer/UndefinedBehaviorSanitizer diagnostics.

**Preserve:** Refactor test organization without replacing these mathematical
oracles by constructor-only or snapshot tests. In particular, keep
state-recomputed Taylor tests for every nonlinear or transformation-sensitive
decision path.

### S-013 — The build has a real deal.II-free verification mode and little incidental machinery

**Evidence:** The 59-line root `CMakeLists.txt` uses a C++17 interface target,
keeps deal.II behind `NMOPT_ENABLE_DEAL_II`, and registers the dense and
semantic tests independently of that backend. With the option disabled the
repository configures, builds, and tests without special stubs. There is no
packaging, generated-code, plugin, or application layer to preserve.

**Preserve:** Keep a fast backend-neutral developer gate and do not add
packaging or a generalized build framework solely for architectural
appearance. Build changes should address explicit configuration, diagnostic,
or compile-cost failures identified below.

## Compiler and deal.II review surface

Evidence wave 3 reviewed the complete compiler product, capability, and
binding headers; traced `DealiiCompiler` validation, target selection,
binding, ownership, and manifest construction; and inspected the public
actions, assembly boundaries, solve policies, and retained state of the v0
model and all five v1-only targets. The deal.II test was read by target and by
diagnostic/manifest assertion rather than treated as one opaque passing
executable.

Authoritative records used in this wave were the
[v1 compiler record](semantic-v1-compiler.md), the
[v0 lowerer record](dealii-v0-lowerer.md), the compiler/lowering sections of
the [system blueprint](system-blueprint.md), the manifest requirements in the
[implementation-readiness review](implementation-readiness-review.md), and
the Chapter 6 provenance requirements in the
[numerical-methods guide](chapter-6-numerical-methods-guide.md).

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

### RF-006 — Contract, semantic, and compiler negative tests are too coarse for refactoring

**Classification:** Test debt.

**Evidence:** `tests/semantic_v1_contract.cc` has 12 negative assertions and
`tests/dealii_diffusion_contract.cc` has 13 compiler-negative assertions that
use only `ValidationReport::has_category()`. Neither asserts the available
`Diagnostic::component_id` or `Diagnostic::capability`, so an unrelated error
in the same category can satisfy a test. The compiler alone has 64 diagnostic
insertion sites, making category-only assertions especially weak.
`tests/reduced_dto_contract.cc` checks many successful numerical identities
but has no expected-failure case for layout mismatch, invalid partition,
solver callback, metric, or constraint construction. All current tests pass,
so this is a containment and localization weakness rather than a failing
baseline.

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
component ID, and capability; make negative semantic and compiler cases
table-driven where that improves readability. Assert whether a rejected
compile returns a diagnostic or throws, because RF-010 shows that both paths
currently exist. Add `expect_contract_error` tests for public contract
preconditions and separate the dense test's algebra, DTO, metric, constraint,
and solver cases into focused functions. CTest executable granularity can be
decided in the broader test/build wave.

**Verification:** Demonstrate that each test fails if its expected capability
or component is changed. Run the backend-neutral suite first, then the full
deal.II suite, and preserve all current numerical checks.

**Tradeoff:** Exact assertions require intentional updates when diagnostic
names change. That maintenance cost is the mechanism that prevents silent
diagnostic drift.

## Architecture and code findings: compiler and deal.II realization

### RF-007 — The component registry is a whitelist in front of a whole-problem dispatcher

**Classification:** Architectural debt and a direct gap against the project’s
composition mission.

**Evidence:** `DealiiLowererRegistryV1` is a final, stateless class whose six
methods return booleans from enum switches
(`include/nmopt/compiler/v1/dealii_capabilities.hpp:9-74`). It stores no
lowering functions and cannot be extended without editing the class. After
those per-kind checks, `validate_registered_graph()` accepts exact residual,
data, loss, observation, metric, and constraint sets selected by whole-graph
feature flags (`include/nmopt/compiler/v1/dealii_compiler.hpp:1047-1371`).
`compile()` then derives ten booleans and selects one of six mutually
exclusive complete models at lines 54-378: Neumann boundary, controlled
Dirichlet, H1 control, coefficient identification, fixed/subdomain assembled
v1, or the direct v0 model.

The five v1-only target headers contain 3,171 lines; the direct v0 target adds
542. All six implement the complete executable model, assembly, objective,
state/adjoint solves, and layouts. Some specialization is mathematically
necessary, but common scalar diffusion-reaction assembly, objective
bookkeeping, variable checks, and the same CG helper are repeated. A kind may
therefore be “registered” while no independently composable lowerer for that
kind exists.

**Authority:** The root mission, [composition boundaries](composition-boundaries.md),
and [system blueprint](system-blueprint.md) require residual, observation,
loss, metric, constraint, transformation, and discretization components to be
lowered without a new problem class for every combination.

**Consequence:** Adding tensor coefficients or transport to P5.1, an energy
or weighted observation to P5.2, or a generalized lifting to P5.4 currently
requires modifying exact-graph validation, the central dispatcher, manifest
boolean logic, and at least one complete target. Combining two individually
supported additions tends toward another target class. The public registry
communicates a stronger compositional architecture than the implementation
provides.

A clean warning-enabled build gives this coupling a practical cost: compiling
the one 1,502-line deal.II test translation unit, which reaches every target
through the public compiler aggregate, contributed to a 54.00 s clean build
with 1.49 GiB peak resident memory; its object file was 14 MiB. These local
numbers are not performance requirements, but they resolve `OBS-BLD-03`:
including every target implementation through one public template aggregate
is already material to the edit-build-test loop.

**Scope relevance:** Before affected features. This is a prerequisite for
component-heavy P5.1/P5.2/P5.4 and the P6.2 advection/stabilization target. It
does not block isolated P6.1 solver-policy work, and it must not pre-implement
mixed Stokes, KKT, or complementarity.

**Action tier:** Before the first selected feature that would add a residual,
observation/loss, or transformation combination.

**Recommendation:** Introduce one validated, ID-resolved compiler view and a
small scalar lowering plan. Registered component handlers should contribute
operators, data requirements, objective pieces, metric/constraint services,
and provenance to that plan; registration should hold those handlers or be
honestly named a capability whitelist. Extract shared scalar FE context,
assembly utilities, and solve services. Keep reconstruction maps, boundary
trace maps, mean-zero saddle policy, and parameter-dependent assembly as
explicit composed strategies because their coordinate and transpose behavior
is genuinely distinct. Do not introduce inheritance from PDE families, a
runtime plugin framework, or a universal builder.

This also resolves `OBS-ARC-03`: file length alone is not the finding, but the
955-line semantic validator and 1,579-line compiler each build indexes,
switch on kinds, and issue 70 and 64 diagnostics respectively. A resolved
structural graph should be produced once by the semantic layer and consumed
by focused lowerability/planning code, preserving the category boundary.

**Verification:** Characterize every current graph before extraction. Compile
the same variants through the new plan, compare manifests structurally, and
retain all current value/JVP/VJP, reconstruction, metric, constraint, and
nullspace checks. Add at least one orthogonal recombination test showing that
changing an observation does not require or alter residual assembly, as the
current subdomain test already partially demonstrates.

**Tradeoff:** A full generic weak-form engine would exceed the documented
endpoint. The useful seam is a modest shared scalar assembly plan plus
specialized policies; selected mixed or all-at-once products can add separate
formulation plans later.

### RF-008 — The compilation manifest is incomplete, stringly typed, and wrong for one current target

**Classification:** Provenance correctness defect and architectural debt.

**Evidence:** `CompilationManifest` is a flat aggregate of free-form strings
and component-ID vectors (`include/nmopt/compiler/v1/compiled_problem.hpp:17-46`).
It has no explicit test, adjoint, or observation spaces; mesh identity or
selection; formulation kind; state and adjoint solve tolerances; or binding
provenance. Those fields are required by the implementation-readiness review
and Chapter 6 numerical-methods guide. The manifest builder accepts nine
feature booleans and reconstructs descriptions in nested conditionals
(`include/nmopt/compiler/v1/dealii_compiler.hpp:1419-1539`) instead of being
populated by the lowerers that made the choices.

There is also a current factual error. Coefficient identification constructs
the metric and box with identifier `l2_cellwise_parameter`
(`include/nmopt/compiler/v1/dealii_coefficient_identification.hpp:103-129`),
and the metric manifest names that Riesz map, but
`constraint_realisation` falls through to the generic
`l2_cellwise` description at compiler lines 1486-1490. Existing manifest tests
mostly use substring searches and do not detect the mismatch.

**Authority:** `docs/implementation-readiness-review.md:523-528` requires
semantic IDs, FE and mesh selections, pairings/dual representation, policies,
provenance, and assumptions. `docs/chapter-6-numerical-methods-guide.md:74-78`
also requires state, test, control, adjoint, and observation spaces;
trial/test relation; formulation; and state, adjoint, metric, and KKT
tolerances.

**Consequence:** Two materially different compiled calculations can have
manifests that omit the distinguishing mesh, solve, or data-binding
provenance, while a parameter-box result can actively report the wrong
projection identifier. Free-form substring tests allow wording and semantic
content to drift independently.

**Scope relevance:** Fix the factual error as common maintenance. Manifest
structure is required before publishing Chapter 6 comparisons or adding P6.2
formulation provenance. It is useful but not a reason to invent unselected
OTD/KKT fields now.

**Action tier:** Correct the existing target before new numerical results;
perform the structured refactor before P6.2 or benchmark/reproduction work.

**Recommendation:** Derive a structured manifest from the resolved lowering
plan and compiled services, not from parallel feature booleans. Use typed
subrecords for spaces and pairings, discretization/mesh summary, formulation,
metric/constraint realization, and solver policies; retain human-readable
rendering as an output view. Record the exact projection metric ID. Where a
deal.II `Function` or generated mesh cannot describe itself, require a
caller-supplied provenance label or documented fingerprint policy rather
than pretending the current object is reproducible. Keep the manifest
descriptive: it must not become a second input configuration.

**Verification:** Compare exact structured fields for every current target,
including the coefficient parameter box. Compile two meshes, two solve
policies, and two formulations and require the relevant manifest fields to
differ. Add serialization/rendering round-trip tests only if a persisted
format is selected.

**Tradeoff:** Typed fields require an explicit compatibility/version policy.
That is preferable to parsing prose in tests or downstream experiment tools.

### RF-009 — The canonical v0/v1 comparison shares its assembly implementation

**Classification:** Test-independence debt and overstated documentation.

**Evidence:** The homogeneous fallback branch aliases `DirectModel` to
`dealii_backend::ScalarDiffusionReactionModel<dim>` and constructs it directly
(`include/nmopt/compiler/v1/dealii_compiler.hpp:335-358`). The comparison test
first constructs that same concrete class at
`tests/dealii_diffusion_contract.cc:1194-1202`, then compiles the canonical
graph through the fallback and compares the two instances at lines
1429-1467. The instances are separate and can catch wrong binding, selection,
packaging, layout, metric, or constraint wiring, but the residual, objective,
derivative, and assembly code are identical.

**Authority:** The [v1 compiler record](semantic-v1-compiler.md) describes
this as a v0/v1 comparison guarantee, while the assessment plan requires
independent verification or an explicit record of shared implementation.

**Consequence:** An assembly sign, quadrature, boundary-row, or objective
error in `ScalarDiffusionReactionModel` appears on both sides and cannot be
found by the equality check. Exact equality at `1e-12` looks stronger than the
actual oracle. The other derivative-consistency tests remain valuable but
mostly establish internal calculus identities rather than the intended weak
form independently.

**Scope relevance:** Common characterization concern before RF-007 changes
the compiler. It does not imply maintaining duplicate production assemblers.

**Action tier:** Before implementation of the compiler refactor.

**Recommendation:** Relabel the current assertion as a compiler wiring and
packaging regression. Add an independent small oracle for the homogeneous
weak form: hand-assembled local/global matrices on a tiny mesh, a manufactured
residual/objective with analytically known integrals, or a separate test-only
assembly path. Do not copy the full production model merely to manufacture
independence.

**Verification:** Seed a temporary sign or coefficient fault in each
production contribution and confirm an independent test fails on only the
faulty path. Retain the current compiled/direct equality check for binding and
port regressions.

**Tradeoff:** A tiny oracle covers fewer meshes and policies than the
production test. Its purpose is semantic independence, while the existing
larger tests retain integration breadth.

### RF-010 — Compiler input failures are split unpredictably between diagnostics and exceptions

**Classification:** Compiler API correctness and diagnostic debt.

**Evidence:** `compile()` returns `CompilationResultT` and explicitly reports
missing diffusion/lifting data, bound representation, positivity, and several
mesh/policy failures as lowerability diagnostics
(`include/nmopt/compiler/v1/dealii_compiler.hpp:86-177`). It does not validate
ordinary diffusion, reaction, or regularisation values, bound vector length,
or bound ordering before construction. Those caller-supplied errors throw
`ContractError` from model and constraint constructors instead: for example,
`ScalarDiffusionReactionModel` requires positive diffusion and
regularisation at `include/nmopt/dealii/scalar_diffusion_reaction.hpp:82-91`,
and `CellwiseBoxConstraint` checks mesh-dependent vector sizes and ordering at
`include/nmopt/dealii/cellwise_box_constraint.hpp:33-44`. Facewise bounds have
the same split. Existing tests expect a diagnostic for a nonpositive
coefficient lower bound but do not cover neighboring invalid inputs.

**Authority:** The semantic/compiler records require unsupported bindings and
realizations to be diagnosed, while `ContractError` is the executable
contract’s invariant-failure mechanism. User-provided compilation data should
have a predictable boundary between those two channels.

**Consequence:** Code consuming `CompilationResultT::succeeded()` can still
be terminated by a foreseeable bad binding. As P5 adds coefficient tensors,
boundary data, weights, and observation data, the inconsistent error surface
will grow and category-only tests may hide it.

**Scope relevance:** Common prerequisite for adding binding kinds and useful
for experiment tooling.

**Action tier:** Before the next compiler data-binding extension.

**Recommendation:** Define and document the failure boundary. Validate
user-supplied scalar values and representation shape early; after FE layouts
exist, run a binding-to-layout validation pass that appends exact diagnostics
before constructing constraints or services. Reserve `ContractError` for
violated internal invariants and direct low-level constructor misuse. If the
project intentionally chooses exceptions for invalid concrete data, make
`compile()` consistently translate or document them rather than mixing both
case by case.

**Verification:** Add table-driven cases for nonfinite/negative coefficients,
zero regularisation, mismatched cellwise and facewise vector sizes, reversed
and nonfinite bounds, empty meshes, and incompatible boundary IDs. Assert the
exact chosen error channel, category, component, and capability.

**Tradeoff:** Some binding validation needs assembled layout information, so
compilation may do limited discretization work before returning a diagnostic.
That cost is preferable to exposing partially constructed products or
unpredictable throws.

### RF-011 — Compiled-product lifetime safety is incidental and does not include the mesh

**Classification:** Ownership and lifetime architectural risk.

**Evidence:** `CompiledProblemT` owns the model, metric, constraint, and solver
function objects, but `make_reduced_dto()` passes a reference to the model
(`include/nmopt/compiler/v1/compiled_problem.hpp:60-137`). The current
`DealiiCompiler` callbacks capture their concrete model `shared_ptr`, so that
specific returned DTO keeps its referenced model alive indirectly. The
public backend-generic `CompiledProblemT` constructor does not require such a
capture, so the guarantee is neither represented by the type nor documented.

More importantly, every concrete target constructs one or more deal.II
`DoFHandler`s from the caller’s `Triangulation&`; for example the v0 target
does so at `include/nmopt/dealii/scalar_diffusion_reaction.hpp:64-80`.
`DealiiCompiler::compile()` accepts that borrowed triangulation and the
compiled product retains no mesh owner or lifetime token. The bound
`Function` data are evaluated during assembly and are not retained, which is
good, but the discretization lifetime remains external.

**Consequence:** A generic compiled product can produce a dangling detached
DTO if its callbacks do not happen to own the model. Any compiled deal.II
product can outlive its triangulation at the C++ type level, and mesh mutation
can invalidate the relation among DoF handlers, layouts, and assembled
operators. This weakens the claim that the result is an immutable compiled
discrete problem and is especially risky for Chapter 6 mesh sequences and
experiment drivers.

**Scope relevance:** Resolve before exposing compiled artifacts to longer-lived
drivers, adaptive/refinement workflows, or stored solver sessions. It need
not block local semantic-only work.

**Action tier:** Before Chapter 6 experiment infrastructure or any public API
that detaches/stores compiled services.

**Recommendation:** Choose an explicit ownership model. A deal.II compilation
context can own the triangulation/discretization and be shared by the target,
or a compiled session can own all backend resources and return only
lifetime-bound views. Make an owned reduced formulation service hold the
executable explicitly; retain the existing non-owning `ReducedDTOT`
constructor for short-lived direct use if its lifetime precondition is
documented. Do not scatter `shared_ptr` through mathematical ports without a
clear owner.

**Verification:** Under ASan, return a compiled reduced service from a helper
after the wrapper’s local variables die and evaluate it safely. Destroying or
mutating the mesh while compiled services exist must either be impossible by
ownership, rejected deterministically, or documented and covered as a
precondition.

**Tradeoff:** Owning a mesh changes how callers construct multiple products
on one triangulation. A shared discretization context can preserve reuse while
making invalidation and recompilation explicit.

### RF-012 — Projection compatibility can be spoofed with a metric string

**Classification:** Architectural correctness risk in the metric/constraint
contract.

**Evidence:** Dense, cellwise deal.II, and facewise deal.II box constraints
implement `supports_projection_in()` by checking layout compatibility and
`MetricT::id()` equality
(`include/nmopt/contract/metric_constraint.hpp:159-171`,
`include/nmopt/dealii/cellwise_box_constraint.hpp:75-89`, and
`include/nmopt/dealii/facewise_box_constraint.hpp:75-89`). `MassMetric` has a
public constructor accepting any sparse square matrix and any nonempty ID
(`include/nmopt/dealii/mass_metric.hpp:34-61`). A caller can therefore give a
non-diagonal SPD metric the ID `l2_cellwise`; the box advertises support and
performs coefficient clipping even though that is generally not the metric
projection. For example, with metric matrix
`[[2,1],[1,2]]`, box `[0,1]^2`, and point `(-1,1)`, clipping returns `(0,1)`
while the metric projection is `(0,0.5)`.

The current `DealiiCompiler` constructs honest metric/constraint pairs from
the same target, so no present reference graph triggers the defect. The
public constructors and backend-generic `CompiledProblemT` nevertheless make
the string an unverified mathematical capability token.

**Authority:** The executable contract and implementation-readiness review
require projection to be qualified by the actual selected metric and
explicitly forbid treating coefficient clipping as an H1 or global metric
projection.

**Consequence:** A typo, reused semantic ID, or future custom metric can turn
a false projection into a reported capability. Layout compatibility cannot
establish the separability property that makes clipping exact.

**Scope relevance:** Before adding H-1, trace, obstacle, state, or
complementarity constraints in P5.2/P5.4/P5.5/P6.5. The current compiler path
can remain operational while the contract is strengthened.

**Action tier:** Before the next metric/constraint compatibility extension.

**Recommendation:** Stop using caller-controlled display IDs as the proof of
projection semantics. Have the owning lowerer construct a coupled projection
service with an opaque/strong realization identity, or otherwise provide a
non-spoofable compatibility witness tied to the actual metric operator.
Continue exposing a human-readable metric ID for provenance, but separate it
from capability identity. Keep projection owned by the constraint, as the
interface specification requires.

**Verification:** Construct a compatible-layout, same-string non-diagonal
metric and require rejection. Retain positive tests for diagonal DGQ(0) cell
mass, face mass, and the dense diagonal reference; add an explicit negative
H1 case.

**Tradeoff:** Runtime type checks alone would couple constraints to concrete
metric classes and still over-accept arbitrary matrices. An opaque witness or
coupled factory adds a small amount of plumbing but expresses the actual
compiled relationship.

### RF-013 — State and adjoint solve policy is hard-coded and repeated across targets

**Classification:** Reproducibility, architecture, and code-quality debt.

**Evidence:** The v0 model and all five v1-only targets contain near-identical
CG helpers using identity preconditioning, maximum iterations
`max(100, 10 * matrix.m())`, and tolerance
`max(1e-14, 1e-12 * ||rhs||)`. The six copies occur in
`scalar_diffusion_reaction.hpp:500-512`,
`dealii_fixed_dirichlet.hpp:638-652`,
`dealii_dirichlet_control.hpp:636-650`,
`dealii_neumann_boundary.hpp:578-590`,
`dealii_h1_control.hpp:467-479`, and
`dealii_coefficient_identification.hpp:561-575`.
`DealiiDiscretisationPolicy` exposes parameters only for the control metric
inverse, and the manifest states only prose about state/adjoint CG or direct
solution, not the iterative tolerances.

**Authority:** The implementation-readiness review requires solve policy and
the Chapter 6 guide requires separate state, adjoint, metric-inverse, and KKT
tolerances in the manifest and numerical reports.

**Consequence:** Numerical experiments cannot select or faithfully record
state/adjoint accuracy, a solver failure is exposed through deal.II rather
than a uniform compiled-service diagnostic, and changing the default requires
six edits. P6.1 line-search and gradient comparisons can be contaminated by
unreported inner accuracy.

**Scope relevance:** Required before Chapter 6 numerical comparisons. It can
be deferred while doing semantic-only Chapter 5 work, but shared solve
extraction is naturally part of RF-007.

**Action tier:** Before P6.1 results are treated as comparable and before P6.2
or later formulation studies.

**Recommendation:** Add a small typed SPD solve policy and one shared serial
deal.II solve service that returns convergence information. Permit state and
adjoint policies to be selected/recorded separately even when they share a
default. Keep the pure-Neumann direct saddle solve as a distinct declared
policy. Do not build a universal linear-solver framework before a selected
target needs additional algorithms.

**Verification:** Test that changing maximum iterations and tolerances changes
the invoked policy and manifest; cover zero RHS, convergence, deliberate
nonconvergence, and state/adjoint residual criteria. Re-run derivative tests
at sufficiently tight declared tolerances.

**Tradeoff:** Exposing inner-solve failure forces callers and solvers to decide
how to stop or report it. That explicit behavior is necessary for the
project’s software-engineering and numerical-methods goals.

## Formulation, solver, test, and build review surface

Evidence wave 4 fully reviewed the only current optimizer,
`include/nmopt/solvers/reduced_gradient.hpp`, and its reduced-formulation
service in `include/nmopt/contract/reduced_dto.hpp`; inventoried every logical
test scenario and its numerical checks; and reviewed the complete root CMake
configuration and repository-local developer tooling. Relevant authority came
from the reduced and verification sections of the
[executable contract](executable-contract-v0.md),
[implementation-readiness review](implementation-readiness-review.md),
[implementation roadmap](implementation-roadmap.md), and
[Chapter 6 numerical-methods guide](chapter-6-numerical-methods-guide.md).

There is no separate formulation package today. That is appropriate for the
implemented slice: `StateControlPartitionT` plus `ReducedDTOT` is explicitly
one-state/one-decision, DTO, and first-order. It should not be generalized to
represent P6.2 OTD or P6.3 KKT products. Those features need separately named
products if selected. The current refactor question is therefore limited to
how reduced algorithms request state/objective/adjoint work and how they
record it.

### Current verification map

| Logical scenario | Strongest current checks | Important gap or characterization need |
| --- | --- | --- |
| Dense linear-quadratic DTO | Residual JVP/VJP, residual and objective finite differences, state residual, reduced derivative, metric identity, box clipping, unconstrained/projected convergence | No direct adjoint-equation residual; solver failure/limit paths and exact accepted Armijo data are absent |
| Alternate dense backend | Backend-parametric pairing and covector subtraction | Only a very small algebra instantiation; no DTO/metric/solver instantiation |
| Semantic v1 graph | Positive reference graphs and twelve rejected mutations | Category-only rejection assertions; RF-002 through RF-005 are uncharacterized |
| Baseline deal.II volume control | State residual, residual JVP/VJP and finite difference, metric action, reduced derivative, projected/unconstrained solves, v0/v1 comparison, representative compiler diagnostics and manifest | v0/v1 assembly is shared (RF-009); only one mesh and the v0 solver success path |
| Fixed Dirichlet | Manufactured objective, reconstructed state residual, JVP/VJP, residual/objective derivatives, state-recomputed Taylor ratio, changed-data invalidation, manifest | Reconstruction is tested through the whole target rather than an independently owned transformation port |
| Controlled Dirichlet | Physical-state reconstruction, state residual, composed JVP/VJP, reduced Taylor ratio, trace metric, unsupported partial boundary, manifest | White-box `dynamic_cast` to a compiler `detail` target; no standalone lifting value/JVP/VJP contract |
| Subdomain observation | Changing the material mask leaves the state/residual unchanged and changes objective/adjoint RHS; exact manifest string | Differential behavior only; no independent restricted-integral value/derivative or reduced Taylor check |
| Neumann boundary control | State residual, boundary coupling JVP/VJP, trace-objective derivative, reduced Taylor ratio, face metric and clipping, manifest | No independent face-integral oracle or degenerate empty/incorrect marker case at the compiled boundary |
| H1 regularisation/metric | Unsupported combinations, stiffness contribution, L2 versus H1 direction distinction, H1 metric identity, state residual, pairing, reduced Taylor ratio, manifest | Absolute tolerances on one mesh; no solver behavior under the H1 metric |
| Coefficient identification | Positive-bound diagnostics, state reassembly, nonlinear residual finite difference and pairing, objective derivative, reduced Taylor ratio, metric/constraint | No solver iteration or reassembly/work-count report; manifest misses the constraint bug in RF-008 |
| Pure Neumann | Incompatible data/control rejection, repeatability, state/adjoint mean-zero gauge, state residual, manifest | White-box cast for gauge inspection; no direct saddle residual or adjoint equation residual |

The existing tests are correctness contracts, not Chapter 6 reproduction
drivers. That separation is good. Mesh/parameter sweeps, convergence trends,
hardware, timings, and published iteration comparisons should remain outside
CTest and should be added only for selected examples with a manifest-bearing
run record. Their present absence is not refactor debt.

## Architecture and engineering findings: solver, tests, and build

### RF-014 — The reduced evaluation protocol performs an adjoint solve for every objective-only trial

**Classification:** Solver/formulation interface debt with immediate P6.1
performance and extensibility impact.

**Evidence:** `ReducedDTOT` exposes one operation, `evaluate(control)`, which
always solves the state, evaluates the full objective derivative, solves the
adjoint, applies the residual VJP, and only then returns the objective value
(`include/nmopt/contract/reduced_dto.hpp:127-170`). The Armijo loop calls that
full operation before deciding whether a trial objective is acceptable
(`include/nmopt/solvers/reduced_gradient.hpp:151-181`). The dense and deal.II
tests explicitly assert
`state_solve_count == adjoint_solve_count == line_search_trial_count + 1`, so
this is established behavior, not a speculative cost.

**Authority:** Armijo acceptance needs a trial objective and state solve, not
a trial derivative. P6.1 then adds policies with different information needs:
exact quadratic search needs a Hessian action, while Wolfe conditions need a
trial derivative. The executable contract should make those costs explicit.

**Consequence:** Every rejected backtracking step pays for an unused adjoint;
the current `initial_step_length = 20` deal.II test deliberately exercises
that waste. More importantly, a P6.1 implementation would either keep
over-solving or add algorithm-name branches around the monolithic evaluation.
The reported solve counts describe calls to `evaluate`, not work requested by
the chosen line search.

**Scope relevance:** Required before implementing selected P6.1 line searches
or treating solver work counts as numerical results. It does not block
semantic-only Chapter 5 work or require RF-007 first.

**Action tier:** Target-dependent pre-feature refactor for P6.1.

**Recommendation:** Split the reduced service by requested information, for
example into a value record containing state/full point/objective and a
derivative operation that augments an accepted value record with adjoint and
reduced covector. Let the iteration driver request value-only Armijo trials
and compute a derivative after acceptance. Add small direction and
line-search policy concepts only for methods actually selected under P6.1;
do not create a universal algorithm registry. Preserve the current
`ReducedDTOT::evaluate()` as a convenience composition if that reduces
migration risk.

**Verification:** Use counting state/adjoint callbacks. A rejected Armijo
trial must increment only the state count; every point at which a stopping
norm or new direction is computed must have one matching adjoint. Add exact
quadratic and Wolfe tests only when those policies are implemented, each
asserting the work it requested.

**Tradeoff:** A split record introduces an accepted-state cache and therefore
more lifetime/layout checks. It pays for itself by making algorithm costs and
derivative validity explicit, and it remains much smaller than a generalized
formulation engine.

### RF-015 — Solver reporting cannot yet audit a Chapter 6 run or carry its compilation provenance

**Classification:** Reproducibility and failure-contract debt; partial P6.1
feature gap.

**Evidence:** `ReducedGradientResultT` records the final control, accepted
objective history, stopping-norm history, four aggregate counts, and one of
three reasons (`gradient_tolerance`, `maximum_iterations`, or
`line_search_failure`) at
`include/nmopt/solvers/reduced_gradient.hpp:11-47`. It does not record
accepted step lengths, relative norm, objective change, actual directional
decrease, metric solves, inner linear iterations/residuals, or a snapshot of
solver parameters. Non-descent, infeasible projection output, bad metric
actions, and state/adjoint callback failures escape as exceptions rather than
typed terminal outcomes. The compilation manifest is accessible from
`CompiledProblemT`, but the solver consumes a detached DTO and its result has
no run-level association with that manifest.

**Authority:** P6.1 requires uniform step, stop, line-search, and separated
work counts. The implementation-readiness review requires compilation
provenance to accompany solver diagnostics, and the numerical-examples guide
requires solver policies, tolerances, mesh policy, and hardware/run metadata
for a reproduction.

**Consequence:** Objective monotonicity can be checked, but the current result
cannot independently reconstruct whether each accepted step met its declared
Armijo inequality. A saved result can be separated from the discretization
and policy that produced it. As new directions add restart or curvature
failure, using `ContractError` for all non-success behavior would confuse
invalid program wiring with normal algorithm termination.

**Scope relevance:** Required before publishing or comparing Chapter 6
numerical results and naturally implemented with P6.1. Current v0 regression
tests can continue using the smaller result until then.

**Action tier:** Target-dependent pre-feature refactor for P6.1/results; not a
common prerequisite for Chapter 5.

**Recommendation:** Define a typed optimization report with an initial/final
and relative stationarity measure, per-accepted-iteration objective and step,
line-search trials, actual descent pairing, separated state/adjoint/metric
work, and an explicit terminal status. Keep invalid layouts and violated
metric/constraint contracts as exceptions; represent expected algorithmic
outcomes such as failed line search, curvature reset/failure, nonfinite trial,
or inner-solve failure in the report when recovery or comparison is
meaningful. Build a small experiment/run envelope outside the generic solver
that pairs the compilation manifest, solver policy snapshot, report, and
environment metadata. Do not make the backend-neutral solver depend on a
deal.II manifest type.

**Verification:** Assert every accepted inequality from reported fields;
check relative and absolute stopping, maximum-iteration and forced
line-search failure paths, nonfinite trial handling, and exact work counts.
Serialize or print one run envelope and show that formulation, mesh, state,
adjoint, metric, and outer-solver policies can all be identified without
inspecting live objects.

**Tradeoff:** Per-iteration records consume more memory than two scalar
histories. Chapter 6 experiments need that evidence; a configurable compact
mode can be considered only if a selected large run demonstrates a problem.

### RF-016 — Eleven logical test scenarios are hidden behind three fail-fast process entries

**Classification:** Test architecture and characterization debt.

**Evidence:** `nmopt_contract_tests` contains two logical functions,
`nmopt_semantic_v1_contract_test` one, and the 1,502-line deal.II executable
eight templated scenario functions. Each `main()` runs them sequentially
inside one catch block, so the first failure prevents all later cases from
running. CTest exposes only three entries and assigns no label or timeout.
The controlled-Dirichlet and pure-Neumann cases `dynamic_cast` the public
executable port to compiler `detail` model types to inspect reconstruction or
gauge behavior. All repetitions passed, so the finding is diagnosis and
coupling, not observed flakiness.

Coverage is strong but incomplete at important refactor seams: there is no
direct adjoint-equation residual check, no solver test for maximum-iteration
or line-search-failure results, no exact accepted Armijo record, and no tests
for callback failure or nonfinite values. Absolute `require_close` tolerances
are repeated locally and are suitable for the fixed meshes but will not scale
to refinement/parameter sweeps. RF-001 through RF-006, RF-008, RF-010, and
RF-012 identify additional missing regressions.

**Authority:** The assessment plan requires failure localization, independent
oracles, degenerate cases, and a practical refactor loop. The Chapter 6 guide
requires direction, acceptance, stopping, and work-count tests independently.

**Consequence:** A failure in the first deal.II scenario can hide seven later
results, focused execution is unavailable through CTest, and a test-only need
to inspect a specialized policy leaks through the same aggregate used for
black-box compiler contracts. Splitting the 1,502-line source into many
executables immediately would instead multiply the already expensive public
template compilation.

**Scope relevance:** Characterization work is a common Stage B prerequisite;
full reorganization can accompany RF-007 and RF-014.

**Action tier:** Add missing defect characterizations before production
refactors; improve process granularity early in Stage B.

**Recommendation:** First add exact regression cases to the existing binaries.
Then let the deal.II executable dispatch one named scenario from a command
argument and register each scenario as a separately named/labeled CTest,
without creating eight translation units. Extract deterministic value,
pairing, Taylor, diagnostic, and solver-report helpers. Keep black-box compiler
port tests separate in intent from white-box realization-policy tests; after
RF-007 extracts lifting/gauge components, test those components directly
instead of casting a compiled executable to `detail`. A third-party test
framework is not necessary to obtain these benefits.

**Verification:** CTest must list and run each logical scenario separately,
one failing case must not suppress unrelated results, and labels must permit
backend-neutral versus deal.II selection. Mutate an expected diagnostic,
transpose sign, and accepted step to confirm the corresponding focused test
fails. Keep the aggregate full run within a documented developer budget.

**Tradeoff:** More CTest processes may repeat deal.II startup. Reusing one
binary and measuring the resulting suite is the low-cost first step; only
split translation units where dependency boundaries or parallelism outweigh
compile cost.

### RF-017 — Enabling deal.II can silently produce a build with no deal.II target or test

**Classification:** Build correctness and dependency-discovery debt.

**Evidence:** `NMOPT_ENABLE_DEAL_II` defaults to `ON`. CMake first calls
`find_package(deal.II CONFIG QUIET)`, then falls back to independent
`find_path` and `find_library` calls. If neither path succeeds, configuration
reaches the end without a status, warning, or error and simply does not create
`nmopt_dealii_contract` or its CTest (`CMakeLists.txt:22-59`). The manual
fallback links one library and include directory but does not call
`deal_ii_setup_target` or propagate the package's dependency libraries,
compile definitions, feature flags, and configuration selection. It was not
executed here because the package configuration was available.

**Authority:** The build convention permits a deal.II-free environment but
requires skipped backend checks to be reported. An option described as
“Build the serial deal.II reference lowerer” should not silently mean “try to
build it.”

**Consequence:** A contributor can request/default to deal.II support, see a
successful configure and green two-test run, and reasonably mistake it for a
full verification. The fallback may configure successfully but fail later or
use incomplete deal.II settings, producing a harder-to-diagnose path than a
clear dependency error.

**Scope relevance:** Common developer-workflow fix; independent of which
Chapter 5/6 features are selected.

**Action tier:** Early Stage B build hardening, before relying on automated
gates.

**Recommendation:** Prefer the official deal.II config package as the one
supported path. If `NMOPT_ENABLE_DEAL_II=ON` and it is unavailable, fail with
an actionable message; users who intentionally need the backend-neutral mode
already have `NMOPT_ENABLE_DEAL_II=OFF`. Remove the manual fallback unless a
real supported environment can verify all required imported-target behavior.
Always print a concise configuration summary including deal.II test
availability.

**Verification:** Configure with the package available, with the option
explicitly off, and with the option on while package discovery is disabled.
The first must register three test groups, the second two, and the third must
fail at configure time with the documented remedy.

**Tradeoff:** Removing an unverified fallback can reduce apparent portability.
For deal.II, a reliable explicit package requirement is preferable to a path
that cannot reproduce the dependency's CMake usage requirements.

### RF-018 — The developer build profile leaves generator, configuration, warnings, and test intent implicit

**Classification:** Developer-experience and reproducibility debt.

**Evidence:** The prescribed command fixes Ninja and the directory `build/`
but omits `CMAKE_BUILD_TYPE`. In this checkout, reusing an older Makefiles
cache failed with CMake's generator-mismatch error; a fresh configure let
deal.II warn and force the empty build type to Debug. The repository has no
CMake presets, project warning option, CI configuration, or CTest labels and
timeouts. Backend-neutral code is clean under `-Wall -Wextra -Wpedantic
-Wconversion -Wshadow`, while the deal.II build exposes the checked-size
boundary in RF-019. Address/undefined-behavior checks are cheap and pass for
the backend-neutral targets when environment-incompatible leak detection is
disabled.

`OBS-BLD-01` is therefore not a source-state or Git risk: generated build
trees are ignored and a commit/tag preserves all authored state. It is a
recoverable CMake cache conflict. The workflow weakness is that the canonical
command provides neither an alternate build-directory convention nor a
documented recovery path. `OBS-BLD-02` is reproducibility debt because an
implicit Debug build is unsuitable as the basis for Chapter 6 performance
observations.

**Authority:** The build convention asks for reproducible out-of-source Ninja
checks and visible failures. The numerical-examples guide separates
correctness from performance and requires the actual configuration to be
recorded.

**Consequence:** Two contributors can run the same documented command but use
different cached settings, or fail before configuration because of an old
generator. Backend-neutral warnings depend on manually supplied flags.
Performance-sensitive experiments risk being run in an implicitly forced
Debug configuration, while the absence of labels makes “fast core” versus
“deal.II integration” an informal distinction.

**Scope relevance:** Common workflow improvement, but not a reason to delay
all semantic fixes. It becomes important before sustained implementation and
mandatory before numerical benchmarking.

**Action tier:** Early Stage B workflow batch; explicit Release/benchmark
profile before any Chapter 6 reported timings.

**Recommendation:** Provide one canonical explicit Debug developer profile
and one backend-neutral profile, either as small CMake presets or fully
specified commands using distinct generated directories. Document generator
cache recovery without deleting a user's build tree automatically. Apply a
project-owned warning interface to project targets without making warnings
from external headers fatal. Add CTest labels and proportionate timeouts.
Use the cheap backend-neutral warning and sanitizer configurations as local or
automated gates. A minimal backend-neutral CI job would address the concrete
absence of any automatic gate; a deal.II CI job should be added only where a
stable package image is actually maintainable. Packaging, installation,
coverage hosting, and a formatter gate are not currently justified.

**Verification:** A new checkout must run the documented profile without an
implicit build-type warning; the backend-neutral profile must not search for
deal.II; test labels must select the intended subsets; and a stale different-
generator tree must have a safe documented recovery. Record configuration in
any experiment run envelope from RF-015.

**Tradeoff:** Presets and extra configurations add a small maintenance matrix.
Keeping only profiles that are actually exercised avoids turning build
polish into a parallel product.

### RF-019 — Contract dimensions and serial deal.II vector sizes meet through unchecked narrowing conversions

**Classification:** Backend adapter safety and warning-hygiene debt.

**Evidence:** Layout dimensions and backend-neutral APIs use `std::size_t`.
On the audited deal.II 9.5.1 build, `dealii::Vector<double>::size_type` is
`unsigned int`. `SerialBackend::zeros(size_t)` passes the value directly to
the vector constructor (`include/nmopt/dealii/serial_backend.hpp:17-21`), and
the two box constraints allocate and index deal.II vectors with `size_t`.
`-Wconversion` reports the same boundary throughout target assembly helpers,
including vector construction and sparse-matrix row/column insertion.

**Consequence:** Present two-dimensional tests are far below the native limit,
so no current result is wrong. The API nevertheless has no checked point at
which an oversized backend-neutral dimension is rejected; a large value can
truncate before later layout/vector agreement checks. The volume of warnings
also prevents a strict conversion profile from distinguishing newly
introduced conversions.

**Scope relevance:** Not a blocker for the bounded small serial Chapter 5/6
targets, but a concrete safety cleanup and prerequisite for making conversion
warnings actionable. It must not be misrepresented as distributed or
large-scale backend support.

**Action tier:** Opportunistic Stage B adapter cleanup, after the correctness
characterizations and before enabling the strict warning profile routinely.

**Recommendation:** Centralize checked conversion from contract dimension to
the backend's native size/index type in `SerialBackend` or a small deal.II
adapter utility. Use native deal.II index types inside vector/matrix loops and
convert back to `size_t` only through widening paths. Apply the helper first at
all allocation and sparse-index boundaries; do not scatter unchecked casts to
silence warnings.

**Verification:** Add a backend adapter test that accepts the maximum
representable native size symbolically without allocating it and rejects a
larger `size_t` when the platform permits one. Rebuild under the strict warning
profile and confirm project-owned conversions are either checked or absent.

**Tradeoff:** The bounded serial examples will never allocate billions of
entries, so the runtime defect risk is low. A centralized boundary still has
learning value and removes noise without redesigning layout dimensions or
claiming a distributed backend.

## Resolution of prior cross-layer observations

| Observation | Resolution |
| --- | --- |
| `OBS-BLD-01`: Ninja conflicts with the existing Makefiles cache | Resolved by RF-018 as a recoverable ignored-cache conflict plus a missing documented/profiled recovery path, not a source or Git-preservation risk |
| `OBS-BLD-02`: deal.II forces an unspecified build type to Debug | Resolved by RF-018 as reproducibility debt; explicit Debug and benchmark profiles are needed |
| `OBS-TST-01`: all deal.II behavior is one CTest | Resolved by RF-016; eight logical scenarios need focused registration without multiplying expensive translation units |
| `OBS-TST-02`: tests cast to compiler `detail` targets | Resolved by RF-016; separate black-box compiler intent from white-box policy intent and test extracted policies directly after RF-007 |
| `OBS-BLD-03`: the compiler aggregate includes all large targets | Resolved by RF-007 with measured build evidence; decompose lowering first and measure before choosing explicit instantiation or additional binary-library machinery |
| `OBS-INT-01`: manifests should accompany solver diagnostics | Resolved by RF-015; use a run envelope outside the backend-neutral solver rather than coupling the solver to compiler manifest types |

## Pending review waves

1. Documentation consistency, conventions, and agent guidance.
2. Final synthesis, action tiers, dependency-ordered refactor batches,
   deferred work, and user decisions.
