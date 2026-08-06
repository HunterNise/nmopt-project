# Pre-Chapter 5/6 refactor assessment

## Assessment status

**In progress — evidence wave 1 of 5 complete.**

This document records the assessment defined by the
[refactor assessment plan](pre-chapter-5-6-refactor-plan.md). The current
contents establish the reproducible baseline and the Chapter 5/6 scope
matrix. Architecture, code, test-quality, build-design, and documentation
findings remain pending until their corresponding review waves are complete.

No code, tests, CMake, conventions, or agent instructions may be changed while
this assessment is in progress.

## Executive verdict

**Pending.** The build and test baseline is healthy, but it is too early to
decide whether a cross-cutting refactor is warranted. The documented endpoint
contains several independent, target-dependent capability families; the
assessment must first determine whether current implementation structure
supports them compositionally or repeats whole-problem machinery.

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

## Findings register

No refactor finding is confirmed in evidence wave 1. The following baseline
observations are queued for later classification:

| Observation | Evidence needed before classification |
| --- | --- |
| `OBS-BLD-01`: the prescribed `build/` Ninja command conflicts with an existing Makefiles cache | Determine whether this is only local stale state or a documentation/workflow weakness reproducible for contributors |
| `OBS-BLD-02`: an unspecified build type is silently forced to `Debug` by deal.II after a warning | Review CMake policy, expected developer configurations, and whether performance-sensitive examples need an explicit choice |
| `OBS-TST-01`: all deal.II behavior is registered as one CTest executable | Review internal test organization, failure handling, runtime, and whether granularity impedes diagnosis or focused refactoring |

## Pending review waves

1. Typed contracts, dense reference behavior, semantic graph, and validation.
2. Compiler ownership, manifests, capability diagnostics, v0/v1 independence,
   and deal.II realizations.
3. Formulations, solvers, tests, numerical verification, CMake, and tooling.
4. Documentation consistency, conventions, and agent guidance.
5. Final synthesis, action tiers, dependency-ordered refactor batches,
   deferred work, and user decisions.
