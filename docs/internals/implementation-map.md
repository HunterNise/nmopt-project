# Implementation map

## Purpose and authority

This document maps the current `nmopt` implementation architecture onto the source tree.
It is for maintainers who already understand the project model and need to know where
a responsibility is implemented, which direction dependencies run, and which tests
provide executable evidence.

Current source and focused tests are authoritative for implementation behavior. This
document explains their structure; `docs/reference/` remains the authority for the
public programming interface.

It is not a public API reference and it is not a replacement for the design record.
For user-facing workflows, start with the [manual overview](../manual/overview/README.md)
and the [public reference](../reference/). For the design principles that guided the
project, see the [architecture record](../design/architecture.md) and
[PDE/formulation/solver boundary](../design/pde-solver-boundary.md).

The source tree implements two peer producer paths:

```text
user or application code
          │
          ▼
     choose producer path
          │
          ├──────────────────────────────┐
          ▼                              ▼
semantic/compiler producer        application-owned producer
      ProblemSpec                  existing PDE application
          │                              │
          ▼                              ▼
 semantic resolution              thin numerical adapters
          │                              │
          ▼                              │
    deal.II compiler                     │
          │                              │
          └──────────────┬───────────────┘
                         ▼
          shared numerical contracts
       model / solves / metric / constraint
                         │
                         ▼
               formulation products
      reduced DTO / supplied OTD / KKT / PDAS
                         │
                         ▼
                       solvers
                         │
                         ▼
          application / experiment evidence
```

The paths converge at numerical and formulation contracts. They do **not** converge
at `ProblemSpec`, `DealiiCompiler`, or `CompiledProblemT`.

## Reusable numerical core

### `include/nmopt/contract/`

This is the backend-parametric numerical vocabulary used by formulations and
algorithms.

The lowest-level pieces are:

- `linalg.hpp` – the reference dense backend and contract failures;
- `layout.hpp` – block layouts, space identifiers, `PrimalBlockT`,
  `CovectorBlockT`, pairing, and block extraction;
- `linear_solve.hpp` – backend-neutral solve results and convergence evidence.

The executable PDE boundary is in:

- `executable_model.hpp`;
- `callback_executable_model.hpp`;
- `metric_constraint.hpp`;
- `reduced_hessian.hpp`.

`ExecutableModelT<Backend>` exposes residual, residual JVP, residual VJP, objective,
and objective derivative actions. State and adjoint inversion are deliberately
separate services. A derivative remains a covector until a `MetricT<Backend>` maps it
to a primal search direction.

Formulation-specific contracts then build on that surface:

- `reduced_dto.hpp` – reduced state/adjoint orchestration;
- `supplied_otd.hpp` and `supplied_otd_kkt.hpp` – supplied all-at-once/OTD systems
  and their validated KKT bridge;
- `quadratic_kkt.hpp` and `quadratic_kkt_solver.hpp` – equality-constrained
  quadratic KKT products and solver policy;
- `complementarity.hpp` and `pdas.hpp` – box complementarity and PDAS.

These products deliberately have different required actions. There is no universal
`Formulation` base class.

### `include/nmopt/solvers/`

This directory owns reusable optimization algorithms that consume the numerical
contracts rather than semantic/compiler objects.

The reduced path is split into:

```text
reduced search policy / directions
              │
              ▼
        line-search solver

reduced trust-region policy
              │
              ▼
       trust-region solver
```

The solvers do not need to know whether their numerical services came from the
semantic compiler or an existing application.

### `include/nmopt/dealii/`

This directory contains reusable deal.II-backed numerical services, not a second
problem-description system.

Important examples are:

- `MassMetric`;
- `Hminus1Metric`;
- `TraceHhalfMetric`;
- `CellwiseBoxConstraint`;
- `FacewiseBoxConstraint`;
- `IndependentStateCoordinates`;
- serial SPD, KKT, and PDAS solve support.

These components are usable numerical building blocks. They do not own semantic
problem selection, compiler orchestration, or application configuration.

`serial_backend.hpp` supplies the deal.II vector backend required by the generic
contract templates.

## Semantic producer

### `include/nmopt/semantic/v1/`

The semantic layer describes a problem independently of deal.II objects.

`types.hpp` owns the typed graph vocabulary. `problem_spec.hpp` is the public
aggregate for that vocabulary. `problem_library.hpp` contains reusable semantic
problem factories and deltas used by tests and applications.

Validation is physically split as:

```text
validation.hpp
      │
      ▼
detail/validation_detail.hpp
```

`validation.hpp` is the small public facade. The inline detail file contains the
ordered validation phases. It remains one implementation unit because the phases and
helpers are tightly coupled.

`resolved_problem.hpp` adds `SemanticResolver` and `ResolvedProblemView`. A resolved
view indexes every stable semantic ID once and borrows the source `ProblemSpec`; it is
intended for one validation/compilation operation rather than long-term storage.

Semantic validation checks graph structure and declared policies. Backend
lowerability and formulation/product capability are compiler responsibilities.

## deal.II compiler producer

### `include/nmopt/compiler/v1/`

The compiler owns the transition from a valid semantic graph plus concrete runtime
bindings to deal.II numerical products.

The main source responsibilities are:

```text
dealii_types.hpp
    runtime bindings, discretisation policy, closed compiler request,
    target-family decisions, compilation session

dealii_scalar_plan.hpp
    bounded component-lowering registry and scalar plans

dealii_* target headers
    typed numerical realizations for specialized target families

dealii_compiler.hpp
    orchestration, validation, target selection, construction,
    provenance finalization, product packaging

compilation_manifest.hpp
    typed compilation/provenance schema

compiled_products.hpp
    owning compiled products

compiled_application_view.hpp
    optional compiler-produced application sidecar

compiled_problem.hpp
    compatibility aggregate
```

`compiled_problem.hpp` is intentionally only an aggregate after the manifest/product
split. New implementation code can include the narrower owning headers.

The compiler is described in detail in [Compiler implementation](compiler.md).

### Why `dealii_compiler.hpp` remains large

`dealii_compiler.hpp` contains recognizable physical regions, but those regions still
share the closed request, semantic lookups, lowerability state, realization data, and
product construction. Splitting them only to reduce file length would add
cross-file navigation and repeated template parsing without establishing a stronger
architectural boundary.

A future extraction is justified when a compiler phase has a stable typed input/output
boundary, narrower dependencies, independently focused tests, or measured build-time
benefit. The current source-navigation comments and typed seams are therefore the
intended implementation map.

## Application-owned producer

The peer producer path begins with an application that already owns its mesh, finite
element spaces, operators, solves, and output.

The reusable boundary is still `include/nmopt/contract/`. An application may create
an `ExecutableModelT` implementation or callback model, provide state/adjoint solve
services and a metric, and construct the formulation required by its algorithm.

The worked external deal.II path lives under:

```text
apps/external-dealii/step-4/
├── source/
├── integration/
├── minimal/
├── evaluation/
├── verification/
└── diagnostics/
```

The `minimal/` consumers are the canonical application-facing examples. Evaluation,
comparison, and verification machinery is intentionally outside that path.

The minimal consumers have a transitive source-level dependency on the optional
Step-4 instrumentation type through integration headers, but they do not instantiate
instrumentation at runtime, execute instrumentation-recording branches, or depend on
the native reference optimizer/oracle/comparison machinery.

The operational integration route is documented in
[External deal.II solver integration](../reference/external-dealii-solver-integration.md).

## Application and experiment layer

### `include/nmopt/application/`

This layer owns reusable application-facing composition above the numerical core.

Important generic pieces are:

- `recipe.hpp` – typed builders for semantic `ProblemSpec`s;
- `scenario.hpp` – typed bundles of problem, compile, solver, and experiment
  options;
- `catalog.hpp` and `metadata.hpp` – application discovery metadata;
- `harness.hpp` – benchmark identity and artifact finalization;
- `runner.hpp` – generic `HeadlessBenchmarkRunnerT`.

Chapter-specific headers in the same tree are concrete project applications, not new
core numerical abstractions.

`HeadlessBenchmarkRunnerT` is generic in both the problem builder and execution
adapter. Its builder return type is inferred; the reusable runner has no intrinsic
requirement that the builder return `ProblemSpec`. The repository Chapter-6 runner
chooses a semantic/compiler builder because that is what B1/B2 use.

### `include/nmopt/experiment/`

`reduced_envelope.hpp` associates a compilation manifest, solver-policy snapshot,
solver report, and caller-supplied environment record.

The envelope owns values. It intentionally does not retain `CompiledProblemT` or
`ReducedDTOT`, so detached evidence does not silently extend executable-service
lifetime.

The module remains a separate top-level area because it owns detached experiment
evidence and no stronger ownership boundary currently justifies folding that
responsibility into the application or solver layers.

## Repository applications

### `apps/nmopt-runner/`

This is the configured Chapter-6 command-line application.

Its current private structure separates:

```text
parameter parsing / binding
          │
          ▼
      RunSetPlan
          │
          ▼
CLI and run lifecycle
          │
          ▼
Chapter-6 execution registration
          │
          ▼
B1 / B2 typed scenarios and adapters
          │
          ▼
artifacts / native output / run manifest
```

`chapter6_execution.hpp` is intentionally an application-specific execution seam. It
was extracted so tests can call production execution logic without including
`main.cc`; it was not intended to become a universal execution framework.

See [Runner implementation](runner.md).

### `apps/external-dealii/`

This area demonstrates the application-owned producer route. It is kept under
`apps/` because concrete consumers and showcases belong there. The `apps/`
tree remains shallow because the current set is too small for another taxonomy layer
to improve navigation.

## Configuration, evidence, and tools

### `parameters/`

Tracked parameter files are executable inputs.

The important distinction is:

```text
docs/studies/    source interpretation and study records
parameters/      versioned executable configuration
runs/            generated evidence
```

`authoritative.prm` means the repository's authoritative reproduction input, not
that every contained value is necessarily a directly transcribed publication fact.

### `tools/`

Two major responsibilities are visible:

- `nmopt_artifacts/` – dependency-free parsing of persisted artifact records;
- `nmopt_postprocess/` – field/post-processing configuration and rendering.

Chapter-6 report/workflow scripts remain project-specific orchestration.

Persisted artifact history is parsed through the shared dependency-free record
implementation so report and plotting code do not maintain different
interpretations of malformed or non-finite evidence.

## Tests as an executable architecture map

Tests are organized primarily by responsibility:

```text
tests/contract/          numerical and formulation contracts
tests/semantic/          semantic graph and resolution behavior
tests/compiler/          compiler planning and deal.II compiler capability
tests/dealii/            lower-level deal.II backend/service behavior
tests/application/       application, runner, and integration boundaries
tests/tools/             persisted-evidence and post-processing behavior
tests/support/           shared fixtures and algebraic reference models
```

The deal.II compiler suite intentionally uses one heavy translation unit:

```text
tests/compiler/dealii/dealii_compiler_contract.cc
        │
        └── scenarios/*.hpp
```

This keeps capability groups readable without recompiling the header-heavy compiler
and deal.II template surface in many translation units.

`tests/support/reference_models/` contains algebraic verification/oracle models. They
were deliberately removed from the public `include/nmopt/reference/` surface.

The scenario-discovery CMake support preserves fine-grained CTest scenario identity
even when several scenarios share one C++ translation unit.

## Build boundary

The project is header-heavy and exposes in-tree interface targets:

```text
nmopt_contract
nmopt_dealii_contract
```

`nmopt_contract` adds the project include tree and C++17 requirement.
`nmopt_dealii_contract` composes that target with the discovered deal.II package.

The current CMake project does **not** install/export an `nmopt` package for downstream
`find_package(nmopt)` use. Do not infer an installed-package workflow from the
interface-target names.

Only the exceptional Step-4 registration lives in
`cmake/ExternalStep4.cmake`. Ordinary project build/test registration remains in the
root `CMakeLists.txt`, while `build.sh` and `CMakePresets.json` define the supported
local profile workflow.

## Ownership and lifetime landmarks

The most important lifetime boundaries are:

- `ResolvedProblemView` borrows its source `ProblemSpec`;
- `DealiiDataBindings<dim>` borrows caller-owned deal.II function objects;
- `DealiiCompilationSession<dim>` exclusively owns a triangulation moved into it;
- compiled products retain the session lifetime token when the owned-session overload
  is used;
- detached reduced and supplied-OTD services retain that token;
- `CompiledApplicationViewT` retains its numerical model through callbacks, but its
  native-output callbacks may still borrow compile-time data bindings such as forcing
  and desired state;
- experiment envelopes own copied values and do not own executable services.

An owned compilation session therefore solves mesh/model lifetime. It does not
magically own every data object referenced by application callbacks.

## Evidence layers

Different layers retain different evidence. They should not be collapsed into one
generic result object.

```text
semantic/compiler producer
    ↓
CompilationManifest
    selected semantic/compiler decision and realized numerical services

solver/formulation
    ↓
reduced/KKT/PDAS result or report
    convergence, work, histories, and retained numerical evaluation

experiment
    ↓
ReducedExperimentEnvelopeT
    manifest + solver-policy snapshot + solver report + run environment

application harness
    ↓
BenchmarkArtifactT
    benchmark identity + diagnostics + envelope + measurements

repository runner
    ├──► artifact.kv / solver-trace.csv / native output
    └──► run-manifest.json
         run-set inventory, configuration provenance, and artifact outcomes
```

The envelope owns values and does not retain the executable compiler product. By
contrast, detached numerical services may retain a compilation-session lifetime token
when their callbacks need backend state. This is why numerical lifetime and persisted
evidence lifetime are documented separately.

## Where to make a change

| Goal | Start here | Usually also inspect |
| --- | --- | --- |
| Change semantic component structure | `include/nmopt/semantic/v1/types.hpp`, `validation.hpp` | `detail/validation_detail.hpp`, `problem_library.hpp`, semantic tests |
| Add a supported semantic family | `include/nmopt/semantic/v1/problem_library.hpp` | problem-authoring reference, semantic/compiler tests |
| Add a deal.II compiler capability | `include/nmopt/compiler/v1/dealii_types.hpp`, `dealii_compiler.hpp` | `dealii_scalar_plan.hpp` or the owning target header, compiler tests, compiler reference |
| Change compiled-product ownership | `compiled_products.hpp` | formulation contracts, manifest, lifetime tests |
| Change manifest/provenance | `compilation_manifest.hpp`, manifest projection in `dealii_compiler.hpp` | manifest test helpers and artifact consumers |
| Add a generic optimization method | `include/nmopt/solvers/` | optimization reference and contract tests |
| Adapt an existing PDE code directly | application code + `CallbackExecutableModelT` | external-integration reference and Step-4 minimal consumers |
| Add an `nmopt`-native application family | `include/nmopt/application/` | application-authoring reference and application tests |
| Register a repository benchmark | `apps/nmopt-runner/benchmark_registry.hpp` | parameter schema/binder, execution registration, parameters, runner tests |
| Change run-set/output policy | `apps/nmopt-runner/` | application-execution and parameter-file references |
| Change persisted-artifact interpretation | `tools/nmopt_artifacts/`, `tools/nmopt_postprocess/` | tool contracts and reporting/post-processing consumers |
| Change build/test registration | root `CMakeLists.txt`, `cmake/` | `.agents/build.md` and the relevant CTest inventory |

## Where to start reading

For a numerical-contract change:

```text
include/nmopt/contract/
        ↓
tests/contract/
        ↓
include/nmopt/dealii/ when a concrete deal.II service is involved
```

For a semantic/compiler feature:

```text
include/nmopt/semantic/v1/types.hpp
        ↓
validation.hpp / detail/validation_detail.hpp
        ↓
resolved_problem.hpp
        ↓
include/nmopt/compiler/v1/dealii_types.hpp
        ↓
dealii_scalar_plan.hpp or one target-specific realization
        ↓
dealii_compiler.hpp
        ↓
tests/compiler/
```

For Chapter-6 execution:

```text
include/nmopt/application/chapter6.hpp
        ↓
apps/nmopt-runner/parameter_files.hpp
        ↓
benchmark_binders.hpp
        ↓
chapter6_execution.hpp
        ↓
include/nmopt/application/dealii/chapter6_b1.hpp or chapter6_b2.hpp
        ↓
tests/application/
```

For an existing external PDE application:

```text
include/nmopt/contract/
        ↓
docs/reference/external-dealii-solver-integration.md
        ↓
apps/external-dealii/step-4/integration/
        ↓
apps/external-dealii/step-4/minimal/
```

The repository-structure decisions behind this current implementation map are preserved in the
[structure-refactor closure](../history/reviews/human-readability-audit/10-repository-structure-refactor-closure.md).
