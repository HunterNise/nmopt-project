# Repository scope — code-derived classification

**Audit baseline:** [`f53b7f009e5c`](https://github.com/HunterNise/nmopt-project/commit/f53b7f009e5c418ec4f3855db29f7eb924faf6f1)  
**Status:** complete for the documentation-refactor decision. This file records the implemented scope reconstructed from the audit baseline.

## 1. Why scope needs to be reconstructed

The repository currently contains:

- backend-neutral mathematical contracts;
- reduced optimization algorithms;
- alternate all-at-once/KKT/active-set products;
- a semantic problem-description graph;
- a substantial deal.II compiler and many registered numerical realizations;
- reusable deal.II numerical services;
- Chapter 5/6 recipe and scenario machinery;
- a headless benchmark runner and artifact layer;
- an external-application integration case study;
- dense/reference mathematical oracles;
- extensive tests and post-processing/reproduction tooling.

These are not all the same kind of product. A human reader should not have to infer their relative importance from the fact that many of them live under `include/nmopt/`.

The goal of this file is therefore to distinguish **capability**, **producer path**, **client/application infrastructure**, and **evidence machinery** before choosing the final project narrative.

## 2. Physical build boundary

The root build establishes a simple first cut.

[`CMakeLists.txt`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/CMakeLists.txt) defines:

```text
nmopt_contract
    INTERFACE
    public include tree
    C++17
    no deal.II dependency

nmopt_dealii_contract
    INTERFACE
    depends on nmopt_contract
    used only when deal.II is enabled
```

Backend-neutral contract, semantic, application-boundary, and runner-contract tests can therefore compile without deal.II.

The deal.II profile adds:

- `nmopt_runner`;
- deal.II compiler/backend tests;
- Chapter 6 deal.II application tests;
- external-integration tests.

The external Step-4 upstream, stripped, and adapted executables are deliberately created without the nmopt interface target. This is strong evidence that an application can exist independently and only its separate nmopt consumer needs the framework.

### Provisional implication

The reusable project cannot be described accurately as “the compiler/runner system.” Those are important paths, but the build graph proves that the numerical contracts and an external application integration have independent existence.

## 3. Public include areas

At the audit baseline the public tree has these top-level namespaces/directories:

```text
include/nmopt/
    application/
    compiler/
    contract/
    dealii/
    experiment/
    reference/
    semantic/
    solvers/
```

The resulting classification follows.

| Area | What the inspected code actually does | Provisional scope class |
| --- | --- | --- |
| `contract/` | Mathematical layouts, primal/covector roles, executable operations, metrics, constraints, solve evidence, reduced formulation, KKT/complementarity/PDAS/supplied-OTD products | **Core numerical interfaces/formulations** |
| `solvers/` | Reduced-space direction, line-search, trust-region, stopping and reporting algorithms over contract types | **Core optimization algorithms** |
| `semantic/v1/` | Backend-neutral problem graph, requirements/policies, structural validation and resolution | **Major optional authoring/front-end subsystem** |
| `compiler/v1/` | Turns a resolved semantic request plus concrete deal.II bindings/policies into executable numerical products, provenance, native views and solver services | **Major semantic-path producer subsystem** |
| `dealii/` | deal.II vector policy and reusable numerical realizations/services such as coordinate maps, metrics, constraints, SPD/KKT/PDAS solve helpers | **Backend/numerical realization support** |
| `application/` | Builds semantic graphs, binds scenarios/configuration, benchmark identity/evidence and headless orchestration; also contains Chapter-specific clients and deal.II adapters | **Mixed: generic client/orchestration + project-specific application code** |
| `experiment/` | Detached association of compilation manifest, solver policy/report and run environment | **Evidence/provenance support** |
| `reference/` | Dense/reference mathematical systems used as independent oracles and contract fixtures | **Verification/reference infrastructure** |

This table deliberately does **not** claim that eight human-facing architectural layers are needed.

## 4. Core mathematical/numerical capability

### 4.1 Representation and contracts

[`layout.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/contract/layout.hpp) establishes one of the project's strongest architectural choices:

```text
BlockLayout
    named spaces + dimensions

PrimalBlockT<Backend>
    values/directions/states/controls/seeds

CovectorBlockT<Backend>
    residuals/derivatives/pullbacks

pair(covector, primal)
    declared coefficient pairing
```

A derivative is not silently identified with a primal optimization vector. That conversion belongs to a metric.

[`metric_constraint.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/contract/metric_constraint.hpp) then separates:

- `MetricT`: primal-to-dual action and inverse Riesz-like action;
- `ConstraintT`: feasibility and metric-compatible projection.

This is not merely a deal.II implementation detail; it is part of the backend-neutral contract.

### 4.2 Residual/objective executable model

[`executable_model.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/contract/executable_model.hpp) supplies:

```text
E(x)
E'(x) dx
E'(x)^* p
J(x)
J'(x)
```

with explicit variable/test layouts.

[`callback_executable_model.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/contract/callback_executable_model.hpp) proves that this interface can be supplied directly by application callbacks. It contains no compiler or deal.II-specific problem-description logic.

### 4.3 Reduced formulation

[`reduced_dto.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/contract/reduced_dto.hpp) is the narrow implemented reduced formulation:

```text
control
  -> state solve
  -> objective

accepted/current value
  -> objective derivative
  -> adjoint solve
  -> residual VJP
  -> reduced control covector
```

It currently assumes:

- two variable blocks;
- one state block;
- one decision/control block;
- one residual test block.

This is an important scope limitation that should be visible to humans. It is a selected reduced formulation, not a universal PDE formulation object.

### 4.4 Reduced optimization

[`solvers/reduced_gradient.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/solvers/reduced_gradient.hpp) and adjacent solver headers own:

- direction policies;
- gradient/metric operations;
- constraints/projection;
- line-search trial construction and acceptance;
- trust-region/search policy;
- stopping;
- work/evidence records.

The implementation does not branch on PDE family or compiler target.

## 5. The project also has non-reduced products

The public contract layer is broader than the reduced path.

It includes:

- [`supplied_otd.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/contract/supplied_otd.hpp): supplied all-at-once state/adjoint/control systems;
- [`quadratic_kkt.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/contract/quadratic_kkt.hpp): equality-constrained quadratic KKT product;
- [`complementarity.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/contract/complementarity.hpp): metric-aware box complementarity and active-set representation;
- [`pdas.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/contract/pdas.hpp): PDAS active-set/KKT machinery.

The compiler mirrors this distinction. [`compiled_problem.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/compiler/v1/compiled_problem.hpp) contains separate products for:

```text
CompiledProblemT                reduced DTO path
CompiledQuadraticKKTProblemT    quadratic KKT path
CompiledPDASProblemT            PDAS path
CompiledSuppliedOTDProblemT     supplied OTD path
```

### Scope question still open

These capabilities are real and tested, but their *narrative centrality* is not yet established. They may be:

1. central evidence that nmopt supports multiple formulations; or
2. bounded research extensions around a primarily reduced-space project.

The final overview should not decide this from file count alone. Later runtime/application inspection must establish how much of the project actually uses these products.

## 6. Semantic authoring/front-end subsystem

[`semantic/v1/types.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/semantic/v1/types.hpp) defines `ProblemSpec` as a **composition root**, not a PDE model class.

It contains named:

- regions;
- spaces and pairings;
- variables and fixed data;
- transformations;
- residual terms and equations;
- observations and losses;
- metrics and constraints;
- requirement policies;
- formulation selection;
- optional supplied-OTD declaration.

The graph is backend-neutral in the strong software sense: it contains no deal.II objects.

However, it is **not purely continuous mathematics**. It also carries explicit backend-neutral discrete realization selections for difficult cases such as traces, negative metrics, transposition formulations, partial Dirichlet policies, and fractional trace metrics.

[`validation.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/semantic/v1/validation.hpp) makes the intended boundary explicit:

- semantic validation owns graph structure and stated policies;
- backend lowerability and formulation capability belong to the compiler.

### Provisional classification

The semantic system appears to be one substantial way to author a problem, not a prerequisite for the numerical solver core. The external Step-4 build separation corroborates that interpretation.

## 7. Compiler and deal.II realization

The compiler is a large subsystem.

The inspected code shows three different concerns that should not be collapsed in human documentation:

### 7.1 Compiler resolution/orchestration

[`dealii_compiler.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/compiler/v1/dealii_compiler.hpp) performs:

```text
semantic resolve
    -> closed compilation request
    -> lowerability/formulation/product validation
    -> registered numerical realization
```

### 7.2 Compiler product/provenance

[`compiled_problem.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/compiler/v1/compiled_problem.hpp) packages executable products and extensive typed compilation provenance.

Its `NativeApplicationViewT` is retained separately from `ExecutableModelT`, so dimensions/output/application information does not pollute the solver boundary.

### 7.3 Concrete deal.II numerical realization

`include/nmopt/dealii/` provides smaller reusable services such as:

- serial vector algebra policy;
- independent state coordinates;
- mass, negative-order and trace metrics;
- box constraints;
- serial SPD/KKT/PDAS solve helpers.

These are conceptually different from the compiler itself.

### Important implementation limitation

The current v1 compiler is not a fully generic independent-component lowerer, but it is also not purely a table of whole-problem implementations. [`dealii_scalar_plan.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/compiler/v1/dealii_scalar_plan.hpp) contains component-level residual, observation, loss, metric, constraint, and transformation handlers that contribute to a `ScalarLoweringPlan`.

The bounded scalar realization later closes that plan against explicitly supported compositions. For example, the residual assembly plan recognizes exact registered term sets for a diffusion/source/control case and for a tensor/transport/Robin case. Specialized compiler targets remain separate.

[`dealii_types.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/compiler/v1/dealii_types.hpp) defines a closed `ResolvedTargetFamily`, and `ResolvedCompilationRequest` records the chosen target family plus detailed realization facts.

Therefore a human-facing description should say something like:

> the semantic graph is compositional as a problem description; the current compiler performs component-level planning for a bounded scalar slice and then realizes only explicitly registered combinations and specialized targets.

This is more precise than either “fully compositional compiler” or “whole-problem switch statement.”

It should not imply arbitrary component recombination is already implemented.

## 8. “Backend” is narrower than “deal.II realization”

[`serial_backend.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/dealii/serial_backend.hpp) is just a vector/storage algebra policy:

```text
Vector
zeros
size
value/set_value
dot
add_scaled
scale
```

By contrast:

- [`independent_state_coordinates.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/dealii/independent_state_coordinates.hpp) realizes FE coordinate reconstruction/embedding/pullback;
- [`mass_metric.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/dealii/mass_metric.hpp) realizes an SPD Riesz map and inverse solve;
- [`serial_spd_solver.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/dealii/serial_spd_solver.hpp) realizes solve policy and returns backend-neutral solve evidence.

### Firm scope conclusion

The final docs should avoid using “backend” as shorthand for the entire deal.II PDE implementation. It is too narrow in the code.

## 9. Application and research infrastructure

The generic pieces of `application/` explicitly say what they do **not** own.

[`recipe.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/application/recipe.hpp):

> a recipe builds a semantic graph and owns no backend, mesh, solver, or experiment.

[`scenario.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/application/scenario.hpp):

> a scenario binds typed problem/compile/solver/experiment choices but does not compile or execute.

[`runner.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/application/runner.hpp):

> the runner owns orchestration; problem building and execution are supplied by callbacks; no compiler or optimizer is reproduced.

This strongly suggests an application/client layer **above** the numerical library.

At the same time, the same namespace contains Chapter 5/6 catalogues and deal.II adapters, and `application.hpp` aggregates both generic and Chapter-specific headers. This mixture will need a dedicated human-API/readability assessment.

## 10. Verification and evidence infrastructure

Several repository areas appear important for trust but not for the minimal runtime model:

- `reference/`: dense/reference mathematical oracles;
- `tests/contract`: backend-neutral behavior and mathematical invariants;
- `tests/semantic`: semantic validation/resolution;
- `tests/dealii`: compiler/numerical realization;
- `tests/application`: client/runner/external integration;
- `tools/`: reproduction, comparison, post-processing and report generation;
- `experiment/`: detached provenance/report association.

These should be visible in a project overview as **how the implementation is verified**, but they probably should not be presented as equal runtime layers.

## 11. Provisional project-scope statement

The implementation supports the following bounded description:

> `nmopt` is a C++17, backend-parametric set of mathematical contracts, formulations, and optimization algorithms for discretized PDE-constrained optimization, with a substantial deal.II realization/compiler path and an independent callback path for existing numerical applications. The repository additionally contains research-application assembly, benchmark/reproduction infrastructure, and alternate KKT/active-set formulation products.

This is deliberately broader and more precise than “a deal.II framework” while not yet claiming that every subsystem has equal status.

The audit will refine this statement after tracing real compiler, external, and Chapter 5/6 execution paths.
