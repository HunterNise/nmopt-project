# Architecture map — derived from code

**Audit baseline:** [`f53b7f009e5c`](https://github.com/HunterNise/nmopt-project/commit/f53b7f009e5c418ec4f3855db29f7eb924faf6f1)  
**Status:** complete code-derived architecture synthesis. The human-facing documentation structure is proposed separately in `07-human-information-architecture.md`.

## 1. Method

This map is derived from responsibility, dependency direction, object ownership, and actual call surfaces. Directory names are supporting evidence, not the decomposition criterion.

The central question is:

> Which pieces need to understand which other pieces for a computation to happen?

## 2. First-pass dependency picture

The current code supports this structure:

```text
                         AUTHORING / PRODUCER PATHS

 semantic ProblemSpec
        |
        v
 semantic validation
        |
        v
 deal.II compiler + closed lowering decision
        |
        v
 typed numerical realization ---------------------+
                                                   |
 existing numerical application                    |
        |                                          |
        +-- application-owned operations ----------+
                                                   v
                                      MATHEMATICAL / FORMULATION PORTS
                                      layouts, E/J actions, solves,
                                      metric, optional capabilities
                                                   |
                         +-------------------------+----------------------+
                         |                                                |
                         v                                                v
                  reduced formulation                           alternate products
                  state/adjoint DTO                    supplied OTD / KKT / PDAS
                         |
                         v
                 reduced optimization
             direction/globalization/stopping
                         |
                         v
                      result

 APPLICATION / RESEARCH CLIENTS
 recipes, scenarios, runner, manifests/artifacts, native output,
 Chapter 5/6 execution and post-processing sit around selected producer
 and solver paths; they are not required by the external application path.
```

A separate representation concern cuts horizontally across numerical layers:

```text
Backend policy
    native vector storage + primitive algebra only
```

The deal.II compiler and FE numerical services are **not** equivalent to the backend policy.

## 3. Responsibility A — mathematical representation and invariant checking

### Primary code

- [`contract/layout.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/contract/layout.hpp)
- [`contract/linalg.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/contract/linalg.hpp)

### Owns

- runtime space identity;
- block dimensions;
- primal versus covector role;
- backend-generic block storage;
- coefficient pairing;
- checked algebra between compatible layouts.

### Does not own

- PDE meaning;
- FE spaces/meshes;
- objective semantics;
- state/adjoint solve policy;
- optimization iteration policy.

### Why this is probably human-visible

Many later design choices only make sense once the reader understands that a derivative/covector is not automatically a primal vector. This is a mathematical convention, not incidental type ceremony.

## 4. Responsibility B — generic mathematical capability contracts

### Primary code

- [`contract/executable_model.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/contract/executable_model.hpp)
- [`contract/metric_constraint.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/contract/metric_constraint.hpp)
- [`contract/linear_solve.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/contract/linear_solve.hpp)
- [`contract/reduced_hessian.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/contract/reduced_hessian.hpp)

### Owns

Capabilities algorithms may request:

```text
ExecutableModel:
    E(x)
    E'(x) dx
    E'(x)^* p
    J(x)
    J'(x)

Metric:
    G v
    G^-1 xi

Constraint:
    feasibility
    projection in a compatible metric

Solve result:
    primal solution
    convergence/work evidence
```

### Does not own

- how a PDE is assembled;
- how a state or adjoint linear system is inverted;
- compiler provenance;
- application output.

### Important nuance

`ExecutableModelT` is a **complete model capability interface**, not the minimum interface of every formulation. The reduced first-order formulation currently uses only a subset of it.

This is one likely place where the project is mathematically clean but external authoring feels verbose.

## 5. Responsibility C — formulations

The code contains more than one formulation/product family.

### C1. Reduced discretize-then-optimize formulation

Primary code: [`contract/reduced_dto.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/contract/reduced_dto.hpp)

Owns the orchestration:

```text
u
 -> solve state y
 -> J(y,u)

(y,u)
 -> J'(y,u)
 -> solve adjoint p from state derivative
 -> E'(y,u)^* p
 -> j'(u) = J_u - E_u^* p
```

The formulation owns **when** these operations are required, not **how** they are numerically implemented.

It currently assumes exactly:

- one state block;
- one decision/control block;
- one residual test block.

The code also deliberately separates value evaluation from derivative augmentation so a line search can evaluate rejected trials without adjoint solves.

### C2. Supplied all-at-once OTD

Primary code: [`contract/supplied_otd.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/contract/supplied_otd.hpp)

Owns a three-block supplied state/adjoint/control system with residual/JVP/VJP/solve actions.

It is a distinct product, not an extension of `ReducedDTOT`.

### C3. Quadratic KKT and PDAS

Primary code:

- [`contract/quadratic_kkt.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/contract/quadratic_kkt.hpp)
- [`contract/complementarity.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/contract/complementarity.hpp)
- [`contract/pdas.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/contract/pdas.hpp)

These represent:

- equality-constrained quadratic products;
- explicit domain/range pairings;
- multiplier/adjoint conversion;
- symmetry/rank/kernel assumptions;
- metric-aware box complementarity;
- active-set restriction and PDAS iteration machinery.

### Human-facing question

It is not yet clear whether “formulation” should be one overview layer containing these branches, or whether the final narrative should introduce the reduced path as the core and place KKT/PDAS/OTD under “additional implemented formulations.”

## 6. Responsibility D — optimization algorithms

### Primary code

- [`solvers/reduced_gradient.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/solvers/reduced_gradient.hpp)
- `solvers/reduced_line_search.hpp`
- `solvers/reduced_search.hpp`
- `solvers/reduced_trust_region.hpp`

### Owns

- search direction policy;
- metric-based gradient conversion/norms;
- trial construction;
- projection through an optional constraint;
- Armijo/Wolfe/exact or other line-search policy;
- trust-region/search policy;
- stopping;
- detailed work and acceptance records.

### Dependency direction

The inspected reduced solver consumes:

```text
ReducedDTOT
MetricT
optional ConstraintT
direction policy
line-search policy
```

It has no semantic/compiler/deal.II/PDE-family dependency.

### Architectural consequence

The optimizer appears to be one of the cleanest reusable layers in the project. Its human documentation should be able to explain it almost entirely from the reduced mathematical problem.

## 7. Responsibility E — semantic problem description

### Primary code

- [`semantic/v1/types.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/semantic/v1/types.hpp)
- [`semantic/v1/resolved_problem.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/semantic/v1/resolved_problem.hpp)
- [`semantic/v1/validation.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/semantic/v1/validation.hpp)

### Owns

Backend-neutral declarations of:

- regions;
- spaces and pairings;
- variables/data;
- transformations;
- residual terms/equations;
- observations/losses;
- metric/constraint declarations;
- requirement policies;
- formulation selection.

### Does not own

- deal.II objects;
- sparse matrices;
- concrete solve execution;
- optimizer iteration.

### Important qualification

The semantic graph contains explicit **discrete realization policy selections** where the mathematics does not determine a unique numerical treatment.

It is therefore better described as:

> backend-neutral problem semantics plus explicit non-inferable realization commitments

than as a purely continuous PDE DSL.

### Semantic validator boundary

`SemanticValidator` validates graph structure and stated policies. Compiler lowerability and formulation capability are intentionally left to the compiler.

## 8. Responsibility F — compiler and compilation products

This responsibility is internally heterogeneous.

### F1. Closed request resolution

[`compiler/v1/dealii_compiler.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/compiler/v1/dealii_compiler.hpp) resolves a semantic graph into a compiler-owned request, then validates lowerability and product capability.

### F2. Component planning plus registered realization closure

[`compiler/v1/dealii_scalar_plan.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/compiler/v1/dealii_scalar_plan.hpp) introduces a bounded component-planning layer. Residual terms, observations, losses, metrics, constraints, and transformations have handlers that contribute typed records to a `ScalarLoweringPlan`.

That plan does not imply arbitrary recombination. Its residual-assembly slice still recognizes exact supported registrations, and specialized targets live outside this common scalar path.

[`compiler/v1/dealii_types.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/compiler/v1/dealii_types.hpp) exposes `ResolvedTargetFamily` and a large `ResolvedCompilationRequest`.

Together these files show a hybrid design:

```text
validated semantic components
    -> bounded component planning
    -> closed registration / target-family decision
    -> concrete deal.II realization
```

Current lowering therefore has real local component structure without claiming arbitrary independent component composition.

### F3. Product/provenance packaging

[`compiler/v1/compiled_problem.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/compiler/v1/compiled_problem.hpp) packages:

- executable mathematical ports;
- metric/constraint/solve services;
- optional Hessian;
- lifetime owner;
- typed compilation manifest;
- optional native application view.

Separate product classes exist for reduced, KKT, PDAS, and supplied-OTD results.

### F4. Native application seam

[`compiler/v1/native_application_view.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/compiler/v1/native_application_view.hpp) retains:

- physical/independent dimensions;
- objective components;
- native output callback.

This data deliberately stays outside `ExecutableModelT`.

### Architectural consequence

The compiler is best understood as a **producer of numerical services and typed provenance**, not as the universal object model consumed by algorithms.

## 9. Responsibility G — backend representation versus numerical realization

This needs especially careful wording.

### G1. Backend policy

[`dealii/serial_backend.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/dealii/serial_backend.hpp)

Owns only:

- concrete native vector type;
- size/value access;
- zero construction;
- dot;
- add/scale.

It makes generic contract/solver templates work over `dealii::Vector<double>`.

### G2. FE/numerical services

Examples:

- [`dealii/independent_state_coordinates.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/dealii/independent_state_coordinates.hpp)
- [`dealii/mass_metric.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/dealii/mass_metric.hpp)
- [`dealii/serial_spd_solver.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/dealii/serial_spd_solver.hpp)

These own real numerical behavior such as:

- state reconstruction `P z + ell`;
- tangent embedding `P dz`;
- covector pullback `P^T q`;
- sparse metric application;
- iterative inverse metric action;
- convergence policy/reporting.

### Firm distinction

```text
backend policy != FE discretization != compiler
```

A human overview should preserve those distinctions while avoiding an explosion into three equal “layers” if a simpler explanation suffices.

## 10. Responsibility H — application/client orchestration

### Generic surfaces

- [`application/recipe.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/application/recipe.hpp)
- [`application/scenario.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/application/scenario.hpp)
- [`application/runner.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/application/runner.hpp)

These explicitly avoid owning PDE solves or optimization algorithms.

A useful dependency picture is:

```text
recipe parameters -> ProblemSpec

scenario
    = problem parameters
    + compile options
    + solver options
    + experiment options

headless benchmark runner
    calls supplied problem builder
    calls supplied execution adapter
    packages evidence/artifact
```

### Project-specific content

The same namespace also contains Chapter 5/6 definitions and deal.II execution adapters. This makes the namespace useful to the research project but potentially less clear as a generic library API.

This will be revisited in the public-API and source-readability phases.

## 11. Responsibility I — evidence, verification and research tooling

### `experiment/`

[`experiment/reduced_envelope.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/experiment/reduced_envelope.hpp) associates values only:

```text
CompilationManifest
solver-policy snapshot
detached solver report
run-environment record
```

It explicitly does not retain executable services.

### `reference/`

Dense/reference models and systems serve as independent mathematical oracles.

### Tests

The top-level test split (`contract`, `semantic`, `dealii`, `application`, `tools`) mirrors evidence categories reasonably well, even if it should not be copied literally into an architecture diagram.

### Tools

Python/shell tooling handles persisted run processing, comparisons, reports, and external-integration checks.

### Human-facing treatment

These areas are critical for **why we trust the implementation**, but most are not runtime layers a library user should have to learn first.

## 12. Two producer paths already visible from code

Even before tracing full runtime calls, the code supports this important convergence:

```text
SEMANTIC PATH

ProblemSpec
   -> semantic validation/resolution
   -> deal.II compiler
   -> compiled mathematical services
   -> ReducedDTOT / other compiled product


EXTERNAL APPLICATION PATH

existing application
   -> CallbackExecutableModel + solve callbacks + metric
   -> ReducedDTOT


COMMON REDUCED CONSUMER

ReducedDTOT
   -> ReducedSearchSolver / ReducedTrustRegionSolver
```

The CMake treatment of Step-4 independently supports this interpretation.

## 13. Initial architectural tensions to audit further

These are observations, not recommendations.

### 13.1 Complete-model contract versus formulation-specific requirements

`ReducedDTOT` does not call `residual()` or `residual_jvp()`, yet an external callback model must provide them because `ExecutableModelT` requires all five operations.

This is real interface-segregation friction, but the external experiment demonstrated that it does not force cross-layer ownership or compiler adoption.

### 13.2 Header size and conceptual locality

Several public semantic/compiler headers are very large:

- `semantic/v1/validation.hpp` ~178 KB;
- `compiler/v1/dealii_compiler.hpp` ~341 KB;
- `compiler/v1/dealii_fixed_dirichlet.hpp` ~96 KB;
- numerous target realizations ~30–50 KB.

This can obstruct human source navigation even if the dependency architecture is sound.

The audit must determine whether the appropriate remedy is:

- better overview/source maps;
- local comments;
- file splitting/reorganization;
- or no code change.

### 13.3 Generic and project-specific application code share one public namespace

`application.hpp` aggregates generic recipe/scenario/runner facilities together with Chapter-specific code.

This may be harmless for the research repository, or it may obscure the reusable boundary for external users. Later phases will judge this from actual use.

### 13.4 Multiple formulation families complicate the “one sentence” project identity

Reduced DTO, supplied OTD, quadratic KKT and PDAS are all real code.

The final overview needs a principled way to explain their relationship without making a new reader understand every formulation before understanding the project.

## 14. Working minimal mental model

At this stage the smallest accurate model appears to be:

```text
1. Describe or supply mathematics.
   - semantic graph + compiler, OR
   - existing numerical application callbacks

2. Realize numerical capabilities.
   - vectors/layouts
   - PDE/objective actions
   - state/adjoint or all-at-once solves
   - metric/constraints

3. Choose a formulation/product.
   - reduced DTO
   - selected alternate KKT/OTD/PDAS products

4. Run a numerical algorithm.
   - reduced search/trust region or corresponding product solver

5. Optionally wrap this in project/application infrastructure.
   - recipes/scenarios
   - runner
   - artifacts/provenance
   - post-processing/reproduction
```

This is a working audit model only. It will be tested against the compiler, Step-4, and Chapter 5/6 runtime paths before becoming documentation architecture.
