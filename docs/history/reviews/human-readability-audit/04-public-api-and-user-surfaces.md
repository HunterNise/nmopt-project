# Public API and user surfaces

**Audit baseline:** [`f53b7f009e5c`](https://github.com/HunterNise/nmopt-project/commit/f53b7f009e5c418ec4f3855db29f7eb924faf6f1)  
**Status:** complete for the documentation-refactor decision. External application, semantic/compiler, algorithm/expert, and research/application surfaces are mapped; remaining questions are explicitly deferred API-maintainability questions.

## 1. There is not one user surface

The code exposes several ways to use the repository. Treating all public headers as one API would misrepresent the project.

At least four user profiles are already visible:

1. **existing numerical application user** — already owns a PDE implementation and wants nmopt formulation/optimization;
2. **semantic/compiler user** — wants to describe a supported problem and have nmopt realize it with deal.II;
3. **algorithm/expert user** — works directly with contracts, formulations, metrics, Hessians, KKT products, or solver policies;
4. **research/application runner user** — selects Chapter 5/6 scenarios, parameter matrices, artifacts and reproduction runs.

These profiles overlap, but they should not be forced to learn the same concepts.

## 2. External application surface

The minimal Step-4 Problem B consumer uses public nmopt types from:

```text
contract/
dealii/serial_backend.hpp
solvers/reduced_gradient.hpp
```

It does not use:

```text
semantic/
compiler/
application/
experiment/
reference/
nmopt-runner
```

### Concepts the binding currently requires

The explicit binding constructs or references:

```text
BlockLayout
SpaceId
PrimalBlockT / CovectorBlockT
CallbackExecutableModelT
StateControlPartitionT
StateAdjointSolversT
FormulationSolveResultT
LinearSolveReport
MetricT subclass
ReducedDTOT
SerialBackend
ReducedSearchSolverT
ReducedSolverParameters
```

### Which of these are mathematically meaningful?

Clearly meaningful to a user:

- state/control/test space identity;
- primal versus covector distinction;
- residual/objective actions;
- state and adjoint solves;
- control metric;
- reduced formulation;
- optimizer settings.

More representation-oriented:

- manually constructing `BlockLayout`;
- manually wrapping/unwrapping one-block vectors;
- translating native solve evidence into `LinearSolveReport`;
- explicit `StateControlPartitionT` construction for a fixed two-block convention;
- constructing the DTO service from several small containers.

This distinction is important. The audit should not count mathematically necessary concepts as “boilerplate” merely because they take lines of C++.

### Current interface-segregation observation

`CallbackExecutableModelT` requires:

```text
residual
residual JVP
residual VJP
objective
objective derivative
```

but the current first-order `ReducedDTOT` execution path calls only:

```text
objective
objective derivative
residual VJP
```

plus separately supplied state and adjoint solves.

Therefore residual/JVP are broader model capabilities rather than minimum reduced-first-order requirements.

This should be explained accurately to humans. Whether to add a narrower convenience surface is a later recommendation question.

## 3. Semantic/compiler surface

A semantic user deals with a different set of concepts:

```text
ProblemSpec
semantic component records
requirement/realization policies
SemanticValidator / resolver (directly or through compiler)
DealiiCompiler
DealiiCompilationSession
DealiiDataBindings
DealiiDiscretisationPolicy
CompilationProduct
CompilationResultT / CompiledProblemT
```

The compiler then supplies the downstream numerical services.

### Benefit

The user does not manually implement residual/JVP/VJP, state/adjoint solves, metric realization, coordinate reconstruction, or native field output for supported compiled targets.

### Cost

The semantic graph and compiler policy surface are substantial, and support is bounded rather than arbitrary. A human user needs a reliable way to answer:

- what can I describe semantically?
- what can the compiler actually realize?
- which discrete policy must I select explicitly?
- what concrete runtime data must I bind?
- which compilation product is available?

The current API may be technically precise but difficult to discover from headers alone because capability closure is spread across semantic enums, validation, capability registries, scalar plans, target families, and specialized realization headers.

## 4. Algorithm/expert surface

Expert users can work directly with:

- `MetricT`;
- `ConstraintT`;
- `ReducedHessianT`;
- reduced search direction and globalization policies;
- trust-region solver;
- quadratic KKT products;
- complementarity/PDAS;
- supplied OTD systems.

This is useful and demonstrates a real backend-neutral mathematical core.

However, the final human documentation should probably introduce this surface **after** a reader has the project mental model. Starting with the full contract directory would make the project appear substantially more complicated than the most common reduced use path.

## 5. Application/research surface

The generic application types are comparatively small:

```text
ProblemRecipeT
ScenarioT
HeadlessBenchmarkRunnerT
catalog / metadata / artifacts
```

But the same public namespace also contains:

- Chapter 5 recipes;
- Chapter 6 scenarios;
- Chapter 6 deal.II execution adapters;
- benchmark-specific data definitions.

And `application.hpp` aggregates them.

### Initial concern

There is a mismatch between **namespace/public-header placement** and **generality**:

```text
generic application orchestration
and
project-specific Chapter 5/6 research clients
```

are both presented under `include/nmopt/application/`.

This may be entirely acceptable for a research repository, but it makes “what is the reusable library?” less obvious to a human.

No move/split recommendation is made yet.

## 6. Umbrella-header observations

### `semantic/v1/problem_spec.hpp`

This compatibility aggregate includes not just core semantic types but also:

- reference specs;
- resolution;
- validation.

The comment says new users may include only what they need.

A future usability question is whether this file is intended as the obvious user entry point. If so, its inclusion of a very large reference-spec header may make source discovery noisier than necessary.

### `application/application.hpp`

This aggregate includes generic application interfaces and Chapter-specific content.

That is convenient for the project runner but may be a poor “start here” public header for external users.

## 7. Human API documentation should probably be path-oriented

Instead of one exhaustive “public API” document, humans likely need entry points such as:

```text
I already have a deal.II application
    -> external integration guide/reference

I want nmopt to construct a supported PDE optimization problem
    -> semantic/compiler guide/reference

I want to implement/change an optimization algorithm
    -> formulation/solver contracts

I want to reproduce the project's Chapter 6 experiments
    -> application/runner guide
```

The audit will test this against the remaining runtime and documentation structure.

## 8. Ergonomics findings

| Observation | Classification now |
| --- | --- |
| External binding must express mathematical spaces/actions/solves/metric | intrinsic or justified public concepts |
| Layout/block wrapping is explicit | possible construction ceremony |
| State/adjoint solve reports are explicit | justified evidence contract, but translation may be ceremony |
| Metric requires subclassing rather than callback bundle | possible ergonomics asymmetry |
| Metric inverse exposes no solve report | current contract limitation/choice; not yet a redesign request |
| Reduced formulation consumes less than full executable model requires | interface-segregation friction |
| Semantic compiler capability is distributed across several records/registries | discoverability/documentation issue, possibly source organization |
| Generic and Chapter-specific application APIs share an aggregate surface | human-boundary clarity issue |

## 9. Deferred API/refactor questions

These questions are not blockers for the documentation refactor. Answering them well would require a dedicated API/source-maintainability audit rather than inference from documentation friction alone.

- Which headers should a human user be told to include for each path?
- Is there a stable “supported public API” concept in the repository today, or merely public file placement?
- How much source inspection is required to construct a new semantic problem without copying an existing recipe?
- Are the exact lifecycle constraints of borrowed callbacks/views sufficiently visible at declaration sites?
- Does any repeated external construction justify a helper, or is documentation/example extraction enough?
- Should Chapter-specific application APIs remain in `include/nmopt/`?
