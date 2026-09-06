# PDE–solver boundary refactor assessment

## Assessment status

**Complete architecture audit; implementation not started.**

This document records the evidence and conclusions used to define the
PDE–solver boundary refactor. It is an audit archive, not the mutable work
queue. The [roadmap](roadmap.md) owns execution order and status; the
[boundary design](../../../design/pde-solver-boundary.md) owns the accepted
long-lived architectural direction.

The audit baseline is the pre-refactor application milestone at
`codex/applications` commit `e8a7975db69a4d13ea8a32c3a94ac2c7295486c6`.
The refactor proceeds on a separate branch from that milestone.

## Executive verdict

The current repository already has strong mathematical and formulation
contracts. The principal problem is not the generic solver API. It is the
native numerical middle layer and compiler composition around it.

The semantic layer describes a compositional graph, and the scalar lowering
plan preserves part of that decomposition. Concrete lowering then tends to
collapse the result into a complete target-specific `*Model` that owns FE
spaces, operators, state/adjoint solves, objective pieces, metrics,
constraints, optional Hessian operations, diagnostics, and output. The compiler
then erases that rich typed object to `ExecutableModelT`, while native
application code later recovers the concrete type with `dynamic_cast` to reach
output and diagnostic operations.

The central corrective direction is therefore:

```text
semantic declaration
        |
        v
closed lowering decision
        |
        v
typed numerical realization
       / \
      /   \
solver-facing ports   typed/native application view
      |                         |
      v                         v
formulation/optimizer      output/diagnostics
```

The typed numerical realization should be composed from independently varying
responsibilities where the current implementations demonstrate real reuse. It
should not be replaced by a new universal PDE hierarchy.

The refactor should be deletion-first. Preserve small contracts and focused
components; migrate unique mathematical behavior and tests; delete duplicate
target identities, repeated service wiring, premature type-erasure recovery,
model-owned output, compatibility state, and redundant complete
implementations when their consumers have moved.

## Audit method

The audit used a repeated `map -> descend -> trace -> summarize -> update map`
method.

For each subsystem or numerical feature it asked:

- what mathematical or architectural object it represents;
- how it is implemented and constructed;
- who owns it and for how long;
- which non-test consumers use it;
- which exact types cross the boundary;
- where policy and mutable state live;
- which tests protect it; and
- whether it should survive, simplify, move, or disappear.

Three views were maintained throughout:

1. dependency map;
2. runtime call graph; and
3. ownership graph.

The audit traced both a native semantic/compiler execution and the intended
external-application integration path. It also compared sibling numerical
realizations across distributed control, Neumann control, Dirichlet control,
coefficient identification, observations, metrics, and second-order
capabilities.

## Current repository architecture

The detailed navigational version of these maps lives in the
[architecture map](architecture-map.md). This assessment keeps only the
evidence needed to support the findings below; the map is the preferred
starting point for reconstructing the whole system.

### Layer map

| Layer | Representative location | Current useful responsibility | Main pressure observed |
| --- | --- | --- | --- |
| Mathematical contracts | `include/nmopt/contract/` | layouts, primal/covector roles, executable operations, reduced DTO, metrics, constraints, KKT/PDAS products | mostly sound; avoid broadening |
| Optimization algorithms | `include/nmopt/solvers/` | search direction, line search, trust region, convergence/reporting | sound separation from PDE details |
| deal.II backend services | `include/nmopt/dealii/` | backend vector policy, metrics, constraints, linear/KKT/PDAS helpers | direct v0 model and some mixed responsibilities remain |
| Semantic graph | `include/nmopt/semantic/v1/` | backend-neutral problem graph, validation, policies, resolved lookup view | large files, mutation-heavy reference builders, some ad-hoc sum-type encoding |
| Compiler/lowering | `include/nmopt/compiler/v1/` | runtime bindings, lowerability, scalar plan, concrete realization, compiled products, manifest | repeated target identity, complete-model dispatch, early erasure |
| Application assembly | `include/nmopt/application/` | recipes, scenarios, harness records, backend-neutral configuration | generally healthy; some duplicated projection/configuration facts |
| deal.II application adapters | `include/nmopt/application/dealii/` | mesh/data binding, compile, optimize, diagnostics/output for B1/B2 | overloaded composition roots; concrete-model downcasts |
| Experiment records | `include/nmopt/experiment/` | detached results/provenance envelope | must avoid retaining live numerical/compiler state |
| Runner | `apps/nmopt-runner/` | CLI, parameter binding, run-set/filesystem orchestration | large but mostly outer-layer concern |
| Persisted tools | `tools/` | artifact/native-output reading and reporting | should remain downstream of stable persisted data |

### Compile-time dependency sketch

```text
contract -> solvers ------------------------------\
    |                                               \
    +-> deal.II services --\                         +-> application
                           +-> compiler ------------/        |
semantic -----------------/        |                         v
    |                              +-> experiment -------> runner
    +-----------------------------------------------------/   |
                                                              v
                                                   persisted output -> tools
```

This direction is broadly healthy. The main architectural problem is not an
obvious inverted include dependency; it is the concentration of unrelated
runtime responsibilities inside concrete numerical targets and the repeated
representation of compiler decisions.

## Runtime call graph

A representative native reduced run currently follows this shape:

```text
CLI / parameter file
    -> runner registration and run-set selection
    -> Scenario
    -> ProblemRecipe
    -> ProblemSpec
    -> runtime deal.II functions + mesh/session
    -> DealiiCompiler
    -> complete concrete target model
    -> CompiledProblemT
    -> ReducedDTOT
    -> reduced optimizer
         -> evaluate_value(control)
              -> solve_state callback
              -> objective
         -> augment_derivative(...)
              -> objective_derivative
              -> solve_adjoint callback
              -> residual_vjp
         -> metric.inverse_apply
         -> line search / trust-region policy
    -> application adapter
         -> recover typed model in some paths
         -> diagnostics/native output
    -> experiment envelope
    -> artifact writer / tools
```

At the formulation boundary the important runtime path is already clean:

```text
solve(u0)
  -> ReducedDTOT::evaluate(...)
       -> solve_state callback
       -> objective
       -> objective_derivative
       -> solve_adjoint callback
       -> residual_vjp
  -> direction policy
       -> metric.inverse_apply
  -> globalization
       -> evaluate_value / augment_derivative
```

The optimizer sees no deal.II mesh, compiler target, DoF handler, or PDE-family
branch.

## Ownership graph

The current native path has several sound lifetime decisions:

- compilation sessions may own a triangulation;
- DoF handlers and concrete model objects depend on that mesh lifetime;
- `CompiledProblemT` owns or retains executable and capability objects;
- model-capturing solve lambdas keep concrete numerical state alive;
- `ReducedDTOT` can retain an opaque lifetime owner for the compiler path.

The problematic ownership transition is conceptual rather than memory-unsafe:

```text
rich concrete numerical object
        |
        v
erased ExecutableModelT
        |
        +--> generic formulation/optimizer       [correct consumer]
        |
        +--> native application needs rich data
                 |
                 v
             dynamic_cast back                  [boundary too early]
```

The target keeps the typed numerical realization alive explicitly for native
application consumers while separately exposing erased solver ports.

## What is already good

### Generic executable and formulation contracts

`ExecutableModelT<Backend>` is intentionally small: variable/test layouts,
residual, residual JVP, residual VJP, objective, and objective derivative. It
expresses the nonlinear mathematical operations required by formulations
without assuming a fixed matrix form.

This is validated by coefficient identification, where
$E(y,m) = A(m)y - f$; a universal $A/B/Q/R$ interface would not be sufficient.

`ReducedDTOT` correctly means **Discretize Then Optimize**. It owns the
first-order reduced orchestration rather than the PDE implementation. Its split
value/derivative evaluation supports globalization without unnecessary adjoint
work.

`StateAdjointSolversT` is already the desired inversion-of-control boundary for
state and adjoint solves.

### Backend-neutral metric and constraint contracts

`MetricT` correctly separates a decision covector from a primal search
direction. `ConstraintT` exposes feasibility/projection behavior without
requiring the optimizer to know the discretization.

The concrete `Hminus1Metric` and `TraceHhalfMetric` demonstrate that metric
application may intrinsically require linear solves without justifying a
universal linear-solver interface.

### Generic KKT product

`EqualityConstrainedQuadraticKKTProductT<Backend>` expresses KKT actions through
callbacks and remains backend-neutral. The serial deal.II packing/solve adapter
is a good example of a typed backend adapter around a generic mathematical
product.

### Semantic representation and validation

`ProblemSpec` is a component graph rather than a PDE-family enum. Residual
terms, observations, losses, metrics, constraints, transformations, and
policies are represented independently. The semantic validator checks graph
meaning and stated policies rather than hard-coding named complete problem IDs.

`ResolvedProblemView` is a useful validated/indexed borrowing view. It resolves
stable IDs for one validation/compilation operation without performing
numerical lowering or target dispatch.

### Focused numerical components already exist

Several current implementation pieces already point toward the desired typed
middle layer:

- `VolumeObservationAssembly`;
- `NeumannControlRealisation` with facewise and continuous implementations;
- backend metric/constraint objects;
- serial SPD/KKT helpers;
- explicit state reconstruction logic in `ScalarComponentModel` and
  `DirichletControlLiftingModel`.

The refactor should reuse and narrow demonstrated components rather than design
a new hierarchy from scratch.

### Recipes, scenarios, and experiment envelope

`ProblemRecipeT` remains a small semantic builder without backend/solver state.
`ScenarioT` is a composition/configuration record rather than an execution
engine. `ReducedExperimentEnvelopeT` detaches recorded experiment evidence from
live compiled objects. These are useful boundaries to preserve.

## Findings

### PSR-001 — The solver-facing contracts are not the central refactor problem

**Observed.** `ExecutableModelT`, `StateAdjointSolversT`, `MetricT`, optional
`ConstraintT`, optional `ReducedHessianT`, and `ReducedDTOT` already separate
generic formulations/optimizers from deal.II details.

**Good.** They are small, mathematically meaningful, and exercised by multiple
concrete targets.

**Direction.** Preserve them. Do not replace them with a universal PDE model,
linear solver, residual hierarchy, observation hierarchy, or output hierarchy.

### PSR-002 — Existing applications lack a callback-backed executable bridge

**Observed.** State/adjoint inversion is already callback-based, but an
application that does not derive from `ExecutableModelT` has no equally small
bridge for residual/objective operations.

**Impact.** The desired SUNDIALS-like inversion-of-control path is incomplete:
an existing PDE application should not be required to adopt `ProblemSpec` or
inherit an nmopt PDE interface.

**Direction.** Add one small backend-neutral callable-to-`ExecutableModelT`
bridge and prove it with a tutorial-style deal.II application before relying on
the boundary for larger native cleanup.

### PSR-003 — Complete concrete `*Model` types are numerical composition roots

**Observed.** `ScalarDiffusionReactionModel`, `ScalarComponentModel`,
`ContinuousControlModel`, `NeumannBoundaryControlModel`,
`DirichletControlLiftingModel`, and `CoefficientIdentificationModel` combine
several independently varying concerns.

Typical responsibilities include FE/DoF state, constraints, operator assembly,
decision coordinates/coupling, observations, objective regularization,
state/adjoint solves, metric/constraint factories, optional Hessian actions,
diagnostics, and native output.

**Impact.** New semantic combinations tend to require a new complete target or
new target-specific branches even when many numerical pieces are identical.

**Direction.** Treat complete models as migration containers, not target
architecture. Extract only demonstrated reusable responsibilities and delete
or narrow the containing complete types when their consumers have moved.

### PSR-004 — Semantic composition survives only partway through lowering

**Observed.** The semantic graph represents residual terms independently.
`ScalarLoweringPlan` also stores component contributions such as
`ScalarResidualContribution`, observations, losses, metric, constraint, and
transformation information.

`ScalarComponentModel` then becomes the complete numerical target and consumes
the plan as one large implementation object.

**Impact.** The semantic/compiler representation advertises more component
structure than the typed numerical layer retains.

**Direction.** Preserve a closed component-aware lowering decision farther into
typed numerical construction without requiring arbitrary graph composition or
one C++ class per semantic node.

### PSR-005 — Volume and Neumann control expose complementary halves of the target design

**Observed.** Volume control is componentized in the scalar compiler plan but
its concrete matrix is owned directly by `ScalarComponentModel`. Neumann
control bypasses that component planner and triggers a whole target family, but
its typed numerics are more modular through `NeumannControlRealisation` and
`coupling_action` / `coupling_transpose_action`.

**Impact.** The repository already contains both useful design halves, but in
different targets.

**Direction.** Move toward component-aware compiler decisions plus narrow typed
numerical operations. Do not infer that all controls require one public
`DecisionRealisation` hierarchy.

### PSR-006 — Decision variables have three different PDE-influence mechanisms

**Observed.** Current targets demonstrate:

1. affine forcing controls: distributed and Neumann control;
2. state-reconstruction control: Dirichlet lifting; and
3. operator parameters: coefficient identification.

**Impact.** A universal fixed coupling matrix abstraction would fail for
Dirichlet lifting and coefficient-dependent operators.

**Direction.** Keep the universal solver contract at $E$, $E'v$, $E'^*p$.
Compiler plans may use a closed decision-realization variant; fixed-linear
coupling helpers are internal specializations only where they delete code.

### PSR-007 — State coordinates are independent of the PDE operator

**Observed.** Similar scalar PDEs use different state-coordinate conventions:
the direct v0 and continuous-control paths use full FE coordinates with
constrained/identity rows, while `ScalarComponentModel` uses independent
coordinates and an affine reconstruction $P \hat{y} + \ell$.
`DirichletControlLiftingModel` extends the same mechanism with $L_D u$.

**Impact.** Boundary/state-coordinate policy is duplicated inside complete
models and cannot be described as part of one PDE-family class.

**Direction.** Treat state-coordinate/reconstruction machinery as an
independent typed numerical responsibility. A concrete reusable component is
justified if migration of at least two current consumers deletes duplicated
logic.

### PSR-008 — Operator structure, coordinate policy, gauge, and solve policy are separate axes

**Observed.** Fixed symmetric diffusion, fixed nonsymmetric transport,
mean-zero pure Neumann, and parameter-dependent diffusion use different solve
requirements while sometimes sharing state FE/coordinate machinery.

**Impact.** Complete models currently mix mathematical operator facts with
concrete solver selection.

**Direction.** Typed numerics may expose symmetry/nullspace/point-dependence as
facts. State/adjoint service construction chooses the concrete solve policy.
Do not add a universal public linear-solver API.

### PSR-009 — Observation maps and compiled tracking forms have different consumers

**Observed.** Volume and $H^{1}$ tracking are naturally compiled to quadratic data
$Q, q, c$. Point sensors, boundary traces, and normal flux additionally retain
explicit measurement value/JVP/VJP data for verification and diagnostics.

**Impact.** Treating every observation as either only a quadratic form or only
an explicit sampled vector would lose useful structure or add unnecessary
materialization.

**Direction.** Keep observation semantics separate from loss semantics. Retain
explicit typed observation ports only where another consumer needs them;
compile efficient objective forms for solver use.

### PSR-010 — Objective regularization, metric, and reduced Hessian are distinct roles

**Observed.** Continuous and Dirichlet control permit independent regularization
norm and search metric choices. The same mass/stiffness/trace operator may be
used by both. `ReducedHessianT` composes PDE sensitivity, observation curvature,
and regularization curvature.

**Impact.** Current classes sometimes make shared numerical data look like one
semantic responsibility. `TraceHhalfMetric`, for example, can also be used as
the operator implementing an $H^{1/2}$ objective norm.

**Direction.** Share immutable numerical ingredients, not semantic ownership.
Keep regularization, metric, and optional reduced Hessian separately composed.
A neutral operator extraction is optional and must prove net simplification.

### PSR-011 — `NeumannControlRealisation` is useful but too broad

**Observed.** The abstraction legitimately supports facewise and continuous
runtime realizations and provides coupling/transpose actions. It also owns or
provides regularization, metric creation, and native output. The containing
model downcasts it to the facewise concrete type for a coefficientwise box
constraint.

**Impact.** One useful polymorphic component has accumulated unrelated
consumers because they share the same boundary discretization data.

**Direction.** Preserve demonstrated control-realization polymorphism but
review whether its stable core should be decision coordinates/coupling and
shared numerical operators. Metric, constraint, regularization, and output
should be capabilities/consumers rather than promises of every realization.

### PSR-012 — Type erasure occurs before native application consumers are finished

**Observed.** The compiler constructs a rich concrete model, stores it as
`ExecutableModelT`, then B1/B2 application adapters use `dynamic_cast` to
recover complete target types for dimensions, objective components, and native
output. Compiler-side manifest utilities also erase to the executable surface
and later rediscover target-specific facts.

**Impact.** The solver boundary becomes an accidental application ownership
boundary and creates explicit erase-then-recover cycles.

**Direction.** Compiler/native construction should return or retain both a
solver-facing erased view and a typed/native realization or closure. Do not add
output/dimension methods to `ExecutableModelT`.

### PSR-013 — Compiler target identity is represented several times

**Observed.** Target family, private compiled target kind, `uses_*` booleans,
predicates, construction branches, and string/provenance labels re-express
substantially overlapping classification state. Actual dispatch often uses
booleans even after a target family has been resolved.

**Impact.** Adding a supported composition requires updating several parallel
representations and makes it difficult to identify the authoritative decision.

**Direction.** Resolve supported numerical choices once into a closed
compilation/lowering decision and consume that decision directly. Do not
replace enum duplication with string dispatch or a registry of complete target
classes.

### PSR-014 — The scalar plan is compositional but complete registrations re-collapse it

**Observed.** Scalar handlers lower semantic terms into component contributions,
but registration logic also recognizes exact complete sets of terms and maps
them to named supported variants.

**Impact.** Some exact-combination validation is legitimate, but component
semantics are repeatedly reclassified as whole targets.

**Direction.** Keep explicit support/lowerability validation while avoiding a
second whole-problem ontology. A closed plan may reject unsupported
combinations without naming every supported Cartesian product as a model type.

### PSR-015 — Compiler service packaging is repeated

**Observed.** Target branches repeatedly perform the same mechanical work:
construct a complete model, obtain metric/constraint, wrap state/adjoint member
functions in lambdas, assign the executable, and expose optional Hessian data.

**Impact.** New targets copy composition code and policy wiring.

**Direction.** Centralize mechanical service packaging only after the typed
result shape is clear. Do not create a common complete-model base class merely
to share wrapping code.

### PSR-016 — `CompiledProblemT` is good compiler-path packaging, not a universal integration type

**Observed.** It validates compatible layouts and packages executable, metric,
optional constraint/Hessian, solve callbacks, manifest, lifetime owner, and
box data. It can construct `ReducedDTOT` cleanly.

**Good.** The package is useful for the native compiler path.

**Direction.** Preserve or simplify it as compiler-path packaging. External
applications should converge on solver ports without adopting compiler
provenance or `ProblemSpec`.

### PSR-017 — Manifest/provenance data has duplicate owners

**Observed.** `ResolvedCompilationDecision`, `CompilationManifest`, legacy flat
fields, compatibility views, strings, and projections repeat realized facts.

**Impact.** Compiler cleanup is constrained by copied representations and
artifact compatibility can accidentally become an excuse to keep duplicate
in-memory state.

**Direction.** Keep one structured in-memory owner per realized fact and render
legacy/persisted forms once at the serialization edge. Preserve persisted
fields only when a current consumer requires them.

### PSR-018 — Native output is misplaced in numerical model types

**Observed.** Several complete targets expose filesystem/native writer methods.
B1/B2 recover concrete models to call those methods.

**Impact.** Numerical executable objects become application/filesystem objects,
and type erasure cannot remain one-way.

**Direction.** Retain typed discretization/native data in the application path
and move output ownership outward. Share concrete writer code only where two
current applications prove a stable common operation.

### PSR-019 — B1/B2 execution adapters are legitimate composition roots but overloaded

**Observed.** They bind runtime data, compile, create initial points, perform
derivative/Hessian evidence, choose/run optimizers, extract native diagnostics,
write fields, and build experiment artifacts.

**Good.** Application-level composition is the correct place for backend data,
solver choice, experiment evidence, and output.

**Pressure.** Numerical compilation, verification evidence, optimizer execution,
and output are interleaved in large target-specific adapters.

**Direction.** First remove the typed-model recovery cycle. Further sharing or
splitting should occur only when it deletes demonstrated duplication; do not
create a benchmark framework merely to shorten the files.

### PSR-020 — The direct v0 scalar model is a migration oracle and probable deletion target

**Observed.** `ScalarDiffusionReactionModel` combines assembly, executable
operations, state/adjoint solves, Hessian, OTD/KKT helpers, metrics,
constraints, and output. The newer compiler path duplicates much of its
functionality.

**Good.** It contains valuable hand-verified numerical oracles and reference
coverage.

**Direction.** Transfer unique tests/behavior to the surviving path before
removing the redundant complete implementation. Do not decompose v0 into a new
public subsystem solely to preserve it.

### PSR-021 — Scalar-specific KKT deserves a consumer audit

**Observed.** A backend-neutral quadratic KKT product and serial adapter already
exist, while a scalar-specific deal.II KKT realization remains tied to the
legacy scalar model.

**Direction.** Trace non-test consumers. If it contributes no independent
behavior after v0 migration, delete it and retain the generic product/adapter.
If independent behavior remains, document that behavior before deciding.

### PSR-022 — Semantic representation is sound; syntax/construction can be simplified later

**Observed.** `RequirementPolicySpec` encodes mutually exclusive typed policy
payloads through many `std::optional` members. Reference builders often clone a
complete graph and mutate components by stable string ID.

**Impact.** Invalid combinations remain representable in C++, and reference
construction can hide semantic deltas behind vector mutation.

**Direction.** This is not the blocking numerical-boundary problem. Later
cleanup may use actual sum types, selectively tagged IDs, named semantic deltas,
and physical file partitioning. Avoid a symbolic PDE DSL or fluent builder that
merely hides the same mutation.

### PSR-023 — Large files are symptoms, not refactor goals

**Observed.** `dealii_compiler.hpp`, semantic validation/reference builders,
runner code, and complete numerical models are large.

**Impact.** Large context burdens humans and agents, but splitting files alone
does not remove duplicated decisions or ownership.

**Direction.** Reorganize files when a responsibility has moved or disappeared.
Use smaller files as a consequence of decomposition, not as the primary
success metric.

### PSR-024 — Documentation authority must be updated deliberately

**Observed.** Existing design documents already describe a component graph and
small executable ports, while implementation/reference/planning documents still
contain current v0/v1 and application-specific facts that will become stale as
cleanup proceeds.

**Direction.** Establish the new boundary and audit now. During each
implementation unit update only the authoritative document whose public or
implemented behavior changed. Perform a final repository-wide stale-reference
sweep after the architecture stabilizes.

## Cross-layer numerical conclusions

### State side

The state side decomposes into at least these independent responsibilities:

```text
StateSpace
StateCoordinates / reconstruction
PhysicalStateOperator
Gauge / nullspace policy
State/adjoint solve realization
```

`ScalarComponentModel` already demonstrates the useful coordinate actions
$P$, $P^T$, and fixed lifting. `DirichletControlLiftingModel` adds $L_D$ and
its pullback. Full-coordinate identity-row schemes are alternative coordinate
policies, not separate PDE families.

The physical operator may be fixed symmetric, fixed nonsymmetric, or
point-dependent. Therefore a universal typed numerical operator must not expose
only a persistent sparse matrix.

### Decision side

Current decision mechanisms are:

```text
affine forcing:
    volume control
    Neumann control

state reconstruction:
    Dirichlet control

operator parameter:
    coefficient identification
```

The common compiler concept is a selected decision realization. The common
runtime solver concept remains the full residual/JVP/VJP contract.

### Observation/objective side

The same physical state can feed volume, subdomain, $H^{1}$, trace, sensor, or flux
observations. Efficient quadratic tracking often lowers to $Q, q, c$; explicit
measurement actions are retained only where another consumer needs them.

Objective regularization belongs to the objective even when its numerical
operator is also used by the optimizer metric.

### Optimization geometry and second order

Metrics are algorithmic primal-dual maps and may be changed independently of
the objective. Constraints are optional capabilities supported by selected
decision realizations. Reduced Hessians are optional composed formulation
capabilities and must not be inferred from first-order ports.

## Current-to-target responsibility map

| Current concept | Assessment | Target direction |
| --- | --- | --- |
| `ExecutableModelT` | keep | stable solver-facing facade |
| `StateAdjointSolversT` | keep | state/adjoint IoC boundary |
| `ReducedDTOT` | keep | first-order reduced orchestration |
| `MetricT`, `ConstraintT`, `ReducedHessianT` | keep | generic optional solver/formulation capabilities |
| `CompiledProblemT` | keep/simplify | compiler-path package, not universal app API |
| `ScalarLoweringPlan` | keep/evolve | evidence for component-aware closed plan |
| `VolumeObservationAssembly` | keep/reuse | focused typed tracking component |
| `NeumannControlRealisation` | narrow/evolve | preserve demonstrated runtime realization polymorphism; remove unrelated roles if migration pays |
| complete compiler `*Model` classes | migrate/narrow/delete | composition containers, not target architecture |
| model-owned native output | move/delete | application/native view |
| B1/B2 concrete executable downcasts | delete | retain typed realization explicitly |
| duplicate target enums/booleans/predicates | delete/simplify | one authoritative closed decision |
| duplicate manifest/compatibility state | delete/simplify | one structured owner + edge projection |
| direct v0 scalar implementation | probable delete | preserve unique oracle/behavior first |
| scalar-specific KKT | audit then likely delete | generic KKT product + serial adapter if sufficient |
| semantic graph/validator | keep | optional later syntax cleanup |
| recipes/scenarios/envelope | keep | current boundaries are useful |

## Naming observations

Naming changes should follow responsibility changes rather than precede them.
The following names deserve review during their migration units:

- `*Model` for types whose actual role becomes a discretization, reconstruction,
  operator bundle, or compiler composition result;
- `VolumeObservationAssembly` if the surviving responsibility is specifically a
  quadratic tracking-form assembly rather than an observation map;
- `TraceHhalfMetric` if a neutral trace quotient/Riesz operator is eventually
  shared directly by both objective regularization and metric roles;
- internal compiler names that encode whole target families after the closed
  plan has become component-oriented;
- `CompiledProblemT` only if its surviving compiler-path responsibility becomes
  materially narrower or broader.

Do not schedule standalone naming commits unless they remove misleading
ownership or accompany a real responsibility migration.

## Decisions deliberately deferred

The audit fixes architectural direction but does not fix every C++ shape.
These decisions belong to the unit with enough current-tree evidence:

- whether state-coordinate data should be a concrete value, concrete class, or
  a few free functions around shared matrices;
- whether volume and Neumann affine coupling share a reusable concrete helper;
- whether regularization needs any named internal type beyond callbacks/values;
- whether a neutral H1/2 operator should be extracted from the existing metric;
- whether a generic quadratic-form value type deletes enough repeated code;
- the exact shape of the compiler's closed lowering-plan variant;
- how much fixed-linear reduced-Hessian code is genuinely reusable;
- whether full-coordinate identity-row state representation survives after
  reconstruction reuse is expanded;
- how application-native typed output state is represented after downcasts are
  removed; and
- final replacement names for complete `*Model` types that disappear or narrow.

Each roadmap unit must state decision criteria rather than guessing these
answers in advance.

## Documentation implications

The new documentation set has five distinct roles:

- `docs/design/pde-solver-boundary.md`: long-lived accepted architecture and
  ownership invariants;
- `architecture-map.md`: detailed current-system topology, runtime/compiler
  flow, typed-numerics axes, ownership, and current-to-target correspondence;
- this assessment: evidence, findings, judgments, and deferred decisions;
- `deletion-ledger.md`: operational migration/deletion state;
- `roadmap.md`: bounded executable work units and mutable handoff.

Existing `system-blueprint.md` and `composition-boundaries.md` are broadly
compatible with the target and should not be rewritten merely to repeat this
audit. Implementation/reference docs should change when their corresponding
code contract changes. A final documentation sweep should search for deleted
type names, stale v0/v1 path descriptions, old target identities, concrete
output ownership, and obsolete roadmap handoffs.

## Acceptance criteria for the audit

The architecture-audit phase is complete when:

1. the long-lived boundary decision is accepted;
2. every major complete numerical target has an entry in the deletion ledger;
3. current dependency, runtime, ownership, typed-numerics, and
   current-to-target maps are recorded in the architecture atlas;
4. the roadmap references findings instead of duplicating assessment prose;
5. unresolved implementation choices have explicit decision criteria;
6. no new public hierarchy is justified only by aesthetics or file size;
7. the external-application proof and native compiler cleanup are both covered;
8. documentation authority and later stale-doc cleanup are explicitly routed;
   and
9. no production code has changed as part of the audit itself.
