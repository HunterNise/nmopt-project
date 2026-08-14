# C1/C2 preparation review and remediation handoff

## Status and authority

This document records the review of the inclusive commit range `f467d49`
through `349d5f1` on 2026-08-12. All commits in the range were authored on
2026-08-09. The range prepared and then declared complete the conditional C1
compiler-product gate and C2 bounded component-lowering gate before P5.1.

This is a static review and implementation handoff, not a second status
ledger. The [implementation roadmap](../../implementation-roadmap.md) remains the
sole owner of mutable status and acceptance state. The
[Stage B roadmap](../pre-ch5-ch6/stage-b-roadmap.md),
[refactor assessment](../pre-ch5-ch6/assessment.md),
[interface specification](../../../design/interface-specification.md),
[implementation-readiness review](../../../implementation/implementation-readiness-review.md),
and [v1 semantic compiler contract](../../../implementation/v1/semantic-compiler.md)
remain authoritative for the C1/C2 contracts.

The preparation established several useful foundations, and the initial
review recorded five issues. The bounded remediation is now
acceptance-complete for both gates; the findings and their historical review
context remain below for traceability:

1. `C1-R1` – the public compilation-result boundary still throws for a null
   session and does not validate ordinary `Function` field shape;
2. `C1-R2` – the structured manifest is reconstructed from target enums,
   roles, and prose and does not identify all executable inputs;
3. `C1-R3` – compiled observation-space dimensions are inferred from input
   coordinate layouts rather than the realized maps;
4. `C2-R1` – the resolver is invoked, but most compiler decisions still rescan
   the raw graph and rebuild whole-target flags; and
5. `C2-R2` – scalar handlers contribute metadata while central dispatch and
   model constructors still select complete recipes.

The final range commit, `349d5f1`, marked C1 and C2 complete while naming the
first P5.1 registered target as the next task. C2's own exit required a
selected new scalar recombination before completion. Existing history must not
be rewritten without explicit authorization; the roadmap must instead show
the current remediation state and the repaired commits must satisfy the
original gates.

## Commit classification

### C1 compiler inputs, products, lifetime, solves, and provenance

- `f467d49` – shared serial SPD solve reporting and migration of the initial
  target families;
- `b9235a0` – stable-ID semantic resolution;
- `2f2e323` – reduced solve reports and owned reduced-service lifetime;
- `5a20f06` – migration of the assembled scalar target to the shared solve;
- `c42a0b1` – `Function` provenance labels;
- `b0ca1ca` – independently selected state and adjoint solve policies;
- `aff82b7` – runtime compilation-input diagnostics;
- `f357f97` – owned deal.II compilation sessions;
- `e681f4c` – honest naming of the capability-only registry;
- `effcb3b` – typed manifest records; and
- `87c9796` – compiler diagnostics, lifetime, solve-report, and manifest tests.

### C2 bounded component lowering

- `1894521` – resolved scalar contribution plans and stored handler records;
  and
- `ed3929f` – integration of the plan into the assembled scalar target and
  manifest.

### Status documentation

- `349d5f1` – v1 compiler and roadmap claims that C1 and C2 are complete.

No commit in the range is incidental to the requested review. The status
commit is evidence for the claimed exit state rather than implementation
evidence.

## Retained behavior

The remediation must preserve the following behavior that static inspection
found consistent with the selected contracts:

- `SemanticResolver` validates before constructing a borrowing
  `ResolvedProblemView` and indexes every semantic component by stable ID;
- the resolver and scalar-plan tests establish declaration-order independence
  for their covered graphs;
- `SPDLinearSolvePolicy` separates state and adjoint tolerances and the shared
  serial service reports convergence, iterations, requested tolerance, and
  achieved residual;
- the pure-Neumann augmented direct solve remains a separate typed policy
  rather than being forced into the SPD service;
- `ReducedDTOT` can own its executable and a lifetime token, and detached
  compiled DTO services retain the model they call;
- `DealiiCompilationSession` owns a moved static triangulation and is retained
  by both the compiled problem and detached DTO;
- scalar diffusion, reaction, regularisation, bound ordering/layout, solve
  policies, empty meshes, and requested boundary IDs have focused
  lowerability diagnostics;
- the metric/constraint projection relationship uses the non-spoofable
  realization witness introduced before C1 rather than a display ID;
- `DealiiCapabilityRegistryV1` is honestly documented as a capability ledger,
  not a lowerer registry;
- the first scalar plan records stable component and handler IDs and rejects
  specialized Neumann graphs outside its bounded path;
- the existing hand-integrated Q1/DGQ0 oracle remains independent of the
  compiled/direct packaging comparison; and
- `349d5f1` records a clean `debug-dealii` build observation of 81.3 seconds
  elapsed and 1.53 GiB peak RSS, satisfying the required C2 cost
  remeasurement without justifying speculative template splitting.

No build or test command was run during this static review. The findings come
from the historical diffs, the code at `349d5f1`, the current implementation,
the committed tests, and the authoritative contracts.

## Required reading for the implementing agent

Before changing code, follow the repository routing in
[`conventions/README.md`](../../../conventions/README.md). Read only the bounded
authorities needed for the selected work unit:

- [Stage B C1 gate](../pre-ch5-ch6/stage-b-roadmap.md#c1--make-compiler-inputs-products-and-provenance-explicit);
- [Stage B C2 gate](../pre-ch5-ch6/stage-b-roadmap.md#c2--introduce-bounded-component-lowering);
- [assessment C1 evidence and exit](../pre-ch5-ch6/assessment.md#c1--make-compiler-inputs-products-and-provenance-explicit);
- [assessment C2 evidence and exit](../pre-ch5-ch6/assessment.md#c2--replace-whole-target-dispatch-with-bounded-component-lowering);
- [interface compiler obligations](../../../design/interface-specification.md#62-compiler-obligations);
- [composition boundaries](../../../design/composition-boundaries.md#what-changes-when-one-component-changes);
- [selected lowering policy](../../../implementation/implementation-readiness-review.md#5-lowering-and-compiler-capabilities);
- [required verification and provenance](../../../implementation/implementation-readiness-review.md#12-required-verification-and-provenance);
- [v1 validation and diagnostics](../../../implementation/v1/semantic-compiler.md#validation-and-diagnostics);
- [v1 registered deal.II realization](../../../implementation/v1/semantic-compiler.md#registered-dealii-realization);
- [P5.1 data-placement repair](p5.1-remediation-review.md#p51-r1--coefficient-and-robin-data-have-missing-or-wrong-regions);
- [P5.1 typed-policy repair](p5.1-remediation-review.md#p51-r2--boundary-and-conormal-policies-are-prose-only-selections);
- [P5.2 realized-dimension repair](p5.2-remediation-review.md#p52-r3--weighted-trace-dimensions-fall-through-to-state-dimensions);
- [P5.3 realized-dimension repair](p5.3-remediation-review.md#p53-r3--incorrect-compiled-normal-flux-observation-dimension); and
- [P5.4 registration repair](p5.4-remediation-review.md#p54-r1--unregistered-section-511-cross-products-can-compile) and
  [realized-dimension repair](p5.4-remediation-review.md#p54-r3--transformed-state-observation-dimensions-are-misreported).

Do not turn this remediation into a runtime plugin framework, universal weak-
form engine, public symbolic DSL, PDE-family class hierarchy, distributed
backend, or all-target rewrite. The selected endpoint remains one modest
scalar contribution path plus explicit specialized strategies whose
coordinate, trace, reconstruction, nullspace, or parameter dependence is
genuinely different.

## C1-R1 – the compilation diagnostic boundary is incomplete

**Severity:** P1 – the documented public compiler contract can still accept or
throw on predictable caller inputs instead of returning exact lowerability
diagnostics.

### Evidence

`aff82b7` added useful early checks for scalar coefficients, regularisation,
bound representation/layout/order, solve policies, mesh emptiness, and
boundary presence. The v1 compiler document consequently states that
caller-provided compilation data use a predictable diagnostic boundary and
that `ContractError` is reserved for direct low-level constructor misuse and
internal invariants after validated lowering.

Two public inputs remain outside that boundary:

1. The `compile()` overload taking
   `std::shared_ptr<DealiiCompilationSession<dim>>` calls
   `contract::require(session)` before constructing a `CompilationResultT`.
   A null public compilation input therefore throws rather than returning a
   `lowerability` diagnostic.
2. `DealiiDataBindings` stores forcing, desired-state, and optional fixed-
   Dirichlet values as base `dealii::Function<dim>` references. The compiler
   does not compare `n_components` with the scalar semantic `DataSpec`. A
   multi-component `Function` is accepted and the backend evaluates component
   zero through `value(point)`. Later weighted-trace and general-scalar paths
   added local scalar-shape checks, which makes the missing common validation
   explicit rather than intentional.

This is not a request to sample analytic data and prove mathematical
properties. Component count and a null session handle are structural runtime
facts available before model construction.

### Required outcome

Every predictable public compilation input must use the documented result
channel. Build one resolved binding-validation pass that matches each supplied
binding to its semantic datum and selected realization. At minimum it must
check:

- required/forbidden binding presence;
- scalar/vector/tensor `Function` component shape;
- exact role, semantic data ID, space, and region expected by the resolved
  lowerer;
- scalar coefficient domains already checked by `aff82b7`;
- bound representation, layout, finiteness, and order;
- solve-policy validity; and
- session presence and mesh-selection validity.

A null session passed to the public compiler should produce one exact
`lowerability` diagnostic. If the project instead chooses to make the session
overload a low-level preconditioned API, add a separately named checked public
entry point and document the distinction. Do not leave two overloads named
`compile()` with different invalid-input channels.

The validator must consume the one resolved compilation request required by
`C2-R1`; it must not reproduce another role-based switch over raw graph
vectors.

### Implementation instructions

1. Introduce a backend-neutral resolved binding request containing semantic
   data ID, required field shape, space/region, runtime representation, and
   selected evaluation rule.
2. Populate it while resolving the lowerer/registration, before deal.II model
   construction.
3. Validate ordinary forcing, desired-state, and fixed-data `Function` component
   counts through the same mechanism used for weighted and general-scalar
   `Function` objects.
4. Translate a null public session to one exact diagnostic such as
   `compilation_session_presence`.
5. Preserve `ContractError` for direct model/metric/constraint constructor
   misuse and impossible internal states after successful validation.
6. Remove duplicated target-local shape checks only after the common pass
   covers the same exact diagnostics and all callers.
7. Update the v1 diagnostic table with the final checked/throwing boundary.

### Required tests

Add table-driven compiler cases for:

- a null owned session;
- two-component forcing;
- two-component desired state;
- two-component fixed-Dirichlet data;
- the existing weighted-trace scalar weight and general-scalar coefficient
  shape failures routed through the common pass;
- a binding attached to the wrong resolved semantic data ID or region; and
- direct low-level constructor misuse that must still throw.

Each compiler rejection must assert diagnostic category, component ID, and
capability. Prove that no backend model constructor or quadrature evaluation
runs after a common binding diagnostic.

## C1-R2 – manifest records are reconstructed and do not identify the run

**Severity:** P1 – blocks the C1 exit that exact structured manifests identify
every current target and distinguish materially different compiled products.

### Evidence

`effcb3b` introduced useful typed subrecords, but
`DealiiCompiler::make_manifest()` still reconstructs provenance after the
executable is chosen. At `349d5f1`, and more extensively in the current code,
it receives a raw `ProblemSpec`, a central `CompiledTargetKind`, a constraint
enum, selected regions, caller bindings, and completed services. It then
rebuilds target booleans and switches on semantic roles.

Several records are not identifying:

- `CompiledBindingRecord` contains only semantic ID, role, representation,
  and a provenance string. It omits the datum's semantic space/region and the
  resolved evaluation rule.
- diffusion, reaction, and regularisation values are converted with
  `std::to_string`, which is a display conversion rather than a lossless typed
  value record.
- lower/upper bound records say only `caller-supplied compiled bound data` and
  a scalar/vector representation. Different values with the same
  representation produce the same structured records.
- the borrowed mesh overload always records `caller-owned triangulation` plus
  dimension and active-cell count. Meshes with different geometry, material
  IDs, or boundary partitions can therefore share the same mesh record.
- human-readable fields such as `data_rule`, `observation_realisation`,
  `lifting_realisation`, and `state_adjoint_solve_policy` are independently
  reconstructed from target flags. They are not renderings of closed typed
  records and can contradict the actual plan.
- `lowering_handler_records` are strings appended from the scalar plan or
  target-specific branches; they do not establish that those decisions
  constructed the returned services.

The test helper added in `87c9796` primarily checks that records are nonempty,
that some roles occur, and that selected display strings contain expected
substrings. It does not compare complete structured records for every target
or prove that changed inputs produce changed provenance.

The later P5.1 through P5.4 reviews exposed concrete consequences of this root
design: wrong or absent binding regions, prose-only policies, unclosed target
cross-products, and observation dimensions reconstructed from unrelated
layouts.

### Required outcome

Create one immutable resolved compilation-decision record before model
construction. The record must be the source for backend construction and for
the manifest. The manifest is a descriptive projection of that record plus
facts reported by the realized services; it must not rediscover choices from
`ProblemSpec`, `ProblemSpec::id`, `CompiledTargetKind`, or display strings.

The closed record needs typed entries for:

- structural registration ID and formulation;
- semantic and realized spaces/pairings;
- mesh/discretization identity and lifetime policy;
- every binding's semantic ID, field shape, space, region, runtime
  representation, evaluation rule, provenance, and reproducible value or
  digest where applicable;
- every transformation, residual, observation, loss, metric, constraint, and
  solve realization actually selected;
- exact state, adjoint, metric-inverse, and specialized solve policies;
- nullspace, trace, reconstruction, boundary, and other typed policies; and
- user-assumed versus compiler-checked assumptions.

Do not serialize opaque `Function` contents by sampling them. Continue
requiring caller provenance for analytic `Function` objects. For scalar and
vector bound data, record typed values or a deterministic exact-layout digest.
For floating-point display, use a lossless representation or retain the value
as `double` in the typed schema; do not use `std::to_string` as identity.

For meshes, require caller provenance on the borrowed path and add a
deterministic structural summary or fingerprint covering the selected
geometry/topology policy, active cells, material IDs, boundary IDs, and
relevant refinement state. A user label alone may remain provenance, but it
must not be presented as compiler-verified mesh identity.

### Implementation instructions

1. Define a `ResolvedCompilationDecision` or equivalently named internal
   record owned by the compiler layer. Document the architectural decision in
   the v1 compiler contract.
2. Populate it from the resolved graph, typed realization policies, concrete
   binding validation, and selected service factories.
3. Remove target booleans from manifest construction. A target/registration ID
   may be recorded, but it must not be used to reconstruct component choices.
4. Extend `CompiledBindingRecord` with space, region, field shape, evaluation
   realization, typed/digested value provenance, and checked/assumed status.
5. Extend the mesh record with caller provenance and a deterministic
   structural identity summary. Require explicit provenance for the borrowed
   overload.
6. Add typed transformation, residual/observation, pairing, boundary-policy,
   and assumption records as the P5 repairs supply their missing contracts.
7. Generate all compatibility prose from those typed records in one renderer.
8. Bump `schema_version` and document compatibility because the structured
   record meaning changes materially.
9. Keep the manifest descriptive. Do not feed rendered manifest data back into
   execution.

### Required tests

For every registered target retained during the work unit, compare exact
structured records, not substrings. Add differential cases proving that the
appropriate record changes for:

- two scalar coefficient values that collide under six-decimal formatting;
- scalar versus exact-layout vector bounds and two different bound values;
- two meshes with equal active-cell counts but different boundary/material
  partitions;
- borrowed versus owned mesh lifetime and provenance;
- independent state and adjoint solve policies;
- a changed observation or transformation with an unchanged residual; and
- arbitrary label/display-text edits with unchanged typed execution.

Add a consistency check that renders the human-readable compatibility fields
from the typed records and cannot construct contradictory prose independently.

## C1-R3 – realized output spaces are inferred from input coordinates

**Severity:** P1 – the C1 manifest already misdescribed a registered boundary
observation at the point the gate was marked complete.

### Evidence

`compiled_space_dimension()` special-cases semantic state, decision, and test
spaces. For an observation output it returns the decision dimension when the
input is the decision and otherwise returns the state-variable dimension.

At `349d5f1`, the registered Neumann target included a state boundary-trace
observation realized by face evaluation and a face pairing. Its output was
therefore recorded with the volume state-coordinate dimension, even though no
realized trace map reported that dimension. Fixed reconstruction likewise
distinguishes independent state coordinates from the physical field consumed
by its observation, but the manifest had no realized map record with which to
state the distinction.

Later point-sensor work added a target-specific dimension special case, while
weighted trace, normal flux, and transformed observations retained or exposed
the same fallback. The P5.2, P5.3, and P5.4 reviews recorded those later
instances. The root defect predates them: semantic source-variable identity is
not a compiled output-space realization.

`CompiledSpaceRecord` also says only `lowered observation coefficients`; it
does not record output ordering, pairing weights, trace/projection rule, or the
value/JVP/VJP object that owns the dimension. Fused observation-loss assembly
therefore leaves manifest generation to guess.

### Required outcome

Every compiled map must publish its realized source and target spaces from the
same object that supplies its executable actions. Introduce a common realized-
map record containing at least:

- semantic component and source/target space IDs;
- realization ID;
- source and output dimensions;
- layout and stable ordering;
- pairing/quadrature realization;
- transformation chain, if any; and
- value, JVP, and transpose/VJP provenance.

An observation fused with a quadratic loss may retain fused assembly for
efficiency, but it must expose an equivalent bounded map and pairing whose
composition is tested against the fused objective and covector. Do not infer
an observation dimension from its input block.

Coordinate one implementation with `P5.2-R3`, `P5.3-R3`, and `P5.4-R3`.
Cover the baseline unweighted Neumann boundary trace as well as the later
weighted trace, normal flux, point sensors, ordinary volume observations, and
transformed physical-state observations.

### Implementation instructions

1. Add backend-neutral `CompiledRealizedSpaceRecord` and
   `CompiledRealizedMapRecord` types, or equivalent closed types.
2. Make each bounded observation/transformation lowerer construct the record
   while constructing its executable map.
3. Choose and document an actual output representation for each fused map:
   FE coefficients with an explicit mass pairing, ordered quadrature samples,
   sensor values, or another proved-equivalent realization.
4. Obtain manifest dimensions only from the realized record. Delete the
   observation input-layout fallback.
5. Generate observation prose and pairing metadata from the same record.
6. Preserve semantic infinite-dimensional spaces separately from their
   compiled finite-dimensional realizations.

### Required tests

Use meshes for which state coordinates, physical state DoFs, boundary trace
samples, facewise controls, and sensor counts are distinct. Check:

- value, JVP, and VJP/pairing for every realized map;
- equality between explicit map-plus-loss composition and fused assembly;
- dimension and ordering reported by the realized record;
- transformed physical-state versus independent-coordinate dimensions;
- baseline and weighted boundary traces;
- normal flux and point sensors; and
- structural manifest equality for unaffected residual, metric, constraint,
  and solve records.

## C2-R1 – the resolved view is not the compiler source of truth

**Severity:** P1 – blocks the C2 requirement for one ID-resolved closed view
and allows validation, dispatch, planning, and provenance to drift.

### Evidence

`b9235a0` added `ResolvedProblemView`, and `ed3929f` changed `validate()` and
`compile_impl()` to invoke `SemanticResolver`. The resolved view is then used
only by `DealiiScalarLoweringPlanner` for the already-selected assembled path.

Before planning, `compile_impl()` rescans the raw `ProblemSpec` to derive
booleans such as fixed reconstruction, controlled Dirichlet, Neumann control,
mean-zero gauge, $H^{1}$ control, coefficient identification, selected tracking
region, and complete-target family. Lowerability and formulation validation
perform further raw-vector searches. Model construction and manifest
generation repeat more selectors. The current compiler has extended the same
pattern into a long `CompiledTargetKind` selection across P5.1 through P5.4.

The scalar planner itself still iterates raw `specification.observations`,
`losses`, and `requirement_policies` for several relations rather than asking
the resolved view for closed selected ports and policies. `ResolvedProblemView`
is therefore an ID index, not yet the one resolved compiler request required
by C2.

This split enabled several later defects: a semantic policy can be accepted
while a raw helper or backend uses another region; manifest decisions can be
reconstructed from a neighboring flag; and unregistered cross-products can
fall through a target-kind ternary.

### Required outcome

After semantic validation, construct exactly one closed, immutable compiler
request. It must resolve:

- formulation state, decision, equation, metric, and optional constraint;
- every equation term and its exact typed variable/data/region ports;
- every observation/loss edge;
- transformations and physical-coordinate chains;
- selected typed requirement policies with status and scope;
- binding requests;
- candidate bounded registration and its explicit exclusions; and
- backend/discretization capabilities required by every component.

All later lowerability checks, mesh checks, scalar planning, specialized
strategy selection, service construction, and manifest records must consume
that request. Do not retain parallel raw-graph feature predicates for the same
decision.

This does not require copying every semantic object. Stable references or IDs
owned by a compilation-scoped record are sufficient, provided lifetime and
closure are explicit.

### Implementation instructions

1. Extend resolution with selected-edge accessors or introduce a compiler-
   layer `ResolvedCompilationRequest` built from `ResolvedProblemView`.
2. Move formulation, tracking-region, boundary-policy, transformation,
   metric/constraint, and target-registration resolution into that one pass.
3. Make diagnostics name the exact resolved component and missing capability.
4. Change `validate_lowerability()`, mesh validation, scalar planning,
   specialized strategy matching, and manifest construction to accept the
   resolved request instead of raw `ProblemSpec`.
5. Delete `uses_*`, `has_*`, and `selected_*` helpers only as their resolved
   equivalents become the sole source. Do not perform a blind mechanical
   rewrite without focused tests.
6. Keep `ProblemSpec::id` and labels descriptive; neither may select a
   registration.
7. Update the v1 compiler flow diagram and exact exclusion language.

### Required tests

Add compiler-level order-independence tests, not only resolver tests:

- reorder every semantic declaration vector and require identical resolved
  decisions, executable behavior, and structured manifest;
- rename problem/component labels while retaining IDs and typed ports;
- change a selected policy region and prove that mesh validation, model
  construction, and manifest all see the same changed region or reject it;
- mutate one observation without changing the residual request;
- mutate one metric without changing objective/residual requests; and
- present a structurally unregistered cross-product and require one exact
  registration diagnostic before backend construction.

## C2-R2 – handlers do not yet lower independent executable contributions

**Severity:** P1 – the C2 plan records handler names but the complete target
recipe remains selected and assembled elsewhere.

### Evidence

`1894521` introduced stored residual, observation, loss, metric, constraint,
and transformation handlers. Their `contribute()` functions append enum/data
records and strings to `ScalarLoweringPlan`. `ed3929f` integrated the plan only
after raw whole-graph booleans selected `uses_assembled_v1_target`.

At the gate boundary, `ScalarComponentModel` required the presence of exactly
diffusion-reaction, source, volume-control, tracking, and control-
regularisation contributions. Its inherited assembly loop always assembled
that complete recipe. Metric and constraint services were still constructed
by the central compiler branch. A handler record therefore proved that a kind
was recognized, not that the handler constructed or parameterized the
executable action.

The current implementation expands the same pattern:

- one constructor requires the complete diffusion-reaction recipe;
- another requires the complete tensor/transport/reaction/Robin recipe;
- `assemble_physical_operators()` branches on whether a general-scalar binding
  pointer is non-null rather than iterating plan contributions;
- observation behavior is selected by model booleans;
- metric/constraint services remain centrally constructed; and
- a long `CompiledTargetKind` ternary still decides the target family and
  later reconstructs provenance.

The tests in the reviewed range checked plan counts, handler strings, fixed
transformation selection, rejection of a Neumann graph, and the two existing
assembled targets. They did not compile a new combination of independently
selected residual and observation components. `349d5f1` nevertheless marked
C2 complete and instructed the next agent to implement the first P5.1 target,
although the C2 exit required that selected recombination as its proof.

Later P5.1 work did add a general scalar/Robin target without a new problem
class and supplied useful term-isolation tests. It still added another exact
constructor recipe, central target selection, and hard-coded binding placement.
It demonstrates the value of the shared model but does not close the original
C2 component-lowering contract.

### Required outcome

Make the bounded scalar plan causally responsible for the executable
contributions it records. The smallest acceptable design is:

- one shared scalar FE/discretization context;
- typed residual contribution descriptors consumed by cell/face assembly;
- typed data-placement and evaluation requests;
- typed observation and loss descriptors consumed by objective assembly;
- transformation strategy/value/JVP/VJP selection from the plan;
- metric and constraint service factories selected from the plan and coupled
  by their realization witness; and
- provenance emitted from those same resolved descriptors and factories.

Handlers may contribute stable enums plus typed payloads to efficient fused
assembly; they do not need heap-polymorphic kernels or one matrix per term.
The assembly loop may fuse contributions, but it must iterate or otherwise
consume the selected contribution records. A nonselected contribution must
not be assembled, and every selected contribution must carry its exact ports,
data placement, region, policy, and transpose behavior.

Keep genuinely distinct strategies specialized: controlled-Dirichlet
reconstruction, independent boundary-control layouts, continuous controls,
coefficient-dependent state matrices, mean-zero saddle systems, and declared
transposition policies. Select those strategies through closed structural
registrations, not problem names or an order-dependent target ternary.

### Implementation instructions

1. Split scalar plan payloads from display provenance. Give each contribution
   its exact resolved variables, data bindings, region, pairing, realization,
   and handler/factory ID.
2. Refactor volume and face assembly so selected residual contributions drive
   local accumulation. Preserve the existing weak signs and exact transpose.
3. Refactor observation/loss assembly so changing an observation leaves the
   residual plan and matrix unchanged.
4. Construct the scalar metric and optional box from the resolved metric and
   constraint contribution records rather than `CompiledTargetKind`.
5. Replace base/general-scalar constructor recipe assertions with a closed
   bounded registration matcher that validates supported combinations before
   model construction.
6. Retain specialized strategies behind separately typed registration
   records. Do not migrate all targets in one flag day.
7. Produce handler/realization manifest records from the objects that actually
   assembled or constructed the services.
8. Coordinate boundary/data payloads with `P5.1-R1` and `P5.1-R2`, and use the
   same closed registration mechanism required by `P5.4-R1`.

### Required tests

Add one explicit recombination matrix whose rows vary only one component
family at a time. At minimum include:

- the retained diffusion-reaction residual with full-volume and material-
  subdomain observations;
- fixed reconstruction combined with the material-subdomain observation;
- the general scalar/Robin residual with an otherwise unchanged compatible
  observation/loss/metric/constraint set;
- optional box present/absent with identical residual/objective actions; and
- one unsupported cross-product rejected before model construction.

For every supported row, verify:

- independent residual value, JVP, and VJP contribution oracles;
- observation value/JVP/VJP and objective derivative;
- reduced Taylor remainder;
- residual invariance when only observation or metric changes;
- objective/adjoint-source invariance when only the search metric changes;
- metric/constraint realization coupling; and
- structured plan/manifest equality for every unchanged component record.

Add a fault-injection or test registry variant showing that removing one
handler prevents that contribution from lowering rather than silently using a
central fallback.

## Ordered implementation plan

### Work unit 1 – one resolved compilation request and binding boundary

**Outcome:** semantic resolution produces one closed compiler request, and all
predictable public binding/session errors use exact diagnostics.

**Boundary:** `C2-R1` and `C1-R1`. Do not change numerical assembly, target
registrations, or manifest schema in this unit.

**Likely files:**

- `docs/implementation/v1/semantic-compiler.md`;
- `include/nmopt/semantic/v1/resolved_problem.hpp`;
- `include/nmopt/compiler/v1/dealii_types.hpp`;
- `include/nmopt/compiler/v1/dealii_compiler.hpp`;
- `tests/semantic_v1_contract.cc`; and
- `tests/dealii_diffusion_contract.cc`.

**Focused verification:** semantic resolution, compiler diagnostics, owned
session, weighted-trace binding, and general-scalar binding scenarios.

**Prospective commit:**
`refactor(compiler): resolve one checked compilation request`

### Work unit 2 – plan-owned scalar residual and data assembly

**Outcome:** selected scalar residual/data contribution records, not a
base/general constructor flag, drive volume and Robin assembly.

**Boundary:** residual/data part of `C2-R2`. Complete P5.1 truthful data spaces
and typed boundary policies before or in directly prerequisite commits. Do not
change observations, metrics, constraints, or specialized target strategies.

**Likely files:**

- `docs/implementation/v1/semantic-compiler.md`;
- `include/nmopt/compiler/v1/dealii_scalar_plan.hpp`;
- `include/nmopt/compiler/v1/dealii_fixed_dirichlet.hpp`;
- `include/nmopt/compiler/v1/dealii_compiler.hpp`;
- `tests/semantic_v1_contract.cc`; and
- `tests/dealii_diffusion_contract.cc`.

**Focused verification:** base scalar, general scalar/Robin, C5.6 composition,
term isolation, nonsymmetric transpose, and scalar-plan negative scenarios.

**Prospective commit:**
`refactor(lowering): assemble scalar residuals from contributions`

### Work unit 3 – plan-owned objective and service recombination

**Outcome:** observation/loss, metric, optional constraint, and transformation
records construct their bounded services and an explicit recombination proves
the C2 exit.

**Boundary:** remaining `C2-R2`. Preserve specialized boundary-control,
continuous-control, coefficient, gauge, and transposition strategies.

**Likely files:**

- `docs/implementation/v1/semantic-compiler.md`;
- `include/nmopt/compiler/v1/dealii_scalar_plan.hpp`;
- `include/nmopt/compiler/v1/dealii_fixed_dirichlet.hpp`;
- `include/nmopt/compiler/v1/dealii_compiler.hpp`;
- `include/nmopt/compiler/v1/compiled_problem.hpp`;
- `tests/semantic_v1_contract.cc`; and
- `tests/dealii_diffusion_contract.cc`.

**Focused verification:** fixed reconstruction plus subdomain tracking, general
scalar residual with compatible observation variants, metric/box independence,
handler removal, and reduced Taylor scenarios.

**Prospective commit:**
`refactor(lowering): compose scalar objective and services`

### Work unit 4 – lossless structured compilation provenance

**Outcome:** one resolved decision and the realized services populate a
versioned manifest that distinguishes all executable inputs; prose is rendered
from typed records.

**Boundary:** `C1-R2`. Consume the typed P5.1–P5.4 policy repairs that have
landed; do not invent policies for unselected targets.

**Likely files:**

- `docs/implementation/implementation-readiness-review.md`;
- `docs/implementation/v1/semantic-compiler.md`;
- `include/nmopt/compiler/v1/compiled_problem.hpp`;
- `include/nmopt/compiler/v1/dealii_types.hpp`;
- `include/nmopt/compiler/v1/dealii_compiler.hpp`;
- the bounded strategy headers that report realized choices; and
- `tests/dealii_diffusion_contract.cc`.

**Focused verification:** exact manifest table for every retained target plus
coefficient, bound, mesh, solve-policy, label-independence, and renderer
consistency differentials.

**Prospective commit:**
`fix(manifest): record resolved compilation decisions`

### Work unit 5 – common realized map and space records

**Outcome:** every bounded transformation/observation reports its executable
source/output realization and actual dimension through one path.

**Boundary:** `C1-R3`, `P5.2-R3`, `P5.3-R3`, and `P5.4-R3`. Do not change the
selected residual, objective, or metric formulas.

**Likely files:**

- `docs/implementation/v1/semantic-compiler.md`;
- `include/nmopt/compiler/v1/compiled_problem.hpp`;
- `include/nmopt/compiler/v1/dealii_compiler.hpp`;
- the bounded observation/transformation model headers; and
- `tests/dealii_diffusion_contract.cc`.

**Focused verification:** ordinary volume observation, baseline and weighted
boundary trace, fixed/controlled transformed state, point sensors, normal
flux, and all P5.4 transformed-state registrations.

**Prospective commit:**
`fix(compiler): record realized map spaces`

## Verification gate

For each unit, run its focused named scenarios first. Then run the required
repository profiles:

```bash
cmake --preset debug-neutral
cmake --build --preset debug-neutral
ctest --preset debug-neutral

cmake --preset debug-dealii
cmake --build --preset debug-dealii --parallel 1
ctest --preset debug-dealii
```

Run the existing backend-neutral sanitizer profile for ownership and resolved-
record lifetime changes:

```bash
cmake --preset sanitize-neutral
cmake --build --preset sanitize-neutral
ctest --preset sanitize-neutral
```

Do not add another permanent build profile merely for this remediation. If a
new deal.II ownership change needs instrumented verification, first document
the concrete gap and dependency compatibility, then follow the build
convention's policy for adding a demonstrated profile.

Before requesting a commit, inspect the staged diff and statistics, confirm
that unrelated remediation units are unstaged, and run `git diff --check`.

## Acceptance checklist

### C1 acceptance

C1 can return to `completed` only when:

- [x] every predictable public compiler input uses the documented diagnostic
      channel;
- [x] ordinary and specialized `Function` bindings validate field shape from
      one resolved binding request;
- [x] owned/detached session and solve-report behavior remains covered;
- [x] state and adjoint policy changes alter actual invocation and exact
      structured provenance;
- [x] the manifest is populated from resolved decisions and realized services,
      not target booleans or prose;
- [x] coefficient, bound, `Function`, and mesh provenance distinguishes
      materially different compiled products;
- [x] every realized map owns its source/output dimension, ordering, and
      pairing record;
- [x] compatibility prose is rendered from typed records;
- [x] exact structured manifest tests cover every retained target; and
- [x] the neutral, deal.II, and applicable sanitizer gates pass.

### C2 acceptance

C2 can return to `completed` only when:

- [x] one resolved closed request drives validation, registration, planning,
      model construction, mesh checks, and provenance;
- [x] raw graph flags no longer duplicate those decisions;
- [x] scalar handlers contribute typed executable payloads, not only handler
      names;
- [x] selected residual/data contributions drive fused assembly;
- [x] selected observations/losses, transformations, metrics, and constraints
      construct their bounded services from the plan;
- [x] specialized strategies remain explicit and do not enter the scalar plan
      by fallback;
- [x] one tested recombination changes residual or observation independently
      without a new complete problem class;
- [x] unsupported cross-products receive one exact registration diagnostic;
- [x] current value/JVP/VJP, reconstruction, nullspace, metric, constraint,
      independent-oracle, transpose, and Taylor tests pass;
- [x] manifests are structurally equal for every unchanged component; and
- [x] the roadmap alone is updated to mark the C1/C2 gates complete.

**Post-audit closure evidence:** the pre-construction decision owns the typed
region, requirement, component-inventory, and map-skeleton records. Its
finalizer consumes the closed request and realized services without receiving
the raw semantic graph, target enum, or a duplicate registration, and
`make_manifest()` remains a pure decision projection. Full compatibility is
unchanged under declaration-order permutations and compiler-policy display-
prose edits. The final `debug-neutral` (12/12), `debug-dealii` (41/41), and
`sanitize-neutral` (12/12) gates pass.
