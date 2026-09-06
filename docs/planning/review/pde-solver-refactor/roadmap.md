# PDE–solver boundary refactor roadmap

## Status and scope

This is the operative implementation plan derived from the
[assessment](assessment.md) and governed by the
[PDE–solver boundary](../../../design/pde-solver-boundary.md).

The roadmap owns only the PDE/formulation/compiler numerical-boundary cleanup
and the directly required application/documentation migrations. It does not
restart paused Chapter 5/6 feature work, expand the benchmark catalogue, or
redesign the runner and post-processing stack.

This programme also does not introduce distributed/MPI backend support,
standardize every linear solver or preconditioner, or make the compiler lower
every semantically valid component combination. Such work requires a separate
selected requirement.

The refactor starts from the frozen pre-refactor application milestone. The
current work branch is `codex/refactor/pde-solver-boundary`. Update this
document's current handoff as units are accepted; do not copy mutable status
into design or assessment documents.

## Target outcome

The completed system has one small solver-facing boundary and two independent
producers:

```text
existing deal.II application -> callback adapter -------\
                                                       +-> ExecutableModelT
ProblemSpec -> validation -> lowering -> typed numerics-/   + solve callbacks
                                                           + metric/capabilities
                                                                    |
                                                                    v
                                                               formulation
                                                                    |
                                                                    v
                                                               optimizer
```

The native compiler path additionally retains typed numerical/application state
for reconstruction, diagnostics, dimensions, and output. It does not recover
that state through `dynamic_cast` from the erased executable interface.

## How to execute this roadmap with Codex

### Read the smallest sufficient context

For one unit, read only:

1. `.agents/README.md` and the routed workflow/Git/code/build/documentation
   instructions required by that unit;
2. the [boundary design](../../../design/pde-solver-boundary.md) sections named
   by the unit;
3. only the [architecture-map](architecture-map.md) sections named by the unit;
4. this roadmap's unit;
5. only the assessment findings named by the unit;
6. only the deletion-ledger entries named by the unit;
7. the current production files and focused tests listed by the unit or found
   by a current consumer search.

Do **not** read the complete assessment by default during implementation. The
assessment is an evidence archive. Follow links/findings when the current unit
exposes a tradeoff or contradiction.

### Prepare a current-tree execution brief before editing

At the start of each implementation unit, inspect the current branch and write
a compact working brief in the assistant response or handoff containing:

```text
Unit:
Current branch/head:
Findings being addressed:
Current non-test consumers:
Files expected to change:
Behavior/invariants preserved:
Implementation decision still open:
Decision criteria:
Deletion enabled by this unit:
Focused verification:
Full verification gate:
Prospective commit:
```

The brief is working context, not normally a permanent repository document.
If current evidence invalidates the roadmap boundary or reveals a new public
architectural decision, stop before editing, update the plan, and obtain review
where required by `.agents/workflow.md`.

### Decision rule

Roadmap units intentionally separate **settled architecture** from **local
implementation choices**.

A unit must not invent a public hierarchy merely because the roadmap names a
responsibility such as “state coordinates” or “typed native realization.”
Choose among a value, concrete class, free functions, callback bundle, or
existing type using current evidence.

A new abstraction is justified only when at least one of these is true:

- two demonstrated production consumers need the same responsibility; or
- the abstraction enables a named deletion larger than itself.

Prefer the smallest representation that leaves ownership explicit and makes
the containing complete-model/dispatch code shrink.

### Verification order

For each unit:

1. add/retain focused characterization where the migration could otherwise
   erase an oracle;
2. run the smallest focused scenario(s);
3. run the applicable full Debug pipeline;
4. run sanitizer/other broad checks when ownership or memory lifetime changes;
5. inspect the production diff separately from tests/documentation; and
6. record production lines added/deleted and the named deletion achieved or
   enabled.

Use the repository's current `.agents/build.md` commands rather than copying
stale commands from this roadmap if the build interface changes.

## Global invariants

Every unit must preserve these unless the unit explicitly changes a documented
contract:

1. Lagrangian sign and primal/covector conventions.
2. Exact residual value/JVP/VJP mathematical actions.
3. `ReducedDTOT` first-order DTO orchestration and evaluation-token safety.
4. Objective regularization remains distinct from search metric.
5. State/adjoint solve reports remain explicit.
6. Formulation/optimizer headers remain deal.II-independent.
7. External applications are not required to use `ProblemSpec`, compiler
   manifests, recipes, scenarios, or the runner.
8. Native field output remains owned by an application/typed backend boundary,
   not a generic optimizer/formulation contract.
9. No universal PDE-family, output, observation, residual, or linear-solver
   hierarchy is added.
10. Existing supported semantic graphs remain valid unless a separately
    accepted correctness change says otherwise.

## Global simplification rules

- Prefer deletion over relocation.
- A temporary code-positive compatibility seam must name the unit that removes
  it.
- Do not split a large class/header unless responsibilities actually move or
  disappear.
- Preserve fused FE assembly when useful; component ownership does not imply
  separate quadrature loops.
- Keep one in-memory owner for each realized compiler decision or fact.
- Keep serialization compatibility at the artifact edge.
- Transfer unique tests/oracles before deleting a legacy implementation.
- Report production C++/Python, tests, and documentation separately.
- The completed cleanup is expected to be net-negative in production code.
- Moving code between files, renaming types, or introducing factories does not
  count as simplification unless it removes duplicated responsibility, state,
  dispatch, or ownership.
- A permanent compatibility facade requires a demonstrated non-test/external
  consumer. Otherwise it must have a named removal gate.

## Comments and documentation during implementation

Documentation changes follow the responsibility being changed.

- Update a nearby comment or documentation comment in the same unit when the
  code change would otherwise make it false or materially misleading.
- Add comments only where they preserve information that is not obvious from
  the code: mathematical meaning or sign, ownership/lifetime, numerical
  assumptions, responsibility/non-responsibility, or the reason for an unusual
  implementation choice.
- Do not narrate straightforward control flow or restate type/function names.
- Do not perform repository-wide comment or docstring normalization during a
  numerical refactor unit.
- Broader stale-comment and documentation cleanup belongs to F2 after the
  implementation shape has stabilized.
- If a comment policy proves generally useful beyond this refactor, consider
  promoting it to the repository code conventions during F2 rather than
  changing agent-wide policy in advance.

## Work sequence

The roadmap is organized into phases so a Codex session can operate on one
bounded concern without loading the whole program.

```text
R0  architecture audit documentation
  |
  v
A   prove the external solver-facing boundary
  A1 callback executable
  A2 unchanged deal.II application proof
  |
  v
B   repair native typed ownership
  B1 retain native typed result; remove application downcasts/output recovery
  B2 centralize mechanical service packaging where now demonstrably common
  |
  v
C   reduce complete numerical targets by demonstrated components
  C1 state coordinates/reconstruction
  C2 decision realization/coupling
  C3 observation/objective/regularization ingredients
  C4 solve-policy composition where target code remains duplicated
  |
  v
D   simplify compiler decisions and provenance
  D1 one closed lowering/realization decision
  D2 one manifest/provenance owner
  |
  v
E   retire redundant implementations
  E1 direct v0 scalar path
  E2 scalar-specific KKT/remaining compatibility paths
  |
  v
F   close documentation and integration proof
```

C1–C4 are evidence-gated. Do not create a component merely to satisfy the
roadmap label. If a proposed extraction does not serve multiple current
consumers or enable a larger deletion, record the negative decision and leave
the responsibility local.

D1 may begin earlier than some C units if current compiler edits are clearly
independent. Do not reorder B1 after E1: direct-model retirement must not happen
while native applications still depend on concrete executable recovery.

---

## R0 — Record the architecture audit

**Status:** documentation preparation/review.

**Outcome:** The refactor has one accepted long-lived boundary, one detailed architecture
atlas, one evidence archive, one deletion ledger, and one executable roadmap.

**Findings:** PSR-001 through PSR-024.

**Changes:**

- add `docs/design/pde-solver-boundary.md`;
- add this review package under
  `docs/planning/review/pde-solver-refactor/`;
- route the package from `docs/README.md` and
  `docs/planning/review/README.md`;
- add `architecture-map.md` as the durable R0 implementation mental model;
- do not mass-update implementation/reference/application documents yet.

**Preserve:** Existing design authorities remain authoritative for the
mathematical/semantic model; this document adds the implementation ownership
boundary rather than replacing the interface specification.

**Do not:** edit C++, tests, build configuration, runner behavior, or generated
run evidence.

**Verification:** Markdown rendering/link review and `git diff --check` after
actual repository insertion.

**Exit:** User accepts the architecture direction and roadmap scope.

**Prospective commit:** `docs(refactor): record PDE solver boundary audit`

---

# Phase A — prove the external boundary

## A1 — Add a callback-backed executable adapter

**Outcome:** An application can provide the existing executable operations
through callables without making its PDE class derive from an nmopt interface.

**Findings:** PSR-001, PSR-002.

**Required reading:** boundary sections “Stable solver-facing surface” and
“Existing deal.II application path”; ledger entry `ExecutableModelT`.

**Expected production surface:** `include/nmopt/contract/` plus focused contract
tests. Re-inspect current paths before editing.

**Preserve:**

- exact `ExecutableModelT` public virtual surface;
- layout validation and primal/covector roles;
- ordinary callable capture lifetime semantics;
- existing `ReducedDTOT` constructors/lifetime behavior.

**Required implementation behavior:**

- accept variable and test layouts;
- forward residual, JVP, and VJP callables;
- forward objective value and derivative callables;
- validate required callables/layouts at construction or first use according to
  current contract conventions;
- allow residual and objective implementation to be replaced independently.

**Local decision:** Use the smallest backend-neutral representation: one small
concrete `ExecutableModelT` implementation, a factory around it, or a compact
callback aggregate. Do not add another abstract interface.

**Decision criteria:**

- complete behavior readable in one small implementation surface;
- no deal.II/compiler/application types;
- no builder hierarchy;
- no output or linear-solver callback;
- construction validation is testable;
- no duplicated ownership protocol beyond ordinary captures.

**Deletion enabled:** later complete-model inheritance/adapter code may be
removed after producer paths migrate.

**Focused verification:** construction validation; forwarding of all five
operations; JVP/VJP pairing on a simple dense/reference model; complete reduced
evaluation; objective replacement without residual change.

**Full gate:** backend-neutral Debug pipeline.

**Prospective commit:** `feat(contract): add callback executable adapter`

---

## A2 — Prove an unchanged deal.II application can use nmopt

**Outcome:** A small tutorial-style deal.II PDE remains independently runnable
and is adapted to an existing nmopt reduced formulation/optimizer without using
the semantic compiler.

**Findings:** PSR-002, PSR-016, PSR-018.

**Required reading:** boundary sections “Existing deal.II application path”,
“State and adjoint solves”, and “Application and output boundary”.

**Expected files:** smallest suitable existing deal.II test fixture/executable
plus one integration translation unit. Re-inspect current test layout before
choosing a new executable.

**Fixture requirements:**

- own mesh, FE/DoFs, assembly, state/adjoint solves, and native output;
- original PDE source contains no nmopt solver/compiler/semantic/application
  dependency and preferably no nmopt include;
- standalone forward solve remains testable independently;
- nmopt integration occurs in one readable adapter/application unit;
- use A1, `StateAdjointSolversT`, an existing metric, `ReducedDTOT`, and an
  existing optimizer.

**Required mathematical proof:**

- state residual at solved state;
- residual JVP finite-difference check;
- JVP/VJP pairing;
- state-recomputed reduced Taylor derivative check;
- optimization convergence for the selected simple problem;
- objective replacement without PDE-source modification;
- native output still written by the fixture/application.

**Local decision:** Reuse an existing test executable if readability remains
high; create a dedicated executable only when isolation materially improves the
proof.

**Do not:** introduce `ProblemSpec` conversion, a tutorial library framework,
a universal PDE interface, or a second runner.

**Deletion enabled:** validates the external boundary required before legacy
integration assumptions are removed.

**Full gate:** deal.II Debug pipeline plus backend-neutral Debug if shared
contract code changed.

**Prospective commit:** `test(dealii): prove external application adapter`

---

# Phase B — repair native typed ownership

## B1 — Retain typed native realization and remove application recovery downcasts

**Outcome:** B1/B2 native application code receives the typed information needed
for dimensions, diagnostics, and output directly from compilation/application
composition. No application path recovers private concrete targets through
`dynamic_cast` from `ExecutableModelT`.

**Findings:** PSR-012, PSR-016, PSR-018, PSR-019.

**Ledger entries:** `CompiledProblemT`, B1/B2 downcasts, model-owned native
output, B1/B2 execution adapters.

**Required reading:** boundary sections “The missing middle layer” and
“Application and output boundary”. Read only B1/B2/compiler files participating
in current downcast/output flow.

**Characterization before migration:**

- inventory every native `dynamic_cast` from an executable model;
- record the exact information each cast recovers;
- identify whether each fact is solver-port data, compiler provenance, typed
  numerical state, or application-only output state;
- record current native filenames, field identities, topology, dimension
  evidence, and objective-component evidence used by tests/artifacts.

**Target behavior:**

```text
compiler/lowerer constructs typed realization
        |\
        | +-> native/app information retained explicitly
        |
        +-> ExecutableModelT erased solver view
```

**Local decision:** Choose the smallest mechanism that retains typed/native
state. Candidate shapes include a compiler-private typed result, application
closures built before erasure, or a product carrying generic services plus a
native view. Do not preselect a universal `NumericalRealisation` public class.

**Decision criteria:**

- all current downcasts disappear;
- `ExecutableModelT` gains no output/dimension/diagnostic methods;
- application/compiler lifetime remains explicit;
- current compiler path still produces `ReducedDTOT` cleanly;
- mechanism can represent at least B1 and B2 without target-string dispatch;
- production code does not increase unless the temporary seam names deletion in
  B/C/E.

**Output migration:** Move filesystem/native writer responsibility outward only
as far as required to eliminate executable-model recovery. Share writer code
only where current B1/B2 consumers prove the same operation.

**Preserve:** output filenames/topology/field semantics and current experiment
artifact values unless they were merely private implementation labels.

**Focused verification:** B1/B2 application contract scenarios, native output
checks, compiler lifetime scenario.

**Full gate:** full deal.II Debug pipeline; sanitizer coverage if ownership or
shared lifetime changes.

**Deletion target:** all B1/B2 executable-model downcasts and model-owned writer
methods whose callers migrate in this unit.

**Prospective commit:** `refactor(application): retain typed numerical ownership`

---

## B2 — Centralize mechanical compiler service packaging

**Outcome:** Concrete lowering results are converted to executable/solve/metric
capability packaging by one mechanical path where behavior is identical.

**Findings:** PSR-015, PSR-016.

**Ledger entries:** `StateAdjointSolversT`, `CompiledProblemT`, complete model
entries.

**Precondition:** B1 has made typed/native ownership explicit enough that a
helper will not re-create early type erasure.

**Current audit:** enumerate repeated target branches that:

- assign executable;
- capture model in state/adjoint solve lambdas;
- attach metric/constraint/Hessian;
- attach lifetime owner.

**Local decision:** Prefer a private free/template helper or compact value over
a common complete-model base class. If target-specific solve signatures/policy
are materially different, centralize only the repeated outer packaging and keep
specialized callables local.

**Decision criteria:**

- adding a compatible lowering target no longer copies the same lambda/package
  boilerplate;
- no new public abstraction;
- no loss of typed/native application state;
- net-negative production C++ for the completed unit.

**Focused verification:** representative scalar volume, continuous control,
Dirichlet control, Neumann control, coefficient identification, plus compiler
lifetime.

**Full gate:** deal.II Debug pipeline.

**Deletion target:** repeated model-capturing solve/service packaging blocks.

**Prospective commit:** `refactor(compiler): centralize compiled service wiring`

---

# Phase C — reduce complete numerical targets

Phase C is deliberately incremental. Each sub-unit must prove that a numerical
responsibility is shared by real current consumers or enables a named deletion.
A negative decision is acceptable: document that the current specializations
remain local and continue.

## C1 — Unify demonstrated state-coordinate/reconstruction machinery

**Outcome:** Fixed-Dirichlet independent-state reconstruction has one reusable
implementation wherever at least two current targets use equivalent $P$, fixed
lifting, tangent embedding, and dual pullback semantics.

**Findings:** PSR-003, PSR-007, PSR-008.

**Primary candidates:** `ScalarComponentModel` and
`DirichletControlLiftingModel`; compare full-coordinate targets only after the
first extraction is understood.

**Characterization:**

- record state physical and solver-facing dimensions;
- record construction of $P$, $P^T$, fixed lifting, independent DoF maps;
- record handling of hanging/Dirichlet constraints;
- record residual/objective/native-output callers;
- compare the two implementations mathematically before sharing code.

**Local decision:** Select concrete value/class/free-function-owned data based on
current common state. No public `StateModel`/`StateCoordinatesT` hierarchy is
required by the roadmap.

**Decision criteria:**

- at least two production consumers;
- one owner for reconstruction matrices/maps;
- value, tangent, and dual pullback semantics are explicit;
- Dirichlet control can compose $L_D$ without forcing every state coordinate
  representation to know about a decision variable;
- containing target classes shrink;
- net deletion after both consumers migrate.

**Do not:** force coefficient identification, pure Neumann gauge, or every
full-coordinate target into this representation merely for uniformity.

**Focused verification:** reconstruction value/JVP/VJP identities, state solve,
objective derivative for both migrated targets, native reconstruction output.

**Full gate:** relevant deal.II target scenarios + full deal.II Debug.

**Deletion target:** duplicated reconstruction/independent-DoF machinery in the
migrated complete targets.

**Naming decision:** choose a responsibility-based name only after the concrete
surviving API is known.

**Prospective commit:** `refactor(dealii): share state reconstruction machinery`

---

## C2 — Separate decision discretization from PDE influence where it deletes duplication

**Outcome:** Current control/parameter implementations expose reusable numerical
ingredients without inventing one universal runtime decision hierarchy.

**Findings:** PSR-005, PSR-006, PSR-011.

**Primary comparison:** DG0 volume control, continuous volume control, facewise
Neumann control, continuous Neumann control. Treat Dirichlet lifting and
coefficient parameter as distinct mechanisms.

**Characterization:** for each candidate record:

- optimizer coordinates/layout;
- physical/backend coordinate embedding;
- PDE influence action and transpose/pullback;
- mass/stiffness/trace operators assembled;
- constraint support;
- regularization and metric consumers;
- output/diagnostic consumers.

**Target direction:**

- affine forcing controls may share an internal `apply`/transpose coupling
  concept if it removes duplicated code;
- Neumann runtime realization polymorphism remains justified;
- Dirichlet control remains a lifting/reconstruction mechanism;
- coefficient identification remains parameterized-operator machinery.

**Local decisions:**

1. whether volume-control coupling deserves a concrete reusable component;
2. how far to narrow `NeumannControlRealisation`;
3. whether shared decision-space mass/stiffness operators should be retained as
   immutable numerical ingredients rather than exposed through metric methods.

**Decision criteria:** every new component has multiple production consumers or
names a larger deletion; complete target/model APIs lose methods/state;
constraint capability remains optional and realization-specific.

**Do not:** create a public `DecisionRealisation` base class covering coupling,
lifting, and parameter derivatives.

**Focused verification:** residual/JVP/VJP pairing for each migrated influence;
box projection where supported; control/parameter layout identities.

**Full gate:** representative deal.II scenarios + full deal.II Debug.

**Deletion targets:** direct matrix/coupling duplication in complete targets;
unrelated regularization/metric/output methods removed from a narrowed Neumann
realization when callers have migrated.

**Prospective commit:** choose a scope matching the actual bounded extraction;
do not combine unrelated decision families merely to fit C2.

---

## C3 — Separate objective ingredients from complete PDE targets

**Outcome:** State tracking and decision regularization are composed as objective
ingredients rather than being inseparable responsibilities of complete PDE
models, where current reuse justifies extraction.

**Findings:** PSR-009, PSR-010.

**Primary evidence:** `VolumeObservationAssembly`, embedded scalar tracking
assembly, point/boundary/flux observation ports, continuous/Dirichlet control
regularization.

**Characterization:**

- identify explicit observation value/JVP/VJP consumers outside objective
  evaluation;
- identify compiled quadratic $Q, q, c$ owners;
- identify regularization operators and their underlying decision-space
  matrices;
- identify reduced-Hessian consumers of $Q$ and regularization curvature.

**Target direction:**

```text
physical state observation -> tracking/loss realization -> J/J'
decision operators ---------> regularization -----------> J/J'
shared decision operators ---> metric ------------------> optimizer
```

**Local decisions:**

- reuse `VolumeObservationAssembly` as-is versus a modest generalization;
- whether a small quadratic-form value/helper deletes enough repeated code;
- whether regularization remains free functions/callbacks or deserves a named
  internal value;
- whether any $H^{1/2}$ neutral-operator extraction is justified now.

**Decision criteria:** do not materialize observation vectors where quadratic
assembly is the natural representation; do not add a public `ObservationBase`
or `RegularizationT`; retain explicit measurement ports when diagnostics/tests
need them.

**Assembly rule:** preserve fused FE traversals when efficient. Separate
component ownership from traversal scheduling.

**Focused verification:** objective value/derivative; explicit sensor/trace/flux
JVP/VJP where migrated; regularization/metric independence; current Hessian
checks if shared ingredients move.

**Full gate:** relevant deal.II scenarios + full deal.II Debug.

**Deletion target:** duplicated embedded tracking/regularization code and
complete-model objective methods/state that become pure composition.

**Prospective commit:** scope to the actual extracted objective component(s).

---

## C4 — Move target-specific solver choice to the state/adjoint service seam where useful

**Outcome:** Concrete state/adjoint solve callbacks own inversion policy without
forcing complete PDE models to also be solver-policy containers.

**Findings:** PSR-008, PSR-015.

**Precondition:** C1–C3/B2 have made numerical ingredients accessible without
requiring a complete model object for every operation.

**Characterization:** compare:

- SPD fixed operator solve;
- nonsymmetric direct/transpose solve;
- mean-zero augmented solve;
- point-dependent coefficient solve.

**Local decision:** extract only repeated policy/composition logic. Specialized
mean-zero and parameter-dependent behavior may remain concrete callables.

**Decision criteria:**

- formulation continues to see only `StateAdjointSolversT`;
- state discretization exposes mathematical/numerical facts rather than an
  optimizer policy;
- no universal linear-solver base class;
- target classes shrink or branches disappear.

**Focused verification:** solve reports, state residuals, adjoint identities,
nonsymmetric transpose path, mean-zero compatibility, coefficient positivity.

**Full gate:** full deal.II Debug; sanitizer if solve-service lifetime changes.

**Deletion target:** repeated complete-model solver-selection branches that have
become composition logic.

**Prospective commit:** `refactor(dealii): compose state and adjoint solve services`

---

# Phase D — simplify compiler decisions and provenance

## D1 — Replace parallel target identities with one closed lowering decision

**Outcome:** Semantic closure/lowerability/construction consume one authoritative
realization decision rather than parallel target enums, booleans, predicates,
and conversion switches.

**Findings:** PSR-004, PSR-013, PSR-014.

**Ledger entries:** `ResolvedTargetFamily`, `CompiledTargetKind`, compiler
`uses_*` booleans, `ScalarLoweringPlan`.

**Characterization:** before editing, list every current field/function that
classifies the complete target and every consumer in validation, construction,
manifest, application, and tests.

**Settled direction:** the closed decision should preserve independently
meaningful numerical axes rather than create another enum of complete target
class names.

**Local decision:** exact C++ representation of the closed plan. Possible tools
include nested values/variants and existing `ScalarLoweringPlan`. Do not assume
one variant per semantic node or one class per target.

**Decision criteria:**

- graph/realization is classified once;
- support/lowerability rejection remains explicit;
- later construction does not reconstruct a target family from booleans;
- no string dispatch or complete-target registry replaces enum duplication;
- all current supported target diagnostics remain explainable;
- production diff is net-negative or names immediate deletion enabled in D2/E.

**Focused verification:** semantic lowerability, all current compiler target
construction scenarios, manifest target/provenance assertions.

**Full gate:** backend-neutral semantic pipeline + full deal.II Debug.

**Deletion target:** `CompiledTargetKind`, conversion switches, redundant
whole-target predicates/booleans, duplicate registration state made obsolete by
the closed plan.

**Prospective commit:** `refactor(compiler): resolve one lowering decision`

---

## D2 — Keep one manifest/provenance owner

**Outcome:** Every realized compiler fact has one in-memory owner; human-readable
and persisted compatibility fields are projections at the serialization edge.

**Findings:** PSR-017, PSR-024.

**Ledger entries:** `ResolvedCompilationDecision`, `CompilationManifest`,
`CompiledCompatibilityView`, artifact manifest fields.

**Consumer audit required:**

- application adapters;
- experiment envelope;
- artifact writer;
- runner output;
- Python post-processing/readers;
- focused manifest tests.

For each legacy field record whether it is:

```text
authoritative semantic/realized fact
current persisted contract
human-readable convenience
internal-only compatibility copy
unused
```

**Local decision:** retain `ResolvedCompilationDecision` as source of truth or
replace it only if D1 produced a clearly smaller structured owner.

**Decision criteria:**

- no realized value stored twice solely for compatibility;
- persisted compatibility remains stable for active consumers;
- projection happens once at artifact/serialization boundary;
- equality/copy checks that existed only to synchronize representations are
  deleted.

**Do not:** break current artifact readers merely to clean internal C++.
Persisted field removal requires its own consumer-backed decision.

**Focused verification:** manifest contract scenarios, B1/B2 artifact contract,
post-processing reader tests for touched fields.

**Full gate:** backend-neutral + deal.II Debug plus relevant Python/tool tests.

**Deletion target:** compatibility view/flat fields/projection code proven
redundant by consumer audit.

**Prospective commit:** `refactor(manifest): keep one realized provenance owner`

---

# Phase E — retire redundant implementations

## E1 — Retire the direct v0 scalar implementation

**Outcome:** The direct v0 deal.II scalar complete model is no longer a parallel
production implementation. Its unique mathematical/reference behavior lives in
the surviving numerical path and tests.

**Findings:** PSR-003, PSR-020.

**Ledger entry:** `ScalarDiffusionReactionModel<dim>`.

**Preconditions:**

- B1 removed application downcast/output dependence;
- surviving native/compiler services cover canonical scalar execution;
- C/D migrations required for unique v0 behavior are complete;
- external application proof A2 is accepted.

**Mandatory oracle transfer before deletion:**

- hand/reference weak-form values;
- residual value/JVP/VJP tests;
- state/adjoint solve identity;
- reduced derivative Taylor test;
- metric/constraint behavior;
- reduced Hessian oracle used by current supported methods;
- supplied OTD behavior still required by non-test/test consumers.

**Consumer audit:** search every production include/construction/downcast and
classify each as migrate, delete, or bounded retain decision.

**Local decisions:**

- exact destination for supplied-OTD construction;
- whether any v0-only convenience still supplies independent public value.

**Stop condition:** If one v0-only capability has a real non-test consumer and
cannot migrate with net simplification, stop and present that capability as a
bounded retain/delete decision. Do not keep the entire model by default.

**Deletion target:**

- `include/nmopt/dealii/scalar_diffusion_reaction.hpp` if gate passes;
- direct-volume compiler branch and v0-only construction helpers;
- concrete-model output/downcast remnants;
- comparison tests whose sole purpose was dual-path equivalence after the
  independent oracle has moved.

**Focused verification:** canonical scalar volume, supplied OTD, Hessian,
metric/constraint, B1, compiler lifetime.

**Full gate:** backend-neutral + deal.II Debug + backend-neutral sanitizer; add
other profiles required by current build instructions.

**Prospective commit:** `refactor(dealii)!: retire direct scalar reference path`

---

## E2 — Remove scalar-specific KKT and remaining compatibility paths where redundant

**Outcome:** KKT and other legacy compatibility code survives only when it has
independent demonstrated behavior.

**Findings:** PSR-021 plus residual PSR-013/017/020 findings after E1.

**Consumer audit:** `ScalarDiffusionReactionKKT`, v0 comparison adapters,
version-specific target labels, compatibility facades introduced during prior
units.

**Target direction:** prefer the backend-neutral quadratic KKT product and
serial adapter when they cover the current behavior.

**Local decision:** retain a scalar-specific implementation only if current
non-test consumers or an independent numerical capability cannot be expressed
through the generic product without losing clarity/correctness.

**Deletion gate:** transfer unique KKT oracle/solve tests before deletion.

**Focused verification:** KKT product/solve, supplied OTD/KKT if still linked,
PDAS consumers where relevant.

**Full gate:** full applicable Debug pipelines.

**Prospective commit:** scope to the actual deletion, for example
`refactor(kkt): remove redundant scalar KKT realization`.

---

# Phase F — close the boundary and documentation

## F1 — Add the final integration guide

**Outcome:** An external deal.II user can understand the integration boundary
without reading compiler internals or historical roadmap documents.

**Inputs:** the accepted A2 fixture and final public API.

**Changes:** add one concise guide under `docs/guides/` showing:

- standalone PDE ownership;
- callback executable construction;
- state/adjoint solves;
- metric/optional capability;
- `ReducedDTOT` and optimizer invocation;
- application-owned output;
- lifetime rules;
- derivative verification expectations.

**Do not:** document speculative adapters that were not implemented.

**Verification:** links/code snippets match the tested API.

**Prospective commit:** `docs(dealii): document external solver integration`

---

## F2 — Repository documentation consistency sweep

**Outcome:** Current documentation describes the implemented post-refactor
architecture without retaining stale parallel narratives.

**Finding:** PSR-024.

**Search targets:**

- deleted/narrowed `*Model` names;
- `dynamic_cast`-based output descriptions;
- old v0/v1 production path claims;
- removed target enums/booleans/compatibility views;
- old manifest ownership claims;
- obsolete roadmap handoffs/status tables;
- documentation that claims component lowering behaves differently from final
  code;
- output ownership and external-integration instructions;
- C++ comments/documentation comments that describe deleted model ownership,
  obsolete v0/v1 paths, old target identity, or model-owned output;
- comments whose type/responsibility terminology no longer matches the
  surviving numerical decomposition.

**Authority rule:**

- design docs describe long-lived boundaries;
- implementation docs describe exact implemented capability;
- reference docs describe current public/application contracts;
- planning roadmap owns mutable work status;
- review assessment remains historical evidence and is not rewritten to pretend
  it observed the final state.

**Local decision:** update, supersede, or delete stale docs according to their
current authority and non-historical value. Do not preserve duplicate status
ledgers.

**Verification:** repository Markdown links/rendering, `git diff --check`, and a
manual authority pass through `docs/README.md`.

**Prospective commit:** `docs(refactor): align documentation with solver boundary`

---

## F3 — Final architecture proof and deletion accounting

**Outcome:** The program closes with measured evidence rather than only passing
tests.

**Verify:**

- external application and native compiler both reach unchanged formulation/
  optimizer contracts;
- no application concrete-model recovery from `ExecutableModelT`;
- objective callbacks and state/adjoint solve callbacks remain independently
  replaceable without changing the PDE implementation;
- a first-order integration does not need to provide `ReducedHessianT` or
  other second-order capability;
- no model-owned filesystem output remains unless explicitly retained with a
  documented reason;
- one authoritative compiler realization decision;
- one in-memory owner for each realized manifest fact;
- no unplanned temporary compatibility seam remains;
- unique mathematical oracles survive;
- full current build/test gates pass.

**Record in roadmap handoff/review:**

```text
production C++ added/deleted
production Python added/deleted
tests added/deleted
documentation added/deleted
complete model/types/files deleted
dispatch branches deleted
duplicate representations deleted
remaining dynamic_casts and reason
remaining version-specific production code and reason
remaining deferred cleanup candidates
```

The production cleanup is expected to be net-negative. If it is not, review
which abstractions failed to enable their named deletions before declaring the
program complete.

**Prospective commit:** normally combined with the final documentation/accounting
unit unless an independently meaningful cleanup remains.

## Current handoff

```text
Completed: architecture exploration and draft audit package
Current:   R0 documentation review before repository commit
Next:      accept/revise R0, then A1 callback executable
Blocked:   no implementation blocker; R0 requires user review
```
