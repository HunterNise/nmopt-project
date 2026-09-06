# PDE–solver refactor deletion ledger

## Purpose

This ledger records what current production types, branches, and stored
representations are expected to survive, narrow, migrate, or disappear during
the PDE–solver boundary refactor.

It is not permission to delete an item merely because it is large or old. Every
deletion requires:

1. a non-test consumer check;
2. identification of unique behavior and mathematical oracles;
3. a migration destination for behavior that survives;
4. focused verification after migration; and
5. satisfaction of the stated deletion gate.

The [assessment](assessment.md) owns the evidence and finding IDs. The
[roadmap](roadmap.md) owns execution order and current status.

## Disposition vocabulary

| Disposition | Meaning |
| --- | --- |
| **Keep** | Current responsibility is already useful and appropriately placed. |
| **Keep / simplify** | Concept survives but duplicate state or wiring should shrink. |
| **Narrow** | Type or component has a useful core but owns unrelated responsibilities. |
| **Migrate then delete** | Unique behavior must move before the current container disappears. |
| **Audit then decide** | Consumer/behavior evidence is still required. |
| **Move outward** | Behavior survives at an application/serialization boundary rather than in the current numerical type. |

## Stable contracts and services

### `ExecutableModelT<Backend>`

- **Location:** `include/nmopt/contract/executable_model.hpp`
- **Disposition:** **Keep**.
- **Current role:** solver-facing layouts, residual, residual JVP/VJP,
  objective, objective derivative.
- **Non-test consumers:** reduced formulations/compiler products and current
  application/compiler execution paths.
- **Unique behavior:** nonlinear-capable mathematical facade; coefficient
  identification proves that fixed-matrix alternatives are insufficient.
- **Target:** unchanged public boundary; may gain a callback-backed concrete
  implementation/factory without widening the interface.
- **Deletion gate:** none.
- **Related findings:** PSR-001, PSR-002.

### `StateAdjointSolversT<Backend>`

- **Location:** formulation/executable contract layer.
- **Disposition:** **Keep**.
- **Current role:** callback boundary for state and adjoint solves plus reports.
- **Target:** remain the formulation inversion-of-control seam; simplify
  repeated compiler construction around it.
- **Deletion gate:** none.
- **Related findings:** PSR-001, PSR-008, PSR-015.

### `ReducedDTOT<Backend>`

- **Location:** `include/nmopt/contract/reduced_dto.hpp`
- **Disposition:** **Keep**.
- **Current role:** Discretize Then Optimize first-order reduced orchestration.
- **Unique behavior:** split value evaluation and derivative augmentation,
  evaluation-token safety, lifetime-owner support.
- **Target:** formulation remains independent of deal.II/compiler target types.
- **Deletion gate:** none.
- **Related findings:** PSR-001.

### `MetricT`, `ConstraintT`, `ReducedHessianT`

- **Disposition:** **Keep**.
- **Target:** generic optional capabilities; construction moves away from
  complete-model ownership where useful.
- **Do not:** merge objective regularization and metric roles.
- **Related findings:** PSR-010.

### serial backend/KKT/PDAS helpers

- **Representative locations:** `include/nmopt/dealii/serial_backend.hpp`,
  `serial_spd_solver.hpp`, serial KKT/PDAS helpers.
- **Disposition:** **Keep / simplify**.
- **Current role:** narrow concrete deal.II services around generic contracts.
- **Target:** continue as concrete backend services; deduplicate only proven
  repeated solve mechanics.
- **Related findings:** PSR-001, PSR-008.

## Compiler-path packaging

### `CompiledProblemT`

- **Location:** `include/nmopt/compiler/v1/compiled_problem.hpp`
- **Disposition:** **Keep / simplify**.
- **Current role:** compiler-path package of executable, solves, metric,
  optional constraint/Hessian, manifest, lifetime owner, and box data.
- **Good:** validates layout compatibility and constructs `ReducedDTOT` cleanly.
- **Target:** compiler-specific package, not universal external-application API.
  It may be reshaped once typed/native realization ownership is explicit.
- **Decision gate:** choose the smallest change that removes native
  erase/recover cycles without turning compiler provenance into an external
  integration requirement.
- **Related findings:** PSR-012, PSR-016.

### `CompilationResultT` nullable product bundle

- **Disposition:** **Audit then decide**.
- **Pressure:** conceptually behaves like a sum of compiled product variants but
  is represented through multiple nullable pointers.
- **Target option:** a real variant/sum type if it removes invalid states and
  branch duplication after higher-priority compiler cleanup.
- **Do not:** prioritize this before target identity and typed realization.
- **Related findings:** PSR-013, PSR-016.

## Compiler decision and manifest state

### `ResolvedTargetFamily`

- **Disposition:** **Keep / simplify or replace once**.
- **Current role:** resolved target-family identity.
- **Target:** one authoritative closed realization decision. It may remain as a
  field/variant within that decision if it still has independent value.
- **Deletion gate:** do not delete until all construction/lowerability/manifest
  consumers use the replacement directly.
- **Related findings:** PSR-013, PSR-014.

### private `CompiledTargetKind`

- **Disposition:** **Migrate then delete**.
- **Current role:** second representation of substantially the same target
  cross-product.
- **Migration:** make construction/provenance consume the retained closed
  decision directly.
- **Deletion gate:** no non-test consumer requires a distinct identity and all
  conversion switches/predicates are gone.
- **Related findings:** PSR-013.

### compiler `uses_*` target booleans and target predicates

- **Disposition:** **Migrate then delete or reduce to local derived queries**.
- **Current role:** repeated target classification and dispatch.
- **Target:** component facts may remain in a closed plan, but whole-target
  booleans should not be parallel authorities.
- **Deletion gate:** construction and lowerability no longer reconstruct target
  identity from boolean combinations.
- **Related findings:** PSR-013, PSR-014.

### `ResolvedCompilationDecision`

- **Disposition:** **Keep / simplify**.
- **Current role:** structured realized compiler decision.
- **Target:** likely source of truth for realized provenance unless the compiler
  plan refactor yields a smaller better owner.
- **Decision gate:** preserve actual persisted consumer needs, not legacy
  internal compatibility for its own sake.
- **Related findings:** PSR-017.

### `CompilationManifest` duplicate structured/flat fields

- **Disposition:** **Narrow**.
- **Current role:** structured decisions plus repeated rendering/compatibility
  fields.
- **Target:** one structured owner; artifact/human-readable fields projected at
  serialization edge.
- **Deletion gate:** audit application/artifact/tool non-test consumers for each
  persisted field before removing compatibility output.
- **Related findings:** PSR-017.

### `CompiledCompatibilityView`

- **Disposition:** **Audit then likely delete/narrow**.
- **Target:** no mirrored in-memory representation solely for convenience.
- **Deletion gate:** current non-test consumer audit confirms required values can
  be projected directly from the authoritative decision.
- **Related findings:** PSR-017.

## Semantic layer

### `ProblemSpec`, semantic component records, `SemanticValidator`

- **Disposition:** **Keep**.
- **Current role:** backend-neutral component graph and semantic/policy
  validation.
- **Target:** unchanged architectural owner; optional syntax cleanup after the
  core numerical refactor.
- **Related findings:** PSR-022.

### `ResolvedProblemView`

- **Disposition:** **Keep**.
- **Current role:** validated indexed borrowing view over one semantic graph.
- **Target:** may gain only demonstrated reusable read-only graph queries.
- **Related findings:** PSR-022.

### `RequirementPolicySpec` optional payload set

- **Disposition:** **Audit then simplify later**.
- **Pressure:** mutually exclusive typed policy selections encoded as many
  independent `std::optional` members.
- **Target option:** actual sum type/variant if migration reduces invalid states
  and validation code.
- **Priority:** below numerical boundary cleanup.
- **Related findings:** PSR-022.

### `reference_specs.hpp` clone-and-mutate builders

- **Disposition:** **Keep then simplify later**.
- **Pressure:** complete graph cloning and stable-string mutation obscure the
  semantic delta.
- **Target option:** named semantic deltas/helpers for demonstrated variations;
  no symbolic PDE DSL.
- **Related findings:** PSR-022, PSR-023.

## Typed numerical targets

### `ScalarDiffusionReactionModel<dim>`

- **Location:** `include/nmopt/dealii/scalar_diffusion_reaction.hpp`
- **Disposition:** **Migrate then delete** unless a final consumer audit proves
  independent production value.
- **Current role:** direct v0 scalar reference implementation containing FE
  state/control, assembly, executable operations, state/adjoint solves,
  metrics, constraints, Hessian, OTD/KKT helpers, and output.
- **Unique behavior to preserve first:** hand-computed weak-form/reference
  oracle, residual/JVP/VJP checks, reduced derivative checks, metric/constraint
  checks, Hessian evidence, supplied-OTD behavior still used by current tests.
- **Migration destination:** surviving compiler/native scalar path and generic
  contracts/components.
- **Deletion gate:** every v0-only non-test consumer is migrated or explicitly
  accepted as a bounded exception; unique oracle coverage runs against the
  surviving path; application output no longer downcasts to this type.
- **Related findings:** PSR-003, PSR-020.

### `ScalarComponentModel<dim>`

- **Disposition:** **Migrate/narrow; final class survival undecided**.
- **Current role:** broad v1 scalar assembled target containing state
  reconstruction, residual assembly, control coupling, observation/loss
  assembly, solves, metrics, constraints, Hessian, diagnostics/output.
- **Valuable internal pieces:** independent state coordinates $P$, $P^T$, fixed
  lifting; component-driven residual assembly; physical/reduced operator
  distinction.
- **Target:** source of reusable typed numerical responsibilities, not a new
  public model hierarchy.
- **Deletion/narrowing gate:** at least two current consumers use extracted
  components and the containing class loses corresponding responsibilities
  with net production deletion.
- **Related findings:** PSR-003, PSR-004, PSR-005, PSR-007, PSR-009.

### `ContinuousControlModel<dim>`

- **Disposition:** **Migrate/narrow; final class survival undecided**.
- **Current role:** complete continuous distributed-control target with state
  and control spaces, coupling, objective, solves, L2/H1/H-1 metrics,
  Hessian, output.
- **Valuable behavior:** independent continuous control coordinates; distinct
  regularization/search metric selections; current LQ reduced-Hessian oracle.
- **Migration candidates:** decision-space operators, affine coupling, objective
  regularization, metric construction, shared state-coordinate/operator
  machinery if evidence supports it.
- **Deletion gate:** no application/compiler consumer requires the complete type
  as one object.
- **Related findings:** PSR-003, PSR-006, PSR-010.

### `NeumannBoundaryControlModel<dim>`

- **Disposition:** **Migrate/narrow; final class survival undecided**.
- **Current role:** state discretization, fixed/mean-zero gauge, transport,
  selected Neumann control realization, observations, solves, objective,
  diagnostics/output.
- **Valuable behavior:** composition with `NeumannControlRealisation`, mean-zero
  solve path, boundary observation maps.
- **Migration candidates:** gauge/solve composition, typed native output,
  observation/tracking components.
- **Deletion gate:** compiler and B2 no longer require the complete model type;
  specialized gauge behavior remains tested.
- **Related findings:** PSR-003, PSR-005, PSR-008, PSR-012.

### `DirichletControlLiftingModel<dim>`

- **Disposition:** **Migrate/narrow; final class survival undecided**.
- **Current role:** independent state coordinates, controlled trace map,
  physical reconstruction, solves, tracking, L2/H1/H1/2 control norms/metrics,
  diagnostics/output.
- **Valuable behavior:** $P$, $L_D$, fixed lifting and dual pullbacks; trace
  metric/regularization evidence.
- **Migration destination:** reusable state-coordinate/reconstruction and
  decision-lifting components if they delete duplication.
- **Deletion gate:** unique chain-rule and trace-norm coverage survives outside
  the complete model.
- **Related findings:** PSR-003, PSR-006, PSR-007, PSR-010.

### `CoefficientIdentificationModel<dim>`

- **Disposition:** **Migrate/narrow; retain specialized parameter-dependent
  machinery unless a smaller composition proves itself**.
- **Current role:** DG0 parameter coordinates, point-dependent diffusion
  operator assembly, parameter JVP/VJP derivative, objective/metric/constraint,
  state/adjoint solves.
- **Valuable behavior:** nonlinear counterexample that prevents overfitting the
  architecture to fixed $A/B/Q/R$ systems.
- **Migration candidates:** parameter coordinate/mass/constraint ingredients;
  point-dependent state-operator service; objective regularization.
- **Deletion gate:** only if replacement retains parameter-dependent residual
  actions and solve behavior without a larger abstraction.
- **Related findings:** PSR-001, PSR-006, PSR-008.

## Focused numerical components

### `NeumannControlRealisation<dim>`

- **Disposition:** **Narrow / keep**.
- **Good core:** demonstrated polymorphism between facewise and continuous
  trace realizations; layout/coordinates; coupling and transpose actions.
- **Pressure:** also owns regularization, L2 metric creation, and native output;
  facewise box construction requires concrete recovery.
- **Target:** preserve only responsibilities shared by actual realizations;
  expose reusable boundary mass/geometry data to regularization, metric,
  constraint, and output consumers where this removes code.
- **Decision gate:** narrowing must delete existing methods/caller branching,
  not merely add indirection.
- **Related findings:** PSR-005, PSR-011.

### `VolumeObservationAssembly<dim>`

- **Disposition:** **Keep / reuse**.
- **Current role:** focused assembly of volume/subdomain quadratic tracking data.
- **Target:** reuse consistently where it replaces embedded equivalent assembly.
- **Naming note:** review whether the surviving responsibility is better named
  as a tracking-form assembly only if the API changes materially.
- **Related findings:** PSR-009.

### `MassMetric`

- **Disposition:** **Keep / minor deduplication**.
- **Pressure:** local inline CG implementation differs from shared SPD helpers.
- **Target:** use shared numerical solve helper if it genuinely reduces code and
  preserves behavior; no public solver hierarchy.
- **Related findings:** PSR-010.

### `Hminus1Metric`, `TraceHhalfMetric`

- **Disposition:** **Keep**.
- **Current role:** mathematically intrinsic metric operators and inverses.
- **Possible later refinement:** neutral reusable trace quotient operator if
  objective regularization and metric consumers demonstrate enough duplicated
  semantics/code.
- **Do not:** extract a neutral operator solely because the current name is
  imperfect.
- **Related findings:** PSR-010.

## KKT and second-order paths

### generic `EqualityConstrainedQuadraticKKTProductT<Backend>`

- **Disposition:** **Keep**.
- **Target:** preferred backend-neutral KKT product.
- **Related findings:** PSR-021.

### `ScalarDiffusionReactionKKT<dim>`

- **Disposition:** **Audit then likely delete**.
- **Current role:** scalar-specific deal.II KKT realization tied to the direct
  scalar path.
- **Deletion gate:** no current non-test consumer requires behavior not
  expressible through the generic KKT product + serial adapter; unique KKT
  oracle tests migrate first.
- **Related findings:** PSR-020, PSR-021.

### complete-model `ReducedHessianT` inheritance

- **Disposition:** **Migrate/narrow**.
- **Current role:** some complete LQ targets directly implement optional reduced
  Hessian capability.
- **Target:** compose optional Hessian provider from compatible typed numerical
  ingredients where shared implementation is demonstrated.
- **Do not:** require nonlinear/Dirichlet/parameter targets to implement a fixed
  LQ formula.
- **Related findings:** PSR-010.

## Application layer

### B1/B2 `dynamic_cast` recovery from `ExecutableModelT`

- **Disposition:** **Delete**.
- **Current role:** recover dimensions, objective components, native output, and
  target-specific diagnostics after compiler type erasure.
- **Migration destination:** explicitly retained typed/native realization or
  small application-owned closures/results created before erasure.
- **Deletion gate:** all native output/evidence consumers have typed ownership;
  no new methods are added to `ExecutableModelT` to compensate.
- **Related findings:** PSR-012, PSR-018, PSR-019.

### model-owned `write_native_output` / filesystem methods

- **Disposition:** **Move outward then delete from numerical model APIs**.
- **Target:** application/native output code with access to typed FE state.
- **Deletion gate:** filenames/topology/field identity remain regression-tested;
  application path no longer downcasts.
- **Related findings:** PSR-018.

### B1/B2 execution adapters

- **Disposition:** **Keep / simplify only after boundary repair**.
- **Current role:** legitimate application composition root plus verification,
  optimization, diagnostics, output, and experiment packaging.
- **Target:** first remove numerical ownership violations. Split/share further
  only when current duplication is demonstrated and the resulting unit is
  easier to review.
- **Related findings:** PSR-019.

### `ProblemRecipeT`, `ScenarioT`, `HeadlessBenchmarkRunnerT`

- **Disposition:** **Keep**.
- **Target:** no refactor unless a core-caused binding change requires it.
- **Related findings:** general audit good boundary.

## Output and persisted evidence

### native output files and writers

- **Disposition:** **Move outward / keep behavior**.
- **Target:** application-owned, typed backend output.
- **Compatibility:** current application filenames and field identities remain
  stable unless a separate accepted output-contract change says otherwise.
- **Related findings:** PSR-018.

### artifact manifest compatibility fields

- **Disposition:** **Audit then project at edge**.
- **Target:** preserve persisted compatibility only for active readers; remove
  mirrored in-memory compiler state.
- **Related findings:** PSR-017.

## Documentation

### `docs/design/pde-solver-boundary.md`

- **Disposition:** **Keep as new authority**.
- **Role:** long-lived ownership and solver-boundary decision.

### this review package

- **Disposition:** **Keep**.
- **Role:** assessment evidence, deletion ledger, and bounded execution roadmap.

### existing `system-blueprint.md` and `composition-boundaries.md`

- **Disposition:** **Keep; update only when implementation truth changes or a
  contradiction appears**.
- **Reason:** current component/port direction is broadly compatible with the
  audit.

### old implementation/application roadmaps and references

- **Disposition:** **Do not mass-edit during R0**.
- **Target:** each implementation unit updates the authoritative docs it
  invalidates; final program performs a stale-reference sweep.
- **Related findings:** PSR-024.

## Deferred cleanup candidates

These are intentionally outside the first numerical-boundary sequence unless a
current unit touches them directly:

- convert semantic policy payload optionals to a sum type;
- introduce selective tagged semantic IDs;
- replace clone-and-mutate reference builders with named semantic deltas;
- physically split large semantic validation/reference files after
  responsibilities are cleaner;
- reconsider `CompilationResultT` nullable product representation;
- neutralize naming of shared $H^{1/2}$ numerical operator if reuse proves it;
- introduce any shared quadratic-form helper only after duplicate code count is
  measured; and
- outer runner/post-processing simplification unrelated to a changed core
  contract.

## Final deletion accounting

At the end of the program, record separately:

```text
production C++ added / deleted
production Python added / deleted
tests added / deleted
documentation added / deleted
complete production types deleted
runtime dispatch branches deleted
duplicate stored representations deleted
remaining concrete-model downcasts
remaining model-owned filesystem operations
remaining v0/version-specific production paths and their reasons
```

Net-negative production code is a design-pressure target for the completed
cleanup, not permission to sacrifice clarity or mathematical verification in
one intermediate unit.
