# deal.II compiler implementation

## Purpose and boundary

This document describes how the current semantic/deal.II compiler is implemented.

It is not the public authoring guide and it is not a second capability catalogue.
For supported authoring and compilation workflows, use
[Problem authoring](../reference/problem-authoring.md) and
[Compiler reference](../reference/compiler.md). This page explains the internal
stages, ownership, provenance construction, product packaging, and test boundaries
behind those interfaces.

The compiler path is one producer of numerical/formulation contracts:

```text
ProblemSpec
    │
    ▼
SemanticResolver
    │
    ▼
ResolvedProblemView
    │
    ▼
closed compiler request
    │
    ├──────────────────────────┐
    ▼                          ▼
bounded scalar planner    closed target strategy
    │                          │
    ▼                          ▼
ScalarLoweringPlan        target-specific model
    │                          │
    ▼                          │
ScalarComponentModel           │
    │                          │
    └────────────┬─────────────┘
                 ▼
       deal.II numerical services
                 │
                 ▼
         typed formulation product
```

An application that already owns its PDE realization may bypass this entire path and
produce the downstream numerical contracts directly.

## Compiler entry point and ownership modes

`include/nmopt/compiler/v1/dealii_compiler.hpp` exposes `DealiiCompiler`.

Its constructor accepts:

- `DealiiCapabilityRegistryV1`;
- `DealiiScalarLowererRegistryV1`.

The main operations are:

```cpp
ValidationReport validate(
  const ProblemSpec &,
  const DealiiDiscretisationPolicy &,
  CompilationProduct = CompilationProduct::reduced_dto) const;
```

and two `compile()` families:

```text
caller-owned triangulation
    compile(specification, triangulation, bindings, ...)

owned compilation session
    compile(specification, shared_ptr<DealiiCompilationSession<dim>>, bindings, ...)
```

The triangulation-reference overload is a borrowed lifetime path. The session overload
owns the mesh lifetime and is the preferred choice when compiled services may outlive
the compilation call.

The compiler currently targets `dealii_backend::SerialBackend`.

## Semantic resolution comes first

`SemanticValidator` checks the graph structure and declared semantic policies.
`SemanticResolver` then creates a `ResolvedProblemView` only when that report is
valid.

A `ResolvedProblemView`:

- borrows the source `ProblemSpec`;
- indexes regions, spaces, pairings, variables, data, transformations, residual
  terms, equations, observations, losses, metrics, constraints, and requirements by
  stable ID;
- removes repeated string-search logic from later compiler stages.

The view is intentionally scoped to one resolution/compilation operation.

Semantic validity is not deal.II lowerability. A semantically coherent graph can
still fail because the compiler has no registered realization for a selected
combination, policy, formulation product, or runtime binding.

## Validation order

`DealiiCompiler::validate()` and the first part of `compile_impl()` follow the same
ordering:

```cpp
auto resolution = SemanticResolver().resolve(specification);
if (!resolution.diagnostics.valid())
  return;

auto request = resolve_compilation_request(*resolution.problem);
auto registration =
  resolve_dirichlet_control_registration(*resolution.problem, request);

close_compilation_request(*resolution.problem, request, registration);

validate_lowerability(specification, request, policy, report);
validate_formulation_capability(specification, report);
validate_supplied_otd_capability(specification, request, product, report);
validate_dirichlet_control_registration(*resolution.problem, request, report);
validate_product_capability(specification, request, policy, product, report);
```

This ordering matters.

Semantic errors are reported before backend capability checks. The compiler then
closes its own request before checking lowerability so later phases reason from one
resolved branch choice rather than repeatedly rediscovering a cross-product of
semantic predicates.

## The closed compiler request

`ResolvedCompilationRequest` in `dealii_types.hpp` is the compiler's internal
decision record before numerical construction.

It owns compiler-layer facts such as:

- `ResolvedTargetFamily`;
- `ResolvedDirichletRegistration`;
- whether fixed reconstruction, Dirichlet control, Neumann control, continuous
  control, coefficient identification, general scalar assembly, special observations,
  or special metrics are required;
- selected boundary/trace/transposition/metric policies;
- resolved semantic region and data IDs;
- expected concrete data-binding ports.

These are **not** new semantic component kinds. They are a closed interpretation of a
validated graph for this compiler.

`ResolvedDataBindingRequest` is particularly important. It records:

```text
semantic datum
    │
    ├── role and kind
    ├── semantic space / region
    ├── compiler binding port
    ├── evaluation realization
    └── expected runtime representation
```

That prevents later binding validation and manifest construction from guessing which
runtime datum corresponds to which semantic node.

## Runtime data are not stored in `ProblemSpec`

`DealiiDataBindings<dim>` supplies concrete values only after semantic validation.

The core bundle contains:

- forcing and desired-state `dealii::Function<dim>` references;
- optional constant diffusion;
- reaction and regularization values;
- optional fixed-Dirichlet data;
- provenance strings.

Optional sub-bundles cover:

- general scalar tensor/vector/scalar coefficients;
- weighted trace data;
- conservative transport data for the narrow Neumann-convection route;
- additive natural-boundary source data.

Cellwise and facewise box bounds are separate binding types because their coefficient
topologies differ.

This separation keeps the semantic graph deal.II-independent while still allowing the
compiler to validate concrete shape, provenance, and evaluation requirements.

## Discretisation policy

`DealiiDiscretisationPolicy` owns backend choices that are not semantic problem
meaning, including:

- state degree;
- assembled versus matrix-free execution selection;
- optional transport-boundary realization;
- optional volume-observation realization;
- state and adjoint solve policies;
- control-metric solve policy;
- PDAS and inner KKT policy.

The existence of an enum value does not imply lowerability. For example, the current
compiler supports assembled execution; unsupported selections are rejected during
compiler validation.

## Two realization mechanisms

The compiler does not lower every supported problem through one mechanism.

### Bounded scalar component planning

`dealii_scalar_plan.hpp` implements a registry-driven component planner for a bounded
scalar subset.

`DealiiScalarLowererRegistryV1` has handlers for registered:

- residual-term kinds;
- observation kinds;
- loss kinds;
- metrics;
- constraints;
- transformations.

Each handler contributes typed entries and a stable handler ID to
`ScalarLoweringPlan`.

The resulting plan records:

```text
semantic problem/state/control/equation IDs
residual contributions
observation contributions
loss contributions
metric / constraint / transformation choices
typed data placements
boundary selections
tracking/observation geometry
handler provenance
```

Two narrower projections are derived from it:

- `ScalarResidualAssemblyPlan` – residual contributions and typed data placements
  needed by equation assembly;
- `ScalarServicePlan` – observation/loss/metric/constraint/transformation information
  used to construct executable services.

The planner is deliberately bounded. Missing handlers produce lowerability
diagnostics; a specialized semantic graph is not silently forced through the generic
scalar plan.

### Closed target-specific strategies

Other capabilities remain clearer as closed target strategies.

The main implementation families include:

- Neumann/boundary control in `dealii_neumann_boundary.hpp`;
- Dirichlet control/lifting in `dealii_dirichlet_control.hpp`;
- continuous control in `dealii_continuous_control.hpp`;
- coefficient identification in `dealii_coefficient_identification.hpp`;
- fixed-Dirichlet/state-coordinate support in `dealii_fixed_dirichlet.hpp`;
- specialized observation/realization support in the corresponding compiler headers.

The compiler selects one target family in the closed request, validates its typed
policies, then constructs that family directly. In particular, Dirichlet-control
realizations are selected by closed structural signatures: independently supported
tracking, loss, metric, or boundary choices do not imply an automatic cross-product or
nearest-match fallback.

This hybrid is intentional:

```text
resolved request
      │
      ├──────────────────────────┐
      ▼                          ▼
bounded compositional      closed specialized
scalar subset              target
      │                          │
      ▼                          ▼
ScalarLoweringPlan         target-specific model
      │                          │
      ▼                          │
ScalarComponentModel            │
      │                          │
      └────────────┬─────────────┘
                   ▼
          common numerical contracts
```

A future feature should join the component planner only when it has a genuinely
independent contribution contract. It should not be forced into that registry merely
to make the implementation look uniform.

## Main compilation flow

After the common validation sequence succeeds, `compile_impl()` performs the
following stages.

### 1. Resolve target-specific preconditions

The compiler derives selected flags and checks mesh/reference-cell, typed boundary,
transposition, fractional metric, target-data, box-data, and related realization
requirements.

Failures are returned as typed diagnostics rather than allowing construction to
proceed with a guessed policy.

### 2. Resolve concrete runtime bindings

The closed data-binding requests are checked against `DealiiDataBindings<dim>` and
optional box bindings.

The compiler validates concrete shape and required provenance at this boundary.

### 3. Build the bounded scalar plan when applicable

For scalar component targets, `DealiiScalarLoweringPlanner` turns the resolved
semantic graph into the plan described above.

The plan is both executable input and provenance: the exact handler IDs used for each
component are later recorded rather than inferred from display strings.

### 4. Create the typed compilation decision

Before numerical model construction is finalized, the compiler creates a
`ResolvedCompilationDecision`.

This is the typed source of truth for compilation provenance. It records the closed
semantic/compiler decision, including:

- target and formulation;
- mesh policy and lifetime mode;
- regions, spaces, pairings, and bindings;
- residual, observation, loss, and transformation realizations;
- selected semantic assumptions/policies;
- realized-map skeletons;
- formulation-specific records.

Model construction consumes the same closed request and scalar plan used to build
this decision.

### 5. Construct numerical services

Depending on the selected target, the compiler constructs a concrete numerical model
and the services needed by the selected formulation.

The reduced path may construct:

- an `ExecutableModelT<SerialBackend>`;
- `MetricT<SerialBackend>`;
- optional `ConstraintT<SerialBackend>`;
- `StateAdjointSolversT<SerialBackend>`;
- optional `ReducedHessianT<SerialBackend>`;
- optional shared compiled box data;
- optional `CompiledApplicationViewT<SerialBackend>`.

The supplied-OTD path additionally constructs a `SuppliedOTDSystemT`.

### 6. Finalize the compilation decision

`finalize_resolved_decision()` projects realized dimensions, solve policies,
operators, maps, metric/constraint records, and other constructed facts into the
typed decision.

This finalization is deliberately not a second target-selection pass. It also reads
the `ProblemSpec` inventory to render semantic identifiers and records, but the target
has already been selected in `ResolvedCompilationDecision` and
`ResolvedCompilationRequest`. It does not use those raw graph inventories to choose a
second target.

### 7. Build optional secondary formulation products

For `CompilationProduct::quadratic_kkt` or `CompilationProduct::pdas`, the compiler
constructs the corresponding KKT product after the numerical services exist.

A supplied-OTD system can be bridged to the canonical quadratic KKT product only when
its typed validity declaration establishes the required affine/block/pairing/sign
properties.

### 8. Package exactly the selected product family

`CompilationResultT<Backend>` contains diagnostics plus separate slots for:

```text
CompiledProblemT
CompiledQuadraticKKTProblemT
CompiledPDASProblemT
CompiledSuppliedOTDProblemT
```

`CompilationResultT::succeeded()` requires valid diagnostics and at least one product.

These slots are intentionally distinct. `CompiledProblemT` is not a universal
container for every formulation. `compiled_problem.hpp` remains a compatibility
aggregate that includes the manifest and compiled-product headers; it does not define
a second product hierarchy.

## Reduced compiled product

`CompiledProblemT<Backend>` owns the compiler-path reduced service bundle:

- erased executable model;
- metric;
- optional constraint;
- state/adjoint solve services;
- compilation manifest;
- lifetime owner;
- optional reduced Hessian;
- optional shared box data;
- optional compiled application view.

`make_reduced_dto()` constructs `ReducedDTOT<Backend>` from the executable model,
state/control partition, solve services, and lifetime owner.

The reduced DTO does not need `CompiledApplicationViewT`.

## KKT, PDAS, and supplied OTD products

`CompiledQuadraticKKTProblemT` owns a quadratic KKT product, manifest, and optional
lifetime owner.

`CompiledPDASProblemT` additionally owns:

- complementarity;
- metric;
- shared compiled box data;
- projection constraint;
- PDAS policy;
- inner KKT solver policy.

The constructor checks that the metric, complementarity, projection constraint, and
shared box token describe one compatible realization.

`CompiledSuppliedOTDProblemT` owns the supplied OTD system, manifest, and lifetime
owner. It intentionally does not expose the reduced DTO service bundle.

This product separation is the implementation reason documentation must not describe
one universal compiled formulation interface.

## `CompiledApplicationViewT`

`CompiledApplicationViewT` is an optional compiler-produced sidecar for application
needs that should not be added to `ExecutableModelT`.

It exposes:

- physical/independent state and control dimensions;
- realized observation dimension;
- native-output callback;
- optional objective-component callback.

Current compiler construction attaches a view to:

- the direct-volume scalar target;
- continuous-control targets built by `ContinuousControlModel`;
- Neumann/boundary-control targets.

The current Dirichlet-control and coefficient-identification branches do not attach
one. General scalar component targets only attach the volume view when the resolved
target family is `direct_volume`.

The Neumann view supplies objective components; the generic volume view supplies
dimensions and native output but no objective-component callback.

KKT, PDAS, and supplied-OTD wrapper products do not expose this reduced-product
sidecar.

### Lifetime caveat

The view retains the numerical model captured by its callbacks. Native-output
callbacks may also borrow compile-time data bindings such as forcing and desired-state
`Function` objects.

Keeping an owned compilation session alive through the compiled product therefore
protects the mesh/model lifetime, but callers must still keep any borrowed data
objects alive while invoking the view.

## Compilation sessions and detached services

`DealiiCompilationSession<dim>` exclusively owns a triangulation moved into it. The
compiler alone gets mutable access during lowering.

The session overload records mesh lifetime as `owned_session` and passes the session
as the lifetime token retained by compiled products.

This enables patterns such as:

```text
owned session
    │
    ▼
compiled problem
    │
    ├── make_reduced_dto()
    │          │
    │          └── detached reduced service retains lifetime token
    │
    └── manifest copied independently
```

The compiler tests explicitly destroy the local compilation/session handles and then
exercise the detached reduced and supplied-OTD/KKT services.

The triangulation-reference overload remains available for source compatibility. Its
manifest records a borrowed immutable mesh lifetime; callers own the obligation to
keep that mesh alive.

## Manifest and provenance

Every successful compiler product carries `CompilationManifest`, currently schema
version 4.

Its primary record is `ResolvedCompilationDecision`. The manifest includes typed
records for:

- mesh identity, provenance, and lifetime;
- semantic regions and spaces;
- concrete bindings;
- pairings;
- residual/observation/loss/transformation realizations;
- realized spaces and maps;
- state and adjoint solve policies;
- metric and constraint realization;
- supplied OTD, KKT, and PDAS products when present;
- typed model-author assumptions and selected policies.

`CompiledCompatibilityView` is a rendered artifact-facing projection. It is not a
second executable configuration source.

The manifest separates:

```text
typed compilation decision
          │
          ▼
descriptive compatibility rendering
```

rather than parsing descriptive strings to recover executable choices.

Binding records can include concrete field shape, runtime representation, checked
status, scalar values, and value digests. Mesh records keep caller provenance
separate from a structural identity derived from the realized mesh.

## Diagnostics versus contract failures

Caller-correctable semantic, backend, binding, and product failures should normally be
returned as `Diagnostic` entries. Examples include unsupported semantic composition,
unsupported execution/product selection, missing or incompatible runtime data, empty
or unsupported meshes, missing bounds, and unsupported KKT/PDAS combinations.

`ContractError` remains the lower-level invariant boundary for direct misuse of
contract constructors and states that should be unreachable after successful compiler
validation. Maintainers should preserve this distinction: turning an ordinary
unsupported request into a contract failure makes the compiler boundary less
inspectable, while silently substituting another realization would change mathematical
meaning.

## Why the compiler composition root remains one large header

`dealii_compiler.hpp` remains a large composition root. It has readable regions for:

- request construction;
- binding validation;
- orchestration;
- lowerability/capability checks;
- semantic lookup;
- numerical realization;
- metadata/provenance;
- manifest projection;
- KKT/PDAS construction.

Those regions are not yet independent compiler phases. They share request state and
construction details heavily enough that a `detail/*.hpp` split would mostly move
navigation boundaries rather than create typed architecture.

The current stopping point is therefore coarse source-navigation comments plus
stronger typed seams around validation, scalar planning, manifest records, compiled
products, and target-specific models. A future physical split should follow a real
typed phase boundary rather than file size alone.

The same reasoning applies to `semantic/v1/detail/validation_detail.hpp`: the public
facade is already narrow, while the private validation phases share indexing, ordering,
and policy context closely enough that physical fragmentation would not create a
stronger interface.

## Verification structure

Compiler evidence is split by responsibility.

Backend-neutral planning is covered by:

```text
tests/compiler/scenarios/scalar_lowering_plan.hpp
tests/semantic/semantic_v1_contract.cc
```

The deal.II capability suite uses:

```text
tests/compiler/dealii/dealii_compiler_contract.cc
tests/compiler/dealii/dealii_compiler_support.hpp
tests/compiler/dealii/scenarios/
```

One heavy compiler translation unit is intentional. Scenario bodies are grouped by
capability so source navigation improves without repeatedly compiling the
header-heavy compiler/deal.II template surface.

The scenario groups cover, among other responsibilities:

- compiled formulation products;
- control metrics;
- diagnostics;
- Dirichlet lowering;
- natural-boundary/Neumann lowering;
- observations;
- session ownership.

Application-level evidence for the compiler-produced sidecar lives in
`tests/application/compiled_application_view_dealii_contract.cc`.

Lower-level deal.II service behavior belongs under `tests/dealii/`, not in the
compiler capability suite.

The verification pattern deliberately combines independent value checks with
derivative and formulation identities: focused scenarios use finite differences for
JVPs, JVP/VJP pairing checks, and state-recomputed reduced-derivative checks where
applicable. Exact scenario counts and assertions belong in the tests rather than in a
second documentation ledger.

## Where to change the compiler

For a new semantic component kind, begin in the semantic types/validation layer. Do
not add a compiler branch that depends on an undeclared semantic convention.

For a new bounded scalar contribution:

```text
semantic kind
    │
    ▼
DealiiScalarLowererRegistryV1 handler
    │
    ▼
ScalarLoweringPlan contribution
    │
    ▼
ScalarResidualAssemblyPlan or ScalarServicePlan
    │
    ▼
ScalarComponentModel realization
    │
    ▼
compiler scenario + manifest provenance
```

For a specialized target, prefer an explicit target realization when the feature owns
coupled spaces, lifting, gauge, trace, or parameterization behavior that is not an
independent scalar contribution.

For a new formulation product, add the backend-neutral contract first, then compiler
capability validation and a distinct compiled product when its required actions differ
from existing products.

For application-only output or diagnostics, do not expand `ExecutableModelT` merely
because one compiled target needs more information. Extend the optional compiled
application seam only when the capability is genuinely compiler-produced and can
respect the required lifetimes.
