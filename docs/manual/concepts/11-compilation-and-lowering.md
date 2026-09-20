# Compilation and lowering

## From an accepted semantic request to an executable numerical product

The previous two chapters stopped at an accepted request.

We had a mathematical problem, translated it into a `ProblemSpec`, validated the
graph, resolved its stable IDs, closed the compiler request, and checked that the
requested realization and product were supported.

Nothing numerical had been assembled yet.

This chapter follows the next transition:

```text
accepted semantic problem
        │
        ▼
compiler planning
        │
        ▼
typed lowering plan / registered target
        │
        ▼
runtime data + mesh/discretization policy
        │
        ▼
finite-element model
        │
        ▼
generic executable services
        │
        ▼
CompiledProblemT
```

The running example is again the scalar distributed-control problem

```math
\begin{aligned}
-\nabla\cdot(\kappa\nabla y)+cy &= f+u
&& \text{in }\Omega,\\
y &= 0
&& \text{on }\partial\Omega,
\end{aligned}
```

with objective

```math
J(y,u)
=
\frac12\lVert y-y_{\mathrm d}\rVert_{L^{2}(\Omega)}^{2}
+
\frac{\beta}{2}\lVert u\rVert_{L^{2}(\Omega)}^{2}.
```

For this baseline problem, the semantic graph already says that:

```text
state equation
    diffusion_reaction
    volume_source
    volume_control

objective
    full-volume state tracking
    full-volume control regularisation

geometry
    control L2 metric

formulation
    reduced DTO
```

Compilation must turn those descriptions into actual spaces, operators, solves,
metrics, and formulation services.

The essential question is no longer

> What does the problem mean?

It is now

> Which executable object realizes each declared semantic role, and how are those
> objects wired together?

## 1. Lowering changes representation, not mathematical intent

The word *lowering* is borrowed from compiler terminology.

Here it means replacing a higher-level semantic declaration by a more concrete
numerical representation while preserving the intended mathematical role.

For example,

```text
ResidualTermKind::diffusion_reaction
```

may eventually become the assembled action

```math
\underbrace{
\int_{\Omega}
\kappa\nabla y\cdot\nabla v
+
cyv
\mathrm{d}x
}_{\text{concrete FE bilinear contribution}}.
```

Likewise,

```text
MetricKind::l2
```

may become a mass-matrix Riesz map

$$
G\mathbf u=M\mathbf u.
$$

The important direction is

```text
semantic meaning
    ↓
selected numerical realization
```

and not the reverse.

A lowerer does not inspect an arbitrary matrix and guess that it "looks like" an
$L^{2}$ metric. The semantic graph has already declared the metric; the compiler
constructs the registered realization of that declaration.

## 2. The public compiler call joins four kinds of input

At the application boundary, `DealiiCompiler::compile()` receives four conceptually
different inputs:

```text
ProblemSpec
    what problem is requested

Triangulation or DealiiCompilationSession
    which concrete mesh is used and who owns it

DealiiDataBindings
    concrete forcing / target / coefficient values

DealiiDiscretisationPolicy
    registered discretization and solve choices
```

Optional bound data and the requested `CompilationProduct` complete the call.

Schematically:

```cpp
auto result =
  compiler.compile(specification,
                   triangulation,
                   data,
                   policy,
                   bounds,
                   facewise_bounds,
                   CompilationProduct::reduced_dto);
```

The arguments should not be collapsed into one configuration object because they
belong to different layers.

`ProblemSpec` is semantic and deal.II-free.

The mesh and data bindings are concrete numerical inputs.

The discretization policy selects implementation choices within the registered
compiler capability.

The compilation product says which solver-facing surface should be packaged at the
end.

## 3. Runtime data bind semantic ports to concrete values

Chapter 9 introduced declarations such as

```text
DataSpec
    id   = forcing
    role = forcing
    kind = function
```

That declaration does not contain the actual function $f(x)$.

Compilation receives the concrete value through `DealiiDataBindings`.

For the baseline problem the main runtime bindings are conceptually

```text
semantic port             concrete runtime object
---------------------------------------------------------
forcing                    dealii::Function<dim>
desired_state              dealii::Function<dim>
diffusion                  scalar coefficient
reaction                   scalar coefficient
regularisation_weight      scalar beta
```

The binding also carries provenance labels.

This produces a useful separation:

```text
ProblemSpec
    "there exists forcing data called forcing"

DealiiDataBindings
    "for this compilation, forcing is this concrete Function object"
```

A parameter study can therefore keep the semantic graph fixed while changing the
concrete forcing or target.

### 3.1 Binding shape is part of lowering

Some problems require more structured runtime data.

Examples include:

```text
tensor diffusion
    dealii::TensorFunction<2, dim>

transport field
    dealii::TensorFunction<1, dim>

weighted boundary observation
    separate scalar Function weight

fixed Dirichlet reconstruction
    separate boundary Function

natural boundary source
    separate boundary Function
```

The compiler does not reinterpret a generic vector of numbers as whichever datum is
currently missing.

The resolved request identifies the semantic port, expected runtime representation,
evaluation location, and whether the binding is required.

That information later appears in compilation diagnostics and provenance records.

## 4. `ResolvedDataBindingRequest` bridges semantics and runtime data

Before model construction, the compiler derives explicit binding requests.

A request records information of the form

```text
semantic_id
    forcing

role
    DataRole::forcing

kind
    DataKind::function

space_id / region_id
    where the datum belongs semantically

evaluation_realisation
    volume_quadrature

runtime_representation
    dealii::Function<dim>

required
    true
```

This is more precise than carrying only the string `"forcing"`.

The same semantic `DataKind::function` can be consumed in different ways:

```text
volume source
    evaluate at volume quadrature

boundary weight
    evaluate at face quadrature

fixed Dirichlet datum
    interpolate on boundary DoFs
```

The compiler resolves those consumption rules before numerical construction.

## 5. The compiler has two lowering styles

The current deal.II compiler deliberately uses two styles side by side.

### 5.1 Bounded component planning

The baseline scalar volume problem and several related scalar targets use a
compositional planner.

Their semantic components are translated independently into a typed
`ScalarLoweringPlan`.

This path is used for the direct scalar volume family and for several assembled scalar
variants such as fixed reconstruction, subdomain/H1 observations, point sensors,
normal flux, and the general scalar Robin composition.

### 5.2 Registered target strategies

Other problem families retain dedicated strategies.

Examples include:

```text
continuous volume control
Neumann boundary control
Dirichlet control
coefficient identification
specialized trace-metric targets
```

These targets still consume the same validated semantic vocabulary and closed compiler
request, but their numerical realization is specialized enough that one dedicated
model strategy is clearer than forcing every detail through the scalar component
planner.

The compiler therefore looks like

```text
ResolvedCompilationRequest
        │
        ├── scalar component target
        │       ↓
        │   ScalarLoweringPlan
        │
        │
        └── specialized target
                ↓
            registered target strategy
```

Both paths converge again at the shared executable contracts and
`CompiledProblemT`.

## 6. `ScalarLoweringPlan` – the worked compiler IR

For the scalar component path, `ScalarLoweringPlan` is the central intermediate
representation.

It is neither the semantic graph nor the final finite-element model.

It records the compiler's selected operators and services.

Its main fields can be grouped as

```text
identity
    semantic_problem_id
    state_variable_id
    decision_variable_id
    equation_id

residual
    residual_terms[]

objective
    observations[]
    losses[]

optimization services
    metric
    constraint
    transformation

realization selections
    boundary / transposition / partial-boundary /
    fractional-metric / boundary-H1 selections

data placement
    data_placements[]

geometric selections
    Dirichlet / Robin / tracking / flux / sensor regions

provenance
    semantic component <- lowering handler
```

That makes the plan a useful compiler IR: it is concrete enough to drive numerical
construction, but still independent of the final matrices and vectors.

## 7. One semantic residual term becomes one typed compiler contribution

Consider the semantic term

```text
volume_control
    kind        = volume_control
    variable    = control
    equation    = state_equation
```

The scalar registry maps its semantic kind to

```text
ScalarResidualOperatorKind::volume_control
```

and assigns a stable handler ID.

The plan contribution retains

```text
component_id
operator_kind
variable_ids
data_ids
region_id
handler_id
```

so the translation remains inspectable.

Conceptually:

```text
ResidualTermSpec
    "volume_control"
        │
        ▼
ScalarResidualTermHandler
        │
        ▼
ScalarResidualContribution
    operator_kind = volume_control
    variable_ids  = {control}
    handler_id    = dealii.scalar.residual.volume_control
```

No matrix has been assembled at this point.

The compiler has only decided which registered numerical operator will realize the
semantic term.

## 8. Observations and losses are lowered separately

The objective path follows the same pattern.

For full-domain state tracking,

```math
\underbrace{
\mathcal O_{y}(y)=y
}_{\texttt{volume\_restriction}}
```

becomes a `ScalarObservationContribution` with

```text
operator_kind = volume_restriction
input_variable_id = state
output_space_id = state_observation_space
```

The tracking loss

```math
\frac12
\lVert
\mathcal O_{y}(y)-y_{\mathrm d}
\rVert^{2}
```

becomes a separate `ScalarLossContribution`.

This separation matters because the observation can change without changing the
mathematical form of the loss.

For example,

```math
\mathcal O_{y}(y)=y
```

can be replaced by

```math
\mathcal O_{\mathrm s}(y)
=
[y(x_{1}),\ldots,y(x_{m})]^{\mathsf T}
```

while the outer objective remains a quadratic tracking loss.

The plan records those two decisions independently:

```text
ObservationSpec
    point_sensor
        │
        ▼
ScalarObservationContribution
    point_sensor

LossSpec
    quadratic_tracking
        │
        ▼
ScalarLossContribution
    quadratic_tracking
```

This is the compiler counterpart of the observation/loss distinction introduced in
Chapter 9.

## 9. Data placement records where a datum is actually consumed

A semantic datum such as

```text
desired_state
```

has a role and a semantic space.

The lowering plan adds another fact:

```text
where/how will this datum be evaluated by the numerical realization?
```

`ScalarDataPlacement` records

```text
semantic_id
role
kind
space_id
region_id
evaluation
handler_id
```

Current scalar evaluations distinguish at least

```text
volume_quadrature
boundary_face_quadrature
```

The baseline forcing and desired state are volume data.

A Robin source belongs to a boundary region and is evaluated at face quadrature.

The plan therefore preserves the semantic identity while adding compiler-specific
placement:

```text
DataSpec
    id = robin_source
    role = robin_source
        │
        ▼
ScalarDataPlacement
    region = robin_boundary
    evaluation = boundary_face_quadrature
    handler = ...
```

The actual function value is still supplied separately by the runtime binding.

## 10. The handler registry is a translation table, not a numerical model

`DealiiScalarLowererRegistryV1` stores handlers for supported semantic kinds.

Conceptually, its mapping looks like:

| Semantic kind | Scalar operator kind | Role of handler |
| --- | --- | --- |
| `diffusion_reaction` | `diffusion_reaction` | add diffusion/reaction residual contribution |
| `volume_source` | `volume_source` | add source contribution |
| `volume_control` | `volume_control` | add distributed-control contribution |
| `volume_restriction` | `volume_restriction` | add full/subdomain observation |
| `quadratic_tracking` | `quadratic_tracking` | add state-tracking objective contribution |
| `quadratic_control_regularisation` | `quadratic_control_regularisation` | add control penalty |
| `l2` metric | `cellwise_l2` | select cellwise mass metric realization |
| `cellwise_box` | `cellwise_box` | select coefficientwise box service |
| `fixed_dirichlet_reconstruction` | `fixed_dirichlet_reconstruction` | select physical-state reconstruction |

The handler does not own a `SparseMatrix`.

Its job is to translate a validated semantic component into a typed plan record.

For a residual term, the implementation is essentially:

```text
check semantic kind
    ↓
append ScalarResidualContribution
    ↓
record component_id <- handler_id provenance
```

The numerical model consumes the resulting plan later.

## 11. Provenance is recorded during lowering, not reconstructed afterward

Each scalar handler appends a record of the form

```text
semantic component
    <- lowering handler
```

For the baseline residual, the plan can retain facts such as

```text
diffusion_reaction
    <- dealii.scalar.residual.diffusion_reaction

volume_source
    <- dealii.scalar.residual.volume_source

volume_control
    <- dealii.scalar.residual.volume_control
```

This matters because once the compiler has produced a matrix, several semantic
contributions may be mixed inside it.

Reconstructing provenance from the matrix afterward would be unreliable.

The lowering plan records the decision at the moment it is made.

## 12. The baseline plan is easy to read as a compiler object

For the running problem, the essential `ScalarLoweringPlan` is conceptually

```text
problem
    scalar_diffusion_reaction

state
    state

decision
    control

equation
    state_equation

residual_terms
    diffusion_reaction
    volume_source
    volume_control

observations
    state volume_restriction
    control volume_restriction

losses
    quadratic_tracking
    quadratic_control_regularisation

metric
    cellwise_l2

constraint
    none
    or cellwise_box

transformation
    none

data placements
    forcing        -> volume quadrature
    desired_state  -> volume quadrature
```

The mathematical problem has not changed.

What changed is that each semantic role now points to a registered numerical
operator/service kind.

This is the compiler's last convenient representation before finite-element
construction begins.

## 13. `ScalarResidualAssemblyPlan` narrows the plan to the PDE residual

The full scalar plan contains objective, metric, constraint, transformation, and
residual information.

The residual assembler should not need to inspect all of it.

The compiler therefore projects

```text
ScalarLoweringPlan
        │
        ▼
ScalarResidualAssemblyPlan
```

containing only

```text
residual_terms
data_placements
Robin boundary IDs
```

This is a useful design pattern.

Instead of passing the entire compiler state into every constructor, each numerical
subsystem receives the smallest closed plan needed for its job.

### 13.1 The residual plan recognizes closed registrations

The residual slice currently recognizes two principal registrations.

The baseline registration is exactly

```text
diffusion_reaction
volume_source
volume_control
```

with each contribution appearing once.

The general scalar registration is exactly

```text
tensor_diffusion
conservative_transport
advective_transport
reaction
volume_source
volume_control
robin_bilinear
robin_source
```

The registration logic is intentionally strict.

A plan containing

```text
diffusion_reaction
volume_source
volume_control
robin_source
```

is not silently treated as the baseline plus "one extra thing".

It needs a registered residual realization that explains the whole combination.

## 14. `ScalarServicePlan` narrows the plan to objective and optimization services

The complementary projection is

```text
ScalarLoweringPlan
        │
        ▼
ScalarServicePlan
```

which retains

```text
observations
losses
metric
constraint
transformation
selected tracking / flux / sensor data
```

This keeps objective and optimization-service construction separate from PDE residual
assembly.

That separation has a concrete mathematical consequence.

Changing

$$
\mathcal O_{y}(y)=y
$$

to a point or subdomain observation should change the objective/adjoint-source
services without automatically changing the state PDE residual.

Likewise, changing the search metric from

$$
G_{L^{2}}
$$

to another registered metric should not rewrite the scalar objective.

The split plan preserves those independent axes.

## 15. Planning is still pure compiler logic

At this stage, the plan can be built from the resolved semantic graph without
assembling deal.II operators.

The sequence is

```text
ResolvedProblemView
        │
        ▼
DealiiScalarLoweringPlanner
        │
        ├── residual handlers
        ├── observation handlers
        ├── loss handlers
        ├── metric handler
        ├── constraint handler
        └── transformation handler
        │
        ▼
ScalarLoweringPlan
```

The planner returns

```text
ScalarPlanResult
    diagnostics
    optional<ScalarLoweringPlan>
```

so an unsupported component can fail before model construction.

This is useful because a missing lowerer is a compiler-planning problem, not a
finite-element assembly exception.

## 16. The plan is created only for scalar-component targets

`DealiiCompiler` first resolves the target family.

It then decides whether the selected family uses the scalar component planner.

In the current implementation,

```text
direct scalar volume
or
assembled scalar v1 target
```

selects the scalar planning path.

The compiler invokes

```cpp
auto planned = scalar_planner_.plan(*resolution.problem);
```

and merges the plan diagnostics into the main compilation report.

If planning succeeds, the resulting `ScalarLoweringPlan` is retained for the scalar
model constructor.

Specialized target families bypass this specific IR but still obey the same earlier
semantic-resolution and capability checks.

## 17. Numerical construction starts where the plan meets the mesh

The plan still contains only semantic IDs, selected operator kinds, regions, data
placements, and realization choices.

Finite-element construction needs concrete geometry.

For the baseline path, the next layer must determine objects such as

```text
state finite element
state/test DoFHandler
independent state coordinates
control representation
quadrature
Dirichlet constraints
mass and stiffness operators
control coupling
observation operators
```

Those objects depend on

```text
ScalarLoweringPlan
+
Triangulation
+
DealiiDiscretisationPolicy
+
DealiiDataBindings
```

This is the point at which lowering becomes a concrete deal.II realization.

Section 36 returns to the resulting coordinate maps, observations, metrics, and
solve services after the compiler packaging path is established.

## 18. Semantic roles become realized spaces and maps

The compiler ultimately needs to turn declarations like

```text
state_space
control_space
state_observation_space
```

into concrete dimensions, layouts, and maps.

The compiled manifest therefore records realized spaces separately from semantic
spaces.

A `CompiledRealizedSpaceRecord` contains facts such as

```text
semantic_id
realization_id
dimension
layout
ordering
pairing_id
```

and a `CompiledRealizedMapRecord` records maps between them.

This distinction is important for constrained coordinates.

Semantically we may have one state variable $y$.

Numerically there can be

```text
independent state coordinates
        │
        ▼
reconstruction map
        │
        ▼
physical FE state
```

with different dimensions.

The compiled representation preserves both instead of pretending that a semantic
space ID is already a concrete vector layout.

## 19. The finite-element model stays private to the compiler

Once the target-specific or scalar-plan construction succeeds, the compiler has a
concrete numerical model.

For the scalar component path, this is an internal assembled model.

That internal class is not the public application abstraction.

The important boundary is

```text
private compiler model
    concrete deal.II spaces
    matrices
    coordinate maps
    assembly state
        │
        ▼
common contract adapters
        │
        ├── ExecutableModelT
        ├── MetricT
        ├── optional ConstraintT
        ├── state/adjoint solve services
        └── formulation-specific products
```

This prevents every solver from depending on the details of each lowerer.

A new compiler target can have a different private implementation as long as it
supplies the accepted contract surfaces.

## 20. The executable model packages actions, not the assembly recipe

By the time a numerical target becomes an `ExecutableModelT`, the solver does not
need to know how the residual was assembled.

For the running problem it needs operations corresponding to

```math
E_{h}(y,u)
```

and its derivatives:

```math
D E_{h}(y,u)[\delta y,\delta u],
```

```math
D E_{h}(y,u)^{\ast} r.
```

It also needs objective value/derivative services.

Thus the information flow is

```text
semantic residual terms
        ↓
ScalarResidualContribution[]
        ↓
FE assembly / operator construction
        ↓
ExecutableModelT actions
```

The semantic IDs and lowering provenance remain available in the manifest, while the
runtime solver sees the operational interface.

## 21. State and adjoint solves are packaged beside the executable model

Chapter 5, [Reduced state–adjoint formulation](05-reduced-state-adjoint-formulation.md), distinguished operator actions from solve services.

Compilation preserves that distinction.

The executable model can apply

```text
residual
JVP
VJP
```

while state/adjoint solver services perform inverse operations such as

```math
D_{y}E_{h}(y,u)\delta y=b
```

or

```math
D_{y}E_{h}(y,u)^{\ast} p = r_{y}.
```

`CompiledProblemT` therefore retains a `StateAdjointSolversT` alongside the executable
model.

This is important for reduced optimization:

```text
ExecutableModelT
    gives E, J, JVP, VJP

StateAdjointSolversT
    gives solve_state / solve_adjoint

ReducedDTOT
    composes both into reduced value/derivative operations
```

The compiler does not force inverse solves into the residual/JVP interface.

## 22. The metric is a separate compiled service

The same separation holds for the optimization metric.

For the baseline control,

$$
G:U_{h}\longrightarrow U_{h}^{\ast}
$$

is realized as the registered $L^{2}$ metric service.

`CompiledProblemT` stores the metric separately from the executable PDE model.

That means

```text
state equation
    does not own the optimization metric

objective
    does not own the optimization metric

solver
    can request metric apply/inverse through MetricT
```

This is why changing the search metric can remain independent of the PDE residual
when the compiler has a registered realization for the new metric.

## 23. Constraints remain optional compiled services

An unconstrained problem has no compiled constraint object.

If the semantic graph selects a registered cellwise or facewise box and concrete
bound data are supplied, compilation additionally creates the corresponding
`ConstraintT` service.

For a cellwise box,

```math
\ell_{i}\leq u_{i}\leq r_{i},
```

the compiled problem may also retain the concrete bound data as one shared source of
truth.

That source can then be reused by

```text
projection
complementarity
PDAS
manifest / provenance
```

rather than reconstructing bounds independently in several subsystems.

The constraint is therefore a separate service attached to the compiled product, not
a modification of the residual operator.

## 24. `CompiledProblemT` is the convergence point of the compiler path

After the numerical model and services are ready, the generic reduced/compiler
product contains the main pieces:

```text
CompiledProblemT
├── ExecutableModelT
├── MetricT
├── optional ConstraintT
├── StateAdjointSolversT
├── optional ReducedHessianT
├── optional compiled box data
├── optional CompiledApplicationViewT
├── CompilationManifest
└── retained lifetime owner
```

The constructor checks that these pieces are compatible.

For example, the control layout seen by the metric must match the control block of the
executable model.

If a constraint is present, its layout must match that same control block.

This is a final contract check after lowering: the compiler has produced several
services independently, but they must describe one coherent discrete problem.

## 25. Reduced DTO is created from the compiled services

For the reduced path, `CompiledProblemT` can construct a `ReducedDTOT`.

Conceptually the factory performs

```text
compiled executable model
        +
state/control partition
        +
state/adjoint solves
        +
lifetime owner
        │
        ▼
ReducedDTOT
```

The resulting formulation then implements the state/adjoint sequence developed in
Chapter 5, [Reduced state–adjoint formulation](05-reduced-state-adjoint-formulation.md).

This is where the compiler path joins the generic formulation layer:

```text
semantic graph
    ↓
compiler
    ↓
CompiledProblemT
    ↓
ReducedDTOT
    ↓
reduced optimizer
```

The reduced solver does not need to know whether the model originated from
`ScalarLoweringPlan`, a Dirichlet-control strategy, or another compiler target.

## 26. Other compilation products package different contract surfaces

The compiler can also package products whose solver-facing structure is not reduced
DTO.

For selected registered targets:

```text
CompilationProduct::quadratic_kkt
    ↓
CompiledQuadraticKKTProblemT
    ↓
EqualityConstrainedQuadraticKKTProductT

CompilationProduct::pdas
    ↓
CompiledPDASProblemT
    ↓
KKT product + complementarity + metric + box data + PDAS policies
```

Supplied OTD is represented through the semantic formulation/provenance and its
dedicated compiled records/actions.

The important point is that the numerical lowerer is reused where possible, but the
final packaging is formulation-specific.

Compilation is therefore not

```text
build one universal object
```

but

```text
build one coherent numerical realization
        │
        ├── package reduced services
        ├── package KKT services
        └── package PDAS / supplied-OTD services
```

subject to the capability checks from Chapter 10.

## 27. The manifest records what was actually built

Compilation produces more than executable callbacks.

It also records a typed `CompilationManifest`.

The manifest contains a `ResolvedCompilationDecision`, which can retain:

```text
semantic problem / formulation / target identity

mesh
    dimension
    active cells
    provenance
    lifetime policy
    structural identity

regions / spaces / pairings

runtime bindings
    semantic role
    concrete representation
    evaluation realization
    provenance

realized residuals / observations / losses / transformations

realized spaces and maps

state and adjoint solve policies

metric and constraint realization

formulation-specific KKT / PDAS / supplied-OTD records

typed boundary / transposition / metric selections

declared assumptions
```

This is not a second executable configuration.

It is evidence of the choices already made by compilation.

### 27.1 Provenance follows the semantic object into the compiled record

Consider the forcing datum.

Its path is roughly

```text
DataSpec
    forcing
        │
        ▼
ResolvedDataBindingRequest
    forcing
    volume_quadrature
    dealii::Function<dim>
        │
        ▼
concrete DealiiDataBindings
        │
        ▼
CompiledBindingRecord
    semantic_id = forcing
    runtime representation
    provenance
    value/digest status
```

The compiler can therefore explain not only that a state operator exists, but which
semantic input and concrete runtime source participated in its construction.

## 28. The manifest distinguishes semantic spaces from realized spaces

A semantic space says what role/topology was requested.

A compiled space record adds concrete realization facts such as

```text
finite element
dimension
runtime role
region
```

The manifest can also record realized coordinate maps separately.

For a constrained state this distinction becomes

```text
semantic state space
        │
        ▼
independent realized state space
        │
        ▼
reconstruction map
        │
        ▼
physical realized state space
```

That record is important because two compilations can share the same semantic problem
while differing in mesh or finite-element realization.

The semantic graph alone cannot distinguish those computations.

## 29. Mesh identity belongs to compiled evidence

The mesh is not stored in `ProblemSpec`.

It enters compilation later.

The compiled manifest records both human-facing mesh provenance and a structural mesh
identity.

This distinction prevents two meshes with the same descriptive label from being
treated as the same discrete realization.

Conceptually:

```text
semantic problem
    unchanged

mesh A
    ↓
compiled realization A

mesh B
    ↓
compiled realization B
```

Both runs may represent the same mathematical problem but not the same discrete
problem.

## 30. Borrowed mesh and owned compilation session are different lifetime modes

The compiler supports two entry patterns.

One takes a caller-owned `Triangulation`.

The other takes an owned/shared `DealiiCompilationSession`.

These imply different mesh lifetime policies.

The compiled manifest records the distinction, while the compiled product retains a
lifetime owner when one is needed.

The deeper ownership model appears in Chapter 13, but one consequence matters here:

> compilation may produce generic type-erased solver services without erasing the
> lifetime obligations of the concrete numerical objects behind them.

## 31. `CompiledApplicationViewT` preserves a narrow typed application seam

The generic solver-facing interfaces intentionally hide deal.II-specific details.

A compiled application may still need operations such as

```text
report physical / independent dimensions
split objective components
write native FE output
```

`CompiledApplicationViewT` provides this optional seam beside the erased solver view.

It is deliberately small:

```text
dimensions
    physical_state
    independent_state
    physical_control
    independent_control
    realized_observation

objective_components(full_point)

write_native_output(...)
```

This does not turn `ExecutableModelT` into an application/output API.

Instead:

```text
solver-facing generic view
    ExecutableModelT / MetricT / formulation products

application-facing typed seam
    CompiledApplicationViewT
```

both refer to the same compiled numerical realization.

## 32. One complete trace through the baseline problem

The entire transition can now be followed without introducing another abstraction.

Start from the state residual:

```math
E(y,u)[v]
=
\int_{\Omega}
\kappa\nabla y\cdot\nabla v
+
cyv
\mathrm{d}x
-
\int_{\Omega}
fv
\mathrm{d}x
-
\int_{\Omega}
uv
\mathrm{d}x.
```

### 32.1 Semantic graph

The graph says

```text
state_equation
    diffusion_reaction
    volume_source
    volume_control
```

with data ports

```text
diffusion
reaction
forcing
```

and variable ports

```text
state
control
```

### 32.2 Resolved semantic graph

`ResolvedProblemView` turns those string references into validated component lookups.

No numerical objects are created.

### 32.3 Closed compiler request

`ResolvedCompilationRequest` identifies the scalar target family and concrete binding
requirements.

For the baseline problem, the forcing and desired state become required runtime data
ports evaluated at volume/observation quadrature.

### 32.4 Scalar lowering plan

The planner translates the semantic terms:

```text
diffusion_reaction
    ↓
ScalarResidualOperatorKind::diffusion_reaction

volume_source
    ↓
ScalarResidualOperatorKind::volume_source

volume_control
    ↓
ScalarResidualOperatorKind::volume_control
```

and records data placements and handler provenance.

### 32.5 Residual assembly slice

The residual plan recognizes the exact three-term registration

```text
diffusion_reaction
volume_source
volume_control
```

and exposes only the data needed by this PDE assembly.

### 32.6 Service slice

The service plan independently selects

```text
state volume observation
control volume observation
quadratic tracking
quadratic control regularisation
cellwise L2 metric
optional cellwise box
```

### 32.7 Numerical model

The scalar lowerer constructs the concrete FE spaces, coordinate maps, matrices,
observation actions, metric, and state/adjoint solve services.

At this point the mathematics has become executable.

### 32.8 Contract adaptation

The private model is exposed through

```text
ExecutableModelT
MetricT
optional ConstraintT
StateAdjointSolversT
```

and optionally the application-native view.

### 32.9 Product packaging

For a reduced request,

```text
CompiledProblemT
    ↓
make_reduced_dto()
    ↓
ReducedDTOT
```

produces the solver-facing formulation.

### 32.10 Evidence

The manifest records the target, mesh, spaces, binding provenance, realized operators,
metric, solve policies, and semantic-to-handler decisions used to produce that
formulation.

The same problem has therefore passed through several representations:

```text
mathematical residual
        ↓
semantic residual terms
        ↓
typed lowering contributions
        ↓
finite-element operators
        ↓
generic executable actions
        ↓
reduced formulation
```

Each representation serves a different layer; none is merely a renamed copy of the
previous one.

## 33. A specialized target follows the same outer pipeline

Consider a Dirichlet-control problem.

Its physical state may satisfy schematically

```math
y\rvert_{\Gamma_{C}}=u.
```

The semantic graph declares a controlled-Dirichlet lifting and the required boundary
spaces, loss, metric, and realization policies.

The compiler still performs

```text
semantic validation
    ↓
semantic resolution
    ↓
ResolvedCompilationRequest
    ↓
lowerability / capability checks
```

but after target resolution it selects a dedicated Dirichlet-control model strategy
rather than the scalar component plan used by the baseline volume problem.

The outer architecture therefore remains the same:

```text
same semantic language
same diagnostic vocabulary
same compiler request closure
different private realization strategy
same contract/formulation boundary
```

This is why target-specific implementation does not imply target-specific public
solver APIs.

## 34. Compilation does not consume semantic display text

Several objects contain labels or descriptive prose for humans.

Those strings are useful in documentation and artifacts but are not meant to drive
the numerical branch selection.

Lowering consumes typed facts such as

```text
ResidualTermKind
ObservationKind
MetricKind
RequirementKind
typed realization selections
region IDs
space IDs
```

and compiler-owned target enums.

The same principle appeared in Chapter 10 for `RequirementPolicySpec`:

> descriptive text explains a decision; typed fields determine it.

Compilation extends that rule to the whole lowering path.

## 35. `CompiledProblemT` is not the universal nmopt integration type

The compiler path ends in `CompiledProblemT` because the compiler has responsibility
for constructing and packaging the numerical realization.

An external application may already own all of those objects.

That path can bind its existing residual/objective actions, solves, and metric
directly to the same contract vocabulary without manufacturing a `ProblemSpec` or
`CompiledProblemT`.

Thus

```text
semantic/compiler path
    ProblemSpec
        ↓
    CompiledProblemT
        ↓
    formulation contracts

native application path
    existing numerical objects
        ↓
    binding/adapters
        ↓
    formulation contracts
```

The convergence point is the numerical/formulation contract layer, not
`CompiledProblemT` itself.

Chapter 13 develops that second path in detail.

## 36. What the compiled realization actually contains

The compiler story ends only when the selected semantic roles have become concrete
numerical actions.

For the baseline scalar problem, the important coordinate relation is

```math
\mathbf y_{\mathrm{phys}}
=
P\mathbf z+\boldsymbol\ell,
```

where $\mathbf z$ contains independent state coordinates, $P$ embeds them into the
physical finite-element vector, and $\boldsymbol\ell$ carries fixed boundary data.

The corresponding tangent and covector maps are

```math
\delta\mathbf y_{\mathrm{phys}}
=
P\delta\mathbf z,
```

and

```math
\mathbf r
=
P^{\mathsf T}\mathbf r_{\mathrm{phys}}.
```

The compiled scalar realization exposes these roles through the concrete
`IndependentStateCoordinates` service:

```text
reconstruct(z)
    P z + ell

embed(z)
    P z

pullback(r_phys)
    P^T r_phys
```

This closes the transition introduced in Chapter 2: the semantic state variable and
the physical finite-element state need not use the same coordinates.

### 36.1 Physical operators become independent-coordinate operators

For the linear residual

```math
\mathbf E_{\mathrm{phys}}
=
A\mathbf y_{\mathrm{phys}}
-
B_{\mathrm{phys}}\mathbf u
-
\mathbf f,
```

substituting the reconstruction and pulling the residual back gives

```math
\mathbf E(\mathbf z,\mathbf u)
=
\underbrace{P^{\mathsf T}AP}_{\widehat A}\mathbf z
-
\underbrace{P^{\mathsf T}B_{\mathrm{phys}}}_{\widehat B}\mathbf u
-
P^{\mathsf T}(\mathbf f-A\boldsymbol\ell).
```

Thus the compiled state and control operators are naturally

```math
\widehat A
=
P^{\mathsf T}AP,
\qquad
\widehat B
=
P^{\mathsf T}B_{\mathrm{phys}}.
```

The compiler-owned scalar model constructs these reduced-coordinate actions from the
physical operators rather than assuming that constrained DoFs can simply be dropped.

### 36.2 Observations also need transpose realizations

A discrete observation has the form

$$
O_{h}:
\mathbb R^{n_{\mathrm{phys}}}
\longrightarrow
\mathbb R^{m}.
$$

For a quadratic loss,

```math
J_{\mathrm{obs}}
=
\frac12
\lVert
O_{h}\mathbf y_{\mathrm{phys}}-\mathbf d
\rVert_{W}^{2},
```

the independent-state covector is

```math
P^{\mathsf T}
O_{h}^{\mathsf T}
W
(O_{h}(P\mathbf z+\boldsymbol\ell)-\mathbf d).
```

That formula explains why an observation realization needs both forward and
transpose actions.

For point sensors,

```math
\mathcal O_{\mathrm s}(y)
=
[y(x_{1}),\ldots,y(x_{m})]^{\mathsf T},
```

the compiler stores one evaluation vector per sensor. Forward evaluation takes dot
products with the physical state; the VJP accumulates the corresponding weighted
evaluation vectors and then pulls the physical covector back with $P^{\mathsf T}$.

The same pattern applies to boundary traces and normal-flux observations, with
different evaluation operators.

### 36.3 Metrics are compiled Riesz services

The metric surface is not restricted to one stored matrix.

For an $L^{2}$ control metric,

$$
G=M,
$$

so `apply()` is a mass action and `inverse_apply()` solves

$$
Mg=r.
$$

For the registered $H^{-1}$ realization,

```math
G
=
MK^{-1}M,
```

so even the forward Riesz action contains an internal solve.

For the trace $H^{1/2}$ realization,

```math
G
=
A_{BB}
-
A_{BI}A_{II}^{-1}A_{IB},
```

and `apply()` performs the minimum-energy interior extension rather than explicitly
forming a dense Schur complement.

All three satisfy the same `MetricT` interface even though their numerical
realizations are substantially different.

### 36.4 State and adjoint solves remain separate services

Residual evaluation computes

```math
\widehat A\mathbf z
-
\widehat B\mathbf u
-
\widehat{\mathbf f}.
```

A state solve instead finds

```math
\widehat A\mathbf z
=
\widehat{\mathbf f}
+
\widehat B\mathbf u.
```

The adjoint solve uses

```math
\widehat A^{\mathsf T}\mathbf p
=
\mathbf r_{y}.
```

For symmetric diffusion–reaction problems, state and adjoint solves may reuse the same
SPD operator and solver implementation. Once transport makes the state Jacobian
nonsymmetric, the adjoint service must use the transpose operator.

`CompiledProblemT` therefore keeps residual/JVP/VJP actions separate from
state/adjoint solve services. The formulation can rely on the mathematical
distinction without depending on whether the underlying implementation reuses one
matrix or one solver object.

### 36.5 Backend algebra is deliberately narrower

The serial deal.II backend, `SerialBackend`, supplies primitive vector operations such as

```text
zeros
size
dot
add_scaled
scale
```

on `dealii::Vector<double>`.

It does not decide what a state, control, observation, metric, or boundary condition
means.

Those meanings belong to the numerical realization assembled by the compiler or
provided by an external application.

This distinction is why the same backend can support several finite-element
realizations without becoming a PDE model itself.

## 37. Following compilation through the source

The compiler entry point is:

- [`include/nmopt/compiler/v1/dealii_compiler.hpp`](../../../include/nmopt/compiler/v1/dealii_compiler.hpp)

A useful reading sequence is:

```text
validate(...)
compile(...)
resolve_compilation_request(...)
close_compilation_request(...)
compile_impl(...)
```

Inside `compile_impl()`, look for the target-family decision and the
`uses_scalar_component_target` branch.

For the compiler-side request and runtime binding vocabulary, read:

- [`include/nmopt/compiler/v1/dealii_types.hpp`](../../../include/nmopt/compiler/v1/dealii_types.hpp)

The central types for this chapter are:

```text
ResolvedCompilationRequest
ResolvedDataBindingRequest
DealiiDataBindings
DealiiDiscretisationPolicy
CompilationProduct
```

For the bounded component path, read:

- [`include/nmopt/compiler/v1/dealii_scalar_plan.hpp`](../../../include/nmopt/compiler/v1/dealii_scalar_plan.hpp)

Follow:

```text
DealiiScalarLowererRegistryV1
        ↓
DealiiScalarLoweringPlanner
        ↓
ScalarLoweringPlan
        ├── ScalarResidualAssemblyPlan
        └── ScalarServicePlan
```

The generic compiler product and typed manifest live in:

- [`include/nmopt/compiler/v1/compiled_problem.hpp`](../../../include/nmopt/compiler/v1/compiled_problem.hpp)

The key distinction there is between

```text
executable runtime services
```

and

```text
CompilationManifest / ResolvedCompilationDecision
```

which records evidence about those services.

The optional application-facing seam is:

- [`include/nmopt/compiler/v1/compiled_application_view.hpp`](../../../include/nmopt/compiler/v1/compiled_application_view.hpp)

For the exact registered target families and their bounded capabilities, use:

- [`docs/internals/compiler/semantic-compiler.md`](../../internals/compiler/semantic-compiler.md)

That ledger remains authoritative for which complete semantic signatures the current
compiler accepts.

## 38. Part III closes at the numerical contract boundary

The three Part III chapters now describe one complete producer path:

```text
Chapter 9
    mathematical problem
        ↓
    ProblemSpec

Chapter 10
    ProblemSpec
        ↓
    validated / resolved / supported request

Chapter 11
    accepted request
        ↓
    lowering plan / registered target
        ↓
    discrete operators / solves / metrics
        ↓
    CompiledProblemT
```

At this point the compiler path has reached the same common numerical/formulation
contracts used by a native application.

Part IV turns both producer paths into user-facing execution stories. Chapter 12,
**Authoring and using compiled problems**, shows how to drive the compiler path from
recipes and scenarios through solving. Chapter 13, **Integrating an existing PDE
application**, shows how an application can supply the same numerical contracts while
retaining ownership of its mesh, matrices, solvers, and output model.

## Read later

Useful existing documents are:

- [Describing and compiling a problem](../overview/semantic-compiler.md), for the
  shorter architecture view.
- [Numerical realization](../overview/numerical-realization.md), for the
  shorter project-wide view of the numerical layer.
- [v1 semantic graph and deal.II compiler](../../internals/compiler/semantic-compiler.md),
  for the exact registered capability and target ledger.
- [Validation, resolution, and capabilities](10-validation-resolution-and-capabilities.md),
  for the acceptance boundaries that precede this chapter.
- [Semantic problem model](09-semantic-problem-model.md), for the semantic vocabulary
  lowered here.

Compilation is the translation layer between an accepted semantic request and the
concrete numerical services studied next.
