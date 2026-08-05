# V1 semantic graph and deal.II compiler

## Status and boundary

This is the first public semantic-to-compiler path. It is deliberately named
**v1** and exists alongside the direct v0 model:

```text
v0 direct reference
ScalarDiffusionReactionModel<dim>

v1 semantic/compiler path
semantic::v1::ProblemSpec
  -> SemanticValidator + compiler diagnostics
  -> compiler::v1::CompiledProblemT<SerialBackend>
  -> generic executable/metric/constraint/DTO services
```

The v0 lowerer is not modified or replaced by this path. The homogeneous v1
reference graph privately constructs a separate v0 executable instance and
packages it behind generic executable ports, so v0 and v1 can be compared at
the same coefficients. A graph that declares fixed-Dirichlet reconstruction,
material-subdomain state tracking, or Neumann boundary control selects a
separate v1-only assembled target; it never extends the v0 direct model. All
targets expose the same generic compiled ports.

The normative semantic protocol remains the
[interface specification](interface-specification.md). This document records
the concrete, intentionally narrow v1 realization of that protocol.

## Public semantic graph

The compatibility aggregate `include/nmopt/semantic/v1/problem_spec.hpp`
includes the focused deal.II-free headers `types.hpp`, `validation.hpp`, and
`reference_specs.hpp`. The last contains
`make_scalar_diffusion_reaction_problem()`, the homogeneous comparison graph.
`make_fixed_dirichlet_scalar_diffusion_reaction_problem()` adds the first
declared physical-field transformation, while
`make_subdomain_tracking_scalar_diffusion_reaction_problem()` adds a named
material-id observation region. `make_neumann_boundary_control_problem()`
declares a facewise Neumann control and state boundary trace.
`make_pure_neumann_boundary_control_problem()` is its separate zero-reaction
mean-constraint variant. `make_h1_regularised_scalar_diffusion_reaction_problem()`
is the separate continuous-control objective variant.
`make_coefficient_identification_problem()` adds a physical diffusion
parameter. The compatibility table below describes the homogeneous graph; the
named graph sections record the added features.

```text
Region       one full volume region, named material-id volume subregions, and Dirichlet boundary ids
Space        scalar H1 state/test; scalar L2 control and observations
Pairing      explicit coefficient pairings for state, test, control, and observations
Variable     one state and one control; the state may name a physical-field transformation
Data         forcing, desired state, fixed Dirichlet lifting, diffusion, reaction,
             regularisation, and optional lower/upper cellwise bounds
Transformation optional fixed-Dirichlet reconstruction of the physical state
Residual     diffusion-reaction, volume source, and volume control
Observation  full-volume control restriction; state restriction on the full volume or one material subregion
Loss         quadratic tracking and quadratic control regularisation
Metric       cellwise L2 control metric
Constraint   optional cellwise L2 box
Formulation  one-state/one-control reduced DTO
```

Concrete values are not semantic objects. `DealiiDataBindings<dim>` binds the
forcing and desired-state `Function` objects plus scalar coefficients only
after validation; its constant-diffusion binding is optional and is required
only by a graph that declares constant diffusion. A graph with the
reconstruction also requires its optional `fixed_dirichlet_data` binding. The
selected data rule interpolates that
`Function` at the declared Dirichlet boundary DoFs to form the lifting; it is
not inferred from the forcing or target. If the graph declares its optional
box, the compiler also requires `CellwiseBoxDataBindings`: both bounds must
be scalar constants or both must be exact-layout `FE_DGQ(0)` coefficient
vectors.

The tracking target has its own required selected policy. The current
registered realization is an analytic desired-state `Function` evaluated at
the compiler-selected `QGauss` volume quadrature points, on exactly the
tracking observation region. The validator rejects a missing selected policy
or one declared on a different region. There is no implicit nodal
interpretation. A future FE projection or interpolation must be represented
by a new declared policy and lowerer with its own binding and cache boundary.

### Material-subdomain tracking

`RegionSpec::material_ids` names the active-cell material ids of a non-full
volume region. The selected state observation owns that region; its observation
space and desired-state data declaration must name the same region. The v1
assembled target builds the tracking mass, desired-state load, and target norm
only on matching cells. It continues to assemble the state matrix, forcing,
control coupling, control mass, and state/adjoint solvers over the full mesh.
Consequently, changing the tracking region changes the objective derivative
and adjoint right-hand side without changing the state equation or solver
interfaces.

### $H^{1}$ control regularisation and search metric

The first half of P2.3 adds a distinct
`quadratic_h1_control_regularisation` loss. It changes the control space to
continuous scalar `FE_Q` on the state mesh and assembles

```math
J_{\mathrm{control}}(u_{h})=\frac{\alpha}{2}u_{h}^{T}(M_{u}+K_{u})u_{h},
```

where $M_{u}$ is the control mass matrix and $K_{u}$ its Laplace stiffness
matrix. Thus the objective derivative receives
$\alpha(M_{u}+K_{u})u_{h}$. `make_h1_regularised_scalar_diffusion_reaction_problem()`
deliberately retains the `l2_continuous` search metric, which maps covectors
through $M_{u}^{-1}$.

`make_h1_metric_scalar_diffusion_reaction_problem()` selects a separate
`MetricKind::h1` realization on the same continuous control layout. It uses
$G=M_{u}+K_{u}$ as the Riesz map, so the positive mass term makes the map
coercive without a separate boundary or mean condition. Its CG inverse is
used only to form a search direction: selecting it does not change the
objective, residual, state equation, or adjoint equation. The first target
supports the homogeneous full-domain volume-control graph only and selects no
box: the existing coefficientwise boxes apply only to the discontinuous
cellwise or facewise control layouts.

### Coefficient identification

`make_coefficient_identification_problem()` uses the existing binary reduced
DTO decision port for a `VariableRole::parameter` named
`diffusion_parameter`. It has a cellwise `FE_DGQ(0)` parameter space and a
separate parameter observation, $L^{2}$ metric, and cellwise box. The selected box
must bind a strictly positive lower bound, so the physical coefficient is
positive without an implicit $m=\exp(q)$ transformation.

The graph replaces the constant-diffusion data port and source-control term by
the `parameter_diffusion_reaction` residual term:

```math
\langle E(y,m),v\rangle=(m\nabla y,\nabla v)+(c y,v)-(f,v).
```

Its private `dealii_coefficient_identification.hpp` target reassembles the
parameter-dependent state matrix for each residual, state solve, and adjoint
solve. It supplies $`D_{y}E(y,m)\delta y`$, $`D_{m}E(y,m)\delta m`$, and their
exact transpose actions. The parameter regularisation is
$`\frac{\alpha}{2} m_{h}^{T}M_{m}m_{h}`$; the parameter $L^{2}$ Riesz map and coefficientwise
box use the `l2_cellwise_parameter` identifier. The target is limited to the
homogeneous, full-domain volume graph with a physical positive parameter;
continuous parameter spaces, logarithmic parameterisation, and parameter
spaces on other meshes remain unsupported.

### Neumann control and boundary tracking

The registered P2.1 realization assigns one scalar control coefficient to
each active state-mesh boundary face with an id in `control_boundary`. The
compiler assembles with `FEFaceValues`:

```math
r_{h}(y,u)=A_{h}y-f_{h}-C_{\Gamma}u,
\qquad
(C_{\Gamma}u)_{i}=\sum_{F\subset\Gamma_{c}} u_{F}\int_{F}\phi_{i}.
```

$C_{\Gamma}^{\ast}$ is used directly for the residual VJP. A separate marked
`observation_boundary` assembles the trace tracking mass, target load, and
target norm with face quadrature. The desired-state `Function` is evaluated
at the selected `QGauss` face quadrature points; no nodal trace or implicit
surface projection is inferred.

The control face indicators have disjoint support, so their boundary $L^{2}$
mass matrix is diagonal with face measures. `MassMetric` receives the
`l2_facewise` identifier, while the separate `FacewiseBoxConstraint` clips
only these face coefficients. `FacewiseBoxDataBindings` is intentionally
separate from `CellwiseBoxDataBindings`; both bounds are scalar constants or
exact face-layout vectors. The compiler rejects a control boundary that
overlaps the fixed homogeneous Dirichlet ids.

### Pure Neumann mean constraint

The P2.2 graph retains the singular natural-boundary residual but declares
`mean_zero_multiplier` on the full volume instead of fixed Dirichlet rows.
The private v1 target solves state and adjoint systems with

```math
\begin{bmatrix}A_{h}&m_{h}\\m_{h}^{T}&0\end{bmatrix}
\begin{bmatrix}y_{h}\\\lambda \end{bmatrix}=
\begin{bmatrix}f_{h}+C_{\Gamma} u\\0\end{bmatrix},
\qquad (m_{h})_{i}=\int_{\Omega}\phi_{i}.
```

It requires zero reaction and checks the discrete constant-mode pairing of
forcing at compilation and of every boundary-control state load at solve time.
The adjoint uses the same augmented system, so state and adjoint have the
recorded zero mean without pinning a DoF. The manifest names the multiplier
policy and serial `SparseDirectUMFPACK` saddle solve. The pure-Neumann graph
does not select a box: the current generic projection cannot preserve its
affine load-compatibility condition.

### Fixed essential reconstruction

The fixed-data graph represents independent state coordinates and a physical
field explicitly:

```math
y_{\mathrm{phys}}=P_{h}\widehat y_{h}+\ell_{0,h}.
```

$P_{h}$ is constructed from homogeneous `AffineConstraints`; the separate
constraint object with the bound `Function` yields $`\ell_{0,h}`$. Residual
and tracking assembly evaluate $`y_{\mathrm{phys}}`$; the state covector,
objective derivative, residual JVP, and residual VJP apply the matching
$`P_{h}^{\ast}`$ pullback. State and adjoint CG solves use the compiled
independent-coordinate system. The compiled model is immutable: compiling
with a different lifting or other data binding assembles a distinct product,
so no data-dependent cached field is shared across compilations.

## Validation and diagnostics

`SemanticValidator` validates semantic structure and declared policies.
`compiler::v1::DealiiCompiler::validate()` appends compiler-specific checks
to the same `ValidationReport`.

| Category | Produced by | Examples in v1 |
| --- | --- | --- |
| `structural` | `SemanticValidator` | missing ports, absent equation test space, wrong term inputs |
| `analytical_policy` | `SemanticValidator` | missing selected fixed-Dirichlet or cellwise-bound policy |
| `lowerability` | `DealiiCompiler` and its `DealiiLowererRegistryV1` | matrix-free execution, zero `FE_Q` degree, unregistered node kind, missing bound or fixed-lifting binding |
| `formulation_capability` | `DealiiCompiler` | all-at-once formulation or a multi-block DTO request |

`CompilationResultT<Backend>` returns the report and only contains a
`CompiledProblemT<Backend>` when it is valid. Unsupported choices are reported
as diagnostics; the compiler does not substitute another residual,
observation, metric, constraint, or formulation.

## Registered deal.II realization

The compatibility aggregate
`include/nmopt/compiler/v1/dealii_scalar_diffusion_reaction.hpp` exposes the
focused compiler headers: `compiled_problem.hpp`, `dealii_types.hpp`,
`dealii_capabilities.hpp`, and `dealii_compiler.hpp`. The compiler's private
`dealii_fixed_dirichlet.hpp` target owns the separate v1 physical-state
assembly used for fixed reconstruction and material-subdomain tracking. The
private `dealii_neumann_boundary.hpp` target owns the distinct Neumann
residual, facewise control layout, boundary trace tracking, facewise metric,
facewise box realization, and pure-Neumann mean-zero saddle realization. The
private `dealii_h1_control.hpp` target owns the distinct continuous-control
$H^{1}$ loss and $H^{1}$ or $L^{2}$ search metrics. The private
`dealii_coefficient_identification.hpp` target owns the cellwise positive
diffusion-parameter residual, its nonlinear first-order actions, parameter
metric, and parameter box. The registry otherwise supports only the listed
volume terms, full-domain control observation, full-domain or material-id
state restriction, quadratic losses, `L2` metric, optional cellwise box, and
fixed-Dirichlet reconstruction. Its selected discrete policies are assembled
serial scalar `FE_Q` state/test with degree at least one and reduced DTO:
`FE_DGQ(0)` volume control on the state mesh, or one facewise-constant
Neumann coefficient for every marked boundary face, or continuous `FE_Q`
volume control for the registered $H^{1}$ loss. Pure Neumann is limited to
zero reaction and compatible forcing/control loads. Coefficient identification
instead selects a cellwise `FE_DGQ(0)` physical parameter with a strictly
positive box and reassembled state/adjoint matrices.

`CompiledProblemT<Backend>::executable_model()`, `metric()`, `constraint()`,
and `make_reduced_dto()` expose only backend-neutral ports and formulation
services. The homogeneous private v0 target, v1 assembled target, mass metric,
and box constraint stay inside `DealiiCompiler`. Solvers have no v1 branches.

Every successful compiled product also carries a `CompilationManifest`. It
records semantic component identities, FE spaces, quadrature, the
dual-coefficient representation, target-data rule, material observation
realization, metric solve tolerances, constraint/lifting/nullspace policies,
DTO provenance, declared transformations, and assumptions. The manifest is
descriptive provenance; it does not create a second configuration channel.

## Comparison guarantee

`tests/dealii_diffusion_contract.cc` creates one direct v0 model and one
homogeneous v1 compiled model on the same triangulation, functions, constants,
and bounds. It verifies equal assembled residual, objective, objective
derivative, and DTO reduced derivative. The same test confirms the optional
compiled box constraint and classifies unsupported matrix-free and all-at-once
requests. It also compiles a nonzero manufactured fixed-Dirichlet state and
checks its physical residual/objective, reconstruction JVP/VJP pairing, and
reduced Taylor remainder; recompilation with changed lifting data must change
the compiled result. The same test builds a mesh with distinct Dirichlet,
Neumann-control, and observation boundary ids. It checks the Neumann
residual pairing, trace-loss derivative, facewise metric and box realization,
and a state-recomputed reduced Taylor remainder.

`tests/semantic_v1_contract.cc` independently verifies the semantic validator
for the canonical and material-subdomain graphs, including missing or
region-mismatched target-data policies. The deal.II test compiles two material
masks on the same mesh and requires the same state solution and residual JVP,
while requiring a different tracking objective and state-objective derivative.
The same test compiles the coefficient-identification graph, rejects a
nonpositive lower bound, verifies reassembly at changed parameter values, and
checks parameter JVP/VJP, finite-difference, and reduced-Taylor identities.

## Exclusions

This v1 registration does not broaden the v0 executable mathematics. Beyond
the selected fixed-data reconstruction, material-id state tracking, and
marked-face Neumann control with boundary tracking, it does not compile
arbitrary geometric or overlapping subdomain restrictions, FE target
projection/interpolation, Robin terms, controlled Dirichlet liftings,
continuous-control bounds, continuous or transformed coefficient parameters,
matrix-free execution, all-at-once/OTD, multiple equations, or multiple
optimisation variables. Each requires its own semantic declaration, registered
lowerer, capability diagnostic, and value/JVP/VJP/reduced tests.
