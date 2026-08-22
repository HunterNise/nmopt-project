# V1 semantic graph and deal.II compiler

## Status and boundary

This is the first public semantic-to-compiler path. It is deliberately named
**v1** and exists alongside the direct v0 model:

```text
v0 direct reference
ScalarDiffusionReactionModel<dim>

v1 semantic/compiler path
semantic::v1::ProblemSpec
  -> SemanticValidator + SemanticResolver + compiler diagnostics
  -> bounded component lowering plan or registered target strategy
  -> compiler::v1::CompiledProblemT<SerialBackend>
  -> generic executable/metric/constraint/DTO services
```

The v0 lowerer is not modified or replaced by this path. The homogeneous v1
reference graph privately constructs a separate v0 executable instance and
packages it behind generic executable ports, so v0 and v1 can be compared at
the same coefficients. A graph that declares fixed-Dirichlet reconstruction,
controlled-Dirichlet lifting, material-subdomain state tracking, Neumann
boundary control, the C5.6 Neumann/transport/subdomain composition, or the
registered general scalar elliptic/Robin composition selects a separate
v1-only assembled target; it never extends the v0 direct model. All targets
expose the same generic compiled ports.

The normative semantic protocol remains the
[interface specification](../../design/interface-specification.md). This document records
the concrete, intentionally narrow v1 realization of that protocol.

## Registered capabilities

This section is the sole detailed capability ledger for the current v1
compiler. The [implementation roadmap](../../planning/implementation-roadmap.md)
owns completion and handoff status; design and policy documents link here
instead of reproducing this table.

Every factory below is validated in
[`semantic_v1_contract.cc`](../../../tests/semantic_v1_contract.cc). The last
column names its focused CTest scenario backed by either
[`dealii_diffusion_contract.cc`](../../../tests/dealii_diffusion_contract.cc)
or
[`dealii_trace_hhalf_metric_contract.cc`](../../../tests/dealii_trace_hhalf_metric_contract.cc).

| Registered semantic graph | Selected implementation | Bounded capability | Focused CTest scenario |
| --- | --- | --- | --- |
| `make_scalar_diffusion_reaction_problem()` | `ScalarDiffusionReactionModel<dim>` direct v0 reference, packaged through v1 ports | Homogeneous fixed Dirichlet, full-volume tracking, `FE_DGQ(0)` volume control, $L^{2}$ metric, and optional cellwise box | `nmopt.dealii.canonical_volume_control` |
| `make_l2_state_tracking_continuous_control_problem()` | `ContinuousControlModel<dim>` | Full-volume $L^{2}$ tracking and regularisation with independent homogeneous-Dirichlet conforming Lagrange volume control (`FE_Q` on hypercubes or `FE_SimplexP` on simplices), an $L^{2}$ metric, and no box | `nmopt.dealii.l2_tracking_continuous_control` |
| `make_fixed_dirichlet_scalar_diffusion_reaction_problem()` | `ScalarComponentModel<dim>` built from `ScalarLoweringPlan` | Fixed-data reconstruction with independent coordinates and optional cellwise box | `nmopt.dealii.fixed_dirichlet` |
| `make_subdomain_tracking_scalar_diffusion_reaction_problem()` | `ScalarComponentModel<dim>` built from `ScalarLoweringPlan` | State tracking on one material-id set while retaining the full-domain state equation | `nmopt.dealii.subdomain_observation` |
| `make_h1_state_tracking_scalar_diffusion_reaction_problem()` | `ScalarComponentModel<dim>` built from `ScalarLoweringPlan` | Full-domain $H^{1}_{0}$ state observation with mass-plus-stiffness tracking and an unchanged control $L^{2}$ metric | `nmopt.dealii.h1_state_observation` |
| `make_point_sensor_scalar_diffusion_reaction_problem()` | `ScalarComponentModel<dim>` built from `ScalarLoweringPlan` | C5.10 finite point-set observation with immutable physical coordinates, finite-dimensional `FE_Q` evaluation, and an assembled very-weak transpose point load | `nmopt.dealii.point_sensor` |
| `make_normal_flux_scalar_diffusion_reaction_problem()` | `ScalarComponentModel<dim>` built from `ScalarLoweringPlan` | C5.8 strong-state normal-flux tracking on selected boundary IDs, outward `FEFaceValues` normal derivatives at face quadrature, and the assembled very-weak boundary-source transpose | `nmopt.dealii.normal_flux` |
| `make_l2_metric_h1_state_tracking_continuous_control_problem()` and `make_hminus1_metric_h1_state_tracking_scalar_diffusion_reaction_problem()` | `ContinuousControlModel<dim>` | Same independent homogeneous-Dirichlet conforming Lagrange control layout and energy-tracking graph with separately selected $L^{2}$ or $H^{-1}$ metric; no box | `nmopt.dealii.hminus1_compilation` |
| `make_dirichlet_control_scalar_diffusion_reaction_problem()` | `DirichletControlLiftingModel<dim>` | Complete-exterior-boundary nodal trace lifting, trace $L^{2}$ metric, and no trace box | `nmopt.dealii.dirichlet_control` |
| `make_l2_dirichlet_laplace_control_problem()` | `DirichletControlLiftingModel<dim>` through the declared conforming-trace equivalence | Chapter 5.11.2 $L^{2}$ state/control transposition parent, $H^{2}\cap H^{1}_{0}$ tests, complete-boundary $U_{h}=\mathrm{tr}_{\Gamma}V_{h}\subset H^{1/2}(\Gamma)$, normalized Laplacian, and no trace box | `nmopt.dealii.l2_dirichlet_transposition` |
| `make_hhalf_dirichlet_laplace_control_problem()` | `DirichletControlLiftingModel<dim>` with `TraceHhalfMetric` | Section 5.11.1 option 1: normalized Laplacian, $L^{2}$ state tracking, and the minimum-extension $H^{1/2}$ action for both control loss and search metric | `nmopt.dealii.section_5_11_compilation` |
| `make_h1_tracking_hhalf_dirichlet_laplace_control_problem()` | `DirichletControlLiftingModel<dim>` with independently selected loss and metric | Section 5.11.1 option 2: physical $H^{1}$ state tracking, boundary $L^{2}$ control loss, minimum-extension $H^{1/2}$ search metric, and no trace box | `nmopt.dealii.section_5_11_compilation` |
| `make_h1_dirichlet_laplace_control_problem()` | `DirichletControlLiftingModel<dim>` with tangential face assembly | Section 5.11.3: normalized Laplacian and boundary $M_{\Gamma,h}+K^{\tau}_{\Gamma,h}$ for the separately declared $H^{1}$ control loss and metric | `nmopt.dealii.section_5_11_compilation` |
| `make_partial_dirichlet_control_scalar_diffusion_reaction_problem()` | `DirichletControlLiftingModel<dim>` with fixed and controlled trace maps | Complete fixed/controlled boundary partition, nonzero nodal fixed lifting, fixed-data precedence at interface DoFs, relative-interior controlled trace, trace $L^{2}$ metric, and no trace box | `nmopt.dealii.partial_dirichlet_control` |
| `make_neumann_boundary_control_problem()` | `NeumannBoundaryControlModel<dim>` | Marked-face Neumann control, boundary trace tracking, `FE_Q` or `FE_SimplexP` state/test coordinates, facewise $L^{2}$ metric, and optional facewise box | `nmopt.dealii.neumann_boundary` |
| `make_neumann_convection_subdomain_tracking_problem()` | `NeumannBoundaryControlModel<dim>` with conservative transport and volume observation | Mixed fixed/Neumann boundary partition, conservative-transport weak form, material-id state tracking, `FE_Q` or `FE_SimplexP` state/test coordinates, facewise $L^{2}$ control metric, and optional facewise box | `nmopt.dealii.neumann_convection_subdomain` |
| `make_weighted_boundary_trace_neumann_control_problem()` | `NeumannBoundaryControlModel<dim>` with fixed weight data | Marked-face Neumann control and the explicit map $y\mapsto h\gamma y$ with an unchanged facewise $L^{2}$ metric | `nmopt.dealii.weighted_boundary_trace` |
| `make_pure_neumann_boundary_control_problem()` | `NeumannBoundaryControlModel<dim>` with mean-zero gauge | Zero-reaction pure Neumann state and adjoint with compatible forcing and controls; no box | `nmopt.dealii.pure_neumann` |
| `make_h1_regularised_scalar_diffusion_reaction_problem()` and `make_h1_metric_scalar_diffusion_reaction_problem()` | `ContinuousControlModel<dim>` | Continuous conforming Lagrange control with $H^{1}$ loss and separately selected $L^{2}$ or $H^{1}$ metric; no box | `nmopt.dealii.h1_control` |
| `make_coefficient_identification_problem()` | `CoefficientIdentificationModel<dim>` | Positive cellwise physical diffusion parameter, reassembled state/adjoint operators, parameter $L^{2}$ metric, and cellwise box | `nmopt.dealii.coefficient_identification` |
| `make_general_scalar_elliptic_robin_problem()` | `ScalarComponentModel<dim>` built from `ScalarLoweringPlan` | Tensor diffusion, conservative and advective transport, reaction, volume source/control, and Robin bilinear/source terms with homogeneous fixed Dirichlet data | `nmopt.dealii.general_scalar_robin` |

### Closed Dirichlet registration table

The compiler resolves the following structural signatures before constructing
`DirichletControlLiftingModel<dim>`. The resolver consumes the resolved
component graph and typed policy selections; it does not inspect a factory or
problem identifier.

| Registration | State tracking | Control observation and loss | Search metric | Boundary realization |
| --- | --- | --- | --- | --- |
| Complete nodal lifting | volume $L^{2}$ | nodal boundary $L^{2}$ | $L^{2}$ | complete controlled boundary |
| Partial nodal lifting | volume $L^{2}$ | nodal boundary $L^{2}$ | $L^{2}$ | disjoint fixed/controlled partition with fixed-data precedence |
| Conforming transposition | continuous $L^{2}$ parent | boundary $L^{2}$ | $L^{2}$ | complete conforming trace subspace and typed transposition equivalence |
| Section 5.11.1 option 1 | volume $L^{2}$ | boundary $H^{1/2}$ | $H^{1/2}$ | complete controlled boundary and minimum-extension Schur metric |
| Section 5.11.1 option 2 | volume $H^{1}$ | boundary $L^{2}$ | $H^{1/2}$ | complete controlled boundary and zero-trace target assumption |
| Section 5.11.3 | volume $L^{2}$ | tangential boundary $H^{1}$ | tangential $H^{1}$ | complete controlled boundary and projected ambient gradients |

The state/test/control/observation topologies, residual formulation,
observation pairings, loss source, metric variable, and required typed
realizations are matched as one closed signature. A structurally valid graph
outside this table receives the lowerability diagnostic
`section_5_11_registered_signature`; it is never assembled using the nearest
target. In particular, changing state tracking, control-loss topology, or
search metric independently does not create an implicit cross-product.

The fixed-reconstruction, subdomain-tracking, $H^{1}_{0}$ state-tracking, and
general scalar Robin
registrations use the bounded component path. `SemanticResolver` turns a valid graph into stable-ID
lookup tables, and `DealiiScalarLoweringPlanner` invokes stored residual,
observation, loss, metric, constraint, and transformation handlers to build a
typed `ScalarLoweringPlan`. The remaining registrations deliberately retain
target-specific strategies. `DealiiCapabilityRegistryV1` is only a capability
ledger used by diagnostics; it is not described as a lowerer. Supporting an
individual kind therefore does not imply arbitrary recombination outside a
registered plan or target strategy.

The bounded scalar path projects that full plan into a
`ScalarResidualAssemblyPlan` before constructing its FE model. This executable
slice contains only the selected residual contributions, their resolved
`ScalarDataPlacement` requests, and the selected Robin boundary IDs. Volume
diffusion/reaction, transport, source, and control actions and Robin
bilinear/source actions are accumulated from those records independently;
unselected records do not trigger coefficient or forcing evaluation. The
plan-owned residual slice also selects the nonsymmetric direct solve when a
transport contribution is present. Objective, metric, constraint, and
specialized-target selection remain separate bounded services.

The scalar compiler also projects the plan into a `ScalarServicePlan`. Its
observation and loss records select which tracking and regularisation actions
are assembled, while its metric, constraint, and transformation records select
the cellwise Riesz service, optional box factory, and fixed-data reconstruction
port. The compiled scalar target consumes these records for its objective and
service wiring; changing a selected observation therefore leaves the residual
assembly unchanged, and changing the search metric does not change the
objective or adjoint source.

The general scalar Robin graph declares six coefficient/data spaces explicitly:
bounded-function tensor, vector, and scalar spaces on the full volume, plus a
bounded Robin coefficient space and a boundary $L^{2}$ source space on
`robin_boundary`. Semantic validation checks those typed ports, including
ownership of each Robin datum by its residual term. The scalar lowering plan
copies each resolved datum into a `ScalarDataPlacement` with semantic ID, role,
field kind, space, region, evaluation realization, and handler ID. General
scalar model construction consumes those placements before assembling volume
or boundary operators, and the compiled manifest preserves the same fields
alongside runtime provenance.

Its boundary policy is likewise resolved as a
`BoundaryRealisationSelection`: fixed Dirichlet is `dirichlet_boundary`, Robin
and natural transport outflow are `robin_boundary`, Neumann and transport
inflow are explicitly empty, and the selected realization is the outward
`(A grad(y) - b y)` conormal with an FE_Q state trace and face QGauss. The
selection is copied into `ScalarLoweringPlan`, checked by general scalar model
construction, and retained as a structured manifest record. The
`selected_policy` strings remain descriptive renderings only.

The Chapter 6 B2 selection extends this vocabulary with an explicit
ordinary-normal transport boundary form. It assembles the equivalent
state-dependent boundary correction and the diffusion-scaled Neumann control
coupling, and is retained in the typed semantic selection and compilation
manifest. The deal.II policy still permits an explicit total-conormal override
for diagnostics; that override does not alter the frozen B2 semantic policy or
runner default.

## Public semantic graph

The compatibility aggregate `include/nmopt/semantic/v1/problem_spec.hpp`
includes the focused deal.II-free headers `types.hpp`, `validation.hpp`,
`resolved_problem.hpp`, and `reference_specs.hpp`. The last contains
`make_scalar_diffusion_reaction_problem()`, the homogeneous comparison graph.
`make_fixed_dirichlet_scalar_diffusion_reaction_problem()` adds the first
declared physical-field transformation, while
`make_dirichlet_control_scalar_diffusion_reaction_problem()` adds an explicit
complete-boundary controlled-Dirichlet lifting.
`make_l2_dirichlet_laplace_control_problem()` instead declares the distinct
continuous transposition residual: an $L^{2}(\Omega)$ state, an
$L^{2}(\Gamma)$ control, $H^{2}(\Omega)\cap H^{1}_{0}(\Omega)$ tests, and the
boundary normal-test-derivative action. Its selected discrete policy lowers
only the conforming trace subspace through the equivalent variational lifting;
the semantic graph itself contains no continuous $H^{1}$ transformation.
`make_hhalf_dirichlet_laplace_control_problem()`,
`make_h1_tracking_hhalf_dirichlet_laplace_control_problem()`, and
`make_h1_dirichlet_laplace_control_problem()` declare the remaining scalar
Section 5.11 variants separately. They share a normalized variational
Laplacian and explicit Dirichlet lifting, but distinguish the control space,
state observation, control loss, and search metric. In particular, Section
5.11.1 option 2 has an $L^{2}(\Gamma)$ control loss and an
$H^{1/2}(\Gamma)$ search metric; those components are not aliases.
`make_partial_dirichlet_control_scalar_diffusion_reaction_problem()` adds a
separate fixed-data input, complete fixed/controlled partition, and selected
fixed-precedence interface policy. The
`make_subdomain_tracking_scalar_diffusion_reaction_problem()` adds a named
material-id observation region.
`make_h1_state_tracking_scalar_diffusion_reaction_problem()` instead selects
the full-domain energy observation and its $H^{1}_{0}$ pairing while retaining
the cellwise $L^{2}$ control metric.
`make_l2_state_tracking_continuous_control_problem()` retains the baseline
$L^{2}$ state loss and changes only the control realization to independent
homogeneous-Dirichlet conforming Lagrange coordinates. Hypercube meshes select
`FE_Q` with `QGauss`; simplex meshes select `FE_SimplexP` with
`QGaussSimplex`. The same state-element and quadrature choice is registered for
the Neumann-boundary-control target family. Mixed-reference-cell meshes and
simplex meshes sent to any other target receive lowerability diagnostics. The
separate
`make_l2_metric_h1_state_tracking_continuous_control_problem()` and
`make_hminus1_metric_h1_state_tracking_scalar_diffusion_reaction_problem()`
factories share one independent homogeneous-Dirichlet continuous-control graph
and change only its selected metric. `make_neumann_boundary_control_problem()`
declares a facewise Neumann control and state boundary trace.
`make_neumann_convection_subdomain_tracking_problem()` instead composes the
same natural control with conservative volume transport and a material-id
state restriction; its desired-state policy selects volume quadrature.
`make_weighted_boundary_trace_neumann_control_problem()` replaces that state
observation by a `weighted_boundary_trace` with an explicit immutable
`boundary_weight` data port; it does not alter the residual or control metric.
`make_pure_neumann_boundary_control_problem()` is its separate zero-reaction
mean-constraint variant. `make_h1_regularised_scalar_diffusion_reaction_problem()`
is the separate continuous-control objective variant.
`make_coefficient_identification_problem()` adds a physical diffusion
parameter. `make_general_scalar_elliptic_robin_problem()` recombines the
registered scalar residual terms and adds one natural Robin region. The
compatibility table below describes the homogeneous graph; the
named graph sections record the added features.

`make_point_sensor_scalar_diffusion_reaction_problem()` selects the C5.10
point-set observation. Its region owns immutable physical coordinates and its
finite-dimensional observation space records the point count. The selected
deal.II plan evaluates `FE_Q` shape functions at those physical coordinates,
assembles the observation transpose as $C_{h}^{\mathsf{T}}$, and records the
very-weak point-load policy in the manifest. No nearest-node, quadrature-point
coincidence, or general transposition fallback is implied.

`make_normal_flux_scalar_diffusion_reaction_problem()` selects the C5.8
strong-state normal-flux observation. Its boundary region owns the selected
boundary IDs, while the semantic policies declare the outward normal,
$Y=H^{2}(\Omega)\cap H^{1}_{0}(\Omega)$, the domain regularity assumption,
and the very-weak transposition. The deal.II plan evaluates
$\nabla y\mathbin\cdot n_{\mathrm{out}}$ using `FEFaceValues` at selected face
quadrature points and assembles the same face map's weighted transpose into
the tracking objective and adjoint source. An $H(\mathrm{div})$ or general
transposition lowerer is not registered.

The P5.3 factories keep the fixed and observed boundary selections explicit:
the point-sensor factory accepts `fixed_dirichlet_boundary_ids` (default
`{0}`), while the normal-flux factory accepts `normal_flux_boundary_ids`
(default `{1}`) and `fixed_dirichlet_boundary_ids` (default `{0, 1}`).
Semantic validation rejects duplicate IDs and requires every normal-flux ID to
belong to the fixed region. The registered point-sensor and normal-flux
lowerers additionally require every exterior mesh face to carry a declared
fixed ID; mesh-dependent failure is reported as the lowerability capability
`p53_complete_fixed_dirichlet_boundary`.

For P5.3 finite observations, the structured `state_observation_space`
manifest record uses the realized output dimension, not the state DoF count:
point sensors report one coefficient per declared coordinate, while
normal-flux observations report one coefficient per selected face-quadrature
sample. The value, JVP, realized-map, and manifest-space dimensions are one
contract.

The P5.3 selections are represented by a typed transposition payload rather
than by the policy description. Point-sensor and normal-flux graphs bind the
same scalar Dirichlet-Laplacian operator contract to different registered
very-weak discrete maps. The Chapter 5.11.2 graph also binds its $L^{2}$
continuous parent, $H^{1/2}$ conforming trace subspace, variational
equivalence, and discrete conormal policy. Semantic validation checks the
space, observation-source, regularity, equivalence references, and the
diffusion/reaction data ports defining $T=-\kappa\Delta+rI$; the normalized
Chapter 5.11.2 Laplacian leaves those ports empty. The
compiler copies the selection into the resolved decision and manifest and
rejects an unsupported operator, isomorphism, observation source, or discrete
realization.

The bounded C1/C2 and P5.2--P5.4 remediation registrations are
acceptance-verified:
weighted traces, normal fluxes, and transformed state observations obtain their
dimensions from realized map/space records; P5.3 requires complete fixed
Dirichlet coverage on the selected meshes; and P5.4 cross-products are limited
to the closed registration set. The scalar residual/data and objective/service
records are consumed by the bounded execution path, and the retained
specialized strategies remain explicit. These records describe the supported
slices only; broader unregistered recombination remains outside the accepted
scope.

Every enum-bearing semantic aggregate defaults to an explicit `unspecified`
value. `SemanticValidator` rejects that sentinel structurally, so default or
incrementally populated aggregates remain safe validation inputs rather than
containing indeterminate enum state. Components that expose a label port must
also provide a non-empty human-readable label.

The current dual-coefficient pairing rule is intentionally exact: a
`PairingSpec` must name the same semantic space ID on its primal and covector
ports, and every equation, observation, loss, and metric checks both ports.
This v1 rule does not encode a nontrivial dual-space relation. A future such
pairing must declare that relation explicitly before relaxing exact identity.

Reference variants are named feature deltas over shared declarations. Their
lookup, replacement, and removal helpers require exactly one matching stable
ID; declaration-vector order is not used as semantic identity.

`ObservationSpec::data_ids` is the explicit immutable-data input port for an
observation map. Existing restrictions and unweighted traces require it to be
empty. The registered weighted trace requires exactly one scalar Function
datum with role `observation_weight` in its output observation space.

```text
Region       one full volume region, named material-id volume subregions, and marked boundary ids
Space        scalar H1 state/test; selected L2-state/H2-test transposition parent; scalar L2 volume or facewise control; selected H1/2 or H1 trace control
Pairing      explicit coefficient pairings for state, test, control, and observations
Variable     one state and one control; the state may name a physical-field transformation
Data         forcing, desired state, fixed Dirichlet lifting, scalar/tensor diffusion,
             conservative/advective transport, reaction, Robin coefficient/source,
             observation weight, regularisation, and optional lower/upper bounds
Transformation optional fixed-Dirichlet reconstruction or complete/partial controlled-Dirichlet physical-state lifting
Residual     diffusion-reaction or normalized Laplacian, conservative transport, volume source/control, Neumann control, and the selected transposition Laplace/boundary-control actions; selected general scalar/Robin terms
Observation  volume/control restrictions; full-volume H1_0 state restriction; boundary trace; fixed-data weighted boundary trace
Loss         quadratic tracking and quadratic L2, H1/2, H1, or parameter regularisation
Metric       selected volume, facewise, trace, parameter, H1/2, H1, or Hminus1 metric
Constraint   optional cellwise or facewise L2 box where registered
Formulation  one-state/one-control reduced DTO
```

Concrete values are not semantic objects. `DealiiDataBindings<dim>` binds the
forcing and desired-state `Function` objects plus scalar coefficients only
after validation; its constant-diffusion binding is optional and is required
only by a graph that declares constant diffusion. A graph with the fixed-data
reconstruction also requires its optional `fixed_dirichlet_data` binding. The
selected data rule interpolates that
`Function` at the declared Dirichlet boundary DoFs to form the lifting; it is
not inferred from the forcing or target. If the graph declares its optional
box, the compiler also requires `CellwiseBoxDataBindings`: both bounds must
be scalar constants or both must be exact-layout `FE_DGQ(0)` coefficient
vectors.

A partial controlled-Dirichlet graph additionally requires the fixed
Dirichlet `Function` binding and provenance used to assemble
$`\ell_{0,h}`$. A graph with `weighted_boundary_trace` additionally requires
`DealiiWeightedTraceDataBindings<dim>`. Its scalar `Function` is independent
of the desired-state binding and carries separate provenance. Both Functions
are evaluated at the selected boundary face quadrature; missing, empty-
provenance, or multi-component weight bindings are lowerability diagnostics.

The P5.1 target additionally requires
`DealiiGeneralScalarDataBindings<dim>`. Rank-specific deal.II
`TensorFunction` ports bind tensor diffusion and the two vector transport
coefficients; one-component `Function` ports bind reaction, Robin coefficient,
and Robin source data. Every binding has separate provenance. A missing
binding, empty provenance, or multi-component scalar coefficient is a
lowerability diagnostic rather than an interpretation of component zero.

The C5.6 target instead requires the narrow
`DealiiConservativeTransportDataBindings<dim>` port. It binds only the vector
coefficient consumed by the conservative term and records its own provenance;
the target does not require unused advective or Robin data.

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

### $H^{1}_{0}$ state observation

`make_h1_state_tracking_scalar_diffusion_reaction_problem()` replaces only the
state observation topology and pairing. The selected full-domain observation
uses identity value, JVP, and coefficient transpose actions on the physical
`FE_Q` state. The quadratic loss lowers the declared pairing to the assembled
mass-plus-stiffness Riesz operator; the concrete target fuses that loss
derivative with the identity observation pullback. For desired state $z_{d}$,
the tracking loss is

```math
J_{\mathrm{state}}(y_{h})=
\frac{1}{2}\left(
\lVert y_{h}-z_{d}\rVert_{L^{2}}^{2}
+\lVert\nabla(y_{h}-z_{d})\rVert_{L^{2}}^{2}
\right).
```

The bound target `Function` supplies both value and gradient at the selected
volume quadrature, and the manifest records that rule. The first realization
requires homogeneous fixed Dirichlet data and the full volume. It retains the
cellwise $L^{2}$ control metric and regularisation; subdomain energy tracking
remains separate. The registered $H^{-1}$ metric is instead composed with the
continuous-control companion described below.

The H1-state registration also carries a typed target-data membership
assumption. It binds `desired_state` to `state_observation_space` and the
selected fixed/controlled Dirichlet boundary, and records H1 value/weak-gradient availability
and zero trace as a continuous-semantic `user_assumed` condition. The compiler
copies this selection into the resolved decision and manifest and renders its
status explicitly; it does not claim to prove the property for the bound
deal.II `Function`. The same assumption is inherited by the L2 and Hminus1
continuous-control comparison factories. A general H1 target remains a
separate semantic registration.

### $H^{-1}$ control metric

The P5.2 comparison factories expose the same independent control coordinates
$P_h\subset H^1_0(\Omega)$, energy observation, volume residual, $L^{2}$
control loss, state solve, and adjoint solve. The companion metric is the
control mass matrix $M_h$. The negative-norm metric is

```math
G_h=M_hK_h^{-1}M_h,
\qquad
G_h^{-1}=M_h^{-1}K_hM_h^{-1},
```

where $K_h$ is the homogeneous-Dirichlet control Laplacian. The boundary
condition makes $K_h$ coercive without a mean constraint. `Hminus1Metric`
performs the Laplacian and mass inverse actions with the selected
identity-preconditioned serial-CG tolerances. The manifest records the control
space, operator, boundary/no-mean choice, preconditioner, and tolerances. The
semantic graph also records a typed metric-realization selection: the control
primal/dual space, mass and Laplacian pairings, fixed boundary region, both
solve-policy IDs, operator and inverse sequences, and nullspace policy. The
compiler passes that resolved selection to the backend and rejects an absent
or mismatched selection; it does not infer the operation from the metric ID or
display formula.

This bounded target retains the existing positive $L^{2}$ control loss; it
verifies metric composition and does not claim the source catalogue's
unregularized infinite-dimensional problem as a separate executable target.
The focused contract proves exact metric action/inversion, identical reduced
objective and covector under $L^{2}$ and $H^{-1}$ selection, distinct search
directions, and quadratic state-recomputed Taylor remainders for both.

### Weighted boundary trace observation

`make_weighted_boundary_trace_neumann_control_problem()` is an ID-based delta
from the registered Neumann graph. Its state observation consumes the physical
state and exactly one immutable `observation_weight` datum and realizes

```math
O_h(y_h)=h\,\gamma y_h.
```

The graph records the model-author assumption
$h\in L^{\infty}(\Gamma_o)$ separately from the selected quadrature rule.

The selected map policy names the source and output spaces, observation
boundary, `boundary_weight` datum, face quadrature, pairing, and the common
transpose realization. The compiler passes this typed selection to the
backend, and the manifest records its stable ID and input/data references.

For both the baseline and weighted Neumann targets, the backend owns one
ordered evaluation table for the selected boundary-face quadrature points.
The realized observation dimension is the number of entries in that table,
and the same entries provide the value and JVP actions. Its recorded
quadrature weights define the output pairing used by the VJP, so the realized
map record, executable ports, and fused tracking objective share one sample
ordering rather than independently counting shape-function contributions.

The value and JVP multiply the state trace by the fixed scalar weight, while
the transpose action multiplies the boundary observation covector by the same
weight before the state-trace pullback. Consequently the quadratic tracking
term assembles $h^2\phi_i\phi_j$ and the target load assembles
$h z_d\phi_i$ at the same `QGauss` face points. The target norm remains
$\lVert z_d\rVert^2$, so the objective is exactly
$\frac12\lVert h\gamma y_h-z_d\rVert^2$. The Neumann residual, facewise
control layout, regularisation, and `l2_facewise` metric are unchanged. The
manifest records both Function bindings, their provenance, the shared face
quadrature rule, and the weighted-observation handler.

### General scalar elliptic and Robin composition

`make_general_scalar_elliptic_robin_problem()` replaces the legacy combined
constant diffusion-reaction contribution by independently handled terms. For
trial state $y$ and test function $v$, the assembled bilinear action is

```math
(A\nabla y,\nabla v)_{\Omega}
-(y,b\mathbin\cdot\nabla v)_{\Omega}
+(c\mathbin\cdot\nabla y,v)_{\Omega}
+(r y,v)_{\Omega}
+(h y,v)_{\Gamma_{R}}.
```

The residual also contains the existing negative volume forcing/control
loads and the negative Robin source load. The declared conormal convention is
$(A\nabla y-b y)\mathbin\cdot n$ with the outward normal. The first registered
boundary partition is one non-empty homogeneous fixed-Dirichlet region and
one disjoint non-empty Robin/transport-outflow region; together they must
cover every exterior face of the compiled mesh.

Uniform ellipticity, coefficient regularity, and coercivity are explicit
model-author assumptions recorded in the manifest, not properties inferred
by sampling the bound Functions. The state matrix is generally nonsymmetric.
The selected serial target factors it with `SparseDirectUMFPACK` and uses the
factorization's exact transpose solve for the adjoint. Residual JVP and VJP
continue to use the assembled matrix and its transpose directly.

### $H^{1}$ control regularisation and search metric

The first half of P2.3 adds a distinct
`quadratic_h1_control_regularisation` loss. It changes the control space to a
continuous scalar conforming Lagrange space on the state mesh and assembles

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

The Neumann residual policy carries a typed control-realization selection. It
binds the control variable, $L^{2}(\Gamma_{c})$ parent space, controlled
boundary, and search metric, and selects either facewise constants or a
continuous nodal trace. Both declarations are semantically valid. The current
deal.II registry lowers only facewise constants; selecting the continuous
trace produces the stable `continuous_neumann_control_lowerer` lowerability
diagnostic. A facewise box combined with the continuous selection is rejected
semantically rather than reinterpreted coefficientwise.

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

### Neumann control with transport and subdomain tracking

`make_neumann_convection_subdomain_tracking_problem()` retains the facewise
Neumann coupling but replaces boundary tracking by a material-id volume
restriction and adds the conservative weak term

```math
-(y_{h} b,\nabla v_{h})_{\Omega}.
```

The generic C5.6 target defaults to the conormal flux of this form. The B2
application selection explicitly replaces that default with the source's
ordinary-normal form and records the choice in the semantic selection and
manifest. Fixed and Neumann regions must form a complete disjoint exterior
partition, while the observed material ids must exist on the mesh. The graph
records the Lipschitz, coercivity, and $b\mathbin\cdot n\leq0$ assumptions
separately. The selected target uses `SparseDirectUMFPACK` for the
nonsymmetric state operator and its exact transpose for the adjoint. Its
independently integrated Q1 contract checks the transport value and
material-subdomain loss; the facewise $L^{2}$ metric and optional box are
unchanged.

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

### Controlled essential lifting

`make_dirichlet_control_scalar_diffusion_reaction_problem()` introduces the
second transformation input explicitly through
`TransformationSpec::control_variable_id`. Its first registered target in
`dealii_dirichlet_control.hpp` has

$$
y_{\mathrm{phys}}=P_{h}\widehat y_{h}+L_{D,h}u_{h}.
$$

It chooses $`\ell_{0,h}=0`$ and one shared nodal trace coefficient for every
state `FE_Q` DoF on the complete exterior controlled boundary. The boundary
$L^{2}$ mass assembles both the control regularisation and the
`l2_dirichlet_trace` metric. The residual and full-volume tracking objective
use $`y_{\mathrm{phys}}`$; their state and control pullbacks are
$`P_{h}^{\ast}`$ and $`L_{D,h}^{\ast}`$. State and adjoint solves use the
independent-coordinate system.

The Chapter 5.11.2 registration does not reuse the preceding graph's
continuous $H^{1}$ declaration. Its parent residual is

```math
\langle E_{\mathrm{tr}}(y,u;f),\psi\rangle
=(y,-\Delta\psi)_{\Omega}-(f,\psi)_{\Omega}
+(u,\partial_{n}\psi)_{\Gamma}.
```

The bounded deal.II strategy nevertheless reuses
`DirichletControlLiftingModel<dim>` with unit diffusion and zero reaction,
because its selected nodal trace space is contained in
$H^{1/2}(\Gamma)$. The manifest records the continuous parent spaces, absence
of a semantic lifting transformation, conforming trace subspace, variational-
transposition equivalence, domain assumption, normalized Laplacian, and
conormal policy. This registration does not accept a facewise or discontinuous
$L^{2}$ Dirichlet control and is not a general transposition lowerer.
`discrete_conormal_covector()` realizes that policy as

```math
q_{h}=L_{D,h}^{\ast}\left(A_{h}^{\mathsf T}P_{h}p_{h}
      -(M_{h}y_{h}-z_{d,h})\right),
```

so the compiled reduced covector is
$\beta M_{\Gamma,h}u_{h}-q_{h}$. It deliberately does not differentiate an
`FE_Q` adjoint pointwise on the boundary.

The remaining Section 5.11 registrations use the ordinary normalized
variational Laplacian with the same explicit complete-boundary lifting. For
Section 5.11.1, `TraceHhalfMetric` realizes
$G_{1/2,h}=A_{BB}-A_{BI}A_{II}^{-1}A_{IB}$ with
$A=M_{\Omega,h}+K_{\Omega,h}$: apply performs the minimum-extension interior
solve and inverse apply uses a full-volume $A$ solve with trace-supported
right-hand side. Option 1 also uses that action in the control loss. Option 2
keeps the loss action $M_{\Gamma,h}$ and uses $G_{1/2,h}$ only for direction
formation, while its state loss assembles value and gradient tracking.
Section 5.11.3 instead assembles
$M_{\Gamma,h}+K^{\tau}_{\Gamma,h}$ from projected tangential shape gradients.
All three produce distinct compiler IDs and structured metric/loss provenance.

The first P5.4 extension selects

$$
y_{\mathrm{phys}}=P_{h}\widehat y_{h}+\ell_{0,h}+L_{D,h}u_{h}
$$

on a complete fixed/controlled boundary partition. Fixed data own every
shared corner or interface DoF; the control is the relative-interior nodal
trace with zero endpoint extension. Both maps use explicit transpose
pullbacks, and the boundary $L^{2}$ mass is assembled only on the independent
control trace. Missing, overlapping, or incomplete partitions are compiler
diagnostics. Hanging-node trace relations, alternate corner policies,
fractional or tangential metrics on partial/nonmatching traces, and trace boxes
remain unregistered. The complete-boundary Section 5.11 Sobolev metrics are
the registered specializations described above.

The P5.4 policy payloads make the currently registered choices explicit. The
$H^{1/2}$ trace metric selects the volume mass-plus-stiffness Schur operator,
minimum-$H^{1}$ extension, full-volume inverse action, and control metric
solve. The boundary $H^{1}$ metric selects boundary mass plus projected
tangential stiffness with a positive-mass nullspace policy. The partial
Dirichlet graph selects fixed-data precedence, relative-interior nodal control
with zero endpoint extension, and unsupported hanging-node handling. These
typed selections are copied into the scalar lowering plan, resolved decision,
and compilation manifest; human-readable policy strings remain explanatory
provenance only.

## Validation and diagnostics

`SemanticValidator` validates semantic structure and declared policies.
`compiler::v1::DealiiCompiler::validate()` appends compiler-specific checks
to the same `ValidationReport`.

The supplied-OTD route has an explicit `SuppliedOTDDeclaration` on the
semantic `ProblemSpec`. It names state, adjoint, and control-stationarity
blocks; each block carries its semantic variable/residual spaces, trial and
test pairings, runtime layout identifiers, and discretisation/action
provenance. The declaration also carries the multiplier convention and
conversion, typed DTO comparison status, and action-level provenance. The
canonical scalar lowerer consumes this declaration and the compiled manifest
retains the typed declaration alongside its runtime layout projection. A DTO
graph with only its formulation kind and provenance changed is therefore
rejected before lowering.

The canonical supplied-OTD-to-quadratic-KKT adapter has a second, typed
`SuppliedOTDQuadraticKKTValidity` declaration on the supplied system. It must
explicitly identify the selected blocks, the affine residual and
point-independent JVP guarantee, canonical block signs, D/D-transpose and
full-KKT pairings, symmetry evidence, rank and kernel assumptions, and the
`lambda=-p` multiplier conversion. The adapter does not sample a callback to
infer constancy or symmetry. An undeclared or incompatible validity record is
a formulation contract error; a declared nonsymmetric record is retained as
GMRES-only and cannot advertise MINRES compatibility.

| Category | Produced by | Examples in v1 |
| --- | --- | --- |
| `structural` | `SemanticValidator` | incomplete nodes, missing labels or ports, incompatible pairings, orphan/duplicate term edges, and variable-space mismatches |
| `analytical_policy` | `SemanticValidator` | missing selected fixed/controlled-Dirichlet or cellwise-bound policy |
| `lowerability` | `DealiiCompiler`, `DealiiCapabilityRegistryV1`, and the selected lowering planner or target strategy | matrix-free execution, zero conforming Lagrange degree, unregistered or unhandled node kind, missing bound or fixed-lifting binding, incomplete controlled boundary |
| `formulation_capability` | `DealiiCompiler` | all-at-once formulation or a multi-block DTO request |

`CompilationResultT<Backend>` returns the report and contains only the selected
compiled product when it is valid. Unsupported choices are reported as
diagnostics; the compiler does not substitute another residual, observation,
metric, constraint, or formulation.

The explicit `CompilationProduct::quadratic_kkt` request is a separate
compiler product boundary. In the current v1 registration it is accepted only
for the unconstrained fixed-Dirichlet scalar DTO target and returns
`CompiledQuadraticKKTProblemT<Backend>`, never a reduced DTO or supplied-OTD
product. Its product owns the canonical primal, multiplier, adjoint,
stationarity, and equality layouts, action provenance, multiplier conversion,
rank/kernel declarations, exact assembled-block provenance,
symmetric-indefinite policy, and identity preconditioner baseline. A direct scalar
DTO, constrained target, or supplied
OTD request receives the stable `compiled_quadratic_kkt`
`formulation_capability` diagnostic before product construction; the compiler
does not infer a KKT product from DTO or supplied-OTD labels alone. The common
KKT contract requires complete one-to-one primal/stationarity and
multiplier/equality block pairings, plus declared evidence for the
$D^{\mathsf T}$ and full-transpose actions, before a symmetric product can
advertise MINRES compatibility. The serial solver packs range blocks through
those pairings rather than relying on equal flattened dimensions or incidental
layout order. Solver validation checks the GMRES restart basis only when GMRES
is selected; MINRES does not inherit an irrelevant GMRES-basis restriction.
Serial KKT results retain linear termination/report separately from
independently recomputed stationarity and equality residual fields.

The explicit `CompilationProduct::pdas` request is a second distinct
formulation boundary. It is registered only for the scalar distributed-control
DTO and canonical supplied-OTD targets with a cellwise `FE_DGQ(0)` box and the
positive-diagonal cellwise $L^{2}$ metric. It returns
`CompiledPDASProblemT<Backend>`, which owns the common KKT product, typed box
complementarity, metric owner, active-set policy, and inner KKT policy. Its
manifest records the selected bound source and control block, multiplier
representation, classification and residual tolerances, active-row rank and
kernel assumptions, inner solver, and explicit exclusions. The supplied-OTD
route is lowered through the canonical supplied-OTD KKT bridge before the box
service is constructed. Continuous controls, facewise or quadrature-point
bounds, unsupported constraint families, and non-cellwise multiplier metrics
receive stable `compiled_pdas_*` formulation diagnostics before construction.

Caller-provided compilation data use the same predictable boundary. Missing
or nonfinite scalar bindings, nonpositive diffusion/regularisation, missing
Function provenance labels, invalid solve policies, empty or unsupported
meshes, absent requested boundary ids, and nonfinite, reversed, or
layout-mismatched bound data are `lowerability` diagnostics. `ContractError`
is reserved for direct low-level constructor misuse and violated internal
invariants after validated lowering.

The compiler creates one resolved binding request after semantic resolution.
Each request record names the semantic datum, its declared space and region,
the runtime binding port and representation, and the selected evaluation
realization. Ordinary
forcing, desired-state, and fixed-Dirichlet `Function` bindings, together with
registered observation and scalar-coefficient Functions, pass through this
request before model construction. Their scalar component count is checked
against the resolved port; a multi-component Function is never silently
reduced to component zero. A null owned compilation session is likewise a
`lowerability` diagnostic (`compilation_session_presence`). The public
compiler therefore returns predictable diagnostics for caller-provided
binding/session errors, while direct low-level constructor misuse retains the
`ContractError` boundary above.

## Registered deal.II realization

The compatibility aggregate
`include/nmopt/compiler/v1/dealii_scalar_diffusion_reaction.hpp` exposes the
focused compiler headers: `compiled_problem.hpp`, `dealii_types.hpp`,
`dealii_capabilities.hpp`, `dealii_scalar_plan.hpp`, and
`dealii_compiler.hpp`. The compiler first resolves the semantic graph once.
For the bounded assembled scalar path it then builds a typed contribution plan;
otherwise it selects one of the registered private target strategies.
`dealii_fixed_dirichlet.hpp` owns `ScalarComponentModel`, the v1 physical-state
assembly that consumes the current fixed-reconstruction,
material-subdomain-tracking, $H^{1}_{0}$ state-tracking, and P5.1 general
scalar/Robin plans. The
private `dealii_neumann_boundary.hpp` target owns the distinct Neumann
residual, facewise control layout, unweighted or fixed-data weighted boundary
trace tracking, the C5.6 conservative-transport/material-subdomain
composition, facewise metric, facewise box realization, and pure-Neumann
mean-zero saddle realization. The
private `dealii_dirichlet_control.hpp` target owns the distinct controlled
essential reconstruction, the complete and selected partial fixed/controlled
trace policies, nodal trace layout, boundary mass metric, and both
transformation pullbacks. It also owns the conforming Galerkin realization of
the Chapter 5.11.2 transposition graph without changing that graph's continuous
spaces or residual. The private `dealii_continuous_control.hpp` target owns
continuous conforming Lagrange control for the registered $L^{2}$ and $H^{1}$ state losses,
$L^{2}$ or $H^{1}$ control loss, and $H^{-1}$, $H^{1}$, or $L^{2}$ search
metrics. The private `dealii_coefficient_identification.hpp` target owns the
cellwise positive
diffusion-parameter residual, its nonlinear first-order actions, parameter
metric, and parameter box. The registry otherwise supports only the listed
volume and Robin terms, full-domain control observation, full-domain or material-id
$L^{2}$ state restriction, full-domain $H^{1}_{0}$ state restriction,
quadratic losses, `L2` metric, optional cellwise box, and
fixed or controlled Dirichlet reconstruction. Its selected discrete policies are
assembled serial scalar conforming Lagrange state/test with degree at least
one and reduced DTO. Continuous-control and Neumann-boundary-control targets
select `FE_Q` on hypercubes or `FE_SimplexP` on simplices; the other registered
targets remain `FE_Q` on hypercubes. Their control realizations are
`FE_DGQ(0)` volume control on the state mesh, one facewise-constant Neumann
coefficient for every marked boundary face, one nodal Dirichlet trace
coefficient for every state DoF on the complete exterior controlled boundary,
one relative-interior nodal trace for the selected partial lifting, or
continuous `FE_Q` or `FE_SimplexP` volume control for a registered
continuous-control graph.
Pure Neumann is limited to
zero reaction and compatible forcing/control loads. Coefficient identification
instead selects a cellwise `FE_DGQ(0)` physical parameter with a strictly
positive box and reassembled state/adjoint matrices.

The P5.1 realization adds selected tensor/vector/scalar coefficient Functions,
one Robin boundary region, and a direct nonsymmetric state/transpose-adjoint
solve while retaining the same `FE_Q`/`FE_DGQ(0)` layouts and reduced DTO
ports.

The C5.6 realization reuses that conservative-transport convention with the
Neumann target's facewise control layout and material-id observation. It binds
only the conservative coefficient and likewise records direct state and
exact-transpose adjoint solves.

`CompiledProblemT<Backend>::executable_model()`, `metric()`, `constraint()`,
`reduced_hessian()`, and `make_reduced_dto()` expose only backend-neutral ports
and formulation services. The Hessian accessor is optional: it is present for
the exact linear-quadratic direct, assembled scalar, and continuous-control
targets and absent for target families without the declared capability. The
homogeneous private v0 target, v1 assembled target, mass metric, and box
constraint stay inside
`DealiiCompiler`. Each registered box is
constructed with the actual positive-diagonal `MassMetric` selected by the
same lowerer and retains that metric's opaque realization witness. Metric
display IDs remain descriptive provenance and cannot grant clipping support to
another operator. Solvers have no v1 branches.

The direct, assembled scalar, and continuous-control targets expose the exact
reduced-Hessian action only because their residual operators, state-tracking
operators, and control regularisation are assembled and linear-quadratic. The
action uses the reduced tangent-state and transpose-adjoint solves; this does
not extend the capability to generic nonlinear residuals or objectives. The
continuous-control native-output path reconstructs its constrained conforming
Lagrange coefficient vector before writing the state, adjoint, control,
forcing, and target fields on the shared volume mesh.

`DealiiCompilationSession<dim>` exclusively owns a static triangulation moved
into it and supplies the lifetime token retained by both the compiled problem
and every detached reduced service. This is the preferred long-lived API and
prevents caller mesh mutation. The source-compatible triangulation-reference
overload remains an explicitly borrowed, immutable-lifetime path recorded as
such in the manifest.

The supplied-OTD lowerer propagates that same owner token into the executable
system. The canonical supplied-OTD KKT bridge retains the token on the
detached quadratic product as well as in its callback-held system copy, so
residual, transpose, conversion, and serial-solve actions remain valid after
the compilation result and session handles are released. Callback-held
resources are destroyed before the owner token.

All iterative state and adjoint targets use the shared serial SPD solve
service with independently selected typed policies. Its result reports
convergence, iterations, requested tolerance, and achieved residual. The
pure-Neumann mean-zero saddle solve remains a separate typed
`SparseDirectUMFPACK` policy rather than being forced through the SPD path.
The P5.1 nonsymmetric target likewise records direct state and exact-transpose
adjoint solves; it does not misreport its operator as SPD.

Every successful compiled product also carries a versioned, structured
`CompilationManifest` (schema version 4). A `ResolvedCompilationDecision` is
created after semantic/lowerability checks and before model construction. It
captures the closed request, typed semantic inventory, regions, unresolved map
skeletons, and selected policy payloads. Model construction consumes the same
closed request and scalar plan recorded by that decision; it does not repeat
target selection from the raw graph. After construction,
`finalize_resolved_decision()` fills only executable dimensions and realized
service records from the decision, closed request, and constructed services.
Its interface does not receive the raw `ProblemSpec`, target enum, or a second
Dirichlet registration.

The finalized decision is therefore the single typed provenance record for the
selected target, formulation, regions, spaces, bindings, pairings, residuals,
observations, losses, transformations, policies, and realized services. The
manifest constructor is a pure descriptive projection of that finalized
decision, including its compatibility rendering; it does not rescan the raw
semantic graph or target registration. Compiler-selected requirement
provenance is rendered from typed IDs, kinds, statuses, scopes, and regions.
Only a `user_assumed` requirement retains its model-author declaration as
free-form semantic input; compiler capability and compatibility rendering do
not depend on selected-policy prose.

Binding records retain their semantic field shape, concrete runtime
representation, checked status, scalar value where applicable, and an exact
FNV-1a identity for vector bounds. Mesh records retain caller/owned lifetime
and provenance separately from a structural fingerprint covering topology,
coordinates, material IDs, and boundary IDs. Existing human-readable fields
are a rendered compatibility view, not a configuration or test-parsing
channel. The compiler selects a typed internal constraint realization
alongside the constructed service, then renders it into that view.

Bounded observations and physical-field transformations additionally publish
`CompiledRealizedMapRecord` entries. Each entry records its source and output
semantic spaces, actual dimensions, layout and ordering, pairing, value/JVP/
VJP provenance, and any transformation chain. Point sensors expose ordered
sensor values, normal flux and boundary traces expose ordered face-quadrature
values, and volume or energy observations expose the realized FE coefficient
map. The corresponding `CompiledRealizedSpaceRecord` entries are the sole
source for observation output dimensions; they do not fall back to the input
variable layout or the semantic space declaration.
Component-planned products additionally record the exact handler provenance
used to lower every contribution; tests assert those records rather than
inferring lowering from display text.
In particular, coefficient identification records
`l2_cellwise_parameter`, while volume and facewise controls record their own
distinct clipping realizations.

## Comparison guarantee

The `nmopt.dealii.canonical_volume_control` scenario contains two deliberately
different checks. A hand-integrated Q1/DGQ0 calculation on a $2\times2$ mesh
checks one residual coefficient and the objective independently of production
assembly. The compiled/direct comparison then constructs the direct model and
the homogeneous compiled product with identical mesh, data, coefficients, and
bounds and compares residual, objective, objective derivative, and reduced
derivative. That second check is a compiler wiring and packaging regression,
not an independent assembly oracle.

The same scenario also exercises the separate supplied-OTD product at a
nonzero manufactured point. It checks every residual block, centered finite
differences for all JVP blocks, the full residual JVP/VJP pairing with an
independent seed, and state/adjoint/stationarity actions reconstructed from
the DTO residual, objective derivative, and pullback under the declared
multiplier conversion. The supplied-OTD manifest records this as evidence
established by the named scenario; it does not imply that every compilation
performed the comparison.

Seven backend-neutral semantic scenarios validate every registered factory,
representative structural or policy failures, incomplete aggregates,
whole-graph closure, two-sided pairing compatibility, and order-independent
reference deltas. They also verify order-independent resolution and the exact
bounded scalar contribution plan, including rejection of a specialized graph.
The seventeen deal.II target scenarios named in the
[capability table](#registered-capabilities) exercise their selected target,
diagnostics, manifest, metric or constraint where applicable, forward and
transpose actions, and state-recomputed reduced derivative. Negative semantic
and compiler checks match diagnostic category, component ID, and capability.
Twelve further deal.II scenarios separately cover metric and continuous-control
component actions, projection-realization compatibility, compiler diagnostics,
owned-session lifetime, serial SPD solve reporting, and the backend's
native-size conversion boundary. The
backend-neutral contract suite separately exercises detached owned-service
lifetime under sanitizers. Exact assertions remain in the test sources rather
than being duplicated as a second inventory here.

## Exclusions

This v1 registration does not broaden the v0 executable mathematics. Beyond
the selected fixed-data reconstruction, complete-boundary nodal
Dirichlet-control lifting, the selected partial fixed/controlled nodal lifting,
material-id $L^{2}$ state tracking, full-domain $H^{1}_{0}$ state tracking,
marked-face Neumann control with unweighted or fixed-data weighted boundary
tracking, the selected C5.6 transport/subdomain composition, and the selected
general scalar/Robin target, it does not compile arbitrary geometric or
overlapping subdomain restrictions, FE target projection/interpolation, Robin
partitions beyond the registered homogeneous fixed/Robin split, alternate
partial controlled-Dirichlet interface policies, fractional or tangential
metrics on partial/nonmatching trace spaces, trace boxes,
continuous-control bounds, continuous or transformed coefficient parameters,
matrix-free execution, non-canonical all-at-once/OTD products, multiple
equations, or multiple optimisation variables. The canonical scalar
supplied-OTD product is registered separately as an all-at-once formulation
and exposed through `CompiledSuppliedOTDProblemT`; other all-at-once/OTD
choices still require their own semantic declaration, registered lowerer,
capability diagnostic, and value/JVP/VJP/reduced tests.
