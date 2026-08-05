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
the same coefficients. A graph that declares fixed-Dirichlet reconstruction
or material-subdomain state tracking selects a separate v1-only assembled
target; it never extends the v0 direct model. Both targets expose the same
generic compiled ports.

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
material-id observation region:

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
after validation. A graph with the reconstruction also requires its optional
`fixed_dirichlet_data` binding. The selected data rule interpolates that
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
registry otherwise supports only the listed volume terms, full-domain control
observation, full-domain or material-id state restriction, quadratic losses,
`L2` metric, optional cellwise box, and fixed-Dirichlet reconstruction. Its
selected discrete policy is assembled serial scalar `FE_Q` state/test with
degree at least one, `FE_DGQ(0)` control on the same mesh, fixed Dirichlet
boundary ids, and reduced DTO.

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
the compiled result.

`tests/semantic_v1_contract.cc` independently verifies the semantic validator
for the canonical and material-subdomain graphs, including missing or
region-mismatched target-data policies. The deal.II test compiles two material
masks on the same mesh and requires the same state solution and residual JVP,
while requiring a different tracking objective and state-objective derivative.

## Exclusions

This v1 registration does not broaden the v0 executable mathematics. Beyond
the selected fixed-data reconstruction and material-id state tracking, it does
not compile boundary observations, arbitrary geometric or overlapping
subdomain restrictions, FE target projection/interpolation, Neumann/Robin
terms, controlled Dirichlet liftings, continuous-control bounds,
non-$L^{2}$ metrics, matrix-free execution, all-at-once/OTD, multiple
equations, or multiple optimisation variables. Each requires its own semantic
declaration, registered lowerer, capability diagnostic, and value/JVP/VJP/
reduced tests.
