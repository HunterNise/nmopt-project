# Semantic problem model

## Describing the problem before constructing its numerical realization

Parts I and II assumed that the numerical objects already existed. We had residual
maps, objective derivatives, metrics, state and adjoint solves, KKT products, and
active-set subproblems, and we asked how those objects should behave once they reached
the formulation and solver layers.

The semantic path begins one stage earlier.

Its question is:

> What must be said about a PDE-constrained optimization problem before deciding how
> deal.II will realize it?

For the running example, the mathematical problem is familiar:

```math
\begin{aligned}
\min_{y,u}\quad
&
\frac12
\lVert y-y_{\mathrm d}\rVert_{L^{2}(\Omega)}^{2}
+
\frac{\beta}{2}
\lVert u\rVert_{L^{2}(\Omega)}^{2},
\\
\text{subject to}\quad
&
-\nabla\cdot(\kappa\nabla y)+cy
=
f+u
\quad\text{in }\Omega,
\\
&
y=0
\quad\text{on }\partial\Omega.
\end{aligned}
```

A finite-element implementation eventually needs a mesh, finite-element spaces,
quadrature, matrices, state and adjoint solvers, a metric realization, and perhaps a
box constraint. But those concrete objects are not the first things that distinguish
this problem from another one.

Before constructing them, we can already say that:

- there is a volume domain and a fixed Dirichlet boundary;
- the state and test functions have $H^{1}$ semantics;
- the control has $L^{2}$ semantics;
- the PDE residual contains diffusion–reaction, forcing, and volume-control terms;
- the objective observes the state and control through specified observation maps;
- one loss tracks a desired state and another regularizes the control;
- the control uses an $L^{2}$ optimization metric;
- the requested formulation is reduced DTO;
- selected analytical or discrete-realization assumptions must hold.

`ProblemSpec` is the current composition root for that description.

It is not an executable PDE class and it is not a symbolic algebra system. Its nodes
describe the **roles and relationships** that must survive until a numerical
realization is selected.

This chapter builds the running problem one layer at a time. Validation, resolution,
and the question of whether the current compiler can lower the resulting graph are
covered in Chapter 10.

## 1. Semantic description separates problem meaning from numerical ownership

The semantic graph is backend-neutral in a specific sense: it does not own deal.II
objects such as

```text
Triangulation
DoFHandler
FE_Q / FE_DGQ objects
SparseMatrix
Vector
SolverCG
Function objects used at runtime
```

Those belong to later numerical realization and runtime binding stages.

At the same time, the semantic graph is **not discretization-free mathematics**.
Some choices affect the meaning of the requested numerical problem and cannot be
safely invented later.

For example:

```text
volume control
vs
boundary control

full-domain tracking
vs
subdomain tracking

point observations
vs
L2 observations

fixed Dirichlet data
vs
controlled Dirichlet data

L2 control metric
vs
H1/2 trace metric
```

Those distinctions determine which weak terms, spaces, observations, pairings, or
realization requirements are needed.

The intended boundary is therefore:

```text
semantic problem
    what mathematical/numerical problem is being requested

runtime data
    which concrete functions/constants instantiate its data ports

compiler / realization
    how the requested problem is built with deal.II objects

solver
    how the resulting formulation product is iterated on
```

The separation is useful because the first layer can remain stable while the later
ones vary.

A parameter study may reuse one semantic graph with different forcing or target
functions. A scenario may use the same problem recipe with a different mesh
refinement or reduced solver. Conversely, changing a volume observation into a point
sensor changes the semantic problem itself, even if the same mesh happens to be used
afterward.

## 2. `ProblemSpec` – the composition root

The mathematical problem is one object, but `nmopt` does not encode it as one
monolithic C++ class. `ProblemSpec` is a composition root whose fields point to
smaller declarations:

```text
ProblemSpec
    id                         // stable graph key for the whole problem
    label                      // human-readable problem name

    regions                    // vector<RegionSpec>
    spaces                     // vector<SpaceSpec>
    pairings                   // vector<PairingSpec>
    variables                  // vector<VariableSpec>
    data                       // vector<DataSpec>
    transformations            // vector<TransformationSpec>
    residual_terms             // vector<ResidualTermSpec>
    equations                  // vector<EquationBlockSpec>
    observations               // vector<ObservationSpec>
    losses                     // vector<LossSpec>
    metrics                    // vector<MetricSpec>
    constraints                // vector<ConstraintSpec>
    requirement_policies       // vector<RequirementPolicySpec>

    formulation                // FormulationSpec
    supplied_otd_declaration   // optional<SuppliedOTDDeclaration>
```

The exact field definitions live in
[`include/nmopt/semantic/v1/types.hpp`](../../../include/nmopt/semantic/v1/types.hpp).
The current reference graphs are built in
[`include/nmopt/semantic/v1/problem_library.hpp`](../../../include/nmopt/semantic/v1/problem_library.hpp).

The useful mental model is not "one C++ struct per equation". It is a translation
through layers:

```text
mathematics
    what is the PDE-constrained problem?
        │
        ▼
semantic graph
    which regions, spaces, maps, losses, and formulation roles exist?
        │
        ▼
C++ construction
    fill those typed declarations and connect them by IDs
        │
        ▼
validation / resolution / lowering
    decide whether the graph is coherent and realizable
        │
        ▼
numerical realization
    construct FE spaces, operators, solves, metrics, ...
```

## 3. `RegionSpec` – where a mathematical object lives

A region answers the geometric question

> On which part of the physical problem is this space, term, observation, or policy
> defined?

For the running problem we need the full volume $\Omega$ and the fixed Dirichlet
boundary $`\Gamma_{D}`$:

```math
\Omega
\qquad\text{and}\qquad
\Gamma_{D}\subset\partial\Omega.
```

The exact schema is:

```text
RegionSpec
    id                  // stable graph key, e.g. "domain"
    label               // human-readable description
    kind                // RegionKind
    is_full_domain      // bool; true when this volume region is all of Omega
    boundary_ids        // vector<unsigned int>; selected boundary markers
    material_ids        // vector<unsigned int>; selected cell/material markers
    point_coordinates   // physical coordinates; used only by point_set
```

Current region kinds are:

| `RegionKind` | Mathematical / semantic picture |
| --- | --- |
| `volume` | a volume subset $\omega\subseteq\Omega$; `is_full_domain = true` represents $\omega=\Omega$ |
| `boundary` | a boundary subset $\Gamma\subseteq\partial\Omega$, selected by boundary IDs |
| `point_set` | a finite set $`X_{\mathrm s}=\{x_{1},\ldots,x_{m}\}\subset\overline\Omega`$ |

A `volume` or `boundary` region is selected through material or boundary IDs.
A `point_set` region stores physical point coordinates.

The semantic declaration stops at this geometric selection. It does not own a
`Triangulation` or decide which cells exist.

## 4. `SpaceSpec` – topology and semantic role

A space declaration answers two questions that should remain separate:

1. what mathematical topology does the quantity belong to?
2. what role does that space play in this problem?

For the baseline problem,

```math
\underbrace{Y}_{\text{state}}
=
H_{0}^{1}(\Omega),
\qquad
\underbrace{Z}_{\text{test}}
=
H_{0}^{1}(\Omega),
\qquad
\underbrace{U}_{\text{control}}
=
L^{2}(\Omega).
```

State and test may have the same topology without having the same semantic role.

The exact `SpaceSpec` schema is:

```text
SpaceSpec
    id          // stable graph key, e.g. "control_space"
    label       // human-readable description
    region_id   // RegionSpec id on which the space lives
    topology    // SpaceTopology
    role        // SpaceRole
    is_scalar   // bool; false for vector/tensor-valued spaces or data
    dimension   // positive only for explicitly finite outputs, e.g. m sensors
```

Current topologies are:

| `SpaceTopology` | Mathematical picture |
| --- | --- |
| `h1` | $H^{1}$-type functions; e.g. conforming state/test roles |
| `h2` | $H^{2}$-type regularity; used by strong/transposition auxiliary spaces |
| `hhalf` | $H^{1/2}$-type trace topology |
| `l2` | $L^{2}$-type fields or observations |
| `bounded_function` | bounded coefficient/data fields; deliberately distinct from an $L^{2}$ data assumption |

Current semantic roles are:

| `SpaceRole` | Typical mathematical role |
| --- | --- |
| `state` | state variable space $Y$ |
| `test` | residual test space $Z$ |
| `control` | optimization control space $U$ |
| `parameter` | parameter/inversion variable space |
| `observation` | codomain of an observation map $\mathcal O$ |
| `data` | semantic space occupied by prescribed coefficient/data fields |
| `auxiliary` | support space used by a formulation or realization, but not a primary variable |

For the control, the semantic statement is simply

```text
space_id = control_space
region = domain
topology = l2
role = control
```

Nothing here says `FE_DGQ(0)`, `FE_Q(1)`, or how many degrees of freedom the eventual
space has.

The exception is a genuinely finite semantic output such as an $m$-sensor
observation space, for which

$$
\mathcal O(y)\in\mathbb R^{m}
$$

and `dimension = m` is already part of the semantic meaning.

## 5. `PairingSpec` – how primal and dual roles are paired

A pairing describes a dual relationship between a primal space and the space in
which a covector is represented.

Its schema is deliberately small:

```text
PairingSpec
    id                  // stable graph key, e.g. "control_pairing"
    label               // human-readable description
    primal_space_id     // SpaceSpec used for primal values v
    covector_space_id   // SpaceSpec used to represent covectors r
```

Mathematically, it gives meaning to expressions of the form

$$
\langle r,v\rangle,
\qquad
r\in V^{\ast},
\quad
v\in V.
$$

This is not yet an optimization metric.

A metric additionally provides a Riesz map

$$
G:V\longrightarrow V^{\ast}.
$$

That distinction is why the baseline graph has pairings for the state, test,
control, and observation spaces even though only the control receives an optimization
metric.

## 6. `VariableSpec` and `DataSpec` – unknowns versus supplied inputs

The PDE contains both unknown fields and externally supplied quantities:

```math
-\nabla\cdot(
\underbrace{\kappa}_{\text{data}}
\nabla
\underbrace{y}_{\text{state}}
)
+
\underbrace{c}_{\text{data}}
\underbrace{y}_{\text{state}}
=
\underbrace{f}_{\text{data}}
+
\underbrace{u}_{\text{control}}.
```

`VariableSpec` is for optimization/state variables:

```text
VariableSpec
    id                          // stable graph key, e.g. "state"
    label                       // human-readable description
    role                        // VariableRole
    space_id                    // SpaceSpec id containing the variable
    physical_field_transform_id // TransformationSpec id; empty if already physical
```

Current variable roles are:

| `VariableRole` | Mathematical role |
| --- | --- |
| `state` | PDE state, e.g. $y$ |
| `control` | optimization control, e.g. $u$ |
| `parameter` | coefficient/parameter to be inferred or optimized, e.g. $m$ |

For the state, this means simply

```text
id = state
role = state
space_id = state_space
```

`DataSpec` is for prescribed inputs:

```text
DataSpec
    id        // stable graph key, e.g. "forcing"
    label     // human-readable description
    kind      // DataKind: shape/storage category of the supplied datum
    role      // DataRole: mathematical meaning in the problem
    space_id  // optional SpaceSpec id; current scalar-constant ports often leave it empty
```

Current data kinds are:

| `DataKind` | Mathematical / runtime picture |
| --- | --- |
| `function` | scalar field $x\mapsto f(x)$ |
| `vector_function` | vector field $x\mapsto \mathbf b(x)$ |
| `tensor_function` | tensor field $x\mapsto K(x)$ |
| `scalar_constant` | one scalar such as $\beta$ or a constant coefficient |
| `cellwise_bound` | bound values attached to cellwise control coefficients |
| `facewise_bound` | bound values attached to facewise control coefficients |

Current data roles are:

| `DataRole` | Typical mathematical meaning |
| --- | --- |
| `forcing` | source $f$ in the state equation |
| `desired_state` | target $`y_{\mathrm d}`$ or target observation data |
| `fixed_dirichlet_lifting` | prescribed boundary/lifting data $g$ used in physical-state reconstruction |
| `diffusion` | scalar or tensor diffusion coefficient such as $\kappa$ or $K$ |
| `conservative_transport` | conservative transport field entering a divergence/flux term |
| `advective_transport` | advection field $\mathbf b$ in a term such as $(\mathbf b\cdot\nabla y)v$ |
| `reaction` | reaction coefficient $c$ |
| `robin_coefficient` | Robin coefficient multiplying a boundary state term |
| `robin_source` | prescribed source in a Robin boundary condition |
| `natural_boundary_source` | prescribed Neumann/natural-boundary source |
| `regularisation_weight` | scalar weight such as $\beta$ in a quadratic penalty |
| `lower_bound` | lower admissible value $\ell$ |
| `upper_bound` | upper admissible value $r$ |
| `observation_weight` | fixed weight $w$ used inside a weighted observation |

The enum describes the semantic vocabulary. Whether a particular combination of
`kind`, `role`, space, and formulation is currently lowerable is a later compiler
question.

For the forcing, this means

```text
id = forcing
kind = function
role = forcing
space_id = state_test_space
```

The actual deal.II `Function` object is bound later. Changing its values does not
necessarily change the `ProblemSpec`.

## 7. `TransformationSpec` – from optimization coordinates to the physical field

A variable is not always the field seen by the PDE.

For fixed Dirichlet coordinates, Part I used

$$
\underbrace{\mathbf y_{\mathrm{phys}}}_{\text{physical state}}
=
\underbrace{P\mathbf z}_{\text{independent coordinates}}
+
\underbrace{\boldsymbol\ell}_{\text{fixed lift}}.
$$

The semantic graph represents this relationship with `TransformationSpec`:

```text
TransformationSpec
    id                   // stable graph key
    label                // human-readable description
    kind                 // TransformationKind
    input_variable_id    // independent-coordinate variable, e.g. z
    output_space_id      // semantic space of the reconstructed physical field
    fixed_data_id        // prescribed lifting/boundary datum g when required
    control_variable_id  // second primal input u for controlled Dirichlet lifting
```

Current transformation kinds are:

| `TransformationKind` | Mathematical picture |
| --- | --- |
| `fixed_dirichlet_reconstruction` | $`y_{\mathrm{phys}}=\mathcal T(z;g)`$: independent state coordinates plus prescribed boundary/lifting data |
| `dirichlet_control_lifting` | $`y_{\mathrm{phys}}=\mathcal T(z;u)`$: independent state coordinates plus a boundary-control input |

The two kinds correspond schematically to

```math
\underbrace{y_{\mathrm{phys}}}_{\text{output field}}
=
\underbrace{\mathcal T(z;g)}_{\text{fixed Dirichlet reconstruction}}
```

and

```math
\underbrace{y_{\mathrm{phys}}}_{\text{output field}}
=
\underbrace{\mathcal T(z;u)}_{\text{Dirichlet-control lifting}}.
```

A variable points to the transformation through
`physical_field_transform_id`. An empty ID means that the variable is already the
physical field.

The semantic graph records the relationship. The actual embedding $P$, lifting
vector, constrained DoFs, and FE interpolation belong to numerical realization.

## 8. `ResidualTermSpec` and `EquationBlockSpec` – from a weak equation to named pieces

For the running problem, the weak residual is

```math
\underbrace{
\int_{\Omega}
\kappa\nabla y\cdot\nabla v
\mathrm{d}x
+
\int_{\Omega}
c yv
\mathrm{d}x
}_{\texttt{diffusion\_reaction}}
-
\underbrace{
\int_{\Omega}
fv
\mathrm{d}x
}_{\texttt{volume\_source}}
-
\underbrace{
\int_{\Omega}
uv
\mathrm{d}x
}_{\texttt{volume\_control}}
=
0.
```

The mathematics is one residual equation. The semantic graph decomposes it into
registered term kinds and then groups those terms into one equation block.

`ResidualTermSpec` has the exact schema

```text
ResidualTermSpec
    id            // stable graph key, e.g. "volume_control"
    label         // human-readable description
    kind          // ResidualTermKind
    equation_id   // EquationBlockSpec this term contributes to
    variable_ids  // variables consumed by the term, e.g. {"state"} or {"control"}
    data_ids      // prescribed inputs consumed by the term
    region_id     // optional RegionSpec; empty for full-volume terms when implicit
```

while `EquationBlockSpec` is

```text
EquationBlockSpec
    id                 // stable graph key, e.g. "state_equation"
    label              // human-readable description
    test_space_id      // SpaceSpec for v
    test_pairing_id    // PairingSpec used for the residual/test duality
    residual_term_ids  // vector/list of ResidualTermSpec ids composing the equation
```

```cpp
ResidualTermSpec control_term;
control_term.kind = ResidualTermKind::volume_control;
control_term.equation_id = "state_equation";
control_term.variable_ids = {"control"};   // the u in -∫_Omega u v
```

Current residual-term kinds represent the following weak-form components. The sign
with which a term enters an equation is fixed by the selected residual convention and
lowerer.

| `ResidualTermKind` | Schematic weak-form / operator picture |
| --- | --- |
| `diffusion_reaction` | $`\int_{\Omega}\kappa\nabla y\cdot\nabla v\mathrm{d}x+\int_{\Omega}cyv\mathrm{d}x`$ |
| `tensor_diffusion` | $`\int_{\Omega}(K\nabla y)\cdot\nabla v\mathrm{d}x`$ |
| `conservative_transport` | conservative/divergence transport involving the flux $\mathbf b y$ |
| `advective_transport` | $`\int_{\Omega}(\mathbf b\cdot\nabla y)v\mathrm{d}x`$ |
| `reaction` | $`\int_{\Omega}cyv\mathrm{d}x`$ |
| `parameter_diffusion_reaction` | diffusion/reaction contribution whose coefficient depends on a declared parameter variable |
| `laplacian` | Laplace contribution, typically represented through $`\int_{\Omega}\nabla y\cdot\nabla v\mathrm{d}x`$ |
| `transposition_laplacian` | very-weak/transposed Laplace action in which derivatives are transferred to the test side |
| `dirichlet_transposition_control` | Dirichlet boundary-control contribution in a transposition/very-weak formulation |
| `volume_source` | source functional $`\int_{\Omega}fv\mathrm{d}x`$ |
| `volume_control` | distributed-control functional $`\int_{\Omega}uv\mathrm{d}x`$ |
| `neumann_control` | boundary-control functional $`\int_{\Gamma}uv\mathrm{d}s`$ |
| `robin_bilinear` | Robin bilinear term such as $`\int_{\Gamma}\rho yv\mathrm{d}s`$ |
| `robin_source` | Robin boundary source $`\int_{\Gamma}gv\mathrm{d}s`$ |
| `natural_boundary_source` | generic natural-boundary functional $`\int_{\Gamma}gv\mathrm{d}s`$ |

These names are not intended to form a universal symbolic PDE language. They name the
typed weak-form components currently represented by the semantic model.

The equation block then says, in effect,

```text
state_equation
    test with state_test_space
    pair through state_test_pairing
    sum:
        diffusion_reaction
        volume_source
        volume_control
```

The later compiler is responsible for realizing the corresponding integrals and
assembling or applying the resulting residual.

## 9. `ObservationSpec` – what quantity the objective actually sees

An observation is a map

$$
\mathcal O:
\text{input variable space}
\longrightarrow
\text{observation space}.
$$

Its exact schema is:

```text
ObservationSpec
    id                 // stable graph key, e.g. "state_observation"
    label              // human-readable description
    kind               // ObservationKind
    input_variable_id  // variable q to which O is applied
    region_id          // RegionSpec over which O acts
    output_space_id    // SpaceSpec containing O(q)
    output_pairing_id  // PairingSpec used on observation values/covectors
    data_ids           // immutable inputs consumed by O, e.g. a weight w
```

Current observation kinds are:

| `ObservationKind` | Mathematical picture |
| --- | --- |
| `volume_restriction` | $`\mathcal O(y)=y\rvert_{\omega}`$ in an $L^{2}$-type output space |
| `h1_state_restriction` | $`\mathcal O(y)=y\rvert_{\omega}`$, measured in an $H^{1}$-type output topology |
| `boundary_trace` | $`\mathcal O(y)=\gamma_{\Gamma}y`$ |
| `boundary_restriction` | $`\mathcal O(q)=q\rvert_{\Gamma}`$ for an already boundary-valued quantity $q$ |
| `weighted_boundary_trace` | $`\mathcal O(y)=w\gamma_{\Gamma}y`$ |
| `point_sensor` | $`\mathcal O(y)=[y(x_{1}),\ldots,y(x_{m})]^{\mathsf T}`$ |
| `normal_flux` | schematically $`\mathcal O(y)=\partial_{n} y`$ or the declared conormal variant |

The last row is intentionally schematic: the exact normal/conormal convention is
represented by requirement policies rather than being hidden inside the observation
name.

### 9.1 The baseline observation

Full-domain state tracking uses

$$
\underbrace{\mathcal O_{y}(y)}_{\texttt{volume\_restriction}}
=
y
\quad\text{on }\Omega.
$$

A minimal builder snippet is:

```cpp
ObservationSpec observation;
observation.kind = ObservationKind::volume_restriction;
observation.input_variable_id = "state";   // O_y takes y
observation.region_id = "domain";          // restriction is over Omega
```

The exact output-space and pairing IDs are part of the schema above.

### 9.2 Point sensors

For sensor locations

$$
X_{\mathrm s}
=
\{x_{1},\ldots,x_{m}\},
$$

the observation is simply

```math
\underbrace{\mathcal O_{\mathrm s}(y)}_{\texttt{point\_sensor}}
=
\begin{bmatrix}
y(x_{1})\\
\vdots\\
y(x_{m})
\end{bmatrix}
\in\mathbb R^{m}.
```

The semantic graph therefore needs:

```text
RegionSpec
    kind = point_set
    point_coordinates = {x_1, ..., x_m}

SpaceSpec
    role = observation
    dimension = m

ObservationSpec
    kind = point_sensor
    input_variable_id = state
    region_id = sensor region
    output_space_id = m-dimensional observation space
```

That is enough to understand what `point_sensor` means. How an FE field is evaluated
at those points, and which analytical/discrete policies make that operation legal,
belong to validation and realization.

## 10. `LossSpec` – how an observation contributes to the objective

`LossSpec` is deliberately downstream of `ObservationSpec`.

Its schema is:

```text
LossSpec
    id                     // stable graph key, e.g. "state_tracking"
    label                  // human-readable description
    kind                   // LossKind
    source_observation_id  // ObservationSpec producing the quantity being penalized
    data_id                // target or scalar weight, depending on LossKind
    pairing_id             // PairingSpec defining the quadratic pairing/norm
```

For the running problem, the objective can be read as a composition:

```math
J(y,u)
=
\underbrace{
\frac12
\left\lVert
\underbrace{\mathcal O_{y}(y)}_{\texttt{ObservationSpec}}
-
\underbrace{y_{\mathrm d}}_{\texttt{DataSpec}}
\right\rVert^{2}
}_{\texttt{quadratic\_tracking}}
+
\underbrace{
\frac{\beta}{2}
\left\lVert
\underbrace{\mathcal O_{u}(u)}_{\texttt{ObservationSpec}}
\right\rVert^{2}
}_{\texttt{quadratic\_control\_regularisation}}.
```

The division of responsibility is:

```text
ObservationSpec
    what quantity is seen?

LossSpec
    how does that quantity enter J?
```

Current loss kinds have the following mathematical interpretations:

| `LossKind` | Mathematical picture |
| --- | --- |
| `quadratic_tracking` | $`\frac12\lVert \mathcal O(q)-d\rVert_{P}^{2}`$ |
| `quadratic_control_regularisation` | $`\frac{\beta}{2}\lVert \mathcal O(u)\rVert_{P}^{2}`$ |
| `quadratic_hhalf_control_regularisation` | $`\frac{\beta}{2}\lVert u\rVert_{H^{1/2}}^{2}`$ |
| `quadratic_h1_control_regularisation` | $`\frac{\beta}{2}\lVert u\rVert_{H^{1}}^{2}`$ |
| `quadratic_parameter_regularisation` | schematically $`\frac{\beta}{2}\lVert \mathcal O(m)-m_{\mathrm{ref}}\rVert_{P}^{2}`$ for a declared parameter observation/reference |

Here $P$ denotes the pairing/topology selected by the surrounding semantic graph.
The last row is intentionally generic because the exact parameter observation and
data port are supplied by the problem builder.

These formulas explain the semantic intent of the names. They do not imply that every
mathematically conceivable combination is registered by the current compiler.

## 11. `MetricSpec` – optimization geometry, not objective regularization

A metric answers

> How is a derivative covector identified with a primal search direction?

Its schema is:

```text
MetricSpec
    id           // stable graph key, e.g. "control_l2_metric"
    label        // human-readable description
    kind         // MetricKind
    variable_id  // VariableSpec whose geometry is being declared
    pairing_id   // PairingSpec used by the Riesz map
```

Current metric kinds are:

| `MetricKind` | Mathematical picture |
| --- | --- |
| `l2` | $`(u,v)_{L^{2}}`$ and the corresponding $L^{2}$ Riesz map |
| `hhalf` | $H^{1/2}$ trace geometry; realization may require an extension/Schur construction |
| `h1` | $`(u,v)_{H^{1}}`$-type geometry, typically mass plus gradient terms |
| `hminus1` | $H^{-1}$ geometry; the current registered realization is based on a mass–Laplacian-inverse–mass construction |

Mathematically they select a Riesz map

$$
\underbrace{G}_{\texttt{MetricSpec}}
:
U_{h}
\longrightarrow
U_{h}^{\ast}.
$$

The running problem happens to use both an $L^{2}$ regularization loss and an $L^{2}$
optimization metric, but these are different semantic nodes:

```math
\underbrace{
\frac{\beta}{2}\lVert u\rVert_{L^{2}}^{2}
}_{\texttt{LossSpec}}
\qquad\neq\qquad
\underbrace{
G_{L^{2}}:U_{h}\to U_{h}^{\ast}
}_{\texttt{MetricSpec}}.
```

The former changes the objective. The latter changes gradient/search geometry.

## 12. `ConstraintSpec` – the admissible set

A box constraint has the mathematical form

$$
\underbrace{
\ell_{i}\leq u_{i}\leq r_{i}
}_{\texttt{ConstraintSpec}}.
$$

Its exact schema is:

```text
ConstraintSpec
    id                   // stable graph key, e.g. "control_box"
    label                // human-readable description
    kind                 // ConstraintKind
    variable_id          // constrained VariableSpec, usually "control"
    lower_bound_data_id  // DataSpec carrying ell
    upper_bound_data_id  // DataSpec carrying r
```

Current constraint kinds are:

| `ConstraintKind` | Mathematical / discrete picture |
| --- | --- |
| `cellwise_box` | coefficientwise bounds $`\ell_{K}\leq u_{K}\leq r_{K}`$ on a cellwise control representation |
| `facewise_box` | coefficientwise bounds $`\ell_{F}\leq u_{F}\leq r_{F}`$ on a facewise/boundary control representation |

The lower and upper values are separate `DataSpec` nodes, with data kinds
`cellwise_bound` or `facewise_bound`.

For a cellwise control box, the semantic mapping is

```text
kind = cellwise_box
variable_id = control
lower_bound_data_id = lower_bound
upper_bound_data_id = upper_bound
```

This node says what the admissible set is. It does not decide whether the eventual
algorithm uses projection or PDAS.

## 13. `RequirementPolicySpec` – assumptions and selected realizations

Some information is neither a variable nor a weak-form term.

Examples include:

```text
the state boundary is fixed Dirichlet
a transport problem is assumed coercive
a normal-flux observation uses a declared conormal convention
point evaluation uses a selected realization
an H1/2 metric uses a selected trace construction
```

Those statements are carried by `RequirementPolicySpec`.

The common schema is:

```text
RequirementPolicySpec
    id               // stable graph key
    subject_id       // id of the variable/equation/metric/etc. being constrained
    kind             // RequirementKind
    status           // RequirementStatus: how the statement is supplied
    scope            // RequirementScope: continuous, discrete, or both
    selected_policy  // human-readable statement of the selected assumption/policy
    region_id        // optional RegionSpec to which the policy applies

    typed_selection                      // optional BoundaryRealisationSelection
    typed_trace_selection                // optional TraceRealisationSelection
    typed_neumann_control_selection      // optional NeumannControlRealisationSelection
    typed_metric_selection               // optional Hminus1MetricRealisationSelection
    typed_transposition_selection        // optional TranspositionRealisationSelection
    typed_partial_boundary_selection     // optional PartialDirichletBoundarySelection
    typed_fractional_metric_selection    // optional FractionalTraceMetricRealisationSelection
    typed_boundary_h1_metric_selection   // optional BoundaryH1MetricRealisationSelection
    typed_h1_target_data_membership_selection // optional H1TargetDataMembershipSelection
```

Status says where the evidence comes from:

| `RequirementStatus` | Meaning |
| --- | --- |
| `provided` | the problem description supplies the required fact/data directly |
| `user_assumed` | an analytical assumption is declared by the author rather than proved by `nmopt` |
| `selected_discrete_realisation` | the problem selects one concrete registered discrete interpretation |

Scope says where the requirement matters:

| `RequirementScope` | Meaning |
| --- | --- |
| `continuous_semantics` | mathematical/analytical problem statement |
| `discrete_compilation` | compiler/discretization realization only |
| `both` | the requirement constrains both levels |

Current requirement kinds are:

| `RequirementKind` | Mathematical / numerical picture |
| --- | --- |
| `fixed_dirichlet` | fixed data such as $y=g$ on $`\Gamma_{D}`$ |
| `controlled_dirichlet` | controlled boundary value such as $y=u$ on $`\Gamma_{C}`$ |
| `mean_zero_multiplier` | normalization such as $`\int_{\Omega}\lambda\mathrm{d}x=0`$ |
| `boundary_trace` | existence/selection of a trace $`\gamma_{\Gamma}y`$ |
| `analytic_quadrature_evaluation` | prescribed data evaluated as $`f(x_{q})`$ at selected quadrature points |
| `discrete_cellwise_bounds` | coefficientwise bounds $`\ell_{K}\leq u_{K}\leq r_{K}`$ |
| `discrete_facewise_bounds` | coefficientwise boundary bounds $`\ell_{F}\leq u_{F}\leq r_{F}`$ |
| `uniform_ellipticity` | assumption such as $`\xi^{\mathsf T}K(x)\xi\geq\kappa_{0}\lVert\xi\rVert^{2}`$ |
| `coefficient_regularity` | declared regularity class for $K$, $c$, $\mathbf b$, or another coefficient |
| `coercivity` | assumption such as $a(v,v)\geq c\lVert v\rVert^{2}$ |
| `conormal_flux` | selected definition/sign of a normal or conormal flux $`q_{n}`$ |
| `boundary_partition` | decomposition such as $`\partial\Omega=\Gamma_{D}\cup\Gamma_{R}\cup\Gamma_{N}`$ |
| `transport_boundary_trace` | selected inflow/outflow trace convention for transport terms |
| `transposition_formulation` | very-weak form obtained by transferring derivatives to test functions |
| `domain_regularity` | regularity assumption on $\Omega$ needed by the selected formulation |
| `conforming_trace_subspace` | declaration that the chosen discrete trace lies in the required conforming subspace |
| `fractional_trace_realisation` | selected realization of an $H^{1/2}$ trace operator/metric |
| `tangential_gradient_realisation` | selected boundary gradient $`\nabla_{\tau}`$ used by an $H^{1}$ boundary metric |
| `target_data_membership` | membership such as $`y_{\mathrm d}\in H^{1}`$ together with any required trace condition |
| `metric_realisation` | selected discrete apply/inverse realization for the declared metric $G$ |

The table mixes analytical assumptions and discrete realization commitments on
purpose; `status` and `scope` distinguish those uses.

The important progression is:

```text
mathematical statement
    "assume the elliptic operator is coercive"
        │
        ▼
RequirementPolicySpec
    kind = coercivity
    status = user_assumed
    scope = continuous_semantics
```

or, for a discrete choice,

```text
mathematical/numerical requirement
    "evaluate point sensors with the registered FE realization"
        │
        ▼
RequirementPolicySpec
    status = selected_discrete_realisation
    scope = discrete_compilation
```

Some requirements also carry typed selection structs for choices that should not be
encoded only in free text. Their precise fields are listed in `types.hpp`; Chapter 10
explains how validation and resolution interpret them.

## 14. `FormulationSpec` – selecting the formulation roles

After declaring the PDE, objective, metric, and optional constraint, the graph must
say which components form the requested numerical product.

The current schema is:

```text
FormulationSpec
    id                   // stable graph key, e.g. "reduced_dto"
    kind                 // FormulationKind
    provenance           // FormulationProvenance
    state_variable_id    // state VariableSpec, e.g. "state"
    control_variable_id  // primary decision VariableSpec; may also be a parameter
    equation_id          // state/equality EquationBlockSpec
    metric_id            // MetricSpec on the decision variable
    constraint_id        // optional ConstraintSpec; empty for unconstrained problems
```

Current formulation kinds are:

| `FormulationKind` | Mathematical organization |
| --- | --- |
| `reduced_dto` | eliminate the state through $`y=S_{h}(u)`$ and optimize the reduced objective $`j_{h}(u)`$ |
| `all_at_once` | retain coupled first-order variables/blocks in one optimality system |

Provenance records how those equations arose:

| `FormulationProvenance` | Meaning |
| --- | --- |
| `dto` | differentiate the already-discretized objective/residual structure |
| `supplied_otd` | consume an application-supplied optimize-then-discretize first-order system |

For the baseline reduced problem,

```math
\underbrace{
u
\mapsto
y=S_{h}(u)
\mapsto
j_{h}(u)
}_{\texttt{FormulationKind::reduced\_dto}}
```

the semantic selection is essentially:

```cpp
specification.formulation.kind = FormulationKind::reduced_dto;
specification.formulation.provenance = FormulationProvenance::dto;
specification.formulation.state_variable_id = "state";      // y
specification.formulation.control_variable_id = "control";  // u
specification.formulation.equation_id = "state_equation";   // E(y,u)=0
specification.formulation.metric_id = "control_l2_metric";   // G on U_h
```

These assignments do not rebuild the PDE. They select previously declared nodes for
the roles required by the formulation.

## 15. `SuppliedOTDDeclaration` – declaring an application-supplied first-order system

A supplied OTD system starts from a different mathematical object.

For the linear running problem, the application may supply the complete first-order
residual

```math
F_{h}(y,p,u)
=
\begin{bmatrix}
\underbrace{Ay-f-Bu}_{\text{state block}}\\
\underbrace{A^{\mathsf T}p-M_{y} y+q}_{\text{adjoint block}}\\
\underbrace{B^{\mathsf T}p+\beta N_{u} u}_{\text{control-stationarity block}}
\end{bmatrix}.
```

The semantic graph must therefore say more than "all at once".

First, the formulation carries the provenance:

```cpp
specification.formulation.kind = FormulationKind::all_at_once;
specification.formulation.provenance =
  FormulationProvenance::supplied_otd;   // F_h is supplied by the application
```

Then `SuppliedOTDDeclaration` records the block structure:

```text
SuppliedOTDDeclaration
    id                           // stable graph key
    state_variable_id            // semantic state variable y
    adjoint_variable_id          // semantic adjoint variable p
    control_variable_id          // semantic control/decision variable u
    state_block                  // SuppliedOTDBlockSpec for E=0
    adjoint_block                // SuppliedOTDBlockSpec for state stationarity
    control_stationarity_block   // SuppliedOTDBlockSpec for control stationarity

    multiplier_convention        // SuppliedOTDMultiplierConvention
    multiplier_conversion        // SuppliedOTDMultiplierConversion

    value_action_provenance      // where supplied residual/value actions come from
    jvp_action_provenance        // provenance of the supplied linearization action
    vjp_action_provenance        // provenance of the supplied transpose action
    solve_provenance             // provenance of the coupled solve action

    comparison_status            // SuppliedOTDComparisonStatus
    comparison_evidence          // human-readable evidence for DTO/OTD comparison
```

Each block is a `SuppliedOTDBlockSpec`:

```text
SuppliedOTDBlockSpec
    id                         // stable graph key, e.g. "adjoint_block"
    label                      // human-readable description
    role                       // state / adjoint / control_stationarity
    variable_id                // variable occupying this block
    residual_id                // residual-equation identifier for this block
    variable_space_id          // semantic variable space
    residual_space_id          // semantic residual/covector space

    runtime_variable_space_id  // BlockLayout-facing runtime identifier
    runtime_residual_space_id  // BlockLayout-facing runtime residual identifier

    trial_pairing_id           // pairing used on the variable/trial side
    test_pairing_id            // pairing used on the residual/test side
    discretisation_provenance  // how this block was discretized
    action_provenance          // where its executable action comes from
```

The block-role options are:

| `SuppliedOTDBlockRole` | Mathematical block |
| --- | --- |
| `state` | state feasibility equation, e.g. $Ay-f-Bu=0$ |
| `adjoint` | state-stationarity / adjoint equation |
| `control_stationarity` | control first-order stationarity equation |

Multiplier convention says how the supplied multiplier relates to the framework
adjoint:

| `SuppliedOTDMultiplierConvention` | Meaning |
| --- | --- |
| `framework_adjoint` | supplied multiplier uses the framework adjoint sign convention |
| `negative_framework_adjoint` | supplied multiplier is the negative of the framework adjoint convention |

Conversion says what map is needed between the two representations:

| `SuppliedOTDMultiplierConversion` | Conversion |
| --- | --- |
| `identity` | $\lambda=p$ |
| `negate` | $\lambda=-p$ |

Comparison status records the relationship to the corresponding DTO construction:

| `SuppliedOTDComparisonStatus` | Meaning |
| --- | --- |
| `not_compared` | no DTO/OTD equivalence comparison is claimed |
| `different` | the supplied system is intentionally not equivalent to the compared DTO system |
| `equivalent_under_declared_conversion` | equivalence has been established after applying the declared multiplier conversion |

The transition to the executable layer is then:

```text
mathematical supplied system F_h
        │
        ▼
SuppliedOTDDeclaration
    block roles + spaces + pairings + provenance
        │
        ▼
compiler / binding
        │
        ▼
SuppliedOTDSystemT
    residual / JVP / VJP / solve
```

The application-facing declaration records **what system is supplied**. The later
runtime object provides **executable actions** for that system.

For an exact source example, see
`make_scalar_diffusion_reaction_supplied_otd_problem()` in
[`problem_library.hpp`](../../../include/nmopt/semantic/v1/problem_library.hpp).
The API/reference material should remain the authority for the complete field-by-field
construction.

## 16. Completed `ProblemSpec` – one graph, several mathematical layers

For the baseline problem, the completed graph can now be read almost directly from
the mathematics:

```text
Omega, Gamma_D
    ↓
RegionSpec

Y, Z, U and observation spaces
    ↓
SpaceSpec + PairingSpec

y, u
    ↓
VariableSpec

f, y_d, kappa, c, beta
    ↓
DataSpec

weak PDE terms
    ↓
ResidualTermSpec
    ↓
EquationBlockSpec

O_y(y), O_u(u)
    ↓
ObservationSpec

tracking + regularization
    ↓
LossSpec

G on U_h
    ↓
MetricSpec

optional ell <= u <= r
    ↓
ConstraintSpec

assumptions / realization commitments
    ↓
RequirementPolicySpec

reduced DTO or supplied OTD
    ↓
formulation declaration
```

The graph is more explicit than a three-line strong PDE statement because it preserves
the information needed by later compilation.

It is still much less concrete than a deal.II application. No mesh, finite-element
object, matrix, solver tolerance, line-search parameter, output directory, or
benchmark manifest belongs here.

## 17. `ProblemRecipeT` and `ScenarioT` – authoring above `ProblemSpec`

The application layer adds two useful levels above the semantic graph.

`ProblemRecipeT<Parameters>` is a typed factory:

```text
typed problem-family parameters
        │
        ▼
ProblemRecipeT
        │
        ▼
ProblemSpec
```

For the source Chapter 5 distributed-control family, one parameter selects whether the
semantic control graph represents the cellwise or continuous distributed-control
variant.

A short usage snippet is enough:

```cpp
auto recipe = make_scalar_distributed_recipe();

ScalarDistributedControlParameters parameters;
parameters.with_cellwise_box = true;  // request a semantic box constraint

ProblemSpec specification = recipe.build(parameters);
```

The recipe still owns no mesh, backend, optimizer, or experiment.

`ScenarioT` sits one level higher and groups choices from several layers:

```text
ScenarioT
├── problem parameters
├── compile options
├── solver options
└── experiment options
```

This is why a scenario is not another semantic specification.

A useful placement rule is:

```text
changes mathematical roles/relationships
    -> ProblemSpec or recipe parameter selecting a ProblemSpec

changes discretization/lowering
    -> compile options

changes optimization iteration
    -> solver options

changes experiment/reproduction organization
    -> experiment options
```

## 18. Reusing `ProblemSpec` components without creating one giant problem type

The reference builders illustrate another benefit of a compositional graph.

A subdomain-tracking problem can begin from the scalar baseline and change the
observation region and observation node.

A conservative-transport Neumann-control problem can reuse the natural-boundary
control structure while adding transport terms and analytical assumptions.

A fixed-Dirichlet problem can add a physical-field reconstruction without changing
the idea of state tracking.

These variations reuse semantic pieces because the changed concepts are explicit.

### 18.1 Reuse should still preserve problem meaning

Copying a baseline graph and modifying a few nodes is an implementation convenience,
not the semantics itself.

The resulting `ProblemSpec` must stand on its own.

A reader or validator should be able to understand the final graph without knowing
which factory it was copied from.

That is why IDs, roles, regions, pairings, and requirement selections remain explicit
after composition.

### 18.2 Reuse changes only the nodes whose mathematics changed

The point-sensor example from `ObservationSpec` illustrates the rule. Replacing
full-domain tracking by sensor tracking changes the observation region, observation
space, observation kind, target-data placement, and related requirements.

It does not require a different state equation merely because the objective now sees

$$
\mathcal O_{\mathrm s}(y)
=
[y(x_{1}),\ldots,y(x_{m})]^{\mathsf T}.
$$

This locality is the practical benefit of the graph: a mathematical change should
have a correspondingly local semantic footprint.

## 19. Limits of `ProblemSpec` – structured, but intentionally not universal

It is worth being explicit about what the current design does **not** promise.

`ProblemSpec` is not:

```text
a general symbolic weak-form language

a universal PDE ontology

a replacement for deal.II

a guarantee that every internally coherent graph is compiler-supported

the object used by every nmopt integration path
```

The current enums and node kinds describe the problem families the project has chosen
to represent.

New mathematical capabilities may require new semantic node kinds or new typed
realization selections.

This is preferable to pretending that an unsupported concept can be expressed by
stuffing it into a label string.

### 19.1 The external-application path does not need `ProblemSpec`

A mature PDE application can already own the correct mesh, spaces, operators, state
and adjoint solves, and output machinery.

In that case it can bind directly to the common numerical contracts developed in
Parts I and II.

The semantic/compiler path is one producer path:

```text
ProblemSpec
    ↓
validation / resolution
    ↓
compiler
    ↓
numerical formulation product
```

The native-integration path is another:

```text
existing numerical application
    ↓
binding/adapters
    ↓
numerical formulation product
```

The semantic graph is valuable when structured problem authoring and compiler reuse
are useful. It is not a compulsory wrapper around every external application.

## 20. `ProblemSpec` correspondence for the running problem

| Mathematical idea | Semantic declaration | Baseline example |
| --- | --- | --- |
| domain or boundary subset | `RegionSpec` | `domain`, `dirichlet_boundary` |
| function-space role | `SpaceSpec` | state $H^{1}$, control $L^{2}$ |
| dual relationship | `PairingSpec` | `state_test_pairing`, `control_pairing` |
| optimization unknown | `VariableSpec` | `state`, `control` |
| prescribed input | `DataSpec` | forcing, target, $\kappa$, $c$, $\beta$ |
| physical-field reconstruction | `TransformationSpec` | absent in homogeneous baseline |
| weak PDE contribution | `ResidualTermSpec` | diffusion–reaction, source, control |
| weak equation | `EquationBlockSpec` | `state_equation` |
| quantity seen by objective | `ObservationSpec` | full-domain state/control restriction |
| objective contribution | `LossSpec` | tracking, control regularization |
| optimization geometry | `MetricSpec` | control $L^{2}$ metric |
| admissible set | `ConstraintSpec` | optional cellwise box |
| analytical/discrete commitment | `RequirementPolicySpec` | fixed boundary, target quadrature |
| requested product | `FormulationSpec` | reduced DTO |
| supplied first-order blocks | `SuppliedOTDDeclaration` | absent in DTO baseline |

The table is a correspondence, not a construction order.

In code, the graph can be assembled in whichever order is convenient provided the
references are coherent when the graph is validated.

## 21. Source orientation for the semantic model

The central type declarations live in:

- [`include/nmopt/semantic/v1/types.hpp`](../../../include/nmopt/semantic/v1/types.hpp)

For this chapter, a useful reading order is:

```text
RegionSpec
SpaceSpec
PairingSpec
VariableSpec
DataSpec
TransformationSpec
ResidualTermSpec
EquationBlockSpec
ObservationSpec
LossSpec
MetricSpec
ConstraintSpec
RequirementPolicySpec
FormulationSpec
SuppliedOTDDeclaration
ProblemSpec
```

`FormulationSpec` selects either the reduced or all-at-once formulation kind; a
genuine supplied OTD request additionally carries `SuppliedOTDDeclaration`.

### 21.1 Read the problem library after the node types

Then read:

- [`include/nmopt/semantic/v1/problem_library.hpp`](../../../include/nmopt/semantic/v1/problem_library.hpp)

Start with

```text
make_scalar_diffusion_reaction_problem
```

because it is almost a literal construction of the running example in this chapter.

Then compare one or two deltas rather than reading the whole file linearly:

```text
fixed Dirichlet reconstruction
subdomain tracking
point sensors
Neumann boundary control
```

The contrast makes the purpose of transformations, observations, regions, and
requirements much easier to see.

### 21.2 The application recipe boundary

Read:

- [`include/nmopt/application/recipe.hpp`](../../../include/nmopt/application/recipe.hpp)
- [`include/nmopt/application/chapter5.hpp`](../../../include/nmopt/application/chapter5.hpp)

`ProblemRecipeT` is intentionally tiny.

`chapter5.hpp` shows how typed problem-family parameters choose among semantic
builders without taking ownership of compilation or execution.

### 21.3 The scenario boundary

Then read:

- [`include/nmopt/application/scenario.hpp`](../../../include/nmopt/application/scenario.hpp)

The class is small because the distinction matters more than the machinery: a
scenario groups problem, compile, solver, and experiment choices, but does not itself
compile or execute.

### 21.4 Semantic contract tests

Finally, the large semantic test file is useful as an executable catalogue:

- [`tests/semantic/semantic_v1_contract.cc`](../../../tests/semantic/semantic_v1_contract.cc)

It exercises the baseline graph and many structured variations, including continuous
controls, fixed-Dirichlet reconstruction, Dirichlet control, point sensors, normal
fluxes, and typed realization policies.

Chapter 10 uses those tests differently: not as a catalogue of graph
shapes, but as evidence for what validation and resolution are responsible for.

## 22. Continue with validation, resolution, and capabilities

At this point we have a graph.

We have **not** yet established that it is valid, fully resolved, or supported by the
current compiler.

Those are three different questions.

Chapter 10, **Validation, resolution, and capabilities**, distinguishes:

```text
structural validity
    do references, roles, regions, and pairings form a coherent graph?

analytical-policy validity
    are required assumptions / selections present and meaningful?

resolution
    have semantically relevant choices been closed before lowering?

compiler lowerability
    does the current registered compiler know how to realize this graph?

formulation capability
    can the requested numerical product be produced?
```

Keeping those questions separate is what prevents `ProblemSpec` from becoming a
hard-coded mirror of the current deal.II compiler.

## Read later

Useful existing documents are:

- [Describing and compiling a problem](../overview/semantic-compiler.md), for the
  shorter project-wide view of the semantic/compiler path.
- [Project architecture](../overview/project-architecture.md), for the relationship
  between semantic authoring and direct native integration.
- [Theoretical formalism](../../design/theoretical-formalism.md), for the mathematical
  spaces, residuals, observations, metrics, and formulation conventions represented
  by the graph.
- [Chapter 5 elliptic optimal-control guide](../../guides/chapter-5-elliptic-control.md),
  for the source problem families that motivated many of the semantic components.
- [v1 semantic/compiler capability](../../implementation/v1/semantic-compiler.md),
  for the exhaustive current capability ledger.

This chapter explains what a semantic problem **is**. The capability ledger remains
authoritative for what the current compiler can actually construct.
