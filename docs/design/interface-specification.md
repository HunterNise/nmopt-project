# Interface specification for composable PDE-constrained optimisation

## Status and normative words

This is an implementation-neutral specification.  It defines interfaces,
ports, and communication protocols; it does not prescribe C++ classes,
storage layouts, or deal.II calls.

The words **MUST**, **MUST NOT**, **SHOULD**, and **MAY** are normative.

The framework describes declared formulations.  It MUST NOT claim that a
declared formulation is well posed, that a trace exists, or that a selected
discretisation is convergent.  The model author supplies those analysis
assumptions through requirements and policies.

## 1. System boundary

The semantic model is

$$
\min_{x\in X_{\mathrm{ad}}}J(x;d)
\qquad\text{subject to}\qquad
E_{a}(x;d)=0\in Z_{a}^{\ast},
\qquad a\in\mathcal A.
$$

Here

$$
X=\prod_{i\in\mathcal I}X_{i},
\qquad
E_{a}:X\times D\to Z_{a}^{\ast}.
$$

A variable block $`x_{i}`$ can be a state, control, parameter, flux, initial
state, or any auxiliary field.  Data $d\in D$ are fixed during an
optimisation solve.

The semantic model is a directed composition graph:

$$
\text{variables and data}
\longrightarrow
\text{transformations}
\longrightarrow
\begin{cases}
\text{equation blocks},\\
\text{observations}
\end{cases}
\longrightarrow
\text{losses, metrics, constraints}.
$$

The compiler lowers this graph to discrete operator actions.  A formulation
builder constructs reduced or all-at-once solver inputs from those actions.

The framework MUST NOT introduce a problem class for a PDE family, control
placement, observation type, or solver method.

A strong differential expression MAY be retained as documentation.  It MUST
NOT be sufficient input to compilation: compilation begins only after trial
and test spaces, boundary interpretation, pairings, and residual terms have
been declared.

## 2. Common rules for every component

Every semantic component MUST declare:

| Declaration | Meaning |
|---|---|
| Identity | Stable identifier and human-readable label |
| Input ports | Named references to spaces, regions, variables, data, or maps |
| Output port | The type, space, and duality status of its result |
| Requirements | Required trace, product, pairing, nullspace, temporal, or discrete policy |
| Derivative contract | Required derivative actions when the component depends on variables |
| Compilation contract | What must be chosen to obtain a discrete realization |

A component MAY reference another component only through a declared port.
It MUST NOT inspect an optimiser, a global matrix, or an unrelated PDE term.

A port connection is valid only if all of the following hold:

1. The source and target spaces agree, or a declared transformation connects
   them.
2. Field shape, domain dimension, region kind, and component count agree.
3. A dual value is paired with its declared primal space.
4. A trace, product, point evaluation, or fractional operation has a declared
   capability and policy.
5. The source component provides every derivative action required downstream.

No implicit conversion from a derivative in $U^{\ast}$ to a vector in $U$ is
allowed.  Such a conversion requires a metric.

## 3. Primitive semantic components

### 3.1 Region

A region identifies a geometric set.  It MUST declare:

| Field | Requirement |
|---|---|
| Kind | volume, boundary, interface, point set, time interval, or space-time set |
| Parent | containing region or geometry |
| Dimension and measure kind | required for integration and traces |
| Identity | stable name usable by all components |
| Relation | optional disjoint partition, subset, orientation, or pairing relation |

A region owns no PDE term and no finite-element policy.

For a point-set region, the declaration additionally owns an ordered finite
set of immutable physical coordinates. A point-sensor observation maps the
declared state to a finite-dimensional observation space with one coefficient
per coordinate. The coordinate evaluation rule and the transpose/very-weak
adjoint policy remain explicit requirements of the observation and compiler;
they are not inferred from mesh nodes or quadrature points.

A boundary partition is a region declaration, for example

$$
\Gamma=\Gamma_{D}\mathbin{\dot\cup}\Gamma_{N}\mathbin{\dot\cup}\Gamma_{R}.
$$

The compiler MUST reject a request for boundary integration or trace evaluation
on a region that is not declared as a compatible boundary set.

When a selected realization depends on boundary roles, the semantic graph MUST
carry a typed boundary selection in addition to its human-readable policy
rendering. The selection identifies fixed-Dirichlet, Robin, Neumann,
transport-inflow, and transport-outflow regions, together with conormal form,
normal orientation, state-trace realization, and face-quadrature realization.
An empty role is represented explicitly by an empty region selection. The
compiler MUST consume this typed selection; it MUST NOT parse policy prose to
choose the executed boundary interpretation.

A Neumann-control residual likewise MUST carry a typed discrete-control
selection. The selection identifies the control variable, its mathematical
space, the controlled boundary, and the selected metric, and distinguishes a
facewise-constant realization from a continuous nodal-trace realization. Both
are discrete subspaces of the declared $L^{2}(\Gamma_{c})$ control space; the
topology of the semantic parent space does not select one implicitly. A
facewise coefficient box applies only to the facewise-constant realization.
The compiler MUST diagnose a selected realization for which it has no
registered coupling, metric, or constraint lowerer.

### 3.2 Space and pairing

A space describes a mathematical source or target.  It MUST declare:

| Field | Requirement |
|---|---|
| Base region | the region on which functions, traces, or values live |
| Field shape | scalar, vector, tensor, or product-block shape |
| Topology | for example $L^{2}$, $H^{1}$, $H(\mathrm{div})$, a trace space, a bounded-function coefficient space, or a Bochner space |
| Role | trial, test, observation, data, control, parameter, or auxiliary |
| Dual pairing | explicit primal-dual pairing or reference to one |
| Capabilities | declared trace, product, derivative, restriction, and temporal-trace capabilities |

A space descriptor is semantic.  It MUST NOT contain a finite element, DoF
handler, matrix, or vector.

A pairing is a separate declaration

$$
\langle\cdot,\cdot\rangle_{S^{\ast},S}:S^{\ast}\times S\to\mathbb R.
$$

A compiler MUST lower every pairing used by a residual, loss, metric, or
transpose action.  It MUST NOT substitute a mass matrix unless that is the
chosen realization of the pairing.

### 3.3 Variable block

A variable block MUST declare:

| Field | Requirement |
|---|---|
| Identifier and role | state, control, parameter, flux, or auxiliary |
| Space | one semantic primal space $`X_{i}`$ |
| Physical-field relation | identity or a declared transformation output |
| Admissible-set reference | optional constraint component |
| Differentiability status | differentiable, nonsmooth, or externally supplied derivative policy |

A variable block owns neither an equation nor an objective.  The same variable
MAY feed multiple residual terms, transformations, observations, and losses.

### 3.4 Data

Data are immutable within one solve.  A data component MUST declare its space,
region, field shape, and value source.  It MAY feed terms, transformations,
observations, losses, or constraints.

The v1 general scalar elliptic/Robin slice uses an explicit bounded-function
topology for tensor, vector, and scalar coefficient fields. Its tensor
diffusion, conservative transport, advective transport, and reaction data
live on the full volume region; its Robin coefficient and source live on the
declared Robin boundary. The Robin source uses boundary $L^{2}$ topology for
the selected trace pairing. These semantic placements are distinct from the
backend's quadrature or Function objects and are carried into the resolved
lowering plan.

A data component MUST NOT expose a derivative block.  To estimate a datum, the
author MUST replace it with a variable block and add all required residual,
loss, metric, and constraint connections.

For a tracking loss declared on an $H^{1}_{0}$ observation space, the graph
MUST carry a typed target-data membership policy. It MUST identify the
`desired_state` datum, the observation space, the selected fixed-Dirichlet
boundary, the availability of values and weak gradients, and zero-trace
membership on that boundary. The current registration marks this condition as
a continuous-semantic `user_assumed` model-author assertion. Validation and
the compiler record the assertion; they MUST NOT present evaluation of an
arbitrary analytic Function as a runtime proof of its trace.

### 3.5 General map

A general map is the universal compositional primitive:

$$
T:X_{1}\times\cdots\times X_{n}\to Y.
$$

It MUST declare source ports and target space.  If variable-dependent, it
MUST provide:

$$
T(x),
\qquad
T'(x)\delta x,
\qquad
T'(x)^{\ast}\eta,
\quad \eta\in Y^{\ast}.
$$

The last action returns one covector in each variable source dual:

$$
T'(x)^{\ast}\eta\in X_{1}^{\ast}\times\cdots\times X_{n}^{\ast}.
$$

For a linear map, the implementation MAY provide a constant action.  It MUST
still provide the same semantic ports.

### 3.6 Transformation

A transformation is a general map whose output is used as a physical field
before residual or observation evaluation.  It MUST declare the physical field
it reconstructs and the stage at which it applies.

Examples are

$$
y_{\mathrm{phys}}=\widehat y+\ell_{D}(g_{D}),
$$

$$
y_{\mathrm{phys}}=\widehat y+\ell_{0}+L_{D}u,
$$

$$
m=\exp(q).
$$

A transformation MUST provide value, JVP, and transpose-JVP.  Its transpose
action is how the chain rule enters a reduced derivative.

A transformation MUST NOT decide whether its output is a state equation,
tracking quantity, or optimisation variable.  Those are downstream
connections.

### 3.7 Residual term

A residual term contributes to exactly one equation block.  For equation test
space $`Z_{a}`$, it represents

$$
E_{a,t}:X_{a,t}\to Z_{a}^{\ast},
$$

through its tested action

$$
e_{a,t}(x;z)
=\langle E_{a,t}(x),z\rangle_{Z_{a}^{\ast},Z_{a}}.
$$

It MUST provide:

$$
e_{a,t}(x;z),
$$

$$
D_{x}e_{a,t}(x;z)\delta x,
$$

$$
D_{x}E_{a,t}(x)^{\ast}p,
\qquad p\in Z_{a}.
$$

The transpose action returns covectors for every variable input of the term.
A term MAY depend on data and transformations.  It MUST NOT choose an
objective, metric, constraint, state solver, or adjoint solver.

Typical term instances are diffusion, transport, reaction, source, boundary
flux, Robin, mixed coupling, time derivative, stabilization, and
control-to-residual coupling.

### 3.8 Equation block

An equation block owns a test space $`Z_{a}`$ and an ordered sum of residual
terms:

$$
E_{a}(x;d)=\sum_{t\in\mathcal T_{a}}E_{a,t}(x;d)\in Z_{a}^{\ast}.
$$

It MUST expose accumulated actions:

$$
E_{a}(x;d),
\qquad
E_{a}'(x;d)\delta x,
\qquad
E_{a}'(x;d)^{\ast}p_{a}.
$$

An equation block MAY be nonlinear, mixed, time dependent, or
Petrov–Galerkin after compilation.  It MUST NOT be named or selected by a
PDE-family enum.

### 3.9 Observation

An observation is a general map

$$
O_{k}:X_{\mathrm{phys}}\times D\to Q_{k}.
$$

It MUST provide value, JVP, and transpose-JVP.  It MUST NOT contain a loss.
Its semantic declaration MUST name every immutable datum consumed by the map
on explicit data-input ports; observation data MUST NOT be inferred from a
downstream loss or residual.

Typical instances are restriction to a subdomain, boundary trace, normal flux,
sensor array, time sampling, and actuator output.

For a sensor array, the selected discrete map may evaluate finite-element shape
functions at the declared physical coordinates and assemble its transpose as a
finite-dimensional point load. This is a discrete realization policy, not an
implicit nearest-node approximation or a claim that the continuous Dirac map
is available on every state space.

For a normal-flux observation, the declaration MUST also identify the outward
normal convention, the boundary region, the state regularity or
$H(\mathrm{div})$ capability that makes the flux meaningful, and the selected
transposition policy for its lower-regularity adjoint. A face-quadrature
`FE_Q` normal derivative and its assembled transpose are one possible
realization; they MUST NOT be inferred from an ordinary boundary trace.

The selected P5.3 transposition realization MUST be structured. It names the
subject equation, strong test space $Y$, operator range, isomorphism,
residual codomain, multiplier space, observation output and source space,
domain-regularity policy, and discrete realization. Point sensors and
normal-flux observations select their respective `FE_Q` very-weak maps. The
scalar diffusion-reaction realization additionally names the diffusion and
reaction data ports that define $T=-\kappa\Delta+rI$; a normalized Laplacian
realization leaves both coefficient ports empty. The
Chapter 5.11.2 realization additionally names the continuous parent space,
conforming trace space, variational-equivalence policy, and discrete conormal
policy. These fields select the executable map; display prose only explains
the selected contract.

A weighted boundary trace is likewise a distinct map realization. Its semantic
policy MUST name the source and output spaces, boundary region, immutable weight
datum, face quadrature rule, pairing rule, and transpose realization. The
value, JVP, and transpose-JVP MUST use that same selected face rule; a prose
label or an untyped weight datum is not sufficient to select the map.

### 3.10 Loss and objective

A loss is a scalar map

$$
\Phi_{k}:Q_{k}\times D\to\mathbb R.
$$

It MUST provide

$$
\Phi_{k}(q;d),
\qquad
D_{q}\Phi_{k}(q;d)\in Q_{k}^{\ast}.
$$

An objective is a sum of loss compositions:

$$
J(x;d)=\sum_{k\in\mathcal K}\Phi_{k}(O_{k}(x;d);d).
$$

The objective component MUST assemble its derivative by the chain rule:

$$
D_{x}J(x;d)
=\sum_{k} O_{k}'(x;d)^{\ast}D_{q}\Phi_{k}(O_{k}(x;d);d).
$$

A tracking norm is a loss and its pairing.  It MUST NOT be represented by an
unqualified string such as “L2 tracking.”

### 3.11 Metric

A metric provides an explicit primal-dual identification for an algorithm.
It MUST declare primal search space $P$, dual space $P^{\ast}$, and actions

$$
G:P\to P^{\ast},
\qquad
G^{-1}:P^{\ast}\to P,
$$

where inverse-apply MAY be unavailable if no selected solver needs it.

If the reduced derivative is naturally in $U^{\ast}$ but the search space is
$P$ with an injection $\iota:P\hookrightarrow U$, the metric protocol MUST
also provide or reference

$$
\iota^{\ast}:U^{\ast}\to P^{\ast},
\qquad
\nabla_{G}j=G^{-1}\iota^{\ast}j'.
$$

A metric MUST NOT alter a residual or objective.  An $H^{1}$ regularisation is
a loss; an $H^{1}$ search geometry is a metric.

The first selected $H^{-1}$ realization makes the discrete policy explicit.
Its search space is the independent homogeneous-Dirichlet coefficient space
$P_h=\operatorname{span}\{\phi_i\}\subset H^1_0(\Omega)$. With the control
mass and Dirichlet-Laplacian matrices

$$
(M_h)_{ij}=(\phi_j,\phi_i)_{L^2},
\qquad
(K_h)_{ij}=(\nabla\phi_j,\nabla\phi_i)_{L^2},
$$

the metric is the pulled-back negative norm

$$
G_h=M_hK_h^{-1}M_h,
\qquad
G_h^{-1}=M_h^{-1}K_hM_h^{-1}.
$$

The fixed-Dirichlet policy removes the constant nullspace, so this realization
has no mean constraint. Both elliptic and mass inverse actions MUST use and
record the selected metric-solve tolerances and preconditioner. Selecting this
metric changes only direction formation; observations, losses, residuals, and
adjoint equations remain unchanged.

The H^{-1} choice MUST be carried by a typed metric-realization policy that
binds the primal and dual spaces, mass and Laplacian pairings, fixed-Dirichlet
boundary region, both solve policies, operator sequence, inverse sequence, and
nullspace policy. The continuous control's fixed boundary policy MUST select
the same complete exterior boundary as the state realization. Human-readable
metric or formula text may explain the choice, but MUST NOT substitute for
these structured references.

The selected P5.4 metric and boundary policies MUST likewise be typed. A
fractional trace metric names its control and volume spaces, trace inclusion,
volume operator, minimum-extension apply action, full-volume inverse action,
and solve policy. The boundary $H^{1}$ metric names its control space and
boundary region, selects boundary mass plus tangential stiffness, identifies
the projected ambient tangential gradient, and declares its positive-mass
nullspace policy. A partial fixed/controlled Dirichlet partition names the
state, lifting transformation, fixed and controlled regions, complete and
disjoint-partition requirements, fixed-data interface ownership, relative-
interior nodal trace rule, and hanging-node status. An alternate realization
is a different semantic selection, not a relabeling of the descriptive policy
string.

### 3.12 Constraint

A constraint acts on one or more variable blocks.  It MUST declare its source
space and supported operations.  It MAY provide any of:

$$
\mathrm{is\_feasible}(x),
\qquad
\Pi_{X_{\mathrm{ad}}}(x),
\qquad
N_{X_{\mathrm{ad}}}(x),
\qquad
\text{multiplier or complementarity relation}.
$$

For a box constraint,

```math
U_{\mathrm{ad}}
=\left\{u\in U:u_{a}\leq u\leq u_{b}\ \text{a.e.}\right\},
```

the constraint owns the projection or normal-cone operation.  It MUST NOT
modify a residual term.

### 3.13 Requirement policy

A requirement policy records a condition that cannot be inferred safely.  It
MUST state:

| Field | Requirement |
|---|---|
| Subject | component and operation requiring the condition |
| Kind | trace, product, point evaluation, fractional norm, nullspace, endpoint, or discrete-only meaning |
| Status | provided, user-assumed, or selected discrete realization |
| Policy | concrete gauge, lifting, trace realization, sensor rule, or norm realization |
| Scope | continuous semantics, discrete compilation, or both |

For transposition, fractional metrics, boundary $H^{1}$ metrics, and partial
Dirichlet interfaces, a selected discrete realization MUST also be carried by
the typed requirement-policy payload. Implementations MAY expose a readable
policy description, but validation and lowering MUST use the structured
identifiers and enum selections.

The validator MAY reject missing policies.  It MUST NOT turn a provided or
user-assumed policy into a claim of mathematical proof.

## 4. Boundary and temporal protocols

### 4.1 Boundary protocol

There is intentionally no universal boundary-condition residual class.

| Boundary statement | Required composition |
|---|---|
| Fixed Dirichlet data | state-space restriction plus fixed-data lifting transformation |
| Dirichlet control | control space, lifting transformation, physical-state reconstruction, residual and observation composition |
| Neumann data or control | boundary residual term using a declared trace pairing |
| Robin data | boundary bilinear residual term plus boundary functional term |
| Periodicity or hanging relation | semantic identification metadata plus discrete constraint realization |
| Pure Neumann problem | residual terms plus nullspace, compatibility, and gauge policies |
| Inflow/outflow transport | selected transport residual form plus oriented boundary-region and trace policy |

The decisive distinction is

$$
\text{Dirichlet control: } y_{\mathrm{phys}}=\widehat y+L_{D}u,
$$

versus

$$
\text{Neumann control: }
\langle E_{u}(u),v\rangle=-\langle u,\gamma v\rangle.
$$

The first is a transformation and chain-rule problem.  The second is a
residual-term problem.  They MUST NOT share a generic boundary-load interface.

### 4.2 Temporal protocol

Time dependence is represented by spaces, regions, data, residual terms, and
losses.  It MUST NOT require a separate parabolic or hyperbolic problem base
class.

For a parabolic state,

$$
Y=L^{2}(0,T;V)\cap H^{1}(0,T;V^{\ast}),
\qquad
Z=L^{2}(0,T;V),
$$

the time term is a residual term

$$
e_{\mathrm{time}}(y;v)
=\int_{0}^{T}\langle\dot y,v\rangle_{V^{\ast},V}\mathrm{d}t.
$$

Initial data MUST be represented by an affine state-space transformation or a
separate equation/constraint block.  Terminal losses MUST be ordinary loss
components.  The selected integration-by-parts convention, initial policy,
and terminal policy MUST be recorded because they determine the adjoint
endpoint condition.

## 5. Problem assembly and validation protocol

### 5.1 Problem specification

The composition root, called here a problem specification, MUST contain only:

1. regions and their relations;
2. spaces, pairings, variables, and data;
3. transformations;
4. equation blocks and residual terms;
5. observations, losses, and objectives;
6. metrics and constraints;
7. requirement policies; and
8. references to a discretisation policy.

It MUST NOT contain solver-specific PDE branches.

### 5.2 Semantic resolution

Resolution proceeds in this order:

1. resolve region identities and partitions;
2. resolve spaces and pairings;
3. resolve variable and data ports;
4. resolve transformations and their physical-field outputs;
5. resolve residual-term inputs and equation test spaces;
6. resolve observations and losses;
7. resolve metrics and constraints;
8. collect unsatisfied requirements.

A successful resolution produces a typed graph.  It does not produce a
theorem about that graph.

### 5.3 Structural validation

The validator MUST reject:

- a residual term without a target equation test space;
- a loss without an observation or variable source of the declared space;
- a transpose action whose declared codomain is not the input dual;
- a natural boundary term without its boundary region and trace pairing;
- a Dirichlet control attached as a load without a reconstruction;
- a metric used to identify a derivative without a declared primal-dual map;
- a pure-Neumann residual without a nullspace policy; and
- a point or fractional operation without a continuous or discrete-only
  policy.

The validator SHOULD report, rather than prove, stated regularity, product,
and compatibility assumptions.

## 6. Discretisation and compilation protocol

### 6.1 Discretisation policy

A discretisation policy selects, for every referenced semantic space,

$$
X_{i}\rightsquigarrow X_{i,h},
\qquad
Z_{a}\rightsquigarrow Z_{a,h},
\qquad
Q_{k}\rightsquigarrow Q_{k,h}.
$$

It MUST declare:

| Selection | Required content |
|---|---|
| Discrete spaces | finite-element family, field components, mesh relation, and trial/test distinction |
| Pairings | quadrature and discrete dual representation |
| Transformations | discrete lifting, trace, restriction, transfer, and parameterisation realization |
| Constraints | affine constraints, periodicity, hanging relations, and gauge realization |
| Execution | assembled or matrix-free operator realization |
| Exceptional policies | point sensors, fractional norms, stabilization, and discrete-only objects |

Galerkin is one policy with $`Y_{h}=Z_{h}`$.  Petrov–Galerkin is a policy with
separately selected $`Y_{h}`$ and $`Z_{h}`$.  The latter MUST preserve the
declared residual target $`Z_{h}^{\ast}`$.

### 6.2 Compiler obligations

The compiler MUST lower each semantic component independently and then
compose the resulting actions.  It MUST provide:

| Semantic component | Required executable output |
|---|---|
| Space and pairing | vector layout, dual layout, pairing application |
| Transformation | value, JVP, transpose-JVP |
| Residual term | residual contribution, JVP, transpose-JVP |
| Equation block | accumulated residual, JVP, transpose-JVP |
| Observation | value, JVP, transpose-JVP |
| Loss | value and dual derivative |
| Metric | apply and required inverse-apply |
| Constraint | selected feasibility/projection/normal-cone/multiplier operation |

The compiler MUST preserve the declared sign convention and pairings.  It
MUST NOT use a raw storage transpose when that differs from the required
pairing-aware transpose.

### 6.3 Executable model

Compilation produces an executable model with:

$$
E_{h}(x_{h}),\qquad
E_{h}'(x_{h})\delta x_{h},\qquad
E_{h}'(x_{h})^{\ast}p_{h},
$$

$$
J_{h}(x_{h}),\qquad
J_{h}'(x_{h}),
$$

and compiled metric and constraint actions.

For every component that declares a transpose action, the implementation MUST
support the adjoint-consistency identity

```math
\langle E_{h}'(x_{h})\delta x_{h},p_{h}\rangle_{Z_{h}^{\ast},Z_{h}}
=
\langle E_{h}'(x_{h})^{\ast}p_{h},\delta x_{h}\rangle_{X_{h}^{\ast},X_{h}}.
```

The test suite MUST exercise this identity for each term and composed block.

## 7. Formulation and solver protocol

### 7.1 Discretize then optimize

The default executable route is

$$
(E,J)
\longrightarrow
(E_{h},J_{h})
\longrightarrow
\text{discrete derivatives}.
$$

For a state-control split $x=(y,u)$, the DTO formulation builder uses

$$
E_{h}(y_{h},u_{h})=0,
$$

$$
D_{y}E_{h}(y_{h},u_{h})^{\ast}p_{h}=D_{y}J_{h}(y_{h},u_{h}),
$$

$$
j_{h}'(u_{h})
=D_{u}J_{h}(y_{h},u_{h})-D_{u}E_{h}(y_{h},u_{h})^{\ast}p_{h}.
$$

A reduced-space solver requires state solution, adjoint solution, reduced
covector, metric, and constraint operations.  An all-at-once solver requires
primal, adjoint, and stationarity residuals with block JVPs and transpose
actions.

A solver MUST NOT request the PDE family or boundary-condition type.

### 7.2 Optimize then discretize

The OTD route is

$$
(E,J)
\longrightarrow
\text{continuous first-order system}
\longrightarrow
\text{compiled optimality system}.
$$

It requires a formulation builder that creates semantic equation blocks for

$$
D_{y}E(y,u)^{\ast}p-D_{y}J(y,u)=0,
$$

$$
D_{u}J(y,u)-D_{u}E(y,u)^{\ast}p=0,
$$

plus constraints as appropriate.

The formulation provenance and executable block shape are separate protocol
choices.  `DTO` identifies a discrete residual/objective pair whose
derivatives are taken after compilation.  `supplied OTD` identifies an
application-provided discrete optimality system whose state, adjoint, and
control-stationarity blocks are supplied directly.  `reduced` and
`all-at-once` describe how those blocks are executed; an all-at-once shape
does not by itself establish OTD provenance.

The selected supplied-OTD interface is deliberately an execution boundary,
not an automatic continuous-adjoint derivation facility.  Its product MUST
declare:

- the primal, adjoint, and control block layouts and their trial/test
  pairings;
- the quadrature and discretisation provenance for each supplied weak block;
- value and block-linearisation actions, including every requested transpose
  action;
- the sign convention and any conversion between the supplied multiplier and
  the framework adjoint; and
- the comparison status relative to a DTO product.

The executor MAY assemble, linearize, solve, and report the supplied blocks,
but MUST NOT reconstruct them from a strong residual or present them as an
exact discrete DTO derivative without an explicit equivalence result.  A
supplied-OTD product is therefore a distinct formulation product and MUST NOT
be represented as a `ReducedDTO` instance.  A request whose supplied adjoint
space, pairing, or sign conversion does not match the declared product MUST
return a formulation diagnostic.

The executable result MUST be marked as OTD and record its discretisation
policy.  It MUST NOT be presented as the exact discrete adjoint unless the
author has established the equivalence.

If DTO produces

$$
A_{h}^{\mathsf T}p_{h}=b_{h}
$$

and OTD produces

$$
\widetilde A_{h}p_{h}=\widetilde b_{h},
$$

the framework MUST allow

$$
\widetilde A_{h}\ne A_{h}^{\mathsf T}.
$$

This difference is expected for some Petrov–Galerkin, stabilized,
time-stepping, quadrature, observation, and lifting choices.

For the selected canonical supplied-OTD shape, a reusable adapter MAY lower
the supplied blocks to the common quadratic KKT boundary. It takes the
negative supplied adjoint-equation JVP as the primal quadratic action, the
state-equation JVP as the equality action, and the supplied residual VJP as
the equality transpose. The adapter MUST retain the supplied block selection
and the declared conversion $`\lambda_{h}=-p_{h}`$; it MUST NOT relabel the
result as a DTO derivative. Non-canonical supplied-OTD block shapes require
their own declared adapter or formulation diagnostic.

### 7.3 Equality-constrained quadratic KKT products

An all-at-once equality-constrained quadratic product is a formulation-level
composition of a primal quadratic operator and a linear equality operator. Let
the primal space be $`X_{h}`$, the residual-multiplier space be
$`\Lambda_{h}`$, and the equality residual space be $`Z_{h}^{\ast}`$. The
product declares

$$
Q_{h}:X_{h}\to X_{h}^{\ast},
\qquad
D_{h}:X_{h}\to Z_{h}^{\ast},
$$

with the quadratic and equality data represented by $`h_{h}\in X_{h}^{\ast}`$
and $`f_{h}\in Z_{h}^{\ast}`$. Its KKT residual is

```math
\begin{bmatrix}
Q_{h}x_{h}-h_{h}+D_{h}^{\ast}\lambda_{h}\\
D_{h}x_{h}-f_{h}
\end{bmatrix}=0,
```

and its block action on an increment is

```math
\begin{bmatrix}
Q_{h} & D_{h}^{\ast}\\
D_{h} & 0
\end{bmatrix}
\begin{bmatrix}
\delta x_{h}\\
\delta\lambda_{h}
\end{bmatrix}.
```

The product MUST declare typed layouts for the primal, multiplier, equality
residual, and KKT residual blocks. It MUST expose the actions of $`Q_{h}`$,
$`D_{h}`$, and $`D_{h}^{\ast}`$, together with the KKT block action and the
corresponding transpose action. An assembled matrix MAY be a lowerer
optimisation; it is not the product contract.

For the selected serial scalar lowerer, an assembled realization MAY expose
the KKT operator as a named three-by-three block matrix with ordering
`[state, control, multiplier]` and row roles
`[state_stationarity, control_stationarity, equality]`. The realization MUST
keep the source and sign of each nonzero block explicit: the state and
control objective Hessian blocks come from the objective derivative, while
the multiplier columns and equality row come from the residual JVP/VJP
pair. Zero objective cross-blocks and the zero multiplier block remain named
blocks. This is lowerer provenance and verification evidence; it does not
add deal.II or assembled-storage requirements to the generic product.

The product MUST record the multiplier convention and any conversion between
the KKT multiplier and the framework adjoint. For the selected scalar target,
the symmetric book multiplier is $`\lambda_{h}=-p_{h}`$; this is a declared
conversion, not a solver-wide sign rule. DTO and supplied-OTD formulation
builders MAY produce the same KKT product shape, but their provenance and
construction evidence remain distinct in the manifest.

The formulation MUST declare the assumptions needed by the selected solve
policy, including the rank or compatibility condition for $D_{h}$ and the
positive-definiteness of $Q_{h}$ on $\ker(D_{h})$. If those assumptions are not
declared or the requested product cannot provide the required pairing, the
compiler MUST return a formulation diagnostic rather than select a
Schur-complement or Krylov policy by implication.

KKT feasibility, stationarity, multiplier conversion, and linear-solve
termination are separate reported quantities. The product MUST NOT identify
the state or control with a particular PDE mass matrix, and a solver MUST NOT
request a PDE family or control placement. Box complementarity, active-set
selection, and multiplier classification belong to the later PDAS product;
they are not part of this quadratic KKT contract.

For a KKT action that is symmetric in its declared pairing, MINRES MAY be
selected only with a compatible symmetric positive-definite preconditioner.
GMRES MUST be used for a declared nonsymmetric action and MAY also be used
for a symmetric action. Variable-work
preconditioners and flexible GMRES belong to the separate preconditioning
extension and MUST not be inferred from a generic KKT product name.

A KKT solver policy MUST name MINRES or GMRES, declare iteration and
tolerance limits, and be checked against the product before execution.
MINRES MUST be rejected for a product that does not declare symmetric-
indefinite compatibility; GMRES remains available for either declared
symmetry. A solve report MUST retain the linear-solver termination separately
from the KKT stationarity and equality/feasibility residual norms, and the
multiplier-to-adjoint conversion remains owned by the product.

### 7.4 Typed box complementarity and primal-dual active sets

A box-complementarity product is a separate formulation service from the
projection operation in Section 3.12. It acts on a declared primal control
space $`U_{h}`$ and its compatible dual $`U_{h}^{\ast}`$, and MUST record the
lower and upper bound representations, the multiplier representation, and the
operations used to compare primal and dual quantities. A box constraint MUST
NOT be treated as a multiplier-bearing complementarity constraint merely
because it already provides coefficientwise projection.

For the selected first realization, $`U_{h}`$ is the cellwise-discontinuous
piecewise-constant volume-control space on the state mesh. The bounds are
coefficientwise values in that same layout, and the declared pairing is the
cellwise $`L^{2}`$ pairing with a positive diagonal Riesz map

$$
R_{u,h}:U_{h}\to U_{h}^{\ast}.
$$

The box multiplier $`\mu_{h}`$ is a dual covector. Active-set classification
MUST first apply the declared dual-to-primal representative
$`R_{u,h}^{-1}\mu_{h}`$; it MUST NOT compare raw dual coefficients with primal
control coefficients. The classification parameter $`c>0`$ is interpreted in
that selected representative. For the framework Lagrangian sign, with
control stationarity $`j_{h}'(u_{h})+\mu_{h}=0`$, define

```math
\begin{aligned}
\widehat\mu_{h} &= R_{u,h}^{-1}\mu_{h},\\
\mathcal A^{+} &= \{i:\widehat\mu_{h,i}+c(u_{h,i}-b_{h,i})>0\},\\
\mathcal A^{-} &= \{i:\widehat\mu_{h,i}+c(u_{h,i}-a_{h,i})<0\},\\
\mathcal I &= \{i:\widehat\mu_{h,i}+c(u_{h,i}-b_{h,i})\leq 0 \leq\widehat\mu_{h,i}+c(u_{h,i}-a_{h,i})\}.
\end{aligned}
```

The three sets MUST be disjoint and exhaustive. A typed selection service
MUST expose the selected representation and provide restriction and
prolongation actions for the free primal variables and the active control
coordinates. It MUST reject incompatible layouts, missing bound data,
nonpositive classification parameters, and any representation that has no
declared primal/dual conversion. Continuous finite-element controls,
facewise controls, quadrature-point inequalities, and state-derived
constraints require separate policies and are not covered by this
realization.

A PDAS service MUST consume an equality-constrained quadratic KKT product
through the contract in Section 7.3. For a selected active set it fixes the
active control coordinates at their declared lower or upper values, composes
the corresponding restricted quadratic and equality actions, and solves the
resulting equality-constrained subproblem through the P6.3 KKT service. The
full primal point is reconstructed before reclassification. The original KKT
multiplier-to-framework-adjoint conversion remains unchanged; the box
multiplier is recovered from the control stationarity covector on active
coordinates and is zero on inactive coordinates.

Every iteration MUST report primal feasibility, dual feasibility,
complementarity, control stationarity, equality feasibility, active-set
changes, and the underlying KKT linear-solve result separately. PDAS MAY stop
only when the active sets are unchanged and the declared full-KKT residual
tolerances are satisfied. Stable active sets alone MUST NOT be reported as
convergence. Failure to construct or solve a restricted KKT product MUST be a
formulation or solve diagnostic, not a silent fallback to projected reduced
optimization.

## 8. Cases that remain cross-cutting

The following cannot be made local to one term.  The stated protocol is
required.

| Case | Why it crosses interfaces | Required procedure |
|---|---|---|
| Spaces and pairings | every derivative, trace, observation, metric, and transpose depends on them | make them first-class typed ports |
| Dirichlet control | control changes the physical state before residual and observation evaluation | use a transformation and compile the full chain rule |
| Very-weak formulation | moving derivatives changes state space, test space, data pairing, and boundary terms | declare a new residual formulation; do not infer it from a strong operator |
| Pure Neumann problem | kernel affects state solve, adjoint, metric, and preconditioner | require compatibility and gauge/nullspace policy |
| Time dependence | temporal spaces, endpoint data, adjoint endpoint, and compiled transpose are coupled | use temporal terms and explicit endpoint policy |
| Point, flux, and fractional objects | continuous map may be unavailable or noncanonical | require a continuous assumption or discrete-only realization |
| Stabilization and nonconformity | compiled residual determines the discrete adjoint | define and test the actual discrete residual and transpose |

These are not exceptions that justify inheritance.  They are explicit
multi-component protocols.

## 9. Laplace binding example

### 9.1 Baseline graph

For

$$
-\Delta y=f+u\quad\text{in }\Omega,
\qquad
y=0\quad\text{on }\Gamma,
$$

use

$$
V=H_{0}^{1}(\Omega),
\qquad
U=L^{2}(\Omega),
\qquad
Q=L^{2}(\omega_{o}).
$$

The registered components are:

| Component | Instance |
|---|---|
| Region | $\Omega$, $\Gamma$, $`\omega_{o}`$ |
| Variables | $y\in V$, $u\in U$ |
| Data | $f\in V^{\ast}$, $`y_{d}\in Q`$, $\alpha>0$ |
| Equation block | test space $V$ |
| Residual terms | $`(\nabla y,\nabla v)_{\Omega}`$, $-\langle f,v\rangle$, $`-(u,v)_{\Omega}`$ |
| Observation | $`C:y\mapsto y\vert_{\omega_{o}}`$ |
| Losses | $`\tfrac{1}{2}\lVert Cy-y_{d}\rVert_{Q}^{2}`$, $`\tfrac{\alpha}{2}\lVert u\rVert_{U}^{2}`$ |
| Metric | selected $`G_{U}`$ |
| Constraint | absent, or a box constraint on $u$ |

The compiled operators are

$$
r_{h}(y_{h},u_{h})=A_{h}y_{h}-f_{h}-B_{h}u_{h},
$$

$$
A_{h}^{\mathsf T}p_{h}
=C_{h}^{\mathsf T}W_{h}(C_{h}y_{h}-d_{h}),
$$

$$
j_{h}'(u_{h})
=\alpha R_{\mathrm{reg},h}u_{h}+B_{h}^{\mathsf T}p_{h}.
$$

### 9.2 Neumann-control delta

Replace the source-control term by

$$
-\langle u,\gamma_{\Gamma_{c}}v\rangle,
\qquad
u\in L^{2}(\Gamma_{c}),
\qquad
\Gamma_{c}\subseteq\Gamma_{N}.
$$

Add a boundary region, boundary-control space, trace map, and boundary
residual term.  The only compiled coupling change is

$$
B_{h}\longrightarrow B_{\Gamma,h},
\qquad
(B_{\Gamma,h})_{j\ell}
=\langle\xi_{\ell},\gamma_{\Gamma_{c}}\psi_{j}\rangle.
$$

The solver protocol is unchanged.

### 9.3 Dirichlet-control delta

Replace the physical state by

$$
y_{\mathrm{phys}}=\widehat y+\ell_{0}+L_{D}u.
$$

Replace the residual by

$$
\langle E_{D}(\widehat y,u),v\rangle
=(\nabla y_{\mathrm{phys}},\nabla v)_{\Omega}-\langle f,v\rangle.
$$

Add a lifting transformation and make both residual and observation consume
$`y_{\mathrm{phys}}`$.  Compilation produces

$$
r_{D,h}(\widehat y_{h},u_{h})
=A_{h}\widehat y_{h}+A_{\mathrm{ext},h}L_{D,h}u_{h}+b_{\ell_{0},h}-f_{h}.
$$

This is not a replacement of $`B_{h}`$ by a boundary load matrix.

The first v1 realization selects the explicit discrete policy
$`\ell_{0,h}=0`$ with one shared nodal trace coefficient for every state DoF on
the complete exterior controlled boundary. It evaluates residual and tracking
on $`y_{\mathrm{phys}}`$ and uses $`P_{h}^{\ast}`$ and $`L_{D,h}^{\ast}`$ for
the two pullbacks. Partial boundary controls, mixed controlled/fixed corners,
interfaces, hanging-node trace relations, and trace-box projections require a
separately declared lifting policy; they are not inferred by this realization.

### 9.4 Coefficient-identification delta

Promote $m$ to a variable and use

$$
\langle E_{m}(y,m),v\rangle
=(m\nabla y,\nabla v)_{\Omega}-\langle f,v\rangle.
$$

Add the parameter derivative port

$$
\langle D_{m}E_{m}(y,m)\delta m,v\rangle
=(\delta m\nabla y,\nabla v)_{\Omega}.
$$

No inverse-problem solver type is required.  The existing equation block,
variable, loss, metric, constraint, and formulation protocols suffice.
The first binary reduced formulation may name its decision port with a
parameter variable rather than a control variable; this changes the semantic
role and residual derivatives, not the generic two-block DTO protocol.

## 10. Required implementation sequence

The first implementation MUST support the following interfaces even if some
advanced instances return an unsupported-capability diagnostic:

1. regions, scalar spaces, pairings, variables, data, and requirement policy;
2. transformations, residual terms, equation blocks, observations, losses,
   metrics, and box constraints;
3. a conforming Galerkin compiler with pairings, affine constraints, assembled
   or matrix-free actions, and exact transpose actions;
4. DTO reduced state-adjoint-gradient formulation;
5. term-level adjoint-consistency and Taylor-remainder tests.

The first concrete instances SHOULD be scalar stationary diffusion-reaction,
fixed Dirichlet data, volume and Neumann controls, distributed and boundary
tracking, $L^{2}$ and $H^{1}$ metrics, and box constraints.

Dirichlet control, coefficient identification, time terms, Petrov–Galerkin,
and OTD SHOULD be added only as new instances or formulation builders.  If an
extension requires a PDE-specific solver branch or a new problem subclass, the
port specification has been violated.
