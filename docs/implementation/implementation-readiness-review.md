# Implementation-readiness review and default policies

## Status and purpose

This document records the outcome of the architecture review. It confirms that
the residual-based formalism is suitable for the project and closes the
implementation ambiguities that would otherwise lead to incompatible
interpretations in a deal.II implementation.

It supplements the architecture record and interface specification. When those
documents permit several choices, the **first default** in this document is
the required choice for the first executable implementation. A later
implementation may select another listed choice only by declaring it in its
discretisation or formulation policy and by providing the required tests.

The result is deliberately conservative:

- the semantic model remains broad;
- the first compiler accepts a small, explicitly supported subset; and
- unsupported combinations fail with a capability diagnostic rather than
  receiving an implicit or mathematically dubious realization.

## 1. Review conclusion

The central model

$$
  \min_{x\in X_{\mathrm{ad}}} J(x)
  \quad\text{subject to}\quad
  E_{a}(x)=0\in Z_{a}^{\ast}
$$

is correct for the intended family of PDE-constrained problems. In
particular, the following decisions are retained:

- residuals target declared test-space duals;
- derivatives are covectors, not silently identified with functions;
- a metric is an algorithmic primal–dual identification, distinct from an
  objective regularisation term;
- transformations carry liftings, parameterisations, restrictions, and their
  chain rule;
- natural and essential boundary conditions have different representations;
- the default solver route is discretise-then-optimise (DTO); and
- PDE terms must not contain solver-specific branches.

The formalism is therefore feasible. What was missing was not another
mathematical primitive, but an unambiguous contract between semantic objects,
discrete linear algebra, and formulation builders.

## 2. Decision rule

Every feature decision follows this order:

1. State the continuous source and target spaces and pairings.
2. Choose a discrete realization, including its coordinate representation.
3. State every action required by the selected formulation and solver.
4. Reject the feature if no registered lowerer and no selected policy can
   provide those actions.
5. Test the compiled actions, rather than relying on a continuous formula to
   imply their correctness.

Semantic validation records assumptions; it does not prove regularity,
well-posedness, inf–sup stability, or convergence.

## 3. Discrete algebra and pairings

### 3.1 The missing distinction

At the discrete level, a primal FE coefficient vector and a covector happen
to use the same backend storage type in many deal.II configurations. They are
nevertheless different mathematical objects. Treating them as interchangeable
is the main route to missing mass matrices and incorrect adjoints.

For every semantic space $S$, compilation MUST create the following typed
objects:

$$
  \mathrm{Primal}(S_{h}),\qquad \mathrm{Covector}(S_{h}),\qquad
  \mathrm{pair}_{S}:\mathrm{Covector}(S_{h})\times\mathrm{Primal}(S_{h})\to\mathbb R.
$$

The wrappers may contain the same `Vector` implementation, but public operator
interfaces MUST preserve the distinction. Products of spaces use typed block
primal vectors and typed block covectors.

### 3.2 Candidate representations and the first default

| Choice | Description | Decision |
|---|---|---|
| Dual-coefficient representation | A residual vector stores $`r_{j}=\langle E,\psi_{j}\rangle`$, and `pair(r,p)` is the coefficient dot product $r^{\mathsf T}p$. | **First default.** It matches ordinary FE assembly and makes a tested residual an actual covector. |
| Riesz-representative representation | A dual element is stored as a primal vector $\bar r$, with $`\langle E,p\rangle=\bar r^{\mathsf T}M_{S}p`$. | Allowed later, but it MUST expose the mass/pairing map explicitly. |
| Untyped backend vector | The meaning of a vector depends on the caller. | Forbidden. |

Under the default, an assembled residual Jacobian has test rows and trial
columns:

$$
  J_{ji}=\langle E'(x)\varphi_{i},\psi_{j}\rangle.
$$

The JVP returns a covector in $`Z_{h}^{\ast}`$, and the pullback of a test vector
returns a covector in $`X_{h}^{\ast}`$. Its coordinate action is the storage
transpose $J^{\mathsf T}p$ *because* the default dual representation has
been declared. This is not a general permission to use a raw storage
transpose.

Every compiled space also records its ownership/ghost layout, constraints,
and backend. Mixed fields use a block layout; they MUST NOT be flattened
without a reversible block map.

### 3.3 Mandatory transpose test

For every term, map, and composed equation block, test random compatible
directions and seeds:

```math
 \mathrm{pair}_{Z}(E_{h}'(x_{h})\delta x_{h},p_{h})
 =
 \mathrm{pair}_{X}(E_{h}'(x_{h})^{\ast}p_{h},\delta x_{h}).
```

The test is performed after applying the same constraint and distribution
rules used in production. For deterministic assembled double-precision tests,
the target relative discrepancy is near floating-point assembly error; MPI
and matrix-free tests may use a documented looser tolerance.

## 4. Composition graph and derivative ownership

### 4.1 Candidate designs

| Choice | Consequence | Decision |
|---|---|---|
| Each residual term handles all upstream transformations itself | Simple local code, but chain rules are duplicated and a transformation can be accidentally omitted from an observation. | Not the default. Allowed only for an opaque custom compiled term. |
| A central typed directed acyclic graph evaluates values, JVPs, and pullbacks | One composition rule works for residuals, observations, and liftings; shared values can be cached. | **First default.** |
| Symbolic expression system with automatic lowering | Attractive future front end, but it is a separate language and cannot be inferred from a strong PDE label. | Future extension only. |

### 4.2 Required compiled-node protocol

Every compiled map owns direct input ports and one output port. It provides:

~~~text
evaluate(inputs, context) -> output
jvp(inputs, input_tangents, context) -> output_tangent
vjp(inputs, output_seed, input_covector_accumulator, context)
~~~

VJP means vector–Jacobian product / pullback. The name is preferred to
“transpose-JVP” because the seed has a type:

- for a transformation or observation $T:X\to Y$, the seed is in $`Y_{h}^{\ast}`$
  and the result accumulates into $`X_{h}^{\ast}`$;
- for a residual $E:X\to Z^{\ast}$, the seed is an adjoint variable in $`Z_{h}`$
  and the result accumulates into $`X_{h}^{\ast}`$.

Residual-term signs belong in the term's value and local derivative. The
global Lagrangian sign belongs only in a formulation builder.

The graph MUST be acyclic after data and variable inputs are distinguished.
It MUST reject ambiguous writers to a physical-field port. A field needed by
both a residual and an observation is represented by one transformation output
with fan-out, not by duplicated lifting code.

Evaluation context is keyed by the evaluation point, mesh revision, data
revision, quadrature policy, and relevant time-step identifier. Reusing a
cached value after any of those changes is forbidden.

### 4.3 Objective assembly rule

Every loss consumes one map output. A direct loss on a variable is resolved to
an explicit identity map during semantic resolution; it is not a second,
special objective path. A loss involving several variables first receives a
product/concatenation map output. This keeps the rule

$$
  J'(x)=\sum_{k} O_{k}'(x)^{\ast}\Phi_{k}'(O_{k}(x))
$$

literal and removes the ambiguity between “loss source is an observation” and
“loss source is a variable.”

## 5. Lowering and compiler capabilities

A semantic declaration is not itself executable. Each component kind needs a
registered lowerer that can produce the selected discrete actions.

| Choice | Consequence | Decision |
|---|---|---|
| Accept arbitrary semantic terms and attempt generic compilation | Impossible without an expression language or user code. | Forbidden. |
| Small built-in lowerer registry | Predictable compilation and testing; feature growth is explicit. | **First default.** |
| User-provided C++ compiled term | Extensible, but the author owns all requested actions and tests. | Allowed after the built-in path works. |
| Symbolic/code-generated lowerer | Possible future front end. | Deferred. |

The compiler performs four separate checks:

1. **Structural:** ports, region kinds, field shapes, dual status, graph
   acyclicity, and declared transformations agree.
2. **Analytical-policy:** a required trace, product, sensor, nullspace, or
   fractional-norm policy is supplied. This produces metadata, not a proof.
3. **Lowerability:** every node has a lowerer for the selected mesh, FE,
   backend, and execution mode.
4. **Formulation capability:** the selected solver has every requested
   derivative, metric, constraint, and linear-solve operation.

Failure in any category is a diagnostic with the component, missing
capability, and candidate remedies. The compiler may fuse logically
independent cell or face contributions into one loop; independent ownership
does not require inefficient separate traversal or separate matrices.

## 6. Formulation ownership and derivative order

### 6.1 Reduced formulation

A reduced solver cannot infer what should be eliminated merely from a variable
role named state. It MUST receive a reduced-formulation specification:

~~~text
implicit variable blocks       state / flux / auxiliary blocks to eliminate
equation blocks                residual blocks defining those variables
decision variable blocks       controls and/or parameters retained by the optimizer
state solve policy             linear or nonlinear solve capability and tolerances
adjoint solve policy           transpose solve capability and tolerances
nullspace/uniqueness policy    gauge, quotient, or multiplier choice when needed
constraint policy              selected feasible-set operation
~~~

The compiler checks that the selected discrete state Jacobian has compatible
block dimensions and that its requested solve policy is available. It does
not claim to prove invertibility. A state solve may be nonlinear; the
formulation records the solve tolerance because inexact state/adjoint solves
change the reduced derivative.

### 6.2 All-at-once and second order

First-order residual ports are sufficient for state, adjoint, and reduced
gradient methods. They are not sufficient for a nonlinear Newton step on the
KKT equations: differentiating the adjoint equation requires terms such as

$$
  D_{x}\bigl(E_{h}'(x_{h})^{\ast}p_{h}\bigr)[\delta x_{h}]
  \quad\text{and}\quad
  J_{h}''(x_{h})\delta x_{h}.
$$

| Choice | Decision |
|---|---|
| Exact second-order action supplied by terms/maps | Future general Newton/SQP contract. |
| Automatic differentiation of local kernels, with an independently tested pullback | Permitted future implementation strategy. |
| Gauss–Newton or quasi-Newton approximation | Permitted only when the formulation declares the approximation. |
| First-order L-BFGS reduced solver | **First default for nonlinear problems.** |
| All-at-once KKT Newton for linear-quadratic problems | Allowed after the first reduced DTO path. |
| General nonlinear all-at-once Newton/SQP | Unsupported until the second-order contract exists. |

The future second-order protocol should expose a Lagrangian Hessian-vector
action, rather than force storage of a tensor. Nonsmooth losses and
constraints must advertise that this action is unavailable.

### 6.3 DTO and OTD

DTO is the only first-default provenance:

$$
 (E,J)\longrightarrow(E_{h},J_{h})\longrightarrow\text{discrete derivatives}.
$$

Every executable model records provenance = DTO, the exact quadrature and
constraint policy, and its data/lifting realization. The direct v0 baseline
selects DTO only. If OTD is registered, it is a separate formulation builder
whose output records provenance = OTD; it may never be labelled the exact
discrete adjoint without an explicit equivalence test.

## 7. Metrics, regularisation, and constraints

### 7.1 Keep these notions separate

An objective term

$$
  \frac{\alpha}{2}\lVert u\rVert_{H^{1}}^{2}
$$

is a regularisation loss and changes $J'(u)$. A search metric maps the
already-computed covector to a direction and does not change $J'(u)$, the
state equation, or the adjoint equation. No document or API may use
“regularisation metric” for both.

A conforming discrete realization may use the full $H^{1}$ matrix
$M_{u}+K_{u}$ in the objective term and its derivative. Its regularisation
factory must remain separate from either an $L^{2}$ mass-matrix Riesz map or
an $H^{1}$ metric with $G=M_{u}+K_{u}$; the positive mass term makes the
latter coercive, and its inverse changes only the search direction.

### 7.2 Search-metric choices

| Choice | Meaning | Decision |
|---|---|---|
| $L^{2}$ metric | $`G=M_{U}`$ in the selected control realization. | **First default.** |
| $H^{1}$ Sobolev metric | Select a search space $P\subseteq U$, injection $\iota:P\to U$, and coercive Riesz map $G:P\to P^{\ast}$. Direction formation solves $Gg=\iota^{\ast}j'$. | Allowed with an explicit coercivity policy and a separately declared compatible constraint, if any. |
| $H^{-1}$-type metric | Must state the actual Hilbert space and operator. If $`(v,w)_{-1}=(A^{-1}v,w)_{L^{2}}`$, then the metric operator is $A^{-1}$ and its inverse is $A$; this is not the same operation as an $H^{1}$ Sobolev-gradient solve. | No generic default. The named P5.2 realization on independent homogeneous-Dirichlet `FE_Q` coordinates uses $G_h=M_hK_h^{-1}M_h$ and is supported only with its explicit boundary and solve policy. |
| Fractional metric | Requires a named discrete realization and spectral/extension/auxiliary problem policy. | Unsupported initially. |

An $H^{1}$ metric includes a positive zero-order term or an explicit
mean/boundary condition; the seminorm alone is not invertible. A metric
inverse is a solver operation and must report its tolerance and
preconditioner policy.

### 7.3 Box constraints

Pointwise box constraints have no FE-independent projection rule.

| Discrete set realization | Consequence | Decision |
|---|---|---|
| Cellwise constants on the same volume mesh | Coefficient clipping exactly represents cellwise bounds; $L^{2}$ metric projection is local. | **First default for volume controls.** |
| Nodal bounds for continuous Lagrange controls | Nodal clipping need not imply bounds between nodes for all polynomial degrees. | Allowed only with a stated nodal, not $a.e.$, semantics. |
| Quadrature-point inequalities | Faithful sampled constraint but projection is a constrained optimization problem. | Future extension. |
| Exact FE obstacle/set projection | Metric-dependent global problem. | Future extension. |

The initial volume-control space is therefore discontinuous piecewise
constant (`FE_DGQ(0)`) on the state mesh, and lower/upper bound data must be
piecewise constant on that mesh. The initial box solver uses the declared
$L^{2}$-metric projection. It MUST NOT claim that coefficient clipping is an
$H^{1}$-metric projection. A solver requesting a projection supplies the
actual metric service. The box accepts only the positive-diagonal metric
realization coupled by its owning lowerer; a human-readable metric identifier
is provenance, not compatibility proof. Absence of that coupled capability is
a failure, not permission to substitute clipping.

## 8. Boundary conditions, liftings, and nullspaces

### 8.1 Fixed essential data

All compiled state operators use independent unknown coordinates. With
$`P_{h}`$ the homogeneous/hanging/periodic reconstruction, the physical field
is represented as

$$
  y_{\mathrm{phys}}=P_{h}\widehat y_{h}+\ell_{0,h}.
$$

Residuals and observations evaluate the physical field. Pullbacks use the
matching $`P_{h}^{\ast}`$. `AffineConstraints` is an implementation mechanism for
$`P_{h}`$ and must be applied consistently in residual, observation, JVP, and
VJP paths.

The first default supports homogeneous or fixed, time-independent Dirichlet
data only. Data discretisation is declared; no implicit nodal interpolation is
permitted.

### 8.2 Dirichlet control

The required representation is

$$
 y_{\mathrm{phys}}=P_{h}\widehat y_{h}+\ell_{0,h}+L_{D,h}u_{h}.
$$

| Choice | Decision |
|---|---|
| Treat controlled data as a boundary load | Forbidden. |
| Rebuild affine constraints for each control and hide the dependence | Not a public derivative contract; not the default. |
| Expose an explicit lifting/reconstruction $`L_{D,h}`$, including JVP and VJP | Required whenever this capability is registered. |

$`L_{D,h}`$ must state its boundary-control discretisation, interior
extension, corner/interface compatibility, and behavior on fixed Dirichlet
portions. The discrete problem records the lifting choice because different
liftings can produce different discrete intermediate systems. Dirichlet
control is excluded from the direct v0 baseline rather than simulated by a
Neumann-like coupling. The [v1 capability table](v1/semantic-compiler.md#registered-capabilities)
owns the status and bounds of any registered realization.

#### 8.2.1 $L^{2}(\Gamma)$ Dirichlet control by a conforming trace subspace

The selected continuous parent for the first Chapter 5.11.2 specialization is
the unconstrained full-boundary Laplace problem with
$f,z_{d}\in L^{2}(\Omega)$, $\beta>0$, state $y\in L^{2}(\Omega)$, and control
$u\in L^{2}(\Gamma)$. It declares the transposition residual

```math
\langle E_{\mathrm{tr}}(y,u;f),\psi\rangle
=(y,-\Delta\psi)_\Omega-(f,\psi)_\Omega
+(u,\partial_{n}\psi)_\Gamma,
\qquad
\psi\in Y=H^{2}(\Omega)\cap H^{1}_{0}(\Omega).
```

The model author must declare a domain on which
$-\Delta:Y\rightarrow L^{2}(\Omega)$ has the required isomorphism property;
convex and $C^{2}$ domains are the source cases. A claim that the exact
optimal control is in $H^{1/2}(\Gamma)$ and its state is in $H^{1}(\Omega)$
uses the stronger $C^{2}$ hypothesis of Proposition 5.16. A convex polygonal
test case may rely on transposition well-posedness and the regularity of its
selected discrete subspace, but must not claim that stronger continuous
bootstrap without a separate result.

The first discrete policy does not introduce an $H^{2}$-conforming state or
test element. It selects the existing complete-boundary nodal trace space

```math
U_{h}=\mathrm{tr}_{\Gamma}V_{h}
\subset H^{1/2}(\Gamma)\subset L^{2}(\Gamma)
```

and the boundary $L^{2}$ pairing independently for the control loss and
control metric. Both lower to the same boundary mass matrix but remain
separate semantic components. Because each $u_{h}$ has an $H^{1/2}$ trace,
its continuous variational and transposition states coincide; the lifted
$H^{1}$ finite-element solve is the conforming Galerkin approximation of that
state. The compiled product must record the continuous parent space, the
conforming discrete subspace, the equivalence used by the lowerer, and the
domain-regularity declaration. It must reject a discontinuous or facewise
Dirichlet control under this policy.

For the pure Laplace bilinear form $a$, define the discrete outward conormal
covector using the same lifting as the physical reconstruction:

```math
\langle q_{h},v_{h}\rangle
=a(L_{D,h}v_{h},p_{h})
-(y_{h}-z_{d},L_{D,h}v_{h})_{\Omega}.
```

The discrete reduced covector is

```math
j_{h}'(u_{h})v_{h}
=\beta(u_{h},v_{h})_{\Gamma}-\langle q_{h},v_{h}\rangle.
```

For a smooth exact adjoint, $q_{h}$ represents the pullback of
$\partial_{n}p$; it is not defined by pointwise differentiation of an
`FE_Q` field. The formula is independent of the interior extension chosen by
$L_{D,h}$ when the discrete adjoint equation is satisfied. A focused contract
must compare this conormal covector with the existing composed lifting
pullback and verify the stationarity sign. A general transposition lowerer
remains a separate requirement for nonconforming $L^{2}$ boundary controls,
normal-flux observations, and Dirac-source adjoints.

### 8.3 Natural boundary data and controls

For the first boundary-control extension, $`\Gamma_{c}`$ is a marked collection
of faces on the same static triangulation and the control is facewise
constant. The boundary term is assembled with `FEFaceValues`; its coupling and
pullback are tested by the ordinary adjoint identity.

Continuous trace FE controls, separate surface meshes, and nonmatching
control meshes are valid later choices, but each needs an explicit trace or
transfer map. They are not inferred from a boundary ID.

### 8.4 Pure Neumann policy

| Choice | Consequence | Decision |
|---|---|---|
| Pin one state DoF | Easy but mesh- and point-dependent gauge. | Not the default. |
| Work in a quotient/mean-zero subspace | Clean mathematically, but needs a robust discrete realization. | Equivalent permitted realization. |
| Augment with one mean constraint/multiplier | Makes the gauge and compatibility visible in the algebra. | **First default.** |

This policy requires the state and adjoint realizations to use the same
visible gauge, check all discrete compatibility conditions, and record the
gauge and solve policy in provenance. The
[v1 capability table](v1/semantic-compiler.md#registered-capabilities) owns
the exact bounds of any registered pure-Neumann target.

## 9. Mesh, data, geometry, and time

### 9.1 Mesh and data defaults

| Area | First default | Deferred alternatives |
|---|---|---|
| Geometry | Fixed geometry and static triangulation throughout one optimisation solve. | Shape/topology optimisation and moving meshes. |
| Mesh relation | One state mesh; volume controls and observations use it; boundary controls use marked faces of it. | Nonmatching meshes and explicit transfer maps. |
| FE state | Scalar conforming `FE_Q` of selected degree. | Vector, mixed, DG, and Petrov–Galerkin systems. |
| Execution | Assembled residual/Jacobian actions. | Matrix-free actions after assembled equivalence tests. |
| Data rule | Analytic data evaluated at declared quadrature points, or an explicitly named projection/interpolation map. | Implicit sampling/interpolation. |

Every data item has a discrete data rule. Changing data, mesh, quadrature, or
a mesh-dependent projection invalidates compiled caches. Geometry parameters
are not ordinary coefficient parameters: shape derivatives, mesh-motion maps,
and geometric conservation policies are outside the first scope and must be
introduced as a separate transformation family.

### 9.2 Time dependence

Time dependence preserves the semantic residual model, but it needs a
dedicated execution policy.

| Choice | Decision |
|---|---|
| Global space–time residual | Valid future option. |
| Fixed-grid one-step residual with the complete trajectory represented in $`E_{h}`$ | **First temporal default.** Use fixed-step backward Euler. |
| Opaque forward time integrator with an independently coded backward adjoint | Forbidden as a DTO implementation unless it exposes the exact residual transpose. |
| Adaptive time steps, events, remeshing, or nonlinear stopping-dependent step decisions | Unsupported until replay and differentiated-control policies exist. |

The temporal compiler owns trajectory storage/checkpointing and replay. It
must return the transpose of the actual time-discrete residual, including
initial and terminal policies. Spatial terms remain reusable kernels, but time
is not merely another cell integral.

## 10. Solver and backend policies

The executable operator contract remains independent of PDE family, but it
does need numerical capabilities. A solve policy is therefore separate from
residual terms and declares:

~~~text
linear/nonlinear solve kind
operator traits actually relied upon (for example SPD, nonsymmetric)
preconditioner and backend
absolute/relative stopping rules
nullspace treatment
whether an approximate solve is permitted in derivative evaluation
~~~

The first implementation uses assembled operators and a single deal.II
backend suitable for serial unit tests. The public typed-vector contract must
not depend on that backend, so a distributed Trilinos/PETSc backend can be
added later. Matrix-free mode is a separate lowerer capability and must pass
the same value/JVP/VJP tests as the assembled mode before it is usable in an
optimizer.

## 11. Baseline slice and extension policies

The selected direct v0 baseline is intentionally narrower than the complete
semantic language:

~~~text
fixed geometry; one static mesh; scalar stationary diffusion-reaction
conforming FE_Q state; homogeneous or fixed Dirichlet data
FE_DGQ(0) volume control on the same mesh
distributed tracking and quadratic L2 control regularisation
assembled DTO residual/JVP/VJP; L2 search metric
unconstrained reduced gradient, then L2-projected box gradient
~~~

This is a policy baseline, not a mutable release ledger. The exact current v1
graphs, target implementations, bounds, tests, and exclusions live in the
[v1 capability table](v1/semantic-compiler.md#registered-capabilities) and its
[exclusions](v1/semantic-compiler.md#exclusions). Task order and completion
state live only in the [implementation roadmap](../planning/implementation-roadmap.md).

Every extension must declare its own lowering boundary, capability diagnostic,
and tests. Fixed reconstruction precedes controlled-essential lifting;
regularisation remains distinct from the search metric; nullspace policy is
explicit; and coefficient transformations are never inferred from positivity
bounds. No extension may add a PDE-name branch to a solver.

## 12. Required verification and provenance

Before an executable component is accepted, tests cover:

1. residual value against a manufactured or independently assembled case;
2. JVP finite-difference/Taylor remainder for every differentiable term;
3. VJP/adjoint consistency for every term, transformation, observation, and
   composed block;
4. objective directional derivative;
5. reduced-derivative Taylor remainder, including an actual state and adjoint
   solve;
6. feasibility and optimality behavior for every supported discrete
   constraint realization; and
7. constraint/lifting/nullspace behavior after the exact same
   `AffineConstraints` and distribution path used in production.

Each output records a compilation manifest containing semantic identifiers,
FE and mesh selections, quadrature, data rules, pairings/dual
representation, constraint and lifting realization, nullspace policy, solver
policy, provenance (DTO or OTD), and declared assumptions. This manifest is
emitted with solver diagnostics so results can be reproduced and DTO/OTD or
assembled/matrix-free comparisons cannot be mistaken for the same discrete
problem.

## 13. Acceptance criterion for future features

A proposed feature is accepted only if its design answers all of the
following:

1. What typed primal and dual ports does it consume and produce?
2. Which graph node owns its value, JVP, and VJP?
3. Which lowerer realizes it for the chosen FE spaces, mesh relation, and
   execution mode?
4. Does it require a user analysis assumption or a selected discrete policy?
5. Which formulation capabilities and derivative order does it need?
6. How are constraints, liftings, nullspaces, data rules, and cached values
   represented?
7. Which value, adjoint-consistency, and Taylor/KKT tests demonstrate the
   compiled behavior?

If any answer is missing, the correct result is an unsupported-capability
diagnostic, not an implicit fallback.
