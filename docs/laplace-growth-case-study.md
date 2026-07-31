# Laplace source control as a composition growth test

## Purpose

This note tests the proposed architecture against one concrete problem.  We
start with the standard Laplace distributed-control example, state what the
program would actually receive, lower it to the objects required by a
numerical solver, and then add one variation at a time.

The test is successful only if a variation changes its owning components and
their declared dependencies, without introducing classes such as
`NeumannBoundaryControlProblem` or `ParabolicTrackingSolver`.

The pseudo-input below describes semantic objects and connections.  It is not
proposed user syntax.

## 1. Baseline: distributed source control for Laplace's equation

Let $\Omega$ be a bounded domain.  The textbook problem is

$$
  \min_{u\in U_\mathrm{ad}}\quad
  J(y,u)=\frac12\int_{\omega_o}(y-y_d)^2
           +\frac\alpha2\int_\Omega u^2
$$

subject to

$$
  -\Delta y=f+u\quad\text{in }\Omega,
  \qquad y=0\quad\text{on }\Gamma=\partial\Omega. \quad\text{(1)}
$$

Here $`\omega_o`$ is an observation subdomain, possibly all of $\Omega$, and
$`U_\mathrm{ad}`$ is initially all of $L^2(\Omega)$.

### 1.1 What (1) does *not* give to the program

The strong equation does not specify a residual target, test space, treatment
of inhomogeneous data, discrete trial/test spaces, or a gradient geometry.
The program cannot safely infer these.  In particular, the expression
$-\Delta y$ must not be the sole input from which the framework guesses a weak
form and boundary treatment.

### 1.2 The semantic input

For the usual weak formulation, choose

$$
  Y=Z=V:=H^1_0(\Omega),\qquad U=L^2(\Omega),\qquad E:Y\times U\to V^{*}=H^{-1}(\Omega).
$$

The program-facing semantic graph is:

```text
regions
  Omega             : volume
  Gamma             : boundary (Dirichlet)
  omega_o           : volume observation region

spaces and pairings
  state Y           : H^1_0(Omega), scalar field
  state test Z      : H^1_0(Omega)
  control U         : L^2(Omega)
  observation Q     : L^2(omega_o)
  residual pairing  : V* x V
  loss pairings     : L^2(omega_o), L^2(Omega)

data
  f                 : source in a declared data space
  y_d               : target in Q
  alpha > 0         : scalar regularisation weight

variables
  y : state, space Y
  u : control, space U, admissible set U_ad=U

state equation block (test space Z)
  diffusion term     : (y,v) |-> integral_Omega grad(y).grad(v)
  source-data term   : (f,v) |-> - integral_Omega f v
  source-control term: (u,v) |-> - integral_Omega u v

objective
  observation C     : y |-> y restricted to omega_o, into Q
  tracking loss     : q |-> 1/2 ||q-y_d||_Q^2
  control loss      : u |-> alpha/2 ||u||_U^2

algorithmic metric (only if a gradient method requires one)
  G_U : U -> U*     : chosen L^2 Riesz map
```

The resulting mathematical objects are

$$
  \langle E(y,u),v\rangle=(\nabla y,\nabla v)_\Omega
    -(f,v)_\Omega-(u,v)_\Omega, \quad\text{(2)}
$$

and

$$
  J(y,u)=\tfrac12\|Cy-y_d\|_Q^2+
           \tfrac\alpha2\|u\|_U^2. \quad\text{(3)}
$$

The initial architecture therefore needs no `LaplaceProblem` object.  It only
needs one equation block containing three residual terms, one observation,
two losses, and one optional metric.

### 1.3 Downstream effect of every baseline input

| Semantic input | Compiled contribution | Consumed by solver |
|---|---|---|
| $Y$, $Z$, field shape | state/test DoF layouts and trial/test FE spaces | state and adjoint vector spaces |
| $U$ | control DoF layout/FE space | control vector and admissible-set operations |
| diffusion term | state residual and its linearisation $`A_h`$ | state solve and adjoint transpose $`A_h^{*}`$ |
| source data $f$ | load vector/function $`f_h`$ | state residual value only |
| source-control term | coupling $`B_h:U_h \rightarrow Z_h^{*}`$ | state solve and reduced derivative through $`B_h^{*}`$ |
| observation $C$ | restriction/interpolation $`C_h:Y_h \rightarrow Q_h`$ | objective and adjoint right-hand side through $`C_h^{*}`$ |
| target $`y_d`$ | observation data $`d_h`$ | objective and adjoint right-hand side |
| tracking loss | observation-space weight/pairing $`W_h`$ | objective and adjoint right-hand side |
| control loss | covector map $`\alpha R_{\mathrm{reg},h}`$ | reduced derivative |
| algorithmic metric $`G_U`$ | $`G_h:U_h \rightarrow U_h^{*}`$, possibly a linear solve | conversion of a reduced derivative to a search direction |
| $`U_\mathrm{ad}`$ | identity/projection/normal-cone operation | stopping condition and constrained update |

The regularisation map $`R_{\mathrm{reg},h}`$ and the algorithmic metric $`G_h`$ are
intentionally separate.  In this basic example both may be the $L^2$ mass
matrix, but they represent different choices and will diverge in later cases.

## 2. What compilation produces

Choose a conforming Galerkin policy, for example

$$
  Y_h=Z_h\subset H^1_0(\Omega),\qquad U_h\subset L^2(\Omega).
$$

The compiler lowers the graph into an `ExecutableModel`.  It need not expose
matrices: the following are operator contracts, usable either assembled or
matrix-free.

```text
ExecutableModel
  layouts
    state: Y_h, control: U_h, state-test: Z_h, observation: Q_h

  residual
    r_h(y_h,u_h) in Z_h*
    apply_jacobian((dy_h,du_h)) in Z_h*
    apply_transpose(p_h) in Y_h* x U_h*

  objective
    J_h(y_h,u_h)
    derivative_y_h in Y_h*, derivative_u_h in U_h*

  transformations and observations
    evaluate / JVP / transpose-JVP for each compiled map

  optimisation geometry
    G_h and (where supported) G_h^{-1}
    U_ad,h projection, normal-cone, or multiplier operations

  metadata/policies
    pairings, quadrature, constraints/liftings, nullspace policy,
    and whether an object is continuous or explicitly discrete-only
```

For this linear quadratic example, an assembled realization may be written
in coordinates as

$$
  r_h(y,u)=A_hy-f_h-B_hu, \quad\text{(4)}
$$

$$
  J_h(y,u)=\tfrac12(C_hy-d_h)^{T} W_h(C_hy-d_h)
            +\tfrac\alpha2u^{T}R_{\mathrm{reg},h}u. \quad\text{(5)}
$$

$`A_h^{T}`$, $`B_h^{T}`$, and $`C_h^{T}`$ in formulas below mean the transpose that
represents the declared discrete pairings.  An implementation must not assume
that a raw storage transpose is correct if its vector representations use
non-Euclidean pairings.

The **final output passed to a solver** is not the original strong PDE.  It is
one of the following generic formulations built from `ExecutableModel`.

```text
Reduced formulation
  solve_state(u_h) -> y_h satisfying r_h(y_h,u_h)=0
  solve_adjoint(y_h,u_h) -> p_h
  reduced_derivative(u_h) in U_h*
  metric/constraint operations

All-at-once KKT formulation
  unknown block (y_h, u_h, p_h [, multipliers])
  primal residual, adjoint residual, stationarity residual
  block Jacobian/JVP/transpose actions and preconditioner hooks
```

Neither formulation needs to know that the residual term was “Laplace.”

## 3. Baseline optimality systems and the two workflows

With the global convention

$$
  \mathcal L_h(y,u,p)=J_h(y,u)-p^{T}r_h(y,u),
$$

the discrete equations for (4)–(5) are

$$
  A_hy=f_h+B_hu, \quad\text{(6a)}
$$

$$
  A_h^{T}p=C_h^{T}W_h(C_hy-d_h), \quad\text{(6b)}
$$

$$
  j_h'(u)=\alpha R_{\mathrm{reg},h}u+B_h^{T}p\in U_h^{*}. \quad\text{(6c)}
$$

An unconstrained gradient method uses

$$
  \nabla_{G_h}j_h=G_h^{-1}j_h'(u). \quad\text{(7)}
$$

For box constraints, (6c) becomes a variational inequality or a projected
gradient/active-set condition; the state and adjoint equations do not change.

### 3.1 Discretize then optimize (DTO)

DTO is the direct route from the compiled model:

```text
semantic graph (E, J)
       -> compile chosen residual and objective (E_h, J_h)
       -> differentiate E_h and J_h
       -> solve the discrete reduced or KKT system
```

For the baseline this is exactly (6).  DTO requires the compiler to provide:

1. the discrete residual value and its discrete Jacobian action;
2. the exact transpose of that discrete Jacobian under the selected discrete
   pairings;
3. the discrete objective and its discrete derivative; and
4. discrete metric and constraint operations.

This is the safest default for an implementation library because every solver
acts on the problem it actually evaluates.  It supports assembly, matrix-free
execution, Petrov–Galerkin test spaces, stabilization, time stepping, and
discrete-only observations without assuming that any continuous and discrete
operations commute.

### 3.2 Optimize then discretize (OTD)

OTD first differentiates the continuous selected formulation:

```text
semantic graph (E, J)
       -> continuous state, adjoint, and stationarity relations
       -> choose and compile a discretisation of those relations
       -> solve the resulting discrete optimality system
```

For the baseline weak form,

$$
  (\nabla y,\nabla v)=(f+u,v),
$$

$$
  (\nabla w,\nabla p)=(Cy-y_d,Cw)_Q,
$$

$$
  j'(u)\delta u=\alpha(u,\delta u)_U+(p,\delta u)_\Omega. \quad\text{(8)}
$$

The formal adjoint in (8) comes from the declared residual and its continuous
transpose port; no solver needs a hand-coded rule for Laplace's equation.
The OTD compiler then selects discrete spaces and forms for all equations in
(8).

OTD and DTO may coincide for this simple conforming Galerkin case when the
forms, quadrature, constraints, and pairings are made compatible.  They must
not be assumed to coincide in general.  Differences can arise from
stabilisation, nonlinear linearisation choices, mass lumping, inexact
quadrature, nonconforming/Petrov–Galerkin test spaces, time stepping,
discrete-only observations, or boundary/lifting treatments.

### 3.3 Architectural consequence

The semantic core must support both provenance paths:

| Path | Object differentiated | Solver's adjoint action |
|---|---|---|
| DTO | compiled $`E_h`$, $`J_h`$ | exact transpose of $`E_h'(x_h)`$ |
| OTD | semantic $E$, $J$, then compiled optimality relations | discretisation of the declared continuous adjoint relation |

The first executable slice should implement DTO and test it rigorously.  OTD
can be represented as a separate **formulation builder** over the same
semantic residual/objective ports.  It must record its provenance and never
masquerade as the exact discrete adjoint unless an equivalence policy is
declared and verified.

This is a small addition to the architecture, not a second PDE hierarchy:
`ReducedFormulation` and `KKTFormulation` consume generic ports; an
`OTDFormulationBuilder` consumes their continuous counterparts.

## 4. Progressive extensions

Each row begins from the preceding baseline graph.  “Output change” says what
must be newly available to a solver after compilation.

### 4.1 Change the fixed source or coefficient data

| Add/change | Semantic input | Downstream effect | Output change |
|---|---|---|---|
| New $f$ | replace `Data(f)` | residual value/load changes | new $`f_h`$; no new operator port |
| Reaction coefficient $c$ | data plus reaction residual term $(cy,v)$ | state Jacobian/adjoint includes reaction | $`A_h`$ action changes |
| Transport $b$ | data plus selected transport residual term | residual and transpose reflect the selected weak form | $`A_h`$, $`A_h^{*}`$; no solver branch |

This establishes the simplest composition rule: PDE physics is a residual
term.  It affects state and adjoint operator actions, but never changes the
objective, metric, or optimisation algorithm interface.

### 4.2 Replace distributed control with Neumann control

Let $`\Gamma_c`$ be a named subset of the Neumann boundary.  Replace

$$
  U=L^2(\Omega),\qquad -(u,v)_\Omega
$$

by

$$
  U_\Gamma=L^2(\Gamma_c),\qquad
  -\langle u,\mathrm{tr}_{\Gamma_c}v\rangle. \quad\text{(9)}
$$

| Component | New input | Downstream effect |
|---|---|---|
| Region | $`\Gamma_c`$ boundary region | boundary quadrature/DoF access |
| Control variable | $`u \in U_\Gamma`$ | boundary control layout and bounds |
| Residual coupling | Neumann control term (9) | $`B_{\Gamma,h}:U_{\Gamma,h} \rightarrow Z_h^{*}`$ |
| Requirement | declared trace $`Z \rightarrow L^2(\Gamma_c)`$ or another selected pairing | validator requires an explicit trace realization |
| Control loss/metric | boundary-space versions | $`R_{\mathrm{reg},h}`$, $`G_h`$ on boundary control DoFs |

The compiled adjoint derivative becomes

$$
  j_h'(u)=\alpha R_{\mathrm{reg},h}u+B_{\Gamma,h}^{T}p.
$$

The state solve, generic adjoint workflow, and outer solver do not change;
only the control-to-residual coupling and its transpose do.

### 4.3 Change where and how the state is observed

Replace the volume restriction $C$ by one of:

| Desired feature | Semantic component replaced/added | Required policy | Compiled solver effect |
|---|---|---|---|
| Boundary tracking | trace observation $`C_\Gamma:Y \rightarrow Q_\Gamma`$ | trace and boundary metric | replace $`C_h^{T} W_h(...)`$ in adjoint RHS |
| Flux tracking | normal-flux observation | flux regularity and trace realization | same observation-transpose port, different implementation |
| Point sensors | $`C_{pt}:Y \rightarrow R^m`$ | regularity or explicitly discrete-only status | finite sensor residual and $`C_{pt,h}^{T}`$ adjoint forcing |
| Energy tracking | change loss to a declared bilinear form | positivity/semidefiniteness policy | different objective derivative/adjoint RHS |

No state residual term changes.  The only solver-visible change is the
objective derivative with respect to the state, supplied through the
observation transpose.

### 4.4 Change the regularisation or the gradient geometry

There are two independent modifications.

| Modification | Input owner | Baseline replacement | Solver effect |
|---|---|---|---|
| $H^1$ control regularisation | control `Space` (or loss domain) plus `Loss` | $`\frac{\alpha}{2}\lVert u\rVert_{L^2}^2`$ becomes $`\frac{\alpha}{2}\lVert u\rVert_{H^1}^2`$ | replace $`\alpha R_{\mathrm{reg},h} u`$ in (6c); compile an $H^1$-conforming control realization |
| $H^1$ or $H^{-1}$ gradient | algorithmic `Metric` | $`G_h=M_h`$ becomes chosen Riesz/auxiliary-solve map | replace only $`G_h^{-1}`$ in (7) |

For example, an $H^1$ regularisation produces an objective derivative
$`\alpha(M_h+K_h)u`$ under its selected boundary convention.  It normally
requires an $H^1$ control space (or, more generally, a declared loss domain
inside the ambient control space) and changes the optimal control problem.
Choosing an $H^1$ search metric for the same $L^2$-regularised objective leaves
(6a)–(6c) unchanged and only preconditions the search direction.  This
separation is essential for correct reduced-space solvers.

### 4.5 Add box constraints

For bounds $`u_a \leq u \leq u_b`$, add one constraint component:

```text
Constraint(Box)
  variable       : u
  lower/upper    : data or parameter fields in U
  operations     : feasible-set test, projection and/or normal-cone relation
```

There is no change to $E$, $`A_h`$, $`B_h`$, $`C_h`$, or the adjoint equation.  The
compiled output adds a discrete admissible-set operation.  A projected
gradient solver needs $`\mathrm{projection}_{U_{\mathrm{ad},h}}`$ compatible with its chosen metric;
a primal-dual active-set/KKT solver instead requests complementarity or
multiplier operations.  This is exactly why constraints belong outside PDE
terms.

### 4.6 Make the control Dirichlet boundary data

This is the important intentionally non-local case.  Let $`\Gamma_c`$ be part
of the Dirichlet boundary and choose a control space and lifting, for example

$$
  U_\Gamma=H^{1/2}(\Gamma_c),\qquad
  L_D:U_\Gamma\longrightarrow H^1(\Omega).
$$

Use a homogeneous unknown and reconstruct the physical state:

$$
  y_\mathrm{phys}=\widehat y+\ell_0+L_Du,\qquad \widehat y\in H^1_{\Gamma_D}(\Omega).
  \quad\text{(10)}
$$

| Component | New input | Downstream effect |
|---|---|---|
| Regions/spaces | controlled and fixed Dirichlet portions; control trace space | state/control layouts and trace requirements |
| Transformation | lifting $`L_D`$, value/JVP/transpose-JVP | physical state used by residual and observation |
| Residual | diffusion applied to reconstructed state | discrete $`A_{\mathrm{ext},h} L_{D,h}`$ coupling, not boundary load $`B_h`$ |
| Observation | reads $`y_\mathrm{phys}`$, not merely $\hat y$ | direct chain-rule contribution through $`L_{D,h}`$ |
| Objective/metric/constraint | boundary-control versions | corresponding control-space operations |

With $`r_h(\widehat y,u)=A_h\widehat y+A_{\mathrm{ext},h}L_{D,h}u-f_h`$, the
reduced covector has the generic form

$$
  j_h'(u)=D_uJ_h-\bigl(A_{\mathrm{ext},h}L_{D,h}\bigr)^{T}p_h. \quad\text{(11)}
$$

Here $`D_uJ_h`$ includes any direct dependence of the objective on the
reconstructed physical state.  Equation (11) is deliberately left as a
composed transpose action: turning it into a simplified boundary formula
requires the selected lifting and additional analysis.  The framework remains
generic precisely by retaining this chain rule.

### 4.7 Estimate a diffusion coefficient instead

Promote a coefficient from immutable data to a parameter variable $m$:

$$
  \langle E(y,m),v\rangle=\int_\Omega m\nabla y\cdot\nabla v-\int_\Omega fv,
  \qquad m\in M_\mathrm{ad}. \quad\text{(12)}
$$

The only new PDE-side object is a parameter residual derivative:

$$
  \langle D_mE(y,m)\delta m,v\rangle
    =\int_\Omega\delta m\nabla y\cdot\nabla v. \quad\text{(13)}
$$

| New component | Downstream effect | Solver output needed |
|---|---|---|
| parameter variable $m \in M$ | residual is nonlinear in full $(y,m)$ | nonlinear residual/JVP/transpose-JVP |
| coefficient residual term | state and parameter Jacobian blocks | $`D_mE_h`$ and $`D_mE_h^{*}`$ |
| parameter loss/metric/constraint | parameter derivative/search/feasibility | parameter Riesz and constraint operations |
| product/positivity requirement | model validation/policy only | recorded policy, not an automatic proof |

The reduced and KKT solvers still receive the same generic ports.  They do
not receive a special “inverse problem” type.

### 4.8 Add time dependence

For the heat equation, replace the stationary state space and add a temporal
residual term:

$$
  Y=L^2(0,T;V)\cap H^1(0,T;V^{*}),\quad Z=L^2(0,T;V),
$$

$$
  \langle E(y,u),v\rangle=
  \int_0^{T}\left[\langle\dot y,v\rangle+a(y,v)-(f+u,v)\right]dt,
  \qquad y(0)=y_0. \quad\text{(14)}
$$

| Component change | Downstream effect |
|---|---|
| time-indexed spaces | space-time or per-step layouts for state, test, control, observation |
| time-derivative residual term | causal state Jacobian and transpose/backward adjoint action |
| initial data/trace policy | state residual or affine state reconstruction |
| time-dependent observation/loss | time-integrated or terminal objective derivative |
| temporal discretisation policy | exact discrete time transpose, checkpointing/trajectory access policy |

All spatial terms from the baseline are reused unchanged.  A time-stepping
compiler produces a global or sequential residual $`E_h`$; DTO differentiates
that actual residual.  This is sufficient for a reduced state/adjoint solver
without introducing a parabolic-specific outer optimiser.

## 5. What a solver method needs—and what it does not

| Method | Required executable ports | It must not need |
|---|---|---|
| Reduced gradient / L-BFGS | state solve, adjoint solve or transpose action, reduced covector, metric inverse, line search, constraint projection if needed | PDE family, boundary-condition type, observation type |
| Newton/SQP in reduced space | above plus derivative of reduced covector/Hessian action or incremental state-adjoint actions | hard-coded Laplace/transport branches |
| All-at-once Newton/Krylov | primal/adjoint/stationarity residuals and block JVP/transpose/preconditioner hooks | a separate class per control placement |
| Primal-dual active set | KKT residual plus box/complementarity operations | details of residual terms beyond generic blocks |
| Gauss–Newton inverse solver | observation derivative/transpose, state/incremental-state actions, parameter metric/constraints | “coefficient identification” solver branch |

The framework need not supply every method at first.  It must ensure that
each method can be written against these executable contracts.

## 6. Robustness assessment

### What separates cleanly

- PDE physics: residual terms and their transpose actions.
- Control placement as a source or Neumann flux: a control coupling and its
  transpose.
- Observation location/type: observation and loss, feeding only the adjoint
  right-hand side.
- Objective regularisation: a loss on control/parameter variables.
- Search geometry: a metric, separate from the objective.
- Box constraints: feasibility operations, separate from the PDE.
- Galerkin versus Petrov–Galerkin: a discrete trial/test policy, provided the
  residual formulation itself is unchanged.

### What must cross components

- Dirichlet control: control space, lifting, state reconstruction, residual,
  observation, and chain rule.
- Very-weak formulations: state/test spaces, data pairings, and residual all
  change together.
- Pure Neumann/nullspace cases: residual, gauge/constraint, metric, and
  linear-solver policy must agree.
- Time evolution: time spaces, endpoint data, temporal residual, observations,
  and exact compiled transpose must agree.
- Fractional norms and point/flux observations: space/trace analysis and a
  concrete continuous or discrete-only realization are inseparable.

These are manageable cross-cutting interfaces, not evidence that an
inheritance hierarchy is needed.  They are precisely the reasons that
transformations, policies, and space/pairing declarations must be first-class
components.

## 7. Feasible implementation path exposed by the case study

1. Implement the baseline DTO path: scalar conforming Galerkin, residual
   terms, one observation/loss, $L^2$ metric, and state/adjoint/reduced-gradient
   interfaces.  Verify the three actions by residual, adjoint-consistency, and
   Taylor-remainder tests.
2. Add Neumann control and boundary observation.  This tests region, trace,
   coupling, and observation composition without changing the solver.
3. Add box constraints.  This tests the separation of KKT/optimisation logic
   from the PDE residual.
4. Add fixed and then controlled Dirichlet liftings.  This validates the most
   important transformation and chain-rule boundary.
5. Add coefficient estimation.  This validates nonlinear variable blocks and
   parameter derivative ports.
6. Add the temporal residual adapter and a single time-stepping policy.  This
   validates the actual-discrete-adjoint contract before pursuing hyperbolic
   or space-time variants.
7. Add OTD only as an explicit formulation-builder path, and compare it with
   DTO on cases where equivalence is expected.  Record differences as a
   property of the selected formulations, never as a solver surprise.

The baseline and its extensions require a modest set of generic interfaces.
The difficult cases are visible early, yet none requires redesigning the
residual/objective/metric/constraint composition model.
