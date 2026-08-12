# Chapter 5 elliptic optimal-control implementation guide

## Purpose

This guide turns Chapter 5 of *Optimal Control of Partial Differential
Equations* (Manzoni, Quarteroni, Salsa; book pages 103–162) into declared
formulations for this project. It is a problem catalogue and implementation
handoff, not a claim that every continuous problem has a currently registered
deal.II compiler target.

The source chapter contains eleven application families. Sections 5.3 and 5.4
are shared theory: the first is the abstract linear-quadratic template and
the second gives the scalar elliptic variational forms used by the application
families. Section 5.9 supplies the transposition formulation needed by the
low-regularity cases; it is not an optimisation problem by itself.

Every implementation must preserve the global convention

$$
\mathcal L(x,p)=J(x)-\langle p,E(x)\rangle.
$$

For a state-control split, the reduced covector is therefore

$$
j'(u)=J_{u}'(y,u)-E_{u}'(y,u)^{\ast}p.
$$

The symbol “gradient” in the source is always qualified here by its chosen
control metric. An $L^{2}$, $H^{1}$, fractional, or $H^{-1}$ gradient is not
interchangeable with the reduced covector.

## How to read a catalogue entry

Each entry identifies the following semantic components.

| Item | Implementation meaning |
| --- | --- |
| State equation | One or more equation blocks, test spaces, residual terms, and state-solve policy. |
| Control | A variable block, its declared primal space, and its constraint if present. |
| Observation and loss | A map from the physical field to an explicit observation space, followed by a quadratic loss. |
| Metric | A separate Riesz map used to turn a control covector into a search direction. |
| Transformation | A declared physical-field map, required for fixed or controlled Dirichlet data. |
| Requirements | Trace, regularity, nullspace, compatibility, inf-sup, or discretisation policies that the compiler must record rather than infer. |

Every differentiable addition must have a residual value test, a JVP Taylor
test, a VJP pairing test, an objective-derivative test, and a state-recomputed
reduced Taylor test. A constraint or metric additionally needs its own
feasibility/projection or apply/inverse-apply test.

## Capability routing and stable feature deltas

Consult the [v1 capability table](../implementation/v1/semantic-compiler.md#registered-capabilities)
for exact compiler support and the
[implementation roadmap](../planning/implementation-roadmap.md#chapter-5-feature-requests)
for completion state. The table below records stable mathematical and
architectural deltas required by the source families; it is not a second
release-status table.

| Source family | Required distinct capability |
| --- | --- |
| 5.1–5.2 volume heat source | Advection residual and transpose actions; any requested continuous-control box needs its own discrete constraint policy. |
| 5.5.1 Robin volume problem | General scalar elliptic coefficients plus Robin boundary residual and source terms. |
| 5.5.2 energy tracking | $H^{1}$ state observation and a separately declared $H^{-1}$ control metric. |
| 5.6 Neumann control with volume tracking | Composition of a Neumann control residual with a volume observation. |
| 5.7 Neumann control with boundary tracking | Weighted trace observation in addition to the unweighted trace case. |
| 5.8 flux observation | Normal-flux observation and its strong or very-weak adjoint policy. |
| 5.10 point observations | Point-set regions, sensor evaluation, and Dirac or very-weak adjoints. |
| 5.11 Dirichlet control | General partial or mixed controlled-boundary transformations and the requested trace metric. |
| 5.12 state constraints | State-derived constraints and their regularised KKT policy. |
| 5.13 Stokes control | Mixed vector state/test blocks, pressure gauge, and saddle-system policies. |

The required reusable extensions are recorded as P5.1–P5.6 in the
[implementation roadmap](../planning/implementation-roadmap.md#chapter-5-feature-requests).

## Shared mathematical declarations

### Abstract linear-quadratic template — source Section 5.3

Declare Hilbert spaces $V$, $U$, and $Z$, a closed convex admissible set
$U_{\mathrm{ad}}\subseteq U$, a continuous coercive bilinear form $a$, a
control operator $B:U\rightarrow V^{\ast}$, and a continuous observation
operator $C:V\rightarrow Z$. The state and objective are

```math
\begin{aligned}
a(y,\phi) &= \langle F+B u,\phi\rangle_{V^{\ast},V}
  && \forall \phi\in V, \\
J(y,u) &= \frac{1}{2}\lVert C y-z_{d}\rVert_{Z}^{2}
  + \frac{\beta}{2}\lVert u\rVert_{U}^{2}.
\end{aligned}
```

The adjoint and minimum principle are

```math
\begin{aligned}
a(\psi,p) &= (C y-z_{d},C\psi)_{Z}
  && \forall \psi\in V, \\
(\beta u,v-u)_{U}+\langle B(v-u),p\rangle_{V^{\ast},V} &\geq 0
  && \forall v\in U_{\mathrm{ad}}.
\end{aligned}
```

The first line is an equation block. The second is obtained from the
residual VJP, the observation VJP, and the declared control loss. It is not a
separate PDE-specific class. For $\beta>0$, the source assumptions guarantee
a unique optimal control. For $\beta=0$, existence requires a bounded
admissible set, and uniqueness additionally requires the relevant
control-to-observation map to be injective.

### Scalar elliptic term palette — source Section 5.4

The chapter's general scalar operator is

```math
E y=-\mathrm{div}(A\nabla y)+\mathrm{div}(b y)+c\mathbin\cdot\nabla y+r y.
```

Its weak volume action is

```math
a(y,\phi)=
(A\nabla y,\nabla\phi)_\Omega
-(y,b\mathbin\cdot\nabla\phi)_\Omega
+(c\mathbin\cdot\nabla y,\phi)_\Omega
+(r y,\phi)_\Omega.
```

The conormal flux is $(A\nabla y-b y)\mathbin\cdot n$. Dirichlet data are
essential, so they require a state restriction or physical-field
transformation. Neumann data/control are boundary functionals on a declared
trace. Robin conditions contribute both a boundary bilinear term and a
boundary load. A mixed boundary partition must name each disjoint region and
the corresponding interpretation.

The source assumes uniform ellipticity and bounded coefficients. Individual
formulations state their coercivity, trace, regularity, or compatibility
requirements. The semantic model records these as requirement policies; the
validator must not represent their presence as a proof of well-posedness.

## Application catalogue

### C5.1 — Unconstrained distributed heat source

The problem is full-volume temperature tracking with distributed source
control:

```math
\begin{aligned}
E y &= -\Delta y+b\mathbin\cdot\nabla y-f-u=0 && \text{in }\Omega, \\
y &= 0 && \text{on }\Gamma, \\
J(y,u) &= \frac{1}{2}\lVert y-z_{d}\rVert_{L^{2}(\Omega)}^{2}
  +\frac{\beta}{2}\lVert u\rVert_{L^{2}(\Omega)}^{2},
  \qquad \beta>0.
\end{aligned}
```

- **Spaces and components:** state/test and adjoint are $H^{1}_{0}(\Omega)$;
  control is $L^{2}(\Omega)$; the control residual term is a volume source;
  the observation is full-volume state restriction.
- **Data and analysis policy:** $f,z_{d}\in L^{2}(\Omega)$; the book assumes
  $b\in C^{1}$ and $\mathrm{div} b=0$, giving coercivity of the stated
  form. The adjoint strong operator is $-\Delta p-b\mathbin\cdot\nabla p$.
- **Optimality:** $p+\beta u=0$ in $L^{2}(\Omega)$. In the framework this is
  the zero reduced covector under the selected $L^{2}$ metric.
- **Framework delta:** the transport term needs P5.1. Do not encode transport
  by modifying the diffusion-reaction lowerer in place.

### C5.2 — Box-constrained distributed heat source

This retains C5.1 but replaces the control domain by

```math
U_{\mathrm{ad}}=
\{u\in L^{2}(\Omega):u_{a}\leq u\leq u_{b}\ \text{a.e. in }\Omega\}.
```

The bounds may be constants or $L^{\infty}$ data. The book allows
$\beta\geq0$. The state and adjoint equations are unchanged, while the
control relation is

```math
(p+\beta u,v-u)_{L^{2}(\Omega)}\geq0
\qquad\forall v\in U_{\mathrm{ad}}.
```

For positive $\beta$, the declared $L^{2}$ projection is pointwise clipping
of $-p/\beta$. For zero $\beta$, use a normal-cone/KKT representation; a
bang-bang control can occur where the adjoint is nonzero. The existing
cellwise box is an exact discrete realisation only for `FE_DGQ(0)` controls
and the `l2_cellwise` metric. Continuous control boxes are a different
feature.

### C5.5.1 — Distributed control and subdomain tracking with Robin data

The state uses a reaction-diffusion residual and a Robin boundary term:

```math
\begin{aligned}
-\mathrm{div}(A\nabla y)+r y &= f+u && \text{in }\Omega, \\
\partial_{n_{E}}y+h y &= g && \text{on }\Gamma, \\
J(y,u) &= \frac{1}{2}\lVert y-z_{d}\rVert_{L^{2}(\Omega_{0})}^{2}
 +\frac{\beta}{2}\lVert u\rVert_{L^{2}(\Omega)}^{2}.
\end{aligned}
```

- **Spaces:** $V=H^{1}(\Omega)$, $U=L^{2}(\Omega)$, and observation space
  $L^{2}(\Omega_{0})$, where $\Omega_{0}$ is a declared volume region.
- **Data and analysis policy:** $f\in L^{2}(\Omega)$,
  $g\in L^{2}(\Gamma)$, $z_{d}\in L^{2}(\Omega_{0})$, and
  $0\leq h\leq h_{0}$ with $h\in L^{\infty}(\Gamma)$. The source's first
  coercive case assumes $r\geq c_{0}>0$.
- **Residual decomposition:** volume diffusion, volume reaction, volume
  source/control, Robin bilinear boundary term, and Robin boundary source.
- **Adjoint:** it has the adjoint volume form, zero Robin datum, and source
  $(y-z_{d})\chi_{\Omega_{0}}$. The control covector is $p+\beta u$.
- **Variants:** no constraint, nonnegative controls, and two-sided boxes. If
  $\beta=0$, the bounded control set gives existence but the subdomain
  observation need not identify a unique control.
- **Framework delta:** the Robin and general-coefficient residual terms need
  P5.1.

### C5.5.2 — Energy tracking with distributed control

The state is $-\Delta y=f+u$ with homogeneous Dirichlet data. The loss is
energy tracking without an explicit control penalty:

```math
J(y,u)=\frac{1}{2}\lVert y-z_{d}\rVert_{H^{1}_{0}(\Omega)}^{2}.
```

The primary variant has $f\in L^{2}(\Omega)$ and
$z_{d}\in H^{1}_{0}(\Omega)$, and uses bounded, closed, convex controls in
$L^{2}(\Omega)$. The adjoint source is the $H^{1}_{0}$ Riesz derivative,
$(I-\Delta)(y-z_{d})$ in $H^{-1}(\Omega)$, and the $L^{2}$ control gradient is
$p$. For an $L^{2}$ ball, the variational inequality has the standard radial
projection solution.

The alternate variant instead takes control and data in $H^{-1}(\Omega)$.
It makes the same quadratic form coercive even for an unbounded admissible
set, but the $H^{-1}$ gradient is $-\Delta p$. This changes the control
metric/Riesz map, not the residual. P5.2 adds the state energy observation and
gives the $H^{-1}$ metric its own selected realization.

The registered discrete realization uses independent homogeneous-Dirichlet
continuous `FE_Q` control coordinates $P_h$ and
$G_h=M_hK_h^{-1}M_h$, with $K_h$ the Dirichlet Laplacian. Its inverse is
$M_h^{-1}K_hM_h^{-1}$, so this is not the existing $H^{1}$ Sobolev-gradient
solve. The fixed boundary removes the constant mode and the compiler records
identity-preconditioned CG tolerances for every inverse action. An $L^{2}$
companion factory uses the same graph and control coordinates. Their compiled
reduced objectives and covectors agree while their search directions differ.
The bounded v1 comparison retains its positive $L^{2}$ control loss and does
not claim to lower the catalogue's unregularized admissible-set variant.

### C5.6 — Subdomain tracking with Neumann boundary control

This is a mixed-boundary convection-diffusion problem:

```math
\begin{aligned}
-\Delta y+\mathrm{div}(b y)&=f && \text{in }\Omega, \\
y&=0 && \text{on }\Gamma_{D}, \\
\partial_{n} y-(b\mathbin\cdot n)y&=u && \text{on }\Gamma_{N}, \\
J(y,u)&=\frac{1}{2}\lVert y-z_{d}\rVert_{L^{2}(\Omega_{0})}^{2}
+\frac{\beta}{2}\lVert u\rVert_{L^{2}(\Gamma_{N})}^{2}.
\end{aligned}
```

Declare separate volume, fixed-Dirichlet, Neumann-control, and observation
regions. State/test space is $H^{1}_{\Gamma_{D}}(\Omega)$; control is
$L^{2}(\Gamma_{N})$. The control residual is a Neumann trace pairing and its
VJP is the adjoint trace. The boundary reduced covector is
$p|_{\Gamma_{N}}+\beta u$. The source assumes a coercivity policy involving
$b\mathbin\cdot n$ and $\mathrm{div} b$ : $b$ has Lipschitz components,
$b\mathbin\cdot n\leq0$ on $\Gamma_{N}$, and
$\mathrm{div} b\geq0$ in $\Omega$. It takes
$f\in L^{2}(\Omega)$ and $z_{d}\in L^{2}(\Omega_{0})$.

This source combination requires composition of volume tracking with boundary
Neumann control, plus P5.1 transport/mixed-boundary support.

### C5.7 — Boundary tracking with Neumann boundary control

The state is

```math
-\Delta y+y=f\quad\text{in }\Omega,
\qquad \partial_{n} y=u\quad\text{on }\Gamma.
```

Control is $L^{2}(\Gamma)$ and the observation is the weighted trace
$h\gamma y$ in $L^{2}(\Gamma)$. The objective is boundary tracking plus
$\frac{\beta}{2}\lVert u\rVert_{L^{2}(\Gamma)}^{2}$. The adjoint has no
volume source and Neumann boundary datum $h(h\gamma y-z_{d})$; the gradient is
$\gamma p+\beta u$. With zero regularisation, a bounded admissible set gives
existence and nonzero $h$ supplies the source's uniqueness condition. Require
$f\in L^{2}(\Omega)$, $z_{d}\in L^{2}(\Gamma)$, and
$h\in L^{\infty}(\Gamma)$; the reaction term $+y$ supplies coercivity.

The registered P5.2 target declares multiplication by boundary data $h$ as an
explicit `weighted_boundary_trace` observation with its own immutable data
port and face-quadrature policy. It reuses the Neumann residual and facewise
$L^{2}$ control metric unchanged; the tracking pullback realizes the boundary
datum $h(h\gamma y-z_d)$.

### C5.8 — Normal-flux tracking with distributed control

The problem tracks a normal derivative on a boundary subset:

```math
\begin{aligned}
-\Delta y&=f+u && \text{in }\Omega, \\
y&=0 && \text{on }\Gamma, \\
J(y,u)&=\frac{1}{2}\lVert \partial_{n} y-z_{d}\rVert_{L^{2}(\Gamma_{0})}^{2}
+\frac{\beta}{2}\lVert u\rVert_{L^{2}(\Omega)}^{2}.
\end{aligned}
```

It requires a convex or sufficiently smooth domain so the state belongs to
$H^{2}(\Omega)\cap H^{1}_{0}(\Omega)$ and its normal trace is meaningful.
The formal adjoint is harmonic but has $L^{2}$ Dirichlet data on the boundary,
so its correct formulation is a very weak solution in $L^{2}(\Omega)$.
The volume reduced covector remains $p+\beta u$. Require
$f\in L^{2}(\Omega)$, $z_{d}\in L^{2}(\Gamma_{0})$, a smooth boundary subset
$\Gamma_{0}$, and $\beta>0$.

An implementation must select either a strong-state normal-flux policy or an
$H(\mathrm{div})$ trace policy, and must lower the adjoint through the
declared transposition policy. P5.3 is required; treating a normal derivative
as an ordinary boundary trace is incorrect.

The first registered C5.8 slice selects the strong-state policy
$Y=H^{2}(\Omega)\cap H^{1}_{0}(\Omega)$ on the declared convex or sufficiently
smooth domain. Its deal.II realization evaluates
$\partial_{n}y=\nabla y\mathbin\cdot n_{\mathrm{out}}$ with `FEFaceValues` at
selected boundary-face quadrature points and assembles the transpose of the
same face map. The objective uses this face-quadrature $L^{2}$ pairing and the
adjoint is recorded as the corresponding very-weak boundary source. An
$H(\mathrm{div})$ realization and alternate flux policies remain unregistered.

### C5.9 — Transposition policy for low-regularity Dirichlet data

For $-\Delta w=f$ with boundary datum $g\in L^{2}(\Gamma)$, the source uses
the test space $Y=H^{2}(\Omega)\cap H^{1}_{0}(\Omega)$ and defines the very
weak solution $w\in L^{2}(\Omega)$ through equations (5.130)–(5.134):

```math
-(w,\Delta\psi)_\Omega=(f,\psi)_\Omega
-(g,\partial_{n}\psi)_\Gamma
\qquad\forall\psi\in Y.
```

Equivalently, declare the residual $E_{\mathrm{tr}}(w,g;f)\in Y^{\ast}$ by

```math
\langle E_{\mathrm{tr}}(w,g;f),\psi\rangle
=(w,-\Delta\psi)_\Omega-(f,\psi)_\Omega
+(g,\partial_{n}\psi)_\Gamma.
```

This is a formulation choice. Proposition 5.11 applies transposition with
$T=-\Delta:Y\rightarrow L^{2}(\Omega)$ a continuous isomorphism under the
declared domain-regularity assumptions. The residual therefore lies in
$Y^{\ast}$, while the Lagrangian seed or multiplier paired with it lies in
$Y$ after the Hilbert-space identification. If $g\in H^{1/2}(\Gamma)$, the
very weak solution gains $H^{1}(\Omega)$ regularity and equals the ordinary
variational solution; this equivalence is the important implementation route
used below.

The transposition policy is needed by C5.8, C5.10, and the continuous
$L^{2}$ Dirichlet-control model in C5.11. It must never be inferred from a
standard $H^{1}$ residual. A discrete subspace contained in
$H^{1/2}(\Gamma)$ may nevertheless use the equivalent variational state
solve, provided that subspace and equivalence are declared and controls
outside it are rejected.

### C5.10 — Finite point-sensor tracking

For $d=2,3$ and a convex or $C^{2}$ domain, choose
$Y=H^{2}(\Omega)\cap H^{1}_{0}(\Omega)$ so that point evaluation is
continuous. Take $r\in L^{\infty}(\Omega)$, $r\geq0$,
$f\in L^{2}(\Omega)$, and $z_{d}\in C(\overline\Omega)$.
The problem is

```math
\begin{aligned}
-\Delta y+r y&=f+u, \qquad y=0, \\
J(y,u)&=\frac{1}{2}\sum_{j=1}^{N}
  (y(\xi_{j})-z_{d}(\xi_{j}))^{2}
+\frac{\beta}{2}\lVert u\rVert_{L^{2}(\Omega)}^{2},
\qquad \beta>0.
\end{aligned}
```

The control is a closed-convex subset of $L^{2}(\Omega)$. The adjoint source
is a weighted sum of Dirac measures at the declared point-set region, hence
the adjoint is very weak and belongs to $L^{2}(\Omega)$ in the source's
setting. The reduced covector is $p+\beta u$.

P5.3 must define sensor locations, their discrete evaluation rule, and the
Dirac/transpose action. In particular, neither nearest-node evaluation nor
quadrature-point coincidence is an acceptable undeclared policy.

The first registered point-sensor slice selects the finite point-evaluation
alternative. Its semantic
point-set region stores finite, unique physical coordinates and its observation
space has one coefficient per coordinate. The continuous policy declares
$Y=H^{2}(\Omega)\cap H^{1}_{0}(\Omega)$ and
$T=-\Delta+rI:Y\rightarrow L^{2}(\Omega)$ under a convex-or-$C^{2}$ domain
assumption. The deal.II realization evaluates `FE_Q` shape functions at each
physical coordinate to form the finite-dimensional map $C_{h}$ and assembles
$C_{h}^{\mathsf{T}}(C_{h}y-z_{d})$ as the point-load transpose for the
very-weak adjoint solve. This is the selected evaluation and transpose policy;
nearest-node, quadrature-coincidence, and general transposition alternatives
remain unregistered. The sibling C5.8 slice uses the separately declared
strong-state normal-flux policy above.

### C5.11 — Dirichlet boundary control and control-space variants

All variants use

```math
-\Delta y=f\quad\text{in }\Omega,
\qquad y=u\quad\text{on }\Gamma,
```

with an unconstrained control and positive regularisation. Dirichlet control
is a physical-state transformation, not a boundary load. The source takes
$f\in L^{2}(\Omega)$ and $z_{d}$ in the selected observation space $Z$:

```math
y_{\mathrm{phys}}=P\widehat y+\ell_{0}+L_{D} u.
```

The source's variants are:

1. **Control space $H^{1/2}(\Gamma)$.** With $L^{2}(\Omega)$ state
   tracking and $H^{1/2}(\Gamma)$ control regularisation, use a lifting and
   the fractional Riesz map in the Euler equation. With $L^{2}(\Gamma)$
   regularisation instead, use $H^{1}(\Omega)$ state tracking to retain
   coercivity; the Euler equation is a boundary duality involving
   $\partial_{n} p$.
2. **Control space $L^{2}(\Gamma)$.** The state is initially very weak and
   needs the C5.9 transposition policy. More precisely, Section 5.11.2 takes
   $Z=L^{2}(\Omega)$, $U=U_{0}=L^{2}(\Gamma)$, and

   ```math
   \begin{aligned}
   J(y,u)&=\frac{1}{2}\lVert y-z_{d}\rVert_{L^{2}(\Omega)}^{2}
   +\frac{\beta}{2}\lVert u\rVert_{L^{2}(\Gamma)}^{2}, \\
   \langle E_{\mathrm{tr}}(y,u;f),\psi\rangle
   &=(y,-\Delta\psi)_\Omega-(f,\psi)_\Omega
   +(u,\partial_{n}\psi)_\Gamma=0
   \quad\forall\psi\in H^{2}(\Omega)\cap H^{1}_{0}(\Omega).
   \end{aligned}
   ```

   Here $f,z_{d}\in L^{2}(\Omega)$ and $\beta>0$. With the project convention
   $\mathcal L=J-\langle p,E\rangle$, the adjoint and Euler equations are

   ```math
   -\Delta p=y-z_{d}\quad\text{in }\Omega,
   \qquad p=0\quad\text{on }\Gamma,
   \qquad \beta u-\partial_{n}p=0\quad\text{on }\Gamma.
   ```

   These signs agree with source equations (5.171)–(5.174) and Proposition
   5.16. Remark 5.18 prints a plus sign in its reduced-gradient formula; that
   sign is inconsistent with the preceding state residual and Euler
   equations and is not selected here.

   Transposition well-posedness requires a convex or sufficiently smooth
   domain. Proposition 5.16 uses a $C^{2}$ domain for the stronger conclusion:
   $p\in H^{2}(\Omega)\cap H^{1}_{0}(\Omega)$ implies
   $\partial_{n}p\in H^{1/2}(\Gamma)$, the Euler equation then gives
   $u\in H^{1/2}(\Gamma)$, and the optimal state is the ordinary
   $H^{1}(\Omega)$ variational solution.

   For the first discrete target, choose a conforming trace-control subspace
   $U_{h}\subset H^{1/2}(\Gamma)\subset L^{2}(\Gamma)$ with the boundary
   $L^{2}$ mass metric. For every discrete datum, the associated continuous
   transposition state has the equivalent $H^{1}$ variational formulation, so
   a lifted conforming Galerkin solve is valid at every discrete iterate, not
   only at the discrete optimum. This reuses the declared
   controlled-Dirichlet transformation while retaining the continuous
   $L^{2}(\Gamma)$ parent problem. A discontinuous or merely facewise
   $L^{2}$ Dirichlet control is outside this shortcut and still requires a
   direct transposition lowerer.
3. **Control space $H^{1}(\Gamma)$.** Use tangential-gradient
   regularisation and a surface weak Euler equation. This avoids a
   fractional norm and the very-weak state, but requires a declared surface
   gradient and its metric.

The first discrete realization of item 1 uses the quotient trace norm induced
by the volume $H^{1}$ energy. For $U_{h}=\mathrm{tr}_{\Gamma}V_{h}$, let
$E_{h}u_{h}$ minimize the discrete $H^{1}(\Omega)$ norm among functions with
trace $u_{h}$ and set

```math
\langle G_{1/2,h}u_{h},v_{h}\rangle
=(E_{h}u_{h},E_{h}v_{h})_{H^{1}(\Omega)}.
```

Thus $G_{1/2,h}$ is the Schur complement of
$M_{\Omega,h}+K_{\Omega,h}$ with respect to the interior coordinates. This
choice supplies the fractional Riesz map required by Section 5.11.1 without
introducing a separate boundary eigendecomposition. Its two source options
remain distinct: the first uses the fractional action for both control loss
and metric with $L^{2}$ state tracking; the second uses an $L^{2}$ control
loss, $H^{1}$ state tracking, and the fractional action only for direction
formation.

The first discrete realization of item 3 uses
$G_{1,h}=M_{\Gamma,h}+K_{\Gamma,h}$, with $K_{\Gamma,h}$ assembled from
tangential shape gradients. It supplies both the selected tangential
$H^{1}(\Gamma)$ control loss and search metric while preserving their
separate semantic identities.

The exact registered Dirichlet-control boundary is recorded in the
[v1 capability table](../implementation/v1/semantic-compiler.md#registered-capabilities).
The complete-boundary factories
`make_hhalf_dirichlet_laplace_control_problem()`,
`make_h1_tracking_hhalf_dirichlet_laplace_control_problem()`, and
`make_h1_dirichlet_laplace_control_problem()` implement items 1 and 3 with
the choices above. Their common contract verifies the reduced stationarity
composition, including that item 1's second option retains an $L^{2}$ control
loss while using the fractional action only as its search metric. P5.4 still
owns alternate partial/nonmatching trace realizations. P5.3 and a general
transposition lowerer remain necessary for a nonconforming $L^{2}$-control
variant.

### C5.12 — State-constrained distributed control

The source problem is on a bounded Lipschitz domain, with
$y_{a},y_{b}\in L^{\infty}(\Omega)$ and a feasible control. It is

```math
\begin{aligned}
-\Delta y+y&=u && \text{in }\Omega, \\
\partial_{n} y&=0 && \text{on }\Gamma, \\
y_{a}\leq y&\leq y_{b} && \text{a.e. in }\Omega, \\
J(y,u)&=\frac{1}{2}\lVert y-z_{d}\rVert_{L^{2}(\Omega)}^{2}
+\frac{\beta}{2}\lVert u\rVert_{L^{2}(\Omega)}^{2},
\qquad \beta>0.
\end{aligned}
```

The dimensions are restricted to two or three so the selected state has a
continuous representative. Original state-constraint multipliers are
measures, which is intentionally outside the current first-order contract.
The source therefore chooses a Lavrentiev regularisation,

```math
y_{a}\leq y+\lambda u\leq y_{b},
\qquad \lambda>0.
```

It introduces the transformed decision $v=(\lambda I+S)u$, where $S$ is the
control-to-state map. This gives ordinary box constraints in $v$ and
$L^{2}$ KKT multipliers. A valid first implementation must compile that
regularised formulation only, record $\lambda$, prove the multiplier and
complementarity identities, and make the original measure-multiplier problem
a diagnostic rather than a silent approximation. P5.5 records the feature.

### C5.13.1 — Distributed velocity control for Stokes flow

The state has velocity and pressure blocks:

```math
\begin{aligned}
-\nu\Delta v+\nabla\pi&=f+u, \\
\mathrm{div} v&=0.
\end{aligned}
```

Velocity has fixed Dirichlet data on $\Gamma_{D}$ and traction data on
$\Gamma_{N}$. Use $X=H^{1}_{\Gamma_{D}}(\Omega)^d$ for velocity,
$Q=L^{2}(\Omega)$ with the selected pressure gauge, and $U=L^{2}(\Omega)^d$
for an unconstrained force control. The objective is

```math
J(v,u)=\frac{1}{2}\lVert v-z_{d}\rVert_{L^{2}(\Omega)^{d}}^{2}
+\frac{\tau}{2}\lVert u\rVert_{L^{2}(\Omega)^{d}}^{2}.
```

The residual consists of momentum and incompressibility equation blocks. Its
well-posedness depends on the declared inf-sup and pressure-gauge policies.
The source takes $f\in L^{2}(\Omega)^{d}$ and traction data in
$L^{2}(\Gamma_{N})^{d}$.
The adjoint is a Stokes system driven by $v-z_{d}$; its velocity component $z$
gives reduced covector $z+\tau u$. P5.6 must extend the formulation and
compiler paths to mixed blocks before this is implemented.

### C5.13.2 — Boundary velocity control with vorticity loss

This two-dimensional Stokes case controls velocity on a subset of a body
boundary. Fixed inflow, wall, and outflow conditions coexist with the
controlled Dirichlet trace. The source selects a normal velocity control in
$H^{1}_{0}(\Gamma_{c})^{2}$ and minimises

```math
J(v,u)=\frac{1}{2}\lVert \nabla\mathbin\times v\rVert_{L^{2}(\Omega)}^{2}
+\frac{\tau}{2}\lVert \nabla_{\Gamma} u\rVert_{L^{2}(\Gamma_{c})^{2}}^{2}.
```

It uses a mixed Stokes formulation with a boundary multiplier for the
Dirichlet relation, rather than an unstated lifting. The adjoint is another
Stokes system with curl-curl forcing. The Euler condition is a surface
Poisson balance between the $H^{1}_{0}(\Gamma_{c})$ regularisation and the
adjoint boundary multiplier.

This needs P5.4 for the trace/surface metric pieces and P5.6 for vector
mixed equations, traction, curl observation, and boundary multiplier blocks.
The implementation must retain all of those as separate components; it must
not add a `StokesVorticityControlProblem` class.

## Implementation sequence for one selected application

For any one catalogue entry, an agent should proceed in this order:

1. Select the exact continuous variant, including every region, state/control
   space, observation space, control metric, admissible set, and source
   regularity assumption. Do not merge variants solely because their strong
   PDE looks similar.
2. State the physical-field transformation before declaring residual or
   observation terms. This is mandatory for every fixed or controlled
   Dirichlet condition.
3. Declare each residual term and its test space. A boundary contribution
   must identify its trace pairing and boundary region. A mixed field must
   use distinct variable and equation blocks.
4. Declare observations independently of losses. In particular, normal flux,
   point sensors, weighted traces, and energy observations require their own
   output spaces and policies.
5. Declare the loss and the algorithmic metric separately. A control
   regularisation changes the objective derivative; a search metric only
   maps the resulting covector to a direction.
6. Supply all requirement policies and compiler bindings, then implement
   value, JVP, and VJP paths before a state/adjoint solve or optimisation
   loop.
7. Add the required contract and reduced-Taylor tests. Add an explicit
   capability diagnostic for every unregistered space, boundary partition,
   trace, or formulation choice.

The roadmap requests below are the only intended route for extending the
registered compiler. They keep the Chapter 5 catalogue compositional
and make the difference between a new general capability and a new application
configuration explicit.
