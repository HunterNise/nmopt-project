# Chapter 6 numerical-methods implementation guide

## Purpose and scope

This guide turns Chapter 6 of *Optimal Control of Partial Differential
Equations* (Manzoni, Quarteroni, Salsa; book pages 167–218) into an
implementation handoff. It catalogues the numerical methods in Sections
6.2–6.4, 6.6, and 6.8, including their discrete assumptions and solver
services. It does not introduce a class per textbook method or PDE.

The numerical experiments in Sections 6.5, 6.7, and 6.9 are deliberately
kept out of this guide. Their PDE data, algorithms, parameters, and reported
quantities are in the separate
[Chapter 6 numerical-examples reference](chapter-6-numerical-examples.md).
They are future reproduction targets, not current acceptance tests.

Sections 6.10 and 6.11, respectively a-priori and a-posteriori error
estimates, are excluded at the user's request. Section 6.12 exercises is also
not an implementation method. References to SQP and nonlinear OCPs belong to
later chapters; this guide records only the capabilities Chapter 6 needs.

Consult the [roadmap handoff](../planning/implementation-roadmap.md#current-handoff-state)
for implemented solver and formulation services and the
[v1 capability table](../implementation/v1/semantic-compiler.md#registered-capabilities)
for exact compiler targets. The required reusable Chapter 6 extensions are
P6.1–P6.5 in the
[implementation roadmap](../planning/implementation-roadmap.md#chapter-6-feature-requests).

## Framework convention and discrete notation

All implementations retain the project-wide Lagrangian convention

$$
\mathcal L_{h}(x_{h},p_{h})=
J_{h}(x_{h})-\langle p_{h},r_{h}(x_{h})\rangle.
$$

For the linear-quadratic scalar model, write the discrete residual and
objective as

```math
\begin{aligned}
r_{h}(y,u) &= A y-f-Bu, \\
J_{h}(y,u) &= \frac{1}{2}(C y-z_{d})^{\mathsf T}M_{z}(C y-z_{d})
 +\frac{\beta}{2}u^{\mathsf T}N u, \\
W &= C^{\mathsf T}M_{z}C.
\end{aligned}
```

Here $A$ is the state operator, $B$ maps control coefficients into the state
residual dual, $C$ is the declared observation, $M_{z}$ is its pairing, and
$N$ is the control-loss pairing. They must remain separate: full-volume,
subdomain, boundary, and point observations do not justify identifying them.

The framework adjoint and reduced covector are

```math
\begin{aligned}
A^{\mathsf T}p &= W y-C^{\mathsf T}M_{z}z_{d}, \\
j_{h}'(u) &= \beta N u+B^{\mathsf T}p.
\end{aligned}
```

Chapter 6 uses the opposite adjoint sign:
$p_{\mathrm{book}}=-p$. Its relations
$A^{\mathsf T}p_{\mathrm{book}}=-M(y-z_{d})$ and
$\beta N u-B^{\mathsf T}p_{\mathrm{book}}=0$ are equivalent to those above.
An all-at-once implementation may use the book's symmetric KKT multiplier
$\lambda=-p$, but the conversion must be explicit at the formulation
boundary; it must not change `ReducedDTO` signs.

A compiler manifest for every Chapter 6 target must record its state, test,
control, adjoint, and observation spaces; trial/test relation; quadrature and
target-data rules; constraint representation; formulation choice; and state,
adjoint, metric-inverse, and KKT tolerances. A numerical method works with
the compiled discrete problem, not merely with the continuous PDE name.

## Method catalogue

### C6.2 — Optimise then discretise and discretise then optimise

The chapter compares two formulation orders:

- **Optimise then discretise (OtD):** derive continuous state, adjoint, and
  control optimality conditions, then discretise every equation.
- **Discretise then optimise (DtO):** lower $r_{h}$ and $J_{h}$, then
  differentiate that exact finite-dimensional problem. The adjoint is the
  residual VJP and the reduced covector is the DTO covector above.

For a Galerkin linear-quadratic problem using the same finite-element space
for state and adjoint, the routes produce the same KKT system. DTO is the
selected default route. It permits different state and control spaces, so a
boundary-control $B$ may be rectangular.

If OtD chooses a different adjoint space, it can yield a non-symmetric reduced
operator that is not the stationarity system of the lowered discrete
objective. It is not a `ReducedDTO` result: it needs a distinct OTD
formulation product with its own residual, transpose, and comparison
diagnostics.

#### Dominant-advection stabilisation

Section 6.2.4 treats

```math
\begin{aligned}
E y &= -\mathrm{div}(\nu\nabla y)+b\mathbin\cdot\nabla y=u
 && \text{in }\Omega, \\
y &= 0 && \text{on }\Gamma_{D}, \\
\nu\partial_{n}y &= 0 && \text{on }\Gamma_{N}.
\end{aligned}
```

The outflow is $\Gamma_{N}$, $\mathrm{div} b\leq0$, and the local Péclet
number is large. The observation is $k y$ on $\Omega_{0}$ and the control is
penalised in $L^{2}(\Omega)$. The Galerkin–least-squares (GLS) state term is

```math
s_{h}(y_{h},\phi_{h})=
\sum_{K\in\mathcal T_{h}}
\int_{K}\delta_{K}\bigl(Ey_{h}-u_{h}\bigr)E\phi_{h}.
```

The source distinguishes three results that an implementation must not
conflate:

1. Stabilising the continuous state and adjoint separately before
   discretisation gives strongly consistent equations, but not generally the
   stationarity system of a discrete objective.
2. Differentiating the GLS-stabilised state residual gives the exact DTO
   adjoint and control covector. The adjoint and stationarity relation acquire
   all GLS VJP contributions, so they need not be strongly consistent with
   the continuous optimality system.
3. The source's consistent stabilised Lagrangian augments the Lagrangian by
   an element sum of state residual times adjoint residual. Its differentiated
   equations retain strong consistency but contain explicit state–adjoint and
   control–adjoint stabilisation couplings.

The first vertical slice is item 2: model every stabilisation contribution in
the discrete residual and verify its JVP/VJP pairing. Items 1 and 3 require
formulation provenance and are P6.2 work. A state-solver option named
“GLS” is insufficient because it hides whether the discrete objective was
differentiated.

#### Stokes discretisation

Section 6.2.5 uses the steady Stokes state $Y=(v,\pi)$ with distributed vector
control $u$ and a stable velocity-pressure pair:

```math
A_{S}Y=f_{S}+V u,
\qquad
A_{S}=\begin{bmatrix}E&B^{\mathsf T}\\B&0\end{bmatrix}.
```

The source's all-at-once preconditioner assumes a positive pressure tracking
term $\delta M_{\pi}$. A reduced target may track velocity only, but that
specific full-KKT Schur complement then becomes singular. The framework needs
multiple residual/test blocks, a pressure gauge, and an inf-sup policy before
this compiles; those are P5.6 prerequisites extended by P6.3–P6.4.

### C6.3 — Reduced-space methods for unconstrained OCPs

Eliminate the state through the declared state solve, $y=y(u)$, and minimise
$j_{h}(u)=J_{h}(y(u),u)$. Every reduced iteration must:

1. solve the state for the accepted control;
2. solve the DTO adjoint and assemble the reduced covector;
3. apply the declared inverse control metric to obtain a primal gradient;
4. construct a primal direction, select a step, and evaluate trials through
   fresh state solves; and
5. report accepted state, adjoint, covector, objective, step, and solve
   counts.

The book writes directions in Euclidean coefficient coordinates. In this
framework, $g=G^{-1}j_{h}'(u)$ is the gradient for the selected metric $G$.
The portable steepest direction is therefore $d=-g$, not bare coefficient
negation.

The chapter introduction also names trust-region methods. The selected scalar
slice now includes an unconstrained matrix-free Cauchy trust-region solver
with explicit reduction-ratio diagnostics. It remains a separate reduced-
Hessian globalization boundary, not an alternate PDE formulation or a
line-search option.

| Method | Required direction service | First general policy |
| --- | --- | --- |
| Steepest descent | $d=-G^{-1}j_{h}'$ | Existing reduced Armijo solver. |
| Nonlinear conjugate gradient | Gradient history and selected Fletcher–Reeves or Polak–Ribière update | The selected slice defaults to metric-aware PR+, exposes Fletcher–Reeves, verifies exact-search quadratic equivalence, and includes a strict classical quadratic-CG policy. |
| Trust region | Quadratic model, Hessian-vector action, and radius update | Unconstrained metric Cauchy step with actual/predicted reduction diagnostics. |
| Newton / truncated Newton | Hessian-vector action and inner linear solve | The selected slice uses capability-gated Newton; explicit truncated-Newton termination remains an extension. |
| BFGS / L-BFGS | Secant history and metric-aware pairings | Start with limited memory; declare memory, curvature test, reset, and initial inverse-metric policy. |

For a linear-quadratic target, the reduced Hessian is constant:

```math
H=\beta N+B^{\mathsf T}A^{-\mathsf T}W A^{-1}B.
```

Its action on a direction $w$ requires no explicit $H$: solve
$A\delta y=Bw$, solve $A^{\mathsf T}\delta p=W\delta y$, and return
$\beta Nw+B^{\mathsf T}\delta p$. This is the matrix-free Hessian action for
CG/Newton and the source's reduced-Hessian system.

For a positive-curvature linear-quadratic direction, exact line search is

$$
\tau=\frac{-j_{h}'(u)[d]}{d^{\mathsf T}H d}.
$$

For general targets, use declared Wolfe conditions or backtracking Armijo:

$$
j_{h}(u+\tau d)\leq
j_{h}(u)+\sigma\tau j_{h}'(u)[d].
$$

The weak Wolfe curvature condition is

$$
j_{h}'(u+\tau d)[d]\geq c_{2}j_{h}'(u)[d]
$$

The strong Wolfe variant
instead uses the absolute value of the trial slope. Both are selectable
policies in the reduced solver, and both pair slopes with the actual returned
displacement.

Reject a non-descent direction before line search. The reduced result records
the final accepted state, adjoint, covector, objective, absolute and relative
gradient histories, actual metric step norms, objective changes, maximum
iteration/line-search failures, and state/adjoint/metric/Hessian-solve counts
as separate fields. Relative-gradient, objective-change, and step tolerances
are independently configurable and disabled by zero values. For didactical
comparisons, `ReducedStoppingCriterion` can explicitly select one of these
criteria or the absolute gradient norm; its default automatic mode preserves
the legacy behavior in which positive optional tolerances are additional
stops. Finite differences are derivative checks, not a high-dimensional
gradient fallback.

The selected slice does not exhaust the alternatives in Chapter 3. The
relative/objective/step stopping policies and Fletcher–Reeves are now
available alongside the default PR+ policy. The exact-search PR+/Fletcher–Reeves
equivalence, strict classical quadratic-CG recurrence, and unconstrained
matrix-free Cauchy trust-region globalization are verified on the
linear-quadratic target. Generic nonlinear second-order actions and projected
nonlinear-CG/L-BFGS or trust-region directions are lower-priority tracks that
require stronger derivative or active-set contracts; they are not inferred
from the current first-order or projection interfaces.

The line-search boundary is typed around a trial-control builder and evaluator.
Armijo, exact-quadratic, weak-Wolfe, and strong-Wolfe policies all evaluate
acceptance with the actual returned displacement, so a projected trial cannot
reuse the nominal unprojected direction in its pairing. Exact search requires
an explicit
positive-curvature reduced-Hessian capability; first-order DTO ports do not
provide that capability implicitly.

### C6.4 — Projected reduced-space methods for box controls

For $a\leq u\leq b$, use the same covector and a metric-qualified projection:

```math
u^{+}=\Pi_{U_{\mathrm{ad}}}(u+\tau d),
\qquad
r_{\Pi}=u-\Pi_{U_{\mathrm{ad}}}\bigl(u-G^{-1}j_{h}'(u)\bigr).
```

Measure the projected residual in the declared metric. Projected Armijo must
use the actual displacement $u^{+}-u$, not the unprojected direction. The
book's coefficient clipping is exact only under a selected discretisation and
projection policy. It is exact for cellwise `FE_DGQ(0)` $L^{2}$ boxes;
continuous FE, trace, and $H^{1}$ controls need their own projection contract.

The selected slice implements and verifies projected steepest descent only.
Projected nonlinear CG, projected L-BFGS, and active-set-aware trust-region
steps are possible later extensions, but they require restart, active-set
transition, and metric-coupling tests before they can be enabled.

### C6.6 — All-at-once methods for unconstrained OCPs

All-at-once methods retain state and control as simultaneous optimisation
variables. With $x=(y,u)$, define

```math
Q=\begin{bmatrix}W&0\\0&\beta N\end{bmatrix},
\qquad
D=\begin{bmatrix}A&-B\end{bmatrix}.
```

For the symmetric book multiplier $\lambda=-p$, the equality-constrained
quadratic-program KKT system is

```math
\begin{bmatrix}Q&D^{\mathsf T}\\D&0\end{bmatrix}
\begin{bmatrix}x\\\lambda\end{bmatrix}
=
\begin{bmatrix}h\\f\end{bmatrix}.
```

This generic operator, not a scalar-PDE KKT class, is the unit to compile.
The compiler must construct its exact action from objective derivatives and
residual JVP/VJP actions; an assembled matrix is a lowerer optimisation. The
solve result must retain conversion from $\lambda$ to the framework adjoint.
The formulation must record that $D$ has the required rank and that $Q$ is
positive definite on $\ker(D)$, or return a diagnostic instead of selecting a
Schur-complement preconditioner that assumes a nonsingular system.

The source calls direct simultaneous solution **full SAND**. It also notes
**reduced SAND** through null-space or range-space KKT factorisations, but
does not specialise an implementation beyond the reduced systems/preconditioners
described below. Represent those as declared KKT block eliminations, never as
separate PDE solver types.

The KKT matrix is symmetric indefinite. Use MINRES only with a symmetric
positive-definite preconditioner, GMRES for nonsymmetric products, and
flexible GMRES when variable inner Krylov work makes the effective
preconditioner non-stationary.

The chapter identifies domain decomposition and multigrid as the two broad
preconditioning families, then develops multigrid-based block
preconditioners. The first framework targets may follow that focus; a domain
decomposition policy is a later approximate-solve service, not a different
KKT formulation.

#### Block preconditioner families

With $S=DQ^{-1}D^{\mathsf T}$, the source proposes these reusable families:

- **Block diagonal:** $P_{d}=\mathrm{diag}(Q,S)$. For full-volume scalar
  control with $B=N=W=M$, the exact Schur block is
  $S=\beta^{-1}M+A M^{-1}A^{\mathsf T}$. Approximate mass inverses with
  lumping, symmetric Gauss–Seidel, or Chebyshev iteration, and approximate
  the state inverse with fixed AMG or geometric-MG cycles. The source's
  small-$\beta$ approximation drops $\beta^{-1}M$ from this block. It is
  mesh-robust under its assumptions but not generally robust as $\beta$
  decreases.
- **Block triangular / Bramble–Pasciak:** use approximate primal and Schur
  blocks and CG in an explicitly declared non-standard inner product. The
  scaling that makes it positive definite is part of the policy.
- **Constraint preconditioner / PPCG:** retain exact constraint blocks in
  $P_{C}=\begin{bmatrix}G&D^{\mathsf T}\\D&0\end{bmatrix}$. PPCG starts from
  a feasible $x_{0}$ satisfying $Dx_{0}=f$ and keeps iterates in that
  equality-constraint manifold.

Boundary control merely makes $B$ rectangular. It preserves the generic
$Q$, $D$, and Schur-complement construction. For Stokes, take
$D_{S}=[A_{S}\ -V]$. The source's nested approximate Stokes solve can use a
stationary Uzawa method; a variable Krylov inner solve instead requires
flexible GMRES outside.

The scalar source also eliminates the unconstrained optimality equation,
giving the state-adjoint system with blocks
$\begin{bmatrix}M&A^{\mathsf T}\\A&-\beta^{-1}M\end{bmatrix}$ and a
factorised block-diagonal preconditioner with
$M+\beta^{1/2}A$ and $\beta^{-1}M+\beta^{-1/2}A^{\mathsf T}$. Treat this as
one declared KKT block-elimination/preconditioner policy, not a distinct
application formulation.

### C6.8 — Primal-dual active sets and regularised mixed constraints

PDAS is the chapter's all-at-once method for box controls. Let $\mu$ be the
multiplier for $a\leq u\leq b$, represented in a declared control dual. With
the framework adjoint sign, stationarity is

$$
j_{h}'(u)+\mu=0.
$$

For a positive PDAS parameter $c$, classify sets under a selected discrete
evaluation policy:

```math
\begin{aligned}
\mathcal A^{+}&=\{\mu+c(u-b)>0\}, \\
\mathcal A^{-}&=\{\mu+c(u-a)<0\}, \\
\mathcal I&=\{\mu+c(u-b)\leq0\leq\mu+c(u-a)\}.
\end{aligned}
```

At every iteration solve the coupled KKT system with upper/lower values on
active sets and zero multiplier on the inactive set; reclassify; and stop only
when all sets are unchanged and the declared KKT residual is small. PDAS is a
regularised semismooth Newton method. Stable sets do not replace feasibility,
duality, complementarity, and stationarity checks.

The source's nodal formula assumes a particular representation of control and
multiplier coefficients. The framework needs typed primal control, dual
multiplier, an explicit primal/dual conversion for classification, and
selection/restriction operators. It may begin with exact cellwise boxes; it
must not classify a continuous-FE dual vector pointwise without a declared
Riesz or interpolation policy.

The source also applies PDAS to the Lavrentiev-regularised mixed constraint

$$
y_{a}\leq y+\varepsilon u\leq y_{b},
\qquad \varepsilon>0.
$$

It is the observation $O_{c}(y,u)=y+\varepsilon u$, not a control box. The
active equalities are $O_{c}=y_{a}$ or $O_{c}=y_{b}$ and stationarity gains
the $\varepsilon\mu$ covector. It reuses P5.5 semantics but needs P6.5 for
the generic PDAS solve.

## Implementation sequence and verification

1. Completed: extend the reduced solver to selected reduced directions and
   line searches (P6.1), while retaining the baseline one-state/one-decision
   DTO boundary.
2. Next: add formulation, trial/test, and stabilisation provenance (P6.2). A
   compiler must say whether it differentiates a discrete residual, lowers an
   OTD system, or differentiates a stabilised Lagrangian.
3. Generalise executable algebra to equality-constrained quadratic KKT
   products (P6.3), then add mixed Stokes blocks through P5.6.
4. Add block preconditioner and Krylov policies (P6.4) as deterministic
   operator actions before an assembly-only hierarchy.
5. Add complementarity, selection, and PDAS (P6.5), then build the
   regularised mixed constraint.

| Capability | Required verification |
| --- | --- |
| DTO residual/objective | Value, JVP Taylor, VJP pairing, objective derivative, state residual, and reduced Taylor tests. |
| Reduced optimisation | Directional descent, line-search acceptance, metric pairing, stop reason, and solve-count tests. |
| Hessian action | Symmetry in the declared pairing and finite-difference reduced-covector comparison. |
| KKT product | Block action/transpose signs, feasibility, stationarity, multiplier conversion, and agreement with DTO on a small LQ case. |
| Preconditioner | Layout compatibility, deterministic diagnostics, and mesh/parameter sweeps separated from correctness tests. |
| PDAS | Active-set update, feasibility, complementarity, stationarity, stable-set termination, and reduced/KKT agreement on an inactive-box case. |
| Stabilisation | Exact residual VJP plus explicit OtD/DtO/stabilised formulation provenance and comparison tests. |

The example reference provides later end-to-end settings. Its iteration
counts and timings are performance observations, not portable unit-test
tolerances.
