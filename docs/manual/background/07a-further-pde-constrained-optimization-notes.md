# 07a · Further notes on PDE-constrained optimization

**Background navigation:** [Index](README.md) \
Main chapter: [07 · PDE-constrained optimization](07-pde-constrained-optimization.md) \
Previous: [07 · PDE-constrained optimization](07-pde-constrained-optimization.md) \
Next: [07b · Extensions of PDE-constrained optimization](07b-extensions-of-pde-constrained-optimization.md)

## Purpose

The main [07 · PDE-constrained optimization](07-pde-constrained-optimization.md) chapter develops the linear-quadratic distributed-control setting far enough to explain state equations, sensitivities, adjoints, reduced derivatives, all-at-once systems, and DTO/OTD discretization. Those ideas survive in much broader problems, but some of the first generalizations would interrupt the prerequisite path if they were developed there.

This companion stays close to that path. It asks what changes when the state equation is nonlinear, when second derivatives are needed, when the control acts on the boundary, or when the optimization variable is a coefficient inside the PDE operator. These extensions are especially useful because they clarify which parts of the state–adjoint construction are structural and which parts were consequences of the simple Poisson model.

The treatment is intentionally selective. We derive enough of each extension to make the new mathematical object visible and usable. Full well-posedness theory for nonlinear PDEs, the approximation theory behind low-regularity boundary formulations, global convergence of Newton/SQP methods, state constraints, time-dependent control, and nonsmooth optimization are deferred to the sources or to [07b · Extensions of PDE-constrained optimization](07b-extensions-of-pde-constrained-optimization.md).

## Before you start

You should be comfortable with:

- the main [07 · PDE-constrained optimization](07-pde-constrained-optimization.md) derivation of state sensitivities and adjoints;
- Fréchet derivatives, duality, and adjoints from [02 · Duality, derivatives, and adjoints](02-duality-derivatives-and-adjoints.md);
- Newton and quasi-Newton ideas from [05 · Unconstrained numerical optimization](05-unconstrained-optimization.md);
- KKT systems and the SQP viewpoint from [06 · Constrained optimization](06-constrained-optimization.md) and its companion.

The boundary-control section also uses traces from [01 · Function spaces and weak PDEs](01-function-spaces-and-weak-pdes.md) and the finite-element viewpoint from [03 · Finite elements](03-finite-elements.md).

## What you will be able to do

After these notes, you should be able to:

- explain why a nonlinear state equation usually gives only a local differentiable control-to-state map;
- derive the linearized state equation for a semilinear elliptic problem;
- derive the corresponding nonlinear adjoint equation and reduced derivative;
- construct a reduced Hessian action using an incremental state and incremental adjoint solve;
- see why second derivatives appear when an all-at-once optimality system is linearized;
- relate a Newton step for the optimality system to an SQP subproblem;
- formulate a Neumann boundary control as a bounded control operator into a dual space;
- explain why Dirichlet boundary control is analytically more delicate and formulate it both through an $H^{1/2}(\Gamma)$ lifting and, at a rougher $L^{2}(\Gamma)$ level, through a very-weak transposition formulation;
- derive the sensitivity and adjoint derivative for a coefficient-identification problem;
- formulate subdomain, boundary-trace, point-sensor, and boundary-flux observation operators and identify the regularity each one requires;
- explain why singular or high-order observations can lower adjoint regularity even when the state is smooth;
- distinguish observation operators, regularization, and optimization metrics in inverse problems.

## Roadmap

We first replace the linear Poisson equation by a semilinear one and follow the same chain as in [07 · PDE-constrained optimization](07-pde-constrained-optimization.md): state map, sensitivity, adjoint, and reduced derivative. We then differentiate once more to see how Hessian-vector products arise without assembling a dense reduced Hessian. The same second-order terms reappear when the full state–adjoint–control system is linearized, which gives the bridge to Newton and SQP methods.

The second half changes the way the control enters the PDE. Neumann boundary control turns the trace map into part of the control operator. Dirichlet control shows why essential data are different: at the natural energy level one introduces a lifting, while rougher boundary controls can force a very-weak transposition formulation. Coefficient identification then places the optimization variable inside the PDE operator itself. A final section makes observation operators concrete, including point sensors and boundary-flux measurements, because these examples show particularly clearly how the observation map determines the adjoint's right-hand side and regularity.

## 1. Nonlinear state equations change the control-to-state map

Consider the semilinear elliptic problem

```math
\begin{aligned}
\min_{y,u}\quad
J(y,u)
&:=
\frac{1}{2}\lVert y-y_{d}\rVert_{L^{2}(\Omega)}^{2}
+
\frac{\beta}{2}\lVert u\rVert_{L^{2}(\Omega)}^{2},\\
\text{subject to}\quad
-\Delta y+y^{3}
&=
f+u
\quad\text{in }\Omega,\\
y&=0
\quad\text{on }\partial\Omega,
\end{aligned}
```

with $\beta>0$. Assume, for this model, that the dimension and regularity are such that the cubic term defines a bounded element of $V^{\ast}$; for the standard energy setting, $d\leq3$ is a convenient case. Let

$$
V:=H_{0}^{1}(\Omega),
\qquad
U:=L^{2}(\Omega).
$$

The weak residual is

$$
E:V\times U\to V^{\ast},
$$

defined by

$$
\langle E(y,u),v\rangle
:=
\int_{\Omega}\nabla y\cdot\nabla v\mathrm{d}x
+
\int_{\Omega}y^{3}v\mathrm{d}x
-
\int_{\Omega}(f+u)v\mathrm{d}x.
$$

The equation is nonlinear because of $y^{3}$. Even if the objective is quadratic, the map

$$
u
\mapsto y
$$

need no longer be affine.

### 1.1 The implicit-function viewpoint is local

Suppose

$$
E(\bar y,\bar u)=0.
$$

The derivative of the residual with respect to the state is the linear operator

$$
D_{y}E(\bar y,\bar u):V\to V^{\ast}
$$

with action

$$
\langle D_{y}E(\bar y,\bar u)z,v\rangle
=
\int_{\Omega}\nabla z\cdot\nabla v\mathrm{d}x
+
\int_{\Omega}3\bar y^{2}zv\mathrm{d}x.
$$

The control derivative is simpler:

$$
\langle D_{u}E(\bar y,\bar u)h,v\rangle
=
-
\int_{\Omega}hv\mathrm{d}x.
$$

If $E$ is sufficiently smooth near $(\bar y,\bar u)$ and

$$
D_{y}E(\bar y,\bar u):V\to V^{\ast}
$$

is an isomorphism, the Banach-space implicit-function theorem gives a neighborhood of $\bar u$ in which the state can be written as a differentiable function

$$
y=S(u).
$$

The word **local** matters. For a nonlinear PDE, uniqueness of the state can fail globally even when a particular solution branch is locally well behaved. A reduced method only needs a single-valued state map on the part of control space that it actually visits, but that is now an analytical and numerical assumption rather than an automatic consequence of linearity.

### 1.2 Differentiate the nonlinear state equation

Let

$$
z:=S'(u)h
$$

be the state perturbation induced by a control perturbation $h$. Differentiating

$$
E(S(u),u)=0
$$

gives

$$
D_{y}E(y,u)z
+
D_{u}E(y,u)h
=
0.
$$

For the semilinear example this becomes

$$
\int_{\Omega}\nabla z\cdot\nabla v\mathrm{d}x
+
\int_{\Omega}3y^{2}zv\mathrm{d}x
=
\int_{\Omega}hv\mathrm{d}x
\qquad
\forall v\in V.
$$

Thus the sensitivity is obtained from a **linear** PDE, but the linearized operator depends on the current nonlinear state $y$.

This is one of the recurring patterns in nonlinear PDE optimization:

```text
nonlinear state solve
        ↓
current state y
        ↓
linearized state operator E_y(y,u)
        ↓
sensitivity and adjoint solves
```

The expensive nonlinear solve creates the point at which the linearized optimization calculations are performed.

### 1.3 Nonlinearity can destroy convexity after reduction

The reduced functional is still

$$
j(u)=J(S(u),u).
$$

The displayed objective $J(y,u)$ is convex in the pair $(y,u)$. That does **not** imply that $j$ is convex in $u$ once $S$ is nonlinear.

The composition

$$
u
\mapsto S(u)
\mapsto
\frac{1}{2}\lVert S(u)-y_{d}\rVert^{2}
$$

can introduce several local stationary points. Consequently,

$$
j'(u)=0
$$

becomes only a first-order necessary condition for a local optimum unless additional structure is available. This is one reason second-order conditions and globalization become more important in nonlinear PDE control than in the linear-quadratic model.

## 2. The nonlinear adjoint still removes the state sensitivity

The direct reduced derivative is

```math
j'(u)[h]
=
D_{y}J(y,u)[z]
+
D_{u}J(y,u)[h],
```

where $z=S'(u)h$ solves the linearized state equation.

For the running example,

$$
D_{y}J(y,u)[z]
=
(y-y_{d},z)_{L^{2}(\Omega)},
$$

and

$$
D_{u}J(y,u)[h]
=
\beta(u,h)_{L^{2}(\Omega)}.
$$

Introduce an adjoint $p\in V$ through the same convention as [07 · PDE-constrained optimization](07-pde-constrained-optimization.md):

$$
D_{y}E(y,u)^{\ast}p
=
D_{y}J(y,u).
$$

The state linearization in this example is symmetric under the weak pairing, so the adjoint equation is

$$
\int_{\Omega}\nabla p\cdot\nabla v\mathrm{d}x
+
\int_{\Omega}3y^{2}pv\mathrm{d}x
=
\int_{\Omega}(y-y_{d})v\mathrm{d}x
\qquad
\forall v\in V.
$$

In strong notation, when the regularity is sufficient,

$$
-\Delta p+3y^{2}p
=
y-y_{d}.
$$

The adjoint operator therefore depends on the current state for the same reason as the sensitivity operator.

Using the adjoint identity and the linearized state equation,

$$
D_{y}J(y,u)[z]
=
\langle D_{y}E(y,u)^{\ast}p,z\rangle
=
\langle p,D_{y}E(y,u)z\rangle
=
-
\langle p,D_{u}E(y,u)h\rangle.
$$

Hence

$$
\boxed{
j'(u)
=
D_{u}J(y,u)
-
D_{u}E(y,u)^{\ast}p.
}
$$

Nothing in this formula required linearity of the original PDE. What changed was the point-dependent operator used in the state and adjoint solves.

For the semilinear source-control problem,

$$
D_{u}E=-I,
$$

so

$$
\boxed{
j'(u)=\beta R_{L^{2}}u+R_{L^{2}}p.
}
$$

With the standard $L^{2}$ Riesz identification, the reduced gradient is simply

$$
\nabla j(u)=\beta u+p.
$$

The formula looks exactly like the linear Poisson case, but computing $p$ now requires the linearization at the current nonlinear state.

## 3. A reduced Hessian action needs incremental state and adjoint equations

A Newton method for the reduced problem needs information about

$$
j''(u).
$$

For a PDE problem, explicitly assembling the full reduced Hessian is usually unnecessary. What a Krylov or truncated-Newton method needs is the action

$$
h
\mapsto
j''(u)h.
$$

The semilinear example makes the structure visible without introducing a large amount of abstract notation.

### 3.1 First incremental solve: the state variation

Given $h\in U$, solve for

$$
z=S'(u)h
$$

from

$$
\int_{\Omega}\nabla z\cdot\nabla v\mathrm{d}x
+
\int_{\Omega}3y^{2}zv\mathrm{d}x
=
\int_{\Omega}hv\mathrm{d}x.
$$

This is the same sensitivity equation as before.

### 3.2 Differentiate the adjoint equation

Let

$$
q:=p'(u)h
$$

be the adjoint variation. The adjoint equation is

$$
A_{y}p=y-y_{d},
$$

where

$$
A_{y}v:=-\Delta v+3y^{2}v.
$$

The state changes by $z$, so the derivative of the coefficient $3y^{2}$ is

$$
6yz.
$$

Differentiating the adjoint equation gives

$$
A_{y}q
+
6yzp
=
z.
$$

In weak form,

```math
\int_{\Omega}\nabla q\cdot\nabla v\mathrm{d}x
+
\int_{\Omega}3y^{2}qv\mathrm{d}x
=
\int_{\Omega}\left(z-6yzp\right)v\mathrm{d}x
\qquad
\forall v\in V.
```

This is an **incremental adjoint equation**. It uses the same principal linearized state operator $A_{y}$ but a new right-hand side containing both the state variation and the second derivative of the nonlinear PDE term.

### 3.3 Differentiate the reduced gradient

The $L^{2}$ gradient was

$$
g(u)=\beta u+p.
$$

Therefore

$$
\boxed{
j''(u)h
\quad\leftrightarrow\quad
\beta h+q
}
$$

under the $L^{2}$ Riesz identification.

A reduced Hessian-vector product can therefore be organized as

```text
control direction h
      ↓
incremental state solve for z
      ↓
incremental adjoint solve for q
      ↓
β h + q
```

This is the second-order analogue of the adjoint-gradient idea: the algorithm works with **operator actions and solves** rather than with an explicitly assembled dense Hessian.

The exact equations change for a different nonlinear PDE or objective, but the pattern survives. Second derivatives of $E$ and $J$ enter through the right-hand sides and couplings of the incremental equations.

## 4. Linearizing the full optimality system leads to Newton–KKT systems

Instead of reducing to the control, retain the state, adjoint, and control as independent unknowns. With the global sign convention

$$
\mathcal L(y,u,p)
=
J(y,u)
-
\langle p,E(y,u)\rangle,
$$

the smooth unconstrained first-order system is

```math
F(y,p,u)
:=
\begin{bmatrix}
E(y,u)\\
D_{y}J(y,u)-D_{y}E(y,u)^{\ast}p\\
D_{u}J(y,u)-D_{u}E(y,u)^{\ast}p
\end{bmatrix}
=0.
```

For a nonlinear problem, $F$ is itself nonlinear. Newton's method asks for a correction

$$
(\delta y,\delta p,\delta u)
$$

satisfying

$$
F'(y,p,u)
\begin{bmatrix}
\delta y\\
\delta p\\
\delta u
\end{bmatrix}
=
-F(y,p,u).
$$

The derivative of the first block contains only the first derivatives

$$
D_{y}E,
\qquad
D_{u}E.
$$

The derivative of the adjoint and stationarity blocks differentiates quantities such as

$$
D_{y}E(y,u)^{\ast}p.
$$

That is why **second derivatives of the original PDE residual appear in a Newton linearization of the first-order optimality system**. They are not optional decorations added by an advanced algorithm; they enter because the first-order equations already contain first derivatives.

### 4.1 The SQP interpretation

The same Newton step can be understood from constrained optimization. At the current pair $(y,u)$, linearize the PDE constraint:

$$
E(y,u)
+
D_{y}E(y,u)\delta y
+
D_{u}E(y,u)\delta u
=0.
$$

Then build a quadratic model of the Lagrangian in the primal variables. Schematically,

```math
\frac{1}{2}
\begin{bmatrix}
\delta y\\
\delta u
\end{bmatrix}^{\mathsf T}
H_{\mathcal L}
\begin{bmatrix}
\delta y\\
\delta u
\end{bmatrix}
+
\left\langle
\nabla_{y,u}\mathcal L,
\begin{bmatrix}
\delta y\\
\delta u
\end{bmatrix}
\right\rangle,
```

subject to the linearized state equation.

Solving this equality-constrained quadratic subproblem produces the classical **sequential quadratic programming** viewpoint. Under the usual smoothness and regularity assumptions, exact-Hessian SQP and Newton applied to the KKT equations are two descriptions of the same local second-order mechanism.

A complete SQP method needs more than this local algebra: merit functions or filters, globalization, constraint qualifications, and inexact linear solves all matter in practice. Those topics belong to full nonlinear-programming treatments. The important PDE-control connection is that each SQP/Newton step again exposes a large structured KKT system, so the block linear algebra from [04 · Numerical linear algebra for PDE and optimization systems](04-numerical-linear-algebra.md) remains relevant.

## 5. The control operator can act through the boundary

The distributed source model used

$$
E(y,u)=Ay-Bu-f
$$

with a volume control operator $B$. Nothing in the adjoint derivation requires $B$ to be a volume embedding. A boundary control simply changes what $B$ means and which control space is natural.

### 5.1 Neumann boundary control

Consider the coercive model

```math
\begin{aligned}
-\Delta y+\sigma y&=f
&&\text{in }\Omega,\\
\partial_{n}y&=u
&&\text{on }\Gamma,
\end{aligned}
```

with $\sigma>0$. Let

$$
V:=H^{1}(\Omega),
\qquad
U:=L^{2}(\Gamma).
$$

The weak equation is

$$
\int_{\Omega}\nabla y\cdot\nabla v\mathrm{d}x
+
\sigma\int_{\Omega}yv\mathrm{d}x
=
\int_{\Omega}fv\mathrm{d}x
+
\int_{\Gamma}u\gamma v\mathrm{d}s
$$

for every $v\in V$, where

$$
\gamma:H^{1}(\Omega)\to L^{2}(\Gamma)
$$

denotes the trace map, using the continuous trace embedding appropriate for a bounded Lipschitz domain.

The control operator is therefore characterized by

$$
\langle Bu,v\rangle_{V^{\ast},V}
=
\int_{\Gamma}u\gamma v\mathrm{d}s.
$$

The trace theorem is doing real work here: it makes the right-hand side a bounded functional of $v$ and therefore makes

$$
B:L^{2}(\Gamma)\to V^{\ast}
$$

a bounded linear map.

### 5.2 What does the adjoint control action become?

Suppose the objective contains

$$
\frac{\beta}{2}\lVert u\rVert_{L^{2}(\Gamma)}^{2}.
$$

For the residual convention

$$
E(y,u)=Ay-Bu-f,
$$

the reduced derivative is

$$
j'(u)
=
D_{u}J-D_{u}E^{\ast}p
=
\beta R_{U}u+B^{\ast}p.
$$

The adjoint control operator is characterized by

$$
\langle B^{\ast}p,h\rangle_{U^{\ast},U}
=
\langle Bh,p\rangle_{V^{\ast},V}
=
\int_{\Gamma}h\gamma p\mathrm{d}s.
$$

Thus, under the standard $L^{2}(\Gamma)$ Riesz identification,

$$
\boxed{
\nabla j(u)
=
\beta u+\gamma p.
}
$$

The control gradient is obtained from the **boundary trace of the adjoint**. This is the boundary-control analogue of the volume formula $\beta u+p$.

### 5.3 Dirichlet boundary control is more delicate

If the control prescribes

$$
y=u
\quad\text{on }\Gamma,
$$

then the control is part of the **essential boundary data**. This is qualitatively different from Neumann control. A Neumann datum enters the weak equation as a bounded boundary functional, whereas a Dirichlet datum determines the affine space in which the state itself lives.

There are two useful levels at which to see the issue.

### 5.3.1 Energy-level control: lift the boundary value

The trace operator

$$
\gamma:H^{1}(\Omega)\to H^{1/2}(\Gamma)
$$

suggests the natural energy-level control space

$$
U=H^{1/2}(\Gamma).
$$

Assume for simplicity that

```math
\begin{aligned}
-\Delta y&=f &&\text{in }\Omega,\\
y&=u &&\text{on }\Gamma.
\end{aligned}
```

Choose a continuous right inverse of the trace,

$$
\mathcal E:H^{1/2}(\Gamma)\to H^{1}(\Omega),
\qquad
\gamma(\mathcal E u)=u,
$$

and write

$$
y=z+\mathcal E u,
\qquad
z\in H_{0}^{1}(\Omega).
$$

The unknown $z$ now satisfies a standard homogeneous variational problem:

$$
\int_{\Omega}\nabla z\cdot\nabla v\mathrm{d}x
=
\int_{\Omega}fv\mathrm{d}x
-
\int_{\Omega}\nabla(\mathcal E u)\cdot\nabla v\mathrm{d}x
\qquad
\forall v\in H_{0}^{1}(\Omega).
$$

This is the same lifting idea used for prescribed nonhomogeneous Dirichlet data in [01 · Function spaces and weak PDEs](01-function-spaces-and-weak-pdes.md). The optimization-specific difference is that the lifting now changes whenever the control changes.

Suppose, for example,

$$
J(y,u)
=
\frac{1}{2}\lVert y-y_{d}\rVert_{L^{2}(\Omega)}^{2}
+
\frac{\beta}{2}\lVert u\rVert_{U}^{2}.
$$

The adjoint still solves the interior equation

```math
\begin{aligned}
-\Delta p&=y-y_{d} &&\text{in }\Omega,\\
p&=0 &&\text{on }\Gamma.
\end{aligned}
```

A boundary perturbation $h\in H^{1/2}(\Gamma)$ is represented in the state by a lifting $\mathcal E h$. After using the adjoint equation and Green's identity, the reduced derivative has the form

$$
j'(u)[h]
=
\beta (u,h)_{U}
-
\left\langle
\partial_{n}p,h
\right\rangle_{H^{-1/2}(\Gamma),H^{1/2}(\Gamma)}.
$$

The normal derivative here should first be understood as a boundary **covector**,

$$
\partial_{n}p\in H^{-1/2}(\Gamma),
$$

not automatically as an ordinary $L^{2}$ function. If $U=H^{1/2}(\Gamma)$, the gradient is obtained only after applying the Riesz inverse for that chosen boundary metric. This is a concrete instance of the derivative/gradient distinction from [02 · Duality, derivatives, and adjoints](02-duality-derivatives-and-adjoints.md).

The exact lifting $\mathcal E$ is a device for constructing the state. The final derivative is determined by the boundary value $h$, not by which particular interior extension was chosen.

### 5.3.2 Rougher $L^{2}$ boundary controls lead to a very-weak state

The choice $H^{1/2}(\Gamma)$ is mathematically natural but can be inconvenient numerically. One may want to optimize directly over

$$
U=L^{2}(\Gamma).
$$

Now a generic control $u\in L^{2}(\Gamma)$ need not be the trace of an $H^{1}(\Omega)$ function. The standard energy formulation therefore cannot simply impose $y=u$.

A useful alternative is the **transposition**, or **very-weak**, formulation. Assume enough elliptic regularity of the domain that

$$
Y:=H^{2}(\Omega)\cap H_{0}^{1}(\Omega)
$$

is an appropriate test space and that the normal derivative of $\psi\in Y$ belongs to $L^{2}(\Gamma)$. Starting from a smooth solution and integrating twice by parts gives

$$
-
\int_{\Omega}y\Delta\psi\mathrm{d}x
=
\int_{\Omega}f\psi\mathrm{d}x
-
\int_{\Gamma}u\partial_{n}\psi\mathrm{d}s.
$$

This identity still makes sense when

$$
y\in L^{2}(\Omega),
\qquad
u\in L^{2}(\Gamma),
$$

because all derivatives have been moved onto the smoother test function $\psi$. We therefore **define** the very-weak state by

$$
-
(y,\Delta\psi)_{L^{2}(\Omega)}
=
(f,\psi)_{L^{2}(\Omega)}
-
(u,\partial_{n}\psi)_{L^{2}(\Gamma)}
\qquad
\forall\psi\in Y.
$$

The idea is worth pausing over. Ordinary weak formulations reduce the derivative requirement on the state by moving one derivative to the test function. Transposition moves still more of the differential operator onto a deliberately smoother test space, allowing the state itself to have lower regularity.

For the $L^{2}$ tracking functional

$$
J(y,u)
=
\frac{1}{2}\lVert y-y_{d}\rVert_{L^{2}(\Omega)}^{2}
+
\frac{\beta}{2}\lVert u\rVert_{L^{2}(\Gamma)}^{2},
$$

the adjoint has the stronger regularity role. Formally, and under the same regularity assumptions,

```math
\begin{aligned}
-\Delta p&=y-y_{d} &&\text{in }\Omega,\\
p&=0 &&\text{on }\Gamma,
\end{aligned}
```

with $p\in Y$. The reduced derivative becomes

$$
\boxed{
j'(u)[h]
=
(\beta u-\partial_{n}p,h)_{L^{2}(\Gamma)}.
}
$$

Thus the rough state formulation trades regularity between the state and the multiplier: the state may live only in $L^{2}(\Omega)$, while the adjoint must be regular enough for its normal derivative to define the boundary gradient.

This is not a free simplification. The transposition theory uses stronger assumptions on the domain than the basic $H^{1}$ variational formulation, and its finite-element discretization needs special care. The full well-posedness and approximation theory is therefore best read in the dedicated treatment in Manzoni–Quarteroni–Salsa §5.9 and §5.11.

### 5.3.3 The important lesson is where the control enters the formulation

For source control, the control appears in the volume residual. For Neumann control, it appears as a natural boundary functional. For Dirichlet control, it changes the state trace itself.

That distinction determines:

- the natural control space;
- whether the state lives in a linear or affine energy space;
- whether a lifting or a transposition formulation is useful;
- which boundary quantity of the adjoint represents the reduced derivative.

So the phrase "boundary control" hides several mathematically different cases. The control operator must be defined together with the state space and pairing; it is not merely a different location at which the same vector is inserted.

## 6. Coefficient identification makes the state map nonlinear in a different way

A second important generalization is to optimize a coefficient inside the PDE operator. Consider

```math
\begin{aligned}
\min_{y,a}\quad
J(y,a)
&:=
\frac{1}{2}\lVert y-y_{d}\rVert_{L^{2}(\Omega)}^{2}
+
\frac{\alpha}{2}\lVert a-a_{\mathrm{ref}}\rVert_{L^{2}(\Omega)}^{2},\\
\text{subject to}\quad
-\nabla\cdot(a\nabla y)&=f
\quad\text{in }\Omega,\\
y&=0
\quad\text{on }\partial\Omega.
\end{aligned}
```

Assume the admissible coefficients belong to $L^{\infty}(\Omega)$ and satisfy uniform positive bounds

$$
0<a_{\min}
\leq
a(x)
\leq
a_{\max}
$$

almost everywhere. These bounds are not merely optimization constraints: the lower bound preserves ellipticity of the state equation. At the energy level, coefficient perturbations $h$ are naturally taken in $L^{\infty}(\Omega)$ so that products such as $h\nabla y$ define bounded weak-form terms. Writing an $L^{2}$ gradient density later is therefore a representation step that may require additional regularity or a chosen coefficient-space setting.

For fixed $a$, the PDE is linear in $y$. As a map

$$
a\mapsto y(a),
$$

however, the problem is nonlinear because the coefficient changes the state operator itself.

### 6.1 Derive the coefficient sensitivity

The weak residual is

$$
\langle E(y,a),v\rangle
=
\int_{\Omega}a\nabla y\cdot\nabla v\mathrm{d}x
-
\int_{\Omega}fv\mathrm{d}x.
$$

Perturb the coefficient by $h$ and let

$$
z=S'(a)h.
$$

Differentiating the state equation gives

$$
\int_{\Omega}a\nabla z\cdot\nabla v\mathrm{d}x
+
\int_{\Omega}h\nabla y\cdot\nabla v\mathrm{d}x
=
0.
$$

Therefore

$$
\boxed{
\int_{\Omega}a\nabla z\cdot\nabla v\mathrm{d}x
=
-
\int_{\Omega}h\nabla y\cdot\nabla v\mathrm{d}x.
}
$$

The right-hand side now depends on the current state gradient. This is the characteristic sensitivity structure of coefficient identification.

### 6.2 Use the adjoint to remove the sensitivity

With the project sign convention,

$$
\mathcal L(y,a,p)
=
J(y,a)
-
\left(
\int_{\Omega}a\nabla y\cdot\nabla p\mathrm{d}x
-
\int_{\Omega}fp\mathrm{d}x
\right).
$$

Stationarity with respect to the state gives

$$
\int_{\Omega}a\nabla p\cdot\nabla v\mathrm{d}x
=
\int_{\Omega}(y-y_{d})v\mathrm{d}x
\qquad
\forall v\in H_{0}^{1}(\Omega).
$$

Now differentiate the Lagrangian with respect to the coefficient in direction $h$:

$$
D_{a}\mathcal L(y,a,p)[h]
=
\alpha(a-a_{\mathrm{ref}},h)_{L^{2}(\Omega)}
-
\int_{\Omega}h\nabla y\cdot\nabla p\mathrm{d}x.
$$

At a feasible state and its adjoint this is the reduced derivative, so formally the $L^{2}$ gradient density is

$$
\boxed{
\nabla j(a)
=
\alpha(a-a_{\mathrm{ref}})
-
\nabla y\cdot\nabla p.
}
$$

The sign of the last term depends on the Lagrangian/residual convention. The formula above follows the global convention used throughout these background notes. Sources using $J+\langle p,E\rangle$ obtain the opposite adjoint sign and an equivalent final optimality system after the multiplier is changed consistently.

### 6.3 Why this example matters beyond inverse problems

The coefficient example separates two ideas that are easy to conflate:

- a PDE can be linear in its **state** for fixed data;
- the control-to-state map can still be nonlinear if the **control changes the operator**.

This is why coefficient identification leads naturally to nonlinear reduced optimization and nonlinear all-at-once KKT systems even when every individual forward solve is a linear elliptic solve.

## 7. Observation operators make the inverse-problem structure explicit

Tracking the full state in $L^{2}(\Omega)$ is only one measurement model. Let

$$
C:Y\to O
$$

be an observation operator into a Hilbert observation space $O$, and consider

$$
J(y,u)
=
\frac{1}{2}\lVert Cy-z_{d}\rVert_{O}^{2}
+
\frac{\alpha}{2}\lVert u-u_{\mathrm{ref}}\rVert_{U}^{2}.
$$

Then

$$
D_{y}J(y,u)[z]
=
(Cy-z_{d},Cz)_{O}
=
\left\langle
C^{\ast}R_{O}(Cy-z_{d}),z
\right\rangle.
$$

The observation residual is therefore pushed back into the adjoint equation by

$$
C^{\ast}R_{O}.
$$

This is more than compact notation. Different observations require different state regularity, and the adjoint inherits the corresponding dual object.

### 7.1 Subdomain observations are close to ordinary $L^{2}$ tracking

Let $\omega\subset\Omega$ and define

$$
Cy=y\rvert_{\omega},
\qquad
O=L^{2}(\omega).
$$

For $y\in L^{2}(\Omega)$ this is a bounded restriction operator. The tracking derivative is

$$
D_{y}J(y,u)[z]
=
\int_{\omega}(y-z_{d})z\mathrm{d}x.
$$

Equivalently, the adjoint receives the observation mismatch extended by zero outside $\omega$. Analytically, this behaves much like full-domain tracking; only the support of the adjoint source changes.

### 7.2 Boundary traces are bounded observations at the energy level

If $Y=H^{1}(\Omega)$, the trace theorem gives a bounded map

$$
\gamma:Y\to H^{1/2}(\Gamma)
$$

and, on a bounded Lipschitz boundary, a continuous embedding into $L^{2}(\Gamma)$. A boundary-state observation on $\Gamma_{\mathrm{obs}}\subset\Gamma$ can therefore be written schematically as

$$
Cy
=
\gamma y\rvert_{\Gamma_{\mathrm{obs}}}.
$$

The adjoint forcing is no longer a volume source: $C^{\ast}$ turns the boundary mismatch into a functional acting on state test functions through their traces. In a weak PDE formulation this typically appears as a natural boundary contribution.

This case is still comparatively mild because the trace is already a bounded operator on the energy space.

### 7.3 Point sensors require a smoother state space

A point sensor at $\xi\in\Omega$ asks for

$$
y\mapsto y(\xi).
$$

For a generic $H^{1}(\Omega)$ state this is not a bounded functional in dimensions two and three. Point evaluation becomes well defined only after stronger regularity is available. For a smooth or convex domain, an elliptic problem with sufficiently regular data may give

$$
Y=H^{2}(\Omega)\cap H_{0}^{1}(\Omega).
$$

When $d\leq3$, Sobolev embedding gives

$$
Y\hookrightarrow C(\overline\Omega),
$$

so evaluation is continuous. For sensor locations $\xi_{1},\ldots,\xi_{m}$ define

$$
C:Y\to\mathbb R^{m},
\qquad
Cy=
\begin{bmatrix}
y(\xi_{1})\\
\vdots\\
y(\xi_{m})
\end{bmatrix}.
$$

If the sensor residual is $r\in\mathbb R^{m}$, then

$$
C^{\ast}r
=
\sum_{j=1}^{m}r_{j}\delta_{\xi_{j}},
$$

where $\delta_{\xi_{j}}$ is the Dirac point-evaluation functional. The adjoint equation therefore contains singular sources. In a Poisson-type model it has the formal form

$$
-\Delta p
=
\sum_{j=1}^{m}
\bigl(y(\xi_{j})-z_{d,j}\bigr)\delta_{\xi_{j}}.
$$

The state had to become smoother so that the observation was meaningful, yet the **adjoint becomes rougher** because the transpose observation operator contains Dirac masses. A standard $H_{0}^{1}$ adjoint formulation may no longer be appropriate; a very-weak/transposition interpretation is natural.

This reversal is a recurring principle:

```text
more singular observation
        ↓
stronger regularity needed for the state
        ↓
more singular C* applied to the data mismatch
        ↓
lower regularity for the adjoint
```

### 7.4 Boundary-flux observations require still more care

Suppose the state solves

```math
\begin{aligned}
-\Delta y&=f+u &&\text{in }\Omega,\\
y&=0 &&\text{on }\Gamma,
\end{aligned}
```

and the measured quantity is the outgoing flux on $\Gamma_{0}\subset\Gamma$:

$$
Cy
=
\partial_{n}y\rvert_{\Gamma_{0}}.
$$

The minimal energy regularity $y\in H_{0}^{1}(\Omega)$ is not enough to regard $\partial_{n}y$ as an $L^{2}$ boundary function. If elliptic regularity gives

$$
y\in Y:=H^{2}(\Omega)\cap H_{0}^{1}(\Omega),
$$

then the normal derivative has a trace and the cost

$$
\frac{1}{2}
\left\lVert
\partial_{n}y-z_{d}
\right\rVert_{L^{2}(\Gamma_{0})}^{2}
$$

is meaningful.

Its state derivative is

$$
D_{y}J(y,u)[z]
=
\int_{\Gamma_{0}}
\left(
\partial_{n}y-z_{d}
\right)
\partial_{n}z
\mathrm{d}s.
$$

This expression already reveals the difficulty: the right-hand side is not a bounded functional of a generic $z\in H_{0}^{1}(\Omega)$ because $\partial_{n}z$ is not defined there.

Using the strong residual and integrating by parts twice leads formally to an adjoint that is harmonic in the interior and whose **Dirichlet boundary data are the flux mismatch**. With the sign convention used in these notes, one convenient form is

```math
\begin{aligned}
-\Delta p&=0 &&\text{in }\Omega,\\
p&=-\chi_{\Gamma_{0}}
\left(
\partial_{n}y-z_{d}
\right)
&&\text{on }\Gamma.
\end{aligned}
```

The datum on the boundary is only $L^{2}$ in this model, so it need not be the trace of an $H^{1}$ function. The adjoint is therefore naturally interpreted as a **very-weak solution**, possibly only in $L^{2}(\Omega)$. For a distributed source control, the reduced derivative then has the familiar volume form

$$
j'(u)[h]
=
(\beta u+p,h)_{L^{2}(\Omega)}.
$$

Boundary-flux observation is a useful warning against treating $C$ as an innocent matrix. The choice of observation can change the correct state space, the multiplier space, and even which formulation of the adjoint PDE is mathematically meaningful.

### 7.5 Observation, regularization, and metric remain different objects

The decomposition

```text
unknown parameter/control
        ↓
PDE state map
        ↓
observation operator C
        ↓
predicted data
        ↓
misfit against measurements
```

is especially useful in inverse problems. The observation operator says **what is measured**. The regularization term says **which unknowns are preferred or penalized** and often stabilizes an ill-posed inverse problem. The optimization metric says **how a derivative covector is converted into a search direction**.

Even if a mass or stiffness operator appears numerically in more than one of these roles, the roles should not be merged conceptually.

## 8. Where these extensions stop being small variations

The topics above still reuse the same core architecture:

```text
state equation
    ↓
linearized state equation
    ↓
adjoint equation
    ↓
reduced derivative
    ↓
optional second-order linearizations
```

Several nearby subjects introduce genuinely new analytical machinery and deserve more than a few extra formulas:

- **time-dependent control** introduces vector-valued function spaces, weak time derivatives, initial/terminal conditions, and forward/backward evolution;
- **nonsmooth objectives** replace ordinary derivatives by subdifferentials or generalized derivatives;
- **state constraints** can produce multipliers with lower regularity, including measure-valued multipliers;
- **variational-inequality state equations** make the control-to-state map itself nonsmooth;
- **shape optimization** changes the domain rather than a field defined on a fixed domain.

These topics are collected as an oriented continuation in [07b · Extensions of PDE-constrained optimization](07b-extensions-of-pde-constrained-optimization.md). That companion derives representative structures where they can be explained compactly and explicitly defers the theory that cannot be made reliable in a short background treatment.

## Where to go next

- **For the current project:** return to the `nmopt` concept manual. In particular:
  - [manual 05 · Reduced state–adjoint formulation](../concepts/05-reduced-state-adjoint-formulation.md) uses the same state/adjoint composition at the executable formulation boundary;
  - [manual 07 · Optimality systems and KKT](../concepts/07-optimality-systems-and-kkt.md) is the natural continuation for nonlinear/all-at-once first-order systems and Newton-like linearizations.
- **For a broader map of PDE optimal control:** [07b · Extensions of PDE-constrained optimization](07b-extensions-of-pde-constrained-optimization.md) introduces evolution problems, nonsmooth control, state constraints, Pontryagin/SQH ideas, more complicated PDEs, OCP error analysis, and shape optimization.
- **For implementation prerequisites:** continue with [08 · Modern C++ for scientific code](08-modern-cpp-for-scientific-code.md) and [09 · The deal.II finite-element workflow](09-dealii-finite-element-workflow.md) if you are moving from the mathematics into the code base.

## References and further reading

Publication details, stable links, and access notes are collected in the [Background reference catalogue](references.md).

- **Andrea Manzoni, Alfio Quarteroni, and Sandro Salsa, *Optimal Control of Partial Differential Equations*.**
  - Chapter 5 for elliptic controls and observations.
  - §§5.8–5.11 for boundary-flux observations, transposition/very-weak solutions, point observations, and Dirichlet boundary control.
  - Chapter 9 for Banach/Hilbert-space nonlinear OCPs, semilinear state equations, coefficient estimation, second-order conditions, iterative methods, and SQP.
- **Luca Heltai, NMOPT course.**
  - Lecture 8, **General Optimality Systems for PDE-Constrained Optimization**, for the implicit-function, nonlinear sensitivity/adjoint, KKT, and second-order viewpoint.
  - Lecture 14, **Boundary Control and Parameter Estimation**, for trace-based boundary controls, liftings, and coefficient/parameter sensitivities.
  - Lecture 15, **One-Shot KKT and PDAS for Inverse Poisson Coefficient Identification**, for nonlinear coefficient identification with Newton/active-set ideas.
- **General nonlinear-programming references in the catalogue.**
  - Use these for full Newton–KKT and SQP theory; this companion only develops the PDE-specific structure needed to recognize those methods.
