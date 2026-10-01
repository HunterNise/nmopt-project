# 07b · Extensions of PDE-constrained optimization

**Background navigation:** [Index](README.md) \
Main chapter: [07 · PDE-constrained optimization](07-pde-constrained-optimization.md) \
Previous: [07a · Further notes on PDE-constrained optimization](07a-further-pde-constrained-optimization-notes.md) \
Next: [08 · Modern C++ for scientific code](08-modern-cpp-for-scientific-code.md)

## Purpose

The prerequisite path ends with [07 · PDE-constrained optimization](07-pde-constrained-optimization.md). The companion [07a · Further notes on PDE-constrained optimization](07a-further-pde-constrained-optimization-notes.md) stays close to that chapter by extending its machinery to nonlinear elliptic equations, second-order methods, boundary controls, and coefficient identification.

This second companion has a different purpose. It collects major directions that are treated in broader PDE optimal-control courses but require enough new mathematics that they should not be compressed into the main prerequisite chapter. The goal is to make these directions **recognizable and conceptually connected**, not to replace their full treatment in the cited sources.

Each section therefore starts from a familiar state–adjoint idea, identifies the genuinely new ingredient, and derives one representative formula or optimality structure. When the missing theory is substantial – for example Bochner-space evolution equations, measure-valued multipliers, generalized derivatives, or shape calculus – the notes state clearly where the short treatment stops.

None of the material in this file is required to understand the current project. It is useful when asking how the same numerical architecture might extend to new PDEs, controls, constraints, or objective models.

## Before you start

You should be comfortable with:

- the main [07 · PDE-constrained optimization](07-pde-constrained-optimization.md) state–sensitivity–adjoint derivation;
- the nearby nonlinear and inverse-problem extensions in [07a · Further notes on PDE-constrained optimization](07a-further-pde-constrained-optimization-notes.md);
- variational inequalities, normal cones, and active sets from [06 · Constrained optimization](06-constrained-optimization.md);
- weak PDE formulations from [01 · Function spaces and weak PDEs](01-function-spaces-and-weak-pdes.md).

For the evolution section, elementary ODE initial-value problems are useful. No prior Bochner-space theory is assumed; the essential vocabulary is introduced here, while proofs are deferred.

## What you will be able to do

After these notes, you should be able to:

- distinguish the main structural features of parabolic and hyperbolic control problems;
- interpret Bochner spaces, the Gelfand triple, and the standard parabolic energy space at a working level;
- explain why parabolic control produces a forward state equation and a backward adjoint equation;
- explain why a discrete time adjoint is the transpose of the discrete state evolution rather than merely a time-reversed continuous equation;
- define the convex subdifferential, explain why it replaces an ordinary derivative at nonsmooth points, and formulate the first-order condition for an $L^{1}$-regularized sparse control;
- explain why state constraints are harder than control bounds and why their multipliers may leave $L^{2}$;
- recognize obstacle and variational-inequality state equations as a source of nonsmooth control-to-state maps;
- derive the Hamiltonian state–costate equations for a basic dynamical problem and explain what the Pontryagin minimum principle adds beyond smooth stationarity;
- connect Newton–KKT, SQP, and sequential quadratic Hamiltonian viewpoints without treating them as project prerequisites;
- identify what changes when the state PDE is Stokes, Navier–Stokes, wave-like, or coupled;
- distinguish OCP-specific error analysis and shape optimization from the core [07 · PDE-constrained optimization](07-pde-constrained-optimization.md) machinery.

## Roadmap

We begin with evolution problems because time dependence changes both the function spaces and the causal structure of the state/adjoint pair. We first compare representative parabolic and hyperbolic equations, then introduce Bochner spaces as the natural way to describe time-dependent fields. We next turn to nonsmooth optimization, starting from the convex subdifferential itself before applying it to sparse $L^{1}$ control, state constraints, and variational-inequality state equations.

The Hamiltonian part then derives the state–costate equations from the Lagrangian, explains the pointwise control condition in the Pontryagin principle, and compares classical SQP with the SQH viewpoint used later in the course. The final sections give concrete PDE examples – Stokes, Navier–Stokes, waves, and coupled reaction–diffusion/ODE systems – and identify exactly which new analytical or numerical difficulty each one adds.

## 1. Time dependence changes both the function spaces and the optimization geometry

A stationary PDE-constrained problem asks for fields on a spatial domain $\Omega$. A time-dependent problem asks for **trajectories of fields** on a space-time cylinder

$$
Q_{T}\coloneqq\Omega\times(0,T).
$$

The control may also vary in time, so instead of one $u(x)$ we optimize a history $u(x,t)$. Objectives may measure tracking over the whole interval,

$$
\frac{1}{2}
\int_{0}^{T}
\lVert y(t)-y_{d}(t)\rVert^{2}
\mathrm{d}t,
$$

at the final time,

$$
\frac{1}{2}
\lVert y(T)-y_{T}\rVert^{2},
$$

or both.

Two broad evolution classes already show why "add time" is not a single modification.

### 1.1 Parabolic examples: diffusion with memory of an initial state

A linear heat or diffusion-reaction equation has the form

$$
\begin{aligned}
y_{t}
-\nabla\cdot(\kappa\nabla y)
+c y
&=Bu+f
&&\text{in }Q_{T},\\
y&=0
&&\text{on }\partial\Omega\times(0,T),\\
y(0)&=y_{0}
&&\text{in }\Omega.
\end{aligned}
$$

An advection-diffusion model adds transport,

$$
y_{t}-\nu\Delta y+b\cdot\nabla y=Bu+f.
$$

These are **parabolic** problems. Diffusion damps high-frequency components and has a regularizing effect as time advances. The state is determined by an initial value and evolves forward.

The cost may observe the whole trajectory or only its endpoint. For example,

$$
J(y,u)
= \frac{1}{2}
\int_{0}^{T}
\lVert y(t)-y_{d}(t)\rVert_{L^{2}(\Omega)}^{2}
\mathrm{d}t
+\frac{\alpha}{2}
\int_{0}^{T}
\lVert u(t)\rVert_{U}^{2}
\mathrm{d}t
+\frac{\beta}{2}
\lVert y(T)-y_{T}\rVert_{L^{2}(\Omega)}^{2}.
$$

The terminal term is important because it becomes a **final condition for the adjoint**.

### 1.2 Hyperbolic example: the wave equation carries two initial data

A model wave-control problem is

$$
\begin{aligned}
y_{tt}-c^{2}\Delta y&=Bu+f
&&\text{in }Q_{T},\\
y&=0
&&\text{on }\partial\Omega\times(0,T),\\
y(0)&=y_{0},\\
y_{t}(0)&=v_{0}.
\end{aligned}
$$

The wave equation is **hyperbolic**, not parabolic. It propagates information at finite speed and does not have the same dissipative smoothing mechanism as the heat equation. Two initial conditions are required because the equation is second order in time.

This changes both analysis and numerics. The natural energy contains the displacement and the velocity, time integration must respect wave stability, and terminal conditions for the adjoint involve two quantities as well. Exact controllability of waves also has a geometric/time-of-propagation character that is very different from diffusion.

The shared optimal-control pattern survives, but the state solver hidden inside it is fundamentally different.

### 1.3 Bochner spaces are ordinary $L^{p}$ spaces with values in a function space

Writing

$$
y(t)\in H_{0}^{1}(\Omega)
$$

for every time is not enough. We need a space that measures the whole trajectory.

Let $X$ be a Banach space. The **Bochner space**

$$
L^{p}(0,T;X)
$$

consists, roughly, of strongly measurable functions

$$
t\mapsto v(t)\in X
$$

whose $X$-norm is $p$-integrable:

$$
\int_{0}^{T}
\lVert v(t)\rVert_{X}^{p}
\mathrm{d}t
\lt\infty.
$$

For $p=2$,

$$
\lVert v\rVert_{L^{2}(0,T;X)}^{2}
= \int_{0}^{T}\lVert v(t)\rVert_{X}^{2}\thinspace\mathrm{d}t.
$$

If $X=\mathbb R^{n}$ this is just the familiar space of square-integrable vector-valued time signals. In PDEs, $X$ itself is a spatial function space, so one should think of a Bochner-space element as a time-indexed family of fields with an integrated spatial norm.

For parabolic equations one commonly uses a **Gelfand triple**

$$
V\hookrightarrow H\hookrightarrow V^{\ast},
$$

with continuous dense embeddings and $H$ identified with its dual. For the heat equation,

$$
V=H_{0}^{1}(\Omega),
\qquad
H=L^{2}(\Omega),
\qquad
V^{\ast}=H^{-1}(\Omega).
$$

The natural energy space is

$$
W(0,T)
\coloneqq
\left\lbrace
y\in L^{2}(0,T;V):
y_{t}\in L^{2}(0,T;V^{\ast})
\right\rbrace.
$$

This choice is not decorative notation. The spatial elliptic operator naturally maps $V$ into $V^{\ast}$, so the equation

$$
y_{t}+Ay=Bu+f
$$

can hold in $V^{\ast}$ for almost every time even when neither $y_{t}$ nor $Ay$ is an $L^{2}(\Omega)$ function.

A fundamental parabolic result gives the continuous embedding

$$
W(0,T)\hookrightarrow C([0,T];H).
$$

That is what makes the endpoint values $y(0)$ and $y(T)$ meaningful. The proof belongs to the full evolution-PDE theory; for the present orientation, the important point is how the spaces fit the equation.

### 1.4 The weak parabolic state equation

Let $A:V\to V^{\ast}$ be induced by a continuous coercive bilinear form $a$. The weak equation is

$$
\langle y_{t}(t),v\rangle_{V^{\ast},V}
+a(y(t),v)
= \langle Bu(t)+f(t),v\rangle_{V^{\ast},V}
$$

for every $v\in V$ and almost every $t$.

The initial value

$$
y(0)=y_{0}
$$

selects the trajectory. Numerically, a time-stepping method therefore advances the state from early times to later times.

This causal direction is the first structural difference from a stationary elliptic solve.

## 2. Integration by parts in time produces the backward adjoint

The backward direction of the adjoint is not an arbitrary convention. It follows from moving the time derivative off the state variation.

Use the residual

$$
E(y,u)
\coloneqq
y_{t}+Ay-Bu-f
$$

and the same Lagrangian convention as the rest of these notes:

$$
\mathcal L(y,u,p)
= J(y,u)
-\int_{0}^{T}\langle p,E(y,u)\rangle\thinspace\mathrm{d}t.
$$

Let $z$ be a state variation satisfying

$$
z(0)=0,
$$

because the initial condition is fixed. The time-derivative contribution to the variation is

$$
-\int_{0}^{T}\langle p,z_{t}\rangle\thinspace\mathrm{d}t.
$$

The weak integration-by-parts identity gives

$$
-\int_{0}^{T}\langle p,z_{t}\rangle\thinspace\mathrm{d}t
= -(p(T),z(T))_{H}
+\int_{0}^{T}\langle p_{t},z\rangle\thinspace\mathrm{d}t,
$$

because $z(0)=0$.

The terminal cost contributes

$$
\beta(y(T)-y_{T},z(T))_{H}.
$$

For the terminal variation to vanish for every $z(T)$, the adjoint must satisfy

$$
\boxed{
p(T)=\beta(y(T)-y_{T}).
}
$$

For the interior variation to vanish, we obtain

$$
\boxed{
-p_{t}+A^{\ast}p
= R_{H}(y-y_{d}).
}
$$

Thus the adjoint is a **final-value problem**. It is solved backward from $T$ toward $0$.

The control derivative is

$$
j'(u)
= \alpha R_{U}u+B^{\ast}p.
$$

The complete computational pattern becomes

```text
control trajectory u(t)
        ↓
state equation forward from t=0
        ↓
state history y(t)
        ↓
adjoint equation backward from t=T
        ↓
control derivative over the whole time interval
```

The forward/backward structure is the evolution analogue of the two elliptic solves in [07 · PDE-constrained optimization](07-pde-constrained-optimization.md).

### 2.1 Why the state history matters

The adjoint right-hand side depends on $y(t)$. When integrating the adjoint backward, the state at the corresponding times must therefore be available.

A simple implementation may store every state time level. Large simulations may instead use checkpointing or recomputation strategies. This is an algorithmic consequence of the mathematics: the reverse adjoint sweep needs information generated during the forward state sweep.

## 3. A discrete adjoint belongs to the discrete state evolution

Time discretization makes the DTO/OTD distinction particularly concrete.

Suppose implicit Euler with time step $\tau$ gives the discrete state equation

$$
(M+\tau K)Y^{k}
= MY^{k-1}
+\tau B U^{k}
+\tau F^{k},
\qquad
k=1,\ldots,N.
$$

Here $M$ is a spatial mass matrix and $K$ the discrete elliptic operator.

If we first define this discrete state evolution and then differentiate the resulting discrete optimization problem, the adjoint recurrence has the transpose structure

$$
(M+\tau K)^{\mathsf T}P^{k}
= M^{\mathsf T}P^{k+1}
+\tau M(Y^{k}-Y_{d}^{k}),
$$

with the appropriate terminal contribution at $k=N$.

The important lesson is not the exact indexing convention. It is the structural statement:

> the discrete adjoint is the transpose of the **discrete state evolution** with respect to the declared discrete pairings.

Merely deriving a continuous adjoint and then informally reversing a time loop does not guarantee this property. Time quadrature, mass matrices, the placement of controls within time intervals, and the terminal objective all affect the transpose relation.

A full treatment of adjoint-consistent time discretization requires more space than this orientation. The NMOPT evolution lectures and Manzoni–Quarteroni–Salsa Chapters 7–8 develop this subject systematically.

## 4. Nonsmooth objectives need a replacement for the derivative

The smooth chapters repeatedly used the first-order condition

$$
j'(u)=0.
$$

That statement assumes that a single linear functional describes the first-order change of $j$ at $u$. A convex function can fail to have such a derivative and still possess a precise first-order geometry.

The basic object is the **convex subdifferential**.

### 4.1 A subgradient is a supporting linear functional

Let $X$ be a Banach space and let

$$
\phi:X\to\mathbb R\cup\lbrace+\infty\rbrace
$$

be a proper convex functional. A covector

$$
\xi\in X^{\ast}
$$

is a **subgradient** of $\phi$ at $x$ if

$$
\boxed{
\phi(v)
\geq
\phi(x)
+\langle\xi,v-x\rangle
\qquad
\forall v\in X.
}
$$

The set of all such subgradients is the **subdifferential**

$$
\partial\phi(x)
\coloneqq
\left\lbrace
\xi\in X^{\ast}:
\phi(v)\geq\phi(x)+\langle\xi,v-x\rangle
\ \forall v
\right\rbrace.
$$

Geometrically, every $\xi\in\partial\phi(x)$ defines an affine supporting hyperplane below the graph of the convex function.

If $\phi$ is Fréchet differentiable at $x$, convexity implies

$$
\partial\phi(x)=\lbrace\phi'(x)\rbrace.
$$

So the subdifferential does not replace derivatives when derivatives work; it extends the same first-order idea to corners and flat set-valued slopes.

The key optimality fact is equally simple:

$$
\boxed{
0\in\partial\phi(\bar x)
\quad\Longleftrightarrow\quad
\bar x\text{ minimizes the convex functional }\phi.
}
$$

This is the nonsmooth analogue of ``gradient equals zero.''

### 4.2 The absolute value shows why the subdifferential is set-valued

For

$$
\phi(s)=\lvert s\rvert,
$$

one finds

$$
\partial \lvert s\rvert
= \begin{cases}
\lbrace1\rbrace, & s\gt0,\\
[-1,1], & s=0,\\
\lbrace-1\rbrace, & s\lt0.
\end{cases}
$$

Away from zero this agrees with the ordinary derivative. At zero there is no single tangent slope, but every slope between $-1$ and $1$ defines a supporting line.

This interval is exactly what allows an optimum to remain at zero even when the smooth part of an objective has a nonzero derivative of moderate size.

### 4.3 $L^{1}$ regularization produces pointwise subgradient conditions

Consider

$$
J(y,u)
= \frac{1}{2}\lVert y-y_{d}\rVert_{L^{2}(\Omega)}^{2}
+\frac{\alpha}{2}\lVert u\rVert_{L^{2}(\Omega)}^{2}
+\beta\lVert u\rVert_{L^{1}(\Omega)}.
$$

Let $j_{\mathrm{sm}}$ denote the differentiable reduced tracking plus quadratic part. The first-order condition becomes

$$
\boxed{
0
\in
j_{\mathrm{sm}}'(u)
+\beta\partial\lVert u\rVert_{L^{1}}.
}
$$

In an $L^{2}$ representation on a bounded domain, this means that there exists a function $\lambda$ such that

$$
j_{\mathrm{sm}}'(u)+\beta\lambda=0,
$$

with the pointwise conditions

$$
\lambda(x)
\in
\begin{cases}
\lbrace1\rbrace, & u(x)\gt0,\\
[-1,1], & u(x)=0,\\
\lbrace-1\rbrace, & u(x)\lt0.
\end{cases}
$$

Suppose the smooth reduced gradient is $g$. Then at a point where

$$
\lvert g(x)\rvert\lt\beta,
$$

the inclusion can be satisfied with

$$
u(x)=0.
$$

This is the mechanism behind **sparsity promotion**: a whole interval of smooth-gradient values is compatible with an exactly zero control.

The proximal map of the scalar absolute value makes the same mechanism explicit. For $\tau\gt0$,

$$
\mathrm{prox}_{\tau\beta\lvert\cdot\rvert}(s)
= \mathrm{sign}(s)
\max\lbrace\lvert s\rvert-\tau\beta,0\rbrace,
$$

the familiar soft-thresholding rule.

### 4.4 Constraints also fit the subdifferential language

Let $K\subset X$ be a closed convex set and define its indicator functional

$$
I_{K}(x)
= \begin{cases}
0, & x\in K,\\
+\infty, & x\notin K.
\end{cases}
$$

Then

$$
\partial I_{K}(x)=N_{K}(x),
$$

where $N_{K}(x)$ is the normal cone from [06 · Constrained optimization](06-constrained-optimization.md). Thus the constrained problem

$$
\min_{x\in K} f(x)
$$

can be written as the unconstrained nonsmooth problem

$$
\min_{x} f(x)+I_{K}(x),
$$

with first-order condition

$$
0
\in
f'(x)+N_{K}(x).
$$

This unifies two topics that can otherwise look unrelated: box constraints and $L^{1}$ sparsity are both naturally expressed through set-valued first-order operators.

### 4.5 Why ordinary Newton is no longer enough

A classical Newton method differentiates a smooth equation. The inclusion

$$
0\in g(u)+\beta\partial\lVert u\rVert_{1}
$$

is set-valued and changes branch when components cross zero.

Useful numerical approaches include proximal methods, active-set formulations, complementarity reformulations, and semismooth Newton methods. The exact choice depends on how the nonsmooth term and the PDE coupling are represented.

The NMOPT nonsmooth lecture develops these connections further. A complete treatment requires generalized derivatives and convergence theory, which are beyond this compendium.

## 5. State constraints are not just another set of box bounds

A control constraint acts directly on the optimization variable:

$$
u_{a}
\leq
u
\leq
u_{b}.
$$

A pointwise state constraint instead asks, for example,

$$
y(u)(x)
\leq
y_{b}(x)
$$

throughout some region.

After reduction this becomes

$$
S(u)
\leq
y_{b}.
$$

The inequality is therefore filtered through a PDE solution operator. A local change in the control can affect the state elsewhere, and the multiplier acts on the **state constraint space**, not directly on the control coefficients.

### 5.1 Why multiplier regularity becomes difficult

For simple finite-dimensional inequalities, the multiplier is a vector. For pointwise state constraints in PDE problems, the natural dual of the space in which pointwise inequalities are imposed may contain objects rougher than functions.

In classical elliptic state-constrained control, the multiplier can be a **measure** rather than an $L^{2}$ function. Then expressions that looked like ordinary multiplier-weighted residuals in [06 · Constrained optimization](06-constrained-optimization.md) must be interpreted through measure/function pairings.

This has several consequences:

- the adjoint equation may receive a measure-valued source;
- the adjoint can lose regularity;
- straightforward finite-element approximation of the multiplier becomes delicate;
- constraint qualifications are more subtle than for control bounds.

The Manzoni–Quarteroni–Salsa treatment uses this difficulty to motivate regularized mixed control–state constraints, including Laurentiev-type regularization. That route can recover optimization problems more amenable to active-set methods, but it belongs to a full state-constraint treatment rather than a prerequisite appendix.

## 6. The PDE itself can be nonsmooth: obstacle and variational-inequality states

So far the state equation has been an equality

$$
E(y,u)=0.
$$

In contact, obstacle, or complementarity models, the state may instead satisfy a variational inequality. A schematic obstacle problem is

$$
\begin{aligned}
y&\geq\psi,\\
Ay-f-u&\geq0,\\
(y-\psi)(Ay-f-u)&=0.
\end{aligned}
$$

At each point either the obstacle is inactive and the PDE equality holds, or the state touches the obstacle and a reaction force becomes active.

The active region can change when the control changes. Consequently the solution map

$$
u
\mapsto y(u)
$$

is generally not Fréchet differentiable in the ordinary sense across switching events.

This breaks the smooth sensitivity derivation from [07 · PDE-constrained optimization](07-pde-constrained-optimization.md) at its first step. One then needs weaker differentiability concepts, complementarity formulations, or semismooth/generalized Newton methods.

The conceptual connection with PDAS is strong: both problems identify active and inactive regions and solve different linear equations on the current branches. The analytical theory is nevertheless richer because the nonsmoothness now lies inside the state system itself.

## 7. The Hamiltonian exposes the local structure of dynamical optimal control

The Lagrangian formulation used throughout these notes already contains the state, adjoint, and control equations. For an evolution problem it is often useful to reorganize the same information around a function called the **Hamiltonian**.

Start with a finite-dimensional dynamical system

$$
\dot y(t)=F(y(t),u(t)),
\qquad
y(0)=y_{0},
$$

and a cost

$$
J(y,u)
= \Phi(y(T))
+\int_{0}^{T}\ell(y(t),u(t))\thinspace\mathrm{d}t.
$$

Using the residual

$$
E(y,u)=\dot y-F(y,u)
$$

and the same minus-sign convention as the rest of the background notes,

$$
\mathcal L(y,u,p)
= J(y,u)
-\int_{0}^{T}p^{\mathsf T}E(y,u)\thinspace\mathrm{d}t.
$$

Expanding the residual gives

$$
\mathcal L
= \Phi(y(T))
+\int_{0}^{T}
\left[
\ell(y,u)
+p^{\mathsf T}F(y,u)
-p^{\mathsf T}\dot y
\right]
\mathrm{d}t.
$$

This motivates the Hamiltonian

$$
\boxed{
H(y,u,p)
\coloneqq
\ell(y,u)+p^{\mathsf T}F(y,u).
}
$$

The Lagrangian can now be read as

$$
\mathcal L
= \Phi(y(T))
+\int_{0}^{T}
\left[
H(y,u,p)-p^{\mathsf T}\dot y
\right]
\mathrm{d}t.
$$

### 7.1 State and costate equations come from the same variation as before

Variation with respect to $p$ recovers

$$
\dot y=F(y,u)=H_{p}.
$$

For the state variation, integrate the term containing $\delta\dot y$ by parts in time. Since the initial state is fixed, $\delta y(0)=0$. The terminal contribution gives

$$
p(T)=\Phi_{y}(y(T)),
$$

while the interior contribution gives

$$
-\dot p=H_{y}(y,u,p).
$$

Thus the familiar forward state/backward adjoint pair becomes

$$
\begin{aligned}
\dot y&=H_{p},\\
-\dot p&=H_{y},\\
p(T)&=\Phi_{y}(y(T)).
\end{aligned}
$$

The word **costate** used in optimal-control theory refers to the same mathematical role that we have called the adjoint or equality-constraint multiplier.

### 7.2 The control condition can be stronger than $H_{u}=0$

If the control is unconstrained and the Hamiltonian is smooth, stationarity gives

$$
H_{u}(y,u,p)=0.
$$

If the admissible control set $U_{\mathrm{ad}}$ is convex, the first-order condition becomes

$$
H_{u}(y,\bar u,p)(v-\bar u)
\geq0
\qquad
\forall v\in U_{\mathrm{ad}}.
$$

Under the standard hypotheses of the Pontryagin principle, a general admissible control set leads to a pointwise optimization of the Hamiltonian:

$$
\boxed{
H(\bar y(t),\bar u(t),\bar p(t))
= \min_{v\in U_{\mathrm{ad}}}
H(\bar y(t),v,\bar p(t))
}
$$

for almost every $t$, under the present minimization/sign convention.

This is one reason the Hamiltonian viewpoint is useful: it separates the expensive global propagation of the state and costate from a control condition that may be local in time, and sometimes local in space as well.

Different books define the Hamiltonian or the Lagrangian with the opposite sign and consequently state a **maximum** principle instead. The invariant content is the same after the costate sign is changed consistently. In particular, NMOPT Lecture 17 uses the opposite plus-sign Lagrangian convention and presents the corresponding Hamiltonian condition as a maximization statement. When comparing its formulas with the convention used here, change the adjoint sign consistently rather than mixing the two sign systems.

### 7.3 A scalar linear-quadratic example makes the Hamiltonian condition concrete

Consider

$$
\dot y=ay+bu,
$$

with

$$
J(y,u)
= \frac{q_{T}}{2}(y(T)-y_{T})^{2}
+\int_{0}^{T}
\left(
\frac{q}{2}y^{2}
+\frac{r}{2}u^{2}
\right)
\mathrm{d}t,
$$

where $r\gt0$. The Hamiltonian is

$$
H(y,u,p)
= \frac{q}{2}y^{2}
+\frac{r}{2}u^{2}
+p(ay+bu).
$$

The state equation is recovered from

$$
H_{p}=ay+bu.
$$

The costate equation is

$$
-\dot p
= H_{y}
= qy+ap,
$$

with terminal condition

$$
p(T)=q_{T}(y(T)-y_{T}).
$$

For an unconstrained control,

$$
H_{u}=ru+bp=0,
$$

so the Hamiltonian minimizer is

$$
\boxed{
u=-\frac{b}{r}p.
}
$$

If instead

$$
u_{a}\leq u\leq u_{b},
$$

the pointwise minimization gives

$$
\boxed{
u
= P_{[u_{a},u_{b}]}
\left(-\frac{b}{r}p\right).
}
$$

This is the dynamical analogue of the projection formulas from [06 · Constrained optimization](06-constrained-optimization.md) and the control constraints discussed in [07 · PDE-constrained optimization](07-pde-constrained-optimization.md). The state trajectory determines the costate, while the costate in turn determines the locally optimal control. The three objects are coupled even though the Hamiltonian control update itself is pointwise.

### 7.4 What the Pontryagin principle adds conceptually

The smooth KKT derivation says that an optimal trajectory satisfies state feasibility, adjoint stationarity, and a control optimality condition. Pontryagin's principle packages these into a dynamical structure and emphasizes that the control is chosen by optimizing the Hamiltonian along the state–costate trajectory.

That viewpoint is particularly useful when:

- the control set is not described conveniently by smooth equality/inequality constraints;
- the control acts pointwise and a local Hamiltonian minimization can be solved explicitly;
- bang-bang or switching controls appear;
- one wants algorithms that update the control through the Hamiltonian rather than through an assembled reduced Hessian.

The full theorem is substantially richer than the derivation above. A rigorous treatment must address measurability, endpoint constraints, constraint qualifications, possible abnormal multipliers, and the precise function spaces for the state and control. Those details are deliberately deferred to the optimal-control sources.

### 7.5 The same idea extends formally to PDEs

For a PDE evolution

$$
y_{t}=\mathcal F(y,u),
$$

the Hamiltonian becomes a functional involving the appropriate spatial pairing,

$$
\mathcal H(y,u,p)
= \ell(y,u)
+\langle p,\mathcal F(y,u)\rangle.
$$

The derivatives $\mathcal H_{y}$ and $\mathcal H_{u}$ are then function-space covectors. The resulting state–adjoint–control structure is the same one already derived through the PDE Lagrangian; the Hamiltonian notation highlights its dynamical organization rather than introducing a different optimization problem.

## 8. Newton–KKT, SQP, and SQH are related local-model viewpoints

Second-order algorithms for nonlinear PDE control are often described with several names. It helps to separate the underlying approximations.

### 8.1 Classical SQP linearizes the constraints and quadratizes the Lagrangian

Consider the abstract equality-constrained problem

$$
\min_{w} J(w)
\qquad
\text{subject to}
\qquad
E(w)=0,
$$

where $w$ may already collect a PDE state and control. Let

$$
\mathcal L(w,p)=J(w)-\langle p,E(w)\rangle.
$$

At an iterate $(w_{k},p_{k})$, an exact-Hessian SQP step $d$ solves the local quadratic problem

$$
\begin{aligned}
\min_{d}\quad&
J'(w_{k})[d]
+\frac{1}{2}
\mathcal L_{ww}''(w_{k},p_{k})[d,d],\\
\text{subject to}\quad&
E(w_{k})+E'(w_{k})d=0.
\end{aligned}
$$

The nonlinear PDE constraint has been replaced by its **linearized state equation**, while the objective curvature is represented by the Hessian of the Lagrangian.

The KKT equations of this quadratic subproblem are exactly the Newton linearization of the original nonlinear KKT system when the same Hessian is used. Thus

```text
Newton applied to nonlinear KKT equations
                ↕
exact-Hessian SQP subproblem
```

are two views of the same local mechanism.

In a PDE problem, $d$ usually contains a state increment and a control increment. Solving one SQP step therefore means solving another structured PDE-constrained quadratic problem, often a large saddle-point system.

### 8.2 Reduced Newton eliminates the state; all-at-once SQP keeps it

A reduced Newton method first eliminates the state and applies Newton to

$$
j(u)=J(S(u),u).
$$

Its Hessian action can be formed through incremental state and adjoint solves, as shown in [07a · Further notes on PDE-constrained optimization](07a-further-pde-constrained-optimization-notes.md).

An all-at-once SQP method instead keeps $(y,u)$ as independent primal variables and imposes the linearized PDE constraint inside the quadratic subproblem. Neither formulation is universally superior. The trade-off is the familiar one:

- reduced methods use repeated PDE solves but optimize in the smaller control space;
- all-at-once methods expose a larger block system but make the coupling explicit and can exploit KKT preconditioners.

Globalization, inexactness, and constraint handling remain essential in a practical SQP method and are not developed here.

### 8.3 SQH builds the local model through the Hamiltonian

A **Sequential Quadratic Hamiltonian** method keeps the state–costate organization explicit and constructs a local regularized model of the Hamiltonian with respect to the control.

A typical outer iteration has the form

```text
current control
      ↓
state solve
      ↓
adjoint/costate solve
      ↓
quadratic or regularized local Hamiltonian model
      ↓
local/pointwise control update
      ↓
trial state solve and acceptance test
```

This resembles SQP because both methods replace the nonlinear problem by a simpler second-order local model. The difference is organizational: SQP is naturally described as a quadratic **constrained optimization subproblem** in primal increments, whereas SQH is naturally described as a local **Hamiltonian control problem** along the current state–costate pair.

When the Hamiltonian separates pointwise in the control, the SQH update can be extremely local. The NMOPT course exploits this for a Poisson control example: state and adjoint remain global finite-element solves, but the control update is obtained from local Hamiltonian optimization and then projected back into the discrete control space.

The exact equivalence or difference between a particular SQH scheme, reduced Newton, and SQP depends on the chosen Hamiltonian approximation, regularization, and globalization. A complete convergence treatment therefore belongs to the dedicated sources rather than to this orientation.

## 9. More complicated PDEs reveal which parts of the machinery are truly reusable

The best way to read a new PDE optimal-control problem is to ask two questions:

1. what is new in the **state equation** before optimization is added?
2. how does that new structure propagate into the sensitivity, adjoint, and KKT systems?

The examples below are not full derivations. They give enough of the equations to make the additional difficulty concrete.

### 9.1 Stokes control: a saddle-point PDE inside a saddle-point optimization problem

For an incompressible viscous flow, the steady Stokes state can be written schematically as

$$
\begin{aligned}
-\nu\Delta v+\nabla\pi&=f+Bu
&&\text{in }\Omega,\\
\nabla\cdot v&=0
&&\text{in }\Omega.
\end{aligned}
$$

The state is already a pair:

$$
(v,\pi),
$$

velocity and pressure. The pressure acts as a multiplier enforcing incompressibility. A conforming finite-element discretization therefore needs compatible velocity and pressure spaces satisfying an inf-sup condition, and the state matrix has the mixed form

$$
\begin{bmatrix}
A & D^{\mathsf T}\\
D & 0
\end{bmatrix},
$$

where $D$ is the discrete divergence coupling. This symbol is kept distinct from the control operator $B$ in the state equation.

Optimization adds an adjoint velocity and adjoint pressure, plus the control. The outer optimality system therefore contains a **PDE saddle point inside the optimization saddle point**.

What survives unchanged is the logic:

```text
state solve
→ linearized state action
→ adjoint solve
→ control derivative
```

What becomes harder is the linear algebra: every state or adjoint application is itself a block mixed solve, pressure nullspaces or normalization may matter, and preconditioners must respect both incompressibility and optimization coupling.

### 9.2 Navier–Stokes control adds nonlinearity to the mixed structure

The steady incompressible Navier–Stokes equations add convection:

$$
\begin{aligned}
-\nu\Delta v
+(v\cdot\nabla)v
+\nabla\pi
&=f+Bu,\\
\nabla\cdot v&=0.
\end{aligned}
$$

Now the control-to-state map is nonlinear. Linearizing about a current velocity $v$ produces an Oseen-type operator containing

$$
(v\cdot\nabla)\delta v
+(\delta v\cdot\nabla)v.
$$

The adjoint contains the transpose of these convection terms, so it is not obtained by simply reusing the forward operator. Several complications arrive at once:

- nonlinear state solves may have only local uniqueness;
- the reduced objective may be nonconvex;
- the Jacobian is generally nonsymmetric;
- state and adjoint solves remain mixed velocity–pressure systems;
- Newton/SQP methods need globalization when the initial iterate is not already close to a solution.

The advanced sources treat both steady and unsteady Navier–Stokes control. Unsteady boundary control combines this nonlinear mixed structure with the forward/backward time organization from Sections 1–3 and with the boundary-control issues from [07a · Further notes on PDE-constrained optimization](07a-further-pde-constrained-optimization-notes.md).

### 9.3 Wave control is evolution without parabolic smoothing

The wave equation

$$
y_{tt}-c^{2}\Delta y=Bu+f
$$

requires both

$$
y(0)=y_{0},
\qquad
y_{t}(0)=v_{0}.
$$

It can be rewritten as a first-order system by introducing the velocity

$$
w=y_{t},
$$

but the resulting dynamics remain hyperbolic. Compared with heat control:

- information propagates at finite speed;
- the natural energy tracks both displacement and velocity;
- there is no analogous diffusion-driven smoothing;
- the time interval and geometry matter strongly for controllability;
- stable time discretization has different requirements.

The adjoint still propagates backward and the discrete-adjoint principle still applies, but now to the discrete wave propagator rather than to a dissipative parabolic step.

### 9.4 Coupled reaction–diffusion/ODE systems add multiple state types and nonlinear coupling

The monodomain model used in cardiac electrophysiology is a representative coupled system:

$$
\begin{aligned}
v_{t}
-\nabla\cdot(\sigma\nabla v)
+I_{\mathrm{ion}}(v,w)
&=I_{e},\\
w_{t}&=g(v,w).
\end{aligned}
$$

Here $v$ is a spatially diffusing electrical potential while $w$ is a recovery/gating variable governed locally by an ODE. The nonlinear functions $I_{\mathrm{ion}}$ and $g$ couple the two.

This introduces a product state space with components of different analytical character. The Jacobian contains PDE–ODE coupling blocks, the adjoint has one component for each state field, and time integration must evolve the coupled system consistently. If the optimization variable also includes a terminal time, as in time-optimal formulations, the time horizon itself becomes part of the optimization problem.

The example is useful because it shows that "more fields" does not merely mean a larger vector. Different blocks can live in different spaces and obey different evolution laws.

### 9.5 A useful complexity map

| State model | New difficulty before optimization | Consequence for the OCP |
| --- | --- | --- |
| scalar elliptic diffusion | coercive spatial solve | baseline state/adjoint pair |
| Stokes | incompressibility and mixed inf-sup structure | nested saddle-point blocks and pressure variables |
| Navier–Stokes | nonlinear convection plus mixed structure | nonlinear/nonconvex OCP, nonsymmetric Jacobians, Newton/SQP |
| parabolic diffusion | first-order time evolution | Bochner spaces, forward state, backward adjoint |
| wave equation | second-order hyperbolic evolution | two initial data, energy propagation, different controllability/stability |
| reaction–diffusion/ODE coupling | heterogeneous state blocks and nonlinear coupling | product spaces, coupled Jacobians/adjoints, multiphysics time stepping |

The point of the table is not that every future PDE fits one row. It is that the reusable optimization layer sits on top of a state operator whose own mathematical structure must first be understood. The more complicated the state problem becomes, the more important it is to preserve explicit spaces, pairings, block roles, and transpose actions.

## 10. Error analysis for an OCP couples state, adjoint, and control errors

[03 · Finite elements](03-finite-elements.md) separated finite-element approximation from the exact variational problem. In an optimal-control problem there are now several coupled discrete objects:

$$
y,
\qquad
p,
\qquad
u.
$$

An error in the state affects the objective and the adjoint source. An error in the adjoint affects the computed control derivative. The control error then changes the state again.

For a linear-quadratic problem, one therefore studies estimates of the schematic form

$$
\lVert y-y_{h}\rVert_{Y}
+\lVert p-p_{h}\rVert_{P}
+\lVert u-u_{h}\rVert_{U}
\leq
\text{approximation terms},
$$

under problem-dependent regularity and stability assumptions.

The exact norms, rates, and coupling terms depend on:

- the state and control spaces;
- whether the control is discretized explicitly;
- regularity of the state and adjoint;
- control constraints;
- whether DTO and OTD produce the same discrete optimality system.

### 10.1 A posteriori estimates support adaptive OCP discretization

A posteriori estimators use the computed discrete solution to estimate where the current mesh or discretization is inadequate. For an OCP, useful indicators may contain residual information from the state equation, adjoint equation, and control optimality condition.

This can drive adaptive refinement toward regions important for the **optimization objective**, not merely regions where the forward PDE state has a large local error.

Developing reliable estimators requires substantial finite-element error theory. Manzoni–Quarteroni–Salsa Sections 6.10–6.11 provide a dedicated linear-quadratic treatment; [03a · Further finite-element notes](03a-further-finite-element-notes.md) and the NMPDE sources provide the underlying finite-element approximation background.

## 11. Shape optimization changes what the optimization variable is

In all previous sections the domain $\Omega$ was fixed. Shape optimization instead asks to vary the domain itself:

$$
\Omega
\mapsto
J(\Omega,y(\Omega)).
$$

Ordinary vector addition

$$
\Omega+h
$$

has no useful meaning. To define a derivative, introduce a deformation field $V$ and perturb the domain through maps such as

$$
T_{t}(x)=x+tV(x).
$$

The perturbed domain is

$$
\Omega_{t}=T_{t}(\Omega).
$$

A **shape derivative** studies

$$
\frac{\mathrm{d}}{\mathrm{d}t}
J(\Omega_{t})
\Big\rvert_{t=0}.
$$

The PDE state must also be transported between changing domains before it can be differentiated consistently. Surface geometry, normal variations, and tangential calculus then enter naturally.

The adjoint method remains valuable because it can again eliminate expensive state sensitivities from the final derivative. However, the derivative now measures sensitivity to a **domain deformation**, not to an additive control field on a fixed function space.

This is why shape optimization should not be presented as just another control example. Manzoni–Quarteroni–Salsa devote a separate chapter to the subject, including domain deformations, shape derivatives, optimality conditions, and numerical approximation. That level of machinery is beyond a self-contained appendix to the current project background.

## 12. Choosing the next theory from the problem feature

The extensions can be organized by asking what feature of the original [07 · PDE-constrained optimization](07-pde-constrained-optimization.md) model has changed.

| If the new problem introduces… | The next mathematical tool is usually… |
| --- | --- |
| nonlinear state dependence | implicit-function/sensitivity theory, second derivatives, Newton or SQP |
| boundary actuation | trace spaces and boundary control operators |
| coefficient or parameter identification | nonlinear parameter-to-state maps and regularization |
| time dependence | Bochner spaces, evolution equations, forward/backward adjoints |
| sparse or switching controls | subdifferentials, proximal or semismooth methods |
| pointwise state inequalities | stronger constraint theory and possibly measure-valued multipliers |
| obstacle/contact state equations | variational inequalities and generalized derivatives |
| fluid constraints | mixed spaces and nested saddle-point solvers |
| domain design | shape calculus |
| adaptive accuracy requirements | OCP-specific a posteriori error analysis |

This table is not a classification of isolated topics. Each row starts from the same question:

> what part of the state–control–adjoint framework no longer behaves like the simple elliptic source-control model?

That question is a useful way to navigate a full optimal-control text without treating every advanced chapter as an unrelated new subject.

## Where to go next

- **For the current project:** this file is an extension map rather than a prerequisite checkpoint. Return to the [background index](README.md) or to the `nmopt` manual once you have the orientation you need.
- **If you are continuing through the software background:** [08 · Modern C++ for scientific code](08-modern-cpp-for-scientific-code.md) begins the implementation route, followed by [08a · Scientific software workflow: shell, compilation, CMake, and tests](08a-scientific-software-workflow.md) and [09 · The deal.II finite-element workflow](09-dealii-finite-element-workflow.md).
- **For deeper mathematical study:** use the topic-specific sources below. Evolution well-posedness, state constraints, generalized derivatives, Pontryagin principles, error estimates, and shape calculus each require enough additional theory that a full treatment is better learned in its native source.

## References and further reading

Publication details, stable links, and access notes are collected in the [Background reference catalogue](references.md).

- **Andrea Manzoni, Alfio Quarteroni, and Sandro Salsa, *Optimal Control of Partial Differential Equations*.**
  - Chapters 7–8 for linear evolution OCPs, including parabolic problems, time-dependent Stokes, the wave equation, and discrete DTO/OTD.
  - Chapter 5, §5.12 for state constraints.
  - Chapter 6, §§6.10–6.11 for a priori and a posteriori error analysis of linear-quadratic OCPs.
  - Chapter 9, especially §§9.6 and 9.8, for nonlinear numerical approximation, Newton/SQP methods, and control constraints.
  - Chapter 10 for advanced Navier–Stokes and cardiac-electrophysiology applications.
  - Chapter 11 for shape optimization.
- **Luca Heltai, NMOPT course.**
  - Lecture 10, **From ODE-Constrained Optimization to Bochner Spaces**, for the transition to evolution function spaces.
  - Lecture 12, **Time Discretization for Parabolic Optimal Control Problems**, for implicit Euler, discrete Lagrangians, discrete adjoints, and DTO/OTD in time.
  - Lecture 13, **deal.II Laboratory: Forward-Backward Projected Gradient for Parabolic Control**, for a concrete forward-state/backward-adjoint implementation.
  - Lecture 16, **Nonsmooth PDE-Constrained Optimization**, for sparse $L^{1}$ control, state constraints, variational-inequality states, complementarity, and semismooth ideas.
  - Lecture 17, **Pontryagin Principle and Sequential Quadratic Hamiltonian Methods**, for the Hamiltonian viewpoint and SQH.
  - Lecture 18, **Implementing SQH for a Poisson Control Problem**, for the finite-element implementation of a callback-based SQH method with a nonsmooth control term.
- **Juan Carlos De los Reyes, *Numerical PDE-Constrained Optimization*.**
  - A natural continuation for nonsmooth equations, semismooth Newton methods, and PDE-optimization algorithms beyond the smooth linear-quadratic setting.
