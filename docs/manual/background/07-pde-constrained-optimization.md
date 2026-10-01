# 07 · PDE-constrained optimization

**Background navigation:** [Index](README.md) \
Previous: [06 · Constrained optimization](06-constrained-optimization.md) \
Next: [08 · Modern C++ for scientific code](08-modern-cpp-for-scientific-code.md) \
Further notes: [07a · Further notes on PDE-constrained optimization](07a-further-pde-constrained-optimization-notes.md) \
Extensions: [07b · Extensions of PDE-constrained optimization](07b-extensions-of-pde-constrained-optimization.md)

## Purpose

A PDE-constrained optimization problem combines two tasks that have so far been treated separately: a differential equation determines a **state**, while an optimization problem chooses a **control** so that the resulting state has desirable properties. The control cannot be optimized independently of the PDE, because changing the control changes the state, and the objective usually depends on both.

This chapter develops the mathematical structure needed to understand that coupling. The main ideas are the control-to-state map, reduced objectives, direct sensitivities, adjoint equations, reduced derivatives, constrained controls, all-at-once optimality systems, and the distinction between discretize-then-optimize and optimize-then-discretize.

The central example is a distributed control of the Poisson equation. It is deliberately simple enough that every step can be derived explicitly, while still exhibiting the same state–adjoint–control structure that appears in much larger PDE optimization problems.

The chapter assumes the weak-PDE, duality, optimization, and finite-element background developed earlier. It does **not** attempt a general existence theory for optimal controls, a complete theory of nonlinear PDE-constrained optimization, or the analysis of state constraints, time-dependent controls, boundary controls of low regularity, or nonsmooth objectives. A short nonlinear bridge is included because it clarifies which parts of the state–adjoint derivation are special to linear problems and which are not.

## Before you start

You should be comfortable with:

- [01 · Function spaces and weak PDEs](01-function-spaces-and-weak-pdes.md), especially weak formulations and Lax–Milgram;
- [02 · Duality, derivatives, and adjoints](02-duality-derivatives-and-adjoints.md), especially Fréchet derivatives, adjoints, and Riesz maps;
- [05 · Unconstrained numerical optimization](05-unconstrained-optimization.md);
- [06 · Constrained optimization](06-constrained-optimization.md) for variational inequalities, projections, and KKT ideas;
- [03 · Finite elements](03-finite-elements.md) for the discretization sections.

[04 · Numerical linear algebra for PDE and optimization systems](04-numerical-linear-algebra.md) is useful when we discuss coupled KKT systems and work accounting.

## What you will be able to do

After this chapter, you should be able to:

- formulate a PDE-constrained optimization problem using state, control, objective, and residual notation;
- explain why the state problem must be well posed before a control-to-state map can be used;
- distinguish states, controls, fixed parameters, forcing data, and observation data;
- derive a reduced objective by eliminating the state;
- derive the linearized state or sensitivity equation;
- explain why direct sensitivities become expensive when the control space is large;
- derive an adjoint equation from the chain rule and an adjoint identity;
- derive the reduced derivative without explicitly computing the state sensitivity for every control direction;
- distinguish the reduced derivative from a control-space gradient under a chosen metric;
- state the constrained first-order condition as a variational inequality and recover the box-projection formula;
- compare reduced and all-at-once formulations;
- derive the state–adjoint–control optimality system from a Lagrangian;
- relate the adjoint sign convention to the conventional KKT multiplier sign;
- explain the difference between discretize-then-optimize and optimize-then-discretize and why they may or may not commute;
- derive the finite-element matrices for a simple distributed-control problem;
- account for the dominant state, adjoint, and trial-state solves in an optimization algorithm;
- verify a reduced derivative with finite differences and adjoint consistency checks.

## Roadmap

We begin by writing a PDE-constrained problem as an objective subject to a residual equation and identifying the roles played by states, controls, and data. The state equation then leads to the control-to-state map, which lets us eliminate the state and obtain a reduced optimization problem.

Differentiating the reduced objective first produces a direct sensitivity equation. That derivation is important because the adjoint method does not replace the chain rule; it reorganizes it. We then introduce the adjoint so that the state sensitivity disappears from the final derivative formula, and we specialize the result to the Poisson tracking problem.

The second half of the chapter compares reduced and all-at-once viewpoints, derives the coupled first-order system from the Lagrangian, and then asks what happens when optimization and finite-element discretization are performed in different orders. We finish with a concrete finite-element realization, computational work accounting, derivative checks, and a scope frontier.

## 1. The anatomy of a PDE-constrained optimization problem

Let

$$
y\in Y
$$

be a state and

$$
u\in U
$$

be a control. A broad class of problems can be written as

$$
\begin{aligned}
\min_{y,u}\quad & J(y,u),\\
\text{subject to}\quad & E(y,u)=0.
\end{aligned}
$$

The objective is a scalar functional

$$
J:Y\times U\to\mathbb R,
$$

while the PDE is represented by a residual

$$
E:Y\times U\to Z^{\ast}.
$$

Here $Z$ is a test space. Writing the residual in $Z^{\ast}$ means that

$$
E(y,u)=0
$$

is shorthand for

$$
\langle E(y,u),v\rangle_{Z^{\ast},Z}=0
\qquad
\forall v\in Z.
$$

This is the natural language for weak PDE formulations.

### 1.1 State and control play different mathematical roles

The state is not normally a freely selectable quantity. Once the control and the fixed data are given, the PDE determines which states are feasible.

The control is the quantity the optimization algorithm is allowed to change. It may be:

- a distributed source;
- a boundary value or boundary flux;
- a coefficient in the PDE;
- a forcing amplitude;
- a geometric or material parameter.

Not every quantity appearing in a PDE model is therefore a control. It is useful to distinguish four roles:

- **state** – determined by the PDE;
- **control** – changed by the optimizer;
- **parameter** – fixed during one optimization problem but possibly varied between studies;
- **data** – prescribed forcing, boundary data, observations, or targets.

The same mathematical quantity could play different roles in different problems. A diffusion coefficient, for example, may be fixed data in one model and an unknown control in an identification problem.

### 1.2 A running distributed-control problem

Let $\Omega\subset\mathbb R^{d}$ be a bounded Lipschitz domain. Consider

$$
\begin{aligned}
\min_{y,u}\quad
J(y,u)
&\coloneqq
\frac{1}{2}\lVert y-y_{d}\rVert_{L^{2}(\Omega)}^{2} +
\frac{\beta}{2}\lVert u\rVert_{L^{2}(\Omega)}^{2},\\
\text{subject to}\quad
-\Delta y
&=
f+u
\quad\text{in }\Omega,\\
y&=0
\quad\text{on }\partial\Omega,
\end{aligned}
$$

with $\beta\gt0$.

The roles are:

- $y$ – state;
- $u$ – distributed control;
- $f$ – fixed forcing;
- $y_{d}$ – desired or target state;
- $\beta$ – regularization weight.

The tracking term rewards states close to $y_{d}$. The control term discourages arbitrarily large controls and, in this linear-quadratic problem, also contributes strict convexity.

The weak state problem is: find

$$
y\in V\coloneqq H_{0}^{1}(\Omega)
$$

such that

$$
a(y,v) = (f+u,v)_{L^{2}(\Omega)}
\qquad
\forall v\in V,
$$

where

$$
a(y,v)\coloneqq
\int_{\Omega}\nabla y\cdot\nabla v\thinspace\mathrm{d}x.
$$

Equivalently, define

$$
E:V\times L^{2}(\Omega)\to V^{\ast}
$$

by

$$
\langle E(y,u),v\rangle\coloneqq
a(y,v)-(f+u,v)_{L^{2}(\Omega)}.
$$

Then the weak PDE is simply

$$
E(y,u)=0
\quad\text{in }V^{\ast}.
$$

This residual convention will be used throughout the chapter.

### 1.3 Observation operators generalize state tracking

The objective need not compare the whole state directly with a target in $L^{2}(\Omega)$. A useful general form is

$$
J(y,u) = \frac{1}{2}\lVert Cy-z_{d}\rVert_{O}^{2} +
\frac{\beta}{2}\lVert u\rVert_{U}^{2},
$$

where

$$
C:Y\to O
$$

is an observation operator into an observation space $O$. Depending on the problem, $C$ may represent observation on a subregion, a trace on part of the boundary, an average, or another measured quantity.

If $O$ is a Hilbert space with Riesz map $R_{O}$, the state derivative of the tracking term is

$$
D_{y}J(y,u) = C^{\ast}R_{O}(Cy-z_{d})\in Y^{\ast}.
$$

The running example corresponds to the simplest case where $O=L^{2}(\Omega)$ and $C$ is the natural embedding of the state into $L^{2}(\Omega)$. Introducing $C$ explicitly is useful because it separates the PDE from the question of what part or feature of the state is actually observed.

## 2. The state equation comes before the optimization

The notation

$$
y=S(u)
$$

looks innocent, but it contains an analytical assumption: for a chosen control $u$, the state equation must determine a state in a sufficiently reliable way.

For the Poisson example, Lax–Milgram gives exactly what is needed. For each

$$
u\in L^{2}(\Omega),
$$

there is a unique

$$
y\in H_{0}^{1}(\Omega)
$$

satisfying the weak state equation, and the state depends continuously on the forcing.

We may therefore define the **control-to-state map**

$$
S:U\to Y,
\qquad
u\mapsto y=S(u).
$$

For the running problem,

$$
U=L^{2}(\Omega),
\qquad
Y=H_{0}^{1}(\Omega).
$$

Because the PDE is linear but contains the fixed forcing $f$, $S$ is generally affine rather than linear. We can write

$$
S(u)=y_{f}+S_{0}u,
$$

where $y_{f}$ solves the equation with $u=0$ and $S_{0}$ maps a control perturbation to the corresponding state perturbation.

### 2.1 Why uniqueness matters for reduction

Suppose the same control could produce several unrelated states. Then the notation $S(u)$ would not identify a unique object, and the expression

$$
J(S(u),u)
$$

would be ambiguous unless an additional rule selected one state.

This does not mean every PDE-constrained method requires a globally defined control-to-state map. All-at-once methods can sometimes work in settings where reduction is inconvenient. But a **reduced** formulation assumes that the state can be selected as a function of the control at the points visited by the algorithm.

### 2.2 A small nonlinear bridge

The same idea survives beyond linear PDEs. Let

$$
E:Y\times U\to Z^{\ast}
$$

be continuously Fréchet differentiable between Banach spaces, and suppose

$$
E(\bar y,\bar u)=0.
$$

If the partial derivative

$$
D_{y}E(\bar y,\bar u):Y\to Z^{\ast}
$$

is a bounded isomorphism, the implicit-function theorem gives, locally around $\bar u$, a differentiable map

$$
y=S(u)
$$

with

$$
E(S(u),u)=0.
$$

Differentiating that identity gives the same sensitivity equation used below.

This is enough nonlinear theory for the present chapter. In genuinely nonlinear PDE control, existence of states, global uniqueness, differentiability of $S$, nonconvexity of the reduced objective, second-order conditions, and globalization of nonlinear state solves all require additional analysis. The important point here is that the state–sensitivity–adjoint structure is not restricted to linear equations.

## 3. Eliminate the state: the reduced objective

Once the state can be written as

$$
y=S(u),
$$

the constrained problem becomes

$$
\min_{u\in U} j(u),
$$

where

$$
j(u)\coloneqq J(S(u),u).
$$

For the running example,

$$
j(u) = \frac{1}{2}\lVert S(u)-y_{d}\rVert_{L^{2}(\Omega)}^{2} +
\frac{\beta}{2}\lVert u\rVert_{L^{2}(\Omega)}^{2}.
$$

This is the **reduced formulation**.

The number of optimization variables has been reduced from the state–control pair $(y,u)$ to the control $u$ alone. The PDE has not disappeared computationally. Evaluating $j(u)$ still requires solving

$$
E(y,u)=0
$$

for $y$.

That distinction is fundamental:

```text
mathematical variable elimination
        does not mean
computational PDE elimination
```

The reduced formulation hides the state solve inside the map $S$.

### 3.1 Why the running linear-quadratic problem has a unique minimizer

A complete existence theory for PDE-constrained optimization is outside this chapter, but the running problem is a useful place to connect the finite-dimensional optimization ideas with the direct method discussed in the companion [05a · Further notes on unconstrained optimization](05a-further-unconstrained-optimization-notes.md).

The regularization term gives

$$
j(u)
\geq
\frac{\beta}{2}\lVert u\rVert_{L^{2}(\Omega)}^{2}.
$$

Thus a minimizing sequence cannot escape to infinite $L^{2}$ norm. Since $L^{2}(\Omega)$ is reflexive, a bounded minimizing sequence has a weakly convergent subsequence. The affine continuous state map and the squared norms make the reduced objective weakly lower semicontinuous, so the weak limit is a minimizer.

Moreover, the term

$$
\frac{\beta}{2}\lVert u\rVert_{L^{2}(\Omega)}^{2}
$$

is strictly convex when $\beta\gt0$. Adding the convex tracking term preserves strict convexity, so the minimizer is unique.

This argument is special in its simplicity. For nonlinear state equations the reduced map can be nonconvex even when the original tracking functional looks convex, and first-order conditions then cease to characterize a global minimizer.

## 4. Differentiate the state equation: direct sensitivities

Let

$$
h\in U
$$

be a control perturbation. The induced state perturbation is

$$
s\coloneqq S'(u)[h].
$$

The letter $s$ is used here for a **sensitivity state**; $p$ will later be reserved for the adjoint.

Start from the identity

$$
E(S(u),u)=0.
$$

Differentiate in the direction $h$. The chain rule gives

$$
D_{y}E(y,u)[s] + D_{u}E(y,u)[h] = 0,
$$

where

$$
y=S(u).
$$

Therefore the sensitivity satisfies the **linearized state equation**

$$
\boxed{
D_{y}E(y,u)[s] = -D_{u}E(y,u)[h].
}
$$

If $D_{y}E(y,u)$ is invertible, this may be written formally as

$$
S'(u)[h] = -D_{y}E(y,u)^{-1}D_{u}E(y,u)[h].
$$

The inverse notation describes the mathematics; a numerical method normally solves a linearized PDE rather than constructing an inverse operator.

### 4.1 Sensitivity equation for the Poisson problem

For

$$
\langle E(y,u),v\rangle =
a(y,v)-(f+u,v)_{L^{2}(\Omega)},
$$

we have

$$
\langle D_{y}E(y,u)[s],v\rangle =
a(s,v)
$$

and

$$
\langle D_{u}E(y,u)[h],v\rangle =
-(h,v)_{L^{2}(\Omega)}.
$$

Hence the sensitivity equation is

$$
a(s,v) =
(h,v)_{L^{2}(\Omega)}
\qquad
\forall v\in V.
$$

In strong notation this is

$$
-\Delta s=h,
\qquad
s=0\text{ on }\partial\Omega.
$$

Because the state equation is linear, this sensitivity does not depend on the current state or control. For a nonlinear PDE, the linearized operator would generally depend on $(y,u)$.

## 5. Direct differentiation of the reduced objective

Apply the chain rule to

$$
j(u)=J(S(u),u).
$$

For a control perturbation $h$,

$$
j'(u)[h] = D_{y}J(y,u)[s] + D_{u}J(y,u)[h],
$$

where

$$
s=S'(u)[h].
$$

For the running objective,

$$
D_{y}J(y,u)[s] = (y-y_{d},s)_{L^{2}(\Omega)}
$$

and

$$
D_{u}J(y,u)[h] = \beta(u,h)_{L^{2}(\Omega)}.
$$

Therefore

$$
j'(u)[h] = (y-y_{d},s)_{L^{2}(\Omega)} +
\beta(u,h)_{L^{2}(\Omega)}.
$$

This formula is correct. The difficulty is that $s$ depends on the direction $h$ through a PDE solve.

### 5.1 One direction is cheap; every direction is not

If we need the derivative only in one chosen direction $h$, direct sensitivity analysis can be perfectly reasonable:

1. solve the state equation for $y$;
2. solve one linearized state equation for $s=S'(u)[h]$;
3. evaluate the two terms in $j'(u)[h]$.

But an optimization method usually wants the entire derivative

$$
j'(u)\in U^{\ast},
$$

or enough information to construct a gradient. If $U_{h}$ is a discrete control space with basis

$$
\lbrace\psi_{1},\ldots,\psi_{N_{u}}\rbrace,
$$

constructing the derivative coefficient by coefficient with direct sensitivities would require solving

$$
D_{y}E(y,u)[s_{k}] = -D_{u}E(y,u)[\psi_{k}],
\qquad
k=1,\ldots,N_{u}.
$$

That is potentially one linearized PDE solve per control degree of freedom.

For a distributed control, $N_{u}$ can be comparable to the number of mesh unknowns. The adjoint method avoids this direction-by-direction cost.

## 6. The adjoint reorganizes the same chain rule

The purpose of the adjoint is to remove $s$ from the final expression for the reduced derivative.

We have

$$
D_{y}J(y,u)\in Y^{\ast}
$$

and

$$
D_{y}E(y,u):Y\to Z^{\ast}.
$$

Its adjoint is

$$
D_{y}E(y,u)^{\ast}:Z\to Y^{\ast},
$$

defined by

$$
\left\langle D_{y}E(y,u)[\delta y],p\right\rangle_{Z^{\ast},Z} =
\left\langle D_{y}E(y,u)^{\ast}p,\delta y\right\rangle_{Y^{\ast},Y}.
$$

Choose the adjoint state

$$
p\in Z
$$

to solve

$$
\boxed{
D_{y}E(y,u)^{\ast}p = D_{y}J(y,u).
}
$$

Now use the adjoint equation in the state-dependent term of the reduced derivative:

$$
\begin{aligned}
D_{y}J(y,u)[s]
&=
\left\langle D_{y}E(y,u)^{\ast}p,s\right\rangle\\
&=
\left\langle D_{y}E(y,u)[s],p\right\rangle.
\end{aligned}
$$

The sensitivity equation gives

$$
D_{y}E(y,u)[s] = -D_{u}E(y,u)[h].
$$

Hence

$$
\begin{aligned}
D_{y}J(y,u)[s]
&=
-\left\langle D_{u}E(y,u)[h],p\right\rangle\\
&=
-\left\langle D_{u}E(y,u)^{\ast}p,h\right\rangle.
\end{aligned}
$$

Substituting into the chain rule yields

$$
\boxed{
j'(u) = D_{u}J(y,u) - D_{u}E(y,u)^{\ast}p,
}
$$

with

$$
\boxed{
D_{y}E(y,u)^{\ast}p = D_{y}J(y,u).
}
$$

The sensitivity $s$ has disappeared. One adjoint solve supplies the state-mediated contribution to the derivative for **all** control directions at once.

### 6.1 Why this is a reverse accumulation

The direct method propagates a chosen control perturbation forward:

```text
control direction h
        ↓
linearized state solve
        ↓
state sensitivity s
        ↓
objective change
```

The adjoint method starts from the scalar objective derivative and propagates its state-side covector backward:

```text
state objective covector D_y J
        ↓
adjoint solve
        ↓
adjoint p
        ↓
D_u E^* p
        ↓
control covector j'(u)
```

This is the same forward-versus-reverse distinction seen in [02 · Duality, derivatives, and adjoints](02-duality-derivatives-and-adjoints.md), now applied to a PDE solve.

## 7. The running Poisson adjoint and reduced derivative

For the Poisson residual,

$$
\langle D_{y}E[s],v\rangle = a(s,v).
$$

The adjoint equation asks for $p\in V$ such that

$$
a(\delta y,p) = (y-y_{d},\delta y)_{L^{2}(\Omega)}
\qquad
\forall \delta y\in V.
$$

Because the Poisson bilinear form is symmetric,

$$
a(\delta y,p)=a(p,\delta y),
$$

so this can also be written

$$
a(p,v) = (y-y_{d},v)_{L^{2}(\Omega)}
\qquad
\forall v\in V.
$$

In strong notation,

$$
-\Delta p=y-y_{d}
\quad\text{in }\Omega,
\qquad
p=0
\quad\text{on }\partial\Omega.
$$

The equality of the state and adjoint differential operators here is a consequence of symmetry. For advection, mixed systems, and many nonlinear PDEs, the adjoint operator is different from the state operator.

For a control direction $h$,

$$
\left\langle D_{u}E[h],p\right\rangle = -(h,p)_{L^{2}(\Omega)}.
$$

Therefore

$$
\begin{aligned}
j'(u)[h]
&=
\beta(u,h)_{L^{2}(\Omega)}-
\left[-(h,p)_{L^{2}(\Omega)}\right]\\
&=
(\beta u+p,h)_{L^{2}(\Omega)}.
\end{aligned}
$$

Thus

$$
\boxed{
j'(u)[h] = (\beta u+p,h)_{L^{2}(\Omega)}.
}
$$

This formula replaces one sensitivity solve per direction with one adjoint solve.

## 8. The derivative is still not automatically the gradient

The reduced derivative is intrinsically a covector,

$$
j'(u)\in U^{\ast}.
$$

A gradient depends on a chosen control-space inner product. Let

$$
G:U\to U^{\ast}
$$

be the corresponding Riesz map. The gradient $g\in U$ is defined by

$$
Gg=j'(u).
$$

Equivalently,

$$
j'(u)[h]
=(g,h)_{G}
\qquad
\forall h\in U.
$$

### 8.1 The $L^{2}$ gradient in the running example

If the control metric is the ordinary $L^{2}(\Omega)$ inner product, then

$$
j'(u)[h] = (\beta u+p,h)_{L^{2}(\Omega)}
$$

shows directly that

$$
\boxed{
\nabla_{L^{2}}j(u)=\beta u+p.
}
$$

If a different metric is chosen, the derivative does not change. Only its primal representative changes:

$$
g=G^{-1}j'(u).
$$

This matters numerically. Choosing a metric can change the scale and smoothness of optimization directions even though the underlying objective and first-order derivative remain the same.

### 8.2 Regularization and optimization metric are different choices

The term

$$
\frac{\beta}{2}\lVert u\rVert_{L^{2}}^{2}
$$

belongs to the objective. It changes which control is optimal.

The metric $G$ belongs to the optimization geometry. It changes how the derivative is represented as a primal direction.

Both may involve the same mass operator after finite-element discretization, but they play different mathematical roles.

## 9. What one reduced optimization iteration actually costs

Suppose an optimization algorithm is at a control $u_{k}$ and needs both the objective and derivative.

A typical reduced evaluation is:

```text
u_k
 ↓
state solve
E(y_k,u_k)=0
 ↓
objective J(y_k,u_k)
 ↓
state objective derivative D_y J
 ↓
adjoint solve
D_y E(y_k,u_k)^* p_k = D_y J(y_k,u_k)
 ↓
reduced derivative
j'(u_k)=D_u J-D_u E^* p_k
```

The expensive operations are usually the PDE solves, not the vector algebra around them.

For a first-order reduced method, an accepted iteration often requires approximately:

- one state solve;
- one adjoint solve;
- one metric solve or inverse metric application if the chosen gradient representation requires it.

The exact count depends on the optimization method.

### 9.1 Line searches can require extra state solves

A trial control

$$
u_{\mathrm{trial}}=u_{k}+\alpha d_{k}
$$

cannot usually be assigned an objective value without solving its state equation. Therefore each rejected line-search trial may cost another state solve.

An important implementation pattern follows:

```text
trial control
    ↓
state solve + objective
    ↓
acceptance test
    ├── reject → no adjoint needed
    └── accept → compute adjoint and derivative
```

There is no reason to pay for an adjoint at a trial point that is rejected using objective information alone.

This is why optimization iteration counts can be misleading for PDE problems. Two algorithms with the same number of accepted iterations may perform very different numbers of state and adjoint solves.

### 9.2 Newton and Hessian actions may add more PDE solves

A reduced Hessian action generally contains differentiated state and adjoint information. Matrix-free Newton–CG methods may therefore require additional linearized state and adjoint solves inside the inner Krylov iteration.

The details depend strongly on the problem and are beyond this chapter. The main accounting principle is simple:

> measure work in expensive PDE and linearized-PDE solves, not only in outer optimization iterations.

## 10. Control constraints change stationarity, not the state equation

Suppose the control is required to lie in a nonempty closed convex set

$$
U_{\mathrm{ad}}\subset U.
$$

The reduced problem becomes

$$
\min_{u\in U_{\mathrm{ad}}}j(u).
$$

At an optimal control $\bar u$, the first-order condition is the variational inequality

$$
\boxed{
j'(\bar u)[v-\bar u]\geq0
\qquad
\forall v\in U_{\mathrm{ad}}.
}
$$

Using the adjoint formula,

$$
\left\langle D_{u}J(\bar y,\bar u) -
D_{u}E(\bar y,\bar u)^{\ast}\bar p, v-\bar u\right\rangle\geq0.
$$

The state and adjoint equations have the same form as before. Only the control stationarity condition changes.

### 10.1 Box constraints in the running example

Let

$$
U_{\mathrm{ad}} = \left\lbrace
u\in L^{2}(\Omega):
 u_{a}(x)\leq u(x)\leq u_{b}(x)
\text{ a.e. in }\Omega
\right\rbrace.
$$

Then

$$
(\beta\bar u+\bar p,v-\bar u)_{L^{2}(\Omega)}
\geq0
\qquad
\forall v\in U_{\mathrm{ad}}.
$$

For $\beta\gt0$, the projection characterization from [06 · Constrained optimization](06-constrained-optimization.md) gives

$$
\boxed{
\bar u = P_{U_{\mathrm{ad}}}\left(-\frac{1}{\beta}\bar p\right).
}
$$

For pointwise box constraints this is simply clipping:

$$
\bar u(x) = \min\left\lbrace
u_{b}(x),
\max\left\lbrace
u_{a}(x),-\frac{1}{\beta}\bar p(x)\right\rbrace
\right\rbrace
$$

for almost every $x\in\Omega$.

The equivalent normal-cone and complementarity formulations were developed in [06 · Constrained optimization](06-constrained-optimization.md). Those forms are especially useful for active-set and primal-dual methods.

## 11. The Lagrangian produces the full optimality system

Reduction is not the only way to organize the first-order conditions. Introduce an adjoint or multiplier

$$
p\in Z
$$

and define the Lagrangian

$$
\boxed{
\mathcal L(y,u,p)\coloneqq J(y,u) -
\langle E(y,u),p\rangle_{Z^{\ast},Z}.
}
$$

The minus sign is a convention. Once it is chosen, the state, adjoint, and control signs follow from differentiation.

For an unconstrained control, stationarity of $\mathcal L$ with respect to the three variables gives the complete first-order system.

### 11.1 Variation with respect to the adjoint: state feasibility

For any $\delta p\in Z$,

$$
D_{p}\mathcal L(y,u,p)[\delta p] = -\langle E(y,u),\delta p\rangle.
$$

Requiring this to vanish for every $\delta p$ gives

$$
E(y,u)=0.
$$

### 11.2 Variation with respect to the state: adjoint equation

For any $\delta y\in Y$,

$$
\begin{aligned}
D_{y}\mathcal L(y,u,p)[\delta y]
&=
D_{y}J(y,u)[\delta y] -
\langle D_{y}E(y,u)[\delta y],p\rangle\\
&=
\left\langle
D_{y}J(y,u)-D_{y}E(y,u)^{\ast}p,
\delta y
\right\rangle.
\end{aligned}
$$

Hence

$$
D_{y}E(y,u)^{\ast}p = D_{y}J(y,u).
$$

### 11.3 Variation with respect to the control: reduced stationarity

For any $h\in U$,

$$
D_{u}\mathcal L(y,u,p)[h] =
\left\langle D_{u}J(y,u)-D_{u}E(y,u)^{\ast}p,h\right\rangle.
$$

For an unconstrained control,

$$
D_{u}J(y,u)-D_{u}E(y,u)^{\ast}p
=0
\quad\text{in }U^{\ast}.
$$

This is exactly the reduced derivative evaluated at a feasible state and its adjoint.

Thus the reduced and all-at-once derivations are not competing first-order theories. They organize the same equations differently.

## 12. Reduced and all-at-once formulations

The **reduced** route eliminates the state at each control:

```text
control u
    ↓
state solve
    ↓
adjoint solve
    ↓
reduced derivative
    ↓
optimization update
```

The **all-at-once** route keeps state, control, and adjoint as simultaneous unknowns and solves the coupled first-order conditions:

```text
unknown (y,u,p)
        ↓
┌────────────────────────────┐
│ state equation             │
│ adjoint equation           │
│ control stationarity       │
└────────────────────────────┘
        ↓
coupled solve
```

A compact comparison is:

| Aspect | Reduced | All-at-once |
| --- | --- | --- |
| Optimization unknown | control | state, control, adjoint |
| State feasibility during outer iterations | enforced by each state solve | may be violated at intermediate iterates |
| Main expensive operation | repeated PDE/adjoint solves | larger coupled solves |
| Natural linear algebra | state and adjoint operators, reduced Hessian actions | block/KKT systems and block preconditioners |
| Useful when | a reliable PDE solver already exists | coupled structure can be exploited effectively |

Neither organization is universally superior. A mature PDE solver can make reduction very attractive. A well-designed block preconditioner can make an all-at-once system attractive, especially for problems where nested solves are expensive or tightly coupled.

## 13. Linear-quadratic problems reveal the KKT structure

For the running problem, the all-at-once first-order conditions are

$$
\begin{aligned}
-\Delta y &= f+u,\\
-\Delta p &= y-y_{d},\\
\beta u+p &= 0,
\end{aligned}
$$

with homogeneous Dirichlet conditions on $y$ and $p$.

The first equation enforces feasibility, the second is the adjoint equation, and the third is control stationarity.

After discretization these equations form a saddle-point system. The same structure can be understood already at the abstract optimization level.

Let

$$
x\coloneqq \begin{bmatrix}
y\\
u
\end{bmatrix}
$$

collect the primal variables, and write the equality constraint abstractly as

$$
Dx=d.
$$

A quadratic objective has the form

$$
\varphi(x) = \frac{1}{2}\langle Qx,x\rangle -
\langle c,x\rangle + c_{0}.
$$

The conventional quadratic-programming Lagrangian is often written

$$
\widehat{\mathcal L}(x,\lambda) = \varphi(x) +
\langle\lambda,Dx-d\rangle.
$$

Its first-order conditions are

$$
Qx-c+D^{\ast}\lambda=0,
$$

and

$$
Dx-d=0.
$$

In finite-dimensional coordinates this yields the canonical KKT system

$$
\begin{bmatrix}
Q & D^{\mathsf T}\\
D & 0
\end{bmatrix}
\begin{bmatrix}
x\\
\lambda
\end{bmatrix} =
\begin{bmatrix}
c\\
d
\end{bmatrix}.
$$

### 13.1 The KKT multiplier and the PDE adjoint may differ by a sign

Our PDE Lagrangian convention is

$$
\mathcal L = J-\langle E,p\rangle,
$$

and for the equality residual

$$
E=Dx-d.
$$

The conventional KKT Lagrangian above uses a plus sign. The two conventions agree when

$$
\boxed{
\lambda=-p.
}
$$

There is no contradiction. Flipping the multiplier sign changes the appearance of the off-diagonal blocks but not the primal optimum. [06a · Further notes on constrained optimization](06a-further-constrained-optimization-notes.md) discusses these sign conventions in more detail.

### 13.2 Why the KKT system is a saddle-point problem

The multiplier has no quadratic objective term, so the multiplier–multiplier block is zero. The resulting operator is not positive definite on the full product space even when the primal objective is strictly convex.

This is why PDE optimality systems naturally lead to the symmetric-indefinite and block-preconditioning ideas introduced in [04 · Numerical linear algebra for PDE and optimization systems](04-numerical-linear-algebra.md).

## 14. Discretize then optimize and optimize then discretize

There are two logically different routes from a continuous PDE-constrained problem to a finite-dimensional optimality system.

### 14.1 Discretize then optimize: DTO

In **discretize then optimize**, first choose finite-dimensional spaces and build a discrete objective and discrete PDE residual:

$$
J_{h}(y_{h},u_{h}),
\qquad
E_{h}(y_{h},u_{h})=0.
$$

Then derive the reduced derivative, adjoint, or KKT system of this finite-dimensional optimization problem.

Schematically,

```text
continuous state equation + objective
        ↓ discretize
discrete objective J_h and residual E_h
        ↓ differentiate / optimize
discrete adjoint and optimality system
```

The discrete adjoint is therefore, by construction, the transpose of the derivative of the **actual discrete residual** under the declared discrete pairings.

### 14.2 Optimize then discretize: OTD

In **optimize then discretize**, first derive the continuous optimality system:

```text
continuous objective + PDE
        ↓ differentiate / optimize
continuous state–adjoint–stationarity system
        ↓ discretize each equation
discrete optimality system
```

The adjoint PDE is obtained at the continuous level and then receives its own numerical discretization.

### 14.3 The two routes can coincide

For the linear-quadratic elliptic problem, a compatible Galerkin discretization can make the two routes produce the same algebraic equations. In particular, if the state and adjoint use compatible trial/test spaces and the discrete objective uses the matching inner products, the discrete adjoint matrix is the transpose required by DTO.

That agreement is important, but it is a property of the discretization, not part of the definitions of DTO and OTD.

### 14.4 The two routes can differ

Optimization and discretization need not commute. Differences can appear when, for example:

- the state and adjoint use different spaces;
- the state equation is stabilized and the continuous adjoint is stabilized independently;
- a Petrov–Galerkin method uses different trial and test spaces;
- quadrature changes the actual discrete residual or objective;
- boundary terms are approximated differently;
- nonlinear terms are linearized or discretized in different orders;
- constraints or transformations are imposed differently at the continuous and discrete levels.

A useful diagnostic question is therefore not merely

> What is the adjoint PDE?

but also

> Adjoint of **which** continuous or discrete operator, under **which** pairings?

### 14.5 Quadrature is part of the distinction in an implemented method

An ideal Galerkin discretization is often written as though all element integrals were exact. An implementation normally evaluates them by quadrature.

If DTO differentiates the arrays actually assembled with a chosen quadrature rule, the resulting discrete adjoint is tied to that numerical realization. An OTD derivation that discretizes the continuous adjoint with different quadrature choices need not give the exact algebraic transpose.

This does not imply that one route is inherently correct and the other incorrect. It means that the numerical realization must be understood precisely when comparing them. [03 · Finite elements](03-finite-elements.md) and [03a · Further finite-element notes](03a-further-finite-element-notes.md) discuss quadrature and variational crimes in more detail.

## 15. Finite-element discretization of the running problem

Choose finite-dimensional spaces

$$
V_{h}\subset H_{0}^{1}(\Omega)
$$

for the state and adjoint, and

$$
U_{h}\subset L^{2}(\Omega)
$$

for the control.

The spaces need not have the same dimension. Let

$$
\lbrace\phi_{1},\ldots,\phi_{N_{y}}\rbrace
$$

be a basis of $V_{h}$ and

$$
\lbrace\psi_{1},\ldots,\psi_{N_{u}}\rbrace
$$

be a basis of $U_{h}$.

Expand

$$
y_{h} = \sum_{j=1}^{N_{y}} y_{j}\phi_{j},
\qquad
u_{h} = \sum_{k=1}^{N_{u}} u_{k}\psi_{k}.
$$

### 15.1 The discrete state equation

Testing with $\phi_{i}$ gives

$$
\sum_{j=1}^{N_{y}}y_{j}a(\phi_{j},\phi_{i}) =
(f,\phi_{i})_{L^{2}} +
\sum_{k=1}^{N_{u}}u_{k}(\psi_{k},\phi_{i})_{L^{2}}.
$$

Define

$$
A_{ij}\coloneqq a(\phi_{j},\phi_{i}),
$$

$$
B_{ik}\coloneqq (\psi_{k},\phi_{i})_{L^{2}(\Omega)},
$$

and

$$
f_{i}\coloneqq (f,\phi_{i})_{L^{2}(\Omega)}.
$$

Then

$$
\boxed{
Ay=f+Bu.
}
$$

The coupling matrix

$$
B\in\mathbb R^{N_{y}\times N_{u}}
$$

is rectangular when the state and control spaces have different dimensions.

### 15.2 The discrete objective

Define the state and control mass matrices

$$
(M_{y})_{ij}\coloneqq (\phi_{j},\phi_{i})_{L^{2}(\Omega)},
$$

and

$$
(M_{u})_{k\ell}\coloneqq (\psi_{\ell},\psi_{k})_{L^{2}(\Omega)}.
$$

If the target $y_{d}$ is not itself represented by a coefficient vector in $V_{h}$, define the load-like vector

$$
(q_{d})_{i}\coloneqq (y_{d},\phi_{i})_{L^{2}(\Omega)}.
$$

Then, up to a constant independent of $(y,u)$,

$$
J_{h}(y,u) = \frac{1}{2}y^{\mathsf T}M_{y}y - q_{d}^{\mathsf T}y +
\frac{\beta}{2}u^{\mathsf T}M_{u}u + \text{constant}.
$$

This form avoids pretending that an arbitrary continuous target automatically has a nodal coefficient vector.

### 15.3 The discrete adjoint equation

The discrete residual is

$$
E_{h}(y,u)\coloneqq Ay-Bu-f.
$$

With the same minus-sign Lagrangian convention,

$$
\mathcal L_{h}(y,u,p) = J_{h}(y,u) -
p^{\mathsf T}(Ay-Bu-f).
$$

State stationarity gives

$$
M_{y}y-q_{d}-A^{\mathsf T}p=0,
$$

or

$$
\boxed{
A^{\mathsf T}p=M_{y}y-q_{d}.
}
$$

### 15.4 Discrete control stationarity

The direct control derivative is

$$
\beta M_{u}u.
$$

Because

$$
D_{u}E_{h}=-B,
$$

the adjoint contribution is

$$
-B^{\mathsf T}p.
$$

Therefore

$$
\boxed{
\beta M_{u}u+B^{\mathsf T}p=0.
}
$$

The complete three-block system is

$$
\begin{bmatrix}
A & 0 & -B\\
M_{y} & -A^{\mathsf T} & 0\\
0 & B^{\mathsf T} & \beta M_{u}
\end{bmatrix}
\begin{bmatrix}
y\\
p\\
u
\end{bmatrix} =
\begin{bmatrix}
f\\
q_{d}\\
0
\end{bmatrix}.
$$

The same equations may be reordered and the multiplier sign changed to obtain the conventional symmetric KKT layout when the primal quadratic form and pairings are symmetric.

### 15.5 Same mesh does not imply same algebraic role

If state and control use the same $L^{2}$ basis, the coupling matrix $B$ may numerically equal a mass matrix. That numerical equality does not make the roles identical:

- $B$ maps control coefficients into the state residual;
- $M_{u}$ represents the control inner product in the objective;
- $M_{y}$ represents the state tracking inner product.

Likewise, the control metric used by an optimization method may also be represented by a mass matrix without being the same mathematical object as the regularization term.

### 15.6 A cellwise control naturally gives a different space

A common choice is a continuous finite-element state together with a discontinuous, cellwise-constant control. Then

$$
N_{u}\neq N_{y}
$$

in general, and the coupling matrix is genuinely rectangular.

This is not an inconvenience to be hidden. It reflects the fact that state and control are different functions with different approximation spaces.

## 16. Reduced algebra from the discrete state equation

Suppose $A$ is invertible. Then

$$
y=A^{-1}(f+Bu).
$$

Substituting this into the discrete objective gives a reduced quadratic function of $u$.

The reduced derivative can be obtained by the state–adjoint sequence

```text
u
 ↓
Ay=f+Bu
 ↓
y
 ↓
A^T p=M_y y-q_d
 ↓
p
 ↓
β M_u u+B^T p
```

For the linear-quadratic problem, the reduced Hessian is

$$
H_{\mathrm{red}} = \beta M_{u} +
B^{\mathsf T}A^{-\mathsf T}M_{y}A^{-1}B.
$$

The formula is useful conceptually, but there is usually no need to assemble the dense operator represented by the inverse factors. A Hessian action can be computed by PDE solves and sparse operator applications.

This expression also makes positivity visible. For any control increment $d$,

$$
\begin{aligned}
d^{\mathsf T}H_{\mathrm{red}}d
&=
\beta d^{\mathsf T}M_{u}d +
(A^{-1}Bd)^{\mathsf T}M_{y}(A^{-1}Bd).
\end{aligned}
$$

If $M_{u}$ is positive definite and $\beta\gt0$, then the reduced Hessian is positive definite. This is the discrete reflection of strict convexity of the linear-quadratic reduced problem.

## 17. Derivative and adjoint verification

An adjoint formula is easy to derive with one sign error and still obtain plausible-looking numerical output. Verification should therefore be part of the numerical workflow.

### 17.1 Check the state solve first

Before testing derivatives, check that the computed state satisfies

$$
\lVert E(y,u)\rVert
$$

to the intended solve tolerance. A derivative check cannot diagnose the adjoint cleanly if the base state is substantially infeasible.

### 17.2 Check the adjoint identity

Choose representative perturbations $\delta y$ and adjoint seeds $p$. Verify numerically that

$$
\left\langle
D_{y}E(y,u)[\delta y],p
\right\rangle
\approx
\left\langle
D_{y}E(y,u)^{\ast}p,\delta y
\right\rangle.
$$

This isolates the transpose implementation from the rest of the optimization derivation.

### 17.3 Check the reduced directional derivative

Choose a control direction $h$ and compare

$$
j'(u)[h]
$$

with finite differences. A centered check is

$$
D_{t}\coloneqq \frac{j(u+t h)-j(u-t h)}{2t}.
$$

For a sufficiently smooth problem and sufficiently accurate state solves,

$$
D_{t}-j'(u)[h] = O(t^{2})
$$

until floating-point error and solve tolerances dominate.

A first-order Taylor remainder is also useful:

$$
R(t)\coloneqq \left\lvert j(u+t h)-j(u)-t j'(u)[h]\right\rvert.
$$

For a correct derivative,

$$
R(t)=O(t^{2})
$$

in the asymptotic range.

### 17.4 A failed derivative test does not identify the bug by itself

Possible causes include:

- the state solve is too loose;
- the adjoint equation uses the wrong sign or boundary condition;
- the transpose action is inconsistent with the forward derivative;
- the objective derivative is wrong;
- the control metric has been confused with the derivative pairing;
- quadrature or constrained-coordinate treatment differs between forward and adjoint code;
- the finite-difference step is too large or too small.

Testing residuals, adjoint identities, and the complete reduced derivative separately makes these failures much easier to localize.

## 18. Scope frontier

This chapter develops the first-order state–adjoint structure needed for the standard PDE-constrained problems considered here. Several nearby subjects are important but belong beyond the prerequisite path. The closest generalizations – nonlinear elliptic state equations, reduced Hessian actions, boundary controls, and coefficient identification – are developed in [07a · Further notes on PDE-constrained optimization](07a-further-pde-constrained-optimization-notes.md). Broader extensions that need substantially new machinery are oriented in [07b · Extensions of PDE-constrained optimization](07b-extensions-of-pde-constrained-optimization.md).

We stop before:

- a general direct-method existence theory for nonlinear optimal controls;
- detailed second-order sufficient conditions for nonlinear PDE problems;
- nonlinear state equations beyond the local implicit-function viewpoint;
- state constraints and measure-valued multipliers;
- nonsmooth objectives such as $L^{1}$ sparsity terms;
- time-dependent PDE control and Bochner-space formulations;
- low-regularity Dirichlet boundary control and very-weak solutions;
- shape optimization;
- Pontryagin maximum principles and Hamiltonian methods;
- specialized reduced-Hessian and KKT preconditioners beyond the orientation given in [04a · Further numerical linear algebra for PDE-constrained optimization](04a-further-linear-algebra-for-pde-optimization.md).

The key transferable structure is already visible:

```text
control
  ↓
state equation
  ↓
objective
  ↓
adjoint equation
  ↓
reduced derivative
```

or, equivalently,

```text
state + control + adjoint
          ↓
full first-order optimality system
```

Everything later changes details of the spaces, equations, discretization, or solver strategy without removing this basic architecture.

## Where to go next

- **Default continuation for understanding `nmopt`:** return to the concept manual. The most useful sequence is:
  - [manual 05 · Reduced state–adjoint formulation](../concepts/05-reduced-state-adjoint-formulation.md), which maps the state/adjoint derivation to the project's reduced formulation services;
  - [manual 06 · Reduced optimization methods](../concepts/06-reduced-optimization-methods.md), which consumes reduced values and derivatives in line-search and trust-region algorithms;
  - [manual 07 · Optimality systems and KKT](../concepts/07-optimality-systems-and-kkt.md), which keeps state, adjoint, and control together in all-at-once systems;
  - [manual 08 · Complementarity and PDAS](../concepts/08-complementarity-and-pdas.md), which adds box constraints and active-set iterations.
- **For nearby mathematical generalizations:** [07a · Further notes on PDE-constrained optimization](07a-further-pde-constrained-optimization-notes.md) develops nonlinear state equations, reduced Hessian actions, Newton/SQP structure, boundary controls, coefficient identification, and observation operators.
- **For broader extensions beyond the current project:** [07b · Extensions of PDE-constrained optimization](07b-extensions-of-pde-constrained-optimization.md) gives a guided introduction to time-dependent control, nonsmooth objectives, state constraints, Hamiltonian methods, more complicated PDEs, OCP error analysis, and shape optimization.
- **For the software route:** continue with [08 · Modern C++ for scientific code](08-modern-cpp-for-scientific-code.md), then [08a · Scientific software workflow: shell, compilation, CMake, and tests](08a-scientific-software-workflow.md) as needed, and finally [09 · The deal.II finite-element workflow](09-dealii-finite-element-workflow.md). [manual overview · Numerical realization](../overview/numerical-realization.md) is the bridge back from those prerequisites to `nmopt`'s coefficient vectors, pairings, metrics, and solve services.

## References and further reading

Publication details, stable links, and access notes are collected in the [Background reference catalogue](references.md).

- **Andrea Manzoni, Alfio Quarteroni, and Sandro Salsa, *Optimal Control of Partial Differential Equations*.**
  - Chapter 5, especially §§5.1–5.3, for linear-quadratic elliptic control, reduced formulations, adjoints, multipliers, box constraints, projections, and the abstract variational framework.
  - Chapter 6, especially §§6.1–6.2, for reduced/all-at-once viewpoints and DTO/OTD.
  - Chapter 6, §6.6 for the all-at-once saddle-point viewpoint.
- **Luca Heltai, NMOPT course.**
  - Lecture 4 for the reduced elliptic state–adjoint formulation.
  - Lecture 5 for control constraints.
  - Lecture 6 for finite-element discretization.
  - Lectures 7–8 for all-at-once and more general optimality systems.
- **NMPDE course material.**
  - The weak-PDE and finite-element ingredients used here are developed in the sources cited by [01 · Function spaces and weak PDEs](01-function-spaces-and-weak-pdes.md) and [03 · Finite elements](03-finite-elements.md).
- **General numerical-optimization references.**
  - The optimization algorithms applied to the reduced problem are covered by the sources cited in [05 · Unconstrained numerical optimization](05-unconstrained-optimization.md) and [06 · Constrained optimization](06-constrained-optimization.md).
