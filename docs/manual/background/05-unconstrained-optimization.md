# 05 · Unconstrained numerical optimization

**Background navigation:** [Index](README.md) \
Previous: [04 · Numerical linear algebra for PDE and optimization systems](04-numerical-linear-algebra.md) \
Next: [06 · Constrained optimization](06-constrained-optimization.md) \
Companion: [05a · Further notes on unconstrained optimization](05a-further-unconstrained-optimization-notes.md)

## Purpose

An optimization problem asks for a point that makes an objective function small. Even in finite dimensions, that description leaves several numerical questions unanswered. How do we turn a derivative into a direction? How far should we move along that direction? When is curvature information worth its cost? How do we keep a fast local method from taking a disastrous step far from the solution? And how do we decide that the computation has done enough work?

This chapter develops the smooth unconstrained optimization background needed later for reduced PDE-constrained optimization. The emphasis is not on cataloguing algorithms. It is on the common structure behind the methods that recur in the project: steepest descent, line searches, nonlinear conjugate gradients, Newton and Newton–CG ideas, BFGS and L-BFGS, and trust regions.

The central distinction from [02 · Duality, derivatives, and adjoints](02-duality-derivatives-and-adjoints.md) remains in force: the derivative is a covector, while a gradient is a primal vector obtained after choosing an inner product or metric. Most formulas are first written in Euclidean coordinates because that is where the classical algorithms are easiest to see. We then state the metric-aware interpretation whenever it matters.

The scope is deliberately narrower than a general numerical-optimization course. We derive the local models and update formulas that explain the algorithms, and we state the practical conditions that make them useful. A companion, [05a · Further notes on unconstrained optimization](05a-further-unconstrained-optimization-notes.md), collects the mathematical existence/convergence background and a few deeper algorithmic derivations that are useful but not prerequisites for following the main chapter. We do not develop derivative-free or stochastic methods, or attempt a comprehensive convergence theory. Constraints are the subject of [06 · Constrained optimization](06-constrained-optimization.md).

## Before you start

This chapter assumes:

- multivariable calculus, including gradients and Hessians in ordinary Euclidean coordinates;
- [02 · Duality, derivatives, and adjoints](02-duality-derivatives-and-adjoints.md), especially the distinction between derivatives and gradients;
- elementary linear algebra;
- [04 · Numerical linear algebra for PDE and optimization systems](04-numerical-linear-algebra.md) for SPD systems, conditioning, and conjugate gradients.

## What you will be able to do

After this chapter, you should be able to:

- distinguish local and global minimizers and explain when convexity turns a local first-order condition into a global one;
- derive first- and second-order necessary conditions and a standard second-order sufficient condition;
- distinguish a derivative covector from a gradient under a chosen metric;
- recognize and verify a descent direction;
- explain why a locally sensible direction still needs globalization;
- derive exact line search for an SPD quadratic and interpret it as an idealized special case;
- apply and interpret Armijo, Wolfe, and strong Wolfe conditions;
- explain why steepest descent slows down on poorly conditioned quadratics;
- distinguish linear conjugate gradients from nonlinear conjugate-gradient optimization;
- write the Fletcher–Reeves and Polak–Ribière+ updates and explain why restarts are useful;
- derive the Newton step from a quadratic model and explain the Newton–CG idea;
- derive the secant equation from gradient differences;
- explain the BFGS curvature condition and the role of L-BFGS storage;
- understand the purpose of the L-BFGS two-loop recursion;
- interpret trust-region acceptance through actual versus predicted reduction;
- combine gradient, step, objective, and work information into sensible stopping criteria.

## Roadmap

We begin by separating local from global minimization and by deriving the first- and second-order conditions that numerical methods try to satisfy. Convexity is introduced at this point because it is the main situation in which first-order stationarity has a global meaning.

A descent method then separates into two decisions: choose a direction and choose a step length. The SPD quadratic is used when a concrete model is helpful because its geometry makes conditioning, exact line search, and conjugacy visible with almost no extra machinery.

From there we develop practical globalization through Armijo and Wolfe conditions, including the reasoning behind the inequalities rather than treating them as recipes. We then distinguish linear CG from nonlinear conjugate-gradient methods and motivate the Fletcher–Reeves and Polak–Ribière formulas from the quadratic conjugacy picture. Newton's method introduces exact curvature; quasi-Newton methods replace it with information inferred from accepted steps, leading to BFGS and L-BFGS.

Trust-region methods provide a second globalization philosophy: instead of committing to a direction and then shortening the step, they restrict where the local model is trusted. We finish with stopping criteria, scaling, and work accounting. The companion notes collect existence via compact level sets and coercivity, convergence rates of sequences, stronger convexity assumptions, and selected deeper derivations and large-scale variants.

## 1. Smooth unconstrained minimization: local, global, and stationary points

Consider

$$
\min_{x\in\mathbb R^{n}} f(x),
$$

with $f:\mathbb R^{n}\to\mathbb R$ differentiable.

A point $x_{\ast}$ is a **local minimizer** if there exists $r>0$ such that

$$
f(x_{\ast})\leq f(x)
\qquad
\text{whenever }\lVert x-x_{\ast}\rVert<r.
$$

It is a **strict local minimizer** if the inequality is strict for every nearby $x\neq x_{\ast}$. It is a **global minimizer** if

$$
f(x_{\ast})\leq f(x)
\qquad
\text{for every }x\in\mathbb R^{n}.
$$

These notions should not be conflated with stationarity. Numerical methods based on local derivatives normally seek a stationary point with properties consistent with a local minimum. Without additional structure, they cannot certify that this point is the global minimizer.

### 1.1 Convexity is the main bridge from local to global

A function $f$ is **convex** if

$$
f((1-t)x+ty)
\leq
(1-t)f(x)+t f(y)
$$

for every $x,y$ and $t\in[0,1]$. If the inequality is strict for $x\neq y$ and $t\in(0,1)$, the function is strictly convex.

For a differentiable convex function, the tangent hyperplane lies below the graph:

$$
f(y)
\geq
f(x)+\nabla f(x)^{\mathsf T}(y-x).
$$

This inequality explains why convexity matters so much in optimization. If

$$
\nabla f(x_{\ast})=0,
$$

then immediately

$$
f(y)\geq f(x_{\ast})
\qquad
\text{for every }y,
$$

so $x_{\ast}$ is a global minimizer. If $f$ is strictly convex, that minimizer is unique.

For a nonconvex objective, stationarity has only local meaning. A stationary point can be a minimum, a maximum, or a saddle point.

### 1.2 First-order stationarity is necessary at an unconstrained local minimum

Take any direction $d\in\mathbb R^{n}$ and define the one-dimensional restriction

$$
\phi(t):=f(x_{\ast}+td).
$$

If $x_{\ast}$ is a local minimizer, then $t=0$ is a local minimizer of $\phi$. Therefore

$$
\phi'(0)=0.
$$

By the chain rule,

$$
\phi'(0)
=
f'(x_{\ast})[d].
$$

Because the direction $d$ was arbitrary,

$$
f'(x_{\ast})=0
\qquad
\text{in }(\mathbb R^{n})^{\ast}.
$$

In Euclidean coordinates this is

$$
\nabla f(x_{\ast})=0.
$$

The argument is worth remembering: an unconstrained local minimum must look like a one-dimensional minimum along **every** line through the point.

### 1.3 Second-order necessary and sufficient conditions

Assume now that $f$ is twice differentiable near $x_{\ast}$. Along the same line,

$$
\phi''(0)
=
d^{\mathsf T}\nabla^{2}f(x_{\ast})d.
$$

A twice-differentiable one-dimensional function has nonnegative second derivative at a local minimum. Therefore every local minimizer satisfies

$$
d^{\mathsf T}\nabla^{2}f(x_{\ast})d\geq0
\qquad
\text{for every }d,
$$

so

$$
\nabla^{2}f(x_{\ast})
\quad\text{is positive semidefinite}.
$$

This is a **necessary** second-order condition. It is still not sufficient: if the Hessian is singular, higher-order terms may decide the nature of the stationary point. In one dimension, $x^{4}$, $-x^{4}$, and $x^{3}$ all have zero first and second derivative at the origin, but the origin is respectively a local minimum, a local maximum, and neither.

A standard sufficient condition is stronger. Suppose

$$
\nabla f(x_{\ast})=0
$$

and

$$
d^{\mathsf T}\nabla^{2}f(x_{\ast})d>0
\qquad
\text{for every nonzero }d.
$$

Taylor expansion gives, for a small displacement $s$,

$$
f(x_{\ast}+s)
=
f(x_{\ast})
+
\frac{1}{2}s^{\mathsf T}\nabla^{2}f(x_{\ast})s
+
o(\lVert s\rVert^{2}).
$$

Positive definiteness makes the quadratic term uniformly positive relative to $\lVert s\rVert^{2}$, while the remainder is smaller order. Hence sufficiently small nonzero $s$ satisfy

$$
f(x_{\ast}+s)>f(x_{\ast}),
$$

and $x_{\ast}$ is a strict local minimizer.

The companion notes develop the related existence questions, compact sublevel sets, minimizing sequences, strong convexity, and convergence rates. The main chapter now turns from characterizing minimizers to constructing iterations that try to reach them.

## 2. Derivative and gradient still depend on different choices

[02 · Duality, derivatives, and adjoints](02-duality-derivatives-and-adjoints.md) emphasized that the derivative is intrinsically a covector. Let

$$
r:=f'(x)\in X^{\ast}
$$

for a finite-dimensional space $X$.

If we choose an inner product represented by an SPD Riesz map

$$
G:X\to X^{\ast},
$$

the corresponding gradient is the primal vector $g\in X$ satisfying

$$
Gg=r.
$$

Thus

$$
g=G^{-1}r.
$$

The directional derivative along $d\in X$ is

$$
f'(x)[d]
=
\langle r,d\rangle
=
\langle Gg,d\rangle.
$$

In ordinary Euclidean coordinates, $G=I$, so the derivative coordinates and the gradient coordinates coincide. Classical optimization texts usually work in this setting, and most formulas below use

$$
g_{k}:=\nabla f(x_{k}).
$$

The underlying rule is nevertheless metric-dependent. When a later application uses a non-Euclidean metric, the primal gradient is obtained from the derivative before a direction such as $-g_{k}$ is formed.

## 3. A local model tells us which directions look promising

At the current iterate $x_{k}$, Taylor expansion gives

$$
f(x_{k}+s)
=
f(x_{k})
+
g_{k}^{\mathsf T}s
+
\frac{1}{2}s^{\mathsf T}H_{k}s
+
o(\lVert s\rVert^{2}),
$$

where

$$
g_{k}=\nabla f(x_{k}),
\qquad
H_{k}=\nabla^{2}f(x_{k}).
$$

The first-order model is

$$
f(x_{k})+g_{k}^{\mathsf T}s.
$$

A vector $d_{k}$ is a **descent direction** if

$$
g_{k}^{\mathsf T}d_{k}<0.
$$

Why is that enough to predict decrease? Set

$$
\phi(\alpha):=f(x_{k}+\alpha d_{k}).
$$

Then

$$
\phi'(0)=g_{k}^{\mathsf T}d_{k}<0.
$$

Therefore, for sufficiently small positive $\alpha$,

$$
f(x_{k}+\alpha d_{k})<f(x_{k}).
$$

The derivative tells us that a small step should work. It does **not** tell us that a unit step will work. This is the basic reason for globalization.

### 3.1 Direction choice and step choice are separate

A line-search method has the generic form

$$
x_{k+1}=x_{k}+\alpha_{k}d_{k}.
$$

There are two logically different decisions:

1. choose $d_{k}$ from local derivative/curvature information;
2. choose $\alpha_{k}>0$ so that the proposed movement is acceptable.

Steepest descent, nonlinear CG, Newton, and BFGS mainly differ in step 1. Armijo and Wolfe rules mainly address step 2.

Keeping these roles separate makes it much easier to understand why the same line search can globalize several different direction formulas.

## 4. Steepest descent is the reference method

In Euclidean geometry, choose

$$
d_{k}:=-g_{k}.
$$

Then

$$
g_{k}^{\mathsf T}d_{k}
=
-\lVert g_{k}\rVert_{2}^{2}<0
$$

whenever $g_{k}\ne0$.

Under a general metric $G$, the same statement becomes

$$
d_{k}=-G^{-1}f'(x_{k}),
$$

and

$$
f'(x_{k})[d_{k}]
=
-\lVert d_{k}\rVert_{G}^{2}.
$$

Steepest descent uses only first-order information. It is therefore cheap and robust enough to serve as a reference method, but it can converge very slowly when the objective has strongly different curvatures in different directions.

## 5. The SPD quadratic makes conditioning visible

Consider

$$
q(x)
=
\frac{1}{2}x^{\mathsf T}H x-b^{\mathsf T}x,
\qquad
H=H^{\mathsf T}>0.
$$

Its gradient is

$$
g(x)=Hx-b,
$$

so the unique stationary point satisfies

$$
Hx_{\ast}=b.
$$

Thus minimizing an SPD quadratic and solving an SPD linear system are the same mathematical problem viewed from two directions.

Let

$$
e_{k}:=x_{k}-x_{\ast}.
$$

Since $Hx_{\ast}=b$,

$$
g_{k}=He_{k}.
$$

The Hessian eigenvectors are the principal axes of the quadratic level sets. If the eigenvalues vary greatly, the level sets are elongated. The condition number

$$
\kappa_{2}(H)
=
\frac{\lambda_{\max}(H)}{\lambda_{\min}(H)}
$$

measures that disparity.

A steepest-descent step points normal to the current level set. On a long narrow ellipse, that normal tends to point across the narrow direction rather than along the valley toward the minimizer. Successive gradients therefore alternate across the valley: the familiar zig-zag pattern.

This is the optimization version of the conditioning phenomenon from [04 · Numerical linear algebra for PDE and optimization systems](04-numerical-linear-algebra.md). A poor Hessian condition number means very different local scales, and a single scalar step length cannot simultaneously be ideal for all eigendirections.

## 6. Exact line search is an illuminating idealization

Given a descent direction $d_{k}$, exact line search asks for

$$
\alpha_{k}
\in
\mathop{\mathrm{argmin}}_{\alpha\ge0}
f(x_{k}+\alpha d_{k}).
$$

For a general nonlinear objective, solving this one-dimensional optimization problem very accurately can cost more than it saves. But for the SPD quadratic it can be derived exactly.

Define

$$
\phi(\alpha)
:=
q(x_{k}+\alpha d_{k}).
$$

Then

$$
\phi'(\alpha)
=
d_{k}^{\mathsf T}
\left(
H(x_{k}+\alpha d_{k})-b
\right)
=
d_{k}^{\mathsf T}g_{k}
+
\alpha d_{k}^{\mathsf T}H d_{k}.
$$

Set the derivative to zero:

$$
\alpha_{k}
=
-\frac{g_{k}^{\mathsf T}d_{k}}
{d_{k}^{\mathsf T}H d_{k}}.
$$

Because $H$ is SPD, the denominator is positive. If $d_{k}$ is a descent direction, the numerator is negative, so $\alpha_{k}>0$.

For steepest descent, $d_{k}=-g_{k}$, hence

$$
\alpha_{k}
=
\frac{g_{k}^{\mathsf T}g_{k}}
{g_{k}^{\mathsf T}H g_{k}}.
$$

Exact line search is useful here because it exposes the geometry cleanly. In a nonlinear problem, practical line searches usually seek an **acceptable** step rather than the exact one-dimensional minimizer.

## 7. Globalization: why a good direction can still need a shorter step

Suppose $d_{k}$ is a descent direction. The first-order model predicts

$$
f(x_{k}+\alpha d_{k})
\approx
f(x_{k})
+
\alpha g_{k}^{\mathsf T}d_{k}.
$$

For very small $\alpha$, this is reliable. For a large $\alpha$, higher-order terms can dominate and the trial point may increase the objective.

A **globalization strategy** modifies or rejects steps so that an algorithm with good local behavior can also make progress from points farther from the minimizer.

Two major strategies will appear in this chapter:

- **line search:** choose a direction first, then regulate its length;
- **trust region:** choose a step by minimizing a local model only inside a region where that model is trusted.

The word *globalization* does not mean that the method finds a global minimizer. It means that local models are embedded in a strategy intended to improve convergence from a wider range of initial guesses.

## 8. Armijo backtracking asks for a controlled fraction of the predicted decrease

Let $d_{k}$ be a descent direction and define

$$
\phi(\alpha):=f(x_{k}+\alpha d_{k}).
$$

At the current point,

$$
\phi'(0)=g_{k}^{\mathsf T}d_{k}<0.
$$

The first-order Taylor model therefore predicts

$$
\phi(\alpha)
\approx
\phi(0)+\alpha\phi'(0).
$$

It would be unreasonable to demand that the nonlinear objective reproduce this linear prediction exactly. Instead choose

$$
c_{1}\in(0,1)
$$

and accept $\alpha>0$ when

$$
f(x_{k}+\alpha d_{k})
\leq
f(x_{k})
+
c_{1}\alpha g_{k}^{\mathsf T}d_{k}.
$$

Because $g_{k}^{\mathsf T}d_{k}<0$, the right-hand side is below the current objective value. More specifically, the linear model predicts a decrease

$$
-\alpha g_{k}^{\mathsf T}d_{k}>0,
$$

while Armijo asks the actual decrease to be at least the fraction $c_{1}$ of that prediction. A small $c_{1}$ is therefore permissive: it asks for definite progress without insisting that the local linear model remain quantitatively accurate over the whole step.

The condition is not arbitrary. Differentiability gives

$$
\phi(\alpha)
=
\phi(0)
+
\alpha\phi'(0)
+
o(\alpha).
$$

Subtract the Armijo right-hand side. We need

$$
(1-c_{1})\alpha\phi'(0)+o(\alpha)\leq0.
$$

The leading term is strictly negative because $c_{1}<1$ and $\phi'(0)<0$. Hence Armijo must hold for all sufficiently small positive $\alpha$. This is why backtracking from a trial step is a sensible algorithm rather than a blind heuristic.

### 8.1 Backtracking

A simple practical strategy is:

1. choose an initial trial step $\alpha_{0}>0$, often $1$;
2. choose a contraction factor $\beta\in(0,1)$;
3. while Armijo fails, replace $\alpha$ by $\beta\alpha$;
4. accept the first successful step.

In pseudocode:

```text
alpha = alpha_0
while f(x_k + alpha d_k) > f(x_k) + c1 alpha g_k^T d_k
    alpha = beta alpha
accept alpha
```

The preceding Taylor argument explains why this loop eventually reaches an acceptable local regime under the usual smoothness assumptions.

Armijo answers mainly the question **is the step too large to realize enough of the predicted decrease?** It does not prevent an unnecessarily tiny step, because a very small step usually satisfies the inequality almost automatically.

## 9. Wolfe conditions also prevent stopping too early

Armijo controls function decrease, but consider what an exact one-dimensional minimizer would satisfy. If

$$
\alpha_{\ast}>0
$$

is an interior minimizer of $\phi(\alpha)=f(x_{k}+\alpha d_{k})$, then

$$
\phi'(\alpha_{\ast})=0.
$$

At the starting point, by contrast,

$$
\phi'(0)=g_{k}^{\mathsf T}d_{k}<0.
$$

So a line search that stops almost immediately can satisfy Armijo while the objective is still descending nearly as steeply as it was at the start. Wolfe's curvature condition asks for evidence that we have moved far enough along the line for that negative slope to relax.

Choose

$$
0<c_{1}<c_{2}<1.
$$

The **Wolfe conditions** require Armijo sufficient decrease together with

$$
\nabla f(x_{k}+\alpha d_{k})^{\mathsf T}d_{k}
\geq
c_{2}g_{k}^{\mathsf T}d_{k}.
$$

Both quantities on the right are negative. The accepted directional derivative must therefore be closer to zero than the initial one by a controlled amount. This is the sense in which the curvature condition prevents an excessively short step.

The **strong Wolfe condition** replaces the curvature inequality by

$$
\left|
\nabla f(x_{k}+\alpha d_{k})^{\mathsf T}d_{k}
\right|
\leq
c_{2}
\left|
g_{k}^{\mathsf T}d_{k}
\right|.
$$

Now the accepted slope must be small in magnitude. Besides avoiding premature stopping, this also limits severe overshoot past a line minimizer. That extra control is particularly useful when the next direction reuses the previous one, as in nonlinear conjugate-gradient methods.

### 9.1 Why Armijo and Wolfe appear together

The two tests regulate different failures of the local model:

- Armijo asks whether the step produced enough objective decrease relative to the initial slope;
- Wolfe asks whether the step moved far enough for the slope along the line to change substantially;
- strong Wolfe also limits how large that slope may become with the opposite sign.

A simple backtracking implementation often checks only Armijo. More elaborate line searches seek both inequalities because curvature information is useful for nonlinear CG and for maintaining good quasi-Newton curvature pairs.

## 10. Linear CG and nonlinear CG are related but different algorithms

[04 · Numerical linear algebra for PDE and optimization systems](04-numerical-linear-algebra.md) introduced **linear conjugate gradients** for the SPD system

$$
Hx=b.
$$

Equivalently, linear CG minimizes the quadratic

$$
q(x)
=
\frac{1}{2}x^{\mathsf T}H x-b^{\mathsf T}x.
$$

Its directions satisfy exact $H$-conjugacy in exact arithmetic:

$$
d_{i}^{\mathsf T}H d_{j}=0,
\qquad i\ne j.
$$

With exact line searches, the method reaches the exact solution in at most $n$ steps in exact arithmetic.

For a nonlinear $f$, there is no single fixed Hessian $H$ governing every iteration. **Nonlinear conjugate gradient** keeps the short recurrence

$$
d_{k}
=
-g_{k}+\beta_{k}d_{k-1},
$$

but the exact linear-CG conjugacy relations no longer persist automatically.

The method should therefore be understood as carrying useful directional memory from one nonlinear step to the next, not as applying the linear CG theorem to a changing nonlinear objective.

## 11. Fletcher–Reeves and Polak–Ribière+ come from the quadratic conjugacy picture

Nonlinear CG keeps the recurrence

$$
d_{k}
=
-g_{k}+\beta_{k}d_{k-1}.
$$

The formulas for $\beta_{k}$ are easier to remember if we first recover one of them from the SPD quadratic. Suppose

$$
q(x)=\frac{1}{2}x^{\mathsf T}H x-b^{\mathsf T}x,
\qquad H=H^{\mathsf T}>0,
$$

and exact line search is used. Consecutive CG gradients are then orthogonal,

$$
g_{k}^{\mathsf T}g_{k-1}=0,
$$

and the new direction is required to be $H$-conjugate to the previous one:

$$
d_{k}^{\mathsf T}H d_{k-1}=0.
$$

Insert $d_{k}=-g_{k}+\beta_{k}d_{k-1}$ and solve for $\beta_{k}$:

$$
\beta_{k}
=
\frac{g_{k}^{\mathsf T}H d_{k-1}}
{d_{k-1}^{\mathsf T}H d_{k-1}}.
$$

For a quadratic objective the gradient change along the previous accepted step is exact:

$$
g_{k}-g_{k-1}
=
\alpha_{k-1}H d_{k-1}.
$$

Using this relation together with the orthogonality identities produced by exact line search gives

$$
\beta_{k}^{\mathrm{FR}}
=
\frac{g_{k}^{\mathsf T}g_{k}}
{g_{k-1}^{\mathsf T}g_{k-1}}.
$$

This is the **Fletcher–Reeves** coefficient. It can be read as retaining an amount of the previous direction proportional to the new squared gradient norm relative to the old one.

The **Polak–Ribière** coefficient keeps the gradient-difference information visible instead of reducing it using the exact quadratic orthogonality identities:

$$
\beta_{k}^{\mathrm{PR}}
=
\frac{
g_{k}^{\mathsf T}(g_{k}-g_{k-1})
}
{g_{k-1}^{\mathsf T}g_{k-1}}.
$$

This has a useful nonlinear interpretation. The difference

$$
y_{k-1}:=g_{k}-g_{k-1}
$$

is a measured change in the derivative, and for a smooth objective it behaves locally like a Hessian action on the accepted displacement. PR therefore lets the update react to how the gradient **changed**, rather than only to the ratio of gradient magnitudes. On an exact quadratic with exact line search, the relevant orthogonality reduces PR to FR; away from that ideal setting they need not agree.

The commonly used **Polak–Ribière+** variant is

$$
\beta_{k}^{\mathrm{PR+}}
=
\max\left\{0,\beta_{k}^{\mathrm{PR}}\right\}.
$$

If the raw PR coefficient is negative, the old direction is not retained and

$$
d_{k}=-g_{k}.
$$

This is a restart to steepest descent. The truncation is a pragmatic safeguard: the exact conjugacy logic that justified the recurrence has been lost, so a negative mixing coefficient is treated as evidence that the stored direction is no longer useful.

### 11.1 Why restarts are normal

Nonlinear CG inherits its structure from an exact quadratic method, but nonlinear curvature changes from step to step and floating-point arithmetic destroys exact conjugacy even for quadratics.

A restart is therefore not merely an emergency response. It is a way to discard stale directional information when it is no longer helping. Common restart triggers include:

- $\beta_{k}\leq0$ in a PR+ rule;
- a direction that is no longer descent;
- a fixed iteration interval;
- a loss-of-orthogonality or gradient-alignment criterion.

The exact policy varies by implementation. The important idea is that nonlinear conjugacy is a useful memory mechanism, not a permanent invariant.

### 11.2 Metric-aware form

If the derivative covector is $r_{k}$ and the metric gradient is

$$
g_{k}=G^{-1}r_{k},
$$

then Euclidean dot products of gradients are replaced by the corresponding primal-dual pairings. For example,

$$
\lVert g_{k}\rVert_{G}^{2}
=
\langle r_{k},g_{k}\rangle.
$$

The conceptual derivation is unchanged: the metric determines the primal representation in which the search directions live, while derivative changes remain covectors.

## 12. Newton's method minimizes the local quadratic model

Return to the second-order Taylor model

$$
m_{k}(d)
:=
f(x_{k})
+
g_{k}^{\mathsf T}d
+
\frac{1}{2}d^{\mathsf T}H_{k}d.
$$

If $H_{k}$ is positive definite, $m_{k}$ has a unique minimizer. Differentiate with respect to $d$:

$$
\nabla m_{k}(d)
=
g_{k}+H_{k}d.
$$

Setting this to zero gives the **Newton equation**

$$
H_{k}d_{k}=-g_{k}.
$$

The next iterate is then

$$
x_{k+1}=x_{k}+\alpha_{k}d_{k},
$$

where a globalized method may use $\alpha_{k}<1$ until the local quadratic model is reliable enough for full steps.

### 12.1 Why Newton can be fast

Near a sufficiently regular strict minimizer, the Hessian captures the local scaling and coupling among variables. Instead of using one scalar step length to compensate for different curvatures, Newton solves for a direction that incorporates the curvature matrix itself.

Under standard local assumptions, full Newton steps can converge quadratically.

### 12.2 Why Newton is not automatically a descent method

If $H_{k}$ is SPD, then

$$
g_{k}^{\mathsf T}d_{k}
=
-d_{k}^{\mathsf T}H_{k}d_{k}<0
$$

for $d_{k}\ne0$.

If $H_{k}$ is indefinite, the Newton direction need not be descent. Far from a minimizer, this is common in nonconvex problems. Practical Newton methods therefore combine curvature modification, line search, or trust-region logic with the raw Newton equation.

## 13. Newton–CG avoids forming or factoring a dense Hessian

For a large problem, forming $H_{k}$ explicitly may be expensive or impossible, while applying the Hessian to a vector may still be feasible.

The Newton equation

$$
H_{k}d=-g_{k}
$$

can then be solved approximately by an iterative method using Hessian-vector products

$$
v\longmapsto H_{k}v.
$$

When $H_{k}$ is SPD, CG is a natural inner solver. A **Newton–CG** method therefore has two nested levels:

```text
outer nonlinear iteration
        ↓
build current gradient and Hessian action
        ↓
approximately solve H_k d_k = -g_k with CG
        ↓
globalize the resulting step
```

One rarely needs to solve the Newton equation to machine precision far from the solution. An inexact inner tolerance can save substantial work while still producing a useful outer step.

If negative curvature is detected, a truncated-CG trust-region method can stop the inner iteration and use that information rather than insisting on an SPD interpretation. Detailed truncated-CG theory belongs beyond this chapter.

## 14. Secant information comes from accepted steps

Newton uses the true Hessian. Quasi-Newton methods try to infer useful curvature from changes already observed during the iteration.

Suppose an accepted step gives

$$
s_{k}
:=
x_{k+1}-x_{k}
$$

and a gradient change

$$
y_{k}
:=
g_{k+1}-g_{k}.
$$

Taylor expansion of the gradient gives

$$
g(x_{k}+s_{k})
=
g(x_{k})
+
H_{k}s_{k}
+
o(\lVert s_{k}\rVert).
$$

Hence

$$
y_{k}
\approx
H_{k}s_{k}.
$$

This suggests requiring a Hessian approximation $B_{k+1}$ to satisfy the **secant equation**

$$
B_{k+1}s_{k}=y_{k}.
$$

If we work instead with an inverse-Hessian approximation $C_{k+1}$, the equivalent requirement is

$$
C_{k+1}y_{k}=s_{k}.
$$

One pair $(s_{k},y_{k})$ does not determine an entire $n\times n$ matrix. A quasi-Newton update combines the new secant relation with the previous approximation and additional structural requirements such as symmetry and positive definiteness.

## 15. BFGS updates curvature while preserving useful structure

The secant equation asks the next Hessian approximation to map the observed displacement $s_{k}$ to the observed gradient change $y_{k}$:

$$
B_{k+1}s_{k}=y_{k}.
$$

One equation does not determine a matrix, so an update needs more design principles. BFGS keeps three properties that are especially useful for optimization:

- symmetry, because a smooth Hessian is symmetric;
- the new secant relation;
- positive definiteness when the observed curvature is positive, so that solving with the approximation still gives a descent direction.

The formula can be understood as a two-part correction. Start from $B_{k}$ and first remove its old action along $s_{k}$:

$$
B_{k}
-
\frac{B_{k}s_{k}s_{k}^{\mathsf T}B_{k}}
{s_{k}^{\mathsf T}B_{k}s_{k}}.
$$

Multiplying this expression by $s_{k}$ gives zero. Now add a symmetric rank-one term that maps $s_{k}$ to the measured value $y_{k}$:

$$
\frac{y_{k}y_{k}^{\mathsf T}}
{y_{k}^{\mathsf T}s_{k}}.
$$

The result is the **BFGS Hessian update**

$$
B_{k+1}
=
B_{k}
-
\frac{B_{k}s_{k}s_{k}^{\mathsf T}B_{k}}
{s_{k}^{\mathsf T}B_{k}s_{k}}
+
\frac{y_{k}y_{k}^{\mathsf T}}
{y_{k}^{\mathsf T}s_{k}}.
$$

The secant equation is now immediate:

```math
\begin{aligned}
B_{k+1}s_{k}
&=
B_{k}s_{k}
-
B_{k}s_{k}
+
y_{k}
\\
&=y_{k}.
\end{aligned}
```

This does not uniquely characterize BFGS among all conceivable secant updates; a deeper derivation views it as a particular minimal-change member of the Broyden family. For the present purpose, the important heuristic is concrete: **replace the old curvature prediction in the direction just explored by the curvature actually observed, while retaining a symmetric positive-definite model elsewhere as much as possible**.

### 15.1 Why the curvature condition matters

The update requires

$$
y_{k}^{\mathsf T}s_{k}>0.
$$

For the exact Hessian, this quantity is locally

$$
y_{k}^{\mathsf T}s_{k}
\approx
s_{k}^{\mathsf T}H_{k}s_{k},
$$

so positivity is precisely what one expects from positive curvature along the accepted displacement.

There is also a direct algebraic reason. For any vector $z$,

```math
\begin{aligned}
z^{\mathsf T}B_{k+1}z
={}&
z^{\mathsf T}B_{k}z
-
\frac{(z^{\mathsf T}B_{k}s_{k})^{2}}
{s_{k}^{\mathsf T}B_{k}s_{k}}
+
\frac{(z^{\mathsf T}y_{k})^{2}}
{y_{k}^{\mathsf T}s_{k}}.
\end{aligned}
```

If $B_{k}$ is SPD, the first two terms are nonnegative by Cauchy–Schwarz in the $B_{k}$ inner product. The final term is nonnegative, and the strict positivity of $y_{k}^{\mathsf T}s_{k}$ prevents the update from becoming singular in the removed $s_{k}$ direction. Thus BFGS preserves positive definiteness.

That property matters algorithmically. If $B_{k}$ is SPD and the quasi-Newton direction solves

$$
B_{k}d_{k}=-g_{k},
$$

then

$$
g_{k}^{\mathsf T}d_{k}
=
-d_{k}^{\mathsf T}B_{k}d_{k}<0.
$$

So the direction is descent. Wolfe-type line searches are useful here because their curvature condition can imply a positive value of $y_{k}^{\mathsf T}s_{k}$ for a descent step. This is why line search and BFGS are mathematically coupled rather than two unrelated implementation choices.

### 15.2 Inverse BFGS

Instead of solving with $B_{k}$, we can maintain an approximation

$$
C_{k}\approx H_{k}^{-1}
$$

and set

$$
d_{k}=-C_{k}g_{k}.
$$

With

$$
\rho_{k}
:=
\frac{1}{y_{k}^{\mathsf T}s_{k}},
$$

the inverse update is

$$
C_{k+1}
=
\left(I-\rho_{k}s_{k}y_{k}^{\mathsf T}\right)
C_{k}
\left(I-\rho_{k}y_{k}s_{k}^{\mathsf T}\right)
+
\rho_{k}s_{k}s_{k}^{\mathsf T}.
$$

Its corresponding secant relation is

$$
C_{k+1}y_{k}=s_{k}.
$$

Full BFGS therefore carries an $n\times n$ matrix. That becomes unattractive when $n$ is large, which motivates the limited-memory form.

## 16. L-BFGS stores curvature pairs instead of a dense matrix

Limited-memory BFGS avoids storing the full inverse-Hessian approximation. It keeps only the most recent $m$ pairs

$$
(s_{i},y_{i}),
$$

with $m$ small compared with $n$.

The memory requirement becomes

$$
O(mn)
$$

rather than

$$
O(n^{2}).
$$

The key point is that the BFGS inverse updates can be **applied** to a vector using those stored pairs without ever forming the dense matrix they implicitly represent.

### 16.1 The two-loop recursion

Suppose we want

$$
d_{k}=-C_{k}g_{k}.
$$

For each stored pair define

$$
\rho_{i}
=
\frac{1}{y_{i}^{\mathsf T}s_{i}}.
$$

A standard two-loop application is:

```text
q = g_k
for stored pairs from newest to oldest
    alpha_i = rho_i s_i^T q
    q = q - alpha_i y_i

r = C_k^(0) q

for stored pairs from oldest to newest
    beta_i = rho_i y_i^T r
    r = r + s_i (alpha_i - beta_i)

d_k = -r
```

Here $C_{k}^{(0)}$ is a simple initial inverse model for the part of the space not represented by the stored history. A common scalar choice is

$$
C_{k}^{(0)}
=
\gamma_{k}I,
\qquad
\gamma_{k}
=
\frac{s_{k-1}^{\mathsf T}y_{k-1}}
{y_{k-1}^{\mathsf T}y_{k-1}},
$$

when the latest curvature pair is valid.

The two passes are not a separate optimization principle. They are an efficient algebraic way to apply the inverse-BFGS approximation represented by a short history.

### 16.2 What L-BFGS forgets

When a new pair arrives and the memory is full, the oldest pair is discarded. Thus L-BFGS does not reproduce full BFGS exactly after many iterations. It deliberately keeps only recent curvature information.

This is especially attractive when objective/gradient evaluations are expensive but storing a dense matrix would be prohibitive.

## 17. Trust regions globalize the model instead of only the step length

Line search asks:

> Given a direction, how far should we move?

A trust-region method asks instead:

> In what neighborhood do we trust the local model enough to choose a step from it?

At iteration $k$, form a quadratic model

$$
m_{k}(s)
=
f(x_{k})
+
g_{k}^{\mathsf T}s
+
\frac{1}{2}s^{\mathsf T}B_{k}s,
$$

where $B_{k}$ may be the exact Hessian or a symmetric approximation.

Then solve, exactly or approximately,

$$
\min_{\lVert s\rVert\leq\Delta_{k}}m_{k}(s).
$$

The radius $\Delta_{k}$ limits how far we are willing to trust the model.

### 17.1 Actual versus predicted reduction

The model predicts the reduction

$$
\mathrm{pred}_{k}
:=
m_{k}(0)-m_{k}(s_{k}).
$$

The objective delivers the actual reduction

$$
\mathrm{ared}_{k}
:=
f(x_{k})-f(x_{k}+s_{k}).
$$

Compare them through

$$
\rho_{k}
:=
\frac{\mathrm{ared}_{k}}
{\mathrm{pred}_{k}}.
$$

Interpretation:

- $\rho_{k}\approx1$: the local model predicted the objective well;
- $\rho_{k}$ small but positive: the model was useful but optimistic;
- $\rho_{k}\le0$: the proposed step failed to reduce the objective.

A typical policy rejects poor steps and shrinks $\Delta_{k}$, accepts sufficiently good steps, and may enlarge the radius when a boundary step was predicted accurately.

The exact threshold values are implementation choices, not mathematical constants.

### 17.2 The Cauchy step shows that the subproblem need not be solved exactly

Even if $B_{k}$ is complicated or indefinite, we can obtain a useful baseline step by minimizing the quadratic model along the steepest-descent ray

$$
s=-\tau g_{k}
$$

subject to

$$
\lVert s\rVert\leq\Delta_{k}.
$$

This **Cauchy step** guarantees a controlled amount of model reduction under mild conditions. Practical trust-region methods often solve the subproblem more accurately, for example by dogleg or truncated CG, but the central principle is already visible: the subproblem only needs to produce a step that makes the model sufficiently better.

## 18. Line search and trust region solve the same reliability problem differently

Both strategies begin with local information and both need to decide whether that information can be trusted far enough to justify a step.

A line-search method typically looks like

```text
choose direction
      ↓
try step length
      ↓
objective/derivative acceptance test
      ↓
shorten if necessary
```

A trust-region method typically looks like

```text
choose model + radius
      ↓
solve model inside radius
      ↓
compare actual and predicted reduction
      ↓
accept/reject and update radius
```

Neither family is uniformly superior. Line searches are often simpler and pair naturally with BFGS and nonlinear CG. Trust regions are particularly attractive when the Hessian is indefinite, when negative curvature is informative, or when one wants to regulate the reliability of a second-order model directly.

## 19. Stopping criteria measure different symptoms

An iterative optimization method must stop before it reaches exact mathematical stationarity. Several quantities are commonly monitored.

### 19.1 Gradient or stationarity norm

For an unconstrained problem, the most direct first-order measure is

$$
\lVert g_{k}\rVert.
$$

A small gradient is evidence of approximate first-order stationarity.

Under a non-Euclidean metric, the norm should be consistent with the chosen primal geometry rather than silently using coefficient Euclidean norm.

A small gradient does not distinguish a local minimum from a saddle point, and poor scaling can make a raw absolute tolerance misleading.

### 19.2 Step norm

A small update

$$
\lVert x_{k+1}-x_{k}\rVert
$$

can indicate that the iterates have stabilized. It can also indicate that a line search is repeatedly choosing tiny steps or that numerical tolerances are preventing progress.

### 19.3 Objective change

A small relative change such as

$$
\frac{|f(x_{k+1})-f(x_{k})|}
{\max\{1,|f(x_{k})|\}}
$$

can be useful, but a flat objective value alone does not guarantee stationarity.

### 19.4 Iteration and work limits

A maximum number of iterations or expensive evaluations is a safety bound, not an optimality condition.

In applications where one function or gradient evaluation contains a large simulation, iteration counts alone can be misleading. A method that converges in fewer outer iterations may still cost more if each iteration requires many additional solves or rejected trial evaluations.

## 20. Scaling affects the geometry seen by every method

Suppose two variables naturally differ by many orders of magnitude. In raw coordinates, a Euclidean gradient may be dominated by the coordinate scaling rather than by physically meaningful sensitivity.

Changing variables

$$
x=D z
$$

with a nonsingular scaling matrix $D$ changes the coordinate representation of the derivative and Hessian. Equivalently, choosing a nontrivial metric changes which primal vector represents the derivative.

This is why scaling, preconditioning, and metric choice repeatedly appear around optimization algorithms. They are not identical concepts, but all address the fact that raw coordinates may present an unnecessarily distorted geometry to the iteration.

## 21. Work per iteration matters as much as the formula

The update equation alone does not determine computational cost.

A rough comparison is:

| Method | Main information | Typical memory beyond vectors | Extra inner work |
| --- | --- | --- | --- |
| Steepest descent | gradient | low | line search |
| Nonlinear CG | gradient + one direction history | low | line search |
| Full BFGS | gradient + dense inverse/Hessian approximation | $O(n^{2})$ | line search, dense matrix-vector work |
| L-BFGS | gradient + recent curvature pairs | $O(mn)$ | line search, two-loop recursion |
| Newton | gradient + Hessian | potentially high | linear solve |
| Newton–CG | gradient + Hessian-vector action | low to moderate | iterative linear solve |
| Trust region | gradient + Hessian/model action | method dependent | trust-region subproblem |

This table is only structural. In a PDE-constrained problem, for example, a “gradient evaluation” may itself contain state and adjoint solves. Later chapters will make that cost model explicit.

## 22. Scope frontier

This chapter has intentionally stopped after the optimization machinery needed to understand the project's reduced methods.

The companion notes retain a few optional topics that are directly useful for interpreting these methods, including existence, convergence rates, stronger convexity assumptions, and deeper BFGS/Newton–CG remarks. Broader optimization references remain the right place for:

- full global/local convergence proofs and complexity bounds;
- derivative-free methods;
- stochastic and online optimization;
- nonlinear least-squares specializations;
- advanced trust-region subproblem algorithms;
- line-search interpolation strategies in detail;
- large-scale preconditioned Newton theory;
- nonconvex quasi-Newton safeguards in depth;
- accelerated/proximal first-order methods;
- all constrained-optimization machinery.

[06 · Constrained optimization](06-constrained-optimization.md) adds constraints without turning into a general nonlinear-programming survey.

## Where to go next

- **Default continuation:** [06 · Constrained optimization](06-constrained-optimization.md) adds feasible directions, projections, Lagrange multipliers, KKT conditions, complementarity, and active sets.
- **For deeper optimization theory:** [05a · Further notes on unconstrained optimization](05a-further-unconstrained-optimization-notes.md) collects optional existence/convergence background, infinite-dimensional orientation, and selected large-scale methods.
- **For the PDE-control route:** [07 · PDE-constrained optimization](07-pde-constrained-optimization.md) shows why evaluating a reduced objective and derivative can require state and adjoint PDE solves.
- **Into the `nmopt` manual:** two pages are especially useful:
  - [manual 06 · Reduced optimization methods](../concepts/06-reduced-optimization-methods.md), for the project's direction policies, globalization, stopping rules, and work accounting;
  - [manual 04 · Metrics, gradients, and constraints](../concepts/04-metrics-gradients-and-constraints.md), for the discrete metric that turns a derivative covector into the primal gradient used by those algorithms.

## References and further reading

Publication details, stable links, and access notes are collected in the [Background reference catalogue](references.md).

- **Andrea Manzoni, Alfio Quarteroni, and Sandro Salsa, *Optimal Control of Partial Differential Equations*.**
  - Chapter 2, §§2.1–2.2 for minimizers, convexity, and first- and second-order conditions.
  - Chapter 3, especially §§3.1–3.2, for descent directions, nonlinear CG, Newton and quasi-Newton methods, line search, Wolfe conditions, BFGS, and trust regions.
- **Luca Heltai, NMOPT course, lecture *Numerical Optimization Toolbox for Reduced OCPs*.**
  - Teaching progression from gradient descent and line search through nonlinear CG, BFGS/L-BFGS, and trust regions.
- **Jorge Nocedal and Stephen J. Wright, *Numerical Optimization*.**
  - Fuller treatment of line-search, trust-region, conjugate-gradient, Newton, and quasi-Newton methods.
