# 05a · Further notes on unconstrained optimization

**Background navigation:** [Index](README.md) \
Main chapter: [05 · Unconstrained numerical optimization](05-unconstrained-optimization.md) \
Previous: [05 · Unconstrained numerical optimization](05-unconstrained-optimization.md) \
Next: [06 · Constrained optimization](06-constrained-optimization.md)

## Purpose

The main [05 · Unconstrained numerical optimization](05-unconstrained-optimization.md) chapter focuses on the algorithms needed later for reduced PDE-constrained optimization. This companion collects mathematical and algorithmic details that help explain why those methods behave as they do but are not required every time the algorithms are used.

The emphasis remains selective. These notes are not a general optimization course. They develop existence and compactness arguments in finite dimensions, convergence language for iterative sequences, stronger convexity assumptions, and a few large-scale ideas that are especially relevant when one objective or derivative evaluation contains a PDE solve. A short bridge to Banach and Hilbert spaces shows which parts of the finite-dimensional theory survive unchanged and which parts require weak compactness or an explicit Riesz map.

## Before you start

You should be comfortable with:

- the main [05 · Unconstrained numerical optimization](05-unconstrained-optimization.md) material;
- elementary real analysis, including sequences and continuity;
- eigenvalues and SPD matrices from [04 · Numerical linear algebra for PDE and optimization systems](04-numerical-linear-algebra.md).

## What you will be able to do

After these notes, you should be able to:

- explain why continuity alone does not guarantee that a minimizer exists on all of $\mathbb R^{n}$;
- use compact sublevel sets or coercivity to obtain existence in finite dimensions;
- interpret a minimizing sequence and distinguish existence from algorithmic convergence;
- distinguish linear, superlinear, and quadratic convergence of iterates;
- relate strong convexity and Lipschitz gradients to Hessian bounds and conditioning;
- recognize how stationarity, convexity, existence, and second-order conditions extend to Banach and Hilbert spaces;
- understand why inexact and truncated Newton ideas are natural for large PDE-control problems;
- recognize what the Rosenbrock function is actually illustrating when it appears in optimization courses.

## Roadmap

We first separate **existence of a minimizer** from **convergence of an algorithm** and refine the convexity picture used in the main chapter. A short infinite-dimensional bridge then shows how the same ideas look in Banach and Hilbert spaces, especially the roles of weak compactness, dual derivatives, and Riesz maps. We next introduce standard convergence-rate language and selected large-scale algorithmic variants. A final short example explains why Rosenbrock is a useful nonlinear stress test without making it part of the prerequisite thread.

## 1. Existence is a different question from stationarity

The equation

$$
\nabla f(x)=0
$$

is a local necessary condition. It says nothing by itself about whether a minimizer exists.

For example,

$$
f(x)=e^{x}
$$

is smooth and bounded below by zero, but it never attains that lower bound on $\mathbb R$. The sequence

$$
x_{k}=-k
$$

satisfies

$$
f(x_{k})\to0,
$$

but there is no $x_{\ast}\in\mathbb R$ with $f(x_{\ast})=0$.

This motivates the notion of a **minimizing sequence**. Let

$$
m:=\inf_{x\in K} f(x).
$$

A sequence $\{x_{k}\}\subset K$ is minimizing if

$$
f(x_{k})\to m.
$$

Every optimization problem with finite infimum admits minimizing sequences by the definition of infimum. The real question is whether one can extract a subsequence converging to a feasible point where the infimum is attained.

### 1.1 Compactness plus continuity gives existence

Suppose $K\subset\mathbb R^{n}$ is nonempty and compact and $f:K\to\mathbb R$ is continuous. Take a minimizing sequence $\{x_{k}\}\subset K$.

Compactness gives a convergent subsequence

$$
x_{k_{j}}\to x_{\ast}\in K.
$$

Continuity then gives

$$
f(x_{\ast})
=
\lim_{j\to\infty}f(x_{k_{j}})
=
m.
$$

So the minimum is attained. This is the finite-dimensional Weierstrass argument.

The logic matters later in infinite-dimensional optimization because the corresponding compactness step becomes much more delicate. Here, however, ordinary compactness in $\mathbb R^{n}$ is enough.

### 1.2 Coercivity replaces a globally bounded search set

The domain $\mathbb R^{n}$ is not compact. A common replacement is **coercivity at infinity**:

$$
f(x)\to+\infty
\qquad
\text{as }\lVert x\rVert\to\infty.
$$

Fix any point $x_{0}$ and consider the sublevel set

$$
L_{0}
:=
\{x:f(x)\leq f(x_{0})\}.
$$

Coercivity makes this set bounded. If $f$ is continuous, it is also closed, hence compact in finite dimensions. Since no global minimizer can have value larger than the trial value $f(x_{0})$, the search can be restricted to $L_{0}$ and the compactness argument applies.

This is why **compact level sets** appear so often in convergence assumptions. They prevent an improving sequence from escaping to infinity.

## 2. Convexity, strong convexity, and smoothness constants

The main chapter used convexity to turn stationarity into global optimality. Two stronger quantitative notions are also useful for understanding rates and conditioning.

### 2.1 Strong convexity gives uniform positive curvature

A differentiable function is $m$-strongly convex, with $m>0$, if

$$
f(y)
\geq
f(x)
+
\nabla f(x)^{\mathsf T}(y-x)
+
\frac{m}{2}\lVert y-x\rVert^{2}
$$

for all $x,y$.

If $f$ is twice continuously differentiable, a sufficient and, on a convex domain, equivalent condition is

$$
\nabla^{2}f(x)\succeq mI
$$

for every $x$ in the domain.

Strong convexity implies a unique global minimizer and rules out arbitrarily flat directions.

### 2.2 A Lipschitz gradient bounds curvature from above

Suppose

$$
\lVert\nabla f(x)-\nabla f(y)\rVert
\leq
L\lVert x-y\rVert.
$$

Then the gradient is $L$-Lipschitz. For a twice differentiable function this corresponds to the upper curvature bound

$$
\nabla^{2}f(x)\preceq LI
$$

when the Hessian is symmetric.

The pair

$$
mI\preceq\nabla^{2}f(x)\preceq LI
$$

is the nonlinear analogue of bounding the eigenvalues of an SPD matrix. The ratio

$$
\frac{L}{m}
$$

plays a conditioning role similar to the condition number of the quadratic Hessian.

### 2.3 The descent lemma gives a quantitative local model

A Lipschitz gradient implies

$$
f(x+s)
\leq
f(x)
+
\nabla f(x)^{\mathsf T}s
+
\frac{L}{2}\lVert s\rVert^{2}.
$$

This is often called the **descent lemma**. It turns the qualitative statement

> sufficiently small descent steps decrease the objective

into a quantitative bound. For steepest descent, $s=-\alpha g$, it gives

$$
f(x-\alpha g)
\leq
f(x)
-
\alpha\left(1-\frac{L\alpha}{2}\right)\lVert g\rVert^{2}.
$$

Thus any

$$
0<\alpha<\frac{2}{L}
$$

produces decrease. This is one analytical explanation for why step-size restrictions depend on curvature.

## 3. What survives in Banach and Hilbert spaces

The finite-dimensional theory above is not tied to coordinates as strongly as its notation may suggest. Let $X$ be a Banach space and let

$$
F:X\to\mathbb R
$$

be a nonlinear functional. If $u_{\ast}$ is an unconstrained local minimizer and $F$ is Fréchet differentiable, then

$$
F'(u_{\ast})=0
\qquad\text{in }X^{\ast}.
$$

This is exactly the same first-order statement as in $\mathbb R^{n}$. What changes is the type of the derivative: $F'(u)$ is naturally a covector in the dual space $X^{\ast}$.

If $X$ is a Hilbert space, a chosen inner product gives a Riesz map

$$
G:X\to X^{\ast},
$$

and one may define the gradient by

$$
G\nabla F(u)=F'(u).
$$

Thus the distinction between derivative and gradient from [02 · Duality, derivatives, and adjoints](02-duality-derivatives-and-adjoints.md) becomes even more important in infinite dimensions. The equation $F'(u)=0$ is metric-independent; the vector called $\nabla F(u)$ depends on the Hilbert-space inner product used to identify $X$ with $X^{\ast}$.

### 3.1 Convexity still turns first-order conditions into global conditions

For a differentiable convex functional on a convex set,

$$
F(v)
\geq
F(u)+\langle F'(u),v-u\rangle_{X^{\ast},X}.
$$

Consequently, in the unconstrained case,

$$
F'(u_{\ast})=0
$$

still implies that $u_{\ast}$ is a global minimizer. Strict convexity still gives uniqueness. These facts do not depend on finite dimensionality.

The constrained version is equally revealing. If $K\subset X$ is closed and convex and $u_{\ast}\in K$ minimizes $F$, then the natural first-order condition is

$$
\langle F'(u_{\ast}),v-u_{\ast}\rangle_{X^{\ast},X}
\geq0
\qquad
\forall v\in K.
$$

This is the same variational inequality developed in [06 · Constrained optimization](06-constrained-optimization.md). A Hilbert structure is useful when we want to rewrite it using gradients, orthogonal projections, or normal cones represented in the primal space, but the derivative-level statement itself only needs the dual pairing.

### 3.2 Existence is where infinite dimensions change the argument

In $\mathbb R^{n}$, a closed and bounded set is compact. That fact powered the finite-dimensional existence argument in Section 1. In an infinite-dimensional Banach space, closed and bounded sets are generally not strongly compact, so the same proof breaks at the subsequence-extraction step.

A standard replacement is to work with **weak convergence**. In a reflexive Banach space, bounded sequences admit weakly convergent subsequences. An existence argument can therefore proceed if:

- a minimizing sequence is bounded, for example by coercivity;
- the feasible set is weakly sequentially closed;
- the functional is weakly sequentially lower semicontinuous.

Convex closed sets and convex lower-semicontinuous functionals fit this framework particularly well. This is the infinite-dimensional version of the compact-level-set argument rather than an unrelated theorem.

### 3.3 Second-order sufficiency has the same shape

Suppose $F$ is twice Fréchet differentiable, $F'(u_{\ast})=0$, and there is a constant $\alpha>0$ such that

$$
\langle F''(u_{\ast})w,w\rangle_{X^{\ast},X}
\geq
\alpha\lVert w\rVert_{X}^{2}
\qquad
\forall w\in X.
$$

A second-order Taylor expansion then gives a strict local minimum, exactly as positive-definite Hessian curvature does in finite dimensions. For constrained problems the coercivity requirement is imposed only on appropriate feasible or critical directions; [06 · Constrained optimization](06-constrained-optimization.md) and its companion develop the finite-dimensional version of that idea.

### 3.4 Why nonlinear PDE control makes this bridge useful

For a PDE-constrained problem, one often writes a reduced objective

$$
j(u)=J(S(u),u),
$$

where $S$ is the control-to-state map. Even when the original cost $J$ is convex in its arguments, a nonlinear map $S$ can make $j$ nonconvex. Then the first-order condition remains necessary, but it is no longer generally sufficient for global optimality.

This observation is useful even when the numerical algorithm is ultimately applied after discretization. It separates two issues that can otherwise be conflated: the function-space calculus that justifies derivatives and adjoints, and the finite-dimensional algorithms used to solve the resulting discrete optimization problem.

## 4. Convergence of iterates and convergence rates

Suppose an algorithm generates $x_{k}\to x_{\ast}$. Define the error

$$
e_{k}:=\lVert x_{k}-x_{\ast}\rVert.
$$

The terminology used in optimization describes how rapidly $e_{k}$ shrinks once the iterates are sufficiently close to the limit.

### 4.1 Linear convergence

The convergence is **linear** if, eventually,

$$
e_{k+1}\leq q e_{k}
$$

for some fixed

$$
0<q<1.
$$

Each iteration removes a roughly fixed fraction of the remaining error.

### 4.2 Superlinear convergence

The convergence is **superlinear** if

$$
\frac{e_{k+1}}{e_{k}}\to0.
$$

The fraction of error removed improves as the solution is approached. BFGS is famous for superlinear local behavior under suitable assumptions.

### 4.3 Quadratic convergence

The convergence is **quadratic** if, eventually,

$$
e_{k+1}\leq C e_{k}^{2}
$$

for some $C>0$.

Once $e_{k}$ is small, squaring it can reduce the error extremely rapidly. Newton's method has this local behavior under the usual nonsingularity and smoothness assumptions.

These are statements about asymptotic local behavior, not about total wall-clock cost. A method with a better local rate may have a much more expensive iteration.

## 5. Why line-search assumptions mention bounded level sets

Suppose a descent algorithm accepts only steps satisfying

$$
f(x_{k+1})\leq f(x_{k}).
$$

Then every iterate remains in the initial sublevel set

$$
L(x_{0})
:=
\{x:f(x)\leq f(x_{0})\}.
$$

If that set is compact, the sequence cannot diverge to infinity and must have accumulation points. A convergence theorem then tries to show that every such accumulation point is stationary, often using assumptions on gradient continuity and the line search.

This explains why compact level-set assumptions occur even though a line-search formula itself contains no compactness. The step rule gives monotonicity; compactness supplies somewhere for the iterates to accumulate.

A full proof of global convergence for Armijo/Wolfe methods requires additional arguments and is intentionally left to the optimization references.

## 6. Inexact Newton methods are especially natural when a Hessian action is expensive

The exact Newton equation is

$$
H_{k}d_{k}=-g_{k}.
$$

For a large problem, solving this system to machine precision at every outer iteration is usually wasteful. Far from the solution, the quadratic model itself is only an approximation, so an extremely accurate inner linear solve does not necessarily produce a meaningfully better outer step.

An **inexact Newton** method instead computes $d_{k}$ so that

$$
\lVert H_{k}d_{k}+g_{k}\rVert
\leq
\eta_{k}\lVert g_{k}\rVert,
$$

where $\eta_{k}$ is a forcing tolerance. Early iterations can use a loose tolerance; the inner solve can become more accurate as the nonlinear residual decreases.

For PDE-constrained optimization this is particularly important because a Hessian-vector product may itself require tangent and incremental-adjoint solves. The useful unit of work is not merely an outer iteration but the number and accuracy of those PDE solves.

### 6.1 Truncated Newton and negative curvature

If CG is used inside Newton and the Hessian ceases to look positive definite, one may encounter a direction $v$ with

$$
v^{\mathsf T}H_{k}v\leq0.
$$

For minimization, this is not a reason to force CG to continue as though the problem were SPD. In a trust-region framework, the inner iteration can stop and use the detected negative-curvature direction to move toward the trust-region boundary.

This is the basic idea behind **truncated Newton** or **truncated-CG trust-region** methods. The main chapter already contains the conceptual ingredients: Hessian actions, CG, negative curvature, and a trust radius. Detailed convergence theory and sophisticated preconditioning remain beyond the background track.

## 7. What the BFGS update is trying to preserve

The main chapter derived the BFGS rank-two correction by first removing the old predicted curvature along $s_{k}$ and then inserting the measured curvature $y_{k}$.

There is a deeper viewpoint. The secant equation

$$
B_{k+1}s_{k}=y_{k}
$$

defines many symmetric matrices. Quasi-Newton theory asks which one should be chosen while changing the previous model as little as possible according to a suitable matrix norm. Different choices produce members of the **Broyden family**; BFGS is a particularly successful member because it combines the secant relation with symmetry and positive-definiteness preservation.

For the project background, the full variational derivation of the Broyden family would add more machinery than insight. The operational lesson is more important:

```text
accepted displacement s_k
        +
observed derivative change y_k
        ↓
one new curvature sample
        ↓
update the inverse/Hessian model without rebuilding it from scratch
```

L-BFGS keeps the same logic but remembers only a short list of those samples.

## 8. Methods especially useful for large reduced PDE objectives

Large reduced PDE objectives favor methods that can be expressed through vectors, derivative evaluations, and operator actions rather than dense matrix construction. That makes a few families particularly natural.

**L-BFGS** is attractive when gradients are available but storing or assembling a reduced Hessian is unrealistic. **Newton–CG and inexact Newton** become attractive when Hessian-vector products can be computed by tangent and incremental-adjoint solves. **Truncated-CG trust-region methods** add a useful response to negative curvature and allow the inner solve to stop before an exact Newton step has been computed. **Nonlinear CG** remains relevant because its storage cost is minimal and each accepted iteration can be built from a new gradient plus short history.

These methods differ substantially in local convergence and robustness, but they share the structural property that matters most for PDE optimization: none requires a dense reduced Hessian. More specialized material such as the full Broyden family, elaborate line-search interpolation, dogleg constructions, or exact trust-region eigensolvers is better treated in a general optimization reference unless a later application actually needs it.

## 9. Rosenbrock is a stress test, not a model problem for PDE control

The two-variable Rosenbrock function is

$$
f(x_{1},x_{2})
=
100(x_{2}-x_{1}^{2})^{2}+(1-x_{1})^{2}.
$$

Its global minimizer is

$$
(1,1),
$$

where the objective vanishes. The function is useful in optimization courses because the minimizer lies in a narrow **curved** valley. This stresses several features at once:

- steepest descent can make slow zig-zagging progress;
- a locally good direction can still need globalization;
- changing curvature challenges methods that reuse old directional information;
- quasi-Newton methods can learn the local valley geometry without an exact Hessian.

It is not theoretically fundamental and it is not a surrogate for a PDE-constrained objective. Its value is diagnostic: it makes nonlinear behavior visible in two dimensions. That is why the NMOPT optimization lecture uses it for finite-dimensional demonstrations, while the main prerequisite chapter relies instead on derivations that transfer directly to large reduced problems.

## 10. What we still leave to broader optimization texts

Even these companion notes stop before:

- full proofs of Zoutendijk-type line-search convergence results;
- detailed nonlinear-CG convergence theory;
- the full Broyden class and matrix minimal-change derivations;
- trust-region subproblem theory, dogleg geometry, and Lanczos methods;
- self-concordance and interior-point complexity;
- stochastic, nonsmooth, derivative-free, and proximal optimization;
- complexity bounds for first- and second-order methods.

Those topics can be valuable, but they are not needed to understand the optimization machinery used by the project.

## Where to go next

- **Return to the main route:** continue with [06 · Constrained optimization](06-constrained-optimization.md) for constraints, multipliers, projection, and KKT conditions.
- **For infinite-dimensional optimization:** keep the Banach/Hilbert-space material here as orientation, then return to it in [07 · PDE-constrained optimization](07-pde-constrained-optimization.md) and [07a · Further notes on PDE-constrained optimization](07a-further-pde-constrained-optimization-notes.md), where the derivative and Hessian structures acquire PDE meaning.
- **Into the `nmopt` manual:** [manual 06 · Reduced optimization methods](../concepts/06-reduced-optimization-methods.md) is the project-specific continuation for nonlinear CG, BFGS/L-BFGS, Newton, line search, and trust regions.

## References and further reading

Publication details, stable links, and access notes are collected in the [Background reference catalogue](references.md).

- **Andrea Manzoni, Alfio Quarteroni, and Sandro Salsa, *Optimal Control of Partial Differential Equations*.**
  - Chapter 2, §§2.1–2.2 for finite-dimensional existence, convexity, admissible directions, and second-order conditions.
  - Chapter 3, §§3.1–3.2 for descent, convergence terminology, BFGS, line searches, and trust regions.
  - Chapter 9, §§9.2–9.3 for optimization in Banach/Hilbert spaces and nonlinear-PDE control.
  - Chapter 6, §6.3 for the cost model of reduced unconstrained PDE optimal control.
- **Luca Heltai, NMOPT course, lecture *Numerical Optimization Toolbox for Reduced OCPs*.**
  - Gradient descent, nonlinear CG, BFGS/L-BFGS, trust regions, and Rosenbrock demonstrations.
- **Jorge Nocedal and Stephen J. Wright, *Numerical Optimization*.**
  - Recommended for full convergence analyses and large-scale variants only sketched here.
