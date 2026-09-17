# Reduced optimization methods

## From a reduced derivative to an optimization iteration

The previous chapter ended with a reduced problem

$$
\min_{u\in U_h} j_h(u)
$$

and a procedure for evaluating

$$
j_h(u)
\qquad\text{and}\qquad
j_h'(u)\in U_h^{\ast}.
$$

That is enough to ask whether a proposed direction is a descent direction, but it does
not yet tell us which direction to choose, how far to move, when to stop, or how to
react when a local model is unreliable.

Those are optimization-algorithm questions. The current reduced solver layer separates
them into two broad families:

```text
line-search methods
    choose a direction d
    then choose a step length α

trust-region methods
    build a local quadratic model
    and choose a step inside ||s||_G ≤ Δ
```

Within the line-search family, the direction itself is another replaceable policy. The
current implementation includes steepest descent, several conjugate-gradient updates,
full and limited-memory BFGS, and Newton directions when an explicit reduced Hessian
action is available.

This chapter develops those methods from the reduced mathematical problem rather than
from the solver class hierarchy. The recurring objects are:

```text
u_k                 current control
r_k = j_h'(u_k)     reduced derivative covector
G                   chosen metric / Riesz map
g_k = G^-1 r_k      metric gradient
d_k                 primal search direction
α_k                 line-search step length
s_k                 actual accepted update
H_k                 exact or approximate reduced Hessian
```

The metric conventions were developed in [Metrics, gradients, and
constraints](metrics-gradients-and-constraints.md), while the state/adjoint work hidden
inside each reduced evaluation was developed in [Reduced state–adjoint
formulation](reduced-state-adjoint-formulation.md).

## 1. One generic line-search iteration

Suppose the current reduced evaluation contains

$$
j_h(u_k)
\qquad\text{and}\qquad
r_k:=j_h'(u_k).
$$

A line-search iteration can be written schematically as

```text
current reduced evaluation
        │
        ▼
choose direction d_k
        │
        ▼
check r_k[d_k] < 0
        │
        ▼
choose trial step α
        │
        ▼
build trial control
u_trial = u_k + α d_k
(or its projected version)
        │
        ▼
state solve + objective
        │
        ▼
acceptance test
   ├── reject → change α and try again
   └── accept
          │
          ▼
      augment derivative
          │
          ▼
      next reduced evaluation
```

The key quantity before globalization is the directional derivative

$$
r_k[d_k]
=
\langle r_k,d_k\rangle.
$$

A descent direction satisfies

$$
r_k[d_k]<0.
$$

For sufficiently small positive $\alpha$,

$$
j_h(u_k+\alpha d_k)
=
j_h(u_k)
+
\alpha r_k[d_k]
+
o(\alpha),
$$

so a negative directional derivative predicts local decrease.

The line search exists because that first-order prediction is only local. A direction
can be a perfectly valid descent direction while a unit step is far too large.

### 1.1 Direction choice and globalization solve different problems

It is useful to keep two questions separate. The **direction policy** asks:

> Which primal vector should represent useful local movement?

The **line-search policy** asks:

> How much of that movement should be accepted?

For example, steepest descent with an Armijo search and Newton with an Armijo search
share the same globalization rule but construct $d_k$ very differently.

Conversely, the same steepest-descent direction can be paired with a fixed step, Armijo
search, weak Wolfe search, or strong Wolfe search.

This is the mathematical reason the current `ReducedSearchSolverT` is parameterized by a
direction policy and a line-search policy rather than encoding each combination as a
separate algorithm.

## 2. Every direction starts from the derivative–metric distinction

The reduced derivative is a covector

$$
r_k\in U_h^{\ast}.
$$

The metric supplies the primal gradient

$$
g_k
:=
G^{-1}r_k.
$$

Its metric norm is

$$
\lVert g_k\rVert_G
=
\sqrt{
\langle Gg_k,g_k\rangle
}.
$$

Since

$$
Gg_k=r_k,
$$

the same norm satisfies

$$
\lVert g_k\rVert_G^2
=
\langle r_k,g_k\rangle.
$$

The helper `make_metric_gradient()` performs exactly these steps:

```text
r_k
 │ inverse_apply
 ▼
g_k
 │ apply
 ▼
G g_k
 │ pair with g_k
 ▼
||g_k||_G
```

This gradient norm is also the basic unconstrained stationarity measure used by the
line-search and trust-region solvers.

## 3. Steepest descent is the reference direction

The simplest direction is

$$
d_k:=-g_k.
$$

Then

$$
r_k[d_k]
=
-\lVert g_k\rVert_G^2
\leq0.
$$

Away from a stationary point the inequality is strict.

This method uses no history and no second-order model. Each iteration asks only for the
current reduced derivative and one inverse metric application.

That makes steepest descent an important reference even when it is not the final
algorithm one intends to use.

Several more elaborate direction policies in the current implementation explicitly fall
back to $-g_k$ when their history becomes unreliable or fails a curvature/descent test.

### 3.1 Why a more elaborate direction can help

Steepest descent chooses the direction of fastest local decrease under the selected
metric, but it does not account for how rapidly the objective bends in different
directions.

For a quadratic reduced objective

$$
j(u)
=
\frac12\langle Hu,u\rangle
-
\langle b,u\rangle,
$$

the level sets can be highly elongated if the Hessian has very different curvatures in
different modes. A steepest path then tends to zig-zag across the narrow direction.

Conjugate-gradient, quasi-Newton, and Newton methods use progressively richer
information about previous derivatives or the Hessian to reduce that behavior.

## 4. Conjugate-gradient directions try not to undo previous progress

The motivation for conjugate gradient is easiest to see for a strictly convex quadratic
reduced objective

$$
j(u)
:=
\frac12\langle Hu,u\rangle
-
\langle b,u\rangle,
$$

with a symmetric positive-definite Hessian

$$
H:
U_h
\longrightarrow
U_h^{\ast}.
$$

Steepest descent repeatedly follows the current metric gradient. On an elongated
quadratic, successive gradients tend to point across the narrow valley in alternating
directions, so one iteration can partially undo progress made by the previous one.

Classical conjugate gradient changes the question. Instead of requiring successive
directions to be geometrically orthogonal in the metric, it constructs directions that
are **conjugate with respect to the Hessian**:

$$
\langle H d_i,d_j\rangle
=
0
\qquad
i\neq j.
$$

This removes the quadratic cross-coupling between directions.

To see why that matters, suppose an exact line search has already minimized the
quadratic along $d_i$, and a later update moves along an $H$-conjugate direction $d_j$.
In the quadratic term, the interaction between those two components contains

$$
\langle H d_i,d_j\rangle.
$$

Because this pairing is zero, movement along $d_j$ does not reintroduce second-order
coupling with the already-treated $d_i$ direction. Informally, each new conjugate
direction can remove a new component of the error without zig-zagging back through the
curvature direction handled previously.

### 4.1 The quadratic method builds an expanding Krylov space

The metric turns the Hessian into a primal-to-primal operator

$$
A
:=
G^{-1}H:
U_h
\longrightarrow
U_h.
$$

Starting from the initial metric gradient $g_0$, classical preconditioned CG explores
the Krylov spaces

$$
\mathcal K_m(A,g_0)
:=
\mathrm{span}
\left\{
g_0,
Ag_0,
A^2g_0,
\ldots,
A^{m-1}g_0
\right\}.
$$

After $m$ iterations, the quadratic minimizer is sought in an affine space of the form

$$
u_0+\mathcal K_m(A,g_0).
$$

The previous direction is therefore not retained merely because "momentum" is useful. It
participates in a recurrence that represents an expanding curvature-informed search
space without storing the whole Krylov basis.

In exact arithmetic, an SPD quadratic in an $n$-dimensional control space can be solved
by classical CG in at most $n$ steps. The practical attraction is that this uses Hessian
actions and short recurrences rather than a dense Hessian factorization.

### 4.2 Nonlinear CG borrows that memory without having a fixed quadratic

For a nonlinear reduced objective there is no single constant Hessian $H$ whose
conjugacy relations remain exact from one iterate to the next.

Nonlinear CG nevertheless keeps the same structural idea:

$$
d_k
=
-g_k+\beta_k d_{k-1}.
$$

The previous direction contains information about the orientation of the local valley
encountered by the last step. The coefficient $\beta_k$ tries to retain the useful part
of that information while adapting it to the new derivative.

This should be understood as a nonlinear analogue of the quadratic Krylov recurrence,
not as exact conjugacy for a changing Hessian. The implementation supports three related
updates:

```text
Polak–Ribière+
Fletcher–Reeves
classical quadratic CG specialization
```

All are expressed with the declared metric and explicit primal/covector pairings.

### 4.3 Fletcher–Reeves in the declared metric

Because

$$
\langle r_k,g_k\rangle
=
\lVert g_k\rVert_G^2,
$$

the Fletcher–Reeves coefficient can be written

$$
\beta_k^{\mathrm{FR}}
:=
\frac{
\langle r_k,g_k\rangle
}{
\langle r_{k-1},g_{k-1}\rangle
}.
$$

Then

$$
d_k
=
-g_k+\beta_k^{\mathrm{FR}}d_{k-1}.
$$

This is the usual Fletcher–Reeves formula with Euclidean gradient dot products replaced
by the primal/dual pairing induced by the selected metric gradient.

### 4.4 Polak–Ribière+

The implemented Polak–Ribière numerator is

```math
\langle r_k,g_k\rangle
-
\langle r_{k-1},g_k\rangle
=
\langle r_k-r_{k-1},g_k\rangle.
```

Thus

$$
\beta_k^{\mathrm{PR}}
:=
\frac{
\langle r_k-r_{k-1},g_k\rangle
}{
\langle r_{k-1},g_{k-1}\rangle
}.
$$

The `+` variant used by the code does not continue with a nonpositive or nonfinite
coefficient. It restarts instead. Conceptually that means replacing

$$
d_k=-g_k+\beta_k d_{k-1}
$$

by

$$
d_k=-g_k
$$

when the conjugacy history is no longer judged useful.

### 4.5 Restarts are part of the method, not exceptional failure

The current nonlinear-CG policy also restarts after a configurable number of directions.
A zero configuration selects an interval equal to the total control dimension.

History is also discarded when the curvature denominator becomes too small or invalid,
or when the proposed direction ceases to be descent.

This is important because the clean conjugacy properties of linear CG hold for a
quadratic objective with exact arithmetic and suitable line searches. A nonlinear
reduced objective does not preserve those ideal conditions automatically.

The restart turns the method back into steepest descent for one iteration and then
allows conjugate history to build again.

### 4.6 The classical quadratic specialization

The third policy is intended for the quadratic/exact-line-search route.

It uses the same norm-ratio coefficient as the Fletcher–Reeves expression, but treats
curvature/descent violations more strictly rather than silently converting them into a
nonlinear restart path.

It is therefore best understood in the context of the exact quadratic line search
described later in this chapter, not as another generic nonlinear-CG recommendation.

## 5. Quasi-Newton methods infer Hessian information from derivative changes

Newton's method uses the reduced Hessian directly. For a large PDE control space,
forming that Hessian explicitly would be expensive even when a matrix-free Hessian
action exists.

If the control dimension is $n_u$, assembling all columns of the reduced Hessian by
probing basis directions would require up to $n_u$ Hessian actions. A single reduced
Hessian action may itself contain tangent-state, incremental-adjoint, or other PDE work.

Quasi-Newton methods instead try to capture useful local curvature information from work
that the optimization iteration has already performed. Suppose two accepted iterates
give the primal displacement

$$
s_k
:=
u_{k+1}-u_k
\in U_h
$$

and the derivative change

$$
\Delta r_k
:=
r_{k+1}-r_k
\in U_h^{\ast}.
$$

The derivative itself has the first-order expansion

$$
j_h'(u_k+s)
=
j_h'(u_k)
+
H_k s
+
o(\lVert s\rVert),
$$

where

$$
H_k
:=
j_h''(u_k).
$$

Therefore an accepted step supplies the observed relation

$$
\Delta r_k
\approx
H_k s_k.
$$

This is one sample of how the Hessian acts, obtained without applying the Hessian to a
basis of the whole control space.

### 5.1 The secant condition turns that observation into an approximation

A quasi-Newton Hessian approximation $B_{k+1}$ can be required to reproduce the observed
change exactly:

$$
\Delta r_k
=
B_{k+1}s_k.
$$

This is the **secant condition**. If we work instead with an inverse-Hessian
approximation

$$
C_{k+1}:
U_h^{\ast}
\longrightarrow
U_h,
$$

the equivalent condition is

$$
C_{k+1}\Delta r_k
=
s_k.
$$

The point is not that one secant pair identifies the entire Hessian. It does not.
Rather, each accepted pair constrains the approximation along a direction the algorithm
has actually explored.

The remaining degrees of freedom are inherited from the previous approximation.

### 5.2 What BFGS adds to the secant idea

Many updates can satisfy one secant equation. BFGS—the Broyden–Fletcher–Goldfarb–Shanno
update—selects a particular rank-two correction of the previous approximation. Its
important properties are:

- it satisfies the new secant condition;
- it preserves symmetry;
- if the previous approximation is positive definite and
  $\langle\Delta r_k,s_k\rangle>0$, the updated approximation remains positive
  definite.

The last property matters because a positive-definite inverse approximation tends to
produce a descent direction

$$
d_k
\approx
-C_k r_k.
$$

So BFGS is not merely "remember the previous iterate". It progressively reshapes an
approximate inverse Hessian so that observed derivative changes are reproduced while
retaining a well-behaved curvature model.

For the `nmopt` solver boundary, those structural properties are more important than the
closed-form rank-two matrix formula. The implementation acts through stored secant
information and inverse-approximation actions rather than assembling a dense approximate
Hessian; a general optimization background chapter can derive the classical formula in
detail.

### 5.3 Why the secant pair has mixed primal/dual types

The displacement

$$
s_k\in U_h
$$

is primal. The derivative change

$$
\Delta r_k\in U_h^{\ast}
$$

is a covector. Therefore

$$
\langle \Delta r_k,s_k\rangle
$$

is a natural dual pairing.

This is exactly how the current BFGS implementations store their history: a primal
displacement and a covector derivative difference, not two untyped coefficient vectors.

The type distinction developed in Part I is doing real algorithmic work here.

## 6. Limited-memory BFGS stores only recent secant information

A full inverse-Hessian approximation is dense even when the underlying PDE
discretization is sparse. L-BFGS avoids storing such a matrix explicitly. Instead it
retains a small number of recent secant pairs

$$
(s_i,\Delta r_i)
$$

and applies the implicit inverse-Hessian approximation through a two-loop recursion.

So neither the exact reduced Hessian nor a dense approximate inverse Hessian needs to be
assembled. The algorithm reconstructs the action of the approximation on the current
derivative from the stored secant information.

The current parameters include:

```text
memory size
curvature tolerance
initial inverse-Hessian scaling
```

The default memory size is finite, so when a new pair arrives after the history is full,
the oldest pair is discarded.

### 6.1 Curvature determines whether a secant pair is retained

For a candidate pair,

$$
s_k=u_k-u_{k-1},
\qquad
\Delta r_k=r_k-r_{k-1},
$$

the relevant curvature is

$$
\rho_k^{-1}
:=
\langle \Delta r_k,s_k\rangle.
$$

If this pairing is not sufficiently positive and finite, the current L-BFGS
implementation clears its history and falls back to steepest descent. That behavior
protects the inverse approximation from incorporating a secant pair that is incompatible
with the positive-curvature model BFGS relies on.

The algorithm records such an event as a curvature reset rather than silently continuing
with a corrupted history.

### 6.2 The metric inverse supplies the baseline inverse approximation

The two-loop recursion needs an initial inverse action after the backward pass through
the history. The default choice is the declared metric inverse:

$$
G^{-1}.
$$

So before any BFGS history has been learned, the method naturally reduces to a
metric-gradient step. This is a useful relationship:

```text
no secant history
    ↓
metric inverse
    ↓
steepest descent

accepted secant history
    ↓
BFGS corrections around that baseline
    ↓
quasi-Newton direction
```

The metric is therefore not discarded when BFGS is enabled. It remains the baseline
primal–dual identification on which the quasi-Newton correction is built.

### 6.3 Optional scalar secant scaling

The implementation can also scale the baseline inverse using the most recent secant
information. If

$$
\widetilde g_k
=
G^{-1}\Delta r_k,
$$

then

$$
\langle \Delta r_k,\widetilde g_k\rangle
=
\langle \Delta r_k,G^{-1}\Delta r_k\rangle.
$$

The scalar scaling used by the current implementation is

$$
\gamma_k
:=
\frac{
\langle \Delta r_k,s_k\rangle
}{
\langle \Delta r_k,G^{-1}\Delta r_k\rangle
}.
$$

The baseline inverse action becomes approximately

$$
\gamma_k G^{-1}
$$

before the two-loop corrections are applied. Again, if the denominator is not suitably
positive, the history is reset rather than forcing the scaling through an invalid
curvature configuration.

## 7. Full BFGS keeps all accepted secant pairs

The project also includes a full-memory BFGS direction policy. "Full" here does **not**
mean that it constructs a dense inverse-Hessian matrix. It retains every accepted secant
pair and applies the corresponding inverse-BFGS history with the same kind of two-loop
recursion.

This has two consequences. First, the storage grows with the number of accepted
iterations rather than being capped by a fixed memory size.

Second, the implementation remains backend-parametric: it only needs primal/covector
updates, pairings, and the metric inverse. It does not require the backend to expose
dense matrix storage.

### 7.1 The same reset logic protects descent

Full BFGS checks

$$
\langle \Delta r_k,s_k\rangle>0
$$

up to its configured curvature tolerance. A bad pair clears the accumulated history and
causes the current direction to fall back to steepest descent.

Even with valid secant pairs, the final proposed direction is checked through

$$
r_k[d_k].
$$

If it is not a finite negative value, the history is cleared and a metric-gradient
direction is used instead. This is a recurring pattern in the reduced solver layer:

> more sophisticated history is allowed to accelerate the search, but a usable
> direction still has to satisfy the basic descent contract before globalization.

## 8. Newton minimizes a local second-order model

The derivative gives the linear part of the local objective change. Newton's method also
uses the leading change of that derivative itself: the Hessian.

Around the current control $u_k$, a second-order Taylor model is

```math
m_k(s)
:=
j_h(u_k)
+
r_k[s]
+
\frac12
\langle H_k s,s\rangle,
```

where

$$
H_k
:=
j_h''(u_k):
U_h
\longrightarrow
U_h^{\ast}.
$$

The first-order term asks whether $s$ points downhill. The quadratic term says how the
slope is expected to change as we move in that direction.

To find the stationary point of this local quadratic model, differentiate it with
respect to an arbitrary perturbation $v\in U_h$:

```math
D m_k(s)[v]
=
r_k[v]
+
\langle H_k s,v\rangle.
```

Requiring this to vanish for every $v$ gives

$$
H_k s
=
-r_k.
$$

The Newton direction is therefore the step that would solve the local quadratic model
exactly.

If the true reduced objective is itself a strictly convex quadratic and the Newton
equation is solved exactly, this step reaches its minimizer in one iteration. For a
nonlinear objective, the quadratic model is only local, but near a nondegenerate
solution it can become accurate enough that Newton steps converge much faster than
first-order directions.

This is the rationale for bringing in second-order information: not merely to scale the
gradient, but to predict how the derivative will change across the step.

### 8.1 What a reduced Hessian action contains

For a PDE-constrained reduced objective, the Hessian is not just the Hessian of the
explicit control regularization term. Changing the control perturbs the state; that
perturbs the objective derivative and the adjoint; those changes contribute to

$$
j_h''(u_k)[w].
$$

A concrete reduced-Hessian provider may therefore perform tangent-state solves,
incremental-adjoint solves, or exploit problem-specific linear-quadratic
simplifications. The optimization layer does not reconstruct those second-order terms
from the first-order DTO interface. It receives an explicit action

$$
w
\longmapsto
H_k w
\in U_h^{\ast}
$$

through `ReducedHessianT`. This keeps the distinction between "first-order model
available" and "Newton-capable reduced model available" explicit.

### 8.2 The current Newton direction is computed by metric-preconditioned CG

The implementation does not form a dense reduced Hessian. Instead it runs an inner
conjugate-gradient solve using the supplied Hessian action. The inner residual begins as

$$
q_0
:=
-r_k.
$$

The metric inverse provides the preconditioned residual

$$
z_0
:=
G^{-1}q_0.
$$

The inner residual norm is measured as

$$
\lVert q\rVert_{G^{-1}}
:=
\sqrt{
\langle q,G^{-1}q\rangle
}.
$$

This is the natural dual norm induced by the same metric that identifies reduced
derivatives with primal gradients. The inner iteration repeatedly requests

$$
H_k p_i
$$

for its current primal CG search direction $p_i$. No Hessian matrix is required.

### 8.3 Positive curvature is required by this Newton–CG path

For the current inner CG solve, each search direction must satisfy a positive curvature
condition

$$
\langle H_k p_i,p_i\rangle>0
$$

up to a numerical tolerance scaled by the metric norm of $p_i$.

If the reduced Hessian has nonpositive curvature in the explored direction, this Newton
direction policy does not convert the event into an indefinite Newton step. It rejects
the inner solve.

This is an important distinction from the trust-region truncated-CG solver discussed
later, which has an explicit negative-curvature path.

### 8.4 Inner and outer convergence are different

The inner CG solve stops when its preconditioned residual reaches

$$
\max
\left(
\tau_{\mathrm{abs}},
\tau_{\mathrm{rel}}\lVert q_0\rVert_{G^{-1}}
\right).
$$

This tolerance controls how accurately the Newton equation is solved. It is not the
outer optimization stopping criterion. The outer solver still decides convergence from
quantities such as the reduced gradient norm, relative gradient norm, objective change,
or step norm.

The result records the inner iteration count and initial/final residual norms separately
from the outer histories.

## 9. A direction is not enough: globalization controls the step

All direction policies ultimately provide a primal $d_k$ and its current directional
derivative

$$
r_k[d_k].
$$

The solver requires this quantity to be negative before entering a line search. The next
question is how far to move. A fixed local model only tells us that sufficiently small
positive steps should decrease the objective. It does not say that

$$
u_k+d_k
$$

is safe. Globalization policies turn local derivative/curvature information into an
acceptance rule. The current line-search layer includes:

```text
fixed step
Armijo backtracking
weak Wolfe
strong Wolfe
exact quadratic step
```

They do not all require the same amount of trial-point information.

## 10. Armijo backtracking protects a local descent model from overshooting

A descent direction only guarantees improvement for **sufficiently small** positive
steps. To see why, restrict the reduced objective to the line

$$
\phi(\alpha)
:=
j_h(u_k+\alpha d_k).
$$

At the current point,

$$
\phi'(0)
=
r_k[d_k]
<
0.
$$

A local quadratic model along the line is

$$
\phi(\alpha)
\approx
\phi(0)
+
\alpha\phi'(0)
+
\frac12\alpha^2\kappa,
$$

where, when second-order information is available,

$$
\kappa
:=
\langle H_k d_k,d_k\rangle.
$$

If $\kappa>0$, the linear term initially drives the objective downward, but the
quadratic curvature term eventually dominates. The minimizer of this one-dimensional
quadratic model is

$$
\alpha_{\ast}
=
-\frac{\phi'(0)}{\kappa},
$$

and the quadratic model returns to its starting value at

$$
\alpha
=
-\frac{2\phi'(0)}{\kappa}.
$$

A trial step substantially larger than that can therefore **shoot past the local minimum
and increase the objective**, even though $d_k$ is a valid descent direction.

For a nonlinear objective we do not generally know the correct quadratic model over a
large step. Backtracking is a cheap way to retreat from an overambitious initial trial
until the step lies in a region where the observed objective decrease agrees with the
local descent information.

It does not try to locate the exact one-dimensional minimum.

### 10.1 The Armijo condition asks for sufficient decrease

For an unconstrained straight-line trial,

$$
u(\alpha)
=
u_k+\alpha d_k.
$$

The classical Armijo condition is

$$
j_h(u_k+\alpha d_k)
\leq
j_h(u_k)
+
c_1\alpha r_k[d_k],
$$

with

$$
0<c_1<1.
$$

Because $r_k[d_k]<0$, the right-hand side lies below the current objective for positive
$\alpha$.

The condition therefore asks the actual objective decrease to be at least a small
fraction of the decrease predicted by the first-order model.

### 10.2 Backtracking

The current Armijo policy begins from an initial step length and repeatedly multiplies
it by a factor

$$
0<\rho<1
$$

until the sufficient-decrease condition holds or the trial limit/minimum step is
reached.

Because $d_k$ is a descent direction and $j_h$ is differentiable, the first-order term
dominates for sufficiently small positive $\alpha$. Geometric reduction is therefore a
simple way of searching for a scale on which the local model becomes trustworthy without
solving a separate one-dimensional optimization problem.

A unit step is especially natural for Newton and quasi-Newton methods. Far from the
solution it may be too aggressive; close to a well-behaved solution the local quadratic
model improves, and Armijo can often accept the full step rather than continuing to damp
it.

Schematically:

```text
α = initial α
    │
    ▼
evaluate j(u + α d)
    │
    ├── sufficient decrease → accept
    │
    └── otherwise
           α ← ρ α
           repeat
```

The value/derivative split from the previous chapter matters here. A rejected Armijo
trial needs:

```text
state solve
objective value
```

but no trial adjoint. Only after the objective test succeeds does the policy augment the
retained trial with its reduced derivative.

### 10.3 The implementation uses the actual update

For an unconstrained straight-line step,

$$
s_k
:=
u_{\mathrm{trial}}-u_k
=
\alpha d_k.
$$

Then

$$
r_k[s_k]
=
\alpha r_k[d_k].
$$

The current policy writes the sufficient-decrease bound using the **actual update**

$$
s_k
$$

rather than reconstructing $\alpha d_k$ separately:

$$
j_h(u_{\mathrm{trial}})
\leq
j_h(u_k)
+
c_1 r_k[s_k].
$$

For an unconstrained line search the two forms are identical.

This formulation becomes useful for the current projected steepest-descent path, where
the projected trial need not lie exactly on the straight ray $u_k+\alpha d_k$.

## 11. A fixed step is the simplest globalization policy

The fixed-step policy uses one prescribed positive step length. It still:

- constructs the trial control;
- evaluates the reduced value;
- rejects a nonfinite objective;
- augments the derivative for a finite accepted trial.

It does not test sufficient decrease.

So "fixed step" should not be read as "no trial evaluation". The state still has to be
solved at the proposed control, because the next iteration needs a feasible reduced
evaluation.

The difference is simply that no adaptive acceptance inequality is used to choose the
step length.

## 12. Wolfe conditions also constrain the trial slope

Armijo controls objective decrease. Wolfe conditions additionally ask whether the
derivative at the trial point has changed enough. For a straight-line unconstrained path

$$
\phi(\alpha)
:=
j_h(u_k+\alpha d_k),
$$

we have

$$
\phi'(0)
=
r_k[d_k]
<0.
$$

A weak Wolfe curvature condition asks

$$
\phi'(\alpha)
\geq
c_2\phi'(0),
$$

where

$$
c_1<c_2<1.
$$

Since $\phi'(0)$ is negative, this says that the slope has become less negative: we have
moved far enough along the direction that the initial steep descent is being used rather
than taking an unnecessarily tiny step.

The strong Wolfe condition instead asks

$$
|\phi'(\alpha)|
\leq
c_2|\phi'(0)|.
$$

This prevents the trial slope from remaining too negative **or** becoming too positive
in magnitude.

### 12.1 Wolfe trials need trial derivatives

Unlike Armijo, either Wolfe condition needs

$$
j_h'(u_{\mathrm{trial}}).
$$

Therefore every tested trial is derivative-augmented:

```text
trial control
    ↓
state solve + objective
    ↓
adjoint solve + reduced derivative
    ↓
sufficient decrease + curvature test
```

A rejected Wolfe trial can therefore cost both a state and an adjoint solve. This is a
concrete example of how a mathematically stronger acceptance condition can change the
PDE work performed inside globalization.

### 12.2 The current implementation again uses the actual update

The current weak/strong Wolfe policies form

$$
s_k
:=
u_{\mathrm{trial}}-u_k
$$

and compare

$$
r_k[s_k]
\qquad\text{and}\qquad
r_{\mathrm{trial}}[s_k].
$$

For a straight line $s_k=\alpha d_k$, the common positive factor $\alpha$ cancels from
the curvature inequality, so this is equivalent to the usual directional-slope form.

For the projected steepest-descent specialization, it instead evaluates the slope along
the actual feasible update produced by the projection.

## 13. Exact quadratic line search uses the reduced Hessian

For an exactly quadratic reduced objective and a straight-line direction $d$,

$$
j(u+\alpha d)
=
j(u)
+
\alpha r[d]
+
\frac12\alpha^2\langle Hd,d\rangle.
$$

Differentiate with respect to $\alpha$:

$$
\frac{d}{d\alpha}
j(u+\alpha d)
=
r[d]
+
\alpha\langle Hd,d\rangle.
$$

If

$$
r[d]<0
$$

and

$$
\langle Hd,d\rangle>0,
$$

the minimizing step is

$$
\boxed{
\alpha_{\ast}
=
-\frac{r[d]}{\langle Hd,d\rangle}.
}
$$

The current exact-quadratic policy implements this formula with one reduced-Hessian
action.

### 13.1 Why the trial must remain on the straight line

The formula for $\alpha_{\ast}$ was derived under

$$
u_{\mathrm{trial}}
=
u+\alpha_{\ast}d.
$$

A projection or another nonlinear trial transformation would invalidate that
one-dimensional quadratic derivation. The implementation therefore explicitly checks
that the trial builder returned the expected straight-line update.

This is also why the current projected reduced solver does not list exact quadratic line
search among its supported constrained combinations.

### 13.2 It still validates the resulting objective

Even with exact quadratic curvature, the implementation evaluates the trial reduced
value and requires a finite objective that does not exceed the current objective beyond
a small tolerance.

The Hessian formula determines the requested step. The actual reduced evaluation remains
the executable evidence that the trial is usable.

## 14. Projected reduced search changes the actual step geometry

For a constrained control set $C$, Part I introduced the metric projection

$$
P_C^G.
$$

The current projected reduced search forms a trial by

$$
u_{\mathrm{trial}}
=
P_C^G(u_k+\alpha d_k).
$$

The actual update is therefore

$$
s_k
:=
u_{\mathrm{trial}}-u_k.
$$

In general,

$$
s_k\neq \alpha d_k.
$$

This is why the projected solver measures the accepted step with the metric norm of the
**actual update** and evaluates the current derivative on that same update.

### 14.1 Projected stationarity is checked before the line search

For the currently supported constrained direction, steepest descent gives

$$
d_k=-g_k.
$$

The solver constructs the unit projected point

$$
\widehat u_k
:=
P_C^G(u_k-g_k)
$$

and the projected update

$$
s_k^{\mathrm{proj}}
:=
\widehat u_k-u_k.
$$

Its stationarity measure is

$$
\lVert s_k^{\mathrm{proj}}\rVert_G.
$$

At a first-order constrained stationary point,

$$
u_k
=
P_C^G(u_k-g_k),
$$

so this norm vanishes even if the unconstrained gradient itself does not. The associated
descent measure is

$$
r_k[s_k^{\mathrm{proj}}].
$$

The line-search solver uses these projected quantities for stopping and descent
validation when a constraint is present.

### 14.2 The current projected capability is intentionally narrower

The generic unconstrained line-search machinery supports several direction policies. The
current constrained specialization accepts only:

```text
steepest descent
+
one of
    Armijo
    fixed step
    weak Wolfe
    strong Wolfe
```

It does not automatically project L-BFGS, full BFGS, nonlinear CG, Newton, or the exact
quadratic search.

That boundary is mathematically sensible to keep visible. Projecting a sophisticated
unconstrained direction does not by itself define the corresponding constrained
algorithm or preserve the assumptions behind its update formulas.

The later complementarity/PDAS chapter will present a different treatment of bounds.

## 15. Stopping criteria answer different notions of "small enough"

Optimization does not have a single universal stopping test. The current line-search
solver records several possibilities.

### Absolute stationarity

For an unconstrained solve,

$$
\lVert g_k\rVert_G
\leq
\tau_g.
$$

For the projected constrained path, the same slot is filled by the projected
stationarity norm

$$
\lVert P_C^G(u_k-g_k)-u_k\rVert_G.
$$

This is the most direct first-order stopping measure.

### Relative stationarity

The solver can compare the current stationarity measure with its initial value:

$$
\frac{\eta_k}{\eta_0}
\leq
\tau_{\mathrm{rel}},
$$

where $\eta_k$ denotes the unconstrained gradient norm or projected norm as appropriate.
This can be useful when the natural scale of the initial derivative varies across
problem instances.

### Objective change

After an accepted step, one may stop when

$$
|j_h(u_{k+1})-j_h(u_k)|
\leq
\tau_J.
$$

A small objective change can indicate practical stagnation, but it is not by itself a
first-order stationarity statement.

### Step norm

Likewise,

$$
\lVert u_{k+1}-u_k\rVert_G
\leq
\tau_s
$$

can detect a search that is no longer moving significantly.

Again, a tiny step can arise because the problem is solved or because globalization has
become restrictive. The stopping reason therefore matters when interpreting the result.

### Objective target

The line-search solver can also stop once

$$
j_h(u_k)
\leq
J_{\mathrm{target}}.
$$

This is an application-level target rather than a stationarity criterion.

### Iteration and globalization limits

Finally, a run can terminate because:

- the maximum accepted-iteration count was reached;
- the line search failed to find an acceptable trial.

Those outcomes are qualitatively different from convergence and are represented by
different stopping reasons.

## 16. "Automatic" stopping combines an absolute baseline with optional tests

The current `automatic` policy preserves the project's historical behavior. The absolute
stationarity tolerance is always active. Positive optional tolerances for

```text
relative gradient/stationarity
objective change
step norm
```

add further stopping conditions. This makes the effective stopping rule a disjunction:

```text
stop if
    absolute stationarity is small enough

or an enabled optional criterion is satisfied
```

Selecting one explicit `ReducedStoppingCriterion` instead gives that criterion its own
required positive tolerance.

The distinction is worth knowing when reproducing a run: the same numerical tolerance
fields can mean "additional automatic stops" or "the selected stop" depending on the
policy.

## 17. The solver records work because PDE iterations are not equal-cost iterations

In ordinary low-dimensional optimization one often reports only the number of outer
iterations. For PDE-constrained optimization that count can hide most of the
computational work.

One outer iteration may contain:

- several state solves for line-search trials;
- one or several adjoint solves;
- metric inverse applications;
- reduced-Hessian actions;
- an inner Newton or trust-region solve.

The reduced solver therefore tracks work explicitly. The common counters include:

```text
state solves
adjoint solves
metric solves
reduced-Hessian actions
direction/history resets
line-search trials
accepted outer iterations
```

### 17.1 Armijo and Wolfe can have the same outer count but different PDE cost

Suppose both methods accept their third trial. An Armijo iteration may perform:

```text
3 state solves
1 adjoint solve
```

because only the accepted value is derivative-augmented. A Wolfe iteration may perform:

```text
3 state solves
3 adjoint solves
```

because every curvature test needs a trial derivative. The accepted-iteration count is
one in both cases. Work counters reveal the difference.

### 17.2 Newton adds inner work

A Newton direction can require:

```text
metric inverse applications
+
several reduced-Hessian actions
```

before the outer line search even begins. The result therefore records
`ReducedHessianSolveDiagnostics` alongside the outer iteration history:

```text
inner iteration count
initial inner residual norm
final inner residual norm
```

This keeps the cost and quality of the direction solve visible.

## 18. Accepted-iteration records retain more than histories of scalars

The solver still exposes compatibility histories such as:

```text
objective history
stationarity-norm history
step-length history
step-norm history
objective-change history
```

but it also builds a record for each accepted iteration. A line-search record includes
information such as:

- objective before and after;
- requested step parameter;
- actual step norm;
- actual descent pairing;
- absolute and relative stationarity;
- number of trials;
- work performed during the iteration;
- cumulative work;
- Hessian inner-solve diagnostics;
- the complete accepted reduced evaluation;
- policy-specific acceptance evidence.

This matters because a single scalar like "step length 0.5" does not fully describe the
accepted update when projection is active, and an objective sequence does not explain
how much PDE work produced it.

### 18.1 Trial records expose rejected work too

The line-search layer also retains trial-level records for policies that populate them.
For Armijo, these include:

```text
trial step
trial objective
actual current slope
sufficient-decrease bound
finite-objective flag
negative-slope flag
accepted flag
```

Rejected trial work is therefore not invisible merely because it did not become an
accepted iterate. This evidence model is one reason the solver code is more elaborate
than the mathematical pseudocode alone would suggest.

## 19. Trust regions globalize a quadratic model instead of a direction

A trust-region method begins from a local quadratic model of the reduced objective. At
$u_k$, write

```math
m_k(s)
:=
j_h(u_k)
+
r_k[s]
+
\frac12
\langle H_k s,s\rangle,
```

where $s\in U_h$ is a candidate step.

Rather than choosing a direction and then a line-search length, the method asks for a
step that approximately minimizes this model inside

$$
\lVert s\rVert_G
\leq
\Delta_k.
$$

The radius $\Delta_k$ expresses how far the algorithm currently trusts the quadratic
model. Here the roles are

$$
r_k\in U_h^{\ast},
\qquad
s\in U_h,
\qquad
H_k:
U_h\longrightarrow U_h^{\ast}.
$$

The trust-region subproblem is therefore

```math
\min_{s\in U_h}
\hspace{0.5em}
r_k[s]
+
\frac12\langle H_k s,s\rangle
\quad
\text{subject to}
\quad
\lVert s\rVert_G\leq\Delta_k.
```

Both terms in the objective are scalars obtained by pairing covectors with the primal
step $s$. The constant $j_h(u_k)$ can be omitted from the subproblem because it does not
affect the minimizer.

### 19.1 Trust regions require explicit reduced-Hessian actions in the current implementation

The current `ReducedTrustRegionSolverT` receives:

```text
ReducedDTOT
MetricT
ReducedHessianT
```

It does not build a Hessian approximation from BFGS history. The Hessian action is used
both to form the local model and, for truncated CG, to solve the subproblem.

The current trust-region implementation is also unconstrained; projection remains a
separate boundary.

## 20. The Cauchy step minimizes the model along the negative gradient ray

The simplest trust-region subproblem method restricts attention to

$$
s=-t g_k,
\qquad
t\geq0.
$$

Recall that

$$
r_k[g_k]
=
\lVert g_k\rVert_G^2.
$$

Substitute $s=-tg_k$ into the quadratic model change:

```math
m_k(-t g_k)-j_h(u_k)
=
-t\langle r_k,g_k\rangle
+
\frac12t^2
\langle H_k g_k,g_k\rangle.
```

Define

$$
a
:=
\langle r_k,g_k\rangle
=
\lVert g_k\rVert_G^2
$$

and

$$
c
:=
\langle H_k g_k,g_k\rangle.
$$

If

$$
a>0,
\qquad
c>0,
$$

the unconstrained minimizer along this ray is

$$
t_{\ast}
=
\frac{a}{c}.
$$

But the trust-region bound requires

$$
\lVert -t g_k\rVert_G
=
t\lVert g_k\rVert_G
\leq
\Delta_k.
$$

Hence

$$
t
\leq
\frac{\Delta_k}{\lVert g_k\rVert_G}.
$$

The Cauchy step uses

$$
t_C
:=
\min
\left(
\frac{a}{c},
\frac{\Delta_k}{\lVert g_k\rVert_G}
\right)
$$

and

$$
s_C
=
-t_C g_k.
$$

The predicted reduction is

$$
\mathrm{pred}_k
=
t_C a
-
\frac12 t_C^2 c.
$$

This is exactly the structure implemented by the current Cauchy subproblem path.

### 20.1 Why positive curvature is required here

The formula

$$
t_{\ast}=a/c
$$

assumes

$$
c>0.
$$

The current Cauchy implementation therefore rejects a nonpositive-curvature model along
the metric-gradient direction. The truncated-CG path treats negative curvature
differently.

## 21. Truncated CG searches the trust-region quadratic more deeply

For a positive-definite quadratic model with no radius boundary, conjugate gradients can
solve

$$
H_k s=-r_k
$$

without forming the Hessian matrix. The trust-region version follows a similar Krylov
process but stops early when one of several events occurs:

```text
inner residual converges
trust-region boundary is reached
negative curvature is detected
inner iteration limit is reached
```

This is commonly called truncated CG or Steihaug-type CG. The current implementation
uses the metric inverse as a preconditioner and measures the subproblem residual in the
associated dual norm.

### 21.1 Negative curvature sends the step to the boundary

Suppose the current CG search direction is $p$ and

$$
\langle H_kp,p\rangle
\leq0.
$$

Moving farther along $p$ does not produce the positive curvature assumed by ordinary CG.
Inside a trust region, that does not have to be treated as an unrecoverable failure.

The implementation instead chooses the point

$$
s+\tau p
$$

where the current ray first reaches

$$
\lVert s+\tau p\rVert_G
=
\Delta_k.
$$

The scalar $\tau$ is found from the quadratic equation

```math
\lVert s+\tau p\rVert_G^2
=
\lVert s\rVert_G^2
+
2\tau(s,p)_G
+
\tau^2\lVert p\rVert_G^2
=
\Delta_k^2.
```

The positive boundary intersection is

```math
\tau
=
\frac{
-(s,p)_G
+
\sqrt{
(s,p)_G^2
-
\lVert p\rVert_G^2
\left(
\lVert s\rVert_G^2-\Delta_k^2
\right)
}
}{
\lVert p\rVert_G^2
}.
```

The same boundary construction is used if an ordinary positive-curvature CG step would
leave the trust region.

This is the major robustness difference from the current Newton–CG line-search
direction: Newton–CG requires positive curvature, while trust-region truncated CG can
return a boundary step when negative curvature is encountered.

## 22. Predicted reduction is compared with actual reduction

Once a trust-region step $s_k$ has been computed, the local model predicts

$$
\mathrm{pred}_k
:=
-
r_k[s_k]
-
\frac12
\langle H_ks_k,s_k\rangle.
$$

The actual reduced evaluation gives

$$
\mathrm{ared}_k
:=
j_h(u_k)-j_h(u_k+s_k).
$$

The agreement ratio is

$$
\rho_k
:=
\frac{
\mathrm{ared}_k
}{
\mathrm{pred}_k
}.
$$

This number asks:

> How well did the local quadratic model predict the objective improvement of the
> step it proposed?

A ratio near one means the model predicted the change well. A small or negative ratio
means the actual objective behaved much worse than the model suggested.

### 22.1 Acceptance and radius updates

The current implementation accepts a trial when

$$
\rho_k
$$

is finite and at least the configured acceptance threshold. Rejected trials shrink the
radius and solve another subproblem around the **same current iterate**.

After an accepted step:

- a high ratio together with a near-boundary step can expand the radius;
- a low ratio can shrink it;
- otherwise the radius is retained.

This gives the trust-region method a feedback loop:

```text
model predicted well
    → trust it farther

model predicted poorly
    → trust it less
```

### 22.2 Rejected trust-region trials can also avoid adjoint solves

A trust-region trial needs the actual objective value to compute

$$
\mathrm{ared}_k.
$$

That requires a state solve. If the ratio rejects the trial, the code does not need the
trial reduced derivative. Only an accepted trial is derivative-augmented.

Thus the retained value/augmentation split from the formulation chapter benefits trust
regions as well as Armijo line searches:

```text
rejected trust-region trial
    state solve + objective

accepted trust-region trial
    state solve + objective
    + adjoint + reduced derivative
```

## 23. Line search and trust region retain different acceptance evidence

The line-search family records evidence appropriate to its policy:

```text
Armijo
    sufficient-decrease bound

Wolfe
    sufficient-decrease + trial slope/curvature condition

exact quadratic
    curvature used for the analytic step
```

Trust-region records instead include:

```text
radius
predicted reduction
actual reduction
reduction ratio
subproblem status
subproblem iteration count
subproblem residual norm
```

These are not merely logging variations. They correspond to different mathematical
questions:

```text
line search
    was this amount of movement along the proposed direction acceptable?

trust region
    did this local quadratic model predict the accepted step accurately enough?
```

## 24. The current capability map

The implemented reduced optimization family can be summarized without treating every
type alias as a separate algorithm.

### Direction policies for the line-search solver

- metric steepest descent;
- nonlinear conjugate gradient with Polak–Ribière+;
- Fletcher–Reeves;
- a classical quadratic-CG specialization;
- limited-memory BFGS;
- full-memory BFGS;
- Newton with an explicit reduced-Hessian action.

### Line-search policies

- Armijo backtracking;
- fixed step;
- weak Wolfe;
- strong Wolfe;
- exact quadratic line search requiring a reduced Hessian and a straight-line trial.

### Constrained projection path

Currently:

```text
steepest descent
+
Armijo / fixed step / weak Wolfe / strong Wolfe
```

provided the constraint can project in the selected metric.

### Trust-region path

Currently:

```text
unconstrained
+
explicit ReducedHessianT
+
Cauchy or truncated-CG subproblem
```

The existence of all these pieces in the repository should not be read as saying that
every Cartesian product of direction, globalization, and constraint policy is
implemented.

The current code has explicit supported combinations and checks them.

## 25. Following one outer iteration through `ReducedSearchSolverT`

With the mathematics in place, the main line-search solver becomes much easier to read.
Its outer loop performs roughly:

```text
current ReducedEvaluation
        │
        ▼
direction_policy.next(...)
        │
        ├── metric work
        └── maybe history/Hessian work
        │
        ▼
stationarity checks
        │
        ▼
line_search_policy.search(...)
        │
        ├── build trial
        ├── evaluate_value
        ├── maybe augment_derivative
        └── accept/fail
        │
        ▼
measure actual accepted step
        │
        ▼
record work + acceptance evidence
        │
        ▼
post-step stopping checks
```

This is almost the entire orchestration of the line-search solver. The mathematical
complexity lives mainly inside the direction and line-search policies.

### 25.1 Direction policies are copied and reset at the start of a solve

History-bearing policies such as nonlinear CG and BFGS maintain state across iterations.
At the beginning of `solve()`, the configured direction policy is copied into a local
instance and reset.

That means history belongs to one solver run, not to a long-lived policy object shared
across independent solves. This source detail matches the mathematical expectation that
secant/conjugacy history is meaningful only for one sequence of iterates.

## 26. Source orientation

A practical reading order is now:

### Shared search vocabulary and direction policies

- [`include/nmopt/solvers/reduced_search.hpp`](../../../include/nmopt/solvers/reduced_search.hpp)

Read:

```text
ReducedSolverParameters
ReducedSearchDirectionT
make_metric_gradient
SteepestDescentDirectionPolicyT
NonlinearConjugateGradientDirectionPolicyT
LimitedMemoryBfgsDirectionPolicyT
FullBfgsDirectionPolicyT
NewtonDirectionPolicyT
```

The file contains the mathematical direction machinery and common evidence records.

### Line-search policies

- [`include/nmopt/solvers/reduced_line_search.hpp`](../../../include/nmopt/solvers/reduced_line_search.hpp)

Read Armijo first because it exposes the value/derivative split most clearly. Then
compare fixed step, weak/strong Wolfe, and exact quadratic search.

### Outer line-search solver

- [`include/nmopt/solvers/reduced_gradient.hpp`](../../../include/nmopt/solvers/reduced_gradient.hpp)

Despite the historical filename, the central class is the more general
`ReducedSearchSolverT`. Follow:

```text
initial reduced evaluation
direction
stopping
trial builder
line search
accepted update
audit record
```

### Trust-region solver

- [`include/nmopt/solvers/reduced_trust_region.hpp`](../../../include/nmopt/solvers/reduced_trust_region.hpp)

Read the Cauchy subproblem before truncated CG, then return to the outer ratio/radius
loop.

### Formulation boundary

When a solver call to `evaluate`, `evaluate_value`, or `augment_derivative` becomes
opaque, return to:

- [`include/nmopt/contract/reduced_dto.hpp`](../../../include/nmopt/contract/reduced_dto.hpp)

The optimization layer intentionally treats the PDE evaluation mechanism through that
reduced interface.

## 27. A compact correspondence

| Optimization idea | Mathematical object | Current code role |
| --- | --- | --- |
| reduced derivative | $r=j_h'(u)\in U_h^{\ast}$ | `ReducedEvaluationT::reduced_derivative` |
| metric gradient | $g=G^{-1}r$ | `make_metric_gradient` |
| steepest direction | $d=-g$ | `SteepestDescentDirectionPolicyT` |
| nonlinear CG | $d_k=-g_k+\beta_kd_{k-1}$ | nonlinear-CG direction policies |
| secant data | $s_k,\ \Delta r_k:=r_{k+1}-r_k$ | BFGS history |
| Newton equation | $H_kd=-r_k$ | `NewtonDirectionPolicyT` |
| line-search trial | $u+\alpha d$ or projected update | trial-control builder |
| Armijo | sufficient decrease | `ArmijoLineSearchPolicyT` |
| Wolfe | decrease + slope condition | weak/strong Wolfe policies |
| exact quadratic step | $-r[d]/\langle Hd,d\rangle$ | `ExactQuadraticLineSearchPolicyT` |
| trust-region model | $r[s]+\frac12\langle Hs,s\rangle$ | `ReducedTrustRegionSolverT` |
| trust radius | $\lVert s\rVert_G\leq\Delta$ | metric trust-region bound |
| model agreement | $\rho=\mathrm{ared}/\mathrm{pred}$ | trust-region acceptance/radius update |

## 28. What this completes, and what comes next

The two reduced chapters now form one complete story. [Reduced state–adjoint
formulation](reduced-state-adjoint-formulation.md) explained how a control produces

$$
j_h(u)
\qquad\text{and}\qquad
j_h'(u).
$$

This chapter explained how reduced optimization methods repeatedly consume those
objects:

```text
formulation
    control → reduced value/derivative
          │
          ▼
optimization
    direction + globalization + stopping
          │
          ▼
    next control
```

The next conceptual shift is substantial. Reduced methods eliminate the PDE state
through an inner state solve.

The forthcoming **Optimality systems and KKT** chapter will instead retain state,
control, and multiplier variables together and derive the coupled first-order system
that an all-at-once method solves.

That chapter will therefore return to the Lagrangian introduced in the reduced
derivation, but it will stop eliminating the state.

## Read later

Useful existing documents are:

- [Theoretical formalism](../../design/theoretical-formalism.md), for the formal
  reduced derivative, Hessian, metric, and optimization conventions.
- [Reduced optimization](../../overview/reduced-optimization.md), for the shorter
  project-wide runtime view.
- [Chapter 6 numerical methods](../../guides/chapter-6-numerical-methods.md), for the
  source-text numerical-method context.
- [Chapter 6 numerical examples](../../guides/chapter-6-numerical-examples.md), for
  the benchmark/problem families on which selected reduced paths are exercised.

The later **Verification and evidence** chapter will revisit iteration records, work
counts, solve evidence, line-search acceptance evidence, and Hessian diagnostics as part
of a broader hierarchy of numerical evidence.
