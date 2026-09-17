# Complementarity and PDAS

## Bounds turn stationarity into a switching condition

The previous chapter ended with an equality-constrained quadratic problem

```math
\begin{aligned}
\min_x\quad&
\frac12\langle Qx,x\rangle
-
\langle c_x,x\rangle,
\\
\text{subject to}\quad&
Dx=d,
\end{aligned}
```

and the corresponding KKT equations

```math
\begin{bmatrix}
Q & D^{\mathsf T}\\
D & 0
\end{bmatrix}
\begin{bmatrix}
x\\
\lambda
\end{bmatrix}
=
\begin{bmatrix}
c_x\\
d
\end{bmatrix}.
```

Now suppose one block of the primal variable is a control $u$ subject to componentwise
bounds

$$
\ell_i
\leq
u_i
\leq
r_i.
$$

The equality equations remain. What changes is the control stationarity condition. For a
free control coefficient, the ordinary stationarity residual should vanish. For a
coefficient fixed at its lower bound, the derivative need not vanish: it may point out
of the admissible interval.

Likewise at an upper bound, a nonzero derivative can be compatible with optimality if it
points out through the upper side. So one equation,

$$
\text{control stationarity}=0,
$$

is replaced by a switching rule:

```text
inside the box
    stationarity must vanish

at the lower bound
    one multiplier sign is allowed

at the upper bound
    the opposite multiplier sign is allowed
```

This is the complementarity structure behind the project's primal-dual active-set (PDAS)
formulation. The chapter develops that structure in four stages:

```text
box-constrained first-order condition
        ↓
normal-cone / multiplier form
        ↓
active, lower-active, upper-active, and free coefficients
        ↓
repeat:
    fix active coefficients at their bounds
    solve a restricted equality-constrained KKT problem
    reconstruct multipliers
    reclassify
```

The equality-constrained KKT product from the previous chapter is therefore not replaced
by PDAS. It becomes the linear algebraic core solved under a changing active-set
hypothesis.

## What PDAS is trying to accomplish

The most useful way to understand PDAS is to imagine that somebody has already told us
which control coefficients are on their bounds.

Suppose, for example, that a control has three coefficients and the correct pattern were

```text
u_0
    lower active

u_1
    inactive

u_2
    upper active
```

Then two unknowns would disappear immediately:

$$
u_0=\ell_0,
\qquad
u_2=r_2.
$$

Only $u_1$ would remain as a free control unknown. The state variables and equality
multipliers would still be unknown, but the resulting problem would be an ordinary
**linear equality-constrained KKT system** on that reduced set of variables.

So, if the active set were known, the difficult bound-constrained problem would be easy
to state:

```text
fix active coefficients at their bounds
        ↓
solve one ordinary KKT system for the remaining variables
```

The real difficulty is that the active set is not known beforehand. PDAS therefore
treats the active set itself as the unknown combinatorial information. It repeatedly:

```text
guess lower / inactive / upper
        ↓
solve the KKT problem implied by that guess
        ↓
inspect the resulting control and multipliers
        ↓
correct the guess
```

The multiplier is what tells us whether fixing a coefficient at a bound was consistent.
If a coefficient was declared inactive but the new stationarity information strongly
pushes it through a bound, the next classification makes it active.

Conversely, if a coefficient was declared active but the reconstructed multiplier has
the wrong sign, the next classification can release it. This gives the method a concrete
purpose:

> identify the correct active bounds while solving only linear
> equality-constrained KKT problems for each current guess.

Once the active set stops changing and the KKT, feasibility, and complementarity
residuals are small, the current point satisfies the full bound-constrained first-order
conditions to the requested tolerances.

The sign and complementarity formulas in the next sections are the machinery that
answers one practical question:

> after solving the current KKT subproblem, should each control coefficient stay
> lower-active, stay upper-active, or become inactive?

## 1. Start from a bound-constrained first-order condition

To isolate the new idea, first ignore the PDE equality constraint and consider

$$
\min_{u\in C} j(u),
$$

where

$$
C
:=
\left\{
u\in U_h:
\ell_i\leq u_i\leq r_i
\right\}.
$$

At a local minimum $u^\ast$ of a differentiable objective over the convex set $C$,

$$
j'(u^\ast)[v-u^\ast]
\geq0
\qquad
\text{for every }v\in C.
$$

This is the variational inequality already encountered in [Metrics, gradients, and
constraints](metrics-gradients-and-constraints.md). It says that no feasible first-order
perturbation decreases the objective.

### 1.1 What this means for one coefficient

Consider one scalar coefficient

$$
\ell\leq u\leq r.
$$

Let the derivative coefficient be $\rho$. If

$$
\ell<u<r,
$$

then both positive and negative small perturbations are feasible. The variational
inequality can hold for both signs only if

$$
\rho=0.
$$

At the lower bound,

$$
u=\ell,
$$

only nonnegative perturbations are locally feasible. Therefore

$$
\rho\geq0
$$

is compatible with optimality. At the upper bound,

$$
u=r,
$$

only nonpositive perturbations are locally feasible. Therefore

$$
\rho\leq0
$$

is compatible with optimality. This is one common derivative-sign convention.

The PDAS contract uses an equivalent **box multiplier** convention in which the
multiplier is added to the unconstrained stationarity residual. We now derive that sign
carefully.

## 2. Introduce a box multiplier through the normal cone

Let

$$
r_u\in U_h^{\ast}
$$

denote the control stationarity covector before adding the box constraint. For the
quadratic KKT product, this is the control block of

$$
Qx+D^{\mathsf T}\lambda-c_x.
$$

The constrained stationarity equation is written

$$
r_u+\mu=0,
$$

where

$$
\mu\in U_h^{\ast}
$$

is the box multiplier. The sign convention is:

```text
lower bound
    mu_i <= 0

free/interior
    mu_i = 0

upper bound
    mu_i >= 0
```

Why these signs? Because

$$
\mu=-r_u.
$$

At the lower bound the unconstrained derivative can satisfy

$$
r_{u,i}\geq0,
$$

so

$$
\mu_i\leq0.
$$

At the upper bound,

$$
r_{u,i}\leq0,
$$

so

$$
\mu_i\geq0.
$$

This is the convention used by the current complementarity and PDAS code.

### 2.1 The normal-cone statement

Define the normal cone with the sign convention

$$
N_C(u)
:=
\left\{
\mu\in U_h^{\ast}:
\langle\mu,v-u\rangle\leq0
\quad
\text{for every }v\in C
\right\}.
$$

Then the first-order condition is

$$
0\in r_u+N_C(u).
$$

Equivalently, there exists

$$
\mu\in N_C(u)
$$

such that

$$
r_u+\mu=0.
$$

For a scalar interval $[\ell,r]$, this normal cone is

```math
N_{[\ell,r]}(u)
=
\begin{cases}
(-\infty,0], & u=\ell,\\
\{0\}, & \ell<u<r,\\
[0,\infty), & u=r.
\end{cases}
```

That is exactly the sign pattern above.

## 3. Complementarity expresses "bound active or multiplier zero"

The normal-cone condition can be written with lower and upper nonnegative multipliers.
Define

$$
\mu_i^{+}
:=
\max(\mu_i,0),
$$

and

$$
\mu_i^{-}
:=
\max(-\mu_i,0).
$$

Then

$$
\mu_i
=
\mu_i^{+}-\mu_i^{-},
$$

with

$$
\mu_i^{+}\geq0,
\qquad
\mu_i^{-}\geq0.
$$

The box conditions become

```math
\begin{aligned}
\ell_i\leq u_i\leq r_i,
\\
\mu_i^{-}(u_i-\ell_i)=0,
\\
\mu_i^{+}(r_i-u_i)=0.
\end{aligned}
```

These are the complementary-slackness relations. They say:

```text
lower side
    either u_i > ell_i and mu_i^- = 0
    or     u_i = ell_i and mu_i^- may be positive

upper side
    either u_i < r_i and mu_i^+ = 0
    or     u_i = r_i and mu_i^+ may be positive
```

A coefficient cannot carry a nonzero bound multiplier while remaining strictly inside
the box.

### 3.1 The project's signed multiplier form

The current PDAS path does not store separate lower and upper multiplier vectors. It
stores one signed dual multiplier

$$
\mu\in U_h^{\ast}.
$$

For classification and complementarity diagnostics, that covector is converted to the
declared primal representative

$$
m:=G^{-1}\mu.
$$

The implementation then measures the two scalar products

```math
\max(m_i,0)(r_i-u_i)
```

and

```math
\min(m_i,0)(u_i-\ell_i).
```

At an exact complementary point both vanish. The maximum absolute violation over all
control coefficients becomes the reported complementarity residual. In the declared
primal multiplier representation, the signs mean:

```text
m_i < 0
    lower-side multiplier

m_i = 0
    free coefficient

m_i > 0
    upper-side multiplier
```

The dual covector $\mu$ remains the quantity added to the KKT stationarity residual. The
representative $m$ is the quantity compared coefficientwise with bounds.

## 4. A dual multiplier cannot be compared directly with a primal bound distance

The box distances

$$
u_i-\ell_i
\qquad\text{and}\qquad
u_i-r_i
$$

are primal quantities. The multiplier

$$
\mu\in U_h^{\ast}
$$

is a covector. It would be mathematically sloppy to add these objects directly. The
complementarity contract therefore carries a **multiplier representation** that maps
between dual multipliers and primal representatives.

For the metric-based representation used by the current implementation,

$$
G:
U_h
\longrightarrow
U_h^{\ast}
$$

is the chosen metric/Riesz map, and the primal representative is

$$
m
:=
G^{-1}\mu
\in U_h.
$$

Conversely,

$$
\mu
=
Gm.
$$

This is the role of `BoxMultiplierRepresentationT`.

### 4.1 Why the representation owns the metric

The conversion is not just a compatible vector cast. It depends on the particular metric
realization. The helper

```text
make_metric_multiplier_representation(metric)
```

therefore retains:

```text
primal layout
dual layout
dual-to-primal action
primal-to-dual action
owned metric
metric realization witness
```

The witness is checked against the owned metric. This is the same realization-identity
idea introduced for metric-dependent projection in Part I. The classification algorithm
can then work with the primal representative $m$ while the stationarity equation retains
the true dual multiplier $\mu$.

## 5. The active-set classification combines position and multiplier information

Let

$$
m:=G^{-1}\mu.
$$

For a positive classification parameter

$$
c>0,
$$

the current `BoxComplementarityT::classify()` evaluates

```math
m_i+c(u_i-r_i)
```

and

```math
m_i+c(u_i-\ell_i).
```

The coefficient is classified as:

```text
upper active
    if m_i + c (u_i - r_i) > 0

lower active
    else if m_i + c (u_i - ell_i) < 0

inactive/free
    otherwise
```

This is a standard primal-dual active-set classification.

### 5.1 Why the upper test uses $u_i-r_i$, not $r_i-u_i$

Earlier, complementarity used the **upper slack**

$$
r_i-u_i.
$$

For a feasible control this quantity is nonnegative. It is exactly the right quantity
for complementary slackness:

$$
\mu_i^{+}(r_i-u_i)=0.
$$

The classification test has a different job. It needs one signed expression whose sign
decides whether the upper bound should be treated as active.

Write the upper classification expression as

$$
m_i+c(u_i-r_i)
=
m_i-c(r_i-u_i).
$$

Now the roles are visible. The first term,

$$
m_i,
$$

is the multiplier signal. A positive value supports upper activity. The second term,

$$
-c(r_i-u_i),
$$

penalizes that hypothesis when the coefficient lies well inside the feasible box. Thus:

```text
at the upper bound
    r_i - u_i = 0
    -> the multiplier sign decides

well below the upper bound
    r_i - u_i > 0
    -> the negative slack term opposes upper activity

above the upper bound
    u_i - r_i > 0
    -> the signed distance reinforces upper activity
```

So the two expressions are not inconsistent:

```text
r_i - u_i
    nonnegative slack used in complementary slackness

u_i - r_i
    signed bound residual used in the active-set classifier
```

The lower test works analogously. Its feasible slack is

$$
u_i-\ell_i\geq0,
$$

while the classifier asks whether the negative multiplier signal is strong enough that

$$
m_i+c(u_i-\ell_i)<0.
$$

### 5.2 Why this identifies an exact upper-active coefficient

At an exact upper-active coefficient,

$$
u_i=r_i,
\qquad
m_i>0.
$$

Then

$$
m_i+c(u_i-r_i)
=
m_i
>
0.
$$

So the coefficient is classified upper active for every $c>0$.

### 5.3 Why it identifies an exact lower-active coefficient

At an exact lower-active coefficient,

$$
u_i=\ell_i,
\qquad
m_i<0.
$$

Then

$$
m_i+c(u_i-\ell_i)
=
m_i
<
0.
$$

So it is classified lower active.

### 5.4 Why it leaves an interior complementary coefficient free

At an exact interior coefficient,

$$
\ell_i<u_i<r_i,
\qquad
m_i=0.
$$

Therefore

$$
m_i+c(u_i-r_i)
=
c(u_i-r_i)
<
0,
$$

while

$$
m_i+c(u_i-\ell_i)
=
c(u_i-\ell_i)
>
0.
$$

Neither active inequality is satisfied. The coefficient is therefore inactive/free.

### 5.5 What the parameter $c$ does

At an exact complementary solution, any positive $c$ produces the same classification.
Away from the solution, $c$ changes the balance between:

- the multiplier signal $m_i$;
- the geometric distance to a bound.

So $c$ affects transient active-set guesses. It does not change the box itself. The
current `PDASPolicy` therefore treats the classification parameter as an algorithmic
positive finite parameter rather than as problem data.

## 6. The classification is closely related to a projection fixed point

Part I characterized constrained stationarity through

$$
u
=
P_C^G(u-g),
$$

where

$$
g
=
G^{-1}r_u
$$

is the metric gradient of the unconstrained stationarity covector. Since

$$
r_u+\mu=0,
$$

we have

$$
G^{-1}\mu
=
-G^{-1}r_u,
$$

or

$$
m=-g.
$$

Thus the same first-order condition can be viewed through either:

```text
projection
    u = P_C^G(u - g)

or multiplier stationarity
    r_u + mu = 0
    together with box complementarity
```

The two viewpoints encode the same constrained first-order idea, but the numerical
algorithms use them differently. A projected-gradient method repeatedly projects trial
points.

PDAS instead guesses which coefficients should be fixed at a bound, solves the remaining
equality-constrained problem, and then revises that guess.

### 6.1 Why the implementation needs an explicit multiplier representation

For a simple diagonal metric, the projection and active-set pictures are especially
transparent because

$$
m_i
=
(G^{-1}\mu)_i
$$

depends only on the corresponding multiplier coordinate. For a more general multiplier
representation, the classification should be understood as being defined **in that
declared representation**.

`BoxComplementarityT` therefore does not pretend that there is a representation-free
coefficientwise comparison between a dual multiplier and a primal bound distance. The
producer supplies the conversion that gives those comparisons their numerical meaning.

This is analogous to the earlier principle:

> a derivative is intrinsic to the problem, while a gradient or primal
> representative depends on a chosen primal–dual identification.

## 7. `BoxBoundsT` and `ActiveSetSelectionT` keep the box algebra small

The complementarity contract separates three pieces of information.

### Bounds

`BoxBoundsT` stores

$$
\ell
\qquad\text{and}\qquad
r
$$

as primal values on one control layout and checks

$$
\ell_i\leq r_i
$$

componentwise. It also provides a direct feasibility test.

### Activity labels

Each control coefficient receives one label:

```text
lower
inactive
upper
```

These are represented by `BoxActivity`.

### Selection

`ActiveSetSelectionT` turns those labels into two index sets:

$$
\mathcal I
:=
\{i:\text{inactive}\},
$$

and

$$
\mathcal A
:=
\{i:\text{lower or upper active}\}.
$$

It records the actual lower/upper label for active coordinates and provides simple
restriction/prolongation operations. The current selection supports exactly one control
block.

### 7.1 Restriction and prolongation

Given a full control vector $u$, the selection can extract the inactive and active
coordinates,

```text
u_I
u_A
```

by index. Conversely, given inactive and active subvectors, it can reconstruct the full
coefficient vector in the original order. If, for example,

```text
activities
    [ upper | lower | free | free ]
```

then

```text
inactive indices
    [ 2, 3 ]

active indices
    [ 0, 1 ]
```

and prolongation inserts the two supplied active values back into positions 0 and 1. The
backend-neutral complementarity test checks exactly this behavior.

## 8. What PDAS assumes during one iteration

Suppose the current classification is

$$
(\mathcal A^{-},\mathcal I,\mathcal A^{+}),
$$

where $\mathcal A^{-}$ is lower-active, $\mathcal I$ is inactive/free, and $\mathcal
A^{+}$ is upper-active.

PDAS temporarily assumes that this classification is correct. The point of the iteration
is to solve the easiest problem consistent with that assumption and then use the result
to decide whether the assumption should be kept.

Under that assumption:

```text
i in A^-
    set u_i = ell_i

i in A^+
    set u_i = r_i

i in I
    leave u_i unknown
    and impose ordinary control stationarity
```

So the active control coefficients stop being unknowns in the current KKT solve. Only
the inactive control coordinates remain in the stationarity system.

This converts the bound-constrained problem into an equality-constrained problem on a
smaller primal space.

### 8.1 Why this is powerful for a quadratic problem

For a quadratic objective and linear equality constraints, once the active set is fixed,
the remaining subproblem is again a quadratic equality-constrained problem.

There is no nonlinear objective left to approximate. The subproblem can therefore be
represented by the same KKT vocabulary as the previous chapter:

```text
restricted Q
restricted D
restricted D^T
restricted stationarity rhs
shifted equality rhs
```

The active-set iteration is nonlinear only because the identity of the active
coefficients is not known in advance. Within one fixed classification, the solve is a
linear KKT solve.

## 9. Fixing active variables changes the right-hand side

This point is easy to overlook. Suppose the base primal vector is split abstractly into
inactive and active coordinates,

$$
x
=
\begin{bmatrix}
x_I\\
x_A
\end{bmatrix},
$$

and the active values are fixed to

$$
x_A=\bar x_A.
$$

Consider a linear equality

$$
Dx=d.
$$

Partition $D$ accordingly:

$$
D
=
\begin{bmatrix}
D_I & D_A
\end{bmatrix}.
$$

Then

```math
D_Ix_I+D_A\bar x_A=d.
```

Therefore the inactive-coordinate problem must satisfy

$$
D_Ix_I
=
d-D_A\bar x_A.
$$

The active values do not simply disappear. Their contribution becomes an affine shift of
the equality right-hand side. The same effect occurs in the stationarity equations
whenever $Q$ couples free and active coordinates.

### 9.1 A tiny quadratic example

Take

$$
Q
=
\begin{bmatrix}
q_{II} & q_{IA}\\
q_{AI} & q_{AA}
\end{bmatrix},
$$

and fix $x_A=\bar x_A$. The inactive-coordinate stationarity component contains

$$
q_{II}x_I
+
q_{IA}\bar x_A
+
D_I^{\mathsf T}\lambda
-
c_I
=
0.
$$

Move the known active contribution to the right:

$$
q_{II}x_I
+
D_I^{\mathsf T}\lambda
=
c_I-q_{IA}\bar x_A.
$$

So the active-set KKT problem is not obtained by merely deleting rows and columns. The
fixed values also shift the affine terms.

## 10. `ActiveSetKKTSubproblemT` performs exactly that restriction

The current `ActiveSetKKTSubproblemT` receives:

```text
base quadratic KKT product
box complementarity object
active-set selection
index of the control block
assumptions for the restricted KKT product
```

It first determines the bound value of every active control coefficient:

```text
lower active
    use lower bound

upper active
    use upper bound
```

Those values form the active vector

$$
u_A.
$$

### 10.1 The inactive control block becomes smaller

If some control coordinates remain inactive, the restricted primal layout keeps the
non-control blocks unchanged and replaces the control block by a smaller block of
dimension

$$
|\mathcal I|.
$$

The matching control-stationarity block is reduced in the same way. If every control
coefficient is active, the control block disappears entirely from the restricted
primal/stationarity layouts.

The state and equality-multiplier blocks remain.

### 10.2 The base point with only active values supplies the affine shift

To obtain the shifted right-hand sides, the implementation constructs a base primal
point whose:

```text
non-control/inactive unknowns
    are zero

active control coordinates
    contain their selected bounds
```

It evaluates the **base KKT residual** at that point. Those residual values are
precisely the affine contributions caused by the fixed active variables.

The restricted product negates the appropriate residual blocks to obtain its new
stationarity and equality right-hand sides. This is the implementation form of the
algebra

$$
d-D_A\bar x_A
$$

and

$$
c_I-Q_{IA}\bar x_A.
$$

### 10.3 Actions are restricted without rebuilding the original operator

For a restricted inactive-coordinate primal direction, the subproblem:

1. expands it to the base primal layout;
2. inserts zeros on active perturbation coordinates;
3. applies the base $Q$ or $D$ action;
4. restricts the resulting stationarity covector back to inactive coordinates.

The multiplier action similarly applies the base $D^{\mathsf T}$ and restricts the
stationarity result. The full KKT transpose is obtained through the same
expand–apply–restrict strategy.

So `ActiveSetKKTSubproblemT` is primarily an operator restriction and affine-shift
adapter around the original quadratic KKT product.

## 11. Why active perturbations are zero even though active values are not

This distinction is central. The current active value may be, for example,

$$
u_i=r_i.
$$

That nonzero value contributes to the affine right-hand side as described above. But
during the restricted solve, the active coefficient is **fixed**. Therefore its
variation is

$$
\delta u_i=0.
$$

That is why restricted operator actions expand a free direction using zeros on active
coordinates. The two roles are:

```text
active base value
    contributes to affine shifts

active perturbation
    equals zero in restricted Q / D actions
```

Confusing them would either lose the effect of the bound value or incorrectly allow a
fixed coordinate to move.

## 12. Solving the restricted KKT problem determines state, free control, and equality multiplier

Once the active-set subproblem has been constructed, PDAS calls the supplied KKT solve
action on that product. Conceptually the unknown is

```text
state / other primal blocks
inactive control coordinates
equality multiplier
```

while active control coefficients are already known. The solve returns a point on the
**restricted** layout. `to_base_point()` then reconstructs the full primal vector by
prolonging:

```text
free coordinates
+
stored active bound values
    ↓
full control block
```

and leaves the equality multiplier unchanged. The result is once again a point in the
base KKT product. The PDAS contract test checks both the reduced free-control dimension
and the reconstructed bound value in the full solution.

## 13. After the KKT solve, the missing box multiplier is reconstructed

The restricted KKT solve does not solve explicitly for one box multiplier per active
control coefficient. Instead, after reconstructing the full base point, PDAS evaluates
the base KKT residual.

Let $r_u$ be the control-stationarity block of that residual. On an inactive coordinate, the
restricted KKT problem has imposed

$$
r_{u,i}=0
$$

up to the KKT solve tolerance. The box multiplier should therefore be

$$
\mu_i=0
\qquad
i\in\mathcal I.
$$

On an active coordinate, ordinary control stationarity was removed from the restricted
problem. The multiplier is chosen to restore constrained stationarity:

$$
r_{u,i}+\mu_i=0.
$$

Hence

$$
\mu_i=-r_{u,i}
\qquad
i\in\mathcal A.
$$

That is exactly what the current `make_box_multiplier()` helper does:

```text
inactive index
    box multiplier = 0

active index
    box multiplier = - base control-stationarity residual
```

The resulting object remains a covector on the control layout.

## 14. The new multiplier produces the next active-set guess

After the restricted solve, PDAS now has:

```text
new full primal point
new equality multiplier
new box multiplier
```

The box multiplier is converted to its declared primal representative,

$$
m=G^{-1}\mu,
$$

and the classification rule is applied again. Thus one iteration has the form:

```text
current active-set guess
        │
        ▼
fix active controls at bounds
        │
        ▼
build restricted KKT product
        │
        ▼
solve restricted KKT system
        │
        ▼
reconstruct full primal point
        │
        ▼
recover box multiplier
        │
        ▼
classify again
        │
        ├── changed → repeat
        └── stable → check full residual conditions
```

The classification is therefore not a one-time preprocessing step. It is the nonlinear
part of the algorithm.

## 15. Why active-set stability alone is not enough

Suppose the classification stops changing. That is encouraging, but it does not prove
that the current point satisfies the bound-constrained KKT conditions accurately.

A stable active set could coexist with:

- a poorly solved restricted KKT system;
- a control slightly outside its bounds;
- a multiplier with the wrong sign;
- a nonzero complementarity product;
- a nonzero constrained stationarity residual;
- a nonzero equality residual.

The current solver therefore requires both:

```text
active set stable
+
numerical KKT/complementarity residuals converged
```

before returning `PDASStoppingReason::converged`. This is a useful design choice: the
discrete active labels and the numerical first-order equations provide different
evidence.

## 16. The PDAS report separates five first-order residuals

The current iteration report records several quantities rather than compressing
everything into one norm.

### Primal feasibility

For each coefficient, the violation is measured against the interval:

$$
\max
\left(
\ell_i-u_i,
u_i-r_i
\right).
$$

The report takes the maximum over all coefficients. At a feasible point this quantity is
nonpositive before the outer maximum with the initial zero, so the stored violation is
zero.

### Dual feasibility

Using the primal multiplier representative $m_i$, the required signs are:

```text
at lower bound
    m_i <= 0

at upper bound
    m_i >= 0

strictly inside
    m_i = 0
```

The implementation measures the corresponding positive violation:

```math
\begin{cases}
\max(m_i,0), & u_i\leq\ell_i,\\
\max(-m_i,0), & u_i\geq r_i,\\
|m_i|, & \ell_i<u_i<r_i.
\end{cases}
```

and keeps the maximum over coefficients.

### Complementarity residual

The signed multiplier form uses

$$
\max(m_i,0)(r_i-u_i)
$$

for the upper side and

$$
\min(m_i,0)(u_i-\ell_i)
$$

for the lower side. The maximum absolute value over both sides and all coefficients is
recorded. This quantity is small only when a nonzero multiplier is attached to a nearly
active bound.

### Constrained stationarity residual

The base KKT stationarity residual contains the unconstrained control block $r_u$.

PDAS adds the box multiplier,

$$
r_u+\mu,
$$

and measures the norm of the **whole stationarity residual**, including unchanged
non-control stationarity blocks. So this check asks whether the complete KKT
stationarity equation, augmented by the box multiplier on the constrained control block,
is satisfied.

### Equality residual

The equality residual remains the base KKT equality block

$$
Dx-d.
$$

PDAS measures its norm independently.

### 16.1 Why keep these residuals separate?

Different failures have different interpretations. For example:

```text
small stationarity, large primal violation
    algebraic KKT balance but control outside the box

small primal violation, wrong dual sign
    feasible control but multiplier inconsistent with the active bound

small complementarity, large equality residual
    bound conditions look plausible but the PDE/equality constraint is unsolved
```

A single combined scalar would hide those distinctions.

## 17. The KKT solve report and the PDAS residual report are also distinct

Each active-set iteration solves a quadratic KKT product. That solve has its own report:

```text
linear solve status
stationarity residual
equality residual
```

PDAS then reconstructs the full point and box multiplier and computes the
**bound-constrained** residuals described above. These two levels are not redundant. The
restricted KKT solve report answers:

> Was the equality-constrained subproblem for this active-set hypothesis solved
> successfully?

The PDAS report answers:

> Does the reconstructed full point satisfy the box-constrained first-order
> conditions, and did the classification stabilize?

The solver returns convergence only when the restricted KKT solve converged and the
outer PDAS conditions also pass.

## 18. One complete PDAS iteration

We can now write the algorithm without hiding any of its mathematical roles. Start from:

- a feasible control;
- an equality-constrained KKT point;
- an initial box multiplier;
- a positive classification parameter $c$.

Then:

```text
1. classify
       (u, mu) -> lower / free / upper

2. freeze active controls
       lower-active u_i = ell_i
       upper-active u_i = r_i

3. form restricted KKT product
       inactive control coordinates remain unknown
       active contributions shift affine rhs terms

4. solve restricted KKT problem
       state / inactive control / equality multiplier

5. prolong to full primal coordinates
       reinsert active bound values

6. evaluate base KKT residual

7. recover box multiplier
       inactive: mu_i = 0
       active: mu_i = -r_{u,i}

8. reclassify using the new (u, mu)

9. measure
       primal feasibility
       dual feasibility
       complementarity
       constrained stationarity
       equality residual

10. stop if
       active set is stable
       AND restricted KKT solve converged
       AND all outer residual tolerances pass
```

If the active set changes, the process repeats with the new classification. If the
restricted KKT solve fails, the solver reports `kkt_solve_failed`. If the iteration
limit is reached first, it reports `maximum_iterations`.

## 19. Why PDAS can converge rapidly once the active set is correct

For the quadratic problem considered here, the only nonlinear uncertainty is which
control coefficients belong to which activity class. Suppose an iteration guesses the
final active set correctly.

Then:

- the correct active variables are fixed to their final bounds;
- the correct free variables remain unknown;
- the remaining subproblem is a linear equality-constrained quadratic KKT system.

Solving that restricted system therefore gives the exact solution associated with that
active-set pattern, up to the numerical KKT solve tolerance. If the reconstructed
multipliers confirm the same classification, the active set is stable and no further
combinatorial change is needed.

This explains the characteristic behavior of active-set methods:

```text
early iterations
    identify which inequalities should bind

late iteration
    once classification is right,
    solve the corresponding smooth/equality-constrained problem
```

The difficult part is not repeatedly minimizing a nonlinear quadratic objective. It is
discovering the correct switching pattern.

## 20. The relation to semismooth Newton

The classification formulas can also be viewed as a nonsmooth equation whose pieces
change when a coefficient crosses an activity boundary. For one coefficient, the
conditions switch among:

```text
u_i = ell_i
u_i = r_i
mu_i = 0
```

depending on the signs of the classification expressions. Within one fixed activity
pattern, the equations are smooth and, in the present quadratic case, linear.

Crossing from one pattern to another changes which equation is active. This
piecewise-smooth structure is why PDAS is closely related to semismooth Newton methods
for complementarity/projection equations.

For this manual, however, the active-set interpretation is the more useful one:

```text
classify
restrict
solve KKT
reclassify
```

The implementation exposes exactly those objects.

A future optimization background chapter can derive the semismooth-Newton equivalence in
more detail without making the project concept chapter depend on that theory.

## 21. Initial feasibility is a real precondition of the current solver

`PDASSolverT::solve()` requires the initial control to satisfy the box bounds. It
rejects an infeasible initial control before beginning the active-set loop. This differs
from algorithms that allow infeasible iterates and drive feasibility residuals to zero
from outside the admissible set.

The current PDAS path starts from a feasible box control and keeps active coordinates on
their bounds through the restricted subproblem construction. The report still measures
primal feasibility after every solve because numerical or contract errors should remain
visible rather than being assumed away.

## 22. The current PDAS contract is deliberately narrower than general inequality KKT

The mathematical KKT theory for inequalities can represent:

- many inequality families;
- coupled inequalities;
- nonlinear constraints;
- complementarity functions;
- separate lower/upper multiplier spaces.

The current contract implements a specific useful slice:

```text
one control block
componentwise finite lower/upper bounds
a declared dual-to-primal multiplier representation
quadratic equality-constrained KKT base product
active-set restriction of that control block
```

`ActiveSetSelectionT` explicitly supports exactly one control block, and the PDAS solver
is constructed with the index of that block. This narrowness is useful to keep visible.

It is an implemented formulation, not a claim that every inequality-constrained PDE
problem reduces to the same contract.

## 23. `BoxComplementarityT` is the bridge between bounds and active-set logic

The complementarity object packages:

```text
BoxBoundsT
BoxMultiplierRepresentationT
optional box-data token
```

and exposes the two operations PDAS needs:

```text
multiplier_to_primal
classify
```

It also exposes the reverse primal-to-dual conversion, which is useful when a primal
multiplier representation needs to be turned back into the covector form used by
stationarity.

The object does not solve a KKT problem. It supplies the bound geometry and multiplier
representation needed to decide which coordinates should be active.

### 23.1 Classification and KKT restriction are separate responsibilities

This separation is important. `BoxComplementarityT` answers:

> Given the current control and multiplier, which coordinates look lower-active,
> free, or upper-active?

`ActiveSetKKTSubproblemT` answers:

> Given that classification, what equality-constrained KKT product should be solved?

`PDASSolverT` orchestrates the iteration between them. The three levels are:

```text
complementarity
    classify

active-set KKT adapter
    restrict and shift

PDAS solver
    iterate until classification + residuals converge
```

That is a more precise description than treating PDAS as one monolithic solver routine.

## 24. The restricted product needs its own structural assumptions

The base quadratic KKT product carries declarations such as:

- equality rank condition;
- positivity on the equality kernel;
- transpose consistency;
- symmetry.

Fixing active control coordinates changes the primal space and can change the effective
equality/KKT operator. The current `PDASPolicy` therefore carries
`active_set_assumptions` for the restricted products.

This is not automatically copied as a theorem from the base product. The producer
declares that the active-set subproblems satisfy the assumptions required by the KKT
solver.

For example, an active-set policy may declare that:

```text
base equality rows
+
the current coordinate restrictions

retain the required rank structure
```

and that the quadratic objective remains positive on the restricted equality kernel. The
contract records those declarations; it does not prove them analytically.

## 25. PDAS inherits the solver compatibility of each restricted KKT product

Each active-set subproblem is itself an `EqualityConstrainedQuadraticKKTProductT`. It
retains the base product's symmetry category and rebuilds the primal/stationarity
pairing on the restricted free-coordinate layouts.

Therefore the supplied KKT solve action can use the same solver-policy machinery
developed in the previous chapter, provided the restricted product satisfies the
necessary assumptions.

Conceptually:

```text
base KKT product
    symmetric-indefinite, for example
        │
        ▼
active-set restriction
        │
        ▼
restricted KKT product
    with inactive control coordinates
        │
        ▼
MINRES/GMRES-compatible solve action
    according to the restricted product
```

PDAS itself does not implement a new Krylov method. It repeatedly supplies a changing
KKT product to an existing KKT solve capability.

## 26. The test problem shows the active-set restriction concretely

The backend-neutral PDAS contract test uses a primal layout with:

```text
state block      dimension 2
control block    dimension 2
```

and a state equality relation whose action is conceptually

$$
y-u.
$$

One test selects:

```text
control coefficient 0
    upper active

control coefficient 1
    free
```

The active-set subproblem therefore keeps:

```text
state block      dimension 2
inactive control dimension 1
```

rather than the original two-dimensional control block. The first control coefficient is
fixed at its upper bound. The test then checks that:

- the restricted equality right-hand side contains the contribution of that fixed
  active value;
- the restricted control-stationarity right-hand side contains the corresponding
  affine shift;
- solving the restricted KKT system and prolonging it returns the active control
  exactly at the bound;
- the inactive control coordinate is solved rather than fixed.

This is a compact executable example of Sections 9–12.

## 27. Complementarity tests isolate the classification contract

The complementarity contract test uses a four-dimensional control with a diagonal metric
representation. It constructs a primal control and a multiplier representative such that
the expected classification is:

```text
upper
lower
free
free
```

The test verifies:

- bound feasibility;
- metric dual-to-primal round-trip behavior;
- the lower/free/upper labels;
- active/free index extraction;
- restriction of full vectors;
- prolongation back to the original ordering;
- rejection of invalid classification parameters and nonfinite values;
- retention and validation of the metric realization witness.

These tests are intentionally independent of the KKT solver. That separation makes it
easier to diagnose whether a failure comes from:

```text
classification / multiplier representation

or

restricted KKT construction / solve
```

## 28. PDAS reports the active-set history explicitly

Each `PDASIterationReportT` records:

```text
iteration index
selection used for the current subproblem
number of activity changes
whether the next selection is identical
primal violation
dual violation
complementarity residual
stationarity residual
equality residual
feasibility/convergence flags
underlying KKT solve report
```

The complete `PDASSolveResultT` then retains:

```text
final base KKT point
final box multiplier
all iteration reports
stopping reason
```

The stopping reasons are currently:

```text
converged
maximum_iterations
kkt_solve_failed
```

This gives the caller enough evidence to distinguish:

- a successful stable active set;
- an iteration cap;
- failure of the linear/KKT subproblem solve.

The report structure mirrors the actual algorithm rather than reducing every outcome to
a boolean.

## 29. Projection and PDAS solve the same first-order problem in different ways

It is worth returning once more to the connection with Part I. For a box-constrained
reduced problem, projected steepest descent uses

$$
\widehat u
=
P_C^G(u-\alpha g).
$$

The projection enforces feasibility after each trial update. PDAS instead carries a box
multiplier and repeatedly solves equality-constrained subproblems with a guessed active
set.

The two routes can be summarized as:

```text
projection method
    derivative
        ↓
    metric gradient
        ↓
    unconstrained trial
        ↓
    project onto box
        ↓
    repeat

PDAS
    KKT stationarity + box multiplier
        ↓
    classify bounds
        ↓
    fix active coordinates
        ↓
    solve restricted KKT problem
        ↓
    update multiplier / classification
        ↓
    repeat
```

Both target the same first-order box-constrained conditions. They exploit them
differently.

### 29.1 Why PDAS belongs beside KKT rather than beside line search

The current PDAS implementation is built around a quadratic KKT product, not around a
reduced objective/gradient callback. Its inner operation is:

```text
solve an equality-constrained KKT system
```

rather than:

```text
take a projected search step and globalize it
```

So although both methods handle the same type of box constraint, their natural software
homes are different. This is why the manual treats projected reduced optimization in the
reduced solver chapter and PDAS in the all-at-once/KKT part.

## 30. A compact mathematical correspondence

| Bound-constrained idea | Mathematical form | Current object |
| --- | --- | --- |
| finite bounds | $\ell\leq u\leq r$ | `BoxBoundsT` |
| box multiplier | $\mu\in U_h^{\ast}$ | multiplier covector |
| primal multiplier representative | $m=G^{-1}\mu$ | `BoxMultiplierRepresentationT` |
| constrained stationarity | $r_u+\mu=0$ | PDAS stationarity check |
| lower activity | $m_i+c(u_i-\ell_i)<0$ | `BoxActivity::lower` |
| upper activity | $m_i+c(u_i-r_i)>0$ | `BoxActivity::upper` |
| free activity | neither active inequality | `BoxActivity::inactive` |
| inactive/active index sets | $\mathcal I,\mathcal A$ | `ActiveSetSelectionT` |
| fixed active values | $u_i=\ell_i$ or $r_i$ | `active_values()` |
| restricted KKT system | solve only inactive primal coordinates | `ActiveSetKKTSubproblemT` |
| active multiplier recovery | $\mu_i=-r_{u,i}$ | `make_box_multiplier` |
| active-set iteration | classify → restrict → solve → reclassify | `PDASSolverT` |

## 31. Following the complementarity/PDAS path through the source

A useful reading order is:

### Bounds, multiplier representation, and classification

- [`include/nmopt/contract/complementarity.hpp`](../../../include/nmopt/contract/complementarity.hpp)

Read:

```text
BoxActivity
BoxMultiplierRepresentationT
BoxBoundsT
ActiveSetSelectionT
BoxComplementarityT
```

The `classify()` method is much easier to understand after Sections 1–6 of this chapter.

### Restricted KKT product

Then read the first part of:

- [`include/nmopt/contract/pdas.hpp`](../../../include/nmopt/contract/pdas.hpp)

Focus on `ActiveSetKKTSubproblemT`. The key operations are:

```text
make_active_values
initialize_restriction_layouts
expand_primal
restrict_stationarity
make_active_product
to_base_point
```

These implement the fixed-active/inactive-coordinate algebra developed above.

### PDAS iteration

Continue in the same file with:

```text
PDASPolicy
PDASIterationReportT
PDASSolveResultT
PDASSolverT
```

Then read `solve()` as the ten-step iteration from Section 18. Finally inspect:

```text
make_box_multiplier
make_report
```

to see the multiplier recovery and residual definitions.

### Contract tests

Compare:

- [`tests/contract/complementarity_contract.cc`](../../../tests/contract/complementarity_contract.cc)
- [`tests/contract/pdas_contract.cc`](../../../tests/contract/pdas_contract.cc)

The first isolates the bound/multiplier/classification layer. The second exercises
active KKT restriction, reconstruction, PDAS convergence, and failure/stopping behavior.

## 32. Part II now contains two distinct formulation families

The formulation/optimization part of the manual has now developed two major ways of
organizing the same PDE-constrained first-order information.

### Reduced family

```text
state eliminated by solve
        ↓
reduced derivative
        ↓
line search / trust region
        ↓
optional metric projection for bounds
```

### All-at-once family

```text
state + adjoint + control retained
        ↓
optimality residual
        ↓
quadratic KKT specialization
        ↓
optional complementarity / PDAS for bounds
```

Neither should be treated as a secondary appendix to the other. They use shared
numerical language—spaces, covectors, transpose actions, metrics, solve evidence—but
organize the optimization problem differently.

This distinction will matter in Part III: a semantic/compiler request does not merely
ask for "an optimizer". It selects a formulation product whose required numerical
realization depends on the chosen route.

## 33. What comes next

The next part of the concept manual changes perspective. So far, we have assumed that
executable residuals, objectives, metrics, KKT products, and solve services already
exist.

Part III asks:

> How does a problem description become those executable numerical objects?

The first chapter in that part, **Semantic problem model**, will begin again from a
concrete PDE-constrained problem and construct its backend-neutral `ProblemSpec`
progressively.

After that, **Validation, resolution, and capabilities** will explain what can be
checked before numerical realization, and **Compilation and lowering** will follow one
resolved problem into the registered deal.II realization path.

## Read later

Useful existing documents are:

- [Theoretical formalism](../../design/theoretical-formalism.md), for the project's
  variational-inequality, multiplier, and first-order conventions.
- [Optimality systems and KKT](optimality-systems-and-kkt.md), for the equality KKT
  product on which the current PDAS subproblems are built.
- [Metrics, gradients, and constraints](metrics-gradients-and-constraints.md), for
  metric projection and the primal–dual identification reused by multiplier
  representations.
- [v1 semantic/compiler capability](../../implementation/v1/semantic-compiler.md),
  for the exact currently registered complementarity/PDAS capability ledger.

A future optimization background chapter can go deeper into normal cones,
complementarity functions, semismooth Newton theory, and convergence results without
turning this project concept chapter into a general nonlinear-optimization textbook.
