# Optimality systems and KKT

## What changes when the state is no longer eliminated

The reduced formulation treated the control as the optimization variable and hid the
state equation inside the control-to-state map

$$
y=S_{h}(u).
$$

That led to a reduced objective

$$
j_{h}(u)
:=
J_{h}(S_{h}(u),u)
$$

and an adjoint construction for its derivative. There is another way to use exactly the
same first-order mathematics.

Instead of solving the state equation first and then optimizing over $u$, retain the
state, control, and adjoint as simultaneous unknowns and solve the coupled optimality
conditions directly.

For a state/control problem, the unknown becomes

$$
(y,p,u)
\in
Y_{h}\times Z_{h}\times U_{h}.
$$

The equations are no longer just the PDE residual. They contain three roles:

```text
state equation
    physical feasibility

adjoint equation
    stationarity with respect to the state

control stationarity
    stationarity with respect to the decision variable
```

This is the basic all-at-once viewpoint.

For a nonlinear problem these equations form a nonlinear residual system. For a linear
PDE with a quadratic objective, the same first-order system becomes a linear
saddle-point problem of KKT type.

This chapter develops that progression:

```text
constrained optimization problem
        ↓
Lagrangian
        ↓
first-order optimality system
        ↓
all-at-once residual in (state, adjoint, control)
        ↓
linear-quadratic specialization
        ↓
canonical equality-constrained KKT system
        ↓
Q / D / D^T actions, pairings, assumptions, solver policy
```

Only after that derivation will we map the objects to `SuppliedOTDSystemT` and
`EqualityConstrainedQuadraticKKTProductT`. Bounds and complementarity are deliberately
deferred. They add another layer of first-order conditions and lead naturally to PDAS in
Chapter 8.

## 1. Return to the Lagrangian

Consider the discrete equality-constrained problem

```math
\begin{aligned}
\min_{y,u}\quad & J_{h}(y,u),\\
\text{subject to}\quad & E_{h}(y,u)=0
\quad\text{in }Z_{h}^{\ast}.
\end{aligned}
```

The project uses the Lagrangian convention

$$
\mathcal L_{h}(y,u,p)
:=
J_{h}(y,u)
-
\langle E_{h}(y,u),p\rangle_{Z_{h}^{\ast},Z_{h}}.
$$

The adjoint variable is

$$
p\in Z_{h}.
$$

The minus sign is the same convention used in the reduced state–adjoint chapter. The
first-order conditions come from requiring the derivative of $`\mathcal L_{h}`$ to vanish
with respect to each unknown role.

### 1.1 Variation with respect to the adjoint gives the state equation

Perturb $p$ by

$$
\delta p\in Z_{h}.
$$

Then

```math
D_{p}\mathcal L_{h}(y,u,p)[\delta p]
=
-
\langle E_{h}(y,u),\delta p\rangle.
```

For this to vanish for every $\delta p$,

$$
E_{h}(y,u)=0
\quad\text{in }Z_{h}^{\ast}.
$$

So the multiplier/adjoint variation simply recovers primal feasibility.

### 1.2 Variation with respect to the state gives the adjoint equation

Now perturb the state by

$$
\delta y\in Y_{h}.
$$

Then

```math
D_{y}\mathcal L_{h}(y,u,p)[\delta y]
=
D_{y}J_{h}(y,u)[\delta y]
-
\left\langle
D_{y}E_{h}(y,u)[\delta y],
p
\right\rangle.
```

Use the transpose identity

$$
\left\langle
D_{y}E_{h}(y,u)[\delta y],
p
\right\rangle
=
\left\langle
D_{y}E_{h}(y,u)^{\ast}p,
\delta y
\right\rangle.
$$

Therefore

```math
D_{y}\mathcal L_{h}(y,u,p)[\delta y]
=
\left\langle
D_{y}J_{h}(y,u)
-
D_{y}E_{h}(y,u)^{\ast}p,
\delta y
\right\rangle.
```

Stationarity for every $\delta y$ gives

$$
D_{y}J_{h}(y,u)
-
D_{y}E_{h}(y,u)^{\ast}p
=
0
\quad\text{in }Y_{h}^{\ast}.
$$

Equivalently,

$$
D_{y}E_{h}(y,u)^{\ast}p
=
D_{y}J_{h}(y,u).
$$

This is exactly the adjoint equation that appeared in the reduced formulation.

### 1.3 Variation with respect to the control gives control stationarity

Finally, perturb

$$
u
\mapsto
u+\varepsilon\delta u.
$$

The control variation is

```math
D_{u}\mathcal L_{h}(y,u,p)[\delta u]
=
\left\langle
D_{u}J_{h}(y,u)
-
D_{u}E_{h}(y,u)^{\ast}p,
\delta u
\right\rangle.
```

For an unconstrained control, first-order stationarity requires

$$
D_{u}J_{h}(y,u)
-
D_{u}E_{h}(y,u)^{\ast}p
=
0
\quad\text{in }U_{h}^{\ast}.
$$

The same covector was the reduced derivative:

$$
j_{h}'(u)
=
D_{u}J_{h}
-
D_{u}E_{h}^{\ast}p.
$$

The difference is not in the equation. The difference is in **how it is used**.

## 2. Reduced and all-at-once methods use the same conditions differently

The reduced route organizes the first-order conditions as a sequence of nested
operations:

```text
given u
    ↓
solve state equation for y
    ↓
solve adjoint equation for p
    ↓
form reduced derivative
    ↓
optimization update for u
```

The all-at-once route instead regards all three equations as one residual system:

```text
unknown (y, p, u)
        │
        ▼
┌──────────────────────────────┐
│ state equation               │
│ adjoint equation             │
│ control stationarity         │
└──────────────────────────────┘
        │
        ▼
solve all blocks together
```

Define the optimality residual

$$
F_{h}(y,p,u)
:=
\begin{bmatrix}
E_{h}(y,u)\\
D_{y}J_{h}(y,u)-D_{y}E_{h}(y,u)^{\ast}p\\
D_{u}J_{h}(y,u)-D_{u}E_{h}(y,u)^{\ast}p
\end{bmatrix}.
$$

Its codomain is the product of the three equation-dual spaces. An all-at-once solution
satisfies

$$
F_{h}(y,p,u)=0.
$$

### 2.1 The state need not be feasible at every intermediate iterate

This is one of the most important algorithmic differences. In a reduced method, every
accepted control is paired with a state satisfying

$$
E_{h}(y,u)=0
$$

to the accuracy of the supplied state solve. An all-at-once nonlinear iteration may
instead visit intermediate points for which

$$
E_{h}(y,u)\neq0.
$$

The state equation is one residual block being driven to zero together with the adjoint
and stationarity blocks.

This can be useful when a coupled Newton or Krylov method has a good block
preconditioner and one does not want to perform a fully converged PDE solve inside every
outer optimization iteration.

It also makes the linear algebra larger and more strongly coupled. Neither organization
is automatically preferable.

### 2.2 What is gained and what is traded

A reduced method has a smaller optimization variable space, but each reduced evaluation
hides state and adjoint solves. An all-at-once method exposes a larger coupled system,
but that exposure creates opportunities for:

- coupled Newton steps;
- block preconditioners;
- simultaneous treatment of state, adjoint, and decision variables;
- Krylov methods that see interactions among all first-order equations directly.

The trade is that the solver must now understand and control the whole saddle-point
structure rather than relying on a nested state-solve abstraction.

This is why reduced and all-at-once formulations are distinct products in the repository
rather than two names for the same runtime object.

## 3. A nonlinear optimality system is itself an executable residual

Once

$$
F_{h}:
Y_{h}\times Z_{h}\times U_{h}
\longrightarrow
R_{h}^{\ast}
$$

has been defined, where $`R_{h}`$ denotes the product test space for the optimality
equations, the same first-order vocabulary from Part I applies again.

At a point

$$
w
:=
(y,p,u),
$$

the linearization is

$$
F_{h}'(w)[\delta w].
$$

Its transpose action is

$$
F_{h}'(w)^{\ast}q.
$$

A Newton-type all-at-once method asks for a correction

$$
\delta w
$$

satisfying

$$
F_{h}'(w)[\delta w]
=
-F_{h}(w).
$$

So the JVP of the optimality system is the operator used by a matrix-free Newton–Krylov
solve.

There is an important derivative-order consequence. The residual $`F_{h}`$ already contains
first derivatives of $`J_{h}`$ and $`E_{h}`$, so differentiating $`F_{h}`$ generally introduces
second derivatives of the original optimization problem. For example, linearizing

$$
D_{y}E_{h}(y,u)^{\ast}p
$$

with respect to $y$ or $u$ differentiates a first derivative of the original PDE
residual.

A supplied OTD system can provide this optimality JVP directly; it need not be
reconstructed from a first-order DTO model whose interface stops at $`J_{h}'`$ and $`E_{h}'`$.

The VJP supplies the transpose action needed for transpose consistency tests,
diagnostics, or algorithms that require the transposed operator. The key point is that
this JVP/VJP pair belongs to the **optimality residual** $`F_{h}`$, not to the original PDE
residual $`E_{h}`$.

## 4. DTO and OTD differ in where discretization occurs

The reduced chapters worked with an already-discrete objective and residual:

$$
J_{h},
\qquad
E_{h}.
$$

The adjoint and reduced derivative were then obtained by differentiating those discrete
maps. That is a discretize-then-optimize, or DTO, construction. An
optimize-then-discretize, or OTD, route changes the order:

```text
continuous PDE-constrained problem
        ↓
derive continuous first-order optimality system
        ↓
choose spaces/discretization for that system
        ↓
discrete state/adjoint/stationarity equations
```

The two routes may produce closely related systems, and for carefully compatible
discretizations they may even coincide algebraically. They are not identical by
definition.

Choices made during discretization—quadrature, stabilization, boundary treatment,
Petrov–Galerkin testing, transformations, or discrete objective realization—can make
"differentiate the discrete system" differ from "discretize the differentiated
continuous system".

### 4.1 Why provenance matters

If an application supplies an OTD system, the framework should not silently reinterpret
it as though it had been generated by differentiating a DTO `ExecutableModelT`.

Its equations may have been derived and discretized independently. That is why the
project has a **supplied OTD** formulation product. The product says, in effect:

> Here is the discrete optimality system I intend you to solve, together with its
> residual, derivative actions, block roles, and solve operation.

It preserves that formulation provenance rather than reconstructing a different
first-order system from a DTO model.

## 5. `SuppliedOTDSystemT` represents one selected three-block system

The current supplied-OTD contract is deliberately narrower than the abstract all-at-once
idea. It expects three variable blocks:

```text
state
adjoint
control
```

and three residual blocks:

```text
state equation
adjoint equation
control stationarity
```

`SuppliedOTDBlockSelection` records which block index plays each role. The default
selected form can be pictured as:

```text
variable point
[ state | adjoint | control ]
        │
        ▼
supplied OTD residual
        │
        ▼
[ state equation | adjoint equation | control stationarity ]
```

The space IDs themselves remain producer-owned. The contract records the block roles and
checks that the selected indices are distinct and within their layouts.

### 5.1 The supplied system exposes residual, JVP, VJP, and solve

`SuppliedOTDSystemT` contains four executable operations:

```text
residual(point)

residual_jvp(point, tangent)

residual_vjp(point, residual_seed)

solve(initial_point)
```

The first three are the same value/linearization/transpose pattern developed in
[Operators, derivatives, and adjoints](03-operators-derivatives-and-adjoints.md), now
applied to the whole first-order system.

The fourth operation is formulation-specific: it solves for the complete
state/adjoint/control point rather than just applying the linearized operator. The
result includes the same backend-neutral solve-report shape used elsewhere in the
contract layer.

### 5.2 Block accessors expose equation roles without redefining the residual

The supplied system also provides convenience views of the residual blocks:

```text
state_residual(point)
adjoint_residual(point)
control_stationarity(point)
```

These are selections from the complete residual. They do not define three independent
models. The complete system remains the object whose JVP, VJP, and solve actions must be
coherent.

## 6. The running elliptic-control problem produces a linear all-at-once system

Return to the discrete linear-quadratic problem from the first chapter. The state
equation is

$$
Ay-Bu-f=0.
$$

The objective is

```math
J_{h}(y,u)
:=
\frac12
\left(
y^{\mathsf T}M_{y} y
-
2q^{\mathsf T}y
+
c
\right)
+
\frac{\beta}{2}
u^{\mathsf T}N_{u} u.
```

The project Lagrangian is

$$
\mathcal L_{h}(y,u,p)
:=
J_{h}(y,u)
-
p^{\mathsf T}(Ay-Bu-f).
$$

Differentiate with respect to each variable.

### 6.1 State equation

Variation with respect to $p$ gives

$$
Ay-Bu-f=0.
$$

### 6.2 Adjoint equation

The state derivative of the objective is

$$
D_{y}J_{h}
=
M_{y} y-q.
$$

The state derivative of the residual is $A$. Therefore

$$
M_{y} y-q-A^{\mathsf T}p=0.
$$

Equivalently,

$$
A^{\mathsf T}p
=
M_{y} y-q.
$$

### 6.3 Control stationarity

The direct control derivative is

$$
D_{u}J_{h}
=
\beta N_{u} u.
$$

Because

$$
D_{u}E_{h}=-B,
$$

the transpose contribution is

$$
D_{u}E_{h}^{\ast}p
=
-B^{\mathsf T}p.
$$

Hence control stationarity is

$$
\beta N_{u} u+B^{\mathsf T}p=0.
$$

### 6.4 The three equations can be solved together

Ordering the unknowns as

$$
(y,p,u),
$$

the optimality system is

```math
\begin{bmatrix}
A & 0 & -B\\
M_{y} & -A^{\mathsf T} & 0\\
0 & B^{\mathsf T} & \beta N_{u}
\end{bmatrix}
\begin{bmatrix}
y\\
p\\
u
\end{bmatrix}
=
\begin{bmatrix}
f\\
q\\
0
\end{bmatrix}.
```

This is a linear all-at-once system.

It contains the same state and adjoint equations used by the reduced formulation, but
they are now coupled to the stationarity equation and solved simultaneously.

The matrix is not yet written in the standard two-by-two KKT form. To reach that form,
combine state and control into one primal variable and account for the project's adjoint
sign convention.

## 7. The canonical quadratic problem combines state and control into one primal variable

Define

$$
x
:=
\begin{bmatrix}
y\\
u
\end{bmatrix}.
$$

The quadratic objective can be written

$$
\varphi(x)
:=
\frac12
x^{\mathsf T}Qx
-
c_{x}^{\mathsf T}x
+
c_{0},
$$

with

```math
Q
:=
\begin{bmatrix}
M_{y} & 0\\
0 & \beta N_{u}
\end{bmatrix},
\qquad
c_{x}
:=
\begin{bmatrix}
q\\
0
\end{bmatrix}.
```

The state equation becomes one linear equality constraint

$$
Dx=d,
$$

where

```math
D
:=
\begin{bmatrix}
A & -B
\end{bmatrix},
\qquad
d:=f.
```

The optimization problem is therefore

```math
\begin{aligned}
\min_{x}\quad&
\frac12 x^{\mathsf T}Qx-c_{x}^{\mathsf T}x+c_{0},
\\
\text{subject to}\quad&
Dx=d.
\end{aligned}
```

This is an equality-constrained quadratic program.

### 7.1 The KKT multiplier uses the opposite sign from the framework adjoint

For the canonical quadratic program it is conventional to use the Lagrangian

$$
\widehat{\mathcal L}(x,\lambda)
:=
\varphi(x)
+
\lambda^{\mathsf T}(Dx-d).
$$

The project PDE convention instead used

$$
\mathcal L_{h}
=
J_{h}-p^{\mathsf T}E_{h}.
$$

Since the equality residual is

$$
E_{h}=Dx-d,
$$

the two conventions agree when

$$
\boxed{
\lambda=-p.
}
$$

This sign conversion is not cosmetic.

It explains why the quadratic KKT product keeps a distinct multiplier role and an
explicit multiplier-to-adjoint conversion rather than simply renaming the adjoint block.

## 8. Derive the canonical KKT equations

Differentiate the canonical quadratic Lagrangian

$$
\widehat{\mathcal L}(x,\lambda)
=
\frac12x^{\mathsf T}Qx
-
c_{x}^{\mathsf T}x
+
\lambda^{\mathsf T}(Dx-d)
+
c_{0}.
$$

Assume $Q$ is symmetric under the declared primal/stationarity pairing. Variation with
respect to $x$ gives

$$
Qx-c_{x}+D^{\mathsf T}\lambda=0.
$$

Variation with respect to $\lambda$ gives

$$
Dx-d=0.
$$

Together,

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
c_{x}\\
d
\end{bmatrix}.
```

This is the canonical equality-constrained KKT system.

### 8.1 Recover the state/adjoint/control equations

Insert

$$
x=
\begin{bmatrix}
y\\
u
\end{bmatrix},
\qquad
\lambda=-p.
$$

The stationarity equation becomes

```math
\begin{bmatrix}
M_{y} & 0\\
0 & \beta N_{u}
\end{bmatrix}
\begin{bmatrix}
y\\
u
\end{bmatrix}
-
\begin{bmatrix}
q\\
0
\end{bmatrix}
+
\begin{bmatrix}
A^{\mathsf T}\\
-B^{\mathsf T}
\end{bmatrix}
(-p)
=
0.
```

So its state component is

$$
M_{y} y-q-A^{\mathsf T}p=0,
$$

and its control component is

$$
\beta N_{u} u+B^{\mathsf T}p=0.
$$

The equality equation is

$$
Ay-Bu=f.
$$

Thus the canonical two-block KKT system and the three-block state/adjoint/control system
encode the same linear-quadratic first-order conditions, provided the multiplier
conversion

$$
\lambda=-p
$$

is respected.

## 9. Why the KKT matrix has saddle-point structure

The canonical KKT operator is

$$
K
:=
\begin{bmatrix}
Q & D^{\mathsf T}\\
D & 0
\end{bmatrix}.
$$

The zero block in the multiplier–multiplier position is not an omission. The objective
has curvature in the primal variable $x$ but not in the Lagrange multiplier $\lambda$.

The multiplier exists to enforce the equality constraint. This creates a saddle-point
system rather than an ordinary positive-definite Hessian system.

### 9.1 Why the matrix is not positive definite

For a nonzero multiplier vector $\lambda$,

$$
\begin{bmatrix}
0\\
\lambda
\end{bmatrix}^{\mathsf T}
K
\begin{bmatrix}
0\\
\lambda
\end{bmatrix}
=
0.
$$

So the KKT matrix cannot be positive definite on the full primal/multiplier product
space. Under the standard regularity and positivity assumptions discussed below, the
system can nevertheless be nonsingular.

For the symmetric case it is then a typical **symmetric-indefinite** linear system. That
distinction matters because solver choices such as CG, MINRES, and GMRES have different
structural requirements.

## 10. Full row rank prevents redundant equality directions

The equality operator is

$$
D:
X_{h}
\longrightarrow
Y_{h}^{\ast},
$$

where $`X_{h}`$ is the KKT primal space and $`Y_{h}^{\ast}`$ is the equality-residual space. A
full-row-rank condition means, in finite-dimensional matrix language,

$$
\mathrm{rank}(D)
=
\dim Y_{h}^{\ast}.
$$

Equivalently,

$$
D^{\mathsf T}\lambda=0
\quad\Longrightarrow\quad
\lambda=0.
$$

So the transpose is injective on multiplier coordinates. Why is this useful?

If equality rows are redundant, the same primal feasible set can be represented by
multiple dependent constraints. The associated multiplier is then not uniquely
determined.

Full row rank rules out that redundancy at the KKT boundary.

The contract does not prove the rank condition from matrix entries. It requires the
producer to declare the condition and its policy/evidence.

## 11. Positivity is needed only on feasible directions

The quadratic Hessian $Q$ does not have to be positive definite on the entire primal
space for the equality-constrained problem to have a well-behaved KKT system.

The relevant directions are those that preserve the equality constraint to first order.
A primal perturbation $z$ is tangent to the feasible affine space when

$$
Dz=0.
$$

Such directions belong to

$$
\ker(D).
$$

The important curvature condition is therefore

$$
\langle Qz,z\rangle
>
0
\qquad
\text{for every nonzero }z\in\ker(D).
$$

This is **positivity on the kernel of the equality operator**.

### 11.1 Why kernel positivity is the right condition

Suppose

$$
K
\begin{bmatrix}
x\\
\lambda
\end{bmatrix}
=
0.
$$

Then

```math
\begin{aligned}
Qx+D^{\mathsf T}\lambda&=0,\\
Dx&=0.
\end{aligned}
```

Pair the first equation with $x$:

```math
\langle Qx,x\rangle
+
\langle D^{\mathsf T}\lambda,x\rangle
=
0.
```

Use transpose consistency:

$$
\langle D^{\mathsf T}\lambda,x\rangle
=
\langle \lambda,Dx\rangle.
$$

But

$$
Dx=0.
$$

Therefore

$$
\langle Qx,x\rangle=0.
$$

Since $x\in\ker(D)$, positivity on the kernel implies

$$
x=0.
$$

The first KKT equation then reduces to

$$
D^{\mathsf T}\lambda=0.
$$

Full row rank of $D$ makes $D^{\mathsf T}$ injective, so

$$
\lambda=0.
$$

Thus the homogeneous KKT system has only the trivial solution. This short argument
explains why the contract records **both** a rank condition and a kernel-positivity
condition.

They protect different parts of the saddle-point problem.

## 12. In the contract, primal and stationarity spaces are paired explicitly

The matrix notation

$$
Qx+D^{\mathsf T}\lambda
$$

can make it look as though all vectors simply live in one Euclidean coordinate space.
The contract keeps more structure. The quadratic KKT product carries five layouts:

```text
primal
multiplier
adjoint
stationarity
equality
```

The key maps are

```math
Q:
X_{h}
\longrightarrow
\Sigma_{h}^{\ast},
```

```math
D:
X_{h}
\longrightarrow
Y_{h}^{\ast},
```

and

```math
D^{\mathsf T}:
\Lambda_{h}
\longrightarrow
\Sigma_{h}^{\ast}.
```

Here:

- $`X_{h}`$ is the primal space;
- $`\Sigma_{h}^{\ast}`$ is the stationarity-covector space;
- $`Y_{h}^{\ast}`$ is the equality-residual space;
- $`\Lambda_{h}`$ is the multiplier space paired with the equality space.

In the simplest matrix example, stationarity coordinates have the same dimension as the
primal coordinates and multiplier coordinates have the same dimension as the equality
residuals.

The contract does not reduce those relationships to "same length". It records explicit
block pairings.

### 12.1 Why pairings are separate from layouts

Suppose the primal has two blocks:

```text
state
control
```

while the stationarity residual has two corresponding blocks:

```text
state stationarity
control stationarity
```

The two layouts may use different `SpaceId` values because they represent primal and
dual roles. A pairing record says which primal block is paired with which stationarity
block.

Likewise, a multiplier block is paired with an equality-residual block.
`QuadraticKKTBlockPairing` stores:

```text
pairing id
domain block map
range block map
pairing ids for each matched pair
```

The constructor requires the maps to:

- cover every block;
- be one-to-one;
- match block dimensions;
- use nonempty unique pairing identifiers.

This is runtime metadata for a mathematical relationship that ordinary vector lengths
cannot express.

## 13. The quadratic product exposes actions rather than a monolithic matrix

`EqualityConstrainedQuadraticKKTProductT` does not require one assembled KKT matrix. It
asks for several actions.

### Quadratic action

$$
x
\longmapsto
Qx.
$$

Runtime method:

```text
apply_q(primal)
```

### Equality action

$$
x
\longmapsto
Dx.
$$

Runtime method:

```text
apply_d(primal)
```

### Multiplier/transpose action

$$
\lambda
\longmapsto
D^{\mathsf T}\lambda.
$$

Runtime method:

```text
apply_d_transpose(multiplier)
```

### KKT action

Given

$$
(x,\lambda),
$$

the product forms

```math
\begin{bmatrix}
Qx+D^{\mathsf T}\lambda\\
Dx
\end{bmatrix}.
```

Runtime method:

```text
apply_kkt(point)
```

### KKT residual

Subtract the stored right-hand sides:

```math
\begin{bmatrix}
Qx+D^{\mathsf T}\lambda-c_{x}\\
Dx-d
\end{bmatrix}.
```

Runtime method:

```text
residual(point)
```

This decomposition allows the producer to realize $Q$, $D$, and $D^{\mathsf T}$ with
native sparse matrices, matrix-free operators, or other callback services.

The KKT abstraction describes the block mathematics without prescribing storage.

## 14. The whole KKT transpose is also an explicit action

The product also exposes

```text
apply_kkt_transpose(seed)
```

rather than assuming that symmetry makes a transpose method unnecessary. Let the seed
contain a stationarity-space primal vector $s$ and an equality-space primal vector $q$.

The KKT action satisfies the pairing identity

```math
\left\langle
K
\begin{bmatrix}
x\\
\lambda
\end{bmatrix},
\begin{bmatrix}
s\\
q
\end{bmatrix}
\right\rangle
=
\left\langle
K^{\mathsf T}
\begin{bmatrix}
s\\
q
\end{bmatrix},
\begin{bmatrix}
x\\
\lambda
\end{bmatrix}
\right\rangle.
```

In ordinary compatible matrix coordinates,

```math
K^{\mathsf T}
=
\begin{bmatrix}
Q^{\mathsf T} & D^{\mathsf T}\\
D & 0
\end{bmatrix}.
```

If $Q$ is symmetric and the supplied $D^{\mathsf T}$ action really is the transpose of
$D$ under the declared pairings, then

$$
K^{\mathsf T}=K.
$$

But the contract does not infer this merely from class names. The producer supplies a
transpose action and declares the relevant consistency evidence.

The backend-neutral KKT contract test checks the pairing identity numerically on a small
dense system.

## 15. Symmetry is an executable property with solver consequences

The product records one of two symmetry categories:

```text
symmetric_indefinite
nonsymmetric
```

A symmetric-indefinite KKT declaration requires more than the dimensions to match. The
current constructor requires:

- complete primal/stationarity pairings;
- complete multiplier/equality pairings;
- declared $D$-transpose consistency;
- declared full-KKT transpose consistency;
- a nonempty policy describing that consistency.

This is deliberate.

A KKT-shaped block matrix is not automatically symmetric if one of its "transpose"
actions was assembled with incompatible quadrature, signs, coordinate maps, or block
orderings.

### 15.1 Why MINRES cares

MINRES is designed for symmetric linear systems and is especially useful for
symmetric-indefinite saddle-point problems. The current `QuadraticKKTSolverPolicy`
therefore permits MINRES only when the product reports that it supports the
symmetric-indefinite requirements.

GMRES does not require symmetry, so it is also available for nonsymmetric KKT products.
The distinction is:

```text
symmetric-indefinite KKT + declared transpose consistency
    → MINRES-compatible

nonsymmetric KKT
    → use a nonsymmetric Krylov policy such as GMRES
```

This is not a ranking of the solvers. It is a compatibility relationship between
operator structure and Krylov method.

## 16. Solver convergence requires both linear and block-residual evidence

A Krylov solver can report that its internal linear residual criterion was satisfied.
For the KKT product, the project also records the resulting physical block residuals:

```text
stationarity residual
equality residual
```

`QuadraticKKTSolveReport::converged()` requires both:

```text
linear solve converged
+
KKT residual blocks passed their convergence check
```

This prevents a solver result from being reported as converged merely because one
internal linear criterion succeeded while the actual stationarity/equality residual
check failed.

The exact numerical realization of the solve lives outside the bare product interface,
but the report shape makes the evidence expected at the boundary explicit.

## 17. The multiplier and adjoint are related but not identical contract roles

The running derivation gave

$$
\lambda=-p.
$$

It would be tempting to treat the multiplier and adjoint as the same layout and rely on
callers to remember the sign. The KKT contract instead stores:

```text
multiplier layout
adjoint layout
multiplier_to_adjoint
adjoint_to_multiplier
```

The conversion can therefore be explicit and checked. For the canonical supplied-OTD
bridge, both operations are simply sign changes:

```math
p=-\lambda,
\qquad
\lambda=-p.
```

The separate roles have two benefits. First, they preserve the project's PDE adjoint
sign convention while allowing the KKT product to use the canonical `+ D^T lambda`
stationarity form.

Second, they make a formulation conversion visible instead of encoding it as an unstated
assumption inside a solver.

## 18. A supplied OTD system is not automatically a quadratic KKT system

This distinction is one of the most important boundaries in the current formulation
layer. A generic supplied OTD system may be:

- nonlinear;
- affine but nonsymmetric;
- based on different block signs;
- missing a canonical primal/equality interpretation;
- based on pairings that do not match the common KKT product;
- lacking the rank or kernel-positivity assumptions needed by the selected quadratic
  solver.

Therefore the existence of

```text
state equation
adjoint equation
control stationarity
```

does **not** prove that the system is an equality-constrained quadratic KKT product. The
canonical adapter requires explicit formulation-owned evidence before performing that
conversion.

## 19. What the supplied-OTD-to-KKT validity record declares

`SuppliedOTDQuadraticKKTValidity` records the evidence needed by the canonical bridge.
The important declarations include:

### Block selection

The adapter must know which supplied variable blocks are

```text
state
adjoint
control
```

and which residual blocks are

```text
state equation
adjoint equation
control stationarity.
```

The declaration must match the actual supplied system.

### Affine residual and constant JVP

The quadratic KKT product assumes a linear KKT operator plus fixed right-hand sides. The
supplied optimality residual must therefore be declared affine, with a point-independent
JVP.

Otherwise extracting $Q$ and $D$ from one linearization point would not define a global
linear KKT operator.

### Canonical block signs

The adapter relies on a particular sign relationship among the state, adjoint, and
stationarity equations. A supplied system with a different convention may still be
mathematically valid, but the canonical adapter must not silently reinterpret its
blocks.

### Pairing and transpose declarations

The bridge needs declared primal/stationarity and multiplier/equality pairings, plus
consistency of the $D^{\mathsf T}$ and whole-KKT transpose actions. These are what
justify calling the adapted product symmetric when that symmetry is claimed.

### Rank and kernel positivity

The bridge carries forward the declared full-row-rank and positivity-on-$\ker(D)$
conditions required by the quadratic KKT product.

### Multiplier conversion

The canonical bridge requires the known relationship

$$
\lambda=-p.
$$

The declaration records that conversion explicitly. None of these boolean declarations
is a theorem prover. They are formulation-owned evidence that the supplied system
satisfies the structural requirements of the narrower product.

## 20. The canonical bridge selects a two-block KKT problem from the three-block OTD system

Once the validity record has been accepted, the adapter

```text
make_canonical_supplied_otd_kkt_product(...)
```

builds a common quadratic KKT view without rewriting the supplied system. The block
transformation is:

```text
supplied OTD variables
[ state | adjoint | control ]
    │        │         │
    │        │         └──────────┐
    │        │                    │
    │        └── λ = -p           │
    │                             │
    └──────────────┐              │
                   ▼              ▼
KKT primal       [ state | control ]

KKT multiplier  [ lambda ]
```

The residual side is reorganized similarly:

```text
supplied OTD residual
[ state equation | adjoint equation | control stationarity ]
        │                 │                   │
        │                 │                   │
        ▼                 └──────────┬────────┘
KKT equality                         ▼
[ state equation ]             KKT stationarity
                               [ state | control ]
```

The adapter does not claim that every three-block system has this interpretation. It
performs this selection only after the validity record has declared the required
linear-quadratic structure and signs.

### 20.1 The equality action comes from the state equation

To evaluate

$$
Dx,
$$

the bridge constructs a supplied-OTD tangent whose:

```text
state block    = KKT primal state component
adjoint block  = 0
control block  = KKT primal control component
```

It applies the supplied JVP and selects the **state-equation residual block**. Because
the supplied residual is declared affine with a constant JVP, this one action represents
the global linear equality operator $D$.

### 20.2 The quadratic action comes from stationarity linearization

The same state/control tangent also produces linearized adjoint-equation and
control-stationarity blocks. Under the canonical supplied sign policy, the bridge uses:

```text
negative supplied adjoint-equation block
+
supplied control-stationarity block
```

as the two KKT stationarity blocks of

$$
Qx.
$$

This sign selection is the implementation counterpart of the state/control stationarity
equations derived earlier. It is deliberately encoded in the adapter rather than left
implicit.

### 20.3 The multiplier action comes from the supplied adjoint variable

To evaluate

$$
D^{\mathsf T}\lambda,
$$

the bridge creates a supplied tangent with:

```text
state block    = 0
adjoint block  = lambda
control block  = 0
```

and applies the supplied JVP. The appropriate adjoint-equation and control-stationarity
components are then selected with the canonical sign conversion.

This is possible because, for an affine optimality system, the JVP contains the constant
block operator relating the supplied adjoint variable to the stationarity equations.

### 20.4 The full KKT transpose is built from the supplied VJP

For a KKT seed consisting of stationarity and equality components, the bridge embeds
those seeds into the supplied residual layout with the required sign adjustment, calls

```text
supplied.residual_vjp(...)
```

and then extracts:

```text
primal transpose result
multiplier transpose result
```

The multiplier component is converted with the same canonical sign convention. So the
bridge reuses the supplied system's own transpose action. It does not synthesize a
transpose by assuming the forward blocks were assembled symmetrically.

## 21. The bridge preserves supplied-OTD provenance

A useful way to understand the adapter is:

```text
supplied OTD system
    owns the original first-order equations
        │
        │ declared linear-quadratic compatibility
        ▼
canonical KKT view
    selects/repackages actions needed by KKT solvers
```

The direction is important. The KKT product is a **view of a supplied OTD system whose
compatibility has been declared**. The adapter does not reverse-engineer a DTO objective
and PDE residual and then claim that the supplied equations came from differentiating
them.

The source comment in `supplied_otd_kkt.hpp` makes this boundary explicit: the adapter
preserves supplied-OTD provenance while selecting the affine blocks needed by the common
quadratic KKT boundary.

This is one of the places where formulation provenance has real mathematical meaning,
not just documentation value.

## 22. Why keep a generic supplied OTD system if a KKT product exists?

Because the two abstractions solve different problems. `SuppliedOTDSystemT` can
represent a broader selected optimality system with:

- a residual evaluated at a point;
- a point-dependent JVP;
- a point-dependent VJP;
- a coupled solve.

That vocabulary can cover nonlinear all-at-once systems. The quadratic KKT product is
narrower:

```text
affine residual
constant block operator
explicit Q / D / D^T structure
quadratic assumptions
Krylov-policy compatibility
```

The narrower structure is valuable because it allows algorithms to exploit more. For
example, a solver can reason about:

- symmetry;
- saddle-point structure;
- rank and kernel positivity;
- MINRES compatibility;
- equality/stationarity residuals.

Those conclusions would be unjustified for an arbitrary nonlinear supplied OTD system.
So the project keeps both levels rather than making the broad type carry false quadratic
promises.

## 23. The quadratic KKT product still does not own every solver detail

The product describes the mathematical operator and its declared assumptions. A separate
`QuadraticKKTSolverPolicy` selects currently supported Krylov behavior:

```text
MINRES
GMRES
```

with numerical controls such as:

```text
maximum iterations
relative tolerance
absolute tolerance
GMRES basis size
```

The policy is validated against the product. A symmetric-indefinite product with the
required transpose declarations can support MINRES. A nonsymmetric product cannot. GMRES
is available without that symmetry requirement.

This separation is analogous to earlier parts of the manual:

```text
numerical product
    what operator exists and what structure is declared

solver policy
    how that operator will be iterated on
```

The product should not pretend that one solver choice is part of the mathematical KKT
definition.

## 24. A small dense KKT test makes the contract concrete

The backend-neutral KKT contract test uses a quadratic action

$$
Q\in\mathbb R^{3\times3}
$$

and an equality action

$$
D\in\mathbb R^{2\times3}.
$$

It constructs the explicit actions

```text
Q x
D x
D^T lambda
KKT transpose
```

and checks that the full action and transpose satisfy the pairing identity. The same
test also checks that:

- a symmetric product reports MINRES compatibility;
- a nonsymmetric product does not;
- multiplier-to-adjoint conversion has the declared sign;
- undeclared rank or kernel conditions are rejected;
- incomplete/incompatible block pairings are rejected;
- missing transpose-consistency declarations are rejected.

These tests are useful because they isolate the formulation contract from any particular
PDE or deal.II realization. The mathematics is small enough to inspect directly.

## 25. The supplied-OTD tests check a different level

The supplied-OTD contract test deliberately uses a more general three-block residual
with arbitrary dense block matrices. Its purpose is not to reproduce one PDE optimality
system.

It checks that:

- residual blocks are selected correctly;
- the supplied JVP agrees with the affine residual;
- the supplied VJP satisfies the pairing identity;
- the coupled solve preserves its solve report;
- invalid variable/residual layouts are rejected;
- lifetime ownership is retained across compiled/supplied products.

Separate tests then exercise the canonical OTD-to-KKT bridge. This test split reflects
the conceptual split in the code:

```text
supplied OTD correctness
        is not the same test as
quadratic KKT structure
```

## 26. A useful comparison with the reduced formulation

The running linear-quadratic problem can now be viewed through two execution
organizations.

### Reduced DTO

```text
control u
    ↓
solve A y = f + B u
    ↓
solve A^T p = M_y y - q
    ↓
reduced derivative
β N_u u + B^T p
    ↓
outer optimization method
```

### All-at-once KKT

```text
unknown (y, u, lambda)
        │
        ▼
[ Q   D^T ] [ x      ] = [ c_x ]
[ D    0  ] [ lambda ]   [ d   ]
        │
        ▼
one coupled saddle-point solve
```

with

$$
x=
\begin{bmatrix}
y\\
u
\end{bmatrix},
\qquad
\lambda=-p.
$$

The mathematics is consistent between the two views. What changes is the numerical
organization:

```text
reduced
    eliminate state/adjoint through nested solves

all-at-once
    retain all first-order unknowns in one coupled system
```

This is a formulation choice, not merely a different linear-solver API.

## 27. What the chapter deliberately does not cover

The KKT system derived here contains only equality constraints. A box-constrained
control introduces additional first-order structure. The simple stationarity equation

$$
D_{u}J_{h}-D_{u}E_{h}^{\ast}p=0
$$

is replaced by a constrained stationarity condition involving, depending on the
representation:

- a variational inequality;
- a normal cone;
- lower/upper multipliers;
- complementary slackness;
- active and free sets.

The quadratic KKT product in this chapter is therefore the linear algebraic core on
which the project's PDAS subproblems are built, not the complete description of a
bound-constrained problem.

That distinction is the subject of **Complementarity and PDAS**.

## 28. Following the optimality-system path through the source

A focused source tour is now possible.

### Supplied all-at-once system

Read:

- [`include/nmopt/contract/supplied_otd.hpp`](../../../include/nmopt/contract/supplied_otd.hpp)

A useful order is:

```text
SuppliedOTDBlockSelection
SuppliedOTDLayout
SuppliedOTDSystemT
SuppliedOTDQuadraticKKTValidity
```

The first three explain the broad three-block supplied system. The validity record is
narrower and should be read only after the general product is understood.

### Common quadratic KKT boundary

Then read:

- [`include/nmopt/contract/quadratic_kkt.hpp`](../../../include/nmopt/contract/quadratic_kkt.hpp)

Look for:

```text
QuadraticKKTAssumptions
QuadraticKKTBlockPairing
QuadraticKKTLayoutT
EqualityConstrainedQuadraticKKTProductT
```

The important actions are:

```text
apply_q
apply_d
apply_d_transpose
apply_kkt
residual
apply_kkt_transpose
```

### Supplied-OTD bridge

Next read:

- [`include/nmopt/contract/supplied_otd_kkt.hpp`](../../../include/nmopt/contract/supplied_otd_kkt.hpp)

The file first validates the explicit compatibility declarations and only then
constructs the KKT layouts/actions. Follow the three forward actions before reading the
full transpose adapter.

### Solver policy

Then read:

- [`include/nmopt/contract/quadratic_kkt_solver.hpp`](../../../include/nmopt/contract/quadratic_kkt_solver.hpp)

This file is intentionally small. It contains the MINRES/GMRES compatibility rule and
the solve-report shape.

### Contract tests

Finally compare:

- [`tests/contract/supplied_otd_contract.cc`](../../../tests/contract/supplied_otd_contract.cc)
- [`tests/contract/quadratic_kkt_contract.cc`](../../../tests/contract/quadratic_kkt_contract.cc)

The former tests the broad supplied-system boundary. The latter tests the structured
quadratic saddle-point boundary.

## 29. A compact correspondence

| Mathematical idea | Role | Current contract object |
| --- | --- | --- |
| all-at-once unknown | $(y,p,u)$ | supplied variable layout |
| all-at-once residual | state + adjoint + stationarity equations | `SuppliedOTDSystemT::residual` |
| optimality JVP | $`F_{h}'(w)[\delta w]`$ | `residual_jvp` |
| optimality VJP | $`F_{h}'(w)^{\ast}q`$ | `residual_vjp` |
| coupled optimality solve | $`F_{h}(w)=0`$ | `SuppliedOTDSystemT::solve` |
| quadratic primal | $x$ | KKT primal layout |
| equality multiplier | $\lambda$ | KKT multiplier layout |
| objective curvature | $Qx$ | `apply_q` |
| equality action | $Dx$ | `apply_d` |
| equality transpose | $D^{\mathsf T}\lambda$ | `apply_d_transpose` |
| KKT residual | $`(Qx+D^{\mathsf T}\lambda-c_{x},\ Dx-d)`$ | `residual` |
| PDE/KKT sign conversion | $\lambda=-p$ | multiplier conversion |
| KKT transpose | $K^{\mathsf T}$ action | `apply_kkt_transpose` |
| saddle-point assumptions | rank + positivity on $\ker(D)$ | `QuadraticKKTAssumptions` |
| symmetric Krylov compatibility | declared symmetric-indefinite structure | `supports_minres()` |

## 30. What to carry into complementarity and PDAS

The chapter has established the equality-constrained first-order system. For a quadratic
problem, the core solve is

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
c_{x}\\
d
\end{bmatrix}.
```

Chapter 8 adds bounds. That changes the problem in two linked ways. First,
control stationarity is no longer simply zero in every coefficient.

Second, each bound-constrained coefficient can switch between:

```text
free
active at the lower bound
active at the upper bound
```

PDAS repeatedly classifies those roles and solves a KKT system restricted by the current
active-set hypothesis. So the quadratic KKT product developed here is not discarded.

It becomes the linear algebraic core of the active-set subproblem.

## Read later

Useful existing documents are:

- [Theoretical formalism](../../design/mathematical-model.md), for the project's
  Lagrangian sign, dual-pairing, DTO/OTD, and first-order conventions.
- [Project architecture](../overview/project-architecture.md), for the relationship
  among reduced, supplied OTD, KKT, and PDAS formulation products.
- [Semantic compiler](../overview/semantic-compiler.md), for the high-level
  compiler path that can construct selected formulation products.
- [v1 semantic/compiler capability](../../internals/compiler/semantic-compiler.md),
  for the exact currently registered formulation/capability ledger.

The next concept chapter will connect the equality-constrained KKT system to bound
constraints, complementarity, and primal-dual active-set iterations.
