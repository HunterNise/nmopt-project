# Metrics, gradients, and constraints

## The derivative tells us how the objective changes, not how to move

The previous chapter ended with a reduced derivative

$$
j_{h}'(u)\in U_{h}^{\ast}.
$$

That object answers a precise question:

> If the control changes by $`\delta u\in U_{h}`$, what is the first-order change in
> the reduced objective?

The answer is

$$
j_{h}'(u)[\delta u].
$$

In coordinates,

$$
j_{h}'(u)[\delta u]
=
\mathbf r^{\mathsf T}\delta\mathbf u,
$$

where $\mathbf r$ contains the dual coordinates of the derivative.

An optimization algorithm needs something slightly different. It needs a **primal
direction** in which to move the control.

That is the role of a metric or, more precisely, a Riesz map.

The distinction is easiest to see if we keep the two questions separate:

```text
derivative
    how does the objective act on a proposed direction?

gradient
    which primal vector represents that derivative
    under a chosen inner product?
```

This chapter develops that distinction carefully and then follows it into the actual
metric and constraint implementations in `nmopt`.

The narrative is:

```text
reduced derivative in U_h*
        ↓
choose an inner product on U_h
        ↓
Riesz map G : U_h → U_h*
        ↓
gradient g = G^{-1} j'(u)
        ↓
search direction / norm / stopping
        ↓
constraint projection in the same geometry
```

The key point is that the metric is not a decorative norm attached to the solver.
It is the operator that identifies primal control directions with dual derivatives.

## 1. The same derivative has different gradients in different geometries

Let $`U_{h}`$ be a finite-dimensional control space.

Suppose the reduced derivative at $u$ is the covector

$$
r
:=
j_{h}'(u)
\in
U_{h}^{\ast}.
$$

Choose an inner product

$$
(\cdot,\cdot)_{G}
$$

on $`U_{h}`$.

The **gradient with respect to that inner product** is the unique
$`g\in U_{h}`$ satisfying

$$
r[\delta u]
=
(g,\delta u)_{G}
\qquad
\text{for every }\delta u\in U_{h}.
$$

This is the finite-dimensional form of the Riesz representation theorem.

The important word is **chosen**.

The derivative $r$ is fixed by the objective and PDE model.

The gradient depends on the geometry used to identify $`U_{h}`$ with its dual.

### 1.1 Coordinate form

Let

$$
G:
U_{h}
\longrightarrow
U_{h}^{\ast}
$$

denote the Riesz map associated with the chosen inner product.

By definition,

$$
(v,w)_{G}
:=
\langle Gv,w\rangle.
$$

If $\mathbf G$ is its matrix in the chosen basis, then

$$
(v,w)_{G}
=
\mathbf v^{\mathsf T}\mathbf G^{\mathsf T}\mathbf w.
$$

For the symmetric positive-definite realizations used here,

$$
\mathbf G^{\mathsf T}
=
\mathbf G.
$$

The gradient condition becomes

$$
\mathbf r^{\mathsf T}\delta\mathbf u
=
\mathbf g^{\mathsf T}\mathbf G\delta\mathbf u
\qquad
\text{for every }\delta\mathbf u.
$$

Therefore

$$
\mathbf G\mathbf g
=
\mathbf r,
$$

and hence

$$
\mathbf g
=
\mathbf G^{-1}\mathbf r.
$$

This is the operation implemented by `MetricT::inverse_apply()`.

### 1.2 Euclidean coordinates are one geometry, not the default truth

If

$$
\mathbf G=I,
$$

then

$$
\mathbf g=\mathbf r.
$$

That is the familiar Euclidean rule.

It is easy to internalize this case so deeply that one begins to identify derivatives
and gradients in general.

Finite-element optimization is precisely where that shortcut becomes dangerous.

The coefficient vector representing a derivative is often **not** the coefficient
vector of the physically meaningful gradient.

## 2. Why the negative metric gradient is steepest descent

The phrase *steepest descent* is also metric-dependent.

Suppose we consider all unit directions in the chosen norm,

$$
\lVert d\rVert_{G}
:=
\sqrt{(d,d)_{G}}.
$$

Among those directions, we want the one that produces the most negative first-order
change in the objective:

```math
\min_{\lVert d\rVert_{G}=1}
j_{h}'(u)[d].
```

Because

$$
j_{h}'(u)[d]
=
(g,d)_{G},
$$

the Cauchy–Schwarz inequality gives

$$
(g,d)_{G}
\geq
-\lVert g\rVert_{G}\lVert d\rVert_{G}
=
-\lVert g\rVert_{G}.
$$

Equality is attained for

$$
d
=
-\frac{g}{\lVert g\rVert_{G}}.
$$

So the negative metric gradient really is the direction of steepest local decrease
**in the chosen geometry**.

Changing the metric therefore changes what "steepest" means.

### 2.1 The unnormalized direction used by an algorithm

A line-search method usually does not normalize the gradient.

It uses

$$
d=-g
$$

and lets the step length $\alpha$ determine how far to move:

$$
u_{\mathrm{trial}}
=
u+\alpha d.
$$

The directional derivative is

```math
j_{h}'(u)[d]
=
-\lVert g\rVert_{G}^{2}.
```

Indeed,

```math
\begin{aligned}
j_{h}'(u)[-g]
&=
-\langle Gg,g\rangle\\
&=
-\lVert g\rVert_{G}^{2}.
\end{aligned}
```

This is exactly the descent relation that the reduced search implementation checks
after obtaining the metric gradient.

## 3. The finite-element $L^{2}$ metric

Return to the distributed-control example.

Let the control be

$$
u_{h}
=
\sum_{i=1}^{n}u_{i}\psi_{i},
$$

with coefficient vector $\mathbf u$.

The finite-element $L^{2}$ inner product is

$$
(u_{h},v_{h})_{L^{2}(\Omega)}
=
\int_{\Omega}u_{h}v_{h}.
$$

In coordinates,

$$
=
\mathbf u^{\mathsf T}M\mathbf v,
$$

where the mass matrix is

$$
M_{ij}
:=
\int_{\Omega}\psi_{j}\psi_{i}.
$$

Therefore the discrete $L^{2}$ Riesz map is

$$
G=M.
$$

If the reduced derivative has covector coordinates $\mathbf r$, then the
$L^{2}$ gradient is found from

$$
M\mathbf g
=
\mathbf r.
$$

This is not merely a preconditioning trick added for numerical convenience. It is the
coordinate form of the mathematical identification between $`U_{h}`$ and
$`U_{h}^{\ast}`$ induced by the $L^{2}$ inner product.

### 3.1 The norm comes from the same map

Once $g$ is known,

$$
\lVert g\rVert_{L^{2}}^{2}
=
\mathbf g^{\mathsf T}M\mathbf g.
$$

Because

$$
M\mathbf g=\mathbf r,
$$

the same quantity can also be written as

$$
\mathbf r^{\mathsf T}\mathbf g.
$$

This identity appears naturally in the solver:

```text
gradient = metric.inverse_apply(derivative)
metric_gradient = metric.apply(gradient)
norm^2 = pair(metric_gradient, gradient)
```

The three operations correspond exactly to

```math
g=G^{-1}r,
\qquad
Gg,
\qquad
\langle Gg,g\rangle.
```

## 4. Regularization and optimization geometry are different choices

The running objective contains

$$
\frac{\beta}{2}
\lVert u\rVert_{L^{2}}^{2}.
$$

Its derivative contributes

$$
\beta Mu.
$$

The same mass matrix $M$ may also be used as the optimization metric.

That coincidence can make it look as though the regularization term somehow
*defines* the gradient geometry.

It does not.

These are two separate uses of $M$.

### Regularization

The objective contains

$$
R(u)
:=
\frac{\beta}{2}
u^{\mathsf T}Mu.
$$

Therefore

$$
R'(u)
=
\beta Mu.
$$

Changing this term changes the optimization problem itself.

### Metric

The solver chooses

$$
G:=M
$$

to identify a derivative covector with an $L^{2}$ gradient.

Changing $G$ changes the representation of the gradient and the search geometry, but
does not by itself change the objective functional.

The distinction can be pictured as:

```text
objective definition
    β/2 <M u, u>
        │
        ▼
derivative contribution
    β M u
        │
        └───────────────┐
                        │
other derivative terms  │
        │               │
        └──────┬────────┘
               ▼
        total covector r
               │
               │ choose optimization metric G
               ▼
        gradient g = G^-1 r
```

Using $G=M$ is common and often natural, but it is still a solver/modeling choice
separate from the regularization coefficient and term.

## 5. Mesh dependence makes the distinction practically important

Suppose one ignores the mass matrix and interprets the derivative coefficients
$\mathbf r$ directly as a primal gradient.

That uses the Euclidean coefficient norm

$$
\lVert \mathbf v\rVert_{2}^{2}
=
\sum_{i}v_{i}^{2}.
$$

This norm depends strongly on the coordinate representation.

Under mesh refinement, the number and scaling of basis functions change. A fixed
physical field can therefore have very different Euclidean coefficient norms on
different meshes.

By contrast, the finite-element mass norm approximates the physical $L^{2}$ norm:

$$
\lVert v_{h}\rVert_{L^{2}}^{2}
=
\mathbf v^{\mathsf T}M\mathbf v.
$$

This is one reason Sobolev or mass-matrix metrics are useful in PDE optimization:
they express geometry in terms of the represented field rather than the arbitrary
raw coefficient scaling.

The [Numerical realization](../overview/numerical-realization.md) overview places these metric choices in the broader
discretization layer. Here the practical message is that a metric is part of the
numerical meaning of a gradient.

## 6. What `MetricT` represents

The contract layer packages the primal–dual identification as `MetricT<Backend>`.

Its essential operations are

```cpp
apply(primal)
inverse_apply(covector)
```

with the mathematical interpretation

```math
\begin{aligned}
\mathrm{apply}(v)
&=
Gv
\in
U_{h}^{\ast},
\\
\mathrm{inverse\_apply}(r)
&=
G^{-1}r
\in
U_{h}.
\end{aligned}
```

The C++ return types make the two directions explicit:

```text
PrimalBlockT
    --apply-->
CovectorBlockT

CovectorBlockT
    --inverse_apply-->
PrimalBlockT
```

This is the runtime version of the Riesz-map diagram.

### 6.1 The current interface is a fixed linear metric on one layout

`MetricT` stores no evaluation point.

The current abstraction therefore represents a fixed primal–dual map for one runtime
layout, not a general point-dependent Riemannian metric of the form

$$
G(u).
$$

That is enough for the mass, negative-order, and trace metrics currently realized by
the project and for the reduced algorithms that consume them.

The word *metric* in this code should therefore be read as something close to

> a constructed Hilbert/Riesz realization that maps between one primal coefficient
> space and its dual.

That is more specific than the many meanings of "metric" used elsewhere in numerical
analysis.

### 6.2 `apply()` can be cheap while `inverse_apply()` is a solve

The interface does not imply that both directions have the same cost.

For a sparse SPD mass matrix,

$$
Gv=Mv
$$

is one matrix-vector product.

Applying the inverse,

$$
G^{-1}r=M^{-1}r,
$$

requires solving a linear system.

For more involved metrics, even `apply()` may contain a solve.

This is why the reduced solver records metric-solve work separately from state and
adjoint solves: obtaining a gradient can itself be a nontrivial numerical operation.

## 7. `MassMetric`: the most direct finite-element realization

`nmopt::dealii_backend::MassMetric` realizes a one-block metric from an explicitly
supplied symmetric positive-definite sparse matrix.

For a standard $L^{2}$ control metric that matrix is the finite-element mass matrix
$M$.

The two operations are literal:

```text
apply(v)
    M v

inverse_apply(r)
    solve M g = r
```

The implementation uses deal.II's sparse matrix `vmult` for the first operation and
conjugate gradients for the inverse.

The layout fixes which coefficient space the matrix acts on.

### 7.1 Why a solve tolerance appears in a metric

Mathematically,

$$
g=M^{-1}r
$$

is exact.

Numerically, an iterative linear solver returns an approximation.

So a metric realization needs numerical policy such as:

- maximum iterations;
- relative tolerance;
- absolute tolerance.

This is an instructive example of the difference between the **abstract metric** and
its **numerical realization**.

The abstract statement is

$$
G^{-1}:U_{h}^{\ast}\to U_{h}.
$$

The deal.II realization must decide how accurately to apply that inverse.

If the metric solve is too inaccurate, the resulting vector is not a sufficiently
accurate representation of the chosen gradient. That can affect descent tests,
stopping criteria, and higher-level optimization behavior.

## 8. Cellwise-constant controls produce a special mass matrix

The cellwise-constant control realization introduced in the first chapter has an
important simplification.

Suppose the control basis contains one indicator-like basis function $`\psi_{K}`$ for
each mesh cell $K$.

Different basis functions have disjoint support, so

$$
\int_{\Omega}\psi_{K}\psi_{L}=0
\qquad
K\neq L.
$$

The mass matrix is therefore diagonal.

For the simplest unscaled cell indicator basis,

$$
M_{KK}
=
|K|,
$$

the measure of the cell.

Thus

$$
M
=
\mathrm{diag}
\left(
|K_{1}|,\ldots,|K_{n}|
\right).
$$

The $L^{2}$ gradient is then obtained coefficient by coefficient:

$$
g_{K}
=
\frac{r_{K}}{|K|}.
$$

This is a useful case because several ideas coincide cleanly:

- each coefficient represents one cell value;
- the $L^{2}$ metric is diagonal;
- the Riesz inverse is a diagonal scaling;
- coefficientwise box constraints have a simple metric projection.

The last point becomes important later in the chapter.

### 8.1 Continuous finite elements do not have the same simplification

For a continuous nodal space such as $`\mathbb P_{1}`$ or $`\mathbb Q_{1}`$, neighboring
basis functions overlap.

The mass matrix therefore has off-diagonal entries.

So even though both discretizations represent an $L^{2}$ control and both use a
finite-element mass matrix, their coefficient geometry differs:

```text
cellwise constant
    diagonal mass matrix

continuous nodal
    sparse coupled mass matrix
```

This is a concrete reason not to equate the phrase "$L^{2}$ metric" with
"coefficientwise scaling".

## 9. The metric need not be an $L^{2}$ metric

The Riesz-map construction works with any chosen Hilbert inner product that has a
suitable discrete realization.

For a control space carrying more or less regularity than $L^{2}$, another geometry
may be more natural.

The repository contains two especially instructive examples:

- a discrete $H^{-1}$ metric;
- a quotient $H^{1/2}$ trace metric.

These are useful not only because they extend the supported problems. They show that
`MetricT` really represents an **operator-defined primal–dual identification**, not a
thin wrapper around a norm formula.

## 10. A discrete $H^{-1}$ metric

The negative-order metric is easier to understand if we first separate the two
finite-element matrices that enter its definition.

### 10.1 Mass and stiffness matrices measure different things

Let

$$
V_{h}
:=
\mathrm{span}\{\phi_{1},\ldots,\phi_{n}\}
\subset H_{0}^{1}(\Omega).
$$

For coefficient vectors $\mathbf v$ and $\mathbf w$, write the corresponding fields
as $`v_{h}`$ and $`w_{h}`$.

The **mass matrix** is

$$
M_{ij}
:=
\int_{\Omega}\phi_{j}\phi_{i}.
$$

It represents the $L^{2}$ pairing:

$$
\mathbf v^{\mathsf T}M\mathbf w
=
\int_{\Omega}v_{h}w_{h}.
$$

The **stiffness matrix** for the homogeneous-Dirichlet Laplacian is

$$
K_{ij}
:=
\int_{\Omega}
\nabla\phi_{j}\cdot\nabla\phi_{i}.
$$

It represents the Dirichlet energy pairing:

$$
\mathbf v^{\mathsf T}K\mathbf w
=
\int_{\Omega}
\nabla v_{h}\cdot\nabla w_{h}.
$$

So $M$ measures the size of the fields themselves, whereas $K$ measures their
spatial variation. On $`H_{0}^{1}(\Omega)`$, the zero boundary condition removes the
constant nullspace, and the Laplace stiffness matrix is positive definite under the
usual assumptions.

Both matrices map primal coefficient vectors to covectors, but they encode different
bilinear forms.

### 10.2 Why a potential appears in an $H^{-1}$ norm

The notation $H^{-1}$ means the dual of $`H_{0}^{1}`$, so its norm is naturally defined
through how a source acts on $`H_{0}^{1}`$ test functions.

Using the Dirichlet-energy norm for the test space, one may think of

```math
\lVert u\rVert_{H^{-1}}
=
\sup_{v\in H_{0}^{1}(\Omega),\ v\neq0}
\frac{\langle u,v\rangle}
{\left(\int_{\Omega}|\nabla v|^{2}\right)^{1/2}}.
```

When the source $u$ is represented by an ordinary $L^{2}$ function, its action is

$$
\langle u,v\rangle
=
\int_{\Omega}uv.
$$

Instead of evaluating the supremum directly, introduce the Riesz representative
$`w\in H_{0}^{1}(\Omega)`$ defined by

$$
\int_{\Omega}\nabla w\cdot\nabla v
=
\int_{\Omega}uv
\qquad
\text{for every }v\in H_{0}^{1}(\Omega).
$$

In strong-form language this is the Poisson problem

$$
-\Delta w=u,
\qquad
w=0\text{ on }\partial\Omega.
$$

The field $w$ is often called the **potential** generated by the source $u$.
Nothing mysterious is intended by that word: $w$ is simply the solution obtained by
applying the inverse Dirichlet Laplacian to $u$.

Why is this useful for the norm? Set $v=w$ in the weak equation:

```math
\int_{\Omega}|\nabla w|^{2}
=
\int_{\Omega}u w.
```

The Riesz-representation argument gives

$$
\lVert u\rVert_{H^{-1}}^{2}
=
\int_{\Omega}|\nabla w|^{2}
=
\int_{\Omega}u w.
$$

So the negative norm measures the source through the energy of the potential it
creates.

This also explains qualitatively why the norm is called *negative order*. Solving an
elliptic equation smooths the source before measuring it. Rapidly oscillating source
components are therefore weighted differently from low-frequency components.

### 10.3 Discretizing the potential equation

Let $\mathbf u$ contain the coefficients of the discrete source field and let
$\mathbf w$ contain the coefficients of its potential.

Testing the weak potential equation with every basis function gives

$$
K\mathbf w
=
M\mathbf u.
$$

The right-hand side is $M\mathbf u$, rather than just $\mathbf u$, because
$\mathbf u$ contains **primal field coefficients**. The load covector with entries

$$
\int_{\Omega}u_{h}\phi_{i}
$$

is exactly $M\mathbf u$.

Solving for the potential gives

$$
\mathbf w
=
K^{-1}M\mathbf u.
$$

Now evaluate its Dirichlet energy:

```math
\begin{aligned}
\mathbf w^{\mathsf T}K\mathbf w
&=
\mathbf w^{\mathsf T}M\mathbf u
\\
&=
\mathbf u^{\mathsf T}M\mathbf w
\\
&=
\mathbf u^{\mathsf T}
M K^{-1}M
\mathbf u.
\end{aligned}
```

The first equality uses the discrete potential equation
$K\mathbf w=M\mathbf u$; the second uses symmetry of $M$.

Hence the project uses the discrete negative-order Riesz map

$$
G
:=
M K^{-1}M.
$$

The norm induced by that map is

$$
\lVert \mathbf u\rVert_{G}^{2}
=
\mathbf u^{\mathsf T}G\mathbf u.
$$

The current `Hminus1Metric` implements this operator exactly in the factorized form
shown above rather than assembling a matrix product $M K^{-1}M$ explicitly.

### 10.4 Applying the $H^{-1}$ metric

To compute

$$
G\mathbf u
=
M K^{-1}M\mathbf u,
$$

the implementation follows the algebra from right to left:

```text
primal source u
      │
      ▼
form its L2 load covector
      M u
      │
      ▼
solve the potential problem
      K w = M u
      │
      ▼
convert the potential to an L2 covector
      M w
      │
      ▼
metric covector G u
```

The middle solve is the discrete Poisson solve that creates the potential.

The final mass action is needed because `MetricT::apply()` must return a **covector**.
The potential coefficients $\mathbf w$ are primal coefficients; $M\mathbf w$ is the
covector that represents their $L^{2}$ action.

Unlike the simple mass metric, applying $G$ itself therefore requires a linear solve.

### 10.5 Applying the inverse metric

From

$$
G
=
M K^{-1}M,
$$

we obtain

$$
G^{-1}
=
M^{-1}K M^{-1}.
$$

One can verify the factorization directly:

```math
\begin{aligned}
G G^{-1}
&=
M K^{-1}M
M^{-1}K M^{-1}
\\
&=
M K^{-1}K M^{-1}
\\
&=
I.
\end{aligned}
```

So for a derivative covector $\mathbf r$, the metric gradient

$$
\mathbf g
=
G^{-1}\mathbf r
$$

is obtained as follows:

```text
covector r
    │
    ▼
first mass solve
    M a = r
    │
    ▼
Laplacian action
    q = K a
    │
    ▼
second mass solve
    M g = q
    │
    ▼
primal gradient g
```

The first solve converts the covector $\mathbf r$ into its $L^{2}$ primal
representative $\mathbf a$. The stiffness action then applies the Dirichlet-energy
operator to that field. The final mass solve converts the resulting covector back to
a primal coefficient vector.

The implementation checks convergence of each SPD solve before returning the result.

### 10.6 Why a weaker metric changes the gradient

The derivative covector is the same no matter which metric is chosen.

With an $L^{2}$ metric,

$$
\mathbf g_{L^{2}}
=
M^{-1}\mathbf r.
$$

With the discrete $H^{-1}$ metric,

$$
\mathbf g_{H^{-1}}
=
M^{-1}K M^{-1}\mathbf r.
$$

The extra stiffness action changes the relative weight of spatial modes. The precise
spectral effect depends on the discretization, but the conceptual point is already
clear: choosing a negative-order metric changes which primal perturbations count as
large or small and therefore changes the meaning of steepest descent.

It does not change $j'(u)$ itself.

## 11. A trace $H^{1/2}$ metric arises from minimum extension

Boundary control creates a different geometric problem. The control coefficients
live only on a selected part of the boundary, while a familiar $H^{1}$ energy is
naturally defined for functions in the **volume**.

The $H^{1/2}$ trace norm connects those two settings through extension.

### 11.1 From a boundary trace to a volume norm

A trace field $b$ on the boundary can usually be extended into the volume in many
ways. If $v$ is one such extension, then all we know is

$$
v|_{\Gamma}=b.
$$

Different choices of the interior values of $v$ can have very different $H^{1}$
energies.

The quotient or minimum-extension viewpoint defines the trace norm by choosing the
least expensive extension:

```math
\lVert b\rVert_{H^{1/2}(\Gamma)}^{2}
\sim
\min_{v|_{\Gamma}=b}
\lVert v\rVert_{H^{1}(\Omega)}^{2}.
```

The exact continuous norm-equivalence constants are not important here. What matters
for the implementation is the construction: **hold the trace fixed and let the
interior adjust to minimize the volume energy**.

This explains why we must distinguish boundary and interior degrees of freedom.
The boundary coefficients are the actual trace/control variables. The interior
coefficients are auxiliary extension variables whose only job is to realize the
minimum-energy continuation of that trace into the volume.

If we simply filled the interior with zeros, we would obtain one arbitrary extension,
not the quotient norm.

### 11.2 The volume $H^{1}$ matrix

Let

$$
A
$$

be the symmetric positive-definite matrix representing the chosen volume $H^{1}$
inner product. In the simplest case it is a mass-plus-stiffness matrix:

$$
A
=
M+K,
$$

so that

$$
\mathbf v^{\mathsf T}A\mathbf v
\approx
\int_{\Omega}
\left(
|v_{h}|^{2}+|\nabla v_{h}|^{2}
\right).
$$

Now reorder the volume degrees of freedom into two groups:

- $B$: the selected **boundary/trace** DoFs;
- $I$: the remaining **interior** DoFs.

With this ordering,

```math
A
=
\begin{bmatrix}
A_{BB} & A_{BI}\\
A_{IB} & A_{II}
\end{bmatrix}.
```

The subscripts describe both the input coefficients and the covector rows:

```text
A_BB
    trace coefficients  → trace covector contribution

A_BI
    interior coefficients → trace covector contribution

A_IB
    trace coefficients → interior covector contribution

A_II
    interior coefficients → interior covector contribution
```

Because $A$ is symmetric,

$$
A_{IB}
=
A_{BI}^{\mathsf T}.
$$

Let

$$
\mathbf b
$$

denote the fixed trace coefficients and

$$
\mathbf i
$$

denote the interior coefficients that we are free to choose.

The full volume extension is therefore

$$
\mathbf e
:=
\begin{bmatrix}
\mathbf b\\
\mathbf i
\end{bmatrix}.
$$

### 11.3 Minimize over the interior coefficients

The extension energy is

```math
\mathcal E(\mathbf b,\mathbf i)
:=
\frac{1}{2}
\begin{bmatrix}
\mathbf b\\
\mathbf i
\end{bmatrix}^{\mathsf T}
A
\begin{bmatrix}
\mathbf b\\
\mathbf i
\end{bmatrix}.
```

Expanding the block product gives

```math
\mathcal E(\mathbf b,\mathbf i)
=
\frac{1}{2}\mathbf b^{\mathsf T}A_{BB}\mathbf b
+
\mathbf b^{\mathsf T}A_{BI}\mathbf i
+
\frac{1}{2}\mathbf i^{\mathsf T}A_{II}\mathbf i.
```

The trace $\mathbf b$ is fixed. We minimize only over $\mathbf i$.

For an interior perturbation $\delta\mathbf i$,

```math
D_{\mathbf i}\mathcal E(\mathbf b,\mathbf i)
[\delta\mathbf i]
=
\delta\mathbf i^{\mathsf T}
\left(
A_{IB}\mathbf b+A_{II}\mathbf i
\right).
```

At the minimum this must vanish for every $\delta\mathbf i$, so

$$
A_{IB}\mathbf b
+
A_{II}\mathbf i
=
0.
$$

Hence the minimum-energy interior coefficients are

$$
\mathbf i_{\ast}(\mathbf b)
=
-A_{II}^{-1}A_{IB}\mathbf b.
$$

This is the interior solve performed by the trace metric.

### 11.4 The Schur complement is the trace Riesz map

Substitute the minimizing interior coefficients into the boundary block of the volume
action:

```math
\begin{aligned}
A_{BB}\mathbf b
+
A_{BI}\mathbf i_{\ast}
&=
A_{BB}\mathbf b
-
A_{BI}A_{II}^{-1}A_{IB}\mathbf b
\\
&=
\left(
A_{BB}
-
A_{BI}A_{II}^{-1}A_{IB}
\right)\mathbf b.
\end{aligned}
```

Define

$$
G
:=
A_{BB}
-
A_{BI}A_{II}^{-1}A_{IB}.
$$

This is the Schur complement obtained by eliminating the interior extension
coefficients.

It is also the trace Riesz map.

To see the energy interpretation explicitly, substitute
$`\mathbf i_{\ast}`$ into the quadratic form. The minimum extension energy becomes

$$
\mathcal E(\mathbf b,\mathbf i_{\ast})
=
\frac{1}{2}
\mathbf b^{\mathsf T}G\mathbf b.
$$

So $G$ measures the boundary trace by the $H^{1}$ energy of its best volume
extension.

### 11.5 Applying $G$ without assembling the Schur complement

The formula

$$
G
=
A_{BB}-A_{BI}A_{II}^{-1}A_{IB}
$$

might suggest constructing a separate Schur-complement matrix.
`TraceHhalfMetric` does not need to do that.

Given trace coefficients $\mathbf b$:

1. place $\mathbf b$ in the trace positions of a full volume vector;
2. solve
   $$
   A_{II}\mathbf i_{\ast}
   =
   -A_{IB}\mathbf b;
   $$
3. fill the interior positions with $`\mathbf i_{\ast}`$;
4. apply the full volume matrix $A$;
5. retain only the trace rows.

Why does that produce $G\mathbf b$?

The full action on the minimum extension is

```math
A
\begin{bmatrix}
\mathbf b\\
\mathbf i_{\ast}
\end{bmatrix}
=
\begin{bmatrix}
A_{BB}\mathbf b+A_{BI}\mathbf i_{\ast}\\
A_{IB}\mathbf b+A_{II}\mathbf i_{\ast}
\end{bmatrix}.
```

The interior block is zero by the minimum-extension equation. The trace block is
exactly

$$
G\mathbf b.
$$

Therefore

```math
A
\begin{bmatrix}
\mathbf b\\
\mathbf i_{\ast}
\end{bmatrix}
=
\begin{bmatrix}
G\mathbf b\\
0
\end{bmatrix}.
```

Restricting the full action to the trace DoFs gives the metric covector without ever
forming $G$ as a standalone matrix.

This is the strategy used by `TraceHhalfMetric::apply_vector()`.

### 11.6 Applying $G^{-1}$ with one full-volume solve

Now suppose we are given a trace covector

$$
\mathbf r_{B}
$$

and want the primal trace vector

$$
\mathbf x_{B}
=
G^{-1}\mathbf r_{B}.
$$

One possibility would be:

1. assemble the Schur complement $G$;
2. solve $`G\mathbf x_{B}=\mathbf r_{B}`$.

The implementation uses a more direct strategy. Embed the trace covector into the
full volume dual space by setting the interior right-hand side to zero, and solve

```math
A
\begin{bmatrix}
\mathbf x_{B}\\
\mathbf x_{I}
\end{bmatrix}
=
\begin{bmatrix}
\mathbf r_{B}\\
0
\end{bmatrix}.
```

Write this as two block equations:

```math
\begin{aligned}
A_{BB}\mathbf x_{B}
+
A_{BI}\mathbf x_{I}
&=
\mathbf r_{B},
\\
A_{IB}\mathbf x_{B}
+
A_{II}\mathbf x_{I}
&=
0.
\end{aligned}
```

The second equation implies

$$
\mathbf x_{I}
=
-A_{II}^{-1}A_{IB}\mathbf x_{B}.
$$

Substitute this into the first equation:

```math
\begin{aligned}
\mathbf r_{B}
&=
A_{BB}\mathbf x_{B}
+
A_{BI}
\left(
-A_{II}^{-1}A_{IB}\mathbf x_{B}
\right)
\\
&=
\left(
A_{BB}
-
A_{BI}A_{II}^{-1}A_{IB}
\right)
\mathbf x_{B}
\\
&=
G\mathbf x_{B}.
\end{aligned}
```

Therefore the trace part of the full-volume solution is precisely

$$
\mathbf x_{B}
=
G^{-1}\mathbf r_{B}.
$$

The full solve has performed the Schur-complement elimination implicitly.

This is why `TraceHhalfMetric::inverse_apply_vector()` can:

```text
trace covector r_B
        │
        ▼
embed as [r_B, 0] in the volume
        │
        ▼
solve the full volume H1 system A x = rhs
        │
        ▼
restrict x to the trace DoFs
        │
        ▼
trace primal G^-1 r_B
```

The two trace-metric directions therefore use the same volume geometry in different
ways:

```text
apply G
    fix the trace
    solve only for its minimum-energy interior extension
    apply A and restrict the boundary covector

inverse_apply G^-1
    prescribe a trace-supported covector
    solve the full volume Riesz problem
    retain the boundary part of the primal solution
```

The solver-facing object remains a map on the trace space, even though its
realization uses auxiliary interior unknowns.

## 12. Geometry and preconditioning are related but not identical ideas

In finite dimensions, the update

$$
d
=
-G^{-1}r
$$

looks exactly like a preconditioned gradient step.

That observation is useful, but two viewpoints should not be collapsed completely.

A **solver metric** in the current reduced algorithms determines:

- how the derivative covector is converted into a primal gradient;
- how gradient and step norms are measured;
- what "steepest descent" means;
- which projection is valid for a constrained problem.

A **linear-solver preconditioner** may instead be an internal device used to
approximately solve some equation faster. It need not define the optimization norm or
the interpretation of a gradient.

The same SPD operator can sometimes play both roles, but the surrounding semantics
are different.

This distinction matters when reading the code. `MetricT` is visible to the
optimization algorithm. A preconditioner buried inside a CG solve is an
implementation detail of one inverse application.

## 13. Box constraints introduce a second metric-dependent operation

Suppose the control is constrained by lower and upper bounds.

At the continuous level one might write

$$
u_{a}(x)
\leq
u(x)
\leq
u_{b}(x).
$$

For a cellwise-constant control, the discrete admissible set has a particularly
simple coefficient representation:

$$
C
:=
\left\{
\mathbf u\in\mathbb R^{n}
:
\ell_{i}\leq u_{i}\leq r_{i}
\right\}.
$$

Testing feasibility is coefficientwise: every entry must lie in its interval.

Projection asks a different question. If a trial point $\mathbf x$ lies outside the
box, which feasible point should replace it?

Once an optimization metric $G$ has been chosen, the natural answer is the **nearest
feasible point in the $G$-norm**:

$$
P_{C}^{G}(\mathbf x)
:=
\mathop{\mathrm{argmin}}_{\mathbf z\in C}
\frac{1}{2}
\lVert \mathbf z-\mathbf x\rVert_{G}^{2}.
$$

Because

$$
\lVert \mathbf z-\mathbf x\rVert_{G}^{2}
=
(\mathbf z-\mathbf x)^{\mathsf T}
G
(\mathbf z-\mathbf x),
$$

this is a small convex quadratic optimization problem whenever $G$ is symmetric
positive definite and $C$ is a closed convex set.

The phrase "project onto the box" is therefore incomplete unless the geometry is
known. Euclidean projection, mass-metric projection, and projection in a coupled
Sobolev metric need not produce the same point.

### 13.1 The projection variational inequality

The minimization problem above has a useful first-order characterization that will
also explain projected-gradient stationarity later.

Let

$$
\mathbf z
:=
P_{C}^{G}(\mathbf x).
$$

Take any other feasible point $\mathbf v\in C$. Because $C$ is convex, the line
segment

$$
\mathbf z(t)
:=
\mathbf z+t(\mathbf v-\mathbf z),
\qquad
0\leq t\leq1,
$$

stays inside $C$.

Define the squared-distance objective along this feasible line:

$$
q(t)
:=
\frac{1}{2}
\lVert
\mathbf z(t)-\mathbf x
\rVert_{G}^{2}.
$$

Since $t=0$ is a minimum over feasible positive $t$, the one-sided derivative must
satisfy

$$
q'(0)
\geq0.
$$

Differentiate:

```math
\begin{aligned}
q'(0)
&=
(\mathbf v-\mathbf z)^{\mathsf T}
G
(\mathbf z-\mathbf x)
\\
&=
\left\langle
G(\mathbf z-\mathbf x),
\mathbf v-\mathbf z
\right\rangle.
\end{aligned}
```

Hence

$$
\left\langle
G(\mathbf z-\mathbf x),
\mathbf v-\mathbf z
\right\rangle
\geq0
\qquad
\text{for every }\mathbf v\in C.
$$

For an SPD metric and a closed convex set, this condition characterizes the metric
projection uniquely.

Geometrically, the metric covector

$$
G(\mathbf z-\mathbf x)
$$

points outward from the feasible set in the sense that every feasible displacement
from $\mathbf z$ has nonnegative pairing with it.

This variational inequality is the bridge between projection and constrained
first-order optimality.

## 14. Why clipping works for a positive diagonal metric

Let

$$
G
:=
\mathrm{diag}(g_{1},\ldots,g_{n}),
\qquad
g_{i}>0.
$$

The projection objective becomes

```math
\frac{1}{2}
(\mathbf z-\mathbf x)^{\mathsf T}
G
(\mathbf z-\mathbf x)
=
\frac{1}{2}
\sum_{i=1}^{n}
g_{i}(z_{i}-x_{i})^{2}.
```

There are no cross-terms involving $`z_{i}z_{j}`$ for $i\neq j$. The minimization
therefore splits into $n$ independent one-dimensional problems:

$$
\min_{\ell_{i}\leq z_{i}\leq r_{i}}
\frac{g_{i}}{2}(z_{i}-x_{i})^{2}.
$$

Because $`g_{i}>0`$, multiplying the one-dimensional objective by $`g_{i}`$ does not
change its minimizer.

Three cases remain:

```text
x_i < ell_i
    nearest feasible value is ell_i

ell_i <= x_i <= r_i
    x_i is already feasible

x_i > r_i
    nearest feasible value is r_i
```

Thus

$$
\left(P_{C}^{G}(\mathbf x)\right)_{i}
=
\min
\left(
r_{i},
\max(\ell_{i},x_{i})
\right).
$$

This is ordinary coefficientwise clipping.

The positive diagonal weights affect the metric distance but not which point in each
interval is nearest. This explains why the cellwise-constant $L^{2}$ case is so
convenient:

- each coefficient has a direct local bound;
- the mass matrix is positive diagonal;
- the metric projection onto the box separates coefficient by coefficient.

## 15. Why clipping fails for a coupled metric

Now take the same box but a coupled SPD metric:

$$
G
:=
\begin{bmatrix}
2&1\\
1&2
\end{bmatrix}.
$$

Project

$$
\mathbf x
:=
\begin{bmatrix}
2\\
0
\end{bmatrix}
$$

onto

$$
C=[0,1]^{2}.
$$

Euclidean clipping would give

$$
\mathbf z_{\mathrm{clip}}
=
\begin{bmatrix}
1\\
0
\end{bmatrix}.
$$

But the off-diagonal entries in $G$ couple the two coordinates. Changing $`z_{2}`$ can
reduce the metric distance caused by the error in $`z_{1}`$.

Consider the upper boundary $`z_{1}=1`$ and write

$$
\mathbf z(t)
:=
\begin{bmatrix}
1\\
t
\end{bmatrix},
\qquad
0\leq t\leq1.
$$

Then

$$
\mathbf z(t)-\mathbf x
=
\begin{bmatrix}
-1\\
t
\end{bmatrix}.
$$

Compute the quadratic form explicitly:

```math
\begin{aligned}
G(\mathbf z(t)-\mathbf x)
&=
\begin{bmatrix}
2&1\\
1&2
\end{bmatrix}
\begin{bmatrix}
-1\\
t
\end{bmatrix}
\\
&=
\begin{bmatrix}
-2+t\\
-1+2t
\end{bmatrix},
\end{aligned}
```

so

```math
\begin{aligned}
(\mathbf z(t)-\mathbf x)^{\mathsf T}
G(\mathbf z(t)-\mathbf x)
&=
(-1)(-2+t)
+t(-1+2t)
\\
&=
2-2t+2t^{2}.
\end{aligned}
```

Differentiate with respect to $t$:

$$
-2+4t=0,
$$

which gives

$$
t=\frac{1}{2}.
$$

The candidate metric projection is therefore

$$
\mathbf z_{\ast}
:=
\begin{bmatrix}
1\\
1/2
\end{bmatrix}.
$$

We can verify that this is the projection using the variational inequality from
Section 13.1.

First,

```math
G(\mathbf z_{\ast}-\mathbf x)
=
G
\begin{bmatrix}
-1\\
1/2
\end{bmatrix}
=
\begin{bmatrix}
-3/2\\
0
\end{bmatrix}.
```

For any feasible

$$
\mathbf v=(v_{1},v_{2})\in[0,1]^{2},
$$

we have

```math
\begin{aligned}
\left\langle
G(\mathbf z_{\ast}-\mathbf x),
\mathbf v-\mathbf z_{\ast}
\right\rangle
&=
-\frac{3}{2}(v_{1}-1)
\\
&\geq0,
\end{aligned}
```

because $`v_{1}\leq1`$.

Thus $`\mathbf z_{\ast}`$ satisfies the projection condition and, by strict convexity,
is the unique metric projection.

The result is

$$
P_{C}^{G}(\mathbf x)
=
\begin{bmatrix}
1\\
1/2
\end{bmatrix},
$$

not the clipped point.

The example exposes the role of the off-diagonal entries: the metric regards the two
coordinate errors as coupled, so the best feasible correction in one coordinate
depends on the other.

This is why **a box constraint does not by itself justify coefficientwise
clipping**. The metric realization must make the projection separable.

## 16. The current cellwise and facewise box realizations exploit separability

The deal.II implementation contains two closely related coefficientwise box
constraints:

- `CellwiseBoxConstraint` for cellwise-constant volume coefficients;
- `FacewiseBoxConstraint` for facewise-constant boundary coefficients.

In both cases, each optimization coefficient has a direct local interpretation and
the intended projection metric is a `MassMetric`.

Their constructors do not accept an arbitrary mass matrix and hope that clipping is
correct.

They ask the metric whether

```cpp
supports_coefficientwise_box_projection()
```

and `MassMetric` returns true only when:

- every off-diagonal matrix entry is zero;
- every diagonal entry is strictly positive.

That is exactly the condition needed by the separability argument above.

### 16.1 Why this excludes a standard continuous mass matrix

A consistent mass matrix for continuous nodal finite elements normally has nonzero
off-diagonal entries because neighboring basis functions overlap.

Therefore the current coefficientwise box classes do **not** treat clipping as a
valid projection in that metric.

This is not a limitation of box constraints mathematically. The metric projection
still exists.

What is missing is a realization of the coupled bound-constrained quadratic
projection problem.

The current classes intentionally implement the simpler diagonal case.

### 16.2 Facewise and cellwise constraints share algebra but not physical meaning

Both classes clamp one coefficient at a time.

The coefficients nevertheless represent different things:

```text
CellwiseBoxConstraint
    one coefficient per selected volume cell

FacewiseBoxConstraint
    one coefficient per selected boundary face
```

This is another example of a recurring theme in the manual: identical-looking
coefficient algebra does not imply identical numerical semantics.

## 17. `ConstraintT` separates feasibility from projection

At the contract level, a constraint exposes:

```cpp
is_feasible(primal)
supports_projection_in(metric)
project_in(primal, metric)
```

These methods answer different questions.

### Feasibility

$$
u\in C?
$$

For a coefficient box, this is a direct bounds check.

### Projection capability

Does this implementation know how to compute

$$
P_{C}^{G}
$$

for the supplied metric realization?

A constraint can recognize feasibility without supporting projection in every
possible metric.

### Projection

If the capability is supported, compute the actual projected primal value.

This decomposition becomes natural once projection has been defined as a
metric-dependent optimization problem rather than as generic "clamping".

## 18. Why the constraint remembers a metric realization, not only a metric name

`MetricT` exposes a human-readable `id()`.

That ID is useful for provenance and diagnostics.

It is not strong enough to establish that a particular projection algorithm is valid.

Imagine two constructed metric objects with:

- the same layout;
- the same string ID;
- different matrices.

If a box constraint was validated against one positive diagonal mass matrix, it
should not automatically accept the other merely because somebody reused the same
name.

For this reason `MetricT` also owns an opaque `MetricRealisationWitness`.

A constraint constructed for a particular metric stores that witness.

Later,

```cpp
supports_projection_in(metric)
```

requires both:

1. a compatible layout;
2. a matching realization witness.

The witness therefore means roughly:

> this is the metric realization whose numerical properties were checked when the
> projection capability was constructed.

It is intentionally stronger than a display label and much smaller than exposing the
whole matrix through the generic contract.

## 19. Projected-gradient stationarity is a projection fixed point

The projection machinery is not used only to repair infeasible trial points. It also
gives a natural first-order stationarity measure for a constrained problem.

Let

$$
C\subset U_{h}
$$

be a closed convex admissible set, and let

$$
r
:=
j'(u)
\in
U_{h}^{\ast}
$$

be the reduced derivative at a feasible point $u\in C$.

Choose the metric $G$ and define the metric gradient

$$
g
:=
G^{-1}r.
$$

Equivalently,

$$
Gg=r.
$$

We now compare two first-order statements.

### 19.1 First-order stationarity of the constrained objective

If $u$ is a local minimizer over the convex set $C$, then every feasible direction
from $u$ must have nonnegative first-order objective change.

Take any $v\in C$. Convexity implies that

$$
u(t)
:=
u+t(v-u)
$$

is feasible for $0\leq t\leq1$.

Define

$$
\varphi(t)
:=
j(u(t)).
$$

A local minimum at $t=0$ requires

$$
\varphi'(0)
\geq0.
$$

Using the chain rule,

$$
\varphi'(0)
=
j'(u)[v-u]
=
\langle r,v-u\rangle.
$$

Thus every constrained local minimizer satisfies the variational inequality

$$
\langle r,v-u\rangle
\geq0
\qquad
\text{for every }v\in C.
$$

This is a **first-order stationarity condition**. For a general nonlinear objective
it is necessary but not sufficient for local optimality. If $j$ is convex, the same
condition is also sufficient for global optimality over $C$.

### 19.2 The projection condition for $`u=P_{C}^{G}(u-g)`$

Now form the point

$$
x
:=
u-g.
$$

Ask whether projecting $x$ back onto $C$ returns $u$:

$$
u
\stackrel{?}{=}
P_{C}^{G}(u-g).
$$

From Section 13.1, a point $z\in C$ is the metric projection of $x$ if and only if

$$
\left\langle
G(z-x),
v-z
\right\rangle
\geq0
\qquad
\text{for every }v\in C.
$$

Substitute

$$
z=u,
\qquad
x=u-g.
$$

First simplify the displacement inside the metric:

```math
\begin{aligned}
z-x
&=
u-(u-g)\\
&=g.
\end{aligned}
```

Apply the metric:

$$
G(z-x)
=
Gg
=
r.
$$

Also,

$$
v-z
=
v-u.
$$

Therefore the projection variational inequality becomes

```math
\left\langle
r,
v-u
\right\rangle
\geq0
\qquad
\text{for every }v\in C.
```

But this is exactly the constrained first-order stationarity condition derived in
Section 19.1.

We have therefore shown the equivalence

```math
\boxed{
\begin{aligned}
\langle j'(u),v-u\rangle
&\geq0
\quad\text{for every }v\in C
\\
&\Longleftrightarrow
\\
u
&=
P_{C}^{G}(u-G^{-1}j'(u)).
\end{aligned}
}
```

The fixed point is not an unrelated projection trick. It is another form of the same
first-order variational inequality.

### 19.3 The projected-gradient residual

Define

$$
\widehat u
:=
P_{C}^{G}(u-g)
$$

and

$$
s
:=
\widehat u-u.
$$

Then

$$
s=0
$$

if and only if $u$ satisfies the constrained first-order stationarity condition.

This makes

$$
\lVert s\rVert_{G}
$$

a natural stationarity measure.

If $u$ is well inside the feasible set and $u-g\in C$, the projection does nothing:

$$
\widehat u=u-g,
$$

so

$$
s=-g.
$$

The projected measure then reduces to the ordinary gradient measure.

At an active bound, however, the unconstrained gradient need not vanish. The
projection removes the components that try to move outside the admissible set, and
$s=0$ can hold even while $g\neq0$.

### 19.4 Why the solver uses a unit metric-gradient step here

More generally, for any scalar $\tau>0$, the stationarity condition is also
equivalent to

$$
u
=
P_{C}^{G}(u-\tau g).
$$

Indeed, the projection variational inequality would contain

$$
G\left(u-(u-\tau g)\right)
=
\tau Gg
=
\tau r,
$$

and multiplying the variational inequality by the positive scalar $\tau$ does not
change its sign.

The current reduced solver uses the unit choice $\tau=1$ to define its projected
stationarity mapping. The actual line-search step length is a separate quantity and
is handled later when trial controls are built.

## 20. The reduced solver uses this projected update as a stationarity measure

The current reduced search code follows the preceding construction quite closely.

First, the direction policy computes the metric gradient

$$
g
=
G^{-1}j'(u)
$$

and the unconstrained steepest direction

$$
d=-g.
$$

For a constrained solve, the solver then forms the unit projected point

$$
\widehat u
=
P_{C}^{G}(u+d)
=
P_{C}^{G}(u-g).
$$

The projected update is

$$
s
=
\widehat u-u.
$$

The code measures

$$
\lVert s\rVert_{G}
$$

and uses this projected norm in place of the ordinary gradient norm for stopping.

It also computes

$$
j'(u)[s]
$$

as the relevant descent measure.

That is substantially more meaningful than testing $`\lVert g\rVert_{G}`$ alone:
at an active constrained optimum the unconstrained gradient need not vanish, while
the projected-gradient mapping does.

### 20.1 Trial steps are projected too

For a line-search step length $\alpha$, the trial builder starts from

$$
u+\alpha d
$$

and, when a constraint is present, projects the result back into the admissible set:

$$
u_{\mathrm{trial}}
=
P_{C}^{G}(u+\alpha d).
$$

The projection therefore serves two roles:

- it keeps line-search trial controls feasible;
- its unit-step fixed-point residual provides the constrained stationarity measure.

Chapter 6 explains how this interacts with Armijo
and other line searches.

## 21. How the reduced search code uses the metric

The metric enters the reduced search in several distinct places.

### Gradient construction

The helper `make_metric_gradient()` computes

$$
g
=
G^{-1}r
$$

using `metric.inverse_apply(reduced_derivative)`.

It then evaluates

$$
Gg
$$

with `metric.apply(gradient)` and forms

$$
\lVert g\rVert_{G}^{2}
=
\langle Gg,g\rangle
$$

using the primal/covector pairing.

The resulting object contains both:

- the primal gradient;
- its metric norm.

### Direction construction

The steepest-descent policy negates the gradient:

$$
d=-g.
$$

The descent measure is not a Euclidean dot product between two primal vectors. It is

$$
j'(u)[d]
=
\langle r,d\rangle.
$$

The solver therefore calls

```cpp
pair(reduced_derivative, direction)
```

using the derivative covector and primal direction.

### Step norms

After a trial step has been accepted, the actual update

$$
s
=
u_{\mathrm{new}}-u
$$

is measured through

$$
\lVert s\rVert_{G}
=
\sqrt{\langle Gs,s\rangle}.
$$

So the same metric controls the meaning of both gradient size and step size.

### Projected stationarity

For constrained steepest descent, the same metric is also passed to the projection,
as developed in Sections 19–20.

This gives one coherent geometry to the whole first-order search.

## 22. One metric inverse application can hide several linear solves

The reduced solver records metric work at the level of calls to the metric-gradient
protocol.

That level should not be confused with the internal cost of one metric realization.

For `MassMetric`,

```text
inverse_apply
    ≈ one mass-matrix CG solve
```

For `Hminus1Metric`,

```text
inverse_apply
    first mass solve
    + Laplace matrix action
    + second mass solve
```

For `TraceHhalfMetric`,

```text
inverse_apply
    one full volume H1 solve
    + trace restriction
```

Thus two metrics can present exactly the same generic operation

$$
G^{-1}r
$$

while having very different costs.

This is one reason the metric abstraction is useful: the optimizer can reason in
terms of geometry while the realization owns the required numerical machinery.

The detailed accounting of those operations belongs to the reduced-optimization and
verification chapters.

## 23. The metric implementations form a useful progression

The concrete realizations discussed in this chapter can be compared compactly.

| Realization | Riesz map $G$ | `apply` | `inverse_apply` |
| --- | --- | --- | --- |
| `DiagonalMetric` | positive diagonal $D$ | diagonal multiplication | diagonal division |
| `MassMetric` | supplied SPD matrix $M$ | sparse matrix action | CG solve with $M$ |
| `Hminus1Metric` | $M K^{-1}M$ | mass action + Laplace solve + mass action | mass solve + Laplace action + mass solve |
| `TraceHhalfMetric` | $`A_{BB}-A_{BI}A_{II}^{-1}A_{IB}`$ | minimum-extension solve + volume action + trace restriction | full volume solve + trace restriction |

The first is a backend-neutral dense/reference realization.

The remaining three are serial deal.II realizations.

The table is not a hierarchy of "better" metrics. Each operator represents a
different geometry appropriate to a different numerical space or modeling choice.

## 24. Metric IDs and realization witnesses answer different questions

A metric object carries two forms of identity.

### Display/provenance identity

`id()` returns a string.

This can answer questions such as:

> Which declared metric did this compiled problem report?

or

> Which metric name should appear in diagnostics?

### Executable realization identity

`realisation_witness()` returns an opaque token tied to the constructed metric
realization.

This can answer a stronger question:

> Is this the same numerical metric realization for which another capability was
> validated?

The box constraints use the second kind.

This distinction is subtle but useful beyond constraints. Numerical compatibility
sometimes depends on the actual constructed operator, not merely on a human-readable
name or matching coefficient dimensions.

## 25. A constraint is not itself a metric

It is also useful to separate two mathematical objects that meet in projected
optimization.

The constraint describes the admissible set

$$
C\subset U_{h}.
$$

The metric describes the geometry

$$
G:U_{h}\to U_{h}^{\ast}.
$$

Feasibility,

$$
u\in C,
$$

is meaningful without choosing $G$.

Metric projection,

$$
P_{C}^{G}(u),
$$

depends on both.

This explains the shape of `ConstraintT`:

```text
layout
is_feasible

supports_projection_in(metric)
project_in(primal, metric)
```

The interface does not force every constraint to implement projection in every
metric.

The optional `box_data_token()` visible in the current contract serves a different
purpose: it lets later complementarity/PDAS machinery retain access to the box data
associated with a constraint realization. Chapters 7 and 8 return to that through KKT and complementarity/PDAS rather than
mixing active-set semantics into the present projection discussion.

## 26. What happens if the metric and constraint disagree?

Several kinds of mismatch are possible.

### Different coefficient spaces

A control constraint built for one layout cannot sensibly project a vector from an
unrelated layout.

The contract checks layout compatibility.

### Same space, unsupported geometry

A coefficient box might be perfectly meaningful in a continuous control space while
the chosen mass metric is coupled.

Feasibility can still be checked, but coefficientwise clipping is not the metric
projection.

The current deal.II box implementations reject this combination.

### Same layout and name, different numerical realization

Even if two metrics have the same dimensions and display ID, a projection algorithm
validated for one constructed matrix should not silently accept the other.

The realization witness distinguishes them.

These cases show why "the dimensions match" is too weak a notion of compatibility
for some numerical operations.

## 27. Where the current constrained reduced solver stops

The current projection-based reduced search intentionally has a narrower capability
than the unconstrained direction framework.

When a `ConstraintT` is supplied, the implementation currently accepts projected
steepest descent with:

- Armijo line search;
- fixed step;
- weak Wolfe line search;
- strong Wolfe line search.

It does not simply combine the projection mechanism with every quasi-Newton,
conjugate-gradient, or Newton direction policy.

This is useful context when navigating the solver code. The presence of a generic
`ConstraintT` does not imply that every reduced algorithm has a mathematically
implemented constrained counterpart.

Chapter 6, **Reduced optimization methods**, explains the unconstrained and
projected algorithm families in more detail.

## 28. Following the metric path through the source

The source can now be read in a natural order.

### Generic metric and constraint contracts

Start with:

- [`include/nmopt/contract/metric_constraint.hpp`](../../../include/nmopt/contract/metric_constraint.hpp)

Read `MetricT` first:

```text
id
layout
realisation_witness
apply
inverse_apply
```

Then read `ConstraintT`:

```text
is_feasible
supports_projection_in
project_in
```

The dense `DiagonalMetric` and dense `CellwiseBoxConstraint` in the same header are
small reference implementations of the ideas developed above.

### Finite-element mass metric

Continue with:

- [`include/nmopt/dealii/mass_metric.hpp`](../../../include/nmopt/dealii/mass_metric.hpp)

Relate the two methods directly to

$$
Mv
$$

and

$$
M^{-1}r.
$$

Then inspect `supports_coefficientwise_box_projection()` and notice that it checks
for a positive diagonal matrix rather than trusting the name "mass metric".

### Negative-order metric

Read:

- [`include/nmopt/dealii/hminus1_metric.hpp`](../../../include/nmopt/dealii/hminus1_metric.hpp)

Follow the operator chain

$$
M K^{-1}M
$$

in `apply()` and

$$
M^{-1}K M^{-1}
$$

in `inverse_apply()`.

This file is particularly helpful for seeing that a metric action can itself contain
PDE-like solves.

### Trace metric

Read:

- [`include/nmopt/dealii/trace_hhalf_metric.hpp`](../../../include/nmopt/dealii/trace_hhalf_metric.hpp)

The comments at the top state the Schur-complement formula. Then follow:

```text
trace vector
→ minimum interior extension
→ volume H1 action
→ restrict to trace
```

for `apply_vector()`.

For the inverse, follow the trace-supported volume right-hand side through the full
volume solve.

### Coefficientwise box realizations

Read:

- [`include/nmopt/dealii/cellwise_box_constraint.hpp`](../../../include/nmopt/dealii/cellwise_box_constraint.hpp)
- [`include/nmopt/dealii/facewise_box_constraint.hpp`](../../../include/nmopt/dealii/facewise_box_constraint.hpp)

Both files become straightforward once the diagonal projection argument is clear.

### Solver consumption

Finally read:

- [`include/nmopt/solvers/reduced_search.hpp`](../../../include/nmopt/solvers/reduced_search.hpp)
- [`include/nmopt/solvers/reduced_gradient.hpp`](../../../include/nmopt/solvers/reduced_gradient.hpp)

Look first for:

```text
make_metric_gradient
make_steepest_descent_direction
metric_norm
evaluate_projected_gradient
```

Those helpers are the direct code counterparts of the derivations in this chapter.

## 29. A final correspondence

The main objects can now be placed side by side.

| Mathematical object | Coordinate realization | Contract/runtime role |
| --- | --- | --- |
| derivative $`r\in U_{h}^{\ast}`$ | covector coefficients $\mathbf r$ | `CovectorBlockT` |
| Riesz map $`G:U_{h}\to U_{h}^{\ast}`$ | SPD operator/matrix | `MetricT::apply` |
| metric gradient $g=G^{-1}r$ | solve/apply inverse operator | `MetricT::inverse_apply` |
| norm $`\lVert v\rVert_{G}^{2}`$ | $\mathbf v^{\mathsf T}G\mathbf v$ | `pair(metric.apply(v), v)` |
| admissible set $C$ | bounds or other constraint data | `ConstraintT` |
| metric projection $`P_{C}^{G}`$ | constrained quadratic minimization | `project_in(primal, metric)` |
| projected stationarity | $`P_{C}^{G}(u-g)-u`$ | projected-gradient update/norm |

The important movement is

```text
covector derivative
        │
        │ inverse Riesz map
        ▼
primal gradient
        │
        │ negate / choose search policy
        ▼
primal direction
        │
        │ project if constrained
        ▼
feasible primal update
```

That sequence is the bridge from differentiation to optimization.

## 30. What to carry into the formulation chapters

Part I of the concept manual has now established the common numerical language:

1. how a PDE-constrained problem becomes a finite-dimensional problem;
2. how spaces, coordinates, primal values, and covectors are represented;
3. how residual derivatives and transpose actions operate;
4. how a metric identifies derivatives with primal gradients and interacts with
   constraints.

We can therefore discuss formulations without repeatedly interrupting the derivation
to explain these foundations.

Chapter 5, **Reduced state–adjoint formulation**, focuses on one coherent
question:

> How do the state equation, objective derivative, adjoint solve, residual VJP, and
> control metric compose into a reduced optimization problem?

After that, **Reduced optimization methods** can focus on search directions,
globalization, stopping, Hessian actions, and work accounting.

## Read later

Useful existing documents for the ideas in this chapter are:

- [Theoretical formalism](../../design/mathematical-model.md), for the project's
  formal primal/dual, metric, constraint, and projection conventions.
- [Numerical realization](../overview/numerical-realization.md), for the architectural
  location of finite-element metric services.
- [Reduced optimization](../overview/reduced-optimization.md), for the high-level
  use of metrics in the current reduced solver path.
- [Chapter 5 elliptic control](../../studies/chapter-5/source-catalogue.md), for
  application families that motivate $L^{2}$, negative-order, and boundary control
  geometries.

The [Numerical realization](../overview/numerical-realization.md) overview returns to the assembly and
deal.II details behind the mass, Laplace, and trace operators. Chapters 7 and 8
develop box constraints again from the KKT/complementarity viewpoint rather than the
projection viewpoint used here.
