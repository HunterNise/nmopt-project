# Anatomy of a discrete PDE-constrained problem

## Why start here

The easiest way to become lost in `nmopt` is to begin from its C++ interfaces.

Names such as `ProblemSpec`, `BlockLayout`, `ExecutableModelT`, `MetricT`, or
`ReducedDTOT` are meaningful only after one understands the numerical problem they
are trying to represent. They are not the mathematics itself. They are pieces of a
software representation built around a particular view of discretized
PDE-constrained optimization.

This chapter therefore starts one level lower than the framework.

We will take a simple elliptic distributed-control problem, derive its weak and
discrete forms, and follow the objects that survive into an executable numerical
problem. Only near the end will we look at how those objects appear in the
repository.

The aim is not to give a complete finite-element or optimization course. The aim is
to make the transition

```text
PDE-constrained problem
        ↓
weak formulation
        ↓
finite-dimensional problem
        ↓
numerical operations
        ↓
nmopt representation
```

concrete enough that later chapters can discuss the framework without requiring the
reader to accept its abstractions on faith.

For the project-wide picture, read
[Project overview and architecture](../overview/project-architecture.md) first or in
parallel. For the project's formal mathematical convention, the authoritative
design note is
[Theoretical formalism](../../design/mathematical-model.md).

## 1. A problem we can carry through the whole chapter

A useful running example is the distributed elliptic control problem used by the
source Chapter 6 numerical examples:

```math
\begin{aligned}
\min_{y,u}\quad &
\frac{1}{2}\lVert y-z_{d}\rVert_{L^{2}(\Omega)}^{2}
+
\frac{\beta}{2}\lVert u\rVert_{L^{2}(\Omega)}^{2},
\\
-\Delta y &= f+u
\qquad \text{in }\Omega,
\\
y &= 0
\qquad \text{on }\partial\Omega.
\end{aligned}
```

The source example takes $\Omega=(0,1)^{2}$ and

$$
z_{d}(x)
:=
10x_{1}(1-x_{1})x_{2}(1-x_{2}).
$$

The project records this source problem as E6.5.1 and uses the B1 scenario family to
exercise a corresponding reduced optimization path. The exact reproduction status
and source omissions are documented separately in
[Chapter 6 numerical examples](../../studies/chapter-6/numerical-examples.md); here we
care about the mathematical structure, not reproduction fidelity.

At first glance the problem contains only two unknown fields:

- the **state** $y$, whose value is constrained by the PDE;
- the **control** $u$, which the optimization algorithm is free to vary.

The forcing $f$, desired state $`z_{d}`$, and regularization parameter $\beta$ are fixed
data.

That distinction already matters. If $f$ were instead an unknown coefficient to be
identified, then it would no longer be immutable data: it would become another
optimization variable, with its own space, derivatives, regularization, and possibly
constraints.

For now, however, the control is the only decision variable after the state equation
has been solved.

### 1.1 What the optimization is trying to balance

The first term,

$$
\frac{1}{2}\lVert y-z_{d}\rVert_{L^{2}(\Omega)}^{2},
$$

penalizes failure to match the desired state. If that were the only term, the
optimizer would use the control as aggressively as the PDE permits.

The second term,

$$
\frac{\beta}{2}\lVert u\rVert_{L^{2}(\Omega)}^{2},
$$

penalizes control effort. Increasing $\beta$ makes large controls more expensive;
decreasing it permits the optimizer to track $`z_{d}`$ more closely at the price of a
larger control.

Nothing about this statement yet tells us how many coefficients will represent
$y$ or $u$, which finite elements will be used, or how the integrals will be
evaluated. Those are discretization choices.

This separation between the continuous problem and its discrete realization is one
of the most important ideas to keep in mind when reading the repository.

### 1.2 One continuous problem can have several discrete realizations

Before comparing the source example with the project, it helps to attach a little
meaning to the finite-element vocabulary that will recur throughout the manual.

A **mesh** $`\mathcal T_{h}`$ partitions the physical domain into simpler cells – for
example triangles or quadrilaterals in two dimensions. The subscript $h$ conventionally
refers to a characteristic mesh scale: smaller $h$ means a finer spatial
discretization. A finite-element space then
specifies what kind of function is allowed on each cell and how neighboring cells are
coupled.

On a triangular mesh, the standard $`\mathbb P_{1}`$ space consists of functions that
are affine on each triangle. When the space is continuous, neighboring triangles
share nodal values, so the global field is continuous across their common edge.
The familiar degrees of freedom are therefore values associated with mesh vertices.

A **cellwise-constant** space makes a different choice: the function is constant on
each cell and may jump across cell boundaries. Its natural degrees of freedom are
cell values rather than shared nodal values. In finite-element notation this is a
piecewise degree-zero space; in the project it appears as one of the available
volume-control realizations.

These few sentences are only enough to read the present example. The
[Numerical realization](../overview/numerical-realization.md) overview places these
choices in the project architecture; this manual focuses on the distinctions that
matter to optimization rather than teaching general finite-element mechanics.

The book's E6.5.1 example reports continuous $`\mathbb P_{1}`$ state, adjoint, and
control spaces on a triangular mesh. The current project supports more than one
realization of the same distributed-control problem family.

The source Chapter 5 recipe behind B1 can request either:

- a cellwise-constant volume control; or
- a continuous homogeneous-Dirichlet control.

The second phrase means that the discrete control is represented by a continuous
finite-element field whose trace is constrained to zero on the boundary. That is a
project realization choice; it is not implied merely by the continuous statement
$u\in L^{2}(\Omega)$.

The default B1 scenario currently selects the first through the recipe defaults.
The simplex-mesh path, on the other hand, is restricted to the continuous-control
realization.

This is not a contradiction. It exposes the distinction between a mathematical
problem and a numerical representation of that problem:

```text
the mathematical problem family
        │
        ├─► one discrete realization
        │
        └─► another discrete realization
```

The PDE says “$u$ is a distributed control over the volume”. It does not by itself
say whether the discrete control is represented by one value per cell, a continuous
nodal field, or some other finite-dimensional space. Even the phrase “continuous
control” is incomplete until a mesh, local polynomial family, boundary treatment,
and coordinate basis have been selected.

Part III explains how the semantic/compiler path records and realizes those
choices. For now, the important point is simpler: **the finite-dimensional problem is
not determined until the spaces and coordinate representations have been chosen.**

When the manual refers to a mathematical finite-element family, it will use notation
such as $`\mathbb P_{1}`$ or $`\mathbb Q_{1}`$. Backticks are reserved for literal
implementation names such as `FE_Q`, `cellwise_constant`, or
`homogeneous_dirichlet_continuous`.

## 2. The strong PDE is not yet the numerical equation

The state equation is printed as

$$
-\Delta y=f+u.
$$

A finite-element code does not normally enforce this identity pointwise. It works
with a variational statement.

Let

$$
V := H_{0}^{1}(\Omega).
$$

Multiplying by a test function $v\in V$ and integrating by parts gives

$$
\int_{\Omega}\nabla y\cdot\nabla v
=
\int_{\Omega}fv
+
\int_{\Omega}uv.
$$

Equivalently, define the residual action

```math
\langle E(y,u),v\rangle
:=
\int_{\Omega}\nabla y\cdot\nabla v
-
\int_{\Omega}fv
-
\int_{\Omega}uv.
```

The state equation becomes

$$
\langle E(y,u),v\rangle = 0
\qquad
\text{for every }v\in V.
$$

This is the first conceptual shift that later framework interfaces are built around.

The object $E(y,u)$ is not best thought of as “a vector of PDE values”. It is a
linear functional on the test space. In functional-analytic notation,

$$
E(y,u)\in V^{\ast}.
$$

The residual becomes a coefficient vector only after we choose a finite-dimensional
test space and a basis.

### 2.1 Why this viewpoint is useful

The strong form is compact and physically recognizable, but it hides choices that a
numerical method must make.

For a more complicated PDE we may need to decide:

- which terms are integrated by parts;
- whether trial and test spaces coincide;
- how essential boundary values are represented;
- whether boundary contributions belong in the residual;
- which weak interpretation is used for low-regularity data;
- whether a control acts in the volume, through a boundary trace, or through a
  separate actuator map.

By treating the tested residual as the executable mathematical object, these choices
become explicit rather than implicit in a PDE string.

The same principle applies to the objective. What matters numerically is not merely
the phrase “$L^{2}$ tracking”, but the actual integral and how it is evaluated on the
chosen discrete spaces.

## 3. Choose finite-dimensional spaces

Now choose a finite-element state/test space

$$
V_{h}
:=
\mathrm{span}\{\varphi_{1},\ldots,\varphi_{n_{y}}\}
\subset V.
$$

Write

$$
y_{h}
:=
\sum_{j=1}^{n_{y}} y_{j}\varphi_{j}.
$$

The state is now represented by a coefficient vector

$$
\mathbf y
:=
(y_{1},\ldots,y_{n_{y}})^{\mathsf T}.
$$

The control requires its own discrete space. Let

$$
U_{h}
:=
\mathrm{span}\{\psi_{1},\ldots,\psi_{n_{u}}\}.
$$

Then

$$
u_{h}
:=
\sum_{k=1}^{n_{u}}u_{k}\psi_{k},
\qquad
\mathbf u
:=
(u_{1},\ldots,u_{n_{u}})^{\mathsf T}.
$$

Nothing requires $`U_{h}`$ to equal $`V_{h}`$.

If the control is continuous $`\mathbb P_{1}`$ on the same triangular mesh, the two spaces may be closely
related and may even have the same dimension after boundary treatment. If the control
is cellwise constant, however, $`U_{h}`$ has one basis function per control cell and
its dimension and topology are different.

This distinction will later motivate the framework's explicit treatment of spaces
and layouts. For the moment, it is enough to notice that

```text
state coefficients   y ∈ R^(n_y)
control coefficients u ∈ R^(n_u)
```

are not interchangeable merely because both are stored as vectors of floating-point
numbers.

## 4. The weak residual becomes a matrix equation

Test the weak equation with each state basis function $`\varphi_{i}`$.

The diffusion matrix is

$$
A_{ij}
:=
\int_{\Omega}
\nabla\varphi_{j}\cdot\nabla\varphi_{i}.
$$

The forcing vector is

$$
(\mathbf f)_{i}
:=
\int_{\Omega} f\varphi_{i}.
$$

The control-to-state coupling is

$$
B_{ik}
:=
\int_{\Omega}\psi_{k}\varphi_{i}.
$$

The finite-element state equation is therefore

$$
A\mathbf y
=
\mathbf f+B\mathbf u.
$$

Move everything to the left and define the discrete residual

$$
\mathbf E(\mathbf y,\mathbf u)
:=
A\mathbf y-B\mathbf u-\mathbf f.
$$

Then

$$
\mathbf E(\mathbf y,\mathbf u)=0.
$$

This is recognizably the finite-dimensional version of the weak residual.

The three matrices/vectors have different origins:

```text
A      state-to-test operator
B      control-to-test coupling
f      fixed test-space load
```

Even in this very simple linear example, the distinction between the state space,
control space, and test space already survives in the shapes of the matrices:

$$
A\in\mathbb R^{n_{y}\times n_{y}},
\qquad
B\in\mathbb R^{n_{y}\times n_{u}}.
$$

If $`U_{h}`$ is cellwise constant, $B$ is generally rectangular. If the control uses
the same continuous nodal basis as the state, $B$ becomes a mass-like square
coupling, but its *role* is still control-to-test coupling.

That role matters more than the accidental matrix shape.

## 5. The objective also becomes a finite-dimensional object

The discrete state is a coefficient vector, but the tracking term is still defined
through a physical-space integral:

$$
\frac{1}{2}\int_{\Omega}(y_{h}-z_{d})^{2}.
$$

Expanding $`y_{h}`$ gives

```math
\frac{1}{2}
\left(
\mathbf y^{\mathsf T}M_{y}\mathbf y
-
2\mathbf q^{\mathsf T}\mathbf y
+
c
\right),
```

where

```math
\begin{aligned}
(M_{y})_{ij}
&:=
\int_{\Omega}\varphi_{j}\varphi_{i},
\\
q_{i}
&:=
\int_{\Omega}z_{d}\varphi_{i},
\\
c
&:=
\int_{\Omega}z_{d}^{2}.
\end{aligned}
```

This form is worth pausing over because it exposes a detail that is easy to hide in
notation such as $`\lVert \mathbf y-\mathbf z_{d}\rVert_{M}^{2}`$.

The target $`z_{d}`$ need not be represented by a state-space coefficient vector at
all. It may remain an analytic function evaluated at quadrature points. In that case
the implementation naturally assembles the linear term $\mathbf q$ and scalar
constant $c$ from the function, rather than first interpolating $`z_{d}`$ into
$`V_{h}`$.

The current B1 path does exactly this kind of thing: the desired state is retained as
a runtime `dealii::Function`, with the B1 polynomial expression supplied by the
application layer. The semantic graph describes a state observation and a quadratic
tracking loss; the concrete compiler decides how that declared observation/loss is
realized numerically.

That is a small but important example of the distinction between:

```text
mathematical data
    z_d(x)

one possible coordinate representation
    vector of interpolated nodal values

the numerical action actually required
    integrate z_d against basis/test functions at quadrature
```

Those three things should not be identified automatically.

### 5.1 The control term introduces its own mass matrix

For the regularization term,

$$
\frac{\beta}{2}
\int_{\Omega}u_{h}^{2},
$$

define the control mass matrix

$$
(N_{u})_{k\ell}
:=
\int_{\Omega}\psi_{\ell}\psi_{k}.
$$

Then

$$
J_{h}(\mathbf y,\mathbf u)
:=
\frac{1}{2}
\left(
\mathbf y^{\mathsf T}M_{y}\mathbf y
-
2\mathbf q^{\mathsf T}\mathbf y
+
c
\right)
+
\frac{\beta}{2}
\mathbf u^{\mathsf T}N_{u}\mathbf u.
$$

The use of a separate symbol $`N_{u}`$ is intentional. If the state and control use
different finite-element spaces, there is no reason for the two mass matrices to be
the same.

For a cellwise-constant control, $`N_{u}`$ is particularly simple: on a standard
elementwise basis it is diagonal, with entries related to cell measures. For a
continuous nodal control, it is the usual sparse finite-element mass matrix.

The continuous objective did not change. The algebra needed to evaluate it did.

## 6. What the executable problem looks like now

After discretization, the optimization problem has become

```math
\begin{aligned}
\min_{\mathbf y,\mathbf u}\quad &
J_{h}(\mathbf y,\mathbf u),
\\
\text{subject to}\quad &
A\mathbf y-B\mathbf u-\mathbf f=0.
\end{aligned}
```

At this stage, one can forget the strong PDE temporarily and ask what operations a
numerical optimization method actually needs.

Given a candidate control $\mathbf u$, the state is determined by

$$
A\mathbf y=\mathbf f+B\mathbf u.
$$

So a reduced evaluation begins with a **state solve**.

Once $\mathbf y$ is available, the objective can be evaluated.

For a first-order optimization method, we also need the derivative of the reduced
objective with respect to $\mathbf u$. Computing that derivative efficiently leads
to the adjoint.

This gives the basic numerical lifecycle:

```text
candidate control u
      │
      ▼
assemble / form state right-hand side
      │
      ▼
solve A y = f + B u
      │
      ▼
evaluate J_h(y,u)
      │
      ├──────────────► objective value
      │
      └─ if derivative is required:
             │
             ▼
          adjoint solve
             │
             ▼
       reduced derivative
```

This is already close to the runtime structure described in the overview, but now
every box has a concrete mathematical meaning.

## 7. Why an adjoint appears

Because the state is constrained by the PDE, it changes whenever the control changes.
After eliminating the state conceptually, the reduced objective is

$$
j_{h}(\mathbf u)
:=
J_{h}(\mathbf y(\mathbf u),\mathbf u).
$$

Suppose we perturb the control by a direction $\delta\mathbf u$. The corresponding
state perturbation $\delta\mathbf y$ is not arbitrary: it must satisfy the
linearized state equation. For the present linear problem,

$$
A\delta\mathbf y
=
B\delta\mathbf u.
$$

The directional derivative of the reduced objective is therefore

```math
j_{h}'(\mathbf u)[\delta\mathbf u]
=
\mathbf d_{y}^{\mathsf T}\delta\mathbf y
+
\mathbf d_{u}^{\mathsf T}\delta\mathbf u,
```

where we introduce the coordinate vectors

```math
\begin{aligned}
\mathbf d_{y}
&:=
M_{y}\mathbf y-\mathbf q,\\
\mathbf d_{u}
&:=
\beta N_{u}\mathbf u.
\end{aligned}
```

These vectors represent the partial derivatives of $`J_{h}`$ with respect to the state
and control coordinates. Chapters 2 and 3 make the primal/dual interpretation
of such coefficient vectors more precise.

For a **single chosen direction** $\delta\mathbf u$, there is nothing inherently
wrong with the direct sensitivity approach. We can form $B\delta\mathbf u$, solve

$$
A\delta\mathbf y=B\delta\mathbf u,
$$

and substitute the resulting $\delta\mathbf y$ into the directional derivative.

The difficulty appears when an optimization method needs the *whole reduced
derivative*, meaning a representation that can act on arbitrary control directions.
A direct basis-by-basis construction would use the control basis vectors
$`\mathbf e_{1},\ldots,\mathbf e_{n_{u}}`$ and solve

$$
A\delta\mathbf y_{k}
=
B\mathbf e_{k},
\qquad
k=1,\ldots,n_{u}.
$$

That is one linearized state solve for each control degree of freedom. Writing

$$
A^{-1}B
$$

is compact algebra, but numerically it hides exactly this collection of solves (or an
equivalent multiple-right-hand-side computation). Even if a factorization of $A$ can
be reused, constructing or applying the full sensitivity map becomes unattractive
when $`n_{u}`$ is large.

The adjoint avoids building that control-to-state sensitivity map. It rearranges the
same chain-rule term so that the expensive inverse of the state operator is applied
once to a state-side quantity rather than separately to every control direction.

### 7.1 Introduce the adjoint once

Define the adjoint $\mathbf p$ as the solution of

$$
A^{\mathsf T}\mathbf p
=
\mathbf d_{y}.
$$

The transpose is important even though the Laplace stiffness matrix in this example is
symmetric. The same derivation must also work for nonsymmetric PDE operators, where
$A^{\mathsf T}$ is genuinely different from $A$.

Now start from the state-dependent part of the directional derivative. Because

$$
\delta\mathbf y
=
A^{-1}B\delta\mathbf u,
$$

we have

```math
\begin{aligned}
\mathbf d_{y}^{\mathsf T}\delta\mathbf y
&=
\mathbf d_{y}^{\mathsf T}A^{-1}B\delta\mathbf u\\
&=
\left(A^{-\mathsf T}\mathbf d_{y}\right)^{\mathsf T}
B\delta\mathbf u\\
&=
\mathbf p^{\mathsf T}B\delta\mathbf u\\
&=
\left(B^{\mathsf T}\mathbf p\right)^{\mathsf T}
\delta\mathbf u.
\end{aligned}
```

Substituting this back into the chain rule gives

```math
j_{h}'(\mathbf u)[\delta\mathbf u]
=
\left(
\mathbf d_{u}
+
B^{\mathsf T}\mathbf p
\right)^{\mathsf T}
\delta\mathbf u.
```

Hence the coefficient vector representing the reduced derivative is

$$
\mathbf r_{u}
:=
\mathbf d_{u}
+
B^{\mathsf T}\mathbf p
=
\beta N_{u}\mathbf u+B^{\mathsf T}\mathbf p.
$$

The sign is a consequence of the residual convention

$$
\mathbf E(\mathbf y,\mathbf u)
:=
A\mathbf y-B\mathbf u-\mathbf f.
$$

Since the control appears in the residual with derivative $-B$, the project's
Lagrangian convention produces the $+B^{\mathsf T}\mathbf p$ contribution above.

The computational saving is now visible. After the state solve, one adjoint solve
produces a vector $B^{\mathsf T}\mathbf p$ that represents the effect of state
dependence on *every* control direction. We no longer need one sensitivity solve per
control basis vector.

The cost of a reduced derivative therefore has the characteristic form:

```text
one state solve
+
objective/partial-derivative evaluation
+
one adjoint solve
+
one transpose control coupling
```

rather than a family of state-sensitivity solves indexed by the control degrees of
freedom.

Chapter 3, **Operators, derivatives, and adjoints**, returns to this derivation
without assuming a linear PDE. Chapter 5, **Reduced state–adjoint formulation**,
then shows how the same structure appears in `ReducedDTOT`.

For this chapter, the key point is that the executable problem is already more than a
function from a control vector to a scalar objective. It contains a state equation,
its linearized and transposed actions, and solve operations that can be composed to
differentiate the reduced objective efficiently.

## 8. The derivative is still not a search direction

The vector

$$
\mathbf r_{u}
=
\beta N_{u}\mathbf u+B^{\mathsf T}\mathbf p
$$

encodes the reduced derivative through

$$
j_{h}'(\mathbf u)[\delta\mathbf u]
=
\mathbf r_{u}^{\mathsf T}\delta\mathbf u.
$$

This object tells us how the objective changes in every control direction. It does
not yet say which *primal control direction* should be called the gradient.

That extra step is the finite-dimensional form of the **Riesz representation**.

Choose an inner product on the control space. The gradient $\mathbf g$ associated
with that inner product is defined by the requirement

$$
j_{h}'(\mathbf u)[\delta\mathbf u]
=
(\mathbf g,\delta\mathbf u)_{U}
\qquad
\text{for every }\delta\mathbf u.
$$

If the chosen geometry is the finite-element $L^{2}$ inner product, then

$$
(\mathbf g,\delta\mathbf u)_{L^{2}}
=
\mathbf g^{\mathsf T}N_{u}\delta\mathbf u.
$$

Comparing the two representations for every $\delta\mathbf u$ gives

$$
N_{u}^{\mathsf T}\mathbf g
=
\mathbf r_{u}.
$$

The mass matrix is symmetric, so this becomes

$$
N_{u}\mathbf g
=
\mathbf r_{u},
$$

or equivalently

$$
\mathbf g
=
N_{u}^{-1}\mathbf r_{u}.
$$

So the appearance of $`N_{u}^{-1}`$ is not an ad hoc rescaling of the derivative. The
mass matrix is the coordinate representation of the $L^{2}$ Riesz map

$$
R_{U}:U_{h}\longrightarrow U_{h}^{\ast},
$$

and solving the mass-matrix system applies $`R_{U}^{-1}`$ to the derivative covector.

If we instead chose the Euclidean coefficient inner product, the Riesz map would be
the identity matrix and the coefficient vector $`\mathbf r_{u}`$ would itself be the
gradient. If we chose an $H^{1}$-type or negative-order metric, a different operator
would appear.

This is also why it is useful to keep two ideas separate even when the same mass
matrix occurs in both:

- $`N_{u}`$ inside the objective contributes to the **regularization derivative**;
- $`N_{u}`$ used as a Riesz map defines the **optimization geometry**.

They happen to coincide for this common $L^{2}$ choice, but they answer different
questions.

A steepest-descent step in the selected geometry uses $-\mathbf g$.

We will defer the full discussion to **Metrics, gradients, and constraints**, where
the same idea will be developed for mass, Sobolev, negative-order, and trace metrics.
For now the chain is:

```text
state/control coordinates
        ↓
objective + residual
        ↓
adjoint calculation
        ↓
reduced derivative in U*
        ↓
Riesz map / chosen control metric
        ↓
primal gradient in U
        ↓
search algorithm
```

This is the mathematical reason the later framework representation distinguishes
derivative-like objects from primal vectors and places a metric operation between
differentiation and optimization.

## 9. The same problem seen at three levels

We can now place three descriptions of the same problem side by side.

### Continuous level

```math
\begin{aligned}
-\Delta y &= f+u,\\
J(y,u) &:=
\frac{1}{2}\lVert y-z_{d}\rVert_{L^{2}}^{2}
+
\frac{\beta}{2}\lVert u\rVert_{L^{2}}^{2}.
\end{aligned}
```

This level expresses the mathematical problem.

### Variational level

```math
\langle E(y,u),v\rangle
:=
\int_{\Omega}\nabla y\cdot\nabla v
-
\int_{\Omega}fv
-
\int_{\Omega}uv.
```

Together with the objective integrals, this level specifies the weak operations that
a finite-element method must realize.

### Discrete algebraic level

```math
\begin{aligned}
A\mathbf y-B\mathbf u-\mathbf f &= 0,\\
J_{h}(\mathbf y,\mathbf u)
&=
\frac{1}{2}
\left(
\mathbf y^{\mathsf T}M_{y}\mathbf y
-
2\mathbf q^{\mathsf T}\mathbf y
+
c
\right)
+
\frac{\beta}{2}
\mathbf u^{\mathsf T}N_{u}\mathbf u.
\end{aligned}
```

This level is close to what numerical kernels and solvers execute.

No one of these descriptions should replace the others.

The continuous form explains the physical/mathematical model. The weak form explains
the discretization contract. The algebraic form explains the actual numerical
operations.

A large part of `nmopt` exists to preserve enough information to move between these
levels without conflating them.

The correspondence can be summarized compactly:

| Concept | Continuous/model view | Variational/FE view | Discrete coefficient view |
| --- | --- | --- | --- |
| State | $y$ | $`y\in V:=H_{0}^{1}(\Omega)`$, later $`y_{h}\in V_{h}`$ | $`\mathbf y\in\mathbb R^{n_{y}}`$ |
| Control | $u$ | $u$ belongs to a chosen control space; later $`u_{h}\in U_{h}`$ | $`\mathbf u\in\mathbb R^{n_{u}}`$ |
| PDE constraint | $-\Delta y=f+u$ | $`\langle E(y,u),v\rangle:=\int_{\Omega}\nabla y\cdot\nabla v-\int_{\Omega}fv-\int_{\Omega}uv`$ | $\mathbf E(\mathbf y,\mathbf u):=A\mathbf y-B\mathbf u-\mathbf f$ |
| State tracking | $`\frac{1}{2}\lVert y-z_{d}\rVert_{L^{2}}^{2}`$ | $`\frac{1}{2}\int_{\Omega}(y_{h}-z_{d})^{2}`$ | $`\frac{1}{2}(\mathbf y^{\mathsf T}M_{y}\mathbf y-2\mathbf q^{\mathsf T}\mathbf y+c)`$ |
| Control penalty | $`\frac{\beta}{2}\lVert u\rVert_{L^{2}}^{2}`$ | $`\frac{\beta}{2}\int_{\Omega}u_{h}^{2}`$ | $`\frac{\beta}{2}\mathbf u^{\mathsf T}N_{u}\mathbf u`$ |
| Reduced derivative | $j'(u)\in U^{\ast}$ | $j'(u)[\delta u]$ acts on a control direction | $`\mathbf r_{u}^{\mathsf T}\delta\mathbf u`$ |

The table is deliberately not an identification of the columns. Each column exposes a
different aspect of the same problem. In particular, a coefficient vector is a
representation of a field or functional in a chosen basis, not the field or
functional itself.

## 10. From the derived problem to one `nmopt` realization

Sections 1–9 followed one mathematical/numerical derivation. From here the chapter
changes perspective.

The remaining sections use the current B1 path as an **orientation tour**: given the
objects we have just derived, where do the corresponding choices and constructions
appear in the repository? The purpose is not to explain the semantic model, compiler,
or deal.II realization completely. The rest of the manual develops those layers in turn.
Here we only want the transition from equations on the page to source code to stop
feeling arbitrary.

The discrete problem contains several distinct pieces of information:

- which fields are variables and which quantities are fixed data;
- which spaces those variables inhabit;
- which weak terms make up the PDE residual;
- which observation is used by the objective;
- which loss is applied to that observation;
- which metric is intended for the control;
- which finite-element spaces, quadrature rules, and meshes realize those declarations;
- which state and adjoint solves are available; and
- which optimization formulation will consume the resulting numerical operations.

A monolithic application could hard-code all of that in one class. `nmopt` instead
separates the declarations from the realization and the optimization algorithm.

The B1 path is a useful place to see the separation because the mathematical problem
is simple enough that the layers remain recognizable.

### 10.1 The application recipe names the problem family

The source Chapter 6 B1 scenario reuses the source Chapter 5 scalar distributed-control recipe.

At the application level, the recipe is selected by the ID

```text
chapter-5.scalar-diffusion-reaction-volume
```

and the associated parameters choose, among other things, the control
discretization:

```text
cellwise_constant
or
homogeneous_dirichlet_continuous
```

The recipe does not create a mesh or solve a PDE. It constructs a semantic problem
description.

That semantic description is produced by
`make_scalar_diffusion_reaction_problem(...)` for the cellwise case, or by the
corresponding continuous-control builder for the continuous case.

A useful way to read this path is:

```text
application recipe
    says which semantic problem family is being requested

semantic builder
    describes the mathematical/structural problem

compiler
    chooses and constructs a registered numerical realization
```

The names matter less at first than seeing how the three stages relate.

### 10.2 The semantic graph resembles the derivation we just made

The scalar diffusion-reaction semantic builder declares:

- a state variable;
- a control variable;
- a state test space;
- forcing and desired-state data;
- diffusion, reaction, and regularization data;
- residual terms for diffusion/reaction, volume forcing, and volume control;
- one state equation made from those residual terms;
- a state observation;
- a control observation;
- a quadratic tracking loss;
- a quadratic control-regularization loss;
- a control $L^{2}$ metric; and
- a reduced DTO formulation request.

That list is much easier to understand after deriving the discrete problem.

For example, the three residual terms correspond to the weak residual

```math
\langle E(y,u),v\rangle
=
\underbrace{\int_{\Omega}\nabla y\cdot\nabla v}_{\text{diffusion/reaction family}}
-
\underbrace{\int_{\Omega}fv}_{\text{volume source}}
-
\underbrace{\int_{\Omega}uv}_{\text{volume control}},
```

with reaction omitted in the pure Laplace case.

Likewise, “state observation + quadratic tracking loss” is a structured way of
describing the operation that eventually yields the tracking integral, while
“control observation + quadratic regularization loss” yields the regularization
term.

The semantic graph is therefore not a second mathematical theory layered on top of
the PDE. It is a software decomposition of the information needed to realize the
selected formulation.

Chapter 9, [Semantic problem model](09-semantic-problem-model.md), examines why the graph is decomposed
in exactly this way and how other problems force additional nodes.

## 11. The B1 scenario adds concrete choices that are not part of the PDE

A problem family is still not executable. B1 must also choose concrete data and
numerical policies.

The current B1 scenario records:

- the unit square $[0,1]^{2}$;
- diffusion coefficient $1$;
- reaction coefficient $0$;
- the desired-state polynomial from the Chapter 6 source;
- a regularization sweep including $10^{-1}$, $10^{-2}$, $10^{-3}$, and
  $10^{-6}$;
- an assembled reduced-DTO product;
- a state and adjoint solve policy;
- steepest descent or limited-memory BFGS;
- Armijo globalization; and
- experiment/provenance information.

It also chooses a concrete forcing definition.

This last point is worth discussing because it shows why the application/scenario
layer exists.

### 11.1 The source problem does not fully determine the experiment

The Chapter 6 source record states the equation with a forcing $f$, but the available
source material does not specialize $f$ numerically.

The project therefore cannot recover a unique experiment from the printed PDE alone.
The B1 scenario currently defaults to a *manufactured zero forcing* and records that
choice and its provenance explicitly.

So the relationship is not

```text
book equation
    = exact executable scenario
```

but rather

```text
book problem family
        │
        ├─ source facts
        │
        ├─ unresolved source details
        │
        └─ declared project choices
                  │
                  ▼
          executable scenario
```

This distinction becomes important whenever one is reading benchmark code. A
numerical result can only be interpreted correctly if one knows which details come
from the source and which were supplied by the project.

The manual uses B1 here because it is concrete, but the same principle applies to
non-benchmark applications: an abstract PDE model rarely determines every mesh,
quadrature, solver, and data-representation choice needed for execution.

## 12. Runtime data and discretization meet only at compilation

B1 keeps the analytic forcing and desired-state functions alive as runtime
`dealii::Function` objects.

The application adapter packages them as `DealiiDataBindings`, together with numeric
data such as diffusion, reaction, and regularization values.

Separately, a compilation session owns or provides the mesh, and a discretization
policy records choices such as the state polynomial degree and linear-solve
tolerances.

Schematically:

```text
semantic problem
    "there is a forcing"
    "there is a desired state"
    "there is a volume control"
          │
          │
runtime data                    discretization/session
f(x), z_d(x), beta              mesh, FE degree, solve policy
          │                             │
          └────────────┬────────────────┘
                       ▼
                    compiler
                       │
                       ▼
             concrete numerical problem
```

This explains a design choice that can otherwise look needlessly indirect.

The semantic graph does not carry a `dealii::Function<2>` because that would make the
problem description itself deal.II-specific. But the compiler cannot assemble
$\mathbf f$ or $\mathbf q$ without the actual functions. Runtime bindings provide the
missing concrete data when numerical realization happens.

The same separation allows the same semantic description to be paired with a
different mesh or compatible data realization without rewriting the mathematical
graph.

## 13. What the compiler must eventually construct

For the simple problem in this chapter, the compiler's job can be understood in very
ordinary finite-element terms.

It must create or identify enough numerical machinery to realize operations such as:

```text
state coordinates y
control coordinates u

A y
B u
f

J_h(y,u)
D_y J_h
D_u J_h

solve A y = f + B u
solve A^T p = D_y J_h

apply the selected control metric
```

The actual implementation contains more bookkeeping because it must also preserve
space identities, provenance, solve reports, native views, optional constraints, and
the ability to package several formulation products.

But the core numerical meaning is the one we have just derived.

This is an important orientation trick when reading a large compiler header: do not
start by trying to memorize its records and helper types. Ask instead:

> Which mathematical or finite-element operation from the discrete problem is this
> code constructing?

That question usually collapses a large amount of implementation detail into a
recognizable role.

## 14. From concrete matrices to operation interfaces

The derivation above used matrices because they make the finite-dimensional
structure easy to see:

$$
A\mathbf y-B\mathbf u-\mathbf f.
$$

The solver-facing interfaces do not require every implementation to expose those
matrices directly.

Instead, later layers work with operations:

```text
residual at (y,u)
linearized residual action
transpose residual action
objective value
objective derivative

state solve
adjoint solve
metric application / inverse application
```

For an assembled deal.II realization, these operations may internally use sparse
matrices $A$, $B$, $`M_{y}`$, and $`N_{u}`$.

For an external PDE application, they may call methods on an existing application
object.

A future matrix-free realization could provide the same mathematical actions without
materializing every matrix in the same form.

This is the bridge from the algebra in this chapter to the numerical contract layer
described in the architecture overview.

Chapters 2 and 3 slow down at this boundary:

- **Spaces, coordinates, and duality** explains what information must accompany
  the coefficient vectors.
- **Operators, derivatives, and adjoints** explains how the matrix formulas
  generalize into residual/JVP/VJP actions.

## 15. A first source tour, following the problem rather than the directory tree

With the running example in mind, the relevant source path becomes much easier to
navigate.

### The source problem and project scenario

Start with:

- [`docs/studies/chapter-6/numerical-examples.md`](../../studies/chapter-6/numerical-examples.md)
  for the E6.5.1 source record;
- [`include/nmopt/application/chapter6.hpp`](../../../include/nmopt/application/chapter6.hpp)
  for B1 problem/scenario options and source-specific choices;
- [`include/nmopt/application/dealii/chapter6_b1.hpp`](../../../include/nmopt/application/dealii/chapter6_b1.hpp)
  for the deal.II execution adapter.

This tells you what experiment is being requested.

### The reusable problem family

Then read:

- [`include/nmopt/application/chapter5.hpp`](../../../include/nmopt/application/chapter5.hpp)
  for the distributed-control recipe;
- [`include/nmopt/semantic/v1/problem_library.hpp`](../../../include/nmopt/semantic/v1/problem_library.hpp)
  for the actual semantic graph produced by that recipe.

This tells you how the mathematical problem family is represented before numerical
realization.

### Numerical realization

Only after that is it useful to enter:

- [`include/nmopt/compiler/v1/dealii_compiler.hpp`](../../../include/nmopt/compiler/v1/dealii_compiler.hpp)
  and the smaller deal.II/compiler helper headers;
- [`include/nmopt/dealii/`](../../../include/nmopt/dealii/)
  for reusable numerical services.

This is where abstract roles turn into finite-element coordinates, matrices,
quadrature actions, metrics, constraints, and solve policies.

### Formulation and optimization

Finally, follow the compiled numerical operations into:

- [`include/nmopt/contract/reduced_dto.hpp`](../../../include/nmopt/contract/reduced_dto.hpp)
  for the reduced state–adjoint formulation;
- [`include/nmopt/solvers/`](../../../include/nmopt/solvers/)
  for reduced optimization algorithms.

This is where the state solve, adjoint solve, reduced derivative, metric, and search
policy are composed into an optimization run.

The useful mental route is therefore:

```text
source problem
    ↓
scenario choices
    ↓
semantic problem family
    ↓
compiler / FE realization
    ↓
formulation
    ↓
optimization algorithm
```

The repository directory tree is not itself the explanation; it is the place where
the successive stages live.

## 16. What changes when the problem becomes less simple

The distributed Laplace example is useful precisely because many complications are
absent.

Changing the problem reveals why later parts of the framework exist.

### Boundary control

If the control acts on a boundary $`\Gamma_{c}`$ rather than in the volume, then

$$
U=L^{2}(\Gamma_{c})
$$

or perhaps a trace space. The coupling matrix $B$ is assembled from boundary
integrals rather than volume integrals, and the appropriate control metric may be a
boundary mass matrix or a fractional trace metric.

The algebraic pattern survives, but the realization changes substantially.

### Parameter identification

If the decision variable is a coefficient inside the PDE operator, then $A$ itself
depends on the parameter:

$$
A(\mathbf m)\mathbf y=\mathbf f.
$$

The residual derivative with respect to $\mathbf m$ is no longer a fixed coupling
matrix $-B$. The JVP/VJP viewpoint becomes more useful than hard-coding the linear
matrix picture.

### Nonlinear PDEs

For a nonlinear residual,

$$
\mathbf E(\mathbf y,\mathbf u)=0,
$$

the state matrix in the sensitivity and adjoint equations is the Jacobian

$$
D_{\mathbf y}\mathbf E(\mathbf y,\mathbf u).
$$

The structure of the adjoint argument is unchanged even though the operators now
depend on the current point.

### Constraints

If

$$
u_{a}\leq u\leq u_{b},
$$

then a reduced derivative alone does not characterize the optimum. Projection,
variational inequalities, KKT conditions, or complementarity enter the picture.

This eventually leads to the project's box-constraint and PDAS machinery.

### All-at-once formulations

The reduced method eliminates the state by repeatedly solving the PDE.

An all-at-once method instead treats state, adjoint/multiplier, and control variables
as one coupled system. The KKT and supplied-OTD paths therefore expose a different
numerical product even though they originate from related optimality conditions.

Chapters 7 and 8 develop those alternatives from the equations
rather than presenting them as a collection of C++ types.

## 17. What to carry forward

The important result of this chapter is not a list of framework classes.

It is the following picture.

A PDE-constrained optimization problem begins with fields and equations, but after a
weak formulation and discretization it becomes a structured collection of
finite-dimensional spaces and operations:

```text
state coefficients y          control coefficients u
        │                              │
        └─────────────┬────────────────┘
                      ▼
               discrete residual
               E_h(y,u)
                      │
            ┌─────────┴─────────┐
            ▼                   ▼
        state solve        derivative actions
            │                   │
            ▼                   ▼
          state y            adjoint solve
            │                   │
            └─────────┬─────────┘
                      ▼
                 objective /
               reduced derivative
                      │
                      ▼
              optimization geometry
                      │
                      ▼
                search algorithm
```

Several of those objects may be stored as vectors of doubles, but they do not
represent the same space or the same mathematical role.

That observation creates the question taken up in Chapter 2:

> Once a PDE has been discretized, what information about spaces, coordinates, and
> duality must survive at runtime so that these numerical objects can be composed
> correctly?

That is where the framework-specific ideas of layouts, primal values, covectors, and
pairings will finally enter.

## Read later

The following documents deepen particular parts of this chapter without replacing its
narrative:

- [Theoretical formalism](../../design/mathematical-model.md) gives the project's
  normative abstract formulation, derivative, adjoint, and metric conventions.
- [Chapter 5 elliptic control](../../studies/chapter-5/source-catalogue.md) develops the
  mathematical application families from the source text.
- [Chapter 6 numerical methods](../../studies/chapter-6/numerical-methods.md) records
  the source numerical-method context.
- [Chapter 6 numerical examples](../../studies/chapter-6/numerical-examples.md)
  distinguishes source facts from project choices for B1/B2 and the other examples.
- [Numerical realization](../overview/numerical-realization.md) returns to the same
  subject from the project-wide architectural viewpoint.
- [Reduced optimization](../overview/reduced-optimization.md) gives the high-level
  runtime path developed in detail by the formulation and solver chapters.
