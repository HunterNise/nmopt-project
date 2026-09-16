# Operators, derivatives, and adjoints

## From a residual value to a residual operator

The first two concept chapters built the objects that meet at the executable-model
boundary.

We now have:

- a discrete variable space $X_{h}$, possibly a product such as
  $Y_{h}\times U_{h}$;
- a test space $Z_{h}$;
- primal coordinate representations for elements of those spaces;
- covector coordinate representations for elements of their duals; and
- a residual

$$
E_{h}(x)\in Z_{h}^{\ast}.
$$

The remaining question is how the residual changes when the variables change.

That question is the source of the two operations that later appear in the code as
`residual_jvp()` and `residual_vjp()`.

This chapter develops those operations before introducing the acronyms. We will begin
with the distributed-control problem from
[Anatomy of a discrete PDE-constrained problem](discrete-problem-anatomy.md), then
move to a nonlinear example where the distinction between an operator action and an
assembled matrix becomes more important.

The main thread is:

```text
residual E(x)
      │
      │ perturb x in direction δx
      ▼
linearized residual E'(x) δx
      │
      │ transpose with respect to dual pairings
      ▼
adjoint action E'(x)* p
      │
      ▼
JVP / VJP executable operations
```

No optimization metric is needed for this story. The transpose is defined by the
pairings between spaces and their duals, not by the choice of gradient geometry.

## 1. Begin with the linear distributed-control residual

For the running example, let

$$
X_{h}
:=
Y_{h}\times U_{h},
$$

and write a point as

$$
x
:=
(y,u).
$$

The discrete residual is

$$
E_{h}(y,u)
:=
Ay-Bu-f
\in
Z_{h}^{\ast}.
$$

For the simple Galerkin Laplace example, $Y_{h}$ and $Z_{h}$ may use the same finite
element space, but we keep their roles distinct.

Now perturb the point by

$$
\delta x
:=
(\delta y,\delta u).
$$

Because the residual is affine,

```math
\begin{aligned}
E_{h}(y+\varepsilon\delta y,u+\varepsilon\delta u)
&=
A(y+\varepsilon\delta y)
-
B(u+\varepsilon\delta u)
-
f\\
&=
E_{h}(y,u)
+
\varepsilon
\left(
A\delta y-B\delta u
\right).
\end{aligned}
```

The linearized residual is therefore

$$
E_{h}'(y,u)[\delta y,\delta u]
=
A\delta y-B\delta u.
$$

In this particular problem the derivative does not depend on $(y,u)$ because the
residual is linear in the variables.

That will not remain true in general.

### 1.1 The derivative is a map between spaces

The important object is not the formula alone but its domain and codomain:

$$
E_{h}'(x):
X_{h}
\longrightarrow
Z_{h}^{\ast}.
$$

At a fixed point $x$, the derivative is a **linear map in the perturbation**
$\delta x$.

The residual itself may be nonlinear in $x$.

This distinction is basic but worth stating explicitly:

```text
E_h
    generally nonlinear map
    X_h → Z_h*

E_h'(x)
    linear map obtained after fixing x
    X_h → Z_h*
```

A numerical interface that exposes the derivative action therefore needs two inputs:

1. the point $x$ at which the derivative is evaluated;
2. the tangent direction $\delta x$ on which the linearized map acts.

That is the mathematical shape behind a Jacobian-vector product.

## 2. Directional derivatives before matrices

In finite dimensions it is tempting to identify the derivative immediately with a
Jacobian matrix.

Suppose coordinates have been chosen:

$$
\mathbf x\in\mathbb R^{n},
\qquad
\mathbf E(\mathbf x)\in\mathbb R^{m}.
$$

Then the Jacobian matrix

$$
J_{E}(\mathbf x)
\in
\mathbb R^{m\times n}
$$

satisfies

$$
E_{h}'(x)[\delta x]
\quad\longleftrightarrow\quad
J_{E}(\mathbf x)\delta\mathbf x.
$$

This matrix picture is useful, but it is not the most general executable
representation.

A code can compute the action

$$
\delta\mathbf x
\longmapsto
J_{E}(\mathbf x)\delta\mathbf x
$$

without ever assembling or storing the whole matrix.

That can happen because:

- the action has an analytic matrix-free implementation;
- automatic differentiation propagates the tangent through residual evaluation;
- an application already exposes a linearized operator;
- the matrix exists conceptually but is too expensive or unnecessary to materialize.

For that reason, nmopt treats the **action** of the derivative as the primary
numerical capability.

The matrix is one possible realization.

## 3. The Jacobian-vector product

The phrase **Jacobian-vector product**, abbreviated JVP, refers to

$$
E_{h}'(x)[\delta x].
$$

In coordinates,

$$
\mathbf r_{\mathrm{jvp}}
:=
J_{E}(\mathbf x)\delta\mathbf x.
$$

The result belongs to the residual codomain:

$$
E_{h}'(x)[\delta x]
\in
Z_{h}^{\ast}.
$$

That codomain is easy to lose sight of when both input and output are ordinary native
vectors.

For the distributed-control residual,

$$
E_{h}'(y,u)[\delta y,\delta u]
=
A\delta y-B\delta u.
$$

So the JVP combines a state perturbation and a control perturbation into a **residual
covector**.

### 3.1 Product-space structure is already visible

With

$$
X_{h}
=
Y_{h}\times U_{h},
$$

the derivative can be viewed blockwise as

$$
E_{h}'(y,u)
=
\begin{bmatrix}
D_{y}E_{h} & D_{u}E_{h}
\end{bmatrix}.
$$

For the linear example,

$$
D_{y}E_{h}=A,
\qquad
D_{u}E_{h}=-B.
$$

Therefore

```math
E_{h}'(y,u)[\delta y,\delta u]
=
D_{y}E_{h} \delta y+D_{u}E_{h} \delta u.
```

This block form becomes important in reduced optimization because the state
linearization and control linearization play different roles:

- $D_{y}E_{h}$ defines the state sensitivity and adjoint systems;
- $D_{u}E_{h}$ carries the effect of the control into the residual.

The generic executable model does not hard-code this state/control split. It exposes
the full derivative action on the variable product. The reduced formulation later
selects the relevant blocks.

## 4. A nonlinear example changes the point, not the derivative concept

Consider a semilinear state equation

$$
-\Delta y+c y^{3}=f+u,
$$

with $c>0$.

Its weak residual is

```math
\langle E(y,u),v\rangle
:=
\int_{\Omega}\nabla y\cdot\nabla v
+
\int_{\Omega}c y^{3}v
-
\int_{\Omega}fv
-
\int_{\Omega}uv.
```

Perturb $y$ by $\delta y$ and $u$ by $\delta u$. Differentiating gives

```math
\langle E'(y,u)[\delta y,\delta u],v\rangle
=
\int_{\Omega}\nabla\delta y\cdot\nabla v
+
\int_{\Omega}3c y^{2}\delta y v
-
\int_{\Omega}\delta u v.
```

Now the linearized state operator depends on the current state $y$ through the
coefficient $3cy^{2}$.

After finite-element discretization, the Jacobian is no longer one fixed matrix
$[A,-B]$. Its state block changes with the evaluation point.

Yet the executable operation has exactly the same abstract form:

$$
(x,\delta x)
\longmapsto
E_{h}'(x)[\delta x].
$$

This is why the interface accepts both `variables` and `variable_tangent`, even
though a linear test problem might make the first input appear redundant.

## 5. How can we tell whether a JVP is correct?

A derivative implementation should agree with the change in the underlying residual.

For a sufficiently smooth residual,

```math
E_{h}(x+\varepsilon\delta x)
=
E_{h}(x)
+
\varepsilon E_{h}'(x)[\delta x]
+
\mathcal O(\varepsilon^{2}).
```

Rearranging gives the first-order finite-difference check

```math
\frac{
E_{h}(x+\varepsilon\delta x)-E_{h}(x)
}{\varepsilon}
\longrightarrow
E_{h}'(x)[\delta x]
```

as $\varepsilon\to0$.

In coefficient form, one can compare the implemented JVP with

$$
\frac{
\mathbf E(\mathbf x+\varepsilon\delta\mathbf x)
-
\mathbf E(\mathbf x)
}{\varepsilon}.
$$

The point of the test is not that finite differences are a better way to compute the
derivative. They are usually worse. They provide an independent numerical check that
the residual and linearized action agree.

### 5.1 A Taylor remainder gives a stronger diagnostic

Instead of dividing by $\varepsilon$, consider the remainder

$$
R(\varepsilon)
:=
E_{h}(x+\varepsilon\delta x)
-
E_{h}(x)
-
\varepsilon E_{h}'(x)[\delta x].
$$

For a smooth residual,

$$
\lVert R(\varepsilon)\rVert
=
\mathcal O(\varepsilon^{2}).
$$

So halving $\varepsilon$ should reduce the remainder by roughly a factor of four
until floating-point or solver errors dominate.

The later verification chapter will discuss such tests systematically. Here they
serve a conceptual purpose: the JVP is the **first-order change predicted by the
residual model**.

## 6. The transpose action is defined by a pairing identity

The JVP moves a perturbation forward through the derivative:

$$
\delta x
\longmapsto
E_{h}'(x)[\delta x]
\in Z_{h}^{\ast}.
$$

Adjoint methods need the derivative in the opposite direction.

Take a test-space element

$$
p\in Z_{h}.
$$

Because the JVP is a residual covector, it can act on $p$ through the dual pairing:

$$
\left\langle
E_{h}'(x)[\delta x],
p
\right\rangle_{Z_{h}^{\ast},Z_{h}}.
$$

For fixed $x$ and $p$, this expression is linear in $\delta x$. It therefore defines
a covector on the variable space.

We call that covector the transpose or adjoint action

$$
E_{h}'(x)^{\ast}p
\in X_{h}^{\ast},
$$

defined by the identity

```math
\left\langle
E_{h}'(x)[\delta x],
p
\right\rangle_{Z_{h}^{\ast},Z_{h}}
=
\left\langle
E_{h}'(x)^{\ast}p,
\delta x
\right\rangle_{X_{h}^{\ast},X_{h}}
```

for every $\delta x\in X_{h}$.

This identity is the central definition.

It tells us:

- what space the seed $p$ belongs to;
- what space the result belongs to;
- what "transpose" means even if no matrix has been assembled.

### 6.1 This transpose does not require an optimization metric

The notation $E'(x)^{\ast}$ can suggest a Hilbert-space adjoint defined through inner
products. That is not the construction being used at this boundary.

Here

$$
E_{h}'(x):
X_{h}\to Z_{h}^{\ast}
$$

already lands in a dual space. Its transpose is defined by the natural dual pairings

$$
Z_{h}^{\ast}\times Z_{h}
\quad\text{and}\quad
X_{h}^{\ast}\times X_{h}.
$$

No Riesz map is needed to define

$$
E_{h}'(x)^{\ast}:Z_{h}\to X_{h}^{\ast}.
$$

A metric enters only if we later want to identify a covector in $X_{h}^{\ast}$ with
a primal vector in $X_{h}$.

That is why the executable model can expose a VJP without knowing the optimization
metric.

## 7. The matrix transpose reappears after coordinates are chosen

Let the JVP have coordinate representation

$$
J_{E}(\mathbf x)\delta\mathbf x.
$$

Let $\mathbf p$ contain the primal coordinates of the test seed.

Using the coefficient pairings from
[Spaces, coordinates, and duality](spaces-coordinates-and-duality.md),

```math
\begin{aligned}
\left\langle
E_{h}'(x)[\delta x],
p
\right\rangle
&=
\left(
J_{E}(\mathbf x)\delta\mathbf x
\right)^{\mathsf T}
\mathbf p\\
&=
\delta\mathbf x^{\mathsf T}
J_{E}(\mathbf x)^{\mathsf T}
\mathbf p.
\end{aligned}
```

Therefore the covector coordinates of the transpose action are

$$
\mathbf r_{\mathrm{vjp}}
:=
J_{E}(\mathbf x)^{\mathsf T}\mathbf p.
$$

This is the **vector-Jacobian product** or VJP.

The terminology can be remembered as follows:

```text
JVP
    Jacobian acts on a variable tangent
    J δx

VJP
    test seed acts back through the Jacobian
    J^T p
```

The implementation need not literally construct $J$ or $J^{\mathsf T}$, but the
pairing identity must hold.

## 8. The linear distributed-control VJP

Return to

$$
E_{h}'(y,u)
=
\begin{bmatrix}
A & -B
\end{bmatrix}.
$$

Let $p\in Z_{h}$.

Then

```math
E_{h}'(y,u)^{\ast}p
=
\begin{bmatrix}
A^{\mathsf T}p\\
-B^{\mathsf T}p
\end{bmatrix}
\in
Y_{h}^{\ast}\times U_{h}^{\ast}.
```

The two components have immediate interpretations:

```text
A^T p
    state-space covector

-B^T p
    control-space covector
```

Check the defining identity:

```math
\begin{aligned}
\left(
A\delta y-B\delta u
\right)^{\mathsf T}p
&=
\delta y^{\mathsf T}A^{\mathsf T}p
-
\delta u^{\mathsf T}B^{\mathsf T}p\\
&=
\left\langle
\begin{bmatrix}
A^{\mathsf T}p\\
-B^{\mathsf T}p
\end{bmatrix},
\begin{bmatrix}
\delta y\\
\delta u
\end{bmatrix}
\right\rangle.
\end{aligned}
```

This is exactly the transpose structure that appeared in the reduced derivative in
the first chapter.

The adjoint state $p$ is chosen so that the state component
$A^{\mathsf T}p$ matches the state part of the objective derivative. Once that is
true, the control component $-B^{\mathsf T}p$ supplies the state-dependence
contribution to the reduced derivative.

### 8.1 Symmetry can hide the transpose

For the simple Laplace operator with a standard conforming Galerkin discretization,

$$
A^{\mathsf T}=A.
$$

It is therefore easy to write an implementation that uses `vmult` in both the state
and adjoint paths and never notice that an adjoint operation has been assumed.

That shortcut does not survive a nonsymmetric operator.

For example, adding transport can produce a discrete state operator

$$
K=A+C,
$$

where the transport contribution $C$ is generally nonsymmetric.

The state linearization uses

$$
K\delta y,
$$

while the transpose action uses

$$
K^{\mathsf T}p.
$$

A code path that silently substitutes $Kp$ for $K^{\mathsf T}p$ may still work on a
symmetric test case and fail on the first transport-dominated problem.

This is one reason the Step-4 integration explicitly provides both state-operator and
transpose-state-operator actions even though its current stiffness matrix is
symmetric.

## 9. Transposing the semilinear weak form

The semilinear example from Section 4 also illustrates that a VJP can be derived
without first writing a Jacobian matrix.

We had

```math
\langle E'(y,u)[\delta y,\delta u],p\rangle
=
\int_{\Omega}\nabla\delta y\cdot\nabla p
+
\int_{\Omega}3c y^{2}\delta y p
-
\int_{\Omega}\delta u p.
```

Regroup the expression by the perturbations:

```math
\begin{aligned}
\langle E'(y,u)[\delta y,\delta u],p\rangle
&=
\left[
\int_{\Omega}\nabla\delta y\cdot\nabla p
+
\int_{\Omega}3c y^{2}\delta y p
\right]\\
&\quad
-
\int_{\Omega}\delta u p.
\end{aligned}
```

The bracketed term is a linear functional of $\delta y$; the last term is a linear
functional of $\delta u$.

Together they define

$$
E'(y,u)^{\ast}p
\in
Y^{\ast}\times U^{\ast}.
$$

For this particular reaction term the state linearization is symmetric in
$\delta y$ and $p$, so the weak transpose has the same-looking integrals.

That is a property of the operator, not a general rule.

The useful lesson is that VJP derivation can remain at the weak-form level:

```text
write the tested linearized residual
        ↓
hold the test seed fixed
        ↓
collect terms by variable perturbation
        ↓
read off the covector components
```

This viewpoint scales better to mixed, boundary, and nonlinear operators than
memorizing matrix transposes.

## 10. Why reverse actions are useful in optimization

The JVP and VJP contain the same linearized information viewed in opposite
directions, but their computational usefulness depends on what one wants to compute.

Suppose

$$
J_{E}
\in
\mathbb R^{m\times n}.
$$

A JVP accepts one direction in the $n$-dimensional variable space and produces its
effect on all $m$ residual coordinates.

A VJP accepts one seed in the $m$-dimensional test space and produces its effect on
all $n$ variable coordinates.

This is the same forward-versus-reverse distinction that appears in automatic
differentiation.

For PDE-constrained optimization, the adjoint method is a reverse operation of
exactly this kind.

After choosing one adjoint seed $p$, the VJP

$$
E'(x)^{\ast}p
$$

returns covector contributions for every variable block at once.

For a state/control product, one call conceptually gives

$$
\left(
D_{y}E_{h}(x)^{\ast}p,
D_{u}E_{h}(x)^{\ast}p
\right).
$$

The reduced method can then use the control component without solving one state
sensitivity problem for every control basis direction.

So the VJP is not an API convenience added after the adjoint theory. It is the
numerical operation that makes the adjoint theory executable.

## 11. Objective derivatives fit the same dual picture

The objective

$$
J_{h}:X_{h}\to\mathbb R
$$

has derivative

$$
J_{h}'(x)\in X_{h}^{\ast}.
$$

For a perturbation $\delta x$,

$$
J_{h}'(x)[\delta x]
=
\left\langle
J_{h}'(x),
\delta x
\right\rangle.
$$

In coordinates,

$$
=
\mathbf r_{J}^{\mathsf T}\delta\mathbf x.
$$

Because the codomain of the objective is already scalar, there is usually no reason
for the executable interface to expose a separate "objective JVP" method: the
returned derivative covector can be paired with any tangent direction.

For the quadratic objective from the first chapter,

```math
J_{h}(y,u)
:=
\frac{1}{2}
\left(
y^{\mathsf T}M_{y}y
-
2q^{\mathsf T}y
+
c
\right)
+
\frac{\beta}{2}
u^{\mathsf T}N_{u}u,
```

the derivative covector is

```math
J_{h}'(y,u)
=
\begin{bmatrix}
M_{y}y-q\\
\beta N_{u}u
\end{bmatrix}.
```

Pairing it with $(\delta y,\delta u)$ gives the directional derivative.

Again, this is a covector. No Riesz inverse is applied merely because the API calls
the operation `objective_derivative()`.

## 12. Objective derivative checks mirror residual derivative checks

A first-order Taylor expansion gives

```math
J_{h}(x+\varepsilon\delta x)
=
J_{h}(x)
+
\varepsilon J_{h}'(x)[\delta x]
+
\mathcal O(\varepsilon^{2}).
```

So the scalar remainder

$$
R_{J}(\varepsilon)
:=
J_{h}(x+\varepsilon\delta x)
-
J_{h}(x)
-
\varepsilon
\left\langle
J_{h}'(x),
\delta x
\right\rangle
$$

should decay quadratically for a smooth objective.

This test checks a different relationship from the JVP/VJP pairing:

```text
Taylor / finite-difference test
    residual or objective
    agrees with its derivative

transpose-pairing test
    JVP
    agrees with the VJP
```

Both are useful because a JVP and VJP can be mutually transposed yet both be the
transpose of the **wrong** derivative. The pairing test alone cannot establish that
they differentiate the intended residual.

## 13. The executable-model interface now has a direct mathematical reading

At this point the five main methods of `ExecutableModelT` are no longer an arbitrary
API list.

The interface represents the maps

```math
\begin{aligned}
E_{h}(x)
&\in
Z_{h}^{\ast},
\\
E_{h}'(x)[\delta x]
&\in
Z_{h}^{\ast},
\\
E_{h}'(x)^{\ast}p
&\in
X_{h}^{\ast},
\\
J_{h}(x)
&\in
\mathbb R,
\\
J_{h}'(x)
&\in
X_{h}^{\ast}.
\end{aligned}
```

The corresponding methods are:

```cpp
residual(x)
residual_jvp(x, dx)
residual_vjp(x, p)
objective(x)
objective_derivative(x)
```

The variable and test layouts provide the runtime coordinate-space identities
developed in the previous chapter.

The type flow can be summarized as:

| Operation | Input | Output | Mathematical meaning |
| --- | --- | --- | --- |
| `residual` | primal in $X_{h}$ | covector in $Z_{h}^{\ast}$ | $E_{h}(x)$ |
| `residual_jvp` | point and tangent in $X_{h}$ | covector in $Z_{h}^{\ast}$ | $E_{h}'(x)[\delta x]$ |
| `residual_vjp` | point in $X_{h}$, seed in $Z_{h}$ | covector in $X_{h}^{\ast}$ | $E_{h}'(x)^{\ast}p$ |
| `objective` | primal in $X_{h}$ | scalar | $J_{h}(x)$ |
| `objective_derivative` | primal in $X_{h}$ | covector in $X_{h}^{\ast}$ | $J_{h}'(x)$ |

The table also explains one detail that can look odd on first inspection:
`residual_vjp()` takes a **primal** test seed.

The seed $p$ is an element of the test space $Z_{h}$. The result is the dual object
in $X_{h}^{\ast}$.

### 13.1 Why the interface exposes actions rather than a Jacobian object

An alternative design could have returned a Jacobian matrix from something like

```text
jacobian(x)
```

and then expected callers to perform matrix-vector and transpose-matrix-vector
products.

The action-oriented interface avoids committing the generic boundary to:

- an assembled sparse-matrix type;
- a particular storage scheme;
- the assumption that JVP and VJP are derived from one materialized matrix;
- a deal.II-specific operator object.

An assembled realization can still implement the callbacks with `vmult` and `Tvmult`.
A matrix-free or application-owned realization can provide equivalent actions by
other means.

The mathematical invariant is the action and its transpose pairing, not the storage
format.

## 14. `CallbackExecutableModelT` adapts native operations to those roles

`CallbackExecutableModelT` is a small concrete implementation of the interface.

It stores five `std::function` objects:

```text
ResidualAction
LinearizedAction
TransposeAction
ObjectiveAction
ObjectiveDerivativeAction
```

The class adds two kinds of structure around the native callbacks:

1. it fixes the variable and test layouts;
2. it checks that inputs and returned values have the expected layouts.

The callbacks themselves remain free to call native application methods, compiler
objects, assembled operators, or other numerical services.

This is easiest to see for the derivative actions.

The JVP adapter receives:

```text
point             Primal(variable_layout)
variable tangent  Primal(variable_layout)
```

and requires the callback to return

```text
Covector(test_layout)
```

The VJP adapter receives:

```text
point      Primal(variable_layout)
test seed  Primal(test_layout)
```

and requires

```text
Covector(variable_layout)
```

So the callback wrapper does not derive the mathematics. It gives already-existing
native numerical operations a checked generic signature.

### 14.1 The point argument can be redundant for a linear model

In the Step-4 Problem B implementation, the residual is affine:

$$
E(z,u)=Kz-b_{F}-Bu.
$$

Its derivative is therefore independent of the current point.

The native method can simply expose

```cpp
residual_jvp(state_tangent, control_tangent)
```

without receiving $(z,u)$.

When that method is adapted to `CallbackExecutableModelT`, the generic callback still
has the signature

```cpp
(point, variable_tangent)
```

and ignores `point`.

That is not wasted information in the general interface. A nonlinear realization
would need the point to assemble or apply $E'(x)$.

The same observation applies to the VJP.

## 15. Step-4 Problem B as an assembled operator example

The native Problem B implementation makes the matrix interpretation unusually easy
to see.

Its residual is conceptually

$$
E(z,u)
:=
Kz-b_{F}-Bu,
$$

where:

- $z$ contains independent state coordinates;
- $K$ is the state operator in those coordinates;
- $B$ couples the full control field into the independent state equations;
- $b_{F}$ is the fixed right-hand side after boundary treatment.

The actual code does not spell this as one matrix expression. It composes native
operations:

```text
apply_state_operator(state)
-
free_system_rhs
-
coupling_apply(control)
```

That is the residual.

### 15.1 JVP

Because the residual is affine, its JVP removes the fixed load and applies the same
linear operators to the tangent:

```text
apply_state_operator(state_tangent)
-
coupling_apply(control_tangent)
```

Mathematically,

$$
E'(z,u)[\delta z,\delta u]
=
K\delta z-B\delta u.
$$

The method does not need the current $(z,u)$.

### 15.2 VJP

Given a test seed $p$, the transpose action is

```math
E'(z,u)^{\ast}p
=
\begin{bmatrix}
K^{\mathsf T}p\\
-B^{\mathsf T}p
\end{bmatrix}.
```

The native implementation therefore calls:

```text
apply_transpose_state_operator(test_seed)
coupling_transpose_apply(test_seed)
```

and multiplies the control result by $-1$.

The state transpose path uses deal.II's `Tvmult`, while the control contribution uses
the explicitly transposed coupling action.

This is a very literal realization of the pairing derivation from Sections 6–8.

### 15.3 The binding adds the space roles

The native Problem B methods accept and return ordinary deal.II vectors.

`ProblemBBinding` then wraps:

```text
native residual / JVP
    → Covector(test_layout)

native VJP state + control pieces
    → Covector(variable_layout)
```

The executable-model interface therefore sees exactly the mathematical roles derived
above even though the application itself remains written in native vector terms.

## 16. Objective derivatives in Step-4 also cross a coordinate map

Problem B's objective uses the physical state, not the 225-dimensional independent
state vector directly.

The native code reconstructs

$$
y_{\mathrm{phys}}
=
Pz+\ell
$$

and forms the physical mass action

$$
M y_{\mathrm{phys}}.
$$

But the objective derivative returned to nmopt must act on an independent state
perturbation $\delta z$.

As derived in the previous chapter, the corresponding state covector is pulled back
through the coordinate map.

For the simple fixed-index coordinate class, that pullback is implemented by
restriction to the free entries. In the general matrix notation it is

$$
P^{\mathsf T}M y_{\mathrm{phys}}.
$$

The control derivative remains the full-space mass action

$$
M u.
$$

So the objective derivative has the product-space form

```math
J'(z,u)
=
\begin{bmatrix}
P^{\mathsf T}M(Pz+\ell)\\
M u
\end{bmatrix}
```

for the Step-4 quadratic objective.

This is a useful example because it combines two different transformations:

```text
primal state
    independent → physical
    P z + ell

state derivative
    physical dual → independent dual
    P^T r
```

The executable interface sees only the final independent-coordinate covector.

## 17. The pairing test is the defining transpose test

For any point $x$, tangent $\delta x$, and test seed $p$, a correct JVP/VJP pair must
satisfy

```math
\left\langle
E'(x)[\delta x],
p
\right\rangle
=
\left\langle
E'(x)^{\ast}p,
\delta x
\right\rangle.
```

After coordinate representation this becomes

$$
\mathrm{pair}(
\mathrm{JVP}(x,\delta x),
p)
=
\mathrm{pair}(
\mathrm{VJP}(x,p),
\delta x).
$$

The backend-neutral executable-model contract test checks exactly this identity.

Its synthetic residual is

```math
E(y_{0},y_{1},u)
=
\begin{bmatrix}
y_{0}-u\\
y_{1}-2u
\end{bmatrix}.
```

Therefore

```math
E'
=
\begin{bmatrix}
1&0&-1\\
0&1&-2
\end{bmatrix},
```

and

```math
E'^{\mathsf T}p
=
\begin{bmatrix}
p_{0}\\
p_{1}\\
-p_{0}-2p_{1}
\end{bmatrix}.
```

The test supplies arbitrary tangent and seed vectors and compares the two pairings.

That test is small, but it captures the essential contract independently of
deal.II.

### 17.1 Why not just compare matrix entries?

For a simple assembled example, comparing a stored transpose matrix to the original
Jacobian may be possible.

The pairing identity is more general:

- it works for matrix-free actions;
- it works when JVP and VJP are assembled through different code paths;
- it works across block/product spaces;
- it tests the actual callable operations seen by the formulation.

It is therefore the natural executable definition of transpose consistency.

## 18. Pairing consistency and derivative consistency are complementary

The two main derivative checks answer different questions.

### Residual/JVP consistency

Does the JVP differentiate the residual?

Check

```math
E(x+\varepsilon\delta x)
-
E(x)
\approx
\varepsilon E'(x)[\delta x].
```

### JVP/VJP consistency

Is the VJP the transpose of that JVP under the declared pairings?

Check

```math
\langle E'(x)[\delta x],p\rangle
=
\langle E'(x)^{\ast}p,\delta x\rangle.
```

If only the second test is present, one could implement a wrong JVP and the matching
transpose of that wrong JVP; the pairing test would still pass.

If only the first test is present, the VJP could be wrong while the JVP is correct.

A serious derivative implementation therefore benefits from both kinds of evidence.

The later verification chapter will discuss tolerances, Taylor rates, randomized
directions, and solver-induced error in more detail.

## 19. Three uses of the word "adjoint"

PDE-constrained optimization uses the word *adjoint* in several closely related
ways. Keeping them separate makes the source much easier to read.

### The adjoint or transpose action

This is the map

$$
E'(x)^{\ast}:
Z_{h}
\longrightarrow
X_{h}^{\ast}.
$$

Given a seed $p$, it returns the variable-space covector

$$
E'(x)^{\ast}p.
$$

In nmopt this is the role of the VJP action.

### The adjoint variable

In a reduced state/control problem, the **adjoint variable** $p$ is a particular
element of the test space chosen to satisfy an adjoint equation.

If

$$
x=(y,u),
$$

the state block of the residual derivative is

$$
D_{y}E_{h}(y,u):
Y_{h}
\longrightarrow
Z_{h}^{\ast}.
$$

The adjoint equation has the form

$$
D_{y}E_{h}(y,u)^{\ast}p
=
D_{y}J_{h}(y,u)
\in
Y_{h}^{\ast}
$$

under the project's Lagrangian sign convention.

The adjoint variable is therefore a *primal element of the test space* whose image
through the transpose state derivative matches the state objective derivative.

### The adjoint solve

The **adjoint solve** is the numerical process that finds this $p$.

For a linear matrix problem it may look like

$$
A^{\mathsf T}p=d_{y}.
$$

For a nonlinear PDE it uses the transpose of the state Jacobian evaluated at the
current state/control point.

These three notions are connected:

```text
transpose action
    defines D_y E_h(x)* p for a supplied p

adjoint equation
    asks for p such that D_y E_h(x)* p = D_y J_h

adjoint solve
    numerically computes that p
```

The distinction explains why a VJP callback and an adjoint-solve callback are separate
capabilities in the reduced formulation.

## 20. Applying an operator is not the same as solving with it

This separation is easy to miss in a matrix example.

If an application provides

$$
p
\longmapsto
A^{\mathsf T}p,
$$

it can apply the adjoint operator.

That does not automatically mean it knows how to solve

$$
A^{\mathsf T}p=b.
$$

A solve may require:

- an assembled matrix;
- a preconditioner;
- nullspace information;
- boundary-condition treatment;
- an iterative tolerance;
- a nonlinear iteration in more general settings.

The executable model therefore describes **model evaluations and derivative
actions**. State and adjoint solves are supplied separately by the formulation's
numerical services.

For the reduced path this separation looks conceptually like:

```text
ExecutableModelT
    E(x)
    E'(x) δx
    E'(x)* p
    J(x)
    J'(x)

StateAdjointSolversT
    control → state satisfying E(y,u)=0
    state RHS → p satisfying adjoint equation
```

The first group describes operators.

The second group provides inverse/solve capabilities needed by one particular
formulation.

This becomes useful for external applications: a mature PDE code may already own a
carefully tuned state solver and adjoint solver, while exposing derivative actions
through different native methods.

## 21. Re-derive the reduced derivative using the VJP language

We can now revisit the derivation from the first chapter in the language of operator
actions.

Let

$$
x=(y,u)
$$

be a feasible state/control point.

The objective derivative is the product covector

$$
J'(x)
=
\left(
D_{y}J,
D_{u}J
\right)
\in
Y_{h}^{\ast}\times U_{h}^{\ast}.
$$

Solve the adjoint equation

$$
D_{y}E_{h}(x)^{\ast}p
=
D_{y}J_{h}.
$$

Now evaluate the **full** residual VJP:

$$
E'(x)^{\ast}p
=
\left(
D_{y}E_{h}(x)^{\ast}p,
D_{u}E_{h}(x)^{\ast}p
\right).
$$

The first component is, by construction, $D_{y}J$.

The second component gives the control contribution from the state equation.

With the project's convention

$$
\mathcal L(x,p)
=
J(x)-\langle p,E(x)\rangle,
$$

the control derivative of the reduced objective is

$$
j'(u)
=
D_{u}J_{h}
-
D_{u}E_{h}(x)^{\ast}p.
$$

For the running residual

$$
E(y,u)=Ay-Bu-f,
$$

we have

$$
D_{u}E_{h}(x)^{\ast}p
=
-B^{\mathsf T}p,
$$

and therefore

$$
j'(u)
=
D_{u}J+B^{\mathsf T}p.
$$

This formulation shows why the generic VJP returns **all variable blocks** rather
than only a control derivative.

The executable model describes the whole derivative

$$
E'(x)^{\ast}p\in X_{h}^{\ast}.
$$

The reduced formulation knows which block is state and which is control and extracts
the pieces relevant to the adjoint and reduced derivative.

## 22. Why the current reduced path can use only part of the executable port

The current `ExecutableModelT` describes a complete first-order model:

```text
residual
JVP
VJP
objective
objective derivative
```

The current first-order reduced state–adjoint path does not need all five actions
during every optimization run.

Once separate state and adjoint solves are supplied, reduced value evaluation needs
the objective, while derivative augmentation uses the objective derivative and VJP.

In the current Step-4 reduced path, the generic `residual()` and `residual_jvp()`
operations are present in the model contract but are not needed by the optimizer's
main first-order execution path.

This does not change their mathematical meaning. It reflects that the executable
model is a broader first-order port than the minimal set of operations consumed by
that one formulation.

Other formulations, diagnostics, or future algorithms can make different use of the
same actions.

The important distinction for a reader is:

```text
what a model is capable of exposing
        is not identical to
what one algorithm happens to call
```

## 23. Why not compute the VJP from the JVP automatically?

If a full Jacobian matrix were always available, one could implement

```text
JVP = J * dx
VJP = J^T * p
```

from the same stored object.

An action-oriented abstraction cannot generally recover the transpose action from a
black-box forward action.

Given only a function that computes

$$
\delta x
\mapsto
E'(x)[\delta x],
$$

there is no cheap generic procedure for constructing

$$
p
\mapsto
E'(x)^{\ast}p
$$

in a high-dimensional problem.

One could probe every basis vector to reconstruct the Jacobian, but that destroys the
point of a matrix-free/action-only interface.

The producer must therefore provide the transpose capability explicitly, or use a
technology such as reverse-mode automatic differentiation that can produce it.

This is why `CallbackExecutableModelT` asks for both `residual_jvp` and
`residual_vjp` callbacks.

The pairing test then checks that the two independently supplied actions agree.

## 24. A compact map of the first-order model

At this point the relationships among the objects can be summarized as:

| Object/action | Map | Coordinate picture | nmopt operation |
| --- | --- | --- | --- |
| Residual | $E:X_{h}\to Z_{h}^{\ast}$ | $\mathbf E(\mathbf x)$ | `residual(x)` |
| Residual derivative | $E'(x):X_{h}\to Z_{h}^{\ast}$ | $J_{E}(\mathbf x)\delta\mathbf x$ | `residual_jvp(x, dx)` |
| Residual transpose | $E'(x)^{\ast}:Z_{h}\to X_{h}^{\ast}$ | $J_{E}(\mathbf x)^{\mathsf T}\mathbf p$ | `residual_vjp(x, p)` |
| Objective | $J:X_{h}\to\mathbb R$ | $J(\mathbf x)$ | `objective(x)` |
| Objective derivative | $J'(x)\in X_{h}^{\ast}$ | $\mathbf r_{J}$ | `objective_derivative(x)` |
| State solve | solve $E(y,u)=0$ for $y$ | inverse/nonlinear solve | separate solve service |
| Adjoint solve | solve $D_{y}E_{h}(x)^{\ast}p=D_{y}J_{h}$ | transpose-system solve | separate solve service |

The table also marks the boundary of this chapter.

We have explained how first-order information is represented and applied.

We have not yet selected a Riesz map to turn an objective/reduced derivative into a
gradient direction. That belongs to the next chapter.

## 25. Following one action through the source

A useful way to navigate the relevant code is to choose one mathematical action and
follow it from the native implementation to the generic contract.

Take the Step-4 VJP.

### Native numerical implementation

Start in:

- [`apps/external-dealii/step-4/integration/problem_b.hpp`](../../../apps/external-dealii/step-4/integration/problem_b.hpp)

Find `ProblemB::residual_vjp()`.

It:

1. applies the transpose state operator with `Tvmult`;
2. applies the transpose control coupling;
3. applies the residual's minus sign to the control component;
4. returns the two native vectors as a small application-owned record.

At this level the code knows native matrices and coordinate maps.

### Generic binding

Then read:

- [`apps/external-dealii/step-4/minimal/problem_b_binding.hpp`](../../../apps/external-dealii/step-4/minimal/problem_b_binding.hpp)

The VJP callback receives a `Primal` test seed, calls the native method, and wraps the
two returned vectors as a `Covector` with the variable layout.

At this level the code supplies the generic mathematical roles.

### Abstract model

Finally read:

- [`include/nmopt/contract/executable_model.hpp`](../../../include/nmopt/contract/executable_model.hpp)

The abstract signature contains no reference to Step-4, sparse matrices, deal.II
`Tvmult`, or the state/control details of Problem B.

It only says:

```text
point in variable space
+
seed in test space
→
covector in variable space
```

The three files are therefore three views of one operation.

## 26. The backend-neutral contract test is worth reading after the interface

After `executable_model.hpp`, read:

- [`include/nmopt/contract/callback_executable_model.hpp`](../../../include/nmopt/contract/callback_executable_model.hpp)
- [`tests/contract/executable_model_contract.cc`](../../../tests/contract/executable_model_contract.cc)

The contract test uses a deliberately tiny model:

```math
E(y_{0},y_{1},u)
=
\begin{bmatrix}
y_{0}-u\\
y_{1}-2u
\end{bmatrix}.
```

Because every derivative can be written by inspection, it is an effective place to
see:

- the expected block layouts;
- residual/JVP/VJP return roles;
- objective derivative blocks;
- layout-rejection behavior;
- the JVP/VJP pairing identity;
- consumption by the reduced DTO formulation.

This test is much easier to understand after the mathematics than before it.

## 27. What changes for compiler-created models?

The semantic/compiler path produces the same kinds of solver-facing actions, but the
native implementation is assembled from declared residual terms and realization
services rather than hand-written in one application class.

For example, a diffusion term, volume-control term, or transport term contributes to:

- the residual value;
- its linearization with respect to relevant variable blocks;
- the corresponding transpose actions.

The compiler's job is not merely to construct a sparse matrix named "Jacobian". It
must realize the set of actions required by the selected formulation and package
them with compatible spaces, data, and solve services.

The later **Semantic problem model**, **Compilation and lowering**, and
**Finite-element realization** chapters will trace that construction in detail.

For now, the important commonality is:

```text
external application
        │
        └─ native residual / JVP / VJP actions
                         │
                         ▼
                 ExecutableModelT

semantic/compiler path
        │
        └─ compiled residual / JVP / VJP actions
                         │
                         ▼
                 ExecutableModelT
```

The producer paths differ. The first-order mathematical vocabulary at the generic
model boundary is the same.

## 28. What to carry into the next chapter

We can now distinguish four operations that are often collapsed in informal
numerical code:

```text
evaluate the residual
    E(x)

differentiate the residual forward
    E'(x) δx

apply the transpose derivative
    E'(x)* p

solve an equation involving one of those operators
    E(y,u)=0
    or
    D_y E_h(x)* p = rhs
```

Likewise, we can distinguish

```text
objective derivative
    J'(x) ∈ X*

from

gradient
    an element of X obtained after choosing a Riesz map
```

The first three concept chapters have therefore brought us to a natural boundary.

We know what the finite-dimensional fields and coordinates represent. We know how
first-order model information acts between their spaces. We know how adjoint
variables arise from transpose actions and solves.

What we have not yet explained is how a covector becomes a primal gradient or search
direction, why different control metrics produce different gradients, or why
constraint projection can depend on the chosen metric realization.

Those are the subjects of **Metrics, gradients, and constraints**.

## Read later

Useful existing documents for this chapter are:

- [Theoretical formalism](../../design/theoretical-formalism.md), for the project's
  formal derivative, transpose, and Lagrangian conventions.
- [Reduced optimization](../../overview/reduced-optimization.md), for the high-level
  state–adjoint execution path.
- [External applications](../../overview/external-applications.md), for the role of
  callback model actions in an independently owned PDE application.
- [Step-4 external integration overview](../../../apps/external-dealii/step-4/external-integration-overview.md),
  for the complete Problem B case study.
- [External deal.II solver integration](../../reference/external-dealii-solver-integration.md),
  for exact public integration signatures.

The forthcoming **Verification and evidence** chapter will return to finite
differences, Taylor tests, transpose-pairing tests, and solve evidence as a coherent
testing strategy rather than treating them only as derivative definitions.
