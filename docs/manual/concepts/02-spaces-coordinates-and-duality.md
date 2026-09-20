# Spaces, coordinates, and duality

## Why a vector is not yet a numerical object

At the end of
[Anatomy of a discrete PDE-constrained problem](01-discrete-problem-anatomy.md), our
continuous control problem had become a collection of finite-dimensional vectors and
operators:

```text
state coefficients
control coefficients
test directions
residual coefficients
objective derivatives
adjoint coefficients
```

At first sight this seems to simplify everything. Once all fields have been expanded
in finite-element bases, do we not simply have vectors in Euclidean spaces?

The difficulty is that the same storage type can represent very different objects.

A `dealii::Vector<double>` with 225 entries might represent an independent state. A
different 225-entry vector might represent a test function. Another might encode a
residual functional acting on that test space. Even when two vectors have identical
length and both consist of ordinary doubles, adding them or passing one where the
other is expected can be mathematically meaningless.

There is a second complication. The coefficient vector manipulated by an optimizer
need not even contain all coefficients of the physical finite-element field. Essential
boundary conditions and other affine constraints can reduce the independent
coordinates, while output and objective evaluation may still use the full physical
field.

This chapter develops the distinctions that make those situations precise. Only
after the mathematical picture is clear will we introduce nmopt's `BlockLayout`,
`PrimalBlockT`, and `CovectorBlockT`.

The central progression is:

```text
function space
    ↓ choose a finite-dimensional subspace
discrete space
    ↓ choose basis / independent coordinates
coordinate vector
    ↓ attach mathematical role
primal or dual coordinate representation
    ↓ retain selected runtime identity
nmopt block value
```

The framework types appear at the bottom of this chain. They are representations of
choices made above them, not substitutes for the underlying mathematics.

## 1. Four things that are easy to call "the space"

It is useful to distinguish several levels that numerical PDE code often discusses
with the same shorthand.

Consider again the state space from the previous chapter.

At the continuous level we used

$$
V:=H_{0}^{1}(\Omega).
$$

After choosing a mesh and finite-element family we obtained a finite-dimensional
subspace

$$
V_{h}\subset V.
$$

After choosing a basis
$`\{\varphi_{1},\ldots,\varphi_{n}\}`$, every $`v_{h}\in V_{h}`$ could be represented by
coordinates

$$
v_{h}
=
\sum_{i=1}^{n}v_{i}\varphi_{i}
\quad\longleftrightarrow\quad
\mathbf v
:=
(v_{1},\ldots,v_{n})^{\mathsf T}\in\mathbb R^{n}.
$$

Finally, a C++ program stores those coordinates in some native vector type.

These are related, but they are not the same object:

| Level | Example | What it tells us |
| --- | --- | --- |
| Continuous space | $`H_{0}^{1}(\Omega)`$ | regularity, boundary meaning, topology |
| Discrete FE space | $`V_{h}`$ on a chosen mesh and element family | basis functions and admissible discrete fields |
| Coordinate space | $\mathbb R^{n}$ with a chosen basis interpretation | finite coefficients used by algebra |
| Native storage | `dealii::Vector<double>` | how those coefficients are stored and manipulated |

The last two rows are where confusion most often enters.

Once a field is represented by $\mathbf v\in\mathbb R^{n}$, ordinary linear algebra
does not remember that the coefficients came from an $H^{1}$ state space rather than
an $L^{2}$ control space, a boundary trace space, or a multiplier space. The
interpretation survives only because the surrounding program remembers it.

### 1.1 Equal dimension does not imply equal space

Suppose two finite-element spaces happen to have the same number of coefficients:

$$
\dim Y_{h}=\dim U_{h}=289.
$$

Then both state and control coordinates are elements of $\mathbb R^{289}$ as plain
arrays.

That does not identify $`Y_{h}`$ with $`U_{h}`$.

The state basis might consist of nodal basis functions satisfying one set of
constraints, while the control basis could describe a different field, different
boundary behavior, or simply a different semantic role. Even if the basis functions
happen to be numerically identical, the state and control remain different factors
of the optimization problem.

The Step-4 external integration contains both kinds of example:

- Problem A uses 289 state and 289 control coefficients;
- Problem B uses 225 independent state coefficients and 289 control coefficients.

Problem B makes a dimension mismatch visible immediately. Problem A is more subtle:
dimension alone cannot distinguish the roles at all.

This is one reason nmopt's runtime checks use more than vector length.

### 1.2 A coordinate vector depends on the chosen basis

The statement

$$
\mathbf v\in\mathbb R^{n}
$$

has meaning only relative to a basis or coordinate map.

If the basis changes, the coefficient vector changes while the represented physical
field may remain exactly the same.

This is familiar in elementary linear algebra, but it matters especially in finite
elements because several coordinate systems can coexist around one field:

- full nodal coefficients;
- free coefficients after essential constraints;
- coefficients on a boundary trace;
- coefficients transferred from another mesh;
- coefficients in a reduced or parameterized basis.

The next section gives a concrete example that appears directly in the repository.

## 2. Physical state fields and independent coordinates

Essential boundary conditions create a useful distinction between the physical
finite-element field and the coordinates that are actually free to vary.

Suppose a finite-element discretization has $N$ physical degrees of freedom. If some
of those coefficients are fixed by Dirichlet data, the admissible state fields do not
fill all of $\mathbb R^{N}$.

Let

$$
\mathbf z\in\mathbb R^{n},
\qquad
n<N,
$$

contain only the independent state coordinates.

The full physical coefficient vector can be reconstructed in affine form as

$$
\mathbf y_{\mathrm{phys}}
:=
P\mathbf z+\boldsymbol\ell.
$$

Here:

- $P\in\mathbb R^{N\times n}$ embeds homogeneous independent coordinates into the
  full finite-element vector;
- $\boldsymbol\ell\in\mathbb R^{N}$ is a fixed lifting that carries the prescribed
  inhomogeneous values.

This formula is more informative than saying merely that boundary values are
"applied".

It separates two different pieces of the physical state:

```text
P z
    the part that changes when the optimizer/state solver changes z

ell
    the fixed affine contribution imposed by the physical constraints
```

### 2.1 The simplest case: fixed boundary nodes

If the only constraints are fixed Dirichlet coefficients, $P$ can be understood as
an insertion matrix.

For a toy vector with five physical coefficients, suppose entries 1 and 4 are fixed.
The free coordinates are

$$
\mathbf z
=
\begin{bmatrix}
z_{0}\\
z_{1}\\
z_{2}
\end{bmatrix}.
$$

One possible reconstruction is

```math
P
:=
\begin{bmatrix}
1&0&0\\
0&0&0\\
0&1&0\\
0&0&1\\
0&0&0
\end{bmatrix},
\qquad
\boldsymbol\ell
:=
\begin{bmatrix}
0\\
g_{1}\\
0\\
0\\
g_{4}
\end{bmatrix}.
```

Then

```math
P\mathbf z+\boldsymbol\ell
=
\begin{bmatrix}
z_{0}\\
g_{1}\\
z_{1}\\
z_{2}\\
g_{4}
\end{bmatrix}.
```

The matrix contains no PDE physics. It is a coordinate map between the independent
state representation and the physical finite-element coefficient vector.

### 2.2 Step-4 Problem B makes this distinction concrete

In the 2D Step-4 case study, the physical deal.II field has 289 degrees of freedom.
After the prescribed boundary values are removed from the independent state
coordinates, Problem B has 225 free state coefficients.

The application therefore distinguishes

$$
\mathbf z\in\mathbb R^{225}
$$

from

$$
\mathbf y_{\mathrm{phys}}\in\mathbb R^{289}.
$$

Its local `ProblemBCoordinates` class stores the free indices and the fixed lifting.
`reconstruct(z)` inserts the free coefficients into a copy of the lifting, while
`embed_free(z)` performs the homogeneous embedding without adding the fixed boundary
values.

That second operation is needed because perturbations do not carry the affine
lifting.

If

$$
\mathbf y_{\mathrm{phys}}
=
P\mathbf z+\boldsymbol\ell,
$$

then a perturbation satisfies

$$
\delta\mathbf y_{\mathrm{phys}}
=
P\delta\mathbf z.
$$

The fixed vector $\boldsymbol\ell$ disappears under differentiation.

This small observation will matter repeatedly when we discuss objective derivatives,
residual linearizations, and adjoints.

## 3. The embedding can encode more than fixed boundary values

The insertion-matrix picture is useful, but it is not the most general finite-element
case.

deal.II can impose homogeneous affine constraints such as hanging-node relations.
Then a constrained physical coefficient may depend linearly on several independent
coefficients.

In that situation, a row of $P$ need not contain either one `1` or all zeros. It can
contain interpolation weights.

The reusable `nmopt::dealii_backend::IndependentStateCoordinates` class is built for
this more general situation.

It receives:

- the physical coefficient dimension;
- a set of homogeneous `dealii::AffineConstraints<double>`;
- a fixed lifting vector.

It finds the unconstrained degrees of freedom and constructs the reconstruction
matrix $P$ column by column. Conceptually, for every independent coordinate it:

1. creates the corresponding physical standard basis vector;
2. distributes the homogeneous affine constraints;
3. records the resulting constrained physical vector as one column of $P$.

The resulting operations are exactly the ones suggested by the mathematics:

```text
embed(z)        → P z
reconstruct(z)  → P z + ell
pullback(r)     → P^T r
```

The last operation is our first explicit encounter with duality.

## 4. Why dual objects transform with the transpose

The transpose pullback becomes completely transparent if we distinguish the scalar
quantity defined on **physical coefficients** from the scalar quantity obtained after
we re-express the state through independent coordinates.

Let

$$
R(\mathbf z)
:=
P\mathbf z+\boldsymbol\ell
$$

be the affine reconstruction map from independent coordinates to the physical
finite-element vector. Suppose

$$
\Phi_{\mathrm{phys}}:
\mathbb R^{N}
\longrightarrow
\mathbb R
$$

is a scalar quantity evaluated from the physical state coefficients. The quantity
seen as a function of the independent coordinates is the composition

$$
\widehat\Phi(\mathbf z)
:=
\Phi_{\mathrm{phys}}(R(\mathbf z))
=
\Phi_{\mathrm{phys}}(P\mathbf z+\boldsymbol\ell).
$$

This distinction matters because a perturbation
$\delta\mathbf z\in\mathbb R^{n}$ is **not** an admissible argument of
$`D\Phi_{\mathrm{phys}}`$ by itself. The derivative of $`\Phi_{\mathrm{phys}}`$ acts on
physical perturbations in $\mathbb R^{N}$.

Assume that, at

$$
\mathbf y_{\mathrm{phys}}
:=
R(\mathbf z),
$$

the physical-coordinate derivative is represented by a covector
$`\mathbf r_{\mathrm{phys}}`$:

$$
D\Phi_{\mathrm{phys}}(\mathbf y_{\mathrm{phys}})
[\delta\mathbf y_{\mathrm{phys}}]
=
\mathbf r_{\mathrm{phys}}^{\mathsf T}
\delta\mathbf y_{\mathrm{phys}}.
$$

Now perturb the **independent** coordinates by $\delta\mathbf z$. The reconstruction
map sends that perturbation to

$$
DR(\mathbf z)[\delta\mathbf z]
=
P\delta\mathbf z.
$$

There is no contribution from $\boldsymbol\ell$ because the lifting is fixed.

Applying the chain rule to the composed function $`\widehat\Phi=\Phi_{\mathrm{phys}}\circ R`$
gives

```math
\begin{aligned}
D\widehat\Phi(\mathbf z)[\delta\mathbf z]
&=
D\Phi_{\mathrm{phys}}(R(\mathbf z))
\left[
DR(\mathbf z)[\delta\mathbf z]
\right]\\
&=
D\Phi_{\mathrm{phys}}(\mathbf y_{\mathrm{phys}})
[P\delta\mathbf z]\\
&=
\mathbf r_{\mathrm{phys}}^{\mathsf T}
P\delta\mathbf z\\
&=
\left(P^{\mathsf T}\mathbf r_{\mathrm{phys}}\right)^{\mathsf T}
\delta\mathbf z.
\end{aligned}
```

The $P$ has therefore not disappeared from the derivative argument. It first pushes
the independent perturbation into physical coordinates. Only after that step do we
move it algebraically to the other side of the pairing, where it appears as
$P^{\mathsf T}$ acting on the covector.

The derivative of the **pulled-back scalar function** $\widehat\Phi$ is consequently
represented in independent coordinates by

$$
\mathbf r_{\mathrm{ind}}
:=
P^{\mathsf T}\mathbf r_{\mathrm{phys}}.
$$

It is useful to picture the two transformations together:

```text
primal perturbation
    δz  ──────────►  P δz
    independent      physical

covector
    P^T r_phys  ◄──  r_phys
    independent      physical
```

Primal perturbations are pushed forward through the coordinate map. Covectors are
pulled back through its transpose.

This is not a finite-element peculiarity. It is the standard chain-rule transformation
of a linear functional under a change of coordinates. The repository's reusable
coordinate class reflects the same derivation directly: `pullback()` applies
`reconstruction_.Tvmult(...)`.

### 4.1 A derivative and a state vector can have the same shape without being the same kind of object

This point is worth making explicit.

Both

$$
\mathbf z
\in
\mathbb R^{n}
$$

and

$$
\mathbf r_{\mathrm{ind}}
\in
\mathbb R^{n}
$$

contain $n$ numbers.

But $\mathbf z$ represents a state-like object: it tells us which linear combination
of basis functions to form.

The derivative vector represents a functional:

$$
\delta\mathbf z
\longmapsto
\mathbf r_{\mathrm{ind}}^{\mathsf T}\delta\mathbf z.
$$

Its natural operation is therefore to **act on** a perturbation.

That distinction is what the project calls the distinction between a **primal**
value and a **covector**.

The terminology is broader than the usual phrase "primal optimization variable".
An adjoint variable, for example, can still be represented as a primal value in its
own test space. "Primal" here means an element of a numerical vector space, as opposed
to an element of its dual.

## 5. Coordinate vectors and covectors

To make that distinction more systematic, let

$$
X_{h}
:=
\mathrm{span}\{\phi_{1},\ldots,\phi_{n}\}.
$$

A primal field

$$
x_{h}
=
\sum_{i=1}^{n}x_{i}\phi_{i}
$$

is represented by the coordinate vector

$$
\mathbf x
:=
(x_{1},\ldots,x_{n})^{\mathsf T}.
$$

Now let

$$
\lambda\in X_{h}^{\ast}
$$

be a linear functional.

Its coordinate representation can be defined by evaluating it on the basis:

$$
\lambda_{i}
:=
\lambda(\phi_{i}).
$$

Collect those values into

$$
\boldsymbol\lambda
:=
(\lambda_{1},\ldots,\lambda_{n})^{\mathsf T}.
$$

Then for

$$
x_{h}
=
\sum_{i=1}^{n}x_{i}\phi_{i},
$$

linearity gives

```math
\begin{aligned}
\lambda(x_{h})
&=
\sum_{i=1}^{n}x_{i}\lambda(\phi_{i})\\
&=
\boldsymbol\lambda^{\mathsf T}\mathbf x.
\end{aligned}
```

This is the coefficient dual pairing.

No mass matrix has appeared.

That is not because we have secretly chosen the Euclidean $L^{2}$ inner product. It
is because $\boldsymbol\lambda$ already contains the coefficients of a **functional**
in the dual coordinate representation.

### 5.1 Where the mass matrix does appear

Suppose instead that we begin with a primal field

$$
g_{h}
=
\sum_{i=1}^{n}g_{i}\phi_{i}
$$

and want the functional

$$
\lambda(v_{h})
:=
(g_{h},v_{h})_{L^{2}(\Omega)}.
$$

Writing

$$
v_{h}
=
\sum_{j=1}^{n}v_{j}\phi_{j},
$$

we obtain

```math
\lambda(v_{h})
=
\mathbf g^{\mathsf T}M\mathbf v,
```

where

$$
M_{ij}
:=
\int_{\Omega}\phi_{j}\phi_{i}.
$$

Therefore the covector coordinates of $\lambda$ are

$$
\boldsymbol\lambda
=
M^{\mathsf T}\mathbf g
=
M\mathbf g
$$

for the symmetric mass matrix.

The pairing is still simply

$$
\boldsymbol\lambda^{\mathsf T}\mathbf v.
$$

The mass matrix appeared earlier, when the primal field $\mathbf g$ was converted into
the dual representation $\boldsymbol\lambda$.

This is exactly the Riesz-map distinction introduced at the end of the previous
chapter:

```text
primal coefficients g
        │
        │ L2 Riesz map
        ▼
covector coefficients lambda = M g
        │
        │ dual pairing with v
        ▼
scalar lambda^T v
```

That separation is central to nmopt's representation.

### 5.2 Why `dot(covector, primal)` is not an $L^{2}$ inner product

The contract layer eventually computes coefficient pairings using the backend dot
product.

Conceptually,

$$
\mathrm{pair}(\boldsymbol\lambda,\mathbf v)
=
\boldsymbol\lambda^{\mathsf T}\mathbf v.
$$

Calling the native vector dot product is therefore appropriate at this stage.

It would be incorrect to insert another mass matrix automatically:

$$
\boldsymbol\lambda^{\mathsf T}M\mathbf v
$$

would apply an extra identification that does not belong to the dual pairing.

This explains a source-code detail that can otherwise be surprising. In
`layout.hpp`, `pair(covector, primal)` loops over blocks and calls `Backend::dot`.
The mathematical meaning comes from the fact that one operand has already been
declared and constructed as a covector.

The backend dot product supplies coefficient contraction, not the optimization
metric.

## 6. Residuals are covectors too

Return to the weak residual from the first chapter.

For a test space

$$
Z_{h}
:=
\mathrm{span}\{\zeta_{1},\ldots,\zeta_{m}\},
$$

the discrete residual is a functional

$$
E_{h}(x_{h})
\in
Z_{h}^{\ast}.
$$

Its coefficient representation is naturally

$$
(\mathbf r_{E})_{i}
:=
\langle E_{h}(x_{h}),\zeta_{i}\rangle.
$$

Given a test function

$$
p_{h}
=
\sum_{i=1}^{m}p_{i}\zeta_{i},
$$

the residual acts as

$$
\langle E_{h}(x_{h}),p_{h}\rangle
=
\mathbf r_{E}^{\mathsf T}\mathbf p.
$$

So the residual coefficient array and the test coefficient array may have the same
length, yet occupy opposite sides of the pairing:

```text
test vector p            element of Z_h
residual r_E             element of Z_h*
pairing                  r_E^T p
```

This becomes especially important in adjoint calculations.

A transpose residual action takes a **primal test-space seed**

$$
p_{h}\in Z_{h}
$$

and returns a **covector in the variable space**

$$
E_{h}'(x_{h})^{\ast}p_{h}
\in
X_{h}^{\ast}.
$$

That is why, later in the source, an adjoint/test seed can be represented by
`PrimalBlockT` even though the word "adjoint" may suggest "dual" in another context.

The mathematical question is not whether a variable is called state, control, or
adjoint. It is which space it belongs to and whether it is currently being used as an
element of that space or as a linear functional on it.

Chapter 3, [Operators, derivatives, and adjoints](03-operators-derivatives-and-adjoints.md), develops this point through JVP and VJP actions.

## 7. Trial and test spaces may coincide without being the same role

For the Poisson example,

$$
Y_{h}=Z_{h}
$$

is a natural Galerkin choice: state and test functions come from the same finite
element space.

It is still useful to distinguish the two roles.

The state coefficient vector answers:

> Which state field are we evaluating?

The test coefficient vector answers:

> Against which test direction are we evaluating the residual?

In a Petrov–Galerkin method the distinction becomes structural rather than merely
role-based:

$$
Y_{h}\neq Z_{h}.
$$

The two spaces may have different bases or even different dimensions.

Keeping a separate test-space identity at the contract boundary therefore avoids
baking a Galerkin-specific coincidence into a more general numerical interface.

This explains why the executable model exposes both a variable layout and a test
layout even when a particular application happens to give them equal dimensions.

## 8. Product spaces explain why nmopt uses blocks

The optimization problem rarely contains only one field.

For the distributed-control example,

$$
X_{h}
:=
Y_{h}\times U_{h}.
$$

A point in the full variable space is

$$
x_{h}
=
(y_{h},u_{h}).
$$

After choosing coordinates in each factor, one could flatten everything into a single
long vector,

$$
\mathbf x
=
\begin{bmatrix}
\mathbf y\\
\mathbf u
\end{bmatrix}.
$$

That is perfectly valid algebraically. But the block structure carries useful
information that would otherwise need to be reconstructed from offsets:

```text
block 0
    state coordinates in Y_h

block 1
    control coordinates in U_h
```

The derivative of a scalar functional on this product space has the corresponding
dual decomposition

$$
X_{h}^{\ast}
=
Y_{h}^{\ast}\times U_{h}^{\ast}.
$$

So a full objective derivative can be viewed as

$$
J_{h}'(x_{h})
=
\left(
D_{y}J_{h},
D_{u}J_{h}
\right).
$$

In coordinates,

$$
\mathbf r_{J}
=
\begin{bmatrix}
\mathbf r_{y}\\
\mathbf r_{u}
\end{bmatrix},
$$

where each block is a covector for the corresponding primal block.

This product-space viewpoint is why a block-oriented runtime representation is
natural for the project.

The particular class `BlockLayout`, however, is still an **implementation choice**.
The mathematics gives us product spaces. nmopt chooses to preserve selected
information about those products as ordered runtime blocks.

### 8.1 Blocks preserve relationships that offsets alone hide

Suppose

$$
\dim Y_{h}=225,
\qquad
\dim U_{h}=289.
$$

A flattened coordinate array would have dimension 514.

Without additional metadata, the number 514 tells us nothing about where the state
ends and the control begins, nor what either range means.

A block representation carries that structure explicitly:

```text
variables
├── state       225 coefficients
└── control     289 coefficients
```

This becomes even more useful for all-at-once formulations, where state, control,
adjoint, multipliers, or other fields may coexist in one product. Chapters 7 and 8
return to richer block structures.

For the reduced formulation, the two-block state/control example is enough to see the
idea.

## 9. Three notions of "space" now coexist in the repository

Before introducing `BlockLayout`, it is worth separating it from two richer notions
that appear elsewhere in the project.

### Mathematical or semantic space

At the problem-description level, a space may carry information such as:

- topology: $H^{1}$, $L^{2}$, trace-like, finite-dimensional Euclidean, and so on;
- role: state, test, control, observation, auxiliary;
- associated region;
- analytical or realization requirements.

The semantic model needs this information because it participates in deciding what
problem has been declared and whether a supported realization can be built.

### Concrete discrete space

After numerical realization, a finite-element space involves much more:

- mesh or triangulation;
- finite-element family and polynomial degree;
- degree-of-freedom numbering;
- affine constraints;
- boundary partition;
- basis functions;
- quadrature policy;
- possibly transfer or trace machinery.

This is the level at which `DoFHandler`, `FiniteElement`, and
`AffineConstraints` live in deal.II.

### Contract layout

At the solver/formulation boundary, most of that detail is no longer needed.

The generic numerical code mainly needs to know:

- how many blocks a value has;
- which runtime space each block claims to represent;
- how many coefficients belong to each block;
- the block ordering.

That smaller record is what `BlockLayout` supplies.

The relationship can be summarized as:

```text
semantic space declaration
    rich problem meaning
            │
            ▼
concrete discrete realization
    mesh + FE + constraints + coordinates + operators
            │
            ▼
BlockLayout
    small runtime identity used at numerical boundaries
```

A layout is therefore not a finite-element space object, and it is not a proof that
two independently constructed spaces are mathematically identical.

It is a producer-assigned runtime description used by the already-constructed
numerical services to keep values from obviously incompatible spaces from being
composed accidentally.

That limited role is important to understand.

## 10. `BlockLayout` as a runtime descriptor

The current implementation stores three pieces of information:

```cpp
std::string              label_;
std::vector<SpaceId>     spaces_;
std::vector<std::size_t> dimensions_;
```

For every block, the layout has one runtime `SpaceId` and one positive dimension.

The block order matters.

Two layouts are considered compatible when their ordered space IDs and dimensions
are equal:

```cpp
return spaces_ == other.spaces_ &&
       dimensions_ == other.dimensions_;
```

The human-readable label is not part of compatibility.

### 10.1 Why a space ID is needed in addition to dimension

Consider two one-block layouts:

```text
state
    SpaceId "state"
    dimension 289

control
    SpaceId "control"
    dimension 289
```

A size check alone would accept either vector wherever the other is expected.

The space IDs keep them distinct.

This is especially valuable because native storage types do not encode such meaning.
Both may be `dealii::Vector<double>`.

A `SpaceId` is deliberately small: currently it is essentially a string token. It
does not know the mesh, basis, topology, or metric.

Its purpose is identity at the numerical contract boundary, not semantic
introspection.

### 10.2 Compatibility is stronger than shape, weaker than mathematical equivalence

This gives `BlockLayout::compatible_with()` a very specific meaning.

If two layouts are compatible, nmopt knows that:

- they have the same number of blocks;
- the corresponding blocks carry the same producer-assigned space IDs;
- the corresponding coefficient dimensions agree.

It does **not** independently verify that:

- the bases are the same;
- the meshes are the same;
- the DoF numbering is the same;
- the affine constraints are the same;
- the two services came from the same numerical realization.

Those deeper facts must already be guaranteed by the producer that created the
layout and associated operations.

Chapters 11 and 13 use additional evidence for stronger notions of identity when it
is needed.

### 10.3 A layout label is for humans

The Step-4 Problem B binding constructs a variable layout labeled

```text
external_step4_minimal_problem_b_variables
```

and a test layout labeled

```text
external_step4_minimal_problem_b_test
```

The labels make diagnostics understandable.

Inside the variable layout, however, the structural IDs are simply:

```text
state
control
```

with dimensions 225 and 289.

The test layout contains:

```text
state_test
```

with dimension 225.

State and test happen to have the same dimension, but the distinct IDs keep the
contract roles separate.

### 10.4 `single_block()` creates a view of one factor

A formulation often needs to talk about the control space alone even when the
executable model uses the product layout

$$
X_{h}=Y_{h}\times U_{h}.
$$

`BlockLayout::single_block()` creates a one-block layout that preserves the selected
space ID and dimension while giving the new view its own descriptive label.

The reduced state/control partition uses this mechanism to obtain state and control
layouts from the full variable layout.

Conceptually:

```text
variables [ state | control ]
             │         │
             ▼         ▼
        state view   control view
```

No coordinate transformation occurs here. The operation only selects the identity
record for one factor of the product.

## 11. `PrimalBlockT` and `CovectorBlockT` attach role to storage

A layout tells us *which discrete spaces* the blocks represent.

We still need to distinguish an element of those spaces from a functional acting on
them.

The contract layer uses two types:

```text
PrimalBlockT<Backend>
CovectorBlockT<Backend>
```

Both ultimately contain native backend vectors and a `LayoutPtr`.

Their storage can therefore look identical.

The distinction is in the C++ type and in how the value is allowed to participate in
the numerical interfaces.

### 11.1 "Primal" does not mean "the optimization primal variable"

This naming deserves care.

In these contracts, a primal value is an element of a declared numerical space.

Examples include:

- a state;
- a control;
- a perturbation direction;
- a test vector;
- an adjoint represented in the test space.

So an adjoint $`p\in Z_{h}`$ can be a `PrimalBlockT` with the test layout.

A covector represents an element of the corresponding dual space.

Examples include:

- a residual in $`Z_{h}^{\ast}`$;
- an objective derivative in $`X_{h}^{\ast}`$;
- a reduced derivative in $`U_{h}^{\ast}`$.

This terminology becomes much less confusing if one asks "element or functional?"
rather than "primal variable or dual variable?"

### 11.2 One layout can describe both a space and its dual coordinate shape

The primal space $`X_{h}`$ and its dual $`X_{h}^{\ast}`$ have the same finite dimension.

nmopt therefore uses the same `BlockLayout` to describe the block identities and
dimensions of both a `PrimalBlockT` and a `CovectorBlockT`.

The distinction between $`X_{h}`$ and $`X_{h}^{\ast}`$ is carried by the wrapper type,
not by a second "dual layout".

That matches the coordinate mathematics:

```text
primal
    x ∈ X_h
    coordinates x ∈ R^n

covector
    lambda ∈ X_h*
    coordinates lambda ∈ R^n
```

Same coefficient count, different transformation and operational meaning.

### 11.3 Construction checks shape against the layout

The common `BlockValuesT` base checks that:

- the layout exists;
- the number of native vectors matches the number of blocks;
- every native vector has the dimension declared for its block.

So the wrapper prevents a 224-entry state vector from being presented as a
225-dimensional state block.

The stronger role distinction then comes from using separate primal and covector
types in signatures.

## 12. Pairing makes the primal/covector distinction operational

The most direct place where the two types meet is

```cpp
pair(covector, primal)
```

Mathematically this represents

$$
\langle \lambda,x\rangle_{X_{h}^{\ast},X_{h}}.
$$

The implementation first requires compatible layouts. It then contracts the
corresponding backend coefficient vectors block by block.

For a two-block product space,

$$
X_{h}
=
Y_{h}\times U_{h},
$$

with

$$
\lambda
=
(\lambda_{y},\lambda_{u}),
\qquad
x
=
(y,u),
$$

the pairing is

$$
\langle \lambda,x\rangle
=
\langle \lambda_{y},y\rangle
+
\langle \lambda_{u},u\rangle.
$$

In coefficient form,

$$
=
\boldsymbol\lambda_{y}^{\mathsf T}\mathbf y
+
\boldsymbol\lambda_{u}^{\mathsf T}\mathbf u.
$$

This is why the block loop in `pair()` is mathematically natural.

### 12.1 Pairing is coordinate contraction, not geometry

The distinction from the metric chapter can now be stated precisely.

If $\lambda$ is already a covector, the dual pairing with $x$ is simply

$$
\boldsymbol\lambda^{\mathsf T}\mathbf x.
$$

If instead we start from a primal vector $g$ and want a functional defined by an
inner product, we must first apply a Riesz map:

$$
g
\overset{R}{\longmapsto}
\lambda
\overset{\text{pair with }x}{\longmapsto}
\lambda(x).
$$

For the finite-element $L^{2}$ geometry:

$$
\boldsymbol\lambda
=
M\mathbf g.
$$

Then

$$
\mathrm{pair}(\boldsymbol\lambda,\mathbf x)
=
\mathbf g^{\mathsf T}M\mathbf x.
$$

The mass matrix belongs to the first arrow, not the second.

That is the conceptual reason `pair()` itself does not know about metrics.

## 13. The Step-4 binding shows the representation in one place

The external Step-4 integration is useful because the native application mathematics
and the nmopt contract wrapping are both visible.

Problem B exposes:

```text
state_dimension()   = 225
control_dimension() = 289
```

The binding constructs the variable layout as

```cpp
BlockLayout(
  "external_step4_minimal_problem_b_variables",
  {SpaceId{"state"}, SpaceId{"control"}},
  {problem.state_dimension(), problem.control_dimension()});
```

and the test layout as

```cpp
BlockLayout(
  "external_step4_minimal_problem_b_test",
  {SpaceId{"state_test"}},
  {problem.state_dimension()});
```

This small piece of code now has a clear mathematical interpretation.

The variable product is

$$
X_{h}
=
Y_{h}^{\mathrm{ind}}
\times
U_{h},
$$

where

$$
\dim Y_{h}^{\mathrm{ind}}=225,
\qquad
\dim U_{h}=289.
$$

The residual test space has dimension 225.

The native `ProblemB` methods still receive ordinary `dealii::Vector<double>`
objects. The binding is where those native vectors acquire the framework's explicit
space and primal/covector roles.

### 13.1 Residual values become test-space covectors

The native application computes a residual vector from state and control coordinates.

The binding wraps that returned vector as

```text
Covector(test_layout, ...)
```

because the residual represents an element of

$$
Z_{h}^{\ast}.
$$

Likewise, the JVP returns another residual covector in the same test-space dual.

The native array alone did not carry that information. The wrapping operation does.

### 13.2 Test seeds are primal values

The VJP callback receives a `Primal` with the test layout.

That represents

$$
p\in Z_{h}.
$$

The native application computes the transpose action and returns state and control
components. The binding then wraps them as a covector with the full variable layout:

$$
E_{h}'(x)^{\ast}p
\in
X_{h}^{\ast}.
$$

So one callback visibly crosses

```text
test primal
    ↓ transpose derivative action
variable covector
```

Chapter 3 derives that operation in detail.

### 13.3 Objective derivatives use the same variable layout as the point

The objective consumes a primal variable point

$$
x=(z,u)\in X_{h}.
$$

Its derivative lies in

$$
X_{h}^{\ast}.
$$

The binding therefore uses the same variable layout but changes the wrapper type:

```text
Primal(variable_layout)
    point (z,u)

Covector(variable_layout)
    objective derivative (D_z J, D_u J)
```

This is a concise example of why "same layout" and "same mathematical role" are
different statements.

## 14. Physical coordinates can disappear before the generic solver sees them

The Step-4 binding's state layout has dimension 225, not 289.

The optimizer therefore never needs to know that the physical state field has 289
coefficients.

That reconstruction remains inside Problem B:

```text
nmopt state block z (225)
        │
        ▼
Problem B coordinate map
        │
        ├─ embed perturbation: P z
        │
        └─ reconstruct field:  P z + ell
        ▼
native Step-4 physical vectors (289)
```

This is a useful example of abstraction that removes *coordinates*, not mathematics.

The physical boundary condition has not disappeared. It has been built into the
coordinate map.

The generic reduced formulation can work entirely in independent coordinates because
Problem B translates every operation appropriately.

### 14.1 The weak equation also transforms under the coordinate map

Suppose, for illustration, that the full physical residual is

$$
\mathbf R_{\mathrm{phys}}
(\mathbf y_{\mathrm{phys}},\mathbf u)
:=
A\mathbf y_{\mathrm{phys}}
-
\mathbf b
-
M\mathbf u.
$$

Admissible test perturbations are also generated by the homogeneous embedding $P$.
If

$$
\mathbf w\in\mathbb R^{n}
$$

contains independent test coordinates, then the physical test vector is

$$
P\mathbf w.
$$

The weak residual equation requires

$$
(P\mathbf w)^{\mathsf T}
\mathbf R_{\mathrm{phys}}
=
0
\qquad
\text{for every }\mathbf w.
$$

Equivalently,

$$
P^{\mathsf T}
\mathbf R_{\mathrm{phys}}
=
0.
$$

Substitute

$$
\mathbf y_{\mathrm{phys}}
=
P\mathbf z+\boldsymbol\ell:
$$

```math
P^{\mathsf T}
\left(
A(P\mathbf z+\boldsymbol\ell)
-
\mathbf b
-
M\mathbf u
\right)
=
0.
```

Grouping terms gives the independent-coordinate system

```math
K\mathbf z
=
\mathbf b_{F}
+
B\mathbf u,
```

with the conceptual definitions

```math
\begin{aligned}
K &:= P^{\mathsf T}AP,\\
B &:= P^{\mathsf T}M,\\
\mathbf b_{F}
&:=
P^{\mathsf T}
\left(
\mathbf b-A\boldsymbol\ell
\right).
\end{aligned}
```

The exact assembly details of an existing application can differ – Step-4 already
contains its own boundary-treated system and lifting corrections – but this reduction
explains the algebraic roles recorded in the case study:

- $P$ maps independent state coordinates to the physical field;
- $P^{\mathsf T}$ pulls physical residual/derivative information back to independent
  coordinates;
- $B=P^{\mathsf T}M$ maps the full control field into the independent state
  equations.

This one calculation connects coordinate reconstruction, weak testing, and dual
pullback.

### 14.2 Why output uses a different representation from optimization

When optimization finishes, the retained state is still an independent vector
$\mathbf z$.

VTK output, physical norms, and native application visualization expect the physical
field

$$
P\mathbf z+\boldsymbol\ell.
$$

So output reconstructs the full state.

This is not merely a data-format conversion. It reflects two legitimate
representations of the same constrained field:

```text
independent coordinates
    best suited to solving/optimization

physical FE coefficients
    best suited to field evaluation/output
```

Chapter 11, [Compilation and lowering](11-compilation-and-lowering.md), revisits this distinction on the compiler-created path.

## 15. The reusable coordinate class generalizes the Step-4 example

`ProblemBCoordinates` is intentionally simple: it knows a set of fixed boundary
indices and performs restriction/insertion.

The reusable `IndependentStateCoordinates` implementation in
`include/nmopt/dealii/independent_state_coordinates.hpp` is more general.

Its homogeneous `AffineConstraints` can encode linear relations among DoFs. The
reconstruction matrix it builds therefore represents the actual constrained
coordinate embedding, not only omission of fixed entries.

This distinction matters on adaptively refined meshes, where hanging-node
constraints are common.

The mathematical pattern remains the same:

```math
\mathbf y_{\mathrm{phys}}
=
P\mathbf z+\boldsymbol\ell.
```

What changes is the structure of $P$.

For simple fixed boundary values, $P$ looks like an insertion matrix.

For general homogeneous affine constraints, constrained physical entries can be
linear combinations of independent coordinates.

The same transpose pullback remains correct in both cases.

This is one reason it is preferable to think in terms of a coordinate map $P$ rather
than memorize "remove boundary entries".

## 16. Product blocks and coordinate maps solve different problems

At this point two structural devices have appeared:

- the coordinate map $P$;
- the runtime `BlockLayout`.

They should not be conflated.

The coordinate map answers:

> How does one representation of a field map into another representation of that
> same field?

For example,

$$
\mathbf z
\longmapsto
P\mathbf z+\boldsymbol\ell.
$$

The layout answers:

> Which runtime spaces and product factors do these native coefficient vectors claim
> to represent?

For example,

```text
variables
├── state     225
└── control   289
```

A coordinate map can change dimension and coefficients.

A layout does not transform values. It describes them.

This distinction is useful when navigating the code:

```text
IndependentStateCoordinates
    numerical coordinate transformation

BlockLayout
    runtime identity metadata

PrimalBlockT / CovectorBlockT
    role-bearing containers
```

The first belongs to the concrete deal.II realization.

The latter two sit at the generic numerical contract boundary.

## 17. What the layout deliberately forgets

The small size of `BlockLayout` is easiest to appreciate by listing what it does *not*
contain.

Given a control block with

```text
SpaceId "control"
dimension 289
```

the layout does not tell us whether the coefficients represent:

- a volume $`\mathbb P_{1}`$ field;
- a $`\mathbb Q_{1}`$ field on quadrilaterals;
- a boundary trace;
- a cellwise-constant field;
- a parameter vector unrelated to finite elements.

It also does not tell us:

- the mesh;
- the basis ordering;
- the mass matrix;
- the boundary constraints;
- the metric;
- the control-to-state coupling.

Those facts live in the producer's numerical realization.

Why keep such a shallow descriptor at all?

Because generic formulation and solver code still benefits from catching errors such
as:

```text
225-state coefficients passed where 289-control coefficients are expected

a test-space vector passed where a control-space vector is expected

a one-block value passed where a two-block variable point is expected
```

The layout carries enough identity for that purpose without dragging deal.II objects
into generic solver interfaces.

This is a recurring theme in the project: a boundary object often carries only the
information required by the consumer at that boundary.

## 18. A compact correspondence

The main distinctions in this chapter can be summarized as follows:

| Question | Mathematical answer | Concrete numerical object | nmopt-facing representation |
| --- | --- | --- | --- |
| What field is this? | $`x_{h}\in X_{h}`$ | coefficients in a chosen basis | `PrimalBlockT` |
| What functional is this? | $`\lambda\in X_{h}^{\ast}`$ | basis evaluations of the functional | `CovectorBlockT` |
| How do they interact? | $\langle\lambda,x\rangle$ | $\boldsymbol\lambda^{\mathsf T}\mathbf x$ | `pair(covector, primal)` |
| Which product factor is this? | $`X_{h}=Y_{h}\times U_{h}`$ | block ordering and sizes | `BlockLayout` |
| How is a constrained state represented? | $`y_{\mathrm{phys}}=Pz+\ell`$ | reconstruction matrix + lifting | `IndependentStateCoordinates` or application-owned equivalent |
| How does a physical covector move to independent coordinates? | pullback through $P$ | $P^{\mathsf T}r$ | `pullback()` / equivalent native operation |

The table compresses the chapter, but the direction is important: the rightmost
column should be read as a representation of the mathematical/numerical objects to
its left.

## 19. Following these ideas through the source

A useful source tour now has a much narrower purpose than the tour in the previous
chapter.

### Runtime block representation

Read:

- [`include/nmopt/contract/layout.hpp`](../../../include/nmopt/contract/layout.hpp)

This one header contains:

- `SpaceId`;
- `BlockLayout`;
- `PrimalBlockT`;
- `CovectorBlockT`;
- `pair(...)`;
- block extraction helpers.

The file is small enough to read in one sitting once the ideas in this chapter are
clear.

### Generic deal.II independent coordinates

Read:

- [`include/nmopt/dealii/independent_state_coordinates.hpp`](../../../include/nmopt/dealii/independent_state_coordinates.hpp)

Look for the three operations:

```text
embed
reconstruct
pullback
```

and relate them directly to

```math
Pz,\qquad
Pz+\ell,\qquad
P^{\mathsf T}r.
```

The code that constructs the columns of `reconstruction_` by distributing
`AffineConstraints` is the implementation detail that turns deal.II constraint
relations into the matrix $P$.

### A simpler application-owned coordinate map

Read:

- [`apps/external-dealii/step-4/integration/problem_b_coordinates.hpp`](../../../apps/external-dealii/step-4/integration/problem_b_coordinates.hpp)

This is the fixed-boundary special case. Its `restrict`, `reconstruct`, and
`embed_free` operations make the coordinate distinction particularly easy to see.

Then read:

- [`apps/external-dealii/step-4/integration/problem_b.hpp`](../../../apps/external-dealii/step-4/integration/problem_b.hpp)

Notice where physical reconstruction is required for the objective and where free
coordinates are used for state solves and residual actions.

### Where native values acquire contract meaning

Finally read:

- [`apps/external-dealii/step-4/minimal/problem_b_binding.hpp`](../../../apps/external-dealii/step-4/minimal/problem_b_binding.hpp)

The layout construction near the top of `ProblemBBinding` is the point where native
state/control/test coordinates receive the IDs and block dimensions used by generic
nmopt code.

The callback definitions then show native vectors being wrapped as either primal
values or covectors according to their mathematical role.

## 20. What this chapter has not yet explained

We now have enough language to state the next question precisely.

A point

$$
x\in X_{h}
$$

is represented by a primal block value.

A residual

$$
E_{h}(x)\in Z_{h}^{\ast}
$$

is represented by a covector in the test layout.

But an optimizer and adjoint method need more than residual values. They need the
linearized map

$$
E_{h}'(x):X_{h}\longrightarrow Z_{h}^{\ast}
$$

and its transpose

$$
E_{h}'(x)^{\ast}:Z_{h}\longrightarrow X_{h}^{\ast}.
$$

Why are these maps represented as JVP and VJP actions rather than matrices? What
exact pairing defines the transpose? How should their dimensions and block structure
be read? How do finite-difference and transpose tests verify them? And how does the
same language cover nonlinear models?

Those are the subjects of
**Operators, derivatives, and adjoints**.

## Read later

Several existing documents cover related material from different viewpoints:

- [Theoretical formalism](../../design/mathematical-model.md) states the project's
  abstract space, duality, derivative, transformation, and adjoint conventions.
- [Integrating an existing PDE application](../overview/external-applications.md)
  places the Step-4 coordinate example in the broader external-application path.
- [Numerical realization](../overview/numerical-realization.md) gives the high-level
  role of coordinate maps and deal.II services.
- [Step-4 external integration overview](../../../apps/external-dealii/step-4/external-integration-overview.md)
  follows Problem B through its native application and optimizer connection.
- [External deal.II solver integration](../../reference/external-dealii-solver-integration.md)
  gives the exact public contract signatures and lifetime rules.

[Numerical realization](../overview/numerical-realization.md) returns to
`IndependentStateCoordinates` from the deal.II side, while Chapter 4,
[Metrics, gradients, and constraints](04-metrics-gradients-and-constraints.md),
develops the Riesz-map side of the primal/covector distinction.
