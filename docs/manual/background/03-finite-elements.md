# 03 · Finite elements

**Background navigation:** [Index](README.md) \
Previous: [02 · Duality, derivatives, and adjoints](02-duality-derivatives-and-adjoints.md) \
Next: [04 · Numerical linear algebra for PDE and optimization systems](04-numerical-linear-algebra.md) \
Companion: [03a · Further finite-element notes](03a-further-finite-element-notes.md)

## Purpose

A weak PDE is still an infinite-dimensional problem. The finite-element method turns it into a finite-dimensional one without discarding the variational structure that made the weak formulation useful in the first place.

The central idea is simple: replace the infinite-dimensional trial and test spaces by carefully constructed finite-dimensional spaces, then require the weak equation to hold for every discrete test function. The practical machinery behind that sentence is what this chapter develops. We will see how local polynomial pieces become a global function space, how degrees of freedom identify a discrete function, how element integrals become local matrices, and how local contributions are assembled into one sparse algebraic system.

The chapter also separates several notions that are often collapsed in informal discussion. **Galerkin approximation** is the act of restricting a variational problem to finite-dimensional spaces. A **finite element** is a local construction used to build those spaces. A **mesh** organizes the local cells. A **degree of freedom** is a linear functional that identifies coefficients. A coefficient vector is a representation of a finite-element function, not the function itself.

The scope is deliberately narrower than a full finite-element course. We derive the mechanics needed to understand conforming scalar finite elements, mass and stiffness matrices, numerical quadrature, local-to-global assembly, Dirichlet and affine constraints, and different state/control discretizations. We will state and derive the basic Galerkin best-approximation result, but we stop before detailed interpolation estimates, Bramble–Hilbert arguments, a posteriori error estimation, discontinuous Galerkin methods, mixed finite elements, and inf-sup theory. The companion [03a · Further finite-element notes](03a-further-finite-element-notes.md) develops selected extensions that are especially useful later, including interpolation and projection of continuous data, quadrature formulas and consistency errors, mesh regularity, adaptivity, and mixed formulations.

## Before you start

This chapter assumes:

- [01 · Function spaces and weak PDEs](01-function-spaces-and-weak-pdes.md), especially weak formulations and the spaces $L^{2}(\Omega)$, $H^{1}(\Omega)$, and $H_{0}^{1}(\Omega)$;
- elementary polynomial interpolation in one variable;
- ordinary matrix and vector algebra.

The vector/covector distinctions from [02 · Duality, derivatives, and adjoints](02-duality-derivatives-and-adjoints.md) are useful, especially when matrices are interpreted as coordinate representations of operators, but they are not required for the first half of the chapter.

## What you will be able to do

After this chapter, you should be able to:

- explain the Galerkin approximation of a weak problem;
- distinguish a finite-dimensional approximation space from the local finite elements used to construct it;
- define a local finite element through a cell, a local function space, and degrees of freedom;
- explain what unisolvence means and why it produces a local basis;
- distinguish $\mathbb P_{k}$ and $\mathbb Q_{k}$ polynomial families at a working level;
- explain why shared degrees of freedom produce globally continuous finite-element functions;
- distinguish continuous nodal spaces from discontinuous cellwise spaces;
- map basis functions and integrals between reference and physical cells;
- derive local stiffness and mass matrices for the one-dimensional $\mathbb P_{1}$ element;
- explain what a quadrature rule approximates and how quadrature enters element assembly;
- derive local-to-global assembly and explain the resulting sparsity pattern;
- explain how prescribed Dirichlet values and more general affine constraints reduce the independent coordinates;
- distinguish mass, stiffness, and rectangular coupling matrices by the bilinear forms and spaces they represent;
- reconstruct a physical finite-element field from its coefficients and avoid identifying the field with the coefficient vector;
- explain, at a high level, why a conforming Galerkin method converges when its spaces approximate the continuous solution space increasingly well.

## Roadmap

We begin by restricting the weak Poisson problem to a finite-dimensional subspace. A one-dimensional piecewise-linear example makes basis functions and coefficients concrete. We then step backward and explain how such spaces are built locally: finite-element triplets, Lagrange degrees of freedom, $\mathbb P_{k}$ and $\mathbb Q_{k}$ families, meshes, conformity, and local-to-global degree-of-freedom numbering.

Once the spaces are clear, we return to computation. Reference-cell maps let us reuse the same local formulas on every physical cell. Quadrature turns element integrals into finite sums. Local matrices and vectors are then accumulated into one global sparse system. The final sections treat affine constraints, different state and control spaces, coefficient-versus-field distinctions, and the basic Galerkin convergence argument.

The running PDE is the homogeneous-Dirichlet Poisson problem from [01 · Function spaces and weak PDEs](01-function-spaces-and-weak-pdes.md). A cellwise-constant control space is introduced later because it exposes a second important finite-element pattern without requiring discontinuous Galerkin theory.

## 1. Galerkin approximation: restrict the weak problem

Consider the weak Poisson problem on a bounded domain $\Omega$:

> find $y\in V:=H_{0}^{1}(\Omega)$ such that
>
> $$
> a(y,v)=F(v)
> \qquad
> \text{for every }v\in V,
> $$
>
> where
>
> $$
> a(y,v)
> :=
> \int_{\Omega}\nabla y\cdot\nabla v \mathrm{d}x,
> \qquad
> F(v)
> :=
> \int_{\Omega}fv \mathrm{d}x.
> $$

The unknown $y$ belongs to an infinite-dimensional function space. A Galerkin method chooses a finite-dimensional subspace

$$
V_{h}\subset V
$$

and asks for a discrete function $y_{h}\in V_{h}$ satisfying the same variational equation against every discrete test function:

$$
a(y_{h},v_{h})
=
F(v_{h})
\qquad
\text{for every }v_{h}\in V_{h}.
$$

The subscript $h$ usually refers to a spatial resolution parameter. For a mesh-based method, $h$ is typically related to the largest cell diameter. Smaller $h$ means that the discrete space is built on a finer mesh.

Two changes have happened and only two:

```text
continuous trial space V      → finite-dimensional trial space V_h
continuous test space V       → finite-dimensional test space V_h
```

The bilinear form $a$ and functional $F$ have not yet been approximated. In this ideal Galerkin statement, the same integrals are used as in the continuous weak problem.

Because $V_{h}\subset V$, the exact solution $y$ also satisfies

$$
a(y,v_{h})=F(v_{h})
\qquad
\text{for every }v_{h}\in V_{h}.
$$

Subtracting the discrete equation gives

$$
a(y-y_{h},v_{h})=0
\qquad
\text{for every }v_{h}\in V_{h}.
$$

This is **Galerkin orthogonality**. It is the structural relation behind the basic error estimate derived near the end of the chapter.

At this stage nothing specifically says “finite element.” Any finite-dimensional subspace could be used. The finite-element method is a particular way to construct useful spaces $V_{h}$ from local pieces.

## 2. A first discrete space: piecewise-linear functions in one dimension

Take

$$
\Omega=(0,1)
$$

and divide the interval into cells

$$
K_{e}
:=
[x_{e},x_{e+1}],
\qquad
0=x_{0}<x_{1}<\cdots<x_{N}=1.
$$

For the moment assume a uniform mesh,

$$
h:=x_{e+1}-x_{e}=\frac{1}{N}.
$$

Define

```math
V_{h}
:=
\left\{
 v_{h}\in C^{0}([0,1])
 :
 v_{h}\rvert_{K_{e}}\text{ is affine for every }K_{e},
 \quad
 v_{h}(0)=v_{h}(1)=0
\right\}.
```

Every function in $V_{h}$ is affine on each cell but need not have the same slope on neighboring cells. The function is continuous at the cell interfaces, while its derivative may jump there.

The interior nodes

$$
x_{1},\ldots,x_{N-1}
$$

provide a natural set of coordinates. For each interior node $x_{i}$ define the hat function $\varphi_{i}\in V_{h}$ by

$$
\varphi_{i}(x_{j})=\delta_{ij}.
$$

Each $\varphi_{i}$ is piecewise affine, equals one at $x_{i}$, equals zero at every other mesh node, and vanishes outside the two cells that meet at $x_{i}$.

Any $v_{h}\in V_{h}$ can then be written uniquely as

$$
v_{h}(x)
=
\sum_{i=1}^{N-1}v_{i}\varphi_{i}(x).
$$

Because the basis is nodal,

$$
v_{i}=v_{h}(x_{i}).
$$

The coefficient vector is

$$
\mathbf v
:=
(v_{1},\ldots,v_{N-1})^{\mathsf T}
\in\mathbb R^{N-1}.
$$

This is the first concrete finite-element representation:

```text
piecewise-polynomial function v_h
        ↕ choose nodal basis
coefficient vector v
```

The function and vector contain equivalent information once the mesh and basis are fixed, but they are not the same mathematical object. The function can be evaluated at arbitrary physical points and differentiated cell by cell. The vector is a coordinate representation relative to a particular basis.

### Why continuity matters for an $H^{1}$-conforming method

The continuous Poisson problem is posed in

$$
V=H_{0}^{1}(0,1).
$$

A discrete space $V_{h}$ is called **conforming** for this problem when

$$
V_{h}\subset V.
$$

That inclusion matters because the discrete problem is then literally the continuous variational problem restricted to a smaller space. In particular, every discrete trial and test function is admissible in the continuous weak problem, the exact solution can be tested against functions in $V_{h}$, and the Galerkin orthogonality derived in the previous section follows without introducing an additional consistency error.

For the present Poisson problem, conformity therefore means

$$
V_{h}\subset H_{0}^{1}(0,1).
$$

The functions in $V_{h}$ are affine inside each cell, so they are classically differentiable there. Their slopes may jump from one cell to the next. We now show why continuity of the function itself is enough to obtain an $H^{1}$ function.

Let $v_{h}$ be continuous and affine on every cell. For any test function $\phi\in C_{c}^{\infty}(0,1)$, integrate by parts cell by cell:

```math
\begin{aligned}
-\int_{0}^{1}v_{h}\phi'\mathrm{d}x
&=
-\sum_{e=0}^{N-1}
\int_{x_{e}}^{x_{e+1}}v_{h}\phi'\mathrm{d}x
\\
&=
\sum_{e=0}^{N-1}
\int_{x_{e}}^{x_{e+1}}v_{h}'\phi\mathrm{d}x
-
\sum_{e=0}^{N-1}
\left[v_{h}\phi\right]_{x_{e}}^{x_{e+1}}.
\end{aligned}
```

The endpoint terms at $0$ and $1$ vanish because $\phi$ has compact support. At an interior node $x_{i}$, the contribution from the cell on the left and the contribution from the cell on the right combine to

$$
\left(
v_{h}(x_{i}^{+})-v_{h}(x_{i}^{-})
\right)\phi(x_{i}).
$$

Continuity gives

$$
v_{h}(x_{i}^{+})=v_{h}(x_{i}^{-}),
$$

so every interface contribution vanishes. Hence

$$
-\int_{0}^{1}v_{h}\phi'\mathrm{d}x
=
\int_{0}^{1}v_{h}'\phi\mathrm{d}x,
$$

where on the right $v_{h}'$ denotes the ordinary derivative taken separately on each cell. This identity is exactly the definition of a weak derivative, so the piecewise classical derivative is the weak derivative of $v_{h}$.

We still have to check that this weak derivative has the regularity required by $H^{1}$. Because $v_{h}$ is affine on every cell, $v_{h}'$ is piecewise constant. The mesh contains finitely many cells, so

$$
\int_{0}^{1}\lvert v_{h}'(x)\rvert^{2}\mathrm{d}x
=
\sum_{e=0}^{N-1}
\int_{K_{e}}\lvert v_{h}'(x)\rvert^{2}\mathrm{d}x
<\infty.
$$

Thus $v_{h}'\in L^{2}(0,1)$ and therefore

$$
v_{h}\in H^{1}(0,1).
$$

Finally, the definition of $V_{h}$ imposes

$$
v_{h}(0)=v_{h}(1)=0.
$$

Using the trace characterization of $H_{0}^{1}(0,1)$ from [01 · Function spaces and weak PDEs](01-function-spaces-and-weak-pdes.md), we conclude that

$$
v_{h}\in H_{0}^{1}(0,1).
$$

Therefore the whole discrete space is conforming:

$$
V_{h}\subset H_{0}^{1}(0,1).
$$

If $v_{h}$ jumped at an interface, the term

$$
\left(
v_{h}(x_{i}^{+})-v_{h}(x_{i}^{-})
\right)\phi(x_{i})
$$

would remain in the cellwise integration-by-parts formula. A general discontinuous piecewise-polynomial function therefore does not belong to $H^{1}$, even though it has an ordinary derivative inside every cell. Discontinuous finite-element methods compensate for those interface jumps explicitly rather than pretending that the space is $H^{1}$-conforming.

## 3. From a global space to local finite elements

The one-dimensional space above was easy to describe globally. In two and three dimensions, constructing a basis directly over the whole domain would be awkward. Finite elements reverse the perspective: define simple polynomial spaces and degrees of freedom on one cell, then glue copies of that local construction across a mesh.

A standard abstract definition describes a finite element by a triplet

$$
(K,\mathcal P,\Sigma).
$$

Here:

- $K$ is the **cell** or element domain, such as an interval, triangle, quadrilateral, tetrahedron, or hexahedron;
- $\mathcal P$ is a finite-dimensional space of functions on $K$, usually polynomials;
- $\Sigma=\{\sigma_{1},\ldots,\sigma_{m}\}$ is a set of $m=\dim\mathcal P$ linearly independent **degrees of freedom**, each $\sigma_{i}$ being a linear functional on $\mathcal P$.

The degrees of freedom identify a function in $\mathcal P$ by a finite list of numbers:

$$
p
\longmapsto
\left(
\sigma_{1}(p),\ldots,\sigma_{m}(p)
\right).
$$

For this identification to be useful, the degrees of freedom must determine $p$ uniquely. The finite element is called **unisolvent** when

$$
\sigma_{i}(p)=0
\quad
\text{for every }i
$$

implies

$$
p=0.
$$

Because $\mathcal P$ and the list of degrees of freedom have the same finite dimension, this is equivalent to saying that arbitrary admissible degree-of-freedom values determine one and only one function in $\mathcal P$.

Unisolvence gives a canonical local basis $\{\phi_{1},\ldots,\phi_{m}\}$ characterized by

$$
\sigma_{i}(\phi_{j})=\delta_{ij}.
$$

This is precisely the primal/dual basis relation from [02 · Duality, derivatives, and adjoints](02-duality-derivatives-and-adjoints.md). The shape functions $\phi_{j}$ live in the local function space; the degrees of freedom $\sigma_{i}$ are linear functionals that act on those functions.

For any $p\in\mathcal P$,

$$
p
=
\sum_{j=1}^{m}\sigma_{j}(p)\phi_{j}.
$$

The finite-element coefficients are therefore degree-of-freedom values whenever this nodal basis is used.

### Lagrange elements

The most familiar finite elements use point evaluation as their degrees of freedom. Choose support points

$$
a_{1},\ldots,a_{m}\in K
$$

and define

$$
\sigma_{i}(p):=p(a_{i}).
$$

The corresponding basis functions satisfy

$$
\phi_{j}(a_{i})=\delta_{ij}.
$$

These are **Lagrange finite elements**. Their coefficients are nodal values.

Not every finite element uses point values. Other families use moments, normal components, tangential components, or derivatives as degrees of freedom. Those cases matter in more advanced PDE discretizations, but point-value Lagrange elements are sufficient for the scalar $H^{1}$ problems studied here.

## 4. Polynomial families: $\mathbb P_{k}$ and $\mathbb Q_{k}$

Two polynomial families occur constantly in finite-element notation.

On a simplex – an interval in one dimension, triangle in two dimensions, tetrahedron in three – $\mathbb P_{k}$ denotes polynomials of **total degree** at most $k$.

For example, on a triangle with coordinates $(x,y)$,

$$
\mathbb P_{1}
=
\mathrm{span}\{1,x,y\},
$$

while

$$
\mathbb P_{2}
=
\mathrm{span}\{1,x,y,x^{2},xy,y^{2}\}.
$$

On a tensor-product reference cell such as a square or cube, $\mathbb Q_{k}$ denotes polynomials of degree at most $k$ **in each coordinate separately**. On a square,

$$
\mathbb Q_{1}
=
\mathrm{span}\{1,x,y,xy\}.
$$

The term $xy$ is why a $\mathbb Q_{1}$ function on a square is usually called **bilinear** rather than affine.

### The reference $\mathbb P_{1}$ triangle

On the reference triangle

$$
\widehat K
:=
\left\{(\xi,\eta):\xi\ge0,\ \eta\ge0,\ \xi+\eta\le1\right\},
$$

use the three vertices

$$
(0,0),\qquad(1,0),\qquad(0,1)
$$

as support points. The nodal basis is

$$
\widehat\phi_{1}(\xi,\eta)=1-\xi-\eta,
$$

$$
\widehat\phi_{2}(\xi,\eta)=\xi,
\qquad
\widehat\phi_{3}(\xi,\eta)=\eta.
$$

Each basis function is affine and takes the value one at one vertex and zero at the other two.

### The reference $\mathbb Q_{1}$ square

On

$$
\widehat K=[0,1]^{2},
$$

use the four vertices as support points. The tensor-product basis is

```math
\begin{aligned}
\widehat\phi_{1}(\xi,\eta)
&=(1-\xi)(1-\eta),
\\
\widehat\phi_{2}(\xi,\eta)
&=\xi(1-\eta),
\\
\widehat\phi_{3}(\xi,\eta)
&=\xi\eta,
\\
\widehat\phi_{4}(\xi,\eta)
&=(1-\xi)\eta.
\end{aligned}
```

Higher-order Lagrange elements add support points on edges and, eventually, in cell interiors. Their construction follows the same principle: the number and placement of the support points must make the interpolation problem unisolvent.

For the rest of the main chapter, first-order elements are enough to expose the mechanics without hiding them under high-order indexing.

## 5. Meshes and global conformity

Let $\mathcal T_{h}$ be a mesh of $\Omega$: a finite collection of cells $K$ whose interiors do not overlap and whose union covers the computational domain.

For a simple conforming mesh, neighboring cells meet along complete lower-dimensional entities – a common vertex, edge, or face – rather than overlapping arbitrarily. We denote the diameter of a cell by $h_{K}$ and often write

$$
h:=\max_{K\in\mathcal T_{h}}h_{K}.
$$

A local polynomial space on each cell does not yet define a global $H^{1}$ finite-element space. We must specify how neighboring cell functions are glued together.

For continuous Lagrange elements, degrees of freedom on shared mesh entities are **identified**. If two triangles share a vertex, the local point-evaluation degree of freedom at that vertex receives one global index, not two. For higher-order elements, the same happens to edge support points shared by neighboring cells.

This identification forces the two local polynomial pieces to take the same values at enough points on their common interface. Because the trace of a polynomial is again a polynomial of the appropriate degree, equality of all interface degrees of freedom forces the two traces to agree along the entire interface. The resulting global function is continuous.

For first-order triangles the mechanism is particularly transparent. The trace of a $\mathbb P_{1}$ polynomial on an edge is a one-dimensional affine function. Two such traces are equal everywhere on the edge if they agree at the two edge vertices. Sharing those vertex degrees of freedom is therefore enough to make the global function continuous across the edge.

A typical conforming scalar space can be written schematically as

```math
V_{h}
:=
\left\{
 v_{h}\in C^{0}(\overline\Omega)
 :
 v_{h}\rvert_{K}\in\mathcal P(K)
 \text{ for every }K\in\mathcal T_{h},
 \text{ with the required boundary values}
\right\}.
```

For $\mathbb P_{1}$ on a triangular mesh,

$$
\mathcal P(K)=\mathbb P_{1}(K).
$$

For a mapped $\mathbb Q_{1}$ construction on a quadrilateral mesh, the local functions are obtained from $\mathbb Q_{1}$ functions on the reference square through the cell map introduced below.

The continuity condition is what makes such piecewise-polynomial spaces conforming subspaces of $H^{1}(\Omega)$.

### A note on mesh regularity

Not every sequence of geometrically valid cells gives a numerically useful family of spaces. If triangles become arbitrarily thin or distorted as the mesh is refined, constants in approximation estimates can deteriorate. Finite-element analysis therefore imposes conditions such as **shape regularity**.

The detailed geometric definitions and the resulting interpolation estimates belong to the error-analysis theory left outside the main path. Here we assume a reasonable family of meshes and focus on the algebra and function-space structure built on them.

## 6. Continuous and discontinuous spaces encode different mathematics

Sharing degrees of freedom across cells is a design choice, not an inevitable property of polynomial approximation.

Suppose we place one constant polynomial on each cell and do **not** identify values between neighboring cells. Then

```math
U_{h}^{0}
:=
\left\{
 u_{h}
 :
 u_{h}\rvert_{K}\text{ is constant for every }K\in\mathcal T_{h}
\right\}.
```

A convenient basis contains one cell indicator per cell:

$$
\chi_{K}(x)
:=
\begin{cases}
1, & x\in K,\\
0, & x\notin K.
\end{cases}
$$

Then

$$
u_{h}
=
\sum_{K\in\mathcal T_{h}}u_{K}\chi_{K}.
$$

Each coefficient $u_{K}$ is one cell value. Different cells retain independent values, so jumps are allowed at cell interfaces.

Such a space is not generally a subspace of $H^{1}(\Omega)$, but it is a perfectly natural finite-dimensional subspace of $L^{2}(\Omega)$. That makes cellwise constants especially useful for variables whose continuous space is $L^{2}$ rather than $H^{1}$ – distributed controls are an important example.

This distinction should not be confused with discontinuous Galerkin methods. A discontinuous finite-element space can be used simply to represent an $L^{2}$ variable. Discontinuous Galerkin methods go further: they use discontinuous trial/test spaces for the PDE itself and add interface terms to recover a suitable weak formulation. That theory is outside the present chapter.

### The cellwise-constant mass matrix

For the cell-indicator basis,

$$
M_{KL}
:=
\int_{\Omega}\chi_{L}\chi_{K} \mathrm{d}x.
$$

Different indicators have disjoint interiors, so for $K\ne L$,

$$
M_{KL}=0.
$$

For $K=L$,

$$
M_{KK}
=
\int_{K}1 \mathrm{d}x
=
|K|.
$$

Therefore

$$
M
=
\mathrm{diag}
\left(
|K_{1}|,\ldots,|K_{N_{K}}|
\right).
$$

The diagonal structure is not a generic property of mass matrices. It comes from this particular discontinuous basis and its nonoverlapping supports.

## 7. Reference cells and physical cells

It would be wasteful to derive and implement a new set of basis functions independently on every physical cell. Finite-element codes instead define shape functions on a small set of **reference cells** and map them to the physical mesh.

The one-dimensional case shows the idea with almost no notation. Let

$$
\widehat K=[0,1]
$$

with reference basis

$$
\widehat\phi_{1}(\xi)=1-\xi,
\qquad
\widehat\phi_{2}(\xi)=\xi.
$$

Let a physical cell be

$$
K=[x_{a},x_{b}],
\qquad
h_{K}:=x_{b}-x_{a}.
$$

The affine map

$$
F_{K}:\widehat K\to K,
\qquad
x=F_{K}(\xi):=x_{a}+h_{K}\xi
$$

maps the reference interval onto the physical interval. Define the physical basis by composition with the inverse map:

$$
\phi_{i}^{K}(x)
:=
\widehat\phi_{i}(F_{K}^{-1}(x)).
$$

The change of variables gives

$$
\mathrm{d}x=h_{K} \mathrm{d}\xi,
$$

and the chain rule gives

$$
\frac{\mathrm{d}\phi_{i}^{K}}{\mathrm{d}x}
=
\frac{1}{h_{K}}
\frac{\mathrm{d}\widehat\phi_{i}}{\mathrm{d}\xi}.
$$

These two scaling factors – one from the volume element, one from the derivative – will immediately determine the local mass and stiffness matrices.

### Affine maps in several dimensions

For an affine simplex map, write

$$
x
=
F_{K}(\widehat x)
:=
B_{K}\widehat x+b_{K},
$$

where $B_{K}$ is an invertible $d\times d$ matrix. Then

$$
\mathrm{d}x
=
|\det B_{K}| \mathrm{d}\widehat x.
$$

If

$$
\phi^{K}(x)
=
\widehat\phi(F_{K}^{-1}(x)),
$$

then the multivariable chain rule gives

$$
\nabla_{x}\phi^{K}(x)
=
B_{K}^{-\mathsf T}
\nabla_{\widehat x}\widehat\phi(\widehat x).
$$

Hence an integral such as

$$
\int_{K}\nabla\phi_{j}^{K}\cdot\nabla\phi_{i}^{K} \mathrm{d}x
$$

can be evaluated on the reference cell using the same reference shape functions and cell-specific geometric factors.

For general quadrilaterals, hexahedra, curved cells, or higher-order geometry mappings, the Jacobian need not be constant over the cell. The same principle still holds, but the geometric factors must be evaluated at the integration points. Isoparametric mappings, curved geometry, and specialized transformations for $H(\mathrm{div})$ or $H(\mathrm{curl})$ elements are beyond the present scope.

## 8. Local stiffness, mass, and load terms

Return to the one-dimensional Poisson equation. On one physical cell $K=[x_{a},x_{b}]$, the local weak contribution is

$$
\int_{K}y_{h}'v_{h}' \mathrm{d}x.
$$

Using the two local $\mathbb P_{1}$ basis functions, define the **local stiffness matrix**

$$
A_{ij}^{K}
:=
\int_{K}
(\phi_{j}^{K})'
(\phi_{i}^{K})'
 \mathrm{d}x,
\qquad
i,j\in\{1,2\}.
$$

On the reference cell,

$$
\widehat\phi_{1}'=-1,
\qquad
\widehat\phi_{2}'=1.
$$

Using the mapping formulas,

```math
\begin{aligned}
A_{ij}^{K}
&=
\int_{0}^{1}
\left(
\frac{1}{h_{K}}\widehat\phi_{j}'
\right)
\left(
\frac{1}{h_{K}}\widehat\phi_{i}'
\right)
h_{K} \mathrm{d}\xi
\\
&=
\frac{1}{h_{K}}
\int_{0}^{1}
\widehat\phi_{j}'\widehat\phi_{i}'
 \mathrm{d}\xi.
\end{aligned}
```

Therefore

$$
A^{K}
=
\frac{1}{h_{K}}
\begin{bmatrix}
1&-1\\
-1&1
\end{bmatrix}.
$$

The corresponding **local mass matrix** represents the $L^{2}$ bilinear form:

$$
M_{ij}^{K}
:=
\int_{K}\phi_{j}^{K}\phi_{i}^{K} \mathrm{d}x.
$$

Changing variables,

$$
M_{ij}^{K}
=
h_{K}
\int_{0}^{1}
\widehat\phi_{j}\widehat\phi_{i} \mathrm{d}\xi.
$$

The three required integrals are

$$
\int_{0}^{1}(1-\xi)^{2} \mathrm{d}\xi
=
\frac{1}{3},
$$

$$
\int_{0}^{1}\xi(1-\xi) \mathrm{d}\xi
=
\frac{1}{6},
$$

$$
\int_{0}^{1}\xi^{2} \mathrm{d}\xi
=
\frac{1}{3}.
$$

Thus

$$
M^{K}
=
\frac{h_{K}}{6}
\begin{bmatrix}
2&1\\
1&2
\end{bmatrix}.
$$

For a source $f$, the **local load vector** is

$$
f_{i}^{K}
:=
\int_{K}f\phi_{i}^{K} \mathrm{d}x.
$$

Unlike the stiffness and mass matrices for these simple affine elements, this integral depends on the actual source data. It is typically evaluated numerically.

The local formulas already show why “mass matrix” and “stiffness matrix” refer to mathematical roles, not just sparse arrays. They are coordinate representations of different bilinear forms:

$$
M
\leftrightarrow
(u,v)_{L^{2}},
\qquad
A
\leftrightarrow
(\nabla u,\nabla v)_{L^{2}}.
$$

## 9. Numerical quadrature turns cell integrals into finite sums

Even after the finite-dimensional space has been chosen, matrix and vector entries are defined by integrals. A computer normally evaluates those integrals through a **quadrature rule**.

On a reference cell $\widehat K$, a quadrature rule consists of points

$$
\widehat x_{1},\ldots,\widehat x_{n_{q}}
$$

and weights

$$
\widehat w_{1},\ldots,\widehat w_{n_{q}}
$$

such that

$$
\int_{\widehat K}g(\widehat x) \mathrm{d}\widehat x
\approx
\sum_{q=1}^{n_{q}}
\widehat w_{q}g(\widehat x_{q}).
$$

A rule has a certain **degree of exactness** if the equality is exact for every polynomial up to that degree. The rule should be chosen accurately enough for the products that occur in the local integrals, together with any non-polynomial coefficients or data.

Under a cell map $F_{K}$, the physical-cell integral becomes

$$
\int_{K}g(x) \mathrm{d}x
=
\int_{\widehat K}
 g(F_{K}(\widehat x))
 |\det J_{K}(\widehat x)|
 \mathrm{d}\widehat x,
$$

where $J_{K}$ is the Jacobian of the map. Quadrature therefore gives

$$
\int_{K}g(x) \mathrm{d}x
\approx
\sum_{q=1}^{n_{q}}
\widehat w_{q}
 g(F_{K}(\widehat x_{q}))
 |\det J_{K}(\widehat x_{q})|.
$$

For an affine map, $J_{K}$ is constant. For a more general mapped cell it varies with the quadrature point.

### Example: assembling a diffusion entry

Let $\kappa:\Omega\to\mathbb R$ denote the scalar **diffusion coefficient** multiplying the gradient term in the weak form. A diffusion matrix entry on one cell has the form

$$
A_{ij}^{K}
=
\int_{K}
\kappa(x)
\nabla\phi_{j}^{K}(x)
\cdot
\nabla\phi_{i}^{K}(x)
\mathrm{d}x.
$$

A quadrature approximation is

```math
A_{ij}^{K}
\approx
\sum_{q=1}^{n_{q}}
 w_{q}^{K}
 \kappa(x_{q}^{K})
 \nabla\phi_{j}^{K}(x_{q}^{K})
 \cdot
 \nabla\phi_{i}^{K}(x_{q}^{K}),
```

where $x_{q}^{K}$ and $w_{q}^{K}$ are the mapped quadrature points and weights.

This formula exposes the data a finite-element assembly loop actually needs at each quadrature point:

```text
coefficient value κ(x_q)
basis values φ_i(x_q) when needed
basis gradients ∇φ_i(x_q) when needed
geometric Jacobian information
quadrature weight
```

For the affine one-dimensional $\mathbb P_{1}$ element derived above, this can be seen immediately. With a constant diffusion coefficient, the local stiffness integrand is constant on the reference cell, whereas each mass-matrix integrand $\widehat\phi_{i}\widehat\phi_{j}$ is quadratic. A rule exact only for constants is therefore enough for the former but not for the latter. When coefficients or data are non-polynomial, exactness is generally replaced by choosing a rule accurate enough that quadrature error is negligible at the intended discretization order.

The companion [03a · Further finite-element notes](03a-further-finite-element-notes.md) goes one step further: it distinguishes nodal interpolation, $L^{2}$ projection, cell averaging, and direct evaluation of continuous data, then relates these choices to Gaussian and tensor-product quadrature formulas. Those details matter when a source, coefficient, desired state, observation, or control is not already represented in the finite-element space.

The ideal Galerkin problem assumes the defining integrals are evaluated exactly. Numerical quadrature introduces another approximation stage. When the quadrature is not exact, the assembled bilinear form and functional differ slightly from their ideal Galerkin counterparts. Detailed analysis of this consistency error is part of the theory sometimes discussed under **variational crimes** and is outside the main path here.

## 10. Local-to-global degree-of-freedom numbering

Every cell has local degrees of freedom. The global finite-element space has one global coefficient vector. We therefore need a map from local indices to global indices.

Suppose a cell $K$ has $m$ local basis functions

$$
\phi_{1}^{K},\ldots,\phi_{m}^{K}.
$$

Let

$$
I_{K}(\alpha)
$$

be the global index associated with local degree of freedom $\alpha$ on cell $K$.

For continuous elements, neighboring cells can map different local indices to the same global index whenever those local degrees of freedom live on a shared mesh entity.

A cell-local coefficient vector is therefore obtained by gathering entries from the global vector:

$$
(\mathbf y^{K})_{\alpha}
=
\mathbf y_{I_{K}(\alpha)}.
$$

Conversely, a local matrix contribution is **scattered** into the global matrix using the same index map.

The numbering itself is not part of the mathematical finite-element space. Renumbering the global degrees of freedom merely permutes the coefficient vector and the corresponding rows and columns of assembled matrices. Good numberings can improve memory locality or solver performance, but they do not change the represented functions.

### Why the assembly formula is a sum over cells

For the Poisson bilinear form,

$$
a(y_{h},v_{h})
=
\int_{\Omega}\nabla y_{h}\cdot\nabla v_{h} \mathrm{d}x.
$$

Because the mesh partitions the domain,

$$
a(y_{h},v_{h})
=
\sum_{K\in\mathcal T_{h}}
\int_{K}\nabla y_{h}\cdot\nabla v_{h} \mathrm{d}x.
$$

Write the restrictions to one cell in its local basis:

$$
y_{h}\rvert_{K}
=
\sum_{\beta=1}^{m}
 y_{I_{K}(\beta)}\phi_{\beta}^{K},
$$

$$
v_{h}\rvert_{K}
=
\sum_{\alpha=1}^{m}
 v_{I_{K}(\alpha)}\phi_{\alpha}^{K}.
$$

Then

```math
\begin{aligned}
a(y_{h},v_{h})
&=
\sum_{K}
\sum_{\alpha=1}^{m}
\sum_{\beta=1}^{m}
 v_{I_{K}(\alpha)}
 \left(
 \int_{K}
 \nabla\phi_{\beta}^{K}
 \cdot
 \nabla\phi_{\alpha}^{K}
 \mathrm{d}x
 \right)
 y_{I_{K}(\beta)}
\\
&=
\sum_{K}
\sum_{\alpha,\beta}
 v_{I_{K}(\alpha)}
 A_{\alpha\beta}^{K}
 y_{I_{K}(\beta)}.
\end{aligned}
```

Thus the global matrix is obtained by adding

$$
A_{\alpha\beta}^{K}
$$

to the entry

$$
A_{I_{K}(\alpha),I_{K}(\beta)}.
$$

In algorithmic form:

```text
for each cell K
    compute local matrix A^K
    get local-to-global indices I_K
    for local rows α
        for local columns β
            A[I_K(α), I_K(β)] += A^K[α,β]
```

The same pattern applies to a local load vector:

$$
f_{I_{K}(\alpha)}
\mathrel{+}=
f_{\alpha}^{K}.
$$

### A complete one-dimensional assembly

Take four uniform cells on $(0,1)$, so

$$
h=\frac{1}{4}.
$$

The global nodes are

$$
x_{0},x_{1},x_{2},x_{3},x_{4},
$$

with $x_{0}$ and $x_{4}$ fixed by homogeneous Dirichlet conditions. The independent finite-element basis therefore contains the three interior hat functions associated with $x_{1},x_{2},x_{3}$.

Each local stiffness matrix is

$$
A^{K}
=
\frac{1}{h}
\begin{bmatrix}
1&-1\\
-1&1
\end{bmatrix}.
$$

If we first assemble on the unconstrained nodal basis including the two boundary nodes, the resulting matrix is

```math
A_{\mathrm{full}}
=
\frac{1}{h}
\begin{bmatrix}
 1&-1& 0& 0& 0\\
-1& 2&-1& 0& 0\\
 0&-1& 2&-1& 0\\
 0& 0&-1& 2&-1\\
 0& 0& 0&-1& 1
\end{bmatrix}.
```

Restricting the test and trial functions to the three free interior degrees of freedom gives

$$
A
=
\frac{1}{h}
\begin{bmatrix}
2&-1&0\\
-1&2&-1\\
0&-1&2
\end{bmatrix}.
$$

The tridiagonal pattern is not put in by hand. It follows from local support: two global basis functions interact only if their supports overlap on at least one cell.

In higher dimensions the pattern is less visually regular because it depends on mesh connectivity and global numbering, but the same locality makes the matrix sparse.

## 11. From the Galerkin equation to the global algebraic system

Let

$$
V_{h}
=
\mathrm{span}\{\varphi_{1},\ldots,\varphi_{n}\}
$$

and write

$$
y_{h}
=
\sum_{j=1}^{n}y_{j}\varphi_{j}.
$$

The discrete weak problem requires

$$
a(y_{h},v_{h})=F(v_{h})
\qquad
\text{for every }v_{h}\in V_{h}.
$$

Because both sides are linear in the test function, it is enough to enforce the equation for every basis function $\varphi_{i}$:

$$
a(y_{h},\varphi_{i})
=
F(\varphi_{i}),
\qquad
i=1,\ldots,n.
$$

Substituting the expansion of $y_{h}$,

$$
\sum_{j=1}^{n}
y_{j}a(\varphi_{j},\varphi_{i})
=
F(\varphi_{i}).
$$

Define

$$
A_{ij}
:=
a(\varphi_{j},\varphi_{i}),
$$

and

$$
f_{i}
:=
F(\varphi_{i}).
$$

Then

$$
A\mathbf y=\mathbf f.
$$

For Poisson,

$$
A_{ij}
=
\int_{\Omega}
\nabla\varphi_{j}\cdot\nabla\varphi_{i}
 \mathrm{d}x.
$$

This is the global **stiffness matrix**.

Several facts follow immediately from the variational problem.

If the bilinear form is symmetric,

$$
a(v,w)=a(w,v),
$$

then

$$
A^{\mathsf T}=A.
$$

If the bilinear form is coercive on $V_{h}$, then for every nonzero coefficient vector $\mathbf z$ representing

$$
z_{h}
=
\sum_{j}z_{j}\varphi_{j},
$$

we have

$$
\mathbf z^{\mathsf T}A\mathbf z
=
a(z_{h},z_{h})
>0.
$$

Thus the matrix is positive definite. These algebraic properties are inherited from the variational form; they are not accidental features of assembly.

## 12. Essential boundary conditions and independent coordinates

Homogeneous Dirichlet conditions were easy above because we simply omitted the boundary hat functions from $V_{h}\subset H_{0}^{1}(\Omega)$. Nonhomogeneous Dirichlet data and more general linear constraints make the coordinate structure more visible.

Suppose the full finite-element coefficient vector has dimension $N$ and is denoted

$$
\mathbf y_{\mathrm{phys}}
\in
\mathbb R^{N}.
$$

Some coefficients are fixed or linearly dependent, leaving only $n$ independent coordinates

$$
\mathbf z\in\mathbb R^{n}.
$$

As in [02 · Duality, derivatives, and adjoints](02-duality-derivatives-and-adjoints.md), write the affine reconstruction as

$$
\mathbf y_{\mathrm{phys}}
=
P\mathbf z+\boldsymbol\ell.
$$

Here $P\in\mathbb R^{N\times n}$ embeds homogeneous independent variations into the full coefficient vector, while $\boldsymbol\ell$ carries the fixed affine part.

A perturbation satisfies

$$
\delta\mathbf y_{\mathrm{phys}}
=
P\delta\mathbf z.
$$

### Fixed nodal Dirichlet values

For simple nodal Lagrange elements, $P$ can be an insertion matrix. Suppose a five-entry physical vector has its second and fifth entries prescribed:

$$
\mathbf y_{\mathrm{phys}}
=
\begin{bmatrix}
z_{1}\\
g_{1}\\
z_{2}\\
z_{3}\\
g_{4}
\end{bmatrix}.
$$

Then

```math
P
=
\begin{bmatrix}
1&0&0\\
0&0&0\\
0&1&0\\
0&0&1\\
0&0&0
\end{bmatrix},
\qquad
\boldsymbol\ell
=
\begin{bmatrix}
0\\
g_{1}\\
0\\
0\\
g_{4}
\end{bmatrix}.
```

### The reduced linear system

Suppose the unconstrained weak residual, expressed in the full nodal coordinates, is

$$
A\mathbf y_{\mathrm{phys}}-\mathbf f.
$$

Only variations in the independent subspace are admissible. Testing the residual against all independent variations $P\delta\mathbf z$ means

$$
(P\delta\mathbf z)^{\mathsf T}
\left(
A(P\mathbf z+\boldsymbol\ell)-\mathbf f
\right)
=0
$$

for every $\delta\mathbf z$. Therefore

$$
P^{\mathsf T}
\left(
A(P\mathbf z+\boldsymbol\ell)-\mathbf f
\right)
=0,
$$

or

$$
P^{\mathsf T}AP\mathbf z
=
P^{\mathsf T}(\mathbf f-A\boldsymbol\ell).
$$

This formula contains both familiar boundary-condition operations:

- $\boldsymbol\ell$ shifts the right-hand side because the fixed physical values contribute to the equation;
- $P^{\mathsf T}$ restricts the residual to the independent test directions.

For simple elimination of boundary nodes, $P^{\mathsf T}AP$ is just the free/free submatrix. The formula also handles more general affine relations.

### Hanging-node-style constraints

Local mesh refinement can create a point on one cell edge that lies in the middle of a coarser neighboring edge. For a continuous piecewise-linear field, the value at such a constrained point may have to satisfy a relation such as

$$
y_{c}
=
\frac{1}{2}y_{1}
+
\frac{1}{2}y_{2}.
$$

Then the corresponding row of $P$ contains interpolation weights rather than a single $1$. The constrained physical coefficient is reconstructed from neighboring independent coordinates.

The geometric details of adaptive meshes are outside this chapter, but the algebraic lesson is general: **the physical finite-element coefficient vector need not be the vector of independent unknowns**.

## 13. Mass, stiffness, and coupling matrices come from spaces and bilinear forms

Matrix names are useful only if their underlying spaces and forms remain visible.

### Mass matrix

For one discrete space

$$
V_{h}
=
\mathrm{span}\{\varphi_{1},\ldots,\varphi_{n}\},
$$

the $L^{2}$ mass matrix is

$$
M_{ij}
:=
\int_{\Omega}\varphi_{j}\varphi_{i} \mathrm{d}x.
$$

It represents the $L^{2}$ pairing of two finite-element functions:

$$
(v_{h},w_{h})_{L^{2}}
=
\mathbf v^{\mathsf T}M\mathbf w.
$$

### Stiffness matrix

For Poisson,

$$
A_{ij}
:=
\int_{\Omega}
\nabla\varphi_{j}\cdot\nabla\varphi_{i}
 \mathrm{d}x.
$$

It represents the energy bilinear form

$$
a(v_{h},w_{h})
=
\mathbf v^{\mathsf T}A\mathbf w.
$$

### Coupling between different spaces

The same assembly mechanism also explains what happens when a PDE contains another field rather than only fixed data. Consider the Poisson equation with a distributed control,

$$
-\Delta y=f+u
\qquad
\text{in }\Omega.
$$

We are not yet formulating an optimization problem; for the present purpose, $u$ is simply another field that enters the PDE. Its weak form contains

$$
\int_{\Omega}\nabla y\cdot\nabla v\mathrm{d}x
=
\int_{\Omega}fv\mathrm{d}x
+
\int_{\Omega}uv\mathrm{d}x.
$$

Let the state/test space be

$$
V_{h}
=
\mathrm{span}\{\varphi_{1},\ldots,\varphi_{n_{y}}\},
$$

with

$$
y_{h}
=
\sum_{j=1}^{n_{y}}y_{j}\varphi_{j},
$$

and let the control use a possibly different discrete space

$$
U_{h}
=
\mathrm{span}\{\psi_{1},\ldots,\psi_{n_{u}}\},
$$

with

$$
u_{h}
=
\sum_{k=1}^{n_{u}}u_{k}\psi_{k}.
$$

Testing the weak equation with $v_{h}=\varphi_{i}$ gives

```math
\begin{aligned}
\sum_{j=1}^{n_{y}}
y_{j}
\int_{\Omega}
\nabla\varphi_{j}\cdot\nabla\varphi_{i}
\mathrm{d}x
&=
\int_{\Omega}f\varphi_{i}\mathrm{d}x
+
\sum_{k=1}^{n_{u}}
u_{k}
\int_{\Omega}
\psi_{k}\varphi_{i}
\mathrm{d}x.
\end{aligned}
```

Define

$$
A_{ij}
:=
\int_{\Omega}
\nabla\varphi_{j}\cdot\nabla\varphi_{i}
\mathrm{d}x,
$$

$$
f_{i}
:=
\int_{\Omega}f\varphi_{i}\mathrm{d}x,
$$

and

$$
B_{ik}
:=
\int_{\Omega}
\psi_{k}\varphi_{i}
\mathrm{d}x.
$$

The discrete equation is therefore

$$
A\mathbf y
=
\mathbf f+B\mathbf u.
$$

The difference between $\mathbf f$ and $B\mathbf u$ is structural. The forcing $f$ is fixed data, so after testing it against the basis functions its contribution is a known vector. The control $u_{h}$ is represented by variable coefficients $\mathbf u$, so testing it produces a linear map from control coefficients to state-test covectors. The matrix $B$ is the coordinate representation of that map.

If the forcing itself depended on an unknown finite-dimensional coefficient vector, the same reasoning would produce another coupling matrix for it. The distinction is therefore not inherent in the words *forcing* and *control*. It comes from which quantities are fixed data and which quantities carry variable coordinates.

Since

$$
B
\in
\mathbb R^{n_{y}\times n_{u}},
$$

there is no reason for $n_{y}$ and $n_{u}$ to agree. If $U_{h}=V_{h}$ with the same basis, $B$ numerically coincides with the mass matrix. If $U_{h}$ is cellwise constant while $V_{h}$ is continuous piecewise linear, $B$ is generally rectangular. Its role remains **control-to-test coupling** in both cases.

### A one-dimensional $\mathbb P_{1}$–$\mathbb P_{0}$ coupling

Let the state/test space use continuous piecewise-linear basis functions and let the control have one constant value $u_{e}$ on each cell $K_{e}$.

On one cell, the two local state basis functions satisfy

$$
\int_{K_{e}}\phi_{1}^{K_{e}} \mathrm{d}x
=
\int_{K_{e}}\phi_{2}^{K_{e}} \mathrm{d}x
=
\frac{h_{e}}{2}.
$$

The cell control therefore contributes

$$
\frac{h_{e}}{2}u_{e}
$$

to each of the two local test equations associated with that cell.

A single control coefficient can influence several state equations, and a single state test function can overlap several control cells. The resulting matrix is sparse because the coupling is still local, but it need not be square.

This example is the finite-element origin of many rectangular control-to-state operators in PDE-constrained optimization.

## 14. Coefficient vectors are not physical fields

For a chosen basis,

$$
y_{h}(x)
=
\sum_{j=1}^{n}y_{j}\varphi_{j}(x).
$$

The vector

$$
\mathbf y
=
(y_{1},\ldots,y_{n})^{\mathsf T}
$$

is a set of coordinates. The physical field is the function obtained after combining those coordinates with the basis.

This distinction has several practical consequences.

### Evaluating the field is a basis operation

At a physical point $x$,

$$
y_{h}(x)
=
\boldsymbol\varphi(x)^{\mathsf T}\mathbf y,
$$

where

$$
\boldsymbol\varphi(x)
:=
\begin{bmatrix}
\varphi_{1}(x)\\
\vdots\\
\varphi_{n}(x)
\end{bmatrix}.
$$

Only basis functions whose support contains $x$ contribute, so evaluation is naturally local.

### A coefficient is not always a nodal value

For the nodal Lagrange bases used above,

$$
y_{j}=y_{h}(x_{j}).
$$

That is a property of this basis, not a universal definition of a finite-element coefficient. In other finite elements, degrees of freedom may be moments or fluxes rather than point values.

### Data need not have finite-element coefficients

A source $f(x)$ or desired field $y_{d}(x)$ may be known analytically or through some external representation. To assemble

$$
f_{i}
=
\int_{\Omega}f\varphi_{i} \mathrm{d}x,
$$

we only need to evaluate $f$ where the quadrature requires it. There is no mathematical requirement to first interpolate $f$ into $V_{h}$.

Likewise, a tracking term

$$
\frac{1}{2}
\int_{\Omega}
(y_{h}-y_{d})^{2}
 \mathrm{d}x
$$

can be evaluated directly at quadrature points even when $y_{d}$ has no state-space coefficient vector.

When a discrete representation of external data is actually needed, several different constructions are possible. Nodal interpolation matches point values, an $L^{2}$ projection matches integral moments, and a cellwise-constant representation may use cell averages. These are different mathematical operations and should not be treated as interchangeable ways of “putting a function on the mesh.” The companion [03a · Further finite-element notes](03a-further-finite-element-notes.md) derives these alternatives.

### The same field has different coordinates in different bases

Changing the finite-element basis changes $\mathbf y$ while leaving $y_{h}$ unchanged. This is the finite-element instance of the coordinate distinction developed in [02 · Duality, derivatives, and adjoints](02-duality-derivatives-and-adjoints.md).

Raw Euclidean operations on coefficient arrays therefore acquire mathematical meaning only after the basis, space, and intended pairing or metric are known.

## 15. A compact assembly workflow

The mathematical steps above can be summarized without committing to any particular finite-element library.

For a scalar Poisson problem:

```text
mesh
    ↓
choose local finite element (K, P, Σ)
    ↓
identify global DoFs and constraints
    ↓
for each cell
    map reference quadrature to the physical cell
    evaluate basis values / gradients
    accumulate local matrix and load vector
    scatter local contributions to global indices
    ↓
apply / eliminate affine constraints
    ↓
solve the resulting sparse linear system
    ↓
reconstruct the physical finite-element field
```

Each line corresponds to a mathematical operation introduced in this chapter. A finite-element library packages and automates much of this bookkeeping, but the library objects remain easier to understand when this underlying flow is clear.

[09 · The deal.II finite-element workflow](09-dealii-finite-element-workflow.md) will map this workflow onto the common deal.II objects used to implement it.

## 16. What convergence means: Galerkin orthogonality and Céa's lemma

The finite-element construction gives a sequence of spaces $V_{h}$. Why should their solutions approach the true weak solution as the mesh is refined?

For the coercive setting of [01 · Function spaces and weak PDEs](01-function-spaces-and-weak-pdes.md), the first part of the answer is almost entirely algebraic.

Assume the bilinear form is bounded:

$$
|a(v,w)|
\le
M
\lVert v\rVert_{V}
\lVert w\rVert_{V},
$$

and coercive:

$$
a(v,v)
\ge
\alpha
\lVert v\rVert_{V}^{2},
$$

for constants $M>0$ and $\alpha>0$.

Let $y\in V$ solve the continuous problem and $y_{h}\in V_{h}\subset V$ solve the Galerkin problem. We already derived

$$
a(y-y_{h},v_{h})=0
\qquad
\text{for every }v_{h}\in V_{h}.
$$

Take any $w_{h}\in V_{h}$. Since $w_{h}-y_{h}\in V_{h}$,

$$
a(y-y_{h},w_{h}-y_{h})=0.
$$

Hence

```math
\begin{aligned}
\alpha
\lVert y-y_{h}\rVert_{V}^{2}
&\le
 a(y-y_{h},y-y_{h})
\\
&=
 a(y-y_{h},y-w_{h})
+
 a(y-y_{h},w_{h}-y_{h})
\\
&=
 a(y-y_{h},y-w_{h})
\\
&\le
 M
 \lVert y-y_{h}\rVert_{V}
 \lVert y-w_{h}\rVert_{V}.
\end{aligned}
```

If $y\ne y_{h}$, cancel one factor of $\lVert y-y_{h}\rVert_{V}$ to obtain

$$
\lVert y-y_{h}\rVert_{V}
\le
\frac{M}{\alpha}
\lVert y-w_{h}\rVert_{V}.
$$

Because this holds for every $w_{h}\in V_{h}$,

$$
\lVert y-y_{h}\rVert_{V}
\le
\frac{M}{\alpha}
\inf_{w_{h}\in V_{h}}
\lVert y-w_{h}\rVert_{V}.
$$

This is the basic form of **Céa's lemma**.

It says that the Galerkin solution is, up to the stability factor $M/\alpha$, as good as the best approximation to $y$ available inside the chosen discrete space.

Céa's lemma does **not** by itself provide a power of $h$. To conclude

$$
\lVert y-y_{h}\rVert_{V}\to0
\qquad
\text{as }h\to0,
$$

we still need the spaces to approximate every relevant $V$ function increasingly well:

$$
\inf_{w_{h}\in V_{h}}
\lVert y-w_{h}\rVert_{V}
\longrightarrow0.
$$

Quantifying that best-approximation error requires interpolation theory, assumptions on mesh regularity, and regularity of the exact solution. Those ingredients lead to estimates such as powers of $h$, but deriving them is deliberately outside the main background needed by the `nmopt` manual.

There is one further caveat. The argument above assumes the discrete bilinear form and right-hand side are exact restrictions of the continuous ones. Numerical quadrature, nonconforming spaces, approximate geometry, and other modifications can break exact Galerkin orthogonality. Their analysis requires additional consistency estimates.

## 17. Scope frontier

The chapter has developed the finite-element mechanics needed by the rest of this background collection:

- conforming Galerkin restriction;
- local finite-element triplets and unisolvence;
- Lagrange $\mathbb P_{1}$ and $\mathbb Q_{1}$ elements;
- discontinuous cellwise-constant spaces;
- meshes and shared degrees of freedom;
- reference-to-physical mappings;
- quadrature;
- local matrices and local-to-global assembly;
- affine coordinate constraints;
- mass, stiffness, and cross-space coupling matrices;
- the distinction between coefficient vectors and finite-element fields;
- the basic best-approximation meaning of Galerkin convergence.

Several important subjects are intentionally deferred. A guided first pass through the main deferred topics is collected in [03a · Further finite-element notes](03a-further-finite-element-notes.md), so the main chapter can keep its conforming scalar path intact.

**Detailed approximation estimates.** Bramble–Hilbert theory, interpolation constants, duality arguments for $L^{2}$ errors, and precise convergence rates are useful numerical-PDE theory but are not required to understand the project manual.

**Adaptive methods and a posteriori error estimation.** Local error indicators, refinement strategies, hanging-node generation, and convergence of adaptive loops deserve their own treatment. We only used a hanging-node-style relation to motivate affine constraints.

**Discontinuous Galerkin methods.** A discontinuous space for an $L^{2}$ control is simple. Using discontinuous spaces for the PDE solution requires interface fluxes, penalties, and a different consistency/stability analysis.

**Mixed and saddle-point finite elements.** Stokes, Darcy, and mixed formulations introduce several fields, spaces with different conformity requirements, and discrete inf-sup conditions. [04 · Numerical linear algebra for PDE and optimization systems](04-numerical-linear-algebra.md) discusses saddle-point matrices algebraically, but not the finite-element stability theory that produces them.

**Nonconforming and Petrov–Galerkin methods.** The trial space need not be a subspace of the continuous space, and the test space need not equal the trial space. Those generalizations are important, but the conforming Galerkin setting is the right baseline for the current prerequisites.

**Advanced geometry mappings.** Curved cells, isoparametric geometry, Piola transforms, and orientation-sensitive finite elements become important for more general PDEs. [09 · The deal.II finite-element workflow](09-dealii-finite-element-workflow.md) only needs the simpler mapping picture developed here to explain the deal.II workflow.

## Where to go next

- **Default continuation:** [04 · Numerical linear algebra for PDE and optimization systems](04-numerical-linear-algebra.md) starts from the sparse operators produced here and develops conditioning, Krylov methods, preconditioning, Schur complements, and saddle-point systems.
- **For deeper finite-element theory:** [03a · Further finite-element notes](03a-further-finite-element-notes.md) develops approximation estimates, quadrature consistency, adaptivity, discontinuous methods, and mixed/inf-sup ideas.
- **For the PDE-control route:** [07 · PDE-constrained optimization](07-pde-constrained-optimization.md) uses distinct state and control spaces, rectangular control couplings, and finite-element discretization of the state–adjoint system.
- **For implementation:** [09 · The deal.II finite-element workflow](09-dealii-finite-element-workflow.md) maps meshes, finite elements, DoFs, quadrature, constraints, assembly, and solvers onto deal.II objects.
- **Into the `nmopt` manual:** the most relevant project-specific continuations are:
  - [manual 01 · Anatomy of a discrete PDE-constrained problem](../concepts/01-discrete-problem-anatomy.md), for the weak-to-discrete OCP pipeline;
  - [manual 02 · Spaces, coordinates, and duality](../concepts/02-spaces-coordinates-and-duality.md), for coefficient spaces and constrained coordinates;
  - [manual 04 · Metrics, gradients, and constraints](../concepts/04-metrics-gradients-and-constraints.md), for mass/stiffness-based optimization geometry;
  - [manual 11 · Compilation and lowering](../concepts/11-compilation-and-lowering.md), for compiler-owned FE realization;
  - [manual 13 · Integrating an existing PDE application](../concepts/13-integrating-an-existing-pde-application.md), for application-owned FE realization.

## References and further reading

Publication details, stable links, and access notes are collected in the [Background reference catalogue](references.md).

- **Stefano Mancini, Andrea Benvenuti, and Luca Heltai, *Numerical Methods for Partial Differential Equations*.**
  - §1.2.5, **Céa's lemma**, and §1.2.6, **The one-dimensional case**, for the basic Galerkin error relation and the hat-function example.
  - §1.3, especially §§1.3.1–1.3.5, for finite-dimensional spaces, bases, dual/nodal functionals, degrees of freedom, and the local finite-element triplet.
  - §1.4, especially §§1.4.1–1.4.5, for Lagrange elements, polynomial families, support points, conformity, triangulations, and DoF enumeration.
  - §2.1 for the affine reference-to-physical-cell scaling argument.
  - §4.1 for Strang's lemmas and consistency errors such as numerical quadrature.
- **Manzoni, Quarteroni, and Salsa, *Optimal Control of Partial Differential Equations*.**
  - Appendix B, §B.2.1 for Galerkin approximation, Galerkin orthogonality, and Céa's inequality.
  - Appendix B, §B.2.2 for the algebraic form of a Galerkin problem.
  - Appendix B, §B.4 for finite-element spaces, Lagrange interpolation, triangulations, and the qualitative origin of sparsity.
  - Appendix B, §B.5 for the detailed a priori error theory omitted here.
  - Chapter 6, §6.2 for distinct state/control discretization spaces and for quadrature as an additional discretization stage.
- **Dario A. Bini, Paola Boito, and Beatrice Meini, *Appunti di Istituzioni di Analisi Numerica*.**
  - Chapter 2, especially §§2.1–2.5, for interpolatory quadrature, convergence, Newton–Cotes, Fejér/Clenshaw–Curtis, and Gaussian formulas.
  - This is useful supplementary reading when the quadrature rule itself needs more detail; the notes are course-local rather than a public prerequisite.
- **Alexandre Ern and Jean-Luc Guermond, *Theory and Practice of Finite Elements*.**
  - **Finite Element Interpolation** and **Quadratures, Assembling, and Storage** are particularly close to this chapter.
- **Susanne C. Brenner and L. Ridgway Scott, *The Mathematical Theory of Finite Element Methods*.**
  - **The Construction of a Finite Element Space**, **Polynomial Approximation Theory in Sobolev Spaces**, and **Variational Crimes** provide deeper theory.
- **Philippe G. Ciarlet, *The Finite Element Method for Elliptic Problems*.**
  - Classical reference for conforming elliptic finite elements.
- **Daniele Boffi, Franco Brezzi, and Michel Fortin, *Mixed Finite Element Methods and Applications*.**
  - A deeper continuation when the mixed and inf-sup material deferred to [03a · Further finite-element notes](03a-further-finite-element-notes.md) becomes central.
