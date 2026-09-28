# 03a · Further finite-element notes

**Background navigation:** [Index](README.md) \
Main chapter: [03 · Finite elements](03-finite-elements.md) \
Previous: [03 · Finite elements](03-finite-elements.md) \
Next: [04 · Numerical linear algebra for PDE and optimization systems](04-numerical-linear-algebra.md)

## Purpose

[03 · Finite elements](03-finite-elements.md) develops the conforming Galerkin mechanics needed by the main background track: local polynomial spaces, degrees of freedom, assembly, quadrature, boundary constraints, and the basic best-approximation interpretation of the method.

This companion note looks just beyond that boundary. Its goal is not to turn the prerequisite material into a full finite-element course. Instead, it gives enough structure to understand what lies behind statements such as “the error is order $h$,” “the mesh must be shape regular,” “interpolate or project this data,” “this quadrature rule is exact to degree $r$,” “quadrature introduces a consistency error,” “refine where the residual is large,” or “this mixed pair must satisfy an inf-sup condition.”

The main chapter remains sufficient for the rest of the core background sequence. This file is for readers who want to understand the deferred theory rather than treating it as a black box.

## Before you start

You should already be comfortable with:

- weak formulations and $H^{1}$ spaces from [01 · Function spaces and weak PDEs](01-function-spaces-and-weak-pdes.md), with higher-order Sobolev notation and seminorms in [01a · Further notes on function spaces and weak PDEs](01a-further-function-space-notes.md);
- duality and bounded operators from [02 · Duality, derivatives, and adjoints](02-duality-derivatives-and-adjoints.md);
- the finite-element construction, quadrature, local-to-global assembly, and Céa estimate from [03 · Finite elements](03-finite-elements.md).

We continue to use scalar second-order elliptic problems as the main reference point. More specialized finite elements for vector-valued spaces, detailed convergence proofs for adaptive algorithms, and full discontinuous-Galerkin or mixed-method theory remain outside the scope of this note.

## What you will be able to do

After this note, you should be able to:

- explain how interpolation estimates turn Céa's lemma into a convergence rate;
- distinguish the mesh size $h_{K}$ of a cell from a measure of its shape quality;
- explain why increasingly thin or distorted cells can spoil mesh-independent estimates;
- distinguish $h$-refinement from increasing the polynomial degree;
- distinguish nodal interpolation, $L^{2}$ projection, cell averaging, and direct evaluation of external data;
- explain how quadrature nodes, weights, and polynomial exactness enter finite-element integration, including the role of Gaussian and tensor-product rules;
- explain why numerical quadrature changes the discrete variational problem and how a consistency term enters the error estimate;
- interpret a residual-based a posteriori estimator and the solve–estimate–mark–refine loop;
- explain what discontinuous Galerkin methods must add when the discrete solution space is not contained in $H^{1}$;
- recognize Petrov–Galerkin and mixed formulations and state the role of a discrete inf-sup condition.

## Roadmap

We begin with the missing link in the main chapter's convergence discussion: how well smooth functions can be approximated by finite-element spaces. That leads naturally to reference-cell scaling and mesh regularity. We then separate refinement in mesh size from refinement in polynomial degree. Before discussing quadrature error, we make a practical distinction that is easy to blur in code: interpolating a function, projecting it into a discrete space, averaging it over cells, and merely evaluating it at quadrature points are different operations. This leads to a short treatment of quadrature formulas and polynomial exactness. We then examine consistency errors and adaptive mesh refinement. Finally, we look at two broader classes of methods – discontinuous/nonconforming methods and mixed/Petrov–Galerkin methods – to see which parts of the simple conforming theory survive and which must be replaced.

## 1. From best approximation to an actual convergence rate

The main chapter ended with Céa's estimate. For a bounded and coercive bilinear form $a$ on $V$, the exact solution $u\in V$, and the conforming Galerkin solution $u_{h}\in V_{h}\subset V$,

$$
\lVert u-u_{h}\rVert_{V}
\le
\frac{M}{\alpha}
\inf_{v_{h}\in V_{h}}
\lVert u-v_{h}\rVert_{V}.
$$

This result says that Galerkin is almost as good as the best function available in $V_{h}$. It does **not** yet say how the error depends on the mesh size. To obtain a rate, we need to construct a particular $v_{h}$ whose approximation error we can estimate.

### 1.1 Interpolation provides a concrete competitor

For a Lagrange finite-element space, suppose first that $u$ is regular enough that its nodal values are meaningful. The interpolation operator

$$
I_{h}u\in V_{h}
$$

is defined by requiring the finite-element degrees of freedom of $I_{h}u$ to equal those of $u$. For a nodal basis,

$$
I_{h}u
=
\sum_{i=1}^{n}u(x_{i})\varphi_{i}.
$$

Since $I_{h}u$ is one admissible member of $V_{h}$,

$$
\inf_{v_{h}\in V_{h}}
\lVert u-v_{h}\rVert_{V}
\le
\lVert u-I_{h}u\rVert_{V}.
$$

Therefore Céa immediately gives

$$
\lVert u-u_{h}\rVert_{V}
\le
\frac{M}{\alpha}
\lVert u-I_{h}u\rVert_{V}.
$$

The finite-element error problem has now been reduced to an interpolation problem.

### 1.2 What the familiar powers of $h$ mean

Consider continuous piecewise-affine elements. On a cell $K$ of diameter $h_{K}$, a sufficiently smooth function is well approximated by its affine interpolant. The basic scaling is

$$
\lVert u-I_{h}u\rVert_{L^{2}(K)}
\lesssim
h_{K}^{2}
\lvert u\rvert_{H^{2}(K)},
$$

while differentiating loses one power of length,

$$
\lVert \nabla(u-I_{h}u)\rVert_{L^{2}(K)}
\lesssim
h_{K}
\lvert u\rvert_{H^{2}(K)}.
$$

The symbol $\lesssim$ means that the left-hand side is bounded by a constant times the right-hand side, with the important intention that the constant not deteriorate as the mesh is refined within an admissible family.

Summing the $H^{1}$ contribution over the mesh gives, schematically,

$$
\lVert u-I_{h}u\rVert_{H^{1}(\Omega)}
\lesssim
h
\lvert u\rvert_{H^{2}(\Omega)},
$$

where

$$
h:=\max_{K\in\mathcal T_{h}}h_{K}.
$$

Combining this with Céa gives the familiar first-order energy-norm estimate for piecewise-affine conforming elements,

$$
\lVert u-u_{h}\rVert_{H^{1}(\Omega)}
\lesssim
h
\lvert u\rvert_{H^{2}(\Omega)},
$$

under the assumptions needed for the interpolation estimate and for the PDE solution to belong to $H^{2}(\Omega)$.

The same idea extends to degree-$k$ elements. If the solution is sufficiently regular, one expects an $H^{1}$ interpolation error with $k$ powers of $h$,

$$
\lVert u-I_{h}u\rVert_{H^{1}(\Omega)}
\lesssim
h^{k}
\lvert u\rvert_{H^{k+1}(\Omega)}.
$$

The precise theorem requires assumptions on the mesh family and the local polynomial construction. The important mechanism is already visible: **Céa transfers approximation quality of the space into solution quality of the Galerkin method.**

### 1.3 Regularity is a ceiling on the rate

The expression

$$
\lvert u\rvert_{H^{k+1}(\Omega)}
$$

is meaningful only if the exact solution actually has that regularity. Increasing the polynomial degree cannot manufacture derivatives that the solution does not possess.

This is why convergence rates depend on more than the formal polynomial degree. They also depend on the PDE coefficients, the forcing, the boundary conditions, and the geometry of the domain, because those ingredients determine how regular the solution can be.

For example, a geometric corner or nonsmooth datum may limit the Sobolev regularity of the solution. In that case, uniform higher-order elements eventually stop delivering the rate one would predict from polynomial degree alone.

## 2. Why mesh shape enters the estimates

The main chapter mapped basis functions from a reference cell $\widehat K$ to a physical cell $K$. Error analysis uses the same map to compare norms on cells of different sizes and shapes.

### 2.1 Two geometric scales

For a cell $K$, let

$$
h_{K}
$$

be its diameter, the largest distance between two points of the cell. We also need a measure of how large a ball can fit inside the cell. Denote by

$$
\rho_{K}
$$

a characteristic inscribed scale, proportional to the diameter of a largest ball contained in $K$.

The ratio

$$
\frac{h_{K}}{\rho_{K}}
$$

measures distortion. A reasonably isotropic cell has $h_{K}$ and $\rho_{K}$ of comparable size. A very thin triangle can have a moderate diameter but a tiny $\rho_{K}$, producing a large ratio.

A family of meshes is called **shape regular** when there is a constant $C$ independent of refinement such that

$$
\frac{h_{K}}{\rho_{K}}
\le C
\qquad
\text{for every cell }K
$$

in every mesh of the family.

### 2.2 Where this ratio comes from

Suppose an affine map sends a reference cell to a physical cell,

$$
F_{K}(\widehat x)
=
B_{K}\widehat x+b_{K}.
$$

The matrix $B_{K}$ stretches reference-cell directions into physical-cell directions. Very roughly,

$$
\lVert B_{K}\rVert
\sim h_{K},
$$

while

$$
\lVert B_{K}^{-1}\rVert
\sim \rho_{K}^{-1}.
$$

The first relation says that the forward map scales lengths by the size of the cell. The second says that mapping back from a very thin cell can amplify some directions strongly.

This matters for derivatives. If

$$
\widehat v=v\circ F_{K},
$$

then the gradient transformation contains $B_{K}^{-1}$:

$$
\nabla v(x)
=
B_{K}^{-\mathsf T}
\widehat\nabla\widehat v(\widehat x).
$$

Thus bounds for derivatives on physical cells inherit factors involving $\lVert B_{K}^{-1}\rVert$. If the cells become increasingly flat or needle-like as the mesh is refined, those constants can grow even while $h$ decreases.

Shape regularity prevents that degeneration. It is not a cosmetic mesh-quality rule; it is one of the hypotheses that makes statements such as

$$
\lVert u-I_{h}u\rVert_{H^{1}}
\lesssim h^{k}
$$

meaningful with a constant that stays controlled across the mesh family.

### 2.3 Shape regular is not the same as quasi-uniform

A shape-regular mesh may contain both small and large cells, provided none of them is arbitrarily distorted. This is essential for local refinement.

A **quasi-uniform** family imposes a stronger requirement: the cell sizes are all comparable to the global mesh scale $h$. Quasi-uniformity simplifies many global estimates, but adaptive finite-element methods deliberately violate it by concentrating small cells only where they are needed.

## 3. $h$-refinement, $p$-refinement, and the regularity trade-off

There are two basic ways to enrich a polynomial finite-element space.

### 3.1 $h$-refinement

In **$h$-refinement**, the polynomial degree stays fixed while cells are subdivided. The local approximation spaces do not become individually richer; instead, the method gains more cells and therefore more local pieces.

For piecewise-affine elements, replacing a cell of diameter $h_{K}$ by smaller cells reduces the local interpolation error because the exact function is being approximated over shorter distances.

### 3.2 $p$-refinement

In **$p$-refinement**, the mesh is kept fixed while the local polynomial degree is increased,

$$
\mathbb P_{1}
\to
\mathbb P_{2}
\to
\mathbb P_{3}
\to\cdots.
$$

This can be very effective for smooth solutions, because higher-degree polynomials can reproduce increasingly rich local behavior. It also increases the number of local degrees of freedom and usually makes element matrices denser at the local level.

The gain is limited by solution regularity. If $u$ has only modest Sobolev regularity, raising $p$ beyond that regularity does not unlock the ideal high-order estimate.

### 3.3 Why local refinement is attractive

A solution may be smooth over most of the domain but have a localized layer, corner singularity, or sharp source feature. Uniform refinement spends degrees of freedom everywhere. Local $h$-refinement instead aims to reduce $h_{K}$ only where the current approximation appears inaccurate.

This motivates a posteriori error estimation, discussed below.

## 4. Discretizing functions and evaluating integrals

Finite-element code repeatedly encounters functions that are not themselves unknown finite-element fields: source terms, diffusion coefficients, desired states, observations, initial guesses, boundary data, and sometimes controls. Saying that such a function is “discretized” is ambiguous. Several different operations are possible, and they serve different purposes.

It is useful to separate four common cases:

```text
continuous or external data
        │
        ├──► nodal interpolation
        │       match selected degrees of freedom
        │
        ├──► projection
        │       match integral moments
        │
        ├──► cell averaging
        │       one representative value per cell
        │
        └──► direct evaluation
                keep the original representation
                and evaluate it where assembly needs it
```

The last option is often the simplest. If an assembly formula only needs the value $f(x_{q})$ at quadrature points, there is no reason to create a finite-element coefficient vector for $f$ first.

### 4.1 Nodal interpolation matches degrees of freedom

Let

$$
V_{h}
=
\mathrm{span}
\{\varphi_{1},\ldots,\varphi_{n}\}
$$

be a nodal Lagrange space with nodes $x_{i}$. If a function $f$ is regular enough for the point values $f(x_{i})$ to be meaningful, its nodal interpolant is

$$
I_{h}f
=
\sum_{i=1}^{n}
f(x_{i})\varphi_{i}.
$$

By construction,

$$
(I_{h}f)(x_{i})
=
f(x_{i}).
$$

More generally, if the finite element has degrees of freedom $\sigma_{i}$ rather than point evaluations, interpolation is defined by

$$
\sigma_{i}(I_{h}f)
=
\sigma_{i}(f).
$$

This definition already reveals a limitation. A generic function in $L^{2}(\Omega)$ is an equivalence class defined only almost everywhere, so point values need not be meaningful. Nodal interpolation is therefore not a universal way to discretize arbitrary data.

### 4.2 An $L^{2}$ projection matches integral moments

Suppose instead that

$$
f\in L^{2}(\Omega).
$$

The $L^{2}$ projection $P_{h}f\in V_{h}$ is defined by

$$
\int_{\Omega}
(P_{h}f-f)v_{h}
 \mathrm{d}x
=
0
\qquad
\text{for every }v_{h}\in V_{h}.
$$

Write

$$
P_{h}f
=
\sum_{j=1}^{n}c_{j}\varphi_{j}.
$$

Testing with $\varphi_{i}$ gives

$$
\sum_{j=1}^{n}
\left(
\int_{\Omega}
\varphi_{j}\varphi_{i}
 \mathrm{d}x
\right)c_{j}
=
\int_{\Omega}
f\varphi_{i}
 \mathrm{d}x.
$$

Thus the coefficient vector $\mathbf c$ satisfies

$$
M\mathbf c=\mathbf b,
$$

where $M$ is the mass matrix and

$$
b_{i}
=
\int_{\Omega}f\varphi_{i} \mathrm{d}x.
$$

Interpolation therefore determines coefficients from degrees of freedom, whereas an $L^{2}$ projection determines them by solving a mass-matrix system. The two coincide only in special cases.

### 4.3 Cell averages are the $L^{2}$ projection onto cellwise constants

The distinction becomes especially transparent for the discontinuous cellwise-constant space

$$
U_{h}
=
\{v_{h}:v_{h}\rvert_{K}\text{ is constant on every cell }K\}.
$$

Let

$$
(P_{h}f)\rvert_{K}=c_{K}.
$$

Testing the $L^{2}$ projection with the indicator of the cell $K$ gives

$$
c_{K}|K|
=
\int_{K}f \mathrm{d}x,
$$

so

$$
c_{K}
=
\frac{1}{|K|}
\int_{K}f \mathrm{d}x.
$$

The projected coefficient is the **cell average** of $f$. Sampling $f$ at the cell center would be a different approximation. For a smooth function and a symmetric cell the two may be close, but they are not the same construction.

This distinction is useful for distributed controls and coefficient fields. A cellwise control coefficient may represent an actual piecewise-constant unknown, while a cellwise approximation of a known coefficient may have been obtained by projection, averaging, or sampling.

### 4.4 Direct evaluation avoids an unnecessary discrete field

Suppose the load vector is

$$
b_{i}
=
\int_{\Omega}
f\varphi_{i} \mathrm{d}x.
$$

If $f$ is available as an analytic function, a data file, or a callback that can evaluate $f(x)$, finite-element assembly can approximate this integral directly:

$$
b_{i}
\approx
\sum_{K}
\sum_{q=1}^{n_{q}}
w_{q}^{K}
f(x_{q}^{K})
\varphi_{i}(x_{q}^{K}).
$$

No intermediate finite-element approximation of $f$ is required.

The same idea applies to a desired state $y_{d}$ in a tracking term,

$$
\frac{1}{2}
\int_{\Omega}
(y_{h}-y_{d})^{2}
 \mathrm{d}x,
$$

or to a variable coefficient $\kappa(x)$ in

$$
\int_{\Omega}
\kappa
\nabla u_{h}\cdot\nabla v_{h}
 \mathrm{d}x.
$$

Whether data should first be projected into a finite-element space is therefore a modeling and numerical decision, not an automatic step in the finite-element method.

### 4.5 A quadrature rule is a weighted evaluation formula

On a one-dimensional reference interval $\widehat K$, a quadrature rule has the form

$$
Q(g)
:=
\sum_{q=1}^{n_{q}}
\widehat w_{q}g(\widehat x_{q})
\approx
\int_{\widehat K}
g(\widehat x)
 \mathrm{d}\widehat x.
$$

The points $\widehat x_{q}$ are the **quadrature nodes** and the numbers $\widehat w_{q}$ are the **weights**.

The **degree of exactness** is the largest integer $r$ such that

$$
Q(p)
=
\int_{\widehat K}p(\widehat x) \mathrm{d}\widehat x
$$

for every polynomial $p$ of degree at most $r$.

This is more informative than the number of quadrature points alone. Different formulas with the same number of points can have different exactness properties.

### 4.6 Gaussian quadrature uses its points efficiently

For integration with unit weight on an interval, an $n_{q}$-point Gauss–Legendre rule is exact for every polynomial of degree at most

$$
2n_{q}-1.
$$

This explains why Gaussian formulas are common in finite-element assembly: polynomial basis products can often be integrated exactly with relatively few evaluations.

For example, on an affine one-dimensional cell with constant coefficients:

- a $\mathbb P_{k}$ mass entry contains a product of two degree-$k$ basis functions and therefore has degree at most $2k$;
- a $\mathbb P_{k}$ stiffness entry contains derivatives of degree at most $k-1$, so their product has degree at most $2k-2$.

Thus exact integration of the mass matrix requires a rule exact through degree $2k$, while the constant-coefficient stiffness matrix only requires exactness through degree $2k-2$.

For $\mathbb P_{1}$ elements this reproduces the observation from the main chapter: the stiffness integrand is constant, while the mass integrand is quadratic.

### 4.7 Tensor-product cells naturally use tensor-product rules

On a rectangular reference cell such as

$$
\widehat K=[-1,1]^{d},
$$

a common construction takes a one-dimensional quadrature rule in each coordinate direction and forms their tensor product.

In two dimensions,

$$
\int_{\widehat K}
g(\xi,\eta)
 \mathrm{d}\xi \mathrm{d}\eta
$$

is approximated by

$$
\sum_{q=1}^{n_{q}}
\sum_{r=1}^{n_{q}}
w_{q}w_{r}
g(\xi_{q},\eta_{r}).
$$

This matches the tensor-product structure of $\mathbb Q_{k}$ elements. Simplices use different multidimensional rules, but the same concepts of nodes, weights, and polynomial exactness apply.

### 4.8 Mapping the rule to a physical cell

Let

$$
F_{K}:\widehat K\to K
$$

be the reference-to-physical map. Then

$$
\int_{K}g(x) \mathrm{d}x
=
\int_{\widehat K}
g(F_{K}(\widehat x))
|\det J_{K}(\widehat x)|
 \mathrm{d}\widehat x.
$$

Applying the reference quadrature rule gives

$$
\int_{K}g(x) \mathrm{d}x
\approx
\sum_{q=1}^{n_{q}}
\widehat w_{q}
g(F_{K}(\widehat x_{q}))
|\det J_{K}(\widehat x_{q})|.
$$

For affine cells the Jacobian determinant is constant. For curved or more general mapped cells it varies with the quadrature point, so even a polynomial physical integrand need not remain a polynomial of the same degree after pullback.

Variable coefficients produce the same issue. If

$$
g(x)
=
\kappa(x)
\nabla\phi_{j}(x)\cdot\nabla\phi_{i}(x)
$$

and $\kappa$ is non-polynomial, no finite polynomial exactness degree makes the rule exactly integrate every such entry. The practical requirement is then accuracy sufficient for the intended finite-element approximation.

### 4.9 Finite-element approximation and numerical integration are separate choices

It is useful to keep two approximations conceptually separate:

```text
continuous problem
    ↓ choose finite-dimensional trial/test spaces
ideal Galerkin problem
    ↓ choose quadrature / approximate geometry / approximate data
implemented discrete problem
```

The first step changes the function spaces. The second changes how the variational forms are evaluated.

If quadrature is exact for the integrands that occur, the two discrete problems coincide. Otherwise the numerical rule changes the bilinear form or right-hand side. The next section explains how that additional consistency error enters the finite-element error analysis.

## 5. Quadrature and other consistency errors

The ideal Galerkin problem uses the exact bilinear form and exact right-hand side:

$$
a(u_{h},v_{h})
=
F(v_{h})
\qquad
\text{for every }v_{h}\in V_{h}.
$$

In an implementation, integrals are normally replaced by quadrature. More generally, geometry, coefficients, or operators may also be approximated. The computer may therefore solve

$$
a_{h}(u_{h},v_{h})
=
F_{h}(v_{h})
\qquad
\text{for every }v_{h}\in V_{h},
$$

where $a_{h}$ and $F_{h}$ are numerical approximations of $a$ and $F$.

The phrase **variational crime** is traditional terminology for such a modification of the ideal variational problem. The name sounds severe, but the point is simply that the error now has two sources:

```text
finite-dimensional approximation
+
consistency error in the discrete form
```

### 5.1 Why Galerkin orthogonality changes

For the exact Galerkin method,

$$
a(u-u_{h},v_{h})=0.
$$

With an approximate discrete form, the exact solution $u$ satisfies

$$
a(u,v_{h})=F(v_{h}),
$$

while the numerical solution satisfies

$$
a_{h}(u_{h},v_{h})=F_{h}(v_{h}).
$$

Subtracting no longer produces zero. Instead,

```math
\begin{aligned}
a_{h}(u-u_{h},v_{h})
&=
a_{h}(u,v_{h})-F_{h}(v_{h})
\\
&=
\bigl(a_{h}(u,v_{h})-a(u,v_{h})\bigr)
+
\bigl(F(v_{h})-F_{h}(v_{h})\bigr).
\end{aligned}
```

The right-hand side measures the **consistency error** introduced by replacing the continuous form with the numerical one.

### 5.2 A Strang-type error decomposition

Assume the discrete bilinear form is stable in the sense that, for some $\alpha_{h}>0$,

$$
a_{h}(w_{h},w_{h})
\ge
\alpha_{h}
\lVert w_{h}\rVert_{V}^{2}
\qquad
\text{for every }w_{h}\in V_{h}.
$$

Choose any comparison function $v_{h}\in V_{h}$. By the triangle inequality,

$$
\lVert u-u_{h}\rVert_{V}
\le
\lVert u-v_{h}\rVert_{V}
+
\lVert v_{h}-u_{h}\rVert_{V}.
$$

Discrete coercivity controls the second term:

$$
\alpha_{h}
\lVert v_{h}-u_{h}\rVert_{V}^{2}
\le
 a_{h}(v_{h}-u_{h},v_{h}-u_{h}).
$$

Using the discrete equation for $u_{h}$ and adding and subtracting the exact form shows that this difference is driven by two kinds of terms:

- how far $v_{h}$ is from $u$;
- how far $a_{h}$ and $F_{h}$ are from $a$ and $F$ when acting on relevant functions.

After taking a supremum over discrete test functions, one obtains a bound of the schematic form

```math
\lVert u-u_{h}\rVert_{V}
\lesssim
\inf_{v_{h}\in V_{h}}
\lVert u-v_{h}\rVert_{V}
+
\sup_{w_{h}\in V_{h}}
\frac{
\lvert a_{h}(u,w_{h})-F_{h}(w_{h})\rvert
}{
\lVert w_{h}\rVert_{V}
}.
```

The first term is approximation error. The second is consistency error. This is the central message of Strang's lemmas.

If the quadrature is exact for all integrands that occur in the assembled problem, then $a_{h}=a$ and $F_{h}=F$ on the discrete spaces and the consistency term disappears. When exactness is impossible, the quadrature should be accurate enough that this extra error does not dominate the finite-element approximation error.

## 6. A posteriori error estimation and adaptive refinement

The interpolation estimates above are **a priori**: they predict the error from assumptions about the exact solution, even though the exact solution is precisely what we do not know.

An **a posteriori** estimator instead uses the computed solution and known problem data to estimate where the current approximation is inaccurate.

### 6.1 The residual of the computed solution

For the Poisson problem

$$
-\Delta u=f,
$$

consider a conforming finite-element solution $u_{h}$. Inside a cell $K$, define the strong cell residual

$$
r_{K}
:=
f+\Delta u_{h}.
$$

For piecewise-affine elements, $\Delta u_{h}=0$ inside each cell, so the cell residual reduces to the forcing there. This does **not** mean the finite-element error is determined only by $f$. The missing information appears at cell interfaces.

The gradient of $u_{h}$ is generally different on the two sides of an interior face $F$. If $n^{+}$ and $n^{-}$ denote the outward normals of its neighboring cells, a normal-flux jump can be written as

$$
j_{F}
:=
\nabla u_{h}^{+}\cdot n^{+}
+
\nabla u_{h}^{-}\cdot n^{-}.
$$

For an exact sufficiently smooth solution of the homogeneous-diffusion Poisson equation, the corresponding interior flux contributions match, so this jump vanishes. A nonzero discrete jump therefore carries information about the error.

### 6.2 Why the two residuals appear

The origin of these terms is another cellwise integration by parts. For a test function $v$,

```math
\sum_{K}
\int_{K}
\nabla(u-u_{h})\cdot\nabla v
\mathrm{d}x
=
\sum_{K}
\left[
\int_{K}
(f+\Delta u_{h})v
\mathrm{d}x
-
\int_{\partial K}
(\nabla u_{h}\cdot n_{K})v
\mathrm{d}s
\right].
```

When the cell-boundary terms are collected globally, external boundary contributions are treated by the boundary conditions, while an interior face receives contributions from both neighboring cells. Those two contributions combine into the flux jump $j_{F}$.

A typical residual indicator therefore has the structure

$$
\eta_{K}^{2}
=
h_{K}^{2}
\lVert r_{K}\rVert_{L^{2}(K)}^{2}
+
\sum_{F\subset\partial K\cap\Omega}
h_{F}
\lVert j_{F}\rVert_{L^{2}(F)}^{2},
$$

up to convention-dependent constants and boundary terms. Here $h_{F}$ is a characteristic size of the face.

The powers of $h_{K}$ and $h_{F}$ compensate for the different dimensions and derivative orders of the residual terms. A complete proof that such an estimator is reliable and locally efficient requires interpolation and localization arguments and is beyond this note.

### 6.3 The adaptive loop

Once local indicators are available, an adaptive method follows a loop such as

```text
SOLVE
  ↓
ESTIMATE local errors η_K
  ↓
MARK cells with significant contribution
  ↓
REFINE those cells
  ↓
repeat
```

Local refinement destroys quasi-uniformity by design, but it can preserve shape regularity.

It can also create **hanging nodes**: fine cells may introduce support points on an interface that do not exist on a neighboring coarse cell. For a conforming method, the corresponding fine-grid degree of freedom cannot remain independent. It must satisfy an interpolation relation such as

$$
y_{c}
=
\frac12 y_{1}
+
\frac12 y_{2},
$$

which is exactly the kind of affine coordinate constraint introduced in the main chapter.

## 7. What changes when the discrete space is discontinuous

A cellwise-discontinuous space can be perfectly natural for a quantity whose continuous space is $L^{2}(\Omega)$. That is why piecewise-constant controls were unproblematic in the main chapter.

The situation changes if we try to approximate an $H^{1}$ PDE solution with functions that are allowed to jump across cell faces. Then

$$
V_{h}\not\subset H^{1}(\Omega),
$$

so the method is no longer a conforming restriction of the original weak problem.

### 7.1 Broken derivatives are not enough

A discontinuous piecewise-polynomial function still has an ordinary gradient inside each cell. We can therefore define the **broken gradient** by

$$
\nabla_{h}v_{h}\rvert_{K}
:=
\nabla(v_{h}\rvert_{K}).
$$

The volume expression

$$
\sum_{K}
\int_{K}
\nabla_{h}u_{h}\cdot\nabla_{h}v_{h}
\mathrm{d}x
$$

is well defined. But it does not by itself communicate how neighboring cells should interact. If two cells share a face, their traces are independent and may disagree.

Define the scalar jump across an interior face by choosing the two traces $v_{h}^{+}$ and $v_{h}^{-}$ and setting

$$
[v_{h}]
:=
v_{h}^{+}-v_{h}^{-}.
$$

For a continuous finite-element function this jump is zero. For a discontinuous function it becomes additional information that the method must control.

### 7.2 Interface terms restore communication

Interior-penalty discontinuous Galerkin methods augment the cellwise gradient form with face terms. Schematically, a symmetric interior-penalty form contains

```math
\begin{aligned}
a_{h}(u_{h},v_{h})
={}&
\sum_{K}
\int_{K}
\nabla u_{h}\cdot\nabla v_{h}
\mathrm{d}x
\\
&-
\sum_{F}
\int_{F}
\{\partial_{n}u_{h}\}[v_{h}]
\mathrm{d}s
-
\sum_{F}
\int_{F}
\{\partial_{n}v_{h}\}[u_{h}]
\mathrm{d}s
\\
&+
\sum_{F}
\frac{\sigma}{h_{F}}
\int_{F}
[u_{h}][v_{h}]
\mathrm{d}s.
\end{aligned}
```

Here $\{\partial_{n}u_{h}\}$ denotes a suitable average of the two neighboring normal derivatives, and $\sigma>0$ is a penalty parameter.

The first face terms encode the flux coupling produced by integration by parts. The penalty term controls the jump itself. With an appropriate scaling and sufficiently large penalty, the resulting method can be stable even though the discrete functions are not members of the continuous $H^{1}$ space.

This is a qualitatively different use of a discontinuous space from a cellwise-constant $L^{2}$ control. In the latter case discontinuity is already compatible with the continuous space. In the PDE case the numerical formulation has to compensate for the loss of $H^{1}$ conformity.

## 8. Petrov–Galerkin and mixed finite elements

The conforming scalar Galerkin problem used the same space for trial and test functions. More general problems may naturally use different spaces.

Let $V$ be the trial space and $W$ the test space. A variational problem can have the form

$$
a(u,w)=F(w)
\qquad
\text{for every }w\in W,
$$

with

$$
a:V\times W\to\mathbb R.
$$

A discrete method chooses

$$
V_{h}\subset V,
\qquad
W_{h}\subset W,
$$

and seeks $u_{h}\in V_{h}$ such that

$$
a(u_{h},w_{h})=F(w_{h})
\qquad
\text{for every }w_{h}\in W_{h}.
$$

When the trial and test spaces differ, this is commonly called a **Petrov–Galerkin** method.

### 8.1 Coercivity is replaced by an inf-sup condition

When $V=W$ and $a$ is coercive, testing with the same function controls its norm directly:

$$
a(v,v)
\ge
\alpha
\lVert v\rVert_{V}^{2}.
$$

If $V$ and $W$ are different, the expression $a(v,v)$ may not even be meaningful. Stability is instead measured by whether every nonzero trial direction can be detected strongly enough by some test direction.

At the discrete level, a representative inf-sup condition is

$$
\inf_{v_{h}\in V_{h},v_{h}\ne0}
\sup_{w_{h}\in W_{h},w_{h}\ne0}
\frac{
\lvert a(v_{h},w_{h})\rvert
}{
\lVert v_{h}\rVert_{V}
\lVert w_{h}\rVert_{W}
}
\ge
\beta_{h}>0.
$$

For a stable family of discretizations, one wants a lower bound

$$
\beta_{h}\ge\beta_{0}>0
$$

that does not collapse as $h\to0$.

This is the analogue of requiring mesh-independent coercivity. A discrete problem can be perfectly finite dimensional and still become unstable if its inf-sup constant tends to zero under refinement.

### 8.2 Mixed formulations produce block systems

A common mixed problem introduces two unknown fields $v\in V$ and $\pi\in Q$ and asks for

```math
\begin{aligned}
a(v,w)+b(w,\pi) &= F(w)
&&\text{for every }w\in V,
\\
b(v,q) &= G(q)
&&\text{for every }q\in Q.
\end{aligned}
```

Choose bases

$$
V_{h}=\mathrm{span}\{\varphi_{1},\ldots,\varphi_{n}\},
\qquad
Q_{h}=\mathrm{span}\{\psi_{1},\ldots,\psi_{m}\}.
$$

The two bilinear forms produce matrices

$$
A_{ij}=a(\varphi_{j},\varphi_{i}),
$$

and

$$
B_{kj}=b(\varphi_{j},\psi_{k}).
$$

The discrete equations have the block form

```math
\begin{bmatrix}
A & B^{\mathsf T}\\
B & 0
\end{bmatrix}
\begin{bmatrix}
\mathbf v\\
\boldsymbol{\pi}
\end{bmatrix}
=
\begin{bmatrix}
\mathbf f\\
\mathbf g
\end{bmatrix}.
```

This matrix is generally indefinite. Its well-posedness is not determined by the quality of $V_{h}$ and $Q_{h}$ separately. The **pair** of spaces matters through a discrete inf-sup condition for $b$.

This is why mixed finite-element literature speaks of stable pairs. For incompressible-flow problems, for example, Taylor–Hood velocity/pressure spaces are a standard stable choice. Other apparently reasonable combinations may require stabilization or may be unsuitable.

The numerical linear-algebra consequences of this block structure – indefinite systems, Schur complements, MINRES/GMRES, and block preconditioning – belong to [04 · Numerical linear algebra for PDE and optimization systems](04-numerical-linear-algebra.md).

## 9. What these extensions change in the mental model

The simple conforming scalar method can be summarized as

```text
choose V_h ⊂ V
        ↓
restrict the same variational problem
        ↓
assemble exact or sufficiently accurate cell contributions
        ↓
solve
```

The extensions in this note modify different parts of that picture:

- **interpolation theory** explains how approximation quality depends on $h$, degree, regularity, and mesh shape;
- **quadrature and variational crimes** change the bilinear form or functional, so approximation error is joined by consistency error;
- **adaptivity** changes the mesh nonuniformly in response to information extracted from the computed solution;
- **discontinuous or nonconforming methods** use spaces outside the original energy space and therefore add interface or consistency machinery;
- **Petrov–Galerkin and mixed methods** change the relation between trial and test spaces and replace simple coercivity by inf-sup stability.

These are not unrelated complications. Each one answers a specific limitation of the baseline method.

## Where to go next

- **Return to the main route:** continue with [04 · Numerical linear algebra for PDE and optimization systems](04-numerical-linear-algebra.md), which studies the sparse and block systems produced by finite-element discretizations.
- **For implementation:** [09 · The deal.II finite-element workflow](09-dealii-finite-element-workflow.md) shows where quadrature, hanging-node constraints, assembly, and solver choices appear in deal.II code.
- **For later advanced applications:** the mixed, adaptive, and nonconforming ideas introduced here become relevant to more complicated PDE systems and discretizations, especially those surveyed in [07b · Extensions of PDE-constrained optimization](07b-extensions-of-pde-constrained-optimization.md).

## References and further reading

Publication details, stable links, and access notes are collected in the [Background reference catalogue](references.md).

- **Stefano Mancini, Andrea Benvenuti, and Luca Heltai, *Numerical Methods for Partial Differential Equations*.**
  - §2.1 for reference-to-physical-cell scaling and the role of cell geometry.
  - §2.2, **Bramble–Hilbert lemma**, and §2.3, **$H^{k}$ error estimates**, for interpolation estimates, shape regularity, and the passage from Céa's lemma to powers of $h$.
  - §3.1 for Scott–Zhang interpolation and projection operators.
  - §§3.2–3.3 for a posteriori error analysis, residual estimators, and adaptivity.
  - §4.1 for Strang's lemmas and consistency error.
  - §§4.2–4.3 for discontinuous Galerkin and Nitsche-type methods.
  - §§5.1–5.3 for BNB/inf-sup conditions, Petrov–Galerkin approximation, and mixed problems.
- **Manzoni, Quarteroni, and Salsa, *Optimal Control of Partial Differential Equations*.**
  - Appendix B, §B.2.1 for Galerkin orthogonality and Céa's estimate.
  - Appendix B, §§B.2.3–B.2.4 for saddle-point problems, discrete inf-sup conditions, and block algebra.
  - Appendix B, §B.4 for finite-element spaces and interpolation.
  - Appendix B, §B.5 for regularity-dependent a priori convergence estimates.
- **Alexandre Ern and Jean-Luc Guermond, *Theory and Practice of Finite Elements*.**
  - **Finite Element Interpolation** continues the approximation material.
  - **Quadratures, Assembling, and Storage** expands the numerical-integration discussion from an implementation-aware viewpoint.
- **Susanne C. Brenner and L. Ridgway Scott, *The Mathematical Theory of Finite Element Methods*.**
  - The chapters on polynomial approximation, adaptive meshes, variational crimes, and mixed methods align with the corresponding sections of this note.
- **Daniele Boffi, Franco Brezzi, and Michel Fortin, *Mixed Finite Element Methods and Applications*.**
  - Recommended when saddle-point and inf-sup theory becomes a primary topic rather than orientation.
- **Philippe G. Ciarlet, *The Finite Element Method for Elliptic Problems*.**
  - Classical reference for the conforming elliptic theory underlying the approximation sections.
- **Bini, Boito, and Meini, *Appunti di Istituzioni di Analisi Numerica*.**
  - Chapter 2, §§2.1–2.5, for a fuller treatment of interpolatory quadrature, convergence, Newton–Cotes, Fejér/Clenshaw–Curtis, and Gaussian formulas.
  - This complements the finite-element-oriented discussion here; the notes are course-local and mainly useful to University of Pisa students.
