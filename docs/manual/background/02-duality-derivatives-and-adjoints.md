# 02 · Duality, derivatives, and adjoints

**Background navigation:** [Index](README.md) \
Previous: [01 · Function spaces and weak PDEs](01-function-spaces-and-weak-pdes.md) \
Next: [03 · Finite elements](03-finite-elements.md)

## Purpose

A vector and a derivative can be stored as the same list of numbers and still play different mathematical roles. A perturbation is something we can *apply* to a state or parameter. A derivative of a scalar functional is something that *acts on* such a perturbation. An operator sends primal objects forward; its transpose or adjoint sends dual information backward. An inner product can identify a dual object with a primal representative, but that identification depends on the chosen geometry.

These distinctions are easy to hide in Euclidean coordinates because the standard dot product identifies $\mathbb{R}^{n}$ with its dual almost invisibly. They become unavoidable in function spaces, finite elements, and PDE-constrained optimization.

This chapter develops the underlying mathematics before any project-specific representation is introduced. We begin with linear functionals and dual coordinates, derive how vectors and covectors transform differently, introduce Riesz maps, and then extend ordinary differentiation to maps between normed spaces. The chain rule leads naturally to pullbacks and adjoint actions. Finally, we distinguish derivatives from gradients, introduce second derivatives as operator actions, and derive the numerical checks that later let us verify derivative and adjoint implementations.

The scope is deliberately practical. We use Banach- and Hilbert-space language where it clarifies the operations that appear later, but we do not develop general operator theory, automatic differentiation, or the full differential calculus of nonlinear operators on arbitrary Banach spaces.

## Before you start

This chapter assumes:

- the function-space and bounded-functional vocabulary from [01 · Function spaces and weak PDEs](01-function-spaces-and-weak-pdes.md), or equivalent familiarity;
- ordinary finite-dimensional vectors, matrices, bases, and transposes;
- multivariable differentiation and the finite-dimensional chain rule.

No prior familiarity with covectors, Riesz maps, Fréchet derivatives, or functional adjoints is assumed. All vector spaces in this chapter are real; complex inner-product spaces require conjugation conventions that are deliberately outside the present scope.

## What you will be able to do

After this chapter, you should be able to:

- distinguish a vector from a linear functional acting on vectors;
- represent vectors and covectors in a basis and explain why their coordinates transform differently;
- distinguish a dual pairing from an inner product;
- explain what the Riesz map does and why it depends on the chosen inner product;
- interpret a Fréchet derivative as a bounded linear first-order model;
- recover directional derivatives from a Fréchet derivative;
- derive the chain rule for compositions of maps between normed spaces;
- define a dual or adjoint action through a pairing identity;
- explain why the derivative of a scalar functional is naturally a covector, while a gradient is a metric-dependent primal representative;
- interpret Hessian actions without requiring an explicitly assembled Hessian matrix;
- use Taylor remainders and pairing identities to check derivative and adjoint calculations numerically.

## Roadmap

We first revisit linear maps and linear functionals, then make the vector/covector distinction concrete using bases and coordinate changes. An inner product gives additional structure, leading to the Riesz map and to the familiar matrix representation $Gg=r$. We next introduce directional and Fréchet derivatives, prove the chain rule at the level needed later, and use it to explain push-forwards and pullbacks. Scalar functionals then provide the cleanest setting for separating derivative from gradient. The last part of the chapter introduces adjoint derivative actions, second derivatives, Hessian actions, and practical derivative checks.

The later manual uses all of these ideas in discrete spaces. This chapter stops before project-specific contracts, finite-element assembly, optimization algorithms, and automatic-differentiation implementation.

## 1. Linear maps and linear functionals

Let $X$ and $Y$ be normed vector spaces. A map

$$
A\colon X\to Y
$$

is **linear** if

$$
A(\alpha x+\beta z)=\alpha Ax+\beta Az
$$

for all $x,z\in X$ and scalars $\alpha,\beta$.

For finite-dimensional spaces every linear map is continuous. In infinite-dimensional spaces this is no longer automatic, so we explicitly require the linear maps used in analysis to be **bounded**: there exists a constant $C\geq 0$ such that

$$
\lVert Ax\rVert_{Y}\leq C\lVert x\rVert_{X}
\qquad
\text{for every }x\in X.
$$

We write

$$
\mathcal L(X,Y)
$$

for the space of bounded linear maps from $X$ to $Y$. Its standard operator norm is

$$
\lVert A\rVert_{\mathcal L(X,Y)}\coloneqq
\sup_{x\ne0}
\frac{\lVert Ax\rVert_{Y}}{\lVert x\rVert_{X}}=
\sup_{\lVert x\rVert_{X}=1}
\lVert Ax\rVert_{Y}.
$$

Boundedness and continuity are equivalent for linear maps. One direction is immediate: if $A$ is bounded, then

$$
\lVert Ax-Az\rVert_{Y}=\lVert A(x-z)\rVert_{Y}\leq
\lVert A\rVert_{\mathcal L(X,Y)}
\lVert x-z\rVert_{X},
$$

so $A$ is in fact Lipschitz continuous.

For the converse, suppose $A$ is continuous at the origin. Taking $\varepsilon=1$, there exists $\delta\gt 0$ such that

$$
\lVert z\rVert_{X}\lt\delta
\quad\Longrightarrow\quad
\lVert Az\rVert_{Y}\lt 1.
$$

For any nonzero $x$, choose

$$
z\coloneqq
\frac{\delta}{2\lVert x\rVert_{X}}x.
$$

Then $\lVert z\rVert_{X}=\delta/2\lt\delta$, so $\lVert Az\rVert_{Y}\lt 1$. By linearity,

$$
\frac{\delta}{2\lVert x\rVert_{X}}\lVert Ax\rVert_{Y}\lt 1,
$$

and therefore

$$
\lVert Ax\rVert_{Y}
\lt\frac{2}{\delta}\lVert x\rVert_{X}.
$$

Thus continuity at one point already forces a linear map to be bounded everywhere.

A **linear functional** on $X$ is a linear map

$$
\ell\colon X\to\mathbb{R}.
$$

When it is bounded, it belongs to the **dual space**

$$
X^{\ast}\coloneqq\mathcal L(X,\mathbb{R}).
$$

An element of $X$ is often called a **primal vector** or simply a vector. An element of $X^{\ast}$ is often called a **covector**. The word *covector* emphasizes the role of the object: it consumes a vector and returns a scalar. A row vector is only one coordinate representation of a covector; the mathematical object is the functional itself.

The action of $\ell\in X^{\ast}$ on $x\in X$ is written

$$
\langle \ell,x\rangle_{X^{\ast},X}\coloneqq\ell(x).
$$

This scalar evaluation is the **dual pairing**. It requires no inner product.

For example, on $X=\mathbb{R}^{2}$ the rule

$$
\ell(x_{1},x_{2})\coloneqq 3x_{1}-2x_{2}
$$

is a linear functional. In the standard basis we may store it as the coefficient array

$$
r =
\begin{bmatrix}
3\\
-2
\end{bmatrix},
$$

so that

$$
\ell(x)=r^{\mathsf T}x.
$$

The same storage shape as a vector does not make $r$ the same mathematical kind of object as $x$. The vector is an input to the functional; the covector records how the functional acts on such inputs.

## 2. Bases, dual bases, and coordinate changes

The distinction becomes sharper when the basis changes.

Let $X$ be an $n$-dimensional vector space with basis

$$
\lbrace\phi_{1},\ldots,\phi_{n}\rbrace.
$$

Every vector $x\in X$ has coordinates $x_{j}$ defined by

$$
x =
\sum_{j=1}^{n}x_{j}\phi_{j}.
$$

The associated **dual basis**

$$
\lbrace\phi^{1},\ldots,\phi^{n}\rbrace
\subset X^{\ast}
$$

is defined by

$$
\phi^{i}(\phi_{j})=\delta_{ij},
$$

where $\delta_{ij}$ is the Kronecker delta: it is $1$ when $i=j$ and $0$ otherwise.

Every covector $\ell\in X^{\ast}$ can therefore be expanded as

$$
\ell=\sum_{i=1}^{n}r_{i}\phi^{i}.
$$

Applying it to $x$ gives

$$
\begin{aligned}
\ell(x)
&=
\sum_{i=1}^{n}r_{i}\phi^{i}
\left(
\sum_{j=1}^{n}x_{j}\phi_{j}
\right)
\\
&=
\sum_{i=1}^{n}
\sum_{j=1}^{n}
r_{i}x_{j}\delta_{ij}
\\
&=
\sum_{i=1}^{n}r_{i}x_{i}.
\end{aligned}
$$

Thus the familiar coordinate formula

$$
\ell(x)=r^{\mathsf T}x
$$

is not itself an inner product. It is the coordinate expression of the natural pairing between a dual basis and its primal basis.

### A basis change

Suppose a second basis $\lbrace\widetilde\phi_{1},\ldots,\widetilde\phi_{n}\rbrace$ is related to the first by an invertible matrix $P$:

$$
\widetilde\phi_{i} =
\sum_{j=1}^{n}\phi_{j}P_{ji}.
$$

The columns of $P$ contain the new basis vectors written in the old basis.

Let $x$ and $\widetilde x$ be the coordinate arrays of the same vector in the old and new bases. Since

$$
\sum_{j=1}^{n}x_{j}\phi_{j} =
\sum_{i=1}^{n}\widetilde x_{i}\widetilde\phi_{i},
$$

substituting the definition of $\widetilde\phi_{i}$ gives

$$
x=P\widetilde x,
$$

hence

$$
\widetilde x=P^{-1}x.
$$

Now let $r$ and $\widetilde r$ be the corresponding coordinate arrays of a fixed covector $\ell$. The scalar pairing cannot depend on which basis we use:

$$
r^{\mathsf T}x =
\widetilde r^{\mathsf T}\widetilde x.
$$

Using $x=P\widetilde x$,

$$
r^{\mathsf T}P\widetilde x =
\widetilde r^{\mathsf T}\widetilde x
$$

for every $\widetilde x$. Therefore

$$
\widetilde r =
P^{\mathsf T}r.
$$

The same change of basis thus acts differently on the two coordinate arrays:

```text
primal coordinates
    x  ──►  P^{-1} x

covector coordinates
    r  ──►  P^{T} r
```

A small example makes the distinction concrete. In $\mathbb{R}^{2}$, let

$$
\widetilde e_{1}=e_{1}+e_{2},
\qquad
\widetilde e_{2}=e_{2}.
$$

Then

$$
P=
\begin{bmatrix}
1&0\\
1&1
\end{bmatrix}.
$$

The vector

$$
x=e_{1}+2e_{2}
$$

has old coordinates

$$
x=
\begin{bmatrix}
1\\
2
\end{bmatrix}
$$

and new coordinates

$$
\widetilde x =
P^{-1}x =
\begin{bmatrix}
1\\
1
\end{bmatrix}.
$$

Now take the covector

$$
\ell(x_{1},x_{2})=3x_{1}+4x_{2},
$$

whose old dual coordinates are

$$
r=
\begin{bmatrix}
3\\
4
\end{bmatrix}.
$$

Its new dual coordinates are

$$
\widetilde r =
P^{\mathsf T}r =
\begin{bmatrix}
7\\
4
\end{bmatrix}.
$$

The pairing is unchanged:

$$
r^{\mathsf T}x
=11,
\qquad
\widetilde r^{\mathsf T}\widetilde x
=11.
$$

This is the basic reason primal perturbations and derivative covectors cannot be treated as interchangeable arrays merely because they have equal length.

## 3. Dual pairing versus inner product

A dual pairing and an inner product both produce scalars, but they answer different questions.

For $\ell\in X^{\ast}$ and $x\in X$,

$$
\langle \ell,x\rangle_{X^{\ast},X}
$$

means: **apply this functional to this vector**.

An inner product on a real vector space $H$ is instead a map

$$
(\cdot,\cdot)_{H}\colon H\times H\to\mathbb{R}
$$

that is bilinear, symmetric, and positive definite. It equips the primal space with geometry: lengths, angles, orthogonality, and projections.

In $\mathbb{R}^{n}$ the standard Euclidean inner product is

$$
(x,z)_{2}=x^{\mathsf T}z.
$$

That formula looks identical to the coordinate pairing $r^{\mathsf T}x$. The roles are different:

- in $r^{\mathsf T}x$, $r$ stores a covector and $x$ stores a vector;
- in $x^{\mathsf T}z$, both $x$ and $z$ store vectors, and the Euclidean inner product supplies the identification that lets one vector act like a covector.

Other inner products produce different identifications. Let $G\in\mathbb{R}^{n\times n}$ be symmetric positive definite. Then

$$
(x,z)_{G}
\coloneqq
x^{\mathsf T}Gz
$$

is also an inner product. The vector $x$ now acts on $z$ through the covector with coordinate array $Gx$, because

$$
(x,z)_{G} =
(Gx)^{\mathsf T}z.
$$

Thus there is no basis-independent rule saying that the numbers stored in a vector are automatically the numbers stored in its corresponding covector. Such a rule appears only after a particular inner product has been chosen.

## 4. Riesz representation and the Riesz map

Let $H$ be a real Hilbert space. We will use the Riesz representation theorem without proving it; its proof belongs to the functional-analysis material deliberately left outside this chapter. The theorem states that every bounded linear functional $\ell\in H^{\ast}$ can be represented uniquely by an element $g\in H$ such that

$$
\langle \ell,v\rangle_{H^{\ast},H} =
(g,v)_{H}
\qquad
\text{for every }v\in H.
$$

The theorem is important because it gives a precise way to move between a dual object and a primal representative once an inner product has been fixed.

For later use, define the **Riesz map** in the primal-to-dual direction,

$$
R_{H}\colon H\to H^{\ast},
$$

by

$$
\langle R_{H}g,v\rangle_{H^{\ast},H}
\coloneqq
(g,v)_{H}.
$$

The Riesz theorem says that $R_{H}$ is a one-to-one and onto isometry. Its inverse

$$
R_{H}^{-1}\colon H^{\ast}\to H
$$

assigns to each functional its unique Riesz representative.

Some texts call the inverse map $H^{\ast}\to H$ the Riesz map instead. The convention is not mathematically significant as long as the direction is stated. Here we use the primal-to-dual direction because it leads directly to the relation used later for gradients:

$$
R_{H}g=\ell.
$$

### Coordinate form: the Gram matrix

Let $H_{h}$ be finite dimensional with basis $\lbrace\phi_{1},\ldots,\phi_{n}\rbrace$. Define the Gram matrix of the inner product by

$$
G_{ij}
\coloneqq
(\phi_{i},\phi_{j})_{H}.
$$

For

$$
g =
\sum_{i=1}^{n}g_{i}\phi_{i},
\qquad
v =
\sum_{j=1}^{n}v_{j}\phi_{j},
$$

we obtain

$$
\begin{aligned}
(g,v)_{H}
&=
\sum_{i=1}^{n}
\sum_{j=1}^{n}
g_{i}v_{j}(\phi_{i},\phi_{j})_{H}
\\
&=
g^{\mathsf T}Gv
\\
&=
(Gg)^{\mathsf T}v,
\end{aligned}
$$

where the last equality uses the symmetry of $G$.

Therefore the covector coordinates of $R_{H}g$ are

$$
r=Gg.
$$

Conversely, if a covector has coordinates $r$, its Riesz representative is obtained from

$$
Gg=r.
$$

The Euclidean case is only the special choice $G=I$.

### Function-space example

In $H=L^{2}(\Omega)$,

$$
(g,v)_{L^{2}} =
\int_{\Omega}gv \thinspace\mathrm{d}x.
$$

Thus an $L^{2}$ function $g$ represents the functional

$$
v
\longmapsto
\int_{\Omega}gv \thinspace\mathrm{d}x.
$$

For any finite-dimensional function basis $\lbrace\phi_{i}\rbrace$, the $L^{2}$ Gram matrix has entries

$$
M_{ij}
\coloneqq
\int_{\Omega}\phi_{i}\phi_{j} \thinspace\mathrm{d}x.
$$

In finite-element terminology, this is the **mass matrix**. The same abstract Riesz map that is invisible in Euclidean coordinates therefore becomes a nontrivial operator after discretization.

## 5. Push-forwards, pullbacks, and adjoint linear maps

Let

$$
A\colon X\to Y
$$

be a bounded linear map. A primal vector moves **forward** through $A$:

$$
x\in X
\quad\longmapsto\quad
Ax\in Y.
$$

A covector on $Y$ moves in the opposite direction by composition. Given $q\in Y^{\ast}$, define

$$
A^{\ast}q\in X^{\ast}
$$

by

$$
\langle A^{\ast}q,x\rangle_{X^{\ast},X}
\coloneqq
\langle q,Ax\rangle_{Y^{\ast},Y}
\qquad
\text{for every }x\in X.
$$

The map

$$
A^{\ast}\colon Y^{\ast}\to X^{\ast}
$$

is the **dual**, **transpose**, or **Banach adjoint** of $A$. The terminology varies by source; the pairing identity is the definition that matters.

In coordinates, let $A$ be represented by a matrix $\mathbf A$. If $q$ has dual-coordinate vector $\mathbf q$, then

$$
\begin{aligned}
\langle q,Ax\rangle
&=
\mathbf q^{\mathsf T}\mathbf A\mathbf x
\\
&=
(\mathbf A^{\mathsf T}\mathbf q)^{\mathsf T}\mathbf x.
\end{aligned}
$$

Therefore

$$
A^{\ast}q
\quad\longleftrightarrow\quad
\mathbf A^{\mathsf T}\mathbf q.
$$

The matrix transpose is not an extra rule imposed on coordinates. It is what the basis representation of the pairing identity becomes.

### The Hilbert adjoint

If $X$ and $Y$ are Hilbert spaces, their Riesz maps allow us to turn the dual operator into a primal-to-primal map. Define

$$
A^{\dagger}
\coloneqq
R_{X}^{-1}A^{\ast}R_{Y}\colon Y\to X.
$$

Then

$$
\begin{aligned}
(Ax,y)_{Y}
&=
\langle R_{Y}y,Ax\rangle
\\
&=
\langle A^{\ast}R_{Y}y,x\rangle
\\
&=
(R_{X}^{-1}A^{\ast}R_{Y}y,x)_{X}
\\
&=
(x,A^{\dagger}y)_{X}.
\end{aligned}
$$

Thus $A^{\dagger}$ is the usual Hilbert-space adjoint. Many texts denote both the dual operator and the Hilbert adjoint by $A^{\ast}$. Keeping the Riesz maps visible makes clear which object is being used.

No inner product is required to define the dual pullback $A^{\ast}\colon Y^{\ast}\to X^{\ast}$. Inner products enter only when we additionally want a primal representative in $X$ or $Y$.

## 6. Directional derivatives and Fréchet derivatives

A nonlinear map is not itself a linear operator, but near a point we may approximate its change by a linear map.

Let $X$ and $Y$ be normed spaces, let $U\subset X$ be open, and let

$$
F\colon U\to Y.
$$

The openness assumption is practical: if $x\in U$, then sufficiently small perturbations $x+h$ remain inside the domain, so it makes sense to ask how $F$ changes in arbitrary nearby directions.

### Directional derivative

For a direction $d\in X$, the **directional derivative** of $F$ at $x$ in direction
$d$ is

$$
D F(x)[d]
\coloneqq
\lim_{t\to0}
\frac{F(x+td)-F(x)}{t},
$$

provided the limit exists in $Y$.

Here $d$ specifies the direction, while $td$ is the actual perturbation of $x$.
This convention will be useful below: we use $d$ for a fixed direction and $h$ for
a generic small displacement, with $h=td$ when we restrict that displacement to a
line.

A directional derivative examines one line through $x$ at a time. The existence of
directional derivatives in every direction does not by itself guarantee that those
directional changes assemble into one bounded linear first-order model.

### Fréchet derivative

The map $F$ is **Fréchet differentiable** at $x$ if there exists a bounded linear operator

$$
F'(x)\in\mathcal L(X,Y)
$$

such that

$$
F(x+h) =
F(x)+F'(x)h+r(h),
$$

where the remainder satisfies

$$
\frac{\lVert r(h)\rVert_{Y}}{\lVert h\rVert_{X}}
\longrightarrow 0
\qquad
\text{as }\lVert h\rVert_{X}\to0.
$$

The notation

$$
r(h)=o(\lVert h\rVert_{X})
$$

is shorthand for exactly this limit. We will write $F'(x)h$ and $F'(x)[h]$ interchangeably for the action of the derivative on a direction.

The Fréchet derivative is unique. If two bounded linear maps $A$ and $B$ both satisfied the definition, then

$$
(A-B)h=o(\lVert h\rVert_{X}).
$$

Set $h=td$ for a fixed direction $d$, divide by $\lvert t\rvert$, and let $t\to0$. Linearity gives

$$
\lVert (A-B)d\rVert_{Y}=0,
$$

so $Ad=Bd$ for every $d\in X$, hence $A=B$. The Fréchet derivative is therefore the unique bounded linear map that captures the entire first-order change of $F$ for all sufficiently small perturbations at once.

It automatically produces every directional derivative. Set $h=t d$ for a fixed direction $d$. Then

$$
F(x+td)-F(x) =
tF'(x)d+r(td).
$$

Dividing by $t\ne0$ gives

$$
\frac{F(x+td)-F(x)}{t} =
F'(x)d+
\frac{r(td)}{t}.
$$

Moreover,

$$
\left\lVert
\frac{r(td)}{t}
\right\rVert_{Y} =
\frac{\lVert r(td)\rVert_{Y}}{\lvert t\rvert\lVert d\rVert_{X}}
\lVert d\rVert_{X}
\longrightarrow0.
$$

Hence

$$
D F(x)[d] =
F'(x)d.
$$

So Fréchet differentiability is stronger than merely being able to differentiate along individual lines.

The term **Gâteaux derivative** is often used for a directional derivative that depends linearly and boundedly on the direction. It is useful in infinite-dimensional analysis, but the distinction between Gâteaux and Fréchet differentiability can be subtle. For the main path here, we use Fréchet differentiability whenever a genuine first-order linear model is needed.

### A nonlinear example

Consider

$$
F\colon\mathbb{R}^{2}\to\mathbb{R}^{2},
\qquad
F(x_{1},x_{2}) =
\begin{bmatrix}
x_{1}^{2}+x_{2}\\
\sin x_{1}
\end{bmatrix}.
$$

Perturb $x$ by $h=(h_{1},h_{2})$. Then

$$
\begin{aligned}
F(x+h)-F(x)
&=
\begin{bmatrix}
(x_{1}+h_{1})^{2}+x_{2}+h_{2}-x_{1}^{2}-x_{2}\\
\sin(x_{1}+h_{1})-\sin x_{1}
\end{bmatrix}
\\
&=
\begin{bmatrix}
2x_{1}h_{1}+h_{2}\\
(\cos x_{1})h_{1}
\end{bmatrix}
+
\begin{bmatrix}
h_{1}^{2}\\
\sin(x_{1}+h_{1})-\sin x_{1}-(\cos x_{1})h_{1}
\end{bmatrix}.
\end{aligned}
$$

The first vector is linear in $h$. The second is of smaller order than $\lVert h\rVert$: its first component is quadratic, and the scalar Taylor expansion of sine gives a quadratic remainder in $h_{1}$. Therefore

$$
F'(x) =
\begin{bmatrix}
2x_{1}&1\\
\cos x_{1}&0
\end{bmatrix}.
$$

At a fixed point $x$, the derivative is a linear operator in the perturbation $h$, even though the original map $F$ is nonlinear in $x$.

## 7. The chain rule

Let

$$
F\colon U\subset X\to Y,
\qquad
G\colon V\subset Y\to Z,
$$

with $F(U)\subset V$. Suppose $F$ is Fréchet differentiable at $x$ and $G$ is Fréchet differentiable at $F(x)$. Then

$$
(G\circ F)'(x) =
G'(F(x))F'(x).
$$

The formula is the same as in finite-dimensional calculus, but it is worth seeing why the bounded-linear remainder definition produces it.

Write

$$
F(x+h) =
F(x)+Ah+r_{F}(h),
\qquad
A\coloneqq F'(x),
$$

with

$$
\lVert r_{F}(h)\rVert_{Y}
=o(\lVert h\rVert_{X}).
$$

Set

$$
k
\coloneqq
Ah+r_{F}(h).
$$

Because $r_{F}(h)=o(\lVert h\rVert_{X})$, there is a neighborhood of the origin in which

$$
\lVert r_{F}(h)\rVert_{Y}
\leq
\lVert h\rVert_{X}.
$$

For such $h$, boundedness of $A$ gives

$$
\lVert k\rVert_{Y}
\leq
\left(\lVert A\rVert_{\mathcal L(X,Y)}+1\right)
\lVert h\rVert_{X}.
$$

Thus $k$ tends to zero at least proportionally to $h$.

Now expand $G$ around $F(x)$:

$$
G(F(x)+k) =
G(F(x))+Bk+r_{G}(k),
\qquad
B\coloneqq G'(F(x)),
$$

where

$$
\lVert r_{G}(k)\rVert_{Z}
=o(\lVert k\rVert_{Y}) =
o(\lVert h\rVert_{X}).
$$

Substituting $k=Ah+r_{F}(h)$ gives

$$
\begin{aligned}
G(F(x+h))-G(F(x))
&=
B\left(Ah+r_{F}(h)\right)+r_{G}(k)
\\
&=
BAh
+
Br_{F}(h)
+
r_{G}(k).
\end{aligned}
$$

Since $B$ is bounded,

$$
\lVert Br_{F}(h)\rVert_{Z}
=o(\lVert h\rVert_{X}),
$$

and the same is true of $r_{G}(k)$. The first-order linear part is therefore $BAh$, which proves the chain-rule formula.

### An affine coordinate map and a pulled-back covector

A common situation is that the independent coordinates of a variable are not the
same as the coordinates used to represent the corresponding physical field.

Let

$$
z\in\mathbb{R}^{n}
$$

denote independent coordinates and let

$$
y_{\mathrm{phys}}\in\mathbb{R}^{N},
\qquad
n\leq N,
$$

denote physical coordinates. Consider the affine reconstruction map

$$
R\colon\mathbb{R}^{n}\to\mathbb{R}^{N},
\qquad
R(z)\coloneqq Pz+\ell,
$$

where

$$
P\in\mathbb{R}^{N\times n},
\qquad
\ell\in\mathbb{R}^{N}.
$$

Thus

$$
y_{\mathrm{phys}} =
R(z) =
Pz+\ell.
$$

The fixed vector $\ell$ contributes to the physical value but not to its
perturbations. If the independent coordinates change by $\delta z$, then

$$
R'(z)[\delta z] =
P\delta z,
$$

so

$$
\delta y_{\mathrm{phys}} =
P\delta z.
$$

Now let

$$
\Phi_{\mathrm{phys}}\colon\mathbb{R}^{N}\to\mathbb{R}
$$

be a differentiable scalar functional of the physical coordinates, and define its
expression in independent coordinates by

$$
\widehat\Phi\colon\mathbb{R}^{n}\to\mathbb{R},
\qquad
\widehat\Phi(z)
\coloneqq
\Phi_{\mathrm{phys}}(R(z)).
$$

The chain rule gives

$$
\widehat\Phi'(z)[\delta z] =
\Phi_{\mathrm{phys}}'(y_{\mathrm{phys}})
[P\delta z].
$$

Suppose the derivative of $\Phi_{\mathrm{phys}}$ at $y_{\mathrm{phys}}$ has
physical-coordinate covector

$$
r_{\mathrm{phys}}\in\mathbb{R}^{N}.
$$

Then

$$
\begin{aligned}
\widehat\Phi'(z)[\delta z]
&=
r_{\mathrm{phys}}^{\mathsf T}P\delta z
\\
&=
\left(P^{\mathsf T}r_{\mathrm{phys}}\right)^{\mathsf T}
\delta z.
\end{aligned}
$$

Therefore the same derivative, expressed as a covector on the independent coordinate
space, has coordinates

$$
r_{\mathrm{ind}} = P^{\mathsf T} r_{\mathrm{phys}}
$$

The two transformations run in opposite directions:

```text
primal perturbation
    δz  ──────────►  P δz
    independent      physical

covector
    Pᵀ r_phys  ◄───  r_phys
    independent      physical
```

Primal perturbations are pushed forward by $P$. Covectors are pulled back by
$P^{\mathsf T}$. This is simply the chain rule applied to an affine coordinate map;
no inner product or Riesz identification is involved.

## 8. Scalar functionals: derivative versus gradient

Let

$$
J\colon U\subset X\to\mathbb{R}
$$

be Fréchet differentiable. Since the codomain is $\mathbb{R}$, its derivative is a bounded linear functional:

$$
J'(x)\in X^{\ast}.
$$

For a perturbation $h\in X$,

$$
J(x+h) =
J(x)+J'(x)[h]+o(\lVert h\rVert_{X}).
$$

The derivative therefore answers the question:

> what is the first-order change in $J$ if we perturb $x$ by $h$?

Nothing in this definition turns $J'(x)$ into a primal vector.

If $X=H$ is a Hilbert space, a chosen inner product supplies the Riesz map

$$
R_{H}\colon H\to H^{\ast}.
$$

The **gradient of $J$ with respect to that inner product** is the vector

$$
\nabla_{H}J(x)
\in H
$$

satisfying

$$
R_{H}\nabla_{H}J(x) =
J'(x),
$$

or equivalently,

$$
J'(x)[h] =
(\nabla_{H}J(x),h)_{H}
\qquad
\text{for every }h\in H.
$$

The derivative is fixed by the functional. The gradient changes when the inner product changes.

### The same derivative under two inner products

Consider

$$
J(x_{1},x_{2}) =
\frac{1}{2}x_{1}^{2}+2x_{2}^{2}.
$$

Expanding $J(x+h)-J(x)$ gives

$$
\begin{aligned}
J(x+h)-J(x)
&=
x_{1}h_{1}+4x_{2}h_{2}
+
\frac{1}{2}h_{1}^{2}+2h_{2}^{2}.
\end{aligned}
$$

Hence

$$
J'(x)[h] =
x_{1}h_{1}+4x_{2}h_{2}.
$$

Its covector coordinates are

$$
r =
\begin{bmatrix}
x_{1}\\
4x_{2}
\end{bmatrix}.
$$

Under the Euclidean inner product, $G=I$, so

$$
\nabla_{2}J(x) =
\begin{bmatrix}
x_{1}\\
4x_{2}
\end{bmatrix}.
$$

Now choose instead

$$
(x,z)_{G}
\coloneqq
x^{\mathsf T}
\begin{bmatrix}
1&0\\
0&4
\end{bmatrix}
z.
$$

The gradient $g_{G}$ must satisfy

$$
Gg_{G}=r.
$$

Therefore

$$
g_{G} =
\begin{bmatrix}
x_{1}\\
x_{2}
\end{bmatrix}.
$$

Both vectors represent the same derivative under different inner products:

$$
J'(x)[h] =
(\nabla_{2}J(x),h)_{2} =
(g_{G},h)_{G}.
$$

The derivative did not change. Only its primal representative changed.

### A function-space example

Let $H=L^{2}(\Omega)$ and

$$
J(u)
\coloneqq
\frac{1}{2}
\lVert u-z\rVert_{L^{2}(\Omega)}^{2}.
$$

For a perturbation $h$,

$$
\begin{aligned}
J(u+h)-J(u)
&= \frac{1}{2}\int_{\Omega}\left((u-z)+h\right)^{2} -
(u-z)^{2}\thinspace\mathrm{d}x
\\
&=
\int_{\Omega}(u-z)h\thinspace\mathrm{d}x
+
\frac{1}{2}\int_{\Omega}h^{2}\thinspace\mathrm{d}x.
\end{aligned}
$$

Thus

$$
J'(u)[h] =
\int_{\Omega}(u-z)h \thinspace\mathrm{d}x.
$$

With the $L^{2}$ inner product, the Riesz representative is simply

$$
\nabla_{L^{2}}J(u)=u-z.
$$

If the same functional is considered on a different Hilbert space with a different inner product, the derivative formula remains the same functional on admissible directions, but its gradient is obtained by solving the corresponding Riesz problem.

## 9. Derivative actions and adjoint actions

Let

$$
F\colon X\to Y
$$

be differentiable. At a fixed point $x$, its derivative is a bounded linear operator

$$
F'(x)\colon X\to Y.
$$

A perturbation moves forward through this linearization:

$$
h
\longmapsto
F'(x)h.
$$

Given a covector $q\in Y^{\ast}$, we may pull it back through the derivative:

$$
F'(x)^{\ast}q
\in X^{\ast},
$$

where

$$
\langle F'(x)^{\ast}q,h\rangle_{X^{\ast},X} =
\langle q,F'(x)h\rangle_{Y^{\ast},Y}
$$

for every $h\in X$.

This is the adjoint or transpose action of the derivative. It is defined entirely by the dual pairings.

### Returning to the nonlinear example

For

$$
F(x_{1},x_{2}) =
\begin{bmatrix}
x_{1}^{2}+x_{2}\\
\sin x_{1}
\end{bmatrix},
$$

we found

$$
F'(x) =
\begin{bmatrix}
2x_{1}&1\\
\cos x_{1}&0
\end{bmatrix}.
$$

For a perturbation

$$
h=
\begin{bmatrix}
h_{1}\\
h_{2}
\end{bmatrix},
$$

the forward derivative action is

$$
F'(x)h =
\begin{bmatrix}
2x_{1}h_{1}+h_{2}\\
(\cos x_{1})h_{1}
\end{bmatrix}.
$$

For a dual seed

$$
q=
\begin{bmatrix}
q_{1}\\
q_{2}
\end{bmatrix},
$$

the transpose action is

$$
F'(x)^{\ast}q
\quad\longleftrightarrow\quad
F'(x)^{\mathsf T}q =
\begin{bmatrix}
2x_{1}q_{1}+(\cos x_{1})q_{2}\\
q_{1}
\end{bmatrix}.
$$

The defining equality can be checked directly:

$$
\begin{aligned}
q^{\mathsf T}F'(x)h
&=
q_{1}(2x_{1}h_{1}+h_{2})
+q_{2}(\cos x_{1})h_{1}
\\
&=
\left(2x_{1}q_{1}+(\cos x_{1})q_{2}\right)h_{1}
+q_{1}h_{2}
\\
&=
\left(F'(x)^{\mathsf T}q\right)^{\mathsf T}h.
\end{aligned}
$$

### When the residual already takes values in a dual space

Variational PDE residuals often have the form

$$
E\colon X\to Z^{\ast}.
$$

Then

$$
E'(x)\colon X\to Z^{\ast}.
$$

For a test-space element $p\in Z$, the residual covector $E'(x)h\in Z^{\ast}$ can act directly on $p$. This defines a covector on $X$ by

$$
\langle E'(x)^{\ast}p,h\rangle_{X^{\ast},X}
\coloneqq
\langle E'(x)h,p\rangle_{Z^{\ast},Z}.
$$

Thus

$$
E'(x)^{\ast}p
\in X^{\ast}.
$$

This is the transpose construction adapted to a map whose codomain is already a dual space: the defining pairing above is all that is needed for the operation. No optimization metric or Riesz map is needed. A metric becomes relevant only if we later want to convert the resulting covector in $X^{\ast}$ into a primal vector in $X$.

This distinction is one of the main reasons to keep pairings and inner products separate from the beginning.

## 10. Second derivatives and Hessian actions

Suppose a scalar functional

$$
J\colon U\subset X\to\mathbb{R}
$$

is Fréchet differentiable and its derivative map

$$
J'\colon U\to X^{\ast}
$$

is itself Fréchet differentiable. The **second derivative** at $x$ is then

$$
J''(x)
\in
\mathcal L(X,X^{\ast}).
$$

Applied to a direction $h\in X$, it produces a covector

$$
J''(x)h
\in X^{\ast}.
$$

That covector can act on another direction $k\in X$:

$$
J''(x)[h,k]
\coloneqq
\langle J''(x)h,k\rangle_{X^{\ast},X}.
$$

If $J$ is twice continuously Fréchet differentiable in a neighborhood of $x$, this bilinear form is symmetric:

$$
J''(x)[h,k] =
J''(x)[k,h].
$$

For the numerical uses that follow, the important operation is often not “form the Hessian matrix,” but rather

$$
h
\longmapsto
J''(x)h.
$$

This is a **Hessian action**.

### Quadratic example

Let

$$
J(x) =
\frac{1}{2}x^{\mathsf T}Qx-b^{\mathsf T}x,
$$

where $Q$ is symmetric. Expanding $J(x+h)-J(x)$ gives

$$
\begin{aligned}
J(x+h)-J(x)
&=
\frac{1}{2}(x+h)^{\mathsf T}Q(x+h)
-b^{\mathsf T}(x+h)
-\frac{1}{2}x^{\mathsf T}Qx
+b^{\mathsf T}x
\\
&=
(Qx-b)^{\mathsf T}h
+
\frac{1}{2}h^{\mathsf T}Qh.
\end{aligned}
$$

Therefore

$$
J'(x)
\quad\longleftrightarrow\quad
Qx-b,
$$

and differentiating the derivative gives

$$
J''(x)h
\quad\longleftrightarrow\quad
Qh.
$$

The Hessian action does not depend on $x$ because the functional is quadratic. For a nonlinear functional it generally does.

If a Hilbert-space gradient rather than a covector is desired, a Riesz map can again be applied after the second derivative. That additional identification depends on the chosen metric; the second derivative itself does not.

## 11. Taylor remainders and numerical checks

The definition of the Fréchet derivative already gives a first-order remainder:

$$
F(x+h)-F(x)-F'(x)h =
o(\lVert h\rVert_{X}).
$$

Along a fixed direction $d$, set $h=\varepsilon d$. Then

$$
F(x+\varepsilon d) =
F(x)
+
\varepsilon F'(x)d
+
o(\lvert\varepsilon\rvert).
$$

This immediately justifies a finite-difference derivative check:

$$
\frac{F(x+\varepsilon d)-F(x)}{\varepsilon}
\longrightarrow
F'(x)d
$$

as $\varepsilon\to0$.

For a twice Fréchet-differentiable map $F\colon X\to Y$, the second derivative $F''(x)$ is a bounded bilinear map from $X\times X$ to $Y$. If $F$ is twice continuously Fréchet differentiable near $x$, Taylor's theorem strengthens the expansion to

$$
F(x+h) =
F(x)
+
F'(x)h
+
\frac{1}{2}F''(x)[h,h]
+
o(\lVert h\rVert_{X}^{2}).
$$

A proof of the general Banach-space Taylor theorem is outside the scope of this prerequisite chapter; it is the infinite-dimensional analogue of the familiar finite-dimensional theorem. Its consequence for verification is simple and useful.

Define the first-order remainder along $d$ by

$$
R_{1}(\varepsilon)
\coloneqq
F(x+\varepsilon d)
-F(x)
-\varepsilon F'(x)d.
$$

Under the twice-differentiable assumptions,

$$
\lVert R_{1}(\varepsilon)\rVert_{Y} =
O(\varepsilon^{2}).
$$

Here $O(\varepsilon^{2})$ means that there are constants $C\gt 0$ and $\varepsilon_{0}\gt 0$ such that

$$
\lVert R_{1}(\varepsilon)\rVert_{Y}
\leq
C\lvert\varepsilon\rvert^{2}
$$

whenever $\lvert\varepsilon\rvert\lt\varepsilon_{0}$. Thus halving $\varepsilon$ should reduce the remainder by approximately a factor of four until roundoff or other numerical errors dominate.

### An adjoint check

A derivative implementation and its transpose implementation can also be checked independently of finite differences. For arbitrary compatible $h$ and $q$, the defining identity requires

$$
\langle q,F'(x)h\rangle_{Y^{\ast},Y} =
\langle F'(x)^{\ast}q,h\rangle_{X^{\ast},X}.
$$

If, in chosen bases, the derivative $F'(x)$ is represented by a matrix
$\mathbf A(x)$, the same identity becomes

$$
q^{\mathsf T}\mathbf A(x)h =
\left(\mathbf A(x)^{\mathsf T}q\right)^{\mathsf T}h.
$$

If the two sides disagree beyond numerical tolerance, the forward derivative action, the transpose action, the chosen pairings, or the coordinate transformations are inconsistent.

### A second-order check

For a scalar functional, a Hessian action can be compared against the change in the first derivative:

$$
\frac{J'(x+\varepsilon h)-J'(x)}{\varepsilon}
\longrightarrow
J''(x)h
$$

in $X^{\ast}$ as $\varepsilon\to0$.

This is the second-derivative analogue of the first-order finite-difference check. It verifies the action $J''(x)h$ without requiring an assembled Hessian matrix.

## 12. What changes in infinite dimensions

The formulas in this chapter look deliberately similar to finite-dimensional calculus, but several assumptions that are automatic in $\mathbb{R}^{n}$ become meaningful in function spaces.

First, **norms matter**. Fréchet differentiability is defined using the norms of the domain and codomain. The same algebraic formula can define a differentiable map between one pair of function spaces and fail to do so between another pair. Pointwise nonlinearities are a common example: whether products or compositions remain in the desired space depends on integrability and regularity.

Second, **directional derivatives are weaker than a full derivative**. A map may possess derivatives along every direction without those directional derivatives forming one bounded linear approximation that is uniform over small perturbations. This is why the Fréchet definition contains an explicit remainder estimate.

Third, **a dual space is not automatically the same space as the primal one**. Hilbert spaces admit a Riesz identification, but the identification depends on the inner product. In a general Banach space there is no canonical inner-product-based gradient.

Fourth, **adjoints depend on the spaces and pairings**. A formula such as “take the transpose” is incomplete until the domain, codomain, and representations of primal and dual objects are known. In coordinates the transpose matrix is a consequence of those choices, not a substitute for them.

Finally, **operator actions need not imply assembled matrices**. In infinite dimensions there may be no finite matrix at all. After discretization, a matrix may exist conceptually but still be unnecessary to store if the action $F'(x)h$, $F'(x)^{\ast}q$, or $J''(x)h$ can be evaluated directly.

These observations are enough for the main route. General reflexivity theory, bidual spaces, subtle examples separating Gâteaux and Fréchet differentiability, complex Hilbert-space adjoints, and automatic-differentiation machinery are intentionally left to specialist references.

## Where to go next

- **Default continuation for optimization:** [05 · Unconstrained numerical optimization](05-unconstrained-optimization.md) uses the derivative/gradient distinction, Hessian actions, and pullbacks developed here to build concrete optimization methods.
- **If finite-element coordinates are still unfamiliar:** read [03 · Finite elements](03-finite-elements.md) before moving into PDE-constrained optimization.
- **For PDE-constrained optimization:** [07 · PDE-constrained optimization](07-pde-constrained-optimization.md) combines this chapter with [01 · Function spaces and weak PDEs](01-function-spaces-and-weak-pdes.md) to derive state sensitivities, adjoint equations, and reduced derivatives.
- **Into the `nmopt` manual:** the same distinctions are specialized to project contracts in:
  - [manual 02 · Spaces, coordinates, and duality](../concepts/02-spaces-coordinates-and-duality.md), for physical versus independent coordinates and primal/dual runtime roles;
  - [manual 03 · Operators, derivatives, and adjoints](../concepts/03-operators-derivatives-and-adjoints.md), for residual JVP/VJP actions and transpose consistency;
  - [manual 04 · Metrics, gradients, and constraints](../concepts/04-metrics-gradients-and-constraints.md), for discrete Riesz maps and metric gradients;
  - [manual 05 · Reduced state–adjoint formulation](../concepts/05-reduced-state-adjoint-formulation.md), for the project reduced-derivative construction.

## References and further reading

Publication details, stable links, and access notes are collected in the [Background reference catalogue](references.md).

- **Stefano Mancini, Andrea Benvenuti, and Luca Heltai, *Numerical Methods for Partial Differential Equations*.**
  - §1.1.5, **Dual spaces**, for bounded linear functionals and dual norms.
  - §1.1.6, **Riesz representation theorem**, for the Hilbert-space primal/dual identification.
  - §1.1.7, **Bilinear forms and operators**, for the operator view of variational forms.
- **Manzoni, Quarteroni, and Salsa, *Optimal Control of Partial Differential Equations*.**
  - Appendix A, §A.2.2, **Functionals, Dual space and Riesz Theorem**.
  - Appendix A, §A.2.4, **Adjoint Operators**.
  - Appendix A, §A.7.1, **The Fréchet Derivative**, including the chain rule and derivative/gradient distinction.
  - Appendix A, §A.7.2, **The Gâteaux Derivative**, for the weaker directional notion and its relation to Fréchet differentiability.
  - Appendix A, §A.7.4, **Second Order Derivatives**, for Hessian actions and the Banach-space Taylor expansion.

The bibliographies collected in Appendices A.2 and A.7 of Manzoni–Quarteroni–Salsa provide routes into deeper functional analysis and Banach-space differential calculus.
