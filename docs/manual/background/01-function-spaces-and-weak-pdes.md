# 01 · Function spaces and weak PDEs

**Background navigation:** [Index](README.md) \
Previous: [Index](README.md) \
Next: [02 · Duality, derivatives, and adjoints](02-duality-derivatives-and-adjoints.md) \
Companion: [01a · Further notes on function spaces and weak PDEs](01a-further-function-space-notes.md)

## Purpose

Many PDEs are first written as pointwise differential equations. For example,

$$
-\Delta y = f
\qquad \text{in }\Omega,
$$

with a boundary condition such as

$$
y=0
\qquad \text{on }\partial\Omega.
$$

That notation is compact, but it hides several assumptions. What does it mean for the second derivatives in $\Delta y$ to exist? In what sense is the boundary value imposed? What class of functions contains the solution? If a numerical method replaces the infinite-dimensional problem by a finite-dimensional one, what exactly is being approximated?

Weak, or variational, formulations answer these questions by changing the point of view. Instead of requiring the differential equation to hold pointwise, they ask the solution to satisfy an integral identity against a family of test functions. This reduces the differentiability demanded of the solution, makes the function spaces part of the mathematical statement, and leads directly to the Galerkin viewpoint used by finite-element methods.

This chapter develops the minimum functional-analytic language needed to read such formulations with understanding. The main example is the Poisson equation, because it exposes the essential ideas without hiding them behind a complicated PDE. The same ideas recur throughout elliptic PDEs and later in PDE-constrained optimization.

The goal is not to build functional analysis from first principles. We will use standard results from Lebesgue and Sobolev-space theory, derive the weak-form calculations that carry the main conceptual load, state rather than prove the trace and Lax–Milgram theorems, and postpone the deeper study of duality, Riesz maps, derivatives, and adjoints to the next chapter. A companion note, [01a · Further notes on function spaces and weak PDEs](01a-further-function-space-notes.md), develops several ideas that are useful later but would interrupt the main path here.

## Before you start

This chapter assumes ordinary multivariable calculus, integration by parts, basic linear algebra, and the usual notation for partial derivatives and the Laplacian. Some familiarity with differential equations is useful, but no prior course in functional analysis or finite elements is assumed.

## What you will be able to do

After this chapter, you should be able to:

- explain why classical differentiability is often too restrictive for PDE solutions;
- interpret $L^{2}(\Omega)$, $H^{1}(\Omega)$, and $H_{0}^{1}(\Omega)$ operationally;
- define a weak derivative from integration by parts;
- explain why boundary values of Sobolev functions require a trace operator;
- derive the weak Poisson problem from its strong form;
- distinguish essential and natural boundary conditions;
- interpret a weak residual as a bounded linear functional on a test space;
- state the Lax–Milgram theorem and check its hypotheses for the model Poisson problem at a working level.

## Roadmap

We begin with the difficulty in asking for classical solutions and introduce just enough normed-space language to discuss convergence and stability. We then move through $L^{2}(\Omega)$, weak derivatives, Sobolev spaces, boundary traces, and linear functionals. With those pieces in place, we derive the weak Poisson equation, rewrite it as an abstract variational problem, and use Lax–Milgram to understand existence, uniqueness, and stability. A mixed-boundary example shows why Dirichlet and Neumann conditions enter a weak formulation in different ways.

The finite-element method itself is deliberately deferred to [03 · Finite elements](03-finite-elements.md). The present chapter stops once the continuous variational problem is mathematically clear. When a trace theorem, compactness result, Sobolev embedding, or regularity statement is mentioned but not developed, the companion [01a · Further notes on function spaces and weak PDEs](01a-further-function-space-notes.md) provide a first orientation and precise references.

## 1. Why a classical solution may ask for too much

Throughout this chapter, a **domain** $\Omega\subset\mathbb R^{d}$ means an open connected set, and $\partial\Omega$ denotes its boundary.

Consider the one-dimensional boundary-value problem

$$
-y'' = f
\qquad \text{in }(0,1),
$$

with $y(0)=y(1)=0$. If we interpret this equation classically, then $y''(x)$ must exist at every interior point where the equation is asserted.

In several dimensions, the Poisson equation

$$
-\Delta y=f
$$

likewise appears to require second partial derivatives of $y$, since the Laplacian is

$$
\Delta y
:=
\sum_{i=1}^{d}
\frac{\partial^{2}y}{\partial x_{i}^{2}}.
$$

For very smooth data and domains, sufficiently smooth solutions may exist. But numerical PDEs should not be built around that fortunate case. Coefficients can jump, forcing terms can be rough, domains can have corners, and even a perfectly meaningful physical solution may fail to possess all derivatives required by the strong equation in the classical sense.

A simple one-dimensional example already shows why a weaker notion of derivative is useful. The function

$$
y(x)=|x|
\qquad \text{for }x\in(-1,1)
$$

is continuous but not classically differentiable at $x=0$. Away from the origin its derivative is

$$
y'(x)=
\begin{cases}
-1, & x<0,\\
1, & x>0.
\end{cases}
$$

The failure at a single point does not prevent this piecewise-constant function from describing the derivative in an integral sense. Section 4 will make that statement precise and verify it directly.

There is a second reason to move away from pointwise equations. Numerical methods such as finite elements do not normally enforce the PDE independently at every point of the domain. They approximate a **weak**, or **variational**, problem posed in a function space. Here *variational* means that the equation is characterized by how it acts against admissible test functions, which represent possible variations or directions in which we probe the equation. It does not, by itself, mean that we are solving an optimization problem.

To understand the discrete method later, we therefore need to know what that continuous variational problem is and why its function spaces are chosen as they are.

## 2. Normed, Banach, and Hilbert spaces

The word *space* in PDE analysis means more than a collection of functions. We need a way to measure the size of a function and to say when a sequence of functions converges.

A **normed vector space** $X$ is a vector space equipped with a norm $\lVert\cdot\rVert_{X}$. The norm satisfies the familiar properties

$$
\lVert x\rVert_{X}\ge 0,
\qquad
\lVert x\rVert_{X}=0 \Longleftrightarrow x=0,
$$

$$
\lVert \alpha x\rVert_{X}=|\alpha|\lVert x\rVert_{X},
\qquad
\lVert x+z\rVert_{X}\le \lVert x\rVert_{X}+\lVert z\rVert_{X}.
$$

The norm lets us define convergence: $x_{n}\to x$ in $X$ when

$$
\lVert x_{n}-x\rVert_{X}\to 0.
$$

A sequence $\{x_{n}\}$ is **Cauchy** if its elements eventually become arbitrarily close to one another, that is,

$$
\lVert x_{n}-x_{m}\rVert_{X}\to0
\qquad
\text{as }n,m\to\infty.
$$

A **Banach space** is a normed space in which every Cauchy sequence converges to an element of the space. Completeness matters because analytical arguments often construct a solution as the limit of approximations. If the space were incomplete, that limiting object could fall outside the space in which the problem was posed.

A **Hilbert space** is a complete inner-product space. Its norm comes from an inner product,

$$
\lVert x\rVert_{H}^{2}=(x,x)_{H}.
$$

The distinction will become important later because inner products give Hilbert spaces additional structure: orthogonality, projections, and the Riesz representation theorem. For this chapter, the most important examples are $L^{2}(\Omega)$ and $H^{1}(\Omega)$, both Hilbert spaces.

We do not need the general theory of Banach spaces before continuing. The point is simply that a PDE is posed in a space whose norm encodes what kind of functions are admissible and what notion of convergence and stability we care about.

## 3. The space $L^{2}(\Omega)$

Let $\Omega\subset\mathbb R^{d}$ be a bounded domain. We will not construct Lebesgue measure or the Lebesgue integral here; that belongs to a course in real analysis. We only need the resulting notion of integrability. All ordinary continuous and piecewise-continuous functions used in the examples below are Lebesgue measurable, so this abstraction does not obstruct the calculations we perform.

Informally, $L^{2}(\Omega)$ contains measurable functions whose square is integrable:

$$
L^{2}(\Omega)
=
\left\{
q:\Omega\to\mathbb R:
\int_{\Omega}|q(x)|^{2}\mathrm{d}x<\infty
\right\},
$$

where functions that agree **almost everywhere** are identified. A property holds almost everywhere, abbreviated a.e., if it can fail only on a set of Lebesgue measure zero. Thus two functions that differ at a single point, or on any other measure-zero set, represent the same element of $L^{2}(\Omega)$.

The norm and inner product are

$$
\lVert q\rVert_{L^{2}(\Omega)}
:=
\left(
\int_{\Omega}|q|^{2}\mathrm{d}x
\right)^{1/2},
$$

and

$$
(q,r)_{L^{2}(\Omega)}
:=
\int_{\Omega}qr\mathrm{d}x.
$$

With this inner product, $L^{2}(\Omega)$ is a Hilbert space.

Two facts are especially important for numerical PDEs.

First, an $L^{2}$ function need not be continuous. Integrability is the requirement; pointwise smoothness is not. For example,

$$
q(x)=
\begin{cases}
1, & 0<x<1/2,\\
0, & 1/2<x<1
\end{cases}
$$

belongs to $L^{2}(0,1)$ even though it jumps at $x=1/2$.

Second, point values are not generally intrinsic pieces of information for an $L^{2}$ element. If we redefine the preceding function at $x=1/2$, its $L^{2}$ norm, inner products, and equivalence class do not change. This is why a notation such as $q(x)$ should not trick us into assuming that every function-space element carries meaningful pointwise data.

The Cauchy–Schwarz inequality in $L^{2}$,

$$
\left|
\int_{\Omega}qr\mathrm{d}x
\right|
\le
\lVert q\rVert_{L^{2}(\Omega)}
\lVert r\rVert_{L^{2}(\Omega)},
$$

will be one of our basic tools. In particular, it immediately tells us when products such as $qr$ are integrable and when an integral expression depends continuously on one of its arguments.

## 4. Weak derivatives from integration by parts

Suppose for the moment that $y$ and a test function $\varphi$ are smooth on an interval and that $\varphi$ vanishes at the boundary. Integration by parts gives

$$
\int y'\varphi\mathrm{d}x
=
-
\int y\varphi'\mathrm{d}x.
$$

Notice that the right-hand side only needs $y$ itself, not its derivative. This suggests turning the identity around and using it as the *definition* of a derivative for functions that may not possess a classical one.

Before doing so, two pieces of notation are worth making explicit.

The space $L_{\mathrm{loc}}^{1}(\Omega)$ consists of functions that are integrable on every compact subset of $\Omega$. In Euclidean space, you may read this as saying that the function is integrable on every closed bounded region that stays strictly inside the domain. We ask only for local integrability because test functions will vanish outside such an interior region.

The notation $C_{c}^{\infty}(\Omega)$ means the set of infinitely differentiable functions with **compact support** inside $\Omega$: each test function is exactly zero outside some compact subset of the domain. Many PDE texts write the same space as $C_{0}^{\infty}(\Omega)$. The compact-support condition is what removes boundary terms in the definition below.

Let $y\in L_{\mathrm{loc}}^{1}(\Omega)$. A function $w\in L_{\mathrm{loc}}^{1}(\Omega)$ is the **weak partial derivative** of $y$ with respect to $x_{i}$ if

$$
\int_{\Omega}w\varphi\mathrm{d}x
=
-
\int_{\Omega}y
\frac{\partial\varphi}{\partial x_{i}}
\mathrm{d}x
\qquad
\text{for every }\varphi\in C_{c}^{\infty}(\Omega).
$$

If $y$ is classically differentiable, ordinary integration by parts shows that its classical derivative satisfies this identity. The weak definition therefore extends the classical derivative rather than replacing it with a conflicting notion on smooth functions.

### Checking the derivative of $|x|$

Return to

$$
y(x)=|x|
\qquad
\text{on }(-1,1),
$$

and define

$$
w(x)=
\begin{cases}
-1, & x<0,\\
1, & x>0.
\end{cases}
$$

The value assigned to $w(0)$ is irrelevant because a single point has measure zero. Let $\varphi\in C_{c}^{\infty}(-1,1)$. Then

```math
\begin{aligned}
-
\int_{-1}^{1}|x|\varphi'(x)\mathrm{d}x
&=
\int_{-1}^{0}x\varphi'(x)\mathrm{d}x
-
\int_{0}^{1}x\varphi'(x)\mathrm{d}x
\\
&=
\left[x\varphi(x)\right]_{-1}^{0}
-
\int_{-1}^{0}\varphi(x)\mathrm{d}x
-
\left[x\varphi(x)\right]_{0}^{1}
+
\int_{0}^{1}\varphi(x)\mathrm{d}x.
\end{aligned}
```

Because $\varphi$ has compact support in $(-1,1)$, it vanishes near $-1$ and $1$. The factors of $x$ also make the two contributions at $x=0$ vanish. Hence

$$
-
\int_{-1}^{1}|x|\varphi'(x)\mathrm{d}x
=
-
\int_{-1}^{0}\varphi(x)\mathrm{d}x
+
\int_{0}^{1}\varphi(x)\mathrm{d}x.
$$

But the right-hand side is exactly

$$
\int_{-1}^{1}w(x)\varphi(x)\mathrm{d}x.
$$

Therefore $w$ is the weak derivative of $|x|$. We did not ignore the failure of classical differentiability at the origin; the integral identity simply shows that this isolated failure does not prevent an $L^{1}_{\mathrm{loc}}$ derivative from existing.

Weak derivatives are defined only up to almost-everywhere equality, just like elements of $L^{p}$ spaces. We can therefore ask not only whether weak derivatives exist, but also whether they have a specified integrability. That leads to Sobolev spaces.

The general theory of distributions pushes this idea further and allows derivatives of objects more singular than locally integrable functions. A short orientation appears in [01a · Further notes on function spaces and weak PDEs](01a-further-function-space-notes.md); distribution theory itself is not needed for the main path of this chapter.

## 5. The Sobolev space $H^{1}(\Omega)$

The space $H^{1}(\Omega)$ consists of square-integrable functions whose first weak derivatives are also square-integrable:

$$
H^{1}(\Omega)
:=
\left\{
y\in L^{2}(\Omega):
\frac{\partial y}{\partial x_{i}}
\in L^{2}(\Omega)
\text{ weakly for }i=1,\ldots,d
\right\}.
$$

The weak derivatives form the **weak gradient**

$$
\nabla y
=
\left(
\frac{\partial y}{\partial x_{1}},
\ldots,
\frac{\partial y}{\partial x_{d}}
\right).
$$

The standard $H^{1}$ inner product is

$$
(y,v)_{H^{1}(\Omega)}
:=
\int_{\Omega}yv\mathrm{d}x
+
\int_{\Omega}\nabla y\cdot\nabla v\mathrm{d}x,
$$

and the induced norm is

$$
\lVert y\rVert_{H^{1}(\Omega)}^{2}
=
\lVert y\rVert_{L^{2}(\Omega)}^{2}
+
\lVert\nabla y\rVert_{L^{2}(\Omega)^{d}}^{2}.
$$

With this inner product, $H^{1}(\Omega)$ is a Hilbert space. The completeness statement is a theorem from Sobolev-space theory; we use it here rather than prove it.

The definition is worth reading operationally. Membership in $H^{1}(\Omega)$ says that both the field and its first weak derivatives have finite $L^{2}$ size. It does **not** say that $y$ is continuously differentiable, or even that it has a meaningful value at every point. The topology of the space is tied to integral control of the function and its weak gradient.

This is precisely the regularity that the weak Poisson form will need: its main integral contains $\nabla y$, but no second derivative of $y$. The earlier example $y(x)=|x|$ belongs to $H^{1}(-1,1)$ because both $|x|$ and its weak derivative are square-integrable, even though the classical derivative fails at the origin.

More generally, Sobolev spaces $W^{k,p}(\Omega)$ contain functions whose weak derivatives through order $k$ belong to $L^{p}(\Omega)$. When $p=2$, the usual notation is $H^{k}(\Omega)$. We will need almost none of that generality in the main chapter. The companion [01a · Further notes on function spaces and weak PDEs](01a-further-function-space-notes.md) make the higher-dimensional notation and Sobolev seminorms explicit before continuing to embeddings and higher regularity, because all of these ideas recur in finite-element estimates and more delicate PDE examples.

## 6. Boundary values, traces, and $H_{0}^{1}(\Omega)$

At first sight, imposing

$$
y=0
\qquad \text{on }\partial\Omega
$$

looks straightforward. For a continuous function we simply restrict $y$ to the boundary. For an arbitrary $H^{1}(\Omega)$ element, however, pointwise boundary values are not part of the definition, so we need a boundary-value operation that is compatible with the $H^{1}$ topology.

For the rest of the chapter we assume that $\Omega$ is a **bounded Lipschitz domain**. Roughly, Lipschitz regularity means that near every boundary point, after rotating coordinates if necessary, the boundary can be represented as the graph of a function whose slope is uniformly bounded. We do not need the precise geometric definition. We need the consequences: this regularity is sufficient for the trace theorem and the integration-by-parts formulas used below. Boundedness will also be used in the Poincaré inequality.

The **trace theorem** states that there is a bounded linear map

$$
\gamma:H^{1}(\Omega)\to H^{1/2}(\partial\Omega)
$$

that agrees with ordinary boundary restriction for smooth functions. In particular, there is a constant $C_{\mathrm{tr}}>0$ such that

$$
\lVert\gamma y\rVert_{L^{2}(\partial\Omega)}
\le
C_{\mathrm{tr}}
\lVert y\rVert_{H^{1}(\Omega)}.
$$

The fractional space $H^{1/2}(\partial\Omega)$ is the natural trace space of $H^{1}(\Omega)$. The main chapter does not develop fractional Sobolev spaces, but one concrete characterization already makes the notation useful:

$$
H^{1/2}(\partial\Omega)
=
\gamma\bigl(H^{1}(\Omega)\bigr),
$$

with a norm measuring the smallest $H^{1}$ norm among all interior functions having the prescribed trace. The companion [01a · Further notes on function spaces and weak PDEs](01a-further-function-space-notes.md) develops this point and explains why $H^{1/2}$ matters for boundary control and nonhomogeneous Dirichlet data.

The space with homogeneous Dirichlet boundary values is

$$
H_{0}^{1}(\Omega)
:=
\overline{C_{c}^{\infty}(\Omega)}^{H^{1}(\Omega)}.
$$

On a bounded Lipschitz domain, this is equivalently the kernel of the trace operator:

$$
H_{0}^{1}(\Omega)
=
\left\{
y\in H^{1}(\Omega):\gamma y=0
\right\}.
$$

The closure definition is more fundamental than the informal phrase “functions in $H^{1}$ that vanish on the boundary,” because it does not assume ordinary pointwise boundary values.

A second crucial result is the **Poincaré inequality**: there is a constant $C_{P}>0$ such that

$$
\lVert y\rVert_{L^{2}(\Omega)}
\le
C_{P}
\lVert\nabla y\rVert_{L^{2}(\Omega)^{d}}
\qquad
\text{for every }y\in H_{0}^{1}(\Omega).
$$

Why does the zero boundary condition matter? On the full space $H^{1}(\Omega)$, the gradient alone cannot be a norm because every nonzero constant $c$ satisfies

$$
\nabla c=0
$$

while $c\ne0$. Homogeneous Dirichlet data remove those constant modes. Poincaré then tells us that the gradient controls the entire $H^{1}$ size:

```math
\begin{aligned}
\lVert y\rVert_{H^{1}(\Omega)}^{2}
&=
\lVert y\rVert_{L^{2}(\Omega)}^{2}
+
\lVert\nabla y\rVert_{L^{2}(\Omega)^{d}}^{2}
\\
&\le
\left(C_{P}^{2}+1\right)
\lVert\nabla y\rVert_{L^{2}(\Omega)^{d}}^{2}.
\end{aligned}
```

Therefore

$$
\lVert y\rVert_{V}
:=
\lVert\nabla y\rVert_{L^{2}(\Omega)^{d}}
$$

is a genuine norm on $H_{0}^{1}(\Omega)$, equivalent to the standard $H^{1}$ norm. This is exactly the norm in which the Poisson bilinear form will become coercive.

Nonhomogeneous Dirichlet data fit the same structure. If we want

$$
\gamma y=g
$$

and $g\in H^{1/2}(\partial\Omega)$, the trace theorem guarantees at least one **lifting** $\widetilde g\in H^{1}(\Omega)$ with

$$
\gamma\widetilde g=g.
$$

Writing

$$
y=z+\widetilde g
$$

reduces the unknown part to $z\in H_{0}^{1}(\Omega)$. This simple change of variables will reappear later whenever essential boundary data are separated from independent unknowns.

## 7. Linear functionals and a first look at dual spaces

A **linear functional** on a vector space $V$ is a linear map

$$
F:V\to\mathbb R.
$$

If $V$ is normed, a functional is **continuous** when convergence of its input implies convergence of its output:

$$
v_{n}\to v\text{ in }V
\quad\Longrightarrow\quad
F(v_{n})\to F(v).
$$

For a linear functional, continuity is equivalent to **boundedness**: there exists a constant $C\ge0$ such that

$$
|F(v)|
\le
C\lVert v\rVert_{V}
\qquad
\text{for every }v\in V.
$$

We will use this equivalence without proving it. The smallest admissible bound is encoded by the dual norm

$$
\lVert F\rVert_{V^{\ast}}
:=
\sup_{v\ne0}
\frac{|F(v)|}{\lVert v\rVert_{V}}.
$$

The collection of all bounded linear functionals on $V$ is the **dual space**, denoted $V^{\ast}$. The action of $F\in V^{\ast}$ on $v\in V$ is often written as the dual pairing

$$
\langle F,v\rangle_{V^{\ast},V}.
$$

For the moment, the distinction between an element of $V$ and an element of $V^{\ast}$ is more important than the full structure of the dual space. The next background chapter develops that distinction, the Riesz map, and adjoint operators carefully.

A familiar integral gives a useful example. Let

$$
V=H_{0}^{1}(\Omega)
$$

with the gradient norm, and let $f\in L^{2}(\Omega)$. Define

$$
F(v)
:=
\int_{\Omega}fv\mathrm{d}x.
$$

This map is linear. By Cauchy–Schwarz and Poincaré,

```math
\begin{aligned}
|F(v)|
&\le
\lVert f\rVert_{L^{2}(\Omega)}
\lVert v\rVert_{L^{2}(\Omega)}
\\
&\le
C_{P}
\lVert f\rVert_{L^{2}(\Omega)}
\lVert v\rVert_{V}.
\end{aligned}
```

So $F$ is bounded and hence continuous. In particular,

$$
\lVert F\rVert_{V^{\ast}}
\le
C_{P}\lVert f\rVert_{L^{2}(\Omega)}.
$$

This is the precise sense in which an $L^{2}$ forcing term naturally defines a functional acting on the test space.

The dual of $H_{0}^{1}(\Omega)$ is commonly denoted

$$
H^{-1}(\Omega)
:=
\left(H_{0}^{1}(\Omega)\right)^{\ast}.
$$

The notation does not mean “ordinary functions with minus one derivative.” It records a dual-space relationship. The important idea for now is that a weak PDE residual is naturally an object that **acts on test functions**, and therefore lives in a dual space.

## 8. Bilinear forms and abstract variational problems

Suppose $V$ is a Hilbert space. A **bilinear form** is a map

$$
a:V\times V\to\mathbb R
$$

that is linear in each argument separately. It is **bounded**, also called continuous, if there is a constant $M>0$ such that

$$
|a(y,v)|
\le
M\lVert y\rVert_{V}\lVert v\rVert_{V}
\qquad
\text{for all }y,v\in V.
$$

It is **coercive** if there is a constant $\alpha>0$ such that

$$
a(v,v)
\ge
\alpha\lVert v\rVert_{V}^{2}
\qquad
\text{for all }v\in V.
$$

Given such a form and a functional $F\in V^{\ast}$, an **abstract variational problem** has the form

> Find $y\in V$ such that
> $$
> a(y,v)=\langle F,v\rangle_{V^{\ast},V}
> \qquad
> \text{for every }v\in V.
> $$

The word *variational* here refers to the fact that the equation is tested against arbitrary admissible variations $v$. It does not require the problem to have arisen from minimizing an energy, although symmetric coercive problems often admit such an equivalent minimization interpretation.

This format separates the structure of the problem from the particular PDE that produced it. The space $V$ says which trial and test functions are admissible. The form $a$ says how the unknown interacts with a test function. The functional $F$ contains the forcing data.

In the Poisson example the unknown and test functions use the same space. More general formulations may use different trial and test spaces or several coupled spaces. Those cases lead to Petrov–Galerkin and mixed formulations and require additional well-posedness tools such as inf-sup conditions. The companion [01a · Further notes on function spaces and weak PDEs](01a-further-function-space-notes.md) give a brief orientation; the main chapter stays with the coercive case.

The same problem can also be written as a residual equation. Define $E(y)\in V^{\ast}$ by

$$
\langle E(y),v\rangle_{V^{\ast},V}
:=
a(y,v)-\langle F,v\rangle_{V^{\ast},V}.
$$

Then

$$
E(y)=0
\qquad \text{in }V^{\ast}
$$

means exactly that

$$
\langle E(y),v\rangle_{V^{\ast},V}=0
\qquad
\text{for every }v\in V.
$$

This is the viewpoint that later transfers naturally to PDE-constrained optimization: a weak residual is characterized by its action on test functions, not fundamentally by a list of pointwise PDE errors.

## 9. From strong Poisson to weak Poisson

We now have the pieces needed to derive the model weak problem rather than simply write it down.

Let $\Omega\subset\mathbb R^{d}$ be a bounded Lipschitz domain and consider

```math
\begin{aligned}
-\Delta y &= f
&& \text{in }\Omega,\\
y &= 0
&& \text{on }\partial\Omega.
\end{aligned}
```

Assume first that $y$ and $f$ are smooth enough for ordinary calculus identities to hold. Choose a smooth test function $v$ whose boundary value is zero. Multiply the PDE by $v$ and integrate:

$$
-
\int_{\Omega}(\Delta y)v\mathrm{d}x
=
\int_{\Omega}fv\mathrm{d}x.
$$

For smooth $y$, define the outward normal derivative by

$$
\frac{\partial y}{\partial n}
:=
\nabla y\cdot n,
$$

where $n$ is the outward unit normal on $\partial\Omega$. Green's identity gives

$$
-
\int_{\Omega}(\Delta y)v\mathrm{d}x
=
\int_{\Omega}\nabla y\cdot\nabla v\mathrm{d}x
-
\int_{\partial\Omega}
\frac{\partial y}{\partial n}v
\mathrm{d}s.
$$

Because $v$ vanishes on $\partial\Omega$, the boundary term disappears. We obtain

$$
\int_{\Omega}\nabla y\cdot\nabla v\mathrm{d}x
=
\int_{\Omega}fv\mathrm{d}x.
$$

The key change is now visible: the strong equation involved second derivatives of $y$, while the tested identity contains only first derivatives.

> **Why $H_{0}^{1}(\Omega)$ is the right space here.** The choice is dictated by the terms we have just derived, not by notation. If $y,v\in H_{0}^{1}(\Omega)$, then their weak gradients belong to $L^{2}(\Omega)^{d}$, so Cauchy–Schwarz gives
> $$
> \left|
> \int_{\Omega}\nabla y\cdot\nabla v\mathrm{d}x
> \right|
> \le
> \lVert\nabla y\rVert_{L^{2}(\Omega)^{d}}
> \lVert\nabla v\rVert_{L^{2}(\Omega)^{d}}.
> $$
> Thus the left-hand side is finite. If $f\in L^{2}(\Omega)$, then
> ```math
> \begin{aligned}
> \left|
> \int_{\Omega}fv\mathrm{d}x
> \right|
> &\le
> \lVert f\rVert_{L^{2}(\Omega)}
> \lVert v\rVert_{L^{2}(\Omega)}
> \\
> &\le
> C_{P}
> \lVert f\rVert_{L^{2}(\Omega)}
> \lVert\nabla v\rVert_{L^{2}(\Omega)^{d}},
> \end{aligned}
> ```
> so the right-hand side is also finite and continuous in the test function. Finally, $H_{0}^{1}(\Omega)$ encodes the homogeneous Dirichlet condition through its zero trace. The smaller space $H^{2}(\Omega)$ would demand second weak derivatives that no longer appear in the variational identity; the larger space $L^{2}(\Omega)$ would not provide the $L^{2}$ gradients required by the left-hand side. In this sense, $H_{0}^{1}(\Omega)$ is the natural **energy space** for this problem: the quantity
> $$
> a(v,v)=\int_{\Omega}|\nabla v|^{2}\mathrm{d}x
> $$
> is finite there and, by Poincaré, controls the full $H^{1}$ size of $v$.

We can therefore enlarge the class of admissible solutions from classically twice-differentiable functions to

$$
V:=H_{0}^{1}(\Omega).
$$

The **weak Poisson problem** is:

> Given a forcing functional $F\in V^{\ast}$, find $y\in V$ such that
> $$
> \int_{\Omega}\nabla y\cdot\nabla v\mathrm{d}x
> =
> \langle F,v\rangle_{V^{\ast},V}
> \qquad
> \text{for every }v\in V.
> $$

When $f\in L^{2}(\Omega)$, the functional is

$$
\langle F,v\rangle_{V^{\ast},V}
=
\int_{\Omega}fv\mathrm{d}x.
$$

The bilinear form is

$$
a(y,v)
:=
\int_{\Omega}\nabla y\cdot\nabla v\mathrm{d}x,
$$

and the residual is

$$
\langle E(y),v\rangle_{V^{\ast},V}
:=
\int_{\Omega}\nabla y\cdot\nabla v\mathrm{d}x
-
\langle F,v\rangle_{V^{\ast},V}.
$$

The weak state equation is simply

$$
E(y)=0
\qquad \text{in }V^{\ast}.
$$

This is a mathematical problem in its own right, not merely a numerical approximation of the strong equation. A finite-element method will later approximate this variational problem by restricting the trial and test spaces to finite-dimensional subspaces.

## 10. Essential and natural boundary conditions

The homogeneous Dirichlet condition disappeared from the boundary integral in Section 9 because we deliberately chose trial and test functions with zero trace. It did not disappear from the problem: it is built into the admissible space. This is why Dirichlet conditions are called **essential boundary conditions** in variational methods.

Neumann conditions enter differently. Split the boundary into two portions, $\Gamma_{D}$ and $\Gamma_{N}$, and consider

```math
\begin{aligned}
-\Delta y &= f
&& \text{in }\Omega,\\
y &= 0
&& \text{on }\Gamma_{D},\\
\frac{\partial y}{\partial n} &= g
&& \text{on }\Gamma_{N}.
\end{aligned}
```

Use the space

$$
V
:=
\left\{
v\in H^{1}(\Omega):
\gamma v=0\text{ on }\Gamma_{D}
\right\}.
$$

For the standard coercive mixed problem, assume that the Dirichlet portion $\Gamma_{D}$ is large enough—for example, that it has positive boundary measure—so that a Poincaré-type inequality holds on $V$. If $\Gamma_{D}$ is empty, the problem is purely Neumann and requires the separate compatibility discussion below.

For smooth functions, Green's identity gives

```math
\begin{aligned}
\int_{\Omega}\nabla y\cdot\nabla v\mathrm{d}x
-
\int_{\Gamma_{D}}
\frac{\partial y}{\partial n}v\mathrm{d}s
-
\int_{\Gamma_{N}}
\frac{\partial y}{\partial n}v\mathrm{d}s
=
\int_{\Omega}fv\mathrm{d}x.
\end{aligned}
```

The trace of $v$ vanishes on $\Gamma_{D}$, so the first boundary term is zero. Substituting the prescribed Neumann datum on $\Gamma_{N}$ gives

$$
\int_{\Omega}\nabla y\cdot\nabla v\mathrm{d}x
=
\int_{\Omega}fv\mathrm{d}x
+
\int_{\Gamma_{N}}gv\mathrm{d}s.
$$

Thus the Dirichlet condition enters through the definition of $V$, whereas the Neumann condition appears as part of the forcing functional. For this reason Neumann data are called **natural boundary conditions** for this variational form.

The trace theorem makes the boundary functional precise. If $g\in L^{2}(\Gamma_{N})$, then

```math
\begin{aligned}
\left|
\int_{\Gamma_{N}}gv\mathrm{d}s
\right|
&\le
\lVert g\rVert_{L^{2}(\Gamma_{N})}
\lVert\gamma v\rVert_{L^{2}(\Gamma_{N})}
\\
&\le
C_{\mathrm{tr}}
\lVert g\rVert_{L^{2}(\Gamma_{N})}
\lVert v\rVert_{H^{1}(\Omega)}.
\end{aligned}
```

So the boundary integral is not merely formal: under these assumptions it defines a bounded functional on the test space.

### Robin conditions

A Robin condition mixes the value of the state with its normal derivative. Suppose, for simplicity, that on a boundary portion $\Gamma_{R}$ we prescribe

$$
\frac{\partial y}{\partial n}
+
\kappa y
=
g,
$$

where $\kappa\in L^{\infty}(\Gamma_{R})$. Equivalently,

$$
\frac{\partial y}{\partial n}
=
g-\kappa y.
$$

Before inserting the condition, integration by parts gives

$$
\int_{\Omega}\nabla y\cdot\nabla v\mathrm{d}x
-
\int_{\Gamma_{R}}
\frac{\partial y}{\partial n}v\mathrm{d}s
=
\int_{\Omega}fv\mathrm{d}x,
$$

where any essential boundary portion has already been absorbed into the test space. Substitute the Robin condition:

```math
\begin{aligned}
\int_{\Omega}\nabla y\cdot\nabla v\mathrm{d}x
-
\int_{\Gamma_{R}}(g-\kappa y)v\mathrm{d}s
&=
\int_{\Omega}fv\mathrm{d}x.
\end{aligned}
```

Rearranging gives the weak form

$$
\int_{\Omega}\nabla y\cdot\nabla v\mathrm{d}x
+
\int_{\Gamma_{R}}\kappa yv\mathrm{d}s
=
\int_{\Omega}fv\mathrm{d}x
+
\int_{\Gamma_{R}}gv\mathrm{d}s.
$$

The new boundary term belongs on the left because it depends bilinearly on the unknown $y$ and the test function $v$. The trace estimate explains why the assumption $\kappa\in L^{\infty}(\Gamma_{R})$ is convenient:

```math
\begin{aligned}
\left|
\int_{\Gamma_{R}}\kappa yv\mathrm{d}s
\right|
&\le
\lVert\kappa\rVert_{L^{\infty}(\Gamma_{R})}
\lVert\gamma y\rVert_{L^{2}(\Gamma_{R})}
\lVert\gamma v\rVert_{L^{2}(\Gamma_{R})}
\\
&\le
C_{\mathrm{tr}}^{2}
\lVert\kappa\rVert_{L^{\infty}(\Gamma_{R})}
\lVert y\rVert_{H^{1}(\Omega)}
\lVert v\rVert_{H^{1}(\Omega)}.
\end{aligned}
```

So the Robin boundary contribution is a bounded bilinear form on $H^{1}(\Omega)$.

> **Natural does not mean automatically well posed.** If the entire boundary carries a pure Neumann condition, constants have zero gradient, so the Poisson bilinear form is not coercive on all of $H^{1}(\Omega)$. The weak equation tested with $v=1$ also forces the compatibility condition
> $$
> 0
> =
> \int_{\Omega}f\mathrm{d}x
> +
> \int_{\partial\Omega}g\mathrm{d}s.
> $$
> When that condition holds, the solution is determined only up to an additive constant unless one imposes an additional normalization. This is a concrete example of why the function space and the well-posedness assumptions must be checked together rather than inferred from the differential equation alone.

## 11. Lax–Milgram: when the variational problem is well posed

Deriving a weak formulation proves that the expressions make sense under the chosen assumptions. It does **not** by itself prove that there is exactly one solution. In this chapter, **well posed** means that a solution exists, is unique, and depends continuously on the data in the chosen norms. For coercive linear problems, the standard result establishing these properties is the **Lax–Milgram theorem**.

Let $V$ be a real Hilbert space. Assume:

1. $a:V\times V\to\mathbb R$ is bilinear and bounded, so there is an $M>0$ with
   $$
   |a(y,v)|
   \le
   M\lVert y\rVert_{V}\lVert v\rVert_{V};
   $$
2. $a$ is coercive, so there is an $\alpha>0$ with
   $$
   a(v,v)
   \ge
   \alpha\lVert v\rVert_{V}^{2};
   $$
3. $F\in V^{\ast}$ is a bounded linear functional.

Under **all three assumptions**, there exists a unique $y\in V$ such that

$$
a(y,v)
=
\langle F,v\rangle_{V^{\ast},V}
\qquad
\text{for every }v\in V.
$$

Moreover,

$$
\lVert y\rVert_{V}
\le
\frac{1}{\alpha}
\lVert F\rVert_{V^{\ast}}.
$$

The theorem therefore gives existence, uniqueness, and stability in one statement. The proof of **existence** is the substantial functional-analytic part and is outside this prerequisite chapter. Two other parts of the logic are cheap enough to see directly.

### Why coercivity gives uniqueness

Suppose $y_{1}$ and $y_{2}$ both solve the same problem. Their difference $w:=y_{1}-y_{2}$ satisfies

$$
a(w,v)=0
\qquad
\text{for every }v\in V.
$$

Choose $v=w$. Then

$$
0
=
a(w,w)
\ge
\alpha\lVert w\rVert_{V}^{2}.
$$

Hence $\lVert w\rVert_{V}=0$, so $w=0$ and $y_{1}=y_{2}$.

### Where the stability estimate comes from

Choose $v=y$ in the variational equation. Coercivity and the definition of the dual norm give

```math
\begin{aligned}
\alpha\lVert y\rVert_{V}^{2}
&\le
a(y,y)
\\
&=
\langle F,y\rangle_{V^{\ast},V}
\\
&\le
\lVert F\rVert_{V^{\ast}}
\lVert y\rVert_{V}.
\end{aligned}
```

If $y\ne0$, divide by $\lVert y\rVert_{V}$ to obtain

$$
\lVert y\rVert_{V}
\le
\frac{1}{\alpha}
\lVert F\rVert_{V^{\ast}}.
$$

If $y=0$, the same estimate is trivially true. This calculation shows exactly where the coercivity constant enters the stability bound.

Boundedness of $a$ and $F$ plays a different role: it says that the two sides of the variational equation depend continuously on the arguments in the topology of $V$. Coercivity supplies a uniform lower bound that prevents nonzero directions from becoming invisible to the bilinear form. Lax–Milgram combines these properties to obtain existence as well.

### Checking the homogeneous-Dirichlet Poisson problem

Take

$$
V=H_{0}^{1}(\Omega)
$$

with

$$
\lVert v\rVert_{V}
:=
\lVert\nabla v\rVert_{L^{2}(\Omega)^{d}}.
$$

Poincaré tells us that this really is a norm equivalent to the usual $H^{1}$ norm. For

$$
a(y,v)
=
\int_{\Omega}\nabla y\cdot\nabla v\mathrm{d}x,
$$

Cauchy–Schwarz gives

$$
|a(y,v)|
\le
\lVert y\rVert_{V}\lVert v\rVert_{V},
$$

so $a$ is bounded with $M=1$. Also,

$$
a(v,v)
=
\int_{\Omega}|\nabla v|^{2}\mathrm{d}x
=
\lVert v\rVert_{V}^{2},
$$

so $a$ is coercive with $\alpha=1$.

For $f\in L^{2}(\Omega)$, Section 7 showed that

$$
F(v)=\int_{\Omega}fv\mathrm{d}x
$$

is a bounded linear functional on $V$. Every Lax–Milgram hypothesis is therefore satisfied. We conclude that there is one and only one weak solution

$$
y\in H_{0}^{1}(\Omega).
$$

For this concrete forcing, we can obtain a slightly more explicit estimate. Testing with $v=y$ gives

$$
\lVert\nabla y\rVert_{L^{2}(\Omega)^{d}}^{2}
=
\int_{\Omega}fy\mathrm{d}x.
$$

Then

```math
\begin{aligned}
\lVert\nabla y\rVert_{L^{2}(\Omega)^{d}}^{2}
&\le
\lVert f\rVert_{L^{2}(\Omega)}
\lVert y\rVert_{L^{2}(\Omega)}
\\
&\le
C_{P}
\lVert f\rVert_{L^{2}(\Omega)}
\lVert\nabla y\rVert_{L^{2}(\Omega)^{d}}.
\end{aligned}
```

Hence

$$
\lVert\nabla y\rVert_{L^{2}(\Omega)^{d}}
\le
C_{P}\lVert f\rVert_{L^{2}(\Omega)}.
$$

This is a stability statement for the state equation: the solution depends continuously on the forcing in these norms.

## 12. Weak solutions, strong solutions, and regularity

A weak formulation asks less differentiability of the solution than the strong PDE. It is therefore useful to separate two directions of implication.

If a classical solution exists, it is also a weak solution: multiply the strong equation by a test function and integrate by parts, exactly as in Section 9.

Conversely, suppose a weak solution is smooth enough that $\Delta y$ is an ordinary integrable function. For every $v\in C_{c}^{\infty}(\Omega)$, the weak equation and integration by parts give

```math
\begin{aligned}
0
&=
\int_{\Omega}\nabla y\cdot\nabla v\mathrm{d}x
-
\int_{\Omega}fv\mathrm{d}x
\\
&=
\int_{\Omega}(-\Delta y-f)v\mathrm{d}x.
\end{aligned}
```

A standard fundamental result for test functions says that if an integrable function $q$ satisfies

$$
\int_{\Omega}qv\mathrm{d}x=0
\qquad
\text{for every }v\in C_{c}^{\infty}(\Omega),
$$

then $q=0$ almost everywhere. Applying it to $q=-\Delta y-f$ recovers

$$
-\Delta y=f
\qquad
\text{a.e. in }\Omega.
$$

So weak and strong formulations agree whenever the weak solution has enough regularity to justify the strong derivatives and boundary interpretation.

The important point is that this extra regularity is **not automatic from the definition of weak solution**. Initially, Lax–Milgram gives only

$$
y\in H_{0}^{1}(\Omega).
$$

Whether $y$ actually belongs to $H^{2}(\Omega)$, is continuous, or has still higher derivatives depends on the coefficients, the data, the geometry of the domain, and the boundary conditions. Establishing such improvements is the subject of **elliptic regularity theory**.

For some smooth elliptic problems on sufficiently regular or convex domains, one can prove estimates of the schematic form

$$
\lVert y\rVert_{H^{2}(\Omega)}
\le
C\lVert f\rVert_{L^{2}(\Omega)}.
$$

On domains with re-entrant corners, with rough coefficients, or with rougher data, such an estimate may fail or need modification. The weak formulation remains meaningful even when this stronger regularity does not hold.

This distinction matters numerically. The weak formulation tells us what problem is well defined at the energy-space level. Additional regularity may justify pointwise observations, stronger norms, or sharper approximation estimates, but it is a separate analytical fact and should not be silently assumed.

A useful mental hierarchy is

```text
strong PDE statement
        ↓ derive by testing and integration by parts
weak / variational problem
        ↓ prove well-posedness in a function space
weak solution
        ↓ add separate regularity assumptions and theorems
possibly smoother / strong solution
```

The arrows indicate logical steps with hypotheses, not unconditional equivalences.

## 13. Scope frontier

The main path has now introduced the functional-analysis machinery needed to understand the first weak PDE and finite-element examples. Several nearby subjects are important but not required before continuing.

We have not constructed Lebesgue measure or the $L^{p}$ spaces from first principles. We have not proved the trace theorem, Poincaré inequality, Sobolev embeddings, or Lax–Milgram. We have not developed compactness and weak convergence, general distribution theory, detailed elliptic regularity, or the inf-sup theory required by mixed problems.

Those deferrals are deliberate, but they need not be black boxes. [01a · Further notes on function spaces and weak PDEs](01a-further-function-space-notes.md) give a compact first look at:

- the trace space $H^{1/2}(\partial\Omega)$ and liftings;
- Sobolev embeddings and when point values become meaningful;
- weak convergence and compactness;
- distributions and singular data;
- elliptic regularity;
- what replaces coercivity for mixed and saddle-point problems.

That companion file is optional: it is there to deepen the surrounding picture without interrupting the prerequisite chain of this chapter.

We have also only touched the dual space. The fact that residuals and derivatives live in dual spaces becomes central once we discuss coordinates, gradients, and adjoints. That is the subject of the next background chapter.

## Where to go next

- **Default continuation:** [02 · Duality, derivatives, and adjoints](02-duality-derivatives-and-adjoints.md) develops the dual-space, derivative, Riesz-map, and adjoint language needed by optimization.
- **If your next goal is discretization:** [03 · Finite elements](03-finite-elements.md) turns the weak Poisson problem into a finite-dimensional Galerkin system and introduces basis functions, degrees of freedom, assembly, and quadrature.
- **For deeper function-space theory:** [01a · Further notes on function spaces and weak PDEs](01a-further-function-space-notes.md) adds higher Sobolev spaces, embeddings, compactness, distributions, elliptic regularity, and the first inf-sup ideas without blocking the main route.
- **Into the `nmopt` manual:** these are the closest project-specific continuations:
  - [manual 01 · Anatomy of a discrete PDE-constrained problem](../concepts/01-discrete-problem-anatomy.md) uses weak residuals and boundary conditions as the starting point for a discrete OCP;
  - [manual 02 · Spaces, coordinates, and duality](../concepts/02-spaces-coordinates-and-duality.md) continues from function spaces to discrete coordinates, primal values, and covectors.

## References and further reading

Publication details, stable links, and access notes are collected in the [Background reference catalogue](references.md).

- **Stefano Mancini, Andrea Benvenuti, and Luca Heltai, *Numerical Methods for Partial Differential Equations*.**
  - §1.1.1 for Lebesgue spaces and §1.1.2 for weak derivatives.
  - §§1.1.3–1.1.4 for Sobolev spaces and $H_{0}^{1}$.
  - §1.1.5 for dual spaces and §§1.1.7–1.1.8 for bilinear forms and traces.
  - §§1.2.1–1.2.3 for the strong Poisson problem, weak formulation, and Lax–Milgram.
  - §§1.2.4–1.2.5 are a useful bridge through minimization and Céa's lemma toward [03 · Finite elements](03-finite-elements.md).
- **Andrea Manzoni, Alfio Quarteroni, and Sandro Salsa, *Optimal Control of Partial Differential Equations: Analysis, Approximation, and Applications*.**
  - Appendix A.1 for Banach and Hilbert spaces.
  - Appendix A.2.1–A.2.2 for bounded operators, functionals, dual spaces, and Riesz representation.
  - Appendix A.4.1 for Lax–Milgram.
  - Appendix A.5.1 and A.5.3–A.5.7 for Lebesgue spaces, weak derivatives, $H^{1}$, $H_{0}^{1}$, and traces.
  - §5.1.1 for the weak state equation in an elliptic control problem and §5.4 for Dirichlet, Neumann, Robin, and mixed boundary formulations.
