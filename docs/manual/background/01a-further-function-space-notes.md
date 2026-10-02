# 01a · Further notes on function spaces and weak PDEs

**Background navigation:** [Index](README.md) \
Main chapter: [01 · Function spaces and weak PDEs](01-function-spaces-and-weak-pdes.md) \
Previous: [01 · Function spaces and weak PDEs](01-function-spaces-and-weak-pdes.md) \
Next: [02 · Duality, derivatives, and adjoints](02-duality-derivatives-and-adjoints.md)

## Purpose

[01 · Function spaces and weak PDEs](01-function-spaces-and-weak-pdes.md) keeps a deliberately narrow main path: enough $L^{2}$, Sobolev, trace, variational, and well-posedness language to derive and understand the weak Poisson problem. Several nearby ideas become important in more advanced PDE-constrained optimization problems, but developing them in the middle of that derivation would obscure the prerequisite chain.

This companion note gives those ideas a first home. It is not a second functional-analysis course. Its purpose is to make later notation recognizable, explain what problem each concept solves, and indicate where a full treatment begins.

## Before you start

Read [01 · Function spaces and weak PDEs](01-function-spaces-and-weak-pdes.md), or be comfortable with $L^{2}(\Omega)$, $H^{1}(\Omega)$, $H_{0}^{1}(\Omega)$, weak derivatives, traces at the level of the trace theorem, and coercive variational problems.

## What you will be able to do

After this note, you should be able to:

- read the notation $W^{k,p}(\Omega)$ and $H^{k}(\Omega)$ in several spatial dimensions;
- distinguish a Sobolev norm from its highest-order seminorm and explain why the distinction matters;
- give a concrete interpretation of the trace space $H^{1/2}(\partial\Omega)$;
- explain why an $H^{1}$ function need not have meaningful point values in several dimensions;
- distinguish strong convergence from weak convergence and explain why compactness matters;
- recognize distributions and singular data such as Dirac point sources as a genuine extension of weak-derivative ideas;
- explain what elliptic regularity adds beyond Lax–Milgram;
- recognize when a variational problem falls outside the coercive Lax–Milgram setting and requires an inf-sup condition instead.

## Roadmap

We first make trace spaces and liftings more concrete. We then fill in two pieces of Sobolev notation that the main chapter only mentioned briefly: higher-order spaces in several dimensions and Sobolev seminorms. Sobolev embeddings then answer when integral regularity implies continuity. Weak convergence and compactness explain how infinite-dimensional existence arguments recover convergent subsequences. A short section on distributions shows how weak differentiation extends to singular data. Finally, elliptic regularity and mixed variational problems indicate two important directions beyond the coercive Poisson model.

## 1. Trace spaces and the meaning of $H^{1/2}(\partial\Omega)$

For a bounded Lipschitz domain, the trace theorem gives a bounded surjective map

$$
\gamma\colon H^{1}(\Omega)\to H^{1/2}(\partial\Omega).
$$

A useful way to read this statement is to define the boundary space by the traces that can actually occur:

$$
H^{1/2}(\partial\Omega)\coloneqq\gamma\bigl(H^{1}(\Omega)\bigr).
$$

This tells us what its elements are, but we also need a norm. A boundary value $g$ may have many different interior extensions $v\in H^{1}(\Omega)$ satisfying $\gamma v=g$. The natural quotient norm is

$$
\lVert g\rVert_{H^{1/2}(\partial\Omega)}\coloneqq
\inf_{\substack{v\in H^{1}(\Omega)\\\gamma v=g}}
\lVert v\rVert_{H^{1}(\Omega)}.
$$

Thus the $H^{1/2}$ size of a boundary field measures how cheaply it can be realized as the trace of an $H^{1}$ interior field. With this norm,

$$
\lVert\gamma v\rVert_{H^{1/2}(\partial\Omega)}
\leq
\lVert v\rVert_{H^{1}(\Omega)}.
$$

The exponent $1/2$ should not be read as a pointwise “half derivative.” Fractional Sobolev spaces can be characterized intrinsically by nonlocal difference quotients, but the trace characterization above is often the most useful one when first meeting boundary control and Dirichlet data.

The kernel of the trace is

$$
\ker\gamma=H_{0}^{1}(\Omega).
$$

This gives a compact structural picture:

```text
H¹(Ω)
  │ trace γ
  ▼
H¹ᐟ²(∂Ω)

kernel of γ = H₀¹(Ω)
```

### Liftings

If $g\in H^{1/2}(\partial\Omega)$, surjectivity of the trace means that there is at least one $\widetilde g\in H^{1}(\Omega)$ with

$$
\gamma\widetilde g=g.
$$

Such a function is a **lifting** or extension of the boundary datum. A nonhomogeneous Dirichlet problem can then be rewritten as

$$
y=z+\widetilde g,
\qquad
z\in H_{0}^{1}(\Omega).
$$

This is more than a proof trick. Numerically, the same decomposition separates prescribed physical values from independent unknown coordinates.

There is also a dual boundary space,

$$
H^{-1/2}(\partial\Omega)\coloneqq\left(H^{1/2}(\partial\Omega)\right)^{\ast},
$$

which is the natural home for weak normal fluxes in many elliptic problems. A full treatment requires more duality than the first background chapter develops, so this notation is best revisited after [02 · Duality, derivatives, and adjoints](02-duality-derivatives-and-adjoints.md).

## 2. Higher-order Sobolev spaces and seminorms

The main chapter defined $H^{1}(\Omega)$ directly because first weak derivatives are all the Poisson weak form needs. Finite-element error estimates and more regular PDEs often use higher orders, so it is useful to make the general notation explicit.

Let

$$
\Omega\subset\mathbb{R}^{d}
$$

and let

$$
\alpha=(\alpha_{1},\ldots,\alpha_{d})
$$

be a **multi-index**, meaning that each $\alpha_{i}$ is a nonnegative integer. Its order is

$$
\lvert\alpha\rvert\coloneqq\alpha_{1}+\cdots+\alpha_{d}.
$$

The notation

$$
D^{\alpha}u
$$

means that we differentiate $u$ $\alpha_{1}$ times with respect to $x_{1}$, $\alpha_{2}$ times with respect to $x_{2}$, and so on, with derivatives understood weakly whenever classical derivatives are unavailable. For example, in two dimensions,

$$
D^{(1,0)}u=\frac{\partial u}{\partial x_{1}},
\qquad
D^{(0,1)}u=\frac{\partial u}{\partial x_{2}},
$$

while the derivatives of total order two are

$$
D^{(2,0)}u,
\qquad
D^{(1,1)}u,
\qquad
D^{(0,2)}u.
$$

For an integer $k\geq 0$ and $1\leq p\lt\infty$, the Sobolev space $W^{k,p}(\Omega)$ is

$$
W^{k,p}(\Omega)\coloneqq\left\lbrace
 u\in L^{p}(\Omega):
 D^{\alpha}u\in L^{p}(\Omega)
 \text{ for every }\lvert\alpha\rvert\leq k
\right\rbrace.
$$

A standard norm is

$$
\lVert u\rVert_{W^{k,p}(\Omega)}\coloneqq
\left(
\sum_{\lvert\alpha\rvert\leq k}
\lVert D^{\alpha}u\rVert_{L^{p}(\Omega)}^{p}
\right)^{1/p}.
$$

When $p=2$, we write

$$
H^{k}(\Omega)\coloneqq W^{k,2}(\Omega).
$$

Thus $H^{2}(\Omega)$ in two dimensions asks for $u$, both first weak derivatives, and all second weak derivatives to lie in $L^{2}(\Omega)$. Schematically,

```text
H² in two dimensions
    u
    first order:   ∂₁u, ∂₂u
    second order:  ∂₁₁u, ∂₁₂u, ∂₂₂u
    all square integrable
```

### Seminorms

Often we want to measure only derivatives of one particular order. The **Sobolev seminorm** of order $k$ is

$$
\lvert u\rvert_{W^{k,p}(\Omega)}\coloneqq
\left(
\sum_{\lvert\alpha\rvert=k}
\lVert D^{\alpha}u\rVert_{L^{p}(\Omega)}^{p}
\right)^{1/p}.
$$

For $p=2$ this becomes $\lvert u\rvert_{H^{k}(\Omega)}$. In particular,

$$
\lvert u\rvert_{H^{1}(\Omega)}
=\lVert\nabla u\rVert_{L^{2}(\Omega)^{d}}.
$$

The word *seminorm* matters. A norm can vanish only on the zero element, whereas a seminorm may vanish on nonzero functions. For example, every constant $c$ satisfies

$$
\lvert c\rvert_{H^{1}(\Omega)}=0,
$$

because $\nabla c=0$, even though $c$ need not be the zero function. More generally, the order-$k$ seminorm cannot detect polynomials whose derivatives of order $k$ vanish.

The full Sobolev norm restores the lower-order information. For $H^{1}$,

$$
\lVert u\rVert_{H^{1}(\Omega)}^{2}
=\lVert u\rVert_{L^{2}(\Omega)}^{2}+\lvert u\rvert_{H^{1}(\Omega)}^{2}.
$$

There are important subspaces on which a seminorm becomes a genuine norm. [01 · Function spaces and weak PDEs](01-function-spaces-and-weak-pdes.md) already encountered the basic example: Poincaré's inequality gives

$$
\lVert u\rVert_{L^{2}(\Omega)}\leq C_{P}\lvert u\rvert_{H^{1}(\Omega)}
\qquad
\text{for }u\in H_{0}^{1}(\Omega).
$$

Hence $\lvert u\rvert_{H^{1}(\Omega)}$ is a norm on $H_{0}^{1}(\Omega)$ and is equivalent there to the full $H^{1}$ norm.

This notation becomes especially common in finite-element analysis. An estimate such as

$$
\lVert u-I_{h}u\rVert_{H^{1}(\Omega)}\leq C h\lvert u\rvert_{H^{2}(\Omega)}
$$

separates two different roles: the left-hand side measures the approximation error in the full $H^{1}$ norm, while the right-hand side measures the second-order variation of the exact function. The derivation of such estimates belongs to finite-element approximation theory; the point here is simply to make the notation and its meaning explicit before it is used.

## 3. Sobolev embeddings: when do point values make sense?

The $H^{1}$ norm controls a function and its first weak derivatives in an integral sense. It does not automatically make the function continuous.

This matters because some operations are pointwise. An observation such as

$$
y(x_{0})
$$

only makes sense as an intrinsic quantity if the function-space element has a well-defined continuous representative, or if some other theory gives the point evaluation a meaning.

Sobolev embedding theorems relate the number of weak derivatives, their integrability, and the spatial dimension. A useful schematic rule is that sufficiently many derivatives relative to the dimension force continuity. For instance, on a regular bounded domain, if

$$
kp\gt d,
$$

then $W^{k,p}(\Omega)$ embeds into a space of continuous functions under the standard hypotheses.

This immediately explains a dimension-dependent fact that often surprises new readers. In two dimensions,

$$
H^{1}(\Omega)=W^{1,2}(\Omega)
$$

does **not** satisfy $1\cdot2\gt 2$, so $H^{1}$ regularity alone does not guarantee continuity. On the other hand,

$$
H^{2}(\Omega)=W^{2,2}(\Omega)
$$

has $2\cdot2\gt 2$ and, under the usual domain assumptions, embeds continuously into $C^{0}(\overline\Omega)$.

This is why a weak state in $H_{0}^{1}(\Omega)$ may be perfectly adequate for integral observations such as

$$
\int_{\Omega}qy\thinspace\mathrm{d}x
$$

while a point observation $y(x_{0})$ can require additional regularity.

The exact embedding theorems have borderline cases and domain assumptions that matter. The purpose here is not to memorize the full catalogue; it is to understand why “the state exists in $H^{1}$” and “the state can be evaluated pointwise” are different analytical statements.

## 4. Weak convergence and compactness

In a finite-dimensional space, every bounded sequence has a convergent subsequence. That fact is so familiar that it is easy to use it unconsciously. It fails in infinite-dimensional normed spaces.

For example, let $\lbrace e_{k}\rbrace$ be an orthonormal sequence in an infinite-dimensional Hilbert space. Then

$$
\lVert e_{j}-e_{k}\rVert^{2}=2
\qquad
\text{for }j\ne k,
$$

so no subsequence can converge strongly in norm.

A weaker notion of convergence restores some compactness. In a Hilbert space $H$, we say

$$
u_{k}\rightharpoonup u
\qquad
\text{weakly in }H
$$

if

$$
(u_{k},v)_{H}\to(u,v)_{H}
\qquad
\text{for every }v\in H.
$$

Strong convergence implies weak convergence, but not conversely. An orthonormal sequence satisfies

$$
e_{k}\rightharpoonup0
$$

although every $e_{k}$ has norm one.

Why is this useful? Hilbert spaces are reflexive, and bounded sequences in reflexive spaces admit weakly convergent subsequences. This is one of the standard tools behind existence proofs in optimization: first obtain a bounded minimizing sequence, then extract a weakly convergent subsequence, and finally show that the objective and constraints behave well enough under that convergence.

Sometimes weak convergence is not enough because a nonlinear term needs strong convergence. **Compact embeddings** bridge the gap. For a bounded Lipschitz domain, the embedding

$$
H^{1}(\Omega)\hookrightarrow L^{2}(\Omega)
$$

is compact: a sequence bounded in $H^{1}$ has a subsequence converging strongly in $L^{2}$. This compactness result is much stronger than the elementary continuous estimate

$$
\lVert v\rVert_{L^{2}(\Omega)}\leq\lVert v\rVert_{H^{1}(\Omega)}.
$$

The first background chapter does not need these ideas to establish Poisson well-posedness, but they become central when proving existence of optimal controls or passing to limits in nonlinear PDEs.

## 5. Distributions and singular data

The weak derivative definition already hints at a broader philosophy: instead of asking an object to possess pointwise derivatives, define differentiation through how it acts on test functions.

A **distribution** is, roughly, a continuous linear functional on the test-function space $C_{c}^{\infty}(\Omega)$. Every locally integrable function $q$ defines a distribution through

$$
\varphi
\mapsto
\int_{\Omega}q\varphi\thinspace\mathrm{d}x.
$$

But distributions also include objects that are not ordinary functions. The Dirac delta at a point $x_{0}$ is defined by

$$
\langle\delta_{x_{0}},\varphi\rangle\coloneqq\varphi(x_{0}).
$$

Its derivative is then defined by moving the derivative onto the test function:

$$
\left\langle
\frac{\partial\delta_{x_{0}}}{\partial x_{i}},
\varphi
\right\rangle\coloneqq -\left\langle
\delta_{x_{0}},
\frac{\partial\varphi}{\partial x_{i}}
\right\rangle.
$$

This is the same integration-by-parts principle used for weak derivatives, now taken as the fundamental definition.

The extension is useful because PDEs can have singular forcing or observations. A point source formally written as

$$
-\Delta y=\delta_{x_{0}}
$$

cannot be treated as an $L^{2}$ forcing problem. Likewise, an objective depending on $y(x_{0})$ introduces a point-evaluation functional whose adjoint forcing may involve a Dirac measure. Such problems are mathematically legitimate, but the correct state and residual spaces can differ from the simple $H_{0}^{1}$/$H^{-1}$ setting of the main chapter.

A full distribution theory is outside the prerequisite scope. The useful lesson is that weak formulations are not merely a trick for slightly nonsmooth functions; they are part of a framework that can accommodate genuinely singular data when the spaces are chosen carefully.

## 6. Elliptic regularity: what happens after existence?

Lax–Milgram answers an energy-space question. For homogeneous Dirichlet Poisson with suitable forcing, it gives

$$
y\in H_{0}^{1}(\Omega)
$$

and a stability estimate in the $H^{1}$-type norm.

It does not say that $y$ has second derivatives in $L^{2}$, is continuous, or possesses a classical normal derivative. Those are **regularity** questions.

Under stronger assumptions on the domain and data, elliptic regularity can upgrade the solution. A common model result is

$$
f\in L^{2}(\Omega)
\quad\Longrightarrow\quad
y\in H^{2}(\Omega)
$$

with an estimate

$$
\lVert y\rVert_{H^{2}(\Omega)}\leq C\lVert f\rVert_{L^{2}(\Omega)},
$$

provided the domain and boundary conditions satisfy the required regularity hypotheses. Convex polygonal domains often enjoy useful $H^{2}$ regularity for Poisson; re-entrant corners can create singularities that reduce it.

Regularity matters numerically for at least three reasons. It controls the approximation rates available to finite elements, it determines whether pointwise or boundary quantities are meaningful, and it affects how smooth an adjoint or sensitivity equation can be expected to be.

The habit to keep is simple: **well posed in $H^{1}$** and **smooth enough for the next operation** are separate claims. Each needs its own assumptions.

## 7. When coercivity is not the right well-posedness condition

The Poisson problem used in the main chapter has one trial space, one test space, and a coercive bilinear form. Many important PDEs do not fit that pattern.

The stationary Stokes equations, for example, couple a velocity $v$ with a pressure $\pi$. A typical weak formulation has the schematic structure

$$
\begin{aligned}
a(v,w)+b(w,\pi)&=F(w),\\
b(v,q)&=0,
\end{aligned}
$$

for every velocity test function $w$ and pressure test function $q$.

The associated block problem is a **mixed** or **saddle-point** variational problem. The full bilinear form is not coercive on the whole product space, so the ordinary Lax–Milgram theorem is not the right tool. Well posedness is instead controlled by coercivity on a suitable kernel together with an **inf-sup**, or Ladyzhenskaya–Babuška–Brezzi, condition of the form

$$
\inf_{q\ne 0}
\sup_{v\ne 0}
\frac{b(v,q)}
{\lVert v\rVert_{V}\lVert q\rVert_{Q}}
\geq
\beta\gt 0.
$$

The inequality says, informally, that no nonzero multiplier direction $q$ is invisible to all admissible primal directions $v$. It plays a role analogous to a stability condition for the coupling.

We do not need to develop this theory before reading the current core `nmopt` material. It becomes important for mixed PDEs such as Stokes, for some Petrov–Galerkin discretizations, and for broader saddle-point formulations. The numerical linear-algebra background chapter will revisit the corresponding block systems from a computational point of view.

## Where to go next

- **Return to the main route:** continue with [02 · Duality, derivatives, and adjoints](02-duality-derivatives-and-adjoints.md) for duality and adjoints, or with [03 · Finite elements](03-finite-elements.md) if your immediate goal is discretization.
- **Where these deferred ideas reappear:** traces and $H^{1/2}$ return in boundary control, compactness in existence arguments, distributions in point observations, regularity in stronger outputs, and inf-sup conditions in mixed systems. [07a · Further notes on PDE-constrained optimization](07a-further-pde-constrained-optimization-notes.md) and [07b · Extensions of PDE-constrained optimization](07b-extensions-of-pde-constrained-optimization.md) show several of those uses in PDE optimal control.

## References and further reading

Publication details, stable links, and access notes are collected in the [Background reference catalogue](references.md).

- **Mancini, Benvenuti, and Heltai, *Numerical Methods for Partial Differential Equations*.**
  - §1.1.8 for trace spaces and §1.1.9 for Sobolev embeddings.
  - The later mixed/Petrov–Galerkin material for inf-sup theory.
  - The approximation chapters for the effect of regularity on finite-element errors.
- **Manzoni, Quarteroni, and Salsa, *Optimal Control of Partial Differential Equations*.**
  - Appendix A.3 for weak convergence and compactness.
  - Appendix A.4.2 for saddle-point problems and inf-sup conditions.
  - Appendix A.5.7 for traces, Appendix A.5.8 for compactness, and Appendix A.5.11 for Sobolev embeddings.
  - §5.11 for trace spaces in Dirichlet boundary control.
  - §§5.9–5.10 for singular/pointwise observations and lower-regularity formulations.

Use the cited source sections when precise hypotheses or full proofs become important; this companion is intended as orientation and working background rather than a replacement for those treatments.
