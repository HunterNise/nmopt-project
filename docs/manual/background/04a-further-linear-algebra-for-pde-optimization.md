# 04a · Further numerical linear algebra for PDE-constrained optimization

**Background navigation:** [Index](README.md) \
Main chapter: [04 · Numerical linear algebra for PDE and optimization systems](04-numerical-linear-algebra.md) \
Previous: [04 · Numerical linear algebra for PDE and optimization systems](04-numerical-linear-algebra.md) \
Next: [05 · Unconstrained numerical optimization](05-unconstrained-optimization.md)

## Purpose

[04 · Numerical linear algebra for PDE and optimization systems](04-numerical-linear-algebra.md) develops the solver vocabulary needed by the main background track: matrix structure, conditioning, Krylov methods, preconditioning, matrix-free actions, block systems, Schur complements, and KKT saddle-point systems.

This companion note goes one step further in the direction most relevant to PDE-constrained optimization. The goal is not to provide another general numerical-linear-algebra course. Instead, it explains a small set of ideas that repeatedly appear when PDE discretizations and optimization are combined:

- the linear system is usually one member of a family indexed by mesh size or another discretization parameter;
- a useful preconditioner should remain effective as that family changes;
- state and adjoint solvers often become building blocks inside larger KKT preconditioners;
- block preconditioners are guided by exact factorizations but use cheap approximate solves;
- robustness may be required with respect to both mesh refinement and optimization parameters;
- inner solves used inside a preconditioner need not always be exact.

The main [04 · Numerical linear algebra for PDE and optimization systems](04-numerical-linear-algebra.md) chapter remains sufficient for the core prerequisite sequence. This file is for readers who want to understand why PDE-optimization papers discuss **mesh-independent**, **spectrally equivalent**, **block**, or **parameter-robust** preconditioners.

## Before you start

You should already be comfortable with:

- the main [04 · Numerical linear algebra for PDE and optimization systems](04-numerical-linear-algebra.md) chapter, especially conditioning, Krylov methods, preconditioning, Schur complements, and saddle-point systems;
- [03 · Finite elements](03-finite-elements.md) for the fact that mesh refinement produces a sequence of discrete operators;
- the distinction between primal and dual actions from [02 · Duality, derivatives, and adjoints](02-duality-derivatives-and-adjoints.md).

We use symmetric positive-definite blocks when that makes an idea transparent. General KKT systems can require weaker assumptions, nonsymmetric solvers, nullspace treatment, or problem-specific analysis.

## What you will be able to do

After this note, you should be able to:

- explain why one should study a family $A_{h}$ rather than one fixed matrix when a PDE mesh is refined;
- define spectral equivalence for SPD operators and relate it to mesh-independent conditioning;
- interpret an approximate PDE solve as a preconditioner action rather than as an explicit inverse matrix;
- derive the ideal block factorization of a simple equality-constrained KKT system;
- explain why practical block preconditioners approximate both the primal block and its Schur complement;
- distinguish mesh robustness from robustness with respect to an optimization parameter such as a regularization weight;
- explain why nested iterative solves can make the effective preconditioner vary from one outer iteration to the next;
- recognize multigrid and Schwarz/domain-decomposition methods as scalable approximate-inverse strategies, without needing their full construction.

## Roadmap

We first replace the idea of “the matrix” by a family of discretized systems $A_{h}$. Spectral equivalence then gives a compact way to state that a preconditioner controls this family uniformly. We next return to equality-constrained KKT systems and use an exact block factorization to identify the algebraic ingredients that an effective preconditioner should imitate. The final sections discuss robustness with respect to mesh and optimization parameters, followed by the practical issue of inexact inner solves. We then give a short orientation to multigrid and Schwarz/domain-decomposition ideas, because these names occur constantly in PDE solver discussions even when their detailed analysis is not part of the prerequisite track.

## 1. PDE discretization produces families of linear systems

A finite-element matrix depends on the discrete space. If the mesh changes, so does the matrix:

$$
A_{h}x_{h}=b_{h}.
$$

The subscript $h$ should be read as a reminder that this is one member of a family. As $h$ decreases,

- the dimension usually grows;
- the sparsity pattern changes;
- the spectrum can spread;
- the condition number can deteriorate;
- the cost of one matrix-vector product changes.

A solver that performs well for one coarse matrix is therefore not automatically a good PDE solver.

For standard elliptic discretizations, the unpreconditioned condition number often grows as the mesh is refined. The exact power of $h$ depends on the operator, norm, basis, and discretization, but the qualitative point is enough here:

> Refining the approximation space can make the algebraic problem harder even while it makes the discretization more accurate.

This is why iteration counts are often reported against mesh size. The important question is not merely whether a Krylov method converges for each fixed $h$, but whether the number of iterations remains controlled as the discrete problem becomes larger.

## 2. Spectral equivalence expresses uniform quality for SPD problems

Suppose $A_{h}$ and $P_{h}$ are SPD matrices. We say they are **spectrally equivalent**, uniformly in $h$, if there exist constants

$$
0<c\le C<\infty
$$

independent of $h$ such that

$$
c v^{\mathsf T}P_{h}v
\le
v^{\mathsf T}A_{h}v
\le
C v^{\mathsf T}P_{h}v
$$

for every vector $v$.

Consider the generalized eigenvalue problem

$$
A_{h}v
=
\lambda P_{h}v.
$$

Pairing with $v$ gives

$$
\lambda
=
\frac{v^{\mathsf T}A_{h}v}{v^{\mathsf T}P_{h}v}.
$$

The spectral-equivalence inequalities therefore imply

$$
c\le\lambda\le C.
$$

Equivalently, the symmetrically preconditioned operator

$$
P_{h}^{-1/2}A_{h}P_{h}^{-1/2}
$$

has all eigenvalues in $[c,C]$, so

$$
\kappa_{2}
\left(
P_{h}^{-1/2}A_{h}P_{h}^{-1/2}
\right)
\le
\frac{C}{c}.
$$

If $c$ and $C$ do not deteriorate as $h\to0$, neither does this bound.

This is the mathematical content behind phrases such as **mesh-independent preconditioning** or **optimal preconditioning** in the SPD setting. It does not mean that the total cost of a solve is independent of mesh size: applying the operator and preconditioner still becomes more expensive as the problem grows. It means that the outer iteration is not becoming harder merely because the discretization was refined.

### 2.1 The preconditioner is judged as an action

Writing $P_{h}^{-1}$ does not require forming an inverse matrix. In practice,

$$
z=P_{h}^{-1}r
$$

means “apply an algorithm that approximately solves a problem represented by $P_{h}$.”

That action could be one or a few multigrid cycles, an incomplete factorization, a sparse direct solve on a subproblem, a few stationary iterations, or another problem-specific approximate inverse.

The main chapter deliberately leaves the construction of these methods out of scope. For PDE-constrained optimization, the important point is that they can be treated as reusable solve actions inside larger block algorithms.

## 3. State and adjoint solves are natural preconditioning ingredients

A reduced PDE-constrained optimization method repeatedly solves state and adjoint equations. An all-at-once method exposes those same PDE operators as blocks of a larger coupled system.

Suppose, schematically, that the state Jacobian is represented by

$$
A_{h}.
$$

The state equation needs an action approximating

$$
A_{h}^{-1},
$$

while the adjoint equation needs an action approximating

$$
A_{h}^{-\mathsf T}.
$$

If efficient PDE solvers for these actions already exist, they are natural building blocks for a KKT preconditioner.

This is an important change of viewpoint:

```text
standalone PDE solver
        ↓
approximate inverse action
        ↓
building block inside a larger optimization solve
```

The optimization algorithm does not need to discard the numerical PDE knowledge encoded by the existing solver. A good block method tries to reuse it.

For a symmetric elliptic operator,

$$
A_{h}^{\mathsf T}=A_{h},
$$

so the state and adjoint blocks may share the same algebraic operator and often the same preconditioning machinery. Their right-hand sides and mathematical roles are still different.

## 4. Exact KKT factorization identifies the ideal block preconditioner

Consider the equality-constrained quadratic KKT matrix

$$
K
=
\begin{bmatrix}
Q & D^{\mathsf T}\\
D & 0
\end{bmatrix},
$$

and assume for this section that $Q$ is SPD and $D$ has full row rank. Then the Schur complement below is SPD as well.

Define the Schur complement

$$
S
:=
DQ^{-1}D^{\mathsf T}.
$$

Then

```math
K
=
\begin{bmatrix}
I & 0\\
DQ^{-1} & I
\end{bmatrix}
\begin{bmatrix}
Q & 0\\
0 & -S
\end{bmatrix}
\begin{bmatrix}
I & Q^{-1}D^{\mathsf T}\\
0 & I
\end{bmatrix}.
```

This factorization is exact. It immediately identifies two difficult inverse actions:

$$
Q^{-1},
\qquad
S^{-1}.
$$

An **ideal** block method would apply both exactly. That would make the outer system extremely easy, but it would usually defeat the purpose: exact application of the Schur-complement inverse can be as difficult as solving the original KKT system.

Practical preconditioners therefore replace the exact blocks by cheaper approximations:

$$
\widehat Q\approx Q,
\qquad
\widehat S\approx S.
$$

For example, a block-diagonal preconditioner may use

$$
P_{\mathrm{BD}}
=
\begin{bmatrix}
\widehat Q & 0\\
0 & \widehat S
\end{bmatrix},
$$

while a triangular construction imitates more of the exact factorization.

The derivation explains why block preconditioners are not arbitrary collections of subsolvers. They are approximations of an exact algebraic decomposition.

### 4.1 The Schur complement is normally applied indirectly

The exact Schur action is

$$
S\lambda
=
D
\left(
Q^{-1}
(D^{\mathsf T}\lambda)
\right).
$$

Even if $S$ is never assembled, its action can be described by three operations:

```text
multiplier vector λ
        ↓
apply D^T
        ↓
solve approximately with Q
        ↓
apply D
        ↓
Schur-complement action
```

This is the same operator-action viewpoint used throughout the main chapter.

A practical $\widehat S$ may be much simpler than this exact composition. Its design depends on the origin of the blocks. For mixed finite elements, mass matrices or pressure-space operators can approximate a Schur complement. For PDE-constrained optimization, approximations often combine mass-like blocks with state or adjoint solve actions.

### 4.2 Solver compatibility still matters

If the KKT matrix is symmetric and an outer MINRES method is desired, the preconditioning strategy must preserve the symmetry assumptions required by that method. In the usual SPD-preconditioned MINRES setting, the preconditioner itself must define an SPD action.

A block-triangular left preconditioner generally destroys ordinary Euclidean symmetry of the preconditioned matrix. GMRES can still be applicable because it does not require symmetry.

Thus an algebraically appealing factorization does not automatically determine the outer Krylov method. The representation of the preconditioned operator matters too.

## 5. Mesh robustness and parameter robustness are different goals

PDE-constrained optimization introduces parameters in addition to the mesh size. A common example is a regularization weight

$$
\beta>0.
$$

After discretization, a KKT system may therefore belong to a two-parameter family,

$$
K_{h,\beta}.
$$

A preconditioner can behave well as $h\to0$ but deteriorate as $\beta\to0$, or vice versa.

It is useful to keep two statements distinct:

- **mesh robust**: the relevant spectral or iteration bounds remain controlled as the mesh is refined;
- **parameter robust**: the bounds remain controlled as an optimization or physical parameter varies over the intended range.

A method can satisfy one without satisfying the other.

This distinction matters when comparing algorithms. An iteration table in which the count is almost constant under mesh refinement does not by itself show robustness with respect to regularization strength. Likewise, a preconditioner tuned for one $\beta$ may hide a severe mesh dependence.

For the prerequisite track, the important habit is simply to ask:

> With respect to which changing quantities is this solver claimed to be robust?

The detailed proof of parameter-robust bounds belongs to specialized numerical analysis, not to the main background sequence.

## 6. Inexact block solves can make the preconditioner variable

An exact block formula often contains inverse actions such as

$$
Q^{-1}r.
$$

In practice, that action may itself be approximated by an iterative solver stopped after a tolerance or a fixed number of iterations.

If the same fixed linear operation is applied every time, the outer method still sees a stationary preconditioner. But if the inner solve changes from one outer step to the next – for example because its stopping criterion depends on the current residual – then the effective preconditioner can vary with the outer iteration.

This matters because classical Krylov recurrences are derived for a fixed linear operator. Flexible variants, most notably flexible GMRES, are designed to tolerate a changing preconditioning action.

The practical lesson is not that every nested solve requires a flexible method. It is that the phrase “apply the preconditioner” may itself hide a numerical solve, and the behavior of that inner solve is part of the outer algorithm.

### 6.1 Inner solves should be accurate enough, not automatically exact

Suppose the outer iteration only needs a useful approximation of

$$
P^{-1}r.
$$

Driving the inner solve to nearly machine precision at every outer iteration can waste work. On the other hand, making the inner solve too crude can destroy the usefulness of the preconditioner.

This creates another accuracy–cost balance:

```text
more accurate inner solve
    ├── more expensive preconditioner application
    └── usually better outer progress

cheaper inner solve
    ├── less work per outer step
    └── potentially more outer iterations
```

Choosing that balance well is a solver-design problem. The background only needs the conceptual consequence: a solve action inside a formulation or preconditioner is a numerical operation with its own tolerance and evidence, not a symbolic inverse.

## 7. Multigrid and domain decomposition: the level of detail we need

Two scalable preconditioning families appear so often in finite-element and PDE-optimization work that it is useful to recognize their basic mechanism, even if we do not develop the algorithms in detail.

### 7.1 Multigrid attacks different error scales on different levels

A stationary relaxation such as Jacobi or Gauss–Seidel can reduce some error components efficiently while leaving slowly varying components largely untouched. Multigrid combines such a **smoother** with a hierarchy of coarser problems. The coarse levels represent error components that are expensive to remove on the fine grid but become easier to correct after restriction to a lower-dimensional space.

A typical multigrid cycle therefore has the schematic form

```text
fine-level residual
        ↓
pre-smoothing
        ↓
restrict residual to a coarser level
        ↓
coarse correction
        ↓
prolong correction to the fine level
        ↓
post-smoothing
```

The important point for this background track is not the exact V-cycle formula. It is that one or a few multigrid cycles can act as an inexpensive approximation to

$$
A_{h}^{-1}.
$$

That makes multigrid useful both as a standalone solver and, very often, as a preconditioner inside CG, MINRES, GMRES, or a block KKT preconditioner.

There are two common ways to obtain the hierarchy. **Geometric multigrid** uses the mesh hierarchy and finite-element transfer operators explicitly. **Algebraic multigrid (AMG)** constructs a hierarchy from the discrete matrix and related algebraic information. The distinction matters in implementation, but both serve the same high-level purpose here: provide an approximate inverse whose quality can remain controlled as the mesh is refined.

deal.II's tutorial sequence makes this role very explicit. Step-16 introduces multigrid as a preconditioner for the Laplace equation; step-50 compares geometric and algebraic multigrid on parallel adaptive meshes; and step-56 uses multigrid inside a block preconditioner for a saddle-point Stokes system. These are useful implementation continuations, not prerequisites for the main manual.

### 7.2 Domain decomposition and Schwarz methods solve local subproblems

A second family splits the discrete space or physical domain into smaller, often overlapping pieces and solves local problems on those pieces. In a Schwarz method, one may write schematically

$$
V_{h}
=
\sum_{j=1}^{J} V_{j},
$$

where each $V_{j}$ is a local subspace. The preconditioner combines approximate inverse actions on these subspaces. In an **additive** Schwarz method, local corrections can be computed independently and summed; in a **multiplicative** Schwarz method, later local solves can use corrections produced by earlier ones.

At this level, the connection with familiar methods is useful: point Jacobi can be interpreted as an additive method with one scalar unknown per subproblem, while Gauss–Seidel has the corresponding sequential character. Block and overlapping variants replace scalar subproblems by patches of degrees of freedom. deal.II's step-63 discusses this Schwarz viewpoint explicitly while constructing multigrid smoothers.

The phrase **domain decomposition** is also used in parallel finite-element software for partitioning the mesh among processes. That decomposition is related to scalable solver design, but it is not automatically the same thing as a Schwarz preconditioner. When reading documentation, distinguish

```text
mesh partitioning for ownership/load balance

from

subdomain solves used as a preconditioner.
```

deal.II's parallel tutorials use the first meaning extensively, while Schwarz-oriented solver components use the second.

### 7.3 Why these ideas matter for PDE-constrained optimization

For `nmopt`, the useful mental model is compositional. A state or adjoint block may be preconditioned by multigrid or a domain-decomposition method; that approximate PDE inverse can then appear inside a larger block or Schur-complement preconditioner for the optimization system.

Thus the hierarchy is often

```text
local / multilevel PDE preconditioner
        ↓
approximate state or adjoint inverse
        ↓
block preconditioner for a coupled KKT system
        ↓
outer Krylov method
```

This is enough background to understand why multigrid, Schwarz methods, AMG, and block preconditioners appear together in PDE-constrained-optimization papers and software. Their detailed construction is a numerical-PDE topic rather than a prerequisite for understanding the `nmopt` architecture or formulations.

## 8. Scope frontier

This note stops before the construction and analysis of particular scalable PDE preconditioners.

In particular, we do not develop:

- convergence proofs or implementation details for geometric or algebraic multigrid;
- convergence theory, coarse-space design, or overlap choices for domain decomposition and Schwarz methods;
- incomplete LU or incomplete Cholesky factorization;
- sparse-direct ordering and fill-reduction algorithms;
- augmented-Lagrangian or constraint preconditioners;
- detailed eigenvalue bounds for block KKT preconditioners;
- operator preconditioning in function-space norms;
- parallel and distributed-memory preconditioning.

These topics become important when solver scalability itself is the research or implementation problem. For the `nmopt` manual, it is enough to recognize their role and understand the block/operator structure they are designed to exploit.

## Where to go next

- **Return to the main route:** continue with [05 · Unconstrained numerical optimization](05-unconstrained-optimization.md); when equality constraints become central, [06 · Constrained optimization](06-constrained-optimization.md) and [07 · PDE-constrained optimization](07-pde-constrained-optimization.md) provide the optimization and PDE origins of the block systems discussed here.
- **Into the `nmopt` manual:** [manual 07 · Optimality systems and KKT](../concepts/07-optimality-systems-and-kkt.md) shows how the project represents KKT operators separately from the numerical solve policy, while [manual 08 · Complementarity and PDAS](../concepts/08-complementarity-and-pdas.md) repeatedly solves restricted KKT systems under changing active sets.
- **For implementation practice:** [09 · The deal.II finite-element workflow](09-dealii-finite-element-workflow.md) and the deal.II tutorial references below show where multigrid, Schwarz, and block preconditioners enter concrete solver code.

## References and further reading

Publication details, stable links, and access notes are collected in the [Background reference catalogue](references.md).

- **Manzoni, Quarteroni, and Salsa, *Optimal Control of Partial Differential Equations*.**
  - Appendix B, §B.1.3 for saddle-point spectra, Schur complements, and ideal block factorizations.
  - Chapter 6, §6.6 for block-diagonal and block-triangular preconditioners in PDE-constrained optimization, including mesh and regularization-parameter robustness.
- **Michele Benzi, Gene H. Golub, and Jörg Liesen, “Numerical solution of saddle point problems.”**
  - Broad survey of saddle-point structure, Krylov solvers, Schur complements, and block preconditioning.
- **Yousef Saad, *Iterative Methods for Sparse Linear Systems*, second edition.**
  - Chapters 9–10 for preconditioned Krylov methods and concrete preconditioners.
  - Chapter 13 for multigrid.
  - Chapter 14 for domain decomposition.
- **Susanne C. Brenner and L. Ridgway Scott, *The Mathematical Theory of Finite Element Methods*.**
  - **Finite Element Multigrid Methods** and **Additive Schwarz Preconditioners** provide the finite-element route into the scalable PDE preconditioners introduced here.
- **Daniele Boffi, Franco Brezzi, and Michel Fortin, *Mixed Finite Element Methods and Applications*.**
  - **Algebraic Aspects of Saddle Point Problems** connects mixed finite elements with the block linear algebra used here.
- **deal.II tutorial programs.**
  - step-16 for geometric multigrid preconditioning of Laplace.
  - step-50 for geometric versus algebraic multigrid.
  - step-63 for additive and multiplicative Schwarz smoothers in multigrid.
  - step-31 and step-56 for Schur-complement and multigrid ideas inside block saddle-point solvers.
- **University of Pisa, *Scientific Computing / Calcolo Scientifico*.**
  - §6.8 for GMRES preconditioning inside the broader Chapter 6 Krylov framework.
  - The bibliography points to Demmel, Saad, and Benzi–Golub–Liesen for deeper treatments.
