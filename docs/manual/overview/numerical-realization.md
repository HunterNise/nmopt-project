# Numerical realization

A PDE-constrained optimization algorithm does not operate on a symbolic PDE. It
operates on a finite-dimensional numerical realization: vectors, operators, solves,
coordinate maps, and objective evaluations.

This layer is where mathematical choices become deal.II computations.

It is also the point where many apparently small implementation decisions acquire
mathematical meaning. Which degrees of freedom are independent? How are fixed
boundary values represented? Does a control coefficient live on cells, nodes, or
boundary faces? Is the optimization gradient measured in coefficient Euclidean norm
or in an $L^{2}$ metric?

The purpose of this page is to make those choices visible without descending into
the compiler's implementation details.

## 1. From continuous spaces to discrete coordinates

Suppose the continuous optimization problem uses a state space $Y$ and control space
$U$. A finite-element implementation replaces them with discrete spaces
$`Y_{h}`$ and $`U_{h}`$ and then represents fields by coefficient vectors.

That sounds straightforward, but the coefficient vector stored in memory is not
always the physical field seen by the PDE.

A common case is a state with prescribed Dirichlet data. Only the free degrees of
freedom are independent. The physical field can be written as

$$
y_{\mathrm{phys}} = Pz+\ell,
$$

where

- $z$ contains the independent state coefficients;
- $P$ injects them into the full finite-element vector;
- $\ell$ carries the fixed boundary values.

A perturbation uses only the homogeneous part,

$$
\delta y_{\mathrm{phys}}=P \delta z.
$$

This distinction appears in adjoints and derivatives, not just in output.

## 2. Example: adapted deal.II Step-4

`step-4` is a deal.II tutorial program for a small finite-element Poisson problem.
The repository keeps that application as the native PDE solver and builds two
optimal-control examples around it. The richer distributed-control example is called
**Problem B** in the integration study.

Problem B makes the coordinate issue concrete.

The control is represented in the full continuous finite-element basis with 289
coefficients. The state has the same physical mesh, but prescribed boundary values
leave 225 independent state coefficients.

The control enters the weak equation as a distributed forcing. Let $M$ be the
finite-element mass matrix. Restricting the weak equations to the free state
coordinates gives a rectangular coupling

$$
B=P^{\mathsf T}M.
$$

The discrete state equation can then be written schematically as

$$
Kz=b_{F}+Bu,
$$

with $K$ the free-state operator and $`b_{F}`$ the boundary-treated load.

The objective used in this integration is

```math
J(z,u)
=
\frac{1}{2}(Pz+\ell)^{\mathsf T}M(Pz+\ell)
+
\frac{1}{2}u^{\mathsf T}Mu.
```

Several project concepts are visible in this one example:

- state and control can have different dimensions;
- physical state reconstruction is not the same as the optimization coordinate;
- the control-to-state coupling can be rectangular;
- the same mass matrix may participate in several mathematical roles.

These distinctions belong to the numerical realization, not to the generic
optimizer.

## 3. Residual evaluation and solving are different operations

For a discrete problem, a residual evaluator answers:

$$
\text{given }(y,u),\quad \text{what is }E(y,u)?
$$

A state solver answers a different question:

$$
\text{given }u,\quad \text{find }y\text{ such that }E(y,u)=0.
$$

The first may be a matrix-vector operation. The second may involve a sparse direct
solve, CG, GMRES, nonlinear iteration, preconditioning, tolerances, and convergence
evidence.

`nmopt` keeps these responsibilities separate because a formulation may need both.
For example:

- derivative verification may evaluate residual and Jacobian actions;
- reduced optimization needs a state solve;
- adjoint construction needs a transpose solve;
- solver reports are useful evidence and should not be confused with residual
  values.

An external application is therefore free to retain its native solve policy while
exposing only the solve operation required by the formulation.

## 4. The executable mathematical operations

At the common numerical boundary, the formulation works with a small family of
operations.

Let $X$ denote the discrete primal variable space – for a state-control problem it
contains the coefficient representation of $(y,u)$ – and let $Z$ denote the
discrete residual test space. Their dual spaces are $X^{\ast}$ and $Z^{\ast}$.
For a point $x\in X$, the relevant operations are conceptually:

```math
\begin{aligned}
E(x) &\in Z^{\ast}, \\
E'(x) \delta x &\in Z^{\ast}, \\
E'(x)^{\ast}p &\in X^{\ast}, \\
J(x) &\in \mathbb{R}, \\
J'(x) &\in X^{\ast}.
\end{aligned}
```

The Jacobian-vector product propagates a primal perturbation forward. The adjoint
Jacobian action pulls a test/adjoint variable back to covectors on the variables.

The formulation does not need to know whether these operations came from assembled
matrices, matrix-free evaluation, callbacks into an application, or compiler-created
deal.II services.

## 5. The coefficient pairing is not the optimization metric

If a covector $r\in U^{\ast}$ acts on a perturbation $\delta u\in U$, the directional
change is

$$
\langle r,\delta u\rangle.
$$

In a coefficient representation this may be the ordinary coefficient dot product
between the stored covector and primal arrays. That pairing expresses the dual
action.

A metric is different. It is an operator

$$
G:U\rightarrow U^{\ast}
$$

used to measure the control space and to convert a derivative into a primal gradient:

$$
Gg=r.
$$

For an $L^{2}$ finite-element geometry, $G$ is naturally represented by a mass
matrix. For a coefficient-space Euclidean geometry, it may be the identity.

This separation is important because the optimizer should not silently insert a mass
matrix into every pairing. The objective, derivative representation, and chosen
optimization geometry are separate mathematical decisions.

## 6. Regularization and metric may use the same operator without being the same thing

Consider

$$
J(y,u)
=
J_{\mathrm{state}}(y)
+
\frac{\beta}{2}u^{\mathsf T}Mu.
$$

The control derivative contains the objective term

$$
J_{u}'=\beta Mu.
$$

If the optimization metric is also $G=M$, the corresponding gradient contribution
is

$$
G^{-1}J_{u}'=\beta u.
$$

The fact that $M$ appears in both places does not mean the regularization and metric
are one subsystem.

- Regularization changes **what objective is minimized**.
- The metric changes **how a covector is represented as a primal gradient and how
  steps are measured**.

Keeping those roles distinct makes algorithmic experiments meaningful.

## 7. Control realization can change the algebraic problem

A “distributed control” is not one unique discrete object.

It may be represented by:

- cellwise constants;
- continuous nodal finite elements;
- another independent finite-element space.

A boundary control may similarly be:

- one value per face;
- a continuous trace field;
- some other boundary discretization.

Changing the control realization changes:

- the control dimension;
- the coupling operator;
- the metric;
- possibly the admissible constraint representation.

This is why control discretization is an explicit part of the semantic/compiler
request and of the Chapter 6 scenario records.

## 8. Observations are numerical operators too

Objectives often compare the state with data only on part of the domain or through
a derived quantity.

A volume observation over a subregion $\omega\subset\Omega$ may contribute

$$
\frac{1}{2}\int_{\omega}(y-y_{\mathrm{d}})^{2} \mathrm{d}x.
$$

Its discrete realization needs to know:

- which cells belong to the observation region;
- the quadrature rule;
- how the desired state is evaluated;
- how the state derivative is accumulated.

These decisions are neither generic optimizer behavior nor incidental I/O. They are
part of the numerical realization of the objective.

The repository's **B2** scenario, derived from the Chapter 6 Graetz-flow
boundary-control example in Manzoni–Quarteroni–Salsa, makes this visible by
selecting a downstream material region and an explicit observation evaluation
policy.

## 9. Boundary conditions are part of the coordinate model

Fixed Dirichlet data is the clearest case, but boundary controls add another layer.

A boundary-control problem must distinguish:

- the part of the boundary with fixed state data;
- the controlled boundary;
- outflow or natural boundary portions;
- the sign and normal/conormal convention of the weak form.

For the current B2 problem, the source-oriented ordinary-normal convention is

$$
\partial_{n} y-(b\cdot n)y=u
$$

on the controlled boundary.

A different conormal convention is not a harmless code refactor; it defines a
different mathematical realization. The compiler therefore records such selections
explicitly rather than inferring them from an implementation branch.

## 10. State and adjoint solves reuse numerical structure

For many linear-quadratic problems, the adjoint system uses the transpose of the
state Jacobian.

If the state equation is

$$
A(u)y=f(u),
$$

the state solve uses the state operator, while the adjoint relation has the form

$$
A(u)^{\mathsf T}p=J_{y}'.
$$

The implementation may reuse matrix structure, constraints, sparsity, or solver
configuration, but the right-hand side and mathematical role are different.

For a symmetric elliptic operator such as the Step-4 example, the same stiffness
matrix can serve both systems. The formulation still treats state and adjoint solves
as separate services because symmetry is not a universal assumption.

## 11. Constraints live in the optimization coordinates

A box constraint such as

$$
u_{\min}\leq u\leq u_{\max}
$$

must be interpreted in the chosen discrete control representation.

For cellwise-discontinuous control, coefficientwise bounds naturally correspond to
cellwise bounds. For a continuous nodal control, coefficient bounds have a different
physical meaning.

This is why primal-dual active-set (PDAS)/complementarity products and projected
reduced methods depend on a
well-defined control layout and metric rather than on an abstract bounded-control flag.

## 12. Backend policy versus numerical realization

The generic solver needs primitive vector operations:

```text
copy      scale      axpy      dot      norm-related primitives
```

Those are backend-policy concerns.

The finite-element realization needs:

```text
mesh
FE spaces
coordinate maps
assembled operators
boundary treatment
observations
state/adjoint solves
metrics
constraints
native output
```

Those are numerical-realization concerns.

The serial deal.II backend happens to use the same native vector type as the
finite-element realization, which makes integration efficient, but the conceptual
boundary should remain clear.

## 13. Compiler-owned and application-owned realization

The same numerical responsibilities can be owned in two ways.

### Compiler-owned path

The semantic/compiler stack chooses and constructs supported deal.II realization
components, then packages them into a compiled problem.

### Application-owned path

An existing numerical application keeps those components and supplies callable
operations to the formulation.

The reduced optimizer does not care which owner produced them.

This is one of the most useful architectural facts for a new contributor: when
reading a numerical service, ask **what mathematical operation it realizes**, not
only which producer path constructed it.

## 14. Native output remains outside the optimizer

A finite-element result usually needs more than an optimization vector.

Writing a final state may require:

- reconstructing fixed boundary values;
- mapping independent coefficients into a physical field;
- attaching control values as cell or point data;
- using the application's `DoFHandler` and triangulation.

Those objects should remain available in the application/native view.

The optimizer only needs the retained numerical state associated with the final
control. The outer application then reconstructs and writes the physical result.

## 15. Read later: authoritative sources

For the detailed numerical boundary, use
[PDE, formulation, and solver boundary](../../design/pde-solver-boundary.md).

For the direct-application path, first read
[Integrating an existing PDE application](external-applications.md). For the concrete
Step-4 realization, use:

- [Step-4 integration overview](../../../apps/external-dealii/step-4/external-integration-overview.md)
  for the worked case-study narrative;
- [External deal.II solver integration](../../reference/external-dealii-solver-integration.md)
  for exact public contracts.

For compiler-owned realizations, the
[Compiler reference](../../reference/compiler.md) records the current public
capability boundary, while
[Compiler implementation](../../internals/compiler.md) explains the lowering
and realization mechanics.

For implementation source, the main areas are:

```text
include/nmopt/dealii/       reusable deal.II numerical services
include/nmopt/compiler/v1/  construction of compiled realizations
apps/external-dealii/       application-owned integration examples
```

Return to [Project overview and architecture](project-architecture.md) for the
whole-system picture or continue with
[Reduced state–adjoint optimization](reduced-optimization.md) to see how these
operations are consumed.
