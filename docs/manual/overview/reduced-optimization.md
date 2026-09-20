# Reduced state–adjoint optimization

Reduced optimization treats the control as the variable seen by the optimizer and
eliminates the PDE state through a state solve.

For PDE-constrained optimization this is a natural approach when a reliable state
solver already exists: the optimizer works in the control space, while the
formulation repeatedly asks the PDE application for state and adjoint solves.

This page follows one reduced evaluation and one search step in enough detail to make
the runtime concrete.

## 1. From constrained problem to reduced objective

Start from

$$
\min_{y,u} J(y,u)
\qquad\text{subject to}\qquad
E(y,u)=0.
$$

Here $y$ is the state, $u$ the control, $J$ the objective, and $E$ the discrete
state equation or residual. Assume that, for controls of interest, solving the PDE
defines a state map $S$ with $y=S(u)$.

The reduced objective is

$$
j(u)=J(S(u),u).
$$

A naive derivative calculation could differentiate the state solution with respect
to every control direction. For a control with many finite-element degrees of
freedom, that is prohibitively expensive.

The adjoint method reorganizes the calculation so that one adjoint solve produces
the full reduced derivative.

## 2. The adjoint relation

Use the Lagrangian

$$
\mathcal{L}(y,u,p)
=
J(y,u)-\langle p,E(y,u)\rangle.
$$

At a state satisfying $E(y,u)=0$, choose $p$ from

$$
E_{y}'(y,u)^{\ast}p
=
J_{y}'(y,u).
$$

Then the reduced derivative is

$$
j'(u)
=
J_{u}'(y,u)-E_{u}'(y,u)^{\ast}p.
$$

The sign follows from the Lagrangian convention above.

The important computational point is the cost structure:

```text
one control u
   │
   ├─ 1 state solve
   │
   └─ 1 adjoint solve
        │
        ▼
full reduced derivative in the control space
```

The number of control coefficients does not imply one PDE solve per coefficient.

## 3. What one value evaluation does

Suppose the optimizer asks for the value at a control $`u_{k}`$.

The reduced formulation:

1. calls the state solver to obtain $`y_{k}`$;
2. records the state-solve convergence evidence;
3. evaluates $`J(y_{k},u_{k})`$;
4. retains the state together with the control and objective value.

Conceptually:

```text
u_k
 │
 ├─► solve E(y_k,u_k)=0
 │
 ▼
(y_k,u_k)
 │
 └─► J(y_k,u_k)
       │
       ▼
 retained value object
```

Retaining the state matters. If the optimizer later needs a derivative at the same
control, the formulation should not solve the state equation again.

## 4. Derivative augmentation reuses the retained state

Given a retained value at $`(y_{k},u_{k})`$, derivative augmentation proceeds as

1. evaluate the state and control parts of $`J'(y_{k},u_{k})`$;
2. solve the adjoint equation for $`p_{k}`$;
3. evaluate the residual transpose action;
4. combine the control covectors into $`j'(u_{k})`$.

In symbols,

```math
\begin{aligned}
E_{y}'(y_{k},u_{k})^{\ast}p_{k}
&=
J_{y}'(y_{k},u_{k}),\\
j'(u_{k})
&=
J_{u}'(y_{k},u_{k})
-
E_{u}'(y_{k},u_{k})^{\ast}p_{k}.
\end{aligned}
```

The returned derivative remains a covector.

## 5. Example: the Step-4 distributed-control case

Here **Step-4** means the adapted deal.II tutorial application used for the
external-integration study, and Problem B is the distributed-control problem built
around it.

For that example,

$$
Kz=b_{F}+Bu
$$

and

```math
J(z,u)
=
\frac{1}{2}(Pz+\ell)^{\mathsf T}M(Pz+\ell)
+
\frac{1}{2}u^{\mathsf T}Mu.
```

The objective derivatives are

```math
\begin{aligned}
J_{z}
&=
P^{\mathsf T}M(Pz+\ell),\\
J_{u}
&=
Mu.
\end{aligned}
```

The adjoint solves

$$
K^{\mathsf T}p=J_{z}.
$$

For this residual convention, the control part of the residual transpose action is

$$
E_{u}'(z,u)^{\ast}p=-B^{\mathsf T}p.
$$

Therefore

$$
j'(u)=Mu+B^{\mathsf T}p.
$$

This is the exact kind of algebra the reduced formulation exists to organize. The
application supplies $K$, $B$, $M$, reconstruction, and solves. The formulation
supplies the state–adjoint recipe.

## 6. A derivative is not yet a search direction

The reduced derivative $j'(u)$ belongs to the dual control space.

Choose a metric

$$
G:U\rightarrow U^{\ast}.
$$

The corresponding primal gradient solves

$$
Gg=j'(u).
$$

For this Step-4 Problem B, $G=M$, so

$$
Mg=j'(u).
$$

A steepest-descent direction is then

$$
d=-g.
$$

Other direction policies, such as quasi-Newton methods, use the derivative and metric
history differently. The formulation does not choose the algorithm.

## 7. Why value and derivative evaluation are separate

Consider an Armijo line search from $`u_{k}`$ along a direction $`d_{k}`$.

A trial control is

$$
u_{\mathrm{trial}}=u_{k}+\alpha d_{k}.
$$

To decide whether the trial has sufficient decrease, the algorithm needs its reduced
objective value. That requires

```text
trial control
   │
   ├─ state solve
   ▼
trial state
   │
   └─ objective value
```

It does **not** yet need an adjoint.

If the trial is rejected, the algorithm shortens $\alpha$ and tries another state
solve. Only after a trial is accepted does it augment that retained state with an
adjoint and derivative.

For expensive PDE solves this separation is significant:

```text
rejected trial  → state solve + objective
accepted point  → retained state + adjoint + reduced derivative
```

The project models this distinction directly rather than presenting “evaluate
objective and gradient” as one inseparable callback.

## 8. The Armijo condition uses the covector pairing

A sufficient-decrease test has the form

$$
j(u_{k}+\alpha d_{k})
\leq
j(u_{k})
+
c \alpha\langle j'(u_{k}),d_{k}\rangle.
$$

Equivalently, using the actual trial displacement

$$
\delta u=u_{\mathrm{trial}}-u_{k},
$$

the directional model uses

$$
\langle j'(u_{k}),\delta u\rangle.
$$

This is a dual pairing. It should not be replaced by a metric norm or by an
arbitrary coefficient-space interpretation.

The metric enters when constructing gradients, directions, and metric-dependent
norms. The derivative pairing expresses the first-order change of the objective.

## 9. What the reduced optimizer controls

Once the formulation can provide values and derivatives, the reduced solver owns the
optimization logic.

Depending on the selected policy, that includes:

- steepest-descent or quasi-Newton directions;
- line-search trial generation;
- Armijo acceptance;
- trust-region logic for supported paths;
- projection for supported constraints;
- stopping criteria;
- iteration and work records.

The solver sees a numerical optimization problem, not a PDE.

A simplified steepest-descent loop is

```text
evaluate and differentiate u₀

repeat
    convert derivative to primal gradient
    choose direction
    propose trial control
    evaluate trial value
    accept/reject according to globalization
    if accepted:
        augment derivative at retained trial state
    test stopping criteria
```

This pseudocode is intentionally formulation-level. The state and adjoint solves are
hidden behind the reduced evaluation service.

## 10. Quasi-Newton methods change the direction model, not the PDE contract

L-BFGS or full BFGS approximates inverse-curvature information from accepted
iterations.

The PDE realization does not need a new state equation merely because the direction
policy changes. It still supplies reduced objective values and derivatives through
the same state–adjoint formulation.

That separation makes algorithm comparisons meaningful: the numerical problem can
remain fixed while the direction and globalization policies change.

One repository benchmark, B1 – the distributed Laplace-control example derived
from Chapter 6 of Manzoni–Quarteroni–Salsa – uses this idea to compare steepest
descent with limited-memory BFGS on the same compiled problem family. It is a
demonstration of the solver boundary rather than a definition of the framework's
scope.

## 11. Newton-type methods require more information

A genuine Newton direction requires a reduced Hessian action or an equivalent
second-order product.

That capability is not inferred automatically from the five first-order executable
operations.

The project therefore treats a supplied reduced Hessian as an additional capability.
Selecting a Newton-like method may impose stronger requirements on the numerical
problem than selecting steepest descent or BFGS.

This is another example of the general rule: **algorithm choice may request extra
mathematical services, but it should not alter the meaning of the existing ones.**

## 12. Constraints and projection

For an admissible control set $`U_{\mathrm{ad}}`$, a reduced method may need to project
a trial point:

$$
u_{\mathrm{trial}}
=
\Pi_{U_{\mathrm{ad}}}^{G}(u+\alpha d).
$$

The superscript emphasizes that projection may depend on the selected metric.

The project keeps constraints and metrics separate from the objective so that the
same objective can be optimized under different admissible sets or geometries.

For strongly structured box-constrained problems, the project also has
complementarity and primal-dual active-set (PDAS) products rather than forcing every
constrained formulation
through a projected reduced method.

## 13. Stopping is measured in the chosen mathematical geometry

A solver may stop based on quantities such as:

- gradient norm;
- step norm;
- objective decrease;
- iteration or work limits.

When a gradient is defined through $Gg=j'(u)$, its natural norm depends on the metric:

$$
\lVert g\rVert_{G}^{2}
=
\langle Gg,g\rangle
=
\langle j'(u),g\rangle.
$$

This is why the metric is a first-class input to the solver rather than an internal
implementation detail of the PDE application.

## 14. The same reduced solver can receive two different producer paths

The producer paths converge before the reduced solver:

```text
ProblemSpec ──► compiler ───────────┐
                                    │
existing application ──► binding ───┤
                                    ▼
                         reduced formulation
                                    │
                                    ▼
                           reduced solver
```

The reduced solver does not branch on which path produced the formulation.

This convergence is one of the main tests of the architecture: problem construction
and optimization orchestration are reusable independently.

## 15. What the optimizer returns

A useful solver result needs more than a final control vector.

The current architecture can retain:

- the final reduced evaluation;
- the state associated with the accepted final control;
- solver status and iteration information;
- work/evaluation evidence.

That retained state allows the outer application to write the final physical field
without solving the PDE once more solely for output.

The research layer can then associate the result with compilation and environment
provenance.

## 16. Current reduced-formulation scope

The present `ReducedDTOT` path is intentionally specialized to the common
one-state/one-decision structure:

```text
variables:  state + control
residual:   one test block
product:    state-eliminated reduced objective
```

The mathematical interfaces elsewhere in the repository are broader, but a reader
should not assume that arbitrary multi-state or multi-control graph structures are
already supported by this particular runtime.

## 17. Read later: authoritative sources

For the mathematical convention behind the adjoint sign, reduced derivative, and
metric distinction, use
[Theoretical formalism](../../design/mathematical-model.md).

For the ownership boundary between state/adjoint solves and optimization, use
[PDE, formulation, and solver boundary](../../design/pde-solver-boundary.md).

For a concrete reduced integration with an existing deal.II code, use the
[Step-4 integration overview](../../../apps/external-dealii/step-4/external-integration-overview.md)
and the
[external integration reference](../../reference/external-dealii-solver-integration.md).

For exact current C++ contract signatures, inspect the public headers under
`include/nmopt/contract/` and `include/nmopt/solvers/` together with their focused
tests. The existing
[interface specification](../../design/interface-specification.md) provides the
normative design vocabulary but should not be read as a substitute for implemented
signatures.

To understand how the numerical objects used here are built, read
[Numerical realization](numerical-realization.md).
