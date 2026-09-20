# Reduced state–adjoint formulation

## Eliminating the state without losing its derivative

The first four concept chapters developed the numerical language shared by several
formulations in `nmopt`. We can now use that language to follow one complete formulation
from its mathematical motivation to its runtime structure.

The reduced state–adjoint formulation begins from a constrained problem

```math
\begin{aligned}
\min_{y,u}\quad & J_{h}(y,u),\\
\text{subject to}\quad & E_{h}(y,u)=0.
\end{aligned}
```

Here

$$
y\in Y_{h},
\qquad
u\in U_{h},
\qquad
E_{h}(y,u)\in Z_{h}^{\ast}.
$$

The central idea is to treat the control $u$ as the optimization variable and the state
$y$ as something determined by the PDE. If the state equation can be solved for a given
control, we may write

$$
y=S_{h}(u),
$$

where $`S_{h}`$ is the discrete **control-to-state map**. The original constrained problem
then becomes the unconstrained reduced problem

$$
\min_{u} j_{h}(u),
\qquad
j_{h}(u)
:=
J_{h}(S_{h}(u),u).
$$

The state has disappeared from the list of optimization variables, but it has not
disappeared computationally. Every evaluation of $`j_{h}(u)`$ still requires a state solve.
Every derivative of $`j_{h}`$ must still account for how that state changes with the
control.

That is where the adjoint enters.

This chapter follows that argument from the equations to the current `ReducedDTOT`
implementation. The optimization algorithm itself—line searches, BFGS, Newton
directions, trust regions, stopping, and work accounting—is developed in Chapter 6.
Here we stop once the reduced objective and reduced derivative have been constructed.

The main path is:

```text
control u
    │
    ▼
state solve
E_h(y,u) = 0
    │
    ▼
feasible point (y,u)
    │
    ├──────────────► objective value j_h(u)
    │
    └─ if derivative is needed:
            │
            ▼
        J_h'(y,u)
            │
            ▼
        adjoint solve
            │
            ▼
        E_h'(y,u)* p
            │
            ▼
        reduced derivative j_h'(u)
```

The split between the left and right branches will become important later: a trial
control can often be rejected using only its state and objective value, without paying
for an adjoint solve.

## 1. What does it mean to eliminate the state?

The notation

$$
y=S_{h}(u)
$$

compresses a numerical solve into one mathematical symbol. For the linear
distributed-control problem from the first chapter,

$$
E_{h}(y,u)
:=
Ay-Bu-f,
$$

the state equation is

$$
Ay=f+Bu.
$$

Assuming $A$ is invertible,

$$
S_{h}(u)
=
A^{-1}(f+Bu).
$$

Writing the inverse explicitly is useful for derivations, but a real implementation does
not normally construct $A^{-1}$. It solves the linear system.

For a nonlinear PDE,

$$
E_{h}(y,u)=0
$$

may require a Newton method or some other nonlinear state solver. The notation $`S_{h}(u)`$
still makes sense locally if the state problem has a well-defined solution, but the
computational object behind it is now a solve procedure rather than a single matrix
inverse.

This distinction is reflected directly in the formulation interface. The reduced
formulation is not given a matrix called $`S_{h}`$. It is given a callable **state solve
service**

```text
control u
    ↓
solve_state
    ↓
state y
```

together with the executable residual/objective model.

### 1.1 The reduced formulation assumes a state can be produced for each admissible control

The reduction only makes sense when the state equation can actually be solved, at least
for the controls visited by the algorithm. Analytically, this may follow from
well-posedness or a local implicit-function argument. Numerically, it also depends on
the chosen discretization and solve policy.

`ReducedDTOT` does not prove those facts. It receives a state solver from the numerical
realization and requires that the solve report says the state solve converged.

So the symbol $`S_{h}`$ should be read as

> the state selected by the supplied numerical state-solve service for this control,

not as a symbolic inverse owned by the reduced formulation.

## 2. Differentiating the reduced objective directly

The reduced objective is

$$
j_{h}(u)
:=
J_{h}(S_{h}(u),u).
$$

Take a control perturbation

$$
\delta u\in U_{h}.
$$

The state changes by

$$
\delta y
=
S_{h}'(u)[\delta u].
$$

By the chain rule,

```math
j_{h}'(u)[\delta u]
=
D_{y} J_{h}(y,u)[\delta y]
+
D_{u} J_{h}(y,u)[\delta u],
```

where

$$
y=S_{h}(u).
$$

The difficulty is the first term: we need the state sensitivity $\delta y$.
Differentiate the state equation

$$
E_{h}(y,u)=0.
$$

Because $y$ depends on $u$,

```math
D_{y} E_{h}(y,u)[\delta y]
+
D_{u} E_{h}(y,u)[\delta u]
=
0.
```

Therefore the sensitivity satisfies

$$
D_{y} E_{h}(y,u)[\delta y]
=
-
D_{u} E_{h}(y,u)[\delta u].
$$

If the state Jacobian can be inverted, then formally

$$
\delta y
=
-
\left(D_{y} E_{h}(y,u)\right)^{-1}
D_{u} E_{h}(y,u)[\delta u].
$$

Substituting into the chain rule gives

```math
\begin{aligned}
j_{h}'(u)[\delta u]
&=
-
D_{y} J_{h}(y,u)
\left[
\left(D_{y} E_{h}(y,u)\right)^{-1}
D_{u} E_{h}(y,u)[\delta u]
\right]
\\
&\quad
+
D_{u} J_{h}(y,u)[\delta u].
\end{aligned}
```

This formula is correct, but it is a poor way to construct the whole reduced derivative
when the control dimension is large.

### 2.1 The direct sensitivity cost

For one chosen direction $\delta u$, the direct method can be entirely reasonable:

1. apply $`D_{u} E_{h}`$ to $\delta u$;
2. solve one linearized state equation for $\delta y$;
3. evaluate $`D_{y}J_{h}[\delta y]`$.

The problem appears when we want a representation of

$$
j_{h}'(u)\in U_{h}^{\ast}
$$

that can act on **every** control direction. 

Let $`\{\psi_{1},\ldots,\psi_{n_{u}}\}`$ be a
basis of $`U_{h}`$. A basis-by-basis sensitivity construction would require

$$
D_{y}E_{h}(y,u)[\delta y_{k}]
=
-
D_{u}E_{h}(y,u)[\psi_{k}],
\qquad
k=1,\ldots,n_{u}.
$$

That is one linearized state solve per control basis direction. For PDE controls, $`n_{u}`$
can be large.

The adjoint method reorganizes the same chain rule so that the state-side inverse is
applied once instead, independently of the number of control directions.

## 3. Introduce the adjoint through the pairing identity

The partial state derivative of the objective is a covector

$$
D_{y}J_{h}(y,u)\in Y_{h}^{\ast}.
$$

The state derivative of the residual is

$$
D_{y}E_{h}(y,u):
Y_{h}
\longrightarrow
Z_{h}^{\ast}.
$$

Its transpose is

$$
D_{y}E_{h}(y,u)^{\ast}:
Z_{h}
\longrightarrow
Y_{h}^{\ast}.
$$

Choose the adjoint $`p\in Z_{h}`$ to satisfy

$$
D_{y}E_{h}(y,u)^{\ast}p
=
D_{y}J_{h}(y,u).
$$

This is the **adjoint equation** used by the current reduced formulation. Now apply the
transpose identity to the state-sensitivity term:

```math
\begin{aligned}
D_{y}J_{h}(y,u)[\delta y]
&=
\left\langle
D_{y}E_{h}(y,u)^{\ast}p,
\delta y
\right\rangle
\\
&=
\left\langle
D_{y}E_{h}(y,u)[\delta y],
p
\right\rangle.
\end{aligned}
```

The linearized state equation tells us that

$$
D_{y}E_{h}(y,u)[\delta y]
=
-
D_{u}E_{h}(y,u)[\delta u].
$$

Therefore

```math
\begin{aligned}
D_{y}J_{h}(y,u)[\delta y]
&=
-
\left\langle
D_{u}E_{h}(y,u)[\delta u],
p
\right\rangle
\\
&=
-
\left\langle
D_{u}E_{h}(y,u)^{\ast}p,
\delta u
\right\rangle.
\end{aligned}
```

Substitute this into the chain rule:

```math
j_{h}'(u)[\delta u]
=
\left\langle
D_{u}J_{h}(y,u)
-
D_{u}E_{h}(y,u)^{\ast}p,
\delta u
\right\rangle.
```

Hence

$$
\boxed{
j_{h}'(u)
=
D_{u}J_{h}(y,u)
-
D_{u}E_{h}(y,u)^{\ast}p
}
$$

with

$$
\boxed{
D_{y}E_{h}(y,u)^{\ast}p
=
D_{y}J_{h}(y,u).
}
$$

The state sensitivity $\delta y$ has disappeared from the final formula. One adjoint
solve produces a test-space vector $p$ whose transpose residual action contains the
state-dependence contribution for all control directions at once.

## 4. The sign comes from the residual convention

For the running linear residual,

$$
E_{h}(y,u)
:=
Ay-Bu-f,
$$

we have

$$
D_{y}E_{h}=A,
\qquad
D_{u}E_{h}=-B.
$$

The adjoint equation is

$$
A^{\mathsf T}p
=
D_{y}J_{h}.
$$

The control part of the residual transpose is

$$
D_{u}E_{h}^{\ast}p
=
-B^{\mathsf T}p.
$$

Therefore

```math
\begin{aligned}
j_{h}'(u)
&=
D_{u}J_{h}
-
\left(-B^{\mathsf T}p\right)
\\
&=
D_{u}J_{h}+B^{\mathsf T}p.
\end{aligned}
```

This is the same formula derived in [Anatomy of a discrete PDE-constrained
problem](01-discrete-problem-anatomy.md), now written in the general derivative/VJP
language of [Operators, derivatives, and
adjoints](03-operators-derivatives-and-adjoints.md).

There is no independent "adjoint sign" to memorize. Once the residual convention and
Lagrangian convention are fixed, the sign follows from the derivative.

### 4.1 The project's Lagrangian convention

The project writes

$$
\mathcal L(y,u,p)
:=
J_{h}(y,u)
-
\langle p,E_{h}(y,u)\rangle.
$$

Stationarity with respect to the state gives

$$
D_{y}J_{h}-D_{y}E_{h}^{\ast}p=0,
$$

or

$$
D_{y}E_{h}^{\ast}p=D_{y}J_{h}.
$$

Stationarity with respect to the control gives

$$
D_{u}J_{h}-D_{u}E_{h}^{\ast}p.
$$

So the reduced-derivative formula above is also the control derivative of the Lagrangian
evaluated at a feasible state and its adjoint. This equivalence is useful because it
connects the reduced formulation to the all-at-once optimality-system viewpoint
developed later.

## 5. The formulation needs solves as well as operator actions

Part I separated four notions:

```text
evaluate E_h(y,u)
apply D E_h(y,u)
apply D E_h(y,u)*

solve an equation involving those operators
```

The reduced state–adjoint formulation uses both categories. The executable model
supplies:

- $`J_{h}(y,u)`$;
- $`J_{h}'(y,u)`$;
- $`E_{h}'(y,u)^{\ast}p`$ through the VJP.

The numerical realization separately supplies:

- a state solver for $`E_{h}(y,u)=0`$;
- an adjoint solver for
  $`D_{y}E_{h}(y,u)^{\ast}p=D_{y}J_{h}(y,u)`$.

This separation is visible in `StateAdjointSolversT`. Conceptually its two callbacks are

```text
solve_state
    u
    ↓
    y satisfying E_h(y,u)=0

solve_adjoint
    (y,u), D_y J_h
    ↓
    p satisfying D_y E_h(y,u)* p = D_y J_h
```

The formulation does not prescribe whether the state solve uses CG, Newton, a direct
factorization, or an application-specific native solver. It asks for the result and
enough solve evidence to know whether the required solve succeeded.

### 5.1 Why the adjoint solver receives the full point

For a linear state equation, the transpose state operator may be independent of the
current state and control. For a nonlinear PDE it generally is not.

Recall the semilinear example

$$
-\Delta y+c y^{3}=f+u.
$$

Its state linearization contains

$$
3cy^{2}\delta y.
$$

Therefore the adjoint operator depends on the current feasible state $y$. The callback
signature

```text
full point (y,u)
+
state objective covector D_y J_h
→
adjoint p
```

contains exactly the information needed to construct and solve the point-dependent
adjoint equation.

### 5.2 The adjoint right-hand side is a state covector

`solve_adjoint` does not receive the entire objective derivative. The full derivative is

$$
J_{h}'(y,u)
=
\left(
D_{y}J_{h},
D_{u}J_{h}
\right)
\in
Y_{h}^{\ast}\times U_{h}^{\ast}.
$$

Only the state component belongs on the right-hand side of

$$
D_{y}E_{h}(y,u)^{\ast}p
=
D_{y}J_{h}.
$$

The reduced formulation extracts that state covector before calling the adjoint solver.
This is one place where the state/control product decomposition from [Spaces,
coordinates, and duality](02-spaces-coordinates-and-duality.md) becomes operational rather
than descriptive.

## 6. The current state/control partition is deliberately narrow

The mathematical reduced argument is not limited to one state field and one control
field. The current `StateControlPartitionT`, however, implements a first narrow
contract.

It expects:

```text
variable layout
    exactly two blocks

test layout
    exactly one block
```

and records which variable block is the state and which is the control. From the full
variable layout it constructs one-block views:

```text
X_h = Y_h × U_h

full variable layout
├── state block
└── control block

selected views
├── state layout
└── control layout
```

The partition then provides three simple operations:

- compose a state and control into a full point;
- extract the state component of a full covector;
- extract the control component of a full covector.

This is enough for the current reduced DTO formulation.

It should not be read as a mathematical claim that reduced PDE optimization always has
exactly two fields. Mixed states, multiple controls, coupled equations, or other
partitions would require a richer formulation boundary.

### 6.1 Why the partition belongs between the model and formulation

`ExecutableModelT` only knows that its variable space is some product $`X_{h}`$.

It does not know that block 0 is to be eliminated as a state while block 1 is to be
optimized as a control.

`StateControlPartitionT` supplies that interpretation to the reduced formulation. So the
layers are:

```text
ExecutableModelT
    generic variable/test spaces and first-order actions
            │
            ▼
StateControlPartitionT
    identifies state and control factors
            │
            ▼
ReducedDTOT
    applies the reduced state–adjoint construction
```

That separation allows the same first-order model vocabulary to participate in
formulations with different block interpretations.

## 7. A reduced value evaluation starts with one state solve

Given a control $u$, the reduced objective value is

$$
j_{h}(u)=J_{h}(S_{h}(u),u).
$$

`ReducedDTOT::evaluate_value()` follows that definition directly. The sequence is:

```text
control u
    │
    ▼
solve_state(u)
    │
    ▼
state y
    │
    ▼
compose x = (y,u)
    │
    ▼
objective(x)
    │
    ▼
j_h(u)
```

The returned `ReducedValueEvaluationT` retains:

- the state $y$;
- the full point $(y,u)$;
- the control $u$;
- the objective value;
- the state-solve report.

It also carries private provenance/lifetime information that we will discuss later.

### 7.1 A value evaluation is already a substantial numerical result

The word *value* can make this object sound like a scalar cache. It is more useful to
think of it as a retained **feasible reduced evaluation**.

The scalar objective

$$
j_{h}(u)
$$

is meaningful only after the state associated with $u$ has been computed. That state is
often one of the most expensive products of the evaluation.

Retaining it means that if the derivative is requested afterward, the formulation can
reuse the same state and full point rather than solve the PDE again.

### 7.2 The state solve must report convergence

A state callback returns `FormulationSolveResultT`, which contains:

```text
solution
+
LinearSolveReport
```

For compiled iterative solves, the report can include:

- algorithm;
- preconditioner;
- iteration limit;
- iterations used;
- requested and achieved tolerances/residuals;
- convergence or failure.

`evaluate_value()` checks that the report says the state solve converged before using
the state to evaluate the objective. The reduced formulation therefore distinguishes

```text
a vector returned by a solver

from

a state accepted as the result of a converged solve
```

without taking ownership of the native solver policy.

## 8. Derivative augmentation reuses the retained state

Once a value evaluation exists, the reduced derivative can be added without repeating
the state solve. `ReducedDTOT::augment_derivative()` implements the adjoint derivation
in almost the same order as the mathematics.

Begin with the retained full point

$$
x=(y,u).
$$

Evaluate the full objective derivative:

$$
J_{h}'(x)
=
\left(
D_{y}J_{h},
D_{u}J_{h}
\right).
$$

Extract the state component

$$
D_{y}J_{h}
$$

and solve the adjoint equation

$$
D_{y}E_{h}(x)^{\ast}p
=
D_{y}J_{h}.
$$

Then evaluate the full residual VJP

$$
E_{h}'(x)^{\ast}p
=
\left(
D_{y}E_{h}(x)^{\ast}p,
D_{u}E_{h}(x)^{\ast}p
\right).
$$

Extract its control component and subtract it from the direct control objective
derivative:

$$
j_{h}'(u)
=
D_{u}J_{h}
-
D_{u}E_{h}(x)^{\ast}p.
$$

The runtime flow is therefore:

```text
retained (y,u)
      │
      ▼
J_h'(y,u)
      │
      ├── state component D_y J_h
      │           │
      │           ▼
      │      adjoint solve
      │           │
      │           ▼
      │           p
      │           │
      │           ▼
      │      E_h'(y,u)* p
      │           │
      │           └── control component
      │
      └── control component D_u J_h
                      │
                      ▼
              subtract pullback
                      │
                      ▼
                 j_h'(u)
```

The result is a `ReducedEvaluationT` containing:

- state;
- adjoint;
- full point;
- reduced derivative;
- objective value;
- state-solve report;
- adjoint-solve report.

### 8.1 The first component of the VJP is not discarded conceptually

The reduced derivative only needs the control component of

$$
E_{h}'(x)^{\ast}p.
$$

But the full VJP returns both state and control covectors. The state component should
satisfy

$$
D_{y}E_{h}(x)^{\ast}p
=
D_{y}J_{h}
$$

because $p$ was obtained from the adjoint solve. So the VJP remains a generic
full-variable transpose action. The reduced formulation uses the control block that
survives after state elimination.

This is an example of the distinction established in Part I:

```text
model capability
    full first-order action

formulation use
    selected block of that action
```

## 9. `evaluate()` is simply the complete two-stage evaluation

For callers that need the derivative immediately, `ReducedDTOT` also exposes

```cpp
evaluate(control)
```

which is conceptually

```text
augment_derivative(
    evaluate_value(control)
)
```

or mathematically:

```text
control
    ↓
state solve
    ↓
objective
    ↓
objective derivative
    ↓
adjoint solve
    ↓
residual VJP
    ↓
reduced derivative
```

The split methods and the combined method therefore describe the same formulation. The
difference is when the second stage is requested. That timing becomes important as soon
as an optimization algorithm evaluates trial points that may be rejected.

## 10. Why the two-stage evaluation matters to line searches

A line search may inspect several trial controls before accepting one. For each trial

$$
u_{\mathrm{trial}}
=
u+\alpha d,
$$

the reduced objective value requires a state solve:

$$
y_{\mathrm{trial}}
=
S_{h}(u_{\mathrm{trial}}).
$$

But not every acceptance policy needs the derivative at every rejected trial.

An Armijo backtracking test, for example, compares the trial objective against a
sufficient-decrease bound using the slope already known at the current iterate. A
rejected trial therefore needs:

```text
trial control
    ↓
state solve
    ↓
objective value
    ↓
reject
```

There is no reason to solve an adjoint for that rejected point. The retained value
interface makes that execution pattern possible:

```text
evaluate_value(trial)
    ↓
inspect objective
    ├── reject → discard retained value
    └── accept → augment_derivative(retained value)
```

The contract tests explicitly check that an Armijo search performs one state solve per
trial but only one adjoint solve, for the accepted trial.

### 10.1 Not every line search can defer derivative work

The split does not imply that every globalization policy can be value-only on rejected
trials.

A Wolfe condition also tests derivative information at the trial point. Such a policy
may need to augment every trial evaluation with an adjoint solve before it can decide
whether the step is acceptable.

So the important capability is not

> rejected trials never compute derivatives,

but rather

> value and derivative work are separable, so each optimization policy can request
> only the stage its acceptance test actually needs.

That distinction is visible in the tests: Armijo defers derivative work until
acceptance, whereas Wolfe searches may augment every trial.

## 11. A retained value is a coherent snapshot, not an arbitrary cache

If `augment_derivative()` accepted any object containing a state, control, and
objective, several subtle inconsistencies would become possible. For example, one could
accidentally combine:

- a state computed by one reduced service;
- a control with a compatible layout but different coefficients;
- a full point assembled in another session;
- an objective value associated with a different numerical realization.

The resulting objects might have the right dimensions and still not describe one
coherent reduced evaluation. `ReducedValueEvaluationT` therefore carries private
provenance information in addition to its public numerical data.

### 11.1 Evaluation-service provenance

Each `ReducedDTOT` constructs a private evaluation token. A value returned by
`evaluate_value()` retains that token. `augment_derivative()` checks that the token
matches the service performing the augmentation.

Thus a value produced by

```text
ReducedDTO A
```

cannot be passed to

```text
ReducedDTO B
```

merely because both use compatible layouts. The contract test checks this case
explicitly and verifies that the foreign value is rejected before an adjoint solve is
invoked.

### 11.2 The full point must still contain the retained state and control

The formulation also checks that the retained

```text
state
control
full_point
```

are mutually consistent. Conceptually,

$$
x
=
\mathrm{compose}(y,u)
$$

must still hold.

This matters because the derivative stage evaluates $`J_{h}'(x)`$ and $`E_{h}'(x)^{\ast}p`$ at
that full point. Reusing a state with a different control would break the mathematical
premise of the adjoint derivation.

### 11.3 Lifetime provenance is a separate issue

A compiled reduced service may capture or own backend objects that must remain alive for
later derivative augmentation. `ReducedValueEvaluationT` therefore also retains the
reduced service's lifetime owner.

The provenance check answers

> did this evaluation come from this reduced service?

The lifetime-owner check answers

> is it associated with the same native lifetime context?

Those are distinct from layout compatibility.

Chapter 13, **Integrating an existing PDE application**, examines these lifetime
relationships in detail. Here they matter because a retained value can outlive the
immediate function call that created it.

## 12. Solve reports make the state and adjoint part of the evaluation record

A reduced evaluation contains not only $y$, $p$, $`j_{h}(u)`$, and $`j_{h}'(u)`$, but also the
reports from the state and adjoint solves.

This is useful because the mathematical equations

$$
E_{h}(y,u)=0
$$

and

$$
D_{y}E_{h}(y,u)^{\ast}p=D_{y}J_{h}
$$

are only satisfied numerically to the accuracy of their solve policies.
`LinearSolveReport` records backend-neutral evidence such as:

```text
algorithm
preconditioner
iteration limit
iterations used
requested tolerance
achieved residual
termination status
```

The reduced formulation itself only requires convergence before continuing. Higher
layers can then retain the evidence for:

- work accounting;
- diagnostics;
- reproducibility;
- accepted-iteration records.

The solve policy remains with the producer that knows how the numerical operator is
realized.

### 12.1 Exact reference solves use the same result shape

Backend-neutral tests often have an exact or direct solve implemented by a tiny
reference model. `FormulationSolveResultT` still wraps that primal solution in the same
result type, with a synthetic report indicating a caller-supplied exact solve.

So the formulation sees the same conceptual object whether the producer used:

```text
exact dense algebra
CG
a direct sparse solve
an application-owned native solver
```

The distinction appears in the report rather than in the formulation API.

## 13. The external Step-4 path shows the composition directly

The Step-4 Problem B binding is a compact concrete example because the existing
application owns the state operator, control coupling, objective, state solve, adjoint
solve, and mass metric.

The binding assembles the reduced formulation from those independent native operations.
Schematically:

```text
Problem B native application
│
├── residual / JVP / VJP
│      ↓
│   CallbackExecutableModelT
│
├── solve_state
├── solve_adjoint
│      ↓
│   StateAdjointSolversT
│
├── state/control layouts
│      ↓
│   StateControlPartitionT
│
└──────────────────────────┐
                           ▼
                       ReducedDTOT
```

The constructor in `ProblemBBinding` makes that composition visible:

```text
model
partition
metric
reduced(model, partition, solvers)
```

The metric is beside the reduced formulation rather than part of the state/adjoint
derivative construction itself.

### 13.1 The state callback wraps the native solve result

Problem B's native state solver accepts the control coefficients and returns:

```text
native state vector
+
native solve evidence
```

The binding wraps the vector as a primal state block and converts the evidence into a
`LinearSolveReport`. Nothing about the reduced formulation needs to know that the
underlying native solve uses the adapted Step-4 matrix and CG implementation.

### 13.2 The adjoint callback does the same for the transpose system

The native adjoint solve receives the state-objective right-hand side. The binding wraps
its result as a primal value in the test layout, again carrying a solve report.

The reduced formulation then combines that adjoint with the VJP callback already exposed
through the executable model. This demonstrates a useful separation:

```text
native application
    owns how the adjoint equation is solved

executable model
    owns how E_h'(x)* p is evaluated

reduced formulation
    owns how those two capabilities are composed
```

## 14. The compiler-created path produces the same formulation ingredients differently

The semantic/compiler path reaches the same `StateAdjointSolversT` boundary from a
different direction. The compiler constructs a typed numerical model and then creates
solver callbacks that call methods such as

```text
solve_state_with_report(control, state_policy)

solve_adjoint_with_report(
    full_point,
    state_rhs,
    adjoint_policy)
```

Those callbacks are packaged into `StateAdjointSolversT`. So the producer paths differ:

```text
external application
    hand-written binding around native operations

semantic/compiler path
    registered lowerer constructs typed numerical services
```

but the reduced formulation consumes the same roles:

```text
ExecutableModelT
StateControlPartitionT
StateAdjointSolversT
```

This is the convergence point described in the architecture overview.

## 15. Reduced DTO means differentiation of the discrete problem

The class is named `ReducedDTOT`. Here DTO refers to **discretize then optimize**. The
formulation begins from an already-discrete executable problem:

$$
E_{h}(y,u)=0,
\qquad
J_{h}(y,u).
$$

Its objective derivative and residual VJP are derivatives of those discrete numerical
maps.

The alternative optimize-then-discretize viewpoint begins by deriving continuous
optimality equations before choosing a discretization. The project can also represent
supplied OTD systems, but through a different formulation product.

The distinction is important because the two approaches are related without being
identical implementations. `ReducedDTOT` is specifically the reduced state–adjoint
construction over the discrete model boundary developed in Part I.

### 15.1 The generic first-order model is broader than the reduced formulation

`ExecutableModelT` exposes:

```text
residual
residual JVP
residual VJP
objective
objective derivative
```

The first-order reduced path principally uses:

```text
objective
objective derivative
residual VJP
state solve
adjoint solve
```

during its main reduced evaluations. The complete executable-model port remains useful
because it describes a fuller first-order numerical model than this one formulation
happens to consume.

This also explains why an external application may have to provide residual/JVP
callbacks that the current reduced search does not call in its ordinary first-order
path.

That is an interface characteristic of the current implementation, not a new
mathematical requirement of the adjoint derivation.

## 16. The reduced derivative is still a covector

After derivative augmentation, the formulation returns

$$
j_{h}'(u)\in U_{h}^{\ast}.
$$

It does **not** have to choose a steepest-descent direction. That distinction was the
subject of [Metrics, gradients, and constraints](04-metrics-gradients-and-constraints.md).
Given a selected metric

$$
G:U_{h}\to U_{h}^{\ast},
$$

one may identify the derivative with the primal metric gradient

$$
g
=
G^{-1}j_{h}'(u).
$$

`ReducedDTOT` contains a small `gradient_direction()` helper that performs this inverse
metric application after checking layouts.

The current reduced search algorithms use their own direction-policy machinery to
perform the same primal/dual identification as part of constructing a search direction.

The mathematical boundary remains:

```text
ReducedDTOT
    produces j_h'(u) in U_h*

metric / direction policy
    turns that covector into a primal direction

globalization
    decides how far to move
```

Keeping these stages conceptually separate becomes especially useful for BFGS, Newton,
trust-region, and constrained methods.

## 17. Constraints do not alter the state–adjoint derivative formula

Suppose the control is restricted to an admissible set

$$
u\in C.
$$

The reduced objective and its derivative are still defined by

$$
j_{h}(u)=J_{h}(S_{h}(u),u)
$$

and

$$
j_{h}'(u)
=
D_{u}J_{h}
-
D_{u}E_{h}^{\ast}p.
$$

The constraint changes the **optimization problem over that reduced objective**. For
projected first-order methods, the solver combines the derivative with the metric and
projection machinery developed in Part I.

For active-set/KKT formulations, the same bound information participates in
complementarity conditions instead. So `ReducedDTOT` itself does not own a
`ConstraintT`. The constraint belongs to the algorithm/formulation layer that decides
how feasible controls are generated and how constrained stationarity is measured.

## 18. Second-order reduced information is an additional capability

The first-order reduced formulation gives us

$$
j_{h}'(u).
$$

Newton-type methods require the action of the reduced Hessian

$$
j_{h}''(u)[w]
\in
U_{h}^{\ast}
$$

for a primal control direction $`w\in U_{h}`$. That action is not implied merely by having:

- residual JVP;
- residual VJP;
- objective derivative;
- state and adjoint solves.

A genuine reduced-Hessian action may require:

- a tangent-state solve;
- second derivatives of the residual/objective;
- an incremental adjoint solve;
- problem-specific simplifications.

The project therefore represents it as a separate `ReducedHessianT` capability:

```text
(control u, direction w)
        ↓
reduced Hessian provider
        ↓
j_h''(u)[w] in U_h*
```

The provider owns whatever extra first- or second-order work is needed internally. This
keeps the first-order formulation honest: an implementation does not become a
Newton-capable model merely because it can compute first derivatives.

Chapter 6 examines how Newton and trust-region policies consume this
capability.

## 19. A complete reduced evaluation in equations and code

It is useful to place the mathematical construction and runtime methods side by side.

| Stage | Mathematics | Runtime operation |
| --- | --- | --- |
| Control | $`u\in U_{h}`$ | control `PrimalBlockT` |
| State | solve $`E_{h}(y,u)=0`$ | `solve_state(u)` |
| Full point | $x:=(y,u)$ | `partition.compose(y,u)` |
| Objective | $`j_{h}(u):=J_{h}(x)`$ | `model.objective(x)` |
| Full derivative | $`J_{h}'(x)`$ | `model.objective_derivative(x)` |
| Adjoint RHS | $`D_{y}J_{h}(x)`$ | `partition.state_component(...)` |
| Adjoint | solve $`D_{y}E_{h}(x)^{\ast}p=D_{y}J_{h}`$ | `solve_adjoint(x, state_rhs)` |
| Residual pullback | $`E_{h}'(x)^{\ast}p`$ | `model.residual_vjp(x,p)` |
| Reduced derivative | $`D_{u}J_{h}-D_{u}E_{h}^{\ast}p`$ | control-block subtraction |
| Metric gradient | $`G^{-1}j_{h}'(u)`$ | solver/metric layer |

The table is short because the preceding derivation has already supplied the meaning of
each row.

## 20. Following the reduced formulation through the source

The relevant source path is now compact.

### The formulation itself

Read:

- [`include/nmopt/contract/reduced_dto.hpp`](../../../include/nmopt/contract/reduced_dto.hpp)

A useful reading order inside the file is:

```text
StateControlPartitionT
StateAdjointSolversT
ReducedValueEvaluationT
ReducedEvaluationT
ReducedDTOT::evaluate_value
ReducedDTOT::augment_derivative
ReducedDTOT::evaluate
```

The implementation closely follows Sections 5–9 of this chapter.

### Solve evidence

Then read:

- [`include/nmopt/contract/linear_solve.hpp`](../../../include/nmopt/contract/linear_solve.hpp)

This file explains the common

```text
solution + LinearSolveReport
```

shape returned by state and adjoint solve callbacks.

### External application composition

Read:

- [`apps/external-dealii/step-4/minimal/problem_b_binding.hpp`](../../../apps/external-dealii/step-4/minimal/problem_b_binding.hpp)

The constructor of `ProblemBBinding` shows the model, partition, metric, and reduced
formulation assembled around application-owned numerical operations. The
`make_solvers()` helper shows how native state/adjoint solves become
`StateAdjointSolversT`.

### Compiler-created solve services

For the compiler path, the relevant construction currently lives in:

- [`include/nmopt/compiler/v1/dealii_compiler.hpp`](../../../include/nmopt/compiler/v1/dealii_compiler.hpp)

Search for `make_state_adjoint_solvers`. It wraps the typed compiled model's
`solve_state_with_report()` and `solve_adjoint_with_report()` methods into the same
formulation-facing callbacks. Part III explains how those typed
numerical models are built.

### Contract tests

Finally read:

- [`tests/contract/reduced_dto_contract.cc`](../../../tests/contract/reduced_dto_contract.cc)

This large test exercises far more than the narrow formulation, but several cases are
especially relevant here:

- value evaluation performs a state solve and objective evaluation without an
  adjoint solve or objective derivative;
- derivative augmentation reuses the retained value;
- foreign retained evaluations are rejected;
- line-search policies exercise the value/augmentation split;
- state and adjoint solve evidence is retained.

These tests make the intended runtime decomposition concrete.

## 21. What belongs to optimization algorithms

This chapter has answered:

> Given a control, how does nmopt construct the reduced objective and reduced
> derivative?

It has deliberately not answered:

> Once $`j_{h}(u)`$ and $`j_{h}'(u)`$ are available, how should the control be updated?

That next question introduces another layer of decisions:

- steepest descent versus nonlinear CG, BFGS, L-BFGS, or Newton;
- metric gradients and history;
- Armijo, Wolfe, exact quadratic, or fixed-step policies;
- trust regions;
- projected constrained steps;
- stopping criteria;
- Hessian inner solves;
- work counts and iteration evidence.

Those belong to **Reduced optimization methods**. The separation is useful because the
reduced state–adjoint formulation remains the same while the optimization strategy above
it changes.

## Read later

Useful existing documents for this chapter are:

- [Theoretical formalism](../../design/mathematical-model.md), for the project's
  Lagrangian, DTO, derivative, and adjoint conventions.
- [Reduced optimization](../overview/reduced-optimization.md), for the high-level
  runtime view.
- [External applications](../overview/external-applications.md), for the Step-4
  path as a peer producer of the same reduced formulation.
- [PDE, formulation, and solver boundary](../../design/pde-solver-boundary.md), for
  the accepted architectural separation among model actions, solve services, and
  optimization.

Chapter 7, **Optimality systems and KKT**, revisits the same first-order
conditions without eliminating the state, which makes the contrast between reduced and
all-at-once formulations explicit.
