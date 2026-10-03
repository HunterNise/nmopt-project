# Possible future extensions

> **Status: continuation ideas, not an active roadmap.**
>
> This document collects substantial capabilities that were left unimplemented,
> deferred, conditional, or incomplete when the project was closed. It does not
> commit the project to implementing them, assign priority, or imply that the old
> roadmaps can be resumed without re-evaluation.

The current implemented boundary is summarized in
[Implemented scope](../manual/overview/implemented-scope.md).

The items below are distilled from the retained
[implementation roadmap](implementation-roadmap.md),
[application roadmap](application-roadmap.md),
[Chapter 5 problem-library roadmap](chapter-5-problem-library-roadmap.md),
[Chapter 6 benchmark-suite roadmap](chapter-6-benchmark-suite-roadmap.md), and
the historical
[application/run/GUI draft](../history/application-and-gui-layer-draft.md).

Those documents record the state and assumptions of their development periods.
Any continuation should begin by checking the current source, references, and tests
rather than treating an old work unit as a ready implementation instruction.

## Application-level continuations

### Complete the Chapter 5 application recipe library

The compiler supports a broader family of scalar problems than the current
application catalog exposes as polished reusable recipes.

The retained application roadmap leaves completion of the Chapter 5 recipe layer as
an unimplemented continuation candidate.

A continuation could complete user-facing recipes for the registered:

- scalar distributed-control families;
- general scalar/Robin families;
- Neumann boundary-control families;
- specialized observation families; and
- registered Dirichlet-control families.

This is primarily an **application-authoring and usability extension**, not a request
for another numerical lowerer. Recipes should continue to produce semantic
specifications and reuse the existing compiler, formulation, and solver services.

Each advertised recipe should still have:

- typed parameters and defaults;
- explicit provenance;
- a fast manufactured configuration;
- focused derivative or KKT evidence;
- clear unsupported-combination diagnostics; and
- at least one comparison against direct semantic composition.

The existing Chapter 5 compiler registrations should be reassessed before deciding
which families are worth promoting into the public application catalog.

### Add the B3 and B4 constrained Chapter 6 studies

The retained benchmark programme selected B3 and B4 as the next application-scale
studies after B1/B2.

#### B3 — symmetric box-constrained Laplace control

The proposed first benchmark uses cellwise-discontinuous control, where the
implemented coefficientwise $L^{2}$ box semantics and PDAS path are already
well-defined.

A complete benchmark would add:

- a frozen source/replacement contract;
- a regularization and bound matrix;
- inactive-box agreement with the unconstrained solution;
- active-set stabilization;
- primal and dual feasibility;
- complementarity and stationarity evidence;
- active-set histories;
- native fields and report output; and
- comparison with projected reduced optimization where meaningful.

This would provide application-scale evidence for the KKT/PDAS machinery already
present in the framework.

#### B4 — asymmetric/spatially varying bounds

B4 would reuse the same KKT/PDAS path with spatially varying lower and upper bound
data.

The retained plan treats the two-dimensional case as sufficient for the first
framework benchmark; a three-dimensional reproduction was explicitly optional.

B4 should not require a new solver path. Its purpose would be to exercise richer
bound data and the corresponding study/evidence layer.

### Add the B5 and B6 all-at-once studies

The retained benchmark roadmap makes B5/B6 conditional later work rather than part
of the completed project slice.

#### B5 — all-at-once Laplace control

B5 would exercise the existing quadratic KKT formulation as an application-scale
all-at-once solve.

The proposed evidence includes:

- KKT action checks;
- agreement with the corresponding reduced solution;
- feasibility and stationarity residuals;
- multiplier-conversion evidence;
- solver diagnostics; and
- separate accounting of outer and inner work if preconditioning is introduced.

#### B6 — diffusion-reaction follow-up

B6 is intended as a regression/application reuse study rather than a new
formulation. It would replace the Laplace state action with diffusion–reaction while
keeping the all-at-once application structure.

These studies should only be activated after a new benchmark-selection decision and
a frozen source/replacement contract. Their historical presence does not make them
required completion work.

## Client and GUI layer

The historical application/GUI draft proposed a headless runner first and a Slint
desktop client second.

Most of the **headless foundation has since been implemented**:

- typed scenarios and catalogs;
- parameter files;
- the repository runner;
- run manifests;
- generated run directories;
- native mesh and field output;
- solver traces;
- artifact records;
- post-processing; and
- deterministic reports.

A future GUI should therefore consume the **current** application/execution
contracts. It should not recreate the earlier draft's runner or introduce another
solver path.

The genuinely unimplemented part is the desktop/client layer.

A possible GUI could provide:

1. target/scenario selection from the application catalog;
2. schema-backed parameter editing and validation;
3. review of the resolved run configuration;
4. launch of the headless runner as a child process;
5. live display of terminal or structured progress;
6. cancellation/status handling;
7. browsing of generated artifacts and reports; and
8. handoff of native field output to an external viewer such as ParaView.

The original proposal used Slint, but the technology choice should be reconsidered
if this work is ever reopened. The architectural constraint matters more than the
GUI toolkit:

> the GUI is a client of the headless application layer, not a second numerical
> application.

An embedded mesh viewer, remote/batch execution, polished packaging, resume
semantics, and broad multi-application discovery were explicitly later ideas in the
historical draft and should remain optional unless a concrete use case justifies
them.

## Broader PDE and discretization support

### Time-dependent optimal control

The implementation roadmap reserved a fixed-step temporal compiler but deliberately
did not activate it.

Its proposed first vertical slice was a backward-Euler heat equation represented as
one global residual over the trajectory, with an exact transpose of that discrete
residual.

Reopening this direction would require explicit policies for:

- trajectory storage or replay;
- temporal state reconstruction;
- differentiated controls;
- state and adjoint solve evidence;
- checkpointing if memory becomes relevant; and
- the relation between time discretization and optimization derivatives.

Adaptive time stepping and event handling should not be treated as part of the first
slice.

### Mixed and multi-equation PDE systems

The current reduced path is centered on one state block, one decision block, and one
test block. The retained roadmap defers general mixed/multi-equation support.

A concrete first extension was a steady Stokes vertical slice with:

- velocity and pressure fields;
- multiple equation/test blocks;
- a pressure gauge;
- a declared stable velocity/pressure finite-element pair;
- coupled state and adjoint solves;
- viscous, pressure-gradient, divergence, traction, and related operators; and
- exact blockwise transpose evidence.

A later boundary-control Stokes problem could then reuse those mixed-PDE pieces.

This work would be a genuine generalization of the formulation and compiler
boundaries, not another registered scalar target.

### Matrix-free and distributed backends

The retained implementation roadmap proposes two distinct stages:

1. establish serial matrix-free equivalence for a bounded realization; then
2. add distributed Trilinos/PETSc-style vector backends with explicit ownership and
   ghost policies.

The present project is serial and does not provide MPI scalability evidence.

Any distributed extension would need to preserve the existing distinctions between
primal/covector roles, model actions, solve services, metrics, and formulation
products rather than bypassing those contracts with backend-specific application
logic.

Scalability claims would require their own benchmark and measurement programme.

## State constraints and broader complementarity

### Regularized state constraints

Chapter 5 Section 5.12 was deliberately excluded from the implemented problem
library.

The retained roadmap proposes a first bounded extension based on a regularized
state-control observation

$$
w = y + \lambda u,
\qquad
\lambda \gt 0,
$$

with box constraints on $w$.

That path would need explicit:

- transformed decision/observation maps;
- regularization provenance;
- multipliers;
- complementarity;
- feasibility and stationarity diagnostics; and
- continuation in positive regularization parameters.

The original unregularized state-constraint problem has measure-valued multipliers
and remains a separate functional-analytic problem. It should not be represented
merely by reusing the current cellwise control-box machinery.

### Broader control-bound semantics

The current compiled PDAS path is intentionally tied to its registered
cellwise-discontinuous $L^{2}$ control representation.

Possible separate extensions include:

- continuous-control bound semantics;
- additional facewise/boundary complementarity representations; and
- constrained second-order or quasi-Newton reduced methods.

Each requires a declared relationship between primal coefficients, the dual
multiplier, the selected metric, and the mathematical bound. Raw coefficientwise
classification should not be generalized without such a policy.

## Broader formulation and solver support

### Automatic OTD derivation

The current framework executes **supplied** OTD blocks. It does not derive an
optimize-then-discretize adjoint or optimality system automatically from a strong PDE
description.

Automatic derivation would require a much stronger symbolic/weak-form contract than
the current supplied-OTD executor.

This should remain distinct from simply adding more supplied OTD registrations.

### Stabilization and GLS variants

The retained Chapter 6 scope explicitly excludes:

- GLS comparisons;
- stabilized Lagrangians; and
- general stabilized OTD/DTO variants.

These are possible future formulation families, but should be added only with
explicit discrete weak forms, spaces, pairings, derivative/transpose contracts, and
evidence showing which discrete formulation is being solved.

They should not appear as silent numerical switches on an existing formulation.

### Nonlinear KKT and SQP

The implemented KKT product is quadratic and the strongest reduced second-order path
assumes an explicit reduced-Hessian capability.

The long-term implementation roadmap leaves open:

- generic nonlinear KKT Newton methods;
- sequential quadratic programming;
- broader second-order constrained optimization; and
- nonlinear PDE optimal-control realizations that provide the required
  linearizations/Hessian information.

This would be a larger formulation/solver expansion and should not be inferred from
the existence of the current reduced Newton consumer.

### Block preconditioners

Reusable block preconditioners were intentionally left **conditional**, not required.

The retained plan identifies possible:

- block-diagonal Schur approximations;
- triangular/Bramble–Pasciak-style compositions;
- constraint preconditioners;
- approximate mass/PDE inverses;
- multigrid-backed block actions; and
- later Stokes/Uzawa variants if mixed PDE support is added.

A first preconditioner should be activated only by a concrete KKT application where
the existing serial solve path is inadequate.

Any robustness claim should be backed by mesh/parameter sweeps with separate inner
and outer iteration accounting.

## Compiler generality

The current deal.II compiler uses closed registered realization families.

One possible continuation is to broaden that registered surface with additional
carefully defined combinations. Another, much larger direction would be to introduce
more generally compositional lowering.

The latter should not be treated as an automatic goal. The current design
deliberately rejects unsupported combinations rather than guessing a nearby
realization.

Useful future compiler work should therefore begin from a concrete missing
mathematical family and define:

- the semantic policy;
- numerical spaces;
- residual/objective actions;
- transpose behavior;
- solve services;
- metric and constraint semantics;
- provenance; and
- focused verification.

Additional generality is valuable only where those contracts remain inspectable.

## Other long-term numerical directions

The retained implementation roadmap also mentions possible work on:

- nonmatching meshes;
- shape optimization; and
- additional specialized metric/observation policies.

Some related metric and observation capabilities originally listed alongside these
ideas were implemented later, so the historical list should not be copied verbatim
as an outstanding backlog.

Nonmatching meshes and shape optimization remain materially different problem
classes and would require their own explicit mathematical and numerical policies
before being added to the framework.

## Downstream library distribution

The current CMake project is used directly from the repository. It does not install
or export an `nmopt` package for downstream

```cmake
find_package(nmopt)
```

use.

If the project were ever intended to become a reusable external library rather than
a course-project repository, a separate distribution effort could add:

- installed public headers;
- exported CMake targets;
- package configuration/version files;
- dependency propagation;
- examples built as downstream consumers; and
- a stable statement of which headers/interfaces are public.

This is a packaging/productization extension, not a numerical capability, and is not
required for the current completed-project role.

## Extensions that should not be inferred from the old plans

The retained plans and historical drafts contain many ideas that later became
implemented. They should not be reopened simply because they still appear in old
documents.

In particular, the following already exist in the current project in some form:

- a headless application runner;
- typed scenarios and recipe boundaries;
- parameter-file configuration;
- run manifests and provenance;
- deterministic artifact records;
- native deal.II field output;
- solver histories and work evidence;
- post-processing and comparison plots;
- B1/B2 application execution;
- reduced Newton and trust-region support;
- supplied OTD execution;
- quadratic KKT products and serial solves;
- typed complementarity and PDAS;
- specialized point and normal-flux observations;
- $H^{-1}$ and $H^{1/2}$-related metric realizations; and
- the external-application integration boundary.

Any future work should start from the current implementation and
[Implemented scope](../manual/overview/implemented-scope.md), not from the
chronological state recorded by an older roadmap.
