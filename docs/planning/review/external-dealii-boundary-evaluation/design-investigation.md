# External deal.II boundary evaluation: design investigation

Status: review context; non-authoritative. Prepared on 2026-09-09 before E1.

This document preserves the architectural reasoning, external precedents,
source-supported facts, and open hypotheses that motivated the controlled
Step-4 evaluation. It defines no interfaces, implementation requirements, or
accepted architectural decisions. The
[evaluation roadmap](../../external-dealii-boundary-evaluation.md) owns the
E0–G1 protocol, sequence, freeze, gates, and phase status. This note does not
amend that protocol or add E1–E5 work.

The [PDE–solver boundary](../../../design/pde-solver-boundary.md) remains design
authority, and the [integration reference](../../../reference/external-dealii-solver-integration.md)
describes the implemented API. Durable conclusions require G1 evidence and an
accepted subsequent decision before promotion to design or public reference
documentation. An early G1 obstruction report is possible; it does not imply
completion of the numerical comparison.

Read this note when interpreting G1 or considering later architecture changes.
It is not required implementation reading for E1–E5.

## 1. Why the evaluation exists

The completed PDE–solver refactor established semantic/compiler and external
applications as separate producers of common downstream operations. Its
[external Poisson fixture](../../../../tests/application/external_application_dealii_contract.cc)
exercises that API, including derivative checks, optimization, objective
replacement, and native output. This is concrete functional evidence. A
fixture written for the contract provides limited evidence about the cost of
adapting a pre-existing application with different access and lifecycle rules.

The next attempt used authentic upstream Step-4. Its effort combined at least
six causes: opening a closed forward application for reuse; defining control,
objective, and adjoint mathematics; adding verification; constructing nmopt
layouts and callbacks; satisfying possibly unnecessary formulation
capabilities; and any deeper coupling between layers. The total wrapper size
could not distinguish those causes. Neither could a native wrapper that had
already been shaped around the expected nmopt construction surface.

An earlier review asserted that the refactor had "stopped one layer too early."
That assertion was withdrawn as an architectural conclusion in the subsequent
review synthesis: unused capabilities were visible, but material integration
cost had not been established. It was not experimentally refuted; the
interpretation was returned to open-hypothesis status.

The resulting experiment preserves the forward application, defines the added
problem separately, verifies an independent native reduced path, and compares
the same mathematics and selected numerical policies through current nmopt.
Native reduced and optimizer orchestration are services to account for, not
work that can be omitted when comparing the binding's size with the native
path. Problem A removes new FE coupling assembly from the first comparison;
it does not remove the need for a coherent adjoint and objective.

### Evidence provenance and vocabulary

The source snapshots are the completed refactor
`67104fd376a9b1d25d7e6cdaca664a7fe1aa98ed` and the tutorial attempt
`b4dce25ab03d893df48ab1d9eee2e4dc7cf9183d`. The protocol was committed as
`e9ae1ac`; these documentation changes do not alter the reviewed numerical
sources. Future reviewers should use the recorded revisions when files move
or the API changes.

Historical reasoning was reconstructed from the local
`.codex/drafts/Project_Architecture_Overview.md` transcript, whose messages are
separated by `> ---`, and the subsequent review synthesis beginning "Thanks.
Your two reviews changed how I think we should proceed." The latter explicitly
rejects drawing the early architectural conclusion. These are review-history
sources, not implementation authority. This note preserves the relevant
reasoning without requiring future readers to recover those local materials.

- **Source-supported fact:** a statement checked against a named code snapshot.
- **External precedent:** documented behavior of another library, with a primary
  source. Its relevance to nmopt is an interpretation, identified as such.
- **Historical reasoning:** why a question or experimental control was chosen.
- **Open hypothesis:** a possible explanation awaiting discriminating evidence.
- **Observed result:** evidence from an actual run, as defined by the roadmap.
  This note introduces no E1–E5 run results. Expected counts remain derived
  schedule invariants, and acceptance thresholds remain frozen protocol.

## 2. External architectural precedents

These sources were checked for this note on 2026-09-09. Versioned links identify
the examples examined; rolling documentation is background, not a dependency
pin for the experiment. The comparisons concern responsibilities, not API size
rankings or measured integration effort.

### SUNDIALS: mathematics, orchestration, algebra, and solves

KINSOL receives the nonlinear system function and can pass application data
through `KINSetUserData()`. Its selected nonlinear and linear solver choices
determine additional requirements. For supported cases, Jacobian or
Jacobian-vector approximations can be supplied by internal difference
quotients rather than mandatory application implementations. This is a
specific example of optional or derived operations, not a claim that every
SUNDIALS capability has a default. [KINSOL usage, v7.7.0](https://sundials.readthedocs.io/en/v7.7.0/kinsol/Usage/index.html).

`N_Vector` separates numerical operations from implementation-specific
contents, constructors, and destruction. Custom implementations provide a
required operation set, with additional optional operations. This separates
algebra/storage from application equations; it does not eliminate vector
compatibility and lifecycle obligations. [NVECTOR API, v7.7.0](https://sundials.readthedocs.io/en/v7.7.0/nvectors/NVector_API_link.html).

The `SUNLinearSolver` API distinguishes matrix-based, matrix-free, and
matrix-embedded services. A matrix-embedded service can keep matrix formation
inside its solve implementation while the external `SUNMatrix` argument is
null. Support and compatibility depend on the consuming package. The useful
precedent is that a solve need not expose its matrix across every boundary;
this does not make a native solve automatically compatible with every
SUNDIALS integration. [SUNLinearSolver API, v7.7.0](https://sundials.readthedocs.io/en/v7.7.0/sunlinsol/SUNLinSol_API_link.html).

CVODES supplies backward integration infrastructure, forward-state recovery,
checkpointing, and interpolation. The application still defines the backward
problem, its endpoint data, and its RHS, including adjoint mathematics when
that is the chosen backward problem. Library-owned adjoint orchestration does
not mean automatic derivation of an application's adjoint equations.
[CVODES adjoint usage](https://sundials.readthedocs.io/en/latest/cvodes/Usage/ADJ.html).

**Interpretation for nmopt:** capability selection and ownership can be
separated. This motivates measuring which operations the chosen formulation
uses and which construction requires. It does not establish that nmopt should
copy a C callback table, introduce a new solve hierarchy, or derive missing
derivatives numerically in this experiment.

### deal.II KINSOL: short wiring can coexist with application changes

Step-77 connects vector initialization, residual evaluation, Jacobian setup,
and linear solution to KINSOL. Its discussion also explains why separating
residual and Jacobian assembly from Step-15 requires care with nonzero
Dirichlet data: the initial state carries prescribed values, while updates
use homogeneous boundary data. Thus a small callback connection can sit on
top of mathematically significant application restructuring.
[deal.II Step-77, v9.5.0](https://dealii.org/9.5.0/doxygen/deal.II/step_77.html).

**Interpretation for Step-4:** count the restructuring independently of the
callback connection. Step-77 supplies no numerical estimate of Step-4's
adaptation cost and does not justify importing its boundary-coordinate policy
into frozen Problem A.

### PETSc/TAO: a different division of services

TAO accepts objective and gradient callbacks using PETSc vectors, with further
capabilities required by the selected method. Its solver environment uses
PETSc algebra and linear-solver facilities. An application can put an already
reduced objective inside those callbacks; the reduced state/adjoint
orchestration then belongs behind them rather than being supplied by that
generic optimization callback interface. [TAO manual](https://petsc.org/release/manual/tao/).

**Interpretation:** comparing a short TAO connection to a complete nmopt
reduced formulation without accounting for the application-owned reduction
would omit a service. Deeper algebra integration is another architectural
choice, not evidence that nmopt's existing backend contract is wrong.

### ROL and Thyra: mathematical semantics and capability negotiation

ROL's `Objective_SimOpt` distinguishes simulation and optimization arguments
and their derivative blocks. Its abstract `Vector` includes `dual()` and
`apply()` as well as an inner product, making primal/dual interpretation
explicit. These are precedents for retaining state/control and geometry
semantics when simplifying construction.
[ROL Objective_SimOpt, Trilinos 13.2.0](https://raw.githubusercontent.com/trilinos/Trilinos/trilinos-release-13-2-0/packages/rol/src/function/simopt/ROL_Objective_SimOpt.hpp),
[ROL Vector, Trilinos 13.2.0](https://raw.githubusercontent.com/trilinos/Trilinos/trilinos-release-13-2-0/packages/rol/src/vector/ROL_Vector.hpp).

Thyra's `ModelEvaluatorBase` exposes input/output support queries, derivative
representations, and operator/solve-related outputs. This is a richer form of
capability negotiation. **Interpretation:** such expressiveness also creates
more concepts and compatibility cases for producers and consumers to manage.
That complexity observation is not a measured defect or a reason to reject
Thyra in its own setting.
[Thyra ModelEvaluatorBase, Trilinos 13.2.0](https://raw.githubusercontent.com/trilinos/Trilinos/trilinos-release-13-2-0/packages/thyra/core/src/interfaces/nonlinear/model_evaluator/fundamental/Thyra_ModelEvaluatorBase_decl.hpp).

The precedents support investigating the division of responsibilities. They
do not establish a universal rule that every mature library exposes exactly
the smallest possible capability set to every consumer.

## 3. Source-supported nmopt facts

The following facts describe the reviewed code; their material consequences
remain to be measured.

1. [`CallbackExecutableModelT`](../../../../include/nmopt/contract/callback_executable_model.hpp)
   requires nonempty residual, residual JVP, residual VJP, objective, and
   objective-derivative callbacks at construction. Providing only the callbacks
   used by one runtime path is not valid construction.
2. [`ReducedDTOT`](../../../../include/nmopt/contract/reduced_dto.hpp), in
   `evaluate_value()` and `augment_derivative()`, calls state solve, objective,
   objective derivative, adjoint solve, and full residual VJP. It does not call
   residual or residual JVP. It uses only the control component of the full
   residual pullback when forming the reduced derivative. The full callback
   contract still requires a valid state component for arbitrary seeds.
3. [`MetricT`](../../../../include/nmopt/contract/metric_constraint.hpp) requires
   primal-to-covector action and inverse action, along with identity/layout
   information. It does not require an assembled mass matrix. The existing
   `MassMetric` is one realization. Objective regularization and metric
   geometry have distinct roles in the
   [accepted design](../../../design/pde-solver-boundary.md#regularization-metric-and-reduced-hessian).
4. `StateAdjointSolversT` accepts application callbacks returning solution
   blocks and reports. Its signature does not transfer matrix or solver-object
   ownership. `ReducedDTOT` checks convergence and layout compatibility; it
   does not choose CG or UMFPACK. The
   [solve-result convenience constructor](../../../../include/nmopt/contract/linear_solve.hpp)
   declares an exact successful solve without inspecting native convergence.
5. The current [`ReducedSearchSolverT`](../../../../include/nmopt/solvers/reduced_gradient.hpp)
   takes a concrete `ReducedDTOT<Backend>` reference and a metric; its policy
   parameters customize direction and line search. It uses the DTO's staged
   value/derivative evaluations, including service-identity checks inside the
   DTO. It is not constructed directly from an arbitrary reduced-objective
   callable. This is a formulation dependency, not proof of inappropriate
   compiler or PDE coupling.
6. That optimizer supplies trial construction, stopping checks, state reuse
   across accepted trial augmentation, histories, and work reporting. These
   services belong in the native-versus-nmopt comparison. Policy templates and
   [optional capabilities](../../../../include/nmopt/solvers/reduced_search.hpp)
   already provide some separation; the architecture is not uniformly maximal.
7. The [external fixture](../../../../tests/application/external_application_dealii_contract.cc)
   explicitly checks objective replacement without changing PDE residual
   callbacks. Its [application](../../../../tests/dealii/external_poisson_fixture.hpp)
   already exposes the required numerical operations. **Interpretation:** it
   supports functional composition more directly than it supports claims
   about adapting a closed upstream application. No tests were rerun for this
   documentation note.

Layout compatibility and solve-report details are recorded in the integration
reference; repeating their full construction examples here would create a
second API guide. None of these facts alone establishes excessive cost.

## 4. What the previous Step-4 attempt could not isolate

The [pinned upstream source](../../../../apps/external-dealii/step-4/source/upstream/step-4.cc)
has public construction and `run()`, private numerical state, CG with identity
preconditioning, and VTK output under the field name `solution`. Making its
assembled system independently callable is application-reuse work, even for a
consumer with no nmopt dependency.

The tutorial-attempt files below refer specifically to `b4dce25`. They are
scheduled for replacement or removal by E2, so their Git revision and symbols
are the historical locators. Paths are under `apps/external-dealii/step-4/`.

| Source at `b4dce25` | Relevant evidence |
| --- | --- |
| `step-4.cc`, `Step4` | Added preparation and mesh/FE/DoF/matrix/RHS views; original solve/output remain private |
| `tutorial_application.hpp`, `TutorialApplication` | Native residual/JVP/full-VJP, objective/partials, state/adjoint solves, metric accessor, and output surface |
| `tutorial_application.cc`, `Impl::assemble_metric()` | New FE mass assembly used by the objective and exposed as a control metric |
| `tutorial_application.cc`, `solve_state()` and `solve_adjoint()` | New `SparseDirectUMFPACK` construction and initialization on each solve |
| `tutorial_application.cc`, `write_native_output()` | Separate `DataOut` path writing VTU with field `state` |
| `tutorial_binding_smoke.cc`, `main()` | Native zero-control solve, residual/shape/output checks; no nmopt DTO or optimizer |

Retrieve a retired source, for example, with:

```bash
git show b4dce25:apps/external-dealii/step-4/tutorial_application.cc
```

The wrapper includes no nmopt library headers, but its native surface already
mirrors most executable-model obligations. Independence of includes therefore
does not by itself establish an independently chosen minimal native reference.
This is a comparison-design concern, not proof that the author intended to
inflate the wrapper.

On the wrapper path, UMFPACK takes the place of Step-4's existing CG machinery,
and native output is implemented a second time with different format, filename,
and field identity. These changes alter numerical and output policies before
nmopt has entered the comparison. The standalone Step-4 path still retains its
original solve/output; saying the entire application was converted to UMFPACK
would overstate the change.

The wrapper's mathematics combines identity load-vector control with a
mass-weighted objective. With $A,b$ denoting the boundary-eliminated Step-4
system and $M$ the newly assembled mass matrix, its definitions are:

```math
E(y,u)=Ay-b-u,
\qquad
J(y,u)=\frac{1}{2}y^{\mathsf T}My+
       \frac{0.2}{2}u^{\mathsf T}Mu.
```

Mass assembly follows that objective/metric choice. The residual does not
apply a mass or FE coupling matrix to the control. Adding $M$ to the objective
does not make identity forcing into ordinary FE distributed control. Since
the control acts on the already-eliminated full system, it also perturbs
boundary rows. This can define a coherent algebraic problem, but it must be
named and interpreted accordingly.

The frozen Problem A deliberately keeps that full-vector algebraic meaning
while selecting a simpler objective and identity metric, and preserving CG and
original output machinery. Those choices define a new controlled probe; its
results must not be presented as solving the old mass-weighted problem. The
roadmap alone owns its equations and constants. Problem B later introduces
genuine FE control/coupling and coordinate treatment as application mathematics.

Earlier draft/hybrid attempts help explain the concern, but their provisional
line counts and suggested replacement interfaces are not a reproducible cost
baseline. The causal ledger instead records reuse seams, mathematical
operations, verification, native orchestration, construction obligations, and
runtime work with multiple consumer flags. An operation introduced for
verification can still be mandatory for nmopt construction.

## 5. Open hypotheses and discriminating evidence

These are review questions, not additional acceptance gates or implementation
instructions. Existing requirements for composition, explicit lifetimes,
application-owned solves/output, and derivative/metric separation are not
themselves open hypotheses. H4–H6 ask whether the current implementation
delivers them effectively. Support or falsification is bounded to the tested
case; Problem A cannot settle general claims about nonlinear, nonsymmetric,
constrained, or distributed applications.

### H1 — Complete executable construction causes material extra obligations

Status: **open before E1**.

Plausible because residual/JVP construction and a complete VJP exceed the
selected runtime use. Support would be a causal inventory showing unavoidable
extra mathematics, data exposure, maintenance, or repeated operator work under
the valid frozen API. A smaller binding alone would not remove such an
obligation. The material-cost claim weakens if valid operations are already
available or cheap and the resulting binding is modest. The mere existence
of an unused callback establishes breadth, not its importance.

### H2 — Most integration friction is mechanical construction

Status: **open before E1**.

Plausible because layouts, packing, captures, and report translation repeat
existing information. Support would be repeated explicit construction with
clear ownership and no missing mathematical capability. A helper could then
be considered after the freeze. Evidence against helper sufficiency would be
mandatory new mathematics, discarded runtime work, or a correctness
obstruction that wrapping the same API cannot remove. H1 and H2 may both hold
for different ledger rows; neither needs to explain every cost.

### H3 — Verification and runtime capabilities merit different construction status

Status: **open before E1**.

Plausible because residual/JVP checks have consumers outside the chosen
optimization schedule. Support would distinguish the independently needed
verification operations from obligations imposed on every runtime binding,
and identify a useful consumer lacking some of them. The hypothesis weakens
if separation adds lifecycle/validation complexity with negligible benefit.
The experiment keeps valid callbacks and shared mathematical definitions;
separate counter scopes do not make construction obligations disappear.

### H4 — Current solve contracts preserve native policy and ownership adequately

Status: **open before E1**.

Opaque native solve callbacks are already permitted by design and signature.
Support would be faithful supplied-RHS CG reuse, actual success reports, native
exception propagation, and no framework-required matrix ownership transfer.
Counter-evidence would be a minimal reproducer showing a valid native solve
cannot be represented without changing policy, ownership, or report meaning.
The old wrapper's voluntary UMFPACK replacement is not such evidence. Step-4
cannot establish adequacy for every nonlinear or distributed solve lifecycle.

### H5 — The existing producer/formulation convergence point is sufficient

Status: **open before E1**.

Separate compiler and external producers already exist; the question concerns
the adequacy of their common operations. An authentic binding needing only
documented public contracts supports the current placement. Required compiler
products, private realization access, or reconstruction from compiler
provenance would argue against it. Concrete dependence on `ReducedDTOT`
justifies examining the formulation boundary if it obstructs an intended
consumer, but not immediately replacing it. Problem A provides external-path
evidence only; compiler-versus-external equivalence remains deferred.

### H6 — Existing composition seams preserve independent mathematical choices

Status: **open before E1**.

Objective/PDE separation and derivative-versus-metric semantics are accepted
requirements. The fixture's objective-replacement check makes successful
composition plausible. Evidence would be that relevant independent choices
remain localized to their mathematical owners. Counter-evidence would be a
necessary change to unrelated PDE or solver code when only an objective or
metric changes. Problem A's fixed quadratic objective and identity metric
provide limited coverage; this hypothesis does not authorize extra variants
or let identity geometry erase primal/dual distinctions from the review.

### H7 — Friction reflects deeper structural coupling

Status: **open before E1**.

This is a competing explanation, not a deduction from wrapper size. Support
would require traced dependencies showing independent changes propagate
across supposedly unrelated layers, ideally in more than one consumer, and
cannot be isolated by a local boundary correction. Examples to investigate
if observed include an objective change forcing compiler/backend/PDE changes,
or output requiring recovery of hidden compiler types. The hypothesis weakens
if integration remains local, numerical services compose correctly, and costs
are attributable to reuse, mathematics, or explicit construction. An
experiment bug, failed tolerance, or one awkward adapter is insufficient.

## 6. Conclusions not established

The available evidence does not justify concluding that:

- the refactor necessarily stopped one layer too early;
- `ReducedDTOT` should be replaced by `ReducedFirstOrderPorts` or any other
  particular new interface;
- a builder/helper necessarily removes the significant costs;
- SUNDIALS's API structure should be copied;
- Problem A generalizes to realistic FE optimal-control problems;
- source line count alone demonstrates architectural failure;
- a mass matrix is intrinsically required by external integration;
- callback or solve-count differences are timing or allocation measurements;
- the authentic integration is successful merely because the fixture passes,
  or failed merely because unused construction requirements exist.

A reproducible frozen-API obstruction can justify an early G1 finding before
E4/E5 complete. Such a finding must state which comparison could not run and
why. Conversely, a successful Problem A comparison would establish a bounded
success, not universal ergonomic sufficiency.

## 7. Use and promotion at G1

Interpret the roadmap's evidence and causal ledger alongside these hypotheses.
For each stable ID, record **supported for the tested case**,
**weakened/refuted for the tested case**, or **unresolved**, with evidence
locations, counter-evidence, and limits. An untested generalization stays
unresolved. An early obstruction report leaves unperformed checks unperformed.

Experiment-specific facts and attribution belong in the G1 report at
`docs/planning/review/external-dealii-boundary-evaluation/g1-report.md`. Keep raw traces
and working estimates in the roadmap's ignored evidence directories. This
reviewed reasoning record is distinct from those raw artifacts and does not
make speculative explanations authoritative.

Only an accepted durable conclusion should lead to a subsequent design or
public-reference change. Helper proposals, boundary changes, deeper refactors,
and Problem B each need their own accepted scope. Preserve this note as the
historical question set, recording later resolutions rather than rewriting
the initial hypotheses as if the outcome had been known before E1.
