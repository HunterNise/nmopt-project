# V0 executable contract

## Scope

This is the first code-level realization of the default policies in the
[implementation-readiness review](../implementation-readiness-review.md). It is
backend-neutral on purpose: the [v1 deal.II compiler](../v1/semantic-compiler.md)
lowers its supported semantic graph to this contract, while the contract
itself does not contain a DoFHandler, finite element, matrix type, or solver
selection.

The current reference implementation is dense and algebraic. It validates the
contract and DTO signs; it is not presented as an FE discretisation or a
replacement for deal.II.

The core types are parameterised by a backend. A backend supplies a copyable
vector type and five operations: construct zeros of a requested size, report
size, dot, scaled addition, and scaling. The dense backend is the tested
reference instance. A deal.II backend will provide the same operations for its
owned and ghosted vector policy; it will not alter the model, metric, or DTO
signatures.

## Exact v0 types

The public executable and solver headers provide:

| Type | Meaning |
|---|---|
| `SpaceId` and `BlockLayout` | Stable semantic space identifiers and compatible block dimensions. |
| `PrimalBlock` | A typed block of coefficients in a declared discrete primal space, with read-only vector access and checked algebraic updates. |
| `CovectorBlock` | A typed block of coefficients in the dual of that layout, with the same layout-preserving storage boundary. |
| `pair` | The declared dual-coefficient pairing. It is the only primitive that pairs a covector with a primal vector. |
| `ExecutableModel` | Residual value, residual JVP, residual VJP, objective value, and objective derivative. |
| `Metric` | Explicit primal-to-dual action and inverse action, plus an opaque realization witness distinct from its display ID. The serial deal.II backend supplies `dealii_backend::MassMetric` for one-block sparse mass matrices. |
| `Constraint` | Feasibility and projection coupled to an actual compatible metric realization. |
| `LinearSolveReport` and `FormulationSolveResultT` | Backend-neutral state/adjoint convergence, tolerance, and work evidence paired with a solved primal block. |
| `ReducedDTO` | The state–adjoint–reduced-covector workflow for one state and one binary decision block. |
| `contract::ReducedHessianT` | An explicit reduced-Hessian capability that applies $`H(u)w`$ as a typed reduced covector; first-order DTO ports do not imply this capability. |
| `solvers::ReducedSearchDirectionT` | A typed primal search direction with its metric gradient norm, reduced-covector directional derivative, and action counts. |
| `solvers::ReducedLineSearchResultT` | A typed line-search result carrying the accepted trial DTO, actual step length, trial count, and any Hessian-action count. |
| `solvers::NonlinearConjugateGradientDirectionPolicyT` | The default metric-aware Polak–Ribière+ direction policy with typed primal/dual history and deterministic restarts. |
| `solvers::FletcherReevesDirectionPolicyT` | The metric-aware Fletcher–Reeves direction policy with the same typed history and deterministic restart protocol. |
| `solvers::QuadraticConjugateGradientDirectionPolicyT` | The strict classical quadratic-CG recurrence with metric-gradient history, periodic dimension restarts, and contract-checked descent. |
| `solvers::LimitedMemoryBfgsDirectionPolicyT` | The metric-aware limited-memory BFGS direction policy with bounded typed secant history and explicit curvature resets. |
| `solvers::NewtonDirectionPolicyT` | The explicit-Hessian Newton direction consumer using metric-preconditioned inner conjugate gradients. |
| `solvers::ArmijoLineSearchPolicyT` | Backtracking Armijo acceptance using the declared pairing and the actual returned trial displacement. |
| `solvers::ExactQuadraticLineSearchPolicyT` | One-step exact line search for a positive-curvature explicit reduced Hessian. |
| `solvers::WolfeLineSearchPolicyT` | Strong Wolfe acceptance using actual trial slopes and declared sufficient-decrease/curvature fractions. |
| `solvers::ReducedSolverResultT` | The shared reduced-solver report containing accepted objectives, gradient norms, accepted steps, objective changes, stopping state, formulation/metric/Hessian action counts, and explicit direction-reset counts. |
| `solvers::ReducedTrustRegionResultT` | The trust-region report containing accepted objectives, gradient histories, radius/step histories, actual and predicted reductions, ratios, acceptance flags, and action counts. |
| `solvers::ReducedGradientSolverT` | Backend-parametric unconstrained or projected reduced Armijo method consuming `ReducedDTOT`, `MetricT`, and an optional `ConstraintT`. |
| `solvers::ReducedConjugateGradientSolverT` | The same reduced execution loop configured with `NonlinearConjugateGradientDirectionPolicyT`. |
| `solvers::ReducedFletcherReevesSolverT` | The same reduced execution loop configured with `FletcherReevesDirectionPolicyT`. |
| `solvers::ReducedExactConjugateGradientSolverT` | The reduced PR+ execution loop combined with `ExactQuadraticLineSearchPolicyT`. |
| `solvers::ReducedExactFletcherReevesSolverT` | The reduced Fletcher–Reeves execution loop combined with `ExactQuadraticLineSearchPolicyT`. |
| `solvers::ReducedQuadraticConjugateGradientSolverT` | The strict classical quadratic-CG recurrence combined with `ExactQuadraticLineSearchPolicyT`. |
| `solvers::ReducedLimitedMemoryBfgsSolverT` | The same reduced execution loop configured with `LimitedMemoryBfgsDirectionPolicyT` for the unconstrained first registration. |
| `solvers::ReducedNewtonSolverT` | The same unconstrained reduced execution loop configured with `NewtonDirectionPolicyT` and an explicit `ReducedHessianT`. |
| `solvers::ReducedExactNewtonSolverT` | The reduced Newton loop combined with `ExactQuadraticLineSearchPolicyT` for positive-curvature quadratic targets. |
| `solvers::ReducedWolfeGradientSolverT` | The reduced steepest-descent loop combined with `WolfeLineSearchPolicyT`. |
| `solvers::ReducedTrustRegionSolverT` | The unconstrained matrix-free Cauchy trust-region solver consuming `ReducedDTOT`, `MetricT`, and an explicit `ReducedHessianT`. |

The unsuffixed public aliases select the dense reference backend. The
corresponding types with a T suffix are backend-parametric, for example
`ExecutableModelT`, `PrimalBlockT`, `CovectorBlockT`, `MetricT`, and `ReducedDTOT`.
Block vectors are supplied at construction and exposed read-only afterward.
The retained `add_scaled_block()` and `scale_block()` operations cannot replace
storage, and scaled updates reject a source whose dimension differs from the
declared block layout.

`dealii_backend::SerialBackend` is the checked adapter between contract
dimensions (`std::size_t`) and the native serial deal.II vector size type. It
accepts the native maximum symbolically without allocation and rejects a
larger contract dimension before vector construction. The cellwise and
facewise box adapters use the same conversion and native index type.

An `ExecutableModel` has the following exact signatures in mathematical form:

$$
  r_{h}(x_{h})\in Z_{h}^{\ast},\qquad
  r_{h}'(x_{h})\delta x_{h}\in Z_{h}^{\ast},\qquad
  r_{h}'(x_{h})^{\ast}p_{h}\in X_{h}^{\ast},\quad p_{h}\in Z_{h},
$$

$$
  J_{h}(x_{h})\in\mathbb R,\qquad J_{h}'(x_{h})\in X_{h}^{\ast}.
$$

The DTO builder uses the global convention

$$
  \mathcal L_{h}(x_{h},p_{h})=J_{h}(x_{h})-\langle p_{h},r_{h}(x_{h})\rangle.
$$

For $`x_{h}=(y_{h},u_{h})`$, it therefore computes

$$
 E_{y}'(x_{h})^{\ast}p_{h}=J_{y}'(x_{h}),\qquad
 j_{h}'(u_{h})=J_{u}'(x_{h})-E_{u}'(x_{h})^{\ast}p_{h}.
$$

No method accepts a PDE name, boundary-condition enum, or untyped vector.

## Chosen representation

V0 uses dual coefficients. If

$$
  r_{j}=\langle E_{h},\psi_{j}\rangle,
$$

then pairing a residual with a test vector of coefficients $p$ is
$r^{\mathsf T}p$. Consequently an assembled Jacobian with test rows and
trial columns has coordinate VJP $J^{\mathsf T}p$. A future backend using
Riesz representatives must provide a different pairing implementation; it
cannot reuse this rule implicitly.

## V0 reduced formulation boundary

The contract deliberately supports only:

~~~text
two variable blocks: one eliminated state and one decision/control-or-parameter
one residual test block
an externally supplied state solve
an externally supplied adjoint solve
DTO only
first derivatives only
~~~

The state/adjoint solves are formulation services, not residual-term methods.
They must operate on the compiled model and honor its constraint, nullspace,
and tolerance policy. Each service returns a `FormulationSolveResultT` whose
typed report states the algorithm, preconditioner, requested tolerances,
iteration limit, work, achieved residual, and termination. `ReducedDTOT`
rejects a nonconverged result before using its primal block. The constructor
validates all block-layout connections; no implicit state/control selection is
made.

The direct `ReducedDTOT` constructor borrows its executable for deliberately
short-lived reference use. Its compiled constructor owns the executable and
an optional backend-session lifetime token, so a reduced service detached from
its `CompiledProblemT` remains valid.

Mixed partitions, multiple equation blocks, nonlinear all-at-once Newton, OTD,
and generic nonlinear second-order actions are intentionally outside this v0
contract. The explicit `ReducedHessianT` capability is the narrow exception
for selected models that provide their own exact or declared Hessian action.
The historical `StateControlPartitionT` name denotes this binary decision
port; a v1 coefficient lowerer may give its layout the semantic identifier
`parameter` without changing the backend-neutral value/JVP/VJP core.

## First reduced solver

`solvers::ReducedGradientSolverT` is the first generic optimizer.
`ReducedSearchDirectionT` stores the primal direction $`d`$, its metric
gradient norm, $`j_{h}'[d]`$, and the metric/Hessian action counts used to form
it. The initial steepest-descent policy forms
$`g=G^{-1}j_{h}'`$ and supplies $`d=-g`$ to the solver. An unconstrained
trial is therefore $`u^{+}=u+\alpha d`$ and is accepted only when it satisfies

$$
  j_{h}(u+\alpha d) \leq j_{h}(u)+c\langle j_{h}',\alpha d\rangle.
$$

The constraint-qualified constructor accepts a `ConstraintT` only when its
layout matches the metric and it declares projection support for that exact
metric realization. Display identifiers are retained for provenance and do
not grant projection capability.
It requires a feasible initial control. With projection $`\Pi_{\mathcal U}`$,
the stopping measure and trial control are

$$
  r_{\Pi}=u-\Pi_{\mathcal U}(u-G^{-1}j_{h}'),\qquad
  u_{\alpha}=\Pi_{\mathcal U}(u-\alpha G^{-1}j_{h}').
$$

It stops when the declared metric norm of $`r_{\Pi}`$ reaches the configured
tolerance. A projected trial is accepted only when

$$
  j_{h}(u_{\alpha}) \leq
  j_{h}(u)+c\langle j_{h}',u_{\alpha}-u\rangle.
$$

The shared `ReducedSolverResultT` retains the final accepted
`ReducedEvaluationT`, accepted-objective and absolute/relative stopping-norm
histories, one nominal step length, metric step norm, and objective change for
each accepted iteration, accepted-iteration and line-search-trial counts,
separate state, adjoint, metric inverse-action, and Hessian-action counts, and
one stopping reason. The optional relative-gradient, objective-change, and
step tolerances are disabled when set to zero; enabled criteria report
`relative_gradient_tolerance`, `objective_change_tolerance`, or
`step_tolerance`, alongside `gradient_tolerance`, `maximum_iterations`, and
`line_search_failure`. Each objective trial is evaluated through
`ReducedDTOT::evaluate`, so both reported formulation solve counts increase
once for the initial point and once for every line-search trial. The metric
count records each direction-forming inverse metric action; the Hessian count
records each explicit provider application. Backend-specific inner iterations
remain in the metric or Hessian realization policy.

The default nonlinear-CG policy uses the metric gradient
$`g_{k}=G^{-1}j_{h}'(u_{k})`$ and the Polak–Ribière+ coefficient

$$
\beta_{k}=\max\left(0,
\frac{\langle j_{h}'(u_{k})-j_{h}'(u_{k-1}),g_{k}\rangle}
     {\langle j_{h}'(u_{k-1}),g_{k-1}\rangle}\right),
\qquad
d_{k}=-g_{k}+\beta_{k}d_{k-1}.
$$

The selectable `FletcherReevesDirectionPolicyT` instead uses

$$
\beta_{k}=\frac{\langle j_{h}'(u_{k}),g_{k}\rangle}
               {\langle j_{h}'(u_{k-1}),g_{k-1}\rangle},
\qquad
d_{k}=-g_{k}+\beta_{k}d_{k-1}.
$$

Both policies restart on the first direction, after the declared interval
(one full layout dimension by default), when the previous curvature
denominator is too small or non-finite, when the selected coefficient is
nonpositive, or when the resulting direction is not descending. Previous
covectors, gradients, and directions retain their layouts and are checked
before every update.

The selected limited-memory BFGS policy stores at most the configured number
of typed secant pairs. For a new iterate it forms the primal displacement
`$s_{k}=u_{k}-u_{k-1}$` and covector difference
`$y_{k}=j_{h}'(u_{k})-j_{h}'(u_{k-1})$`, accepts the pair only when
`$\langle y_{k},s_{k}\rangle$` is finite and larger than the configured
curvature tolerance, and applies the standard two-loop recursion. Its initial
inverse-Hessian action is the declared metric inverse, so the resulting
direction remains metric-aware without identifying covectors and primal
vectors. A rejected pair clears the secant history, returns the current
steepest direction, and increments the typed direction-reset count; the
policy also exposes whether the latest update was initial, accepted, or a
curvature reset. Layout mismatches are rejected rather than silently
discarded. The first registered `ReducedLimitedMemoryBfgsSolverT` integration
is the unconstrained one-state/one-decision reduced DTO with the mass metric.

The Newton policy requires a `ReducedHessianT` capability and solves
`$H(u)d=-j_{h}'(u)$` with metric-preconditioned inner conjugate gradients. It
rejects an absent capability, incompatible layouts, non-positive Hessian
curvature, and an unconverged inner solve. Each Hessian-vector application
and metric inverse action is returned in the shared solver report. The
reference linear-quadratic model supplies the exact action by solving the
tangent-state and incremental-adjoint systems; no dense reduced Hessian is
assembled.

The line-search protocol receives a typed current evaluation, a direction, a
trial-control builder, and a trial evaluator. The builder is responsible for
the caller's feasibility/projection policy; every acceptance condition then
uses the actual displacement `$u_{\mathrm{trial}}-u$`, not the nominal scalar
step times the unprojected direction. Armijo backtracks until sufficient
decrease holds, exact quadratic search uses
`$\tau=-j_{h}'[d]/\langle H d,d\rangle$` for positive curvature, and Wolfe
also checks the actual trial slope against its declared curvature fraction.
These policies return failure rather than silently accepting a non-descent or
non-finite trial.

`ReducedSearchSolverT` delegates trial construction and evaluation to the
selected line-search policy. Its default Armijo policy is initialized from
the legacy `ReducedSolverParameters` fields, while explicit policy instances
select combinations such as exact Newton or Wolfe steepest descent. Trial
counts are accumulated from the policy result and state/adjoint counts are
incremented by the evaluator callback, keeping reporting consistent across
all combinations.

The exact-search nonlinear-CG aliases combine the typed PR+ or
Fletcher–Reeves direction policy with the explicit positive-curvature
quadratic line search. For the linear-quadratic reference DTO, exact line
search makes the two coefficient updates equivalent in exact arithmetic; the
contract test compares their accepted objective histories and final controls.
The separate `QuadraticConjugateGradientDirectionPolicyT` keeps the classical
recurrence explicit: it uses the Fletcher–Reeves coefficient, periodically
restarts at the declared layout dimension, and rejects invalid curvature or a
non-descent recurrence instead of applying nonlinear-CG fallback behaviour.
`ReducedQuadraticConjugateGradientSolverT` composes it with exact quadratic
search.

`ReducedTrustRegionSolverT` is a separate globalization boundary rather than
another line-search policy. It forms the metric Cauchy step
$`s=-\tau G^{-1}j_{h}'`$, with $`\tau`$ clipped by the current radius, and
uses the matrix-free quadratic model

$$
m(s)=j_{h}(u)+\langle j_{h}',s\rangle+
\frac{1}{2}\langle H(u)s,s\rangle.
$$

The solver evaluates the actual reduction, computes the actual/predicted ratio,
accepts or rejects the trial against its threshold, and shrinks or expands the
radius according to the declared thresholds. Its report keeps one diagnostic
record per trial, including rejected trials. This first realization is
unconstrained and uses the explicit Hessian capability; projected trust-region
steps and a full Newton/truncated-CG trust-region subproblem remain separate
extensions.

## Metric and constraint boundary

The backend-neutral concrete metric supplied is a positive diagonal metric.
The serial deal.II backend also supplies `dealii_backend::MassMetric`, which
uses a one-block sparse mass matrix for its primal-to-dual action and a
CG inverse apply with explicit iteration and relative/absolute tolerance
parameters. The concrete dense and serial deal.II constraints are cellwise
boxes. The v0 control lowerer uses identifier `l2_cellwise`; the separately
owned v1 coefficient target uses `l2_cellwise_parameter`. Each lowerer creates
its box only for matching `FE_DGQ(0)` coefficients and an exact-layout scalar
or vector bound representation. This prevents accidental coefficient clipping
in continuous controls or an $H^{1}$ geometry.

## Reference model and verification

The reference linear-quadratic model implements

$$
  r(y,u)=Ay-f-Bu,
$$

$$
  J(y,u)=\tfrac{1}{2}(Cy-d)^{\mathsf T}W(Cy-d)
         +\tfrac{\alpha}{2}u^{\mathsf T}Ru.
$$

The `CTest` scenarios verify:

1. residual JVP/VJP pairing;
2. residual finite-difference derivative;
3. objective directional derivative;
4. state solve residual;
5. reduced DTO derivative;
6. metric inverse/apply consistency;
7. the selected cellwise $L^{2}$ box projection;
8. rejection of a non-diagonal metric that reuses a compatible metric's
   display identifier;
9. read-only block storage, rejected dimension-changing updates, and preserved
   pairing under checked algebraic updates;
10. exact `ContractError` messages for representative partition, callback,
   metric, constraint, and projected-solver precondition failures;
11. exact reduced-Hessian finite-difference and symmetry identities, plus the
    explicit-capability Newton convergence path;
12. exact quadratic, actual-displacement Armijo, and Wolfe line-search
    acceptance policies, including exact-Newton and Wolfe solver combinations;
13. dense and deal.II unconstrained/projected Armijo convergence, including
   active-bound and projected-stationarity checks; and
14. checked acceptance/rejection at the serial deal.II native-size boundary;
    and
15. detached owned reduced-service lifetime and typed state/adjoint solve
    reports, including sanitizer coverage in the backend-neutral profile.

This establishes the small executable algebra that a deal.II compiler must
produce. The first serial scalar diffusion-reaction compiler now exists; its
exact finite-element scope and exclusions are recorded in the
[deal.II v0 lowerer](dealii-lowerer.md).
