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
| `PrimalBlock` | A typed block of coefficients in a declared discrete primal space. |
| `CovectorBlock` | A typed block of coefficients in the dual of that layout. |
| `pair` | The declared dual-coefficient pairing. It is the only primitive that pairs a covector with a primal vector. |
| `ExecutableModel` | Residual value, residual JVP, residual VJP, objective value, and objective derivative. |
| `Metric` | Explicit primal-to-dual action and inverse action. The serial deal.II backend supplies `dealii_backend::MassMetric` for one-block sparse mass matrices. |
| `Constraint` | Feasibility and a metric-specific projection capability. |
| `ReducedDTO` | The state–adjoint–reduced-covector workflow for one state and one binary decision block. |
| `solvers::ReducedGradientSolverT` | Backend-parametric unconstrained or projected reduced Armijo method consuming `ReducedDTOT`, `MetricT`, and an optional `ConstraintT`. |

The unsuffixed public aliases select the dense reference backend. The
corresponding types with a T suffix are backend-parametric, for example
`ExecutableModelT`, `PrimalBlockT`, `CovectorBlockT`, `MetricT`, and `ReducedDTOT`.

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
and tolerance policy. The constructor validates all block-layout connections;
no implicit state/control selection is made.

Mixed partitions, multiple equation blocks, nonlinear all-at-once Newton,
OTD, and Hessian-vector actions are intentionally outside this v0 contract.
The historical `StateControlPartitionT` name denotes this binary decision
port; a v1 coefficient lowerer may give its layout the semantic identifier
`parameter` without changing the backend-neutral value/JVP/VJP core.

## First reduced solver

`solvers::ReducedGradientSolverT` is the first generic optimizer. Given the
DTO covector $`j_{h}'`$ and a declared metric $`G`$, it forms
$`g=G^{-1}j_{h}'`$ and attempts the unconstrained update
$`u^{+}=u-\alpha g`$. A trial is accepted only when it satisfies

$$
  j_{h}(u-\alpha g) \leq j_{h}(u)-c\alpha\langle j_{h}',g\rangle.
$$

The constraint-qualified constructor accepts a `ConstraintT` only when its
layout matches the metric and it declares projection support for that metric.
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

The solver returns accepted-objective and stopping-norm histories,
accepted-iteration and line-search-trial counts, state/adjoint solve counts,
and one stopping reason: `gradient_tolerance`, `maximum_iterations`, or
`line_search_failure`. Each objective trial is evaluated through
`ReducedDTOT::evaluate`, so both reported solve counts increase once for the
initial point and once for every line-search trial.

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

The `CTest` executable verifies:

1. residual JVP/VJP pairing;
2. residual finite-difference derivative;
3. objective directional derivative;
4. state solve residual;
5. reduced DTO derivative;
6. metric inverse/apply consistency; and
7. the selected cellwise $L^{2}$ box projection;
8. exact `ContractError` messages for representative partition, callback,
   metric, constraint, and projected-solver precondition failures; and
9. dense and deal.II unconstrained/projected Armijo convergence, including
   active-bound and projected-stationarity checks.

This establishes the small executable algebra that a deal.II compiler must
produce. The first serial scalar diffusion-reaction compiler now exists; its
exact finite-element scope and exclusions are recorded in the
[deal.II v0 lowerer](dealii-lowerer.md).
