# V0 executable contract

## Scope

This is the first code-level realization of the default policies in the
[implementation-readiness review](implementation-readiness-review.md). It is
backend-neutral on purpose: a future deal.II compiler lowers a semantic
problem to this contract, while the contract itself does not contain a
DoFHandler, finite element, matrix type, or solver selection.

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
| `ReducedDTO` | The state–adjoint–reduced-covector workflow for one state and one control block. |
| `solvers::ReducedGradientSolverT` | Backend-parametric, unconstrained reduced Armijo method consuming only `ReducedDTOT` and `MetricT`. |

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
two variable blocks: one eliminated state and one decision/control
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
They will extend, rather than alter, the value/JVP/VJP core.

## First reduced solver

`solvers::ReducedGradientSolverT` is the first generic optimizer. Given the
DTO covector $`j_{h}'`$ and a declared metric $`G`$, it forms
$`g=G^{-1}j_{h}'`$ and attempts the unconstrained update
$`u^{+}=u-\alpha g`$. A trial is accepted only when it satisfies

$$
  j_{h}(u-\alpha g) \leq j_{h}(u)-c\alpha\langle j_{h}',g\rangle.
$$

The solver returns accepted-objective and metric-gradient-norm histories,
accepted-iteration and line-search-trial counts, state/adjoint solve counts,
and one stopping reason: `gradient_tolerance`, `maximum_iterations`, or
`line_search_failure`. Each objective trial is evaluated through
`ReducedDTOT::evaluate`, so both reported solve counts increase once for the
initial point and once for every line-search trial. It implements no
constraint projection; the cellwise box extension owns that behavior.

## Metric and constraint boundary

The backend-neutral concrete metric supplied is a positive diagonal metric.
The serial deal.II backend also supplies `dealii_backend::MassMetric`, which
uses a one-block sparse mass matrix for its primal-to-dual action and a
CG inverse apply with explicit iteration and relative/absolute tolerance
parameters. The only concrete constraint supplied is a cellwise box
constraint, and it explicitly supports projection only in a metric whose
identifier is `l2_cellwise`. This models the selected `FE_DGQ(0)`
volume-control policy and prevents accidental coefficient clipping in an
$H^{1}$ geometry.

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
7. the selected cellwise $L^{2}$ box projection.

This establishes the small executable algebra that a deal.II compiler must
produce. The first serial scalar diffusion-reaction compiler now exists; its
exact finite-element scope and exclusions are recorded in the
[deal.II v0 lowerer](dealii-v0-lowerer.md).
