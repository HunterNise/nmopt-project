# Optimization reference

This reference describes how to consume `nmopt` formulation products with the
current reduced-search, trust-region, KKT, and PDAS interfaces.

For the mathematical background, use:

- [Reduced optimization methods](../manual/concepts/06-reduced-optimization-methods.md);
- [Optimality systems and KKT](../manual/concepts/07-optimality-systems-and-kkt.md);
- [Complementarity and PDAS](../manual/concepts/08-complementarity-and-pdas.md).

This page stays at the programming boundary: what object you need, what a
parameter means mathematically, which combinations are implemented, and what
result evidence to inspect.

## Start from the formulation product you actually have

The ordinary reduced path consumes:

```text
ReducedDTOT
   │
   ├── MetricT
   ├── optional ConstraintT
   └── optional ReducedHessianT
          │
          ▼
   reduced optimization solver
```

The DTO may come from the semantic compiler:

```cpp
const auto &compiled = *compilation_result.problem;

auto reduced = compiled.make_reduced_dto();
const auto &metric = compiled.metric();
```

or from a directly adapted PDE application:

```cpp
const auto &reduced = binding.reduced();
const auto &metric = binding.metric();
```

See [Compiler](compiler.md) and
[External deal.II solver integration](external-dealii-solver-integration.md)
for those producer paths.

Once the contracts are the same, optimization code does not need to know how
they were produced.

Quadratic KKT, supplied OTD, and PDAS are different products. They have their
own sections below rather than being forced through the reduced interface.

## Reduced derivative, metric gradient, and search direction

The reduced DTO returns a covector $j'(u)$. The metric $G$ identifies that
covector with a primal gradient $g$ through

$$
G g = j'(u).
$$

In code, the same operation is:

```cpp
const auto evaluation = reduced.evaluate(control);

const auto gradient =
  metric.inverse_apply(evaluation.reduced_derivative);
```

The metric norm used by the reduced solvers is

$$
\lVert g \rVert_{G}
=
\sqrt{\langle Gg,g\rangle}
=
\sqrt{\langle j'(u),g\rangle}.
$$

The helper used internally by the solvers performs exactly this conversion and
norm calculation.

For steepest descent, the search direction is

$$
d=-g.
$$

This distinction explains two API details that otherwise look arbitrary:

- `ReducedEvaluationT::reduced_derivative` is a **covector**;
- `ReducedDTOT::gradient_direction()` returns the metric gradient, not the
  negative search direction.

Direction policies own the sign and any further transformation of the
gradient.

## Working path: steepest descent with Armijo

The default reduced solver is steepest descent with Armijo backtracking:

```cpp
#include "nmopt/solvers/reduced_gradient.hpp"

using Backend = nmopt::dealii_backend::SerialBackend;
using Primal = nmopt::contract::PrimalBlockT<Backend>;
using Solver = nmopt::solvers::ReducedGradientSolverT<Backend>;
```

### Configure the outer solve

```cpp
nmopt::solvers::ReducedSolverParameters parameters;

parameters.maximum_iterations = 200;
parameters.maximum_line_search_trials = 25;

parameters.gradient_tolerance = 1.0e-8;
parameters.stopping_criterion =
  nmopt::solvers::ReducedStoppingCriterion::gradient_norm;

parameters.initial_step_length = 1.0;
parameters.minimum_step_length = 0.0;

parameters.armijo_fraction = 1.0e-4;
parameters.backtracking_factor = 0.5;
```

The defaults most relevant to this path are:

| Field | Default | Contract |
| --- | ---: | --- |
| `maximum_iterations` | `100` | Maximum accepted outer iterations |
| `maximum_line_search_trials` | `20` | Maximum Armijo trials per outer iteration |
| `gradient_tolerance` | `1e-8` | Absolute metric-gradient stopping tolerance |
| `initial_step_length` | `1.0` | First trial step parameter |
| `armijo_fraction` | `1e-4` | Armijo sufficient-decrease constant $c_{1}$ |
| `backtracking_factor` | `0.5` | Multiplicative reduction $\rho$ for rejected trials |

### Understand the Armijo fields

For an unconstrained straight-line trial

$$
u_{\mathrm{trial}}=u+\alpha d,
$$

Armijo accepts when

$$
j(u+\alpha d)
\le
j(u)+c_{1}\alpha j'(u)[d].
$$

The code fields are:

```text
initial_step_length  = initial α
armijo_fraction      = c1
backtracking_factor  = ρ
```

and each rejection updates

$$
\alpha \leftarrow \rho\alpha.
$$

The implementation records the actual update

$$
s=u_{\mathrm{trial}}-u
$$

and checks

$$
j(u_{\mathrm{trial}})
\le
j(u)+c_{1} j'(u)[s].
$$

That form matters for projected search, where the realized update need not be
exactly $\alpha d$ after projection.

`minimum_step_length` is a floor for Armijo backtracking. A value of zero
disables that floor.

### Build the initial control

The initial point is a primal in the control/metric layout:

```cpp
dealii::Vector<double> values(
  metric.layout()->dimension(0));
values = 0.0;

Primal initial_control(
  metric.layout(),
  {std::move(values)});
```

If a projection constraint is active, the initial control must already be
feasible. The solver does not silently project the starting point.

### Solve

```cpp
Solver solver(
  reduced,
  metric,
  parameters);

const auto result = solver.solve(initial_control);
```

The main evaluation loop is:

```text
u_k
 │
 ▼
ReducedDTOT::evaluate
 │  state + objective + adjoint + reduced derivative
 ▼
j'(u_k)
 │
 ▼
G^-1
 │
 ▼
g_k
 │
 ▼
direction policy
 │
 ▼
d_k
 │
 ▼
line-search trials
 │
 ▼
u_{k+1}
```

### Inspect the stopping reason and retained evaluation

```cpp
std::cout
  << nmopt::solvers::reduced_stopping_reason_name(
       result.stopping_reason)
  << '\n'
  << "accepted iterations: "
  << result.accepted_iterations
  << '\n';
```

The returned result keeps:

```cpp
result.control;
result.final_evaluation;
result.objective_history;
result.gradient_norm_history;
result.step_length_history;
result.iteration_records;
result.line_search_trials;
```

`maximum_iterations` and `line_search_failure` are valid termination states.
They are not convergence aliases. Application code should decide explicitly
which stopping reasons count as successful completion.

## Stopping criteria

`ReducedStoppingCriterion::automatic` preserves the compatibility behavior:
absolute gradient stopping is always active, and each positive optional
tolerance adds another stopping condition.

For absolute metric-gradient stopping,

$$
\lVert g_{k}\rVert_{G}
\le
\texttt{gradient\_tolerance}.
$$

For relative gradient stopping,

$$
\frac{\lVert g_{k}\rVert_{G}}{\lVert g_{0}\rVert_{G}}
\le
\texttt{relative\_gradient\_tolerance}.
$$

Configure it with:

```cpp
parameters.stopping_criterion =
  nmopt::solvers::ReducedStoppingCriterion::relative_gradient_norm;
parameters.relative_gradient_tolerance = 1.0e-6;
```

Objective-change stopping uses the accepted reduction

$$
j(u_{k})-j(u_{k+1})
\le
\texttt{objective\_change\_tolerance},
$$

while step stopping uses the metric norm

$$
\lVert u_{k+1}-u_{k}\rVert_{G}
\le
\texttt{step\_tolerance}.
$$

A positive optional tolerance enables the corresponding condition in automatic
mode; zero disables it.

`objective_target` is an independent absolute threshold on the objective.

The special `stationary` result is reserved for a metric-gradient norm at the
implementation's machine-scale stationarity threshold. It is distinct from
meeting a user-specified gradient tolerance.

## Direction policies

The line-search solver separates the direction policy from the globalization
policy. The main implemented directions are:

| Direction | Public policy | Additional data |
| --- | --- | --- |
| Steepest descent | `SteepestDescentDirectionPolicyT` | Metric |
| Nonlinear CG, PR+ | `NonlinearConjugateGradientDirectionPolicyT` | Previous derivative/gradient/direction |
| Fletcher–Reeves | `FletcherReevesDirectionPolicyT` | Previous derivative/gradient/direction |
| Classical quadratic CG | `QuadraticConjugateGradientDirectionPolicyT` | Positive gradient curvature |
| L-BFGS | `LimitedMemoryBfgsDirectionPolicyT` | Accepted secant history |
| Full BFGS | `FullBfgsDirectionPolicyT` | All accepted secant pairs |
| Newton | `NewtonDirectionPolicyT` | Explicit `ReducedHessianT` |

### Nonlinear conjugate gradient

The current default nonlinear-CG update is Polak–Ribière+ in the declared
metric.

Write

$$
g_{k}=G^{-1}j'(u_{k}).
$$

The denominator used by the implementation is

$$
\langle j'(u_{k-1}),g_{k-1}\rangle.
$$

For Fletcher–Reeves,

$$
\beta_{k}^{\mathrm{FR}}
=
\frac{\langle j'(u_{k}),g_{k}\rangle}
     {\langle j'(u_{k-1}),g_{k-1}\rangle},
$$

while the PR+ numerator is

$$
\langle j'(u_{k}),g_{k}\rangle
-
\langle j'(u_{k-1}),g_{k}\rangle,
$$

with the final coefficient clipped below by zero.

The policy fields are:

```cpp
nmopt::solvers::NonlinearConjugateGradientParameters cg;

cg.restart_interval = 0;
cg.curvature_tolerance = 1.0e-14;
```

`restart_interval == 0` selects one interval equal to the total coefficient
count of the current reduced layout.

`curvature_tolerance` guards the denominator. If the denominator is non-finite
or too small relative to its scale, the PR+/FR policies restart to steepest
descent rather than use an unstable coefficient.

```cpp
nmopt::solvers::NonlinearConjugateGradientDirectionPolicyT<Backend>
  direction(cg);

nmopt::solvers::ReducedConjugateGradientSolverT<Backend>
  solver(
    reduced,
    metric,
    parameters,
    direction);
```

### L-BFGS

For consecutive accepted controls, define the typed secant pair

$$
s_{k}=u_{k}-u_{k-1},
\qquad
 y_{k}=j'(u_{k})-j'(u_{k-1}).
$$

Here $s_{k}$ is primal and $y_{k}$ is a covector. The pair is accepted only when

$$
\langle y_{k},s_{k}\rangle
>
\texttt{curvature\_tolerance}.
$$

A failed curvature test clears the stored history and returns to a metric
steepest-descent direction.

```cpp
nmopt::solvers::LimitedMemoryBfgsParameters lbfgs;

lbfgs.memory_size = 10;
lbfgs.curvature_tolerance = 1.0e-14;
lbfgs.initial_inverse_hessian_scaling =
  nmopt::solvers::LimitedMemoryBfgsInitialScaling::metric_inverse;
```

`memory_size` is the maximum number of secant pairs retained by the two-loop
recursion.

With `metric_inverse`, the initial inverse-Hessian action is $G^{-1}$.

With `scalar_secant`, `nmopt` scales that metric inverse by

$$
\gamma_{k}
=
\frac{\langle y_{k},s_{k}\rangle}
     {\langle y_{k},G^{-1}y_{k}\rangle}.
$$

That is what the field `initial_inverse_hessian_scaling` selects; it is not a
Euclidean-only scalar rule hidden behind the metric abstraction.

```cpp
nmopt::solvers::LimitedMemoryBfgsDirectionPolicyT<Backend>
  direction(lbfgs);

nmopt::solvers::ReducedLimitedMemoryBfgsSolverT<Backend>
  solver(
    reduced,
    metric,
    parameters,
    direction);
```

The result's `direction_reset_count` is useful when diagnosing frequent
curvature failures.

### Full BFGS

Full BFGS uses the same typed secant pair and curvature test, but retains every
accepted pair:

```cpp
nmopt::solvers::FullBfgsParameters bfgs;
bfgs.curvature_tolerance = 1.0e-14;

nmopt::solvers::FullBfgsDirectionPolicyT<Backend>
  direction(bfgs);
```

Use it when retaining the full secant history is intentional; the public
interface otherwise matches the same reduced-search framework.

### Newton

Newton needs an explicit reduced Hessian. The first-order DTO does not imply
one.

The direction solves

$$
H(u_{k})d_{k}=-j'(u_{k}).
$$

The current policy uses metric-preconditioned CG. Its residual is

$$
r=H d+j'(u),
$$

and the reported residual norm is

$$
\lVert r\rVert_{G^{-1}}
=
\sqrt{\langle r,G^{-1}r\rangle}.
$$

The inner solve target is

$$
\lVert r\rVert_{G^{-1}}
\le
\max\left(
  \texttt{absolute\_tolerance},
  \texttt{relative\_tolerance}\lVert r_{0}\rVert_{G^{-1}}
\right).
$$

Configure it with:

```cpp
nmopt::solvers::ReducedNewtonParameters newton;

newton.maximum_inner_iterations = 100;
newton.relative_tolerance = 1.0e-10;
newton.absolute_tolerance = 1.0e-12;
newton.curvature_tolerance = 1.0e-14;

nmopt::solvers::NewtonDirectionPolicyT<Backend>
  direction(hessian, newton);
```

`curvature_tolerance` checks the Hessian curvature used by inner CG. For an
inner search direction $p$, `nmopt` requires

$$
\langle Hp,p\rangle
>
\texttt{curvature\_tolerance}
\lVert p\rVert_{G}^{2}.
$$

A default-constructed Newton policy has no Hessian capability and cannot
produce a Newton direction.

The solver result retains `hessian_solve_history`, `hessian_action_count`, and
`metric_solve_count`, so inner work is not hidden inside the outer iteration
count.

## Line-search policies

### Armijo

Armijo was described in the working path. The explicit policy record is:

```cpp
nmopt::solvers::ArmijoLineSearchParameters line;

line.maximum_trials = 30;
line.initial_step_length = 1.0;
line.armijo_fraction = 1.0e-4;
line.backtracking_factor = 0.5;
line.minimum_step_length = 0.0;
```

Use an explicit policy when line-search settings should be independent of the
compatibility fields in `ReducedSolverParameters`.

### Fixed step

A fixed-step policy declares one positive step parameter:

```cpp
nmopt::solvers::FixedStepLineSearchParameters line;
line.step_length = 0.25;

nmopt::solvers::FixedStepLineSearchPolicyT<Backend>
  line_search(line);
```

It evaluates one trial and accepts a finite objective. There is no Armijo or
curvature predicate.

For a straight-line unprojected update,

$$
u_{k+1}=u_{k}+\alpha d_{k},
\qquad
\alpha=\texttt{step\_length}.
$$

### Strong Wolfe

The strong-Wolfe policy has:

```cpp
nmopt::solvers::WolfeLineSearchParameters wolfe;

wolfe.maximum_trials = 30;
wolfe.initial_step_length = 1.0;
wolfe.backtracking_factor = 0.5;
wolfe.sufficient_decrease_fraction = 1.0e-4;
wolfe.curvature_fraction = 0.9;
```

With trial update $s=u_{\mathrm{trial}}-u$, the sufficient-decrease condition
is

$$
j(u_{\mathrm{trial}})
\le
j(u)+c_{1} j'(u)[s],
$$

where

```text
sufficient_decrease_fraction = c1.
```

The strong curvature condition is

$$
\left|j'(u_{\mathrm{trial}})[s]\right|
\le
c_{2}\left|j'(u)[s]\right|,
$$

where

```text
curvature_fraction = c2.
```

The implementation requires

$$
0<c_{1}<c_{2}<1.
$$

Because Wolfe evaluates the trial derivative, derivative-bearing trials can
incur adjoint solves. That cost difference is visible in the retained solve
counts.

### Weak Wolfe

The weak-Wolfe policy uses the same parameter record, but replaces the strong
curvature condition with

$$
j'(u_{\mathrm{trial}})[s]
\ge
c_{2} j'(u)[s].
$$

Construct it as:

```cpp
nmopt::solvers::WeakWolfeLineSearchPolicyT<Backend>
  line_search(wolfe);
```

### Exact quadratic line search

For a quadratic reduced model and a direction $d$, the exact step is

$$
\alpha_{\ast}
=
-
\frac{j'(u)[d]}
     {\langle H d,d\rangle}.
$$

The policy therefore requires the same explicit reduced Hessian capability as
Newton:

```cpp
nmopt::solvers::ExactQuadraticLineSearchParameters exact;
exact.curvature_tolerance = 1.0e-14;
exact.objective_tolerance = 1.0e-12;

nmopt::solvers::ExactQuadraticLineSearchPolicyT<Backend>
  line_search(hessian, exact);
```

`curvature_tolerance` requires

$$
\langle Hd,d\rangle
>
\texttt{curvature\_tolerance}.
$$

After taking the exact straight-line step, the implementation still checks
that the trial objective has not increased beyond

$$
j(u)
+
\texttt{objective\_tolerance}
\max(1,|j(u)|).
$$

That second field is therefore a numerical acceptance tolerance, not another
curvature parameter.

A typical exact-Newton composition is:

```cpp
nmopt::solvers::NewtonDirectionPolicyT<Backend>
  direction(hessian, newton);

nmopt::solvers::ExactQuadraticLineSearchPolicyT<Backend>
  line_search(hessian, exact);

nmopt::solvers::ReducedExactNewtonSolverT<Backend>
  solver(
    reduced,
    metric,
    parameters,
    direction,
    line_search);
```

## Projected reduced search

A projection constraint is supplied separately from the DTO:

```cpp
const auto *constraint = compiled.constraint();
if (constraint == nullptr)
  throw std::runtime_error("projection constraint required");

nmopt::solvers::ReducedGradientSolverT<Backend>
  solver(
    reduced,
    metric,
    *constraint,
    parameters);
```

The constraint must support projection in the **actual metric realization**:

```cpp
constraint->supports_projection_in(metric);
```

Layout compatibility alone is not enough.

For a direction $d$ and step parameter $\alpha$, the trial builder has the
form

$$
u_{\mathrm{trial}}
=
P_{G}(u+\alpha d),
$$

where $P_{G}$ is the constraint projection in the declared metric. That is why
Armijo/Wolfe acceptance uses the actual update

$$
s=u_{\mathrm{trial}}-u
$$

rather than assuming $s=\alpha d$.

The current projected line-search combinations are deliberately limited to
steepest descent with:

```text
Armijo
fixed step
strong Wolfe
weak Wolfe
```

Passing a projection constraint to nonlinear CG, BFGS, Newton, or exact
quadratic search is currently rejected. Trust region is also unconstrained.

For compiler-produced cellwise boxes, use the compiled constraint. It carries
the metric-realization witness and box-data identity required to keep
projection semantics consistent with the lowered problem.

## Trust-region optimization

Trust region consumes:

```text
ReducedDTOT
MetricT
ReducedHessianT
```

and is a separate solver rather than a line-search policy.

The local quadratic model is

$$
m_{k}(s)
=
j(u_{k})
+
j'(u_{k})[s]
+
\frac12\langle H_{k} s,s\rangle,
$$

with metric trust region

$$
\lVert s\rVert_{G}\le\Delta_{k}.
$$

Configure it with:

```cpp
nmopt::solvers::ReducedTrustRegionParameters trust;

trust.maximum_iterations = 100;
trust.maximum_trials_per_iteration = 10;
trust.gradient_tolerance = 1.0e-8;

trust.initial_radius = 1.0;
trust.minimum_radius = 1.0e-12;
trust.maximum_radius = 1.0e6;

trust.acceptance_threshold = 0.1;
trust.shrink_threshold = 0.25;
trust.expansion_threshold = 0.75;
trust.shrink_factor = 0.25;
trust.expansion_factor = 2.0;
```

For a trial step $s$, `nmopt` records

$$
\text{ared}
=
j(u_{k})-j(u_{k}+s),
$$

$$
\text{pred}
=
m_{k}(0)-m_{k}(s),
$$

and

$$
\rho
=
\frac{\text{ared}}{\text{pred}}.
$$

The field mapping is direct:

```text
acceptance_threshold  -> accept when ρ >= this value
shrink_threshold      -> accepted step shrinks radius when ρ is below this
expansion_threshold   -> expand when ρ is above this and step is near boundary
shrink_factor         -> Δ <- factor * Δ
expansion_factor      -> Δ <- factor * Δ
```

Rejected trials always shrink the radius.

The default subproblem method is the Cauchy step. For truncated CG:

```cpp
trust.subproblem_method =
  nmopt::solvers::ReducedTrustRegionSubproblemMethod::
    truncated_conjugate_gradient;

trust.maximum_subproblem_iterations = 100;
trust.subproblem_relative_tolerance = 1.0e-10;
trust.subproblem_absolute_tolerance = 1.0e-12;
```

The truncated-CG residual uses the same metric-preconditioned norm convention
as the Newton inner solve.

The trust-region result retains the radius, predicted/actual reduction, ratio,
subproblem status, subproblem iteration count, and acceptance decision for
every trial.

## Work accounting and accepted-iteration evidence

Reduced results retain work explicitly:

```cpp
result.state_solve_count;
result.adjoint_solve_count;
result.metric_solve_count;
result.hessian_action_count;
result.direction_reset_count;
```

Line-search results also retain:

```cpp
result.iteration_records;
result.line_search_trials;
```

An accepted-iteration record carries objective before/after, actual change,
requested step parameter, realized step norm, descent pairing, stationarity,
per-iteration work, cumulative work, accepted evaluation, and policy-specific
acceptance evidence.

This is the correct source when comparing methods. Objective-history length
alone cannot tell you how many rejected value trials, derivative trials,
adjoint solves, or Hessian actions were required.

## Quadratic KKT products

A quadratic KKT product is not a reduced DTO.

For a problem of the form

$$
\min_{x}
\frac12\langle Qx,x\rangle
-
\langle b,x\rangle
\quad\text{subject to}\quad
Dx=c,
$$

the KKT equations are represented abstractly as

$$
\begin{bmatrix}
Q & D^{\ast} \\
D & 0
\end{bmatrix}
\begin{bmatrix}
x\\
\lambda
\end{bmatrix}
=
\begin{bmatrix}
b\\
c
\end{bmatrix}.
$$

The `nmopt` product stores typed primal/multiplier layouts, the $Q$, $D$,
$D^{\ast}$ and transpose actions, right-hand sides, pairings, multiplier
conversion, and declared structural assumptions.

From the compiler:

```cpp
const auto &compiled_kkt =
  *compilation_result.kkt_problem;

const auto &product =
  compiled_kkt.product();
```

### Select the iterative policy

```cpp
nmopt::contract::QuadraticKKTSolverPolicy policy;

policy.method =
  nmopt::contract::QuadraticKKTSolverMethod::minres;

policy.maximum_iterations = 0;
policy.relative_tolerance = 1.0e-10;
policy.absolute_tolerance = 1.0e-12;
```

`MINRES` is accepted only when the product declares the required symmetric
indefinite capability, complete pairings, and transpose-consistency evidence.

`GMRES` is the nonsymmetric-capable alternative and adds the
`gmres_maximum_basis` policy field.

Validate the product/policy combination before solving:

```cpp
nmopt::contract::validate(product, policy);
```

### Serial deal.II solve

For the current serial deal.II backend:

```cpp
#include "nmopt/dealii/serial_kkt_solver.hpp"

const auto kkt_result =
  nmopt::dealii_backend::solve_serial_quadratic_kkt(
    product,
    policy);
```

A zero `maximum_iterations` selects the adapter's dimension-based automatic
limit; it does not mean zero iterations.

The returned report contains the linear-solve report plus stationarity and
equality residuals.

For another backend or a different preconditioner, implement an adapter that
returns the same `QuadraticKKTSolveResultT<Backend>` contract. The typed KKT
product remains the source of operator semantics; packing, preconditioning,
and native iterative algebra belong in the backend adapter.

## PDAS products

The compiled PDAS product combines:

```text
quadratic KKT product
        +
box complementarity
        +
metric-backed multiplier representation
        +
shared box data
        +
PDAS / KKT policies
```

```cpp
const auto &compiled_pdas =
  *compilation_result.pdas_problem;

const auto &product = compiled_pdas.product();
const auto &complementarity = compiled_pdas.complementarity();
const auto &metric = compiled_pdas.metric();
const auto &box = compiled_pdas.box_data();

const auto &pdas_policy = compiled_pdas.pdas_policy();
const auto &kkt_policy = compiled_pdas.kkt_solver_policy();
```

The multiplier in complementarity is a covector. Classification obtains the
corresponding primal representative through the declared metric; it is not
correct to compare raw dual coefficients to primal bound distances as if all
representations were Euclidean-identical.

For the serial deal.II KKT subproblem:

```cpp
#include "nmopt/dealii/serial_kkt_solver.hpp"

auto solver =
  compiled_pdas.make_solver(
    [&](const auto &active_product) {
      return
        nmopt::dealii_backend::solve_serial_quadratic_kkt(
          active_product,
          kkt_policy);
    });
```

The initial KKT point must use the product layouts and the initial control must
already satisfy the box. The initial box multiplier must use the
complementarity layout.

The solver then returns:

```cpp
result.solution;
result.box_multiplier;
result.iterations;
result.stopping_reason;
```

with per-iteration active-set, primal/dual feasibility, complementarity,
stationarity, equality, and inner-KKT evidence.

The current stopping reasons are:

```text
converged
maximum_iterations
kkt_solve_failed
```

Active-set stability by itself is therefore not the complete convergence
criterion.

For direct construction of the current serial cellwise-$L^{2}$
complementarity object, see
[`SerialCellwiseBoxComplementarity`](../../include/nmopt/dealii/serial_pdas.hpp).
When the compiler already produced a PDAS bundle, use that bundle rather than
reconstructing its metric/box identity separately.

## Supplied OTD

A supplied optimize-then-discretize product is also distinct from the reduced
DTO path:

```cpp
const auto &compiled_otd =
  *compilation_result.supplied_otd_problem;

const auto &system = compiled_otd.system();
```

The supplied system exposes its all-at-once state/adjoint/control equations and
its own solve service.

If the semantic supplied-OTD formulation is explicitly compiled to a
compatible quadratic KKT product, consume the resulting KKT product through
the KKT path instead. There is no `CompilationProduct::supplied_otd` selector;
the semantic formulation determines that supplied-OTD product.

## Choosing the public path

Use the formulation product and available capabilities to choose the next
object:

```text
ReducedDTOT + MetricT
    │
    ├── first-order only
    │      ├── steepest descent
    │      ├── nonlinear CG
    │      ├── L-BFGS
    │      └── full BFGS
    │
    ├── + projection ConstraintT
    │      └── projected steepest descent
    │          with Armijo/fixed/Wolfe
    │
    └── + ReducedHessianT
           ├── Newton + line search
           └── unconstrained trust region

QuadraticKKTProduct
    └── backend KKT solve adapter

CompiledPDASProblemT
    └── active-set solver + KKT subproblem adapter

CompiledSuppliedOTDProblemT
    └── supplied all-at-once system service
```

This routing is more reliable than choosing a solver name first and then
trying to coerce an incompatible formulation product into it.

Optimization ends with a numerical result and solver evidence. Run-set
expansion, parameter-file precedence, artifact naming, and plotting policy are
handled by [Application execution](application-execution.md) and
[Parameter files](parameter-files.md).
