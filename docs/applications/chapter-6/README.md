# Chapter 6 application and benchmark scenarios

Chapter 6 scenarios freeze a configuration of a Chapter 5 recipe. They are
not a second problem library and they do not contain PDE assembly or solver
implementations. The [benchmark roadmap](../../planning/chapter-6-benchmark-suite-roadmap.md)
owns selection and status, the [numerical-examples reference](../../guides/chapter-6-numerical-examples.md)
owns source equations and data, and the [numerical-methods guide](../../guides/chapter-6-numerical-methods.md)
owns reusable solver contracts.

The public execution path is:

```text
typed scenario
  -> Chapter 5 problem recipe
  -> ProblemSpec + runtime data bindings
  -> owned deal.II compilation session
  -> selected reduced/KKT/PDAS product
  -> solver report + compilation manifest
  -> detached experiment envelope and benchmark artifacts
```

## Scenario contract

Every frozen scenario must record these fields, even when a source example
left one of them unspecified:

| Group | Frozen contents |
| --- | --- |
| Identity | Scenario ID, recipe ID/version, source section/equation, and source revision. |
| Geometry | Dimension, domain, mesh generator or input, refinement, boundary labels, and material IDs. |
| Semantic choices | State/test/control/observation spaces, residual terms, observation region, loss, metric, constraint, and formulation product. |
| Runtime data | Forcing, target, coefficients, fixed data, transport/weight data, bounds, and provenance strings. |
| Discretization | Finite-element degrees, quadrature/evaluation policies, execution realization, and boundary orientation. |
| Algorithm | Initial control, solver family, direction, line search or fixed step, iteration limits, and stopping tolerances. |
| Evidence | Objective components, state/adjoint/metric/KKT histories, feasibility, reduced-gradient or complementarity records, and run environment. |

If the source does not specify a value, the scenario must name a recovered
value or a manufactured replacement and explain why. It must not silently use
a library default when that changes the mathematical or numerical comparison.

## Activation matrix

The current benchmark order is deliberately narrow:

| ID | Source example | Recipe family | Decision | Product/method |
| --- | --- | --- | --- | --- |
| B0 | Common harness | Any selected recipe | Required harness preparation | Deterministic manifest/report artifact path |
| B1 | E6.5.1 distributed Laplace control | `scalar-diffusion-reaction-volume` | Selected first benchmark | Reduced DTO; steepest descent versus L-BFGS |
| B2 | E6.5.2 Graetz-flow boundary control | `scalar-neumann-convection-subdomain` | Selected second benchmark | Reduced DTO; BFGS or declared recovered Armijo/fixed-step policy |
| B3 | E6.9.1 symmetric box Laplace | Distributed scalar with cellwise box | Desirable follow-up | Scalar KKT/PDAS; cellwise-discontinuous control |
| B4 | E6.9.2 asymmetric box Laplace | Distributed scalar with spatially varying cellwise bounds | Desirable follow-up | Reuse B3 KKT/PDAS product |
| B5 | E6.7.1 all-at-once Laplace | Distributed scalar | Conditional | Quadratic KKT and optional preconditioning |
| B6 | E6.7.2 diffusion-reaction all-at-once | Reuse B5 recipe/scenario shape | Follow-up | Reuse B5 implementation |

B1 and B2 are the current benchmark activation targets after the minimum
recipe boundary and B0 harness exist. B3/B4 are not prerequisites for those
reduced-space benchmarks. B5/B6 do not activate automatically; preconditioning
is conditional on measured need. Stokes, stabilization comparisons, generic
continuous-control boxes, measure-valued state constraints, and automatic OtD
derivation remain excluded.

## B0 — common harness

The harness consumes a frozen scenario and the public compiled product. It may
measure time and memory, serialize reports, and select output fields, but it
must not contain a second lowerer or optimization loop.

At minimum, B0 must produce:

- a deterministic scenario/recipe identity record;
- the `CompilationManifest` and complete diagnostic report;
- the solver policy snapshot and typed solver report;
- objective, tracking, gradient/KKT, feasibility, and solve-count histories;
- mesh, parameter, and artifact provenance; and
- a caller-supplied `RunEnvironmentRecord`.

The in-memory
`experiment::ReducedSearchExperimentEnvelopeT<Backend>` is the association
boundary for reduced runs. It owns values only; artifact serialization is a
harness concern. `application::benchmark::BenchmarkHarnessT<Scenario>` adds
the deterministic scenario identity, complete validation diagnostics,
measurements, and selected output fields around that envelope. It deliberately
does not execute a solver or write files. Use `release-dealii` for reproduction
runs and smaller development meshes for local iteration.

`application::benchmark::BenchmarkArtifactWriter` renders the captured record
as deterministic escaped `key=value` lines to a caller-owned stream. The
`application::benchmark::HeadlessBenchmarkRunnerT<Scenario>` connects the
remaining steps: it calls a problem builder with `scenario.problem`, passes the
resulting `ProblemSpec` and scenario to the execution adapter, captures runner
wall time, finalizes the harness artifact, and renders it with the writer.
The execution adapter supplies backend compilation, solver reports, detached
experiment data, and non-wall-clock measurements. Path selection and scenario
artifact-directory creation remain outside this API.

## B1 — distributed Laplace control

### Scenario selection

Use the distributed scalar recipe with homogeneous Dirichlet data, full-domain
$L^2$ state tracking, a cellwise `FE_DGQ(0)` volume control, and the positive
cellwise $L^2$ metric. The source uses continuous `P1` control coordinates;
the first framework-native scenario must record the chosen discrete control
representation explicitly rather than claiming source discretization parity
without evidence.

The source target is

```text
z_d(x) = 10 x_1(1-x_1) x_2(1-x_2),  x in (0,1)^2.
```

The source does not fully specify `f`. Select either the recovered source
forcing or a manufactured replacement and store that choice in the scenario
record. Absolute objective values from different choices are not comparable.
The default `make_b1_scenario()` selects the explicitly named manufactured-zero
forcing replacement; pass `B1ForcingSelection::recovered_source` only when the
recovered function and provenance are available.

The semantic graph is assembled through the Chapter 5 recipe rather than by
repeating graph construction in the benchmark adapter:

```cpp
const auto scenario =
  nmopt::application::chapter6::make_b1_scenario(method);
const auto specification =
  nmopt::application::chapter6::make_b1_problem_spec(scenario);
```

`make_b1_problem_spec(...)` only returns the backend-neutral `ProblemSpec`.
The execution adapter still supplies forcing and target functions, scalar
coefficients, the owned mesh session, and the selected solver product.

### Frozen inputs

| Input | Scenario value |
| --- | --- |
| Domain | $(0,1)^2$ |
| State equation | $-\Delta y=f+u$, homogeneous Dirichlet boundary |
| Forcing | `B1ForcingSelection::manufactured_zero` by default, or an explicitly named recovered source; provenance is mandatory |
| Diffusion/reaction | Diffusion `1`, reaction `0` unless the recovered source says otherwise |
| Regularization sweep | $\beta=10^{-1},10^{-2},10^{-3}$; retain $10^{-6}$ only for the source field illustration |
| Initial control | Zero control in the selected control layout |
| Methods | Steepest descent and L-BFGS from the same initial control |
| Source line search | Armijo $\rho=0.7$, $\sigma=10^{-5}$, at most five rescalings; encode any unavailable minimum-step rule explicitly |
| Source stopping | Steepest-descent tolerance $10^{-3}$; the current factory declares $10^{-8}$ for L-BFGS and records it separately from the source rule |
| Mesh | Source square/triangular mesh when available; development refinement must be labeled as such |

The current public `ReducedSolverParameters` defaults do not equal all source
line-search values. A B1 scenario must therefore supply the intended solver
policy as typed options and preserve the policy snapshot in its artifact.

### Required evidence

- objective and tracking reduction for every regularization value;
- iteration, state-solve, adjoint-solve, metric-solve, and line-search counts;
- reduced-gradient histories and terminal stopping reason;
- finite-difference verification of the linear-quadratic Hessian action; and
- the qualitative source trend that L-BFGS is substantially faster than
  steepest descent, without treating runtime as a portable correctness value.

### B1 backend execution adapter

The deal.II-specific adapter is declared in
`include/nmopt/application/dealii/chapter6_b1.hpp`. For the manufactured
forcing choice, the complete assembly is:

```cpp
#include "nmopt/application/dealii/chapter6_b1.hpp"

const auto scenario =
  nmopt::application::chapter6::make_b1_scenario(method);
const auto specification =
  nmopt::application::chapter6::make_b1_problem_spec(scenario);

nmopt::application::chapter6::dealii::B1ManufacturedDataT<2> data;
const auto runtime =
  nmopt::application::chapter6::dealii::make_b1_manufactured_runtime_data(
    scenario, data);
const auto session =
  nmopt::application::chapter6::dealii::make_b1_compilation_session<2>(
    scenario);
nmopt::application::chapter6::dealii::B1ReducedExecutionAdapterT<2> execute{
  beta, runtime, session, environment};

using Runner =
  nmopt::application::benchmark::HeadlessBenchmarkRunnerT<
    decltype(scenario)>;
const auto result = Runner(scenario).run(
  [](const auto &parameters) {
    return nmopt::application::chapter6::make_b1_problem_spec(parameters);
  },
  execute);
```

The adapter binds forcing, the polynomial desired state, diffusion, reaction,
and the per-run regularization value; creates the owned square-domain mesh
session; compiles the assembled reduced DTO; and dispatches either steepest
descent or limited-memory BFGS. The returned detached envelope contains the
compiler manifest, typed solver report, policy snapshot, and environment.
Additional artifact fields record the B1 method, regularization, histories,
stopping reason, and solve counts. Recovered forcing uses a caller-owned
`B1RuntimeDataT<dim>` instead of `B1ManufacturedDataT<dim>`.

## B2 — Graetz-flow boundary control

### Scenario selection

Use the `scalar-neumann-convection-subdomain` recipe, seeded by
`make_neumann_convection_subdomain_tracking_problem(observed_material_id,
with_facewise_box=false)`. The public B2 helper adds the nonzero fixed-data
port and lifting policy before compilation, so an adapter only needs to
assemble the typed scenario and bind the declared runtime functions:

```cpp
const auto scenario =
  nmopt::application::chapter6::make_b2_scenario(
    nmopt::application::chapter6::GraetzCase::observation_wings_constant_target);
const auto specification =
  nmopt::application::chapter6::make_b2_problem_spec(scenario);
```

The resulting state has scalar conservative transport, fixed temperature on
the declared Dirichlet boundary, controlled Neumann flux on a marked
boundary, insulated outflow, and volume observation on a downstream material
region. The control is facewise and the metric is `l2_facewise`. The helper
declares `fixed_dirichlet_data` as a `fixed_dirichlet_lifting` function and
connects it through `fixed_dirichlet_reconstruction`; passing that function
to a graph without the declared port remains an error.

The four scenarios are the Cartesian product of two observation regions and
two targets:

| Case | Observation region | Target |
| --- | --- | --- |
| a | $x_1>1$, with $x_2<0.3$ or $x_2>0.7$ | $z_d=2$ |
| b | $x_1>1$ | $z_d=2$ |
| c | $x_1>1$, with $x_2<0.3$ or $x_2>0.7$ | $z_d=4x_2(1-x_2)$ |
| d | $x_1>1$ | $z_d=4x_2(1-x_2)$ |

The material-ID realization of these regions must be frozen in the mesh
record. It must not be replaced by a runtime geometric predicate that is
absent from the semantic graph.

The four case records are discoverable through `make_catalog()`. The default
case keeps the stable ID `chapter-6.b2.graetz-flow`; the other cases append
`full-constant`, `wings-parabolic`, or `full-parabolic` to that ID.

### Frozen inputs and bindings

| Input | Scenario value |
| --- | --- |
| Domain | $(0,1+l)\times(0,1)$ with $l=3$ |
| Diffusion | $\mu=0.1$ |
| Transport | $b(x)=(1.5x_2(1-x_2),0)$, bound through `conservative_transport` |
| Reaction | `0` for the stated Graetz equation |
| Regularization | $\beta=10^{-3}$ |
| Boundary data | Fixed temperature `1` on `dirichlet_boundary`; zero natural outflow on the declared outflow region |
| Control | Facewise Neumann flux on `control_boundary`; no facewise box for B2 |
| Initial control | Zero facewise control |
| Method | Source BFGS; if a fixed step is recovered, record it; otherwise use a declared Armijo policy |
| Discretization | Source `P1` triangular mesh when available; development refinement is labeled |

The runtime binding must include forcing, desired state, diffusion, reaction,
regularization, fixed Dirichlet data, conservative transport, and provenance
for each. The graph and mesh must agree on fixed, controlled, outflow, and
observation labels. A transport field must not be passed through the general
scalar bundle with a different semantic role.

### Required evidence

- state/control dimensions and the compiled manifest;
- transport, boundary partition, and inflow/outflow policy;
- residual JVP/VJP and reduced-Taylor checks for the boundary control;
- objective and relative-gradient reduction for all four cases; and
- state/control field behavior when observation region or target changes.

The first reproduction is the stated Galerkin formulation. GLS and other
stabilization policies are outside B2. If the unstabilized discretization is
inadequate, record the limitation instead of silently changing the scenario.

## Later scenario contracts

B3 and B4 must use cellwise-discontinuous control for the first framework
native PDAS runs. The coefficientwise box, positive-diagonal `l2_cellwise`
metric, KKT product, PDAS complementarity, and active-set evidence must all
refer to the same control layout. A continuous `Q1` box reproduction is a
separate feature and is not implied by B3.

B5 and B6 are all-at-once products. They must preserve the source's objective
scaling, especially the $\beta\lVert u\rVert^2$ convention in E6.7.1, and must
record multiplier conversion, KKT residuals, and inner/outer work. Do not
activate a preconditioner merely because the source used one; measure whether
the selected scalar run needs it first.

## Scenario handoff checklist

Before a benchmark run is enabled, an agent should be able to answer without
reading compiler code:

1. Which recipe ID/version creates the graph?
2. Which mesh labels become semantic region IDs and material IDs?
3. Which runtime binding ports are required and what are their provenance
   strings?
4. Which control layout, loss, metric, bounds, and product are selected?
5. Which solver policy and stopping rule reproduce the stated comparison?
6. Which source omissions were recovered or manufactured?
7. Which evidence fields and artifacts determine pass/fail versus merely
   provide a performance observation?
