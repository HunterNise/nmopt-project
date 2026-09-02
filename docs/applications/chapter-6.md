# Chapter 6 application scenarios

This document catalogues Chapter 6 application scenarios, each of which
describes how a Chapter 5 recipe is configured. It is not a second problem
library and does not contain PDE assembly or solver implementations. The
[benchmark roadmap](../planning/chapter-6-benchmark-suite-roadmap.md)
owns selection and status, the [numerical-examples reference](../guides/chapter-6-numerical-examples.md)
owns source equations and data, and the [numerical-methods guide](../guides/chapter-6-numerical-methods.md)
owns reusable solver contracts.

The public execution path is:

```text
typed scenario
  -> Chapter 5 problem recipe
  -> ProblemSpec + runtime data bindings
  -> owned deal.II compilation session
  -> selected reduced/KKT/PDAS product
  -> solver report + compilation manifest
  -> detached experiment envelope
```

## Scenario contract

Every scenario must record these fields, even when a source example
left one of them unspecified:

| Group | Frozen contents |
| --- | --- |
| Identity | Scenario ID, recipe ID/version, source section/equation, and source revision. |
| Geometry | Dimension, domain, mesh generator or input, refinement, boundary labels, and material IDs. |
| Semantic choices | State/test/control/observation spaces, residual terms, observation region, loss, metric, constraint, and formulation product. |
| Runtime data | Forcing, target, coefficients, fixed data, transport/weight data, bounds, and provenance strings. |
| Discretization | Finite-element degrees, quadrature/evaluation policies, execution realization, and boundary orientation. |
| Algorithm | Initial control, solver family, direction, line search or fixed step, iteration limits, and stopping tolerances. |

If the source does not specify a value, the scenario must name a recovered
value or a manufactured replacement and explain why. It must not silently use
a library default when that changes the mathematical or numerical comparison.
Benchmark-specific evidence, artifact formats, and acceptance criteria belong
to the [Chapter 6 benchmark specification](../benchmarks/chapter-6.md).

## Activation matrix

The catalogue covers the current B0–B2 scenarios and the later B3–B6
application shapes. The benchmark roadmap owns their activation order and
dependencies.

| ID | Source example | Recipe family | Product/method |
| --- | --- | --- | --- |
| B0 | Common harness | Any selected recipe | Deterministic manifest/report artifact path |
| B1 | E6.5.1 distributed Laplace control | `scalar-diffusion-reaction-volume` | Reduced DTO; steepest descent versus L-BFGS |
| B2 | E6.5.2 Graetz-flow boundary control | `scalar-neumann-convection-subdomain` | Reduced DTO; full-BFGS direction with declared Armijo or fixed-step globalization |
| B3 | E6.9.1 symmetric box Laplace | Distributed scalar with cellwise box | Scalar KKT/PDAS; cellwise-discontinuous control |
| B4 | E6.9.2 asymmetric box Laplace | Distributed scalar with spatially varying cellwise bounds | Reuse B3 KKT/PDAS product |
| B5 | E6.7.1 all-at-once Laplace | Distributed scalar | Quadratic KKT and optional preconditioning |
| B6 | E6.7.2 diffusion-reaction all-at-once | Reuse B5 recipe/scenario shape | Reuse B5 implementation |

Stokes, stabilization comparisons, generic continuous-control boxes,
measure-valued state constraints, and automatic OtD derivation remain outside
this application catalogue.

## Benchmark boundary

The B0 harness, benchmark artifact schema, run path, output inventory, and
acceptance evidence are defined in the [Chapter 6 benchmark
specification](../benchmarks/chapter-6.md#common-b0-contract). Application
scenarios provide the typed scenario, recipe, runtime bindings, and selected
compiled product consumed by that boundary; they do not redefine the harness
or artifact contract.

## B1 — distributed Laplace control

### Scenario selection

Use the distributed scalar recipe with homogeneous Dirichlet data and
full-domain $L^{2}$ state tracking. The authoritative scenario selects a
cellwise `FE_DGQ(0)` volume control and its positive cellwise $L^{2}$ metric.
Development scenarios may instead select independent homogeneous-Dirichlet
continuous `FE_Q` control with its assembled $L^{2}$ metric. The scenario
records this discrete control representation and exposes the runtime target
and forcing selections; the
[source definition](../guides/chapter-6-numerical-examples.md) and the
[current benchmark freeze](../benchmarks/chapter-6.md) are maintained in their
respective documents.

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

### B1 backend execution adapter

The deal.II-specific adapter is declared in
`include/nmopt/application/dealii/chapter6_b1.hpp`. For a registered forcing
choice, the complete assembly is:

```cpp
#include "nmopt/application/dealii/chapter6_b1.hpp"

const auto scenario =
  nmopt::application::chapter6::make_b1_scenario(method);
const auto specification =
  nmopt::application::chapter6::make_b1_problem_spec(scenario);

nmopt::application::chapter6::dealii::B1SelectedDataT<2> data(
  scenario.problem.forcing, scenario.problem.desired_state);
const auto runtime =
  nmopt::application::chapter6::dealii::make_b1_runtime_data(
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

The adapter binds the selected forcing and desired-state scalar definitions,
diffusion, reaction, and the per-run regularization value; creates the owned
square-domain mesh session; compiles the assembled reduced DTO; and dispatches
either steepest descent or limited-memory BFGS. The mesh selection may be the original refined
hypercube, a uniformly subdivided triangular mesh, or a triangular mesh in
which a deterministic subset of base triangles is split at its centroid. The
last selection exposes the base subdivision count, number of splits, and
selection seed so that a candidate can match reported cell and vertex counts
without presenting the source's omitted connectivity as recovered fact.
Simplex selections currently require the registered continuous homogeneous-
Dirichlet control target.

The returned detached envelope contains the compiler manifest, typed solver
report, policy snapshot, and environment. The authoritative family selects the
source-oriented constant-half replacement; manufactured zero and the
constant-one hypothesis inferred from the Figure 6.2 extrema remain checked
development definitions. Parameter families may also supply any finite
constant or a validated scalar expression through the same backend-neutral
function record. Arbitrary caller-owned deal.II functions can still use
`B1RuntimeDataT<dim>` directly.

`parameters/chapter-6/b1/development/continuous-control.prm` is the checked
candidate motivated by the source's continuous linear control-space statement
and equal reported state/control/adjoint counts. On the default quadrilateral
mesh it realizes the conforming `Q1` analogue, not the source's undisclosed
triangular `P1` mesh, so it is comparison evidence rather than source parity.
The simplex generators make the corresponding `P1` experiments executable;
`continuous-control-structured-simplex.prm` records the nearly count-matched
regular-grid hypothesis and `continuous-control-count-matched-simplex.prm`
records the exact-count deterministic centroid-split hypothesis. Neither file
claims to recover the omitted source connectivity.
The companion `continuous-control-constant-one.prm` changes only the forcing
to the Figure 6.2 constant-one hypothesis while retaining the same continuous
control and recovered Figure 6.3 solver policies.

The benchmark-specific B1 freeze and acceptance evidence are defined in the
[Chapter 6 benchmark specification](../benchmarks/chapter-6.md).

## B2 — Graetz-flow boundary control

### Scenario selection

Use the `scalar-neumann-convection-subdomain` recipe, seeded by
`make_neumann_convection_subdomain_tracking_problem(observed_material_id,
with_facewise_box=false, control_discretisation)`. The public B2 helper adds
the nonzero fixed-data port and lifting policy before compilation, so an
adapter only needs to assemble the typed scenario and bind the declared
runtime functions:

```cpp
const auto scenario =
  nmopt::application::chapter6::make_b2_scenario(
    nmopt::application::chapter6::B2ObservationRegion::wings,
    "constant");
const auto specification =
  nmopt::application::chapter6::make_b2_problem_spec(scenario);
```

The resulting state has scalar conservative transport, fixed temperature on
the declared Dirichlet boundary, controlled Neumann flux on a marked
boundary, insulated outflow, and volume observation on a downstream material
region. The frozen default uses facewise constants and `l2_facewise`.
Development scenarios may instead select a continuous degree-one nodal trace
and `l2_neumann_trace`; that representation rejects the facewise box. The
helper declares `fixed_dirichlet_data` as a `fixed_dirichlet_lifting` function
and connects it through `fixed_dirichlet_reconstruction`; passing that
function to a graph without the declared port remains an error.

The scenario solver record also selects `ReducedGlobalization::armijo` or
`ReducedGlobalization::fixed_step`. Armijo is the frozen default. Fixed-step
uses `initial_step_length` as the declared positive step, evaluates one trial,
and does not add a sufficient-decrease condition or alter the selected
stopping criterion. B1 remains Armijo-only until its separate execution
adapter registers another policy.

The B2 compile record requires `VolumeObservationOptions`. Its frozen value
uses order `3` and analytic desired-state evaluation at observation quadrature
points. Development scenarios may instead select a different positive order
and `state-fe-interpolation`, which interpolates the desired-state function
into the scalar state finite-element space before observation assembly. B1
rejects this B2-only option rather than silently ignoring it.

For B2, equation (6.65) fixes the source boundary form as the ordinary-normal
condition
$\partial_{n} y-(b\mathbin\cdot n)y=u$ on the control boundary and zero on the
outflow. `make_b2_scenario(...)` also accepts
`TransportBoundaryForm::total_conormal` as an explicit diagnostic choice. That
choice realizes the outward diffusion-minus-transport conormal
$\mu\partial_{n} y-(b\mathbin\cdot n)y$ on the same boundary partition; it does
not change the frozen benchmark and must not be used as source-replication
evidence. The compiler derives its deal.II realization from the typed boundary
selection and records the selected form in the manifest. The [source
definition](../guides/chapter-6-numerical-examples.md)
and [frozen benchmark policy](../benchmarks/chapter-6.md) are maintained in
their respective documents.

The public catalog exposes the four B2 records through the Cartesian product
of the independent observation-region and target-profile selections in
`make_catalog()`. Production parameter files carry the selected observation
ID and its scalar indicator definition; the application owns the stable IDs
and recipe/runtime construction, while the benchmark specification owns the
observation/target values and frozen run matrix.

### B2 backend execution adapter

The deal.II-specific adapter is declared in
`include/nmopt/application/dealii/chapter6_b2.hpp`. It supplies the two
dimensional rectangle mesh, selected as either the frozen framework-native
quadrilateral realization, a structured triangular realization with per-axis
subdivisions, or a deterministic centroid-split variant of that triangular
mesh. All paths use the same material-ID observation realization and boundary
labels. Each centroid split replaces one selected base triangle by three,
adding one vertex and two cells; the seed makes this topology candidate
repeatable but does not identify the source's omitted connectivity. The
adapter also supplies zero forcing, the selected target function, the selected
fixed Dirichlet function, and the selected conservative transport field. The
complete manufactured-data path is:

```cpp
#include "nmopt/application/dealii/chapter6_b2.hpp"

const auto scenario =
  nmopt::application::chapter6::make_b2_scenario(
    nmopt::application::chapter6::B2ObservationRegion::wings);
const auto specification =
  nmopt::application::chapter6::make_b2_problem_spec(scenario);

nmopt::application::chapter6::dealii::B2ManufacturedDataT<2> data{
  nmopt::application::chapter6::b2_manufactured_wings_observation(),
  scenario.problem.fixed_dirichlet_data,
  scenario.problem.forcing,
  nmopt::application::chapter6::b2_target_definition(
    scenario.problem.target_catalog),
  scenario.problem.conservative_transport};
const auto runtime =
  nmopt::application::chapter6::dealii::make_b2_manufactured_runtime_data(
    scenario, data);
const auto session =
  nmopt::application::chapter6::dealii::make_b2_compilation_session<2>(
    scenario);
nmopt::application::chapter6::dealii::B2ReducedExecutionAdapterT<2> execute{
  runtime, session, environment};

using Runner =
  nmopt::application::benchmark::HeadlessBenchmarkRunnerT<
    decltype(scenario)>;
const auto result = Runner(scenario).run(
  [](const auto &parameters) {
    return nmopt::application::chapter6::make_b2_problem_spec(parameters);
  },
  execute);
```

The adapter binds every declared port, including `fixed_dirichlet_data` and
`conservative_transport`, compiles the assembled reduced DTO, and dispatches
full BFGS with the selected Armijo or fixed-step globalization from zero
control in the selected layout. It also maps `VolumeObservationOptions` to the
compiler policy and checks the effective target realization and quadrature
order against the compilation manifest. The declared discretization and
effective solver policy are retained in the artifact and solver snapshot.
The backend-neutral `B2TargetParameters` record carries two
`ScalarFunctionDefinition` values for the constant and parabolic target
profiles. Its defaults are the source value `2` and expression
`4.0*x1*(1.0-x1)`; parameter-file development candidates may change either
definition while retaining explicit provenance. The deal.II adapter realizes
the selected definition as a constant or parsed coordinate expression. The
artifact records the selected definition, kind, value, and expression so
target-transcription experiments remain distinguishable from source data.
`B2RuntimeDataT<dim>` is
the extension point for recovered forcing, targets, or transport fields; the
caller owns all referenced Function objects for the duration of compilation
and solving. Native boundary output stores facewise constants as cell data and
continuous trace controls as point data on connected line cells. The
postprocessor renders either representation; point-data values are averaged at
each line segment for the line color while retaining the nodal extrema for the
color scale.

### Runtime bindings

The runtime binding must include forcing, desired state, diffusion, reaction,
regularization, fixed Dirichlet data, conservative transport, and provenance
for each. The graph and mesh must agree on fixed, controlled, outflow, and
observation labels. A transport field must not be passed through the general
scalar bundle with a different semantic role.

The benchmark-specific frozen B2 inputs, stabilization boundary, and required
evidence are defined in the [Chapter 6 benchmark
specification](../benchmarks/chapter-6.md).

## Later scenario contracts

B3 and B4 must use cellwise-discontinuous control for the first framework
native PDAS runs. The coefficientwise box, positive-diagonal `l2_cellwise`
metric, KKT product, PDAS complementarity, and active-set evidence must all
refer to the same control layout. A continuous `Q1` box reproduction is a
separate feature and is not implied by B3.

B5 and B6 are all-at-once products. Their source scaling and KKT details are
recorded in the [numerical-examples reference](../guides/chapter-6-numerical-examples.md);
their scenario construction belongs here, and their frozen run/evidence
contracts belong in the benchmark specification. Do not activate a
preconditioner merely because the source used one; measure whether the
selected scalar run needs it first.

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
7. Which benchmark-specific evidence fields and artifacts must this scenario
   expose for acceptance and performance reporting?
