# deal.II compiler reference

This reference covers the public v1 deal.II compiler: validation, runtime data
binding, mesh lifetime, discretization policy, product selection, diagnostics,
and the compiled products returned to application code.

It assumes that a valid `semantic::v1::ProblemSpec` already describes the
mathematical problem. See [Problem authoring](problem-authoring.md) for that
programming interface. For the conceptual compiler pipeline, see
[Compilation and lowering](../manual/concepts/11-compilation-and-lowering.md)
and [Authoring and using compiled problems](../manual/concepts/12-authoring-and-using-compiled-problems.md).

The compiler combines semantic meaning with a concrete numerical realization.
It does not choose an optimization algorithm and it does not silently repair
unsupported semantic requests.

## Public headers and namespaces

Typical compiled-path code includes:

```cpp
#include "nmopt/compiler/v1/dealii_compiler.hpp"
#include "nmopt/compiler/v1/dealii_types.hpp"
#include "nmopt/semantic/v1/problem_spec.hpp"

#include <deal.II/grid/grid_generator.h>
```

The examples use:

```cpp
namespace sem = nmopt::semantic::v1;
namespace comp = nmopt::compiler::v1;

using Backend = nmopt::dealii_backend::SerialBackend;
```

## Compiler workflow

The public path is:

```text
                       ProblemSpec
                           │
                           ▼
                  DealiiCompiler::validate
                           │
             diagnostics ──┴── accepted request
                                  │
          ┌───────────────────────┼────────────────────────┐
          │                       │                        │
          ▼                       ▼                        ▼
DealiiCompilationSession   DealiiDataBindings   DealiiDiscretisationPolicy
          │                       │                        │
          └───────────────────────┼────────────────────────┘
                                  ▼
                         DealiiCompiler::compile
                                  │
                                  ▼
                         CompilationResultT
                                  │
            ┌──────────┬──────────┼──────────────────┐
            ▼          ▼          ▼                  ▼
        problem  kkt_problem  pdas_problem  supplied_otd_problem
```

Only one product member is expected for a successful compilation.

The semantic formulation and the requested `CompilationProduct` jointly decide
which member is returned. In particular, supplied OTD is a semantic
formulation, not a fourth `CompilationProduct` enum value.

## Working path: compile one problem end to end

This section continues the scalar distributed-control example from
[Problem authoring](problem-authoring.md). The semantic graph is already
decided; the compiler now needs the concrete numerical information that was
intentionally absent from `ProblemSpec`.

The four inputs have different jobs:

```text
ProblemSpec                 what mathematical problem?
    │
    ├──── runtime bindings  what concrete data values?
    ├──── mesh/session      on what discrete domain and with whose lifetime?
    ├──── policy            with what registered numerical realization?
    └──── product           which solver-facing product should be built?
```

Keeping those concerns separate is the central rule of the compiled path.

### Step 1 – recover the semantic request

```cpp
#include "nmopt/application/chapter5.hpp"
#include "nmopt/compiler/v1/dealii_compiler.hpp"

#include <deal.II/base/function_lib.h>
#include <deal.II/grid/grid_generator.h>

#include <iostream>
#include <memory>
#include <utility>

namespace sem = nmopt::semantic::v1;
namespace comp = nmopt::compiler::v1;
namespace ch5 = nmopt::application::chapter5;

using Backend = nmopt::dealii_backend::SerialBackend;

constexpr int dim = 2;

void
compile_example()
{
  // This is the semantic family chosen in the authoring layer.
  //
  // The recipe only constructs ProblemSpec. It does not know the mesh,
  // forcing Function, finite-element degree, or outer optimizer.
  const auto recipe = ch5::make_scalar_distributed_recipe();

  ch5::ScalarDistributedControlParameters problem_parameters;
  problem_parameters.discretisation =
    ch5::DistributedControlDiscretisation::cellwise_constant;
  problem_parameters.with_cellwise_box = false;

  const sem::ProblemSpec spec = recipe(problem_parameters);
```

**At this point:** `spec` tells the compiler what the PDE/control problem
means, but none of its semantic spaces has a concrete coefficient dimension.

### Step 2 – bind concrete runtime data

The graph declared ports such as `"forcing"`, `"desired_state"`,
`"diffusion"`, and `"regularisation_weight"`. Now supply the actual values for
this run.

```cpp
  // These deal.II Functions remain application objects. DealiiDataBindings
  // stores references, so they must outlive every compiled callback that may
  // evaluate them.
  dealii::Functions::ZeroFunction<dim> forcing;
  dealii::Functions::ConstantFunction<dim> desired_state(1.0);

  // DealiiDataBindings is the bridge from semantic data ports to concrete
  // backend values.
  //
  // The order here is:
  //   forcing Function
  //   desired-state Function
  //   optional constant diffusion
  //   constant reaction
  //   regularisation weight
  //   provenance for the common Function ports
  //
  // The semantic graph decides whether each port is meaningful for the
  // selected target. These numeric values do not alter ProblemSpec itself.
  comp::DealiiDataBindings<dim> data{
    forcing,
    desired_state,
    1.0,    // diffusion coefficient used by this scalar target
    0.0,    // reaction coefficient
    1e-2,   // control regularisation weight
    {
      "generated zero forcing",
      "constant target 1",
      ""    // no fixed Dirichlet datum in this graph
    }
  };
```

The distinction here is important:

```text
ProblemSpec::DataSpec
    "desired_state"
          │
          │ declares meaning / expected role
          ▼
DealiiDataBindings
    desired_state Function
          │
          │ supplies values
          ▼
compiler realization
```

Changing the values represented by `desired_state` does not require rebuilding
the semantic graph. Changing the *role or observation structure* of that datum
does.

### Step 3 – create the mesh and decide who owns it

```cpp
  auto triangulation =
    std::make_unique<dealii::Triangulation<dim>>();

  dealii::GridGenerator::hyper_cube(
    *triangulation,
    0.0,
    1.0);

  triangulation->refine_global(4);

  // Move the triangulation into an owned compilation session.
  //
  // This is the normal long-lived compiled path: the session becomes the
  // lifetime owner of the mesh and can be retained by compiled callbacks.
  // The caller no longer has a mutable Triangulation handle after the move.
  auto session =
    std::make_shared<comp::DealiiCompilationSession<dim>>(
      std::move(triangulation),
      "unit square; GridGenerator::hyper_cube; global refinement 4");
```

Use a meaningful provenance string. It is retained separately from the mesh's
structural identity; two different meshes can come from the same named source.

**At this point:** the application has a concrete mesh and concrete runtime
data, but no finite-element spaces, operators, metric, or state/adjoint solver
services have been constructed.

### Step 4 – choose a registered numerical realization

```cpp
  comp::DealiiDiscretisationPolicy policy;

  // Degree one is already the default, but setting it here makes the example's
  // FE request explicit. The semantic graph said "H1 state"; the compiler
  // policy decides the registered discrete degree used to realize it.
  policy.state_degree = 1;

  // The enum also contains matrix_free, but the current v1 compiler accepts
  // only assembled execution. Asking for matrix_free produces a lowerability
  // diagnostic instead of silently reverting to assembled.
  policy.execution =
    comp::DealiiDiscretisationPolicy::Execution::assembled;

  // reduced_dto asks for the ordinary compiled reduced service for this
  // reduced semantic formulation.
  const auto product =
    comp::CompilationProduct::reduced_dto;
```

The policy is not where you choose steepest descent, BFGS, or a trust-region
method. It controls compiler/backend realization.

### Step 5 – validate capability before construction

```cpp
  comp::DealiiCompiler compiler;

  const sem::ValidationReport validation =
    compiler.validate(spec, policy, product);

  if (!validation.valid()) {
    for (const auto &d : validation.diagnostics()) {
      std::cerr
        << "component: " << d.component_id << '\n'
        << "capability: " << d.capability << '\n'
        << "remedy: " << d.remedy << "\n\n";
    }

    // `validate()` has not needed the mesh or Functions above; in a larger
    // application this check can be performed earlier, before expensive
    // backend construction.
    return;
  }
```

This answers a different question from `SemanticValidator`:

```text
SemanticValidator
    “Is this graph coherent?”
             │
             ▼
DealiiCompiler::validate
    “Can this compiler realize this coherent graph
     with this policy and requested product?”
```

For example, selecting `matrix_free` leaves the graph semantically valid but
fails compiler lowerability.

### Step 6 – compile the concrete problem

```cpp
  const auto result = compiler.compile(
    spec,
    session,
    data,
    policy,
    std::nullopt,  // no cellwise bound values
    std::nullopt,  // no facewise bound values
    product);

  if (!result.succeeded()) {
    // Compilation also validates concrete inputs: mesh, bound Functions,
    // scalar values, optional bounds, and other runtime realization data.
    //
    // A request may therefore pass compiler.validate() and still fail here
    // because the concrete bindings or mesh are invalid.
    for (const auto &d : result.diagnostics.diagnostics()) {
      std::cerr
        << "component: " << d.component_id << '\n'
        << "capability: " << d.capability << '\n'
        << "remedy: " << d.remedy << "\n\n";
    }
    return;
  }
```

A successful compile has now created the finite-element realization, bound the
runtime data, constructed the executable model, selected the metric, and built
the state/adjoint services appropriate to this target.

### Step 7 – take the solver-facing product

For this semantic formulation plus `CompilationProduct::reduced_dto`,
`result.problem` is the selected product.

```cpp
  const auto &compiled = *result.problem;

  // make_reduced_dto() packages the executable model and the compiler's state
  // and adjoint solve services behind the common reduced formulation
  // interface. Application code should normally use this instead of manually
  // rebuilding the partition and callbacks from compiler internals.
  auto reduced = compiled.make_reduced_dto();

  // The metric is a separate solver input because a reduced derivative is a
  // covector and optimization methods need the selected Riesz map to obtain
  // a primal search gradient/direction.
  const auto &metric = compiled.metric();

  // `reduced` and `metric` are now ready for the optimization layer.
}
```

**At this point:** the semantic graph has become an executable numerical
service. From here onward, concrete dimensions and realized FE choices belong
to the compiled product and its manifest, not to `ProblemSpec`.

The important object progression is:

```text
ProblemSpec
   │ backend-independent request
   ▼
CompilationResultT
   │ validated concrete construction
   ▼
CompiledProblemT
   │ compiler-owned executable bundle
   ├── executable model
   ├── metric
   ├── optional constraint / Hessian / box data
   ├── state-adjoint solves
   └── manifest
           │
           ▼
     make_reduced_dto()
           │
           ▼
       ReducedDTOT
```

## Variations on the working example

The shortest way to learn the compiler API is to change one layer at a time
while keeping the rest of the preceding example fixed.

### Variation: add a cellwise box

First change the semantic family:

```cpp
problem_parameters.with_cellwise_box = true;
const sem::ProblemSpec bounded_spec =
  recipe(problem_parameters);
```

That makes the lower/upper bound ports and `ConstraintSpec` part of the
semantic request.

Then supply the runtime values separately:

```cpp
comp::CellwiseBoxDataBindings bounds{
  -0.5,  // same lower value for every realized cellwise control coefficient
   0.5   // same upper value for every realized cellwise control coefficient
};

const auto bounded_result = compiler.compile(
  bounded_spec,
  session,
  data,
  policy,
  bounds,        // concrete values for the semantic cellwise box
  std::nullopt,  // no facewise box
  comp::CompilationProduct::reduced_dto);
```

After success:

```cpp
const auto &bounded = *bounded_result.problem;

const auto *constraint = bounded.constraint();
if (constraint == nullptr)
  throw std::runtime_error(
    "bounded compilation did not produce its projection constraint");

auto bounded_reduced = bounded.make_reduced_dto();
```

The important split is:

```text
ProblemSpec              declares “this problem has a cellwise box”
CellwiseBoxDataBindings  supplies “these are the lower/upper values”
CompiledProblemT         returns the compatible realized constraint
```

Do not pass `CellwiseBoxDataBindings` to an unbounded graph and expect the
compiler to invent a semantic constraint.

### Variation: bind fixed nonzero Dirichlet data

Use a graph that declares fixed-data reconstruction:

```cpp
const sem::ProblemSpec fixed_spec =
  sem::make_fixed_dirichlet_scalar_diffusion_reaction_problem();
```

Create the additional application-owned function:

```cpp
dealii::Functions::ConstantFunction<dim> fixed_boundary_value(2.0);
```

Then fill the optional fixed-data binding in `DealiiDataBindings`:

```cpp
comp::DealiiDataBindings<dim> fixed_data{
  forcing,
  desired_state,
  1.0,
  0.0,
  1e-2,
  {
    "generated zero forcing",
    "constant target 1",
    "constant fixed boundary value 2"
  },
  std::cref(fixed_boundary_value)
};
```

The semantic transformation says **why** fixed data is needed and where it
acts. The binding supplies the concrete `Function` that realizes that port.

### Variation: request a quadratic KKT product

Keep the compatible canonical unconstrained scalar graph and change only the
requested product:

```cpp
const auto kkt_result = compiler.compile(
  spec,
  session,
  data,
  policy,
  std::nullopt,
  std::nullopt,
  comp::CompilationProduct::quadratic_kkt);

if (!kkt_result.succeeded())
  return;

const auto &compiled_kkt = *kkt_result.kkt_problem;
const auto &kkt_product = compiled_kkt.product();
```

Do not call `result.problem->make_reduced_dto()` in this branch:
`kkt_problem` is a different solver-facing product.

If the same request is made for a reduced graph outside the registered
quadratic KKT capability, compilation returns a `formulation_capability`
diagnostic.

### Variation: semantic supplied OTD

For supplied OTD, change the semantic graph rather than inventing a new
`CompilationProduct` value:

```cpp
const sem::ProblemSpec otd_spec =
  sem::make_scalar_diffusion_reaction_supplied_otd_problem();

const auto otd_result = compiler.compile(
  otd_spec,
  session,
  data,
  policy,
  std::nullopt,
  std::nullopt,
  comp::CompilationProduct::reduced_dto);

if (!otd_result.succeeded())
  return;

const auto &compiled_otd =
  *otd_result.supplied_otd_problem;

const auto &system =
  compiled_otd.system();
```

Here `reduced_dto` means “do not additionally request KKT or PDAS.” The
semantic formulation itself selects the supplied-OTD execution product.

## Validate before expensive construction

`DealiiCompiler::validate()` takes:

```cpp
sem::ValidationReport validate(
  const sem::ProblemSpec &specification,
  const comp::DealiiDiscretisationPolicy &policy,
  comp::CompilationProduct product =
    comp::CompilationProduct::reduced_dto) const;
```

It performs semantic resolution first. If the semantic graph is invalid, the
compiler returns those diagnostics without attempting lowerability analysis.

For a valid semantic graph it then checks:

- whether the complete graph matches a registered lowering path;
- the selected finite-element/discretization policy;
- formulation capability;
- supplied-OTD declaration capability;
- requested KKT or PDAS product capability;
- typed realization requirements.

Validation does not require a triangulation or concrete runtime `Function`
objects. Use it when selecting among problem families or products before
building expensive backend state.

Compilation repeats the relevant validation and additionally checks the
concrete mesh and runtime bindings.

## Reading diagnostics as instructions

A diagnostic should tell you **which layer to edit**. A few representative
cases make the distinction clearer.

### Example: unsupported execution policy

```cpp
policy.execution =
  comp::DealiiDiscretisationPolicy::Execution::matrix_free;

const auto report =
  compiler.validate(spec, policy, product);
```

The current compiler reports capability `assembled_execution` in the
`lowerability` category.

The semantic graph is not the problem. Change the compiler policy:

```cpp
policy.execution =
  comp::DealiiDiscretisationPolicy::Execution::assembled;
```

### Example: broken semantic reference

Suppose a hand-authored equation names a missing test space:

```cpp
spec.equations.front().test_space_id =
  "state_test_space_typo";
```

This is a `structural` problem. Do not change deal.II policy or bindings; fix
the stable ID relationship in the graph.

### Example: box declared but values omitted

If `spec.formulation.constraint_id` names a registered cellwise box but the
compile call supplies:

```cpp
std::nullopt  // cellwise bounds
```

the failure belongs to concrete lowerability. The graph correctly says that
bounds exist; the application failed to supply their numerical realization.

### Example: request the wrong product

Requesting:

```cpp
comp::CompilationProduct::quadratic_kkt
```

for a semantically valid graph that is not one of the registered KKT targets
produces a `formulation_capability` diagnostic.

Do not reinterpret the graph as quadratic or switch products silently. Either
choose `reduced_dto` for that application or add the missing KKT capability
deliberately.

These cases summarize the editing rule:

```text
structural             → edit ProblemSpec
analytical_policy      → edit semantic requirement declarations
lowerability           → edit mesh, bindings, compiler policy, or supported graph
formulation_capability → edit formulation/product choice or implement the capability
```

## Diagnostic categories

All validation layers use `semantic::v1::ValidationReport`.

The category tells you which layer owns the correction:

- `structural` – fix the semantic graph;
- `analytical_policy` – declare the required assumption or realization;
- `lowerability` – change the concrete compiler request, bindings, mesh, or
  discretization;
- `formulation_capability` – select a supported formulation/product
  combination.

Each `Diagnostic` also records:

```cpp
component_id
capability
remedy
```

Treat `remedy` as the supported correction, not as permission to substitute a
nearby component silently.

Expected unsupported-input failures are returned as diagnostics. Direct
low-level contract misuse and violated invariants may throw
`nmopt::contract::ContractError`.

## The compiler input boundary

Before looking at individual binding types, it is useful to classify values by
where they belong.

If changing a value changes the **meaning or composition** of the mathematical
problem, it belongs in `ProblemSpec` or recipe parameters.

If it supplies the **concrete value of a declared semantic datum**, it belongs
in `DealiiDataBindings` or one of its specialized binding records.

If it changes the **discrete realization or backend solve policy**, it belongs
in `DealiiDiscretisationPolicy`.

If it changes the **outer optimization algorithm**, it does not belong in the
compiler at all.

For example:

```text
“use a Neumann control”              → semantic graph
“velocity b(x) = (...)”              → runtime binding
“state FE degree = 2”                → compiler policy
“use L-BFGS with memory 10”          → optimization layer
```

This classification is more useful than memorizing constructor argument lists.

## Runtime data bindings

`ProblemSpec::data` declares semantic ports. `DealiiDataBindings<dim>` supplies
the concrete values used by the current deal.II compiler.

The base constructor is:

```cpp
comp::DealiiDataBindings<dim> data{
  forcing,
  desired_state,
  diffusion,
  reaction,
  regularisation_weight,
  provenance,
  fixed_dirichlet_data,
  general_scalar_data,
  weighted_trace_data,
  conservative_transport_data,
  natural_boundary_source_data
};
```

Only the first six arguments are required by the C++ constructor. The semantic
graph determines which ports are actually required by a particular
compilation.

### Base bindings

The common ports are:

| Binding | Concrete type | Typical semantic role |
| --- | --- | --- |
| `forcing` | `const dealii::Function<dim> &` | volume forcing |
| `desired_state` | `const dealii::Function<dim> &` | target evaluated by the selected observation |
| `diffusion` | `std::optional<double>` | constant diffusion for coefficient-bound scalar targets |
| `reaction` | `double` | constant reaction coefficient |
| `regularisation_weight` | `double` | scalar objective regularization weight |
| `fixed_dirichlet_data` | optional `reference_wrapper<const Function<dim>>` | fixed nonzero essential data |

The concrete functions are referenced, not copied. Keep every bound object
alive for as long as the compiled callbacks may evaluate it.

`diffusion` is optional because not every target receives its diffusion through
the constant scalar port. Coefficient identification, for example, takes the
physical diffusion from its decision variable.

### Provenance

`DealiiBindingProvenance` records labels for the common function inputs:

```cpp
comp::DealiiBindingProvenance provenance{
  "forcing from input/source.dat",
  "target generated analytically",
  "Dirichlet data from boundary.json"
};
```

A `dealii::Function` cannot describe where its values came from. Provenance is
therefore a separate compiler input and is retained in the compilation
manifest.

Do not use provenance strings as capability selectors. They are descriptive
records.

### General scalar coefficients

The registered general scalar elliptic/Robin family binds spatially varying
coefficient functions through `DealiiGeneralScalarDataBindings<dim>`:

```cpp
comp::DealiiGeneralScalarDataBindings<dim> general{
  diffusion_tensor,
  conservative_transport,
  advective_transport,
  reaction_function,
  robin_coefficient,
  robin_source,
  {
    "diffusion tensor source",
    "conservative transport source",
    "advective transport source",
    "reaction source",
    "Robin coefficient source",
    "Robin source"
  }
};

comp::DealiiDataBindings<dim> data{
  forcing,
  desired_state,
  std::nullopt,
  0.0,
  beta,
  common_provenance,
  std::nullopt,
  general
};
```

The types distinguish tensor, vector, and scalar fields at the compiler
boundary. The compiler does not reinterpret a scalar function as a tensor
coefficient or silently take component zero from a multi-component function.

### Weighted trace data

A weighted boundary observation uses:

```cpp
comp::DealiiWeightedTraceDataBindings<dim> weighted_trace{
  boundary_weight,
  "boundary-weight source"
};
```

This supplies the immutable weight named by the corresponding semantic
observation.

### Conservative transport data

The narrower conservative-transport Neumann family can bind only the vector
coefficient it needs:

```cpp
comp::DealiiConservativeTransportDataBindings<dim> transport{
  velocity,
  {"transport source"}
};
```

This is intentionally different from supplying the full general-scalar
coefficient bundle.

### Natural boundary source

An additive immutable natural-boundary source is bound as:

```cpp
comp::DealiiNaturalBoundarySourceBinding<dim> natural_source{
  source_function,
  "natural boundary source provenance"
};
```

The semantic residual term decides where this datum enters the state load.

## Mesh lifetime

The compiler has two public `compile()` overload families.

### Owned session

Prefer:

```cpp
auto session =
  std::make_shared<comp::DealiiCompilationSession<dim>>(
    std::move(triangulation),
    mesh_provenance);

auto result =
  compiler.compile(spec, session, data, policy, ...);
```

The constructor requires:

- a non-null `std::unique_ptr<dealii::Triangulation<dim>>`;
- nonempty mesh provenance.

Ownership is transferred into the session. The compiler may mutate the
triangulation while constructing its finite-element realization; callers
cannot mutate it through the public session API afterward.

Compiled products retain the session as an opaque lifetime owner where needed.
The manifest records the mesh lifetime as `owned_session`.

This is the safest default when a compiled `ReducedDTOT`, KKT product, or
supplied OTD system may be detached from the local compilation variables.

### Borrowed triangulation

The alternate overload accepts:

```cpp
dealii::Triangulation<dim> &triangulation
```

directly.

```cpp
auto result = compiler.compile(
  spec,
  triangulation,
  data,
  policy,
  bounds,
  facewise_bounds,
  product);
```

This is an explicitly borrowed immutable-lifetime path. The caller must keep
the triangulation and all callback dependencies alive and must not invalidate
them while compiled services can still execute.

Use it only when the surrounding application already has a clear lifetime
owner that is stronger than the compiler session.

## Bound data

Bounds are runtime realization data and therefore are passed separately from
the semantic `ConstraintSpec`.

### Cellwise bounds

For a cellwise box:

```cpp
comp::CellwiseBoxDataBindings bounds{
  -1.0,   // lower
   1.0    // upper
};

auto result = compiler.compile(
  spec,
  session,
  data,
  policy,
  bounds,
  std::nullopt,
  comp::CompilationProduct::reduced_dto);
```

Each side is:

```cpp
std::variant<double, dealii::Vector<double>>
```

A scalar means one constant bound for every coefficient. A vector provides the
coefficientwise values.

Both bounds must be finite and ordered. Vector bounds must match the realized
control layout.

### Facewise bounds

Boundary facewise controls use the distinct type:

```cpp
comp::FacewiseBoxDataBindings bounds{
  lower,
  upper
};
```

Pass it in the second optional bound slot:

```cpp
auto result = compiler.compile(
  spec,
  session,
  data,
  policy,
  std::nullopt,
  bounds,
  comp::CompilationProduct::reduced_dto);
```

Facewise and cellwise bounds are not interchangeable. They correspond to
different discrete control orderings.

### One compiled box realization

For a compiled cellwise box, the compiler creates one shared
`CompiledCellwiseBoxDataT` and derives the projection constraint and PDAS
complementarity data from it. The objects retain a common token, metric
realization, layout signature, bounds digest, and provenance.

Application code should consume `compiled.constraint()`,
`compiled.box_data()`, or the corresponding PDAS product rather than
reconstructing a second nominally identical box by hand.

## Discretization policy

`DealiiDiscretisationPolicy` is the main compiler option record:

```cpp
comp::DealiiDiscretisationPolicy policy;

policy.state_degree = 2;
policy.execution =
  comp::DealiiDiscretisationPolicy::Execution::assembled;

policy.state_solve.maximum_iterations = 500;
policy.state_solve.relative_tolerance = 1e-10;
policy.state_solve.absolute_tolerance = 1e-12;
```

The fields with general application are:

| Field | Default | Behavior |
| --- | --- | --- |
| `state_degree` | `1` | conforming scalar Lagrange degree; zero is rejected |
| `execution` | `assembled` | current v1 lowerers accept only assembled execution |
| `state_solve` | CG policy, automatic iteration cap | state solve policy for iterative SPD targets |
| `adjoint_solve` | CG policy, automatic iteration cap | adjoint solve policy for iterative SPD targets |
| `control_metric_solve` | 1000, `1e-12`, `1e-14` | iterative inverse for applicable mass-like metrics |
| `pdas` | typed PDAS defaults | outer active-set policy for a compiled PDAS product |
| `pdas_kkt_solver` | KKT policy, maximum 200 | inner KKT solver policy for compiled PDAS |

`SPDLinearSolvePolicy::maximum_iterations == 0` selects the current
dimension-dependent rule

```math
\max(100, 10 n),
```

where $n$ is the realized operator dimension. Relative and absolute
tolerances must be positive and finite.

`control_metric_solve.maximum_iterations` must be positive.

### Current execution support

Although the policy enum contains `matrix_free`, v1 currently rejects every
non-assembled request with the `assembled_execution` lowerability diagnostic.

Select:

```cpp
policy.execution =
  comp::DealiiDiscretisationPolicy::Execution::assembled;
```

unless the compiler implementation is extended with a registered matrix-free
lowerer.

### Target-specific policy

Some fields are meaningful only for particular registered targets.

`volume_observation` optionally selects:

```cpp
comp::VolumeObservationDiscretisationPolicy{
  quadrature_order,
  comp::VolumeObservationTargetRealisation::analytic_quadrature
};
```

or:

```cpp
comp::VolumeObservationTargetRealisation::state_fe_interpolation
```

The current conservative-transport/subdomain target requires a positive
quadrature order when this option is supplied. If omitted, the effective
policy uses `state_degree + 2` and analytic quadrature.

`transport_boundary_realisation` is likewise a specialized boundary policy.
Do not set specialized options speculatively; validation checks whether the
selected target accepts them.

## Compilation products

`CompilationProduct` currently contains:

```cpp
comp::CompilationProduct::reduced_dto
comp::CompilationProduct::quadratic_kkt
comp::CompilationProduct::pdas
```

It is a request for the downstream numerical product, not a semantic
formulation label.

The result shape is:

```cpp
template <typename Backend>
struct CompilationResultT
{
  sem::ValidationReport diagnostics;

  std::shared_ptr<const CompiledProblemT<Backend>>
    problem;

  std::shared_ptr<const CompiledQuadraticKKTProblemT<Backend>>
    kkt_problem;

  std::shared_ptr<const CompiledPDASProblemT<Backend>>
    pdas_problem;

  std::shared_ptr<const CompiledSuppliedOTDProblemT<Backend>>
    supplied_otd_problem;

  bool succeeded() const;
};
```

The following lookup is the important one:

| Semantic formulation | Requested product | Successful result member |
| --- | --- | --- |
| reduced DTO | `reduced_dto` | `problem` |
| reduced DTO | `quadratic_kkt` | `kkt_problem`, only for a registered canonical KKT target |
| reduced DTO with supported cellwise box | `pdas` | `pdas_problem`, only for the registered PDAS target |
| supplied OTD without a box | `reduced_dto` | `supplied_otd_problem` |
| supplied OTD without a box | `quadratic_kkt` | `kkt_problem` when the canonical OTD-to-KKT bridge is valid |
| supplied OTD with supported cellwise box | `pdas` | `pdas_problem` when the canonical OTD-to-KKT/PDAS path is valid |

There is no `CompilationProduct::supplied_otd`. The semantic formulation
selects supplied OTD; the product argument optionally asks the compiler to
bridge that system into the registered KKT or PDAS product.

Unsupported combinations produce `formulation_capability` diagnostics rather
than falling back to another product.

## Consuming a reduced compiled problem

For the ordinary reduced path:

```cpp
const auto &compiled = *result.problem;
```

`CompiledProblemT<Backend>` is the executable bundle returned by the compiler.

Its principal operations are:

```cpp
compiled.executable_model()
compiled.metric()
compiled.constraint()
compiled.box_data()
compiled.reduced_hessian()
compiled.state_adjoint_solvers()
compiled.make_reduced_dto()
compiled.compiled_application_view()
compiled.manifest()
```

The first two are always present on a valid `CompiledProblemT`.

The optional capabilities have precise meanings:

- `constraint()` is non-null when a projection constraint was compiled;
- `box_data()` is non-null when the compiler constructed shared cellwise box
  data;
- `reduced_hessian()` is non-null only when the selected numerical target
  exposes the explicit reduced-Hessian action;
- `compiled_application_view()` is non-null only for compiler paths that expose
  the small application-output seam.

Check the pointer. Do not infer an optional capability from the semantic metric
name, objective type, or problem label.

### Create the reduced formulation

The normal handoff to reduced optimization is:

```cpp
auto reduced = compiled.make_reduced_dto();

nmopt::solvers::ReducedGradientSolverT<Backend> solver(
  reduced,
  compiled.metric(),
  parameters);

auto optimization_result = solver.solve(initial_control);
```

When a compatible projection constraint exists:

```cpp
const auto *constraint = compiled.constraint();

if (constraint == nullptr)
  throw std::runtime_error("this solve requires the compiled box");

nmopt::solvers::ReducedGradientSolverT<Backend> solver(
  reduced,
  compiled.metric(),
  *constraint,
  parameters);
```

For the meaning and relationships of the reduced optimization methods, see
[Reduced optimization methods](../manual/concepts/06-reduced-optimization-methods.md).
The solver API and policy combinations are covered by the dedicated
optimization reference.

### Use the executable model only when needed

`executable_model()` exposes the backend-neutral residual/objective contract
used by formulations:

```cpp
const auto &model = compiled.executable_model();

auto residual = model.residual(full_point);
auto derivative = model.objective_derivative(full_point);
```

Application code normally does not need to call these operations directly
when it is using `make_reduced_dto()`. They are useful for verification,
custom formulations, and numerical diagnostics.

### State and adjoint services

`state_adjoint_solvers()` returns the exact services used by the reduced
formulation:

```cpp
const auto &solves = compiled.state_adjoint_solvers();
```

They carry `LinearSolveReport` evidence and the compiler-selected lifetime
dependencies.

Prefer `make_reduced_dto()` to reassembling the compiler's model, partition,
and solvers manually unless a custom formulation actually needs the separate
ports.

## Consuming a compiled KKT problem

A successful `quadratic_kkt` request returns:

```cpp
const auto &compiled_kkt = *result.kkt_problem;
const auto &product = compiled_kkt.product();
const auto &manifest = compiled_kkt.manifest();
```

The product is the common
`EqualityConstrainedQuadraticKKTProductT<Backend>` interface. It exposes the
KKT actions, pairings, multiplier/adjoint conversion, assumptions, and
transpose capability required by the KKT solver layer.

For the mathematical block structure, multiplier conventions, and distinction
between reduced and all-at-once systems, see
[Optimality systems and KKT](../manual/concepts/07-optimality-systems-and-kkt.md).

The current compiler accepts this request only for:

- the canonical unconstrained scalar DTO target; or
- the canonical scalar supplied-OTD target with its declared conversion.

A generic reduced problem is not automatically a quadratic KKT product.

## Consuming a compiled PDAS problem

A successful `pdas` request returns:

```cpp
const auto &compiled_pdas = *result.pdas_problem;

const auto &product = compiled_pdas.product();
const auto &complementarity = compiled_pdas.complementarity();
const auto &metric = compiled_pdas.metric();
const auto &box = compiled_pdas.box_data();

const auto &outer_policy = compiled_pdas.pdas_policy();
const auto &inner_policy = compiled_pdas.kkt_solver_policy();
```

The compiled object keeps the KKT product, metric, projection constraint,
box-complementarity representation, and shared box data consistent.

To construct the generic PDAS solver, provide the concrete KKT solve action:

```cpp
auto solver = compiled_pdas.make_solver(
  [&](const auto &active_product) {
    return solve_active_kkt(active_product, inner_policy);
  });
```

The exact serial KKT solve adapter belongs to the optimization layer. For the
algorithmic roles of the KKT and active-set pieces, see
[Optimality systems and KKT](../manual/concepts/07-optimality-systems-and-kkt.md)
and [Complementarity and PDAS](../manual/concepts/08-complementarity-and-pdas.md).

The current compiler PDAS path requires the registered cellwise box and a
positive-diagonal cellwise $L^{2}$ metric. Continuous controls, facewise boxes,
and unrelated metric realizations are not silently converted to this path.

## Consuming a supplied OTD problem

With an all-at-once supplied-OTD semantic graph and the default
`reduced_dto` product request, the successful result is:

```cpp
const auto &compiled_otd = *result.supplied_otd_problem;
const auto &system = compiled_otd.system();
const auto &manifest = compiled_otd.manifest();
```

The result does not contain `CompiledProblemT` and therefore does not expose
`make_reduced_dto()`.

The system provides its own residual, JVP, VJP, and solve operations according
to the supplied-OTD block declaration.

If the caller requests `quadratic_kkt` or `pdas`, the compiler first validates
the canonical supplied-OTD KKT evidence and returns the corresponding KKT/PDAS
product instead.

## Compiled application view

Some compiled reduced targets expose a `CompiledApplicationViewT<Backend>`
beside the erased numerical model.

Check it explicitly:

```cpp
if (const auto *view = compiled.compiled_application_view()) {
  const auto dimensions = view->dimensions();

  // After solving:
  view->write_native_output(
    output_directory,
    state,
    control,
    adjoint,
    uncontrolled_state_or_null);
}
```

The view may also expose objective-component evaluation.

This is intentionally **not** part of `ExecutableModelT`. The numerical
formulation does not gain filesystem or native finite-element ownership merely
because the compiler can retain an application-facing output callback.

## Inspect the compiled service before solving

A first integration should inspect the product once after compilation. This
catches incorrect assumptions about dimensions and optional capabilities
before an optimizer starts calling it repeatedly.

```cpp
const auto &compiled = *result.problem;

const auto &model = compiled.executable_model();
const auto &metric = compiled.metric();

std::cout
  << "variable blocks: "
  << model.variable_layout()->n_blocks() << '\n'
  << "test blocks: "
  << model.test_layout()->n_blocks() << '\n'
  << "metric: "
  << metric.id() << '\n';

if (compiled.constraint())
  std::cout << "projection constraint: available\n";

if (compiled.reduced_hessian())
  std::cout << "reduced Hessian action: available\n";

if (compiled.compiled_application_view())
  std::cout << "native application output view: available\n";
```

Do not hard-code optional capability assumptions from the recipe name. Query
the compiled product that will actually be passed downstream.

For exact finite-element names, realized dimensions, data provenance, and
compiler-selected policies, inspect the manifest rather than trying to infer
them from the semantic graph.

## Compilation manifest

Every successful compiled product carries a `CompilationManifest`.

Use the typed `resolved_decision` records when code needs to inspect what was
actually realized:

```cpp
const auto &manifest = compiled.manifest();
const auto &decision = manifest.resolved_decision;

std::cout
  << "problem: " << decision.semantic_problem_id << '\n'
  << "target: " << decision.target_id << '\n'
  << "mesh cells: " << decision.mesh_record.active_cells << '\n'
  << "mesh provenance: " << decision.mesh_record.provenance << '\n';

for (const auto &space : decision.spaces) {
  std::cout
    << space.semantic_id << ": "
    << space.finite_element << ", dimension "
    << space.dimension << '\n';
}

for (const auto &binding : decision.bindings) {
  std::cout
    << binding.semantic_id << " <- "
    << binding.provenance << '\n';
}
```

The manifest records, among other things:

- semantic problem/formulation identity;
- mesh dimension, active cells, provenance, lifetime mode, and structural
  identity;
- realized semantic spaces;
- runtime data bindings and their provenance;
- pairings;
- residual/observation/loss/transformation realizations;
- realized maps and output dimensions;
- state and adjoint solve policies;
- metric realization;
- constraint realization;
- KKT/PDAS records when selected;
- typed boundary/transposition/trace policies and assumptions.

`manifest.compatibility` is a human-readable compatibility projection of the
typed decision. It is useful for artifacts and diagnostics, but it is not a
second executable configuration.

Do not reconstruct what the compiler did from a problem factory name. Inspect
the manifest when the realized FE, dimension, policy, or provenance matters.

## Current deal.II realization at a glance

The compiler is intentionally registered rather than universally
compositional. The following summary is useful when deciding whether an
existing semantic graph is near a supported path.

### Volume-control families

The baseline scalar component path supports homogeneous or fixed-data
Dirichlet state realization, full-volume or material-subdomain tracking,
$H^{1}$ state tracking, registered scalar residual contributions, cellwise
`FE_DGQ(0)` control, an $L^{2}$ metric, and an optional cellwise box.

The general scalar Robin path additionally supports tensor diffusion,
conservative/advective transport, reaction, Robin bilinear/source terms, and
their typed boundary realization.

### Continuous volume control

Registered continuous-control families use independent homogeneous-Dirichlet
conforming Lagrange coefficients. Hypercube meshes use `FE_Q`; simplex meshes
use `FE_SimplexP`.

Current registered combinations cover:

- $L^{2}$ tracking with an $L^{2}$ metric;
- $H^{1}$ state tracking with $L^{2}$ or $H^{-1}$ metric;
- $H^{1}$ control regularization with $L^{2}$ or $H^{1}$ metric.

These paths do not acquire a cellwise box merely because the control has the
same scalar value type.

### Boundary control

Neumann control supports:

- one facewise-constant coefficient per selected boundary face; or
- a continuous degree-one nodal trace.

The metric and optional box must match that realization. The optional box is
registered for the facewise path, not for the continuous trace path.

Dirichlet-control families use nodal trace coordinates and closed
registrations for the complete/partial boundary and the supported
$L^{2}$/$H^{1/2}$/$H^{1}$ variants.

### Point and flux observations

Point-sensor and normal-flux tracking are specialized registrations with typed
transposition policies. They are not generic observation plugins that can be
attached to every target family.

### Coefficient identification

The registered coefficient-identification path treats a positive cellwise
diffusion parameter as the decision variable, reassembles the state/adjoint
operators, uses a parameter $L^{2}$ metric, and requires its cellwise box.

### Mesh reference cells

Continuous-control and Neumann-boundary-control registrations support the
registered `FE_Q` hypercube and `FE_SimplexP` simplex realizations. Other
current targets remain hypercube `FE_Q` paths.

Mixed reference-cell meshes are not a general fallback.

## Registered factory lookup

For a new application it is often faster to start from the closest complete
semantic family than to assemble a graph from isolated registered kinds.

Useful current factories include:

```cpp
sem::make_scalar_diffusion_reaction_problem(...)
sem::make_fixed_dirichlet_scalar_diffusion_reaction_problem(...)
sem::make_subdomain_tracking_scalar_diffusion_reaction_problem(...)
sem::make_h1_state_tracking_scalar_diffusion_reaction_problem()
sem::make_l2_state_tracking_continuous_control_problem()
sem::make_l2_metric_h1_state_tracking_continuous_control_problem()
sem::make_hminus1_metric_h1_state_tracking_scalar_diffusion_reaction_problem()
sem::make_h1_regularised_scalar_diffusion_reaction_problem()
sem::make_h1_metric_scalar_diffusion_reaction_problem()
sem::make_coefficient_identification_problem()
sem::make_neumann_boundary_control_problem(...)
sem::make_neumann_convection_subdomain_tracking_problem(...)
sem::make_weighted_boundary_trace_neumann_control_problem(...)
sem::make_pure_neumann_boundary_control_problem()
sem::make_dirichlet_control_scalar_diffusion_reaction_problem()
sem::make_partial_dirichlet_control_scalar_diffusion_reaction_problem()
sem::make_l2_dirichlet_laplace_control_problem()
sem::make_hhalf_dirichlet_laplace_control_problem()
sem::make_h1_tracking_hhalf_dirichlet_laplace_control_problem()
sem::make_h1_dirichlet_laplace_control_problem()
sem::make_point_sensor_scalar_diffusion_reaction_problem(...)
sem::make_normal_flux_scalar_diffusion_reaction_problem(...)
sem::make_general_scalar_elliptic_robin_problem(...)
```

The factory list is a lookup aid, not a guarantee of arbitrary composition.
Run `validate()` on the final graph and policy.

## Capability checks that commonly reject a request

### Matrix-free execution

Current v1 compilation is assembled only.

```cpp
policy.execution =
  comp::DealiiDiscretisationPolicy::Execution::matrix_free;
```

produces a lowerability diagnostic.

### State degree zero

`state_degree == 0` is rejected. Current state/test registrations require a
conforming scalar Lagrange degree of at least one.

### Wrong runtime binding shape

The compiler checks the concrete representation expected by each resolved
semantic datum. Supplying a vector-valued `Function` to a scalar port is not
silently accepted by taking one component.

### Missing bound data

Declaring a semantic box but omitting the corresponding runtime bounds is a
lowerability error.

Supplying bounds to an unrelated graph does not create a semantic constraint.

### Unsupported KKT request

`CompilationProduct::quadratic_kkt` is not a universal conversion from every
DTO. The current compiler accepts only its registered canonical scalar DTO or
canonical supplied-OTD target.

### Unsupported PDAS request

PDAS additionally needs the registered box and metric realization. A facewise
box or continuous control is not converted to the cellwise PDAS path.

### Specialized graph outside its closed signature

Some Dirichlet-control, transposition, point/flux, and boundary targets match a
complete registered signature. A graph may be structurally valid yet fail
lowerability if one space, pairing, observation, metric, or typed realization
moves it outside that signature.

This is intentional: the compiler does not choose the “closest” implementation.

## Lifetime checklist

For the owned-session path:

```text
DealiiCompilationSession
        │
        ├── triangulation
        │
        └── retained compiler resources
                   │
                   ▼
             compiled product
                   │
                   ├── callbacks
                   ├── metric
                   ├── solves
                   └── detached formulation products
```

The session protects mesh lifetime, but it does not automatically own arbitrary
objects referenced by `DealiiDataBindings`. Keep bound deal.II functions and
other referenced runtime inputs alive while callbacks may still access them.

For the borrowed-triangulation overload, the caller additionally owns the mesh
lifetime and immutability requirement.

When using `make_reduced_dto()`, the returned service retains the compiler
lifetime owner supplied by `CompiledProblemT`. It is safe to detach it from the
local compilation result with respect to compiler-owned session state. That
does not manufacture ownership of unrelated raw references captured elsewhere.

## Minimal failure-handling pattern

A practical application should treat diagnostics as normal control flow:

```cpp
template <typename Result>
void
print_diagnostics(const Result &result)
{
  for (const auto &d : result.diagnostics.diagnostics()) {
    std::cerr
      << "[" << static_cast<int>(d.category) << "] "
      << d.component_id << '\n'
      << "  capability: " << d.capability << '\n'
      << "  remedy: " << d.remedy << '\n';
  }
}

const auto result = compiler.compile(
  spec,
  session,
  data,
  policy,
  bounds,
  facewise_bounds,
  product);

if (!result.succeeded()) {
  print_diagnostics(result);
  return;
}
```

Do not continue by choosing whichever product pointer happens to be present
after a failed result.

## Compilation checkpoint

Before handing a product to optimization or experiment code, you should know:

- which `CompilationResultT` member is the selected product;
- who owns the triangulation and bound runtime data;
- whether the product has a projection constraint;
- whether an explicit reduced Hessian is available;
- which metric realization was compiled;
- whether state and adjoint solves are iterative or direct;
- where to read the realized dimensions and provenance in the manifest.

You should not need to know which private compiler class assembled a residual
term or which internal planner branch selected the target. Those are
implementation details unless you are extending the compiler itself.

## Where the compiler boundary ends

A successful compilation has already decided and recorded:

- concrete FE spaces and coefficient orderings;
- runtime data placement;
- metric and optional constraint realization;
- state/adjoint solve services;
- optional reduced-Hessian capability;
- requested KKT/PDAS/supplied-OTD product;
- lifetime retention;
- provenance manifest.

It has **not** chosen:

- steepest descent versus BFGS versus Newton;
- line search versus trust region;
- stopping tolerances for the outer optimization;
- experiment run-set organization;
- output-directory policy.

Those decisions belong to the optimization and application/experiment
interfaces. See
[Reduced optimization methods](../manual/concepts/06-reduced-optimization-methods.md)
for the algorithmic background and
[Application execution](application-execution.md) for the current run,
artifact, and experiment interfaces.
