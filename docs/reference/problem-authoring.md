# Problem authoring reference

This reference is for code that constructs `nmopt`'s semantic problem description.
It starts from mathematical ingredients that have already been chosen and shows
how to represent them as a `semantic::v1::ProblemSpec`, validate the graph, and
package reusable families as `ProblemRecipeT`.

The semantic layer describes **what the problem is**. It does not own a mesh,
deal.II `Function` objects, finite-element degrees, linear-solver tolerances, or
optimization algorithms. Those values enter later through the compiler and
application layers.

For the conceptual model behind this API, see
[Semantic problem model](../manual/concepts/09-semantic-problem-model.md) and
[Validation, resolution, and capabilities](../manual/concepts/10-validation-resolution-and-capabilities.md).
This page focuses on what to write and on the relationships the API checks.

## Public surface

Most authoring code needs:

```cpp
#include "nmopt/semantic/v1/problem_spec.hpp"
```

This compatibility aggregate includes the public semantic types, validation,
resolution, and the ready-made problem-library factories.

Reusable application families additionally use:

```cpp
#include "nmopt/application/recipe.hpp"
```

The examples below use:

```cpp
namespace sem = nmopt::semantic::v1;
namespace app = nmopt::application;
```

## Choose the authoring boundary first

There are three common ways to enter the compiled path.

```text
existing registered family
        │
        ▼
ProblemRecipeT<Parameters>
        │
        │ build(parameters)
        ▼
     ProblemSpec
        │
        ├──────────── hand-authored or composed ProblemSpec
        │
        ▼
semantic validation
        │
        ▼
deal.II compiler
```

Use a `ProblemRecipeT` when several applications or experiments should share
one semantic family while varying a small set of meaningful parameters.

Construct or compose a `ProblemSpec` directly when defining a new semantic
family, extending a registered one, or writing focused tests.

Do not construct a `ProblemSpec` merely to wrap an already existing PDE
application that owns its numerical realization. That integration path starts
from numerical callbacks instead; see
[Integrating an existing PDE application](../manual/concepts/13-integrating-an-existing-pde-application.md)
for the conceptual path and
[External deal.II solver integration](external-dealii-solver-integration.md)
for the current programming interface.

## Working path: obtain, inspect, vary, and validate a problem

Most application code should not begin by filling every vector in
`ProblemSpec`. Start from the closest registered problem family, expose the
semantic variations your application actually needs, and inspect the resulting
graph before moving to compilation.

The example in this section stays with one scalar distributed-control problem
from beginning to end.

### Step 1 – choose the closest semantic family

```cpp
#include "nmopt/application/chapter5.hpp"
#include "nmopt/semantic/v1/problem_spec.hpp"

#include <iostream>

namespace ch5 = nmopt::application::chapter5;
namespace sem = nmopt::semantic::v1;

// A recipe is a typed builder for one semantic problem family.
//
// No mesh, deal.II Function, finite-element object, or solver is created here.
// We are choosing the mathematical structure that the compiler will later
// realize numerically.
const auto recipe = ch5::make_scalar_distributed_recipe();
```

`make_scalar_distributed_recipe()` represents the family

```text
state equation: scalar diffusion–reaction
decision:       distributed volume control
objective:      state tracking + control regularization
formulation:    reduced DTO
```

The recipe has already encoded the regions, spaces, weak residual terms,
observations, losses, metric, formulation, and compiler policies that define
that family.

This is the first useful authoring question:

> **Does an existing recipe already describe the mathematics I need?**

If yes, use its typed parameters instead of copying its `ProblemSpec` and
editing component IDs manually.

### Step 2 – choose semantic variations

The distributed-control recipe accepts
`ScalarDistributedControlParameters`.

```cpp
ch5::ScalarDistributedControlParameters problem_parameters;

// Keep the default cellwise-constant distributed control.
//
// This choice changes the semantic control realization. It is not an
// optimizer choice and it is not the value of the control.
problem_parameters.discretisation =
  ch5::DistributedControlDiscretisation::cellwise_constant;

// Do not request a box constraint in the first run.
//
// If this becomes true, the recipe adds the lower/upper bound data ports,
// the ConstraintSpec, the required discrete box policy, and the formulation's
// constraint reference. The numerical bound values are still supplied later.
problem_parameters.with_cellwise_box = false;
```

A recipe parameter belongs here when changing it changes the semantic problem
family. The actual forcing function, target values, regularization coefficient,
mesh refinement, state-solve tolerance, and optimization stopping criterion do
not belong in this record.

### Step 3 – materialize the semantic graph

```cpp
// Calling the recipe builds an ordinary backend-independent ProblemSpec.
//
// From this point onward `spec` is the complete semantic request. It contains
// stable IDs and references among the problem components, but still contains
// no deal.II objects and no realized vector dimensions.
const sem::ProblemSpec spec = recipe(problem_parameters);
```

**At this point:** the mathematical/compositional request exists, but no
numerical problem has been built.

That distinction is useful in application architecture. A recipe and its
parameters can be selected, inspected, validated, or catalogued without first
creating a mesh.

### Step 4 – inspect the pieces that matter to the application

You do not normally need to print or traverse every component. Inspect the
parts whose identity affects later application decisions.

```cpp
std::cout
  << "problem id: " << spec.id << '\n'
  << "problem label: " << spec.label << '\n'
  << "formulation: " << spec.formulation.id << '\n'
  << "decision variable: "
  << spec.formulation.control_variable_id << '\n'
  << "metric: " << spec.formulation.metric_id << '\n';

if (spec.formulation.constraint_id.empty())
  std::cout << "constraint: none\n";
else
  std::cout << "constraint: "
            << spec.formulation.constraint_id << '\n';
```

For this recipe the primary decision variable is named `"control"`. The field
is called `control_variable_id` for historical reasons; coefficient
identification uses the same formulation port for a semantic parameter.

When writing generic application code, prefer the IDs stored in the
formulation and component records over assumptions such as “the second
variable must be the control.”

### Step 5 – validate the graph before touching the backend

```cpp
const sem::ValidationReport semantic_report =
  sem::SemanticValidator{}.validate(spec);

if (!semantic_report.valid()) {
  for (const sem::Diagnostic &diagnostic :
       semantic_report.diagnostics()) {
    std::cerr
      << "component: " << diagnostic.component_id << '\n'
      << "capability: " << diagnostic.capability << '\n'
      << "remedy: " << diagnostic.remedy << "\n\n";
  }

  // Do not proceed to backend binding with a structurally invalid graph.
  return;
}
```

A successful semantic validation means that references, roles, pairings,
policies, and formulation structure are internally coherent.

It does **not** yet mean that the current deal.II compiler can realize the
graph. That second question belongs to `DealiiCompiler::validate()` and needs
a discretization policy.

### Step 6 – make a supported semantic variation through the recipe

Suppose the same application now needs a cellwise box.

Change the semantic choice at the recipe boundary:

```cpp
problem_parameters.with_cellwise_box = true;

const sem::ProblemSpec bounded_spec =
  recipe(problem_parameters);
```

The recipe changes several graph components together. Conceptually the delta
is:

```text
unconstrained graph
      │
      ├── add lower-bound data port
      ├── add upper-bound data port
      ├── add cellwise-box constraint
      ├── add box-realization requirement
      └── set formulation.constraint_id
              │
              ▼
        bounded semantic graph
```

The recipe does **not** choose numerical lower and upper values. During
compilation you will additionally provide `CellwiseBoxDataBindings`.

This is why recipe parameters are preferable to scattered mutations when a
variation is already registered: one typed choice updates the whole semantic
contract coherently.

### Step 7 – understand incompatible semantic choices early

The same family also supports a continuous homogeneous-Dirichlet control:

```cpp
problem_parameters.discretisation =
  ch5::DistributedControlDiscretisation::
    homogeneous_dirichlet_continuous;
```

That representation does not support the cellwise box. Therefore this is an
invalid family request:

```cpp
problem_parameters.with_cellwise_box = true;

// The recipe rejects this combination instead of constructing a graph whose
// control and box realizations disagree.
const sem::ProblemSpec invalid =
  recipe(problem_parameters);
```

The recipe throws `std::invalid_argument` for that incompatible typed
combination.

This is different from compiler lowerability. A recipe can reject a
nonsensical combination of its own parameters before a `ProblemSpec` exists;
the semantic validator checks graph coherence after construction; the compiler
then checks whether a coherent graph has a registered numerical realization.

### Step 8 – hand the semantic request to the compiler

For the original unconstrained example, restore:

```cpp
problem_parameters.discretisation =
  ch5::DistributedControlDiscretisation::cellwise_constant;
problem_parameters.with_cellwise_box = false;

const sem::ProblemSpec compilable_spec =
  recipe(problem_parameters);
```

The next layer will combine `compilable_spec` with four things that are
deliberately absent here:

```text
                         ProblemSpec
                             │
          ┌──────────────────┼────────────────────┐
          │                  │                    │
          ▼                  ▼                    ▼
        mesh         concrete runtime data   discretization policy
          │                  │                    │
          └──────────────────┼────────────────────┘
                             ▼
                       DealiiCompiler
```

The [compiler reference](compiler.md) continues this same example.

## When to author a `ProblemSpec` directly

The recipe-first path is intentionally the shortest route. Write or compose a
`ProblemSpec` directly when you are:

- defining a new reusable semantic family;
- adding a component combination for which no recipe exists;
- implementing a new compiler capability and its focused tests;
- constructing a very small graph for contract testing.

In that case, the full example below is the useful starting point. It shows
the relationships that a recipe normally hides.

## A complete scalar problem

The following function constructs the same kind of graph as the canonical
scalar diffusion-reaction volume-control problem. It is intentionally written
out rather than hidden behind a helper so the connections among components are
visible.

```cpp
sem::ProblemSpec
make_volume_control_problem()
{
  sem::ProblemSpec spec;

  // Give the whole graph a stable programmatic identity and a separate
  // human-readable label. Other components do not refer to this label.
  spec.id = "my.scalar-volume-control";
  spec.label = "Scalar diffusion-reaction volume control";

  // Regions answer "where does this object live?" independently of the
  // finite-element mesh. Boundary and material IDs are semantic selectors
  // that the deal.II compiler later resolves against the concrete mesh.
  spec.regions = {
    {"domain",
     "Full volume domain",
     sem::RegionKind::volume,
     true,
     {},
     {},
     {}},

    {"dirichlet_boundary",
     "Homogeneous Dirichlet boundary",
     sem::RegionKind::boundary,
     false,
     {0},
     {},
     {}}
  };

  // Spaces describe mathematical/semantic roles, not concrete FE_Q/FE_DGQ
  // objects. Their dimensions are normally unknown until compilation.
  spec.spaces = {
    {"state_space",
     "State",
     "domain",
     sem::SpaceTopology::h1,
     sem::SpaceRole::state},

    {"state_test_space",
     "State test",
     "domain",
     sem::SpaceTopology::h1,
     sem::SpaceRole::test},

    {"control_space",
     "Cellwise control",
     "domain",
     sem::SpaceTopology::l2,
     sem::SpaceRole::control},

    {"state_observation_space",
     "State observation",
     "domain",
     sem::SpaceTopology::l2,
     sem::SpaceRole::observation},

    {"control_observation_space",
     "Control observation",
     "domain",
     sem::SpaceTopology::l2,
     sem::SpaceRole::observation}
  };

  // Pairings make primal/covector compatibility explicit. They do not
  // choose a Riesz map; the optimization metric is declared separately.
  spec.pairings = {
    {"state_pairing",
     "State coefficient pairing",
     "state_space",
     "state_space"},

    {"state_test_pairing",
     "State-test coefficient pairing",
     "state_test_space",
     "state_test_space"},

    {"control_pairing",
     "Control coefficient pairing",
     "control_space",
     "control_space"},

    {"state_observation_pairing",
     "State-observation coefficient pairing",
     "state_observation_space",
     "state_observation_space"},

    {"control_observation_pairing",
     "Control-observation coefficient pairing",
     "control_observation_space",
     "control_observation_space"}
  };

  // Variables connect the formulation unknowns to semantic spaces.
  // The empty transformation ID means these coordinates already represent
  // the physical field consumed by this baseline problem.
  spec.variables = {
    {"state",
     "State",
     sem::VariableRole::state,
     "state_space",
     ""},

    {"control",
     "Control",
     sem::VariableRole::control,
     "control_space",
     ""}
  };

  // Data nodes declare runtime input ports. They name what must later be
  // bound, but they deliberately contain no dealii::Function or numeric
  // coefficient value.
  spec.data = {
    {"forcing",
     "Volume forcing",
     sem::DataKind::function,
     sem::DataRole::forcing,
     "state_test_space"},

    {"desired_state",
     "Desired state",
     sem::DataKind::function,
     sem::DataRole::desired_state,
     "state_observation_space"},

    {"diffusion",
     "Diffusion coefficient",
     sem::DataKind::scalar_constant,
     sem::DataRole::diffusion,
     ""},

    {"reaction",
     "Reaction coefficient",
     sem::DataKind::scalar_constant,
     sem::DataRole::reaction,
     ""},

    {"regularisation_weight",
     "Control regularisation",
     sem::DataKind::scalar_constant,
     sem::DataRole::regularisation_weight,
     ""}
  };

  // Residual terms are the compositional pieces of the weak equation.
  // Each term explicitly names the variables and immutable data it consumes.
  spec.residual_terms = {
    {"diffusion_reaction",
     "Diffusion and reaction",
     sem::ResidualTermKind::diffusion_reaction,
     "state_equation",
     {"state"},
     {"diffusion", "reaction"},
     ""},

    {"volume_source",
     "Volume source",
     sem::ResidualTermKind::volume_source,
     "state_equation",
     {},
     {"forcing"},
     ""},

    {"volume_control",
     "Volume control",
     sem::ResidualTermKind::volume_control,
     "state_equation",
     {"control"},
     {},
     ""}
  };

  // The equation closes the residual by naming the test space, pairing,
  // and exact ordered set of terms that belong to this weak block.
  spec.equations = {
    {"state_equation",
     "State residual",
     "state_test_space",
     "state_test_pairing",
     {"diffusion_reaction", "volume_source", "volume_control"}}
  };

  // Observations are objective-facing maps. Keeping them separate from the
  // state residual allows the same PDE to be combined with different
  // tracking operators without rewriting the equation.
  spec.observations = {
    {"state_observation",
     "Full-domain state restriction",
     sem::ObservationKind::volume_restriction,
     "state",
     "domain",
     "state_observation_space",
     "state_observation_pairing",
     {}},

    {"control_observation",
     "Full-domain control restriction",
     sem::ObservationKind::volume_restriction,
     "control",
     "domain",
     "control_observation_space",
     "control_observation_pairing",
     {}}
  };

  // Losses consume observation outputs and immutable data. The scalar
  // regularisation weight remains a data port until compilation.
  spec.losses = {
    {"state_tracking",
     "Quadratic tracking",
     sem::LossKind::quadratic_tracking,
     "state_observation",
     "desired_state",
     "state_observation_pairing"},

    {"control_regularisation",
     "Quadratic control regularisation",
     sem::LossKind::quadratic_control_regularisation,
     "control_observation",
     "regularisation_weight",
     "control_observation_pairing"}
  };

  // The search metric is a separate semantic choice from objective
  // regularisation. Here both happen to be L2-related, but the API does not
  // identify those two concepts.
  spec.metrics = {
    {"control_l2_metric",
     "Cellwise L2 metric",
     sem::MetricKind::l2,
     "control",
     "control_pairing"}
  };

  // Requirement policies carry assumptions or discrete realization
  // choices that cannot be recovered from graph connectivity alone.
  spec.requirement_policies = {
    {"state_fixed_dirichlet",
     "state",
     sem::RequirementKind::fixed_dirichlet,
     sem::RequirementStatus::selected_discrete_realisation,
     sem::RequirementScope::discrete_compilation,
     "homogeneous full-vector Dirichlet rows",
     "dirichlet_boundary"},

    {"desired_state_quadrature_policy",
     "desired_state",
     sem::RequirementKind::analytic_quadrature_evaluation,
     sem::RequirementStatus::selected_discrete_realisation,
     sem::RequirementScope::discrete_compilation,
     "analytic Function evaluated at selected volume quadrature",
     "domain"}
  };

  // Finally, select how the graph is exposed to numerical consumers.
  // The reduced DTO formulation names the state, primary decision variable,
  // state equation, search metric, and optional constraint.
  spec.formulation = {
    "reduced_dto",
    sem::FormulationKind::reduced_dto,
    sem::FormulationProvenance::dto,
    "state",
    "control",
    "state_equation",
    "control_l2_metric",
    ""
  };

  return spec;
}
```

The important property of this code is not the aggregate-initializer syntax.
It is that every relation is explicit.

```text
region ───────────────► space
                         │
                         ├────────► variable
                         │
                         ├────────► data
                         │
                         └────────► pairing
                                      │
variable + data ──► residual term ────┤
                         │            │
                         ▼            │
                      equation ◄──────┘
                         │
variable ──► observation ──► loss
   │                          │
   ├────────► metric          │
   └────────► constraint?     │
                              │
                 formulation ◄┘
```

A compiler does not recover these edges from names or labels. It consumes the
IDs stored in the graph.

## How to read an unfamiliar `ProblemSpec`

When opening application code that returns a `ProblemSpec`, do not start by
reading every vector in declaration order. A faster inspection sequence is:

1. read `spec.formulation`;
2. resolve the state and decision variables named by it;
3. inspect the referenced equation and metric;
4. follow the equation's residual-term IDs;
5. inspect losses and their observations;
6. inspect any constraint and requirement policies;
7. only then inspect auxiliary spaces/data needed by those paths.

For example:

```cpp
const auto &formulation = spec.formulation;

std::cout
  << "state: " << formulation.state_variable_id << '\n'
  << "decision: " << formulation.control_variable_id << '\n'
  << "equation: " << formulation.equation_id << '\n'
  << "metric: " << formulation.metric_id << '\n';
```

For repeated lookup, resolve the graph once:

```cpp
const sem::SemanticResolution resolution =
  sem::SemanticResolver{}.resolve(spec);

if (!resolution.succeeded())
  return;

const sem::ResolvedProblemView &problem = *resolution.problem;

// Follow the formulation outward instead of searching vectors manually.
const auto &state =
  problem.variable(formulation.state_variable_id);

const auto &decision =
  problem.variable(formulation.control_variable_id);

const auto &equation =
  problem.equation(formulation.equation_id);

const auto &metric =
  problem.metric(formulation.metric_id);
```

This is also a useful pattern for tooling and agents: the formulation is the
entry point into the graph; stable IDs are the edges.

## Stable IDs are the graph API

Every reusable component has an `id` and usually a human-readable `label`.
Treat them differently:

- `id` is the programmatic identity used by references from other components;
- `label` is presentation text and may change without rewiring the graph.

For example, this equation refers to the **IDs** of its test space, pairing,
and residual terms:

```cpp
sem::EquationBlockSpec equation{
  "state_equation",
  "State residual",
  "state_test_space",
  "state_test_pairing",
  {"diffusion_reaction", "volume_source", "volume_control"}
};
```

Changing the label `"State residual"` is harmless. Changing
`"state_test_space"` without changing the corresponding `SpaceSpec::id`
disconnects the graph and produces a structural diagnostic.

Use stable, descriptive IDs in reusable code. Avoid generating IDs from vector
indices, pointer values, or display labels.

## Component responsibilities and references

The semantic graph has several component kinds because different relationships
need to remain explicit. The following table is intended as a lookup while
writing a graph; it omits fields that are merely descriptive.

| Component | What you choose | References it carries |
| --- | --- | --- |
| `RegionSpec` | volume, boundary, or point set and its selector data | none |
| `SpaceSpec` | topology, semantic role, region | `region_id` |
| `PairingSpec` | declared primal/covector pairing | two space IDs |
| `VariableSpec` | state/control/parameter role and coordinate space | `space_id`, optional physical-field transformation |
| `DataSpec` | runtime datum kind and mathematical role | optional `space_id` |
| `TransformationSpec` | reconstruction/lifting from optimization coordinates | variable, output space, fixed datum, optional control |
| `ResidualTermSpec` | one contribution to an equation | equation, variables, data, optional region |
| `EquationBlockSpec` | test block and its residual terms | test space, test pairing, term IDs |
| `ObservationSpec` | map used by a loss | variable, region, output space/pairing, optional data |
| `LossSpec` | objective contribution | observation, datum, pairing |
| `MetricSpec` | search geometry for one decision variable | variable and pairing |
| `ConstraintSpec` | admissible set | variable and bound-data IDs |
| `RequirementPolicySpec` | assumption or selected realization | subject, optional region, optional typed realization |
| `FormulationSpec` | solver-facing interpretation of the graph | state, decision, equation, metric, optional constraint |

`ProblemSpec` owns vectors of these components plus one `FormulationSpec` and
an optional supplied-OTD declaration. It is a composition root, not a
polymorphic PDE base class.

## Regions

`RegionSpec` has three current `RegionKind` values:

- `volume`;
- `boundary`;
- `point_set`.

A full-domain volume region normally sets `is_full_domain = true`. A boundary
region names deal.II boundary IDs in `boundary_ids`. A material subdomain uses
`material_ids`.

A point-set region stores physical coordinates directly:

```cpp
sem::RegionSpec sensors{
  "sensor_points",
  "Observation sensor locations",
  sem::RegionKind::point_set,
  false,
  {},
  {},
  {{0.25, 0.5}, {0.75, 0.5}}
};
```

The semantic graph stores those coordinates without a deal.II point type. The
compiler later checks their dimension against the concrete mesh.

Do not use a region label as a selector. The compiler consumes the typed
boundary IDs, material IDs, or point coordinates.

## Spaces, pairings, and variables

`SpaceSpec` records semantic topology and role, not a concrete finite element.

The currently represented topologies are `h1`, `h2`, `hhalf`, `l2`, and
`bounded_function`. The role states why the space exists: state, test,
control, parameter, observation, data, or auxiliary.

Most spaces leave `dimension` at zero because the dimension is determined only
after discretization. A finite-dimensional point-sensor observation space is
the notable current case that declares a positive semantic dimension.

A `PairingSpec` records which semantic primal and covector roles are paired.
It is not a Riesz map and it does not choose an optimization metric. A metric
is declared separately through `MetricSpec`.

A `VariableSpec` binds an optimization variable to a semantic space:

```cpp
sem::VariableSpec state{
  "state",
  "State",
  sem::VariableRole::state,
  "state_space",
  ""
};
```

The final field, `physical_field_transform_id`, is empty when the optimization
coordinates already represent the physical field. It names a
`TransformationSpec` when the physical field must first be reconstructed.

## Physical-field transformations

Use `TransformationSpec` when the variable stored by the optimization problem
is not itself the complete physical field consumed by the residual or
observation.

For fixed nonzero Dirichlet data, the semantic relationship is:

```text
independent state coordinates
          │
          ├──────────────┐
          │              │ fixed datum
          ▼              ▼
     fixed Dirichlet reconstruction
                  │
                  ▼
           physical state field
```

The corresponding declaration has this shape:

```cpp
spec.data.push_back({
  "fixed_dirichlet_data",
  "Fixed Dirichlet data",
  sem::DataKind::function,
  sem::DataRole::fixed_dirichlet_lifting,
  "state_space"
});

spec.transformations.push_back({
  "fixed_dirichlet_reconstruction",
  "Fixed Dirichlet reconstruction",
  sem::TransformationKind::fixed_dirichlet_reconstruction,
  "state",                  // input_variable_id
  "state_space",            // output_space_id
  "fixed_dirichlet_data",   // fixed_data_id
  ""                        // no control input
});

state.physical_field_transform_id = "fixed_dirichlet_reconstruction";
```

Controlled Dirichlet lifting uses the other registered transformation kind,
`dirichlet_control_lifting`, and fills `control_variable_id` as the second
primal input.

Do not represent an essential-control lifting as an unrelated boundary source
term. If the physical state depends on a transformation of optimization
coordinates, make that map explicit.

## Residual terms and equations

A residual term declares one semantic contribution. Its ports are explicit:

```cpp
sem::ResidualTermSpec term{
  "tensor_diffusion",
  "Tensor diffusion",
  sem::ResidualTermKind::tensor_diffusion,
  "state_equation",
  {"state"},
  {"diffusion_tensor"},
  ""
};
```

`variable_ids` names unknowns consumed by the term. `data_ids` names immutable
runtime inputs. `region_id` is used when the contribution belongs to a
particular boundary or subdomain.

`EquationBlockSpec` owns the ordered list of residual-term IDs that constitute
one weak equation block. A term that exists in `ProblemSpec::residual_terms`
but is not listed by an equation is not implicitly assembled into that
equation.

The current deal.II capability registry recognizes these residual kinds:

- `laplacian`;
- `diffusion_reaction`;
- `tensor_diffusion`;
- `conservative_transport`;
- `advective_transport`;
- `reaction`;
- `parameter_diffusion_reaction`;
- `transposition_laplacian`;
- `dirichlet_transposition_control`;
- `volume_source`;
- `volume_control`;
- `neumann_control`;
- `robin_bilinear`;
- `robin_source`;
- `natural_boundary_source`.

Recognition of a kind does **not** mean every arbitrary mixture of recognized
terms can be compiled. The compiler accepts bounded component compositions and
a set of closed specialized registrations. Use `DealiiCompiler::validate()`
to test the complete graph and requested product.

## Observations and losses

An `ObservationSpec` is the map whose result is consumed by an objective loss.
This separation lets the same state equation be combined with different
tracking maps without disguising the observation as part of the residual.

A standard full-volume observation is:

```cpp
sem::ObservationSpec observation{
  "state_observation",
  "Full-domain state restriction",
  sem::ObservationKind::volume_restriction,
  "state",
  "domain",
  "state_observation_space",
  "state_observation_pairing",
  {}
};
```

Weighted observations name immutable data explicitly:

```cpp
observation.data_ids = {"boundary_weight"};
```

The current compiler recognizes volume restriction, $H^{1}$ state
restriction, boundary trace/restriction, weighted boundary trace, point
sensor, and normal-flux observations.

A `LossSpec` connects an observation to the datum against which it is compared
or to the scalar coefficient that weights it:

```cpp
sem::LossSpec tracking{
  "state_tracking",
  "Quadratic tracking",
  sem::LossKind::quadratic_tracking,
  "state_observation",
  "desired_state",
  "state_observation_pairing"
};
```

Current registered loss kinds are quadratic tracking, ordinary quadratic
control regularization, $H^{1/2}$ control regularization, $H^{1}$ control
regularization, and quadratic parameter regularization.

The semantic graph does not store the numerical value of the desired state or
regularization coefficient. It declares their ports. Concrete values are
bound during compilation.

## Metrics and constraints

A `MetricSpec` identifies the search metric attached to the primary decision
variable:

```cpp
sem::MetricSpec metric{
  "control_l2_metric",
  "Control L2 metric",
  sem::MetricKind::l2,
  "control",
  "control_pairing"
};
```

The current deal.II compiler recognizes `l2`, `hhalf`, `h1`, and `hminus1`
metric kinds. The metric choice is independent of the objective loss. A graph
may therefore use one topology for regularization and another registered
metric for optimization.

A box constraint needs both semantic structure and runtime bound values. The
semantic part is:

```cpp
spec.data.push_back({
  "lower_bound",
  "Control lower bound",
  sem::DataKind::cellwise_bound,
  sem::DataRole::lower_bound,
  "control_space"
});

spec.data.push_back({
  "upper_bound",
  "Control upper bound",
  sem::DataKind::cellwise_bound,
  sem::DataRole::upper_bound,
  "control_space"
});

spec.constraints.push_back({
  "control_box",
  "Cellwise control box",
  sem::ConstraintKind::cellwise_box,
  "control",
  "lower_bound",
  "upper_bound"
});

spec.requirement_policies.push_back({
  "control_box_policy",
  "control_box",
  sem::RequirementKind::discrete_cellwise_bounds,
  sem::RequirementStatus::selected_discrete_realisation,
  sem::RequirementScope::discrete_compilation,
  "FE_DGQ(0) coefficientwise clipping in l2_cellwise",
  "domain"
});

spec.formulation.constraint_id = "control_box";
```

The corresponding lower and upper values are **not** placed in `ProblemSpec`.
They are supplied later as `CellwiseBoxDataBindings`.

`facewise_box` is a distinct constraint kind for facewise boundary controls.
Do not substitute cellwise and facewise bounds merely because both happen to
be coefficient vectors.

## Requirement policies

Some graph properties cannot be represented only by component connectivity.
`RequirementPolicySpec` records those assumptions or selected realizations.

The status is significant:

- `provided` records evidence represented directly by the model;
- `user_assumed` records an assumption the model author is asserting;
- `selected_discrete_realisation` records an explicit realization required for
  compilation.

The scope is also explicit: continuous semantics, discrete compilation, or
both.

A free-form `selected_policy` string is descriptive. It is not a substitute
for typed realization data when a capability has a typed selection.

For example, the general scalar Robin family records its boundary realization
as a `BoundaryRealisationSelection`:

```cpp
sem::RequirementPolicySpec boundary_policy{
  "scalar_boundary_partition",
  "state",
  sem::RequirementKind::boundary_partition,
  sem::RequirementStatus::selected_discrete_realisation,
  sem::RequirementScope::both,
  "fixed Dirichlet plus Robin/transport outflow partition",
  ""
};

boundary_policy.typed_selection = sem::BoundaryRealisationSelection{
  "scalar_boundary_partition",
  "state",
  "dirichlet_boundary",
  "robin_boundary",
  {},  // no Neumann regions
  {},  // no separate transport-inflow regions
  "robin_boundary",
  sem::ConormalForm::diffusion_minus_transport,
  sem::NormalOrientation::outward,
  sem::TraceEvaluationRealisation::fe_q_state_trace,
  sem::FaceQuadratureRealisation::qgauss_face
};

spec.requirement_policies.push_back(boundary_policy);
```

Other current capabilities have dedicated typed fields on
`RequirementPolicySpec`, including:

- `typed_trace_selection`;
- `typed_neumann_control_selection`;
- `typed_metric_selection`;
- `typed_transposition_selection`;
- `typed_partial_boundary_selection`;
- `typed_fractional_metric_selection`;
- `typed_boundary_h1_metric_selection`;
- `typed_h1_target_data_membership_selection`.

Use the dedicated typed record when the selected realization affects
lowerability. Do not encode its semantics only in `selected_policy`.

## Formulation selection

`FormulationSpec` says how the semantic graph is intended to be consumed.

The ordinary reduced path has:

```cpp
spec.formulation = {
  "reduced_dto",
  sem::FormulationKind::reduced_dto,
  sem::FormulationProvenance::dto,
  "state",
  "control",
  "state_equation",
  "control_l2_metric",
  ""  // optional constraint ID
};
```

`control_variable_id` is the name retained from the first control-only slice.
It is the primary decision-variable port and may identify a semantic
`parameter` in coefficient-identification problems.

The current all-at-once path uses:

```cpp
sem::FormulationKind::all_at_once
sem::FormulationProvenance::supplied_otd
```

and additionally requires `ProblemSpec::supplied_otd_declaration`.

Changing the formulation labels on a reduced DTO graph is not enough to create
a supplied OTD system. The declaration must describe the state, adjoint, and
control-stationarity blocks, their semantic and runtime space IDs, pairings,
action provenance, multiplier convention, and comparison status.

For the registered scalar supplied-OTD family, prefer the existing factory:

```cpp
sem::ProblemSpec spec =
  sem::make_scalar_diffusion_reaction_supplied_otd_problem();
```

Hand-author a `SuppliedOTDDeclaration` only when adding a genuinely new
registered formulation, because the compiler validates its block structure and
conventions as a closed capability.

## Ready-made semantic families

`problem_library.hpp` provides complete graphs used by applications and
contract tests. They are useful starting points and, more importantly, show
the exact semantic signatures that the current compiler understands.

The principal families are:

- scalar diffusion-reaction volume control, with optional cellwise box;
- continuous distributed control with $L^{2}$ or $H^{1}$ control geometry;
- $H^{1}$ state tracking with $L^{2}$ or $H^{-1}$ control metric;
- coefficient identification with a positive cellwise parameter;
- fixed Dirichlet reconstruction;
- complete and partial Dirichlet control;
- the registered $L^{2}$, $H^{1/2}$, and $H^{1}$ Dirichlet-control variants;
- Neumann boundary control with facewise or continuous trace coordinates;
- conservative-transport Neumann control with subdomain tracking;
- weighted boundary-trace tracking;
- pure-Neumann control with mean constraint;
- point-sensor tracking;
- normal-flux tracking;
- general scalar elliptic/Robin composition.

Representative factories include:

```cpp
sem::make_scalar_diffusion_reaction_problem(...)
sem::make_l2_state_tracking_continuous_control_problem()
sem::make_hminus1_metric_h1_state_tracking_scalar_diffusion_reaction_problem()
sem::make_coefficient_identification_problem()
sem::make_neumann_boundary_control_problem(...)
sem::make_neumann_convection_subdomain_tracking_problem(...)
sem::make_dirichlet_control_scalar_diffusion_reaction_problem()
sem::make_partial_dirichlet_control_scalar_diffusion_reaction_problem()
sem::make_hhalf_dirichlet_laplace_control_problem()
sem::make_h1_dirichlet_laplace_control_problem()
sem::make_point_sensor_scalar_diffusion_reaction_problem(...)
sem::make_normal_flux_scalar_diffusion_reaction_problem(...)
sem::make_general_scalar_elliptic_robin_problem(...)
```

A ready-made factory is not an inheritance point. It returns an ordinary
`ProblemSpec`. You may compose a new graph from its pieces, but validation is
performed on the resulting graph and the compiler does not promise that every
cross-product of individually registered pieces is lowerable.

## Recipes package semantic variation

Use `ProblemRecipeT<Parameters>` when callers should select meaningful problem
choices without editing graph nodes directly.

The Chapter 5 distributed-control recipe follows this pattern:

```cpp
enum class ControlDiscretisation
{
  cellwise_constant,
  homogeneous_dirichlet_continuous
};

struct ProblemParameters
{
  bool with_cellwise_box = false;
  ControlDiscretisation control =
    ControlDiscretisation::cellwise_constant;
};

app::ProblemRecipeT<ProblemParameters> recipe{
  {
    "my.scalar-control",
    "Scalar control family",
    "Scalar diffusion-reaction volume control",
    "application",
    {"scalar", "volume-control"}
  },

  [](const ProblemParameters &p) {
    switch (p.control) {
      case ControlDiscretisation::cellwise_constant:
        return sem::make_scalar_diffusion_reaction_problem(
          p.with_cellwise_box);

      case ControlDiscretisation::homogeneous_dirichlet_continuous:
        if (p.with_cellwise_box)
          throw std::invalid_argument(
            "continuous control does not support the cellwise box");

        return sem::make_l2_state_tracking_continuous_control_problem();
    }

    throw std::invalid_argument("unknown control discretisation");
  }
};

ProblemParameters parameters;
parameters.with_cellwise_box = true;

sem::ProblemSpec spec = recipe(parameters);
```

A recipe should expose **semantic family choices**:

- control or parameter representation;
- optional semantic constraint;
- selected boundary/subdomain IDs;
- other choices that change the problem graph.

It should not own:

- a triangulation;
- concrete forcing or target `Function` objects;
- state/adjoint solver tolerances;
- optimizer iteration limits;
- output directories.

Those belong to later layers.

## Validate the semantic graph

Use `SemanticValidator` when you want to check only semantic structure and
declared policies:

```cpp
sem::SemanticValidator validator;
const sem::ValidationReport report = validator.validate(spec);

if (!report.valid()) {
  for (const sem::Diagnostic &d : report.diagnostics()) {
    std::cerr
      << "component: " << d.component_id << '\n'
      << "capability: " << d.capability << '\n'
      << "remedy: " << d.remedy << "\n\n";
  }
}
```

Semantic validation checks, among other things:

- nonempty stable IDs and labels;
- duplicate IDs;
- required enum selections;
- referenced regions, spaces, variables, data, pairings, terms, and
  observations;
- role and topology compatibility;
- equation/term ownership;
- observation/loss compatibility;
- metric and constraint attachment;
- formulation references;
- supplied-OTD declaration structure;
- required analytical/discrete policies.

A valid semantic graph is not necessarily compilable by the current deal.II
backend. Backend lowerability and formulation-product checks are added by
`DealiiCompiler::validate()`.

`DiagnosticCategory` separates the source of a problem:

- `structural` – the graph is internally incomplete or inconsistent;
- `analytical_policy` – a required assumption or selected policy is missing;
- `lowerability` – the concrete compiler cannot realize the valid request;
- `formulation_capability` – the requested formulation/product is unsupported.

The first two are produced by semantic validation. The compiler appends the
latter two to the same report type.

## Resolve IDs for graph-processing code

Application code normally passes `ProblemSpec` directly to the compiler.
Code that needs repeated typed lookup inside a validated graph can use
`SemanticResolver`:

```cpp
sem::SemanticResolution resolution =
  sem::SemanticResolver{}.resolve(spec);

if (!resolution.succeeded())
  return;

const sem::ResolvedProblemView &problem = *resolution.problem;

const auto &state = problem.variable("state");
const auto &space = problem.space(state.space_id);
const auto &equation =
  problem.equation(spec.formulation.equation_id);
```

`ResolvedProblemView` borrows the original `ProblemSpec`; it does not copy or
own it. Keep the specification alive for the full lifetime of the view.

Do not use the resolver as a second authoring mechanism. Change the
`ProblemSpec`, then resolve it again.

## Registered semantic kinds versus registered compositions

The v1 deal.II kind registry currently recognizes the following individual
semantic kinds:

| Family | Registered kinds |
| --- | --- |
| residual | laplacian; diffusion-reaction; tensor diffusion; conservative and advective transport; reaction; parameter diffusion-reaction; transposition Laplacian; Dirichlet transposition control; volume source/control; Neumann control; Robin bilinear/source; natural boundary source |
| observation | volume restriction; $H^{1}$ state restriction; boundary trace/restriction; weighted boundary trace; point sensor; normal flux |
| loss | quadratic tracking; $L^{2}$, $H^{1/2}$, and $H^{1}$ control regularization; parameter regularization |
| metric | $L^{2}$; $H^{1/2}$; $H^{1}$; $H^{-1}$ |
| constraint | cellwise box; facewise box |
| transformation | fixed Dirichlet reconstruction; Dirichlet-control lifting |

This table answers only: “does the compiler know this semantic kind?”

It does **not** answer: “can these kinds be combined arbitrarily?” Several
families are closed registrations with specific spaces, pairings, policies,
and boundary realizations. In particular, the Dirichlet-control variants,
point/flux observations, supplied OTD, KKT, and PDAS paths have additional
whole-graph requirements.

Always validate the complete graph with the compilation product you intend to
request:

```cpp
nmopt::compiler::v1::DealiiCompiler compiler;
nmopt::compiler::v1::DealiiDiscretisationPolicy policy;

const auto diagnostics = compiler.validate(
  spec,
  policy,
  nmopt::compiler::v1::CompilationProduct::reduced_dto);
```

The compiler reports an unsupported composition; it does not silently replace
a component with the nearest registered alternative.

## Common authoring mistakes

### Putting runtime values in `ProblemSpec`

Do not store a `dealii::Function`, triangulation, matrix, or mesh-dependent
vector in the semantic graph. Declare the corresponding `DataSpec` and bind
the concrete value at compilation time.

### Using labels as references

This is wrong in reusable code:

```cpp
equation.test_space_id = "State test";  // label, not ID
```

Use the stable component ID:

```cpp
equation.test_space_id = "state_test_space";
```

### Treating metric and regularization as one choice

A `LossSpec` determines an objective contribution. A `MetricSpec` determines
the search geometry. They may intentionally differ.

### Adding a bound without a constraint policy

A pair of lower/upper `DataSpec` objects does not by itself create an
admissible set. Add the `ConstraintSpec`, the required realization policy, and
the formulation's `constraint_id`.

### Copying a factory and changing only enum labels

A registered target is defined by a complete structural signature. Changing
one enum or formulation label can produce a semantically valid but currently
unlowerable graph. Use compiler validation as part of authoring.

### Encoding a typed realization only as prose

When `RequirementPolicySpec` has a dedicated typed selection for the chosen
capability, fill it. `selected_policy` is useful description, not executable
evidence.

## Authoring checkpoint

Before leaving the semantic layer, you should be able to answer these
questions from the graph itself:

- Which variable is the state?
- Which variable or parameter is the primary optimization decision?
- Which weak equation determines the state?
- Which immutable data ports must be supplied later?
- Which observations and losses make up the objective?
- Which metric defines the decision-space search geometry?
- Is there a semantic constraint?
- Which assumptions or discrete realization policies are explicitly declared?
- Which formulation is requested?

You should **not** yet need answers to:

- How many DoFs does the state have?
- Which deal.II finite element was constructed?
- Which concrete `Function` object supplies the target?
- What CG tolerance will the state solve use?
- What vector should the optimizer start from?

Those are downstream numerical/application choices.

## From authoring to compilation

At the end of this layer you should have:

```cpp
sem::ProblemSpec spec = recipe(parameters);

const auto semantic_report =
  sem::SemanticValidator{}.validate(spec);

if (!semantic_report.valid())
  // Fix the graph before binding backend data.
```

The next layer supplies the mesh, concrete runtime data, discretization policy,
and desired compilation product:

```text
ProblemSpec
    │
    │ semantic meaning
    ▼
DealiiCompiler
    ▲
    │ concrete mesh + data + policy + bounds
    │
runtime application
```

See [Compiler reference](compiler.md) for that boundary.
