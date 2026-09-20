# Application authoring reference

This reference explains how to author a **reusable nmopt-native application
family** so that downstream users can select a scenario, supply a `.prm`
configuration, and run it without reconstructing the semantic graph or backend
integration by hand.

It fills the layer between the semantic/compiler APIs and the execution/config
APIs:

```text
user supplies .prm / CLI selections
        │
        ▼
runner / configuration integration
        │
        ▼
typed application scenario
        │
        ├── problem parameters ──► ProblemRecipeT ──► ProblemSpec
        │                                      │
        └── compile / solver / experiment ─────┤
                                               ▼
                                    backend execution adapter
                                               │
                                               ▼
                                    compiler + optimizer
```

Use [Problem authoring](problem-authoring.md) when defining the mathematical
`ProblemSpec` itself. Use [Compiler](compiler.md) for the deal.II lowering API.
Use [Application execution](application-execution.md) once the application has
already been authored and should be run. If you already have a complete PDE
application and want to keep it application-owned, use the peer path in
[External deal.II solver integration](external-dealii-solver-integration.md)
instead.

## What belongs in an application family

A reusable application family normally has these pieces:

```text
semantic family
    typed recipe parameters
    ProblemRecipeT

application family
    problem/runtime parameter record
    compile option record
    solver option record
    experiment option record
    ScenarioT alias
    scenario validation
    scenario factory
    ProblemSpec builder

backend realization
    runtime-data adapter
    mesh/session builder
    compiler-policy mapper
    execution adapter

repository runner integration (optional)
    parameter schema adapter
    parameter -> scenario binder
    benchmark metadata registration
    artifact-coordinate planner
    execution registration
```

Those layers should remain separable. In particular, the semantic recipe
should not know about `.prm` syntax, and the parameter-file parser should not
know how to assemble finite-element operators.

## Start with the semantic recipe

`ProblemRecipeT<Parameters>` is the reusable boundary between typed semantic
choices and a `ProblemSpec`:

```cpp
#include "nmopt/application/recipe.hpp"

struct MyRecipeParameters
{
  bool with_cellwise_box = false;
};

using MyRecipe =
  nmopt::application::ProblemRecipeT<MyRecipeParameters>;
```

Construct the recipe with metadata and one builder:

```cpp
MyRecipe
make_my_recipe()
{
  return MyRecipe{
    {
      "my.scalar-volume-control",
      "My scalar volume-control problem",
      "Reusable scalar diffusion-reaction control family",
      "my-applications",
      {"scalar", "volume-control"}
    },
    [](const MyRecipeParameters &parameters) {
      return nmopt::semantic::v1::
        make_scalar_diffusion_reaction_problem(
          parameters.with_cellwise_box);
    }
  };
}
```

The builder should contain **semantic** choices only.

For example, `chapter5::ScalarDistributedControlParameters` chooses whether
the control is cellwise or continuous and whether a compatible cellwise box is
part of the problem. It does not carry a triangulation, a deal.II `Function`,
linear-solver tolerances, or output directories.

If your application family matches an existing Chapter 5 recipe, reuse it
instead of creating a wrapper recipe with identical semantics.

## Decide what the application problem record owns

A scenario-level problem record usually combines:

```text
recipe parameters
        +
runtime-data definitions/provenance
        +
application-specific semantic selections
```

The Chapter 6 B1 record is a good example:

```cpp
struct B1ProblemParameters
{
  chapter5::ScalarDistributedControlParameters recipe;
  ScalarRuntimeDataOptions                     data;
  std::vector<double>                          regularisation_sweep;
  ScalarFunctionDefinition                     forcing;
  ScalarFunctionDefinition                     desired_state;
};
```

The important boundary is that `forcing` and `desired_state` here are still
**typed definitions**, not deal.II `Function` objects.

For the running `MyScenario` examples below, use the same pattern:

```cpp
struct MyProblemParameters
{
  MyRecipeParameters recipe;

  nmopt::application::chapter6::ScalarRuntimeDataOptions data;

  nmopt::application::ScalarFunctionDefinition forcing;
  nmopt::application::ScalarFunctionDefinition desired_state;
};
```

The `chapter6::ScalarRuntimeDataOptions` reuse is only convenient because this
example has the same scalar diffusion/reaction/regularisation ports. A new
application should define its own runtime-data record when its ports differ.

That lets the same scenario be:

- validated before deal.II construction;
- built from a parameter file;
- listed in a catalog;
- compared or transformed by tests;
- instantiated by a backend-specific runtime-data adapter later.

If an application needs to augment an existing recipe, do that explicitly in a
problem-spec builder rather than smuggling the extra semantics into runtime
bindings.

Chapter 6 B2 demonstrates this pattern:

```cpp
semantic::v1::ProblemSpec
make_b2_problem_spec(const B2ProblemParameters &parameters)
{
  auto spec =
    chapter5::make_neumann_convection_recipe()(
      parameters.recipe);

  // Application-family semantic augmentation:
  // - concrete boundary-region policy
  // - optional natural source port
  // - fixed Dirichlet reconstruction
  // - application-specific requirement selections

  return spec;
}
```

The result is still an ordinary backend-neutral `ProblemSpec`.

## Define the option records separately

`ScenarioT` does not prescribe the option record types:

```cpp
template <
  typename ProblemParameters,
  typename CompileOptions,
  typename SolverOptions,
  typename ExperimentOptions>
struct ScenarioT;
```

Use that freedom to keep each layer typed.

If an application's option surface matches an existing reusable record, reuse
it directly. For example, a scalar Chapter-6-style family can use:

```cpp
using MyCompileOptions =
  nmopt::application::chapter6::CompileOptions;

using MySolverOptions =
  nmopt::application::chapter6::SolverOptions;

using MyExperimentOptions =
  nmopt::application::chapter6::ExperimentOptions;
```

Those records already separate mesh/compiler choices, reduced-solver choices,
and experiment/output choices.

If a new application needs a genuinely different surface, define its own small
records with the same responsibility split rather than extending the Chapter 6
types with unrelated fields. The exact records are application API: their job
is to represent the choices a caller is allowed to make, not every field
exposed by lower layers.

That distinction is important. The low-level optimizer may support Newton,
Wolfe, PDAS, and trust region, while one application scenario may intentionally
expose only assembled reduced DTO + Armijo + L-BFGS.

Do not widen the scenario option record merely because a lower-level enum has
more values.

## Define the scenario type

```cpp
using MyScenario =
  nmopt::application::ScenarioT<
    MyProblemParameters,
    MyCompileOptions,
    MySolverOptions,
    MyExperimentOptions>;
```

A `ScenarioT` contains:

```text
metadata
problem
compile
solver
experiment
```

and nothing executable.

The scenario is the typed application contract that a parameter binder or C++
caller should produce.

## Give the scenario stable metadata

```cpp
nmopt::application::ScenarioMetadata metadata{
  "my-app.case-a",
  "My application – case A",
  "Volume-control scenario used by the application",
  "my-applications",
  "my.scalar-volume-control",
  {"assembled", "serial-dealii"}
};
```

`recipe_id` should name the semantic family actually used by the problem-spec
builder.

`ApplicationCatalog` later uses this metadata for discovery; do not overload
the scenario ID with filesystem coordinates or run-slot numbers.

## Write scenario validation at the application boundary

`ScenarioT::validate()` checks common metadata only. Application invariants
belong in a family validator:

```cpp
void
validate_my_scenario(const MyScenario &scenario)
{
  scenario.validate();

  validate_my_compile_options(scenario.compile);
  validate_my_solver_options(scenario.solver);
  validate_my_runtime_data(scenario.problem.data);

  if (scenario.compile.execution !=
      nmopt::application::chapter6::ExecutionSelection::assembled)
    throw std::invalid_argument(
      "my application currently requires assembled execution");

  if (scenario.compile.product !=
      nmopt::application::chapter6::ProductSelection::reduced_dto)
    throw std::invalid_argument(
      "my application currently exposes the reduced DTO product");
}
```

This is where combinations of otherwise valid lower-level values become
application-invalid.

Examples from the current Chapter 6 families include:

- B1 rejecting a cellwise box with continuous homogeneous-Dirichlet control;
- B1 restricting simplex meshes to the continuous-control target;
- B2 requiring volume-observation policy;
- B2 requiring distinct fixed/control/outflow boundary IDs;
- B1 exposing only steepest descent or L-BFGS;
- B2 exposing BFGS while allowing selected globalization choices;
- both freezing assembled reduced-DTO compilation.

Keeping those rules in the scenario validator gives parameter files, C++
callers, and tests the same application contract.

## Provide a scenario factory

A scenario factory should create one complete, internally consistent default
case:

```cpp
MyScenario
make_my_scenario()
{
  MyProblemParameters problem;
  problem.recipe.with_cellwise_box = false;

  problem.forcing =
    nmopt::application::chapter6::b1_manufactured_zero_forcing();
  problem.desired_state =
    nmopt::application::chapter6::b1_manufactured_desired_state();

  problem.data.diffusion = 1.0;
  problem.data.reaction = 0.0;
  problem.data.regularisation_weight = 1.0e-2;
  problem.data.forcing_provenance = problem.forcing.provenance;
  problem.data.desired_state_provenance = problem.desired_state.provenance;

  MyCompileOptions compile;
  compile.mesh.dimension = 2;
  compile.mesh.refinement = 4;

  MySolverOptions solver;
  solver.method =
    nmopt::application::chapter6::ReducedMethod::steepest_descent;
  solver.parameters.gradient_tolerance = 1.0e-8;

  MyExperimentOptions experiment;
  experiment.scenario_output_id = "my-app.case-a";
  experiment.source_reference = "my application definition";
  experiment.source_revision = "1";
  experiment.build_profile = "release-dealii";

  MyScenario scenario{
    {
      "my-app.case-a",
      "My application – case A",
      "Volume-control scenario used by the application",
      "my-applications",
      "my.scalar-volume-control",
      {"assembled", "serial-dealii"}
    },
    std::move(problem),
    std::move(compile),
    std::move(solver),
    std::move(experiment)
  };

  validate_my_scenario(scenario);
  return scenario;
}
```

Downstream configuration code should start from this typed object and override
declared fields. It should not construct an unvalidated aggregate independently
in every call site.

## Provide one `ProblemSpec` builder for the scenario family

If the scenario uses an existing recipe directly:

```cpp
nmopt::semantic::v1::ProblemSpec
make_my_problem_spec(
  const MyProblemParameters &parameters)
{
  return make_my_recipe()(parameters.recipe);
}
```

Also provide a scenario overload when useful:

```cpp
nmopt::semantic::v1::ProblemSpec
make_my_problem_spec(const MyScenario &scenario)
{
  validate_my_scenario(scenario);
  return make_my_problem_spec(scenario.problem);
}
```

If the application augments the base graph, this function is the place to do
it. Keep runtime values out of the graph unless they are genuinely semantic
choices.

For detailed graph authoring rules, use
[Problem authoring](problem-authoring.md).

## Add the application to the metadata catalog

If the family should be discoverable:

```cpp
nmopt::application::ApplicationCatalog
make_catalog()
{
  nmopt::application::ApplicationCatalog catalog;

  catalog.add(make_my_recipe().metadata());
  catalog.add(make_my_scenario().metadata);

  return catalog;
}
```

The catalog is metadata-only. It should be usable without constructing a mesh,
compiler, runtime function, or optimizer.

This is why catalog registration is separate from executable runner
registration.

## Add the backend realization

The application-facing scenario is still backend-neutral. For a deal.II-backed
application, put the concrete realization in a backend-specific layer, such as:

```text
include/nmopt/application/dealii/my_application.hpp
```

The current Chapter 6 code uses this layer for four distinct tasks.

### Convert typed runtime definitions to deal.II objects

For example, a scalar function definition can be selected in the scenario and
realized later as a `dealii::Function<dim>`.

The current Chapter 6 helpers follow the pattern:

```cpp
MyRuntimeDataT<dim>
make_my_runtime_data(
  const MyScenario &scenario,
  const MySelectedDataT<dim> &selected)
{
  validate_my_scenario(scenario);

  // Check that selected concrete functions match the scenario definitions.
  // Return references/owners required by compiler bindings.
}
```

The scenario remains the source of requested identity/provenance; the backend
adapter creates the concrete objects.

### Build compiler data bindings

```cpp
nmopt::compiler::v1::DealiiDataBindings<dim>
make_my_data_bindings(
  const MyProblemParameters &parameters,
  const MyRuntimeDataT<dim> &runtime)
{
  return {
    runtime.forcing,
    runtime.desired_state,
    parameters.data.diffusion,
    parameters.data.reaction,
    parameters.data.regularisation_weight,
    {
      parameters.data.forcing_provenance,
      parameters.data.desired_state_provenance,
      ""
    }
  };
}
```

The exact constructor/optional ports depend on the semantic graph.

### Build the compilation session

```cpp
template <int dim>
std::shared_ptr<
  nmopt::compiler::v1::DealiiCompilationSession<dim>>
make_my_compilation_session(
  const MyScenario &scenario)
{
  validate_my_scenario(scenario);

  auto mesh =
    std::make_unique<dealii::Triangulation<dim>>();

  // Map scenario.compile.mesh to the application-supported mesh policy.

  return std::make_shared<
    nmopt::compiler::v1::DealiiCompilationSession<dim>>(
      std::move(mesh),
      scenario.compile.mesh.mesh_provenance);
}
```

Mesh generation belongs here because it is a numerical realization, not a
semantic recipe operation.

### Map application compile options to compiler policy

```cpp
nmopt::compiler::v1::DealiiDiscretisationPolicy
make_my_discretisation_policy(
  const MyCompileOptions &options)
{
  validate_my_compile_options(options);

  nmopt::compiler::v1::DealiiDiscretisationPolicy policy;

  policy.state_degree = options.state_degree;
  policy.execution =
    nmopt::compiler::v1::DealiiDiscretisationPolicy::
      Execution::assembled;

  policy.state_solve = {
    options.state_solve.maximum_iterations,
    options.state_solve.relative_tolerance,
    options.state_solve.absolute_tolerance
  };

  return policy;
}
```

This mapper is intentionally explicit. Scenario enums are application-facing;
compiler enums are compiler-facing. The mapping is where unsupported
application combinations should be rejected, not silently reinterpreted.

## Write one execution adapter

The backend execution adapter turns:

```text
ProblemSpec + Scenario
```

into detached benchmark execution evidence.

The Chapter 6 pattern is:

```cpp
using Envelope =
  nmopt::experiment::ReducedSearchExperimentEnvelopeT<Backend>;

using Evidence =
  nmopt::application::benchmark::
    BenchmarkExecutionEvidenceT<Envelope>;

class MyReducedExecutionAdapter final
{
public:
  Evidence
  operator()(
    const nmopt::semantic::v1::ProblemSpec &spec,
    const MyScenario &scenario) const
  {
    validate_my_scenario(scenario);

    const auto bindings =
      make_my_data_bindings(
        scenario.problem,
        runtime_);

    const auto policy =
      make_my_discretisation_policy(
        scenario.compile);

    nmopt::compiler::v1::DealiiCompiler compiler;

    const auto compilation =
      compiler.compile(
        spec,
        session_,
        bindings,
        policy,
        std::nullopt,
        std::nullopt,
        nmopt::compiler::v1::CompilationProduct::reduced_dto);

    if (!compilation.succeeded())
      return failure_evidence(compilation.diagnostics);

    auto reduced =
      compilation.problem->make_reduced_dto();

    const auto report =
      solve_selected_method(
        reduced,
        compilation.problem->metric(),
        scenario.solver);

    return make_detached_evidence(
      compilation.problem->manifest(),
      report,
      scenario,
      environment_);
  }
};
```

The helper names in this sketch are application-local. The key is the
responsibility split:

```text
scenario
    says what is requested

execution adapter
    realizes that request using compiler + optimizer

execution evidence
    detaches what needs to survive the numerical objects
```

The current concrete examples are:

- [`B1ReducedExecutionAdapterT`](../../include/nmopt/application/dealii/chapter6_b1.hpp);
- [`B2ReducedExecutionAdapterT`](../../include/nmopt/application/dealii/chapter6_b2.hpp).

## Keep artifact fields out of the core numerical adapter when possible

The execution adapter should return the numerical envelope, validation
diagnostics, measurements, selected fields, and application evidence needed by
the generic harness.

Repository-specific path coordinates or parameter-matrix labels can be added
by the runner layer after execution:

```cpp
auto evidence = execute(specification, scenario);

add_parameter_artifact_fields(
  evidence,
  configuration,
  combination);
```

That keeps one backend execution adapter reusable across multiple run-set
layouts.

## Make the application configurable from `.prm`

Everything up to this point can be used directly from C++.

To let another user supply only `.prm`/CLI configuration, the repository
runner adds an application-private schema/binding layer.

This part currently lives under:

```text
apps/nmopt-runner/
```

rather than `include/nmopt/` because it is runner policy, not a universal
framework API.

### Register the parameter schema

`ParameterSchemaAdapter` says which benchmark/recipe IDs select the schema,
which matrix axes are valid, and which scalar-definition catalogues can be
discovered.

Conceptually:

```cpp
ParameterSchemaAdapter my_adapter{
  "my-app",
  {"my-app"},
  {"my-app."},
  {"my.scalar-volume-control"},
  {},
  scalar_definition_schemas
};

append_matrix_axis(my_adapter, "method");
append_matrix_axis(my_adapter, "regularisation");
```

Do not add arbitrary keys to a generic catch-all parser. Add the application
schema deliberately so unsupported configuration fails early.

See [Parameter files](parameter-files.md) for the user-facing syntax.

### Bind one resolved combination to the typed scenario

The binder should start from the scenario factory and override typed fields:

```cpp
auto scenario = make_my_scenario();

bind_my_scenario(
  scenario,
  file,
  combination,
  scenario_id);

validate_my_scenario(scenario);
```

The current Chapter 6 binders map `.prm` values into:

```text
scenario.problem
scenario.compile
scenario.solver
scenario.experiment
```

and leave compiler/backend object creation to the execution adapter.

This is the central reason for having a typed scenario layer: parameter files
never need to know `DealiiDataBindings` constructor order or solver template
arguments.

## Register the runnable benchmark

The current repository runner has two separate registrations.

Metadata registration:

```cpp
struct BenchmarkRegistration
{
  std::string_view id;
  std::string_view parameter_benchmark_id;
  std::string_view default_parameter_file;
};
```

Execution registration:

```cpp
struct BenchmarkExecutionRegistration
{
  const BenchmarkRegistration *metadata;
  BenchmarkArtifactPlanner     artifact_planner;
  BenchmarkExecutionCallback   execute;
};
```

For a new runnable family, add both.

The artifact planner turns one resolved matrix combination into stable path
components. The execution callback receives the already validated
`RunSetPlan`, parameter file, run configuration, and `RunSetManifest`.

The current B1 execution flow is representative:

```text
resolved matrix combination
        │
        ▼
make_b1_scenario(...)
        │
        ▼
bind_b1_scenario(...)
        │
        ▼
make runtime data + compilation session
        │
        ▼
B1ReducedExecutionAdapterT
        │
        ▼
HeadlessBenchmarkRunnerT::run
        │
        ▼
artifact.kv + native output + solver trace
        │
        ▼
RunSetManifest::record_success/failure
```

The concrete orchestration is in
[`chapter6_execution.hpp`](../../apps/nmopt-runner/chapter6_execution.hpp).

## Decide where a new choice belongs

When extending an application, use the layer that owns the meaning of the
choice.

A change belongs in the **recipe/ProblemSpec** when it changes mathematical
composition, for example:

```text
cellwise vs continuous control
new observation operator
new constraint
new residual term
new semantic boundary partition
```

A change belongs in **scenario compile options** when it changes numerical
realization without changing the semantic problem, for example:

```text
mesh generation/refinement
FE degree
state/adjoint solve tolerances
supported execution mode
```

A change belongs in **scenario solver options** when it changes optimization
policy:

```text
steepest descent vs L-BFGS
Armijo vs fixed step
stopping tolerances
initial control
```

A change belongs in **experiment options** when it changes run/evidence
behavior:

```text
source revision
build profile
retain fields
measure timings
output identity
```

A change belongs in **parameter schema/binding** only when users should be able
to choose an already-supported typed option from `.prm`.

Adding a parameter-file key before the typed application layer supports the
choice reverses the ownership direction and usually creates stringly typed
logic in the runner.

## A complete authoring route

For a new nmopt-native application family, the shortest maintainable route is:

1. Reuse or define a semantic recipe and typed recipe parameters.
2. Define the application problem record with runtime definitions/provenance.
3. Define the allowed compile, solver, and experiment option records.
4. Create a `ScenarioT` alias.
5. Write one application validator that rejects unsupported combinations.
6. Provide a scenario factory with useful, internally consistent defaults.
7. Provide one `ProblemSpec` builder from problem parameters/scenario.
8. Add metadata catalog entries if the family should be discoverable.
9. Implement backend runtime-data/session/policy mappers.
10. Implement one execution adapter returning detached evidence.
11. Only then add `.prm` schema/binding and runner registration if CLI-driven
    execution is desired.

After step 10, the family is usable from C++. After step 11, a downstream user
can operate it through the higher-level workflow documented in
[Application execution](application-execution.md) and
[Parameter files](parameter-files.md) without touching semantic/compiler code.
