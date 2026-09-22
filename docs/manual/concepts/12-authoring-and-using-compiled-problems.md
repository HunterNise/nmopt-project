# Authoring and using compiled problems

## The compiled path is a complete user path, not only an internal architecture

Parts III described how `ProblemSpec` is validated, resolved, lowered, and packaged.

A reader who wants to use that machinery still needs one more view:

> Which objects do I create, which choices belong together, and what do I do with the
> compiled product afterward?

The compiled path has two common entry points.

The first uses a ready-made problem family:

```text
recipe / scenario already exists
        ↓
select parameters
        ↓
bind runtime data
        ↓
compile
        ↓
select formulation product
        ↓
solve
```

The second authors a new family:

```text
mathematical problem
        ↓
typed parameters
        ↓
ProblemRecipeT
        ↓
ProblemSpec
        ↓
runtime/compiler/solver integration
```

Both paths use the same layers. The difference is how much of the setup has already
been packaged for the user.

This chapter uses the source Chapter 5 scalar distributed-control recipe and the Chapter 6
B1 scenario as concrete examples.

## 1. Start by deciding which layer you intend to change

Several objects that look like "problem options" belong to different parts of the
system.

The most important separation is:

| Choice | Owning layer | Example |
| --- | --- | --- |
| mathematical/semantic structure | recipe parameters | cellwise versus continuous control, add a box constraint |
| runtime mathematical data | runtime/problem data | forcing function, desired state, coefficient values, regularization weight |
| numerical realization | compile options | mesh, FE degree, state/adjoint solve policy |
| optimization algorithm | solver options | steepest descent, L-BFGS, stopping tolerances |
| experiment/output policy | experiment options | retained fields, output directory, timing collection |

Changing one row should not require rebuilding the others unless the new choice
changes a declared capability.

For example, changing the desired state from one function to another need not change
the semantic graph.

Changing a cellwise control to a continuous control does change the semantic problem
and may select a different compiler target.

## 2. `ProblemRecipeT` owns problem construction, not execution

A recipe is a typed mapping

```math
\text{parameters}
\longmapsto
\texttt{ProblemSpec}.
```

In code, the shape is:

```text
ProblemRecipeT<Parameters>
    metadata
    builder : Parameters -> ProblemSpec
```

It does not own

```text
mesh
deal.II Function objects
compiler
solver
experiment output
```

That narrow responsibility is useful because the same semantic problem family can be
compiled on different meshes, with different data, or solved by different algorithms.

The source Chapter 5 scalar distributed-control recipe has the parameter type

```text
ScalarDistributedControlParameters
    with_cellwise_box
    discretisation
```

where the current discretization choices are

```text
cellwise_constant
homogeneous_dirichlet_continuous
```

The builder dispatches to the corresponding semantic specification.

A small use therefore looks like:

```cpp
auto recipe = chapter5::make_scalar_distributed_recipe();

chapter5::ScalarDistributedControlParameters parameters;
parameters.with_cellwise_box = true;  // semantic choice

ProblemSpec specification = recipe.build(parameters);
```

The result is still only a semantic graph.

No mesh, matrix, or optimization iteration has been created.

## 3. A ready-made recipe is the shortest way to request a supported family

When a repository recipe already describes the mathematics you need, prefer changing
its typed parameters over reconstructing `ProblemSpec` node by node.

For the scalar distributed recipe, the user can currently select:

```text
cellwise constant control
    optional registered cellwise box

or

continuous homogeneous-Dirichlet control
    no cellwise box
```

The recipe rejects the unsupported combination

```text
continuous control + cellwise box
```

before the compiler is involved.

This is a useful division of responsibility.

The recipe knows which semantic variants belong to one named problem family.

The compiler later decides whether the resulting complete graph has a registered
numerical realization.

## 4. A scenario binds problem, compile, solver, and experiment choices

`ScenarioT` is a typed aggregate:

```text
ScenarioT
    metadata
    problem
    compile
    solver
    experiment
```

It does not compile or execute.

A scenario says:

> For this named application case, these are the problem data and these are the
> compile, solve, and experiment choices that should travel together.

The distinction between recipe and scenario is therefore

```text
recipe
    parameters -> semantic problem

scenario
    one application case
        recipe/problem parameters
        + compile options
        + solver options
        + experiment options
```

The source Chapter 6 B1 scenario uses exactly this pattern.

## 5. The B1 scenario separates semantic choice from runtime data

`B1ProblemParameters` contains several kinds of problem-facing data:

```text
recipe
    ScalarDistributedControlParameters

data
    diffusion
    reaction
    regularisation_weight
    provenance strings

forcing
    scalar function definition

desired_state
    scalar function definition

regularisation_sweep
    values used by the benchmark/application case
```

Only `recipe` is sent to the source Chapter 5 `ProblemRecipeT` to determine the semantic
graph.

The forcing, desired state, and coefficient values are bound later as runtime data.

This means the same `ProblemSpec` can represent several runs with different numerical
values while preserving the same mathematical structure.

## 6. `make_b1_scenario()` is a useful ready-made starting point

The repository provides

```cpp
auto scenario = chapter6::make_b1_scenario();
```

which selects a complete B1 application configuration.

Its defaults include a reduced DTO product, assembled execution, a mesh configuration,
solver parameters, runtime data definitions, and experiment metadata.

The scenario is not immutable.

A user can start from it and change the layer they actually intend to study.

For example:

```cpp
auto scenario = chapter6::make_b1_scenario(
  chapter6::ReducedMethod::limited_memory_bfgs);

scenario.compile.mesh.refinement = 6;
scenario.solver.initial_control_value = 0.25;
```

The important point is not the exact values above.

It is that runtime problem data and solver initialization are changed in their own
records instead of being encoded into the semantic graph.

## 7. Scenario validation checks cross-layer application assumptions

A generic `ScenarioT::validate()` validates only scenario metadata.

Application families can add stronger validation.

B1 checks, among other things:

```text
mesh options are coherent

runtime diffusion/reaction/regularization values are finite

selected control discretization is recognized

continuous control is not combined with the cellwise box

the selected solver method is allowed by the B1 application case

the selected product remains assembled reduced DTO

forcing and desired-state definitions agree with their provenance
```

These checks are not substitutes for semantic/compiler validation.

They protect assumptions of the named application case.

The sequence is:

```text
scenario validation
    application-case consistency
        ↓
ProblemSpec construction
        ↓
semantic/compiler validation
    framework capability
```

## 8. Building the `ProblemSpec` from a scenario remains explicit

B1 provides

```cpp
ProblemSpec specification =
  chapter6::make_b1_problem_spec(scenario);
```

Internally, this validates the scenario and calls the source Chapter 5 recipe with

```text
scenario.problem.recipe
```

rather than copying the entire scenario into `ProblemSpec`.

That distinction keeps compile policy and optimizer policy outside the semantic graph.

The resulting graph can therefore still be inspected, validated, or compiled
independently of the scenario object that produced it.

## 9. Runtime function definitions become concrete deal.II objects later

B1 stores forcing and target as backend-neutral scalar function definitions.

The deal.II application adapter later creates the actual

```text
dealii::Function<dim>
```

objects.

The flow is:

```text
ScalarFunctionDefinition
        ↓
B1SelectedDataT
        ↓
B1RuntimeDataT
        ↓
DealiiDataBindings
```

This is another boundary worth preserving.

A recipe/scenario can be discovered and validated without constructing deal.II
objects.

Only the backend adapter needs to turn the selected data definition into a concrete
runtime function.

## 10. `DealiiDataBindings` supplies values to semantic data ports

For B1, the concrete binding contains

```text
forcing Function

desired-state Function

diffusion coefficient

reaction coefficient

regularization weight

provenance
```

The adapter creates it with a function such as

```cpp
auto bindings =
  make_b1_data_bindings(scenario.problem,
                        regularisation,
                        runtime);
```

The binding does not redefine what `"forcing"` means.

That meaning was already declared in `ProblemSpec`.

The binding says which concrete value supplies the declared port for this
compilation.

## 11. Mesh ownership is explicit through `DealiiCompilationSession`

The ready-made B1 path builds an owned compilation session:

```text
mesh options
    ↓
Triangulation
    ↓
DealiiCompilationSession<dim>
```

The session owns the triangulation and records mesh provenance.

This is the simplest lifetime arrangement when the compiled services may retain
objects related to the mesh.

A caller that already owns a triangulation can use the borrowed-triangulation
compiler overload instead, but then the caller must enforce the required lifetime.

The choice is about ownership, not PDE mathematics.

## 12. Backend-neutral compile options are translated to compiler policy

Source Chapter 6 `CompileOptions` does not expose the full compiler implementation type.

The backend adapter translates it to

```text
DealiiDiscretisationPolicy
```

including:

```text
state FE degree

execution realization

state solve policy

adjoint solve policy

control-metric solve policy
```

For B1:

```cpp
auto policy =
  make_b1_discretisation_policy(scenario.compile);

auto product =
  make_b1_compilation_product(scenario.compile);
```

The scenario-level records are convenient application configuration.

The compiler-level records remain the authoritative inputs to
`DealiiCompiler`.

## 13. Compilation joins semantic meaning and runtime realization

At this point the inputs are:

```text
ProblemSpec
    requested semantic problem

DealiiCompilationSession
    owned mesh

DealiiDataBindings
    concrete data values

DealiiDiscretisationPolicy
    numerical realization choices

CompilationProduct
    requested solver-facing product
```

The B1 adapter performs the essential call:

```cpp
DealiiCompiler compiler;

auto compilation =
  compiler.compile(specification,
                   session,
                   bindings,
                   policy,
                   std::nullopt,   // no cellwise bounds in this run
                   std::nullopt,   // no facewise bounds
                   product);
```

A user-level path should treat the result as fallible:

```cpp
if (!compilation.succeeded() || !compilation.problem)
  // inspect diagnostics and stop
```

Do not continue by guessing a nearby supported target.

## 14. Validation can be used before expensive construction

`DealiiCompiler::validate()` accepts

```text
ProblemSpec
DealiiDiscretisationPolicy
CompilationProduct
```

without requiring runtime mesh/data construction.

That makes it useful when an application lets the user choose among several semantic
or product variants.

Conceptually:

```cpp
auto report =
  compiler.validate(specification, policy, product);

if (!report.valid())
  // report diagnostics before creating the full run
```

Compilation still performs the checks again and can additionally fail on concrete
mesh/data conditions that static validation cannot see.

## 15. A successful reduced compilation is immediately usable

For a reduced DTO request:

```cpp
auto &compiled = *compilation.problem;

auto reduced = compiled.make_reduced_dto();
```

The compiled problem also exposes the metric and state/adjoint services needed by
the reduced formulation.

The user does not need to retrieve the compiler's private matrices and assemble a
new formulation object manually.

The compiled product is already the bridge from compiler-owned numerical realization
to the generic optimization layer.

## 16. The initial control must come from the compiled layout

The optimizer operates on the realized control space, not on an arbitrary vector of a
guessed dimension.

B1 obtains the control layout from the compiled model through a state/control
partition.

Conceptually:

```text
compiled executable model
        ↓
StateControlPartitionT
        ↓
control_layout()
        ↓
initial PrimalBlockT
```

A uniform initial value can then be written into the realized control coefficients.

This matters because changing the control discretization can change

```text
dimension
space identity
ordering
```

even when the continuous variable is still called $u$.

## 17. Solver choice is made after compilation

For the reduced product, B1 can select different reduced algorithms without changing
the compiler call.

The current application adapter dispatches between, for example,

```text
ReducedGradientSolverT
    steepest descent + selected globalization

ReducedLimitedMemoryBfgsSolverT
    L-BFGS direction + same reduced formulation/metric
```

Both receive

```text
reduced formulation
compiled metric
solver parameters
```

and solve from the same layout-compatible initial control.

That demonstrates the intended separation:

```text
compiler
    constructs the numerical problem

solver
    chooses how to optimize it
```

## 18. A complete reduced-use skeleton is short

Removing application-specific benchmarking, the essential shape is:

```cpp
auto scenario = chapter6::make_b1_scenario();

auto specification =
  chapter6::make_b1_problem_spec(scenario);

chapter6::dealii::B1SelectedDataT<2> selected_data(
  scenario.problem.forcing,
  scenario.problem.desired_state);

auto runtime =
  chapter6::dealii::make_b1_runtime_data(scenario, selected_data);

auto session =
  chapter6::dealii::make_b1_compilation_session<2>(scenario);

auto bindings =
  chapter6::dealii::make_b1_data_bindings(
    scenario.problem,
    scenario.problem.data.regularisation_weight,
    runtime);

auto policy =
  chapter6::dealii::make_b1_discretisation_policy(scenario.compile);

DealiiCompiler compiler;
auto result =
  compiler.compile(specification,
                   session,
                   bindings,
                   policy,
                   std::nullopt,
                   std::nullopt,
                   CompilationProduct::reduced_dto);

// inspect result, then consume result.problem
```

The full B1 execution adapter adds benchmark-specific validation, method dispatch,
evidence, and output handling around this core.

## 19. `CompiledProblemT` is useful beyond `make_reduced_dto()`

A successful reduced compilation exposes several related services:

```text
executable_model()
    residual/objective/JVP/VJP actions

metric()
    optimization geometry

constraint()
    optional projection constraint

state_adjoint_solvers()
    PDE solve services

reduced_hessian()
    optional compiled second-order action

compiled_application_view()
    compiler-path FE/application seam when available

manifest()
    typed record of realized decisions
```

Use the highest-level product that matches the task.

For ordinary reduced optimization, that is normally

```text
make_reduced_dto()
+
metric()
+
optional constraint()
```

The lower-level pieces are useful for diagnostics, specialized algorithms, or
application output, not as a requirement to reconstruct the formulation.

## 20. Other compilation products have different consumers

`CompilationProduct` currently includes

```text
reduced_dto

quadratic_kkt

pdas
```

The result members are product-specific.

Conceptually:

```text
CompilationProduct::reduced_dto
    -> CompiledProblemT

CompilationProduct::quadratic_kkt
    -> compiled quadratic KKT product

CompilationProduct::pdas
    -> compiled PDAS product
```

A supplied OTD formulation can also produce its own executable product when the
semantic declaration and target support it.

Do not compile one product and reinterpret it as another.

The choice should be made before compilation and checked by the capability layer.

## 21. Bounds are runtime realization data as well as semantic structure

A semantic box constraint declares the existence and role of lower/upper bounds.

The actual coefficient values are passed to compilation separately.

For a cellwise control:

```cpp
CellwiseBoxDataBindings bounds{lower, upper};

auto result =
  compiler.compile(specification,
                   session,
                   data,
                   policy,
                   bounds,
                   std::nullopt,
                   CompilationProduct::reduced_dto);
```

A facewise control uses `FacewiseBoxDataBindings`.

The two are not interchangeable because they belong to different realized control
layouts.

## 22. Using a ready-made problem effectively means changing the narrowest layer

Consider several common changes.

### Change the forcing function

Keep:

```text
recipe
ProblemSpec structure
compile policy
solver
```

and change the runtime function definition/binding.

### Change the desired state

Usually the same: change runtime data while keeping the graph.

### Change the regularization weight

For B1, this is runtime data. The graph still has the same regularization loss.

### Add a cellwise box

Change the recipe parameter so the graph declares the constraint, and supply concrete
bound data at compilation.

### Change from cellwise to continuous control

Change the recipe parameter.

This changes the semantic spaces/metric realization and therefore requires a new
compiled numerical problem.

### Change steepest descent to L-BFGS

Keep the compiled problem and select a different reduced solver, provided the required
capabilities are already present.

This layer-by-layer reasoning is the practical payoff of the architecture developed
in Parts I–III.

## 23. A scenario may intentionally freeze choices

A named scenario is not necessarily a general configuration object.

B1, for example, currently freezes

```text
assembled execution
reduced DTO product
selected solver families/globalization
specific application assumptions
```

That is appropriate for a benchmark or reproduction case.

If a user wants to explore a product outside those constraints, the right move may be
to use the underlying recipe and compiler directly rather than weakening the named
scenario until it no longer represents its original case.

A scenario should remain meaningful as a named application configuration.

## 24. Authoring a new compiled-path problem starts from mathematics

A new recipe should begin with the mathematical problem, not with an enum or compiler
branch.

Write down at least:

```text
state and decision variables

state equation / weak residual terms

boundary regions and conditions

observations

objective losses

control metric

constraints

formulation

analytical/discrete requirements
```

Then decide whether those concepts already exist in the semantic vocabulary from
Chapter 9.

If they do, a new problem family may require only a new composition.

If they do not, extending the semantic/compiler capability is a larger design change
than writing a recipe.

## 25. The recipe parameter type should expose meaningful family variations

Suppose a new scalar problem family has two genuine semantic variations:

```text
observation region

bounded versus unbounded control
```

A suitable parameter type could be:

```cpp
struct MyProblemParameters
{
  unsigned int observed_material_id = 1;
  bool with_cellwise_box = false;
};
```

The parameters should express choices meaningful to the problem family.

They should not include:

```text
CG tolerance

output directory

mesh filename

line-search constant
```

because those belong to later layers.

## 26. The recipe builder should compose a complete `ProblemSpec`

The recipe builder can use existing semantic factories, mutate a registered reference
spec deliberately, or construct the graph directly.

Its contract remains

```cpp
ProblemRecipeT<MyProblemParameters>{
  metadata,
  [](const MyProblemParameters &p) {
    ProblemSpec spec;
    // declare the complete semantic graph
    return spec;
  }
};
```

For a new family, the returned graph should stand on its own.

A downstream compiler should not need to know which C++ builder function produced it.

That is why semantic identity lives in the graph rather than in the recipe object's
type name.

## 27. Prefer composition over copying a ready-made spec blindly

Starting from a reference spec can be useful when the mathematical relationship is
clear.

For example:

```text
baseline scalar diffusion-reaction problem
    +
replace observation region
    +
add one registered requirement
```

can be a sensible derivation.

But copying a large `ProblemSpec` and changing string IDs until validation passes is a
poor authoring strategy.

Each modification should correspond to a mathematical statement.

The semantic graph is intended to preserve that relationship.

## 28. Runtime binders are a separate authoring task

Once the recipe exists, the application still needs a mapping from its runtime data to
the compiler bindings.

For each semantic `DataSpec`, decide:

```text
where does the concrete value come from?

what deal.II type realizes it?

who owns it?

what provenance should be recorded?
```

The B1 adapter's forcing and target functions are one example.

A general scalar target may additionally need tensor/vector coefficient bindings.

A boundary target may need fixed Dirichlet or observation-weight functions.

This is application/backend integration, not recipe construction.

## 29. Compile-option translation should be explicit

If the application has its own backend-neutral compile options, write one adapter that
translates them to the concrete compiler policy.

For B1:

```text
chapter6::CompileOptions
        ↓
make_b1_discretisation_policy
        ↓
DealiiDiscretisationPolicy
```

This is preferable to letting several call sites each translate state degree,
tolerances, and execution choices differently.

The adapter becomes the single boundary between application configuration and
compiler policy.

## 30. Solver-option translation should be equally explicit

The same principle applies after compilation.

The application may expose

```text
ReducedMethod
ReducedGlobalization
SolverOptions
```

while the actual solver constructors require

```text
ReducedSolverParameters
direction-policy objects
line-search-policy objects
```

The execution adapter owns that translation.

The semantic recipe should not.

## 31. The scenario is optional for reusable library code

A recipe is useful whenever a typed problem family exists.

A `ScenarioT` is useful when a named application case should keep together:

```text
problem selections
compile selections
solver selections
experiment metadata
```

Not every library user needs to create a scenario.

A small application can use:

```text
recipe
    ↓
ProblemSpec
    ↓
compiler
    ↓
solver
```

directly.

The scenario layer becomes more valuable when the repository needs discovery,
reproduction, benchmark identity, or several coordinated option groups.

## 32. The benchmark runner is orchestration, not another execution model

`HeadlessBenchmarkRunnerT` receives:

```text
scenario

problem builder

execution adapter
```

and performs

```text
build problem
    ↓
execute scenario
    ↓
collect measurements/evidence
    ↓
finalize artifact
```

It does not reproduce compilation or optimization logic.

The B1 execution adapter still performs the real compiler and solver calls.

This is useful when reading application code: the runner belongs to the experiment
shell, not to the numerical formulation.

## 33. After solving, keep the result at the right level

A reduced solver report contains optimization-level information such as

```text
final control
final retained evaluation
objective history
gradient-norm history
step history
solve counts
stopping reason
```

The final evaluation retains the state/adjoint information computed by the
formulation.

That means application output should normally reuse the retained result rather than
triggering a fresh PDE solve merely to reconstruct fields.

The compiler-path `CompiledApplicationViewT` exists for exactly this kind of
application-facing work.

## 34. `CompiledApplicationViewT` is the compiler path back to finite-element output

When present, the native view exposes information such as

```text
physical / independent state dimensions

physical / independent control dimensions

realized observation dimension

objective components

native output action
```

B1 uses it after optimization to write the retained final state, control, and adjoint.

This is intentionally separate from the generic optimizer interface.

The solver does not become responsible for VTK, reconstruction, or deal.II output.

## 35. The manifest answers "what did I actually compile?"

The semantic graph records the requested problem.

The compilation manifest records the realized discrete problem.

After a successful compile, retain

```text
compiled.manifest()
```

with the run/report when the exact realization matters.

The manifest can distinguish, for example:

```text
mesh identity

realized spaces and dimensions

binding provenance

selected formulation/product

metric realization

state/adjoint solve policies

constraint realization
```

This is particularly important when two runs share one semantic problem but differ in
mesh or numerical policy.

## 36. A useful decision tree for the compiled path

When approaching a new task, ask:

```text
Does a ready-made scenario already represent the run?
    yes -> start from the scenario

no:
Does a ready-made recipe represent the mathematical family?
    yes -> use the recipe directly and provide your own runtime/compile/solver setup

no:
Can the problem be expressed with existing semantic components
and registered compiler capabilities?
    yes -> author a new ProblemRecipeT and backend adapter

no:
the work is a semantic/compiler capability extension,
not merely application authoring
```

This prevents application code from becoming an accidental substitute for missing
framework capability.

## 37. The compiled path in one complete picture

```text
ready-made or new typed parameters
        ↓
ProblemRecipeT
        ↓
ProblemSpec
        │
        ├── runtime data definitions
        ├── mesh/session
        └── compile policy
        ↓
DealiiCompiler
        ↓
CompilationResultT
        │
        ├── diagnostics
        ├── CompiledProblemT / KKT / PDAS / OTD product
        └── manifest
        ↓
formulation product
        │
        ├── metric
        ├── optional constraint
        └── initial realized control
        ↓
solver
        ↓
report + retained final evaluation
        ↓
compiled application view / output / experiment shell
```

The semantic/compiler path is therefore usable at several levels.

A reader can consume a frozen ready-made scenario, reuse a recipe with different
runtime/solver choices, or author a new problem family without rewriting the
optimization algorithms.

## 38. Following the ready-made path through the source

The generic application types are:

- [`include/nmopt/application/recipe.hpp`](../../../include/nmopt/application/recipe.hpp)
- [`include/nmopt/application/scenario.hpp`](../../../include/nmopt/application/scenario.hpp)

The source Chapter 5 recipe family is:

- [`include/nmopt/application/chapter5.hpp`](../../../include/nmopt/application/chapter5.hpp)

The source Chapter 6 problem/scenario records and B1/B2 builders are:

- [`include/nmopt/application/chapter6.hpp`](../../../include/nmopt/application/chapter6.hpp)

The concrete deal.II B1 path is:

- [`include/nmopt/application/dealii/chapter6_b1.hpp`](../../../include/nmopt/application/dealii/chapter6_b1.hpp)

That file contains the complete sequence:

```text
selected data
runtime data
bindings
compilation session
compiler policy
compilation
reduced DTO
initial control
solver dispatch
native output
artifact/evidence envelope
```

The benchmark orchestration shell is:

- [`include/nmopt/application/runner.hpp`](../../../include/nmopt/application/runner.hpp)

For exact current programming interfaces, use:

- [Problem authoring](../../reference/problem-authoring.md);
- [Compiler](../../reference/compiler.md);
- [Application authoring](../../reference/application-authoring.md);
- [Optimization](../../reference/optimization.md).

For implementation mechanics, use the
[compiler implementation map](../../internals/compiler.md). Current supported
public compiler paths and rejection boundaries belong to the
[compiler reference](../../reference/compiler.md).

## Read next

The compiled path assumes that `nmopt` constructs the supported numerical
realization.

[Integrating an existing PDE application](13-integrating-an-existing-pde-application.md)
starts from the opposite situation: the mesh, matrices, solves, and output already
exist, and the task is to expose only the operations required by the same formulation
and solver layer.
