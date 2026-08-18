# Application API reference

This is the agent-facing reference for assembling and running a supported
application. It describes the public boundary between a typed application
description, the backend compiler, an optimization product, and experiment
provenance. The intended workflow is to select a recipe, fill its parameters,
select options, and call the documented interfaces; reading a lowerer should
not be necessary.

This document is an exact reference for the current public surface. The
mathematical and semantic contracts remain authoritative in the
[interface specification](../design/interface-specification.md), while the
implemented v1 capability and exclusion record is in the
[v1 semantic compiler record](../implementation/v1/semantic-compiler.md).
Concrete Chapter 5/6 recipes will be documented separately under
[`docs/applications/`](../applications/README.md).

## The assembly boundary

An application supplies values at these layers:

| Layer | Public boundary | Application responsibility |
| --- | --- | --- |
| Recipe | `application::ProblemRecipeT<Parameters>` | Map typed problem parameters to a semantic `ProblemSpec`. |
| Semantic graph | `semantic::v1::ProblemSpec` | Declare stable component IDs, composition, roles, spaces, and policies. |
| Runtime data | `compiler::v1::DealiiDataBindings<dim>` | Bind deal.II functions, scalar coefficients, and provenance after graph validation. |
| Discretization | `compiler::v1::DealiiDiscretisationPolicy` | Select finite-element degree, execution, linear-solve, and PDAS policies. |
| Product | `compiler::v1::CompilationProduct` | Select `reduced_dto`, `quadratic_kkt`, or `pdas`. |
| Executable service | `compiler::v1::CompiledProblemT<dealii_backend::SerialBackend>` | Consume the compiled reduced product and expose metric, constraint, and state/adjoint services. |
| Optimization | `solvers::Reduced*SolverT<Backend>` | Select solver family and stopping/line-search parameters. |
| Provenance | `experiment::ReducedExperimentEnvelopeT<Policy, Report>` | Detach the report, policy snapshot, compilation manifest, and environment from the service. |

The semantic graph is backend-neutral. Meshes, deal.II functions, solver
policies, and experiment metadata do not belong in `ProblemSpec`.

## Canonical workflow

The following sequence is the complete application path. A recipe or scenario
can package the inputs, but it does not change their ownership or ordering.

1. Select a recipe and construct its typed problem parameters.
2. Build the `ProblemSpec`. Every reusable semantic component needs a stable,
   unique ID; references use those IDs rather than display labels.
3. Create the deal.II data bindings and record their source provenance.
4. Create an owned `DealiiCompilationSession<dim>` for the mesh.
5. Choose a `DealiiDiscretisationPolicy` and a `CompilationProduct`.
6. Call `DealiiCompiler::validate` before compilation when a fast capability
   check is useful.
7. Call the session overload of `DealiiCompiler::compile` and inspect both
   `CompilationResultT::succeeded()` and its diagnostics.
8. Retrieve the product selected in step 5, configure the matching solver,
   solve from a layout-compatible initial control, and retain the returned
   report together with the manifest and policy snapshot.

The owned-session path is the default application path because the compiled
service may retain the mesh and bound callbacks. The borrowed triangulation
overload is available for callers that can enforce the required immutable
lifetime themselves.

### Minimal compilation outline

The following is an integration outline. `dim`, `forcing`, `desired_state`,
and `initial_control` are application-owned objects whose concrete types come
from the selected recipe and backend setup.

```cpp
using Backend = nmopt::dealii_backend::SerialBackend;
namespace v1 = nmopt::compiler::v1;

const nmopt::semantic::v1::ProblemSpec specification = recipe(parameters);

v1::DealiiDiscretisationPolicy policy;
const auto product = v1::CompilationProduct::reduced_dto;

v1::DealiiCompiler compiler;
const auto validation = compiler.validate(specification, policy, product);
if (!validation.valid()) {
  // Report every validation.diagnostics() entry and select a supported
  // recipe/product combination; do not guess at an unsupported lowering.
}

auto triangulation =
  std::make_unique<dealii::Triangulation<dim>>();
// Generate or read the mesh before transferring ownership to the session.
// dealii::GridGenerator::hyper_cube(*triangulation, 0.0, 1.0);

auto session = std::make_shared<v1::DealiiCompilationSession<dim>>(
  std::move(triangulation), "mesh-file-or-generation-record");

v1::DealiiDataBindings<dim> data{
  forcing,
  desired_state,
  1.0,                         // constant diffusion, when required
  0.0,                         // reaction coefficient
  parameters.regularisation,  // control regularisation weight
  {"forcing-source", "desired-state-source", "fixed-data-source"}
};

const auto result = compiler.compile(
  specification,
  session,
  data,
  policy,
  std::nullopt,
  std::nullopt,
  product);

if (!result.succeeded()) {
  // Inspect result.diagnostics().diagnostics() and stop.
}

const auto &compiled = *result.problem;
const auto reduced = compiled.make_reduced_dto();
```

The `DealiiDataBindings` constructor accepts, in order, forcing and desired
state functions, an optional constant diffusion value, reaction and
regularization scalars, and `DealiiBindingProvenance`. Optional ports are
then available for fixed Dirichlet data, general scalar data, weighted trace
data, and conservative transport data. Supply an optional port only when the
semantic graph declares the corresponding datum.

For a scalar reference problem, the existing semantic factories such as
`semantic::v1::make_scalar_diffusion_reaction_problem` can provide the graph
while recipe work is being completed. A final application recipe should hide
that construction behind `ProblemRecipeT` so an agent selects parameters
rather than reconstructing graph nodes by hand.

## Recipes, scenarios, and discovery

`ProblemRecipeT<Parameters>` is a typed builder. It owns only metadata and a
callable that returns a `ProblemSpec`; it does not own a mesh, compiler,
solver, or run output.

```cpp
struct ProblemParameters {
  double regularisation = 1.0e-2;
  bool bounded_control = false;
};

nmopt::application::ProblemRecipeT<ProblemParameters> recipe{
  {"scalar-diffusion-reaction",
   "Scalar diffusion-reaction control",
   "Chapter 5 scalar elliptic control recipe",
   "chapter-5",
   {"assembled", "serial", "scalar"}},
  [](const ProblemParameters &parameters) {
    return nmopt::semantic::v1::make_scalar_diffusion_reaction_problem(
      parameters.bounded_control);
  }
};

const auto specification = recipe(problem_parameters);
```

`ScenarioT<ProblemParameters, CompileOptions, SolverOptions,
ExperimentOptions>` binds one recipe ID to typed choices for all four option
groups. It validates metadata but intentionally does not compile or execute.
`ApplicationCatalog` stores metadata-only `RecipeMetadata` and
`ScenarioMetadata` entries for discovery and rejects duplicate IDs. Discovery
should therefore identify the recipe and scenario before any backend object is
constructed.

## Compilation options and products

### Discretization policy

`DealiiDiscretisationPolicy` is the compiler-facing option group. Its public
members are:

| Option | Default | Meaning |
| --- | --- | --- |
| `state_degree` | `1` | Finite-element degree used for the state and compiler-selected compatible spaces. |
| `execution` | `Execution::assembled` | Select assembled or matrix-free execution where the registered product supports it. |
| `control_metric_solve` | default `MassMetricSolveParameters` | Linear-solve policy for the control metric. |
| `state_solve` | default `SPDLinearSolvePolicy` | Linear-solve policy for state equations. |
| `adjoint_solve` | default `SPDLinearSolvePolicy` | Linear-solve policy for adjoint equations. |
| `pdas` | registered default `PDASPolicy` | Active-set policy for a PDAS product. |
| `pdas_kkt_solver` | default policy with v1 maximum iterations set to `200` | Equality-constrained KKT solve policy used by PDAS. |

Default construction is the intended starting point. Override a member only
when the recipe or benchmark contract names that option. A policy value is
not evidence that a product is supported: call `validate` and honor its
diagnostics.

### Product selection

| `CompilationProduct` | Result member | Consumer |
| --- | --- | --- |
| `reduced_dto` | `result.problem` | `CompiledProblemT`, reduced-space solvers, and trust-region methods. |
| `quadratic_kkt` | `result.kkt_problem` | Equality-constrained quadratic KKT contract. |
| `pdas` | `result.pdas_problem` | Box-constrained PDAS contract; requires the registered metric and bounds realization. |

The result also has `supplied_otd_problem` for a separately executable
supplied-OTD formulation when the selected semantic formulation and compiler
path produce that product. Supplied OTD is not a reduced DTO and must be
consumed through its own system contract. Never cast between result members.

Bounds are passed independently of `ProblemSpec`:

```cpp
std::optional<v1::CellwiseBoxDataBindings> bounds =
  v1::CellwiseBoxDataBindings{lower_value, upper_value};

const auto result = compiler.compile(
  specification, session, data, policy, bounds, std::nullopt, product);
```

Use `FacewiseBoxDataBindings` for a facewise boundary-control realization.
Cellwise and facewise bounds are distinct because they have different
discrete layouts and cannot be substituted for one another.

### Compiled reduced service

When `result.problem` is present, `CompiledProblemT` exposes:

- `executable_model()` for residual, objective, and derivative actions;
- `metric()` for the reduced covector-to-primal identification;
- `constraint()` when the compiled product has a compatible projection
  constraint;
- `box_data()` when bounds were compiled;
- `reduced_hessian()` when the selected product supplies one;
- `state_adjoint_solvers()` for the compiled state and adjoint services;
- `make_reduced_dto()` for the backend-neutral reduced optimization port; and
- `manifest()` for compilation and discretization provenance.

The returned reduced DTO retains the compiled lifetime owner. Keep the DTO or
compiled product alive for as long as a solver may call its executable model.

## Reduced solver options

The standard reduced-space entry point is
`solvers::ReducedGradientSolverT<Backend>`. It is the steepest-descent plus
Armijo specialization of `ReducedSearchSolverT`.

```cpp
nmopt::solvers::ReducedSolverParameters solver_parameters;
solver_parameters.maximum_iterations = 100;
solver_parameters.gradient_tolerance = 1.0e-8;

nmopt::solvers::ReducedGradientSolverT<Backend> solver(
  reduced, compiled.metric(), solver_parameters);
const auto report = solver.solve(initial_control);
```

For a compiled projection constraint, use the constrained constructor and
pass `*compiled.constraint()`:

```cpp
nmopt::contract::require(compiled.constraint() != nullptr,
                         "this scenario requires a compiled constraint");
nmopt::solvers::ReducedGradientSolverT<Backend> solver(
  reduced,
  compiled.metric(),
  *compiled.constraint(),
  solver_parameters);
```

`ReducedSolverParameters` has these defaults:

| Option | Default | Meaning |
| --- | --- | --- |
| `maximum_iterations` | `100` | Outer optimization iteration limit. |
| `maximum_line_search_trials` | `20` | Trial limit per line search. |
| `gradient_tolerance` | `1e-8` | Primary gradient stopping tolerance. |
| `stopping_criterion` | `automatic` | Built-in stopping policy selection. |
| `relative_gradient_tolerance` | `0` | Optional relative gradient criterion; zero disables it. |
| `objective_change_tolerance` | `0` | Optional objective-change criterion; zero disables it. |
| `step_tolerance` | `0` | Optional step criterion; zero disables it. |
| `initial_step_length` | `1` | Initial line-search step. |
| `armijo_fraction` | `1e-4` | Sufficient-decrease fraction for Armijo. |
| `backtracking_factor` | `0.5` | Multiplicative reduction after a rejected trial. |

Other direction and line-search policies are selected through the
`ReducedSearchSolverT` template parameters and constructor overloads. A
recipe should expose those choices as typed solver options only when its
contract requires them; it should not encode solver policy into the semantic
graph.

The initial control must have a layout compatible with the compiled metric.
For a constrained solve it must also be feasible. These are runtime contract
checks, not compiler option defaults.

## Diagnostics and provenance

Both `DealiiCompiler::validate` and `CompilationResultT::diagnostics` expose a
`semantic::v1::ValidationReport`. Check `valid()` first, then inspect every
`Diagnostic` in `diagnostics()` when it is false. Each diagnostic records:

- a category, such as structural, lowerability, or formulation capability;
- the affected component or formulation ID;
- the capability or policy that failed; and
- a remedy describing the supported correction.

An agent should use the remedy to change the recipe, binding, policy, or
product selection. It should not silently replace a requested component,
drop a bound, or switch formulations.

`CompiledProblemT::manifest()` is the authoritative record of what was
actually lowered. It includes the semantic problem identity, mesh lifetime and
provenance, data-binding provenance, realized formulation, solve policies,
and any constraint/product records. It is the input to an experiment
envelope, not a substitute for the solver report.

For reduced search reports, create a policy snapshot with
`experiment::make_reduced_search_policy_snapshot(report)`. Then construct a
`ReducedSearchExperimentEnvelopeT<Backend>` with the compilation manifest,
snapshot, report, and a caller-supplied `RunEnvironmentRecord`:

```cpp
const auto solver_policy =
  nmopt::experiment::make_reduced_search_policy_snapshot(report);

nmopt::experiment::RunEnvironmentRecord environment{
  source_revision,
  build_profile,
  compiler_name,
  compiler_version,
  standard_library,
  operating_system,
  architecture,
  hardware
};

nmopt::experiment::ReducedSearchExperimentEnvelopeT<Backend> envelope{
  compiled.manifest(),
  solver_policy,
  report,
  environment
};
```

The envelope owns values only and does not retain the compiled executable
service. This permits a run record to be serialized or archived after the
solver lifetime ends.

For Chapter 6 benchmark runs, wrap the detached envelope with
`application::benchmark::BenchmarkHarnessT<Scenario>`. The harness projects
scenario metadata into a deterministic `BenchmarkIdentity` and
`finalize(...)` creates a `BenchmarkArtifactT<Envelope>` containing:

- the scenario and recipe IDs, source reference/revision, build profile, and
  artifact directory;
- the complete compiler `ValidationReport`;
- the detached experiment envelope;
- optional wall-clock, CPU, and peak-memory measurements; and
- the explicitly selected output fields.

The harness validates identity and measurement shape but does not execute a
solver, lower a PDE, or serialize files. Those operations belong to the
headless runner and artifact writer. This keeps B0 from becoming a second
compiler or optimizer.

`application::benchmark::BenchmarkArtifactWriter` is the deterministic text
writer for the captured boundary. It emits canonical `key=value` lines with
stable ordering, escaped values, diagnostics, measurements, selected fields,
and caller-supplied artifact fields. It writes to a caller-owned stream; path
selection and directory creation remain orchestration concerns.

`application::benchmark::HeadlessBenchmarkRunnerT<Scenario>` provides that
orchestration boundary. Its `run(build_problem, execute)` call performs four
steps in order:

1. invoke `build_problem(scenario.problem)` to obtain the public `ProblemSpec`;
2. invoke `execute(specification, scenario)` to obtain
   `BenchmarkExecutionEvidenceT<Envelope>`;
3. capture runner wall time when the scenario requests timing measurements and
   finalize the evidence through `BenchmarkHarnessT`; and
4. render the detached artifact with `BenchmarkArtifactWriter`.

The execution callback owns backend compilation, solver invocation, and
construction of the detached envelope. It also supplies validation
diagnostics, non-wall-clock measurements, selected fields, and any additional
artifact fields. The runner does not choose a filesystem path, create a
directory, lower a PDE, or implement an optimization algorithm. An agent can
therefore assemble the typed scenario and provide the two application-specific
callbacks without reading the harness implementation.

The selected B1 deal.II integration is a separate backend-specific header,
`include/nmopt/application/dealii/chapter6_b1.hpp`. It supplies the
application-owned target function, mesh-session factory, runtime data-binding
factory, and reduced-solver execution adapter. It is intentionally not
included by the backend-neutral `application.hpp` umbrella header.

The B2 semantic and deal.II boundaries are assembled with
`chapter6::make_b2_scenario(GraetzCase)` and
`chapter6::make_b2_problem_spec(scenario)`. The four values in
`chapter6::b2_case_order` cover wings/full observation and constant/parabolic
targets; `chapter6::make_catalog()` registers each case with a unique stable
scenario ID. The helper adds the declared `fixed_dirichlet_data` lifting port
and `fixed_dirichlet_reconstruction` transformation to the Chapter 5
convection recipe. The deal.II adapter in
`include/nmopt/application/dealii/chapter6_b2.hpp` binds that function, the
conservative transport field, and the remaining scalar runtime data; it must
not patch a homogeneous graph at execution time. Its owned session realizes
the rectangle, boundary labels, and material-ID observation region, and its
execution adapter dispatches the selected full-BFGS reduced run.

## Headless executable boundary

The repository's `nmopt_runner` executable is the headless orchestration
boundary for benchmark runs. It owns command-line selection, the output root,
artifact-directory creation, and file writing; it does not own a second
compiler or optimization loop. The current boundary exposes `--list` for
metadata discovery, `--output DIRECTORY` for the runner-owned artifact root,
and `--run-kind` for selecting the run policy. `reproduction` is the default
and requires the `release-dealii` build profile. `--refinement` is an optional
benchmark mesh override; when it is absent, the selected benchmark supplies
its mesh default. A benchmark may instead select an adaptive strategy or
another mesh definition. The artifact retains the resolved mesh counts and
structural identity, so the run is not classified by an assumed universal
refinement number. The run kind does not yet select a separate generated-
output layout; that is the next runner-organization unit. Benchmark execution
commands are registered by their benchmark integration units and must consume
the typed scenario and execution-adapter interfaces described above.

## Agent checklist

Before implementation, an agent should be able to answer these questions from
the recipe and scenario documents alone:

1. What recipe ID and typed problem parameters produce the `ProblemSpec`?
2. Which semantic component IDs require runtime data bindings, and what are
   the provenance strings for those data sources?
3. Which dimension, mesh source, and mesh lifetime policy are required?
4. Which compile policy members and product are selected?
5. Are cellwise or facewise bounds required, and what layout do they have?
6. Which solver family, parameters, and initial-control feasibility rule apply?
7. Which manifest, report, policy snapshot, and environment fields must be
   retained as run evidence?

If one of these answers requires reading a lowerer, the corresponding recipe,
scenario, or reference documentation is incomplete and should be extended
before implementing the application.
