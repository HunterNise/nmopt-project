# Representative runtime paths

**Audit baseline:** [`f53b7f009e5c`](https://github.com/HunterNise/nmopt-project/commit/f53b7f009e5c418ec4f3855db29f7eb924faf6f1)  
**Status:** complete for the architecture decision. The external reduced path and semantic/compiler reduced path are traced in detail; the research runner is traced to the level needed to classify its responsibilities. Non-reduced product centrality is established from compiler/build/test evidence rather than a second full end-to-end runtime trace.

## Why runtime traces matter

Static type and directory structure can suggest boundaries that disappear at runtime. This file follows actual construction and call sequences to answer:

- who owns state;
- which subsystem calls which;
- where type/representation adaptation happens;
- what infrastructure is optional;
- which layers own mathematical decisions versus numerical orchestration.

## Path A — external Step-4 Problem B to reduced optimization

### 1. Native application preparation

Primary source:

- [`apps/external-dealii/step-4/source/adapted/step-4.cc`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/apps/external-dealii/step-4/source/adapted/step-4.cc)

The adapted `Step4<2>` remains an ordinary deal.II application. It owns:

```text
Triangulation
FE_Q
DoFHandler
SparseMatrix
system RHS
native solution
boundary values
CG solve
VTK output
```

The adaptation adds reusable numerical seams:

```text
prepare_for_external_use()
system_matrix_view()
system_rhs_view()
dof_handler_view()
boundary_values_view()
solve(rhs, solution)
output_results(state, filename)
```

`prepare_for_external_use()` simply exposes the original preparation sequence:

```text
make_grid()
setup_system()
assemble_system()
```

No nmopt type appears in this file.

### 2. Application-owned optimal-control problem

Primary source:

- [`apps/external-dealii/step-4/integration/problem_b.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/apps/external-dealii/step-4/integration/problem_b.hpp)

`ProblemB` is also outside nmopt. It owns the OCP mathematics built on Step-4:

```text
free/full state coordinates
mass matrix and rectangular control coupling
residual
residual JVP
residual VJP
objective
objective derivative
state solve
adjoint solve
```

Its state/control dimensions are intentionally different:

```text
state     = independent/free state coordinates
control   = full FE coefficient vector
```

The residual is implemented in the expected form:

```text
state operator
- restricted base RHS
- control coupling
```

The state solve:

1. starts from Step-4's already boundary-treated RHS;
2. adds the embedded distributed-control load;
3. calls Step-4's native `solve`;
4. restricts the full solution to independent state coordinates.

The adjoint solve similarly embeds the state-objective covector into the full Step-4 system and reuses the native solve.

### 3. Application-owned metric

Primary source:

- [`apps/external-dealii/step-4/integration/problem_b_metric.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/apps/external-dealii/step-4/integration/problem_b_metric.hpp)

`ProblemBMetric` owns the FE mass geometry and its inverse solve:

```text
apply(v)          = M v
inverse_apply(r)  = CG solve M g = r
```

It records native convergence evidence.

This is important because the external application already has the numerical operation nmopt needs; the binding does not assemble or reinterpret the metric mathematics.

### 4. nmopt binding

Primary source:

- [`apps/external-dealii/step-4/minimal/problem_b_binding.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/apps/external-dealii/step-4/minimal/problem_b_binding.hpp)

The binding introduces the nmopt representation boundary.

Construction sequence:

```text
ProblemB
   |
   v
BlockLayout(variable = state + control)
BlockLayout(test = state-test)
   |
   v
CallbackExecutableModel
   |
   +-- residual callback -----------------> ProblemB::residual
   +-- residual JVP callback -------------> ProblemB::residual_jvp
   +-- residual VJP callback -------------> ProblemB::residual_vjp
   +-- objective callback ----------------> ProblemB::objective
   +-- objective derivative callback -----> ProblemB::objective_derivative
   |
   v
StateControlPartition
   |
   +-- state layout
   +-- control layout
   |
   +--> state solve callback -------------> ProblemB::solve_state
   +--> adjoint solve callback -----------> ProblemB::solve_adjoint
   |
   v
metric adapter ---------------------------> ProblemBMetric
   |
   v
ReducedDTOT
```

The binding's main responsibilities are therefore:

- declare runtime layout identity;
- wrap native vectors in primal/covector blocks;
- translate native solve evidence into `LinearSolveReport`;
- expose the native metric through `MetricT`;
- construct the reduced formulation service.

It does **not**:

- create a semantic `ProblemSpec`;
- invoke the compiler;
- create a compilation manifest;
- use application recipes/scenarios;
- use the project runner;
- take ownership of the Step-4 mesh or PDE matrices.

### 5. Reduced solver

Primary source:

- [`apps/external-dealii/step-4/minimal/problem_b.cc`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/apps/external-dealii/step-4/minimal/problem_b.cc)

The executable constructs:

```text
Step4 application
 -> prepare
ProblemB
ProblemBBinding
initial control
ReducedSearchSolver(binding.reduced(), binding.metric(), parameters)
 -> solve(initial)
```

The solver sees only:

- a reduced formulation service;
- a metric;
- solver parameters;
- a control primal block.

The application/compiler distinction has disappeared at this point.

### 6. One initial reduced evaluation

From `ReducedDTOT` and the Problem B callbacks, the actual call flow is:

```text
solver
  |
  v
ReducedDTOT::evaluate(u)
  |
  +--> evaluate_value(u)
  |      |
  |      +--> ProblemB::solve_state(u)
  |      |      |
  |      |      +--> Step4::solve(controlled RHS)
  |      |
  |      +--> ProblemB::objective(state, u)
  |
  +--> augment_derivative(value)
         |
         +--> ProblemB::objective_derivative(state, u)
         |
         +--> ProblemB::solve_adjoint(J_state)
         |      |
         |      +--> Step4::solve(adjoint RHS)
         |
         +--> ProblemB::residual_vjp(adjoint)
         |
         +--> reduced covector = J_control - E_control^* p
```

Notably:

- `ProblemB::residual()` is not called in this optimization sequence;
- `ProblemB::residual_jvp()` is not called in this optimization sequence;
- the full VJP returns state and control components, although `ReducedDTOT` ultimately extracts the control component after using the state derivative separately for the adjoint solve.

These operations remain valuable for verification and other formulations, but they are broader than the minimum first-order reduced runtime need.

### 7. Trial evaluations and accepted steps

The reduced search implementation separates trial value work from derivative work.

For a line-search trial:

```text
trial control
 -> ReducedDTOT::evaluate_value
 -> state solve
 -> objective
```

Only after acceptance:

```text
accepted retained value
 -> augment_derivative
 -> objective derivative
 -> adjoint solve
 -> VJP
```

This is a real formulation/algorithm service: rejected trials do not pay for adjoints merely because the application exposes them.

### 8. Final output returns to the application

The final result retains the accepted state in reduced/free coordinates.

`problem_b.cc` then performs:

```text
result.final_evaluation.state
 -> ProblemB coordinates reconstruct
 -> full physical FE state
 -> Step4::output_results
```

Output ownership therefore returns to the native application. The optimizer does not know how a finite-element field is reconstructed or written.

## Ownership summary

```text
Step4
  owns mesh, FE, DoFHandler, stiffness matrix, base RHS, native solve, VTK

ProblemB
  owns OCP coordinates, control coupling, mass-based objective,
  derivatives, state/adjoint RHS construction

ProblemBMetric
  owns native FE control metric and inverse solve

ProblemBBinding
  owns nmopt representation/adaptation objects
  borrows ProblemB

ReducedDTOT
  owns reduced state/adjoint orchestration
  borrows/copies the supplied model/solve services according to constructor

ReducedSearchSolver
  owns optimization orchestration and work/evidence records
  borrows ReducedDTOT and Metric

main
  owns lifetime/order of all of the above
```

## What this path proves about repository scope

The following subsystems are **not on the runtime path** for an existing application using reduced optimization:

```text
semantic/v1
compiler/v1
application recipes/scenarios
experiment envelope
nmopt-runner
Chapter 5/6 catalog
post-processing/report tooling
```

This does not make those subsystems unimportant. It establishes that they are not prerequisites for the core external numerical-library use case.

## Human-readability observations from this path

### Good

- ownership stays local and visible;
- application mathematics is not hidden inside nmopt;
- the optimizer boundary is genuinely independent of deal.II/compiler concepts;
- the final state can remain application-owned/output using native code;
- callback adaptation is explicit rather than inheritance-based.

### Friction

The binding must understand several nmopt representation concepts:

```text
BlockLayout
SpaceId
PrimalBlockT / CovectorBlockT
CallbackExecutableModelT
StateControlPartitionT
StateAdjointSolversT
FormulationSolveResultT / LinearSolveReport
MetricT subclass
ReducedDTOT
```

Some of this is mathematically meaningful; some is construction ceremony. Later `04-public-api-and-user-surfaces.md` will distinguish the two.

The metric adapter exposes another asymmetry: the native metric inverse returns solve evidence, but `MetricT::inverse_apply` returns only a primal block. The minimal adapter checks convergence and discards the detailed solve record. This is an interface characteristic to document, not yet a request to alter the API.

## Trace coverage decision

The detailed external trace establishes the minimal existing-application path.
Path B below establishes that semantic compilation rejoins the same reduced
formulation/optimizer spine. The Chapter 6 runner is inspected far enough to
separate research execution, parameter binding, provenance, and artifact
responsibilities from the reusable numerical core.

A second full end-to-end trace for KKT/PDAS/supplied-OTD was not necessary for
the documentation decision: dedicated contract tests, compiler integration
tests, and build targets establish that those products are implemented and
independent enough to be presented as secondary capabilities.

The dense/reference paths are verification infrastructure and likewise do not
change the user-facing architecture.

# Path B — semantic Chapter 6 B1 to the same reduced optimizer

Primary sources:

- [`include/nmopt/application/chapter6.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/application/chapter6.hpp)
- [`include/nmopt/application/dealii/chapter6_b1.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/application/dealii/chapter6_b1.hpp)
- [`include/nmopt/compiler/v1/dealii_compiler.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/compiler/v1/dealii_compiler.hpp)
- [`include/nmopt/compiler/v1/compiled_problem.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/compiler/v1/compiled_problem.hpp)

This path adds a substantial front end around the same downstream reduced services.

## 1. Scenario/configuration remains backend-neutral

`chapter6.hpp` defines typed B1 problem, compilation, solver, mesh, and experiment records without including deal.II.

The B1 scenario shape is:

```text
B1Scenario
    problem     B1ProblemParameters
    compile     CompileOptions
    solver      SolverOptions
    experiment  ExperimentOptions
```

The compile record includes mesh/discretization/solve-policy selections. The solver record includes reduced method/globalization/stopping settings. Experiment settings include source references, build profile, output identity, and retention policy.

This is application/research configuration, not the numerical formulation itself.

## 2. The deal.II execution adapter turns the scenario into concrete runtime data

`B1ReducedExecutionAdapterT` receives:

```text
regularisation value
runtime deal.II data functions
owned DealiiCompilationSession
run environment
optional native-output directory
```

On invocation it receives a `ProblemSpec` and the B1 scenario.

It then:

1. validates that the scenario and semantic problem match the expected B1 family;
2. converts application compile options into `DealiiDiscretisationPolicy`;
3. creates concrete deal.II data bindings;
4. calls `DealiiCompiler::compile(...)`.

The adapter is therefore a **client that bridges project-level scenario choices into the general compiler interface**.

## 3. Compilation produces a reduced product

The B1 path requests `CompilationProduct::reduced_dto`.

Conceptually:

```text
ProblemSpec
+ DealiiCompilationSession
+ concrete data bindings
+ discretisation policy
+ requested product
        |
        v
DealiiCompiler::compile
        |
        v
CompilationResultT
        |
        +-- diagnostics
        |
        +-- CompiledProblemT
```

`CompiledProblemT` packages:

```text
ExecutableModelT
MetricT
optional ConstraintT
StateAdjointSolversT
optional ReducedHessianT
CompilationManifest
optional NativeApplicationViewT
lifetime owner
```

The compiler-specific typed realization stays behind those ports.

## 4. The path rejoins the external runtime spine

The adapter immediately performs:

```text
compiled problem
    |
    +--> make_reduced_dto()
    |
    +--> metric()
    |
    +--> executable_model() -- used here to derive the control layout
    |
    v
ReducedGradientSolver / ReducedLimitedMemoryBfgsSolver
```

The optimizer receives the same categories of objects as the external Step-4 consumer:

```text
ReducedDTOT
MetricT
control PrimalBlock
solver parameters
```

At this point it does not know that its services came from `ProblemSpec` or `DealiiCompiler`.

## 5. Compiler-specific application information remains beside the solver boundary

After optimization, B1 asks:

```text
compiled_problem.native_application_view()
```

for:

- physical/independent state/control dimensions;
- native field output.

This is deliberately separate from `ExecutableModelT`.

The runtime therefore has two parallel outputs from compilation:

```text
typed numerical realization
     |\
     | +--> solver-facing mathematical services --> ReducedDTOT --> optimizer
     |
     +----> NativeApplicationView --> dimensions/native output
```

That is a concrete implementation of the principle that solver type erasure is not the ownership boundary for all numerical state.

## 6. Experiment/provenance packaging is post-solver infrastructure

After solving, B1 builds:

- Hessian verification evidence;
- selected artifact fields;
- solver-policy snapshot;
- compilation manifest;
- run-environment record.

These are combined into a `ReducedSearchExperimentEnvelope`.

The envelope intentionally stores **values only**, not executable services.

This is important for project scope:

> provenance/artifact infrastructure surrounds a numerical run; it is not part of the formulation or optimizer contract needed to perform that run.

## 7. Comparison with the external path

```text
EXTERNAL STEP-4                     SEMANTIC B1

native PDE/OCP operations           ProblemSpec
       |                                |
       |                                v
       |                          DealiiCompiler
       |                                |
       v                                v
callback model + solves          CompiledProblemT
       |                                |
       +---------------+----------------+
                       |
                       v
                  ReducedDTOT
                       |
                       v
                reduced optimizer
```

The two paths differ substantially **before** the formulation:

- external path: application already owns the numerical realization;
- semantic path: nmopt owns the authoring model and compiler realizes it.

They converge before optimization.

## 8. What compilation adds

Relative to the minimal external path, the semantic route adds legitimate responsibilities:

- declarative problem construction;
- semantic structural validation;
- explicit lowerability/capability checks;
- FE/discretization realization;
- lifetime-managed compiled numerical state;
- typed compilation provenance;
- optional native application/output view;
- project scenario/configuration bridging.

It does **not** add a second reduced optimization algorithm.

## 9. Human-facing implication

The project should probably be explained with the formulation/optimizer as a shared center and the semantic/compiler system as a producer path, rather than as:

```text
semantic -> compiler -> application -> solver
```

for every user.

The latter would be false for the external integration and would make the repository look more vertically coupled than the implementation actually is.

# Path C — Chapter 6 runner and research-execution wrapper

The executable `apps/nmopt-runner` adds another shell outside Path B.

Its source set includes:

```text
main.cc                  ~54 KB
parameter_files.hpp      ~53 KB
runner.hpp               ~25 KB
parameter_binding.hpp    ~13 KB
run_set_plan.hpp
benchmark registry/binders
```

The top of `main.cc` handles:

- command-line parsing;
- parameter-file selection;
- benchmark catalog selection;
- run-set coordinates;
- artifact path decisions;
- framework revision/build-profile provenance;
- Chapter-specific B1/B2 dispatch.

This strongly suggests that `nmopt-runner` should be treated as **research/application execution infrastructure**, not as the primary entry point to the numerical library.

The relevant separation is visible in the implemented layering:

```text
recipe/scenario
    |
execution adapter
    |
HeadlessBenchmarkRunnerT
    |
nmopt_runner CLI / parameter matrices / run-set planning
    |
artifacts, traces, manifests, native output, post-processing
```

The numerical solve is delegated through the execution adapter; the outer
runner owns experiment selection, parameterization, provenance, filesystem
layout, and artifact policy. A more exhaustive line-by-line runner trace would
be useful only for runner maintenance, not for deciding the project's
human-facing architecture.
