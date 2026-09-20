# Experiments and replication

The reusable numerical core answers questions such as:

- how do I represent a derivative?
- how do I combine state and adjoint solves into a reduced objective?
- how does a line search obtain trial values?

The project also needs to answer a different class of questions:

- which published or manufactured problem are we running?
- which mesh, target, control discretization, and solver policy were selected?
- how can the run be reproduced?
- which outputs and manifests belong to that run?

Those concerns belong to the experiment and replication layer.

The repository's recurring **Chapter 5** and **Chapter 6** labels refer to the
corresponding chapters of *Optimal Control of Partial Differential Equations* by
Manzoni, Quarteroni, and Salsa. **B1** and **B2** are repository scenario IDs for
two selected Chapter 6 numerical examples: distributed Laplace control and
Graetz-flow boundary control, respectively.

This layer is important to the repository, but it sits **around** the numerical
library rather than defining the library's core architecture.

## 1. Why an experiment and replication layer exists

A numerical optimization result is difficult to interpret if the run does not record
the choices that produced it.

For example, two **B1** runs could differ in:

- mesh topology;
- control discretization;
- forcing;
- desired state;
- regularization parameter;
- direction policy;
- stopping tolerance.

The experiment layer makes these choices explicit and associates them with solver and
compiler evidence.

The intended result is not only

```text
final control vector
```

but something closer to

```text
problem identity
+ declared numerical choices
+ effective compiler choices
+ solver policy and report
+ output artifacts
+ environment/provenance
```

## 2. Recipe, scenario, benchmark, and parameter file are different layers

These names are easy to blur because they often appear in the same execution path.

| Concept | Owns | Typical reuse |
| --- | --- | --- |
| Recipe | A typed builder for one reusable semantic problem family | Shared by multiple scenarios |
| Scenario | One concrete problem/data/discretization/solver selection | One named application case or a family of benchmark runs |
| Benchmark contract | What should be executed, recorded, compared, and accepted | One scientific/reproduction evaluation |
| Parameter file | A serialized way to choose supported typed scenario/run options | Convenient run configuration; not a new mathematical API |

A recipe should not own meshes, solvers, or files. A scenario should not implement a
second compiler. A benchmark should not redefine the PDE because a result is hard to
reproduce. A parameter file should not make a semantically unsupported combination
possible merely because it can express a string.

## 3. The main execution layers

The Chapter 6 path can be understood as six layers:

```text
scenario
   │
   ▼
problem recipe
   │
   ▼
semantic problem + runtime data
   │
   ▼
compilation session
   │
   ▼
execution adapter
   │
   ▼
headless runner
   │
   ▼
detached evidence + output artifacts
```

Each layer answers a different question.

### Scenario – which experiment?

A scenario gives a stable identity to one numerical experiment and freezes the
choices that should be reproducible.

It can include:

- geometry selection;
- PDE and target variant;
- control realization;
- formulation product;
- optimization method;
- stopping and globalization policy.

### Recipe – which problem family?

A recipe constructs the semantic problem structure shared by related scenarios.

For example, current B1 scenarios reuse the same distributed scalar-control family
while changing selected data or solver method.

### Runtime data – which concrete functions and coefficients?

Runtime bindings provide forcing, desired states, transport fields, fixed boundary
data, coefficients, and provenance.

### Compilation session – which owned numerical context?

The session owns or coordinates the deal.II realization needed to compile the
selected problem.

### Execution adapter – how is this compiled product run?

The adapter connects the compiled numerical product to the selected algorithm or
formulation solver and translates the result into the application-level evidence
model.

### Headless runner – how is the run orchestrated?

The runner provides the generic project workflow around problem construction,
execution, and detached reporting.

## 4. B1 as a complete example

The B1 scenario represents distributed Laplace-type control.

At the application level, a run can choose, among other things:

- a mesh realization;
- a distributed control discretization;
- a forcing definition;
- a desired-state definition;
- regularization;
- steepest descent or L-BFGS.

The scenario and recipe build the semantic problem. Runtime data supplies the
selected functions. A compilation session realizes the problem with deal.II. The
execution adapter obtains the compiled reduced formulation and runs the chosen
solver.

The important architectural point is that the runner is not assembling the PDE
itself.

```text
runner
  │
  ├─ asks application code to build problem/scenario
  ├─ asks adapter to compile and execute
  └─ records detached result/evidence
```

The compiler and optimizer remain reusable below that orchestration layer.

## 5. B2 shows why provenance must be explicit

The B2 Graetz-flow boundary-control case has more choices whose meaning matters to a
replication attempt.

Examples include:

- the observation region;
- the desired-state profile;
- the boundary-control representation;
- the transport boundary convention;
- quadrature/evaluation policy;
- mesh topology candidate;
- Armijo versus fixed-step globalization.

Several of these choices are not interchangeable implementation details. They can
define a different numerical problem.

The execution layer therefore records them rather than relying on “whatever default
the library currently uses.”

This is especially important when the source literature omits details. A development
candidate can be recorded as a hypothesis without being presented as recovered
source truth.

## 6. Benchmark coverage is not framework capability

The current executable benchmark programme is intentionally narrower than the
framework's implemented formulation surface.

B1 and B2 received the most complete end-to-end application work as the selected
Chapter 6 reduced-space demonstrations. Other formulation families do not have
equivalent named reproduction campaigns. That does **not** mean their support is
absent: selected KKT, supplied-OTD, complementarity, and PDAS products exist in the
compiler/contract layers with focused tests.

The experiment layer should therefore be read as a set of demonstrations and
scientific evaluations over the framework, not as the framework's capability
registry.

## 7. Parameter files are scenario inputs, not the architecture

Checked parameter files under `parameters/` provide convenient, reviewable run
configurations.

They are useful because they let a reader see a concrete experiment without
recompiling the application code for every parameter change.

But a parameter file does not replace the typed scenario or semantic model.

The flow is conceptually

```text
parameter file
     │
     ▼
typed application parameters
     │
     ▼
scenario / runtime selections
     │
     ▼
normal semantic + compiler + solver path
```

This keeps text parsing at the outer edge of the project.

## 8. Detached evidence keeps provenance out of numerical ownership

Compilation and solving produce evidence that should survive after the numerical
objects themselves are gone.

The project therefore uses detached records for things such as:

- compilation manifest;
- selected/effective policies;
- solver report;
- environment information.

These records can be written into artifacts without forcing the optimization solver
to own the mesh or compiler session.

This separation is useful for both software design and scientific reproducibility.

## 9. Output is owned by the native application side

The numerical solver may return a final control and a retained state, but writing a
meaningful finite-element result still requires application knowledge.

For example, B2 boundary control may need to be written differently depending on
whether it is facewise constant or a continuous trace field.

The execution adapter/native view has access to that information. The generic
optimizer does not.

The result path is therefore roughly

```text
solver result
   │
   ├─ final control
   └─ retained final state
          │
          ▼
native reconstruction / output
          │
          ▼
VTK or other application artifacts
```

## 10. Run sets and artifact directories

A research run often belongs to a family rather than standing alone.

The runner layer can organize:

- one run;
- a parameter sweep;
- method comparisons;
- benchmark/reproduction matrices.

Artifacts then give each run a stable location for its report, manifest, numerical
outputs, and post-processing inputs.

The important point is that run-directory structure is project infrastructure. It
does not leak into the reduced formulation or compiler contracts.

## 11. Post-processing is downstream of numerical execution

Plots and comparison summaries consume completed run artifacts.

They should not change the mathematical definition of the run.

This makes it possible to regenerate presentation outputs from recorded artifacts
without solving the PDE again, and it makes missing or changed post-processing code
less dangerous to the scientific meaning of the numerical result.

## 12. The current Chapter 6 application status

The current active application/parameter path is concentrated on B1 and B2.

B1 is the distributed scalar-control family used for the main reduced-method
comparisons.

B2 is the Graetz-flow boundary-control family with observation-region and boundary
realization choices.

The repository also contains later scenario contracts and formulation capabilities
for KKT/PDAS/all-at-once work, but those should not be read as if every scenario has
the same current reproduction status or runner coverage.

This is one reason the application, benchmark, and planning documentation remains
separate from the reusable numerical overviews.

## 13. How this layer relates to a library user

A user integrating their own application does **not** need to adopt:

- Chapter 6 scenario IDs;
- repository parameter-file schemas;
- benchmark run directories;
- project post-processing scripts.

They can use the numerical contracts and solvers directly.

A user who wants reproducible project-style experiments may reuse some of the
application/experiment machinery, but it is an outer convenience layer rather than a
required dependency of the reduced solver.

## 14. Source orientation

The relevant repository areas are:

```text
include/nmopt/application/
    generic orchestration plus project application/scenario support

include/nmopt/experiment/
    detached provenance/evidence records

apps/nmopt-runner/
    command-line project runner

parameters/
    checked B1/B2 parameter families

docs/studies/
    chapter-5/   source catalogue and application recipes
    chapter-6/   numerical methods, examples, scenarios, benchmarks, and replication
```

This mixed structure reflects the project's history: reusable orchestration support
and repository-specific application code currently live near each other.

A reader should distinguish them by responsibility rather than assume everything
under `application/` is required for library integration.

## 15. A practical reading path for reproducing a run

If your goal is to reproduce one of the repository's numerical experiments:

1. identify the scenario in the application documentation;
2. inspect the corresponding checked parameter file;
3. understand the problem recipe and runtime data it selects;
4. run through `nmopt_runner`;
5. inspect the solver report and compilation manifest;
6. inspect native numerical output;
7. run the post-processing/comparison tools if needed.

If your goal is instead to understand or reuse the underlying optimization
machinery, return to [Project architecture](project-architecture.md) and
[Reduced optimization](reduced-optimization.md); the runner is not the best entry
point for that task.

## 16. Why this layer is separate

The separation can be summarized as:

```text
numerical core:
    "How is this PDE optimization problem evaluated and solved?"

experiment and replication layer:
    "Exactly which problem/run did we execute, and what evidence belongs to it?"
```

Both are necessary for the repository's goals, but they should not be confused.

The numerical core aims at reusable mathematical software. The experiment and
replication layer selects, executes, records, and compares concrete uses of that
software. Benchmarks are therefore clients and evidence producers for the framework,
not the framework's architectural center.

## Read later: authoritative sources

Use the following documents after this overview according to the task:

- [Application execution](../../reference/application-execution.md) owns the exact run,
  artifact, manifest, and output contracts.
- [Parameter files](../../reference/parameter-files.md) owns the supported parameter
  schema and checked file conventions.
- [Application assembly API](../../reference/application-api.md) owns the construction
  interfaces between scenarios, semantic problems, runtime data, and compilation.
- [Chapter 5 application recipes](../../studies/chapter-5/recipes.md) and
  [Chapter 6 application scenarios](../../studies/chapter-6/scenarios.md) own the concrete
  application-level records.
- [Chapter 6 benchmark specification](../../studies/chapter-6/benchmarks.md) owns benchmark
  acceptance and evidence requirements.
- [Chapter 6 numerical-examples reference](../../studies/chapter-6/numerical-examples.md)
  owns what the source book actually states and what it omits.

Planning and review files may explain how those contracts were reached, but they
should not be used as the default source for current run semantics.
