# External deal.II tutorial integration roadmap

## Status and authority

**Superseded for current implementation work.** Follow the
[external deal.II boundary evaluation roadmap](external-dealii-boundary-evaluation.md)
for the current protocol, work units, gates, and phase status.

This document preserves the historical T0–T5 proposal following the completed
[PDE–solver boundary refactor](review/pde-solver-refactor/roadmap.md).
The intent, selection rules, distributed-control/mass-metric path, work
sequence, and acceptance matrix below are historical planning, not current
implementation instructions. Their preservation does not assert that their
acceptance criteria were completed. The refactor roadmap remains complete.

The proposal's authoritative inputs were:

- [PDE–solver and application boundary](../design/pde-solver-boundary.md);
- [external deal.II integration reference](../reference/external-dealii-solver-integration.md);
- [A2 external-application proof](review/pde-solver-refactor/roadmap.md#a2--prove-an-unchanged-dealii-application-can-use-nmopt);
- the selected upstream deal.II tutorial source and its license/provenance
  record, once T0 is complete.

## Intent

Download one upstream deal.II tutorial step verbatim, preserve its standalone
forward behavior, and add the smallest project-owned binding needed to run a
real reduced optimization through nmopt.

The result should demonstrate that an application which was not designed
around nmopt can retain ownership of:

- the mesh, finite elements, DoF handlers, constraints, and assembled
  operators;
- state and adjoint solves and their diagnostics;
- control realization, objective data, and physical reconstruction; and
- native deal.II output and application lifetime.

nmopt should receive only the existing public solver-facing contracts:
`CallbackExecutableModelT`, `StateAdjointSolversT`, a metric and any selected
capabilities, `ReducedDTOT`, and a reduced optimizer.

## Initial selection rule

Start with a fixed-mesh Laplace/Poisson tutorial step, with deal.II step-4 as
the initial candidate. Fixed dimensions keep the reduced layouts stable while
the first external tutorial integration is being validated. Adaptive
refinement, nonlinear PDEs, distributed execution, and a second tutorial
family are follow-on work rather than prerequisites.

The final source selection is a T0 decision. It must be recorded before code
is copied into the repository.

## Source and provenance policy

The upstream tutorial remains the numerical reference. The imported source
must retain its upstream attribution and licensing information. T0 records:

- the official source URL;
- the exact deal.II release, tag, or commit;
- retrieval date and a content hash where practical;
- the applicable license/attribution notice; and
- every local patch or binding seam introduced after import.

The upstream baseline and the nmopt binding should remain visibly separate.
If the tutorial's original executable structure prevents reuse, introduce a
small, explicit seam or façade and record it as a local patch. Do not silently
rewrite the tutorial into an nmopt-specific application or hide the numerical
diff inside a generic adapter.

## Target layout

The sequence should produce two related executables:

1. an independently runnable upstream tutorial target, preserving the
   tutorial's forward solve and native output; and
2. an nmopt integration target that uses the same tutorial-owned numerical
   realization, adds the application-owned control/objective binding, and
   runs a small reduced optimization.

The first target proves that the imported source remains meaningful. The
second proves that the public external boundary is sufficient for a real
application rather than only for a purpose-built test fixture.

## Work sequence

### T0 — Select and pin the upstream tutorial

Choose the fixed-mesh tutorial step and record its source URL, version,
license, attribution, retrieval date, and source hash. Import the source
without nmopt dependencies and establish the standalone build/run command.

**Done when:** the source is identifiable and attributable, the standalone
forward target builds, and its output/residual behavior is recorded as the
baseline. No binding design should proceed against an unpinned source.

Prospective commit:

```text
docs(dealii): pin external tutorial source
```

### T1 — Preserve the standalone tutorial boundary

Package the imported tutorial so its original forward path remains directly
runnable. Keep the upstream numerical code and the project-owned integration
seam separate. If a minimal access seam is required, expose only the
application-owned assembled data and solve operations needed by the binding;
do not introduce an nmopt base class or compiler dependency into the tutorial
core.

**Done when:** the original tutorial target still runs independently, the
local patch list is complete, and the application still owns mesh, assembly,
solves, reconstruction, and native output.

Prospective commit:

```text
test(dealii): preserve upstream tutorial application
```

### T2 — Add the application-owned control and objective binding

Add the smallest tutorial-specific binding for a fixed-dimensional distributed
control and tracking objective. The binding should expose:

- residual, residual JVP, and residual VJP;
- objective and objective derivative;
- state and adjoint solves through `StateAdjointSolversT`;
- a compatible application-owned mass metric; and
- physical reconstruction and native output.

Control realization, objective assembly, and adjoint signs remain local to the
tutorial application. The generic nmopt contracts must not gain tutorial or
PDE-family knowledge.

**Done when:** the binding has one readable construction path, preserves the
tutorial's native state representation, and has no compiler or `ProblemSpec`
prerequisite.

Prospective commit:

```text
test(dealii): bind tutorial application to nmopt
```

### T3 — Compose the reduced optimization path

Construct the public integration boundary in a separate driver or adapter:

```text
tutorial application
    -> CallbackExecutableModelT
    -> StateControlPartitionT + StateAdjointSolversT + metric
    -> ReducedDTOT
    -> ReducedGradientSolverT
```

Keep all captured tutorial objects alive for the full lifetime of the DTO and
optimizer. Run a small optimization and write output through the tutorial's
native output path.

**Done when:** the integration target compiles without `CompiledProblemT`,
reaches the selected stopping criterion, reduces the objective, and emits
application-owned output.

Prospective commit:

```text
test(dealii): run reduced optimization from tutorial application
```

### T4 — Close the mathematical and ownership proof

Add focused checks for the imported application and its binding:

- standalone state residual at the solved state;
- residual JVP centered finite difference;
- residual JVP/VJP pairing;
- reduced derivative finite difference and Taylor remainder;
- state and adjoint solve reports;
- objective reduction and optimizer stopping evidence;
- unchanged objective/PDE callback separation; and
- native output written by the tutorial application.

Run the focused tutorial target and the applicable `debug-dealii` pipeline.
Run the neutral or sanitizer pipeline when shared backend-neutral code is
changed. A release run is an additional validation gate, not a reason to
weaken the focused mathematical checks.

**Done when:** the tutorial-specific proof passes, the upstream standalone
target remains green, and the ownership/lifetime boundary is documented from
actual code rather than inferred from snippets.

Prospective commit:

```text
test(dealii): verify external tutorial integration
```

### T5 — Publish the tutorial handoff

Update the external integration reference with the exact tutorial source,
build/run commands, binding seam, local patches, and verification results.
Add a concise tutorial entry or README next to the example. Keep the
reference's generic API guidance separate from tutorial-specific numerical
choices.

**Done when:** a new user can identify the upstream source, run it unchanged,
run the nmopt-integrated target, understand the ownership split, and reproduce
the derivative and optimization checks.

Prospective commit:

```text
docs(dealii): document external tutorial integration
```

## Boundaries and explicit non-goals

- Do not reopen the completed PDE–solver refactor roadmap.
- Do not add a tutorial framework, PDE-family inheritance hierarchy, or
  universal application interface.
- Do not make `CompiledProblemT`, semantic compiler recipes, or manifest
  provenance prerequisites for the external tutorial path.
- Do not move mesh, assembly, solve, objective, reconstruction, or output
  ownership into nmopt.
- Do not begin with adaptive mesh changes, MPI, nonlinear optimization, or a
  second tutorial family.
- Do not claim out-of-tree package/install support unless a separate packaging
  task validates the exported consumer surface.

## Acceptance matrix

| Boundary | Evidence required |
| --- | --- |
| Upstream fidelity | Pinned source, attribution, local patch list, and independently runnable forward target |
| Application ownership | Mesh, assembly, solves, reconstruction, and native output remain in the tutorial application |
| Public API sufficiency | Integration uses callback/model, solve, metric, DTO, and optimizer contracts only |
| Mathematical correctness | Residual, JVP/VJP, reduced derivative, Taylor, solve, and optimization checks pass |
| Lifetime correctness | Captured application state outlives DTO and optimizer operations |
| Documentation | Exact source provenance, commands, seam, and verification are published |

## Historical handoff and supersession

As checked on 2026-09-09, the previous tutorial branch
`codex/external-dealii-tutorial` at
`b4dce25ab03d893df48ab1d9eee2e4dc7cf9183d` already contains the
[pinned upstream Step-4 source](../../apps/external-dealii/step-4/upstream/step-4.cc),
an [adapted copy](../../apps/external-dealii/step-4/step-4.cc), and the native
wrapper/smoke attempt (`tutorial_application.hpp`, `tutorial_application.cc`,
and `tutorial_binding_smoke.cc` under `apps/external-dealii/step-4/`).
The earlier statement that T0 source selection had not started is stale.

That snapshot does not contain the controlled native-versus-current-nmopt
evaluation now planned. The existing Poisson contract fixture exercises API
functionality; it does not settle authentic adaptation cost or ergonomic
sufficiency. Current work and status are owned by the
[evaluation roadmap](external-dealii-boundary-evaluation.md), which preserves
both the completed refactor and the previous tutorial attempt as baselines.
