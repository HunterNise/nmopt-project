# Minimal Step-4 consumer assessment

Status: MC1–MC6 complete through `2ba749b`, reviewed and closed on 2026-09-12.
This report compares the two runnable consumers; the
[closure report](closure-report.md) assesses the full phase against its
original goals. G1 and G2 remain completed evaluations.

## Result

The current public boundary is sufficient to wire both existing Step-4
problems into the reduced optimizer without changing shared nmopt, the
compiler, or the deal.II backend. The consumers are explicit application
examples: each author must assemble the public layouts, executable callbacks,
state/control partition, native solve services, metric, reduced DTO, and
optimizer configuration.

The examples demonstrate functional wiring for these two linear, symmetric,
serial Step-4 cases. They are not claims of globally smallest integration,
production readiness, newcomer usability, or applicability to nonlinear,
nonsymmetric, constrained, distributed, compiler, or packaged applications.
No unfamiliar-user exercise was performed. Repeated construction is an
observed authoring limitation, not merely an absence of user testing.

Problem A and Problem B use the same public composition shape. B additionally
has genuine application-specific work: it borrows a prepared `Step4<2>`
application through `ProblemB`, uses free state and adjoint coordinates with
full control coordinates, reconstructs the physical state for output, and
adapts a nonidentity mass metric. These differences are not generic framework
ceremony.

## Ownership map

The following map separates reused numerical capability, functional consumer
wiring, and accessory evaluation support.

| Area | Problem A | Problem B | Role |
| --- | --- | --- | --- |
| Step-4 discretization and writer | `source/adapted/step-4.cc`; `integration/problem_a.hpp` owns its prepared `Step4<2>` application | `source/adapted/step-4.cc`; `minimal/problem_b.cc` prepares `Step4<2>` and `integration/problem_b.hpp` borrows it | Existing application capability |
| PDE/OCP operations | `integration/problem_a.hpp` owns residual, derivatives, objective, adjoint, and control pullback | `integration/problem_b.hpp`, `problem_b_coordinates.hpp`, `problem_b_mass.hpp`, and `problem_b_metric.hpp` own residual, derivatives, objective, coordinates, coupling, and mass operations | Existing application formulation |
| Public binding | `minimal/problem_a_binding.hpp` supplies layouts, five callbacks, solve wrappers, partition, and identity metric | `minimal/problem_b_binding.hpp` supplies the same public pieces plus B's mass metric adapter and coordinate-aware formulation calls | Required functional wiring |
| Entry point | `minimal/problem_a.cc` constructs `ProblemA`, binding, zero control, policy, solver, and output | `minimal/problem_b.cc` constructs/prepares the application, `ProblemB`, binding, zero control, policy, solver, reconstruction, and output | Required consumer behavior |
| Evaluated reference | Native A reduced evaluation/optimizer and `integration/nmopt_binding.hpp` | Native B reduced evaluation/optimizer and `integration/nmopt_problem_b_binding.hpp` | Accessory comparison support |
| Independent checks | `tests/application/external_step4_minimal_problem_a_contract.cc` | `tests/application/external_step4_minimal_problem_b_contract.cc` | Accessory verification |
| Evidence and diagnostics | `evaluation/`, `verification/`, `diagnostics/`, and ignored `runs/` artifacts | Same directories and artifact role | Accessory; not needed by the consumer |

`ProblemABinding` and `ProblemBBinding` are local containers for the borrowed
objects and callbacks. Their class names and convenience accessors are not
new public framework abstractions. After MC4, solve services are constructed
at the `ReducedDTO` boundary rather than retained as a redundant member.

The lifetime obligations are functional: A's `ProblemA` owns its prepared
Step-4 application and outlives its binding and solver; B's entry point keeps
the prepared application alive before `ProblemB`, its binding, metric, and
solver. The binding callbacks borrow the corresponding problem and must not
outlive it.

## Minimum functional wiring

Here “minimum” means the smallest complete path demonstrated by these
examples, not a proof that every external application needs exactly these
lines.

### Shared public-contract work

Each consumer author supplies:

- a variable layout containing `state` and `control` blocks and a residual
  test layout;
- all five `CallbackExecutableModelT` operations: residual, residual JVP,
  residual VJP, objective, and objective derivative;
- a `StateControlPartitionT` and state/adjoint solve callbacks;
- actual native iterative evidence translated into `LinearSolveReport` values;
- a `MetricT`, a `ReducedDTOT`, an initial control, and a solver invocation;
  and
- application-owned output if a physical field or checkpoint is required.

The current executable-model contract requires all five callbacks even though
the selected first-order reduced optimization does not call residual or JVP
on its successful path. State and adjoint callbacks must preserve native
convergence and failure behavior; they cannot label a failed native solve as
an exact solve.

### A-specific work

A uses the full 289-entry algebraic control and a local identity metric.
`ProblemA` provides the native operations and owns the prepared Step-4
application, so its retained solver state can be sent directly to the writer.

### B-specific work

B uses 225 free state/adjoint entries and 289 full control entries. It does
not use a free-coordinate control. `ProblemB` owns the free/full coordinate
maps, lifting, coupling, objective, derivatives, and native mass service. The
entry point reconstructs the full state before output.

The local `ProblemBMassMetric` adapts the native mass service to `MetricT`.
Its inverse checks the native CG result and returns the primal block required
by the public metric interface. The public metric interface does not return a
metric-solve report. State and adjoint solve wrappers do return public solve
reports because those reports are part of `FormulationSolveResultT`.

## Functional versus accessory structure

The following distinctions prevent the example's implementation choices from
being mistaken for universal nmopt obligations.

| Concern | Required for this consumer path | Local or evaluation-only choice |
| --- | --- | --- |
| Layouts, five callbacks, partition, solve services | Yes; the current public contracts require them | The `ProblemABinding`/`ProblemBBinding` class wrappers and their names are local packaging |
| State/adjoint solve reports | Yes; translate actual native CG evidence to `LinearSolveReport` | Optional counters, serialized solve records, and native optimization traces |
| Metric | Yes; A uses identity and B uses its mass metric | The metric class names and retained-object arrangement are local implementation choices |
| `ReducedDTO` and optimizer policy | The DTO and a complete policy are needed to run the selected example | Zero initialization, frozen tolerances, stopping choice, CLI parsing, and report formatting are example policy |
| Output | Application-owned if the consumer must write a field | VTK comparison, timestamp stripping, output logs, and ignored run allocation are test support |
| Diagnostics | No instrumentation is instantiated by either consumer | Diagnostic headers may remain transitive dependencies of reused application headers; purging them is a separate cleanup |

The focused tests contain much more code than the consumers. They create
separate native and public instances, compare vectors and histories, recompute
fresh state/adjoint/covector values, use dense A/B oracles, audit final
stationarity, validate executable reports and VTK payloads, and write ignored
evidence. None of that is needed for an external application to run the
reduced solver.

## A and B construction comparison

| Construction concern | Problem A | Problem B |
| --- | --- | --- |
| Layouts | Two variable blocks (`state`, `control`) plus one test block; both problem vectors have 289 entries | The same block structure; state/adjoint have 225 free entries and control has 289 full entries |
| Five callbacks | Block vectors are passed to `ProblemA`'s full-coordinate residual, derivative, objective, and VJP operations | The same block conversions feed `ProblemB`; its formulation applies free/full coordinate and coupling operations |
| Solve wrappers | State and adjoint callbacks call `ProblemA` and wrap native vectors and actual CG reports | State and adjoint callbacks call `ProblemB` and perform the corresponding free-coordinate wrapping and report translation |
| Metric adaptation | Local identity `MetricT` copies the control block in both directions | Local `ProblemBMassMetric` forwards mass apply/inverse, checks native convergence, and returns the primal block; it exposes no metric report |
| Partition and DTO | One state/control partition; solve services are passed directly while constructing `ReducedDTO` | Same partition/DTO shape, with B's metric and coordinate-aware solve services |
| Lifetime | The binding borrows `ProblemA`, which owns its prepared application | The binding borrows `ProblemB`, which borrows the prepared application; application → problem → binding → metric/DTO/solver must remain alive |
| Output | Retained state is written directly through `ProblemA` | Free state is reconstructed to a full state before the existing Step-4 writer is called |

This comparison identifies repeated public-contract construction and the
real B-specific numerical/ownership additions without assigning subjective
percentages or treating the local class layout as required API surface.

## Source-size observation

These counts are nonblank, non-comment physical lines, omitting lines whose
first non-whitespace characters are `//` and retaining code-bearing lines.
They are descriptive observations, not acceptance thresholds. The counts
reflect MC4's consumer cleanup and MC5's test-only executable validation;
existing application mathematics and adapted source are not counted as new
consumer wiring.

| Scope | Problem A | Problem B |
| --- | ---: | ---: |
| Minimal binding | 206 | 222 |
| Minimal entry point | 107 | 112 |
| Minimal functional source total | 313 | 334 |
| Evaluated nmopt binding | 235 | 336 |
| Focused minimal contract test | 451 | 705 |
| Native reduced evaluation and optimizer | 327 | 459 |

B's functional consumer is 21 code-bearing lines larger in this snapshot,
reflecting its coordinate and mass-metric requirements. The larger B contract
test contains dense/reference and executable checks; those lines do not
belong to the functional consumer.

The [implementation report](../../../../apps/external-dealii/step-4/integration-report.md#4-source-accounting)
extends this inventory to native OCP code, shared Step-4 adaptation, and
comparison support. Its [binding-region breakdown](../../../../apps/external-dealii/step-4/integration-report.md#5-inside-the-minimal-bindings)
reconciles the 206/222 totals to disjoint construction responsibilities.

## Validation and usability limits

MC5 extended the existing A/B contract-test registries so their executable
scenarios launch the actual minimal programs. Each scenario checks completion,
stopping reason, relevant dimensions, iteration counts, objective, gradient
norm, and generated VTK data against a freshly audited native reference. A's
final audit independently recomputes state before adjoint and covector
checks. The focused scenarios and the required Debug pipelines passed at
`f5f2fe5`:

- `debug-dealii`: 195/195;
- `debug-neutral`: 67/67.

The executable and contract outputs remain ignored run artifacts below
`runs/external-dealii/step-4/minimal/`; no tracked run evidence is part of
this follow-up. Successful compilation and execution demonstrate
functionality, not ease of authoring for an unfamiliar user.

## Decision

Keep the two explicit consumers and their ownership separation. The tested
public boundary is functionally sufficient, and no shared nmopt change,
generic helper, control-only pullback API, diagnostic cleanup, or further PDE
experiment is selected by this evidence.

The repeated construction can motivate a separately scoped ergonomic
experiment if a future decision asks for it. That experiment need not be
conditioned on first adding another application or collecting performance
measurements; broader generality and performance are separate questions. Any
such work would require its own scope and acceptance. This follow-up ends at
the current decision gate.
