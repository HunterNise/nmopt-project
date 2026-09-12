# Minimal Step-4 consumer assessment

Status: complete MC3 assessment for the separately accepted minimal-consumer
follow-up. This report compares the two runnable consumers at commit
`99bef20` and does not reopen the completed G1 or G2 evaluations.

## Result

The current public boundary is sufficient to wire both existing Step-4
problems into the reduced optimizer without changing shared nmopt, the
compiler, or the deal.II backend. The functional code is explicit and
application-owned. It is not a trivial one-line connection: the application
author must construct layouts, all five executable callbacks, the state/control
partition, native solve reports, a metric, a reduced DTO, and the optimizer
policy.

Problem A and Problem B share almost all of that public-contract construction.
The B-specific additions are real application choices, not generic framework
ceremony: B uses a prepared `Step4<2>` borrowed by `ProblemB`, maps between
free state and full state coordinates, and adapts its nonidentity mass metric.
The examples therefore answer the usability question without justifying a
helper or a public API change.

The examples demonstrate functional wiring for these two linear, symmetric,
serial Step-4 cases. They do not establish global minimality, newcomer
usability, performance materiality, or applicability to nonlinear,
nonsymmetric, constrained, distributed, compiler, or packaged applications.

## Ownership map

The following is the boundary a user needs to understand. Existing numerical
code is listed to distinguish reused application capability from newly authored
consumer wiring.

| Area | Problem A | Problem B | Functional or accessory? |
| --- | --- | --- | --- |
| Step-4 discretization and output | `source/adapted/step-4.cc`, reused by `integration/problem_a.hpp` | `source/adapted/step-4.cc`, borrowed by `integration/problem_b.hpp` | Existing application capability |
| PDE/OCP operations | `integration/problem_a.hpp` owns residual, JVP, VJP, objective, derivative, and native solves | `integration/problem_b.hpp` plus `problem_b_coordinates.hpp`, `problem_b_mass.hpp`, and `problem_b_metric.hpp` own the same operations, lifting, coupling, and mass geometry | Existing application formulation |
| Public binding | `minimal/problem_a_binding.hpp` | `minimal/problem_b_binding.hpp` | Required functional wiring |
| Application entry point | `minimal/problem_a.cc` constructs `ProblemA`, zero control, policy, solver, and output | `minimal/problem_b.cc` constructs `Step4<2>`, prepares it, constructs `ProblemB`, zero control, policy, solver, state reconstruction, and output | Required consumer behavior |
| Evaluation baseline | `integration/nmopt_binding.hpp` and native A evaluation | `integration/nmopt_problem_b_binding.hpp` and native B evaluation | Accessory comparison support |
| Independent checks | `tests/application/external_step4_minimal_problem_a_contract.cc` | `tests/application/external_step4_minimal_problem_b_contract.cc` | Accessory verification |
| Evidence and diagnostics | `evaluation/`, `verification/`, `diagnostics/`, and ignored `runs/` artifacts | Same directories and artifact role | Accessory; not needed by the consumer |

The minimal bindings deliberately borrow the existing `ProblemA` or
`ProblemB`; they do not own an evaluation reference, instrumentation, oracle,
comparison ledger, or run manifest. `ProblemA` internally owns its prepared
Step-4 application. `ProblemB` borrows the application supplied by its caller,
so the B entry point makes the application → problem → binding → solver
lifetime order visible. In both cases the application and problem outlive the
callbacks and output operation.

## Minimum functional wiring

The minimum here means the smallest complete path demonstrated by these
examples, not a claim that every external application will need exactly these
lines.

### Shared A/B work

Each consumer author writes:

- a variable layout with `state` and `control` blocks and a residual-test
  layout;
- all five `CallbackExecutableModelT` operations: residual, residual JVP,
  residual VJP, objective, and objective derivative;
- a `StateControlPartitionT` and state/adjoint solve callbacks;
- translation of native CG evidence into `LinearSolveReport` values;
- a `ReducedDTOT`, zero initial control, the frozen reduced-solver policy, and
  the solver invocation; and
- application-owned final-state output.

The reduced first-order path does not call every callback during a successful
optimization, but the current executable-model contract still requires all
five. The state and adjoint callbacks must preserve native convergence and
failure behavior; they cannot report a failed native solve as an exact solve.

### Problem A additions

The A consumer uses the full 289-entry algebraic control and a local identity
`MetricT`. `ProblemA` provides the native operations and retains the adapted
Step-4 application internally. Final output can use the retained state
directly.

### Problem B additions

The B consumer uses 225 free state/adjoint entries and 289 full control
entries. Its binding retains the existing B mass/coupling and coordinate
objects through `ProblemB`, and supplies a local `MetricT` adapter whose
inverse calls the native mass CG solve and returns actual convergence evidence.
The entry point reconstructs the full physical state before calling the
existing Step-4 writer. These are B's mathematical and ownership requirements,
not generic work that should be removed from the example.

## What belongs only to evaluation

The focused contract tests intentionally contain much more code than the
consumers. They create separate native and public instances, run matched
optimization, compare histories and vectors, recompute fresh state/adjoint and
gradient values, inspect dense operators, perform the independent dense mass
stationarity audit and KKT oracle check, compare output payloads, and write
ignored evidence. None of that is needed for an external application to run
the reduced solver.

The older evaluated bindings remain useful as comparison fixtures, but they
are not the minimal path. They add optional instrumentation and, especially
for B, hide application construction inside the evaluated binding. Native
optimization loops, dense verification, operation attribution, exception
evidence, and run artifacts belong to the evaluation layers. Removing them
from `minimal/` is the meaningful cleanup; deleting the historical evaluation
layers would erase a different responsibility.

No tracked run evidence is part of this follow-up. The ignored `runs/` tree
and ignored build-profile artifacts are recreated by the focused tests and
manual commands documented in the [minimal consumer README](../../../../apps/external-dealii/step-4/minimal/README.md).

## Source-size observation

These counts use nonblank, non-comment physical lines, omitting lines whose
first non-whitespace characters are `//` and retaining code-bearing lines.
They are descriptive observations, not acceptance thresholds. Existing
application mathematics and the adapted source are not counted as new
consumer wiring.

| Scope | Problem A | Problem B |
| --- | ---: | ---: |
| Minimal binding | 216 | 232 |
| Minimal entry point | 85 | 90 |
| Minimal functional source total | 301 | 322 |
| Evaluated nmopt binding | 235 | 336 |
| Focused minimal contract test | 294 | 517 |
| Native reduced evaluation and optimizer | 327 | 459 |

The B functional consumer is 21 code-bearing lines larger in this snapshot.
That small difference should not be read as a general scaling law: B's
existing formulation is substantially different, and the line count does not
measure conceptual or runtime cost. The B evaluated binding is 101 lines
larger than A's evaluated binding, reflecting B-specific coordinate and metric
adaptation in addition to its application ownership and instrumentation. The
focused B test is larger because it contains the independent dense and output
checks, not because the B consumer needs those checks.

## Concepts, coupling, and documentation

The public concepts that an external author must learn are:

1. block layouts identify ordered spaces and dimensions; raw vector sizes are
   not enough;
2. callbacks use block primal/covector wrappers and must return compatible
   layouts;
3. the reduced DTO separates the full executable model from state/control
   partition and native state/adjoint solves;
4. solve reports carry native iterative evidence and failure status;
5. a metric maps the reduced covector to the optimizer's search direction; and
6. callback captures impose an application/problem lifetime obligation.

The repeated mechanical conversions are the vector-to-block wrappers,
native solve evidence to public solve reports, native derivative structs to
block covectors, and (for B) free/full coordinate reconstruction. The full
residual VJP is also a repeated public runtime operation in the B path, but
its performance materiality was not measured here. These obligations are
visible in the examples and are not evidence of a correctness defect.

The [public integration reference](../../../reference/external-dealii-solver-integration.md)
documents the common callback, layout, solve, metric, reduced-DTO, and
lifetime contracts. The minimal README adds the concrete A/B ownership map and
commands. A future documentation improvement could add a short example of a
borrowed problem plus an application-specific iterative metric adapter to the
public reference; the current examples are already functional and this is not
an accepted API or implementation task.

No unfamiliar-user exercise was performed. Successful compilation and
execution are evidence of functionality, not evidence that the authoring
surface is easy for a newcomer.

## Decision

Keep the two explicit consumers and their ownership separation. Close the
minimal-consumer follow-up with no generic helper, control-only pullback API,
shared-nmopt change, diagnostic cleanup, or further PDE experiment. Any
mechanical helper or performance study requires a new, separately scoped
question supported by another authentic application or direct measurements.

The report, README, and roadmap now make the minimum functional source and
the accessory evaluation source explicit. The phase has no next implementation
unit under this plan.
