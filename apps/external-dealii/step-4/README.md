# External deal.II Step-4 integration

This case study connects the authentic deal.II `v9.5.1` Step-4 tutorial to
nmopt's existing reduced optimizer. Problem A adds algebraic RHS control;
Problem B adds FE distributed control with fixed state boundary data and a
mass metric. Both retain the application's native numerical and output
policies. The boundary evaluation and minimal-consumer follow-up are closed.

## Reading paths

| Purpose | Start here |
| --- | --- |
| Understand how the application, formulation, backend, and optimizer work together | [Explanatory overview](external-integration-overview.md) |
| Inspect the Step-4 adaptations, mathematics/API mapping, and source counts | [Implementation report](integration-report.md) |
| Build and run the complete minimal consumers | [Minimal consumer commands](minimal/README.md) |
| Look up exact public types, callbacks, solve reports, and lifetimes | [External API reference](../../../docs/reference/external-dealii-solver-integration.md) |
| Review objectives, evidence, conclusions, and limits | [Closure audit](../../../docs/planning/review/external-dealii-boundary-evaluation/closure-report.md) |
| Inspect the frozen experimental design | [Evaluation protocol](../../../docs/planning/external-dealii-boundary-evaluation.md) and [Problem B protocol](../../../docs/planning/review/external-dealii-boundary-evaluation/problem-b-protocol.md) |

## Source and reproduction

The exact [upstream source](source/upstream/step-4.cc), numerically unchanged
[stripped baseline](source/baseline/step-4-stripped.cc), and reusable
[adapted source](source/adapted/step-4.cc) are retained separately. The adapted
standalone program still runs Step-4's original 2D/3D forward sequence.

Use the [implementation report's reproduction commands](integration-report.md#7-evidence-and-reproduction)
for token/forward fidelity and evaluation checks. The minimal consumers have
their own [build and run commands](minimal/README.md#build-and-run). Generated
VTK output and evaluation evidence belong under ignored `runs/` directories.
The comparison machinery is separate from the minimal consumers' runtime;
the report identifies the shared native code and its optional diagnostics.
