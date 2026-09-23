# External deal.II Step-4 integration

This self-contained case study asks a narrow question: can an existing deal.II
application remain the owner of its mesh, finite-element discretization, linear
solves, boundary treatment, and output while using `nmopt` as an optimization
library?

For the two tested Step-4 optimal-control problems, the answer is yes. The
adapted tutorial retains its numerical implementation, the OCP mathematics stays
application-owned, and the canonical `nmopt` consumers are confined to explicit
bindings under `minimal/`. The much larger comparison, verification, and evidence
machinery is kept outside that consumer path.

Problem A is a deliberately simple algebraic RHS-control case. Problem B is the
stronger test: it adds free/full state coordinates, fixed boundary lifting, a
rectangular finite-element control coupling, and a nonidentity mass metric without
requiring a shared `nmopt` API, compiler, formulation, or optimizer change.

## Reading paths

| Purpose | Start here |
| --- | --- |
| Review the experiment, its result, evidence, and limits | [External integration overview](external-integration-overview.md) |
| Understand the wider `nmopt` architecture and producer paths | [Project overviews](../../../docs/manual/overview/README.md) |
| Inspect exact source responsibilities, LOC accounting, and reproduction details | [Implementation report](integration-report.md) |
| Build and run the canonical external-consumer examples | [Minimal consumer commands](minimal/README.md) |
| Look up exact public callbacks, solve reports, metrics, and lifetime rules | [External API reference](../../../docs/reference/external-dealii-solver-integration.md) |
| Inspect the historical evaluation decisions and numerical evidence | [Closure audit](../../../docs/history/reviews/external-dealii-boundary-evaluation/closure-report.md) |

## Experiment layout

The directory roles are intentionally separated:

```text
source/        preserved upstream, stripped baseline, and adapted Step-4
integration/   application-owned OCP mathematics and native numerical services
minimal/       canonical external-consumer nmopt bindings and executables
evaluation/    instrumented nmopt bindings and independent native references
verification/  derivative/equation/oracle checks
diagnostics/   optional counters and evidence records
```

`integration/adapted_step4.hpp` is the single private reuse seam that consumes the
preserved adapted tutorial without its standalone `main()`. No `nmopt` include or
type is introduced into `source/adapted/step-4.cc`.

Generated VTK output and evaluation evidence remain ignored below `runs/`. The
comparison and verification layers exist to establish confidence in the experiment;
they are not dependencies of the minimal consumers.

A fresh reproduction at `170c9f1` re-established the historical A/B numerical
results and forward VTK behavior without overwriting the earlier evidence. The
[closure audit](../../../docs/history/reviews/external-dealii-boundary-evaluation/closure-report.md)
records the comparison.
