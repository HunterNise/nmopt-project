# External deal.II Problem B boundary evaluation: G2 report

Status: successful Problem B comparison and completed G2 attribution review,
including independent stationarity and runtime-attribution corrections.
Generated run artifacts remain ignored and are recreated by the commands below;
no run evidence is tracked.

Date: 2026-09-12  
Evaluated source revision: `46aaa7b` (`test(dealii): retain Problem B operation
attribution`)
Protocol: [`problem-b-protocol.md`](problem-b-protocol.md)  
Earlier decision: [`g1-report.md`](g1-report.md)

Final phase disposition and consumer follow-up: [closure report](closure-report.md).
The numerical observations below retain their original evaluated revision.

## Question and bounded result

Problem B tested whether a genuine finite-element distributed control, fixed
essential state data, a rectangular mass coupling, and a nonidentity control
metric can remain application-owned while the same problem is evaluated and
optimized through the current public nmopt contracts.

The answer is yes for this tested Step-4 case. The native reference and the
current public-nmopt path represented the same 2D linear problem, matched on
the reduced evaluations and optimization trace, passed the independent final
acceptance checks, and produced the same retained-state output. No shared
nmopt, compiler, private-header, upstream-source, or stripped-baseline change
was needed.

This is a bounded result. It covers one linear, symmetric, serial, fixed-mesh,
unconstrained distributed-control problem with the Step-4 discretization. It
does not establish the same conclusion for nonlinear or nonsymmetric systems,
MPI, adaptive meshes, constrained controls, boundary controls, alternate
metrics, or compiler-produced external applications.

## Minimum functional wiring

The B path that makes the public connection functional is:

```text
source/adapted/step-4.cc
        -> integration/problem_b_coordinates.hpp
        -> integration/problem_b_mass.hpp
        -> integration/problem_b.hpp
        -> integration/problem_b_metric.hpp
        -> integration/nmopt_problem_b_binding.hpp
        -> existing nmopt public contracts
```

The ownership of that path is deliberately narrow:

| Location | Functional responsibility | Needed by B wiring? |
| --- | --- | --- |
| `source/adapted/step-4.cc` | Existing mesh, FE, assembly, native solve, and VTK output seams. | Yes |
| `integration/problem_b_coordinates.hpp` | Free-state indices, boundary lifting, restriction, reconstruction, and homogeneous adjoint embedding. | Yes |
| `integration/problem_b_mass.hpp` | Full consistent FE mass and free-row rectangular coupling, including applications and transpose applications. | Yes |
| `integration/problem_b.hpp` | B residual, objective, derivatives, and native state/adjoint solves. | Yes |
| `integration/problem_b_metric.hpp` | Application-owned consistent-mass metric and fresh serial CG inverse. | Yes |
| `integration/nmopt_problem_b_binding.hpp` | Public layouts, callbacks, solve services, metric adapter, and reduced DTO construction. | Yes |
| `evaluation/native_problem_b_reduced.hpp` | Native reduced reference used to compare the public path. | No |
| `evaluation/native_problem_b_optimization.hpp` | Native Armijo reference and frozen optimization arithmetic. | No |
| `verification/problem_b_verification.hpp` | Dense KKT oracle and independent equation, derivative, and metric audits. | No |
| `diagnostics/instrumentation.hpp` | Optional counters and solve evidence for the evaluation. | No |
| `tests/dealii/` and `tests/application/` | Contract, native, comparison, and acceptance drivers. | No |
| `source/upstream/` and `source/baseline/` | Provenance and fidelity fixtures. | No |

The integration headers accept nullable instrumentation pointers, so the
functional binding can be constructed without diagnostics. The headers still
have a source dependency on the diagnostics type; removing that dependency is
not necessary to establish the B boundary and is not recommended as an
automatic follow-up.

The B-specific mathematics is not generic wiring. In particular, the
application owns the full control basis, the free-state coordinate map, the
lifting contribution to the objective, the consistent mass matrix and
rectangular coupling, the residual/objective derivatives, and the native
linear solves. nmopt receives those operations through the existing model,
solve-service, metric, layout, and reduced-DTO contracts.

## Execution ownership

The native and public paths have different owners but the same numerical
meaning:

| Concern | Native reference | Public-nmopt path |
| --- | --- | --- |
| FE application and B mathematics | `ProblemB` and the adapted Step-4 application | The same `ProblemB` operations owned by the binding's application instance |
| Reduced evaluation | `NativeProblemBReduced` | Existing `ReducedDTOT` |
| Control gradient | Direct application control pullback $`-B^{\mathsf T}p`$ | Full residual VJP followed by extraction of its control block |
| Metric | `ProblemBMetric` | `ProblemBNmoptMetric`, forwarding to `ProblemBMetric` |
| Optimization orchestration | `NativeProblemBArmijoSolver` | Existing `ReducedSearchSolverT` |
| Output | Adapted Step-4 VTK writer | Same writer policy applied to the retained reconstructed state |

The callback model still requires residual, residual JVP, residual VJP,
objective, and objective-derivative operations at construction. The matched
optimization used no residual value or JVP callbacks; the public path used the
full residual VJP for the initial derivative and each accepted derivative
update. This is a distinction between construction requirements and runtime
consumers, not a missing native operation.

## Evidence and observed results

PB1 through PB5 were completed from the frozen B protocol. The implementation
revision sequence is recorded in the ignored unit plans. The relevant
regenerated artifacts are:

- Native derivative and metric checks:
  `runs/external-dealii/step-4/problem-b/native-verification/derivative-metric-1789201388478789/`.
- Native/public reduced-evaluation samples:
  `runs/external-dealii/step-4/problem-b/reduced-evaluation/comparison-1789201522707939/`.
- Corrected native/public matched optimization and final acceptance:
  `runs/external-dealii/step-4/problem-b/optimization/paired-1789207421441367/`.
- Corrected standalone native optimization attribution:
  `runs/external-dealii/step-4/problem-b/optimization/native-1789207286761907/`.

The five reduced-evaluation controls—zero, constant, ramp, alternating, and
repeated ramp—matched in objective, state, adjoint, reduced covector, and
mass-metric gradient. The largest recorded full-state difference was
`3.3e-15`; the other recorded comparison errors were zero for that run.

The matched optimization started both paths from zero control. Both stopped
by `gradient_tolerance` after five accepted iterations and five line-search
trials, with `first_divergence none`. Each path recorded six state solves and
six adjoint solves. Native recorded six direct control pullbacks; the public
path recorded six full residual VJPs. Each path recorded six metric inverses
and eleven metric applications. The corrected evidence retains six successful
`gradient_norm` metric-solve records per path; both paths have iterations
`29, 27, 25, 21, 16, 12`, with the corresponding initial and final residuals
in `metric-solve-records.csv`.

The corrected runtime matrix-action ledger separates the explicit actions by
their purpose. The counts below exclude separately labeled verification rows:

| Matrix action and purpose | Native | Public nmopt |
| --- | ---: | ---: |
| Coupling apply — state solve | 6 | 6 |
| Mass apply — objective | 12 | 12 |
| Mass apply — objective derivative | 12 | 12 |
| Coupling transpose apply — control VJP | 6 | 0 |
| Stiffness transpose apply — residual VJP | 0 | 6 |
| Coupling transpose apply — residual VJP | 0 | 6 |
| Mass apply — metric apply | 11 | 11 |

The native state and adjoint stiffness work is internal to its existing linear
solve and is not presented as an explicit action count. Likewise, the metric
ledger retains CG iterations and residuals but does not infer unobservable
internal sparse-matrix actions. These are work-attribution limits, not zero
work claims.

For each optimizer output, final state, adjoint, and covector were recomputed
with a fresh verification instance. The final stationarity norm below was then
computed by an independent dense mass factorization, rather than from the
optimizer's stored CG metric result. The independent final acceptance
reported the following common values:

| Quantity | Native | Public nmopt |
| --- | ---: | ---: |
| Final mass-gradient norm from dense audit | `4.8845615376171785e-08` | `4.8845615376171785e-08` |
| Mass-norm distance to dense oracle control | `4.6930791164347435e-08` | `4.6930791164347435e-08` |
| Signed objective gap | `3.5527136788005009e-15` | `3.5527136788005009e-15` |
| Recomputed state residual | `1.1100296398130756e-14` | `1.1100296398130756e-14` |
| Recomputed adjoint residual | `1.0059078255252107e-13` | `1.0059078255252107e-13` |

The dense KKT oracle independently reported state stationarity
`5.7964881954658515e-15`, control stationarity
`6.34455338340255e-18`, feasibility
`1.2986067586327786e-15`, and oracle gradient norm
`2.2993555877574972e-15`. The two retained VTK
payloads were each 26450 bytes after removing only the generated timestamp
header, and their numerical payloads were equal. The output and later audits
left the optimization runtime counters unchanged.

## Attribution

The incremental B obligations fall into four bounded categories:

| Category | Observed responsibility | Interpretation |
| --- | --- | --- |
| Application formulation | Coordinates, lifting, consistent mass, rectangular coupling, B objective/residual, and native solves. | Capability/formulation obligation owned by the application. |
| Mechanical public adaptation | Two variable blocks, one test block, callback packaging, solve-report translation, and the mass-metric adapter. | Local construction work required to express B through the existing public contracts. |
| Native/public runtime difference | Six native direct control pullbacks versus six public full residual VJPs; six metric inverses and eleven metric applications per path; explicit matrix actions are retained by purpose in the corrected ledger. | Measured incidence for this run. The state part of the public VJP and the internal work of CG are not performance measurements. |
| Evaluation support | Native reference loop, independent dense oracle, residual/derivative audits, counters, trace comparison, output comparison, and failure evidence. | Verification and diagnostics, not functional binding code. |

The PB1–PB5 implementation diff from the recorded pre-B starting revision
`0bb1307` through `514c2d9` contains 5841 insertions and 10 deletions across
13 files. Most of that expansion is native contract testing, independent
verification, and comparison/evidence drivers. The G2 correction commits
through `46aaa7b` add only verification and diagnostic plumbing in the
experiment-local layers; they do not expand the minimum functional wiring or
change the public boundary. These source-size observations do not define a
minimum adapter size or establish an external application cost; the functional
wiring is the smaller ownership slice listed above.

No timing, allocation, cache, or internal sparse-matrix work claim follows
from these counts. The full residual VJP is a confirmed repeated operation in
the public path, but its performance materiality remains unresolved.

## Limits and decision

Problem B strengthens the tested-case conclusion from G1: fixed essential
state data, a rectangular FE coupling, and a nonidentity mass metric did not
require a shared API change or reveal a correctness obstruction in the current
public boundary. The required code stayed in the external application's
integration/evaluation/verification and test layers.

No generic helper, control-only-pullback API, or boundary redesign was selected
by G2. Closure reconciliation on 2026-09-12 removes the earlier requirement
to add another authentic application before considering a construction helper:
the subsequent [minimal-consumer comparison](minimal-consumers-report.md)
already identifies repeated construction. That can motivate a separately
scoped ergonomic experiment; broader generality and performance are different
questions. Timing the repeated full-VJP work would likewise require its own
scope. Neither implementation is selected by this closure.

The Problem B evaluation is therefore closed. The shared-nmopt freeze remains
active, and no follow-on implementation unit is accepted without a new,
explicitly scoped question.

## Reproducible verification

The final checks were run from `/root/nmopt-project` on the source tree at
`46aaa7b`:

```bash
./build.sh build debug-dealii \
  --target nmopt_external_step4_native_contract_test \
  --target nmopt_external_step4_problem_b_optimization_contract_test
build/debug-dealii/bin/nmopt_external_step4_native_contract_test \
  native_problem_b_optimization_evidence
build/debug-dealii/bin/nmopt_external_step4_problem_b_optimization_contract_test \
  problem_b_nmopt_matched_optimization
./build.sh pipeline debug-dealii
./build.sh pipeline debug-neutral
```

The focused scenario passed. The complete deal.II Debug pipeline passed
`191/191` tests, including the Problem B matched-optimization scenario. The
backend-neutral Debug pipeline passed `67/67` tests. The run paths above are
ignored artifacts recreated by these tests; no tracked evidence file is
required.
