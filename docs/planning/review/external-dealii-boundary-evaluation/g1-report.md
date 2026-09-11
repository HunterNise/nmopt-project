# External deal.II boundary evaluation: G1 report

Status: successful Problem A comparison and completed EC5 failure-evidence
corrections retained on 2026-09-11. Generated run artifacts
remain ignored and are recreated by the
commands recorded below; no run output is copied into tracked documentation.
The original E0–E5 evidence remains identified by its evaluated revision, and
the corrected closure evidence is identified separately.

Date: 2026-09-11
Original evaluated revision: `277fbf4` (`test(dealii): compare native and nmopt optimization paths`)
Corrected numerical path: `d44dded1400365979d1633c9eda8e561a776c32c`
Pre-EC5 closure verification: `6a9d1a565bdf0b21ad72889e6f3e06bcdf5a3282`
EC5 implementation: `0357d4f`, `e8050ff`
Final native trace follow-up: `0bb1307` plus the tested change to
`tests/dealii/external_step4_native_contract.cc`, pending commit
Roadmap: [external deal.II boundary evaluation](../../external-dealii-boundary-evaluation.md)
Review context: [design investigation](design-investigation.md)

## Post-G1 reconciliation

The original report was evaluated at `277fbf4`, before the evidence and
ownership corrections. The following bounded changes were subsequently
reviewed and committed:

- C1 (`2ff6eec`) uses deal.II's solver-control monitored residual as the final
  solve residual, removes the extra post-solve matrix-vector product, and
  keeps failed-solve records separate from successful records. The
  `NoConvergence` catch preserves the solver's last step and residual; the
  finite non-convergence branch was not runtime-triggered because the frozen
  system converges and non-finite inputs are rejected earlier by deal.II.
- C2 (`303e00e`) instruments the local identity metric callbacks directly. The
  latest matched artifact recorded 1,657 metric `apply` calls, 829
  `inverse_apply` calls, and 829 solver-reported metric solves. The latest
  reduced-evaluation artifact recorded zero metric callbacks. The direct
  counts are distinct from the solver-reported summary and are not timing or
  allocation measurements.
- O1 (`8ad5495`) separates source, integration, evaluation, verification, and
  diagnostics ownership. O2 (`b1e83a6`) makes diagnostics optional for the
  functional binding path and prevents the binding from being copied or
  moved because its callbacks capture owned application state.

- EC1 (`a7f6cf7`) retains explicitly recorded comparison failures and solve
  records, and prevents runtime tests from rewriting the shared attribution
  ledger.
- EC2 (`b235a85`, `21d775e`, `6193519`, `d3bdaa0`, and `d44dded`) removes the
  extra solve action, independently audits residuals and final gradients,
  enforces the absolute oracle-distance gate, rejects nonfinite acceptance
  data, and closes the runtime-count checks.
- EC3 (`6a9d1a5`) makes the forward comparator reject nonfinite geometry and
  field values while retaining its failure report.
- EC5 (`0357d4f`, `e8050ff`) records original exception diagnostics at explicit
  catch boundaries and serializes each completed native/current-nmopt trace
  in the paired driver before later checks can fail. Direct-throw and
  post-result-failure retries retain failed status, diagnostics, available
  traces, counters, and solve records.

The corrected reduced comparison is recreated at
`runs/external-dealii/step-4/reduced-evaluation/reduced-1789148523091774/comparison.csv`;
the corrected matched-optimization summary and counters are recreated at
`runs/external-dealii/step-4/optimization/paired-1789148523884470/comparison/summary.txt`
and its `counters.csv`. The corrected working ledger is at
`runs/external-dealii/step-4/working/attribution.csv`. These paths are ignored
run artifacts, not tracked evidence files. The corrections preserve the
original G1 decision: the current public boundary is adequate for this tested
external case, and no shared helper or API change is implied.

### Failure-evidence qualification after EC5

Review at `ed450bd` found that
[`EvidenceGuard`](../../../../tests/dealii/external_step4_evidence.hpp)
could replace an uncaught exception's original message with the generic
`scenario terminated before completion` during destruction. A temporary C++17
probe reproduced this without changing repository sources. The
[paired optimization driver](../../../../tests/application/external_step4_optimization_contract.cc)
also delayed trace serialization until both solvers, convergence checks, and
the oracle returned; a later failure could discard a completed earlier result.

EC5 closes both delayed-evidence cases within the existing Step-4 test
boundary. Each guarded driver catches while its `EvidenceGuard` and
instrumentation remain alive, calls `fail_current_exception()`, and rethrows;
the guard no longer has to infer an exception message during destruction. The
paired optimization driver writes the native trace immediately after the
native solver returns and the current-nmopt trace immediately after the nmopt
solver returns, before later finite, convergence, oracle, output, or
comparison checks.

Review of `0bb1307` found that the standalone native driver still delayed its
trace until after all acceptance audits. The final follow-up separates that
trace from the oracle-dependent summary and writes it immediately after the
solver returns. Its new `native_trace_failure` scenario invokes the actual
driver and throws immediately after serialization, before acceptance audits.
Both retries preserve completed trial/accepted records, the original exception,
counters, and solve records, while later summary/audit files remain absent.

The final focused native optimization/trace-failure selection passed `2/2`.
The full required pipelines passed `177/177` for `debug-dealii` and `67/67`
for `debug-neutral`, including the earlier exception and paired trace tests.
These checks used `0bb1307` plus the native-driver correction named above;
all other numerical sources were unchanged. Generated evidence stays in
unique ignored run directories. No tracked run evidence is required or added.

The correction remains bounded: trials held only inside a solver that throws
before returning a result are not exposed by the current solver interface.
Such a failed run is incomplete evidence; EC5 does not claim full partial
trace retention and does not introduce a progress callback or shared API.
Existing successful results, numerical audits, and the G1 boundary conclusion
are unaffected.

The EC5 checks are reproducible with the existing profiles:

```bash
./build.sh build debug-dealii --target nmopt_external_step4_native_contract_test
./build.sh test debug-dealii \
  --regex '^nmopt\.external_tutorial_step_4\.(native_optimization|native_trace_failure)$'
ctest --test-dir build/debug-dealii --output-on-failure \
  -R '^nmopt\.external\.tutorial_step_4\.failure_evidence$'
ctest --test-dir build/debug-dealii --output-on-failure \
  -R '^(nmopt\.external\.tutorial_step_4\.(matched_optimization|completed_trace_failure))$'
./build.sh pipeline debug-dealii debug-neutral
```

## Decision

The frozen Problem A comparison is successful. The authentic Step-4 source can
be opened through small application-owned seams, its control mathematics can be
verified independently, and the same native and current-nmopt reduced
optimization reaches the same result under the frozen policy. The existing
public boundary is therefore adequate for this tested external case.

No shared nmopt change, new public helper, or boundary redesign is recommended
from this evidence alone. Keep the Step-4 binding and Problem A mathematics
local to the application evaluation. A future capability-selection or
control-only-pullback proposal would need its own reproducible question and
scope; this report does not turn the observed work difference into a timing or
allocation claim.

The conclusion is deliberately bounded. The experiment covers one linear,
symmetric, square 2D Problem A with full-vector load control, identity metric,
and no constraints. It does not establish ergonomic sufficiency for nonlinear,
nonsymmetric, distributed-control, boundary-control, constrained, alternate-
metric, or compiler-produced applications.

## Scope and provenance

The pinned upstream source is deal.II `v9.5.1`, with SHA-256
`be9e694f5f3c9177b7cd18200ff8173337c2b16e1ee72d45ab1ba7c6e105be5f`. The
comment-stripped baseline has SHA-256
`b21212764c50401089612c6ac2bb196395e3ac9c261120e0513be341825399b2`, and the
current adapted reusable source has SHA-256
`c85e04027681d007498ee9b5b885e4b30f9b020c659abee26e703c71024be71b`.
The upstream and stripped token stream was checked by the bounded stripper;
the adapted source was compared separately and is not used to regenerate the
baseline. The obsolete wrapper and smoke implementation from the previous
attempt were removed in `9a339e7` and remain recoverable in Git history.

The closure verification recreated the forward comparisons for upstream versus
stripped at
`runs/external-dealii/step-4/forward-comparison/20260911T174353Z-201717a3/comparison.txt`
and upstream versus adapted at
`runs/external-dealii/step-4/forward-comparison/adapted/20260911T174403Z-82ee1d64/comparison.txt`.
Both reports show identical stdout, 2D output with 1,024 points and 256 cells,
and 3D output with 32,768 points and 4,096 cells, including zero numeric array
differences. The paths are ignored run artifacts and are recreated by the
documented commands.

The pre-EC5 closure verification used the existing Debug profiles and the
machine's configured build limits:

~~~bash
python3 tools/external_dealii/strip_comments.py \
  --input apps/external-dealii/step-4/source/upstream/step-4.cc \
  --check apps/external-dealii/step-4/source/baseline/step-4-stripped.cc
ctest --test-dir build/debug-dealii --output-on-failure \
  -R '^nmopt\.external_tutorial_step_4\.(forward_comparator_contract|forward_comparison|adapted_forward_comparison)$'
./build.sh pipeline debug-dealii
./build.sh pipeline debug-neutral
~~~

The focused forward selection passed `3/3`; the initial complete Debug
deal.II pipeline passed `175/175`, and the backend-neutral pipeline passed
`67/67`. These results predate EC5; the final pipeline counts after EC5 are
recorded in the failure-evidence section above.
The machine was Linux x86_64 under WSL2, using GCC 13.3.0, CMake 3.28.3,
Ninja 1.11.1, and deal.II 9.5.1 from `/usr/share/cmake/deal.II`. No release
timing, allocation count, or internal uninstrumented operator count is
claimed.

## Mathematical and integration evidence

The native Problem A owns the frozen equations

~~~math
E(y,u)=Ay-b-u,
\qquad
J(y,u)=\frac{1}{2}y^{\mathsf T}y+\frac{1}{2}u^{\mathsf T}u.
~~~

The native staged evaluator uses Step-4's assembled matrix, right-hand side,
CG policy, and output machinery. Its independent dense oracle solved
$(I+A^{\mathsf T}A)y=A^{\mathsf T}b$ and derived $u=Ay-b$. The oracle reported
system residual `3.2811052213435387e-16`, stationarity residual
`1.4056590879085018e-15`, and zero matrix-symmetry error. Off-solution
residual/JVP/full-VJP pairings, objective directional derivatives, centered
finite differences, and Taylor remainders passed; the observed Taylor ratios
were within `3.9999999994`–`4.0000000004`.

The E4 paired reduced-evaluation run covered all seven prescribed controls,
including the repeated-control sequence. Every paired state, objective,
adjoint, and reduced-gradient error was zero. Both paths recorded seven state
solves and seven adjoint solves. Native used seven direct control pullbacks;
current nmopt used seven full residual VJPs and seven explicit transpose
actions. The result is preserved in
[the recreated reduced-evaluation artifact](../../../../runs/external-dealii/step-4/reduced-evaluation/reduced-1789148523091774/comparison.csv).

The E5 matched optimization run used the same zero initial control and frozen
steepest-descent/Armijo policy. Both paths stopped by gradient tolerance after
828 accepted iterations and 6,025 total line-search trials. Both recorded:

- one assembly;
- 6,026 state solves and 829 adjoint solves;
- 6,026 objective values and 829 objective-derivative augmentations;
- the same scalar trial/accepted trace and stopping reason;
- zero final control and state difference;
- identical retained-state VTK output.

The final gradient norm was `9.5036543162094535e-7`. Native used 829 direct
control VJPs. Current nmopt used 829 residual VJPs and 829 explicit matrix
transpose actions, plus 829 reported identity-metric solves and no Hessian
actions. The corrected paired summary is
[here](../../../../runs/external-dealii/step-4/optimization/paired-1789148523884470/comparison/summary.txt);
its counters, residual audits, gradient audits, and complete traces are
retained beside it in the ignored run directory.

## Causal attribution

The working ledger is [here](../../../../runs/external-dealii/step-4/working/attribution.csv);
it is intentionally ignored and is recreated or reviewed alongside the run
artifacts rather than copied into tracked documentation.
The categories below are causal classifications, not mutually exclusive
claims about source files.

| Consumer and operation | Observed evidence | G1 classification |
| --- | --- | --- |
| Step-4 preparation and supplied-RHS CG solve (`S01`, `S02`) | The adapted source exposes preparation, matrix/RHS views, and an in/out solve while retaining original assembly and CG. | Application reuse; reused numerical machinery, not nmopt overhead. |
| Objective and control pullback (`M01`, `M02`) | Problem A defines and verifies the identity quadratic objective and $u+p$ reduced derivative. | Capability/formulation obligation: added control mathematics owned by the application. |
| Residual and residual JVP (`V01`, `V02`) | Required and checked for the complete nmopt executable model; neither callback was used by the matched reduced runtime. | Verification capability plus frozen construction obligation; not repeated E5 runtime work. |
| State component of full residual VJP (`N01`) | Current nmopt performed 829 full VJPs and explicit matrix transposes; native used only 829 direct control pullbacks. | Repeated work required by the frozen API. The count is measured; its time and memory impact and performance materiality are unknown. |
| Full-variable block composition (`N02`) | One state block, one control block, and one test block were constructed with checked compatible layouts and borrowed lifetimes. | Mechanical nmopt adaptation. |
| Identity metric realization (`N03`) | A local identity `MetricT` was required; nmopt reported 829 metric solves, with no iterative metric or Hessian work. | Mechanical adaptation and solver service, not a new PDE or metric formulation. |
| Trial and accepted-step orchestration (`O01`) | Native loop and current nmopt matched 6,025 trials and 828 accepted steps, including value-only rejected trials and accepted-state reuse. | Optimizer orchestration service. The native implementation is reference instrumentation, not an nmopt boundary defect. |

The ledger intentionally does not estimate copies, allocations, timing, or
unobserved CG matrix actions. The runtime test also does not collect Git
metadata or source hashes.

## Source-size accounting

Convention: count physical lines containing code after excluding blank lines
and lines consisting only of C++ comments, while counting code-bearing lines
that contain inline comments. Python support counts nonblank, non-comment
physical lines. These are descriptive measurements, not acceptance thresholds.
The line counts describe the original E5 snapshot at `277fbf4`; the paths
below use the current ownership names so the historical measurements remain
locatable after the O1 reorganization.

| Scope | Files or source span | Lines |
| --- | --- | ---: |
| Pinned/derived baseline | `source/upstream/step-4.cc` and `source/baseline/step-4-stripped.cc` | 191 each |
| Adapted reusable Step-4 | `source/adapted/step-4.cc`, including baseline-derived code | 240 |
| Actual reuse seam delta | baseline-to-adapted physical diff | +74 / −15 |
| Native mathematical operations | `verification/scenario.hpp`, `integration/problem_a.hpp`, `evaluation/native_reduced.hpp` | 346 |
| Native optimization orchestration | `evaluation/optimization_policy.hpp`, `evaluation/native_optimization.hpp` | 270 |
| Experiment instrumentation | `diagnostics/instrumentation.hpp` | 59 |
| nmopt application binding | `integration/nmopt_binding.hpp` | 221 |
| Independent verification implementation | `verification/verification.hpp` | 208 |
| Contract/evidence drivers | three Step-4 test drivers | 1,861 |
| Reusable support tools | `tools/external_dealii/*.py` | 376 |

The adapted 240 lines include the 191-line baseline-derived source, so those
rows must not be added. Likewise, test-driver assertions are reported
separately rather than being presented as binding or library implementation.
The source-size result shows that the experiment contains substantial
verification and evidence-driver code; it does not show that the public
boundary is architecturally wrong.

## Hypothesis dispositions

The dispositions are limited to this tested external path. “Supported” means
the evidence supports the bounded case, not the corresponding universal claim.

| ID | Disposition | Evidence and limit |
| --- | --- | --- |
| H1 — Complete executable construction causes material extra obligations | Required extra obligations/work established; materiality unresolved | Residual/JVP/full-VJP construction was required, and the full-VJP state block caused 829 repeated explicit transpose actions. The run measures incidence and count, not timing, allocation, or performance materiality. |
| H2 — Most integration friction is mechanical construction | Unresolved | Block/layout/metric/report wiring is mechanical, but the single binding has no helper comparison and the case also requires real mathematics and repeated full-VJP work. |
| H3 — Verification and runtime capabilities merit different construction status | Supported for the tested case | Residual/JVP were verified and required for construction but had zero matched E5 runtime calls; full VJP remained both a construction and runtime consumer. |
| H4 — Current solve contracts preserve native policy and ownership adequately | Supported for the tested case | Supplied-RHS CG, actual solve reports, zero initialization, retained states, and native exception behavior were preserved. The case is linear and symmetric. |
| H5 — Existing producer/formulation convergence point is sufficient | Supported for the tested external path | The binding used documented public contracts only and required no compiler, private-header, or shared implementation change. Compiler-versus-external equivalence remains untested. |
| H6 — Existing composition seams preserve independent mathematical choices | Supported for the tested case | Problem A objective, control pullback, identity metric, native solves, and output remained application-owned without changing the PDE or shared solver. Alternate formulations remain untested. |
| H7 — Friction reflects deeper structural coupling | Weakened for the tested case | All required changes stayed in the Step-4 evaluation, test, and direct CMake registration layers; no cross-layer correctness obstruction appeared. This does not refute the hypothesis for other applications. |

## Recommendation and limits

Retain the successful E0–E5 Problem A evaluation with the completed EC5
evidence correction. The corrected run outputs remain ignored and reproducible
from the documented commands. Do not add a generic helper or alter the shared
nmopt API as an automatic follow-up.

If future authentic applications reproduce the same block/layout/metric
construction pattern, a separately scoped mechanical helper may be evaluated.
If repeated full-VJP work becomes a demonstrated performance concern, formulate
a separate capability or formulation experiment with direct measurements. Do
not infer that need from the present Debug run.

PB0 subsequently froze the
[Problem B execution protocol](problem-b-protocol.md).
Its [review reasoning](design-investigation.md#8-post-g1-review-and-problem-b-candidate)
remains non-authoritative. B compares another OCP on Step-4; it does not provide
an independent sample of adapting another external application. No B numerical
result is claimed. Nonlinear or nonsymmetric systems, control constraints,
compiler integration, and package/install behavior remain outside that scope.
