# External deal.II boundary evaluation: G1 report

Status: G1 report reconciled with the bounded post-G1 cleanup on 2026-09-11.
The original E0–E5 evidence remains identified by its evaluated revision;
the reconciliation below records later evidence corrections without changing
the public nmopt boundary.

Date: 2026-09-11  
Evaluated revision: `277fbf4` (`test(dealii): compare native and nmopt optimization paths`)  
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

The corrected paired artifacts are the
[reduced-evaluation comparison](../../../../runs/external-dealii/step-4/reduced-evaluation/1789129214466132/comparison.csv)
and the [matched-optimization summary](../../../../runs/external-dealii/step-4/optimization/1789129275006264/comparison/summary.txt).
These changes preserve the original G1 decision: the current public boundary
is adequate for this tested external case, and no shared helper or API change
is implied.

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
adapted reusable source has SHA-256
`aec5a0e25578ad9b8c762bf74b2f7f333f23dabdd777c77c4e06e3b0be377bc8`.
The upstream and stripped token stream was checked by the bounded stripper;
the adapted source was compared separately and is not used to regenerate the
baseline. The obsolete wrapper and smoke implementation from the previous
attempt were removed in `9a339e7` and remain recoverable in Git history.

The latest recorded forward comparison passed for upstream versus stripped in
[this artifact](../../../../runs/external-dealii/step-4/forward-comparison/20260911T093700Z-ef4e676c/comparison.txt),
and for upstream versus adapted in
[this artifact](../../../../runs/external-dealii/step-4/forward-comparison/adapted/20260911T093708Z-f89b903f/comparison.txt).
Both report identical stdout, 2D output with 1,024 points and 256 cells, and
3D output with 32,768 points and 4,096 cells, including zero numeric array
differences.

The measured commands used the existing `debug-dealii` profile and one build
job where a deal.II build was required:

~~~bash
python3 tools/external_dealii/strip_comments.py \
  --input apps/external-dealii/step-4/source/upstream/step-4.cc \
  --check apps/external-dealii/step-4/source/baseline/step-4-stripped.cc
./build.sh build debug-dealii --target nmopt_external_step4_optimization_contract_test
ctest --test-dir build/debug-dealii --output-on-failure \
  -R '^nmopt\.external\.tutorial_step_4\.matched_optimization$'
./build.sh pipeline debug-dealii
~~~

The focused matched-optimization selection passed `1/1`; the complete Debug
deal.II pipeline passed `172/172`. No release timing, allocation count, or
internal uninstrumented operator count is claimed.

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
[the reduced-evaluation artifact](../../../../runs/external-dealii/step-4/reduced-evaluation/1789119304531046/comparison.csv).

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
actions. The paired summary is
[here](../../../../runs/external-dealii/step-4/optimization/1789119247586875/comparison/summary.txt);
both complete traces are retained beside it.

## Causal attribution

The working ledger is [here](../../../../runs/external-dealii/step-4/working/attribution.csv).
The categories below are causal classifications, not mutually exclusive
claims about source files.

| Consumer and operation | Observed evidence | G1 classification |
| --- | --- | --- |
| Step-4 preparation and supplied-RHS CG solve (`S01`, `S02`) | The adapted source exposes preparation, matrix/RHS views, and an in/out solve while retaining original assembly and CG. | Application reuse; reused numerical machinery, not nmopt overhead. |
| Objective and control pullback (`M01`, `M02`) | Problem A defines and verifies the identity quadratic objective and $u+p$ reduced derivative. | Capability/formulation obligation: added control mathematics owned by the application. |
| Residual and residual JVP (`V01`, `V02`) | Required and checked for the complete nmopt executable model; neither callback was used by the matched reduced runtime. | Verification capability plus frozen construction obligation; not repeated E5 runtime work. |
| State component of full residual VJP (`N01`) | Current nmopt performed 829 full VJPs and explicit matrix transposes; native used only 829 direct control pullbacks. | Frozen API capability with small justified repeated overhead in this case. The count is measured; its time and memory impact are unknown. |
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
| H1 — Complete executable construction causes material extra obligations | Supported for the tested case | Residual/JVP/full-VJP construction was required, and the full-VJP state block caused 829 repeated explicit transpose actions. No timing or allocation materiality is claimed. |
| H2 — Most integration friction is mechanical construction | Unresolved | Block/layout/metric/report wiring is mechanical, but the single binding has no helper comparison and the case also requires real mathematics and repeated full-VJP work. |
| H3 — Verification and runtime capabilities merit different construction status | Supported for the tested case | Residual/JVP were verified and required for construction but had zero matched E5 runtime calls; full VJP remained both a construction and runtime consumer. |
| H4 — Current solve contracts preserve native policy and ownership adequately | Supported for the tested case | Supplied-RHS CG, actual solve reports, zero initialization, retained states, and native exception behavior were preserved. The case is linear and symmetric. |
| H5 — Existing producer/formulation convergence point is sufficient | Supported for the tested external path | The binding used documented public contracts only and required no compiler, private-header, or shared implementation change. Compiler-versus-external equivalence remains untested. |
| H6 — Existing composition seams preserve independent mathematical choices | Supported for the tested case | Problem A objective, control pullback, identity metric, native solves, and output remained application-owned without changing the PDE or shared solver. Alternate formulations remain untested. |
| H7 — Friction reflects deeper structural coupling | Weakened for the tested case | All required changes stayed in the Step-4 evaluation, test, and direct CMake registration layers; no cross-layer correctness obstruction appeared. This does not refute the hypothesis for other applications. |

## Recommendation and limits

Close the E0–E5 evaluation for the tested Problem A case. Do not perform
cleanup, add a generic helper, or alter the shared nmopt API as an automatic
follow-up. Preserve the current attribution ledger and raw artifacts as
working evidence.

If future authentic applications reproduce the same block/layout/metric
construction pattern, a separately scoped mechanical helper may be evaluated.
If repeated full-VJP work becomes a demonstrated performance concern, formulate
a separate capability or formulation experiment with direct measurements. Do
not infer that need from the present Debug run.

Problem B, nonlinear or nonsymmetric systems, constraints, FE-coupled controls,
alternate metrics, compiler integration, and package/install behavior remain
outside this decision. Any of them requires an accepted follow-up scope.
