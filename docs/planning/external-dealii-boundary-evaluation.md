# External deal.II boundary evaluation roadmap

This roadmap is the execution contract for evaluating the current external Step-4 boundary. E0 is adopted; its protocol was committed as `e9ae1ac`, promoting the reviewed local handoff at `.codex/plans/external/nmopt-external-dealii-boundary-evaluation-handoff.md`. This tracked roadmap is the current authority; the local handoff is supporting source material. Subsequent implementation does not require that handoff or the review conversation.

## Phase status

| Field | Status |
| --- | --- |
| Phase | External deal.II boundary evaluation on `codex/evaluate/external-dealii-boundary` |
| Current unit | E2.b — Obsolete Step-4 wrapper retirement and final E2 evidence; implementation complete, review pending |
| Last completed/adopted gate | E2.a — Minimal reusable Step-4 seams and native reuse contract, verified and committed as `f7e44e7` |
| Next unit | E3 — Native Problem A, after the E2.b cleanup commit |
| Shared-nmopt freeze | Active through G1 |

Maintain progress status here. Future units also update their required evidence, attribution records, and runnable documentation. Explicitly accepted protocol amendments must be recorded with their rationale; centralizing status does not prohibit those updates.

### Protocol and evidence labels

- **Frozen protocol** — choices adopted by E0, including numerical constants and acceptance thresholds. They remain fixed for the first comparison unless a documented protocol amendment is accepted. They are not experimental observations.
- **Derived schedule invariant** — an operation/count relationship implied by the prescribed successful execution schedule. It must be tested; deriving it does not make it an observed result.
- **Observed result** — evidence actually produced by a run. Record the revision, environment, command, and evidence location. Prior recorded results remain historical until reproduced for this comparison.
- **Open hypothesis** — an architectural interpretation awaiting evidence. Whether the current boundary is too broad, ergonomic helpers would suffice, or deeper coupling exists remains open through G1.

Sections 3–8 define the frozen protocol, except for explicitly identified derived schedule invariants and evidence-recording guidance. Numerical constants and acceptance thresholds are **Frozen protocol; not yet validated experimentally.** No E1–E5 result is asserted by E0.

All code paths and paths in backticks below are relative to the repository root unless a table gives a narrower base. Markdown links are relative to this document. Repository agent instructions continue to govern execution. This roadmap defines experimental scope and decisions; it does not authorize switching branches, committing, or changing machine configuration independently of the user's implementation instructions.

## 1. Objective, scope, and history

Evaluate the current external deal.II first-order reduced-optimization boundary using upstream Step-4, while distinguishing application reuse, added control mathematics, verification, native reduced and optimizer orchestration, nmopt construction obligations, nmopt runtime work, and any cross-layer coupling.

The first completed comparison is Problem A through native and current-nmopt paths. It must establish equivalent mathematics and matched selected numerical policies before attributing integration burden. A correct, modestly sized binding is a valid outcome. A reproducible inability to integrate under the frozen API is also a valid diagnostic outcome, but not successful completion of the numerical comparison. Name the two comparisons "reduced evaluation" and "matched optimization"; reserve Problem A/Problem B for mathematical problems to avoid confusing the two axes.

Baselines confirmed during review and rechecked for E0 on 2026-09-09:

- Completed refactor: `codex/refactor/pde-solver-boundary`, `67104fd376a9b1d25d7e6cdaca664a7fe1aa98ed`.
- Previous tutorial attempt: `codex/external-dealii-tutorial`, `b4dce25ab03d893df48ab1d9eee2e4dc7cf9183d`.
- Phase branch: `codex/evaluate/external-dealii-boundary`, created from the recorded tutorial tip `b4dce25ab03d893df48ab1d9eee2e4dc7cf9183d`.

Verify these references before implementation. If the current tree or references have advanced, record and reconcile the difference; never reset user work to make the references match. Keep the completed refactor baseline and previous attempt available in Git. Do not merge the experiment backward into the refactor baseline. Commit boundaries below are prospective until commits are explicitly authorized.

### Scope through G1

Permitted changes:

- The Step-4 experiment directory, its new focused tests, and directly necessary target/test registration in `CMakeLists.txt`.
- Evaluation documentation and factual documentation corrections listed in E0.
- Small experiment-local data structures, an identity metric implementing the existing public contract, and explicit instrumentation.
- Replacing the obsolete tutorial wrapper/smoke target as specified in E2, preserving its evidence in Git.

Frozen:

- All shared nmopt library implementation and interfaces, including `include/nmopt/` and `src/` if present.
- Compiler/lowering, generic solver algorithms, shared backend/storage implementations, existing generic helpers, and existing test-support implementations.
- The exact upstream Step-4 file and, after E1, the stripped baseline.

Non-goals:

- New public helpers or formulation ports; a general external-application framework; PDE-family inheritance; generic solver/output hierarchies.
- Recipes, `ProblemSpec`, compiler products, manifests, parameter parsers, a second runner, or a configurable experiment framework.
- MPI, adaptive meshes, nonlinear PDEs, new optimization methods, package/install support, benchmark replication, or numerical timing claims.
- Reducing line counts by changing the mathematical problem or moving obligations into hidden wrappers.

Necessary local fixes to experiment code are allowed. A shared-library bug or missing capability must be reduced to evidence and taken to the decision gate; it is not an exception authorizing a shared-library fix during the measurement.

## 2. Authority and required reading

Before every unit, read [agent routing](../../.agents/README.md), follow its routing, inspect current repository state when required, then read this roadmap's unit and phase status. Reuse already-read, unchanged instructions within a session. Use the [documentation map](../README.md) for routing instead of reading all design documents.

Relevant action instructions are `.agents/workflow.md`, `.agents/git.md`, `.agents/code.md`, `.agents/documentation.md`, `.agents/build.md`, and `.agents/run.md` when their actions apply. This experiment is not a Chapter 6 run and does not adopt Chapter 6 parameter or manifest schemas merely because run instructions mention them.

Authority split:

- This roadmap owns current work sequence, frozen experiment conventions, gates, and status.
- The [PDE–solver boundary](../design/pde-solver-boundary.md) and [v0 executable contract](../implementation/v0/executable-contract.md) remain the architecture and current executable-contract authorities. Evaluate the current contract faithfully; its ergonomic sufficiency is under investigation.
- The [external integration reference](../reference/external-dealii-solver-integration.md) describes the existing public API and tested reference consumer, with the factual qualifications below.
- The [tutorial roadmap](external-dealii-tutorial-roadmap.md) is a superseded historical plan for this work. Its original intended sequence remains historical planning; this roadmap owns current evaluation work.
- Raw traces, working attribution, and speculative explanations remain ignored working evidence. At G1, promote the reviewed factual report and attribution; long-lived design changes require a subsequent accepted decision.

E0 factual corrections:

1. The old roadmap's stale assertion that source selection had not started is replaced with a dated historical qualification: `b4dce25` contains the upstream source, an adapted copy, and a native wrapper/smoke target; it does not contain the planned authentic Step-4 nmopt comparison. Do not retroactively mark all old acceptance criteria complete.
2. The integration reference's opening next-step paragraph points to this evaluation roadmap. The existing fixture exercises API functionality, while authentic adaptation cost and ergonomic sufficiency are under evaluation.
3. Qualify the reference statement that nmopt receives only operations needed by the selected formulation: construction currently requires all five executable operations, although first-order reduced evaluation does not call residual or JVP and consumes only the control component of the full VJP.
4. Correct the layout paragraph: `BlockLayout::compatible_with()` compares ordered space IDs and dimensions. Separate layout objects, and different display labels, can be compatible. Object pointer identity is not required. Equal raw dimensions alone are insufficient.
5. Clarify solve reporting: the one-argument `FormulationSolveResultT` constructor fabricates a converged exact-solve report; it does not verify a native solve. Iterative success should carry actual evidence. A failure must never be labeled converged; returning a failed report or propagating the native exception both prevent reduced evaluation from using a successful result. This experiment uses exception propagation.

No long-lived design rewrite is justified by the reviewed evidence. Do not edit design text merely because a design choice is being tested.

Unit-specific reading uses contract-header names relative to `include/nmopt/contract/`, solver-header names relative to `include/nmopt/solvers/`, and `dealii/serial_backend.hpp` relative to `include/nmopt/`. Step-4 paths are under `apps/external-dealii/step-4/`; `CMakeLists.txt` is at the repository root. These are current-source references, not proposed replacement interfaces.

| Unit | Documents and exact implementation areas |
| --- | --- |
| E0 | `docs/README.md`; current tutorial roadmap; external integration reference; boundary document sections on decision, external application, solves, and output; v0 contract sections on types and reduced formulation. Inspect `layout.hpp`, `callback_executable_model.hpp`, `linear_solve.hpp`, and `reduced_dto.hpp` to support factual corrections. |
| E1 | Step-4 `README.md` and exact upstream source; existing Step-4 targets in `CMakeLists.txt`; build and documentation instructions. |
| E2 | Upstream and stripped sources; current adapted source and obsolete wrapper; boundary sections on solves and output; installed deal.II `matrix_tools.h`, `solver_control.h`, and `solver_cg.h` only as needed to verify initialization and reporting. |
| E3 | Frozen Problem A and policy sections here; E2 seams; boundary mathematical sign/metric conventions; installed deal.II dense factorization API used by the verification oracle. |
| E4 | External integration reference; v0 executable contract; `callback_executable_model.hpp`, `layout.hpp`, `linear_solve.hpp`, `metric_constraint.hpp`, `reduced_dto.hpp`, `dealii/serial_backend.hpp`; existing `tests/application/external_application_dealii_contract.cc` as API evidence. |
| E5 | `reduced_search.hpp`, `reduced_line_search.hpp`, `reduced_gradient.hpp`; staged-evaluation and work-count tests in `tests/contract/reduced_dto_contract.cc`; frozen optimization policy below. |
| G1 | Completed E0–E5 evidence, ledger, final diff, and accepted boundary; [design investigation](review/external-dealii-boundary-evaluation/design-investigation.md) as non-authoritative review context. Read other subsystems only to investigate a named finding. |

## 3. Source organization and ownership

Retain the current directory root `apps/external-dealii/step-4/` for the
tutorial source fixtures. Reusable external-deal.II tooling lives under
`tools/external_dealii/`; do not create an additional application or numerical
library hierarchy to realize the table.

| Path under that directory | Responsibility |
| --- | --- |
| `upstream/step-4.cc` | Exact pinned source; never edited. |
| `baseline/step-4-stripped.cc` | Mechanical comment-stripped baseline, retaining legal notice and all numerical tokens, including both standalone dimensions. Frozen after E1. |
| `step-4.cc` | Adapted tutorial, directly comparable to the stripped baseline. Original numerical machinery plus reuse seams only. |
| `evaluation/scenario.hpp` | Plain native constants and deterministic sample-vector construction. No nmopt includes, parser, recipe, registry, or selection framework. |
| `evaluation/problem_a.hpp` | Native Problem A operations and native vector/solve-result types. Owns control meaning, residual/derivative formulas, objective, adjoint interpretation, and control pullback. |
| `evaluation/native_reduced.hpp` | Native value evaluation and derivative augmentation, with explicit state reuse. No nmopt types. |
| `evaluation/native_armijo.hpp` | E5-only bounded native steepest-descent/Armijo loop and its trace. No policy templates for other algorithms. |
| `evaluation/nmopt_binding.hpp` | Public-API construction, explicit callbacks, report translation, identity `MetricT<SerialBackend>`, and binding lifetime. All nmopt-specific experiment code belongs here or in its test driver. |
| `evaluation/instrumentation.hpp` | Plain counters and evidence records. No framework, allocator replacement, or numerical policy decisions. |
| `evaluation/verification.hpp` | Native oracle, independent residual/derivative checks, result comparison, and evidence serialization. No nmopt numerical dependency. |
| `tools/external_dealii/strip_comments.py` | Bounded reusable source transformation and token-preservation check. Preserves literal contents, token separation, and legal notice. Never rewrites upstream or unrelated files. |
| `tools/external_dealii/check_forward.py` | Focused reusable harness running supplied forward executables in separate directories and comparing stdout/numerical VTK content. No experiment configuration system. |
| `README.md` | Source lineage, exact runnable commands, mathematical name, output locations, and interpretation limits. |

New test drivers:

- `tests/dealii/external_step4_native_contract.cc`: E2 reuse checks, E3 native mathematical checks, and E5 native optimizer scenario.
- `tests/application/external_step4_nmopt_contract.cc`: E4 paired evaluation and E5 paired optimization.

The headers may contain small inline/native definitions. Do not introduce PImpl or a new public façade solely to split files. For this tutorial, `problem_a.hpp` may include the adapted template source behind a tutorial-local `STEP4_NO_MAIN` guard, exactly once per translation unit. This is a recorded source-packaging seam, not an nmopt dependency. Keep the standalone `main()` behind the same guard; do not add nmopt names to the adapted tutorial.

Native targets must not include nmopt library headers or link nmopt numerical/contract targets. The standard-only test-discovery harness may be used by test drivers; it is test infrastructure, not a native numerical dependency. Do not use `nmopt_add_dealii_executable` for native targets because that helper links `nmopt_dealii_contract`. Declare those targets locally, use `deal_ii_setup_target`, existing build flags where needed, and existing scenario discovery. Do not change the shared helper. This local construction follows the phase's explicit independence requirement.

Keep the existing raw-upstream target. Add separate stripped and adapted standalone targets. Keep new contract-test scenarios in the established discovery mechanism. A focused Python forward-comparison CTest may invoke those unmodified standalone programs because they cannot implement scenario discovery without altering the baselines.

In E2, retire `tutorial_application.hpp`, `tutorial_application.cc`, `tutorial_binding_smoke.cc`, and their old binding-smoke target when the replacement native reuse checks are present. Their prior experiment remains in Git. Do not retain FE accessors or a compatibility wrapper merely to keep the obsolete mass/UMFPACK experiment running.

Ownership details:

- Step-4 owns mesh, FE, DoFs, original assembly/boundary elimination, assembled matrix/RHS, CG machinery, and VTK field writing.
- Problem A owns mathematics. Residual/JVP/VJP formulas remain mathematical operations even when verification is their first consumer; the verification harness owns the checks, not a second production definition of those formulas.
- The native reduced evaluator owns solve/objective/adjoint/pullback ordering. The nmopt path delegates that ordering to `ReducedDTOT`.
- Scenario/driver owns initial guesses for new controlled solves, optimization constants, run directories, and output timing.
- The binding owns conversion and explicit borrowed lifetimes. Construct numerical application/problem, model, partition, metric, DTO, and solver in an order that outlives every callback. Keep the binding immovable/noncopyable if callbacks capture its members by address. No lifetime token is required for a deliberately scoped borrowed binding; do not borrow a temporary or return a service whose application has died.

## 4. Frozen Problem A

The control experiments use Step-4 in 2D only: its existing square $[-1,1]^{2}$, four global refinements, `FE_Q(1)`, and original quadrature/forcing/essential-data assembly. The expected full dimension is 289. Both baseline and adapted standalone programs must still run the original 2D and 3D cases.

Let $A,b$ denote the matrix and RHS after Step-4's original boundary elimination. All control and state entries use those full algebraic coordinates.

```math
E(y,u)=Ay-b-u=0,
\qquad y,u\in\mathbb{R}^{n},
\qquad J(y,u)=\frac{1}{2} y^{\mathsf T}y+\frac{1}{2} u^{\mathsf T}u.
```

Zero target, regularization coefficient 1, identity control influence, Euclidean metric, no optimization constraints, initial control zero. This is an algebraic load-vector architecture probe. It is neither a physical distributed-control nor a Dirichlet-control discretization. Perturbing already-eliminated boundary rows is intentional; do not clamp the controlled state afterward or zero boundary components of controls/derivatives.

Mathematical operations:

```math
\begin{aligned}
E'(y,u)(v_{y},v_{u})&=Av_{y}-v_{u},\\
E'(y,u)^{\ast}q&=(A^{\mathsf T}q,-q),\\
J_{y}(y,u)&=y,\qquad J_{u}(y,u)=u,\\
Ay&=b+u,\\
A^{\mathsf T}p&=y,\\
j'(u)&=u+p,\qquad g(u)=j'(u),\qquad d(u)=-g(u).
\end{aligned}
```

The full VJP must be valid for arbitrary supplied $q$, not only for an adjoint returned by a solve. Do not replace its state block with $`J_{y}`$ merely because the exact adjoint equation would identify the two at one particular point.

For this pinned assembled symmetric operator, the native adjoint can reuse the same CG solve. Verify symmetry once before relying on that reuse. This is not a generic adjoint rule.

Independent optimum oracle:

```math
(I+A^{\mathsf T}A)y^{\ast}=A^{\mathsf T}b,
\qquad u^{\ast}=Ay^{\ast}-b.
```

Build and factor the small dense oracle system using deal.II dense/LAPACK facilities. Keep its storage, factorization, and work in verification. Do not call native reduced evaluation, nmopt, or a reduced-Hessian implementation to produce the oracle. Do not form a dense inverse. Independently evaluate its system residual and the stationarity expression $`y_{\mathrm{oracle}} + A^{\mathsf T}u_{\mathrm{oracle}}`$; normalize the latter by $`\max(1, \lVert y_{\mathrm{oracle}}\rVert_{2}, \lVert A^{\mathsf T}u_{\mathrm{oracle}}\rVert_{2})`$. The oracle's matrix extraction/transpose/indexing code should not reuse a production reduced-gradient implementation.

Since the reduced Hessian is $I + A^{-\mathsf T} A^{-1}$, its smallest eigenvalue is at least one. This gives a useful control-error interpretation of the reduced gradient norm, subject to solve and oracle accuracy.

## 5. Frozen solve, initialization, output, and execution policies

### Original standalone initialization

`setup_system()` zero-initializes/reinitializes the solution vector. The original call to `MatrixTools::apply_boundary_values()` then presets essential entries to the prescribed data, leaving unconstrained entries zero. The original `solve()` must pass that same current vector to CG. Do not zero it, substitute an analytic full field, or reconstruct it through a new coordinate map before the standalone solve.

### Reusable solve seam

Prefer the narrow in/out form:

```cpp
SolveInfo solve_system(const Vector &rhs, Vector &solution) const;
```

`solution` is correctly sized on entry and contains the initial guess; on successful return it contains the solved vector. This is an application-local operation, not a new nmopt interface. No algorithm, tolerance, preconditioner, transpose-mode, cache-policy, or variadic options are accepted. Document nonaliasing of RHS and solution and preserve the input RHS.

Inside the seam, retain the original `SolverControl(1000, 1e-12)`, `SolverCG<Vector<double>>`, default CG additional data, and `PreconditionIdentity`. Construct the solver/control per call, as the original code does. Reuse the assembled matrix for the entire scenario. Do not reassemble, factorize, or switch algorithms per control.

`SolveInfo` contains only native convergence evidence needed by callers, such as iterations and the solver-control final residual value. Original `solve()` invokes this seam using its existing `system_rhs` and `solution`, and prints the original iteration message. The seam itself does not print an optimization trace.

For every controlled state solve and every adjoint solve, allocate/reinitialize a correctly sized vector to zero and pass it to the same seam. Never warm-start from a previous control, accepted state, previous adjoint, or the original boundary-initialized standalone vector. This deterministic policy is identical in native and nmopt paths.

Consequently, zero-control solutions must agree with the standalone solution within numerical tolerance, but their CG iteration counts need not equal standalone counts. Native and nmopt zero-control counts should agree with each other. Standalone upstream/stripped/adapted counts should agree with each other on the same build/environment.

On CG failure, retain native exception propagation. The top-level experiment records the last iteration/residual if available, marks failure, and does not continue with the vector. Do not retry using UMFPACK, relax tolerances, or relabel failure as a rejected Armijo trial.

Translate successful solves into `LinearSolveReport` with algorithm CG, identity preconditioner, maximum iterations 1000, actual iterations, relative tolerance 0, absolute/requested tolerance `1e-12`, solver-control final residual, and converged termination. State explicitly that `achieved_residual` is the solver's monitored residual. Independently recomputed residuals belong to verification and must not be confused with that monitor.

### Other frozen policies

| Concern | Policy |
| --- | --- |
| Precision/backend | `double`, native `dealii::Vector<double>`; nmopt uses the existing `SerialBackend` unchanged. |
| Discretization | Exact pinned Step-4 setup; no adaptive refinement, renumbering, extra boundary elimination, or matrix rescaling. |
| Assembly lifetime | Once per application instance. Separate deterministic native/nmopt instances for comparison; verify their $A,b$ agree before evaluating. |
| Caching | No cross-control state/adjoint cache and no factorization cache. Preserve one accepted trial's computed state for derivative augmentation. |
| Execution | Sequential comparisons in one recorded compiler/deal.II/build environment. Record thread settings and numeric compilation flags; avoid fast-math or path-specific numerical flags. |
| Field writer | Original `DataOut` patch construction, VTK format, and field name `solution`. A supplied-state writer may take an output stream; original standalone filename selection remains intact. |
| Field output | Standalone writes its original files. Controlled runs write the final accepted state once, from the retained evaluation, through Step-4's writer. No extra state solve for output. |
| Trace output | Buffer simple scalar records; serialize outside the measured numerical calls. Same numeric precision in evidence, preferably `max_digits10`. |
| Verification work | Separate counter scope; no residual audits or oracle construction counted as optimization runtime. |

## 6. Matched optimization policy

**Frozen protocol; not yet validated experimentally.** Use these first-comparison constants unchanged on both paths. They are not measured convergence claims.

| `ReducedSolverParameters` field | Value |
| --- | --- |
| `maximum_iterations` | 5000 accepted iterations |
| `maximum_line_search_trials` | 30 total trials per iteration, including the initial trial |
| `gradient_tolerance` | `1e-6` |
| `stopping_criterion` | `gradient_norm` explicitly |
| `relative_gradient_tolerance` | `0.0` |
| `objective_change_tolerance` | `0.0` |
| `step_tolerance` | `0.0` |
| `objective_target` | unset |
| `initial_step_length` | `1.0`, restarted on every iteration |
| `minimum_step_length` | `0.0` |
| `armijo_fraction` | `1e-4` |
| `backtracking_factor` | `0.5` |

The `1e-6` stopping threshold is a deliberate first-probe choice to reduce cancellation-driven line-search ambiguity close to a nonzero objective minimum. It still has a direct control-error interpretation for this strongly convex problem. Do not silently substitute the library default `1e-8`. If these constants prevent completion, preserve the failure and review the numerical protocol before changing either path.

Native loop behavior:

1. Evaluate state, objective, adjoint, and reduced derivative at zero control once.
2. At the beginning of each iteration, form $g=j'(u)$, $d=-g$, and $\sqrt{g^{\mathsf T}g}$ using the native vector dot product. If this norm is at most $10^{-6}$, stop successfully, including at iteration zero or after the last permitted accepted step.
3. Check finite values and descent. If the accepted-iteration budget has been reached without convergence, return iteration-limit failure.
4. Start the trial step $t$ at 1. Build each trial by copying the current control and using the native scaled-add operation with $td$.
5. Solve only the trial state and evaluate only its objective. Compute the actual update $`\delta u=u_{\mathrm{trial}}-u`$ by subtracting the current control from the computed trial control, then compute $s=\langle j'(u),\delta u\rangle$.
6. Accept exactly when $s<0$, the trial objective is finite, and $`J_{\mathrm{trial}}\leq J_{\mathrm{current}}+10^{-4}s`$. A nonfinite slope is a failure. A nonfinite objective is recorded and rejected, as in current nmopt. Native solve failures abort rather than becoming objective infinities.
7. On rejection, halve the step and try again, up to 30 total trials. No adjoint, objective derivative, or metric gradient is evaluated at a rejected trial. Each retry is formed from the original current control, not from the preceding rejected trial.
8. On acceptance, retain the already-solved trial state and objective, compute the objective derivatives and adjoint once, and form the reduced derivative. Do not solve the state again. Increment accepted iterations and return to the stopping check.
9. If no trial is accepted, report line-search failure. Do not alter stopping criteria or invoke another method.

The actual-update Armijo slope matches `reduced_line_search.hpp`. Preserve this operation order; `dot` in the following pseudocode denotes the native vector dot product:

```text
actual_update = trial_control - current_control
slope = dot(reduced_derivative, actual_update)
```

Do not substitute `step * dot(reduced_derivative, direction)` for the slope calculation. Although algebraically equivalent, its rounding can change acceptance near a threshold. The native loop need not reproduce nmopt's generic policy classes, layout checks, optional histories, or unnecessary copies. Record the extra services/work instead.

## 7. Instrumentation and attribution

Use plain externally owned counters and compact event records. Keep instrumentation out of adapted Step-4 and out of shared nmopt. Native mathematical-operation entries and nmopt callback entries have separate counters; do not sum both as though they were two mathematical evaluations.

Instrument these observable operations:

- Assembly; state solve; adjoint solve; objective; objective state/control derivatives.
- Residual, JVP, state VJP, control VJP, and nmopt full-VJP callback.
- Explicit matrix `vmult` and `Tvmult` calls outside CG, including the otherwise unused state-VJP calculation.
- Metric `apply` and `inverse_apply` callbacks.
- Trial evaluation, accepted derivative augmentation, and field output.
- CG calls, iterations, monitored residuals, and failures, distinguished by state/adjoint role.

Use phases `setup`, `runtime`, `verification`, and `output`; add path `native` or `nmopt`, scenario ID, control/evaluation ID, accepted-iteration index, and trial index where applicable. Reset/snapshot counters at phase boundaries without resetting or mutating numerical state.

**Derived schedule invariant:** one successful `evaluate_value` plus one successful `augment_derivative` in a fresh counter scope implies one state solve, one objective evaluation, one adjoint solve, and no repeated state solve. On the nmopt side it implies zero residual/JVP runtime callbacks and one full-VJP callback. The native reference needs only the control pullback. The full VJP must really compute its state component; count that matrix-transpose action.

**Derived schedule invariants:** for successful optimization with $K$ accepted steps and $R$ rejected trials, before output/audits:

| Operation | Count |
| --- | --- |
| State solves | $1+K+R$ |
| Adjoint solves | $1+K$ |
| Objective values | $1+K+R$ |
| Objective-derivative evaluations | $1+K$ |
| nmopt full-VJP callbacks | $1+K$ |
| nmopt residual callbacks | $0$ |
| nmopt JVP callbacks | $0$ |

These identities must be tested against observed results. Repair local implementation bugs that violate the prescribed schedule, retaining the failing evidence. Never manipulate counters or add/remove work merely to reproduce the expected counts, and never change shared nmopt to make the identities hold. Report actual observed counts and explain discrepancies, including any unresolved difference taken to G1.

If objective state/control partials are separate native calls, record each as $1+K$ rather than miscounting them as duplicate full derivative evaluations. Record the exact grouping used.

Compare measured external counts with nmopt's reported counts. `metric_solve_count` counts inverse-metric requests for this solver; an identity inverse has no iterative solve. Report it as such. Metric applications used for nmopt's gradient/step diagnostics can legitimately exceed native work and should be measured rather than artificially added to the native loop.

Instrumentation limits are part of the result:

- Ordinary wrappers cannot observe every copy, allocation, or operator call inside unmodified nmopt/deal.II.
- Count explicit experiment-owned copies at named sites. Audit core copies at fixed source locations and label their counts source-derived, with assumptions about calls and copy elision. Never label such totals measured.
- Record copied vector lengths and coefficient-volume estimates separately from actual allocations. An allocation count is unknown unless directly measured.
- CG iteration counts are measured. Do not relabel them as exact internal matrix-application counts.
- Do not replace `SerialBackend`, vector storage, CG matrix types, or global allocation to make instrumentation easier in the primary experiment. Optional profiling is a later diagnostic if needed.

### Minimal working ledger

Maintain one CSV, one row per logical operation/seam/conversion, not per source line. Initial location: `runs/external-dealii/step-4/working/attribution.csv`. Never delete that working directory as part of a test cleanup.

```text
id,source_site,operation,first_requiring_consumer,other_consumers,native_runtime_required,nmopt_construction_required,verification_required,cost_category,phase,frequency,explicit_copy_volume,operator_or_solve_work,ownership_or_lifetime_obligation,evidence_kind,evidence_ref,notes
```

Allowed evidence kinds: `measured`, `source-derived`, `estimated`, `unknown`.

Use `first_requiring_consumer` as a causal explanation, not simply the earliest unit to implement the operation. The independent requirement flags are essential: verification introduced in E3 must not conceal that nmopt construction would still demand the operation if verification were absent.

Suggested initial rows:

| ID | Operation | Native runtime | nmopt construction | Verification | Attribution |
| --- | --- | --- | --- | --- | --- |
| S01 | Separate prepare from run | yes | indirectly | yes | Application reuse |
| S02 | Solve supplied RHS | yes | yes | yes | Application reuse; reused numerical machinery |
| M01 | Objective and partial derivatives | yes | yes | yes | Added control mathematics |
| M02 | Control VJP $-q$ | yes | yes | yes | Added control mathematics |
| V01 | Residual evaluation | no | yes | yes | Shared verification plus frozen API obligation |
| V02 | Residual JVP | no | yes | yes | Shared verification plus frozen API obligation |
| N01 | State component of full VJP | no | yes | yes | Additional frozen API capability/runtime work |
| N02 | Full-variable block composition | no | yes | no | nmopt representation and validation |
| N03 | Identity metric contract wrapper | no wrapper needed | yes | optional check | nmopt mechanical adaptation |
| O01 | Value-only rejected trial handling | yes | supplied by solver | checked | Optimizer orchestration service |

Also keep small CSV traces for solves, evaluations, and trials. Fixed headers and scalar values are sufficient; no general serializer is required. Every run records both baseline hashes, experiment source revision and dirty-state description, scenario constants, compiler/deal.II versions, profile, and commands. Keep a final reviewed ledger/report in `docs/planning/review/external-dealii-boundary-evaluation.md` at G1, or attach a small CSV beside it. Do not promote speculative diagnoses or raw field files merely to complete the report.

Add a small source-size summary to the final report using one stated counting convention, such as nonblank, noncomment physical source lines. Separate original/derived baseline source, actual reuse delta, native mathematical operations, native orchestration, binding, and verification. Use disjoint source spans for additive totals; shared operations are counted once with multiple consumer flags. Do not sum overlapping ledger symbols, count an included template definition once per translation unit, or compare a stripped file's total size with an adapter's incremental size. Source-size measurements supplement the causal ledger; they are not acceptance thresholds.

## 8. Verification definitions and tolerances

**Frozen protocol; not yet validated experimentally.** All tolerances below are acceptance thresholds for this fixed 2D problem, not observed successes. Log actual errors. Do not relax a threshold to make an unexplained failure pass.

For paired vectors and paired scalar objectives, respectively, use:

```math
\begin{aligned}
\lVert a-b\rVert_{2}&\leq 10^{-11}+10^{-10}\max(\lVert a\rVert_{2},\lVert b\rVert_{2}),\\
|a-b|&\leq 10^{-12}+10^{-11}\max(|a|,|b|).
\end{aligned}
```

Matrix/RHS assembly should be identical in the paired environment; investigate any difference before reduced evaluation. For symmetry use:

```math
\frac{\lVert A-A^{\mathsf T}\rVert_{\mathrm F}}{\max(1,\lVert A\rVert_{\mathrm F})}\leq 10^{-14}.
```

Require successful CG under its unchanged `1e-12` monitored-residual policy. Independently check:

```math
\frac{\lVert Ax-\mathrm{rhs}\rVert_{2}}{\max(1,\lVert\mathrm{rhs}\rVert_{2})}\leq 10^{-10}.
```

This is a verification tolerance, not a replacement solver stopping rule. Log both residuals and their definitions.

For $n=289$, define native vectors by index $i=0,\ldots,n-1$, avoiding platform-dependent random generation:

```math
\begin{aligned}
c_{i}&=\frac{1}{\sqrt{n}},\\
r_{i}&=\frac{i/(n-1)-1/2}{\sqrt{n}},\\
a_{i}&=\begin{cases}
-1/\sqrt{n},&i\text{ odd},\\
1/\sqrt{n},&i\text{ even}.
\end{cases}
\end{aligned}
```

Reduced evaluation controls: $0$, $0.1c$, $0.1r$, $0.1a$. Include a fresh repeated evaluation of $0.1r$ after another control to detect leaked mutable solve state. Normalize directions $r$ and $a$ to unit Euclidean norm for reduced derivative checks. Test an off-solution full point with nonzero state/control and independent nonzero state/control tangents and test seed.

Checks:

- Residual JVP centered finite difference and full JVP/VJP dual pairing, including separate state-only and control-only tangents. Relative/scaled error target $10^{-8}$ for finite differences, $10^{-12}$ for direct linear pairing.
- Objective directional derivative at an off-solution point, target $10^{-8}$ scaled error.
- Reduced centered finite differences with states recomputed at every perturbation, using $h=10^{-2},10^{-3},10^{-4},10^{-5},10^{-6}$. Require agreement within $10^{-7}\max(1,|\langle j'(u),v\rangle|)$ for at least two adjacent usable steps, where $v$ is the normalized direction. Save the full table, including the small-step noise floor.
- Reduced first-order Taylor remainder at $h=0.1,0.05,0.025$. Require positive remainders above ten times the observed repeat-evaluation objective variation and consecutive halving ratios between 3.5 and 4.5. If numerical noise defeats this fixed range, report and review it; do not silently filter inconvenient points.
- Dense oracle system relative residual and independent stationarity residual at most $10^{-10}$ under their recorded normalizations. Compare the native reduced gradient evaluated at oracle control against zero with norm at most $10^{-8}$.
- For optimization, require both paths to stop by gradient norm and independently audit final gradients at most $1.1\cdot10^{-6}$. Require control distance to oracle at most $2\cdot10^{-6}$, contingent on the oracle checks passing. Retain objective gaps as evidence as well.

Paired evaluation should pass at every chosen control before E5. During E5, initially expect matching trial/accepted-step sequences in the fixed environment, but bitwise equality is not the mathematical acceptance criterion. On a divergence, stop automatic acceptance, preserve the earliest differing controls/values and the signed Armijo margin, and explain whether the cause is arithmetic, solve accuracy, different schedules, or a defect. Continue only after that explanation is accepted or a local experiment bug is repaired. Never call unexplained trajectory differences architectural evidence.

## 9. Work units, artifacts, and gates

### E0 — Protocol and documentation

Outcome: one current execution authority with all choices above explicit.

Files: `docs/planning/external-dealii-boundary-evaluation.md`, `docs/README.md`, `docs/planning/external-dealii-tutorial-roadmap.md`, and `docs/reference/external-dealii-solver-integration.md`. No C++, CMake, test, or design changes.

Checks: source-check every factual correction; inspect Markdown/math/links; confirm that historical notes and speculative interfaces do not become implementation instructions. Confirm baseline refs and branch-start instructions.

Artifacts: adopted protocol, documented constants, attribution CSV schema/sample rows, status identifying E1 as next. E0 should not assert that future tests pass.

Prospective commit: `docs(dealii): define external boundary evaluation protocol`.

Gate: E0 is adopted only after its documentation/source verification and the authorized documentation commit. Until then its status is E0 ready for review. E1 begins only from that coherent adopted protocol. Branch creation/switching and the E0 commit were explicitly authorized for this unit; future actions remain subject to the user's scope and repository instructions.

### E1 — Comments-stripped baseline

Outcome: a reproducible stripped numerical reference independent of adaptation.

Files: new baseline, bounded stripping/forward-comparison scripts, stripped standalone target and forward-comparison registration, README/provenance evidence.

The stripper must preserve literals, token boundaries, preprocessor directives, and legal notice. Compare non-comment token sequences; compilation/output agreement alone is insufficient proof of a comment-only transformation. Do not use a regular expression that can consume comments inside literals. Regeneration for verification writes to a temporary file and compares it to the tracked baseline. Follow repository rules for generators/rewrites; E1 creates a new file and must not bulk-rewrite existing sources.

Checks: upstream hash remains `be9e694f5f3c9177b7cd18200ff8173337c2b16e1ee72d45ab1ba7c6e105be5f`; token equivalence; original 2D and 3D counts, output, and numerical VTK content. Recorded prior counts are 289/4913 DoFs and 26/30 CG iterations, but rerun both sources in the same current environment before asserting equivalence. Upstream and stripped runs use distinct working directories.

Artifacts: source hashes, transform command/script, token-check result, forward logs and output comparisons. Parsed numeric VTK comparison must cover field values and mesh/connectivity, not file existence alone.

Prospective commit: `test(dealii): establish stripped Step-4 baseline`.

Gate: freeze the stripped source only after fidelity passes. No adaptation before this gate.

### E2 — Minimal reusable Step-4

Outcome: reusable original numerical machinery with no optimal-control or nmopt concepts.

Files: adapted `step-4.cc`, adapted standalone target, initial native reuse contract test, README, removal of obsolete wrapper/smoke files and registration once replacement checks exist.

Seams: preparation; const assembled matrix/RHS views; narrow supplied-RHS/in-out-vector solve with native evidence; supplied-state VTK writer; main guard. Prefer the existing matrix/RHS view names. Do not expose mesh/FE/DoF objects merely because Problem B may need them later. No retained optimization-specific wrapper compatibility.

Checks: upstream/stripped/adapted standalone equivalence in both dimensions, including original initialization and stdout; matrix symmetry in 2D; supplied-RHS solves; original RHS and assembled matrix unchanged by those solves; zero-RHS zero-initial solve succeeds without forced iterations; zero-control zero-initial solution matches standalone within tolerance; supplied-state writer uses original field identity/format. Native target remains independent of nmopt library headers/targets.

Artifacts: stripped-to-adapted diff, seam ledger, native solve evidence, standalone and reusable-path output comparisons. Record any iteration-count difference caused by the explicitly different controlled initialization.

Commits: E2.a was adopted as `f7e44e7` (`refactor(dealii): expose reusable
Step-4 operations`); E2.b is pending review with the prospective message
`refactor(dealii): retire obsolete Step-4 wrapper`.

Gate: all reuse/fidelity checks pass; otherwise fix adaptation before adding control mathematics.

E2.b evidence (2026-09-10): with `HEAD` `f7e44e7` plus the uncommitted
retirement changes, the superseded `tutorial_application.hpp`,
`tutorial_application.cc`, and `tutorial_binding_smoke.cc` files were deleted.
Their prior implementation remains recoverable in Git from the historical
tutorial attempt. The old target and CTest registration had already been
removed in E2.a; repository search now finds no active source or build
reference to the wrapper or smoke target, while historical planning documents
continue to mention them intentionally. The complete Debug deal.II pipeline
passed all `163/163` tests:

```text
./build.sh pipeline debug-dealii
100% tests passed, 0 tests failed out of 163
```

E2 is ready to close after review and the prospective cleanup commit. E3 may
start only after that commit is adopted.

### E3 — Native Problem A

Outcome: verified native problem and reduced evaluator, including separate value and derivative stages.

Files: scenario, Problem A, native reduced evaluation, minimal instrumentation, native verification/oracle, native test scenarios. No nmopt numerical dependency, new public metric, or optimizer yet.

Checks: frozen equations at arbitrary/off-solution inputs, state and adjoint residuals, symmetry assumption, deterministic sample evaluations, repeated-control independence, derivative/pairing/Taylor checks, oracle residual/stationarity, native output from retained state, and exact state-reuse solve counts.

Artifacts: native operation inventory, required capability flags including future frozen-nmopt obligations, sample evaluation/solve tables, derivative/Taylor table, oracle errors, updated ledger.

Prospective commit: `test(dealii): add native algebraic control reference`.

Gate: native mathematical proof passes before E4. Do not implement a second native formula to make nmopt agree. This unit is cohesive; split oracle-only verification into a follow-up commit only if size warrants it, and keep the E3 gate closed until both are complete.

### E4 — Current nmopt Problem A

Outcome: complete valid construction through the frozen public API and explained agreement at the reduced-evaluation level.

Files: explicit nmopt binding, local identity metric, callback/report instrumentation, paired-evaluation test driver, target registration. No changes to native problem equations or Step-4 seams merely for convenience; any new required seam must be explicit, causally attributed, and separately reverified.

Use one state block, one control block, one test block, distinct space identities, and actual dimensions obtained from the native problem. Supply residual, JVP, full VJP, objective, and objective derivative. Do not use throwing placeholders, fabricated derivative blocks, a dummy executable, or callbacks that bypass `ReducedDTOT` by calling the native reduced evaluator.

Checks: each callback checked against native operations; off-solution full VJP; correct report translation; all frozen sample controls and repeated-control case; staged evaluation counts; metric gradient and negative search direction; valid borrowed lifetimes; no compiler/private-header dependence. Verify that distinct compatible layout objects are not assumed incompatible.

Artifacts: paired result/error table, construction capability inventory, runtime callback/native-operation counts, source-derived copy audit, header/link dependency evidence, updated ledger.

Prospective commit: `test(dealii): compare native and nmopt reduced evaluations`.

Gate: E5 cannot begin with unexplained mismatches. A genuine frozen-API obstruction produces a minimal reproducer and an early G1 diagnostic report; it does not authorize modifying the API.

### E5 — Matched optimization

Outcome: the bounded native loop and existing nmopt steepest-descent/Armijo solver solve the same problem under the frozen selected policy.

Files: native Armijo loop, native and paired optimization test scenarios, scalar trace/evidence output. Reuse E3 mathematics and E4 binding unchanged except for diagnosed local defects, which require repeating affected evaluation checks.

Checks: exact selected policy and stopping precedence, value-only rejection, accepted-state reuse, observed counts against the derived schedule invariants, both stopping reasons, final independently audited optimality/oracle agreement, recorded CG work, and final output from retained states. Compare source-level copies/operator actions without inflating native work to imitate unnecessary framework work.

Artifacts: both traces, first-divergence evidence if any, solve/iteration totals, final errors/gaps, output comparison, updated service/attribution ledger. Do not use only objective reduction as a convergence criterion.

Prospective commit: `test(dealii): compare native and nmopt optimization paths`.

Gate: every numerical or schedule discrepancy must be resolved or explicitly accepted as explained numerical variation before a successful comparison is claimed. If the selected fixed constants fail, preserve the failure; propose a symmetric protocol amendment instead of tuning nmopt and native independently.

### G1 — Attribution and decision

Outcome: reviewed factual account and a bounded next decision.

Files: final report under `docs/planning/review/external-dealii-boundary-evaluation.md`, compact reviewed attribution if separate, this roadmap's current status, README/reference links where factual outcomes warrant them. No library or design refactor in this unit.

Report mathematical success or failure, fidelity, extra required capabilities, one-time/repeated costs, library-provided services, source-derived versus measured evidence, and limits. Distinguish incidence in native versus nmopt paths from universal claims about all PDE-control applications. Do not force every row into an exclusive category when it has multiple consumers.

Classify findings as mechanical ergonomics, capability/formulation obligation, structural cross-layer coupling, small justified overhead, or unresolved. A local programming defect or numerical-policy failure is not itself a boundary diagnosis. Recommend a helper only if explicit repeated mechanical construction motivates it; recommend boundary investigation only with a concrete unnecessary obligation or correctness obstruction that a helper cannot remove.

Prospective commit: `docs(dealii): report external boundary evaluation findings`.

Gate: stop implementation after presenting the recommendation. Helper/API changes and Problem B require their own accepted scope. When possible, prefer running Problem B against the same frozen baseline before changing the boundary, so its incremental comparison remains interpretable.

## 10. Verification execution and artifacts

Follow `.agents/build.md` for required pipelines. In particular, use existing `build.sh` profiles and the machine's existing `build.local.conf`; do not invent configuration or job limits. Documentation-only E0/G1 needs no build unless new code changes require it. C++/CMake units run focused new scenarios and the required Debug pipelines; reuse valid results within a unit and rerun affected checks after fixes.

Use the existing helper for configuring/building/testing. If a manual deal.II build is necessary, it must use `--parallel 1`. Do not configure or rebuild release without the explicit authorization required by build instructions. This phase uses work counts, not timing claims, so release timing is not a gate. Do not create sanitizer-deal.II profiles or weaken warning settings.

Choose explicit finite CTest timeouts appropriate to forward/native checks and the bounded optimization; an infrastructure timeout is a failed run, not architectural evidence. New scenarios use existing discovery. Add only direct target-specific construction necessary to keep native executables independent.

Ignored run-artifact root:

```text
runs/external-dealii/step-4/
  working/attribution.csv
  forward-comparison/<run-id>/{upstream,stripped,adapted,comparison}/
  reduced-evaluation/<run-id>/{native,nmopt,comparison}/
  optimization/<run-id>/{native,nmopt,comparison}/
```

Use isolated directories and unique run IDs so retries do not destroy failing
traces or overwrite native with nmopt output. Directory names describe the
comparison or scenario, not the roadmap unit that produced them. This is a
simple artifact convention, not a new runner schema. Record commands,
revisions, constants, compiler/deal.II/build information, test outcomes,
source hashes, and numerical errors in plain text/CSV. Evidence files are
outputs, not mutable configuration.

At each unit handoff state: completed artifacts and checks, actual failures, protocol deviations, updated attribution, exact next unit, and prospective commit boundary. Do not mark a unit complete because it compiles or emits files. Final scope check must show no shared nmopt implementation changes and no edits to upstream/frozen stripped source.

## 11. Hard stops and unresolved questions

Stop dependent work and preserve evidence on:

- A changed/unreconciled baseline or conflicting authoritative instructions.
- Upstream token/provenance loss, standalone numerical-policy drift, or unexplained output/iteration differences between equivalent standalone paths.
- Failure of native mathematical verification or the oracle.
- Unexplained paired-evaluation or optimization divergence.
- Need to change frozen shared code, replace a solve policy, import compiler/private interfaces, or supply invalid callbacks.
- Inability to observe a requested cost honestly: record it as unknown rather than inventing a number; only stop the whole unit if that observation is essential to its gate.
- Scope growth into a generic framework, unrelated cleanup, additional algorithms, or new PDE families.

The implementation model may diagnose and repair in-scope local bugs without requesting a new architectural decision for every edit. It must not silently alter mathematical conventions, acceptance tolerances, numerical constants, or the shared-code freeze. Proposed amendments must include the observed failure and apply equally to both paths.

**Frozen protocol:** full-vector Problem A first; 2D control experiments; both-dimensional standalone fidelity; exact signs/objective/metric; original versus controlled initial guesses; narrow in/out solve seam; exception behavior/report semantics; Armijo expression and stopping precedence; no cross-control caching; explicit instrument limitations; working versus promoted evidence; treatment of obsolete wrapper files.

**Open hypotheses through G1:** whether construction obligations are materially costly; whether mechanical helpers would suffice; whether unused full-VJP work matters; whether source-derived copy estimates justify profiling; whether a local boundary issue or deeper structural coupling exists; which boundary change, if any, should be investigated; and how much of Problem A generalizes.

Problem B's control space, coupling assembly, $L^{2}$ operators, free-coordinate/lifting policy, and metric realization are later decisions. A nonsymmetric probe, rectangular injection, alternate metric, compiler-versus-external equivalence, and package/install evaluation are optional follow-on experiments tied to specific unanswered questions. Do not implement them during E0–G1.
