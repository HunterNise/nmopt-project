# Step-4 external integration: implementation report

## Findings

This report records the final implementation state of the external deal.II Step-4
experiment at `170c9f1` (`refactor(step4): isolate adapted application reuse`). It is
an implementation/accounting report, not a claim that the measured source sizes are
universal lower bounds.

The experiment was motivated by a repository-level concern: a predominantly
Codex-authored codebase had grown large, so it was important to determine whether
that size reflected an intrinsically heavyweight external-application boundary or a
large amount of evaluation machinery around a smaller functional path.

The final source attribution separates those two things clearly:

| Scope | Code-bearing lines | Needed by the canonical minimal consumers? |
| --- | ---: | --- |
| Functional Step-4/OCP/nmopt path, using the adapted source once | **1,722** | Yes, by role |
| Of that: minimal nmopt binding headers | **426** | Yes |
| Of that: minimal executable entry points | **219** | Yes, but largely ordinary program/policy code |
| Evaluation, comparison, verification, diagnostics, and test/evidence support | **10,256** | **No** |

The large experiment footprint therefore should not be read as the amount of code an
external application must adopt. The evaluation layers exist to test fidelity,
mathematical correctness, native/nmopt equivalence, attribution, and output.

The tested boundary is nevertheless not effortless. The minimal bindings still
require several public concepts: layouts, five executable-model callbacks, two solve
services with truthful reports, a metric adapter, a state/control partition, and a
`ReducedDTOT`. The experiment establishes that these obligations are localized and
stable across the tested A/B cases; it does not establish globally minimal authoring
cost.

## 1. Source roles and ownership

The final experiment layout is:

```text
source/
  upstream/       authentic deal.II provenance source
  baseline/       comment-stripped numerical baseline
  adapted/        reusable Step-4 source; still a standalone program

integration/
  adapted_step4.hpp       one private include seam for the adapted source
  problem_a.hpp           Problem A OCP operations and native solves
  problem_b.hpp           Problem B OCP operations and native solves
  problem_b_coordinates.hpp
  problem_b_mass.hpp
  problem_b_metric.hpp    Problem B coordinates, mass/coupling, and metric

minimal/
  problem_a_binding.hpp
  problem_b_binding.hpp   canonical nmopt-facing adapters
  problem_a.cc
  problem_b.cc            complete minimal consumers

evaluation/
  nmopt_problem_*_binding.hpp
  native_*                instrumented nmopt bindings and independent native paths

verification/             equation/derivative/oracle checks
diagnostics/              optional counters and evidence records
```

The important ownership rule is that the application remains the numerical owner.
`source/adapted/step-4.cc` owns its mesh, FE space, matrix, RHS, native CG solve, and
writer. The OCP layer owns the control model, objective, derivatives, coordinate
maps, and metric. The nmopt binding owns only contract adaptation and reduced-service
composition.

The optional diagnostics used by the experiment are source-level dependencies of
some OCP/evaluation headers, but the minimal consumers instantiate no
`Instrumentation`, native reference optimizer, dense oracle, or comparison runner.

## 2. The starting Step-4 application

The experiment starts from the authentic deal.II `step-4` tutorial rather than from
a framework-shaped mock application. In its original form, Step-4 is a small
finite-element Poisson solver. It builds a hypercube mesh, globally refines it, uses
continuous first-order `FE_Q` elements, assembles the stiffness matrix and load
vector, applies prescribed Dirichlet data, solves the resulting linear system with
conjugate gradients and identity preconditioning, and writes the finite-element
field to VTK.

At the continuous level, the tutorial has the familiar form

```math
\begin{aligned}
-\Delta y &= f && \text{in } \Omega, \\
y &= g && \text{on } \partial\Omega.
\end{aligned}
```

The tutorial itself chooses the forcing and boundary values. For this integration
study, the important fact is not their particular formula but the ownership model:
Step-4 already has a complete mesh/assembly/solve/output lifecycle before `nmopt` is
introduced.

With the continuous first-order `FE_Q` space $V_{h}$ and its homogeneous-test
subspace $V_{h,0}$, the finite-element problem is: \
find $y_{h}\in V_{h}$ with the prescribed discrete boundary values such that

```math
\int_{\Omega} \nabla y_{h} \cdot \nabla v_{h} \mathrm{d}x
=
\int_{\Omega} f v_{h} \mathrm{d}x
\qquad
\text{for all } v_{h}\in V_{h,0}.
```

After assembly and Step-4's existing Dirichlet treatment, this becomes the algebraic
system

```math
A y = b.
```

The symbols $A$ and $b$ used below refer to these already assembled and
boundary-treated Step-4 objects. The OCP experiment does not replace this
discretization: Problem A modifies its algebraic right-hand side, while Problem B
builds a free-state coordinate system and finite-element control coupling around the
same native realization.

Its original public interface is essentially one monolithic `run()` operation:

```text
make grid
  -> distribute DoFs / allocate system
  -> assemble PDE
  -> solve native linear system
  -> write VTK output
```

The standalone `main()` executes that forward problem in both 2D and 3D. In the 2D
configuration used by Problems A and B, four global refinements produce 256 active
cells and 289 degrees of freedom. Problem B later distinguishes the 225 free state
coefficients from the full 289-entry physical/control representation, but that
coordinate split belongs to the OCP layer rather than to the original tutorial.

This starting point matters for interpreting the experiment. Step-4 is a forward PDE
application, not an optimal-control application and not an `nmopt` client. Turning it
into an OCP necessarily requires application mathematics such as a control,
objective, derivatives, adjoint, and metric. Those additions should not be counted
as framework binding overhead. The first question is therefore only whether the
existing numerical application can be made reusable without being rewritten around
`nmopt`.

## 3. Step-4 before and after adaptation

The upstream and stripped baseline both contain **191 code-bearing lines**. The
adapted tutorial contains **244**, a net increase of **53**.

The adaptation makes an originally monolithic forward example reusable without
turning it into an nmopt application. It exposes:

| Reuse seam | Purpose |
| --- | --- |
| preparation/assembly entry | build the mesh/system once before repeated solves |
| matrix and RHS views | allow application-owned OCP operators to reuse the assembled system |
| `solve(rhs, solution)` | reuse Step-4's native CG policy for caller-supplied right-hand sides and return monitored evidence |
| `output_results(state, path)` | reuse the original VTK writer for a caller-supplied state |
| DoF-handler and boundary-data views | let Problem B construct its coordinate map and FE mass/coupling |
| `STEP4_NO_MAIN` guard | permit the preserved translation unit to be reused without its standalone `main()` |

The preprocessor protocol is centralized in
[`integration/adapted_step4.hpp`](integration/adapted_step4.hpp), a four-line private
reuse header. Consumer code includes that named seam rather than repeating a
`#define`/include/`#undef` sequence.

No nmopt header, type, callback, metric, or solver enters
`source/adapted/step-4.cc`. Its original forward `run()` remains available, and the
standalone adapted target remains a direct source build.

## 4. Application-owned OCP mathematics

The OCP implementation is deliberately separated from nmopt-specific binding code.
These mathematical choices would still have to exist if a different optimization
library were used.

### 4.1 Problem A – algebraic RHS control

[`integration/problem_a.hpp`](integration/problem_a.hpp) contains **214 code-bearing
lines** and owns the prepared Step-4 application plus the A-specific residual,
objective, derivatives, state/adjoint solves, and output forwarding.

Problem A keeps the same assembled matrix $A$ and turns the existing algebraic
right-hand side into the control channel. Its primal equation, adjoint, objective,
and reduced derivative are

```math
\left\{
\begin{aligned}
Ay &= b+u,\\
A^{\mathsf T}p &= y
\end{aligned}
\right.
\qquad
\begin{aligned}
J(y,u) &= \frac{1}{2}y^{\mathsf T}y
        + \frac{1}{2}u^{\mathsf T}u,\\
j'(u) &= u+p
\end{aligned}
```

For the public model callback, the same state equation is
represented by the residual $E(y,u)=Ay-b-u$.

State and control are both 289-entry algebraic vectors. A is intentionally simple
and does not represent physical FE volume forcing with preserved state boundary
values.

### 4.2 Problem B – finite-element distributed control

Problem B's native OCP implementation contains **615 code-bearing lines**:

| Source | Lines | Responsibility |
| --- | ---: | --- |
| `integration/problem_b.hpp` | 233 | residual, objective, derivatives, state/adjoint solve wrapping |
| `integration/problem_b_coordinates.hpp` | 118 | free/full state map, boundary lifting, restriction/reconstruction |
| `integration/problem_b_mass.hpp` | 183 | consistent FE mass and rectangular control coupling |
| `integration/problem_b_metric.hpp` | 81 | native mass metric apply/inverse service |
| **Total** | **615** | application-owned B OCP implementation |

The full control has 289 continuous FE coefficients while the state/test spaces use
225 free coordinates. With $P$ the free-to-full injection, $\ell$ the fixed boundary
lifting, and $M$ the full consistent mass matrix, define

```math
\begin{aligned}
y_{\mathrm{phys}} &= Pz+\ell, \\
K &= P^{\mathsf T}AP, \quad
b_{F} = P^{\mathsf T}b, \quad
B = P^{\mathsf T}M.
\end{aligned}
```

The resulting primal equation, adjoint, objective,
and reduced derivative are

```math
\left\{
\begin{aligned}
Kz &= b_F+Bu, \\
K^{\mathsf T}p &= P^{\mathsf T}M y_{\mathrm{phys}}
\end{aligned}
\right.
\qquad
\begin{aligned}
J(z,u) &= \frac{1}{2}y_{\mathrm{phys}}^{\mathsf T}
          M y_{\mathrm{phys}}
          + \frac{1}{2}u^{\mathsf T}Mu,\\
j'(u) &= Mu+B^{\mathsf T}p
\end{aligned}
```

For the public model callback, the primal equation is equivalently represented by
$E(z,u)=Kz-b_{F}-Bu$. The application owns both the OCP operators and the mass
metric $G=M$. State and adjoint solves reuse the native Step-4 linear solver;
physical output reconstructs $y_{\mathrm{phys}}$ before calling Step-4's writer.

The growth from A to B is therefore mainly OCP/numerical work: coordinate semantics,
FE coupling, a physical objective, and a nonidentity metric. It is not growth in the
nmopt optimizer or compiler.

## 5. Mathematical operations and the public connection

The canonical minimal A and B bindings expose the same public shape.

| Mathematical operation | Native owner | nmopt representation |
| --- | --- | --- |
| residual $E$ | `ProblemA` / `ProblemB` | test-layout covector callback |
| residual JVP $E'\delta x$ | native problem | test-layout covector callback |
| residual VJP $E'^{\ast}p$ | native problem | full state/control covector callback |
| objective $J$ | native problem | scalar callback |
| objective derivative $J'$ | native problem | full state/control covector callback |
| state solve | native problem/Step-4 | `StateAdjointSolversT::solve_state` |
| adjoint solve | native problem/Step-4 | `StateAdjointSolversT::solve_adjoint` |
| control metric | identity for A; native mass metric for B | local `MetricT` adapter |
| reduced derivative | composed from the operations above | `ReducedDTOT` |

The binding wraps/unwraps native `dealii::Vector<double>` values while preserving
state/control/test layouts and primal/covector roles. It does not reimplement the
PDE algebra.

A successful first-order reduced evaluation uses:

```text
ReducedSearchSolverT::solve(control)
  ReducedDTOT::evaluate_value(control)
    native state solve
    native objective
  ReducedDTOT::augment_derivative(retained value)
    native objective derivative
    native adjoint solve
    residual VJP
    reduced control covector
  MetricT::inverse_apply
  line search / stopping / retained result
```

This separation is also where current API limitations are visible: all five model
callbacks are mandatory at construction even though this successful first-order
runtime does not call residual or JVP; the full VJP computes a state component that
the reduced path later discards when extracting the control contribution; and
`MetricT::inverse_apply` does not expose a public metric-solve report.

## 6. Source accounting

### 6.1 Counting method

The experiment's authoritative accounting counts physical C++ lines that remain
after blank and comment-only lines (including leading block-comment licenses) are
removed. Includes, declarations, and brace-only lines count. These numbers describe
source ownership; they are not development-time estimates or universal minima.

The conventional gross size of the repository is less informative here because the
same experiment intentionally contains both the consumer path and a much larger,
independent validation harness.

### 6.2 Functional source

All paths below are relative to `apps/external-dealii/step-4/` unless noted.

| Responsibility | A | B | Shared | Minimal runtime? |
| --- | ---: | ---: | ---: | --- |
| Preserved stripped forward baseline | – | – | 191 | Existing application |
| Net native reuse adaptation | – | – | 53 | Yes, for this originally monolithic tutorial |
| Adapted-source reuse helper | – | – | 4 | Yes, experiment-local include seam |
| OCP operations/native solve wrapping | 214 | 233 | – | Yes, application/OCP work |
| Free/full state coordinates and lifting | – | 118 | – | Yes, B OCP work |
| FE mass assembly and control coupling | – | 183 | – | Yes, B OCP work |
| Native mass metric | – | 81 | – | Yes, B OCP work |
| Minimal nmopt binding | 206 | 220 | – | **Yes, framework-specific** |
| Consumer entry point | 107 | 112 | – | Yes, mixed nmopt use + ordinary program policy |

Counting the complete adapted source once rather than the baseline plus delta gives a
**functional gross of 1,722 code-bearing lines**:

```text
adapted Step-4 source             244
adapted-source reuse helper         4
application-owned OCP             829
minimal nmopt bindings            426
minimal consumer entry points     219
                                -----
functional gross                1,722
```

Relative to the 191-line stripped forward baseline, the experiment adds **1,531
functional code-bearing lines**. They divide more usefully by responsibility:

```text
native reuse delta + helper        57
application-owned OCP             829
minimal nmopt binding             426
mixed consumer entry points       219
                                -----
net functional addition         1,531
```

The 829 OCP lines are not nmopt boilerplate: they define the optimal-control problem
that the original forward solver did not contain. The 219 executable lines are also
not pure framework burden; they include frozen example policy, output allocation,
CLI handling, reporting, and error handling.

### 6.3 Minimal nmopt binding breakdown

The two canonical bindings total **426 code-bearing lines**. Their explicit A/B
separation is intentional: it shows that the richer B mathematics does not require a
different integration pattern.

| Binding responsibility | Problem A | Problem B | What it is for |
| --- | ---: | ---: | --- |
| Metric adapter | 51 | 58 | implement `MetricT`; A copies through identity, B forwards native mass apply/inverse |
| Layout construction | 9 | 9 | identify state, control, and residual-test spaces |
| Model callback factory | 39 | 39 | expose residual, JVP, VJP, objective, and objective derivative |
| Native solve-report translation | 15 | 15 | preserve the actual CG convergence/work evidence |
| State/adjoint solve factory | 23 | 24 | adapt native solves to `StateAdjointSolversT` |
| Partition / reduced composition | 7 | 7 | build state/control partition, metric connection, and `ReducedDTOT` |
| Packaging, aliases, accessors, storage | 62 | 68 | local C++ ownership/lifetime and readable example packaging |
| **Total** | **206** | **220** | |

The first six rows are the current contract construction. The final row is not a
claim of unavoidable framework API: it includes local class organization, aliases,
accessors, identifiers, member storage, and deleted copy/move operations.

The example does not hide repeated A/B construction behind a custom helper merely to
produce a smaller LOC number. Such a helper could shorten the files while leaving
the same public obligations in place.

### 6.4 Evaluation and evidence support

The confidence-building machinery is much larger than the consumer path and is
intentionally separate:

| Responsibility | Lines | Minimal runtime? |
| --- | ---: | --- |
| Evaluated/instrumented `nmopt` bindings | 571 | No |
| Independent native reduced/optimizer paths and frozen policy | 811 | No |
| Verification/oracles | 959 | No |
| Diagnostics | 207 | No |
| Contract/acceptance test drivers | 7,430 | No |
| Failure/evidence support | 278 | No |
| **Total evaluation/evidence overhead** | **10,256** | **No** |

The current file-level attribution is:

| Source | Lines | Why it exists |
| --- | ---: | --- |
| `evaluation/nmopt_problem_a_binding.hpp` | 235 | instrumented A binding used by the comparison harness |
| `evaluation/nmopt_problem_b_binding.hpp` | 336 | instrumented B binding with comparison/metric evidence |
| `evaluation/native_reduced.hpp` | 81 | independent native reduced evaluation for A |
| `evaluation/native_optimization.hpp` | 246 | independent native A Armijo optimization path |
| `evaluation/native_problem_b_reduced.hpp` | 145 | independent native reduced evaluation for B |
| `evaluation/native_problem_b_optimization.hpp` | 315 | independent native B optimization path |
| `evaluation/optimization_policy.hpp` | 24 | frozen policy shared by matched native/public comparisons |
| `verification/verification.hpp` | 294 | A derivative/equation checks and dense reference logic |
| `verification/problem_b_verification.hpp` | 589 | B equation, derivative, mass, boundary, and dense KKT audits |
| `verification/scenario.hpp` | 76 | deterministic verification inputs |
| `diagnostics/instrumentation.hpp` | 207 | optional counters and solve/matrix/metric records |
| `tests/dealii/external_step4_native_contract.cc` | 2,725 | native application/OCP contracts, failure probes, and operation attribution |
| `tests/application/external_step4_minimal_problem_a_contract.cc` | 452 | minimal A binding/executable fidelity against audited native behavior |
| `tests/application/external_step4_minimal_problem_b_contract.cc` | 706 | minimal B binding/executable fidelity, reconstruction, and metric checks |
| `tests/application/external_step4_nmopt_contract.cc` | 943 | evaluated A public-binding and reduced-evaluation contracts |
| `tests/application/external_step4_optimization_contract.cc` | 888 | evaluated A matched optimization and failure behavior |
| `tests/application/external_step4_problem_b_nmopt_contract.cc` | 700 | evaluated B public-binding and reduced-evaluation contracts |
| `tests/application/external_step4_problem_b_optimization_contract.cc` | 1,016 | evaluated B matched optimization, metric, and trace checks |
| `tests/dealii/external_step4_evidence.hpp` | 278 | experiment evidence writing/retention support |

This is the main answer to the repository-size concern. The experiment is expensive
to **verify**, because it maintains independent native paths, dense audits, detailed
operation attribution, failure evidence, and executable/output comparisons. That
verification cost is evidence for the library boundary; it is not a dependency of
an external application.

### 6.5 Reproduce the source counts

The role-based accounting above is more informative than a single gross C++ total.
The underlying per-file counts can be reproduced from the repository root with the
same comment-stripping utility used by the experiment:

```bash
python3 - <<'PY'
from pathlib import Path
import runpy

strip = runpy.run_path('tools/external_dealii/strip_comments.py')


def count(text):
    clean, _ = strip['strip_comment_lines'](text)
    start, license_block = strip['leading_block'](clean.splitlines(keepends=True))
    lines = clean.splitlines()[start if license_block else 0:]
    return sum(bool(line.strip()) for line in lines)


base = Path('apps/external-dealii/step-4')
files = {path for path in base.rglob('*') if path.suffix in ('.cc', '.hpp')}
files.update(Path('tests/dealii').glob('external_step4*.cc'))
files.update(Path('tests/dealii').glob('external_step4*.hpp'))
files.update(Path('tests/application').glob('external_step4*.cc'))

for path in sorted(files):
    print(f'{count(path.read_text()):5d} {path}')
PY
```

These are physical code-bearing lines after the documented transformation, not an
estimate of development effort and not a claim that every line is an unavoidable
framework obligation.

## 7. What another application would need

For an application that already has callable assembly/solve/output and already owns
the chosen OCP mathematics, the current reduced path requires the following nmopt
adaptation:

1. compatible state/control/test layouts;
2. five model callbacks;
3. state and adjoint solve services with truthful reports;
4. a control metric adapter;
5. a state/control partition and `ReducedDTOT`;
6. an initial control and selected optimization policy.

An application that does **not** already contain an OCP must additionally define its
control coordinates/coupling, objective, derivatives, adjoint, and metric. Those are
mathematical/application obligations, not evidence that nmopt rewrites the PDE.

The current reduced DTO supports one state block, one control block, and one
residual-test block. Application-specific state reconstruction and output remain
outside the generic solver.

## 8. Evidence and validation

The numerical evaluation is preserved in the historical reports, where revisions,
tolerances, operation counts, and artifact locations are frozen. The current
implementation still retains the same separation of native and public paths.

| Property | Evidence |
| --- | --- |
| Forward fidelity | upstream, stripped, and adapted programs compare original 2D/3D output |
| Mathematical correctness | off-solution residual/JVP/VJP checks, reduced derivative checks, fresh equation audits, independent dense optimum oracles |
| Native/nmopt equivalence | matched reduced evaluations and matched optimization traces for A and B |
| FE metric verification | B final stationarity checked independently with dense mass algebra |
| Consumer fidelity | actual minimal executables checked against audited native reports and VTK payloads |
| Post-refactor reproduction | fresh Problem A and B reproduction at `170c9f1` matched the historical numerical evidence |
| Current routine gate | `debug-dealii` pipeline 178/178 at `170c9f1` with deal.II build jobs = 1 |

A dedicated post-closure reproduction on 2026-09-19 preserved the old artifacts and
wrote fresh evidence to unique directories. The selected reproduction covered 11
Problem A scenarios and 14 Problem B scenarios, including native/public comparisons,
minimal executable checks, independent audits, and the recorded metric/operation
evidence.

The central numerical comparison was exact:

| Quantity | Historical | Fresh at `170c9f1` |
| --- | ---: | ---: |
| A accepted iterations | 828 | 828 |
| A line-search trials | 6,025 | 6,025 |
| A final objective | 54.376840518174902 | 54.376840518174902 |
| A final optimizer gradient norm | $9.5036543162094535\times10^{-7}$ | $9.5036543162094535\times10^{-7}$ |
| B accepted iterations | 5 | 5 |
| B line-search trials | 5 | 5 |
| B final objective | 3.4971143909160936 | 3.4971143909160936 |
| B final optimizer gradient norm | $4.8845615376102665\times10^{-8}$ | $4.8845615376102665\times10^{-8}$ |
| B independently audited mass-gradient norm | $4.8845615376171785\times10^{-8}$ | $4.8845615376171785\times10^{-8}$ |

The remaining reproduced quantities also matched: state/control dimensions, solve
counts, recomputed state and adjoint residuals, dense-oracle quantities, first
divergence, Problem B metric iteration history, operation ledgers, traces, summaries,
audits, counters, and solve records. Problem A structured evidence and the analogous
Problem B evidence were byte-identical to their historical counterparts. Forward
upstream-versus-stripped and upstream-versus-adapted comparisons matched in both 2D
and 3D; fresh VTK payloads matched historical payloads after ignoring generated
timestamp headers.

The routine 178-test gate remains separate from the reproduction result and is not
directly comparable with the historical 195/195 count because the repository's later
test policy excludes extended and reproduction labels from the routine pipeline.
The
[closure report](../../../docs/history/reviews/external-dealii-boundary-evaluation/closure-report.md)
records the original evidence and this post-closure revalidation separately.

The final code-hygiene unit also verified that:

- deal.II compilation used one build job;
- `git diff --check` passed;
- `STEP4_NO_MAIN` appears in code only in the preserved adapted guard and the private
  reuse helper; and
- the only direct source-level include of `source/adapted/step-4.cc` is from
  `integration/adapted_step4.hpp`.

The first validation build encountered one transient 30-second CMake
scenario-discovery timeout for the Problem B optimization executable. Running its
`--list-scenarios` operation directly succeeded, and the subsequent build and full
routine pipeline passed. No source/compiler failure was observed.

For historical numerical details, see the
[A report](../../../docs/history/reviews/external-dealii-boundary-evaluation/g1-report.md),
[B report](../../../docs/history/reviews/external-dealii-boundary-evaluation/problem-b-report.md),
and [minimal-consumer assessment](../../../docs/history/reviews/external-dealii-boundary-evaluation/minimal-consumers-report.md).

## 9. Reproduction

The canonical runnable consumers and focused commands are documented in
[minimal/README.md](minimal/README.md).

For source fidelity, the three standalone targets remain:

```text
nmopt_external_tutorial_step_4
nmopt_external_tutorial_step_4_stripped
nmopt_external_tutorial_step_4_adapted
```

They use deal.II directly and do not link nmopt. From the repository root:

```bash
python3 tools/external_dealii/strip_comments.py \
  --input apps/external-dealii/step-4/source/upstream/step-4.cc \
  --check apps/external-dealii/step-4/source/baseline/step-4-stripped.cc

./build.sh configure debug-dealii
./build.sh build debug-dealii --target nmopt_external_tutorial_step_4
./build.sh build debug-dealii --target nmopt_external_tutorial_step_4_stripped
./build.sh build debug-dealii --target nmopt_external_tutorial_step_4_adapted
```

Compare the upstream and stripped programs, then the upstream and adapted program:

```bash
python3 tools/external_dealii/check_forward.py \
  --upstream-executable build/debug-dealii/bin/nmopt_external_tutorial_step_4 \
  --stripped-executable build/debug-dealii/bin/nmopt_external_tutorial_step_4_stripped \
  --output-root runs/external-dealii/step-4/forward-comparison \
  --file solution-2d.vtk \
  --file solution-3d.vtk

python3 tools/external_dealii/check_forward.py \
  --upstream-executable build/debug-dealii/bin/nmopt_external_tutorial_step_4 \
  --stripped-executable build/debug-dealii/bin/nmopt_external_tutorial_step_4_adapted \
  --stripped-label adapted \
  --output-root runs/external-dealii/step-4/forward-comparison \
  --file solution-2d.vtk \
  --file solution-3d.vtk
```

The comparator checks stdout plus VTK geometry, connectivity, cell types, and numeric
arrays in separate run directories.

The complete Step-4 CTest selection remains:

```bash
ctest --test-dir build/debug-dealii --output-on-failure \
  -R '^nmopt\.(external_tutorial_step_4|external\.tutorial_step_4)\.'
```

Generated output/evidence belongs below ignored `runs/external-dealii/step-4/`
directories. Exact historical reproduction protocols and accepted numerical values
remain in `docs/history/reviews/external-dealii-boundary-evaluation/`.

## 10. Limits and conclusion

The experiment demonstrates that the tested Step-4 application can act as an
independently owned numerical producer for nmopt's reduced formulation and optimizer.
No semantic compiler, recipe, manifest, project runner, shared optimizer change, or
framework-specific PDE base class is required by the minimal path.

It does **not** establish equivalent integration cost for nonlinear/nonsymmetric
PDEs, MPI, adaptivity, constrained or boundary controls, unrelated applications,
package installation, newcomer authoring time, or performance-sensitive production
use. It also does not prove that 206/220 lines are universal lower bounds.

The observed remaining API ergonomics are concrete and bounded: mandatory model
callbacks exceed what this first-order run consumes, full VJP work includes a state
component the reduced path later discards, and metric inversion exposes no public
solve report. None blocked correctness for A or B.

The practical result is therefore narrower but useful: the very large experiment is
large mainly because it was designed to challenge and verify the boundary. The code
an external consumer actually needs is localized, inspectable, and remains nearly
the same size when moving from the simple A probe to the substantially richer FE
Problem B.
