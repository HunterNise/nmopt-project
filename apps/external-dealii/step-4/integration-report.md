# Step-4 external integration: implementation report

## Findings

The implementation preserves Step-4's native assembly, CG solve, and VTK
output while adding two application-owned optimal-control problems. Both
problems run through the existing public nmopt contracts and agree with
independently audited native references. Shared nmopt, compiler, and backend
code did not change during the evaluation.

The final minimal bindings contain 206 code-bearing lines for A and 222 for
B. The native OCP implementations contain 216 and 614 lines respectively.
These counts describe implemented responsibilities, including checks and
packaging; they are not lower bounds on mathematical or framework work.
The tables below identify every counted scope and separate the functional
consumers from the comparison machinery.

Source snapshot: `2ba749b`. The evaluation is closed. Start with the
[explanatory overview](external-integration-overview.md) for the architecture;
use the [API reference](../../../docs/reference/external-dealii-solver-integration.md)
for exact contracts and the [closure audit](../../../docs/planning/review/external-dealii-boundary-evaluation/closure-report.md)
for objective dispositions and remaining hypotheses.

## 1. Step-4 before and after adaptation

The pinned [upstream source](source/upstream/step-4.cc) is deal.II `v9.5.1`.
In the evaluated 2D configuration, the forward problem is:

```math
\begin{aligned}
-\Delta y&=4(x_{1}^{4}+x_{2}^{4}) &&\text{in }\Omega=[-1,1]^{2},\\
y&=x_{1}^{2}+x_{2}^{2} &&\text{on }\partial\Omega.
\end{aligned}
```

Step-4 uses continuous $`Q_{1}`$ elements, four global refinements, 256 cells,
and 289 DoFs. It assembles stiffness and load using its original two-point
Gauss rule in each coordinate direction. Dirichlet elimination produces the
matrix $A$ and RHS $b$. The resulting symmetric positive-definite system is
solved by CG with identity preconditioning and
`SolverControl(1000, 1e-12)`, then written as VTK. The standalone program runs
both 2D and 3D; optimization uses 2D only.

The original class has this structure (shortened declaration):

```cpp
template <int dim>
class Step4
{
public:
  Step4();
  void run();

private:
  void make_grid();
  void setup_system();
  void assemble_system();
  void solve();
  void output_results() const;

  Triangulation<dim> triangulation;
  FE_Q<dim> fe;
  DoFHandler<dim> dof_handler;
  SparsityPattern sparsity_pattern;
  SparseMatrix<double> system_matrix;
  Vector<double> solution;
  Vector<double> system_rhs;
};
```

The [stripped baseline](source/baseline/step-4-stripped.cc) preserves the
license and non-comment token stream. The [adapted source](source/adapted/step-4.cc)
opens these native operations:

| Adapted member | Purpose | Consumer |
| --- | --- | --- |
| `prepare_for_external_use()` | Perform grid/setup/assembly without immediately solving and writing | A and B |
| `system_matrix_view()`, `system_rhs_view()` | Inspect the existing assembled system | A and B |
| `solve(rhs, solution)` | Reuse CG with supplied RHS and in/out solution; return native monitored evidence | A and B |
| `output_results(state, path)` | Write a supplied state through the original VTK machinery | A and B |
| `dof_handler_view()`, `boundary_values_view()` | Supply native discretization and boundary data for mass assembly and coordinates | B |

The adapted `run()` still performs the original forward computation.
`STEP4_NO_MAIN` permits reuse of the tutorial translation unit. No nmopt
include or type enters Step-4. The adapted source has 244 code-bearing lines
against 191 in the stripped baseline: a shared net increase of 53 lines.
This is net source growth, not a count of every edited line or a claim that
all 53 lines would be needed by an already reusable application.

## 2. Application-owned OCP mathematics

### A: algebraic RHS control

[ProblemA](integration/problem_a.hpp) owns its prepared Step-4 instance.
State and control are full 289-entry coefficient vectors, with identity
control coupling and metric:

```math
\begin{aligned}
E(y,u)&=Ay-b-u,\\
J(y,u)&=\frac{1}{2}y^{\mathsf T}y+\frac{1}{2}u^{\mathsf T}u,\\
Ay&=b+u,\qquad A^{\mathsf T}p=y,\qquad r=j'(u)=u+p.
\end{aligned}
```

Control acts on the already boundary-treated algebraic equations, including
the boundary rows. A deliberately avoids FE coupling and mass assembly; it
is not distributed forcing with the original fixed physical boundary data.

### B: FE distributed control

[ProblemB](integration/problem_b.hpp) borrows a prepared Step-4 application.
The full continuous $`Q_{1}`$ control has 289 coefficients; the state and test
spaces use 225 free coordinates. Let $P$ inject free coordinates, $\ell$
contain fixed boundary values and zeros on free entries, and $M$ be the full
consistent mass matrix assembled on the original mesh:

```math
\begin{aligned}
y_{\mathrm{phys}}&=Pz+\ell,\qquad
K=P^{\mathsf T}AP,\quad b_{F}=P^{\mathsf T}b,\quad B=P^{\mathsf T}M,\\
E(z,u)&=Kz-b_{F}-Bu,\\
J(z,u)&=\frac{1}{2}(Pz+\ell)^{\mathsf T}M(Pz+\ell)
       +\frac{1}{2}u^{\mathsf T}Mu,\\
K^{\mathsf T}p&=P^{\mathsf T}M(Pz+\ell),\qquad
r=j'(u)=Mu+B^{\mathsf T}p,\quad g=M^{-1}r.
\end{aligned}
```

B is FE volume forcing, $-\Delta y=f+u$, with zero desired state and unit
regularization. The original $b$ already contains the lifting correction.
Boundary-node control coefficients contribute to volume forcing; they do
not change the prescribed state boundary values. The mass matrix is neither
lumped nor boundary-eliminated.

B supplies state RHS $b+PBu$ and adjoint RHS $`PJ_{z}`$ to the original full
solver, then restricts its output. Reusing CG for the transpose is justified
by the verified symmetry. Each OCP solve starts from zero and retains
Step-4's CG policy and exception behavior. The native mass inverse uses a
fresh CG solve with identity preconditioning, at most 1000 iterations, and
threshold $`\max(10^{-14},10^{-12}\lVert r\rVert_{2})`$. Both paths use that
same metric service. The [B protocol](../../../docs/planning/review/external-dealii-boundary-evaluation/problem-b-protocol.md)
contains the complete frozen discretization and policy.

## 3. Mathematical operations and the public connection

The minimal [A](minimal/problem_a_binding.hpp) and
[B](minimal/problem_b_binding.hpp) bindings expose the same contract shape.
The B mapping is:

| Operation | Native implementation | Public representation |
| --- | --- | --- |
| $E(z,u)$ | `ProblemB::residual` | Test-layout covector from residual callback |
| $E'(z,u)(v,w)=Kv-Bw$ | `ProblemB::residual_jvp` | Test-layout covector from JVP callback |
| $`E'(z,u)^{\ast}q=(K^{\mathsf T}q,-B^{\mathsf T}q)`$ | `ProblemB::residual_vjp` | Full state/control covector from VJP callback |
| $J(z,u)$ | `ProblemB::objective` | Scalar objective callback |
| $`(J_{z},J_{u})`$ | `ProblemB::objective_derivative` | Full state/control covector from derivative callback |
| $`Kz=b_{F}+Bu`$ | `ProblemB::solve_state` | `StateAdjointSolversT::solve_state` |
| $`K^{\mathsf T}p=J_{z}`$ | `ProblemB::solve_adjoint` | `StateAdjointSolversT::solve_adjoint` |
| $Mu$ and $`M^{-1}r`$ | `ProblemBMetric::apply`, `inverse_apply` | Local `MetricT` adapter |
| $`r=J_{u}-E_{u}^{\ast}p`$ | Composed from the preceding operations | `ReducedDTOT::augment_derivative` |

State/adjoint callbacks wrap native vectors in one-block layouts and return
actual `LinearSolveReport` evidence. The metric inverse checks its native CG
result but returns only a primal block; public `MetricT` has no solve-report
return. Application objects and matrices remain outside the generic solver.

A successful value/derivative evaluation follows these function calls:

```text
ReducedSearchSolverT::solve(initial_control)
  ReducedDTOT::evaluate_value(control)
    state-solve callback -> ProblemB::solve_state -> Step4::solve
    objective callback -> ProblemB::objective
  ReducedDTOT::augment_derivative(value)
    objective-derivative callback -> ProblemB::objective_derivative
    adjoint-solve callback -> ProblemB::solve_adjoint -> Step4::solve
    VJP callback -> ProblemB::residual_vjp
    form the reduced control covector
  direction policy -> metric adapter -> ProblemBMetric::inverse_apply
  line-search policy -> further value evaluations
  accepted trial -> derivative augmentation using the retained state
```

A changes the dimensions, mathematical operations, and metric implementation;
it does not change this construction pattern. B's
[main](minimal/problem_b.cc) reconstructs
`result.final_evaluation.state.block(0)` and calls Step-4's writer. A's
[main](minimal/problem_a.cc) writes its full retained state through `ProblemA`.
The [overview](external-integration-overview.md#7-what-changed-from-a-to-b)
compares the mathematical changes and unchanged framework roles.

## 4. Source accounting

Counts below are physical lines containing C++ code after removing blank and
comment-only lines, including block-comment licenses. Includes, declarations,
and brace-only lines count. Inline comments do not remove a code-bearing
line. The existing comment-strip utility supplies the transformation; a
reproduction command appears below. Source ranges in the next section refer
to the original files at `2ba749b`, before comment removal.

This is a snapshot inventory of the implemented paths, not development time,
a cumulative diff, or a universal lower bound. Files may combine mathematics,
checks, result structures, and optional instrumentation. Their role identifies
ownership, not a claim that every contained line is mathematically necessary.

### Native application and functional consumers

All paths in this table are relative to this directory. Shared source is
counted once; A and B columns contain separate files.

| Responsibility | A | B | Shared | Owning source |
| --- | ---: | ---: | ---: | --- |
| Forward application baseline | – | – | 191 | `source/baseline/step-4-stripped.cc` establishes the baseline |
| Net native reuse adaptation | – | – | 53 | `source/adapted/step-4.cc` is 244 total |
| OCP operations and native solve wrapping | 216 | 232 | – | `integration/problem_a.hpp`, `integration/problem_b.hpp` |
| Free/full state coordinates and lifting | – | 118 | – | `integration/problem_b_coordinates.hpp` |
| FE mass assembly and control coupling | – | 183 | – | `integration/problem_b_mass.hpp` |
| Native mass metric service | – | 81 | – | `integration/problem_b_metric.hpp` |
| Public nmopt binding, including metric adaptation | 206 | 222 | – | `minimal/problem_a_binding.hpp`, `minimal/problem_b_binding.hpp` |
| Consumer, solver policy, reporting, and output paths | 107 | 112 | – | `minimal/problem_a.cc`, `minimal/problem_b.cc` |
| **Functional source inventory** | **529** | **948** | **244** | **1,721 combined, counting shared Step-4 once** |

The role totals are:

- **Native OCP implementation:** 830 lines (216 A and 614 B). Including the
  shared 53-line net reuse increase gives 883 lines of application-side source
  growth over the forward baseline.
- **Minimal nmopt-specific bindings:** 428 lines (206 A and 222 B).
- **Mixed consumer entry points:** 219 lines. These contain both nmopt setup
  and ordinary program policy, error handling, reporting, and filesystem work.
  They are not counted again as binding or OCP lines.

A alone uses 773 lines including the adapted application; B alone uses 1,192.
Their sum would count Step-4 twice. The combined source increase over the
191-line forward baseline is 1,530, split into the three roles above. These
are source-footprint statements, not an estimate of how little code a new
application could require. Optional instrumentation remains in the reused
OCP headers, although the minimal consumers instantiate none.

### Comparison and verification support

These scopes are outside the minimal consumers' execution path. Alternative
evaluated bindings are listed separately because they are functional bindings
used for comparison, not pure diagnostic code.

| Responsibility | Lines | Exact source scope |
| --- | ---: | --- |
| Alternative evaluated bindings | 571 | `integration/nmopt_binding.hpp` (235), `integration/nmopt_problem_b_binding.hpp` (336) |
| Native reduced references and optimizers | 786 | `evaluation/native_reduced.hpp` (81), `native_optimization.hpp` (246), `native_problem_b_reduced.hpp` (145), `native_problem_b_optimization.hpp` (314) |
| Shared frozen comparison policy | 24 | `evaluation/optimization_policy.hpp` |
| Native verification and scenarios | 959 | `verification/verification.hpp` (294), `problem_b_verification.hpp` (589), `scenario.hpp` (76) |
| Optional diagnostic records | 207 | `diagnostics/instrumentation.hpp` |
| Native contract driver and failure-evidence support | 2,976 | `tests/dealii/external_step4_native_contract.cc` (2,713), `external_step4_evidence.hpp` (263), relative to repository root |
| A application contract drivers | 2,275 | `tests/application/external_step4_nmopt_contract.cc` (939), `external_step4_optimization_contract.cc` (885), `external_step4_minimal_problem_a_contract.cc` (451) |
| B application contract drivers | 2,423 | `tests/application/external_step4_problem_b_nmopt_contract.cc` (701), `external_step4_problem_b_optimization_contract.cc` (1,017), `external_step4_minimal_problem_b_contract.cc` (705) |
| **Support subtotal, excluding alternative bindings** | **9,650** | The preceding support scopes are disjoint |
| **Support plus alternative bindings** | **10,221** | Separate from the 1,721-line functional inventory |

This explains the large experiment footprint: it includes an independently
implemented native optimization schedule, equation and derivative audits,
dense oracles, comparisons, work attribution, and failure-evidence tests.
An application integration does not need to reproduce that harness. Native
validation appropriate to its mathematics is still necessary.

The inventory counts C++ sources only. Documentation, generated run data,
Python source/fidelity tools and their tests, and CMake target registration
are outside these totals. The two additional unchanged provenance copies
(upstream and stripped) are also excluded from the functional total.
All counts describe the present files, including earlier corrections; they
are not a measure of lines authored during only the minimal-consumer follow-up.

## 5. Inside the minimal bindings

The following disjoint regions reconcile exactly to each binding's total.
Ranges are inclusive physical line numbers in the linked
[A](minimal/problem_a_binding.hpp) and [B](minimal/problem_b_binding.hpp)
files at `2ba749b`; counts exclude blank/comment-only lines within each range.

| Region | A lines / range | B lines / range | Responsibility |
| --- | --- | --- | --- |
| Metric class | 51 / 24–81 | 58 / 29–93 | Identity or native mass forwarding, layout checks, ID, and local storage |
| Layout construction | 9 / 97–105 | 9 / 112–120 | Variable and test layouts |
| Contract member initialization | 7 / 106–112 | 7 / 121–127 | Model, partition, metric adapter, and reduced DTO construction |
| Model callback factory | 39 / 155–194 | 39 / 170–209 | Five wrappers and factory scaffolding |
| Solve-report translation | 15 / 139–153 | 15 / 154–168 | Actual native policy and convergence evidence |
| Solve-service factory | 23 / 196–218 | 24 / 211–234 | Two wrappers and factory scaffolding |
| Remaining packaging | 62 / complement | 70 / complement | Includes, namespace/class declarations, aliases, constructor setup, accessors, deleted copy/move operations, and members |
| **Total** | **206** | **222** | Every code-bearing line appears once |

The metric class count includes its own aliases, checks, and members; those
lines are not counted again as packaging. B's remaining constructor setup
also constructs its native metric. These regions expose the implementation
shape without classifying every brace or declaration as a separate API
obligation.

Layouts, correctly typed callbacks, native solve results, metric operations,
and reduced construction are required by the current path. The class names,
member organization, and accessors are local choices. The repeated wrappers
are visible in both bindings; the metric and mathematical operations they
wrap differ. A future convenience could package repeated construction, but
removing residual/JVP requirements or requesting only a control VJP would
change the underlying capability boundary and requires a separate decision.

### Reproduce the counts

Run from the repository root. This reads source only, reuses the existing
comment utilities, and prints the inventory and region reconciliation:

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
files = {p for p in base.rglob('*') if p.suffix in ('.cc', '.hpp')}
files.update(Path('tests/dealii').glob('external_step4*.cc'))
files.update(Path('tests/dealii').glob('external_step4*.hpp'))
files.update(Path('tests/application').glob('external_step4*.cc'))
for path in sorted(files):
    print(count(path.read_text()), path)

regions = {
    'a': [(24, 81), (97, 105), (106, 112), (155, 194), (139, 153), (196, 218)],
    'b': [(29, 93), (112, 120), (121, 127), (170, 209), (154, 168), (211, 234)],
}
for problem, spans in regions.items():
    text = (base / 'minimal' / f'problem_{problem}_binding.hpp').read_text()
    lines = text.splitlines(keepends=True)
    parts = [count(''.join(lines[first - 1:last])) for first, last in spans]
    packaging = count(text) - sum(parts)
    print(problem, 'regions', parts, 'packaging', packaging, 'total', count(text))
PY
```

## 6. What another application would need

The source roles identify where analogous work would belong, not its expected
size. An application with callable assembly/solve/output can reuse those
interfaces directly. If it already implements the chosen OCP and adjoint,
those operations can also be reused. Otherwise their mathematical definition
and validation are application work, regardless of the optimizer selected.

Connecting those operations to the present nmopt path requires compatible
state/control/test layouts, five model callbacks, two solve services, a
metric, and reduced/optimizer construction. The current DTO supports one
state block, one control block, and one residual-test block. Output uses the
application's physical reconstruction. The
[overview's adoption sequence](external-integration-overview.md#9-connecting-another-application)
and [API reference](../../../docs/reference/external-dealii-solver-integration.md)
show how these responsibilities fit together.

## 7. Evidence and reproduction

| Property | Evidence |
| --- | --- |
| Forward fidelity | Upstream, stripped, and adapted programs agree in original 2D/3D output comparisons. |
| Mathematical correctness | Off-solution residual/JVP/VJP checks, reduced derivatives, fresh equation audits, and independent dense optimum oracles pass. |
| Native/nmopt equivalence | The same problem and policies produce matched reduced evaluations and optimization traces. |
| FE metric verification | B's independent dense mass audit checks final stationarity separately from runtime CG inversion. |
| Consumer fidelity | Actual minimal executables reproduce audited native reports and VTK fields. |
| Regression evidence | Implementation handoff passed 195/195 deal.II and 67/67 neutral tests. |

The [A report](../../../docs/planning/review/external-dealii-boundary-evaluation/g1-report.md),
[B report](../../../docs/planning/review/external-dealii-boundary-evaluation/problem-b-report.md),
and [minimal-consumer assessment](../../../docs/planning/review/external-dealii-boundary-evaluation/minimal-consumers-report.md)
retain exact revisions, tolerances, counts, and artifact locations. These are
existing implementation results; this documentation reorganization does not
claim a new numerical experiment. A and B use different objectives and
metrics, so their different iteration counts are not a performance comparison.

For the runnable consumers, follow the [minimal README](minimal/README.md).
For source fidelity, the three standalone targets are
`nmopt_external_tutorial_step_4`, `nmopt_external_tutorial_step_4_stripped`,
and `nmopt_external_tutorial_step_4_adapted`. They use deal.II directly without
linking nmopt. From the repository root:

```bash
python3 tools/external_dealii/strip_comments.py \
  --input apps/external-dealii/step-4/source/upstream/step-4.cc \
  --check apps/external-dealii/step-4/source/baseline/step-4-stripped.cc
./build.sh configure debug-dealii
./build.sh build debug-dealii --target nmopt_external_tutorial_step_4
./build.sh build debug-dealii --target nmopt_external_tutorial_step_4_stripped
./build.sh build debug-dealii --target nmopt_external_tutorial_step_4_adapted
python3 tools/external_dealii/check_forward.py \
  --upstream-executable build/debug-dealii/bin/nmopt_external_tutorial_step_4 \
  --stripped-executable build/debug-dealii/bin/nmopt_external_tutorial_step_4_stripped \
  --output-root runs/external-dealii/step-4/forward-comparison \
  --file solution-2d.vtk \
  --file solution-3d.vtk
```

For adapted fidelity, replace the stripped executable with the adapted target
and add `--stripped-label adapted`. The comparator checks stdout, geometry,
connectivity, cell types, and numeric VTK arrays in separate run directories.
The upstream tag, retrieval details, hashes, and legal attribution remain in
the [evaluation roadmap](../../../docs/planning/external-dealii-boundary-evaluation.md).

Once all external test targets are built, select the complete Step-4 checks:

```bash
ctest --test-dir build/debug-dealii --output-on-failure \
  -R '^nmopt\.(external_tutorial_step_4|external\.tutorial_step_4)\.'
```

Generated evidence stays ignored under `runs/external-dealii/step-4/`,
relative to the command's working directory, with build-profile artifacts
under the corresponding build tree. The minimal consumers instantiate no
instrumentation, reference optimizer, oracle, or comparison runner. Reused
application headers retain optional diagnostic types as source dependencies.

The evaluated cases are linear, symmetric, serial, fixed-mesh, and
unconstrained in the control. Mandatory residual/JVP construction and unused
state-VJP work remain known limitations. Generality, newcomer effort,
performance materiality, and any future helper are separate questions; their
full disposition is in the [closure audit](../../../docs/planning/review/external-dealii-boundary-evaluation/closure-report.md).
