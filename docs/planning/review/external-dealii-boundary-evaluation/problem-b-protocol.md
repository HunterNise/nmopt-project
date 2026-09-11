# External deal.II Problem B evaluation protocol

Status: PB0 protocol frozen on 2026-09-11 following the request to make the
review corrections and documentation updates. No Problem B implementation or
numerical result is claimed. The first implementation unit is PB1 below.

This is the execution contract for the next bounded experiment, not a new
nmopt interface or a universal FE-control policy. The
[Problem A roadmap](../../external-dealii-boundary-evaluation.md) retains A's frozen
protocol and history. Its [G1 report](g1-report.md)
records the successful comparison; the
[design investigation](design-investigation.md#8-post-g1-review-and-problem-b-candidate)
preserves the reasoning that led to this selection.

## 1. Question, baseline, and scope

Test whether genuine FE distributed control, fixed essential state data,
rectangular coupling, and a nonidentity metric remain application-owned when
the same problem is evaluated and optimized natively and through current
public nmopt contracts. Attribute additional work to application access,
OCP mathematics, native orchestration, binding construction, and runtime
services separately. Problem B is another OCP on Step-4, not a second sample
of adapting an unrelated external application.

Retain the completed refactor at `67104fd376a9b1d25d7e6cdaca664a7fe1aa98ed`,
the exact upstream/stripped Step-4 sources, and corrected Problem A. The
starting committed evaluation revision is `0bb1307`; the subsequent native
trace-order correction and this protocol must be identified by their actual
commit or dirty source state in implementation handoffs. Stay on
`codex/evaluate/external-dealii-boundary`; do not merge changes backward.

Shared nmopt implementation, contracts, compiler/lowering, backend/storage,
generic helpers, and shared test support remain frozen. Existing public
contracts may be used by the nmopt binding. The native path has no nmopt
numerical or contract dependency. Do not add a builder, reduced port, generic
optimizer, or progress framework to prepare the comparison.

The selected case is linear, symmetric, serial, fixed-mesh, and unconstrained
in the control. Fixed state boundary data are part of the PDE discretization,
not an optimizer constraint. Nonlinearity, nonsymmetry, MPI, adaptivity,
additional optimization algorithms, timing/profiling, package/install work,
and a compiler-versus-external comparison are outside scope.

## 2. Frozen discretization and coordinates

Use original Step-4 in two dimensions: $\Omega=[-1,1]^{2}$, four global
refinements, continuous $`Q_{1}`$ state FE, original DoF numbering, forcing
$`f(x)=4(x_{1}^{4}+x_{2}^{4})`$, and prescribed boundary values
$`y(x)=x_{1}^{2}+x_{2}^{2}`$. Preserve the original forcing/stiffness assembly
and its two-point tensor Gauss rule; do not reassemble a more accurately
integrated forward PDE for B. Both standalone dimensions remain fidelity
checks; optimization remains 2D.

The control uses the same full continuous $`Q_{1}`$ basis and coefficient
ordering. It represents a volume field with no boundary restriction. There
are 289 control/physical-state coefficients and 225 free state coordinates;
derive index maps from the application's boundary data and verify these
counts, rather than hardcoding an index pattern.

Let $F$ be the free indices, ordered by increasing native global DoF index,
and let $P$ inject those coordinates into a full vector. Let $\ell$ contain
the original prescribed boundary coefficients and zero entries on $F$.
With $A,b$ denoting Step-4's already boundary-eliminated matrix and RHS:

```math
y_{\mathrm{phys}}=Pz+\ell,\qquad
K=P^{\mathsf T}AP,\qquad b_{F}=P^{\mathsf T}b.
```

The lifting correction is already present in $b$; do not subtract it again.
State/test coordinates have dimension 225, control coordinates dimension 289.
State tangents and adjoints embed using $P$ with zero boundary entries. Native
output and physical-state tracking use the reconstructed field.

Assemble one full consistent FE mass matrix with two Gauss points per
coordinate direction on each original cell:

```math
M_{ij}=\int_{\Omega}\phi_{i}\phi_{j}\mathrm{d}x,
\qquad B=P^{\mathsf T}M.
```

Use the same quadrature realization for tracking and regularization. No mass
lumping, boundary elimination of $M$, zeroing of control columns, or additional
control DoF handler is selected. $M$ is symmetric positive definite; $B$ is
rectangular. Store/apply the row restriction locally; explicit storage of $B$
versus applying $M$ followed by restriction is an implementation choice whose
work and copy consequences must be recorded.

## 3. Frozen optimal-control problem

Select zero desired state, $`y_{d}=0`$, and $\alpha=1$:

```math
\begin{aligned}
E(z,u)&=Kz-b_{F}-Bu,\\
J(z,u)&=\frac{1}{2}(Pz+\ell)^{\mathsf T}M(Pz+\ell)
       +\frac{\alpha}{2}u^{\mathsf T}Mu.
\end{aligned}
```

This is FE distributed forcing, $-\Delta y=f+u$. The objective includes the
inhomogeneous physical boundary lifting, including its cross terms. Controls
associated with boundary nodes remain volume coefficients; they do not change
the prescribed state values or the eliminated boundary equations.

Use Euclidean coefficient dual pairing for covectors. The complete executable
operations, including arbitrary off-solution test seeds $q$, are:

```math
\begin{aligned}
E'(z,u)(v,w)&=Kv-Bw,\\
E'(z,u)^{\ast}q&=(K^{\mathsf T}q,-B^{\mathsf T}q),\\
J_{z}&=P^{\mathsf T}M(Pz+\ell),\qquad J_{u}=\alpha Mu,\\
Kz&=b_{F}+Bu,\qquad K^{\mathsf T}p=J_{z},\\
r=j'(u)&=\alpha Mu+B^{\mathsf T}p,\qquad g=M^{-1}r.
\end{aligned}
```

The native reduced path needs the control pullback $`-B^{\mathsf T}p`$.
The public model still requires a valid full VJP. Verify residual/JVP/full
VJP separately and count the actual runtime consumers; do not substitute
zero or an on-solution-only expression for an unused block.

## 4. Ownership and numerical services

Step-4 retains its mesh, FE/DoF objects, boundary prescription, forward
assembly, supplied-RHS CG service, and native writer. Add only native read-only
access to the DoF/FE and boundary information needed for B assembly and
coordinates. Do not expose nmopt types or put OCP mass assembly into Step-4.
Retain its standalone `run()` behavior and both-dimensional output fidelity.

The B control problem owns its free/full maps, mass/coupling operations,
objective and derivatives, adjoint interpretation, and native metric service.
Its lifetime must cover borrowed Step-4 data and any sparse-matrix sparsity
patterns. The binding maps these operations, dimensions, and actual solve
reports into existing public contracts. It owns no alternate PDE assembly.

State solve: pass $`b+PBu`$ to the existing full-system solve, then restrict
the returned state to $z$. Adjoint solve: pass $`PJ_{z}`$ with zero boundary
RHS entries, solve the same full system, and restrict to $p$. Reusing CG for
the transpose is justified only by the verified symmetry of this fixed case.
Do not clamp a corrupted state after a solve to make boundary checks pass.

For both state and adjoint use zero full-vector initial guesses, a fresh
`SolverControl(1000, 1e-12)`, original CG and `PreconditionIdentity` on every
call. Retain actual monitored convergence evidence and native exception
propagation. No direct solver, factorization cache, cross-control solve cache,
or revised PDE tolerance is allowed. Original standalone initialization is
unchanged. Assemble once per application instance; retain an accepted trial's
state for derivative augmentation and output without another state solve.

The control metric is $M$, independently of objective regularization. Native
metric application is one $M$ action. For metric inversion use serial CG with
identity preconditioning, zero initial vector, maximum 1000 iterations, and
absolute stopping threshold:

```math
\tau_{M}(r)=\max(10^{-14},10^{-12}\lVert r\rVert_{2}).
```

Construct solver/control afresh per inverse application and propagate failure.
Both paths call the same application-native metric operations through direct
calls or a thin `MetricT` adapter. This keeps numerical policy and measured
metric work identical without giving the native reference an nmopt dependency.
The public `MassMetric` demonstrates existing support but is not a reason to
introduce a second numerical implementation into this comparison.

This shared-space case remains favorable: $`B^{\mathsf T}=MP`$, so in exact
arithmetic $`g=\alpha u+Pp`$. The comparison deliberately exercises the declared
metric service on both paths; neither path uses that cancellation to bypass
the mass solve. Attribute this solve to the selected numerical policy, not to
an unavoidable cost of every FE control problem. This limitation does not
justify adding another control space during the experiment.

## 5. Matched optimization policy

Use double-precision native deal.II vectors and the unchanged serial nmopt
backend, in the same recorded compiler/build environment. Run comparisons
sequentially. Do not introduce path-specific flags, fast-math, or threading.

Use steepest descent and Armijo, with zero initial control, at most 5000
accepted iterations and 30 trials per iteration. Reset the trial step to 1,
halve on rejection, use Armijo fraction $10^{-4}$, and minimum step 0.
Select `gradient_norm` stopping with tolerance $10^{-6}$; set relative-gradient,
objective-change, and step tolerances to zero and leave objective target unset.
These values retain A's algorithm policy while changing the declared geometry.
They are frozen choices, not predicted iteration counts or observed successes.

At every stopping check compute $g$ through the declared mass solve, explicitly
apply $M$ to $g$, and use $`\sqrt{g^{\mathsf T}Mg}`$. Direction is $d=-g$.
After an accepted step, measure its norm by explicitly applying $M$ to the
actual coefficient update. Match current nmopt's arithmetic order rather than
replacing the stopping expression by $`\sqrt{r^{\mathsf T}g}`$, which can
differ when the mass solve is inexact. Armijo uses the covector pairing:

```math
\delta u=u_{\mathrm{trial}}-u,\qquad
s=r^{\mathsf T}\delta u,\qquad
J_{\mathrm{trial}}\leq J_{\mathrm{current}}+10^{-4}s.
```

Form trial controls from the retained current control, and compute the actual
update by subtraction; do not replace $s$ by a step-scaled direction pairing.
Check gradient convergence before declaring the accepted-iteration limit.
Require finite gradients/slopes and strict descent. Native solve exceptions
abort; they are not rejected Armijo trials. Record and reject nonfinite trial
objectives as the current solver does. If no trial is accepted, retain a
line-search failure. Do not tune either path after observing a mismatch.

Rejected trials compute only state and objective. Accepted trials retain
those results and add one adjoint/derivative evaluation. The native reference
implements these policies directly, without nmopt calls or a generic optimizer
abstraction. Preserve the A reference as its own fixed experiment.

## 6. Independent verification and acceptance

All numerical acceptance checks reject nonfinite inputs/results. Verification
uses separate instrumentation after snapshotting runtime counts. Successful
native/nmopt agreement alone does not verify shared OCP mathematics.

### Assembly, coordinates, and derivative checks

- Separate native/nmopt instances must assemble identical $A,b,M$ and free
  index maps before comparison; investigate any difference. Verify derived
  dimensions, $`P^{\mathsf T}P=I`$, reconstruction, homogeneous tangent/adjoint
  boundary values, and state boundary error at most $10^{-11}$ in maximum norm.
  Check the raw full-vector CG result before restriction and its agreement
  with reconstruction under the paired vector tolerance. Reconstruction alone
  must not hide an incorrectly controlled boundary row.
- Verify symmetry of $A,K,M$ with Frobenius relative error at most $10^{-14}$,
  normalized by $`\max(1,\lVert\cdot\rVert_{\mathrm{F}})`$. Verify positive
  definiteness of $K,M$ using independent dense Cholesky in verification.
- Check $`\boldsymbol{1}^{\mathsf T}M\boldsymbol{1}=4`$ within $10^{-12}$.
  Independently evaluate mass pairings, control weak actions against embedded
  free test fields, and physical tracking by cell quadrature with three Gauss
  points per coordinate. Compare to assembled actions/objectives with the
  scalar tolerance below; this check must not reuse the assembled matrices.
- Use A's deterministic control samples $0,0.1c,0.1r,0.1a$ in the common
  289-entry ordering, including a repeated ramp after an intervening control.
  For free-state/test vectors, apply the same index formulas at dimension
  225. Use nonzero independent off-solution state/control tangents and seed.
- Verify residual JVP and objective directional derivative by centered
  differences at $h=10^{-6}$ with scaled error at most $10^{-8}$. Verify full,
  state-only, and control-only JVP/VJP pairings with scaled error at most
  $10^{-12}$. Scaled error divides by the maximum of 1 and the compared norms
  or scalar magnitudes, as appropriate.
- At all four controls, use ramp and alternating control directions normalized
  in the mass norm. For reduced centered differences, use
  $h=10^{-2},10^{-3},10^{-4},10^{-5},10^{-6}$ and recompute perturbed states.
  Require two adjacent steps to satisfy error at most
  $`10^{-7}\max(1,|r^{\mathsf T}v|)`$. Retain the complete table.
- At zero control, check first-order Taylor remainders for both directions at
  $h=0.1,0.05,0.025$: positive, above ten times repeat-objective variation,
  with consecutive halving ratios in $[3.5,4.5]$. No silent point filtering.

### Solve, metric, and paired checks

Require successful native CG reports under the selected stopping policies.
Independently recompute state, adjoint, and metric equation residuals, using
the actual RHS of each equation, and require:

```math
\frac{\lVert Lx-q\rVert_{2}}{\max(1,\lVert q\rVert_{2})}\leq10^{-10}.
```

Audit both restricted and expanded state/adjoint equations at sample controls
and final iterates. Retain monitored and recomputed residuals separately.
For metric verification, use the deterministic ramp/alternating vectors as
primal and covector inputs: check inverse/apply recovery under the vector
tolerance below, positive metric pairings, and symmetry of cross-pairings.
Reject invalid or failed mass solves rather than treating them as exact.

Paired coefficient vectors (including $z,u,p,r,g$) and scalars use:

```math
\begin{aligned}
\lVert a-b\rVert_{2}&\leq10^{-11}+10^{-10}\max(\lVert a\rVert_{2},\lVert b\rVert_{2}),\\
|a-b|&\leq10^{-12}+10^{-11}\max(|a|,|b|).
\end{aligned}
```

Compare sample objectives, state/adjoint solutions, reduced covectors, metric
gradients, and solve evidence. During optimization compare stopping reasons,
trial/accepted sequences, objectives, steps, slopes, gradient/step norms, and
work schedules. Bitwise identity is not required. Unexplained first divergence
blocks acceptance; retain its signed Armijo margin and diagnose before changing
the protocol. Native output must preserve the same reconstructed field and
writer policy; compare parsed numerical VTK data, ignoring generated timestamps.

### Independent optimum oracle and final acceptance

Use a verification-only dense factorization of the full KKT system, assembled
directly from $K,B,M,P,\ell$. Define $`Q=P^{\mathsf T}MP`$ and
$`c=P^{\mathsf T}M\ell`$:

```math
\begin{bmatrix}
Q&0&-K^{\mathsf T}\\
0&\alpha M&B^{\mathsf T}\\
-K&B&0
\end{bmatrix}
\begin{bmatrix}z^{\ast}\\u^{\ast}\\p^{\ast}\end{bmatrix}
=\begin{bmatrix}-c\\0\\-b_{F}\end{bmatrix}.
```

Use dense LU with pivoting; do not form an inverse or call native/nmopt reduced
evaluation to build or solve this oracle. Shared assembled operators are
checked independently by the assembly tests above. Audit all three KKT block
residuals normalized by the maximum of 1 and the norms of the constituent
terms; each must be at most $10^{-10}$. Evaluate the native reduced covector
at the oracle control and check its dual mass norm at most $10^{-8}$.

Both optimizers must stop by the selected gradient tolerance. Recompute state,
adjoint, and covector at each final control using a fresh native verification
instance; compare with returned values. Compute the independent final
stationarity norm by a dense mass solve, not the optimizer's stored metric
gradient, and require it at most $1.1\cdot10^{-6}$. Require:

```math
\lVert u-u^{\ast}\rVert_{M}\leq2\cdot10^{-6}.
```

This is an absolute mass-norm distance, not A's Euclidean coefficient distance
or a relative bound. With $\alpha=1$, strong convexity in this metric motivates
the bound; passing still requires the independent oracle and gradient audits.
Retain signed objective gaps. Any failed threshold is a failed gate requiring
diagnosis or an explicit, justified protocol amendment on both paths.

## 7. Attribution and evidence

Use unique B run directories under `runs/external-dealii/step-4/problem-b/`,
with `native-verification/`, `reduced-evaluation/`, and `optimization/` roles.
Reuse existing standalone fidelity artifact roots; preserve A artifacts.
Raw runs and reviewed working attribution remain ignored. Tracked closure
documentation records findings, source state, environment, and commands;
runtime tests do not gather Git metadata or rewrite reviewed attribution.

Separate setup, runtime, verification, and output. Record state/adjoint CG
calls, iterations, residuals and failures; objective/derivative calls;
residual/JVP/full-VJP versus control pullbacks; metric apply/inverse calls and
CG work; explicit stiffness, mass, and coupling actions. Attribute a mass
action to its actual purpose, so objective mass work is not counted as metric
work. Record coordinate reconstruction/restriction and known copy volume as
source-derived where not instrumented. Unobservable allocations or internal
CG actions remain unknown; callback counts are not timing measurements.

For $N$ accepted iterations and $T$ total trials, the prescribed successful
schedule predicts $1+T$ state/objective evaluations, $1+N$ adjoint/derivative
evaluations and control pullbacks, $1+N$ metric inversions, and $1+2N$ metric
applications (gradient norms plus accepted-step norms) on each path. Native
uses a control pullback; current nmopt uses a full VJP. Residual/JVP runtime
calls should be zero. Assert and explain actual counts rather than modifying
work to force these identities. A's bookkeeping fields need not be reused.

Create evidence destinations before work can fail, record original exceptions
while instrumentation is alive, and serialize each available completed result
before subsequent checks or the other path can fail. Preserve failed status
and retry isolation. Internal partial histories of a throwing solver remain
unavailable unless already exposed; do not expand the shared API to obtain
them or label incomplete evidence a successful comparison.

## 8. Implementation units and decision gate

Paths below are relative to the repository root. New B files belong alongside
A under existing Step-4 `integration/`, `evaluation/`, `verification/`, and
`diagnostics/` directories, with distinct B names. Keep A's mathematics,
binding, and native reference intact. Native and nmopt comparison tests use
separate executables; direct CMake target/scenario registration is permitted.
No application runner or new general test-support layer is selected.

| Unit | Outcome and boundary | Gate and prospective commit |
| --- | --- | --- |
| PB1 – Native FE problem | Minimal Step-4 access seams; B coordinates, mass/coupling, residual/objective, solves and metric; native contract tests. | Standalone and A regressions, assembly/coordinate/operation checks. `feat(dealii): define native Step-4 distributed control` |
| PB2 – Native verification | B native reduced evaluator, independent oracle and derivative/metric audits; no nmopt binding. | All native mathematical gates pass. `test(dealii): verify the native distributed-control reference` |
| PB3 – Native optimization | B-specific native metric steepest-descent/Armijo reference and evidence. | Convergence, oracle, schedule, and failure gates. `test(dealii): add native distributed-control optimization` |
| PB4 – Public binding | Explicit B binding and sampled native/nmopt comparison using current API. | Paired values, derivatives, metrics, solves, and ownership pass. `test(dealii): bind distributed control through the public API` |
| PB5 – Matched optimization | Paired B optimization, traces, output and attribution. | Matched schedule and independent final acceptance. `test(dealii): compare distributed-control optimization paths` |
| G2 – Attribution review | Factual report in the existing external-boundary review folder; update phase status. | Bounded interpretation and next decision. `docs(dealii): report distributed-control boundary findings` |

These are work boundaries, not a requirement to finish a large unit in one
commit. Split a unit locally if its implementation grows beyond one coherent
reviewable change; preserve the gate and native-before-binding order. Follow
[agent routing](../../../../.agents/README.md), including focused tests and required
Debug pipelines, existing machine limits, and no release builds for this
correctness experiment. Add exact runnable commands as executables exist.

Current status: PB0 complete; PB1 is the next implementation unit and has not
started. Implementation authorization is separate from this documentation
update. At each handoff update the ignored unit record and the tracked current
status, naming tested revisions/dirty inputs and actual evidence locations.

G2 compares incremental responsibilities and required operations, not an
adapter line-count budget. Repeated mechanical construction can motivate a
scoped helper experiment. A boundary change needs a reproducible valid
capability/ownership obstruction or a justified unnecessary obligation that
a helper cannot remove. A deeper diagnosis needs traced propagation across
otherwise independent responsibilities. Do not infer a universal external
integration cost from two OCPs on Step-4, or promote hypotheses to design
merely because both optimizations converge.
