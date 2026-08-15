# Implementation roadmap and agent handoff

## Purpose

This is the ranked continuation plan after the v0 executable contract and the
first serial deal.II lowerer. It is intentionally dependency-ordered: finish
the useful vertical slice before broadening the semantic language, and broaden
the semantic language before adding advanced PDE variants.

The fast baseline before and after every task is the explicit backend-neutral
Debug profile:

~~~bash
cmake --preset debug-neutral
cmake --build --preset debug-neutral
ctest --preset debug-neutral
~~~

Run `debug-dealii` with `--parallel 1` as well for compiler, lowerer, backend,
or numerical changes. The [build convention](../../conventions/build.md) owns
the complete profile matrix and cache-recovery guidance.

## Current handoff state

Stage B common stabilization is complete on
`codex/ch5-ch6-development`. Batches R0
(`RF-020`), R1 (`RF-006`, `RF-009`, and `RF-016`), R2a (`RF-001` and the
relevant `RF-006` characterization), R2b (`RF-002` through `RF-005` and the
relevant `RF-006` cases), R2c (the current factual defect in `RF-008` plus
`RF-012`), and R3 (`RF-016` through `RF-019`) are complete. P4.1 and P4.2 are
ignored for the current ordered implementation run because their scope is too
broad. The conditional C1 and C2 preparation implementations and remediation
gates are acceptance-complete.
Their remediation slices now include one checked resolved request, plan-owned
scalar residual/data and objective/service lowering, a versioned
resolved-decision manifest, and common realized-map/space records. Independent
scalar recombination checks and the neutral, deal.II-enabled, and sanitizer
verification profiles close the documented gates. Their stable-ID index, owned
compilation session,
shared solve reporting, projection witness, independent scalar oracle, and
build-cost remeasurement remain retained foundations. P5.1 is acceptance-
complete after its coefficient-placement and typed-boundary remediation gates
passed. P5.2 is acceptance-complete for its selected bounded registrations.
The typed trace and negative-metric policies, explicit $H^{1}_{0}$ target-data
assumption, control-boundary realization, and realized observation-space
dimensions are covered by the semantic and deal.II contracts. The reviewed
operator formulas, exact-transpose solve, energy/weighted observation
assembly, and negative-metric formulas remain the selected bounded behavior.
P5.3's bounded C5.8 and C5.10 implementations are likewise
acceptance-complete. The typed transposition contract, complete fixed-boundary
coverage, custom-ID manufactured checks, realized observation dimensions, and
point/flux dual-residual, directional, and reduced-Taylor checks pass their
documented gates. The outward-normal face evaluation, immutable point
coordinates, and explicit exclusions remain the selected bounded policies.
The C5.6-style Neumann composition is complete. All selected scalar Section
5.11/P5.4 Dirichlet-control slices are acceptance-complete for their registered
combinations. The typed shared policies, closed registration matching, and
realized transformed-observation dimensions pass the documented semantic and
deal.II contracts. The partial fixed-precedence lifting, selected fractional
and tangential operators, conforming-trace equivalence, and explicit
exclusions remain the bounded policies for these registered slices. A general
nonconforming transposition lowerer remains unselected.

The following pieces exist and are tested:

| Layer | Existing artifact | Meaning |
|---|---|---|
| Typed algebra | `include/nmopt/contract/layout.hpp` | `PrimalBlockT` and `CovectorBlockT` are distinct typed wrappers, even when a backend uses one vector storage type. Their block storage is read-only after construction; checked algebraic updates preserve the declared dimensions. |
| V1 semantic graph | `include/nmopt/semantic/v1/{types,validation,reference_specs}.hpp` | Deal.II-free selected graph with safe incomplete states, whole-graph closure checks, explicit two-sided pairings, structural/policy diagnostics, and ID-based reference deltas. |
| V1 compiler | `include/nmopt/compiler/v1/{compiled_problem,dealii_compiler,dealii_scalar_plan}.hpp` | Backend-generic compiled package, typed manifest container, and stored scalar handler plan. The bounded scalar path consumes the resolved request and typed residual/service records; specialized registrations retain the explicit strategies listed in the [v1 capability table](../implementation/v1/semantic-compiler.md#registered-capabilities). |
| Operator contract | `include/nmopt/contract/executable_model.hpp` | Residual, JVP, VJP, objective, and objective derivative. |
| DTO workflow | `include/nmopt/contract/reduced_dto.hpp` | One state block, one decision block (control or parameter), one test block, externally supplied state/adjoint solves. |
| Formulation solves and lifetime | `include/nmopt/contract/linear_solve.hpp`, `include/nmopt/dealii/serial_spd_solver.hpp`, and `include/nmopt/compiler/v1/dealii_types.hpp` | Typed state/adjoint solve reports, one shared serial SPD policy/service for symmetric targets, recorded direct and exact-transpose solves for the P5.1 nonsymmetric target, an owned static-mesh compilation session, and detached reduced services that retain executable/session lifetime. |
| Reference oracle | `include/nmopt/reference/linear_quadratic_model.hpp` | Dense linear-quadratic model used to test signs and derivatives independently of deal.II. |
| deal.II backend | `include/nmopt/dealii/serial_backend.hpp` | Serial Vector backend with a checked conversion from contract dimensions to the native deal.II size type. |
| Direct deal.II v0 lowerer | `include/nmopt/dealii/scalar_diffusion_reaction.hpp` | Preserved assembled scalar `FE_Q` diffusion-reaction reference with `FE_DGQ(0)` volume control, full-domain tracking, homogeneous Dirichlet data, and DTO solves. |
| deal.II metrics | `include/nmopt/dealii/{mass_metric,hminus1_metric,trace_hhalf_metric}.hpp` | One-block sparse SPD Riesz actions for the registered volume, boundary, trace, parameter, fractional-trace, and negative-norm layouts. The selected $H^{-1}$ realization applies $M_hK_h^{-1}M_h$; the $H^{1/2}$ realization applies the minimum-volume-$H^{1}$ Schur complement without forming a dense fractional matrix. All inverse actions use recorded serial-CG policies and operator-bound realization witnesses. |
| deal.II constraints | `include/nmopt/dealii/{cellwise,facewise}_box_constraint.hpp` | Coefficientwise boxes coupled to the actual positive-diagonal cellwise-volume or facewise-boundary $L^{2}$ metric realization, never to its display string. |
| Reduced solver | `include/nmopt/solvers/reduced_gradient.hpp` | Backend-parametric reduced search loop over `ReducedDTOT`, `MetricT`, and optional `ConstraintT`, with typed direction policies, explicit Hessian/Newton support, configurable Armijo/exact/Wolfe line searches, and uniform action reporting. |
| Build/test workflow | `CMakePresets.json` and `CMakeLists.txt` | Explicit neutral/deal.II Debug, neutral sanitizer, and deal.II Release profiles; requested dependency failures; target-scoped warnings; and labeled, time-bounded scenarios. |
| Tests | `tests/{reduced_dto_contract,semantic_v1_contract,dealii_diffusion_contract,dealii_trace_hhalf_metric_contract}.cc` | Four binaries expose forty-one independently named, labeled, and time-bounded CTest scenarios: five dense/backend contract cases, seven semantic graph/resolution/lowering-plan cases, and twenty-nine deal.II compiler/lowering/adapter cases. Negative checks identify exact diagnostics or contract failures, including block/layout, graph-closure, coefficient shape, boundary partition, binding, solve-policy, manifest-realization, lifetime, projection-coupling, observation topology/region, point-sensor mesh placement, normal-flux orientation and face-transpose behavior, and native-size invariants. Independent oracles additionally cover the exact trace Schur complement, tangential stiffness, loss/metric separation, all three remaining Section 5.11 stationarity compositions, and their reduced Taylor tests. |

The public v1 semantic path is deliberately not a general component compiler
yet. It resolves valid graphs by stable ID and has one bounded scalar
component-planning path; specialized graphs still select one of a bounded set
of target-specific implementations. The direct
deal.II v0 class remains a concrete reference lowerer rather than a public
problem hierarchy.

## Accepted Chapter 5/6 scope

The current project scope is a scalar, reduced-space DTO framework with
selected volume and boundary controls. The reusable Chapter 5 problem recipes
and the frozen Chapter 6 benchmark scenarios are owned by separate roadmaps;
this section records only the feature dependencies and implementation
priority.

The accepted scope is:

- keep the C1/C2 remediation gates closed; P5.2–P5.4 are acceptance-complete
  for their registered bounded slices;
- keep the selected scalar Section 5.11 and Neumann boundary-control slices;
- skip book Sections 5.12 and 5.13, corresponding to the roadmap's P5.5
  regularised state constraints and P5.6 Stokes/mixed-block work;
- retain the implemented selected P6.1 reduced-space policies, including
  nonlinear conjugate gradients, L-BFGS, line-search policies, and a
  Hessian-vector service with Newton/truncated-Newton support where the
  available derivative contract is sufficient;
- add a bounded P6.2 interface for executing an explicitly supplied OtD
  optimality system, without automatic continuous adjoint derivation;
- implement the scalar equality-constrained KKT product and PDAS services in
  P6.3 and P6.5;
- treat P6.4 preconditioning, other continuous-control bound policies, and
  the all-at-once numerical examples as conditional extensions.

The following remain outside the accepted scope: stabilization policies and
GLS-specific P6.2 work, automatic OtD derivation, measure-valued state
constraints, Stokes, and a generic continuous-control box semantics. A
preconditioner may be added when a selected all-at-once benchmark demonstrates
that it is needed; it is not a prerequisite for the reduced-space examples.

## Non-negotiable rules for every task

1. Preserve the global convention

   The global convention is $`\mathcal L(x,p)=J(x)-\langle p,E(x)\rangle`$.

   Thus the reduced covector is always $`J_{u}'-E_{u}'^{\ast}p`$.

2. Residuals and objective derivatives return covectors. A metric alone maps
   the reduced covector to a primal direction.
3. A VJP seed for a residual is primal in its test space. Do not replace it
   with a covector or add mass matrices unless the chosen discrete dual
   representation requires one.
4. New physics belongs in a residual-term lowerer, observations in an
   observation lowerer, and optimizer behavior in a solver. Do not introduce
   classes named after complete PDE/control combinations.
5. Every differentiable new coupling needs a value test, a JVP Taylor test,
   and a VJP pairing test. Every reduced feature needs a reduced-derivative
   Taylor test.
6. Unsupported combinations must yield a capability diagnostic. They must not
   quietly fall back to a nearby but different mathematical problem.

## Ranked work items

### P0.1 — Add a real deal.II $L^{2}$ metric — completed

**Why first:** The lowerer already assembles the control mass matrix
$`M_{u}`$, and the reduced DTO builder already produces $`j_{h}'(u)`$. The
missing operation is

$$
  \nabla_{L^{2}}j=M_{u}^{-1}j_{h}'.
$$

Without it, the result is a correct covector but not yet a usable search
direction for a deal.II optimizer.

**Implement:**

- A generic deal.II mass-metric implementation of `MetricT<SerialBackend>`.
- Apply with a sparse mass-matrix multiplication and inverse-apply with a
  declared linear solver/tolerance.
- Expose the control mass matrix through a narrow compiled-metric factory, not
  through optimizer knowledge of `ScalarDiffusionReactionModel`.

**Primary files:** add `include/nmopt/dealii/mass_metric.hpp`; minimally extend
the current lowerer with a metric factory or compiled-object registry.

**Done when:**

- applying the inverse followed by apply recovers a random covector;
- a metric direction satisfies
  $`\langle M_{u} g,\delta u\rangle=\langle j',\delta u\rangle`$;
- the test uses actual deal.II vectors and a nontrivial control mass matrix.

**Implemented:** `dealii_backend::MassMetric` provides the generic one-block
sparse-mass realization of `MetricT<SerialBackend>`. Its inverse action uses
serial CG with declared iteration and relative/absolute tolerances, and
`ScalarDiffusionReactionModel::control_l2_metric()` is the control-space
compiled-metric factory.

### P0.2 — Add a generic reduced Armijo gradient solver — completed

**Why second:** This closes the current vertical slice: the existing state
solve, adjoint solve, reduced covector, and new metric become a real
optimization loop. The solver must consume only `ReducedDTOT` and `MetricT`.

**Implement:**

- A backend-parametric reduced gradient solver with Armijo backtracking.
- Objective history, gradient norm in the selected metric, state/adjoint solve
  counts, and clear stopping reasons.
- An unconstrained path first; no PDE type checks or casts.

**Primary files:** add `include/nmopt/solvers/reduced_gradient.hpp` and a
deal.II integration test.

**Done when:**

- the objective decreases monotonically under the line-search criterion;
- the norm of the reduced covector/gradient reaches a configured tolerance on
  a manufactured linear-quadratic case;
- the same solver passes the dense reference model test unchanged.

**Implemented:** `solvers::ReducedGradientSolverT` evaluates every initial and
trial control with `ReducedDTOT::evaluate`, converts the DTO covector through
the supplied metric, and applies an unconstrained Armijo update. Its result
records objective and metric-gradient-norm histories, accepted iterations,
line-search trials, state/adjoint solve counts, and explicit stopping reason.

### P0.3 — Add the selected cellwise $L^{2}$ box constraint — completed

**Why now:** `FE_DGQ(0)` was chosen precisely because the $L^{2}$ box projection
is unambiguous: clipping one control coefficient per cell represents the
declared cellwise-constant admissible set.

**Implement:**

- A deal.II implementation of `ConstraintT<SerialBackend>` for coefficientwise
  lower/upper vectors on the control layout.
- Projected Armijo gradient updates, with a projected-gradient stopping
  measure.
- Data rules for bounds: initially constants or `FE_DGQ(0)` vectors on the
  control mesh only.

**Do not do:** reuse this projection for continuous controls, nodal bounds,
or an $H^{1}$ metric.

**Done when:** a test reaches a bound on a manufactured problem, preserves
feasibility every iteration, and satisfies a discrete projected-stationarity
criterion.

**Implemented:** `dealii_backend::CellwiseBoxConstraint` is the serial
coefficientwise `ConstraintT<SerialBackend>` realization. The lowerer exposes
only `FE_DGQ(0)` control factories for scalar constants or exact-layout
coefficient vectors. The constraint-qualified `ReducedGradientSolverT`
projects every trial in `l2_cellwise` and stops on the metric norm of the
projected-gradient residual.

### P1.1 — Build the narrow v1 semantic-to-compiler path — completed

**Why before additional PDE features:** The current lowerer is intentionally
concrete. Adding Neumann, boundary observation, or coefficients directly to
it would recreate the forbidden problem-class pattern.

**Implement only the current slice of a public semantic graph:**

~~~text
Region:       volume and boundary ids
Space:        scalar state/test/control/observation declarations
Variable:     state and control
Data:         forcing, target, constants, bounds
ResidualTerm: diffusion-reaction, volume source, volume control
Observation:  volume restriction (initially full domain)
Loss:         quadratic tracking and quadratic control regularisation
Metric:       L2
Constraint:   optional cellwise box
~~~

Add a semantic validator with separate structural, policy, lowerability, and
formulation-capability diagnostics. Add a lowerer registry that recognizes
only the listed v0 term kinds and creates the existing deal.II executable
objects.

**Primary files:** `include/nmopt/semantic/v1/{types,validation,reference_specs}.hpp`,
`include/nmopt/compiler/v1/{compiled_problem,dealii_compiler}.hpp`,
`docs/implementation/v1/semantic-compiler.md`, and focused semantic/deal.II
contract tests.

**Done when:** constructing the current problem through a `ProblemSpec` produces
the same residual, objective, and reduced derivative as the direct deal.II
reference setup.

**Implemented (v1):** `semantic::v1::ProblemSpec` and
`SemanticValidator` describe and validate the selected graph without backend
objects. `compiler::v1::DealiiCompiler` appends lowerability and formulation
diagnostics through a small explicit lowerer registry, then produces a
separately owned compiled executable. The v0 direct
`ScalarDiffusionReactionModel` remains unchanged as the reference path; v1
constructs its own instance from the semantic declaration, so both paths can
be compared without overwriting v0.

### P1.2 — Generalize fixed essential conditions through reconstruction — completed

**Why:** The current lowerer correctly handles only homogeneous Dirichlet
data. Fixed inhomogeneous data is the smallest meaningful test of the
transformation/reconstruction contract; it must happen before controlled
Dirichlet data.

**Choices and default:**

- **Default:** represent
  $`y_{\mathrm{phys}}=P_{h}\widehat y_{h}+\ell_{0,h}`$, where $`P_{h}`$ handles all
  affine constraints and $`\ell_{0,h}`$ is a declared fixed lifting.
- Alternative: work only with full vectors plus modified rows. This is the
  current homogeneous implementation and must not be generalized implicitly
  to nonzero data.

**Done when:** residual, observation, JVP, and VJP operate on the physical
field; a nonzero manufactured Dirichlet solution passes the reduced Taylor
test; and changing the lifting/data invalidates all appropriate caches.

**Implemented (v1):** A state variable can declare the
`fixed_dirichlet_reconstruction` transformation, and the compiler then binds
explicit fixed-Dirichlet `Function` data. Its private v1-only compiler target
compiles independent `FE_Q` coordinates with
$`y_{\mathrm{phys}}=P_{h}\widehat y_{h}+\ell_{0,h}`$, evaluates residual and
tracking on the physical field, and applies $`P_{h}^{\ast}`$ for every
state-side covector. The direct v0 homogeneous model is untouched. The
manifest records the transformation, nodal boundary interpolation data rule,
and lifting realization; recompilation is the immutable data-cache boundary.
The deal.II contract test covers a nonzero manufactured state, reconstruction
JVP/VJP and objective derivatives, a reduced Taylor remainder, missing data
diagnostics, and changed lifting data.

### P1.3 — Add subdomain observations and data projection policy — completed

**Why:** It exercises the observation/loss port without changing PDE physics.
It also forces the compiler to make region and data semantics real.

**Implement:**

- Named material/subdomain regions and an observation mask/restriction.
- A tracking mass/load assembled only on the observation region.
- Explicit data rule: analytic `Function` evaluated at quadrature or a declared
  FE projection/interpolation. No implicit nodal interpretation.

**Done when:** changing only the observation region changes the adjoint RHS
and objective, not $`A_{h}`$, $`B_{h}`$, the state solution, or solver code.

**Implemented (v1):** `RegionSpec` now names material-id volume subregions;
the state tracking observation, observation space, desired state, and selected
analytic-quadrature data policy are checked to use the same region. The
private v1 assembled target builds tracking mass/load/norm only on those
cells, leaving the state matrix, forcing, control coupling, state solution,
and generic solver interfaces untouched. `CompilationManifest` records both
the analytic `QGauss` target-data rule and material observation realization.
Focused semantic and deal.II tests reject an absent or region-mismatched
target policy and prove that changing only the material mask changes the
objective and adjoint RHS but not the state solution or residual action. FE
target projection/interpolation remains an explicit future lowerer; no nodal
interpretation is inferred.

### P2.1 — Add Neumann control and boundary tracking — completed

**Why:** This is the first nontrivial boundary composition test while keeping
the state parameterization fixed. It tests face integration, regions, trace
semantics, control layout, and pullbacks.

**Default first discrete choice:** facewise-constant boundary control on
marked boundary faces of the state triangulation. Assemble with `FEFaceValues`.

**Implement:**

- Boundary Region declarations and a boundary control layout.
- A residual term
  $`-\langle u,\gamma v\rangle_{\Gamma_{c}}`$ and its VJP.
- Boundary trace observation/loss on a marked boundary region.
- Independent boundary $L^{2}$ metric and box constraint realization.

**Done when:** the boundary coupling obeys the pairing test and a
finite-difference reduced derivative test. Do not call it a Dirichlet
control or use a generic boundary-load abstraction for both cases.

**Implemented (v1):** `make_neumann_boundary_control_problem()` declares
separate Dirichlet, control, and observation boundary regions; an explicit
`neumann_control` residual term; a state `boundary_trace` observation; and
the trace policies required for both. The private v1 Neumann target assigns
one scalar to each marked active boundary face and uses `FEFaceValues` to
assemble $`-\langle u,\gamma v\rangle_{\Gamma_{c}}`$, its transpose, and
boundary tracking. It has a separate diagonal face-measure `l2_facewise`
metric and `FacewiseBoxConstraint`, with separate facewise bound bindings.
Focused semantic tests reject an absent trace policy or a mismatched control
region; the deal.II contract test verifies the coupling pairing, trace-loss
derivative, facewise projection, manifest, and a reduced Taylor remainder.
The v0 volume-control model remains unchanged.

### P2.2 — Support pure Neumann with the selected mean-constraint policy — completed

**Why:** This is a correctness feature, not a convenience flag. It affects
state uniqueness, adjoints, observation meaning, metrics, and preconditioners.

**Default:** augment state and adjoint systems with one mean-zero
constraint/multiplier. Check compatibility of every volume and boundary load,
including controls.

**Done when:** incompatible data is rejected; compatible problems report the
gauge in their compilation manifest; state and adjoint are invariant under no
unrecorded pinning convention.

**Implemented (v1):** `make_pure_neumann_boundary_control_problem()` keeps
the facewise natural-boundary residual but selects a full-domain
`mean_zero_multiplier` policy. With zero reaction, its private v1 target
solves state and adjoint systems through an explicit one-multiplier saddle
matrix. Compilation rejects nonzero reaction and incompatible forcing; state
solves reject incompatible boundary controls. The manifest records the gauge
and `SparseDirectUMFPACK` solve, and focused contracts check zero means and no
hidden DoF pin. The v0 model remains unchanged.

### P2.3 — Add $H^{1}$ regularisation and $H^{1}$ search geometry separately — completed

**Why:** The current code makes the distinction possible but does not yet
exercise it.

**Implement in two different tasks or commits:**

1. $H^{1}$ **regularisation**: a new control loss and likely an `FE_Q` control
   space; it changes $`J_{u}'`$.
2. $H^{1}$ **metric**: a Riesz map/inverse used only for
   $G^{-1}j'$; it does not change the objective or adjoint.

**Implemented (v1):**
`make_h1_regularised_scalar_diffusion_reaction_problem()` selects a new
`quadratic_h1_control_regularisation` loss and continuous control `FE_Q` space.
Its private v1 target adds $`\alpha(M_{u}+K_{u})u`$ to the control objective
derivative while retaining the `l2_continuous` mass-matrix search metric.
`make_h1_metric_scalar_diffusion_reaction_problem()` selects the separate
`h1_continuous` metric on the same continuous control space. Its Riesz map is
$G=M_{u}+K_{u}$, with the positive mass term providing coercivity; its CG
inverse is used only to form the search direction. It does not alter the
objective, residual, state solve, or adjoint solve, and neither realization
selects a box constraint. Focused contracts compare the two metrics while
checking the stiffness contribution, residual pullback identity, and reduced
Taylor remainder.

The $H^{1}$ metric needs a positive zero-order term or an explicit
boundary/mean policy to be invertible. P5.2 separately registers one named
$H^{-1}$ realization with an explicit homogeneous-Dirichlet search space and
operator; no generic negative-norm default is inferred.

### P3.1 — Add coefficient identification and nonlinear first-order actions — completed

**Why:** This validates the parameter block and nonlinear residual paths while
remaining compatible with reduced DTO and L-BFGS.

**Implement:**

- Parameter variable, bounds/positivity policy, and optionally
  $m=\exp(q)$ transformation.
- Residual
  $(m\nabla y,\nabla v)$, JVPs in $y$ and $m$, and VJPs for both.
- State solve policy that reassembles at each parameter iterate.
- L-BFGS or Gauss–Newton only after reduced derivatives pass Taylor tests.

**Do not do:** introduce an inverse-problem solver type or a general nonlinear
KKT Newton method. The latter needs the separate second-order contract.

**Implemented (v1):**
`make_coefficient_identification_problem()` uses the binary reduced DTO
decision port for a `diffusion_parameter` variable rather than a source
control. The parameter is cellwise `FE_DGQ(0)` and must select a cellwise L2
box whose lower bound is strictly positive; this is the first explicit
positivity policy, with no logarithmic transformation. Its private v1 target
realises $`(m\nabla y,\nabla v)+(c y,v)-(f,v)`$, reassembles the SPD state
matrix for each parameter point, and supplies both JVP and VJP parameter
actions. The parameter L2 metric and box use the separate
`l2_cellwise_parameter` identifier. Focused contracts check positivity
diagnostics, reassembly, residual finite differences, residual pullback,
objective derivative, and a state-recomputed reduced Taylor remainder. The
generic first-order reduced workflow is reused; L-BFGS, Gauss–Newton, and
second-order KKT actions remain separate future solver work.

### P3.2 — Add Dirichlet control through an explicit lifting — completed

**Why later:** It is a transformation problem, not a load term. The fixed
lifting work in P1.2 is its prerequisite.

**Default:** compile an explicit map

$$
  y_{\mathrm{phys}}=P_{h}\widehat y_{h}+\ell_{0,h}+L_{D,h}u_{h}
$$

with forward value, JVP, and VJP. Both residual and observation consume
$`y_{\mathrm{phys}}`$.

**Done when:** the implementation passes a composed-transformation VJP test
and a reduced Taylor test. Reject boundary spaces or corner/interface
configurations with no declared lifting policy.

**Implementation:** `make_dirichlet_control_scalar_diffusion_reaction_problem()`
selects a separate v1-only assembled target with
$`y_{\mathrm{phys}}=P_{h}\widehat y_{h}+L_{D,h}u_{h}`$. The initial registered
trace policy has one shared nodal control coefficient for every `FE_Q` state
DoF on the complete exterior boundary; its $`\ell_{0,h}`$ term is zero. The
compiler rejects partial boundaries, undeclared corner or interface rules,
hanging-node relations, and box constraints rather than choosing a trace
value implicitly. Residual and full-volume tracking assemble on the physical
state, while state and control covectors use $`P_{h}^{\ast}`$ and
$`L_{D,h}^{\ast}`$. The deal.II contract test verifies a nonzero manufactured
state, composed lifting VJP, reduced Taylor remainder, trace metric, manifest,
and incomplete-boundary diagnostic. The v0 lowerer remains unchanged.

### P4.1 — Add the fixed-step temporal compiler — ignored

**Status:** ignored for the current ordered implementation run. Reactivate
only through an explicit user decision.

**Default:** a global residual for a fixed-step backward-Euler heat equation,
with the complete trajectory available for its exact transpose.

**Prerequisites:** stable reconstruction semantics, compilation manifests,
and state/adjoint solve diagnostics.

**Do not do:** hide a forward time integrator behind a separately coded
backward adjoint, or introduce adaptive steps/events before replay and
differentiated-control policies exist.

### P4.2 — Generalize algebra, execution, and formulations — ignored

**Status:** ignored for the current ordered implementation run. Reactivate
individual capabilities only through an explicit user decision.

This long-term group has value, but should not block the preceding vertical
slices:

1. Multiple state/equation blocks, mixed fields, and Petrov–Galerkin spaces.
2. Serial matrix-free equivalence, then distributed Trilinos/PETSc vector
   backends with ownership/ghost policies.
3. Automatic OTD formulation derivation; the selected scope provides only an
   executor for supplied optimality-system blocks.
4. Generic nonlinear KKT Newton/SQP and other second-order constrained
   methods beyond the selected scalar reduced Newton consumer.
5. Fractional metrics, point/flux observations, nonmatching meshes, and shape
   optimization only with their own explicit policies.

## Chapter 5 feature requests

The Chapter 5 implementation guide records the mathematical variants and the
required composition boundary for each request. The requests below are
ordered by reusable capability, not by the source chapter's presentation
order. They must not be implemented as a hierarchy of named textbook problem
classes.

### P5.1 — Compose general scalar elliptic volume and Robin boundary terms — acceptance complete

**Motivation:** C5.1, C5.5.1, and C5.6 use scalar operators beyond the
current diffusion-reaction slice. Their differences are residual terms and
boundary policies, not different optimisation frameworks.

**Declare and implement:**

- Tensor diffusion, conservative transport `div(b y)`, advective transport
  `c · grad(y)`, reaction, volume source, and volume-control residual terms.
- Robin bilinear boundary term and Robin boundary source, with a declared
  conormal-flux convention.
- A boundary-partition policy that distinguishes fixed Dirichlet, Neumann,
  Robin, and transport inflow/outflow regions.
- Requirement policies for uniform ellipticity, coefficient regularity,
  coercivity, and trace interpretation. They record assumptions; they are not
  theorem provers.

**First registered target:** serial scalar `FE_Q` state/test with selected
coefficient `Function` bindings, volume `FE_DGQ(0)` control, and one declared
Robin boundary region. It may initially retain homogeneous fixed Dirichlet
data; partial controlled Dirichlet boundaries belong to P5.4.

**Derivative and test contract:** each variable-dependent residual term must
have value, JVP, and VJP pairing tests. The combined target must pass a
non-symmetric transport adjoint test, a Robin boundary contribution test, and
a reduced Taylor test. Unsupported coefficient shapes, boundary overlaps, or
missing policies must return compiler diagnostics.

**Implemented (v1):**
`make_general_scalar_elliptic_robin_problem()` recombines independently
handled tensor diffusion, conservative transport, advective transport,
reaction, volume source/control, and Robin bilinear/source contributions in
`ScalarLoweringPlan`. Rank-specific deal.II `TensorFunction` bindings and
one-component scalar `Function` bindings carry separate provenance. The
serial `ScalarComponentModel` target retains `FE_Q` state/test and
`FE_DGQ(0)` control layouts, requires a disjoint complete homogeneous
Dirichlet/Robin boundary partition, assembles the declared conormal weak form,
and uses `SparseDirectUMFPACK` for the nonsymmetric state operator and its
exact-transpose adjoint solve. Focused contracts isolate every
state-dependent term's value/JVP/VJP actions, the Robin source, combined
transport transpose, reduced Taylor remainder, manifest, coefficient-shape
diagnostic, and boundary overlap/completeness diagnostics.

**Review status:** the
[P5.1 remediation review](review/chapter-5/p5.1-remediation-review.md) findings
are repaired: coefficient and Robin data have truthful semantic placement,
and the selected boundary/conormal/transport/trace realization is typed and
consumed by lowering. The required neutral and deal.II Debug verification
profiles pass, so P5.1 is acceptance-complete. The assembled weak signs,
nonsymmetric exact transpose, complete mesh partition, and existing
term/Taylor tests remain retained behavior.

### P5.2 — Add energy-volume and weighted-trace observations, then the selected $H^{-1}$ metric separately — acceptance complete

**Motivation:** C5.5.2 needs state tracking in an energy space, while C5.7
needs a boundary trace multiplied by declared data. The alternate C5.5.2
control-space formulation requires a separate $H^{-1}$ primal-dual map.

**First observation slice:**

- Add a full-domain `h1_state_restriction` observation with its explicitly
  declared $H^{1}_{0}$ pairing. Its quadratic loss lowers to selected mass and
  stiffness contributions, and its VJP supplies the corresponding state
  covector.
- Add `weighted_boundary_trace`, a general map from a state trace and fixed
  boundary weight data to a declared boundary observation space. It must not
  be encoded by modifying a tracking loss or Neumann residual.
- Record the target-data and trace quadrature policies in the manifest.

**Separate metric slice:** add an $H^{-1}$ metric only after stating its
control space, Riesz operator, boundary/mean policy, and inverse-solve
tolerances. It changes only the search direction; it does not change the
energy observation, residual, or adjoint equation.

**Done when:** state energy and weighted-trace value/JVP/VJP tests pass, the
metric apply/inverse tests establish the declared pairing, and a reduced
Taylor test distinguishes the $L^{2}$ and $H^{-1}$ search directions without
changing the reduced covector.

**Implemented slices:**
`make_h1_state_tracking_scalar_diffusion_reaction_problem()` selects the new
`h1_state_restriction` kind and an explicit $H^{1}_{0}$ observation pairing.
The bounded scalar component target assembles the mass-plus-stiffness tracking
operator, target value/gradient load and norm, and corresponding state
covector while retaining the cellwise $L^{2}$ control metric. The manifest
records value-and-gradient target quadrature and handler provenance. Focused
contracts cover semantic topology, unsupported subdomain lowering, an
independent polynomial energy value, the observation JVP/VJP derivative,
reduced Taylor convergence, and unchanged metric provenance. The
`make_weighted_boundary_trace_neumann_control_problem()` separately selects
`weighted_boundary_trace`, whose explicit immutable `boundary_weight` data
port is evaluated with the desired target at boundary face quadrature. Its
Neumann target assembles $h^{2}$ in the tracking operator and $h z_d$ in the
target load, while leaving the residual and `l2_facewise` control metric
unchanged. Semantic contracts cover the data port and quadrature policy;
deal.II contracts cover missing/provenance/shape diagnostics, exact
constant-weight value and pullback scaling, residual JVP/VJP identity,
unchanged metric action, reduced Taylor convergence, and manifest provenance.

`make_l2_metric_h1_state_tracking_continuous_control_problem()` and
`make_hminus1_metric_h1_state_tracking_scalar_diffusion_reaction_problem()`
select the same independent homogeneous-Dirichlet continuous `FE_Q` control
coordinates. The first retains $G_h=M_h$; the second selects
$G_h=M_hK_h^{-1}M_h$, where $K_h$ is the Dirichlet Laplacian, so
$G_h^{-1}=M_h^{-1}K_hM_h^{-1}$ and no mean constraint is needed. Both retain
the same energy observation, volume residual, $L^{2}$ control loss, state
solve, and adjoint solve. The registered target has no coefficientwise box.
The shared identity-preconditioned CG tolerances are recorded for the mass and
Laplacian inverse actions. Independent metric tests establish exact action,
inverse recovery, and symmetry; the compiled comparison establishes identical
reduced objectives/covectors, distinct $L^{2}$/$H^{-1}$ directions, quadratic
reduced Taylor remainders for both directions, and complete manifest
provenance. These implementation slices have landed.

**Acceptance status:** complete for the registered bounded slices. The typed
trace and negative-metric selections, control-boundary realization, explicit
$H^{1}_{0}$ target-data assumption, and realized weighted-observation
dimensions are covered by semantic and deal.II contracts. The mass-plus-
stiffness observation, weighted pullback, $M_{h}K_{h}^{-1}M_{h}$ action/inverse,
and comparison/Taylor tests remain retained behavior.

### P5.3 — Add normal-flux and point-sensor observations through an explicit strong/very-weak policy — acceptance complete

**Motivation:** C5.8 and C5.10 are not ordinary restriction/trace
observations. Their adjoints use boundary data or Dirac sources with lower
regularity. The implemented C5.11.2 boundary-control problem uses a related
transposition principle; this roadmap item registers bounded C5.8 and C5.10
slices while leaving alternate flux/evaluation formulations unselected.

**Declare and implement:**

- Point-set regions with immutable sensor coordinates, a selected discrete
  evaluation rule, and a finite-dimensional observation pairing.
- A normal-flux observation whose input has an explicit strong $H^{2}$ or
  declared $H(\mathrm{div})$ trace capability. The selected policy must state
  the normal orientation and observation boundary region.
- A transposition formulation policy declaring the strong test space $Y$, the
  isomorphism $T:Y\rightarrow H$, the residual codomain, and the multiplier
  space. It must identify which domain-regularity assumptions are supplied by
  the model author.

**First registered targets:** C5.8 strong-state normal-flux observation and
C5.10 point-sensor observation. C5.8 selects
$Y=H^{2}(\Omega)\cap H^{1}_{0}(\Omega)$, the outward normal, selected
boundary-face quadrature, and an assembled very-weak boundary-source
transpose. C5.10 selects physical point evaluation and its assembled
very-weak point-load transpose. Alternate flux/evaluation policies and a
general transposition lowerer remain rejected. Neither adjoint is treated as
an ordinary boundary trace or undeclared nodal approximation.

**Done when:** the selected observation derivatives pass value/JVP/VJP tests,
the low-regularity adjoint passes its dual residual test, and the compiled
product records all regularity, orientation, and evaluation policies.

**Implemented first slices:**
`make_normal_flux_scalar_diffusion_reaction_problem()` registers the C5.8
strong-state normal-flux graph. Its boundary region owns the selected IDs and
its observation space is an $L^{2}$ face-quadrature output. The selected
strong-state parent is
$Y=H^{2}(\Omega)\cap H^{1}_{0}(\Omega)$ with
$T=-\kappa\Delta+rI:Y\rightarrow L^{2}(\Omega)$ using the bound diffusion
and reaction data ports, under a declared convex-or-$C^{2}$ domain assumption.
The deal.II lowerer evaluates the outward normal
derivative with `FEFaceValues` at selected boundary-face quadrature and
assembles the weighted transpose of the same face map; the adjoint solve
therefore uses the declared very-weak boundary source. An $H(\mathrm{div})$
realization and general transposition lowerer remain rejected. The sibling
`make_point_sensor_scalar_diffusion_reaction_problem()` registers the C5.10
finite point-sensor graph with finite, unique physical coordinates, its
physical `FE_Q` evaluation, and its assembled very-weak point-load transpose.

**Acceptance status:** complete for the registered C5.8 and C5.10 slices. The
typed transposition declaration and validation, complete fixed-boundary
coverage, custom-ID manufactured checks, realized normal-flux dimensions, and
point/flux dual-residual, directional, and reduced-Taylor checks pass the
documented verification gate. Outward-normal face evaluation, immutable point
coordinates, and the explicit rejected alternatives remain unchanged. The
post-audit non-unit diffusion regression also solves the manufactured state
with consistently scaled forcing and checks both the state and discrete
residual, so operator provenance is exercised by the realized action.

### P5.4 — Generalize Dirichlet-control transformations and trace metrics — acceptance complete

**Acceptance status:** complete for the registered Section 5.11 slices. The
typed partition/interface/trace/transposition policies, closed registration
matching for supported combinations, and realized transformed-observation
dimensions pass the documented semantic and deal.II verification gate. The
selected fixed-precedence, fractional, tangential, conforming-trace, and
explicit exclusion policies remain bounded; unregistered cross-products stay
rejected.

**Motivation:** C5.11 and C5.13.2 require choices beyond the current complete
boundary nodal $L^{2}$ lifting: mixed fixed/controlled boundaries, nonzero
fixed data, fractional trace geometry, and tangential $H^{1}$ regularisation.

**Declare and implement:**

- A boundary-partition/lifting policy that explicitly resolves every
  fixed-controlled interface, corner, and hanging-node relation in
  $`y_{\mathrm{phys}}=P_{h}\widehat y_{h}+\ell_{0,h}+L_{D,h}u_{h}`$.
- Selected trace-space declarations and metrics for boundary $L^{2}$,
  fractional $H^{1/2}$, and tangential $H^{1}$ choices. A metric is separate
  from the associated control regularisation loss.
- A surface-gradient observation/map for the tangential $H^{1}$ loss, with
  boundary mesh, endpoint, and orientation policies.
- A transposition formulation for Dirichlet data below the ordinary
  $H^{1/2}(\Gamma)$ trace class, with explicit continuous parent spaces,
  domain regularity, discrete subspace, equivalence, and conormal policies.

**First registered target:** one partial scalar controlled boundary disjoint
from one fixed nonzero Dirichlet boundary, with a complete declaration of the
corner policy. Start with an $L^{2}$ trace metric; add fractional and
tangential metrics as separate follow-up targets rather than silently using a
mass matrix for all three.

**Done when:** the physical reconstruction and both state/control pullbacks
pass VJP tests under changed fixed and controlled data; each metric passes its
own pairing test; and incomplete or ambiguous boundary declarations are
rejected.

**Implemented first slice:**
`make_partial_dirichlet_control_scalar_diffusion_reaction_problem()` declares
one fixed-data port, a complete fixed/controlled boundary partition, and
fixed-data precedence at shared corner/interface DoFs. The private nodal target
uses the relative-interior controlled trace with zero endpoint extension,
assembles the independent trace $L^{2}$ metric, and realizes
$`P_{h}\widehat y_{h}+\ell_{0,h}+L_{D,h}u_{h}`$. Focused contracts distinguish
fixed and controlled values, test state/control residual and objective
pullbacks, metric apply/inverse, changed fixed data, and exact diagnostics for
overlapping or absent boundary regions. Sobolev metrics on this partial trace,
alternate interface policies, hanging-node trace relations, and trace boxes
remain P5.4 follow-ups.

**Implemented Chapter 5.11.2 slice:**
`make_l2_dirichlet_laplace_control_problem()` declares
$Y=H^{2}(\Omega)\cap H^{1}_{0}(\Omega)$,
$T=-\Delta:Y\rightarrow L^{2}(\Omega)$, residual codomain $Y^{\ast}$, and
multiplier space $Y$ for the complete-boundary $L^{2}(\Gamma)$ Dirichlet
control. Its selected discrete space
$U_{h}=\mathrm{tr}_{\Gamma}V_{h}\subset H^{1/2}(\Gamma)$ permits the equivalent
lifted $H^{1}$ Galerkin solve. The manifest records the continuous parent,
conforming subspace, equivalence, normalized Laplacian, and exclusions; it
does not claim support for facewise or discontinuous controls or a general
transposition lowerer.

The corresponding discrete conormal is the lifting pullback of the adjoint
residual, not a pointwise `FE_Q` normal derivative. The focused contract
verifies

```math
\langle q_{h},v_{h}\rangle
=a(L_{D,h}v_{h},p_{h})
-(y_{h}-z_{d},L_{D,h}v_{h})_{\Omega}
```

and the reduced stationarity covector
$\beta M_{\Gamma,h}u_{h}-q_{h}$. This fixes the sign from source equations
(5.171)–(5.174), rather than the inconsistent plus sign printed in Remark
5.18. `DirichletControlLiftingModel::discrete_conormal_covector()` and
`nmopt.dealii.l2_dirichlet_transposition` provide the independent contract.

**Implemented remaining Section 5.11 slices:** the Section 5.11.1
$H^{1/2}(\Gamma)$ metric is the minimum-volume-$H^{1}$ extension norm on
$U_{h}=\mathrm{tr}_{\Gamma}V_{h}$, realized by the Schur complement of
$M_{\Omega,h}+K_{\Omega,h}$. Register both source options: fractional control
loss with $L^{2}$ state tracking, and boundary $L^{2}$ control loss with
$H^{1}$ state tracking. Section 5.11.3 uses
$M_{\Gamma,h}+K_{\Gamma,h}$ for its tangential $H^{1}(\Gamma)$ loss and
metric. Loss and search metric remain distinct components in every case, and
none of these registrations has a trace box. `TraceHhalfMetric` applies the
Schur complement through one interior minimum-extension solve and applies its
inverse through one full-volume $H^{1}$ solve; it does not assemble a dense
fractional matrix. `DirichletControlLiftingModel` independently selects the
state tracking, control loss, and search metric, and assembles projected
tangential gradients for Section 5.11.3. The
`nmopt.dealii.section_5_11_compilation` contract verifies all three metric
round trips, reduced Taylor tests, stationarity decompositions, and distinct
structured manifests.

### P5.5 — Add regularised state-observation constraints with KKT provenance

**Status:** deferred by the accepted scope; this is the roadmap item for book
Section 5.12. Its measure-valued state-constraint problem remains excluded.

**Motivation:** C5.12 constrains the state, so it cannot be represented by the
current control-only projection capability. The unregularised problem has
measure multipliers and is outside the present first-order executable scope.

**First supported formulation:** declare a state-control observation

```math
w=O_{c}(y,u)=y+\lambda u,
\qquad \lambda>0,
```

with box constraints on $w$. The compiler must record the selected
Lavrentiev regularisation and may introduce the transformed decision
$`v=(\lambda I+S)u`$ only through declared maps and solves. It must expose the
$L^{2}$ multipliers, complementarity residual, and regularised adjoint relation.

**Do not do:** represent the original state constraint as a cellwise control
box; conceal the control-to-state inverse inside an optimiser; or present an
$L^{2}$ multiplier as the multiplier of the unregularised state constraint.

**Done when:** a feasible manufactured regularised problem satisfies primal
feasibility, dual feasibility, complementarity, stationarity, and a reduced
derivative test. A continuation test must record the sequence of positive
regularisation parameters. The original measure-multiplier problem remains a
capability diagnostic until its own functional-analytic contract exists.

### P5.6 — Generalize to mixed equation blocks and register a Stokes vertical slice

**Status:** deferred by the accepted scope; this is the roadmap item for book
Section 5.13 and the Chapter 6 Stokes examples.

**Motivation:** C5.13 needs vector fields, multiple equation blocks, a
pressure gauge, inf-sup stability, traction, and—later—boundary multipliers.
This is the concrete Chapter 5 instance of P4.2 item 1.

**Declare and implement:**

- General reduced DTO support for multiple state and test blocks and an
  externally supplied coupled state/adjoint solve. The reduced decision port
  remains explicit; do not make it infer a velocity-control problem.
- Vector-valued spaces and reusable viscous, pressure-gradient, divergence,
  traction, curl, and boundary-multiplier terms.
- Requirement policies for the pressure gauge, velocity/pressure FE pair,
  inf-sup assumption, traction trace, and controlled-boundary relation.

**First registered target:** serial steady Stokes with distributed vector
force control, velocity $L^{2}$ tracking, homogeneous Dirichlet and declared
traction boundaries, and a pressure mean policy. Its adjoint is the exact
coupled transpose and its reduced covector is the adjoint velocity plus the
regularised force covector.

**Follow-up target:** boundary velocity control with vorticity observation,
surface $H^{1}$ regularisation, and a declared boundary multiplier. It must
reuse the mixed equation and trace components rather than introduce a
named Stokes-control application type.

**Done when:** the mixed residual has blockwise JVP/VJP pairing tests, the
chosen FE pair/policy is in the manifest, the state and adjoint satisfy their
block residuals, and each target passes a coupled reduced Taylor test.

## Chapter 6 feature requests

The Chapter 6 guides separate reusable numerical-method infrastructure from
the source's scalar, boundary-control, and Stokes examples. These requests
extend formulation and solver layers without adding an optimiser or KKT class
for a named PDE. The accepted scope below supersedes the deferred GLS,
stabilised-Lagrangian, automatic-OtD, and Stokes variants; the separate
problem-library and benchmark-suite roadmaps own recipe and example
selection.

### P6.1 — Generalise reduced-space search strategies and line-search policies

**Status:** implemented for the selected scalar reduced DTO boundary.

**Motivation:** `ReducedGradientSolverT` supplies the first
steepest-descent/Armijo slice, while Section 6.3 also uses nonlinear CG,
Newton, and BFGS. They share reduced covectors, metrics, state/adjoint solves,
and reporting; they are not separate problem formulations.

**Declare and implement:**

- A typed primal search-direction protocol supplied by a reduced covector and
  declared inverse metric.
- Deterministic steepest-descent, nonlinear-CG (PR+ and Fletcher–Reeves
  updates with restart), and limited-memory BFGS policies. Store history with
  primal/dual layout checks and declared curvature/reset behaviour. An explicit
  quadratic-CG policy remains an extension candidate. A trust-region policy
  belongs to the extension track after a Hessian-vector service is available.
- A typed Hessian-vector service for the selected scalar DTO targets. Start
  with the exact linear-quadratic tangent-state/incremental-adjoint action;
  expose Newton or truncated Newton only for models that provide the required
  second-order actions. Do not imply generic nonlinear second-order support
  from a first-order residual JVP/VJP pair.
- Exact quadratic, Armijo, and later Wolfe line-search policies. Acceptance
  must use declared pairings and the actual projected displacement.
- A uniform report with accepted objective, covector/gradient norm, step,
  stop reason, line-search trials, and state/adjoint/metric solve counts.

**First registered target:** the existing one-state/one-decision linear DTO
path, using a mass metric and unconstrained L-BFGS. Keep projected steepest
descent as the only projected method until box transition/restart tests exist.

**Done when:** each selected direction has a descent test, each accepted trial
meets its declared inequality, and gradient/metric identities plus solve-count
accounting are verified. BFGS history must neither mix layouts nor silently
fall back after a failed curvature test. Every selected Hessian-vector path
also passes symmetry and finite-difference reduced-covector tests.

**Implementation closure:** the shared typed direction/result protocol now
covers steepest descent, Polak–Ribiere-plus and Fletcher–Reeves nonlinear CG,
and metric-aware limited-memory BFGS with declared layout, curvature, and
restart behaviour.
The explicit `ReducedHessianT` capability and exact linear-quadratic provider
support metric-preconditioned Newton. Exact quadratic, Armijo, and Wolfe
policies are selectable through the reduced solver and evaluate acceptance
using the actual trial displacement, including projected trials. The selected
scalar one-state/one-decision DTO target is verified across the neutral,
deal.II, sanitizer, and release profiles; projected steepest descent remains
the only projected direction policy in this slice.

#### P6.1 extension ladder

The status above closes the selected scalar P6.1 slice; it does not claim that
every alternative described in Chapter 3 has been implemented. The following
extensions are intentionally recorded so that a later implementation does
not silently broaden the selected slice:

1. **Near-term iterative closure:** relative-gradient, objective-change, and
   step-size stopping policies, final accepted state/adjoint/covector
   reporting, and Fletcher–Reeves are implemented. The remaining candidate is
   either classical quadratic CG or an exact-search PR+ equivalence test.
2. **Globalization extension:** add a matrix-free trust-region policy with a
   quadratic model, actual/predicted reduction ratio, radius update, and
   acceptance diagnostics. It is a separate globalization boundary, not a
   line-search option.
3. **Backend and benchmark parity:** provide the selected exact
   linear-quadratic reduced-Hessian action for the deal.II scalar target and
   exercise the iterative policies in the selected Chapter 6 benchmark
   harness.
4. **Lower-priority generic second order:** support nonlinear DTO models only
   when they provide an explicit Lagrangian/reduced Hessian-vector action,
   incremental-adjoint action, or a declared Gauss–Newton approximation.
   First-order JVP/VJP ports must not imply this capability.
5. **Lower-priority projected directions:** add projected Fletcher–Reeves or
   PR+, projected L-BFGS, and active-set-aware trust-region steps only after
   projection transition, restart, and metric-coupling tests exist. This does
   not authorize generic continuous-control box semantics.

Full BFGS, quadratic/cubic interpolation line searches, and other textbook
variants remain optional alternatives rather than prerequisites for the
selected framework slice.

#### Hessian and Newton boundary

P6.1 has three separable layers:

1. A **Newton optimizer** is a generic consumer. It receives a reduced
   covector and a service that applies $`H(u)w`$, solves the reduced Newton
   system, possibly approximately, and performs globalization. It contains no
   PDE-specific formulas.
2. A **Hessian-vector provider** is a model capability. For the selected
   linear-quadratic DTO problem, it is exact: a tangent-state solve and an
   incremental-adjoint solve produce

   ```math
   \begin{aligned}
   A\,\delta y &= B w,\\
   A^{\ast}\,\delta p &= W\,\delta y,\\
   H w &= \beta N w+B^{\ast}\delta p.
   \end{aligned}
   ```

   This is enough to run Newton or truncated Newton on E6.5.1 without
   assembling a dense reduced Hessian.
3. **Generic nonlinear second-order support** means producing that service
   for arbitrary nonlinear residuals and objectives. It requires derivatives
   of the residual Jacobian and objective derivative, or an explicitly
   supplied incremental-adjoint/Hessian action. The current first-order
   residual JVP/VJP and objective-derivative ports do not provide this
   information automatically.

The selected scope therefore includes the generic Newton consumer and the
linear-quadratic Hessian provider first. A nonlinear model may use Newton only
when it declares the stronger Hessian-vector capability; the optimizer must
reject rather than approximate a missing second-order action with finite
differences unless a separate finite-difference policy is explicitly selected.

### P6.2 — Record formulation provenance and execute supplied OtD systems

**Status:** selected only as a supplied-OtD execution interface. Stabilization,
GLS, stabilized Lagrangians, and automatic OtD derivation are excluded.

**Motivation:** Section 6.2 shows that OtD, DtO, and GLS-stabilised variants
can be distinct discrete systems even when they share a continuous PDE. The
compiler must make its selected formulation reviewable.

**Declare and implement:**

- A formulation policy naming DTO or supplied OTD, plus state/adjoint trial
  and test spaces, pairings, quadrature, and provenance.
- A separately owned supplied-OtD optimality-system product whose state,
  adjoint, and control-stationarity weak blocks are provided by the
  application author. The executor may assemble, linearize, solve, and
  report those blocks, but it does not derive them from a strong PDE.
- An explicit DTO/OTD sign, space, and comparison boundary. A supplied OTD
  product is not a `ReducedDTO` instance and is not presented as an exact
  discrete derivative unless equivalence is established.

**First registered target:** a small scalar supplied-OtD optimality system
whose weak state, adjoint, and control blocks are explicit and can be compared
with the corresponding DTO target. The target should exercise the execution
interface, not introduce a PDE-specific OTD class.

**Done when:** the manifest distinguishes DTO and supplied OTD, the supplied
blocks pass their value/JVP/VJP or block-linearization tests, and a
mismatched-adjoint-space request cannot be reported as a DTO optimum.

### P6.3 — Add reusable equality-constrained quadratic KKT products

**Status:** selected after P6.1 and the supplied-OtD interface; required by
the selected PDAS path.

**Motivation:** All-at-once OCPs are equality-constrained quadratic programs
with a symmetric indefinite KKT operator. PDAS repeatedly solves related KKT
subproblems. A generic block product is needed before problem-specific
preconditioners.

**Declare and implement:**

- A formulation service for $Q$, $D$, and the KKT action
  $\begin{bmatrix}Q&D^{\mathsf T}\\D&0\end{bmatrix}$, assembled from
  objective derivatives and residual JVP/VJP actions.
- Explicit primal, residual-multiplier, and block-layout descriptors; KKT
  multiplier sign conversion to the framework adjoint belongs in the manifest.
  Record rank and kernel-positivity assumptions needed for a nonsingular KKT
  product or return a formulation diagnostic.
- Krylov policies for symmetric-indefinite MINRES and nonsymmetric GMRES, with
  formulation/layout compatibility and residual diagnostics.

**First registered target:** serial scalar linear-quadratic DTO with a
volume-control matrix-free/assembled-equivalence test. Do not require a
scalar PDE name or identify control and state mass matrices in the generic
product.

**Done when:** KKT block action and transpose signs pass, its solution agrees
with reduced DTO on the same target, and it reports feasibility, stationarity,
multiplier conversion, and solver termination independently.

### P6.4 — Compose block preconditioners from declared approximate solves

**Status:** conditional. Do not implement this as a Chapter 6 prerequisite.
Activate a bounded scalar preconditioner only if the selected all-at-once
application benchmark cannot be run or interpreted with the available direct
or basic serial solves.

**Motivation:** Chapter 6 uses diagonal, triangular, and constraint
preconditioners built from mass, PDE, Schur-complement, multigrid, and Uzawa
actions. These must remain reusable approximate operator services.

**Declare and implement:**

- Block-diagonal Schur, Bramble–Pasciak triangular, and constraint
  preconditioner compositions, each with its required symmetry/inner-product
  property and compatible outer Krylov method.
- Approximate mass and PDE inverse actions with fixed work/tolerance policies
  and nested-solve diagnostics. Variable inner Krylov work must select a
  flexible outer method.
- Parameter-robustness and mesh-scaling claims as benchmark metadata, never
  as an unchecked property of a preconditioner name.

**Conditional first target:** the scalar all-at-once product of P6.3 with one
positive-definite block-diagonal preconditioner, mass lumping or a declared
stationary mass approximation, and a fixed-cycle serial multigrid state
approximation. Do not add the full family of preconditioners unless a
benchmark requires it.

**Done when:** preconditioner layouts and symmetry restrictions are validated,
apply diagnostics are deterministic, and a mesh/regularisation sweep records
outer and inner iterations separately. Add Stokes/Uzawa after P5.6 supplies
the mixed state operator and pressure policy only if the deferred mixed-PDE
scope is later reopened.

### P6.5 — Add typed complementarity, selection, and PDAS services

**Status:** selected after P6.3, initially for cellwise $L^{2}$ control boxes.
Continuous-control and other bound semantics remain separate optional
extensions.

**Motivation:** Section 6.8 treats control boxes and regularised mixed
state-control bounds with semismooth primal-dual active sets. The current
control projection cannot express a multiplier, active selection, or a
sequence of KKT subproblems.

**Declare and implement:**

- A primal constraint observation with a compatible dual multiplier,
  complementarity residual, selected primal/dual representation for active
  classification, and active-set restriction/prolongation operators.
- A PDAS/semismooth-Newton service that constructs KKT equality subproblems
  through P6.3, updates sets, and reports stable-set plus full-KKT
  convergence.
- Do not include the regularised observation $O_{c}(y,u)=y+\varepsilon u$
  in the selected path. It belongs to deferred P5.5/Section 5.12 work, and
  the original measure-multiplier state constraint remains out of scope.

**First registered target:** cellwise-discontinuous distributed control with
two-sided boxes and an $L^{2}$ multiplier, beginning with an inactive-box case
where PDAS agrees with the unconstrained KKT solution.

**Done when:** every iteration exposes primal/dual feasibility,
complementarity, stationarity, active-set change, and KKT diagnostics; an
active manufactured case stabilises at correct bounds; and no continuous
control coefficient or multiplier is classified pointwise without its
declared conversion policy.

## Current next-agent sequence

The C5.6 Neumann composition is complete. P5.1 is acceptance-complete. The
remediation review documents are static evidence; this roadmap is the status
ledger. The C1/C2 remediation slices are implemented and acceptance-complete:
one checked resolved request, truthful P5.1 data placements and boundary
realization, plan-owned scalar residual/data and objective/service lowering,
resolved-decision manifest records, common realized-map/space records, and
the documented differential and verification gates. P5.2–P5.4 are likewise
acceptance-complete for their registered bounded slices. The review's “Work
unit” names are cross-references to those documents, not a separate roadmap
vocabulary.

### Completed remediation and acceptance

- **C1/C2 work unit 1:** [one resolved compilation request and binding boundary](review/chapter-5/c1-c2-preparation-remediation-review.md#work-unit-1--one-resolved-compilation-request-and-binding-boundary), including exact public binding/session diagnostics.
- **P5.1 work unit 1:** [truthful P5.1 data spaces and regions](review/chapter-5/p5.1-remediation-review.md#work-unit-1--truthful-p51-data-spaces-and-regions), carried into the resolved request.
- **P5.1 work unit 2:** [typed boundary and conormal selection](review/chapter-5/p5.1-remediation-review.md#work-unit-2--typed-p51-boundary-and-conormal-selection), with shared boundary, orientation, and trace-realization vocabulary and no target-specific policy enums.
- **C1/C2 work units 2–3:** [plan-owned scalar residual/data assembly](review/chapter-5/c1-c2-preparation-remediation-review.md#work-unit-2--plan-owned-scalar-residual-and-data-assembly) and [objective/service recombination](review/chapter-5/c1-c2-preparation-remediation-review.md#work-unit-3--plan-owned-objective-and-service-recombination), including an independently varied recombination.
- **P5.2–P5.4 acceptance closure:** typed trace and negative-metric selections, target-data assumptions, control-boundary realization, typed transposition, fixed-boundary coverage, point/flux checks, closed registration matching, and realized transformed-observation dimensions. The `debug-neutral`, `debug-dealii`, and `sanitize-neutral` gates pass.
- **Supporting evidence:** the closed request and pre-construction decision own target selection, typed regions and requirements, component inventories, and map skeletons; finalization adds only realized service facts, and manifest construction is a pure projection. Full compatibility remains stable under declaration-order and display-prose changes. Realized map/space records cover weighted traces, normal flux, transformed state observations, and the baseline boundary trace. The outward normal, face-quadrature transpose, immutable physical-point evaluation, and current exclusions remain explicit.
- **P6.1 reduced-space solver slice:** the selected direction, Hessian/Newton, line-search, reporting, and configurable-composition contracts are implemented for the scalar reduced DTO boundary. The final verification pass covers the neutral, deal.II Debug, sanitizer, and deal.II Release profiles.

### Future implementation sequence

P4.1 and P4.2 remain ignored for the current ordered implementation run. The
selected P6.1 scalar reduced DTO boundary is complete; future work begins
with P6.2:

1. **P6.2:** implement the supplied-OtD execution interface separately from
   `ReducedDTO`; do not add GLS, stabilization, or automatic continuous
   adjoint derivation.
2. **P6.3/P6.5:** implement the scalar KKT product and diagnostics, then
   complementarity/PDAS for the selected cellwise box representation.
3. **Chapter 6 benchmarks:** use the separate [Chapter 6 benchmark suite roadmap](chapter-6-benchmark-suite-roadmap.md) to run E6.5.1, E6.5.2, and E6.9.1/E6.9.2. Activate bounded P6.4 preconditioner work only if E6.7.1 is selected and basic serial solves are insufficient.
4. **P6.1 extension track:** after the required P6.2/P6.3/P6.5 gates and
   selected benchmark evidence, activate the [P6.1 extension ladder](#p61-extension-ladder)
   in order. Generic nonlinear second order and projected directions remain
   lower priority than the selected scalar formulation and KKT/PDAS paths.

Follow the [Stage B routing protocol](review/pre-ch5-ch6/README.md) for each
future gate. The next conditional gate is P6.2; do not silently activate S1
or any of the remaining Stokes, measure-constraint, stabilization,
automatic-OtD, or broad continuous-bound work. Those remain outside the
current ordered run.
