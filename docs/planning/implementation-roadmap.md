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

Run `debug-dealii` as well for compiler, lowerer, backend, or numerical
changes. The [build convention](../../conventions/build.md) owns the complete
profile matrix and cache-recovery guidance.

## Current handoff state

Stage B common stabilization is complete on
`codex/refactor-ch5-ch6-readiness`. Batches R0
(`RF-020`), R1 (`RF-006`, `RF-009`, and `RF-016`), R2a (`RF-001` and the
relevant `RF-006` characterization), R2b (`RF-002` through `RF-005` and the
relevant `RF-006` cases), R2c (the current factual defect in `RF-008` plus
`RF-012`), and R3 (`RF-016` through `RF-019`) are complete. P4.1 and P4.2 are
ignored for the current ordered implementation run because their scope is too
broad. P5.1 and its conditional C1 (`RF-008` through `RF-013`) and C2 gates
are complete. The P5.2 full-domain $H^{1}_{0}$ state-observation unit is
complete, and the weighted-boundary-trace observation is now complete as its
own review boundary. The separately specified $H^{-1}$ metric is the selected
next vertical slice.

The following pieces exist and are tested:

| Layer | Existing artifact | Meaning |
|---|---|---|
| Typed algebra | `include/nmopt/contract/layout.hpp` | `PrimalBlockT` and `CovectorBlockT` are distinct typed wrappers, even when a backend uses one vector storage type. Their block storage is read-only after construction; checked algebraic updates preserve the declared dimensions. |
| V1 semantic graph | `include/nmopt/semantic/v1/{types,validation,reference_specs}.hpp` | Deal.II-free selected graph with safe incomplete states, whole-graph closure checks, explicit two-sided pairings, structural/policy diagnostics, and ID-based reference deltas. |
| V1 compiler | `include/nmopt/compiler/v1/{compiled_problem,dealii_compiler,dealii_scalar_plan}.hpp` | Backend-generic compiled package and structured manifest; fixed-reconstruction, subdomain-tracking, $H^{1}_{0}$ state-tracking, and general scalar/Robin graphs use stable-ID resolution plus stored component handlers, while specialized registrations retain bounded target strategies listed in the [v1 capability table](../implementation/v1/semantic-compiler.md#registered-capabilities). |
| Operator contract | `include/nmopt/contract/executable_model.hpp` | Residual, JVP, VJP, objective, and objective derivative. |
| DTO workflow | `include/nmopt/contract/reduced_dto.hpp` | One state block, one decision block (control or parameter), one test block, externally supplied state/adjoint solves. |
| Formulation solves and lifetime | `include/nmopt/contract/linear_solve.hpp`, `include/nmopt/dealii/serial_spd_solver.hpp`, and `include/nmopt/compiler/v1/dealii_types.hpp` | Typed state/adjoint solve reports, one shared serial SPD policy/service for symmetric targets, recorded direct and exact-transpose solves for the P5.1 nonsymmetric target, an owned static-mesh compilation session, and detached reduced services that retain executable/session lifetime. |
| Reference oracle | `include/nmopt/reference/linear_quadratic_model.hpp` | Dense linear-quadratic model used to test signs and derivatives independently of deal.II. |
| deal.II backend | `include/nmopt/dealii/serial_backend.hpp` | Serial Vector backend with a checked conversion from contract dimensions to the native deal.II size type. |
| Direct deal.II v0 lowerer | `include/nmopt/dealii/scalar_diffusion_reaction.hpp` | Preserved assembled scalar `FE_Q` diffusion-reaction reference with `FE_DGQ(0)` volume control, full-domain tracking, homogeneous Dirichlet data, and DTO solves. |
| deal.II metrics | `include/nmopt/dealii/mass_metric.hpp` | One-block sparse SPD Riesz actions for the registered volume, boundary, trace, and parameter layouts, with serial CG inverse apply and an operator-bound realization witness. |
| deal.II constraints | `include/nmopt/dealii/{cellwise,facewise}_box_constraint.hpp` | Coefficientwise boxes coupled to the actual positive-diagonal cellwise-volume or facewise-boundary $L^{2}$ metric realization, never to its display string. |
| Reduced solver | `include/nmopt/solvers/reduced_gradient.hpp` | Backend-parametric unconstrained and projected Armijo method over `ReducedDTOT`, `MetricT`, and optional `ConstraintT`. |
| Build/test workflow | `CMakePresets.json` and `CMakeLists.txt` | Explicit neutral/deal.II Debug, neutral sanitizer, and deal.II Release profiles; requested dependency failures; target-scoped warnings; and labeled, time-bounded scenarios. |
| Tests | `tests/{reduced_dto_contract,semantic_v1_contract,dealii_diffusion_contract}.cc` | Three binaries expose twenty-eight independently named, labeled, and time-bounded CTest scenarios: five dense/backend contract cases, seven semantic graph/resolution/lowering-plan cases, and sixteen deal.II compiler/lowering/adapter cases. Negative checks identify exact diagnostics or contract failures, including block/layout, graph-closure, coefficient shape, boundary partition, binding, solve-policy, manifest-realization, lifetime, projection-coupling, observation topology/region, and native-size invariants; the canonical deal.II scenario includes a hand-integrated weak-form oracle separately from its compiled/direct wiring comparison, the energy-observation scenario has an independent polynomial $H^{1}_{0}$ oracle, and the weighted-trace scenario compares against an exact constant-weight scaling oracle while proving the residual and metric unchanged. |

The public v1 semantic path is deliberately not a general component compiler
yet. It resolves valid graphs by stable ID and has one bounded scalar
component-planning path; specialized graphs still select one of a bounded set
of target-specific implementations. The direct
deal.II v0 class remains a concrete reference lowerer rather than a public
problem hierarchy.

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

The metric needs a positive zero-order term or an explicit boundary/mean
policy to be invertible. An $H^{-1}$-type metric remains unsupported until
its exact Hilbert space and discrete operator are stated.

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
3. OTD formulation builder with explicit provenance and DTO comparison tests.
4. Second-order Lagrangian Hessian-vector contract, then nonlinear KKT
   Newton/SQP and active-set methods.
5. Fractional metrics, point/flux observations, nonmatching meshes, and shape
   optimization only with their own explicit policies.

## Chapter 5 feature requests

The Chapter 5 implementation guide records the mathematical variants and the
required composition boundary for each request. The requests below are
ordered by reusable capability, not by the source chapter's presentation
order. They must not be implemented as a hierarchy of named textbook problem
classes.

### P5.1 — Compose general scalar elliptic volume and Robin boundary terms — completed

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

### P5.2 — Add energy-volume and weighted-trace observations, then the selected $H^{-1}$ metric separately

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

**Implemented observation slices:**
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
The separate $H^{-1}$ metric remains pending.

### P5.3 — Add normal-flux and point-sensor observations through an explicit strong/very-weak policy

**Motivation:** C5.8 and C5.10 are not ordinary restriction/trace
observations. Their adjoints use boundary data or Dirac sources with lower
regularity. The $L^{2}$ Dirichlet-control variant in C5.11 uses the same
transposition principle.

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

**First registered target:** scalar Dirichlet Laplace or reaction-diffusion
on a convex or declared $C^{2}$ serial domain. Choose one explicit discrete
point-evaluation policy and one flux policy; reject every alternative until
it has its own lowerer. The corresponding adjoint must use the declared
transpose/very-weak formulation, not an undeclared nodal Dirac approximation.

**Done when:** sensor and flux observation derivatives pass value/JVP/VJP
tests; the low-regularity adjoint passes its dual residual test; and the
compiled product records all regularity, orientation, and evaluation policies.

### P5.4 — Generalize Dirichlet-control transformations and trace metrics

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

**First registered target:** one partial scalar controlled boundary disjoint
from one fixed nonzero Dirichlet boundary, with a complete declaration of the
corner policy. Start with an $L^{2}$ trace metric; add fractional and
tangential metrics as separate follow-up targets rather than silently using a
mass matrix for all three.

**Done when:** the physical reconstruction and both state/control pullbacks
pass VJP tests under changed fixed and controlled data; each metric passes its
own pairing test; and incomplete or ambiguous boundary declarations are
rejected.

### P5.5 — Add regularised state-observation constraints with KKT provenance

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
for a named PDE. P5.6 remains the prerequisite for the Stokes portions.

### P6.1 — Generalise reduced-space search strategies and line-search policies

**Motivation:** `ReducedGradientSolverT` supplies the first
steepest-descent/Armijo slice, while Section 6.3 also uses nonlinear CG,
Newton, and BFGS. They share reduced covectors, metrics, state/adjoint solves,
and reporting; they are not separate problem formulations.

**Declare and implement:**

- A typed primal search-direction protocol supplied by a reduced covector and
  declared inverse metric.
- Deterministic steepest-descent, nonlinear-CG (selected update and restart),
  and limited-memory BFGS policies. Store history with primal/dual layout
  checks and declared curvature/reset behaviour. Add a trust-region policy
  only after a Hessian-vector service is available.
- Exact quadratic, Armijo, and later Wolfe line-search policies. Acceptance
  must use declared pairings and the actual projected displacement.
- A uniform report with accepted objective, covector/gradient norm, step,
  stop reason, line-search trials, and state/adjoint/metric solve counts.

**First registered target:** the existing one-state/one-decision linear DTO
path, using a mass metric and unconstrained L-BFGS. Keep projected steepest
descent as the only projected method until box transition/restart tests exist.

**Done when:** each direction has a descent test, each accepted trial meets
its declared inequality, and gradient/metric identities plus solve-count
accounting are verified. BFGS history must neither mix layouts nor silently
fall back after a failed curvature test.

### P6.2 — Record formulation, trial/test, and stabilisation provenance

**Motivation:** Section 6.2 shows that OtD, DtO, and GLS-stabilised variants
can be distinct discrete systems even when they share a continuous PDE. The
compiler must make its selected formulation reviewable.

**Declare and implement:**

- A formulation policy naming DTO, OTD, or declared stabilised-Lagrangian
  construction, plus state/adjoint trial and test spaces and quadrature.
- Stabilisation terms as residual/objective components with exact value, JVP,
  and VJP ownership. Include selected element residual, local parameter, and
  inflow/outflow policies in the manifest.
- A separately owned OTD optimality-system product. It is not a
  `ReducedDTO` instance unless a comparison proves discrete equivalence.

**First registered target:** scalar advection-diffusion DTO with one GLS
policy and full-volume tracking. Its DTO adjoint is the exact transpose of the
lowered residual. Add the strongly consistent stabilised-Lagrangian policy
only after its coupled derivatives are explicit.

**Done when:** the manifest distinguishes all formulations, stabilised
JVP/VJP and reduced Taylor tests pass, and a mismatched-adjoint-space request
cannot be reported as a DTO optimum.

### P6.3 — Add reusable equality-constrained quadratic KKT products

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

**First registered target:** the scalar all-at-once product of P6.3 with one
positive-definite block-diagonal preconditioner, mass lumping or a declared
stationary mass approximation, and a fixed-cycle serial multigrid state
approximation.

**Done when:** preconditioner layouts and symmetry restrictions are validated,
apply diagnostics are deterministic, and a mesh/regularisation sweep records
outer and inner iterations separately. Add Stokes/Uzawa after P5.6 supplies
the mixed state operator and pressure policy.

### P6.5 — Add typed complementarity, selection, and PDAS services

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
- The regularised observation $O_{c}(y,u)=y+\varepsilon u$ of P5.5 as a
  separate target. Its multiplier VJP contributes to state and control
  stationarity; the original measure-multiplier state constraint remains out
  of scope.

**First registered target:** cellwise-discontinuous distributed control with
two-sided boxes and an $L^{2}$ multiplier, beginning with an inactive-box case
where PDAS agrees with the unconstrained KKT solution.

**Done when:** every iteration exposes primal/dual feasibility,
complementarity, stationarity, active-set change, and KKT diagnostics; an
active manufactured case stabilises at correct bounds; and no continuous
control coefficient or multiplier is classified pointwise without its
declared conversion policy.

## Current next-agent sequence

P5.2 is the selected bounded Chapter 5 vertical slice. P4.1 and P4.2 remain
ignored for the current ordered implementation run. Continue as follows:

1. P5.1 is complete: its component plan and first registered deal.II target
   verify individual term actions, the nonsymmetric adjoint, Robin value/load
   contributions, boundary/shape diagnostics, and the reduced derivative.
2. Reuse the completed C1/C2 compiler and component-lowering boundaries for
   P5.2; do not reopen the broad P4.2 algebra/formulation group.
3. The full-domain `h1_state_restriction` and `weighted_boundary_trace`
   observation units are complete, including declared weight data, trace
   quadrature, and value/JVP/VJP coverage.
4. Implement the selected $H^{-1}$ metric as the next separate review
   boundary. State its control space, Riesz operator, boundary/mean policy,
   and inverse-solve tolerances before changing code.

Follow the [Stage B routing protocol](refactor/README.md) for each gate. Do not
run S1 before P6.1 reaches the front of the ordered implementation run.
