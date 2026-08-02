# Implementation roadmap and agent handoff

## Purpose

This is the ranked continuation plan after the v0 executable contract and the
first serial deal.II lowerer. It is intentionally dependency-ordered: finish
the useful vertical slice before broadening the semantic language, and broaden
the semantic language before adding advanced PDE variants.

The command below is the baseline check before and after every task:

~~~bash
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build
ctest --test-dir build --output-on-failure
~~~

## Current handoff state

The following pieces exist and are tested:

| Layer | Existing artifact | Meaning |
|---|---|---|
| Typed algebra | `include/nmopt/contract/layout.hpp` | `PrimalBlockT` and `CovectorBlockT` are distinct typed wrappers, even when a backend uses one vector storage type. |
| V1 semantic graph | `include/nmopt/semantic/v1/{types,validation,reference_specs}.hpp` | Deal.II-free selected graph, explicit pairings, and structural/policy diagnostics. |
| V1 compiler | `include/nmopt/compiler/v1/{compiled_problem,dealii_compiler}.hpp` | Backend-generic compiled package, manifest, and registered assembled volume slice. |
| Operator contract | `include/nmopt/contract/executable_model.hpp` | Residual, JVP, VJP, objective, and objective derivative. |
| DTO workflow | `include/nmopt/contract/reduced_dto.hpp` | One state block, one control block, one test block, externally supplied state/adjoint solves. |
| Reference oracle | `include/nmopt/reference/linear_quadratic_model.hpp` | Dense linear-quadratic model used to test signs and derivatives independently of deal.II. |
| deal.II backend | `include/nmopt/dealii/serial_backend.hpp` | Serial Vector backend for the backend-parametric contract. |
| deal.II lowerer | `include/nmopt/dealii/scalar_diffusion_reaction.hpp` | Assembled scalar `FE_Q` diffusion-reaction state, `FE_DGQ(0)` volume control, full-domain tracking, homogeneous Dirichlet data, and DTO solves. |
| deal.II metric | `include/nmopt/dealii/mass_metric.hpp` | One-block sparse-mass $L^{2}$ Riesz action with serial CG inverse apply. |
| deal.II constraint | `include/nmopt/dealii/cellwise_box_constraint.hpp` | `FE_DGQ(0)` coefficientwise box projection for the declared `l2_cellwise` metric. |
| Reduced solver | `include/nmopt/solvers/reduced_gradient.hpp` | Backend-parametric unconstrained and projected Armijo method over `ReducedDTOT`, `MetricT`, and optional `ConstraintT`. |
| Tests | `tests/reduced_dto_contract.cc` and `tests/dealii_diffusion_contract.cc` | Pairing, JVP, VJP, residual finite difference, state solve, reduced derivative, deal.II metric, and unconstrained/projected Armijo checks. |

The implementation is deliberately not a general compiler yet. The deal.II
class is a concrete lowerer/reference slice, not the future public
`ProblemSpec` interface.

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

**Implemented:** `dealii_backend::MassMetric` provides the generic one-block
sparse-mass realization of `MetricT<SerialBackend>`. Its inverse action uses
serial CG with declared iteration and relative/absolute tolerances, and
`ScalarDiffusionReactionModel::control_l2_metric()` is the control-space
compiled-metric factory.

**Primary files:** add `include/nmopt/dealii/mass_metric.hpp`; minimally extend
the current lowerer with a metric factory or compiled-object registry.

**Done when:**

- applying the inverse followed by apply recovers a random covector;
- a metric direction satisfies
  $`\langle M_{u} g,\delta u\rangle=\langle j',\delta u\rangle`$;
- the test uses actual deal.II vectors and a nontrivial control mass matrix.

### P0.2 — Add a generic reduced Armijo gradient solver — completed

**Why second:** This closes the current vertical slice: the existing state
solve, adjoint solve, reduced covector, and new metric become a real
optimization loop. The solver must consume only `ReducedDTOT` and `MetricT`.

**Implement:**

- A backend-parametric reduced gradient solver with Armijo backtracking.
- Objective history, gradient norm in the selected metric, state/adjoint solve
  counts, and clear stopping reasons.
- An unconstrained path first; no PDE type checks or casts.

**Implemented:** `solvers::ReducedGradientSolverT` evaluates every initial and
trial control with `ReducedDTOT::evaluate`, converts the DTO covector through
the supplied metric, and applies an unconstrained Armijo update. Its result
records objective and metric-gradient-norm histories, accepted iterations,
line-search trials, state/adjoint solve counts, and explicit stopping reason.

**Primary files:** add `include/nmopt/solvers/reduced_gradient.hpp` and a
deal.II integration test.

**Done when:**

- the objective decreases monotonically under the line-search criterion;
- the norm of the reduced covector/gradient reaches a configured tolerance on
  a manufactured linear-quadratic case;
- the same solver passes the dense reference model test unchanged.

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

**Implemented:** `dealii_backend::CellwiseBoxConstraint` is the serial
coefficientwise `ConstraintT<SerialBackend>` realization. The lowerer exposes
only `FE_DGQ(0)` control factories for scalar constants or exact-layout
coefficient vectors. The constraint-qualified `ReducedGradientSolverT`
projects every trial in `l2_cellwise` and stops on the metric norm of the
projected-gradient residual.

**Done when:** a test reaches a bound on a manufactured problem, preserves
feasibility every iteration, and satisfies a discrete projected-stationarity
criterion.

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

**Implemented (v1):** `semantic::v1::ProblemSpec` and
`SemanticValidator` describe and validate the selected graph without backend
objects. `compiler::v1::DealiiCompiler` appends lowerability and formulation
diagnostics through a small explicit lowerer registry, then produces a
separately owned compiled executable. The v0 direct
`ScalarDiffusionReactionModel` remains unchanged as the reference path; v1
constructs its own instance from the semantic declaration, so both paths can
be compared without overwriting v0.

**Primary files:** `include/nmopt/semantic/v1/{types,validation,reference_specs}.hpp`,
`include/nmopt/compiler/v1/{compiled_problem,dealii_compiler}.hpp`,
`docs/semantic-v1-compiler.md`, and focused semantic/deal.II contract tests.

**Done when:** constructing the current problem through a `ProblemSpec` produces
the same residual, objective, and reduced derivative as the direct deal.II
reference setup.

### P1.2 — Generalize fixed essential conditions through reconstruction

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

### P1.3 — Add subdomain observations and data projection policy

**Why:** It exercises the observation/loss port without changing PDE physics.
It also forces the compiler to make region and data semantics real.

**Implement:**

- Named material/subdomain regions and an observation mask/restriction.
- A tracking mass/load assembled only on the observation region.
- Explicit data rule: analytic `Function` evaluated at quadrature or a declared
  FE projection/interpolation. No implicit nodal interpretation.

**Done when:** changing only the observation region changes the adjoint RHS
and objective, not $`A_{h}`$, $`B_{h}`$, the state solution, or solver code.

### P2.1 — Add Neumann control and boundary tracking

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

### P2.2 — Support pure Neumann with the selected mean-constraint policy

**Why:** This is a correctness feature, not a convenience flag. It affects
state uniqueness, adjoints, observation meaning, metrics, and preconditioners.

**Default:** augment state and adjoint systems with one mean-zero
constraint/multiplier. Check compatibility of every volume and boundary load,
including controls.

**Done when:** incompatible data is rejected; compatible problems report the
gauge in their compilation manifest; state and adjoint are invariant under no
unrecorded pinning convention.

### P2.3 — Add $H^{1}$ regularisation and $H^{1}$ search geometry separately

**Why:** The current code makes the distinction possible but does not yet
exercise it.

**Implement in two different tasks or commits:**

1. $H^{1}$ **regularisation**: a new control loss and likely an `FE_Q` control
   space; it changes $`J_{u}'`$.
2. $H^{1}$ **metric**: a Riesz map/inverse used only for
   $G^{-1}j'$; it does not change the objective or adjoint.

The metric needs a positive zero-order term or an explicit boundary/mean
policy to be invertible. An $H^{-1}$-type metric remains unsupported until
its exact Hilbert space and discrete operator are stated.

### P3.1 — Add coefficient identification and nonlinear first-order actions

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

### P3.2 — Add Dirichlet control through an explicit lifting

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

### P4.1 — Add the fixed-step temporal compiler

**Default:** a global residual for a fixed-step backward-Euler heat equation,
with the complete trajectory available for its exact transpose.

**Prerequisites:** stable reconstruction semantics, compilation manifests,
and state/adjoint solve diagnostics.

**Do not do:** hide a forward time integrator behind a separately coded
backward adjoint, or introduce adaptive steps/events before replay and
differentiated-control policies exist.

### P4.2 — Generalize algebra, execution, and formulations

This final group has value, but should not block the preceding vertical
slices:

1. Multiple state/equation blocks, mixed fields, and Petrov–Galerkin spaces.
2. Serial matrix-free equivalence, then distributed Trilinos/PETSc vector
   backends with ownership/ghost policies.
3. OTD formulation builder with explicit provenance and DTO comparison tests.
4. Second-order Lagrangian Hessian-vector contract, then nonlinear KKT
   Newton/SQP and active-set methods.
5. Fractional metrics, point/flux observations, nonmatching meshes, and shape
   optimization only with their own explicit policies.

## Suggested next-agent sequence

For the next agent, take **P1.2 only** unless explicitly asked for a broader
change:

1. Read this roadmap, the interface specification, the v1 semantic/compiler
   record, the executable contract, and the deal.II lowerer document.
2. Run the baseline `CTest` command.
3. Add an explicit fixed-data lifting/reconstruction to the v1 path without
   altering the v0 homogeneous reference model.
4. Preserve physical-field residual, observation, JVP, and VJP semantics;
   extend the compiler manifest and diagnostics for the selected lifting.
5. Prove a nonzero manufactured Dirichlet state with a reduced Taylor test,
   then run the baseline tests again.

This keeps the public contracts stable and makes each step reviewable.
