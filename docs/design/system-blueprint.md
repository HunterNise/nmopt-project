# System blueprint: from variational problem to executable optimisation

## Purpose

This guide is the shortest path to a working mental model of `nmopt`. It
connects the mathematical model, the semantic specification, the current v0
contracts, the serial deal.II reference lowerer, and the tests.

It is not a second source of authority: the [interface specification](interface-specification.md)
is normative for the intended semantic API; the
[v0 executable contract](../implementation/v0/executable-contract.md) and
[deal.II lowerer record](../implementation/v0/dealii-lowerer.md) define the
direct reference slice; and the
[v1 capability table](../implementation/v1/semantic-compiler.md#registered-capabilities)
owns exact compiler support and exclusions.

## The one-sentence model

`nmopt` does not model a “diffusion-control problem” or a “Neumann-control solver” as a class. It models a graph of independently meaningful maps:

```text
variables + fixed data
          │
          ├── transformations ──> physical fields
          │                         │
          ├── residual terms ──> equation block E(x) in Z* ──> state/adjoint solves
          │                         │
          └── observations ──> losses ──> objective J(x) ──> reduced covector
                                                               │
                                 metric + constraint ─────────┘
                                                               │
                                                         optimisation step
```

The reusable core is not the particular PDE. It is the ability to evaluate a map, apply its Jacobian to a primal direction, and pull a seed back through its adjoint:

```math
F(x), \qquad F'(x)\delta x, \qquad F'(x)^{\ast}q.
```

For a residual, the output is a covector in a *test-space dual*, so its VJP seed is a primal test-space vector. That distinction is the source of much of the design.

## Intended architecture versus present code

```text
PRESENT V1 SEMANTIC/COMPILER PATH

ProblemSpec ──> semantic validation + kind whitelist
                              │
                              └──> whole-graph predicate dispatch
                                           │
                                           └──> one target implementation
                                                       │
                                                       └──> executable ports

TARGET COMPONENT-LOWERING ARCHITECTURE (NOT YET IMPLEMENTED)

ProblemSpec ──> validated component graph ──> independently owned lowerers
                                                     │
                                                     └──> composed executable model


IMPLEMENTED V0 VERTICAL SLICE

LinearQuadraticModel ──────────────────────────────────> ExecutableModelT
                                                       ↗
ScalarDiffusionReactionModel (preserved hand-written deal.II reference)
      │
      └──> ReducedDTOT + supplied state/adjoint solves ──> reduced covector
```

The `ScalarDiffusionReactionModel` is a reference lowerer for one selected
finite-element problem. It is **not** the public `ProblemSpec` and must not
grow into a hierarchy of complete PDE/control combinations. The present v1
kind registry improves diagnostics, but it does not independently lower an
arbitrary combination of whitelisted components.

| Layer | Question answered | Authority | Representative code |
| --- | --- | --- | --- |
| Theory | What mathematical object is solved? | [Formalism](theoretical-formalism.md) | `LinearQuadraticModel` oracle |
| Semantic specification | Which components and ports should exist? | [Interface specification](interface-specification.md) | `semantic::v1::ProblemSpec` |
| Compilation policy | How do spaces, pairings, liftings, and execution become discrete? | [V1 semantic compiler](../implementation/v1/semantic-compiler.md) | `compiler::v1::DealiiCompiler` |
| Executable contract | What may algorithms call after lowering? | [V0 contract](../implementation/v0/executable-contract.md) | `include/nmopt/contract/` |
| Formulation | How do residual and objective become first-order operations? | [Interface specification](interface-specification.md) | `ReducedDTOT` |
| Backend/lowerer | How is one model assembled in deal.II? | [deal.II lowerer](../implementation/v0/dealii-lowerer.md) | `ScalarDiffusionReactionModel` |
| Verification | How do values and derivatives agree? | [Roadmap](../planning/implementation-roadmap.md) | `tests/*_contract.cc` |

## The vocabulary: component cards

The semantic layer is deliberately broader than the implemented v1 slice.
Each card says what the component owns and which kind of port it exposes. For
current implementations of these cards, consult the
[v1 capability table](../implementation/v1/semantic-compiler.md#registered-capabilities).

| Component | Owns | Communicates through |
| --- | --- | --- |
| `Region` | Named volume, boundary, interface, point set, or time set | Identity, dimension, relation |
| `Space` and `Pairing` | Field shape, topology, role, primal/dual pairing | Typed source and target ports |
| `VariableBlock` | State, control, parameter, flux, or auxiliary unknown | One primal space; feeds maps |
| `Data` | Fixed forcing, target, coefficient, lifting, or bound | Read-only ports; never a derivative block |
| `Transformation` | Reconstruction, lifting, parameterisation, restriction, transfer | Value, JVP, VJP |
| `ResidualTerm` | One physical contribution to one equation | Tested value, JVP, VJP |
| `EquationBlock` | Sum of residual terms and its test space | $E$, $E'\delta x$, $E'^{\ast}p$ |
| `Observation` | Map from physical variables to observation space | Value, JVP, VJP |
| `Loss` | One scalar penalty of an observation | Scalar value and output covector |
| `Objective` | Sum of loss compositions | $J$ and $J'$ |
| `Metric` | Algorithmic map $G:P\to P^{\ast}$ | `apply`, `inverse_apply` |
| `Constraint` | Feasibility, projection, normal cone, multipliers | Operations in a named metric |
| `RequirementPolicy` | A non-inferable trace, nullspace, point, or discrete-only choice | Validator metadata |
| `DiscretisationPolicy` | FE family, mesh relation, quadrature, lifting, execution | Input to lowerers |

The ownership rule is the target architecture: new *physics* belongs in a
residual-term lowerer, a new measurement in an observation/loss lowerer, and a
new search geometry in a metric. Until component lowering is introduced, a
feature also needs one explicitly bounded whole-target registration. It must
not become a solver branch keyed by PDE, control placement, or boundary
condition.

## The essential type distinction

At the mathematical level, a state coefficient vector and its derivative are different objects:

```text
primal coefficients:       delta x_h in X_h
covector coefficients:     xi_h      in X_h*
declared pairing:           <xi_h, delta x_h>
```

The default representation stores tested covector coefficients. If $`r_{j}=\langle E_{h},\psi_{j}\rangle`$, then

```math
\langle r,p\rangle_{Z_{h}^{\ast},Z_{h}}=r^{\mathsf T}p.
```

This is why a usual assembled matrix with test rows and trial columns can use the coordinate transpose in v0. It is **not** a general permission to add or remove mass matrices: a different dual representation needs a different pairing and transpose action.

```text
BlockLayout                          named spaces + dimensions
    │
    ├── PrimalBlockT<Backend>         variable / direction / adjoint-test seed
    └── CovectorBlockT<Backend>       residual / objective derivative / pullback

pair(covector, primal)               the only primitive dual pairing
```

See [`layout.hpp`](../../include/nmopt/contract/layout.hpp) and [`linalg.hpp`](../../include/nmopt/contract/linalg.hpp). Layout checks reject a state vector where a control vector is required, even when both use the same storage type and happen to have the same length.

## From a strong PDE to an input graph

A strong form is source material, not sufficient program input. For the current example,

```math
-\nabla\cdot(k\nabla y)+c y=f+u \quad\text{in }\Omega,
\qquad y=0 \quad\text{on }\Gamma_{D},
```

the program-level statement is the selected weak residual

```math
\langle E(y,u),v\rangle
=(k\nabla y,\nabla v)_{\Omega}+(c y,v)_{\Omega}
-(f,v)_{\Omega}-(u,v)_{\Omega}.
```

It requires choices that the displayed strong expression cannot determine:

```text
geometry/regions         Omega and Gamma_D
spaces                   state Y, test Z, control U, observation Q
essential condition      restriction or lifting for Gamma_D
residual terms           diffusion, reaction, forcing, volume control
observation              y restricted to an observation region
losses                   tracking norm and control regularisation
search geometry          a metric, if an algorithm needs a direction
constraints              admissible control set, if any
discrete policy          FEs, quadrature, coefficient representation, constraints
```

Essential and natural boundary conditions intentionally follow different paths:

```text
fixed Dirichlet data:    y_phys = P y_hat + ell_0      transformation/reconstruction
Dirichlet control:       y_phys = P y_hat + ell_0 + L_D u
Neumann data/control:    add a boundary residual term -<u, trace(v)>
```

Treating controlled Dirichlet data as a boundary load loses the chain rule through $`L_{D}`$; treating a Neumann control as a lifting changes a different problem. See the [boundary protocol](interface-specification.md#41-boundary-protocol).

## The executable boundary

After compilation, solvers only need this abstract surface:

```text
ExecutableModelT<Backend>
    residual(x)                 -> E_h(x) in Z_h*
    residual_jvp(x, dx)         -> E_h'(x) dx in Z_h*
    residual_vjp(x, p)          -> E_h'(x)* p in X_h*    (p is in Z_h)
    objective(x)                -> J_h(x)
    objective_derivative(x)     -> J_h'(x) in X_h*

MetricT<Backend>
    apply(g)                    -> G g
    inverse_apply(j_prime)      -> G^-1 j_prime

ConstraintT<Backend>
    is_feasible(u)
    project_in(u, metric)       -> metric-specific projection, when supported
```

The interface is [`executable_model.hpp`](../../include/nmopt/contract/executable_model.hpp); metric and constraint interfaces are in [`metric_constraint.hpp`](../../include/nmopt/contract/metric_constraint.hpp). No method accepts a PDE name or boundary-condition enum. That is what allows an optimizer to be reused once a lowerer has supplied these actions.

## DTO: how a control becomes a reduced covector

V0 uses **discretize then optimize**. It fixes $`E_{h}`$ and $`J_{h}`$ first, then derives discrete actions. Its global convention is

```math
\mathcal L_{h}(x_{h},p_{h})=J_{h}(x_{h})-
\langle p_{h},E_{h}(x_{h})\rangle.
```

For $`x_{h}=(y_{h},u_{h})`$, `ReducedDTOT::evaluate(control)` follows this exact path:

```text
control u_h
   │
   ├── solve_state(u_h) ─────────────> y_h satisfying E_h(y_h,u_h)=0
   │                                     │
   ├── compose(y_h, u_h) ───────────────> full point x_h
   │                                     │
   ├── objective_derivative(x_h) ───────> (J_y', J_u')
   │                                     │
   ├── solve_adjoint(x_h, J_y') ────────> p_h satisfying E_y'* p_h = J_y'
   │                                     │
   └── residual_vjp(x_h, p_h) ──────────> (E_y'*p_h, E_u'*p_h)
                                         │
reduced covector j_h'(u_h) = J_u' - E_u'*p_h
```

`j_h'` is still in $`U_{h}^{\ast}`$, not a control-space vector. A metric creates an algorithmic direction:

```math
g_{h}=G^{-1}j_{h}'.
```

Regularisation and metric must stay separate. Adding $`\frac{\alpha}{2}\lVert u\rVert_{H^{1}}^{2}`$ changes $`J_{u}'`$ and therefore the reduced derivative. Replacing the search metric by an $H^{1}$ metric changes only the last covector-to-direction map.

[`reduced_dto.hpp`](../../include/nmopt/contract/reduced_dto.hpp) implements this narrow workflow: one eliminated state block, one control block, one residual test block, and externally supplied state/adjoint solves.

## Direct v0 deal.II model: theory to assembled objects

`ScalarDiffusionReactionModel<dim>` realizes the scalar, stationary, homogeneous-Dirichlet case. Its discrete model is

```math
r_{h}(y_{h},u_{h})=A_{h}y_{h}-f_{h}-B_{h}u_{h},
```

```math
J_{h}(y_{h},u_{h})=
\frac{1}{2}y_{h}^{\mathsf T}M_{y}y_{h}-q_{y}^{\mathsf T}y_{h}
+\frac{1}{2}\lVert y_{d}\rVert_{L^{2}(\Omega)}^{2}
+\frac{\alpha}{2}u_{h}^{\mathsf T}M_{u}u_{h}.
```

| Theory | Discrete object | `ScalarDiffusionReactionModel` member/action |
| --- | --- | --- |
| Diffusion-reaction residual | $`A_{h}`$ | `system_matrix_` |
| Fixed source | $`f_{h}`$ | `forcing_load_` |
| Volume control coupling | $`B_{h}`$ | `control_coupling_` |
| State tracking | $`M_{y}`$, $`q_{y}`$, target constant | `state_mass_`, `desired_state_load_`, `desired_state_norm_` |
| Control $L^{2}$ regularisation | $`\alpha M_{u}`$ | `regularisation_weight_`, `control_mass_` |
| Fixed zero Dirichlet condition | constrained coordinate policy | `state_constraints_`, `constrained_state_dofs_` |
| Residual/JVP/VJP | $r$, $A\delta y-B\delta u$, $(A^{\mathsf T}p,-B^{\mathsf T}p)$ | `residual`, `residual_jvp`, `residual_vjp` |
| State/adjoint solve | $`A_{h}y=f+B u`$, $`A_{h}^{\mathsf T}p=J_{y}'`$ | `solve_state`, `solve_adjoint` |

The selected policy is scalar `FE_Q` state/test, `FE_DGQ(0)` control on the same active cells, constant $k>0$ and $c\geq0$, data sampled from deal.II `Function` objects at cell quadrature, serial assembled matrices, homogeneous Dirichlet ids, and SPD CG solves. `assemble()` populates cached matrices/vectors cell by cell; `residual*()` and `objective*()` compose those objects without rederiving the weak form.

### A sign check you can do by hand

The residual uses $A y-f-B u$, so its pullback with adjoint $p$ is

```math
E'(y,u)^{\ast}p=(A^{\mathsf T}p,-B^{\mathsf T}p).
```

The objective control derivative is $`J_{u}'=\alpha M_{u}u`$. The fixed DTO formula therefore gives

```math
j_{h}'(u)=\alpha M_{u}u-(-B^{\mathsf T}p)
=\alpha M_{u}u+B^{\mathsf T}p.
```

The plus sign is correct. When debugging a new term, write its residual sign, its VJP sign, and only then apply the single global subtraction in the formulation builder.

## Source-code tour in dependency order

| Read next | Why it exists | Read it when you need to understand or change |
| --- | --- | --- |
| [`linalg.hpp`](../../include/nmopt/contract/linalg.hpp) | Dense reference vector/matrix algebra and `ContractError` | Backend capability minimum and reference solves |
| [`layout.hpp`](../../include/nmopt/contract/layout.hpp) | Named layouts, primal/covector wrappers, pairing | Type and duality errors |
| [`executable_model.hpp`](../../include/nmopt/contract/executable_model.hpp) | Five universal executable actions | What a compiler/lowerer must emit |
| [`metric_constraint.hpp`](../../include/nmopt/contract/metric_constraint.hpp) | Metric and projection boundaries | Search direction versus feasibility |
| [`reduced_dto.hpp`](../../include/nmopt/contract/reduced_dto.hpp) | State/adjoint/reduced-covector orchestration | First-order workflow and signs |
| [`reduced_gradient.hpp`](../../include/nmopt/solvers/reduced_gradient.hpp) | Unconstrained and projected reduced Armijo method | First optimizer and its diagnostics |
| [`linear_quadratic_model.hpp`](../../include/nmopt/reference/linear_quadratic_model.hpp) | Transparent matrix oracle | Check algebra before FE assembly |
| [`serial_backend.hpp`](../../include/nmopt/dealii/serial_backend.hpp) | Adapter from five vector operations to deal.II | Backend parameterisation |
| [`mass_metric.hpp`](../../include/nmopt/dealii/mass_metric.hpp) | Sparse SPD decision Riesz map | deal.II $L^{2}$ and $H^{1}$ control or parameter search directions |
| [`cellwise_box_constraint.hpp`](../../include/nmopt/dealii/cellwise_box_constraint.hpp) | `FE_DGQ(0)` coefficientwise box projection | Feasible deal.II control or parameter updates |
| [`facewise_box_constraint.hpp`](../../include/nmopt/dealii/facewise_box_constraint.hpp) | Facewise-constant coefficientwise box projection | Feasible Neumann boundary-control updates |
| [`scalar_diffusion_reaction.hpp`](../../include/nmopt/dealii/scalar_diffusion_reaction.hpp) | Concrete deal.II lowerer | FE assembly, constraints, solves |
| [`types.hpp`](../../include/nmopt/semantic/v1/types.hpp) | Narrow deal.II-free v1 graph types | Semantic component ports |
| [`validation.hpp`](../../include/nmopt/semantic/v1/validation.hpp) | Structural and policy diagnostics | Semantic validation |
| [`compiled_problem.hpp`](../../include/nmopt/compiler/v1/compiled_problem.hpp) | Backend-generic compiled package and manifest | Solver-facing compiled ports and provenance |
| [`reduced_envelope.hpp`](../../include/nmopt/experiment/reduced_envelope.hpp) | In-memory manifest, policy, report, and environment association | Detached Chapter 6 experiment provenance |
| [`dealii_compiler.hpp`](../../include/nmopt/compiler/v1/dealii_compiler.hpp) | V1 registered deal.II compiler path | Capability checks plus private v0 comparison or v1 assembled targets |
| [`dealii_fixed_dirichlet.hpp`](../../include/nmopt/compiler/v1/dealii_fixed_dirichlet.hpp) | V1 physical-state assembly target | Independent coordinates, fixed lifting, material tracking, and pullbacks |
| [`dealii_dirichlet_control.hpp`](../../include/nmopt/compiler/v1/dealii_dirichlet_control.hpp) | V1 controlled physical-state target | Complete-boundary nodal lifting, trace metric, and state/control pullbacks |
| [`dealii_neumann_boundary.hpp`](../../include/nmopt/compiler/v1/dealii_neumann_boundary.hpp) | V1 natural-boundary assembly target | Facewise Neumann coupling, boundary tracking, pullbacks, and pure-Neumann mean-zero saddle solves |
| [`dealii_continuous_control.hpp`](../../include/nmopt/compiler/v1/dealii_continuous_control.hpp) | V1 continuous-control target | Continuous control loss and explicitly selectable search-metric realizations |
| [`dealii_coefficient_identification.hpp`](../../include/nmopt/compiler/v1/dealii_coefficient_identification.hpp) | V1 positive cellwise coefficient target | Parameter-dependent state/adjoint actions and exact first-order pullbacks |
| [`reduced_dto_contract.cc`](../../tests/reduced_dto_contract.cc) | Contract tests against dense oracle | Minimal executable example |
| [`semantic_v1_contract.cc`](../../tests/semantic_v1_contract.cc) | Deal.II-free graph validation | Structural and analytical-policy diagnostics |
| [`dealii_diffusion_contract.cc`](../../tests/dealii_diffusion_contract.cc) | Compiler and lowerer checks through real deal.II assembly | Exact scenarios are listed in the [v1 capability table](../implementation/v1/semantic-compiler.md#registered-capabilities) |

The dense model has no mesh or FE code. It makes an incorrect formula, type pairing, or DTO sign fail independently of deal.II. The deal.II test then establishes that the same contract survives real `DoFHandler`, quadrature, sparse assembly, `AffineConstraints`, and CG solve operations.

## Verification schematic

Every differentiable feature has four increasingly global checks:

```text
1. value test
   Does the assembled residual/objective match a simple expected value?

2. JVP Taylor test
   [E(x + eps dx) - E(x)] / eps  approximately equals  E'(x) dx

3. VJP pairing test
   <E'(x) dx, p>  approximately equals  <E'(x)* p, dx>

4. reduced Taylor test
   [j(u + eps du) - j(u)] / eps  approximately equals  <j'(u), du>
   where the perturbed state is solved again
```

The dense contract test exercises pairing, residual finite difference,
objective directional derivative, state residual, reduced derivative, metric,
box projection, and unconstrained/projected Armijo convergence. The semantic
test exercises graph and policy validation. The exact v1 deal.II scenarios are
listed only in the [capability table](../implementation/v1/semantic-compiler.md#registered-capabilities).

## Capability boundary

Do not mistake a documented architectural slot or a whitelisted semantic kind
for working functionality. The exact current registrations and unsupported
combinations are maintained in the
[v1 capability table](../implementation/v1/semantic-compiler.md#registered-capabilities)
and [v1 exclusions](../implementation/v1/semantic-compiler.md#exclusions).
The direct reference slice has its own
[v0 exclusions](../implementation/v0/dealii-lowerer.md#explicit-exclusions),
and the [roadmap](../planning/implementation-roadmap.md) alone records task
completion and handoff.

## Blueprint for adding one feature yourself

1. Classify the feature: region, space/pairing, variable/data, transformation, residual term, observation, loss, metric, constraint, or discretisation policy. If it fits none, record an architectural decision first.
2. State source/target spaces, regions, pairings, and requirements. A strong PDE label or norm name is not sufficient.
3. Write the tested value and JVP/VJP on paper, including the residual sign. A transformation must reach both residual and every observation using its physical field.
4. Select FE spaces, quadrature, dual-coefficient pairing, lifting/constraint/nullspace policy, and assembled or matrix-free execution.
5. Put code in the owning layer: PDE physics in residual lowering, tracking in observation/loss, Riesz map in metric, feasible-set operation in constraint.
6. Add value, JVP, and VJP pairing tests. If it affects an optimised variable, add a state-recomputed reduced Taylor test. Test the metric relation separately when a solver needs a direction.
7. Reject unsupported combinations with a capability diagnostic; never use a nearby but mathematically different fallback.

For worked deltas from the baseline, use the [Laplace growth case study](../case-studies/laplace-growth.md) and [formula-delta guide](../case-studies/laplace-interface-formulas.md). They show exactly which component changes for Neumann or Dirichlet control, observation changes, metrics, box constraints, coefficient identification, and time dependence.
