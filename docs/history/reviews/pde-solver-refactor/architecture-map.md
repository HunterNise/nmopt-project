# PDE–solver architecture map

## Purpose and status

This document is the architecture atlas produced by the R0 PDE–solver audit.
It records the implemented topology at the start of the refactor, the durable
numerical axes found in that implementation, and the correspondence between
the current code and the target
[PDE–solver boundary](../../../design/pde-solver-boundary.md).

The boundary document owns architectural rules and invariants. This map is
descriptive review evidence: it gives humans and agents a durable mental model
without requiring them to rediscover the system from large compiler and
numerical-target headers.

The [assessment](assessment.md) owns findings and their rationale. The
[deletion ledger](deletion-ledger.md) owns migration/deletion state. The
[roadmap](roadmap.md) owns execution status.

## How to read this document

Use the first part as a true matryoshka: each level magnifies one box from the
previous diagram.

```text
Level 0  whole system
    ↓ magnify the producer side
Level 1  producer paths
    ↓ magnify native compilation
Level 2  compiler / numerical middle
    ↓ magnify typed numerics
Level 3  typed-numerics decomposition
```

After Level 3, choose only the sibling zoom relevant to the task. The state,
decision, objective, ownership, and runtime sections are **different views of
Level 3**, not deeper universal layers.

| Task | Usually sufficient reading |
| --- | --- |
| Understand the whole architecture | Levels 0–3, current-to-target correspondence |
| Integrate an existing deal.II application | Levels 0–1, runtime zoom, ownership zoom |
| Change compiler/lowering | Levels 1–3, current-to-target correspondence |
| Change state representation | Level 3, state zoom |
| Change control/parameter realization | Level 3, decision zoom |
| Change observations/objective | Level 3, objective and geometry zoom |
| Change metric/constraint/Hessian | Level 3, objective and geometry zoom |
| Change B1/B2 output or diagnostics | Levels 1–3, ownership zoom |
| Locate current code after choosing a box | Code locator |

The boxes name responsibilities. They do **not** imply one C++ class per box.
A responsibility may be implemented by a value, callback bundle, free
function, focused concrete component, or an existing interface.

# Part I — Matryoshka magnification

## Level 0 — Whole-system architecture

At the coarsest level, two producers can feed the same generic mathematical
boundary.

```text
                    ┌──────────────────────────────────────────┐
                    │ PROBLEM / APPLICATION PRODUCER           │
                    │                                          │
                    │ semantic compiler path                   │
                    │          or                              │
                    │ existing deal.II application             │
                    └────────────────────┬─────────────────────┘
                                         │ expose
                                         ▼
                    ┌──────────────────────────────────────────┐
                    │ SOLVER PORTS                             │
                    │                                          │
                    │ E, E'v, E'^*p                            │
                    │ J, J'                                    │
                    │ state / adjoint solves                   │
                    │ metric + optional capabilities           │
                    └────────────────────┬─────────────────────┘
                                         │ orchestrate
                                         ▼
                    ┌──────────────────────────────────────────┐
                    │ FORMULATION                              │
                    │                                          │
                    │ Reduced DTO                              │
                    │ KKT products                             │
                    │ supplied OTD                             │
                    └────────────────────┬─────────────────────┘
                                         │ drive
                                         ▼
                    ┌──────────────────────────────────────────┐
                    │ ALGORITHM                                │
                    │                                          │
                    │ gradient / BFGS / L-BFGS / trust region  │
                    └────────────────────┬─────────────────────┘
                                         │ result / report
                                         ▼
                    ┌──────────────────────────────────────────┐
                    │ APPLICATION                              │
                    │                                          │
                    │ diagnostics / artifacts / native output  │
                    └──────────────────────────────────────────┘
```

The formulation and optimizer do not know which producer supplied the
mathematical operations.

The next level magnifies the producer box.

## Level 1 — Two producer paths

### Native semantic/compiler producer

```text
┌──────────────────────────────────────────┐
│ APPLICATION CONFIGURATION                │
│                                          │
│ CLI / parameter files                    │
│ run-set / benchmark registration         │
│ Scenario / Recipe                        │
└────────────────────┬─────────────────────┘
                     │ build
                     ▼
┌──────────────────────────────────────────┐
│ SEMANTIC                                 │
│                                          │
│ “What mathematical problem do I mean?”   │
│                                          │
│ ProblemSpec                              │
│ regions / variables / equations          │
│ observations / losses / policies         │
└────────────────────┬─────────────────────┘
                     │ compile / lower
                     ▼
┌──────────────────────────────────────────┐
│ TYPED NUMERICS                           │
│                                          │
│ “What concrete discrete thing exists?”   │
│                                          │
│ deal.II FE spaces                        │
│ DoFs / matrices / constraints            │
│ reconstruction / observations            │
│ concrete operators                       │
└──────────────┬─────────────────┬─────────┘
               │                 │
               │ adapt           │ retain
               ▼                 ▼
┌──────────────────────┐   ┌────────────────────────┐
│ SOLVER PORTS         │   │ NATIVE APP VIEW        │
│                      │   │                        │
│ ExecutableModelT     │   │ output                 │
│ solve callbacks      │   │ reconstruction         │
│ metric/capabilities  │   │ dimensions             │
└──────────┬───────────┘   │ diagnostics/evidence   │
           │               └───────────┬────────────┘
           ▼                           │
     formulation                       │
           │                           │
           ▼                           │
       algorithm                       │
           │                           │
           └────── result/report ──────┤
                                       ▼
                                 application /
                                  experiment
```

The target rule is visible here:

> Adapt typed numerics to generic solver ports **while retaining the typed
> native view beside them**.

### Existing deal.II application producer

An external application skips the semantic/compiler stack:

```text
┌──────────────────────────────────────────┐
│ EXISTING deal.II APPLICATION             │
│                                          │
│ mesh / FE / DoFs                         │
│ assembly / native solvers                │
│ original output                          │
└──────────────┬─────────────────┬─────────┘
               │                 │
               │ adapt           │ retain
               ▼                 ▼
┌──────────────────────┐   ┌────────────────────────┐
│ SOLVER PORTS         │   │ ORIGINAL APP VIEW      │
│                      │   │                        │
│ residual callbacks   │   │ standalone solve       │
│ objective callbacks  │   │ native output          │
│ state/adjoint solves │   │ diagnostics            │
│ metric/capabilities  │   │ application state      │
└──────────┬───────────┘   └────────────────────────┘
           │
           ▼
      formulation
           │
           ▼
       algorithm
```

The original PDE class need not inherit an nmopt PDE interface and need not
adopt `ProblemSpec`, compiler manifests, recipes, scenarios, or the runner.

The next level magnifies **native compilation** from the first producer path.

## Level 2 — Compiler and numerical middle

```text
ProblemSpec
    │
    ▼
┌─────────────────────────────────────────────┐
│ SEMANTIC RESOLUTION                         │
│                                             │
│ validate graph                              │
│ resolve IDs                                 │
│ expose stable read-only semantic view       │
└───────────────────┬─────────────────────────┘
                    ▼
             ResolvedProblemView
                    │
                    ▼
┌─────────────────────────────────────────────┐
│ CLOSED COMPILATION DECISION                 │
│                                             │
│ supported component combination             │
│ backend-neutral realization choices         │
│ required policies / regions / selections    │
└───────────────────┬─────────────────────────┘
                    │
                    ├── deal.II runtime data
                    ├── mesh / session
                    └── discretization policy
                    ▼
┌─────────────────────────────────────────────┐
│ BACKEND LOWERABILITY                        │
│                                             │
│ reference cells / mesh facts                │
│ boundary + material IDs                     │
│ runtime bindings                            │
│ sensor/value compatibility                  │
└───────────────────┬─────────────────────────┘
                    ▼
┌─────────────────────────────────────────────┐
│ NUMERICAL REALIZATION                       │
│                                             │
│ FE / DoFs / constraints                     │
│ coordinates / reconstruction                │
│ state + decision operators                  │
│ observations / objective ingredients        │
│ solve services / metrics / capabilities     │
└───────────────────┬─────────────────────────┘
                    ▼
          typed numerical realization
               ╱                 ╲
              ╱                   ╲
             ▼                     ▼
    generic solver view       native/app view
             │                     │
             ▼                     └──► native data / evidence
      CompiledProblemT
             │
             └──► formulation / generic products
```

### Current compiler pressure

The current compiler re-encodes substantially the same whole-target choice in
several forms:

```text
semantic component graph
        ↓
request booleans / uses_* predicates
        ↓
ResolvedTargetFamily
        ↓
private compiled-target identity
        ↓
whole-target construction branch
        ↓
manifest / compatibility projection
```

The target is not necessarily one giant enum. The target is **one authoritative
closed realization decision** whose structured choices are consumed by
construction and provenance rather than reconstructed repeatedly.

One scalar path already shows the useful upper half of this target:

```text
semantic residual components
        ↓
registered scalar lowerers
        ↓
ScalarLoweringPlan
        │
        ├── diffusion / transport / reaction
        ├── source / volume control
        └── observations / losses
        ↓
current complete-registration match
        ↓
current ScalarComponentModel
```

The refactor should preserve component-aware lowering while avoiding an
unnecessary return to complete-problem identity.

The next level magnifies the `NUMERICAL REALIZATION` box.

## Level 3 — Typed-numerics decomposition

The audit found a numerical middle layer with several independent axes, but no
need for one universal `NumericalModel`.

```text
                            ┌─────────────────────┐
                            │ STATE SPACE         │
                            │ FE / DoFs / mesh    │
                            └──────────┬──────────┘
                                       │
                                       ▼
                            ┌─────────────────────┐
                            │ STATE COORDINATES   │
                            │ reconstruction      │
                            │ pullback / gauge    │
                            └──────────┬──────────┘
                                       │
                                       ▼
                              physical state field
                                       │
                ┌──────────────────────┼──────────────────────┐
                │                      │                      │
                ▼                      ▼                      ▼
       ┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐
       │ STATE OPERATOR  │    │ OBSERVATION     │    │ NATIVE VIEW     │
       │                 │    │                 │    │                 │
       │ PDE residual    │    │ physical-state  │    │ output          │
       │ transpose facts │    │ measurement     │    │ diagnostics     │
       └────────┬────────┘    └────────┬────────┘    └─────────────────┘
                │                      │
                │                      ▼
                │             tracking / loss data
                │
                │             ┌───────────────────────────────┐
                │             │ DECISION DISCRETIZATION       │
                │             │                               │
                │             │ layout / coordinates          │
                │             │ M / K / trace ingredients     │
                │             └──────────────┬────────────────┘
                │                            │
                │             ┌──────────────┼──────────────┐
                │             │              │              │
                │             ▼              ▼              ▼
                │       PDE influence   regularization   metric/constraint
                │             │
                └─────────────┼─────────────────────────────────────┐
                              │                                     │
                              ▼                                     ▼
                       executable ops                      solve services /
                         E, E'v, E'^*p                      optional Hessian
                         J, J'
```

| Axis | Demonstrated current variants |
| --- | --- |
| State coordinates | full FE coordinates; independent coordinates with $P$ and lifting; Dirichlet-control reconstruction; mean-zero gauge |
| State operator | fixed symmetric; fixed nonsymmetric; decision-dependent $A(m)$ |
| Decision coordinates | DG0 volume; continuous volume; facewise Neumann; continuous trace; Dirichlet trace; DG0 parameter |
| Decision → PDE influence | affine volume/boundary coupling; state lifting; parameterized operator |
| Observation | full/subdomain volume; $H^{1}$ state; boundary trace; weighted trace; point sensors; normal flux |
| Regularization | $L^{2}$, $H^{1}$, trace $H^{1/2}$ |
| Search metric | $L^{2}$, $H^{1}$, $H^{-1}$, trace $H^{1/2}$ |
| Constraint capability | cellwise box; facewise box; unsupported/optional elsewhere |
| Solve structure | symmetric SPD; nonsymmetric direct; augmented mean-zero; point-dependent SPD |

These are **axes**, not instructions to introduce matching classes.

# Part II — Sibling zooms from typed numerics

## State realization zoom

State representation, physical operator, gauge, and solve realization vary
independently.

### State coordinates

```text
FULL COORDINATES
────────────────────────────────

solver state y
    =
physical FE coefficients


AFFINE INDEPENDENT COORDINATES
────────────────────────────────

independent state y_hat
        │
        ▼
y_phys = P y_hat + ell

tangent:  delta y_hat ───► P delta y_hat
pullback: q_phys ─────────► P^T q_phys


```

`ScalarComponentModel` demonstrates the affine independent-state mechanism.
Dirichlet control extends this coordinate map with a decision-dependent lifting
term; that cross-cutting case is shown once in the
[decision realization zoom](#decision-realization-zoom). These are
coordinate/reconstruction mechanisms, not PDE-family identities.

### Physical state operator and solve facts

```text
FIXED SYMMETRIC
A y
    ↓
SPD state / adjoint solves


FIXED NONSYMMETRIC
A y, A^T p
    ↓
direct or other forward / transpose solves


DECISION-DEPENDENT
A(m) y
D_m A(m)[delta m] y
parameter transpose derivative
    ↓
assemble-at-point + solve
```

Pure-Neumann/mean-zero treatment adds an orthogonal gauge/nullspace decision.

Architectural consequence:

> Do not make the typed state layer assume one persistent sparse matrix or one
> linear solver. The universal downstream language remains $E$, $E'v$,
> $E'^*p$, plus supplied state/adjoint solve callbacks.

## Decision realization zoom

The audit found three different mathematical mechanisms by which a decision
affects the PDE.

```text
1. AFFINE FORCING
────────────────────────────────

u ───────────► B u ───────────► residual

reverse:
p ───────────► B^T p

Examples:
DG0 volume control
continuous volume control
Neumann boundary control


2. STATE RECONSTRUCTION
────────────────────────────────

y_hat ─► P y_hat ─────┐
                      ├──► y_phys ───► residual + observation
u ─────► L_D u ───────┤
ell_0 ────────────────┘

reverse:
q_phys ─► P^T q_phys
       └► L_D^T q_phys

Example:
Dirichlet control


3. OPERATOR PARAMETER
────────────────────────────────

m ───────────► A(m)
               │
               ├──► A(m)y
               ├──► D_m A(m)[delta m] y
               └──► parameter transpose derivative

Example:
coefficient identification
```

Decision coordinates, PDE influence, and optimization geometry are separate:

```text
DecisionDiscretization
        │
        ├── optimizer layout / physical embedding
        ├── reusable M / K / trace ingredients
        └── geometry / coordinates
        │
        ├────────► PDE influence
        ├────────► regularization
        ├────────► search metric
        └────────► optional constraint
```

A broad runtime `DecisionRealisation` interface would conflate orthogonal
concerns and duplicate `ExecutableModelT`. A compiler-only closed realization
variant may still be useful because the compiler must know which mechanism to
construct.

## Observation, objective, and optimization-geometry zoom

Observations act naturally on the **physical state**:

```text
solver state coordinates
        │
        ▼
state reconstruction
        │
        ▼
physical FE state
        │
        ├──► volume / subdomain observation
        ├──► H1 observation
        ├──► boundary / weighted trace
        ├──► point sensors
        └──► normal flux
```

Two numerical representations are useful for different consumers.

```text
EXPLICIT MEASUREMENT
────────────────────────────────

physical state y
      │
      ▼
     O y
      │
      ├──► diagnostics / evidence / output
      └──► transpose action


COMPILED QUADRATIC TRACKING
────────────────────────────────

observation + target + pairing
      │
      ▼
Q = O^* W O
q = O^* W d
c = d^* W d
      │
      ▼
J_track(y) = 1/2 y^T Q y − q^T y + 1/2 c
```

Point sensors, trace, and flux may need the explicit observation itself.
Volume and $H^{1}$ tracking can often use the compiled quadratic form directly.

### Regularization versus metric

```text
decision numerical ingredients
M / K / trace / quotient operators
        │
        ├────────────────────────────┐
        │                            │
        ▼                            ▼
OBJECTIVE REGULARIZATION         SEARCH METRIC
R                                G, G^{-1}
changes J                        changes optimizer geometry
        │                            │
        ▼                            ▼
      J, J'                    search direction
```

One matrix may realize both $R = M$ and $G = M$. Shared numerical data does not
merge semantic ownership.

Demonstrated metric constructions include:

```text
L2:         G = M
H1:         G ≈ M + K
H-1:        G = M K^{-1} M
trace H1/2: G = A_BB − A_BI A_II^{-1} A_IB
```

Solves intrinsic to negative/fractional metric application are part of applying
the mathematical metric. They are not a reason to introduce a universal PDE
linear-solver interface.

### Optional reduced Hessian

For the fixed linear-quadratic affine-control case:

```text
control direction w
        ↓
      B w
        ↓
solve A delta y = B w
        ↓
      Q delta y
        ↓
solve A^T delta p = Q delta y
        ↓
B^T delta p + alpha R w
        ↓
     H_red w
```

`ReducedHessianT` is therefore an optional composed second-order capability,
not a primitive required of every PDE or every first-order integration.

## Runtime formulation zoom

The first-order reduced runtime is already a clean boundary:

```text
control u
    │
    ▼
evaluate_value(u)
    │
    ├──► solve_state(u) ───► state y
    ├──► compose x = (y, u)
    └──► objective J(x)

derivative needed
    │
    ▼
augment_derivative(...)
    │
    ├──► objective_derivative(x)
    ├──► solve_adjoint(x, J'_y)
    ├──► residual_vjp(x, p)
    └──► j'(u) = J'_u − E'_u(x)^* p
                 │
                 ▼
          reduced covector
                 │
                 ▼
       MetricT::inverse_apply
                 │
                 ▼
          search direction
                 │
                 ▼
       globalization / step
```

`ReducedDTOT` means **Discretize Then Optimize**. The formulation owns the call
order. Callback providers own the concrete state/adjoint solves. The algorithm
owns search, globalization, and convergence.

The value/derivative split permits rejected trial points to avoid an adjoint
solve when the selected algorithm does not require one.

## Ownership and lifetime zoom

```text
┌──────────────────────────────────────────────┐
│ APPLICATION / COMPILATION SESSION            │
│                                              │
│ mesh lifetime                                │
│ runtime data                                 │
│ execution / parameters                       │
└───────────────────┬──────────────────────────┘
                    │ owns / outlives
                    ▼
┌──────────────────────────────────────────────┐
│ TYPED NUMERICAL REALIZATION                  │
│                                              │
│ DoF handlers                                 │
│ constraints / coordinate maps                │
│ assembled or point-dependent operators       │
│ observations / numerical ingredients         │
│ reconstruction / native diagnostics          │
└───────────────┬───────────────────┬──────────┘
                │                   │
                │ adapt             │ retain
                ▼                   ▼
┌────────────────────────┐   ┌──────────────────────┐
│ SOLVER-FACING OBJECTS  │   │ NATIVE APP VIEW      │
│                        │   │                      │
│ executable facade      │   │ output               │
│ solve callbacks        │   │ dimensions           │
│ metric/capabilities    │   │ diagnostics          │
└────────────┬───────────┘   └──────────────────────┘
             │
             ▼
        formulation
             │
             ▼
          optimizer
```

Useful current lifetime mechanisms should survive:

- compilation sessions may own triangulations;
- DoF handlers borrow meshes;
- compiler packages retain executable/capability objects;
- callbacks may capture shared typed numerical state;
- `ReducedDTOT` can retain an opaque lifetime owner when detached lifetime is
  required.

The refactor changes **responsibility ownership**, not the fundamental lifetime
facts. It should not copy meshes, DoF handlers, or large operators merely to
hide borrowing.

# Part III — Current-to-target correspondence

## Current structural pressures

| Pressure | Current pattern | Target direction |
| --- | --- | --- |
| Composition collapse | semantic components → component-aware plan → complete target `*Model` | retain meaningful numerical choices farther into lowering |
| Early type erasure | typed target → `ExecutableModelT` → application `dynamic_cast` | erase solver operations only; retain native view separately |
| Role conflation | coupling + regularization + metric + constraint + output in one target | share numerical ingredients while keeping role ownership distinct |
| Repeated compiler authority | `uses_*` → target family → target kind → branch → manifest copy | one structured closed decision |
| Legacy duplication | direct v0 scalar implementation beside newer compiler path | transfer unique oracle/capability, then delete redundant path |

These are refactor pressures, not requested class hierarchies.

## Current-to-target type correspondence

| Current structure | Useful part | Target correspondence |
| --- | --- | --- |
| `ExecutableModelT` | small generic mathematical facade | survives |
| `StateAdjointSolversT` | inversion-of-control solve seam | survives |
| `ReducedDTOT` | first-order reduced orchestration | survives |
| `MetricT`, `ConstraintT`, `ReducedHessianT` | generic capabilities | survive |
| `CompiledProblemT` | compiler-path packaging and lifetime retention | remains compiler-specific; may simplify |
| `ScalarLoweringPlan` | component-aware semantic lowering | evolves toward one closed numerical decision |
| `VolumeObservationAssembly` | focused volume tracking realization | reuse/evolve where it replaces embedded duplication |
| `NeumannControlRealisation` | demonstrated trace-control runtime polymorphism | narrow to demonstrated discretization/coupling responsibility |
| complete compiler `*Model` classes | tested numerical machinery | migrate useful pieces; complete-target identity is not the target architecture |
| model-owned output | working native visualization | move filesystem/output orchestration to application side; retain typed field access |
| B1/B2 downcasts | access to typed data | delete by retaining native view explicitly |
| parallel compiler identities | dispatch evidence | collapse into one authoritative decision |
| manifest compatibility copies | persisted compatibility evidence | project once at serialization edge |
| direct v0 scalar model | hand-written oracle | transfer unique behavior, then likely retire duplicate path |

# Part IV — Code locator

Use this only after selecting the relevant architectural box. These are
representative current locations, not permanent ownership promises.

| Responsibility | Representative current location |
| --- | --- |
| Layout/primal/covector contracts | `include/nmopt/contract/layout.hpp`, `include/nmopt/contract/linalg.hpp` |
| Executable facade | `include/nmopt/contract/executable_model.hpp` |
| Reduced DTO and solve callbacks | `include/nmopt/contract/reduced_dto.hpp` |
| Metrics/constraints contract | `include/nmopt/contract/metric_constraint.hpp` |
| Optional reduced Hessian | `include/nmopt/contract/reduced_hessian.hpp` |
| Generic KKT products | `include/nmopt/contract/quadratic_kkt.hpp` and related contract headers |
| Reduced algorithms | `include/nmopt/solvers/` |
| Backend metrics | `include/nmopt/dealii/mass_metric.hpp`, `hminus1_metric.hpp`, `trace_hhalf_metric.hpp` |
| Backend constraints | `include/nmopt/dealii/cellwise_box_constraint.hpp`, `facewise_box_constraint.hpp` |
| Backend solve helper | `include/nmopt/dealii/serial_spd_solver.hpp` |
| Direct v0 scalar model | `include/nmopt/dealii/scalar_diffusion_reaction.hpp` |
| Semantic graph | `include/nmopt/semantic/v1/types.hpp` |
| Semantic validation/resolution | `include/nmopt/semantic/v1/validation.hpp`, `resolved_problem.hpp` |
| Reference semantic builders | `include/nmopt/semantic/v1/reference_specs.hpp` |
| Compiler composition root | `include/nmopt/compiler/v1/dealii_compiler.hpp` |
| Scalar component plan | `include/nmopt/compiler/v1/dealii_scalar_plan.hpp` |
| General scalar realization | `include/nmopt/compiler/v1/dealii_fixed_dirichlet.hpp` |
| Continuous volume control | `include/nmopt/compiler/v1/dealii_continuous_control.hpp` |
| Neumann boundary target | `include/nmopt/compiler/v1/dealii_neumann_boundary.hpp` |
| Neumann control realizations | `include/nmopt/compiler/v1/dealii_neumann_control_realisation.hpp` |
| Dirichlet control lifting | `include/nmopt/compiler/v1/dealii_dirichlet_control.hpp` |
| Coefficient identification | `include/nmopt/compiler/v1/dealii_coefficient_identification.hpp` |
| Volume observation assembly | `include/nmopt/compiler/v1/dealii_volume_observation.hpp` |
| Compiler package | `include/nmopt/compiler/v1/compiled_problem.hpp` |
| Chapter 5/6 options | `include/nmopt/application/chapter5.hpp`, `chapter6.hpp` |
| B1 native composition | `include/nmopt/application/dealii/chapter6_b1.hpp` |
| B2 native composition | `include/nmopt/application/dealii/chapter6_b2.hpp` |
| Experiment envelope | `include/nmopt/experiment/` |
| CLI / run-set composition | `apps/nmopt-runner/` |
| Artifact / post-processing consumers | `tools/` |

## Architectural reading summary

The strong generic core already exists. The principal refactor need is in the
middle:

```text
semantic meaning
      ↓
closed numerical choices
      ↓
typed numerical pieces
      ├──► narrow generic mathematical ports
      └──► retained native/application access
```

When an implementation choice is unclear, descend through Levels 0–3, choose
the relevant sibling zoom, inspect the current consumers in the code locator,
and prefer the smallest representation that preserves the boundary while
enabling demonstrated deletion or reuse.
