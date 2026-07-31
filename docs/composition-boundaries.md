# Composition boundaries: an actionable architecture summary

## Decision in one sentence

Do not model a “maximally general PDE-control problem.”  Model a typed graph
of small components that contribute residuals, transformations, observations,
metrics, constraints, and discretisation choices.  The only shared language is
spaces, pairings, regions, and derivative actions.

The problem shell should be a container of connected components, not a base
class with a growing list of optional PDE, boundary-condition, control, and
solver fields.

```text
 spaces + pairings + named regions
        |          |             |
 variable/data -- transformations -- residual equation blocks
                                      |
                                  state equations

 variable/data -- observations -- losses/objective
        |
   admissible-set constraints

 all semantic components -- discretisation policy --> executable operators
```

## Minimal component set

| Component | It owns | Its ports | It must not own |
|---|---|---|---|
| `Region` | domain, subdomain, boundary, time interval identity | geometry and named-set references | PDE or FE policy |
| `Space` / `Pairing` | field shape, continuous trial/test space, dual pairing, trace/product requirements | used by every map | a particular finite element |
| `Variable` | unknown block, role, continuous space, admissible-set reference | residual, objective, transformation | state/adjoint solver logic |
| `Data` | immutable coefficient, source, target, boundary/initial datum | terms, observations, transformations | derivative or optimisation behaviour |
| `Transformation` | a map such as lifting, trace, restriction, parameterisation, or actuator | value, Jacobian, transpose/chain rule | an equation or an objective by itself |
| `ResidualTerm` | one tested contribution to an equation | value, JVP, transpose-JVP | global adjoint orchestration or optimiser choices |
| `EquationBlock` | sum of terms and its test-space dual target | residual and linearised actions | a PDE family name |
| `Observation` | map from physical variables to an observation space | value, JVP, transpose-JVP | loss or PDE residual |
| `Loss` | scalar functional on an observation/control/parameter space | value and derivative | the observation map or search metric |
| `Metric` | a declared primal-dual identification and, if available, inverse | pairing, Riesz/duality actions | the PDE or objective definition |
| `Constraint` | admissible-set, projection, normal-cone, active-set, or multiplier relation | variable block and solver/formulation | residual terms |
| `DiscretisationPolicy` | finite-element trial/test realizations, quadrature, constraints, lifting realization, assembled/matrix-free choice | lowering of semantic ports | continuous meaning or well-posedness |

The public `ProblemSpec` is just the composition root: it registers these
objects, connects their named ports, and validates that every referenced port
exists.  It should contain no `if (pde_type == ...)` or `if (solver == ...)`
branches.

## The universal executable ports

Every nonlinear semantic map needed by an algorithm reduces to three actions:

$$
  F(x),\qquad F'(x)\delta x,\qquad F'(x)^{\ast}q.
$$

For a residual, $F(x)$ lies in a declared test-space dual.  For an observation
or transformation, its source and target are declared directly.  For a loss,
the derivative lies in the dual of its argument space.  This is enough to
compose the state residual, adjoint right-hand sides, reduced derivatives, and
KKT blocks without adding PDE-specific solver code.

## What changes when one component changes

The table is deliberately directional.  “Affects” means that the component
must be recomposed or recompiled; it does not mean that it takes ownership of
the affected concern.

| Change one thing | Local owner | Necessarily affected | Deliberately unchanged |
|---|---|---|---|
| Change a volume source $f$ | `Data` | residual values | spaces, adjoints’ structure, optimiser |
| Add diffusion/reaction/transport | `ResidualTerm` | residual and its derivative/transpose | objective, constraints, outer optimiser |
| Replace volume control by Neumann control | control `Variable` + residual coupling term | control space, boundary trace requirement, reduced derivative | equation/adjoint workflow, optimiser |
| Change the target region for tracking | `Observation` | loss derivative and therefore adjoint forcing | state residual, control constraint, optimiser |
| Change a tracking norm | `Loss` and its observation-space pairing | objective derivative and adjoint forcing | state residual |
| Change only search geometry from $L^2$ to $H^1$ or $H^{-1}$ | algorithmic `Metric` | dual-to-search-direction conversion | objective, residual, adjoint equation |
| Add a box constraint | `Constraint` | optimality condition and constrained optimisation method | PDE terms, observation, adjoint equation |
| Change scalar to vector field | `Space` field shape plus compatible terms | tensor contractions and FE component layout | residual/adjoint protocol, optimiser |
| Add a mixed flux variable | new `Variable` and `EquationBlock` terms | block residual and block adjoint | objective/constraint architecture |
| Switch Galerkin to Petrov–Galerkin | `DiscretisationPolicy` | discrete trial/test spaces and discrete transpose | continuous semantic problem |
| Add stabilization/upwinding | residual term and sometimes test-space policy | discrete/continuous residual and transpose | objective and outer optimiser |
| Replace a coefficient datum by an estimated coefficient | promote `Data` to `Variable` + parameter residual term | residual nonlinearity, parameter derivative, parameter constraints | adjoint orchestration, optimiser protocol |

Two distinctions prevent common accidental coupling:

1. An **objective norm** is part of a `Loss`, hence changing it changes the
   adjoint source.  An **algorithmic gradient metric** maps a dual derivative
   to a primal search direction, hence changing it does not change the
   objective or adjoint equation.
2. A **natural boundary contribution** is a residual term.  An **essential
   boundary condition** changes a state-space parameterisation and is a
   transformation/constraint.  They must not share a generic “boundary load”
   interface.

## Cross-cutting aspects that cannot be fully isolated

Some interactions are mathematical, not architectural accidents.  The design
must make them explicit rather than hide them in a component.

| Cross-cutting aspect | Why it crosses boundaries | Required explicit contract |
|---|---|---|
| Spaces and dual pairings | Every residual, observation, trace, metric, and derivative is typed by them | source, target, pairing, regularity/trace/product requirement |
| Regions and geometry | Terms, controls, observations, and boundary conditions all refer to physical sets | named-region identity and dimension/measure kind |
| Essential boundary conditions | They alter the physical state before both residual and observation evaluation | homogeneous space plus lifting/reconstruction and chain rule |
| Dirichlet boundary control | It simultaneously determines control trace space, lifting, residual dependence, and reduced derivative | $U \rightarrow \mathrm{state}$ lifting with derivative and discrete realization |
| Time dependence | It changes spaces, residual endpoint terms, observations, controls, and the discrete transpose | time-space declaration plus initial/terminal policy |
| Nullspaces and compatibility | A residual’s kernel affects solvability, adjoints, metrics, and linear solvers | explicit gauge/mean/compatibility/nullspace policy |
| Very-weak formulations | Moving derivatives changes state/test spaces, data pairings, and admissible observations | full residual formulation; never an automatic flag |
| Discrete-only objects | Point sensors, fractional norms, and some traces may lack the desired continuous map | an explicit discrete-only policy and lowering |
| Nonconforming/stabilized schemes | The actual discrete residual, not just its spaces, determines the correct discrete adjoint | compiled residual and exact transpose action |

These are the cases where a validator should require a policy or a supplied
map.  It may not silently choose one.

## Boundary conditions: ownership rules

| Boundary statement | Semantic owner | Additional components touched |
|---|---|---|
| Fixed Dirichlet datum | lifting/transformation of state | state space, compiler’s affine constraints |
| Dirichlet control | control-to-state lifting | control space, residual, objective through physical state, chain rule |
| Neumann datum/control | boundary residual term | trace pairing; control variable if optimised |
| Robin datum | boundary bilinear residual term plus boundary functional | coefficient/data and region |
| Periodic/hanging relation | discrete constraint realization | matching semantic region/identification metadata |
| Pure Neumann | residual terms plus nullspace policy | constraints, metric, linear solver/preconditioner |
| Transport inflow/outflow | selected transport residual formulation | region orientation, trace and stabilization requirements |

## Time: a structured extension, not a separate problem hierarchy

An evolution equation should not lead to `ParabolicProblem` and
`HyperbolicProblem` base classes.  It adds the following components to an
otherwise ordinary residual graph:

- time-indexed spaces for state, test, control, and observations;
- one or more temporal residual terms ($\dot y$, or first-order blocks
  $\dot y-v$ and $\dot v+A y$);
- initial-trace data or a corresponding equation block;
- terminal loss and/or terminal test policy when present; and
- a temporal discretisation policy, whose compiled transpose determines the
  discrete backward adjoint.

The spatial PDE terms, control couplings, observations, losses, constraints,
and metric interfaces remain the same components.

## The non-negotiable composition rule for inputs

The program should accept a selected variational residual, not infer one from
a strong differential expression.  A component may expose a strong-form
label for documentation, but its executable semantic contract is its tested
action and derivative ports.

For example, the textbook equation

$$
  -\Delta y=f+u,\qquad y|_\Gamma=0
$$

becomes the small graph

```text
state variable:       y in V = H^1_0(Omega)
test space:           Z = V
residual term:        (grad y, grad v) - (f, v) - (u, v)
control coupling:     U=L^2(Omega) -> V*  [volume source]
observation + loss:   y -> L^2(omega_o) -> tracking scalar
control loss:         u -> L^2(Omega) -> regularisation scalar
metric (optional):    U <-> U* for a search direction
constraint (optional): box constraint on U
```

Changing $u$ to Neumann control replaces only the control-coupling edge by a
boundary trace-adjoint edge.  Changing it to Dirichlet control instead adds a
control-to-state lifting before residual and observation evaluation; this is
one of the intentional cross-cutting cases above.

## Decisions to make now

1. Fix the typed ports: spaces, pairings, map value/JVP/transpose-JVP,
   transformations, residual blocks, observations/losses, metrics, and
   constraints.
2. Fix one global Lagrangian sign and require every compiled operator to
   expose its exact transpose action.
3. Make regions, field shapes, and trial/test spaces first-class references;
   never use PDE-specific string switches as a substitute.
4. Treat liftings/reconstructions and discrete lowering as first-class
   components from day one, even if the first implementation supports only
   simple homogeneous Dirichlet data.
5. Add a requirements/policy mechanism for the cross-cutting cases rather
   than attempting automatic analysis.
6. Keep the first executable slice narrow.  Its interfaces, not its feature
   list, need to be general.

## Fast architecture review checklist

When proposing a feature, ask:

1. Is it a residual term, transformation, observation, loss, metric,
   constraint, or discretisation policy?  If none, is a new primitive truly
   required?
2. What are its source and target spaces, pairings, and regions?
3. Does it supply a value, JVP, and transpose-JVP action?
4. Does it modify the physical state through a transformation, especially at
   an essential boundary?
5. Does it require an explicit cross-cutting policy (trace, nullspace,
   endpoint, point sensor, fractional norm, or discrete-only meaning)?
6. Can an optimiser consume it solely through the generic executable ports?

If the final answer is no because a solver needs to know the PDE name, the
proposed boundary is too coupled.
