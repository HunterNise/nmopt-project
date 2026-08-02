# V1 semantic graph and deal.II compiler

## Status and boundary

This is the first public semantic-to-compiler path. It is deliberately named
**v1** and exists alongside the direct v0 model:

```text
v0 direct reference
ScalarDiffusionReactionModel<dim>

v1 semantic/compiler path
semantic::v1::ProblemSpec
  -> SemanticValidator + compiler diagnostics
  -> compiler::v1::CompiledProblem<dim>
  -> separately constructed ScalarDiffusionReactionModel<dim>
```

The v0 lowerer is not modified or replaced by this path. The v1 compiler owns
a separate executable instance produced from a validated `ProblemSpec`; v0
and v1 can therefore be evaluated side by side at the same coefficients.
This first registration intentionally reuses the already-tested v0 assembly
as its executable target. A later independently assembled v1 lowerer may be
registered without changing the semantic API or the v0 reference.

The normative semantic protocol remains the
[interface specification](interface-specification.md). This document records
the concrete, intentionally narrow v1 realization of that protocol.

## Public semantic graph

`include/nmopt/semantic/v1/problem_spec.hpp` contains deal.II-free component
descriptions and `make_scalar_diffusion_reaction_problem()`. The current
factory declares only this graph:

```text
Region       one full volume region and homogeneous Dirichlet boundary ids
Space        scalar H1 state/test; scalar L2 control and observations
Pairing      explicit coefficient pairings for state, test, control, and observations
Variable     one state and one control
Data         forcing, desired state, diffusion, reaction, regularisation,
             and optional lower/upper cellwise bounds
Residual     diffusion-reaction, volume source, and volume control
Observation  full-volume restriction of state and control
Loss         quadratic tracking and quadratic control regularisation
Metric       cellwise L2 control metric
Constraint   optional cellwise L2 box
Formulation  one-state/one-control reduced DTO
```

Concrete values are not semantic objects. `DealiiDataBindings<dim>` binds the
forcing and desired-state `Function` objects plus scalar coefficients only
after validation. If the graph declares its optional box, the compiler also
requires `CellwiseBoxDataBindings`: both bounds must be scalar constants or
both must be exact-layout `FE_DGQ(0)` coefficient vectors.

## Validation and diagnostics

`SemanticValidator` validates semantic structure and declared policies.
`compiler::v1::DealiiCompiler::validate()` appends compiler-specific checks
to the same `ValidationReport`.

| Category | Produced by | Examples in v1 |
| --- | --- | --- |
| `structural` | `SemanticValidator` | missing ports, absent equation test space, wrong term inputs |
| `analytical_policy` | `SemanticValidator` | missing selected fixed-Dirichlet or cellwise-bound policy |
| `lowerability` | `DealiiCompiler` and its `DealiiLowererRegistryV1` | matrix-free execution, zero `FE_Q` degree, unregistered node kind, missing bound binding |
| `formulation_capability` | `DealiiCompiler` | all-at-once formulation or a multi-block DTO request |

`CompilationResult<dim>` returns the report and only contains a
`CompiledProblem<dim>` when it is valid. Unsupported choices are reported as
diagnostics; the compiler does not substitute another residual, observation,
metric, constraint, or formulation.

## Registered deal.II realization

`include/nmopt/compiler/v1/dealii_scalar_diffusion_reaction.hpp` registers
only the listed volume terms, full-domain volume restriction, quadratic
losses, `L2` metric, and optional cellwise box. Its selected discrete policy
is the v0 assembled serial realization: scalar `FE_Q` state/test with degree
at least one, `FE_DGQ(0)` control on the same mesh, homogeneous Dirichlet
boundary ids, and reduced DTO.

`CompiledProblem<dim>::executable_model()` exposes the separately owned
executable for formulation construction. `control_l2_metric()` and
`control_constraint()` expose only the metric and optional constraint that
were declared by the graph. Solvers continue to consume the backend-neutral
executable/DTO/metric/constraint contracts; they have no v1 branches.

## Comparison guarantee

`tests/dealii_diffusion_contract.cc` creates one direct v0 model and one v1
compiled model on the same triangulation, functions, constants, and bounds.
It verifies equal assembled residual, objective, objective derivative, and
DTO reduced derivative. The same test also confirms the optional compiled box
constraint and classifies unsupported matrix-free and all-at-once requests.

`tests/semantic_v1_contract.cc` independently verifies the semantic validator
for the canonical graph and for structural and analytical-policy failures.

## Exclusions

This v1 registration does not broaden the v0 executable mathematics. It does
not compile subdomain or boundary observations, Neumann/Robin terms, lifting
transformations, continuous-control bounds, non-$L^{2}$ metrics,
matrix-free execution, all-at-once/OTD, multiple equations, or multiple
optimisation variables. Each requires its own semantic declaration,
registered lowerer, capability diagnostic, and value/JVP/VJP/reduced tests.
