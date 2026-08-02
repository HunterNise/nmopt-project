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
  -> compiler::v1::CompiledProblemT<SerialBackend>
  -> generic executable/metric/constraint/DTO services
```

The v0 lowerer is not modified or replaced by this path. The v1 compiler
privately constructs a separate v0 executable instance from a validated
`ProblemSpec` and packages it behind generic executable ports. V0 and v1 can
therefore be evaluated side by side at the same coefficients, while a later
independently assembled v1 lowerer can replace the private target without
changing compiler consumers or the v0 reference.

The normative semantic protocol remains the
[interface specification](interface-specification.md). This document records
the concrete, intentionally narrow v1 realization of that protocol.

## Public semantic graph

The compatibility aggregate `include/nmopt/semantic/v1/problem_spec.hpp`
includes the focused deal.II-free headers `types.hpp`, `validation.hpp`, and
`reference_specs.hpp`. The last contains
`make_scalar_diffusion_reaction_problem()`, the current reference graph:

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

`CompilationResultT<Backend>` returns the report and only contains a
`CompiledProblemT<Backend>` when it is valid. Unsupported choices are reported
as diagnostics; the compiler does not substitute another residual,
observation, metric, constraint, or formulation.

## Registered deal.II realization

The compatibility aggregate
`include/nmopt/compiler/v1/dealii_scalar_diffusion_reaction.hpp` exposes the
focused compiler headers: `compiled_problem.hpp`, `dealii_types.hpp`,
`dealii_capabilities.hpp`, and `dealii_compiler.hpp`. The last registers only
the listed volume terms, full-domain volume restriction, quadratic losses,
`L2` metric, and optional cellwise box. Its selected discrete policy is the
v0 assembled serial realization: scalar `FE_Q` state/test with degree at
least one, `FE_DGQ(0)` control on the same mesh, homogeneous Dirichlet
boundary ids, and reduced DTO.

`CompiledProblemT<Backend>::executable_model()`, `metric()`, `constraint()`,
and `make_reduced_dto()` expose only backend-neutral ports and formulation
services. The concrete v0 model, mass metric, and box constraint stay inside
`DealiiCompiler`. Solvers have no v1 branches.

Every successful compiled product also carries a `CompilationManifest`. It
records semantic component identities, FE spaces, quadrature, the
dual-coefficient representation, data rule, metric solve tolerances,
constraint/lifting/nullspace policies, DTO provenance, and declared
assumptions. The manifest is descriptive provenance; it does not create a
second configuration channel.

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
