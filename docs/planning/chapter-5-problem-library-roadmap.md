# Chapter 5 problem library roadmap

## Purpose

This roadmap owns the reusable problem recipes derived from Chapter 5. It is
separate from the [Chapter 5 elliptic-control guide](../guides/chapter-5-elliptic-control.md),
which records the mathematical catalogue, and from the
[Chapter 6 benchmark suite roadmap](chapter-6-benchmark-suite-roadmap.md),
which records frozen numerical experiments.

The problem library is intended to showcase the framework and provide
parameterized standard problems for exploration. A library entry is therefore
not a new PDE/control class. It is a named convenience recipe that creates a
compositional semantic specification from parameters such as:

- coefficients, forcing, target data, and regularization;
- mesh, finite-element, boundary-region, and observation-region choices;
- control placement and control metric;
- optional constraints supported by the selected representation;
- state, adjoint, and metric solve policies.

The recipe must lower through the same residual, objective, observation,
metric, constraint, and solver services used by an explicitly composed
problem. It may provide friendly defaults and validation, but it must not
duplicate PDE assembly or optimization logic.

## Library versus benchmark

The same recipe can produce both an exploratory application and a benchmark
configuration, but the two uses have different contracts.

| Problem-library recipe | Benchmark scenario |
| --- | --- |
| Parameters are intentionally changeable. | Parameters are frozen and versioned. |
| Demonstrates composition and usability. | Tests reproducibility and regression behavior. |
| May expose many valid combinations supported by the compiler. | Selects one exact discretization, data set, and algorithm. |
| Success means a valid compiled problem and useful diagnostics. | Success means derivative/KKT correctness plus agreed numerical trends. |
| Exploratory plots and parameter sweeps are welcome. | Outputs and tolerances must be declared before the run. |

Benchmark reproduction is not the only correctness evidence: value, JVP,
VJP, reduced-Taylor, and KKT contract tests remain authoritative for local
mathematical correctness. Benchmarks provide system-level validation that the
compiled components work together and that selected source results are
reproduced.

## Accepted library families

The initial library covers the selected scalar Chapter 5 scope:

| Family | Chapter 5 variants | Role in the library |
| --- | --- | --- |
| Scalar distributed control | C5.1, C5.2, and C5.5 variants | Baseline elliptic, transport, reaction, Robin, tracking, and selected cellwise-box problems. |
| Neumann boundary control | C5.6 and C5.7 | Standard facewise control with volume or boundary observations. |
| Observation variants | C5.8 and C5.10 | Normal-flux and finite point-sensor problems under their explicit policies. |
| Dirichlet boundary control | C5.11 variants | Complete or selected fixed/controlled trace transformations with the registered control losses and metrics. |
| General scalar/Robin composition | C5.5.1 and C5.5.2 | Parameterized tensor diffusion, transport, reaction, Robin data, and energy tracking. |

Sections 5.12 and 5.13 are excluded from this library. They would require
state-constraint/multiplier and mixed Stokes problem families, respectively.

## Parameterization policy

Each recipe should expose parameters at the level at which users naturally
experiment:

- scalar coefficients and transport fields;
- forcing and target functions;
- named boundary and observation regions;
- control boundary or volume region;
- loss and search metric choices that are registered for the family;
- regularization and, where supported, cellwise or facewise bounds;
- mesh refinement and solver tolerances.

Parameters that change the mathematical formulation must be part of the
semantic specification and manifest. Parameters that only affect an
experiment, such as iteration limits or output frequency, belong to the
solver or experiment configuration and must not alter the compiled PDE.

Every recipe should have:

- a default manufactured configuration suitable for fast contract tests;
- a documented parameter sweep for exploratory use;
- explicit diagnostics for unsupported combinations;
- a comparison against direct composition for at least one configuration.

## Implementation sequence

### L0 — Recipe boundary and configuration records

Define the public boundary between a standard-problem recipe and the semantic
compiler. A recipe returns a `ProblemSpec` or equivalent declarative request;
it does not return a PDE-specific executable model. Keep scenario parameters,
source provenance, and mesh/data bindings structured rather than encoded in
display strings.

### L1 — Scalar baseline recipes

Provide parameterized distributed-control and general scalar/Robin recipes.
Use them to exercise diffusion, transport, reaction, subdomain tracking,
energy tracking, and selected cellwise boxes without creating one factory per
cross-product.

### L2 — Neumann boundary-control recipes

Expose volume-tracking and boundary-tracking Neumann controls, including the
convection/subdomain composition used by the Graetz benchmark. The recipe
must make control boundary, orientation, transport convention, and observation
region explicit.

### L3 — Dirichlet-control recipes

Expose the registered complete and partial fixed/controlled trace variants,
including the selected Section 5.11 losses and metrics. Keep the physical
lifting/transformation in the compiled formulation rather than in the recipe
or optimizer.

### L4 — Exploration and documentation layer

Add small parameter-sweep drivers and field/objective reporting that consume
compiled problem services. The drivers must not become a second numerical
lowerer. Link every standard recipe to the relevant Chapter 5 entry and to
any Chapter 6 benchmark that freezes it.

## Library acceptance

The library is ready for a family when:

1. its default request validates and compiles through the selected path;
2. changing an allowed parameter changes only the declared semantic/data
   components;
3. residual/objective/reduced derivative tests pass for the default;
4. unsupported cross-products fail with a capability diagnostic;
5. at least one recipe configuration is consumed unchanged by the benchmark
   harness.

The library is not required to expose every theoretical Chapter 5 variant.
Its purpose is a coherent, inspectable set of standard scalar problems that
can be recombined and varied without proliferating application classes.
