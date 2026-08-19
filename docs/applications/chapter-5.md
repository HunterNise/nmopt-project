# Chapter 5 application recipes

This document turns the selected Chapter 5 catalogue into application-facing
recipe contracts. It is an authoring aid: it names the semantic graph seed,
the required runtime bindings, the discrete control representation, and the
allowed option seams. The [Chapter 5 guide](../guides/chapter-5-elliptic-control.md)
remains authoritative for the mathematics, the [v1 capability table](../implementation/v1/semantic-compiler.md#registered-capabilities)
remains authoritative for realized compiler registrations, and the
[problem-library roadmap](../planning/chapter-5-problem-library-roadmap.md)
owns recipe implementation order and status.

A recipe is a typed builder of `semantic::v1::ProblemSpec`. It is not a PDE
class and it does not own deal.II functions, a mesh, a compiler, a solver, or
run output. Those are supplied through the separate public boundaries in the
[application assembly API reference](../reference/application-api.md).

The family IDs in this document are stable application-level names for the
recipe records being added. The named `make_*` functions are current semantic
graph seeds; they are not themselves the final typed recipe registry.

## Selection rules

Choose the smallest graph that expresses the requested mathematical problem.
Use the following routing table before writing a recipe:

| Need | Graph seed | Control realization | Important restriction |
| --- | --- | --- | --- |
| Full-volume scalar tracking/control | `make_scalar_diffusion_reaction_problem(with_cellwise_box)` | `FE_DGQ(0)` volume control with `l2_cellwise` metric | A cellwise box is exact only for this representation. |
| Tensor/transport/reaction/Robin terms | `make_general_scalar_elliptic_robin_problem(fixed_ids, robin_ids, with_cellwise_box)` | Same cellwise volume control | Coefficient ranks, boundary partition, and data ports are closed by the registered signature. |
| Tracking on a material subdomain | `make_subdomain_tracking_scalar_diffusion_reaction_problem(material_id, with_cellwise_box)` | Cellwise volume control | The observation region is a declared material-ID region, not an arbitrary runtime predicate. |
| Energy state tracking | `make_h1_state_tracking_scalar_diffusion_reaction_problem()` | Cellwise volume control with unchanged `l2_cellwise` search metric | The desired state has the selected $H^1_0$ target-data assumption. |
| Different continuous-control search geometry | `make_l2_metric_h1_state_tracking_continuous_control_problem()` or `make_hminus1_metric_h1_state_tracking_scalar_diffusion_reaction_problem()` | Independent homogeneous-Dirichlet `FE_Q` control | No coefficientwise box is registered for this continuous control. |
| Neumann boundary control | `make_neumann_boundary_control_problem(with_facewise_box)` | Facewise boundary control with `l2_facewise` metric | Bounds, when present, are facewise rather than cellwise. |
| Graetz/conservative transport boundary control | `make_neumann_convection_subdomain_tracking_problem(material_id, with_facewise_box)` | Facewise Neumann control | Transport and fixed/inflow/outflow policies must be bound explicitly. |
| Point sensors or normal flux | `make_point_sensor_scalar_diffusion_reaction_problem(points, fixed_ids)` or `make_normal_flux_scalar_diffusion_reaction_problem(flux_ids, fixed_ids)` | Cellwise volume control | These are explicit very-weak/strong observation registrations, not ordinary traces. |
| Dirichlet boundary control | One of the complete, partial, transposition, fractional, or tangential factories | Registered trace realization | Do not add a trace box or infer a metric from the loss. |

If no row matches, stop at the recipe boundary and report the unsupported
capability. Do not combine the nearest rows by changing a label or by
inserting a backend-specific branch into a generic solver.

## Common recipe contract

Every Chapter 5 recipe should expose only the parameters that change the
semantic problem. A useful parameter record normally contains:

- scalar coefficients, regularization, and formulation flags;
- stable region IDs and boundary/material IDs;
- control placement and the registered loss/metric choice;
- observation-region or sensor data that is part of the semantic graph; and
- an explicit constraint choice when the selected representation supports it.

The recipe may provide defaults, but it must not hide a mathematical choice
in a display string. In particular:

| Concern | Put it in | Do not put it in |
| --- | --- | --- |
| Residual terms, spaces, observations, loss, metric, transformations | `ProblemSpec` built by the recipe | `DealiiDataBindings` or solver code |
| Forcing/target functions and coefficient runtime values | `DealiiDataBindings` | Semantic component IDs or a PDE-specific executable class |
| Mesh, finite-element degree, linear-solve policies | `DealiiCompilationSession` and `DealiiDiscretisationPolicy` | The semantic graph |
| Cellwise or facewise numerical bounds | `CellwiseBoxDataBindings` or `FacewiseBoxDataBindings` | A generic `ConstraintT` implementation or an unrelated metric |
| Iteration, line search, direction, and stopping policy | Solver option record | `ProblemSpec` |
| Source revision, hardware, and run history | Experiment envelope | Recipe metadata or compilation policy |

The semantic graph must have stable IDs for every component that a binding,
manifest record, or scenario refers to. A recipe should return a complete
graph or a graph that is intentionally rejected by semantic validation; it
should not rely on a lowerer to fill missing structure.

## Recipe family contracts

The following matrix is the application-level view of the scalar family
contracts. “Runtime bindings” names the public `DealiiDataBindings` port, not
the semantic data IDs alone. A recipe record may add a declared semantic
delta to a graph seed when the source application needs a supported variant;
it must not pass a runtime binding that the graph does not declare.

| Recipe family | Chapter 5 use | Runtime bindings | Metric/constraint | Default product |
| --- | --- | --- | --- | --- |
| `scalar-diffusion-reaction-volume` | C5.1/C5.2 baseline | `forcing`, `desired_state`, `diffusion`, `reaction`, `regularisation_weight` | `l2_cellwise`; optional cellwise box | `reduced_dto` |
| `scalar-elliptic-robin-volume` | C5.5.1 general scalar/Robin composition | Baseline ports plus `general_scalar` tensor/vector/function bundle | `l2_cellwise`; optional cellwise box | `reduced_dto` |
| `scalar-subdomain-tracking` | C5.5.1 observation variant | Baseline ports; target is evaluated on the declared material region | `l2_cellwise`; optional cellwise box | `reduced_dto` |
| `scalar-energy-tracking` | C5.5.2 state energy observation | Baseline ports; desired state must satisfy the selected target-data assumption | `l2_cellwise`; no metric change | `reduced_dto` |
| `scalar-continuous-control-l2` | C5.5.2 metric comparison | Baseline ports; independent homogeneous-Dirichlet control coordinates | continuous-control `l2`; no box | `reduced_dto` |
| `scalar-continuous-control-hminus1` | C5.5.2 negative metric | Baseline ports; independent homogeneous-Dirichlet control coordinates | continuous-control `hminus1`; no box | `reduced_dto` |
| `scalar-coefficient-identification` | Physical diffusion identification | Forcing, target, reaction, regularization; diffusion is the decision variable | parameter `l2`; cellwise parameter box | `reduced_dto` |
| `scalar-neumann-boundary` | C5.7 boundary tracking/control | Baseline ports; facewise control and boundary target | `l2_facewise`; optional facewise box | `reduced_dto` |
| `scalar-neumann-convection-subdomain` | C5.6/Graetz composition | Baseline ports plus `conservative_transport`; a declared fixed-data port when nonzero | `l2_facewise`; optional facewise box | `reduced_dto` |
| `scalar-weighted-boundary-trace` | C5.7 weighted observation | Neumann ports plus `weighted_trace` | unchanged `l2_facewise`; optional facewise box | `reduced_dto` |
| `scalar-pure-neumann` | Pure-Neumann gauge variant | Neumann ports with compatible forcing/target data | `l2_facewise`; no box | `reduced_dto` |
| `scalar-point-sensor` | C5.10 finite point observations | Baseline ports; sensor coordinates are semantic region data | `l2_cellwise`; no point-specific constraint | `reduced_dto` |
| `scalar-normal-flux` | C5.8 strong normal-flux observation | Baseline ports; flux and fixed boundary IDs are semantic region data | `l2_cellwise`; no flux-specific constraint | `reduced_dto` |
| `scalar-dirichlet-trace` | C5.11 complete/partial trace control | Baseline ports; fixed data when a fixed partition is declared | Registered trace metric; no trace box | `reduced_dto` |

The default product is a routing recommendation, not a promise that every
combination accepts every `CompilationProduct`. Validate the graph, policy,
and product together before compiling.

## Runtime binding matrix

The graph determines which ports are required. The following rules keep the
runtime binding minimal and truthful:

| Graph feature | Required binding | Binding rule |
| --- | --- | --- |
| Volume forcing and target | `forcing`, `desired_state` | The function objects must outlive the compiled service and receive stable provenance strings. |
| Constant diffusion/reaction | `diffusion`, `reaction` | Use the optional diffusion scalar only when the selected graph consumes it; do not use it to stand in for tensor data. |
| Control loss | `regularisation_weight` | This scalar changes the objective, not the search metric. |
| Fixed Dirichlet reconstruction | `fixed_dirichlet_data` | Supply it only when the graph declares the fixed-data port; fixed-data precedence is a graph policy. |
| General scalar/Robin | `general_scalar` | Bind rank-specific tensor/vector/function objects for all declared coefficient ports, with matching provenance. Missing or wrong-rank data is a capability error. |
| Weighted boundary trace | `weighted_trace` | Bind the immutable boundary weight separately from the desired target and Neumann residual. |
| Conservative transport | `conservative_transport` | Bind the transport field separately; do not encode transport as a reaction or change the sign convention in the runner. |
| Cellwise box | `CellwiseBoxDataBindings` | Use scalar constants or exact-layout `dealii::Vector<double>` values with `l2_cellwise`. |
| Facewise box | `FacewiseBoxDataBindings` | Use scalar constants or exact-layout `dealii::Vector<double>` values with `l2_facewise`. |

The compiler manifest should make the same semantic ID, region, runtime
representation, and provenance visible after lowering. A recipe author must
not make a binding “optional” merely because one lowerer happens to ignore it.

## Source and mathematical authority

The [Chapter 5 guide](../guides/chapter-5-elliptic-control.md) is authoritative
for the mathematical forms, regularity assumptions, source-specific policies,
and registered catalogue variants. Use its [application catalogue](../guides/chapter-5-elliptic-control.md#application-catalogue)
when a recipe needs a source-level definition or a new mathematical choice.

This document records only the application boundary: which graph seed is
selected, which semantic delta and runtime ports are required, which discrete
control realization is supported, and which product is exposed. A new
mathematical variant must be added to the guide and capability/roadmap records
before it is introduced as a recipe family here.

## Exclusions

The current recipe library does not include state-constrained Section 5.12,
Stokes/mixed-block Section 5.13, automatic continuous-adjoint derivation,
GLS/stabilized Chapter 6 variants, or generic continuous-control box
semantics. A request for one of these must become an explicit scope decision,
not an inferred extension of a scalar recipe.

## Recipe author checklist

Before handing a recipe to an application or benchmark author, verify:

1. The recipe parameters identify one mathematical graph and have stable
   defaults.
2. The returned graph validates without backend objects.
3. Every runtime datum has a declared semantic ID, space/region, and binding
   port.
4. The control layout, loss pairing, search metric, and bound layout agree.
5. The selected product and all unsupported products are documented.
6. A manufactured default has value, JVP, VJP, and reduced-Taylor evidence.
7. The recipe can be consumed unchanged by a frozen Chapter 6 scenario.
