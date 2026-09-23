# Validation, resolution, and capabilities

## A valid problem is not yet an executable problem

The previous chapter built a `ProblemSpec` from mathematical ingredients:

```text
PDE + objective + spaces + observations + metric + constraints
        │
        ▼
semantic graph
```

At that point the graph may look complete, but several different questions remain:

- Does every reference point to something that exists?
- Are the declared spaces, roles, observations, and policies mutually coherent?
- Have all choices needed by later compilation been made explicitly?
- Does the current deal.II compiler know how to realize the requested combination?
- Can it construct the particular formulation product being requested?

Those questions belong to different layers. Treating all of them as one vague
"validation" step would make the architecture much harder to reason about.

The current path is better read as

```text
ProblemSpec
    │
    ▼
semantic validation
    structural + analytical-policy checks
    │
    ▼
semantic resolution
    validated graph + stable-ID lookup view
    │
    ▼
compiler request resolution / closure
    choose the registered target family and realization facts
    │
    ▼
lowerability checks
    can this graph be realized by the current deal.II lowerers?
    │
    ▼
formulation / product capability checks
    can the requested numerical product be constructed?
    │
    ▼
numerical realization
```

The boundary considered here stops before the accepted request is lowered into
finite-element objects; that construction is the subject of Chapter 11.

## 1. Four diagnostic categories separate four kinds of failure

`nmopt` reports validation and compiler failures through one shared diagnostic
vocabulary:

```text
DiagnosticCategory
    structural
    analytical_policy
    lowerability
    formulation_capability
```

These names are more useful when attached to concrete questions.

| Category | Question being answered | Typical failure |
| --- | --- | --- |
| `structural` | Is the semantic graph internally coherent? | missing/duplicate ID, invalid reference, mismatched region/space/role |
| `analytical_policy` | Has a required mathematical or realization assumption been declared? | point evaluation or transposition policy missing |
| `lowerability` | Can the current deal.II compiler realize this valid request? | unsupported registered signature, incompatible metric realization, missing concrete binding |
| `formulation_capability` | Can the selected formulation/product be constructed from that realization? | unsupported KKT/PDAS shape, unsupported constraint/metric combination |

A single `ValidationReport` carries diagnostics from all four categories.

Its logical contract is simple:

```text
ValidationReport
    diagnostics : vector<Diagnostic>

Diagnostic
    category       // one of the four categories above
    component_id   // semantic/compiler component responsible for the problem
    capability     // stable name for the missing/violated capability
    remedy         // human-readable corrective action
```

The report is valid exactly when it contains no diagnostics.

### 1.1 One report does not mean one validator owns every category

The semantic validator emits only semantic-level failures.

The compiler appends lowerability and formulation/product diagnostics later.

This is an important architecture boundary:

```text
SemanticValidator
    structural
    analytical_policy

DealiiCompiler
    lowerability
    formulation_capability
```

The categories share one report so callers can handle diagnostics uniformly, not
because all checks belong to one layer.

## 2. `SemanticValidator` – is the semantic graph coherent?

`SemanticValidator` is deal.II-free.

It consumes only a `ProblemSpec` and asks whether the semantic description is
self-consistent.

Conceptually:

```cpp
SemanticValidator validator;
ValidationReport report = validator.validate(specification);

if (!report.valid())
  // do not proceed to semantic resolution or compiler lowering
```

The validator first constructs temporary ID indexes for the graph components, then
checks labels, required enum selections, regions, spaces, pairings, variables, data,
transformations, equations, residual terms, observations, losses, metrics,
constraints, formulation references, supplied-OTD declarations, and requirement
policies.

It does not assemble a matrix and it does not ask whether deal.II has an
implementation for the graph.

### 2.1 Structural validity begins with identity

Every semantic component has a stable ID.

Therefore the first invariants are mundane but fundamental:

```text
id is non-empty
id is unique within its component kind
label is non-empty
required enum-valued fields are not "unspecified"
```

Without those conditions, references such as

```text
state.space_id = "state_space"
loss.source_observation_id = "state_observation"
formulation.metric_id = "control_l2_metric"
```

cannot be interpreted reliably.

A duplicate ID is not a compiler limitation. It is an ill-formed semantic graph, so
the correct category is `structural`.

## 3. Structural validation – references must agree with mathematical roles

A graph can contain only valid IDs and still be incoherent.

Suppose an observation says

```text
input_variable_id = state
output_space_id   = sensor_space
region_id         = sensor_region
```

The validator must check that those references exist and that their roles make sense
together.

For a point observation

```math
\mathcal O_{\mathrm s}(y)
=
\begin{bmatrix}
y(x_{1})\\
\vdots\\
y(x_{m})
\end{bmatrix},
```

the semantic structure contains at least

```text
point-set region
    X_s = {x_1, ..., x_m}

observation space
    dimension = m

point_sensor observation
    state → observation space
```

The condition

$$
\dim(\text{observation space})
=
\lvert X_{\mathrm s}\rvert
$$

is a semantic consistency condition. It does not depend on how deal.II will later
evaluate the finite-element function at the points.

### 3.1 Structural checks can involve several nodes

The validator therefore checks relationships, not just individual fields.

Examples include:

```text
SpaceSpec.region_id
    must refer to a suitable RegionSpec

VariableSpec.space_id
    must refer to an existing SpaceSpec

PairingSpec
    must refer to existing primal/covector spaces

ResidualTermSpec.equation_id
    must refer to the equation containing the term

EquationBlockSpec.residual_term_ids
    must name residual terms that belong to that equation

ObservationSpec
    must connect a valid variable, region, output space, and pairing

LossSpec
    must connect a valid observation, datum, and pairing

MetricSpec
    must apply to the selected variable with a coherent pairing

ConstraintSpec
    must bind bounds to the constrained variable

FormulationSpec
    must point to the declared state/decision/equation/metric/constraint
```

The same principle applies to supplied OTD: the state, adjoint, and
control-stationarity block declarations must refer to coherent semantic spaces and
pairings before the compiler is allowed to ask whether that OTD system is supported.

## 4. `analytical_policy` – the graph may need assumptions beyond its nodes

Some mathematical statements cannot be inferred merely from the existence of a
space, term, or observation.

A point observation is again a useful example.

Writing

$$
\mathcal O_{\mathrm s}(y)
=
[y(x_{1}),\ldots,y(x_{m})]^{\mathsf T}
$$

does not, by itself, explain the analytical/discrete assumptions under which point
evaluation and the corresponding transpose action are to be interpreted.

The reference point-sensor graph therefore carries requirement policies for both the
point evaluation and the transposition realization.

If those policies are absent, the semantic graph is rejected with
`analytical_policy` diagnostics.

### 4.1 Structural failure and policy failure are intentionally different

Consider two broken point-sensor graphs.

In the first, the transposition selection refers to no diffusion data at all:

```text
transposition selection
    diffusion_data_id = ""
```

The graph contains an incomplete typed relationship. That is structural.

In the second, the graph contains a valid point-sensor observation but omits the
required point-evaluation policy entirely.

The graph shape is understandable, but a required assumption/realization commitment
has not been supplied. That is an analytical-policy failure.

The distinction is:

```text
structural
    "the declaration does not form a coherent graph"

analytical_policy
    "the graph is coherent enough to understand,
     but a required assumption or declared realization is missing"
```

### 4.2 `user_assumed` is a declaration, not a proof

A requirement with

```text
status = user_assumed
```

means that the model author is asserting an analytical condition.

For example,

```math
a(v,v)
\geq
c\lVert v\rVert^{2}
```

may be declared as a coercivity assumption.

`SemanticValidator` can check that the required policy is present, correctly scoped,
and attached to the right subject.

It does not prove coercivity from the coefficients.

That responsibility remains with the mathematical model author or external evidence.

## 5. Typed selections carry semantics; display strings do not

`RequirementPolicySpec` contains both descriptive text and, for important choices,
typed selection records.

This distinction matters.

A policy may include

```text
selected_policy
    human-readable explanation
```

together with a typed field such as

```text
typed_transposition_selection
typed_fractional_metric_selection
typed_partial_boundary_selection
...
```

The validator and compiler reason from the typed selection.

They do not parse the prose in `selected_policy` to decide what the problem means.

The semantic contract tests exercise this explicitly: clearing the display text from
the point-sensor policies does not invalidate an otherwise complete typed selection.

That gives the project a useful rule:

> prose explains a choice; typed fields determine it.

## 6. `SemanticResolver` – validation plus stable-ID lookup

After semantic validation succeeds, `SemanticResolver` constructs a
`ResolvedProblemView`.

The current resolver is intentionally small:

```text
ProblemSpec
    │
    ├── SemanticValidator
    │       │
    │       └── invalid --> diagnostics only
    │
    ▼
ResolvedProblemView
    validated graph + stable-ID lookup tables
```

`ResolvedProblemView` indexes the semantic component vectors by ID:

```text
region(id)
space(id)
pairing(id)
variable(id)
datum(id)
transformation(id)
residual_term(id)
equation(id)
observation(id)
loss(id)
metric(id)
constraint(id)
requirement(id)
```

The point is not to invent new semantics.

It is to replace repeated searches such as

```text
"find the one SpaceSpec whose id is control_space"
```

with a validated lookup.

### 6.1 Resolution does not copy or own the problem

The view borrows the original `ProblemSpec`.

That means

```text
ProblemSpec
    owns semantic component storage

ResolvedProblemView
    stores pointers/indexes into that ProblemSpec
```

The view is intentionally confined to one validation/compilation operation.

This matters for lifetime reasoning: semantic resolution is not a hidden ownership
transfer.

### 6.2 `SemanticResolution` combines diagnostics and the optional view

The resolver returns

```text
SemanticResolution
    diagnostics : ValidationReport
    problem     : optional<ResolvedProblemView>
```

and succeeds only when

```text
diagnostics.valid()
and
problem.has_value()
```

A failed semantic validation therefore never produces a resolved view.

## 7. There are two distinct meanings of "resolution" in the current path

The word *resolution* appears at two different boundaries, and distinguishing them
prevents confusion.

### 7.1 Semantic resolution resolves graph references

`SemanticResolver` gives the compiler a validated, indexable semantic graph.

It answers questions such as

```text
what SpaceSpec is named by this ID?
what RegionSpec does this observation use?
what MetricSpec belongs to the formulation?
```

It does not yet choose a deal.II target family.

### 7.2 Compiler request resolution closes realization choices

Once semantic resolution succeeds, `DealiiCompiler` derives a
`ResolvedCompilationRequest`.

That object is compiler-specific.

It records facts such as

```text
target family
Dirichlet-control registration
whether fixed reconstruction is used
whether the control is Neumann / Dirichlet / continuous volume
whether the metric is H1, H1/2, H-1, or L2
whether the observation is point-sensor / normal-flux / subdomain / ...
selected boundary regions
selected typed realization records
required concrete data-binding ports
```

The progression is therefore more precisely

```text
ProblemSpec
    │
    ▼
SemanticResolver
    │
    ▼
ResolvedProblemView
    semantic graph closed under stable-ID lookup
    │
    ▼
compiler request resolution
    │
    ▼
ResolvedCompilationRequest
    compiler-specific realization facts closed before lowering
```

The first object belongs to the semantic namespace.

The second belongs to the compiler namespace.

## 8. `ResolvedCompilationRequest` – a closed compiler-layer request

A compiler should not rediscover the same cross-product of semantic predicates every
time it constructs a numerical object.

`ResolvedCompilationRequest` collects the decisions that have already been made from
the validated graph.

Its central discriminator is

```text
ResolvedTargetFamily
```

with current families such as

```text
direct_volume
assembled_volume
neumann_boundary
weighted_boundary_trace
pure_neumann
dirichlet_control
l2_dirichlet_transposition
hhalf_dirichlet_control
h1_tracking_hhalf_dirichlet_control
h1_dirichlet_control
h1_control_l2_metric
h1_control_h1_metric
hminus1_control_metric
continuous_control_l2_metric
coefficient_identification
general_scalar_robin
point_sensor
normal_flux
```

The request also retains independent realization facts instead of encoding everything
in one target-name string.

For example,

```text
uses_point_sensor = true
transposition_selection = ...
point_sensor_evaluation_policy = ...
```

can remain explicit even after `target_family = point_sensor` has been chosen.

### 8.1 Why close the request before construction?

Suppose the semantic graph requests an $H^{1/2}$ boundary metric.

The compiler should not reach a late matrix-construction function and decide, based
on whatever data happen to be available, which fractional metric realization it will
use.

Instead, the relevant typed selection is resolved before lowering:

```text
semantic metric
    H1/2

requirement policy
    selected fractional realization

ResolvedCompilationRequest
    fractional_metric_selection
        │
        ▼
lowerability check
        │
        ▼
metric construction
```

This keeps the effective numerical experiment inspectable.

## 9. Compiler validation starts only after semantic resolution succeeds

The public compiler validation path follows this sequence:

```text
SemanticResolver::resolve
        │
        ├── invalid --> return semantic diagnostics
        │
        ▼
resolve_compilation_request
        │
        ▼
resolve / close registered target choices
        │
        ▼
validate_lowerability
        │
        ▼
validate_formulation_capability
        │
        ▼
validate supplied-OTD / target / product capability
```

The early return is important.

A malformed semantic graph does not continue into compiler matching, where a missing
ID might otherwise be misreported as an unsupported deal.II feature.

Semantic errors stay semantic.

Compiler errors stay compiler errors.

## 10. `lowerability` – can the current compiler realize the valid graph?

A semantically valid request can still be outside the current deal.II compiler.

That is what `lowerability` means.

The question is not

> Is this mathematical problem meaningful?

It is

> Does the current registered compiler know how to turn this specific semantic
> request into concrete numerical objects?

### 10.1 A semantic kind can exist before every combination is supported

Chapter 9 introduced semantic enums such as

```text
ObservationKind::point_sensor
MetricKind::hhalf
ResidualTermKind::neumann_control
```

Those names belong to the semantic vocabulary.

Their existence does not imply that every possible combination involving those kinds
is executable.

For example, a graph might combine

```text
point-sensor state tracking
H1/2 boundary-control metric
Neumann control
PDAS product
```

in a way that is individually understandable at the semantic level.

The current compiler is not obligated to invent a realization for that cross-product.

### 10.2 The capability registry is a kind-level ledger, not a universal composer

`DealiiCapabilityRegistryV1` records whether the compiler has a registered lowerer
capability for individual semantic kinds:

```text
has_residual_term_lowerer(...)
has_observation_lowerer(...)
has_loss_lowerer(...)
has_metric_lowerer(...)
has_constraint_lowerer(...)
has_transformation_lowerer(...)
```

For example, the current registry knows about all of the semantic observation kinds
listed in Chapter 9, including point sensors and normal flux.

That still does not mean

```text
registered observation kind
+
registered metric kind
+
registered control kind
=
automatically supported problem
```

The combination must also match a coherent lowering plan or a registered target
strategy.

This distinction is central to the current compiler design.

## 11. Bounded composition – registered pieces, registered combinations

The scalar compiler demonstrates the idea particularly clearly.

For the baseline distributed-control residual,

```math
E(y,u)[v]
=
\underbrace{
\int_{\Omega}
\kappa\nabla y\cdot\nabla v
+
cyv
\mathrm{d}x
}_{\texttt{diffusion\_reaction}}
-
\underbrace{
\int_{\Omega}fv\mathrm{d}x
}_{\texttt{volume\_source}}
-
\underbrace{
\int_{\Omega}uv\mathrm{d}x
}_{\texttt{volume\_control}},
```

the bounded scalar planner can produce residual contributions for the three semantic
terms.

The residual assembly layer recognizes the closed registration

```text
diffusion_reaction
volume_source
volume_control
```

as one supported scalar residual family.

A more general registered family contains

```text
tensor_diffusion
conservative_transport
advective_transport
reaction
volume_source
volume_control
robin_bilinear
robin_source
```

and is recognized separately.

A random subset or recombination of those eight entries is not silently interpreted
as "close enough".

If no registered residual signature applies, the request is not lowered through the
nearest available target.

### 11.1 Composition is real, but bounded

This is why the compiler is neither

```text
one giant switch over whole applications
```

nor

```text
an unrestricted symbolic weak-form compiler.
```

It can reuse registered handlers for

```text
residual terms
observations
losses
metrics
constraints
transformations
```

while still requiring the completed plan to match a known coherent realization.

## 12. A closed signature prevents accidental cross-products

The Dirichlet-control registrations give another concrete example.

The current compiler recognizes several distinct combinations of

```text
state tracking topology
control observation/loss
search metric
boundary realization
```

such as

```text
L2 state tracking
+
L2 boundary control loss
+
L2 metric
+
complete nodal controlled boundary
```

or

```text
H1 state tracking
+
L2 boundary control loss
+
H1/2 metric
+
declared trace realization
```

These are separate registered signatures.

Changing just one axis does not automatically create a new compiler capability.

For example, replacing the second signature's $H^{1/2}$ metric by an arbitrary other
metric is not accepted merely because that other metric exists elsewhere in the
compiler.

The failure is a lowerability failure:

```text
semantic graph
    coherent

individual semantic kinds
    known

combined registered realization
    absent
```

That is exactly the case that a closed capability table is meant to expose.

## 13. Typed realization selections make lowerability explicit

Some semantic choices need an additional typed realization before numerical
construction can proceed.

The $H^{1/2}$ metric is a good example.

At the mathematical level, the semantic graph asks for

$$
G_{1/2}:U_{h}\longrightarrow U_{h}^{\ast}.
$$

Chapter 4 showed one realization through a minimum-energy extension and Schur
complement.

The current compiler does not infer such a construction merely from

```text
MetricKind::hhalf
```

alone.

It also expects the matching typed fractional-trace selection.

The lowerability question is therefore:

```text
semantic metric kind
    hhalf
        │
        ▼
typed realization selection
    operator realization
    apply realization
    inverse realization
        │
        ▼
registered compiler signature?
```

For the current registered path, the compiler expects the selected
minimum-extension / Schur-complement realization and its corresponding inverse.

A different typed choice may still be mathematically meaningful, but it is not
automatically lowerable.

### 13.1 Point sensors require both semantic and compiler-side closure

For a point observation

$$
\mathcal O_{\mathrm s}(y)
=
[y(x_{1}),\ldots,y(x_{m})]^{\mathsf T},
$$

semantic validation checks the graph and required policies.

Compiler lowerability then checks that the resolved transposition selection matches
the registered point-sensor realization.

The progression is:

```text
point_sensor ObservationSpec
        │
        ├── point-set region
        ├── finite observation space
        ├── evaluation policy
        └── transposition policy
        │
        ▼
semantic validation succeeds
        │
        ▼
ResolvedCompilationRequest
    point-sensor target
    typed transposition selection
        │
        ▼
lowerability
    registered point-sensor realization matches?
```

Thus the same mathematical concept appears at several levels without being duplicated:

```text
math
    O_s(y)

semantic graph
    what observation is requested?

policy
    what assumptions/realization commitments accompany it?

compiler request
    which registered realization was selected?

lowerer
    how is the FE operator actually built?
```

## 14. Runtime bindings can introduce lowerability failures too

Not every compiler decision can be made from `ProblemSpec` alone.

The semantic graph may say

```text
forcing
    DataKind::function

diffusion
    DataKind::tensor_function

observation weight
    DataKind::function
```

but the compiler eventually receives concrete deal.II bindings.

At that point it can check facts such as

```text
is the required binding present?
does it have the expected number of components?
does it carry provenance?
does its concrete runtime representation match the resolved port?
```

A missing or malformed concrete binding is still a lowerability problem.

The semantic graph was coherent; the current numerical realization cannot be
constructed from the concrete inputs that were supplied.

### 14.1 Static compiler validation and compilation are not identical

`DealiiCompiler::validate(...)` receives

```text
ProblemSpec
DealiiDiscretisationPolicy
CompilationProduct
```

but not the full concrete mesh and runtime data bindings used by `compile(...)`.

Therefore static validation can reject unsupported semantic/policy combinations, but
some errors necessarily appear only when compilation sees the actual session, mesh,
or bound data.

For example, compilation can reject:

```text
null compilation session
missing required Function binding
empty binding provenance
multi-component Function where one scalar component is required
unsupported concrete reference-cell / target combination
```

The diagnostic category remains `lowerability` because the failure concerns concrete
realization, not semantic meaning.

## 15. `formulation_capability` – can the requested numerical product be built?

Lowerability answers whether the numerical ingredients can be realized.

Formulation capability asks whether those ingredients can be assembled into the
requested solver-facing formulation product.

These are related but different questions.

Suppose a compiler can realize

```text
state space
control space
PDE residual
L2 metric
cellwise box
```

individually.

It does not follow that every product built from those ingredients is implemented.

### 15.1 The base DTO capability remains intentionally narrow

For the current generic DTO formulation contract, the compiler expects one state,
one primary decision variable, and one equation block.

Schematically:

```text
state       y
decision    u or parameter m
equation    E(y,u)=0
```

A valid semantic graph with several coupled state equations could therefore be
structurally meaningful while lying outside this formulation capability.

The diagnostic is not "bad PDE".

It is

```text
current executable DTO product
    does not support this block shape
```

### 15.2 Constraints also interact with formulation capability

A constraint may be lowerable by itself but unsupported in a particular formulation.

For the current projected-gradient surface, the registered box types are the
coefficientwise cellwise or facewise $L^{2}$ boxes.

So the distinction is:

```text
ConstraintSpec is coherent
        │
        ▼
constraint lowerer exists
        │
        ▼
selected formulation knows how to use this constraint?
```

The last question belongs to `formulation_capability`.

## 16. Compiler products are a second axis beside semantic formulation

The semantic graph contains a `FormulationSpec`.

The compiler also receives a `CompilationProduct`.

Current compiler product requests are:

| `CompilationProduct` | Requested compiled surface |
| --- | --- |
| `reduced_dto` | reduced DTO formulation services |
| `quadratic_kkt` | quadratic equality-constrained KKT product |
| `pdas` | bound-constrained PDAS product built on the registered KKT machinery |

This enum should not be confused with `FormulationKind`.

The semantic formulation says what first-order organization/provenance the problem
declares.

The compilation product says which concrete solver-facing product the caller wants
the compiler to package.

### 16.1 Supplied OTD is represented through semantic provenance

A supplied OTD graph declares

```text
FormulationKind::all_at_once
FormulationProvenance::supplied_otd
SuppliedOTDDeclaration
```

rather than appearing as a separate `CompilationProduct` enum value.

The compiler then checks whether the supplied declaration is compatible with the
requested product.

This is why the compiler validation pipeline has a dedicated
supplied-OTD capability step in addition to the generic formulation/product checks.

## 17. KKT capability is narrower than "a KKT matrix can be written"

Part II derived the abstract equality-constrained quadratic system

```math
\begin{bmatrix}
Q & D^{\mathsf T}\\
D & 0
\end{bmatrix}
\begin{bmatrix}
x\\
\lambda
\end{bmatrix}
=
\begin{bmatrix}
c\\
d
\end{bmatrix}.
```

Many numerical problems can be written in this form.

The current compiler does not claim to construct an arbitrary KKT product from any
such semantic graph.

Its registered quadratic-KKT product is tied to the supported canonical scalar DTO
or supplied-OTD target and to the assumptions encoded by that path.

Thus

```text
abstract KKT mathematics exists
```

is not the same statement as

```text
current compiler can construct this KKT product.
```

That latter claim is what `formulation_capability` tracks.

## 18. PDAS adds another layer of capability conditions

PDAS needs more than a lowerable box constraint.

The current compiled PDAS product requires a specific combination:

```text
registered distributed scalar target
cellwise L2 box
positive-diagonal cellwise L2 metric
valid PDAS policy
declared active-row rank assumption
declared kernel-positivity assumption
valid inner KKT solver policy
```

The reason is visible in the mathematics from Chapter 8.

The active-set method repeatedly constructs restricted equality-constrained KKT
systems. It also needs a primal representation of the box multiplier, which depends
on the metric realization.

If the metric were replaced by an arbitrary coupled $H^{1/2}$ metric, merely having a
valid Riesz map would not make the current coefficientwise PDAS multiplier conversion
correct.

Therefore a request can pass

```text
semantic validation
lowerability of metric
lowerability of box
```

and still fail

```text
compiled_pdas_metric_realisation
```

at the formulation-capability layer.

This is an example of a failure that would be obscured if every compiler rejection
were simply called "unsupported".

## 19. Supplied OTD has its own capability checks

Chapter 9 showed the supplied first-order residual

```math
F_{h}(y,p,u)
=
\begin{bmatrix}
Ay-f-Bu\\
A^{\mathsf T}p-M_{y} y+q\\
B^{\mathsf T}p+\beta N_{u} u
\end{bmatrix}.
```

A structurally valid `SuppliedOTDDeclaration` can name all three blocks correctly and
still be outside the current supplied-OTD lowerer.

The compiler currently checks, among other things, that the registered scalar target,
spaces, pairings, and multiplier convention match the supplied-OTD capability.

For the canonical registered path, the multiplier convention is the framework
adjoint with identity conversion.

So these are three different statements:

```text
the application has supplied an OTD declaration

the declaration is structurally coherent

the current compiler has a lowerer/product for that supplied declaration
```

Only the last one is a formulation/compiler capability claim.

## 20. One diagnostic report supports a useful failure ladder

The separation between categories is visible by following one problem as its request
becomes progressively more complete.

Consider point tracking.

### 20.1 Incomplete graph

Suppose the point observation names a missing output space.

Result:

```text
category
    structural

meaning
    the graph itself cannot be interpreted coherently
```

Fix the reference.

### 20.2 Missing analytical policy

Now the point observation is structurally complete, but the required point-evaluation
policy is absent.

Result:

```text
category
    analytical_policy

meaning
    the mathematical/realization assumption has not been declared
```

Add the required policy.

### 20.3 Unsupported typed realization

Now the graph and policies are valid, but the typed transposition realization does not
match the one registered for the current point-sensor lowerer.

Result:

```text
category
    lowerability

meaning
    the current compiler cannot realize this valid request
```

Select a registered realization.

### 20.4 Unsupported solver-facing product

Finally suppose the point-sensor target is lowerable, but the caller asks for a
compiled product that is not registered for this target.

Result:

```text
category
    formulation_capability

meaning
    numerical realization exists, but not in the requested formulation/product shape
```

The four outcomes preserve which layer must change instead of collapsing every
failure into one generic exception.

## 21. Capability support is explicit, not "nearest match"

A bounded compiler must be predictable when it encounters a graph outside its
registered set.

The current rule is:

> reject an unsupported signature rather than substituting a nearby realization.

For example, if one registered Dirichlet-control signature uses

```text
H1 state tracking
L2 control loss
H1/2 search metric
```

then changing the search metric creates a different request.

The compiler does not silently keep the old $H^{1/2}$ construction because it looks
similar.

Likewise, if a residual assembly registration expects

```text
tensor diffusion
conservative transport
advective transport
reaction
source
volume control
Robin bilinear term
Robin source
```

then deleting one term does not mean "use the general scalar target with that term
turned off" unless such a combination is explicitly represented by the registered
planning path.

This behavior is important for reproducibility.

A compiler that silently chooses the nearest implementation would make the actual
numerical problem depend on undocumented matching heuristics.

## 22. The capability ledger and the capability registry answer different questions

The project has two related but different sources of capability information.

`DealiiCapabilityRegistryV1` is a small code-level registry used during diagnostics.
It answers questions such as

```text
is there any registered residual lowerer for this ResidualTermKind?
is there any registered observation lowerer for this ObservationKind?
```

The current [Compiler reference](../../reference/compiler.md) records
supported realization families, factory routes, products, and common rejection
boundaries. For one concrete graph, policy, and product,
`DealiiCompiler::validate()` is the public whole-request capability check.

Those are therefore the better answers to

> Can the current compiler build this whole problem?

The registry is useful evidence for a narrower question:

> Is this individual semantic kind known to the compiler at all?

### 22.1 Individual support does not imply Cartesian-product support

Suppose the registry reports support for

```text
point_sensor
hhalf metric
facewise_box
```

individually.

It is invalid to infer that

```text
point_sensor × hhalf × facewise_box
```

is automatically a registered target.

The valid combinations are constrained by the lowering plans, target-family
registrations, typed realization selections, and formulation/product checks.

## 23. Capability is not the same as benchmark coverage

A compiler capability may exist even when it has not received the same end-to-end
benchmark treatment as another capability.

The current compiler surface includes reduced DTO, supplied OTD, quadratic KKT, and
PDAS-related products for selected registered problems.

The repository's named benchmark/reproduction programme exercises those surfaces
unevenly.

Therefore keep two questions separate:

```text
framework capability
    can the compiler construct and verify this product?

benchmark coverage
    has this product been exercised as a named application/reproduction campaign?
```

A missing benchmark does not erase an accepted compiler capability.

Conversely, a benchmark for one target does not prove that neighboring semantic
combinations are supported.

## 24. Diagnostics should identify the failed layer and the remedy

Each `Diagnostic` carries

```text
category
component_id
capability
remedy
```

The `capability` field is especially useful for tests and tooling because it is more
stable than matching an entire prose error message.

For example, a diagnostic may identify

```text
component_id
    control_l2_metric

capability
    compiled_pdas_metric_realisation

remedy
    select the registered positive-diagonal cellwise L2 metric
```

The diagnostic tells the reader three things:

```text
where?
    the selected metric

what?
    the PDAS metric realization capability

how to proceed?
    choose the registered realization
```

That is more actionable than a generic "compilation failed".

### 24.1 Diagnostics are evidence, not automatic repair

The compiler reports a remedy but does not rewrite the semantic graph to make it
compile.

If an $H^{1/2}$ metric is requested and only a different metric would make one target
lowerable, the compiler does not silently change the metric.

The user or problem recipe must make that mathematical/numerical decision.

## 25. What each layer is allowed to know

The boundaries can be summarized compactly.

| Layer | It knows | It should not decide |
| --- | --- | --- |
| `ProblemSpec` | mathematical roles, regions, spaces, terms, observations, metrics, constraints, declared policies | deal.II objects, solver tolerances, target-matching heuristics |
| `SemanticValidator` | graph coherence and required semantic/analytical policies | whether a deal.II implementation exists |
| `ResolvedProblemView` | validated IDs and direct component lookup | target family or FE realization |
| `ResolvedCompilationRequest` | compiler-specific target/realization facts derived from the graph | optimization algorithm policy |
| lowerability validation | registered compiler signatures, concrete realization restrictions, binding/mesh requirements | changing the requested mathematics to fit |
| formulation/product validation | supported reduced/KKT/PDAS/supplied-OTD product shapes | how the optimizer iterates after product construction |

This division prevents two common architectural mistakes.

The first is pushing compiler limitations into semantic validity:

```text
"the current deal.II compiler cannot do it"
        ≠
"the semantic problem is meaningless"
```

The second is allowing a semantically valid graph to reach numerical construction
without first closing the choices on which that construction depends.

## 26. A practical authoring workflow

When adding or modifying a semantic problem family, the progression is:

```text
1. write the mathematical problem

2. build the ProblemSpec
       regions / spaces / variables / terms / observations / ...

3. make semantic validation pass
       structural
       analytical_policy

4. resolve the graph
       stable IDs become safe lookups

5. inspect compiler lowerability
       target family
       typed realization selections
       registered component/target signatures

6. select the requested CompilationProduct

7. make formulation/product capability checks pass

8. bind runtime data and compile
       remaining mesh/data lowerability checks happen here
```

A failure at step 3 should generally be fixed in the semantic graph.

A failure at step 5 may require either

```text
choosing an already registered realization
```

or

```text
implementing a new compiler capability.
```

Those are very different changes and should not be conflated.

## 27. Three examples across the layers

### 27.1 Distributed volume control

Mathematics:

```math
-\nabla\cdot(\kappa\nabla y)+cy=f+u.
```

Semantic structure:

```text
diffusion_reaction
volume_source
volume_control
L2 control metric
reduced DTO formulation
```

Semantic validation checks the graph connections and policies.

Compiler resolution identifies the direct/assembled scalar target and required data
bindings.

Lowerability checks whether the selected control realization, metric, mesh policy,
and bindings fit the registered target.

The reduced product is then checked against the formulation capability.

### 27.2 Point-sensor tracking

Mathematics:

```math
\mathcal O_{\mathrm s}(y)
=
[y(x_{1}),\ldots,y(x_{m})]^{\mathsf T}.
```

Semantic structure:

```text
point-set region
m-dimensional observation space
point_sensor observation
tracking loss
point-evaluation policy
transposition policy
```

Semantic validation checks both the graph and required policies.

Compiler request resolution selects the point-sensor target and typed transposition
selection.

Lowerability checks that the registered FE point-evaluation and very-weak transpose
realization match that selection.

The observation does not become executable merely because the enum
`ObservationKind::point_sensor` exists.

### 27.3 PDAS for a cellwise box

Mathematics:

```math
\ell_{i}
\leq
u_{i}
\leq
r_{i}.
```

Semantic structure:

```text
cellwise_box constraint
lower/upper bound data
```

Lowerability establishes the concrete control/metric/bound realization.

Formulation capability additionally requires the registered PDAS-compatible scalar
target, metric conversion, rank/kernel assumptions, and inner KKT policy.

The box itself is only one ingredient.

## 28. Following the boundary through the source

For semantic diagnostics, start with:

- [`include/nmopt/semantic/v1/validation.hpp`](../../../include/nmopt/semantic/v1/validation.hpp)

The entry point is `SemanticValidator::validate()`.

Its call sequence is useful because it reveals the validator's scope: identity and
labels first, then regions/spaces/pairings, variables/data/transformations, equations
and terms, observations/losses, metrics/constraints, formulation/OTD declarations,
and finally policies.

For semantic resolution, read:

- [`include/nmopt/semantic/v1/resolved_problem.hpp`](../../../include/nmopt/semantic/v1/resolved_problem.hpp)

`ResolvedProblemView` and `SemanticResolver` are both small enough to read in one
sitting.

For compiler-side request closure and diagnostics, use:

- [`include/nmopt/compiler/v1/dealii_types.hpp`](../../../include/nmopt/compiler/v1/dealii_types.hpp)
- [`include/nmopt/compiler/v1/dealii_capabilities.hpp`](../../../include/nmopt/compiler/v1/dealii_capabilities.hpp)
- [`include/nmopt/compiler/v1/dealii_compiler.hpp`](../../../include/nmopt/compiler/v1/dealii_compiler.hpp)

`ResolvedCompilationRequest` is the key type to keep in mind.

For the bounded scalar planning path, read:

- [`include/nmopt/compiler/v1/dealii_scalar_plan.hpp`](../../../include/nmopt/compiler/v1/dealii_scalar_plan.hpp)

The distinction between `ScalarLoweringPlan`,
`ScalarResidualAssemblyPlan`, and `ScalarServicePlan` is the starting point for
Chapter 11.

For executable examples of semantic failures, see:

- [`tests/semantic/semantic_v1_contract.cc`](../../../tests/semantic/semantic_v1_contract.cc)

The point-sensor tests are particularly useful because they distinguish structural
errors, missing analytical policies, and typed transposition requirements in one
problem family.

## 29. Continue with Chapter 11

At the end of this chapter, the request has crossed the main acceptance boundaries:

```text
ProblemSpec
    validated
        │
        ▼
ResolvedProblemView
    indexed
        │
        ▼
ResolvedCompilationRequest
    closed
        │
        ▼
registered capability checks
    accepted
```

The next question is no longer

> Is this request valid and supported?

It is

> How does the compiler turn that accepted request into executable numerical
> services?

Chapter 11, **Compilation and lowering**, follows that process through component
handlers, `ScalarLoweringPlan`, data placements, registered target strategies,
runtime bindings, and `CompiledProblemT`.

## Read later

Useful existing documents are:

- [Describing and compiling a problem](../overview/semantic-compiler.md), for the
  short architecture-level description of the path.
- [Compiler reference](../../reference/compiler.md), for current public
  capabilities, products, diagnostics, and rejection boundaries.
- [Compiler implementation](../../internals/compiler.md), for request closure,
  lowering/planning mechanics, realization ownership, and verification.
- [Semantic problem model](09-semantic-problem-model.md), for the `ProblemSpec` vocabulary
  assumed by this chapter.
- [Project architecture](../overview/project-architecture.md), for the relationship
  between compiler-produced and externally supplied numerical services.

The compiler reference records the supported public paths; validation of a
concrete request determines whether the current compiler accepts it.
