# Concept chapters

The concept chapters are the explanatory core of the manual.

They sit between the project overviews and the exact design/reference material. Their
job is to build enough mathematical and software understanding that a reader can
follow the implementation rather than merely recognize its names.

The chapters are organized as a progression. They are not intended to mirror the
directory tree or enumerate public classes.

## Part I – The numerical language of nmopt

### 1. Anatomy of a discrete PDE-constrained problem

[Anatomy of a discrete PDE-constrained problem](discrete-problem-anatomy.md) starts
from a concrete distributed elliptic control problem and follows it from the strong
PDE to a weak residual, finite-element spaces, coefficient vectors, matrices, state
solves, objective evaluation, and the adjoint calculation.

This is the recommended first concept chapter. It explains what sort of numerical
object the framework is ultimately trying to construct before introducing the
project-specific runtime representation.

### 2. Spaces, coordinates, and duality

[Spaces, coordinates, and duality](spaces-coordinates-and-duality.md) begins where the
first chapter leaves off: several coefficient vectors may all be arrays of numbers,
yet represent different spaces, different coordinate systems, and different
mathematical roles.

It develops physical versus independent coordinates, affine state reconstruction,
primal and dual coordinate representations, pairings, product spaces, and finally
`BlockLayout`, `PrimalBlockT`, and `CovectorBlockT` as the project's deliberately
small runtime representation of those distinctions.

### 3. Operators, derivatives, and adjoints

[Operators, derivatives, and adjoints](operators-derivatives-and-adjoints.md) develops
the residual as a nonlinear map, its directional derivative, Jacobian-vector and
transpose-Jacobian actions, and the pairing identity that defines an adjoint action.

It then connects those maps to `ExecutableModelT`, callback adapters, derivative
tests, and the Step-4 implementation without assuming that an explicit Jacobian
matrix must exist.

### 4. Metrics, gradients, and constraints

[Metrics, gradients, and constraints](metrics-gradients-and-constraints.md) starts
from a reduced derivative in $U_{h}^{\ast}$ and derives the Riesz representation that
turns it into a primal gradient.

It then compares several concrete finite-element metric realizations, explains why
regularization and optimization geometry are distinct choices, and develops
metric-dependent projection for box constraints before mapping the ideas to
`MetricT`, `ConstraintT`, and the reduced search code.

## Part II – Formulations and optimization algorithms

### 5. Reduced state–adjoint formulation

[Reduced state–adjoint formulation](reduced-state-adjoint-formulation.md) derives the
control-to-state reduction, compares direct sensitivities with the adjoint
elimination, and then follows the same construction through state/adjoint solve
services, retained value evaluations, and derivative augmentation.

### 6. Reduced optimization methods

[Reduced optimization methods](reduced-optimization-methods.md) follows the reduced
optimization loop after $j_h(u)$ and $j_h'(u)$ are available. It develops metric
steepest descent, nonlinear conjugate gradient, full and limited-memory BFGS,
Newton–CG, line-search globalization, projected steepest descent, trust regions,
stopping criteria, and numerical work/evidence records.

### 7. Optimality systems and KKT

[Optimality systems and KKT](optimality-systems-and-kkt.md) returns to the Lagrangian
without eliminating the state, derives the coupled state/adjoint/stationarity system,
distinguishes supplied OTD from DTO-derived formulations, and then develops the
equality-constrained quadratic KKT product, its pairings, assumptions, transpose
structure, and solver compatibility.

### 8. Complementarity and PDAS

[Complementarity and PDAS](complementarity-and-pdas.md) starts from the
box-constrained variational inequality, introduces the signed box multiplier and its
primal representation, derives lower/inactive/upper classification, and then follows
PDAS through active-coordinate restriction, shifted KKT subproblems, multiplier
recovery, residual checks, and active-set convergence.

## Part III – From problem descriptions to executable models

### 9. Semantic problem model

[Semantic problem model](semantic-problem-model.md) shows how a mathematical
PDE-constrained problem is translated into the semantic vocabulary accepted by
`nmopt`. Each major node is introduced through its mathematical role, an annotated
schema, a table of current semantic options, and—where useful—a small commented C++
fragment. The chapter also develops observation/loss formulas, supplied-OTD block
declarations, recipes, scenarios, and the boundary to later validation and numerical
realization.

### 10. Validation, resolution, and capabilities

[Validation, resolution, and capabilities](validation-resolution-and-capabilities.md)
separates structural and analytical-policy validation from semantic resolution,
compiler request closure, lowerability, and formulation/product capability. It
explains `ValidationReport`, `SemanticResolver`, `ResolvedProblemView`,
`ResolvedCompilationRequest`, bounded capability registration, and why a valid
semantic graph need not be executable by the current deal.II compiler.

### 11. Compilation and lowering

[Compilation and lowering](compilation-and-lowering.md) follows an accepted scalar
distributed-control problem through runtime data ports, `ScalarLoweringPlan`,
residual/service projections, finite-element model construction, coordinate
reconstruction and pullback, observation transposes, metric and state/adjoint solve
services, common executable contracts, `CompiledProblemT`, formulation-product
packaging, and the typed compilation manifest.

## Part IV – Using and integrating nmopt

### 12. Authoring and using compiled problems

[Authoring and using compiled problems](authoring-and-using-compiled-problems.md)
follows both the ready-made and authoring paths through recipes, scenarios, runtime
data, compilation sessions, compiler policies, compiled products, reduced solvers,
native output, and manifests. It also separates changes to semantic structure from
runtime data, numerical realization, solver policy, and experiment configuration.

### 13. Integrating an existing PDE application

[Integrating an existing PDE application](integrating-an-existing-pde-application.md)
starts from a numerical application that already owns its mesh, finite-element
operators, solves, and output. The Step-4 Problem B integration shows how layouts,
callback model actions, state/adjoint solve adapters, a native metric, ownership, and
lifetimes connect that application directly to the same reduced formulation and
optimizer contracts used by the compiler path.

## How the chapters should read

A chapter should normally introduce project-specific machinery only after the reader
has encountered the problem that machinery solves.

For example, the manual should not begin a discussion of `BlockLayout` by assuming
that “layouts” are self-evident mathematical objects. It should first establish that
several coefficient vectors with identical storage shapes may represent distinct
discrete spaces, and then explain why nmopt retains some of that identity at runtime.

Likewise, JVP/VJP terminology should follow a derivation of the linearized residual
and its transpose; `MetricT` should follow the derivative-versus-gradient problem;
and compiler planning types should follow the need to turn a semantic description
into one of a bounded set of supported numerical realizations.

The manual may revisit the same example in several chapters. The repetition is
intentional when a later chapter reveals another layer of the same object.
