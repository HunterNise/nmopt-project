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

Planned. A concrete problem will be rebuilt progressively as regions, spaces,
variables, data, residual terms, observations, losses, transformations, metrics,
constraints, requirements, and a formulation declaration.

### 10. Validation, resolution, and capabilities

Planned. This chapter will distinguish semantic validity, resolution, lowerability,
compiler capability, and formulation/product capability before numerical realization.

### 11. Compilation and lowering

Planned. One resolved problem will be followed through component planning,
`ScalarLoweringPlan`, runtime data bindings, registered realization strategies, and
compiled-product packaging.

### 12. Finite-element realization

Planned. This chapter will examine the concrete deal.II constructions behind
coordinates, liftings, residuals, observations, metrics, constraints, and PDE solve
services.

## Part IV – Integration and evidence

### 13. Native integration and ownership

Planned. This chapter will explain how an independently owned PDE application can
provide operations to nmopt without surrendering its native mesh, matrices, solves,
or output model.

### 14. Verification and evidence

Planned. This chapter will connect derivative identities, algebraic contract tests,
compiler integration tests, solve reports, benchmark evidence, and reproduction
artifacts.

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
