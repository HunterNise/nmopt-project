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

Planned. This chapter will derive the reduced formulation, compare direct state
sensitivity with the adjoint elimination, and then explain state/adjoint solve
services, retained value evaluations, and derivative augmentation.

### 6. Reduced optimization methods

Planned. This chapter will follow the actual optimization loop: metric gradients,
steepest descent, L-BFGS/BFGS, line search, reduced Hessian actions, trust regions,
stopping criteria, and work accounting.

### 7. All-at-once, KKT, and active-set formulations

Planned. This chapter will develop supplied OTD, quadratic KKT systems,
complementarity, and PDAS as related but distinct formulation products.

## Part III – Constructing executable problems

### 8. Semantic problem model

Planned. A concrete problem will be rebuilt progressively as regions, spaces,
variables, data, residual terms, observations, losses, transformations, metrics,
constraints, requirements, and a formulation declaration.

### 9. Compilation and lowering

Planned. One semantic problem will be followed through validation, resolution,
planning, registered lowering, deal.II realization, and compiled-product packaging.

### 10. Finite-element realization

Planned. This chapter will examine the concrete deal.II constructions behind
coordinates, liftings, observations, metrics, constraints, and PDE solve services.

## Part IV – Integration and evidence

### 11. Native integration and ownership

Planned. This chapter will explain how an independently owned PDE application can
provide operations to nmopt without surrendering its native mesh, matrices, solves,
or output model.

### 12. Verification and evidence

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
