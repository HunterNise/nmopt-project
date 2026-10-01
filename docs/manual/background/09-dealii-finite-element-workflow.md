# 09 · The deal.II finite-element workflow

**Background navigation:** [Index](README.md) \
Previous: [08 · Modern C++ for scientific code](08-modern-cpp-for-scientific-code.md) \
Next: [manual overview · Numerical realization](../overview/numerical-realization.md)

## Purpose

[03 · Finite elements](03-finite-elements.md) developed finite elements as mathematics: a mesh, local polynomial spaces, degrees of freedom, basis functions, quadrature, local element contributions, global assembly, constraints, and a linear system. This chapter answers a different question:

> How are those concepts represented and connected in deal.II?

The goal is not to reproduce the deal.II tutorial or catalogue the library API. The current deal.II documentation already does that more comprehensively and stays synchronized with the library. Instead, this chapter builds a **conceptual map** between the finite-element language used earlier and the small set of deal.II objects that recur throughout `nmopt`'s finite-element realization.

The central workflow is

```text
Triangulation
    ↓
FiniteElement + Mapping + Quadrature
    ↓
DoFHandler
    ↓
constraints + sparsity + vectors/matrices
    ↓
cell/face loops with FEValues / FEFaceValues
    ↓
global algebraic system
    ↓
linear solver
    ↓
finite-element coefficient vector
    ↓
DataOut / application output
```

A compact Poisson problem will run through the whole chapter. It is deliberately close to deal.II's early tutorial programs because those tutorials are the best next source after this orientation.

The chapter does **not** teach adaptive refinement, parallel distributed meshes, matrix-free operator evaluation, mixed finite elements, hp methods, multigrid internals, or the full deal.II type hierarchy. Those can be learned when a concrete application requires them.

## Before you start

You should be comfortable with:

- [03 · Finite elements](03-finite-elements.md), especially basis functions, degrees of freedom, assembly, quadrature, and constrained coefficients;
- [04 · Numerical linear algebra for PDE and optimization systems](04-numerical-linear-algebra.md), especially sparse systems, CG, and preconditioning;
- [08 · Modern C++ for scientific code](08-modern-cpp-for-scientific-code.md), especially templates, references, object lifetime, and class interfaces;
- [08a · Scientific software workflow: shell, compilation, CMake, and tests](08a-scientific-software-workflow.md) if CMake/build terminology is unfamiliar.

This chapter uses deal.II class and method names as code symbols. Their exact signatures can evolve, so the current deal.II documentation is authoritative when adapting an example to an installed version.

## What you will be able to do

After this chapter, you should be able to:

- explain the role of `Triangulation<dim>` in a finite-element program;
- distinguish a mesh from a finite-element space and from its global degree-of-freedom numbering;
- explain what `FE_Q<dim>` describes;
- explain what `DoFHandler<dim>` adds to a triangulation and finite element;
- identify where reference-to-physical-cell mappings enter computation;
- choose and interpret a quadrature object at a basic level;
- explain why `FEValues<dim>` combines a mapping, finite element, quadrature rule, and requested update data;
- read a standard cell assembly loop and connect every index to the formulas of [03 · Finite elements](03-finite-elements.md);
- understand `JxW(q)` as the mapped quadrature weight;
- distinguish `FEValues` from `FEFaceValues`;
- explain why sparse matrix structure can be built before matrix entries are assembled;
- understand the role of `AffineConstraints<double>` in eliminating or relating constrained degrees of freedom;
- recognize the roles of deal.II vectors, sparse matrices, `SolverControl`, iterative solvers, and preconditioners;
- understand how `Function<dim>`-like objects supply spatial data;
- understand the purpose of `DataOut<dim>`;
- follow the dependencies among the main finite-element objects and avoid obvious lifetime mistakes;
- trace a Poisson problem from weak form to a deal.II program skeleton;
- see where `nmopt` adds semantic, formulation, and optimization layers above the finite-element workflow.

## Roadmap

We first map the mathematical objects from [03 · Finite elements](03-finite-elements.md) to deal.II classes without assembling anything. We then focus on the core operation of a Galerkin code: reinitializing finite-element data on one cell, evaluating shape functions and gradients at quadrature points, computing a local matrix/vector, and scattering those local contributions into global algebraic objects.

After assembly, we treat constraints, sparse linear algebra, solution of the discrete system, and output. The final sections walk through the complete Poisson pipeline and explain how `nmopt` uses deal.II as a numerical realization beneath more abstract optimization interfaces.

## 1. One finite-element program contains several distinct structures

Return to the weak Poisson problem

$$
\text{find }y\in H_{0}^{1}(\Omega)
$$

such that

$$
\int_{\Omega}\nabla y\cdot\nabla v\thinspace\mathrm{d}x
= \int_{\Omega}fv\thinspace\mathrm{d}x
\qquad
\forall v\in H_{0}^{1}(\Omega).
$$

A conforming finite-element discretization chooses a space

$$
V_{h}
= \mathop{\mathrm{span}}\lbrace\varphi_{1},\ldots,\varphi_{N}\rbrace
$$

and seeks

$$
y_{h}
= \sum_{j=1}^{N}Y_{j}\varphi_{j}.
$$

The discrete system is

$$
AY=F,
$$

with

$$
A_{ij}
= \int_{\Omega}\nabla\varphi_{j}\cdot\nabla\varphi_{i}\thinspace\mathrm{d}x,
\qquad
F_{i}
= \int_{\Omega}f\varphi_{i}\thinspace\mathrm{d}x.
$$

A program must represent at least four conceptually different things:

```text
geometry
    cells and their connectivity

local finite element
    shape functions and local DoFs on a reference cell

global finite-element space
    one global numbering that identifies shared DoFs

algebraic data
    coefficient vectors and sparse matrix entries
```

A common source of confusion is to call all of these “the mesh.” deal.II keeps them in different objects precisely because they carry different information.

## 2. `Triangulation<dim>` represents the computational mesh

A `Triangulation<dim>` stores the cells and mesh topology. In a simple program one might write

```cpp
dealii::Triangulation<dim> triangulation;
dealii::GridGenerator::hyper_cube(triangulation, 0.0, 1.0);
triangulation.refine_global(3);
```

The mathematical counterpart is the triangulation

$$
\mathcal T_{h}
= \lbrace K\rbrace.
$$

The triangulation knows which cells exist and how they meet. It does not by itself say whether the approximation uses piecewise linear, quadratic, continuous, discontinuous, scalar, or vector-valued finite elements.

### 2.1 Cell iterators let code traverse the mesh

A typical loop has the form

```cpp
for (const auto &cell : dof_handler.active_cell_iterators())
  {
    // work on one active cell
  }
```

An **active cell** is a cell not currently refined into active children. On a uniformly unrefined mesh this distinction is nearly invisible; it becomes essential with adaptive refinement.

The cell iterator provides geometric/topological information and, when obtained through a `DoFHandler`, access to finite-element degree-of-freedom indices associated with that cell.

## 3. A `FiniteElement` object describes the local approximation

For a standard continuous scalar Lagrange element of polynomial degree $r$, deal.II commonly uses

```cpp
dealii::FE_Q<dim> fe(r);
```

For example,

```cpp
dealii::FE_Q<dim> fe(1);
```

represents a continuous first-degree scalar element on the reference-cell family supported by `FE_Q`.

The important conceptual role is:

```text
FiniteElement object
    local shape functions
    local DoF placement/meaning
    transformation rules
    number of local DoFs
```

It does **not** enumerate all global basis functions over the mesh. That is the job of the `DoFHandler` once the local element has been associated with the triangulation.

### 3.1 `fe.n_dofs_per_cell()` is the local algebraic dimension

For one cell $K$, let

$$
\lbrace\varphi_{\alpha}^{K}\rbrace_{\alpha=1}^{n_{\text{loc}}}
$$

be its local basis. Then the code quantity corresponding to $n_{\text{loc}}$ is conceptually

```cpp
const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
```

This determines the size of a cell matrix:

```cpp
dealii::FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
```

and cell right-hand side:

```cpp
dealii::Vector<double> cell_rhs(dofs_per_cell);
```

The Greek local indices $\alpha,\beta$ from [03 · Finite elements](03-finite-elements.md) become integer loop indices such as `i` and `j` in assembly code.

## 4. `DoFHandler<dim>` turns local finite elements into a globally numbered space

The mathematical finite-element space is not just the direct collection of independent cell polynomials. Continuous finite elements identify local degrees of freedom that represent the same global entity.

A `DoFHandler<dim>` combines a triangulation with a finite-element description and enumerates the resulting global degrees of freedom.

A common setup is

```cpp
dealii::DoFHandler<dim> dof_handler(triangulation);
dof_handler.distribute_dofs(fe);
```

After distribution,

```cpp
dof_handler.n_dofs()
```

is the global dimension $N$ of the discrete space.

The deal.II documentation describes `distribute_dofs()` as enumerating the basis functions/DoFs needed by the finite element over the triangulation. The exact global numbering order is an implementation detail unless explicitly renumbered.

### 4.1 Local indices must be mapped to global indices

On one cell, local basis function number `i` corresponds to some global basis-function number $I_{K}(i)$.

The code obtains these indices with a container such as

```cpp
std::vector<dealii::types::global_dof_index>
  local_dof_indices(dofs_per_cell);

cell->get_dof_indices(local_dof_indices);
```

Then

```cpp
local_dof_indices[i]
```

is the global index associated with local basis function `i` on that cell.

This is exactly the local-to-global map used in the assembly formulas of [03 · Finite elements](03-finite-elements.md).

### 4.2 Mesh vertices are not the same thing as DoFs

For $Q_{1}$ scalar elements on a simple conforming mesh, one can easily get the impression that “there is one DoF per vertex.” That intuition fails as soon as one uses:

- higher polynomial degree;
- vector-valued elements;
- discontinuous elements;
- edge/face/interior DoFs;
- mixed finite elements.

The robust mental model is:

```text
Triangulation
    geometric/topological entities

FiniteElement
    says which local DoFs live on those entities

DoFHandler
    gives the resulting global enumeration
```

This is why `DoFHandler` deserves to be a separate object rather than an integer array attached casually to the mesh.

## 5. The mapping connects the reference cell to a physical cell

Finite-element shape functions and quadrature rules are naturally described on a reference cell. Integrals, gradients, and geometry live on the physical cell $K$.

[03 · Finite elements](03-finite-elements.md) introduced a map

$$
F_{K}:\widehat K\to K.
$$

In deal.II, a `Mapping` object represents this geometric transformation machinery.

For simple straight-sided meshes, a common explicit choice is

```cpp
dealii::MappingQ1<dim> mapping;
```

More general mappings are required when curved geometry or higher-order geometric approximation matters.

### 5.1 Why the mapping matters even for a simple Laplace problem

The reference gradient of a shape function is not in general the physical gradient used in

$$
\int_{K}\nabla\varphi_{i}\cdot\nabla\varphi_{j}\thinspace\mathrm{d}x.
$$

The gradient must be transformed by the cell map, and the integration measure picks up the Jacobian determinant.

These transformations are routine but error-prone to implement manually. `FEValues` coordinates the finite element, mapping, and quadrature so assembly code can request physical-cell values and gradients directly.

## 6. A quadrature object represents numerical integration points and weights

For cell integrals a common rule is

```cpp
dealii::QGauss<dim> quadrature(fe.degree + 1);
```

`QGauss<dim>` is the tensor-product Gauss rule in `dim` dimensions for tensor-product cells. The integer argument gives the number of points in each one-dimensional direction.

The mathematical approximation is

$$
\int_{K}g(x)\thinspace\mathrm{d}x
\approx
\sum_{q=1}^{n_{q}}
w_{q}^{K}g(x_{q}^{K}).
$$

### 6.1 Quadrature degree is part of the numerical method

Choosing a finite-element space does not automatically fix the integration rule. If the integrand is polynomial and the rule is sufficiently exact, quadrature reproduces the intended Galerkin integral. With variable coefficients, curved mappings, nonlinear terms, or nonpolynomial data, exactness arguments change.

The purpose of this chapter is to show where quadrature appears in code, not to repeat its analysis. See [03a · Further finite-element notes](03a-further-finite-element-notes.md) for quadrature formulas, exactness, mapped integration, and variational crimes.

### 6.2 Cell and face quadrature have different dimensions

A volume integral over $K$ uses a `Quadrature<dim>` object. A boundary/face integral uses points on a $(d-1)$-dimensional reference face, commonly represented by something like

```cpp
dealii::QGauss<dim - 1> face_quadrature(...);
```

This distinction later appears in `FEValues` versus `FEFaceValues`.

## 7. `FEValues<dim>` is the bridge from symbolic FE formulas to cell data

A typical construction is

```cpp
dealii::FEValues<dim> fe_values(
  mapping,
  fe,
  quadrature,
  dealii::update_values |
  dealii::update_gradients |
  dealii::update_quadrature_points |
  dealii::update_JxW_values);
```

The object is told four things:

```text
mapping
    how reference and physical cells are related

finite element
    which shape functions are used

quadrature
    where integrals are sampled

update flags
    which derived data will actually be needed
```

### 7.1 Update flags are a computational request

If assembly only needs shape gradients and mapped quadrature weights, asking for Hessians or physical quadrature points would compute and store unnecessary information.

For the Poisson stiffness matrix, the essential requests are

```cpp
dealii::update_gradients |
dealii::update_JxW_values
```

If the right-hand side $f(x)$ must be evaluated at physical quadrature points, add

```cpp
dealii::update_values |
dealii::update_quadrature_points
```

because we also need $\varphi_{i}(x_{q})$ and $x_{q}$.

### 7.2 `reinit(cell)` specializes the object to one physical cell

Before accessing cell-dependent data:

```cpp
fe_values.reinit(cell);
```

This causes `FEValues` to compute or expose the requested mapped information for that particular cell.

The pattern is therefore

```text
construct FEValues once
        ↓
for each cell
    reinit(cell)
    use data for this cell
```

not “construct a brand-new finite element for every quadrature point.”

## 8. `shape_value`, `shape_grad`, and `JxW` implement the quadrature formula

For a local shape function `i` and quadrature point `q`, assembly code can use

```cpp
fe_values.shape_value(i, q)
```

for

$$
\varphi_{i}^{K}(x_{q}^{K}),
$$

and

```cpp
fe_values.shape_grad(i, q)
```

for the physical gradient

$$
\nabla\varphi_{i}^{K}(x_{q}^{K}).
$$

The quantity

```cpp
fe_values.JxW(q)
```

is the quadrature weight transformed to the physical cell. Conceptually,

$$
\mathrm{JxW}_{q}
= \left\lvert\det DF_{K}(\widehat x_{q})\right\rvert\widehat w_{q}
$$

for an ordinary full-dimensional cell map.

Thus a cell contribution such as

$$
A_{ij}^{K}
= \int_{K}\nabla\varphi_{j}^{K}\cdot\nabla\varphi_{i}^{K}\thinspace\mathrm{d}x
$$

becomes

```cpp
cell_matrix(i, j) +=
  fe_values.shape_grad(i, q) *
  fe_values.shape_grad(j, q) *
  fe_values.JxW(q);
```

The code is now almost a transcription of the variational formula.

## 9. A cell assembly loop mirrors the mathematics

For the Poisson problem, a complete local assembly core looks schematically like

```cpp
for (const auto &cell : dof_handler.active_cell_iterators())
  {
    fe_values.reinit(cell);

    cell_matrix = 0;
    cell_rhs = 0;

    for (unsigned int q = 0; q < quadrature.size(); ++q)
      {
        const auto x_q = fe_values.quadrature_point(q);
        const double f_q = rhs.value(x_q);

        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              cell_matrix(i, j) +=
                fe_values.shape_grad(i, q) *
                fe_values.shape_grad(j, q) *
                fe_values.JxW(q);

            cell_rhs(i) +=
              f_q *
              fe_values.shape_value(i, q) *
              fe_values.JxW(q);
          }
      }

    cell->get_dof_indices(local_dof_indices);

    constraints.distribute_local_to_global(
      cell_matrix,
      cell_rhs,
      local_dof_indices,
      system_matrix,
      system_rhs);
  }
```

The last call performs the local-to-global scatter while respecting the affine constraints already declared for the discrete space. In an unconstrained formulation, the same operation could be written as explicit additions into the global matrix and vector; using the constraint object keeps the algebraic coordinate reduction in one place.

Every nested loop has a mathematical meaning:

```text
cell loop
    sum integrals over K in the mesh

quadrature loop q
    approximate each cell integral

local i loop
    test-function index

local j loop
    trial/basis-function index
```

The local matrix is dense even though the global matrix is sparse, because all local basis functions on one cell may interact.

## 10. Global sparsity comes from local support

The global stiffness matrix has

$$
A_{IJ}\neq0
$$

only when the corresponding global basis functions overlap on at least one cell, aside from constraint-related modifications.

Therefore the sparsity graph can be determined from the mesh, DoF layout, and constraints before computing numerical matrix entries.

A common pattern is conceptually

```cpp
dealii::DynamicSparsityPattern dsp(dof_handler.n_dofs());
dealii::DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);

sparsity_pattern.copy_from(dsp);
system_matrix.reinit(sparsity_pattern);
```

The dynamic pattern is convenient while entries are being discovered. The compressed `SparsityPattern` is then used by a `SparseMatrix<double>`.

### 10.1 Sparsity pattern and matrix values are separate objects

This separation is useful because the nonzero *positions* often remain fixed while the numerical entries change, for example across nonlinear iterations or parameter updates on a fixed discretization.

A sparse matrix is therefore not just a dense matrix with many zeros omitted after assembly. Its allowed nonzero structure is part of its storage design.

## 11. Scattering maps local element arrays into global arrays

Without constraints, the conceptual assembly is

$$
A_{I_{K}(i),I_{K}(j)}
\mathrel{+}=
A_{ij}^{K},
$$

$$
F_{I_{K}(i)}
\mathrel{+}=
F_{i}^{K}.
$$

One could code those index additions explicitly. With constrained DoFs, however, the local-to-global map must also account for elimination/redistribution of constrained contributions.

This motivates using `AffineConstraints<double>` as part of assembly rather than treating boundary or hanging-node corrections as an unrelated postprocessing step.

## 12. `AffineConstraints<double>` represents linear relations among DoFs

A constrained coefficient may satisfy an affine relation

$$
Y_{i}
= \sum_{j}c_{ij}Y_{j}+b_{i}.
$$

This general form covers several important cases:

- homogeneous Dirichlet value: $Y_{i}=0$;
- inhomogeneous prescribed value: $Y_{i}=b_{i}$;
- hanging-node relation: one constrained DoF is a linear combination of neighboring unconstrained DoFs.

A typical object is

```cpp
dealii::AffineConstraints<double> constraints;
```

The workflow is conceptually:

```text
create constraint object
    ↓
add hanging-node relations and/or boundary conditions
    ↓
close constraints
    ↓
use constraints while building sparsity and assembling
    ↓
solve for algebraic unknowns
    ↓
distribute constrained values into final coefficient vector
```

### 12.1 Constraint-aware assembly keeps independent and physical coefficients distinct

[03 · Finite elements](03-finite-elements.md) wrote a physical coefficient vector schematically as

$$
y_{\mathrm{phys}}
= Pz+\ell,
$$

where $z$ contains independent coordinates and $\ell$ carries fixed values.

`AffineConstraints` expresses the same kind of algebraic relationship operationally. It allows local cell contributions to be condensed into the independent system while preserving enough information to reconstruct the full finite-element vector afterward.

This is why constraints are mathematically relevant to derivatives and adjoints as well as to boundary-value bookkeeping.

### 12.2 Hanging nodes are the canonical reason constraints become more than fixed boundary values

On a locally refined conforming mesh, a fine-grid node may lie on the interior of a coarse-grid face. Its value is not an independent degree of freedom if global continuity is to be preserved. The required interpolation relation is naturally linear.

The deal.II step-6 tutorial is the early tutorial that introduces adaptive refinement and hanging-node handling. For the current `nmopt` background, understanding the *existence and role* of such affine relations is enough; adaptivity itself remains out of scope.

## 13. `FEFaceValues` provides the analogous machinery on faces

Boundary integrals require values on a cell face rather than over the cell volume. A common construction is

```cpp
dealii::FEFaceValues<dim> fe_face_values(
  mapping,
  fe,
  face_quadrature,
  dealii::update_values |
  dealii::update_quadrature_points |
  dealii::update_normal_vectors |
  dealii::update_JxW_values);
```

The reinitialization now specifies both a cell and one of its faces:

```cpp
fe_face_values.reinit(cell, face_no);
```

This supports terms such as

$$
\int_{\Gamma_{N}}g_{N}v\thinspace\mathrm{d}s
$$

or boundary-control couplings.

### 13.1 Normals belong to face geometry

If an integrand contains

$$
n(x)
$$

or a flux involving the outward normal, request

```cpp
dealii::update_normal_vectors
```

and access the normal at a quadrature point through the face-values object.

This is one place where sign conventions from the weak form become concrete implementation choices. The mathematical derivation should determine which normal or conormal quantity the code must assemble; the library only provides the geometric normal.

## 14. `Function<dim>`-like objects supply spatial data

A right-hand side $f(x)$, boundary datum, target state, or coefficient often needs to be evaluated at quadrature points.

A deal.II-style interface is based on objects that can answer queries such as

```cpp
rhs.value(point)
```

for a physical `Point<dim>`.

Library classes such as `Function<dim>` and derived/helper function classes provide this abstraction.

The useful separation is:

```text
assembly code
    asks for value at physical point x_q

data object
    decides how that value is represented or computed
```

The assembly loop therefore does not need to know whether $f$ is a formula, constant, parsed expression, interpolated dataset, or application-specific callback, provided the appropriate interface is supplied.

### 14.1 Data evaluation and finite-element coefficients are not the same operation

A prescribed smooth function may be evaluated directly at quadrature points without first being represented as a finite-element coefficient vector.

Conversely, if a control or state is itself a discrete finite-element field, its values at quadrature points are obtained from its coefficient vector and shape functions.

This is the implementation counterpart of the distinction made in [03 · Finite elements](03-finite-elements.md) among interpolation, projection, cell averages, and direct quadrature-point evaluation.

## 15. Vectors and matrices store coefficient-space algebra

For a serial problem, simple deal.II programs often use types conceptually like

```cpp
dealii::Vector<double> solution;
dealii::Vector<double> system_rhs;
dealii::SparseMatrix<double> system_matrix;
```

After DoF distribution,

```cpp
solution.reinit(dof_handler.n_dofs());
system_rhs.reinit(dof_handler.n_dofs());
```

sets the vector dimensions.

The coefficient vector

```text
solution
```

is not the finite-element function itself in an abstract mathematical sense. It is the coordinate representation from which the finite-element field is reconstructed using the global basis associated with the `DoFHandler`.

That distinction is essential in `nmopt`, where physical fields, independent coordinates, primal vectors, and covectors may deliberately have different representations.

## 16. `SolverControl` describes when an iterative solver should stop

A solver needs a policy for acceptable convergence and a limit on work. A simple configuration might be

```cpp
dealii::SolverControl solver_control(1000, 1e-12);
```

followed by

```cpp
dealii::SolverCG<dealii::Vector<double>> solver(solver_control);
```

The two roles are different:

```text
SolverControl
    convergence / stopping policy and report state

SolverCG
    iterative algorithm generating candidate solutions
```

The exact constructors and available control classes depend on deal.II version and desired criterion, so consult the current API for production code.

### 16.1 Solver choice follows matrix structure

For the symmetric positive-definite Poisson stiffness system, CG is a natural choice. That is a mathematical statement from [04 · Numerical linear algebra for PDE and optimization systems](04-numerical-linear-algebra.md), not a deal.II-specific fact.

For nonsymmetric or indefinite systems, different algorithms are required. The library supplies many solver classes; choosing one should follow the operator structure rather than whichever class name is easiest to instantiate.

### 16.2 A preconditioner is passed as another object

Schematically:

```cpp
solver.solve(system_matrix,
             solution,
             system_rhs,
             preconditioner);
```

The solver algorithm and preconditioner remain separate roles. This reflects the linear-algebra distinction developed in [04 · Numerical linear algebra for PDE and optimization systems](04-numerical-linear-algebra.md).

For the simplest orientation one can use an identity preconditioner. Real PDE problems commonly need a stronger preconditioner for scalable iteration counts; [04a · Further numerical linear algebra for PDE-constrained optimization](04a-further-linear-algebra-for-pde-optimization.md) explains why.

## 17. After solving, distribute constrained values

If the algebraic system was assembled in a constraint-aware condensed form, the solution must be expanded back into the full coefficient representation:

```cpp
constraints.distribute(solution);
```

Conceptually this applies the affine reconstruction relationships, restoring fixed and dependent DoFs.

The resulting vector can then be interpreted as the physical finite-element coefficient vector for evaluation and output.

This step is easy to treat as implementation housekeeping. It is actually the operational form of the coordinate map between independent algebraic unknowns and the full finite-element field.

## 18. `DataOut<dim>` turns finite-element coefficients into visualization data

After solving, a common output sequence is

```cpp
dealii::DataOut<dim> data_out;
data_out.attach_dof_handler(dof_handler);
data_out.add_data_vector(solution, "solution");
data_out.build_patches();
```

and then a format-specific write operation, for example VTU output.

The `DoFHandler` is needed because a coefficient vector by itself does not say where basis functions live geometrically. Output needs the finite-element-space interpretation of the coefficients.

### 18.1 Output is not the numerical solution algorithm

Writing a `.vtu` file can be important for scientific inspection, but it should remain conceptually separate from:

- state assembly;
- linear solution;
- convergence checks;
- optimization metrics.

A pretty visualization is not evidence that a residual is small or an optimizer converged. Conversely, a mathematically correct solve may use no visualization at all.

## 19. Object dependencies and lifetime are part of the program design

A typical problem class might contain members in this conceptual order:

```cpp
dealii::Triangulation<dim> triangulation;
dealii::FE_Q<dim> fe;
dealii::DoFHandler<dim> dof_handler;
dealii::AffineConstraints<double> constraints;
dealii::SparsityPattern sparsity_pattern;
dealii::SparseMatrix<double> system_matrix;
dealii::Vector<double> solution;
dealii::Vector<double> system_rhs;
```

The objects are related:

```text
DoFHandler
    uses the triangulation and distributed FE description

sparsity pattern
    derived from DoF connectivity + constraints

matrix/vector dimensions
    derived from n_dofs()

FEValues during assembly
    combines mapping + FE + quadrature + current cell
```

These dependencies matter to lifetime. If an object holds or uses a reference to another library object, the referenced object must remain valid for the period required by the API.

### 19.1 Prefer learning the dependency graph over memorizing constructor order

A robust question is:

> Which mathematical object must already exist for this object to be meaningful?

For example:

- a `DoFHandler` needs a triangulation;
- global vectors cannot be correctly sized until the global DoF count is known;
- an assembly loop cannot evaluate physical shape gradients until finite element, mapping, quadrature, and current cell are available.

C++ declaration order in a class then follows those semantic dependencies.

### 19.2 Check current ownership semantics rather than inferring them from names

deal.II objects do not all store their constructor arguments in the same way. For example, current `DoFHandler::distribute_dofs()` documentation states that the finite element is copied into the handler, whereas the handler itself is associated with a triangulation.

When a lifetime detail becomes important, consult the class documentation instead of assuming that every `const T &` constructor parameter is stored as a reference forever or that every object owns a deep copy.

The mental model should be precise enough to trigger the right question, not replace the API contract.

## 20. A compact Poisson program has a recognizable phase structure

A deal.II tutorial-style problem class often separates the computation into functions such as

```text
make_grid
setup_system
assemble_system
solve
output_results
run
```

This decomposition is pedagogically useful because it mirrors mathematical stages.

### 20.1 `make_grid`: choose the discrete geometry

```cpp
void make_grid()
{
  dealii::GridGenerator::hyper_cube(triangulation, 0.0, 1.0);
  triangulation.refine_global(3);
}
```

The output of this phase is the mesh $\mathcal T_{h}$.

### 20.2 `setup_system`: build the finite-dimensional coordinate structure

Conceptually:

```cpp
dof_handler.distribute_dofs(fe);

constraints.clear();
// add boundary and hanging-node constraints as appropriate
constraints.close();

// build sparsity from dof_handler + constraints
// initialize matrix and vectors
```

At the end of this phase, the program knows the discrete-space dimension and which global matrix entries may occur.

No cell integrals need to have been evaluated yet.

### 20.3 `assemble_system`: evaluate the weak form

This is where the variational problem becomes matrix and vector entries.

For each cell:

```text
reinitialize FEValues
        ↓
evaluate local quadrature sum
        ↓
obtain local-to-global DoF indices
        ↓
constraint-aware scatter
```

The result is the global discrete equation.

### 20.4 `solve`: apply numerical linear algebra

The finite-element method has now produced

$$
AY=F.
$$

The solver stage chooses an algorithm appropriate to $A$, applies a stopping criterion, and records convergence information.

It is conceptually downstream of finite-element assembly even though advanced matrix-free methods may blur the distinction by applying the operator without explicitly assembling all entries.

### 20.5 `output_results`: interpret coordinates as a field again

The coefficient vector becomes meaningful for visualization only together with the DoFHandler/finite-element space. `DataOut` performs this reconstruction for an output representation.

### 20.6 `run`: organize the phases

A simple `run()` function makes the progression visible:

```cpp
void run()
{
  make_grid();
  setup_system();
  assemble_system();
  solve();
  output_results();
}
```

Real applications may repeat some phases, perform nonlinear iterations, refine the mesh, or solve many right-hand sides. The phase separation remains a useful way to reason about what must be recomputed when something changes.

## 21. What must be rebuilt when the problem changes?

The object map helps answer performance questions.

Suppose the mesh and finite-element space remain fixed but only the right-hand side $f$ changes. Then:

```text
mesh topology                    unchanged
DoF numbering                    unchanged
sparsity pattern                 unchanged
Poisson stiffness matrix         unchanged if coefficients unchanged
right-hand-side vector           changes
solution                         changes
```

Therefore rebuilding the mesh and sparsity structure for every new forcing would be unnecessary.

If the diffusion coefficient changes but the mesh does not, the sparsity pattern may remain fixed while matrix entries must be reassembled.

If the mesh changes, DoF numbering and sparsity usually need reconstruction as well.

This dependency reasoning becomes important in PDE-constrained optimization, where state and adjoint equations are solved repeatedly with related operators.

## 22. State and adjoint assembly often reuse the same finite-element infrastructure

In a linear distributed-control problem, the state equation might have matrix form

$$
AY=f+Bu.
$$

The adjoint equation may use

$$
A^{\mathsf T}P=r_{y}.
$$

For a symmetric Poisson operator,

$$
A^{\mathsf T}=A,
$$

so state and adjoint solves can reuse the same assembled stiffness matrix and often the same preconditioning infrastructure.

The right-hand sides and mathematical roles are still different.

For nonsymmetric PDEs, the transpose operator must be treated explicitly. A library implementation may assemble a transpose matrix, provide transpose application, or construct the adjoint weak form directly, depending on the formulation.

The deal.II object model does not decide this mathematical question. It supplies the finite-element mechanisms needed to realize the chosen state and adjoint operators.

## 23. Control spaces may use a different `DoFHandler`

PDE-constrained optimization often has

$$
Y_{h}\neq U_{h}.
$$

For example, the state may use continuous $Q_{1}$ elements while a distributed control uses cellwise constants or a discontinuous space.

Then one should expect separate objects such as conceptually

```text
state finite element
state DoFHandler

control finite element
control DoFHandler
```

and possibly rectangular coupling matrices.

This is the implementation form of the [03 · Finite elements](03-finite-elements.md) relation

$$
AY=f+BU,
$$

where $B$ need not be square.

A single physical mesh can support several different discrete spaces and therefore several independent DoF numberings.

## 24. Where `nmopt` sits above deal.II

A direct deal.II application can own all finite-element decisions itself:

```text
application code
    chooses mesh
    chooses FE spaces
    assembles operators
    solves systems
    writes output
```

`nmopt` introduces additional layers because it wants optimization formulations to work through explicit numerical contracts rather than depending directly on one application implementation.

A useful picture is

```text
structured semantic problem        existing deal.II application
          │                                   │
          ▼ compile                           ▼ adapt
compiler-owned deal.II realization  application-owned realization
          │                                   │
          └──────────────┬────────────────────┘
                         ▼
              common numerical services
        layouts, residual/JVP/VJP, solves,
              metrics, constraints
                         ▼
                    formulation
          reduced, KKT, complementarity, PDAS
                         ▼
                     algorithm
```

The deal.II classes in this chapter live primarily in the numerical-realization boxes. They realize finite-dimensional spaces and operations; they do not by themselves encode whether a vector is a state, control, adjoint, derivative, or optimization gradient.

### 24.1 Lowering preserves semantic roles while choosing deal.II objects

The `nmopt` compiler path can take a semantic declaration such as a scalar diffusion-reaction residual, an $L^{2}$ control metric, and a particular control realization, then construct the deal.II spaces and operators that implement those roles.

This is why the project documentation distinguishes semantic declarations from runtime mesh/data bindings and discretization policy. The same mathematical role may admit several numerical realizations.

The `nmopt` manual chapter [manual 11 · Compilation and lowering](../concepts/11-compilation-and-lowering.md) develops that pipeline. This chapter supplies the deal.II vocabulary needed to understand what the lowerer eventually creates.

## 25. A map from finite-element terminology to deal.II

The following compact map can be used while reading code:

| Finite-element concept | deal.II object or operation |
| --- | --- |
| triangulation $\mathcal T_{h}$ | `Triangulation<dim>` |
| local Lagrange element | `FE_Q<dim>` |
| reference-to-physical geometry map | `Mapping...<dim>` |
| global FE basis / DoF numbering | `DoFHandler<dim>` + `distribute_dofs()` |
| quadrature nodes and weights | `QGauss<dim>` or another `Quadrature<dim>` |
| values/gradients on one cell | `FEValues<dim>` after `reinit(cell)` |
| values/normals on one face | `FEFaceValues<dim>` |
| mapped weight $\lvert\det DF_{K}\rvert w_{q}$ | `JxW(q)` |
| local matrix/vector | `FullMatrix<double>`, `Vector<double>` |
| local-to-global DoF map | `cell->get_dof_indices(...)` |
| affine DoF relations | `AffineConstraints<double>` |
| global sparsity | `DynamicSparsityPattern`, `SparsityPattern` |
| assembled sparse operator | `SparseMatrix<double>` or another backend matrix |
| coefficient vector | deal.II/backend vector type |
| iterative stopping policy | `SolverControl` |
| CG algorithm | `SolverCG<...>` |
| prescribed spatial data | `Function<dim>`-style interface |
| visualization/output reconstruction | `DataOut<dim>` |

The table is an orientation map, not a promise that every application uses exactly these concrete classes. deal.II offers multiple backends, hp collections, matrix-free facilities, parallel vectors, different mappings, and many other alternatives.

## 26. A useful reading order for unfamiliar deal.II code

When opening a finite-element application, do not begin by decoding every template argument. A productive order is:

1. identify the main problem class and its member objects;
2. find the mesh creation function;
3. find where `distribute_dofs()` is called;
4. identify the finite element and quadrature rule;
5. find the `FEValues`/`FEFaceValues` construction and update flags;
6. inspect the local assembly formula;
7. see how local indices are scattered under constraints;
8. identify the solver and preconditioner;
9. find where constrained values are reconstructed;
10. inspect output only after the numerical path is understood.

This order follows the mathematical data flow rather than the order in which C++ declarations happen to appear in the file.

## 27. Scope frontier

The current project requires a much smaller deal.II vocabulary than the full library. We stop before:

- adaptive refinement algorithms and error estimators;
- solution transfer across changing meshes;
- hp finite elements;
- discontinuous Galerkin formulations beyond recognizing discontinuous control spaces;
- mixed finite elements and inf-sup implementation details;
- distributed meshes and MPI;
- PETSc/Trilinos distributed linear algebra in detail;
- matrix-free `MatrixFree`/`FEEvaluation` programming;
- multigrid implementation;
- manifold descriptions and high-order geometry beyond the mapping concept;
- automatic differentiation and symbolic differentiation facilities;
- parameter-handler frameworks;
- advanced postprocessing and data interpretation;
- mesh generators beyond what is needed to orient a simple example.

The deal.II tutorial is designed precisely so that these topics can be learned incrementally after the basic pipeline is understood.

## Where to go next

- **Default continuation into `nmopt`:** use the manual according to what you are trying to understand:
  - [manual overview · Numerical realization](../overview/numerical-realization.md), for how FE coordinates, operators, metrics, and solves fit into the project architecture;
  - [manual 11 · Compilation and lowering](../concepts/11-compilation-and-lowering.md), if `nmopt` is constructing the deal.II realization from a semantic problem;
  - [manual 12 · Authoring and using compiled problems](../concepts/12-authoring-and-using-compiled-problems.md), for the user workflow around compiled problems;
  - [manual 13 · Integrating an existing PDE application](../concepts/13-integrating-an-existing-pde-application.md), if an application already owns its mesh, assembly, solves, and output.
- **For more deal.II practice:** the natural upstream tutorial sequence is step-1 (mesh creation), step-2 (DoF distribution and sparsity), step-3 (a complete Laplace solve), step-4 (dimension-independent programming and nonconstant data), and step-6 when hanging-node constraints and adaptive refinement become relevant.
- **For an interactive route:** Wolfgang Bangerth's MATH 676 lectures cover the same early workflow and later discuss quadrature, solvers, and preconditioning. The recordings can use older APIs, so prefer the current deal.II tutorial whenever an interface differs.

## References and further reading

Publication details, stable links, and access notes are collected in the [Background reference catalogue](references.md).

- **deal.II current tutorial.**
  - step-1 for grids.
  - step-2 for degrees of freedom and sparsity.
  - step-3 for the complete Laplace assembly/solve/output pipeline.
  - step-4 for dimension-independent programming and nonconstant data.
  - step-6 for hanging-node constraints and adaptivity.
- **deal.II class documentation.**
  - Consult `DoFHandler`, `FiniteElement`/`FE_Q`, `FEValues`, `FEFaceValues`, `AffineConstraints`, `SolverControl`, solver classes, `Function`, and `DataOut` when exact signatures or lifetime guarantees matter.
- **Wolfgang Bangerth, MATH 676 video lectures.**
  - Lectures 4, 9, 10, 12, 13, and 16 are particularly close to this chapter's early progression.
  - Lectures 33.5, 34, 35, 37, and 38 provide later orientation on quadrature, solvers, and preconditioners.
  - The video page warns that older recordings can differ from current deal.II APIs; use the current tutorial as the present-day authority.
