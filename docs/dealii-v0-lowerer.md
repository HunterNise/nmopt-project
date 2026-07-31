# deal.II v0 scalar diffusion-reaction lowerer

## Implemented slice

The header `include/nmopt/dealii/scalar_diffusion_reaction.hpp` is the first
deal.II lowerer for the executable contract. It compiles the concrete model

$$
  -\nabla\cdot(k\nabla y)+c y=f+u
  \quad\text{in }\Omega,\qquad
  y=0\quad\text{on }\Gamma_{D},
$$

with the objective

$$
  J(y,u)=\frac{1}{2}\lVert y-y_{d}\rVert_{L^{2}(\Omega)}^{2}+
         \frac{\alpha}{2}\lVert u\rVert_{L^{2}(\Omega)}^{2}.
$$

The emitted DTO objects are

$$
  r_{h}(y_{h},u_{h})=A_{h}y_{h}-f_{h}-B_{h}u_{h},
$$

$$
  J_{h}(y_{h},u_{h})=\frac{1}{2}y_{h}^{\mathsf T}M_{y}y_{h}-q_{y}^{\mathsf T}y_{h}+
               \frac{1}{2}\lVert y_{d}\rVert^{2}_{L^{2}(\Omega)}+
               \frac{\alpha}{2}u_{h}^{\mathsf T}M_{u}u_{h}.
$$

The source, target, and sign convention are thus exactly those required by
the backend-parametric `ExecutableModelT` contract.

## Selected discrete policies

| Concern | V0 deal.II policy |
|---|---|
| Geometry and execution | One static serial `Triangulation`; assembled `SparseMatrix` and `Vector` operations. |
| State and test | Scalar `FE_Q` state/test space, degree at least one. |
| Volume control | `FE_DGQ(0)` control space on the same active cells. |
| PDE coefficients | Constant scalar diffusion $k>0$ and reaction $c\geq0$. |
| Data rule | Forcing and target are deal.II `Function` objects evaluated at declared cell quadrature points. |
| Observation/loss | Distributed state tracking assembled as state mass matrix, target load, and target norm. |
| Essential boundary | Selected homogeneous Dirichlet boundary ids, defaulting to id zero. |
| State and adjoint solves | Serial CG with identity preconditioner; valid only because this v0 operator is symmetric positive definite. |
| Provenance | DTO. |

The control mass matrix is the $L^{2}$ metric and regularisation realization.
Because `FE_DGQ(0)` is cellwise constant, a future deal.II cellwise box
constraint may use coefficient clipping only with that declared $L^{2}$
metric.

## Affine constraints and coordinate policy

The lowerer creates and retains `AffineConstraints` for the homogeneous
Dirichlet conditions. Its current executable coordinate policy represents
the constrained state space with full serial vectors, zero constrained
entries, zero residual/objective coupling on constrained rows and columns,
and identity residual rows for essential values. This makes the same
coordinate convention visible to:

1. residual evaluation;
2. JVP and VJP;
3. state solve;
4. adjoint solve; and
5. objective derivative.

This is a valid fixed homogeneous-Dirichlet realization of the declared
restriction. It deliberately rejects hanging-node, periodic, and any other
non-Dirichlet affine constraint. Supporting them requires the more general
reconstruction/pullback policy described in the implementation-readiness
review; no attempt is made to silently treat them as fixed zero DoFs.

## Verified properties

The deal.II `CTest` case assembles a refined two-dimensional mesh and verifies:

1. the state solve residual;
2. residual JVP/VJP adjoint consistency;
3. residual finite-difference JVP consistency; and
4. the reduced DTO derivative by a state-recomputed Taylor difference.

The test runs through the real deal.II `DoFHandler`, `FEValues`, `FEFace`-independent
volume assembly, `AffineConstraints`, sparse matrices, and CG solves.

## Explicit exclusions

The lowerer returns an unsupported capability rather than attempting any of
the following:

- inhomogeneous, controlled, periodic, or hanging essential conditions;
- Neumann/Robin terms, boundary controls, and boundary observations;
- vector, mixed, DG state, nonmatching meshes, or distributed MPI vectors;
- variable/nonlinear coefficients, nonlinear state solves, and second
  derivatives;
- $H^{1}$, $H^{-1}$, fractional, or non-diagonal metric realizations;
- box constraints in a deal.II vector backend;
- time dependence, matrix-free execution, OTD, or KKT Newton.

These are extensions of the contract, not branches to insert into this
lowerer or its reduced DTO solver.
