# Background for the `nmopt` manual

This directory collects the mathematical, numerical, and software background that the `nmopt` manual assumes but does not teach in full. It is not a second manual and it is not intended to be read linearly from beginning to end. Use it as a prerequisite route: start where your background becomes uncertain, follow the default continuation for your goal, and return to the `nmopt` manual once the missing concepts are comfortable.

The numbered main chapters form the core prerequisite material. Chapter links keep the file number in their visible label (for example, `03 · Finite elements`) so that cross-references remain recognizable even when the manual and background sequences are open at the same time. Files with an `a` suffix are optional companions that deepen nearby theory without blocking the main route. [07b · Extensions of PDE-constrained optimization](07b-extensions-of-pde-constrained-optimization.md) has a different role: it surveys broader directions that are useful for orientation but are not prerequisites for the current project. [08a · Scientific software workflow: shell, compilation, CMake, and tests](08a-scientific-software-workflow.md) is a practical companion rather than a deeper C++ theory chapter.

## Choose a starting point

Use the first row that describes something you would not yet be comfortable explaining or using.

| If this is unfamiliar | Start here |
| --- | --- |
| $L^2$, $H^1$, weak derivatives, test functions, weak boundary conditions | [01 · Function spaces and weak PDEs](01-function-spaces-and-weak-pdes.md) |
| dual spaces, covectors, Riesz maps, Fréchet derivatives, JVPs/VJPs, adjoints | [02 · Duality, derivatives, and adjoints](02-duality-derivatives-and-adjoints.md) |
| meshes, finite elements, degrees of freedom, quadrature, assembly, mass/stiffness matrices | [03 · Finite elements](03-finite-elements.md) |
| conditioning, CG/MINRES/GMRES, preconditioning, Schur complements, saddle-point systems | [04 · Numerical linear algebra for PDE and optimization systems](04-numerical-linear-algebra.md) |
| line search, Armijo/Wolfe conditions, nonlinear CG, Newton, BFGS/L-BFGS, trust regions | [05 · Unconstrained numerical optimization](05-unconstrained-optimization.md) |
| feasible directions, projections, KKT conditions, complementarity, active sets, PDAS | [06 · Constrained optimization](06-constrained-optimization.md) |
| control-to-state maps, sensitivities, adjoint gradients, reduced versus all-at-once PDE optimization | [07 · PDE-constrained optimization](07-pde-constrained-optimization.md) |
| references, pointers, ownership, templates, virtual interfaces, lambdas, callbacks, `std::function` | [08 · Modern C++ for scientific code](08-modern-cpp-for-scientific-code.md) |
| shell commands, compilation/linking, Make/Ninja, CMake, build directories, tests | [08a · Scientific software workflow: shell, compilation, CMake, and tests](08a-scientific-software-workflow.md) |
| `Triangulation`, `DoFHandler`, `FEValues`, `AffineConstraints`, deal.II assembly and solvers | [09 · The deal.II finite-element workflow](09-dealii-finite-element-workflow.md) |

## Reading paths

There is no single required order. The most common routes are:

```text
PDE / finite-element route
    01  Function spaces and weak PDEs
        → 03  Finite elements
        → 04  Numerical linear algebra

Optimization route
    01  Function spaces and weak PDEs
        → 02  Duality, derivatives, and adjoints
        → 05  Unconstrained numerical optimization
        → 06  Constrained optimization
        → 07  PDE-constrained optimization

Software route
    08  Modern C++ for scientific code
        → 08a Scientific software workflow
        → 09  The deal.II finite-element workflow
```

The routes meet rather than replace one another. For example, the PDE-constrained optimization chapter assumes the derivative/adjoint language from the optimization route and uses the weak-PDE language from the PDE route. A reader interested in the discretized examples will also benefit from the finite-element and numerical-linear-algebra chapters.

## Main chapters and companions

| Core chapter | Optional continuation |
| --- | --- |
| [01 · Function spaces and weak PDEs](01-function-spaces-and-weak-pdes.md) | [01a · Further notes on function spaces and weak PDEs](01a-further-function-space-notes.md): higher Sobolev spaces, embeddings, weak convergence, distributions, regularity, inf-sup ideas |
| [02 · Duality, derivatives, and adjoints](02-duality-derivatives-and-adjoints.md) | no separate companion; later chapters reuse these concepts directly |
| [03 · Finite elements](03-finite-elements.md) | [03a · Further finite-element notes](03a-further-finite-element-notes.md): approximation estimates, quadrature consistency, adaptivity, DG/nonconforming and mixed methods |
| [04 · Numerical linear algebra for PDE and optimization systems](04-numerical-linear-algebra.md) | [04a · Further numerical linear algebra for PDE-constrained optimization](04a-further-linear-algebra-for-pde-optimization.md): mesh/parameter robustness, multigrid, domain decomposition, and block/KKT preconditioning |
| [05 · Unconstrained numerical optimization](05-unconstrained-optimization.md) | [05a · Further notes on unconstrained optimization](05a-further-unconstrained-optimization-notes.md): existence/convergence background, infinite-dimensional orientation, and additional large-scale methods |
| [06 · Constrained optimization](06-constrained-optimization.md) | [06a · Further notes on constrained optimization](06a-further-constrained-optimization-notes.md): constraint qualifications, sign conventions, second-order conditions, and nearby methods |
| [07 · PDE-constrained optimization](07-pde-constrained-optimization.md) | [07a · Further notes on PDE-constrained optimization](07a-further-pde-constrained-optimization-notes.md) for nearby generalizations; [07b · Extensions of PDE-constrained optimization](07b-extensions-of-pde-constrained-optimization.md) for broader directions beyond the current project |
| [08 · Modern C++ for scientific code](08-modern-cpp-for-scientific-code.md) | [08a · Scientific software workflow: shell, compilation, CMake, and tests](08a-scientific-software-workflow.md): practical build and command-line background |
| [09 · The deal.II finite-element workflow](09-dealii-finite-element-workflow.md) | continue through the official deal.II tutorials as needed |

The shared [background reference catalogue](references.md) records the main books, course material, official documentation, and further-reading sources used across the chapters.

## From the `nmopt` manual to the right background

The `nmopt` manual is organized around the current project, so its numbering is independent of the background sequence. Cross-references therefore use document titles rather than bare chapter numbers.

| If you are reading the `nmopt` manual about... | Background that is most useful |
| --- | --- |
| [manual 01 · Anatomy of a discrete PDE-constrained problem](../concepts/01-discrete-problem-anatomy.md) | [01 · Function spaces and weak PDEs](01-function-spaces-and-weak-pdes.md); [03 · Finite elements](03-finite-elements.md) |
| [manual 02 · Spaces, coordinates, and duality](../concepts/02-spaces-coordinates-and-duality.md) | [02 · Duality, derivatives, and adjoints](02-duality-derivatives-and-adjoints.md); [03 · Finite elements](03-finite-elements.md) |
| [manual 03 · Operators, derivatives, and adjoints](../concepts/03-operators-derivatives-and-adjoints.md) | [02 · Duality, derivatives, and adjoints](02-duality-derivatives-and-adjoints.md) |
| [manual 04 · Metrics, gradients, and constraints](../concepts/04-metrics-gradients-and-constraints.md) | [02 · Duality, derivatives, and adjoints](02-duality-derivatives-and-adjoints.md); [06 · Constrained optimization](06-constrained-optimization.md) as needed |
| [manual 05 · Reduced state–adjoint formulation](../concepts/05-reduced-state-adjoint-formulation.md) | [07 · PDE-constrained optimization](07-pde-constrained-optimization.md) |
| [manual 06 · Reduced optimization methods](../concepts/06-reduced-optimization-methods.md) | [05 · Unconstrained numerical optimization](05-unconstrained-optimization.md); [04 · Numerical linear algebra for PDE and optimization systems](04-numerical-linear-algebra.md) |
| [manual 07 · Optimality systems and KKT](../concepts/07-optimality-systems-and-kkt.md) | [06 · Constrained optimization](06-constrained-optimization.md); [04 · Numerical linear algebra for PDE and optimization systems](04-numerical-linear-algebra.md); [07 · PDE-constrained optimization](07-pde-constrained-optimization.md) |
| [manual 08 · Complementarity and PDAS](../concepts/08-complementarity-and-pdas.md) | [06 · Constrained optimization](06-constrained-optimization.md) and its [06a · Further notes on constrained optimization](06a-further-constrained-optimization-notes.md) |
| [manual 11 · Compilation and lowering](../concepts/11-compilation-and-lowering.md) | [03 · Finite elements](03-finite-elements.md); [08 · Modern C++ for scientific code](08-modern-cpp-for-scientific-code.md); [09 · The deal.II finite-element workflow](09-dealii-finite-element-workflow.md) |
| [manual 12 · Authoring and using compiled problems](../concepts/12-authoring-and-using-compiled-problems.md) | [08 · Modern C++ for scientific code](08-modern-cpp-for-scientific-code.md); [08a · Scientific software workflow: shell, compilation, CMake, and tests](08a-scientific-software-workflow.md) |
| [manual 13 · Integrating an existing PDE application](../concepts/13-integrating-an-existing-pde-application.md) | [07 · PDE-constrained optimization](07-pde-constrained-optimization.md); [08 · Modern C++ for scientific code](08-modern-cpp-for-scientific-code.md); [09 · The deal.II finite-element workflow](09-dealii-finite-element-workflow.md) |

## When to return to the manual

You do not need to finish every companion before continuing. The main chapters are written to provide competence endpoints; the companions are there when a later topic exposes a gap or when you want more theory.

For the mathematical route, [07 · PDE-constrained optimization](07-pde-constrained-optimization.md) is the final core background chapter. Once its state–sensitivity–adjoint and reduced/all-at-once viewpoints are comfortable, continue with the `nmopt` manual's [manual 05 · Reduced state–adjoint formulation](../concepts/05-reduced-state-adjoint-formulation.md), [manual 06 · Reduced optimization methods](../concepts/06-reduced-optimization-methods.md), and [manual 07 · Optimality systems and KKT](../concepts/07-optimality-systems-and-kkt.md).

For implementation, [09 · The deal.II finite-element workflow](09-dealii-finite-element-workflow.md) is the bridge back to the manual's [manual overview · Numerical realization](../overview/numerical-realization.md), [manual 11 · Compilation and lowering](../concepts/11-compilation-and-lowering.md), and [manual 13 · Integrating an existing PDE application](../concepts/13-integrating-an-existing-pde-application.md).
