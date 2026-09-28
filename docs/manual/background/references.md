# Background reference catalogue

**Background navigation:** [Index](README.md)

## Purpose

This file collects recurring references used across the background chapters. Chapter-local **References and further reading** sections use one primary list item per source and nested items for the chapters, sections, lectures, or documentation pages most relevant to that chapter. Publication data, access links, and source-family context live here so that the chapter bibliographies can stay compact and consistent.

The catalogue keeps general books and papers separate from the NMPDE and NMOPT course ecosystems, University of Pisa course-local notes, and software/tool documentation.

## Note on provenance

These notes were drafted and revised with LLM assistance. The bibliography lists sources explicitly consulted or supplied for this documentation; it is not a complete provenance record for information present in the model's training data.

## Public finite-element and PDE references

### `ciarlet-1978`

Philippe G. Ciarlet, *The Finite Element Method for Elliptic Problems*, North-Holland, 1978.

- Classical reference for conforming finite elements for elliptic problems.
- Most relevant to the theory behind [03 · Finite elements](03-finite-elements.md) and its companion: finite-element construction, interpolation, approximation, and elliptic error analysis.
- Listed among the principal references of the NMPDE course.

### `ern-guermond-2004`

Alexandre Ern and Jean-Luc Guermond, *Theory and Practice of Finite Elements*, Applied Mathematical Sciences 159, Springer, 2004.

- Springer: <https://link.springer.com/book/10.1007/978-1-4757-4355-5>
- Particularly useful for this background track:
  - **Finite Element Interpolation**;
  - **Approximation in Banach Spaces by Galerkin Methods**;
  - **Coercive Problems** and **Mixed Problems**;
  - **Quadratures, Assembling, and Storage**;
  - **Linear Algebra**.
- Combines finite-element theory with implementation-oriented material and is a useful continuation of [03 · Finite elements](03-finite-elements.md) and [04 · Numerical linear algebra for PDE and optimization systems](04-numerical-linear-algebra.md).

### `brenner-scott-2008`

Susanne C. Brenner and L. Ridgway Scott, *The Mathematical Theory of Finite Element Methods*, third edition, Texts in Applied Mathematics 15, Springer, 2008.

- Springer: <https://link.springer.com/book/10.1007/978-0-387-75934-0>
- Particularly relevant:
  - **The Construction of a Finite Element Space**;
  - **Polynomial Approximation Theory in Sobolev Spaces**;
  - **Finite Element Multigrid Methods**;
  - **Additive Schwarz Preconditioners**;
  - **Adaptive Meshes**;
  - **Variational Crimes**;
  - **Mixed Methods**.
- The current NMPDE bibliography identifies the same third edition through its later softcover publication metadata; the Springer page gives copyright 2008.

### `boffi-brezzi-fortin-2013`

Daniele Boffi, Franco Brezzi, and Michel Fortin, *Mixed Finite Element Methods and Applications*, Springer Series in Computational Mathematics 44, Springer, 2013.

- Springer: <https://link.springer.com/book/10.1007/978-3-642-36519-5>
- Particularly relevant after the main finite-element prerequisite:
  - **Variational Formulations and Finite Element Methods**;
  - **Function Spaces and Finite Element Approximations**;
  - **Algebraic Aspects of Saddle Point Problems**;
  - **Saddle Point Problems in Hilbert Spaces**;
  - **Approximation of Saddle Point Problems**.
- Deeper reading for mixed methods, inf-sup theory, and the connection between mixed FEM and block saddle-point algebra.

### `quarteroni-valli-1994`

Alfio Quarteroni and Alberto Valli, *Numerical Approximation of Partial Differential Equations*, 1994.

- A broad numerical-PDE reference that appears in both the NMPDE and University of Pisa numerical-analysis bibliographies.
- Useful as a continuation when approximation theory and PDE discretization need more depth than the prerequisite chapters provide.

## Public numerical linear algebra and scientific-computing references

### `demmel-1997`

James W. Demmel, *Applied Numerical Linear Algebra*, SIAM, 1997.

- SIAM: <https://epubs.siam.org/doi/book/10.1137/1.9781611971446>
- Particularly relevant: Chapter 2, **Linear Equation Solving**, and Chapter 6, **Iterative Methods for Linear Systems**.
- Good general reference for conditioning, direct solvers, iterative methods, sparse problems, and the relation between algorithms and computer architecture.

### `trefethen-bau-1997`

Lloyd N. Trefethen and David Bau III, *Numerical Linear Algebra*, SIAM, 1997.

- Author page: <https://people.maths.ox.ac.uk/trefethen/text.html>
- DOI route: <https://doi.org/10.1137/1.9780898719574>
- Particularly relevant: Part III, **Conditioning and Stability**; Part IV, **Systems of Equations**; and Part VI, **Iterative Methods**, including Arnoldi, GMRES, Lanczos, conjugate gradients, and preconditioning.

### `saad-2003`

Yousef Saad, *Iterative Methods for Sparse Linear Systems*, second edition, SIAM, 2003.

- Author-hosted edition: <https://www-users.cse.umn.edu/~saad/IterMethBook_2ndEd.pdf>
- SIAM: <https://epubs.siam.org/doi/book/10.1137/1.9780898718003>
- Particularly relevant: Chapter 3 on sparse matrices, Chapter 5 on projection methods, Chapters 6–7 on Krylov methods, Chapters 9–10 on preconditioning, Chapter 13 on multigrid, and Chapter 14 on domain decomposition.
- Best continuation when [04 · Numerical linear algebra for PDE and optimization systems](04-numerical-linear-algebra.md) becomes a primary topic.

### `bini-capovani-menchi-1988`

Dario Bini, Milvio Capovani, and Ornella Menchi, *Metodi Numerici per l'Algebra Lineare*, Zanichelli, 1988.

- Public author-hosted PDF: <https://people.dm.unipi.it/bini/Metodi_Numerici_per_l%27Algebra_Lineare.pdf>
- Particularly relevant: vector and matrix norms, direct linear solvers, stationary iterations, relaxation methods, and conjugate gradients.
- Useful Italian-language reference and a recurring source in University of Pisa numerical-analysis courses.

### `benzi-golub-liesen-2005`

Michele Benzi, Gene H. Golub, and Jörg Liesen, “Numerical solution of saddle point problems,” *Acta Numerica* 14 (2005), 1–137.

- DOI: <https://doi.org/10.1017/S0962492904000212>
- Advanced reference for saddle-point systems, Schur complements, iterative solvers, and block preconditioning.
- Especially relevant to [04a · Further numerical linear algebra for PDE-constrained optimization](04a-further-linear-algebra-for-pde-optimization.md) and later KKT material.

### `barrett-et-al-1994`

Richard Barrett, Michael Berry, Tony F. Chan, James Demmel, June Donato, Jack Dongarra, Victor Eijkhout, Roldan Pozo, Charles Romine, and Henk van der Vorst, *Templates for the Solution of Linear Systems: Building Blocks for Iterative Methods*, SIAM, 1994.

- Public Netlib edition: <https://www.netlib.org/templates/Templates.html>
- DOI: <https://doi.org/10.1137/1.9781611971538>
- Compact algorithm-oriented reference for classical iterative solver families and preconditioners.

### `leveque-2007`

Randall J. LeVeque, *Finite Difference Methods for Ordinary and Partial Differential Equations: Steady-State and Time-Dependent Problems*, SIAM, 2007.

- DOI: <https://doi.org/10.1137/1.9780898717839>
- Chapters 3–4 connect elliptic PDE discretization to sparse linear systems and iterative solution.
- More PDE-oriented than the general linear-algebra books above.

## Public numerical-optimization references

### `nocedal-wright-1999`

Jorge Nocedal and Stephen J. Wright, *Numerical Optimization*, first edition, Springer, 1999.

- Springer first-edition record: <https://link.springer.com/book/10.1007/b98874>
- The University of Pisa optimization course lists this first edition; chapter and section locators in these notes therefore refer to the first-edition structure.
- Verified first-edition topics especially relevant to the next background chapters include **Fundamentals of Unconstrained Optimization**, **Line Search Methods**, **Trust-Region Methods**, **Conjugate Gradient Methods**, **Practical Newton Methods**, **Quasi-Newton Methods**, **Theory of Constrained Optimization**, **Fundamentals of Algorithms for Nonlinear Constrained Optimization**, and **Quadratic Programming**.
- This is the main general-purpose continuation for [05 · Unconstrained numerical optimization](05-unconstrained-optimization.md) and [06 · Constrained optimization](06-constrained-optimization.md).

### `bazaraa-sherali-shetty-1993`

Mokhtar S. Bazaraa, Hanif D. Sherali, and C. M. Shetty, *Nonlinear Programming: Theory and Algorithms*, second edition, Wiley, 1993.

- Listed among the principal references of the University of Pisa optimization course.
- Particularly relevant at topic level for convexity, optimality conditions, constraint qualifications, Lagrangian duality, unconstrained methods, and feasible-direction methods.
- References to this book are kept at topic level here rather than assigning unverified chapter or section numbers.

### `bertsekas-1999`

Dimitri P. Bertsekas, *Nonlinear Programming*, second edition, Athena Scientific, 1999.

- Listed among the principal references of the University of Pisa optimization course; some course pages cite later reprints/printings, so the edition should be stated whenever an exact locator is used.
- Particularly relevant at topic level for unconstrained optimization, gradient and Newton methods, conjugate directions, multiplier theory, constrained optimization, and duality.
- References to this book are kept at topic level here because chapter numbering varies across editions and printings.

### `beck-2017`

Amir Beck, *First-Order Methods in Optimization*, MOS-SIAM Series on Optimization 25, SIAM, 2017.

- SIAM: <https://epubs.siam.org/doi/book/10.1137/1.9781611974997>
- Especially useful when projection, nonsmooth convex models, or proximal viewpoints become relevant.
- The publisher table of contents verifies Chapter 8, **Primal and Dual Projected Subgradient Methods**, and Chapter 10, **The Proximal Gradient Method**. These are deeper continuations rather than prerequisites for the main unconstrained chapter.

## Public PDE-constrained-optimization references

### `manzoni-quarteroni-salsa-2021`

Andrea Manzoni, Alfio Quarteroni, and Sandro Salsa, *Optimal Control of Partial Differential Equations: Analysis, Approximation, and Applications*, Applied Mathematical Sciences 207, Springer, 2021.

- DOI: <https://doi.org/10.1007/978-3-030-77226-0>
- Appendix A is especially useful for functional analysis and derivatives.
- Appendix B provides matrix properties, saddle-point systems, Galerkin approximation, and finite elements.
- Chapters 3–6 cover optimization, optimality systems, elliptic optimal control, discretization, and numerical algorithms.

### `troeltzsch-2010`

Fredi Tröltzsch, *Optimal Control of Partial Differential Equations: Theory, Methods and Applications*, AMS, 2010.

- One of the primary references listed by the public NMOPT course.
- Useful for continuous optimal-control theory, elliptic/parabolic control, optimality conditions, and control constraints.

### `de-los-reyes-2015`

Juan Carlos De los Reyes, *Numerical PDE-Constrained Optimization*, Springer, 2015.

- One of the primary references listed by the public NMOPT course.
- Particularly relevant when numerical algorithms, semismooth/active-set methods, and PDE-constrained optimization are treated together.

### `hinze-et-al-2009`

Michael Hinze, René Pinnau, Michael Ulbrich, and Stefan Ulbrich, *Optimization with PDE Constraints*, Springer, 2009.

- Appears in later NMOPT lecture references.
- Useful continuation for reduced methods, discretization, control constraints, and PDE-optimization algorithms.

## NMPDE course ecosystem

These entries record the public NMPDE course context from which several prerequisites originate.

### `nmpde-course-page`

Luca Heltai, **Numerical Methods for Partial Differential Equations**, University of Pisa.

- Teaching page: <https://luca-heltai.github.io/numerical-methods-for-pdes/>
- Public course book: <https://luca-heltai.github.io/nmpde/>
- GitHub repository: <https://github.com/luca-heltai/nmpde>

The public course description emphasizes finite-element-space construction, polynomial approximation in Sobolev spaces, convergence for elliptic problems, variational crimes, mixed methods, adaptivity, and deal.II practice.

### `mancini-benvenuti-heltai-2025`

Stefano Mancini, Andrea Benvenuti, and Luca Heltai, *Numerical Methods for Partial Differential Equations*, University of Pisa lecture notes, 2025.

- Public source repository: <https://github.com/luca-heltai/nmpde-notes>
- Closely aligned with [01 · Function spaces and weak PDEs](01-function-spaces-and-weak-pdes.md), [02 · Duality, derivatives, and adjoints](02-duality-derivatives-and-adjoints.md), and [03 · Finite elements](03-finite-elements.md): Sobolev spaces, weak formulations, Galerkin approximation, finite elements, interpolation/error analysis, and mixed formulations.

### `nmpde-reference-set`

The current public NMPDE bibliography collects, among others:

- Ciarlet, *The Finite Element Method for Elliptic Problems*;
- Ern and Guermond, *Theory and Practice of Finite Elements*;
- Brenner and Scott, *The Mathematical Theory of Finite Element Methods*;
- Boffi, Brezzi, and Fortin, *Mixed Finite Element Methods and Applications*;
- Quarteroni and Valli, *Numerical Approximation of Partial Differential Equations*.

The course also points to Wolfgang Bangerth's video lectures for **MATH 676: Finite element methods in scientific computing**:

- <https://www.math.colostate.edu/~bangerth/videos.html>

These are course-recommended continuations, not additional prerequisites for the `nmopt` background track.

## NMOPT course ecosystem

These entries record the NMOPT optimal-control course separately from the NMPDE course, even where the two share prerequisites.

### `nmopt-course-page`

Luca Heltai, **Numerical Methods for Optimal Control**, University of Pisa.

- Teaching page: <https://luca-heltai.github.io/numerical-methods-for-optimal-control/>
- Public course material: <https://luca-heltai.github.io/nmopt/>
- GitHub repository: <https://github.com/luca-heltai/nmopt>

The public course description covers optimal-control formulations, gradient methods, direct and indirect approaches for ODE/PDE-constrained control, numerical techniques, implementations, and case studies.

### `nmopt-primary-reference-set`

The public NMOPT course repository lists three primary references:

- Fredi Tröltzsch, *Optimal Control of Partial Differential Equations*, AMS, 2010;
- Andrea Manzoni, Alfio Quarteroni, and Sandro Salsa, *Optimal Control of Partial Differential Equations*, Springer, 2021;
- Juan Carlos De los Reyes, *Numerical PDE-Constrained Optimization*, Springer, 2015.

Individual lectures add more specialized references, including Hinze–Pinnau–Ulbrich–Ulbrich for optimization with PDE constraints and classical sources for time-dependent and boundary-control topics.

This section records the course context and a public continuation route; the chapter-local bibliographies point to the underlying books or papers when a specific source is more useful.

## General numerical-optimization course ecosystem

These materials provide general finite-dimensional optimization background for [05 · Unconstrained numerical optimization](05-unconstrained-optimization.md) and [06 · Constrained optimization](06-constrained-optimization.md).

### `unipi-ottimizzazione-non-lineare-2026-27`

University of Pisa, **Ottimizzazione Non Lineare**, Massimo Pappalardo, academic year 2026/27.

The course covers classification and equivalent transformations of optimization problems, convex functions and sets, tangent/linearized/feasible/critical cones, constraint regularity, first- and second-order optimality conditions, duality, and algorithms for unconstrained and constrained nonlinear optimization. The listed reference text is:

- A. Carpignani and M. Pappalardo, *Theory and Methods of Optimization*, Springer, 2025.

### `unipi-tmo-reference-set`

An earlier University of Pisa optimization-course bibliography remains useful as supplementary reading:

- <https://pages.di.unipi.it/bigi/dida/tmo/2021/bibliografia.html>
- J. Nocedal and S. J. Wright, *Numerical Optimization*, Springer, 1999;
- M. S. Bazaraa, H. D. Sherali, and C. M. Shetty, *Nonlinear Programming: Theory and Algorithms*, Wiley, 1993;
- D. P. Bertsekas, *Nonlinear Programming*, Athena Scientific;
- A. Beck, *First-Order Methods in Optimization*, SIAM, 2017.

## University of Pisa course-local notes

These notes are useful to fellow University of Pisa students because their level and notation are close to the courses from which this project grew. They should not be treated as public prerequisites: availability may depend on the University e-learning platform or institutional access.

### `unipi-analisi-numerica-2025`

Dario A. Bini, Beatrice Meini, and Leonardo Robol, *Dispense del corso Analisi Numerica con Laboratorio*, University of Pisa, updated 21 November 2025.

- Relevant material: its Chapter 3 on numerical linear algebra and conditioning, Chapter 4 on direct methods for linear systems, and Chapter 5 on iterative methods.
- Its bibliography points, among other sources, to Bini–Capovani–Menchi, Golub–Van Loan, and Higham.

### `unipi-calcolo-scientifico-2025-26`

*Scientific Computing / Calcolo Scientifico*, University of Pisa lecture notes, academic year 2025/26.

- The notes identify Demmel's *Applied Numerical Linear Algebra* as the reference book.
- Especially relevant to [04 · Numerical linear algebra for PDE and optimization systems](04-numerical-linear-algebra.md): Chapter 6 on Krylov methods for linear systems, §6.8 on GMRES preconditioning, and §6.10 on Lanczos, MINRES, and conjugate gradients.
- The bibliography also points to Benzi–Golub–Liesen for saddle-point systems and Saad for iterative methods.

### `unipi-ian-2026`

Dario A. Bini, Paola Boito, and Beatrice Meini, *Appunti di Istituzioni di Analisi Numerica*, University of Pisa, updated 13 May 2026.

- Particularly relevant to [03 · Finite elements](03-finite-elements.md): Chapter 2 develops interpolatory quadrature, convergence, Newton–Cotes, Fejér/Clenshaw–Curtis, and Gaussian formulas; §7.8 introduces variational methods.
- Its bibliography includes Barrett et al., Bini–Capovani–Menchi, and Quarteroni–Valli among the sources useful for numerical PDE work.

## Suggested C++ references — not vetted for this background sequence

The references in this section are **optional suggestions and have not been vetted as part of this background sequence**. They are offered as possible continuations, not as prerequisites.

### `stroustrup-tour-2022`

Bjarne Stroustrup, *A Tour of C++*, third edition, Addison-Wesley, 2022.

- Author page: <https://www.stroustrup.com/tour3.html>
- A compact overview of modern C++ for readers who already have programming experience.
- Useful continuation for [08 · Modern C++ for scientific code](08-modern-cpp-for-scientific-code.md) when a broader view of the language and standard library is desired.
- The third edition discusses C++20, while `nmopt` currently requires only the C++17 subset used in its interfaces.

### `stroustrup-ppp-2024`

Bjarne Stroustrup, *Programming: Principles and Practice Using C++*, third edition, Addison-Wesley, 2024.

- Author page: <https://www.stroustrup.com/programming.html>
- More introductory and pedagogical than *A Tour of C++*; potentially useful when the reader's general programming background is weaker.
- It covers substantially more programming material than is needed to read `nmopt`.

### `cpp-core-guidelines`

Bjarne Stroustrup and Herb Sutter, eds., *C++ Core Guidelines*.

- <https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines>
- A living design-guideline document rather than a sequential language tutorial.
- Particularly relevant to [08 · Modern C++ for scientific code](08-modern-cpp-for-scientific-code.md) for interfaces, parameter passing, ownership, RAII, class hierarchies, const-correctness, lambdas, and resource management.
- The guidelines explicitly cover modern C++ including advice applicable to C++17.

### `cppreference`

cppreference, C++ language and standard-library reference.

- <https://en.cppreference.com/cpp>
- Useful as a precise lookup reference for language rules and standard-library facilities such as smart pointers, `std::move`, lambdas, `std::function`, `std::optional`, and `std::variant`.
- Better used after a concept has been introduced than as a first tutorial.

## Scientific-software and build-tool references

These references support [08a · Scientific software workflow: shell, compilation, CMake, and tests](08a-scientific-software-workflow.md). Unlike the suggested C++ books above, most are official tool manuals; they are authoritative for tool behavior but are not intended to be read cover to cover as a course.

### `cmake-current`

CMake, current official documentation.

- Getting started: <https://cmake.org/getting-started/>
- Documentation index: <https://cmake.org/documentation/>
- The step-by-step CMake Tutorial is the main introductory route.
- The `cmake-buildsystem(7)`, `cmake-presets(7)`, and CTest manuals are useful once targets, presets, and registered tests become relevant.

### `bash-manual`

GNU Bash Reference Manual.

- <https://www.gnu.org/software/bash/manual/>
- Authoritative reference for shell command execution, quoting, expansions, variables, redirection, pipes, and exit status.
- Used as a reference for [08a · Scientific software workflow: shell, compilation, CMake, and tests](08a-scientific-software-workflow.md) rather than as a required shell course.

### `gnu-make-manual`

GNU Make manual.

- <https://www.gnu.org/software/make/manual/>
- Useful for understanding targets, prerequisites, recipes, dependency updating, and variables when Make appears as the build executor beneath a generated build system.

### `ninja-manual`

Ninja build system manual.

- <https://ninja-build.org/manual.html>
- Useful for understanding Ninja as a low-level build executor designed to run generated dependency graphs efficiently.
- In this background track, CMake remains the project-description layer; generated `build.ninja` files are not intended for manual editing.

### `bangerth-math676-videos`

Wolfgang Bangerth, **MATH 676: Finite element methods in scientific computing**, video lectures.

- <https://www.math.colostate.edu/~bangerth/videos.html>
- Especially relevant software/workflow lectures include Lecture 2.9 on the command line, Lecture 2.91 on compiling programs, Lecture 8.01 on direct CMake project setup, Lecture 12 on C++ templates, Lecture 18 on Debug versus optimized mode, and Lectures 42–43 on scientific-computing workflows and large-software development.
- For deal.II itself, Lectures 4, 9, 10, 13, and 16 accompany the early tutorial progression; later lectures discuss quadrature, solvers, and preconditioners.
- The video page warns that deal.II evolves and older recordings can show APIs or tutorial code that have since changed. The current deal.II documentation should resolve such discrepancies.

## Software and library documentation

### `dealii-current`

The deal.II finite-element library, current documentation.

- Main documentation: <https://dealii.org/current/doxygen/deal.II/index.html>
- Tutorial: <https://dealii.org/current/doxygen/deal.II/Tutorial.html>
- For [09 · The deal.II finite-element workflow](09-dealii-finite-element-workflow.md), the most useful early tutorial path is step-1 (grids), step-2 (DoFs and sparsity), step-3 (Laplace assembly, solve, and output), step-4 (dimension-independent programming and nonconstant data), and step-6 when hanging-node constraints/adaptivity become relevant.
- Class documentation for `DoFHandler`, `FE_Q`, `FEValues`, `FEFaceValues`, `AffineConstraints`, `SolverControl`, solver classes, `Function`, and `DataOut` should be consulted when exact APIs or lifetime/ownership details matter.
- The deal.II documentation is used as an implementation reference rather than as a substitute for [03 · Finite elements](03-finite-elements.md) and [03a · Further finite-element notes](03a-further-finite-element-notes.md).

For linear solvers and scalable PDE preconditioning, useful tutorial continuations include:

- step-5: preconditioned conjugate gradients;
- step-16: geometric multigrid preconditioning for Laplace;
- step-31: Schur-complement block preconditioning with AMG building blocks;
- step-40: parallel elliptic solves with algebraic multigrid;
- step-50: comparison of matrix-based geometric, matrix-free geometric, and algebraic multigrid;
- step-56: geometric multigrid inside a block Stokes preconditioner;
- step-63: Schwarz-type block smoothers for a nonsymmetric multigrid problem.

These examples are implementation continuations for [04 · Numerical linear algebra for PDE and optimization systems](04-numerical-linear-algebra.md) and [09 · The deal.II finite-element workflow](09-dealii-finite-element-workflow.md). They are not additional prerequisites for the main background sequence.

## Additional numerical-analysis references

These broader references are not required by the current background chapters, but they are useful for deeper study of matrix computations, numerical stability, and general numerical analysis.

- Gene H. Golub and Charles F. Van Loan, *Matrix Computations*, fourth edition, Johns Hopkins University Press, 2013 – broad reference for matrix factorizations and numerical linear algebra.
- Nicholas J. Higham, *Accuracy and Stability of Numerical Algorithms*, second edition, SIAM, 2002 – deeper treatment of floating-point error, backward stability, and conditioning. DOI: <https://doi.org/10.1137/1.9780898718027>.
- Alfio Quarteroni, Riccardo Sacco, Fausto Saleri, and Paola Gervasio, *Matematica Numerica*, Springer, 2014 – broad Italian-language numerical-analysis reference.
