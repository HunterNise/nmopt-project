# Describing and compiling a problem

The semantic/compiler path is the part of `nmopt` that turns a structured
description of a PDE optimal-control problem into a concrete numerical product.

Its purpose is not to replace finite-element mathematics with configuration. It is
to separate **what problem is being requested** from **how that request is realized
with deal.II objects**.

This is useful when many related problems share the same implementation machinery.
Instead of hand-writing a new application class for every combination of diffusion,
control type, observation, metric, or formulation, the project records those choices
in a semantic graph and lets a bounded compiler assemble a supported realization.

## 1. What the semantic description contains

A semantic problem description records concepts such as:

- named state, control, test, and observation spaces;
- variables and their roles;
- PDE residual terms;
- boundary and region selections;
- observations and losses;
- the optimization metric;
- constraints;
- the requested formulation product;
- discretization or realization commitments that cannot be inferred later.

The description is backend-neutral in the sense that it does not own
`Triangulation`, `DoFHandler`, sparse matrices, or linear solvers.

It is not completely abstract mathematics, however. Some numerical choices are part
of the problem meaning and therefore need to be stated explicitly. A control may be
cellwise constant rather than continuous, an observation may act only on a material
region, or a boundary control may use one specific normal/conormal convention.
Those are semantic commitments even though they affect the eventual discretization.

## 2. A concrete example: Chapter 6 B1

Here **Chapter 6** refers to the numerical-methods chapter of *Optimal Control of
Partial Differential Equations* by Manzoni, Quarteroni, and Salsa. **B1** is the
repository's stable scenario ID for its distributed Laplace-control example.

The B1 application is a distributed scalar-control problem. At a conceptual level,
the scenario requests something like

```math
\begin{aligned}
-\nabla\cdot(\kappa\nabla y) + c y &= f + u
&& \text{in } \Omega, \\
y &= 0
&& \text{on } \partial\Omega,
\end{aligned}
```

Here $y$ is the state, $u$ the distributed control, $\kappa$ the diffusion
coefficient, $c$ a reaction coefficient, and $f$ the forcing. The desired state
$`y_{\mathrm{d}}`$ and regularization weight $\beta$ enter an objective of the form

```math
J(y,u)
=
\frac{1}{2}\lVert y-y_{\mathrm{d}}\rVert_{L^{2}(\Omega)}^{2}
+
\frac{\beta}{2}\lVert u\rVert_{L^{2}(\Omega)}^{2}.
```

The semantic request needs to say more than those equations alone. It also chooses,
for example, the control representation and its metric, the observation region, the
reduced formulation, and the relevant realization policies.

The Chapter 6 application code then follows this shape:

```text
typed scenario
     │
     ▼
problem recipe
     │
     ▼
backend-neutral problem description
     │
     ├─ validation
     ▼
resolved request
     │
     ├─ forcing / target / coefficients
     ├─ mesh and FE policy
     └─ compiler policy
     ▼
owned deal.II compilation session
     │
     ▼
compiled reduced product
```

The typed scenario belongs to the project application layer. The semantic graph
itself is the reusable problem description.

## 3. Why validation is split into two stages

A useful distinction is the difference between **semantic validity** and
**compiler support**.

Suppose a problem says:

- there is one scalar state and one scalar control;
- the PDE has a diffusion term and a distributed control term;
- the objective observes the state over a declared region;
- the control uses a positive metric;
- the requested product is reduced optimization.

The semantic validator can check whether those declarations are internally
consistent. It can catch things such as a missing variable, a region reference that
does not exist, or an incompatible role connection.

But semantic consistency does not prove that the current deal.II compiler knows how
to construct the requested numerical realization.

The compiler may still reject the request because, for example:

- a particular control discretization is not registered for that PDE family;
- an observation realization is unavailable;
- a requested boundary form is unsupported;
- the chosen formulation product has no lowering strategy for that combination.

This split prevents the semantic layer from becoming a disguised list of current
deal.II implementation cases.

## 4. Resolution closes the request before numerical construction

Validation answers whether the graph makes sense. Resolution then makes choices that
must be fixed before lowering.

The important property is that the compiler should receive a **closed request**,
rather than repeatedly guessing defaults while constructing numerical objects.

That makes the compilation process easier to reason about:

```text
open description
   │
   ├─ declared choices
   ├─ validated defaults/policies
   ▼
resolved request
   │
   └─ all semantically relevant choices fixed
   ▼
compiler planning and lowering
```

This is especially important for research/reproduction work. If the control
discretization, boundary convention, observation policy, or target realization can
silently change during compilation, two runs that look similar at the scenario level
may not represent the same numerical experiment.

## 5. Runtime data is separate from the semantic graph

A semantic graph can say that the problem has

- a forcing function;
- a desired state;
- a diffusion coefficient;
- fixed Dirichlet data;
- a transport field;
- a regularization parameter.

It should not need to own the actual deal.II `Function` objects used to evaluate
those quantities.

Runtime bindings connect the semantic ports to concrete numerical data.

For a Chapter 6 run, the outer application can therefore choose a target profile,
forcing definition, or transport field and bind it to the same semantic structure.
The compiler sees both the resolved semantic meaning and the concrete data required
to assemble operators.

This separation is useful for parameter studies: the graph may stay fixed while the
runtime values change.

## 6. The compiler is bounded and compositional

The current compiler is neither of the following extremes:

- it is **not** an unrestricted symbolic compiler able to lower arbitrary weak forms;
- it is **not** merely a switch over a few whole pre-written applications.

Instead, it composes registered pieces within a bounded capability set.

For example, the compiler can recognize a supported scalar PDE family, select the
realization of its terms, construct coordinates and operators, assemble observation
and metric objects, and package the requested formulation product.

That means two levels of reuse coexist:

1. **component reuse** – the same mass metric, coordinate map, observation mechanism,
   or solve service can participate in multiple problems;
2. **registered strategy reuse** – only combinations with a known coherent lowering
   path are accepted.

This is why the implementation has a capability ledger. The architecture is
compositional, but support remains explicit.

## 7. Numerical construction happens after semantic choices are fixed

Once a request is resolved and runtime data is present, the deal.II compiler creates
the numerical realization.

Depending on the problem, this may include:

- mesh and finite-element spaces;
- state/control coordinate maps;
- stiffness, mass, coupling, and observation operators;
- boundary lifting;
- state and adjoint solve services;
- optimization metrics;
- constraints;
- native output support.

These objects are ordinary typed numerical C++ objects. The semantic layer does not
replace them; it determines which ones should be constructed and how they should be
connected.

The [Numerical realization](numerical-realization.md) overview explains these
objects from the finite-element side.

## 8. What a compiled product contains

For the main reduced path, the useful result of compilation is not “a matrix.”
It is a coherent numerical problem containing the services needed by the formulation
and by the surrounding application.

Conceptually:

```text
compiled problem
├── solver-facing numerical services
│   ├── executable residual/objective operations
│   ├── state and adjoint solves
│   ├── metric
│   └── selected formulation product
│
├── typed compiled application view
│   ├── dimensions and coordinates
│   ├── application-native numerical objects
│   └── output/reconstruction support
│
└── compilation evidence
    ├── effective policies
    └── provenance / manifest information
```

The solver-facing part can be generic or type-erased where that simplifies the
optimization layer. The native view remains typed so that the application can still
do finite-element-specific work after compilation.

This is why a compiled problem should not be interpreted as the universal
application interface. It is the product of one producer path.

## 9. The compiler is formulation-plural

The semantic/compiler path is not a reduced-only front end. Its registered product
surface includes several formulations with different numerical solve structures.

For selected registered problem families, the compiler can produce distinct
formulation products:

```text
resolved semantic problem
           │
           ▼
        compiler
           │
     ┌─────┼──────────┬────────────┐
     ▼     ▼          ▼            ▼
 reduced   supplied   quadratic    complementarity /
 product   OTD        KKT          PDAS product
```

The supplied optimize-then-discretize (OTD), quadratic Karush–Kuhn–Tucker (KKT),
and primal-dual active-set (PDAS) surfaces have their own contract and deal.II
verification. They represent genuine framework capability.

The repository's executable benchmark programme is narrower than that product surface.
B1/B2 exercise reduced formulations extensively; PDAS, supplied-OTD, and KKT support
is instead evidenced mainly by focused compiler, contract, and numerical tests rather
than equivalent named end-to-end benchmark campaigns.

That distinction is useful when reading the project:

> **Compiler support tells you what the framework can construct; benchmark coverage
> tells you what has been exercised as a named reproduction/application campaign.**

Do not infer one from the other.

## 10. The semantic path does not own the optimizer

Compilation ends when the numerical formulation product is ready.

The compiler does not decide:

- whether the reduced solver uses steepest descent or L-BFGS;
- what Armijo constants are chosen;
- how many iterations are allowed;
- how a benchmark run directory is named.

Those are optimization or application/execution decisions.

This boundary matters because the same compiled reduced problem can be paired with
different solver policies without recompiling the PDE meaning from scratch.

## 11. How this differs from the external-application path

An external application can bypass the semantic/compiler layers entirely.

If an application already has the correct numerical realization, there is little
value in rebuilding it through `ProblemSpec`. The application can expose the needed
residual/objective actions, solves, and metric directly.

The two paths differ mainly in who constructs the numerical services:

```text
semantic ProblemSpec ──► validation / compiler ──┐
                                                 │
existing PDE application ──► thin binding ───────┤
                                                 ▼
                                   common numerical services
                                                 │
                              ┌──────────────────┼──────────────────┐
                              ▼                  ▼                  ▼
                           reduced           KKT / OTD           PDAS
                         formulation         products          formulation
```

From the common numerical boundary onward, formulation and algorithm code need not
know which producer path created the services. This does not mean every formulation
consumes one identical service interface: reduced, supplied-OTD, KKT, and PDAS
products have distinct formulation-specific actions and assumptions built on the
shared contract vocabulary.

This is an important practical choice for users. The semantic path is a convenience
and composition mechanism, not a compulsory front door.

## 12. Where the source code lives

The main implementation areas are:

- `include/nmopt/semantic/v1/` – semantic types, graph construction, validation,
  resolution;
- `include/nmopt/compiler/v1/` – compiler requests, planning/lowering, compiled
  products;
- `include/nmopt/dealii/` – reusable deal.II realization services;
- `include/nmopt/application/` – project recipes and scenarios that construct
  semantic problems for the Chapter 5/6 application layer.

The current operational capability description lives in the
[Compiler reference](../../reference/compiler.md). The corresponding
implementation stages, ownership, and verification map live in
[Compiler implementation](../../internals/compiler.md).

## 13. When to use this path

Use the semantic/compiler path when:

- the problem belongs to a reusable family that should be expressed independently of
  one executable application;
- you want spaces, PDE terms, observations, metrics, constraints, and formulation
  choices to be explicit and inspectable;
- you want `nmopt` to construct the supported deal.II realization;
- several scenarios should reuse the same problem-family builder.

Use a recipe when several application scenarios share that semantic structure.
Use a scenario when one concrete experiment needs to freeze data, discretization,
and solver choices. Use the experiment/replication layer only when you also need
run organization, manifests, benchmark evidence, and post-processing.

Prefer [direct external integration](external-applications.md) when a mature PDE code
already owns the numerical realization and expressing the entire application through
the compiler would add indirection without meaningful reuse.

The two entry paths solve different integration problems. Neither is a fallback for
the other.

## Read later: authoritative sources

For deeper detail, use the document whose authority matches the question:

- [Compiler reference](../../reference/compiler.md) records the current public
  capability and diagnostic boundary.
- [Compiler implementation](../../internals/compiler.md) records current
  lowering, realization, product, provenance, and verification mechanics.
- [Problem authoring](../../reference/problem-authoring.md),
  [Compiler](../../reference/compiler.md), and
  [Application authoring](../../reference/application-authoring.md) record the
  current programming interfaces around semantic problems, compilation, and
  reusable application composition.
- [Composition boundaries](../../design/composition-boundaries.md) records the
  long-lived rules governing what belongs in semantics, runtime bindings, compiler
  policy, and solver code.
- [Chapter 5 elliptic optimal-control guide](../../studies/chapter-5/source-catalogue.md)
  records the source problem catalogue drawn from the book.
- [Chapter 6 numerical-methods guide](../../studies/chapter-6/numerical-methods.md) and
  [numerical-examples reference](../../studies/chapter-6/numerical-examples.md) record
  the source methods and experiments behind the B-series application scenarios.

The overview explains the compiler's role; those documents own the exhaustive
capability and source records.
