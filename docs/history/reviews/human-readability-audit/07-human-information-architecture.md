# Proposed human information architecture

**Audit baseline:** [`f53b7f009e5c`](https://github.com/HunterNise/nmopt-project/commit/f53b7f009e5c418ec4f3855db29f7eb924faf6f1)  
**Status:** recommended documentation structure for the refactor, based on the code-derived architecture and documentation audit. Exact wording and routing remain subject to review as the documentation commits are prepared.

## 1. Design goal

A technically competent reader should be able to answer, in order:

1. What is nmopt?
2. What is the reusable numerical core?
3. How can mathematics reach that core?
4. What is deal.II's role?
5. What is optional project/research infrastructure?
6. Where do I go for exact APIs or historical evidence?

without opening a roadmap or an implementation review.

## 2. Canonical mental model

The first architecture page should center on the convergence already present in the code:

```text
                         semantic authoring path
ProblemSpec
    |
    v
validation / resolution
    |
    v
deal.II compiler + numerical realization
    |
    +-----------------------------------+
                                        |
existing numerical application          |
    |                                   |
    +--> callback/solve/metric binding -+
                                        |
                                        v
                         mathematical/formulation boundary
                                        |
             +--------------------------+----------------------+
             |                                                 |
             v                                                 v
       ReducedDTOT                                  additional products
             |                                  supplied OTD / KKT / PDAS
             v
   reduced search / trust region
             |
             v
          result
```

Underneath the numerical values:

```text
Backend policy = vector/storage algebra
```

Around selected project runs:

```text
recipes/scenarios -> runner -> artifacts -> post-processing
```

The outer application shell is a client of the semantic path, not a mandatory layer between compiler and solver.

## 3. Add `docs/overview/`

Recommended files:

```text
docs/overview/
  project-architecture.md
  semantic-compiler.md
  numerical-realization.md
  reduced-optimization.md
  research-execution.md
```

Five files are justified by the current repository. Fewer would force the same “one giant blueprint” problem to recur.

## 4. `project-architecture.md`

### Purpose

Primary human entry point. Roughly 6–10 pages, not an exhaustive API document.

### It should answer

- the mathematical problem family;
- project scope and tested limitations;
- core versus producer/client infrastructure;
- two producer paths;
- formulation/optimizer convergence;
- backend versus numerical realization;
- additional KKT/OTD/PDAS capabilities;
- verification architecture;
- source map.

### Recommended sections

```text
1. What nmopt is
2. The mathematical center
3. The architecture in one diagram
4. Two ways to supply numerical mathematics
5. The common numerical/formulation boundary
6. Reduced optimization as the primary runtime path
7. Additional implemented formulation products
8. deal.II and numerical realization
9. Optional research/application execution shell
10. Verification and evidence
11. Source map
12. Where to read next
```

### Scope statement

A more accurate current description is:

> nmopt is a C++17 collection of backend-parametric mathematical contracts, formulations, and optimization algorithms for discretized PDE-constrained optimization. Supported problems may be realized by its deal.II semantic/compiler path or supplied by an existing numerical application through callbacks and solve services. The repository also contains bounded all-at-once/KKT/active-set products and a research execution layer for Chapter 5/6 scenarios.

This statement makes deal.II important without making the compiler/runner universal.

## 5. `semantic-compiler.md`

### Purpose

Explain the authoring/compiler producer path without duplicating the v1 capability table.

### Key distinctions

```text
ProblemSpec
    backend-neutral graph + explicit non-inferable realization policies

SemanticValidator / Resolver
    graph closure and stated-policy consistency

ResolvedCompilationRequest
    compiler-owned closed interpretation

ScalarLoweringPlan / specialized registration
    bounded component planning, not arbitrary combination

Concrete realization
    FE spaces, coordinates, operators, solves, metrics, constraints

Compiled* product
    solver-facing ports + typed provenance + optional native view
```

### Must say explicitly

- semantic graph != concrete deal.II objects;
- semantic graph != purely continuous mathematics;
- component planning exists;
- arbitrary independent component recombination is **not** implemented;
- compiler is optional for an existing application;
- `CompiledProblemT` is not a universal solver integration type.

### Deep links

- normative semantics: `design/interface-specification.md`;
- current boundary: `design/pde-solver-boundary.md`;
- exact capabilities: `implementation/v1/semantic-compiler.md`;
- exact application assembly API: `reference/application-api.md`.

## 6. `numerical-realization.md`

### Purpose

Explain the layer that is easiest to misunderstand from names alone.

### Core message

```text
backend policy != finite-element realization != compiler
```

### Explain

`SerialBackend`:

- native vector type;
- primitive algebra.

deal.II numerical services:

- state coordinates;
- reconstruction and pullback;
- mass/negative/trace metrics;
- constraints;
- state/adjoint/KKT solve services;
- assembled operators.

compiler:

- decides and composes supported numerical realizations.

external application:

- may already own equivalent realization operations.

This page should use small examples such as:

```text
y_phys = P y_hat + ell
G u = M u
G^-1 xi = solve(M g = xi)
```

and link to exact backend/compiler code.

## 7. `reduced-optimization.md`

### Purpose

Explain the computational center shared by semantic and external paths.

### Main runtime

```text
u
 -> solve state
 -> objective
 -> objective derivative
 -> solve adjoint
 -> residual pullback
 -> reduced derivative
 -> metric inverse
 -> search direction
 -> globalization
 -> stopping
```

Explain the value/derivative split:

```text
rejected trial:
  state solve + objective

accepted trial:
  augment with adjoint + derivative
```

### Ownership

| Concern | Owner |
| --- | --- |
| PDE/objective numerical operations | producer/application |
| when state/adjoint operations are needed | formulation |
| primal-dual search geometry | metric |
| search direction/globalization/stopping | optimizer |
| native field output | application/native view |

### Limitations

State explicitly:

- current reduced DTO is one state block + one decision block + one test block;
- full `ExecutableModelT` requires residual/JVP even though selected first-order reduced execution does not call them;
- serial deal.II backend is the tested backend;
- additional KKT/PDAS products are separate, not hidden reduced modes.

## 8. `research-execution.md`

### Purpose

Prevent the size of `application/`, `apps/nmopt-runner`, parameters, artifacts, and tools from making the numerical library look more coupled than it is.

### Explain

```text
ProblemRecipe
    parameters -> ProblemSpec

Scenario
    problem + compile + solver + experiment choices

execution adapter
    scenario + runtime data -> compiler + solver -> evidence

HeadlessBenchmarkRunner
    orchestration and artifact finalization

nmopt_runner
    CLI, matrix expansion, run-set directories, manifests

tools/
    persisted-evidence post-processing/reporting
```

### Explicit statement

> None of this shell is required by the external Step-4 library-use path.

### Deep links

- `reference/application-execution.md`;
- `reference/parameter-files.md`;
- Chapter 6 application/benchmark docs.

## 9. No separate overview for KKT/PDAS initially

The audit found that these capabilities are real:

- dedicated neutral contract tests;
- compiled KKT/PDAS deal.II integration scenarios.

But they are not the primary application runner path.

Therefore:

- introduce them in `project-architecture.md`;
- reference exact numerical-contract documentation;
- add a dedicated overview later only if they become a primary user workflow.

## 10. Rewrite `docs/README.md` as a router

Suggested shape:

```text
# Documentation

## Understand the project
- Project architecture
- Semantic/compiler path
- Numerical realization
- Reduced optimization
- Research execution

## Use nmopt
- Existing deal.II application integration
- Semantic/compiler application API
- Current compiler capability table
- Numerical contracts

## Reproduce project applications
- Chapter 6 execution
- Parameter files
- benchmark contracts

## Design and mathematical foundations
- theoretical formalism
- interface specification
- PDE/solver boundary
- composition boundaries

## Evidence and development history
- external Step-4 report/closure
- implementation/history records
- planning/review archive

## Coding-agent workflow
- AGENTS.md / .agents/
```

No giant “choose by task” table is necessary once the main paths are obvious. A short table for common tasks is fine.

## 11. Rewrite root `README.md` routing

Keep:

- concise project description;
- build/test;
- basic application execution;
- repository layout.

Change:

### Scope

Replace active-development framing with:

- implemented scope;
- tested boundaries;
- explicitly deferred work.

### Where to start

Route:

1. project architecture;
2. choose either external application or semantic/compiler path;
3. exact API/reference only afterward.

### Agent material

Keep a short final section, but do not make agent workflow part of the project identity.

## 12. Current versus historical material

Every major document should be classifiable without reading several pages.

Recommended labels:

```text
Overview
Design authority
Reference
Implementation record
Guide
Case study
Review/evidence
Historical roadmap
Draft
```

Do not rewrite old reviews to sound current.

Instead add a short status banner where needed:

> Historical development record. Preserved for rationale/evidence; not current API or project-status authority.

## 13. Source map for the primary overview

Recommended table:

| Responsibility | Primary source |
| --- | --- |
| Typed mathematical values/layouts | `include/nmopt/contract/` |
| Formulations / KKT / complementarity | `include/nmopt/contract/` |
| Reduced algorithms | `include/nmopt/solvers/` |
| Semantic graph/validation | `include/nmopt/semantic/v1/` |
| Compiler and compiled products | `include/nmopt/compiler/v1/` |
| deal.II numerical services | `include/nmopt/dealii/` |
| Generic application/evidence shell | `include/nmopt/application/`, `include/nmopt/experiment/` |
| Headless project runner | `apps/nmopt-runner/` |
| External integration example | `apps/external-dealii/step-4/` |
| Dense/reference oracles | `include/nmopt/reference/` |
| Verification | `tests/` |
| Post-processing/reproduction | `tools/`, `parameters/` |

## 14. Human reading paths

### “I want to understand the project”

```text
README
 -> overview/project-architecture
 -> overview/reduced-optimization
 -> design/pde-solver-boundary
```

### “I already have a deal.II application”

```text
overview/project-architecture
 -> Step-4 external overview
 -> external integration reference
 -> minimal A/B consumers
```

### “I want nmopt to build a supported problem”

```text
overview/project-architecture
 -> overview/semantic-compiler
 -> reference/application-api
 -> implementation/v1/semantic-compiler capability table
```

### “I want to modify algorithms”

```text
overview/reduced-optimization
 -> reference/numerical-contracts
 -> include/nmopt/contract + include/nmopt/solvers
```

### “I want to reproduce the research runs”

```text
overview/research-execution
 -> reference/application-execution
 -> reference/parameter-files
 -> benchmark/application contract
```

## 15. What not to do

Do not:

- make one enormous replacement blueprint;
- duplicate the v1 capability table into overview docs;
- move all review/history files merely for aesthetics;
- introduce generated API docs as a second source of authority;
- describe the semantic/compiler route as mandatory;
- call the deal.II FE realization simply “the backend”;
- make the Chapter 6 runner define the reusable library boundary.
