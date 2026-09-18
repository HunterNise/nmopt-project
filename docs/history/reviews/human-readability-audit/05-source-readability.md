# Source readability audit

**Audit baseline:** [`f53b7f009e5c`](https://github.com/HunterNise/nmopt-project/commit/f53b7f009e5c418ec4f3855db29f7eb924faf6f1)  
**Status:** complete as a bounded readability pass for the documentation-refactor decision. This is not a dedicated repository-wide source-maintainability or refactor audit.

## 1. The code is not uniformly hard to read

A useful correction to the initial concern is that many low-level/public contract headers already contain good local explanations.

Examples include:

- [`contract/linear_solve.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/contract/linear_solve.hpp): explains why solve policy belongs to backend/compiler while formulations consume evidence;
- [`compiler/v1/native_application_view.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/compiler/v1/native_application_view.hpp): explains why native application access stays beside the erased solver view and states borrowing/lifetime assumptions;
- [`dealii/serial_backend.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/dealii/serial_backend.hpp): clearly identifies itself as a vector policy and excludes distributed ownership;
- [`dealii/independent_state_coordinates.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/dealii/independent_state_coordinates.hpp): explains the reconstruction `P` and fixed lifting role;
- [`dealii/mass_metric.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/dealii/mass_metric.hpp): states the SPD/Riesz interpretation and inverse-solve policy;
- [`application/recipe.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/application/recipe.hpp), [`scenario.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/application/scenario.hpp), and [`runner.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/application/runner.hpp): explicitly say which responsibilities they do not own.

So the final humanization work should not blanket the code with comments. It should preserve and extend the strongest existing pattern: explain mathematical meaning, ownership, lifetime, and non-obvious architectural boundaries.

## 2. Public source currently mixes enduring semantics with development-history vocabulary

This is the clearest source-documentation problem found so far.

Examples in public headers refer to internal project phases or benchmark labels such as:

```text
P5.1
P6.2
C5.6
Section 5.11.2
B1 / B2
selected target
selected product
historical policy
```

Representative files include:

- [`compiler/v1/dealii_types.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/compiler/v1/dealii_types.hpp);
- [`contract/supplied_otd.hpp`](https://github.com/HunterNise/nmopt-project/blob/f53b7f009e5c418ec4f3855db29f7eb924faf6f1/include/nmopt/contract/supplied_otd.hpp);
- application/Chapter-specific headers by design.

### Why this matters

There are two legitimate categories:

1. **Chapter/application code**, where `B1`, `B2`, Chapter 5/6 and source-study names are natural domain identifiers.
2. **Reusable public contracts/compiler infrastructure**, where comments should preferably describe the mathematical or software condition directly.

For example, a comment on a generic compiler binding is more durable if it says:

> “This binding is used only by general scalar tensor/transport/Robin realizations.”

than if understanding it requires finding what “P5.1” meant in a roadmap.

### Preliminary recommendation class

**Source comment cleanup**, not architectural change.

Historical labels may remain in audit/provenance documents and test names where traceability is the purpose.

## 3. Header-only organization exposes implementation scale directly to readers

The root CMake target `nmopt_contract` is an `INTERFACE` target. Much of the project is implemented in headers.

Several inspected headers are very large:

| File | Approximate size at baseline | Role |
| --- | ---: | --- |
| `semantic/v1/validation.hpp` | 178 KB | semantic validation |
| `compiler/v1/dealii_compiler.hpp` | 341 KB | compiler orchestration + many realization paths |
| `compiler/v1/dealii_fixed_dirichlet.hpp` | 96 KB | specialized deal.II realization |
| `compiler/v1/dealii_neumann_boundary.hpp` | 50 KB | specialized realization |
| `application/dealii/chapter6_b2.hpp` | 52 KB | project application adapter |
| `apps/nmopt-runner/main.cc` | 54 KB | CLI/run orchestration |
| `apps/nmopt-runner/parameter_files.hpp` | 53 KB | parameter/schema handling |

Large files are not automatically bad, especially for template-heavy C++ and an implementation that intentionally avoids a large compiled library.

But they affect human navigation:

- “go to the compiler header” is not a useful orientation instruction when that header is hundreds of KB;
- public implementation details and architectural entry points occupy the same physical files;
- GitHub browsing becomes less effective;
- responsibility boundaries that are clear conceptually may be hard to see locally.

### Audit question

The later recommendation must decide among:

- overview/source maps only;
- stronger section-level comments;
- internal helper/header extraction;
- moving implementation detail to `detail/`;
- introducing `.cc` compilation units where templates do not require headers;
- leaving code unchanged because source movement risk outweighs readability gain.

No file split is recommended yet.

## 4. Some namespace/directory names are broader or narrower than their contents

### `dealii_backend`

`SerialBackend` is genuinely a backend vector/algebra policy.

The same namespace/directory also contains:

- FE independent-state coordinates;
- sparse metrics;
- box constraints;
- serial SPD/KKT/PDAS solve services.

To a human, “backend” could incorrectly imply that all of these are primitive storage operations.

The architecture audit currently distinguishes:

```text
backend representation policy
vs
deal.II numerical realization/services
```

The final docs may solve this without renaming code. But the name/content mismatch should be acknowledged.

### `application`

The namespace contains both:

```text
generic recipe/scenario/runner abstractions
and
Chapter 5/6 project-specific clients/adapters
```

This weakens the ability to infer reusable scope from source placement.

Again, documentation may be enough; later phases will decide.

## 5. Umbrella headers optimize convenience more than conceptual locality

### `semantic/v1/problem_spec.hpp`

This aggregate includes:

- core semantic types;
- reference specs;
- resolution;
- validation.

The file comments acknowledge that new code may include smaller pieces.

Potential human issue:

A reader looking for the definition of `ProblemSpec` through the obvious filename is routed into an aggregate whose dependencies include large reference/example material. The actual type definition lives in `types.hpp`.

### `application/application.hpp`

This aggregate includes generic application facilities and Chapter-specific content.

That is convenient for `nmopt-runner` but less useful as a “public library entry point” for an external human user.

### Preliminary recommendation class

**API discoverability / include-surface documentation** first. Consider physical restructuring only if later inspection shows real maintenance cost.

## 6. Lifetime and ownership are often handled carefully, but the pattern is sophisticated

Several interfaces have nontrivial lifetime behavior:

- `ReducedDTOT` supports borrowed and owned executable models/lifetime tokens;
- retained reduced values carry private tokens to prevent misuse across evaluations;
- compiled products retain a lifetime owner;
- `NativeApplicationViewT` can borrow compile-time data bindings even while owning a retained model;
- supplied-OTD products deliberately declare owner members before callback members so reverse destruction order is safe;
- external bindings rely on object construction/member order.

Positive finding:

The code often comments these rules exactly where they matter.

Human risk:

A user may still need architecture/reference documentation before they know **which** lifetime model applies to their use case.

### Likely documentation need

One concise ownership/lifetime diagram per major producer path is likely more useful than adding comments everywhere.

## 7. Mathematical naming is generally strong

The inspected core uses terms such as:

```text
Primal
Covector
Pairing
Metric
Residual JVP
Residual VJP
StateAdjointSolvers
ReducedHessian
QuadraticKKT
Complementarity
```

These names communicate mathematical roles better than generic software names such as “data,” “handler,” or “manager.”

This is worth preserving.

A humanization refactor should **not** simplify mathematically precise language merely to make the API look smaller. It should explain the terms once and use them consistently.

## 8. Compiler naming becomes more implementation-oriented

Within the compiler, names such as:

```text
ResolvedCompilationRequest
ResolvedTargetFamily
ScalarLoweringPlan
ScalarResidualAssemblyPlan
DealiiDataBindings
DealiiCompilationSession
CompilationManifest
```

are understandable individually, but the relationship among them is difficult to infer without following the compile function.

This is primarily an **explanatory/source-map problem** at present.

A short compiler source map could explain:

```text
ProblemSpec
 -> ResolvedProblemView
 -> ResolvedCompilationRequest
 -> ScalarLoweringPlan / specialized registration
 -> numerical realization
 -> Compiled* product + manifest
```

That may remove more human friction than renaming these types.

## 9. Project application code is appropriately more concrete

The Step-4 integration is comparatively readable because responsibilities are physically separated:

```text
adapted Step-4
ProblemB
ProblemB coordinates
ProblemB mass/coupling
ProblemB native metric
minimal nmopt binding
minimal executable
```

This organization is a useful example for the rest of the repo:

- files correspond to human-recognizable responsibilities;
- application mathematics can be read without nmopt;
- nmopt binding can be read without evaluation instrumentation;
- evidence machinery lives elsewhere.

The final source-readability recommendations should ask where the main library could benefit from similar responsibility locality.

## 10. Runner/application implementation is a likely readability hotspot

`apps/nmopt-runner` is already split into several helper headers, but:

- `main.cc` is still ~54 KB;
- parameter-file machinery is ~53 KB;
- Chapter 6 deal.II adapters are large;
- runner code necessarily combines CLI, scenario selection, provenance, filesystem layout, and benchmark-specific dispatch.

This may be acceptable because it is **application code**, not a reusable numerical core.

The key humanization question is therefore not automatically “make the runner small.” It is:

> Can a human understand that this complexity belongs to the research execution shell and is not required to use nmopt as a numerical library?

If the documentation makes that boundary clear, aggressive runner refactoring may be unnecessary.

## 11. Comment/docstring strategy

A likely final policy, subject to full audit, is:

### Public reusable types

Document when one of these is non-obvious:

- mathematical role;
- ownership;
- lifetime/borrowing;
- invariants;
- representation conventions;
- why a capability exists;
- important exclusions.

### Internal implementation

Comment:

- non-obvious mathematical transformations;
- policy reasons;
- subtle lifetime/cache behavior;
- why a branch exists.

Do not comment:

- ordinary C++ syntax;
- obvious field assignments;
- roadmap history that does not explain current behavior.

### Doxygen/API generation

Do not add generated API documentation in this refactor.

The inspected files mainly use ordinary `//` comments rather than a systematic
generated-API documentation style. The most important missing explanations are
architectural, mathematical, and ownership-oriented rather than signature
listings. A generated surface would add another authority layer and would
require a substantial annotation pass before it became useful.

Reconsider Doxygen or another generator only if nmopt later becomes a
separately packaged library whose users need exhaustive symbol browsing.

## 12. Intervention classification

| Finding | Current intervention level |
| --- | --- |
| Good mathematical/lifetime comments in core | preserve |
| Development-phase vocabulary in reusable comments | cheap candidate, but defer from the documentation-only branch |
| Huge header-only compiler/validation files | document/source-map first; deeper restructuring requires a dedicated source audit |
| `dealii_backend` concept broader in directory than vector backend meaning | docs terminology first |
| generic + Chapter-specific `application/` content | clarify information architecture; source reorganization is not justified by this audit |
| aggregate headers obscure minimal include surface | docs/API discoverability first |
| external binding ceremony | future API ergonomics audit; no redesign decision |
| runner size | likely project-shell complexity; avoid refactoring merely for LOC |


## 13. Boundary of the source-refactor conclusion

This pass was intentionally coupled to architecture reconstruction: it looked
for source-level problems that were obvious enough, cheap enough, or
architecturally important enough to affect the human-readability refactor.

It supports the following bounded conclusion:

> No broad code, API, namespace, or file-organization refactor is justified as
> part of the documentation phase from the evidence inspected here.

It does **not** establish that deeper restructuring could never be beneficial.
Systematic header decomposition, include-cost reduction, API simplification,
namespace/directory reorganization, and long-term maintenance hotspots would
require a separate source-refactor audit with evidence such as dependency
structure, compile cost, change frequency, review friction, and repeated
maintenance patterns.

For `codex/docs-refactor`, implementation code should therefore remain frozen
unless a correctness issue is discovered. The small comment-only terminology
cleanup identified above is deliberately deferred so that the branch remains
strictly documentation-only.
