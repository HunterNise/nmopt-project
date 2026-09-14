# nmopt human-readability and architecture audit

**Audit baseline:** [`f53b7f009e5c`](https://github.com/HunterNise/nmopt-project/commit/f53b7f009e5c418ec4f3855db29f7eb924faf6f1) (`codex/main` at the start of the audit)  
**Baseline commit message:** `docs(dealii): close the external boundary evaluation`  
**Status:** closed for the documentation-refactor decision. Refactor outcomes will be recorded in a later `closure-report.md`.

## Purpose

This audit reconstructs the current architecture from implementation evidence and
then asks how the repository should be presented to human readers, users,
reviewers, and maintainers.

Existing documentation is treated as evidence to reconcile rather than assumed
to be current. Public headers, build dependencies, tests, representative runtime
paths, and implemented application/compiler behavior are the primary sources of
truth.

## Audit package

The audit contains ten documents: this index, eight focused audit reports, and
the original audit plan.

| File | Purpose | Final state |
| --- | --- | --- |
| [`audit-plan.md`](audit-plan.md) | Original methodology, scope, and stopping criteria | historical plan |
| `00-audit-index.md` | Audit entry point and final synthesis | complete |
| [`01-repository-scope.md`](01-repository-scope.md) | Reconstruct current project scope and classify core versus optional/project infrastructure | complete |
| [`02-architecture-map.md`](02-architecture-map.md) | Reconstruct responsibilities, ownership, and dependency boundaries from code | complete |
| [`03-runtime-paths.md`](03-runtime-paths.md) | Trace representative external, semantic/compiler, and research-execution paths | complete for the architecture decision |
| [`04-public-api-and-user-surfaces.md`](04-public-api-and-user-surfaces.md) | Separate the major user surfaces and record API ergonomics observations | complete for the documentation refactor |
| [`05-source-readability.md`](05-source-readability.md) | Assess source-level readability and identify whether cheap code cleanup is justified | complete as a bounded readability pass |
| [`06-documentation-audit.md`](06-documentation-audit.md) | Reconcile current, stale, historical, reference, and planning documentation | complete for the refactor scope |
| [`07-human-information-architecture.md`](07-human-information-architecture.md) | Propose the human-facing documentation/navigation structure | recommended structure |
| [`08-recommendations.md`](08-recommendations.md) | Define the documentation-refactor sequence and explicit deferrals | final handoff |

## Final architectural synthesis

The implementation supports a compact current mental model:

```text
semantic ProblemSpec -> compiler/numerical realization --+
                                                         |
existing numerical application -> adapters/solve ports --+--> shared formulation contracts
                                                                  |
                                                                  +--> primary reduced optimization
                                                                  |
                                                                  +--> secondary OTD/KKT/PDAS products

research/application runner -------------------------------------> outer execution/provenance shell
```

The important conclusions are:

1. The semantic/compiler path and the external-application path are independent
   producers that converge before optimization.
2. `CompiledProblemT` is a compiler product, not the universal integration type
   for all users.
3. The primary shared runtime is the state-adjoint reduced formulation plus
   metric and reduced solver.
4. Supplied OTD, quadratic KKT, complementarity, and PDAS are real implemented
   secondary capabilities and should remain visible without dominating
   onboarding.
5. Backend vector/storage policy is narrower than FE/numerical realization.
6. Typed numerical ownership can remain beside solver-facing type erasure;
   native output and application dimensions need not be forced through the
   erased solver API.
7. The application/runner layer is an optional research/reproduction shell,
   not a prerequisite for using the numerical library.

## Main human-readability finding

The repository's dominant readability problem is not missing documentation or
an obviously broken architecture. It is **inspectability and authority
presentation**:

- no small canonical current overview layer;
- root/documentation routing that emphasizes implementation history and agent
  workflow;
- current API material still routed through historical `v0` documentation;
- stale status and repository-layout claims;
- overlapping design, implementation, planning, and reference documents whose
  authority is not obvious to a new reader.

The recommended phase is therefore documentation-only.

## Boundary of the code-refactor conclusion

The source-readability pass inspected representative public contracts,
semantic/compiler code, deal.II realization services, application/runner code,
large header-only units, naming, ownership/lifetime explanations, and external
integration ergonomics.

That evidence is sufficient to say:

> No obvious or inexpensive code, API, namespace, or file-organization refactor
> is justified as part of this human-readability phase.

This is deliberately **not** a claim that deeper source restructuring could
never be useful. A future decision to decompose large headers, simplify public
APIs, reorganize namespaces/directories, or optimize include/build structure
would require a dedicated source-refactor audit using maintenance, dependency,
compile-cost, edit-frequency, and reviewer-workflow evidence.

Even the small comment-only cleanup candidates found during this audit are
deferred so that `codex/docs-refactor` remains strictly documentation-only.

## Refactor handoff

The recommended order is:

1. preserve this audit under `docs/planning/review/`;
2. add a small canonical `docs/overview/` layer;
3. add a current numerical-contract reference;
4. reroute root and documentation entry points;
5. reconcile stale/current/historical authority markers;
6. run a documentation-wide consistency and link sweep;
7. add `closure-report.md` recording what was actually changed, deferred, and
   verified.

Detailed recommendations are in [`08-recommendations.md`](08-recommendations.md).

## Evidence convention

Claims in the focused audit files are tied where practical to exact source,
test, build, or documentation paths at the baseline commit. Baseline-pinned
links are intentional: these files are review evidence, not mutable current
architecture documentation.
