# PDE–solver boundary refactor review

This directory records the architecture audit and bounded implementation plan
for the PDE/formulation/compiler cleanup.

Read the smallest document that answers the current task:

- [PDE–solver boundary](../../../design/pde-solver-boundary.md) — accepted
  long-lived ownership and integration rules.
- [Architecture map](architecture-map.md) — detailed current architecture,
  runtime/compiler flows, typed-numerics decomposition, ownership, and
  current-to-target correspondence.
- [Assessment](assessment.md) — audit evidence, finding IDs, judgments, and
  deferred decisions. Read only assigned findings during implementation.
- [Deletion ledger](deletion-ledger.md) — current types/paths to keep, narrow,
  migrate, audit, or delete and the gates that protect unique behavior.
- [Roadmap](roadmap.md) — executable work units, verification, decision
  criteria, and current handoff.

## Agent routing

For one implementation unit:

1. read the repository agent instructions routed for that action;
2. read the roadmap unit;
3. read only the boundary sections, architecture-map sections, assessment
   findings, and ledger entries named by that unit;
4. inspect the current non-test consumers and focused tests; and
5. prepare the unit execution brief required by the roadmap before editing.

Do not read the complete assessment by default and do not treat it as mutable
project status. If implementation evidence changes a planned boundary, update
the roadmap rather than rewriting historical audit evidence.

## Scope

This review owns PDE numerical realization, compiler composition, solver-facing
integration, directly required application/native-output migration, and
redundant legacy-path retirement.

It does not restart unrelated runner, GUI, benchmark-catalogue,
post-processing, or Chapter 5/6 feature work.
