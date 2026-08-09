# Workflow conventions

These conventions govern how repository work is scoped, sequenced, executed,
and handed off. [`git.md`](git.md) governs repository state and history;
action-specific convention files govern the content of each change.

## Work units and review boundaries

- Before substantial work, divide the objective into ordered units with one
  concrete outcome, an explicit boundary, and focused verification for each
  unit. A roadmap item or architectural gate may require several units.
- An assistant turn may complete several units when that is the most efficient
  use of time and context. Propose the independent review boundaries before
  substantial implementation and, when commits are requested, map them to
  separate commits according to [`git.md`](git.md). Do not spend a separate
  turn on a small follow-up unless it needs a user decision.
- Keep each unit self-contained and leave the repository consistent, buildable,
  and testable where practical. Include the tests and authoritative
  documentation needed to understand the behavior changed by the unit.
- If a unit grows beyond its stated boundary, stop and repartition it before
  absorbing adjacent cleanup, another feature, or the next roadmap item.
- Several units may proceed in one turn, but pause before a new architectural
  decision, a material expansion of scope or risk, or a review boundary
  requested by the user.
- When correctness requires tightly coupled changes that cannot be separated
  into honest intermediate states, explain the coupling, affected interfaces,
  expected size, and verification scope before implementation. A genuinely
  coupled unit may remain one detailed commit even when it is large. Still
  separate preparatory, mechanical, and documentation work when it remains
  independently meaningful.

## Efficient execution and context

- Read the routed authoritative documents and affected files needed for the
  current unit; do not scan unrelated files speculatively. Stop discovery once
  the relevant contract, boundary, and verification path are established.
- Batch independent read-only discovery and verification when practical. Reuse
  valid build results and run checks in proportion to the change, reserving
  broader suites for integration boundaries and final handoffs.
- Keep a compact running handoff of decisions, assumptions, completed commits,
  verification, and exact remaining work. Use it after context compaction or a
  resumed session instead of rediscovering unchanged repository state.
- Make progress updates decision-dense: report outcomes, changed assumptions,
  blockers, and the next boundary rather than narrating routine commands.
- Optimize for review and revert value rather than the number of turns or
  commits. Do not create artificial checkpoints that add cost without giving
  the user a meaningful opportunity to assess risk or change direction.
