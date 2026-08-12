# Pre-Chapter 5/6 review routing

This directory contains the operational guidance and exhaustive evidence for
the pre-Chapter 5/6 review. An agent executing one batch should not read the
full assessment by default.

## Authority and status ownership

- [Stage B roadmap](stage-b-roadmap.md) defines static batch scope,
  prerequisites, exit checks, and feature gates.
- [Assessment](assessment.md) owns the evidence and tradeoffs behind findings
  `RF-001` through `RF-020`.
- [Assessment plan](assessment-plan.md) records the completed Stage A audit
  method. It is historical context, not required execution reading.
- [Implementation roadmap](../../implementation-roadmap.md) remains the sole
  owner of mutable project status, completed tasks, and the current handoff.

The integration branch is `codex/ch5-ch6-development`. The immutable
behavioral comparison point is tag `pre-refactor-ch5-ch6` at `7c2496b`.

## Reading protocol for one batch

Every agent must:

1. Read root `AGENTS.md`, `conventions/README.md`, and every convention routed
   for the planned actions.
2. Read the global rules and only the assigned batch in the
   [Stage B roadmap](stage-b-roadmap.md).
3. Read every finding assigned to that batch in the
   [assessment](assessment.md), including its evidence, verification, and
   tradeoff.
4. Use the root [documentation map](../../../README.md) to select the authoritative
   contract or guide for the affected layer.
5. Inspect the affected implementation and tests directly. Do not rely on the
   assessment as a substitute for current code.

Read adjacent findings only when the batch declares a dependency, verification
reveals a cross-layer interaction, or the proposed change would alter a
documented boundary. Read the full assessment only for a new architecture or
scope decision.

## Batch reading map

| Batch | Findings | Primary documents after this file |
| --- | --- | --- |
| R0 | RF-020 | [Implementation roadmap](../../implementation-roadmap.md), [v1 compiler record](../../../implementation/v1/semantic-compiler.md), and status-bearing documents named by RF-020 |
| R1 | RF-006, RF-009, RF-016 | Affected executable/semantic/compiler contracts and current tests |
| R2a | RF-001, RF-006 | [V0 executable contract](../../../implementation/v0/executable-contract.md) and typed contract tests |
| R2b | RF-002 through RF-006 | [Interface specification](../../../design/interface-specification.md), [v1 compiler record](../../../implementation/v1/semantic-compiler.md), and semantic tests |
| R2c | RF-008, RF-012 | Selected policies, compiler record, and deal.II contract tests |
| R3 | RF-016 through RF-019 | Build conventions, root CMake configuration, and test registration |
| C1 | RF-008 through RF-013 | V0/v1 contracts, selected policies, compiler code, and numerical provenance requirements |
| C2 | RF-007, RF-009 | [Composition boundaries](../../../design/composition-boundaries.md), compiler record, selected feature guide, and current lowerers |
| S1 | RF-014, RF-015 | [V0 executable contract](../../../implementation/v0/executable-contract.md), Chapter 6 P6.1 sections, solver code, and contract tests |

## Handoff requirements

A completed batch handoff states:

- batch and finding IDs;
- behavior preserved and deliberate contract changes;
- files changed;
- focused and full verification performed;
- comparison with the tagged baseline when numerical behavior could drift;
- documentation/status updates made in their authoritative owner;
- the next conditional gate, without silently selecting a feature.
