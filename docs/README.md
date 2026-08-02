# Documentation map

This directory is organized by purpose and audience. Start with the smallest
document that answers the question, then follow links to the deeper material.

## Authority and status

- `interface-specification.md` is the normative, implementation-neutral
  contract for semantic interfaces and component protocols.
- `implementation-readiness-review.md` selects the required default policies
  for the first executable implementation when the specification allows
  multiple choices.
- `executable-contract-v0.md` and `dealii-v0-lowerer.md` describe the current
  implemented contracts, supported slice, and explicit exclusions.
- `semantic-v1-compiler.md` describes the implemented v1 semantic graph and
  the separately owned compiler path that is compared with the v0 reference.
- `architecture.md` records the long-term design rationale and boundaries.
- `implementation-roadmap.md` records current implementation state, task
  order, acceptance checks, and agent handoff guidance.
- Project-wide working conventions live in the
  [project conventions](../conventions/README.md); read the applicable
  convention before changing files.
- The remaining documents provide mathematical background, examples, and
  implementation analysis; they do not silently override the normative
  contracts above.

## Choose by task

| Audience or task | Start with | Then consult |
| --- | --- | --- |
| New contributor or agent | `system-blueprint.md` | `../AGENTS.md`, `implementation-roadmap.md`, `interface-specification.md` |
| Understand the whole system and its code correspondence | `system-blueprint.md` | `interface-specification.md`, `executable-contract-v0.md` |
| Changing semantic interfaces or ports | `interface-specification.md` | `architecture.md`, `implementation-readiness-review.md` |
| Changing the backend-neutral executable API | `executable-contract-v0.md` | `implementation-readiness-review.md`, `implementation-roadmap.md` |
| Changing deal.II code or the compiler/lowerer | `dealii-v0-lowerer.md` | `executable-contract-v0.md`, `implementation-roadmap.md` |
| Changing the v1 semantic graph/compiler | `semantic-v1-compiler.md` | `interface-specification.md`, `implementation-readiness-review.md` |
| Deciding component ownership or layer boundaries | `composition-boundaries.md` | `architecture.md`, `interface-specification.md` |
| Checking mathematical signs or formulas | `theoretical-formalism.md` | `laplace-interface-formulas.md`, `laplace-growth-case-study.md` |
| Understanding the concrete end-to-end example | `laplace-growth-case-study.md` | `laplace-interface-formulas.md` |
| Choosing the next implementation task | `implementation-roadmap.md` | The task-specific contract listed above |
| Reviewing why the design is shaped this way | `architecture.md` | `theoretical-formalism.md` |
| Editing Markdown or LaTeX | `../conventions/documentation.md` | The document being changed |

## Document index

- [Architecture record](architecture.md)
- [System blueprint](system-blueprint.md)
- [Composition boundaries](composition-boundaries.md)
- [Theoretical formalism](theoretical-formalism.md)
- [Interface specification](interface-specification.md)
- [Implementation-readiness review](implementation-readiness-review.md)
- [V0 executable contract](executable-contract-v0.md)
- [deal.II v0 lowerer](dealii-v0-lowerer.md)
- [V1 semantic graph and compiler](semantic-v1-compiler.md)
- [Laplace growth case study](laplace-growth-case-study.md)
- [Laplace interface formulas](laplace-interface-formulas.md)
- [Implementation roadmap](implementation-roadmap.md)
- [Project conventions](../conventions/README.md)
