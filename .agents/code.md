# Code conventions

## Layer boundaries

- Keep semantic problem descriptions, executable contracts, discretization
  backends, and solver/formulation code in separate layers.
- Connect layers through explicit, stable interfaces or ports. Do not bypass
  a layer by adding backend or solver details to a semantic component.
- Do not create subclasses for particular combinations of PDE, objective,
  control, or boundary condition.
- Do not add domain-specific branches to generic solvers or solver-specific
  branches to reusable components.
- Keep implementation details in the narrowest layer that needs them.

## Changes and tests

- Do not invoke `clang-format` or another automated source formatter unless the
  user explicitly approves the exact command and scope for the current work.
  The same approval requirement applies to linter auto-fixes, code generators,
  and bulk-rewrite commands that can modify existing source. Ordinary
  compiler/test diagnostics and read-only search or diff inspection do not
  require this approval.
- Before changing code, identify the affected layer and read its authoritative
  contract or implementation document from [`docs/`](../docs/README.md).
- When a public interface or file boundary changes, update its authoritative
  documentation and add or update focused tests.
- Keep tests close to the contract they verify and cover changed public
  behavior, including forward and reverse actions where applicable.
- Prefer small, coherent changes. Do not broaden an interface merely to avoid
  making a deliberate architectural decision.
- Use compatibility seams or staged migrations when they create honest,
  independently testable intermediate states. Do not add a permanent public
  interface solely to manufacture an artificial split.
