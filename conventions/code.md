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

- Before changing code, identify the affected layer and read its authoritative
  contract or implementation document from [`docs/`](../docs/README.md).
- When a public interface or file boundary changes, update its authoritative
  documentation and add or update focused tests.
- Keep tests close to the contract they verify and cover changed public
  behavior, including forward and reverse actions where applicable.
- Prefer small, coherent changes. Do not broaden an interface merely to avoid
  making a deliberate architectural decision.
