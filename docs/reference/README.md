# Public API reference

This directory contains exact, agent-facing references for the public C++ API
and its configuration surfaces.

It will document:

- semantic `ProblemSpec` assembly and stable component identifiers;
- compiler data bindings, discretization policies, and compilation products;
- reduced, KKT, and PDAS solver options;
- experiment provenance and report records; and
- diagnostics and supported capability combinations.

The generic assembly workflow and current option defaults are documented in
the [application API reference](application-api.md). Concrete capability
combinations and Chapter 5/6 recipe choices will be added as separate
reference units; this document does not replace those application contracts.

The reference documents describe existing public contracts. They do not replace
the normative [interface specification](../design/interface-specification.md),
the [v1 compiler record](../implementation/v1/semantic-compiler.md), or the
[implementation roadmap](../planning/implementation-roadmap.md).
