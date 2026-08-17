# Application recipes

Application documents describe how to assemble concrete problems from the
public API. They are intentionally separate from the mathematical catalogue
in `docs/guides/` and from mutable implementation status in
`docs/planning/`.

Chapter 5 documents will define reusable, parameterized problem recipes.
Chapter 6 documents will define application-level method/scenario choices
that consume those recipes.

Start here for the concrete authoring contracts:

- [Chapter 5 recipe contracts](chapter-5/README.md) route a mathematical
  variant to a registered semantic graph and its required bindings.
- [Chapter 6 scenario contracts](chapter-6/README.md) freeze the selected B0,
  B1, and B2 application runs and identify the later benchmark extensions.

A recipe returns or constructs a semantic `ProblemSpec`. Compilation data,
discretization choices, solver policies, and experiment metadata remain
separate typed inputs.
