# Application recipes

Application documents describe how to assemble concrete problems from the
public API. They are intentionally separate from the mathematical catalogue
in `docs/guides/` and from mutable implementation status in
`docs/planning/`.

Chapter 5 documents will define reusable, parameterized problem recipes.
Chapter 6 documents will define application-level method/scenario choices
that consume those recipes.

A recipe returns or constructs a semantic `ProblemSpec`. Compilation data,
discretization choices, solver policies, and experiment metadata remain
separate typed inputs.
