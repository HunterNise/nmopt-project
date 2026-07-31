# Project guide

## Mission

Build a deal.II-based framework for PDE-constrained optimal control and
inverse problems. The project should make combinations of PDE operators,
controls, observations, norms, constraints, and solvers reusable without
creating a new problem class for every combination.

The documentation map is [docs/README.md](docs/README.md). Use it to select
the documents relevant to a task; do not assume that every document is needed.

## Code and file architecture

- Keep semantic problem descriptions, executable contracts, discretization
  backends, and solver/formulation code in separate layers.
- Connect layers through explicit, stable interfaces or ports. Do not bypass
  a layer by adding backend or solver details to a semantic component.
- Do not create subclasses for particular combinations of PDE, objective,
  control, or boundary condition.
- Do not add PDE-specific branches to generic solvers or solver-specific
  branches to PDE components.
- Keep implementation details in the narrowest layer that needs them.
- When a public interface or file boundary changes, update its authoritative
  documentation and add or update focused tests.
- Use the implementation roadmap for task scope and sequencing. If a task
  requires a new architectural decision, document the decision instead of
  silently extending an existing interface.

## Working rules

- Inspect `git status` before changing files and preserve existing user work.
- Work on the currently selected branch or worktree; do not switch branches
  for convenience.
- Do not merge, push, force-push, or rewrite history unless explicitly asked.
- Prefer small, coherent changes with a clear handoff.
- Do not commit build products, generated files, editor files, or secrets.

## Validation

Run the baseline checks after implementation changes:

```bash
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

If an optional backend is unavailable, run the tests that can be built and
report which checks were skipped.

## Documentation and Markdown style

- Use `$...$` for inline LaTeX and `$$...$$` for display equations. Do not use
  `\(...\)` or `\[...\]`, because the project Markdown viewer does not render
  those delimiters reliably.
- Always brace superscripts and subscripts, including one-character scripts:
  write `x_{i}`, `x^{*}`, and `E'(x)^{*}` rather than `x_i`, `x^*`, or
  `E'(x)^*`. This avoids parser ambiguity around primes and nested scripts.
- For inline formulas containing underscores or other Markdown punctuation, use
  GitHub's backtick-delimited math form:

  ```text
  $`x_{i}`$
  ```

- Do not use `\operatorname`, `\tag`, `\hbox`, `\tt`, or `\rm`; GitHub's math
  pipeline does not handle all of these consistently. Prefer `\mathrm{...}` or
  `\text{...}`, and put equation numbers in `\text{(1)}` or in prose.
- Do not use fine-spacing commands such as `\!` or `\,`. Use ordinary spaces,
  `\;`, or `\quad` when extra spacing is actually needed.
- Use `\lVert u\rVert` and `\rVert` for norm bars rather than constructions
  such as `\|\|u\|\|`.
- Use a literal Unicode en dash `–` in prose. Preserve `--` only where it is
  Markdown table syntax or part of a shell command/options.
- Use inline `$...$` math inside tables. Move longer or display-sized equations
  outside tables rather than putting `$$...$$` in a table cell.
- Use proper LaTeX notation instead of shorthand when precision matters. For
  example, write `$L^2(\Omega)$`, `$\partial\Omega$`, and `$\nabla u$` rather
  than plain-text substitutes.
- Put code symbols, class names, functions, commands, options, filenames,
  paths, branch names, and environment variables in backticks: `ProblemSpec`,
  `docs/architecture.md`, or `cmake --build build`.
- Use fenced code blocks with a language tag, such as a `cpp` block or a `bash`
  block.
- Use relative Markdown links for repository files and give links descriptive
  labels.
- Keep headings hierarchical, use blank lines around lists and code blocks,
  and keep tables limited to compact comparisons.
- Before finishing, inspect rendered-looking Markdown for unmatched backticks,
  broken links, malformed tables, and inconsistent notation.

## Handoff

Report the files changed, architectural decisions made, validation commands
run, skipped checks, known limitations, and useful follow-up work.
