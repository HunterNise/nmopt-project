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
- On `main` or any branch that does not match `codex/*`, do not edit, create,
  delete, or rename repository files unless the user explicitly requests that
  change. Read-only inspection and validation are allowed. If a task might
  require edits but the authorization is unclear, first describe the intended
  changes and ask for permission.
- On a `codex/*` branch, Codex may make in-scope project changes without
  requesting per-file permission. It must still follow the architecture,
  validation, and preservation rules in this file and must not switch branches
  on its own.
- After completing each requested change, or each coherent atomic part of a
  larger change, inspect the diff and suggest a Git commit message in the final
  reply on `codex/*` branches. Do not run `git commit`, ask for commit approval,
  or include unrelated user changes; leave the changes pending so the user can
  review or adjust them before committing.
- Do not merge, push, force-push, or rewrite history unless explicitly asked.
- Prefer small, coherent changes with a clear handoff.
- Do not commit build products, generated files, editor files, or secrets.
- Do not change the Git author or committer identity, or rewrite existing
  authorship, unless the user explicitly requests it.

### Commit messages

Use this Conventional Commits format:

```text
<type>(<scope>): <imperative summary>
```

- Use a lowercase `type` such as `feat`, `fix`, `refactor`, `docs`, `test`,
  `build`, `ci`, `chore`, `perf`, or `revert`.
- Use a short, lowercase `scope` when it clarifies the affected component;
  omit it when it adds no information.
- Keep the summary imperative, specific, and preferably within 72 characters.
  Do not end it with a period.
- Add a body after a blank line when the reason, tradeoff, or compatibility
  impact is not clear from the summary.
- Mark breaking changes with `!` after the type or scope and explain them in
  the body or a `BREAKING CHANGE:` footer.
- Keep each commit focused on one coherent purpose.

Examples:

```text
docs(agents): define Codex edit and commit policy
docs(markdown): fix GitHub math rendering
refactor(spec): separate semantic contracts from backends
```

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

- Use `$...$` for inline LaTeX and `$$...$$` or fenced `math` blocks for
  display equations. Prefer a fenced `math` block when a multiline equation
  contains underscores, asterisks, escaped braces, or environments such as
  `aligned` and `bmatrix`; this prevents Markdown from rewriting the LaTeX.
  Do not use `\(...\)` or `\[...\]`, because the project Markdown viewer does
  not render those delimiters reliably.
- Always brace superscripts and subscripts, including one-character scripts:
  write `x_{i}`, `x^{\ast}`, and `E'(x)^{\ast}` rather than `x_i`, `x^*`, or
  `E'(x)^*`. Use `\ast` instead of a literal `*` for adjoints and duals, since
  the literal asterisk can be consumed as Markdown emphasis before math is
  rendered. This avoids parser ambiguity around primes and nested scripts.
- For inline formulas containing underscores or other Markdown punctuation, use
  GitHub's backtick-delimited math form:

  ```text
  $`x_{i}`$
  ```

- Do not use `\operatorname`, `\tag`, `\hbox`, `\tt`, or `\rm`; GitHub's math
  pipeline does not handle all of these consistently. Prefer `\mathrm{...}` or
  `\text{...}`, and put equation numbers in `\text{(1)}` or in prose.
- Do not use fine-spacing commands such as `\!`, `\,`, or `\;`; some project
  viewers display the punctuation literally. Use ordinary spaces, `\quad`, or
  explicit operators such as `\mathrm{d}t` when extra spacing is needed.
- Use `\lVert u\rVert` for norm bars rather than doubled vertical bars.
- Use a literal Unicode en dash `–` in prose. Preserve `--` only where it is
  Markdown table syntax or part of a shell command/options.
- Use inline `$...$` math inside tables. Move longer or display-sized equations
  outside tables rather than putting `$$...$$` in a table cell.
- Use proper LaTeX notation instead of shorthand when precision matters. For
  example, write `$L^{2}(\Omega)$`, `$\partial\Omega$`, and `$\nabla u$` rather
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
