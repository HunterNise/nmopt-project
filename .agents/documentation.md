# Documentation and Markdown conventions

These conventions apply to Markdown and LaTeX stored in the repository,
especially the project’s Markdown viewer and GitHub-rendered documentation.

They do not govern Codex chat replies. In chat, use `\(...\)` for inline
mathematics and `\[...\]` for display mathematics, because the chat renderer
uses a different Markdown/LaTeX pipeline.

## Mathematics

- Use `$...$` for inline LaTeX and `$$...$$` for display equations. Do not use
  `\(...\)` or `\[...\]`; the project Markdown viewer does not render those
  delimiters reliably.
- For multiline equations containing underscores, asterisks, escaped braces,
  or environments such as `aligned` and `bmatrix`, prefer a fenced `math`
  block so Markdown does not rewrite the LaTeX:

  ```math
  \begin{aligned}
  E(x) &= 0 \\
  J(x) &= \frac{1}{2}\lVert x \rVert^{2}
  \end{aligned}
  ```

- Always brace superscripts and subscripts, including one-character scripts:
  write `$x_{i}$`, `$x^{\ast}$`, and `$E'(x)^{\ast}$` rather than `$x_i$`,
  `$x^*$`, or `$E'(x)^*$`.
- Use `\ast` instead of a literal `*` for adjoints and duals. A literal
  asterisk can be consumed as Markdown emphasis before math is rendered.
- For inline formulas containing underscores or other Markdown punctuation, use
  GitHub’s backtick-delimited math form:

  ```text
  $`x_{i}`$
  ```

- Do not use `\operatorname`, `\tag`, `\hbox`, `\tt`, or `\rm`; prefer
  `\mathrm{...}` or `\text{...}`. Put equation numbers in `\text{(1)}` or in
  prose.
- Do not use fine-spacing commands such as `\!`, `\,`, or `\;`; some viewers
  display the punctuation literally. Use ordinary spaces, `\quad`, or
  explicit operators such as `\mathrm{d}t` when extra spacing is needed.
- Use `\lVert u\rVert` for norm bars rather than doubled vertical bars.
- Use proper LaTeX notation when precision matters. For example, write
  `$L^{2}(\Omega)$`, `$\partial\Omega$`, and `$\nabla u$` rather than
  plain-text substitutes.
- Use inline `$...$` math inside tables. Move longer or display-sized
  equations outside tables rather than putting `$$...$$` in a table cell.
- In Markdown tables, do not use a literal `|` as a mathematical restriction or
  evaluation bar; it is also the table-cell delimiter. Use `\rvert` instead,
  for example `$u\rvert_{\Gamma}$`.

## Markdown

- Put code symbols, class names, functions, commands, options, filenames,
  paths, branch names, and environment variables in backticks, for example
  `ProblemSpec`, `docs/design/architecture.md`, and `cmake --build build`.
- In current prose, headings, table cells, and ordinary blockquotes, write the
  project/library name as `nmopt`. This includes possessive and compound forms
  such as `nmopt`'s compiler and `nmopt`-native application. Do not add an
  extra code span inside fenced blocks, URLs, literal program output, or larger
  code identifiers such as `nmopt::application`, `nmopt_runner`, or
  `nmopt_contract`. The exact repository title `nmopt-project` may remain plain
  as the root document title.
- Treat established project/tool names such as deal.II, CMake, CTest, Ninja,
  Python, and ParaView as normal prose names. Use backticks only for exact
  symbols, configuration values, commands, or identifiers such as
  `DEAL_II_DIR` or `debug-dealii`.
- Use fenced code blocks with a language tag, such as `cpp`, `bash`, `text`,
  `math`, or `mermaid`.
- Use relative Markdown links for repository files and give links descriptive
  labels.
- Keep headings hierarchical and use blank lines around lists and code blocks.
- Keep tables limited to compact comparisons. Move explanatory prose or long
  equations outside tables.
- Choose diagram format deliberately. Prefer fenced `text` for compact local
  flows, directory trees, record/layout sketches, and visuals where exact
  monospace alignment carries meaning. Prefer fenced `mermaid` when graph
  topology is the main information and automatic layout materially improves a
  multi-node dependency, branching/converging flow, cycle, sequence, or state
  diagram. Do not convert a diagram merely because it is large.
- In new or current `text` diagrams, prefer Unicode structural and arrow
  symbols such as `│`, `─`, `├`, `└`, `┬`, `┴`, `┼`, `→`, `←`, `↑`, `↓`, and
  `↔` over diagrammatic ASCII such as `|`, `+---`, `->`, and `<-`. Preserve
  literal code, CLI syntax, serialized formats, and mathematical notation.
  Historical documents do not need cosmetic diagram conversion unless a
  diagram is broken or misleading.
- When a current document deliberately retains rollout, status, or other
  historical context, use a `>` block with a concise bold label such as
  `**Historical context.**` or `**Historical rollout note.**` when that helps
  prevent the material from being mistaken for current authority.
- Use a literal Unicode en dash `–` in prose. Preserve `--` where it is
  Markdown table syntax or part of a shell command or option.
- Preserve historical evidence as historical evidence. Fix broken links,
  rendering defects, and misleading present-day routing, but do not perform
  bulk cosmetic modernization of `docs/history/` merely to match current
  typography.
- Before finishing, inspect the rendered-looking Markdown for unmatched
  backticks, broken links, malformed tables, and inconsistent notation.
