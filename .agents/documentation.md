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

## Markdown

- Put code symbols, class names, functions, commands, options, filenames,
  paths, branch names, and environment variables in backticks, for example
  `ProblemSpec`, `docs/design/architecture.md`, and `cmake --build build`.
- Use fenced code blocks with a language tag, such as `cpp`, `bash`, or `text`.
- Use relative Markdown links for repository files and give links descriptive
  labels.
- Keep headings hierarchical and use blank lines around lists and code blocks.
- Keep tables limited to compact comparisons. Move explanatory prose or long
  equations outside tables.
- Use a literal Unicode en dash `–` in prose. Preserve `--` where it is
  Markdown table syntax or part of a shell command or option.
- Before finishing, inspect the rendered-looking Markdown for unmatched
  backticks, broken links, malformed tables, and inconsistent notation.
