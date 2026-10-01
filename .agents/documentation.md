# Documentation and Markdown conventions

These conventions apply to Markdown and LaTeX stored in the repository,
especially GitHub-rendered documentation, VS Code's built-in Markdown preview,
and other conventional math-aware Markdown readers. The goal is to use a
small, semantically clear TeX subset that renders correctly across those
readers while avoiding syntax that GitHub's Markdown processing can reinterpret
before the math renderer sees it.

When several equivalent spellings exist, prefer a supported named TeX command
over a literal or escaped Markdown punctuation character. Prefer ordinary
portable `$...$` and `$$...$$` math; use GitHub-specific protected inline math
or fenced `math` blocks only as escape hatches when the natural TeX cannot be
made robust without renderer-specific escaping or source distortion.

These conventions do not govern Codex chat replies. In chat, use `\(...\)`
for inline mathematics and `\[...\]` for display mathematics,
because the chat renderer uses a different Markdown/LaTeX pipeline.

## Mathematics

### Delimiters and source layout

- Use `$...$` for ordinary inline mathematics.
- Use `$$...$$` for ordinary display mathematics. Put display delimiters on
  their own lines and keep blank lines around the display block.
- Plain source newlines inside `$...$` or `$$...$$` are TeX whitespace; they
  do not need to become visible mathematical line breaks. Wrap long source for
  readability without introducing `\\`, `aligned`, `array`, or similar
  environments unless the rendered equation genuinely needs multiple lines or
  alignment.
- Preserve readable TeX source wrapping unless the exact source text forms a
  Markdown construct. Do not treat punctuation at the beginning of a source
  line as unsafe by itself. In particular, `-`, `+`, and `*` are list markers
  only when the required following whitespace is present, and `=` or `-`
  acts as a Setext heading underline only when the line has the required
  underline form. Reflow only actual Markdown parsing hazards, while keeping
  the rendered mathematics unchanged.
- Use `\begin{aligned}...\end{aligned}` when the rendered mathematics is
  intentionally multiline or aligned.
- Prefer `cases` for genuine piecewise definitions. Use `array` only when its
  more general layout is actually required.
- GitHub's protected inline form ``$`...`$`` is an escape hatch for an inline
  expression whose natural TeX contains unavoidable Markdown-sensitive syntax.
  Prefer ordinary `$...$` whenever a clean portable spelling exists. Do not
  use protected math for code or identifiers that should instead be code spans.
- A fenced `math` block is an escape hatch for display mathematics that cannot
  be expressed naturally and safely in `$$...$$`. It protects the TeX from
  Markdown parsing on GitHub and VS Code, but other Markdown readers may show it
  as a code block. Do not use fenced `math` merely because an equation is long
  or multiline.
- Do not use `\(...\)` or `\[...\]` in repository Markdown.

### Canonical scripts and symbols

- Always brace superscripts and subscripts, including one-character scripts:
  write `$x_{i}$`, `$x^{2}$`, `$x^{\ast}$`, and `$E'(x)^{\ast}$` rather than
  `$x_i$`, `$x^2$`, `$x^*$`, or `$E'(x)^*$`. Raw `_` and `^` are valid TeX;
  the braces are the canonical project form.
- Check Markdown-sensitive `_` and `*` contextually rather than treating their
  presence as unsafe by itself. Inspect the surrounding paragraph or table cell,
  because Markdown emphasis delimiters can pair across separate `$...$` spans.
  Ordinary scripts such as `$x_{i}$`, `$\Gamma_{D}$`, and
  `$H_{0}^{1}(\Omega)$` do not by themselves require protected math.
- When such a collision occurs, first prefer a mathematically natural rewrite in
  the portable TeX subset. If preserving the natural notation still leaves an
  unavoidable Markdown collision, use GitHub's protected inline ``$`...`$`` form
  as an escape hatch rather than introducing artificial TeX or changing the
  mathematical notation merely to avoid punctuation.
- Prefer supported named mathematical commands over ambiguous punctuation when
  the meaning is the same. In particular:

  | Meaning | Preferred form | Avoid as the project form |
  | --- | --- | --- |
  | Asterisk / adjoint star | `\ast` or `\star` | literal `*`, `\*` |
  | Multiplication | juxtaposition, `\cdot`, or `\times` as appropriate | literal `*` |
  | Less / greater | `\lt`, `\gt` | literal `<`, `>` |
  | Set braces | `\lbrace`, `\rbrace` | `\{`, `\}` |
  | Scalable braces | `\left\lbrace`, `\right\rbrace` | `\left\{`, `\right\}` |
  | Absolute value | `\lvert x\rvert` | raw `|x|` |
  | Norm | `\lVert x\rVert` | `||x||`, `\|x\|` |
  | Set-builder separator | `\mid` | raw `|` |
  | Restriction / evaluation | `\rvert_{\Gamma}` | raw `|_{\Gamma}` |
  | Pairing | `\langle u,v\rangle` | `<u,v>` |
  | Map/type colon | `F\colon X\to Y` | `F:X\to Y` |
  | Definition | `a\coloneqq b` | `a := b` |

- Do not solve a Markdown collision by inserting a Markdown backslash escape
  into TeX. Forms such as `\*` and `\_` can be consumed differently by
  Markdown and can be invalid commands in KaTeX. Use the proper TeX command or
  a protected-math escape hatch instead.
- Do not replace `_` or `^` with `\sb` or `\sp`; use ordinary braced TeX
  scripts.

### Operators, text, and fonts

- Use the standard built-in operator command when one exists, for example
  `\min`, `\max`, `\inf`, `\sup`, `\lim`, `\ker`, `\det`, `\log`, `\exp`,
  `\sin`, and `\cos`. Small spacing differences between renderers are
  acceptable; do not replace a semantic built-in operator merely to force
  pixel-identical output.
- For a custom named operator, use `\mathop{\mathrm{...}}`, for example
  `\mathop{\mathrm{span}}`, `\mathop{\mathrm{diag}}`,
  `\mathop{\mathrm{div}}`, or `\mathop{\mathrm{tr}}`.
- When a custom operator carries side scripts, use `\nolimits` where that is the
  intended mathematical layout and it gives a Markdown-safe spelling, for
  example `$\mathop{\mathrm{tr}}\nolimits_{\Gamma_{D}} v$`.
  Do not add `\nolimits` mechanically to operators without scripts or where
  limits-style placement is mathematically intended.
- Do not use `\operatorname` or `\operatorname*`; GitHub currently rejects
  those macros even though stock MathJax and KaTeX support them.
- Use `\text{...}` for ordinary prose inside mathematics, `\mathrm{...}` for
  roman mathematical labels and symbols, `\mathbf{...}` for bold roman
  mathematical symbols, and `\mathtt{...}` when monospaced mathematical text
  is genuinely intended.
- Write the differential with a roman `d`, for example
  `\thinspace\mathrm{d}x` in explicit integrals.
- Do not use legacy declarations such as `\rm`, `\tt`, or `\bf` in new or
  revised mathematics. They may render, but the scoped modern commands above
  are the canonical project forms.
- Do not use `\mbox`. Do not use `\hbox` for ordinary equation text; use
  `\text{...}` instead. If a true box-level construct is ever required, verify
  it deliberately in every supported renderer before introducing it.

### Spacing

- Fine mathematical spacing is allowed and should be used where it improves
  standard mathematical typography. Do not omit meaningful spacing merely to
  avoid GitHub parsing issues.
- Prefer letter-named or otherwise Markdown-safe forms:

  | Intended spacing | Preferred form |
  | --- | --- |
  | Thin space, approximately `3mu` | `\thinspace` |
  | Negative thin space, approximately `-3mu` | `\negthinspace` |
  | Medium math space, approximately `4mu` | `\mkern4mu` |
  | Thick math space, approximately `5mu` | `\mkern5mu` |
  | Custom math spacing | `\mkern<n>mu` |
  | Large structural space | `\quad` or `\qquad` |

- Do not use the punctuation control-symbol forms `\,`, `\:`, `\>`, `\;`, or
  `\!` in dollar-delimited repository math; GitHub Markdown can consume the
  backslash before the math renderer sees it.
- Do not use `\medspace` or `\thickspace`; they are not accepted consistently
  by the supported GitHub rendering path even though other engines support
  them.
- Do not add spacing mechanically. Evaluate the mathematical role first. In
  particular, use `\thinspace` before an explicit differential in an integral;
  use `\mkern` only when a specific fine spacing is mathematically or
  typographically justified; use `\quad` and `\qquad` only for genuinely large
  structural separation.

### Equation layout and numbering

- For a single rendered equation, source wrapping may span several Markdown
  lines, but do not insert `\\` unless a visible line break is intended.
- For intentionally multiline equations, prefer `aligned` inside the display
  delimiters. Avoid top-level `align` in repository Markdown; its numbering and
  layout behavior varies more across renderers.
- Do not use `\tag`. It is invalid in useful inner environments such as
  `aligned` and its placement can overlap or differ across supported
  renderers. For a short manually numbered equation, place the number as
  ordinary math text such as `\qquad\text{(1)}` when that renders cleanly;
  otherwise put the equation number or reference in prose.
- Do not use `\label`, `\ref`, or `\eqref` for Markdown equation navigation.
  Use prose and Markdown links instead.

### Tables and existing mathematics

- Use inline `$...$` math inside Markdown tables. Move display-sized,
  multiline, fenced, or otherwise complicated mathematics outside the table.
- In Markdown tables, never use a literal `|` as a mathematical restriction,
  evaluation bar, norm, or set separator. Use the named forms above.
- Use proper mathematical notation when precision matters, for example
  `$L^{2}(\Omega)$`, `$\partial\Omega$`, and `$\nabla u$` rather than
  plain-text substitutes.
- During an ordinary local edit, do not rewrite unrelated existing mathematics
  solely for house-style normalization. A repository-wide documentation polish
  task may deliberately normalize current documentation when that normalization
  is part of the requested scope.

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
- When a fenced example contains another fenced block, make the outer fence
  longer than the inner fence rather than escaping the inner backticks.
- Use relative Markdown links for repository files and give links descriptive
  labels. Keep headings that are likely link targets textually simple; avoid
  unnecessary inline math, inline code, or punctuation in such headings.
- Keep headings hierarchical and use blank lines around lists and code blocks.
- For an intentional hard line break, use a terminal backslash rather than two
  trailing spaces. Prefer a normal paragraph break when a paragraph boundary
  is intended.
- Keep tables limited to compact comparisons and inline content. Move
  explanatory prose, nested structures, or long equations outside tables.
- Do not use raw HTML merely to work around Markdown rendering when a portable
  Markdown form exists. In particular, do not rely on custom CSS, `style`,
  `class`, or manually assigned `id` attributes for repository documentation.
- Choose diagram format deliberately. Prefer fenced `text` for compact local
  flows, directory trees, record/layout sketches, and visuals where exact
  monospace alignment carries meaning. Prefer fenced `mermaid` when graph
  topology is the main information and automatic layout materially improves a
  multi-node dependency, branching/converging flow, cycle, sequence, or state
  diagram. Do not convert a diagram merely because it is large.
- In new or current `text` diagrams, use Unicode structure deliberately.
  For branching or converging graphs, use box-drawing geometry such as `│`,
  `─`, `├`, `└`, `┬`, `┴`, and `┼`, with `►` for horizontal branch edges and
  `▼` for vertical branch/progression edges. For one-line relations or
  mappings, use `→`, `←`, or `↔`; use mathematical `↦` for mapsto. For a
  vertical linear flow or pipeline, use `↓`. Avoid diagrammatic ASCII such as
  `|`, `+---`, `->`, `-->`, and `<-`. Preserve literal code, CLI syntax,
  serialized formats, C++ member access, Mermaid edge syntax, and mathematical
  notation that is not being used as diagram geometry. Historical documents do
  not need cosmetic diagram conversion unless a diagram is broken or
  misleading.
- When a current document deliberately retains rollout, status, or other
  historical context, use a `>` block with a concise bold label such as
  `**Historical context.**` or `**Historical rollout note.**` when that helps
  prevent the material from being mistaken for current authority.
- Use a literal Unicode en dash `–` in prose. Preserve `--` where it is
  Markdown table syntax or part of a shell command or option.
- Preserve historical evidence as historical evidence. Scan `docs/history/`
  for broken rendering, broken links, and misleading present-day routing, but
  do not bulk-modernize historical typography merely to match current style.
- Before finishing a Markdown edit, inspect the changed regions for unmatched
  or incorrectly nested backticks and fences, malformed math delimiters,
  Markdown-active continuation lines inside `$$...$$`, literal Markdown
  punctuation used where the preferred named TeX form exists, unsupported or
  discouraged macros, malformed tables, broken relative links or heading
  fragments, accidental raw HTML, and inconsistent notation.
