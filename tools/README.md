# Repository tools

This directory contains repository-local scripts and reusable analysis tooling.
It is not a second application or numerical implementation layer: tools
consume public interfaces and persisted artifacts.

## Entry points

| Tool | Purpose |
| --- | --- |
| `lint_markdown_math.py` | Read-only linter for the repository's portable Markdown/LaTeX conventions and known GitHub parser hazards. |
| `run_chapter6.sh` | Run an already-built Chapter 6 B1/B2 executable, post-process its run set, and generate its report. |
| `postprocess.py` | Profile-driven field rendering and run-root comparison entry point. |
| `chapter6_postprocess.py` | Compatibility wrapper for older Chapter 6 post-processing invocations. |
| `chapter6_report.py` | Deterministic CSV/Markdown benchmark report from a persisted run manifest. |
| `external_dealii/` | Reusable source and forward-output checks for external deal.II tutorials. |

The reusable Python implementation is under
`nmopt_postprocess/`. Its Chapter 6 profile, mesh/field and solver-history
readers, renderers, comparison builders, and output manifests are
implementation details of the public post-processing entry point.

The `external_dealii/` directory contains small command-line tools that
operate on tutorial source files and standalone executable outputs. They do
not own PDE assembly, numerical methods, or experiment configuration.

For commands, output paths, native-file names, supported formats, and agent
verification, read the [application execution and artifact
reference](../docs/reference/application-execution.md). The wrapper does not
build the C++ runner. Build and test it from the repository root with the
preferred root-level helper:

```bash
./build.sh init-config
./build.sh pipeline debug-dealii
```

When only the runner target is needed, use the atomic build action:

```bash
./build.sh build debug-dealii --target nmopt_runner
```

`build.sh` delegates to the checked-in CMake presets and applies the
machine-local configuration from `build.local.conf`, including the maximum
number of jobs for each profile. Direct CMake and CTest commands remain useful
when a workflow needs lower-level control; see the [dependencies and
environment reference](../DEPENDENCIES.md) for equivalent manual commands
and the [agent build instructions](../.agents/build.md) for agent-specific
verification requirements.

## Markdown and LaTeX linting

`lint_markdown_math.py` is a zero-dependency, read-only check for the portable
Markdown/LaTeX subset defined in
[the documentation conventions](../.agents/documentation.md). It catches known
Markdown-versus-TeX collisions and deterministic project-style violations; it
does not modify files or mechanically choose mathematical rewrites.

Run it from the repository root:

```bash
# All tracked Markdown
python3 tools/lint_markdown_math.py

# Only currently modified or untracked Markdown
python3 tools/lint_markdown_math.py --changed

# A selected file or subtree
python3 tools/lint_markdown_math.py docs/manual README.md

# Machine-readable findings
python3 tools/lint_markdown_math.py --changed --json

# Rule catalogue
python3 tools/lint_markdown_math.py --list-rules
```

Warnings and errors fail by default. Use `--fail-on error` to fail only on
errors, `--fail-on info` to include advisory findings, or `--fail-on none` for
a report-only run. Style-only findings under `docs/history/` are suppressed by
default; use `--include-history-style` only when historical material is
intentionally being modernized.

The linter is deliberately narrower than a full GFM parser or TeX renderer.
Known limits are intentional:

- it does not validate general TeX syntax, balance arbitrary braces/environments,
  or prove that an unlisted macro is supported by every renderer;
- ordinary/protected inline math is expected to stay on one physical source line;
  unusual multiline inline constructs require manual review;
- underscore-emphasis analysis models the project hazards conservatively rather
  than implementing every GFM delimiter-run and nested-link edge case; and
- semantic choices such as whether a bar means `\mid` or `\rvert`, or where
  fine spacing is typographically appropriate, remain reviewer decisions.

When renderer compatibility is uncertain, manually inspect GitHub, VS Code, and
another target Markdown reader rather than extending the linter with a blind
mechanical rewrite.

Run its focused contract tests with:

```bash
python3 tests/tools/markdown_math_lint_contract.py
```

## Dependencies and output boundaries

The Markdown/LaTeX linter and report generator use only Python's standard
library. The post-processing tools use `python3`, `meshio`, and `matplotlib`.
These Python dependencies are runtime tooling dependencies, not C++ or CTest
dependencies.

Tools must:

- read authoritative native fields and artifact records without modifying
  them;
- write derived plots, indexes, and reports under the selected run output;
- preserve missing or failed artifact records instead of inferring values; and
- keep numerical assembly, PDE lowering, and optimization algorithms in the
  application and library layers.

Generated Python caches such as `__pycache__/` are ignored. Do not commit
generated tool output; experiment evidence belongs under the ignored
`runs/` tree.
