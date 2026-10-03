# Repository tools

`tools/` contains repository-local command-line utilities around the numerical
applications and their persisted output.

The tools do **not** implement another PDE, compiler, or optimizer layer. They
either orchestrate an already-built application or consume files that an
application has already written.

The main Chapter 6 flow is:

```text
built nmopt_runner
      │
      ▼
run_chapter6.sh
      │
      ├──► persisted run / native fields
      │          │
      │          ├──► postprocess.py ──► figures / comparisons
      │          │
      │          └──► chapter6_report.py ──► summary tables
      │
      └──► generated output under runs/
```

For installation and build profiles, start with [`BUILD.md`](../BUILD.md).
For the exact run-set, manifest, and artifact contracts, use
[Application execution](../docs/reference/application-execution.md).

## Common workflows

Run these commands from the repository root.

### Run and post-process a Chapter 6 development case

Build the deal.II runner first, then:

```bash
tools/run_chapter6.sh --benchmark b1 --refinement 1
```

or:

```bash
tools/run_chapter6.sh --benchmark b2 --refinement 1
```

The wrapper defaults to:

- `build/debug-dealii/bin/nmopt_runner`;
- a `development` run;
- refinement `1`;
- the repository `runs/` output root;
- `python3`; and
- PNG post-processing.

It runs the selected benchmark, finds the run directory allocated by the
runner, post-processes the persisted output, and generates a manifest-aware
report.

Use:

```bash
tools/run_chapter6.sh --help
```

for runner, Python, output, refinement, and format overrides.

The wrapper **does not build the C++ runner**. See [`BUILD.md`](../BUILD.md)
when the executable is missing or the required profile has not been built.

### Re-render an existing run

Post-processing reads persisted artifacts and native fields; it does not rerun
the PDE:

```bash
python3 tools/postprocess.py \
  --input runs/... \
  --output runs/.../postprocess
```

For one artifact:

```bash
python3 tools/postprocess.py \
  --artifact runs/.../artifacts/... \
  --output runs/.../artifacts/.../postprocess
```

The effective plotting policy is normally reconstructed from the configuration
snapshots stored with the run. See
[Parameter files and plotting profiles](../docs/reference/parameter-files.md)
for the tracked configuration and precedence rules.

### Regenerate a benchmark report

Reports consume the persisted run manifest:

```bash
python3 tools/chapter6_report.py \
  --run-manifest runs/.../run-manifest.json \
  --output runs/.../report
```

The report preserves failed, pending, and missing artifact status rather than
silently reporting only successful coordinates.

### Check Markdown and LaTeX

For files currently changed in the working tree:

```bash
python3 tools/lint_markdown_math.py --changed
```

For all tracked Markdown:

```bash
python3 tools/lint_markdown_math.py
```

The linter is read-only. It checks the repository's portable Markdown/LaTeX
conventions and known GitHub parser hazards; it does not rewrite prose or
mathematics.

### Inspect external deal.II integrations

[`external_dealii/`](external_dealii/) contains source and forward-output checks
used by the external deal.II studies. These tools support the evidence and
comparison workflow; they do not participate in the minimal `nmopt` consumer
path.

Start with the
[Step-4 integration README](../apps/external-dealii/step-4/README.md) before
using those scripts.

## Tool map

| Entry point | Use it when... |
| --- | --- |
| [`run_chapter6.sh`](run_chapter6.sh) | you want one command for a B1/B2 development run, post-processing, and report generation |
| [`postprocess.py`](postprocess.py) | you already have persisted numerical output and want figures or comparisons |
| [`chapter6_report.py`](chapter6_report.py) | you already have a run manifest and want deterministic CSV/Markdown summaries |
| [`chapter6_postprocess.py`](chapter6_postprocess.py) | you are reproducing an older invocation that still uses the compatibility wrapper |
| [`lint_markdown_math.py`](lint_markdown_math.py) | you want to check repository Markdown/LaTeX conventions |
| [`external_dealii/`](external_dealii/) | you are running source/output checks for the external deal.II studies |
| [`nmopt_postprocess/`](nmopt_postprocess/) | you are maintaining the reusable Python implementation behind `postprocess.py` |

## Chapter 6 run wrapper

`run_chapter6.sh` is convenience orchestration around the already-built
repository runner.

Its high-level behavior is:

```text
nmopt_runner
   │
   ├──► run-manifest.json
   ├──► artifact.kv
   ├──► solver traces
   └──► native FE output
             │
             ▼
       postprocess.py
             │
             ▼
       chapter6_report.py
```

The wrapper records the current Git revision by default and lets the runner
allocate the concrete development run slot. Generated data stays under the
ignored `runs/` tree.

For direct `nmopt_runner` commands, authoritative versus development runs,
matrix selection, run slots, manifests, and artifact layout, use
[Application execution](../docs/reference/application-execution.md).

For the meaning and structure of checked-in `.prm` and plotting inputs, use
[Parameter files](../docs/reference/parameter-files.md).

## Post-processing

`postprocess.py` is the public profile-driven rendering/comparison entry point.

The reusable implementation lives under
[`nmopt_postprocess/`](nmopt_postprocess/), including:

- native mesh and field readers;
- solver-history readers;
- field renderers;
- comparison builders;
- plotting-profile handling; and
- derived-output manifests.

The postprocessor consumes persisted numerical output. It must not reconstruct
an authoritative finite-element field by rerunning a solve or inventing data
that the application did not write.

The current Chapter 6 workflow uses NumPy, Matplotlib, and meshio. Those are
runtime tooling dependencies, not C++ build dependencies.

`chapter6_postprocess.py` is retained as a compatibility wrapper for older
commands. New workflows should use `postprocess.py`.

## Reporting

`chapter6_report.py` projects one persisted `run-manifest.json` and its
artifacts into deterministic summary output such as:

```text
report/
  summary.csv
  summary.md
```

It does not rerun the application or postprocessor and does not infer missing
numerical values.

This separation is intentional:

```text
run manifest + artifacts ──► report tables

artifact metadata
+ native fields
+ solver histories ────────► plots / comparisons
```

## Markdown and LaTeX linting

`lint_markdown_math.py` checks the portable subset defined by the
[documentation conventions](../.agents/documentation.md).

Useful invocations are:

```bash
# Only modified or untracked Markdown
python3 tools/lint_markdown_math.py --changed

# Selected files or subtrees
python3 tools/lint_markdown_math.py README.md docs/manual

# Machine-readable findings
python3 tools/lint_markdown_math.py --changed --json

# Available rules
python3 tools/lint_markdown_math.py --list-rules
```

Warnings and errors fail by default. Use `--fail-on error` to fail only on
errors or `--fail-on none` for report-only use.

The linter is deliberately narrower than a full Markdown or TeX parser. It
catches deterministic repository-style and renderer hazards while leaving
semantic choices—such as mathematical notation or prose rewrites—to the
reviewer.

Its focused contract test is:

```bash
python3 tests/tools/markdown_math_lint_contract.py
```

## Dependencies and output boundaries

The tools have deliberately different dependency levels:

| Tool family | Runtime dependencies |
| --- | --- |
| Markdown linter | Python standard library |
| Chapter 6 report | Python standard library |
| Field/history post-processing | Python 3, NumPy, Matplotlib, meshio |
| Chapter 6 run wrapper | built `nmopt_runner`, Bash, Git, plus the post-processing dependencies |

See [`DEPENDENCIES.md`](../DEPENDENCIES.md) for installation and environment
setup.

Generated Python caches and generated analysis output are ignored. Persistent
experiment inputs belong under `parameters/`; generated run evidence and
derived plots belong under `runs/`.
