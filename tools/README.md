# Repository tools

This directory contains repository-local scripts and reusable analysis tooling.
It is not a second application or numerical implementation layer: tools
consume public interfaces and persisted artifacts.

## Entry points

| Tool | Purpose |
| --- | --- |
| `run_chapter6.sh` | Run an already-built Chapter 6 B1/B2 executable, post-process its run set, and generate its report. |
| `postprocess.py` | Profile-driven field rendering and run-root comparison entry point. |
| `chapter6_postprocess.py` | Compatibility wrapper for older Chapter 6 post-processing invocations. |
| `chapter6_report.py` | Deterministic CSV/Markdown benchmark report from a persisted run manifest. |

The reusable Python implementation is under
`nmopt_postprocess/`. Its Chapter 6 profile, mesh/field and solver-history
readers, renderers, comparison builders, and output manifests are
implementation details of the public post-processing entry point.

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

## Dependencies and output boundaries

The post-processing tools use `python3`, `meshio`, and `matplotlib`.
The report generator itself reads persisted records and has no third-party
dependencies. These Python dependencies are runtime tooling dependencies, not
C++ or CTest dependencies.

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
