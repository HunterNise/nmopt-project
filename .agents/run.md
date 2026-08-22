# Application run and post-processing conventions

These conventions govern Chapter 6 parameter placement, runner execution,
generated run organization, post-processing, and experiment promotion. The
exact public schemas remain in the
[application execution reference](../docs/reference/application-execution.md)
and [parameter-file reference](../docs/reference/parameter-files.md).

## Required companion instructions

- Read [Build and test](build.md) before configuring, building, or selecting a
  build profile.
- Read [Workflow](workflow.md) and [Git](git.md) before changing tracked
  parameter files or agent instructions.
- Read [Documentation](documentation.md) before changing the execution or
  parameter references.
- Read the relevant benchmark contract before describing a run as a
  reproduction or using its results as benchmark evidence.

Do not duplicate numerical schemas or benchmark facts here. This file owns the
repository workflow for producing and organizing their generated evidence.

## Parameter-file organization

Tracked parameter files are stable, reviewable inputs:

```text
parameters/chapter-6/<benchmark>/authoritative.prm
parameters/chapter-6/<benchmark>/development/<family>.prm
```

Every tracked Chapter 6 `.prm` uses:

```text
set output root = runs
```

The file's `Benchmark/id` and `Run/kind` select the canonical suffix. The
runner, not the parameter file, allocates the concrete run directory:

```text
runs/chapter-6/<benchmark>/authoritative/
runs/chapter-6/<benchmark>/development/<NNN>/
```

`Run/kind = reproduction` maps to the fixed `authoritative` slot.
`Run/kind = development` maps to the next zero-padded development slot. Do not
put `chapter-6/<benchmark>/<run-slot>` into `Run/output root`; doing so repeats
the canonical suffix.

Ignored, disposable experiment files belong directly below:

```text
runs/parameters/<experiment>.prm
```

They use a distinct experiment root:

```text
set output root = runs/experiments/<experiment>
```

Their concrete runs retain the canonical suffix below that root. An ignored
experiment therefore normally produces:

```text
runs/experiments/<experiment>/chapter-6/<benchmark>/development/<NNN>/
```

If an ignored parameter file includes a tracked file, it must override
`Run/output root` explicitly so the inherited stable root is not used.

## Choosing and building the runner

Reuse an existing compatible runner when possible. Do not rebuild merely to
regenerate reports or plots.

- Use `debug-dealii` for smoke runs and ordinary development checks.
- Use `release-dealii` for source-scale evidence and any reported timing.
- Never build a deal.II profile with more than one job. Use `--parallel 1`.
- Do not configure or rebuild `release-dealii` without the authorization
  required by [Build and test](build.md).

The parameter file records its intended build profile. A reproduction run must
use a compatible release runner; do not relabel Debug output as reproduction
evidence.

## Executing a parameter family

Run a tracked stable family with its checked-in parameter path:

```bash
build/<profile>/bin/nmopt_runner \
  --parameter-file parameters/chapter-6/<benchmark>/development/<family>.prm \
  --framework-revision <revision>
```

For the stable reproduction family, use
`parameters/chapter-6/<benchmark>/authoritative.prm` instead.

Run an ignored experiment from `runs/parameters/` in the same way:

```bash
build/<profile>/bin/nmopt_runner \
  --parameter-file runs/parameters/<experiment>.prm \
  --framework-revision <revision>
```

Use `--select AXIS=VALUE` for a focused subset of the declared matrix. Use
`--output` only for an intentional destination-level override; it replaces the
parameter file's output root, not the canonical benchmark/run-slot suffix.
The parameter file owns `Run/kind`, so do not try to override it on the command
line.

Prefer the smallest matrix that answers the current question. A quiet runner
may spend several minutes in a source-sized iterative solve; confirm the live
process before treating absent output as failure.

## Verifying a numerical run

Do not infer success solely from emitted artifact paths. Inspect the concrete
run's `run-manifest.json` and require:

- `status` is `complete`;
- `success_count` equals `expected_artifact_count`;
- `failure_count` and `pending_count` are zero;
- the recorded build profile, framework revision, parameter hash, and resolved
  matrix match the intended run.

Use `artifact.kv` and native field files as numerical evidence. Report timing
only from a release run. Generated manifests and copied parameter snapshots are
provenance records; do not edit them to repair a moved or stale run.

## Post-processing

Post-process the concrete run-set directory, not the parameter output root:

```bash
python3 tools/postprocess.py \
  --input <concrete-run-directory> \
  --format png
```

By default the tool consumes the parameter and plotting-profile snapshots
stored with the run. Supplying a different profile or format is an explicit
derived-output override and must remain visible in post-processing provenance.

After processing, inspect `postprocess/postprocess-index.json`. Require all
artifacts to succeed and no comparison errors. The `--format` option governs
derived plots; the native runner may independently retain a mesh SVG as part
of its native evidence contract.

Do not modify native VTU fields to improve a plot. Treat rendering failures as
post-processing failures unless the native artifact itself is demonstrably
invalid.

## Promoting an experiment

When an ignored experiment becomes useful beyond the investigation:

1. add a clearly named development `.prm` below `parameters/`;
2. change its output root to `runs`;
3. retain comments and provenance that distinguish hypotheses from source
   facts;
4. add or extend a focused parameter-file contract; and
5. update the parameter reference when it inventories the new stable family.

Do not promote generated run directories, copied run snapshots, or derived
plots as configuration. Evidence below `runs/` remains ignored and disposable
unless a benchmark report explicitly records its findings.

## Run handoff

Report the parameter file, concrete run directory, build profile, framework
revision, artifact counts, stopping reasons, and material failures. Link the
run manifest and requested derived outputs. State clearly whether the result
is a smoke run, a development experiment, or reproduction evidence.
