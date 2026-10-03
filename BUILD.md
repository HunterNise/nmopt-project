# Build, test, and run

This page is the normal human workflow for configuring, building, testing, and
running the repository.

Use [Dependencies and environment](DEPENDENCIES.md) for installation and
machine setup. Use [Application execution](docs/reference/application-execution.md) for the
exact runner, run-set, manifest, and artifact contracts. Use
[Parameter files](docs/reference/parameter-files.md) for `.prm` and plotting-profile semantics.

## Quick start

From the repository root:

```bash
./build.sh init-config
./build.sh pipeline debug-neutral
```

If deal.II is installed and you want the full backend tests and repository
runner:

```bash
./build.sh pipeline debug-dealii
```

Then inspect the registered applications:

```bash
build/debug-dealii/bin/nmopt_runner --list
```

A small Chapter 6 development run, including post-processing and report
generation, is:

```bash
tools/run_chapter6.sh --benchmark b1 --refinement 1
```

The first `init-config` call is normally needed only once per checkout.

## Machine-local build configuration

`build.sh` is the preferred entry point. It wraps the checked-in CMake presets
and applies the ignored machine-local `build.local.conf` file.

Create that file with:

```bash
./build.sh init-config
```

It records two machine-specific concerns:

- the maximum build job count for each profile; and
- an optional `NMOPT_DEAL_II_DIR` when normal CMake discovery cannot locate
  `deal.IIConfig.cmake`.

Inspect the effective settings with:

```bash
./build.sh show-config
```

Deal.II builds are usually more memory-intensive than neutral builds, so their
job limits can be lowered independently in `build.local.conf`.

## Choose a build profile

| Profile | Use it for | deal.II |
| --- | --- | --- |
| `debug-neutral` | Fast default build and backend-neutral contract tests | No |
| `debug-dealii` | deal.II/compiler/numerical tests and the repository runner | Yes |
| `sanitize-neutral` | AddressSanitizer/UBSan checks of the backend-neutral code | No |
| `release-dealii` | Optimized verification, source-scale runs, and timing evidence | Yes |

For a first checkout, `debug-neutral` is the fastest verification. Use
`debug-dealii` when you want the finite-element realization, compiler-backed
applications, or the full repository runner.

`release-dealii` is not a routine smoke-test profile. It exists for optimized
verification and runs whose scale or timing depends on an optimized deal.II
build.

List the profiles at any time with:

```bash
./build.sh list
```

## What a pipeline does

A pipeline performs the three normal phases for each selected profile:

```text
configure → build → test
```

For example:

```bash
./build.sh pipeline debug-neutral
./build.sh pipeline debug-dealii
```

Multiple profiles may be selected explicitly:

```bash
./build.sh pipeline debug-neutral sanitize-neutral
```

The complete supported set is:

```bash
./build.sh all
```

Running `./build.sh` without a command is equivalent to:

```bash
./build.sh pipeline debug-neutral
```

## Build or test only what you need

The helper also exposes the phases separately:

```bash
./build.sh configure debug-dealii
./build.sh build debug-dealii
./build.sh test debug-dealii
```

Build one target without rebuilding every executable:

```bash
./build.sh build debug-dealii --target nmopt_runner
```

Filter the test phase when working on a focused area:

```bash
./build.sh test debug-dealii --regex REGEX
./build.sh test debug-dealii --label LABEL
```

Use:

```bash
./build.sh --help
```

for pass-through CMake/build/CTest options, `--jobs`, `--clean-first`, dry-run,
and other helper controls.

For repository changes, the stricter verification conventions under
[`.agents/build.md`](.agents/build.md) describe which profiles are expected
for different categories of changes.

## Run the included applications

The repository runner is built by the deal.II profile:

```text
build/debug-dealii/bin/nmopt_runner
build/release-dealii/bin/nmopt_runner
```

List application/catalog entries with:

```bash
build/debug-dealii/bin/nmopt_runner --list
```

Use the runner's own help for its current CLI surface:

```bash
build/debug-dealii/bin/nmopt_runner --help
```

The current runnable Chapter 6 benchmark registrations are `b1` and `b2`.
Their tracked configuration lives under `parameters/chapter-6/`.

### Quick development run

For an ordinary development run, the repository wrapper is the shortest path:

```bash
tools/run_chapter6.sh --benchmark b1 --refinement 1
```

By default the wrapper:

- uses `build/debug-dealii/bin/nmopt_runner`;
- records the run as `development`;
- uses refinement `1`;
- records the current Git revision;
- writes below `runs/`;
- post-processes the persisted numerical output; and
- generates a manifest-aware report.

The wrapper does **not** build the C++ runner. Build `debug-dealii` first, or
build only the runner target:

```bash
./build.sh build debug-dealii --target nmopt_runner
```

Post-processing uses Python 3 plus NumPy, Matplotlib, and meshio. If those
packages live in another environment, supply its interpreter explicitly:

```bash
tools/run_chapter6.sh \
  --benchmark b1 \
  --refinement 1 \
  --python /path/to/python
```

### Direct runner use

Use `nmopt_runner` directly when you need a tracked parameter family, matrix
selection, output override, or exact control over the run kind.

The two normal configuration routes are:

```text
--benchmark ID
    select the registered default parameter family

--parameter-file FILE
    execute one explicit tracked or experimental .prm family
```

Do not combine those two selectors. Every direct run also records a framework
revision.

The exact options, run-kind rules, output placement, manifest lifecycle, and
examples are documented in
[Application execution](docs/reference/application-execution.md).

## Generated output

Repository-generated output is deliberately separated from tracked inputs:

```text
build/<profile>/    configured and compiled build trees
runs/               application runs, manifests, native fields, plots, reports
parameters/         tracked executable inputs
```

`build/` and `runs/` are ignored. A normal application run should not write
generated evidence back into `parameters/` or the source tree.

For Chapter 6 development runs, the runner allocates numbered directories such
as:

```text
runs/chapter-6/b1/development/001/
```

Reproduction runs use an authoritative slot and require the compatible release
profile. Do not treat a Debug/refinement smoke run as source-scale reproduction
evidence.

## Where to look next

| Need | Read |
| --- | --- |
| Install compilers, deal.II, Python packages, or check the environment | [Dependencies and environment](DEPENDENCIES.md) |
| Understand runner options, run sets, manifests, artifacts, and post-processing | [Application execution](docs/reference/application-execution.md) |
| Understand `.prm` files, matrix selection, plotting profiles, and precedence | [Parameter files](docs/reference/parameter-files.md) |
| Use the wrapper, postprocessor, report generator, or documentation linter | [Repository tools](tools/README.md) |
| Understand repository verification rules while modifying code | [Build and test conventions](.agents/build.md) |
| Understand run/evidence conventions while modifying experiments | [Run and post-process conventions](.agents/run.md) |
