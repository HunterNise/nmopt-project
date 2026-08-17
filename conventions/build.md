# Build and test conventions

## Read-only tasks

For read-only tasks, do not configure, build, or run tests unless the user
explicitly requests verification or the answer depends on test output.
Existing build artifacts may be inspected, but should not be regenerated. If
nothing was edited, skip build and test by default.

## Baseline profiles

Use the checked-in CMake presets. They select Ninja, an explicit build type,
testing, dependency intent, warning policy, and a distinct ignored directory.
The fast baseline before and after every task is:

```bash
cmake --preset debug-neutral
cmake --build --preset debug-neutral
ctest --preset debug-neutral
```

Changes that touch deal.II code, compiler/lowering, numerical behavior, or
backend test registration must also run:

```bash
cmake --preset debug-dealii
cmake --build --preset debug-dealii --parallel 1
ctest --preset debug-dealii
```

The persistent profiles are:

| Preset | Generated directory | Purpose |
| --- | --- | --- |
| `debug-neutral` | `build/debug-neutral/` | Fast backend-neutral Debug loop with strict project warnings. |
| `debug-dealii` | `build/debug-dealii/` | Full deal.II Debug correctness gate with project warnings. |
| `sanitize-neutral` | `build/sanitize-neutral/` | Backend-neutral address/undefined-behavior checks. |
| `release-dealii` | `build/release-dealii/` | Optimized verification and Chapter 6 timing profile. |

## Tooling and profile cost

Use `python3` explicitly for repository Python scripts. The unversioned
`python` command is not guaranteed to exist across Linux distributions and
may refer to a different Python major version when it is provided. Check the
selected interpreter before running a script when the environment is
uncertain:

```bash
which python3
python3 --version
```

Report generation and refinement-1 development runs do not require the
optimized profile. Reuse or build the Debug deal.II runner for those tasks:

```bash
cmake --build --preset debug-dealii --target nmopt_runner
```

The `release-dealii` build is optimized and can take substantially longer on
a local machine. Build it only when refreshing source-scale reproduction
artifacts or the optimized verification gate; do not rebuild it merely to
regenerate reports or inspect existing debug artifacts.

Deal.II-enabled builds must always use one build job (`--parallel 1`). The
deal.II translation units are memory-intensive, and allowing Ninja to compile
multiple large units concurrently can exhaust a constrained development
environment. This restriction applies to `debug-dealii` and
`release-dealii`; backend-neutral profiles may use their normal parallelism.

Warnings are errors in the Debug and sanitizer profiles. Strict conversion and
shadow warnings remain backend-neutral until the audited deal.II adapter seams
are clean; do not silence external-header warnings with global compiler flags.
The sanitizer test preset disables leak detection because LeakSanitizer cannot
run in the ptrace-based agent environment. To request leak detection in a
compatible environment, bypass that test preset after building:

```bash
ASAN_OPTIONS=detect_leaks=1 \
  ctest --test-dir build/sanitize-neutral --output-on-failure
```

Run `release-dealii` before reporting numerical timings. Debug timings are not
benchmark evidence.

Each profile keeps CMake's own state at the profile root, while project outputs
use this layout:

```text
build/debug-neutral/
  bin/                         # executable targets
  lib/                         # library targets
  generated/scenarios/         # generated CTest registrations
    nmopt_reduced_dto_contract_test/
      discovered_scenarios.cmake
      include.cmake
```

The same layout applies to the other profiles. The
`discovered_scenarios.cmake` files are produced after their registered test
executables are built, while the small `include.cmake` wrappers are generated
at configure time. Neither is a source file; do not edit them by hand.

### Unix Makefiles fallback

The checked-in presets intentionally require Ninja. If Ninja is unavailable,
reuse a configure preset's cache variables while overriding its generator and
generated directory locally:

```bash
cmake --preset debug-neutral \
  -G "Unix Makefiles" \
  -B build/debug-neutral-make
cmake --build build/debug-neutral-make
ctest --test-dir build/debug-neutral-make --output-on-failure
```

Substitute another configure preset and a matching `-make` directory when that
profile is required. Its checked-in build and test presets still refer to the
Ninja directory, so use the generator-independent `cmake --build` and
`ctest --test-dir` forms for this fallback. These local directories are not
additional persistent project profiles.

## Focused tests

Every logical scenario is a separate CTest entry with labels and a timeout.
Use labels or names for focused feedback, for example:

```bash
ctest --preset debug-neutral -L semantic
ctest --preset debug-dealii -L dealii
ctest --preset debug-dealii -R coefficient_identification
```

Each test source owns its scenario names, CTest names, labels, timeouts, and
callbacks in one registry. The current convention is one test executable per
test translation unit: `tests/foo.cc` becomes the project target
`nmopt_foo_test`. A translation unit may still contain multiple logical
scenarios.

After linking a test executable, the generic `scenario_discovery_register`
helper queries its `--list-scenarios` manifest and generates the individual
CTest entries. Adding or removing a logical scenario therefore changes only
the test source; do not duplicate its inventory in `CMakeLists.txt`. Separate
processes localize ordinary failures, crashes, and timeouts while allowing
later CTest entries to continue. Running an executable without a scenario
remains an aggregate convenience, and its failure text identifies the active
scenario.

## Dependency intent and cache discipline

`NMOPT_ENABLE_DEAL_II=ON` requires the official deal.II CMake package. Missing
discovery is a configuration error with instructions for `deal.II_DIR` and the
backend-neutral alternative. `NMOPT_ENABLE_DEAL_II=OFF` intentionally performs
no deal.II discovery and registers only backend-neutral tests. Every configure
prints a concise summary of the selected mode.

The root of `build/` is a container only. Never configure or build directly in
`build/`; always use a named subdirectory such as `build/debug-neutral/` or
`build/debug-dealii/`. A CMake cache, `CMakeFiles/`, or build product directly
under `build/` is invalid legacy output and should not be reused. If a named
profile directory contains a stale or different-generator cache, remove only
that generated directory and configure the profile again, or select a fresh
ignored directory when its old contents must be retained temporarily.

## Build boundaries

- Keep generated files and build products under `build/` or another ignored
  build directory.
- Add ordinary backend-neutral executables with the neutral executable helper;
  add tests with its test variant so scenario discovery is registered exactly
  once. Use the corresponding deal.II helpers only inside the deal.II-enabled
  block. When adding a scenario to an existing executable, update only its C++
  registry with useful labels and a proportionate timeout, then rebuild before
  running CTest.
- The deal.II backend is optional only through the explicit neutral mode. Do
  not hide a requested unavailable dependency, failed build, or failed test.
- Do not add packaging, coverage hosting, formatter gates, or additional
  permanent profiles without a demonstrated project need.
