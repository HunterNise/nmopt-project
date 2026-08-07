# Build and test conventions

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
cmake --build --preset debug-dealii
ctest --preset debug-dealii
```

The persistent profiles are:

| Preset | Generated directory | Purpose |
| --- | --- | --- |
| `debug-neutral` | `build/debug-neutral/` | Fast backend-neutral Debug loop with strict project warnings. |
| `debug-dealii` | `build/debug-dealii/` | Full deal.II Debug correctness gate with project warnings. |
| `sanitize-neutral` | `build/sanitize-neutral/` | Backend-neutral address/undefined-behavior checks. |
| `release-dealii` | `build/release-dealii/` | Optimized verification and Chapter 6 timing profile. |

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

## Focused tests

Every logical scenario is a separate CTest entry with labels and a timeout.
Use labels or names for focused feedback, for example:

```bash
ctest --preset debug-neutral -L semantic
ctest --preset debug-dealii -L dealii
ctest --preset debug-dealii -R coefficient_identification
```

The deal.II executable remains one translation unit and dispatches the named
scenario supplied by CTest; adding a scenario does not require another
expensive executable.

## Dependency intent and cache recovery

`NMOPT_ENABLE_DEAL_II=ON` requires the official deal.II CMake package. Missing
discovery is a configuration error with instructions for `deal.II_DIR` and the
backend-neutral alternative. `NMOPT_ENABLE_DEAL_II=OFF` intentionally performs
no deal.II discovery and registers only backend-neutral tests. Every configure
prints a concise summary of the selected mode.

The preset subdirectories can coexist with a legacy cache directly under
`build/`; do not delete or rewrite that cache automatically. If one specific
preset directory contains a stale different-generator cache, inspect it and
either move that exact generated directory aside or choose a fresh ignored
directory before configuring again.

## Build boundaries

- Keep generated files and build products under `build/` or another ignored
  build directory.
- When adding a target or test, update the relevant `CMakeLists.txt`, assign
  useful labels and a proportionate timeout, and run the applicable profiles.
- The deal.II backend is optional only through the explicit neutral mode. Do
  not hide a requested unavailable dependency, failed build, or failed test.
- Do not add packaging, coverage hosting, formatter gates, or additional
  permanent profiles without a demonstrated project need.
