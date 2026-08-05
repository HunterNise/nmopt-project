# Build and test conventions

## Baseline commands

Use an out-of-source `build/` directory with the Ninja generator.

```bash
cmake -S . -B build -G Ninja -DBUILD_TESTING=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

Run the baseline checks after implementation changes. For documentation-only
changes, at least run `git diff --check` and inspect the rendered-looking
Markdown.

## Build boundaries

- Keep generated files and build products under `build/` or another ignored
  build directory.
- When adding a target or test, update the relevant `CMakeLists.txt` and run
  the baseline commands.
- The deal.II backend is optional in environments where deal.II is unavailable.
  Run the backend-neutral tests and report any skipped deal.II checks.
- Do not hide a failed build or test behind a changed command, disabled test,
  or unreported optional dependency.
