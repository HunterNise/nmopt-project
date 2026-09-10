# External deal.II Step-4 fixtures

This directory contains the pinned upstream tutorial, a mechanically
comment-stripped copy, and the application-owned adapted copy used by the
external deal.II evaluation. The evaluation record and roadmap status live in
the [external deal.II boundary roadmap](../../../docs/planning/external-dealii-boundary-evaluation.md).

## Fixtures

- [`upstream/step-4.cc`](upstream/step-4.cc) is the verbatim deal.II `v9.5.1`
  source and must not be edited.
- [`baseline/step-4-stripped.cc`](baseline/step-4-stripped.cc) is generated
  from the upstream source by
  [`strip_comments.py`](../../../tools/external_dealii/strip_comments.py).
  The tool removes complete comment lines, preserves the leading license
  block and code-bearing lines, and checks non-comment token equivalence.
- [`step-4.cc`](step-4.cc) is copied from the stripped baseline and contains
  only the reusable seams needed by the evaluation. It remains an independent
  source file; it is not used to regenerate the baseline.

The raw source is published at [deal.II Step-4
`v9.5.1`](https://github.com/dealii/dealii/blob/v9.5.1/examples/step-4/step-4.cc).
The pinned revision, retrieval details, hashes, and license attribution are
recorded in the roadmap.

Regenerate-check the committed stripped fixture with:

```bash
python3 tools/external_dealii/strip_comments.py \
  --input apps/external-dealii/step-4/upstream/step-4.cc \
  --check apps/external-dealii/step-4/baseline/step-4-stripped.cc
```

## Targets and tests

The standalone targets are:

- `nmopt_external_tutorial_step_4` — upstream source;
- `nmopt_external_tutorial_step_4_stripped` — stripped baseline;
- `nmopt_external_tutorial_step_4_adapted` — adapted source.

They use deal.II directly and do not link nmopt. The native reuse contract is
`nmopt_external_step4_native_contract_test`; it includes the adapted source
and uses only the standard scenario-discovery helper, so it also remains
independent of nmopt headers and targets.

Configure and build the standalone fixtures with:

```bash
./build.sh configure debug-dealii
./build.sh build debug-dealii --target nmopt_external_tutorial_step_4
./build.sh build debug-dealii --target nmopt_external_tutorial_step_4_stripped
./build.sh build debug-dealii --target nmopt_external_tutorial_step_4_adapted
```

Run the complete external Step-4 checks with:

```bash
ctest --test-dir build/debug-dealii --output-on-failure \
  -R '^nmopt\\.external_tutorial_step_4'
```

The focused forward comparator can also be run directly. It executes the
supplied programs in separate temporary run directories and compares stdout,
mesh geometry/connectivity, cell types, and numeric VTK arrays:

```bash
python3 tools/external_dealii/check_forward.py \
  --upstream-executable build/debug-dealii/bin/nmopt_external_tutorial_step_4 \
  --stripped-executable build/debug-dealii/bin/nmopt_external_tutorial_step_4_stripped \
  --output-root runs/external-dealii/step-4/forward-comparison \
  --file solution-2d.vtk \
  --file solution-3d.vtk
```

The adapted comparison uses the same command with the adapted executable and
`--stripped-label adapted`. The generated files are ignored run artifacts.
CTest uses stable working directories under
`runs/external-dealii/step-4/forward-comparison/ctest/`; direct comparisons
use a unique run ID under the selected output root.

## Adapted surface

The adapted source preserves the original 2D/3D forward sequence and exposes
only preparation, const assembled-matrix/RHS views, a supplied-RHS/in-out
vector solve, and a supplied-state VTK writer. Its standalone `main()` is
guarded by `STEP4_NO_MAIN` for native reuse tests. Control, objective, adjoint,
and nmopt binding code belong to later evaluation units.
