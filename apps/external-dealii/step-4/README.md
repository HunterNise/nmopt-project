# External deal.II Step-4 baseline

This directory pins the upstream deal.II Step-4 tutorial and its
comment-line-stripped copy as standalone forward references for the
external-application integration sequence.

## Upstream provenance

- Source: [deal.II Step-4 at `v9.5.1`](https://github.com/dealii/dealii/blob/v9.5.1/examples/step-4/step-4.cc)
- Raw source: [`step-4.cc`](upstream/step-4.cc)
- Upstream revision: `v9.5.1`
- Retrieved: `2026-09-08`
- SHA-256: `be9e694f5f3c9177b7cd18200ff8173337c2b16e1ee72d45ab1ba7c6e105be5f`
- Source license: GNU LGPL 2.1 or later, as stated in the upstream file
- Upstream copyright: deal.II authors; the file credits Wolfgang Bangerth

The stripped baseline is [`step-4-stripped.cc`](baseline/step-4-stripped.cc).
It is generated from the raw source by
[`strip_comments.py`](../../../tools/external_dealii/strip_comments.py). The
tool removes complete comment lines while preserving the leading license
block, code-bearing lines, inline comments, literals, and whitespace. It also
checks non-comment token equivalence. Its current SHA-256 is
`b21212764c50401089612c6ac2bb196395e3ac9c261120e0513be341825399b2`.

The installed project dependency is deal.II `9.5.1`, recorded by the local
`deal.IIConfig.cmake`, so this baseline uses the matching upstream tag.

The imported source is verbatim: it contains no nmopt include, namespace, or
link dependency. The stripped baseline is a separate generated reference;
the adapted [`step-4.cc`](step-4.cc) remains the later application-binding
source and is not used to generate the baseline.

## Standalone target

The raw project target is `nmopt_external_tutorial_step_4`; the stripped
target is `nmopt_external_tutorial_step_4_stripped`. Both use only deal.II's
target setup and do not link `nmopt_contract`.

CTest registers the same standalone executable as
`nmopt.external_tutorial_step_4.forward`, and registers the stripped target as
`nmopt.external_tutorial_step_4.stripped_forward`. Their working directories
are separate subdirectories of the ignored
`runs/external-dealii/step-4/forward-comparison/` directory.

Configure and build it with:

```bash
./build.sh configure debug-dealii
./build.sh build debug-dealii --target nmopt_external_tutorial_step_4
./build.sh build debug-dealii --target nmopt_external_tutorial_step_4_stripped
```

Run and compare both executables with the reusable forward-checking tool:

```bash
python3 tools/external_dealii/check_forward.py \\
  --upstream-executable build/debug-dealii/bin/nmopt_external_tutorial_step_4 \\
  --stripped-executable build/debug-dealii/bin/nmopt_external_tutorial_step_4_stripped \\
  --output-root runs/external-dealii/step-4/forward-comparison \\
  --file solution-2d.vtk \\
  --file solution-3d.vtk
```

Each invocation creates a unique run directory containing separate upstream
and stripped outputs, stdout/stderr logs, and `comparison.txt`. The comparison
checks exact stdout, mesh points, cell connectivity, cell types, and numeric
VTK arrays. It uses no roadmap-unit name in the artifact path.

The tutorial remains a forward-only Laplace solve at this stage. Control,
objective, adjoint, and nmopt bindings are deliberately deferred to the next
roadmap units.

## Historical T0 baseline evidence

The standalone target was configured and built with the repository's Debug
deal.II profile, then run from a dedicated output directory. The observed
output was:

```text
Solving problem in 2 space dimensions.
   Number of active cells: 256
   Total number of cells: 341
   Number of degrees of freedom: 289
   26 CG iterations needed to obtain convergence.
Solving problem in 3 space dimensions.
   Number of active cells: 4096
   Total number of cells: 4681
   Number of degrees of freedom: 4913
   30 CG iterations needed to obtain convergence.
```

The historical run produced root-level files under
`runs/external-dealii/step-4/` without modifying the imported source. Those
ignored files are retained as legacy local artifacts, not current comparison
evidence.

The standalone boundary can be rerun with:

```bash
ctest --test-dir build/debug-dealii --output-on-failure \\
  -R '^nmopt\\.external_tutorial_step_4\\.(forward|stripped_forward|forward_comparison)$'
```

The E1 comparison is accepted only after the current run records matching
2D/3D stdout and numerical VTK content. Historical counts above are not
substituted for that current evidence.

## E1 current forward-comparison evidence

The first current raw-versus-stripped comparison was run on `2026-09-10`
with the Debug deal.II profile, using `HEAD` `0f82566` plus the uncommitted
E1.b source and CMake changes. The required command was:

```bash
./build.sh pipeline debug-dealii
```

The complete pipeline passed all 160 tests. The focused comparison created
the ignored artifact directory
`runs/external-dealii/step-4/forward-comparison/20260910T075148Z-e0a3b3aa/`.
Its `comparison.txt` records identical stdout, identical mesh points, cell
connectivity, cell types, and one numeric VTK array for both output files.
The maximum point and array differences were both zero:

```text
solution-2d.vtk: points=1024, cells=256, arrays=1, max_point_difference=0, max_array_difference=0
solution-3d.vtk: points=32768, cells=4096, arrays=1, max_point_difference=0, max_array_difference=0
```

The observed standalone counts were 289/4913 degrees of freedom and 26/30
CG iterations for 2D/3D respectively on both executables. The raw and
stripped VTK file hashes are not expected to match because their generated
date headers differ; the comparator ignores those headers and compares the
parsed numerical contents.

## T2 application binding

The adapted [`step-4.cc`](step-4.cc) removes tutorial exposition, selects the
2D path, and adds only a preparation entry point plus read-only views of the
tutorial-owned deal.II state. It contains no control or objective code.

[`tutorial_application.cc`](tutorial_application.cc) adds a compact
application-owned binding around that seam. The binding uses a full-coordinate
affine control, the tutorial operator as the residual operator, and one
finite-element mass assembly reused for the state objective and control
metric. It does not build a separate control-coupling matrix or duplicate the
tutorial PDE assembly.

Build and run the binding smoke target with:

```bash
./build.sh build debug-dealii \
  --target nmopt_external_tutorial_binding_smoke
ctest --test-dir build/debug-dealii --output-on-failure \
  -R '^nmopt\\.external_tutorial_step_4\\.binding$'
```

The smoke path verifies the zero-control state residual, objective derivative
layout, metric dimensions, and native VTU output. Both state and control have
dimension `289` for this fixed 2D realization.
