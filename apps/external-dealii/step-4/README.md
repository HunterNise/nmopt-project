# External deal.II Step-4 baseline

This directory pins the upstream deal.II Step-4 tutorial as the standalone
baseline for the external-application integration sequence.

## Upstream provenance

- Source: [deal.II Step-4 at `v9.5.1`](https://github.com/dealii/dealii/blob/v9.5.1/examples/step-4/step-4.cc)
- Raw source: [`step-4.cc`](upstream/step-4.cc)
- Upstream revision: `v9.5.1`
- Retrieved: `2026-09-08`
- SHA-256: `be9e694f5f3c9177b7cd18200ff8173337c2b16e1ee72d45ab1ba7c6e105be5f`
- Source license: GNU LGPL 2.1 or later, as stated in the upstream file
- Upstream copyright: deal.II authors; the file credits Wolfgang Bangerth

The installed project dependency is deal.II `9.5.1`, recorded by the local
`deal.IIConfig.cmake`, so this baseline uses the matching upstream tag.

The imported source is verbatim: it contains no nmopt include, namespace, or
link dependency. Local binding code will be added beside `upstream/` in later
roadmap units. The stripped adapted copy in [`step-4.cc`](step-4.cc) is kept
separate from that reference source and is the only source included by the
application binding.

## Standalone target

The project target is `nmopt_external_tutorial_step_4`. It is configured only
with deal.II's target setup and does not link `nmopt_contract`.

CTest registers the same standalone executable as
`nmopt.external_tutorial_step_4.forward`; its working directory is the
ignored `runs/external-dealii/step-4/` directory.

Configure and build it with:

```bash
./build.sh configure debug-dealii
./build.sh build debug-dealii --target nmopt_external_tutorial_step_4
```

Run it into the ignored external-tutorial run directory because the upstream
tutorial writes `solution-2d.vtk` and `solution-3d.vtk` in the current
directory:

```bash
mkdir -p runs/external-dealii/step-4
cmake -E chdir runs/external-dealii/step-4 \\
  ../../../build/debug-dealii/bin/nmopt_external_tutorial_step_4
```

The tutorial remains a forward-only Laplace solve at this stage. Control,
objective, adjoint, and nmopt bindings are deliberately deferred to the next
roadmap units.

## T0 baseline evidence

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

The run produced `runs/external-dealii/step-4/solution-2d.vtk` and
`runs/external-dealii/step-4/solution-3d.vtk` without modifying the imported
source. These generated files are ignored runtime output, not committed
application source or benchmark evidence.

The standalone boundary can be rerun with:

```bash
ctest --test-dir build/debug-dealii --output-on-failure \\
  -R '^nmopt\\.external_tutorial_step_4\\.forward$'
```

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
