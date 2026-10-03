# Applications

`apps/` contains concrete consumers of the reusable `nmopt` core. Application code
may assemble repository-specific scenarios, bind runtime data, or preserve an
independently owned PDE implementation, but reusable numerical contracts and
algorithms remain under [`include/nmopt/`](../include/nmopt/).

## Repository runner

[`nmopt-runner/`](nmopt-runner/) is the repository's headless application for the
registered Chapter 6 benchmark families. It combines:

```text
tracked parameters
      ↓
scenario + run-set planning
      ↓
compiler/backend execution
      ↓
artifacts + manifest + native fields
```

For the normal build and run workflow, see
[Build, test, and run](../BUILD.md).

For exact runner behavior, run directories, manifests, and artifacts, see
[Application execution](../docs/reference/application-execution.md). Maintainers
working on the runner implementation should use the
[runner implementation map](../docs/internals/runner.md).

## External deal.II integration

[`external-dealii/`](external-dealii/) demonstrates the peer application-owned path:
an existing deal.II code keeps its mesh, assembly, solves, boundary treatment, and
native output while exposing the numerical operations required by `nmopt`
formulations and solvers.

The current case study is
[Step-4](external-dealii/step-4/README.md).

## Inputs and generated output

Tracked executable configuration lives under [`parameters/`](../parameters/).
Generated application evidence belongs under the ignored `runs/` tree.

The application directories are final consumers and examples of the public
boundaries; they should not become a second home for reusable compiler, formulation,
or optimization implementations.
