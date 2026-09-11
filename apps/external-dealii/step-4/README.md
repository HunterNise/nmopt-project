# External deal.II Step-4 boundary fixture

This directory contains one adapted external application and the local code
needed to evaluate its connection to nmopt. The minimum functional path is
kept separate from reference evaluators, verification, and diagnostics:

```text
source/adapted/step-4.cc
        -> integration/problem_a.hpp
        -> integration/nmopt_binding.hpp
        -> existing nmopt public contracts
```

`NmoptBinding` is an experiment-local adapter, not a new application
framework. The evaluation record and roadmap status live in the
[external deal.II boundary roadmap](../../../docs/planning/external-dealii-boundary-evaluation.md).

## Ownership

| Directory or file | Owns | Required by the minimum wiring? |
| --- | --- | --- |
| `source/upstream/step-4.cc` | Verbatim deal.II `v9.5.1` provenance input; never edit. | No |
| `source/baseline/step-4-stripped.cc` | Comment-stripped fidelity baseline. | No |
| `source/adapted/step-4.cc` | Application-owned Step-4 mesh, assembly, solve, and output seams. | Yes |
| `integration/problem_a.hpp` | Problem A residual, objective, derivatives, native solves, and output adapter. Diagnostics are optional. | Yes |
| `integration/nmopt_binding.hpp` | Layouts, callbacks, state/adjoint services, identity metric, and reduced DTO construction using public nmopt contracts. | Yes |
| `evaluation/` | Native reduced and optimization reference paths plus the frozen experiment policy. | No |
| `verification/` | Deterministic scenarios, independent oracle, and comparison checks. | No |
| `diagnostics/` | Evaluation-only counters and solve evidence used to explain the comparison. | No |

The minimum path does not construct an `Instrumentation` object. However,
`problem_a.hpp` and `nmopt_binding.hpp` currently include the diagnostics
header to support their nullable instrumentation pointers. Runtime diagnostic
collection is therefore optional, while the header remains a source
dependency; removing that dependency is outside this closure.

The generic source tools remain under
[`tools/external_dealii/`](../../../tools/external_dealii/). The ignored
`runs/external-dealii/step-4/` tree contains generated evidence only; it is
not an application dependency or a source layout.

## Source fixtures

- [`source/upstream/step-4.cc`](source/upstream/step-4.cc) is the verbatim
  deal.II `v9.5.1` source and must not be edited.
- [`source/baseline/step-4-stripped.cc`](source/baseline/step-4-stripped.cc)
  is generated from the upstream source by
  [`strip_comments.py`](../../../tools/external_dealii/strip_comments.py).
  The tool removes complete comment lines, preserves the leading license
  block and code-bearing lines, and checks non-comment token equivalence.
- [`source/adapted/step-4.cc`](source/adapted/step-4.cc) is copied from the
  stripped baseline and contains only the reusable seams needed by the
  evaluation. It remains independent and does not regenerate the baseline.

The raw source is published at [deal.II Step-4
`v9.5.1`](https://github.com/dealii/dealii/blob/v9.5.1/examples/step-4/step-4.cc).
The pinned revision, retrieval details, hashes, and license attribution are
recorded in the roadmap.

Regenerate-check the committed stripped fixture with:

```bash
python3 tools/external_dealii/strip_comments.py \
  --input apps/external-dealii/step-4/source/upstream/step-4.cc \
  --check apps/external-dealii/step-4/source/baseline/step-4-stripped.cc
```

## Targets and checks

The standalone targets are:

- `nmopt_external_tutorial_step_4` — upstream source;
- `nmopt_external_tutorial_step_4_stripped` — stripped baseline; and
- `nmopt_external_tutorial_step_4_adapted` — adapted source.

They use deal.II directly and do not link nmopt. The native reuse contract is
`nmopt_external_step4_native_contract_test`; it includes the adapted source
and uses only the standard scenario-discovery helper, so it remains
independent of nmopt headers and targets. The application contract and paired
optimization checks consume the `integration/`, `evaluation/`, and
`verification/` files separately.

The three Step-4 C++ test drivers, the forward-comparator contract, and the
`tools/external_dealii/` scripts are verification or evaluation support. They
are not needed to construct or run the minimum binding.

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
  -R '^nmopt\.(external_tutorial_step_4|external\.tutorial_step_4)\.'
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

The evaluators and their reports are verification code, not required pieces
of an external application's nmopt wiring.

## Adapted surface

The adapted source preserves the original 2D/3D forward sequence and exposes
only preparation, const assembled-matrix/RHS views, a supplied-RHS/in-out
vector solve, and a supplied-state VTK writer. Its standalone `main()` is
guarded by `STEP4_NO_MAIN` for native reuse tests. Problem A and the nmopt
binding are application-local integration code; native comparison,
verification, and instrumentation do not belong to the minimum application
path.
