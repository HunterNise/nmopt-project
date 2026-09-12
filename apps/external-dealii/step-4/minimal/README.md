# Minimal Problem A consumer

This directory is a small, complete external application that uses the
current nmopt public contracts to optimize the existing Step-4 Problem A.
It is a usability example, not a replacement for the evaluated integration
or for the Step-4 application itself.

## Functional consumer path

The code needed to wire and run the application is deliberately explicit:

```text
minimal/problem_a.cc
  ├── constructs integration/ProblemA
  ├── constructs minimal/ProblemABinding
  ├── creates zero control and the frozen solver policy
  ├── runs ReducedSearchSolverT
  └── sends the retained final state to ProblemA's VTK writer

minimal/problem_a_binding.hpp
  ├── declares variable/test layouts
  ├── supplies the five CallbackExecutableModel callbacks
  ├── declares the state/control partition and native solve services
  ├── translates native CG evidence into LinearSolveReport
  ├── adapts the identity control metric
  └── constructs the ReducedDTO

integration/problem_a.hpp
  └── owns the existing Step-4 Problem A operations and borrows the adapted app

source/adapted/step-4.cc
  └── owns the Step-4 mesh, FE, assembly, solves, and field output
```

The binding borrows `ProblemA`; the application keeps both objects alive until
the solver and output have finished. It has no instrumentation, native
reference optimizer, oracle, comparison logic, or run manifest. The current
callback contract still requires all five executable operations even though
this reduced evaluation consumes only the objective, objective derivative,
full VJP, and solve services during a successful optimization.

The control is the full 289-entry algebraic vector used by Step-4, and the
consumer uses the identity metric. It therefore does not include Problem B's
free-coordinate control or mass metric.

## Validation-only code

`tests/application/external_step4_minimal_problem_a_contract.cc` is outside
the consumer path. It invokes the new binding and executable behavior while
using the existing native reference and independent dense Problem A oracle to
check matrix/RHS reuse, matched optimization, fresh state/adjoint/gradient
audits, control distance, and output consistency. Its evidence is disposable
and is written below the ignored `runs/` tree. The evaluated `evaluation/`,
`verification/`, and `diagnostics/` directories likewise do not belong to the
minimal consumer.

## Build and run

Build the application and its focused validation target with the existing
deal.II profile:

```bash
./build.sh build debug-dealii \
  --target nmopt_external_step4_minimal_problem_a
./build.sh build debug-dealii \
  --target nmopt_external_step4_minimal_problem_a_contract_test
```

Run the consumer with an optional output path. Parent directories are created
by the application:

```bash
build/debug-dealii/bin/nmopt_external_step4_minimal_problem_a \
  runs/external-dealii/step-4/minimal/manual/problem-a-solution.vtk
```

Run the focused contract and executable checks:

```bash
ctest --test-dir build/debug-dealii --output-on-failure \
  -R '^nmopt\.external\.tutorial_step_4\.minimal_problem_a_binding_contract$'
ctest --test-dir build/debug-dealii --output-on-failure \
  -R '^nmopt\.external\.tutorial_step_4\.minimal_problem_a$'
```

The required profile gates for this deal.II/CMake change are:

```bash
./build.sh pipeline debug-dealii
./build.sh pipeline debug-neutral
```

The public callback, formulation, metric, lifetime, and solver contracts are
described in the [external deal.II integration reference](../../../../docs/reference/external-dealii-solver-integration.md).
The evaluated wiring and its limits remain documented in the [Step-4
README](../README.md) and the [boundary-evaluation roadmap](../../../../docs/planning/external-dealii-boundary-evaluation.md).
