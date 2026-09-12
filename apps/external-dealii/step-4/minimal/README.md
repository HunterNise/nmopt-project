# Minimal Step-4 consumers

This directory contains small, complete external applications that use the
current nmopt public contracts to optimize the existing Step-4 Problems A and
B. They are usability examples, not replacements for the evaluated
integration or for the Step-4 application itself.

## Problem A consumer

The A consumer's functional path is deliberately explicit:

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

The A binding borrows `ProblemA`; the application keeps both objects alive
until the solver and output have finished. It has no instrumentation, native
reference optimizer, oracle, comparison logic, or run manifest. The current
callback contract still requires all five executable operations even though
this reduced evaluation consumes only the objective, objective derivative,
full VJP, and solve services during a successful optimization.

The control is the full 289-entry algebraic vector used by Step-4, and the
consumer uses the identity metric. It therefore does not include Problem B's
free-coordinate control or mass metric.

## Problem B consumer

The B consumer retains the real FE distributed-control differences rather
than hiding them behind the A path:

```text
minimal/problem_b.cc
  ├── prepares Step4<2> and constructs ProblemB
  ├── constructs minimal/ProblemBBinding
  ├── creates zero full control and the frozen solver policy
  ├── runs ReducedSearchSolverT with the mass metric
  ├── reconstructs the physical state from free state coordinates
  └── sends the retained full state to Step-4's VTK writer

minimal/problem_b_binding.hpp
  ├── reuses ProblemB's free/full coordinate maps and FE mass/coupling
  ├── declares 225-state/289-control layouts and five callbacks
  ├── declares the state/control partition and native solve services
  ├── translates native CG evidence into LinearSolveReport
  ├── adapts ProblemB's mass metric through MetricT
  └── constructs the ReducedDTO

integration/problem_b.hpp, problem_b_coordinates.hpp,
integration/problem_b_mass.hpp, and integration/problem_b_metric.hpp
  └── own B's coordinates, lifting, distributed-control mathematics, and metric
```

The B binding borrows `ProblemB`; `ProblemB` in turn borrows the prepared
`Step4<2>` application. The application, problem, binding, metric, and solver
therefore remain alive in that order. The 289 control entries are full
algebraic coefficients, while the state and adjoint use 225 free coordinates.
The mass metric is applied and inverted through the existing native B metric
service; no A identity metric or boundary-coordinate control is substituted.

## Validation-only code

The two `tests/application/external_step4_minimal_problem_*_contract.cc`
files are outside the consumer paths. They invoke the new bindings and
executable behavior while using the existing native references and independent
dense A/B oracles to check matched optimization, fresh state/adjoint/gradient
audits, oracle control distance, and output consistency. Their evidence is
disposable and is written below the ignored `runs/` tree. The evaluated
`evaluation/`, `verification/`, and `diagnostics/` directories likewise do not
belong to the minimal consumers.

## Build and run

Build the application and its focused validation target with the existing
deal.II profile:

```bash
./build.sh build debug-dealii \
  --target nmopt_external_step4_minimal_problem_a
./build.sh build debug-dealii \
  --target nmopt_external_step4_minimal_problem_a_contract_test
./build.sh build debug-dealii \
  --target nmopt_external_step4_minimal_problem_b
./build.sh build debug-dealii \
  --target nmopt_external_step4_minimal_problem_b_contract_test
```

Run the consumer with an optional output path. Parent directories are created
by the application:

```bash
build/debug-dealii/bin/nmopt_external_step4_minimal_problem_a \
  runs/external-dealii/step-4/minimal/manual/problem-a-solution.vtk
build/debug-dealii/bin/nmopt_external_step4_minimal_problem_b \
  runs/external-dealii/step-4/minimal/manual/problem-b-solution.vtk
```

Run the focused contract and executable checks:

```bash
ctest --test-dir build/debug-dealii --output-on-failure \
  -R '^nmopt\.external\.tutorial_step_4\.minimal_problem_a_binding_contract$'
ctest --test-dir build/debug-dealii --output-on-failure \
  -R '^nmopt\.external\.tutorial_step_4\.minimal_problem_a$'
ctest --test-dir build/debug-dealii --output-on-failure \
  -R '^nmopt\.external\.tutorial_step_4\.minimal_problem_b_binding_contract$'
ctest --test-dir build/debug-dealii --output-on-failure \
  -R '^nmopt\.external\.tutorial_step_4\.minimal_problem_b$'
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
The functional/accessory source split and the A/B comparison are recorded in
the [minimal consumer assessment](../../../../docs/planning/review/external-dealii-boundary-evaluation/minimal-consumers-report.md).
