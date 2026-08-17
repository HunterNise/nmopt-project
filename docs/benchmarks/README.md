# Benchmark specifications

Benchmark documents define frozen Chapter 6 experiments. Each benchmark must
identify its problem recipe, data, mesh, discretization, formulation, solver
parameters, tolerances, outputs, provenance, and expected numerical evidence.

Benchmarks are not a second problem library. A benchmark selects one exact
configuration from a reusable recipe and must not introduce another PDE
lowerer or optimizer.

The [Chapter 6 benchmark suite roadmap](../planning/chapter-6-benchmark-suite-roadmap.md)
owns benchmark selection and sequence. These documents will own the detailed
scenario contracts once implementation begins.
