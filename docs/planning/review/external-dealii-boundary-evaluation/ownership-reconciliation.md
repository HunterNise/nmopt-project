# Step-4 ownership and evidence reconciliation

Status: closure record updated for the bounded post-G1 corrections on
2026-09-11.

This note makes the implementation boundary explicit after the original G1
report and records where corrected ignored artifacts are recreated. It does
not change the G1 conclusion, the shared nmopt API, or the frozen Problem A
comparison. The authoritative execution record remains the
[external deal.II boundary evaluation roadmap](../../external-dealii-boundary-evaluation.md).

## Minimum functional path

The smallest application-owned path exercised by the binding is:

```text
source/adapted/step-4.cc
        -> integration/problem_a.hpp
        -> integration/nmopt_binding.hpp
        -> existing nmopt public contracts
```

`source/adapted/step-4.cc` owns the deal.II application machinery. The
application-side `ProblemA` adapter exposes the residual, objective,
derivatives, native state/adjoint solves, and native output. The nmopt binding
constructs layouts, callbacks, solve services, the local identity metric, and
the reduced DTO using existing public contracts only.

Diagnostics are optional at runtime: `ProblemA` and the identity metric accept
an instrumentation pointer, while `NmoptBinding binding;` constructs and
solves without a diagnostics object. Instrumented constructors remain
available to the evaluation tests. `NmoptBinding` is explicitly non-copyable
and non-movable because its callbacks capture its owned `ProblemA`.

## Supporting code

The following code is intentionally outside the minimum wiring path:

| Location | Responsibility |
| --- | --- |
| `source/upstream/` | Verbatim external provenance. |
| `source/baseline/` | Comment-stripped fidelity comparison input. |
| `evaluation/` | Native reduced and optimization reference paths and frozen policy. |
| `verification/` | Deterministic scenarios, independent oracle, and comparison checks. |
| `diagnostics/` | Evaluation-only counters and solve evidence for explaining the comparison. |
| `tools/external_dealii/` | Reusable source and forward-comparison utilities. |
| `tests/dealii/` and `tests/application/` | Evaluation and contract drivers; not binding code. |
| `runs/external-dealii/step-4/` | Ignored generated artifacts only. |

The application directory therefore does not contain multiple functional
versions of Step-4 or a large undifferentiated evaluation folder. The complete
directory map and runnable commands are in the
[Step-4 README](../../../../apps/external-dealii/step-4/README.md).

The minimum path does not construct diagnostics. The current integration
headers still include `diagnostics/instrumentation.hpp` to provide nullable
instrumentation pointers, so collection is optional at runtime while the
header remains a source dependency. Removing that dependency is separate from
this evidence closure.

## Evidence corrections retained from G1 review

- C1 records the deal.II solver-control monitored residual as the final solve
  residual and removes the extra post-solve matrix-vector product. Successful
  and failed solve records remain distinct; native exception evidence is
  preserved when deal.II reports `NoConvergence`.
- C2 instruments identity-metric `apply` and `inverse_apply` callbacks directly.
  The matched run recorded 1,657 `apply` calls, 829 `inverse_apply` calls,
  and 829 solver-reported metric solves. The reduced-evaluation path recorded
  no metric callbacks. The direct counts and their schedule interpretation
  are recorded in the [reconciled G1 report](g1-report.md).
- O1 moved files by ownership and updated CMake, tests, and current links
  without changing the upstream, stripped, or adapted source contents.
- O2 removed the requirement to construct diagnostics for the functional
  path, without adding a shared helper or framework interface.
- EC1–EC3 preserved failure evidence, repaired the numerical acceptance and
  runtime-count checks, and made forward comparison reject nonfinite data.
  Their generated comparisons, summaries, counters, and working attribution
  are recreated below the ignored `runs/external-dealii/step-4/` tree.

## Test organization decision

The three Step-4 drivers remain separate:

- the native deal.II driver owns reuse, mathematical, oracle, and native
  optimization scenarios;
- the application driver owns reduced-evaluation binding checks and their
  compact paired artifact; and
- the paired optimization driver owns native-versus-nmopt trace and output
  comparison.

They retain local repository-root discovery and artifact creation because the
artifact layouts, lifetimes, and evidence responsibilities differ. Extracting
those functions into `tests/support/` would create a Step-4-specific test
utility and hide the distinction between reference, comparison, and oracle
artifacts. No shared test helper was added.

The ignored run tree was not pruned or reorganized. Its semantic folders and
run IDs are artifact storage, not implementation ownership or roadmap-unit
names.

## Verification

The final implementation was checked with the existing `debug-dealii` and
`debug-neutral` profiles. The focused forward selection passed 3/3, the
complete deal.II pipeline passed 175/175, and the neutral pipeline passed
67/67. The G1 report records the source revisions, current source hashes,
commands, environment, and generated run locations; raw traces, counters,
and VTK files remain outside the repository and are recreated by those
commands.
