# September 2026 build and run refresh

> **Historical evidence.** This report reviews the 24 September 2026 campaign
> at execution commit `18dd7ddba440efbe0f2d8e5d43b94e1806b9b5fd`.
> It was written after the runs. A later documentation-only publication commit
> does not change the revision recorded in their manifests or filenames.

Generated paths below refer to ignored `build/`, `runs/`, and
`tmp/old-runs/` trees intended for release assets. They are not part of the
Git source archive containing this report.

## Question and result

The campaign asked whether the refactored repository could be built cleanly,
run the Chapter 6 examples and external Step-4 consumers, and recover the
behavior represented by the earlier output archive. It also gathered the
inputs, logs, and comparisons needed for a release. The old archive is mixed
historical evidence, not a uniform baseline: its 211 Chapter 6 manifests
contain 200 complete runs, 10 failures, and one unfinished run, across 165
Debug and 46 Release profiles and 28 revision strings. The
[B1](../../studies/chapter-6/b1-replication.md) and
[B2](../../studies/chapter-6/b2-replication.md) studies remain the authority for
what the runs mean as benchmark evidence.

The result is strong but qualified. The campaign supplies direct
equivalent-input comparisons, including a focused include-dependency audit,
while several old/new status changes and one stopping-reason difference
prevent a blanket claim that every invocation behaved identically.

## Results at a glance

| Area | Observed result |
| --- | --- |
| Clean builds | Four profiles configured and built; no compilation failure. |
| Final tests | Four suites passed after replay: 74/74 for each neutral profile and 202/202 for each deal.II profile. |
| Step-4 | All five saved executable/check commands passed. |
| Chapter 6 | 213 Release run sets: 206 numerically complete; seven reproduced archived failures. |
| Derived output | 38 complete runs could not plot without retained fields; 30 further runs marked `passed` have comparison-panel errors. |
| Refactor comparison | 47 completed include-bearing pairs have matching checked objectives; four archived failed/running cases now complete, and one B1 stopping reason differs. |

The initial pipeline tests failed on a generated-run fixture that was absent
after the clean build. They all passed when repeated after replay. Numerical
completion, successful plotting, and agreement with the published B1/B2
sources are separate questions throughout this report.

## What was executed

The first-pass driver, now preserved byte for byte at
`runs/refresh_all_runs.py`, invoked `./build.sh pipeline <profile> --clean-first` for `debug-neutral`, `sanitize-neutral`, `debug-dealii`, and
`release-dealii`. It then ran the historical Chapter 6 plan with the Release
runner and checked the external Step-4 executables. The plan at
`runs/refresh-18dd7ddba440/plan.json` maps 211 archived run sets and two
additional tracked profiles to 213 destinations. Every new Chapter 6
manifest records `release-dealii` and the execution commit above.

The first replay met old `.prm` syntax that the current runner could not
parse. The failed tasks were resumed after translating inputs to
self-contained current syntax. The first-pass `summary.json` is therefore
an intermediate state, not the final campaign verdict. The final
`resume-1/summary.json` classifies 168 tasks as passed, 38 as
`postprocess-failed`, and seven as `failed-as-before`.

## Build and test timings

These are phase times reported by the saved clean pipeline and final test
logs, not the elapsed time of the full rerun campaign. “Configure” and
“generate” are CMake’s separate `Configuring done` and `Generating done`
messages; “build” is the helper’s `build elapsed` value after
`--clean-first`; test times are CTest’s `Total Test time (real)`.

| Profile | CMake configure | CMake generate | Clean build | Initial test | Final test after replay |
| --- | ---: | ---: | ---: | ---: | ---: |
| `debug-neutral` | 0.8 s | 0.0 s | 38 s | 8.20 s, 73/74 | 10.68 s, 74/74 |
| `sanitize-neutral` | 0.6 s | 0.0 s | 1 m 19 s | 15.53 s, 73/74 | 11.65 s, 74/74 |
| `debug-dealii` | 1.0 s | 6.4 s | 19 m 52 s | 43.18 s, 177/178 | 242.09 s, 202/202 |
| `release-dealii` | 1.8 s | 1.3 s | 45 m 37 s | 49.60 s, 201/202 | 44.11 s, 202/202 |

All four initial test stages failed only
`nmopt.tools.postprocess_configuration`. The
[postprocessing configuration test](../../../tests/tools/postprocess_configuration_contract.py)
copies `runs/chapter-6/b1/development/008/run-manifest.json` into a
temporary fixture. This hardcoded numbered run did not exist immediately
after the clean build. Once replay regenerated it, all four final test
stages passed. The test therefore depends on a particular generated
`runs/` directory structure and on execution order, rather than creating
its own manifest fixture. A future test cleanup should remove that
dependency; this report does not change the test.

The command logs are under `runs/refresh-18dd7ddba440/logs/` and
`resume-1/logs/`; the first-pass and resumed result ledgers give the command,
exit status, and destination of each task. Build and test elapsed times are
printed in those logs. The two preserved Python drivers describe the original
campaign layout; `runs/README.md` explains their old `tools/` and
`tmp/old-runs/` path assumptions.

## Input provenance and numerical comparison

The old run's `parameters.prm` snapshot and the new run's saved input are the
first comparison boundary. Many old inputs used `include`; the old manifest's
parameter hash covers only its top-level snapshot, not the included file.
The replay wrote self-contained `.prm` files under `runs/parameters/` and
recorded source and prepared-input hashes in
`runs/refresh-18dd7ddba440/parameter-migration.json`. The original named
files remain under `tmp/old-runs/parameters/`. The plan's `old_manifest` and
`destination` fields provide the path correspondence; matching names alone
are insufficient.

The `include-audit.json` examines 50 runs with historical include
provenance. For 47 completed old/new pairs, it checked all corresponding
objectives and found matches; the remaining three were failed on both sides.
The audit also records six ignored dependency edges without a corroborating
same-revision base snapshot, so those inputs have a narrower provenance
claim. Two B1 runs, `archived-062` and `archived-063`, initially inherited the
current relative-gradient policy instead of their historical absolute policy.
After correction, `include-corrections.json` records the original and
corrected input hashes, focused reruns, equal numerical records, and matching
native fields. The corrected outputs, rather than the first attempts, are the
current destinations.

The early B1 numbered runs `001`–`004` had no saved `.prm`; their new inputs
were reconstructed from command-line defaults and cannot establish parameter
syntax equivalence. In `archived-006`, the old and new B1 L-BFGS
`beta-1e-2` cases have effectively the same final objective, but the stopping
reason changes from `line_search_failure` to `gradient_tolerance` at
floating-point precision. The old run was Debug and the new run Release.
These limits matter when interpreting apparent equality or difference.

## Failures and status changes

The seven current numerical failures reproduce the archived failures, with
the same diagnostic families:

- `archived-091`, `092`, `107`, and `166` are B2 conormal or boundary screens
  whose reduced-gradient finite-difference evidence fails.
- `archived-194` and `195` are B1 objective-target screens whose reference
  method does not reach its stopping tolerance.
- `archived-205` is a B1 control-space screen requiring an independent
  control degree of freedom.

These are retained diagnostic attempts. They do not change the status of the
named authoritative B1 and B2 runs, and they should not be presented as new
failures introduced by the refresh. The archive also contains three old
failures and one unfinished run that now complete:

| Task | Historical outcome | Refreshed outcome | Interpretation limit |
| --- | --- | --- | --- |
| `archived-069` | B1 objective-target reference tolerance failure, Debug | Complete, Release | Build profile changed; status alone does not isolate a code change. |
| `archived-085` | B1 run left pending after one of two cases, Debug | Both cases complete, Release | The old record was unfinished. |
| `archived-167` | B2 transition at `x = 0.0` rejected, Release | Complete, Release | The saved input differs only by a redundant Release override; the endpoint acceptance behavior changed across revisions. |
| `archived-209` | B2 observation catalog ID unregistered, Debug | Complete, Release | Both code revision and build profile changed. |

`archived-167` is an observable behavior difference, even though the later
B2 forward-state handoff at
`tmp/old-runs/analysis/b2-forward-state-replication/b2-f3-f4-handoff.md`
describes the lower endpoint as accepted. The run comparison does not by
itself assign that change to the refactor rather than another intervening
correction. It must be excluded from any claim of exact old/new behavior
preservation. The broader [B2 replication result](../../studies/chapter-6/b2-replication.md)
remains framework-verified but not source-replication-verified.

## Derived-output limits

All 38 tasks classified `postprocess-failed` have complete numerical
manifests but `retain fields = false` in their saved input, so the fields
needed for plotting were absent. The resulting plot failures do not undo
the numerical solves; they do limit what can be viewed or compared from those
runs. A further 30 tasks labelled `passed` have nonempty
`comparison_errors` in their `postprocess-index.json`: 15 report a history
panel with no series, and 15 report incomplete comparison matrices. The
driver's `passed` status checks artifact and plot counts, not this separate
comparison-error field. For figure evidence, inspect the individual index;
do not infer a clean comparison from the campaign label alone.

The current and historical run READMEs explain where to start with concrete
B1/B2 plots, numerical `artifact.kv` records, solver traces, and native fields.
They also separate the B2 source-literal run from fitted but non-source
experiments. This refresh does not alter the benchmark conclusions recorded
in the study documents.

## Archive stewardship

Campaign-created backups were moved out of the historical archive. The
`archive-relocations.json` ledger records 209 operations: 144 input-file
moves into `archived-inputs/`, 62 byte-identical input backups deduplicated
against current `runs/parameters/`, and three directory moves into
`archived-attempts/`. Original execution-time paths remain in older ledgers. Separately,
`old-svg-pruning.json` records 100 removed non-mesh SVG plots (798,245,194
bytes), each checked against a same-stem PNG. Native mesh SVGs and two
unpaired report SVGs remain. Some historical postprocessing indexes still
list the pruned SVGs, so the pruning ledger is part of the archive's
provenance.
