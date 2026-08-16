# Chapter 6 integration and benchmark readiness review

## Status and authority

This document records the final static cross-batch review at `942ede6` on
2026-08-16. It synthesizes the completed P6.1, P6.2, P6.3, and P6.5 reviews
without reopening every implementation diff. Those reviews cover 40
implementation commits after `codex/main`. P6.4 remains an intentionally
unimplemented conditional extension and is not treated as a missing
prerequisite for the selected reduced-space or PDAS benchmarks.

This is review evidence and a remediation handoff, not a second status ledger.
The [implementation roadmap](../../implementation-roadmap.md) remains the sole
owner of mutable feature status. The
[Chapter 6 benchmark suite roadmap](../../chapter-6-benchmark-suite-roadmap.md)
owns benchmark selection and acceptance, while the
[Chapter 6 numerical-examples reference](../../../guides/chapter-6-numerical-examples.md),
[Chapter 6 numerical-methods guide](../../../guides/chapter-6-numerical-methods.md),
[interface specification](../../../design/interface-specification.md), and
[v1 semantic compiler](../../../implementation/v1/semantic-compiler.md) remain
authoritative for the numerical, method, formulation, and compiler boundaries.

The individual reviews record 16 open P1 defects and seven open P2 gaps across
the four implemented P6 batches. Those counts are not repeated as new defects
below. The cross-batch findings identify how the existing defects combine at
integration and benchmark boundaries.

Overall Chapter 6 status is **implemented in separate bounded prototypes but
not acceptance-complete or benchmark-ready**. No B0 harness, executable
benchmark scenario, source-sized run, or tracked numerical result exists at
the reviewed head.

| Cross-batch finding | Severity | Status at `942ede6` |
| --- | --- | --- |
| `C6-I1` – no compiler product closes the selected formulation-to-method path | P1 | blocked by P6.2-R2, P6.3-R1, and P6.5-R3 |
| `C6-I2` – supplied OTD, KKT, and PDAS do not preserve one valid typed operator contract | P1 | blocked by P6.2-R1, P6.3-R2/R3, and P6.5-R1/R2 |
| `C6-I3` – executable ownership is not transitive across the composition chain | P1 | blocked by P6.2-R3, P6.3-R4, and P6.5-R4/R5 |
| `C6-I4` – projection and complementarity cannot prove they use the same compiled box | P1 | open; blocks the B3 cross-method comparison |
| `C6-B1` – the benchmark and prerequisite recipe layers are unimplemented | P1 | not started |
| `C6-B2` – no Release benchmark evidence or reproducible result record exists | P2 | absent |
| `C6-I5` – the mutable roadmap completion claims conflict with the review evidence | P2 | open |

## Integrated capability assessment

The implementation has useful pieces at each layer, but the seams remain
incomplete:

| Layer | Available implementation | Integration status |
| --- | --- | --- |
| semantic/compiler | registered scalar DTO graphs, projection constraints, optional reduced Hessian, and a distinct supplied-OTD wrapper | no compiled KKT, complementarity, or PDAS product; supplied-OTD semantics and ownership remain incomplete |
| reduced methods | backend-neutral P6.1 direction, line-search, Newton, trust-region, stopping, and reporting policies | selected scalar services exist, but P6.1 review defects remain and no experiment runner composes them from a frozen scenario |
| supplied OTD | explicit value/JVP/VJP/solve system and canonical scalar deal.II realization | executable provenance is detached, the semantic request is still DTO-shaped, and owned-session teardown is unsafe |
| quadratic KKT | generic $Q$, $D$, $D^{\mathsf T}$, residual, transpose, conversion, and Krylov services | direct adapters only; compiler registration, typed pairing, supplied validity, and transitive ownership remain open |
| complementarity/PDAS | represented-dual classification, active selections, reports, and a serial cellwise realization | direct construction only; active subproblem semantics and algebra, compiler registration, and ownership remain open |
| benchmark infrastructure | benchmark roadmap, numerical source catalogue, and a `release-dealii` CMake preset | no recipe implementation, B0 runner, frozen scenario file, benchmark executable, serializer, measurement layer, or result artifact |

The backend-neutral and deal.II CTest scenarios remain valuable local contract
evidence. They do not satisfy the benchmark roadmap's system-level requirement
that one frozen semantic problem, discretization, algorithm configuration, and
output record pass through the full compiled path.

## C6-I1 – no compiler product closes the selected formulation-to-method path

**Severity:** P1 – KKT and PDAS results cannot be attributed to the same
resolved semantic decision and compilation manifest as the DTO problem they
are intended to compare with.

### Evidence

`CompiledProblemT` owns an executable model, metric, optional projection
constraint, state/adjoint solvers, optional reduced Hessian, manifest, and
lifetime owner. It can produce a reduced DTO service. The distinct
`CompiledSuppliedOTDProblemT` owns a supplied system and manifest.
`CompilationResultT` can return one of those two products.

There is no compiled KKT or PDAS product in that result, no corresponding
resolved-decision record, and no KKT/complementarity/PDAS manifest section.
P6.3 tests construct `ScalarDiffusionReactionKKT` or adapt a supplied system
directly. P6.5 tests then construct the complementarity and `PDASSolverT`
directly over that KKT product. The semantic compiler does not participate in
either path.

This is not merely a missing convenience factory. The omitted compiler seam
is where target eligibility, formulation provenance, pairing, multiplier
conversion, bounds, active-set assumptions, and lifetime must be checked and
recorded. B3 cannot demonstrate inactive-box agreement or active-set behavior
for one compiled problem when its reduced and PDAS paths are assembled from
unrelated public objects. B5 likewise cannot produce the required
source-consistent KKT manifest.

### Required outcome and tests

Keep reduced DTO, supplied OTD, quadratic KKT, and PDAS as distinct product
types, but construct every selected type from one resolved compiler decision
and owned session. Add structured KKT and complementarity/PDAS records rather
than hiding solver choices in display strings or a named PDE adapter.

Compile one canonical scalar request to each supported formulation product and
verify that the shared mesh, data, spaces, objective scaling, multiplier
conversion, bounds, and provenance records agree where they should. Reject an
unsupported formulation before direct object construction. Run the selected
method from the returned owner-bearing product after releasing the original
compiler inputs.

## C6-I2 – supplied OTD, KKT, and PDAS do not preserve one valid typed operator contract

**Severity:** P1 – the all-at-once composition chain manufactures or loses
properties required by the next solver layer, and the final active KKT action
is algebraically inconsistent.

### Evidence

The P6.2 executable supplied system does not carry the mandatory pairings,
quadrature/discretization declaration, multiplier convention, or typed DTO
comparison status as one inseparable contract. Its compiler request does not
actually supply an OTD weak formulation.

The P6.3 bridge accepts that weak input, freezes JVPs at zero, and declares
affinity, rank, kernel positivity, multiplier conversion, and symmetric-
indefinite compatibility without receiving typed evidence for those
properties. The generic KKT product also lacks the domain–range pairing needed
to justify MINRES.

P6.5 inherits that declaration, but its active subproblem does not construct
the specified free-coordinate restriction. It retains the complete primal,
appends active equality rows and multipliers, and drops those multipliers after
the solve. In that augmented product, `apply_d_transpose()` also overwrites a
selected base contribution with the active multiplier instead of adding both
terms, while `apply_kkt_transpose()` performs an addition. The product remains
marked symmetric despite those different actions.

The direct cross-path tests use the same undeclared assumptions or assemble a
dense matrix from the same defective action. They therefore cannot establish
the missing formulation and pairing facts independently.

### Required outcome and tests

Carry one typed declaration from the supplied weak blocks through KKT and
PDAS. It must include block roles, trial/test pairings, discretization and
quadrature provenance, affine/quadratic validity, signs, multiplier
conversion, symmetry pairing, and rank/kernel assumptions. The KKT adapter
must reject absent or incompatible evidence.

Implement P6.5's declared free-coordinate restriction and reconstruction.
Remove the augmented-row alternative or register it as a separate formulation;
if it remains, correct its additive transpose contribution and validate its
own pairing and assumptions. Add independent bilinear transpose identities,
DTO/KKT solution agreement, and lower/upper-active PDAS stationarity checks
through the compiled product.

## C6-I3 – executable ownership is not transitive across the composition chain

**Severity:** P1 – an apparently valid compiled or adapted method object can
outlive the deal.II mesh or callback state it executes.

### Evidence

The ownership defects form one continuous chain:

1. `CompiledSuppliedOTDProblemT` declares its lifetime owner after its system,
   so reverse member destruction releases the owner before destroying the
   callback-bearing system.
2. The P6.3 supplied-to-KKT bridge copies the system but does not retain the
   compiled product's lifetime owner.
3. The generic P6.5 metric representation captures a metric by reference.
4. An active P6.5 product captures its `ActiveSetKKTSubproblemT` through
   `this`, while `PDASSolverT` stores its product and complementarity by
   reference.

Local tests keep triangulation, model, metric, products, and solvers alive in
safe lexical order. That does not prove that the public copyable services are
self-owning. A compiled supplied system adapted to KKT and then PDAS can lose
its owner at multiple independent seams.

### Required outcome and tests

Define one transitive executable owner from compilation session through
supplied system, KKT adapter, complementarity representation, active
subproblem, and PDAS solver. Callbacks must capture owner-bearing immutable
state, not stack builders or borrowed polymorphic references.

Add one sanitizer integration scenario that compiles an owned-session target,
constructs every supported method layer, destroys all source wrappers and
external handles, and then executes value, JVP, VJP, KKT, transpose,
conversion, classification, active-subproblem, solve, and destruction paths.

## C6-I4 – projection and complementarity cannot prove they use the same compiled box

**Severity:** P1 – B3's projected-reduced/PDAS comparison can silently compare
different constraints or multiplier metrics.

### Evidence

The v1 compiled DTO product may own a `ConstraintT` that supports projection
in its metric. Its manifest has a projection-oriented constraint record. P6.5
introduces a separate `BoxBoundsT` and multiplier representation, constructed
from independent lower/upper values and conversion callbacks. There is no
shared compiled box-data product, bound identity/digest, or conversion witness
from which both services are derived.

Consequently a benchmark can compile one projection constraint and manually
construct a similar-looking complementarity object without proving identical
bound coefficients, control layout, cell ordering, metric realization, or
data provenance. P6.1 additionally exposes projected advanced directions that
its own selected scope excludes, so an unsupported projected solver can be
chosen for the comparison.

### Required outcome and tests

Compile one owner-bearing typed box-data product from the semantic constraint
and bound bindings. Derive both coefficientwise projection and cellwise
complementarity from that product, preserving exact layout, ordering, bound
digests, metric identity, and provenance in the manifest. Expose only the
projected reduced methods supported by P6.1.

For B3, assert that both algorithms consume the same compiled box token and
initial control. Compare objective, feasibility, final control, and declared
residuals; reject deliberately changed bound data, metric, ordering, or
control layout before either solve begins.

## C6-B1 – the benchmark and prerequisite recipe layers are unimplemented

**Severity:** P1 – no selected Chapter 6 experiment can currently satisfy the
benchmark contract, even where a local numerical method prototype exists.

### Evidence

Repository search finds Chapter 6 benchmark identifiers only in planning and
guide documents. `CMakeLists.txt` registers contract and compiler tests, but
no benchmark runner or benchmark scenario. There is no tracked implementation
for the Chapter 5 problem-library L0 recipe boundary or its L1 scalar and L2
Neumann recipes, which B0, B1, and B2 require.

The following required B0 artifacts are absent:

- a versioned frozen scenario/configuration record;
- a runner consuming a compiled problem and algorithm policy;
- deterministic manifest and report serialization;
- mesh, parameter-sweep, and initial-control selection;
- objective, gradient, KKT, feasibility, and action-count histories;
- timing and memory measurement definitions;
- field-output and result locations; and
- checks that the runner has not become a second PDE lowerer or optimizer.

The scenario choices are also not frozen. E6.5.1 still needs a declared
forcing replacement or recovered source value; E6.5.2 needs a selected
fixed-step or Armijo policy; and the framework-native E6.9.1 case needs a
declared mesh and an explicit record that its cellwise `DG0` control differs
from the source's continuous `Q1` control. The benchmark roadmap labels B3/B4
“desirable,” while the implementation roadmap's next sequence instructs the
agent to run E6.9.1/E6.9.2; the required queue therefore needs one explicit
selection decision.

### Required outcome and tests

After the relevant method remediation, implement the problem-library L0–L2
boundaries needed by the selected scenarios, then implement B0 as a thin
consumer of compiler and solver services. Freeze source replacements,
discretization differences, algorithms, tolerances, outputs, and provenance
before running B1 or B2. Explicitly decide whether B3 and B4 are current
acceptance targets or remain desirable follow-ups.

Test schema validation, deterministic serialization, manifest/scenario
identity, exact rerun configuration, output collision handling, and rejection
of unsupported formulations. A recipe configuration and its direct semantic
composition must compile to the same resolved decision.

## Benchmark readiness by scenario

| Benchmark | Status at `942ede6` | Blocking boundary |
| --- | --- | --- |
| B0 common harness | not started | no recipe implementation, scenario schema, runner, serialization, or measurements |
| B1 E6.5.1 distributed Laplace | method pieces only | B0/L0/L1 absent; forcing replacement unfrozen; P6.1 action-count, stopping, and compiled-Hessian evidence relevant to the required record remains open |
| B2 E6.5.2 Graetz boundary control | compiler target and method pieces only | B0/L0/L2 absent; four frozen scenarios and BFGS step policy absent; no field/history report |
| B3 E6.9.1 symmetric box Laplace | direct manufactured prototype only | P6.3/P6.5 compiled, pairing, subproblem, algebra, ownership, and shared-box gates; B0 absent; target selection unresolved |
| B4 E6.9.2 asymmetric box Laplace | not started | B3 plus spatially varying bound recipe/scenario; target selection unresolved |
| B5 E6.7.1 all-at-once Laplace | conditional and unselected | P6.3 acceptance and B0 first; P6.4 remains conditional on measured basic-solve inadequacy |
| B6 E6.7.2 diffusion-reaction | not started follow-up | B5 |

P6.4 is therefore correctly left unimplemented. It must not be activated by
this audit or used to postpone B0–B3. Only a selected B5 run that demonstrates
inadequate direct/basic serial behavior can open the bounded preconditioner
work.

## C6-B2 – no Release benchmark evidence or reproducible result record exists

**Severity:** P2 – Debug and sanitizer scenario counts cannot substantiate
numerical reproduction, timings, memory use, or source-sized convergence
claims.

### Evidence

`CMakePresets.json` provides a `release-dealii` configure/build/test profile,
but the source tree contains no benchmark target or measurement code for that
profile. The preset builds the same registered contract tests in optimized
mode. A preset description calling it a Chapter 6 timing profile does not
define a clock, warm-up policy, repeated-run policy, memory measurement,
hardware record, output schema, or source-sized scenario.

There are no tracked benchmark manifests, reports, histories, tables, field
outputs, or Release run summaries. The roadmap's Debug and sanitizer pass
counts demonstrate intended local regression coverage only. They are not B1
through B6 evidence and cannot support performance or source-trend claims.

### Required outcome and tests

Once B0 and a scenario are ready, run correctness profiles first and then the
source-sized configuration with `release-dealii`. Record compiler/dependency
versions, hardware, mesh and dimensions, complete scenario and compilation
manifests, objective/residual histories, method and solve counts, timing
definitions, and memory methodology. Treat timing and iteration results as
observations, not portable correctness tolerances.

Retain a small deterministic development scenario for regression and a
separate source-sized reproduction record. Do not record a Release result
until the executable was built from the reviewed head and its exact manifest
is stored with the report.

## C6-I5 – the mutable roadmap completion claims conflict with the review evidence

**Severity:** P2 – the project's sole mutable status ledger currently directs
benchmark work as though its prerequisite method gates were accepted.

### Evidence

The implementation roadmap marks P6.1 implemented, P6.2 complete, P6.3
complete, and P6.5 complete, records final scenario pass counts, and directs
the next agent to the benchmark gate. The four static implementation reviews
instead record 16 open P1 defects and seven P2 gaps, including compiler,
pairing, algebra, ownership, and selected-method correctness blockers.

The discrepancy is especially consequential because the roadmap declares
itself the sole owner of feature and acceptance status. A future agent reading
only that ledger can start B3 or B5 on products the reviews explicitly reject
for acceptance.

### Required outcome

After this review is accepted, update the roadmap statuses to distinguish
“implementation present” from “acceptance pending remediation,” link the five
review documents, and replace the direct benchmark handoff with the bounded
remediation/selection sequence below. Keep P6.4 conditional and preserve all
recorded exclusions.

## Remediation and benchmark activation order

1. Reconcile the implementation-roadmap status and benchmark selection with
   the accepted review evidence (`C6-I5`, `C6-B1`).
2. Close only the P6.1 findings exercised by B1/B2 first: action-count
   evidence, selected stopping behavior, and compiled Hessian verification.
   Do not make unrelated second-order variants a benchmark prerequisite.
3. Close the P6.2 semantic provenance and owned-session boundary before using
   supplied OTD in any cross-formulation result.
4. Register and validate the P6.3 DTO KKT product, typed pairing, solver, and
   lifetime; then repair the supplied bridge as a separate acceptance slice.
5. Implement the P6.5 restricted active subproblem, correct its algebra, and
   add compiled complementarity/PDAS ownership and diagnostics.
6. Compile one shared box-data product for projection and complementarity
   comparisons (`C6-I4`).
7. Implement the minimum L0–L2 problem recipes and B0 harness, then freeze B1
   and B2. Decide explicitly whether B3/B4 join the current acceptance gate.
8. Run B1 and B2 through Debug correctness and `release-dealii` reproduction;
   if selected, run repaired B3 and then B4 through the same artifact path.
9. Leave B5/B6 and P6.4 conditional. Activate one bounded P6.4 target only
   after a selected B5 measurement demonstrates the need.

This order avoids using benchmark trends as a substitute for unresolved local
mathematics or lifetime correctness. It also avoids broadening the review into
Stokes, stabilization, automatic OTD derivation, continuous-control box
semantics, measure constraints, or other excluded Chapter 6 families.
