# B2 replication findings

## Scope and status

This report records the development investigation of E6.5.2, the Graetz-flow
boundary-control example shown in Table 6.2 and Figures 6.4--6.5 of the source.
It supplements the [frozen B2 benchmark contract](chapter-6.md#b2--e652-graetz-flow-boundary-control)
and the [source catalogue](../guides/chapter-6-numerical-examples.md#e652--graetz-flow-boundary-control).
It separates source facts, deductions from the published counts, framework
replacement choices, completed negative evidence, and hypotheses that still
need experiments.

The historical persisted comparison baseline is Debug deal.II run `003` at framework
revision `923e4b9`, using refinement 6 of the framework-native quadrilateral
mesh. A behavior-neutral rerun at revision `db5a07a` reproduced its objectives
and solver histories and additionally recorded structural,
objective-component, field-extrema, and derivative-norm evidence. Targeted
Debug experiments at revision `fdc75cf` subsequently screened the two boundary
interpretations, native and structured simplex meshes, centroid-split simplex
connectivity, diffusion values, and several observation-objective
realisations. A second campaign at revision `21982c6` screened positive
forcing, regularisation, fixed-step policy, and joint Figure--Table fits on
the source-oriented $160\times40$ simplex mesh with continuous trace control.
These are development experiments rather than reproduction evidence.
Generated outputs remain disposable; the conclusions and change boundaries
are recorded here.

A September 2026 release campaign at framework revision `df50946` then
regenerated the authoritative and retained development profiles, calibrated a
quantitative comparison against the source raster, tested the zero-forcing
target-transcription hypothesis, and completed a deliberately non-source
forcing/target factorial. Across the seven release run sets used in that
campaign comparison, all `37/37` manifest artifacts completed, passed their
derivative evidence, and postprocessed successfully, with no failed or pending
cases. The follow-up Debug audit at revision `038cd7f` added 17 fresh,
profile-labelled forward records across F0–F2. Together these campaigns close
the replication attempt negatively: no source-compatible interpretation
reproduces the Figure 6.5 no-control field. Further optimization work is not
part of this replication attempt.

The current assessment is:

- the stated equation, four observation/target combinations, fixed data,
  transport, and full-BFGS execution path are represented and executable;
- all four current cases reach the 100-iteration limit, with objective values,
  terminal relative gradients, and state ranges far from the published data;
- the largest field discrepancy is already present at zero control: the
  current uncontrolled state lies in $[1,1.11959]$, while Figure 6.5 shows a
  range of approximately $[1,7.22]$;
- boundary-aligned structured triangles and centroid-split triangles change
  this maximum and the initial objectives negligibly, so the missing magnitude
  is not explained by the tested cell shape or connectivity;
- the attached source references [187] and [205] confirm that the B2 volume
  right-hand side is zero and that “forced convection” refers to the prescribed
  velocity field, not a volume source;
- the book's equation (6.65) explicitly selects the ordinary-normal-minus-
  transport boundary condition. The diffusion-weighted total-conormal
  interpretation is therefore not a source ambiguity, even though it remains
  useful as a historical diagnostic screen;
- positive constant-forcing screens can recover some plotted magnitudes and
  objective orders, but they contradict the source and must not be treated as
  B2 replication candidates;
- calibrated source-raster comparison confirms that forcing $f=0.47009$
  closely reproduces the no-control field's displayed range and normalized
  shape, while none of the retained optimized fields reproduces all four
  panels' sign and spatial behavior; this remains a non-source diagnostic;
- the coupled interior-scaling ray reaches the displayed range near $s=8.5$,
  but its normalized raster correlation remains approximately $0.696$ and
  its normalized MAE approximately $0.375$;
- a calibrated outlet-only natural source reaches the displayed range at
  $g=13.0$, but its normalized raster correlation is approximately $0.700$
  and its cross-stream profile MAE is approximately $0.171$; and
- the wall-flux fingerprint is intentionally deferred because a separate
  source-region selector would be a new reusable framework decision;
- the zero-forcing target-transcription gate finds a common positive constant
  target near $15.60$ only after swapping the two printed constant rows, but
  this does not explain the parabolic rows or the no-control field;
- the closest four-row initial-objective reconstruction needs $f=0.64$, a
  constant target $20$, and the same row swap; these are three explicit source
  contradictions rather than a replacement specification;
- the source counts constrain aggregate boundary subdivision under standard
  $P_1$ assumptions, but do not determine interior mesh connectivity;
- the current uniform facewise metric makes the relative metric-gradient and
  coefficient-derivative norms identical, so this realization cannot
  distinguish the book's gradient convention; and
- Table 6.2 and Figure 6.5 contain consistency questions that must be retained
  as source ambiguities rather than fitted silently through a framework
  change.

B2 is therefore **framework-verified but not replication-verified**. The
completed campaign supports a bounded negative conclusion for the tested
homogeneous, coupled-scaling, and outlet-load interpretations. The best image
fit is the explicitly non-source volume-load diagnostic. No BFGS or
regularisation result is promoted as replication evidence, and the replication
attempt is closed unless new source evidence justifies reopening it.

## Source facts and omissions

The source fixes $\Omega=(0,4)\times(0,1)$, diffusion $\mu=0.1$,
regularisation $\beta=10^{-3}$, zero initial control $u_0=0$, fixed
temperature $1$, transport field $b(x)=(1.5x_2(1-x_2),0)$, two downstream
observation regions, two targets, linear triangular finite elements, and
BFGS. It reports 11,028 vertices,
21,653 triangles, 10,907 state and adjoint degrees of freedom, and 243 control
degrees of freedom. Table 6.2 reports initial and final objectives, accepted
iterations, objective reduction, and the terminal-to-initial gradient-norm
ratio for all four cases.

The source does not provide:

- mesh coordinates or connectivity;
- the boundary-node subdivision or the endpoint convention behind the 243
  control coordinates;
- the norm used for the reported discrete gradient;
- the constant BFGS step value, stopping rule, or full update safeguards;
- volume, boundary, and target quadrature rules;
- linear-solver choices and tolerances; or
- numerical arrays and plotting settings behind Figure 6.5.

The book does explicitly provide the volume right-hand side and boundary
operator. Equation (6.65) uses zero volume forcing and
$\partial_{n} y-(b\mathbin\cdot n)y$ on both controlled and outflow boundaries;
the ordinary normal derivative is not diffusion-weighted. References [187]
and [205] support this classification: [187], equation (37), also has zero
volume right-hand side, while [205] separates a distributed heat-source test
(Test 2) from a zero-volume-source boundary-control test (Test 3). These
omissions still prevent coefficient-wise parity, but they do not justify
fitting a nonzero B2 forcing or selecting total conormal as the source form.

At zero control, $b\mathbin\cdot n=0$ on the horizontal control walls, while
$b\mathbin\cdot n>0$ on the outlet. Thus $y\equiv1$ does not satisfy the
outflow condition, and a nontrivial uncontrolled state is compatible with a
zero volume right-hand side. The Figure 6.5 maximum therefore cannot be used
by itself to infer an omitted forcing term.

## Constraints from the published counts

For a conforming triangulation of a simply connected rectangle with $V$
vertices, $T$ triangles, and $B$ boundary vertices,

```math
T = 2V - B - 2.
```

The published $V=11028$ and $T=21653$ imply $B=401$. If the state uses one
nodal $P_1$ coordinate per vertex, the difference $11028-10907=121$ is exactly
the number of strongly constrained Dirichlet vertices. The fixed boundary is
one connected arc from the lower transition point around the inlet to the
upper transition point, so those 121 vertices imply 120 fixed-boundary edges.

Under the additional standard assumption that the control uses independent
continuous $P_1$ traces on the two disconnected downstream wall segments, 243
control nodes imply 241 control edges: each segment contributes one more node
than edge. The remaining boundary then contains $401-120-241=40$ outflow
edges. Thus the reported dimensions are compatible with the aggregate split

```math
(N_D,N_c,N_{\mathrm{out}})=(120,241,40).
```

This deduction depends on nodal endpoint and junction conventions, and the
odd control-edge total permits asymmetric top and bottom subdivisions. It
does not reveal interior coordinates, diagonal choices, local refinement, or
connectivity. Matching these counts with a synthetic mesh would demonstrate
compatibility, not recover the source mesh.

The current refinement-6 realization is structurally different:

| Quantity | Source | Current development baseline |
| --- | ---: | ---: |
| Cell type and count | 21,653 triangles | 4,096 quadrilaterals |
| Physical state vertices | 11,028 | 4,225 |
| Independent state coordinates | 10,907 | 4,128 |
| Control coordinates | 243 | 96 facewise constants |
| Fixed/control/outflow boundary faces | inferred $120/241/40$ | $96/96/64$ |
| Wings/full observation measure | exact $1.8/3$ | $1.78125/3$ |

The wings-measure error is caused by selecting complete cells from their
centres at the $x_2=0.3$ and $x_2=0.7$ interfaces. It is small compared with the
present field and objective discrepancies, but it must be controlled in a
source-oriented mesh.

## Numerical baseline against Table 6.2

The refreshed refinement-6 baseline reproduced run `003` exactly. Here
relative gradient means terminal norm divided by initial norm, matching the
Table 6.2 caption; it is not one minus that ratio.

| Case | Region / target | Iterations, source/current | Initial objective, source/current | Final objective, source/current | Reduction, source/current | Relative gradient, source/current |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| a | wings / constant | $59/100$ (limit) | $316.6661/0.865898$ | $3.5682/0.0241174$ | $98.87\%/97.21\%$ | $0.0250/(3.63\times10^{-8})$ |
| b | full / constant | $54/100$ (limit) | $192.8385/1.45870$ | $2.6368/0.0450412$ | $98.63\%/96.91\%$ | $0.0569/(2.49\times10^{-8})$ |
| c | wings / parabolic | $48/100$ (limit) | $29.2188/0.310643$ | $0.7826/0.0502506$ | $97.32\%/83.82\%$ | $0.0753/(2.37\times10^{-7})$ |
| d | full / parabolic | $87/100$ (limit) | $45.9996/0.315375$ | $0.8464/0.122553$ | $98.16\%/61.14\%$ | $0.0387/(1.79\times10^{-7})$ |

The discrepancy is not a single missing scale: the current-to-source ratio
varies from $0.27\%$ to $1.06\%$ for the initial objectives and from $0.68\%$
to $14.48\%$ for the final objectives. The parabolic cases also fail to reproduce
the source's objective reduction. Scaling the reported functional would not
repair those case-dependent differences.

Figure 6.5 provides a stronger early diagnostic than the optimized
objectives. It shows the source zero-control state ranging from approximately
$1$ to $7.22$; the current zero-control state ranges from $1$ to $1.11959$.
The optimized state ranges also disagree:

| Case | Source state range | Current state range |
| --- | ---: | ---: |
| a | about $[-0.420,5.71]$ | $[1,2.37914]$ |
| b | about $[-1.88,4.45]$ | $[1,2.59184]$ |
| c | about $[-0.588,2.86]$ | $[0.263395,1]$ |
| d | about $[-1.38,3.60]$ | $[0.447711,1]$ |

Because the uncontrolled state is independent of optimization policy,
boundary and PDE discretisation candidates should be tested before BFGS step,
stopping, or update details.

## Targeted discretisation and boundary screens

The first screen used refinement-4 native quadrilaterals and two
boundary-aligned $40\times10$ triangular meshes. The structured mesh divides
each rectangle along one diagonal; the centroid-split mesh adds one interior
vertex per rectangle and connects it to the four corners. All runs used the
ordinary-normal boundary interpretation and one BFGS iteration, because only
the zero-control state and initial objective were relevant.

| Case | Source | Native quadrilateral | Structured triangles | Centroid-split triangles |
| --- | ---: | ---: | ---: | ---: |
| Wings / constant | $316.6661$ | $0.911380$ | $0.874930$ | $0.874927$ |
| Full / constant | $192.8385$ | $1.458507$ | $1.458526$ | $1.458523$ |
| Wings / parabolic | $29.2188$ | $0.312065$ | $0.310995$ | $0.310997$ |
| Full / parabolic | $45.9996$ | $0.315493$ | $0.315441$ | $0.315443$ |
| Uncontrolled state maximum | about $7.22$ | $1.119225$ | $1.119505$ | $1.119507$ |

The boundary-aligned meshes correct the wings observation measure from
$1.875$ on the native refinement-4 mesh to the exact value $1.8$. This
explains the visible but still small change in the wings/constant objective.
Structured and centroid-split triangles differ by less than $8\times10^{-6}$
relatively in every initial objective. The tested cell shape, interface
alignment, and connectivity therefore cannot explain the source discrepancy;
exact count-matched connectivity is not justified by this evidence.

The retained refinement-6 fields also permit observation-objective variants to
be evaluated without rerunning the PDE. Interpolating the target into the
state space changes the objectives by at most about $0.05\%$, and consistent
mass lumping changes them by at most about $0.16\%$. An unweighted sum over
nodal coefficients produces much larger values, but with case-dependent
factors and without representing the stated $L^{2}$ tracking functional.
Neither the defensible target interpolation nor quadrature variants explain
Table 6.2.

The second screen historically varied diffusion and the boundary form for the
full/constant case on the native refinement-4 mesh:

| Boundary interpretation | Diffusion $\mu$ | Uncontrolled maximum | Initial objective |
| --- | ---: | ---: | ---: |
| Ordinary normal | $0.03$ | $1.032906$ | $1.496341$ |
| Ordinary normal | $0.3$ | $1.403171$ | $1.174268$ |
| Ordinary normal | $1$ | $2.261675$ | $0.459786$ |
| Total conormal | $0.3$ | $15.133817$ | $35.187812$ |
| Total conormal | $0.4$ | $7.680102$ | $7.143744$ |
| Total conormal | $0.4125$ | $7.220595$ | $6.062632$ |
| Total conormal | $0.415$ | $7.135224$ | $5.870977$ |
| Total conormal | $0.5$ | $5.110750$ | $2.217154$ |

The value $\mu=0.4125$ was obtained by interpolation to test the Figure 6.5
range; it is not inferred source data. Although it reproduces the plotted
maximum to within $0.001$, its initial objective is smaller than the
full/constant Table 6.2 value $192.8385$ by a factor of about $31.81$. It also
contradicts the stated $\mu=0.1$. This tuned match shows only that one scalar
can fit one plotted output; it does not define a coherent replication.

The attached book source removes the earlier ambiguity behind this screen:
equation (6.65) uses the ordinary-normal-minus-transport form explicitly.
Accordingly, the total-conormal rows above are historical diagnostics and are
rejected as B2 source interpretations. The current ordinary realization in
the compiler is the source-aligned boundary choice; the remaining difference
between its zero-control maximum and Figure 6.5 is an unresolved source or
implementation discrepancy, not a reason to add volume forcing.

At the stated $\mu=0.1$, the total-conormal diagnostic is materially
incompatible in the opposite direction: a small comparison case produced an
uncontrolled maximum about $1393.5$ and an initial objective about $259295$.
The refinement-4 wings/constant derivative check had relative central
finite-difference mismatch about $1.46\times10^{-8}$ and Taylor order about
$2.002$, but its absolute mismatch $0.0137$ exceeded the fixed $10^{-7}$
evidence threshold because the objective and derivative were of order
$10^{5}$ and $10^{6}$. This is a scale-robustness defect in the development
evidence gate, not evidence of an incorrect derivative.

The direct weak-form audit is now complete. The assembled ordinary realization
uses the conservative volume form
$\int_{\Omega}\mu\nabla y\mathbin\cdot\nabla v-yb\mathbin\cdot\nabla v$
and adds $(1-\mu)(b\mathbin\cdot n)y$ on the control and outflow faces. Its
control load is correspondingly scaled by $\mu$, as required when the source
condition is $\partial_{n}y-(b\mathbin\cdot n)y=u$. On the B2 rectangle, a
constant test/state check gives
$(1-0.1)\int_{0}^{1}1.5x_{2}(1-x_{2})\,dx_{2}=0.225$ for the
ordinary-minus-total residual difference. The Debug contract
`nmopt.application.dealii.b2_ordinary_transport_boundary_operator` now locks
this coefficient and the outlet integral; the existing realization comparison
also passes. No forcing or production weak-form change is justified by this
audit.

## Historical forcing and optimization magnitude campaign

The second campaign comprised 17 complete Debug run sets and 92 artifacts at
revision `21982c6`. It used the boundary-aligned $160\times40$ structured
simplex mesh, continuous $P_1$ trace control, order-three analytic observation
quadrature, and the ordinary-normal boundary form unless stated otherwise.
All manifests completed successfully. Debug timings are not benchmark
evidence. Following the source audit, this campaign is retained only as
historical negative evidence: [187], [205], and the book all distinguish
volume forcing from the boundary-controlled B2 formulation, whose volume
right-hand side is zero.

Increasing ordinary-normal diffusion cannot reproduce Figure 6.5: the
uncontrolled maximum approaches a value near $5.07$ by $\mu=100$, still below
the plotted value about $7.22$. Positive constant forcing has the required
effect. Its zero-control response was:

| Forcing $f$ | Uncontrolled state range | Initial $J$, wings/constant | Initial $J$, full/constant | Initial $J$, wings/parabolic | Initial $J$, full/parabolic |
| ---: | ---: | ---: | ---: | ---: | ---: |
| $0$ | $[1,1.119609]$ | $0.875017$ | $1.458699$ | $0.310938$ | $0.315372$ |
| $1$ | $[1,14.09577]$ | $45.9714$ | $75.5286$ | $65.2729$ | $103.266$ |
| $2$ | $[1,27.07193]$ | $206.789$ | $339.998$ | $245.956$ | $396.616$ |

Interpolation of this response gives $f=0.47009$, for which the uncontrolled
maximum is $7.21957$. This is an exact Figure 6.5 range fit but not a coherent
source interpretation: it contradicts the explicit zero forcing in equation
(6.65) and the related formulations in [187] and [205]. Its four initial
objectives are only $7.66096$, $12.5635$, $16.4355$, and $24.9966$.
Optimizing this candidate with the source-stated $\beta=10^{-3}$ and Armijo
globalization drives all four objectives below $0.062$ and reduces them by
more than $99.5\%$, again incompatible with Table 6.2.

A regularisation sweep rules out $\beta$ as the sole missing detail. For zero
forcing, lowering $\beta$ from $10^{-3}$ to $10^{-6}$ raises the constant-case
reductions to about $98.9\%$, but the parabolic reductions remain near
$84.9\%$ on the wings and $62.0\%$ on the full region. With $f=1$, however,
$\beta=10^{-2}$ and Armijo give reductions between $98.43\%$ and $99.00\%$.
The earlier fixed-step calibration identified step $0.05$ as the strongest
common candidate for comparing the published iteration counts.

At those counts, $f=1$, $\beta=10^{-2}$, and fixed step $0.05$ place every
reported quantity within one decimal order of Table 6.2:

| Case | Iteration | Initial $J$, source/candidate | $J$ at iteration, source/candidate | Reduction, source/candidate | Relative metric gradient, source/candidate |
| --- | ---: | ---: | ---: | ---: | ---: |
| a | $59$ | $316.6661/45.9714$ | $3.5682/0.831156$ | $98.87\%/98.19\%$ | $0.0250/0.04874$ |
| b | $54$ | $192.8385/75.5286$ | $2.6368/1.15928$ | $98.63\%/98.47\%$ | $0.0569/0.06185$ |
| c | $48$ | $29.2188/65.2729$ | $0.7826/1.45561$ | $97.32\%/97.77\%$ | $0.0753/0.08570$ |
| d | $87$ | $45.9996/103.266$ | $0.8464/1.05098$ | $98.16\%/98.98\%$ | $0.0387/0.01065$ |

This is the best all-case order-of-magnitude fit found, but both $f=1$ and
$\beta=10^{-2}$ contradict printed source values. A forcing near $0.65$ is a
stronger simultaneous Figure--Table compromise, especially for the parabolic
cases:

| Case | Initial $J$, source/candidate | $J$ at source iteration, source/candidate | Reduction, source/candidate | Relative metric gradient, source/candidate |
| --- | ---: | ---: | ---: | ---: |
| a | $316.6661/17.0244$ | $3.5682/0.328626$ | $98.87\%/98.07\%$ | $0.0250/0.04843$ |
| b | $192.8385/27.9462$ | $2.6368/0.453623$ | $98.63\%/98.38\%$ | $0.0569/0.06171$ |
| c | $29.2188/29.3729$ | $0.7826/0.666825$ | $97.32\%/97.73\%$ | $0.0753/0.08554$ |
| d | $45.9996/45.5753$ | $0.8464/0.470078$ | $98.16\%/98.97\%$ | $0.0387/0.01050$ |

Its uncontrolled maximum is $9.55411$, about $32\%$ above the Figure 6.5
maximum. The parabolic initial objectives differ from Table 6.2 by only
$0.5\%$ and $0.9\%$; the wings final objective and relative gradient differ
by about $15\%$ and $14\%$. The full/parabolic terminal gradient remains the
largest optimization-side mismatch.

The constant rows show a separate, striking reconstruction clue. At $f=0.65$,
recomputing the tracking term with constant target $20$ and swapping the two
constant-region Table labels gives initial objectives $187.788$ and $314.284$,
only $2.6\%$ and $0.75\%$ from the swapped published values. Independently
solving for forcing gives $f=0.621019$ and $0.641707$ for those constant cases,
and $f=0.648059$ and $0.653155$ for the two parabolic cases. This clustering
near $f\approx0.64$ motivates a possible implementation or transcription
hypothesis: nonzero forcing near $0.64$, a constant target near $20$, and
swapped constant-case labels. Every part of that reconstruction contradicts
the printed source, so it is a candidate for further testing rather than a
replacement benchmark contract.

Four tracked development families retain the strongest distinct historical
diagnostics, not source-replication hypotheses:

- `figure-6.5-state-fit.prm` uses $f=0.47009$ to reproduce the uncontrolled
  plotted maximum;
- `table-6.2-order-fit.prm` uses $f=1$ for the strongest all-case decimal-order
  agreement; and
- `figure-6.5-table-6.2-parabolic-fit.prm` uses $f=0.65$ for the strongest
  simultaneous field and parabolic-row agreement; and
- `target-transcription-gate.prm` holds the source PDE and zero forcing fixed
  while comparing targets $2$, $20$, and $4x_{2}(1-x_{2})$ at zero control.

All four use continuous trace control, the source-oriented structured simplex
mesh, fixed step $0.05$, and retained fields. The three fitted optimization
families use $\beta=10^{-2}$ and 100 steps; the target-transcription gate keeps
the source $\beta=10^{-3}$ and retains one step only to serialize the complete
artifact contract. Their manifests and parameter provenance preserve every
source contradiction. They must not be promoted to B2 reproduction evidence.

The retained executions are canonical historical development runs `004`,
`005`, and `006`, respectively. Each manifest is complete with four successful
artifacts, no failures or pending artifacts, the expected parameter hash, and
framework revision `21982c6`. Every artifact retains its volume and boundary
fields. The source-iteration values above were reproduced exactly by these
runs.

The retained extrema already show that none reproduces the four optimized
Figure 6.5 panels, even before rendering:

| Case | Source | Run `004`, $f=0.47009$ | Run `005`, $f=1$ | Run `006`, $f=0.65$ |
| --- | ---: | ---: | ---: | ---: |
| a | about $[-0.420,5.71]$ | $[1,2.62764]$ | $[1,3.19042]$ | $[1,2.82049]$ |
| b | about $[-1.88,4.45]$ | $[1,2.48552]$ | $[1,2.89706]$ | $[1,2.62692]$ |
| c | about $[-0.588,2.86]$ | $[0.124623,1.41008]$ | $[-0.206220,1.93165]$ | $[0.0125864,1.58517]$ |
| d | about $[-1.38,3.60]$ | $[0.141438,1.40618]$ | $[-0.364096,1.90369]$ | $[-0.0299511,1.57425]$ |

Run `005` is closest in sign and magnitude for the optimized parabolic fields,
but its uncontrolled maximum is $14.0958$, nearly twice the plotted value.
Run `006` remains the best Table 6.2/parabolic compromise, not a Figure 6.5
reproduction.

The runner now exposes the constant and parabolic target definitions as
expression-backed scalar-function records. The source defaults remain $2$ and
$4x_{2}(1-x_{2})$; the selected definition, kind, value, and expression are
retained in B2 artifact evidence. This enables the conditional constant-target
transcription hypothesis to be tested without changing the PDE, adding
forcing, or introducing an objective multiplier.

Runs `004`--`006` were generated before the self-contained parameter-file
policy. Their snapshots retain top-level `include` directives, so their
run-set postprocessing remains blocked even though their numerical and native
field evidence is complete. Tracked parameter files are now required to be
self-contained; fresh candidate runs will therefore produce portable
snapshots. Existing generated snapshots must not be edited in place because
doing so would invalidate their provenance.

Versioned reruns `004-v2`, `005-v2`, and `006-v2` were then generated from the
self-contained candidates with the Debug runner at framework revision
`211ccfa`. Each run has four successful artifacts and a complete manifest. The
postprocessor now accepts the continuous trace's point-valued boundary output
as well as facewise cell data; all three reruns report `4/4` processed artifacts,
the eight expected PNG outputs, and a generated Chapter 6 report. These are
fresh diagnostic evidence and do not replace the historical `004`--`006`
snapshots.

### Regenerated run correspondence

The September 2026 profile rebuild regenerated the tracked B2 profiles with
the release runner and descriptive output names. The mapping to the retained
numbered runs is:

| Current named run | Tracked parameter file | Numbered predecessor | Numeric comparison |
| --- | --- | --- | --- |
| `figure-6.5-state-fit` | [figure-6.5-state-fit.prm](../../parameters/chapter-6/b2/development/figure-6.5-state-fit.prm) | `004`, `004-v2` | Identical numeric `artifact.kv` records for all four artifacts. |
| `figure-6.5-table-6.2-parabolic-fit` | [figure-6.5-table-6.2-parabolic-fit.prm](../../parameters/chapter-6/b2/development/figure-6.5-table-6.2-parabolic-fit.prm) | `006`, `006-v2` | Identical numeric `artifact.kv` records for all four artifacts. |
| `forcing-sweep` | [forcing-sweep.prm](../../parameters/chapter-6/b2/development/forcing-sweep.prm) | none | No numbered run contains a matching forcing-sweep parameter snapshot. |
| `table-6.2-order-fit` | [table-6.2-order-fit.prm](../../parameters/chapter-6/b2/development/table-6.2-order-fit.prm) | `005`, `005-v2` | Identical numeric `artifact.kv` records for all four artifacts. |

The new authoritative output is `runs/chapter-6/b2/authoritative/`; its
preserved pre-refactor counterpart is
`runs/chapter-6/b2/development/authoritative-20260818-f1a32e9/`. Those four
authoritative artifacts also have identical numeric records. Runs `001` and
`002` have no parameter-file snapshot, while `003` used the then-authoritative
configuration rather than a current named development profile. The historical
numbered runs and the regenerated named runs use different build profiles in
some cases, but the mapped B2 numeric records remain identical.

## September 2026 release evidence campaign

The release evidence campaign used the rebuilt `release-dealii` runner at
framework revision `df50946`. Its seven retained run sets comprise the
authoritative matrix, the three named fitted profiles, the forcing sweep, the
six-case zero-forcing target-transcription gate, and the 12-case forensic
forcing/target factorial. Their manifests report `37/37` successful artifacts,
zero failures or pending cases, `37/37` passed derivative records, and `37/37`
successful postprocessing records. Every artifact stopped at its configured
iteration limit: 100 accepted iterations in the optimization and forcing
profiles, and one retained iteration in the two initial-objective gates.

### Figure 6.5 raster evidence

The source comparison decodes the tracked `source-page-190.png` image with
SHA-256
`2cddb01d9f53bb7da2eaaa6d17e5a22b8e006704ee491d3301770ec2fe5246f7`.
It calibrates each panel independently against its printed colorbar and
compares the reconstructed raster with native finite-element fields in source
coordinates. The audit examined 144 retained artifact paths, representing 63
unique panel/field comparisons after duplicate field hashes were collapsed.

The no-control panel sharply separates the literal and fitted inputs:

| Candidate | Native state range | Normalized raster correlation | Normalized MAE | Classification |
| --- | ---: | ---: | ---: | --- |
| zero forcing | $[1,1.119597]$ | $0.6969$ | $0.3722$ | source literal |
| $f=0.47009$ | $[1,7.219451]$ | $0.9986$ | $0.0417$ | fitted, non-source |
| $f=0.65$ | $[1,9.553958]$ | $0.9987$ | $0.0399$ | fitted, non-source and too large |

The strong $f=0.47009$ match identifies a positive volume load, or
mathematically equivalent unprinted PDE or boundary data, as the best
descriptive explanation found for this panel. It does not identify source
data because equation (6.65) explicitly has zero right-hand side.

The optimized panels cannot be explained by rescaling the same fields. For
cases a--d, the source-literal normalized correlations are respectively
$0.368$, $-0.092$, $0.246$, and $0.554$. Their signs, peak locations, and
streamwise trends also differ. The printed panel b has about $11.28\%$
negative decoded pixels, while none of the 13 unique retained panel-b fields
crosses zero. The no-control fit is therefore not a hidden four-panel
reproduction.

### Zero-forcing target-transcription gate

The promoted [target-transcription gate](../../parameters/chapter-6/b2/development/target-transcription-gate.prm)
uses the source-oriented $160\mathbin\times40$ simplex mesh, continuous
$P_{1}$ trace control, zero forcing, source coefficients, and zero initial
control. It compares constant targets $2$ and $20$ with the source parabolic
target. The preserved release run completed all six artifacts and
postprocessed all six without failures.

| Target | Wings $J(u_{0})$ | Full $J(u_{0})$ | Full minus wings |
| --- | ---: | ---: | ---: |
| constant $2$ | $0.875017$ | $1.458699$ | $0.583682$ |
| constant $20$ | $324.413590$ | $540.695200$ | $216.281610$ |
| $4x_{2}(1-x_{2})$ | $0.310938$ | $0.315372$ | $0.004433$ |

Every computed pair satisfies the required nesting order. Solving the
constant-target quadratic for each printed objective gives incompatible
positive roots $19.7719$ and $12.3525$ under the printed row association.
Swapping only the two printed constant rows gives roots $15.6520$ and
$15.5438$. Their midpoint, $15.5979$, predicts both swapped objectives within
$0.75\%$. This is a qualified table-only transcription clue: it changes the
association of printed rows, not the selected wings/full scenario, and it does
not explain the parabolic objectives or Figure 6.5.

### Forensic forcing/target factorial

The release factorial tested $f\in\{0,0.64\}$, both regions, and targets
$2$, $20$, and $4x_{2}(1-x_{2})$ in 12 one-step cases. The best tested table
reconstruction uses $f=0.64$, target $20$, and the swapped constant-row
association:

| Compared case | Candidate $J(u_{0})$ | Printed $J(u_{0})$ | Relative error |
| --- | ---: | ---: | ---: |
| wings/constant $20$ to row b | $189.5193$ | $192.8385$ | $-1.72\%$ |
| full/constant $20$ to row a | $317.1569$ | $316.6661$ | $+0.16\%$ |
| wings/parabolic | $28.5555$ | $29.2188$ | $-2.27\%$ |
| full/parabolic | $44.2697$ | $45.9996$ | $-3.76\%$ |

The next-best tested combination has maximum relative error above $64\%$, so
the clustering is real. It is nevertheless explicitly non-source: it jointly
requires nonzero forcing, target $20$, and a row reinterpretation, while its
no-control maximum $9.42435$ is $30.53\%$ above the displayed $7.22$. The
factorial is recorded as contradiction evidence and is not promoted as stable
configuration.

### Source-literal trace checkpoints

The authoritative release histories sampled at the four source iteration
counts remain far from the source. Table 6.2 reports one gradient ratio; the
artifacts retain both the metric-gradient and coefficient-derivative ratios,
which coincide for the uniform facewise control realization.

| Case | Count | Source/current $J(u_{0})$ | Source/current $J$ | Source/current reduction | Source/current metric ratio | Current coefficient ratio |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| a | 59 | $316.6661/0.865904$ | $3.5682/0.0232535$ | $98.87\%/97.31\%$ | $0.0250/4.600\mathbin\times10^{-5}$ | $4.600\mathbin\times10^{-5}$ |
| b | 54 | $192.8385/1.458708$ | $2.6368/0.0437242$ | $98.63\%/97.00\%$ | $0.0569/4.895\mathbin\times10^{-5}$ | $4.895\mathbin\times10^{-5}$ |
| c | 48 | $29.2188/0.310638$ | $0.7826/0.0500038$ | $97.32\%/83.90\%$ | $0.0753/1.843\mathbin\times10^{-4}$ | $1.843\mathbin\times10^{-4}$ |
| d | 87 | $45.9996/0.315369$ | $0.8464/0.122407$ | $98.16\%/61.19\%$ | $0.0387/2.126\mathbin\times10^{-7}$ | $2.126\mathbin\times10^{-7}$ |

These results establish documented non-reproducibility for the completed
joint campaign, but they do not justify ending all B2 investigation. The next
campaign returns to the no-control forward state and admits optimization
evidence only after a candidate matches its range, sign, peak, and normalized
shape.

## Published consistency questions

The source audit used the book's references [187] and [205]. Reference [187]
also writes the Graetz state equation with zero volume right-hand side. In
reference [205], Test 2 is a different distributed heat-source-control
problem, while Test 3 is a boundary-control problem with zero volume source.
Neither reference supports introducing a fitted volume forcing into B2.

Reference [205] does expose a historically relevant difference: its Test 3
uses diffusion-only Neumann data on the control walls and outlet, rather than
the book's ordinary-normal-minus-transport condition. That difference does
not explain Figure 6.5's no-control field. With zero control and forcing, the
constant field $y=1$ satisfies the Test 3 PDE, inlet value, and homogeneous
diffusion-Neumann conditions. In the cross-section-averaged model
$\bar y=C_{1}+C_{2}e^{2.5x_{1}}$, the outlet condition
$\bar y'(4)=0$ forces $C_{2}=0$. The state shown on the retained
[reference page](../guides/assets/chapter-6/reference-205-page-A2336.png) is
an optimized controlled state, not an uncontrolled-state record. The
ancestral natural-Neumann form is therefore analytically rejected as a
global explanation of the displayed range; adding a third production
boundary policy solely to run that zero-control check is not justified.

The wings observation region is a subset of the full downstream region. At
the common zero control and for the same target, a positively weighted
tracking integral must therefore satisfy

```math
J_{h,\mathrm{wings}}(u_0) \leq J_{h,\mathrm{full}}(u_0).
```

Table 6.2 instead reports $316.6661>192.8385$ for the constant target,
while its parabolic cases have the expected ordering. Possible explanations
include a case-label or transcription error, non-nested discrete observation
sets, or an unreported objective convention. The current evidence cannot
choose among them.

There is a second conditional inconsistency. If Figure 6.5's zero-control
colour range $[1,7.22]$ is the actual data range rather than a clipped plotting
range, the stated continuous tracking functional with target $2$ is bounded
above by about $24.52$ on the wings and $40.87$ on the full region. Both
constant-target initial objectives in Table 6.2 exceed those bounds. The
table and figure may therefore use an unreported discrete scaling or may not
describe identical data. Numerical candidates must retain this ambiguity
instead of treating objective magnitude alone as an acceptance target.

## Gradient interpretation

The refreshed artifacts record both the declared $L^{2}$-metric gradient and
the Euclidean norm of the reduced derivative coefficients. On the uniform
facewise control mesh, every control-face mass is equal. Consequently the two
absolute norms differ by the constant factor four, and that factor cancels in
the terminal-to-initial ratio. Both interpretations give the current relative
gradient values in the comparison table.

Table 6.2 therefore cannot distinguish the two interpretations under the
frozen realization. A nonuniform boundary mesh or continuous $P_{1}$ boundary
control gives a non-scalar mass matrix and can make the ratios differ. B2 now
records the selected control and compiled metric explicitly, so candidate runs
can make that comparison without a further artifact-schema change.

## Candidate experiments and change boundaries

| Priority | Candidate | Motivation and decisive evidence | Required change |
| ---: | --- | --- | --- |
| F0 | [205] diffusion-only natural Neumann | **Analytically rejected for the no-control field.** At zero forcing and control, $y=1$ satisfies the PDE and all boundary data, independently of the transport amplitude. | Do not add a production boundary policy for this check. Reconsider only for a separately motivated nonzero-control ancestry audit. |
| F1 | Ancestral inlet-only Dirichlet partition with the book outflow form | **Open as a forward-only forensic check.** Figure 6.4 fixes the source partition, but reference [205] fixes only the inlet and controls the horizontal walls. This combination can change the no-control state because it retains the book outflow condition. | Parameter-only through the upstream transition; label it non-source. |
| F2 | Historically motivated transport amplitudes and direction | **Open only under the book outflow form.** Zero, unit, source $1.5$, and referenced scaled amplitudes can localize whether the field gap enters through volume transport or its coupled outlet term. | Parameter-only through the transport expression; retain exact provenance. |
| 1 | Diffusion-weighted conormal alternative | **Rejected as a source interpretation.** Equation (6.65) explicitly uses the ordinary-normal-minus-transport form; the total-conormal screen is retained only as historical diagnostic evidence. | No further framework change for B2. Any independently scaled boundary-transport coefficient would be a new, explicitly non-source hypothesis. |
| 2 | Source-oriented triangular $P_{1}$ state mesh | **Screened.** Boundary-aligned structured and centroid-split meshes give nearly identical states and objectives; connectivity sensitivity is negligible at this scale. | No further change unless source connectivity becomes available. |
| 3 | Continuous $P_{1}$ boundary control | **Screened and retained in historical diagnostics.** The source states linear finite elements and $N_{u}=243$; the $160\times40$ realization has 242 trace controls and distinguishes the metric from coefficient geometry. It does not by itself repair the field scale. | No further change for the current diagnostic; exact odd source counts would require an asymmetric or imported mesh. |
| 4 | Boundary-aligned observation geometry | **Screened.** The aligned meshes recover exact measure $1.8$ but do not materially reduce the objective discrepancy. | No further change for structured meshes. |
| 5 | Constant-step and stopping candidates | **Screened.** Step $0.05$ with evaluation at the source counts gives the strongest common reductions; it does not repair the full/parabolic terminal gradient. | No further change for fixed-step runs. Initial full-BFGS inverse-Hessian policies remain a possible framework extension. |
| 6 | Metric versus coefficient gradient norm | **Screened with continuous control.** The two relative histories differ, but neither consistently resolves all four published ratios. | No further evidence-schema change; nonuniform or exact source topology remains optional. |
| 7 | Objective quadrature or coefficient scaling | **Partially screened.** Target interpolation and mass lumping are far too small; an unweighted coefficient loss is large but case-dependent and changes the stated objective. | Explicit positive quadrature order and analytic-versus-state-FE target selection are implemented through B2 and the runner; candidate parameter files and runs remain. Arbitrary objective scaling remains excluded. |
| 8 | Exact count-matched connectivity or mesh import | **Deprioritised.** Counts constrain boundary totals but not the interior mesh, and the two tested triangular connectivities are numerically indistinguishable. | Reconsider only with source connectivity or contrary mesh-sensitivity evidence. |
| 9 | Linear-solver tolerances | These are omitted by the source but are unlikely to explain the smooth zero-control field gap. | Parameter-only with the existing solve-policy entries. |

Galerkin discretisation is stated by the source and remains fixed during these
screens. Stabilisation should be introduced only if a source-oriented mesh
shows a numerical inadequacy, and then reported as a project alternative
rather than an inferred source choice.

## Reopened forward-state evidence gate

The completed joint campaign found no coherent literal realization of
equation (6.65), Figure 6.5, and Table 6.2. The positive-forcing family
recovers selected ranges and decimal orders, but remains an invalid
calibration of the printed PDE rather than evidence for omitted forcing.
The target-transcription gate and forensic factorial are now complete; no
objective convention or BFGS change is the next experiment.

B2 is reopened at the more basic no-control forward problem. The next bounded
campaign must use one-step development runs only as a serialization vehicle
and compare the recorded `state_uncontrolled` field before considering the
objective, adjoint, gradient, or optimized state. It begins with the existing
source-literal release baseline, then tests the parameter-expressible
ancestral boundary partition and historically motivated transport choices
under the book boundary form. Analytically constant natural-Neumann cases,
completed mesh/connectivity screens, targets, regularisation, and optimizer
settings are excluded from this first stage.

A forward candidate may advance only if it simultaneously:

- preserves a defensible source or explicitly labelled historical provenance;
- reproduces the displayed $[1,7.22]$ range to the precision supported by the
  source raster;
- matches normalized shape, sign, streamwise trend, and peak location rather
  than only one endpoint; and
- remains stable under one coarser/finer paired mesh check.

Only after this gate passes should the four initial objectives be evaluated
with the source targets and nested-region invariant. Optimization policy,
gradient convention, and control scaling remain a later stage. If the
parameter-only forward tests isolate a missing boundary coefficient or weak
form that the current framework cannot express, that need must be reported as
a separate architectural decision before tracked implementation. A fitted
forcing, independently scaled boundary term, or arbitrary objective
multiplier remains excluded from source-replication evidence unless new
primary evidence supports it.

## Forward-state replication handoff

The bounded forward campaign specified in the reopened plan was provisionally
executed through Units 0–4 on 2026-09-03, then followed by the F0 provenance
repair and F1–F2 Debug diagnostics. The fresh F0–F2 records use the current
`debug-dealii` executable at framework revision
`038cd7f59242567b7df27be6d60efd43029d0a9a`; the Release rebuild was not
repeated. The existing source-literal comparison point remains the historical
release artifact at `runs/chapter-6/b2/authoritative` and was not rebuilt.

The required framework change was to make the lower endpoint of
`Boundary/upstream transition` inclusive while retaining the strict upper
endpoint. This permits the exact inlet-only diagnostic at transition `0.0`:
the left edge remains fixed Dirichlet and the horizontal exterior faces are
controlled. The change, regression contract, and parameter-contract note are
in commit `d403dfc` (`fix(b2): support inlet-only boundary partition`). The
Debug build and full pipeline passed all 148 tests after that change.

### Experiments performed

Every successful new family below declared the six-case matrix
`observation-region=wings,full` by
`target-profile=constant-2,constant-20,parabolic-source`. Each completed
family produced six native artifacts, six derivative-evidence passes, and
six successful PNG post-processing results.

| Unit | Parameter family and run root | Forward change | Result |
| --- | --- | --- | --- |
| 1 | `b2-forward-inlet-only-book-outflow-release` | Source partition changed to inlet-only, transition `0.0`, printed ordinary outlet form | Completed `6/6`; uncontrolled range `[1, 1.119651]`, raster correlation `0.6971`, and normalized MAE `0.3720`; effectively identical to the literal baseline and rejected as the range explanation. |
| 2 | `b2-forward-transport-c0-book-release`, `c1`, `c10`, `c15`, and `cminus1p5` | Source partition restored, transport `c*x1*(1-x1)` with `c` equal to `0`, `1`, `10`, `15`, and `-1.5` | All five completed `6/6`; `c=0` gives `y=1` to `2.5×10^-12`, positive coefficients remain near maximum `1.115`–`1.117`, and `c=-1.5` reverses the trend with range `[0.00376,1]`; no coefficient passes the source gate. |
| 3 | `b2-forward-boundary-ordinary-c1p5-book-release` | Source partition, `c=1.5`, μ=`0.1`, ordinary-normal-minus-transport | Completed `6/6`; range `[1, 1.119607]`, correlation `0.6969`, normalized MAE `0.3722`; it reproduces the literal low branch, not Figure 6.5. |
| 3 | `b2-forward-boundary-total-conormal-c1p5-book-release` | Same data with the implemented total-conormal endpoint diagnostic | All six cases were rejected before serialization by the absolute reduced-gradient finite-difference guard. The existing typed boundary comparison passes and separates the total-conormal zero-control state by more than `10×` at refinement 2; the documented μ=`0.1` diagnostic reaches approximately `1393.5`, so this form is rejected as a source interpretation. |

The existing literal `c=1.5` point is the sixth transport reference. Its
uncontrolled range is `[1,1.119607]`, with source-raster correlation
`0.6969`; it agrees with the new ordinary-form `c=1.5` result. The negative
coefficient reverses the downstream trend as an implementation-sign check,
so the positive branch is not a sign-orientation artifact. None of the
parameter-expressible homogeneous forward cases reaches the acceptance range
`[0.98,1.02]` to `[7.00,7.45]`, normalized correlation `0.98`, and normalized
MAE `0.06` simultaneously. Since no candidate passed the forward gate, no
refinement partner was justified.

### Follow-up F0 — provenance repair

F0 allocated fresh, profile-labelled, one-case parameter files and immutable
run roots. The selected case was wings/constant-2 because the uncontrolled
state is independent of target and observation-region selection. Every fresh
run has one expected and one successful artifact, passed derivative evidence,
and a post-processing index with no comparison errors. The current F0 ledger
contains seven complete Debug records at revision `038cd7f`; earlier
mislabelled or exploratory records remain historical evidence only.

| Replacement | Profile/revision | Run result |
| --- | --- | --- |
| `b2-f0-debug-inlet-only-038cd7f` | `debug-dealii` / `038cd7f` | Range `[1, 1.119651]`, no-flip correlation `0.6971`, normalized MAE `0.3720`; rejected as the range explanation. |
| `b2-f0-debug-transport-{c0,c1,c10,c15,cminus1p5}-038cd7f` | `debug-dealii` / `038cd7f` | Maxima `[1.000000, 1.117482, 1.114968, 1.114573]` for `c=0,1,10,15`; `c=-1.5` gives `[0.003762,1]`. All are provenance-clean exclusions, not source candidates. |
| `b2-f0-debug-ordinary-c1p5-038cd7f` | `debug-dealii` / `038cd7f` | Range `[1, 1.119609]`, no-flip correlation `0.6969`, normalized MAE `0.3722`; it matches the historical Release baseline. |

The F0 analysis script reconstructs the raw native extrema, calibrated raster
metrics, profiles, provenance, and decisions from each manifest and native
field. Its output records seven complete Debug records and zero source-gate
passes. This repairs the evidence classification and supplies the baseline for
the final F4 closure.

The provisional negative result left the coupled-scaling boundary
normalization, missing or nonzero volume data, and inconsistency between the
printed problem and Figure 6.5 as the remaining plausible classes. F1 and F2
now bound those parameter-expressible load classes: coupled scaling fails the
shape gate, the outlet-only source fails the shape/profile gate, and the
constant volume load is the best image fit but contradicts the printed zero
volume forcing. Unit 5 (initial-objective checksum) and Unit 6 (optimizer
campaign design) were not launched because F4 did not pass.

### Follow-up F1 — coupled interior-scaling ray

F1 was run on 2026-09-03 using the existing `debug-dealii` executable at
framework revision `038cd7f`, because rebuilding Release was too costly for
this screen. These are diagnostic results only; they are not promotion-quality
source evidence. The seven self-contained parameter files scale both

$$
\mu=0.1s,
\qquad
b=s[1.5x_{2}(1-x_{2}),0],
$$

with the prescribed values

$$
s\in\{5,7,8,8.5,8.75,9,9.5\}.
$$

Each case uses the source-sized `160,40`
structured-simplex mesh, the ordinary boundary form, and one fixed step.

All seven fresh run roots completed with one expected/successful artifact,
derivative evidence, and a clean PNG post-processing index. The range bracket
is hit at $s=8.5$, but the spatial gate fails throughout the ray:

| $s$ | Native uncontrolled range | No-flip raster correlation | Normalized MAE | Peak $(x,y)$ |
| --- | --- | --- | --- | --- |
| $5$ | $[1, 2.0883]$ | $0.6964$ | $0.3737$ | $(3.976, 0.502)$ |
| $7$ | $[1, 3.5558]$ | $0.6962$ | $0.3746$ | $(3.976, 0.501)$ |
| $8$ | $[1, 5.3956]$ | $0.6960$ | $0.3751$ | $(3.975, 0.502)$ |
| $8.5$ | $[1, 7.2360]$ | $0.6959$ | $0.3753$ | $(3.975, 0.502)$ |
| $8.75$ | $[1, 8.7075]$ | $0.6958$ | $0.3754$ | $(3.975, 0.501)$ |
| $9$ | $[1, 10.9128]$ | $0.6958$ | $0.3756$ | $(3.975, 0.501)$ |
| $9.5$ | $[1, 21.8934]$ | $0.6957$ | $0.3759$ | $(3.975, 0.502)$ |

The streamwise trend remains positive and the peak stays near the source
location, but the normalized shape is far outside the acceptance gate
$\mathrm{corr}\ge 0.98$, $\mathrm{MAE}\le 0.06$. Streamwise, centreline, and source-peak
cross-stream profile correlations are also only approximately
$0.696$, $0.672$, and $0.260$, respectively. Therefore F1 has no candidate
for local or mesh refinement: no complete forward gate passed. A Release
confirmation would still be required even if a Debug case had passed.

The reproducible F1 ledger and analysis are retained as ignored artifacts in
`runs/analysis/b2-forward-state-replication/`:

- `analyze_f1.py` reconstructs provenance, native extrema, raster metrics,
  peak locations, and profile comparisons from the manifests and native fields.
- `b2-f1-coupled-scaling-ledger.csv` is the row-oriented evidence ledger.
- `b2-f1-coupled-scaling-analysis.json` contains the gate decisions and source
  reference metrics.

### Conditional F2 — normalized load fingerprints

Because F1 failed the shape gate, the plan's conditional F2 was audited. The
existing B2 parameter contract can express the normalized interior-volume
forcing diagnostic already recorded in the earlier campaign (best fit
approximately $f=0.47009$, correlation $0.9986$, normalized MAE $0.0417$). It
does not expose a fixed wall-flux offset or an outlet-only boundary load as a
separate source datum. The B2 helper currently distinguishes fixed Dirichlet
data, controlled Neumann flux, and the ordinary/total-conormal boundary forms;
those controls do not provide the missing F2 load fingerprints.

F2 therefore separated the parameter-expressible volume-load diagnostic from
the boundary-load cases. The branch now implements the reusable immutable
natural-boundary-source contract and binds the B2 source to the outflow region
for its outlet-only diagnostic. The wall-flux case remains deferred because
the B2 parameter binding does not yet select a separate wall source region.
The reusable contract provides:

- an explicit boundary source datum and boundary region, separate from
  Dirichlet lifting and the Neumann control, with a declared trace/boundary
  pairing;
- an explicit normal/orientation policy and a choice between ordinary
  $\partial_{n}y-(b\mathbin\cdot n)y$ and total-conormal boundary terms;
- distinct wall and outlet source locations, with provenance, scaling, units,
  and face-partition information persisted in the run manifest; and
- regression coverage for zero-source identity, constant wall source,
  outlet-only source, sign/orientation reversal, disjoint partitions, and
  scale-aware derivative evidence.

The outlet-only source reaches the source range at $g=13.0$ but fails the
normalized shape and profile gates. The volume diagnostic at $f=0.47009$ is
the strongest image fit, with native range $[1,7.2195725]$, raster
correlation $0.998648$, and normalized MAE $0.041690$, but it is explicitly
not a source reconstruction. This closes the forward replication attempt;
Units 5 and 6 remain out of scope unless new source evidence reopens it.

Detailed ignored evidence is consolidated in
`runs/analysis/b2-forward-state-replication/`:

- `b2-forward-ledger.csv` contains the baseline and Units 1–3 provenance,
  hashes, counts, extrema, raster metrics, profiles, and decisions.
- `b2-transport-sweep.json` contains the complete Unit 2 profiles and the
  existing literal sixth-point reference.
- `b2-boundary-form-endpoint.json` contains the Unit 3 ordinary artifact and
  total-conormal failure/endpoint evidence.
- `b2-source-sized-forward-confirmation.json` records the Unit 4 gate and
  negative conclusion.
- `analyze_f0.py` is the reproducible F0 analysis entry point.
- `b2-f0-provenance-ledger.csv` and `b2-f0-provenance-analysis.json` contain
  the fresh F0 provenance, hashes, native extrema, calibrated metrics, and
  classifications.

## Final forward closure

The F0–F2 follow-up was completed at framework revision `038cd7f` with the
existing `debug-dealii` executable. It produced 17 fresh one-case records:
seven F0 provenance checks, seven F1 coupled-scaling checks, and three F2 load
fingerprints. All had complete manifests, one successful artifact, derivative
evidence, native fields, and clean PNG post-processing.

The constant volume-load case with $f=0.47009$ is the most informative
diagnostic. It reaches native range $[1,7.2195725]$, raster correlation
$0.998648$, normalized MAE $0.041690$, and passing streamwise and centreline
profile gates. It is nevertheless not a B2 reconstruction because the source
specifies zero volume forcing. It is promoted as the explicitly labelled
[Figure 6.5 volume-load image-fit diagnostic](../../parameters/chapter-6/b2/development/figure-6.5-volume-load-diagnostic.prm).

The outlet-only source reaches the displayed range at $g=13.0$ but fails the
spatial-shape and profile gates. The coupled-scaling ray reaches the range near
$s=8.5$ but also fails the shape gate. The wall-flux fingerprint is deferred:
running it would require a separate reusable source-region selector, not a
B2-specific switch.

The final verdict is **framework-verified but not replication-verified**. The
tested forward interpretations do not reproduce Figure 6.5, and no optimizer
result is promoted as replication evidence. Units 5 and 6 are not launched.
The reproducible final handoff is retained in the ignored
`runs/analysis/b2-forward-state-replication/b2-f3-f4-handoff.md` record.
