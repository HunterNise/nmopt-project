# B2 replication findings

## Scope and status

This report records the development investigation of E6.5.2, the Graetz-flow
boundary-control example shown in Table 6.2 and Figures 6.4--6.5 of the source.
It supplements the [frozen B2 benchmark contract](chapter-6.md#b2--e652-graetz-flow-boundary-control)
and the [source catalogue](../guides/chapter-6-numerical-examples.md#e652--graetz-flow-boundary-control).
It separates source facts, deductions from the published counts, framework
replacement choices, and hypotheses that still need experiments.

The persisted comparison baseline is Debug deal.II run `003` at framework
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
- the source counts constrain aggregate boundary subdivision under standard
  $P_1$ assumptions, but do not determine interior mesh connectivity;
- the current uniform facewise metric makes the relative metric-gradient and
  coefficient-derivative norms identical, so this realization cannot
  distinguish the book's gradient convention; and
- Table 6.2 and Figure 6.5 contain consistency questions that must be retained
  as source ambiguities rather than fitted silently through a framework
  change.

B2 is therefore **executable but not replication-verified**. The completed
screens rule out several low-cost explanations, but do not identify one
coherent realization of both Figure 6.5 and Table 6.2.

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

Three tracked development families retain the strongest distinct historical
diagnostics, not source-replication hypotheses:

- `figure-6.5-state-fit.prm` uses $f=0.47009$ to reproduce the uncontrolled
  plotted maximum;
- `table-6.2-order-fit.prm` uses $f=1$ for the strongest all-case decimal-order
  agreement; and
- `figure-6.5-table-6.2-parabolic-fit.prm` uses $f=0.65$ for the strongest
  simultaneous field and parabolic-row agreement.

All three use continuous trace control, the source-oriented structured simplex
mesh, $\beta=10^{-2}$, fixed step $0.05$, and retained fields. Their manifests
identify them as Debug development runs and explicitly preserve the source
contradictions. They must not be promoted to B2 reproduction evidence.

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

## Published consistency questions

The source audit used the book's references [187] and [205]. Reference [187]
also writes the Graetz state equation with zero volume right-hand side. In
reference [205], Test 2 is a different distributed heat-source-control
problem, while Test 3 is a boundary-control problem with zero volume source.
Neither reference supports introducing a fitted volume forcing into B2.

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

## Next evidence gate

The low-cost state, objective, and optimization screens are complete. No
literal realization of the printed data gives one coherent account of both
Figure 6.5 and Table 6.2. The positive-forcing family recovered some decimal
orders, but the source audit now classifies that success as an invalid
calibration of the B2 PDE rather than evidence for omitted forcing.

The three four-case Debug run sets have completed and their traces
confirm the screened diagnostic comparisons. The direct weak-form and
outflow-boundary audit against equation (6.65) also passes. Target definitions
are now parameterized as scalar values or expressions without changing the
frozen defaults, so the next evidence gate is the explicit constant-target
transcription experiment, followed by objective conventions or named
full-BFGS initial-scaling policies.
A fitted volume
forcing, independently scaled boundary term, or arbitrary objective multiplier
remains excluded from source-replication evidence.
