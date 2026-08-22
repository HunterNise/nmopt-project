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
realisations. These are development experiments rather than reproduction
evidence. Generated outputs remain disposable; the conclusions and change
boundaries are recorded here.

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
- the diffusion-weighted conormal interpretation at the stated $\mu=0.1$
  makes the zero-control solution orders of magnitude too large, while tuning
  that interpretation to $\mu=0.4125$ reproduces the plotted maximum but not
  the corresponding Table 6.2 objective;
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
- an unambiguous ordinary-normal versus diffusion-weighted conormal
  interpretation;
- the norm used for the reported discrete gradient;
- the constant BFGS step value, stopping rule, or full update safeguards;
- volume, boundary, and target quadrature rules;
- linear-solver choices and tolerances; or
- numerical arrays and plotting settings behind Figure 6.5.

These omissions prevent coefficient-wise parity and make several apparently
simple numerical comparisons convention-dependent.

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

The second screen varied diffusion for the full/constant case on the native
refinement-4 mesh:

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

At the stated $\mu=0.1$, the total-conormal candidate is materially
incompatible in the opposite direction: a small comparison case produced an
uncontrolled maximum about $1393.5$ and an initial objective about $259295$.
The refinement-4 wings/constant derivative check had relative central
finite-difference mismatch about $1.46\times10^{-8}$ and Taylor order about
$2.002$, but its absolute mismatch $0.0137$ exceeded the fixed $10^{-7}$
evidence threshold because the objective and derivative were of order
$10^{5}$ and $10^{6}$. This is a scale-robustness defect in the development
evidence gate, not evidence of an incorrect derivative.

## Published consistency questions

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

The refreshed artifacts record both the declared $L^2$-metric gradient and
the Euclidean norm of the reduced derivative coefficients. On the uniform
facewise control mesh, every control-face mass is equal. Consequently the two
absolute norms differ by the constant factor four, and that factor cancels in
the terminal-to-initial ratio. Both interpretations give the current relative
gradient values in the comparison table.

Table 6.2 therefore cannot distinguish the two interpretations under the
current realization. A nonuniform boundary mesh or continuous $P_1$ boundary
control gives a non-scalar mass matrix and can make the ratios differ. The
evidence fields already support that comparison; no further post-processing
or artifact-schema change is required before those discretisations exist.

## Candidate experiments and change boundaries

| Priority | Candidate | Motivation and decisive evidence | Required change |
| ---: | --- | --- | --- |
| 1 | Diffusion-weighted conormal alternative | **Screened and rejected at the stated $\mu=0.1$.** Tuning $\mu$ can fit Figure 6.5 but fails Table 6.2 and changes stated source data. | No further framework change for the two coherent boundary interpretations. An independently scaled boundary-transport coefficient would be a new hypothesis. |
| 2 | Source-oriented triangular $P_1$ state mesh | **Screened.** Boundary-aligned structured and centroid-split meshes give nearly identical states and objectives; connectivity sensitivity is negligible at this scale. | No further change unless source connectivity becomes available. |
| 3 | Continuous $P_1$ boundary control | The source states linear finite elements and $N_u=243$; the current 96 facewise constants cannot test that space or its gradient metric. | Compiler/framework extension plus parameter files. |
| 4 | Boundary-aligned observation geometry | **Screened.** The aligned meshes recover exact measure $1.8$ but do not materially reduce the objective discrepancy. | No further change for structured meshes. |
| 5 | Constant-step and stopping candidates | The source attributes its high counts to an unstated constant step, while the current Armijo BFGS reaches a much smaller gradient ratio at the iteration limit. Screen fixed iteration counts and relative thresholds first. | Parameter-only for iteration/tolerance screens; framework extension for a genuinely unconditional fixed step. |
| 6 | Metric versus coefficient gradient norm | The book names a discrete gradient but not its Hilbert metric. The current uniform facewise realization makes both relative ratios identical. | Evidence is already present; differentiation needs the continuous or nonuniform discretisation above. |
| 7 | Objective quadrature or coefficient scaling | **Partially screened.** Target interpolation and mass lumping are far too small; an unweighted coefficient loss is large but case-dependent and changes the stated objective. | No framework extension is justified for the tested variants. |
| 8 | Exact count-matched connectivity or mesh import | **Deprioritised.** Counts constrain boundary totals but not the interior mesh, and the two tested triangular connectivities are numerically indistinguishable. | Reconsider only with source connectivity or contrary mesh-sensitivity evidence. |
| 9 | Linear-solver tolerances | These are omitted by the source but are unlikely to explain the smooth zero-control field gap. | Parameter-only with the existing solve-policy entries. |

Galerkin discretisation is stated by the source and remains fixed during these
screens. Stabilisation should be introduced only if a source-oriented mesh
shows a numerical inadequacy, and then reported as a project alternative
rather than an inferred source choice.

## Next evidence gate

The low-cost state and objective screens are complete. Neither literal
boundary interpretation at the stated diffusion, either triangular
connectivity, boundary-aligned observation geometry, target interpolation,
nor mass lumping gives one coherent account of Figure 6.5 and Table 6.2.

Continuous $P_1$ boundary control, gradient interpretation, and BFGS policy
remain necessary to compare optimized results and reported gradient ratios,
but none can alter the zero-control state or initial tracking objective. They
therefore cannot resolve the primary source inconsistency. Before implementing
those larger framework units, the investigation needs an explicit choice:

1. introduce and screen a separately parameterized coefficient multiplying the
   transport term in the natural boundary operator as an acknowledged project
   hypothesis, while keeping $\mu=0.1$; or
2. accept that the published initial data are not jointly reproducible from
   the stated model, retain the current mathematically literal baseline, and
   implement continuous control and solver variants only to study the
   optimization-side omissions.

The first choice requires a typed semantic, compiler, application, and runner
extension and can fit a field scale, but it has no present source provenance
and must not be promoted as the authoritative benchmark. Scaling only the
control load would not answer this question because the reported initial
control is zero. The second choice preserves the stated mathematics and gives
the remaining framework work a narrower, defensible purpose.
