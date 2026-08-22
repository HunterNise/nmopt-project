# B1 replication findings

## Scope and status

This report records the development investigation of E6.5.1, the distributed
Laplace-control example shown in Figures 6.2–6.3 of the source. It supplements
the [frozen B1 benchmark contract](chapter-6.md#b1--e651-distributed-laplace-control)
and the [source catalogue](../guides/chapter-6-numerical-examples.md#e651--distributed-laplace-control).
It records the evidence behind the benchmark's explicit replacement choices;
those choices remain project policy rather than recovered source facts.

The initial evidence was produced with Debug deal.II builds. The decisive
source-sized experiments were repeated with `release-dealii` and record
framework revisions from `4221c67` through `56bd47a`. The selected
source-oriented seven-case reproduction was then completed with
`release-dealii` at framework revision `631537a`. It verifies the frozen
project replacement choices; it is not a book-equivalent recovery of the
omitted source data. Generated outputs below `runs/` remain disposable
evidence.

The current assessment is:

- the mathematical problem, field identities, sign convention, and the
  qualitative Figure 6.3 method comparison are reproduced;
- the authoritative seven-case matrix completes with seven valid artifacts,
  seven successful PNG post-processing records, and no failures;
- continuous `P1` control, constant-forcing, triangular-mesh, and solver-policy
  candidates can be selected independently in parameter files;
- neither the source mesh nor the omitted forcing has been uniquely recovered;
- a relative-gradient threshold of $10^{-3}$ consistently selects the fourth
  L-BFGS iterate and reproduces the Figure 6.2 sign-changing
  $\beta=10^{-6}$ adjoint morphology;
- constants near $f=0.4$–$0.5$ give the strongest aggregate agreement, but
  different extrema favor different constants; and
- the Figure 6.3 iteration counts are effectively independent of the tested
  forcing and control representation, while their objective levels are not.

The benchmark is therefore **reproduction-verified under the frozen project
contract**, while parity with the undisclosed source realization remains
unresolved.

## Source facts and omissions

The source fixes the square domain, target, continuous linear triangular state,
adjoint, and control spaces, a zero initial control, and the broad optimization
methods. It reports 34,320 triangles, 17,361 vertices, and 16,961 independent
coordinates for each field. For Figure 6.3 it also gives the steepest-descent
relative-gradient threshold and the Armijo constants.

The source does not give:

- the forcing $f$;
- mesh coordinates or connectivity;
- quadrature and target-evaluation rules;
- linear-solver tolerances;
- L-BFGS memory, initial scaling, curvature safeguards, or a complete update
  policy;
- an unambiguous stopping rule for the Figure 6.2 field snapshots; or
- numerical arrays behind either figure.

These omissions prevent coefficient-wise or history-wise parity from being
established from the book alone.

## A consistency constraint from optimality

For the framework adjoint sign convention, unconstrained first-order
stationarity in the stated common continuous `P1` control and adjoint space is

```math
p + \beta u = 0.
```

The comparison-only book-convention field is $`p_{\mathrm{book}}=-p`$, so a
converged discrete solution must satisfy
$`p_{\mathrm{book}}=\beta u`$. The mass matrix does not weaken this conclusion:
the same invertible `P1` mass matrix multiplies both terms.

The top row of Figure 6.2 is consistent with this relation. For
$\beta=10^{-3}$, the plotted control maximum is approximately $8.57$ and the
book-convention adjoint maximum is approximately $8.57\times10^{-3}$, as
required by the factor $\beta$.

The bottom row is not consistent with stationarity. Its plotted control is
nonnegative with maximum approximately $8.90$, whereas its adjoint changes
sign over the domain and ranges approximately from $-1.55\times10^{-5}$ to
$3.31\times10^{-5}$. A converged solution of the stated discrete problem
cannot have both fields. At least one of the following must therefore hold:

- the plotted $\beta=10^{-6}$ fields are from an early iterate;
- the adjoint and control snapshots do not come from the same iterate;
- an unreported discrete control or gradient convention was used; or
- the source figure or caption contains a field-identification error.

The experiments strongly support the early-iterate explanation, although
they do not prove it.

## Experimental findings

### Figure 6.3 solver comparison

The recovered policy stops steepest descent at relative gradient norm
$10^{-3}$ and stops L-BFGS when it reaches the corresponding terminal
steepest-descent objective. The original zero-forcing framework-native run
gave steepest-descent counts `64/548/2202`, L-BFGS counts `2/4/4`, and matched
objectives `0.0541866/0.0443562/0.0158240` for
$\beta=10^{-1},10^{-2},10^{-3}$.

Two release families then used the source-oriented continuous `P1` control and
regular 131-subdivision triangular mesh:

| Forcing | $\beta$ | Initial objective | Steepest-descent iterations | L-BFGS iterations | Matched objective |
| ---: | ---: | ---: | ---: | ---: | ---: |
| $0.5$ | $10^{-1}$ | 0.0488993 | 64 | 2 | 0.0476791 |
| $0.5$ | $10^{-2}$ | 0.0488993 | 547 | 4 | 0.0389394 |
| $0.5$ | $10^{-3}$ | 0.0488993 | 2,175 | 4 | 0.0137935 |
| $0.4150674$ | $10^{-1}$ | 0.0500000 | 64 | 2 | 0.0487524 |
| $0.4150674$ | $10^{-2}$ | 0.0500000 | 547 | 4 | 0.0398164 |
| $0.4150674$ | $10^{-3}$ | 0.0500000 | 2,180 | 4 | 0.0141057 |

Both families reproduce the source trend and approximate raster locations:
steepest descent becomes dramatically slower as $\beta$ decreases, while
L-BFGS reaches the matched objective in a few iterations. The forcing changes
the vertical objective levels but not the convergence mechanism. The
objective-matched constant improves the visible common starting level and the
$\beta=10^{-2}$ plateau; $f=0.5$ appears closer for the other two terminal
levels. The raster is not precise enough to select between them.

Changing only the control representation in the earlier runs also preserved
the iteration counts. The Figure 6.3 histories therefore identify the broad
solver policy but not the forcing or control discretisation. The source
explicitly states continuous `P1`, so that remains the appropriate
source-oriented choice.

A final release check replaced the L-BFGS objective target with the same
$10^{-3}$ relative-gradient threshold used by steepest descent. For
$\beta=10^{-1},10^{-2},10^{-3}$ it selected the same `2/4/4` iterates, with
solver traces and final objectives identical to the objective-target runs. It
also selected iteration 4 for $\beta=10^{-6}$, reproducing the constant-half
Figure 6.2 candidate. This establishes numerical equivalence for the tested
family, not equivalence of the two stopping-policy declarations.

### Constant-forcing screen

Constant forcings $f\in\{0,0.5,1\}$ were screened with continuous control on
the refinement-7 framework mesh. The $\beta=10^{-3}$ fields use a relative
gradient tolerance of $10^{-5}$; the $\beta=10^{-6}$ comparison uses the
common screening tolerance $10^{-4}$. Here $`p_{\mathrm{book}}=-p`$.

| $\beta$ | Candidate | State maximum | Control maximum | $`p_{\mathrm{book}}`$ range |
| ---: | --- | ---: | ---: | ---: |
| $10^{-3}$ | Book | about 0.475 | about 8.57 | about $[0,8.57\times10^{-3}]$ |
| $10^{-3}$ | $f=0$ | 0.474663 | 9.04268 | $[0,9.04281\times10^{-3}]$ |
| $10^{-3}$ | $f=0.5$ | 0.482429 | 8.50349 | $[0,8.50359\times10^{-3}]$ |
| $10^{-3}$ | $f=1$ | 0.490195 | 7.96429 | $[0,7.96438\times10^{-3}]$ |
| $10^{-6}$ | Book | about 0.625 | about 8.90 | about $[-1.55,3.31]\times10^{-5}$ |
| $10^{-6}$ | $f=0$ | 0.626008 | 10.2539 | $[0,1.09134\times10^{-5}]$ |
| $10^{-6}$ | $f=0.5$ | 0.625707 | 9.78228 | $[0,1.09214\times10^{-5}]$ |
| $10^{-6}$ | $f=1$ | 0.625774 | 9.27675 | $[0,1.00212\times10^{-5}]$ |

For $\beta=10^{-3}$, zero forcing reproduces the state maximum almost exactly,
while a forcing near $0.5$ better reproduces the control and adjoint magnitude.
For $\beta=10^{-6}$, $f=1$ gives the closest tested control magnitude and all
three candidates reproduce the state maximum closely. None of the later
tolerance-terminated screens reproduces the sign-changing adjoint.

No single tested constant forcing explains both source rows at these later
solutions. In particular, the extrema do not justify promoting $f=1$ from a
motivated hypothesis to a source fact. The later screen remains useful for
separating forcing effects from the much larger stopping-state effect.

### Source-oriented four-iteration Figure 6.2 sweep

The completed source-oriented sweep combined continuous homogeneous-Dirichlet
`P1` control, the regular 131-subdivision triangular mesh, zero initial
control, the book Armijo constants, metric-inverse memory-5 L-BFGS, and a
$10^{-3}$ relative-gradient stop. Every tested constant stopped after exactly
four accepted unit steps.

| $\beta$ | Candidate | State maximum | Control maximum | $`p_{\mathrm{book}}`$ range |
| ---: | --- | ---: | ---: | ---: |
| $10^{-3}$ | Book | about 0.475 | about 8.57 | about $[0,8.57\times10^{-3}]$ |
| $10^{-3}$ | $f=0$ | 0.474561 | 9.03985 | $[0,9.04217\times10^{-3}]$ |
| $10^{-3}$ | $f=0.4150674$ | 0.481010 | 8.59257 | $[0,8.59461\times10^{-3}]$ |
| $10^{-3}$ | $f=0.5$ | 0.482329 | 8.50105 | $[0,8.50303\times10^{-3}]$ |
| $10^{-3}$ | $f=1$ | 0.490098 | 7.96224 | $[0,7.96388\times10^{-3}]$ |
| $10^{-6}$ | Book | about 0.625 | about 8.90 | about $[-1.55,3.31]\times10^{-5}$ |
| $10^{-6}$ | $f=0$ | 0.618022 | 9.60297 | $[-1.734,4.318]\times10^{-5}$ |
| $10^{-6}$ | $f=0.4150674$ | 0.618802 | 9.12043 | $[-1.529,3.866]\times10^{-5}$ |
| $10^{-6}$ | $f=0.5$ | 0.618962 | 9.02305 | $[-1.487,3.774]\times10^{-5}$ |
| $10^{-6}$ | $f=1$ | 0.619897 | 8.46315 | $[-1.241,3.233]\times10^{-5}$ |

For $\beta=10^{-3}$, zero forcing best matches the state, while
$f=0.4150674$ best matches the control and adjoint magnitude. For
$\beta=10^{-6}$, $f=0.5$ gives the most balanced control and adjoint range;
$f=0.4150674$ best matches the negative minimum, and $f=1$ best matches the
positive maximum. All four forcings reproduce the source's central positive
region, four negative lobes, and positive edge regions. The sign changes are
present in the native VTU fields rather than introduced by plotting.

At zero control, the discrete state depends affinely on a constant forcing,
so the initial objective is quadratic. The three computed values at
$f=0,0.5,1$ give

```math
J_{0}(f)=0.0008510022804 f^{2}-0.0137379298020 f+0.055555555556.
```

Matching the approximately `0.05` initial level visible in Figure 6.3 gives
roots `0.4150674` and `15.7282`; only the smaller root is plausible in the
screened range. This is a raster-matched diagnostic, not a recovered source
value.

The earlier run `006` first exposed the four-iteration explanation with zero
forcing and cellwise control. The source-oriented sweep shows that the same
explanation survives the stated continuous triangular discretisation. Fixed
iteration scans at $f=0.5$ and $f=1$ further localize the effect: iteration 3
has adjoint extrema an order of magnitude too large, iteration 4 is the first
to satisfy the $10^{-3}$ relative threshold, iteration 5 is nearly identical,
and iteration 6 begins moving away from the source-like range. The stopping
rule remains inferred because the source states it for steepest descent, not
for the Figure 6.2 BFGS snapshots.

### Mesh counts and connectivity

The published counts constrain the topology but do not determine the mesh.
For a conforming triangulation of a simply connected square with $V$ vertices,
$T$ triangles, and $B$ boundary vertices,

```math
T = 2V - B - 2.
```

The book's $V=17{,}361$ and $T=34{,}320$ therefore imply $B=400$, which also
explains the 16,961 independent homogeneous-Dirichlet coordinates. The counts
do not reveal the interior coordinates, diagonal choices, refinement history,
or connectivity.

| Mesh candidate | Cells | Physical vertices | Independent coordinates |
| --- | ---: | ---: | ---: |
| Book | 34,320 | 17,361 | 16,961 |
| Framework quadrilateral, refinement 7 | 16,384 | 16,641 | 16,129 |
| Regular structured simplex, 131 subdivisions | 34,322 | 17,424 | 16,900 |
| Count-matched centroid-split hypothesis | 34,320 | 17,361 | 16,961 |

The count-matched candidate begins with a 100-by-100 structured triangulation
and splits 7,160 selected triangles at their centroids. It proves that the
three counts can be matched, but its artificial connectivity is not evidence
that the book used that construction. A more natural interpretation is a
400-boundary-vertex mesh with unreported interior refinement or generation.

The regular structured-simplex $f=1$ experiment changed the field extrema by
less than 0.1% relative to the quadrilateral $f=1$ screen:

| $\beta$ | Mesh | State maximum | Control maximum |
| ---: | --- | ---: | ---: |
| $10^{-3}$ | Quadrilateral | 0.490195 | 7.96429 |
| $10^{-3}$ | Structured simplex | 0.490105 | 7.96375 |
| $10^{-6}$ | Quadrilateral | 0.625774 | 9.27675 |
| $10^{-6}$ | Structured simplex | 0.625750 | 9.27058 |

At this resolution, cell shape and the tested regular connectivity do not
explain the extrema discrepancy. Exact source connectivity could still affect
fine contour details, but recovering it would require the original mesh rather
than inference from counts.

### L-BFGS details

The omitted L-BFGS scaling and stopping choices materially affect the later
$\beta=10^{-6}$ solve, but memory does not affect the source-like fourth
iterate. Memory 5 and memory 20 with metric-inverse scaling produced identical
iteration-4 fields in the tested cases. Changing only the initial scaling to
scalar-secant produced the following tradeoff:

| Forcing | Scaling | Relative gradient | State maximum | Control maximum | $`p_{\mathrm{book}}`$ range |
| ---: | --- | ---: | ---: | ---: | ---: |
| $0.5$ | Metric inverse | $7.06\times10^{-4}$ | 0.618962 | 9.02305 | $[-1.487,3.774]\times10^{-5}$ |
| $0.5$ | Scalar secant | $6.00\times10^{-4}$ | 0.620112 | 8.87188 | $[-1.121,3.251]\times10^{-5}$ |
| $1$ | Metric inverse | $6.52\times10^{-4}$ | 0.619897 | 8.46315 | $[-1.241,3.233]\times10^{-5}$ |
| $1$ | Scalar secant | $5.54\times10^{-4}$ | 0.620824 | 8.36652 | $[-0.915,2.825]\times10^{-5}$ |

Scalar-secant scaling improves some scalar extrema, but it weakens the
negative adjoint lobe. Metric-inverse memory 5 therefore remains the more
balanced Figure 6.2 snapshot candidate. This is a candidate selection, not a
recovery of the source configuration.

With strict tolerances, the original memory-5 metric-inverse policy reached
4,814 iterations and a line-search failure for $\beta=10^{-3}$ and the
5,000-iteration limit for $\beta=10^{-6}$. Increasing memory alone did not
remove the long tail.

A scalar-secant initial inverse-Hessian scaling with memory 20 gave robust
screening behavior: four to five iterations for $\beta=10^{-3}$ and 8 or 13
iterations for $\beta=10^{-6}$ at relative tolerances $10^{-4}$ and $10^{-5}$,
respectively. Driving the small-regularisation case to stricter tolerances
sometimes needed more line-search reductions than the five reported by the
source.

These settings are reasonable for obtaining a stable numerical solution, but
they are not necessarily the source settings. Conversely, a $10^{-3}$
relative-gradient snapshot reaches the Figure 6.2-like early field before the
long-tail behavior matters and accepts unit steps in run `006`. Linear-solver
tolerance variants were not run: the observed forcing and scaling effects are
smooth and systematic, so solve noise is not a plausible explanation for the
remaining discrepancies.

### Promoted development families

Four stable development profiles record the two useful forcing hypotheses and
the two figure-specific solver policies:

- [Figure 6.2, constant half](../../parameters/chapter-6/b1/development/figure-6.2-early-stop-constant-half.prm)
- [Figure 6.2, objective matched](../../parameters/chapter-6/b1/development/figure-6.2-early-stop-objective-matched.prm)
- [Figure 6.3, constant half](../../parameters/chapter-6/b1/development/figure-6.3-constant-half.prm)
- [Figure 6.3, objective matched](../../parameters/chapter-6/b1/development/figure-6.3-objective-matched.prm)

These profiles write stable outputs below `runs/chapter-6/b1/development/`.
The $f=0.5$ profiles are the simple balanced development candidates; the
$f=0.4150674$ profiles preserve the independent initial-objective clue. Neither
family is authoritative or establishes the omitted source forcing.

### Complete source matrix

The two source figures contain seven unique method–regularisation
combinations: steepest descent at $\beta=10^{-1},10^{-2},10^{-3}$ and L-BFGS
at those values plus $\beta=10^{-6}$. A release diagnostic of the unreported
steepest-descent $\beta=10^{-6}$ case reached the 5,000-iteration limit after
approximately six and a half minutes without satisfying the relative-gradient
threshold. The parameter schema now supports exact matrix exclusions, and the
selected authoritative family omits that eighth case instead of retaining it
merely to make the matrix rectangular.

### Authoritative release reproduction

The checked-in
[authoritative parameter profile](../../parameters/chapter-6/b1/authoritative.prm)
freezes the selected $f=0.5$ candidate, regular 131-subdivision triangular
mesh, continuous homogeneous-Dirichlet `P1` control, common $10^{-3}$ relative
gradient threshold, and metric-inverse memory-5 L-BFGS policy. The resulting
release run at framework revision `631537a` recorded parameter hash
`fnv1a64:d9f046cb19ac4d2f` and completed all seven intended combinations:

| Method | $\beta$ | Accepted iterations | Final objective | Final gradient norm | Wall time (s) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Steepest descent | $10^{-1}$ | 64 | 0.0476791331 | $1.55386\times10^{-5}$ | 5.842 |
| Steepest descent | $10^{-2}$ | 547 | 0.0389393934 | $1.56861\times10^{-5}$ | 45.077 |
| Steepest descent | $10^{-3}$ | 2,175 | 0.0137935422 | $1.58044\times10^{-5}$ | 189.922 |
| L-BFGS | $10^{-1}$ | 2 | 0.0476791320 | $3.40480\times10^{-6}$ | 1.591 |
| L-BFGS | $10^{-2}$ | 4 | 0.0389393836 | $2.27613\times10^{-8}$ | 1.296 |
| L-BFGS | $10^{-3}$ | 4 | 0.0137934431 | $7.48796\times10^{-7}$ | 1.202 |
| L-BFGS | $10^{-6}$ | 4 | $2.50521163\times10^{-5}$ | $1.11636\times10^{-5}$ | 1.246 |

Every case stopped at the declared relative-gradient threshold. At the three
shared regularisation values, steepest descent and L-BFGS agree in final
objective to within $10^{-7}$ while retaining the expected `64/547/2175`
versus `2/4/4` work contrast. All artifacts report valid diagnostics and
passing finite-difference Hessian checks, with errors near
$3.50\times10^{-13}$ and symmetry errors no larger than
$2.28\times10^{-13}$.

The realized mesh has 34,322 active cells, 17,424 physical state vertices,
16,900 independent state coordinates, and 16,900 control coordinates, as
required by the selected regular-mesh replacement rather than the unavailable
source mesh. PNG-only post-processing completed for all seven artifacts. The
Figure 6.3 comparison contains the six source history cases, and the
$\beta=10^{-6}$ L-BFGS field retains the source-like central positive region,
four negative lobes, and positive edge regions in the book-convention
adjoint. The sparse comparison leaves the unreported steepest-descent
$\beta=10^{-6}$ position empty.

This run closes the source-sized execution and comparison gate for the frozen
B1 contract. It does not change the evidence classification of $f=0.5$, the
mesh, or the L-BFGS details: they remain explicit project replacements for
information omitted from the source.

## Assessment and remaining evidence gap

- **Mesh connectivity – high confidence.** The counts imply 400 boundary
  vertices, not unique connectivity, and regular triangles have negligible
  effect on extrema. Exact parity needs the source mesh and a mesh-import
  framework extension; another synthetic mesh is not informative.
- **Forcing – high confidence.** Constants around $f=0.4$–$0.5$ give the best
  aggregate match, but no tested value explains all extrema. $f=0.5$ is the
  simplest balanced candidate and $f=0.4150674$ matches the visible initial
  objective. A denser constant sweep would fit raster-reading uncertainty
  rather than recover a source fact.
- **Control discretisation – high confidence.** Continuous
  homogeneous-Dirichlet `P1` is the stated source choice; changing from
  cellwise control does not explain solver histories. It is already a
  parameter-only choice.
- **Figure 6.2 stopping – medium-high confidence.** The completed release
  sweep shows that a $10^{-3}$ relative-gradient threshold consistently stops
  after four unit steps and produces the sign-changing $\beta=10^{-6}$
  adjoint across all tested constants. The source still does not state this
  stopping rule for the field snapshots.
- **Figure 6.3 policy – medium-high confidence.** The recovered
  objective-target L-BFGS policy reproduces the published trend and
  approximate iteration scale for both promoted forcing hypotheses. The
  iteration counts do not identify the forcing; exact validation would require
  source arrays.
- **L-BFGS memory and scaling – high confidence.** Metric-inverse memory 5 is
  retained for the Figure 6.2 snapshots because it better preserves the
  negative adjoint lobe. Scalar-secant memory 20 remains the more robust
  project policy for later convergence. Both are parameter choices, not
  recovered source facts. Exact full BFGS would require separate framework
  work and is lower priority because Figure 6.3 labels LM-BFGS.
- **Quadrature and target evaluation – medium confidence.** These remain
  unresolved but are unlikely to explain the large $\beta=10^{-6}$ morphology
  change at the tested resolution. Alternative assembly policies would
  require a scoped framework change and rebuild.
- **Adjoint sign and plotting – high confidence.** Native $p$ and comparison
  field $`p_{\mathrm{book}}=-p`$ are now explicit; the run `006` sign changes are
  numerical rather than a color-map artifact. No further framework change is
  needed.

The planned parameter-only experiments are complete. No remaining sweep can
uniquely resolve the source omissions: more constants would overfit raster
values, solver tolerances are unlikely to explain the smooth discrepancies,
and additional synthetic meshes cannot recover connectivity. The next
decisive evidence would be the original mesh, numerical arrays, or
implementation details. Supporting arbitrary mesh import or alternative
quadrature and target-evaluation policies would require scoped framework
changes and a rebuild.

The selected authoritative candidate combines $f=0.5$, the regular
131-subdivision triangular mesh, continuous homogeneous-Dirichlet `P1`
control, metric-inverse L-BFGS memory 5, and the common $10^{-3}$ relative
gradient threshold. Keep $f=0.4150674$ as an objective-matched diagnostic. The
project has reproduced the model, field identities, early-iterate morphology,
and method trend; it has not uniquely reproduced the omitted source data. The
completed authoritative release run makes B1 reproduction-verified under the
frozen project contract without claiming book-equivalent source parity.
