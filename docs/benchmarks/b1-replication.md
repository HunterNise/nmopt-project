# B1 replication findings

## Scope and status

This report records the development investigation of E6.5.1, the distributed
Laplace-control example shown in Figures 6.2–6.3 of the source. It supplements
the [frozen B1 benchmark contract](chapter-6.md#b1--e651-distributed-laplace-control)
and the [source catalogue](../guides/chapter-6-numerical-examples.md#e651--distributed-laplace-control).
It does not change the frozen benchmark or promote an inferred source detail
to an authoritative choice.

The evidence was produced with the Debug deal.II build through framework
revision `18e005f`. The tracked authoritative profile remains a
framework-native manufactured benchmark; it is stale with respect to this
investigation and must not be described as a book-equivalent reproduction.
Development outputs under `runs/` are diagnostic evidence and may be deleted.

The current assessment is:

- the mathematical problem, field identities, sign convention, and the
  qualitative Figure 6.3 method comparison are reproduced;
- continuous `P1` control, constant-forcing, triangular-mesh, and solver-policy
  candidates can now be selected independently in parameter files;
- neither the source mesh nor the omitted forcing has been uniquely recovered;
- the Figure 6.2 extrema cannot all be explained by one tested constant
  forcing at the later tolerance-terminated solutions; and
- `runs/chapter-6/b1/development/006` currently gives the closest
  Figure 6.2 negative-adjoint morphology for $\beta=10^{-6}$, but it is an
  intentionally truncated four-iteration solution rather than a converged
  reference.

The benchmark is therefore **framework-verified but not
reproduction-verified**.

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

Development run `008` implements the recovered Figure 6.3 policy with zero
forcing and the framework-native mesh. Steepest descent stops at relative
gradient norm $10^{-3}$; L-BFGS stops when it reaches the corresponding final
steepest-descent objective.

| $\beta$ | Steepest-descent iterations | L-BFGS iterations | Matched objective |
| ---: | ---: | ---: | ---: |
| $10^{-1}$ | 64 | 2 | 0.0541866 |
| $10^{-2}$ | 548 | 4 | 0.0443562 |
| $10^{-3}$ | 2,202 | 4 | 0.0158240 |

The histories reproduce the source trend and the approximate locations and
shapes visible in the rasterized plots: steepest descent becomes dramatically
slower as $\beta$ decreases, while L-BFGS reaches the reference objectives in
a few iterations. The source arrays are unavailable, so this is qualitative
and raster-level agreement rather than numerical parity.

Changing only the control from cellwise constant to continuous homogeneous
Dirichlet `P1` in development run `009` preserves the iteration counts
`64/548/2202` and `2/4/4`. The solver trend therefore does not identify the
control representation. The source explicitly states continuous `P1`, so that
is the appropriate source-replication choice even though the frozen
framework-native benchmark uses cellwise control.

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

No single tested constant forcing explains both source rows. In particular,
the extrema do not justify promoting $f=1$ from a motivated hypothesis to a
source fact. Zero forcing remains plausible for Figure 6.3 and the
$\beta=10^{-3}$ state, while a positive source improves other extrema.

### Four-iteration $\beta=10^{-6}$ field

Development run `006` uses zero forcing, cellwise control, and exactly four
accepted L-BFGS iterations. Its $\beta=10^{-6}$ field has:

- state maximum `0.618087`;
- control maximum `9.59824`;
- $`p_{\mathrm{book}}`$ range
  $[-1.68607\times10^{-5},4.32398\times10^{-5}]$; and
- final relative gradient
  $1.27534\times10^{-5}/0.0168617=7.56\times10^{-4}$.

Its book-convention adjoint has the same central positive region, four negative
lobes, and positive edge regions as the source plot. Its range is also much
closer to the source range than any later tolerance-terminated screen. The
similarity is not caused by plotting: the sign changes are present in the
native VTU field.

This makes a relative-gradient stopping threshold near $10^{-3}$ the leading
candidate for the Figure 6.2 BFGS snapshots. Such a policy stops after four
iterations in run `006` and is independently stated by the source for steepest
descent. The source does not state that it was used for these BFGS fields, so
the interpretation remains a hypothesis.

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

The omitted L-BFGS scaling and stopping choices materially affect the
$\beta=10^{-6}$ solve. With strict tolerances, the original memory-5
metric-inverse policy reached 4,814 iterations and a line-search failure for
$\beta=10^{-3}$ and the 5,000-iteration limit for $\beta=10^{-6}$. Increasing
memory alone did not remove the long tail.

A scalar-secant initial inverse-Hessian scaling with memory 20 gave robust
screening behavior: four to five iterations for $\beta=10^{-3}$ and 8 or 13
iterations for $\beta=10^{-6}$ at relative tolerances $10^{-4}$ and $10^{-5}$,
respectively. Driving the small-regularisation case to stricter tolerances
sometimes needed more line-search reductions than the five reported by the
source.

These settings are reasonable for obtaining a stable numerical solution, but
they are not necessarily the source settings. Conversely, a $10^{-3}$
relative-gradient snapshot reaches the Figure 6.2-like early field before the
long-tail behavior matters and accepts unit steps in run `006`.

## Inferences and remaining experiments

- **Mesh connectivity – high confidence.** The counts imply 400 boundary
  vertices, not unique connectivity, and regular triangles have negligible
  effect on extrema. Exact parity needs the source mesh and a mesh-import
  framework extension; another synthetic mesh is not informative.
- **Forcing – high confidence.** Constant $f=0$ and $f=1$ each explain
  different observations, but no tested constant explains all fields.
  Further constants or expressions are parameter-only changes and should not
  be frozen without stronger evidence.
- **Control discretisation – high confidence.** Continuous
  homogeneous-Dirichlet `P1` is the stated source choice; changing from
  cellwise control does not explain solver histories. It is already a
  parameter-only choice.
- **Figure 6.2 stopping – medium confidence.** A relative-gradient threshold
  near $10^{-3}$ plausibly explains the four-iteration, sign-changing
  $\beta=10^{-6}$ adjoint. Testing it requires parameter changes only.
- **Figure 6.3 policy – medium-high confidence.** The recovered
  objective-target L-BFGS policy reproduces the published trend and
  approximate iteration scale. No framework change is needed; exact
  validation would require source arrays.
- **L-BFGS memory and scaling – high confidence.** Scalar-secant memory 20 is
  a robust project policy, not a recovered source fact. It is parameter-only.
  Exact full BFGS would require separate framework work and is lower priority
  because Figure 6.3 labels LM-BFGS.
- **Quadrature and target evaluation – medium confidence.** These remain
  unresolved but are unlikely to explain the large $\beta=10^{-6}$ morphology
  change at the tested resolution. Alternative assembly policies would
  require a scoped framework change and rebuild.
- **Adjoint sign and plotting – high confidence.** Native $p$ and comparison
  field $`p_{\mathrm{book}}=-p`$ are now explicit; the run-006 sign changes are
  numerical rather than a color-map artifact. No further framework change is
  needed.

The highest-information follow-up is a small Figure 6.2 matrix using continuous
`P1` control, the regular 131-subdivision triangular mesh, zero initial
control, the book Armijo constants, and a $10^{-3}$ relative-gradient stop.
It should compare $f=0$, $0.5$, and $1$ for $\beta=10^{-3}$ and $10^{-6}$.
This combines the source's stated discretisation with the run-006 stopping
clue and requires parameter files only.

Until that experiment or original-source material resolves the contradiction,
the correct conclusion is not that one candidate reproduces B1. The project
has reproduced the model and method trend, bracketed the field extrema, ruled
out regular triangular connectivity as the main cause, and isolated the
$\beta=10^{-6}$ stopping state as the most consequential omitted numerical
detail.
