# Chapter 6 numerical-examples reference

## Purpose and use

This is the deferred reproduction catalogue for the numerical examples in
Chapter 6 of *Optimal Control of Partial Differential Equations* (Manzoni,
Quarteroni, Salsa; book pages 187–218).

Each record gives source equation, data, discretisation, method, parameters,
and reported comparison quantities.

The source has a few notation inconsistencies, noted rather than silently
normalised. In particular, it uses both $\alpha$ and $\beta$ for a
regularisation coefficient in some captions and tables. Equation and section
number are authoritative for each record.

## Source catalogue policy

This document is authoritative for the source content of all Chapter 6
numerical examples: equations, data, discretisations, algorithms, reported
tables, figure references, and source omissions. Each example must distinguish
what the book states from what can be inferred mathematically and from what
remains unresolved.

This catalogue does not freeze project implementation choices. A forcing
replacement, mesh generator, finite-element realization, quadrature rule,
solver option, boundary interpretation, or stopping policy that is not an
obvious consequence of the source belongs in the owning application,
benchmark, or implementation contract. Those documents may link here for the
source record but should not silently rewrite it.

## Section 6.5 — Reduced-space examples

### E6.5.1 — Distributed Laplace control

**Source:** equations (6.64), Figures 6.2–6.3, book pages 187–188.

~~~math
\begin{aligned}
\min_{y,u}\quad &
\frac{1}{2}\lVert y-z_{d}\rVert_{L^{2}(\Omega)}^{2}
+\frac{\beta}{2}\lVert u\rVert_{L^{2}(\Omega)}^{2}, \\
-\Delta y &= f+u && \text{in }\Omega, \\
y&=0 && \text{on }\partial\Omega.
\end{aligned}
~~~

Take $\Omega=(0,1)^{2}$ and
$z_{d}(x)=10x_{1}(1-x_{1})x_{2}(1-x_{2})$. The supplied text does not
specialise $f$ beyond equation (6.64); a reproduction must resolve it from
the original source material or declare a manufactured replacement.

- **Discretisation:** continuous linear triangular (`P1`) state, adjoint,
  and distributed-control spaces; 17,361 vertices and 34,320 triangles;
  reported $N_{y}=N_{u}=N_{p}=16{,}961$.
- **Methods:** the text names steepest descent and BFGS with backtracking
  Armijo; Figure 6.3 labels the quasi-Newton curves LM-BFGS. The comparison
  uses $\beta=10^{-1}$, $10^{-2}$, and $10^{-3}$; Figure 6.2 also plots BFGS
  fields for $\beta=10^{-3}$ and $10^{-6}$.
- **Start and stopping:** the discussion compares gradient norm with the
  no-control start $u_{0}=0$. Steepest descent uses tolerance $10^{-3}$.
  BFGS is stopped when its objective reaches the steepest-descent final
  objective, rather than at a separate common tolerance.
- **Armijo data:** $\rho=0.7$, $\sigma=10^{-5}$, and at most five step
  rescalings; minimum step $0.01$ for BFGS and $0.2$ for steepest descent.
- **Published trend:** smaller $\beta$ lowers objective and tracking error;
  steepest descent is much slower than L-BFGS.

Source references:
- [page 187, equations (6.64) and test-case data](assets/chapter-6/source-page-187.png);
- [page 188, Figures 6.2–6.3](assets/chapter-6/source-page-188.png).

**Source completeness.** The source does not provide the forcing $f$, the
coordinates or connectivity of its triangular mesh, the quadrature and
target-data evaluation rules, the linear-solver tolerances, or the exact
BFGS/LM-BFGS memory and update policy. It also publishes Figures 6.2–6.3 as
rendered plots rather than numerical field or iteration-history data.
Consequently, the source alone cannot determine the exact coefficient vectors,
objective histories, or field images.

The source does state that Figure 6.2 includes BFGS fields for
$\beta=10^{-3}$ and $10^{-6}$, while Figure 6.3 contains cost and gradient
histories for $\beta=10^{-1},10^{-2},10^{-3}$. The figures do not provide the
underlying numerical arrays. Any project policy used to fill the missing
forcing, mesh, discretisation, or solver details cannot be recovered from the
published figures alone.

### E6.5.2 — Graetz-flow boundary control

**Source:** equation (6.65), Table 6.2, Figures 6.4–6.5, book pages 188–190.

~~~math
\begin{aligned}
\min_{y,u}\quad &
\frac{1}{2}\lVert y-z_{d}\rVert_{L^{2}(\Omega_{0})}^{2}
+\frac{\beta}{2}\lVert u\rVert_{L^{2}(\Gamma_{c})}^{2}, \\
-\mathrm{div}(\mu\nabla y)+\mathrm{div}(b y)&=0 && \text{in }\Omega, \\
y&=1 && \text{on }\Gamma_{D}, \\
\partial_{n}y-(b\mathbin\cdot n)y&=u && \text{on }\Gamma_{c}, \\
\partial_{n}y-(b\mathbin\cdot n)y&=0 && \text{on }\Gamma_{\mathrm{out}}.
\end{aligned}
~~~

- **Data:** $\mu=0.1$, $\beta=10^{-3}$,
  $\Omega=(0,1+l)\mathbin\times(0,1)$ with $l=3$, and
  $b(x)=(1.5x_{2}(1-x_{2}),0)$. The state space is
  $H^{1}_{\Gamma_{D}}(\Omega)$ and the control space is
  $L^{2}(\Gamma_{c})$.
- **Scenarios:**
  $\Omega_{0}^{1}=\{x:x_{1}>1,\ x_{2}<0.3\ \text{or}\ x_{2}>0.7\}$ and
  $\Omega_{0}^{2}=\{x:x_{1}>1\}$; target choices are $z_{d}^{1}=2$ and
  $z_{d}^{2}=4x_{2}(1-x_{2})$.
- **Boundary geometry (Figure 6.4):** $\Gamma_{D}$ contains the left edge
  and the upstream portions of the top and bottom walls for $0\leq x_{1}\leq 1$;
  $\Gamma_{c}$ contains the downstream wall portions for $1\leq x_{1}\leq 4$;
  and $\Gamma_{\mathrm{out}}$ is the right edge $x_{1}=4$.
- **Discretisation:** `P1` triangular finite elements; 11,028 vertices and
  21,653 triangles; reported $N_{y}=N_{p}=10{,}907$ and $N_{u}=243$.
- **Method:** BFGS for all scenarios. The nearby text attributes high
  iteration counts to a constant step but does not state its value; a
  reproduction must declare a fixed-step or line-search policy.

| Case | Observation | Target | Iterations | $J_{h}(u_{0})$ | $J_{h}(\hat u)$ | Reduction | Relative gradient |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| a | $\Omega_{0}^{1}$ | $z_{d}^{1}$ | 59 | 316.6661 | 3.5682 | 98.87% | 0.0250 |
| b | $\Omega_{0}^{2}$ | $z_{d}^{1}$ | 54 | 192.8385 | 2.6368 | 98.63% | 0.0569 |
| c | $\Omega_{0}^{1}$ | $z_{d}^{2}$ | 48 | 29.2188 | 0.7826 | 97.32% | 0.0753 |
| d | $\Omega_{0}^{2}$ | $z_{d}^{2}$ | 87 | 45.9996 | 0.8464 | 98.16% | 0.0387 |

Source references:
- [page 189, equations (6.65), Figure 6.4 and Table 6.2](assets/chapter-6/source-page-189.png);
- [page 190, Figure 6.5](assets/chapter-6/source-page-190.png).

**Source completeness.** The book does not provide the mesh connectivity,
boundary-node subdivision, or the exact basis and quadrature used for the
$N_{u}=243$ control degrees of freedom. It also does not fully disambiguate
whether the boundary notation $\partial_{n}y-(b\mathbin\cdot n)y$ uses an
ordinary or diffusion-weighted normal derivative. Exact numerical parity
therefore requires the source mesh and implementation. The ordinary-normal
versus diffusion-weighted-conormal choice is unresolved in this source record
and cannot be resolved from the published notation alone.

The source does not state the value of the constant BFGS step mentioned in the
discussion, the stopping tolerance that produced the reported iteration
counts, the linear-solver tolerances, or the precise target and boundary
quadrature policies. Figure 6.5 contains rendered state fields but no source
field data or plotting configuration. These details remain unresolved in the
source record.

### E6.5.3 — Reduced-space Stokes control

**Source:** equation (6.66), Table 6.3, Figure 6.6, book pages 190–192.

~~~math
\begin{aligned}
\min_{v,\pi,u}\quad &
\frac{1}{2}\lVert v-v_{d}\rVert_{L^{2}(\Omega)^{2}}^{2}
+\frac{\alpha}{2}\lVert u\rVert_{L^{2}(\Omega)^{2}}^{2}, \\
-\nu\Delta v+\nabla\pi&=f+u && \text{in }\Omega, \\
\mathrm{div} v&=0 && \text{in }\Omega.
\end{aligned}
~~~

The velocity is zero on $\partial\Omega\setminus\Gamma_{\mathrm{right}}$. On
$\Gamma_{\mathrm{right}}$, the source uses $v\mathbin\cdot t=0$ and
$(T(v,\pi)n)\mathbin\cdot n=0$, with
$T(v,\pi)n=-\pi n+\nu\partial_{n}v$.

- **Data:** $\Omega=(0,1)\mathbin\times(0,2)$,
  $\Gamma_{\mathrm{right}}=\{x:x_{1}=1\}$, $f=-e_{2}$, and
  $v_{d}=x_{2}e_{1}$.
- **Discretisation:** inf-sup compatible `P2–P1` pair on 5,140 vertices and
  10,278 triangles.
- **Method:** steepest descent with fixed $\tau=5$ for
  $\alpha=10^{-2}$, $10^{-3}$, and $10^{-4}$.
- **Naming:** equation (6.66) uses $\alpha$, while the accompanying Table
  6.3 labels its coefficient column $\beta$.

| $\beta$ | $J_{h}(\hat u)$ | Iterations | $\lVert\hat u\rVert_{L^{2}(\Omega)}$ |
| ---: | ---: | ---: | ---: |
| $10^{-2}$ | 0.00317078 | 29 | 0.752571 |
| $10^{-3}$ | 0.000383566 | 70 | 0.827743 |
| $10^{-4}$ | $5.9189\mathbin\cdot10^{-5}$ | 90 | 0.847255 |

The source gives $J_{h}(u_{0})=0.076222$ and reports that smaller
regularisation strengthens control action and lowers cost.

Source references:
- [page 190, equations (6.66) and boundary conditions](assets/chapter-6/source-page-190.png);
- [page 191, Figure 6.6](assets/chapter-6/source-page-191.png);
- [page 192, Table 6.3](assets/chapter-6/source-page-192.png).

**Source completeness.** The source does not state the viscosity $\nu$, a
pressure gauge, the mesh connectivity, the linear-solver tolerances, or the
stopping criterion behind the reported convergence counts. The prose says
that the distance $\lVert\hat v-v_{d}\rVert_{L^{2}(\Omega)^{2}}$ is reported,
but the printed Table 6.3 contains no distance column. Figure 6.6 is a
rendered field figure without source field arrays or plotting settings. These
details remain unresolved in the source record.

## Section 6.7 — All-at-once scalar examples

### E6.7.1 — All-at-once Laplace control

**Source:** equations (6.84)–(6.88), Table 6.4, Figure 6.7, book pages
202–205.

~~~math
\begin{aligned}
\min_{y,u}\quad &
\frac{1}{2}\lVert y-z_{d}\rVert_{L^{2}(\Omega)}^{2}
+\beta\lVert u\rVert_{L^{2}(\Omega)}^{2}, \\
-\Delta y&=u && \text{in }\Omega, \\
y&=z_{d} && \text{on }\partial\Omega.
\end{aligned}
~~~

The source uses $\beta\lVert u\rVert^{2}$, not
$\beta\lVert u\rVert^{2}/2$; its KKT control block is consequently
$2\beta M$. It takes $\Omega=(0,1)^{2}$ and

~~~math
z_{d}(x)=
\begin{cases}
(2x_{1}-1)^{2}(2x_{2}-1)^{2} & x\in(0,1/2]^{2}, \\
0 & \text{otherwise}.
\end{cases}
~~~

The source uses a DtO discretisation with stiffness matrix $K$, mass matrix
$M$, $b=Mz_{d}$, and a boundary-value vector $d$. Its discrete state
equation is $Ky=Mu+d$, and the KKT system is

~~~math
\begin{bmatrix}
2\beta M & 0 & M^{\mathsf T} \\
0 & M & K^{\mathsf T} \\
M & K & 0
\end{bmatrix}
\begin{bmatrix}u \\ y \\ p\end{bmatrix}
=
\begin{bmatrix}0 \\ b \\ d\end{bmatrix}.
~~~

- **Discretisation/formulation:** DtO with conforming bilinear `Q1`
  elements; Table 6.4 uses $h=2^{-2},\ldots,2^{-9}$.
- **Methods:** direct Matlab backslash, MINRES with geometric-MG or AMG
  block-diagonal preconditioners, and PPCG with a constraint preconditioner.
- **Preconditioner settings:** two multigrid V-cycles; two pre- and
  post-smoothing relaxed-Jacobi steps with $\omega=8/9$; Chebyshev or cheap
  lumped-style mass approximations. The block-diagonal Schur approximation
  uses the approximate stiffness and mass solves; the constraint
  preconditioner uses the corresponding approximate $K$ and $M$ actions and
  is also tested with $G=\mathrm{diag}(K)$.
- **Sweep/stopping:** $\beta=10^{-2}$, $5\mathbin\cdot10^{-5}$, $10^{-5}$,
  and $10^{-8}$; tolerance $10^{-8}$. MINRES uses relative Euclidean
  residual; PPCG uses $r^{\mathsf T}g$ relative to its initial value.
- **Source provenance:** the example is adapted from [226] and the book
  footnote points to the authors' [Poisson-control implementation](https://github.com/tyronerees/poisson-control).

Table 6.4 reports timings and outer iterations for the three tabulated
coefficients. The compact projection below retains the iteration counts; each
entry is geometric-MG MINRES / AMG MINRES / PPCG / PPCG with
$G=\mathrm{diag}(K)$. The complete source table, including timings and
$3N_{y}$, is linked below.

| $h$ | $10^{-2}$ | $5\mathbin\cdot10^{-5}$ | $10^{-5}$ |
| --- | --- | --- | --- |
| $2^{-2}$ | 12 / 12 / 3 / 6 | 14 / 20 / 6 / 6 | 18 / 22 / 8 / 6 |
| $2^{-3}$ | 12 / 12 / 3 / 11 | 28 / 32 / 10 / 11 | 44 / 50 / 15 / 10 |
| $2^{-4}$ | 12 / 12 / 3 / 11 | 30 / 34 / 11 / 11 | 48 / 58 / 19 / 11 |
| $2^{-5}$ | 12 / 12 / 3 / 11 | 30 / 36 / 11 / 11 | 51 / 60 / 19 / 10 |
| $2^{-6}$ | 12 / 12 / 3 / 10 | 30 / 36 / 10 / 10 | 51 / 60 / 18 / 10 |
| $2^{-7}$ | 12 / 13 / 3 / 10 | 32 / 36 / 10 / 10 | 52 / 60 / 16 / 10 |
| $2^{-8}$ | 12 / 14 / 3 / 10 | 32 / 38 / 8 / 10 | 54 / 60 / 14 / 10 |
| $2^{-9}$ | 12 / 16 / 3 / 9 | 32 / 38 / 8 / 9 | 54 / 62 / 14 / 9 |

The source reports mesh-independent outer iterations at moderate $\beta$ and
larger counts as $\beta$ decreases. Its runtime table is not portable.

Source references:
- [page 202, problem data and equations (6.84)](assets/chapter-6/source-page-202.png);
- [page 203, discrete KKT system and preconditioners](assets/chapter-6/source-page-203.png);
- [page 204, Figure 6.7 and iteration trends](assets/chapter-6/source-page-204.png);
- [page 205, complete Table 6.4](assets/chapter-6/source-page-205.png).

**Source completeness.** The narrative describes mesh levels through
$h=2^{-8}$, while Table 6.4 includes $h=2^{-9}$. The table’s timing values
are tied to the source’s Matlab/HSL environment and a footnoted 2.2 GHz Intel
Core i7 with 16 GB RAM; the source does not provide the exact mesh
connectivity, software versions, or timing-unit definition. Figure 6.7 labels
the second method family as PCG, while the surrounding text calls it PPCG.

### E6.7.2 — All-at-once diffusion-reaction control

**Source:** equations (6.89)–(6.91), Figure 6.8, book pages 204–206.

This has the same form and data as E6.7.1, including its target, domain,
`Q1` DtO discretisation, preconditioners, parameter sweep, and stopping
criteria. The state equation is instead

~~~math
-\Delta y+y=u \quad\text{in }\Omega,
\qquad y=z_{d}\quad\text{on }\partial\Omega.
~~~

With $b=Mz_{d}$ and the same boundary-value vector $d$, the corresponding
discrete KKT system is

~~~math
\begin{bmatrix}
2\beta M & 0 & M^{\mathsf T} \\
0 & M & (K+M)^{\mathsf T} \\
M & K+M & 0
\end{bmatrix}
\begin{bmatrix}u \\ y \\ p\end{bmatrix}
=
\begin{bmatrix}0 \\ b \\ d\end{bmatrix}.
~~~

The block-diagonal preconditioners have the same form as in E6.7.1, with
approximations of the stiffness solves replaced by approximations of $K+M$;
the source states that the same multigrid solvers are used. Figure 6.8 plots
the MINRES method with $\widehat P_{da}$ and PPCG with $\widehat P_{ca}$ for
$\beta=10^{-2},5\mathbin\cdot10^{-5},10^{-5},10^{-8}$ as functions of the
mesh size $h$. The source reports trends similar to the Laplace case and
does not repeat the timing table: it states that the computational-time
dependence on $h$ and $\beta$ is the same as in Table 6.4.

Source references:
- [page 204, start of the diffusion-reaction case and equation (6.89)](assets/chapter-6/source-page-204.png);
- [page 205, equation (6.91) and comparison with case 1](assets/chapter-6/source-page-205.png);
- [page 206, Figure 6.8](assets/chapter-6/source-page-206.png).

**Source completeness.** The source does not provide the numerical arrays
underlying Figure 6.8, its plot-generation settings, or a separate timing
table. The Figure 6.8 caption labels the right-hand method as PCG, while the
surrounding text and the preconditioner notation call it PPCG.

## Section 6.9 — PDAS examples

### E6.9.1 — Symmetric box-constrained Laplace control

**Source:** equation (6.105), Figures 6.9–6.13, book pages 211–214.

~~~math
\begin{aligned}
\min_{y,u}\quad &
\frac{1}{2}\lVert y-z_{d}\rVert_{L^{2}(\Omega)}^{2}
+\frac{\beta}{2}\lVert u\rVert_{L^{2}(\Omega)}^{2}, \\
-\Delta y&=u && \text{in }\Omega, \\
y&=z_{d} && \text{on }\partial\Omega, \\
-a&\leq u\leq a.
\end{aligned}
~~~

Take $\Omega=(0,1)^{2}$, $a>0$, and
$z_{d}=\sin(2\pi x_{1})\sin(2\pi x_{2})$. The source quotes the
unconstrained exact control
$u_{\mathrm{ex}}=8\pi^{2}\sin(2\pi x_{1})\sin(2\pi x_{2})$ and considers
$\beta>0$.

- **Discretisation:** conforming bilinear `Q1` elements on a quadrilateral
  mesh; mesh size is not stated.
- **Method:** PDAS from $u_{0}=0$ and multiplier $\lambda_{0}=0$.
- **Sweep:** $\beta=10^{-2},10^{-4},10^{-6},10^{-8}$ and
  $a=100,70,50,30$.
- **Published trend:** $a=100$ converges in one iteration because
  $\lVert u_{\mathrm{ex}}\rVert_{L^{\infty}(\Omega)}\leq80$. Tighter bounds
  require more PDAS iterations, and finer meshes increase the iteration count
  for fixed $a$ and $\beta$. Smaller $\beta$ increases the control norm and
  active-region extent while lowering tracking error; larger $\beta$ lowers
  the iteration count and control norm while increasing the distance from the
  target.

Figures 6.9–6.13 show, respectively, the optimal state for $\beta=10^{-8}$
and $a=30,50,70$ plus the unconstrained case; the unconstrained optimal
control for the four $\beta$ values; constrained controls for the three
finite bounds and four $\beta$ values; the iteration count, target distance,
and control norm as functions of $\beta$ for $a=100,70,50,30$; and the
iteration count versus mesh size for $\beta=10^{-6},10^{-8}$ and
$a=30,50,70$.

Source references:
- [page 211, equations (6.105) and test-case data](assets/chapter-6/source-page-211.png);
- [page 212, Figures 6.9–6.10](assets/chapter-6/source-page-212.png);
- [page 213, Figures 6.11–6.12](assets/chapter-6/source-page-213.png);
- [page 214, Figure 6.13](assets/chapter-6/source-page-214.png).

**Source completeness.** The source does not state the mesh size or
connectivity, quadrature and target-evaluation rules, the PDAS switching
parameter, the stopping tolerance or maximum iteration count, the linear
solver and its tolerances, or the numerical arrays and plotting settings
behind Figures 6.9–6.13. One paragraph refers to the regularisation
coefficient as $\alpha$, although equation (6.105), the sweep, and the figure
captions use $\beta$.

### E6.9.2 — Asymmetric box-constrained Laplace control

**Source:** equation (6.106), Figures 6.14–6.15, book pages 214–216.

~~~math
\begin{aligned}
\min_{y,u}\quad &
\frac{1}{2}\lVert y-z_{d}\rVert_{L^{2}(\Omega)}^{2}
+\frac{\beta}{2}\lVert u\rVert_{L^{2}(\Omega)}^{2}, \\
-\Delta y&=u && \text{in }\Omega=(0,1)^{d}, \\
y&=z_{d} && \text{on }\partial\Omega, \\
a(x)&\leq u(x)\leq b(x).
\end{aligned}
~~~

For $d=2$ and $d=3$, use the target and bounds

~~~math
\begin{aligned}
z_{d}(x)&=-x_{1}\exp\left(-\sum_{i=1}^{d}(x_{i}-0.5)^{2}\right), \\
a(x)&=
\begin{cases}
-0.15 & x_{1}<0.5, \\
-0.2 & x_{1}\geq0.5,
\end{cases} \\
b(x)&=-0.01\exp\left(-\sum_{i=1}^{d}x_{i}^{2}\right).
\end{aligned}
~~~

- **Discretisation:** the source states conforming piecewise bilinear `Q1`
  finite elements on a quadrilateral mesh for $d=2$ and a hexahedral mesh for
  $d=3$, respectively.
- **Method/sweep:** PDAS from $u_{0}=0$, $\lambda_{0}=0$ for
  $\beta=10^{-2},10^{-4},10^{-6}$.
- **Published provenance:** the source identifies this test case as
  originally proposed in [228].
- **Published trend:** iteration count, target distance, and control norm vary
  with $\beta$ as in E6.9.1. With the discontinuous lower bound $a(x)$, a
  large regularisation coefficient can limit the control norm, as illustrated
  by the top row of Figure 6.14.

Figures 6.14 and 6.15 show the optimal state, optimal control, and adjoint
state for dimensions $d=2$ and $d=3$, respectively, with
$\beta=10^{-2},10^{-4},10^{-6}$ from top to bottom.

Source references:
- [page 214, equations (6.106), target, and bounds](assets/chapter-6/source-page-214.png);
- [page 215, Figure 6.14](assets/chapter-6/source-page-215.png);
- [page 216, Figure 6.15](assets/chapter-6/source-page-216.png).

**Source completeness.** The source does not state the mesh sizes or
connectivity, quadrature and target-evaluation rules, the PDAS switching
parameter, the stopping tolerance or maximum iteration count, the linear
solver and its tolerances, or the numerical arrays and plotting settings
behind Figures 6.14–6.15. It reports qualitative trends rather than a table
of numerical values.

### E6.9.3 — Box-constrained Stokes control

**Source:** equation (6.107), Table 6.5, Figure 6.16, book pages 215–218.

~~~math
\begin{aligned}
\min_{v,\pi,u}\quad &
\frac{1}{2}\lVert v-v_{d}\rVert_{L^{2}(\Omega)^{2}}^{2}
+\frac{\alpha}{2}\lVert u\rVert_{L^{2}(\Omega)^{2}}^{2}, \\
-\nu\Delta v+\nabla\pi&=f+u && \text{in }\Omega, \\
\mathrm{div} v&=0 && \text{in }\Omega, \\
v&=g && \text{on }\partial\Omega, \\
-550&\leq u_{i}\leq550 && i=1,2.
\end{aligned}
~~~

The domain is $\Omega=(0,1)^{2}$ and
$\Gamma_{\mathrm{right}}=\{x:x_{1}=1\}$. The boundary datum is $g=0$ on
$\Gamma_{\mathrm{right}}$ and $g=-e_{2}$ elsewhere. The target is

~~~math
v_{d}=(x_{2}-0.5)e_{1}-(x_{1}-0.5)e_{2}.
~~~

- The source sets
  $V=H_{0}^{1}(\Gamma_{\mathrm{right}})^{2}\times L^{2}(\Omega)$ and
  $U_{\mathrm{ad}}=\{u\in L^{2}(\Omega)^{2}:-550\leq u_{i}\leq550
  \text{ a.e. in }\Omega,\ i=1,2\}$.
- **Discretisation:** the inf-sup compatible `Q2–Q1` pair; the reported
  dimensions are $N_{v}=51{,}842$ and $N_{\pi}=6{,}561$.
- **Method:** the numerical-examples section presents the case as a PDAS
  computation, but this case does not state the initial control or multiplier.
- **Naming:** equation (6.107) calls regularisation $\alpha$, while Figure
  6.16 and Table 6.5 label the same sweep $\beta$.
- **Published provenance:** the source identifies this test case as
  originally proposed in [229].

| $\beta$ | $J_{h}(\hat u)$ | Iterations | $\lVert\hat u\rVert_{L^{2}(\Omega)}$ | $\lVert\hat v-v_{d}\rVert_{L^{2}(\Omega)}$ |
| ---: | ---: | ---: | ---: | ---: |
| $10^{-2}$ | 0.0579601 | 1 | 0.326406 | 0.338902 |
| $10^{-4}$ | 0.0430436 | 1 | 10.3743 | 0.274453 |
| $10^{-6}$ | 0.0244265 | 1 | 72.2504 | 0.208885 |
| $10^{-8}$ | 0.0176069 | 6 | 308.912 | 0.185093 |

The uncontrolled reference is $J_{h}(u_{0})=0.0585097$. The target velocity
is incompatible with the prescribed state-velocity boundary data, so the
cost cannot be reduced below a certain value. The source observes that the
control moves the vortex from the right side towards the centre, with a more
pronounced shift for smaller $\alpha$. Table 6.5 reports that the iteration
count, target distance, and control norm vary with the regularisation as in
E6.9.1; its caption calls the iterations those of steepest descent, despite
the surrounding numerical-examples section presenting the tests as PDAS.

Figure 6.16 shows the optimal velocity and pressure, adjoint velocity, and
optimal control for $\beta=10^{-2},10^{-4},10^{-6},10^{-8}$, ordered from top
to bottom and left to right. The source notes that pressure peaks are
localized near the boundary and are not clearly visible.

Source references:
- [page 215, equation (6.107), domain, and boundary data](assets/chapter-6/source-page-215.png);
- [page 216, target, state space, and box constraint](assets/chapter-6/source-page-216.png);
- [page 217, discretisation, trends, and Figure 6.16](assets/chapter-6/source-page-217.png);
- [page 218, complete Table 6.5](assets/chapter-6/source-page-218.png).

**Source completeness.** The source does not state the forcing $f$, viscosity
$\nu$, pressure gauge, mesh size or connectivity, quadrature and target-data
evaluation rules, the PDAS initial guesses or switching/stopping policy, the
linear solver and its tolerances, or the numerical arrays and plotting
settings behind Figure 6.16. The Table 6.5 caption’s steepest-descent label
conflicts with the section’s PDAS context, and the source does not resolve
whether its iteration column records PDAS or steepest-descent iterations.
