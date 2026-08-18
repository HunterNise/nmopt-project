# Chapter 6 numerical-examples reference

## Purpose and use

This is the deferred reproduction catalogue for the numerical examples in
Chapter 6 of *Optimal Control of Partial Differential Equations* (Manzoni,
Quarteroni, Salsa; book pages 187–217). It is separate from the
[numerical-methods implementation guide](chapter-6-numerical-methods.md)
so the framework can first be implemented and verified with manufactured
tests.

Each record gives source equation, data, discretisation, method, parameters,
and reported comparison quantities. A reproduction must add its compiler
manifest, mesh generator/version, quadrature and target-data policy, linear
solver/preconditioner tolerances, hardware, and wall-clock definition. Times
and iteration counts are not cross-platform correctness criteria.

Before serialization or comparative-run orchestration exists, the in-memory
[`ReducedExperimentEnvelopeT`](../../include/nmopt/experiment/reduced_envelope.hpp)
binds one copied `CompilationManifest`, one typed solver-policy snapshot, one
typed reduced report, and one caller-supplied `RunEnvironmentRecord`. The
envelope owns values rather than a compiled problem or reduced DTO service, so
the provenance remains available after the executable service is detached.
Timing collection and artifact serialization remain outside this contract.

The source has a few notation inconsistencies, noted rather than silently
normalised. In particular, it uses both $\alpha$ and $\beta$ for a
regularisation coefficient in some captions and tables. Equation and section
number are authoritative for each record.

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
- **Methods:** steepest descent and BFGS/limited-memory BFGS, both with
  backtracking Armijo. The comparison uses $\beta=10^{-1}$, $10^{-2}$, and
  $10^{-3}$; Figure 6.2 also plots BFGS fields for $\beta=10^{-3}$ and
  $10^{-6}$.
- **Start and stopping:** the discussion compares gradient norm with the
  no-control start $u_{0}=0$. Steepest descent uses tolerance $10^{-3}$.
  BFGS is stopped when its objective reaches the steepest-descent final
  objective, rather than at a separate common tolerance.
- **Armijo data:** $\rho=0.7$, $\sigma=10^{-5}$, and at most five step
  rescalings; minimum step $0.01$ for BFGS and $0.2$ for steepest descent.
- **Published trend:** smaller $\beta$ lowers objective and tracking error;
  steepest descent is much slower than L-BFGS.

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

The book does not provide the mesh connectivity, boundary-node subdivision,
or the exact basis and quadrature used for the $N_{u}=243$ control degrees of
freedom. It also does not fully disambiguate whether the boundary notation
$\partial_{n}y-(b\mathbin\cdot n)y$ uses an ordinary or diffusion-weighted
normal derivative. Exact numerical parity therefore requires the source mesh
and implementation, or an explicit framework-native realization policy. The
framework-native B2 policy now selects the ordinary normal derivative
interpretation, while retaining the diffusion-weighted conormal as an explicit
diagnostic alternative.

| Case | Observation | Target | Iterations | $J_{h}(u_{0})$ | $J_{h}(\hat u)$ | Reduction | Relative gradient |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| a | $\Omega_{0}^{1}$ | $z_{d}^{1}$ | 59 | 316.6661 | 3.5682 | 98.87% | 0.0250 |
| b | $\Omega_{0}^{2}$ | $z_{d}^{1}$ | 54 | 192.8385 | 2.6368 | 98.63% | 0.0569 |
| c | $\Omega_{0}^{1}$ | $z_{d}^{2}$ | 48 | 29.2188 | 0.7826 | 97.32% | 0.0753 |
| d | $\Omega_{0}^{2}$ | $z_{d}^{2}$ | 87 | 45.9996 | 0.8464 | 98.16% | 0.0387 |

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
  10,278 triangles. The source does not state viscosity or pressure gauge;
  both must be selected before a reproduction.
- **Method:** steepest descent with fixed $\tau=5$ for
  $\alpha=10^{-2}$, $10^{-3}$, and $10^{-4}$.
- **Naming:** equation (6.66) uses $\alpha$, while the accompanying Table
  6.3 labels its coefficient column $\beta$.

| Regularisation | $J_{h}(\hat u)$ | Iterations | $\lVert\hat u\rVert_{L^{2}(\Omega)}$ |
| ---: | ---: | ---: | ---: |
| $10^{-2}$ | 0.00317078 | 29 | 0.752571 |
| $10^{-3}$ | 0.000383566 | 70 | 0.827743 |
| $10^{-4}$ | $5.9189\mathbin\cdot10^{-5}$ | 90 | 0.847255 |

The source gives $J_{h}(u_{0})=0.076222$ and reports that smaller
regularisation strengthens control action and lowers cost.

## Section 6.7 — All-at-once scalar examples

### E6.7.1 — All-at-once Laplace control

**Source:** equations (6.84)–(6.88), Table 6.4, Figure 6.7, book pages
202–204.

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

- **Discretisation/formulation:** DtO with conforming bilinear `Q1`
  elements; Table 6.4 uses $h=2^{-2},\ldots,2^{-9}$.
- **Methods:** MINRES with geometric-MG or AMG block-diagonal
  preconditioner, and PPCG with a constraint preconditioner.
- **Preconditioner settings:** two multigrid V-cycles; two pre- and
  post-smoothing relaxed-Jacobi steps with $\omega=8/9$; Chebyshev or cheap
  lumped-style mass approximations.
- **Sweep/stopping:** $\beta=10^{-2}$, $5\mathbin\cdot10^{-5}$, $10^{-5}$,
  and $10^{-8}$; tolerance $10^{-8}$. MINRES uses relative Euclidean
  residual; PPCG uses $r^{\mathsf T}g$ relative to its initial value.

Table 6.4 reports the following outer iterations for the three tabulated
coefficients. Each entry is geometric-MG MINRES / AMG MINRES / PPCG / PPCG
with $G=\mathrm{diag}(K)$.

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

### E6.7.2 — All-at-once diffusion-reaction control

**Source:** equations (6.89)–(6.91), Figure 6.8, book pages 204–205.

This repeats E6.7.1, including target, domain, `Q1` DtO discretisation,
preconditioners, parameter sweep, and stopping criteria, replacing the state
equation by

~~~math
-\Delta y+y=u \quad\text{in }\Omega,
\qquad y=z_{d}\quad\text{on }\partial\Omega.
~~~

Replace every stiffness action in the preconditioners by $K+M$. Figure 6.8
reports qualitatively similar iteration behaviour and no second full timing
table.

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

Take $\Omega=(0,1)^{2}$ and
$z_{d}=\sin(2\pi x_{1})\sin(2\pi x_{2})$. The source quotes the
unconstrained exact control
$u_{\mathrm{ex}}=8\pi^{2}\sin(2\pi x_{1})\sin(2\pi x_{2})$.

- **Discretisation:** conforming bilinear `Q1` elements on a quadrilateral
  mesh; mesh size is not stated.
- **Method:** PDAS from $u_{0}=0$ and multiplier $\mu_{0}=0$.
- **Sweep:** $\beta=10^{-2},10^{-4},10^{-6},10^{-8}$ and
  $a=100,70,50,30$.
- **Published trend:** $a=100$ converges in one iteration because
  $\lVert u_{\mathrm{ex}}\rVert_{L^{\infty}(\Omega)}\leq80$. Tighter bounds
  and finer meshes increase PDAS iterations. Smaller $\beta$ increases norm
  and active-region extent while lowering tracking error.

### E6.9.2 — Asymmetric box-constrained Laplace control

**Source:** equation (6.106), Figures 6.14–6.15, book pages 214–215.

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

For $d=2$ and $d=3$, use

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

- **Discretisation:** `Q1` bilinear elements in two dimensions and trilinear
  elements in three dimensions, on quadrilateral/hexahedral meshes.
- **Method/sweep:** PDAS from $u_{0}=0$, $\mu_{0}=0$ for
  $\beta=10^{-2},10^{-4},10^{-6}$.
- **Published trend:** objective, tracking error, and control norm vary with
  $\beta$ as in E6.9.1. The discontinuous lower bound makes the
  large-$\beta$ pattern sensitive to the bound jump.

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

- **Discretisation:** inf-sup compatible `Q2–Q1`; reported
  $N_{v}=51{,}842$ and $N_{\pi}=6{,}561$. The source does not state a
  pressure gauge; select and record one.
- **Method:** PDAS. The chapter introduces zero control/multiplier starts for
  its PDAS tests; preserve or explicitly replace them.
- **Naming:** equation (6.107) calls regularisation $\alpha$, while Figure
  6.16 and Table 6.5 label the same sweep $\beta$.

| Coefficient | $J_{h}(\hat u)$ | Iterations | $\lVert\hat u\rVert_{L^{2}(\Omega)}$ | $\lVert\hat v-v_{d}\rVert_{L^{2}(\Omega)}$ |
| ---: | ---: | ---: | ---: | ---: |
| $10^{-2}$ | 0.0579601 | 1 | 0.326406 | 0.338902 |
| $10^{-4}$ | 0.0430436 | 1 | 10.3743 | 0.274453 |
| $10^{-6}$ | 0.0244265 | 1 | 72.2504 | 0.208885 |
| $10^{-8}$ | 0.0176069 | 6 | 308.912 | 0.185093 |

The uncontrolled reference is $J_{h}(u_{0})=0.0585097$. The target violates
the prescribed velocity boundary data, so exact matching is impossible; the
source observes a stronger vortex shift for smaller regularisation.

## Reproduction protocol

1. Implement the residual, observation, loss, metric, constraint, and solver
   policy with derivative/KKT tests from the numerical-methods guide.
2. Reproduce the source finite-element layout as far as stated data permit.
   Record every source omission and replacement, especially forcing,
   viscosity, pressure gauge, and line-search rule.
3. Compare fields, objective components, feasibility/KKT residuals, outer and
   inner solver histories, and mesh/regularisation trends before iterations or
   time.
4. Treat qualitative trends and mesh-independent preconditioned behaviour as
   the first benchmark gate. Any numerical tolerance is a target-specific,
   versioned decision, not a generic framework test.
