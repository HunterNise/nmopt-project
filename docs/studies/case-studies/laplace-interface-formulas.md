# Laplace variations: formula deltas and interface requirements

## 1. Baseline input

$$
V=H_{0}^{1}(\Omega),\qquad
U=L^{2}(\Omega),\qquad
Q=L^{2}(\omega_{o}).
$$

$$
y\in V,\qquad
u\in U,\qquad
f\in V^{\ast},\qquad
y_{d}\in Q,\qquad
\alpha>0.
$$

$$
\langle E(y,u),v\rangle_{V^{\ast},V}
=(\nabla y,\nabla v)_{\Omega}
-\langle f,v\rangle_{V^{\ast},V}
-(u,v)_{\Omega},
\qquad v\in V.
$$

$$
J(y,u)
=\frac{1}{2}\lVert Cy-y_{d}\rVert_{Q}^{2}
+\frac{\alpha}{2}\lVert u\rVert_{U}^{2},
\qquad
C:V\to Q,\quad Cy=y\vert_{\omega_{o}}.
$$

| Component | Concrete input |
|---|---|
| Region | $\Omega$, $\Gamma=\partial\Omega$, $`\omega_{o}`$ |
| Space and pairing | $V$, $U$, $Q$, $V^{\ast}\times V$ |
| Variables | state $y$, control $u$ |
| Data | $f$, $`y_{d}`$, $\alpha$ |
| Equation block | $E:V\times U\to V^{\ast}$ |
| Residual terms | diffusion, source, source-control |
| Observation | $C:V\to Q$ |
| Losses | tracking and control regularisation |
| Metric | $`G_{U}:U\to U^{\ast}`$, if a gradient is required |

The derivative actions are

$$
\begin{aligned}
\langle D_{y}E(y,u)w,v\rangle
&=(\nabla w,\nabla v)_{\Omega},\\
\langle D_{u}E(y,u)\delta u,v\rangle
&=-(\delta u,v)_{\Omega},\\
D_{y}J(y,u)w
&=(Cy-y_{d},Cw)_{Q},\\
D_{u}J(y,u)\delta u
&=\alpha(u,\delta u)_{U}.
\end{aligned}
$$

With

$$
\mathcal L(y,u,p)=J(y,u)-\langle p,E(y,u)\rangle_{V,V^{\ast}},
$$

the solver equations are

$$
\begin{aligned}
(\nabla y,\nabla v)_{\Omega}
&=\langle f,v\rangle+(u,v)_{\Omega},
&&v\in V,\\
(\nabla w,\nabla p)_{\Omega}
&=(Cy-y_{d},Cw)_{Q},
&&w\in V,\\
j'(u)\delta u
&=\alpha(u,\delta u)_{U}+(p,\delta u)_{\Omega},
&&\delta u\in U.
\end{aligned}
$$

The search direction is

$$
\nabla_{G_{U}}j(u)=G_{U}^{-1}j'(u).
$$

## 2. Baseline compiled output

Choose

```math
Y_{h}=\mathrm{span}\left\{\varphi_{i}\right\},\qquad
Z_{h}=\mathrm{span}\left\{\psi_{j}\right\},\qquad
U_{h}=\mathrm{span}\left\{\chi_{\ell}\right\}.
```

Define

$$
(A_{h})_{ji}=(\nabla\varphi_{i},\nabla\psi_{j})_{\Omega},
\qquad
(B_{h})_{j\ell}=(\chi_{\ell},\psi_{j})_{\Omega},
\qquad
(f_{h})_{j}=\langle f,\psi_{j}\rangle.
$$

For an observation basis $`\left\{\eta_{k}\right\}`$,

$$
(W_{h})_{kr}=(\eta_{k},\eta_{r})_{Q},
\qquad
(R_{\mathrm{reg},h})_{\ell m}=(\chi_{\ell},\chi_{m})_{U}.
$$

Compilation returns

$$
r_{h}(y_{h},u_{h})=A_{h}y_{h}-f_{h}-B_{h}u_{h},
$$

$$
J_{h}(y_{h},u_{h})
=\frac{1}{2}(C_{h}y_{h}-d_{h})^{\mathsf T}W_{h}(C_{h}y_{h}-d_{h})
+\frac{\alpha}{2}u_{h}^{\mathsf T}R_{\mathrm{reg},h}u_{h},
$$

```math
E_{h}'(y_{h},u_{h})
\begin{bmatrix}
\delta y_{h}\\
\delta u_{h}
\end{bmatrix}
=A_{h}\delta y_{h}-B_{h}\delta u_{h},
\qquad
E_{h}'(y_{h},u_{h})^{\ast}p_{h}
=
\begin{bmatrix}
A_{h}^{\mathsf T}p_{h}\\
-B_{h}^{\mathsf T}p_{h}
\end{bmatrix}.
```

The DTO solver receives

$$
\begin{aligned}
A_{h}y_{h}&=f_{h}+B_{h}u_{h},\\
A_{h}^{\mathsf T}p_{h}&=C_{h}^{\mathsf T}W_{h}(C_{h}y_{h}-d_{h}),\\
j_{h}'(u_{h})&=\alpha R_{\mathrm{reg},h}u_{h}+B_{h}^{\mathsf T}p_{h},\\
\nabla_{G_{h}}j_{h}(u_{h})&=G_{h}^{-1}j_{h}'(u_{h}).
\end{aligned}
$$

The mandatory executable ports are

$$
E_{h}(x_{h}),\qquad
E_{h}'(x_{h})\delta x_{h},\qquad
E_{h}'(x_{h})^{\ast}p_{h},\qquad
J_{h}(x_{h}),\qquad
J_{h}'(x_{h}),\qquad
G_{h}^{-1}.
$$

## 3. Formula deltas

### 3.1 Reaction and transport

Replace the diffusion form by

$$
a(y,v)
=(K\nabla y,\nabla v)_{\Omega}
+(b\cdot\nabla y,v)_{\Omega}
+(cy,v)_{\Omega}.
$$

Then

$$
\langle E(y,u),v\rangle
=a(y,v)-\langle f,v\rangle-(u,v)_{\Omega},
$$

$$
\langle D_{y}E(y,u)^{\ast}p,w\rangle=a(w,p),
$$

$$
(A_{h})_{ji}=a(\varphi_{i},\psi_{j}).
$$

Required interfaces: residual term; data $K,b,c$; transport/stabilisation
policy if the chosen form needs one.

### 3.2 Neumann control

Replace the full Dirichlet boundary setting by

```math
\Gamma=\Gamma_{D}\mathbin{\dot\cup}\Gamma_{N},
\qquad
V=H_{\Gamma_{D}}^{1}(\Omega)
=\left\{v\in H^{1}(\Omega):\gamma v=0\ \text{on }\Gamma_{D}\right\},
\qquad
\Gamma_{c}\subseteq\Gamma_{N}.
```

Set

$$
U_{\Gamma}=L^{2}(\Gamma_{c}).
$$

Replace the control contribution by

$$
\langle E_{\Gamma}(y,u),v\rangle
=(\nabla y,\nabla v)_{\Omega}-\langle f,v\rangle
-\langle u,\gamma_{\Gamma_{c}}v\rangle.
$$

Then

$$
\langle D_{u}E_{\Gamma}(y,u)\delta u,v\rangle
=-\langle\delta u,\gamma_{\Gamma_{c}}v\rangle,
$$

$$
j_{\Gamma}'(u)\delta u
=\alpha(u,\delta u)_{U_{\Gamma}}
+\langle\delta u,\gamma_{\Gamma_{c}}p\rangle.
$$

With $`U_{\Gamma,h}=\mathrm{span}\left\{\xi_{\ell}\right\}`$,

$$
(B_{\Gamma,h})_{j\ell}
=\langle\xi_{\ell},\gamma_{\Gamma_{c}}\psi_{j}\rangle,
\qquad
j_{\Gamma,h}'(u_{h})
=\alpha R_{\mathrm{reg},\Gamma,h}u_{h}+B_{\Gamma,h}^{\mathsf T}p_{h}.
$$

Required interfaces: boundary partition and state space, boundary control
space, trace map, boundary residual term, boundary quadrature.

### 3.3 Boundary tracking

For a nontrivial trace observation, use the same boundary partition

$$
\Gamma=\Gamma_{D}\mathbin{\dot\cup}\Gamma_{N},
\qquad
V=H_{\Gamma_{D}}^{1}(\Omega),
\qquad
\Gamma_{o}\subseteq\Gamma_{N}.
$$

Let

$$
Q_{\Gamma}=L^{2}(\Gamma_{o}),
\qquad
C_{\Gamma} y=\gamma_{\Gamma_{o}}y.
$$

Replace the tracking loss by

$$
J_{\Gamma}(y,u)
=\frac{1}{2}\lVert C_{\Gamma} y-y_{d,\Gamma}\rVert_{Q_{\Gamma}}^{2}
+\frac{\alpha}{2}\lVert u\rVert_{U}^{2}.
$$

The state equation is unchanged.  The adjoint becomes

$$
(\nabla w,\nabla p)_{\Omega}
=\langle C_{\Gamma} y-y_{d,\Gamma},C_{\Gamma} w\rangle_{Q_{\Gamma}}.
$$

The compiled adjoint right-hand side is

$$
A_{h}^{\mathsf T}p_{h}
=C_{\Gamma,h}^{\mathsf T}W_{\Gamma,h}
(C_{\Gamma,h}y_{h}-d_{\Gamma,h}).
$$

Required interfaces: boundary partition and state space, trace observation,
boundary loss, boundary quadrature.  Once this state space has been selected,
changing only $`\Gamma_{o}`$ changes only the observation and loss.

### 3.4 Point observations

$$
C_{\mathrm{pt}}y=
\begin{bmatrix}
y(x_{1})\\
\vdots\\
y(x_{m})
\end{bmatrix}
\in\mathbb R^{m},
$$

$$
J_{\mathrm{pt}}(y,u)
=\frac{1}{2}(C_{\mathrm{pt}}y-d)^{\mathsf T}W(C_{\mathrm{pt}}y-d)
+\frac{\alpha}{2}\lVert u\rVert_{U}^{2},
$$

$$
D_{y}J_{\mathrm{pt}}(y,u)
=C_{\mathrm{pt}}^{\ast}W(C_{\mathrm{pt}}y-d).
$$

Required interfaces: observation with target $\mathbb R^{m}$, transpose
action, and a policy declaring continuous validity or discrete-only meaning.

### 3.5 $H^{1}$ regularisation

Replace the control loss by

$$
\frac{\alpha}{2}\left(
\lVert u\rVert_{L^{2}(\Omega)}^{2}
+\ell^{2}\lVert\nabla u\rVert_{L^{2}(\Omega)}^{2}
\right),
\qquad
U=H^{1}(\Omega).
$$

Then

$$
j_{\mathrm{reg}}'(u)\delta u
=\alpha\left(
(u,\delta u)_{\Omega}
+\ell^{2}(\nabla u,\nabla\delta u)_{\Omega}
\right)
+(p,\delta u)_{\Omega},
$$

$$
R_{\mathrm{reg},h}=M_{U,h}+\ell^{2}K_{U,h}.
$$

Required interfaces: control space, control loss, control discretisation, and
possibly a new constraint space.

### 3.6 $H^{1}$ or $H^{-1}$ search metric

Leave the objective unchanged.  To use a Sobolev search space while the
control remains $U=L^{2}(\Omega)$, declare

$$
U_{G}=H^{1}(\Omega),
\qquad
\iota:U_{G}\hookrightarrow U,
\qquad
\iota^{\ast}:U^{\ast}\to U_{G}^{\ast}.
$$

Define

$$
\langle G_{U}g,\delta u\rangle_{U_{G}^{\ast},U_{G}}
=(g,\delta u)_{\Omega}
+\ell^{2}(\nabla g,\nabla\delta u)_{\Omega}.
$$

Then only

$$
\nabla_{G_{U}}j=G_{U}^{-1}\iota^{\ast}j'
$$

changes.  The state, adjoint, and reduced covector $j'$ are unchanged.

Required interface: metric with apply and inverse-apply.  An $H^{-1}$
metric may implement inverse-apply by an auxiliary elliptic solve.  If the
control space itself is changed to $H^{1}(\Omega)$, the injection
$\iota^{\ast}$ is unnecessary.

### 3.7 Box constraints

```math
U_{\mathrm{ad}}
=\left\{u\in U:u_{a}\leq u\leq u_{b}\ \text{a.e. in }\Omega\right\}.
```

Replace unconstrained stationarity by

$$
0\in j'(u)+N_{U_{\mathrm{ad}}}(u),
$$

or

$$
u=\Pi_{U_{\mathrm{ad}}}\bigl(u-\tau\nabla_{G_{U}}j(u)\bigr).
$$

Required interfaces: constraint feasibility, projection, normal cone, or
multiplier operation.  No residual or adjoint formula changes.

### 3.8 Dirichlet control

Let

```math
\Gamma_{D}=\Gamma_{0}\mathbin{\dot\cup}\Gamma_{c},
\qquad
V_{0}=\left\{v\in H^{1}(\Omega):\gamma v=0\ \text{on }\Gamma_{D}\right\},
```

$$
U_{\Gamma}=H^{1/2}(\Gamma_{c}),
\qquad
L_{D}:U_{\Gamma}\to H^{1}(\Omega).
$$

With a fixed-data lifting $`\ell_{0}`$,

$$
y_{\mathrm{phys}}=\widehat y+\ell_{0}+L_{D}u,
\qquad
\widehat y\in V_{0}.
$$

The residual is

$$
\langle E_{D}(\widehat y,u),v\rangle
=(\nabla(\widehat y+\ell_{0}+L_{D}u),\nabla v)_{\Omega}
-\langle f,v\rangle,
\qquad v\in V_{0}.
$$

The objective is

$$
J_{D}(\widehat y,u)
=\frac{1}{2}\lVert C(\widehat y+\ell_{0}+L_{D}u)-y_{d}\rVert_{Q}^{2}
+\frac{\alpha}{2}\lVert u\rVert_{U_{\Gamma}}^{2}.
$$

The changed derivatives are

$$
\begin{aligned}
\langle D_{\widehat y}E_{D}w,v\rangle
&=(\nabla w,\nabla v)_{\Omega},\\
\langle D_{u}E_{D}\delta u,v\rangle
&=(\nabla L_{D}\delta u,\nabla v)_{\Omega},\\
D_{u}J_{D}\delta u
&=\langle Cy_{\mathrm{phys}}-y_{d},CL_{D}\delta u\rangle_{Q}
+\alpha(u,\delta u)_{U_{\Gamma}}.
\end{aligned}
$$

Therefore

$$
(\nabla w,\nabla p)_{\Omega}=D_{\widehat y}J_{D}w,
$$

$$
j_{D}'(u)\delta u
=D_{u}J_{D}\delta u
-(\nabla L_{D}\delta u,\nabla p)_{\Omega}.
$$

The compiled residual is

$$
r_{D,h}(\widehat y_{h},u_{h})
=A_{h}\widehat y_{h}
+A_{\mathrm{ext},h}L_{D,h}u_{h}
+b_{\ell_{0},h}-f_{h}.
$$

Required interfaces: controlled/fixed boundary regions, lifting
transformation with value/JVP/transpose-JVP, residual and observation on the
physical state, and compiler support for $`L_{D,h}`$ and affine constraints.

### 3.9 Coefficient identification

Introduce

$$
m\in M_{\mathrm{ad}},
$$

and use

$$
\langle E_{m}(y,m),v\rangle
=(m\nabla y,\nabla v)_{\Omega}-\langle f,v\rangle.
$$

The parameter derivative is

$$
\langle D_{m}E_{m}(y,m)\delta m,v\rangle
=(\delta m\nabla y,\nabla v)_{\Omega}.
$$

The adjoint and parameter derivative are

$$
(m\nabla w,\nabla p)_{\Omega}=D_{y}J_{m}(y,m)w,
$$

$$
j_{m}'(m)\delta m
=D_{m}J_{m}(y,m)\delta m
-(\delta m\nabla y,\nabla p)_{\Omega}.
$$

Required interfaces: parameter variable and constraint, nonlinear residual
term with parameter JVP/transpose-JVP, parameter loss/metric, product and
positivity policy.

### 3.10 Parabolic source control

$$
V=H_{0}^{1}(\Omega),
\qquad
Y=L^{2}(0,T;V)\cap H^{1}(0,T;V^{\ast}),
\qquad
Z=L^{2}(0,T;V).
$$

With $`y(0)=y_{0}`$,

$$
\langle E_{T}(y,u),v\rangle
=\int_{0}^{T}
\left[
\langle\dot y,v\rangle_{V^{\ast},V}
+(\nabla y,\nabla v)_{\Omega}
-\langle f,v\rangle
-(u,v)_{\Omega}
 \right]\mathrm{d}t.
$$

For

$$
J_{T}(y,u)
=\frac{1}{2}\int_{0}^{T}\lVert Cy(t)-y_{d}(t)\rVert_{Q}^{2}\mathrm{d}t
+\frac{\alpha}{2}\int_{0}^{T}\lVert u(t)\rVert_{U}^{2}\mathrm{d}t,
$$

the adjoint relation is

$$
\int_{0}^{T}
\left[
\langle\dot w,p\rangle_{V^{\ast},V}
+(\nabla w,\nabla p)_{\Omega}
 \right]\mathrm{d}t
=\int_{0}^{T}(Cy-y_{d},Cw)_{Q}\mathrm{d}t.
$$

For $w(0)=0$,

$$
-\dot p-\Delta p=C^{\ast}(Cy-y_{d}),
\qquad
p(T)=0.
$$

The reduced derivative is

$$
j_{T}'(u)\delta u
=\alpha\int_{0}^{T}(u,\delta u)_{U}\mathrm{d}t
+\int_{0}^{T}(p,\delta u)_{\Omega}\mathrm{d}t.
$$

Required interfaces: time-space descriptors, time-derivative residual term,
initial-trace data/constraint, time loss, temporal compiler, and exact
transpose of the compiled time residual.

## 4. DTO and OTD

### Discretize then optimize

Differentiate

$$
(E_{h},J_{h}).
$$

The discrete adjoint is

$$
E_{h}'(x_{h})^{\ast}p_{h}=D_{y}J_{h}(x_{h}).
$$

Required solver inputs:

$$
E_{h}(x_{h}),\quad
E_{h}'(x_{h})\delta x_{h},\quad
E_{h}'(x_{h})^{\ast}p_{h},\quad
J_{h}(x_{h}),\quad
J_{h}'(x_{h}),\quad
G_{h}^{-1},
$$

plus constraint operations.

### Optimize then discretize

Differentiate

$$
(E,J)
$$

first:

$$
E_{y}'(x)^{\ast}p=D_{y}J(x).
$$

Compile that equation separately.  If it gives

$$
\widetilde A_{h}p_{h}=\widetilde b_{h},
$$

there is no general identity

$$
\widetilde A_{h}=A_{h}^{\mathsf T}.
$$

They can agree for compatible conforming Galerkin forms.  They need not agree
for Petrov–Galerkin, stabilisation, time stepping, inexact quadrature,
discrete observations, or different lifting treatments.

Required additional interfaces:

| Component | Required action |
|---|---|
| Formulation builder | build $`E_{y}'(x)^{\ast}p-D_{y}J(x)=0`$ as a semantic equation block |
| Provenance | record DTO or OTD and the lowering policy |

## 5. Minimal interface specification

| Interface | Mathematical contract |
|---|---|
| Space | domain/region, field shape, trial/test role, pairing, trace and product capabilities |
| Variable and data | identifier, space, value, admissible-set reference; derivatives only for variables |
| Map | $T:X\to Y$, $T(x)$, $T'(x)\delta x$, $T'(x)^{\ast}\eta$ |
| Residual term | $`\langle E_{t}(x),z\rangle`$, JVP, transpose-JVP |
| Equation block | test space, sum of terms, residual/JVP/transpose-JVP |
| Transformation | map applied before residual/observation evaluation |
| Observation | map from physical variables to $Q$ |
| Loss | $\Phi:Q\to\mathbb R$, value and $D\Phi(q)\in Q^{\ast}$ |
| Metric | $G:X\to X^{\ast}$, apply and inverse-apply |
| Constraint | feasibility, projection, normal cone, or multiplier operation |
| Requirement policy | trace, product, nullspace, point, fractional, or discrete-only policy |
| Compiler | choose spaces, lower maps, preserve pairing-aware JVP and transpose-JVP |
| Formulation builder | construct reduced/KKT systems and optional OTD equations |
