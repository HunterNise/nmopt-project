# Laplace variations: formula deltas and interface requirements

## 1. Baseline input

$$
V=H_0^1(\Omega),\qquad
U=L^2(\Omega),\qquad
Q=L^2(\omega_o).
$$

$$
y\in V,\qquad
u\in U,\qquad
f\in V^{\ast},\qquad
y_d\in Q,\qquad
\alpha>0.
$$

$$
\langle E(y,u),v\rangle_{V^{\ast},V}
=(\nabla y,\nabla v)_\Omega
-\langle f,v\rangle_{V^{\ast},V}
-(u,v)_\Omega,
\qquad v\in V.
$$

$$
J(y,u)
=\frac12\lVert Cy-y_d\rVert_Q^2
+\frac{\alpha}{2}\lVert u\rVert_U^2,
\qquad
C:V\to Q,\quad Cy=y|_{\omega_o}.
$$

| Component | Concrete input |
|---|---|
| Region | $\Omega$, $\Gamma=\partial\Omega$, $`\omega_o`$ |
| Space and pairing | $V$, $U$, $Q$, $V^{\ast}\times V$ |
| Variables | state $y$, control $u$ |
| Data | $f$, $`y_d`$, $\alpha$ |
| Equation block | $E:V\times U\to V^{\ast}$ |
| Residual terms | diffusion, source, source-control |
| Observation | $C:V\to Q$ |
| Losses | tracking and control regularisation |
| Metric | $`G_U:U\to U^{\ast}`$, if a gradient is required |

The derivative actions are

$$
\begin{aligned}
\langle D_yE(y,u)w,v\rangle
&=(\nabla w,\nabla v)_\Omega,\\
\langle D_uE(y,u)\delta u,v\rangle
&=-(\delta u,v)_\Omega,\\
D_yJ(y,u)w
&=(Cy-y_d,Cw)_Q,\\
D_uJ(y,u)\delta u
&=\alpha(u,\delta u)_U.
\end{aligned}
$$

With

$$
\mathcal L(y,u,p)=J(y,u)-\langle p,E(y,u)\rangle_{V,V^{\ast}},
$$

the solver equations are

$$
\begin{aligned}
(\nabla y,\nabla v)_\Omega
&=\langle f,v\rangle+(u,v)_\Omega,
&&v\in V,\\
(\nabla w,\nabla p)_\Omega
&=(Cy-y_d,Cw)_Q,
&&w\in V,\\
j'(u)\delta u
&=\alpha(u,\delta u)_U+(p,\delta u)_\Omega,
&&\delta u\in U.
\end{aligned}
$$

The search direction is

$$
\nabla_{G_U}j(u)=G_U^{-1}j'(u).
$$

## 2. Baseline compiled output

Choose

$$
Y_{h}=\mathrm{span}\left\{\varphi_{i}\right\},\qquad
Z_{h}=\mathrm{span}\left\{\psi_{j}\right\},\qquad
U_{h}=\mathrm{span}\left\{\chi_{\ell}\right\}.
$$

Define

$$
(A_h)_{ji}=(\nabla\varphi_i,\nabla\psi_j)_\Omega,
\qquad
(B_h)_{j\ell}=(\chi_\ell,\psi_j)_\Omega,
\qquad
(f_h)_j=\langle f,\psi_j\rangle.
$$

For an observation basis $`\{\eta_k\}`$,

$$
(W_h)_{kr}=(\eta_k,\eta_r)_Q,
\qquad
(R_{\mathrm{reg},h})_{\ell m}=(\chi_\ell,\chi_m)_U.
$$

Compilation returns

$$
r_h(y_h,u_h)=A_hy_h-f_h-B_hu_h,
$$

$$
J_h(y_h,u_h)
=\frac12(C_hy_h-d_h)^{\mathsf T}W_h(C_hy_h-d_h)
+\frac{\alpha}{2}u_h^{\mathsf T}R_{\mathrm{reg},h}u_h,
$$

$$
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
$$

The DTO solver receives

$$
\begin{aligned}
A_hy_h&=f_h+B_hu_h,\\
A_h^{\mathsf T}p_h&=C_h^{\mathsf T}W_h(C_hy_h-d_h),\\
j_h'(u_h)&=\alpha R_{\mathrm{reg},h}u_h+B_h^{\mathsf T}p_h,\\
\nabla_{G_h}j_h(u_h)&=G_h^{-1}j_h'(u_h).
\end{aligned}
$$

The mandatory executable ports are

$$
E_h(x_h),\qquad
E_h'(x_h)\delta x_h,\qquad
E_h'(x_h)^{\ast}p_h,\qquad
J_h(x_h),\qquad
J_h'(x_h),\qquad
G_h^{-1}.
$$

## 3. Formula deltas

### 3.1 Reaction and transport

Replace the diffusion form by

$$
a(y,v)
=(K\nabla y,\nabla v)_\Omega
+(b\cdot\nabla y,v)_\Omega
+(cy,v)_\Omega.
$$

Then

$$
\langle E(y,u),v\rangle
=a(y,v)-\langle f,v\rangle-(u,v)_\Omega,
$$

$$
\langle D_yE(y,u)^{\ast}p,w\rangle=a(w,p),
$$

$$
(A_h)_{ji}=a(\varphi_i,\psi_j).
$$

Required interfaces: residual term; data $K,b,c$; transport/stabilisation
policy if the chosen form needs one.

### 3.2 Neumann control

Replace the full Dirichlet boundary setting by

$$
\Gamma=\Gamma_D\mathbin{\dot\cup}\Gamma_N,
\qquad
V=H_{\Gamma_{D}}^1(\Omega)
=\left\{v\in H^1(\Omega):\gamma v=0\ \text{on }\Gamma_{D}\right\},
\qquad
\Gamma_c\subseteq\Gamma_N.
$$

Set

$$
U_\Gamma=L^2(\Gamma_c).
$$

Replace the control contribution by

$$
\langle E_\Gamma(y,u),v\rangle
=(\nabla y,\nabla v)_\Omega-\langle f,v\rangle
-\langle u,\gamma_{\Gamma_c}v\rangle.
$$

Then

$$
\langle D_uE_\Gamma(y,u)\delta u,v\rangle
=-\langle\delta u,\gamma_{\Gamma_c}v\rangle,
$$

$$
j_\Gamma'(u)\delta u
=\alpha(u,\delta u)_{U_\Gamma}
+\langle\delta u,\gamma_{\Gamma_c}p\rangle.
$$

With $`U_{\Gamma,h}=\mathrm{span}\{\xi_\ell\}`$,

$$
(B_{\Gamma,h})_{j\ell}
=\langle\xi_\ell,\gamma_{\Gamma_c}\psi_j\rangle,
\qquad
j_{\Gamma,h}'(u_h)
=\alpha R_{\mathrm{reg},\Gamma,h}u_h+B_{\Gamma,h}^{\mathsf T}p_h.
$$

Required interfaces: boundary partition and state space, boundary control
space, trace map, boundary residual term, boundary quadrature.

### 3.3 Boundary tracking

For a nontrivial trace observation, use the same boundary partition

$$
\Gamma=\Gamma_D\mathbin{\dot\cup}\Gamma_N,
\qquad
V=H_{\Gamma_D}^1(\Omega),
\qquad
\Gamma_o\subseteq\Gamma_N.
$$

Let

$$
Q_\Gamma=L^2(\Gamma_o),
\qquad
C_\Gamma y=\gamma_{\Gamma_o}y.
$$

Replace the tracking loss by

$$
J_\Gamma(y,u)
=\frac12\lVert C_\Gamma y-y_{d,\Gamma}\rVert_{Q_\Gamma}^2
+\frac{\alpha}{2}\lVert u\rVert_U^2.
$$

The state equation is unchanged.  The adjoint becomes

$$
(\nabla w,\nabla p)_\Omega
=\langle C_\Gamma y-y_{d,\Gamma},C_\Gamma w\rangle_{Q_\Gamma}.
$$

The compiled adjoint right-hand side is

$$
A_h^{\mathsf T}p_h
=C_{\Gamma,h}^{\mathsf T}W_{\Gamma,h}
(C_{\Gamma,h}y_h-d_{\Gamma,h}).
$$

Required interfaces: boundary partition and state space, trace observation,
boundary loss, boundary quadrature.  Once this state space has been selected,
changing only $`\Gamma_o`$ changes only the observation and loss.

### 3.4 Point observations

$$
C_\mathrm{pt}y=
\begin{bmatrix}
y(x_1)\\
\vdots\\
y(x_m)
\end{bmatrix}
\in\mathbb R^m,
$$

$$
J_\mathrm{pt}(y,u)
=\frac12(C_\mathrm{pt}y-d)^{\mathsf T}W(C_\mathrm{pt}y-d)
+\frac{\alpha}{2}\lVert u\rVert_U^2,
$$

$$
D_yJ_\mathrm{pt}(y,u)
=C_\mathrm{pt}^{\ast}W(C_\mathrm{pt}y-d).
$$

Required interfaces: observation with target $\mathbb R^m$, transpose
action, and a policy declaring continuous validity or discrete-only meaning.

### 3.5 $H^1$ regularisation

Replace the control loss by

$$
\frac{\alpha}{2}\left(
\lVert u\rVert_{L^2(\Omega)}^2
+\ell^2\lVert\nabla u\rVert_{L^2(\Omega)}^2
\right),
\qquad
U=H^1(\Omega).
$$

Then

$$
j_\mathrm{reg}'(u)\delta u
=\alpha\left(
(u,\delta u)_\Omega
+\ell^2(\nabla u,\nabla\delta u)_\Omega
\right)
+(p,\delta u)_\Omega,
$$

$$
R_{\mathrm{reg},h}=M_{U,h}+\ell^2K_{U,h}.
$$

Required interfaces: control space, control loss, control discretisation, and
possibly a new constraint space.

### 3.6 $H^1$ or $H^{-1}$ search metric

Leave the objective unchanged.  To use a Sobolev search space while the
control remains $U=L^2(\Omega)$, declare

$$
U_{G}=H^1(\Omega),
\qquad
\iota:U_{G}\hookrightarrow U,
\qquad
\iota^{\ast}:U^{\ast}\to U_{G}^{\ast}.
$$

Define

$$
\langle G_Ug,\delta u\rangle_{U_G^{\ast},U_G}
=(g,\delta u)_\Omega
+\ell^2(\nabla g,\nabla\delta u)_\Omega.
$$

Then only

$$
\nabla_{G_U}j=G_U^{-1}\iota^{\ast}j'
$$

changes.  The state, adjoint, and reduced covector $j'$ are unchanged.

Required interface: metric with apply and inverse-apply.  An $H^{-1}$
metric may implement inverse-apply by an auxiliary elliptic solve.  If the
control space itself is changed to $H^1(\Omega)$, the injection
$\iota^{\ast}$ is unnecessary.

### 3.7 Box constraints

$$
U_{\mathrm{ad}}
=\left\{u\in U:u_{a}\leq u\leq u_{b}\ \text{a.e. in }\Omega\right\}.
$$

Replace unconstrained stationarity by

$$
0\in j'(u)+N_{U_\mathrm{ad}}(u),
$$

or

$$
u=\Pi_{U_\mathrm{ad}}\bigl(u-\tau\nabla_{G_U}j(u)\bigr).
$$

Required interfaces: constraint feasibility, projection, normal cone, or
multiplier operation.  No residual or adjoint formula changes.

### 3.8 Dirichlet control

Let

$$
\Gamma_{D}=\Gamma_{0}\mathbin{\dot\cup}\Gamma_{c},
\qquad
V_{0}=\left\{v\in H^1(\Omega):\gamma v=0\ \text{on }\Gamma_{D}\right\},
$$

$$
U_{\Gamma}=H^{1/2}(\Gamma_{c}),
\qquad
L_D:U_\Gamma\to H^1(\Omega).
$$

With a fixed-data lifting $`\ell_0`$,

$$
y_\mathrm{phys}=\widehat y+\ell_0+L_Du,
\qquad
\widehat y\in V_0.
$$

The residual is

$$
\langle E_D(\widehat y,u),v\rangle
=(\nabla(\widehat y+\ell_0+L_Du),\nabla v)_\Omega
-\langle f,v\rangle,
\qquad v\in V_0.
$$

The objective is

$$
J_D(\widehat y,u)
=\frac12\lVert C(\widehat y+\ell_0+L_Du)-y_d\rVert_Q^2
+\frac{\alpha}{2}\lVert u\rVert_{U_\Gamma}^2.
$$

The changed derivatives are

$$
\begin{aligned}
\langle D_{\widehat y}E_Dw,v\rangle
&=(\nabla w,\nabla v)_\Omega,\\
\langle D_uE_D\delta u,v\rangle
&=(\nabla L_D\delta u,\nabla v)_\Omega,\\
D_uJ_D\delta u
&=\langle Cy_\mathrm{phys}-y_d,CL_D\delta u\rangle_Q
+\alpha(u,\delta u)_{U_\Gamma}.
\end{aligned}
$$

Therefore

$$
(\nabla w,\nabla p)_\Omega=D_{\widehat y}J_Dw,
$$

$$
j_D'(u)\delta u
=D_uJ_D\delta u
-(\nabla L_D\delta u,\nabla p)_\Omega.
$$

The compiled residual is

$$
r_{D,h}(\widehat y_h,u_h)
=A_h\widehat y_h
+A_{\mathrm{ext},h}L_{D,h}u_h
+b_{\ell_0,h}-f_h.
$$

Required interfaces: controlled/fixed boundary regions, lifting
transformation with value/JVP/transpose-JVP, residual and observation on the
physical state, and compiler support for $`L_{D,h}`$ and affine constraints.

### 3.9 Coefficient identification

Introduce

$$
m\in M_\mathrm{ad},
$$

and use

$$
\langle E_m(y,m),v\rangle
=(m\nabla y,\nabla v)_\Omega-\langle f,v\rangle.
$$

The parameter derivative is

$$
\langle D_mE_m(y,m)\delta m,v\rangle
=(\delta m\nabla y,\nabla v)_\Omega.
$$

The adjoint and parameter derivative are

$$
(m\nabla w,\nabla p)_\Omega=D_yJ_m(y,m)w,
$$

$$
j_m'(m)\delta m
=D_mJ_m(y,m)\delta m
-(\delta m\nabla y,\nabla p)_\Omega.
$$

Required interfaces: parameter variable and constraint, nonlinear residual
term with parameter JVP/transpose-JVP, parameter loss/metric, product and
positivity policy.

### 3.10 Parabolic source control

$$
V=H_0^1(\Omega),
\qquad
Y=L^2(0,T;V)\cap H^1(0,T;V^{\ast}),
\qquad
Z=L^2(0,T;V).
$$

With $`y(0)=y_0`$,

$$
\langle E_T(y,u),v\rangle
=\int_0^{T}
\left[
\langle\dot y,v\rangle_{V^{\ast},V}
+(\nabla y,\nabla v)_\Omega
-\langle f,v\rangle
-(u,v)_\Omega
\right]dt.
$$

For

$$
J_T(y,u)
=\frac12\int_0^{T}\lVert Cy(t)-y_d(t)\rVert_Q^2dt
+\frac{\alpha}{2}\int_0^{T}\lVert u(t)\rVert_U^2dt,
$$

the adjoint relation is

$$
\int_0^{T}
\left[
\langle\dot w,p\rangle_{V^{\ast},V}
+(\nabla w,\nabla p)_\Omega
\right]dt
=\int_0^{T}(Cy-y_d,Cw)_Qdt.
$$

For $w(0)=0$,

$$
-\dot p-\Delta p=C^{\ast}(Cy-y_d),
\qquad
p(T)=0.
$$

The reduced derivative is

$$
j_T'(u)\delta u
=\alpha\int_0^{T}(u,\delta u)_Udt
+\int_0^{T}(p,\delta u)_\Omega\mathrm{d}t.
$$

Required interfaces: time-space descriptors, time-derivative residual term,
initial-trace data/constraint, time loss, temporal compiler, and exact
transpose of the compiled time residual.

## 4. DTO and OTD

### Discretize then optimize

Differentiate

$$
(E_h,J_h).
$$

The discrete adjoint is

$$
E_h'(x_h)^{\ast}p_h=D_yJ_h(x_h).
$$

Required solver inputs:

$$
E_h(x_h),\quad
E_h'(x_h)\delta x_h,\quad
E_h'(x_h)^{\ast}p_h,\quad
J_h(x_h),\quad
J_h'(x_h),\quad
G_h^{-1},
$$

plus constraint operations.

### Optimize then discretize

Differentiate

$$
(E,J)
$$

first:

$$
E_y'(x)^{\ast}p=D_yJ(x).
$$

Compile that equation separately.  If it gives

$$
\widetilde A_hp_h=\widetilde b_h,
$$

there is no general identity

$$
\widetilde A_h=A_h^{\mathsf T}.
$$

They can agree for compatible conforming Galerkin forms.  They need not agree
for Petrov–Galerkin, stabilisation, time stepping, inexact quadrature,
discrete observations, or different lifting treatments.

Required additional interfaces:

| Component | Required action |
|---|---|
| Formulation builder | build $`E_y'(x)^{\ast}p-D_yJ(x)=0`$ as a semantic equation block |
| Provenance | record DTO or OTD and the lowering policy |

## 5. Minimal interface specification

| Interface | Mathematical contract |
|---|---|
| Space | domain/region, field shape, trial/test role, pairing, trace and product capabilities |
| Variable and data | identifier, space, value, admissible-set reference; derivatives only for variables |
| Map | $T:X\to Y$, $T(x)$, $T'(x)\delta x$, $T'(x)^{\ast}\eta$ |
| Residual term | $`\langle E_t(x),z\rangle`$, JVP, transpose-JVP |
| Equation block | test space, sum of terms, residual/JVP/transpose-JVP |
| Transformation | map applied before residual/observation evaluation |
| Observation | map from physical variables to $Q$ |
| Loss | $\Phi:Q\to\mathbb R$, value and $D\Phi(q)\in Q^{\ast}$ |
| Metric | $G:X\to X^{\ast}$, apply and inverse-apply |
| Constraint | feasibility, projection, normal cone, or multiplier operation |
| Requirement policy | trace, product, nullspace, point, fractional, or discrete-only policy |
| Compiler | choose spaces, lower maps, preserve pairing-aware JVP and transpose-JVP |
| Formulation builder | construct reduced/KKT systems and optional OTD equations |
