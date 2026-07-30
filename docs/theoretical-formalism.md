# Formal mathematical model and the strong-to-variational bridge

## Purpose and boundary

This note fixes the theoretical input language of the project.  It is more
general than the first implementation slice, but deliberately less ambitious
than a general-purpose theorem prover or symbolic PDE system.

The library is intended to compose a *declared formulation* of a
PDE-constrained optimisation or inverse problem and then lower that
formulation to a discrete executable problem.  The user remains responsible
for well-posedness, regularity, identifiability, and for choosing an
appropriate formulation and discretisation.  The library may check that a
component's declared requirements and policies have been supplied; it must
not present that check as a proof of analysis.

In particular, a strong PDE as printed in a textbook is useful source
material, but is not by itself an executable input.  The program needs the
spaces, pairings, boundary interpretation, variational residual, and
derivative actions that make that PDE meaningful in the selected setting.

## Conventions

Let $\Omega$ be a spatial domain, $\Gamma = \partial \Omega$ its boundary, and,
for an evolution problem, $I = (0,T)$.  A boundary may be partitioned into
named, disjoint portions such as

$$
  \Gamma = \overline{\Gamma_D}\mathbin\cup
           \overline{\Gamma_N}\mathbin\cup
           \overline{\Gamma_R}.
$$

The notation

$$
  \langle \xi,z\rangle_{Z^*,Z}
$$

always denotes the dual pairing between a test space $Z$ and its dual.  It is
not silently replaced by an $L^2$ inner product.  When both factors happen to
be $L^2$ functions, the pairing can of course be represented by the $L^2$
integral.

We use the Lagrangian convention

$$
  \mathcal L(x,p;d) = J(x;d)
    - \sum_{a\in\mathcal A}\langle p_a,E_a(x;d)\rangle_{Z_a,Z_a^*}.
$$

Here $d$ denotes fixed data.  The sign is global: individual terms and
solvers must not select a different convention.

## 1. The abstract optimisation problem

### 1.1 Variables, data, equations, and constraints

Let $\mathcal I$ be a finite set of variable blocks.  A block may represent a
state, a control, a parameter, a flux, an initial value, or an auxiliary
mixed variable.  For every $i \in \mathcal I$, declare a topological vector
space $X_i$; the full unknown is

$$
  x=(x_i)_{i\in\mathcal I}\in X:=\prod_{i\in\mathcal I}X_i.
$$

Examples are $H^1(\Omega)$, $L^2(\Omega)^m$, $H(div;\Omega)$,
$L^2(\Gamma_c)$, a finite-dimensional Euclidean space, and a Bochner space
such as $L^2(I;V)$.  Scalar, vector, and tensor fields differ only in the
field shape carried by their space; they are not different kinds of problem.

Fixed coefficients, source functions, desired states, geometry tags, and
initial data belong to data spaces $D_r$.  With

$$
  d=(d_r)_{r\in\mathcal R}\in D:=\prod_{r\in\mathcal R}D_r,
$$

the residual is allowed to depend on data,

$$
  E_a:X\times D\longrightarrow Z_a^*,\qquad a\in\mathcal A,
$$

where $Z_a$ is the declared test space for equation block $a$.  The
mathematical problem is

$$
  \begin{aligned}
    &\min_{x\in X_{\rm ad}(d)} J(x;d),\\
    &E_a(x;d)=0\quad\text{in }Z_a^*,\qquad a\in\mathcal A.
  \end{aligned} \tag{1}
$$

The residual equation means

$$
  \langle E_a(x;d),z_a\rangle_{Z_a^*,Z_a}=0
  \quad\text{for every }z_a\in Z_a. \tag{2}
$$

Thus a mixed system simply has several equation blocks and test spaces.  A
Petrov--Galerkin formulation is already covered: its trial spaces $X_i$ and
test spaces $Z_a$ need not coincide.

$X_{\rm ad}(d)$ contains optimisation constraints, such as box constraints,
affine restrictions, a mean-zero gauge, or positivity of a coefficient.  An
essential boundary condition is normally **not** stored as a load term in
$E$: it restricts or parameterises a physical trial space, as described in
Section 2.3.

Data are not fake variables.  For example, the affine residual

$$
  E(y,u;f) = A y - B u - F(f)
$$

has $f$ as immutable data.  If a source or coefficient is to be estimated,
it is explicitly promoted from $d$ to one block of $x$; then its derivative
and admissible set become part of the problem.

### 1.2 Terms, observations, and objective functionals

At the semantic level, an equation is a sum of independently meaningful
terms:

$$
  E_a(x;d)=\sum_{\ell\in\mathcal T_a}E_{a\ell}(x;d).
$$

For a weak finite-element formulation, a term need only define its tested
action

$$
  (x,d,z_a)\longmapsto
  \langle E_{a\ell}(x;d),z_a\rangle.
$$

It may be a diffusion integral, a reaction integral, a boundary functional,
a mixed coupling, a time derivative, a stabilization term, or a lifting
contribution.  Calling it a `diffusion` or `Robin` term is useful metadata;
its typed map and tested action are the mathematical contract.

Objectives are factored through observations:

$$
  O_k:X\times D\longrightarrow Q_k,
  \qquad
  J(x;d)=\sum_{k\in\mathcal K}\Phi_k(O_k(x;d);d). \tag{3}
$$

The observation space $Q_k$ is explicit.  It can be a volume space, a trace
space, a flux space, a space-time space, or $R^m$ for a finite sensor array.
The loss $\Phi_k$ supplies a scalar value and derivative in $Q_k^*$.  Hence a
tracking term is not intrinsically an $L^2$ mass matrix.  For example,

$$
  \tfrac12\|O(y)-y_d\|_{Q}^2
$$

requires a declared norm or metric on $Q$; a weighted $L2$, an energy norm,
a boundary norm, and a discrete sensor covariance are different choices.

### 1.3 Derivatives, adjoints, and gradients

At a differentiability point, a residual term provides partial linearizations

$$
  D_iE_a(x;d):X_i\longrightarrow Z_a^*,
  \qquad
  D_iE_a(x;d)^*:Z_a\longrightarrow X_i^*. \tag{4}
$$

The latter is the transpose action with respect to the declared dual
pairings.  An objective provides

$$
  D_iJ(x;d)\in X_i^*.
$$

For a state/control split $x=(y,u)$ and adjoints $p_a \in Z_a$, the first
order pattern for an unconstrained state is

$$
  \begin{aligned}
  E_a(y,u;d)&=0,\\
  \sum_a D_yE_a(y,u;d)^*p_a&=D_yJ(y,u;d),\\
  j'(u)&=D_uJ(y,u;d)-\sum_aD_uE_a(y,u;d)^*p_a\in U^*.
  \end{aligned} \tag{5}
$$

With an admissible control set, the last line is replaced by the appropriate
variational inequality, normal-cone relation, projection, or multiplier KKT
condition.  A reduced gradient is not obtained until a primal--dual
identification is chosen.  In a Hilbert control space, a metric/Riesz map

$$
  R_U:U\longrightarrow U^*,\qquad \nabla_U j=R_U^{-1}j'(u) \tag{6}
$$

defines one such identification.  Choosing $L^2$, $H^1$, boundary $L^2$, or an
auxiliary-solve $H^{-1}$ gradient changes (6), not the residual or the adjoint
equation.  A fractional metric requires its own declared realization.

### 1.4 Transformations and physical fields

The variable in (1) need not be the physical field appearing in a PDE.  A
reconstruction or transformation is a map

$$
  T:\widehat X\times D\longrightarrow X_{\rm phys}.
$$

The declared residual is then the composite

$$
  \widehat E(\widehat x;d)=E(T(\widehat x;d);d). \tag{7}
$$

Liftings, trace maps, restrictions to a subdomain, coefficient
parameterizations, control-to-actuator maps, and transfers between meshes all
fit this pattern.  Their derivatives are composed by the chain rule.  This
is the correct home for essential boundary conditions and Dirichlet control.

## 2. From a textbook strong form to an input formulation

### 2.1 The strong form is not enough

A textbook may begin with

$$
  -\nabla\!\cdot(K\nabla y)+b\!\cdot\nabla y+c y
  = f+\chi_{\omega_c}u\quad\hbox{in }\Omega, \tag{8}
$$

followed by boundary equalities.  Equation (8) alone does not determine a
library input.  It leaves open, among other things:

- whether derivatives are interpreted classically, distributionally, or by
  transposition;
- the state and test spaces, and whether the form is coercive, inf--sup, or
  stabilized;
- which boundary values are essential restrictions and which are boundary
  functionals;
- whether transport is integrated by parts and how inflow/outflow is treated;
- the control and data spaces that make the displayed products and traces
  meaningful; and
- the discrete scheme, quadrature, and treatment of nonconforming or
  Petrov--Galerkin test functions.

The library should therefore not accept a string such as `diffusion-
transport-reaction` and infer all of this.  A strong form can be retained as
documentation or as a future symbolic front end, but compilation starts from
the selected residual formulation.

### 2.2 Required bridge steps

To turn a strong statement into an abstract problem, the model author makes
the following mathematical choices.

1. **Declare geometry, regions, fields, and data.**  Name volume and
   boundary regions, choose field shapes, and state which symbols are
   optimised and which are fixed data.
2. **Choose trial and test spaces.**  This includes regularity, traces,
   duals, product spaces, and temporal topology if relevant.
3. **Handle essential conditions.**  Restrict the trial space to homogeneous
   traces or introduce a lifting/reconstruction.
4. **Choose the formulation.**  Multiply by test functions, integrate by
   parts only where intended, and record every surviving volume, boundary,
   initial, and terminal contribution.
5. **Declare each residual term and its linearized/transpose actions.**
   Built-in terms may supply these; a user-defined term must do so.
6. **Declare observations, losses, metrics, and admissible-set operations.**
7. **State unverified analysis requirements and exceptional policies.**
   Examples are a trace/product requirement, a pure-Neumann compatibility
   condition, a point-sensor regularity claim, or a chosen fractional-norm
   realization.

Steps 2--4 are mathematical modelling choices, not defaults that the program
can safely guess.

### 2.3 Essential, natural, and Robin conditions

Suppose $\Gamma_D$ carries a Dirichlet datum $g_D$.  Set

$$
  V:=\{v\in H^1(\Omega):\operatorname{tr}_{\Gamma_D}v=0\}
$$

and choose a lifting $\ell_D(g_D)$ with the desired trace.  The physical
state is represented as

$$
  y=\widehat y+\ell_D(g_D),\qquad \widehat y\in V. \tag{9}
$$

Thus an essential condition modifies the trial variable or its
reconstruction.  It is not an ordinary functional on the test space.

By contrast, a Neumann datum is a natural boundary functional, for example
$v \mapsto \int_{\Gamma_N} g_N v$.  A Robin condition contributes both a
boundary bilinear term and, usually, a boundary functional.  A model with a
pure Neumann operator must additionally provide a compatibility and
nullspace/gauge policy.  Periodicity and hanging-node relations are likewise
constraints on a discrete realization, rather than physical load terms.

### 2.4 A very-weak formulation is an explicit different residual

Even a simple strong equation can lead to genuinely different semantic
models.  For illustration, take

$$
  -\Delta y=f\quad\hbox{in }\Omega,\qquad y=g\quad\hbox{on }\Gamma.
$$

The ordinary weak approach seeks a lifting plus an $H^1_0(\Omega)$ unknown.
If the available data or desired solution regularity instead calls for a
transposition formulation, one possible model seeks

$$
  y\in Y=L^2(\Omega),\qquad
  Z=H^2(\Omega)\cap H^1_0(\Omega),
$$

and declares the residual

$$
  \langle E(y;f,g),z\rangle
   =\int_\Omega y(-\Delta z)-\int_\Omega fz
     +\langle g,\partial_nz\rangle_{\Gamma}=0
  \quad\forall z\in Z.
$$

The final pairing in this formula must be meaningful in the model author's selected
trace spaces.  Domain regularity and the exact test space may also need to be
changed.  This is not a variant the compiler can discover from the strong
equation: $Y$, $Z$, the normal-derivative trace, and the data pairing are
different from the weak $H^1$ residual.  Both fit (1), but only after their
respective maps and requirements are declared.

## 3. Worked recovery of representative cases

The examples below use the same abstract ingredients.  Their substitutions
make visible what changes locally and what must not change in the core.

### 3.1 Stationary diffusion--transport--reaction with volume control

Consider the strong model

$$
\begin{aligned}
 -\nabla\!\cdot(K\nabla y)+b\!\cdot\nabla y+c y
     &= f+\chi_{\omega_c}u &&\text{in }\Omega,\\
 y&=g_D &&\text{on }\Gamma_D,\\
 (K\nabla y)\!\cdot n&=g_N &&\text{on }\Gamma_N,\\
 (K\nabla y)\!\cdot n+\rho y&=g_R &&\text{on }\Gamma_R.
\end{aligned} \tag{10}
$$

One possible primal weak formulation has the following substitutions into
(1):

| Abstract component | This problem |
|---|---|
| Optimisation variables | $x=(\hat y,u)$ |
| State trial space | $Y=V={v \in H^1(\Omega): tr_{\Gamma_D} v=0}$ |
| Control space | $U=L^2(\omega_c)$ |
| Test space | $Z=V$ |
| Fixed data | $d=(K,b,c,f,g_D,g_N,\rho,g_R,y_d,\alpha)$ |
| Reconstruction | $y=\hat y+\ell_D(g_D)$ |
| Equation | One residual $E(\hat y,u;d) \in V^*$ |

If $\Gamma_D=\Gamma$, then $V=H^1_0(\Omega)$ and this residual lies in
$V^*=H^{-1}(\Omega)$.

For $v \in V$, define

$$
\begin{aligned}
 \langle E(\widehat y,u;d),v\rangle
 ={}&\int_\Omega K\nabla(\widehat y+\ell_D)\cdot\nabla v
       +(b\cdot\nabla(\widehat y+\ell_D))v
       +c(\widehat y+\ell_D)v\\
 &+\int_{\Gamma_R}\rho(\widehat y+\ell_D)v
   -\int_\Omega fv-\int_{\omega_c}uv
   -\int_{\Gamma_N}g_Nv-\int_{\Gamma_R}g_Rv .
\end{aligned} \tag{11}
$$

Equation (11), not (10), is the residual input.  It exposes the chosen
transport convention: the transport derivative was left on the state.  An
equally legitimate adjoint-oriented or stabilized formulation is a different
declared residual and may have different trace requirements.

For distributed tracking and $L^2$ control regularisation, select

$$
\begin{aligned}
 Q_\Omega&=L^2(\omega_o),& O_\Omega(\widehat y,u)&=(\widehat y+\ell_D)|_{\omega_o},\\
 Q_u&=L^2(\omega_c),& O_u(\widehat y,u)&=u,\\
 J&=\tfrac12\|O_\Omega-y_d\|_{L^2(\omega_o)}^2+
       \tfrac\alpha2\|u\|_{L^2(\omega_c)}^2. \tag{12}
\end{aligned}
$$

The linearized residual in directions $(w,\delta u)$ is obtained from (11):

$$
  \langle E'(\widehat y,u)(w,\delta u),v\rangle
   =a(w,v)-\int_{\omega_c}\delta u\,v, \tag{13}
$$

where $a$ is the bilinear part of (11).  The adjoint is defined by the
transpose of exactly this form:

$$
  a(w,p)=\int_{\omega_o}(y-y_d)w\quad\hbox{for every }w\in V. \tag{14}
$$

Under the global Lagrangian convention, the reduced derivative is

$$
  j'(u)\delta u=\alpha\int_{\omega_c}u\,\delta u+
                       \int_{\omega_c}p\,\delta u. \tag{15}
$$

Only after selecting, for example, the $L^2(\omega_c)$ Riesz map does (15)
become the represented gradient $\alpha u+p|_{\omega_c}$.

For a vector state, replace $H^1(\Omega)$ by $H^1(\Omega)^m$ and give $K$,
$b$, and $c$ compatible tensor shapes.  The residual contract itself is
unchanged.

### 3.2 Neumann and Dirichlet boundary control are different substitutions

For a Neumann control on $\Gamma_c \subset \Gamma_N$, replace the volume
coupling in (11) by

$$
  B_N(u,v)=\int_{\Gamma_c}u\,\operatorname{tr}v,
  \qquad U=L^2(\Gamma_c). \tag{16}
$$

The residual contains $-B_N(u,v)$.  The model must declare the trace needed
for (16), for example $V \rightarrow L^2(\Gamma_c)$ in the selected setting.  Its
transpose is the trace-adjoint coupling, and (15) becomes

$$
  j'(u)\delta u=\alpha(u,\delta u)_{L^2(\Gamma_c)}+
                  (\operatorname{tr}p,\delta u)_{L^2(\Gamma_c)}.
$$

For Dirichlet control, take a boundary control space such as
$U=H^{1/2}(\Gamma_c)$ and a declared lifting

$$
  L_D:U\longrightarrow H^1(\Omega),
  \qquad \operatorname{tr}_{\Gamma_c}L_Du=u.
$$

With any fixed Dirichlet datum absorbed in $\ell_0$, reconstruct

$$
  y=\widehat y+\ell_0+L_Du,\qquad \widehat y\in V. \tag{17}
$$

The control enters the residual through $a(L_Du,v)$, not through a boundary
load.  Consequently

$$
  D_uE(\widehat y,u)\delta u:v\longmapsto a(L_D\delta u,v),
$$

and the reduced derivative contains the chain-rule term

$$
  j'(u)\delta u=D_uJ\,\delta u-a(L_D\delta u,p). \tag{18}
$$

Any conversion of (18) into a boundary expression requires additional
analysis and a particular lifting/trace realization; it is not a generic
source-control rule.

### 3.3 Mixed formulation and Petrov--Galerkin spaces

The framework does not require a state to be an $H^1$ field.  For instance,
one mixed form of a diffusion equation introduces flux $q$ and scalar state
$y$:

$$
  q+K\nabla y=0,\qquad \nabla\!\cdot q=f+u.
$$

One possible substitution is

$$
  X_q=H(\operatorname{div};\Omega),\quad X_y=L^2(\Omega),
  \qquad Z_1=H(\operatorname{div};\Omega),\quad Z_2=L^2(\Omega).
$$

For test functions $(r,v) \in Z_1 \times Z_2$, residual blocks may be

$$
\begin{aligned}
 \langle E_1(q,y),r\rangle&=(K^{-1}q,r)_\Omega-(y,\nabla\!\cdot r)_\Omega
                              +\text{declared boundary contribution},\\
 \langle E_2(q,y,u),v\rangle&=(\nabla\!\cdot q,v)_\Omega-(f+u,v)_\Omega.
\end{aligned} \tag{19}
$$

The boundary contribution in the first line depends on the chosen boundary
condition and formulation; it must be made explicit rather than inferred.

At compilation, select trial spaces $X_{q,h}, X_{y,h}$ and test spaces
$Z_{1,h}, Z_{2,h}$.  Galerkin is the special case in which corresponding
trial and test spaces coincide.  Petrov--Galerkin uses different spaces or a
test transformation, while preserving the same typed residual form:

$$
  E_h:X_h\longrightarrow Z_h^*,\qquad
  \langle E_h(x_h),z_h\rangle=0\quad\forall z_h\in Z_h. \tag{20}
$$

Stabilisation and upwinding are therefore residual terms and/or a declared
test-space policy, not hidden solver options.

### 3.4 Boundary, flux, point, and energy observations

The same state equation supports different objectives by changing only (3).

| Observation | Map and observation space | Required declaration |
|---|---|---|
| Distributed tracking | $O(y)=y\|_{\omega_o}$, $Q=L^2(\omega_o)$ | restriction and $L^2$ pairing |
| Boundary tracking | $O(y)=tr_{\Gamma_o} y$, $Q=L^2(\Gamma_o)$ or a trace space | available trace and chosen boundary metric |
| Flux tracking | $O(y)=(K \nabla y) \cdot n$, $Q$ a declared boundary/dual space | normal-flux regularity and realization |
| Point sensors | $O(y)=(y(x_1),...,y(x_m))$, $Q=R^m$ | sufficient regularity, or an explicit discrete-only sensor policy |
| Energy tracking | $O(y)=y$ with $\Phi(y)=\tfrac12 a(y-y_d,y-y_d)$ | declared bilinear form and positivity/semidefiniteness policy |

For example, boundary tracking adds

$$
  \tfrac12\|\operatorname{tr}_{\Gamma_o}y-y_{d,\Gamma}\|_{L^2(\Gamma_o)}^2
$$

to (12).  Its derivative is a trace-adjoint contribution to the adjoint
right-hand side; it does not alter the state residual.  A point observation
is not automatically meaningful for $y \in H^1(\Omega)$ in all dimensions.
The program records the author's supplied regularity rationale or accepts a
clearly marked discrete-only observation; it cannot establish either.

### 3.5 Parabolic control: time is part of the space and residual

Consider

$$
  \partial_t y-\nabla\!\cdot(K\nabla y)+c y=f+u
  \quad\hbox{in }I\times\Omega,
  \qquad y(0)=y_0, \tag{21}
$$

with the spatial boundary conditions treated as in Section 3.1.  Let

$$
  V=H^1_{\Gamma_D}(\Omega),\qquad H=L^2(\Omega),
  \qquad
  Y=L^2(I;V)\cap H^1(I;V^*),\qquad Z=L^2(I;V). \tag{22}
$$

One formulation retains the time derivative on the state.  With $a_t$ the
spatial weak form and $B_t$ the control coupling, it declares

$$
  \langle E(y,u),v\rangle=
  \int_0^T\bigl[
    \langle \dot y(t),v(t)\rangle_{V^*,V}+a_t(y(t),v(t))
    -F_t(v(t))-B_t(u(t),v(t))
  \bigr]dt, \tag{23}
$$

together with the initial-trace condition $y(0)=y_0$, represented either by
the affine state space or by a separate residual/constraint block.  Typical
space-time tracking is

$$
  \tfrac12\int_0^T\|C y(t)-y_d(t)\|_Q^2dt
  +\tfrac\alpha2\int_0^T\|u(t)\|_U^2dt,
$$

possibly with a terminal cost.

An integration-by-parts-in-time formulation is also possible, but is a
different declared residual.  For sufficiently regular $v$,

$$
  \int_0^T\langle\dot y,v\rangle dt
  =(y(T),v(T))_H-(y(0),v(0))_H-
    \int_0^T(y,\dot v)_Hdt. \tag{24}
$$

Equation (24) shows exactly why initial and terminal terms cannot be guessed
by an implementation.  The backward adjoint and its terminal condition arise
from the transpose of the chosen time residual and the objective; they are
not a special case coded into an optimiser.

At the discrete level, (23) may be lowered to a space--time method or to a
chosen time-stepping residual.  These are distinct $E_h$ objects.  The
discrete adjoint must be the transpose of the actual compiled residual and
objective; the framework must not assume that discretising a displayed
continuous adjoint gives the same result.

### 3.6 Hyperbolic dynamics by block residuals

For the wave-type strong equation

$$
  \ddot y+A y=f+B u,\qquad y(0)=y_0,\quad\dot y(0)=v_0, \tag{25}
$$

introduce a velocity $v=\dot y$ if a first-order form is preferred.  The
unknown is then $x=(y,v,u)$ and the two equation blocks are, schematically,

$$
  E_1=\dot y-v,\qquad E_2=\dot v+A y-f-Bu. \tag{26}
$$

Each block has its own time-space test space and initial trace.  The spatial
operator $A$ may itself contain vector elasticity, flux, boundary, or mixed
terms.  This is not a new category of optimisation problem: it is a larger
block residual with temporal transformations and endpoint data.  Choice of
energy space, characteristic/inflow boundary treatment, damping, and time
integrator remains a declared formulation and policy.

### 3.7 Coefficient identification

Let a scalar coefficient $m$ be an optimisation variable in

$$
  -\nabla\!\cdot(m\nabla y)=f,
  \qquad m\in M_{\rm ad}. \tag{27}
$$

For example, choose $Y=Z=H^1_{\Gamma_D}(\Omega)$ and a parameter space
$M$ with the regularity and bounds needed for the product.  The residual is

$$
  \langle E(y,m),v\rangle=\int_\Omega m\nabla y\cdot\nabla v-\int_\Omega fv.
  \tag{28}
$$

The parameter derivative is a separate local contribution:

$$
  \langle D_mE(y,m)\,\delta m,v\rangle
  =\int_\Omega\delta m\,\nabla y\cdot\nabla v. \tag{29}
$$

Thus the problem is nonlinear in the full variable $(y,m)$, even though the
state equation is linear when $m$ is fixed.  A parameterisation such as
$m=\exp(q)$ is a transformation $T(q)$ and is handled by (7), including its
chain rule.  Positivity, identifiability, and product regularity are user
analysis obligations that must be declared rather than presumed.

## 4. Continuous semantics versus a discrete executable problem

The continuous semantic model records $X_i$, $Z_a$, residual terms,
transformations, observations, losses, pairings, constraints, and declared
requirements.  It contains no deal.II vectors, finite elements, DoF
handlers, or `AffineConstraints`.

A discretisation policy subsequently chooses, for every relevant space,

$$
  X_{i,h}\quad\hbox{and}\quad Z_{a,h},
$$

plus finite elements, meshes, component layouts, quadrature, trace
realizations, lifting construction, and assembly or matrix-free execution.
It lowers the semantic residual to

$$
  E_h:X_h\times D_h\longrightarrow Z_h^*, \tag{30}
$$

and similarly lowers observations, metrics, constraints, and transformations.
Dirichlet restrictions, periodicity, and hanging nodes are then realized by
the appropriate affine constraints and/or discrete liftings.

The discrete model must expose the value, Jacobian action, and transpose
action of *its own* residual,

$$
  E_h'(x_h)\,\delta x_h,\qquad E_h'(x_h)^*p_h,
$$

along with $J_h$ and its derivative.  This is the basis for exact discrete
adjoint-consistency and Taylor tests.  It makes no blanket claim that a
continuous adjoint, discretised independently, is identical to the discrete
transpose.

## 5. What the program may reasonably require as input

For each model, the public semantic input should be able to state the
following.

1. **Geometry and named regions:** domain, subdomains, boundary pieces, and
   time interval where applicable.
2. **Field and space declarations:** role, scalar/vector/tensor shape,
   continuous trial space, test space, dual pairing, and time structure.
3. **Variables, data, and admissible sets:** which symbols are unknown versus
   immutable data; box/affine/gauge constraints where relevant.
4. **Transformations:** liftings, restrictions, traces, parameterisations,
   and their required chain-rule actions.
5. **Residual blocks and terms:** a tested weak/mixed/very-weak action with
   source and target spaces, plus residual/Jacobian/transpose actions.
6. **Boundary and endpoint semantics:** essential restrictions or liftings;
   natural and Robin residual terms; initial/terminal traces; explicit
   nullspace policy if needed.
7. **Observations and objectives:** their target spaces, losses, target data,
   and derivative actions.
8. **Metrics and algorithmic geometry:** a pairing and any required inverse
   Riesz/duality map, rather than an unqualified norm name.
9. **Requirements and explicit policies:** traces, products, regularity,
   point-observation status, fractional-norm realization, compatibility, and
   any discrete-only semantics.

What it must *not* require is a proof of well-posedness, an automatically
derived weak form, a universally canonical finite-element space, or a PDE-
specific solver branch.  Conversely, what it must *not* accept as sufficient
input is merely a strong differential expression, a boundary-control flag,
or labels such as `Hminus1` and `point tracking` without their spaces,
pairings, and policies.

## 6. Consequences for the first implementation slice

The first executable slice can remain much smaller than the semantics above:
scalar stationary diffusion--reaction; primal conforming Galerkin spaces;
volume and Neumann controls; homogeneous or fixed lifted Dirichlet data;
distributed/boundary observations; $L^2$ and $H^1$ metrics; and box constraints.

Its interfaces should nevertheless already carry separate trial/test spaces,
dual residuals, transformations, observations, and metrics.  Then mixed
systems, Petrov--Galerkin methods, time dependence, Dirichlet control,
hyperbolic dynamics, and coefficient identification become additions of
terms, spaces, transformations, and compilation policies—not redesigns of
the optimisation core.
