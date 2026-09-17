# Architecture maps

This page is the visual companion to the prose overviews.

The first map shows the widest project view. The later maps zoom into individual
parts of the same architecture rather than introducing alternative architectures.
Read them like a nested set of maps: start broad, then move inward only where you
need more detail.

The diagrams are for orientation rather than navigation. Mermaid rendering is not a
reliable place to encode repository links, so each map is followed by prose that
points to the relevant overview, reference, design, or source area.

## Visual legend

The maps use a small visual grammar consistently.

| Visual form | Meaning |
| --- | --- |
| Rectangle | An owned subsystem or application component |
| Rounded node | An adapter, formulation, or process-oriented service |
| Double-bordered node | A shared architectural boundary or product exposed downstream |
| Slanted node | Serialized/external input or persisted artifact |
| Solid arrow | Main construction, data, or execution flow |
| Dashed arrow | Orchestration, provenance, borrowing, or other side-channel relationship |
| Text on an arrow | The action connecting two nouns: validate, compile, adapt, expose, solve, record, and so on |

Where a node contains an example, the example appears after a blank line in smaller
italic text. It illustrates the category rather than defining it.

## 1. Master map

This is the detailed version of the simplified diagram in
[Project overview and architecture](project-architecture.md).

```mermaid
flowchart TB
  subgraph Start["User-facing starting points"]
    direction LR
    PF[/"Parameter files<br/><br/><small><i>checked .prm run configurations</i></small>"/]
    SC["Typed scenario<br/><br/><small><i>B1/B2-style problem, discretization, and solver selections</i></small>"]
    RC(["Recipe<br/><br/><small><i>reusable typed problem-family builder</i></small>"])
    DC(["Direct semantic construction<br/><br/><small><i>library or client code building ProblemSpec</i></small>"])
    EA["Existing PDE/OCP application<br/><br/><small><i>deal.II code with mesh, assembly, native solves, and output</i></small>"]

    PF -->|parse into| SC
    SC -->|selects and configures| RC
  end

  PS[["Semantic ProblemSpec<br/><br/><small><i>spaces · PDE terms · objective · metric · requested formulation</i></small>"]]
  RC -->|builds| PS
  DC -->|builds| PS

  RR(["Resolved semantic request<br/><br/><small><i>validated choices and explicit requirements</i></small>"])
  PS -->|validate + resolve| RR

  CR["Compiler-owned numerical realization<br/><br/><small><i>deal.II coordinates · operators · solves · metric · native view</i></small>"]
  RR -->|compile with runtime data + discretization session| CR

  AB(["Application-owned binding<br/><br/><small><i>model actions · state/adjoint solves · metric adapters</i></small>"])
  EA -->|adapt only required operations| AB

  CN[["Common numerical services<br/><br/><small><i>layouts · primal/covector values · formulation-specific actions · solve reports · metrics · constraints</i></small>"]]
  CR -->|exposes| CN
  AB -->|exposes| CN

  subgraph Products["Formulation products"]
    direction LR
    RED(["Reduced state–adjoint<br/><br/><small><i>state elimination + adjoint reduced derivative</i></small>"])
    OTD(["Supplied OTD<br/><br/><small><i>all-at-once optimality system</i></small>"])
    KKT(["Quadratic KKT<br/><br/><small><i>equality-constrained saddle-point product</i></small>"])
    PDAS(["Complementarity / PDAS<br/><br/><small><i>active-set formulation</i></small>"])
  end

  CN -->|compose as| RED
  CN -->|compose as| OTD
  CN -->|compose as| KKT
  CN -->|compose as| PDAS

  subgraph Algorithms["Algorithms and numerical solves"]
    direction LR
    RS["Reduced search algorithms<br/><br/><small><i>steepest descent · L-BFGS · Newton / trust-region</i></small>"]
    LS["Linear / Krylov solve services<br/><br/><small><i>KKT and supplied-OTD systems</i></small>"]
    AS["Active-set services<br/><br/><small><i>PDAS iterations + complementarity evidence</i></small>"]
  end

  RED -->|optimize with| RS
  OTD -->|solve with| LS
  KKT -->|solve with| LS
  PDAS -->|iterate with| AS

  OUT[["Results and numerical evidence<br/><br/><small><i>accepted state/control · reports · work records</i></small>"]]
  RS -->|returns| OUT
  LS -->|returns| OUT
  AS -->|returns| OUT

  EXP["Experiment and replication orchestration<br/><br/><small><i>nmopt_runner · run sets · manifests · artifact paths</i></small>"]
  ART[/"Persisted artifacts<br/><br/><small><i>native fields · traces · reports · comparisons</i></small>"/]

  SC -.->|benchmark identity + run policy| EXP
  OUT -.->|solver/formulation evidence| EXP
  CR -.->|native output access| EXP
  EA -.->|native output access when externally owned| EXP
  EXP -->|write + post-process| ART
```

The initial subgraph deliberately contains both structured authoring and an existing
application. They are equally valid starting points. Parameter files, scenarios, and
recipes are *outer authoring conveniences* for the semantic path; an existing
application reaches the same common services through an adapter instead.

The experiment and replication layer is shown with dashed incoming edges because it
orchestrates selected uses of the numerical system rather than defining the PDE or
formulation.

Read [Project overview and architecture](project-architecture.md) for the prose
interpretation of the whole map.

## 2. Structured authoring and compiler path

This map zooms into the left side of the master map.

```mermaid
flowchart TB
  subgraph Authoring["Typed problem authoring"]
    direction LR
    PF[/"Parameter file<br/><br/><small><i>serialized supported choices</i></small>"/]
    SC["Scenario<br/><br/><small><i>concrete data · discretization · formulation · solver choices</i></small>"]
    RC(["Recipe<br/><br/><small><i>reusable semantic problem-family builder</i></small>"])
    DC(["Direct semantic construction<br/><br/><small><i>client code</i></small>"])
    PF -->|parse into| SC
    SC -->|configures| RC
  end

  PS[["ProblemSpec<br/><br/><small><i>backend-neutral semantic graph</i></small>"]]
  RC -->|builds| PS
  DC -->|builds| PS

  SV(["Semantic validation<br/><br/><small><i>graph structure · roles · declared policies</i></small>"])
  RR[["Resolved request<br/><br/><small><i>semantically relevant choices closed</i></small>"]]
  PS -->|check| SV
  SV -->|resolve| RR

  RT[/"Runtime data bindings<br/><br/><small><i>forcing · target · coefficients · boundary data</i></small>"/]
  DS["Compilation session / discretization policy<br/><br/><small><i>mesh · FE choices · realization policy</i></small>"]
  CP(["Compiler planning + lowering"])
  RR -->|request| CP
  RT -->|supply concrete data| CP
  DS -->|supply realization context| CP

  NR["deal.II numerical realization<br/><br/><small><i>coordinates · operators · metrics · solves</i></small>"]
  CP -->|construct| NR

  PROD[["Compiled product<br/><br/><small><i>formulation services · native view · manifest</i></small>"]]
  NR -->|package| PROD
```

The recipe/scenario machinery is not the semantic language itself. A recipe is a
typed builder for a reusable problem family; a scenario freezes one concrete set of
problem and execution choices. Direct client code can build the same `ProblemSpec`
without either layer.

For prose, read [Describing and compiling a problem](semantic-compiler.md). The
detailed current capability ledger is
[`docs/implementation/v1/semantic-compiler.md`](../implementation/v1/semantic-compiler.md).

## 3. Existing-application integration path

This map zooms into the other producer path.

```mermaid
flowchart TB
  APP["Existing numerical application<br/><br/><small><i>mesh · FE spaces · assembled operators · native solves · output</i></small>"]

  PM["Application OCP mathematics<br/><br/><small><i>control coordinates · objective · derivative actions · metric meaning</i></small>"]
  APP -.->|provides numerical operators to| PM

  subgraph Bind["Thin nmopt binding"]
    direction LR
    MOD(["Executable model adapters"])
    SOL(["State / adjoint solve adapters"])
    MET(["Metric / constraint adapters"])
    LAY(["Layouts + state/control partition"])
  end

  APP -->|native callbacks and borrowed objects| MOD
  APP -->|native solve services| SOL
  PM -->|mathematical roles| MOD
  PM -->|optimization geometry| MET
  PM -->|coordinate structure| LAY

  CN[["Common numerical services"]]
  MOD -->|exposes| CN
  SOL -->|exposes| CN
  MET -->|exposes| CN
  LAY -->|labels| CN

  FORM(["Selected formulation"])
  CN -->|compose as| FORM

  ALG["Algorithm / formulation solver"]
  FORM -->|consume with| ALG

  RES[["Result"]]
  ALG -->|returns| RES
  RES -.->|reconstruct/write through native application| APP
```

The important ownership rule is visible in the dashed relationships: the binding can
borrow application objects, but it does not become their owner. Output returns to the
native application because physical reconstruction may require its mesh, coordinate
maps, and writers.

For prose, read
[Integrating an existing PDE application](external-applications.md). The worked
deal.II example remains beside its source at
[`apps/external-dealii/step-4/external-integration-overview.md`](../../apps/external-dealii/step-4/external-integration-overview.md).

## 4. Shared foundations and formulation-specific products

This map zooms into the convergence point. The formulations share typed numerical
foundations, but they do not all consume one identical service bundle.

```mermaid
flowchart TB
  subgraph Base["Shared numerical foundations"]
    direction LR
    L["Layouts + primal/covector values"]
    P["Dual pairings + solve reports + lifetime evidence"]
    M["Metrics / constraints where required"]
  end

  subgraph ReducedInputs["Reduced-formulation services"]
    direction LR
    E["Executable model actions<br/><br/><small><i>residual · JVP · VJP · objective · derivative</i></small>"]
    S["State / adjoint solve services"]
    H["Optional reduced Hessian action"]
  end

  subgraph OTDInputs["Supplied-OTD services"]
    direction LR
    OACT["Supplied optimality-system actions<br/><br/><small><i>residual · JVP · VJP</i></small>"]
    OSOL["Supplied system solve action"]
  end

  subgraph KKTInputs["Quadratic-KKT services"]
    direction LR
    QACT["Quadratic / equality / transpose actions<br/><br/><small><i>Q · D · D-transpose · KKT transpose</i></small>"]
    KASS["KKT assumptions + block pairings<br/><br/><small><i>rank · kernel positivity · symmetry evidence</i></small>"]
  end

  RED(["Reduced state–adjoint formulation"])
  OTD(["Supplied OTD formulation"])
  KKT(["Quadratic KKT formulation"])
  PDAS(["Complementarity / PDAS formulation"])

  L --> RED
  P --> RED
  M --> RED
  E --> RED
  S --> RED
  H -.->|Newton-type methods| RED

  L --> OTD
  P --> OTD
  OACT --> OTD
  OSOL --> OTD

  L --> KKT
  P --> KKT
  QACT --> KKT
  KASS --> KKT

  KKT -->|base KKT product| PDAS
  M --> PDAS
  PDIN["Box complementarity + bounds<br/><br/><small><i>classification policy · active/free sets</i></small>"] --> PDAS

  RS["Reduced search<br/><br/><small><i>steepest descent · L-BFGS · Newton / trust-region</i></small>"]
  KS["KKT / all-at-once linear solves"]
  AS["Active-set iteration"]

  RED -->|optimize| RS
  OTD -->|solve| KS
  KKT -->|solve| KS
  PDAS -->|iterate| AS

  ReducedInputs ~~~ Base
  KKTInputs ~~~ OTDInputs
```

The distinction is deliberate. `ExecutableModelT` plus state/adjoint solves is the
current reduced service boundary. A supplied-OTD system instead owns its own
three-block residual/JVP/VJP/solve actions. The quadratic KKT product owns explicit
quadratic, equality, transpose, pairing, and assumption data. PDAS is built around a
KKT product plus box complementarity, metric/bound data, and active-set policy.

The commonality is therefore **typed numerical roles and composition rules**, not a
single universal formulation interface.

For the reduced path, continue with
[Reduced state–adjoint optimization](reduced-optimization.md). For exact
mathematical conventions, use
[Theoretical formalism](../design/theoretical-formalism.md).

## 5. Reduced evaluation lifecycle

This map zooms into one reduced evaluation and makes *derivative augmentation*
explicit.

```mermaid
flowchart TB
  U["Control u"]

  V(["Evaluate value"])
  SS["State solve"]
  J["Objective value"]
  RV[["Retained value<br/><br/><small><i>control + state + objective + solve evidence</i></small>"]]

  U -->|request value| V
  V -->|solve PDE| SS
  SS -->|retain state| RV
  V -->|evaluate objective| J
  J -->|store| RV

  TEST{"Need derivative here?"}
  RV --> TEST
  TEST -->|no: trial value is enough| DONE["Return / acceptance test"]

  AUG(["Augment derivative<br/><br/><small><i>extend retained value without another state solve</i></small>"])
  TEST -->|yes| AUG

  JD["Objective derivatives"]
  ADJ["Adjoint solve"]
  VJP["Residual transpose action"]
  RD[["Reduced covector"]]

  AUG -->|differentiate objective| JD
  AUG -->|solve using state derivative| ADJ
  ADJ -->|pull back through residual| VJP
  JD -->|control contribution| RD
  VJP -->|state-dependence contribution| RD

  MET(["Metric inverse / direction policy"])
  RD -->|convert or update| MET
```

The branch at **Need derivative here?** is the key runtime optimization. A rejected
line-search trial can stop after value evaluation. When a derivative is needed, the
retained state is augmented with adjoint-derived first-order information.

Read [Reduced state–adjoint optimization](reduced-optimization.md) for the equations
behind this map.

## 6. Experiment and evidence flow

This map zooms into the outer execution layer. It is intentionally downstream of the
numerical architecture.

```mermaid
flowchart TB
  BC["Benchmark contract<br/><br/><small><i>what to run · record · compare · accept</i></small>"]
  PF[/"Parameter files<br/><br/><small><i>serialized supported run choices</i></small>"/]
  SC["Typed scenario<br/><br/><small><i>problem + discretization + solver selection</i></small>"]

  PF -->|parse into| SC
  BC -.->|defines run/evidence requirements for| SC

  EX(["Execution adapter"])
  SC -->|configure| EX

  CORE[["Numerical formulation / solver result"]]
  EX -->|invoke framework| CORE

  MAN[/"Detached evidence<br/><br/><small><i>compilation manifest · solver report · environment</i></small>"/]
  NAT[/"Native numerical output<br/><br/><small><i>VTU fields · mesh views · traces</i></small>"/]

  CORE -->|record| MAN
  CORE -.->|reconstruct through native view/application| NAT

  RUN["Run controller / nmopt_runner<br/><br/><small><i>run sets · paths · failures · artifact identity</i></small>"]
  SC -.->|benchmark identity| RUN
  MAN -->|persist under| RUN
  NAT -->|persist under| RUN

  PP(["Post-processing / comparison"])
  RUN -->|read persisted artifacts| PP

  REP[/"Reports and derived plots"/]
  PP -->|write| REP
```

A benchmark contract evaluates framework behavior; it is not the capability registry
for the compiler or formulation layers. This is why unimplemented B3–B6 benchmark
campaigns do not imply absent KKT/PDAS capability.

Read [Experiments and replication](experiments-and-replication.md) for the prose
model and the distinction among recipes, scenarios, parameter files, benchmarks, and
artifacts.

## 7. Source-tree orientation

The maps above describe responsibilities. The source tree is only an approximate
physical projection of those responsibilities:

```text
include/nmopt/contract/       common numerical types and formulations
include/nmopt/solvers/        reduced optimization algorithms
include/nmopt/semantic/v1/    semantic problem language and validation
include/nmopt/compiler/v1/    planning, lowering, and compiled products
include/nmopt/dealii/         vector backend + reusable deal.II numerical services
include/nmopt/application/    recipes, scenarios, adapters, runner support
include/nmopt/experiment/     detached experiment/provenance records

apps/external-dealii/         direct-integration case studies
apps/nmopt-runner/            repository experiment runner
parameters/                   checked run configurations
tests/                        contract, compiler, numerical, and application evidence
```

Do not infer architecture solely from directory names. `application/` mixes generic
orchestration with project-specific scenario code, and `dealii/` contains both the
narrow vector backend and richer numerical realization services. The prose overviews
explain those distinctions where they matter.
