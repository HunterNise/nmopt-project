# Parameter and plotting profiles

## Status

Accepted repository and application-layer decision. The schema fixtures and
the first reusable loader/runner integration exist.

## Context

The headless Chapter 6 runner had constructed benchmark scenarios from C++
defaults and hardcoded matrices. The post-processing pipeline had read a
profile selected on the command line, while some benchmark-specific layout
policy remained inside the generic renderer. This made numerical experiments
harder to vary and made old plots dependent on current tooling defaults.

The project needs two reproducible configuration boundaries:

1. a reviewable description of the numerical experiment and its run matrix;
2. a reviewable description of how persisted fields are presented.

These boundaries must expose source omissions and development hypotheses
without moving PDE assembly or solver logic into a configuration language.

## Decision

### Storage

Versioned configuration inputs live below `parameters/`:

```text
parameters/
  chapter-6/
    b1/
      authoritative.prm
    b2/
      authoritative.prm
      development/
        <named-development-family>.prm
    plotting/
      chapter-6-b1.json
      chapter-6-b2.json
```

The `docs/benchmarks/` documents remain the human-readable frozen scientific
contracts. The files below `parameters/` are executable instantiations of
those contracts or explicitly marked development variations. Generated
evidence remains below `runs/` and is never a source of configuration.

### Experiment files

Numerical experiment families use Deal.II-style `.prm` files parsed through a
typed `ParameterHandler` boundary. The file owns every choice that affects
the numerical run or retained evidence:

- benchmark identity, recipe, source reference, and source revision;
- matrix axes and optional selection filters;
- problem and runtime data;
- registered forcing, target, fixed-data, and transport specifications;
- observation and boundary regions, IDs, orientations, realizations, and
  quadrature;
- geometry, mesh, finite-element degree, refinement, and provenance;
- execution product and compilation-session policy;
- solver method, initial control, stopping, line search, and iteration
  policies;
- run kind, build profile, output destination, harness measurements, and
  retained native fields; and
- the post-processing style reference and comparison-axis binding.

Function entries are declarative inputs to registered function constructors.
They are not arbitrary C++. A constructor may accept a constrained expression
record when its syntax and available functions are declared and validated by
the application capability boundary. B1 uses this mechanism for scalar
forcing definitions; other function ports remain bounded by their registered
constructors.

### Matrix expansion

Only entries under `Matrix` participate in expansion. If no `Selection`
section is present, the runner executes the full Cartesian product of the
declared axes. A concrete artifact is one resolved combination; a run set may
contain many artifacts.

For example, B1 declares four regularisation values and two methods, then
excludes the unreported steepest-descent $\beta=10^{-6}$ coordinate to retain
the seven source figure artifacts. B2 declares two observation regions and two
target profiles, yielding four artifacts. A development forcing family may
declare three forcing values with one fixed beta and yield three artifacts.

Selection is optional and only filters the declared matrix. It must use stable
axis and value IDs, not numeric positions. The resolver rejects unknown
values, empty products, duplicate combinations, and unsupported combinations
before execution.

### Plotting profiles

Plotting profiles use JSON with schema ID `nmopt-plot-v1`. They own reusable
presentation policy:

- persisted field names and logical display names;
- field titles and colorbar labels;
- colormap, normalization, interpolation, and mesh-overlay policy;
- axis labels and output formats; and
- default comparison-axis bindings and axis label/order metadata.

The profile is a style and field contract, not a fixed artifact matrix. The
`.prm` file may override the comparison rows, columns, and grouping for its
declared matrix. Therefore the B2 style can be reused by a three-forcing
development sweep whose comparison is one row by three forcing columns,
without reusing the authoritative 2x2 arrangement.

The generic renderer receives an explicit comparison plan. It must not infer
benchmark layout from the number of panels, and it must not contain B1/B2
branches.

### Precedence and smoke runs

The effective configuration is resolved in this order:

```text
parameter-file defaults
  → Selection section
  → explicit CLI selection filter
  → explicit CLI refinement override
```

The CLI refinement override is retained specifically for smoke runs. It may
change the realized mesh refinement without creating an ad hoc parameter file,
but it must update mesh provenance and be recorded as an explicit override in
the run manifest.

Output-directory relocation may remain a destination-level CLI concern.
Framework revision is provenance of the executable/environment rather than a
mathematical experiment value and may be supplied or discovered by the
runner. All other run-affecting choices belong in the `.prm` file.

### Provenance

Every run manifest must record:

- the parameter-file path and content hash;
- the declared matrix and any selection;
- every resolved artifact combination;
- the plotting-profile path and content hash;
- the resolved comparison plan; and
- all CLI overrides, including refinement.

The runner should copy the resolved parameter and plotting documents into the
run directory. Post-processing uses those snapshots by default. Applying a
different plotting profile is an explicit derived-output override and is
recorded in the post-processing manifest.

## Consequences

This decision provides one authoritative B1 or B2 file for the full benchmark
matrix while retaining per-artifact reproducibility. Development families can
vary missing source components without changing the frozen benchmark file.

The parameter loader resolves matrix values into the existing typed
scenario records; it must not bypass semantic validation or introduce
backend-specific branches into generic solvers. The post-processing loader
must resolve JSON into explicit profile and comparison-policy objects; it must
not reconstruct authoritative numerical values.

The initial schema fixtures are:

- [parameter-file reference](../reference/parameter-files.md);
- [B1 authoritative family](../../parameters/chapter-6/b1/authoritative.prm);
- [B2 authoritative family](../../parameters/chapter-6/b2/authoritative.prm);
- [B2 forcing development family](../../parameters/chapter-6/b2/development/forcing-sweep.prm);
- [B1 plotting profile](../../parameters/plotting/chapter-6-b1.json); and
- [B2 plotting profile](../../parameters/plotting/chapter-6-b2.json).

The current content hash is a labelled deterministic FNV-1a-64 digest used to
detect configuration drift; it is not an authentication checksum.
