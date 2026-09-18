# Executable parameter inputs

The files under `parameters/` are tracked, executable inputs for repository
applications and reproduction workflows. They are distinct from the
interpretive material under `docs/` and from generated evidence under `runs/`.

`authoritative.prm` identifies the accepted repository reproduction input for a
benchmark family. It does not claim that every value in the file is directly
stated by the source publication; values may also encode the repository's
selected numerical and execution policy.

Files under `development/` record explicit investigation or hypothesis
families. They are versioned so that their assumptions and results remain
reviewable, but they are not automatically reproduction evidence.

Plotting JSON files under `parameters/plotting/` are versioned presentation
policies. Generated parameter snapshots, run manifests, native fields, and
derived plots belong under `runs/`, not in this input tree.
