#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage: tools/run_chapter6.sh --benchmark b1|b2 [options]

Run an already-built Chapter 6 runner, then generate post-processing and a
manifest-aware report in the run set allocated by the runner.

Options:
  --benchmark BENCHMARK       required: b1 or b2
  --framework-revision REV    provenance value; defaults to current Git HEAD
  --run-kind KIND              defaults to development
  --refinement N               mesh override; defaults to 1
  --output DIRECTORY           run root; defaults to runs
  --runner PATH                runner executable; defaults to
                              build/debug-dealii/bin/nmopt_runner
  --python PATH                Python interpreter; defaults to python3
  --format FORMAT [FORMAT ...] post-processing formats: png and/or svg;
                              defaults to png
  --help                      show this message

Example:
  tools/run_chapter6.sh --benchmark b2 --refinement 1 --format png svg
EOF
}

die() {
  printf 'error: %s\n' "$1" >&2
  exit 1
}

script_directory="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repository_root="$(cd -- "$script_directory/.." && pwd)"
cd -- "$repository_root"

benchmark=''
framework_revision=''
run_kind='development'
refinement='1'
output_directory='runs'
runner_path='build/debug-dealii/bin/nmopt_runner'
python_command='python3'
output_formats=()

while (($# > 0)); do
  case "$1" in
    --benchmark)
      (($# >= 2)) || die "--benchmark needs a value"
      benchmark="$2"
      shift 2
      ;;
    --framework-revision)
      (($# >= 2)) || die "--framework-revision needs a value"
      framework_revision="$2"
      shift 2
      ;;
    --run-kind)
      (($# >= 2)) || die "--run-kind needs a value"
      run_kind="$2"
      shift 2
      ;;
    --refinement)
      (($# >= 2)) || die "--refinement needs a value"
      refinement="$2"
      shift 2
      ;;
    --output)
      (($# >= 2)) || die "--output needs a value"
      output_directory="$2"
      shift 2
      ;;
    --runner)
      (($# >= 2)) || die "--runner needs a value"
      runner_path="$2"
      shift 2
      ;;
    --python)
      (($# >= 2)) || die "--python needs a value"
      python_command="$2"
      shift 2
      ;;
    --format)
      shift
      (($# > 0)) || die "--format needs at least one value"
      while (($# > 0)) && [[ "$1" != -* ]]; do
        case "$1" in
          png|svg)
            output_formats+=("$1")
            ;;
          *)
            die "unsupported post-processing format '$1'"
            ;;
        esac
        shift
      done
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      die "unknown option '$1' (use --help for usage)"
      ;;
  esac
done

[[ "$benchmark" == b1 || "$benchmark" == b2 ]] || \
  die "--benchmark must be b1 or b2"

if [[ -z "$framework_revision" ]]; then
  framework_revision="$(git rev-parse HEAD)"
fi

if [[ "$runner_path" != */* ]]; then
  runner_path="$(command -v "$runner_path" || true)"
fi
[[ -x "$runner_path" ]] || die "runner is not executable: $runner_path"
command -v "$python_command" >/dev/null 2>&1 || \
  die "Python interpreter not found: $python_command"

runner_arguments=(
  --benchmark "$benchmark"
  --framework-revision "$framework_revision"
  --run-kind "$run_kind"
  --output "$output_directory"
)
if [[ -n "$refinement" ]]; then
  runner_arguments+=(--refinement "$refinement")
fi

runner_output_file="$(mktemp)"
cleanup() {
  rm -f -- "$runner_output_file"
}
trap cleanup EXIT

printf 'running %s\n' "$runner_path"
if ! "$runner_path" "${runner_arguments[@]}" | tee "$runner_output_file"; then
  die "benchmark runner failed"
fi

run_directory="$(awk '
  /^B[12] wrote / {
    path = $0
    sub(/^B[12] wrote /, "", path)
    sub(/\/artifacts\/.*\/artifact\.kv$/, "", path)
    print path
    exit
  }
' "$runner_output_file")"
[[ -n "$run_directory" ]] || \
  die "could not determine the allocated run directory from runner output"

postprocess_arguments=(
  tools/postprocess.py
  --profile chapter6
  --input "$run_directory"
  --output "$run_directory/postprocess"
)
if ((${#output_formats[@]} > 0)); then
  postprocess_arguments+=(--format "${output_formats[@]}")
fi

"$python_command" "${postprocess_arguments[@]}"
"$python_command" tools/chapter6_report.py \
  --run-manifest "$run_directory/run-manifest.json" \
  --output "$run_directory/report"

printf 'completed run set: %s\n' "$run_directory"
