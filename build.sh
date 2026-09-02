#!/usr/bin/env bash

set -euo pipefail

# The atomic actions below stay close to the underlying CMake commands.  The
# pipeline actions add timing and presentation without changing what
# configure, build, and test mean on their own.

usage() {
  cat <<'EOF'
Usage: ./build.sh [COMMAND] [PROFILE ...] [OPTIONS]

Atomic commands:
  configure [PROFILE ...]  Run cmake --preset for each profile.
  build [PROFILE ...]      Run cmake --build for each profile.
  test [PROFILE ...]       Run ctest --preset for each profile.
  pipeline [PROFILE ...]   Configure, time the build, and test profiles.
  all                      Run the pipeline for every supported profile.

Configuration commands:
  init-config              Create the machine-local build configuration.
  show-config              Display the active machine-local configuration.
  list                     List the supported profiles.

Profiles default to debug-neutral when omitted.  The no-argument command is
equivalent to pipeline debug-neutral.  Build-producing commands require an
existing build.local.conf; create it explicitly with init-config.

Options:
  --config FILE             Use FILE instead of build.local.conf.
  --cmake-arg ARG           Add ARG to configure commands; repeatable.
  --build-arg ARG           Add ARG to cmake --build; repeatable.
  --ctest-arg ARG           Add ARG to ctest; repeatable.
  --clean-first             Clean before building.
  --jobs N                  Request at most N jobs for each profile.
  --target TARGET           Build TARGET; only valid for one profile.
  --regex REGEX             Run tests whose names match REGEX.
  --label LABEL             Run tests with LABEL.
  --keep-going              Continue with later profiles after a failure.
  --dry-run                 Print commands without executing them.
  --verbose                 Print commands while executing them.
  --help, -h                Show this help text.

Examples:
  ./build.sh init-config
  ./build.sh configure debug-dealii
  ./build.sh build debug-neutral sanitize-neutral --jobs 4
  ./build.sh test debug-dealii --label dealii
  ./build.sh pipeline debug-neutral sanitize-neutral --clean-first
  ./build.sh all --keep-going --dry-run
EOF
}

die() {
  printf 'error: %s\n' "$1" >&2
  exit 1
}

script_directory="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repository_root="$script_directory"
cd -- "$repository_root"

config_file="$repository_root/build.local.conf"

declare -a supported_profiles=(
  debug-neutral
  sanitize-neutral
  debug-dealii
  release-dealii
)

declare -a cmake_args=()
declare -a build_args=()
declare -a ctest_args=()
declare -a selected_profiles=()

clean_first=false
requested_jobs=''
target_name=''
test_regex=''
test_label=''
keep_going=false
dry_run=false
verbose=false

declare -A profile_max_jobs=()
config_dealii_dir=''

print_profiles() {
  cat <<'EOF'
Supported profiles:
  debug-neutral    Fast backend-neutral Debug build.
  sanitize-neutral Backend-neutral AddressSanitizer/UBSan build.
  debug-dealii     Debug build with the deal.II backend.
  release-dealii   Optimized build with the deal.II backend.
EOF
}

is_supported_profile() {
  local candidate="$1"
  local profile

  for profile in "${supported_profiles[@]}"; do
    if [[ "$candidate" == "$profile" ]]; then
      return 0
    fi
  done
  return 1
}

validate_profiles() {
  local profile

  if ((${#selected_profiles[@]} == 0)); then
    selected_profiles=(debug-neutral)
  fi

  for profile in "${selected_profiles[@]}"; do
    is_supported_profile "$profile" || \
      die "unsupported profile '$profile' (use --help or list)"
  done
}

require_value() {
  local option="$1"
  (($# >= 2)) || die "$option needs a value"
}

parse_config_option() {
  local option="$1"

  require_value "$option" "$@"
  [[ -n "$2" ]] || die "$option needs a non-empty path"
  config_file="$2"
  if [[ "$config_file" != /* ]]; then
    config_file="$repository_root/$config_file"
  fi
}

parse_config_only_arguments() {
  while (($# > 0)); do
    case "$1" in
      --config)
        parse_config_option "$@"
        shift 2
        ;;
      --config=*)
        config_file="${1#*=}"
        [[ -n "$config_file" ]] || die "--config needs a non-empty path"
        if [[ "$config_file" != /* ]]; then
          config_file="$repository_root/$config_file"
        fi
        shift
        ;;
      --help|-h)
        usage
        exit 0
        ;;
      *)
        die "unexpected argument '$1'"
        ;;
    esac
  done
}

detect_processor_count() {
  local count=''

  if command -v nproc >/dev/null 2>&1; then
    count="$(nproc 2>/dev/null || true)"
  fi

  if [[ ! "$count" =~ ^[1-9][0-9]*$ ]] && \
     command -v getconf >/dev/null 2>&1; then
    count="$(getconf _NPROCESSORS_ONLN 2>/dev/null || true)"
  fi

  [[ "$count" =~ ^[1-9][0-9]*$ ]] || \
    die "could not determine the available processor count"

  printf '%s\n' "$count"
}

cache_dealii_dir() {
  local cache_file line candidate

  for cache_file in \
    "$repository_root/build/debug-dealii/CMakeCache.txt" \
    "$repository_root/build/release-dealii/CMakeCache.txt"; do
    [[ -r "$cache_file" ]] || continue

    while IFS= read -r line; do
      if [[ "$line" == deal.II_DIR:* ]]; then
        candidate="${line#*=}"
        if [[ -f "$candidate/deal.IIConfig.cmake" ]]; then
          printf '%s\n' "$candidate"
          return 0
        fi
      fi
    done < "$cache_file"
  done

  return 1
}

infer_dealii_dir() {
  local candidate=''

  for candidate in "${NMOPT_DEAL_II_DIR:-}" "${DEAL_II_DIR:-}"; do
    if [[ -n "$candidate" ]] && \
       [[ -f "$candidate/deal.IIConfig.cmake" ]]; then
      printf '%s\n' "$candidate"
      return 0
    fi
  done

  cache_dealii_dir || true
}

init_config() {
  local parent_directory processor_count dealii_dir

  [[ ! -e "$config_file" ]] || \
    die "configuration already exists: $config_file (edit it instead)"

  parent_directory="$(dirname -- "$config_file")"
  [[ -d "$parent_directory" ]] || \
    die "configuration directory does not exist: $parent_directory"

  processor_count="$(detect_processor_count)"
  dealii_dir="$(infer_dealii_dir)"

  # The generated file is deliberately a small trusted Bash fragment.  It
  # keeps paths readable and lets users edit the two machine-specific values
  # without touching the project presets or this script.
  {
    printf '%s\n' '# Machine-local nmopt build configuration.'
    printf '%s\n' '# This file is created by ./build.sh init-config.'
    printf '%s\n' '#'
    printf '%s\n' '# NMOPT_MAX_JOBS sets an independent maximum for each profile.'
    printf '%s\n' '# Deal.II profiles are often more memory-intensive than neutral ones.'
    printf '%s\n' '# Lower their entries independently if a build exhausts memory.'
    printf '%s\n' 'declare -A NMOPT_MAX_JOBS=('
    for profile in "${supported_profiles[@]}"; do
      printf '  [%q]=%q\n' "$profile" "$processor_count"
    done
    printf '%s\n' ')'
    printf '%s\n' '#'
    printf '%s\n' '# Set NMOPT_DEAL_II_DIR to the directory containing'
    printf '%s\n' '# deal.IIConfig.cmake when normal CMake discovery cannot find it.'
    printf 'NMOPT_DEAL_II_DIR=%q\n' "$dealii_dir"
  } > "$config_file"

  printf 'created build configuration: %s\n' "$config_file"
  printf '%s\n' 'lower the affected profile entry in NMOPT_MAX_JOBS if a build exhausts memory'
  if [[ -n "$dealii_dir" ]]; then
    printf 'inferred deal.II directory: %s\n' "$dealii_dir"
  else
    printf '%s\n' 'deal.II directory was not inferred; CMake will use normal discovery'
  fi
}

load_config() {
  [[ -f "$config_file" ]] || {
    printf 'error: build configuration does not exist: %s\n' "$config_file" >&2
    printf '%s\n' 'run ./build.sh init-config before a build command' >&2
    exit 1
  }

  unset NMOPT_MAX_JOBS NMOPT_DEAL_II_DIR

  # The local file is user-owned and generated in Bash assignment syntax.
  # Only the two allow-listed variables below are consumed after loading it.
  # shellcheck disable=SC1090
  source "$config_file"

  local declaration profile maximum
  declaration="$(declare -p NMOPT_MAX_JOBS 2>/dev/null || true)"
  [[ "$declaration" == 'declare -A NMOPT_MAX_JOBS='* ]] || \
    die "NMOPT_MAX_JOBS must be an associative array in $config_file; run init-config and edit the generated profile entries"

  profile_max_jobs=()
  for profile in "${supported_profiles[@]}"; do
    maximum="${NMOPT_MAX_JOBS[$profile]:-}"
    [[ "$maximum" =~ ^[1-9][0-9]*$ ]] || \
      die "NMOPT_MAX_JOBS[$profile] must be a positive integer in $config_file"
    profile_max_jobs["$profile"]="$maximum"
  done
  config_dealii_dir="${NMOPT_DEAL_II_DIR:-}"

  if [[ -n "$config_dealii_dir" ]] && \
     [[ ! -f "$config_dealii_dir/deal.IIConfig.cmake" ]]; then
    die "NMOPT_DEAL_II_DIR does not contain deal.IIConfig.cmake: $config_dealii_dir"
  fi
}

show_config() {
  load_config
  printf 'configuration: %s\n' "$config_file"
  printf 'maximum jobs by profile:\n'
  for profile in "${supported_profiles[@]}"; do
    printf '  %-18s %s\n' "$profile" "${profile_max_jobs[$profile]}"
  done
  if [[ -n "$config_dealii_dir" ]]; then
    printf 'deal.II directory: %s\n' "$config_dealii_dir"
  else
    printf '%s\n' 'deal.II directory: CMake default discovery'
  fi
}

effective_jobs() {
  local profile="$1"
  local requested="$2"
  local maximum="${profile_max_jobs[$profile]}"

  if ((requested > maximum)); then
    printf '%s\n' "$maximum"
  else
    printf '%s\n' "$requested"
  fi
}

print_command() {
  printf '+'
  printf ' %q' "$@"
  printf '\n'
}

run_command() {
  if [[ "$verbose" == true || "$dry_run" == true ]]; then
    print_command "$@"
  fi

  # A dry run must not configure, compile, test, or create any files.
  [[ "$dry_run" == true ]] || "$@"
}

run_configure() {
  local profile="$1"
  local pipeline_mode="${2:-false}"
  local -a command=(cmake --preset "$profile")
  local argument

  # The full pipeline is intentionally compact.  --verbose restores CMake's
  # normal STATUS messages, while an explicit --cmake-arg can override this
  # default for a particular invocation.
  if [[ "$pipeline_mode" == true && "$verbose" == false ]]; then
    command+=(--log-level=WARNING)
  fi
  if [[ "$profile" == *-dealii ]] && [[ -n "$config_dealii_dir" ]]; then
    command+=("-Ddeal.II_DIR=$config_dealii_dir")
  fi
  for argument in "${cmake_args[@]}"; do
    command+=("$argument")
  done

  run_command "${command[@]}"
}

run_build() {
  local profile="$1"
  local jobs requested
  local argument
  local -a command=(cmake --build --preset "$profile")

  requested="${requested_jobs:-${profile_max_jobs[$profile]}}"
  [[ "$requested" =~ ^[1-9][0-9]*$ ]] || \
    die "--jobs must be a positive integer"
  jobs="$(effective_jobs "$profile" "$requested")"

  command+=(--parallel "$jobs")
  if "$clean_first"; then
    command+=(--clean-first)
  fi
  if [[ -n "$target_name" ]]; then
    command+=(--target "$target_name")
  fi
  for argument in "${build_args[@]}"; do
    command+=("$argument")
  done

  printf '%s: using %s build jobs\n' "$profile" "$jobs"
  run_command "${command[@]}"
}

run_test() {
  local profile="$1"
  local progress="${2:-false}"
  local no_label_summary="${3:-false}"
  local argument
  local -a command=(ctest --preset "$profile")

  if [[ "$progress" == true ]]; then
    command+=(--progress)
  fi
  if [[ "$no_label_summary" == true ]]; then
    command+=(--no-label-summary)
  fi
  if [[ -n "$test_regex" ]]; then
    command+=(-R "$test_regex")
  fi
  if [[ -n "$test_label" ]]; then
    command+=(-L "$test_label")
  fi
  for argument in "${ctest_args[@]}"; do
    command+=("$argument")
  done

  run_command "${command[@]}"
}

print_rule() {
  printf '%s\n' '------------------------------------------------------------------------'
}

print_double_rule() {
  printf '%s\n' '========================================================================'
}

print_pipeline_phase() {
  local profile="$1"
  local phase="$2"

  printf '\n'
  print_rule
  printf '%s: %s\n' "$profile" "$phase"
  print_rule
}

format_elapsed() {
  local total_seconds="$1"
  local hours minutes seconds

  hours=$((total_seconds / 3600))
  minutes=$(((total_seconds % 3600) / 60))
  seconds=$((total_seconds % 60))

  # Keep short builds compact while making long builds readable at a glance.
  if ((hours > 0)); then
    printf '%dh %02dm %02ds' "$hours" "$minutes" "$seconds"
  elif ((minutes > 0)); then
    printf '%dm %02ds' "$minutes" "$seconds"
  else
    printf '%ds' "$seconds"
  fi
}

run_timed_build() {
  local profile="$1"
  local started elapsed elapsed_display status

  if [[ "$dry_run" == true ]]; then
    run_build "$profile"
    printf '%s: build timing skipped (dry run)\n' "$profile"
    return 0
  fi

  started="$SECONDS"
  if run_build "$profile"; then
    status=0
  else
    status=$?
  fi
  elapsed=$((SECONDS - started))
  elapsed_display="$(format_elapsed "$elapsed")"
  printf '%s: build elapsed %s\n' "$profile" "$elapsed_display"
  return "$status"
}

run_pipeline_profile() {
  local profile="$1"

  printf '\n\n'
  print_double_rule
  printf 'PROFILE: %s\n' "$profile"
  print_double_rule

  print_pipeline_phase "$profile" configure
  if ! run_configure "$profile" true; then
    return 1
  fi

  print_pipeline_phase "$profile" build
  if ! run_timed_build "$profile"; then
    return 1
  fi

  print_pipeline_phase "$profile" test
  if ! run_test "$profile" true true; then
    return 1
  fi

  printf '\n'
  print_double_rule
}

run_profiles() {
  local action="$1"
  local profile
  local -a failed_profiles=()

  for profile in "${selected_profiles[@]}"; do
    if "$action" "$profile"; then
      continue
    fi

    failed_profiles+=("$profile")
    if [[ "$keep_going" != true ]]; then
      return 1
    fi
  done

  if ((${#failed_profiles[@]} > 0)); then
    printf 'failed profiles:' >&2
    printf ' %s' "${failed_profiles[@]}" >&2
    printf '\n' >&2
    return 1
  fi

  return 0
}

parse_action_arguments() {
  while (($# > 0)); do
    case "$1" in
      --config|--config=*)
        if [[ "$1" == --config ]]; then
          parse_config_option "$@"
          shift 2
        else
          config_file="${1#*=}"
          [[ -n "$config_file" ]] || die "--config needs a non-empty path"
          if [[ "$config_file" != /* ]]; then
            config_file="$repository_root/$config_file"
          fi
          shift
        fi
        ;;
      --cmake-arg)
        require_value "$1" "$@"
        cmake_args+=("$2")
        shift 2
        ;;
      --cmake-arg=*)
        cmake_args+=("${1#*=}")
        shift
        ;;
      --build-arg)
        require_value "$1" "$@"
        build_args+=("$2")
        shift 2
        ;;
      --build-arg=*)
        build_args+=("${1#*=}")
        shift
        ;;
      --ctest-arg)
        require_value "$1" "$@"
        ctest_args+=("$2")
        shift 2
        ;;
      --ctest-arg=*)
        ctest_args+=("${1#*=}")
        shift
        ;;
      --clean-first)
        clean_first=true
        shift
        ;;
      --jobs)
        require_value "$1" "$@"
        requested_jobs="$2"
        shift 2
        ;;
      --jobs=*)
        requested_jobs="${1#*=}"
        shift
        ;;
      --target)
        require_value "$1" "$@"
        target_name="$2"
        shift 2
        ;;
      --target=*)
        target_name="${1#*=}"
        shift
        ;;
      --regex)
        require_value "$1" "$@"
        test_regex="$2"
        shift 2
        ;;
      --regex=*)
        test_regex="${1#*=}"
        shift
        ;;
      --label)
        require_value "$1" "$@"
        test_label="$2"
        shift 2
        ;;
      --label=*)
        test_label="${1#*=}"
        shift
        ;;
      --keep-going)
        keep_going=true
        shift
        ;;
      --dry-run)
        dry_run=true
        shift
        ;;
      --verbose)
        verbose=true
        shift
        ;;
      --help|-h)
        usage
        exit 0
        ;;
      -*)
        die "unknown option '$1' (use --help for usage)"
        ;;
      *)
        selected_profiles+=("$1")
        shift
        ;;
    esac
  done
}

if (($# == 0)); then
  command_name=pipeline
else
  command_name="$1"
  shift
fi

case "$command_name" in
  init-config)
    parse_config_only_arguments "$@"
    init_config
    ;;
  show-config)
    parse_config_only_arguments "$@"
    show_config
    ;;
  list)
    parse_config_only_arguments "$@"
    print_profiles
    ;;
  configure|build|test)
    parse_action_arguments "$@"
    validate_profiles

    case "$command_name" in
      configure)
        ((${#build_args[@]} == 0)) || \
          die "--build-arg is only valid with build"
        ((${#ctest_args[@]} == 0)) || \
          die "--ctest-arg is only valid with test"
        [[ "$clean_first" == false ]] || \
          die "--clean-first is only valid with build"
        [[ -z "$requested_jobs" ]] || \
          die "--jobs is only valid with build"
        [[ -z "$target_name" ]] || \
          die "--target is only valid with build"
        [[ -z "$test_regex" ]] || \
          die "--regex is only valid with test"
        [[ -z "$test_label" ]] || \
          die "--label is only valid with test"
        if [[ -f "$config_file" ]]; then
          load_config
        fi
        run_profiles run_configure
        ;;
      build)
        ((${#cmake_args[@]} == 0)) || \
          die "--cmake-arg is only valid with configure"
        ((${#ctest_args[@]} == 0)) || \
          die "--ctest-arg is only valid with test"
        ((${#test_regex} == 0)) || \
          die "--regex is only valid with test"
        ((${#test_label} == 0)) || \
          die "--label is only valid with test"
        if [[ -n "$target_name" ]] && ((${#selected_profiles[@]} != 1)); then
          die "--target requires exactly one profile"
        fi
        load_config
        run_profiles run_build
        ;;
      test)
        ((${#cmake_args[@]} == 0)) || \
          die "--cmake-arg is only valid with configure"
        ((${#build_args[@]} == 0)) || \
          die "--build-arg is only valid with build"
        [[ "$clean_first" == false ]] || \
          die "--clean-first is only valid with build"
        [[ -z "$requested_jobs" ]] || \
          die "--jobs is only valid with build"
        [[ -z "$target_name" ]] || \
          die "--target is only valid with build"
        run_profiles run_test
        ;;
    esac
    ;;
  pipeline)
    parse_action_arguments "$@"
    validate_profiles
    ((${#target_name} == 0)) || \
      die "--target is only valid with build"
    load_config
    run_profiles run_pipeline_profile
    ;;
  all)
    parse_action_arguments "$@"
    ((${#selected_profiles[@]} == 0)) || \
      die "all selects every profile; do not provide profile names"
    selected_profiles=("${supported_profiles[@]}")
    ((${#target_name} == 0)) || \
      die "--target is only valid with build"
    load_config
    run_profiles run_pipeline_profile
    ;;
  --help|-h)
    usage
    ;;
  *)
    die "unknown command '$command_name' (use --help for usage)"
    ;;
esac
