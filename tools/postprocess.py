#!/usr/bin/env python3
"""Process native runner fields with a selected post-processing profile."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from nmopt_postprocess import PostprocessError
from nmopt_postprocess.chapter6 import (
    CHAPTER6_PROFILE,
    load_json_profile,
    with_parameter_configuration,
    with_comparison_plan,
)
from nmopt_postprocess.parameters import load_postprocess_configuration
from nmopt_postprocess.pipeline import (
    PostprocessProfile,
    process_artifact,
    process_input_root,
)
from nmopt_postprocess.render import (
    OutputFormat,
    OutputFormats,
)


PROFILES: dict[str, PostprocessProfile] = {
    "chapter6": CHAPTER6_PROFILE,
}
OUTPUT_FORMATS: tuple[OutputFormat, ...] = ("png", "svg")


def build_parser(
    default_profile: str | None = None, description: str | None = None
) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=description or __doc__)
    profile = parser.add_mutually_exclusive_group(required=False)
    if default_profile is None:
        profile.add_argument(
            "--profile",
            choices=tuple(PROFILES),
            help="built-in application profile",
        )
        profile.add_argument(
            "--profile-file",
            type=Path,
            help="JSON plotting profile (normally copied into a run root)",
        )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--artifact",
        type=Path,
        help="one runner artifact directory containing artifact.kv/native",
    )
    source.add_argument(
        "--input",
        type=Path,
        help="run root containing one or more artifact.kv files",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help=(
            "derived output directory; defaults to <artifact>/postprocess for "
            "--artifact or <input>/postprocess for --input"
        ),
    )
    parser.add_argument(
        "--format",
        dest="output_formats",
        nargs="+",
        choices=OUTPUT_FORMATS,
        default=None,
        metavar="FORMAT",
        help="one or more plot formats; defaults to png",
    )
    return parser


def _profile_from_run_snapshot(
    source: Path, fallback: PostprocessProfile, use_snapshot_style: bool = True
) -> PostprocessProfile:
    """Assemble the effective profile from the run snapshot and its ``.prm``."""

    root = source if source.is_dir() else source.parent
    for candidate_root in (root, *root.parents):
        snapshot = candidate_root / "plotting-profile.json"
        parameter_file = candidate_root / "parameters.prm"
        manifest_path = candidate_root / "run-manifest.json"
        if not parameter_file.exists():
            if not snapshot.exists() or not manifest_path.exists():
                continue
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            comparison = manifest.get("plotting", {}).get("resolved_comparison", {})
            return with_comparison_plan(
                fallback,
                str(comparison.get("rows", "")),
                str(comparison.get("columns", "")),
                str(comparison.get("group_by", "")),
            )
        if use_snapshot_style and snapshot.exists():
            fallback = load_json_profile(snapshot)
        configuration = load_postprocess_configuration(parameter_file)
        return with_parameter_configuration(fallback, configuration)
    return fallback


def _effective_output_formats(
    requested: list[OutputFormat] | None, profile: PostprocessProfile
) -> OutputFormats:
    if requested is not None:
        return tuple(requested)
    return profile.output_formats


def main(
    default_profile: str | None = None, description: str | None = None
) -> int:
    parser = build_parser(default_profile, description)
    arguments = parser.parse_args()
    profile_name = default_profile or getattr(arguments, "profile", None)
    explicit_profile = (
        default_profile is not None
        or getattr(arguments, "profile", None) is not None
        or getattr(arguments, "profile_file", None) is not None
    )
    if getattr(arguments, "profile_file", None) is not None:
        profile = load_json_profile(arguments.profile_file.resolve())
    else:
        profile = PROFILES[profile_name or "chapter6"]
    try:
        if arguments.artifact is not None:
            artifact = arguments.artifact.resolve()
            profile = _profile_from_run_snapshot(
                artifact, profile, use_snapshot_style=not explicit_profile
            )
            output_formats = _effective_output_formats(
                arguments.output_formats, profile
            )
            output = (
                arguments.output.resolve()
                if arguments.output is not None
                else artifact / "postprocess"
            )
            manifest = process_artifact(
                artifact,
                output,
                profile,
                output_formats=output_formats,
            )
            print(f"wrote {len(manifest['plots'])} field plot sets to {output}")
            return 0

        input_root = arguments.input.resolve()
        profile = _profile_from_run_snapshot(
            input_root, profile, use_snapshot_style=not explicit_profile
        )
        output_formats = _effective_output_formats(
            arguments.output_formats, profile
        )
        output_root = (
            arguments.output.resolve()
            if arguments.output is not None
            else input_root / "postprocess"
        )
        index = process_input_root(
            input_root,
            output_root,
            profile,
            output_formats=output_formats,
        )
        print(
            f"processed {index['success_count']}/{index['artifact_count']} "
            f"artifacts to {output_root}"
        )
        return (
            0
            if index["failure_count"] == 0
            and not index.get("comparison_errors")
            else 1
        )
    except (OSError, ValueError, PostprocessError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    raise SystemExit(main())
