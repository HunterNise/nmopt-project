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
    effective_profile_record,
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


def _find_run_snapshot(
    source: Path,
) -> tuple[Path, Path | None, Path | None, Path | None] | None:
    root = source if source.is_dir() else source.parent
    for candidate_root in (root, *root.parents):
        snapshot = candidate_root / "plotting-profile.json"
        parameter_file = candidate_root / "parameters.prm"
        manifest_path = candidate_root / "run-manifest.json"
        if parameter_file.exists() or (snapshot.exists() and manifest_path.exists()):
            return (
                candidate_root,
                parameter_file if parameter_file.exists() else None,
                snapshot if snapshot.exists() else None,
                manifest_path if manifest_path.exists() else None,
            )
    return None


def _fnv1a64(path: Path) -> str:
    digest = 1469598103934665603
    for byte in path.read_bytes():
        digest ^= byte
        digest = (digest * 1099511628211) & 0xFFFFFFFFFFFFFFFF
    return f"fnv1a64:{digest:016x}"


def _manifest_section(
    manifest: dict[str, object], section: str
) -> dict[str, object]:
    value = manifest.get(section, {})
    return value if isinstance(value, dict) else {}


def _postprocess_provenance(
    source: Path,
    profile: PostprocessProfile,
    output_formats: OutputFormats,
    *,
    explicit_profile: bool,
    profile_file: Path | None,
    profile_name: str | None,
    requested_formats: list[OutputFormat] | None,
) -> dict[str, object]:
    """Describe the source documents and effective choices used for plotting."""

    snapshot = _find_run_snapshot(source)
    manifest: dict[str, object] = {}
    parameter_snapshot: Path | None = None
    plotting_snapshot: Path | None = None
    manifest_path: Path | None = None
    if snapshot is not None:
        _, parameter_snapshot, plotting_snapshot, manifest_path = snapshot
        if manifest_path is not None:
            document = json.loads(manifest_path.read_text(encoding="utf-8"))
            if isinstance(document, dict):
                manifest = document

    parameters_manifest = _manifest_section(manifest, "parameters")
    parameters: dict[str, object] = {}
    if isinstance(parameters_manifest.get("file"), str):
        parameters["file"] = parameters_manifest["file"]
    elif parameter_snapshot is not None:
        parameters["file"] = str(parameter_snapshot)
    manifest_parameter_hash = parameters_manifest.get("content_hash")
    if parameter_snapshot is not None:
        parameters["content_hash"] = _fnv1a64(parameter_snapshot)
        if (
            isinstance(manifest_parameter_hash, str)
            and manifest_parameter_hash != parameters["content_hash"]
        ):
            parameters["run_manifest_content_hash"] = manifest_parameter_hash
    elif isinstance(manifest_parameter_hash, str):
        parameters["content_hash"] = manifest_parameter_hash
    if parameter_snapshot is not None:
        parameters["snapshot_file"] = str(parameter_snapshot)
    for key in ("selection", "declared_matrix", "resolved_combinations"):
        if key in parameters_manifest:
            parameters[key] = parameters_manifest[key]

    plotting_manifest = _manifest_section(manifest, "plotting")
    plotting: dict[str, object] = {}
    snapshot_profile_file = plotting_manifest.get("profile_file")
    snapshot_profile_hash = plotting_manifest.get("content_hash")
    if explicit_profile and profile_file is not None:
        plotting["source"] = "explicit-profile-file"
        plotting["profile_file"] = str(profile_file)
        plotting["content_hash"] = _fnv1a64(profile_file)
        if isinstance(snapshot_profile_file, str):
            plotting["run_snapshot"] = {
                "profile_file": snapshot_profile_file,
                "content_hash": snapshot_profile_hash,
            }
    elif explicit_profile:
        plotting["source"] = "explicit-built-in-profile"
        plotting["profile"] = profile_name or "chapter6"
        if isinstance(snapshot_profile_file, str):
            plotting["run_snapshot"] = {
                "profile_file": snapshot_profile_file,
                "content_hash": snapshot_profile_hash,
            }
    elif plotting_snapshot is not None:
        plotting["source"] = "run-snapshot"
        plotting["profile_file"] = (
            snapshot_profile_file or str(plotting_snapshot)
        )
        plotting["content_hash"] = _fnv1a64(plotting_snapshot)
        if (
            isinstance(snapshot_profile_hash, str)
            and snapshot_profile_hash != plotting["content_hash"]
        ):
            plotting["run_manifest_content_hash"] = snapshot_profile_hash
        plotting["snapshot_file"] = str(plotting_snapshot)
    elif isinstance(snapshot_profile_file, str):
        plotting["source"] = "run-manifest"
        plotting["profile_file"] = snapshot_profile_file
        plotting["content_hash"] = snapshot_profile_hash
    else:
        plotting["source"] = "built-in-profile"
        plotting["profile"] = profile_name or "chapter6"

    overrides: dict[str, object] = {}
    if explicit_profile:
        overrides["profile"] = (
            str(profile_file) if profile_file is not None else profile_name or "chapter6"
        )
    if requested_formats is not None:
        overrides["formats"] = list(requested_formats)

    return {
        "parameters": parameters,
        "plotting": plotting,
        "effective": effective_profile_record(profile, output_formats),
        "overrides": overrides,
    }


def _profile_from_run_snapshot(
    source: Path, fallback: PostprocessProfile, use_snapshot_style: bool = True
) -> PostprocessProfile:
    """Assemble the effective profile from the run snapshot and its ``.prm``."""

    snapshot = _find_run_snapshot(source)
    if snapshot is None:
        return fallback
    _, parameter_file, plotting_snapshot, manifest_path = snapshot
    if parameter_file is None:
        if manifest_path is None:
            return fallback
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        comparison = manifest.get("plotting", {}).get("resolved_comparison", {})
        return with_comparison_plan(
            fallback,
            str(comparison.get("rows", "")),
            str(comparison.get("columns", "")),
            str(comparison.get("group_by", "")),
        )
    if use_snapshot_style and plotting_snapshot is not None:
        fallback = load_json_profile(plotting_snapshot)
    configuration = load_postprocess_configuration(parameter_file)
    return with_parameter_configuration(fallback, configuration)


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
    selected_profile_file = getattr(arguments, "profile_file", None)
    if selected_profile_file is not None:
        selected_profile_file = selected_profile_file.resolve()
    try:
        if arguments.artifact is not None:
            artifact = arguments.artifact.resolve()
            profile = _profile_from_run_snapshot(
                artifact, profile, use_snapshot_style=not explicit_profile
            )
            output_formats = _effective_output_formats(
                arguments.output_formats, profile
            )
            provenance = _postprocess_provenance(
                artifact,
                profile,
                output_formats,
                explicit_profile=explicit_profile,
                profile_file=selected_profile_file,
                profile_name=profile_name,
                requested_formats=arguments.output_formats,
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
                provenance=provenance,
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
        provenance = _postprocess_provenance(
            input_root,
            profile,
            output_formats,
            explicit_profile=explicit_profile,
            profile_file=selected_profile_file,
            profile_name=profile_name,
            requested_formats=arguments.output_formats,
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
            provenance=provenance,
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
