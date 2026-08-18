#!/usr/bin/env python3
"""Render Chapter 6 native fields with meshio and matplotlib.

The tool reads one runner artifact directory and writes derived plots outside
the authoritative ``native/`` directory.  It accepts the current native file
names and the historical ``fields.vtu``/``control.vtu`` names used by older
development runs.  With ``--input`` it processes every artifact below a run
root into a mirrored output tree and writes an aggregate index.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from nmopt_postprocess import PostprocessError
from nmopt_postprocess.pipeline import (
    FieldSpec,
    PostprocessProfile,
    process_artifact as process_profile_artifact,
    process_input_root as process_profile_input_root,
)


def axes_title(metadata: dict[str, str], field_name: str) -> str:
    output_id = metadata.get("identity.output_id", "Chapter 6 artifact")
    return f"{output_id}: {field_name}"


def comparison_group(metadata: dict[str, str]) -> str:
    scenario = metadata.get("identity.scenario_id", "")
    parts = scenario.split(".")
    if len(parts) >= 3 and parts[0] == "chapter-6":
        return ".".join(parts[:3])
    return metadata.get("identity.recipe_id", "unknown")


def comparison_title(field: FieldSpec) -> str:
    return f"Chapter 6 comparison: {field.output_name}"


CHAPTER6_PROFILE = PostprocessProfile(
    volume_source_names=("fields-volume.vtu", "fields.vtu"),
    boundary_source_names=("control-boundary.vtu", "control.vtu"),
    volume_fields=(
        FieldSpec("state", "state"),
        FieldSpec("adjoint", "adjoint"),
        FieldSpec("control", "control"),
    ),
    boundary_fields=(FieldSpec("control-boundary", "control"),),
    item_label=lambda metadata: metadata.get(
        "identity.output_id", "Chapter 6 artifact"
    ),
    field_title=lambda metadata, field: axes_title(metadata, field.output_name),
    comparison_group=comparison_group,
    comparison_title=comparison_title,
    missing_fields_message="no state, adjoint, or control fields found",
)


def process_artifact(artifact: Path, output: Path) -> dict[str, object]:
    return process_profile_artifact(artifact, output, CHAPTER6_PROFILE)


def process_input_root(input_root: Path, output_root: Path) -> dict[str, object]:
    return process_profile_input_root(input_root, output_root, CHAPTER6_PROFILE)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
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
            "derived output directory; defaults to <artifact>/derived for "
            "--artifact or <input>/postprocessed for --input"
        ),
    )
    arguments = parser.parse_args()
    try:
        if arguments.artifact is not None:
            artifact = arguments.artifact.resolve()
            output = (
                arguments.output.resolve()
                if arguments.output is not None
                else artifact / "derived"
            )
            manifest = process_artifact(artifact, output)
            print(f"wrote {len(manifest['plots'])} field plot sets to {output}")
            return 0

        input_root = arguments.input.resolve()
        output_root = (
            arguments.output.resolve()
            if arguments.output is not None
            else input_root / "postprocessed"
        )
        index = process_input_root(input_root, output_root)
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
