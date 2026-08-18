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
import json
from pathlib import Path

import meshio

from nmopt_postprocess import (
    PostprocessError,
    find_field,
    first_existing,
    read_metadata,
)
from nmopt_postprocess.render import (
    RenderItem,
    plot_boundary_comparison,
    plot_boundary_field,
    plot_volume_comparison,
    plot_volume_field,
)


def axes_title(metadata: dict[str, str], field_name: str) -> str:
    output_id = metadata.get("identity.output_id", "Chapter 6 artifact")
    return f"{output_id}: {field_name}"


def volume_comparison_items(
    artifacts: list[Path], field_name: str
) -> list[RenderItem]:
    items: list[RenderItem] = []
    for artifact in artifacts:
        fields_path = first_existing(artifact, ("fields-volume.vtu", "fields.vtu"))
        if fields_path is None:
            continue
        mesh = meshio.read(fields_path)
        field = find_field(mesh, field_name)
        if field is not None:
            metadata = read_metadata(artifact)
            items.append(
                RenderItem(
                    metadata.get("identity.output_id", "Chapter 6 artifact"),
                    mesh,
                    field,
                )
            )
    return items


def boundary_comparison_items(
    artifacts: list[Path],
) -> list[RenderItem]:
    items: list[RenderItem] = []
    for artifact in artifacts:
        boundary_path = first_existing(
            artifact, ("control-boundary.vtu", "control.vtu")
        )
        if boundary_path is None:
            continue
        mesh = meshio.read(boundary_path)
        field = find_field(mesh, "control")
        if field is not None:
            metadata = read_metadata(artifact)
            items.append(
                RenderItem(
                    metadata.get("identity.output_id", "Chapter 6 artifact"),
                    mesh,
                    field,
                )
            )
    return items


def comparison_group(artifact: Path) -> str:
    metadata = read_metadata(artifact)
    scenario = metadata.get("identity.scenario_id", "")
    parts = scenario.split(".")
    if len(parts) >= 3 and parts[0] == "chapter-6":
        return ".".join(parts[:3])
    return metadata.get("identity.recipe_id", "unknown")


def build_comparisons(
    artifacts: list[Path], output_root: Path
) -> tuple[dict[str, dict[str, list[str]]], dict[str, dict[str, str]]]:
    grouped: dict[str, list[Path]] = {}
    for artifact in artifacts:
        grouped.setdefault(comparison_group(artifact), []).append(artifact)

    generated: dict[str, dict[str, list[str]]] = {}
    errors: dict[str, dict[str, str]] = {}
    for group, group_artifacts in grouped.items():
        comparison_dir = output_root / "comparisons" / group
        comparison_dir.mkdir(parents=True, exist_ok=True)
        group_generated: dict[str, list[str]] = {}
        group_errors: dict[str, str] = {}
        for field_name in ("state", "adjoint", "control"):
            try:
                items = volume_comparison_items(group_artifacts, field_name)
                if len(items) >= 2:
                    group_generated[field_name] = plot_volume_comparison(
                        items,
                        field_name,
                        comparison_dir / f"comparison-{field_name}",
                        title=f"Chapter 6 comparison: {field_name}",
                    )
            except (OSError, ValueError, PostprocessError) as error:
                group_errors[field_name] = str(error)

        try:
            items = boundary_comparison_items(group_artifacts)
            if len(items) >= 2:
                group_generated["control-boundary"] = plot_boundary_comparison(
                    items,
                    comparison_dir / "comparison-control-boundary",
                    title="Chapter 6 comparison: control-boundary",
                )
        except (OSError, ValueError, PostprocessError) as error:
            group_errors["control-boundary"] = str(error)

        if group_generated:
            generated[group] = group_generated
        if group_errors:
            errors[group] = group_errors

    return generated, errors


def process_artifact(artifact: Path, output: Path) -> dict[str, object]:
    metadata = read_metadata(artifact)
    fields_path = first_existing(artifact, ("fields-volume.vtu", "fields.vtu"))
    if fields_path is None:
        raise PostprocessError(
            f"no volume field file found below '{artifact / 'native'}'"
        )

    output.mkdir(parents=True, exist_ok=True)
    volume_mesh = meshio.read(fields_path)
    generated: dict[str, list[str]] = {}
    for field_name in ("state", "adjoint", "control"):
        field = find_field(volume_mesh, field_name)
        if field is not None:
            generated[field_name] = plot_volume_field(
                volume_mesh,
                field,
                axes_title(metadata, field_name),
                output / field_name,
            )

    boundary_path = first_existing(artifact, ("control-boundary.vtu", "control.vtu"))
    if boundary_path is not None:
        boundary_mesh = meshio.read(boundary_path)
        boundary_field = find_field(boundary_mesh, "control")
        if boundary_field is not None:
            generated["control-boundary"] = plot_boundary_field(
                boundary_mesh,
                boundary_field,
                axes_title(metadata, "control-boundary"),
                output / "control-boundary",
            )

    if not generated:
        raise PostprocessError("no state, adjoint, or control fields found")

    manifest = {
        "artifact_directory": str(artifact.resolve()),
        "volume_source": str(fields_path.relative_to(artifact)),
        "boundary_source": (
            str(boundary_path.relative_to(artifact))
            if boundary_path is not None
            else None
        ),
        "plots": generated,
    }
    (output / "postprocess.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest


def process_input_root(input_root: Path, output_root: Path) -> dict[str, object]:
    artifacts = sorted(path.parent for path in input_root.rglob("artifact.kv"))
    if not artifacts:
        raise PostprocessError(f"no artifact.kv files found below '{input_root}'")

    output_root.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, object]] = []
    successful_artifacts: list[Path] = []
    for artifact in artifacts:
        relative = artifact.relative_to(input_root)
        output = output_root / relative
        try:
            manifest = process_artifact(artifact, output)
            records.append(
                {
                    "artifact": str(relative),
                    "status": "ok",
                    "plots": manifest["plots"],
                }
            )
            successful_artifacts.append(artifact)
        except (OSError, ValueError, PostprocessError) as error:
            output.mkdir(parents=True, exist_ok=True)
            (output / "postprocess.json").write_text(
                json.dumps(
                    {
                        "artifact_directory": str(artifact.resolve()),
                        "status": "error",
                        "error": str(error),
                    },
                    indent=2,
                    sort_keys=True,
                )
                + "\n",
                encoding="utf-8",
            )
            records.append(
                {
                    "artifact": str(relative),
                    "status": "error",
                    "error": str(error),
                }
            )

    comparisons, comparison_errors = build_comparisons(
        successful_artifacts,
        output_root,
    )
    index = {
        "input_root": str(input_root.resolve()),
        "output_root": str(output_root.resolve()),
        "artifact_count": len(records),
        "success_count": sum(record["status"] == "ok" for record in records),
        "failure_count": sum(record["status"] == "error" for record in records),
        "artifacts": records,
        "comparisons": comparisons,
    }
    if comparison_errors:
        index["comparison_errors"] = comparison_errors
    (output_root / "postprocess-index.json").write_text(
        json.dumps(index, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return index


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
