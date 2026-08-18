"""Profile-driven artifact and run-root post-processing."""

from __future__ import annotations

import json
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import meshio

from .errors import PostprocessError
from .fields import ScalarField, find_field
from .records import read_metadata
from .render import (
    DEFAULT_OUTPUT_FORMATS,
    OutputFormats,
    RenderItem,
    plot_boundary_comparison,
    plot_boundary_field,
    plot_volume_comparison,
    plot_volume_field,
)
from .sources import first_existing


@dataclass(frozen=True)
class FieldSpec:
    """Map one persisted field name to its derived output name."""

    output_name: str
    mesh_name: str
    colorbar_label: str | None = None


@dataclass(frozen=True)
class PostprocessProfile:
    """Application-specific field and labeling choices for the generic engine."""

    volume_source_names: tuple[str, ...]
    boundary_source_names: tuple[str, ...]
    volume_fields: tuple[FieldSpec, ...]
    boundary_fields: tuple[FieldSpec, ...]
    item_label: Callable[[dict[str, str]], str]
    field_title: Callable[[dict[str, str], FieldSpec], str]
    comparison_group: Callable[[dict[str, str]], str]
    comparison_title: Callable[[FieldSpec], str]
    missing_fields_message: str = "no configured fields found"


def _render_item(
    artifact: Path,
    field: ScalarField,
    profile: PostprocessProfile,
    mesh: meshio.Mesh,
) -> RenderItem:
    metadata = read_metadata(artifact)
    return RenderItem(profile.item_label(metadata), mesh, field)


def _volume_items(
    artifacts: list[Path], field_spec: FieldSpec, profile: PostprocessProfile
) -> list[RenderItem]:
    items: list[RenderItem] = []
    for artifact in artifacts:
        fields_path = first_existing(artifact, profile.volume_source_names)
        if fields_path is None:
            continue
        mesh = meshio.read(fields_path)
        field = find_field(mesh, field_spec.mesh_name)
        if field is not None:
            items.append(_render_item(artifact, field, profile, mesh))
    return items


def _boundary_items(
    artifacts: list[Path], field_spec: FieldSpec, profile: PostprocessProfile
) -> list[RenderItem]:
    items: list[RenderItem] = []
    for artifact in artifacts:
        boundary_path = first_existing(artifact, profile.boundary_source_names)
        if boundary_path is None:
            continue
        mesh = meshio.read(boundary_path)
        field = find_field(mesh, field_spec.mesh_name)
        if field is not None:
            items.append(_render_item(artifact, field, profile, mesh))
    return items


def build_comparisons(
    artifacts: list[Path],
    output_root: Path,
    profile: PostprocessProfile,
    output_formats: OutputFormats = DEFAULT_OUTPUT_FORMATS,
) -> tuple[dict[str, dict[str, list[str]]], dict[str, dict[str, str]]]:
    """Render all configured comparisons grouped by the application profile."""

    grouped: dict[str, list[Path]] = {}
    for artifact in artifacts:
        grouped.setdefault(
            profile.comparison_group(read_metadata(artifact)), []
        ).append(artifact)

    generated: dict[str, dict[str, list[str]]] = {}
    errors: dict[str, dict[str, str]] = {}
    for group, group_artifacts in grouped.items():
        comparison_dir = output_root / "comparisons" / group
        comparison_dir.mkdir(parents=True, exist_ok=True)
        group_generated: dict[str, list[str]] = {}
        group_errors: dict[str, str] = {}

        for field_spec in profile.volume_fields:
            try:
                items = _volume_items(group_artifacts, field_spec, profile)
                if len(items) >= 2:
                    group_generated[field_spec.output_name] = plot_volume_comparison(
                        items,
                        field_spec.mesh_name,
                        comparison_dir / f"comparison-{field_spec.output_name}",
                        title=profile.comparison_title(field_spec),
                        colorbar_label=field_spec.colorbar_label,
                        output_formats=output_formats,
                    )
            except (OSError, ValueError, PostprocessError) as error:
                group_errors[field_spec.output_name] = str(error)

        for field_spec in profile.boundary_fields:
            try:
                items = _boundary_items(group_artifacts, field_spec, profile)
                if len(items) >= 2:
                    group_generated[field_spec.output_name] = plot_boundary_comparison(
                        items,
                        comparison_dir / f"comparison-{field_spec.output_name}",
                        title=profile.comparison_title(field_spec),
                        colorbar_label=field_spec.colorbar_label,
                        output_formats=output_formats,
                    )
            except (OSError, ValueError, PostprocessError) as error:
                group_errors[field_spec.output_name] = str(error)

        if group_generated:
            generated[group] = group_generated
        if group_errors:
            errors[group] = group_errors

    return generated, errors


def process_artifact(
    artifact: Path,
    output: Path,
    profile: PostprocessProfile,
    output_formats: OutputFormats = DEFAULT_OUTPUT_FORMATS,
) -> dict[str, object]:
    """Render one artifact according to a profile and write its manifest."""

    metadata = read_metadata(artifact)
    fields_path = first_existing(artifact, profile.volume_source_names)
    if fields_path is None:
        raise PostprocessError(
            f"no volume field file found below '{artifact / 'native'}'"
        )

    output.mkdir(parents=True, exist_ok=True)
    volume_mesh = meshio.read(fields_path)
    generated: dict[str, list[str]] = {}
    for field_spec in profile.volume_fields:
        field = find_field(volume_mesh, field_spec.mesh_name)
        if field is not None:
            generated[field_spec.output_name] = plot_volume_field(
                volume_mesh,
                field,
                profile.field_title(metadata, field_spec),
                output / field_spec.output_name,
                colorbar_label=field_spec.colorbar_label,
                output_formats=output_formats,
            )

    boundary_path = first_existing(artifact, profile.boundary_source_names)
    if boundary_path is not None:
        boundary_mesh = meshio.read(boundary_path)
        for field_spec in profile.boundary_fields:
            field = find_field(boundary_mesh, field_spec.mesh_name)
            if field is not None:
                generated[field_spec.output_name] = plot_boundary_field(
                    boundary_mesh,
                    field,
                    profile.field_title(metadata, field_spec),
                    output / field_spec.output_name,
                    colorbar_label=field_spec.colorbar_label,
                    output_formats=output_formats,
                )

    if not generated:
        raise PostprocessError(profile.missing_fields_message)

    manifest = {
        "artifact_directory": str(artifact.resolve()),
        "formats": list(output_formats),
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


def process_input_root(
    input_root: Path,
    output_root: Path,
    profile: PostprocessProfile,
    output_formats: OutputFormats = DEFAULT_OUTPUT_FORMATS,
) -> dict[str, object]:
    """Process every artifact below a run root and write the aggregate index."""

    artifacts = sorted(path.parent for path in input_root.rglob("artifact.kv"))
    if not artifacts:
        raise PostprocessError(f"no artifact.kv files found below '{input_root}'")

    output_root.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, object]] = []
    successful_artifacts: list[Path] = []
    for artifact in artifacts:
        relative = artifact.relative_to(input_root)
        output = artifact / "postprocess"
        try:
            manifest = process_artifact(
                artifact, output, profile, output_formats=output_formats
            )
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
                        "formats": list(output_formats),
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
        profile,
        output_formats=output_formats,
    )
    index = {
        "input_root": str(input_root.resolve()),
        "output_root": str(output_root.resolve()),
        "formats": list(output_formats),
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
