"""Profile-driven artifact and run-root post-processing."""

from __future__ import annotations

import json
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import meshio

from .errors import PostprocessError
from .fields import ScalarField, find_field
from .records import read_metadata, read_numeric_history
from .render import (
    DEFAULT_OUTPUT_FORMATS,
    OutputFormats,
    HistoryLine,
    HistoryPanel,
    RenderItem,
    plot_boundary_comparison,
    plot_boundary_field,
    plot_history_figure,
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
class ComparisonPlan:
    """Named matrix axes that determine comparison rows and columns."""

    rows: tuple[str, ...] = ()
    columns: tuple[str, ...] = ()
    group_by: tuple[str, ...] = ()


@dataclass(frozen=True)
class HistorySeriesSpec:
    """Select and style one artifact in an optimization-history figure."""

    selector: tuple[tuple[str, str], ...]
    label: str
    color: str | None = None
    linestyle: str = "-"
    marker: str | None = None


@dataclass(frozen=True)
class HistoryPlotSpec:
    """Select one persisted history for one figure panel."""

    source: str
    title: str
    x_label: str
    y_label: str
    x_scale: str = "linear"
    y_scale: str = "linear"
    iteration_origin: int = 0
    x_limits: tuple[float, float] | None = None
    y_limits: tuple[float, float] | None = None


@dataclass(frozen=True)
class HistoryFigureSpec:
    """A profile-defined multi-panel history comparison."""

    output_name: str
    plots: tuple[HistoryPlotSpec, ...]
    series: tuple[HistorySeriesSpec, ...]


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
    comparison_sort_key: Callable[[dict[str, str]], tuple[object, ...]] | None = None
    missing_fields_message: str = "no configured fields found"
    comparison_plan: ComparisonPlan = ComparisonPlan()
    axis_orders: dict[str, tuple[str, ...]] | None = None
    axis_labels: dict[str, dict[str, str]] | None = None
    history_figures: tuple[HistoryFigureSpec, ...] = ()


def _render_item(
    artifact: Path,
    field: ScalarField,
    profile: PostprocessProfile,
    mesh: meshio.Mesh,
) -> RenderItem:
    metadata = read_metadata(artifact)
    return RenderItem(profile.item_label(metadata), mesh, field, metadata)


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


def _matches_history_series(
    metadata: dict[str, str], series: HistorySeriesSpec
) -> bool:
    return all(
        _axis_value(metadata, axis) == value for axis, value in series.selector
    )


def _history_panels(
    artifacts: list[Path], figure: HistoryFigureSpec
) -> tuple[HistoryPanel, ...]:
    metadata_by_artifact = [
        (artifact, read_metadata(artifact)) for artifact in artifacts
    ]
    panels: list[HistoryPanel] = []
    for plot in figure.plots:
        lines: list[HistoryLine] = []
        for series in figure.series:
            matches = [
                metadata
                for _, metadata in metadata_by_artifact
                if _matches_history_series(metadata, series)
            ]
            if len(matches) > 1:
                selector = ", ".join(
                    f"{axis}={value}" for axis, value in series.selector
                )
                raise PostprocessError(
                    f"history selector '{selector}' matches more than one artifact"
                )
            if not matches:
                continue
            values = read_numeric_history(matches[0], plot.source)
            if not values:
                continue
            lines.append(
                HistoryLine(
                    series.label,
                    values,
                    series.color,
                    series.linestyle,
                    series.marker,
                )
            )
        panels.append(
            HistoryPanel(
                plot.title,
                plot.x_label,
                plot.y_label,
                plot.x_scale,
                plot.y_scale,
                plot.iteration_origin,
                tuple(lines),
                plot.x_limits,
                plot.y_limits,
            )
        )
    return tuple(panels)


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
        if profile.comparison_sort_key is not None:
            group_artifacts.sort(
                key=lambda artifact: profile.comparison_sort_key(
                    read_metadata(artifact)
                )
            )
        comparison_dir = output_root / "comparisons" / group
        comparison_dir.mkdir(parents=True, exist_ok=True)
        group_generated: dict[str, list[str]] = {}
        group_errors: dict[str, str] = {}
        rows, columns = comparison_dimensions(group_artifacts, profile)

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
                        rows=rows,
                        columns=columns,
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
                        rows=rows,
                        columns=columns,
                        output_formats=output_formats,
                    )
            except (OSError, ValueError, PostprocessError) as error:
                group_errors[field_spec.output_name] = str(error)

        for figure in profile.history_figures:
            try:
                group_generated[figure.output_name] = plot_history_figure(
                    _history_panels(group_artifacts, figure),
                    comparison_dir / figure.output_name,
                    output_formats=output_formats,
                )
            except (OSError, ValueError, PostprocessError) as error:
                group_errors[figure.output_name] = str(error)

        if group_generated:
            generated[group] = group_generated
        if group_errors:
            errors[group] = group_errors

    return generated, errors


def _axis_value(metadata: dict[str, str], axis: str) -> str:
    for key in (f"parameters.{axis}", f"benchmark.{axis}", f"b2.{axis}", axis):
        if key in metadata:
            return metadata[key]
    return ""


def _axis_values(
    artifacts: list[Path], profile: PostprocessProfile, axis: str
) -> list[str]:
    declared = list((profile.axis_orders or {}).get(axis, ()))
    observed = {
        value
        for artifact in artifacts
        if (value := _axis_value(read_metadata(artifact), axis))
    }
    return declared + sorted(observed.difference(declared))


def comparison_dimensions(
    artifacts: list[Path], profile: PostprocessProfile
) -> tuple[int, int]:
    """Resolve an explicit row/column plan without using panel-count magic."""

    plan = profile.comparison_plan
    row_count = 1
    column_count = 1
    if plan.rows:
        for axis in plan.rows:
            row_count *= max(1, len(_axis_values(artifacts, profile, axis)))
    if plan.columns:
        for axis in plan.columns:
            column_count *= max(1, len(_axis_values(artifacts, profile, axis)))
    if not plan.columns:
        column_count = max(1, len(artifacts))
    return row_count, column_count


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
        artifacts, output_root, profile, output_formats=output_formats
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
