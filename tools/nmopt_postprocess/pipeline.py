"""Profile-driven artifact and run-root post-processing."""

from __future__ import annotations

import json
from collections.abc import Callable
from dataclasses import dataclass
from itertools import product
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
    matrix_axis_values: dict[str, tuple[str, ...]] | None = None
    matrix_combinations: tuple[dict[str, str], ...] = ()
    output_formats: OutputFormats = DEFAULT_OUTPUT_FORMATS


@dataclass(frozen=True)
class ComparisonGrid:
    """Validated row-major artifact order and dimensions for one comparison."""

    artifacts: tuple[Path, ...]
    rows: int
    columns: int


def effective_profile_record(
    profile: PostprocessProfile, output_formats: OutputFormats
) -> dict[str, object]:
    """Return the JSON-safe plotting choices applied to derived outputs."""

    plan = profile.comparison_plan
    return {
        "fields": {
            "volume": [field.output_name for field in profile.volume_fields],
            "boundary": [field.output_name for field in profile.boundary_fields],
        },
        "matrix_axes": {
            axis: list(values)
            for axis, values in (profile.matrix_axis_values or {}).items()
        },
        "matrix_combinations": [
            dict(combination) for combination in profile.matrix_combinations
        ],
        "comparison": {
            "rows": list(plan.rows),
            "columns": list(plan.columns),
            "group_by": list(plan.group_by),
        },
        "formats": list(output_formats),
    }


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
        metadata = read_metadata(artifact)
        base_group = profile.comparison_group(metadata)
        group_by = profile.comparison_plan.group_by
        if group_by:
            suffix = "__".join(
                f"{axis}={_axis_value(metadata, axis)}" for axis in group_by
            )
            base_group = f"{base_group}/{suffix}"
        grouped.setdefault(base_group, []).append(artifact)

    generated: dict[str, dict[str, list[str]]] = {}
    errors: dict[str, dict[str, str]] = {}
    for group, group_artifacts in grouped.items():
        comparison_dir = output_root / "comparisons" / group
        comparison_dir.mkdir(parents=True, exist_ok=True)
        group_generated: dict[str, list[str]] = {}
        group_errors: dict[str, str] = {}
        try:
            grid = comparison_grid(group_artifacts, profile)
        except (ValueError, PostprocessError) as error:
            errors[group] = {"comparison": str(error)}
            continue
        group_artifacts = list(grid.artifacts)
        rows, columns = grid.rows, grid.columns

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
    """Return observed axis values in profile-declared presentation order."""

    matrix_values = (profile.matrix_axis_values or {}).get(axis)
    if matrix_values is not None:
        return list(matrix_values)

    observed = {
        value
        for artifact in artifacts
        if (value := _axis_value(read_metadata(artifact), axis))
    }
    declared = (profile.matrix_axis_values or {}).get(
        axis, (profile.axis_orders or {}).get(axis, ())
    )
    ordered = []
    seen: set[str] = set()
    for value in declared:
        if value in observed and value not in seen:
            ordered.append(value)
            seen.add(value)
    ordered.extend(sorted(observed.difference(seen)))
    return ordered


def comparison_dimensions(
    artifacts: list[Path], profile: PostprocessProfile
) -> tuple[int, int]:
    """Return dimensions after validating one comparison matrix."""

    grid = comparison_grid(artifacts, profile)
    return grid.rows, grid.columns


def _axis_tuples(
    artifacts: list[Path], profile: PostprocessProfile, axes: tuple[str, ...]
) -> list[tuple[str, ...]]:
    if not axes:
        return [()]
    values = [_axis_values(artifacts, profile, axis) for axis in axes]
    if any(not axis_values for axis_values in values):
        missing = [axis for axis, axis_values in zip(axes, values) if not axis_values]
        raise PostprocessError(
            "comparison plan references axes without active values: "
            + ", ".join(missing)
        )
    return list(product(*values))


def _validate_plan_axes(profile: PostprocessProfile) -> None:
    plan = profile.comparison_plan
    planned = plan.rows + plan.columns + plan.group_by
    if len(set(planned)) != len(planned):
        raise PostprocessError("comparison axes may only appear in one plan section")
    if profile.matrix_axis_values is None:
        return
    declared = set(profile.matrix_axis_values)
    unknown = set(planned).difference(declared)
    if unknown:
        raise PostprocessError(
            "comparison plan references undeclared matrix axes: "
            + ", ".join(sorted(unknown))
        )
    unbound = {
        axis
        for axis, values in profile.matrix_axis_values.items()
        if len(values) > 1 and axis not in planned
    }
    if unbound:
        raise PostprocessError(
            "varying matrix axes must be rows, columns, or group_by: "
            + ", ".join(sorted(unbound))
        )


def comparison_grid(
    artifacts: list[Path], profile: PostprocessProfile
) -> ComparisonGrid:
    """Validate and order one explicit comparison matrix in row-major order."""

    if not artifacts:
        raise PostprocessError("comparison matrix has no artifacts")
    plan = profile.comparison_plan
    if not plan.rows and not plan.columns and not profile.matrix_combinations:
        ordered = list(artifacts)
        if profile.comparison_sort_key is not None:
            ordered.sort(
                key=lambda artifact: profile.comparison_sort_key(
                    read_metadata(artifact)
                )
            )
        return ComparisonGrid(tuple(ordered), 1, len(ordered))
    _validate_plan_axes(profile)
    group_values = {
        axis: _axis_value(read_metadata(artifacts[0]), axis)
        for axis in plan.group_by
    }
    for artifact in artifacts[1:]:
        metadata = read_metadata(artifact)
        if any(
            _axis_value(metadata, axis) != value
            for axis, value in group_values.items()
        ):
            raise PostprocessError("comparison group contains mixed group_by values")

    row_keys = _axis_tuples(artifacts, profile, plan.rows)
    column_keys = _axis_tuples(artifacts, profile, plan.columns)
    expected_keys = {
        (row_key, column_key) for row_key in row_keys for column_key in column_keys
    }
    if profile.matrix_combinations:
        combinations = [
            combination
            for combination in profile.matrix_combinations
            if all(combination.get(axis) == value for axis, value in group_values.items())
        ]
        if not combinations:
            raise PostprocessError("comparison group has no declared matrix combination")
        expected_keys = {
            (
                tuple(combination.get(axis, "") for axis in plan.rows),
                tuple(combination.get(axis, "") for axis in plan.columns),
            )
            for combination in combinations
        }

    coordinates: dict[tuple[tuple[str, ...], tuple[str, ...]], Path] = {}
    for artifact in artifacts:
        metadata = read_metadata(artifact)
        key = (
            tuple(_axis_value(metadata, axis) for axis in plan.rows),
            tuple(_axis_value(metadata, axis) for axis in plan.columns),
        )
        if key in coordinates:
            raise PostprocessError(
                "comparison matrix has duplicate coordinates: "
                + repr(key)
            )
        coordinates[key] = artifact

    missing = expected_keys.difference(coordinates)
    extra = set(coordinates).difference(expected_keys)
    if missing or extra:
        details = []
        if missing:
            details.append("missing " + ", ".join(map(repr, sorted(missing))))
        if extra:
            details.append("unexpected " + ", ".join(map(repr, sorted(extra))))
        raise PostprocessError(
            "comparison matrix is incomplete: " + "; ".join(details)
        )

    ordered = tuple(
        coordinates[(row_key, column_key)]
        for row_key in row_keys
        for column_key in column_keys
        if (row_key, column_key) in expected_keys
    )
    return ComparisonGrid(ordered, len(row_keys), len(column_keys))


def process_artifact(
    artifact: Path,
    output: Path,
    profile: PostprocessProfile,
    output_formats: OutputFormats = DEFAULT_OUTPUT_FORMATS,
    provenance: dict[str, object] | None = None,
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
    if provenance is not None:
        manifest["provenance"] = provenance
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
    provenance: dict[str, object] | None = None,
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
                artifact,
                output,
                profile,
                output_formats=output_formats,
                provenance=provenance,
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
            manifest = {
                "artifact_directory": str(artifact.resolve()),
                "formats": list(output_formats),
                "status": "error",
                "error": str(error),
            }
            if provenance is not None:
                manifest["provenance"] = provenance
            (output / "postprocess.json").write_text(
                json.dumps(
                    manifest,
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
    if provenance is not None:
        index["provenance"] = provenance
    if comparison_errors:
        index["comparison_errors"] = comparison_errors
    (output_root / "postprocess-index.json").write_text(
        json.dumps(index, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return index
