"""Chapter 6 field profile for the reusable post-processing pipeline."""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path
from typing import Any

from .pipeline import ComparisonPlan, FieldSpec, PostprocessProfile


def is_b2_case(metadata: dict[str, str]) -> bool:
    scenario = metadata.get("identity.scenario_id", "")
    return scenario == "chapter-6.b2.graetz-flow" or scenario.startswith(
        "chapter-6.b2.graetz-flow."
    )


def axes_title(metadata: dict[str, str], field_name: str) -> str:
    if is_b2_case(metadata):
        output_id = b2_case_label(metadata)
    else:
        output_id = metadata.get("identity.output_id", "Chapter 6 artifact")
    return f"{output_id}: {field_name}"


def comparison_group(metadata: dict[str, str]) -> str:
    scenario = metadata.get("identity.scenario_id", "")
    parts = scenario.split(".")
    if len(parts) >= 3 and parts[0] == "chapter-6":
        return ".".join(parts[:3])
    return metadata.get("identity.recipe_id", "unknown")


def comparison_title(field: FieldSpec) -> str:
    titles = {
        "state": "B2 state with optimized control",
        "state-no-control": "B2 state without control (u = 0)",
        "target": "B2 target function",
        "forcing": "B2 forcing function",
        "observation-region": "B2 observation domains",
        "adjoint": "B2 adjoint",
        "negative-adjoint": "B2 negative adjoint (-p)",
    }
    return titles.get(field.output_name, f"Chapter 6 comparison: {field.output_name}")


def b2_case_label(metadata: dict[str, str]) -> str:
    region = metadata.get("b2.observation_region", "unknown domain")
    target = metadata.get("b2.target_profile", "unknown target")
    return f"{region} / {target} target"


def comparison_sort_key(metadata: dict[str, str]) -> tuple[object, ...]:
    """Order B1 panels as method rows and beta columns."""

    scenario = metadata.get("identity.scenario_id", "")
    if scenario == "chapter-6.b1.distributed-laplace":
        method_order = {"steepest-descent": 0, "l-bfgs": 1}
        method = metadata.get("benchmark.method", "")
        try:
            beta = -float(metadata.get("benchmark.regularisation", "nan"))
        except ValueError:
            beta = float("inf")
        return (method_order.get(method, 99), beta, method)
    if is_b2_case(metadata):
        target_order = {"constant": 0, "parabolic": 1}
        region_order = {"wings": 0, "full": 1}
        return (
            target_order.get(metadata.get("b2.target_profile", ""), 99),
            region_order.get(metadata.get("b2.observation_region", ""), 99),
            metadata.get("identity.output_id", ""),
        )
    return (0, 0.0, metadata.get("identity.output_id", ""))


CHAPTER6_PROFILE = PostprocessProfile(
    volume_source_names=("fields-volume.vtu", "fields.vtu"),
    boundary_source_names=("control-boundary.vtu", "control.vtu"),
    volume_fields=(
        FieldSpec("state", "state"),
        FieldSpec("state-no-control", "state_uncontrolled", "state (u = 0)"),
        FieldSpec("target", "target", "target function"),
        FieldSpec("forcing", "forcing", "forcing function"),
        FieldSpec("observation-region", "observation_region", "observation domain"),
        FieldSpec("adjoint", "adjoint"),
        FieldSpec("negative-adjoint", "negative_adjoint", "book adjoint (-p)"),
        FieldSpec("control", "control"),
    ),
    boundary_fields=(FieldSpec("control-boundary", "control"),),
    item_label=lambda metadata: (
        b2_case_label(metadata)
        if is_b2_case(metadata)
        else metadata.get("identity.output_id", "Chapter 6 artifact")
    ),
    field_title=lambda metadata, field: axes_title(metadata, field.output_name),
    comparison_group=comparison_group,
    comparison_title=comparison_title,
    comparison_sort_key=comparison_sort_key,
    missing_fields_message=(
        "no configured state, target, forcing, adjoint, or control fields found"
    ),
)


def _axis_value(metadata: dict[str, str], axis: str) -> str:
    for key in (f"parameters.{axis}", f"benchmark.{axis}", f"b2.{axis}", axis):
        if key in metadata:
            return metadata[key]
    return ""


def _format_template(template: str, metadata: dict[str, str]) -> str:
    values: dict[str, str] = {}
    for key, value in metadata.items():
        values[key] = value
        values[key.rsplit(".", 1)[-1]] = value
    for axis in ("method", "regularisation", "forcing", "observation-region", "target-profile"):
        values[axis] = _axis_value(metadata, axis)
    try:
        return template.format_map(values)
    except (KeyError, ValueError):
        return template


def load_json_profile(path: Path) -> PostprocessProfile:
    """Load fields, labels, axes, and comparison policy from ``nmopt-plot-v1``."""

    document: dict[str, Any] = json.loads(path.read_text(encoding="utf-8"))
    if document.get("schema") != "nmopt-plot-v1":
        raise ValueError(f"unsupported plotting profile schema in '{path}'")
    fields = document.get("fields", {})
    volume_fields: list[FieldSpec] = []
    boundary_fields: list[FieldSpec] = []
    field_titles: dict[str, str] = {}
    for output_name, specification in fields.items():
        if not isinstance(specification, dict) or "source" not in specification:
            raise ValueError(f"plotting field '{output_name}' needs a source")
        field = FieldSpec(
            output_name,
            str(specification["source"]),
            specification.get("colorbar_label"),
        )
        field_titles[output_name] = str(specification.get("title", output_name))
        if specification.get("mesh") == "boundary":
            boundary_fields.append(field)
        else:
            volume_fields.append(field)

    axes = document.get("axes", {})
    axis_orders = {
        axis: tuple(str(value) for value in specification.get("order", []))
        for axis, specification in axes.items()
        if isinstance(specification, dict)
    }
    axis_labels = {
        axis: {
            str(value): str(label)
            for value, label in specification.get("labels", {}).items()
        }
        for axis, specification in axes.items()
        if isinstance(specification, dict)
    }
    comparison = document.get("default_comparison", {})
    plan = ComparisonPlan(
        tuple(str(axis) for axis in comparison.get("rows", [])),
        tuple(str(axis) for axis in comparison.get("columns", [])),
        tuple(str(axis) for axis in comparison.get("group_by", [])),
    )
    templates = {
        str(field): str(template)
        for field, template in document.get("title_templates", {}).items()
    }

    def item_label(metadata: dict[str, str]) -> str:
        scenario = metadata.get("identity.scenario_id", "Chapter 6 artifact")
        values = [
            axis_labels.get(axis, {}).get(_axis_value(metadata, axis), _axis_value(metadata, axis))
            for axis in plan.rows + plan.columns
            if _axis_value(metadata, axis)
        ]
        return " / ".join(values) if values else scenario

    def field_title(metadata: dict[str, str], field: FieldSpec) -> str:
        template = templates.get(field.output_name, field_titles.get(field.output_name, field.output_name))
        return _format_template(template, metadata)

    def comparison_title(field: FieldSpec) -> str:
        return field_titles.get(field.output_name, f"Chapter 6 comparison: {field.output_name}")

    def comparison_sort_key(metadata: dict[str, str]) -> tuple[object, ...]:
        values: list[object] = []
        for axis in plan.rows + plan.columns:
            value = _axis_value(metadata, axis)
            order = axis_orders.get(axis, ())
            values.append(order.index(value) if value in order else len(order))
        return tuple(values) + (metadata.get("identity.output_id", ""),)

    return PostprocessProfile(
        volume_source_names=CHAPTER6_PROFILE.volume_source_names,
        boundary_source_names=CHAPTER6_PROFILE.boundary_source_names,
        volume_fields=tuple(volume_fields),
        boundary_fields=tuple(boundary_fields),
        item_label=item_label,
        field_title=field_title,
        comparison_group=comparison_group,
        comparison_title=comparison_title,
        comparison_sort_key=comparison_sort_key,
        missing_fields_message=CHAPTER6_PROFILE.missing_fields_message,
        comparison_plan=plan,
        axis_orders=axis_orders,
        axis_labels=axis_labels,
    )


def with_comparison_plan(
    profile: PostprocessProfile, rows: str, columns: str, group_by: str
) -> PostprocessProfile:
    """Apply the run manifest's comparison override to a loaded style profile."""

    split = lambda value: tuple(item for item in value.split(",") if item and item != "none")
    return replace(
        profile,
        comparison_plan=ComparisonPlan(split(rows), split(columns), split(group_by)),
    )
