"""Chapter 6 field profile for the reusable post-processing pipeline."""

from .pipeline import FieldSpec, PostprocessProfile


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
