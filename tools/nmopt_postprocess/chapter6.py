"""Chapter 6 field profile for the reusable post-processing pipeline."""

from .pipeline import FieldSpec, PostprocessProfile


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
