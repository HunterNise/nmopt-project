"""Read the application-owned post-processing settings from a ``.prm`` file."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from pathlib import Path

from .pipeline import ComparisonPlan
from .render import DEFAULT_OUTPUT_FORMATS, OutputFormats


@dataclass(frozen=True)
class PostprocessConfiguration:
    """Effective post-processing selections declared by one experiment file."""

    matrix_axes: dict[str, tuple[str, ...]]
    matrix_combinations: tuple[dict[str, str], ...]
    fields: tuple[str, ...]
    style_profile: str
    comparison_plan: ComparisonPlan
    output_formats: OutputFormats = DEFAULT_OUTPUT_FORMATS


def _read_values(path: Path) -> dict[tuple[str, ...], str]:
    values: dict[tuple[str, ...], str] = {}
    sections: list[str] = []
    for line_number, line in enumerate(
        path.read_text(encoding="utf-8").splitlines(), start=1
    ):
        text = line.strip()
        if not text or text.startswith("#"):
            continue
        if text.startswith("subsection "):
            sections.append(text.removeprefix("subsection ").strip())
            continue
        if text == "end":
            if not sections:
                raise ValueError(f"unexpected 'end' in '{path}' at line {line_number}")
            sections.pop()
            continue
        if not text.startswith("set ") or "=" not in text:
            raise ValueError(f"unsupported parameter syntax in '{path}' at line {line_number}")
        key, value = text.removeprefix("set ").split("=", 1)
        values[tuple(sections) + (key.strip(),)] = value.strip()
    if sections:
        raise ValueError(f"unterminated subsection in '{path}'")
    return values


def _split(value: str | None) -> tuple[str, ...]:
    if value is None:
        return ()
    return tuple(item.strip() for item in value.split(",") if item.strip())


def _value(values: dict[tuple[str, ...], str], *path: str) -> str | None:
    return values.get(tuple(path))


def _comparison_axes(value: str | None) -> tuple[str, ...]:
    return tuple(axis for axis in _split(value) if axis != "none")


def _excluded_combinations(
    value: str | None, axes: dict[str, tuple[str, ...]], path: Path
) -> tuple[dict[str, str], ...]:
    if value is None or not value.strip():
        return ()

    excluded: list[dict[str, str]] = []
    signatures: set[tuple[tuple[str, str], ...]] = set()
    for entry in value.split(";"):
        text = entry.strip()
        if not text.startswith("[") or not text.endswith("]"):
            raise ValueError(
                f"excluded matrix combination '{text}' in '{path}' "
                "must be enclosed in brackets"
            )
        coordinate: dict[str, str] = {}
        for assignment in text[1:-1].split(","):
            if "=" not in assignment:
                raise ValueError(
                    f"excluded matrix combination '{text}' in '{path}' "
                    "must contain axis=value assignments"
                )
            axis, selected = (part.strip() for part in assignment.split("=", 1))
            if not axis or not selected:
                raise ValueError(
                    f"excluded matrix combination '{text}' in '{path}' "
                    "contains an empty axis or value"
                )
            if axis in coordinate:
                raise ValueError(
                    f"excluded matrix combination '{text}' in '{path}' "
                    f"repeats axis '{axis}'"
                )
            coordinate[axis] = selected

        unknown_axes = set(coordinate).difference(axes)
        missing_axes = set(axes).difference(coordinate)
        if unknown_axes or missing_axes:
            details = []
            if unknown_axes:
                details.append("unknown axes " + ", ".join(sorted(unknown_axes)))
            if missing_axes:
                details.append("missing axes " + ", ".join(sorted(missing_axes)))
            raise ValueError(
                f"excluded matrix combination '{text}' in '{path}' has "
                + " and ".join(details)
            )
        for axis, selected in coordinate.items():
            if selected not in axes[axis]:
                raise ValueError(
                    f"excluded matrix combination '{text}' in '{path}' selects "
                    f"unknown value '{selected}' on axis '{axis}'"
                )

        signature = tuple((axis, coordinate[axis]) for axis in axes)
        if signature in signatures:
            raise ValueError(
                f"excluded matrix combination '{text}' is duplicated in '{path}'"
            )
        signatures.add(signature)
        excluded.append(coordinate)
    return tuple(excluded)


def _matrix(
    values: dict[tuple[str, ...], str], path: Path
) -> tuple[dict[str, tuple[str, ...]], tuple[dict[str, str], ...]]:
    axes = {
        key[-1]: _split(value)
        for key, value in values.items()
        if len(key) == 2 and key[0] == "Matrix"
    }
    if any(not axis_values for axis_values in axes.values()):
        raise ValueError(f"matrix axes in '{path}' must have at least one value")

    selections = {
        key[-1]: _split(value)
        for key, value in values.items()
        if len(key) == 2
        and key[0] == "Selection"
        and key[-1] != "exclude combinations"
    }
    excluded = _excluded_combinations(
        _value(values, "Selection", "exclude combinations"), axes, path
    )
    resolved: dict[str, tuple[str, ...]] = {}
    for axis, axis_values in axes.items():
        selected = selections.get(axis, axis_values)
        unknown = set(selected).difference(axis_values)
        if unknown:
            raise ValueError(
                f"selection for matrix axis '{axis}' contains unknown values: "
                + ", ".join(sorted(unknown))
            )
        if not selected:
            raise ValueError(f"selection for matrix axis '{axis}' is empty")
        resolved[axis] = selected
    unknown_axes = set(selections).difference(axes)
    if unknown_axes:
        raise ValueError(
            "selection names undeclared matrix axes: "
            + ", ".join(sorted(unknown_axes))
        )

    axis_names = tuple(resolved)
    expanded = tuple(
        dict(zip(axis_names, combination))
        for combination in product(*(resolved[axis] for axis in axis_names))
    )
    excluded_signatures = {
        tuple((axis, coordinate[axis]) for axis in axes) for coordinate in excluded
    }
    combinations = tuple(
        combination
        for combination in expanded
        if tuple((axis, combination[axis]) for axis in axes)
        not in excluded_signatures
    )
    if not combinations:
        raise ValueError(f"matrix selection in '{path}' resolves to no combinations")
    return resolved, combinations


def _formats(value: str | None, path: Path) -> OutputFormats:
    formats = _split(value) or DEFAULT_OUTPUT_FORMATS
    unsupported = set(formats).difference({"png", "svg"})
    if unsupported:
        raise ValueError(
            f"unsupported post-processing formats in '{path}': "
            + ", ".join(sorted(unsupported))
        )
    if len(set(formats)) != len(formats):
        raise ValueError(f"post-processing formats in '{path}' must be unique")
    return tuple(formats)  # type: ignore[return-value]


def load_postprocess_configuration(path: Path) -> PostprocessConfiguration:
    """Load matrix and post-processing selections from a Deal.II-style file."""

    values = _read_values(path)
    matrix_axes, combinations = _matrix(values, path)
    fields = _split(
        _value(values, "Postprocessing", "fields")
        or _value(values, "Output", "selected fields")
    )

    rows = _comparison_axes(_value(values, "Postprocessing", "comparison rows"))
    columns = _comparison_axes(
        _value(values, "Postprocessing", "comparison columns")
    )
    group_by = _comparison_axes(
        _value(values, "Postprocessing", "comparison group by")
    )
    style_profile = _value(values, "Postprocessing", "style profile")
    if not style_profile:
        raise ValueError(f"'{path}' must declare a post-processing style profile")

    return PostprocessConfiguration(
        matrix_axes,
        combinations,
        fields,
        style_profile,
        ComparisonPlan(rows, columns, group_by),
        _formats(_value(values, "Postprocessing", "output formats"), path),
    )
