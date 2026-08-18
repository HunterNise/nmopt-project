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
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(tempfile.gettempdir()) / "nmopt-matplotlib"),
)

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import meshio
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.tri import Triangulation


class PostprocessError(RuntimeError):
    """An input artifact cannot be rendered by this post-processor."""


@dataclass(frozen=True)
class ScalarField:
    name: str
    values: np.ndarray
    location: str
    block_index: int | None


def unescape(value: str) -> str:
    decoded: list[str] = []
    index = 0
    escaped = {"\\": "\\", "n": "\n", "r": "\r", "t": "\t", "=": "="}
    while index < len(value):
        if value[index] == "\\" and index + 1 < len(value):
            decoded.append(escaped.get(value[index + 1], value[index + 1]))
            index += 2
        else:
            decoded.append(value[index])
            index += 1
    return "".join(decoded)


def read_metadata(artifact: Path) -> dict[str, str]:
    path = artifact / "artifact.kv"
    if not path.is_file():
        return {}
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line and "=" in line:
            key, value = line.split("=", 1)
            values[key] = unescape(value)
    return values


def first_existing(artifact: Path, names: Iterable[str]) -> Path | None:
    candidates = [artifact / "native" / name for name in names]
    candidates.extend(artifact / name for name in names)
    return next((path for path in candidates if path.is_file()), None)


def scalar_values(values: object, name: str) -> np.ndarray:
    array = np.asarray(values)
    if array.ndim == 1:
        return array.astype(float, copy=False)
    if array.ndim == 2 and array.shape[1] == 1:
        return array[:, 0].astype(float, copy=False)
    raise PostprocessError(f"field '{name}' is not scalar")


def find_field(mesh: meshio.Mesh, name: str) -> ScalarField | None:
    if name in mesh.point_data:
        return ScalarField(
            name,
            scalar_values(mesh.point_data[name], name),
            "point",
            None,
        )

    for block_index, arrays in enumerate(mesh.cell_data.get(name, [])):
        return ScalarField(
            name,
            scalar_values(arrays, name),
            "cell",
            block_index,
        )
    return None


def volume_block(mesh: meshio.Mesh) -> tuple[int, meshio.CellBlock]:
    for index, block in enumerate(mesh.cells):
        if block.type in {"triangle", "quad", "polygon"}:
            return index, block
    raise PostprocessError("volume field file contains no 2D cells")


def triangulate(block: meshio.CellBlock) -> tuple[np.ndarray, np.ndarray]:
    triangles: list[list[int]] = []
    owners: list[int] = []
    for cell_index, cell in enumerate(block.data):
        if len(cell) < 3:
            continue
        for vertex in range(1, len(cell) - 1):
            triangles.append([int(cell[0]), int(cell[vertex]), int(cell[vertex + 1])])
            owners.append(cell_index)
    if not triangles:
        raise PostprocessError("volume mesh contains no renderable triangles")
    return np.asarray(triangles, dtype=int), np.asarray(owners, dtype=int)


def axes_title(metadata: dict[str, str], field_name: str) -> str:
    output_id = metadata.get("identity.output_id", "Chapter 6 artifact")
    return f"{output_id}: {field_name}"


def save_figure(figure: plt.Figure, output: Path) -> list[str]:
    generated: list[str] = []
    for suffix in (".png", ".svg"):
        path = output.with_suffix(suffix)
        figure.savefig(path, dpi=180, bbox_inches="tight")
        generated.append(path.name)
    plt.close(figure)
    return generated


def draw_volume_field(
    axis: plt.Axes,
    mesh: meshio.Mesh,
    field: ScalarField,
    norm: Normalize | None = None,
) -> object:
    block_index, block = volume_block(mesh)
    if field.location == "cell" and field.block_index != block_index:
        raise PostprocessError(
            f"field '{field.name}' is attached to a non-volume cell block"
        )

    triangles, owners = triangulate(block)
    points = np.asarray(mesh.points)[:, :2]
    if field.location == "point":
        if len(field.values) != len(points):
            raise PostprocessError(
                f"point field '{field.name}' has {len(field.values)} values, "
                f"expected {len(points)}"
            )
    else:
        if len(field.values) != len(block.data):
            raise PostprocessError(
                f"cell field '{field.name}' has {len(field.values)} values, "
                f"expected {len(block.data)}"
            )

    triangulation = Triangulation(points[:, 0], points[:, 1], triangles)
    if field.location == "point":
        image = axis.tripcolor(
            triangulation,
            field.values,
            shading="gouraud",
            cmap="viridis",
            norm=norm,
        )
    else:
        image = axis.tripcolor(
            triangulation,
            facecolors=field.values[owners],
            shading="flat",
            cmap="viridis",
            norm=norm,
        )
    axis.triplot(triangulation, color="black", linewidth=0.25, alpha=0.25)
    axis.set_aspect("equal")
    axis.set_xlabel("x")
    axis.set_ylabel("y")
    return image


def plot_volume_field(
    mesh: meshio.Mesh,
    field: ScalarField,
    metadata: dict[str, str],
    output: Path,
) -> list[str]:
    figure, axis = plt.subplots(figsize=(7, 5.5), constrained_layout=True)
    image = draw_volume_field(axis, mesh, field)
    axis.set_title(axes_title(metadata, field.name))
    figure.colorbar(image, ax=axis, label=field.name)
    return save_figure(figure, output)


def boundary_field_block(
    mesh: meshio.Mesh, field: ScalarField
) -> tuple[int, meshio.CellBlock]:
    if field.location != "cell" or field.block_index is None:
        raise PostprocessError("boundary control must be cell data on line cells")
    block = mesh.cells[field.block_index]
    if block.type not in {"line", "line3"}:
        raise PostprocessError("boundary control is not attached to line cells")
    return field.block_index, block


def draw_boundary_field(
    axis: plt.Axes,
    mesh: meshio.Mesh,
    field: ScalarField,
    norm: Normalize | None = None,
) -> object:
    _, block = boundary_field_block(mesh, field)
    if len(field.values) != len(block.data):
        raise PostprocessError(
            f"boundary field has {len(field.values)} values, "
            f"expected {len(block.data)}"
        )
    points = np.asarray(mesh.points)[:, :2]
    segments = points[block.data][:, :, :2]
    collection = LineCollection(
        segments,
        array=field.values,
        cmap="viridis",
        linewidths=3.0,
        norm=norm,
    )
    axis.add_collection(collection)
    axis.autoscale()
    axis.set_aspect("equal")
    axis.set_xlabel("x")
    axis.set_ylabel("y")
    return collection


def plot_boundary_field(
    mesh: meshio.Mesh,
    field: ScalarField,
    metadata: dict[str, str],
    output: Path,
) -> list[str]:
    figure, axis = plt.subplots(figsize=(7, 5.5), constrained_layout=True)
    collection = draw_boundary_field(axis, mesh, field)
    axis.set_title(axes_title(metadata, "control-boundary"))
    figure.colorbar(collection, ax=axis, label="control")
    return save_figure(figure, output)


def comparison_norm(
    items: list[tuple[Path, dict[str, str], meshio.Mesh, ScalarField]]
) -> Normalize:
    finite_values = [
        item[3].values[np.isfinite(item[3].values)] for item in items
    ]
    values = np.concatenate([array for array in finite_values if len(array)])
    if len(values) == 0:
        raise PostprocessError("comparison field contains no finite values")
    lower = float(np.min(values))
    upper = float(np.max(values))
    if lower == upper:
        padding = max(abs(lower) * 0.05, 1.0e-12)
        lower -= padding
        upper += padding
    return Normalize(vmin=lower, vmax=upper)


def comparison_layout(count: int) -> tuple[int, int]:
    columns = min(3, count)
    rows = (count + columns - 1) // columns
    return rows, columns


def plot_volume_comparison(
    items: list[tuple[Path, dict[str, str], meshio.Mesh, ScalarField]],
    field_name: str,
    output: Path,
) -> list[str]:
    norm = comparison_norm(items)
    rows, columns = comparison_layout(len(items))
    figure, axes = plt.subplots(
        rows,
        columns,
        figsize=(6.0 * columns, 4.8 * rows),
        squeeze=False,
        constrained_layout=True,
    )
    image = None
    for index, (_, metadata, mesh, field) in enumerate(items):
        axis = axes.flat[index]
        image = draw_volume_field(axis, mesh, field, norm)
        axis.set_title(
            metadata.get("identity.output_id", "Chapter 6 artifact"),
            fontsize=8,
        )
    for axis in axes.flat[len(items) :]:
        axis.axis("off")
    if image is None:
        raise PostprocessError(f"no volume fields available for '{field_name}'")
    figure.suptitle(f"Chapter 6 comparison: {field_name}")
    figure.colorbar(image, ax=axes.ravel().tolist(), label=field_name)
    return save_figure(figure, output)


def plot_boundary_comparison(
    items: list[tuple[Path, dict[str, str], meshio.Mesh, ScalarField]],
    output: Path,
) -> list[str]:
    norm = comparison_norm(items)
    rows, columns = comparison_layout(len(items))
    figure, axes = plt.subplots(
        rows,
        columns,
        figsize=(6.0 * columns, 4.8 * rows),
        squeeze=False,
        constrained_layout=True,
    )
    collection = None
    for index, (_, metadata, mesh, field) in enumerate(items):
        axis = axes.flat[index]
        collection = draw_boundary_field(axis, mesh, field, norm)
        axis.set_title(
            metadata.get("identity.output_id", "Chapter 6 artifact"),
            fontsize=8,
        )
    for axis in axes.flat[len(items) :]:
        axis.axis("off")
    if collection is None:
        raise PostprocessError("no boundary fields available for comparison")
    figure.suptitle("Chapter 6 comparison: control-boundary")
    figure.colorbar(collection, ax=axes.ravel().tolist(), label="control")
    return save_figure(figure, output)


def volume_comparison_items(
    artifacts: list[Path], field_name: str
) -> list[tuple[Path, dict[str, str], meshio.Mesh, ScalarField]]:
    items: list[tuple[Path, dict[str, str], meshio.Mesh, ScalarField]] = []
    for artifact in artifacts:
        fields_path = first_existing(artifact, ("fields-volume.vtu", "fields.vtu"))
        if fields_path is None:
            continue
        mesh = meshio.read(fields_path)
        field = find_field(mesh, field_name)
        if field is not None:
            items.append((artifact, read_metadata(artifact), mesh, field))
    return items


def boundary_comparison_items(
    artifacts: list[Path],
) -> list[tuple[Path, dict[str, str], meshio.Mesh, ScalarField]]:
    items: list[tuple[Path, dict[str, str], meshio.Mesh, ScalarField]] = []
    for artifact in artifacts:
        boundary_path = first_existing(
            artifact, ("control-boundary.vtu", "control.vtu")
        )
        if boundary_path is None:
            continue
        mesh = meshio.read(boundary_path)
        field = find_field(mesh, "control")
        if field is not None:
            items.append((artifact, read_metadata(artifact), mesh, field))
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
                    )
            except (OSError, ValueError, PostprocessError) as error:
                group_errors[field_name] = str(error)

        try:
            items = boundary_comparison_items(group_artifacts)
            if len(items) >= 2:
                group_generated["control-boundary"] = plot_boundary_comparison(
                    items,
                    comparison_dir / "comparison-control-boundary",
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
                metadata,
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
                metadata,
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
