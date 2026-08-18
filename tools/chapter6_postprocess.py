#!/usr/bin/env python3
"""Render Chapter 6 native fields with meshio and matplotlib.

The tool reads one runner artifact directory and writes derived plots outside
the authoritative ``native/`` directory.  It accepts the current native file
names and the historical ``fields.vtu``/``control.vtu`` names used by older
development runs.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import meshio
import numpy as np
from matplotlib.collections import LineCollection
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


def plot_volume_field(
    mesh: meshio.Mesh,
    field: ScalarField,
    metadata: dict[str, str],
    output: Path,
) -> list[str]:
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
    figure, axis = plt.subplots(figsize=(7, 5.5), constrained_layout=True)
    if field.location == "point":
        image = axis.tripcolor(
            triangulation,
            field.values,
            shading="gouraud",
            cmap="viridis",
        )
    else:
        image = axis.tripcolor(
            triangulation,
            facecolors=field.values[owners],
            shading="flat",
            cmap="viridis",
        )
    axis.triplot(triangulation, color="black", linewidth=0.25, alpha=0.25)
    axis.set_aspect("equal")
    axis.set_xlabel("x")
    axis.set_ylabel("y")
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


def plot_boundary_field(
    mesh: meshio.Mesh,
    field: ScalarField,
    metadata: dict[str, str],
    output: Path,
) -> list[str]:
    _, block = boundary_field_block(mesh, field)
    if len(field.values) != len(block.data):
        raise PostprocessError(
            f"boundary field has {len(field.values)} values, "
            f"expected {len(block.data)}"
        )
    points = np.asarray(mesh.points)[:, :2]
    segments = points[block.data][:, :, :2]
    figure, axis = plt.subplots(figsize=(7, 5.5), constrained_layout=True)
    collection = LineCollection(
        segments,
        array=field.values,
        cmap="viridis",
        linewidths=3.0,
    )
    axis.add_collection(collection)
    axis.autoscale()
    axis.set_aspect("equal")
    axis.set_xlabel("x")
    axis.set_ylabel("y")
    axis.set_title(axes_title(metadata, "control-boundary"))
    figure.colorbar(collection, ax=axis, label="control")
    return save_figure(figure, output)


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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--artifact",
        type=Path,
        required=True,
        help="one runner artifact directory containing artifact.kv/native",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="derived output directory (default: <artifact>/derived)",
    )
    arguments = parser.parse_args()
    artifact = arguments.artifact.resolve()
    output = (
        arguments.output.resolve()
        if arguments.output is not None
        else artifact / "derived"
    )
    try:
        manifest = process_artifact(artifact, output)
    except (OSError, ValueError, PostprocessError) as error:
        parser.error(str(error))
    print(f"wrote {len(manifest['plots'])} field plot sets to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
