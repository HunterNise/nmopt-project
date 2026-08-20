"""Headless Matplotlib renderers for scalar native fields."""

from __future__ import annotations

import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

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

from .errors import PostprocessError
from .fields import ScalarField
from .geometry import boundary_field_block, triangulate, volume_block


OutputFormat = Literal["png", "svg"]
OutputFormats = tuple[OutputFormat, ...]
DEFAULT_OUTPUT_FORMATS: OutputFormats = ("png",)
# The comparison policy requested for the Chapter 6 reproduction is the
# perceptually smoother ``turbo`` map, applied consistently to individual and
# comparison panels.
BOOK_COLORMAP = "turbo"
COLORBAR_TICK_COUNT = 5


@dataclass(frozen=True)
class RenderItem:
    """A field, its mesh, and the label used in a comparison panel."""

    label: str
    mesh: meshio.Mesh
    field: ScalarField


def field_norm(values: np.ndarray) -> Normalize:
    """Create a finite-value scale whose limits are the field extrema."""

    finite_values = values[np.isfinite(values)]
    if len(finite_values) == 0:
        raise PostprocessError("field contains no finite values")
    lower = float(np.min(finite_values))
    upper = float(np.max(finite_values))
    if lower == upper:
        padding = max(abs(lower) * 0.05, 1.0e-12)
        lower -= padding
        upper += padding
    return Normalize(vmin=lower, vmax=upper)


def colorbar_ticks(norm: Normalize) -> np.ndarray:
    """Include both extrema while retaining a compact interior scale."""

    return np.linspace(norm.vmin, norm.vmax, COLORBAR_TICK_COUNT)


def save_figure(
    figure: plt.Figure,
    output: Path,
    output_formats: OutputFormats = DEFAULT_OUTPUT_FORMATS,
) -> list[str]:
    """Write the selected figure formats, then close Matplotlib resources."""

    generated: list[str] = []
    for output_format in output_formats:
        path = output.with_suffix(f".{output_format}")
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
    """Draw a scalar point or cell field on a two-dimensional mesh."""

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
            cmap=BOOK_COLORMAP,
            norm=norm,
        )
    else:
        image = axis.tripcolor(
            triangulation,
            facecolors=field.values[owners],
            shading="flat",
            cmap=BOOK_COLORMAP,
            norm=norm,
        )
    # Keep scalar-field renders free of a mesh overlay so the field can be
    # compared with the source raster figures. The native VTU remains the
    # authoritative place to inspect the realized topology.
    axis.set_aspect("equal")
    axis.set_xlabel("x")
    axis.set_ylabel("y")
    return image


def plot_volume_field(
    mesh: meshio.Mesh,
    field: ScalarField,
    title: str,
    output: Path,
    colorbar_label: str | None = None,
    output_formats: OutputFormats = DEFAULT_OUTPUT_FORMATS,
) -> list[str]:
    """Render one scalar volume field in the selected formats."""

    figure, axis = plt.subplots(figsize=(7, 5.5), constrained_layout=True)
    norm = field_norm(field.values)
    image = draw_volume_field(axis, mesh, field, norm)
    axis.set_title(title)
    figure.colorbar(
        image,
        ax=axis,
        ticks=colorbar_ticks(norm),
        label=colorbar_label or field.name,
    )
    return save_figure(figure, output, output_formats)


def draw_boundary_field(
    axis: plt.Axes,
    mesh: meshio.Mesh,
    field: ScalarField,
    norm: Normalize | None = None,
) -> object:
    """Draw a scalar field stored on boundary line cells."""

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
        cmap=BOOK_COLORMAP,
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
    title: str,
    output: Path,
    colorbar_label: str | None = None,
    output_formats: OutputFormats = DEFAULT_OUTPUT_FORMATS,
) -> list[str]:
    """Render one scalar boundary field in the selected formats."""

    figure, axis = plt.subplots(figsize=(7, 5.5), constrained_layout=True)
    norm = field_norm(field.values)
    collection = draw_boundary_field(axis, mesh, field, norm)
    axis.set_title(title)
    figure.colorbar(
        collection,
        ax=axis,
        ticks=colorbar_ticks(norm),
        label=colorbar_label or field.name,
    )
    return save_figure(figure, output, output_formats)


def comparison_norm(items: list[RenderItem]) -> Normalize:
    """Create one finite-value color scale shared by comparison panels."""

    finite_values = [item.field.values[np.isfinite(item.field.values)] for item in items]
    values = np.concatenate([array for array in finite_values if len(array)])
    if len(values) == 0:
        raise PostprocessError("comparison field contains no finite values")
    return field_norm(values)


def comparison_layout(count: int) -> tuple[int, int]:
    """Choose a compact comparison grid for the active case matrix."""

    if count == 8:
        return 2, 4
    if count == 4:
        return 2, 2
    columns = min(3, count)
    rows = (count + columns - 1) // columns
    return rows, columns


def plot_volume_comparison(
    items: list[RenderItem],
    field_name: str,
    output: Path,
    title: str,
    colorbar_label: str | None = None,
    output_formats: OutputFormats = DEFAULT_OUTPUT_FORMATS,
) -> list[str]:
    """Render volume fields with one color scale across all panels."""

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
    for index, item in enumerate(items):
        axis = axes.flat[index]
        image = draw_volume_field(axis, item.mesh, item.field, norm)
        axis.set_title(item.label, fontsize=8)
    for axis in axes.flat[len(items) :]:
        axis.axis("off")
    if image is None:
        raise PostprocessError(f"no volume fields available for '{field_name}'")
    figure.suptitle(title)
    figure.colorbar(
        image,
        ax=axes.ravel().tolist(),
        ticks=colorbar_ticks(norm),
        label=colorbar_label or field_name,
    )
    return save_figure(figure, output, output_formats)


def plot_boundary_comparison(
    items: list[RenderItem],
    output: Path,
    title: str,
    colorbar_label: str | None = None,
    output_formats: OutputFormats = DEFAULT_OUTPUT_FORMATS,
) -> list[str]:
    """Render boundary fields with one color scale across all panels."""

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
    for index, item in enumerate(items):
        axis = axes.flat[index]
        collection = draw_boundary_field(axis, item.mesh, item.field, norm)
        axis.set_title(item.label, fontsize=8)
    for axis in axes.flat[len(items) :]:
        axis.axis("off")
    if collection is None:
        raise PostprocessError("no boundary fields available for comparison")
    figure.suptitle(title)
    figure.colorbar(
        collection,
        ax=axes.ravel().tolist(),
        ticks=colorbar_ticks(norm),
        label=colorbar_label or items[0].field.name,
    )
    return save_figure(figure, output, output_formats)
