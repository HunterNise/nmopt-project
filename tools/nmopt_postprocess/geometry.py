"""Mesh geometry helpers shared by native-field renderers."""

import meshio
import numpy as np

from .errors import PostprocessError
from .fields import ScalarField


def volume_block(mesh: meshio.Mesh) -> tuple[int, meshio.CellBlock]:
    """Return the first two-dimensional cell block suitable for plotting."""

    for index, block in enumerate(mesh.cells):
        if block.type in {"triangle", "quad", "polygon"}:
            return index, block
    raise PostprocessError("volume field file contains no 2D cells")


def triangulate(block: meshio.CellBlock) -> tuple[np.ndarray, np.ndarray]:
    """Fan-triangulate polygonal cells and return their source cell indices."""

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


def boundary_field_block(
    mesh: meshio.Mesh, field: ScalarField
) -> tuple[int, meshio.CellBlock]:
    """Return the line cell block carrying a boundary cell field."""

    if field.location != "cell" or field.block_index is None:
        raise PostprocessError("boundary control must be cell data on line cells")
    block = mesh.cells[field.block_index]
    if block.type not in {"line", "line3"}:
        raise PostprocessError("boundary control is not attached to line cells")
    return field.block_index, block
