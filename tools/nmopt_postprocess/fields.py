"""Extract scalar point and cell fields from meshio meshes."""

from dataclasses import dataclass

import meshio
import numpy as np

from .errors import PostprocessError


@dataclass(frozen=True)
class ScalarField:
    """A scalar mesh field together with its storage location."""

    name: str
    values: np.ndarray
    location: str
    block_index: int | None


def scalar_values(values: object, name: str) -> np.ndarray:
    """Convert a one-component point or cell array to a flat float array."""

    array = np.asarray(values)
    if array.ndim == 1:
        return array.astype(float, copy=False)
    if array.ndim == 2 and array.shape[1] == 1:
        return array[:, 0].astype(float, copy=False)
    raise PostprocessError(f"field '{name}' is not scalar")


def find_field(mesh: meshio.Mesh, name: str) -> ScalarField | None:
    """Find a named scalar point field or the first matching cell field."""

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
