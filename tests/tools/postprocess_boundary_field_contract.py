#!/usr/bin/env python3
"""Focused contract for facewise and continuous boundary-field rendering."""

from __future__ import annotations

import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import meshio
import numpy as np
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY_ROOT / "tools"))

from nmopt_postprocess.fields import find_field
from nmopt_postprocess.render import draw_boundary_field, field_norm


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> int:
    points = np.asarray(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
        dtype=float,
    )
    cells = [meshio.CellBlock("line", np.asarray([[0, 1], [1, 2]], dtype=int))]

    continuous_mesh = meshio.Mesh(
        points,
        cells,
        point_data={"control": np.asarray([0.0, 2.0, 4.0])},
    )
    continuous_field = find_field(continuous_mesh, "control")
    require(
        continuous_field is not None and continuous_field.location == "point",
        "continuous control was not detected as point data",
    )
    figure, axis = plt.subplots()
    collection = draw_boundary_field(
        axis,
        continuous_mesh,
        continuous_field,
        field_norm(continuous_field.values),
    )
    np.testing.assert_allclose(collection.get_array(), [1.0, 3.0])
    plt.close(figure)

    facewise_mesh = meshio.Mesh(
        points,
        cells,
        cell_data={"control": [np.asarray([5.0, 7.0])]},
    )
    facewise_field = find_field(facewise_mesh, "control")
    require(
        facewise_field is not None and facewise_field.location == "cell",
        "facewise control was not detected as cell data",
    )
    figure, axis = plt.subplots()
    collection = draw_boundary_field(
        axis,
        facewise_mesh,
        facewise_field,
        field_norm(facewise_field.values),
    )
    np.testing.assert_allclose(collection.get_array(), [5.0, 7.0])
    plt.close(figure)

    print("postprocess boundary-field contract passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
