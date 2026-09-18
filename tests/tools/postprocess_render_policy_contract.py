#!/usr/bin/env python3
"""Focused contracts for profile-driven rendering policy and axis lookup."""

from __future__ import annotations

import json
import sys
import tempfile
from dataclasses import replace
from pathlib import Path
from unittest.mock import Mock

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import meshio
import numpy as np
from matplotlib.collections import LineCollection


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY_ROOT / "tools"))

from nmopt_postprocess.chapter6 import load_json_profile
from nmopt_postprocess.fields import find_field
from nmopt_postprocess.pipeline import comparison_grid
from nmopt_postprocess.render import (
    colorbar_ticks,
    draw_volume_field,
    field_norm,
    save_figure,
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> int:
    profile_path = REPOSITORY_ROOT / "parameters/plotting/chapter-6-b1.json"
    profile = load_json_profile(profile_path)
    require(profile.render_policy.colormap == "turbo", "B1 colormap changed")
    require(
        profile.render_policy.normalization == "finite-extrema",
        "B1 normalization changed",
    )
    require(
        profile.render_policy.comparison_normalization == "shared-finite-extrema",
        "B1 comparison normalization changed",
    )
    require(
        profile.render_policy.volume_interpolation == "gouraud",
        "B1 volume interpolation changed",
    )
    require(
        profile.render_policy.volume_mesh_overlay is False,
        "B1 mesh overlay changed",
    )
    require(
        profile.render_policy.colorbar_ticks == "endpoint-inclusive"
        and profile.render_policy.colorbar_tick_count == 5,
        "B1 colorbar policy changed",
    )
    require(profile.render_policy.dpi == 180, "B1 dpi changed")
    require(
        profile.render_policy.coordinate_axis_labels == ("x", "y"),
        "B1 axis labels changed",
    )
    require(profile.output_formats == ("png",), "B1 output format changed")

    with tempfile.TemporaryDirectory(prefix="nmopt-render-policy-") as directory:
        root = Path(directory)
        document = json.loads(profile_path.read_text(encoding="utf-8"))
        document["defaults"].update(
            {
                "colormap": "viridis",
                "volume_mesh_overlay": True,
                "colorbar_tick_count": 3,
                "dpi": 72,
                "axis_labels": ["X1", "X2"],
                "output_formats": ["svg"],
            }
        )
        modified_path = root / "modified.json"
        modified_path.write_text(json.dumps(document), encoding="utf-8")
        modified = load_json_profile(modified_path)
        policy = modified.render_policy
        require(policy.colormap == "viridis", "modified colormap was ignored")
        require(policy.volume_mesh_overlay, "modified overlay was ignored")
        require(policy.colorbar_tick_count == 3, "modified tick count was ignored")
        require(policy.dpi == 72, "modified dpi was ignored")
        require(
            policy.coordinate_axis_labels == ("X1", "X2"),
            "modified axis labels were ignored",
        )
        require(modified.output_formats == ("svg",), "modified format was ignored")

        points = np.asarray(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
            dtype=float,
        )
        mesh = meshio.Mesh(
            points,
            [meshio.CellBlock("quad", np.asarray([[0, 1, 2, 3]], dtype=int))],
            point_data={"state": np.asarray([0.0, 1.0, 2.0, 3.0])},
        )
        field = find_field(mesh, "state")
        require(field is not None, "synthetic volume field was not found")
        figure, axis = plt.subplots()
        norm = field_norm(field.values, policy)
        image = draw_volume_field(axis, mesh, field, norm, policy)
        require(image.get_cmap().name == "viridis", "policy colormap was not rendered")
        require(axis.get_xlabel() == "X1" and axis.get_ylabel() == "X2", "axis labels were not rendered")
        overlays = [
            collection
            for collection in axis.collections
            if isinstance(collection, LineCollection)
        ]
        require(len(overlays) == 1, "volume mesh overlay was not rendered")
        require(
            len(overlays[0].get_segments()) == 4,
            "volume overlay did not retain the original cell edges",
        )
        ticks = colorbar_ticks(norm, policy)
        require(len(ticks) == 3, "policy tick count was not rendered")
        np.testing.assert_allclose(ticks[[0, -1]], [norm.vmin, norm.vmax])

        figure.savefig = Mock()
        generated = save_figure(
            figure, root / "rendered", modified.output_formats, policy
        )
        require(generated == ["rendered.svg"], "policy output format was not saved")
        figure.savefig.assert_called_once()
        require(
            figure.savefig.call_args.kwargs["dpi"] == 72,
            "policy dpi was not passed to savefig",
        )

        unsupported = {
            "normalization": "unsupported",
            "comparison_normalization": "unsupported",
            "volume_interpolation": "unsupported",
            "colorbar_ticks": "unsupported",
        }
        for key, value in unsupported.items():
            invalid_document = json.loads(json.dumps(document))
            invalid_document["defaults"][key] = value
            invalid_path = root / f"invalid-{key}.json"
            invalid_path.write_text(
                json.dumps(invalid_document), encoding="utf-8"
            )
            try:
                load_json_profile(invalid_path)
            except ValueError:
                pass
            else:
                raise RuntimeError(f"unsupported {key} was accepted")

        legacy_profile = replace(
            profile,
            axis_value=lambda metadata, axis: metadata.get(f"legacy.{axis}", ""),
            matrix_axis_values=None,
            matrix_combinations=(),
        )
        artifacts = []
        for method in ("steepest-descent", "l-bfgs"):
            for regularisation in ("1e-1", "1e-2"):
                artifact = root / "artifacts" / method / regularisation
                artifact.mkdir(parents=True)
                (artifact / "artifact.kv").write_text(
                    "\n".join(
                        (
                            "identity.scenario_id=legacy-scenario",
                            f"legacy.method={method}",
                            f"legacy.regularisation={regularisation}",
                        )
                    )
                    + "\n",
                    encoding="utf-8",
                )
                artifacts.append(artifact)
        grid = comparison_grid(artifacts, legacy_profile)
        require(
            (grid.rows, grid.columns) == (2, 2),
            "profile axis resolver was not used for comparison coordinates",
        )

        plt.close(figure)

    print("postprocess render-policy contract passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
