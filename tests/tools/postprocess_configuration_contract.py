#!/usr/bin/env python3
"""Focused contracts for .prm-authoritative post-processing configuration."""

from __future__ import annotations

import json
import shutil
import sys
import tempfile
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY_ROOT / "tools"))

from nmopt_postprocess.chapter6 import (
    CHAPTER6_PROFILE,
    load_json_profile,
    with_parameter_configuration,
)
from nmopt_postprocess.parameters import load_postprocess_configuration
from postprocess import _postprocess_provenance, _profile_from_run_snapshot
from nmopt_postprocess.pipeline import process_input_root


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> int:
    parameter_file = (
        REPOSITORY_ROOT
        / "parameters/chapter-6/b1/development/figure-6.3-book-policy.prm"
    )
    configuration = load_postprocess_configuration(parameter_file)
    require(
        configuration.matrix_axes["regularisation"] == ("1e-1", "1e-2", "1e-3"),
        "matrix values were not read from the experiment file",
    )
    require(
        len(configuration.matrix_combinations) == 6,
        "matrix combinations were not expanded from the experiment file",
    )
    require(configuration.output_formats == ("png",), "PNG is not the .prm format")
    require(
        configuration.comparison_plan.group_by == (),
        "comparison 'none' sentinel was treated as a matrix axis",
    )

    profile = load_json_profile(
        REPOSITORY_ROOT / "parameters/plotting/chapter-6-b1.json"
    )
    effective = with_parameter_configuration(profile, configuration)
    require(
        effective.comparison_plan == configuration.comparison_plan,
        "comparison plan did not come from the .prm file",
    )
    require(
        effective.matrix_axis_values == configuration.matrix_axes,
        "effective profile did not retain the active matrix axes",
    )
    require(
        {field.output_name for field in effective.volume_fields}
        == set(configuration.fields),
        "selected plotting fields did not come from the .prm file",
    )

    b2_configuration = load_postprocess_configuration(
        REPOSITORY_ROOT / "parameters/chapter-6/b2/authoritative.prm"
    )
    b2_profile = load_json_profile(
        REPOSITORY_ROOT / "parameters/plotting/chapter-6-b2.json"
    )
    b2_effective = with_parameter_configuration(b2_profile, b2_configuration)
    require(
        {field.output_name for field in b2_effective.volume_fields}
        == {
            "state",
            "state-no-control",
            "target",
            "forcing",
            "observation-region",
            "adjoint",
            "negative-adjoint",
        },
        "B2 volume field aliases were not resolved from the .prm file",
    )
    require(
        {field.output_name for field in b2_effective.boundary_fields}
        == {"control-boundary"},
        "B2 boundary field selection was not resolved from the .prm file",
    )
    forcing_configuration = load_postprocess_configuration(
        REPOSITORY_ROOT / "parameters/chapter-6/b2/development/forcing-sweep.prm"
    )
    require(
        forcing_configuration.comparison_plan.rows == ()
        and forcing_configuration.comparison_plan.columns == ("forcing",),
        "B2 one-row forcing comparison plan was not resolved",
    )

    with tempfile.TemporaryDirectory(prefix="nmopt-snapshot-") as directory:
        snapshot_root = Path(directory)
        shutil.copyfile(parameter_file, snapshot_root / "parameters.prm")
        shutil.copyfile(
            REPOSITORY_ROOT / "parameters/plotting/chapter-6-b1.json",
            snapshot_root / "plotting-profile.json",
        )
        shutil.copyfile(
            REPOSITORY_ROOT / "runs/chapter-6/b1/development/008/run-manifest.json",
            snapshot_root / "run-manifest.json",
        )
        snapshot_profile = _profile_from_run_snapshot(
            snapshot_root, CHAPTER6_PROFILE
        )
        require(
            snapshot_profile.output_formats == ("png",),
            "run snapshot did not use the .prm output format",
        )
        require(
            snapshot_profile.matrix_axis_values == configuration.matrix_axes,
            "run snapshot did not use the .prm matrix",
        )
        provenance = _postprocess_provenance(
            snapshot_root,
            snapshot_profile,
            ("png",),
            explicit_profile=False,
            profile_file=None,
            profile_name=None,
            requested_formats=None,
        )
        require(
            provenance["parameters"]["content_hash"]
            == "fnv1a64:c118220bc5bd1849",
            "parameter provenance hash was not retained",
        )
        require(
            provenance["plotting"]["content_hash"]
            == "fnv1a64:cdee2422f2c79582",
            "plotting provenance hash was not retained",
        )
        effective = provenance["effective"]
        require(
            effective["comparison"]
            == {"rows": ["method"], "columns": ["regularisation"], "group_by": []},
            "effective comparison plan was not recorded",
        )

        for method in ("steepest-descent", "l-bfgs"):
            for beta in ("1e-1", "1e-2", "1e-3"):
                artifact = snapshot_root / "artifacts" / method / f"beta-{beta}"
                artifact.mkdir(parents=True)
                (artifact / "artifact.kv").write_text(
                    "\n".join(
                        (
                            "identity.scenario_id=chapter-6.b1.distributed-laplace",
                            f"benchmark.method={method}",
                            f"benchmark.regularisation={beta}",
                            "solver.objective_history=0.05,0.04,0.03",
                            "solver.gradient_norm_history=0.01,0.001,0.0001",
                        )
                    )
                    + "\n",
                    encoding="utf-8",
                )
        index = process_input_root(
            snapshot_root,
            snapshot_root / "postprocess-contract",
            snapshot_profile,
            output_formats=("png",),
            provenance=provenance,
        )
        require(
            index["provenance"] == provenance,
            "aggregate post-processing provenance was not written",
        )
        artifact_manifest = json.loads(
            (
                snapshot_root
                / "artifacts/steepest-descent/beta-1e-1/postprocess/postprocess.json"
            ).read_text(encoding="utf-8")
        )
        require(
            artifact_manifest["provenance"] == provenance,
            "artifact post-processing provenance was not written",
        )

    with tempfile.TemporaryDirectory(prefix="nmopt-prm-") as directory:
        path = Path(directory) / "selection.prm"
        path.write_text(
            """subsection Matrix
  set method = steepest-descent, l-bfgs
  set regularisation = 1e-1, 1e-2
end
subsection Selection
  set method = l-bfgs
end
subsection Output
  set selected fields = state
end
subsection Postprocessing
  set style profile = style.json
  set comparison rows = method
  set comparison columns = regularisation
  set comparison group by = none
  set output formats = svg
end
""",
            encoding="utf-8",
        )
        selected = load_postprocess_configuration(path)
        require(
            selected.matrix_combinations
            == ({"method": "l-bfgs", "regularisation": "1e-1"},
                {"method": "l-bfgs", "regularisation": "1e-2"}),
            "matrix selection was not applied before plotting",
        )
        require(selected.output_formats == ("svg",), "explicit SVG format was lost")

    print("postprocess configuration contract passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
