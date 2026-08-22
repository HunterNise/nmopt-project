#!/usr/bin/env python3
"""Focused contracts for grouped and validated comparison matrices."""

from __future__ import annotations

import json
import sys
import tempfile
from dataclasses import replace
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY_ROOT / "tools"))

from nmopt_postprocess.chapter6 import load_json_profile, with_parameter_configuration
from nmopt_postprocess.errors import PostprocessError
from nmopt_postprocess.parameters import load_postprocess_configuration
from nmopt_postprocess.pipeline import (
    ComparisonPlan,
    build_comparisons,
    comparison_grid,
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def write_artifact(path: Path, method: str, beta: str, forcing: str) -> None:
    path.mkdir(parents=True)
    (path / "artifact.kv").write_text(
        "\n".join(
            (
                "identity.scenario_id=chapter-6.b1.distributed-laplace",
                f"benchmark.method={method}",
                f"benchmark.regularisation={beta}",
                f"benchmark.forcing={forcing}",
                "solver.objective_history=0.05,0.04,0.03",
                "solver.gradient_norm_history=0.01,0.001,0.0001",
            )
        )
        + "\n",
        encoding="utf-8",
    )


def main() -> int:
    profile_template = load_json_profile(
        REPOSITORY_ROOT / "parameters/plotting/chapter-6-b1.json"
    )
    with tempfile.TemporaryDirectory(prefix="nmopt-comparison-plan-") as directory:
        root = Path(directory)
        parameter_file = root / "comparison.prm"
        parameter_file.write_text(
            """subsection Matrix
  set method = steepest-descent, l-bfgs
  set regularisation = 1e-1, 1e-2
  set forcing = zero, one
end
subsection Postprocessing
  set fields = state
  set style profile = chapter-6-b1.json
  set comparison rows = method
  set comparison columns = regularisation
  set comparison group by = forcing
  set output formats = png
end
""",
            encoding="utf-8",
        )
        configuration = load_postprocess_configuration(parameter_file)
        profile = with_parameter_configuration(profile_template, configuration)

        artifacts: list[Path] = []
        for forcing in ("zero", "one"):
            for method in ("steepest-descent", "l-bfgs"):
                for beta in ("1e-1", "1e-2"):
                    artifact = (
                        root
                        / "artifacts"
                        / forcing
                        / method
                        / f"beta-{beta}"
                    )
                    write_artifact(artifact, method, beta, forcing)
                    artifacts.append(artifact)

        generated, errors = build_comparisons(
            artifacts, root / "postprocess", profile
        )
        require(not errors, f"grouped comparison failed: {json.dumps(errors)}")
        require(
            set(generated) == {
                "chapter-6.b1.distributed-laplace/forcing=zero",
                "chapter-6.b1.distributed-laplace/forcing=one",
            },
            "group_by did not partition the comparison into one grid per forcing",
        )

        zero_artifacts = [
            artifact
            for artifact in artifacts
            if "forcing=zero" not in str(artifact)
            and "artifacts/zero/" in str(artifact)
        ]
        require(
            comparison_grid(zero_artifacts, profile).rows == 2,
            "grouped comparison did not retain method rows",
        )
        require(
            comparison_grid(zero_artifacts, profile).columns == 2,
            "grouped comparison did not retain regularisation columns",
        )
        legacy = replace(
            profile,
            comparison_plan=ComparisonPlan(),
            matrix_axis_values=None,
            matrix_combinations=(),
        )
        legacy_grid = comparison_grid(zero_artifacts, legacy)
        require(
            (legacy_grid.rows, legacy_grid.columns) == (1, 4),
            "legacy profile did not retain its one-row fallback",
        )

        try:
            comparison_grid(zero_artifacts[:-1], profile)
        except PostprocessError as error:
            require("incomplete" in str(error), "missing coordinate was not reported")
        else:
            raise RuntimeError("incomplete comparison matrix was accepted")

        try:
            comparison_grid(zero_artifacts + [zero_artifacts[0]], profile)
        except PostprocessError as error:
            require("duplicate coordinates" in str(error), "duplicate was not reported")
        else:
            raise RuntimeError("duplicate comparison coordinate was accepted")

        ungrouped = replace(
            profile,
            comparison_plan=ComparisonPlan(
                rows=("method",), columns=("regularisation",), group_by=()
            ),
        )
        try:
            comparison_grid(artifacts, ungrouped)
        except PostprocessError as error:
            require(
                "varying matrix axes" in str(error),
                "unbound matrix axis was not reported",
            )
        else:
            raise RuntimeError("unbound matrix axis was accepted")

    print("postprocess comparison-plan contract passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
