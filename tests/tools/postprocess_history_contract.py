#!/usr/bin/env python3
"""Focused contract checks for profile-driven optimization-history plots."""

from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY_ROOT / "tools"))

from nmopt_postprocess.chapter6 import load_json_profile
from nmopt_postprocess.pipeline import build_comparisons, comparison_dimensions
from nmopt_postprocess.records import read_numeric_history
from nmopt_postprocess.render import comparison_layout


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def write_artifact(path: Path, method: str, beta: str) -> None:
    path.mkdir(parents=True)
    (path / "artifact.kv").write_text(
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


def main() -> int:
    profile = load_json_profile(
        REPOSITORY_ROOT / "parameters/plotting/chapter-6-b1.json"
    )
    require(len(profile.history_figures) == 1, "B1 needs one history figure")
    require(
        len(profile.history_figures[0].plots) == 2,
        "Figure 6.3 needs objective and gradient panels",
    )
    require(
        read_numeric_history({"history": "1,2.5,3"}, "history")
        == (1.0, 2.5, 3.0),
        "numeric history parsing changed",
    )

    with tempfile.TemporaryDirectory(prefix="nmopt-history-") as directory:
        root = Path(directory)
        artifacts = []
        for method in ("steepest-descent", "l-bfgs"):
            for beta in ("1e-1", "1e-2", "1e-3"):
                artifact = root / "artifacts" / method / f"beta-{beta}"
                write_artifact(artifact, method, beta)
                artifacts.append(artifact)
        require(
            comparison_dimensions(artifacts, profile) == (2, 3),
            "active three-beta matrix must produce a 2x3 comparison grid",
        )
        require(
            comparison_layout(6) == (2, 3),
            "fallback comparison layout should remain compact and generic",
        )
        generated, errors = build_comparisons(
            artifacts, root / "postprocess", profile
        )
        require(not errors, f"history rendering failed: {json.dumps(errors)}")
        outputs = generated["chapter-6.b1.distributed-laplace"]["figure-6.3"]
        require(outputs == ["figure-6.3.png"], "history output name changed")
        require(
            (
                root
                / "postprocess"
                / "comparisons"
                / "figure-6.3.png"
            ).is_file(),
            "history figure was not written",
        )
    print("postprocess history contract passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
