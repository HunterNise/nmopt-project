#!/usr/bin/env python3
"""Focused contracts for the external deal.II forward comparator."""

from __future__ import annotations

import stat
import subprocess
import sys
import tempfile
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY_ROOT / "tools" / "external_dealii"))

from check_forward import ComparisonError, compare_vtk  # noqa: E402


VTK_TEMPLATE = """# vtk DataFile Version 3.0
forward comparator contract
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 1 float
{point}
CELLS 1 2
1 0
CELL_TYPES 1
1
POINT_DATA 1
SCALARS field float 1
LOOKUP_TABLE default
{field}
"""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def write_vtk(path: Path, point: str = "0 0 0", field: str = "1") -> None:
    path.write_text(
        VTK_TEMPLATE.format(point=point, field=field), encoding="utf-8"
    )


def expect_comparison_failure(
    upstream: Path, stripped: Path, expected: str
) -> None:
    try:
        compare_vtk(upstream, stripped, 1e-14, 1e-12)
    except ComparisonError as error:
        require(expected in str(error), f"unexpected comparison error: {error}")
    else:
        raise RuntimeError("non-equivalent VTK data was accepted")


def write_runner(path: Path, point: str, field: str) -> None:
    path.write_text(
        "#!/usr/bin/env python3\n"
        "from pathlib import Path\n"
        f"Path('solution.vtk').write_text({VTK_TEMPLATE.format(point=point, field=field)!r}, encoding='utf-8')\n"
        "print('fixture stdout')\n",
        encoding="utf-8",
    )
    path.chmod(path.stat().st_mode | stat.S_IXUSR)


def check_cli_failure() -> None:
    with tempfile.TemporaryDirectory(prefix="nmopt-forward-contract-") as directory:
        root = Path(directory)
        upstream_runner = root / "upstream-fixture"
        stripped_runner = root / "stripped-fixture"
        write_runner(upstream_runner, "0 0 0", "1")
        write_runner(stripped_runner, "0 0 0", "nan")
        output_root = root / "runs"
        result = subprocess.run(
            [
                sys.executable,
                str(REPOSITORY_ROOT / "tools/external_dealii/check_forward.py"),
                "--upstream-executable",
                str(upstream_runner),
                "--stripped-executable",
                str(stripped_runner),
                "--output-root",
                str(output_root),
                "--file",
                "solution.vtk",
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        require(result.returncode != 0, "CLI accepted nonfinite VTK data")
        run_directories = list(output_root.iterdir())
        require(len(run_directories) == 1, "CLI did not create one run directory")
        report = run_directories[0] / "comparison.txt"
        require(report.is_file(), "CLI failure did not retain comparison report")
        report_text = report.read_text(encoding="utf-8")
        require("status: FAIL" in report_text, "failure report lost its status")
        require("nonfinite" in report_text, "failure report lost the root cause")


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="nmopt-forward-data-") as directory:
        root = Path(directory)
        finite_upstream = root / "finite-upstream.vtk"
        finite_stripped = root / "finite-stripped.vtk"
        write_vtk(finite_upstream)
        write_vtk(finite_stripped)
        require(
            "points=1, cells=1, arrays=1"
            in compare_vtk(finite_upstream, finite_stripped, 1e-14, 1e-12),
            "identical finite VTK data no longer compares successfully",
        )

        write_vtk(finite_stripped, field="1.1")
        expect_comparison_failure(finite_upstream, finite_stripped, "numeric VTK data differs")

        for side in ("upstream", "stripped"):
            for location in ("point", "field"):
                for value in ("nan", "inf", "-inf"):
                    upstream_path = root / "nonfinite-upstream.vtk"
                    stripped_path = root / "nonfinite-stripped.vtk"
                    write_vtk(upstream_path)
                    write_vtk(stripped_path)
                    target = upstream_path if side == "upstream" else stripped_path
                    if location == "point":
                        write_vtk(target, point=f"{value} 0 0")
                        expected = f"{side} point"
                    else:
                        write_vtk(target, field=value)
                        expected = f"{side} array"
                    expect_comparison_failure(upstream_path, stripped_path, expected)

    check_cli_failure()
    print("external deal.II forward comparator contract passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
