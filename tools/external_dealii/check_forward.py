#!/usr/bin/env python3
"""Compare two standalone legacy-ASCII VTK forward runs."""

from __future__ import annotations

import argparse
import subprocess
import sys
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


@dataclass(frozen=True)
class VtkData:
    points: tuple[tuple[float, ...], ...]
    cells: tuple[tuple[int, ...], ...]
    cell_types: tuple[int, ...]
    arrays: dict[tuple[str, str], tuple[float, ...]]


class ComparisonError(Exception):
    pass


def tokens(path: Path) -> list[str]:
    return path.read_text(encoding="utf-8").split()


def take(values: list[str], index: int) -> tuple[str, int]:
    if index >= len(values):
        raise ComparisonError("unexpected end of VTK file")
    return values[index], index + 1


def parse_vtk(path: Path) -> VtkData:
    values = tokens(path)
    try:
        dataset = values.index("DATASET")
    except ValueError as error:
        raise ComparisonError(f"{path} has no DATASET section") from error

    dataset_type, index = take(values, dataset + 1)
    if dataset_type != "UNSTRUCTURED_GRID":
        raise ComparisonError(f"{path} uses unsupported dataset {dataset_type!r}")

    points: tuple[tuple[float, ...], ...] = ()
    cells: tuple[tuple[int, ...], ...] = ()
    cell_types: tuple[int, ...] = ()
    arrays: dict[tuple[str, str], tuple[float, ...]] = {}
    data_location = ""
    data_count = 0

    while index < len(values):
        section, index = take(values, index)
        if section == "POINTS":
            count_text, index = take(values, index)
            _, index = take(values, index)
            count = int(count_text)
            raw, index = values[index : index + 3 * count], index + 3 * count
            if len(raw) != 3 * count:
                raise ComparisonError(f"{path} has incomplete POINTS data")
            points = tuple(
                tuple(float(raw[offset + component]) for component in range(3))
                for offset in range(0, len(raw), 3)
            )
        elif section == "CELLS":
            count_text, index = take(values, index)
            _, index = take(values, index)
            cells_list: list[tuple[int, ...]] = []
            for _ in range(int(count_text)):
                size_text, index = take(values, index)
                size = int(size_text)
                cell, index = values[index : index + size], index + size
                if len(cell) != size:
                    raise ComparisonError(f"{path} has incomplete CELLS data")
                cells_list.append(tuple(int(value) for value in cell))
            cells = tuple(cells_list)
        elif section == "CELL_TYPES":
            count_text, index = take(values, index)
            count = int(count_text)
            raw, index = values[index : index + count], index + count
            if len(raw) != count:
                raise ComparisonError(f"{path} has incomplete CELL_TYPES data")
            cell_types = tuple(int(value) for value in raw)
        elif section in {"POINT_DATA", "CELL_DATA"}:
            data_location = section
            count_text, index = take(values, index)
            data_count = int(count_text)
        elif section == "SCALARS":
            name, index = take(values, index)
            _, index = take(values, index)
            components_text, index = take(values, index)
            components = int(components_text)
            lookup, index = take(values, index)
            if lookup != "LOOKUP_TABLE":
                raise ComparisonError(f"{path} has malformed SCALARS metadata")
            _, index = take(values, index)
            raw, index = values[index : index + data_count * components], index + data_count * components
            if len(raw) != data_count * components:
                raise ComparisonError(f"{path} has incomplete SCALARS data")
            arrays[(data_location, name)] = tuple(float(value) for value in raw)
        elif section == "VECTORS":
            name, index = take(values, index)
            _, index = take(values, index)
            raw, index = values[index : index + 3 * data_count], index + 3 * data_count
            if len(raw) != 3 * data_count:
                raise ComparisonError(f"{path} has incomplete VECTORS data")
            arrays[(data_location, name)] = tuple(float(value) for value in raw)
        else:
            raise ComparisonError(f"{path} has unsupported VTK section {section!r}")

    return VtkData(points, cells, cell_types, arrays)


def max_numeric_difference(left: tuple[tuple[float, ...], ...], right: tuple[tuple[float, ...], ...]) -> float:
    if len(left) != len(right):
        raise ComparisonError(f"numeric sequence lengths differ: {len(left)} != {len(right)}")
    difference = 0.0
    for left_row, right_row in zip(left, right):
        if len(left_row) != len(right_row):
            raise ComparisonError("numeric row lengths differ")
        difference = max(
            difference,
            *(abs(left_value - right_value) for left_value, right_value in zip(left_row, right_row)),
        )
    return difference


def compare_vtk(upstream: Path, stripped: Path, abs_tol: float, rel_tol: float) -> str:
    left = parse_vtk(upstream)
    right = parse_vtk(stripped)
    if left.cells != right.cells:
        raise ComparisonError("cell connectivity differs")
    if left.cell_types != right.cell_types:
        raise ComparisonError("cell types differ")
    if left.arrays.keys() != right.arrays.keys():
        raise ComparisonError("VTK data-array names or locations differ")

    point_difference = max_numeric_difference(left.points, right.points)
    value_difference = 0.0
    for key in left.arrays:
        values_left = (left.arrays[key],)
        values_right = (right.arrays[key],)
        value_difference = max(value_difference, max_numeric_difference(values_left, values_right))

    scale = max(
        1.0,
        *(abs(value) for point in left.points for value in point),
        *(abs(value) for array in left.arrays.values() for value in array),
    )
    tolerance = abs_tol + rel_tol * scale
    if point_difference > tolerance or value_difference > tolerance:
        raise ComparisonError(
            f"numeric VTK data differs: points={point_difference:g}, "
            f"arrays={value_difference:g}, tolerance={tolerance:g}"
        )
    return (
        f"points={len(left.points)}, cells={len(left.cells)}, "
        f"arrays={len(left.arrays)}, max_point_difference={point_difference:g}, "
        f"max_array_difference={value_difference:g}"
    )


def run(executable: Path, directory: Path) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(
        [str(executable)],
        cwd=directory,
        capture_output=True,
        text=True,
        check=False,
    )
    (directory / "stdout.txt").write_text(result.stdout, encoding="utf-8")
    (directory / "stderr.txt").write_text(result.stderr, encoding="utf-8")
    if result.returncode != 0:
        raise ComparisonError(f"{executable} exited with {result.returncode}")
    return result


def make_run_directory(root: Path) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    name = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid.uuid4().hex[:8]
    directory = root / name
    directory.mkdir()
    (directory / "upstream").mkdir()
    (directory / "stripped").mkdir()
    return directory


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upstream-executable", required=True, type=Path)
    parser.add_argument("--stripped-executable", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--file", action="append", required=True, dest="files")
    parser.add_argument("--absolute-tolerance", type=float, default=1e-14)
    parser.add_argument("--relative-tolerance", type=float, default=1e-12)
    return parser.parse_args()


def main() -> int:
    options = parse_arguments()
    upstream_executable = options.upstream_executable.resolve()
    stripped_executable = options.stripped_executable.resolve()
    run_directory = make_run_directory(options.output_root.resolve())
    report: list[str] = [
        f"upstream executable: {upstream_executable}",
        f"stripped executable: {stripped_executable}",
        f"run directory: {run_directory}",
    ]

    try:
        upstream = run(upstream_executable, run_directory / "upstream")
        stripped = run(stripped_executable, run_directory / "stripped")
        if upstream.stdout != stripped.stdout:
            raise ComparisonError("stdout differs")
        report.append("stdout: identical")
        for filename in options.files:
            upstream_file = run_directory / "upstream" / filename
            stripped_file = run_directory / "stripped" / filename
            if not upstream_file.is_file() or not stripped_file.is_file():
                raise ComparisonError(f"missing expected output {filename}")
            summary = compare_vtk(
                upstream_file,
                stripped_file,
                options.absolute_tolerance,
                options.relative_tolerance,
            )
            report.append(f"{filename}: identical ({summary})")
        report.insert(0, "status: PASS")
        exit_code = 0
    except (ComparisonError, OSError) as error:
        report.insert(0, f"status: FAIL ({error})")
        exit_code = 1

    (run_directory / "comparison.txt").write_text("\n".join(report) + "\n", encoding="utf-8")
    print("\n".join(report))
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
