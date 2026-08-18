#!/usr/bin/env python3
"""Build deterministic Chapter 6 benchmark summaries and SVG plots.

The report consumes the runner's persisted ``artifact.kv`` files and the
optional ``solver-trace.csv`` sidecars.  It deliberately has no third-party
dependencies so it can be used on the same minimal environment as the
benchmark runner.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Sequence


SUMMARY_FIELDS = (
    "benchmark",
    "output_id",
    "artifact_path",
    "artifact_status",
    "diagnostic",
    "run_directory",
    "method_or_case",
    "regularisation",
    "mesh_refinement",
    "initial_objective",
    "final_objective",
    "initial_gradient",
    "final_gradient",
    "accepted_iterations",
    "line_search_trial_count",
    "stopping_reason",
    "trace_rows",
    "field_outputs",
    "native_outputs",
)


@dataclass(frozen=True)
class ManifestArtifact:
    relative_path: str
    status: str
    diagnostic: str


@dataclass(frozen=True)
class RunManifest:
    path: Path
    status: str
    benchmark: str
    run_kind: str
    build_profile: str
    framework_revision: str
    artifacts: tuple[ManifestArtifact, ...]


@dataclass(frozen=True)
class Run:
    artifact_path: Path
    values: dict[str, str]
    objective_history: tuple[float, ...]
    gradient_history: tuple[float, ...]
    trace: tuple[dict[str, str], ...]
    artifact_status: str = "ok"
    diagnostic: str = ""
    relative_artifact_path: str = ""
    fallback_benchmark: str = ""
    fallback_output_id: str = ""
    fallback_method_or_case: str = "unknown"
    fallback_regularisation: str = "n/a"

    @property
    def directory(self) -> Path:
        return self.artifact_path.parent

    @property
    def scenario(self) -> str:
        return self.values.get("identity.scenario_id") or self.fallback_benchmark or "unknown"

    @property
    def output_id(self) -> str:
        return self.values.get("identity.output_id") or self.fallback_output_id or self.directory.name

    @property
    def method_or_case(self) -> str:
        return (
            self.values.get("benchmark.method")
            or self.values.get("benchmark.graetz_case")
            or self.values.get("solver.method")
            or self.fallback_method_or_case
        )

    @property
    def regularisation(self) -> str:
        return (
            self.values.get("benchmark.regularisation")
            or self.values.get("b1.regularisation_weight")
            or self.fallback_regularisation
        )

    @property
    def field_outputs(self) -> str:
        native_directory = self.directory / "native"
        names = [
            name
            for name in ("fields-volume.vtu", "control-boundary.vtu")
            if (native_directory / name).is_file()
            or (self.directory / name).is_file()
        ]
        return ",".join(names) if names else "none"

    @property
    def native_outputs(self) -> str:
        native_directory = self.directory / "native"
        names = [
            name
            for name in (
                "mesh-volume.vtu",
                "mesh-volume.svg",
                "fields-volume.vtu",
                "control-boundary.vtu",
            )
            if (native_directory / name).is_file()
            or (self.directory / name).is_file()
        ]
        return ",".join(names) if names else "none"


def unescape(value: str) -> str:
    """Decode the escaping used by the C++ deterministic artifact writer."""

    decoded: list[str] = []
    index = 0
    escaped = {"\\": "\\", "n": "\n", "r": "\r", "t": "\t", "=": "="}
    while index < len(value):
        if value[index] == "\\" and index + 1 < len(value):
            decoded.append(escaped.get(value[index + 1], value[index + 1]))
            index += 2
        else:
            decoded.append(value[index])
            index += 1
    return "".join(decoded)


def read_artifact(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key] = unescape(value)
    return values


def read_history(values: dict[str, str], key: str) -> tuple[float, ...]:
    raw = values.get(key, "")
    numbers: list[float] = []
    for item in raw.split(","):
        if not item:
            continue
        try:
            number = float(item)
        except ValueError:
            continue
        if math.isfinite(number):
            numbers.append(number)
    return tuple(numbers)


def read_trace(path: Path) -> tuple[dict[str, str], ...]:
    if not path.is_file():
        return ()
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        return tuple(dict(row) for row in reader if row)


def read_run_manifest(path: Path) -> RunManifest:
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError(f"invalid run manifest '{path}': {error}") from error

    if not isinstance(document, dict):
        raise ValueError(f"run manifest '{path}' must contain a JSON object")
    if document.get("schema") != "nmopt-run-set-v1":
        raise ValueError(
            f"run manifest '{path}' has unsupported schema "
            f"{document.get('schema')!r}"
        )

    status = document.get("status")
    if status not in {"running", "complete", "failed"}:
        raise ValueError(f"run manifest '{path}' has invalid run-set status")

    def required_string(key: str) -> str:
        value = document.get(key)
        if not isinstance(value, str) or not value:
            raise ValueError(f"run manifest '{path}' needs a nonempty '{key}'")
        return value

    raw_artifacts = document.get("artifacts")
    if not isinstance(raw_artifacts, list):
        raise ValueError(f"run manifest '{path}' needs an 'artifacts' array")

    artifacts: list[ManifestArtifact] = []
    seen_paths: set[str] = set()
    for index, raw_artifact in enumerate(raw_artifacts):
        if not isinstance(raw_artifact, dict):
            raise ValueError(
                f"run manifest '{path}' artifact {index} must be an object"
            )
        relative_path = raw_artifact.get("path")
        artifact_status = raw_artifact.get("status")
        diagnostic = raw_artifact.get("error", "")
        if not isinstance(relative_path, str) or not relative_path:
            raise ValueError(
                f"run manifest '{path}' artifact {index} needs a nonempty path"
            )
        relative = Path(relative_path)
        if (
            relative.is_absolute()
            or ".." in relative.parts
            or relative.name != "artifact.kv"
        ):
            raise ValueError(
                f"run manifest '{path}' artifact {index} has an unsafe path"
            )
        if relative_path in seen_paths:
            raise ValueError(
                f"run manifest '{path}' contains duplicate artifact path "
                f"'{relative_path}'"
            )
        if artifact_status not in {"pending", "ok", "error"}:
            raise ValueError(
                f"run manifest '{path}' artifact {index} has invalid status"
            )
        if not isinstance(diagnostic, str):
            raise ValueError(
                f"run manifest '{path}' artifact {index} has a non-string error"
            )
        seen_paths.add(relative_path)
        artifacts.append(ManifestArtifact(relative_path, artifact_status, diagnostic))

    return RunManifest(
        path=path,
        status=status,
        benchmark=required_string("benchmark"),
        run_kind=required_string("run_kind"),
        build_profile=required_string("build_profile"),
        framework_revision=required_string("framework_revision"),
        artifacts=tuple(artifacts),
    )


def fallback_metadata(
    manifest: RunManifest, relative_artifact_path: str
) -> tuple[str, str, str, str]:
    parts = Path(relative_artifact_path).parent.parts
    if parts:
        output_id = "/".join(parts)
    else:
        output_id = relative_artifact_path

    if manifest.benchmark == "b1" and len(parts) >= 2:
        return (
            output_id,
            parts[-2],
            "beta-" + parts[-1].removeprefix("beta-"),
            manifest.benchmark,
        )
    if manifest.benchmark == "b2" and parts:
        return output_id, parts[-1], "n/a", manifest.benchmark
    return output_id, parts[-1] if parts else "unknown", "n/a", manifest.benchmark


def run_from_manifest(
    manifest: RunManifest, artifact: ManifestArtifact
) -> Run:
    artifact_path = manifest.path.parent / artifact.relative_path
    output_id, method_or_case, regularisation, benchmark = fallback_metadata(
        manifest, artifact.relative_path
    )
    status = artifact.status
    diagnostic = artifact.diagnostic
    values: dict[str, str] = {}
    objective_history: tuple[float, ...] = ()
    gradient_history: tuple[float, ...] = ()
    trace: tuple[dict[str, str], ...] = ()

    if status == "ok":
        if not artifact_path.is_file():
            status = "missing"
            diagnostic = "manifest marks artifact as ok, but artifact.kv is absent"
        else:
            try:
                values = read_artifact(artifact_path)
            except OSError as error:
                status = "error"
                diagnostic = f"could not read artifact.kv: {error}"
            else:
                objective_history = read_history(values, "solver.objective_history")
                gradient_history = read_history(
                    values, "solver.gradient_norm_history"
                )
                trace = read_trace(artifact_path.parent / "solver-trace.csv")
    elif status == "pending" and not diagnostic:
        diagnostic = "artifact was still pending when the run manifest was written"
    elif status == "error" and not diagnostic:
        diagnostic = "runner reported an error without a diagnostic"

    return Run(
        artifact_path=artifact_path,
        values=values,
        objective_history=objective_history,
        gradient_history=gradient_history,
        trace=trace,
        artifact_status=status,
        diagnostic=diagnostic,
        relative_artifact_path=artifact.relative_path,
        fallback_benchmark=benchmark,
        fallback_output_id=output_id,
        fallback_method_or_case=method_or_case,
        fallback_regularisation=regularisation,
    )


def load_manifest_runs(manifest: RunManifest) -> list[Run]:
    return [run_from_manifest(manifest, artifact) for artifact in manifest.artifacts]


def load_runs(input_root: Path) -> list[Run]:
    paths = sorted(input_root.rglob("artifact.kv"))
    runs = []
    for path in paths:
        values = read_artifact(path)
        runs.append(
            Run(
                artifact_path=path,
                values=values,
                objective_history=read_history(values, "solver.objective_history"),
                gradient_history=read_history(values, "solver.gradient_norm_history"),
                trace=read_trace(path.parent / "solver-trace.csv"),
            )
        )
    return sorted(runs, key=lambda run: (run.output_id, str(run.artifact_path)))


def value_or_na(values: dict[str, str], key: str) -> str:
    return values.get(key, "n/a")


def summary_row(run: Run, input_root: Path) -> dict[str, str]:
    relative = run.directory.relative_to(input_root)
    artifact_relative = run.artifact_path.relative_to(input_root)
    return {
        "benchmark": run.scenario,
        "output_id": run.output_id,
        "artifact_path": str(artifact_relative),
        "artifact_status": run.artifact_status,
        "diagnostic": run.diagnostic or "n/a",
        "run_directory": str(relative) if str(relative) else ".",
        "method_or_case": run.method_or_case,
        "regularisation": run.regularisation,
        "mesh_refinement": value_or_na(run.values, "benchmark.mesh_refinement"),
        "initial_objective": format_number(run.objective_history[0])
        if run.objective_history
        else "n/a",
        "final_objective": format_number(run.objective_history[-1])
        if run.objective_history
        else value_or_na(run.values, "solver.final_objective"),
        "initial_gradient": format_number(run.gradient_history[0])
        if run.gradient_history
        else "n/a",
        "final_gradient": format_number(run.gradient_history[-1])
        if run.gradient_history
        else "n/a",
        "accepted_iterations": value_or_na(run.values, "solver.accepted_iterations"),
        "line_search_trial_count": value_or_na(
            run.values, "solver.line_search_trial_count"
        ),
        "stopping_reason": value_or_na(run.values, "solver.stopping_reason"),
        "trace_rows": str(len(run.trace)),
        "field_outputs": run.field_outputs,
        "native_outputs": run.native_outputs,
    }


def format_number(number: float) -> str:
    return f"{number:.12g}"


def write_summary_csv(path: Path, rows: Sequence[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def markdown_cell(value: str) -> str:
    return value.replace("|", "\\|").replace("\n", " ")


def write_summary_markdown(
    path: Path,
    input_root: Path,
    runs: Sequence[Run],
    rows: Sequence[dict[str, str]],
    manifest: Optional[RunManifest],
) -> None:
    successful = sum(run.artifact_status == "ok" for run in runs)
    lines = [
        "# Chapter 6 benchmark report",
        "",
        f"- Input root: `{input_root}`",
        f"- Artifact records selected: {len(runs)}",
        f"- Successful artifact records: {successful}",
        "- Data sources: `artifact.kv`, with optional `solver-trace.csv` and VTU sidecars.",
    ]
    if manifest is not None:
        lines.extend(
            [
                f"- Run manifest: `{manifest.path}`",
                f"- Run-set status: `{manifest.status}`",
                f"- Build profile: `{manifest.build_profile}`",
                f"- Framework revision: `{manifest.framework_revision}`",
                f"- Expected artifacts: {len(manifest.artifacts)}",
            ]
        )
    else:
        lines.append(
            "- Run manifest: not supplied; artifacts were recursively discovered "
            "below the input root."
        )
    lines.extend(
        [
            "",
            "The report is a deterministic projection of persisted runner outputs. It does "
            "not rerun a benchmark or infer values that are absent from an artifact.",
            "",
            "## Run summary",
            "",
            "| " + " | ".join(SUMMARY_FIELDS) + " |",
            "| " + " | ".join("---" for _ in SUMMARY_FIELDS) + " |",
        ]
    )
    for row in rows:
        lines.append("| " + " | ".join(markdown_cell(row[key]) for key in SUMMARY_FIELDS) + " |")

    incomplete = [run for run in runs if run.artifact_status != "ok"]
    missing_trace = [
        run.output_id for run in runs if run.artifact_status == "ok" and not run.trace
    ]
    missing_fields = [
        run.output_id
        for run in runs
        if run.artifact_status == "ok" and run.field_outputs == "none"
    ]
    lines.extend(
        [
            "",
            "## Artifact inventory",
            "",
            f"- Objective-history plot: `objective-history.svg`",
            f"- Armijo trial plot: `armijo-trials.svg`",
            f"- Machine-readable table: `summary.csv`",
        ]
    )
    if incomplete:
        lines.extend(
            [
                "",
                "The selected run set contains artifact records that are not successful "
                "`artifact.kv` inputs:",
                "",
                *[
                    f"- `{run.relative_artifact_path}`: `{run.artifact_status}` — "
                    f"{run.diagnostic or 'no diagnostic recorded'}"
                    for run in incomplete
                ],
            ]
        )
    if missing_trace:
        lines.extend(
            [
                "",
                "The following runs have no `solver-trace.csv`; their Armijo trial "
                "diagnostics are unavailable:",
                "",
                *[f"- `{output_id}`" for output_id in missing_trace],
            ]
        )
    if missing_fields:
        lines.extend(
            [
                "",
                "The following runs have no ParaView sidecar fields:",
                "",
                *[f"- `{output_id}`" for output_id in missing_fields],
            ]
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def finite_points(series: Iterable[tuple[float, float]]) -> list[tuple[float, float]]:
    return [
        (x, y)
        for x, y in series
        if math.isfinite(x) and math.isfinite(y)
    ]


def tick_values(low: float, high: float, count: int = 5) -> list[float]:
    if high == low:
        return [low]
    return [low + (high - low) * index / (count - 1) for index in range(count)]


def svg_text(x: float, y: float, text: str, size: int = 13, anchor: str = "start") -> str:
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" font-size="{size}" '
        f'text-anchor="{anchor}" fill="#202124">{html.escape(text)}</text>'
    )


def write_line_plot(
    path: Path,
    title: str,
    y_label: str,
    series: Sequence[tuple[str, Sequence[tuple[float, float]]]],
    log_y: bool = False,
) -> None:
    width, height = 1100, 620
    left, right, top, bottom = 90, 260, 70, 75
    plot_width = width - left - right
    plot_height = height - top - bottom
    prepared = [(label, finite_points(points)) for label, points in series]
    all_points = [point for _, points in prepared for point in points]
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        svg_text(width / 2, 34, title, 20, "middle"),
    ]
    if not all_points:
        svg.append(svg_text(width / 2, height / 2, "No data found", 18, "middle"))
        svg.append("</svg>")
        path.write_text("\n".join(svg) + "\n", encoding="utf-8")
        return

    x_low = min(x for x, _ in all_points)
    x_high = max(x for x, _ in all_points)
    if x_low == x_high:
        x_low -= 1
        x_high += 1
    positive = [(x, y) for x, y in all_points if y > 0] if log_y else all_points
    if not positive:
        log_y = False
        positive = all_points
    y_values = [math.log10(y) for _, y in positive] if log_y else [y for _, y in positive]
    y_low, y_high = min(y_values), max(y_values)
    if y_low == y_high:
        padding = max(abs(y_low) * 0.1, 1.0)
        y_low -= padding
        y_high += padding
    else:
        padding = (y_high - y_low) * 0.08
        y_low -= padding
        y_high += padding

    def x_position(x: float) -> float:
        return left + (x - x_low) / (x_high - x_low) * plot_width

    def y_position(y: float) -> float:
        value = math.log10(y) if log_y else y
        return top + (y_high - value) / (y_high - y_low) * plot_height

    for y_tick in tick_values(y_low, y_high):
        y = top + (y_high - y_tick) / (y_high - y_low) * plot_height
        label = f"10^{y_tick:.1f}" if log_y else format_number(y_tick)
        svg.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left + plot_width}" y2="{y:.2f}" stroke="#e5e7eb"/>')
        svg.append(svg_text(left - 10, y + 4, label, 12, "end"))
    for x_tick in tick_values(x_low, x_high):
        x = x_position(x_tick)
        svg.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top + plot_height}" stroke="#f1f3f4"/>')
        svg.append(svg_text(x, top + plot_height + 24, format_number(x_tick), 12, "middle"))

    svg.extend(
        [
            f'<line x1="{left}" y1="{top + plot_height}" x2="{left + plot_width}" y2="{top + plot_height}" stroke="#202124"/>',
            f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_height}" stroke="#202124"/>',
            svg_text(left + plot_width / 2, height - 20, "sample / iteration", 13, "middle"),
            svg_text(22, top + plot_height / 2, y_label, 13, "middle"),
        ]
    )
    colors = ("#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#17becf")
    legend_x, legend_y = left + plot_width + 25, top + 20
    for index, (label, points) in enumerate(prepared):
        color = colors[index % len(colors)]
        valid = [(x, y) for x, y in points if not log_y or y > 0]
        if valid:
            coordinates = " ".join(f"{x_position(x):.2f},{y_position(y):.2f}" for x, y in valid)
            svg.append(f'<polyline fill="none" stroke="{color}" stroke-width="2" points="{coordinates}"/>')
        legend_line_y = legend_y + index * 24
        svg.append(f'<line x1="{legend_x}" y1="{legend_line_y - 4:.2f}" x2="{legend_x + 18}" y2="{legend_line_y - 4:.2f}" stroke="{color}" stroke-width="3"/>')
        svg.append(svg_text(legend_x + 26, legend_line_y, label, 12))
    svg.append("</svg>")
    path.write_text("\n".join(svg) + "\n", encoding="utf-8")


def objective_series(runs: Sequence[Run]) -> list[tuple[str, list[tuple[float, float]]]]:
    return [
        (
            run.output_id,
            [(float(index), value) for index, value in enumerate(run.objective_history)],
        )
        for run in runs
        if run.objective_history
    ]


def trace_series(runs: Sequence[Run]) -> list[tuple[str, list[tuple[float, float]]]]:
    result: list[tuple[str, list[tuple[float, float]]]] = []
    for run in runs:
        points: list[tuple[float, float]] = []
        for index, row in enumerate(run.trace):
            try:
                objective = float(row["objective_value"])
            except (KeyError, TypeError, ValueError):
                continue
            points.append((float(index), objective))
        if points:
            result.append((run.output_id, points))
    return result


def build_report(
    input_root: Path, output_root: Path, run_manifest_path: Optional[Path] = None
) -> None:
    manifest = None
    if run_manifest_path is None:
        candidate = input_root / "run-manifest.json"
        if candidate.is_file():
            run_manifest_path = candidate
    if run_manifest_path is not None:
        manifest = read_run_manifest(run_manifest_path)
        input_root = manifest.path.parent
        runs = load_manifest_runs(manifest)
    else:
        runs = load_runs(input_root)
    if not runs:
        raise ValueError(f"no artifact.kv files found below {input_root}")
    output_root.mkdir(parents=True, exist_ok=True)
    rows = [summary_row(run, input_root) for run in runs]
    write_summary_csv(output_root / "summary.csv", rows)
    write_summary_markdown(
        output_root / "summary.md", input_root, runs, rows, manifest
    )
    write_line_plot(
        output_root / "objective-history.svg",
        "Chapter 6 objective histories",
        "objective",
        objective_series(runs),
    )
    write_line_plot(
        output_root / "armijo-trials.svg",
        "Chapter 6 Armijo trial objectives",
        "trial objective",
        trace_series(runs),
        log_y=True,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument(
        "--input",
        type=Path,
        help="directory containing runner artifact directories; recursive discovery "
        "is retained for legacy or aggregate inputs",
    )
    selection.add_argument(
        "--run-manifest",
        type=Path,
        help="authoritative run-manifest.json selecting one runner run set",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="directory in which the report files will be written",
    )
    arguments = parser.parse_args()
    try:
        input_root = (
            arguments.input.resolve()
            if arguments.input is not None
            else arguments.run_manifest.resolve().parent
        )
        run_manifest_path = (
            arguments.run_manifest.resolve()
            if arguments.run_manifest is not None
            else None
        )
        build_report(input_root, arguments.output.resolve(), run_manifest_path)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print(f"wrote Chapter 6 report to {arguments.output.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
