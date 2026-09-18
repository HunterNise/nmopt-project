#!/usr/bin/env python3
"""Focused standard-library contract for Chapter 6 report inputs and output."""

from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY_ROOT / "tools"))

from chapter6_report import (
    build_report,
    read_run_manifest,
    run_from_manifest,
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def manifest_document(artifacts: list[dict[str, str]]) -> dict[str, object]:
    return {
        "schema": "nmopt-run-set-v1",
        "status": "complete",
        "benchmark": "b1",
        "run_kind": "development",
        "build_profile": "debug-dealii",
        "framework_revision": "test-revision",
        "artifacts": artifacts,
    }


def write_manifest(path: Path, document: dict[str, object]) -> None:
    path.write_text(json.dumps(document), encoding="utf-8")


def assert_rejected(path: Path, document: dict[str, object], message: str) -> None:
    write_manifest(path, document)
    try:
        read_run_manifest(path)
    except ValueError:
        return
    raise RuntimeError(message)


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="nmopt-report-") as directory:
        root = Path(directory)
        manifest_path = root / "run-manifest.json"
        artifacts = [
            {
                "path": "artifacts/steepest-descent/beta-1e-1/artifact.kv",
                "status": "ok",
            },
            {
                "path": "artifacts/l-bfgs/beta-1e-1/artifact.kv",
                "status": "error",
                "error": "synthetic failure",
            },
            {
                "path": "artifacts/l-bfgs/beta-1e-2/artifact.kv",
                "status": "pending",
            },
        ]
        document = manifest_document(artifacts)
        write_manifest(manifest_path, document)
        manifest = read_run_manifest(manifest_path)
        require(manifest.status == "complete", "manifest status was not preserved")
        require(manifest.benchmark == "b1", "manifest benchmark was not preserved")
        require(
            manifest.run_kind == "development"
            and manifest.build_profile == "debug-dealii"
            and manifest.framework_revision == "test-revision",
            "manifest metadata was not preserved",
        )
        require(
            [artifact.status for artifact in manifest.artifacts]
            == ["ok", "error", "pending"],
            "artifact statuses were not preserved",
        )
        require(
            manifest.artifacts[1].diagnostic == "synthetic failure",
            "artifact diagnostic was not preserved",
        )

        for unsafe_path in (
            "../artifact.kv",
            "/absolute/path/artifact.kv",
            "artifacts/case/not-artifact.txt",
        ):
            unsafe = manifest_document(
                [{"path": unsafe_path, "status": "ok"}]
            )
            assert_rejected(
                root / "unsafe.json",
                unsafe,
                f"unsafe artifact path was accepted: {unsafe_path}",
            )
        duplicate = manifest_document(
            [
                {"path": "artifacts/a/artifact.kv", "status": "ok"},
                {"path": "artifacts/a/artifact.kv", "status": "pending"},
            ]
        )
        assert_rejected(
            root / "duplicate.json",
            duplicate,
            "duplicate artifact path was accepted",
        )

        unsupported_schema = manifest_document([])
        unsupported_schema["schema"] = "unsupported"
        assert_rejected(
            root / "schema.json",
            unsupported_schema,
            "unsupported manifest schema was accepted",
        )
        invalid_run_status = manifest_document([])
        invalid_run_status["status"] = "invalid"
        assert_rejected(
            root / "run-status.json",
            invalid_run_status,
            "invalid run-set status was accepted",
        )
        invalid_artifact_status = manifest_document(
            [{"path": "artifacts/a/artifact.kv", "status": "invalid"}]
        )
        assert_rejected(
            root / "artifact-status.json",
            invalid_artifact_status,
            "invalid artifact status was accepted",
        )

        artifact_path = root / artifacts[0]["path"]
        artifact_path.parent.mkdir(parents=True)
        artifact_path.write_text(
            "\n".join(
                (
                    "identity.scenario_id=chapter-6.b1.distributed-laplace",
                    r"identity.output_id=actual\=output",
                    "benchmark.method=actual-method",
                    "benchmark.regularisation=actual-beta",
                    "solver.objective_history=1.0,0.5,0.25",
                    "solver.gradient_norm_history=0.4,0.2,0.1",
                    "solver.accepted_iterations=3",
                    "solver.stopping_reason=gradient tolerance",
                )
            )
            + "\n",
            encoding="utf-8",
        )
        run = run_from_manifest(manifest, manifest.artifacts[0])
        require(run.artifact_status == "ok", "valid artifact was not accepted")
        require(
            run.output_id == "actual=output"
            and run.method_or_case == "actual-method"
            and run.regularisation == "actual-beta",
            "artifact metadata did not override fallback metadata",
        )
        require(
            run.objective_history == (1.0, 0.5, 0.25)
            and run.gradient_history == (0.4, 0.2, 0.1),
            "valid histories were not parsed",
        )
        require(
            run.values["solver.accepted_iterations"] == "3"
            and run.values["solver.stopping_reason"] == "gradient tolerance",
            "persisted solver metadata was not retained",
        )

        missing_document = manifest_document(
            [{"path": "artifacts/missing/artifact.kv", "status": "ok"}]
        )
        missing_manifest_path = root / "missing-manifest.json"
        write_manifest(missing_manifest_path, missing_document)
        missing_manifest = read_run_manifest(missing_manifest_path)
        missing = run_from_manifest(missing_manifest, missing_manifest.artifacts[0])
        require(
            missing.artifact_status == "missing"
            and missing.diagnostic
            == "manifest marks artifact as ok, but artifact.kv is absent",
            "missing successful artifact was not diagnosed",
        )
        error_run = run_from_manifest(manifest, manifest.artifacts[1])
        require(
            error_run.artifact_status == "error"
            and error_run.diagnostic == "synthetic failure",
            "error artifact diagnostic was not retained",
        )
        pending_run = run_from_manifest(manifest, manifest.artifacts[2])
        require(
            pending_run.artifact_status == "pending"
            and pending_run.diagnostic
            == "artifact was still pending when the run manifest was written",
            "pending artifact diagnostic was not generated",
        )

        output_root = root / "report"
        build_report(root, output_root, manifest_path)
        summary_csv = output_root / "summary.csv"
        summary_md = output_root / "summary.md"
        require(summary_csv.is_file() and summary_md.is_file(), "report files were not created")
        csv_text = summary_csv.read_text(encoding="utf-8")
        markdown = summary_md.read_text(encoding="utf-8")
        for artifact in artifacts:
            require(
                artifact["path"] in csv_text,
                f"artifact record missing from summary.csv: {artifact['path']}",
            )
        require(
            "synthetic failure" in csv_text
            and "pending when the run manifest was written" in csv_text,
            "status diagnostics missing from summary.csv",
        )
        require(
            "debug-dealii" in markdown
            and "test-revision" in markdown
            and "Run-set status: `complete`" in markdown,
            "manifest metadata missing from summary.md",
        )

    print("chapter6 report contract passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
