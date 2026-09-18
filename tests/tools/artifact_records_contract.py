#!/usr/bin/env python3
"""Focused contract for dependency-free persisted artifact parsing."""

from __future__ import annotations

import tempfile
from pathlib import Path

import sys


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY_ROOT / "tools"))

from nmopt_artifacts.records import (
    read_artifact_file,
    read_numeric_history,
    unescape,
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> int:
    require(
        unescape(r"slash\\ newline\n return\r tab\t equals\=")
        == "slash\\ newline\n return\r tab\t equals=",
        "artifact escaping was not decoded",
    )
    require(
        unescape(r"left\=middle\\right\tend") == "left=middle\\right\tend",
        "multiple escaped delimiters were not decoded",
    )

    with tempfile.TemporaryDirectory(prefix="nmopt-artifact-records-") as directory:
        root = Path(directory)
        artifact_path = root / "artifact.kv"
        artifact_path.write_text(
            "\n".join(
                (
                    "",
                    "this line has no separator",
                    r"escaped=hello\nworld\=value",
                    "equals=left=right",
                    "duplicate=first",
                    "duplicate=last",
                )
            )
            + "\n",
            encoding="utf-8",
        )
        values = read_artifact_file(artifact_path)
        require(
            "this line has no separator" not in values,
            "invalid artifact line was not ignored",
        )
        require(
            values["escaped"] == "hello\nworld=value",
            "escaped artifact value was not decoded",
        )
        require(
            values["equals"] == "left=right",
            "artifact parser did not split at the first equals sign",
        )
        require(values["duplicate"] == "last", "duplicate key was not last-write-wins")

        missing_path = root / "missing.kv"
        try:
            read_artifact_file(missing_path)
        except FileNotFoundError:
            pass
        else:
            raise RuntimeError("missing artifact file did not raise FileNotFoundError")

    require(read_numeric_history({}, "history") == (), "missing history was not empty")
    require(
        read_numeric_history({"history": ""}, "history") == (),
        "empty history was not empty",
    )
    require(
        read_numeric_history({"history": "1,2.5,-3e-2"}, "history")
        == (1.0, 2.5, -0.03),
        "valid history was parsed incorrectly",
    )

    invalid_histories = {
        "1,,2": "non-numeric value at index 1",
        "1,abc,2": "non-numeric value at index 1",
        "1,nan,2": "non-finite value at index 1",
        "1,inf,2": "non-finite value at index 1",
        "1,-inf,2": "non-finite value at index 1",
    }
    for text, expected in invalid_histories.items():
        try:
            read_numeric_history({"solver.history": text}, "solver.history")
        except ValueError as error:
            message = str(error)
            require(
                "solver.history" in message and expected in message,
                f"history diagnostic was incomplete for {text!r}: {message}",
            )
        else:
            raise RuntimeError(f"invalid history was accepted: {text!r}")

    print("artifact records contract passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
