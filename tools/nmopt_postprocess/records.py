"""Compatibility adapters for artifact records used by post-processing."""

from pathlib import Path

from nmopt_artifacts.records import (
    read_artifact_file,
    read_numeric_history,
    unescape,
)


def read_metadata(artifact: Path) -> dict[str, str]:
    """Read ``artifact.kv`` from an artifact, if it is present."""

    path = artifact / "artifact.kv"
    if not path.is_file():
        return {}
    return read_artifact_file(path)
