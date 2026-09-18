"""Read persisted runner artifact records without third-party dependencies."""

from collections.abc import Mapping
import math
from pathlib import Path


__all__ = ["unescape", "read_artifact_file", "read_numeric_history"]


def unescape(value: str) -> str:
    """Decode the escaping used by the runner's ``artifact.kv`` format."""

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


def read_artifact_file(path: Path) -> dict[str, str]:
    """Read one escaped ``artifact.kv`` file."""

    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line and "=" in line:
            key, value = line.split("=", 1)
            values[key] = unescape(value)
    return values


def read_numeric_history(
    values: Mapping[str, str], key: str
) -> tuple[float, ...]:
    """Read one finite comma-separated solver history."""

    text = values.get(key, "")
    if not text:
        return ()
    result: list[float] = []
    for index, item in enumerate(text.split(",")):
        try:
            value = float(item)
        except ValueError as error:
            raise ValueError(
                f"history '{key}' has a non-numeric value at index {index}"
            ) from error
        if not math.isfinite(value):
            raise ValueError(
                f"history '{key}' has a non-finite value at index {index}"
            )
        result.append(value)
    return tuple(result)
