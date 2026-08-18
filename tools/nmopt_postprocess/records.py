"""Read the small key-value records emitted for runner artifacts."""

from pathlib import Path


def unescape(value: str) -> str:
    """Decode the escapes used by the runner's ``artifact.kv`` format."""

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


def read_metadata(artifact: Path) -> dict[str, str]:
    """Read ``artifact.kv`` from an artifact, if it is present."""

    path = artifact / "artifact.kv"
    if not path.is_file():
        return {}
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line and "=" in line:
            key, value = line.split("=", 1)
            values[key] = unescape(value)
    return values
