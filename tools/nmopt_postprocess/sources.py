"""Resolve native and legacy files within an artifact."""

from collections.abc import Iterable
from pathlib import Path


def first_existing(artifact: Path, names: Iterable[str]) -> Path | None:
    """Return the first matching native or legacy artifact file."""

    candidates = [artifact / "native" / name for name in names]
    candidates.extend(artifact / name for name in names)
    return next((path for path in candidates if path.is_file()), None)
