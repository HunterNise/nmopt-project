#!/usr/bin/env python3
"""Render Chapter 6 native fields with meshio and matplotlib.

The tool reads one runner artifact directory and writes derived plots outside
the authoritative ``native/`` directory.  It accepts the current native file
names and the historical ``fields.vtu``/``control.vtu`` names used by older
development runs.  With ``--input`` it processes every artifact below a run
root into a mirrored output tree and writes an aggregate index.
"""

from __future__ import annotations

from pathlib import Path

from nmopt_postprocess.chapter6 import CHAPTER6_PROFILE
from nmopt_postprocess.pipeline import (
    process_artifact as process_profile_artifact,
    process_input_root as process_profile_input_root,
)
from postprocess import main as postprocess_main


def process_artifact(artifact: Path, output: Path) -> dict[str, object]:
    return process_profile_artifact(artifact, output, CHAPTER6_PROFILE)


def process_input_root(input_root: Path, output_root: Path) -> dict[str, object]:
    return process_profile_input_root(input_root, output_root, CHAPTER6_PROFILE)


def main() -> int:
    return postprocess_main(default_profile="chapter6", description=__doc__)


if __name__ == "__main__":
    raise SystemExit(main())
