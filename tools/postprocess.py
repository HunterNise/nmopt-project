#!/usr/bin/env python3
"""Process native runner fields with a selected post-processing profile."""

from __future__ import annotations

import argparse
from pathlib import Path

from nmopt_postprocess import PostprocessError
from nmopt_postprocess.chapter6 import CHAPTER6_PROFILE
from nmopt_postprocess.pipeline import (
    PostprocessProfile,
    process_artifact,
    process_input_root,
)
from nmopt_postprocess.render import (
    DEFAULT_OUTPUT_FORMATS,
    OutputFormat,
    OutputFormats,
)


PROFILES: dict[str, PostprocessProfile] = {
    "chapter6": CHAPTER6_PROFILE,
}
OUTPUT_FORMATS: tuple[OutputFormat, ...] = ("png", "svg")


def build_parser(
    default_profile: str | None = None, description: str | None = None
) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=description or __doc__)
    if default_profile is None:
        parser.add_argument(
            "--profile",
            choices=tuple(PROFILES),
            required=True,
            help="application profile used to select fields and comparisons",
        )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--artifact",
        type=Path,
        help="one runner artifact directory containing artifact.kv/native",
    )
    source.add_argument(
        "--input",
        type=Path,
        help="run root containing one or more artifact.kv files",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help=(
            "derived output directory; defaults to <artifact>/postprocess for "
            "--artifact or <input>/postprocess for --input"
        ),
    )
    parser.add_argument(
        "--format",
        dest="output_formats",
        nargs="+",
        choices=OUTPUT_FORMATS,
        default=DEFAULT_OUTPUT_FORMATS,
        metavar="FORMAT",
        help="one or more plot formats; defaults to png",
    )
    return parser


def main(
    default_profile: str | None = None, description: str | None = None
) -> int:
    parser = build_parser(default_profile, description)
    arguments = parser.parse_args()
    profile_name = default_profile or arguments.profile
    profile = PROFILES[profile_name]
    output_formats: OutputFormats = tuple(arguments.output_formats)
    try:
        if arguments.artifact is not None:
            artifact = arguments.artifact.resolve()
            output = (
                arguments.output.resolve()
                if arguments.output is not None
                else artifact / "postprocess"
            )
            manifest = process_artifact(
                artifact,
                output,
                profile,
                output_formats=output_formats,
            )
            print(f"wrote {len(manifest['plots'])} field plot sets to {output}")
            return 0

        input_root = arguments.input.resolve()
        output_root = (
            arguments.output.resolve()
            if arguments.output is not None
            else input_root / "postprocess"
        )
        index = process_input_root(
            input_root,
            output_root,
            profile,
            output_formats=output_formats,
        )
        print(
            f"processed {index['success_count']}/{index['artifact_count']} "
            f"artifacts to {output_root}"
        )
        return (
            0
            if index["failure_count"] == 0
            and not index.get("comparison_errors")
            else 1
        )
    except (OSError, ValueError, PostprocessError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    raise SystemExit(main())
