#!/usr/bin/env python3
"""Strip complete C++ comment lines and verify the remaining token stream."""

from __future__ import annotations

import argparse
import hashlib
import sys
import tempfile
from pathlib import Path


def read_text(path: Path) -> str:
    with path.open("r", encoding="utf-8", newline="") as source:
        return source.read()


def write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as destination:
        destination.write(content)


def block_end(lines: list[str], start: int) -> tuple[int, int] | None:
    for line_number in range(start, len(lines)):
        offset = lines[line_number].find("*/")
        if offset >= 0:
            return line_number, offset + 2
    return None


def leading_block(lines: list[str]) -> tuple[int, str | None]:
    for line_number, line in enumerate(lines):
        if line.strip():
            if not line.lstrip().startswith("/*"):
                return line_number, None
            end = block_end(lines, line_number)
            if end is None:
                raise ValueError("unterminated leading comment block")
            return end[0] + 1, "".join(lines[line_number : end[0] + 1])
    return len(lines), None


def strip_comment_lines(source: str) -> tuple[str, int]:
    lines = source.splitlines(keepends=True)
    first_code_line, header = leading_block(lines)
    output = lines[:first_code_line]
    removed = 0
    index = first_code_line

    while index < len(lines):
        line = lines[index]
        stripped = line.lstrip()

        if stripped.startswith("//"):
            removed += 1
            index += 1
            continue

        if stripped.startswith("/*"):
            end = block_end(lines, index)
            if end is None:
                raise ValueError(f"unterminated block comment at line {index + 1}")
            end_line, end_offset = end
            if lines[end_line][end_offset:].strip():
                output.extend(lines[index : end_line + 1])
            else:
                removed += end_line - index + 1
            index = end_line + 1
            continue

        output.append(line)
        index += 1

    result = "".join(output)
    if header is not None and not result.startswith(header):
        raise ValueError("the leading comment block was not preserved")
    return result, removed


def quoted_end(source: str, start: int) -> int:
    quote = source[start]
    index = start + 1
    while index < len(source):
        if source[index] == "\\":
            index += 2
        elif source[index] == quote:
            return index + 1
        else:
            index += 1
    raise ValueError(f"unterminated literal at offset {start}")


def raw_string_end(source: str, start: int) -> int | None:
    prefixes = ("u8R\"", "uR\"", "UR\"", "LR\"", "R\"")
    prefix = next((value for value in prefixes if source.startswith(value, start)), None)
    if prefix is None:
        return None
    delimiter_start = start + len(prefix)
    opening = source.find("(", delimiter_start)
    if opening < 0:
        raise ValueError(f"invalid raw string literal at offset {start}")
    delimiter = source[delimiter_start:opening]
    closing_marker = ")" + delimiter + '"'
    closing = source.find(closing_marker, opening + 1)
    if closing < 0:
        raise ValueError(f"unterminated raw string literal at offset {start}")
    return closing + len(closing_marker)


def non_comment_tokens(source: str) -> tuple[str, ...]:
    """Compare code tokens without treating comment markers in literals as comments."""

    tokens: list[str] = []
    index = 0
    while index < len(source):
        if source[index].isspace():
            index += 1
        elif source.startswith("//", index):
            newline = source.find("\n", index + 2)
            index = len(source) if newline < 0 else newline + 1
        elif source.startswith("/*", index):
            end = source.find("*/", index + 2)
            if end < 0:
                raise ValueError(f"unterminated block comment at offset {index}")
            index = end + 2
        else:
            raw_end = raw_string_end(source, index)
            if raw_end is not None:
                tokens.append(source[index:raw_end])
                index = raw_end
            elif source[index] in {'"', "'"}:
                end = quoted_end(source, index)
                tokens.append(source[index:end])
                index = end
            elif source[index].isalnum() or source[index] == "_":
                end = index + 1
                while end < len(source) and (
                    source[end].isalnum() or source[end] in {"_", ".", "'"}
                ):
                    end += 1
                tokens.append(source[index:end])
                index = end
            else:
                tokens.append(source[index])
                index += 1
    return tuple(tokens)


def verify(source: str, candidate: str, expected: Path | None = None) -> None:
    if non_comment_tokens(source) != non_comment_tokens(candidate):
        raise ValueError("non-comment token sequence differs")
    if expected is not None and candidate != read_text(expected):
        raise ValueError(f"regenerated output differs from {expected}")


def sha256(content: str) -> str:
    return hashlib.sha256(content.encode("utf-8")).hexdigest()


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    destination = parser.add_mutually_exclusive_group(required=True)
    destination.add_argument("--output", type=Path)
    destination.add_argument("--check", type=Path, metavar="EXPECTED")
    return parser.parse_args()


def main() -> int:
    options = arguments()
    source = read_text(options.input)
    try:
        candidate, removed = strip_comment_lines(source)
        verify(source, candidate, options.check)
    except ValueError as error:
        print(f"strip_comments.py: error: {error}", file=sys.stderr)
        return 1

    if options.output is not None:
        write_text(options.output, candidate)
        print(f"generated {options.output} ({removed} comment lines removed, "
              f"sha256={sha256(candidate)})")
        return 0

    with tempfile.TemporaryDirectory(prefix="nmopt-external-dealii-strip-") as directory:
        temporary = Path(directory) / "candidate.cc"
        write_text(temporary, candidate)
        if read_text(temporary) != read_text(options.check):
            print("strip_comments.py: error: temporary output changed", file=sys.stderr)
            return 1
    print(f"verified {options.check} (sha256={sha256(candidate)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
