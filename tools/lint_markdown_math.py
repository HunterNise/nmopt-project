#!/usr/bin/env python3
"""Lint repository Markdown/LaTeX against nmopt's portable-math conventions.

This is a read-only project linter. It detects known Markdown/TeX interactions
and deterministic house-style violations, but it deliberately does not rewrite
mathematics: many correct fixes require understanding the mathematical role of
the notation.

Usage from the repository root
------------------------------

Scan every tracked Markdown file::

    python3 tools/lint_markdown_math.py

Scan only Markdown that is currently modified or untracked::

    python3 tools/lint_markdown_math.py --changed

Scan one or more files/directories::

    python3 tools/lint_markdown_math.py docs/manual README.md

Emit machine-readable findings for agent/tool consumption::

    python3 tools/lint_markdown_math.py --changed --json

List the implemented rules::

    python3 tools/lint_markdown_math.py --list-rules

The default exit status is nonzero when a warning or error is present. Use
``--fail-on error`` to make only errors fail, ``--fail-on info`` to include
advisory findings, or ``--fail-on none`` for a report-only run.

By default, style-only findings under ``docs/history/`` are suppressed; actual
rendering/parser hazards are still reported there. Pass
``--include-history-style`` when historical material is deliberately being
modernized.

Scope and limits
----------------

The linter is a fast deterministic guard for the project's documented portable
subset. It is not a complete CommonMark/GFM parser, a TeX parser, or a substitute
for rendering unusual expressions in GitHub, VS Code, and another target reader.
In particular, novel macros/environments, deeply unusual Markdown nesting, and
semantic typography choices still require review. A finding identifies a place
to inspect; it does not authorize a mechanical rewrite.
"""

from __future__ import annotations

import argparse
import json
import os
from dataclasses import dataclass
from pathlib import Path
import re
import subprocess
import sys
import unicodedata
from typing import Iterable, Iterator, Sequence


# ---------------------------------------------------------------------------
# Findings and rules
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Finding:
    path: str
    line: int
    column: int
    severity: str
    rule: str
    message: str
    excerpt: str
    style_only: bool = False

    def as_dict(self) -> dict[str, object]:
        return {
            "path": self.path,
            "line": self.line,
            "column": self.column,
            "severity": self.severity,
            "rule": self.rule,
            "message": self.message,
            "excerpt": self.excerpt,
        }


RULES: dict[str, str] = {
    "MATH001": "raw vertical bar in mathematics; use a named delimiter/operator",
    "MATH002": "Markdown-sensitive TeX control-symbol escape",
    "MATH003": "unsupported or discouraged TeX macro",
    "MATH004": "possible GFM emphasis pairing involving math underscore(s)",
    "MATH005": "Markdown block construct inside dollar-delimited display math",
    "MATH006": "noncanonical display-math delimiter/source layout",
    "MATH007": "literal asterisk in mathematics",
    "MATH008": "literal less-than/greater-than in mathematics",
    "MATH009": "literal := definition operator",
    "MATH010": "unbraced superscript or subscript",
    "MATH011": "generic vertical delimiter used as restriction/evaluation bar",
    "MATH012": "custom operator with side script should be reviewed for \\nolimits",
    "MATH013": "integral differential spacing should be reviewed",
    "MATH014": "repository Markdown uses noncanonical math delimiters",
    "MATH015": "display mathematics inside a Markdown table row",
    "MATH016": "top-level align environment is renderer-sensitive",
    "MATH017": "unmatched ordinary dollar math delimiter",
    "MD001": "two-space hard line break; prefer a terminal backslash",
}

SEVERITY_ORDER = {"error": 3, "warning": 2, "info": 1}

# Macros which are either known renderer problems or explicitly outside the
# project's portable subset.  The value is (severity, explanation, style_only).
FORBIDDEN_MACROS: dict[str, tuple[str, str, bool]] = {
    r"\operatorname": (
        "error",
        r"use a built-in operator or \mathop{\mathrm{...}}; GitHub rejects \operatorname",
        False,
    ),
    r"\operatorname*": (
        "error",
        r"use a built-in operator or \mathop{\mathrm{...}}; GitHub rejects \operatorname*",
        False,
    ),
    r"\tag": (
        "warning",
        r"avoid \tag; numbering/layout is renderer-sensitive",
        False,
    ),
    r"\label": (
        "warning",
        r"use prose/Markdown links instead of TeX \label/\ref navigation",
        True,
    ),
    r"\ref": (
        "warning",
        r"use prose/Markdown links instead of TeX \label/\ref navigation",
        True,
    ),
    r"\eqref": (
        "warning",
        r"use prose/Markdown links instead of TeX \eqref navigation",
        True,
    ),
    r"\sb": ("warning", r"use a braced _{...} script", True),
    r"\sp": ("warning", r"use a braced ^{...} script", True),
    r"\mbox": ("warning", r"use \text{...} for ordinary equation text", True),
    r"\hbox": (
        "warning",
        r"use \text{...} for ordinary equation text; reserve \hbox for deliberate box-level needs",
        True,
    ),
    r"\rm": ("warning", r"use scoped \mathrm{...}", True),
    r"\tt": ("warning", r"use scoped \mathtt{...}", True),
    r"\bf": ("warning", r"use scoped \mathbf{...}", True),
    r"\medspace": (
        "error",
        r"use \mkern4mu when a medium math space is actually justified",
        False,
    ),
    r"\thickspace": (
        "error",
        r"use \mkern5mu when a thick math space is actually justified",
        False,
    ),
}

# Markdown can consume the backslash from these TeX control-symbol forms before
# the math renderer sees them.  Named forms are intentionally suggested.
CONTROL_SYMBOLS: dict[str, str] = {
    r"\,": r"use \thinspace",
    r"\:": r"use \mkern4mu when that spacing is justified",
    r"\>": r"use \mkern4mu when that spacing is justified",
    r"\;": r"use \mkern5mu when that spacing is justified",
    r"\!": r"use \negthinspace",
    r"\{": r"use \lbrace",
    r"\}": r"use \rbrace",
    r"\|": r"use \lVert/\rVert or another named vertical delimiter",
    r"\*": r"use \ast or another semantic operator",
    r"\_": r"if a literal underscore is intended, use a protected math escape hatch; otherwise use ordinary _{...} for a script",
    r"\$": r"if a literal dollar is intended, use protected inline math or a fenced math block",
    r"\%": r"if a literal percent sign is intended, use a protected math escape hatch",
    r"\#": r"if a literal hash is intended, use a protected math escape hatch",
    r"\&": r"if a literal ampersand is intended, use a protected math escape hatch",
}


@dataclass(frozen=True)
class MathSegment:
    text: str
    line: int
    column: int  # 1-based column of first character inside math container
    kind: str  # inline, protected, display, fence


@dataclass(frozen=True)
class InlineRange:
    start: int  # inclusive content index in line, 0-based
    end: int  # exclusive content index
    protected: bool = False


# ---------------------------------------------------------------------------
# Small Markdown helpers
# ---------------------------------------------------------------------------


def _is_escaped(text: str, index: int) -> bool:
    backslashes = 0
    i = index - 1
    while i >= 0 and text[i] == "\\":
        backslashes += 1
        i -= 1
    return backslashes % 2 == 1


def _is_unicode_whitespace(ch: str | None) -> bool:
    return ch is None or ch.isspace()


def _is_unicode_punctuation(ch: str | None) -> bool:
    if ch is None:
        return False
    cat = unicodedata.category(ch)
    return cat.startswith("P") or cat.startswith("S")


def underscore_open_close(text: str, index: int) -> tuple[bool, bool]:
    """Approximate CommonMark/GFM delimiter flanking for a single underscore.

    This deliberately models the contextual property Luna was missing: an
    underscore after punctuation (for example ``)_{...}`` or ``}_{...}``) can
    be both an opener and a closer, while ``x_{...}`` is not an opener.
    """

    prev = text[index - 1] if index > 0 else None
    nxt = text[index + 1] if index + 1 < len(text) else None

    left_flanking = (
        not _is_unicode_whitespace(nxt)
        and (
            not _is_unicode_punctuation(nxt)
            or _is_unicode_whitespace(prev)
            or _is_unicode_punctuation(prev)
        )
    )
    right_flanking = (
        not _is_unicode_whitespace(prev)
        and (
            not _is_unicode_punctuation(prev)
            or _is_unicode_whitespace(nxt)
            or _is_unicode_punctuation(nxt)
        )
    )

    can_open = left_flanking and (
        not right_flanking or _is_unicode_punctuation(prev)
    )
    can_close = right_flanking and (
        not left_flanking or _is_unicode_punctuation(nxt)
    )
    return can_open, can_close


def mask_inline_code(line: str) -> str:
    """Replace inline-code content with spaces while preserving line length.

    The implementation handles arbitrary backtick run lengths on one line,
    which is sufficient for linting repository prose and avoids interpreting
    literal TeX examples inside code spans as real mathematics.
    """

    chars = list(line)
    i = 0
    while i < len(line):
        if line[i] != "`":
            i += 1
            continue
        # GitHub protected inline math uses $`...`$.  Its backticks are math
        # delimiters, not ordinary Markdown code-span delimiters.
        if i > 0 and line[i - 1] == "$":
            protected_close = line.find("`$", i + 1)
            if protected_close >= 0:
                i = protected_close + 2
                continue
        j = i
        while j < len(line) and line[j] == "`":
            j += 1
        run = line[i:j]
        close = line.find(run, j)
        if close < 0:
            i = j
            continue
        for k in range(i, close + len(run)):
            if chars[k] not in "\r\n":
                chars[k] = " "
        i = close + len(run)
    return "".join(chars)


def protected_math_ranges(line: str) -> list[InlineRange]:
    """Find GitHub protected inline math of the form $`...`$."""

    ranges: list[InlineRange] = []
    i = 0
    while True:
        start = line.find("$`", i)
        if start < 0:
            break
        if _is_escaped(line, start):
            i = start + 2
            continue
        close = line.find("`$", start + 2)
        if close < 0:
            break
        ranges.append(InlineRange(start + 2, close, protected=True))
        i = close + 2
    return ranges


def mask_ranges(text: str, ranges: Sequence[tuple[int, int]]) -> str:
    chars = list(text)
    for start, end in ranges:
        for i in range(max(0, start), min(len(chars), end)):
            if chars[i] not in "\r\n":
                chars[i] = " "
    return "".join(chars)


def find_inline_math_ranges(line: str) -> list[InlineRange]:
    """Find ordinary $...$ ranges after code/protected spans have been masked."""

    code_masked = mask_inline_code(line)
    protected = protected_math_ranges(code_masked)
    masked = mask_ranges(
        code_masked,
        [(r.start - 2, r.end + 2) for r in protected],
    )

    ranges: list[InlineRange] = list(protected)
    i = 0
    while i < len(masked):
        if masked[i] != "$" or _is_escaped(masked, i):
            i += 1
            continue
        if i + 1 < len(masked) and masked[i + 1] == "$":
            i += 2
            continue
        j = i + 1
        while j < len(masked):
            if masked[j] == "$" and not _is_escaped(masked, j):
                if j + 1 < len(masked) and masked[j + 1] == "$":
                    j += 2
                    continue
                ranges.append(InlineRange(i + 1, j, protected=False))
                i = j + 1
                break
            j += 1
        else:
            i += 1
    ranges.sort(key=lambda r: r.start)
    return ranges


def _inside(index: int, ranges: Sequence[InlineRange], *, protected: bool | None = None) -> bool:
    for r in ranges:
        if r.start <= index < r.end and (protected is None or r.protected == protected):
            return True
    return False


def _looks_like_table_row(line: str) -> bool:
    stripped = line.strip()
    if not stripped or "|" not in stripped:
        return False
    # A conservative heuristic: leading/trailing pipe, or at least two pipes.
    return stripped.startswith("|") or stripped.endswith("|") or stripped.count("|") >= 2


def _excerpt(line: str, limit: int = 180) -> str:
    s = line.rstrip("\r\n")
    if len(s) <= limit:
        return s
    return s[: limit - 1] + "…"


def mask_tex_text_arguments(text: str) -> str:
    r"""Mask literal-text arguments where punctuation is not math syntax.

    ``\text{...}`` and ``\mathtt{...}`` may legitimately contain characters
    such as ``*``, ``|``, ``<`` or ``>`` as text/code-like content. The mask is
    brace-aware and length-preserving so diagnostic columns remain accurate.
    """

    chars = list(text)
    pattern = re.compile(r"\\(?:text|mathtt)\{")
    pos = 0
    while True:
        m = pattern.search(text, pos)
        if not m:
            break
        open_brace = m.end() - 1
        depth = 1
        i = open_brace + 1
        while i < len(text) and depth:
            if text[i] == "{" and not _is_escaped(text, i):
                depth += 1
            elif text[i] == "}" and not _is_escaped(text, i):
                depth -= 1
            if depth:
                chars[i] = " "
            i += 1
        pos = max(i, m.end())

    # In an array preamble, a raw vertical bar is TeX layout syntax rather than
    # a mathematical delimiter, for example ``\begin{array}{c|c}``. Do not
    # diagnose that structural use as MATH001.
    for m in re.finditer(r"\\begin\{array\}\{([^{}]*)\}", text):
        start, end = m.span(1)
        for i in range(start, end):
            chars[i] = " "

    return "".join(chars)


# ---------------------------------------------------------------------------
# TeX/math lint rules
# ---------------------------------------------------------------------------


def lint_math_segment(path: str, seg: MathSegment, source_line: str) -> list[Finding]:
    findings: list[Finding] = []
    text = seg.text
    syntax_text = mask_tex_text_arguments(text)
    base = seg.column

    def add(
        offset: int,
        severity: str,
        rule: str,
        message: str,
        *,
        style_only: bool = False,
    ) -> None:
        findings.append(
            Finding(
                path=path,
                line=seg.line,
                column=max(1, base + offset),
                severity=severity,
                rule=rule,
                message=message,
                excerpt=_excerpt(source_line),
                style_only=style_only,
            )
        )

    # Raw | is never the canonical project form.  \| is handled separately so
    # it does not create duplicate raw-pipe reports.
    for m in re.finditer(r"(?<!\\)\|", syntax_text):
        add(
            m.start(),
            "error",
            "MATH001",
            r"raw `|` in math; use \mid, \lvert/\rvert, \lVert/\rVert, or another semantic named form",
        )

    # Control-symbol escapes that Markdown can consume in ordinary dollar math.
    # Protected inline math and fenced math exist specifically to preserve raw
    # TeX punctuation when no clean named alternative exists.
    literal_escape_hatch_tokens = {r"\_", r"\$", r"\%", r"\#", r"\&"}
    for token, replacement in CONTROL_SYMBOLS.items():
        start = 0
        while True:
            idx = syntax_text.find(token, start)
            if idx < 0:
                break
            if seg.kind in {"inline", "display"}:
                add(idx, "error", "MATH002", f"avoid `{token}` here; {replacement}")
            elif token not in literal_escape_hatch_tokens:
                add(
                    idx,
                    "warning",
                    "MATH002",
                    f"`{token}` is protected from Markdown here, but prefer the project's named portable form when equivalent; {replacement}",
                    style_only=True,
                )
            start = idx + len(token)

    # Longer forbidden macro names first so \operatorname* is not swallowed by
    # the \operatorname check.
    for macro in sorted(FORBIDDEN_MACROS, key=len, reverse=True):
        sev, explanation, style_only = FORBIDDEN_MACROS[macro]
        # Command boundary: do not match \ref inside a longer command.
        suffix_guard = r"(?!\*)" if macro == r"\operatorname" else ""
        pattern = re.compile(re.escape(macro) + suffix_guard + r"(?![A-Za-z])")
        for m in pattern.finditer(text):
            add(m.start(), sev, "MATH003", explanation, style_only=style_only)

    for m in re.finditer(r"(?<!\\)\*", syntax_text):
        add(
            m.start(),
            "warning",
            "MATH007",
            r"literal `*` in math; use \ast, \star, \cdot, \times, or juxtaposition according to meaning",
        )

    for m in re.finditer(r"[<>]", syntax_text):
        add(
            m.start(),
            "warning",
            "MATH008",
            r"literal `<`/`>` in math; prefer \lt/\gt when that is the mathematical relation",
            style_only=True,
        )

    for m in re.finditer(r":=", syntax_text):
        add(
            m.start(),
            "warning",
            "MATH009",
            r"use \coloneqq for a definition operator when `:=` has that meaning",
            style_only=True,
        )

    # Braces are canonical around all scripts.  Ignore escaped \_ and \^; \_
    # is already diagnosed by MATH002.
    for m in re.finditer(r"(?<!\\)([_^])(?!\{)", syntax_text):
        add(
            m.start(),
            "warning",
            "MATH010",
            f"brace the `{m.group(1)}` script as `{m.group(1)}{{...}}`",
            style_only=True,
        )

    for m in re.finditer(r"\\vert\s*[_^]", syntax_text):
        add(
            m.start(),
            "warning",
            "MATH011",
            r"for a restriction/evaluation bar, use \rvert_{...} rather than generic \vert_{...}",
            style_only=True,
        )

    # A grouped custom operator immediately followed by a script is both a
    # semantic/layout review point and a known place where GitHub can pair the
    # punctuation-adjacent underscore with a later underscore.  \nolimits
    # breaks that pattern when side-script placement is intended.
    custom_operator_script = re.compile(
        r"\\mathop\{\\mathrm\{[^{}]+\}\}(?!\\nolimits)(\s*)([_^])"
    )
    for m in custom_operator_script.finditer(syntax_text):
        add(
            m.start(),
            "warning",
            "MATH012",
            r"custom operator has a side script; if side placement is intended, prefer \mathop{\mathrm{...}}\nolimits_{...}",
        )

    # Typography opportunity, deliberately advisory rather than mechanical.
    if r"\int" in syntax_text:
        for m in re.finditer(r"\\mathrm\{d\}", syntax_text):
            prefix = syntax_text[: m.start()].rstrip()
            if not prefix.endswith(r"\thinspace"):
                add(
                    m.start(),
                    "info",
                    "MATH013",
                    r"integral contains \mathrm{d}; review whether an explicit \thinspace before the differential is appropriate",
                    style_only=True,
                )

    for m in re.finditer(r"\\begin\{align\*?\}", syntax_text):
        add(
            m.start(),
            "warning",
            "MATH016",
            r"avoid top-level align/align* in repository Markdown; use aligned inside $$...$$ when rendered alignment is intended",
        )

    return findings


# ---------------------------------------------------------------------------
# File parser/linter
# ---------------------------------------------------------------------------


def lint_cross_line_underscore_pairs(
    lines: Sequence[str], display_path: str
) -> list[Finding]:
    """Find math underscores that can pair across soft line breaks.

    The main pass checks same-line GFM emphasis hazards. CommonMark emphasis can
    also span a soft line break inside one paragraph, so this second lightweight
    pass checks cross-line pairs. It intentionally requires both endpoints to
    occur in ordinary inline math to keep the heuristic conservative.
    """

    findings: list[Finding] = []
    fence_char: str | None = None
    fence_len = 0
    display_open = False
    paragraph: list[tuple[int, int, bool, bool, str]] = []

    def flush() -> None:
        nonlocal paragraph
        for i, (line_a, pos_a, can_open, _can_close, excerpt_a) in enumerate(paragraph):
            if not can_open:
                continue
            for line_b, _pos_b, _open_b, can_close_b, _excerpt_b in paragraph[i + 1 :]:
                if line_b == line_a or not can_close_b:
                    continue
                findings.append(
                    Finding(
                        display_path,
                        line_a,
                        pos_a + 1,
                        "error",
                        "MATH004",
                        "underscores in ordinary inline math can form GFM emphasis across a soft line break; prefer a portable TeX restructuring or protected inline math when unavoidable",
                        excerpt_a,
                        False,
                    )
                )
                paragraph = []
                return
        paragraph = []

    for idx, raw_line in enumerate(lines):
        line_no = idx + 1
        line = raw_line.rstrip("\r\n")

        if fence_char is not None:
            stripped = line.lstrip(" ")
            if (
                re.match(rf"^{re.escape(fence_char)}{{{fence_len},}}\s*$", stripped)
                and len(line) - len(stripped) <= 3
            ):
                fence_char = None
                fence_len = 0
            flush()
            continue

        open_match = FENCE_OPEN.match(line)
        if open_match:
            marker = open_match.group(2)
            fence_char = marker[0]
            fence_len = len(marker)
            flush()
            continue

        masked_code = mask_inline_code(line)
        protected = protected_math_ranges(masked_code)
        protected_containers = [(r.start - 2, r.end + 2) for r in protected]
        work = mask_ranges(masked_code, protected_containers)

        if display_open:
            if "$$" in work:
                display_open = False
            flush()
            continue
        if "$$" in work:
            # Any display delimiter breaks this lightweight paragraph model.
            if work.count("$$") % 2 == 1:
                display_open = True
            flush()
            continue

        if not line.strip():
            flush()
            continue

        ranges = [r for r in find_inline_math_ranges(work) if not r.protected]
        for r in ranges:
            for pos in range(r.start, r.end):
                if work[pos] != "_" or _is_escaped(work, pos):
                    continue
                can_open, can_close = underscore_open_close(work, pos)
                paragraph.append(
                    (line_no, pos, can_open, can_close, _excerpt(raw_line))
                )

    flush()
    return findings


FENCE_OPEN = re.compile(r"^( {0,3})(`{3,}|~{3,})(.*)$")


def lint_file(path: Path, display_path: str, include_history_style: bool) -> list[Finding]:
    text = path.read_text(encoding="utf-8")
    lines = text.splitlines(keepends=True)
    findings: list[Finding] = []

    fence_char: str | None = None
    fence_len = 0
    fence_is_math = False
    display_open = False
    display_open_line = 0

    # Track standalone display delimiters so blank-line checks can be applied
    # after the main pass without complicating the state machine.
    standalone_display_lines: list[tuple[int, str]] = []  # (index, open|close)

    for idx, raw_line in enumerate(lines):
        line_no = idx + 1
        line = raw_line.rstrip("\r\n")

        # Fenced blocks.  Non-math fences are literal examples/code and should
        # not be inspected as repository math.
        if fence_char is not None:
            stripped = line.lstrip(" ")
            close_match = re.match(rf"^{re.escape(fence_char)}{{{fence_len},}}\s*$", stripped)
            if close_match and len(line) - len(stripped) <= 3:
                fence_char = None
                fence_len = 0
                fence_is_math = False
                continue
            if fence_is_math:
                seg = MathSegment(line, line_no, 1, "fence")
                findings.extend(lint_math_segment(display_path, seg, raw_line))
            continue

        open_match = FENCE_OPEN.match(line)
        if open_match:
            marker = open_match.group(2)
            info = open_match.group(3).strip()
            fence_char = marker[0]
            fence_len = len(marker)
            fence_is_math = bool(info) and info.split()[0].lower() == "math"
            continue

        # Mask ordinary inline code before any Markdown/math checks.
        masked_code = mask_inline_code(line)

        # Repository math delimiters \(...\) and \[...\] are noncanonical.
        for token in (r"\(", r"\)", r"\[", r"\]"):
            start = 0
            while True:
                pos = masked_code.find(token, start)
                if pos < 0:
                    break
                findings.append(
                    Finding(
                        display_path,
                        line_no,
                        pos + 1,
                        "warning",
                        "MATH014",
                        rf"repository Markdown should use `$...$` / `$$...$$`, not `{token}` delimiters",
                        _excerpt(raw_line),
                        True,
                    )
                )
                start = pos + len(token)

        # Project Markdown prefers a terminal backslash for an intentional hard
        # break instead of two spaces.  Ignore blank lines.
        content_no_nl = raw_line.rstrip("\r\n")
        if content_no_nl.strip() and content_no_nl.endswith("  "):
            findings.append(
                Finding(
                    display_path,
                    line_no,
                    len(content_no_nl) - 1,
                    "warning",
                    "MD001",
                    "two-space hard line break; use a terminal backslash when a hard break is intentional",
                    _excerpt(raw_line),
                    True,
                )
            )

        # Protected inline math is parsed first so its contents are not treated
        # as ordinary dollar math for Markdown-emphasis collision analysis.
        protected_ranges = protected_math_ranges(masked_code)
        for r in protected_ranges:
            seg = MathSegment(
                line[r.start : r.end],
                line_no,
                r.start + 1,
                "protected",
            )
            findings.extend(lint_math_segment(display_path, seg, raw_line))

        protected_container_ranges = [(r.start - 2, r.end + 2) for r in protected_ranges]
        work = mask_ranges(masked_code, protected_container_ranges)

        # Display math state.  Canonical delimiters are standalone `$$` lines.
        if display_open:
            if re.fullmatch(r"\s*\$\$\s*", work):
                display_open = False
                standalone_display_lines.append((idx, "close"))
                continue

            closing = work.find("$$")
            if closing >= 0:
                # Content before a non-standalone closer is still math.
                content = line[:closing]
                if content:
                    seg = MathSegment(content, line_no, 1, "display")
                    findings.extend(lint_math_segment(display_path, seg, raw_line))
                findings.append(
                    Finding(
                        display_path,
                        line_no,
                        closing + 1,
                        "warning",
                        "MATH006",
                        "put the closing `$$` delimiter on its own line",
                        _excerpt(raw_line),
                        True,
                    )
                )
                display_open = False
                # Any material after the closer is intentionally not reparsed;
                # the layout finding is sufficient and avoids cascading noise.
                continue

            # Exact GFM block constructs matter; punctuation alone does not.
            stripped = work.rstrip()
            hazard: str | None = None
            if re.fullmatch(r" {0,3}(?:=+|-+)[ \t]*", stripped):
                hazard = "Setext heading underline"
            elif re.match(r"^ {0,3}[-+*](?:[ \t]+|$)", stripped):
                hazard = "bullet-list marker"
            elif re.match(r"^ {0,3}\d{1,9}[.)](?:[ \t]+|$)", stripped):
                hazard = "ordered-list marker"
            elif re.match(r"^ {0,3}>($|[ \t])", stripped):
                hazard = "block-quote marker"
            elif re.match(r"^ {0,3}#{1,6}(?:[ \t]+|$)", stripped):
                hazard = "ATX-heading marker"
            elif re.fullmatch(r" {0,3}(?:\*\s*){3,}", stripped) or re.fullmatch(
                r" {0,3}(?:_\s*){3,}", stripped
            ):
                hazard = "thematic break"

            if hazard:
                findings.append(
                    Finding(
                        display_path,
                        line_no,
                        1,
                        "error",
                        "MATH005",
                        f"{hazard} is Markdown-active inside `$$...$$`; reflow the TeX source without changing rendered math",
                        _excerpt(raw_line),
                        False,
                    )
                )

            seg = MathSegment(line, line_no, 1, "display")
            findings.extend(lint_math_segment(display_path, seg, raw_line))
            continue

        # Not currently in display mode.
        if re.fullmatch(r"\s*\$\$\s*", work):
            display_open = True
            display_open_line = line_no
            standalone_display_lines.append((idx, "open"))
            continue

        # Same-line display block, or an opening delimiter with content after it.
        first_dd = work.find("$$")
        if first_dd >= 0:
            second_dd = work.find("$$", first_dd + 2)
            if second_dd >= 0:
                content = line[first_dd + 2 : second_dd]
                seg = MathSegment(content, line_no, first_dd + 3, "display")
                findings.extend(lint_math_segment(display_path, seg, raw_line))
                findings.append(
                    Finding(
                        display_path,
                        line_no,
                        first_dd + 1,
                        "warning",
                        "MATH006",
                        "put display `$$` delimiters on their own lines with blank lines around the block",
                        _excerpt(raw_line),
                        True,
                    )
                )
                if _looks_like_table_row(line):
                    findings.append(
                        Finding(
                            display_path,
                            line_no,
                            first_dd + 1,
                            "error",
                            "MATH015",
                            "do not put display mathematics in a Markdown table cell; use inline math or move the equation outside the table",
                            _excerpt(raw_line),
                            False,
                        )
                    )
                # Mask the whole display block before ordinary inline scanning.
                work = mask_ranges(work, [(first_dd, second_dd + 2)])
            else:
                content = line[first_dd + 2 :]
                if content.strip():
                    seg = MathSegment(content, line_no, first_dd + 3, "display")
                    findings.extend(lint_math_segment(display_path, seg, raw_line))
                findings.append(
                    Finding(
                        display_path,
                        line_no,
                        first_dd + 1,
                        "warning",
                        "MATH006",
                        "put the opening `$$` delimiter on its own line",
                        _excerpt(raw_line),
                        True,
                    )
                )
                display_open = True
                display_open_line = line_no
                work = work[:first_dd]

        # Ordinary inline math after code/protected/display regions are masked.
        inline_ranges = [r for r in find_inline_math_ranges(work) if not r.protected]
        for r in inline_ranges:
            seg = MathSegment(line[r.start : r.end], line_no, r.start + 1, "inline")
            findings.extend(lint_math_segment(display_path, seg, raw_line))

        # Detect leftover ordinary dollar delimiters after complete inline spans
        # have been masked. This catches common typos without pretending to be
        # a full Markdown math parser.
        complete_inline_containers = [(r.start - 1, r.end + 1) for r in inline_ranges]
        leftover = mask_ranges(work, complete_inline_containers)
        for pos, ch in enumerate(leftover):
            if ch != "$" or _is_escaped(leftover, pos):
                continue
            if (pos > 0 and leftover[pos - 1] == "$") or (
                pos + 1 < len(leftover) and leftover[pos + 1] == "$"
            ):
                continue
            findings.append(
                Finding(
                    display_path,
                    line_no,
                    pos + 1,
                    "warning",
                    "MATH017",
                    "unmatched ordinary `$` math delimiter; verify the inline math boundary or escape/protect a literal dollar",
                    _excerpt(raw_line),
                    False,
                )
            )

        # Contextual underscore-emphasis analysis.  This does not search for
        # literal "__".  It applies GFM-style opener/closer eligibility to all
        # underscores on the line and reports a possible pair when at least one
        # endpoint lies inside ordinary, unprotected inline math.
        if inline_ranges:
            code_and_protected_mask = mask_ranges(masked_code, protected_container_ranges)
            positions: list[tuple[int, bool, bool, bool]] = []
            for pos, ch in enumerate(code_and_protected_mask):
                if ch != "_" or _is_escaped(code_and_protected_mask, pos):
                    continue
                can_open, can_close = underscore_open_close(code_and_protected_mask, pos)
                in_math = _inside(pos, inline_ranges, protected=False)
                positions.append((pos, can_open, can_close, in_math))

            reported_pair = False
            for a_i, (a_pos, a_open, _a_close, a_math) in enumerate(positions):
                if not a_open:
                    continue
                for b_pos, _b_open, b_close, b_math in positions[a_i + 1 :]:
                    if not b_close or not (a_math or b_math):
                        continue
                    # At least one endpoint must be in ordinary math, and both
                    # should be on the same Markdown line.  This catches both
                    # within-span and cross-span pairings without claiming that
                    # every underscore is unsafe.
                    span_kind = "across inline-math spans" if a_math and b_math and not any(
                        r.start <= a_pos < r.end and r.start <= b_pos < r.end
                        for r in inline_ranges
                    ) else "involving inline math"
                    findings.append(
                        Finding(
                            display_path,
                            line_no,
                            a_pos + 1,
                            "error",
                            "MATH004",
                            f"underscores can form GFM emphasis {span_kind}; prefer a portable TeX restructuring or protected inline math when unavoidable",
                            _excerpt(raw_line),
                            False,
                        )
                    )
                    reported_pair = True
                    break
                if reported_pair:
                    break

    if display_open:
        findings.append(
            Finding(
                display_path,
                display_open_line,
                1,
                "error",
                "MATH006",
                "unclosed `$$` display-math block",
                _excerpt(lines[display_open_line - 1]) if lines else "",
                False,
            )
        )

    # Blank lines around canonical display blocks.
    for index, kind in standalone_display_lines:
        if kind == "open":
            if index > 0 and lines[index - 1].strip():
                findings.append(
                    Finding(
                        display_path,
                        index + 1,
                        1,
                        "warning",
                        "MATH006",
                        "keep a blank line before a `$$` display block",
                        _excerpt(lines[index]),
                        True,
                    )
                )
        else:
            if index + 1 < len(lines) and lines[index + 1].strip():
                findings.append(
                    Finding(
                        display_path,
                        index + 1,
                        1,
                        "warning",
                        "MATH006",
                        "keep a blank line after a `$$` display block",
                        _excerpt(lines[index]),
                        True,
                    )
                )

    findings.extend(lint_cross_line_underscore_pairs(lines, display_path))

    # Suppress style-only modernization in historical evidence unless explicitly
    # requested.  Rendering/parser hazards still surface there.
    is_history = display_path.replace("\\", "/").startswith("docs/history/")
    if is_history and not include_history_style:
        findings = [f for f in findings if not f.style_only]

    # De-duplicate identical findings and provide deterministic ordering.
    unique: dict[tuple[object, ...], Finding] = {}
    for f in findings:
        key = (f.path, f.line, f.column, f.severity, f.rule, f.message)
        unique[key] = f
    return sorted(
        unique.values(),
        key=lambda f: (f.path, f.line, f.column, -SEVERITY_ORDER[f.severity], f.rule),
    )


# ---------------------------------------------------------------------------
# Repository path discovery and CLI
# ---------------------------------------------------------------------------


def _git_root() -> Path | None:
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"],
            stderr=subprocess.DEVNULL,
            text=True,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return None
    return Path(out) if out else None


def _git_paths(args: Sequence[str], root: Path) -> list[Path]:
    try:
        proc = subprocess.run(
            ["git", *args],
            cwd=root,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        raise RuntimeError(f"git path discovery failed: {exc}") from exc
    items = [p for p in proc.stdout.decode("utf-8", errors="surrogateescape").split("\0") if p]
    return [root / p for p in items]


def discover_markdown_paths(explicit: Sequence[str], changed: bool) -> tuple[Path, list[Path]]:
    root = _git_root() or Path.cwd().resolve()

    if explicit:
        paths: list[Path] = []
        for item in explicit:
            p = Path(item)
            if not p.is_absolute():
                p = (Path.cwd() / p).resolve()
            if p.is_dir():
                paths.extend(sorted(q for q in p.rglob("*.md") if q.is_file()))
            elif p.is_file() and p.suffix.lower() == ".md":
                paths.append(p)
        return root, sorted(set(paths))

    if (root / ".git").exists() or _git_root() is not None:
        if changed:
            tracked = _git_paths(
                ["diff", "--name-only", "-z", "HEAD", "--", "*.md"], root
            )
            untracked = _git_paths(
                ["ls-files", "--others", "--exclude-standard", "-z", "--", "*.md"],
                root,
            )
            return root, sorted(set(tracked + untracked))
        return root, _git_paths(["ls-files", "-z", "--", "*.md"], root)

    return root, sorted(Path.cwd().rglob("*.md"))


def relative_display(path: Path, root: Path) -> str:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return path.as_posix()


def print_text(findings: Sequence[Finding]) -> None:
    for f in findings:
        print(
            f"{f.path}:{f.line}:{f.column}: {f.severity} {f.rule}: {f.message}"
        )
        if f.excerpt:
            print(f"    {f.excerpt}")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "paths",
        nargs="*",
        help="Markdown files/directories; default is all tracked *.md files",
    )
    parser.add_argument(
        "--changed",
        action="store_true",
        help="scan only changed/untracked Markdown files relative to HEAD",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="emit machine-readable JSON findings",
    )
    parser.add_argument(
        "--include-history-style",
        action="store_true",
        help="also emit style-only findings under docs/history/",
    )
    parser.add_argument(
        "--fail-on",
        choices=("error", "warning", "info", "none"),
        default="warning",
        help="minimum severity that produces a non-zero exit status (default: warning)",
    )
    parser.add_argument(
        "--list-rules",
        action="store_true",
        help="print rule identifiers and exit",
    )
    ns = parser.parse_args(argv)

    if ns.list_rules:
        for rule, description in sorted(RULES.items()):
            print(f"{rule}: {description}")
        return 0

    root, paths = discover_markdown_paths(ns.paths, ns.changed)
    findings: list[Finding] = []
    for path in paths:
        if not path.exists():
            continue
        display = relative_display(path, root)
        try:
            findings.extend(lint_file(path, display, ns.include_history_style))
        except UnicodeDecodeError as exc:
            findings.append(
                Finding(
                    display,
                    1,
                    1,
                    "error",
                    "MD001",
                    f"file is not valid UTF-8: {exc}",
                    "",
                    False,
                )
            )

    findings.sort(
        key=lambda f: (f.path, f.line, f.column, -SEVERITY_ORDER[f.severity], f.rule)
    )

    if ns.json:
        json.dump([f.as_dict() for f in findings], sys.stdout, indent=2)
        sys.stdout.write("\n")
    else:
        print_text(findings)
        counts = {
            sev: sum(1 for f in findings if f.severity == sev)
            for sev in ("error", "warning", "info")
        }
        print(
            f"\nScanned {len(paths)} Markdown file(s): "
            f"{counts['error']} error(s), {counts['warning']} warning(s), "
            f"{counts['info']} info finding(s)."
        )

    if ns.fail_on == "none":
        return 0
    threshold = SEVERITY_ORDER[ns.fail_on]
    return 1 if any(SEVERITY_ORDER[f.severity] >= threshold for f in findings) else 0


if __name__ == "__main__":
    raise SystemExit(main())
