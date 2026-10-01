#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "tools" / "lint_markdown_math.py"
spec = importlib.util.spec_from_file_location("lint_markdown_math", SCRIPT)
assert spec and spec.loader
lint = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = lint
spec.loader.exec_module(lint)


class MarkdownMathLintTests(unittest.TestCase):
    def scan(self, source: str, name: str = "docs/test.md"):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "test.md"
            p.write_text(source, encoding="utf-8")
            return lint.lint_file(p, name, include_history_style=False)

    def rules(self, source: str, name: str = "docs/test.md"):
        return [f.rule for f in self.scan(source, name)]

    def test_safe_braced_subscript_does_not_trigger_emphasis(self):
        self.assertNotIn("MATH004", self.rules("Value $x_{i}$ and $\\Gamma_{D}$.\n"))

    def test_grouped_operator_subscript_triggers_contextual_emphasis(self):
        rs = self.rules(r"Value $\mathop{\mathrm{tr}}_{\Gamma_{D}} v$." + "\n")
        self.assertIn("MATH004", rs)
        self.assertIn("MATH012", rs)

    def test_nolimits_operator_is_safe(self):
        rs = self.rules(r"Value $\mathop{\mathrm{tr}}\nolimits_{\Gamma_{D}} v$." + "\n")
        self.assertNotIn("MATH004", rs)
        self.assertNotIn("MATH012", rs)

    def test_underscores_pair_across_separate_math_spans(self):
        rs = self.rules(
            r"Terms $(\nabla y,\nabla v)_{\Omega}$, $-\langle f,v\rangle$, $-(u,v)_{\Omega}$."
            + "\n"
        )
        self.assertIn("MATH004", rs)

    def test_underscores_can_pair_across_soft_line_break(self):
        rs = self.rules(
            r"Terms $(\nabla y,\nabla v)_{\Omega}$," + "\n"
            r"and $-(u,v)_{\Omega}$." + "\n"
        )
        self.assertIn("MATH004", rs)

    def test_protected_math_is_not_treated_as_ordinary_emphasis(self):
        rs = self.rules(r"Value $`(u,v)_{\Omega}`$ and $`(p,q)_{\Omega}`$." + "\n")
        self.assertNotIn("MATH004", rs)

    def test_raw_pipe_in_math_is_reported(self):
        self.assertIn("MATH001", self.rules(r"Use $u|_{\Gamma}$." + "\n"))

    def test_array_column_rule_is_not_treated_as_math_delimiter(self):
        rs = self.rules(r"Use $\begin{array}{c|c} a & b \\ c & d \end{array}$." + "\n")
        self.assertNotIn("MATH001", rs)

    def test_named_restriction_bar_is_safe(self):
        self.assertNotIn("MATH001", self.rules(r"Use $u\rvert_{\Gamma}$." + "\n"))

    def test_minus_without_space_is_not_a_list_hazard(self):
        rs = self.rules(
            "$$\n"
            r"\begin{aligned}" + "\n"
            r"-\nabla\cdot(\kappa\nabla y) &= f," + "\n"
            r"\end{aligned}" + "\n"
            "$$\n"
        )
        self.assertNotIn("MATH005", rs)

    def test_minus_with_space_is_a_list_hazard(self):
        rs = self.rules("$$\n- \\nabla y\n$$\n")
        self.assertIn("MATH005", rs)

    def test_standalone_equals_is_setext_hazard(self):
        self.assertIn("MATH005", self.rules("$$\nA x\n=\nb\n$$\n"))

    def test_equals_with_other_characters_is_not_setext_hazard(self):
        self.assertNotIn("MATH005", self.rules("$$\n=A x\n$$\n"))

    def test_non_math_code_span_is_ignored(self):
        rs = self.rules(r"Literal code: `u|_{\Gamma}` and `\operatorname{foo}`." + "\n")
        self.assertNotIn("MATH001", rs)
        self.assertNotIn("MATH003", rs)

    def test_punctuation_spacing_command_is_reported_in_dollar_math(self):
        self.assertIn("MATH002", self.rules(r"Value $a\,b$." + "\n"))

    def test_literal_dollar_escape_requires_protection(self):
        self.assertIn("MATH002", self.rules(r"Value $\sqrt{\$4}$." + "\n"))

    def test_literal_dollar_escape_is_allowed_in_protected_math(self):
        findings = self.scan(r"Value $`\sqrt{\$4}`$." + "\n")
        self.assertFalse(any(f.rule == "MATH002" and f.severity == "error" for f in findings))

    def test_text_argument_punctuation_is_not_treated_as_math_syntax(self):
        rs = self.rules(r"Value $x \text{ if a*b | c<d}$." + "\n")
        self.assertNotIn("MATH001", rs)
        self.assertNotIn("MATH007", rs)
        self.assertNotIn("MATH008", rs)

    def test_same_line_display_delimiters_are_reported(self):
        self.assertIn("MATH006", self.rules("Before.\n\n$$ a=b $$\n"))

    def test_top_level_align_is_reported(self):
        rs = self.rules("$$\n\\begin{align}\na &= b\n\\end{align}\n$$\n")
        self.assertIn("MATH016", rs)

    def test_unmatched_dollar_is_reported(self):
        self.assertIn("MATH017", self.rules("Broken $x_{i}.\n"))

    def test_hbox_is_reported(self):
        self.assertIn("MATH003", self.rules(r"Value $x\hbox{ if }y$." + "\n"))

    def test_operatorname_star_is_not_double_reported_as_operatorname(self):
        findings = self.scan(r"Value $\operatorname*{argmin}_{x} f(x)$." + "\n")
        msgs = [f.message for f in findings if f.rule == "MATH003" and "operatorname" in f.message]
        self.assertEqual(1, len(msgs))


    def test_standalone_plus_is_context_sensitive_warning_not_error(self):
        findings = self.scan("$$\n+\n$$\n")
        matches = [f for f in findings if f.rule == "MATH005"]
        self.assertEqual(1, len(matches))
        self.assertEqual("warning", matches[0].severity)

    def test_standalone_ordered_marker_is_context_sensitive_warning_not_error(self):
        findings = self.scan("$$\n0.\n$$\n")
        matches = [f for f in findings if f.rule == "MATH005"]
        self.assertEqual(1, len(matches))
        self.assertEqual("warning", matches[0].severity)

    def test_plus_with_following_content_is_definite_list_hazard(self):
        findings = self.scan("$$\n+ x\n$$\n")
        matches = [f for f in findings if f.rule == "MATH005"]
        self.assertEqual(1, len(matches))
        self.assertEqual("error", matches[0].severity)

    def test_quoted_display_math_treats_quote_prefix_as_container(self):
        findings = self.scan(
            "> Quoted equation:\n"
            ">\n"
            "> $$\n"
            "> \\begin{aligned}\n"
            "> x &= y, \\\\n"
            "> z &= w.\n"
            "> \\end{aligned}\n"
            "> $$\n"
            ">\n"
            "> Following text.\n"
        )
        self.assertNotIn("MATH005", [f.rule for f in findings])
        self.assertNotIn("MATH006", [f.rule for f in findings])

    def test_quoted_math_fence_treats_quote_prefix_as_container(self):
        findings = self.scan(
            "> ```math\n"
            "> \\mathrm{is\\_feasible}\n"
            "> ```\n"
        )
        self.assertFalse(any(f.rule == "MATH002" for f in findings))
        self.assertNotIn("MATH006", [f.rule for f in findings])

    def test_extra_quote_marker_inside_quoted_display_is_still_hazard(self):
        findings = self.scan(">\n> $$\n> >0\n> $$\n>\n")
        matches = [f for f in findings if f.rule == "MATH005"]
        self.assertEqual(1, len(matches))
        self.assertEqual("error", matches[0].severity)

    def test_literal_underscore_escape_is_reported_in_ordinary_math(self):
        findings = self.scan(r"Value $\mathrm{is\_feasible}$." + "\n")
        self.assertTrue(any(f.rule == "MATH002" and f.severity == "error" for f in findings))

    def test_literal_underscore_escape_is_allowed_in_protected_math(self):
        findings = self.scan(r"Value $`\mathrm{is\_feasible}`$." + "\n")
        self.assertFalse(any(f.rule == "MATH002" for f in findings))

    def test_literal_underscore_escape_is_allowed_in_math_fence(self):
        findings = self.scan("```math\n\\mathrm{is\\_feasible}\n```\n")
        self.assertFalse(any(f.rule == "MATH002" for f in findings))

    def test_history_suppresses_style_only_but_not_rendering_hazards(self):
        findings = self.scan(
            r"Old $x_i := y|_\Gamma$." + "\n",
            "docs/history/old.md",
        )
        rules = [f.rule for f in findings]
        self.assertNotIn("MATH009", rules)
        self.assertNotIn("MATH010", rules)
        self.assertIn("MATH001", rules)


if __name__ == "__main__":
    unittest.main()
