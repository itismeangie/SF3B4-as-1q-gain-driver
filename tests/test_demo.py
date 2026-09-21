"""Offline tests for the demo and selected, unchanged research functions."""

import csv
from itertools import product
from math import comb
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET

from demo.run_demo import (
    ROOT, analyze, benjamini_hochberg, fisher_exact_one_sided_greater,
    read_counts, render_preview, run_demo,
)
from analyze_1q_pathway_signed_roles import parse_kegg_signed_roles


class StatisticsTests(unittest.TestCase):
    def test_fisher_against_exact_combinatorics(self):
        for a, b, c, d in product(range(5), repeat=4):
            with self.subTest(table=(a, b, c, d)):
                row1, row2, col1 = a + b, c + d, a + c
                if row1 == 0 or row2 == 0:
                    expected = 1.0
                else:
                    numerator = sum(
                        comb(row1, x) * comb(row2, col1 - x)
                        for x in range(a, min(row1, col1) + 1)
                    )
                    expected = numerator / comb(row1 + row2, col1)
                self.assertAlmostEqual(
                    fisher_exact_one_sided_greater(a, b, c, d), expected, places=12
                )

    def test_bh_preserves_input_order(self):
        actual = benjamini_hochberg([0.04, 0.01, 0.03, 0.2])
        for value, expected in zip(actual, [0.16 / 3, 0.04, 0.16 / 3, 0.2]):
            self.assertAlmostEqual(value, expected)
        self.assertEqual(benjamini_hochberg([]), [])
        self.assertEqual(benjamini_hochberg([0.0, 0.0, 1.0]), [0.0, 0.0, 1.0])

    def test_toy_signed_graph(self):
        # Numeric tokens are toy parser IDs, not biological gene assignments.
        kgml = '''<pathway>
          <entry id="1" name="hsa:1" type="gene"/>
          <entry id="2" name="hsa:2" type="gene"/>
          <entry id="3" name="hsa:3" type="gene"/>
          <entry id="4" name="hsa:4" type="gene"/>
          <relation entry1="1" entry2="4"><subtype name="activation"/></relation>
          <relation entry1="2" entry2="4"><subtype name="inhibition"/></relation>
          <relation entry1="3" entry2="1"><subtype name="activation"/></relation>
          <relation entry1="3" entry2="2"><subtype name="inhibition"/></relation>
        </pathway>'''
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "toy.xml"
            path.write_text(kgml, encoding="utf-8")
            self.assertEqual(
                parse_kegg_signed_roles(path),
                ({"1", "2", "3", "4"}, {"1"}, {"2"}, {"3"}),
            )


class DemoTests(unittest.TestCase):
    def test_offline_deterministic_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            with patch("urllib.request.urlopen", side_effect=AssertionError("Network forbidden")):
                results = run_demo(output)
                first = {p.name: p.read_bytes() for p in output.iterdir()}
                run_demo(output)
                second = {p.name: p.read_bytes() for p in output.iterdir()}
            self.assertEqual(first, second)
            self.assertEqual(set(first), {"results.csv", "preview.svg"})
            self.assertEqual(len(results), 3)
            self.assertTrue(all(r["data_type"] == "synthetic_only_not_study_results" for r in results))
            ET.fromstring(first["preview.svg"])
            self.assertIn(b"NOT STUDY RESULTS", first["preview.svg"])

    def test_reference_values(self):
        results = analyze(read_counts(ROOT / "demo/data/synthetic_counts.csv"))
        with (ROOT / "demo/example_output/results.csv").open(newline="", encoding="utf-8") as handle:
            reference = list(csv.DictReader(handle))
        self.assertEqual(len(results), len(reference))
        for actual, expected in zip(results, reference):
            self.assertEqual(actual["scenario"], expected["scenario"])
            for key in ("activator_1q_fraction", "repressor_1q_fraction", "p_one_sided_greater", "q_bh"):
                self.assertAlmostEqual(actual[key], float(expected[key]), places=10)
        self.assertAlmostEqual(results[0]["p_one_sided_greater"], 0.0115070687795, places=10)
        self.assertAlmostEqual(results[0]["q_bh"], 0.0345212063385, places=10)

    def test_invalid_tables(self):
        header = "scenario,activator_1q,activator_other,repressor_1q,repressor_other\n"
        invalid = [
            "wrong,columns\n", header,
            header + "negative,-1,2,3,4\n",
            header + "fraction,1.2,2,3,4\n",
            header + "empty,0,0,3,4\n",
            header + "empty,1,2,0,0\n",
            header + "same,1,2,3,4\nsame,2,3,4,5\n",
            header + ",1,2,3,4\n",
            header + "missing,1,2,3\n",
            header + "extra,1,2,3,4,5\n",
        ]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "invalid.csv"
            for content in invalid:
                with self.subTest(content=content):
                    path.write_text(content, encoding="utf-8")
                    with self.assertRaises(ValueError):
                        read_counts(path)

    def test_svg_escapes_labels(self):
        rows = read_counts(ROOT / "demo/data/synthetic_counts.csv")
        rows[0]["scenario"] = "Toy <A> & B"
        svg = render_preview(analyze(rows))
        ET.fromstring(svg)
        self.assertIn("Toy &lt;A&gt; &amp; B", svg)


if __name__ == "__main__":
    unittest.main()
