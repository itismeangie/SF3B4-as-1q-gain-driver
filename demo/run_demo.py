#!/usr/bin/env python3
"""Run an offline, synthetic-only example using existing research functions."""

from __future__ import annotations

import argparse
import csv
from html import escape
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "01_pan_cancer_1q"))
from analyze_1q_pathway_signed_roles import (  # noqa: E402
    benjamini_hochberg,
    fisher_exact_one_sided_greater,
)

COUNT_FIELDS = (
    "activator_1q", "activator_other", "repressor_1q", "repressor_other"
)
OUTPUT_FIELDS = (
    "scenario", *COUNT_FIELDS, "activator_1q_fraction", "repressor_1q_fraction",
    "p_one_sided_greater", "q_bh", "data_type",
)
FIXTURE = ROOT / "demo/data/synthetic_counts.csv"


def read_counts(path: Path) -> list[dict]:
    """Validate a small count table before passing it to research functions."""
    rows = []
    names = set()
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        expected = ["scenario", *COUNT_FIELDS]
        if reader.fieldnames != expected:
            raise ValueError("Expected columns: " + ", ".join(expected))
        for line, raw in enumerate(reader, start=2):
            if None in raw or any(value is None for value in raw.values()):
                raise ValueError(f"Line {line}: incorrect number of fields")
            name = raw["scenario"].strip()
            if not name or name in names:
                raise ValueError(f"Line {line}: scenario must be non-empty and unique")
            row = {"scenario": name}
            for field in COUNT_FIELDS:
                value = raw[field].strip()
                if not value or not value.isascii() or not value.isdecimal():
                    raise ValueError(f"Line {line}: {field} must be a non-negative integer")
                row[field] = int(value)
            if not row["activator_1q"] + row["activator_other"]:
                raise ValueError(f"Line {line}: activator group must be non-empty")
            if not row["repressor_1q"] + row["repressor_other"]:
                raise ValueError(f"Line {line}: repressor group must be non-empty")
            rows.append(row)
            names.add(name)
    if not rows:
        raise ValueError("The count table must contain at least one scenario")
    return rows


def analyze(rows: list[dict]) -> list[dict]:
    results = []
    for row in rows:
        a, b, c, d = (row[field] for field in COUNT_FIELDS)
        results.append({
            **row,
            "activator_1q_fraction": a / (a + b),
            "repressor_1q_fraction": c / (c + d),
            "p_one_sided_greater": fisher_exact_one_sided_greater(a, b, c, d),
            "data_type": "synthetic_only_not_study_results",
        })
    qvals = benjamini_hochberg([row["p_one_sided_greater"] for row in results])
    for row, qval in zip(results, qvals):
        row["q_bh"] = qval
    return results


def render_preview(results: list[dict]) -> str:
    """A small accessible SVG, deliberately labelled as synthetic throughout."""
    height = 238 + 108 * len(results)
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="{height}" viewBox="0 0 1000 {height}" role="img" aria-labelledby="title desc">',
        '<title id="title">Synthetic count-table demonstration, not study results</title>',
        '<desc id="desc">Invented scenarios compare the fraction on 1q in two toy gene-role groups, with one-sided Fisher p-values and BH-adjusted q-values.</desc>',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<g font-family="Arial, Helvetica, sans-serif" fill="#20303d">',
        '<text x="32" y="38" font-size="13" font-weight="700" fill="#9a531c">SYNTHETIC DEMO · NOT STUDY RESULTS</text>',
        '<text x="32" y="74" font-size="28" font-weight="700">From count tables to statistical summaries</text>',
        '<text x="32" y="102" font-size="15" fill="#52626c">Existing research functions · toy scenarios · no downloads or controlled data</text>',
        '<rect x="264" y="124" width="14" height="14" fill="#206d94"/>',
        '<text x="285" y="136" font-size="14">Activator-like</text>',
        '<rect x="441" y="124" width="14" height="14" fill="#bb6b29"/>',
        '<text x="462" y="136" font-size="14">Repressor-like</text>',
        '<text x="785" y="136" font-size="14">One-sided p / BH q</text>',
    ]
    for index, row in enumerate(results):
        y = 166 + index * 108
        label = escape(row["scenario"])
        svg.append(f'<text x="32" y="{y + 26}" font-size="17" font-weight="700">{label}</text>')
        for field, offset, color in (
            ("activator_1q_fraction", 0, "#206d94"),
            ("repressor_1q_fraction", 30, "#bb6b29"),
        ):
            fraction = row[field]
            width = 480 * fraction
            svg.append(f'<rect x="264" y="{y + offset}" width="480" height="22" rx="3" fill="#edf1f4"/>')
            svg.append(f'<rect x="264" y="{y + offset}" width="{width:.3f}" height="22" rx="3" fill="{color}"/>')
            svg.append(f'<text x="{272 + width:.3f}" y="{y + offset + 16}" font-size="13">{fraction:.0%}</text>')
        svg.append(f'<text x="785" y="{y + 31}" font-size="16">{row["p_one_sided_greater"]:.3g} / {row["q_bh"]:.3g}</text>')
    bottom = 166 + len(results) * 108
    svg.extend([
        f'<text x="264" y="{bottom - 24}" font-size="13" fill="#52626c">Fraction on 1q within each toy group (0–100%)</text>',
        f'<text x="32" y="{bottom + 12}" font-size="13">Alternative: greater representation among activator-like genes; depletion is not tested in that direction.</text>',
        f'<text x="32" y="{bottom + 37}" font-size="13" fill="#9a531c">Illustration only: no biological conclusions or manuscript results are represented.</text>',
        '</g></svg>',
    ])
    return "\n".join(svg) + "\n"


def run_demo(output_dir: Path) -> list[dict]:
    """Always use the bundled invented fixture, never user-supplied research data."""
    results = analyze(read_counts(FIXTURE))
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "results.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS)
        writer.writeheader()
        for row in results:
            writer.writerow({
                key: format(value, ".12g") if isinstance(value, float) else value
                for key, value in row.items()
            })
    (output_dir / "preview.svg").write_text(render_preview(results), encoding="utf-8")
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "demo/output")
    args = parser.parse_args()
    try:
        run_demo(args.output_dir)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print("Synthetic demo only; not study results.")
    print(f"CSV: {args.output_dir / 'results.csv'}")
    print(f"SVG: {args.output_dir / 'preview.svg'}")


if __name__ == "__main__":
    main()
