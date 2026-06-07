#!/usr/bin/env python3
"""Summarize major functional/cancer gene categories across chromosome arms.

Outputs:
- results/arm_category_density.tsv
- results/category_1q_summary.tsv
"""

from __future__ import annotations

import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median
from typing import Dict, List, Set, Tuple

import analyze_1q_developmental_enrichment as base

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
OUT_DIR = BASE_DIR / "results"

ONCOKB_URL = "https://www.oncokb.org/api/v1/utils/allCuratedGenes"
ONCOKB_PATH = DATA_DIR / "oncokb_allCuratedGenes.json"
EXCLUDED_CYTOBAND_STAINS = {"acen", "gvar", "stalk"}

CATEGORY_SPECS: List[dict] = [
    {
        "id": "oncogenes",
        "label": "Oncogenes",
        "source": "oncokb",
        "oncokb_types": {"ONCOGENE", "ONCOGENE_AND_TSG"},
    },
    {
        "id": "tumor_suppressors",
        "label": "Tumor suppressors",
        "source": "oncokb",
        "oncokb_types": {"TSG", "ONCOGENE_AND_TSG"},
    },
    {
        "id": "transcription_factors",
        "label": "Transcription factors",
        "source": "go",
        "go_roots": ["GO:0140110"],  # transcription regulator activity
    },
    {
        "id": "protein_kinases",
        "label": "Protein kinases",
        "source": "go",
        "go_roots": ["GO:0004672"],  # protein kinase activity
    },
    {
        "id": "chromatin_remodeling",
        "label": "Chromatin remodeling",
        "source": "go",
        "go_roots": ["GO:0006338"],
    },
    {
        "id": "epigenetic_regulation",
        "label": "Epigenetic regulation",
        "source": "go",
        "go_roots": ["GO:0040029"],  # regulation of gene expression, epigenetic
    },
    {
        "id": "dna_repair",
        "label": "DNA repair",
        "source": "go",
        "go_roots": ["GO:0006281"],
    },
    {
        "id": "rna_splicing",
        "label": "RNA splicing",
        "source": "go",
        "go_roots": ["GO:0008380"],
    },
    {
        "id": "ubiquitin_ligases",
        "label": "Ubiquitin ligases",
        "source": "go",
        "go_roots": ["GO:0004842"],  # ubiquitin-protein transferase activity
    },
    {
        "id": "cell_surface_receptors",
        "label": "Cell-surface receptors",
        "source": "go",
        "go_roots": ["GO:0004888"],  # transmembrane signaling receptor activity
    },
    {
        "id": "developmental_process",
        "label": "Developmental process",
        "source": "go",
        "go_roots": ["GO:0032502"],
    },
    {
        "id": "stem_cell_programs",
        "label": "Stem-cell programs",
        "source": "go",
        "go_roots": [
            "GO:0019827",  # stem cell population maintenance
            "GO:0048863",  # stem cell differentiation
            "GO:0048864",  # stem cell development
        ],
    },
]


def arm_sort_key(arm: str):
    chrom = arm[:-1]
    arm_letter = arm[-1]
    if chrom == "X":
        c = 23
    elif chrom == "Y":
        c = 24
    else:
        c = int(chrom)
    return (c, 0 if arm_letter == "p" else 1)


def merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        last_s, last_e = merged[-1]
        if start <= last_e:
            merged[-1] = (last_s, max(last_e, end))
        else:
            merged.append((start, end))
    return merged


def load_euchromatin_arm_sizes_and_intervals(cytoband_path: Path):
    arm_intervals: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    with base.open_text(cytoband_path) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) < 5:
                continue
            chrom, start_s, end_s, band, stain = row
            if not chrom.startswith("chr"):
                continue
            chrom_label = chrom[3:]
            if chrom_label not in base.ALLOWED_CHR:
                continue
            if stain in EXCLUDED_CYTOBAND_STAINS:
                continue
            if not band or band[0] not in {"p", "q"}:
                continue
            arm = f"{chrom_label}{band[0]}"
            try:
                start = int(start_s)
                end = int(end_s)
            except ValueError:
                continue
            if end <= start:
                continue
            arm_intervals[arm].append((start, end))

    merged_intervals: Dict[str, List[Tuple[int, int]]] = {}
    arm_sizes_bp: Dict[str, int] = {}
    for arm, intervals in arm_intervals.items():
        merged = merge_intervals(intervals)
        size_bp = sum(end - start for start, end in merged)
        if size_bp <= 0:
            continue
        merged_intervals[arm] = merged
        arm_sizes_bp[arm] = size_bp
    return arm_sizes_bp, merged_intervals


def midpoint_in_intervals(pos: int, intervals: List[Tuple[int, int]]) -> bool:
    for start, end in intervals:
        if start <= pos < end:
            return True
    return False


def load_eligible_gene_to_arm(
    protein_coding: Dict[str, str],
    gene_locs: Dict[str, base.GeneLoc],
    centromere_q_start: Dict[str, int],
    arm_sizes_bp: Dict[str, int],
    euchromatin_intervals: Dict[str, List[Tuple[int, int]]],
):
    gene_to_arm: Dict[str, str] = {}
    all_counts = Counter()
    eligible = set(protein_coding).intersection(gene_locs)
    for gene_id in eligible:
        gloc = gene_locs[gene_id]
        arm = base.assign_arm(gloc, centromere_q_start)
        if not arm or arm not in arm_sizes_bp:
            continue
        intervals = euchromatin_intervals.get(arm)
        if not intervals or not midpoint_in_intervals(gloc.midpoint, intervals):
            continue
        gene_to_arm[gene_id] = arm
        all_counts[arm] += 1
    return gene_to_arm, all_counts


def build_go_category_gene_sets(gene2go_path: Path, category_terms: Dict[str, Set[str]]) -> Dict[str, Set[str]]:
    # term -> categories lookup for single-pass gene2go scan
    term_to_categories: Dict[str, List[str]] = defaultdict(list)
    for cat_id, terms in category_terms.items():
        for term in terms:
            term_to_categories[term].append(cat_id)

    go_category_genes: Dict[str, Set[str]] = {cat_id: set() for cat_id in category_terms}

    with base.open_text(gene2go_path) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if row[0] != "9606":
                continue
            go_id = row[2]
            if go_id not in term_to_categories:
                continue
            qualifier = row[4]
            if "NOT" in qualifier:
                continue
            gene_id = row[1]
            for cat_id in term_to_categories[go_id]:
                go_category_genes[cat_id].add(gene_id)

    return go_category_genes


def load_oncokb_category_gene_sets(path: Path) -> Dict[str, Set[str]]:
    with open(path, "r", encoding="utf-8") as f:
        genes = json.load(f)

    by_type: Dict[str, Set[str]] = defaultdict(set)
    for item in genes:
        gene_type = item.get("geneType")
        try:
            entrez = int(item.get("entrezGeneId", -1))
        except (TypeError, ValueError):
            continue
        if entrez <= 0:
            continue
        by_type[gene_type].add(str(entrez))

    out = {}
    for spec in CATEGORY_SPECS:
        if spec["source"] != "oncokb":
            continue
        cat_ids = spec["oncokb_types"]
        merged = set()
        for t in cat_ids:
            merged |= by_type.get(t, set())
        out[spec["id"]] = merged
    return out


def zscore(values: List[float]) -> List[float]:
    m = sum(values) / len(values)
    var = sum((x - m) ** 2 for x in values) / len(values)
    sd = math.sqrt(var)
    if sd == 0:
        return [0.0 for _ in values]
    return [(x - m) / sd for x in values]


def main() -> int:
    DATA_DIR.mkdir(exist_ok=True)
    OUT_DIR.mkdir(exist_ok=True)

    files = base.url_file_map()
    for key, url in base.URLS.items():
        base.download_if_missing(url, files[key])
    base.download_if_missing(ONCOKB_URL, ONCOKB_PATH)

    accn_to_chr = base.parse_assembly_refseq_to_chr(files["assembly_report"])
    protein_coding = base.load_protein_coding_genes(files["gene_info"])
    gene_locs = base.load_gene_locations(files["gene2refseq"], accn_to_chr)

    centromere_q_start = base.load_centromere_boundaries(files["cytoband"])
    base.load_chrom_sizes(files["chrom_sizes"])  # kept for parity with prior downloads
    arm_sizes_bp, euchromatin_intervals = load_euchromatin_arm_sizes_and_intervals(files["cytoband"])

    gene_to_arm, all_counts = load_eligible_gene_to_arm(
        protein_coding, gene_locs, centromere_q_start, arm_sizes_bp, euchromatin_intervals
    )

    # Build GO descendants for requested categories.
    children = base.parse_go_graph(files["go_basic"])
    go_category_terms: Dict[str, Set[str]] = {}
    for spec in CATEGORY_SPECS:
        if spec["source"] != "go":
            continue
        terms = set()
        for root in spec["go_roots"]:
            terms |= base.descendants(children, root)
        go_category_terms[spec["id"]] = terms

    go_category_genes = build_go_category_gene_sets(files["gene2go"], go_category_terms)
    oncokb_category_genes = load_oncokb_category_gene_sets(ONCOKB_PATH)

    category_to_genes: Dict[str, Set[str]] = {}
    for spec in CATEGORY_SPECS:
        cat_id = spec["id"]
        if spec["source"] == "go":
            category_to_genes[cat_id] = go_category_genes.get(cat_id, set())
        else:
            category_to_genes[cat_id] = oncokb_category_genes.get(cat_id, set())

    arms = sorted(arm_sizes_bp, key=arm_sort_key)

    arm_rows = []
    summary_rows = []

    for spec in CATEGORY_SPECS:
        cat_id = spec["id"]
        cat_label = spec["label"]
        gene_ids = category_to_genes[cat_id]

        # Restrict to protein-coding genes that have assigned chromosome arm.
        eligible_genes = [gid for gid in gene_ids if gid in gene_to_arm]

        counts = Counter(gene_to_arm[gid] for gid in eligible_genes)
        densities = [
            (counts.get(arm, 0) / arm_sizes_bp[arm]) * 1_000_000 for arm in arms
        ]
        zscores = zscore(densities)

        rank_order = sorted(
            ((arm, densities[i]) for i, arm in enumerate(arms)),
            key=lambda x: x[1],
            reverse=True,
        )
        arm_rank = {arm: i + 1 for i, (arm, _) in enumerate(rank_order)}

        for i, arm in enumerate(arms):
            arm_rows.append(
                {
                    "arm": arm,
                    "arm_size_mb": arm_sizes_bp[arm] / 1_000_000,
                    "total_protein_coding_genes_in_arm": all_counts.get(arm, 0),
                    "category_id": cat_id,
                    "category_label": cat_label,
                    "category_gene_count": counts.get(arm, 0),
                    "category_genes_per_mb": densities[i],
                    "category_density_zscore": zscores[i],
                    "arm_rank_within_category": arm_rank[arm],
                    "arms_in_category": len(arms),
                    "category_gene_fraction_in_arm": (
                        counts.get(arm, 0) / all_counts.get(arm, 1)
                        if all_counts.get(arm, 0) > 0
                        else float("nan")
                    ),
                }
            )

        one_q_idx = arms.index("1q")
        summary_rows.append(
            {
                "category_id": cat_id,
                "category_label": cat_label,
                "eligible_category_genes": len(eligible_genes),
                "one_q_gene_count": counts.get("1q", 0),
                "one_q_genes_per_mb": densities[one_q_idx],
                "one_q_density_zscore": zscores[one_q_idx],
                "one_q_rank": arm_rank["1q"],
                "arms_total": len(arms),
                "mean_genes_per_mb": sum(densities) / len(densities),
                "median_genes_per_mb": median(densities),
            }
        )

    arm_out = OUT_DIR / "arm_category_density.tsv"
    with open(arm_out, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "arm",
                "arm_size_mb",
                "total_protein_coding_genes_in_arm",
                "category_id",
                "category_label",
                "category_gene_count",
                "category_genes_per_mb",
                "category_density_zscore",
                "arm_rank_within_category",
                "arms_in_category",
                "category_gene_fraction_in_arm",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in arm_rows:
            writer.writerow(
                {
                    k: (
                        f"{v:.6f}"
                        if isinstance(v, float) and not math.isnan(v)
                        else ("nan" if isinstance(v, float) else v)
                    )
                    for k, v in row.items()
                }
            )

    summary_out = OUT_DIR / "category_1q_summary.tsv"
    with open(summary_out, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "category_id",
                "category_label",
                "eligible_category_genes",
                "one_q_gene_count",
                "one_q_genes_per_mb",
                "one_q_density_zscore",
                "one_q_rank",
                "arms_total",
                "mean_genes_per_mb",
                "median_genes_per_mb",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in summary_rows:
            writer.writerow(
                {
                    k: (
                        f"{v:.6f}"
                        if isinstance(v, float) and not math.isnan(v)
                        else ("nan" if isinstance(v, float) else v)
                    )
                    for k, v in row.items()
                }
            )

    print(f"Wrote {arm_out}")
    print(f"Wrote {summary_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
