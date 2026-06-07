#!/usr/bin/env python3
"""Test whether 1q is an outlier for oncogenic pathways (gene density by arm size).

Pathway source: KEGG REST (human pathways).
Outputs:
- results/oncogenic_pathway_arm_density.tsv
- results/oncogenic_pathway_1q_summary.tsv
"""

from __future__ import annotations

import csv
import math
import re
import urllib.request
from collections import Counter
from pathlib import Path
from typing import Dict, List, Set, Tuple

import analyze_1q_developmental_enrichment as base

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
OUT_DIR = BASE_DIR / "results"

# Canonical oncogenic signaling and cancer-relevant pathways.
PATHWAYS = [
    ("hsa05200", "Pathways in cancer"),
    ("hsa04010", "MAPK signaling pathway"),
    ("hsa04151", "PI3K-Akt signaling pathway"),
    ("hsa04014", "Ras signaling pathway"),
    ("hsa04012", "ErbB signaling pathway"),
    ("hsa04310", "Wnt signaling pathway"),
    ("hsa04330", "Notch signaling pathway"),
    ("hsa04340", "Hedgehog signaling pathway"),
    ("hsa04350", "TGF-beta signaling pathway"),
    ("hsa04390", "Hippo signaling pathway"),
    ("hsa04630", "JAK-STAT signaling pathway"),
    ("hsa04150", "mTOR signaling pathway"),
    ("hsa04110", "Cell cycle"),
    ("hsa04115", "p53 signaling pathway"),
    ("hsa04210", "Apoptosis"),
    ("hsa03440", "Homologous recombination"),
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


def benjamini_hochberg(pvals: List[float]) -> List[float]:
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    qvals = [1.0] * m
    prev = 1.0
    for rank, idx in enumerate(reversed(order), start=1):
        i = m - rank + 1
        p = pvals[idx]
        q = min(prev, (p * m) / i)
        qvals[idx] = q
        prev = q
    return qvals


def fetch_kegg_pathway_text(pathway_id: str, out_path: Path) -> None:
    if out_path.exists() and out_path.stat().st_size > 0:
        return
    url = f"https://rest.kegg.jp/get/{pathway_id}"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req) as r, open(out_path, "wb") as w:
        w.write(r.read())


def parse_kegg_gene_ids(path: Path) -> Set[str]:
    genes: Set[str] = set()
    in_gene = False
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if len(line) < 12:
                continue
            key = line[:12].strip()
            rest = line[12:].rstrip("\n")
            if key:
                in_gene = (key == "GENE")
                if not in_gene:
                    continue
            if not in_gene:
                continue
            # First token is Entrez Gene ID.
            token = rest.strip().split(" ", 1)[0] if rest.strip() else ""
            if re.fullmatch(r"\d+", token):
                genes.add(token)
    return genes


def build_arm_sizes(centromere_q_start: Dict[str, int], chrom_sizes: Dict[str, int]) -> Dict[str, int]:
    arm_sizes_bp: Dict[str, int] = {}
    for chrom, size in chrom_sizes.items():
        q_start = centromere_q_start.get(chrom)
        if q_start is None:
            continue
        label = chrom[3:]
        arm_sizes_bp[f"{label}p"] = q_start
        arm_sizes_bp[f"{label}q"] = size - q_start
    return arm_sizes_bp


def rate_ratio(k1: int, t1: int, k2: int, t2: int):
    if k1 == 0 or k2 == 0 or t1 <= 0 or t2 <= 0:
        return float("nan")
    return (k1 / t1) / (k2 / t2)


def main() -> int:
    DATA_DIR.mkdir(exist_ok=True)
    OUT_DIR.mkdir(exist_ok=True)

    files = base.url_file_map()
    for key, url in base.URLS.items():
        base.download_if_missing(url, files[key])

    # Load reference mapping.
    accn_to_chr = base.parse_assembly_refseq_to_chr(files["assembly_report"])
    protein_coding = base.load_protein_coding_genes(files["gene_info"])
    gene_locs = base.load_gene_locations(files["gene2refseq"], accn_to_chr)
    centromere_q_start = base.load_centromere_boundaries(files["cytoband"])
    chrom_sizes = base.load_chrom_sizes(files["chrom_sizes"])
    arm_sizes_bp = build_arm_sizes(centromere_q_start, chrom_sizes)

    eligible = set(protein_coding).intersection(gene_locs)
    gene_to_arm: Dict[str, str] = {}
    all_counts = Counter()
    for gid in eligible:
        arm = base.assign_arm(gene_locs[gid], centromere_q_start)
        if not arm or arm not in arm_sizes_bp:
            continue
        gene_to_arm[gid] = arm
        all_counts[arm] += 1

    arms = sorted(arm_sizes_bp, key=arm_sort_key)
    total_arm_bp = sum(arm_sizes_bp.values())
    one_q_bp = arm_sizes_bp["1q"]

    arm_rows = []
    summary_rows = []

    pvals = []
    idx_map = []

    for idx, (pid, pname) in enumerate(PATHWAYS):
        pfile = DATA_DIR / f"{pid}.txt"
        fetch_kegg_pathway_text(pid, pfile)
        genes = parse_kegg_gene_ids(pfile)
        mapped = [g for g in genes if g in gene_to_arm]

        counts = Counter(gene_to_arm[g] for g in mapped)
        densities = {arm: (counts.get(arm, 0) / arm_sizes_bp[arm]) * 1_000_000 for arm in arms}

        vals = [densities[a] for a in arms]
        mean = sum(vals) / len(vals)
        var = sum((x - mean) ** 2 for x in vals) / len(vals)
        sd = math.sqrt(var)
        z_1q = (densities["1q"] - mean) / sd if sd > 0 else 0.0

        rank = 1 + sum(1 for a in arms if densities[a] > densities["1q"])

        k1 = counts.get("1q", 0)
        n = sum(counts.values())
        p0 = one_q_bp / total_arm_bp
        p_one_sided = base.binom_tail_geq(k1, n, p0) if n > 0 else 1.0

        expected_1q = n * p0
        observed_over_expected = (k1 / expected_1q) if expected_1q > 0 else float("nan")
        rr = rate_ratio(k1, one_q_bp, n - k1, total_arm_bp - one_q_bp)

        summary_rows.append(
            {
                "pathway_id": pid,
                "pathway_name": pname,
                "pathway_gene_total": len(genes),
                "pathway_gene_mapped": n,
                "one_q_gene_count": k1,
                "one_q_expected_count_by_size": expected_1q,
                "one_q_obs_over_exp": observed_over_expected,
                "one_q_genes_per_mb": densities["1q"],
                "one_q_density_zscore": z_1q,
                "one_q_rank": rank,
                "arms_total": len(arms),
                "one_sided_p_higher": p_one_sided,
                "rate_ratio_1q_vs_rest": rr,
            }
        )
        pvals.append(p_one_sided)
        idx_map.append(idx)

        for arm in arms:
            arm_rows.append(
                {
                    "pathway_id": pid,
                    "pathway_name": pname,
                    "arm": arm,
                    "arm_size_mb": arm_sizes_bp[arm] / 1_000_000,
                    "pathway_gene_count": counts.get(arm, 0),
                    "pathway_genes_per_mb": densities[arm],
                    "one_q_genes_per_mb": densities["1q"],
                }
            )

    qvals = benjamini_hochberg(pvals)
    for i, q in enumerate(qvals):
        summary_rows[i]["fdr_bh"] = q

    summary_rows.sort(key=lambda r: (r["fdr_bh"], -r["one_q_obs_over_exp"]))

    arm_out = OUT_DIR / "oncogenic_pathway_arm_density.tsv"
    with open(arm_out, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "pathway_id",
                "pathway_name",
                "arm",
                "arm_size_mb",
                "pathway_gene_count",
                "pathway_genes_per_mb",
                "one_q_genes_per_mb",
            ],
            delimiter="\t",
        )
        w.writeheader()
        for row in arm_rows:
            w.writerow({
                k: (f"{v:.6f}" if isinstance(v, float) and not math.isnan(v) else v)
                for k, v in row.items()
            })

    sum_out = OUT_DIR / "oncogenic_pathway_1q_summary.tsv"
    with open(sum_out, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "pathway_id",
                "pathway_name",
                "pathway_gene_total",
                "pathway_gene_mapped",
                "one_q_gene_count",
                "one_q_expected_count_by_size",
                "one_q_obs_over_exp",
                "one_q_genes_per_mb",
                "one_q_density_zscore",
                "one_q_rank",
                "arms_total",
                "one_sided_p_higher",
                "fdr_bh",
                "rate_ratio_1q_vs_rest",
            ],
            delimiter="\t",
        )
        w.writeheader()
        for row in summary_rows:
            w.writerow({
                k: (f"{v:.6g}" if isinstance(v, float) and not math.isnan(v) else ("nan" if isinstance(v, float) else v))
                for k, v in row.items()
            })

    print(f"Wrote {arm_out}")
    print(f"Wrote {sum_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
