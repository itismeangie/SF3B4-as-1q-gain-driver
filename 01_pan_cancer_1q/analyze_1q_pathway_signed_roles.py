#!/usr/bin/env python3
"""Pathway-signed role analysis for 1q using KEGG KGML edges.

For each selected pathway, classify genes by pathway-specific signed outgoing edges:
- activator_like: outgoing activation/expression > outgoing inhibition/repression
- repressor_like: outgoing inhibition/repression > outgoing activation/expression
- dual_like: both signs and tied

Then test whether activator-like genes are more represented on 1q than repressor-like genes.
"""

from __future__ import annotations

import csv
import math
import re
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Set, Tuple

import analyze_1q_developmental_enrichment as base
from analyze_1q_oncogenic_pathways import PATHWAYS

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
OUT_DIR = BASE_DIR / "results"

POS_SIGNS = {"activation", "expression"}
NEG_SIGNS = {"inhibition", "repression"}


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


def log_comb(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def fisher_exact_one_sided_greater(a: int, b: int, c: int, d: int) -> float:
    # Table:
    #              1q   non-1q
    # activator     a      b
    # repressor     c      d
    row1 = a + b
    row2 = c + d
    col1 = a + c
    n = row1 + row2
    if row1 == 0 or row2 == 0:
        return 1.0

    lo = max(0, col1 - row2)
    hi = min(row1, col1)

    def log_hypergeom(x: int) -> float:
        return log_comb(row1, x) + log_comb(row2, col1 - x) - log_comb(n, col1)

    logs = [log_hypergeom(x) for x in range(a, hi + 1)]
    m = max(logs)
    return math.exp(m) * sum(math.exp(x - m) for x in logs)


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


def fetch_kgml(pathway_id: str, out_path: Path) -> None:
    if out_path.exists() and out_path.stat().st_size > 0:
        return
    url = f"https://rest.kegg.jp/get/{pathway_id}/kgml"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req) as r, open(out_path, "wb") as w:
        w.write(r.read())


def parse_hsa_gene_ids(text: str) -> Set[str]:
    return set(re.findall(r"hsa:(\d+)", text or ""))


def parse_kegg_signed_roles(kgml_path: Path) -> Tuple[Set[str], Set[str], Set[str], Set[str]]:
    tree = ET.parse(kgml_path)
    root = tree.getroot()

    entry_genes_raw: Dict[str, Set[str]] = {}
    group_components: Dict[str, List[str]] = {}

    for entry in root.findall("entry"):
        eid = entry.attrib.get("id")
        if not eid:
            continue
        etype = entry.attrib.get("type", "")
        name = entry.attrib.get("name", "")
        genes = parse_hsa_gene_ids(name)
        entry_genes_raw[eid] = set(genes)
        if etype == "group":
            comps = [c.attrib.get("id") for c in entry.findall("component") if c.attrib.get("id")]
            group_components[eid] = comps

    memo: Dict[str, Set[str]] = {}

    def resolve_entry_genes(eid: str, seen: Set[str] | None = None) -> Set[str]:
        if eid in memo:
            return memo[eid]
        if seen is None:
            seen = set()
        if eid in seen:
            return set()
        seen.add(eid)

        out = set(entry_genes_raw.get(eid, set()))
        for comp in group_components.get(eid, []):
            out |= resolve_entry_genes(comp, seen)

        memo[eid] = out
        return out

    pos_out: Dict[str, int] = {}
    neg_out: Dict[str, int] = {}
    all_genes_in_pathway: Set[str] = set()

    for rel in root.findall("relation"):
        e1 = rel.attrib.get("entry1")
        e2 = rel.attrib.get("entry2")
        if not e1 or not e2:
            continue

        subtypes = [s.attrib.get("name", "") for s in rel.findall("subtype")]
        has_pos = any(s in POS_SIGNS for s in subtypes)
        has_neg = any(s in NEG_SIGNS for s in subtypes)

        if has_pos and has_neg:
            continue
        if not has_pos and not has_neg:
            continue

        source_genes = resolve_entry_genes(e1)
        target_genes = resolve_entry_genes(e2)

        if not source_genes or not target_genes:
            continue

        all_genes_in_pathway |= source_genes
        all_genes_in_pathway |= target_genes

        for g1 in source_genes:
            if has_pos:
                pos_out[g1] = pos_out.get(g1, 0) + 1
            elif has_neg:
                neg_out[g1] = neg_out.get(g1, 0) + 1

    activator_like: Set[str] = set()
    repressor_like: Set[str] = set()
    dual_like: Set[str] = set()

    for g in all_genes_in_pathway:
        p = pos_out.get(g, 0)
        n = neg_out.get(g, 0)
        if p == 0 and n == 0:
            continue
        if p > n:
            activator_like.add(g)
        elif n > p:
            repressor_like.add(g)
        else:
            dual_like.add(g)

    return all_genes_in_pathway, activator_like, repressor_like, dual_like


def main() -> int:
    DATA_DIR.mkdir(exist_ok=True)
    OUT_DIR.mkdir(exist_ok=True)

    files = base.url_file_map()
    for key, url in base.URLS.items():
        base.download_if_missing(url, files[key])

    accn_to_chr = base.parse_assembly_refseq_to_chr(files["assembly_report"])
    protein_coding = base.load_protein_coding_genes(files["gene_info"])
    gene_locs = base.load_gene_locations(files["gene2refseq"], accn_to_chr)
    centromere_q_start = base.load_centromere_boundaries(files["cytoband"])
    chrom_sizes = base.load_chrom_sizes(files["chrom_sizes"])
    arm_sizes_bp = build_arm_sizes(centromere_q_start, chrom_sizes)

    eligible = set(protein_coding).intersection(gene_locs)
    gene_to_arm: Dict[str, str] = {}
    for gid in eligible:
        arm = base.assign_arm(gene_locs[gid], centromere_q_start)
        if arm and arm in arm_sizes_bp:
            gene_to_arm[gid] = arm

    total_bp = sum(arm_sizes_bp.values())
    p0 = arm_sizes_bp["1q"] / total_bp

    rows: List[Dict[str, object]] = []
    pvals: List[float] = []

    for pid, pname in PATHWAYS:
        kgml_path = DATA_DIR / f"{pid}.kgml"
        fetch_kgml(pid, kgml_path)

        all_genes, act_like, rep_like, dual_like = parse_kegg_signed_roles(kgml_path)

        mapped_all = {g for g in all_genes if g in gene_to_arm}
        mapped_act = {g for g in act_like if g in gene_to_arm}
        mapped_rep = {g for g in rep_like if g in gene_to_arm}
        mapped_dual = {g for g in dual_like if g in gene_to_arm}

        n_act = len(mapped_act)
        n_rep = len(mapped_rep)

        k_act_1q = sum(1 for g in mapped_act if gene_to_arm[g] == "1q")
        k_rep_1q = sum(1 for g in mapped_rep if gene_to_arm[g] == "1q")

        exp_act = n_act * p0
        exp_rep = n_rep * p0

        oe_act = (k_act_1q / exp_act) if exp_act > 0 else float("nan")
        oe_rep = (k_rep_1q / exp_rep) if exp_rep > 0 else float("nan")

        oe_act_adj = ((k_act_1q + 0.5) / (exp_act + 0.5)) if exp_act > 0 else float("nan")
        oe_rep_adj = ((k_rep_1q + 0.5) / (exp_rep + 0.5)) if exp_rep > 0 else float("nan")

        if math.isnan(oe_act_adj) or math.isnan(oe_rep_adj):
            net_log2 = float("nan")
        else:
            net_log2 = math.log2(oe_act_adj / oe_rep_adj)

        if n_act > 0 and n_rep > 0:
            p_dir = fisher_exact_one_sided_greater(
                k_act_1q,
                n_act - k_act_1q,
                k_rep_1q,
                n_rep - k_rep_1q,
            )
        else:
            p_dir = 1.0

        rows.append(
            {
                "pathway_id": pid,
                "pathway_name": pname,
                "pathway_genes_with_signed_edges": len(mapped_all),
                "activator_like_genes": n_act,
                "repressor_like_genes": n_rep,
                "dual_like_genes": len(mapped_dual),
                "activator_1q_count": k_act_1q,
                "repressor_1q_count": k_rep_1q,
                "activator_oe_1q": oe_act,
                "repressor_oe_1q": oe_rep,
                "activator_oe_1q_adj": oe_act_adj,
                "repressor_oe_1q_adj": oe_rep_adj,
                "net_log2_activator_over_repressor": net_log2,
                "directional_p_one_sided": p_dir,
            }
        )
        pvals.append(p_dir)

    qvals = benjamini_hochberg(pvals)
    for i, q in enumerate(qvals):
        rows[i]["directional_fdr_bh"] = q

    rows.sort(
        key=lambda r: (
            r["directional_fdr_bh"],
            -(r["net_log2_activator_over_repressor"] if not math.isnan(r["net_log2_activator_over_repressor"]) else -999),
        )
    )

    out = OUT_DIR / "oncogenic_pathway_signed_role_1q_summary.tsv"
    with open(out, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "pathway_id",
                "pathway_name",
                "pathway_genes_with_signed_edges",
                "activator_like_genes",
                "repressor_like_genes",
                "dual_like_genes",
                "activator_1q_count",
                "repressor_1q_count",
                "activator_oe_1q",
                "repressor_oe_1q",
                "activator_oe_1q_adj",
                "repressor_oe_1q_adj",
                "net_log2_activator_over_repressor",
                "directional_p_one_sided",
                "directional_fdr_bh",
            ],
            delimiter="\t",
        )
        w.writeheader()
        for r in rows:
            w.writerow(
                {
                    k: (f"{v:.6g}" if isinstance(v, float) and not math.isnan(v) else ("nan" if isinstance(v, float) else v))
                    for k, v in r.items()
                }
            )

    print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
