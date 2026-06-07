#!/usr/bin/env python3
"""Assess whether 1q has a disproportionate density of developmental genes.

Data sources (downloaded automatically):
- NCBI gene_info (protein-coding gene catalog)
- NCBI gene2ensembl (gene genomic coordinates)
- NCBI gene2go + GO basic OBO (developmental process descendants)
- UCSC hg38 cytoband + chromosome sizes (arm boundaries and arm sizes)
"""

from __future__ import annotations

import csv
import gzip
import math
import os
import urllib.request
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
OUT_DIR = BASE_DIR / "results"
DEFAULT_REFERENCE_CACHE = Path.home() / ".cache" / "1q-gain-sf3b4-analysis" / "1q_genes"

URLS = {
    "gene_info": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz",
    "gene2refseq": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz",
    "gene2go": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz",
    "go_basic": "https://purl.obolibrary.org/obo/go/go-basic.obo",
    "assembly_report": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt",
    "cytoband": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz",
    "chrom_sizes": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes",
}

ALLOWED_CHR = {str(i) for i in range(1, 23)} | {"X", "Y"}
ROOT_GO = "GO:0032502"  # developmental process


def reference_data_dir() -> Path:
    """Cache large downloaded reference files outside the git repo by default."""
    return Path(os.environ.get("ONE_Q_REFERENCE_CACHE", DEFAULT_REFERENCE_CACHE)).expanduser()


def url_file_map(data_dir: Path | None = None) -> Dict[str, Path]:
    root = data_dir or reference_data_dir()
    return {key: root / Path(url).name for key, url in URLS.items()}


@dataclass
class GeneLoc:
    chrom: str
    start: int
    end: int

    @property
    def midpoint(self) -> int:
        return (self.start + self.end) // 2


def download_if_missing(url: str, path: Path) -> None:
    if path.exists() and path.stat().st_size > 0:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {url} -> {path}")
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req) as r, open(path, "wb") as w:
        w.write(r.read())


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, "rt", encoding="utf-8", newline="")


def parse_assembly_refseq_to_chr(path: Path) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    with open_text(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            sequence_role = parts[1]
            assigned_molecule = parts[2]
            refseq_accn = parts[6]
            if sequence_role != "assembled-molecule":
                continue
            if assigned_molecule in ALLOWED_CHR and refseq_accn != "na":
                mapping[refseq_accn] = assigned_molecule
    return mapping


def load_protein_coding_genes(path: Path) -> Dict[str, str]:
    genes: Dict[str, str] = {}
    with open_text(path) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            # tax_id, GeneID, Symbol, ..., type_of_gene
            if row[0] != "9606":
                continue
            gene_id = row[1]
            symbol = row[2]
            gene_type = row[9]
            if gene_type == "protein-coding":
                genes[gene_id] = symbol
    return genes


def load_gene_locations(path: Path, accn_to_chr: Dict[str, str]) -> Dict[str, GeneLoc]:
    # Aggregate potentially repeated rows per gene and chromosome.
    accum: Dict[str, Dict[str, List[Tuple[int, int]]]] = defaultdict(lambda: defaultdict(list))
    with open_text(path) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 13:
                continue
            if row[0] != "9606":
                continue
            gene_id = row[1]
            genomic_accn = row[7]
            if genomic_accn not in accn_to_chr:
                continue
            chrom = accn_to_chr[genomic_accn]
            if chrom not in ALLOWED_CHR:
                continue
            try:
                s = int(row[9])
                e = int(row[10])
            except ValueError:
                continue
            if s < 0 or e < 0:
                continue
            start, end = (s, e) if s <= e else (e, s)
            accum[gene_id][chrom].append((start, end))

    locs: Dict[str, GeneLoc] = {}
    for gene_id, by_chr in accum.items():
        # Pick chromosome with most supporting transcript rows.
        best_chr = max(by_chr, key=lambda c: len(by_chr[c]))
        starts = [s for s, _ in by_chr[best_chr]]
        ends = [e for _, e in by_chr[best_chr]]
        locs[gene_id] = GeneLoc(best_chr, min(starts), max(ends))
    return locs


def parse_go_graph(path: Path) -> Dict[str, Set[str]]:
    # child adjacency from is_a + part_of relationships
    parents: Dict[str, Set[str]] = defaultdict(set)
    current_id = None
    obsolete = False
    current_parents: Set[str] = set()

    def flush_term() -> None:
        nonlocal current_id, obsolete, current_parents
        if current_id and not obsolete:
            parents[current_id].update(current_parents)
        current_id = None
        obsolete = False
        current_parents = set()

    with open_text(path) as f:
        for raw in f:
            line = raw.strip()
            if line == "[Term]":
                flush_term()
                continue
            if not line:
                continue
            if line.startswith("id: GO:"):
                current_id = line.split("id: ", 1)[1]
            elif line.startswith("is_a: GO:"):
                parent = line.split("is_a: ", 1)[1].split(" ! ", 1)[0]
                current_parents.add(parent)
            elif line.startswith("relationship: part_of GO:"):
                parent = line.split("relationship: part_of ", 1)[1].split(" ! ", 1)[0]
                current_parents.add(parent)
            elif line == "is_obsolete: true":
                obsolete = True
        flush_term()

    children: Dict[str, Set[str]] = defaultdict(set)
    for child, pars in parents.items():
        for p in pars:
            children[p].add(child)
    return children


def descendants(children: Dict[str, Set[str]], root: str) -> Set[str]:
    seen = set()
    stack = [root]
    while stack:
        node = stack.pop()
        if node in seen:
            continue
        seen.add(node)
        stack.extend(children.get(node, ()))
    return seen


def load_developmental_gene_ids(path: Path, dev_go_terms: Set[str]) -> Set[str]:
    gene_ids: Set[str] = set()
    with open_text(path) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if row[0] != "9606":
                continue
            gene_id = row[1]
            go_id = row[2]
            qualifier = row[4]
            category = row[7] if len(row) > 7 else ""
            if category != "Process":
                continue
            if "NOT" in qualifier:
                continue
            if go_id in dev_go_terms:
                gene_ids.add(gene_id)
    return gene_ids


def load_centromere_boundaries(path: Path) -> Dict[str, int]:
    # boundary is q-arm start (first q acen start)
    q_starts: Dict[str, List[int]] = defaultdict(list)
    with open_text(path) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) < 5:
                continue
            chrom, start, _end, band, stain = row
            if stain != "acen":
                continue
            if not chrom.startswith("chr"):
                continue
            label = chrom[3:]
            if label not in ALLOWED_CHR:
                continue
            if band.startswith("q"):
                q_starts[chrom].append(int(start))

    out: Dict[str, int] = {}
    for chrom, starts in q_starts.items():
        if starts:
            out[chrom] = min(starts)
    return out


def load_chrom_sizes(path: Path) -> Dict[str, int]:
    sizes: Dict[str, int] = {}
    with open_text(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) != 2:
                continue
            chrom, size = parts
            if not chrom.startswith("chr"):
                continue
            label = chrom[3:]
            if label in ALLOWED_CHR:
                sizes[chrom] = int(size)
    return sizes


def assign_arm(gene: GeneLoc, centromere_q_start: Dict[str, int]) -> str | None:
    chr_ucsc = f"chr{gene.chrom}"
    q_start = centromere_q_start.get(chr_ucsc)
    if q_start is None:
        return None
    arm = "p" if gene.midpoint < q_start else "q"
    return f"{gene.chrom}{arm}"


def safe_rate(count: int, size_bp: int) -> float:
    return (count / size_bp) * 1_000_000 if size_bp > 0 else float("nan")


def log_binom_pmf(k: int, n: int, p: float) -> float:
    if k < 0 or k > n:
        return float("-inf")
    if p == 0:
        return 0.0 if k == 0 else float("-inf")
    if p == 1:
        return 0.0 if k == n else float("-inf")
    return (
        math.lgamma(n + 1)
        - math.lgamma(k + 1)
        - math.lgamma(n - k + 1)
        + k * math.log(p)
        + (n - k) * math.log(1 - p)
    )


def binom_tail_geq(k: int, n: int, p: float) -> float:
    logs = [log_binom_pmf(i, n, p) for i in range(k, n + 1)]
    m = max(logs)
    if m == float("-inf"):
        return 0.0
    return math.exp(m) * sum(math.exp(x - m) for x in logs)


def rate_ratio_ci(k1: int, t1: float, k2: int, t2: float, z: float = 1.96) -> Tuple[float, float, float]:
    r1 = k1 / t1
    r2 = k2 / t2
    rr = r1 / r2
    se = math.sqrt((1 / k1) + (1 / k2))
    lo = math.exp(math.log(rr) - z * se)
    hi = math.exp(math.log(rr) + z * se)
    return rr, lo, hi


def main() -> int:
    OUT_DIR.mkdir(exist_ok=True)

    files = url_file_map()
    for key, url in URLS.items():
        download_if_missing(url, files[key])

    accn_to_chr = parse_assembly_refseq_to_chr(files["assembly_report"])
    protein_coding = load_protein_coding_genes(files["gene_info"])
    gene_locs = load_gene_locations(files["gene2refseq"], accn_to_chr)

    children = parse_go_graph(files["go_basic"])
    dev_terms = descendants(children, ROOT_GO)
    dev_gene_ids = load_developmental_gene_ids(files["gene2go"], dev_terms)

    centromere_q_start = load_centromere_boundaries(files["cytoband"])
    chrom_sizes = load_chrom_sizes(files["chrom_sizes"])

    arm_sizes_bp: Dict[str, int] = {}
    for chrom, size in chrom_sizes.items():
        q_start = centromere_q_start.get(chrom)
        if q_start is None:
            continue
        label = chrom[3:]
        arm_sizes_bp[f"{label}p"] = q_start
        arm_sizes_bp[f"{label}q"] = size - q_start

    all_counts = Counter()
    dev_counts = Counter()

    eligible_gene_ids = set(protein_coding).intersection(gene_locs)
    for gene_id in eligible_gene_ids:
        arm = assign_arm(gene_locs[gene_id], centromere_q_start)
        if not arm or arm not in arm_sizes_bp:
            continue
        all_counts[arm] += 1
        if gene_id in dev_gene_ids:
            dev_counts[arm] += 1

    rows = []
    for arm, size_bp in arm_sizes_bp.items():
        rows.append(
            {
                "arm": arm,
                "size_mb": size_bp / 1_000_000,
                "protein_coding_genes": all_counts.get(arm, 0),
                "developmental_genes": dev_counts.get(arm, 0),
                "developmental_genes_per_mb": safe_rate(dev_counts.get(arm, 0), size_bp),
                "dev_fraction_of_genes": (
                    dev_counts.get(arm, 0) / all_counts.get(arm, 1)
                    if all_counts.get(arm, 0) > 0
                    else float("nan")
                ),
            }
        )

    rows.sort(key=lambda r: r["developmental_genes_per_mb"], reverse=True)

    out_table = OUT_DIR / "chrom_arm_developmental_density.tsv"
    with open(out_table, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(
            [
                "arm",
                "size_mb",
                "protein_coding_genes",
                "developmental_genes",
                "developmental_genes_per_mb",
                "dev_fraction_of_genes",
            ]
        )
        for r in rows:
            writer.writerow(
                [
                    r["arm"],
                    f"{r['size_mb']:.4f}",
                    r["protein_coding_genes"],
                    r["developmental_genes"],
                    f"{r['developmental_genes_per_mb']:.6f}",
                    f"{r['dev_fraction_of_genes']:.6f}" if not math.isnan(r["dev_fraction_of_genes"]) else "nan",
                ]
            )

    # 1q vs all other arms: one-sided test for higher developmental-gene rate per Mb.
    k1 = dev_counts.get("1q", 0)
    t1 = arm_sizes_bp.get("1q", 0)
    k_total = sum(dev_counts.values())
    t_total = sum(arm_sizes_bp.values())
    k2 = k_total - k1
    t2 = t_total - t1

    if k1 == 0 or k2 == 0 or t1 == 0 or t2 == 0:
        raise RuntimeError("Insufficient data for rate-ratio test.")

    p0 = t1 / (t1 + t2)
    pval = binom_tail_geq(k1, k1 + k2, p0)
    rr, rr_lo, rr_hi = rate_ratio_ci(k1, t1, k2, t2)

    rank_1q = next((i + 1 for i, r in enumerate(rows) if r["arm"] == "1q"), None)
    one_q_density = next(r["developmental_genes_per_mb"] for r in rows if r["arm"] == "1q")
    mean_density = sum(r["developmental_genes_per_mb"] for r in rows) / len(rows)

    summary_path = OUT_DIR / "summary_1q_vs_other_arms.txt"
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("1q developmental-gene density analysis (hg38, protein-coding genes)\n")
        f.write(f"Eligible protein-coding genes with mapped location: {sum(all_counts.values())}\n")
        f.write(f"Eligible developmental genes: {k_total}\n")
        f.write(f"1q developmental genes: {k1}\n")
        f.write(f"1q size (Mb): {t1 / 1_000_000:.4f}\n")
        f.write(f"1q developmental genes per Mb: {one_q_density:.6f}\n")
        f.write(f"Genome-wide mean arm developmental genes per Mb: {mean_density:.6f}\n")
        f.write(f"1q rank by developmental genes per Mb (1=highest): {rank_1q} of {len(rows)}\n")
        f.write(f"Rate ratio (1q vs all other arms): {rr:.4f}\n")
        f.write(f"95% CI for rate ratio: [{rr_lo:.4f}, {rr_hi:.4f}]\n")
        f.write(f"One-sided binomial exact p-value (H1: 1q higher rate): {pval:.6g}\n")

    print(f"Wrote {out_table}")
    print(f"Wrote {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
