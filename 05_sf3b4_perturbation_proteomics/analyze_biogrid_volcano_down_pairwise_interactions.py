#!/usr/bin/env python3
"""Test BioGRID interactions among highlighted downregulated volcano hits.

The volcano labels are the manuscript-facing labels, but BioGRID may store the
same protein under a different official symbol. Matching therefore uses the
volcano gene/feature labels plus UniProt accessions from the responsive-hit
table.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


WORK_DIR = Path(__file__).resolve().parent
DEFAULT_RESULT_DIR = (
    WORK_DIR
    / "dep_shctrl_adjusted_no_rdes_strong_kd_plus_mhhes1_cellline_weighted_78_decontam_results"
)
DEFAULT_BIOGRID_FILE = (
    WORK_DIR
    / "resources"
    / "biogrid"
    / "5.0.258"
    / "BIOGRID-ORGANISM-Homo_sapiens-5.0.258.tab3.txt"
)
BIOGRID_RELEASE = "5.0.258"
BIOGRID_ARCHIVE_URL = (
    "https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/"
    "BIOGRID-5.0.258/BIOGRID-ORGANISM-5.0.258.tab3.zip"
)
BIOGRID_DOWNLOAD_PAGE = (
    "https://downloads.thebiogrid.org/BioGRID/Release-Archive/BIOGRID-5.0.258/"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Find BioGRID pairwise interactions among highlighted downregulated "
            "publication volcano hits."
        )
    )
    parser.add_argument("--result-dir", type=Path, default=DEFAULT_RESULT_DIR)
    parser.add_argument("--biogrid-file", type=Path, default=DEFAULT_BIOGRID_FILE)
    parser.add_argument("--out-dir", type=Path, default=None)
    return parser.parse_args()


def split_terms(value: str | None) -> list[str]:
    if value is None:
        return []
    value = value.strip()
    if not value or value in {"-", "NA", "NaN", "nan"}:
        return []
    parts = re.split(r"[|;,]", value)
    return [part.strip() for part in parts if part.strip() and part.strip() != "-"]


def add_accession_terms(terms: set[str], value: str | None) -> None:
    for token in split_terms(value):
        terms.add(token)
        if "-" in token:
            terms.add(token.split("-", 1)[0])


def read_csv_by_key(path: Path, key: str) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rows[row[key]] = row
    return rows


def build_hit_records(result_dir: Path) -> list[dict[str, object]]:
    publication_dir = result_dir / "tables" / "publication"
    volcano_file = publication_dir / "volcano_single_publication_hit_plot_data.csv"
    responsive_file = publication_dir / "sf3b4_shctrl_adjusted_publication_responsive_hits.csv"

    if not volcano_file.exists():
        raise FileNotFoundError(volcano_file)
    if not responsive_file.exists():
        raise FileNotFoundError(responsive_file)

    responsive_by_feature = read_csv_by_key(responsive_file, "feature_id")
    responsive_by_gene = read_csv_by_key(responsive_file, "gene_symbol")

    records: list[dict[str, object]] = []
    with volcano_file.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            direction = row.get("response_direction", "").lower()
            log_fc = float(row.get("logFC") or row.get("primary_logFC") or row.get("imputed_logFC"))
            if direction != "down" or log_fc >= 0:
                continue

            label = row["gene_symbol"]
            feature_id = row.get("feature_id", "")
            responsive = responsive_by_feature.get(feature_id) or responsive_by_gene.get(label) or {}

            terms = {label, feature_id}
            for col in ("Genes", "gene_symbol", "feature_id"):
                for token in split_terms(responsive.get(col)):
                    terms.add(token)
            add_accession_terms(terms, responsive.get("Protein.Group"))

            records.append(
                {
                    "volcano_symbol": label,
                    "feature_id": feature_id,
                    "terms": {term for term in terms if term},
                    "response_direction": row.get("response_direction", ""),
                    "hit_class": row.get("hit_class", ""),
                    "plot_group": row.get("plot_group", ""),
                    "logFC": row.get("logFC", ""),
                    "adj.P.Val": row.get("adj.P.Val", ""),
                    "Protein.Group": responsive.get("Protein.Group", ""),
                    "Genes": responsive.get("Genes", ""),
                    "Protein.Names": responsive.get("Protein.Names", ""),
                    "First.Protein.Description": responsive.get("First.Protein.Description", ""),
                    "confidence_class": responsive.get("confidence_class", row.get("hit_class", "")),
                }
            )

    return records


def interactor_terms(row: dict[str, str], suffix: str) -> set[str]:
    terms: set[str] = set()
    for col in (
        f"Official Symbol Interactor {suffix}",
        f"Systematic Name Interactor {suffix}",
    ):
        for token in split_terms(row.get(col)):
            terms.add(token)
    for col in (
        f"SWISS-PROT Accessions Interactor {suffix}",
        f"TREMBL Accessions Interactor {suffix}",
        f"REFSEQ Accessions Interactor {suffix}",
    ):
        add_accession_terms(terms, row.get(col))
    return terms


def normalise_header(fieldnames: list[str] | None) -> list[str]:
    if fieldnames is None:
        raise ValueError("BioGRID file has no header")
    return [re.sub(r"^#", "", field) for field in fieldnames]


def md5sum(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def collapse(values: list[str]) -> str:
    clean = sorted({value for value in values if value and value != "-"})
    return ";".join(clean)


def main() -> None:
    args = parse_args()
    result_dir = args.result_dir.resolve()
    biogrid_file = args.biogrid_file.resolve()
    out_dir = (args.out_dir or (result_dir / "tables" / "biogrid")).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    if not biogrid_file.exists():
        raise FileNotFoundError(biogrid_file)

    hit_records = build_hit_records(result_dir)
    hit_order = {record["volcano_symbol"]: index for index, record in enumerate(hit_records)}
    hit_labels = [record["volcano_symbol"] for record in hit_records]

    term_to_hits: dict[str, set[str]] = defaultdict(set)
    hit_terms: dict[str, set[str]] = {}
    for record in hit_records:
        label = str(record["volcano_symbol"])
        terms = set(record["terms"])
        hit_terms[label] = terms
        for term in terms:
            term_to_hits[term].add(label)

    def matching_hits(terms: set[str]) -> set[str]:
        matches: set[str] = set()
        for term in terms:
            matches.update(term_to_hits.get(term, set()))
        return matches

    mapping_official: dict[str, set[str]] = defaultdict(set)
    mapping_terms: dict[str, set[str]] = defaultdict(set)
    all_biogrid_row_counts: dict[str, int] = defaultdict(int)
    evidence_rows: list[dict[str, str]] = []

    required_cols = {
        "BioGRID Interaction ID",
        "Official Symbol Interactor A",
        "Official Symbol Interactor B",
        "Experimental System",
        "Experimental System Type",
        "Author",
        "Publication Source",
        "Organism ID Interactor A",
        "Organism ID Interactor B",
        "Throughput",
        "Score",
        "Modification",
        "Qualifications",
        "Tags",
        "Source Database",
        "SWISS-PROT Accessions Interactor A",
        "SWISS-PROT Accessions Interactor B",
        "TREMBL Accessions Interactor A",
        "TREMBL Accessions Interactor B",
        "REFSEQ Accessions Interactor A",
        "REFSEQ Accessions Interactor B",
        "Organism Name Interactor A",
        "Organism Name Interactor B",
    }

    with biogrid_file.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        reader.fieldnames = normalise_header(reader.fieldnames)
        missing = required_cols.difference(reader.fieldnames)
        if missing:
            raise ValueError(f"BioGRID file is missing required columns: {sorted(missing)}")

        for row in reader:
            terms_a = interactor_terms(row, "A")
            terms_b = interactor_terms(row, "B")
            hits_a = matching_hits(terms_a)
            hits_b = matching_hits(terms_b)

            for label in hits_a:
                mapping_official[label].add(row["Official Symbol Interactor A"])
                mapping_terms[label].update(hit_terms[label].intersection(terms_a))
                all_biogrid_row_counts[label] += 1
            for label in hits_b:
                mapping_official[label].add(row["Official Symbol Interactor B"])
                mapping_terms[label].update(hit_terms[label].intersection(terms_b))
                all_biogrid_row_counts[label] += 1

            if not hits_a or not hits_b:
                continue
            if row["Organism ID Interactor A"] != "9606" or row["Organism ID Interactor B"] != "9606":
                continue

            for hit_a, hit_b in itertools.product(hits_a, hits_b):
                if hit_a == hit_b:
                    continue
                pair = tuple(sorted((hit_a, hit_b), key=lambda label: hit_order[label]))
                evidence_rows.append(
                    {
                        "query_symbol_1": pair[0],
                        "query_symbol_2": pair[1],
                        "query_pair": f"{pair[0]}--{pair[1]}",
                        "query_symbol_interactor_a": hit_a,
                        "query_symbol_interactor_b": hit_b,
                        "biogrid_interaction_id": row["BioGRID Interaction ID"],
                        "official_symbol_interactor_a": row["Official Symbol Interactor A"],
                        "official_symbol_interactor_b": row["Official Symbol Interactor B"],
                        "experimental_system": row["Experimental System"],
                        "experimental_system_type": row["Experimental System Type"],
                        "author": row["Author"],
                        "publication_source": row["Publication Source"],
                        "throughput": row["Throughput"],
                        "score": row["Score"],
                        "modification": row["Modification"],
                        "qualifications": row["Qualifications"],
                        "tags": row["Tags"],
                        "source_database": row["Source Database"],
                        "organism_name_interactor_a": row["Organism Name Interactor A"],
                        "organism_name_interactor_b": row["Organism Name Interactor B"],
                    }
                )

    evidence_rows.sort(
        key=lambda row: (
            hit_order[row["query_symbol_1"]],
            hit_order[row["query_symbol_2"]],
            row["experimental_system_type"],
            row["experimental_system"],
            row["publication_source"],
            row["biogrid_interaction_id"],
        )
    )

    summary_by_pair: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in evidence_rows:
        summary_by_pair[(row["query_symbol_1"], row["query_symbol_2"])].append(row)

    summary_rows: list[dict[str, str | int]] = []
    for pair, rows in summary_by_pair.items():
        summary_rows.append(
            {
                "query_symbol_1": pair[0],
                "query_symbol_2": pair[1],
                "query_pair": f"{pair[0]}--{pair[1]}",
                "n_biogrid_evidence_rows": len(rows),
                "n_physical_evidence_rows": sum(
                    row["experimental_system_type"] == "physical" for row in rows
                ),
                "n_genetic_evidence_rows": sum(
                    row["experimental_system_type"] == "genetic" for row in rows
                ),
                "experimental_system_types": collapse(
                    [row["experimental_system_type"] for row in rows]
                ),
                "experimental_systems": collapse([row["experimental_system"] for row in rows]),
                "throughput": collapse([row["throughput"] for row in rows]),
                "publication_sources": collapse([row["publication_source"] for row in rows]),
                "authors": collapse([row["author"] for row in rows]),
                "biogrid_interaction_ids": collapse(
                    [row["biogrid_interaction_id"] for row in rows]
                ),
            }
        )

    summary_rows.sort(
        key=lambda row: (
            hit_order[str(row["query_symbol_1"])],
            hit_order[str(row["query_symbol_2"])],
        )
    )

    matrix_rows: list[dict[str, str | int | bool]] = []
    for first, second in itertools.combinations(hit_labels, 2):
        rows = summary_by_pair.get((first, second), [])
        matrix_rows.append(
            {
                "query_symbol_1": first,
                "query_symbol_2": second,
                "query_pair": f"{first}--{second}",
                "has_biogrid_interaction": bool(rows),
                "n_biogrid_evidence_rows": len(rows),
                "n_physical_evidence_rows": sum(
                    row["experimental_system_type"] == "physical" for row in rows
                ),
                "n_genetic_evidence_rows": sum(
                    row["experimental_system_type"] == "genetic" for row in rows
                ),
            }
        )

    mapping_rows: list[dict[str, str | int]] = []
    for record in hit_records:
        label = str(record["volcano_symbol"])
        official_symbols = sorted(mapping_official.get(label, set()))
        mapping_rows.append(
            {
                "volcano_symbol": label,
                "feature_id": str(record["feature_id"]),
                "response_direction": str(record["response_direction"]),
                "hit_class": str(record["hit_class"]),
                "confidence_class": str(record["confidence_class"]),
                "logFC": str(record["logFC"]),
                "adj.P.Val": str(record["adj.P.Val"]),
                "protein_group": str(record["Protein.Group"]),
                "genes": str(record["Genes"]),
                "protein_names": str(record["Protein.Names"]),
                "first_protein_description": str(record["First.Protein.Description"]),
                "query_terms_used": collapse(sorted(hit_terms[label])),
                "biogrid_official_symbols": collapse(official_symbols),
                "matched_terms_in_biogrid": collapse(sorted(mapping_terms.get(label, set()))),
                "n_biogrid_rows_with_this_protein": all_biogrid_row_counts.get(label, 0),
                "mapping_status": "mapped" if official_symbols else "unmapped",
            }
        )

    write_csv(out_dir / "volcano_down_highlighted_biogrid_input_mapping.csv", mapping_rows)
    write_csv(out_dir / "volcano_down_highlighted_biogrid_pairwise_evidence.csv", evidence_rows)
    write_csv(out_dir / "volcano_down_highlighted_biogrid_pairwise_summary.csv", summary_rows)
    write_csv(out_dir / "volcano_down_highlighted_biogrid_pairwise_matrix.csv", matrix_rows)
    write_readme(
        out_dir=out_dir,
        result_dir=result_dir,
        biogrid_file=biogrid_file,
        hit_labels=hit_labels,
        mapping_rows=mapping_rows,
        summary_rows=summary_rows,
        matrix_rows=matrix_rows,
    )

    print(f"Wrote BioGRID pairwise outputs to: {out_dir}")
    print(
        f"Found {len(summary_rows)} interacting protein pairs among "
        f"{len(hit_labels)} highlighted downregulated volcano proteins."
    )


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_readme(
    out_dir: Path,
    result_dir: Path,
    biogrid_file: Path,
    hit_labels: list[str],
    mapping_rows: list[dict[str, str | int]],
    summary_rows: list[dict[str, str | int]],
    matrix_rows: list[dict[str, str | int | bool]],
) -> None:
    interacting_pairs = [
        f"{row['query_pair']} ({row['n_biogrid_evidence_rows']} evidence rows)"
        for row in summary_rows
    ]
    unmapped = [str(row["volcano_symbol"]) for row in mapping_rows if row["mapping_status"] != "mapped"]
    mapped_aliases = [
        f"{row['volcano_symbol']} -> {row['biogrid_official_symbols']}"
        for row in mapping_rows
        if row["volcano_symbol"] != row["biogrid_official_symbols"]
    ]
    n_possible = len(matrix_rows)
    n_observed = len(summary_rows)

    lines = [
        "# BioGRID pairwise interactions among highlighted downregulated volcano proteins",
        "",
        f"Generated: {datetime.now(timezone.utc).isoformat()}",
        "",
        "Inputs:",
        "",
        f"- Volcano hits: `{result_dir / 'tables' / 'publication' / 'volcano_single_publication_hit_plot_data.csv'}`",
        f"- Responsive-hit metadata: `{result_dir / 'tables' / 'publication' / 'sf3b4_shctrl_adjusted_publication_responsive_hits.csv'}`",
        f"- BioGRID Homo sapiens Tab 3 file: `{biogrid_file}`",
        f"- BioGRID release: `{BIOGRID_RELEASE}`",
        f"- BioGRID archive URL: {BIOGRID_ARCHIVE_URL}",
        f"- BioGRID download page: {BIOGRID_DOWNLOAD_PAGE}",
        f"- BioGRID file MD5: `{md5sum(biogrid_file)}`",
        "",
        "Protein set:",
        "",
        f"- {len(hit_labels)} highlighted downregulated volcano proteins: `{', '.join(hit_labels)}`",
        "- Self-interactions were excluded from the pairwise interaction summary.",
        "- Only BioGRID rows with both interactors annotated as Homo sapiens (`9606`) were used for pairwise interactions.",
        "",
        "Symbol/accession mapping:",
        "",
        *(f"- {item}" for item in mapped_aliases),
        *(["- Unmapped proteins: `" + "`, `".join(unmapped) + "`"] if unmapped else []),
        "",
        "Pairwise result:",
        "",
        f"- Observed BioGRID interactions in {n_observed} of {n_possible} possible non-self pairs.",
        *(f"- {pair}" for pair in interacting_pairs),
        "",
        "Outputs:",
        "",
        "- `volcano_down_highlighted_biogrid_input_mapping.csv`: volcano labels, accession/symbol terms used, and BioGRID mapping status.",
        "- `volcano_down_highlighted_biogrid_pairwise_evidence.csv`: BioGRID evidence rows for observed interactions between highlighted downregulated proteins.",
        "- `volcano_down_highlighted_biogrid_pairwise_summary.csv`: one row per observed interacting pair.",
        "- `volcano_down_highlighted_biogrid_pairwise_matrix.csv`: all possible non-self pairs with interaction flags and evidence counts.",
    ]
    (out_dir / "README_volcano_down_highlighted_pairwise.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
