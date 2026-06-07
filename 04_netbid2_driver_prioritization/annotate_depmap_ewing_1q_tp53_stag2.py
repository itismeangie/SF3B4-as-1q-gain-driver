#!/usr/bin/env python3
"""Annotate heatmap Ewing sarcoma cell lines for 1q gain, TP53, and STAG2.

Inputs are DepMap Public 26Q1 files cached locally. The script downloads
missing DepMap/CN/mutation files via the DepMap metadata cache when needed.
"""

from __future__ import annotations

import gzip
import json
import os
import pathlib
import re
import subprocess
import urllib.request

import numpy as np
import pandas as pd


TOPIC_DIR = pathlib.Path(__file__).resolve().parent
DEP_MAP_DIR = pathlib.Path(
    os.environ.get("DEPMAP_26Q1_DIR", TOPIC_DIR / "depmap" / "DepMap_Public_26Q1")
)
ANNOTATION_DIR = pathlib.Path(
    os.environ.get("DEPMAP_GENE_ANNOTATION_DIR", TOPIC_DIR / "depmap" / "gene_annotation")
)
TABLE_DIR = pathlib.Path(
    os.environ.get("DEPMAP_TABLE_DIR", TOPIC_DIR / "tables")
)
HEATMAP_DEPMAP_TABLE = pathlib.Path(
    os.environ.get(
        "HEATMAP_DEPMAP_TABLE",
        TABLE_DIR / "top25_driver_multimodal_heatmap_depmap_ewing_gene_effect.tsv",
    )
)

DEPMAP_BASE_URL = "https://depmap.org"
CYTOBAND_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz"
NCBI_GENE_INFO_URL = (
    "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
)

ONEQ_GAIN_GENE_CN_MEAN_THRESHOLD = 1.10

# Literature/resource curation after cross-checking DepMap 26Q1 default calls
# against Cellosaurus, Sanger Cell Model Passports, and Ewing sarcoma genomic
# cell-line tables. These columns are used for figure annotation, while the
# Broad/DepMap-default columns are retained unchanged for auditability.
CURATED_MUTATION_OVERRIDES = {
    "MHH-ES-1": {
        "STAG2": "p.Gln735fs [Cellosaurus CVCL_1411; Tirode et al., Cancer Discovery 2014]"
    },
    "SK-N-MC": {
        "TP53": "p.M1_T125Del/c.170_572del [Cellosaurus CVCL_0530; Tirode et al., Cancer Discovery 2014]",
        "STAG2": "p.M1_R546Del [Tirode et al., Cancer Discovery 2014]",
    },
}

CURATED_ONEQ_STATUS_OVERRIDES = {
    "ES4": ("Gain", "Sanger Cell Model Passports 1q gene-level CNV/GISTIC"),
    "ES5": ("Gain", "Sanger Cell Model Passports 1q gene-level CNV/GISTIC"),
    "ES8": ("Gain", "Sanger Cell Model Passports 1q gene-level CNV/GISTIC"),
    "EW-22": ("Gain", "Sanger Cell Model Passports 1q gene-level CNV/GISTIC"),
    "EW-7": ("Gain", "Sanger Cell Model Passports 1q gene-level CNV/GISTIC"),
    "EW-1": ("No gain", "Sanger Cell Model Passports 1q gene-level CNV/GISTIC"),
    "EW-16": ("No gain", "Sanger Cell Model Passports 1q gene-level CNV/GISTIC"),
    "MC-IXC": ("No gain", "Sanger Cell Model Passports 1q gene-level CNV/GISTIC"),
}

LOW_PARTIAL_ONEQ_CELLS = {
    "CADO-ES1": "DepMap threshold-level gain; Sanger supports only partial/low-level 1q increase",
    "SK-NEP-1": "DepMap threshold-level/fragmented 1q increase; WGS segment mean is not broad-gain level",
    "SK-PN-DW": "DepMap threshold-level gain; Sanger 1q gene-level CNV does not support a clean broad gain",
}


def depmap_download_url(file_name: str) -> str:
    metadata_file = DEP_MAP_DIR / "downloads_metadata.json"
    with metadata_file.open() as handle:
        metadata = json.load(handle)

    items: list[dict] = []

    def walk(x):
        if isinstance(x, dict):
            if x.get("fileName") and x.get("downloadUrl"):
                items.append(x)
            for value in x.values():
                walk(value)
        elif isinstance(x, list):
            for value in x:
                walk(value)

    walk(metadata)
    matches = [
        item
        for item in items
        if item.get("fileName") == file_name and "26q1" in item.get("downloadUrl", "")
    ]
    if not matches:
        raise FileNotFoundError(f"No 26Q1 DepMap download URL found for {file_name}")
    url = matches[0]["downloadUrl"]
    return DEPMAP_BASE_URL + url if url.startswith("/") else url


def ensure_depmap_file(file_name: str) -> pathlib.Path:
    out_file = DEP_MAP_DIR / file_name
    if out_file.exists() and out_file.stat().st_size > 0:
        return out_file

    DEP_MAP_DIR.mkdir(parents=True, exist_ok=True)
    url = depmap_download_url(file_name)
    subprocess.run(
        [
            "curl",
            "-L",
            "-C",
            "-",
            "-A",
            "Mozilla/5.0",
            "--fail",
            "--retry",
            "3",
            "--retry-delay",
            "5",
            "-o",
            str(out_file),
            url,
        ],
        check=True,
    )
    return out_file


def ensure_cytoband() -> pathlib.Path:
    out_file = DEP_MAP_DIR / "hg38_cytoBand.txt"
    if out_file.exists() and out_file.stat().st_size > 0:
        return out_file
    with urllib.request.urlopen(CYTOBAND_URL, timeout=60) as response:
        out_file.write_bytes(gzip.decompress(response.read()))
    return out_file


def ensure_gene_info() -> pathlib.Path:
    ANNOTATION_DIR.mkdir(parents=True, exist_ok=True)
    out_file = ANNOTATION_DIR / "Homo_sapiens.gene_info.gz"
    if out_file.exists() and out_file.stat().st_size > 0:
        return out_file
    with urllib.request.urlopen(NCBI_GENE_INFO_URL, timeout=120) as response:
        out_file.write_bytes(response.read())
    return out_file


def heatmap_models(model_file: pathlib.Path) -> pd.DataFrame:
    cell_order = (
        pd.read_csv(HEATMAP_DEPMAP_TABLE, sep="\t")["CellLineName"].drop_duplicates().tolist()
    )
    model = pd.read_csv(model_file, low_memory=False)
    models = model.loc[
        model["CellLineName"].isin(cell_order),
        ["ModelID", "CellLineName", "StrippedCellLineName", "CCLEName"],
    ].copy()
    models["CellLineName"] = pd.Categorical(
        models["CellLineName"], categories=cell_order, ordered=True
    )
    return models.sort_values("CellLineName").reset_index(drop=True)


def oneq_bounds(cytoband_file: pathlib.Path) -> tuple[int, int]:
    cytoband = pd.read_csv(
        cytoband_file,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "band", "gie"],
    )
    q_arm = cytoband[(cytoband["chrom"] == "chr1") & cytoband["band"].str.startswith("q")]
    return int(q_arm["start"].min()), int(q_arm["end"].max())


def wgs_oneq_summary(segment_file: pathlib.Path, models: pd.DataFrame) -> pd.DataFrame:
    q_start, q_end = oneq_bounds(ensure_cytoband())
    model_ids = set(models["ModelID"])
    chunks = []
    for chunk in pd.read_csv(
        segment_file,
        usecols=[
            "ModelID",
            "Chrom",
            "Start",
            "End",
            "RelCopyRatio",
            "IsDefaultEntryForModel",
        ],
        chunksize=200_000,
        low_memory=False,
    ):
        sub = chunk[
            (chunk["ModelID"].isin(model_ids))
            & (chunk["Chrom"] == "chr1")
            & (chunk["IsDefaultEntryForModel"] == "Yes")
        ].copy()
        if sub.empty:
            continue
        sub["ov_start"] = sub["Start"].clip(lower=q_start)
        sub["ov_end"] = sub["End"].clip(upper=q_end)
        sub = sub[sub["ov_end"] > sub["ov_start"]]
        if not sub.empty:
            chunks.append(sub)

    if not chunks:
        return pd.DataFrame(
            columns=[
                "ModelID",
                "oneq_wgs_relcopy_weighted_mean",
                "oneq_wgs_fraction_relcopy_ge_1_15",
                "oneq_wgs_covered_bp",
                "oneq_wgs_segments",
            ]
        )

    segments = pd.concat(chunks, ignore_index=True)
    rows = []
    for model_id, group in segments.groupby("ModelID"):
        length = (group["ov_end"] - group["ov_start"]).astype(float)
        total = length.sum()
        rows.append(
            {
                "ModelID": model_id,
                "oneq_wgs_relcopy_weighted_mean": (
                    group["RelCopyRatio"] * length
                ).sum()
                / total,
                "oneq_wgs_fraction_relcopy_ge_1_15": length[
                    group["RelCopyRatio"] >= 1.15
                ].sum()
                / total,
                "oneq_wgs_covered_bp": total,
                "oneq_wgs_segments": len(group),
            }
        )
    return pd.DataFrame(rows)


def ncbi_oneq_gene_ids(gene_info_file: pathlib.Path) -> set[str]:
    gene_info = pd.read_csv(
        gene_info_file,
        sep="\t",
        compression="gzip",
        comment="#",
        header=None,
        dtype=str,
        names=[
            "tax_id",
            "GeneID",
            "Symbol",
            "LocusTag",
            "Synonyms",
            "dbXrefs",
            "chromosome",
            "map_location",
            "description",
            "type_of_gene",
            "Symbol_from_nomenclature_authority",
            "Full_name_from_nomenclature_authority",
            "Nomenclature_status",
            "Other_designations",
            "Modification_date",
            "Feature_type",
        ],
    )
    return set(
        gene_info.loc[
            gene_info["map_location"].fillna("").str.startswith("1q"), "GeneID"
        ].astype(str)
    )


def gene_level_oneq_summary(cn_gene_file: pathlib.Path, models: pd.DataFrame) -> pd.DataFrame:
    oneq_ids = ncbi_oneq_gene_ids(ensure_gene_info())
    header = pd.read_csv(cn_gene_file, nrows=0).columns.tolist()
    oneq_cols = []
    for column in header[1:]:
        match = re.search(r"\((\d+)\)$", column)
        if match and match.group(1) in oneq_ids:
            oneq_cols.append(column)

    cn = pd.read_csv(cn_gene_file, usecols=["Unnamed: 0", *oneq_cols])
    cn = cn[cn["Unnamed: 0"].isin(set(models["ModelID"]))].rename(
        columns={"Unnamed: 0": "ModelID"}
    )
    score = cn.set_index("ModelID")[oneq_cols]
    denominator = score.notna().sum(axis=1)
    return pd.DataFrame(
        {
            "ModelID": score.index,
            "oneq_gene_cn_log2_mean": score.mean(axis=1, skipna=True).values,
            "oneq_gene_cn_log2_median": score.median(axis=1, skipna=True).values,
            "oneq_gene_cn_gene_count": denominator.values,
            "oneq_gene_cn_fraction_ge_1_10": (
                score.ge(1.10).sum(axis=1) / denominator
            ).values,
            "oneq_gene_cn_fraction_ge_1_15": (
                score.ge(1.15).sum(axis=1) / denominator
            ).values,
        }
    )


def mutation_summary(mutation_file: pathlib.Path, models: pd.DataFrame) -> pd.DataFrame:
    model_ids = set(models["ModelID"])
    usecols = [
        "ModelID",
        "IsDefaultEntryForModel",
        "HugoSymbol",
        "ProteinChange",
        "DNAChange",
        "VepImpact",
    ]
    chunks = []
    variant_counts = {model_id: 0 for model_id in model_ids}
    for chunk in pd.read_csv(
        mutation_file,
        usecols=usecols,
        chunksize=250_000,
        low_memory=False,
    ):
        default = chunk[
            (chunk["ModelID"].isin(model_ids)) & (chunk["IsDefaultEntryForModel"] == "Yes")
        ]
        counts = default["ModelID"].value_counts()
        for model_id, count in counts.items():
            variant_counts[model_id] += int(count)
        sub = default[default["HugoSymbol"].isin(["TP53", "STAG2"])].copy()
        if not sub.empty:
            chunks.append(sub)

    mutations = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame(columns=usecols)
    mutation_rows_file = TABLE_DIR / "depmap_ewing_tp53_stag2_mutation_rows.tsv"
    mutations.to_csv(mutation_rows_file, sep="\t", index=False)

    if mutations.empty:
        summary = pd.DataFrame(columns=["ModelID", "HugoSymbol"])
    else:
        mutations["mutation_label"] = mutations.apply(
            lambda row: (
                str(row["ProteinChange"])
                if pd.notna(row["ProteinChange"])
                else str(row["DNAChange"])
            )
            + " ["
            + str(row["VepImpact"])
            + "]",
            axis=1,
        )
        summary = (
            mutations.groupby(["ModelID", "HugoSymbol"])
            .agg(
                mutation_count=("HugoSymbol", "size"),
                mutation_details=("mutation_label", lambda x: "; ".join(x.astype(str))),
            )
            .reset_index()
        )

    wide = models.copy()
    for gene in ["TP53", "STAG2"]:
        gene_summary = summary[summary["HugoSymbol"] == gene].rename(
            columns={
                "mutation_count": f"{gene}_mutation_count",
                "mutation_details": f"{gene}_mutation_details",
            }
        )
        wide = wide.merge(
            gene_summary[
                ["ModelID", f"{gene}_mutation_count", f"{gene}_mutation_details"]
            ],
            on="ModelID",
            how="left",
        )
        wide[f"{gene}_mutation_count"] = (
            wide[f"{gene}_mutation_count"].fillna(0).astype(int)
        )
        wide[f"{gene}_mutated"] = wide[f"{gene}_mutation_count"] > 0
        wide[f"{gene}_mutation_details"] = wide[f"{gene}_mutation_details"].fillna("")

    availability = pd.DataFrame(
        {
            "ModelID": list(variant_counts.keys()),
            "mutation_profile_variant_rows": list(variant_counts.values()),
        }
    )
    availability["mutation_profile_available"] = (
        availability["mutation_profile_variant_rows"] > 0
    )
    return wide.merge(availability, on="ModelID", how="left")


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    model_file = ensure_depmap_file("Model.csv")
    mutation_file = ensure_depmap_file("OmicsSomaticMutations.csv")
    segment_file = ensure_depmap_file("OmicsCNSegmentsWGS.csv")
    cn_gene_file = ensure_depmap_file("PortalOmicsCNGeneLog2.csv")

    models = heatmap_models(model_file)
    gene_cn = gene_level_oneq_summary(cn_gene_file, models)
    wgs_cn = wgs_oneq_summary(segment_file, models)
    mutations = mutation_summary(mutation_file, models)

    final = (
        models.merge(gene_cn, on="ModelID", how="left")
        .merge(wgs_cn, on="ModelID", how="left")
        .merge(
            mutations.drop(columns=["StrippedCellLineName", "CCLEName", "CellLineName"]),
            on="ModelID",
            how="left",
        )
    )
    final["oneq_gain_status"] = np.where(
        final["oneq_gene_cn_log2_mean"].isna(),
        "CN missing",
        np.where(
            final["oneq_gene_cn_log2_mean"] >= ONEQ_GAIN_GENE_CN_MEAN_THRESHOLD,
            "Gain",
            "No gain",
        ),
    )
    final["oneq_gain_basis"] = (
        "DepMap PortalOmicsCNGeneLog2 mean across NCBI genes mapped to 1q; "
        f"Gain if mean >= {ONEQ_GAIN_GENE_CN_MEAN_THRESHOLD:.2f}"
    )
    final["oneq_gain_status_curated"] = final["oneq_gain_status"]
    final["oneq_gain_curated_basis"] = final["oneq_gain_basis"]
    for cell_line, (status, basis) in CURATED_ONEQ_STATUS_OVERRIDES.items():
        mask = final["CellLineName"] == cell_line
        final.loc[mask, "oneq_gain_status_curated"] = status
        final.loc[mask, "oneq_gain_curated_basis"] = basis
    for cell_line, basis in LOW_PARTIAL_ONEQ_CELLS.items():
        mask = final["CellLineName"] == cell_line
        final.loc[mask, "oneq_gain_status_curated"] = "Low/partial"
        final.loc[mask, "oneq_gain_curated_basis"] = basis

    for gene in ["TP53", "STAG2"]:
        final[f"{gene}_mutated_curated"] = final[f"{gene}_mutated"]
        final[f"{gene}_mutation_details_curated"] = final[f"{gene}_mutation_details"]
        final[f"{gene}_curated_basis"] = "DepMap Public 26Q1 default mutation profile"

    for cell_line, gene_overrides in CURATED_MUTATION_OVERRIDES.items():
        mask = final["CellLineName"] == cell_line
        for gene, detail in gene_overrides.items():
            detail_col = f"{gene}_mutation_details_curated"
            basis_col = f"{gene}_curated_basis"
            final.loc[mask, f"{gene}_mutated_curated"] = True
            existing = final.loc[mask, detail_col].fillna("")
            final.loc[mask, detail_col] = existing.apply(
                lambda value: "; ".join(part for part in [value, detail] if part)
            )
            final.loc[mask, basis_col] = (
                "Curated external evidence: Cellosaurus/Sanger/literature cross-check"
            )

    concise_columns = [
        "CellLineName",
        "ModelID",
        "oneq_gain_status",
        "oneq_gain_status_curated",
        "oneq_gene_cn_log2_mean",
        "oneq_gene_cn_fraction_ge_1_10",
        "oneq_wgs_relcopy_weighted_mean",
        "oneq_wgs_fraction_relcopy_ge_1_15",
        "oneq_gain_curated_basis",
        "TP53_mutated",
        "TP53_mutation_details",
        "TP53_mutated_curated",
        "TP53_mutation_details_curated",
        "TP53_curated_basis",
        "STAG2_mutated",
        "STAG2_mutation_details",
        "STAG2_mutated_curated",
        "STAG2_mutation_details_curated",
        "STAG2_curated_basis",
        "mutation_profile_available",
        "mutation_profile_variant_rows",
        "oneq_gain_basis",
    ]
    final[concise_columns].to_csv(
        TABLE_DIR / "depmap_ewing_1q_tp53_stag2_status.tsv",
        sep="\t",
        index=False,
    )
    gene_cn.to_csv(TABLE_DIR / "depmap_ewing_1q_gene_cn_scores.tsv", sep="\t", index=False)

    print(
        "Wrote",
        TABLE_DIR / "depmap_ewing_1q_tp53_stag2_status.tsv",
        "for",
        len(final),
        "cell lines.",
    )


if __name__ == "__main__":
    main()
