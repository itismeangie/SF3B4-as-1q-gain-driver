# GitHub Publication Readiness

Audit date: 2026-06-07

## Current Release Shape

This derived folder is a public-facing computational-analysis code package. It is separate from the original working repository and intentionally excludes:

- Sensitive human subject or controlled-access data contents.
- Generated figures and figure-rendering scripts.
- Previously published assay code and outputs.
- Raw proteomics, raw sequencing, raw NetBID2 state, and large public-resource caches.
- Local operating-system, editor, package-cache, and execution-state files.

Filename-only references to controlled inputs are acceptable for provenance.

## Topic Map

| Topic | Included computational entry points | Data boundary |
|---|---|---|
| Pan-cancer 1q and chr1q gene context | `01_pan_cancer_1q/` | TCGA/reference source files stay external; public reference downloads can be cached locally. |
| Ewing CNA and survival | `02_ewing_cna_survival/` | `combined_df.csv` stays outside Git. |
| INFORM transcriptome and stemness | `03_inform_transcriptome_stemness/` | INFORM expression, clinical, CNV, and sample-map inputs stay outside Git. |
| NetBID2 driver prioritization | `04_netbid2_driver_prioritization/` | NetBID2 workspaces and controlled TCGA covariates stay external; DepMap resources are download-on-demand. |
| SF3B4 perturbation proteomics | `05_sf3b4_perturbation_proteomics/` | Raw proteomics and decontaminated DEP tables stay external unless separately approved as derived tables. |
| Xenograft growth | `06_xenograft_growth/` | Source workbooks and animal-level records stay external unless separately approved. |

## Remaining Pre-Push Items

1. Confirm the MIT license in `LICENSE` is approved by the institution and coauthors, or replace it with the approved license.
2. Confirm the author list in `AUTHORS.md`, `CITATION.cff`, and `.zenodo.json` against the accepted manuscript proofs.
3. Run `python3 docs/check_rerun_readiness.py --write docs/rerun_path_audit.tsv` and classify unresolved paths as controlled, public-download, generated, or legacy.
4. Confirm no sensitive files are present using the file checks in `docs/MANUSCRIPT_TOPIC_CODE_BLOCKS.md`.
5. Review `docs/sensitive_data_file_inventory.tsv`; it should remain filename-only.
6. Archive a tagged release with Zenodo and add the DOI to `docs/CODE_AVAILABILITY.md`.
