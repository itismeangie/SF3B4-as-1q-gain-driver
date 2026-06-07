# Rerunning Analyses

Run analysis scripts from their topic folder unless a script-specific README says otherwise. Each top-level topic directory contains its own scripts and local `run.sh`.

## Current State

The repository is not a raw-data archive and is not a plotting-code release. It contains computational scripts only. Scripts that rebuild analyses from raw RNA-seq, TCGA source downloads, raw NetBID2 expression inputs, raw proteomics-derived tables, animal workbooks, or external genome annotations still need those external files.

Install/runtime notes are in `docs/DEPENDENCIES.md`. Required external inputs are summarized in `docs/EXTERNAL_INPUTS.md`.

Use this audit to list hardcoded local paths and whether they resolve:

```bash
python3 docs/check_rerun_readiness.py --write docs/rerun_path_audit.tsv
```

Current unresolved items are expected external or controlled-input cases.

## Useful Environment Variables

- `ONE_Q_REFERENCE_CACHE`: cache directory for large NCBI/UCSC reference downloads used by `01_pan_cancer_1q` Python scripts. Defaults to the script-defined cache path unless overridden.
- `INFORM_ALL_ROOT`: external INFORM analysis root used by controlled-data recomputation scripts.
- `INFORM_TPM_SYMBOL_CSV`: INFORM TPM matrix used by selected 1q-gene expression summaries.
- `COMBINED_COHORT_DF`: controlled combined-cohort CNA/survival input named `combined_df.csv`.
- `INFORM_MANUAL_FINAL_CSV`: controlled INFORM arm-call table named `manual_final_inform.csv`.
- `NETBID_RUN_ROOT`: external NetBID2 full run tree for `04_netbid2_driver_prioritization/run_sjaracne_networks.sh`.
- `NETBID_TCGA_BASE_DIR`: external TCGA/NetBID2 base folder for NetBID2 preparation and adjusted-model scripts.
- `REBUILD_MOUSE_TABLES`: set to `1` when approved local animal source workbooks are available.
- `RUN_PROTEOMICS_TABLES`: set to `1` when approved local decontaminated proteomics tables and BioGRID resources are available.
- `LOCAL_R_LIB`: optional extra R library path for scripts that used to assume a user-specific `.libPaths()` entry.

## Repo-Local Fixes Already Applied

- `01_pan_cancer_1q` Python scripts resolve lightweight `data/` and `results/` relative to their own topic folder, while large downloaded genome/reference tables are cached outside the repository through `ONE_Q_REFERENCE_CACHE`.
- `02_ewing_cna_survival` can use `COMBINED_COHORT_DF` for an external controlled `combined_df.csv` instead of tracking that file in the public repository.
- INFORM scripts accept explicit `--root`, `--expr_matrix`, and output arguments so controlled inputs can stay external.
- NetBID2 scripts expose TCGA/NetBID2 paths as command-line arguments or environment variables.
