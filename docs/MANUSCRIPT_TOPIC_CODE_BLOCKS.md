# Manuscript Topic Code Blocks

Use these blocks from the repository root. The top-level folders are the manuscript topics, and each folder contains its own computational scripts. Plotting, figure-rendering, previously published assay code, generated figures, and sensitive data contents are intentionally excluded.

Sensitive human subject and controlled-access inputs must stay outside Git. Filename-only references are allowed for provenance. Use environment variables to point to local controlled storage when rerunning analyses that require those files.

## Setup

```bash
git clone https://github.com/itismeangie/SF3B4-as-1q-gain-driver.git "SF3B4 as 1q gain driver"
cd "SF3B4 as 1q gain driver"
```

```bash
export CONTROLLED_DATA_DIR="/path/to/controlled-inputs"
export INFORM_ALL_ROOT="${CONTROLLED_DATA_DIR}/inform_all"
export INFORM_TPM_SYMBOL_CSV="${CONTROLLED_DATA_DIR}/INFORM_EWS_BCOR_CIC_merged_TPM_SYMBOL.csv"
export COMBINED_COHORT_DF="${CONTROLLED_DATA_DIR}/combined_df.csv"
export NETBID_TCGA_BASE_DIR="${CONTROLLED_DATA_DIR}/tcga_netbid2"
export NETBID_RUN_ROOT="${NETBID_TCGA_BASE_DIR}/netbid2_runs"
```

## Topic Runner Scripts

```bash
bash 01_pan_cancer_1q/run.sh
RUN_PUBLIC_DOWNLOADS=1 bash 01_pan_cancer_1q/run.sh
COMBINED_COHORT_DF="${COMBINED_COHORT_DF}" bash 02_ewing_cna_survival/run.sh
INFORM_ALL_ROOT="${INFORM_ALL_ROOT}" INFORM_TPM_SYMBOL_CSV="${INFORM_TPM_SYMBOL_CSV}" bash 03_inform_transcriptome_stemness/run.sh
RUN_PUBLIC_DOWNLOADS=1 bash 04_netbid2_driver_prioritization/run.sh
RUN_PROTEOMICS_TABLES=1 bash 05_sf3b4_perturbation_proteomics/run.sh
REBUILD_MOUSE_TABLES=1 bash 06_xenograft_growth/run.sh
```

## Topic 1: Pan-Cancer 1q Gain And Chr1q Gene Context

```bash
cd 01_pan_cancer_1q
export RUN_PUBLIC_DOWNLOADS=1
python3 analyze_arm_gene_categories.py
python3 analyze_1q_pathway_signed_roles.py
```

```bash
cd 01_pan_cancer_1q
Rscript tcga_1q_prognosis_meta_analysis.R \
  --tcga_root="/path/to/TCGA" \
  --reference_root="/path/to/tcga_clinical" \
  --out_root="tcga_1q_prognosis"
```

## Topic 2: Ewing Sarcoma CNA, Survival, And Relapse Enrichment

```bash
cd 02_ewing_cna_survival
COMBINED_COHORT_DF="/path/to/combined_df.csv" \
Rscript compute_CNA_survival_models.R
```

Do not commit `combined_df.csv` or registry/sample-map inputs. Keep only filename-level provenance in documentation.

## Topic 3: INFORM 1q Gain Transcriptome And Stemness

```bash
cd 03_inform_transcriptome_stemness
Rscript dge_1q_gain_gsea.R \
  --root="${INFORM_ALL_ROOT}" \
  --out_root="dge_1q_gain_gsea"
```

```bash
cd 03_inform_transcriptome_stemness
Rscript inform_mRNAsi_scoring.R \
  --root="${INFORM_ALL_ROOT}" \
  --expr_matrix="${INFORM_TPM_SYMBOL_CSV}" \
  --out_dir="stemness_outputs"
```

## Topic 4: NetBID2 Driver Prioritization And DepMap Annotation

```bash
cd 04_netbid2_driver_prioritization
python3 annotate_depmap_ewing_1q_tp53_stag2.py
```

```bash
cd 04_netbid2_driver_prioritization
Rscript netbid2_sensitivity_prepare.R \
  --tcga_base_dir="${NETBID_TCGA_BASE_DIR}"
```

```bash
cd 04_netbid2_driver_prioritization
Rscript netbid2_adjusted_best_models.R \
  --tcga_base_dir="${NETBID_TCGA_BASE_DIR}" \
  --netbid_dir="${NETBID_RUN_ROOT}"
```

```bash
cd 04_netbid2_driver_prioritization
NETBID_RUN_ROOT="${NETBID_RUN_ROOT}" python3 rewrite_sjaracne_scripts.py
NETBID_RUN_ROOT="${NETBID_RUN_ROOT}" bash run_sjaracne_networks.sh
```

Raw NetBID2 workspaces and controlled TCGA/clinical inputs stay outside Git.

## Topic 5: SF3B4 Perturbation And Proteomics

```bash
cd 05_sf3b4_perturbation_proteomics
Rscript identify_decontam_stringent_concordant_downregulated.R
Rscript summarise_publication_responsive_hits.R
python3 analyze_biogrid_volcano_down_pairwise_interactions.py
```

Raw proteomics files, decontaminated DEP output tables, and public-resource caches stay outside Git unless separately approved as deidentified derived tables.

## Topic 6: SF3B4 Xenograft Growth

```bash
cd 06_xenograft_growth
python3 rebuild_mouse_measurements.py
```

Review animal-level identifiers before public release. Keep source workbooks outside the public code package unless separately approved.

## Audit

```bash
python3 docs/check_rerun_readiness.py --write docs/rerun_path_audit.tsv
```

```bash
git check-ignore -v \
  02_ewing_cna_survival/combined_df.csv \
  metadata/inform_target/MARVINextract_Inform_EwingSarcoma_20240807.csv \
  metadata/inform_target/inform_rnaseq_sample_map.csv \
  metadata/inform_target/manual_final_inform.csv
```
