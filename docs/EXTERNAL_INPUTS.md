# External Inputs

This repository does not include sensitive human subject data, controlled-access data, raw data, raw proteomics files, NetBID2 workspaces, animal source workbooks, or generated figures. File names are documented for provenance only.

## 01 Pan-Cancer 1q

Public/reference inputs:

- NCBI `gene_info.gz`
- NCBI `gene2refseq.gz`
- NCBI `gene2go.gz`
- UCSC `cytoBand.txt.gz`
- UCSC `hg38.chrom.sizes`
- OncoKB curated gene JSON
- KEGG pathway resources

The 1q gene scripts download public reference files as needed. Large files should be cached outside Git with `ONE_Q_REFERENCE_CACHE`.

The topic runner skips these downloads unless `RUN_PUBLIC_DOWNLOADS=1` or `RUN_1Q_CONTEXT=1` is set.

TCGA prognosis recomputation requires local TCGA/reference inputs such as:

- `PANCAN_ArmCallsAndAneuploidyScore_092817.txt`
- `infiltration_estimation_for_tcga.csv`
- `TCGA-CDR-SupplementalTableS1.xlsx`
- `Aran_TCGA_purity_Supplementary_Data1.xlsx`
- TCGA mutation manifest and mutation covariate files

## 02 Ewing CNA And Survival

Controlled input:

- `combined_df.csv`

This table contains cohort survival/CNA records and must not be committed publicly.

## 03 INFORM Transcriptome And Stemness

Controlled inputs include:

- `INFORM_EWS_BCOR_CIC_merged_TPM_SYMBOL.csv`
- `manual_final_inform.csv`
- `sample_info.tsv`
- `covariates_full.tsv`
- `mRNAsi_SC_PCBC_stemSig_weights.tsv`
- `tp53_cdkn2a_specimen_status_resolved.tsv`
- `cnv_mapping_loss_breakdown.csv`
- immune-deconvolution and methylation-derived local metadata files, where used

These files may link samples to controlled molecular or clinical data and must stay outside Git.

## 04 NetBID2 Driver Prioritization

Public-download inputs:

- DepMap Public 26Q1 model, copy-number, mutation, and gene-effect tables
- UCSC hg38 cytoband table
- NCBI human gene information

External local inputs:

- TCGA TPM matrices
- TCGA arm-call table
- NetBID2 run workspaces
- SJAracne network directories
- TCGA clinical/covariate files

Raw NetBID2 expression inputs and run-state files stay outside Git.

## 05 SF3B4 Perturbation Proteomics

External local inputs:

- Decontaminated DEP output tables
- Imputed and non-imputed responsive-hit tables
- BioGRID Homo sapiens release resources

Raw proteomics instrument files and BioGRID caches stay outside Git.

## 06 Xenograft Growth

External local input:

- `mouse_documentation_Angelina_090226_B410_Geyer.xlsx`

Review animal-level identifiers before any public data release. This public repository keeps the code only.
