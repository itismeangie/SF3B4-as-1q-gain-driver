# SF3B4 as 1q gain driver

Computational-analysis code package for the manuscript. The repository is organized by manuscript topic and intentionally excludes plotting scripts, generated figures, workbooks, local caches, previously published assay code, and sensitive human subject files.

Current version: `0.1.0-prepublication`.

## Topics

- `01_pan_cancer_1q/`: pan-cancer 1q gain and chr1q gene context.
- `02_ewing_cna_survival/`: Ewing sarcoma CNA, survival, and relapse enrichment.
- `03_inform_transcriptome_stemness/`: INFORM 1q gain transcriptome and stemness.
- `04_netbid2_driver_prioritization/`: NetBID2 driver prioritization and DepMap annotation.
- `05_sf3b4_perturbation_proteomics/`: SF3B4 perturbation and proteomics.
- `06_xenograft_growth/`: SF3B4 xenograft growth.

Each topic folder contains its respective computational scripts, a local `run.sh` wrapper, and a README. Steps requiring controlled or external inputs are skipped unless the corresponding environment variables are set.

## Data Boundary

Sensitive human subject and controlled-access inputs are not included. Public documentation may mention required input filenames, but not their contents. Set the environment variables documented in `docs/MANUSCRIPT_TOPIC_CODE_BLOCKS.md` when rerunning controlled-data analyses locally.

## Publication Metadata

- Version: `VERSION`
- Changelog: `CHANGELOG.md`
- License: `LICENSE`
- Authors: `AUTHORS.md`
- Citation metadata: `CITATION.cff`
- Zenodo metadata template: `.zenodo.json`
- Dependency notes: `docs/DEPENDENCIES.md`
- External input manifest: `docs/EXTERNAL_INPUTS.md`
- Draft code-availability wording: `docs/CODE_AVAILABILITY.md`
- Pre-publication checklist: `docs/RELEASE_CHECKLIST.md`

## Quick Start

Running a topic wrapper without the required environment variables is safe: controlled-data and public-download steps are skipped by default.

```bash
bash 01_pan_cancer_1q/run.sh
RUN_PUBLIC_DOWNLOADS=1 bash 01_pan_cancer_1q/run.sh
COMBINED_COHORT_DF=/path/to/combined_df.csv bash 02_ewing_cna_survival/run.sh
INFORM_ALL_ROOT=/path/to/inform_all INFORM_TPM_SYMBOL_CSV=/path/to/INFORM_EWS_BCOR_CIC_merged_TPM_SYMBOL.csv bash 03_inform_transcriptome_stemness/run.sh
RUN_PUBLIC_DOWNLOADS=1 bash 04_netbid2_driver_prioritization/run.sh
RUN_PROTEOMICS_TABLES=1 bash 05_sf3b4_perturbation_proteomics/run.sh
REBUILD_MOUSE_TABLES=1 bash 06_xenograft_growth/run.sh
```

See `docs/MANUSCRIPT_TOPIC_CODE_BLOCKS.md` for the expanded command blocks and controlled-input environment variables.
