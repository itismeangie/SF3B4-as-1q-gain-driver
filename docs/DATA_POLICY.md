# Data Tracking Policy

This repository tracks computational analysis code and filename-only provenance for required external inputs. It does not track generated figures, sensitive data contents, or raw-data archives.

For a public GitHub release, sensitive human subject or controlled-access data must not be included. It is acceptable to document the names of required controlled input files, but not their row-level contents. This includes patient IDs, internal registry IDs, sample IDs, dates, birth years, treatment-center fields, survival/event fields, sample maps, and any table that can link a sample to clinical outcome or controlled sequencing data.

It does not track:

- Raw or aligned RNA-seq data: `fastq`, `bam`, `cram`, `bigWig`, STAR/RSEM full run outputs, genome indices, and compressed RNA-seq archives.
- Raw TCGA data: mutation lists, MC3 MAF files, probe-level or bin-level source extracts, and raw GDC/Firehose-style source tables.
- Raw proteomics instrument output: Bruker `.d` directories, `analysis.tdf`, `analysis.tdf_bin`, chromatogram SQLite files, and related raw acquisition files.
- Raw NetBID2 expression inputs and intermediate workspaces: `.exp`, `.RData`, `.rds`, and equivalent run-state files.
- Local execution state: caches, package libraries, virtual environments, `.RData`, `.Rhistory`, `__pycache__`, logs, and temporary files.
- Local archive metadata and compressed source dumps: `__MACOSX`, `.DS_Store`, `.gz`, `.zip`, and similar sidecar/archive files.
- Serialized R work objects: `.rds` and `.RData` files are treated as regenerable binary state, not source-controlled analysis records.
- Large binary figure/source files such as TIFF and HEIC image sets unless explicitly curated later.
- Any file larger than the importer limit, currently 25 MiB.
- Human patient/sample-level clinical, survival, registry, molecular-sample-map, or controlled-access metadata. Public documentation may list the file names and expected roles only.

The goal is a reproducible, readable public computational-code repository, not a raw-data archive or controlled-data mirror. Heavy/raw and sensitive assets should stay in approved external storage and be referenced through filename-only manifests.

## Review Rules

Before public release, keep only lightweight source-code and documentation files in this repository. Do not add:

- Controlled or patient/sample-level data contents.
- Raw sequencing, raw proteomics, raw NetBID2, or source animal workbook files.
- Generated figures, local caches, or package libraries.
- Files that can be regenerated from controlled inputs unless explicitly approved as deidentified derived tables.

If a future analysis requires a file currently excluded, add a deidentified aggregate or controlled-access retrieval note rather than committing the raw, linked, or sensitive source file.
