# Pan-Cancer 1q Gain And Chr1q Gene Context

Manuscript panels: Fig. 1a-c.

Code:

- `tcga_1q_prognosis_meta_analysis.R`
- `analyze_1q_developmental_enrichment.py`
- `analyze_1q_oncogenic_pathways.py`
- `analyze_arm_gene_categories.py`
- `analyze_1q_pathway_signed_roles.py`

The 1q gene-context scripts download public reference files to a cache. The local runner skips those downloads unless `RUN_PUBLIC_DOWNLOADS=1` or `RUN_1Q_CONTEXT=1` is set. Set `RUN_TCGA_PROGNOSIS=1` only after providing local TCGA/reference inputs for the prognosis meta-analysis.
