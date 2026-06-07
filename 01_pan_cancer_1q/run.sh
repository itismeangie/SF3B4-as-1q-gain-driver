#!/usr/bin/env bash
set -euo pipefail

topic_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${topic_dir}"

run() {
  printf '\n## %s\n' "$*"
  "$@"
}

skip() {
  printf '\nSkipping %s\n' "$1"
}

if [[ "${RUN_TCGA_PROGNOSIS:-0}" == "1" ]]; then
  run Rscript tcga_1q_prognosis_meta_analysis.R
else
  skip "TCGA prognosis meta-analysis; set RUN_TCGA_PROGNOSIS=1 after providing local TCGA/reference inputs."
fi

if [[ "${RUN_PUBLIC_DOWNLOADS:-0}" == "1" || "${RUN_1Q_CONTEXT:-0}" == "1" ]]; then
  run python3 analyze_arm_gene_categories.py
  run python3 analyze_1q_pathway_signed_roles.py
else
  skip "public-reference 1q context summaries; set RUN_PUBLIC_DOWNLOADS=1 or RUN_1Q_CONTEXT=1 to download public references and run them."
fi
