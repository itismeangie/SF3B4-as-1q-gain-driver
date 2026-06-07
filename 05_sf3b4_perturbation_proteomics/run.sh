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

if [[ "${RUN_PROTEOMICS_TABLES:-0}" == "1" ]]; then
  run Rscript identify_decontam_stringent_concordant_downregulated.R
  run Rscript summarise_publication_responsive_hits.R
  run python3 analyze_biogrid_volcano_down_pairwise_interactions.py
else
  skip "proteomics table recomputation; set RUN_PROTEOMICS_TABLES=1 after providing local decontaminated DEP outputs and BioGRID resources."
fi
