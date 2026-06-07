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

netbid_tcga_args=()
if [[ -n "${NETBID_TCGA_BASE_DIR:-}" ]]; then
  netbid_tcga_args+=(--tcga_base_dir="${NETBID_TCGA_BASE_DIR}")
fi

netbid_adjusted_args=("${netbid_tcga_args[@]}")
if [[ -n "${NETBID_RUN_ROOT:-}" ]]; then
  netbid_adjusted_args+=(--netbid_dir="${NETBID_RUN_ROOT}")
fi

if [[ "${RUN_PUBLIC_DOWNLOADS:-0}" == "1" ]]; then
  run python3 annotate_depmap_ewing_1q_tp53_stag2.py
else
  skip "public DepMap download/update; set RUN_PUBLIC_DOWNLOADS=1 to refresh DepMap annotations."
fi

if [[ "${RUN_NETBID2_SENSITIVITY_PREP:-0}" == "1" ]]; then
  run Rscript netbid2_sensitivity_prepare.R "${netbid_tcga_args[@]}"
else
  skip "NetBID2 sensitivity input preparation; set RUN_NETBID2_SENSITIVITY_PREP=1 after providing local TCGA/NetBID2 inputs."
fi

if [[ "${RUN_NETBID2_ADJUSTED_MODELS:-0}" == "1" ]]; then
  run Rscript netbid2_adjusted_best_models.R "${netbid_adjusted_args[@]}"
else
  skip "NetBID2 adjusted best-model analysis; set RUN_NETBID2_ADJUSTED_MODELS=1 after providing local TCGA/NetBID2 inputs."
fi
