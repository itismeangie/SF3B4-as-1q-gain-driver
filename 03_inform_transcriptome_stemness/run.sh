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

if [[ -n "${INFORM_ALL_ROOT:-}" && -d "${INFORM_ALL_ROOT}" ]]; then
  run Rscript dge_1q_gain_gsea.R \
    --root="${INFORM_ALL_ROOT}" \
    --out_root="dge_1q_gain_gsea"
else
  skip "controlled INFORM GSEA recomputation; set INFORM_ALL_ROOT to the local controlled INFORM analysis root."
fi

if [[ -n "${INFORM_ALL_ROOT:-}" && -d "${INFORM_ALL_ROOT}" && -n "${INFORM_TPM_SYMBOL_CSV:-}" && -f "${INFORM_TPM_SYMBOL_CSV}" ]]; then
  run Rscript inform_mRNAsi_scoring.R \
    --root="${INFORM_ALL_ROOT}" \
    --expr_matrix="${INFORM_TPM_SYMBOL_CSV}" \
    --out_dir="stemness_outputs"
else
  skip "controlled mRNAsi recomputation; set INFORM_ALL_ROOT and INFORM_TPM_SYMBOL_CSV."
fi
