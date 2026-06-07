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

if [[ -n "${COMBINED_COHORT_DF:-}" && -f "${COMBINED_COHORT_DF}" ]]; then
  run env COMBINED_COHORT_DF="${COMBINED_COHORT_DF}" Rscript compute_CNA_survival_models.R
else
  skip "controlled recomputation; set COMBINED_COHORT_DF to the local controlled file named combined_df.csv."
fi
