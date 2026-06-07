#!/usr/bin/env bash

set -euo pipefail

# Ensure the user-level bin directory is available for locally installed tools.
LOCAL_BIN="$HOME/.local/bin"
if [[ -d "$LOCAL_BIN" && ":$PATH:" != *":$LOCAL_BIN:"* ]]; then
  export PATH="$LOCAL_BIN:$PATH"
fi

# Directory containing the NetBID2 per-project run folders.
NETBID_RUN_ROOT=${NETBID_RUN_ROOT:-"netbid2_runs"}

if ! command -v sjaracne >/dev/null 2>&1; then
  echo "Error: sjaracne is not available on PATH." >&2
  exit 1
fi

log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*"
}

declare -a projects

if (( $# > 0 )); then
  projects=("$@")
else
  if [[ ! -d "$NETBID_RUN_ROOT" ]]; then
    log "No NetBID2 run root found at ${NETBID_RUN_ROOT}"
    exit 0
  fi
  mapfile -t projects < <(find "$NETBID_RUN_ROOT" -mindepth 1 -maxdepth 1 -type d -exec basename {} \; | sort)
fi

if (( ${#projects[@]} == 0 )); then
  log "No NetBID2 projects found under ${NETBID_RUN_ROOT}"
  exit 0
fi

for project in "${projects[@]}"; do
  project_dir="${NETBID_RUN_ROOT}/${project}"
  run_dir="${project_dir}/SJAR/${project}"

  if [[ ! -d "$run_dir" ]]; then
    log "Skipping ${project}: no SJAR directory at ${run_dir}"
    continue
  fi

  for mode in tf sig; do
    script="${run_dir}/run_${mode}.sh"
    log_file="${run_dir}/run_${mode}.log"
    label=$(echo "$mode" | tr '[:lower:]' '[:upper:]')

    if [[ ! -f "$script" ]]; then
      log "Skipping ${project} ${label}: missing ${script}"
      continue
    fi

    log "Running ${project} ${label} network (logging to ${log_file})"
    bash "$script" |& tee "$log_file"
    log "Finished ${project} ${label}"
  done
done
