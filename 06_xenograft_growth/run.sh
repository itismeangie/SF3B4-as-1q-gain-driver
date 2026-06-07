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

if [[ "${REBUILD_MOUSE_TABLES:-0}" == "1" ]]; then
  run python3 rebuild_mouse_measurements.py
else
  skip "mouse-table rebuild; set REBUILD_MOUSE_TABLES=1 if local source workbooks are available."
fi
