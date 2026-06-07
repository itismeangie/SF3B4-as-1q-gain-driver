#!/usr/bin/env python3
"""
Rewrite the NetBID2-generated run_tf/run_sig scripts to use the modern
SJARACNe CLI (which requires the 'local' subcommand) instead of the legacy
four-argument wrapper.
"""
from __future__ import annotations

import os
import shlex
from pathlib import Path

ROOT = Path(os.environ.get("NETBID_RUN_ROOT", "netbid2_runs"))


def resolve_new_script(cmd_parts: list[str]) -> str | None:
    """
    Map `sjaracne JOB EXP HUB BASE` to the modern invocation.
    """
    if len(cmd_parts) < 5 or cmd_parts[0] != "sjaracne":
        return None

    _, job_name, exp_file, hub_file, base_dir, *extra = cmd_parts
    base_dir = base_dir.rstrip("/")
    if not base_dir:
        base_dir = "."

    out_dir = f"{base_dir}/{job_name}"
    tmp_dir = f"{out_dir}_tmp"

    extra_args = ""
    if extra:
        extra_args = " \\\n  " + " \\\n  ".join(shlex.quote(arg) for arg in extra)

    new_script = (
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n\n"
        f'sjaracne local \\\n'
        f'  -e "{exp_file}" \\\n'
        f'  -g "{hub_file}" \\\n'
        f'  -o "{out_dir}" \\\n'
        f'  -tmp "{tmp_dir}"{extra_args}\n'
    )
    return new_script


def main() -> None:
    scripts = sorted(ROOT.rglob("run_*.sh"))
    if not scripts:
        raise SystemExit("No run_*.sh files found under netbid2_runs/")

    for script in scripts:
        original = script.read_text().strip()
        if not original:
            print(f"Skipping empty script: {script}")
            continue

        try:
            parts = shlex.split(original)
        except ValueError as exc:  # malformed quoting
            print(f"Skipping {script}: could not parse line ({exc})")
            continue

        new_contents = resolve_new_script(parts)
        if new_contents is None:
            print(f"Skipping {script}: unexpected command format")
            continue

        script.write_text(new_contents)
        os.chmod(script, 0o755)
        print(f"Updated {script}")


if __name__ == "__main__":
    main()
