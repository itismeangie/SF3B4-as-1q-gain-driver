#!/usr/bin/env python3
"""Audit analysis code for repo-local rerun blockers."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]

CODE_SUFFIXES = {".R", ".Rmd", ".py", ".sh"}
LOCAL_PATH_RE = re.compile(
    r"""(?P<quote>["'])(?P<path>(?:/(?:home|media|Users)/[^"']+|~/[^"']+))(?P=quote)"""
)
OUTPUT_CONTEXT_RE = re.compile(
    r"\b(write\.(?:csv|table)|write_csv|saveRDS|save|to_csv)\b"
)

PREFIX_MAPS: dict[str, str] = {}


def iter_code_files() -> list[Path]:
    files: list[Path] = []
    for path in REPO_ROOT.rglob("*"):
        if ".git" in path.parts:
            continue
        if path == Path(__file__).resolve():
            continue
        if path.is_file() and path.suffix in CODE_SUFFIXES:
            files.append(path)
    return sorted(files)


def expand_userish(path: str) -> Path:
    if path.startswith("~/"):
        return Path.home() / path[2:]
    return Path(path)


def mapped_repo_path(path: str) -> Path | None:
    expanded = str(expand_userish(path))
    for old_prefix, repo_prefix in sorted(PREFIX_MAPS.items(), key=lambda item: len(item[0]), reverse=True):
        if expanded == old_prefix or expanded.startswith(old_prefix + "/"):
            suffix = expanded[len(old_prefix) :].lstrip("/")
            return REPO_ROOT / repo_prefix / suffix
    return None


def classify(path: str, line: str = "") -> tuple[str, str]:
    mapped = mapped_repo_path(path)
    if mapped is not None:
        if mapped.exists():
            return "repo_mapped_exists", str(mapped.relative_to(REPO_ROOT))
        external = expand_userish(path)
        if external.exists():
            return "repo_mapped_missing_external_exists", str(external)
        if OUTPUT_CONTEXT_RE.search(line):
            return "generated_output_missing", str(mapped.relative_to(REPO_ROOT))
        return "repo_mapped_missing", str(mapped.relative_to(REPO_ROOT))

    external = expand_userish(path)
    if external.exists():
        return "external_exists", str(external)
    if OUTPUT_CONTEXT_RE.search(line):
        return "generated_output_missing", str(external)
    return "external_missing_or_unmapped", str(external)


def audit() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for code_file in iter_code_files():
        rel = code_file.relative_to(REPO_ROOT)
        for line_no, line in enumerate(code_file.read_text(errors="replace").splitlines(), start=1):
            for match in LOCAL_PATH_RE.finditer(line):
                literal = match.group("path")
                status, resolved = classify(literal, line)
                rows.append(
                    {
                        "file": str(rel),
                        "line": str(line_no),
                        "status": status,
                        "literal": literal,
                        "resolved": resolved,
                    }
                )
    return rows


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    fields = ["file", "line", "status", "literal", "resolved"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", type=Path, help="Optional TSV output path.")
    parser.add_argument("--strict", action="store_true", help="Exit nonzero for missing/unmapped external paths.")
    args = parser.parse_args()

    rows = audit()
    if args.write:
        write_tsv(args.write, rows)
    else:
        fields = ["file", "line", "status", "literal", "resolved"]
        print("\t".join(fields))
        for row in rows:
            print("\t".join(row[field] for field in fields))

    missing = [row for row in rows if row["status"] in {"repo_mapped_missing", "external_missing_or_unmapped"}]
    if args.strict and missing:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
