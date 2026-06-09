#!/usr/bin/env python3
"""
wagtail_archive.py — bundle a finished Wagtail run's stage folders into .zip
archives. Called from rule cleanup (before the intermediate folders are
deleted).

The run directory holds the numbered stage folders:
    0_logs_wagtail  1_marker_id  1_manifest  2_qc_wagtail  3_deblur_wagtail
    4_rep_seqs_wagtail  4_table_wagtail  5_taxonomy_wagtail  6_condensed_wagtail
    7_metadata
Folders are selected by their leading "<digit>_" prefix (any "*_tmp" scratch
folder is always skipped).

Options (any combination):
    --archive-result             zip the 0_, 6_, 7_ folders
                                   → <run-name>_final_results.zip
    --archive-intermediate SPEC  zip selected intermediates (SPEC = "all" or a
                                   comma list of 1-5, e.g. 1,3,5)
                                   → <run-name>_intermediate.zip
    --archive-all                zip every 0_..7_ folder
                                   → <run-name>_all_results.zip

Archives are written at the top level of --run-dir, so they are never recursed
into (they don't start with a "<digit>_" prefix and the running .zip is skipped
explicitly).
"""

from __future__ import annotations

import argparse
import os
import sys
import zipfile
from pathlib import Path

RESULT_DIGITS = {"0", "6", "7"}
ALL_DIGITS    = {"0", "1", "2", "3", "4", "5", "6", "7"}


def parse_args():
    p = argparse.ArgumentParser(description="Archive Wagtail run stage folders into .zip files.")
    p.add_argument("--run-dir",  required=True, help="Directory holding the numbered stage folders.")
    p.add_argument("--run-name", required=True, help="Run name; prefixes the output archive filenames.")
    p.add_argument("--archive-result", action="store_true",
                   help="Archive the 0_, 6_, 7_ folders into <run-name>_final_results.zip")
    p.add_argument("--archive-intermediate", metavar="SPEC", default=None,
                   help="Archive intermediates into <run-name>_intermediate.zip. "
                        "SPEC = 'all' or a comma list of 1-5 (e.g. 1,3,5).")
    p.add_argument("--archive-all", action="store_true",
                   help="Archive every 0_..7_ folder into <run-name>_all_results.zip")
    return p.parse_args()


def intermediate_digits(spec: str) -> set[str]:
    """Resolve an --archive-intermediate SPEC to a set of digit strings (1-5)."""
    spec = spec.strip().lower()
    if spec == "all":
        return {"1", "2", "3", "4", "5"}
    digits = set()
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if not part.isdigit() or not 1 <= int(part) <= 5:
            raise ValueError(f"intermediate selection {part!r} out of range (1-5)")
        digits.add(part)
    if not digits:
        raise ValueError("no valid intermediate stages selected")
    return digits


def dirs_for_digits(run_dir: Path, digits: set[str]) -> list[Path]:
    """Return stage folders whose name starts with '<digit>_' for digit in *digits*.
    Always skips '*_tmp' scratch folders."""
    out = []
    for d in sorted(run_dir.iterdir()):
        if not d.is_dir():
            continue
        if d.name.endswith("_tmp"):
            continue
        if len(d.name) >= 2 and d.name[0] in digits and d.name[1] == "_":
            out.append(d)
    return out


def zip_dirs(run_dir: Path, dirs: list[Path], zip_path: Path) -> None:
    """Zip *dirs* (recursively) into *zip_path*, storing paths relative to run_dir."""
    if not dirs:
        print(f"  [skip] no folders to archive for {zip_path.name}")
        return
    n_files = 0
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED, allowZip64=True) as zf:
        for d in dirs:
            for root, _, files in os.walk(d):
                for fn in files:
                    fp = Path(root) / fn
                    if fp.resolve() == zip_path.resolve():
                        continue   # never archive the archive being written
                    zf.write(fp, fp.relative_to(run_dir))
                    n_files += 1
    folders = ", ".join(d.name for d in dirs)
    print(f"  [ok] {zip_path.name}: {n_files} file(s) from {len(dirs)} folder(s) [{folders}]")


def main() -> int:
    args = parse_args()
    run_dir = Path(args.run_dir)
    if not run_dir.is_dir():
        print(f"--run-dir not found: {run_dir}", file=sys.stderr)
        return 1

    if not (args.archive_result or args.archive_intermediate or args.archive_all):
        print("Nothing to do: pass --archive-result, --archive-intermediate, and/or --archive-all")
        return 0

    print(f"Archiving run '{args.run_name}' in {run_dir}")

    if args.archive_all:
        zip_dirs(run_dir,
                 dirs_for_digits(run_dir, ALL_DIGITS),
                 run_dir / f"{args.run_name}_all_results.zip")

    if args.archive_result:
        zip_dirs(run_dir,
                 dirs_for_digits(run_dir, RESULT_DIGITS),
                 run_dir / f"{args.run_name}_final_results.zip")

    if args.archive_intermediate:
        try:
            digits = intermediate_digits(args.archive_intermediate)
        except ValueError as exc:
            print(f"--archive-intermediate: {exc}", file=sys.stderr)
            return 1
        zip_dirs(run_dir,
                 dirs_for_digits(run_dir, digits),
                 run_dir / f"{args.run_name}_intermediate.zip")

    return 0


if __name__ == "__main__":
    sys.exit(main())
