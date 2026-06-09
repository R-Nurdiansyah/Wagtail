#!/usr/bin/env python3
"""
resolve_reference.py — resolve the reference database for one sample at RUNTIME.

This is the small extractor that replaces the old DAG-time Snakemake helpers
(get_deblur_db / get_mappy_ref / get_tax_ref).  Because marker_id is now a plain
rule (not a checkpoint), the per-sample {sample}_marker.json does not exist when
Snakemake builds the DAG — so reference resolution is deferred to the rule's
shell, which calls this script.

It reads the 'reference' key written by identify_amplicon.py:
    "DB_16S" / "DB_ITS" / …  →  database family "16S" / "ITS" / …
    ""  (empty)              →  no confident call; nothing to resolve (soft-fail)

and maps that family + a --kind to a single file under --db-dir:
    --kind deblur     →  {db_dir}/{marker}_deblur*     (16S uses QIIME2 built-in → empty)
    --kind database   →  {db_dir}/{marker}_database*
    --kind taxonomy   →  {db_dir}/{marker}_taxonomy*

Output (stdout):
    --what path   (default)  →  resolved path, or "" when none applies
    --what marker            →  database marker, or "16S" fallback
    --what both              →  "<marker>\\t<path>"

The script never hard-fails the calling rule: resolution problems (0 or >1
matches) are reported to stderr and yield an empty path, so the per-sample rule
fails/soft-fails for just that sample rather than aborting the whole batch.
Database presence is validated up front by `wagtail.py --check-db`.
"""

from __future__ import annotations

import argparse
import glob
import json
import sys
from pathlib import Path


def read_marker(marker_json: str) -> tuple[str, str]:
    """Return (database_marker, status) from a sample's _marker.json.

    database_marker derives from the 'reference' key (DB_16S → 16S); if that is
    absent it falls back to predicted_marker, else "" (treated as soft-fail).
    """
    try:
        data = json.loads(Path(marker_json).read_text())
    except (OSError, ValueError):
        return "", "ERROR"
    ref = (data.get("reference") or "").strip()
    if ref.startswith("DB_"):
        marker = ref[3:]
    elif ref:
        marker = ref
    else:
        marker = ""   # empty reference = non-OK sample → soft-fail downstream
    return marker, data.get("status", "")


def resolve_path(marker: str, db_dir: str, kind: str) -> str:
    """Resolve the single db_dir/{marker}_{kind}* file, or "" when N/A."""
    if not marker:
        return ""                       # non-OK sample: nothing to resolve
    if kind == "deblur" and marker == "16S":
        return ""                       # denoise-16S uses QIIME2's built-in reference
    matches = sorted(glob.glob(str(Path(db_dir) / f"{marker}_{kind}*")))
    if len(matches) == 1:
        return matches[0]
    if not matches:
        print(f"resolve_reference: no '{marker}_{kind}*' in {db_dir}", file=sys.stderr)
    else:
        print(f"resolve_reference: multiple '{marker}_{kind}*' in {db_dir}: {matches}",
              file=sys.stderr)
    return ""


def main() -> int:
    p = argparse.ArgumentParser(description="Resolve a sample's reference DB at runtime.")
    p.add_argument("--marker-json", required=True, help="Path to the sample's _marker.json")
    p.add_argument("--db-dir",      required=True, help="Database directory to search")
    p.add_argument("--kind", required=True, choices=["deblur", "database", "taxonomy"],
                   help="Which database family to resolve")
    p.add_argument("--what", default="path", choices=["path", "marker", "both"],
                   help="What to print (default: path)")
    args = p.parse_args()

    marker, _status = read_marker(args.marker_json)
    out_marker = marker or "16S"        # safe fallback for deblur's required --marker

    if args.what == "marker":
        print(out_marker)
    elif args.what == "both":
        print(f"{out_marker}\t{resolve_path(marker, args.db_dir, args.kind)}")
    else:
        print(resolve_path(marker, args.db_dir, args.kind))
    return 0


if __name__ == "__main__":
    sys.exit(main())
