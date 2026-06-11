#!/usr/bin/env python3
"""
marker_summary.py — write a one-row-per-sample TSV summarising marker_id.

Columns (tab-separated, with header):
    sample    status    amplicon

  status:
    OK      — confident marker call (auto-detect mode)
    FAIL    — AMBIGUOUS / UNKNOWN / missing / unreadable marker (auto-detect)
    FORCED  — forced-amplicon mode (--forced); no identification was run
  amplicon:
    OK      → the predicted marker (16S/18S/ITS/CO1)
    FAIL    → UNKNOWN
    FORCED  → the forced marker

In auto-detect mode each sample's <marker-dir>/<sample>_marker.json is read.
In --forced mode every sample is reported FORCED/<marker> without reading JSONs.

Usage:
    marker_summary.py --marker-dir DIR --sample-list FILE --output TSV [--forced MARKER]
"""

import argparse
import json
import sys
from pathlib import Path


def auto_row(marker_dir: str, sample: str):
    """Return (status, amplicon) for one sample in auto-detect mode."""
    path = Path(marker_dir) / f"{sample}_marker.json"
    try:
        data = json.loads(path.read_text())
        if not isinstance(data, dict):
            data = {}
    except (OSError, ValueError):
        data = {}
    if data.get("status") == "OK":
        return "OK", (data.get("predicted_marker") or "UNKNOWN")
    return "FAIL", "UNKNOWN"


def main() -> int:
    p = argparse.ArgumentParser(description="Summarise marker_id results into a TSV.")
    p.add_argument("--marker-dir", required=True, help="Directory of <sample>_marker.json files.")
    p.add_argument("--sample-list", required=True, help="One sample/accession per line.")
    p.add_argument("--output", required=True, help="Output TSV path.")
    p.add_argument("--forced", default=None,
                   help="Forced marker. If set (and not 'ALL'), every row is FORCED/<marker>.")
    args = p.parse_args()

    forced = args.forced
    if forced and forced.upper() == "ALL":
        forced = None   # auto-detect

    with open(args.sample_list) as fh:
        samples = [line.strip() for line in fh if line.strip()]

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as fh:
        fh.write("sample\tstatus\tamplicon\n")
        for sample in samples:
            if forced:
                status, amplicon = "FORCED", forced
            else:
                status, amplicon = auto_row(args.marker_dir, sample)
            fh.write(f"{sample}\t{status}\t{amplicon}\n")

    print(f"Wrote marker summary for {len(samples)} sample(s) to {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
