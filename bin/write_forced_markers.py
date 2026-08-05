#!/usr/bin/env python3
"""
write_forced_markers.py — write a trivial OK _marker.json for every sample.

Used by rule marker_id's forced-amplicon fast path (wagtail.py --amplicon
16S/18S/ITS/CO1): no alignment is needed, so each sample gets a confident OK
marker pointing at the forced database family (reference = DB_<marker>).

Usage:
    python write_forced_markers.py --sample-list FILE --output-dir DIR --marker 16S

Writes <output-dir>/<sample>_marker.json for each non-blank line in the sample
list. The JSON keys mirror those written by identify_amplicon.py.
"""

import argparse
import json
import sys
from pathlib import Path


def forced_record(marker: str) -> dict:
    """A confident OK marker record for forced-amplicon mode."""
    return {
        "status":                "OK",
        "predicted_marker":      marker,
        "best_combined_score":   "1.0000",
        "best_mean_identity":    "1.0000",
        "best_hit_fraction":     "1.000",
        "second_marker":         "None",
        "second_combined_score": "NA",
        "second_mean_identity":  "NA",
        "second_hit_fraction":   "NA",
        "score_margin":          "1.0000",
        "reference":             f"DB_{marker}",
    }


def main() -> int:
    p = argparse.ArgumentParser(description="Write forced-amplicon marker JSONs.")
    p.add_argument("--sample-list", required=True, help="One sample/accession per line.")
    p.add_argument("--output-dir",  required=True, help="Directory to write <sample>_marker.json.")
    p.add_argument("--marker",      required=True, help="Forced marker: 16S | 18S | ITS | CO1.")
    args = p.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    text = json.dumps(forced_record(args.marker), indent=2)

    n = 0
    with open(args.sample_list) as fh:
        for line in fh:
            sample = line.strip()
            if not sample:
                continue
            (out_dir / f"{sample}_marker.json").write_text(text)
            n += 1

    print(f"Wrote {n} forced '{args.marker}' marker JSON file(s) to {out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
