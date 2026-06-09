#!/usr/bin/env python3
"""
read_marker_fields.py — print selected fields from a sample's _marker.json,
one per line, for the marker_log rule.

Usage:
    python read_marker_fields.py <path/to/{sample}_marker.json>

Prints six lines in this fixed order:
    status
    predicted_marker
    best_mean_identity
    second_marker
    second_mean_identity
    score_margin

A missing, empty, or corrupt file yields the same placeholder set the marker_log
rule treats as a failure (status = ERROR), so the caller can stay simple.
"""

import json
import sys

# (key, default-when-absent) — order defines the printed line order.
FIELDS = [
    ("status",               "ERROR"),
    ("predicted_marker",     ""),
    ("best_mean_identity",   ""),
    ("second_marker",        "None"),
    ("second_mean_identity", "NA"),
    ("score_margin",         ""),
]


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: read_marker_fields.py <marker.json>", file=sys.stderr)
        return 2

    try:
        with open(sys.argv[1]) as fh:
            data = json.load(fh)
        if not isinstance(data, dict):
            data = {}
    except (OSError, ValueError):
        data = {}

    for key, default in FIELDS:
        print(data.get(key, default))
    return 0


if __name__ == "__main__":
    sys.exit(main())
