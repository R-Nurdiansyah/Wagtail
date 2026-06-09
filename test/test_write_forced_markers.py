#!/usr/bin/env python3
"""
Unit tests for bin/write_forced_markers.py (v1.2 forced-amplicon fast path).

Covers the forced_record() schema and the CLI that writes one _marker.json per
sample listed in --sample-list.
"""

import json
import os
import subprocess
import sys
import tempfile
import unittest

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), "..")] + sys.path

from bin.write_forced_markers import forced_record

_ROOT   = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..")
_SCRIPT = os.path.join(_ROOT, "bin", "write_forced_markers.py")

# Must match identify_amplicon.write_marker_file's schema so downstream readers
# (read_marker_fields.py, resolve_reference.py) work on forced and auto runs alike.
_SCHEMA_KEYS = (
    "status", "predicted_marker", "best_combined_score", "best_mean_identity",
    "best_hit_fraction", "second_marker", "second_combined_score",
    "second_mean_identity", "second_hit_fraction", "score_margin", "reference",
)


class TestForcedRecord(unittest.TestCase):

    def test_status_ok_and_reference_set(self):
        r = forced_record("18S")
        self.assertEqual(r["status"], "OK")
        self.assertEqual(r["predicted_marker"], "18S")
        self.assertEqual(r["reference"], "DB_18S")

    def test_has_full_marker_schema(self):
        r = forced_record("16S")
        for key in _SCHEMA_KEYS:
            self.assertIn(key, r)


class TestWriteForcedMarkersCli(unittest.TestCase):

    def test_writes_one_json_per_nonblank_sample(self):
        with tempfile.TemporaryDirectory() as d:
            sample_list = os.path.join(d, "samples.txt")
            with open(sample_list, "w") as fh:
                fh.write("S1\nS2\n\nS3\n")   # blank line must be ignored
            out_dir = os.path.join(d, "markers")

            subprocess.run(
                [sys.executable, _SCRIPT,
                 "--sample-list", sample_list,
                 "--output-dir",  out_dir,
                 "--marker",      "CO1"],
                check=True, capture_output=True, text=True,
            )

            self.assertEqual(
                sorted(os.listdir(out_dir)),
                ["S1_marker.json", "S2_marker.json", "S3_marker.json"],
            )
            with open(os.path.join(out_dir, "S1_marker.json")) as fh:
                rec = json.load(fh)
            self.assertEqual(rec["reference"], "DB_CO1")
            self.assertEqual(rec["status"], "OK")


if __name__ == "__main__":
    unittest.main()
