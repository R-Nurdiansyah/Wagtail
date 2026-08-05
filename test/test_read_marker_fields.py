#!/usr/bin/env python3
"""
Unit tests for bin/read_marker_fields.py (v1.2 marker_log helper).

Invoked as a CLI script by the marker_log rule, so the tests run it as a
subprocess and check the six fields it prints (one per line).
"""

import json
import os
import subprocess
import sys
import tempfile
import unittest

_ROOT   = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..")
_SCRIPT = os.path.join(_ROOT, "bin", "read_marker_fields.py")


def _run(path):
    out = subprocess.run([sys.executable, _SCRIPT, path],
                         capture_output=True, text=True)
    return out.returncode, out.stdout.splitlines()


class TestReadMarkerFields(unittest.TestCase):

    def test_six_fields_from_valid_json(self):
        record = {
            "status": "OK",
            "predicted_marker": "16S",
            "best_mean_identity": "0.9900",
            "second_marker": "18S",
            "second_mean_identity": "0.3000",
            "score_margin": "0.6900",
        }
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
            json.dump(record, fh)
            path = fh.name
        try:
            rc, lines = _run(path)
        finally:
            os.unlink(path)
        self.assertEqual(rc, 0)
        self.assertEqual(lines, ["OK", "16S", "0.9900", "18S", "0.3000", "0.6900"])

    def test_missing_keys_use_defaults(self):
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
            json.dump({"status": "OK", "predicted_marker": "ITS"}, fh)
            path = fh.name
        try:
            _, lines = _run(path)
        finally:
            os.unlink(path)
        # status, predicted, best_mean_identity, second_marker, second_mean_identity, score_margin
        self.assertEqual(lines, ["OK", "ITS", "", "None", "NA", ""])

    def test_missing_file_yields_error_placeholders(self):
        rc, lines = _run("/nonexistent/marker.json")
        self.assertEqual(rc, 0)   # never hard-fails the marker_log rule
        self.assertEqual(lines, ["ERROR", "", "", "None", "NA", ""])


if __name__ == "__main__":
    unittest.main()
