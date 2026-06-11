#!/usr/bin/env python3
"""
Unit tests for bin/marker_summary.py (v1.2 per-run marker summary TSV).

Covers auto_row() classification and the CLI (auto + forced modes).
"""

import json
import os
import subprocess
import sys
import tempfile
import unittest

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), "..")] + sys.path

from bin.marker_summary import auto_row

_ROOT   = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..")
_SCRIPT = os.path.join(_ROOT, "bin", "marker_summary.py")


def _write_json(directory, sample, data):
    with open(os.path.join(directory, f"{sample}_marker.json"), "w") as fh:
        json.dump(data, fh)


class TestAutoRow(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def test_ok_returns_predicted_marker(self):
        _write_json(self.tmp.name, "S1", {"status": "OK", "predicted_marker": "ITS"})
        self.assertEqual(auto_row(self.tmp.name, "S1"), ("OK", "ITS"))

    def test_ambiguous_is_fail_unknown(self):
        _write_json(self.tmp.name, "S2", {"status": "AMBIGUOUS", "predicted_marker": "16S"})
        self.assertEqual(auto_row(self.tmp.name, "S2"), ("FAIL", "UNKNOWN"))

    def test_unknown_is_fail_unknown(self):
        _write_json(self.tmp.name, "S3", {"status": "UNKNOWN", "predicted_marker": ""})
        self.assertEqual(auto_row(self.tmp.name, "S3"), ("FAIL", "UNKNOWN"))

    def test_missing_file_is_fail_unknown(self):
        self.assertEqual(auto_row(self.tmp.name, "nope"), ("FAIL", "UNKNOWN"))

    def test_empty_file_is_fail_unknown(self):
        open(os.path.join(self.tmp.name, "S4_marker.json"), "w").close()
        self.assertEqual(auto_row(self.tmp.name, "S4"), ("FAIL", "UNKNOWN"))


class TestCli(unittest.TestCase):

    def _run(self, marker_dir, sample_list, out, forced=None):
        cmd = [sys.executable, _SCRIPT,
               "--marker-dir", marker_dir,
               "--sample-list", sample_list,
               "--output", out]
        if forced:
            cmd += ["--forced", forced]
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        with open(out) as fh:
            return [line.rstrip("\n") for line in fh]

    def test_auto_mode_mixed(self):
        with tempfile.TemporaryDirectory() as d:
            mdir = os.path.join(d, "markers"); os.makedirs(mdir)
            _write_json(mdir, "S1", {"status": "OK", "predicted_marker": "16S"})
            _write_json(mdir, "S2", {"status": "AMBIGUOUS"})
            # S3 has no JSON → FAIL/UNKNOWN
            sl = os.path.join(d, "samples.txt")
            with open(sl, "w") as fh:
                fh.write("S1\nS2\nS3\n")
            out = os.path.join(d, "summary.tsv")
            lines = self._run(mdir, sl, out)
            self.assertEqual(lines[0], "sample\tstatus\tamplicon")
            self.assertEqual(lines[1], "S1\tOK\t16S")
            self.assertEqual(lines[2], "S2\tFAIL\tUNKNOWN")
            self.assertEqual(lines[3], "S3\tFAIL\tUNKNOWN")

    def test_forced_mode_overrides_everything(self):
        with tempfile.TemporaryDirectory() as d:
            mdir = os.path.join(d, "markers"); os.makedirs(mdir)
            # Even an OK-looking JSON must be reported as FORCED.
            _write_json(mdir, "S1", {"status": "OK", "predicted_marker": "18S"})
            sl = os.path.join(d, "samples.txt")
            with open(sl, "w") as fh:
                fh.write("S1\nS2\n")
            out = os.path.join(d, "summary.tsv")
            lines = self._run(mdir, sl, out, forced="CO1")
            self.assertEqual(lines[1], "S1\tFORCED\tCO1")
            self.assertEqual(lines[2], "S2\tFORCED\tCO1")

    def test_forced_all_is_treated_as_auto(self):
        with tempfile.TemporaryDirectory() as d:
            mdir = os.path.join(d, "markers"); os.makedirs(mdir)
            _write_json(mdir, "S1", {"status": "OK", "predicted_marker": "16S"})
            sl = os.path.join(d, "samples.txt")
            with open(sl, "w") as fh:
                fh.write("S1\n")
            out = os.path.join(d, "summary.tsv")
            lines = self._run(mdir, sl, out, forced="ALL")
            self.assertEqual(lines[1], "S1\tOK\t16S")


if __name__ == "__main__":
    unittest.main()
