#!/usr/bin/env python3
"""
Unit tests for bin/resolve_reference.py (v1.2 runtime reference resolver).

Covers:
  - read_marker()   – derive (database_marker, status) from a _marker.json
  - resolve_path()  – glob the database dir for {marker}_{kind}* (with the
                      16S-deblur-builtin and empty/non-OK special cases)
"""

import json
import os
import sys
import tempfile
import unittest

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), "..")] + sys.path

from bin.resolve_reference import read_marker, resolve_path


class TestReadMarker(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def _json(self, data):
        path = os.path.join(self.tmp.name, "m.json")
        with open(path, "w") as fh:
            json.dump(data, fh)
        return path

    def test_db_prefix_stripped(self):
        path = self._json({"reference": "DB_16S", "status": "OK"})
        self.assertEqual(read_marker(path), ("16S", "OK"))

    def test_empty_reference_gives_empty_marker(self):
        path = self._json({"reference": "", "status": "AMBIGUOUS"})
        self.assertEqual(read_marker(path), ("", "AMBIGUOUS"))

    def test_missing_reference_key_gives_empty_marker(self):
        # No 'reference' key at all → empty marker (sample soft-fails downstream).
        path = self._json({"status": "UNKNOWN", "predicted_marker": "ITS"})
        marker, status = read_marker(path)
        self.assertEqual(marker, "")
        self.assertEqual(status, "UNKNOWN")

    def test_reference_without_db_prefix_used_verbatim(self):
        path = self._json({"reference": "CO1", "status": "OK"})
        self.assertEqual(read_marker(path), ("CO1", "OK"))

    def test_missing_file_returns_error(self):
        self.assertEqual(read_marker("/nonexistent/m.json"), ("", "ERROR"))

    def test_corrupt_json_returns_error(self):
        path = os.path.join(self.tmp.name, "bad.json")
        with open(path, "w") as fh:
            fh.write("{not valid json")
        self.assertEqual(read_marker(path), ("", "ERROR"))


class TestResolvePath(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.db = self.tmp.name

    def tearDown(self):
        self.tmp.cleanup()

    def _touch(self, name):
        open(os.path.join(self.db, name), "w").close()

    def test_empty_marker_returns_empty(self):
        self._touch("16S_database_x.fasta")
        self.assertEqual(resolve_path("", self.db, "database"), "")

    def test_16s_deblur_is_builtin_returns_empty(self):
        # 16S denoise uses QIIME2's built-in reference → no file expected.
        self._touch("16S_deblur_x.qza")
        self.assertEqual(resolve_path("16S", self.db, "deblur"), "")

    def test_single_match_returned(self):
        self._touch("ITS_database_unite.fasta")
        res = resolve_path("ITS", self.db, "database")
        self.assertEqual(os.path.basename(res), "ITS_database_unite.fasta")

    def test_16s_database_still_resolves(self):
        # 16S only skips the *deblur* reference; database/taxonomy still resolve.
        self._touch("16S_database_mfd.fasta")
        res = resolve_path("16S", self.db, "database")
        self.assertEqual(os.path.basename(res), "16S_database_mfd.fasta")

    def test_no_match_returns_empty(self):
        self.assertEqual(resolve_path("CO1", self.db, "taxonomy"), "")

    def test_multiple_matches_returns_empty(self):
        self._touch("18S_database_a.fasta")
        self._touch("18S_database_b.fasta")
        self.assertEqual(resolve_path("18S", self.db, "database"), "")


if __name__ == "__main__":
    unittest.main()
