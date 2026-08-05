#!/usr/bin/env python3
"""
Unit tests for bin/sql_merge.py — v1.2 _expand_db_args() directory globbing.

cleanup now passes the 0_tmp *directory* (not thousands of file paths) so the
command line stays under ARG_MAX at scale; sql_merge must expand a directory to
its *_log.sql files while still accepting explicit file paths.
"""

import os
import sys
import tempfile
import unittest

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), "..")] + sys.path

from bin.sql_merge import _expand_db_args


def _touch(path):
    open(path, "w").close()
    return path


class TestExpandDbArgs(unittest.TestCase):

    def test_directory_globs_only_log_sql(self):
        with tempfile.TemporaryDirectory() as d:
            _touch(os.path.join(d, "a_log.sql"))
            _touch(os.path.join(d, "b_log.sql"))
            _touch(os.path.join(d, "notes.txt"))      # must be ignored
            result = _expand_db_args([d])
            self.assertEqual([os.path.basename(p) for p in result],
                             ["a_log.sql", "b_log.sql"])

    def test_explicit_file_passes_through(self):
        with tempfile.TemporaryDirectory() as d:
            f = _touch(os.path.join(d, "x_log.sql"))
            self.assertEqual(_expand_db_args([f]), [f])

    def test_mixed_files_and_directory(self):
        with tempfile.TemporaryDirectory() as d:
            sub = os.path.join(d, "tmp")
            os.makedirs(sub)
            _touch(os.path.join(sub, "c_log.sql"))
            f = _touch(os.path.join(d, "x_log.sql"))
            result = _expand_db_args([f, sub])
            self.assertEqual([os.path.basename(p) for p in result],
                             ["x_log.sql", "c_log.sql"])

    def test_empty_directory_yields_nothing(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(_expand_db_args([d]), [])


if __name__ == "__main__":
    unittest.main()
