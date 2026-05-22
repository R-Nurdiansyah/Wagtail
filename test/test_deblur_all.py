#!/usr/bin/env python3
"""
Unit tests for bin/deblur_all.py

Covers:
  - calculate_avg_read_length()  – average read length from gzipped FASTQ
  - deblur()                     – correct QIIME2 subcommand dispatched per marker
                                   (denoise-16S vs denoise-other, presence/absence
                                   of --i-reference-seqs, trim-length arithmetic)
"""

import gzip
import os
import sys
import tempfile
import unittest
from unittest.mock import call, patch

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), "..")] + sys.path

from bin.deblur_all import calculate_avg_read_length, deblur


# ── Shared helpers ────────────────────────────────────────────────────────────

def _write_fastq_gz(path, seqs):
    """Write a gzipped FASTQ file from a list of sequence strings."""
    with gzip.open(path, "wt") as fh:
        for i, seq in enumerate(seqs):
            qual = "I" * len(seq)
            fh.write(f"@read{i}\n{seq}\n+\n{qual}\n")


def _run_deblur(marker, avg_read_length=150, mock_run=None):
    """Helper: call deblur() with fixed dummy paths and return the command string."""
    deblur(
        input="input.qza",
        avg_read_length=avg_read_length,
        threads=4,
        representative="rep.qza",
        table="table.qza",
        stats="stats.qza",
        reference="ref.qza",
        marker=marker,
    )
    assert mock_run.called
    return mock_run.call_args[0][0]   # the shell command string


# ── calculate_avg_read_length ─────────────────────────────────────────────────

class TestCalculateAvgReadLength(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def _path(self, name="reads.fastq.gz"):
        return os.path.join(self.tmp.name, name)

    def test_uniform_length(self):
        path = self._path()
        _write_fastq_gz(path, ["ACGT" * 37] * 10)   # 148 bp each
        self.assertAlmostEqual(calculate_avg_read_length(path), 148.0)

    def test_mixed_lengths(self):
        path = self._path()
        _write_fastq_gz(path, ["A" * 100, "A" * 200])
        self.assertAlmostEqual(calculate_avg_read_length(path), 150.0)

    def test_single_read(self):
        path = self._path()
        _write_fastq_gz(path, ["C" * 250])
        self.assertAlmostEqual(calculate_avg_read_length(path), 250.0)

    def test_empty_file_raises_value_error(self):
        path = self._path()
        with gzip.open(path, "wt") as fh:
            pass   # empty file
        with self.assertRaises(ValueError):
            calculate_avg_read_length(path)


# ── deblur dispatch ───────────────────────────────────────────────────────────

class TestDeblurDispatch(unittest.TestCase):
    """
    deblur() must call the correct QIIME2 sub-command for each marker:
      16S  → denoise-16S   (no --i-reference-seqs)
      18S  → denoise-other (includes --i-reference-seqs)
      ITS  → denoise-other (includes --i-reference-seqs)
      CO1  → denoise-other (includes --i-reference-seqs)
    """

    @patch("bin.deblur_all.subprocess.run")
    def test_16s_uses_denoise_16s(self, mock_run):
        cmd = _run_deblur("16S", mock_run=mock_run)
        self.assertIn("denoise-16S", cmd)

    @patch("bin.deblur_all.subprocess.run")
    def test_16s_has_no_reference_seqs_flag(self, mock_run):
        """denoise-16S uses QIIME2's built-in reference; --i-reference-seqs must be absent."""
        cmd = _run_deblur("16S", mock_run=mock_run)
        self.assertNotIn("--i-reference-seqs", cmd)

    @patch("bin.deblur_all.subprocess.run")
    def test_18s_uses_denoise_other(self, mock_run):
        cmd = _run_deblur("18S", mock_run=mock_run)
        self.assertIn("denoise-other", cmd)

    @patch("bin.deblur_all.subprocess.run")
    def test_18s_includes_reference_seqs(self, mock_run):
        cmd = _run_deblur("18S", mock_run=mock_run)
        self.assertIn("--i-reference-seqs", cmd)
        self.assertIn("ref.qza", cmd)

    @patch("bin.deblur_all.subprocess.run")
    def test_its_uses_denoise_other(self, mock_run):
        cmd = _run_deblur("ITS", mock_run=mock_run)
        self.assertIn("denoise-other", cmd)

    @patch("bin.deblur_all.subprocess.run")
    def test_its_includes_reference_seqs(self, mock_run):
        cmd = _run_deblur("ITS", mock_run=mock_run)
        self.assertIn("--i-reference-seqs", cmd)
        self.assertIn("ref.qza", cmd)

    @patch("bin.deblur_all.subprocess.run")
    def test_co1_uses_denoise_other(self, mock_run):
        cmd = _run_deblur("CO1", mock_run=mock_run)
        self.assertIn("denoise-other", cmd)

    @patch("bin.deblur_all.subprocess.run")
    def test_co1_includes_reference_seqs(self, mock_run):
        cmd = _run_deblur("CO1", mock_run=mock_run)
        self.assertIn("--i-reference-seqs", cmd)

    # ── trim-length arithmetic ─────────────────────────────────────────────────

    @patch("bin.deblur_all.subprocess.run")
    def test_trim_length_is_avg_minus_10(self, mock_run):
        """trim_length must equal avg_read_length - 10."""
        cmd = _run_deblur("16S", avg_read_length=160, mock_run=mock_run)
        self.assertIn("--p-trim-length 150", cmd)

    @patch("bin.deblur_all.subprocess.run")
    def test_trim_length_for_non_16s(self, mock_run):
        cmd = _run_deblur("ITS", avg_read_length=200, mock_run=mock_run)
        self.assertIn("--p-trim-length 190", cmd)

    # ── required QIIME2 flags present ─────────────────────────────────────────

    @patch("bin.deblur_all.subprocess.run")
    def test_16s_required_flags(self, mock_run):
        cmd = _run_deblur("16S", mock_run=mock_run)
        for flag in [
            "--i-demultiplexed-seqs",
            "--p-trim-length",
            "--p-sample-stats",
            "--p-min-reads",
            "--p-jobs-to-start",
            "--o-representative-sequences",
            "--o-table",
            "--o-stats",
        ]:
            self.assertIn(flag, cmd, f"Required flag '{flag}' missing from denoise-16S command")

    @patch("bin.deblur_all.subprocess.run")
    def test_other_required_flags(self, mock_run):
        cmd = _run_deblur("ITS", mock_run=mock_run)
        for flag in [
            "--i-demultiplexed-seqs",
            "--i-reference-seqs",
            "--p-trim-length",
            "--p-sample-stats",
            "--p-min-reads",
            "--p-jobs-to-start",
            "--o-representative-sequences",
            "--o-table",
            "--o-stats",
        ]:
            self.assertIn(flag, cmd, f"Required flag '{flag}' missing from denoise-other command")

    # ── subprocess called exactly once ────────────────────────────────────────

    @patch("bin.deblur_all.subprocess.run")
    def test_subprocess_called_once(self, mock_run):
        _run_deblur("16S", mock_run=mock_run)
        self.assertEqual(mock_run.call_count, 1)

    @patch("bin.deblur_all.subprocess.run")
    def test_subprocess_check_true(self, mock_run):
        """subprocess.run must be called with check=True so errors propagate."""
        _run_deblur("16S", mock_run=mock_run)
        _, kwargs = mock_run.call_args
        self.assertTrue(kwargs.get("check", False))


if __name__ == "__main__":
    unittest.main()
