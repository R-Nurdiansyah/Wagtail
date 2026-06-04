#!/usr/bin/env python3
"""
Unit tests for bin/identify_amplicon.py

Covers:
  - reads_from_filemap()   – TSV lookup helpers
  - subsample_reads()      – reservoir-sampling FASTQ reader
  - identify_marker()      – ranking, confidence gates, 16S priority rule
"""

import csv
import gzip
import os
import sys
import tempfile
import unittest

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), "..")] + sys.path

try:
    from bin.identify_amplicon import (
        AlignmentStats,
        identify_marker,
        reads_from_filemap,
        subsample_reads,
    )
    _MAPPY_AVAILABLE = True
except ImportError:
    _MAPPY_AVAILABLE = False


# ── Shared helpers ────────────────────────────────────────────────────────────

def _make_stats(marker, hit_fraction, mean_identity, n_reads=1000):
    """Construct an AlignmentStats with uniform per-hit identity values."""
    n_hits = int(hit_fraction * n_reads)
    return AlignmentStats(
        marker=marker,
        n_reads=n_reads,
        n_hits=n_hits,
        hit_fraction=hit_fraction,
        mean_identity=mean_identity,
        median_identity=mean_identity,
        identities=[mean_identity] * n_hits,
    )


def _write_fastq(fh, n_reads, seq_len=150):
    """Write *n_reads* FASTQ records to an open text file handle."""
    seq  = "A" * seq_len
    qual = "I" * seq_len
    for i in range(n_reads):
        fh.write(f"@read{i}\n{seq}\n+\n{qual}\n")


# ── reads_from_filemap ────────────────────────────────────────────────────────

@unittest.skipUnless(_MAPPY_AVAILABLE, "mappy not installed in this environment")
class TestReadsFromFilemap(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def _make_filemap(self, rows):
        path = os.path.join(self.tmp.name, "filemap.tsv")
        with open(path, "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            for row in rows:
                writer.writerow(row)
        return path

    def test_first_row_found(self):
        fmap = self._make_filemap([
            ["ACC001", "/data/ACC001.fastq.gz"],
            ["ACC002", "/data/ACC002.fastq.gz"],
        ])
        self.assertEqual(reads_from_filemap(fmap, "ACC001"), "/data/ACC001.fastq.gz")

    def test_non_first_row_found(self):
        fmap = self._make_filemap([
            ["ACC001", "/data/ACC001.fastq.gz"],
            ["ACC002", "/data/ACC002.fastq.gz"],
        ])
        self.assertEqual(reads_from_filemap(fmap, "ACC002"), "/data/ACC002.fastq.gz")

    def test_not_found_raises_value_error(self):
        fmap = self._make_filemap([["ACC001", "/data/ACC001.fastq.gz"]])
        with self.assertRaises(ValueError):
            reads_from_filemap(fmap, "ACC999")

    def test_strips_whitespace_in_accession_and_path(self):
        """Leading/trailing whitespace in the TSV cells should be ignored."""
        path = os.path.join(self.tmp.name, "filemap.tsv")
        with open(path, "w") as fh:
            fh.write("  ACC001  \t  /data/ACC001.fastq.gz  \n")
        self.assertEqual(reads_from_filemap(path, "ACC001"), "/data/ACC001.fastq.gz")


# ── subsample_reads ───────────────────────────────────────────────────────────

@unittest.skipUnless(_MAPPY_AVAILABLE, "mappy not installed in this environment")
class TestSubsampleReads(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def _make_fastq(self, n_reads, seq_len=150, gzipped=False):
        suffix = ".fastq.gz" if gzipped else ".fastq"
        path = os.path.join(self.tmp.name, f"reads{suffix}")
        opener = gzip.open if gzipped else open
        with opener(path, "wt") as fh:
            _write_fastq(fh, n_reads, seq_len)
        return path

    def test_returns_all_when_file_smaller_than_n(self):
        path = self._make_fastq(50)
        result = subsample_reads(path, 1000)
        self.assertEqual(len(result), 50)

    def test_returns_exactly_n_when_file_larger(self):
        path = self._make_fastq(5000)
        result = subsample_reads(path, 100)
        self.assertEqual(len(result), 100)

    def test_returns_all_when_file_equals_n(self):
        path = self._make_fastq(100)
        result = subsample_reads(path, 100)
        self.assertEqual(len(result), 100)

    def test_deterministic_with_same_seed(self):
        path = self._make_fastq(5000)
        r1 = subsample_reads(path, 200, seed=42)
        r2 = subsample_reads(path, 200, seed=42)
        self.assertEqual(r1, r2)

    def test_gzipped_fastq_reads_correctly(self):
        path = self._make_fastq(200, gzipped=True)
        result = subsample_reads(path, 50)
        self.assertEqual(len(result), 50)

    def test_missing_file_returns_empty_list(self):
        result = subsample_reads("/nonexistent/reads.fastq", 100)
        self.assertEqual(result, [])

    def test_returns_sequences_not_headers_or_quality(self):
        """Reservoir must return bare DNA sequences, not FASTQ header/quality lines."""
        path = self._make_fastq(20, seq_len=150)
        result = subsample_reads(path, 20)
        for seq in result:
            self.assertFalse(seq.startswith("@"), "Header line returned instead of sequence")
            self.assertFalse(seq.startswith("+"), "Quality separator returned instead of sequence")
            self.assertTrue(len(seq) > 0)


# ── identify_marker ───────────────────────────────────────────────────────────

@unittest.skipUnless(_MAPPY_AVAILABLE, "mappy not installed in this environment")
class TestIdentifyMarker(unittest.TestCase):
    """
    Tests for identify_marker() covering:
      - 16S priority rule (fires / does not fire / exact boundary)
      - UNKNOWN gate (best hit fraction below min_confidence)
      - AMBIGUOUS gate (identity margin below min_margin)
      - OK clear winner
      - Correct second-best reporting
    """

    CONF   = 0.50   # min_confidence
    MARGIN = 0.05   # min_margin
    FRAC   = 0.90   # min_16s_fraction

    def _call(self, stats_list, conf=CONF, margin=MARGIN, frac=FRAC):
        stats = {s.marker: s for s in stats_list}
        return identify_marker(stats, conf, margin, frac)

    # ── 16S priority rule ─────────────────────────────────────────────────────

    def test_16s_priority_fires_even_when_18s_has_higher_identity(self):
        """
        Regression: V6V8 sample where 18S edged 16S by 0.0007 mean identity.
        With 16S fraction = 0.996 (> 0.90), result must be OK + 16S.
        """
        stats = [
            _make_stats("16S", hit_fraction=0.996, mean_identity=0.9931),
            _make_stats("18S", hit_fraction=0.999, mean_identity=0.9938),
            _make_stats("ITS", hit_fraction=0.000, mean_identity=0.0),
            _make_stats("CO1", hit_fraction=0.000, mean_identity=0.0),
        ]
        status, predicted, best, _ = self._call(stats)
        self.assertEqual(status, "OK")
        self.assertEqual(predicted, "16S")
        self.assertEqual(best.marker, "16S")

    def test_16s_priority_at_exact_threshold(self):
        """hit_fraction == threshold must trigger the rule (>=, not >)."""
        stats = [
            _make_stats("16S", hit_fraction=0.90, mean_identity=0.99),
            _make_stats("18S", hit_fraction=0.91, mean_identity=0.995),
            _make_stats("ITS", hit_fraction=0.00, mean_identity=0.0),
            _make_stats("CO1", hit_fraction=0.00, mean_identity=0.0),
        ]
        status, predicted, _, _ = self._call(stats)
        self.assertEqual(status, "OK")
        self.assertEqual(predicted, "16S")

    def test_16s_priority_not_triggered_below_threshold(self):
        """16S fraction just below threshold must fall through to normal gates."""
        stats = [
            _make_stats("16S", hit_fraction=0.899, mean_identity=0.99),
            _make_stats("18S", hit_fraction=0.000, mean_identity=0.0),
            _make_stats("ITS", hit_fraction=0.000, mean_identity=0.0),
            _make_stats("CO1", hit_fraction=0.000, mean_identity=0.0),
        ]
        # Normal gates: 16S fraction 0.899 > 0.6 conf, and margin vs zero is large → OK
        status, predicted, _, _ = self._call(stats)
        self.assertEqual(status, "OK")
        self.assertEqual(predicted, "16S")

    def test_16s_priority_second_is_highest_non_16s(self):
        """When priority rule fires, second-best must be the top non-16S marker."""
        stats = [
            _make_stats("16S", hit_fraction=0.99, mean_identity=0.993),
            _make_stats("18S", hit_fraction=0.99, mean_identity=0.994),
            _make_stats("ITS", hit_fraction=0.50, mean_identity=0.950),
            _make_stats("CO1", hit_fraction=0.00, mean_identity=0.0),
        ]
        _, _, _, second = self._call(stats)
        self.assertIsNotNone(second)
        self.assertEqual(second.marker, "18S")   # highest non-16S by combined score

    # ── UNKNOWN gate ──────────────────────────────────────────────────────────

    def test_unknown_when_best_fraction_below_min_confidence(self):
        """No marker clears the confidence threshold → UNKNOWN."""
        stats = [
            _make_stats("16S", hit_fraction=0.50, mean_identity=0.95),
            _make_stats("18S", hit_fraction=0.10, mean_identity=0.80),
            _make_stats("ITS", hit_fraction=0.05, mean_identity=0.70),
            _make_stats("CO1", hit_fraction=0.00, mean_identity=0.0),
        ]
        status, _, _, _ = self._call(stats)
        self.assertEqual(status, "UNKNOWN")

    # ── AMBIGUOUS gate ────────────────────────────────────────────────────────

    def test_ambiguous_when_margin_too_small(self):
        """Two markers with near-identical combined scores → AMBIGUOUS.
        ITS: 0.90×0.9500=0.8550, 18S: 0.88×0.9490=0.8351 → margin 0.0199 < 0.05"""
        stats = [
            _make_stats("ITS", hit_fraction=0.90, mean_identity=0.9500),
            _make_stats("18S", hit_fraction=0.88, mean_identity=0.9490),
            _make_stats("16S", hit_fraction=0.05, mean_identity=0.80),
            _make_stats("CO1", hit_fraction=0.00, mean_identity=0.0),
        ]
        status, _, _, _ = self._call(stats)
        self.assertEqual(status, "AMBIGUOUS")

    def test_ok_when_margin_equals_threshold(self):
        """Combined-score margin exactly equal to min_margin is NOT ambiguous (>= threshold).
        ITS: 1.00×0.900=0.900, 18S: 1.00×0.850=0.850 → margin exactly 0.05"""
        stats = [
            _make_stats("ITS", hit_fraction=1.00, mean_identity=0.9000),
            _make_stats("18S", hit_fraction=1.00, mean_identity=0.8500),  # score margin = 0.05
            _make_stats("16S", hit_fraction=0.05, mean_identity=0.80),
            _make_stats("CO1", hit_fraction=0.00, mean_identity=0.0),
        ]
        status, predicted, _, _ = self._call(stats)
        self.assertEqual(status, "OK")
        self.assertEqual(predicted, "ITS")

    # ── OK clear winner ───────────────────────────────────────────────────────

    def test_ok_clear_16s_sample(self):
        stats = [
            _make_stats("16S", hit_fraction=0.99, mean_identity=0.993),
            _make_stats("18S", hit_fraction=0.00, mean_identity=0.0),
            _make_stats("ITS", hit_fraction=0.00, mean_identity=0.0),
            _make_stats("CO1", hit_fraction=0.00, mean_identity=0.0),
        ]
        status, predicted, _, _ = self._call(stats)
        self.assertEqual(status, "OK")
        self.assertEqual(predicted, "16S")

    def test_ok_clear_its_sample(self):
        stats = [
            _make_stats("ITS", hit_fraction=0.95, mean_identity=0.97),
            _make_stats("16S", hit_fraction=0.05, mean_identity=0.80),
            _make_stats("18S", hit_fraction=0.05, mean_identity=0.79),
            _make_stats("CO1", hit_fraction=0.00, mean_identity=0.0),
        ]
        status, predicted, _, _ = self._call(stats)
        self.assertEqual(status, "OK")
        self.assertEqual(predicted, "ITS")

    def test_ok_clear_co1_sample(self):
        stats = [
            _make_stats("CO1", hit_fraction=0.92, mean_identity=0.96),
            _make_stats("16S", hit_fraction=0.00, mean_identity=0.0),
            _make_stats("18S", hit_fraction=0.00, mean_identity=0.0),
            _make_stats("ITS", hit_fraction=0.00, mean_identity=0.0),
        ]
        status, predicted, _, _ = self._call(stats)
        self.assertEqual(status, "OK")
        self.assertEqual(predicted, "CO1")

    # ── Edge cases ────────────────────────────────────────────────────────────

    def test_second_is_none_with_single_marker_in_stats(self):
        """When only one marker was aligned (early-exit), second must be None."""
        stats = [_make_stats("16S", hit_fraction=0.99, mean_identity=0.99)]
        _, _, _, second = self._call(stats)
        self.assertIsNone(second)

    def test_custom_16s_fraction_threshold(self):
        """Raising min_16s_fraction above the actual 16S fraction must not trigger rule."""
        stats = [
            _make_stats("16S", hit_fraction=0.92, mean_identity=0.993),
            _make_stats("18S", hit_fraction=0.93, mean_identity=0.996),
            _make_stats("ITS", hit_fraction=0.00, mean_identity=0.0),
            _make_stats("CO1", hit_fraction=0.00, mean_identity=0.0),
        ]
        # With frac=0.95, 16S fraction 0.92 should NOT trigger the rule
        _, _, best, _ = self._call(stats, frac=0.95)
        # 18S has higher mean identity and rule didn't fire → 18S wins (AMBIGUOUS or OK)
        self.assertEqual(best.marker, "18S")


if __name__ == "__main__":
    unittest.main()
