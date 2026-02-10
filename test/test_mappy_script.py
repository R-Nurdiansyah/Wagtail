#!/usr/bin/env python3
"""Unit tests for bin/mappy_script.py"""

import os
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

import polars as pl

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')] + sys.path

mock_mappy = types.ModuleType('mappy')
mock_mappy.Aligner = None
mock_mappy.fastx_read = None
sys.modules.setdefault('mappy', mock_mappy)

from bin.mappy_script import process_alignment

DATA_DIR = Path(__file__).resolve().parents[1] / 'test' / 'test_data'


class FakeHit:
    def __init__(self, contig, mapq=0, primary=True):
        self.ctg = contig
        self.mapq = mapq
        self.is_primary = primary


class FakeAligner:
    def __init__(self, hits_map):
        self._hits_map = hits_map

    def map(self, sequence):
        return self._hits_map.get(sequence, [])


def _read_fastq_records(path):
    records = []
    with open(path, 'r', encoding='utf-8') as handle:
        while True:
            header = handle.readline().strip()
            if not header:
                break
            sequence = handle.readline().strip()
            handle.readline()  # plus line
            quality = handle.readline().strip()
            records.append((header.lstrip('@'), sequence, quality))
    return records


class TestMappyScript(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_process_alignment_writes_alignment_and_metadata(self):
        alignment_output = os.path.join(self.temp_dir.name, 'alignment.tsv')
        metadata_output = os.path.join(self.temp_dir.name, 'metadata.tsv')
        records = _read_fastq_records(DATA_DIR / 'mock_reads.fastq')

        hits_map = {
            records[0][1]: [FakeHit('contigA', mapq=42)],
            records[1][1]: [FakeHit(None, mapq=0, primary=False)],
        }

        with mock.patch('bin.mappy_script.mp.Aligner', return_value=FakeAligner(hits_map)) as mock_aligner:
            with mock.patch('bin.mappy_script.mp.fastx_read', return_value=records):
                process_alignment(
                    str(DATA_DIR / 'mock_reads.fastq'),
                    str(DATA_DIR / 'mock_reference.mmi'),
                    alignment_output,
                    metadata_output,
                    'SAMPLE123',
                )

        mock_aligner.assert_called_once()
        alignment_df = pl.read_csv(alignment_output, separator='\t')
        self.assertIn('contig', alignment_df.columns)
        self.assertTrue(alignment_df['contig'].str.contains('contigA').any())

        metadata_df = pl.read_csv(metadata_output, separator='\t')
        self.assertEqual(metadata_df['mapped reads'][0], 1)
        self.assertEqual(metadata_df['unmapped reads'][0], 1)
        self.assertAlmostEqual(metadata_df['mapping percentage'][0], 50.0)
        self.assertAlmostEqual(metadata_df['average mapq'][0], 42.0)


if __name__ == '__main__':
    unittest.main()
