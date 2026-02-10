#!/usr/bin/env python3
"""Unit tests for bin/extract_taxonomy.py"""

import os
import sys
import tempfile
import unittest
from pathlib import Path

import polars as pl

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')] + sys.path

from bin.extract_taxonomy import extract_taxonomy

DATA_DIR = Path(__file__).resolve().parents[1] / 'test' / 'test_data'


class TestExtractTaxonomy(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_extract_taxonomy_with_reference_database(self):
        output_path = os.path.join(self.temp_dir.name, 'gg2_taxonomy.tsv')
        extract_taxonomy(
            str(DATA_DIR / 'gg2_primary.tsv'),
            str(DATA_DIR / 'gg2_table.tsv'),
            str(DATA_DIR / 'gg2_reference.tsv'),
            'SAMPLE_A',
            output_path,
        )

        result = pl.read_csv(output_path, separator='\t')
        self.assertIn('taxonomy', result.columns)
        self.assertIn('coverage', result.columns)
        self.assertTrue(all(value.startswith('Root; ') for value in result['taxonomy']))
        self.assertTrue((result['coverage'] > 0).all())
        self.assertTrue((result['sample'] == 'SAMPLE_A').all())

    def test_extract_taxonomy_without_reference_database(self):
        output_path = os.path.join(self.temp_dir.name, 'danica_taxonomy.tsv')
        extract_taxonomy(
            str(DATA_DIR / 'danica_primary.tsv'),
            str(DATA_DIR / 'danica_table.tsv'),
            None,
            'DANICA_SAMPLE',
            output_path,
        )

        result = pl.read_csv(output_path, separator='\t')
        self.assertIn('taxonomy', result.columns)
        self.assertTrue((result['sample'] == 'DANICA_SAMPLE').all())
        self.assertTrue(result['taxonomy'].str.starts_with('Root').all())


if __name__ == '__main__':
    unittest.main()
