#!/usr/bin/env python3
"""Unit tests for qiime2 export helper scripts."""

import os
import shutil
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')] + sys.path

# Provide lightweight stand-ins for optional dependencies so the import succeeds.
sys.modules.setdefault('qiime2', types.ModuleType('qiime2'))
sys.modules.setdefault('biom', types.ModuleType('biom'))

from bin.qiime2_biom_export import run_qiime_export as run_biom_export
from bin.qiime2_seqs_export import run_qiime_export as run_seqs_export

DATA_DIR = Path(__file__).resolve().parents[1] / 'test' / 'test_data'


class TestQiimeExports(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_biom_export_copies_feature_table(self):
        output_dir = os.path.join(self.temp_dir.name, 'biom_export')
        source_biom = DATA_DIR / 'feature-table.biom'
        exported_biom = os.path.join(output_dir, 'renamed_feature_table.biom')

        def _fake_run(*_, **__):
            os.makedirs(output_dir, exist_ok=True)
            shutil.copyfile(source_biom, os.path.join(output_dir, 'feature-table.biom'))

        with mock.patch('bin.qiime2_biom_export.subprocess.run', side_effect=_fake_run) as mocked:
            run_biom_export('input_table.qza', output_dir, 'renamed_feature_table.biom')

        mocked.assert_called_once()
        self.assertTrue(os.path.exists(exported_biom))
        with open(source_biom, 'r', encoding='utf-8') as src, open(exported_biom, 'r', encoding='utf-8') as dest:
            self.assertEqual(src.read(), dest.read())

    def test_seqs_export_runs_command_and_creates_directory(self):
        output_dir = os.path.join(self.temp_dir.name, 'seqs_export')
        with mock.patch('bin.qiime2_seqs_export.subprocess.run') as mocked:
            run_seqs_export('rep_seqs.qza', output_dir)

        mocked.assert_called_once()
        self.assertTrue(os.path.isdir(output_dir))


if __name__ == '__main__':
    unittest.main()
