#!/usr/bin/env python3
"""Unit tests for bin/deblur_all.py"""

import gzip
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')] + sys.path

from bin.deblur_all import calculate_avg_read_length, deblur, run_qiime_export

DATA_DIR = Path(__file__).resolve().parents[1] / 'test' / 'test_data'


class TestDeblurAll(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_run_qiime_export_returns_fastq_path(self):
        output_dir = os.path.join(self.temp_dir.name, 'export')
        exported_fastq = os.path.join(output_dir, 'mock.fastq.gz')

        def _fake_run(*_, **__):
            os.makedirs(output_dir, exist_ok=True)
            with open(exported_fastq, 'wb') as handle:
                handle.write(b'FASTQ')

        with mock.patch('bin.deblur_all.subprocess.run', side_effect=_fake_run) as mocked:
            result = run_qiime_export('input.qza', output_dir)

        mocked.assert_called_once()
        self.assertEqual(result, exported_fastq)

    def test_run_qiime_export_raises_without_fastq(self):
        output_dir = os.path.join(self.temp_dir.name, 'empty_export')
        os.makedirs(output_dir, exist_ok=True)

        with mock.patch('bin.deblur_all.subprocess.run'):
            with self.assertRaises(ValueError):
                run_qiime_export('input.qza', output_dir)

    def test_calculate_avg_read_length_uses_fastq_content(self):
        gz_path = os.path.join(self.temp_dir.name, 'reads.fastq.gz')
        with open(DATA_DIR / 'mock_reads.fastq', 'r', encoding='utf-8') as src, gzip.open(gz_path, 'wt') as dest:
            dest.write(src.read())

        avg_length = calculate_avg_read_length(gz_path)
        # mock_reads.fastq contains reads of length 10 and 6
        self.assertAlmostEqual(avg_length, 8.0)

    def test_deblur_invokes_subprocess_with_trim_length(self):
        with mock.patch('bin.deblur_all.subprocess.run') as mocked:
            deblur('input.qza', 55.2, '4', 'rep.qza', 'table.qza', 'stats.qza')

        called_command = mocked.call_args.args[0]
        self.assertIn('--p-trim-length 45', called_command)
        self.assertIn('--p-jobs-to-start 4', called_command)


if __name__ == '__main__':
    unittest.main()
