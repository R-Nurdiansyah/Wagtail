#!/usr/bin/env python3
"""Unit tests for bin/error_annotation.py"""

import os
import subprocess
import sys
import unittest
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / 'test' / 'test_data'
SCRIPT_PATH = PROJECT_ROOT / 'bin' / 'error_annotation.py'


class TestErrorAnnotation(unittest.TestCase):
    def _run_script(self, log_file, rule):
        result = subprocess.run(
            [sys.executable, str(SCRIPT_PATH), str(DATA_DIR / log_file), rule],
            capture_output=True,
            text=True,
            check=True,
        )
        return result.stdout.strip()

    def test_mappy_polars_hint(self):
        output = self._run_script('log_mappy_polars.txt', 'mappy')
        self.assertIn('Sample not mapped to used database', output)

    def test_quality_control_filtered_out_hint(self):
        output = self._run_script('log_quality_control.txt', 'quality_control')
        self.assertIn('Sample may not 16S', output)

    def test_deblur_empty_sequence_hint(self):
        output = self._run_script('log_deblur.txt', 'deblur')
        self.assertIn('Sample may not 16S', output)

    def test_upstream_failure_silently_exits(self):
        output = self._run_script('log_upstream.txt', 'mappy')
        self.assertEqual(output, '')

    def test_generic_message_for_unknown_rule(self):
        output = self._run_script('log_generic.txt', 'unknown_rule')
        self.assertEqual(output, 'Check STDERR in log')


if __name__ == '__main__':
    unittest.main()
