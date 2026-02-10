#!/usr/bin/env python3
"""Unit tests for bin/wagtail_metadata_qiimes.py"""

import os
import sys
import tempfile
import unittest
from unittest import mock

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')] + sys.path

from bin.wagtail_metadata_qiimes import deblur_stat_export, qc_export


class TestWagtailMetadataQiimes(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_qc_export_invokes_qiime_tools(self):
        output_dir = os.path.join(self.temp_dir.name, 'qc')
        with mock.patch('bin.wagtail_metadata_qiimes.subprocess.run') as mocked:
            qc_export('input_filter.qza', output_dir)

        mocked.assert_called_once()
        command = mocked.call_args.args[0]
        self.assertIn('qiime tools export', command)
        self.assertIn('input_filter.qza', command)
        self.assertTrue(os.path.isdir(output_dir))

    def test_deblur_export_invokes_qiime_tools(self):
        output_dir = os.path.join(self.temp_dir.name, 'deblur')
        with mock.patch('bin.wagtail_metadata_qiimes.subprocess.run') as mocked:
            deblur_stat_export('deblur.qza', output_dir)

        mocked.assert_called_once()
        command = mocked.call_args.args[0]
        self.assertIn('deblur.qza', command)
        self.assertTrue(os.path.isdir(output_dir))


if __name__ == '__main__':
    unittest.main()
