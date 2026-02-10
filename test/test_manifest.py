#!/usr/bin/env python3
#script to test manifest part of wagtail

import unittest
import os
import sys
import polars as pl
import logging
import tempfile

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

from bin.create_manifest_wagtail import create_manifest

class TestCreateManifest(unittest.TestCase):
    
    def setUp(self):
        # Set up temporary directory for manifest output and shared file map
        self.temp_dir = tempfile.TemporaryDirectory()
        self.output_dir = self.temp_dir.name
        self.data_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'test_data')
        self.file_map_path = os.path.join(self.data_dir, 'manifest_file_map.tsv')
        
    def tearDown(self):
        self.temp_dir.cleanup()

    def test_create_manifest_successful(self):
        # Test if manifest is created when accession exists in file map
        input_accession = 'ACC_TEST'
        create_manifest(input_accession, self.file_map_path, self.output_dir)
        
        # Check if manifest file is created correctly
        expected_manifest_path = os.path.join(self.output_dir, input_accession + "_manifest.csv")
        self.assertTrue(os.path.exists(expected_manifest_path))
        
        # Verify content of manifest file
        manifest_df = pl.read_csv(expected_manifest_path)
        self.assertEqual(manifest_df['sample-id'][0], input_accession)
        expected_path = os.path.abspath(os.path.join(self.data_dir, 'sample_forward.fastq'))
        self.assertEqual(manifest_df['absolute-filepath'][0], expected_path)
        self.assertEqual(manifest_df['direction'][0], "forward")

    def test_create_manifest_not_found(self):
        # Test if no manifest file is created when accession is not in file map
        input_accession = 'ACC_UNKNOWN'
        with self.assertRaises(SystemExit):
            create_manifest(input_accession, self.file_map_path, self.output_dir)

        expected_manifest_path = os.path.join(self.output_dir, input_accession + "_manifest.csv")
        self.assertFalse(os.path.exists(expected_manifest_path))

if __name__ == '__main__':
    unittest.main()