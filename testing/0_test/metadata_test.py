#!/usr/bin/env python3
#script to test mappy part of wagtail

import unittest
import os
import sys
import polars as pl
import logging
import tempfile

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

#append the path to the the related data, database, and script in its respective directory
script_dir = os.path.abspath("../0_scripts")
database_dir = os.path.abspath("../0_database")
data_dir = os.path.abspath("test_data")

#append the script directory to the sys.path so it can be imported
if script_dir not in sys.path:
    sys.path.insert(0, script_dir)

from wagtail_metadata_meta_combine import combine_metadata

class TestMappyAlignment(unittest.TestCase):
    
    def setUp(self):
        # Set up temporary directory for output
        self.temp_dir = tempfile.TemporaryDirectory()
        self.output_dir = self.temp_dir.name 
        
    def tearDown(self):
        self.temp_dir.cleanup()

    def test_mappy_alignment_successful(self):
        # Test if alignment is successful and alignment file and metadata files are created
        input_qc = os.path.join(data_dir, 'test_qc.csv')
        input_deblur = os.path.join(data_dir, 'test_deblur.csv')
        input_mappy = os.path.join(data_dir, 'test_metadata.tsv')
        run_name = 'test_sample'
        combine_metadata(input_qc, input_deblur, input_mappy, self.output_dir, run_name)
        
        # Check if manifest file is created correctly
        expected_output_path = os.path.join(self.output_dir, 'test_sample_metadata.tsv')
        self.assertTrue(os.path.exists(expected_output_path))
        
        # Verify content of metadata file
        metadata_df = pl.read_csv(expected_output_path, separator='\t')
        #check if the output content has 25 columns
        self.assertEqual(len(metadata_df.columns), 25)
        #check if the first row of the output has the sample-id as the sample-id from all inputs
        #take the content of the first row of the first column of the input qc file
        qc_df = pl.read_csv(input_qc)
        sample_accession = qc_df['sample-id'][0]
        self.assertEqual(metadata_df['sample-id'][0], sample_accession)
        #check each row to have non-empty values
        for i in range(len(metadata_df.columns)):
            self.assertTrue(metadata_df[metadata_df.columns[i]].is_not_null().all())

if __name__ == '__main__':
    unittest.main()