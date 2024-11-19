#!/usr/bin/env python3
#script to test mappy part of wagtail

import unittest
import os
import sys
import polars as pl
import logging
import tempfile

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

#setting the directory path for database and test data
data_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)),'test_data')

from bin.extract_taxonomy_danica import extract_taxonomy

class TestExtractTaxon(unittest.TestCase):
    
    def setUp(self):
        # Set up temporary directory for output
        self.temp_dir = tempfile.TemporaryDirectory()
        self.output_dir = self.temp_dir.name
        # Set up metadata output in temporary directory with '_condensed' suffix
        self.condensed_output = os.path.join(self.output_dir, 'test_sample_condensed.tsv') 
        
    def tearDown(self):
        self.temp_dir.cleanup()

    def test_mappy_alignment_successful(self):
        # Test if alignment is successful and alignment file and metadata files are created
        input_align = os.path.join(data_dir, 'test_alignment.tsv')
        input_table = os.path.join(data_dir, 'test_table.tsv')
        sample_accession = 'test_sample'
        extract_taxonomy(input_align, input_table, sample_accession, self.condensed_output)
        
        # Check if manifest file is created correctly
        expected_output_path = os.path.join(self.output_dir, 'test_sample_condensed.tsv')
        self.assertTrue(os.path.exists(expected_output_path))
        
        # Verify content of condensed file
        condensed_df = pl.read_csv(expected_output_path, separator='\t')
        self.assertEqual(condensed_df['sample'][0], sample_accession)
        #make sure the content of coverage is float
        self.assertTrue(isinstance(condensed_df['coverage'][0], float))
        #make sure the content of taxonomy have 'Root;'
        self.assertTrue(condensed_df['taxonomy'][0].startswith('Root; '))
        #check each row to have non-empty values
        for i in range(len(condensed_df.columns)):
            self.assertTrue(condensed_df[condensed_df.columns[i]].is_not_null().all())

if __name__ == '__main__':
    unittest.main()