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
database_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)),'../database')
data_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)),'test_data')

from bin.mappy_script import process_alignment

class TestMappyAlignment(unittest.TestCase):
    
    def setUp(self):
        # Set up temporary directory for output
        self.temp_dir = tempfile.TemporaryDirectory()
        self.output_dir = self.temp_dir.name
        # Set up metadata output in temporary directory with appropriate suffix
        self.meta_output = os.path.join(self.output_dir, 'test_sample_metadata.tsv') 
        self.alignment_output = os.path.join(self.output_dir, 'test_sample_alignment.tsv')
        
    def tearDown(self):
        self.temp_dir.cleanup()

    def test_mappy_alignment_successful(self):
        # Test if alignment is successful and alignment file and metadata files are created
        input = os.path.join(data_dir, 'dna-sequences.fasta')
        ref = os.path.join(database_dir, 'danica.mmi')
        sample_accession = 'test_sample'
        process_alignment(input, ref, self.alignment_output, self.meta_output, sample_accession)
        
        # Check if output files are created correctly
        expected_metadata_path = os.path.join(self.output_dir, 'test_sample_metadata.tsv')
        expected_alignment_path = os.path.join(self.output_dir, 'test_sample_alignment.tsv')
        self.assertTrue(os.path.exists(expected_metadata_path))
        self.assertTrue(os.path.exists(expected_alignment_path))
        
        # Verify content of alignment file
        alignment_df = pl.read_csv(expected_alignment_path, separator='\t')
        self.assertEqual(alignment_df['read'][0], 'fcb15987386cce65c60466dbc4d440b8')
        self.assertEqual(alignment_df['contig'][0], 'FLASV381465.1373;tax=d:Bacteria,p:Desulfobacterota,c:Desulfovibrionia,o:Desulfovibrionales,f:Desulfovibrionaceae,g:MFD_g_381465,s:MFD_s_381465;')
        #check each row to have non-empty values
        for i in range(len(alignment_df.columns)):
            self.assertTrue(alignment_df[alignment_df.columns[i]].is_not_null().all())

        # Verify content of metadata file
        metadata_df = pl.read_csv(expected_metadata_path, separator='\t')
        self.assertEqual(metadata_df['sample-id'][0], sample_accession)
        #check each row to have non-empty values
        for i in range(len(metadata_df.columns)):
            self.assertTrue(metadata_df[metadata_df.columns[i]].is_not_null().all())

if __name__ == '__main__':
    unittest.main()