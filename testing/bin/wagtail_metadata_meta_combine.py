#script to create running metadata for wagtail
#running metadata-> combination of qc-stats, deblur-stats, and mappy metadata

import argparse
import logging
import sys
import subprocess
import shutil
import polars as pl
from pathlib import Path
import os

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

#edit metadata
def combine_metadata(qc_file, deblur_file, mappy_file, output_path, run_name):
    #read the qc, deblur, mappy files in their respective output directory
    qc = pl.read_csv(qc_file)
    deblur = pl.read_csv(deblur_file)
    mappy = pl.read_csv(mappy_file, separator='\t', has_header=True)
    #sum the number on column total-input-reads and total-retained-reads
    qc = qc.with_columns((pl.col("total-retained-reads") / pl.col("total-input-reads") * 100).alias("retained-reads-percentage"))
    #combine the metadata into one file based on the sample_name
    metadata = qc.join(deblur, on='sample-id').join(mappy, on='sample-id')

    # append the metadata to the created metadata file from def empty_metadata
    metadata.write_csv(f"{output_path}/{run_name}_metadata.tsv", separator='\t')
    
if __name__ == "__main__":
    parent_parser = argparse.ArgumentParser()
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--qc-file', help='directory of the qiime qc-stat input (filter-stats.qza)', required=True)
    parent_parser.add_argument('--deblur-file', help='directory of the qiime deblur-stat input (deblur-stats.qza)', required=True)
    parent_parser.add_argument('--mappy-file', help='mappy metadata input file', required=True)
    parent_parser.add_argument('--output-path', help='directory of combined metadata output', required=True)
    parent_parser.add_argument('--run-name', help='name of the run', required=True)

    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')

    # Run the Qiime export command
    combine_metadata(args.qc_file, args.deblur_file, args.mappy_file, args.output_path, args.run_name)