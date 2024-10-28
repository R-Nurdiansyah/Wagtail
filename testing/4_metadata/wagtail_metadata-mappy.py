#script to create running metadata for wagtail
#running metadata-> combination of qc-stats, deblur-stats, and mappy metadata

import argparse
import logging
import sys
import subprocess
import shutil
import qiime2
import samtools
import polars as pl
from pathlib import Path
import os

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

#run export for qc-stats
def qc_export(input_qc, output_qc, sample_name):
    # Create the output directory
    os.makedirs(output_qc, exist_ok=True)

    # Qiime export command
    qiime_command = f"qiime tools export --input-path {input_qc} --output-path {output_qc}"

    # Run the Qiime command
    subprocess.run(qiime_command, shell=True, check=True)

    #manipulate the csv data inside the output directory
    for item in os.listdir(output_qc):
        if item.endswith('.csv'):
            #read the csv file using polars
            df = pl.read_csv(os.path.join(output_qc, item))
            #sum the number on column total-input-reads and total-retained-reads
            df = df.with_columns((pl.col("total-retained-reads") / pl.col("total-input-reads") * 100).alias("retained-reads-percentage"))
            #save the result to a csv file with sample_name + _qc_stats.csv in output_qc directory
            df.write_csv(os.path.join(output_qc, sample_name, '_qc_stats.csv'))
    
    qc_file = [f for f in os.listdir(output_qc) if f.endswith("_qc_stats.csv")]
    return qc_file

#run export for deblur-stats
def deblur_stat_export(input_deblur, output_deblur, sample_name):
    # Create the output directory
    os.makedirs(output_deblur, exist_ok=True)

    # Qiime export command
    qiime_command = f"qiime tools export --input-path {input_deblur} --output-path {output_deblur}"

    # Run the Qiime command
    subprocess.run(qiime_command, shell=True, check=True)

    #manipulate the csv data inside the output directory
    for item in os.listdir(output_deblur):
        if item.endswith('.csv'):
            #read the csv file using polars
            df = pl.read_csv(os.path.join(output_deblur, item))
            #save the result to a csv file with sample_name + _deblur_stats.csv in output_qc directory
            df.write_csv(os.path.join(output_deblur, sample_name, '_deblur_stats.csv'))
    
    deblur_file = [f for f in os.listdir(output_deblur) if f.endswith("_deblur_stats.csv")]
    return deblur_file

def combine_metadata(qc_file, deblur_file, mappy_file, output_path, run_name):
    #read the qc, deblur, mappy files in their respective output directory
    qc_files = pl.read_csv(qc_file)
    deblur_files = pl.read_csv(deblur_file)
    mappy_files = pl.read_csv(mappy_file)
    #combine the metadata into one file based on the sample_name
    metadata = qc_files.join(deblur_files, on='sample_name').join(mappy_files, on='sample_name')

    # append the metadata to the created metadata file from def empty_metadata
    metadata.write_csv(os.path.join(output_path, run_name+'_metadata.txt'), separator='\t')
    
if __name__ == "__main__":
    parent_parser = argparse.ArgumentParser()
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--input-qc', help='directory of the qiime qc-stat input (filter-stats.qza)', required=True)
    parent_parser.add_argument('--output-qc', help='directory of the qiime qc-stat output', required=True)
    parent_parser.add_argument('--input-deblur', help='directory of the qiime deblur-stat input (deblur-stats.qza)', required=True)
    parent_parser.add_argument('--output-deblur', help='directory of the qiime deblur-stat output', required=True)
    parent_parser.add_argument('--input-mappy', help='mappy metadata input file', required=True)
    parent_parser.add_argument('--output-path', help='directory of combined metadata output', required=True)
    parent_parser.add_argument('--sample-name', help='name of the sample', required=True)

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
    qc_file = qc_export(args.input_qc, args.output_qc, args.sample_name)
    deblur_file = deblur_stat_export(args.input_deblur, args.output_deblur, args.sample_name)
    combine_metadata(qc_file, deblur_file, args.input_mappy, args.output_path, args.sample_name)