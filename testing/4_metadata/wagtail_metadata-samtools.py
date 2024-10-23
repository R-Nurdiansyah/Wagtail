#script to create running metadata for wagtail
#content: retained reads after qc + trimming, mapping percentage, number of reads and save it in a txt file

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
#function to create empty metadata in txt format
def empty_metadata(output_path, run_name):
    file_path = os.path.join(output_path, run_name+'_metadata.txt')
    if not file_path.exists():
        file_path.touch()
    else:
        logging.info(f'{file_path} already exists')

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

#run samtools stats to get the mapped reads by minimap2
def run_samtools_stats(input_sam, output_sam, sample_name):
    # Create the output directory
    os.makedirs(output_sam, exist_ok=True)

    # Samtools stats command
    output_name = os.path.join(output_sam, sample_name, "_samtools_stats.txt")
    samtools_command = f"samtools flagstats -O tsv {input_sam} > {output_name}"

    # Run the Samtools command
    subprocess.run(samtools_command, shell=True, check=True)

    #read the txt to polars
    df = pl.read_csv(output_name, separator= '\t', has_header=False)
    index_to_keep = [0, 1, 2, 3, 7, 8, 9]
    index_to_delete = [1]
    #take only the needed rows
    df = df.filter(pl.Series(range(len(df))).is_in(index_to_keep))
    #transpose and make the last column as the header
    df = df.select([pl.col("column_1"), pl.col("column_2")]).transpose()
    df.columns = df["column_3"].to_list()
    #delete the unnecessary row
    df = df.filter(~pl.Series(range(len(df))).is_in(index_to_delete))
    #add one more column with the column name as sample-id and the value is as sample_name
    df = df.with_columns(pl.Series([sample_name for i in range(len(df))]).alias("sample-id"))
    #save the result to a csv file with sample_name + _samtools_stats.csv in output_sam directory
    df.write_csv(os.path.join(output_sam, sample_name, '_samtools_edited_stats.csv'))

    sam_file = [f for f in os.listdir(output_sam) if f.endswith("_samtools_edited_stats.csv")]
    return sam_file

def combine_metadata(qc_file, deblur_file, sam_file, output_path, run_name):
    #read the qc, deblur, samtools files in their respective output directory
    qc_files = pl.read_csv(qc_file)
    deblur_files = pl.read_csv(deblur_file)
    sam_files = pl.read_csv(sam_file)
    #combine the metadata into one file based on the sample_name
    metadata = qc_files.join(deblur_files, on='sample_name').join(sam_files, on='sample_name')

    # append the metadata to the created metadata file from def empty_metadata
    metadata.write_csv(os.path.join(output_path, run_name+'_metadata.txt'), separator='\t')
    

if __name__ == "__main__":
    parent_parser = argparse.ArgumentParser()
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--input-qc', help='location of the qiime species artifact file', required=True)
    parent_parser.add_argument('--input-sam', help='location of the sam file', required=True)
    parent_parser.add_argument('--output-path', help='directory of biom output', required=True)
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
    run_qc_export(args.input_qc, args.output_path, args.sample_name)
    run_samtools_stats(args.input_sam, args.output_path, args.sample_name)