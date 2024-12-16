#script to run the exporter and deblur in one go
#usage: python read_seqs.py -i <artifact.qza> -o <output_dir> -t <threads> -r <representative_seqs.qza> -a <table.qza> -s <stats.qza>

#!/usr/bin/env python3
###############################################################################
#
#    Copyright (C) 2024 Rizky Nurdiansyah
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

__author__ = "Rizky Nurdiansyah"
__copyright__ = "Copyright 2024"
__credits__ = ["Rizky Nurdiansyah"]
__license__ = "GPL3"
__maintainer__ = "Rizky Nurdiansyah"
__email__ = "nurdiansyah.rizky.bio near gmail.com"
__status__ = "Development"

import gzip
import argparse
import logging
import sys
import os
import subprocess

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

def run_qiime_export(input, output_path):
    # Create the output directory
    os.makedirs(output_path, exist_ok=True)

    # Qiime export command
    qiime_command = f"qiime tools export --input-path {input} --output-path {output_path}"

    # Run the Qiime command
    subprocess.run(qiime_command, shell=True, check=True)
    #check any file with .fastq.gz extension in the output_path
    #if found, return the file name as input variable for calculate_avg_read_length
    #if not found, raise ValueError
    if any(fname.endswith('.fastq.gz') for fname in os.listdir(output_path)):
        file = [fname for fname in os.listdir(output_path) if fname.endswith('.fastq.gz')][0]
        #create the full path of the file
        extracted = os.path.join(output_path, file)
        return extracted
    else:
        raise ValueError("No FASTQ file found in the output directory")

def calculate_avg_read_length(extracted):
    total_length = 0
    read_count = 0
    line_counter = 0
    current_block = []

    def process_block(block):
        # validate and process the block
        if len(block) >= 2:
            seq = block[1].strip()
            return len(seq)
        return 0

    with gzip.open(extracted, 'rt') as f:
        for line in f:
            line_counter += 1
            current_block.append(line.strip())

            # If read 4 lines or the end of the block, process the block
            if len(current_block) == 4:
                total_length += process_block(current_block)
                read_count += 1
                current_block = []

        # process the last block
        if current_block:
            total_length += process_block(current_block)
            read_count += 1

    if read_count == 0:
        raise ValueError("Empty FASTQ file")

    avg_read_length = total_length / read_count
    logging.info(f"Average read length: {avg_read_length:.2f} basepairs")
    return avg_read_length

def deblur(input, avg_read_length, threads, representative, table, stats):
# Create the output directory
    #calculate trim length and make it integer
    trim_length = int(avg_read_length - 10)
    # Qiime export command
    qiime_command = f"qiime deblur denoise-16S --i-demultiplexed-seqs {input} --p-trim-length {trim_length} --p-sample-stats --p-min-reads 0 --p-jobs-to-start {threads} --o-representative-sequences {representative} --o-table {table} --o-stats {stats}"

    # Run the Qiime command
    subprocess.run(qiime_command, shell=True, check=True)

if __name__ == "__main__":
    parent_parser = argparse.ArgumentParser()
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--input', '-i', help='data input', required=True)
    parent_parser.add_argument('--output-path', '-o', help='directory of fastq.gz output', required=True)
    parent_parser.add_argument('--threads', '-t', help='threads number for multithreading', required=True)
    parent_parser.add_argument('--representative', '-r', help='output: representative seqs from deblur', required=True)
    parent_parser.add_argument('--table', '-a', help='output: frequency table from deblur', required=True)
    parent_parser.add_argument('--stats', '-s', help='output: deblur stats for metadata', required=True)
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
    extracted = run_qiime_export(args.input, args.output_path)
    avg_read_length = calculate_avg_read_length(extracted)
    deblur(args.input, avg_read_length, args.threads, args.representative, args.table, args.stats)
