#script to create running metadata for wagtail
#running metadata-> combination of qc-stats, deblur-stats, and mappy metadata

import argparse
import logging
import sys
import subprocess
import qiime2
from pathlib import Path
import os

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

#run export for qc-stats
def qc_export(input_qc, output_qc):
    # Create the output directory
    os.makedirs(output_qc, exist_ok=True)

    # Qiime export command
    qiime_command = f"qiime tools export --input-path {input_qc} --output-path {output_qc}"

    # Run the Qiime command
    subprocess.run(qiime_command, shell=True, check=True)

#run export for deblur-stats
def deblur_stat_export(input_deblur, output_deblur):
    # Create the output directory
    os.makedirs(output_deblur, exist_ok=True)

    # Qiime export command
    qiime_command = f"qiime tools export --input-path {input_deblur} --output-path {output_deblur}"

    # Run the Qiime command
    subprocess.run(qiime_command, shell=True, check=True)
    
if __name__ == "__main__":
    parent_parser = argparse.ArgumentParser()
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--input-qc', help='directory of the qiime qc-stat input (filter-stats.qza)', required=True)
    parent_parser.add_argument('--output-qc', help='directory of the qiime qc-stat output', required=True)
    parent_parser.add_argument('--input-deblur', help='directory of the qiime deblur-stat input (deblur-stats.qza)', required=True)
    parent_parser.add_argument('--output-deblur', help='directory of the qiime deblur-stat output', required=True)
    
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
    qc_file = qc_export(args.input_qc, args.output_qc)
    deblur_file = deblur_stat_export(args.input_deblur, args.output_deblur)