#script to export qc_stats from qiime2 deblur-stat artifact and then copy the stats.csv to a new directory
#usage: python export_biom.py <artifact.qza> <output_dir>

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

import argparse
import logging
import sys
import subprocess
import shutil
import qiime2
import biom
import os

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

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

if __name__ == "__main__":
    parent_parser = argparse.ArgumentParser()
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")
    parent_parser.add_argument('--input-path', help='location of the deblur stat artifact file', required=True)
    parent_parser.add_argument('--output-path', help='directory of output', required=True)
    parent_parser.add_argument('--sample-name', help='sample name', required=True)
    
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
    run_qiime_export(args.input_path, args.output_path, args.sample_name)
    
    print("Process completed.")