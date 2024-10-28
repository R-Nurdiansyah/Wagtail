#script to run mappy for both single-end as part of Wagtail pipeline

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
import os
import polars as pl
import mappy as mp
import time
from collections import Counter

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

#run the mappy for both paired and single-end
def mappy_run(index, read1, read2, single, output):
    #load or build index
    ref = mp.Aligner(index, preset = 'sr')
    if not ref: raise Exception("ERROR: failed to load/build index")
    #print to console if index is loaded from the file
    logging.info("Index loaded from the file {}".format(index))
   
    #give error if the output is not a directory
    #if not os.path.isdir(output):
    #    logging.error("Output directory is not exist")
    #    sys.exit(1)

    # check the number of sample reads in the input to ensure how much subdirectory is needed
    #if single:
    #    num_reads = len(read1)
    #else:
    #    num_reads = len(read1) + len(read2)
    
    
    #create the output directory (user's defined) if it does not exist
    if not os.path.exists(output):
        os.makedirs(output)
    
    #create subdirectory based on the input filename and make sure the output is saved in the subdirectory 
    #and the subdirectory do not exceed 1000 directories in one parent directory
    #if it is more than 1000 then make a new parent directory
    #parse the output to the lowest subdirectory
    #if num_reads > 1000:
    #    parent_dir = os.path.join(output, "output_{}".format(num_reads // 1000))
    #    if not os.path.exists(parent_dir):
    #        os.makedirs(parent_dir)
    #    output = parent_dir
    #    logging.info("Output is saved in {}".format(output))
    #else:
    #    output = output
    #    logging.info("Output is saved in {}".format(output))

    start_time = time.time()

    #check the input, if --single is used, then it will run the single-end alignment
    if single:
        #print to console if single-end alignment is used
        logging.info("Running single-end alignment")
        for i in read1:
            logging.info("Running analysis for {}".format(i))
            result_per_i = {}
            #Loop through single-end reads
            for name, seq, qual in mp.fastx_read(i):
                barcode = os.path.basename(i).split('.')[0]
                names = "{}_barcode+{}".format(barcode, name)
            #Perform single-end alignment
                for hits in ref.map(seq, cs=True, MD=True):
                    if hits.is_primary:
                        value = Counter(hits.ctg)
                        result_per_i[hits.ctg] = (value)
            #record the result per iteration and export it as tsv to output directory
            df = pl.from_dict(result_per_i)
            df = df.transpose(include_header=True)
            df.columns = ['read_name', 'ctg']
            df = df.with_columns([pl.col("read_name").str.splitn("_barcode+", 2).struct.rename_fields(["barcode", "read"]).alias("taxonomy")]).unnest("taxonomy")
            df = df.drop("read_name")
            df = df.select(["barcode", "read", "ctg"])
            #export to csv with barcode as the filename
            df.write_csv("{}/{}.tsv".format(output, barcode), separator='\t')
    
    else:
        #print to console if paired-end alignment is used
        logging.info("Running paired-end alignment")
        for i1, i2 in zip(read1, read2):
            logging.info("Running analysis for {} + {}".format(i1, i2))
            result_per_i = {}
            #Loop through paired-end reads
            for (name_f, seq_f, qual_f), (name_r, seq_r, qual_r) in zip(mp.fastx_read(i1), mp.fastx_read(i2)):
                #take the basename of the input file for barcode 
                barcode_f = os.path.basename(i1).split('.')[0]
                barcode_r = os.path.basename(i2).split('.')[0]
                names = "{}_paired_{}_barcode+{}".format(barcode_f, barcode_r, name_f)
            #Perform paired-end alignment
                for hits in ref.map(seq_f, seq_r, cs=True, MD=True):
                    if hits.is_primary:
                        result_per_i[names] = (hits.ctg)
            #record the result per iteration and export it as tsv to output directory
            df = pl.from_dict(result_per_i)
            df = df.transpose(include_header=True)
            df.columns = ['read_name', 'ctg']
            df = df.with_columns([pl.col("read_name").str.splitn("_barcode+", 2).struct.rename_fields(["barcode", "read"]).alias("taxonomy")]).unnest("taxonomy")
            df = df.drop("read_name")
            df = df.select(["barcode", "read", "ctg"])
            #export to csv with barcode as the filename
            df.write_csv("{}/{}+{}.tsv".format(output, barcode_f, barcode_r), separator='\t')
    
    #print to console if the analysis is done and time elapsed
    logging.info("Analysis done")
    elapsed_time = int(time.time() - start_time)
    #print elapsed time in h:mm:ss
    logging.info("Elapsed time: {:02d}:{:02d}:{:02d}".format(elapsed_time // 3600, (elapsed_time % 3600 // 60), elapsed_time % 60))
    logging.info("Output is saved in {}".format(output))


if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(description='Run mappy for both paired and single-end')
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")
    parent_parser.add_argument('-i', '--index', help='Path to the index file', required=True)
    #add nargs='+' to accept multiple input files for both r1 and r2
    parent_parser.add_argument('-1', '--read1', help='Path to the first/forward read file(s) in fasta/q or in fastq.gz format', required=True, nargs='+')
    parent_parser.add_argument('-2', '--read2', help='Path to the second/reverse read file(s) in fasta/q or in fastq.gz format, if paired-end is used', nargs='+')
    parent_parser.add_argument('-o', '--output', help='Path to the output directory', required=True)
    parent_parser.add_argument('-s', '--single', help='Run single-end alignment only', action="store_true")

    args = parent_parser.parse_args()
    
    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')

    mappy_run(args.index, args.read1, args.read2, args.single, args.output)
