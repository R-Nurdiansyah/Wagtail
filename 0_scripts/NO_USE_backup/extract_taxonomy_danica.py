#script to extract taxonomy information from Greengenes2 database and match it with the hit from the previous result (primary_silva.tsv)
#usable for KMA software only
#!/usr/bin/env python3

###############################################################################
#
#    Copyright (C) 2023 Rizky Nurdiansyah
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
__copyright__ = "Copyright 2023"
__credits__ = ["Rizky Nurdiansyah"]
__license__ = "GPL3"
__maintainer__ = "Rizky Nurdiansyah"
__email__ = "nurdiansyah.rizky.bio near gmail.com"
__status__ = "Development"

import polars as pl
import argparse
import logging
import os
import sys

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

#GREENGENES2 = 'greengenes2'

def extract_taxonomy(input_primaryID, sample_name, output_taxonomy_list):
    #1. read the input from the previous result (primary_[database].tsv) based on input
    #infer schema length set to 0 to ensure the read table is read as string -> easier that way
    df = pl.read_csv(input_primaryID, separator=';', has_header=False, infer_schema_length=0)
    
    #2. delete column_3 and rename column_1 to ID and column_2 to Taxonomy
    df = df.drop('column_3')
    df = pl.DataFrame.rename(df, {"column_1": "ID", "column_2": "taxonomy"})
    
    #delete the string "tax=" from the Taxonomy column
    df = df.with_columns(pl.col("taxonomy").str.replace("tax", ""))
    
    #replace the d: to d__ and p: to p__ in the Taxonomy column
    df = df.with_columns(pl.col("taxonomy").str.replace_all("=d:", "; d__")) #just to make sure we don't replace the wrong thing
    df = df.with_columns(pl.col("taxonomy").str.replace_all(",p:", "; p__"))
    df = df.with_columns(pl.col("taxonomy").str.replace_all(",c:", "; c__"))
    df = df.with_columns(pl.col("taxonomy").str.replace_all(",o:", "; o__"))
    df = df.with_columns(pl.col("taxonomy").str.replace_all(",f:", "; f__"))
    df = df.with_columns(pl.col("taxonomy").str.replace_all(",g:", "; g__"))
    df = df.with_columns(pl.col("taxonomy").str.replace_all(",s:", "; s__"))
    
    #count the duplicates and keep the ID and Taxonomy columns. put the count in a new column called count and sort by count
    df = df.groupby("taxonomy").agg(pl.count("taxonomy").alias("coverage"))
    df = df.sort("coverage", descending=True)

    #add "Root" at the start of each data in taxonomy column
    df = df.with_columns(pl.col("taxonomy").str.replace("^", "Root"))

    #add sample column on the first column and fill the value with user's input (sample name)
    df = df.with_columns(pl.Series("sample", [sample_name]*len(df)))

    #reorder the column to match the biobox script input while delete the unnecessary columns
    df = df[["sample", "coverage", "taxonomy"]]

     #6. write the result to a new file
    df.write_csv(output_taxonomy_list, separator='\t', has_header=True)

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--input-primaryID', help='input primaryID path', required=True)
    #parent_parser.add_argument('--input-taxonomy-reference', help='input taxonomy database path', required=True)
    parent_parser.add_argument('--sample-name', help='sample name for this data', required=True)
    parent_parser.add_argument('--output-taxonomy-list', help='output taxonomy list path in csv/tsv', required=True)
    #parent_parser.add_argument('--db', help='database name', required=True, choices=[SILVA, GREENGENES2])

    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    extract_taxonomy(args.input_primaryID, args.sample_name, args.output_taxonomy_list)