#script to extract taxonomy from combination of table and rep-seqs
#!/usr/bin/env python3

###############################################################################
#
#    Copyright (C) 2025 Rizky Nurdiansyah
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
__copyright__ = "Copyright 2025"
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

def extract_taxonomy(input_primaryID, input_table, input_taxonomy_reference, sample_name, output_taxonomy_list):
    #1a. read the input from the previous results (primary_.tsv)
    df = pl.read_csv(input_primaryID, separator='\t', has_header=False)
    
    #1b. read the table input
    df_table = pl.read_csv(input_table, separator='\t', has_header=False)

    #2. combine the primaryID and table
    joined = df.join(df_table, on="column_1", how="inner")

    #3. check any duplicates in the primaryID column
    duplicate_counts = joined["column_1"].value_counts()

    #4. merge the duplicate count with the joined table
    joined = joined.join(duplicate_counts, on="column_1", how="left")

    #5. adjust the column_2 based on the duplicate count -> divide by the count and round up
    joined = joined.with_columns([pl.col("column_2_right") / pl.col("count")])
    joined = joined.with_columns([pl.col("column_2_right").round()])

    #6. check if there is any duplicate in column_2_right and if any, aggregate the value in column_2
    joined = joined.group_by("column_2").agg(pl.sum("column_2_right").alias("column_2_right"))

    #7. edit the data to condensed format
    joined = pl.DataFrame.rename(joined, {"column_2_right": "coverage", "column_2": "taxonomy"})

    #split here based on the input reference
    #input reference is provided -> greengenes2
    #input reference is not provided -> danica
    if input_taxonomy_reference:
        #read the GG2 database
        ref = pl.read_csv(input_taxonomy_reference, separator='\t')
        logging.info(f"Read input {len(ref)} taxonomy entry from the Database")
        col1 = ref.columns[0] #take the name of the first column for merging
        #Assume greenegenes2 need to manipulate the ID a little bit
        joined = joined.with_columns(pl.col("taxonomy").str.replace("_.*", ""))
        final = joined.join(ref, left_on="taxonomy", right_on=col1, how="left") 
        #delete the unnecessary column after merger
        final = final.drop("taxonomy", "Confidence")
        #rename the Taxon column to taxonomy
        final = final.rename({"Taxon": "taxonomy"})
         # recheck for duplicates and then sort by coverage
        final = final.group_by("taxonomy").agg(pl.sum("coverage").alias("coverage"))
        final = final.sort("coverage", descending=True)
        #replace spaces with underscores in the taxonomy column
        final = final.with_columns(pl.col("taxonomy").str.replace_all(" ", "_").alias("taxonomy"))
        final = final.with_columns(pl.col("taxonomy").str.replace_all(r";_", "; ").alias("taxonomy"))

        #add "Root" at the start of each data in taxonomy column
        final = final.with_columns(pl.col("taxonomy").str.replace("^", "Root; "))

    else:
        #take only the second list after exploding taxonomy column by ';'
        final = joined.with_columns(pl.col("taxonomy").str.split(by=";").list.get(1).alias("taxonomy"))
        
        #delete the string "tax=" from the Taxonomy column
        final = final.with_columns(pl.col("taxonomy").str.replace("tax", ""))
        
        #replace the d: to d__ and p: to p__ in the Taxonomy column
        #just to make sure we don't replace the wrong thing
        final = final.with_columns(pl.col("taxonomy").str.replace_all("=d:", "; d__")) 
        final = final.with_columns(pl.col("taxonomy").str.replace_all(",p:", "; p__"))
        final = final.with_columns(pl.col("taxonomy").str.replace_all(",c:", "; c__"))
        final = final.with_columns(pl.col("taxonomy").str.replace_all(",o:", "; o__"))
        final = final.with_columns(pl.col("taxonomy").str.replace_all(",f:", "; f__"))
        final = final.with_columns(pl.col("taxonomy").str.replace_all(",g:", "; g__"))
        final = final.with_columns(pl.col("taxonomy").str.replace_all(",s:", "; s__"))
        
        # recheck for duplicates and then sort by coverage
        final = final.group_by("taxonomy").agg(pl.sum("coverage").alias("coverage"))
        final = final.sort("coverage", descending=True)

        #add "Root" at the start of each data in taxonomy column
        final = final.with_columns(pl.col("taxonomy").str.replace("^", "Root"))

    #add sample column on the first column and fill the value with user's input (sample name)
    final = final.with_columns(pl.Series("sample", [sample_name]*len(final)))

    #reorder the column to match the biobox script input while delete the unnecessary columns
    final = final[["sample", "coverage", "taxonomy"]]

    #8. write the result to a new file
    final.write_csv(output_taxonomy_list, separator='\t', include_header=True)

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")
    parent_parser.add_argument('--input-primaryID', '-i', help='input primaryID path', required=True)
    parent_parser.add_argument('--input-table','-t', help='input qiime2 table tsv output path', required=True)
    parent_parser.add_argument('--input-taxonomy-reference', '-r', help='input taxonomy database path, use this for GG2 database')
    parent_parser.add_argument('--sample-name', '-s', help='sample name for this data', required=True)
    parent_parser.add_argument('--output-taxonomy-list', '-o', help='output taxonomy list path in csv/tsv', required=True)
    
    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    extract_taxonomy(args.input_primaryID, args.input_table, args.input_taxonomy_reference, args.sample_name, args.output_taxonomy_list)