#script to extract taxonomy from combination of table and rep-seqs
#!/usr/bin/env python3

###############################################################################
#
#    Copyright (C) 2026 Rizky Nurdiansyah
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

def extract_taxonomy(input_alignment, input_table, input_taxonomy_reference, sample_name, output_taxonomy_list):
    df  = pl.read_csv(input_alignment, separator='\t', has_header=False)
    df_table = pl.read_csv(input_table, separator='\t', has_header=False)

    joined = df.join(df_table, on="column_1", how="inner")

    hit_counts = joined["column_1"].value_counts()
    joined = joined.join(hit_counts, on="column_1", how="left")
    joined = joined.with_columns(
        (pl.col("column_2_right") / pl.col("count")).round().alias("coverage"))
    per_taxonomy = joined.group_by("column_2").agg(pl.sum("coverage")).rename({"column_2": "id"})

    if input_taxonomy_reference is not None:
        # External reference database (e.g. GG2): join contig IDs to taxonomy strings.
        # Use inner join so only contigs with a known taxonomy are retained.
        ref = pl.read_csv(input_taxonomy_reference, separator='\t', has_header=False)
        logging.info(f"Read {len(ref)} taxonomy entries from reference database")
        combined = (per_taxonomy
                    .join(ref, left_on="id", right_on="column_1", how="inner")
                    .group_by("column_2").agg(pl.sum("coverage").alias("coverage"))
                    .sort("coverage", descending=True)
                    .with_columns(pl.col("column_2").str.replace("^", "Root;"))
                    )
        final = combined.rename({"column_2": "taxonomy"})
    else:
        # No external reference: taxonomy is embedded in the contig names themselves.
        combined = (per_taxonomy
                    .sort("coverage", descending=True)
                    .with_columns(pl.col("id").str.replace("^", "Root; "))
                    )
        final = combined.rename({"id": "taxonomy"})

    final = final.with_columns(pl.Series("sample", [sample_name] * len(final)))
    final = final[["sample", "coverage", "taxonomy"]]
    final.write_csv(output_taxonomy_list, separator='\t', include_header=True)

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")
    parent_parser.add_argument('--input-alignment', '-i', help='input primaryID path', required=True)
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
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%d/%m/%Y %I:%M:%S %p')

    extract_taxonomy(args.input_alignment, args.input_table, args.input_taxonomy_reference, args.sample_name, args.output_taxonomy_list)