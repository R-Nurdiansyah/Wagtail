import mappy as mp
import polars as pl
import os
import argparse
import sys
import logging
from collections import Counter

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

def process_alignment(input, reference, alignment_output, sample_id):
    # Load reference database
    ref = mp.Aligner(reference, preset='sr')
    if not ref:
        logging.error(f"ERROR: failed to load/build index from {reference}")

    # Create empty lists to store alignment data and metadata counts
    total_mapq = 0
    alignment_counter = Counter()

    # Process each input sequence in the input file (assuming FASTA/FASTQ format)
    for i in input:
        result_per_i = {}
        for name, seq, qual in mp.fastx_read(i, read_comment=False):
            is_mapped = False  # To track duplicates
            # Perform alignment and store results
            for hits in ref.map(seq):
                if hits.is_primary:
                    result_per_i[name] = (hits.ctg)
                    alignment_counter['mapped'] += 1
                    total_mapq += hits.mapq
                    is_mapped = True
                    break # Only consider primary alignment
                else:
                    alignment_counter['unmapped'] += 1
        #count the number of hits using Counter
        ctg_counts = Counter(result_per_i.values())

        #record the ctg_counts dictionary to a tsv file using polars
        df = pl.from_dict(ctg_counts)
        df = df.transpose(include_header=True)
        #add sample column on the first column and fill the value with user's input (sample name)
        df.columns = ['contig', 'count']
        df = df.with_columns(pl.Series("sample", [sample_id]*len(df)))
        #reorder the column to match the biobox script input while delete the unnecessary columns
        df = df[["sample", "contig", "count"]]
        #sort based on count
        df = df.sort("count", descending=True)
        df.write_csv(f"{alignment_output}/{sample_id}_alignment.tsv", separator='\t')

        # Calculate metadata
        total_reads = alignment_counter["mapped"] + alignment_counter["unmapped"]
        mapping_percentage = (alignment_counter["mapped"] / total_reads) * 100 if total_reads > 0 else 0
        average_mapq = total_mapq / alignment_counter["mapped"] if alignment_counter["mapped"] > 0 else 0

    # Create metadata DataFrame using Polars
    metadata = {
    "Total Reads": total_reads,
    "Mapped Reads": alignment_counter["mapped"],
    "Unmapped Reads": alignment_counter["unmapped"],
    "Mapping Percentage": f"{mapping_percentage:.2f}%",
    "Average Mapping Quality": f"{average_mapq:.2f}"
    }

    # Return metadata to be used in another function
    return metadata
    
def save_metadata_to_txt(metadata, metadata_output, sample_id):
    filename = f"{metadata_output}/{sample_id}_metadata.txt"
    try:
        with open(filename, "w") as f:
            for key, value in metadata.items():
                f.write(f"{key}: {value}\n")
        logging.info(f"Metadata successfully saved to {filename}")
    except Exception as e:
        print(f"Error saving metadata: {e}")

if __name__ == "__main__":
    parent_parser = argparse.ArgumentParser(description='Run mappy for both paired and single-end')
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")
    parent_parser.add_argument('-r', '--reference', help='Path to the index/reference file', required=True)
    parent_parser.add_argument('-i', '--input', help='Path to the first/forward read file(s) in fasta/q or in fastq.gz format', required=True, nargs='+')
    parent_parser.add_argument('-a', '--alignment_output', help='Path to the alignment output directory', required=True)
    parent_parser.add_argument('-m', '--metadata_output', help='Path to the metadata output directory', required=True)
    parent_parser.add_argument('-s', '--sample_id', help='Unique sample ID', required=True)
    args = parent_parser.parse_args()
    
    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')
    
    # Process alignment and generate output
    process_alignment(args.input, args.reference, args.alignment_output, args.sample_id)
    metadata = process_alignment(args.input, args.reference, args.alignment_output, args.sample_id)
    save_metadata_to_txt(metadata, args.metadata_output, args.sample_id)
