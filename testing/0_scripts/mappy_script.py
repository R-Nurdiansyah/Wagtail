import polars as pl
import mappy as mp
import os
import argparse
import sys
import logging
from collections import Counter

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

def process_alignment(input, reference, alignment_output, metadata_output, sample_id):
    # Load reference database
    ref = mp.Aligner(reference, preset = 'sr') # short-read preset Illumina
    if not ref: raise Exception("ERROR: failed to load/build index")

    # Create empty lists to store alignment data and metadata counts
    total_mapq = 0
    alignment_counter = Counter()

    # Process each input sequence in the input file (assuming FASTA/FASTQ format)
    result_per_i = {}
    for name, seq, qual in mp.fastx_read(input, read_comment=False):
        # Perform alignment and store results
        for hits in ref.map(seq):
            if hits.is_primary:
                result_per_i[name] = (hits.ctg)
                alignment_counter['mapped'] += 1
                total_mapq += hits.mapq
            elif hits.ctg is None:
                alignment_counter['unmapped'] += 1
        #else non primaru but mapped to something else (secondary, supplementary)
        
    #record the ctg_counts dictionary to a tsv file using polars
    df = pl.from_dict(result_per_i)
    df = df.transpose(include_header=True)
    #add sample column on the first column and fill the value with user's input (sample name)
    df.columns = ['read', 'contig']
    df.write_csv(alignment_output, separator='\t')
    #logging.info(f"Alignment successfully saved to {alignment_output}/{sample_id}_alignment.tsv")

    # Calculate metadata
    total_reads = alignment_counter["mapped"] + alignment_counter["unmapped"]
    mapping_percentage = (alignment_counter["mapped"] / total_reads) * 100 if total_reads > 0 else 0
    average_mapq = total_mapq / alignment_counter["mapped"] if alignment_counter["mapped"] > 0 else 0

    # Create metadata DataFrame using Polars
    metadata = pl.DataFrame({
    'sample-id': [sample_id],
    'total reads': [total_reads],
    'mapped reads': [alignment_counter['mapped']],
    'unmapped reads': [alignment_counter['unmapped']],
    'mapping percentage': [mapping_percentage],
    'average mapq': [average_mapq]
    })

    metadata.write_csv(metadata_output, separator='\t')

if __name__ == "__main__":
    parent_parser = argparse.ArgumentParser(description='Run mappy for both paired and single-end')
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")
    parent_parser.add_argument('-r', '--reference', help='Path to the index/reference file', required=True)
    parent_parser.add_argument('-i', '--input', help='Path to the first/forward read file(s) in fasta/q or in fastq.gz format', required=True)
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
    process_alignment(args.input, args.reference, args.alignment_output, args.metadata_output, args.sample_id)