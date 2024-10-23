#testing mappy

import polars as pl
import mappy as mp
import os
import argparse
import logging
import sys
from collections import Counter
from collections import OrderedDict

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

#dataclass
#create a class to store the input file, base name, and output directory
class AmpliconDataset:
    input_file: str = None
    basename: str = None
    output_directory: str = None

    def __str__(self):
        return f"Dataset: {self.input_file} -> {self.output_directory}"

def check_input(input):
    # Check if input is a file or a directory and list all files in the directory
    input_files=[]
    printout=[]
    if isinstance (input, str):
        input = [input] #convert to list if input is a string

    for i in input:
        if os.path.isfile(i):
            input_files.append(i)
            #logging.info(f"Input is a file: {input_files}")
        else:
            # List files in the directory
            all_files = os.listdir(i)
            input_files.extend(os.path.join(i, f) for f in all_files if os.path.isfile(os.path.join(i, f)))
            printout.append(all_files)
    logging.info(f"Files in the input directory(ies) {input}: {printout}")
    logging.info(f"Input files: {input_files}")

    input_files = sorted(input_files)
    file_names = [os.path.basename(f) for f in input_files]

    datasets = []

    if single:
        #create a list of base file name without the extension for create_output_directories function
        for file in file_names:
            base = file.split('.')[0]

            d = AmpliconDataset()
            d.input_file = file
            d.basename = base
            datasets.append(d)

    return datasets

#def create_output_directories(output_dir, chunk_length, max_files_per_dir, base_file):
def create_output_directories(output_dir, datasets):
    #iterate through the files and create the directories and if the directory is already exist, it will skip the creation
    for d in datasets:
        file = d.basename
        # Extract parent directory name (first 3 digits after "ERR")
        if len(file) < 6:
            logging.error(f"Filename {file} is too short")
            sys.exit(1)
        parent_dir = file[:6]
    
        # Extract the last digit to determine the subdirectory ("000" to "009")
        sub_dir_num = int(file[-3:]) % 1000
    
        # Create the full path for the accession directory
        new_dir = os.path.join(output_dir, f"{parent_dir}", f"{sub_dir_num}", file)
        d.output_directory = new_dir
          
        #Create new directory if not exist
        if not os.path.exists(new_dir):
            os.makedirs(new_dir, exist_ok=True)
            logging.info(f"Created directory: {new_dir}")
        else:
            logging.info(f"Directory already exists: {new_dir}")

#run the mappy
def mappy_run(index, datasets):
    #load or build index
    ref = mp.Aligner(index, preset = 'sr')
    if not ref: raise Exception("ERROR: failed to load/build index")
    #print to console if index is loaded from the file
    logging.info("Index loaded from the file {}".format(index))
    
    #check the input, if --single is used, then it will run the single-end alignment
    if single:
        #print to console if single-end alignment is choosen
        logging.info("Running single-end alignment")
        for da in datasets:
            i = da.input_file
            d = da.output_directory

            logging.info("Running analysis for {}".format(i))
            result_per_i = {}
            ctg_counts = Counter()
            #Loop through single-end reads
            for name, seq, qual in mp.fastx_read(i):
                barcode = os.path.basename(i).split('.')[0]
                names = "{}_barcode+{}".format(barcode, name)
                #Perform single-end alignment
                for hits in ref.map(seq, cs=True, MD=True):
                    if hits.is_primary:
                        result_per_i[names] = (hits.ctg)
                #count the number of contigs using Counter
                ctg_counts.update(result_per_i.values())
            #record the ctg_counts dictionary to a tsv file using polars
            df = pl.from_dict(ctg_counts)
            df = df.transpose(include_header=True)
            #add sample column on the first column and fill the value with user's input (sample name)
            df.columns = ['contig', 'count']
            df = df.with_columns(pl.Series("sample", [barcode]*len(df)))
            #reorder the column to match the biobox script input while delete the unnecessary columns
            df = df[["sample", "contig", "count"]]
            #sort it by count
            df = df.sort("count", descending=True)
            df.write_csv("{}/{}.csv".format(d, barcode), separator=',')
            logging.info(f"Output saved at: {d}/{barcode}.csv")
    
    # else:
    #     #print to console if paired-end alignment is choosen
    #     logging.info("Running paired-end alignment")
    #     #unzip the input_list to read1 and read2 and make it a list
    #     read1, read2 = zip(*input_list)
    #     read1 = list(read1)
    #     read2 = list(read2)
    #     for i1, i2, d in zip(read1, read2, directories_list):
    #         logging.info("Running analysis for {} + {}".format(i1, i2))
    #         result_per_i = {}
    #         ctg_counts = Counter()
    #         #Loop through paired-end reads
    #         for (name_f, seq_f, qual_f), (name_r, seq_r, qual_r) in zip(mp.fastx_read(i1), mp.fastx_read(i2)):
    #             #take the basename of the input file for barcode 
    #             barcode_f = os.path.basename(i1).split('.')[0]
    #             barcode_r = os.path.basename(i2).split('.')[0]
    #             names = "{}_paired_{}_barcode+{}".format(barcode_f, barcode_r, name_f)
    #         #Perform paired-end alignment
    #             for hits in ref.map(seq_f, seq_r, cs=True, MD=True):
    #                 if hits.is_primary:
    #                     result_per_i[names] = (hits.ctg)
    #             #count the number of contigs using Counter
    #             ctg_counts.update(result_per_i.values())
    #         #record the ctg_counts dictionary to a tsv file using polars
    #         df = pl.from_dict(ctg_counts)
    #         df = df.transpose(include_header=True)
    #         #add sample column on the first column and fill the value with user's input (sample name)
    #         df.columns = ['contig', 'count']
    #         df = df.with_columns(pl.Series("sample", [barcode_f]*len(df)))
    #         #reorder the column to match the biobox script input while delete the unnecessary columns
    #         df = df[["sample", "contig", "count"]]
    #         #sort it by count
    #         df = df.sort("count", descending=True)
    #         df.write_csv("{}/{}+{}.csv".format(d, barcode_f, barcode_r), separator=',')
    #         logging.info(f"Output saved at: {d}/{barcode_f}+{barcode_r}.csv")

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(description='Run mappy for both paired and single-end')
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")
    parent_parser.add_argument('-d', '--index', help='Path to the index or reference file', required=True)
    #parent_parser.add_argument('-s', '--single', help='Run single-end alignment mode only', action="store_true") -> no need as the input would be single end from the start
    parent_parser.add_argument('-i', '--input-directories', help='Input directory or input file path', nargs='+', required=True)
    parent_parser.add_argument('-o', '--output-dir', help='Output directory path', required=True)
    
    args = parent_parser.parse_args()
    
    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')
    
    # Run the function
    datasets = check_input(args.input)
    create_output_directories(args.output_dir, datasets)
    mappy_run(args.index, args.single, datasets)
    