import os
import argparse
import logging
import polars as pl

def create_manifest(input, file_map, output):
    #1. read the file map. file map has 2 columns: accession and path. create a dictionary from those 2 columns
    #2. read the input and then match the accession with the dictionary. if the accession is found, then create the manifest file, single end only
    #3. save the manifest file to the output directory

    # Read the file map
    file_map = pl.read_csv(file_map, separator = '\t', has_header=False)
    file_map_dict = dict(zip(file_map['column_1'], file_map['column_2']))

    # Create manifest single end file from input and file map and then save it to the output directory
    if input in file_map_dict:
        forward_file = file_map_dict[input]
        df_manifest = pl.DataFrame({
            "sample-id": [input],
            "absolute-filepath": [os.path.abspath(forward_file)],
            "direction": ["forward"]
        })
        manifest_file = os.path.join(output, input + "_manifest.csv")
        df_manifest.write_csv(manifest_file, has_header=True)
        logging.info("Created manifest file for sample:", input)
    else:
        logging.warning(f"Accession {input} not found in the file map")

if __name__ == "__main__":
    # Setup parser argument
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', help='input accession name', required=True)
    parser.add_argument('--file_map', help='file map location, 2 column of tsv without header. column 1 accession, column 2 absolute path of the data', required=True)
    parser.add_argument('--output', help='output directory', required=True)
    args = parser.parse_args()

    # Call the function
    create_manifest(args.input, args.file_map, args.output)
