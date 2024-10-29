import os
import argparse
import logging

def create_filemap(input, output_filemap, output_accession):
    with open(output_filemap, "w") as f_map, open(output_accession, "w") as f_acc:
        for filename in os.listdir(input):
            if os.path.isfile(os.path.join(input, filename)):
                #name_without_ext = os.path.splitext(filename)[0]  # Nama file tanpa ekstensi
                name_without_ext = filename.split(".")[0]  # Nama file tanpa ekstensi
                absolute_path = os.path.abspath(os.path.join(input, filename))  # Absolute path
                f_map.write(f"{name_without_ext}\t{absolute_path}\n")
                f_acc.write(f"{name_without_ext}\n")

if __name__ == "__main__":
    # Setup parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='input directory', required=True)
    parser.add_argument('-o', '--output-filemap', help='output file in txt', required=True)
    parser.add_argument('-a', '--output-accession', help='output file for accession in txt', required=True)
    args = parser.parse_args()

    # Call function
    create_filemap(args.input, args.output_filemap, args.output_accession)
