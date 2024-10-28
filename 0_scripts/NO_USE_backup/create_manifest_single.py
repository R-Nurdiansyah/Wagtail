import os
import argparse
import logging
import pandas as pd

def create_manifest(input_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for file in os.listdir(input_dir):
        if file.endswith("-1.fastq.gz"):
            #the file is accession-1.fastq.gz, take only the CAMI_number or 1st and second element of the split
            sample_name = file.split("-1")[0]
            forward_file = os.path.join(input_dir, file)
            #reverse_file = os.path.join(input_dir, sample_name + "-2.fastq.gz")

            # Buat dataframe manifest dengan 3 kolom: sample-id, absolute-filepath, dan direction
            df = pd.DataFrame(columns=["sample-id", "absolute-filepath", "direction"])
                
                # Tambahkan baris untuk read maju
            df.loc[0] = [sample_name, os.path.abspath(forward_file), "forward"]
                # Tambahkan baris untuk read mundur
            #    df.loc[1] = [sample_name, os.path.abspath(reverse_file), "reverse"]
                
                # Simpan manifest sebagai file CSV
            manifest_file = os.path.join(output_dir, sample_name + "_manifest.csv")
            df.to_csv(manifest_file, index=False)
                
            print("Created manifest file for sample:", sample_name)
        else:
            logging.warning(f"Not a forward {file}")

if __name__ == "__main__":
    # Setup parser argumen baris perintah
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', help='input directory', required=True)
    parser.add_argument('--output_dir', help='output directory', required=True)
    args = parser.parse_args()

    # Panggil fungsi untuk membuat manifest
    create_manifest(args.input_dir, args.output_dir)
