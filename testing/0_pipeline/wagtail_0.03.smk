#snakemake file for wagtail pipeline version 0.03 (config, mappy, metadata)

###tools import###
import os

###data, scripts and database directory###
parent_dir = os.path.dirname(workflow.basedir)
db_dir = os.path.join(parent_dir, "0_database")
data_dir = os.path.join(parent_dir, "0_data")
script_dir = os.path.join(parent_dir, "0_scripts")
metadata_dir = os.path.join(parent_dir, "0_metadata")
configfile: "config.yaml"

###list of accession to be used###
#take the filename list from config.yaml
filename_list = config["sample_list"]
#take the accession list from the user defined file, make sure it is just the forward end
filenames = [line.strip() for line in open(filename_list)]

#rule all to run all the rules at once
rule all:
    input:
        expand("6_condensed_combine/{filename}_condensed.tsv", filename = filenames),
        "0_metadata/{params.run_name}/metadata.tsv"

rule manifest:
#read the data in filenames and check with file map
#create manifest file for each accession, single end data only
    input:
        #list of filenames and filemap
        filemap = config["file_map"]
    output:
        #create manifest file for each accession in the data_dir
        output_dir = directory(f"{data_dir}/manifest/{{filename}}"),
        manifest = f"{data_dir}/manifest/{{filename}}/{{filename}}_manifest.csv"
    params:
        #defining the used script for this rule
        script = f"{script_dir}/create_manifest_wagtail.py"
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    log:
        "0_logs_combine/{filename}/manifest_creation.log"
    benchmark:
        "0_logs_combine/{filename}/manifest_creation.benchmark.txt"
    threads:
        1
    resources:
        mem_mb = 32000,
        runtime = "96h"
    shell:
        "(python {params.script} --input {input.filename} --file_map {input.filemap} "
        "--output {output.output_dir}| touch {output.manifest}) 2> {log}"

rule qiime2_import: #read the data inside the directory
#import data using manifest into qiime2 artifact and single end data only
#need script to create manifest file??
    input:
        manifest = f"{data_dir}/manifest/{{filename}}/{{filename}}_manifest.csv"
    output:
        "1_import_combine/{filename}-single.qza"
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    log:
        "0_logs_combine/{filename}/qiime2_import.log"
        #log and benchmark files are located in 0_logs folder and the subdirectory with user-defined name
    benchmark:
        "0_logs_combine/{filename}/qiime2_import.benchmark.txt"
    threads:
        1
    resources:
        mem_mb = 32000,
        runtime = "96h"
    shell:
        "(qiime tools import --type 'SampleData[SequencesWithQuality]' "
        "--input-path {input.manifest} "
        "--input-format SingleEndFastqManifestPhred33 "
        "--output-path {output}) 2> {log}"

rule quality_control:
    #run the quality control using qiime2 quality-filter q-score
    input: 
        "1_import_combine/{filename}-single.qza"
    output: 
        filtered = "2_qc_combine/{filename}-filtered.qza",
        stats = "2_qc_combine/{filename}-qc-stats.qza"
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    log:
        "0_logs_combine1/{filename}/quality_control.log"
    benchmark:
        "0_logs_combine1/{filename}/quality_control.benchmark.txt"
    threads:
        1
    resources:
        mem_mb = 32000,
        runtime = "96h"
    shell:
        "(qiime quality-filter q-score --i-demux {input} "
        "--o-filtered-sequences {output.filtered} --o-filter-stats {output.stats}) 2> {log}"

rule deblur:
#demultiplex the data to get the representative sequences, OTU table, and denoising stats
#delete the trunc length as it is only single end data
    input: 
        "2_qc_combine/{filename}-filtered.qza"
    output: 
        representative = "3_deblur_combine/{filename}-rep-seqs.qza",
        table = "3_deblur_combine/{filename}-table.qza",
        stats = "3_deblur_combine/{filename}-deblur-stats.qza"
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    log:
        "0_logs_combine/{filename}/qc.log"
    benchmark:
        "0_logs_combine/{filename}/qc.benchmark.txt"
    threads:
        1
    resources:
        mem_mb = 32000,
        runtime = "96h"
    shell:
        "(qiime deblur denoise-16S --i-demultiplexed-seqs {input} "
        "--p-trim-length -1 --p-sample-stats "
        "--p-min-reads 0 "
        "--o-representative-sequences {output.representative} " 
        "--o-table {output.table} " 
        "--o-stats {output.stats}) 2> {log}"

rule export_seqs:
#export the biom table from the abundance table using the custom script
    input: 
        representative = "3_deblur_combine/{filename}-rep-seqs.qza"
    output:
        output = directory("4_rep_seqs_combine/{filename}"),
        final = "4_rep_seqs_combine/{filename}/dna-sequences.fasta"
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    params:
        #defining the used script for this rule
        script = f"{script_dir}/qiime2_seqs_export.py",
    log:
        "0_logs_combine/{filename}/export_seqs.log"
    benchmark:
        "0_logs_combine/{filename}/export_seqs.benchmark.txt"
    threads:
        1
    resources:
        mem_mb = 32000,
        runtime = "96h"
    shell:
        "(python {params.script} --input-path {input.representative} "
        "--output-path {output.output}) 2> {log}"

rule mappy:
#use mappy script in version 0.03
    input:
        #as the product is one file, we can use single input
        F = "4_rep_seqs_combine/{filename}/dna-sequences.fasta",
    #ref is the indexed database from rule indexing
        ref = f"{db_dir}/danica.mmi"
    output:
        directory = directory("5_taxonomy_combine/{filename}"),
        align = "5_taxonomy_combine/{filename}/{filename}_alignment.tsv",
        meta = "5_taxonomy_combine/{filename}/{filename}_metadata.tsv"
    conda:
        "envs/mappy.yaml"
    params:
        #defining the used script for this rule
        script = f"{script_dir}/mappy.py",
        sample = "{filename}"
    log:
        "0_logs_combine/{filename}/mappy.log"
        #log and benchmark files are located in 0_logs folder and the subdirectory with user-defined name
    benchmark:
        "0_logs_combine/{filename}/mappy.benchmark.txt"
    threads:
        1
    resources:
        mem_mb = 16000,
        runtime = "96h"
    shell:
        "(python {params.script} -r {input.F} "
        "-i {input.ref} -a {output.directory} -m {output.directory} "
        "-s {params.sample} 2> {log}"

rule export_table:
#export the biom table from the abundance table using the custom script
    input: 
        table = "3_deblur_combine/{filename}-table.qza"
    output:
        output = directory("4_table_combine/{filename}_biom"),
        final = "4_table_combine/{filename}_biom/{filename}.biom"
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    params:
        #defining the used script for this rule
        script = f"{script_dir}/qiime2_biom_export.py",
        #parameter to change the filename for the script
        name = "{filename}.biom"
    log:
        "0_logs_combine/{filename}/export_table.log"
    benchmark:
        "0_logs_combine/{filename}/export_biom.benchmark.txt"
    threads:
        1
    resources:
        mem_mb = 32000,
        runtime = "96h"
    shell:
        "(python {params.script} --input-path {input.table} "
        "--output-path {output.output} "
        "--new-filename {params.name}) 2> {log}"

rule biom_to_tsv:
#convert the biom table to tsv table
    input:
    #check the directory from the previous rule and read the biom table inside
        "4_table_combine/{filename}_biom/{filename}.biom"
    output: 
        "4_table_combine/{filename}_table.tsv"
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    log:
        "0_logs_combine/{filename}/biom_to_tsv.log"
    benchmark:
        "0_logs_combine/{filename}/biom_to_tsv.benchmark.txt"
    threads:
        1
    resources:
        mem_mb = 32000,
        runtime = "96h"
    shell:
        "(biom convert -i {input} -o {output} --to-tsv) 2> {log}"

rule edit_table:
#to edit the abundance tsv file for input in the next script, take all from 3rd lines
    input: 
        "4_table_combine/{filename}_table.tsv"
    output: 
        "5_taxonomy_combine/{filename}/{filename}_table_edit.tsv"
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    log:
        "0_logs_combine/{filename}/edit_table.log"
    benchmark:
        "0_logs_combine/{filename}/edit_table.benchmark.txt"
    threads:
        1
    resources:
        mem_mb = 32000,
        runtime = "96h"
    shell:
        "(tail -n +3 {input} > {output}) 2> {log}"

rule extract_taxonomy:
#extract taxonomy data based on the database and format the output for filling the taxonomy
    input:
        primary = "5_taxonomy_combine/{filename}/{filename}_alignment.tsv",
        table = "5_taxonomy_combine/{filename}/{filename}_table_edit.tsv"
    output: 
        "6_condensed_combine/{filename}_condensed.tsv"
    params:
        #sample name for the extract_taxonomy.py script -> needed for the script
        sample = "{filename}",
        #defining the used script for this rule
        script = f"{script_dir}/extract_taxonomy_danica.py"
    conda:
        "envs/mappy.yaml"
    log:
        "0_logs_combine/{filename}/extract_taxonomy.log"
    benchmark:
        "0_logs_combine/{filename}/extract_taxonomy.benchmark.txt"
    threads:
        1
    resources:
        mem_mb = 16000,
        runtime = "96h"
    shell:
        "(python {params.script} -i {input.primary} "
        "-t {input.table} "
        "-o {output} "
        "-s {params.sample}) 2> {log}"

rule qiime_stats:
#wrote the metadata for the run
    input: 
        qc = "2_qc_combine/{filename}-qc-stats.qza",
        deblur = "3_deblur_combine/{filename}-deblur-stats.qza"
    output: 
        output_qc = directory("2_qc_combine/{filename}"),
        output_deblur = directory("3_deblur_combine/{filename}"),
        final_qc = "2_qc_combine/{filename}/stats.csv",
        final_deblur = "3_deblur_combine/{filename}/stats.csv"
    params:
        #defining the used script for this rule
        script = f"{script_dir}/wagtail_metadata-qiimes.py"
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    log:
        "0_logs_combine/{filename}/qiime_stats.log"
    benchmark:
        "0_logs_combine/{filename}/qiime_stats.benchmark.txt"
    threads:
        1
    resources:
        mem_mb = 32000,
        runtime = "96h"
    shell:
        "(python {params.script} --input-qc {input.qc} --output-qc {output.output_qc} "
        "--input-deblur {input.deblur} --output-deblur {output.output_deblur} 2> {log}"

rule metadata_creation:
#wrote the metadata for the run
    input: 
        final_qc = "2_qc_combine/{filename}/stats.csv",
        final_deblur = "3_deblur_combine/{filename}/stats.csv",
        final_mappy = "5_taxonomy_combine/{filename}_metadata.tsv"
    output: 
        directory = directory("0_metadata/{filename}"),
        final = "0_metadata/{filename}/{filename}_metadata.tsv"
    params:
        #defining the used script for this rule
        script = f"{script_dir}/wagtail_metadata-meta-combine.py",
        run_name = "{filename}"
    conda:
        "envs/mappy.yaml"
    log:
        "0_logs_combine/{filename}/metadata.log"
    benchmark:
        "0_logs_combine/{filename}/metadata.benchmark.txt"
    threads:
        1
    resources:
        mem_mb = 32000,
        runtime = "96h"
    shell:
        "(python {params.script} --qc-file {input.final_qc} --deblur-file {input.final_deblur} "
        "--mappy-file {input.final_mappy} --output-path {output.directory} "
        "--run-name {params.run_name} 2> {log}"

rule metadata_combine:
#combine all metadata files into one
    input:
        metadata = expand("0_metadata/{filename}/{filename}_metadata.tsv", filename = filenames)
    output:
        "0_metadata/{params.run_name}/metadata.tsv"
    params:
        run_name = config["run_name"] #use the one in the config file
    conda:
        "envs/mappy.yaml"
    log:
        "0_logs_combine/{params.run_name}/metadata_combine.log"
    benchmark:
        "0_logs_combine/{params.run_name}/metadata_combine.benchmark.txt"
    threads:
        1
    resources:
        mem_mb = 32000,
        runtime = "96h"
    shell:
        "(head -n 1 {input.metadata[0]} > {output} && tail -n +2 -q {input.metadata} >> {output}) 2> {log}"