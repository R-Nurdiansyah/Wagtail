#snakemake file for wagtail pipeline version 0.03 (config, mappy, metadata)

###tools import###
import os
#from scripts.reference_loader import load_reference

###data, scripts and database directory###
parent_dir = os.path.dirname(workflow.basedir)
db_dir = os.path.join(parent_dir, "0_database")
data_dir = os.path.join(parent_dir, "0_data")
script_dir = os.path.join(parent_dir, "0_scripts")
metadata_dir = os.path.join(parent_dir, "0_metadata")
run_dir = os.path.join(parent_dir, "0_run")

###config file###
configfile: "config.yaml"

###variable list###
#take the filename list from config.yaml
filename_list = config["sample_list"]
#take the accession list from the user defined file, make sure it is just the forward end
filenames = [line.strip() for line in open(filename_list)]
#take the run name from config.yaml
run_name = config["run_name"]

#rule all to run all the rules at once
rule all:
    input:
        expand(f"{run_dir}/{run_name}/6_condensed_wagtail/{{filename}}_condensed.tsv", filename = filenames),
        (f"{run_dir}/{run_name}/7_metadata/full_metadata.tsv")

rule manifest:
#read the data in filenames and check with file map
#create manifest file for each accession, single end data only
    input:
        #list of filenames and filemap
        filemap = config["file_map"]
    output:
        #create manifest file for each accession in the data_dir
        output_dir = directory(f"{run_dir}/{run_name}/0_manifest/{{filename}}"),
        manifest = f"{run_dir}/{run_name}/0_manifest/{{filename}}/{{filename}}_manifest.csv"
    wildcard_constraints:
        filename = r"[^\.]+"  # Regex to ensure no '.' in 'filename' wildcard
    group:
        "initialization"
    params:
        #defining the used script for this rule
        script = f"{script_dir}/create_manifest_wagtail.py",
        sample = "{filename}"
    conda:
        "envs/mappy.yaml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/manifest_creation.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/manifest_creation.benchmark.txt"
    threads:
        1
    resources:
        mem_mb = 32000,
        runtime = "96h"
    shell:
        "(python {params.script} --input {params.sample} --file-map {input.filemap} "
        "--output {output.output_dir}| touch {output.manifest}) 2> {log}"

rule qiime2_import: #read the data inside the directory
#import data using manifest into qiime2 artifact and single end data only
    input:
        manifest = f"{run_dir}/{run_name}/0_manifest/{{filename}}/{{filename}}_manifest.csv"
    output:
        f"{run_dir}/{run_name}/1_import_wagtail/{{filename}}-single.qza"
    wildcard_constraints:
        filename = r"[^\.]+"  # Regex to ensure no '.' in 'filename' wildcard
    group:
        "initialization"
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/qiime2_import.log"
        #log and benchmark files are located in 0_logs folder and the subdirectory with user-defined name
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/qiime2_import.benchmark.txt"
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
        f"{run_dir}/{run_name}/1_import_wagtail/{{filename}}-single.qza"
    output: 
        filtered = f"{run_dir}/{run_name}/2_qc_wagtail/{{filename}}-filtered.qza",
        stats = f"{run_dir}/{run_name}/2_qc_wagtail/{{filename}}-qc-stats.qza"
    wildcard_constraints:
        filename = r"[^\.]+"  # Regex to ensure no '.' in 'filename' wildcard 
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/quality_control.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/quality_control.benchmark.txt"
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
        f"{run_dir}/{run_name}/2_qc_wagtail/{{filename}}-filtered.qza"
    output: 
        representative = f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}-rep-seqs.qza",
        table = f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}-table.qza",
        stats = f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}-deblur-stats.qza"
    wildcard_constraints:
        filename = r"[^\.]+"  # Regex to ensure no '.' in 'filename' wildcard 
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/qc.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/qc.benchmark.txt"
    threads:
        1
    resources:
        mem_mb = 32000,
        runtime = "96h"
    shell:
        "(qiime deblur denoise-16S --i-demultiplexed-seqs {input} "
        #"--p-trim-length -1 --p-sample-stats "
        "--p-trim-length 150 --p-sample-stats "
        "--p-min-reads 0 "
        "--o-representative-sequences {output.representative} " 
        "--o-table {output.table} " 
        "--o-stats {output.stats}) 2> {log}"

rule export_seqs:
#export the biom table from the abundance table using the custom script
    input: 
        representative = f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}-rep-seqs.qza"
    output:
        output = directory(f"{run_dir}/{run_name}/4_rep_seqs_wagtail/{{filename}}"),
        final = f"{run_dir}/{run_name}/4_rep_seqs_wagtail/{{filename}}/dna-sequences.fasta"
    wildcard_constraints:
        filename = r"[^\.]+"  # Regex to ensure no '.' in 'filename' wildcard 
    group:
        "mappy_alignment"
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    params:
        #defining the used script for this rule
        script = f"{script_dir}/qiime2_seqs_export.py",
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/export_seqs.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/export_seqs.benchmark.txt"
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
        F = f"{run_dir}/{run_name}/4_rep_seqs_wagtail/{{filename}}/dna-sequences.fasta",
    #ref is the indexed database from rule indexing
        ref = f"{db_dir}/danica.mmi"
    output:
        align = f"{run_dir}/{run_name}/5_taxonomy_wagtail/{{filename}}/{{filename}}_alignment.tsv",
        meta = f"{run_dir}/{run_name}/5_taxonomy_wagtail/{{filename}}/{{filename}}_metadata.tsv"
    wildcard_constraints:
        filename = r"[^\.]+"  # Regex to ensure no '.' in 'filename' wildcard
    group:
        "mappy_alignment" #make sure no whitespace
    conda:
        "envs/mappy.yaml"
    params:
        #defining the used script for this rule
        script = f"{script_dir}/mappy_script.py",
        sample = "{filename}"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/mappy.log"
        #log and benchmark files are located in 0_logs folder and the subdirectory with user-defined name
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/mappy.benchmark.txt"
    threads:
        1
    resources:
        mem_mb = 32000,
        runtime = "96h"
    shell:
        "(python {params.script} -i {input.F} "
        "-r {input.ref} -a {output.align} -m {output.meta} "
        "-s {params.sample}) 2> {log}"

rule export_table:
#export the biom table from the abundance table using the custom script
    input: 
        table = f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}-table.qza"
    output:
        output = directory(f"{run_dir}/{run_name}/4_table_wagtail/{{filename}}_biom"),
        final = f"{run_dir}/{run_name}/4_table_wagtail/{{filename}}_biom/{{filename}}.biom"
    wildcard_constraints:
        filename = r"[^\.]+"  # Regex to ensure no '.' in 'filename' wildcard
    group:
        "table_creation"
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    params:
        #defining the used script for this rule
        script = f"{script_dir}/qiime2_biom_export.py",
        #parameter to change the filename for the script
        name = "{filename}.biom"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/export_table.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/export_biom.benchmark.txt"
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
        f"{run_dir}/{run_name}/4_table_wagtail/{{filename}}_biom/{{filename}}.biom"
    output: 
        f"{run_dir}/{run_name}/4_table_wagtail/{{filename}}_table.tsv"
    wildcard_constraints:
        filename = r"[^\.]+"  # Regex to ensure no '.' in 'filename' wildcard
    group:
        "table_creation"
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/biom_to_tsv.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/biom_to_tsv.benchmark.txt"
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
        f"{run_dir}/{run_name}/4_table_wagtail/{{filename}}_table.tsv"
    output: 
        f"{run_dir}/{run_name}/5_taxonomy_wagtail/{{filename}}/{{filename}}_table_edit.tsv"
    wildcard_constraints:
        filename = r"[^\.]+"  # Regex to ensure no '.' in 'filename' wildcard
    group:
        "table_creation"
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/edit_table.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/edit_table.benchmark.txt"
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
        primary = f"{run_dir}/{run_name}/5_taxonomy_wagtail/{{filename}}/{{filename}}_alignment.tsv",
        table = f"{run_dir}/{run_name}/5_taxonomy_wagtail/{{filename}}/{{filename}}_table_edit.tsv"
    output: 
        f"{run_dir}/{run_name}/6_condensed_wagtail/{{filename}}_condensed.tsv"
    wildcard_constraints:
        filename = r"[^\.]+"  # Regex to ensure no '.' in 'filename' wildcard
    params:
        #sample name for the extract_taxonomy.py script -> needed for the script
        sample = "{filename}",
        #defining the used script for this rule
        script = f"{script_dir}/extract_taxonomy_danica.py"
    conda:
        "envs/mappy.yaml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/extract_taxonomy.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/extract_taxonomy.benchmark.txt"
    threads:
        1
    resources:
        mem_mb = 32000,
        runtime = "96h"
    shell:
        "(python {params.script} -i {input.primary} "
        "-t {input.table} "
        "-o {output} "
        "-s {params.sample}) 2> {log}"

rule qiime_stats:
#wrote the metadata for the run
    input: 
        qc = f"{run_dir}/{run_name}/2_qc_wagtail/{{filename}}-qc-stats.qza",
        deblur = f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}-deblur-stats.qza"
    output: 
        qc_dir = directory(f"{run_dir}/{run_name}/2_qc_wagtail/{{filename}}"),
        deblur_dir = directory(f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}"),
        final_qc = f"{run_dir}/{run_name}/2_qc_wagtail/{{filename}}/stats.csv",
        final_deblur = f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}/stats.csv"
    wildcard_constraints:
        filename = r"[^\.]+"  # Regex to ensure no '.' in 'filename' wildcard
    group:
        "metadatas"
    params:
        #defining the used script for this rule
        script = f"{script_dir}/wagtail_metadata-qiimes.py"
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/qiime_stats.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/qiime_stats.benchmark.txt"
    threads:
        1
    resources:
        mem_mb = 32000,
        runtime = "96h"
    shell:
        "(python {params.script} --input-qc {input.qc} --output-qc {output.qc_dir} "
        "--input-deblur {input.deblur} --output-deblur {output.deblur_dir}) 2> {log}"

rule metadata_creation:
#wrote the metadata for the run
    input: 
        final_qc = f"{run_dir}/{run_name}/2_qc_wagtail/{{filename}}/stats.csv",
        final_deblur = f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}/stats.csv",
        final_mappy = f"{run_dir}/{run_name}/5_taxonomy_wagtail/{{filename}}/{{filename}}_metadata.tsv"
    output: 
        directory = directory(f"{run_dir}/{run_name}/7_metadata/{{filename}}"),
        final = f"{run_dir}/{run_name}/7_metadata/{{filename}}/{{filename}}_metadata.tsv"
    wildcard_constraints:
        filename = r"[^\.]+"  # Regex to ensure no '.' in 'filename' wildcard
    group:
        "metadatas"
    params:
        #defining the used script for this rule
        script = f"{script_dir}/wagtail_metadata-meta-combine.py",
        run = "{filename}"
    conda:
        "envs/mappy.yaml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/metadata.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/metadata.benchmark.txt"
    threads:
        1
    resources:
        mem_mb = 32000,
        runtime = "96h"
    shell:
        "(python {params.script} --qc-file {input.final_qc} --deblur-file {input.final_deblur} "
        "--mappy-file {input.final_mappy} --output-path {output.directory} "
        "--run-name {params.run}) 2> {log}"

rule metadata_combine:
#combine all metadata files into one
    input:
        metadata = expand(f"{run_dir}/{run_name}/7_metadata/{{filename}}/{{filename}}_metadata.tsv", filename = filenames)
    output:
        f"{run_dir}/{run_name}/7_metadata/full_metadata.tsv"
    wildcard_constraints:
        filename = r"[^\.]+"  # Regex to ensure no '.' in 'filename' wildcard
    group:
        "metadatas"
    conda:
        "envs/mappy.yaml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/metadata_wagtail.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/metadata_wagtail.benchmark.txt"
    threads:
        1
    resources:
        mem_mb = 32000,
        runtime = "96h"
    shell:
        "(head -n 1 {input.metadata[0]} > {output} && tail -n +2 -q {input.metadata} >> {output}) 2> {log}"