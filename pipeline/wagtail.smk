#snakemake file for wagtail pipeline version 0.09 (rule reduction, 12 rules instead of 14, grouping of rules for HPC, and local rule declared correctly)

###tools import###
import os
from datetime import datetime

###data, scripts and database directory###
parent_dir = os.path.dirname(workflow.basedir)
db_dir = os.path.join(parent_dir, "database")
data_dir = os.path.join(parent_dir, "data")
script_dir = os.path.join(parent_dir, "bin")
run_dir = os.path.join(parent_dir, "run")

###config file###
configfile: "config.yaml"

###variable list###
#take the filename list from config.yaml
filename_list = config["sample_list"]
#take the accession list from the user defined file, make sure it is just the forward end
filenames = [line.strip() for line in open(filename_list)]
#take the run name from config.yaml
run_name = config["run_name"]

###create timestamp###
def generate_timestamp():
    #get current date
    return datetime.now().strftime('%Y%m%d')

###rules###
#rule all to run all the rules at once
rule all:
    input:
        expand(f"{run_dir}/{run_name}/6_condensed_wagtail/{{filename}}_condensed.tsv", filename = filenames),
        (f"{run_dir}/{run_name}/7_metadata/{run_name}_full_metadata.tsv")
    localrule: True

rule manifest:
    localrule: True
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
    params:
        #defining the used script for this rule
        script = f"{script_dir}/create_manifest_wagtail.py",
        sample = "{filename}",
        error_log = f"{run_dir}/{run_name}/0_logs_wagtail/{generate_timestamp()}_status.log" # To store any error within the pipeline
    conda:
        "envs/mappy.yaml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/manifest_creation.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/manifest_creation.benchmark.txt"
    threads:
        1
    resources:
        mem_mb = 320,
        runtime = "1h"
    shell:
        "python {params.script} --input {params.sample} --file-map {input.filemap} "
        "--output {output.output_dir} &> {log} || echo '{wildcards.filename} Error at manifest' >> {params.error_log}"
   
rule qiime2_import: #read the data inside the directory
#import data using manifest into qiime2 artifact and single end data only
    localrule: True
    input:
        manifest = f"{run_dir}/{run_name}/0_manifest/{{filename}}/{{filename}}_manifest.csv"
    output:
        f"{run_dir}/{run_name}/1_import_wagtail/{{filename}}-single.qza"
    wildcard_constraints:
        filename = r"[^\.]+"  # Regex to ensure no '.' in 'filename' wildcard
    params:
        error_log = f"{run_dir}/{run_name}/0_logs_wagtail/{generate_timestamp()}_status.log" # To store any error within the pipeline
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
        mem_mb = 640,
        runtime = "1h"
    shell:
        "qiime tools import --type 'SampleData[SequencesWithQuality]' "
        "--input-path {input.manifest} "
        "--input-format SingleEndFastqManifestPhred33 "
        "--output-path {output} &> {log} || echo '{wildcards.filename} Error at qiime_import' >> {params.error_log}"

rule quality_control:
    group: "deblur"
    #run the quality control using qiime2 quality-filter q-score
    input: 
        f"{run_dir}/{run_name}/1_import_wagtail/{{filename}}-single.qza"
    output: 
        filtered = f"{run_dir}/{run_name}/2_qc_wagtail/{{filename}}-filtered.qza",
        stats = f"{run_dir}/{run_name}/2_qc_wagtail/{{filename}}-qc-stats.qza"
    wildcard_constraints:
        filename = r"[^\.]+"  # Regex to ensure no '.' in 'filename' wildcard
    params:
        error_log = f"{run_dir}/{run_name}/0_logs_wagtail/{generate_timestamp()}_status.log" # To store any error within the pipeline
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/quality_control.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/quality_control.benchmark.txt"
    threads:
        1
    resources:
        mem_mb = 4000,
        runtime = "1h"
    shell:
        "qiime quality-filter q-score --i-demux {input} "
        "--o-filtered-sequences {output.filtered} --o-filter-stats {output.stats} &> {log} || echo '{wildcards.filename} Error at quality_control' >> {params.error_log}"

rule deblur:
    group: "deblur"
#script version
    input: 
        f"{run_dir}/{run_name}/2_qc_wagtail/{{filename}}-filtered.qza"
    output: 
        path = directory(f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}"),
        representative = f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}-rep-seqs.qza",
        table = f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}-table.qza",
        stats = f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}-deblur-stats.qza"
    wildcard_constraints:
        filename = r"[^\.]+"  # Regex to ensure no '.' in 'filename' wildcard
    params:
        script = f"{script_dir}/deblur_all.py",
        error_log = f"{run_dir}/{run_name}/0_logs_wagtail/{generate_timestamp()}_status.log" # To store any error within the pipeline
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/deblur.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/deblur.benchmark.txt"
    threads:
        2
    resources:
        mem_mb = 4000,
        runtime = "47h"
    shell:
        "python {params.script} -i {input} -o {output.path} -t {threads} "
        "-r {output.representative} -a {output.table} "
        "-s {output.stats} &> {log} || echo '{wildcards.filename} Error at deblur' >> {params.error_log}"

rule export_seqs:
    group: "mappy"
#export the biom table from the abundance table using the custom script
    input: 
        representative = f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}-rep-seqs.qza"
    output:
        output = directory(f"{run_dir}/{run_name}/4_rep_seqs_wagtail/{{filename}}"),
        final = f"{run_dir}/{run_name}/4_rep_seqs_wagtail/{{filename}}/dna-sequences.fasta"
    wildcard_constraints:
        filename = r"[^\.]+"  # Regex to ensure no '.' in 'filename' wildcard 
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    params:
        #defining the used script for this rule
        script = f"{script_dir}/qiime2_seqs_export.py",
        error_log = f"{run_dir}/{run_name}/0_logs_wagtail/{generate_timestamp()}_status.log" # To store any error within the pipeline
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/export_seqs.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/export_seqs.benchmark.txt"
    threads:
        1
    resources:
        mem_mb = 1000,
        runtime = "1h"
    shell:
        "python {params.script} --input-path {input.representative} "
        "--output-path {output.output} &> {log} || echo '{wildcards.filename} Error at export_seqs' >> {params.error_log}"

rule mappy:
    group: "mappy"
    input:
        #as the product is one file, we can use single input
        F = f"{run_dir}/{run_name}/4_rep_seqs_wagtail/{{filename}}/dna-sequences.fasta",
    #this rule is dependent on the reference being loaded
        ref = f"{db_dir}/danica.mmi"
    output:
        align = f"{run_dir}/{run_name}/5_taxonomy_wagtail/{{filename}}/{{filename}}_alignment.tsv",
        meta = f"{run_dir}/{run_name}/5_taxonomy_wagtail/{{filename}}/{{filename}}_metadata.tsv"
    wildcard_constraints:
        filename = r"[^\.]+"  # Regex to ensure no '.' in 'filename' wildcard
    conda:
         "envs/mappy.yaml"
    params:
        #defining the used script for this rule
        script = f"{script_dir}/mappy_script.py",
        sample = "{filename}",
        error_log = f"{run_dir}/{run_name}/0_logs_wagtail/{generate_timestamp()}_status.log" # To store any error within the pipeline
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/mappy.log"
        #log and benchmark files are located in 0_logs folder and the subdirectory with user-defined name
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/mappy.benchmark.txt"
    threads:
        1
    resources:
        mem_mb = 7000,
        runtime = "24h"
    shell:
        "python {params.script} -i {input.F} "
        "-r {input.ref} -a {output.align} -m {output.meta} "
        "-s {params.sample} &> {log} || echo '{wildcards.filename} Error at mappy' >> {params.error_log}"

rule export_and_edit_table:
    localrule: True
#export the biom table from the abundance table using the custom script
    input: 
        table = f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}-table.qza"
    output:
        folder = directory(f"{run_dir}/{run_name}/4_table_wagtail/{{filename}}_biom"),
        biom = f"{run_dir}/{run_name}/4_table_wagtail/{{filename}}_biom/{{filename}}.biom",
        table = f"{run_dir}/{run_name}/4_table_wagtail/{{filename}}_table.tsv",
        edited = f"{run_dir}/{run_name}/5_taxonomy_wagtail/{{filename}}/{{filename}}_table_edit.tsv"
    wildcard_constraints:
        filename = r"[^\.]+"  # Regex to ensure no '.' in 'filename' wildcard
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    params:
        #defining the used script for this rule
        script = f"{script_dir}/qiime2_biom_export.py",
        #parameter to change the filename for the script
        name = "{filename}.biom",
        error_log = f"{run_dir}/{run_name}/0_logs_wagtail/{generate_timestamp()}_status.log" # To store any error within the pipeline
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/export_and_edit_table.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/export_and_edit_table.benchmark.txt"
    threads:
        1
    resources:
        mem_mb = 1000,
        runtime = "1h"
    shell:
        "python {params.script} --input-path {input.table} "
        "--output-path {output.folder} "
        "--new-filename {params.name} && "
        "biom convert -i {output.biom} -o {output.table} --to-tsv && "
        "(tail -n +3 {output.table} > {output.edited}) &> {log} || echo '{wildcards.filename} Error at export_and_edit_table' >> {params.error_log}"

rule extract_taxonomy:
    localrule: True
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
        script = f"{script_dir}/extract_taxonomy_danica.py",
        error_log = f"{run_dir}/{run_name}/0_logs_wagtail/{generate_timestamp()}_status.log" # To store any error within the pipeline
    conda:
        "envs/mappy.yaml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/extract_taxonomy.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/extract_taxonomy.benchmark.txt"
    threads:
        1
    resources:
        mem_mb = 160,
        runtime = "1h"
    shell:
        "python {params.script} -i {input.primary} "
        "-t {input.table} "
        "-o {output} "
        "-s {params.sample} &> {log} || echo '{wildcards.filename} Error at extract_taxonomy' >> {params.error_log}"

rule qiime_stats:
    localrule: True
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
    params:
        #defining the used script for this rule
        script = f"{script_dir}/wagtail_metadata_qiimes.py",
        error_log = f"{run_dir}/{run_name}/0_logs_wagtail/{generate_timestamp()}_status.log" # To store any error within the pipeline
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/qiime_stats.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/qiime_stats.benchmark.txt"
    threads:
        1
    resources:
        mem_mb = 640,
        runtime = "1h"
    shell:
        "python {params.script} --input-qc {input.qc} --output-qc {output.qc_dir} "
        "--input-deblur {input.deblur} --output-deblur {output.deblur_dir} &> {log} || echo '{wildcards.filename} Error at qiime_stats' >> {params.error_log}"

rule metadata_creation:
    localrule: True
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
    params:
        #defining the used script for this rule
        script = f"{script_dir}/wagtail_metadata_meta_combine.py",
        run = "{filename}",
        error_log = f"{run_dir}/{run_name}/0_logs_wagtail/{generate_timestamp()}_status.log" # To store any error within the pipeline
    conda:
        "envs/mappy.yaml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/metadata.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/metadata.benchmark.txt"
    threads:
        1
    resources:
        mem_mb = 320,
        runtime = "1h"
    shell:
        "python {params.script} --qc-file {input.final_qc} --deblur-file {input.final_deblur} "
        "--mappy-file {input.final_mappy} --output-path {output.directory} "
        "--run-name {params.run} &> {log} || echo '{wildcards.filename} Error at metadata' >> {params.error_log}"

rule clean_intermidiates:
    localrule: True
#rule to clean every in-between result after metadata for each filename is created. cleaning target is 1_* until 5_*
    input:
        metadata = expand(f"{run_dir}/{run_name}/7_metadata/{{filename}}/{{filename}}_metadata.tsv", filename = filenames),
        community = expand(f"{run_dir}/{run_name}/6_condensed_wagtail/{{filename}}_condensed.tsv", filename = filenames)
    output:
        touch(f"{run_dir}/{run_name}/cleanup_done.txt")
    wildcard_constraints:
        filename = r"[^\.]+"  # Regex to ensure no '.' in 'filename' wildcard
    conda:
        "envs/mappy.yaml"
    params:
        target = f"{run_dir}/{run_name}/"
    threads:
        1
    resources:
        mem_mb = 160,
        runtime = "1h"
    shell:
        "rm -rf {params.target}/[1-5]_* && touch {output}"

rule metadata_combine:
    localrule: True
#combine all metadata files into one
    input:
        metadata = expand(f"{run_dir}/{run_name}/7_metadata/{{filename}}/{{filename}}_metadata.tsv", filename = filenames),
        cleanup = f"{run_dir}/{run_name}/cleanup_done.txt"
    output:
        f"{run_dir}/{run_name}/7_metadata/{run_name}_full_metadata.tsv"
    wildcard_constraints:
        filename = r"[^\.]+"  # Regex to ensure no '.' in 'filename' wildcard
    conda:
        "envs/mappy.yaml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/full_metadata_wagtail.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/full_metadata_wagtail.benchmark.txt"
    threads:
        1
    resources:
        mem_mb = 160,
        runtime = "1h"
    shell:
        "head -n 1 {input.metadata[0]} > {output} && tail -n +2 -q {input.metadata} >> {output} &> {log}"