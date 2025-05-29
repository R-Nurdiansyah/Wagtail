#snakemake file for wagtail pipeline version 0.13 (sqlite separate per sample and then merge in the end)

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

timestamp = generate_timestamp()
sqlite_db_path = f"{run_dir}/{run_name}/0_logs_wagtail/{run_name}_{timestamp}_log.sql"

# Change: Use per-sample SQLite database for logging, then merge at the end.

###rules###
#rule all to run all the rules at once
rule all:
    input:
        expand(f"{run_dir}/{run_name}/6_condensed_wagtail/{{filename}}_condensed.tsv", filename = filenames),
        f"{run_dir}/{run_name}/7_metadata/{run_name}_full_metadata.tsv",
        # Add merged log as a final target so it is included in the DAG
        f"{run_dir}/{run_name}/0_logs_wagtail/{run_name}_{timestamp}_log.sql"
    localrule: True

rule manifest:
    group: "initialize"
    #read the data in filenames and check with file map
#create manifest file for each accession, single end data only
    input:
        #list of filenames and filemap
        filemap = config["file_map"]
    output:
        #create manifest file for each accession in the data_dir
        output_dir = directory(f"{run_dir}/{run_name}/0_manifest/{{filename}}"),
        manifest = f"{run_dir}/{run_name}/0_manifest/{{filename}}/{{filename}}_manifest.csv",
    wildcard_constraints:
        filename = r"[^\.]+"  # Regex to ensure no '.' in 'filename' wildcard
    params:
        #defining the used script for this rule
        script = f"{script_dir}/create_manifest_wagtail.py",
        sample = "{filename}",
        tmpdir = f"{run_dir}/{run_name}/0_tmp/",
        sqlite_db = f"{run_dir}/{run_name}/0_tmp/{{filename}}_log.sql",
        rule_name = "manifest",
        logger = f"{script_dir}/wagtail_sqlite_logger.py",
        checker = f"{script_dir}/check_sqlite_status.py"
    conda:
        "envs/mappy.yaml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/manifest_creation.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/manifest_creation.benchmark.log"
    threads:
        1
    resources:
        mem_mb = 320,
        runtime = "1h"
    shell:
        r"""
        # Ensure the directory for the per-sample sqlite_db exists before logging
        mkdir -p {params.tmpdir}
        set +e
        python {params.script} --input {params.sample} --file-map {input.filemap} --output {output.output_dir} &> {log}
        status=$?
        set -e
        if [[ $status -ne 0 ]]; then
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log}
            touch {output.manifest}
        else
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} OK {log}
        fi
        sleep 2
        find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
        """

rule qiime2_import:
    group: "initialize"
    input:
        manifest = f"{run_dir}/{run_name}/0_manifest/{{filename}}/{{filename}}_manifest.csv"
    output:
        qza = f"{run_dir}/{run_name}/1_import_wagtail/{{filename}}-single.qza"
    wildcard_constraints:
        filename = r"[^\.]+"
    params:
        sqlite_db = f"{run_dir}/{run_name}/0_tmp/{{filename}}_log.sql",
        rule_name = "qiime2_import",
        upstream_rule = "manifest",
        tmpdir = f"{run_dir}/{run_name}/0_tmp/{{filename}}",
        logger = f"{script_dir}/wagtail_sqlite_logger.py",
        checker = f"{script_dir}/check_sqlite_status.py"
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/qiime2_import.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/qiime2_import.benchmark.log"
    threads:
        1
    resources:
        mem_mb = 640,
        runtime = "1h"
    shell:
        r"""
        sqlite_status=$(python {params.checker} {params.sqlite_db} "{wildcards.filename}" "{params.upstream_rule}")
        if [ "$sqlite_status" = "FAILED" ]; then
            touch {output.qza}
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log}
            find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
            exit 0
        fi
        export TMPDIR={params.tmpdir}
        mkdir -p $TMPDIR
        set +e
        qiime tools import --type 'SampleData[SequencesWithQuality]' --input-path {input.manifest} --input-format SingleEndFastqManifestPhred33 --output-path {output.qza} &> {log}
        status=$?
        set -e
        if [[ $status -ne 0 ]]; then
            touch {output.qza}
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log}
        else
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} OK {log}
        fi
        sleep 2
        find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
        """

rule quality_control:
    group: "initialize"
    #run the quality control using qiime2 quality-filter q-score
    input: 
        qza = f"{run_dir}/{run_name}/1_import_wagtail/{{filename}}-single.qza"
    output:
        filtered = f"{run_dir}/{run_name}/2_qc_wagtail/{{filename}}-filtered.qza",
        stats = f"{run_dir}/{run_name}/2_qc_wagtail/{{filename}}-qc-stats.qza"
    wildcard_constraints:
        filename = r"[^\.]+"  # Regex to ensure no '.' in 'filename' wildcard
    params:
        sqlite_db = f"{run_dir}/{run_name}/0_tmp/{{filename}}_log.sql",
        rule_name = "quality_control",
        upstream_rule = "qiime2_import",
        tmpdir = f"{run_dir}/{run_name}/0_tmp/{{filename}}",
        logger = f"{script_dir}/wagtail_sqlite_logger.py",
        checker = f"{script_dir}/check_sqlite_status.py"
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/quality_control.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/quality_control.benchmark.log"
    threads:
        1
    resources:
        mem_mb = 2000,
        runtime = "1h"
    shell:
        r"""
        sqlite_status=$(python {params.checker} {params.sqlite_db} "{wildcards.filename}" "{params.upstream_rule}")
        if [ "$sqlite_status" = "FAILED" ]; then
            touch {output.filtered} {output.stats}
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log}
            find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
            exit 0
        fi
        export TMPDIR={params.tmpdir}
        set +e
        qiime quality-filter q-score --i-demux {input.qza} --o-filtered-sequences {output.filtered} --o-filter-stats {output.stats} &> {log}
        status=$?
        set -e
        if [[ $status -ne 0 ]]; then
            touch {output.filtered} {output.stats}
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log}
        else
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} OK {log}
        fi
        sleep 2
        find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
        """

rule deblur:
    group: "initialize"
#script version
    input: 
        filtered = f"{run_dir}/{run_name}/2_qc_wagtail/{{filename}}-filtered.qza"
    output:
        path = directory(f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}"),
        representative = f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}-rep-seqs.qza",
        table = f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}-table.qza",
        stats = f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}-deblur-stats.qza"
    wildcard_constraints:
        filename = r"[^\.]+"  # Regex to ensure no '.' in 'filename' wildcard
    params:
        script = f"{script_dir}/deblur_all.py",
        sqlite_db = f"{run_dir}/{run_name}/0_tmp/{{filename}}_log.sql",
        rule_name = "deblur",
        upstream_rule = "quality_control",
        tmpdir = f"{run_dir}/{run_name}/0_tmp/{{filename}}",
        logger = f"{script_dir}/wagtail_sqlite_logger.py",
        checker = f"{script_dir}/check_sqlite_status.py"
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/deblur.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/deblur.benchmark.log"
    threads:
        1
    resources:
        mem_mb = 2000,
        runtime = "47h"
    shell:
        r"""
        sqlite_status=$(python {params.checker} {params.sqlite_db} "{wildcards.filename}" "{params.upstream_rule}")
        if [ "$sqlite_status" = "FAILED" ]; then
            touch {output.representative} {output.table} {output.stats}
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log}
            find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
            exit 0
        fi
        export TMPDIR={params.tmpdir}
        set +e
        python {params.script} -i {input.filtered} -o {output.path} -t {threads} -r {output.representative} -a {output.table} -s {output.stats} &> {log}
        status=$?
        set -e
        if [[ $status -ne 0 ]]; then
            touch {output.representative} {output.table} {output.stats}
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log}
        else
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} OK {log}
        fi
        sleep 2
        find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
        """

rule export_seqs:
    group: "per_sample"
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
        script = f"{script_dir}/qiime2_seqs_export.py",
        sqlite_db = f"{run_dir}/{run_name}/0_tmp/{{filename}}_log.sql",
        rule_name = "export_seqs",
        upstream_rule = "deblur",
        tmpdir = f"{run_dir}/{run_name}/0_tmp/{{filename}}",
        logger = f"{script_dir}/wagtail_sqlite_logger.py",
        checker = f"{script_dir}/check_sqlite_status.py"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/export_seqs.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/export_seqs.benchmark.log"
    threads:
        1
    resources:
        mem_mb = 2000,
        runtime = "1h"
    shell:
        r"""
        sqlite_status=$(python {params.checker} {params.sqlite_db} "{wildcards.filename}" "{params.upstream_rule}")
        if [ "$sqlite_status" = "FAILED" ]; then
            touch {output.final}
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log}
            find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
            exit 0
        fi
        export TMPDIR={params.tmpdir}
        set +e
        python {params.script} --input-path {input.representative} --output-path {output.output} &> {log}
        status=$?
        set -e
        if [[ $status -ne 0 ]]; then
            touch {output.final}
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log}
        else
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} OK {log}
        fi
        sleep 2
        find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
        """

rule mappy:
    group: "per_sample"
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
        sqlite_db = f"{run_dir}/{run_name}/0_tmp/{{filename}}_log.sql",
        rule_name = "mappy",
        upstream_rule = "export_seqs",
        logger = f"{script_dir}/wagtail_sqlite_logger.py",
        checker = f"{script_dir}/check_sqlite_status.py"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/mappy.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/mappy.benchmark.log"
    threads:
        1
    resources:
        mem_mb = 4000,
        runtime = "24h"
    shell:
        r"""
        sqlite_status=$(python {params.checker} {params.sqlite_db} "{wildcards.filename}" "{params.upstream_rule}")
        if [ "$sqlite_status" = "FAILED" ]; then
            touch {output.align} {output.meta}
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log}
            find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
            exit 0
        fi
        set +e
        python {params.script} -i {input.F} -r {input.ref} -a {output.align} -m {output.meta} -s {params.sample} &> {log}
        status=$?
        set -e
        if [[ $status -ne 0 ]]; then
            touch {output.align} {output.meta}
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log}
        else
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} OK {log}
        fi
        sleep 2
        find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
        """

rule export_and_edit_table:
    group: "per_sample"
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
        sqlite_db = f"{run_dir}/{run_name}/0_tmp/{{filename}}_log.sql",
        rule_name = "export_and_edit_table",
        upstream_rule = "deblur",
        tmpdir = f"{run_dir}/{run_name}/0_tmp/{{filename}}",
        logger = f"{script_dir}/wagtail_sqlite_logger.py",
        checker = f"{script_dir}/check_sqlite_status.py"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/export_and_edit_table.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/export_and_edit_table.benchmark.log"
    threads:
        1
    resources:
        mem_mb = 1000,
        runtime = "1h"
    shell:
        r"""
        sqlite_status=$(python {params.checker} {params.sqlite_db} "{wildcards.filename}" "{params.upstream_rule}")
        if [ "$sqlite_status" = "FAILED" ]; then
            touch {output.biom} {output.table} {output.edited}
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log}
            find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
            exit 0
        fi
        export TMPDIR={params.tmpdir}
        set +e
        python {params.script} --input-path {input.table} --output-path {output.folder} --new-filename {params.name} && biom convert -i {output.biom} -o {output.table} --to-tsv && (tail -n +3 {output.table} > {output.edited}) &> {log}
        status=$?
        set -e
        if [[ $status -ne 0 ]]; then
            touch {output.biom} {output.table} {output.edited}
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log}
        else
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} OK {log}
        fi
        sleep 2
        find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
        """

rule extract_taxonomy:
    group: "per_sample"
#extract taxonomy data based on the database and format the output for filling the taxonomy
    input:
        primary = f"{run_dir}/{run_name}/5_taxonomy_wagtail/{{filename}}/{{filename}}_alignment.tsv",
        table = f"{run_dir}/{run_name}/5_taxonomy_wagtail/{{filename}}/{{filename}}_table_edit.tsv"
    output:
        condensed = f"{run_dir}/{run_name}/6_condensed_wagtail/{{filename}}_condensed.tsv"
    wildcard_constraints:
        filename = r"[^\.]+"  # Regex to ensure no '.' in 'filename' wildcard
    params:
        sample = "{filename}",
        script = f"{script_dir}/extract_taxonomy_danica.py",
        sqlite_db = f"{run_dir}/{run_name}/0_tmp/{{filename}}_log.sql",
        rule_name = "extract_taxonomy",
        upstream_rule = "mappy",
        logger = f"{script_dir}/wagtail_sqlite_logger.py",
        checker = f"{script_dir}/check_sqlite_status.py"
    conda:
        "envs/mappy.yaml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/extract_taxonomy.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/extract_taxonomy.benchmark.log"
    threads:
        1
    resources:
        mem_mb = 160,
        runtime = "1h"
    shell:
        r"""
        sqlite_status=$(python {params.checker} {params.sqlite_db} "{wildcards.filename}" "{params.upstream_rule}")
        if [ "$sqlite_status" = "FAILED" ]; then
            echo "FAILED" > {output.condensed}
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log}
            find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
            exit 0
        fi
        set +e
        python {params.script} -i {input.primary} -t {input.table} -o {output.condensed} -s {params.sample} &> {log}
        status=$?
        set -e
        if [[ $status -ne 0 ]]; then
            touch {output.condensed}
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log}
        else
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} OK {log}
        fi
        sleep 2
        find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
        """

rule qiime_stats:
    group: "per_sample"
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
        script = f"{script_dir}/wagtail_metadata_qiimes.py",
        sqlite_db = f"{run_dir}/{run_name}/0_tmp/{{filename}}_log.sql",
        rule_name = "qiime_stats",
        upstream_rule = "deblur",
        tmpdir = f"{run_dir}/{run_name}/0_tmp/{{filename}}",
        logger = f"{script_dir}/wagtail_sqlite_logger.py",
        checker = f"{script_dir}/check_sqlite_status.py"
    conda:
        "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/qiime_stats.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/qiime_stats.benchmark.log"
    threads:
        1
    resources:
        mem_mb = 640,
        runtime = "1h"
    shell:
        r"""
        sqlite_status=$(python {params.checker} {params.sqlite_db} "{wildcards.filename}" "{params.upstream_rule}")
        if [ "$sqlite_status" = "FAILED" ]; then
            echo "Upstream stats missing, skipping qiime_stats" > {output.final_qc}
            echo "Upstream stats missing, skipping qiime_stats" > {output.final_deblur}
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log}
            find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
            exit 0
        fi
        # For qiime_stats, check both inputs are present and non-empty
        if [ ! -s {input.qc} ] || [ ! -s {input.deblur} ]; then
            echo "Upstream stats missing, skipping qiime_stats" > {output.final_qc}
            echo "Upstream stats missing, skipping qiime_stats" > {output.final_deblur}
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log}
            find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
            exit 0
        fi
        export TMPDIR={params.tmpdir}
        set +e
        python {params.script} --input-qc {input.qc} --output-qc {output.qc_dir} --input-deblur {input.deblur} --output-deblur {output.deblur_dir} &> {log}
        status=$?
        set -e
        if [[ $status -ne 0 ]]; then
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log}
        else
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} OK {log}
        fi
        sleep 2
        find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
        """

rule metadata_creation:
    group: "per_sample"
#wrote the metadata for the run
    input: 
        final_qc = f"{run_dir}/{run_name}/2_qc_wagtail/{{filename}}/stats.csv",
        final_deblur = f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}/stats.csv",
        final_mappy = f"{run_dir}/{run_name}/5_taxonomy_wagtail/{{filename}}/{{filename}}_metadata.tsv"
    output: 
        directory = directory(f"{run_dir}/{run_name}/7_metadata/{{filename}}"),
        final = f"{run_dir}/{run_name}/7_metadata/{{filename}}/{{filename}}_metadata.tsv",
        sqlite_db = f"{run_dir}/{run_name}/0_tmp/{{filename}}_log.sql"
    wildcard_constraints:
        filename = r"[^\.]+"  # Regex to ensure no '.' in 'filename' wildcard
    params:
        script = f"{script_dir}/wagtail_metadata_meta_combine.py",
        run = "{filename}",
        rule_name = "metadata_creation",
        upstream_rule = "qiime_stats",
        tmpdir = f"{run_dir}/{run_name}/0_tmp/{{filename}}",
        logger = f"{script_dir}/wagtail_sqlite_logger.py",
        checker = f"{script_dir}/check_sqlite_status.py"
    conda:
        "envs/mappy.yaml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/metadata.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/metadata.benchmark.log"
    threads:
        1
    resources:
        mem_mb = 320,
        runtime = "1h"
    shell:
        r"""
        sqlite_status=$(python {params.checker} {output.sqlite_db} "{wildcards.filename}" "{params.upstream_rule}")
        if [ "$sqlite_status" = "FAILED" ]; then
            mkdir -p {output.directory}
            echo "ERROR: metadata_creation failed for {wildcards.filename}" > {output.final}
            python {params.logger} {output.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log}  
            find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
            exit 0
        fi
        set +e
        python {params.script} --qc-file {input.final_qc} --deblur-file {input.final_deblur} --mappy-file {input.final_mappy} --output-path {output.directory} --run-name {params.run} &> {log}
        status=$?
        set -e
        if [[ $status -ne 0 ]]; then
            echo "ERROR: metadata_creation failed for {wildcards.filename}" > {output.final}
            python {params.logger} {output.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log}
        else
            python {params.logger} {output.sqlite_db} "{wildcards.filename}" {params.rule_name} OK {log}
        fi
        sleep 2
        find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
        """

rule merge_logs:
    localrule: True
    input:
        metadata = expand(f"{run_dir}/{run_name}/7_metadata/{{filename}}/{{filename}}_metadata.tsv", filename = filenames),
        community = expand(f"{run_dir}/{run_name}/6_condensed_wagtail/{{filename}}_condensed.tsv", filename = filenames),
        dbs = expand(f"{run_dir}/{run_name}/0_tmp/{{filename}}_log.sql", filename=filenames)
    output:
        merged = sqlite_db_path
    conda:
        "envs/mappy.yaml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/merger.log"
    threads:
        1
    shell:
        r"""
        # Create merged database and logs table if not exists
        touch {input.metadata[0]}  # Ensure the first metadata file exists
        touch {input.community[0]}  # Ensure the first community file exists
        sqlite3 {output.merged} "CREATE TABLE IF NOT EXISTS logs (id TEXT, rule TEXT, outcome TEXT, log TEXT);"
        # Merge all per-sample logs, but skip if db does not exist or is empty
        for db in {input.dbs}; do
            if [ -s "$db" ]; then
                count=$(sqlite3 "$db" "SELECT count(*) FROM sqlite_master WHERE type='table' AND name='logs';")
                if [ "$count" = "1" ]; then
                    nrows=$(sqlite3 "$db" "SELECT count(*) FROM logs;")
                    if [ "$nrows" -gt 0 ]; then
                        # Use a transaction for better performance and reliability
                        sqlite3 {output.merged} "BEGIN;"
                        sqlite3 {output.merged} "ATTACH DATABASE '$db' AS to_merge; INSERT INTO logs SELECT * FROM to_merge.logs; DETACH DATABASE to_merge;"
                        sqlite3 {output.merged} "COMMIT;"
                    fi
                fi
            fi
        done
        # Wait a bit to ensure file system sync (especially on network filesystems)
        sync
        sleep 1
        """

rule clean_intermidiates:
    localrule: True
    #rule to clean every in-between result after metadata for each filename is created. cleaning target is 1_* until 5_*
    input:
        metadata = expand(f"{run_dir}/{run_name}/7_metadata/{{filename}}/{{filename}}_metadata.tsv", filename = filenames),
        community = expand(f"{run_dir}/{run_name}/6_condensed_wagtail/{{filename}}_condensed.tsv", filename = filenames),
        merged = f"{run_dir}/{run_name}/0_logs_wagtail/{run_name}_{timestamp}_log.sql"
    output:
        touch(f"{run_dir}/{run_name}/0_logs_wagtail/cleanup_done.txt")
    wildcard_constraints:
        filename = r"[^\.]+"
    conda:
        "envs/mappy.yaml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/cleaning.log"
    threads:
        1
    resources:
        mem_mb = 160,
        runtime = "1h"
    params:
        target = f"{run_dir}/{run_name}/",
        tmpdir = f"{run_dir}/{run_name}/0_tmp/"
    shell:
        r"""
        # Remove blank or FAILED condensed files
        for f in {run_dir}/{run_name}/6_condensed_wagtail/*_condensed.tsv; do
            # Remove if file is empty
            if [ ! -s "$f" ]; then
                rm -f "$f"
                continue
            fi
            # Remove if file contains only the word FAILED (from soft fail)
            if grep -q "^FAILED" "$f"; then
                rm -f "$f"
            fi
        done
        # Clean up temporary directory and intermediate files
        rm -rf {params.tmpdir}
        rm -rf {params.target}/[1-5]_*
        """ 

rule metadata_combine:
    localrule: True
#combine all metadata files into one
    input:
        metadata = expand(f"{run_dir}/{run_name}/7_metadata/{{filename}}/{{filename}}_metadata.tsv", filename = filenames),
        cleanup = f"{run_dir}/{run_name}/0_logs_wagtail/cleanup_done.txt"
    output:
        f"{run_dir}/{run_name}/7_metadata/{run_name}_full_metadata.tsv"
    wildcard_constraints:
        filename = r"[^\.]+"  # Regex to ensure no '.' in 'filename' wildcard
    conda:
        "envs/mappy.yaml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/full_metadata.log"
    threads:
        1
    resources:
        mem_mb = 160,
        runtime = "1h"
    shell:
        "(head -n 1 {input.metadata[0]} > {output} && tail -n +2 -q {input.metadata} >> {output}) &> {log}"