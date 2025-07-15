#snakemake file for wagtail pipeline version 0.14 (revised the soft fail and adding error annotation)

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
    group: "wagtail"
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
        checker = f"{script_dir}/check_sqlite_status.py",
        annotation = f"{script_dir}/error_annotation.py"
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
            # If the script fails, create the output directory and touch manifest file
            mkdir -p {output.output_dir}
            touch {output.manifest}
            error=$(python {params.annotation} {log} {params.rule_name})
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log} "$error"
        else
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} OK {log} ""
        fi
        sleep 2
        find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
        """

rule qiime2_import:
    group: "wagtail"
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
        checker = f"{script_dir}/check_sqlite_status.py",
        annotation = f"{script_dir}/error_annotation.py"
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
            echo "Upstream fails" > {log}
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log} ""
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
            error=$(python {params.annotation} {log} {params.rule_name})
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log} "$error"
        else
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} OK {log} ""
        fi
        sleep 2
        find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
        """

rule quality_control:
    group: "wagtail"
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
        checker = f"{script_dir}/check_sqlite_status.py",
        annotation = f"{script_dir}/error_annotation.py"
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
            echo "Upstream fails" > {log}
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log} ""
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
            error=$(python {params.annotation} {log} {params.rule_name})
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log} "$error"
        else
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} OK {log} ""
        fi
        sleep 2
        find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
        """

rule deblur:
    group: "wagtail"
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
        checker = f"{script_dir}/check_sqlite_status.py",
        annotation = f"{script_dir}/error_annotation.py"
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
            mkdir -p {output.path}
            touch {output.representative} {output.table} {output.stats}
            echo "Upstream fails" > {log}
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log} ""
            find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
            exit 0
        fi
        export TMPDIR={params.tmpdir}
        set +e
        python {params.script} -i {input.filtered} -o {output.path} -t {threads} -r {output.representative} -a {output.table} -s {output.stats} &> {log}
        status=$?
        set -e
        if [[ $status -ne 0 ]]; then
            mkdir -p {output.path}
            touch {output.representative} {output.table} {output.stats}
            error=$(python {params.annotation} {log} {params.rule_name})
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log} "$error"
        else
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} OK {log} ""
        fi
        sleep 2
        find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
        """

rule export_seqs:
    group: "wagtail"
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
        checker = f"{script_dir}/check_sqlite_status.py",
        annotation = f"{script_dir}/error_annotation.py"
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
            mkdir -p {output.output}
            touch {output.final}
            echo "Upstream fails" > {log}
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log} ""
            find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
            exit 0
        fi
        export TMPDIR={params.tmpdir}
        set +e
        python {params.script} --input-path {input.representative} --output-path {output.output} &> {log}
        status=$?
        set -e
        if [[ $status -ne 0 ]]; then
            mkdir -p {output.output}
            touch {output.final}
            error=$(python {params.annotation} {log} {params.rule_name})
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log} "$error"
        else
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} OK {log} ""
        fi
        sleep 2
        find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
        """

rule mappy:
    group: "wagtail"
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
        checker = f"{script_dir}/check_sqlite_status.py",
        annotation = f"{script_dir}/error_annotation.py"
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
            echo "Upstream fails" > {log}
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log} ""
            find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
            exit 0
        fi
        set +e
        python {params.script} -i {input.F} -r {input.ref} -a {output.align} -m {output.meta} -s {params.sample} &> {log}
        status=$?
        set -e
        if [[ $status -ne 0 ]]; then
            touch {output.align} {output.meta}
            error=$(python {params.annotation} {log} {params.rule_name})
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log} "$error"
        else
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} OK {log} ""
        fi
        sleep 2
        find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
        """

rule export_and_edit_table:
    group: "wagtail"
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
        checker = f"{script_dir}/check_sqlite_status.py",
        annotation = f"{script_dir}/error_annotation.py"
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
            mkdir -p {output.folder}
            touch {output.biom} {output.table} {output.edited}
            echo "Upstream fails" > {log}
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log} ""
            find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
            exit 0
        fi
        export TMPDIR={params.tmpdir}
        set +e
        python {params.script} --input-path {input.table} --output-path {output.folder} --new-filename {params.name} && biom convert -i {output.biom} -o {output.table} --to-tsv && (tail -n +3 {output.table} > {output.edited}) &> {log}
        status=$?
        set -e
        if [[ $status -ne 0 ]]; then
            mkdir -p {output.folder}
            touch {output.biom} {output.table} {output.edited}
            error=$(python {params.annotation} {log} {params.rule_name})
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log} "$error"
        else
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} OK {log} ""
        fi
        sleep 2
        find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
        """

rule extract_taxonomy:
    group: "wagtail"
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
        checker = f"{script_dir}/check_sqlite_status.py",
        annotation = f"{script_dir}/error_annotation.py"
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
            touch {output.condensed}
            echo "Upstream fails" > {log}
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log} ""
            find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
            exit 0
        fi
        set +e
        python {params.script} -i {input.primary} -t {input.table} -o {output.condensed} -s {params.sample} &> {log}
        status=$?
        set -e
        if [[ $status -ne 0 ]]; then
            touch {output.condensed}
            error=$(python {params.annotation} {log} {params.rule_name})
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log} "$error"
        else
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} OK {log} ""
        fi
        sleep 2
        find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
        """

rule qiime_stats:
    group: "wagtail"
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
        checker = f"{script_dir}/check_sqlite_status.py",
        annotation = f"{script_dir}/error_annotation.py"
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
            mkdir -p {output.qc_dir} {output.deblur_dir}
            touch {output.final_qc} {output.final_deblur}
            echo "Upstream fails" > {log}
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log} ""
            find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
            exit 0
        fi
        # For qiime_stats, check both inputs are present and non-empty
        if [ ! -s {input.qc} ] || [ ! -s {input.deblur} ]; then
            mkdir -p {output.qc_dir} {output.deblur_dir}
            touch {output.final_qc} {output.final_deblur}
            echo "Upstream fails" > {log}
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log} ""
            find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
            exit 0
        fi
        export TMPDIR={params.tmpdir}
        set +e
        python {params.script} --input-qc {input.qc} --output-qc {output.qc_dir} --input-deblur {input.deblur} --output-deblur {output.deblur_dir} &> {log}
        status=$?
        set -e
        if [[ $status -ne 0 ]]; then
            mkdir -p {output.qc_dir} {output.deblur_dir}
            touch {output.final_qc} {output.final_deblur}
            error=$(python {params.annotation} {log} {params.rule_name})
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log} "$error"
        else
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} OK {log} ""
        fi
        sleep 2
        find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
        """

rule metadata_creation:
    group: "wagtail"
#wrote the metadata for the run
    input: 
        final_qc = f"{run_dir}/{run_name}/2_qc_wagtail/{{filename}}/stats.csv",
        final_deblur = f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}/stats.csv",
        final_mappy = f"{run_dir}/{run_name}/5_taxonomy_wagtail/{{filename}}/{{filename}}_metadata.tsv"
    output: 
        directory = directory(f"{run_dir}/{run_name}/7_metadata/{{filename}}"),
        final = f"{run_dir}/{run_name}/7_metadata/{{filename}}/{{filename}}_metadata.tsv",
    wildcard_constraints:
        filename = r"[^\.]+"  # Regex to ensure no '.' in 'filename' wildcard
    params:
        script = f"{script_dir}/wagtail_metadata_meta_combine.py",
        sqlite_db = f"{run_dir}/{run_name}/0_tmp/{{filename}}_log.sql",
        run = "{filename}",
        rule_name = "metadata_creation",
        upstream_rule = "qiime_stats",
        tmpdir = f"{run_dir}/{run_name}/0_tmp/{{filename}}",
        logger = f"{script_dir}/wagtail_sqlite_logger.py",
        checker = f"{script_dir}/check_sqlite_status.py",
        annotation = f"{script_dir}/error_annotation.py"
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
        sqlite_status=$(python {params.checker} {params.sqlite_db} "{wildcards.filename}" "{params.upstream_rule}")
        if [ "$sqlite_status" = "FAILED" ]; then
            mkdir -p {output.directory}
            touch {output.final}
            echo "Upstream fails" > {log}
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log} "" 
            find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
            exit 0
        fi
        set +e
        python {params.script} --qc-file {input.final_qc} --deblur-file {input.final_deblur} --mappy-file {input.final_mappy} --output-path {output.directory} --run-name {params.run} &> {log}
        status=$?
        set -e
        if [[ $status -ne 0 ]]; then
            mkdir -p {output.directory}
            echo "ERROR: metadata_creation failed for {wildcards.filename}" > {output.final}
            error=$(python {params.annotation} {log} {params.rule_name})
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} FAILED {log} "$error"
        else
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" {params.rule_name} OK {log} ""
        fi
        sleep 2
        find $(dirname {log}) -type f ! -name "$(basename {log})" ! -name "*.log" ! -name "*.sql" -delete
        """

rule cleanup:
    group: "cleanup"
    input:
        metadata = expand(f"{run_dir}/{run_name}/7_metadata/{{filename}}/{{filename}}_metadata.tsv", filename = filenames),
        community = expand(f"{run_dir}/{run_name}/6_condensed_wagtail/{{filename}}_condensed.tsv", filename = filenames),
    output:
        merged = sqlite_db_path,
        cleanup_done = f"{run_dir}/{run_name}/0_logs_wagtail/cleanup_done.txt",
        full_metadata = f"{run_dir}/{run_name}/7_metadata/{run_name}_full_metadata.tsv"
    wildcard_constraints:
        filename = r"[^\.]+"  # Regex to ensure no '.' in 'filename' wildcard
    params:
        merge_script = f"{script_dir}/sql_merge.py",
        dbs = expand(f"{run_dir}/{run_name}/0_tmp/{{filename}}_log.sql", filename=filenames),
        target = f"{run_dir}/{run_name}/",
        tmpdir = f"{run_dir}/{run_name}/0_tmp/"
    conda:
        "envs/mappy.yaml"
    log:
        merger = f"{run_dir}/{run_name}/0_logs_wagtail/merger.log",
        cleaning = f"{run_dir}/{run_name}/0_logs_wagtail/cleaning.log",
        full_metadata = f"{run_dir}/{run_name}/0_logs_wagtail/full_metadata.log"
    threads:
        1
    resources:
        mem_mb = 2000,
        runtime = "2h"
    shell:
        r"""
        # Step 1: Merge logs
        touch {input.metadata[0]}  # Ensure the first metadata file exists
        touch {input.community[0]}  # Ensure the first community file exists
        python {params.merge_script} {output.merged} {params.dbs} &> {log.merger}
        sleep 2
        
        # Step 2: Clean intermediates
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
        touch {output.cleanup_done}
        echo "Cleanup completed" &> {log.cleaning}
        
        # Step 3: Combine metadata files
        (head -n 1 {input.metadata[0]} > {output.full_metadata} && tail -n +2 -q {input.metadata} >> {output.full_metadata}) &> {log.full_metadata}
        """