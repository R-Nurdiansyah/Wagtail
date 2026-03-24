# =============================================================================
# Wagtail Pipeline v0.15
# 16S amplicon sequencing analysis pipeline (Illumina single-end only)
# Reference database: Microflora Danica
# =============================================================================

import os
from datetime import datetime

# -----------------------------------------------------------------------------
# Directory layout
# -----------------------------------------------------------------------------
parent_dir  = os.path.dirname(workflow.basedir)
db_dir      = os.path.join(parent_dir, "database")
data_dir    = os.path.join(parent_dir, "data")
script_dir  = os.path.join(parent_dir, "bin")
run_dir     = os.path.join(parent_dir, "run")

# -----------------------------------------------------------------------------
# Config
# -----------------------------------------------------------------------------
configfile: "config.yaml"

run_name      = config["run_name"]
filenames     = [l.strip() for l in open(config["sample_list"])]
timestamp     = datetime.now().strftime("%Y%m%d")

# Derived paths used in multiple rules
RUN           = f"{run_dir}/{run_name}"
TMP_ROOT      = f"{RUN}/0_tmp"
LOG_ROOT      = f"{RUN}/0_logs_wagtail"
sqlite_db_path = f"{LOG_ROOT}/{run_name}_{timestamp}_log.sql"

# -----------------------------------------------------------------------------
# Scripts & environments
# -----------------------------------------------------------------------------
LOGGER_SCRIPT     = f"{script_dir}/wagtail_sqlite_logger.py"
CHECKER_SCRIPT    = f"{script_dir}/check_sqlite_status.py"
ANNOTATION_SCRIPT = f"{script_dir}/error_annotation.py"

QIIME_ENV  = "envs/qiime2-amplicon-2023.9-py38-linux-conda.yml"
MAPPY_ENV  = "envs/mappy.yaml"

# -----------------------------------------------------------------------------
# Wildcard constraint shared by all per-sample rules
# -----------------------------------------------------------------------------
WILDCARD_FILENAME = r"[^\.]+"

# -----------------------------------------------------------------------------
# Soft-fail bash helpers (injected into every shell block)
# -----------------------------------------------------------------------------
# soft_fail_upstream: propagates a FAILED status from an upstream rule.
#   Touches output file(s), writes a log entry, then exits 0 so Snakemake
#   considers the rule "done" (soft fail rather than hard abort).
#
# run_and_log_soft_fail: runs a command, catches non-zero exit, annotates the
#   error, writes a SQLite log entry, then exits 0.
#   On success it also writes an OK log entry.
#   Either way it removes all non-log artefacts from the log directory so
#   intermediate clutter is minimised as the pipeline runs.
SOFT_FAIL_HELPERS = r"""
soft_fail_upstream() {
    local sqlite_status="$1"
    local touch_cmd="$2"
    local sqlite_db="$3"
    local sample="$4"
    local rule_name="$5"
    local log_file="$6"
    local logger="$7"

    if [ "$sqlite_status" = "FAILED" ]; then
        eval "$touch_cmd"
        echo "Upstream failed" > "$log_file"
        python "$logger" "$sqlite_db" "$sample" "$rule_name" FAILED "$log_file" ""
        find "$(dirname "$log_file")" -type f \
            ! -name "$(basename "$log_file")" \
            ! -name "*.log" \
            ! -name "*.sql" \
            -delete
        exit 0
    fi
}

run_and_log_soft_fail() {
    local run_cmd="$1"
    local fail_cmd="$2"
    local sqlite_db="$3"
    local sample="$4"
    local rule_name="$5"
    local log_file="$6"
    local logger="$7"
    local annotation="$8"

    set +e
    eval "$run_cmd" &> "$log_file"
    local exit_status=$?
    set -e

    if [ $exit_status -ne 0 ]; then
        eval "$fail_cmd"
        local error_msg
        error_msg=$(python "$annotation" "$log_file" "$rule_name")
        python "$logger" "$sqlite_db" "$sample" "$rule_name" FAILED "$log_file" "$error_msg"
    else
        python "$logger" "$sqlite_db" "$sample" "$rule_name" OK "$log_file" ""
    fi

    sleep 2
    find "$(dirname "$log_file")" -type f \
        ! -name "$(basename "$log_file")" \
        ! -name "*.log" \
        ! -name "*.sql" \
        -delete
}
"""

# =============================================================================
# Rules
# =============================================================================

rule all:
    """Final targets that drive the entire DAG."""
    localrule: True
    input:
        expand(f"{RUN}/6_condensed_wagtail/{{filename}}_condensed.tsv", filename=filenames),
        f"{RUN}/7_metadata/{run_name}_full_metadata.tsv",
        sqlite_db_path


# -----------------------------------------------------------------------------
# Step 0 – Manifest creation
# -----------------------------------------------------------------------------
rule manifest:
    """Create a QIIME 2 single-end manifest CSV for each sample."""
    group: "wagtail"
    input:
        filemap = config["file_map"]
    output:
        output_dir = directory(f"{RUN}/0_manifest/{{filename}}"),
        manifest   = f"{RUN}/0_manifest/{{filename}}/{{filename}}_manifest.csv"
    wildcard_constraints:
        filename = WILDCARD_FILENAME
    conda:
        MAPPY_ENV
    log:
        f"{LOG_ROOT}/{{filename}}/manifest_creation.log"
    benchmark:
        f"{LOG_ROOT}/{{filename}}/manifest_creation.benchmark.log"
    threads: 1
    resources:
        mem_mb  = 320,
        runtime = "1h"
    params:
        script     = f"{script_dir}/create_manifest_wagtail.py",
        sample     = "{filename}",
        tmpdir     = f"{TMP_ROOT}/",
        sqlite_db  = f"{TMP_ROOT}/{{filename}}_log.sql",
        rule_name  = "manifest",
        logger     = LOGGER_SCRIPT,
        checker    = CHECKER_SCRIPT,
        annotation = ANNOTATION_SCRIPT
    shell:
        r"""
        {SOFT_FAIL_HELPERS}
        mkdir -p {params.tmpdir}
        run_and_log_soft_fail \
            "python {params.script} \
                --input    {params.sample} \
                --file-map {input.filemap} \
                --output   {output.output_dir}" \
            "mkdir -p {output.output_dir} && touch {output.manifest}" \
            "{params.sqlite_db}" "{wildcards.filename}" "{params.rule_name}" \
            "{log}" "{params.logger}" "{params.annotation}"
        """


# -----------------------------------------------------------------------------
# Step 1 – QIIME 2 import
# -----------------------------------------------------------------------------
rule qiime2_import:
    """Import demultiplexed FASTQ into a QIIME 2 artifact."""
    group: "wagtail"
    input:
        manifest = f"{RUN}/0_manifest/{{filename}}/{{filename}}_manifest.csv"
    output:
        qza = f"{RUN}/1_import_wagtail/{{filename}}-single.qza"
    wildcard_constraints:
        filename = WILDCARD_FILENAME
    conda:
        QIIME_ENV
    log:
        f"{LOG_ROOT}/{{filename}}/qiime2_import.log"
    benchmark:
        f"{LOG_ROOT}/{{filename}}/qiime2_import.benchmark.log"
    threads: 1
    resources:
        mem_mb  = 640,
        runtime = "1h"
    params:
        sqlite_db     = f"{TMP_ROOT}/{{filename}}_log.sql",
        rule_name     = "qiime2_import",
        upstream_rule = "manifest",
        tmpdir        = f"{TMP_ROOT}/{{filename}}",
        logger        = LOGGER_SCRIPT,
        checker       = CHECKER_SCRIPT,
        annotation    = ANNOTATION_SCRIPT
    shell:
        r"""
        {SOFT_FAIL_HELPERS}
        sqlite_status=$(python {params.checker} {params.sqlite_db} \
            "{wildcards.filename}" "{params.upstream_rule}")
        soft_fail_upstream "$sqlite_status" \
            "touch {output.qza}" \
            "{params.sqlite_db}" "{wildcards.filename}" "{params.rule_name}" \
            "{log}" "{params.logger}"

        export TMPDIR={params.tmpdir}
        mkdir -p "$TMPDIR"
        run_and_log_soft_fail \
            "qiime tools import \
                --type         'SampleData[SequencesWithQuality]' \
                --input-path   {input.manifest} \
                --input-format SingleEndFastqManifestPhred33 \
                --output-path  {output.qza}" \
            "touch {output.qza}" \
            "{params.sqlite_db}" "{wildcards.filename}" "{params.rule_name}" \
            "{log}" "{params.logger}" "{params.annotation}"
        """


# -----------------------------------------------------------------------------
# Step 2 – Quality control
# -----------------------------------------------------------------------------
rule quality_control:
    """Filter reads by Q-score using QIIME 2 quality-filter."""
    group: "wagtail"
    input:
        qza = f"{RUN}/1_import_wagtail/{{filename}}-single.qza"
    output:
        filtered = f"{RUN}/2_qc_wagtail/{{filename}}-filtered.qza",
        stats    = f"{RUN}/2_qc_wagtail/{{filename}}-qc-stats.qza"
    wildcard_constraints:
        filename = WILDCARD_FILENAME
    conda:
        QIIME_ENV
    log:
        f"{LOG_ROOT}/{{filename}}/quality_control.log"
    benchmark:
        f"{LOG_ROOT}/{{filename}}/quality_control.benchmark.log"
    threads: 1
    resources:
        mem_mb  = 2000,
        runtime = "1h"
    params:
        sqlite_db     = f"{TMP_ROOT}/{{filename}}_log.sql",
        rule_name     = "quality_control",
        upstream_rule = "qiime2_import",
        tmpdir        = f"{TMP_ROOT}/{{filename}}",
        logger        = LOGGER_SCRIPT,
        checker       = CHECKER_SCRIPT,
        annotation    = ANNOTATION_SCRIPT
    shell:
        r"""
        {SOFT_FAIL_HELPERS}
        sqlite_status=$(python {params.checker} {params.sqlite_db} \
            "{wildcards.filename}" "{params.upstream_rule}")
        soft_fail_upstream "$sqlite_status" \
            "touch {output.filtered} {output.stats}" \
            "{params.sqlite_db}" "{wildcards.filename}" "{params.rule_name}" \
            "{log}" "{params.logger}"

        export TMPDIR={params.tmpdir}
        run_and_log_soft_fail \
            "qiime quality-filter q-score \
                --i-demux              {input.qza} \
                --o-filtered-sequences {output.filtered} \
                --o-filter-stats       {output.stats}" \
            "touch {output.filtered} {output.stats}" \
            "{params.sqlite_db}" "{wildcards.filename}" "{params.rule_name}" \
            "{log}" "{params.logger}" "{params.annotation}"
        """


# -----------------------------------------------------------------------------
# Step 3 – Deblur denoising
# -----------------------------------------------------------------------------
rule deblur:
    """Denoise with Deblur to produce representative sequences and a table."""
    group: "wagtail"
    input:
        filtered = f"{RUN}/2_qc_wagtail/{{filename}}-filtered.qza"
    output:
        representative = f"{RUN}/3_deblur_wagtail/{{filename}}-rep-seqs.qza",
        table          = f"{RUN}/3_deblur_wagtail/{{filename}}-table.qza",
        stats          = f"{RUN}/3_deblur_wagtail/{{filename}}-deblur-stats.qza"
    wildcard_constraints:
        filename = WILDCARD_FILENAME
    conda:
        QIIME_ENV
    log:
        f"{LOG_ROOT}/{{filename}}/deblur.log"
    benchmark:
        f"{LOG_ROOT}/{{filename}}/deblur.benchmark.log"
    threads: 1
    resources:
        mem_mb  = 2000,
        runtime = "47h"
    params:
        path          = directory(f"{RUN}/3_deblur_wagtail/{{filename}}"),
        script        = f"{script_dir}/deblur_all.py",
        sqlite_db     = f"{TMP_ROOT}/{{filename}}_log.sql",
        rule_name     = "deblur",
        upstream_rule = "quality_control",
        tmpdir        = f"{TMP_ROOT}/{{filename}}",
        logger        = LOGGER_SCRIPT,
        checker       = CHECKER_SCRIPT,
        annotation    = ANNOTATION_SCRIPT
    shell:
        r"""
        {SOFT_FAIL_HELPERS}
        sqlite_status=$(python {params.checker} {params.sqlite_db} \
            "{wildcards.filename}" "{params.upstream_rule}")
        soft_fail_upstream "$sqlite_status" \
            "mkdir -p {params.path} && touch {output.representative} {output.table} {output.stats}" \
            "{params.sqlite_db}" "{wildcards.filename}" "{params.rule_name}" \
            "{log}" "{params.logger}"

        export TMPDIR={params.tmpdir}
        run_and_log_soft_fail \
            "python {params.script} \
                -i {input.filtered} \
                -o {params.path} \
                -t {threads} \
                -r {output.representative} \
                -a {output.table} \
                -s {output.stats}" \
            "mkdir -p {params.path} && touch {output.representative} {output.table} {output.stats}" \
            "{params.sqlite_db}" "{wildcards.filename}" "{params.rule_name}" \
            "{log}" "{params.logger}" "{params.annotation}"
        """


# -----------------------------------------------------------------------------
# Step 4a – Export representative sequences
# -----------------------------------------------------------------------------
rule export_seqs:
    """Export representative FASTA sequences from the Deblur QZA artifact."""
    group: "wagtail"
    input:
        representative = f"{RUN}/3_deblur_wagtail/{{filename}}-rep-seqs.qza"
    output:
        output_dir = directory(f"{RUN}/4_rep_seqs_wagtail/{{filename}}"),
        fasta      = f"{RUN}/4_rep_seqs_wagtail/{{filename}}/dna-sequences.fasta"
    wildcard_constraints:
        filename = WILDCARD_FILENAME
    conda:
        QIIME_ENV
    log:
        f"{LOG_ROOT}/{{filename}}/export_seqs.log"
    benchmark:
        f"{LOG_ROOT}/{{filename}}/export_seqs.benchmark.log"
    threads: 1
    resources:
        mem_mb  = 2000,
        runtime = "1h"
    params:
        script        = f"{script_dir}/qiime2_seqs_export.py",
        sqlite_db     = f"{TMP_ROOT}/{{filename}}_log.sql",
        rule_name     = "export_seqs",
        upstream_rule = "deblur",
        tmpdir        = f"{TMP_ROOT}/{{filename}}",
        logger        = LOGGER_SCRIPT,
        checker       = CHECKER_SCRIPT,
        annotation    = ANNOTATION_SCRIPT
    shell:
        r"""
        {SOFT_FAIL_HELPERS}
        sqlite_status=$(python {params.checker} {params.sqlite_db} \
            "{wildcards.filename}" "{params.upstream_rule}")
        soft_fail_upstream "$sqlite_status" \
            "mkdir -p {output.output_dir} && touch {output.fasta}" \
            "{params.sqlite_db}" "{wildcards.filename}" "{params.rule_name}" \
            "{log}" "{params.logger}"

        export TMPDIR={params.tmpdir}
        run_and_log_soft_fail \
            "python {params.script} \
                --input-path  {input.representative} \
                --output-path {output.output_dir}" \
            "mkdir -p {output.output_dir} && touch {output.fasta}" \
            "{params.sqlite_db}" "{wildcards.filename}" "{params.rule_name}" \
            "{log}" "{params.logger}" "{params.annotation}"
        """


# -----------------------------------------------------------------------------
# Step 4b – Export and edit abundance table
# -----------------------------------------------------------------------------
rule export_and_edit_table:
    """Export BIOM table from Deblur QZA and convert to TSV."""
    group: "wagtail"
    input:
        table = f"{RUN}/3_deblur_wagtail/{{filename}}-table.qza"
    output:
        folder = directory(f"{RUN}/4_table_wagtail/{{filename}}_biom"),
        biom   = f"{RUN}/4_table_wagtail/{{filename}}_biom/{{filename}}.biom",
        table  = f"{RUN}/4_table_wagtail/{{filename}}_table.tsv",
        edited = f"{RUN}/5_taxonomy_wagtail/{{filename}}/{{filename}}_table_edit.tsv"
    wildcard_constraints:
        filename = WILDCARD_FILENAME
    conda:
        QIIME_ENV
    log:
        f"{LOG_ROOT}/{{filename}}/export_and_edit_table.log"
    benchmark:
        f"{LOG_ROOT}/{{filename}}/export_and_edit_table.benchmark.log"
    threads: 1
    resources:
        mem_mb  = 1000,
        runtime = "1h"
    params:
        script        = f"{script_dir}/qiime2_biom_export.py",
        biom_name     = "{filename}.biom",
        sqlite_db     = f"{TMP_ROOT}/{{filename}}_log.sql",
        rule_name     = "export_and_edit_table",
        upstream_rule = "deblur",
        tmpdir        = f"{TMP_ROOT}/{{filename}}",
        logger        = LOGGER_SCRIPT,
        checker       = CHECKER_SCRIPT,
        annotation    = ANNOTATION_SCRIPT
    shell:
        r"""
        {SOFT_FAIL_HELPERS}
        sqlite_status=$(python {params.checker} {params.sqlite_db} \
            "{wildcards.filename}" "{params.upstream_rule}")
        soft_fail_upstream "$sqlite_status" \
            "mkdir -p {output.folder} && touch {output.biom} {output.table} {output.edited}" \
            "{params.sqlite_db}" "{wildcards.filename}" "{params.rule_name}" \
            "{log}" "{params.logger}"

        export TMPDIR={params.tmpdir}
        run_and_log_soft_fail \
            "python {params.script} \
                --input-path   {input.table} \
                --output-path  {output.folder} \
                --new-filename {params.biom_name} \
            && biom convert \
                -i {output.biom} \
                -o {output.table} \
                --to-tsv \
            && tail -n +3 {output.table} > {output.edited}" \
            "mkdir -p {output.folder} && touch {output.biom} {output.table} {output.edited}" \
            "{params.sqlite_db}" "{wildcards.filename}" "{params.rule_name}" \
            "{log}" "{params.logger}" "{params.annotation}"
        """


# -----------------------------------------------------------------------------
# Step 5 – Minimap2 alignment (mappy)
# -----------------------------------------------------------------------------
rule mappy:
    """Align representative sequences against the reference database with Mappy."""
    group: "wagtail"
    input:
        fasta = f"{RUN}/4_rep_seqs_wagtail/{{filename}}/dna-sequences.fasta",
        ref   = f"{db_dir}/danica.mmi"
    output:
        align = f"{RUN}/5_taxonomy_wagtail/{{filename}}/{{filename}}_alignment.tsv",
        meta  = f"{RUN}/5_taxonomy_wagtail/{{filename}}/{{filename}}_metadata.tsv"
    wildcard_constraints:
        filename = WILDCARD_FILENAME
    conda:
        MAPPY_ENV
    log:
        f"{LOG_ROOT}/{{filename}}/mappy.log"
    benchmark:
        f"{LOG_ROOT}/{{filename}}/mappy.benchmark.log"
    threads: 1
    resources:
        mem_mb  = 4000,
        runtime = "24h"
    params:
        script        = f"{script_dir}/mappy_script.py",
        sample        = "{filename}",
        sqlite_db     = f"{TMP_ROOT}/{{filename}}_log.sql",
        rule_name     = "mappy",
        upstream_rule = "export_seqs",
        logger        = LOGGER_SCRIPT,
        checker       = CHECKER_SCRIPT,
        annotation    = ANNOTATION_SCRIPT
    shell:
        r"""
        {SOFT_FAIL_HELPERS}
        sqlite_status=$(python {params.checker} {params.sqlite_db} \
            "{wildcards.filename}" "{params.upstream_rule}")
        soft_fail_upstream "$sqlite_status" \
            "touch {output.align} {output.meta}" \
            "{params.sqlite_db}" "{wildcards.filename}" "{params.rule_name}" \
            "{log}" "{params.logger}"

        run_and_log_soft_fail \
            "python {params.script} \
                -i {input.fasta} \
                -r {input.ref} \
                -a {output.align} \
                -m {output.meta} \
                -s {params.sample}" \
            "touch {output.align} {output.meta}" \
            "{params.sqlite_db}" "{wildcards.filename}" "{params.rule_name}" \
            "{log}" "{params.logger}" "{params.annotation}"
        """


# -----------------------------------------------------------------------------
# Step 6 – Extract taxonomy
# -----------------------------------------------------------------------------
rule extract_taxonomy:
    """Combine alignment and abundance table into a condensed taxonomy profile."""
    group: "wagtail"
    input:
        alignment = f"{RUN}/5_taxonomy_wagtail/{{filename}}/{{filename}}_alignment.tsv",
        table     = f"{RUN}/5_taxonomy_wagtail/{{filename}}/{{filename}}_table_edit.tsv"
    output:
        condensed = f"{RUN}/6_condensed_wagtail/{{filename}}_condensed.tsv"
    wildcard_constraints:
        filename = WILDCARD_FILENAME
    conda:
        MAPPY_ENV
    log:
        f"{LOG_ROOT}/{{filename}}/extract_taxonomy.log"
    benchmark:
        f"{LOG_ROOT}/{{filename}}/extract_taxonomy.benchmark.log"
    threads: 1
    resources:
        mem_mb  = 160,
        runtime = "1h"
    params:
        sample        = "{filename}",
        script        = f"{script_dir}/extract_taxonomy_danica.py",
        sqlite_db     = f"{TMP_ROOT}/{{filename}}_log.sql",
        rule_name     = "extract_taxonomy",
        upstream_rule = "mappy",
        logger        = LOGGER_SCRIPT,
        checker       = CHECKER_SCRIPT,
        annotation    = ANNOTATION_SCRIPT
    shell:
        r"""
        {SOFT_FAIL_HELPERS}
        sqlite_status=$(python {params.checker} {params.sqlite_db} \
            "{wildcards.filename}" "{params.upstream_rule}")
        soft_fail_upstream "$sqlite_status" \
            "touch {output.condensed}" \
            "{params.sqlite_db}" "{wildcards.filename}" "{params.rule_name}" \
            "{log}" "{params.logger}"

        run_and_log_soft_fail \
            "python {params.script} \
                -i {input.alignment} \
                -t {input.table} \
                -o {output.condensed} \
                -s {params.sample}" \
            "touch {output.condensed}" \
            "{params.sqlite_db}" "{wildcards.filename}" "{params.rule_name}" \
            "{log}" "{params.logger}" "{params.annotation}"
        """


# -----------------------------------------------------------------------------
# Step 7a – Export QIIME 2 stats (QC + Deblur)
# -----------------------------------------------------------------------------
rule qiime_stats:
    """Export QC and Deblur statistics from QZA artifacts to CSV."""
    group: "wagtail"
    input:
        qc     = f"{RUN}/2_qc_wagtail/{{filename}}-qc-stats.qza",
        deblur = f"{RUN}/3_deblur_wagtail/{{filename}}-deblur-stats.qza"
    output:
        qc_dir     = directory(f"{RUN}/2_qc_wagtail/{{filename}}"),
        deblur_dir = directory(f"{RUN}/3_deblur_wagtail/{{filename}}"),
        qc_csv     = f"{RUN}/2_qc_wagtail/{{filename}}/stats.csv",
        deblur_csv = f"{RUN}/3_deblur_wagtail/{{filename}}/stats.csv"
    wildcard_constraints:
        filename = WILDCARD_FILENAME
    conda:
        QIIME_ENV
    log:
        f"{LOG_ROOT}/{{filename}}/qiime_stats.log"
    benchmark:
        f"{LOG_ROOT}/{{filename}}/qiime_stats.benchmark.log"
    threads: 1
    resources:
        mem_mb  = 640,
        runtime = "1h"
    params:
        script        = f"{script_dir}/wagtail_metadata_qiimes.py",
        sqlite_db     = f"{TMP_ROOT}/{{filename}}_log.sql",
        rule_name     = "qiime_stats",
        upstream_rule = "deblur",
        tmpdir        = f"{TMP_ROOT}/{{filename}}",
        logger        = LOGGER_SCRIPT,
        checker       = CHECKER_SCRIPT,
        annotation    = ANNOTATION_SCRIPT
    shell:
        r"""
        {SOFT_FAIL_HELPERS}

        # (1) Propagate upstream failure
        sqlite_status=$(python {params.checker} {params.sqlite_db} \
            "{wildcards.filename}" "{params.upstream_rule}")
        soft_fail_upstream "$sqlite_status" \
            "mkdir -p {output.qc_dir} {output.deblur_dir} && touch {output.qc_csv} {output.deblur_csv}" \
            "{params.sqlite_db}" "{wildcards.filename}" "{params.rule_name}" \
            "{log}" "{params.logger}"

        # (2) Guard: both QZA inputs must be non-empty to proceed
        #     (upstream OK but Deblur produced an empty artifact)
        if [ ! -s {input.qc} ] || [ ! -s {input.deblur} ]; then
            echo "Empty input QZA — marking as FAILED" > "{log}"
            python {params.logger} {params.sqlite_db} "{wildcards.filename}" \
                "{params.rule_name}" FAILED "{log}" "Empty input QZA"
            mkdir -p {output.qc_dir} {output.deblur_dir}
            touch {output.qc_csv} {output.deblur_csv}
            exit 0
        fi

        export TMPDIR={params.tmpdir}
        run_and_log_soft_fail \
            "python {params.script} \
                --input-qc      {input.qc} \
                --output-qc     {output.qc_dir} \
                --input-deblur  {input.deblur} \
                --output-deblur {output.deblur_dir}" \
            "mkdir -p {output.qc_dir} {output.deblur_dir} && touch {output.qc_csv} {output.deblur_csv}" \
            "{params.sqlite_db}" "{wildcards.filename}" "{params.rule_name}" \
            "{log}" "{params.logger}" "{params.annotation}"
        """


# -----------------------------------------------------------------------------
# Step 7b – Per-sample metadata
# -----------------------------------------------------------------------------
rule metadata_creation:
    """Combine QC, Deblur and alignment metadata into one per-sample TSV."""
    group: "wagtail"
    input:
        qc_csv     = f"{RUN}/2_qc_wagtail/{{filename}}/stats.csv",
        deblur_csv = f"{RUN}/3_deblur_wagtail/{{filename}}/stats.csv",
        mappy_meta = f"{RUN}/5_taxonomy_wagtail/{{filename}}/{{filename}}_metadata.tsv"
    output:
        directory = directory(f"{RUN}/7_metadata/{{filename}}"),
        metadata  = f"{RUN}/7_metadata/{{filename}}/{{filename}}_metadata.tsv"
    wildcard_constraints:
        filename = WILDCARD_FILENAME
    conda:
        MAPPY_ENV
    log:
        f"{LOG_ROOT}/{{filename}}/metadata.log"
    benchmark:
        f"{LOG_ROOT}/{{filename}}/metadata.benchmark.log"
    threads: 1
    resources:
        mem_mb  = 320,
        runtime = "1h"
    params:
        script        = f"{script_dir}/wagtail_metadata_meta_combine.py",
        run           = "{filename}",
        sqlite_db     = f"{TMP_ROOT}/{{filename}}_log.sql",
        rule_name     = "metadata_creation",
        upstream_rule = "qiime_stats",
        tmpdir        = f"{TMP_ROOT}/{{filename}}",
        logger        = LOGGER_SCRIPT,
        checker       = CHECKER_SCRIPT,
        annotation    = ANNOTATION_SCRIPT
    shell:
        r"""
        {SOFT_FAIL_HELPERS}
        sqlite_status=$(python {params.checker} {params.sqlite_db} \
            "{wildcards.filename}" "{params.upstream_rule}")
        soft_fail_upstream "$sqlite_status" \
            "mkdir -p {output.directory} && touch {output.metadata}" \
            "{params.sqlite_db}" "{wildcards.filename}" "{params.rule_name}" \
            "{log}" "{params.logger}"

        run_and_log_soft_fail \
            "python {params.script} \
                --qc-file     {input.qc_csv} \
                --deblur-file {input.deblur_csv} \
                --mappy-file  {input.mappy_meta} \
                --output-path {output.directory} \
                --run-name    {params.run}" \
            "mkdir -p {output.directory} && touch {output.metadata}" \
            "{params.sqlite_db}" "{wildcards.filename}" "{params.rule_name}" \
            "{log}" "{params.logger}" "{params.annotation}"
        """


# -----------------------------------------------------------------------------
# Step 8 – Cleanup & merge
# -----------------------------------------------------------------------------
rule cleanup:
    """
    Final rule — runs after all per-sample rules are complete:
      1. Merges per-sample SQLite logs into one run-level database.
      2. Combines per-sample metadata TSVs into a single run-level file.
      3. Removes all intermediate directories to free inodes/disk space.
    """
    group: "cleanup"
    input:
        metadata  = expand(f"{RUN}/7_metadata/{{filename}}/{{filename}}_metadata.tsv", filename=filenames),
        community = expand(f"{RUN}/6_condensed_wagtail/{{filename}}_condensed.tsv",    filename=filenames)
    output:
        merged        = sqlite_db_path,
        cleanup_done  = f"{LOG_ROOT}/cleanup_done.txt",
        full_metadata = f"{RUN}/7_metadata/{run_name}_full_metadata.tsv"
    wildcard_constraints:
        filename = WILDCARD_FILENAME
    conda:
        MAPPY_ENV
    log:
        merger        = f"{LOG_ROOT}/merger.log",
        cleaning      = f"{LOG_ROOT}/cleaning.log",
        full_metadata = f"{LOG_ROOT}/full_metadata.log"
    threads: 1
    resources:
        mem_mb  = 2000,
        runtime = "2h"
    params:
        merge_script = f"{script_dir}/sql_merge.py",
        per_sample_dbs = expand(f"{TMP_ROOT}/{{filename}}_log.sql", filename=filenames),
        run_root     = RUN,
        tmpdir       = TMP_ROOT
    shell:
        r"""
        # ------------------------------------------------------------------
        # 1. Merge per-sample SQLite logs
        # ------------------------------------------------------------------
        python {params.merge_script} {output.merged} {params.per_sample_dbs} \
            &> {log.merger}

        # ------------------------------------------------------------------
        # 2. Combine per-sample metadata into one run-level TSV
        #    Strategy: write header from the first file, then append data
        #    rows (skip header) from all files.
        # ------------------------------------------------------------------
        (
            head -n 1 {input.metadata[0]}
            tail -n +2 -q {input.metadata}
        ) > {output.full_metadata} 2> {log.full_metadata}

        # ------------------------------------------------------------------
        # 3. Remove intermediate directories
        #    Targets:
        #      0_manifest/          - per-sample manifest CSVs
        #      0_tmp/               - per-sample SQLite logs
        #      1_import_wagtail/    - raw QZA imports
        #      2_qc_wagtail/        - QC QZAs and exported stats
        #      3_deblur_wagtail/    - Deblur QZAs and exported stats
        #      4_rep_seqs_wagtail/  - exported FASTA files
        #      4_table_wagtail/     - exported BIOM/TSV tables
        #      5_taxonomy_wagtail/  - alignment and intermediate TSVs
        # ------------------------------------------------------------------
        rm -rf \
            {params.run_root}/0_manifest \
            {params.tmpdir} \
            {params.run_root}/1_import_wagtail \
            {params.run_root}/2_qc_wagtail \
            {params.run_root}/3_deblur_wagtail \
            {params.run_root}/4_rep_seqs_wagtail \
            {params.run_root}/4_table_wagtail \
            {params.run_root}/5_taxonomy_wagtail \
            2> {log.cleaning}

        touch {output.cleanup_done}
        echo "Cleanup complete at $(date)" >> {log.cleaning}
        """
