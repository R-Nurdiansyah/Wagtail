# Wagtail pipeline v0.15
# Changes from v0.14:
#   - Bash boilerplate extracted into shared shell functions (wagtail_functions.sh)
#   - Added marker_id rule (mappy-based 16S/18S/ITS/CO1 identification) as pre-flight QC
#   - marker_id runs before manifest; non-16S samples soft-fail early

import os
import glob
import csv
from pathlib import Path
from datetime import datetime

# ── Directories ──────────────────────────────────────────────────────────────
parent_dir = os.path.dirname(workflow.basedir)
db_dir     = os.path.join(parent_dir, "database")
data_dir   = os.path.join(parent_dir, "data")
script_dir = os.path.join(parent_dir, "bin")
run_dir    = os.path.join(parent_dir, "run")

# ── Config ────────────────────────────────────────────────────────────────────
configfile: "config.yaml"

filename_list = config["sample_list"]
filenames     = [line.strip() for line in open(filename_list)]
run_name      = config["run_name"]

# ── Timestamp / SQLite ───────────────────────────────────────────────────
def make_timestamp():
    fmt = "%Y%m%d"
    return datetime.now().strftime(fmt)

timestamp = make_timestamp()
sql_log   = f"{run_dir}/{run_name}/0_logs_wagtail/{run_name}_{timestamp}_log.sql"

# ── Shared path helpers (called from params) ──────────────────────────────────
def log_dir(filename):
    return f"{run_dir}/{run_name}/0_logs_wagtail/{filename}"

def tmp_db(filename):
    return f"{run_dir}/{run_name}/0_tmp/{filename}_log.sql"

def tmp_dir(filename):
    return f"{run_dir}/{run_name}/0_tmp/{filename}"

# Path to the shared bash helper library (generated alongside this .smk)
BASH_LIB = os.path.join(workflow.basedir, "wagtail_functions.sh")

# ── Marker-aware glob helpers ─────────────────────────────────────────────────
def get_marker(wc):
    """Read predicted marker from .marker TSV output.
    Calls checkpoints.marker_id.get() so Snakemake defers DAG construction
    until the checkpoint has run and the file actually exists.
    Returns "16S" as a safe fallback when the file is empty or unreadable
    (e.g. OOM kill) — downstream rules will still soft-fail via SQLite check.
    """
    checkpoints.marker_id.get(filename=wc.filename)
    marker_file = f"{run_dir}/{run_name}/0_marker_id/{wc.filename}.marker"
    try:
        with open(marker_file) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            row = next(reader)
        return row["predicted_marker"]
    except (StopIteration, KeyError, OSError):
        return "16S"

def _resolve(marker: str, db_dir: str, pattern: str) -> str:
    """Resolve a single file matching db_dir/pattern. Raises if 0 or >1 match."""
    matches = glob.glob(f"{db_dir}/{pattern}")
    if not matches:
        raise FileNotFoundError(f"No file matching '{pattern}' found in {db_dir}")
    if len(matches) > 1:
        raise ValueError(f"Multiple files matching '{pattern}' in {db_dir}: {matches}")
    return matches[0]

def get_deblur_db(wc) -> str:
    marker = get_marker(wc)
    if marker == "16S":
        # denoise-16S uses QIIME2's built-in reference — no external database file exists.
        # Return the marker file itself as a harmless dummy so Snakemake can build the DAG;
        # deblur_all.py ignores --reference when marker == 16S.
        return f"{run_dir}/{run_name}/0_marker_id/{wc.filename}.marker"
    return _resolve(marker, db_dir, pattern=f"{marker}_deblur*")

def get_mappy_ref(wc) -> str:
    marker = get_marker(wc)
    return _resolve(marker, db_dir, pattern=f"{marker}_database*")

def get_tax_ref(wc) -> str:
    marker = get_marker(wc)
    return _resolve(marker, db_dir, pattern=f"{marker}_taxonomy*")

DB_16S = _resolve("16S", db_dir, pattern="16S_database*")
DB_18S = _resolve("18S", db_dir, pattern="18S_database*")
DB_ITS = _resolve("ITS", db_dir, pattern="ITS_database*")
DB_CO1 = _resolve("CO1", db_dir, pattern="CO1_database*")

# ── rule all ─────────────────────────────────────────────────────────────────
rule all:
    input:
        expand(f"{run_dir}/{run_name}/6_condensed_wagtail/{{filename}}_condensed.tsv", filename=filenames),
        f"{run_dir}/{run_name}/7_metadata/{run_name}_full_metadata.tsv",
        sql_log
    localrule: True

# ── marker_id ─────────────────────────────────────────────────────────────────
# Pre-flight QC: subsample reads and align against 16S, 18S, ITS, CO1 databases.
# Writes the predicted marker to a .marker file. soft-fails if ambiguous.
checkpoint marker_id:
    group: "wagtail"
    input:
        filemap = config["file_map"],
    output:
        marker = f"{run_dir}/{run_name}/0_marker_id/{{filename}}.marker"
    wildcard_constraints:
        filename = r"[^\.]+"
    params:
        script          = f"{script_dir}/identify_amplicon.py",
        sample          = "{filename}",
        db_16S = DB_16S,
        db_18S = DB_18S,
        db_ITS = DB_ITS,
        db_CO1 = DB_CO1,
        sqlite_db       = lambda wc: tmp_db(wc.filename),
        rule_name       = "marker_id",
        tmpdir          = lambda wc: tmp_dir(wc.filename),
        logger          = f"{script_dir}/wagtail_sqlite_logger.py",
        annotation      = f"{script_dir}/error_annotation.py",
        n_subsample     = 1000,
        min_confidence  = 0.6,
        bash_lib        = BASH_LIB
    conda:
        "envs/mappy.yaml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/marker_id.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/marker_id.benchmark.log"
    threads: 1
    resources:
        mem_mb  = 24000,
        runtime = "1h"
    shell:
        r"""
        source {params.bash_lib}
        LOGGER={params.logger}
        mkdir -p $(dirname {output.marker}) {params.tmpdir}

        set +e
        python {params.script} \
            --file-map  {input.filemap} \
            --sample    {params.sample} \
            --db-16S    {params.db_16S} \
            --db-ITS    {params.db_ITS} \
            --db-CO1    {params.db_CO1} \
            --db-18S    {params.db_18S} \
            --n-subsample   {params.n_subsample} \
            --min-confidence {params.min_confidence} \
            --output    {output.marker} \
            &> {log}
        status=$?
        set -e

        if [[ $status -ne 0 ]]; then
            touch {output.marker}
            wagtail_log_fail {params.sqlite_db} {wildcards.filename} {params.rule_name} {log} {params.annotation}
        else
            marker_status=$(tail -n 1 {output.marker} | cut -f1)
            predicted=$(tail -n 1 {output.marker} | cut -f2)
            second=$(tail -n 1 {output.marker} | cut -f5)
            best_id=$(tail -n 1 {output.marker} | cut -f3)
            second_id=$(tail -n 1 {output.marker} | cut -f6)
            margin=$(tail -n 1 {output.marker} | cut -f8)

            if [[ "$marker_status" == "AMBIGUOUS" ]]; then
                echo "ERROR: {wildcards.filename} marker ambiguous between best ($predicted mean_id=$best_id) and second-best ($second mean_id=$second_id); margin=$margin" >> {log}
                wagtail_log_fail {params.sqlite_db} {wildcards.filename} {params.rule_name} {log} {params.annotation}
            elif [[ "$marker_status" == "UNKNOWN" ]]; then
                echo "ERROR: {wildcards.filename} amplicon is not 16S, 18S, ITS, or CO1 (best hit fraction below confidence threshold)" >> {log}
                wagtail_log_fail {params.sqlite_db} {wildcards.filename} {params.rule_name} {log} {params.annotation}
            else
                wagtail_log_ok {params.sqlite_db} {wildcards.filename} {params.rule_name} {log}
            fi
        fi
        wagtail_cleanup {log}
        exit 0
        """

# ── manifest ──────────────────────────────────────────────────────────────────
rule manifest:
    group: "wagtail"
    input:
        filemap = config["file_map"],
        marker  = f"{run_dir}/{run_name}/0_marker_id/{{filename}}.marker"
    output:
        output_dir = directory(f"{run_dir}/{run_name}/0_manifest/{{filename}}"),
        manifest   = f"{run_dir}/{run_name}/0_manifest/{{filename}}/{{filename}}_manifest.csv"
    wildcard_constraints:
        filename = r"[^\.]+"
    params:
        script    = f"{script_dir}/create_manifest_wagtail.py",
        sample    = "{filename}",
        sqlite_db = lambda wc: tmp_db(wc.filename),
        rule_name = "manifest",
        upstream  = "marker_id",
        tmpdir    = lambda wc: tmp_dir(wc.filename),
        logger    = f"{script_dir}/wagtail_sqlite_logger.py",
        checker   = f"{script_dir}/check_sqlite_status.py",
        annotation= f"{script_dir}/error_annotation.py",
        bash_lib  = BASH_LIB
    conda:
        "envs/mappy.yaml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/manifest_creation.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/manifest_creation.benchmark.log"
    threads: 1
    resources:
        mem_mb  = 320,
        runtime = "1h"
    shell:
        r"""
        source {params.bash_lib}
        LOGGER={params.logger}
        CHECKER={params.checker}
        mkdir -p {params.tmpdir}
        wagtail_check_upstream {params.sqlite_db} {wildcards.filename} {params.upstream} \
            "{output.manifest}" \
            && {{ mkdir -p {output.output_dir}; touch {output.manifest}; }} \
            && wagtail_log_fail {params.sqlite_db} {wildcards.filename} {params.rule_name} {log} {params.annotation} \
            && wagtail_cleanup {log} && exit 0

        set +e
        python {params.script} \
            --input    {params.sample} \
            --file-map {input.filemap} \
            --output   {output.output_dir} \
            &> {log}
        status=$?
        set -e

        if [[ $status -ne 0 ]]; then
            mkdir -p {output.output_dir}; touch {output.manifest}
            wagtail_log_fail {params.sqlite_db} {wildcards.filename} {params.rule_name} {log} {params.annotation}
        else
            wagtail_log_ok {params.sqlite_db} {wildcards.filename} {params.rule_name} {log}
        fi
        wagtail_cleanup {log}
        exit 0
        """

# ── qiime2_import ─────────────────────────────────────────────────────────────
rule qiime2_import:
    group: "wagtail"
    input:
        manifest = f"{run_dir}/{run_name}/0_manifest/{{filename}}/{{filename}}_manifest.csv"
    output:
        qza = f"{run_dir}/{run_name}/1_import_wagtail/{{filename}}-single.qza"
    wildcard_constraints:
        filename = r"[^\.]+"
    params:
        sqlite_db = lambda wc: tmp_db(wc.filename),
        rule_name = "qiime2_import",
        upstream  = "manifest",
        tmpdir    = lambda wc: tmp_dir(wc.filename),
        logger    = f"{script_dir}/wagtail_sqlite_logger.py",
        checker   = f"{script_dir}/check_sqlite_status.py",
        annotation= f"{script_dir}/error_annotation.py",
        bash_lib  = BASH_LIB
    conda:
        "envs/qiime2.yml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/qiime2_import.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/qiime2_import.benchmark.log"
    threads: 1
    resources:
        mem_mb  = 640,
        runtime = "1h"
    shell:
        r"""
        source {params.bash_lib}
        LOGGER={params.logger}
        CHECKER={params.checker}
        wagtail_check_upstream {params.sqlite_db} {wildcards.filename} {params.upstream} \
            "{output.qza}" \
            && touch {output.qza} \
            && wagtail_log_fail {params.sqlite_db} {wildcards.filename} {params.rule_name} {log} {params.annotation} \
            && wagtail_cleanup {log} && exit 0

        export TMPDIR={params.tmpdir}; mkdir -p $TMPDIR
        set +e
        qiime tools import \
            --type 'SampleData[SequencesWithQuality]' \
            --input-path   {input.manifest} \
            --input-format SingleEndFastqManifestPhred33 \
            --output-path  {output.qza} \
            &> {log}
        status=$?
        set -e

        wagtail_handle_status $status {params.sqlite_db} {wildcards.filename} {params.rule_name} {log} {params.annotation} \
            "{output.qza}"
        wagtail_cleanup {log}
        exit 0
        """

# ── quality_control ───────────────────────────────────────────────────────────
rule quality_control:
    group: "wagtail"
    input:
        qza = f"{run_dir}/{run_name}/1_import_wagtail/{{filename}}-single.qza"
    output:
        filtered = f"{run_dir}/{run_name}/2_qc_wagtail/{{filename}}-filtered.qza",
        stats    = f"{run_dir}/{run_name}/2_qc_wagtail/{{filename}}-qc-stats.qza"
    wildcard_constraints:
        filename = r"[^\.]+"
    params:
        sqlite_db = lambda wc: tmp_db(wc.filename),
        rule_name = "quality_control",
        upstream  = "qiime2_import",
        tmpdir    = lambda wc: tmp_dir(wc.filename),
        logger    = f"{script_dir}/wagtail_sqlite_logger.py",
        checker   = f"{script_dir}/check_sqlite_status.py",
        annotation= f"{script_dir}/error_annotation.py",
        bash_lib  = BASH_LIB
    conda:
        "envs/qiime2.yml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/quality_control.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/quality_control.benchmark.log"
    threads: 1
    resources:
        mem_mb  = 2000,
        runtime = "1h"
    shell:
        r"""
        source {params.bash_lib}
        LOGGER={params.logger}
        CHECKER={params.checker}
        wagtail_check_upstream {params.sqlite_db} {wildcards.filename} {params.upstream} \
            "{output.filtered} {output.stats}" \
            && touch {output.filtered} {output.stats} \
            && wagtail_log_fail {params.sqlite_db} {wildcards.filename} {params.rule_name} {log} {params.annotation} \
            && wagtail_cleanup {log} && exit 0

        export TMPDIR={params.tmpdir}
        set +e
        qiime quality-filter q-score \
            --i-demux            {input.qza} \
            --o-filtered-sequences {output.filtered} \
            --o-filter-stats     {output.stats} \
            &> {log}
        status=$?
        set -e

        wagtail_handle_status $status {params.sqlite_db} {wildcards.filename} {params.rule_name} {log} {params.annotation} \
            "{output.filtered} {output.stats}"
        wagtail_cleanup {log}
        exit 0
        """

# ── deblur ────────────────────────────────────────────────────────────────────
rule deblur:
    group: "wagtail"
    input:
        filtered    = f"{run_dir}/{run_name}/2_qc_wagtail/{{filename}}-filtered.qza",
        reference   = lambda wc: get_deblur_db(wc),
    output:
        representative = f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}-rep-seqs.qza",
        table          = f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}-table.qza",
        stats          = f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}-deblur-stats.qza"
    wildcard_constraints:
        filename = r"[^\.]+"
    params:
        path      = directory(f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}"),
        script    = f"{script_dir}/deblur_all.py",
        marker    = lambda wc: get_marker(wc),
        sqlite_db = lambda wc: tmp_db(wc.filename),
        rule_name = "deblur",
        upstream  = "quality_control",
        tmpdir    = lambda wc: tmp_dir(wc.filename),
        logger    = f"{script_dir}/wagtail_sqlite_logger.py",
        checker   = f"{script_dir}/check_sqlite_status.py",
        annotation= f"{script_dir}/error_annotation.py",
        bash_lib  = BASH_LIB
    conda:
        "envs/qiime2.yml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/deblur.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/deblur.benchmark.log"
    threads: 1
    resources:
        mem_mb  = 2000,
        runtime = "47h"
    shell:
        r"""
        source {params.bash_lib}
        LOGGER={params.logger}
        CHECKER={params.checker}
        wagtail_check_upstream {params.sqlite_db} {wildcards.filename} {params.upstream} \
            "{output.representative} {output.table} {output.stats}" \
            && mkdir -p {params.path} \
            && touch {output.representative} {output.table} {output.stats} \
            && wagtail_log_fail {params.sqlite_db} {wildcards.filename} {params.rule_name} {log} {params.annotation} \
            && wagtail_cleanup {log} && exit 0

        export TMPDIR={params.tmpdir}
        set +e
        python {params.script} \
            -i {input.filtered} -o {params.path} -t {threads} \
            -r {output.representative} -a {output.table} -s {output.stats} \
            --marker    {params.marker} \
            --reference {input.reference} \
            &> {log}
        status=$?
        set -e

        if [[ $status -ne 0 ]]; then
            mkdir -p {params.path}
            wagtail_log_fail {params.sqlite_db} {wildcards.filename} {params.rule_name} {log} {params.annotation}
            touch {output.representative} {output.table} {output.stats}
        else
            wagtail_log_ok {params.sqlite_db} {wildcards.filename} {params.rule_name} {log}
        fi
        wagtail_cleanup {log}
        exit 0
        """

# ── export_seqs ───────────────────────────────────────────────────────────────
rule export_seqs:
    group: "wagtail"
    input:
        representative = f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}-rep-seqs.qza"
    output:
        output = directory(f"{run_dir}/{run_name}/4_rep_seqs_wagtail/{{filename}}"),
        final  = f"{run_dir}/{run_name}/4_rep_seqs_wagtail/{{filename}}/dna-sequences.fasta"
    wildcard_constraints:
        filename = r"[^\.]+"
    params:
        script    = f"{script_dir}/qiime2_seqs_export.py",
        sqlite_db = lambda wc: tmp_db(wc.filename),
        rule_name = "export_seqs",
        upstream  = "deblur",
        tmpdir    = lambda wc: tmp_dir(wc.filename),
        logger    = f"{script_dir}/wagtail_sqlite_logger.py",
        checker   = f"{script_dir}/check_sqlite_status.py",
        annotation= f"{script_dir}/error_annotation.py",
        bash_lib  = BASH_LIB
    conda:
        "envs/qiime2.yml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/export_seqs.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/export_seqs.benchmark.log"
    threads: 1
    resources:
        mem_mb  = 2000,
        runtime = "1h"
    shell:
        r"""
        source {params.bash_lib}
        LOGGER={params.logger}
        CHECKER={params.checker}
        wagtail_check_upstream {params.sqlite_db} {wildcards.filename} {params.upstream} \
            "{output.final}" \
            && mkdir -p {output.output} && touch {output.final} \
            && wagtail_log_fail {params.sqlite_db} {wildcards.filename} {params.rule_name} {log} {params.annotation} \
            && wagtail_cleanup {log} && exit 0

        export TMPDIR={params.tmpdir}
        set +e
        python {params.script} \
            --input-path  {input.representative} \
            --output-path {output.output} \
            &> {log}
        status=$?
        set -e

        wagtail_handle_status $status {params.sqlite_db} {wildcards.filename} {params.rule_name} {log} {params.annotation} \
            "{output.final}"
        wagtail_cleanup {log}
        exit 0
        """

# ── mappy ─────────────────────────────────────────────────────────────────────
rule mappy:
    group: "wagtail"
    input:
        F   = f"{run_dir}/{run_name}/4_rep_seqs_wagtail/{{filename}}/dna-sequences.fasta",
        ref = lambda wc: get_mappy_ref(wc),
    output:
        align = f"{run_dir}/{run_name}/5_taxonomy_wagtail/{{filename}}/{{filename}}_alignment.tsv",
        meta  = f"{run_dir}/{run_name}/5_taxonomy_wagtail/{{filename}}/{{filename}}_metadata.tsv"
    wildcard_constraints:
        filename = r"[^\.]+"
    params:
        script    = f"{script_dir}/mappy_script.py",
        sample    = "{filename}",
        sqlite_db = lambda wc: tmp_db(wc.filename),
        rule_name = "mappy",
        upstream  = "export_seqs",
        logger    = f"{script_dir}/wagtail_sqlite_logger.py",
        checker   = f"{script_dir}/check_sqlite_status.py",
        annotation= f"{script_dir}/error_annotation.py",
        bash_lib  = BASH_LIB
    conda:
        "envs/mappy.yaml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/mappy.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/mappy.benchmark.log"
    threads: 1
    resources:
        mem_mb  = 4000,
        runtime = "24h"
    shell:
        r"""
        source {params.bash_lib}
        LOGGER={params.logger}
        CHECKER={params.checker}
        wagtail_check_upstream {params.sqlite_db} {wildcards.filename} {params.upstream} \
            "{output.align} {output.meta}" \
            && touch {output.align} {output.meta} \
            && wagtail_log_fail {params.sqlite_db} {wildcards.filename} {params.rule_name} {log} {params.annotation} \
            && wagtail_cleanup {log} && exit 0

        set +e
        python {params.script} \
            -i {input.F} -r {input.ref} \
            -a {output.align} -m {output.meta} -s {params.sample} \
            &> {log}
        status=$?
        set -e

        wagtail_handle_status $status {params.sqlite_db} {wildcards.filename} {params.rule_name} {log} {params.annotation} \
            "{output.align} {output.meta}"
        wagtail_cleanup {log}
        exit 0
        """

# ── export_and_edit_table ─────────────────────────────────────────────────────
rule export_and_edit_table:
    group: "wagtail"
    input:
        table = f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}-table.qza"
    output:
        folder = directory(f"{run_dir}/{run_name}/4_table_wagtail/{{filename}}_biom"),
        biom   = f"{run_dir}/{run_name}/4_table_wagtail/{{filename}}_biom/{{filename}}.biom",
        table  = f"{run_dir}/{run_name}/4_table_wagtail/{{filename}}_table.tsv",
        edited = f"{run_dir}/{run_name}/5_taxonomy_wagtail/{{filename}}/{{filename}}_table_edit.tsv"
    wildcard_constraints:
        filename = r"[^\.]+"
    params:
        script    = f"{script_dir}/qiime2_biom_export.py",
        name      = "{filename}.biom",
        sqlite_db = lambda wc: tmp_db(wc.filename),
        rule_name = "export_and_edit_table",
        upstream  = "deblur",
        tmpdir    = lambda wc: tmp_dir(wc.filename),
        logger    = f"{script_dir}/wagtail_sqlite_logger.py",
        checker   = f"{script_dir}/check_sqlite_status.py",
        annotation= f"{script_dir}/error_annotation.py",
        bash_lib  = BASH_LIB
    conda:
        "envs/qiime2.yml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/export_and_edit_table.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/export_and_edit_table.benchmark.log"
    threads: 1
    resources:
        mem_mb  = 1000,
        runtime = "1h"
    shell:
        r"""
        source {params.bash_lib}
        LOGGER={params.logger}
        CHECKER={params.checker}
        wagtail_check_upstream {params.sqlite_db} {wildcards.filename} {params.upstream} \
            "{output.biom} {output.table} {output.edited}" \
            && mkdir -p {output.folder} && touch {output.biom} {output.table} {output.edited} \
            && wagtail_log_fail {params.sqlite_db} {wildcards.filename} {params.rule_name} {log} {params.annotation} \
            && wagtail_cleanup {log} && exit 0

        export TMPDIR={params.tmpdir}
        set +e
        python {params.script} \
            --input-path   {input.table} \
            --output-path  {output.folder} \
            --new-filename {params.name} \
            && biom convert -i {output.biom} -o {output.table} --to-tsv \
            && tail -n +3 {output.table} > {output.edited} \
            &> {log}
        status=$?
        set -e

        if [[ $status -ne 0 ]]; then
            mkdir -p {output.folder}
            wagtail_log_fail {params.sqlite_db} {wildcards.filename} {params.rule_name} {log} {params.annotation}
            touch {output.biom} {output.table} {output.edited}
        else
            wagtail_log_ok {params.sqlite_db} {wildcards.filename} {params.rule_name} {log}
        fi
        wagtail_cleanup {log}
        exit 0
        """

# ── extract_taxonomy ──────────────────────────────────────────────────────────
rule extract_taxonomy:
    group: "wagtail"
    input:
        primary     = f"{run_dir}/{run_name}/5_taxonomy_wagtail/{{filename}}/{{filename}}_alignment.tsv",
        table       = f"{run_dir}/{run_name}/5_taxonomy_wagtail/{{filename}}/{{filename}}_table_edit.tsv",
        reference   = lambda wc: get_tax_ref(wc),
    output:
        condensed = f"{run_dir}/{run_name}/6_condensed_wagtail/{{filename}}_condensed.tsv"
    wildcard_constraints:
        filename = r"[^\.]+"
    params:
        script    = f"{script_dir}/extract_taxonomy.py",
        sample    = "{filename}",
        sqlite_db = lambda wc: tmp_db(wc.filename),
        rule_name = "extract_taxonomy",
        upstream  = "mappy",
        logger    = f"{script_dir}/wagtail_sqlite_logger.py",
        checker   = f"{script_dir}/check_sqlite_status.py",
        annotation= f"{script_dir}/error_annotation.py",
        bash_lib  = BASH_LIB
    conda:
        "envs/mappy.yaml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/extract_taxonomy.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/extract_taxonomy.benchmark.log"
    threads: 1
    resources:
        mem_mb  = 160,
        runtime = "1h"
    shell:
        r"""
        source {params.bash_lib}
        LOGGER={params.logger}
        CHECKER={params.checker}
        wagtail_check_upstream {params.sqlite_db} {wildcards.filename} {params.upstream} \
            "{output.condensed}" \
            && touch {output.condensed} \
            && wagtail_log_fail {params.sqlite_db} {wildcards.filename} {params.rule_name} {log} {params.annotation} \
            && wagtail_cleanup {log} && exit 0

        set +e
        python {params.script} \
            -i {input.primary} -t {input.table} -r {input.reference} \
            -o {output.condensed} -s {params.sample} \
            &> {log}
        status=$?
        set -e

        wagtail_handle_status $status {params.sqlite_db} {wildcards.filename} {params.rule_name} {log} {params.annotation} \
            "{output.condensed}"
        wagtail_cleanup {log}
        exit 0
        """

# ── qiime_stats ───────────────────────────────────────────────────────────────
rule qiime_stats:
    group: "wagtail"
    input:
        qc     = f"{run_dir}/{run_name}/2_qc_wagtail/{{filename}}-qc-stats.qza",
        deblur = f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}-deblur-stats.qza"
    output:
        qc_dir      = directory(f"{run_dir}/{run_name}/2_qc_wagtail/{{filename}}"),
        deblur_dir  = directory(f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}"),
        final_qc    = f"{run_dir}/{run_name}/2_qc_wagtail/{{filename}}/stats.csv",
        final_deblur= f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}/stats.csv"
    wildcard_constraints:
        filename = r"[^\.]+"
    params:
        script    = f"{script_dir}/wagtail_metadata_qiimes.py",
        sqlite_db = lambda wc: tmp_db(wc.filename),
        rule_name = "qiime_stats",
        upstream  = "deblur",
        tmpdir    = lambda wc: tmp_dir(wc.filename),
        logger    = f"{script_dir}/wagtail_sqlite_logger.py",
        checker   = f"{script_dir}/check_sqlite_status.py",
        annotation= f"{script_dir}/error_annotation.py",
        bash_lib  = BASH_LIB
    conda:
        "envs/qiime2.yml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/qiime_stats.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/qiime_stats.benchmark.log"
    threads: 1
    resources:
        mem_mb  = 640,
        runtime = "1h"
    shell:
        r"""
        source {params.bash_lib}
        LOGGER={params.logger}
        CHECKER={params.checker}
        wagtail_check_upstream {params.sqlite_db} {wildcards.filename} {params.upstream} \
            "{output.final_qc} {output.final_deblur}" \
            && mkdir -p {output.qc_dir} {output.deblur_dir} \
            && touch {output.final_qc} {output.final_deblur} \
            && wagtail_log_fail {params.sqlite_db} {wildcards.filename} {params.rule_name} {log} {params.annotation} \
            && wagtail_cleanup {log} && exit 0

        # Also guard against empty upstream outputs
        if [[ ! -s {input.qc} ]] || [[ ! -s {input.deblur} ]]; then
            mkdir -p {output.qc_dir} {output.deblur_dir}
            touch {output.final_qc} {output.final_deblur}
            echo "Upstream fails" > {log}
            wagtail_log_fail {params.sqlite_db} {wildcards.filename} {params.rule_name} {log} {params.annotation}
            wagtail_cleanup {log}
            exit 0
        fi

        export TMPDIR={params.tmpdir}
        set +e
        python {params.script} \
            --input-qc     {input.qc}     --output-qc     {output.qc_dir} \
            --input-deblur {input.deblur} --output-deblur {output.deblur_dir} \
            &> {log}
        status=$?
        set -e

        if [[ $status -ne 0 ]]; then
            mkdir -p {output.qc_dir} {output.deblur_dir}
            touch {output.final_qc} {output.final_deblur}
            wagtail_log_fail {params.sqlite_db} {wildcards.filename} {params.rule_name} {log} {params.annotation}
        else
            wagtail_log_ok {params.sqlite_db} {wildcards.filename} {params.rule_name} {log}
        fi
        wagtail_cleanup {log}
        exit 0
        """

# ── metadata_creation ─────────────────────────────────────────────────────────
rule metadata_creation:
    group: "wagtail"
    input:
        final_qc    = f"{run_dir}/{run_name}/2_qc_wagtail/{{filename}}/stats.csv",
        final_deblur= f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}/stats.csv",
        final_mappy = f"{run_dir}/{run_name}/5_taxonomy_wagtail/{{filename}}/{{filename}}_metadata.tsv"
    output:
        directory = directory(f"{run_dir}/{run_name}/7_metadata/{{filename}}"),
        final     = f"{run_dir}/{run_name}/7_metadata/{{filename}}/{{filename}}_metadata.tsv"
    wildcard_constraints:
        filename = r"[^\.]+"
    params:
        script    = f"{script_dir}/wagtail_metadata_meta_combine.py",
        sqlite_db = lambda wc: tmp_db(wc.filename),
        run       = "{filename}",
        rule_name = "metadata_creation",
        upstream  = "qiime_stats",
        tmpdir    = lambda wc: tmp_dir(wc.filename),
        logger    = f"{script_dir}/wagtail_sqlite_logger.py",
        checker   = f"{script_dir}/check_sqlite_status.py",
        annotation= f"{script_dir}/error_annotation.py",
        bash_lib  = BASH_LIB
    conda:
        "envs/mappy.yaml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/metadata.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/metadata.benchmark.log"
    threads: 1
    resources:
        mem_mb  = 320,
        runtime = "1h"
    shell:
        r"""
        source {params.bash_lib}
        LOGGER={params.logger}
        CHECKER={params.checker}
        wagtail_check_upstream {params.sqlite_db} {wildcards.filename} {params.upstream} \
            "{output.final}" \
            && mkdir -p {output.directory} && touch {output.final} \
            && wagtail_log_fail {params.sqlite_db} {wildcards.filename} {params.rule_name} {log} {params.annotation} \
            && wagtail_cleanup {log} && exit 0

        set +e
        python {params.script} \
            --qc-file     {input.final_qc} \
            --deblur-file {input.final_deblur} \
            --mappy-file  {input.final_mappy} \
            --output-path {output.directory} \
            --run-name    {params.run} \
            &> {log}
        status=$?
        set -e

        if [[ $status -ne 0 ]]; then
            mkdir -p {output.directory}
            echo "ERROR: metadata_creation failed for {wildcards.filename}" > {output.final}
            wagtail_log_fail {params.sqlite_db} {wildcards.filename} {params.rule_name} {log} {params.annotation}
        else
            wagtail_log_ok {params.sqlite_db} {wildcards.filename} {params.rule_name} {log}
        fi
        wagtail_cleanup {log}
        exit 0
        """

# ── cleanup ───────────────────────────────────────────────────────────────────
rule cleanup:
    group: "cleanup"
    input:
        metadata  = expand(f"{run_dir}/{run_name}/7_metadata/{{filename}}/{{filename}}_metadata.tsv", filename=filenames),
        community = expand(f"{run_dir}/{run_name}/6_condensed_wagtail/{{filename}}_condensed.tsv",    filename=filenames)
    output:
        merged        = sql_log,
        cleanup_done  = f"{run_dir}/{run_name}/0_logs_wagtail/cleanup_done.txt",
        full_metadata = f"{run_dir}/{run_name}/7_metadata/{run_name}_full_metadata.tsv"
    wildcard_constraints:
        filename = r"[^\.]+"
    params:
        merge_script = f"{script_dir}/sql_merge.py",
        dbs          = expand(f"{run_dir}/{run_name}/0_tmp/{{filename}}_log.sql", filename=filenames),
        target       = lambda w, output: str(Path(output.merged).parent.parent),
        tmpdir       = lambda w, output: str(Path(output.merged).parent.parent / "0_tmp")
    conda:
        "envs/mappy.yaml"
    log:
        merger       = f"{run_dir}/{run_name}/0_logs_wagtail/merger.log",
        cleaning     = f"{run_dir}/{run_name}/0_logs_wagtail/cleaning.log",
        full_metadata= f"{run_dir}/{run_name}/0_logs_wagtail/full_metadata.log"
    threads: 1
    resources:
        mem_mb  = 2000,
        runtime = "2h"
    shell:
        r"""
        # Step 1: Merge per-sample SQLite logs
        touch {input.metadata[0]} {input.community[0]}
        python {params.merge_script} {output.merged} {params.dbs} &> {log.merger}
        sleep 2

        # Step 2: Clean intermediate directories
        rm -rf {params.tmpdir}
        rm -rf {params.target}/[1-3]_* {params.target}/5_*
        touch {output.cleanup_done}
        echo "Cleanup completed" &> {log.cleaning}

        # Step 3: Merge per-sample metadata into one file
        (head -n 1 {input.metadata[0]} > {output.full_metadata} \
            && tail -n +2 -q {input.metadata} >> {output.full_metadata}) \
            &> {log.full_metadata}
        """
