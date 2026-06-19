# Wagtail pipeline v1.2
# Changes from v1.1:
#   - marker_id is now a single BATCHED plain RULE over the whole sample list
#     (no longer a checkpoint): identify_amplicon.py loads each database once and
#     loops every sample, so the multi-GB minimap2 index is built once per batch
#     instead of once per sample. Writes all 1_marker_id/{sample}_marker.json and
#     touches ONE marker_id.done sentinel — not a 5000-file `expand` output —
#     which keeps the cluster submit command under ARG_MAX at 5k+ samples.
#   - reference resolution moved from DAG-time Snakemake helpers to RUNTIME:
#     deblur/mappy/extract_taxonomy each call bin/resolve_reference.py to read
#     their sample's 'reference' key and glob the database dir. This removes the
#     checkpoint dependency that the plain rule can no longer satisfy.
#   - new lightweight per-sample marker_log rule reads each _marker.json (via the
#     sentinel gate) and records OK/FAIL into SQLite under rule_name "marker_id"
#     (manifest's upstream check is unchanged). Grouped with "wagtail".
#   - manifest depends on 1_marker_logged/{sample}.logged instead of the JSON.
#   - the all-samples fan-in is collapsed by a `localrule gather_data` into one
#     `all_samples.done` sentinel; cleanup depends on that single file, not the
#     per-sample outputs. Supports optional result/intermediate/all archiving via
#     wagtail.py --archive-* flags.

import os
import glob
import csv
import json
from functools import lru_cache
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

# Amplicon filter ─────────────────────────────────────────────────────────────
# "all"  (default) — auto-detect per sample via the marker_id checkpoint.
# "16S" | "18S" | "ITS" | "CO1" — skip identification; treat every sample as
# chosen amplicon type. Pass via wagtail.py --amplicon or --config amplicon=ITS.
_AMPLICON_RAW     = config.get("amplicon", "all")
FORCED_AMPLICON   = _AMPLICON_RAW.upper()
_VALID_AMPLICONS  = {"ALL", "16S", "18S", "ITS", "CO1"}
if FORCED_AMPLICON not in _VALID_AMPLICONS:
    raise WorkflowError(
        f"config amplicon={_AMPLICON_RAW!r} is not valid. "
        f"Accepted values: {', '.join(sorted(_VALID_AMPLICONS - {'ALL'}))} "
        f"— or omit to auto-detect."
    )
AUTO_DETECT = FORCED_AMPLICON == "ALL"   # convenience flag

# ── Archive options (set from wagtail.py --archive-* flags) ───────────────────
def _truthy(v) -> bool:
    return str(v).strip().lower() in {"1", "true", "yes", "on"}

ARCHIVE_RESULT       = _truthy(config.get("archive_result", ""))   # zip 0_/6_/7_
ARCHIVE_ALL          = _truthy(config.get("archive_all", ""))      # zip 0_..7_
ARCHIVE_INTERMEDIATE = str(config.get("archive_intermediate", "")).strip()  # "" | "all" | "1,3,5"

def _archive_cli() -> str:
    """Build the wagtail_archive.py flag string for the configured options. Empty = no archiving."""
    parts = []
    if ARCHIVE_ALL:
        parts.append("--archive-all")
    if ARCHIVE_RESULT:
        parts.append("--archive-result")
    if ARCHIVE_INTERMEDIATE:
        parts.append(f"--archive-intermediate {ARCHIVE_INTERMEDIATE}")
    return " ".join(parts)

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

# ── Marker / reference helpers ────────────────────────────────────────────────
# marker_id writes one _marker.json per sample plus a single batch sentinel (marker_done).
# Single batch sentinel produced by rule marker_id (avoids huge sample to overflow ARG_MAX when the job is submitted to the cluster).
marker_done = f"{run_dir}/{run_name}/1_marker_id/marker_id.done"

def marker_json(filename: str) -> str:
    """Path to a sample's _marker.json (read at runtime by resolve_reference.py)."""
    return f"{run_dir}/{run_name}/1_marker_id/{filename}_marker.json"

# Runtime reference resolver — called from deblur / mappy / extract_taxonomy.
RESOLVER = f"{script_dir}/resolve_reference.py"

@lru_cache(maxsize=None)
def _resolve(marker: str, db_dir: str, pattern: str) -> str:
    """Resolve a single file matching db_dir/pattern. Raises if 0 or >1 match.
    Used at parse time for the identification databases below."""
    matches = glob.glob(f"{db_dir}/{pattern}")
    if not matches:
        raise FileNotFoundError(f"No file matching '{pattern}' found in {db_dir}")
    if len(matches) > 1:
        raise ValueError(f"Multiple files matching '{pattern}' in {db_dir}: {matches}")
    return matches[0]

# Identification databases are only needed for auto-detection.
# In forced-amplicon mode these are never passed to identify_amplicon.py,
DB_16S = _resolve("16S", db_dir, pattern="16S_database*") if AUTO_DETECT else ""
DB_18S = _resolve("18S", db_dir, pattern="18S_database*") if AUTO_DETECT else ""
DB_ITS = _resolve("ITS", db_dir, pattern="ITS_database*") if AUTO_DETECT else ""
DB_CO1 = _resolve("CO1", db_dir, pattern="CO1_database*") if AUTO_DETECT else ""

# ── rule all ─────────────────────────────────────────────────────────────────
rule all:
    input:
        expand(f"{run_dir}/{run_name}/6_condensed_wagtail/{{filename}}_condensed.tsv", filename=filenames),
        f"{run_dir}/{run_name}/7_metadata/{run_name}_full_metadata.tsv",
        sql_log
    localrule: True

# ── marker_id (batched) ───────────────────────────────────────────────────────
# Pre-flight QC for the WHOLE sample list in a single job. identify_amplicon.py
# loads each of the 4 databases once and loops over every sample. 
# Writes one _marker.json per sample into 1_marker_id/, then touches `marker_id.done`.

rule marker_id:
    localrule: not AUTO_DETECT
    input:
        filemap     = config["file_map"],
        sample_list = filename_list,
    output:
        sentinel = marker_done,
        summary  = f"{run_dir}/{run_name}/0_logs_wagtail/{run_name}_marker_summary.tsv"
    params:
        script          = f"{script_dir}/identify_amplicon.py",
        forced_script   = f"{script_dir}/write_forced_markers.py",
        summary_script  = f"{script_dir}/marker_summary.py",
        marker_dir      = f"{run_dir}/{run_name}/1_marker_id",
        db_16S          = DB_16S,
        db_18S          = DB_18S,
        db_ITS          = DB_ITS,
        db_CO1          = DB_CO1,
        forced_amplicon = FORCED_AMPLICON,    # "ALL" → auto-detect; else forced
        n_subsample     = 1000,
        min_confidence  = 0.6,
        minimizer_w     = 19,
    conda:
        "envs/mappy.yaml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/marker_id_batch.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/marker_id_batch.benchmark.log"
    threads: 1
    resources:
        # One index resident at a time (~6.5 GB at w=19) + all subsampled reads.
        mem_mb  = 320 if not AUTO_DETECT else 10000,
        # One job for the whole batch — size batches so the auto path fits in 47h.
        runtime = "1h" if not AUTO_DETECT else "47h"
    shell:
        r"""
        mkdir -p {params.marker_dir} $(dirname {log}) $(dirname {output.summary})

        # ── Forced-amplicon fast path ─────────────────────────────────────────
        # No alignment needed — write a trivial OK marker JSON for every sample.
        if [[ "{params.forced_amplicon}" != "ALL" ]]; then
            python {params.forced_script} \
                --sample-list {input.sample_list} \
                --output-dir  {params.marker_dir} \
                --marker      {params.forced_amplicon} \
                &> {log}
            # Summary: every sample reported FORCED/<marker>.
            python {params.summary_script} \
                --marker-dir  {params.marker_dir} \
                --sample-list {input.sample_list} \
                --output      {output.summary} \
                --forced      {params.forced_amplicon} \
                &>> {log}
            touch {output.sentinel}
            exit 0
        fi

        # ── Auto-detection: load each DB once, loop all samples ───────────────
        python {params.script} \
            --sample-list {input.sample_list} \
            --file-map    {input.filemap} \
            --db-16S {params.db_16S} --db-ITS {params.db_ITS} \
            --db-18S {params.db_18S} --db-CO1 {params.db_CO1} \
            --n-subsample {params.n_subsample} \
            --min-confidence {params.min_confidence} \
            --minimizer-w {params.minimizer_w} \
            --output-dir {params.marker_dir} \
            &> {log}
        # Summary: OK/<marker> for confident calls, FAIL/UNKNOWN otherwise.
        python {params.summary_script} \
            --marker-dir  {params.marker_dir} \
            --sample-list {input.sample_list} \
            --output      {output.summary} \
            &>> {log}
        touch {output.sentinel}
        """

# ── marker_log ────────────────────────────────────────────────────────────────
# Lightweight per-sample step: read this sample's _marker.json (produced by the
# marker_id rule) and record OK/FAIL into its SQLite log under rule_name
# "marker_id", so manifest's existing upstream check is unchanged.

rule marker_log:
    group: "wagtail"
    input:
        sentinel = marker_done
    output:
        logged = f"{run_dir}/{run_name}/1_marker_logged/{{filename}}.logged"
    wildcard_constraints:
        filename = r"[^\.]+"
    params:
        marker    = lambda wc: marker_json(wc.filename),
        reader    = f"{script_dir}/read_marker_fields.py",
        sqlite_db = lambda wc: tmp_db(wc.filename),
        rule_name = "marker_id",     # log under marker_id so manifest's check matches
        logger    = f"{script_dir}/wagtail_sqlite_logger.py",
        annotation= f"{script_dir}/error_annotation.py",
        bash_lib  = BASH_LIB
    conda:
        "envs/mappy.yaml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/marker_log.log"
    threads: 1
    resources:
        mem_mb  = 160,
        runtime = "1h"
    shell:
        r"""
        source {params.bash_lib}
        LOGGER={params.logger}
        mkdir -p $(dirname {output.logged}) $(dirname {log}) $(dirname {params.sqlite_db})

        # Empty/missing marker JSON = failed sample (no reads, missing file) → soft-fail
        if [[ ! -s {params.marker} ]]; then
            echo "ERROR: empty/missing marker file for {wildcards.filename}" > {log}
            wagtail_log_fail {params.sqlite_db} {wildcards.filename} {params.rule_name} {log} {params.annotation}
            touch {output.logged}
            wagtail_cleanup {log}
            exit 0
        fi

        # Parse the marker JSON via a helper script (6 fields, one per line).
        _json_vals=$(python {params.reader} {params.marker})
        marker_status=$(echo "$_json_vals" | sed -n '1p')
        predicted=$(echo "$_json_vals"     | sed -n '2p')
        best_id=$(echo "$_json_vals"       | sed -n '3p')
        second=$(echo "$_json_vals"        | sed -n '4p')
        second_id=$(echo "$_json_vals"     | sed -n '5p')
        margin=$(echo "$_json_vals"        | sed -n '6p')

        if [[ "$marker_status" == "OK" ]]; then
            wagtail_log_ok {params.sqlite_db} {wildcards.filename} {params.rule_name} {log}
        elif [[ "$marker_status" == "AMBIGUOUS" ]]; then
            echo "ERROR: {wildcards.filename} marker ambiguous between best ($predicted mean_id=$best_id) and second-best ($second mean_id=$second_id); margin=$margin" > {log}
            wagtail_log_fail {params.sqlite_db} {wildcards.filename} {params.rule_name} {log} {params.annotation}
        else
            echo "ERROR: {wildcards.filename} marker_id status=$marker_status (predicted=$predicted)" > {log}
            wagtail_log_fail {params.sqlite_db} {wildcards.filename} {params.rule_name} {log} {params.annotation}
        fi
        touch {output.logged}
        wagtail_cleanup {log}
        exit 0
        """

# ── manifest ──────────────────────────────────────────────────────────────────
rule manifest:
    group: "wagtail"
    input:
        filemap = config["file_map"],
        logged  = f"{run_dir}/{run_name}/1_marker_logged/{{filename}}.logged"
    output:
        output_dir = directory(f"{run_dir}/{run_name}/1_manifest/{{filename}}"),
        manifest   = f"{run_dir}/{run_name}/1_manifest/{{filename}}/{{filename}}_manifest.csv"
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

# ── import_and_qc ─────────────────────────────────────────────────────────────
# Merged from qiime2_import + quality_control

rule import_and_qc:
    group: "wagtail"
    input:
        manifest = f"{run_dir}/{run_name}/1_manifest/{{filename}}/{{filename}}_manifest.csv"
    output:
        filtered = f"{run_dir}/{run_name}/2_qc_wagtail/{{filename}}-filtered.qza",
        stats    = f"{run_dir}/{run_name}/2_qc_wagtail/{{filename}}-qc-stats.qza"
    wildcard_constraints:
        filename = r"[^\.]+"
    params:
        sqlite_db = lambda wc: tmp_db(wc.filename),
        rule_name = "import_and_qc",
        upstream  = "manifest",
        tmpdir    = lambda wc: tmp_dir(wc.filename),
        logger    = f"{script_dir}/wagtail_sqlite_logger.py",
        checker   = f"{script_dir}/check_sqlite_status.py",
        annotation= f"{script_dir}/error_annotation.py",
        bash_lib  = BASH_LIB
    conda:
        "envs/qiime2.yml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/import_and_qc.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/import_and_qc.benchmark.log"
    threads: 1
    resources:
        mem_mb  = 2000,
        runtime = "1h"
    shell:
        r"""
        source {params.bash_lib}
        LOGGER={params.logger}
        CHECKER={params.checker}
        export TMPDIR={params.tmpdir}; mkdir -p $TMPDIR
        wagtail_check_upstream {params.sqlite_db} {wildcards.filename} {params.upstream} \
            "{output.filtered} {output.stats}" \
            && touch {output.filtered} {output.stats} \
            && wagtail_log_fail {params.sqlite_db} {wildcards.filename} {params.rule_name} {log} {params.annotation} \
            && wagtail_cleanup {log} && exit 0

        single_qza=$TMPDIR/{wildcards.filename}-single.qza

        set +e
        qiime tools import \
            --type 'SampleData[SequencesWithQuality]' \
            --input-path   {input.manifest} \
            --input-format SingleEndFastqManifestPhred33 \
            --output-path  $single_qza \
            &> {log}
        status=$?

        if [[ $status -eq 0 ]]; then
            qiime quality-filter q-score \
                --i-demux              $single_qza \
                --o-filtered-sequences {output.filtered} \
                --o-filter-stats       {output.stats} \
                &>> {log}
            status=$?
        fi
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
    output:
        representative = f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}-rep-seqs.qza",
        table          = f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}-table.qza",
        stats          = f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}-deblur-stats.qza"
    wildcard_constraints:
        filename = r"[^\.]+"
    params:
        path        = directory(f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}"),
        script      = f"{script_dir}/deblur_all.py",
        resolver    = RESOLVER,
        marker_json = lambda wc: marker_json(wc.filename),
        db_dir      = db_dir,
        sqlite_db = lambda wc: tmp_db(wc.filename),
        rule_name = "deblur",
        upstream  = "import_and_qc",
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
        runtime = lambda wc, attempt: 4 * 90 * attempt   # 6h → 12h → 18h on retry
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
        # Resolve marker + deblur reference at runtime from this sample's JSON
        # (16S → empty reference, deblur_all.py uses QIIME2's built-in for 16S).
        MARKER=$(python {params.resolver} --marker-json {params.marker_json} \
            --db-dir {params.db_dir} --kind deblur --what marker)
        REF=$(python {params.resolver} --marker-json {params.marker_json} \
            --db-dir {params.db_dir} --kind deblur --what path)

        set +e
        python {params.script} \
            -i {input.filtered} -o {params.path} -t {threads} \
            -r {output.representative} -a {output.table} -s {output.stats} \
            --marker    "$MARKER" \
            --reference "$REF" \
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

# ── export_all ────────────────────────────────────────────────────────────────
# Merged from export_seqs + export_and_edit_table

rule export_all:
    group: "wagtail"
    input:
        representative = f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}-rep-seqs.qza",
        table          = f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}-table.qza"
    output:
        seqs_dir  = directory(f"{run_dir}/{run_name}/4_rep_seqs_wagtail/{{filename}}"),
        fasta     = f"{run_dir}/{run_name}/4_rep_seqs_wagtail/{{filename}}/dna-sequences.fasta",
        folder    = directory(f"{run_dir}/{run_name}/4_table_wagtail/{{filename}}_biom"),
        biom      = f"{run_dir}/{run_name}/4_table_wagtail/{{filename}}_biom/{{filename}}.biom",
        table_tsv = f"{run_dir}/{run_name}/4_table_wagtail/{{filename}}_table.tsv",
        edited    = f"{run_dir}/{run_name}/5_taxonomy_wagtail/{{filename}}/{{filename}}_table_edit.tsv"
    wildcard_constraints:
        filename = r"[^\.]+"
    params:
        script_seqs = f"{script_dir}/qiime2_seqs_export.py",
        script_biom = f"{script_dir}/qiime2_biom_export.py",
        name        = "{filename}.biom",
        sqlite_db   = lambda wc: tmp_db(wc.filename),
        rule_name   = "export_all",
        upstream    = "deblur",
        tmpdir      = lambda wc: tmp_dir(wc.filename),
        logger      = f"{script_dir}/wagtail_sqlite_logger.py",
        checker     = f"{script_dir}/check_sqlite_status.py",
        annotation  = f"{script_dir}/error_annotation.py",
        bash_lib    = BASH_LIB
    conda:
        "envs/qiime2.yml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/export_all.log"
    benchmark:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/export_all.benchmark.log"
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
            "{output.fasta} {output.biom} {output.table_tsv} {output.edited}" \
            && mkdir -p {output.seqs_dir} {output.folder} \
            && touch {output.fasta} {output.biom} {output.table_tsv} {output.edited} \
            && wagtail_log_fail {params.sqlite_db} {wildcards.filename} {params.rule_name} {log} {params.annotation} \
            && wagtail_cleanup {log} && exit 0

        export TMPDIR={params.tmpdir}; mkdir -p $TMPDIR
        set +e
        python {params.script_seqs} \
            --input-path  {input.representative} \
            --output-path {output.seqs_dir} \
            &> {log}
        status=$?

        if [[ $status -eq 0 ]]; then
            {{ python {params.script_biom} \
                --input-path   {input.table} \
                --output-path  {output.folder} \
                --new-filename {params.name} \
            && biom convert -i {output.biom} -o {output.table_tsv} --to-tsv \
            && tail -n +3 {output.table_tsv} > {output.edited}; }} &>> {log}
            status=$?
        fi
        set -e

        if [[ $status -ne 0 ]]; then
            mkdir -p {output.seqs_dir} {output.folder}
            wagtail_log_fail {params.sqlite_db} {wildcards.filename} {params.rule_name} {log} {params.annotation}
            touch {output.fasta} {output.biom} {output.table_tsv} {output.edited}
        else
            wagtail_log_ok {params.sqlite_db} {wildcards.filename} {params.rule_name} {log}
        fi
        wagtail_cleanup {log}
        exit 0
        """

# ── mappy ─────────────────────────────────────────────────────────────────────
rule mappy:
    group: "wagtail"
    input:
        F   = f"{run_dir}/{run_name}/4_rep_seqs_wagtail/{{filename}}/dna-sequences.fasta",
    output:
        align = f"{run_dir}/{run_name}/5_taxonomy_wagtail/{{filename}}/{{filename}}_alignment.tsv",
        meta  = f"{run_dir}/{run_name}/5_taxonomy_wagtail/{{filename}}/{{filename}}_metadata.tsv"
    wildcard_constraints:
        filename = r"[^\.]+"
    params:
        script      = f"{script_dir}/mappy_script.py",
        resolver    = RESOLVER,
        marker_json = lambda wc: marker_json(wc.filename),
        db_dir      = db_dir,
        sample    = "{filename}",
        sqlite_db = lambda wc: tmp_db(wc.filename),
        rule_name = "mappy",
        upstream  = "export_all",
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
        runtime = "12h"
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

        # Resolve the alignment database at runtime from this sample's JSON.
        REF=$(python {params.resolver} \
            --marker-json {params.marker_json} --db-dir {params.db_dir} \
            --kind database --what path)

        set +e
        python {params.script} \
            -i {input.F} -r "$REF" \
            -a {output.align} -m {output.meta} -s {params.sample} \
            &> {log}
        status=$?
        set -e

        wagtail_handle_status $status {params.sqlite_db} {wildcards.filename} {params.rule_name} {log} {params.annotation} \
            "{output.align} {output.meta}"
        wagtail_cleanup {log}
        exit 0
        """

# ── extract_taxonomy ──────────────────────────────────────────────────────────
rule extract_taxonomy:
    group: "wagtail"
    input:
        align       = f"{run_dir}/{run_name}/5_taxonomy_wagtail/{{filename}}/{{filename}}_alignment.tsv",
        table       = f"{run_dir}/{run_name}/5_taxonomy_wagtail/{{filename}}/{{filename}}_table_edit.tsv",
    output:
        condensed = f"{run_dir}/{run_name}/6_condensed_wagtail/{{filename}}_condensed.tsv"
    wildcard_constraints:
        filename = r"[^\.]+"
    params:
        script      = f"{script_dir}/extract_taxonomy.py",
        resolver    = RESOLVER,
        marker_json = lambda wc: marker_json(wc.filename),
        db_dir      = db_dir,
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

        # Resolve the taxonomy reference at runtime; pass -r only when one exists.
        REF=$(python {params.resolver} \
            --marker-json {params.marker_json} --db-dir {params.db_dir} \
            --kind taxonomy --what path)

        set +e
        if [[ -n "$REF" ]]; then
            python {params.script} \
                -i {input.align} -t {input.table} -r "$REF" \
                -o {output.condensed} -s {params.sample} \
                &> {log}
        else
            python {params.script} \
                -i {input.align} -t {input.table} \
                -o {output.condensed} -s {params.sample} \
                &> {log}
        fi
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

# ── gather_data ───────────────────────────────────────────────────────────────
# Wide fan-in over every sample's two final per-sample outputs, collapsed into a single `all_samples.done` sentinel.
# `localrule` so it runs inside the Snakemake controller.
rule gather_data:
    localrule: True
    input:
        metadata  = expand(f"{run_dir}/{run_name}/7_metadata/{{filename}}/{{filename}}_metadata.tsv", filename=filenames),
        community = expand(f"{run_dir}/{run_name}/6_condensed_wagtail/{{filename}}_condensed.tsv",    filename=filenames)
    output:
        gathered = f"{run_dir}/{run_name}/0_logs_wagtail/all_samples.done"
    shell:
        "mkdir -p $(dirname {output.gathered}) && touch {output.gathered}"

# ── cleanup ───────────────────────────────────────────────────────────────────
rule cleanup:
    group: "cleanup"
    input:
        # Gate on gather_data's single sentinel — NOT the per-sample file lists
        gathered = f"{run_dir}/{run_name}/0_logs_wagtail/all_samples.done"
    output:
        merged        = sql_log,
        cleanup_done  = f"{run_dir}/{run_name}/0_logs_wagtail/cleanup_done.txt",
        full_metadata = f"{run_dir}/{run_name}/7_metadata/{run_name}_full_metadata.tsv"
    params:
        merge_script   = f"{script_dir}/sql_merge.py",
        target         = lambda w, output: str(Path(output.merged).parent.parent),
        tmpdir         = lambda w, output: str(Path(output.merged).parent.parent / "0_tmp"),
        metadata_dir   = f"{run_dir}/{run_name}/7_metadata",
        archive_script = f"{script_dir}/wagtail_archive.py",
        archive_cli    = _archive_cli(),     # "" → no archiving requested
        run_name       = run_name
    conda:
        "envs/mappy.yaml"
    log:
        merger       = f"{run_dir}/{run_name}/0_logs_wagtail/merger.log",
        cleaning     = f"{run_dir}/{run_name}/0_logs_wagtail/cleaning.log",
        archive      = f"{run_dir}/{run_name}/0_logs_wagtail/archive.log",
        full_metadata= f"{run_dir}/{run_name}/0_logs_wagtail/full_metadata.log"
    threads: 1
    resources:
        mem_mb  = 2000,
        runtime = "2h"
    shell:
        r"""

        # Step 1: Merge per-sample SQLite logs (script globs 0_tmp/*_log.sql)
        python {params.merge_script} {output.merged} {params.tmpdir} &> {log.merger}
        sleep 2

        # Step 2: Merge per-sample metadata into one file.
        #   Header: from the first file via `find -print -quit` (NOT `find | sort |
        #   head`, which SIGPIPEs `sort` under `set -o pipefail` once its output
        #   exceeds the pipe buffer at large sample counts → silent cleanup failure).
        #   Bodies: `find | sort | while read` drains every line (no early pipe
        #   close → SIGPIPE-safe) and keeps no path list on a command line.
        (
            first_meta=$(find {params.metadata_dir} -mindepth 2 -name '*_metadata.tsv' -type f -print -quit)
            if [[ -n "$first_meta" ]]; then
                head -n 1 "$first_meta" > {output.full_metadata}
                find {params.metadata_dir} -mindepth 2 -name '*_metadata.tsv' -type f \
                    | sort | while IFS= read -r f; do tail -n +2 -q "$f"; done \
                    >> {output.full_metadata}
            else
                : > {output.full_metadata}
            fi
        ) &> {log.full_metadata}

        # Step 3: Drop the scratch dir, then (optionally) archive. 0_tmp is
        #         removed first so it never lands in a 0_ archive; the 1-5
        #         intermediates are still present here for --archive-intermediate
        #         / --archive-all.
        rm -rf {params.tmpdir}
        if [[ -n "{params.archive_cli}" ]]; then
            python {params.archive_script} \
                --run-dir  {params.target} \
                --run-name {params.run_name} \
                {params.archive_cli} \
                &> {log.archive}
        else
            echo "No archiving requested" &> {log.archive}
        fi

        # Step 4: Delete the intermediate stage directories (1-5)
        rm -rf {params.target}/[1-5]_*
        touch {output.cleanup_done}
        echo "Cleanup completed" &> {log.cleaning}
        """
