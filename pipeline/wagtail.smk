# Wagtail pipeline v1.2
# Changes from v1.1:
#   - marker_id is now a single BATCHED checkpoint over the whole sample list:
#     identify_amplicon.py loads each database once and loops every sample, so
#     the multi-GB minimap2 index is built once per batch instead of once per
#     sample. Fixes the per-sample memory multiplication / OOM and the
#     checkpoint-update storm. Writes all 1_marker_id/{sample}.marker files.
#   - new lightweight per-sample marker_log rule reads each .marker and records
#     OK/FAIL into SQLite under rule_name "marker_id" (manifest's upstream check
#     is unchanged). Grouped with "wagtail"; does no alignment.
#   - manifest now depends on 1_marker_logged/{sample}.logged instead of the
#     marker file directly.
# Changes from v0.15 (carried over from v1.1):
#   - qiime2_import + quality_control merged into import_and_qc
#   - export_seqs + export_and_edit_table merged into export_all
#   - _resolve() @lru_cache; get_marker() cached in _marker_cache
#   - deblur runtime dynamic on retry

import os
import glob
import csv
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
# that amplicon type. Pass via wagtail.py --amplicon or --config amplicon=ITS.
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
_marker_cache: dict = {}

def get_marker(wc):
    """Return the predicted (or forced) amplicon marker for a sample.

    Forced-amplicon mode (--amplicon 16S / 18S / ITS / CO1):
      Returns FORCED_AMPLICON immediately.  The marker_id checkpoint still
      runs (writing trivial marker files) so the marker_log/manifest input
      dependencies are satisfied, but no checkpoint DAG edge is created here.

    Auto-detect mode (default):
      Calls checkpoints.marker_id.get() — now a single batched job over the
      whole sample list — so Snakemake defers DAG construction until that job
      has written every marker file.  Falls back to "16S" when a file is empty
      or unreadable (e.g. failed sample) — downstream rules still soft-fail via
      the SQLite check.
    """
    if not AUTO_DETECT:
        return FORCED_AMPLICON

    if wc.filename in _marker_cache:
        return _marker_cache[wc.filename]

    checkpoints.marker_id.get()
    marker_file = f"{run_dir}/{run_name}/1_marker_id/{wc.filename}.marker"
    try:
        with open(marker_file) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            row = next(reader)
        marker = row["predicted_marker"]
    except (StopIteration, KeyError, OSError):
        marker = "16S"
    _marker_cache[wc.filename] = marker
    return marker

@lru_cache(maxsize=None)
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
        return f"{run_dir}/{run_name}/1_marker_id/{wc.filename}.marker"
    return _resolve(marker, db_dir, pattern=f"{marker}_deblur*")

def get_mappy_ref(wc) -> str:
    marker = get_marker(wc)
    return _resolve(marker, db_dir, pattern=f"{marker}_database*")

def get_tax_ref(wc) -> str:
    marker = get_marker(wc)
    return _resolve(marker, db_dir, pattern=f"{marker}_taxonomy*")

# Identification databases are only needed for auto-detection.
# In forced-amplicon mode these are never passed to identify_amplicon.py,
# so we skip the resolution step (avoids FileNotFoundError for missing DBs).
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
# loads each of the 4 databases once and loops over every sample, so the multi-GB
# minimap2 index is built once per batch instead of once per sample — this is the
# fix for the per-sample memory multiplication / OOM. Writes one .marker file per
# sample into 1_marker_id/. SQLite logging is deferred to the per-sample
# marker_log rule below (keeps this job free of per-sample bookkeeping).
checkpoint marker_id:
    input:
        filemap     = config["file_map"],
        sample_list = filename_list,
    output:
        markers = expand(f"{run_dir}/{run_name}/1_marker_id/{{fn}}.marker", fn=filenames)
    params:
        script          = f"{script_dir}/identify_amplicon.py",
        marker_dir      = f"{run_dir}/{run_name}/1_marker_id",
        db_16S          = DB_16S,
        db_18S          = DB_18S,
        db_ITS          = DB_ITS,
        db_CO1          = DB_CO1,
        forced_amplicon = FORCED_AMPLICON,    # "ALL" → auto-detect; else forced
        n_subsample     = 1000,
        min_confidence  = 0.6,
        minimizer_w     = 19,
        # 10-column header matching identify_amplicon.py's write_marker_file()
        header          = "status\tpredicted_marker\tbest_combined_score\tbest_mean_identity\tbest_hit_fraction\tsecond_marker\tsecond_combined_score\tsecond_mean_identity\tsecond_hit_fraction\tscore_margin",
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
        mkdir -p {params.marker_dir} $(dirname {log})

        # ── Forced-amplicon fast path ─────────────────────────────────────────
        # No alignment needed — write a trivial OK marker for every sample.
        if [[ "{params.forced_amplicon}" != "ALL" ]]; then
            while IFS= read -r s || [[ -n "$s" ]]; do
                [[ -z "$s" ]] && continue
                printf '%s\n' "{params.header}" > {params.marker_dir}/$s.marker
                printf 'OK\t{params.forced_amplicon}\t1.0\t1.0\t1.0\tNone\tNA\tNA\tNA\t1.0\n' \
                    >> {params.marker_dir}/$s.marker
            done < {input.sample_list}
            echo "Forced amplicon {params.forced_amplicon} for all samples" > {log}
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
        """

# ── marker_log ────────────────────────────────────────────────────────────────
# Lightweight per-sample step: read this sample's .marker file (produced by the
# batched checkpoint) and record OK/FAIL into its SQLite log under rule_name
# "marker_id", so manifest's existing upstream check is unchanged. No alignment,
# no DB load — just a file read + one SQLite write. Grouped with "wagtail" so it
# rides along in the same cluster job as manifest/import_and_qc for each sample.
rule marker_log:
    group: "wagtail"
    input:
        marker = f"{run_dir}/{run_name}/1_marker_id/{{filename}}.marker"
    output:
        logged = f"{run_dir}/{run_name}/1_marker_logged/{{filename}}.logged"
    wildcard_constraints:
        filename = r"[^\.]+"
    params:
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

        # Empty/missing marker = failed sample (no reads, missing file) → soft-fail
        if [[ ! -s {input.marker} ]]; then
            echo "ERROR: empty/missing marker file for {wildcards.filename}" > {log}
            wagtail_log_fail {params.sqlite_db} {wildcards.filename} {params.rule_name} {log} {params.annotation}
            touch {output.logged}
            wagtail_cleanup {log}
            exit 0
        fi

        # Column map (10-col): f1 status f2 predicted f4 best_mean_id
        #   f6 second_marker f8 second_mean_id f10 score_margin
        marker_status=$(tail -n 1 {input.marker} | cut -f1)
        predicted=$(tail -n 1 {input.marker} | cut -f2)
        best_id=$(tail -n 1 {input.marker} | cut -f4)
        second=$(tail -n 1 {input.marker} | cut -f6)
        second_id=$(tail -n 1 {input.marker} | cut -f8)
        margin=$(tail -n 1 {input.marker} | cut -f10)

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
# Merged from qiime2_import + quality_control: the intermediate -single.qza
# is kept in TMPDIR and never materialised as a tracked output, eliminating
# one DAG node and one round-trip to shared storage per sample.
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

# ── export_all ────────────────────────────────────────────────────────────────
# Merged from export_seqs + export_and_edit_table: both rules depend only on
# deblur outputs and share the same conda env, so running them sequentially
# in one job eliminates one DAG node and one scheduler round-trip per sample.
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

# ── extract_taxonomy ──────────────────────────────────────────────────────────
rule extract_taxonomy:
    group: "wagtail"
    input:
        align       = f"{run_dir}/{run_name}/5_taxonomy_wagtail/{{filename}}/{{filename}}_alignment.tsv",
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
            -i {input.align} -t {input.table} -r {input.reference} \
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
        rm -rf {params.target}/[1-5]_*
        touch {output.cleanup_done}
        echo "Cleanup completed" &> {log.cleaning}

        # Step 3: Merge per-sample metadata into one file
        (head -n 1 {input.metadata[0]} > {output.full_metadata} \
            && tail -n +2 -q {input.metadata} >> {output.full_metadata}) \
            &> {log.full_metadata}
        """
