# Wagtail RESUME pipeline (v1.2) — DRAFT
# ──────────────────────────────────────────────────────────────────────────────
# Runs the SECOND HALF of Wagtail on already-produced intermediate results:
#   mappy → extract_taxonomy → (qiime_stats) → metadata_creation → gather_data → cleanup
#
# Triggered by:  wagtail.py --resume ...
#
# You provide the upstream outputs; this pipeline picks up from there. For the
# current entry point (mappy, stage 5) the following must already exist under
# run/<run_name>/ for every sample:
#     4_rep_seqs_wagtail/<s>/dna-sequences.fasta        (mappy input)
#     5_taxonomy_wagtail/<s>/<s>_table_edit.tsv         (extract_taxonomy input)
#     2_qc_wagtail/<s>-qc-stats.qza                      (qiime_stats input, for metadata)
#     3_deblur_wagtail/<s>-deblur-stats.qza              (qiime_stats input, for metadata)
# Plus markers for runtime reference resolution, supplied one of two ways:
#     • auto:   1_marker_id/<s>_marker.json provided by you, or
#     • forced: pass --amplicon 16S|18S|ITS|CO1 and prepare_markers regenerates them.
#
# Differences from the main pipeline (intentional):
#   - No marker_id / manifest / import_and_qc / deblur / export_all.
#   - ENTRY rules (mappy, qiime_stats) do NOT seed or check an upstream SQLite
#     status — there is no in-run upstream. Logging simply STARTS at the entry
#     rule. Downstream rules (extract_taxonomy, metadata_creation) keep their
#     normal upstream check, which passes because the entry rule logged this run.
#   - cleanup uses the same logic (merge SQLite + metadata, archive, delete 1_-5_).
#
# TODO (generalisation): make the entry stage configurable (config["from_stage"])
# and auto-wire each rule's provided inputs from the stage number. For now the
# entry point is hard-set to `mappy` (stage 5).

import os
from pathlib import Path
from datetime import datetime

# ── Directories ──────────────────────────────────────────────────────────────
parent_dir = os.path.dirname(workflow.basedir)
db_dir     = os.path.join(parent_dir, "database")
script_dir = os.path.join(parent_dir, "bin")
run_dir    = os.path.join(parent_dir, "run")

# ── Config ───────────────────────────────────────────────────────────────────
configfile: "config.yaml"

filename_list = config["sample_list"]
filenames     = [line.strip() for line in open(filename_list) if line.strip()]
run_name      = config["run_name"]

# Amplicon: "all" → auto (markers must be provided); else forced (regenerated).
_AMPLICON_RAW   = config.get("amplicon", "all")
FORCED_AMPLICON = _AMPLICON_RAW.upper()
_VALID          = {"ALL", "16S", "18S", "ITS", "CO1"}
if FORCED_AMPLICON not in _VALID:
    raise WorkflowError(f"config amplicon={_AMPLICON_RAW!r} invalid. "
                        f"Choose {', '.join(sorted(_VALID - {'ALL'}))} or omit.")

# Entry stage is hard-set to `mappy` (stage 5) for now; see the header TODO for
# the planned config["from_stage"] generalisation.

# ── Archive options (shared with the main pipeline) ──────────────────────────
def _truthy(v) -> bool:
    return str(v).strip().lower() in {"1", "true", "yes", "on"}

ARCHIVE_RESULT       = _truthy(config.get("archive_result", ""))
ARCHIVE_ALL          = _truthy(config.get("archive_all", ""))
ARCHIVE_INTERMEDIATE = str(config.get("archive_intermediate", "")).strip()

def _archive_cli() -> str:
    parts = []
    if ARCHIVE_ALL:
        parts.append("--archive-all")
    if ARCHIVE_RESULT:
        parts.append("--archive-result")
    if ARCHIVE_INTERMEDIATE:
        parts.append(f"--archive-intermediate {ARCHIVE_INTERMEDIATE}")
    return " ".join(parts)

# ── Timestamp / SQLite / paths ───────────────────────────────────────────────
timestamp = datetime.now().strftime("%Y%m%d")
sql_log   = f"{run_dir}/{run_name}/0_logs_wagtail/{run_name}_{timestamp}_log.sql"

def tmp_db(filename):  return f"{run_dir}/{run_name}/0_tmp/{filename}_log.sql"
def tmp_dir(filename): return f"{run_dir}/{run_name}/0_tmp/{filename}"

BASH_LIB = os.path.join(workflow.basedir, "wagtail_functions.sh")

# Runtime reference resolver (same as main): reads <s>_marker.json + globs db_dir.
RESOLVER    = f"{script_dir}/resolve_reference.py"
marker_done = f"{run_dir}/{run_name}/1_marker_id/marker_id.done"

def marker_json(filename: str) -> str:
    return f"{run_dir}/{run_name}/1_marker_id/{filename}_marker.json"

# ── rule all ─────────────────────────────────────────────────────────────────
rule all:
    input:
        expand(f"{run_dir}/{run_name}/6_condensed_wagtail/{{filename}}_condensed.tsv", filename=filenames),
        f"{run_dir}/{run_name}/7_metadata/{run_name}_full_metadata.tsv",
        sql_log
    localrule: True

# ── prepare_markers ──────────────────────────────────────────────────────────
# Ensure every sample has a marker JSON for runtime reference resolution.
# Forced mode regenerates them cheaply; auto mode assumes you provided them.
rule prepare_markers:
    localrule: True
    output:
        sentinel = marker_done
    params:
        forced_script = f"{script_dir}/write_forced_markers.py",
        marker_dir    = f"{run_dir}/{run_name}/1_marker_id",
        forced        = FORCED_AMPLICON,
        sample_list   = filename_list,
    shell:
        r"""
        mkdir -p {params.marker_dir}
        if [[ "{params.forced}" != "ALL" ]]; then
            python {params.forced_script} \
                --sample-list {params.sample_list} \
                --output-dir  {params.marker_dir} \
                --marker      {params.forced}
        fi
        touch {output.sentinel}
        """

# ── mappy (ENTRY — no upstream check) ────────────────────────────────────────
rule mappy:
    group: "wagtail"
    input:
        F       = f"{run_dir}/{run_name}/4_rep_seqs_wagtail/{{filename}}/dna-sequences.fasta",
        markers = marker_done,
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
        sample      = "{filename}",
        sqlite_db   = lambda wc: tmp_db(wc.filename),
        rule_name   = "mappy",
        logger      = f"{script_dir}/wagtail_sqlite_logger.py",
        annotation  = f"{script_dir}/error_annotation.py",
        bash_lib    = BASH_LIB
    conda:
        "envs/mappy.yaml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/mappy.log"
    threads: 1
    resources:
        mem_mb  = 4000,
        runtime = "12h"
    shell:
        r"""
        source {params.bash_lib}
        LOGGER={params.logger}
        mkdir -p $(dirname {output.align}) $(dirname {log}) $(dirname {params.sqlite_db})

        # Resume entry point: logging STARTS here — no upstream to seed or check.
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

# ── extract_taxonomy (unchanged from main; upstream mappy logged this run) ────
rule extract_taxonomy:
    group: "wagtail"
    input:
        align = f"{run_dir}/{run_name}/5_taxonomy_wagtail/{{filename}}/{{filename}}_alignment.tsv",
        table = f"{run_dir}/{run_name}/5_taxonomy_wagtail/{{filename}}/{{filename}}_table_edit.tsv",
    output:
        condensed = f"{run_dir}/{run_name}/6_condensed_wagtail/{{filename}}_condensed.tsv"
    wildcard_constraints:
        filename = r"[^\.]+"
    params:
        script      = f"{script_dir}/extract_taxonomy.py",
        resolver    = RESOLVER,
        marker_json = lambda wc: marker_json(wc.filename),
        db_dir      = db_dir,
        sample      = "{filename}",
        sqlite_db   = lambda wc: tmp_db(wc.filename),
        rule_name   = "extract_taxonomy",
        upstream    = "mappy",
        logger      = f"{script_dir}/wagtail_sqlite_logger.py",
        checker     = f"{script_dir}/check_sqlite_status.py",
        annotation  = f"{script_dir}/error_annotation.py",
        bash_lib    = BASH_LIB
    conda:
        "envs/mappy.yaml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/extract_taxonomy.log"
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

        REF=$(python {params.resolver} \
            --marker-json {params.marker_json} --db-dir {params.db_dir} \
            --kind taxonomy --what path)

        set +e
        if [[ -n "$REF" ]]; then
            python {params.script} -i {input.align} -t {input.table} -r "$REF" \
                -o {output.condensed} -s {params.sample} &> {log}
        else
            python {params.script} -i {input.align} -t {input.table} \
                -o {output.condensed} -s {params.sample} &> {log}
        fi
        status=$?
        set -e

        wagtail_handle_status $status {params.sqlite_db} {wildcards.filename} {params.rule_name} {log} {params.annotation} \
            "{output.condensed}"
        wagtail_cleanup {log}
        exit 0
        """

# ── qiime_stats (ENTRY — no upstream check; needs provided 2_/3_ stats) ───────
rule qiime_stats:
    group: "wagtail"
    input:
        qc     = f"{run_dir}/{run_name}/2_qc_wagtail/{{filename}}-qc-stats.qza",
        deblur = f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}-deblur-stats.qza"
    output:
        qc_dir       = directory(f"{run_dir}/{run_name}/2_qc_wagtail/{{filename}}"),
        deblur_dir   = directory(f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}"),
        final_qc     = f"{run_dir}/{run_name}/2_qc_wagtail/{{filename}}/stats.csv",
        final_deblur = f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}/stats.csv"
    wildcard_constraints:
        filename = r"[^\.]+"
    params:
        script     = f"{script_dir}/wagtail_metadata_qiimes.py",
        sqlite_db  = lambda wc: tmp_db(wc.filename),
        rule_name  = "qiime_stats",
        tmpdir     = lambda wc: tmp_dir(wc.filename),
        logger     = f"{script_dir}/wagtail_sqlite_logger.py",
        annotation = f"{script_dir}/error_annotation.py",
        bash_lib   = BASH_LIB
    conda:
        "envs/qiime2.yml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/qiime_stats.log"
    threads: 1
    resources:
        mem_mb  = 640,
        runtime = "1h"
    shell:
        r"""
        source {params.bash_lib}
        LOGGER={params.logger}
        mkdir -p {output.qc_dir} {output.deblur_dir} $(dirname {log}) $(dirname {params.sqlite_db})

        # Resume entry point: no upstream check. Guard against empty inputs only.
        if [[ ! -s {input.qc} ]] || [[ ! -s {input.deblur} ]]; then
            touch {output.final_qc} {output.final_deblur}
            echo "Provided qc/deblur stats missing or empty" > {log}
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
            touch {output.final_qc} {output.final_deblur}
            wagtail_log_fail {params.sqlite_db} {wildcards.filename} {params.rule_name} {log} {params.annotation}
        else
            wagtail_log_ok {params.sqlite_db} {wildcards.filename} {params.rule_name} {log}
        fi
        wagtail_cleanup {log}
        exit 0
        """

# ── metadata_creation (unchanged from main) ──────────────────────────────────
rule metadata_creation:
    group: "wagtail"
    input:
        final_qc     = f"{run_dir}/{run_name}/2_qc_wagtail/{{filename}}/stats.csv",
        final_deblur = f"{run_dir}/{run_name}/3_deblur_wagtail/{{filename}}/stats.csv",
        final_mappy  = f"{run_dir}/{run_name}/5_taxonomy_wagtail/{{filename}}/{{filename}}_metadata.tsv"
    output:
        directory = directory(f"{run_dir}/{run_name}/7_metadata/{{filename}}"),
        final     = f"{run_dir}/{run_name}/7_metadata/{{filename}}/{{filename}}_metadata.tsv"
    wildcard_constraints:
        filename = r"[^\.]+"
    params:
        script     = f"{script_dir}/wagtail_metadata_meta_combine.py",
        sqlite_db  = lambda wc: tmp_db(wc.filename),
        run        = "{filename}",
        rule_name  = "metadata_creation",
        upstream   = "qiime_stats",
        logger     = f"{script_dir}/wagtail_sqlite_logger.py",
        checker    = f"{script_dir}/check_sqlite_status.py",
        annotation = f"{script_dir}/error_annotation.py",
        bash_lib   = BASH_LIB
    conda:
        "envs/mappy.yaml"
    log:
        f"{run_dir}/{run_name}/0_logs_wagtail/{{filename}}/metadata.log"
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

# ── gather_data (localrule fan-in, same as main) ─────────────────────────────
rule gather_data:
    localrule: True
    input:
        metadata  = expand(f"{run_dir}/{run_name}/7_metadata/{{filename}}/{{filename}}_metadata.tsv", filename=filenames),
        community = expand(f"{run_dir}/{run_name}/6_condensed_wagtail/{{filename}}_condensed.tsv",    filename=filenames)
    output:
        gathered = f"{run_dir}/{run_name}/0_logs_wagtail/all_samples.done"
    shell:
        "mkdir -p $(dirname {output.gathered}) && touch {output.gathered}"

# ── cleanup (same logic as main) ─────────────────────────────────────────────
rule cleanup:
    group: "cleanup"
    input:
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
        archive_cli    = _archive_cli(),
        run_name       = run_name
    conda:
        "envs/mappy.yaml"
    log:
        merger        = f"{run_dir}/{run_name}/0_logs_wagtail/merger.log",
        cleaning      = f"{run_dir}/{run_name}/0_logs_wagtail/cleaning.log",
        archive       = f"{run_dir}/{run_name}/0_logs_wagtail/archive.log",
        full_metadata = f"{run_dir}/{run_name}/0_logs_wagtail/full_metadata.log"
    threads: 1
    resources:
        mem_mb  = 2000,
        runtime = "2h"
    shell:
        r"""
        # Same logic as the main pipeline's cleanup (ARG_MAX-safe: no per-sample
        # file lists on the command line). The intermediate range deleted starts
        # at the provided-input stage and runs through stage 5.
        python {params.merge_script} {output.merged} {params.tmpdir} &> {log.merger}
        sleep 2

        # Header via `find -print -quit` (NOT `find | sort | head`, which SIGPIPEs
        # `sort` under `set -o pipefail` at large sample counts → silent failure).
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

        rm -rf {params.target}/[1-5]_*
        touch {output.cleanup_done}
        echo "Cleanup completed" &> {log.cleaning}
        """
