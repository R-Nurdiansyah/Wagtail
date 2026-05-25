#!/usr/bin/env bash
# wagtail_functions.sh
# Shared bash helpers for the Wagtail Snakemake pipeline.
# Source this at the top of every rule shell block:
#   source {params.bash_lib}
#
# Functions
# ─────────────────────────────────────────────────────────────────────────────
# wagtail_check_upstream  DB SAMPLE RULE TOUCH_FILES
#   Returns exit code 1 (upstream FAILED) or 0 (OK).
#   Caller uses && chaining to handle the FAILED case inline.
#
# wagtail_log_ok          DB SAMPLE RULE LOG
#   Log a successful rule execution.
#
# wagtail_log_fail        DB SAMPLE RULE LOG ANNOTATION_SCRIPT
#   Annotate and log a failed rule execution.
#
# wagtail_handle_status   STATUS DB SAMPLE RULE LOG ANNOTATION_SCRIPT TOUCH_FILES
#   Convenience wrapper: calls log_ok or log_fail based on $status,
#   and touches output files on failure.
#
# wagtail_cleanup         LOG
#   Sleep 2s then delete non-log/non-sql files in the log directory.
# ─────────────────────────────────────────────────────────────────────────────

wagtail_check_upstream() {
    local db="$1" sample="$2" rule="$3"
    local status
    status=$(python "${CHECKER}" "$db" "$sample" "$rule" 2>/dev/null || echo "FAILED")
    if [[ "$status" == "FAILED" ]]; then
        return 0   # upstream FAILED → return 0 so && chain executes for soft-fail handling
    fi
    return 1  # upstream OK → return 1 so && chain short-circuits; set -e safe inside && list
}

wagtail_log_ok() {
    local db="$1" sample="$2" rule="$3" log="$4"
    python "${LOGGER}" "$db" "$sample" "$rule" OK "$log" ""
}

wagtail_log_fail() {
    local db="$1" sample="$2" rule="$3" log="$4" annotation_script="$5"
    local error
    error=$(python "$annotation_script" "$log" "$rule" 2>/dev/null || echo "")
    python "${LOGGER}" "$db" "$sample" "$rule" FAILED "$log" "$error"
    echo "Upstream fails" >> "$log"
}

wagtail_handle_status() {
    local status="$1" db="$2" sample="$3" rule="$4" log="$5"
    local annotation_script="$6"
    shift 6
    local touch_files=("$@")   # remaining args are files to touch on failure

    if [[ $status -ne 0 ]]; then
        for f in "${touch_files[@]}"; do touch $f; done
        wagtail_log_fail "$db" "$sample" "$rule" "$log" "$annotation_script"
    else
        wagtail_log_ok "$db" "$sample" "$rule" "$log"
    fi
}

wagtail_cleanup() {
    local log="$1"
    sleep 2
    find "$(dirname "$log")" -type f \
        ! -name "$(basename "$log")" \
        ! -name "*.log" \
        ! -name "*.sql" \
        -delete
}
