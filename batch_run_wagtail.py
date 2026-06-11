#!/usr/bin/env python3
"""
batch_run_wagtail.py — split an accession list into batches and run the Wagtail
pipeline on each, locally or on a cluster (via mqsub).

Completion is decided by the pipeline's own ``cleanup_done.txt`` sentinel
(written by rule cleanup), not by the queue-exit status: mqsub returns 0 on
*submission* and qstat only tells us a job left the queue — neither proves the
pipeline finished. A batch is "done" iff
``<wagtail_root>/run/<run_name>/0_logs_wagtail/cleanup_done.txt`` exists.
"""

import argparse
import concurrent.futures
import glob
import math
import os
import random
import re
import subprocess
import threading
import time
from dataclasses import dataclass
from datetime import datetime

# Batch files are named batch_<num>.txt; numbers are zero-padded to >= 3 digits
# for nice sorting but parsed width-agnostically so >999 batches still work.
_BATCH_NUM_RE = re.compile(r"^batch_(\d+)\.txt$")
# Serialises appends to the shared log file when running batches in parallel.
_LOG_LOCK = threading.Lock()


def batch_str(num: int) -> str:
    """Zero-padded batch label (>= 3 digits)."""
    return str(num).zfill(3)


# ── Config ────────────────────────────────────────────────────────────────────

@dataclass
class Config:
    input_mode: str
    acc_list_file: str
    input_dir: str
    batch_size: int
    output_dir: str
    log_file: str
    file_map: str
    pipeline: str
    max_workers: int
    run_name_prefix: str
    request_cores: int
    request_mem: int
    request_hours: int
    wait_time: int
    execution_mode: str
    amplicon: str
    conda_prefix: str
    group_components: int
    randomize: bool
    random_seed: int
    start_batch: int
    end_batch: int
    batch_make: bool

    # ── derived paths (single source of truth, fixes scattered string-building) ──
    def run_name(self, bstr: str) -> str:
        return f"{self.run_name_prefix}_batch_{bstr}"

    def batch_file(self, bstr: str) -> str:
        return os.path.join(self.output_dir, f"batch_{bstr}.txt")

    def config_file(self, bstr: str) -> str:
        return os.path.join(self.output_dir, f"config_batch_{bstr}.yaml")

    def status_file(self, bstr: str) -> str:
        return os.path.join(self.output_dir, f"batch_{bstr}_run_status.txt")

    def cleanup_done(self, bstr: str) -> str:
        """Path to the pipeline's completion sentinel for this batch's run."""
        wagtail_root = os.path.dirname(os.path.abspath(self.pipeline))
        return os.path.join(wagtail_root, "run", self.run_name(bstr),
                            "0_logs_wagtail", "cleanup_done.txt")


# ── Argument parsing & validation ─────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Batch run Wagtail pipeline.")
    parser.add_argument("-a", "--acc_list_file", help="Path to the accession list file (for CSV input mode).")
    parser.add_argument("-i", "--input_dir", help="Directory containing batch_*.txt files (for directory input mode).")
    parser.add_argument("-b", "--batch_size", type=int, default=5000, help="Number of accessions per batch (only for CSV input mode).")
    parser.add_argument("-o", "--output_dir", required=True, help="Directory to save batch configs and logs.")
    parser.add_argument("-l", "--log_file", default=None, help="Path to the log file. Required unless --batch-make is used.")
    parser.add_argument("-f", "--file_map", required=True, help="Path to the file map.")
    parser.add_argument("-p", "--pipeline", default=None, help="Path to the Wagtail wrapper (wagtail.py). Required unless --batch-make is used.")
    parser.add_argument("-w", "--max_workers", type=int, default=1, help="Maximum number of parallel workers.")
    parser.add_argument("-r", "--run_name_prefix", required=True, help="Prefix for the run name in the batch configuration.")
    parser.add_argument("-c", "--request_cores", type=int, default=4, help="Number of cores to request per batch job.")
    parser.add_argument("-m", "--request_mem", type=int, default=8, help="Memory (GB) to request per batch job.")
    parser.add_argument("-t", "--request_hours", type=int, default=48, help="Number of hours to request per batch job.")
    parser.add_argument("--start-batch", type=int, default=1, help="Batch number to start from (inclusive, default: 1).")
    parser.add_argument("--end-batch", type=int, default=None, help="Batch number to end at (inclusive, default: last batch).")
    parser.add_argument("--wait-time", type=int, default=30, help="Time to wait between job status checks (seconds).")
    parser.add_argument("--execution-mode", choices=["local", "cluster"], default="local", help="Execution mode: 'local' for local execution, 'cluster' for mqsub.")
    parser.add_argument("--randomize", action="store_true", help="Randomize accessions before batching to balance sample sizes (only for CSV input mode).")
    parser.add_argument("--random-seed", type=int, default=None, help="Random seed for reproducible randomization.")
    parser.add_argument("--amplicon", choices=["16S", "18S", "ITS", "CO1"], default=None,
                        help="Force all batches to a specific amplicon type, skipping auto-detection. "
                             "Omit to auto-detect per sample (default).")
    parser.add_argument("--conda-prefix", default=None,
                        help="Path to the shared conda environment prefix directory "
                             "(passed to Snakemake as --conda-prefix).")
    parser.add_argument("--group-components", type=int, default=30,
                        help="Max per-sample chains merged into one Snakemake 'wagtail' group job "
                             "(--group-components wagtail=N, default 30). Keep this modest: each "
                             "group job re-invokes Snakemake with an explicit target per (rule x "
                             "sample), and a too-large group overflows the cluster submit command "
                             "(OSError: Argument list too long). ~30 leaves wide headroom; raise "
                             "only if submissions stay well under ~128 KB.")
    parser.add_argument("--batch-make", action="store_true",
                        help="Only create batch sample files and config files, then exit. "
                             "Does not run the pipeline. Requires -a (CSV mode) and -r. "
                             "Prints ready-to-use Snakemake commands for each batch.")
    return parser


def config_from_args(args, parser) -> Config:
    """Validate args (single pass, no duplicated checks) and build a Config."""
    # Input mode
    if not args.acc_list_file and not args.input_dir:
        parser.error("Either -a/--acc_list_file (CSV mode) or -i/--input_dir (directory mode) must be specified.")
    if args.acc_list_file and args.input_dir:
        parser.error("Cannot specify both -a/--acc_list_file and -i/--input_dir. Choose one input mode.")

    # --batch-make requires CSV mode
    if args.batch_make and not args.acc_list_file:
        parser.error("--batch-make requires -a/--acc_list_file (CSV mode); directory mode already has pre-made batch files.")

    # Args required for a real run but not for --batch-make
    if not args.batch_make:
        if not args.log_file:
            parser.error("-l/--log_file is required unless --batch-make is used.")
        if not args.pipeline:
            parser.error("-p/--pipeline is required unless --batch-make is used.")
        if os.path.isdir(args.log_file):
            parser.error(f"--log_file must be a file path, not a directory: {args.log_file}")
        if not os.path.exists(args.pipeline):
            parser.error(f"Pipeline script not found: {args.pipeline}")

    input_mode = "csv" if args.acc_list_file else "directory"
    if input_mode == "csv" and not os.path.exists(args.acc_list_file):
        parser.error(f"Accession list file not found: {args.acc_list_file}")
    if input_mode == "directory" and not os.path.exists(args.input_dir):
        parser.error(f"Input directory not found: {args.input_dir}")
    if not os.path.exists(args.file_map):
        parser.error(f"File map not found: {args.file_map}")
    if args.group_components < 1:
        parser.error("--group-components must be >= 1.")

    return Config(
        input_mode=input_mode,
        acc_list_file=args.acc_list_file,
        input_dir=args.input_dir,
        batch_size=args.batch_size,
        output_dir=args.output_dir,
        log_file=args.log_file,
        file_map=args.file_map,
        pipeline=args.pipeline,
        max_workers=args.max_workers,
        run_name_prefix=args.run_name_prefix,
        request_cores=args.request_cores,
        request_mem=args.request_mem,
        request_hours=args.request_hours,
        wait_time=args.wait_time,
        execution_mode=args.execution_mode,
        amplicon=args.amplicon,
        conda_prefix=args.conda_prefix,
        group_components=args.group_components,
        randomize=args.randomize,
        random_seed=args.random_seed,
        start_batch=args.start_batch,
        end_batch=args.end_batch,
        batch_make=args.batch_make,
    )


# ── Small IO helpers ──────────────────────────────────────────────────────────

def log_event(cfg: Config, msg: str) -> None:
    """Timestamped, thread-safe append to the shared log file."""
    with _LOG_LOCK:
        with open(cfg.log_file, "a") as logf:
            logf.write(f"[{datetime.now().isoformat()}] {msg}\n")


def write_text_if_changed(path: str, content: str) -> bool:
    """Write *content* to *path*. If the file exists with different content,
    warn (stale batch/config from changed inputs) and overwrite. Returns True
    if the file was (re)written, False if it was already up to date."""
    if os.path.exists(path):
        with open(path) as fh:
            if fh.read() == content:
                return False
        print(f"  WARNING: overwriting {os.path.basename(path)} — content changed since last run.")
    with open(path, "w") as fh:
        fh.write(content)
    return True


def batch_config_yaml(cfg: Config, bstr: str, sample_file: str) -> str:
    """Render the per-batch config.yaml (single source of truth — used by both
    --batch-make and the live run path)."""
    text = (f"run_name: {cfg.run_name(bstr)}\n"
            f"sample_list: {sample_file}\n"
            f"file_map: {cfg.file_map}\n")
    if cfg.amplicon:
        text += f"amplicon: {cfg.amplicon}\n"
    return text


def read_accessions(path: str) -> list:
    with open(path, "r") as fh:
        return [line.strip() for line in fh if line.strip()]


def parse_batch_number(basename: str):
    """Return the integer batch number from 'batch_<n>.txt', or None.
    Width-agnostic, so >999 batches parse correctly."""
    m = _BATCH_NUM_RE.match(basename)
    return int(m.group(1)) if m else None


# ── Cluster job polling ───────────────────────────────────────────────────────

def check_job_status(job_id: str) -> str:
    """'running' while qstat still knows the job, 'finished' once it leaves the
    queue, 'unknown' if qstat itself errored."""
    try:
        result = subprocess.run(["qstat", "-j", str(job_id)],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return "running" if result.returncode == 0 else "finished"
    except Exception as e:
        print(f"Error checking job status: {e}")
        return "unknown"


def wait_for_job(job_id: str, bstr: str, wait_time: int, max_unknown: int = 5) -> None:
    """Block until the job leaves the queue. Transient 'unknown' qstat results
    are retried (bounded) rather than treated as completion, so a hiccup no
    longer makes us declare the batch done prematurely."""
    print(f"Waiting for batch {bstr} job {job_id} to complete...")
    unknown = 0
    while True:
        status = check_job_status(job_id)
        if status == "finished":
            print(f"Batch {bstr} job {job_id} has left the queue.")
            return
        if status == "running":
            unknown = 0
            print(f"Batch {bstr} job {job_id} still running, waiting {wait_time}s...")
            time.sleep(wait_time)
            continue
        # unknown
        unknown += 1
        if unknown >= max_unknown:
            print(f"Batch {bstr} job {job_id} status unknown {unknown}x; stopping wait "
                  f"(completion will be decided by cleanup_done.txt).")
            return
        print(f"Batch {bstr} job {job_id} status unknown, retrying in {wait_time}s...")
        time.sleep(wait_time)


def parse_job_id(stdout: str):
    """Best-effort job ID from mqsub stdout. Completion no longer depends on
    this (cleanup_done.txt is authoritative), so a miss only affects whether we
    actively poll vs. fall back to a timed wait."""
    for line in stdout.splitlines():
        m = re.search(r"\b(\d{4,})\b", line)
        if m:
            return m.group(1)
    return None


# ── make_batches (--batch-make) ───────────────────────────────────────────────

def make_batches(cfg: Config, accs: list) -> list:
    """Create batch sample files + config files without running anything, and
    print a ready-to-paste Snakemake command per batch."""
    if cfg.randomize:
        if cfg.random_seed is not None:
            random.seed(cfg.random_seed)
            print(f"Randomizing with seed {cfg.random_seed}...")
        else:
            print("Randomizing accessions (no seed — not reproducible)...")
        random.shuffle(accs)

    total_acc = len(accs)
    total_batches = math.ceil(total_acc / cfg.batch_size)
    print(f"\nSplitting {total_acc} accessions into {total_batches} batch(es) "
          f"of up to {cfg.batch_size} each.\n")

    created = []   # (batch_num, batch_str, batch_file, config_file, n_samples)
    for i in range(total_batches):
        bnum = i + 1
        bstr = batch_str(bnum)
        batch_accs = accs[i * cfg.batch_size:(i + 1) * cfg.batch_size]

        sample_file = cfg.batch_file(bstr)
        config_file = cfg.config_file(bstr)
        write_text_if_changed(sample_file, "\n".join(batch_accs) + "\n")
        write_text_if_changed(config_file, batch_config_yaml(cfg, bstr, sample_file))

        created.append((bnum, bstr, sample_file, config_file, len(batch_accs)))
        print(f"  Batch {bstr}: {len(batch_accs):>5} samples  →  {os.path.basename(config_file)}")

    smk = cfg.pipeline if cfg.pipeline else "pipeline/wagtail.smk"
    conda_flag = f" --conda-prefix {cfg.conda_prefix}" if cfg.conda_prefix else ""
    print(f"\n{'─'*70}")
    print("Run each batch in a separate tmux window:\n")
    for _, bstr, _, config_file, _ in created:
        print(f"  # Batch {bstr}")
        print(f"  snakemake -s {smk} \\")
        print(f"    --configfile {config_file} \\")
        print(f"    --profile aqua --jobs 50 \\")
        print(f"    --keep-going --rerun-incomplete \\")
        print(f"    --use-conda{conda_flag}")
        print()
    print(f"{'─'*70}")
    print(f"Total: {total_batches} batch file(s) written to {cfg.output_dir}\n")
    return created


# ── run a single batch ────────────────────────────────────────────────────────

def base_command(cfg: Config, config_file: str) -> list:
    """The wagtail/snakemake command shared by both execution modes."""
    cmd = [
        "python", cfg.pipeline,
        "--profile", "aqua",
        "--configfile", config_file,
        "--keep-going",
        "--conda-frontend", "mamba",
        "--rerun-incomplete",
        "--use-conda",
    ]
    if cfg.conda_prefix:
        cmd += ["--conda-prefix", cfg.conda_prefix]
    return cmd


def run_batch(cfg: Config, batch_num: int, batch_data: dict):
    """Create the per-batch config and run it locally or submit it to the cluster.
    Returns (batch_num, status_string)."""
    bstr = batch_str(batch_num)
    config_file = cfg.config_file(bstr)
    status_file = cfg.status_file(bstr)

    # Resolve the sample list for this batch.
    if batch_data["mode"] == "csv":
        sample_file = cfg.batch_file(bstr)
        write_text_if_changed(sample_file, "\n".join(batch_data["accessions"]) + "\n")
        n_samples = len(batch_data["accessions"])
    else:  # directory mode — use the existing batch file as-is
        sample_file = batch_data["file"]
        n_samples = len(read_accessions(sample_file))

    write_text_if_changed(config_file, batch_config_yaml(cfg, bstr, sample_file))
    log_event(cfg, f"Batch {bstr} STARTED: {n_samples} accessions, config: {config_file}")

    base = base_command(cfg, config_file)

    if cfg.execution_mode == "local":
        wagtail_cmd = base + [
            "--jobs", "50",
            "--local-cores", str(cfg.request_cores),
            "--cores", str(cfg.request_cores * 4),
            "--group-components", f"wagtail={cfg.group_components}",
        ]
        with open(status_file, "w") as sf:
            sf.write("STARTED\n")
        try:
            log_event(cfg, f"Batch {bstr} RUNNING LOCALLY: {' '.join(wagtail_cmd)}")
            print(f"Running batch {bstr} locally...")
            result = subprocess.run(wagtail_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            with open(status_file, "a") as sf:
                sf.write(f"STDOUT:\n{result.stdout}\n")
                sf.write(f"STDERR:\n{result.stderr}\n")
                sf.write(f"EXITCODE={result.returncode}\n")
            if result.returncode == 0:
                status = "COMPLETED successfully"
                with open(status_file, "a") as sf:
                    sf.write("COMPLETED\n")
            else:
                status = f"FAILED with return code {result.returncode}"
                with open(status_file, "a") as sf:
                    sf.write("FAILED\n")
        except Exception as e:
            status = f"FAILED with exception: {e}"
            with open(status_file, "a") as sf:
                sf.write(f"EXCEPTION: {e}\nEXITCODE=1\n")

    else:  # cluster (mqsub)
        wagtail_cmd = base + [
            "--jobs", "50",
            "--local-cores", str(cfg.request_cores),
            "--cores", str(cfg.request_cores * 4),
            "--group-components", f"wagtail={cfg.group_components}",
            "--retries", "3",
        ]
        cmd = [
            "mqsub",
            "-t", str(cfg.request_cores),
            "-m", str(cfg.request_mem),
            "--hours", str(cfg.request_hours),
            "--no-email",
            "--",
        ] + wagtail_cmd

        with open(status_file, "w") as sf:
            sf.write("STARTED\n")
        status = "SUBMITTED to cluster"
        try:
            log_event(cfg, f"Batch {bstr} SUBMITTING TO CLUSTER: {' '.join(cmd)}")
            print(f"Submitting batch {bstr} to cluster...")
            result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            with open(status_file, "a") as sf:
                sf.write(f"MQSUB_STDOUT:\n{result.stdout}\n")
                sf.write(f"MQSUB_STDERR:\n{result.stderr}\n")
                sf.write(f"MQSUB_EXITCODE={result.returncode}\n")

            if result.returncode != 0:
                status = f"FAILED submit (return code {result.returncode})"
                with open(status_file, "a") as sf:
                    sf.write(f"SUBMIT_FAILED\nEXITCODE={result.returncode}\n")
            else:
                job_id = parse_job_id(result.stdout)
                log_event(cfg, f"Batch {bstr} SUBMITTED TO CLUSTER"
                               f"{' job_id=' + job_id if job_id else ''}")
                if job_id:
                    print(f"Batch {bstr} submitted (job {job_id}), waiting for completion...")
                    wait_for_job(job_id, bstr, cfg.wait_time)
                else:
                    print(f"Batch {bstr} submitted (no job ID parsed); waiting {cfg.wait_time}s...")
                    log_event(cfg, f"Batch {bstr} WARNING: could not parse job ID from "
                                   f"mqsub stdout: {result.stdout!r}")
                    time.sleep(cfg.wait_time)

                # Completion is decided by the pipeline's cleanup_done.txt, NOT by
                # the queue-exit: mqsub returns 0 on submission and qstat only says
                # the job left the queue. No sentinel → not finished → rerun later
                # (EXITCODE=2 so Step 2 reruns it rather than skipping).
                if os.path.exists(cfg.cleanup_done(bstr)):
                    status = "COMPLETED (cleanup_done present)"
                    with open(status_file, "a") as sf:
                        sf.write("COMPLETED\nEXITCODE=0\n")
                else:
                    status = "INCOMPLETE (no cleanup_done — will rerun)"
                    with open(status_file, "a") as sf:
                        sf.write("INCOMPLETE\nEXITCODE=2\n")
        except Exception as e:
            status = f"FAILED with exception: {e}"
            with open(status_file, "a") as sf:
                sf.write(f"EXCEPTION: {e}\nEXITCODE=1\n")

    log_event(cfg, f"Batch {bstr} {status}")
    return batch_num, status


# ── Build the batch_info map for a run ────────────────────────────────────────

def build_batch_info(cfg: Config):
    """Return (batch_info, total_batches). CSV mode slices the accession list;
    directory mode discovers existing batch_*.txt files."""
    if cfg.input_mode == "csv":
        print(f"Reading accession list: {cfg.acc_list_file}")
        accs = read_accessions(cfg.acc_list_file)
        total_acc = len(accs)
        total_batches = math.ceil(total_acc / cfg.batch_size)
        print(f"Found {total_acc} accessions, will make {total_batches} batches "
              f"(each {cfg.batch_size}, last batch might be less).")

        if cfg.randomize:
            print("Randomizing accessions to balance batch sizes...")
            if cfg.random_seed is not None:
                random.seed(cfg.random_seed)
                print(f"Using random seed: {cfg.random_seed}")
            random.shuffle(accs)
            print("Accessions randomized.")
        else:
            print("Using original order (no randomization).")

        batch_info = {}
        for i in range(total_batches):
            bnum = i + 1
            batch_info[bnum] = {
                "accessions": accs[i * cfg.batch_size:(i + 1) * cfg.batch_size],
                "batch_str": batch_str(bnum),
                "mode": "csv",
            }
        return batch_info, total_batches

    # directory mode
    print(f"Scanning directory for batch files: {cfg.input_dir}")
    batch_files = sorted(glob.glob(os.path.join(cfg.input_dir, "batch_*.txt")))
    if not batch_files:
        raise ValueError(f"No batch_*.txt files found in {cfg.input_dir}")

    batch_info = {}
    for batch_file in batch_files:
        bnum = parse_batch_number(os.path.basename(batch_file))
        if bnum is None:
            print(f"Warning: Could not extract batch number from {os.path.basename(batch_file)}")
            continue
        batch_info[bnum] = {
            "file": batch_file,
            "batch_str": batch_str(bnum),
            "mode": "directory",
        }
    if not batch_info:
        raise ValueError(f"No valid batch files found in {cfg.input_dir}")

    total_batches = max(batch_info.keys())
    print(f"Found {len(batch_info)} batch files, highest batch number: {total_batches}")
    return batch_info, total_batches


def select_batches_to_run(cfg: Config, batch_info: dict, start_batch: int, end_batch: int):
    """Decide which batches still need running. A batch is skipped if its
    cleanup_done.txt exists (authoritative) or its status file records a
    terminal EXITCODE (0 = done, 1 = sample-level failure not worth rerunning)."""
    batches_to_run = []
    for bnum in range(start_batch, end_batch + 1):
        if bnum not in batch_info:
            print(f"Warning: Batch {batch_str(bnum)} not found, skipping.")
            continue
        bstr = batch_str(bnum)

        # Authoritative completion signal from the pipeline itself.
        if os.path.exists(cfg.cleanup_done(bstr)):
            print(f"Batch {bstr} already completed (cleanup_done present), skipping.")
            continue

        status_file = cfg.status_file(bstr)
        skip = False
        if os.path.exists(status_file):
            with open(status_file) as sf:
                exitcodes = [ln for ln in sf if ln.startswith("EXITCODE=")]
            if exitcodes:
                last_exit = exitcodes[-1].strip().split("=")[-1]
                if last_exit == "0":
                    skip = True
                    print(f"Batch {bstr} already completed (EXITCODE=0), skipping.")
                elif last_exit == "1":
                    skip = True
                    print(f"Batch {bstr} has EXITCODE=1, skipping rerun (likely sample issues).")
                else:
                    print(f"Batch {bstr} has EXITCODE={last_exit}, will rerun.")
            else:
                print(f"Batch {bstr} has no EXITCODE, will rerun.")
        else:
            print(f"Batch {bstr} has no status file, will run.")

        if not skip:
            batches_to_run.append((bnum, batch_info[bnum]))
    return batches_to_run


# ── Logging the run header ────────────────────────────────────────────────────

def log_run_header(cfg: Config) -> None:
    print(f"Input mode: {cfg.input_mode.upper()}")
    print(f"Execution mode: {cfg.execution_mode.upper()}")
    print(f"Amplicon filter: {cfg.amplicon} (forced)" if cfg.amplicon
          else "Amplicon: auto-detect per sample")
    print("Running batches locally for testing..." if cfg.execution_mode == "local"
          else "Running batches on cluster using mqsub...")

    lines = [
        f"\n[{datetime.now().isoformat()}] === Starting batch Wagtail processing ===",
        f"Input mode: {cfg.input_mode.upper()}",
        f"Execution mode: {cfg.execution_mode.upper()}",
        f"Amplicon: {cfg.amplicon if cfg.amplicon else 'auto-detect'}",
    ]
    if cfg.conda_prefix:
        lines.append(f"Conda prefix: {cfg.conda_prefix}")
    if cfg.input_mode == "csv":
        lines.append(f"Accession list: {cfg.acc_list_file}")
        lines.append(f"Batch size: {cfg.batch_size}")
        lines.append(f"Randomize: {cfg.randomize}")
        if cfg.randomize and cfg.random_seed is not None:
            lines.append(f"Random seed: {cfg.random_seed}")
    else:
        lines.append(f"Input directory: {cfg.input_dir}")
    lines.append(f"Output directory: {cfg.output_dir}")
    lines.append(f"Pipeline: {cfg.pipeline}")
    lines.append(f"Group components: wagtail={cfg.group_components}")
    if cfg.execution_mode == "cluster":
        lines.append(f"Cluster resources: {cfg.request_cores} cores, {cfg.request_mem}GB RAM, {cfg.request_hours}h")
    with _LOG_LOCK:
        with open(cfg.log_file, "a") as logf:
            logf.write("\n".join(lines) + "\n")


# ── main ──────────────────────────────────────────────────────────────────────

def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    cfg = config_from_args(args, parser)
    os.makedirs(cfg.output_dir, exist_ok=True)

    # --batch-make: create files + print commands, then exit.
    if cfg.batch_make:
        print(f"Reading accession list: {cfg.acc_list_file}")
        make_batches(cfg, read_accessions(cfg.acc_list_file))
        return 0

    log_run_header(cfg)

    batch_info, total_batches = build_batch_info(cfg)

    start_batch = cfg.start_batch
    end_batch = cfg.end_batch if cfg.end_batch is not None else total_batches
    if start_batch < 1 or start_batch > total_batches:
        raise ValueError(f"--start-batch must be between 1 and {total_batches}")
    if end_batch < start_batch or end_batch > total_batches:
        raise ValueError(f"--end-batch must be between {start_batch} and {total_batches}")

    batches_to_run = select_batches_to_run(cfg, batch_info, start_batch, end_batch)
    print(f"\nWill process {len(batches_to_run)} batches out of {len(batch_info)} available batches.")

    results = []
    finished = 0
    if not batches_to_run:
        print("No batches to process.")
        print(f"\nAll batches processing completed.\nCheck log file for details: {cfg.log_file}")
        return 0

    def record(bn, status):
        nonlocal finished
        print(f"Batch {batch_str(bn)}: {status}")
        results.append((bn, status))
        finished += 1
        print(f"Progress: {finished}/{len(batches_to_run)} "
              f"({finished / len(batches_to_run) * 100:.2f}%)")

    if cfg.max_workers > 1:
        with concurrent.futures.ThreadPoolExecutor(max_workers=cfg.max_workers) as executor:
            future_to_batch = {executor.submit(run_batch, cfg, bnum, data): bnum
                               for bnum, data in batches_to_run}
            for future in concurrent.futures.as_completed(future_to_batch):
                bnum = future_to_batch[future]
                try:
                    bn, status = future.result()
                    record(bn, status)
                except Exception as exc:
                    record(bnum, f"EXCEPTION: {exc}")
    else:
        for bnum, data in batches_to_run:
            try:
                bn, status = run_batch(cfg, bnum, data)
                record(bn, status)
            except Exception as exc:
                record(bnum, f"EXCEPTION: {exc}")

    results.sort(key=lambda x: x[0])
    print("\nAll batches finished.")
    print(f"\nAll batches processing completed.\nCheck log file for details: {cfg.log_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
