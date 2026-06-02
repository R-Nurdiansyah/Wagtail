#!/usr/bin/env python3
import os
import math
import subprocess
import concurrent.futures
import argparse
import glob
import random
import re
import time
from datetime import datetime

# === Argparse for Config ===
parser = argparse.ArgumentParser(description="Batch run Wagtail pipeline.")
parser.add_argument("-a", "--acc_list_file", help="Path to the accession list file (for CSV input mode).")
parser.add_argument("-i", "--input_dir", help="Directory containing batch_*.txt files (for directory input mode).")
parser.add_argument("-b", "--batch_size", type=int, default=5000, help="Number of accessions per batch (only for CSV input mode).")
parser.add_argument("-o", "--output_dir", required=True, help="Directory to save batch configs and logs.")
parser.add_argument("-l", "--log_file", required=True, help="Path to the log file.")
parser.add_argument("-f", "--file_map", required=True, help="Path to the file map.")
parser.add_argument("-p", "--pipeline", required=True, help="Path to the Snakemake pipeline.")
parser.add_argument("-w", "--max_workers", type=int, default=1, help="Maximum number of parallel workers.")
parser.add_argument("-r", "--run_name_prefix", required=True, help="Prefix for the run name in the batch configuration.")
parser.add_argument("-c", "--request_cores", type=int, default=4, help="Number of cores to request per batch job.")
parser.add_argument("-m", "--request_mem", type=int, default=8, help="Memory (GB) to request per batch job.")
parser.add_argument("-t", "--request_hours", type=int, default=48, help="Number of hours to request per batch job.")
parser.add_argument("--start-batch", type=int, default=1, help="Batch number to start from (inclusive, default: 1).")
parser.add_argument("--end-batch", type=int, default=None, help="Batch number to end at (inclusive, default: last batch).")
parser.add_argument("--wait-time", type=int, default=30, help="Time to wait between job status checks (seconds).")
parser.add_argument("--execution-mode", choices=["local", "cluster"], default="local", help="Execution mode: 'local' for local execution, 'cluster' for mqsub (default: cluster).")
parser.add_argument("--randomize", action="store_true", help="Randomize accessions before batching to balance sample sizes (only for CSV input mode).")
parser.add_argument("--random-seed", type=int, default=None, help="Random seed for reproducible randomization.")
parser.add_argument("--amplicon", choices=["16S", "18S", "ITS", "CO1"], default=None,
                    help="Force all batches to a specific amplicon type, skipping auto-detection. "
                         "Omit to auto-detect per sample (default).")
parser.add_argument("--conda-prefix", default=None,
                    help="Path to the shared conda environment prefix directory "
                         "(passed to Snakemake as --conda-prefix).")
args = parser.parse_args()

# === Config ===
# Validate input modes
if not args.acc_list_file and not args.input_dir:
    raise ValueError("Either --acc_list_file (CSV mode) or --input_dir (directory mode) must be specified.")
if args.acc_list_file and args.input_dir:
    raise ValueError("Cannot specify both --acc_list_file and --input_dir. Choose one input mode.")

input_mode = "csv" if args.acc_list_file else "directory"
acc_list_file = args.acc_list_file
input_dir = args.input_dir
batch_size = args.batch_size
output_dir = args.output_dir
os.makedirs(output_dir, exist_ok=True)
log_file = args.log_file
constant_filemap = args.file_map
pipeline = args.pipeline
max_workers = args.max_workers
run_name_prefix = args.run_name_prefix
request_cores = args.request_cores
request_mem = args.request_mem
request_hours = args.request_hours
wait_time = args.wait_time
execution_mode = args.execution_mode
amplicon = args.amplicon           # None → auto-detect; "16S"/"18S"/"ITS"/"CO1" → forced
conda_prefix = args.conda_prefix   # None → Snakemake default

# === Early validation =========================================================
if os.path.isdir(log_file):
    raise ValueError(f"--log_file must be a file path, not a directory: {log_file}")
if input_mode == "csv" and not os.path.exists(acc_list_file):
    raise ValueError(f"Accession list file not found: {acc_list_file}")
if input_mode == "directory" and not os.path.exists(input_dir):
    raise ValueError(f"Input directory not found: {input_dir}")
if not os.path.exists(constant_filemap):
    raise ValueError(f"File map not found: {constant_filemap}")
if not os.path.exists(pipeline):
    raise ValueError(f"Pipeline script not found: {pipeline}")

print(f"Input mode: {input_mode.upper()}")
print(f"Execution mode: {execution_mode.upper()}")
if amplicon:
    print(f"Amplicon filter: {amplicon} (forced)")
else:
    print("Amplicon: auto-detect per sample")
if execution_mode == "local":
    print("Running batches locally for testing...")
else:
    print("Running batches on cluster using mqsub...")

# Write log header now that validation has passed
with open(log_file, "a") as logf:
    logf.write(f"\n[{datetime.now().isoformat()}] === Starting batch Wagtail processing ===\n")
    logf.write(f"Input mode: {input_mode.upper()}\n")
    logf.write(f"Execution mode: {execution_mode.upper()}\n")
    logf.write(f"Amplicon: {amplicon if amplicon else 'auto-detect'}\n")
    if conda_prefix:
        logf.write(f"Conda prefix: {conda_prefix}\n")
    if input_mode == "csv":
        logf.write(f"Accession list: {acc_list_file}\n")
        logf.write(f"Batch size: {batch_size}\n")
        logf.write(f"Randomize: {args.randomize}\n")
        if args.randomize and args.random_seed:
            logf.write(f"Random seed: {args.random_seed}\n")
    else:
        logf.write(f"Input directory: {input_dir}\n")
    logf.write(f"Output directory: {output_dir}\n")
    logf.write(f"Pipeline: {pipeline}\n")
    if execution_mode == "cluster":
        logf.write(f"Cluster resources: {request_cores} cores, {request_mem}GB RAM, {request_hours}h\n")

# === Step 1: Process input based on mode ===
if input_mode == "csv":
    # CSV input mode - read accessions and create batches
    print(f"Reading accession list: {acc_list_file}")
    with open(acc_list_file, "r") as f:
        accs = [line.strip() for line in f if line.strip()]
    
    total_acc = len(accs)
    total_batches = math.ceil(total_acc / batch_size)
    print(f"Found {total_acc} accessions, will make {total_batches} batches (each {batch_size}, last batch might be less).")
    
    # Optional randomization
    if args.randomize:
        print("Randomizing accessions to balance batch sizes...")
        if args.random_seed is not None:
            random.seed(args.random_seed)
            print(f"Using random seed: {args.random_seed}")
        random.shuffle(accs)
        print("Accessions randomized.")
    else:
        print("Using original order (no randomization).")
    
    # Create batch info for CSV mode
    batch_info = {}
    for i in range(total_batches):
        batch_num = i + 1
        start_idx = i * batch_size
        end_idx = min((i + 1) * batch_size, total_acc)
        batch_accs = accs[start_idx:end_idx]
        batch_info[batch_num] = {
            'accessions': batch_accs,
            'batch_str': str(batch_num).zfill(3),
            'mode': 'csv'
        }

else:
    # Directory input mode - find existing batch files
    print(f"Scanning directory for batch files: {input_dir}")
    batch_files = glob.glob(os.path.join(input_dir, "batch_*.txt"))
    batch_files.sort()
    
    if not batch_files:
        raise ValueError(f"No batch_*.txt files found in {input_dir}")
    
    # Create batch info for directory mode
    batch_info = {}
    for batch_file in batch_files:
        basename = os.path.basename(batch_file)
        if basename.startswith("batch_") and basename.endswith(".txt"):
            try:
                batch_str = basename[6:9]  # Extract the 3-digit number
                batch_num = int(batch_str)
                batch_info[batch_num] = {
                    'file': batch_file,
                    'batch_str': batch_str,
                    'mode': 'directory'
                }
            except ValueError:
                print(f"Warning: Could not extract batch number from {basename}")
                continue
    
    if not batch_info:
        raise ValueError(f"No valid batch files found in {input_dir}")
    
    total_batches = max(batch_info.keys())
    print(f"Found {len(batch_info)} batch files, highest batch number: {total_batches}")

# === Batch range selection ===
start_batch = args.start_batch
end_batch = args.end_batch if args.end_batch is not None else total_batches

if start_batch < 1 or start_batch > total_batches:
    raise ValueError(f"--start-batch must be between 1 and {total_batches}")
if end_batch < start_batch or end_batch > total_batches:
    raise ValueError(f"--end-batch must be between {start_batch} and {total_batches}")

# === Helper functions ===
def check_job_status(job_id):
    """Check the status of a job using qstat."""
    try:
        result = subprocess.run(["qstat", "-j", str(job_id)], 
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.returncode == 0:
            return "running"
        else:
            return "finished"
    except Exception as e:
        print(f"Error checking job status: {e}")
        return "unknown"

def wait_for_job(job_id, batch_str, wait_time):
    """Wait for a job to complete."""
    print(f"Waiting for batch {batch_str} job {job_id} to complete...")
    while True:
        status = check_job_status(job_id)
        if status == "finished":
            print(f"Batch {batch_str} job {job_id} has completed.")
            break
        elif status == "running":
            print(f"Batch {batch_str} job {job_id} is still running, waiting {wait_time} seconds...")
            time.sleep(wait_time)
        else:
            print(f"Batch {batch_str} job {job_id} status unknown, assuming finished.")
            break

# === Function to process 1 batch ===
def run_batch(batch_num, batch_data):
    """Process a single batch by creating config files and running either locally or on cluster."""
    batch_str = str(batch_num).zfill(3)
    batch_acc_filename = os.path.join(output_dir, f"batch_{batch_str}.txt")
    config_filename = os.path.join(output_dir, f"config_batch_{batch_str}.yaml")
    run_status_filename = os.path.join(output_dir, f"batch_{batch_str}_run_status.txt")

    # Handle different input modes
    if batch_data['mode'] == 'csv':
        batch_accs = batch_data['accessions']
        # Create batch file if it doesn't exist
        if not os.path.exists(batch_acc_filename):
            with open(batch_acc_filename, "w") as f:
                for acc in batch_accs:
                    f.write(acc + "\n")
    else:  # directory mode
        batch_acc_filename = batch_data['file']  # Use existing file
        with open(batch_acc_filename, "r") as f:
            batch_accs = [line.strip() for line in f if line.strip()]

    # Create config file if it doesn't exist
    if not os.path.exists(config_filename):
        with open(config_filename, "w") as cf:
            cf.write(f"run_name: {run_name_prefix}_batch_{batch_str}\n")
            cf.write(f"sample_list: {batch_acc_filename}\n")
            cf.write(f"file_map: {constant_filemap}\n")
            if amplicon:
                cf.write(f"amplicon: {amplicon}\n")
    
    # Log that the batch is starting
    with open(log_file, "a") as logf:
        logf.write(f"[{datetime.now().isoformat()}] Batch {batch_str} STARTED: {len(batch_accs)} accessions, config: {config_filename}\n")

    # Build the base wagtail command shared by both execution modes.
    # amplicon is written to the config file above, so Snakemake picks it up
    # via config.get("amplicon", "all") without needing an extra CLI flag.
    _base_cmd = [
        "python", pipeline,
        "--profile", "aqua",
        "--configfile", config_filename,
        "--keep-going",
        "--conda-frontend", "mamba",
        "--rerun-incomplete",
        "--use-conda",
    ]
    if conda_prefix:
        _base_cmd += ["--conda-prefix", conda_prefix]

    if execution_mode == "local":
        # Local execution mode
        wagtail_cmd = _base_cmd + [
            "--jobs", "50",
            "--local-cores", str(request_cores),
            "--cores", str(request_cores * 4),
            "--group-components", "wagtail=160",
        ]
        
        # Write a checkpoint file before running
        with open(run_status_filename, "w") as status_file:
            status_file.write("STARTED\n")
        
        try:
            # Log command execution
            with open(log_file, "a") as logf:
                logf.write(f"[{datetime.now().isoformat()}] Batch {batch_str} RUNNING LOCALLY: {' '.join(wagtail_cmd)}\n")
            
            print(f"Running batch {batch_str} locally...")
            
            # Run the job locally
            result = subprocess.run(wagtail_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            
            # Log the output
            with open(run_status_filename, "a") as status_file:
                status_file.write(f"STDOUT:\n{result.stdout}\n")
                status_file.write(f"STDERR:\n{result.stderr}\n")
                status_file.write(f"EXITCODE={result.returncode}\n")
            
            if result.returncode == 0:
                status = "COMPLETED successfully"
                with open(run_status_filename, "a") as status_file:
                    status_file.write(f"COMPLETED\n")
            else:
                status = f"FAILED with return code {result.returncode}"
                with open(run_status_filename, "a") as status_file:
                    status_file.write(f"FAILED\n")
        
        except Exception as e:
            status = f"FAILED with exception: {str(e)}"
            with open(run_status_filename, "a") as status_file:
                status_file.write(f"EXCEPTION: {str(e)}\nEXITCODE=1\n")

    else:
        # Cluster execution mode (mqsub)
        wagtail_cmd = _base_cmd + [
            "--jobs", "50",
            "--local-cores", str(request_cores),
            "--cores", str(request_cores * 4),
            "--group-components", "wagtail=50",
            "--retries", "3",
        ]

        cmd = [
            "mqsub",
            "-t", str(request_cores),
            "-m", str(request_mem),
            "--hours", str(request_hours),
            "--no-email",
            "--"
        ] + wagtail_cmd

        # Write a checkpoint file before running
        with open(run_status_filename, "w") as status_file:
            status_file.write("STARTED\n")
        status = "SUBMITTED to cluster"
        try:
            # Log command execution
            with open(log_file, "a") as logf:
                logf.write(f"[{datetime.now().isoformat()}] Batch {batch_str} SUBMITTING TO CLUSTER: {' '.join(cmd)}\n")
            
            print(f"Submitting batch {batch_str} to cluster...")
            
            # Submit the job
            result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            
            # Log the mqsub output for debugging
            with open(run_status_filename, "a") as status_file:
                status_file.write(f"MQSUB_STDOUT:\n{result.stdout}\n")
                status_file.write(f"MQSUB_STDERR:\n{result.stderr}\n")
                status_file.write(f"MQSUB_EXITCODE={result.returncode}\n")
            
            if result.returncode == 0:
                # Job submitted successfully — parse the job ID from mqsub stdout
                # so we can poll qstat instead of guessing at output files.
                job_id = None
                for line in result.stdout.splitlines():
                    m = re.search(r'\b(\d{4,})\b', line)  # job IDs are typically 4+ digits
                    if m:
                        job_id = m.group(1)
                        break

                with open(log_file, "a") as logf:
                    logf.write(f"[{datetime.now().isoformat()}] Batch {batch_str} SUBMITTED TO CLUSTER"
                               f"{' job_id=' + job_id if job_id else ''}\n")

                if job_id:
                    print(f"Batch {batch_str} submitted (job {job_id}), waiting for completion...")
                    wait_for_job(job_id, batch_str, wait_time)
                    status = "COMPLETED (cluster job left queue)"
                    with open(run_status_filename, "a") as status_file:
                        status_file.write(f"COMPLETED\nEXITCODE=0\n")
                else:
                    # mqsub output did not contain a parseable job ID —
                    # fall back to a simple timed wait, then move on.
                    print(f"Batch {batch_str} submitted (no job ID parsed from mqsub output); "
                          f"waiting {wait_time}s then continuing...")
                    with open(log_file, "a") as logf:
                        logf.write(f"[{datetime.now().isoformat()}] Batch {batch_str} WARNING: "
                                   f"could not parse job ID from mqsub stdout: {result.stdout!r}\n")
                    time.sleep(wait_time)
                    status = "SUBMITTED (job ID unknown — not waited)"
                    with open(run_status_filename, "a") as status_file:
                        status_file.write(f"SUBMITTED\nEXITCODE=0\n")
            else:
                status = f"FAILED pipeline (return code {result.returncode})"
                with open(run_status_filename, "a") as status_file:
                    status_file.write(f"SUBMIT_FAILED\nEXITCODE={result.returncode}\n")
        
        except Exception as e:
            status = f"FAILED with exception: {str(e)}"
            with open(run_status_filename, "a") as status_file:
                status_file.write(f"EXCEPTION: {str(e)}\nEXITCODE=1\n")

    # Final log entry
    with open(log_file, "a") as logf:
        logf.write(f"[{datetime.now().isoformat()}] Batch {batch_str} {status}\n")
    
    return batch_num, status

# === Step 2: Check which batches need to be run ===
batches_to_run = []

for batch_num in range(start_batch, end_batch + 1):
    if batch_num not in batch_info:
        print(f"Warning: Batch {batch_num:03d} not found, skipping.")
        continue
        
    batch_str = str(batch_num).zfill(3)
    run_status_filename = os.path.join(output_dir, f"batch_{batch_str}_run_status.txt")
    
    skip = False
    rerun = False

    if os.path.exists(run_status_filename):
        with open(run_status_filename) as sf:
            lines = sf.readlines()
            exitcodes = [line for line in lines if line.startswith("EXITCODE=")]
            if exitcodes:
                last_exit = exitcodes[-1].strip().split("=")[-1]
                if last_exit == "0":
                    skip = True
                    print(f"Batch {batch_str} already completed (EXITCODE=0), skipping.")
                elif last_exit == "1":
                    skip = True
                    print(f"Batch {batch_str} has EXITCODE=1, skipping rerun (could be sample issues).")
                else:
                    rerun = True
                    print(f"Batch {batch_str} has EXITCODE={last_exit}, will rerun.")
            else:
                rerun = True
                print(f"Batch {batch_str} has no EXITCODE, will rerun.")
    else:
        rerun = True
        print(f"Batch {batch_str} has no status file, will run.")

    if not skip:
        batches_to_run.append((batch_num, batch_info[batch_num]))

print(f"\nWill process {len(batches_to_run)} batches out of {len(batch_info)} available batches.")

# === Step 3: Run batches (sequential to avoid overwhelming the system) ===
results = []
finished_batches = 0

if batches_to_run:
    if max_workers > 1:
        # Parallel execution
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_batch = {executor.submit(run_batch, batch_num, batch_data): batch_num for batch_num, batch_data in batches_to_run}
            for future in concurrent.futures.as_completed(future_to_batch):
                batch_num = future_to_batch[future]
                try:
                    bn, status = future.result()
                    print(f"Batch {bn:03d}: {status}")
                    results.append((bn, status))
                except Exception as exc:
                    print(f"Batch {batch_num:03d} generated an exception: {exc}")
                    results.append((batch_num, f"EXCEPTION: {exc}"))
                finished_batches += 1
                percentage = (finished_batches / len(batches_to_run)) * 100
                print(f"Progress: {finished_batches}/{len(batches_to_run)} ({percentage:.2f}%)")
    else:
        # Sequential execution (recommended for job queue systems)
        for batch_num, batch_data in batches_to_run:
            try:
                bn, status = run_batch(batch_num, batch_data)
                print(f"Batch {bn:03d}: {status}")
                results.append((bn, status))
            except Exception as exc:
                print(f"Batch {batch_num:03d} generated an exception: {exc}")
                results.append((batch_num, f"EXCEPTION: {exc}"))
            finished_batches += 1
            percentage = (finished_batches / len(batches_to_run)) * 100
            print(f"Progress: {finished_batches}/{len(batches_to_run)} ({percentage:.2f}%)")

    results.sort(key=lambda x: x[0])
    print("\nAll batches finished.")
else:
    print("No batches to process.")

print("\nAll batches processing completed.")

print("\nAll batches processing completed.")
print(f"Check log file for details: {log_file}")