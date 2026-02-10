![Wagtail Logo](https://github.com/R-Nurdiansyah/Wagtail/blob/development/wagtail_logo_(29-7-2024).png?raw=true)

# Wagtail
_Version 0.14_

Wagtail is an accurate and scalable tools to analyze 16S dataset with easy-to-swap reference database. As for now, the functionality is **limited to 16S amplicon sequences sequenced from Illumina sequencer only**.\
The tool is a combination of Qiime 2 [Deblur](https://github.com/biocore/deblur/blob/master/README.md) plugin for quality control and Minimap2 aligner to align representative sequences from deblur. The result is a taxonomy profile in the format of condense and running metadata. The combination is weaved using [Snakemake](https://snakemake.github.io/) to allow easy reproducibility, benchmarking, and packaging.\
[Minimap2](https://github.com/lh3/minimap2), in the form of its python interface (Mappy), is chosen to its versatility in using reference database without/with minor modification to the reference database. User only need to give the database in the fasta form or use the Minimap2 to help create the index for faster reading (optional).

KEY FEATURES
- Accurate taxonomic profiling using custom databases
- Scalable batch processing for large datasets
- Quality control and denoising with QIIME2 Deblur
- Fast alignment with Minimap2/Mappy
- Comprehensive logging and metadata generation

## How to use
Wagtail is built on Snakemake so executing the tool should be based on how Snakemake is executed. Make sure Snakemake is installed in the system before using it. 
Please check [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/index.html) for further reading in Snakemake execution.

### Prerequisites
- Snakemake installed
- Conda/Mamba for environment management
- Reference database in FASTA format (e.g., danica.mmi in database/ folder)

### Input
Wagtail needs 4 kind of inputs to work:
1. Database. Before executing the pipeline, please put your desired database in the database in the fasta format or in indexed format.
2. Input data in fastq/fastq.gz, recommended to put inside the data directory. Recommended format: filename/accession-1.fastq.gz **Wagtail only need the first-end of the read to work and does not support paired-end. Should your data is paired-end, then use the forward-end only**
3. file_map, contain 2 column with filename/accession in the first column and the absolute path to the file in the second column
4. sample_list, contain only filename/accession.
Make sure that the no 3 and no 4 have the same filename/accession. We provide the filemap creation script in the bin directory and example input data, filemap, and sample list in data directory for your convenience.

### Directories
In the Wagtail directory, you will have several directories:
1. bin (contain python scripts important to the Wagtail functionality)
2. pipeline (contain Wagtail snakemake (Wagtail.smk), config.yaml, dag and rulegraph of Wagtail, and envs directory containing the .yaml file to create necessary conda environment to Wagtail)
3. test (contain python scripts to test the Wagtail functions and test data)
4. data (contain several example data and to put the user's necessary input data).

Users should create several empty directories:
1. database (to save the database directory needed for alignment part)
2. run (to put the result of Wagtail pipeline).

### Output structure
run/run_name/
- 0_logs_wagtail/             # Execution logs and benchmarks
- 6_condensed_wagtail/        # Taxonomy profiles per sample
  - sample_condensed.tsv      # Final taxonomy results
- 7_metadata/                 # Sample metadata
  - sample_metadata.tsv       # Per-sample metadata
- run_name_full_metadata.tsv  # Combined metadata
- [1-5]_*/                    # Intermediate files (cleaned up)

### Key output files
- sample_condensed.tsv: Taxonomic profile with abundances
- run_name_full_metadata.tsv: Combined metadata for all samples
- run_name_YYYYMMDD_log.sql: SQLite database with execution logs

### Executing Wagtail
Please check the config.yaml to configure on how the Wagtail is executed:
1. run_name (running name for this execution. The results will be saved at run/run_name)
2. sample_list (absolute path to sample list .txt)
3. file_map (absolute path to file map .txt).

#### Single run (No batching)

STEP 1: PREPARE INPUT FILES\
Create a file map (e.g., file_map.txt):
```bash
sample001    /absolute/path/to/sample001.fastq.gz
sample002    /absolute/path/to/sample002.fastq.gz
```

Create a sample list (e.g., sample_list.txt):
```bash
sample001
sample002
```

STEP 2: CONFIGURE\
Edit pipeline/config.yaml:
```bash
run_name: "my_analysis"
sample_list: "/path/to/sample_list.txt"
file_map: "/path/to/file_map.txt"
```

STEP 3: RUN WAGTAIL\
Using the wrapper script (recommended)
```bash
python wagtail.py --use-conda -c 8 --configfile pipeline/config.yaml
```

Or directly with snakemake
```bash
snakemake --use-conda -c 8 -s pipeline/wagtail.smk --configfile pipeline/config.yaml
```

STEP 4: CHECK RESULTS\
Results will be in: run/my_analysis/
- 0_logs_wagtail/[run_name]_[date]_log.sql  # full log of each run
- 6_condensed_wagtail/                      # Taxonomy profiles per sample
- 7_metadata/                               # Metadata and combined results

If it is the first time you execute Wagtail, Snakemake will create the conda environment automatically and it will take a few minutes depending on your internet access speed and system.\
All of the result will be saved in the run directory.

#### Batch processing

WHEN TO USE BATCH PROCESSING
- Large datasets (>1000 samples)
- Cluster computing environments
- Need to resume interrupted runs
- Want to process samples in manageable chunks
- Assuming that you have aqua (QUT HPC) environment

STEP 1: PREPARE ACCESSION LIST\
Create a text file with one accession per line:\
```bash
SRR12345678
SRR12345679
SRR12345680
...
```

STEP 2: RUN BATCH SCRIPT

Local execution (testing/small datasets)
```bash
python batch_run_wagtail.py \
  -a accession_list.txt \
  -b 100 \
  -o /path/to/output \
  -l batch_log.log \
  -f /path/to/file_map.txt \
  -p /path/to/wagtail.py \
  -r "my_batch_run" \
  --execution-mode local \
  -c 4 -m 8
```

Cluster execution (production)
```bash
python batch_run_wagtail.py \
  -a accession_list.txt \
  -b 1000 \
  -o /path/to/output \
  -l batch_log.log \
  -f /path/to/file_map.txt \
  -p /path/to/wagtail.py \
  -r "my_batch_run" \
  --execution-mode cluster \
  -c 8 -m 16 -t 48
```

BATCH SCRIPT OPTIONS
INPUT MODES:
```bash
  -a FILE         Accession list file (CSV mode)
  -i DIR          Directory with existing batch files (Directory mode)
```

EXECUTION:
```bash
  --execution-mode local|cluster    Run locally or submit to cluster
  -b SIZE         Batch size (samples per batch)
  -c CORES        CPU cores per batch job
  -m MEM          Memory (GB) per batch job
  -t HOURS        Time limit (hours) per batch job
```

RANDOMIZATION:
```bash
  --randomize     Shuffle accessions before batching
  --random-seed N Use specific random seed
```

RESUME/RANGE:
```bash
  --start-batch N Start from specific batch number
  --end-batch N   End at specific batch number
```

STEP 3: MONITOR PROGRESS
Check batch status
```bash
tail -f batch_log.log
```

Check individual batch status
```bash
ls output_dir/batch_*_run_status.txt
```

Rerun specific batch range
```bash
python batch_run_wagtail.py \
  -i /path/to/existing/batches \
  --start-batch 5 --end-batch 10 \
  [other options...]
```

#### Code examples
MINIMAL TEST RUN: _Test with 3 samples locally_
```bash
python batch_run_wagtail.py \
  -a test_samples.txt \
  -b 3 \
  -o test_output \
  -l test.log \
  -f file_map.txt \
  -p wagtail.py \
  -r "test_run" \
  --execution-mode local \
  -c 2 -m 4
```

PRODUCTION CLUSTER RUN: _Process 10,000 samples in batches of 500_
```bash
python batch_run_wagtail.py \
  -a large_dataset.txt \
  -b 500 \
  -o production_output \
  -l production.log \
  -f file_map.txt \
  -p wagtail.py \
  -r "production_run" \
  --execution-mode cluster \
  --randomize \
  -c 16 -m 32 -t 72
```

RESUME INTERRUPTED RUN: _Continue from batch 10 onwards_
```bash
python batch_run_wagtail.py \
  -i production_output \
  --start-batch 10 \
  --execution-mode cluster \
  -c 16 -m 32 -t 72
```

### Execution result and logs
After successful execution, all results will be saved in run/run_name directory. It will have several directories:
1. 0_logs_wagtail (contain the logs of each pipeline execution)
2. 6_condesed_wagtail (contain the taxonomy profile for each filename/accession)
3. 7_metadata (contains the running metadata for each filename/accession and combined metadata for all filename/accession.

Should the pipeline execution result in error, there will be intermediate results for debugging and error logs in the 0_logs with the name run_name_status.log to pinpoint which filename/accession and on which part of the pipeline return the error.

COMMON ISSUES
1. "conda environment not found"
   → First run takes time to create environments
   → Ensure conda/mamba is properly installed

2. "File not found" errors
   → Check file paths are absolute paths
   → Verify file_map points to existing files

3. Batch jobs fail
   → Check cluster resource limits
   → Reduce batch size (-b parameter)
   → Check disk space in output directory

4. Empty results
   → Verify input FASTQ files are valid
   → Check reference database format
   → Review log files in 0_logs_wagtail/

GETTING HELP
------------
1. Check log files in output directory
2. Review Snakemake documentation
3. Ensure all input files are properly formatted
4. Test with small dataset first

### Tool limitations
Sequence data\
:x: Paired-end reads (use forward reads only)\
:x: Long-read sequencing (Nanopore, PacBio)\
:x: Non-16S amplicon data\
:x: Whole genome shotgun data\
:white_check_mark: Single-end 16S amplicon (Illumina)\
:white_check_mark: Forward reads from paired-end 16S data\

Platform\
✗ Windows (Linux/Unix only)\
:white_check_mark: Linux with Conda/Mamba\
:white_check_mark: HPC clusters with job schedulers\
:white_check_mark: Local workstations\

Database
- Tested primarily with bacterial/archaeal 16S databases
- Custom databases must be in FASTA format
- No built-in database downloading (manual setup required)

## Code development notice
This tool is created as a part of PhD project for analyzing the global amplicon datasets and integrate them in a machine learning model to understand the connection of global microbiome communities with its environment and vice-versa.\
**This tool is only tested for 16S amplicon sequences with bacteria and archaea database (MFD)**
