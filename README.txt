===============================================================================
                                   WAGTAIL
                    16S Amplicon Analysis Pipeline (Version 0.14)
===============================================================================

OVERVIEW
--------
Wagtail is a Snakemake-based pipeline for analyzing 16S amplicon sequences.
It combines QIIME2 Deblur for quality control and Minimap2 for taxonomic 
classification against custom reference databases.

KEY FEATURES
- Accurate taxonomic profiling using custom databases
- Scalable batch processing for large datasets
- Quality control and denoising with QIIME2 Deblur
- Fast alignment with Minimap2/Mappy
- Comprehensive logging and metadata generation

===============================================================================
                              QUICK START GUIDE
===============================================================================

PREREQUISITES
-------------
- Snakemake installed
- Conda/Mamba for environment management
- Reference database in FASTA format (e.g., danica.mmi in database/ folder)

INPUT REQUIREMENTS
------------------
1. FASTQ files (.fastq or .fastq.gz) - SINGLE-END ONLY
2. File map (2 columns: sample_name, absolute_path_to_file)
3. Sample list (1 column: sample_names matching file map)
4. Reference database (FASTA format in database/ directory)

DIRECTORY STRUCTURE
-------------------
wagtail/
├── wagtail.py          # Main wrapper script
├── pipeline/
│   ├── wagtail.smk     # Snakemake pipeline
│   └── config.yaml     # Configuration file
├── batch_run_wagtail.py  # Batch processing script
├── database/           # Reference databases (create this)
├── data/              # Input data
└── run/               # Output results (create this)

===============================================================================
                        METHOD 1: SINGLE RUN (No Batching)
===============================================================================

STEP 1: PREPARE INPUT FILES
---------------------------
Create a file map (e.g., file_map.txt):
sample001    /absolute/path/to/sample001.fastq.gz
sample002    /absolute/path/to/sample002.fastq.gz

Create a sample list (e.g., sample_list.txt):
sample001
sample002

STEP 2: CONFIGURE
-----------------
Edit pipeline/config.yaml:
run_name: "my_analysis"
sample_list: "/path/to/sample_list.txt"
file_map: "/path/to/file_map.txt"

STEP 3: RUN WAGTAIL
-------------------
# Using the wrapper script (recommended)
python wagtail.py --use-conda -c 8 --configfile pipeline/config.yaml

# Or directly with snakemake
snakemake --use-conda -c 8 -s pipeline/wagtail.smk --configfile pipeline/config.yaml

STEP 4: CHECK RESULTS
---------------------
Results will be in: run/my_analysis/
- 6_condensed_wagtail/    # Taxonomy profiles per sample
- 7_metadata/             # Metadata and combined results

===============================================================================
                        METHOD 2: BATCH PROCESSING
===============================================================================

WHEN TO USE BATCH PROCESSING
-----------------------------
- Large datasets (>1000 samples)
- Cluster computing environments
- Need to resume interrupted runs
- Want to process samples in manageable chunks
- Assuming that you have aqua (QUT HPC) environment

STEP 1: PREPARE ACCESSION LIST
------------------------------
Create a text file with one accession per line:
SRR12345678
SRR12345679
SRR12345680
...

STEP 2: RUN BATCH SCRIPT
------------------------

# Local execution (testing/small datasets)
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

# Cluster execution (production)
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

BATCH SCRIPT OPTIONS
--------------------
INPUT MODES:
  -a FILE         Accession list file (CSV mode)
  -i DIR          Directory with existing batch files (Directory mode)

EXECUTION:
  --execution-mode local|cluster    Run locally or submit to cluster
  -b SIZE         Batch size (samples per batch)
  -c CORES        CPU cores per batch job
  -m MEM          Memory (GB) per batch job
  -t HOURS        Time limit (hours) per batch job

RANDOMIZATION:
  --randomize     Shuffle accessions before batching
  --random-seed N Use specific random seed

RESUME/RANGE:
  --start-batch N Start from specific batch number
  --end-batch N   End at specific batch number

STEP 3: MONITOR PROGRESS
------------------------
# Check batch status
tail -f batch_log.log

# Check individual batch status
ls output_dir/batch_*_run_status.txt

# Rerun specific batch range
python batch_run_wagtail.py \
  -i /path/to/existing/batches \
  --start-batch 5 --end-batch 10 \
  [other options...]

===============================================================================
                              INPUT/OUTPUT DETAILS
===============================================================================

INPUT FILES
-----------
1. FASTQ Files:
   - Single-end reads only
   - .fastq or .fastq.gz format
   - Forward reads only (if paired-end available)

2. File Map:
   - Tab-separated: sample_name<TAB>absolute_path
   - Must match sample names in sample list

3. Sample List:
   - One sample name per line
   - Must match first column in file map

4. Reference Database:
   - FASTA format in database/ directory
   - Pre-indexed .mmi files supported for speed

OUTPUT STRUCTURE
----------------
run/run_name/
├── 0_logs_wagtail/           # Execution logs and benchmarks
├── 6_condensed_wagtail/      # Taxonomy profiles per sample
│   └── sample_condensed.tsv  # Final taxonomy results
├── 7_metadata/               # Sample metadata
│   ├── sample_metadata.tsv   # Per-sample metadata
│   └── run_name_full_metadata.tsv  # Combined metadata
└── [1-5]_*/                  # Intermediate files (cleaned up)

KEY OUTPUT FILES
----------------
- sample_condensed.tsv: Taxonomic profile with abundances
- run_name_full_metadata.tsv: Combined metadata for all samples
- run_name_YYYYMMDD_log.sql: SQLite database with execution logs

===============================================================================
                              TOOL LIMITATIONS
===============================================================================

SEQUENCE DATA LIMITATIONS
--------------------------
✗ Paired-end reads (use forward reads only)
✗ Long-read sequencing (Nanopore, PacBio)
✗ Non-16S amplicon data
✗ Whole genome shotgun data
✓ Single-end 16S amplicon (Illumina)
✓ Forward reads from paired-end 16S data

PLATFORM LIMITATIONS
---------------------
✗ Windows (Linux/Unix only)
✓ Linux with Conda/Mamba
✓ HPC clusters with job schedulers
✓ Local workstations

DATABASE LIMITATIONS
--------------------
- Tested primarily with bacterial/archaeal 16S databases
- Custom databases must be in FASTA format
- No built-in database downloading (manual setup required)

COMPUTATIONAL REQUIREMENTS
---------------------------
- Minimum: 4 GB RAM, 2 CPU cores
- Recommended: 16+ GB RAM, 8+ CPU cores
- Storage: ~10-50 GB per 1000 samples (depending on sample size)
- Network: Required for conda environment setup

===============================================================================
                               TROUBLESHOOTING
===============================================================================

COMMON ISSUES
-------------
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

===============================================================================
                                  EXAMPLES
===============================================================================

MINIMAL TEST RUN
----------------
# Test with 3 samples locally
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

PRODUCTION CLUSTER RUN
----------------------
# Process 10,000 samples in batches of 500
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

RESUME INTERRUPTED RUN
----------------------
# Continue from batch 10 onwards
python batch_run_wagtail.py \
  -i production_output \
  --start-batch 10 \
  --execution-mode cluster \
  -c 16 -m 32 -t 72