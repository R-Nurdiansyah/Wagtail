![Wagtail Logo](https://github.com/R-Nurdiansyah/Wagtail/blob/development/wagtail_logo_(29-7-2024).png?raw=true)

# Wagtail
_Version 1.1_

Wagtail is an accurate and scalable tool to analyze amplicon sequencing datasets with easy-to-swap reference databases. As of version 1.1, Wagtail supports **16S, 18S, ITS, and CO1** amplicon sequences sequenced from the Illumina platform.\
The tool is a combination of Qiime 2 [Deblur](https://library.qiime2.org/plugins/qiime2/q2-deblur/overview) plugin for quality control and [Minimap2](https://github.com/lh3/minimap2) aligner to align representative sequences from Deblur. The result is a taxonomy profile in the format of condensed and running metadata. The combination is weaved using [Snakemake](https://snakemake.github.io/) to allow easy reproducibility, benchmarking, and packaging.\
[Minimap2](https://github.com/lh3/minimap2), in the form of its Python interface (Mappy), is chosen for its versatility in using reference databases without or with minor modification. Users only need to provide the database in FASTA format, or use the Minimap2 index format for faster loading (optional).

KEY FEATURES
- **Automatic marker identification**: pre-flight QC step identifies the amplicon type (16S / 18S / ITS / CO1) directly from raw reads before processing begins
- Multi-marker support: 16S, 18S, ITS, and CO1 in a single pipeline
- Accurate taxonomic profiling using custom, marker-specific databases
- Scalable batch processing for large datasets
- Quality control and denoising with QIIME2 Deblur
- Fast alignment with Minimap2/Mappy
- Soft-fail system: samples that fail at any step are logged and skipped gracefully — the pipeline continues for all remaining samples
- Comprehensive per-sample SQLite logging and metadata generation

## How to use
Wagtail is built on Snakemake so executing the tool should be based on how Snakemake is executed. Make sure Snakemake is installed in the system before using it.
Please check [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/index.html) for further reading on Snakemake execution.

### Prerequisites
- Snakemake installed
- Conda/Mamba for environment management
- Reference databases in FASTA format (see Database Setup below)

### Database Setup

Wagtail v0.15 requires marker-specific databases placed in the `database/` directory. Databases are resolved by filename prefix, so each file must start with the marker name.

**Marker identification databases** (used by the `marker_id` pre-flight step — all four are required):

| Marker | Filename pattern | Example |
|--------|-----------------|---------|
| 16S | `16S_database*` | `16S_database_MFD_NR987.fasta.gz` |
| 18S | `18S_database*` | `18S_database_pr2_v5.fasta.gz` |
| ITS | `ITS_database*` | `ITS_database_sh_refs_qiime.fasta.gz` |
| CO1 | `CO1_database*` | `CO1_database_MIDORI2_CO1.fasta.gz` |

**Deblur reference databases** (used during denoising — required for non-16S markers only; 16S uses QIIME2's built-in reference):

| Marker | Filename pattern | Example |
|--------|-----------------|---------|
| 18S | `18S_deblur*` | `18S_deblur_pr2.fasta.gz` |
| ITS | `ITS_deblur*` | `ITS_deblur_sh_refs.fasta.gz` |
| CO1 | `CO1_deblur*` | `CO1_deblur_MIDORI2.fasta.gz` |

**Taxonomy databases** (used during alignment — one per marker):

| Marker | Filename pattern | Example |
|--------|-----------------|---------|
| 16S | `16S_taxonomy*` | `16S_taxonomy_MFD.fasta.gz` |
| 18S | `18S_taxonomy*` | `18S_taxonomy_pr2.fasta.gz` |
| ITS | `ITS_taxonomy*` | `ITS_taxonomy_sh_refs.fasta.gz` |
| CO1 | `CO1_taxonomy*` | `CO1_taxonomy_MIDORI2.fasta.gz` |

> Each pattern must match **exactly one file** in the `database/` directory.

#### Taxonomy file format

Taxonomy reference files (`*_taxonomy*`) must be **uncompressed, tab-separated** files with two columns: feature ID and semicolon-delimited taxonomy string. No header row.

**Amplicon markers (16S, 18S, CO1)** — 7 core levels with `d__` or `k__` prefix notation:
```
FLASV1.1346	d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Propionibacteriales; f__Propionibacteriaceae; g__Cutibacterium; s__Cutibacterium_acnes
FLASV2.1310	d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Xanthobacteraceae; g__Bradyrhizobium; s__MFD_s_2
FLASV3.1372	d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Bacillaceae; g__Bacillus; s__MFD_s_3
```

**ITS (UNITE)** — 7 core levels starting with `k__` plus an optional `sh__` (species hypothesis; added by [UNITE](https://unite.ut.ee/index.php#main)) rank:
```
SH1227328.10FU_MT153946_refs	k__Fungi; p__Ascomycota; c__Dothideomycetes; o__Abrothallales; f__Abrothallaceae; g__Abrothallus; s__Abrothallus_subhalei; sh__SH1227328.10FU
SH1227742.10FU_JN206177_refs	k__Fungi; p__Mucoromycota; c__Mucoromycetes; o__Mucorales; f__Mucoraceae; g__Mucor; s__Mucor_inaequisporus; sh__SH1227742.10FU
SH1232203.10FU_KY102517_refs	k__Fungi; p__Ascomycota; c__Saccharomycetes; o__Saccharomycetales; f__Saccharomycetales_fam_Incertae_sedis; g__Candida; s__Candida_vrieseae; sh__SH1232203.10FU
```

`--check-db` reports the `sh__` rank explicitly in the level summary, e.g.:
```
Rows: 96606  |  Levels: 8 levels (7 core + sh__)
Optional ranks present: sh__ (e.g. UNITE species hypothesis)
✓ PASS — all 96606 rows OK
```

#### Creating deblur reference `.qza` files

The deblur reference databases for non-16S markers must be QIIME2 artifacts (`.qza`). Create them from a plain FASTA file using the standard [QIIME2](https://docs.qiime2.org/2024.10/tutorials/importing/) import command, for example:

```bash
qiime tools import \
    --type 'FeatureData[Sequence]' \
    --input-path  18S_reference.fasta \
    --output-path 18S_deblur_pr2.qza
```

Place the resulting `.qza` files in the `database/` directory. The filenames must match the `{marker}_deblur*` pattern (e.g. `18S_deblur_pr2.qza`). 16S does not need a deblur reference. QIIME2 [q2-deblur](https://library.qiime2.org/plugins/qiime2/q2-deblur/overview) `denoise-16S` uses a built-in reference automatically.

### Input

Wagtail needs 4 kinds of inputs to work:
1. **Databases** — see Database Setup above.
2. **Input data** in `fastq` / `fastq.gz` format. Wagtail only needs the first (forward) read and does not support paired-end. If your data is paired-end, use the forward read only.
3. **file_map** — a headerless 2-column TSV: accession in column 1, absolute path to the FASTQ file in column 2:
   ```
   sample001    /absolute/path/to/sample001.fastq.gz
   sample002    /absolute/path/to/sample002.fastq.gz
   ```
4. **sample_list** — a plain text file with one accession per line. The accessions must match those in the file_map.

### Directories

In the Wagtail directory, you will have several directories:
1. `bin/` — Python scripts that implement Wagtail functionality
2. `pipeline/` — the Snakemake file (`wagtail.smk`), `config.yaml`, shared bash helpers (`wagtail_functions.sh`), and the `envs/` directory with conda environment definitions
3. `test/` — unit tests for Wagtail scripts and test data
4. `data/` — example input data, filemap, and sample list

Users should create the following empty directories before the first run:
1. `database/` — to hold the reference databases
2. `run/` — to receive pipeline outputs

### Pipeline overview

The pipeline runs the following steps in order for each sample:

```
marker_id → manifest → import_and_qc → deblur
                                           ↓
                              ┌────────────┴──────────────┐
                           export_all                 qiime_stats
                           ↙         ↘                    ↓
                        mappy     (table_edit)    metadata_creation
                          ↓            ↓          ↑
                     extract_taxonomy ─┘          │
                              ↓                   │ (mappy metadata)
                              └───────────────────┘
                                       ↓
                                    cleanup
```

| Step | Description |
|------|-------------|
| `marker_id` | Sub-samples 1,000 reads and competitively aligns against all four databases to identify the amplicon type. Runs as a standalone job before the main pipeline group. Soft-fails if the marker is ambiguous or unrecognised. |
| `manifest` | Creates a QIIME2-format sample manifest from the file_map. |
| `import_and_qc` | Imports the raw FASTQ into a QIIME2 artifact, then immediately quality-filters reads using `qiime quality-filter q-score`. The intermediate import artifact is kept in scratch space and never written to shared storage. |
| `deblur` | Denoises reads: `denoise-16S` for 16S samples (uses QIIME2's built-in reference), `denoise-other` with an explicit reference for 18S / ITS / CO1. Wall-time scales automatically with retries (4 h → 8 h → 12 h). |
| `export_all` | Exports representative ASV sequences as FASTA **and** exports the BIOM frequency table to TSV in a single job, eliminating one scheduler round-trip per sample. |
| `mappy` | Aligns representative sequences to the marker-specific taxonomy database. |
| `extract_taxonomy` | Combines alignment and frequency table into a condensed taxonomy profile. |
| `qiime_stats` | Extracts QC and deblur statistics from QIIME2 artifacts. |
| `metadata_creation` | Combines QC, deblur, and alignment metadata into a per-sample TSV. |
| `cleanup` | Merges all per-sample SQLite logs into one, combines per-sample metadata, and removes intermediate files. |

**Soft-fail behaviour**: if a sample fails at any step (e.g. low quality, wrong marker, OOM), the failure is recorded in a per-sample SQLite log. All downstream steps for that sample are skipped automatically, but the pipeline continues processing all other samples to completion.

### Output structure

```
run/run_name/
├── 0_logs_wagtail/             # Execution logs, benchmarks, and SQLite logs
│   └── {sample}/
│       ├── marker_id.log
│       ├── import_and_qc.log
│       ├── deblur.log
│       ├── export_all.log
│       └── *.benchmark.log
├── 6_condensed_wagtail/        # Taxonomy profiles per sample
│   └── {sample}_condensed.tsv
├── 7_metadata/                 # Sample metadata
│   ├── {sample}/
│   │   └── {sample}_metadata.tsv
│   └── {run_name}_full_metadata.tsv   # Combined metadata for all samples
└── {run_name}_YYYYMMDD_log.sql        # SQLite execution log
```

Intermediate directories (`1_manifest/` through `5_taxonomy_wagtail/`) are cleaned up after the pipeline completes. The intermediate QIIME2 import artifact (previously `1_import_wagtail/`) is now written to per-sample scratch space (`TMPDIR`) and is never materialised on shared storage.

### Key output files
- `{sample}_condensed.tsv` — taxonomic profile with abundances for one sample
- `{run_name}_full_metadata.tsv` — combined metadata for all samples in the run
- `{run_name}_YYYYMMDD_log.sql` — SQLite database with per-sample execution status for every pipeline step

### Executing Wagtail

Check `pipeline/config.yaml` before running:

```yaml
run_name: "my_analysis"
sample_list: "/absolute/path/to/sample_list.txt"
file_map: "/absolute/path/to/file_map.txt"
```

#### Single run (no batching)

**STEP 1: PREPARE INPUT FILES**

Create a file map (`file_map.txt`):
```
sample001    /absolute/path/to/sample001.fastq.gz
sample002    /absolute/path/to/sample002.fastq.gz
```

Create a sample list (`sample_list.txt`):
```
sample001
sample002
```

[!WARNING]
deblur cannot work with sample name with underscores, make sure to avoid underscores in the sample name

**STEP 2: CONFIGURE**

Edit `pipeline/config.yaml`:
```yaml
run_name: "my_analysis"
sample_list: "/path/to/sample_list.txt"
file_map: "/path/to/file_map.txt"
```

**STEP 3: RUN WAGTAIL**

Using the wrapper script (recommended):
```bash
python wagtail.py --use-conda -c 8 --configfile pipeline/config.yaml
```

Or with the GreenGenes2 database variant:
```bash
python wagtail.py --gg2 --use-conda -c 8 --configfile pipeline/config.yaml
```

Or directly with Snakemake:
```bash
snakemake --use-conda -c 8 -s pipeline/wagtail.smk --configfile pipeline/config.yaml
```

**STEP 4: CHECK RESULTS**

Results will be in `run/my_analysis/`:
```
0_logs_wagtail/[run_name]_[date]_log.sql   # full execution log
0_marker_id/                               # marker identification results
6_condensed_wagtail/                       # taxonomy profiles per sample
7_metadata/                                # metadata and combined results
```

If it is the first time you execute Wagtail, Snakemake will create the conda environments automatically, which may take a few minutes depending on your internet speed and system.

#### Batch processing

**WHEN TO USE BATCH PROCESSING**
- Large datasets (>1000 samples)
- Cluster computing environments
- Need to resume interrupted runs
- Want to process samples in manageable chunks
- Assuming that you have an HPC (e.g. QUT Aqua) environment

**STEP 1: PREPARE SAMPLE NAME/ACCESSION LIST**

Create a text file with one sample name/accession per line:
```
SRR12345678
SRR12345679
SRR12345680
```

**STEP 2: RUN BATCH SCRIPT**

Local execution (testing/small datasets):
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

Cluster execution (production):
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

**BATCH SCRIPT OPTIONS**

Input modes:
```
  -a FILE         Accession list file
  -i DIR          Directory with existing batch files (resume mode)
```

Execution:
```
  --execution-mode local|cluster    Run locally or submit to cluster
  -b SIZE         Batch size (samples per batch)
  -c CORES        CPU cores per batch job
  -m MEM          Memory (GB) per batch job
  -t HOURS        Time limit (hours) per batch job
```

Randomization:
```
  --randomize     Shuffle accessions before batching
  --random-seed N Use specific random seed
```

Resume/range:
```
  --start-batch N Start from specific batch number
  --end-batch N   End at specific batch number
```

**STEP 3: MONITOR PROGRESS**
```bash
tail -f batch_log.log
ls output_dir/batch_*_run_status.txt
```

Resume a specific batch range:
```bash
python batch_run_wagtail.py \
  -i /path/to/existing/batches \
  --start-batch 5 --end-batch 10 \
  [other options...]
```

#### Code examples

Minimal test run (3 samples locally):
```bash
python batch_run_wagtail.py \
  -a test_samples.txt -b 3 \
  -o test_output -l test.log \
  -f file_map.txt -p wagtail.py \
  -r "test_run" --execution-mode local \
  -c 2 -m 4
```

Production cluster run (10,000 samples in batches of 500):
```bash
python batch_run_wagtail.py \
  -a large_dataset.txt -b 500 \
  -o production_output -l production.log \
  -f file_map.txt -p wagtail.py \
  -r "production_run" --execution-mode cluster \
  --randomize -c 16 -m 32 -t 72
```

Resume interrupted run from batch 10:
```bash
python batch_run_wagtail.py \
  -i production_output --start-batch 10 \
  --execution-mode cluster \
  -c 16 -m 32 -t 72
```

### Execution result and logs

After successful execution, all results will be saved in `run/run_name/`. It will contain several directories:
1. `0_logs_wagtail/` — execution logs and benchmarks for each step; the merged SQLite database (`run_name_YYYYMMDD_log.sql`) records the status of every sample at every step
2. `6_condensed_wagtail/` — taxonomy profile for each sample
3. `7_metadata/` — running metadata for each sample and a combined metadata file for all samples

Should the pipeline execution result in errors, the SQLite log at `0_logs_wagtail/run_name_YYYYMMDD_log.sql` records exactly which sample failed at which step. The affected sample's downstream steps are skipped while all other samples continue normally.

### Database management

`wagtail.py` includes built-in commands for checking that all required database files are present and correctly formatted before running a batch. These are especially useful when setting up a new environment or adding a new marker.

#### Check databases (`--check-db`)

Scans the `database/` directory for all required files (identification, deblur reference, and taxonomy for each marker), verifies each file exists and is non-empty, and then runs the taxonomy format validator on any `.tsv` taxonomy files it finds.

```bash
python wagtail.py --check-db
```

Override the database directory:
```bash
python wagtail.py --check-db --db-dir /path/to/database
```

Example output:
```
Checking databases in: /path/to/database

  [16S] identification    OK  16S_database_MFD_NR987.fasta.gz
  [16S] deblur_ref        OK  (built-in QIIME2 reference)
  [16S] taxonomy          OK  16S_taxonomy_MFD_fixed.tsv
  [18S] identification    OK  18S_database_pr2_v5.fasta.gz
  [18S] deblur_ref        OK  18S_deblur_pr2.qza
  [18S] taxonomy          OK  18S_taxonomy_pr2_fixed.tsv
  ...

Validating taxonomy format (4 file(s))...
  ✓ PASS — all rows OK
  ...

======================================================================
  PASS — all databases present and taxonomy formats valid
======================================================================
```

Exit code 0 = all checks passed; exit code 1 = at least one problem found. This makes it suitable for use in pre-run CI or cluster submission scripts.

> **Note:** Taxonomy format validation only runs on uncompressed `.tsv` files. If your taxonomy database is in a compressed or non-TSV format, run `bin/check_db_taxonomy.py` directly.

#### List databases (`--list-db`)

Prints a formatted table of every database file detected in the `database/` directory without performing any validation. Useful for a quick sanity check before a run.

```bash
python wagtail.py --list-db
python wagtail.py --list-db --db-dir /path/to/database
```

Example output:
```
Database directory: /path/to/database

  Marker   Type             File
  ------   ----------------   --------------------------------------------------
  16S      identification   16S_database_MFD_NR987.fasta.gz
  16S      deblur_ref       (built-in QIIME2 reference)
  16S      taxonomy         16S_taxonomy_MFD_fixed.tsv
  18S      identification   18S_database_pr2_v5.fasta.gz
  18S      deblur_ref       18S_deblur_pr2.qza
  18S      taxonomy         18S_taxonomy_pr2_fixed.tsv
  ITS      identification   ITS_database_sh_refs.fasta.gz
  ...
```

#### Print version (`--version`)

```bash
python wagtail.py --version
# Wagtail 0.15
```

### Troubleshooting

**AMBIGUOUS marker result**\
If `0_marker_id/{sample}.marker` contains `status=AMBIGUOUS`, the reads aligned to two databases with nearly identical scores. Check whether your sample is truly 16S/18S/ITS/CO1. The identity margin reported in the marker file indicates how close the top two hits were.

**UNKNOWN marker result**\
The best hit fraction was below the confidence threshold (default 0.60). The sample may not be a supported amplicon type, or the reads may be of very low quality.

**COMMON ISSUES**

1. "conda environment not found"\
   → First run takes time to create environments\
   → Ensure conda/mamba is properly installed

2. "No file matching pattern" errors\
   → Run `python wagtail.py --check-db` to get a clear report of what is missing\
   → Verify database filenames match the expected patterns (e.g. `16S_database*`)\
   → Each pattern must match exactly one file in `database/`

3. "File not found" or "Sample not in filemap"\
   → Check that all accessions in `sample_list.txt` have a matching entry in `file_map.txt`\
   → Verify file paths in the filemap are absolute and the files exist

4. Batch jobs fail / OOM\
   → Check cluster resource limits; `marker_id` requires 8 GB in auto-detect mode\
   → `deblur` wall-time scales automatically on retry (4 h → 8 h → 12 h); ensure `--restart-times 2` is set\
   → Reduce batch size (`-b` parameter)\
   → Check disk space in the output directory

5. Empty results\
   → Verify input FASTQ files are valid and non-empty\
   → Check the reference database format\
   → Review log files in `0_logs_wagtail/{sample}/`

GETTING HELP
1. Run `python wagtail.py --check-db` to verify all databases are present and correctly formatted
2. Check log files in `0_logs_wagtail/` and the SQLite run log
3. Review Snakemake documentation
4. Ensure all input files are properly formatted
5. Test with a small dataset first

### Tool limitations

Sequence data\
:white_check_mark: Single-end 16S amplicon (Illumina)\
:white_check_mark: Single-end 18S amplicon (Illumina)\
:white_check_mark: Single-end ITS amplicon (Illumina)\
:white_check_mark: Single-end CO1 amplicon (Illumina)\
:white_check_mark: Forward reads from paired-end data\
:x: Paired-end reads (use forward reads only)\
:x: Long-read sequencing (Nanopore, PacBio)\
:x: Whole genome shotgun data

Platform\
:x: Windows (Linux/Unix only)\
:white_check_mark: Linux with Conda/Mamba\
:white_check_mark: HPC clusters with job schedulers\
:white_check_mark: Local workstations

Database\
- Tested with bacterial/archaeal 16S ([MFD](https://zenodo.org/records/17162544), [GreenGenes2](https://greengenes2.ucsd.edu/)), 18S [(PR2)](https://pr2-database.org/), ITS [(UNITE)](https://unite.ut.ee/index.php#main), and CO1 [(MIDORI2)](https://www.reference-midori.info/) databases\
- Custom databases must be in FASTA format (gzipped or plain)\
- Filename must match the marker prefix pattern (e.g. `16S_database*`)\
- No built-in database downloading or editing. User should procure and edit manually

## Code development notice

This tool is created as part of a PhD project for analysing global amplicon datasets and integrating them into a machine learning model to understand the connection of global microbiome communities with their environment.\
As of v0.15, Wagtail has been tested with 16S (bacteria/archaea, MFD and GreenGenes2 databases), 18S (PR2), ITS (UNITE), and CO1 (MIDORI2) amplicon datasets.
