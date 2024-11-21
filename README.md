![Wagtail Logo](https://github.com/R-Nurdiansyah/Wagtail/blob/development/wagtail_logo_(29-7-2024).png?raw=true)

# Wagtail
_Version 0.05 (21-11-2024)_

Wagtail is an accurate and scalable tools to analyze 16S dataset with easy-to-swap reference database. As for now, the functionality is **limited to 16S amplicon sequences sequenced from Illumina sequencer only**.\
The tool is a combination of Qiime 2 [Deblur](https://github.com/biocore/deblur/blob/master/README.md) plugin for quality control and Minimap2 aligner to align representative sequences from deblur. The result is a taxonomy profile in the format of condense and running metadata. The combination is weaved using [Snakemake](https://snakemake.github.io/) to allow easy reproducibility, benchmarking, and packaging.\
[Minimap2](https://github.com/lh3/minimap2), in the form of its python interface (Mappy), is chosen to its versatility in using reference database without/with minor modification to the reference database. User only need to give the database in the fasta form or use the Minimap2 to help create the index for faster reading (optional).

## How to use
Wagtail is built on Snakemake so executing the tool should be based on how Snakemake is executed. Make sure Snakemake is installed in the system before using it. 
Please check [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/index.html) for further reading in Snakemake execution.

### Directories
In the Wagtail directory, you will have several directories:
1. bin (contain python scripts important to the Wagtail functionality)
2. pipeline (contain Wagtail snakemake (Wagtail.smk), config.yaml, dag and rulegraph of Wagtail, and envs directory containing the .yaml file to create necessary conda environment to Wagtail)
3. test (contain python scripts to test the Wagtail functions and test data)
4. data (contain several example data and to put the user's necessary input data).

Users should create several empty directories:
1. database (to save the database directory needed for alignment part)
2. run (to put the result of Wagtail pipeline).

### Input
Wagtail needs 4 kind of inputs to work:
1. Database. Before executing the pipeline, please put your desired database in the database in the fasta format or in indexed format.
2. Input data in fastq/fastq.gz, recommended to put inside the data directory. Recommended format: filename/accession-1.fastq.gz **Wagtail only need the first-end of the read to work and does not support paired-end. Should your data is paired-end, then use the forward-end only**
3. file_map, contain 2 column with filename/accession in the first column and the absolute path to the file in the second column
4. sample_list, contain only filename/accession.

Make sure that the no 3 and no 4 have the same filename/accession. We provide the filemap creation script in the bin directory and example input data, filemap, and sample list in data directory for your convenience.

### Executing Wagtail

Please check the config.yaml to configure on how the Wagtail is executed:
1. run_name (running name for this execution. The results will be saved at run/run_name)
2. sample_list (absolute path to sample list .txt)
3. file_map (absolute path to file map .txt).

Example execution command:
use conda, 8 cores, and edited the config file
> snakemake --use-conda -c 8 --snakefile path/to/wagtail.smk --configfile=path/to/config.yaml

use conda, 8 cores, and edit the config file in the terminal
> snakemake --use-conda -c 8 --snakefile path/to/wagtail.smk --configfile=path/to/config.yaml --config run_name=user_configured_run_name file_map=path/to/filemap sample_list=path/to/sample_list

If it is the first time you execute Wagtail, Snakemake will create the conda environment automatically and it will take a few minutes depending on your internet access speed and system.\
All of the result will be saved in the run directory.

### Execution result and logs
After successful execution, all results will be saved in run/run_name directory. It will have several directories:
1. 0_logs_wagtail (contain the logs of each pipeline execution)
2. 6_condesed_wagtail (contain the taxonomy profile for each filename/accession)
3. 7_metadata (contains the running metadata for each filename/accession and combined metadata for all filename/accession.

Should the pipeline execution result in error, there will be intermediate results for debugging and error logs in the 0_logs with the name run_name_status.log to pinpoint which filename/accession and on which part of the pipeline return the error.

## Code development notice
This tool is created as a part of PhD project for analyzing the global amplicon datasets and integrate them in a machine learning model to understand the connection of global microbiome communities with its environment and vice-versa.\
**This tool is only tested for 16S amplicon sequences with bacteria and archaea database (MFD)**
