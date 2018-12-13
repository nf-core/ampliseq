# nf-core/ampliseq: Usage

## Table of contents

* [Introduction](#general-nextflow-info)
* [Running the pipeline](#running-the-pipeline)
* [Updating the pipeline](#updating-the-pipeline)
* [Reproducibility](#reproducibility)
* [Main arguments](#main-arguments)
    * [`-profile`](#-profile-single-dash)
    * [`--reads`](#--reads)
    * [`--FW_primer` and `--RV_primer`](#--FW_primer-and---RV_primer)
    * [`--retain_untrimmed`](#--retain_untrimmed)
    * [`--metadata`](#--metadata)
* [Cutoffs](#cutoffs)
    * [`--trunclenf` and `--trunclenr`](#--trunclenf-and---trunclenr)
    * [`--trunc_qmin`](#--trunc_qmin)
    * [`--untilQ2import`](#--untilQ2import)
    * [`--Q2imported`](#--Q2imported)
* [Reference database](#reference-database)
    * [`--classifier`](#--classifier)
    * [`--classifier_removeHash`](#--classifier_removeHash)
* [Statistics](#statistics)
    * [`--metadata_category`](#--metadata_category)
* [Filters](#filters)
    * [`--exclude_taxa`](#--exclude_taxa)
* [Skipping steps](#skipping-steps)
    * [`--skip_fastqc`](#--skip_fastqc)
    * [`--skip_alpha_rarefaction`](#--skip_alpha_rarefaction)
    * [`--skip_taxonomy`](#--skip_taxonomy)
    * [`--skip_barplot`](#--skip_barplot)
    * [`--skip_abundance_tables`](#--skip_abundance_tables)
    * [`--skip_diversity_indices`](#--skip_diversity_indices)
    * [`--skip_ancom`](#--skip_ancom)
* [Job Resources](#job-resources)
* [Automatic resubmission](#automatic-resubmission)
* [Custom resource requests](#custom-resource-requests)
* [AWS batch specific parameters](#aws-batch-specific-parameters)
    * [`-awsbatch`](#-awsbatch)
    * [`--awsqueue`](#--awsqueue)
    * [`--awsregion`](#--awsregion)
* [Other command line parameters](#other-command-line-parameters)
    * [`--outdir`](#--outdir)
    * [`--email`](#--email)
    * [`-name`](#-name-single-dash)
    * [`-resume`](#-resume-single-dash)
    * [`-c`](#-c-single-dash)
    * [`--max_memory`](#--max_memory)
    * [`--max_time`](#--max_time)
    * [`--max_cpus`](#--max_cpus)
    * [`--plaintext_emails`](#--plaintext_emails)
    * [`--sampleLevel`](#--sampleLevel)
    * [`--multiqc_config`](#--multiqc_config)


## General Nextflow info
Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline
The typical command for running the pipeline is as follows:
```bash
nextflow run nf-core/ampliseq --reads 'data' -profile standard,docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline
When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/ampliseq
```

### Reproducibility
It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/ampliseq releases page](https://github.com/nf-core/ampliseq/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.


## Main Arguments

### `-profile`
Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. Note that multiple profiles can be loaded, for example: `-profile standard,docker` - the order of arguments is important!

* `standard`
    * The default profile, used if `-profile` is not specified at all.
    * Runs locally and expects all software to be installed and available on the `PATH`.
* `docker`
    * A generic configuration profile to be used with [Docker](http://docker.com/)
    * Pulls software from dockerhub: [`nfcore/ampliseq`](http://hub.docker.com/r/nfcore/ampliseq/)
* `singularity`
    * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
    * Pulls software from singularity-hub
* `conda`
    * A generic configuration profile to be used with [conda](https://conda.io/docs/)
    * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `awsbatch`
    * A generic configuration profile to be used with AWS Batch.
* `test`
    * A profile with a complete configuration for automated testing
    * Includes links to test data so needs no other parameters
* `none`
    * No configuration at all. Useful if you want to build your own config from scratch and want to avoid loading in the default `base` config profile (not recommended).

### `--reads`
Use this to specify the location of your input paired-end FastQ files. ***Currently, sequencing data originating from multiple sequencing runs cannot be analysed apppropriately by this pipeline.*** For example:

```bash
--reads 'path/to/data/'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The folder must containing gzip compressed Casava 1.8 paired-end demultiplexed fastq files with the naming sheme *_L001_R{1,2}_001.fastq.gz
3. All sequencing data should originate from one sequencing run, because processing relies on run-specific error models that are unreliable when data from several sequencing runs are mixed.

### `--FW_primer` and `--RV_primer`
In Amplicon sequencing methods, PCR with specific primers produces the amplicon of intrest. These primer sequences need to be trimmed from the reads before further processing and are also required for producing an appropriate classifier. For example:

```bash
--FW_primer GTGYCAGCMGCCGCGGTAA --RV_primer GGACTACNVGGGTWTCTAAT
```

### `--retain_untrimmed`
When read sequences are trimmed, routinely untrimmed read pairs are discarded. Use this option to retain untrimmed read pairs. For example:

```bash
--retain_untrimmed
```

### `--metadata`
For performing downstream analysis such as barplots, diversity indices or differential abundance testing, a metadata file is essential. For example:

```bash
--metadata "$PWD/data/Metadata.tsv"
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The metadata sheet has to follow the [QIIME2 specifications](https://docs.qiime2.org/2018.6/tutorials/metadata/)

## Cutoffs

### `--trunclenf` and `--trunclenr`
Read denoising by DADA2 creates an error profile specific to a sequencing run and uses this to correct sequencing errors. This method requires all reads to have the same length and as high quality as possible while maintaining at least 20 bp overlap for merging. One cutoff for the forward read `--trunclenf` and one for the reverse read `--trunclenr` truncate all longer reads at that position and drop all shorter reads. 
These cutoffs are usually chosen visually using `--untilQ2import`, inspecting the quality plots in the "result folder/demux", and resuming analysis with `--Q2imported [Path]`. If not set, these cutoffs will be determined automatically for the position before the mean quality score drops below `--trunc_qmin`. For example:

```bash
--trunclenf 180 --trunclenr 120
```

Please note:

1. Too agressive truncation might lead to insufficient overlap for read merging
2. Too less truncation might reduce denoised reads
3. The code choosing these values automatically cannot take the points above into account, therefore setting `--trunclenf` and `--trunclenr` is preferred

### `--trunc_qmin`
Automatically determine `--trunclenf` and `--trunclenr` before the mean quality score drops below `--trunc_qmin` (default: 25). Setting values for `--trunclenf` and `--trunclenr` is strongly encouraged. For example:

```bash
--trunc_qmin 20
```

Please note:

1. The code choosing `--trunclenf` and `--trunclenr` using `--trunc_qmin` automatically cannot take amplicon length or overlap requirements for merging into account, therefore setting `--trunclenf` and `--trunclenr` is preferred

### `--untilQ2import`
Computes all steps until quality plots aiding the choosing of `--trunclenf` and `--trunclenr`.

### `--Q2imported`
Analysis starting with a QIIME2 artefact with trimmed reads, typically produced before with `--untilQ2import`.

## Reference database
By default, the workflow downloads [SILVA](https://www.arb-silva.de/) [v132](https://www.arb-silva.de/documentation/release-132/) and extracts reference sequences and taxonomy clustered at 99% similarity and trains a Naive Bayes classifier to assign taxonomy to features.

### `--classifier`
If you have trained a compatible classifier before, from sources such as [SILVA](https://www.arb-silva.de/), [Greengenes](http://greengenes.secondgenome.com/downloads) or [RDP](https://rdp.cme.msu.edu/). For example:

```bash
--classifier "FW_primer-RV_primer-classifier.qza"
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The cassifier is a Naive Bayes classifier produced by "qiime feature-classifier fit-classifier-naive-bayes" (e.g. by this pipeline or from [QIIME2 resources](https://docs.qiime2.org/2018.6/data-resources/))
3. The primer pair for the amplicon PCR and the computing of the classifier are exactly the same
4. The classifier has to be trained by the same version of scikit-learn as this version of the pipeline uses (0.19.1).

### `--classifier_removeHash`
Remove all hash signs from taxonomy strings, resolves a rare ValueError during classification (process classifier).


## Statistics

### `--metadata_category`
Here columns in the metadata sheet can be chosen with groupings that are used for diversity indices and differential abundance analysis. By default, all suitable columns in the metadata sheet will be used if this option is not specified. Suitable are columns which are categorical (not numerical) and have multiple different values which are not all unique. For example:

```bash
--metadata_category "treatment1,treatment2"
```

Please note the following requirements:

1. Comma seperated list enclosed in quotes and may not contain whitespace characters
2. Each comma seperated term has to match exactly one column name in the metadata sheet

## Filters

### `--exclude_taxa`
Depending on the primers used, PCR might amplify unwanted or off-target DNA, e.g. originating from mitochondria or chloroplasts. Here you can specify taxa that are excluded from further analysis. For example:

```bash
--exclude_taxa "mitochondria,chloroplast"
```

Please note the following requirements:

1. Comma seperated list enclosed in quotes and may not contain whitespace characters
2. Features that contain one or several of these terms in their taxonomical classification are excluded from further analysis
3. The taxonomy level is not taken into consideration

## Skipping steps:

### `--skip_fastqc`
Skip FastQC, minor time saving.

### `--skip_alpha_rarefaction`
Skip alpha rarefaction, minor time saving.

### `--skip_taxonomy`
Skip taxonomic classification, will essentially truncate the workflow after denoising.

### `--skip_barplot`
Skip producing barplot, minor time saving.

### `--skip_abundance_tables`
Skip producing most relative abundance tables, minor time saving.

### `--skip_diversity_indices`
Skip alpha and beta diversity analysis, large time saving.

### `--skip_ancom`
Skip differential abundance testing, large time saving.

## Job Resources
### Automatic resubmission
Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests
Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files in [`conf`](../conf) for examples.

## AWS Batch specific parameters
Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use the `-awsbatch` profile and then specify all of the following parameters.
### `--awsqueue`
The JobQueue that you intend to use on AWS Batch.
### `--awsregion`
The AWS region to run your job in. Default is set to `eu-west-1` but can be adjusted to your needs.

Please make sure to also set the `-w/--work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't.

## Other command line parameters

### `--outdir`
The output directory where the results will be saved.

### `--email`
Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to speicfy this on the command line for every run.

### `-name`
Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`
Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override defaults. For example, you can specify a config file using `-c` that contains the following:

```nextflow
process.$multiqc.module = []
```

### `--max_memory`
Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'``

### `--max_time`
Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`
Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--plaintext_email`
Set to receive plain-text e-mails instead of HTML formatted.

### `--multiqc_config`
Specify a path to a custom MultiQC configuration file.
