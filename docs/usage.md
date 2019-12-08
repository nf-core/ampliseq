# nf-core/ampliseq: Usage

## Table of contents

* [nf-core/ampliseq: Usage](#nf-coreampliseq-usage)
  * [Table of contents](#table-of-contents)
  * [General Nextflow info](#general-nextflow-info)
  * [Running the pipeline](#running-the-pipeline)
    * [Updating the pipeline](#updating-the-pipeline)
    * [Reproducibility](#reproducibility)
  * [Main arguments](#main-arguments)
    * [-profile](#profile)
    * [--reads](#reads)
    * [--FW_primer and --RV_primer](#fwprimer-and---rvprimer)
    * [--metadata](#metadata)
  * [Other input options](#other-input-options)
    * [--extension](#extension)
    * [--multipleSequencingRuns](#multiplesequencingruns)
    * [--split](#split)
    * [--phred64](#phred64)
  * [Cutoffs](#cutoffs)
    * [--trunclenf and --trunclenr](#trunclenf-and---trunclenr)
    * [--trunc_qmin](#truncqmin)
  * [Other options](#other-options)
    * [--untilQ2import](#untilq2import)
    * [--Q2imported](#q2imported)
    * [--keepIntermediates](#keepintermediates)
      * [Visually choosing sequencing read truncation cutoffs](#visually-choosing-sequencing-read-truncation-cutoffs)
  * [Reference database](#reference-database)
    * [--classifier](#classifier)
    * [--classifier_removeHash](#classifierremovehash)
  * [Statistics](#statistics)
    * [--metadata_category](#metadatacategory)
  * [Filters](#filters)
    * [--retain_untrimmed](#retainuntrimmed)
    * [--exclude_taxa](#excludetaxa)
    * [--min_frequency](#minfrequency)
    * [--min_samples](#minsamples)
  * [Skipping steps](#skipping-steps)
    * [--onlyDenoising](#onlydenoising)
    * [--skip_fastqc](#skipfastqc)
    * [--skip_alpha_rarefaction](#skipalphararefaction)
    * [--skip_taxonomy](#skiptaxonomy)
    * [--skip_barplot](#skipbarplot)
    * [--skip_abundance_tables](#skipabundancetables)
    * [--skip_diversity_indices](#skipdiversityindices)
    * [--skip_ancom](#skipancom)
  * [Job Resources](#job-resources)
    * [Automatic resubmission](#automatic-resubmission)
    * [Custom resource requests](#custom-resource-requests)
  * [AWS Batch specific parameters](#aws-batch-specific-parameters)
    * [--awsqueue](#awsqueue)
    * [--awsregion](#awsregion)
  * [Other command line parameters](#other-command-line-parameters)
    * [--outdir](#outdir)
    * [--email](#email)
    * [-name](#name)
    * [-resume](#resume)
    * [-c](#c)
    * [--custom_config_version](#customconfigversion)
    * [--custom_config_base](#customconfigbase)
    * [--max_memory](#maxmemory)
    * [--max_time](#maxtime)
    * [--max_cpus](#maxcpus)
    * [--plaintext_email](#plaintextemail)
    * [--monochrome_logs](#monochromelogs)
    * [--multiqc_config](#multiqcconfig)

## General Nextflow info

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/ampliseq \
    -profile singularity \
    --reads "data" \
    --FW_primer GTGYCAGCMGCCGCGGTAA \
    --RV_primer GGACTACNVGGGTWTCTAAT \
    --metadata "data/Metadata.tsv"
```

This will launch the pipeline with the `singularity` configuration profile. See below [`-profile`](#-profile) for more information about profiles.

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

## Main arguments

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. Note that multiple profiles can be loaded, for example: `-profile docker` - the order of arguments is important!

If `-profile` is not specified at all the pipeline will be run locally and expects all software to be installed and available on the `PATH`.

* `awsbatch`
  * A generic configuration profile to be used with AWS Batch.
* `conda`
  * A generic configuration profile to be used with [conda](https://conda.io/docs/)
  * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `docker`
  * A generic configuration profile to be used with [Docker](http://docker.com/)
  * Pulls software from dockerhub: [`nfcore/ampliseq`](http://hub.docker.com/r/nfcore/ampliseq/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
  * Pulls software from DockerHub: [`nfcore/ampliseq`](http://hub.docker.com/r/nfcore/ampliseq/)
* `test`
  * A profile with a complete configuration for automated testing
  * Includes links to test data so needs no other parameters

### `--reads`

Use this to specify the location of your input paired-end FastQ files.  

For example:

```bash
--reads 'path/to/data'
```

Example for input data organization from one sequencing run with two samples:

```bash
data
  |-sample1_1_L001_R1_001.fastq.gz
  |-sample1_1_L001_R2_001.fastq.gz
  |-sample2_1_L001_R1_001.fastq.gz
  |-sample2_1_L001_R2_001.fastq.gz
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The folder must contain gzip compressed paired-end demultiplexed fastq files. If the file names do not follow the default (`"/*_R{1,2}_001.fastq.gz"`), please check [`--extension`](#--extension).
3. If your data is scattered, a directory with symlinks to your actual data might be a solution.
4. All sequencing data should originate from one sequencing run, because processing relies on run-specific error models that are unreliable when data from several sequencing runs are mixed. Sequencing data originating from multiple sequencing runs requires additionally the parameter `--multipleSequencingRuns` and a specific folder structure, see [here](#--multipleSequencingRuns).

### `--FW_primer` and `--RV_primer`

In amplicon sequencing methods, PCR with specific primers produces the amplicon of intrest. These primer sequences need to be trimmed from the reads before further processing and are also required for producing an appropriate classifier. For example:

```bash
--FW_primer GTGYCAGCMGCCGCGGTAA --RV_primer GGACTACNVGGGTWTCTAAT
```

### `--metadata`

This is optional, but for performing downstream analysis such as barplots, diversity indices or differential abundance testing, a metadata file is essential. For example:

```bash
--metadata "path/to/metadata.tsv"
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The metadata file has to follow the [QIIME2 specifications](https://docs.qiime2.org/2019.10/tutorials/metadata/)
3. In case of multiple sequencing runs, specific naming of samples are required, see [here](#--multipleSequencingRuns)

## Other input options

### `--extension`

Indicates the naming of sequencing files (default: `"/*_R{1,2}_001.fastq.gz"`).

Please note:

1. The prepended slash (`/`) is required
2. The star (`*`) is the required wildcard for sample names
3. The curly brackets (`{}`) enclose the orientation for paired end reads, seperated by a comma (`,`).
4. The pattern must be enclosed in quotes

For example for one sample (name: `1`) with forward (file: `1_a.fastq.gz`) and reverse (file: `1_b.fastq.gz`) reads in folder `data`:

```bash
--reads "data" --extension "/*_{a,b}.fastq.gz"
```

### `--multipleSequencingRuns`

If samples were sequenced in multiple sequencing runs. Expects one subfolder per sequencing run
in the folder specified by [`--reads`](#--reads) containing sequencing data of the specific run. Also, fastQC and MultiQC are skipped because multiple sequencing runs might create overlapping file names that crash MultiQC.

To prevent overlapping sample names from multiple sequencing runs, sample names obtained from the sequencing files will be renamed automatically by adding the folder name as prefix seperated by a string specified by [`--split`](#--split). Accordingly, the sample name column in the metadata file [`--metadata`](#--metadata) require values following `subfolder-samplename`.

Example for input data organization:

```bash
data
  |-run1
  |  |-sample1_1_L001_R{1,2}_001.fastq.gz
  |  |-sample2_1_L001_R{1,2}_001.fastq.gz
  |
  |-run2
     |-sample3_1_L001_R{1,2}_001.fastq.gz
     |-sample4_1_L001_R{1,2}_001.fastq.gz
```

In this example the first column in the metadata file requires the values `run1-sample1` ... `run2-sample4` (instead of `sample1`, ..., `sample4`).

Example command to analyze this data in one pipeline run:

```bash
nextflow run nf-core/ampliseq \
    -profile singularity \
    --reads "data" \
    --FW_primer GTGYCAGCMGCCGCGGTAA \
    --RV_primer GGACTACNVGGGTWTCTAAT \
    --metadata "data/Metadata.tsv" \
    --multipleSequencingRuns
```

### `--split`

A string that will be used between the prepended run/folder name and the sample name. Only used with [`--multipleSequencingRuns`](#--multipleSequencingRuns) (default: `"-"`).

For example using the string `link`:

```bash
--split "link"
```

Please note:

1. Run/folder names may not contain the string specified by `--split`
2. No underscore(s) allowed
3. Must be enclosed in quotes
4. The metadata sheet has to be adjusted, instead of using `run-sample` in the first column, in this example `runlinksample` is required

### `--phred64`

If the sequencing data has PHRED 64 encoded quality scores (default: PHRED 33)

## Cutoffs

### `--trunclenf` and `--trunclenr`

Read denoising by DADA2 creates an error profile specific to a sequencing run and uses this to correct sequencing errors. This method requires all reads to have the same length and as high quality as possible while maintaining at least 20 bp overlap for merging. One cutoff for the forward read `--trunclenf` and one for the reverse read `--trunclenr` truncate all longer reads at that position and drop all shorter reads.
These cutoffs are usually chosen visually using [`--untilQ2import`](#--untilQ2import), inspecting the quality plots in "results/demux", and resuming analysis with [`--Q2imported`](#--Q2imported). If not set, these cutoffs will be determined automatically for the position before the mean quality score drops below [`--trunc_qmin`](#--trunc_qmin).

For example:

```bash
--trunclenf 180 --trunclenr 120
```

Please note:

1. Too agressive truncation might lead to insufficient overlap for read merging
2. Too less truncation might reduce denoised reads
3. The code choosing these values automatically cannot take the points above into account, therefore setting `--trunclenf` and `--trunclenr` is recommended

### `--trunc_qmin`

Automatically determine [`--trunclenf` and `--trunclenr`](#--trunclenf-and---trunclenr) before the mean quality score drops below `--trunc_qmin` (default: 25). For example:

```bash
--trunc_qmin 35
```

Please note:

1. The code choosing `--trunclenf` and `--trunclenr` using `--trunc_qmin` automatically cannot take amplicon length or overlap requirements for merging into account, therefore setting `--trunclenf` and `--trunclenr` is preferred
2. The default value of 25 is recommended. However, very good quality data with large paired sequence overlap might justify a higher value (e.g. 35). Also, very low quality data might require a lower value.

## Other options

### `--untilQ2import`

Computes all steps until quality plots aiding the choosing of [`--trunclenf` and `--trunclenr`](#--trunclenf-and---trunclenr).

### `--Q2imported`

Analysis starting with a QIIME2 artefact with trimmed reads, typically produced before with [`--untilQ2import`](#--untilQ2import). This is only supported for data from a single sequencing run.

For data from multiple sequencing runs with [`--multipleSequencingRuns`](#--multipleSequencingRuns) the pipeline can be first run with `--untilQ2import` and next run without `--untilQ2import` but with [`-resume`](#-resume-single-dash). For more details see [here](#visually-choosing-sequencing-read-truncation-cutoffs).

### `--keepIntermediates`

Keep additional intermediate files, such as trimmed reads or various QIIME2 archives.

#### Visually choosing sequencing read truncation cutoffs

While [`--untilQ2import`](#--untilQ2import) with `--multipleSequencingRuns` is currently supported, [`--Q2imported`](#--Q2imported) is not. The pipeline can be first run with `--untilQ2import`, than [`--trunclenf` and `--trunclenr`](#--trunclenf-and---trunclenr) are visually chosen, and than the pipeline can be continued without `--untilQ2import` but with `--trunlenf`, `--trunlenr`, and [`-resume`](#-resume).

For example:

(1) To produce quality plots and choose truncation values:

```bash
nextflow run nf-core/ampliseq \
    -profile singularity \
    --reads "data" \
    --FW_primer GTGYCAGCMGCCGCGGTAA \
    --RV_primer GGACTACNVGGGTWTCTAAT \
    --metadata "data/Metadata.tsv" \
    --multipleSequencingRuns \
    --untilQ2import
```

(2) To finish analysis:

```bash
nextflow run nf-core/ampliseq \
    -profile singularity \
    --reads "data" \
    --FW_primer GTGYCAGCMGCCGCGGTAA \
    --RV_primer GGACTACNVGGGTWTCTAAT \
    --metadata "data/Metadata.tsv" \
    --multipleSequencingRuns \
    --trunclenf 200 \
    --trunclenr 180 \
    -resume
```

## Reference database

By default, the workflow downloads [SILVA](https://www.arb-silva.de/) [v132](https://www.arb-silva.de/documentation/release-132/) and extracts reference sequences and taxonomy clustered at 99% similarity and trains a Naive Bayes classifier to assign taxonomy to features.

### `--classifier`

If you have trained a compatible classifier before, from sources such as [SILVA](https://www.arb-silva.de/), [Greengenes](http://greengenes.secondgenome.com/downloads) or [RDP](https://rdp.cme.msu.edu/). For example:

```bash
--classifier "FW_primer-RV_primer-classifier.qza"
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The cassifier is a Naive Bayes classifier produced by "qiime feature-classifier fit-classifier-naive-bayes" (e.g. by this pipeline or from [QIIME2 resources](https://docs.qiime2.org/2019.10/data-resources/))
3. The primer pair for the amplicon PCR and the computing of the classifier are exactly the same
4. The classifier has to be trained by the same version of scikit-learn as this version of the pipeline uses (0.21.2)

### `--classifier_removeHash`

Remove all hash signs from taxonomy strings, resolves a rare ValueError during classification (process classifier).

## Statistics

### `--metadata_category`

Here columns in the metadata sheet can be chosen with groupings that are used for diversity indices and differential abundance analysis. By default, all suitable columns in the metadata sheet will be used if this option is not specified. Suitable are columns which are categorical (not numerical) and have multiple different values which are not all unique. For example:

```bash
--metadata_category "treatment1,treatment2"
```

Please note the following requirements:

1. Comma seperated list enclosed in quotes
2. May not contain whitespace characters
3. Each comma seperated term has to match exactly one column name in the metadata sheet

## Filters

### `--retain_untrimmed`

When read sequences are trimmed, untrimmed read pairs are discarded routinely. Use this option to retain untrimmed read pairs. This is usually not recommended and is only of advantage for specific protocols that prevent sequencing PCR primers. For example:

```bash
--retain_untrimmed
```

### `--exclude_taxa`

Depending on the primers used, PCR might amplify unwanted or off-target DNA. By default sequences originating from mitochondria or chloroplasts are removed.
The taxa specified are excluded from further analysis.

For example to exclude any taxa that contain mitochondria, chloroplast, or archea:

```bash
--exclude_taxa "mitochondria,chloroplast,archea"
```

If you prefer not filtering the data, specify:

```bash
--exclude_taxa "none"
```

Please note the following requirements:

1. Comma seperated list enclosed in quotes
2. May not contain whitespace characters
3. Features that contain one or several of these terms in their taxonomical classification are excluded from further analysis
4. The taxonomy level is not taken into consideration

### `--min_frequency`

Remove entries from the feature table below an absolute abundance threshold (default: 1, meaning filter is disabled). Singletons are often regarded as artifacts, choosing a value of 2 removes sequences with less than 2 total counts from the feature table.

For example to remove singletons choose:

```bash
--min_frequency 2
```

### `--min_samples`

Filtering low prevalent features from the feature table, e.g. keeping only features that are present in at least two samples can be achived by choosing a value of 2 (default: 1, meaning filter is disabled). Typically only used when having replicates for all samples.

For example to retain features that are present in at least two sample:

```bash
--min_samples 2
```

Please note this is independent of abundance.

## Skipping steps

### `--onlyDenoising`

Skip all steps after denoising, produce only sequences and abundance tables on ASV level.

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

Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files hosted at [`nf-core/configs`](https://github.com/nf-core/configs/tree/master/conf) for examples.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-core-invite.herokuapp.com/).

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

Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.

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

Note - you can use this to override pipeline defaults.

### `--custom_config_version`

Provide git commit id for custom Institutional configs hosted at `nf-core/configs`. This was implemented for reproducibility purposes. Default is set to `master`.

```bash
## Download and use config file with following git commid id
--custom_config_version d52db660777c4bf36546ddb188ec530c3ada1b96
```

### `--custom_config_base`

If you're running offline, nextflow will not be able to fetch the institutional config files
from the internet. If you don't need them, then this is not a problem. If you do need them,
you should download the files from the repo and tell nextflow where to find them with the
`custom_config_base` option. For example:

```bash
## Download and unzip the config files
cd /path/to/my/configs
wget https://github.com/nf-core/configs/archive/master.zip
unzip master.zip

## Run the pipeline
cd /path/to/my/data
nextflow run /path/to/pipeline/ --custom_config_base /path/to/my/configs/configs-master/
```

> Note that the nf-core/tools helper package has a `download` command to download all required pipeline
> files + singularity containers + institutional configs in one go for you, to make this process easier.

### `--max_memory`

Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'`

### `--max_time`

Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`

Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--plaintext_email`

Set to receive plain-text e-mails instead of HTML formatted.

### `--monochrome_logs`

Set to disable colourful command line output and live life in monochrome.

### `--multiqc_config`

Specify a path to a custom MultiQC configuration file.
