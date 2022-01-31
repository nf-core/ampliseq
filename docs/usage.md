# nf-core/ampliseq: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/ampliseq/usage](https://nf-co.re/ampliseq/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Input specifications

The input data can be passed to nf-core/ampliseq in three possible ways using the `--input` parameter, either a folder containing zipped FastQ files, a tab-seperated samplesheet, or a fasta file to be taxonomically classified.

Optionally, a metadata sheet can be specified for downstream analysis.

### Direct FASTQ input

The easiest way is to specify directly the path to the folder that contains your input FASTQ files. For example:

```console
--input 'path/to/data/'
```

File names must follow a specific pattern, default is `/*_R{1,2}_001.fastq.gz`, but this can be adjusted with `--extension`.

For example, the following files in folder `data` would be processed as `sample1` and `sample2`:

```console
data
    |-sample1_1_L001_R1_001.fastq.gz
    |-sample1_1_L001_R2_001.fastq.gz
    |-sample2_1_L001_R1_001.fastq.gz
    |-sample2_1_L001_R2_001.fastq.gz
```

All sequencing data should originate from one sequencing run, because processing relies on run-specific error models that are unreliable when data from several sequencing runs are mixed. Sequencing data originating from multiple sequencing runs requires additionally the parameter `--multiple_sequencing_runs` and a specific folder structure, for example:

```console
data
    |-runA
    |   |-sample1_1_L001_R1_001.fastq.gz
    |   |-sample1_1_L001_R2_001.fastq.gz
    |   |-sample2_1_L001_R1_001.fastq.gz
    |   |-sample2_1_L001_R2_001.fastq.gz
    |
    |-runB
        |-sample3_1_L001_R1_001.fastq.gz
        |-sample3_1_L001_R2_001.fastq.gz
        |-sample4_1_L001_R1_001.fastq.gz
        |-sample4_1_L001_R2_001.fastq.gz
```

Where `sample1` and `sample2` were sequenced in one sequencing run and `sample3` and `sample4` in another sequencing run.

Please note the following additional requirements:

* Files names must be unique
* Valid file extensions: `.fastq.gz`, `.fq.gz` (files must be compressed)
* The path must be enclosed in quotes
* `--extension` must have at least one `*` wildcard character
* When using the pipeline with paired end data, the `--extension` must use `{1,2}` (or similar) notation to specify read pairs
* To run single-end data you must additionally specify `--single_end` and `--extension` may not include curly brackets `{}`
* Sample identifiers are extracted from file names, i.e. the string before the first underscore `_`, these must be unique (also across sequencing runs)
* If your data is scattered, produce a sample sheet

### Samplesheet input

The sample sheet file is an alternative way to provide input reads, it must be a tab-separated file ending with `.tsv` that must have two to four columns with the following headers:

| Column | Necessity | Description |
|-|-|-|
| sampleID | required | Unique sample identifiers |
| forwardReads | required | Paths to (forward) reads zipped FastQ files |
| reverseReads | optional | Paths to reverse reads zipped FastQ files, required if the data is paired-end |
| run | optional | If the data was produced by multiple sequencing runs, any string |

```console
--input 'path/to/samplesheet.tsv'
```

For example, the samplesheet may contain:

```console
sampleID    forwardReads    reverseReads    run
sample1    ./data/S1_R1_001.fastq.gz    ./data/S1_R2_001.fastq.gz    A
sample2    ./data/S2_fw.fastq.gz    ./data/S2_rv.fastq.gz    A
sample3    ./S4x.fastq.gz    ./S4y.fastq.gz    B
sample4    ./a.fastq.gz    ./b.fastq.gz    B
```

Please note the following requirements:

* 2 to 4 tab-seperated columns
* Valid file extension: `.tsv`
* Must contain the header `sampleID` and `forwardReads`
* May contain the header `reverseReads` and `run`
* Sample IDs must be unique
* Sample IDs must not contain a dot `.`
* Sample IDs starting with a number are not allowed when using metadata (because these strings will be modified)
* FastQ files must be compressed (`.fastq.gz`, `.fq.gz`)
* Within one samplesheet, only one type of raw data should be specified (same amplicon & sequencing method)

An [example samplesheet](../assets/samplesheet.tsv) has been provided with the pipeline.

### Fasta file input, to taxonomically classify previously produced ASV/OTU sequences

When pointing at a file ending with `.fasta`, `.fna` or `.fa`, the containing sequences will be taxonomically classified. All other pipeline steps will be skipped.

```console
--input 'path/to/amplicon_sequences.fasta'
```

Please note the following requirements:

* Valid file extensions: `.fasta`, `.fna` or `.fa`

### Metadata

The metadata file must be tab-separated with a header line, for more details please see the [nf-core/ampliseq website documentation](https://nf-co.re/ampliseq/parameters).

```console
--metadata "path/to/metadata.tsv"
```

For example:

```console
ID    condition
sample1    control
sample2    treatment
sample3    control
sample4    treatment
```

## Running the pipeline

The typical command for running the pipeline is as follows:

```console
nextflow run nf-core/ampliseq \
    -profile singularity \
    --input "data" \
    --FW_primer GTGYCAGCMGCCGCGGTAA \
    --RV_primer GGACTACNVGGGTWTCTAAT \
    --metadata "data/Metadata.tsv"
```

This will launch the pipeline with the `singularity` configuration profile. See below [`-profile`](#-profile) for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```console
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

> **NB:** If the data originates from multiple sequencing runs, the error profile of each of those sequencing runs needs to be considered separately. Using the `run` column in the samplesheet input or adding `--multiple_sequencing_runs` for Direct FASTQ input will separate certain processes by the sequencing run. Please see the following example:

<p align="center">
    <img src="images/ampliseq_workflow_multiplesequencingruns.png" alt="nf-core/ampliseq workflow overview with --multiple_sequencing_runs" width="40%">
</p>

See the [nf-core/ampliseq website documentation](https://nf-co.re/ampliseq/parameters) for more information about pipeline specific parameters.

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```console
nextflow pull nf-core/ampliseq
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/ampliseq releases page](https://github.com/nf-core/ampliseq/releases) and find the latest version number - numeric only (eg. `2.0.0`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 2.0.0`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below. When using Biocontainers, most of these software packaging methods pull Docker containers from quay.io e.g [FastQC](https://quay.io/repository/biocontainers/fastqc) except for Singularity which directly downloads Singularity images via https hosted by the [Galaxy project](https://depot.galaxyproject.org/singularity/) and Conda which downloads and installs software locally from [Bioconda](https://bioconda.github.io/).

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

* `docker`
    * A generic configuration profile to be used with [Docker](https://docker.com/)
* `singularity`
    * A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
* `podman`
    * A generic configuration profile to be used with [Podman](https://podman.io/)
* `shifter`
    * A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
* `charliecloud`
    * A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
* `conda`
    * A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.
* `test`, `test_multi`, `test_full`, `test_pacbio_its`, `test_iontorrent`, `test_doubleprimers`
    * Profiles with a complete configuration for automated testing
    * Includes links to test data so needs no other parameters

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

For example, if the nf-core/rnaseq pipeline is failing after multiple re-submissions of the `STAR_ALIGN` process due to an exit code of `137` this would indicate that there is an out of memory issue:

```console
[62/149eb0] NOTE: Process `RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137) -- Execution is retried (1)
Error executing process > 'RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)'

Caused by:
    Process `RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137)

Command executed:
    STAR \
        --genomeDir star \
        --readFilesIn WT_REP1_trimmed.fq.gz  \
        --runThreadN 2 \
        --outFileNamePrefix WT_REP1. \
        <TRUNCATED>

Command exit status:
    137

Command output:
    (empty)

Command error:
    .command.sh: line 9:  30 Killed    STAR --genomeDir star --readFilesIn WT_REP1_trimmed.fq.gz --runThreadN 2 --outFileNamePrefix WT_REP1. <TRUNCATED>
Work dir:
    /home/pipelinetest/work/9d/172ca5881234073e8d76f2a19c88fb

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`
```

To bypass this error you would need to find exactly which resources are set by the `STAR_ALIGN` process. The quickest way is to search for `process STAR_ALIGN` in the [nf-core/rnaseq Github repo](https://github.com/nf-core/rnaseq/search?q=process+STAR_ALIGN). We have standardised the structure of Nextflow DSL2 pipelines such that all module files will be present in the `modules/` directory and so based on the search results the file we want is `modules/nf-core/software/star/align/main.nf`. If you click on the link to that file you will notice that there is a `label` directive at the top of the module that is set to [`label process_high`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/modules/nf-core/software/star/align/main.nf#L9). The [Nextflow `label`](https://www.nextflow.io/docs/latest/process.html#label) directive allows us to organise workflow processes in separate groups which can be referenced in a configuration file to select and configure subset of processes having similar computing requirements. The default values for the `process_high` label are set in the pipeline's [`base.config`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L33-L37) which in this case is defined as 72GB. Providing you haven't set any other standard nf-core parameters to __cap__ the [maximum resources](https://nf-co.re/usage/configuration#max-resources) used by the pipeline then we can try and bypass the `STAR_ALIGN` process failure by creating a custom config file that sets at least 72GB of memory, in this case increased to 100GB. The custom config below can then be provided to the pipeline via the [`-c`](#-c) parameter as highlighted in previous sections.

```nextflow
process {
    withName: STAR_ALIGN {
        memory = 100.GB
    }
}
```

> **NB:** We specify just the process name i.e. `STAR_ALIGN` in the config file and not the full task name string that is printed to screen in the error message or on the terminal whilst the pipeline is running i.e. `RNASEQ:ALIGN_STAR:STAR_ALIGN`. You may get a warning suggesting that the process selector isn't recognised but you can ignore that if the process name has been specified correctly. This is something that needs to be fixed upstream in core Nextflow.

### Updating containers

The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. If for some reason you need to use a different version of a particular tool with the pipeline then you just need to identify the `process` name and override the Nextflow `container` definition for that process using the `withName` declaration. For example, in the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline a tool called [Pangolin](https://github.com/cov-lineages/pangolin) has been used during the COVID-19 pandemic to assign lineages to SARS-CoV-2 genome sequenced samples. Given that the lineage assignments change quite frequently it doesn't make sense to re-release the nf-core/viralrecon everytime a new version of Pangolin has been released. However, you can override the default container used by the pipeline by creating a custom config file and passing it as a command-line argument via `-c custom.config`.

1. Check the default version used by the pipeline in the module file for [Pangolin](https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/modules/nf-core/software/pangolin/main.nf#L14-L19)
2. Find the latest version of the Biocontainer available on [Quay.io](https://quay.io/repository/biocontainers/pangolin?tag=latest&tab=tags)
3. Create the custom config accordingly:

    * For Docker:

        ```nextflow
        process {
            withName: PANGOLIN {
                container = 'quay.io/biocontainers/pangolin:3.0.5--pyhdfd78af_0'
            }
        }
        ```

    * For Singularity:

        ```nextflow
        process {
            withName: PANGOLIN {
                container = 'https://depot.galaxyproject.org/singularity/pangolin:3.0.5--pyhdfd78af_0'
            }
        }
        ```

    * For Conda:

        ```nextflow
        process {
            withName: PANGOLIN {
                conda = 'bioconda::pangolin=3.0.5'
            }
        }
        ```

> **NB:** If you wish to periodically update individual tool-specific results (e.g. Pangolin) generated by the pipeline then you must ensure to keep the `work/` directory otherwise the `-resume` ability of the pipeline will be compromised and it will restart from scratch.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```console
NXF_OPTS='-Xms1g -Xmx4g'
```
