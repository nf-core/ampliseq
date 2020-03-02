# ![nf-core/ampliseq](docs/images/nf-core-ampliseq_logo.png)

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)
[![nf-core](https://img.shields.io/badge/nf--core-pipeline-brightgreen.svg)](https://nf-co.re/)
[![DOI](https://zenodo.org/badge/150448201.svg)](https://zenodo.org/badge/latestdoi/150448201)
[![Cite Preprint](https://img.shields.io/badge/Cite%20Us!-Cite%20Preprint-orange)](https://biorxiv.org/cgi/content/short/2019.12.17.880468v1)

[![GitHub Actions CI Status](https://github.com/nf-core/ampliseq/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/ampliseq/actions)
[![GitHub Actions Linting Status](https://github.com/nf-core/ampliseq/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/ampliseq/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/ampliseq.svg)](https://hub.docker.com/r/nfcore/ampliseq)
![Singularity Container available](https://img.shields.io/badge/singularity-available-7E4C74.svg)

[![Joins us on Slack](https://img.shields.io/badge/slack-nfcore/ampliseq-blue.svg)](https://nfcore.slack.com/channels/ampliseq)

## Introduction

**nfcore/ampliseq** is a bioinformatics analysis pipeline used for 16S rRNA amplicon sequencing data.

The workflow processes raw data from FastQ inputs ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)), trims primer sequences from the reads ([Cutadapt](https://journal.embnet.org/index.php/embnetjournal/article/view/200)), imports data into [QIIME2](https://www.nature.com/articles/s41587-019-0209-9), generates amplicon sequencing variants (ASV, [DADA2](https://www.nature.com/articles/nmeth.3869)), classifies features against the [SILVA](https://www.arb-silva.de/) [v132](https://www.arb-silva.de/documentation/release-132/) database, excludes unwanted taxa, produces absolute and relative feature/taxa count tables and plots, plots alpha rarefaction curves, computes alpha and beta diversity indices and plots thereof, and finally calls differentially abundant taxa ([ANCOM](https://www.ncbi.nlm.nih.gov/pubmed/26028277)). See the [output documentation](docs/output.md) for more details of the results.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Documentation

The nf-core/ampliseq pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](https://nf-co.re/usage/installation)
2. Pipeline configuration
    * [Local installation](https://nf-co.re/usage/installation)
    * [Adding your own system config](https://nf-co.re/usage/adding_own_config)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](https://nf-co.re/usage/troubleshooting)

## Quick Start

i. Install [`nextflow`](https://nf-co.re/usage/installation)

ii. Install either [`Docker`](https://docs.docker.com/engine/installation/) or [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) for full pipeline reproducibility (please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))

iii. Download the pipeline and test it on a minimal dataset with a single command

```bash
nextflow run nf-core/ampliseq -profile test,<docker/singularity/conda/institute>
```

> Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.

iv. Start running your own analysis!

```bash
nextflow run nf-core/ampliseq -profile <docker/singularity/conda/institute> --reads 'data' --FW_primer <forward-primer-sequence> --RV_primer <reverse-primer-sequence> --metadata "data/Metadata.tsv"
```

See [usage docs](docs/usage.md) for all of the available options when running the pipeline.

## Credits

nf-core/ampliseq was originally written by Daniel Straub ([@d4straub](https://github.com/d4straub)) and Alexander Peltzer ([@apeltzer](https://github.com/apeltzer)) for use at the [Quantitative Biology Center (QBiC)](http://www.qbic.life) and [Microbial Ecology, Center for Applied Geosciences](http://www.uni-tuebingen.de/de/104325), part of Eberhard Karls Universität Tübingen (Germany).

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on [Slack](https://nfcore.slack.com/channels/ampliseq) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citation

If you use  nf-core/ampliseq for your analysis, please cite it using the following doi:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3568091.svg)](https://doi.org/10.5281/zenodo.3568091)

The pre-print can be cited as follows:

[![DOI](https://img.shields.io/badge/-Cite%20Preprint-orange)](https://www.biorxiv.org/content/10.1101/2019.12.17.880468v1)

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).  
> ReadCube: [Full Access Link](https://rdcu.be/b1GjZ)
