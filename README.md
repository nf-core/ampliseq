# ![nf-core/ampliseq](docs/images/nf-core-ampliseq_logo.png)

**16S rRNA amplicon sequencing analysis workflow using QIIME2**.

[![nf-core](https://img.shields.io/badge/nf--core-pipeline-brightgreen.svg)](https://nf-co.re/)
[![DOI](https://zenodo.org/badge/150448201.svg)](https://zenodo.org/badge/latestdoi/150448201)
[![Cite Preprint](https://img.shields.io/badge/Cite%20Us!-Cite%20Publication-important)](https://doi.org/10.3389/fmicb.2020.550420)

[![GitHub Actions CI Status](https://github.com/nf-core/ampliseq/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/ampliseq/actions)
[![GitHub Actions Linting Status](https://github.com/nf-core/ampliseq/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/ampliseq/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/ampliseq.svg)](https://hub.docker.com/r/nfcore/ampliseq)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23ampliseq-4A154B?logo=slack)](https://nfcore.slack.com/channels/ampliseq)

## Introduction

**nfcore/ampliseq** is a bioinformatics analysis pipeline used for 16S rRNA or ITS amplicon sequencing data (currently supported is Illumina paired end or PacBio).

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Quick Start

1. Install [`nextflow`](https://nf-co.re/usage/installation)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`Podman`](https://podman.io/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```bash
    nextflow run nf-core/ampliseq -profile test,<docker/singularity/podman/conda/institute>
    ```

    > Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.

4. Start running your own analysis!

    ```bash
    nextflow run nf-core/ampliseq -profile <docker/singularity/podman/conda/institute> --input "data" --FW_primer GTGYCAGCMGCCGCGGTAA --RV_primer GGACTACNVGGGTWTCTAAT --metadata "data/Metadata.tsv"
    ```

See [usage docs](https://nf-co.re/ampliseq/usage) and [parameter docs](https://nf-co.re/ampliseq/parameters) for all of the available options when running the pipeline.

## Pipeline Summary

By default, the pipeline currently performs the following:

* Sequencing quality control ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
* Trimming of reads ([Cutadapt](https://journal.embnet.org/index.php/embnetjournal/article/view/200))
* Illumina read processing with [QIIME2](https://www.nature.com/articles/s41587-019-0209-9)
* Infer Amplicon Sequence Variants (ASVs) ([DADA2](https://doi.org/10.1038/nmeth.3869))
* Taxonomical classification based on [SILVA](https://www.arb-silva.de/) [v132](https://www.arb-silva.de/documentation/release-132/) or [UNITE](https://unite.ut.ee/) database
* excludes unwanted taxa, produces absolute and relative feature/taxa count tables and plots, plots alpha rarefaction curves, computes alpha and beta diversity indices and plots thereof ([QIIME2](https://www.nature.com/articles/s41587-019-0209-9))
* Calls differentially abundant taxa ([ANCOM](https://www.ncbi.nlm.nih.gov/pubmed/26028277))
* Overall pipeline run summaries ([MultiQC](https://multiqc.info/))

## Documentation

The nf-core/ampliseq pipeline comes with documentation about the pipeline: [usage](https://nf-co.re/ampliseq/usage) and [output](https://nf-co.re/ampliseq/output).

## Credits

nf-core/ampliseq was originally written by Daniel Straub ([@d4straub](https://github.com/d4straub)) and Alexander Peltzer ([@apeltzer](https://github.com/apeltzer)) for use at the [Quantitative Biology Center (QBiC)](http://www.qbic.life) and [Microbial Ecology, Center for Applied Geosciences](http://www.uni-tuebingen.de/de/104325), part of Eberhard Karls Universität Tübingen (Germany).

We thank the following people for their extensive assistance in the development of this pipeline (in alphabetical order):

* [Daniel Lundin](https://github.com/erikrikarddaniel)
* [Diego Brambilla](https://github.com/DiegoBrambilla)
* [Emelie Nilsson](https://github.com/emnilsson)
* [Jeanette Tångrot](https://github.com/jtangrot)
* [Sabrina Krakau](https://github.com/skrakau)

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#ampliseq` channel](https://nfcore.slack.com/channels/ampliseq) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

If you use `nf-core/ampliseq` for your analysis, please cite the `ampliseq` article as follows:
> Daniel Straub, Nia Blackwell, Adrian Langarica-Fuentes, Alexander Peltzer, Sven Nahnsen, Sara Kleindienst **Interpretations of Environmental Microbial Community Studies Are Biased by the Selected 16S rRNA (Gene) Amplicon Sequencing Pipeline** *Frontiers in Microbiology* 2020, 11:2652 [doi: 10.3389/fmicb.2020.550420](https://doi.org/10.3389/fmicb.2020.550420).

You can cite the `nf-core/ampliseq` zenodo record for a specific version using the following [doi: 10.5281/zenodo.1493841](https://zenodo.org/badge/latestdoi/150448201)

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
> ReadCube: [Full Access Link](https://rdcu.be/b1GjZ)

<!-- TODO In addition, references of tools and data used in this pipeline are as follows: -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->
