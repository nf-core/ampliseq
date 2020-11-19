# nf-core/ampliseq: Output

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/ampliseq/output](https://nf-co.re/ampliseq/output)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [nf-core/ampliseq: Output](#nf-coreampliseq-output)
  * [Pipeline overview](#pipeline-overview)
    * [FastQC](#fastqc) - Read quality control
    * [Cutadapt](#cutadapt) - Primer trimming
    * [MultiQC](#multiqc) - Aggregate report describing results
    * [QIIME2](#qiime2) - Import & quality control
    * [DADA2](#dada2) - Infer Amplicon Sequence Variants (ASVs)
    * [Taxonomic classification](#taxonomic-classification) - Taxonomical classification of ASVs
    * [Exclude taxa](#exclude-taxa) - Remove unwanted ASV based on taxonomy
    * [Relative abundance tables](#relative-abundance-tables) - Exported relative abundance tables
    * [Barplot](#barplot) - Interactive barplot
    * [Alpha diversity rarefaction curves](#alpha-diversity-rarefaction-curves) - Rarefaction curves for quality control
    * [Alpha diversity indices](#alpha-diversity-indices) - Diversity within samples
    * [Beta diversity indices](#beta-diversity-indices) - Diversity between samples (e.g. PCoA plots)
    * [ANCOM](#ancom) - Differential abundance analysis
    * [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution
  * [More help](#more-help)
  * [Citations](#citations)

### FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences.

For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

**Output files:**

* `fastqc/`
  * `*_fastqc.html`: FastQC report containing quality metrics for your untrimmed raw fastq files.
* `fastqc/zips/`
  * `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality.

### Cutadapt

[Cutadapt](https://journal.embnet.org/index.php/embnetjournal/article/view/200) is trimming primer sequences from sequencing reads. Primer sequences are non-biological sequences that often introduce point mutations that do not reflect sample sequences. This is especially true for degenerated PCR primer. If primer trimming would be omitted, artifactual amplicon sequence variants might be computed by the denoising tool or sequences might be lost due to become labelled as PCR chimera.

**Output files:**

* `trimmed/logs/`: directory containing log files with retained reads, trimming percentage, etc. for each sample.

### MultiQC

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarizing all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability.

For more information about how to use MultiQC reports, see [https://multiqc.info](https://multiqc.info).

**Output files:**

* `multiqc/`
  * `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  * `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  * `multiqc_plots/`: directory containing static images from the report in various formats.

### QIIME2

**Quantitative Insights Into Microbial Ecology 2** ([QIIME2](https://qiime2.org/)) is a next-generation microbiome bioinformatics platform and the successor of the widely used [QIIME1](https://www.nature.com/articles/nmeth.f.303). QIIME2 is currently **under heavy development** and often updated, this version of ampliseq uses QIIME2 2019.10. QIIME2 has a wide variety of analysis tools available and has excellent support in its [forum](https://docs.qiime2.org/2019.10/).

At this point of the analysis the trimmed reads are imported into QIIME2 and an interactive quality plot is made.

**Output files:**

* `demux/`  
  * `index.html`: Quality plots that can be viewed in your web browser.
  * `demux.qza` (only when --untilQ2import is true): QIIME2 artefact for imported reads.

All following analysis steps are performed in QIIME2, except DADA2 in the case of pacbio data.

### DADA2

[DADA2](https://www.nature.com/articles/nmeth.3869) performs fast and accurate sample inference from amplicon data with single-nucleotide resolution. It infers exact amplicon sequence variants (ASVs) from amplicon data with fewer false positives than many other methods while maintaining high sensitivity.

DADA2 computes an error model on the sequencing reads (forward and reverse independently), therefore quality filtering or paired read merging may not be performed before. Each sequencing run varies in their error profile and it is recommended that DADA2 runs separately on data from each run individually. It is recommended to use the ampliseq option `--multipleSequencingRuns` to analyse such data.

DADA2 reduces sequence errors and dereplicates sequences by quality filtering, denoising, read pair merging (for paired end Illumina reads only) and PCR chimera removal.

**Output files:**

* `representative_sequences/unfiltered/`  
  * `sequences.fasta`: Fasta file with ASV sequences.
  * `index.html`: ASV IDs, sequences and blast results in an interactive table that can be viewed in your web browser.
  * `rep-seqs.qza`: QIIME2 data artefact.
* `abundance-table/unfiltered/`  
  * `dada_report.txt`: DADA2 verbose output.
  * `dada_stats.tsv`: Tab-separated table of DADA2 statistics.
  * `feature-table.biom`: Abundance table in biom format for importing into downstream analysis tools.
  * `feature-table.tsv`: Tab-separated abundance table for each ASV and each sample.
  * `rel-feature-table.biom`: Relative abundance table in biom format for importing into downstream analysis tools.
  * `rel-feature-table.tsv`: Tab-separated relative abundance table for each ASV and each sample.
  * `table.qza`: QIIME2 data artefact.

### Taxonomic classification

ASV abundance and sequences inferred in DADA2 are informative but routinely taxonomic classifications such as family or genus annotation is desireable. ASV sequences are classified by default against the [SILVA](https://www.arb-silva.de/) [v132](https://www.arb-silva.de/documentation/release-132/) database to add taxonomic information, but a custom database is used if provided. In particular, a [UNITE](https://unite.ut.ee/repository.php) fasta file can be provided to classify fungal ITS sequences.

**Output files:**

* `taxonomy/`
  * `taxonomy.tsv`: Tab-separated table with taxonomic classification for each ASV
  * `index.html`: ASV IDs with taxonomic classification in an interactive table that can be viewed in your web browser

### Exclude taxa

Removes unwanted taxa in DADA2 output sequences and abundance tables by taxonomic classification. Unwanted taxa are often off-targets generated in PCR with primers that are not perfectly specific for the target DNA. For example, PCR with commonly used primers also amplifyies mitrochindrial or chloroplast rRNA genes and therefore leads to non-bacteria products. These mitrochondria or chloroplast amplicons are removed in this step.

All following analysis is based on these filtered tables.

**Output files:**

* `representative_sequences/filtered/`
  * `sequences.fasta`: Fasta file with ASV sequences.
  * `index.html`: ASV IDs, sequences and blast results in an interactive table that can be viewed in your web browser.
  * `rep-seqs.qza`: QIIME2 data artefact.
* `abundance-table/filtered/`
  * `abs-abund-table-2.tsv`: Tab-separated absolute abundance table at phylum level.
  * `abs-abund-table-3.tsv`: Tab-separated absolute abundance table at class level.
  * `abs-abund-table-4.tsv`: Tab-separated absolute abundance table at order level.
  * `abs-abund-table-5.tsv`: Tab-separated absolute abundance table at family level.
  * `abs-abund-table-6.tsv`: Tab-separated absolute abundance table at genus level.
  * `abs-abund-table-7.tsv`: Tab-separated absolute abundance table at species level.
  * `count_table_filter_stats.tsv`: Tab-separated table with information on how much counts were filtered for each sample.
  * `feature-table.biom`: Abundance table in biom format for importing into downstream analysis tools.
  * `feature-table.tsv`: Tab-separated abundance table for each ASV and each sample.
  * `table.qza`: QIIME2 data artefact.

### Relative abundance tables

Absolute abundance tables produced by the previous steps contain count data, but the compositional nature of 16S rRNA amplicon sequencing requires sequencing depth normalisation. This step computes relative abundance tables for various taxonomic levels and a detailed table for all ASVs with taxonomic classification, sequence and relative abundance for each sample. Typically used for in depth investigation of taxa abundances.

**Output files:**

* `rel_abundance_tables/`
  * `rel-table-2.tsv`: Tab-separated relative abundance table at phylum level.
  * `rel-table-3.tsv`: Tab-separated relative abundance table at class level.
  * `rel-table-4.tsv`: Tab-separated relative abundance table at order level.
  * `rel-table-5.tsv`: Tab-separated relative abundance table at family level.
  * `rel-table-6.tsv`: Tab-separated relative abundance table at genus level.
  * `rel-table-7.tsv`: Tab-separated relative abundance table at species level.
  * `rel-table-ASV.tsv`: Tab-separated relative abundance table for all ASVs.
  * `qiime2_ASV_table.tsv`: Tab-separated table for all ASVs with taxonomic classification, sequence and relative abundance.

### Barplot

Produces an interactive abundance plot count tables that aids exploratory browsing the discovered taxa and their abundance in samples and allows sorting for associated meta data.

**Output files:**

* `barplot/`
  * `index.html`: Interactive barplot for taxa abundance per sample that can be viewed in your web browser.

### Alpha diversity rarefaction curves

Produces rarefaction plots for several alpha diversity indices, and is primarily used to determine if the richness of the samples has been fully observed or sequenced. If the slope of the curves does not level out and the lines do not becomes horizontal, this might be because the sequencing depth was too low to observe all diversity or that sequencing error artificially increases sequence diversity and causes false discoveries.

**Output files:**

* `alpha-rarefaction/`
  * `index.html`: Interactive alphararefaction curve for taxa abundance per sample that can be viewed in your web browser.

### Alpha diversity indices

Alpha diversity measures the species diversity within samples. This step calculates alpha diversity using various methods and performs pairwise comparisons of groups of samples.

**Output files:**

* `alpha-diversity`
  * `evenness_vector/index.html`: Pielou’s Evenness.
  * `faith_pd_vector/index.html`: Faith’s Phylogenetic Diversity (qualitiative, phylogenetic).
  * `observed_otus_vector/index.html`: Observed OTUs (qualitative).
  * `shannon_vector/index.html`: Shannon’s diversity index (quantitative).

### Beta diversity indices

Beta diversity measures the species community differences between samples. This step calculates beta diversity distances using various methods and performs pairwise comparisons of groups of samples. Additionally principle coordinates analysis (PCoA) plots are produced that can be visualized with [Emperor](https://biocore.github.io/emperor/build/html/index.html) in your default browser without the need for installation.

**The following methods are used to calculate community dissimilarities:**

* Jaccard distance (qualitative)
* Bray-Curtis distance (quantitative)
* unweighted UniFrac distance (qualitative, phylogenetic)
* weighted UniFrac distance (quantitative, phylogenetic)

**Output files:**

* `beta-diversity/`
  * `<method>_distance_matrix-<treatment>/index.html`
  * `<method>_pcoa_results-PCoA/index.html`
    * method: bray_curtis, jaccard, unweighted_unifrac, weighted_unifrac
    * treatment: depends on your metadata sheet or what metadata categories you have specified

### ANCOM

Analysis of Composition of Microbiomes ([ANCOM](https://www.ncbi.nlm.nih.gov/pubmed/26028277)) is applied to identify features that are differentially abundant across sample groups. A key assumption made by ANCOM is that few taxa (less than about 25%) will be differentially abundant between groups otherwise the method will be inaccurate.

ANCOM is applied to each suitable or specified metadata column for 6 taxonomic levels.

**Output files:**

* `ancom/`
  * `Category-<treatment>-<taxonomic level>/index.html`
    * treatment: depends on your metadata sheet or what metadata categories you have specified
    * taxonomic level: level-2 (phylum), level-3 (class), level-4 (order), level-5 (family), level-6 (genus), ASV

## Pipeline information

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

**Output files:**

* `pipeline_info/`
  * Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  * Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.csv`.
  * Documentation for interpretation of results in HTML format: `results_description.html`.

## More help

QIIME2 is currently **under heavy development** and often updated, this version of ampliseq uses QIIME2 2019.10. QIIME2 has excellent support in its [forum](https://docs.qiime2.org/2019.10/).

## Citations

Besides citing the [pipeline](https://doi.org/10.5281/zenodo.3568091) and its [publication](https://doi.org/10.3389/fmicb.2020.550420), all tools that were used inside the pipeline have to be cited in a publication properly:

* FastQC, "Andrews, Simon. "FastQC: a quality control tool for high throughput sequence data." (2010)."
* Cutadapt "Martin, Marcel. "Cutadapt removes adapter sequences from high-throughput sequencing reads." EMBnet. journal 17.1 (2011): pp-10."
* MultiQC, "Ewels, Philip, et al. "MultiQC: summarize analysis results for multiple tools and samples in a single report." Bioinformatics 32.19 (2016): 3047-3048."
* QIIME2, "Bolyen, Evan, et al. "Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2." Nature Biotechnology 37 (2019): 852–857."
* DADA2, "Callahan, Benjamin J., et al. "DADA2: high-resolution sample inference from Illumina amplicon data." Nature methods 13.7 (2016): 581."
* Matplotlib, "Hunter, John D. "Matplotlib: A 2D graphics environment." Computing in science & engineering 9.3 (2007): 90-95."
* Feature-classifier, "Bokulich, Kaehler, et al. "Optimizing taxonomic classification of marker-gene amplicon sequences with QIIME 2's q2-feature-classifier plugin." Microbiome 6 (2018): 90.
* SILVA database, "Quast, Pruesse, et al. 2013. 'The SILVA ribosomal RNA gene database project: improved data processing and web-based tools', Nucleic Acids Research, 41: D590-D96."
* Mafft, "Katoh, Kazutaka and Standley, Daron M. "MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Molecular biology and evolution 4 (2013): 772-780"
* ANCOM, "Mandal, Siddhartha et al. “Analysis of composition of microbiomes: a novel method for studying microbial composition” Microbial ecology in health and disease vol. 26 27663. 29 May. 2015, doi:10.3402/mehd.v26.27663"
