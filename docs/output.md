# nf-core/ampliseq: Output

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
    * [DADA2](#dada2) - Infer Amplicon Sequence Variants (ASVs) and taxonomic classification
      * [ITSx](#itsx) - Optionally, taxonomic classification can be performed on ITS region only
    * [QIIME2](#qiime2) - Secondary analysis
      * [Taxonomic classification](#taxonomic-classification) - Taxonomical classification of ASVs
      * [Exclude taxa](#exclude-taxa) - Remove unwanted ASV based on taxonomy
      * [Relative abundance tables](#relative-abundance-tables) - Exported relative abundance tables
      * [Barplot](#barplot) - Interactive barplot
      * [Alpha diversity rarefaction curves](#alpha-diversity-rarefaction-curves) - Rarefaction curves for quality control
      * [Alpha diversity indices](#alpha-diversity-indices) - Diversity within samples
      * [Beta diversity indices](#beta-diversity-indices) - Diversity between samples (e.g. PCoA plots)
      * [ANCOM](#ancom) - Differential abundance analysis
    * [Read count report](#Read-count-report) - Report of read counts during various steps of the pipeline
    * [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution
  * [Citations](#citations)

### FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences.

For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

**Output files:**

* `fastqc/`
  * `*_fastqc.html`: FastQC report containing quality metrics for your untrimmed raw fastq files.

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality.

### Cutadapt

[Cutadapt](https://journal.embnet.org/index.php/embnetjournal/article/view/200) is trimming primer sequences from sequencing reads. Primer sequences are non-biological sequences that often introduce point mutations that do not reflect sample sequences. This is especially true for degenerated PCR primer. If primer trimming would be omitted, artifactual amplicon sequence variants might be computed by the denoising tool or sequences might be lost due to become labelled as PCR chimera.

**Output files:**

* `cutadapt/`: directory containing log files with retained reads, trimming percentage, etc. for each sample.
  * `cutadapt_summary.tsv`: Summary of read numbers that pass cutadapt.

### MultiQC

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarizing all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability.

For more information about how to use MultiQC reports, see [https://multiqc.info](https://multiqc.info).

**Output files:**

* `multiqc/`
  * `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  * `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.

### DADA2

[DADA2](https://www.nature.com/articles/nmeth.3869) performs fast and accurate sample inference from amplicon data with single-nucleotide resolution. It infers exact amplicon sequence variants (ASVs) from amplicon data with fewer false positives than many other methods while maintaining high sensitivity.

DADA2 computes an error model on the sequencing reads (forward and reverse independently), therefore quality filtering or paired read merging may not be performed before. Each sequencing run varies in their error profile and it is recommended that DADA2 runs separately on data from each run individually. It is recommended to use the ampliseq option `--multiple_sequencing_runs` to analyse such data.

DADA2 reduces sequence errors and dereplicates sequences by quality filtering, denoising, read pair merging (for paired end Illumina reads only) and PCR chimera removal.

Additionally, DADA2 taxonomically classifies the ASVs using pre-trained databases.

**Output files:**

* `dada2/`
  * `ASV_seqs.fasta`: Fasta file with ASV sequences.
  * `ASV_table.tsv`: Counts for each ASV sequence.
  * `ASV_tax.tsv`: Taxonomic classification for each ASV sequence.
  * `ASV_tax_species.tsv`: Species classification for each ASV sequence.
  * `DADA2_stats.tsv`: Tracking read numbers through DADA2 processing steps, for each sample.
  * `DADA2_table.rds`: DADA2 ASV table as R object.
  * `DADA2_tables.tsv`: DADA2 ASV table.
* `dada2/args/`: Directory containing all parameters for DADA2 steps.
* `dada2/log/`: Directory containing log files for DADA2 steps.
* `dada2/QC/`
  * `*.err.convergence.txt`: Convergence values for DADA2's dada command, should reduce over several magnitudes and approaching 0.
  * `*.err.pdf`: Estimated error rates for each possible transition. The black line shows the estimated error rates after convergence of the machine-learning algorithm. The red line shows the error rates expected under the nominal definition of the Q-score. The estimated error rates (black line) should be a good fit to the observed rates (points), and the error rates should drop with increased quality.
  * `*_qual_stats.pdf`: Overall read quality profiles: heat map of the frequency of each quality score at each base position. The mean quality score at each position is shown by the green line, and the quartiles of the quality score distribution by the orange lines. The red line shows the scaled proportion of reads that extend to at least that position.

#### ITSx

Optionally, the ITS region can be extracted from each ASV sequence using ITSx, and taxonomic classification is performed based on the ITS sequence.

**Output files:**

* `itsx/`
  * `ASV_ITS_seqs.full.fasta`: Fasta file with ITS region from each ASV sequence.
* `dada2/`
  * `ASV_ITS_tax.tsv`: Taxonomic classification with ITS region of each ASV sequence.
  * `ASV_ITS_tax_species.tsv`: Species classification with ITS region of each ASV sequence.

### QIIME2

**Quantitative Insights Into Microbial Ecology 2** ([QIIME2](https://qiime2.org/)) is a next-generation microbiome bioinformatics platform and the successor of the widely used [QIIME1](https://www.nature.com/articles/nmeth.f.303).

ASV sequences and counts as produced before with DADA2 are imported into QIIME2 and further analysed. First, ASVs are taxonomically classified, than filtered (`--exclude_taxa`, `--min_frequency`, `--min_samples`), and abundance tables exported. Following, diversity indices are calculated and testing for differential abundant features between sample groups is performed.

#### Taxonomic classification

ASV abundance and sequences inferred in DADA2 are informative but routinely taxonomic classifications such as family or genus annotation is desireable.

**Output files:**

* `taxonomy/`
  * `taxonomy.tsv`: Tab-separated table with taxonomic classification for each ASV
  * `*-classifier.qza`: QIIME2 artefact of the trained classifier. Can be supplied to other pipeline runs with `--classifier`

#### Exclude taxa

Removes unwanted taxa in DADA2 output sequences and abundance tables by taxonomic classification. Unwanted taxa are often off-targets generated in PCR with primers that are not perfectly specific for the target DNA. For example, PCR with commonly used primers also amplifyies mitrochindrial or chloroplast rRNA genes and therefore leads to non-bacteria products. These mitrochondria or chloroplast amplicons are removed in this step by default (`--exclude_taxa`).

All following analysis is based on these filtered tables.

**Output files:**

* `qiime2/representative_sequences/`
  * `rep-seq.fasta`: Fasta file with ASV sequences.
  * `descriptive_stats.tsv`: Length, mean, etc. of ASV sequences.
  * `seven_number_summary.tsv`: Length of ASV sequences in different quantiles.
* `qiime2/abundance_tables/`
  * `abs-abund-table-2.tsv`: Tab-separated absolute abundance table at phylum level.
  * `abs-abund-table-3.tsv`: Tab-separated absolute abundance table at class level.
  * `abs-abund-table-4.tsv`: Tab-separated absolute abundance table at order level.
  * `abs-abund-table-5.tsv`: Tab-separated absolute abundance table at family level.
  * `abs-abund-table-6.tsv`: Tab-separated absolute abundance table at genus level.
  * `abs-abund-table-7.tsv`: Tab-separated absolute abundance table at species level.
  * `count_table_filter_stats.tsv`: Tab-separated table with information on how much counts were filtered for each sample.
  * `feature-table.biom`: Abundance table in biom format for importing into downstream analysis tools.
  * `feature-table.tsv`: Tab-separated abundance table for each ASV and each sample.

#### Relative abundance tables

Absolute abundance tables produced by the previous steps contain count data, but the compositional nature of 16S rRNA amplicon sequencing requires sequencing depth normalisation. This step computes relative abundance tables for various taxonomic levels and a detailed table for all ASVs with taxonomic classification, sequence and relative abundance for each sample. Typically used for in depth investigation of taxa abundances.

**Output files:**

* `qiime2/rel_abundance_tables/`
  * `rel-table-2.tsv`: Tab-separated relative abundance table at phylum level.
  * `rel-table-3.tsv`: Tab-separated relative abundance table at class level.
  * `rel-table-4.tsv`: Tab-separated relative abundance table at order level.
  * `rel-table-5.tsv`: Tab-separated relative abundance table at family level.
  * `rel-table-6.tsv`: Tab-separated relative abundance table at genus level.
  * `rel-table-7.tsv`: Tab-separated relative abundance table at species level.
  * `rel-table-ASV.tsv`: Tab-separated relative abundance table for all ASVs.
  * `qiime2_ASV_table.tsv`: Tab-separated table for all ASVs with taxonomic classification, sequence and relative abundance.

#### Barplot

Produces an interactive abundance plot count tables that aids exploratory browsing the discovered taxa and their abundance in samples and allows sorting for associated meta data.

**Output files:**

* `qiime2/barplot/`
  * `index.html`: Interactive barplot for taxa abundance per sample that can be viewed in your web browser.

#### Alpha diversity rarefaction curves

Produces rarefaction plots for several alpha diversity indices, and is primarily used to determine if the richness of the samples has been fully observed or sequenced. If the slope of the curves does not level out and the lines do not becomes horizontal, this might be because the sequencing depth was too low to observe all diversity or that sequencing error artificially increases sequence diversity and causes false discoveries.

**Output files:**

* `qiime2/alpha-rarefaction/`
  * `index.html`: Interactive alphararefaction curve for taxa abundance per sample that can be viewed in your web browser.

#### Alpha diversity indices

Alpha diversity measures the species diversity within samples. Diversity calculations are based on sub-sampled data rarefied to the minimum read count of all samples. This step calculates alpha diversity using various methods and performs pairwise comparisons of groups of samples. It is based on a phylogenetic tree of all ASV sequences.

**Output files:**

* `qiime2/phylogenetic_tree/`
  * `tree.nwk`: Phylogenetic tree in newick format.
  * `rooted-tree.qza`: Phylogenetic tree in QIIME2 format.
* `qiime2/diversity/`
  * `*.txt`: File that describes the rarefaction depth (file name and file contant).
* `qiime2/diversity/alpha_diversity/`
  * `evenness_vector/index.html`: Pielou’s Evenness.
  * `faith_pd_vector/index.html`: Faith’s Phylogenetic Diversity (qualitiative, phylogenetic).
  * `observed_otus_vector/index.html`: Observed OTUs (qualitative).
  * `shannon_vector/index.html`: Shannon’s diversity index (quantitative).

#### Beta diversity indices

Beta diversity measures the species community differences between samples. Diversity calculations are based on sub-sampled data rarefied to the minimum read count of all samples. This step calculates beta diversity distances using various methods and performs pairwise comparisons of groups of samples. Additionally principle coordinates analysis (PCoA) plots are produced that can be visualized with [Emperor](https://biocore.github.io/emperor/build/html/index.html) in your default browser without the need for installation. This calculations are based on a phylogenetic tree of all ASV sequences.

**The following methods are used to calculate community dissimilarities:**

* Jaccard distance (qualitative)
* Bray-Curtis distance (quantitative)
* unweighted UniFrac distance (qualitative, phylogenetic)
* weighted UniFrac distance (quantitative, phylogenetic)

**Output files:**

* `qiime2/phylogenetic_tree/`
  * `tree.nwk`: Phylogenetic tree in newick format.
  * `rooted-tree.qza`: Phylogenetic tree in QIIME2 format.
* `qiime2/diversity/`
  * `*.txt`: File that describes the rarefaction depth (file name and file contant).
* `qiime2/diversity/beta_diversity/`
  * `<method>_distance_matrix-<treatment>/index.html`
  * `<method>_pcoa_results-PCoA/index.html`
    * method: bray_curtis, jaccard, unweighted_unifrac, weighted_unifrac
    * treatment: depends on your metadata sheet or what metadata categories you have specified

#### ANCOM

Analysis of Composition of Microbiomes ([ANCOM](https://www.ncbi.nlm.nih.gov/pubmed/26028277)) is applied to identify features that are differentially abundant across sample groups. A key assumption made by ANCOM is that few taxa (less than about 25%) will be differentially abundant between groups otherwise the method will be inaccurate.

ANCOM is applied to each suitable or specified metadata column for 6 taxonomic levels.

**Output files:**

* `qiime2/ancom/`
  * `Category-<treatment>-<taxonomic level>/index.html`
    * treatment: depends on your metadata sheet or what metadata categories you have specified
    * taxonomic level: level-2 (phylum), level-3 (class), level-4 (order), level-5 (family), level-6 (genus), ASV

## Read count report

This report includes information on how many reads per sample passed each pipeline step in which a loss can occur. Specifically, how many read pairs entered cutadapt, were reverse complemented, passed trimming; how many read pairs entered DADA2, were denoised, merged and non-chimeric; and how many counts were lost during excluding unwanted tax and removing low abundance/prevalence sequences in QIIME2.

**Output files:**

* `overall_summary.tsv`

## Pipeline information

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

**Output files:**

* `pipeline_info/`
  * Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  * Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.csv`.
  * Documentation for interpretation of results in HTML format: `results_description.html`.

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
