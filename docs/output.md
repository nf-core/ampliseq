# nf-core/ampliseq: Output

This document describes the output produced by the pipeline.


## Pipeline overview
The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [FastQC](#fastqc) - read quality control
* [Cutadapt](#cutadapt) - primer trimming
* [MultiQC](#multiqc) - aggregate report for FastQC and Cutadapt
* [QIIME2](#qiime2) - import trimmed read into QIIME2, most following steps are wrapped by QIIME2
* [DADA2](#dada2) - denoise reads and infer amplicon sequencing variants
* [Exclude taxa](#exclude-taxa) - exclude unwanted taxa that were PCR off-targets
* [Relative abundance tables](#relative-abundance-tables) - produce relative abundance tables for various taxonomic levels
* [Barplot](#barplot) - produces an interactive abundance plot count tables
* [Alpha diversity rarefaction curves](#alpha-diversity-rarefaction-curves) - plots alpha rarefaction curves
* [Alpha diversity indices](#alpha-diversity-indices) - computes alpha diversity indices and pairwise comparisons
* [Beta diversity indices](#beta-diversity-indices) - computes beta diversity indices and PCoA plots
* [ANCOM](#ancom) - calls differentially abundant taxa
* [More help](#more-help) - more help and further reading
* [Citations](#citations) - Citations for tools involved

## FastQC
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other overrepresented sequences.

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

**Output directory: `results/fastqc`**

* `sample_fastqc.html`
  * FastQC report, containing quality metrics for your untrimmed raw fastq files
* `zips/sample_fastqc.zip`
  * zip file containing the FastQC report, tab-delimited data file and plot images

## Cutadapt
[Cutadapt](https://journal.embnet.org/index.php/embnetjournal/article/view/200) is trimming primer sequences from sequencing reads. Primer sequences are non-biological sequences that often introduce point mutations that do not reflect sample sequences. This is especially true for degenerated PCR primer. If primer trimming would be omitted, artifactual amplicon sequence variants might be computed by the denoising tool or sequences might be lost due to become labelled as PCR chimera.

**Output directory: `results/trimmed/logs`**
  * Log files with retained reads, trimming percentage, etc. for each sample.


## MultiQC
[MultiQC](http://multiqc.info) is a visualisation tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in within the report data directory.

The pipeline has special steps which allow the software versions used to be reported in the MultiQC output for future traceability.

**Output directory: `results/multiqc`**

* `Project_multiqc_report.html`
  * MultiQC report - a standalone HTML file that can be viewed in your web browser
* `Project_multiqc_data/`
  * Directory containing parsed statistics from the different tools used in the pipeline

For more information about how to use MultiQC reports, see http://multiqc.info

## QIIME2
**Quantitative Insights Into Microbial Ecology 2** ([QIIME2](https://qiime2.org/)) is a next-generation microbiome bioinformatics platform and the successor of the widely used [QIIME1](https://www.nature.com/articles/nmeth.f.303). QIIME2 is currently **under heavy development** and often updated, this version of ampliseq uses QIIME2 2018.6. QIIME2 has a wide variety of analysis tools available and has excellent support in its [forum](https://docs.qiime2.org/2018.6/).

At this point of the analysis the trimmed reads are imported into QIIME2 and an interactive quality plot is made.

**Output directory: `results/demux`**
* `index.html`
  * Quality plots that can be viewed in your web browser
* `demux.qza` (only when --untilQ2import is true)
  * QIIME2 artefact for imported reads

All following steps are performed in QIIME2.

## DADA2
[DADA2](https://www.nature.com/articles/nmeth.3869) performs fast and accurate sample inference from amplicon data with single-nucleotide resolution. It infers exact amplicon sequence variants (ASVs) from amplicon data with fewer false positives than many other methods while maintaining high sensitivity.

DADA2 computes an error model on the sequencing reads (forward and reverse independently), therefore quality filtering or paired read merging may not be performed before. Each sequencing run varies in their error profile and it is recommended that DADA2 runs separately on data from each run individually. It is recommended to use the ampliseq option `--multipleSequencingRuns` to analyse such data.

DADA2 reduces sequence errors and dereplicates sequences by quality filtering, denoising, read pair merging and PCR chimera removal.

**Output directory: `results/representative_sequences/unfiltered`**
* `sequences.fasta`
  * Fasta file with ASV sequences
* `index.html`
  * ASV IDs, sequences and blast results in an interactive table that can be viewed in your web browser
* `rep-seqs.qza`
  * QIIME2 data artefact

**Output directory: `results/abundance-table/unfiltered`**
* `dada_report.txt`
  * DADA2 verbose output
* `dada_stats.tsv`
  * Tab-separated table of DADA2 statistics
* `feature-table.biom`
  * Abundance table in biom format for importing into downstream analysis tools
* `feature-table.tsv`
  * Tab-separated abundance table for each ASV and each sample
* `rel-feature-table.biom`
  * Relative abundance table in biom format for importing into downstream analysis tools
* `rel-feature-table.tsv`
  * Tab-separated relative abundance table for each ASV and each sample
* `table.qza`
  * QIIME2 data artefact

## Taxonomic classification
ASV abundance and sequences inferred in DADA2 are informative but routinely taxonomic classifications such as family or genus annotation is desireable. ASV sequences are classified by default against the [SILVA](https://www.arb-silva.de/) [v132](https://www.arb-silva.de/documentation/release-132/) database to add taxonomic information.

**Output directory: `results/taxonomy`**
* `taxonomy.tsv`
  * Tab-separated table with taxonomic classification for each ASV
* `index.html`
  * ASV IDs with taxonomic classification in an interactive table that can be viewed in your web browser

## Exclude taxa
Removes unwanted taxa in DADA2 output sequences and abundance tables by taxonomic classification. Unwanted taxa are often off-targets generated in PCR with primers that are not perfectly specific for the target DNA. For example, PCR with commonly used primers also amplifyies mitrochindrial or chloroplast rRNA genes and therefore leads to non-bacteria products. These mitrochondria or chloroplast amplicons are removed in this step.

All following analysis is based on these filtered tables.

**Output directory: `results/representative_sequences/filtered`**
* `sequences.fasta`
  * Fasta file with ASV sequences
* `index.html`
  * ASV IDs, sequences and blast results in an interactive table that can be viewed in your web browser
* `rep-seqs.qza`
  * QIIME2 data artefact

**Output directory: `results/abundance-table/filtered`**
* `count_table_filter_stats.csv`
  * Comma-separated table with information on how much counts were filtered for each sample
* `feature-table.biom`
  * Abundance table in biom format for importing into downstream analysis tools
* `feature-table.tsv`
  * Tab-separated abundance table for each ASV and each sample
* `table.qza`
  * QIIME2 data artefact

## Relative abundance tables
Absolute abundance tables produced by the previous steps contain count data, but the compositional nature of 16S rRNA amplicon sequencing requires sequencing depth normalisation. This step computes relative abundance tables for various taxonomic levels and a detailed table for all ASVs with taxonomic classification, sequence and relative abundance for each sample. Typically used for in depth investigation of taxa abundances.

**Output directory: `results/rel_abundance_tables`**
* `rel-table-2.tsv`
  * Tab-separated relative abundance table at phylum level
* `rel-table-3.tsv`
  * Tab-separated relative abundance table at class level
* `rel-table-4.tsv`
  * Tab-separated relative abundance table at order level
* `rel-table-5.tsv`
  * Tab-separated relative abundance table at family level
* `rel-table-6.tsv`
  * Tab-separated relative abundance table at genus level
* `rel-table-7.tsv`
  * Tab-separated relative abundance table at species level
* `rel-table-ASV.tsv`
  * Tab-separated relative abundance table for all ASVs
* `qiime2_ASV_table.csv`
  * Tab-separated table for all ASVs with taxonomic classification, sequence and relative abundance

## Barplot
Produces an interactive abundance plot count tables that aids exploratory browsing the discovered taxa and their abundance in samples and allows sorting for associated meta data.

**Output directory: `results/barplot`**
* `index.html`
  * Interactive barplot for taxa abundance per sample that can be viewed in your web browser

## Alpha diversity rarefaction curves
Produces rarefaction plots for several alpha diversity indices, and is primarily used to determine if the richness of the samples has been fully observed or sequenced. If the slope of the curves does not level out and the lines do not becomes horizontal, this might be because the sequencing depth was too low to observe all diversity or that sequencing error artificially increases sequence diversity and causes false discoveries.

**Output directory: `results/alpha-rarefaction`**
* `index.html`
  * Interactive alphararefaction curve for taxa abundance per sample that can be viewed in your web browser

## Alpha diversity indices
Alpha diversity measures the species diversity within samples. This step calculates alpha diversity using various methods and performs pairwise comparisons of groups of samples.

**Output directory: `results/alpha-diversity`** (all *.html files can be viewed in your web browser)
* `evenness_vector/index.html`
  * Pielou’s Evenness
* `faith_pd_vector/index.html`
  * Faith’s Phylogenetic Diversity (qualitiative, phylogenetic)
* `observed_otus_vector/index.html`
  * Observed OTUs (qualitative)
* `shannon_vector/index.html`
  * Shannon’s diversity index (quantitative)


## Beta diversity indices
Beta diversity measures the species community differences between samples. This step calculates beta diversity distances using various methods and performs pairwise comparisons of groups of samples. Additionally principle coordinates analysis (PCoA) plots are produced that can be visualized with [Emperor](https://biocore.github.io/emperor/build/html/index.html) in your default browser without the need for installation.

**The following methods are used to calculate community dissimilarities:**
* Jaccard distance (qualitative)
* Bray-Curtis distance (quantitative)
* unweighted UniFrac distance (qualitative, phylogenetic)
* weighted UniFrac distance (quantitative, phylogenetic)

**Output directory: `results/beta-diversity`** (all *.html files can be viewed in your web browser)
* `<method>_distance_matrix-<treatment>/index.html`
* `<method>_pcoa_results-PCoA/index.html`
  * method: bray_curtis, jaccard, unweighted_unifrac, weighted_unifrac
  * treatment: depends on your metadata sheet or what metadata categories you have specified

## ANCOM
Analysis of Composition of Microbiomes ([ANCOM](https://www.ncbi.nlm.nih.gov/pubmed/26028277)) is applied to identify features that are differentially abundant across sample groups. A key assumption made by ANCOM is that few taxa (less than about 25%) will be differentially abundant between groups otherwise the method will be inaccurate.

ANCOM is applied to each suitable or specified metadata column for 6 taxonomic levels.

**Output directory: `results/ancom`** (all *.html files can be viewed in your web browser)
* `Category-<treatment>-<taxonomic level>/index.html`
  * treatment: depends on your metadata sheet or what metadata categories you have specified
  * taxonomic level: level-2 (phylum), level-3 (class), level-4 (order), level-5 (family), level-6 (genus), ASV

# More help
QIIME2 is currently **under heavy development** and often updated, this version of ampliseq uses QIIME2 2018.6. QIIME2 has excellent support in its [forum](https://docs.qiime2.org/2018.6/).

# Citations
All tools inside the pipeline have to be cited in a publication properly:

* FastQC, "Andrews, Simon. "FastQC: a quality control tool for high throughput sequence data." (2010)."
* Cutadapt "Martin, Marcel. "Cutadapt removes adapter sequences from high-throughput sequencing reads." EMBnet. journal 17.1 (2011): pp-10."
* MultiQC, "Ewels, Philip, et al. "MultiQC: summarize analysis results for multiple tools and samples in a single report." Bioinformatics 32.19 (2016): 3047-3048."
* QIIME2, "Bolyen, Evan, et al. QIIME 2: Reproducible, interactive, scalable, and extensible microbiome data science. No. e27295v1. PeerJ Preprints, 2018."
* DADA2, "Callahan, Benjamin J., et al. "DADA2: high-resolution sample inference from Illumina amplicon data." Nature methods 13.7 (2016): 581."
* Matplotlib, "Hunter, John D. "Matplotlib: A 2D graphics environment." Computing in science & engineering 9.3 (2007): 90-95."
* Feature-classifier, "Bokulich, Kaehler, et al. "Optimizing taxonomic classification of marker-gene amplicon sequences with QIIME 2's q2-feature-classifier plugin." Microbiome 6 (2018): 90. 
* SILVA database, "Quast, Pruesse, et al. 2013. 'The SILVA ribosomal RNA gene database project: improved data processing and web-based tools', Nucleic Acids Research, 41: D590-D96."
* Mafft, "Katoh, Kazutaka and Standley, Daron M. "MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Molecular biology and evolution 4 (2013): 772-780"
* ANCOM, "Mandal, Siddhartha et al. “Analysis of composition of microbiomes: a novel method for studying microbial composition” Microbial ecology in health and disease vol. 26 27663. 29 May. 2015, doi:10.3402/mehd.v26.27663"