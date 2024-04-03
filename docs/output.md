---
output:
  pdf_document: default
  html_document: default
---

# nf-core/ampliseq: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Input](#input) - Input files
- [Pipeline summary report](#pipeline-summary-report) - Overview of pipeline output
- [Preprocessing](#preprocessing)
  - [FastQC](#fastqc) - Read quality control
  - [Cutadapt](#cutadapt) - Primer trimming
  - [MultiQC](#multiqc) - Aggregate report describing results
- [ASV inferrence with DADA2](#asv-inferrence-with-dada2) - Infer Amplicon Sequence Variants (ASVs)
- [Optional ASV filtering](#optional-asv-filtering) - Filter ASVs to optimize downstream analysis
  - [VSEARCH cluster](#vsearch-cluster) - Centroid fasta file, filtered asv table, and stats
  - [Barrnap](#barrnap) - Predict ribosomal RNA sequences and optional filtering
  - [Length filter](#length-filter) - Optionally, ASV can be filtered by length thresholds
  - [Codons](#codons) - Optionally the ASVs can be filtered by presence of stop codons.
  - [ITSx](#itsx) - Optionally, the ITS region can be extracted
- [Taxonomic classification](#taxonomic-classification) - Taxonomic classification of (filtered) ASVs
  - [DADA2](#dada2) - Taxonomic classification with DADA2
  - [assignSH](#assignsh) - Optionally, a UNITE species hypothesis (SH) can be added to the DADA2 taxonomy
  - [SINTAX](#sintax) - Taxonomic classification with SINTAX
  - [Kraken2](#kraken2) - Taxonomic classification with Kraken2
  - [QIIME2](#qiime2) - Taxonomic classification with QIIME2
- [Phlogenetic placement and taxonomic classification](#phylogenetic-placement-and-taxonomic-classification) - Placing ASVs into a phyloenetic tree
- [Multiple region analysis with Sidle](#multiple-region-analysis-with-sidle) - Scaffolding multiple regions along a reference
- [Secondary analysis with QIIME2](#secondary-analysis-with-qiime2) - Visualisations, diversity and differential abundance analysis with QIIME2
  - [Abundance tables](#abundance-tables) - Exported abundance tables
  - [Relative abundance tables](#relative-abundance-tables) - Exported relative abundance tables
  - [Barplot](#barplot) - Interactive barplot
  - [Alpha diversity rarefaction curves](#alpha-diversity-rarefaction-curves) - Rarefaction curves for quality control
  - [Diversity analysis](#diversity-analysis) - High level overview with different diversity indices
  - [ANCOM](#ancom) - Differential abundance analysis
- [PICRUSt2](#picrust2) - Predict the functional potential of a bacterial community
- [SBDI export](#sbdi-export) - Swedish Biodiversity Infrastructure (SBDI) submission file
- [Phyloseq](#phyloseq) - Phyloseq R objects
- [Read count report](#read-count-report) - Report of read counts during various steps of the pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### Input

Samplesheet, ASV fasta, and metadata file are copied into the results folder.

<details markdown="1">
<summary>Output files</summary>

- `input/`
  - `*`: Samplesheet input if specified with `--input`.
  - `*.tsv`: Metadata input if specified with `--metadata`.
  - `*`: ASV sequence input if specified with `--input_fasta`.

</details>

### Pipeline summary report

A summary report for most pipeline results in html format produced by [R Markdown](https://rmarkdown.rstudio.com/). The report gives a general overview of the analysis, includes many tables and visualizations, and links to interactive downstream analysis results, if available.

<details markdown="1">
<summary>Output files</summary>

- `summary_report/`
  - `summary_report.html`: pipeline summary report as standalone HTML file that can be viewed in your web browser.
  - `*.svg*`: plots that were produced for (and are included in) the report.
  - `versions.yml`: software versions used to produce this report.

</details>

### Preprocessing

#### FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics for your untrimmed raw fastq files.

</details>

#### Cutadapt

[Cutadapt](https://journal.embnet.org/index.php/embnetjournal/article/view/200) is trimming primer sequences from sequencing reads. Primer sequences are non-biological sequences that often introduce point mutations that do not reflect sample sequences. This is especially true for degenerated PCR primer. If primer trimming would be omitted, artifactual amplicon sequence variants might be computed by the denoising tool or sequences might be lost due to become labelled as PCR chimera.

<details markdown="1">
<summary>Output files</summary>

- `cutadapt/`: directory containing log files with retained reads, trimming percentage, etc. for each sample.
  - `cutadapt_summary.tsv`: Summary of read numbers that pass cutadapt.
  - `assignTaxonomy.cutadapt.log`: Contains how many expected amplified sequences were extracted from the DADA2 reference taxonomy database. Optional.

</details>

#### MultiQC

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

:::note
The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality.
:::

### ASV inferrence with DADA2

[DADA2](https://www.nature.com/articles/nmeth.3869) performs fast and accurate sample inference from amplicon data with single-nucleotide resolution. It infers exact amplicon sequence variants (ASVs) from amplicon data with fewer false positives than many other methods while maintaining high sensitivity.

DADA2 computes an error model on the sequencing reads (forward and reverse independently), therefore quality filtering or paired read merging may not be performed before. Each sequencing run varies in their error profile and it is recommended that DADA2 runs separately on data from each run individually. It is recommended to use the ampliseq option `--multiple_sequencing_runs` to analyse such data.

DADA2 reduces sequence errors and dereplicates sequences by quality filtering, denoising, read pair merging (for paired end Illumina reads only) and PCR chimera removal.

<details markdown="1">
<summary>Output files</summary>

- `dada2/`
  - `ASV_seqs.fasta`: Fasta file with ASV sequences.
  - `ASV_table.tsv`: Counts for each ASV sequence.
  - `DADA2_stats.tsv`: Tracking read numbers through DADA2 processing steps, for each sample.
  - `DADA2_table.rds`: DADA2 ASV table as R object.
  - `DADA2_table.tsv`: DADA2 ASV table.
- `dada2/args/`: Directory containing files with all parameters for DADA2 steps.
- `dada2/log/`: Directory containing log files for DADA2 steps.
- `dada2/QC/`
  - `*.err.convergence.txt`: Convergence values for DADA2's dada command, should reduce over several magnitudes and approaching 0.
  - `*.err.pdf`: Estimated error rates for each possible transition. The black line shows the estimated error rates after convergence of the machine-learning algorithm. The red line shows the error rates expected under the nominal definition of the Q-score. The estimated error rates (black line) should be a good fit to the observed rates (points), and the error rates should drop with increased quality.
  - `*_qual_stats.pdf`: Overall read quality profiles: heat map of the frequency of each quality score at each base position. The mean quality score at each position is shown by the green line, and the quartiles of the quality score distribution by the orange lines. The red line shows the scaled proportion of reads that extend to at least that position.
  - `*_preprocessed_qual_stats.pdf`: Same as above, but after preprocessing.

</details>

For binned quality scores in NovaSeq data, monotonicity in the fitted error model is enforced when using `--illumina_novaseq`. Consequently, additional QC data is generated.

<details markdown="1">
<summary>Output files</summary>

- `dada2/QC/`
  - `*.md.err.convergence.txt`: Convergence values for DADA2's dada command on monotone decreasing (corrected) quality scores, should reduce over several magnitudes and approaching 0.
  - `*.md.err.pdf`: Estimated error rates for each possible transition on monotone decreasing (corrected) quality scores. The black line shows the estimated error rates after convergence of the machine-learning algorithm. The red line shows the error rates expected under the nominal definition of the Q-score. The estimated error rates (black line) should be a good fit to the observed rates (points), and the error rates should drop with increased quality.

</details>

### Optional ASV filtering

#### VSEARCH cluster

Optionally, VSEARCH can be used to cluster the denoised ASVs. This will be performed before other filtering steps.
This directory will hold the centroid fasta file, the filtered asv count table (after merging non-centroid counts with their respective centroid counts), and a stats table.

<details markdown="1">
<summary>Output files</summary>

- `vsearch_cluster/`
  - `ASV_post_clustering_filtered.fna`: Centroid fasta file.
  - `ASV_post_clustering_filtered.stats.tsv`: Stats table.
  - `ASV_post_clustering_filtered.table.tsv`: ASV table.

</details>

#### Barrnap

Barrnap predicts the location of ribosomal RNA genes in genomes, here it can be used to discriminate rRNA sequences from potential contamination. It supports bacteria (5S,23S,16S), archaea (5S,5.8S,23S,16S), metazoan mitochondria (12S,16S) and eukaryotes (5S,5.8S,28S,18S).

Optionally, ASV sequences can be filtered for rRNA sequences identified by Barrnap with `--filter_ssu` that can take a list of abbreviations of the above supported categories (kingdoms), e.g. `bac,arc,mito,euk`. This filtering takes place after DADA2's ASV computation (i.e. after chimera removal) but _before_ taxonomic classification (also applies to above mentioned taxonomic classification with DADA2, i.e. files `ASV_tax.tsv` & `ASV_tax_species.tsv`).

<details markdown="1">
<summary>Output files</summary>

- `barrnap/`
  - `rrna.<kingdom>.gff`: GFF3 output for rRNA matches per kingdom, where kingdom is one of `bac,arc,mito,euk`.
  - `summary.tsv`: Summary of evalues for each ASV and kingdom
  - `ASV_seqs.ssu.fasta`: Fasta file with filtered ASV sequences, only if `--filter_ssu` is set.
  - `ASV_table.ssu.tsv`: Counts for each filtered ASV sequence, only if `--filter_ssu` is set.
  - `stats.ssu.tsv`: Tracking read numbers through filtering, for each sample, only if `--filter_ssu` is set.

</details>

#### Length filter

Optionally, a length filter can be used to reduce potential contamination after ASV computation. For example with 515f and 806r primers the majority of 16S rRNA amplicon sequences should have a length of 253 bp and amplicons vary significantely are likely spurious.

The minimum ASV length threshold can be set by `--min_len_asv` and the maximum length threshold with `--max_len_asv`. If no threshold is set, the filter (and output) is omitted.

<details markdown="1">
<summary>Output files</summary>

- `asv_length_filter/`
  - `ASV_seqs.len.fasta`: Fasta file with filtered ASV sequences.
  - `ASV_table.len.tsv`: Counts for each filtered ASV sequence.
  - `ASV_len_orig.tsv`: ASV length distribution before filtering.
  - `ASV_len_filt.tsv`: ASV length distribution after filtering.
  - `stats.len.tsv`: Tracking read numbers through filtering, for each sample.

</details>

#### Codons

Optionally, the ASVs can be filtered against the presence of stop codons in the specified open reading frame of the ASV. The filtering step can also filter out ASVs that are not multiple of 3 in length.

Codon filtering can be activated by `--filter_codons`. By default, the codons are calculated and checked from the beginning (`--orf_start 1`) to the end (`--orf_end` not set) of the ASV sequence and the filter also checks if the length is a multiple of 3. By default, the filter screens for the presence of the stop codons `TAA` and `TAG`. To define different stop codons than default, use `--stop_codons`(comma-separated list, e.g. `--stop_codons "TAA,TAG,TGA"`).

<details markdown="1">
<summary>Output files</summary>

- `codon_filter/`
  - `ASV_codon_filtered.fna`: Fasta file of ASV sequences that passes the filter thresholds explained above.
  - `ASV_codon_filtered.table.tsv`: The count table of ASVs that successfully passed through the filter thresholds.
  - `ASV_codon_filtered.list`: List of ASV IDs that pass through the filter thresholds.
  - `codon.filtered.stats.tsv`: Tracking read numbers through filtering, for each sample.

</details>

#### ITSx

Optionally, the ITS region can be extracted from each ASV sequence using ITSx, and taxonomic classification is performed based on the ITS sequence. Only sequences with at minimum 50bp in length are retained.

<details markdown="1">
<summary>Output files</summary>

- `itsx/`
  - `ASV_ITS_seqs.full.fasta`: Fasta file with full ITS region from each ASV sequence.
  - `ASV_ITS_seqs.ITS1.fasta` or `ASV_ITS_seqs.ITS2.fasta`: If using --cut_its "its1" or --cut_its "its2"; fasta file with ITS1 or ITS2 region from each ASV sequence.
  - `ASV_ITS_seqs.full_and_partial.fasta`: If using --its_partial; fasta file with full and partial ITS regions from each ASV sequence.
  - `ASV_ITS_seqs.ITS1.full_and_partial.fasta` or `ASV_ITS_seqs.ITS2.full_and_partial.fasta`: If using --cut_its "its1" or --cut_its "its2" and --its_partial; fasta file with complete and partial ITS1 or ITS2 regions from each ASV sequence.
  - `ASV_ITS_seqs.summary.txt`: Summary information from ITSx.
  - `ITSx.args.txt`: File with parameters passed to ITSx.
  - `ASV_seqs.len.fasta`: Fasta file with filtered ASV sequences.
  - `ASV_len_orig.tsv`: ASV length distribution before filtering.
  - `ASV_len_filt.tsv`: ASV length distribution after filtering.

</details>

### Taxonomic classification

Taxonomic classification of ASVs can be performed with a choice of DADA2, SINTAX, Kraken2 or QIIME2 using supplied databases or user supplied databases (see parameter documentation). By default, DADA2 is used for the classification. The taxonomic classification will be done based on filtered ASV sequences (see above).

#### DADA2

DADA2 is the default method for taxonomy classification. Depending on the reference taxonomy database, sequences can be classified down to species rank. Species classification is reported in columns "Species" using DADA2's assignTaxonomy function or "Species_exact" using DADA2's addSpecies function, the latter only assigns exact sequence matches. Generally, species assignment without exact matches are much less trustworthy than those with exact matches. With short amplicons, e.g. 16S rRNA gene V4 region, the non-exact species annotation is not recommended to be trusted. The longer the ASVs are, the more acceptable is the non-exact species classification, e.g. PacBio (nearly) full length 16S rRNA gene sequences are thought to be trustworthy.

Files when _not_ using ITSx (default):

<details markdown="1">
<summary>Output files</summary>

- `dada2/`
  - `ASV_tax.*.tsv`: Taxonomic classification for each ASV sequence.
  - `ASV_tax_species.*.tsv`: Exact species classification for each ASV sequence.
  - `ref_taxonomy.*.txt`: Information about the used reference taxonomy, such as title, version, citation.

</details>

Files when using ITSx:

<details markdown="1">
<summary>Output files</summary>

- `dada2/`
  - `ASV_ITS_tax.*.tsv`: Taxonomic classification with ITS region of each ASV sequence.
  - `ASV_ITS_tax_species.*.tsv`: Exact species classification with ITS region of each ASV sequence.
  - `ASV_tax.*.tsv`: Taxonomic classification of each ASV sequence, based on the ITS region.
  - `ASV_tax_species.*.tsv`: Exact species classification of each ASV sequence, based on the ITS region.
  - `ref_taxonomy.*.txt`: Information about the used reference taxonomy, such as title, version, citation.

</details>

#### assignSH

Optionally, a UNITE species hypothesis (SH) can be added to the DADA2 taxonomy. In short, the sequences are matched to the selected UNITE database with VSEARCH, using the sequence identity cutoff 98.5%. If a good enough and unique match is found, the sequence is assigned the SH of the matching sequence, and the taxonomy assignment for the sequence is changed to that of the SH. The look-up tables for SH taxonomies are generated as descibed at https://github.com/biodiversitydata-se/unite-shinfo. If no match is found, the taxonomy generated by assignTaxonomy and/or assignSpecies is kept.

<details markdown="1">
<summary>Output files</summary>

- `assignsh/`
  - `ASV_tax_species_SH.tsv`: Taxonomic classification with SH taxonomy added in case of a match.
  - `*.vsearch.txt`: Raw vsearch results.

</details>

#### SINTAX

Alternatively, ASVs can be taxonomically classified using the SINTAX algorithm as implemented in VSEARCH. Depending on the reference taxonomy database (specified with `sintax_ref_taxonomy`), sequences can be classified down to species rank.

Files when _not_ using ITSx (default):

<details markdown="1">
<summary>Output files</summary>

- `sintax/`
  - `ASV_tax_sintax.*.raw.tsv`: Raw output from SINTAX.
  - `ASV_tax_sintax.*.tsv`: Taxonomic classification for each ASV sequence in a format similar to DADA2 output.
  - `ref_taxonomy.*.txt`: Information about the used reference taxonomy, such as title, version, citation.

</details>

Files when using ITSx:

<details markdown="1">
<summary>Output files</summary>

- `sintax/`
  - `ASV_ITS_tax_sintax.*.raw.tsv`: Raw output from SINTAX, when using cut sequences as input.
  - `ASV_tax_sintax.*.tsv`: Taxonomic classification for each ASV sequence, based on the chosen ITS region, in a format similar to DADA2 output.
  - `ref_taxonomy.*.txt`: Information about the used reference taxonomy, such as title, version, citation.

</details>

#### Kraken2

Kraken2 taxonomically classifies ASVs using exact k-mer matches. Kraken2 matches each k-mer within a query sequence to the lowest common ancestor (LCA) of all genomes/sequences containing the given k-mer.

<details markdown="1">
<summary>Output files</summary>

- `kraken2`
  - `ASV_tax.*.kraken2.report.txt`: Kraken2 report file, i.e. taxonomic classification of ASVs (shows number of ASVs matched at any given taxonomic level)
  - `ASV_tax.*.kraken2.keys.tsv`: Tab-separated table with extracted information from the report file
  - `ASV_tax.*.kraken2.classifiedreads.txt`: Classified sequence file, i.e. taxonomic classification (leaf) per ASV
  - `ASV_tax.*.kraken2.complete.tsv`: Tab-separated table with all extracted and parsed information from report and classified sequence file for each ASV
  - `ASV_tax.*.kraken2.tsv`: Tab-separated table with chosen taxonomic ranks per ASV
  - `ASV_tax.*.kraken2.into-qiime2.tsv`: Table with two tab-separated columns, `ASV_ID` and aggregated `taxonomy` (semicolon separated string), input to QIIME2

</details>

#### QIIME2

Taxonomic classification with QIIME2 is based on a classifier trained on sequences extracted with the primers.

<details markdown="1">
<summary>Output files</summary>

- `qiime2/taxonomy/`
  - `taxonomy.tsv`: Tab-separated table with taxonomic classification for each ASV
  - `*-classifier.qza`: QIIME2 artefact of the trained classifier. Can be supplied to other pipeline runs with `--classifier`
  - `ref_taxonomy.txt`: Information about the used reference taxonomy, such as title, version, citation.

</details>

### Phylogenetic placement and taxonomic classification

Phylogenetic placement grafts sequences onto a phylogenetic reference tree and optionally outputs taxonomic annotations. The reference tree is ideally made from full-length high-quality sequences containing better evolutionary signal than short amplicons. It is hence superior to estimating de-novo phylogenetic trees from short amplicon sequences. On providing required reference files, ASV sequences are aligned to the reference alignment with either [HMMER](http://hmmer.org/) (default) or [MAFFT](https://mafft.cbrc.jp/alignment/software/). Subsequently, phylogenetic placement of query sequences is performed with [EPA-NG](https://github.com/Pbdas/epa-ng), and finally a number of summary operations are performed with [Gappa](https://github.com/lczech/gappa). This uses code from [nf-core/phyloplace](https://nf-co.re/phyloplace) in the form of its main [subworkflow](https://github.com/nf-core/modules/tree/master/subworkflows/nf-core/fasta_newick_epang_gappa), therefore its detailed documentation also applies here.

<details markdown="1">
<summary>Output files</summary>

- `pplace/`
  - `*.graft.*.epa_result.newick`: Full phylogeny with query sequences grafted on to the reference phylogeny, in newick format.
  - `*.taxonomy.per_query.tsv`: Tab separated file with taxonomy information per query from classification by `gappa examine examinassign`
  - `*.per_query_unique.tsv`: Tab separated file with taxonomy information as above, but one row per query, by filtering for lowest LWR (likelihood weight ratio)
  - `*.heattree.tree.svg`: Heattree in SVG format from calling `gappa examine heattree`, see [Gappa documentation](https://github.com/Pbdas/epa-ng/blob/master/README.md) for details.
  - `pplace/hmmer/`: Contains intermediatary files if HMMER is used
  - `pplace/mafft/`: Contains intermediatary files if MAFFT is used
  - `pplace/epang/`: Output files from EPA-NG.
  - `pplace/gappa/`: Gappa output described in the [Gappa documentation](https://github.com/Pbdas/epa-ng/blob/master/README.md).

</details>

### Multiple region analysis with Sidle

Instead of relying on one short amplicon, scaffolding multiple regions along a reference can improve resolution over a single region. This method applies [Sidle (SMURF Implementation Done to acceLerate Efficiency)](https://doi.org/10.1101/2021.03.23.436606) within [QIIME2](https://pubmed.ncbi.nlm.nih.gov/31341288/).

Sidle reconstructs taxonomy profiles and abundances of several regions using a taxonomc database, therefore the previous sections about taxonomic classification are not applied.

Additional output includes reconstruction of abundance table and taxonomic information from multiple regions and a phylogenetic tree that can be further analysed. Apart from region-specific sequences, no useful de-novo sequences are generated.

<details markdown="1">
<summary>Output files</summary>

- `sidle/per-region/`
  - `ASV_seqs_region*_<FW_primer>_<RV_primer>.fasta`: ASV sequences per region
  - `ASV_table_region*_<FW_primer>_<RV_primer>.fasta`: ASV abundances per region
  - `DADA2_table_region*_<FW_primer>_<RV_primer>.fasta`: ASV abundances and sequences per region
- `sidle/DB/3_reconstructed/reconstruction_summary/index.html`: Information about the reconstructed reference taxonomy database
- `sidle/reconstructed/`
  - `reconstructed_feature-table.biom`: Unified abundance table in biom format
  - `reconstructed_feature-table.tsv`: Tab-separated unified abundance table
  - `reconstructed_taxonomy.tsv`: Tab-separated unified taxonomy table
  - `reconstructed_merged.tsv`: Tab-separated unified table with merged abundance and taxonomy information
  - `reconstructed_tree.nwk`: Phylogenetic tree
  - `reconstruction_table/index.html`: Information about the unified abundance table

More intermediate output is populated into the results folder when using `--save_intermediates`.

</details>

### Secondary analysis with QIIME2

**Quantitative Insights Into Microbial Ecology 2** ([QIIME2](https://qiime2.org/)) is a next-generation microbiome bioinformatics platform and the successor of the widely used [QIIME1](https://www.nature.com/articles/nmeth.f.303).

ASV sequences, counts, and taxonomic classification as produced before with DADA2 are imported into QIIME2 and further analysed. Optionally, ASVs can be taxonomically classified also with QIIME2 against a database chosen with `--qiime_ref_taxonomy` (but DADA2 taxonomic classification takes precedence). Next, ASVs are filtered (`--exclude_taxa`, `--min_frequency`, `--min_samples`), and abundance tables are exported. Following, diversity indices are calculated and testing for differential abundant features between sample groups is performed.

Intermediate data imported to QIIME2 is saved as QIIME2 fragments, that can be conveniently used for custom QIIME2 analysis.

<details markdown="1">
<summary>Output files</summary>

- `qiime2/input/`
  - `table.qza`: ASV count table.
  - `rep-seqs.qza`: ASV sequences.
  - `taxonomy.qza`: ASV taxonomic classification.

</details>

#### Abundance tables

The abundance tables are the final data for further downstream analysis and visualisations. The tables are based on the computed ASVs and taxonomic classification (in the following priotity: phylogenetic placement [EPA-NG, Gappa], DADA2, QIIME2), but after removal of unwanted taxa. Unwanted taxa are often off-targets generated in PCR with primers that are not perfectly specific for the target DNA (can be specified by `--exclude_taxa`), by default mitrochondria and chloroplast sequences are removed because these are frequent unwanted non-bacteria PCR products.

All following analysis is based on these filtered tables.

<details markdown="1">
<summary>Output files</summary>

- `qiime2/representative_sequences/`
  - `rep-seq.fasta`: Fasta file with ASV sequences.
  - `descriptive_stats.tsv`: Length, mean, etc. of ASV sequences.
  - `seven_number_summary.tsv`: Length of ASV sequences in different quantiles.
  - `filtered-sequences.qza`: QIIME2 fragment.
- `qiime2/abundance_tables/`
  - `abs-abund-table-*.tsv`: Tab-separated absolute abundance table at taxa level `*`, where `*` ranges by default from 2 to 6, specified by the `--tax_agglom_min` and `--tax_agglom_max` parameters.
  - `count_table_filter_stats.tsv`: Tab-separated table with information on how much counts were filtered for each sample.
  - `feature-table.biom`: Abundance table in biom format for importing into downstream analysis tools.
  - `feature-table.tsv`: Tab-separated abundance table for each ASV and each sample.
  - `filtered-table.qza`: QIIME2 fragment.

</details>

#### Relative abundance tables

Absolute abundance tables produced by the previous steps contain count data, but the compositional nature of 16S rRNA amplicon sequencing requires sequencing depth normalisation. This step computes relative abundance tables using TSS (Total Sum Scaling normalisation) for various taxonomic levels and detailed tables for all ASVs with taxonomic classification, sequence and relative abundance for each sample. Typically used for in depth investigation of taxa abundances. If not specified, the tables are based on the computed taxonomic classification (DADA2 classification takes precedence over QIIME2 classifications).

<details markdown="1">
<summary>Output files</summary>

- `qiime2/rel_abundance_tables/`
  - `rel-table-*.tsv`: Tab-separated absolute abundance table at taxa level `*`, where `*` ranges by default from 2 to 6, specified by the `--tax_agglom_min` and `--tax_agglom_max` parameters.
  - `rel-table-ASV.tsv`: Tab-separated relative abundance table for all ASVs.
  - `rel-table-ASV_with-DADA2-tax.tsv`: Tab-separated table for all ASVs with DADA2 taxonomic classification, sequence and relative abundance.
  - `rel-table-ASV_with-QIIME2-tax.tsv`: Tab-separated table for all ASVs with QIIME2 taxonomic classification, sequence and relative abundance.
  - `rel-table-ASV_with-PPLACE-tax.tsv`: Tab-separated table for all ASVs with EPA-NG - Gappa taxonomic classification, sequence and relative abundance.

</details>

#### Barplot

Produces an interactive abundance plot count tables that aids exploratory browsing the discovered taxa and their abundance in samples and allows sorting for associated meta data, DADA2 classification takes precedence over QIIME2 classifications. Optionally, barplots with average relative abundance values are produced for each metadata column chosen by `--metadata_category_barplot`.

<details markdown="1">
<summary>Output files</summary>

- `qiime2/barplot/`
  - `index.html`: Interactive barplot for taxa abundance per sample that can be viewed in your web browser.
- `qiime2/barplot_average/barplot_<category>`
  - `index.html`: Interactive barplot for averaged relative taxa abundance per group that can be viewed in your web browser.

</details>

#### Alpha diversity rarefaction curves

Produces rarefaction plots for several alpha diversity indices, and is primarily used to determine if the richness of the samples has been fully observed or sequenced. If the slope of the curves does not level out and the lines do not becomes horizontal, this might be because the sequencing depth was too low to observe all diversity or that sequencing error artificially increases sequence diversity and causes false discoveries.

<details markdown="1">
<summary>Output files</summary>

- `qiime2/alpha-rarefaction/`
  - `index.html`: Interactive alphararefaction curve for taxa abundance per sample that can be viewed in your web browser.

</details>

#### Diversity analysis

Diversity measures summarize important sample features (alpha diversity) or differences between samples (beta diversity). To do so, sample data is first rarefied to the minimum number of counts per sample. Parameter `--diversity_rarefaction_depth` can increase the rarefaction depth at the cost of excluding low count samples. Also, a phylogenetic tree of all ASVs is computed to provide phylogenetic information.

<details markdown="1">
<summary>Output files</summary>

- `qiime2/diversity/`
  - `Use the sampling depth of * for rarefaction.txt`: File that reports the rarefaction depth in the file name and file content.
- `qiime2/phylogenetic_tree/`
  - `tree.nwk`: Phylogenetic tree in newick format.
  - `rooted-tree.qza`: Phylogenetic tree in QIIME2 format.

</details>

##### Alpha diversity indices

Alpha diversity measures the species diversity within samples. Diversity calculations are based on sub-sampled data rarefied to the minimum read count of all samples. This step calculates alpha diversity using various methods and performs pairwise comparisons of groups of samples. It is based on a phylogenetic tree of all ASV sequences.

<details markdown="1">
<summary>Output files</summary>

- `qiime2/diversity/alpha_diversity/`
  - `evenness_vector/index.html`: Pielou’s Evenness.
  - `faith_pd_vector/index.html`: Faith’s Phylogenetic Diversity (qualitiative, phylogenetic).
  - `observed_otus_vector/index.html`: Observed OTUs (qualitative).
  - `shannon_vector/index.html`: Shannon’s diversity index (quantitative).

</details>

##### Beta diversity indices

Beta diversity measures the species community differences between samples. Diversity calculations are based on sub-sampled data rarefied to the minimum read count of all samples. This step calculates beta diversity distances using various methods and performs pairwise comparisons of groups of samples.
Additionally, principle coordinates analysis (PCoA) plots are produced that can be visualized with [Emperor](https://biocore.github.io/emperor/build/html/index.html) in your default browser without the need for installation. This calculations are based on a phylogenetic tree of all ASV sequences.
Furthermore, ADONIS permutation-based statistical test in vegan-R determine whether groups of samples are significantly different from one another. This test will only be executed when a custom formula is supplied with `--qiime_adonis_formula`, multiple formulas are comma separated. ADONIS computes an R2 value (effect size) which shows the percentage of variation explained by a condition, as well as a p-value to determine the statistical significance. The sequence of conditions in the formula matters, the variance of factors is removed (statistically controlled for) from beginning to end of the formula. E.g. formula "A+B" determines the effect size of "B" while controlling for "A".

**The following methods are used to calculate community dissimilarities:**

- Jaccard distance (qualitative)
- Bray-Curtis distance (quantitative)
- unweighted UniFrac distance (qualitative, phylogenetic)
- weighted UniFrac distance (quantitative, phylogenetic)

<details markdown="1">
<summary>Output files</summary>

- `qiime2/diversity/beta_diversity/`
  - `<method>_distance_matrix-<treatment>/index.html`: Box plots and significance analysis (PERMANOVA).
  - `<method>_pcoa_results-PCoA/index.html`: Interactive PCoA plot.
- `qiime2/diversity/beta_diversity/adonis/`
  - `<method>_distance_matrix-<adonis formula>/index.html`: Interactive (and .tsv) table of metadata feature importance and significance.
    - method: bray_curtis, jaccard, unweighted_unifrac, weighted_unifrac
    - treatment: depends on your metadata sheet or what metadata categories you have specified
    - adonis formula: Applied ADONIS formula

</details>

#### ANCOM

Analysis of Composition of Microbiomes ([ANCOM](https://www.ncbi.nlm.nih.gov/pubmed/26028277)) is applied to identify features that are differentially abundant across sample groups. A key assumption made by ANCOM is that few taxa (less than about 25%) will be differentially abundant between groups otherwise the method will be inaccurate. Parameter `--ancom_sample_min_count` sets the minimum sample counts to retain a sample for ANCOM analysis.

ANCOM is applied to each suitable or specified metadata column for 5 taxonomic levels (2-6).

<details markdown="1">
<summary>Output files</summary>

- `qiime2/ancom/`
  - `Category-<treatment>-<taxonomic level>/index.html`: Statistical results and interactive Volcano plot.
    - treatment: depends on your metadata sheet or what metadata categories you have specified
    - taxonomic level: level-2 (phylum), level-3 (class), level-4 (order), level-5 (family), level-6 (genus), ASV

</details>

### PICRUSt2

PICRUSt2 (Phylogenetic Investigation of Communities by Reconstruction of Unobserved States) is a software for predicting functional abundances based only on marker gene sequences. On demand (`--picrust`), Enzyme Classification numbers (EC), KEGG orthologs (KO) and MetaCyc ontology predictions will be made for each sample.

PICRUSt2 is preferentially applied to filtered data by QIIME2 but will use DADA2 output in case QIIME2 isnt run.

<details markdown="1">
<summary>Output files</summary>

- `PICRUSt2/`
  - `EC_pred_metagenome_unstrat_descrip.tsv`: Predicted quantifications for Enzyme Classification numbers (EC).
  - `KO_pred_metagenome_unstrat_descrip.tsv`: Predicted quantifications for KEGG orthologs (KO).
  - `METACYC_path_abun_unstrat_descrip.tsv`: Predicted quantifications for MetaCyc ontology.
  - `picrust.args.txt`: File containing arguments from the config file
- `PICRUSt2/all_output`

</details>

:::note
Quantifications are not normalized yet, they can be normalized e.g. by the total sum per sample.
:::

### SBDI export

You can use the `--sbdiexport` flag (or `sbdiexport: true` in a nextflow parameter file using `-params-file` in yml format) to generate tab separated files in preparation for submission to the [Swedish Biodiversity Infrastructure (SBDI)](https://biodiversitydata.se/).

Tables are generated from the DADA2 denoising and taxonomy assignment steps.
Each table, except `annotation.tsv`, corresponds to one tab in the [submission template](https://asv-portal.biodiversitydata.se/submit).
See [`docs/usage.md`](docs/usage.md) for further information.
Most of the fields in the template will not be populated by the export process, but if you run nf-core/ampliseq with a sample metadata table (`--metadata`) any fields corresponding to a field in the template will be used.

<details markdown="1">
<summary>Output files</summary>

- `SBDI/`
  - `annotation.tsv`: SBDI specific output for taxonomic reannotation, not used in submission to SBDI.
  - `asv-table.tsv`: asv-table tab of template.
  - `emof.tsv`: emof tab of template.
  - `event.tsv`: event tab of template.
  - `mixs.tsv`: mixs tab of template.

</details>

### Phyloseq

This directory will hold phyloseq objects for each taxonomy table produced by this pipeline. The objects will contain an ASV abundance table and a taxonomy table. If the pipeline is provided with metadata, that metadata will also be included in the phyloseq object. A phylogenetic tree will also be included if the pipeline produces a tree.

<details markdown="1">
<summary>Output files</summary>

- `phyloseq/`
  - `<taxonomy>_phyloseq.rds`: Phyloseq R object.

</details>

## Read count report

This report includes information on how many reads per sample passed each pipeline step in which a loss can occur. Specifically, how many read pairs entered cutadapt, were reverse complemented, passed trimming; how many read pairs entered DADA2, were denoised, merged and non-chimeric; and how many counts were lost during excluding unwanted taxa and removing low abundance/prevalence sequences in QIIME2.

<details markdown="1">
<summary>Output files</summary>

- `overall_summary.tsv`: Tab-separated file with count summary.

</details>

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>
