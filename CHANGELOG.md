# nf-core/ampliseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## nf-core/ampliseq version 2.9.0 - 2024-04-03

### `Added`

- [#700](https://github.com/nf-core/ampliseq/pull/700) - Optional `--save_intermediates` to publish QIIME2 data objects (.qza) and visualisation objects (.qzv)
- [#702](https://github.com/nf-core/ampliseq/pull/702),[#723](https://github.com/nf-core/ampliseq/pull/723),[#728](https://github.com/nf-core/ampliseq/pull/728),[#729](https://github.com/nf-core/ampliseq/pull/729) - Add multiple regions analysis (including 5R / SMURF / q2-sidle)

### `Changed`

- [#719](https://github.com/nf-core/ampliseq/pull/719) - Versions of all (instead of selected) processes are now exported to `pipeline_info/software_versions.yml`

### `Fixed`

- [#697](https://github.com/nf-core/ampliseq/pull/697),[#699](https://github.com/nf-core/ampliseq/pull/699),[#713](https://github.com/nf-core/ampliseq/pull/713) - Template update for nf-core/tools version 2.13.1
- [#711](https://github.com/nf-core/ampliseq/pull/711) - From r207 and onwards Archaea sequences were omitted when parsing GTDB databases. (This did not affect `sbdi-gtdb` databases, only `gtdb`.)
- [#715](https://github.com/nf-core/ampliseq/pull/715) - Fix filtering vsearch clusters for high number of clusters
- [#717](https://github.com/nf-core/ampliseq/pull/717) - Fix edge case for sorting file names by using radix method
- [#718](https://github.com/nf-core/ampliseq/pull/718) - Require a minimum sequence length of 50bp for taxonomic classifcation after using ITSx
- [#721](https://github.com/nf-core/ampliseq/pull/721) - Fix error `unknown recognition error type: groovyjarjarantlr4.v4.runtime.LexerNoViableAltException` caused by a missing `\` in nf-core module `pigz/uncompress` (which had no consequences but was confusing)
- [#722](https://github.com/nf-core/ampliseq/pull/722) - When barrnap detects several genes select the lowest e-value
- [#726](https://github.com/nf-core/ampliseq/pull/726) - Add fallback to `download_pipeline.yml` because the pipeline does not support stub runs ([#2846](https://github.com/nf-core/tools/pull/2846))

### `Dependencies`

- [#720](https://github.com/nf-core/ampliseq/pull/720) - Updated nf-core modules, DADA2, and Phyloseq

| software | previously | now    |
| -------- | ---------- | ------ |
| cutadapt | 3.4        | 4.6    |
| DADA2    | 1.28.0     | 1.30.0 |
| Phyloseq | 1.44.0     | 1.46.0 |

### `Removed`

- [#710](https://github.com/nf-core/ampliseq/pull/710) - Removed Phyloref from DADA2 reference option because it's part of PR2 5.0.0

## nf-core/ampliseq version 2.8.0 - 2024-01-16

### `Added`

- [#666](https://github.com/nf-core/ampliseq/pull/666) - Added Greengenes2 database, version 2022.10, support for QIIME2 taxonomic classification.
- [#667](https://github.com/nf-core/ampliseq/pull/667),[#691](https://github.com/nf-core/ampliseq/pull/691) - Added `--qiime_ref_tax_custom` to permit custom reference database for QIIME2 taxonomic classification
- [#674](https://github.com/nf-core/ampliseq/pull/674) - Add PhytoRef database for DADA2 taxonomy assignment using `--dada_ref_taxonomy phytoref`
- [#675](https://github.com/nf-core/ampliseq/pull/675) - Add the Zehr lab nifH database for DADA2 taxonomy assignment using `--dada_ref_taxonomy zehr-nifh`
- [#681](https://github.com/nf-core/ampliseq/pull/681) - For DADA2, with `--dada_addspecies_allowmultiple` multiple exact species matches are reported and with `--dada_taxonomy_rc` reverse-complement matches are also considered in taxonomic classification

### `Changed`

- [#677](https://github.com/nf-core/ampliseq/pull/677) - Added cut_its information to SDBI export

### `Fixed`

- [#672](https://github.com/nf-core/ampliseq/pull/672),[#688](https://github.com/nf-core/ampliseq/pull/688),[#691](https://github.com/nf-core/ampliseq/pull/691) - Updated documentation
- [#676](https://github.com/nf-core/ampliseq/pull/676) - Phyloseq sometimes only produced one of multiple output files
- [#679](https://github.com/nf-core/ampliseq/pull/679) - Prevent masking low complexity regions by VSEARCH with lower case letters
- [#680](https://github.com/nf-core/ampliseq/pull/680),[#673](https://github.com/nf-core/ampliseq/pull/673) - Improved pipeline summary report & error messages
- [#683](https://github.com/nf-core/ampliseq/pull/683) - Template update for nf-core/tools version 2.11
- [#687](https://github.com/nf-core/ampliseq/pull/687) - Correct conda package for ASV SSU filtering

### `Dependencies`

### `Removed`

## nf-core/ampliseq version 2.7.1 - 2023-11-14

### `Added`

### `Changed`

- [#657](https://github.com/nf-core/ampliseq/pull/657) - Improved parameter descriptions and sequence

### `Fixed`

- [#655](https://github.com/nf-core/ampliseq/pull/655) - Added `NUMBA_CACHE_DIR` to fix downstream analysis with QIIME2 that failed on some systems
- [#656](https://github.com/nf-core/ampliseq/pull/656) - Moved conda-check to script-section and replaced `exit 1` with `error()`
- [#657](https://github.com/nf-core/ampliseq/pull/657) - Corrected inaccurate reporting of QIIME2 taxonomic classifications and ASV length filtering

### `Dependencies`

### `Removed`

## nf-core/ampliseq version 2.7.0 - 2023-10-20

### `Added`

- [#558](https://github.com/nf-core/ampliseq/pull/558),[#619](https://github.com/nf-core/ampliseq/pull/619),[#625](https://github.com/nf-core/ampliseq/pull/625),[#632](https://github.com/nf-core/ampliseq/pull/632),[#644](https://github.com/nf-core/ampliseq/pull/644) - Pipeline summary report
- [#615](https://github.com/nf-core/ampliseq/pull/615) - Phyloseq R object creation
- [#622](https://github.com/nf-core/ampliseq/pull/622) - ASV post-clustering with Vsearch
- [#637](https://github.com/nf-core/ampliseq/pull/637) - Taxonomic classification with Kraken2, parameter `--kraken2_ref_taxonomy`, `--kraken2_ref_tax_custom`, `--kraken2_assign_taxlevels`, `--kraken2_confidence`
- [#639](https://github.com/nf-core/ampliseq/pull/639) - GTDB release 214.1 for taxonomic classification with DADA2, using `--dada_ref_taxonomy gtdb` or `--dada_ref_taxonomy gtdb=R08-RS214`
- [#641](https://github.com/nf-core/ampliseq/pull/641) - Continue analysis even when individual files fail the filtering threshold, added parameter `--ignore_failed_filtering`

### `Changed`

- [#616](https://github.com/nf-core/ampliseq/pull/616) - When using a sample sheet with `--input` containing forward and reverse reads, specifying `--single_end` will only extract forward reads and treat the data as single ended instead of extracting forward and reverse reads.
- [#616](https://github.com/nf-core/ampliseq/pull/616) - `--input` was split into three params: (1) `--input` for samplesheet, (2) `--input_fasta` for ASV/OTU fasta input, (3) `--input_folder` direct FASTQ input

| Param updated | Param old | Accepts                                  |
| ------------- | --------- | ---------------------------------------- |
| input         | input     | samplesheet, .tsv/.csv/.yml/.yaml        |
| input_fasta   | input     | ASV/OTU sequences, .fasta                |
| input_folder  | input     | Folder containing compressed fastq files |

- [#639](https://github.com/nf-core/ampliseq/pull/639) - `--dada_ref_taxonomy gtdb` points towards GTDB release 214.1 instead of GTDB release 207 for taxonomic classification with DADA2
- [#645](https://github.com/nf-core/ampliseq/pull/645) - Updated documentation, including workflow figure

### `Fixed`

- [#605](https://github.com/nf-core/ampliseq/pull/605) - Make `--sbdiexport` compatible with PR2 version 5.0.0
- [#614](https://github.com/nf-core/ampliseq/pull/614),[#620](https://github.com/nf-core/ampliseq/pull/620),[#642](https://github.com/nf-core/ampliseq/pull/642) - Template update for nf-core/tools version 2.10
- [#617](https://github.com/nf-core/ampliseq/pull/617) - Fix database compatibility check for `--sbdiexport`
- [#628](https://github.com/nf-core/ampliseq/pull/628) - Fix edge case for sample sheet input when using specific combinations of sampleID and forwardReads or reverseReads that will forward one file too much to cutadapt
- [#630](https://github.com/nf-core/ampliseq/pull/630) - ASV rRNA (barrnap), length, and codon filter now work with ASV fasta file input
- [#633](https://github.com/nf-core/ampliseq/pull/633) - UNIFRAC in QIIME2_DIVERSITY_CORE is now prevented from using a GPU to avoid errors
- [#643](https://github.com/nf-core/ampliseq/pull/643) - Fix using `--skip_dada_addspecies` without `--dada_ref_tax_custom_sp` which was broken in 2.6.0 & 2.6.1
- [#647](https://github.com/nf-core/ampliseq/pull/647) - Update of credits

### `Dependencies`

- [#646](https://github.com/nf-core/ampliseq/pull/646) - Updated dependencies, see below:

| software | previously | now    |
| -------- | ---------- | ------ |
| FASTQC   | 0.11.9     | 0.12.1 |
| DADA2    | 1.22.0     | 1.28.0 |
| PICRUSt2 | 2.5.0      | 2.5.2  |
| QIIME2   | 2022.11    | 2023.7 |

### `Removed`

## nf-core/ampliseq version 2.6.1 - 2023-06-28

### `Fixed`

- [#603](https://github.com/nf-core/ampliseq/pull/603) - Fix all containers registry

## nf-core/ampliseq version 2.6.0 - 2023-06-27

### `Added`

- [#580](https://github.com/nf-core/ampliseq/pull/580) - Add NF-TEST pipeline end-to-end tests for existing CI tests
- [#591](https://github.com/nf-core/ampliseq/pull/591) - New version of the Unite taxonomy databases: 9.0
- [#596](https://github.com/nf-core/ampliseq/pull/596) - New version of the PR2 taxonomy database: 5.0.0, only available with DADA2 (`--dada_ref_taxonomy`)
- [#564](https://github.com/nf-core/ampliseq/pull/564),[#567](https://github.com/nf-core/ampliseq/pull/567),[#582](https://github.com/nf-core/ampliseq/pull/582) - Added phylogenetic placement
- [#577](https://github.com/nf-core/ampliseq/pull/577) - Added SINTAX for taxonomic classification
- [#575](https://github.com/nf-core/ampliseq/pull/575), [#586](https://github.com/nf-core/ampliseq/pull/586) - Added filtering step for stop codons for ASVs that are of coding regions.
- [#597](https://github.com/nf-core/ampliseq/pull/597) - Samples with less reads than specified with `--min_read_counts` (default: 1) stop the pipeline, previously the threshold was 1KB in size.

### `Changed`

- [#580](https://github.com/nf-core/ampliseq/pull/580) - GitHub Actions CI - pull_request to `dev` tests with NXF_VER `latest-everything` & pull_request to `master` tests with NXF_VER `22.10.1` & `latest-everything`
- [#563](https://github.com/nf-core/ampliseq/pull/563) - Renamed DADA2 taxonomic classification files to include the chosen reference taxonomy abbreviation.
- [#567](https://github.com/nf-core/ampliseq/pull/567) - Renamed `--dada_tax_agglom_min` and `--qiime_tax_agglom_min` to `--tax_agglom_min` and `--dada_tax_agglom_max` and `--qiime_tax_agglom_max` to `--tax_agglom_max`
- [#598](https://github.com/nf-core/ampliseq/pull/598) - Updated Workflow figure with SINTAX and phylogenetic placement
- [#599](https://github.com/nf-core/ampliseq/pull/599) - For exact species assignment (DADA2's addSpecies) PR2 taxonomy database (e.g. `--dada_ref_taxonomy pr2`) now excludes any taxa that end with " sp.".

### `Fixed`

- [#553](https://github.com/nf-core/ampliseq/pull/553) - Handle empty barrnap results files
- [#554](https://github.com/nf-core/ampliseq/pull/554) - Accept taxonomy strings that contain `#`,`'`
- [#569](https://github.com/nf-core/ampliseq/pull/569) - Make header of overall_summary.tsv consistent between input data types
- [#573](https://github.com/nf-core/ampliseq/pull/573) - Avoid parser error for single-end data when an empty read file is detected
- [#578](https://github.com/nf-core/ampliseq/pull/578) - Template update for nf-core/tools version 2.8, including changing `System.exit(1)` to `Nextflow.error()`
- [#594](https://github.com/nf-core/ampliseq/pull/594) - Update metadata documentation
- [#595](https://github.com/nf-core/ampliseq/pull/595) - Closing gaps in rarefaction depth for diversity calculations (`mindepth` in QIIME2_DIVERSITY_CORE)

### `Dependencies`

### `Removed`

## nf-core/ampliseq version 2.5.0 - 2023-03-02

### `Added`

- [#518](https://github.com/nf-core/ampliseq/pull/518),[#534](https://github.com/nf-core/ampliseq/pull/534) - Add COIDB DADA2 reference taxonomy database
- [#521](https://github.com/nf-core/ampliseq/pull/521) - Export svg in addition to pdf files for quality plots from DADA2
- [#538](https://github.com/nf-core/ampliseq/pull/538) - Parameter `--diversity_rarefaction_depth` controls the minimum rarefaction depth for diversity analysis, this allows increasing the rarefaction depth at the cost of excluding low count samples. Parameter `--ancom_sample_min_count` sets the minimum sample counts to retain a sample for ANCOM analysis.

### `Changed`

- [#537](https://github.com/nf-core/ampliseq/pull/537) - Update output generated with option sbdi-export
- [#541](https://github.com/nf-core/ampliseq/pull/541) - Remove adjustments of taxonomic levels for RDP & SILVA & GTDB & UNITE database for DADA2 taxonomic classification, reduced default of `--dada_tax_agglom_max` from 7 to 6
- [#548](https://github.com/nf-core/ampliseq/pull/548) - `--filter_ssu` accepted any barrnap hit to a kingdom (domain) (any occurence in resulting gff) to choose an ASV, now only ASVs with the kingdom (domain) that has the lowest evalue are accepted.

### `Fixed`

- [#513](https://github.com/nf-core/ampliseq/pull/513) - Template update for nf-core/tools version 2.7.2
- [#519](https://github.com/nf-core/ampliseq/pull/519) - Adding the pipeline reference to the MultiQC report
- [#520](https://github.com/nf-core/ampliseq/pull/520),[#530](https://github.com/nf-core/ampliseq/pull/530) - Fix conda packages
- [#531](https://github.com/nf-core/ampliseq/pull/531),[#546](https://github.com/nf-core/ampliseq/pull/546) - Update documentation
- [#535](https://github.com/nf-core/ampliseq/pull/535) - Make sure barrnap runs with fasta input
- [#544](https://github.com/nf-core/ampliseq/pull/544) - Adding module to fix header in fasta input if needed

### `Dependencies`

- [#528](https://github.com/nf-core/ampliseq/pull/528) - Updated QIIME2

| Tool   | Previous version | New version |
| ------ | ---------------- | ----------- |
| QIIME2 | 2022.8           | 2022.11     |

### `Removed`

- [#513](https://github.com/nf-core/ampliseq/pull/513) - Removed parameter `--enable_conda`.

## nf-core/ampliseq version 2.4.1 - 2022-12-07

### `Added`

- [#494](https://github.com/nf-core/ampliseq/pull/494) - `--metadata_category_barplot` accepts a comma separated list of metadata categories and plots for each barplots with average relative abundance.

### `Changed`

- [#492](https://github.com/nf-core/ampliseq/pull/492) - `--qiime_adonis_formula` accepts a comma separated list of formulas.

### `Fixed`

- [#486](https://github.com/nf-core/ampliseq/pull/486) - Fixed typo in error message stating `--skip_classifer` instead of `--classifier`.
- [#487](https://github.com/nf-core/ampliseq/pull/487),[#488](https://github.com/nf-core/ampliseq/pull/488) - Update stale links in usage documentation.
- [#489](https://github.com/nf-core/ampliseq/pull/489) - Reduce linting warnings for nf-core tools version 2.5.1.
- [#491](https://github.com/nf-core/ampliseq/pull/491) - Make output from --addSH match UNITE format by replacing spaces with underscores.
- [#495](https://github.com/nf-core/ampliseq/pull/495) - Template update for nf-core/tools version 2.6
- [#501](https://github.com/nf-core/ampliseq/pull/501) - Check for empty fields in samplesheet column "run" and raise an appropriate error.
- [#503](https://github.com/nf-core/ampliseq/pull/503) - Changed environment for formatting databases.
- [#504](https://github.com/nf-core/ampliseq/pull/504) - Fixed warnings with nextflow 22.10 (and later) about processes that are defined more than once.

### `Dependencies`

| Tool   | Previous version | New version |
| ------ | ---------------- | ----------- |
| QIIME2 | 2021.8           | 2022.8      |

### `Removed`

## nf-core/ampliseq version 2.4.0 - 2022-09-07

### `Added`

- [#456](https://github.com/nf-core/ampliseq/pull/456) - An optional ASV length filter can be activated using `--min_len_asv <int>` and/or `--max_len_asv <int>`.
- [#458](https://github.com/nf-core/ampliseq/pull/458) - Samplesheet, ASV fasta file, and/or metadata sheet is now exported into `<results>/input/`
- [#459](https://github.com/nf-core/ampliseq/pull/459) - MIDORI2 CO1 database with keys `midori2-co1=gb250` and `midori2-co1` for `--dada_ref_taxonomy`
- [#460](https://github.com/nf-core/ampliseq/pull/460) - Taxonomic ranks for DADA2 taxonomic classification can be now adjusted using `--dada_assign_taxlevels <comma separated string>`.
- [#461](https://github.com/nf-core/ampliseq/pull/461) - A custom DADA2 reference taxonomy database can now be used with `--dada_ref_tax_custom` and `--dada_ref_tax_custom_sp`, typically accompanied by `--dada_assign_taxlevels`.
- [#446](https://github.com/nf-core/ampliseq/pull/446),[#467](https://github.com/nf-core/ampliseq/pull/467) - Binned quality scores from Illumina NovaSeq data can be now corrected with `--illumina_novaseq`.
- [#477](https://github.com/nf-core/ampliseq/pull/477) - QC plots of DADA2's plotQualityProfile are now also produced after preprocessing.
- [#478](https://github.com/nf-core/ampliseq/pull/478) - Added GTDB R07-RS207 DADA2 taxonomy reference databases

### `Changed`

- [#444](https://github.com/nf-core/ampliseq/pull/444),[#457](https://github.com/nf-core/ampliseq/pull/457),[#463](https://github.com/nf-core/ampliseq/pull/463),[#465](https://github.com/nf-core/ampliseq/pull/465),[#466](https://github.com/nf-core/ampliseq/pull/466),[#469](https://github.com/nf-core/ampliseq/pull/469) - Updated the documentation.
- [#445](https://github.com/nf-core/ampliseq/pull/445) - The minimum number of total bases to use for error rate learning by default is 1e8 (DADA2, learnErrors, nbases). Previously, samples were read in the provided order until enough reads were obtained (DADA2, learnErrors, randomize=FALSE). Now, samples are picked at random from those provided (DADA2, learnError, randomize=TRUE) and a seed is set.
- [#453](https://github.com/nf-core/ampliseq/pull/453) - Export a few more basic QIIME2 fragments (zipped files) that can be easily imported into the correct QIIME2 version for custom analysis.
- [#464](https://github.com/nf-core/ampliseq/pull/464) - Reported taxonomic classifications on species level based on DADA2's assignTaxonomy (approximations) is now listed in column "Species" while exact matches based on DADA2's addSpecies are now reported in column "Species_exact".

### `Fixed`

- [#448](https://github.com/nf-core/ampliseq/pull/448) - Updated SBDI export scripts to include Unite species hypothesis information if available.
- [#451](https://github.com/nf-core/ampliseq/pull/451) - Pairwise statistics will be now performed on a subset of metadata columns specified with `--metadata_category` instead of ignoring that setting.
- [#451](https://github.com/nf-core/ampliseq/pull/451) - Replace busybox with Ubuntu base image for GCP support.
- [#455](https://github.com/nf-core/ampliseq/pull/455) - Stop with descriptive error when only one of `--trunclenf` and `--trunclenr` is given, earlier it was silently ignored.
- [#474](https://github.com/nf-core/ampliseq/pull/474) - Template update for nf-core/tools version 2.5.1
- [#475](https://github.com/nf-core/ampliseq/pull/475) - Report software versions for DADA2_TAXONOMY

### `Dependencies`

- [#479](https://github.com/nf-core/ampliseq/pull/479) - Updated software

| Tool     | Previous version | New version |
| -------- | ---------------- | ----------- |
| PICRUSt2 | 2.4.2            | 2.5.0       |
| MultiQC  | 1.12             | 1.13a       |

### `Removed`

## nf-core/ampliseq version 2.3.2 - 2022-05-27

### `Added`

- [#429](https://github.com/nf-core/ampliseq/pull/429) - `--cutadapt_min_overlap` sets cutadapt's global minimum overlap (`-O`) and `--cutadapt_max_error_rate` sets cutadapt's global maximum error rate (`-e`) for trimming primer sequences.
- [#431](https://github.com/nf-core/ampliseq/pull/431) - `--skip_dada_quality` allows to skip quality check with DADA2. This is only allowed when `--trunclenf` and `--trunclenr` are set.
- [#434](https://github.com/nf-core/ampliseq/pull/434) - `--addsh` adds UNITE species hypothesis (SH) to the taxonomy. Only available for UNITE databases.

### `Changed`

- [#432](https://github.com/nf-core/ampliseq/pull/432) - The number of records to sample from a fastq file was decreased from 5e+06 to 5e+04 for plotQualityProfile (DADA2_QUALITY), therefore a smaller subset of reads is sampled for determining `--trunlenf` and `--trunclenr`. This should make the process more robust also from larger data sets.

### `Fixed`

- [#428](https://github.com/nf-core/ampliseq/pull/428) - Fixed samplesheet sampleID entries, now allows dashes.
- [#433](https://github.com/nf-core/ampliseq/pull/433) - Fixed typos and improved documentation layout.
- [#437](https://github.com/nf-core/ampliseq/pull/437) - Template update for nf-core/tools version 2.4
- [#439](https://github.com/nf-core/ampliseq/pull/439) - Fixed a bug in DADA2_QUALITY process with large number of nucleotides

### `Dependencies`

### `Removed`

## nf-core/ampliseq version 2.3.1 - 2022-04-05

### `Added`

### `Changed`

### `Fixed`

- [#415](https://github.com/nf-core/ampliseq/pull/415) - ADONIS test was running by default in version 2.2.0 and 2.3.0 on metadata columns that did not always have the required format, now specifying `--qiime_adonis_formula` is required to run that step.

### `Dependencies`

### `Removed`

## nf-core/ampliseq version 2.3.0 - 2022-04-04

### `Added`

- [#385](https://github.com/nf-core/ampliseq/pull/385) - `--skip_cutadapt` allows to skip primer trimmimg.
- [#390](https://github.com/nf-core/ampliseq/pull/390), [#408](https://github.com/nf-core/ampliseq/pull/408) - Add values to option `--cut_its`, specifying which ITS part to use. Also, add option `--its_partial <x>` to allow partial ITS sequences longer than given cutoff.
- [#395](https://github.com/nf-core/ampliseq/pull/395) - `--seed` specifies the random seed.
- [#396](https://github.com/nf-core/ampliseq/pull/396) - Barrnap annotates ASV sequences for SSU's, it can be skipped with `--skip_barrnap`. `--filter_ssu` takes a comma separated list of "bac,arc,mito,euk" and enables SSU filtering depending on Barrnap (default: off).
- [#397](https://github.com/nf-core/ampliseq/pull/397) - Complement README.md with links to the nf-core bytesize 25 (nf-core/ampliseq).

### `Changed`

- [#385](https://github.com/nf-core/ampliseq/pull/385) - `--FW_primer` and `--RV_primer` are not obligatory any more, however primer sequences are still required with cutadapt (i.e. without `--skip_cutadapt`), `--qiime_ref_taxonomy`, and `--cut_dada_ref_taxonomy` (cuts reference with primer sequences).

### `Fixed`

- [#384](https://github.com/nf-core/ampliseq/pull/384) - For QIIME2 beta diversity, make directory before execution.
- [#394](https://github.com/nf-core/ampliseq/pull/394) - Prevent simultaneous usage of `--qiime_ref_taxonomy` and `--classifier`.
- [#402](https://github.com/nf-core/ampliseq/pull/402), [#410](https://github.com/nf-core/ampliseq/pull/410) - Template update for nf-core/tools version 2.3.2
- [#403](https://github.com/nf-core/ampliseq/pull/403) - Limit number of files for DADA2_QUALITY (plotQualityProfile) by read numbers.
- [#405](https://github.com/nf-core/ampliseq/pull/405) - Fix omitting taxonomic filtering (QIIME2_FILTERTAXA).
- [#404](https://github.com/nf-core/ampliseq/pull/404) - Rephrase error message for empty input files and empty files after trimming with cutadapt.
- [#407](https://github.com/nf-core/ampliseq/pull/407) - Fix reformatting script for Unite.

### `Dependencies`

| Tool    | Previous version | New version |
| ------- | ---------------- | ----------- |
| MultiQC | 1.11             | 1.12        |

### `Removed`

## nf-core/ampliseq version 2.2.0 - 2022-01-31

### `Added`

- [#352](https://github.com/nf-core/ampliseq/pull/352), [#372](https://github.com/nf-core/ampliseq/pull/372) - `--skip_dada_addspecies` allows to skip species level classification to reduce memory requirements, incompatible with `--sbdiexport` that expect species annotation.
- [#364](https://github.com/nf-core/ampliseq/pull/364) - Adonis in QIIME2 for testing feature importance in beta diversity distances, `--qiime_adonis_formula` can be set to provide a custom formula.
- [#366](https://github.com/nf-core/ampliseq/pull/366) - New version of the SBDI-GTDB taxonomy database: v. 3. (Fixes problem with `Reverse_` added to some domain strings.)

### `Changed`

- [#354](https://github.com/nf-core/ampliseq/pull/354) - Input files and files after primer trimming with cutadapt are required to be >1KB (i.e. not empty) and either the pipeline will stop if at least one sample file fails or the failing samples will be ignored when using `--ignore_empty_input_files` or `--ignore_failed_trimming`, respectively.
- [#376](https://github.com/nf-core/ampliseq/pull/376) - Forbid sampleIDs starting with a number when also `--metadata` is used, because such strings are unintentionally modified and the metadata will not match any more.

### `Fixed`

- [#377](https://github.com/nf-core/ampliseq/pull/377)- An error message will occur when `--sbdiexport` is used with `--skip_taxonomy` or `--skip_dada_addspecies`
- [#375](https://github.com/nf-core/ampliseq/pull/375)- Updated documentation regarding not using curly brackets in `--extension` with `--single_end`
- [#362](https://github.com/nf-core/ampliseq/pull/362)- Template update for nf-core/tools version 2.2, now requires nextflow version `>= 21.10.3`
- [#374](https://github.com/nf-core/ampliseq/pull/374)- Cutadapt results can be now also viwed in the MultiQC report

### `Dependencies`

| Tool     | Previous version | New version |
| -------- | ---------------- | ----------- |
| Cutadapt | 3.2              | 3.4         |
| DADA2    | 1.20.0           | 1.22.0      |
| QIIME2   | 2021.2           | 2021.8      |
| PICRUSt2 | 2.4.1            | 2.4.2       |
| MultiQC  | 1.10.1           | 1.11        |

### `Removed`

- [#350](https://github.com/nf-core/ampliseq/pull/350) - Remove redundant derepFastq step (has no impact on results)

## nf-core/ampliseq version 2.1.1 - 2021-10-28

### `Added`

- [#336](https://github.com/nf-core/ampliseq/pull/336) - Taxa agglomeration levels with `--dada_tax_agglom_min`, `--dada_tax_agglom_max`, `--qiime_tax_agglom_min`, `--qiime_tax_agglom_max`, with defaults that go to genus level for abundance tables and ANCOM analysis

### `Changed`

- [338](https://github.com/nf-core/ampliseq/pull/338) - Write empty space instead of `NA` for missing values in output files.
- [342](https://github.com/nf-core/ampliseq/pull/342) - Added PICRUSt2 to summary figure.

### `Fixed`

- [#329](https://github.com/nf-core/ampliseq/issues/329) - Improve error message when no data files are found
- [#330](https://github.com/nf-core/ampliseq/issues/330) - Make `--skip_fastqc` usable again
- [#339](https://github.com/nf-core/ampliseq/issues/339) - Fix sample names when using `--double_primer` or `--illumina_pe_its`

### `Dependencies`

### `Removed`

## nf-core/ampliseq version 2.1.0 "Gray Steel Boa" - 2021-09-14

### `Added`

- [#322](https://github.com/nf-core/ampliseq/pull/322) - Export to the Swedish Biodiversity Infrastructure (SBDI)
- [#294](https://github.com/nf-core/ampliseq/pull/294) - New version of the PR2 taxonomy database: 4.14.0, contains also some Bacteria, see [PR2 release 4.14.0 notes](https://github.com/pr2database/pr2database/releases/tag/v4.14.0)
- [#294](https://github.com/nf-core/ampliseq/pull/294) - New version of the Unite taxonomy databases: 8.3
- [#302](https://github.com/nf-core/ampliseq/pull/302) - Pipeline workflow figure in README.md
- [#307](https://github.com/nf-core/ampliseq/pull/307) - Functional predictions with PICRUSt2, on demand with `--picrust`
- [#310](https://github.com/nf-core/ampliseq/pull/310) - Workflow figure in usage.md when using `--multiple_sequencing_runs`
- [#312](https://github.com/nf-core/ampliseq/pull/312) - Added curated GTDB 16S taxonomy: `sbdi-gtdb` as parameter to `--dada_ref_taxonomy`
- [#318](https://github.com/nf-core/ampliseq/pull/318) - Output information about the used reference taxonomy in a separate file in results folder `dada2/` or `qiime2/taxonomy`

### `Changed`

- [#313](https://github.com/nf-core/ampliseq/pull/313) - Relative abundance tables in `qiime2/rel_abundance_tables/` on ASV level were renamed and with either DADA2 (`rel-table-ASV_with-DADA2-tax.tsv`) or QIIME2 classifications (`rel-table-ASV_with-QIIME2-tax.tsv`), if available.

### `Fixed`

- [#306](https://github.com/nf-core/ampliseq/pull/306) - Sample names can now be identical to basenames of read files
- [#299](https://github.com/nf-core/ampliseq/pull/299), [#301](https://github.com/nf-core/ampliseq/pull/301)- Template update for nf-core/tools version 2.1
- [#303](https://github.com/nf-core/ampliseq/pull/303) - Reverse primer of PacBio and IonTorrent reads should be now given to `--RV_primer` in usual direction (before: reverse complement)
- [#305](https://github.com/nf-core/ampliseq/pull/305) - `--max_len` now accepts integers as expected
- [#314](https://github.com/nf-core/ampliseq/pull/314), [#315](https://github.com/nf-core/ampliseq/pull/315) - ASV fasta input via --input fixed (was broken in 2.0.0) and a test profile was added.

### `Dependencies`

- [#299](https://github.com/nf-core/ampliseq/pull/299) - Updated MultiQC to v1.10
- [#319](https://github.com/nf-core/ampliseq/pull/319) - Updated DADA2 from 1.18.0 to 1.20.0

### `Removed`

## nf-core/ampliseq version 2.0.0 "Blue Copper Kangaroo" - 2021-06-29

Re-wrote whole pipeline in nextflow [DSL2](https://www.nextflow.io/docs/latest/dsl2.html) instead of DSL1

### `Added`

- [#229](https://github.com/nf-core/ampliseq/pull/229) - `--single_end` for single-ended Illumina data
- [#229](https://github.com/nf-core/ampliseq/pull/229), [#245](https://github.com/nf-core/ampliseq/pull/245), [#267](https://github.com/nf-core/ampliseq/pull/267) - Taxonomic classification with DADA2
- [#229](https://github.com/nf-core/ampliseq/pull/229) - `--dada_ref_taxonomy` for taxonomic classification with DADA2's assignTaxonomy and addSpecies functions
- [#278](https://github.com/nf-core/ampliseq/pull/278) - `--qiime_ref_taxonomy` for taxonomic classification with QIIME2
- [#239](https://github.com/nf-core/ampliseq/pull/239) - Support of RDP database for DADA2 classification
- [#237](https://github.com/nf-core/ampliseq/pull/237) - Support of UNITE database for DADA2 classification
- [#229](https://github.com/nf-core/ampliseq/pull/229) - `--input` may point (1) at a fasta file ending with `.fasta`/`.fna`/`.fa` that will be taxonomically classified, (2) at a samples sheet ending with `.tsv` that allows analysis of multiple sequencing runs by reading the optional column `run`, or (3) at a folder input
- [#229](https://github.com/nf-core/ampliseq/pull/229) - `--sample_inference`, `--concatenate_reads`, `--illumina_pe_its`; please check the documentation for their function
- [#275](https://github.com/nf-core/ampliseq/pull/275) - Read count summary
- [#274](https://github.com/nf-core/ampliseq/pull/274) - `--skip_qiime` to prevent any steps that are executed with QIIME2
- [#272](https://github.com/nf-core/ampliseq/pull/272) - `--cut_its` to cut ASV sequence to ITS region before performing taxonomic classification with DADA2
- [#280](https://github.com/nf-core/ampliseq/pull/280) - Added support for IonTorrent data
- [#283](https://github.com/nf-core/ampliseq/pull/283) - `--cut_dada_ref_taxonomy` allows extracting expected amplicons from DADA2 reference taxonomy database

### `Changed`

- [#254](https://github.com/nf-core/ampliseq/pull/254) - Updated CamelCase parameters to be lower_case_snake_case:
  - `multipleSequencingRuns` to `multiple_sequencing_runs`
  - `minLen` to `min_len`
  - `maxLen` to `max_len`
  - `maxEE` to `max_ee`
- [#277](https://github.com/nf-core/ampliseq/pull/277) - Requires nextflow version `>= 21.04.0`

### `Fixed`

- [#273](https://github.com/nf-core/ampliseq/pull/273) - Template update for nf-core/tools version 1.14

### `Dependencies`

- [#272](https://github.com/nf-core/ampliseq/pull/272) - New dependency ITSx v1.1.3
- [#229](https://github.com/nf-core/ampliseq/pull/229) - Updated from cutadapt v2.8 to v3.2
- [#229](https://github.com/nf-core/ampliseq/pull/229) - Updated DADA2 from v1.10 to v1.18.0, now not using QIIME2 for ASV generation any more
- [#229](https://github.com/nf-core/ampliseq/pull/229) - Updated QIIME2 to v2021.2

### `Removed`

- [#229](https://github.com/nf-core/ampliseq/pull/229) - `--manifest` is superseeded by `--input` that can now also handle a sample sheet file input (required extension: `.tsv`)
- [#229](https://github.com/nf-core/ampliseq/pull/229) - `--Q2imported` and `untilQ2import` are removed because pausing at that point is not neccessary
- [#229](https://github.com/nf-core/ampliseq/pull/229) - `--split` is no longer supported, therefore all sample IDs have to be unique
- [#229](https://github.com/nf-core/ampliseq/pull/229) - `--classifier_removeHash` and `--qiime_timezone` became unnecessary
- [#229](https://github.com/nf-core/ampliseq/pull/229) - `--onlyDenoising` is deprecated in favour of `--skip_taxonomy` (which does the exact same thing)
- `--taxon_reference` became unnecessary
- [#229](https://github.com/nf-core/ampliseq/pull/229) - `--reference_database` and `--dereplication` are not supported any more. `--qiime_ref_taxonomy` allows now choosing a taxonomic reference

## nf-core/ampliseq version 1.2.0 "Teal Bronze Lion" - 2021-02-04

### `Added`

- [#106](https://github.com/nf-core/ampliseq/issues/106) - Added support for PacBio data
- Added `--taxon_reference` to be able to support both 'silva' and 'unite'
- [#157](https://github.com/nf-core/ampliseq/issues/157) - Added possibility to run double cutadapt steps, `--double_primer`
- [#211](https://github.com/nf-core/ampliseq/issues/211) - Added quality filter option `--maxEE`

### `Fixed`

- [#182](https://github.com/nf-core/ampliseq/issues/182) - Fix input in case there are no underscores in sample IDs
- [#186](https://github.com/nf-core/ampliseq/issues/186) - Update github actions
- [#187](https://github.com/nf-core/ampliseq/issues/187) - Sample ids are in incorrect order in feature-table from PacBio data
- [#201](https://github.com/nf-core/ampliseq/pull/201) - Template update for nf-core/tools version 1.12.1
- [#147](https://github.com/nf-core/ampliseq/issues/147) - Split `make_classifier` in two different processes that can be allocated different resources
- [#183](https://github.com/nf-core/ampliseq/issues/183) - Don't fetch taxonomy/create classifier when run with `--skip_taxonomy`
- [#180](https://github.com/nf-core/ampliseq/issues/180) - MultiQC, cutadapt and fastQC now work with `--multipleSequencingRuns`

### `Dependencies`

- Updated from cutadapt v2.6 to v2.8

### `Deprecated`

## nf-core/ampliseq version 1.1.3 - 2020-11-02

### `Added`

- [#170](https://github.com/nf-core/ampliseq/issues/170) - Cite paper for initial release
- [#111](https://github.com/nf-core/ampliseq/issues/111) - Added parameter for user specified manifest file
- [#118](https://github.com/nf-core/ampliseq/issues/118) - Added social preview images
- [#135](https://github.com/nf-core/ampliseq/issues/135) - Added `--trunc_rmin` to make sure that auto trunc cutoff retaines a certain fraction of reads

### `Fixed`

- [#172](https://github.com/nf-core/ampliseq/pull/172) - Template update for nf-core/tools v1.11
- [#163](https://github.com/nf-core/ampliseq/pull/163) - Template update for nf-core/tools v1.10.2
- [#136](https://github.com/nf-core/ampliseq/issues/136) - Pipeline fails with remote working directory
- [#152](https://github.com/nf-core/ampliseq/issues/152) - Don't fetch taxonomy/create classifier when run with `--onlyDenoising`

### `Dependencies`

- Updated from MultiQC v1.6 to v1.9

### `Deprecated`

- `--reads` is replaced by `--input` due to nf-core/tools v1.10.2

## nf-core/ampliseq version 1.1.2 - 2019-12-19

- No further changes, except a bugfix for the [timezone](https://github.com/nf-core/ampliseq/issues/114) issue found by @marchoeppner
- Specification of `--qiime_timezone` might be required to run the analysis appropriately

## nf-core/ampliseq version 1.1.1 - 2019-12-09

### Pipeline Updates

- Update from QIIME2 v2018.6 to v2019.10, including DADA2 v1.6 to DADA2 v1.10
- Export absolute abundance files into 'results/abundance-table/filtered/' for optional external secondary analysis

### Bugfixes

- [#78](https://github.com/nf-core/ampliseq/issues/78) - All sequenced classifed to the same species

## nf-core/ampliseq version 1.1.0 "Silver Lime Bee" - 2019-07-15

### Pipeline updates

- [#40](https://github.com/nf-core/ampliseq/issues/40) - Added support for data originating from multiple sequencing runs
- [#53](https://github.com/nf-core/ampliseq/issues/53) - DADA2 report is always exported
- [#49](https://github.com/nf-core/ampliseq/issues/49) - Allowed more filtering options
- [#5](https://github.com/nf-core/ampliseq/issues/5) - Introduced check for existence of input files
- Extended parameter sanity check, including [#15](https://github.com/nf-core/ampliseq/issues/15)
- [#61](https://github.com/nf-core/ampliseq/issues/61) - Improved documentation
- [#62](https://github.com/nf-core/ampliseq/pull/62) - Utilize nf-core/configs centrally for this pipeline
- [#63](https://github.com/nf-core/ampliseq/issues/63) - QIIME imports files by using a manifest, giving more freedom with input file names
- [#84](https://github.com/nf-core/ampliseq/issues/84) - Add proper nf-core logo

### Bug fixes

- [#57](https://github.com/nf-core/ampliseq/issues/57) - Indicate exact regex for sequencing file names
- [#60](https://github.com/nf-core/ampliseq/issues/60) - publish demux.qza when --untilQ2import

## nf-core/ampliseq version 1.0.0 "Olive Steel Panda" - 2018-11-23

Initial release of nf-core/ampliseq, created with the [nf-core](http://nf-co.re/) template.
