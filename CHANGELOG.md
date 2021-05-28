# nf-core/ampliseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.0.0dev - [date]

Re-wrote whole pipeline in nextflow [DSL2](https://www.nextflow.io/docs/latest/dsl2.html) instead of DSL1

### `Added`

* [#229](https://github.com/nf-core/ampliseq/pull/229) - `--single_end` for single-ended Illumina data
* [#229](https://github.com/nf-core/ampliseq/pull/229), [#245](https://github.com/nf-core/ampliseq/pull/245), [#267](https://github.com/nf-core/ampliseq/pull/267)  - Taxonomic classification with DADA2
* [#229](https://github.com/nf-core/ampliseq/pull/229) - `--dada_ref_taxonomy` for taxonomic classification with DADA2's assignTaxonomy and addSpecies functions
* [#278](https://github.com/nf-core/ampliseq/pull/278) - `--qiime_ref_taxonomy` for taxonomic classification with QIIME2
* [#239](https://github.com/nf-core/ampliseq/pull/239) - Support of RDP database for DADA2 classification
* [#237](https://github.com/nf-core/ampliseq/pull/237) - Support of UNITE database for DADA2 classification
* [#229](https://github.com/nf-core/ampliseq/pull/229) - `--input` may point (1) at a fasta file ending with `.fasta`/`.fna`/`.fa` that will be taxonomically classified, (2) at a samples sheet ending with `.tsv` that allows analysis of multiple sequencing runs by reading the optional column `run`, or (3) at a folder input
* [#229](https://github.com/nf-core/ampliseq/pull/229) - `--sample_inference`, `--concatenate_reads`, `--illumina_pe_its`; please check the documentation for their function
* [#275](https://github.com/nf-core/ampliseq/pull/275) - Read count summary
* [#274](https://github.com/nf-core/ampliseq/pull/274) - `--skip_qiime` to prevent any steps that are executed with QIIME2
* [#272](https://github.com/nf-core/ampliseq/pull/272) - `--cut_its` to cut ASV sequence to ITS region before performing taxonomic classification with DADA2

### `Changed`

* [#254](https://github.com/nf-core/ampliseq/pull/254) - Updated CamelCase parameters to be lower_case_snake_case:
  * `multipleSequencingRuns` to `multiple_sequencing_runs`
  * `minLen` to `min_len`
  * `maxLen` to `max_len`
  * `maxEE` to `max_ee`
* [#277](https://github.com/nf-core/ampliseq/pull/277) - Requires nextflow version `>= 21.04.0`

### `Fixed`

* [#273](https://github.com/nf-core/ampliseq/pull/273) - Template update for nf-core/tools version 1.14

### `Dependencies`

* [#272](https://github.com/nf-core/ampliseq/pull/272) - New dependency ITSx v1.1.3
* [#229](https://github.com/nf-core/ampliseq/pull/229) - Updated from cutadapt v2.8 to v3.2
* [#229](https://github.com/nf-core/ampliseq/pull/229) - Updated DADA2 from v1.10 to v1.18.0, now not using QIIME2 for ASV generation any more
* [#229](https://github.com/nf-core/ampliseq/pull/229) - Updated QIIME2 to v2021.2

### `Removed`

* [#229](https://github.com/nf-core/ampliseq/pull/229) - `--manifest` is superseeded by `--input` that can now also handle a sample sheet file input (required extension: `.tsv`)
* [#229](https://github.com/nf-core/ampliseq/pull/229) - `--Q2imported` and `untilQ2import` are removed because pausing at that point is not neccessary
* [#229](https://github.com/nf-core/ampliseq/pull/229) - `--split` is no longer supported, therefore all sample IDs have to be unique
* [#229](https://github.com/nf-core/ampliseq/pull/229) - `--classifier_removeHash` and `--qiime_timezone` became unnecessary
* [#229](https://github.com/nf-core/ampliseq/pull/229) - `--onlyDenoising` is deprecated in favour of `--skip_taxonomy` (which does the exact same thing)
* `--taxon_reference` became unnecessary
* [#229](https://github.com/nf-core/ampliseq/pull/229) - `--reference_database` and `--dereplication` are not supported any more. `--qiime_ref_taxonomy` allows now choosing a taxonomic reference

## nf-core/ampliseq version 1.2.0 "Teal Bronze Lion" - 2021

### `Added`

* [#106](https://github.com/nf-core/ampliseq/issues/106) - Added support for PacBio data
* Added `--taxon_reference` to be able to support both 'silva' and 'unite'
* [#157](https://github.com/nf-core/ampliseq/issues/157) - Added possibility to run double cutadapt steps, `--double_primer`
* [#211](https://github.com/nf-core/ampliseq/issues/211) - Added quality filter option `--maxEE`

### `Fixed`

* [#182](https://github.com/nf-core/ampliseq/issues/182) - Fix input in case there are no underscores in sample IDs
* [#186](https://github.com/nf-core/ampliseq/issues/186) - Update github actions
* [#187](https://github.com/nf-core/ampliseq/issues/187) - Sample ids are in incorrect order in feature-table from PacBio data
* [#201](https://github.com/nf-core/ampliseq/pull/201) - Template update for nf-core/tools version 1.12.1
* [#147](https://github.com/nf-core/ampliseq/issues/147) - Split `make_classifier` in two different processes that can be allocated different resources
* [#183](https://github.com/nf-core/ampliseq/issues/183) - Don't fetch taxonomy/create classifier when run with `--skip_taxonomy`
* [#180](https://github.com/nf-core/ampliseq/issues/180) - MultiQC, cutadapt and fastQC now work with `--multipleSequencingRuns`

### `Dependencies`

* Updated from cutadapt v2.6 to v2.8

### `Deprecated`

## nf-core/ampliseq version 1.1.3 - 2020

### `Added`

* [#170](https://github.com/nf-core/ampliseq/issues/170) - Cite paper for initial release
* [#111](https://github.com/nf-core/ampliseq/issues/111) - Added parameter for user specified manifest file
* [#118](https://github.com/nf-core/ampliseq/issues/118) - Added social preview images
* [#135](https://github.com/nf-core/ampliseq/issues/135) - Added `--trunc_rmin` to make sure that auto trunc cutoff retaines a certain fraction of reads

### `Fixed`

* [#172](https://github.com/nf-core/ampliseq/pull/172) - Template update for nf-core/tools v1.11
* [#163](https://github.com/nf-core/ampliseq/pull/163) - Template update for nf-core/tools v1.10.2
* [#136](https://github.com/nf-core/ampliseq/issues/136) - Pipeline fails with remote working directory
* [#152](https://github.com/nf-core/ampliseq/issues/152) - Don't fetch taxonomy/create classifier when run with `--onlyDenoising`

### `Dependencies`

* Updated from MultiQC v1.6 to v1.9

### `Deprecated`

* `--reads` is replaced by `--input` due to nf-core/tools v1.10.2

## nf-core/ampliseq version 1.1.2 - 2019

* No further changes, except a bugfix for the [timezone](https://github.com/nf-core/ampliseq/issues/114) issue found by @marchoeppner
* Specification of `--qiime_timezone` might be required to run the analysis appropriately

## nf-core/ampliseq version 1.1.1 - 2019

### Pipeline Updates

* Update from QIIME2 v2018.6 to v2019.10, including DADA2 v1.6 to DADA2 v1.10
* Export absolute abundance files into 'results/abundance-table/filtered/' for optional external secondary analysis

### Bugfixes

* [#78](https://github.com/nf-core/ampliseq/issues/78) - All sequenced classifed to the same species

## nf-core/ampliseq version 1.1.0 "Silver Lime Bee" - 2019

### Pipeline updates

* [#40](https://github.com/nf-core/ampliseq/issues/40) - Added support for data originating from multiple sequencing runs
* [#53](https://github.com/nf-core/ampliseq/issues/53) - DADA2 report is always exported
* [#49](https://github.com/nf-core/ampliseq/issues/49) - Allowed more filtering options
* [#5](https://github.com/nf-core/ampliseq/issues/5) - Introduced check for existence of input files
* Extended parameter sanity check, including [#15](https://github.com/nf-core/ampliseq/issues/15)
* [#61](https://github.com/nf-core/ampliseq/issues/61) - Improved documentation
* [#62](https://github.com/nf-core/ampliseq/pull/62) - Utilize nf-core/configs centrally for this pipeline
* [#63](https://github.com/nf-core/ampliseq/issues/63) - QIIME imports files by using a manifest, giving more freedom with input file names
* [#84](https://github.com/nf-core/ampliseq/issues/84) - Add proper nf-core logo

### Bug fixes

* [#57](https://github.com/nf-core/ampliseq/issues/57) - Indicate exact regex for sequencing file names
* [#60](https://github.com/nf-core/ampliseq/issues/60) - publish demux.qza when --untilQ2import

## nf-core/ampliseq version 1.0.0 "Olive Steel Panda" - 2018-11-23

Initial release of nf-core/ampliseq, created with the [nf-core](http://nf-co.re/) template.
