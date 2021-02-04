# nf-core/ampliseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
