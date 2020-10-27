# nf-core/ampliseq

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.2.0dev - [date]

### `Added`

* [#170](https://github.com/nf-core/ampliseq/issues/170) - Cite paper for initial release
* [#111](https://github.com/nf-core/ampliseq/issues/111) - Added parameter for user specified manifest file

### `Fixed`

* [#172](https://github.com/nf-core/ampliseq/pull/163) - Template update for nf-core/tools v1.11
* [#163](https://github.com/nf-core/ampliseq/pull/163) - Template update for nf-core/tools v1.10.2
* [#136](https://github.com/nf-core/ampliseq/issues/136) - Pipeline fails with remote working directory

### `Dependencies`

* Updated from MultiQC v1.6 to v1.9

### `Deprecated`

* --reads is replaced by --input due to nf-core/tools v1.10.2

## nf-core/ampliseq version 1.1.2 - 2019

* No further changes, except a bugfix for the [timezone](https://github.com/nf-core/ampliseq/issues/114) issue found by @marchoeppner
* Specification of '--qiime_timezone' might be required to run the analysis appropriately

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
