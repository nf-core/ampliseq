# nf-core/ampliseq

## nf-core/ampliseq version 1.1.1

### Pipeline updates

* Update from QIIME2 v2018.6 to v2019.7, including DADA2 v1.6 to DADA2 v1.10

### Bug fixes

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
