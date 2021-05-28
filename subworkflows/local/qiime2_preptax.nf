/*
 * Training of a classifier with QIIME2
 */

params.options = [:]

include { FORMAT_TAXONOMY_QIIME } from '../../modules/local/format_taxonomy_qiime'
include { QIIME2_EXTRACT        } from '../../modules/local/qiime2_extract' addParams( options: params.options )
include { QIIME2_TRAIN          } from '../../modules/local/qiime2_train'   addParams( options: params.options )

workflow QIIME2_PREPTAX {
    take:
	ch_dada_ref_taxonomy //channel, list of files
	FW_primer //val
	RV_primer //val
    
    main:
	FORMAT_TAXONOMY_QIIME ( ch_dada_ref_taxonomy )

	ch_ref_database = FORMAT_TAXONOMY_QIIME.out.fasta.combine(FORMAT_TAXONOMY_QIIME.out.tax)
	ch_ref_database
		.map {
			db ->
				def meta = [:]
				meta.FW_primer = FW_primer
				meta.RV_primer = RV_primer
				[ meta, db ] }
		.set { ch_ref_database }
	QIIME2_EXTRACT ( ch_ref_database )
	QIIME2_TRAIN ( QIIME2_EXTRACT.out.qza )

    emit:
    classifier      = QIIME2_TRAIN.out.qza
    version         = QIIME2_TRAIN.out.version
}