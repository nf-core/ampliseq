/*
 * Training of a classifier with QIIME2
 */

params.options = [:]

include { QIIME2_EXTRACT } from '../process/qiime2' addParams( options: params.options )
include { QIIME2_TRAIN   } from '../process/qiime2' addParams( options: params.options )

workflow QIIME2_PREPTAX {
    take:
	ch_fasta_to_classifier //channel path
	ch_tax_to_classifier //channel path
	FW_primer //val
	RV_primer //val
    
    main:
	ch_ref_database = ch_fasta_to_classifier.combine(ch_tax_to_classifier)

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