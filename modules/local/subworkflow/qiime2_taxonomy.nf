/*
 * Taxonomic classification with QIIME2
 */

params.options = [:]

include { QIIME2_INSEQ                  } from '../process/qiime2' addParams( options: params.options )
include { QIIME2_CLASSIFY               } from '../process/qiime2' addParams( options: params.options )

workflow QIIME2_TAXONOMY {
    take:
    ch_fasta
	ch_classifier
    
    main:	
	QIIME2_INSEQ ( ch_fasta )
	QIIME2_CLASSIFY ( ch_classifier, QIIME2_INSEQ.out.qza )

    emit:
    qza     = QIIME2_CLASSIFY.out.qza
	tsv     = QIIME2_CLASSIFY.out.tsv
    version = QIIME2_INSEQ.out.version
}