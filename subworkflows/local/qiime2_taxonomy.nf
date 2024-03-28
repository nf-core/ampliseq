/*
 * Taxonomic classification with QIIME2
 */

include { QIIME2_INSEQ                  } from '../../modules/local/qiime2_inseq'
include { QIIME2_CLASSIFY               } from '../../modules/local/qiime2_classify'

workflow QIIME2_TAXONOMY {
    take:
    ch_fasta
    ch_classifier

    main:
    QIIME2_INSEQ ( ch_fasta )
    QIIME2_CLASSIFY ( ch_classifier, QIIME2_INSEQ.out.qza )

    emit:
    qza      = QIIME2_CLASSIFY.out.qza
    tsv      = QIIME2_CLASSIFY.out.tsv
    versions = QIIME2_INSEQ.out.versions.mix(QIIME2_CLASSIFY.out.versions)
}
