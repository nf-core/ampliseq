/*
 * Make average barplot for metadata categories specified with params.metadata_category_barplot
 */

include { QIIME2_INASV as QIIME2_INASV_BPAVG } from '../../modules/local/qiime2_inasv'
include { QIIME2_FEATURETABLE_GROUP          } from '../../modules/local/qiime2_featuretable_group'
include { QIIME2_BARPLOT as QIIME2_BPAVG     } from '../../modules/local/qiime2_barplot'

workflow QIIME2_BARPLOTAVG {
    take:
    ch_metadata
    ch_rel_tsv
    ch_tax
    metadata_category_barplot

    main:
    ch_versions_qiime2_barplotavg = Channel.empty()
    ch_metadata_category_barplot = Channel.fromList(metadata_category_barplot.tokenize(','))

    //Import raltive ASV table
    QIIME2_INASV_BPAVG ( ch_rel_tsv )
    ch_versions_qiime2_barplotavg = ch_versions_qiime2_barplotavg.mix(QIIME2_INASV_BPAVG.out.versions)

    //group by metadata category (ch_metadata_category_barplot)
    QIIME2_FEATURETABLE_GROUP (
        QIIME2_INASV_BPAVG.out.qza
        .combine(ch_metadata)
        .combine(ch_metadata_category_barplot)
    )
    ch_versions_qiime2_barplotavg = ch_versions_qiime2_barplotavg.mix(QIIME2_FEATURETABLE_GROUP.out.versions)

    //Barplot
    QIIME2_BPAVG ( [], QIIME2_FEATURETABLE_GROUP.out.qza, ch_tax, 'average' )
    ch_versions_qiime2_barplotavg = ch_versions_qiime2_barplotavg.mix(QIIME2_BPAVG.out.versions)

    emit:
    versions = ch_versions_qiime2_barplotavg
}
