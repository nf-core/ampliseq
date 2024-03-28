/*
 * Diversity indices with QIIME2
 */

include { QIIME2_FILTERSAMPLES as QIIME2_FILTERSAMPLES_ANCOM } from '../../modules/local/qiime2_filtersamples'
include { QIIME2_ANCOM_TAX                 } from '../../modules/local/qiime2_ancom_tax'
include { QIIME2_ANCOM_ASV                 } from '../../modules/local/qiime2_ancom_asv'

workflow QIIME2_ANCOM {
    take:
    ch_metadata
    ch_asv
    ch_metacolumn_all
    ch_tax
    tax_agglom_min
    tax_agglom_max

    main:
    ch_versions_qiime2_ancom = Channel.empty()

    //Filter ASV table to get rid of samples that have no metadata values
    ch_metadata
        .combine( ch_asv )
        .combine( ch_metacolumn_all )
        .set{ ch_for_filtersamples }
    QIIME2_FILTERSAMPLES_ANCOM ( ch_for_filtersamples )
    ch_versions_qiime2_ancom = ch_versions_qiime2_ancom.mix(QIIME2_FILTERSAMPLES_ANCOM.out.versions)

    //ANCOM on various taxonomic levels
    ch_taxlevel = Channel.of( tax_agglom_min..tax_agglom_max )
    ch_metadata
        .combine( QIIME2_FILTERSAMPLES_ANCOM.out.qza )
        .combine( ch_tax )
        .combine( ch_taxlevel )
        .set{ ch_for_ancom_tax }
    QIIME2_ANCOM_TAX ( ch_for_ancom_tax )
    ch_versions_qiime2_ancom = ch_versions_qiime2_ancom.mix(QIIME2_ANCOM_TAX.out.versions)
    QIIME2_ANCOM_TAX.out.ancom.subscribe { if ( it.baseName[0].toString().startsWith("WARNING") ) log.warn it.baseName[0].toString().replace("WARNING ","QIIME2_ANCOM_TAX: ") }

    QIIME2_ANCOM_ASV ( ch_metadata.combine( QIIME2_FILTERSAMPLES_ANCOM.out.qza.flatten() ) )
    ch_versions_qiime2_ancom = ch_versions_qiime2_ancom.mix(QIIME2_ANCOM_ASV.out.versions)

    emit:
    ancom    = QIIME2_ANCOM_ASV.out.ancom.mix(QIIME2_ANCOM_TAX.out.ancom)
    versions = ch_versions_qiime2_ancom
}
