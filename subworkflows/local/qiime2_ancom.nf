/*
 * Diversity indices with QIIME2
 */

include { QIIME2_FILTERASV                 } from '../../modules/local/qiime2_filterasv'
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
    //Filter ASV table to get rid of samples that have no metadata values
    QIIME2_FILTERASV ( ch_metadata, ch_asv, ch_metacolumn_all )

    //ANCOM on various taxonomic levels
    ch_taxlevel = Channel.from( tax_agglom_min..tax_agglom_max )
    ch_metadata
        .combine( QIIME2_FILTERASV.out.qza.flatten() )
        .combine( ch_tax )
        .combine( ch_taxlevel )
        .set{ ch_for_ancom_tax }
    QIIME2_ANCOM_TAX ( ch_for_ancom_tax )
    QIIME2_ANCOM_TAX.out.ancom.subscribe { if ( it.baseName[0].toString().startsWith("WARNING") ) log.warn it.baseName[0].toString().replace("WARNING ","QIIME2_ANCOM_TAX: ") }

    QIIME2_ANCOM_ASV ( ch_metadata.combine( QIIME2_FILTERASV.out.qza.flatten() ) )
}
