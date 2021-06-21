/*
 * Diversity indices with QIIME2
 */

params.filterasv_options = [:]
params.ancom_tax_options = [:]
params.ancom_asv_options = [:]

include { QIIME2_FILTERASV                 } from '../../modules/local/qiime2_filterasv'  addParams( options: params.filterasv_options )
include { QIIME2_ANCOM_TAX                 } from '../../modules/local/qiime2_ancom_tax'  addParams( options: params.ancom_tax_options )
include { QIIME2_ANCOM_ASV                 } from '../../modules/local/qiime2_ancom_asv'  addParams( options: params.ancom_asv_options )

workflow QIIME2_ANCOM {
    take:
    ch_metadata
    ch_asv
    ch_metacolumn_all
    ch_tax
    
    main:	
    //Filter ASV table to get rid of samples that have no metadata values
    QIIME2_FILTERASV ( ch_metadata, ch_asv, ch_metacolumn_all )

    //ANCOM on various taxonomic levels
    ch_taxlevel = Channel.from( 2, 3, 4, 5, 6 )
    ch_metadata
        .combine( QIIME2_FILTERASV.out.qza.flatten() )
        .combine( ch_tax )
        .combine( ch_taxlevel )
        .set{ ch_for_ancom_tax }
    QIIME2_ANCOM_TAX ( ch_for_ancom_tax )
    QIIME2_ANCOM_TAX.out.ancom.subscribe { if ( it.baseName[0].toString().startsWith("WARNING") ) log.warn it.baseName[0].toString().replace("WARNING ","QIIME2_ANCOM_TAX: ").replace(" ancom_log", "") }

    QIIME2_ANCOM_ASV ( ch_metadata.combine( QIIME2_FILTERASV.out.qza.flatten() ) )
}
