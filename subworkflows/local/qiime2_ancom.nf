/*
 * Diversity indices with QIIME2
 */

include { QIIME2_FILTERSAMPLES as QIIME2_FILTERSAMPLES_ANCOM } from '../../modules/local/qiime2_filtersamples'
include { QIIME2_ANCOM_TAX                 } from '../../modules/local/qiime2_ancom_tax'
include { QIIME2_ANCOM_ASV                 } from '../../modules/local/qiime2_ancom_asv'
include { QIIME2_ANCOMBC_ASV               } from '../../modules/local/qiime2_ancombc_asv'
include { QIIME2_ANCOMBC_TAX               } from '../../modules/local/qiime2_ancombc_tax'
include { QIIME2_ANCOMBC_ASV as ANCOMBC_FORMULA_ASV } from '../../modules/local/qiime2_ancombc_asv'
include { QIIME2_ANCOMBC_TAX as ANCOMBC_FORMULA_TAX } from '../../modules/local/qiime2_ancombc_tax'

workflow QIIME2_ANCOM {
    take:
    ch_metadata
    ch_asv
    ch_metacolumn_all
    ch_tax
    tax_agglom_min
    tax_agglom_max
    ancombc_formula

    main:
    ch_versions_qiime2_ancom = Channel.empty()

    ch_taxlevel = Channel.of( tax_agglom_min..tax_agglom_max )

    //Filter ASV table to get rid of samples that have no metadata values
    ch_metadata
        .combine( ch_asv )
        .combine( ch_metacolumn_all )
        .set{ ch_for_filtersamples }
    QIIME2_FILTERSAMPLES_ANCOM ( ch_for_filtersamples )
    ch_versions_qiime2_ancom = ch_versions_qiime2_ancom.mix(QIIME2_FILTERSAMPLES_ANCOM.out.versions)

    if ( params.ancom ) {
        //ANCOM on various taxonomic levels
        ch_metadata
            .combine( QIIME2_FILTERSAMPLES_ANCOM.out.qza )
            .combine( ch_tax )
            .combine( ch_taxlevel )
            .set{ ch_for_ancom_tax }
        QIIME2_ANCOM_TAX ( ch_for_ancom_tax )
        ch_versions_qiime2_ancom = ch_versions_qiime2_ancom.mix(QIIME2_ANCOM_TAX.out.versions)
        QIIME2_ANCOM_TAX.out.ancom.subscribe { if ( it.baseName[0].toString().startsWith("WARNING") ) log.warn it.baseName[0].toString().replace("WARNING ","QIIME2_ANCOM_TAX: ") }

        //ANCOM on ASVs
        QIIME2_ANCOM_ASV ( ch_metadata.combine( QIIME2_FILTERSAMPLES_ANCOM.out.qza.flatten() ) )
        ch_versions_qiime2_ancom = ch_versions_qiime2_ancom.mix(QIIME2_ANCOM_ASV.out.versions)
    }

    if ( params.ancombc ) {
        //ANCOMBC on various taxonomic levels
        ch_metadata
            .combine( QIIME2_FILTERSAMPLES_ANCOM.out.qza )
            .combine( ch_tax )
            .combine( ch_taxlevel )
            .combine( Channel.fromList([""]) )
            .set{ ch_for_ancombc_tax }
        QIIME2_ANCOMBC_TAX ( ch_for_ancombc_tax )
        ch_versions_qiime2_ancom = ch_versions_qiime2_ancom.mix(QIIME2_ANCOMBC_TAX.out.versions)
        QIIME2_ANCOMBC_TAX.out.da_barplot.subscribe { if ( it.baseName[0].toString().startsWith("WARNING") ) log.warn it.baseName[0].toString().replace("WARNING ","QIIME2_ANCOMBC_TAX: ") }

        //ANCOMBC on ASVs
        QIIME2_ANCOMBC_ASV ( ch_metadata.combine( QIIME2_FILTERSAMPLES_ANCOM.out.qza.flatten() ).combine( Channel.fromList([""]) ) )
        ch_versions_qiime2_ancom = ch_versions_qiime2_ancom.mix(QIIME2_ANCOMBC_ASV.out.versions)
    }

    if ( ancombc_formula ) {
        ch_ancombc_formula = Channel.fromList( ancombc_formula.toString().replace(" ","").tokenize(',') )

        //ANCOMBC with ancombc_formula on various taxonomic levels
        ch_taxlevel = Channel.of( tax_agglom_min..tax_agglom_max )
        ch_metadata
            .combine( ch_asv )
            .combine( ch_tax )
            .combine( ch_taxlevel )
            .combine( ch_ancombc_formula )
            .set{ ch_for_ancombc_tax }
        ANCOMBC_FORMULA_TAX ( ch_for_ancombc_tax )
        ch_versions_qiime2_ancom = ch_versions_qiime2_ancom.mix(ANCOMBC_FORMULA_TAX.out.versions)
        ANCOMBC_FORMULA_TAX.out.da_barplot.subscribe { if ( it.baseName[0].toString().startsWith("WARNING") ) log.warn it.baseName[0].toString().replace("WARNING ","QIIME2_ANCOMBC_TAX: ") }

        //ANCOMBC with ancombc_formula on ASVs
        ANCOMBC_FORMULA_ASV ( ch_metadata.combine( ch_asv ).combine( ch_ancombc_formula ) )
        ch_versions_qiime2_ancom = ch_versions_qiime2_ancom.mix(ANCOMBC_FORMULA_ASV.out.versions)
    }

    emit:
    ancom    = params.ancom ? QIIME2_ANCOM_ASV.out.ancom.mix(QIIME2_ANCOM_TAX.out.ancom) : Channel.empty()
    ancombc  = params.ancombc ? QIIME2_ANCOMBC_ASV.out.da_barplot.mix(QIIME2_ANCOMBC_TAX.out.da_barplot) : Channel.empty()
    ancombc_formula = ancombc_formula ? ANCOMBC_FORMULA_ASV.out.da_barplot.mix(ANCOMBC_FORMULA_TAX.out.da_barplot) : Channel.empty()
    versions = ch_versions_qiime2_ancom
}
