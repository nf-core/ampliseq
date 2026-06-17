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
include { QIIME2_ANCOMBC2_ASV               } from '../../modules/local/qiime2_ancombc2_asv'
include { QIIME2_ANCOMBC2_TAX               } from '../../modules/local/qiime2_ancombc2_tax'
include { QIIME2_ANCOMBC2_ASV as ANCOMBC2_FORMULA_ASV } from '../../modules/local/qiime2_ancombc2_asv'
include { QIIME2_ANCOMBC2_TAX as ANCOMBC2_FORMULA_TAX } from '../../modules/local/qiime2_ancombc2_tax'

workflow QIIME2_ANCOM {
    take:
    ch_metadata
    ch_asv
    ch_metacolumn_all
    ch_tax
    tax_agglom_min
    tax_agglom_max
    ancombc_formula
    ancombc2_formula

    main:
    ch_taxlevel = channel.of( tax_agglom_min..tax_agglom_max )

    //Filter ASV table to get rid of samples that have no metadata values
    ch_metadata
        .combine( ch_asv )
        .combine( ch_metacolumn_all )
        .set{ ch_for_filtersamples }
    QIIME2_FILTERSAMPLES_ANCOM ( ch_for_filtersamples )

    if ( params.ancom ) {
        //ANCOM on various taxonomic levels
        ch_metadata
            .combine( QIIME2_FILTERSAMPLES_ANCOM.out.qza )
            .combine( ch_tax )
            .combine( ch_taxlevel )
            .set{ ch_for_ancom_tax }
        QIIME2_ANCOM_TAX ( ch_for_ancom_tax )
        QIIME2_ANCOM_TAX.out.ancom.subscribe { it -> if ( it.baseName[0].toString().startsWith("WARNING") ) log.warn it.baseName[0].toString().replace("WARNING ","QIIME2_ANCOM_TAX: ") }

        //ANCOM on ASVs
        QIIME2_ANCOM_ASV ( ch_metadata.combine( QIIME2_FILTERSAMPLES_ANCOM.out.qza.flatten() ) )
    }

    if ( params.ancombc ) {
        //ANCOMBC on various taxonomic levels
        ch_metadata
            .combine( QIIME2_FILTERSAMPLES_ANCOM.out.qza )
            .combine( ch_tax )
            .combine( ch_taxlevel )
            .combine( channel.fromList([""]) )
            .set{ ch_for_ancombc_tax }
        QIIME2_ANCOMBC_TAX ( ch_for_ancombc_tax )
        QIIME2_ANCOMBC_TAX.out.da_barplot.subscribe { it -> if ( it.baseName[0].toString().startsWith("WARNING") ) log.warn it.baseName[0].toString().replace("WARNING ","QIIME2_ANCOMBC_TAX: ") }

        //ANCOMBC on ASVs
        QIIME2_ANCOMBC_ASV ( ch_metadata.combine( QIIME2_FILTERSAMPLES_ANCOM.out.qza.flatten() ).combine( channel.fromList([""]) ) )
    }

    if ( ancombc_formula ) {
        ch_ancombc_formula = channel.fromList( ancombc_formula.toString().replace(" ","").tokenize(',') )

        //ANCOMBC with ancombc_formula on various taxonomic levels
        ch_taxlevel = channel.of( tax_agglom_min..tax_agglom_max )
        ch_metadata
            .combine( ch_asv )
            .combine( ch_tax )
            .combine( ch_taxlevel )
            .combine( ch_ancombc_formula )
            .set{ ch_for_ancombc_tax }
        ANCOMBC_FORMULA_TAX ( ch_for_ancombc_tax )
        ANCOMBC_FORMULA_TAX.out.da_barplot.subscribe { it -> if ( it.baseName[0].toString().startsWith("WARNING") ) log.warn it.baseName[0].toString().replace("WARNING ","QIIME2_ANCOMBC_TAX: ") }

        //ANCOMBC with ancombc_formula on ASVs
        ANCOMBC_FORMULA_ASV ( ch_metadata.combine( ch_asv ).combine( ch_ancombc_formula ) )
    }

    if ( params.ancombc2 ) {
        //ANCOMBC2 on various taxonomic levels
        ch_metadata
            .combine( QIIME2_FILTERSAMPLES_ANCOM.out.qza )
            .combine( ch_tax )
            .combine( ch_taxlevel )
            .combine( channel.fromList([""]) )
            .set{ ch_for_ancombc2_tax }
        QIIME2_ANCOMBC2_TAX ( ch_for_ancombc2_tax )
        QIIME2_ANCOMBC2_TAX.out.plot.subscribe { it -> if ( it.baseName[0].toString().startsWith("WARNING") ) log.warn it.baseName[0].toString().replace("WARNING ","QIIME2_ANCOMBC2_TAX: ") }

        //ANCOMBC2 on ASVs
        QIIME2_ANCOMBC2_ASV ( ch_metadata.combine( QIIME2_FILTERSAMPLES_ANCOM.out.qza.flatten() ).combine( channel.fromList([""]) ) )
    }

    if ( ancombc2_formula ) {
        ch_ancombc2_formula = channel.fromList( ancombc2_formula.toString().replace(" ","").tokenize(',') )

        //ANCOMBC2 with ancombc2_formula on various taxonomic levels
        ch_taxlevel = channel.of( tax_agglom_min..tax_agglom_max )
        ch_metadata
            .combine( ch_asv )
            .combine( ch_tax )
            .combine( ch_taxlevel )
            .combine( ch_ancombc2_formula )
            .set{ ch_for_ancombc2_tax }
        ANCOMBC2_FORMULA_TAX ( ch_for_ancombc2_tax )
        ANCOMBC2_FORMULA_TAX.out.plot.subscribe { it -> if ( it.baseName[0].toString().startsWith("WARNING") ) log.warn it.baseName[0].toString().replace("WARNING ","QIIME2_ANCOMBC2_TAX: ") }

        //ANCOMBC2 with ancombc2_formula on ASVs
        ANCOMBC2_FORMULA_ASV ( ch_metadata.combine( ch_asv ).combine( ch_ancombc2_formula ) )
    }

    emit:
    ancom    = params.ancom ? QIIME2_ANCOM_ASV.out.ancom.mix(QIIME2_ANCOM_TAX.out.ancom) : channel.empty()
    ancombc  = params.ancombc ? QIIME2_ANCOMBC_ASV.out.da_barplot.mix(QIIME2_ANCOMBC_TAX.out.da_barplot) : channel.empty()
    ancombc_formula = ancombc_formula ? ANCOMBC_FORMULA_ASV.out.da_barplot.mix(ANCOMBC_FORMULA_TAX.out.da_barplot) : channel.empty()
    ancombc2  = params.ancombc2 ? QIIME2_ANCOMBC2_ASV.out.plot.mix(QIIME2_ANCOMBC2_TAX.out.plot) : channel.empty()
    ancombc2_formula = ancombc2_formula ? ANCOMBC2_FORMULA_ASV.out.plot.mix(ANCOMBC2_FORMULA_TAX.out.plot) : channel.empty()
}
