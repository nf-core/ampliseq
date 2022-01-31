/*
 * Export filtered tables from QIIME2
 */

include { QIIME2_EXPORT_ABSOLUTE                } from '../../modules/local/qiime2_export_absolute'
include { QIIME2_EXPORT_RELASV                  } from '../../modules/local/qiime2_export_relasv'
include { QIIME2_EXPORT_RELTAX                  } from '../../modules/local/qiime2_export_reltax'
include { COMBINE_TABLE as COMBINE_TABLE_QIIME2 } from '../../modules/local/combine_table'
include { COMBINE_TABLE as COMBINE_TABLE_DADA2  } from '../../modules/local/combine_table'

workflow QIIME2_EXPORT {
    take:
    ch_asv
    ch_seq
    ch_tax
    ch_QIIME2_tax_tsv
    ch_DADA2_tax_tsv
    tax_agglom_min
    tax_agglom_max

    main:
    //export_filtered_dada_output (optional)
    QIIME2_EXPORT_ABSOLUTE ( ch_asv, ch_seq, ch_tax, tax_agglom_min, tax_agglom_max )

    //RelativeAbundanceASV (optional)
    QIIME2_EXPORT_RELASV ( ch_asv )

    //RelativeAbundanceReducedTaxa (optional)
    QIIME2_EXPORT_RELTAX ( ch_asv, ch_tax, tax_agglom_min, tax_agglom_max )

    //combine_table.r (optional), similar to DADA2_table.tsv but with additionally taxonomy merged
    COMBINE_TABLE_QIIME2 ( QIIME2_EXPORT_RELASV.out.tsv, QIIME2_EXPORT_ABSOLUTE.out.fasta, ch_QIIME2_tax_tsv, 'rel-table-ASV_with-QIIME2-tax.tsv' )

    //combine_table.r (optional), similar to DADA2_table.tsv but with additionally taxonomy merged
    COMBINE_TABLE_DADA2 ( QIIME2_EXPORT_RELASV.out.tsv, QIIME2_EXPORT_ABSOLUTE.out.fasta, ch_DADA2_tax_tsv, 'rel-table-ASV_with-DADA2-tax.tsv' )

    emit:
    abs_fasta      = QIIME2_EXPORT_ABSOLUTE.out.fasta
    abs_tsv        = QIIME2_EXPORT_ABSOLUTE.out.tsv
    rel_tsv        = QIIME2_EXPORT_RELASV.out.tsv
}
