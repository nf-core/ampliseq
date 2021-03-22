/*
 * Export filtered tables from QIIME2
 */

params.absolute_options = [:]
params.relasv_options = [:]
params.reltax_options = [:]
params.combine_table_options = [:]

include { QIIME2_EXPORT_ABSOLUTE      } from '../../modules/local/qiime2_export_absolute' addParams( options: params.absolute_options      )
include { QIIME2_EXPORT_RELASV        } from '../../modules/local/qiime2_export_relasv'   addParams( options: params.relasv_options        )
include { QIIME2_EXPORT_RELTAX        } from '../../modules/local/qiime2_export_reltax'   addParams( options: params.reltax_options        )
include { COMBINE_TABLE               } from '../../modules/local/combine_table'          addParams( options: params.combine_table_options )

workflow QIIME2_EXPORT {
    take:
    ch_asv
    ch_seq
    ch_tax
    ch_tax_tsv
    
    main:
    //export_filtered_dada_output (optional)
    QIIME2_EXPORT_ABSOLUTE ( ch_asv, ch_seq, ch_tax )
    
    //RelativeAbundanceASV (optional)
    QIIME2_EXPORT_RELASV ( ch_asv )

    //RelativeAbundanceReducedTaxa (optional)
    QIIME2_EXPORT_RELTAX ( ch_asv, ch_tax )

    //combine_table.r (optional), seems similar to DADA2_table.tsv but with additionally taxonomy merged
    COMBINE_TABLE ( QIIME2_EXPORT_RELASV.out.tsv, QIIME2_EXPORT_ABSOLUTE.out.fasta, ch_tax_tsv )
}