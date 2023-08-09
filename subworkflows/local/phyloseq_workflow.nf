/*
 * Create phyloseq objects
 */

include { PHYLOSEQ                                } from '../../modules/local/phyloseq'
include { PHYLOSEQ_INASV                          } from '../../modules/local/phyloseq_inasv'
include { PHYLOSEQ_INTAX as PHYLOSEQ_INTAX_PPLACE } from '../../modules/local/phyloseq_intax'
include { PHYLOSEQ_INTAX as PHYLOSEQ_INTAX_QIIME2 } from '../../modules/local/phyloseq_intax'

workflow PHYLOSEQ_WORKFLOW {
    take:
    ch_dada2_tax
    ch_sintax_tax
    ch_pplace_tax
    ch_qiime2_tax
    ch_tsv
    ch_meta
    ch_tree
    run_qiime2
    
    main:
    if ( params.metadata ) {
        ch_phyloseq_inmeta = ch_meta.first() // The .first() is to make sure it's a value channel
    } else { 
        ch_phyloseq_inmeta = [] 
    }

    ch_phyloseq_intax = Channel.empty()
    if ( !params.skip_dada_taxonomy ) {
        ch_phyloseq_intax = ch_phyloseq_intax.mix ( 
            ch_dada2_tax.map { it = [ "dada2", file(it) ] } 
        )
    }

    if ( params.sintax_ref_taxonomy ) {
        ch_phyloseq_intax = ch_phyloseq_intax.mix ( 
            ch_sintax_tax.map { it = [ "sintax", file(it) ] } 
        )
    }

    if ( params.pplace_tree ) {
        ch_phyloseq_intax = ch_phyloseq_intax.mix ( 
            PHYLOSEQ_INTAX_PPLACE ( 
                ch_pplace_tax
            ).tsv.map { it = [ "pplace", file(it) ] } 
        )

        ch_phyloseq_intree = ch_tree.map { it = it[1] }.first()
    } else {
        ch_phyloseq_intree = []
    }
        
    if ( run_qiime2 ) {
        ch_phyloseq_intax = ch_phyloseq_intax.mix ( 
            PHYLOSEQ_INTAX_QIIME2 ( 
                ch_qiime2_tax
            ).tsv.map { it = [ "qiime2", file(it) ] } 
        )

        if ( params.exclude_taxa != "none" || params.min_frequency != 1 || params.min_samples != 1 ) {
            ch_phyloseq_inasv = PHYLOSEQ_INASV ( ch_tsv ).tsv

        } else { 
            ch_phyloseq_inasv = ch_tsv 
        }
    } else { 
        ch_phyloseq_inasv = ch_tsv
    }

    PHYLOSEQ ( ch_phyloseq_intax, ch_phyloseq_inasv, ch_phyloseq_inmeta, ch_phyloseq_intree )

    emit:
    rds     = PHYLOSEQ.out.rds
    versions= PHYLOSEQ.out.versions
}
