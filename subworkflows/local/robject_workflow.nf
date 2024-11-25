/*
 * Create phyloseq objects
 */

include { PHYLOSEQ                                } from '../../modules/local/phyloseq'
include { PHYLOSEQ_INASV                          } from '../../modules/local/phyloseq_inasv'
include { TREESUMMARIZEDEXPERIMENT                } from '../../modules/local/treesummarizedexperiment'

workflow ROBJECT_WORKFLOW {
    take:
    ch_tax
    ch_tsv
    ch_meta
    ch_tree
    run_qiime2

    main:
    ch_versions_robject_workflow = Channel.empty()

    if ( params.metadata ) {
        ch_robject_inmeta = ch_meta.first() // The .first() is to make sure it's a value channel
    } else {
        ch_robject_inmeta = []
    }

    if ( params.pplace_tree ) {
        ch_robject_intree = ch_tree.map { it = it[1] }.first()
    } else {
        ch_robject_intree = []
    }

    if ( run_qiime2 ) {
        if ( params.exclude_taxa != "none" || params.min_frequency != 1 || params.min_samples != 1 ) {
            ch_robject_inasv = PHYLOSEQ_INASV ( ch_tsv ).tsv
        } else {
            ch_robject_inasv = ch_tsv
        }
    } else {
        ch_robject_inasv = ch_tsv
    }

    if ( params.phyloseq ) {
        PHYLOSEQ ( ch_tax.combine(ch_robject_inasv), ch_robject_inmeta, ch_robject_intree )
        ch_versions_robject_workflow = ch_versions_robject_workflow.mix(PHYLOSEQ.out.versions)
    }

    if ( params.treesummarizedexperiment ) {
        TREESUMMARIZEDEXPERIMENT ( ch_tax.combine(ch_robject_inasv), ch_robject_inmeta, ch_robject_intree )
        ch_versions_robject_workflow = ch_versions_robject_workflow.mix(TREESUMMARIZEDEXPERIMENT.out.versions)
    }

    emit:
    phyloseq                 = params.phyloseq ? PHYLOSEQ.out.rds : []
    treesummarizedexperiment = params.treesummarizedexperiment ? TREESUMMARIZEDEXPERIMENT.out.rds : []
    versions                 = ch_versions_robject_workflow
}
