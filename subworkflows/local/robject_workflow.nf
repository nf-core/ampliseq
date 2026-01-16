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
    ch_robject_intree
    run_qiime2

    main:
    ch_versions_robject_workflow = channel.empty()

    if ( run_qiime2 ) {
        if ( params.exclude_taxa != "none" || params.min_frequency != 1 || params.min_samples != 1 ) {
            ch_robject_inasv = PHYLOSEQ_INASV ( ch_tsv ).tsv
        } else {
            ch_robject_inasv = ch_tsv
        }
    } else {
        ch_robject_inasv = ch_tsv
    }

    ch_tax
        .combine(ch_robject_inasv)
        .combine(ch_meta.ifEmpty('EMPTY'))
        .combine(ch_robject_intree.ifEmpty('EMPTY'))
        .map {
            id, tax, tsv, meta, tree ->
                def meta_new = ( meta != 'EMPTY' ? meta : [] )
                def tree_new = ( tree != 'EMPTY' ? tree : [] )
                [ id, tax, tsv, meta_new, tree_new ] }
        .set { ch_for_r_objects }

    if ( !params.skip_phyloseq ) {
        PHYLOSEQ ( ch_for_r_objects )
        ch_versions_robject_workflow = ch_versions_robject_workflow.mix(PHYLOSEQ.out.versions)
    }

    if ( !params.skip_tse ) {
        TREESUMMARIZEDEXPERIMENT ( ch_for_r_objects )
        ch_versions_robject_workflow = ch_versions_robject_workflow.mix(TREESUMMARIZEDEXPERIMENT.out.versions)
    }

    emit:
    phyloseq = !params.skip_phyloseq ? PHYLOSEQ.out.rds : []
    tse      = !params.skip_tse ? TREESUMMARIZEDEXPERIMENT.out.rds : []
    versions = ch_versions_robject_workflow
}
