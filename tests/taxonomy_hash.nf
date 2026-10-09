include { PHYLOSEQ                 } from '../modules/local/phyloseq/main'
include { TREESUMMARIZEDEXPERIMENT } from '../modules/local/treesummarizedexperiment/main'

workflow TAXONOMY_HASH {
    take:
    ch_input

    main:
    PHYLOSEQ(ch_input)
    TREESUMMARIZEDEXPERIMENT(ch_input)

    emit:
    phyloseq = PHYLOSEQ.out.rds
    tse      = TREESUMMARIZEDEXPERIMENT.out.rds
    versions = PHYLOSEQ.out.versions_phyloseq.mix(TREESUMMARIZEDEXPERIMENT.out.versions_treesummarizedexperiment)
}
