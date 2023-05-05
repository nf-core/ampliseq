/*
 * Phylogenetic placement and taxonomic classification
 */

include { FASTA_NEWICK_EPANG_GAPPA      } from '../../subworkflows/nf-core/fasta_newick_epang_gappa/main'

include { FORMAT_PPLACETAX              } from '../../modules/local/format_pplacetax'

workflow PPLACE_TAXONOMY_WF {
    take:
    ch_fasta

    main:
    ch_versions_pplace_taxonomy = Channel.empty()

    ch_pp_data = ch_fasta.map { it ->
        [ meta: [ id: params.pplace_name ?: 'user_tree' ],
        data: [
            alignmethod:  params.pplace_alnmethod ?: 'hmmer',
            queryseqfile: it,
            refseqfile:   file( params.pplace_aln, checkIfExists: true ),
            hmmfile:      [],
            refphylogeny: file( params.pplace_tree, checkIfExists: true ),
            model:        params.pplace_model,
            taxonomy:     params.pplace_taxonomy ? file( params.pplace_taxonomy, checkIfExists: true ) : []
        ] ]
    }
    FASTA_NEWICK_EPANG_GAPPA ( ch_pp_data )
    ch_versions_pplace_taxonomy = ch_versions_pplace_taxonomy.mix( FASTA_NEWICK_EPANG_GAPPA.out.versions )

    ch_pplace_tax = FORMAT_PPLACETAX ( FASTA_NEWICK_EPANG_GAPPA.out.taxonomy_per_query ).tsv

    emit:
    tax      = ch_pplace_tax
    versions = ch_versions_pplace_taxonomy
}
