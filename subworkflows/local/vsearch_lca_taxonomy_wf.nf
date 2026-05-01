/*
 * Taxonomic classification with VSEARCH usearch_global
 */

include { VSEARCH_USEARCHGLOBAL as VSEARCH_USEARCHGLOBAL_LCA } from '../../modules/nf-core/vsearch/usearchglobal/main'
include { FORMAT_TAXRESULTS_VSEARCH_LCA } from '../../modules/local/format_taxresults_vsearch_lca'

workflow VSEARCH_LCA_TAXONOMY_WF {
    take:
    ch_ref_taxonomy
    ch_fasta
    ch_full_fasta
    vsearch_lca_taxlevels
    vsearch_lca_id

    main:
    if (params.cut_its == "none") {
        ASV_tax_name = "ASV_tax_vsearch_lca"
        ASV_tax_name2 = "ASV_tax_vsearch_lca"
    } else {
        ASV_tax_name = "ASV_ITS_tax_vsearch_lca"
        ASV_tax_name2 = "ASV_tax_vsearch_lca"
    }

    ch_fasta
        .map { fasta ->
            def meta = [:]
            meta.id = ASV_tax_name
            [ meta, fasta ]
        }
        .set { ch_fasta_map }

    VSEARCH_USEARCHGLOBAL_LCA( ch_fasta_map, ch_ref_taxonomy, vsearch_lca_id, 'lcaout', "" )
    FORMAT_TAXRESULTS_VSEARCH_LCA( VSEARCH_USEARCHGLOBAL_LCA.out.lca, ch_full_fasta, ASV_tax_name2 + ".tsv", vsearch_lca_taxlevels )

    emit:
    raw_lca = VSEARCH_USEARCHGLOBAL_LCA.out.lca
    tax     = FORMAT_TAXRESULTS_VSEARCH_LCA.out.tsv
}
