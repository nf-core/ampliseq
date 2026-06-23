/*
 * Taxonomic classification with VSEARCH usearch_global
 */

include { VSEARCH_USEARCHGLOBAL as VSEARCH_USEARCHGLOBAL_LCA } from '../../modules/nf-core/vsearch/usearchglobal/main'
include { FORMAT_TAXONOMY_VSEARCH_LCA } from '../../modules/local/format_taxonomy_vsearch_lca'
include { FORMAT_TAXRESULTS_VSEARCH_LCA } from '../../modules/local/format_taxresults_vsearch_lca'

workflow VSEARCH_LCA_TAXONOMY_WF {
    take:
    ch_ref_taxonomy
    val_ref_taxonomy
    ch_fasta
    ch_full_fasta
    vsearch_lca_taxlevels
    vsearch_lca_id

    main:
    FORMAT_TAXONOMY_VSEARCH_LCA ( ch_ref_taxonomy )
    ch_lca_db = FORMAT_TAXONOMY_VSEARCH_LCA.out.db

    if (params.cut_its == "none") {
        ASV_tax_name = "ASV_tax_vsearch_lca.${val_ref_taxonomy}"
        ASV_tax_name2 = "ASV_tax_vsearch_lca.${val_ref_taxonomy}"
    } else {
        ASV_tax_name = "ASV_ITS_tax_vsearch_lca.${val_ref_taxonomy}"
        ASV_tax_name2 = "ASV_tax_vsearch_lca.${val_ref_taxonomy}"
    }

    ch_fasta_map =
        ch_fasta
            .map { fasta ->
                def meta = [:]
                meta.id = ASV_tax_name
                [ meta, fasta ] }

    VSEARCH_USEARCHGLOBAL_LCA( ch_fasta_map, ch_lca_db, vsearch_lca_id, 'lcaout', "" )
    FORMAT_TAXRESULTS_VSEARCH_LCA( VSEARCH_USEARCHGLOBAL_LCA.out.lca, ch_full_fasta, ASV_tax_name2 + ".tsv", vsearch_lca_taxlevels )

    emit:
    raw_lca = VSEARCH_USEARCHGLOBAL_LCA.out.lca
    tax     = FORMAT_TAXRESULTS_VSEARCH_LCA.out.tsv
}
