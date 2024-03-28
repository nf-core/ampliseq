/*
 * Taxonomic classification with SINTAX
 */

include { VSEARCH_SINTAX                } from '../../modules/nf-core/vsearch/sintax/main'

include { FORMAT_TAXONOMY_SINTAX        } from '../../modules/local/format_taxonomy_sintax'
include { FORMAT_TAXRESULTS_SINTAX      } from '../../modules/local/format_taxresults_sintax'

workflow SINTAX_TAXONOMY_WF {
    take:
    ch_sintax_ref_taxonomy
    val_sintax_ref_taxonomy
    ch_fasta
    ch_full_fasta
    sintax_taxlevels

    main:
    ch_versions_sintax_taxonomy = Channel.empty()

    //format taxonomy file
    FORMAT_TAXONOMY_SINTAX ( ch_sintax_ref_taxonomy )
    ch_versions_sintax_taxonomy = ch_versions_sintax_taxonomy.mix(FORMAT_TAXONOMY_SINTAX.out.versions)
    ch_sintaxdb = FORMAT_TAXONOMY_SINTAX.out.db

    //set file prefix
    if (params.cut_its == "none") {
        ASV_tax_name = "ASV_tax_sintax.${val_sintax_ref_taxonomy}"
        ASV_tax_name2 = "ASV_tax_sintax.${val_sintax_ref_taxonomy}"
    } else {
        ASV_tax_name = "ASV_ITS_tax_sintax.${val_sintax_ref_taxonomy}"
        ASV_tax_name2 = "ASV_tax_sintax.${val_sintax_ref_taxonomy}"
    }

    //search taxonomy database with SINTAX
    ch_fasta
        .map {
            fasta ->
                def meta = [:]
                meta.id = ASV_tax_name + ".raw"
                [ meta, fasta ] }
        .set { ch_fasta_sintax }
    VSEARCH_SINTAX( ch_fasta_sintax, ch_sintaxdb )
    ch_versions_sintax_taxonomy = ch_versions_sintax_taxonomy.mix(VSEARCH_SINTAX.out.versions)

    //convert SINTAX output to DADA2 like taxonomy table
    FORMAT_TAXRESULTS_SINTAX( VSEARCH_SINTAX.out.tsv, ch_full_fasta, ASV_tax_name2 + '.tsv', sintax_taxlevels )
    ch_versions_sintax_taxonomy = ch_versions_sintax_taxonomy.mix(FORMAT_TAXRESULTS_SINTAX.out.versions)
    ch_sintax_tax = FORMAT_TAXRESULTS_SINTAX.out.tsv

    emit:
    tax      = ch_sintax_tax
    versions = ch_versions_sintax_taxonomy
}
