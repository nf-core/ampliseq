/*
 * Taxonomic classification with Kraken2
 */

include { UNTAR                     } from '../../modules/nf-core/untar/main'
include { KRAKEN2_KRAKEN2           } from '../../modules/nf-core/kraken2/kraken2/main'
include { FORMAT_TAXRESULTS_KRAKEN2 } from '../../modules/local/format_taxresults_kraken2'

workflow KRAKEN2_TAXONOMY_WF {
    take:
    ch_kraken2_ref_taxonomy
    val_kraken2_ref_taxonomy
    ch_fasta
    kraken2_taxlevels

    main:
    ch_versions_kraken2_taxonomy = Channel.empty()

    // format taxonomy file
    // TODO: accept folder with database files (not '*.tar.gz')
    UNTAR (
        ch_kraken2_ref_taxonomy.dump(tag:'ch_kraken2_ref_taxonomy')
            .map {
                db ->
                    def meta = [:]
                    meta.id = val_kraken2_ref_taxonomy
                    [ meta, db ] } )
    ch_kraken2db = UNTAR.out.untar.map{ it[1] }

    // search taxonomy database with kraken2
    ch_fasta
        .map {
            fasta ->
                def meta = [:]
                meta.id = "ASV_tax.${val_kraken2_ref_taxonomy}"
                meta.single_end = true
                [ meta, fasta ] }
        .set { ch_fasta_kraken2 }
    KRAKEN2_KRAKEN2( ch_fasta_kraken2, ch_kraken2db, false, true )
    ch_versions_kraken2_taxonomy = ch_versions_kraken2_taxonomy.mix(KRAKEN2_KRAKEN2.out.versions)

    // convert kraken2 output to ASV taxonomy table
    FORMAT_TAXRESULTS_KRAKEN2( KRAKEN2_KRAKEN2.out.report, KRAKEN2_KRAKEN2.out.classified_reads_assignment, kraken2_taxlevels )

    emit:
    qiime2_tsv = FORMAT_TAXRESULTS_KRAKEN2.out.qiime2_tsv
    tax_tsv    = FORMAT_TAXRESULTS_KRAKEN2.out.tsv
    versions   = ch_versions_kraken2_taxonomy
}
