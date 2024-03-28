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
    ch_kraken2_ref_taxonomy
        .branch {
            tar: it.isFile() && ( it.getName().endsWith(".tar.gz") || it.getName().endsWith (".tgz") )
            dir: it.isDirectory()
            failed: true
        }.set { ch_kraken2_ref_taxonomy }
    ch_kraken2_ref_taxonomy.failed.subscribe { error "$it is neither a directory nor a file that ends in '.tar.gz' or '.tgz'. Please review input." }

    UNTAR (
        ch_kraken2_ref_taxonomy.tar
            .map {
                db ->
                    def meta = [:]
                    meta.id = val_kraken2_ref_taxonomy
                    [ meta, db ] } )
    ch_versions_kraken2_taxonomy = ch_versions_kraken2_taxonomy.mix(UNTAR.out.versions)
    ch_kraken2db = UNTAR.out.untar.map{ it[1] }
    ch_kraken2db = ch_kraken2db.mix(ch_kraken2_ref_taxonomy.dir)

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
