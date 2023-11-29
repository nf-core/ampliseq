/*
 * Training of a classifier with QIIME2
 */

include { UNTAR                 } from '../../modules/nf-core/untar/main'
include { GZIP_DECOMPRESS       } from '../../modules/local/gzip_decompress.nf'
include { FORMAT_TAXONOMY_QIIME } from '../../modules/local/format_taxonomy_qiime'
include { QIIME2_EXTRACT        } from '../../modules/local/qiime2_extract'
include { QIIME2_TRAIN          } from '../../modules/local/qiime2_train'

workflow QIIME2_PREPTAX {
    take:
    ch_qiime_ref_taxonomy //channel, list of files
    val_qiime_ref_taxonomy //val
    FW_primer //val
    RV_primer //val

    main:
    ch_qiime2_preptax_versions = Channel.empty()

    if (params.qiime_ref_tax_custom) {
        ch_qiime_ref_taxonomy.view()

        // ch_qiime_ref_taxonomy
        //     .branch {
        //         gzip: it.isFile() && ( it.getName().endsWith(".gz") )
        //         decompressed: it.isFile() && ( it.getName().endsWith(".fna") || it.getName().endsWith (".tax") )
        //         failed: true
        //     }.set { ch_qiime_ref_tax_branched }
        // ch_qiime_ref_tax_branched.failed.subscribe { error "$it is neither a compressed or decompressed sequence or taxonomy file. Please review input." }

        // GZIP_DECOMPRESS(ch_qiime_ref_tax_branched.gzip)
        // ch_qiime2_preptax_versions = ch_qiime2_preptax_versions.mix(GZIP_DECOMPRESS.out.versions)

        // ch_qiime_db_files = GZIP_DECOMPRESS.out.ungzip
        // ch_qiime_db_files = ch_qiime_db_files.mix(ch_qiime_ref_tax_branched.decompressed)

        // ch_ref_database = ch_qiime_db_files.collate(2)
    } else {
        FORMAT_TAXONOMY_QIIME ( ch_qiime_ref_taxonomy.collect() )

        ch_ref_database = FORMAT_TAXONOMY_QIIME.out.fasta.combine(FORMAT_TAXONOMY_QIIME.out.tax)
    }

    ch_ref_database
        .map {
            db ->
                def meta = [:]
                meta.FW_primer = FW_primer
                meta.RV_primer = RV_primer
                [ meta, db ] }
        .set { ch_ref_database }
    QIIME2_EXTRACT ( ch_ref_database )
    QIIME2_TRAIN ( QIIME2_EXTRACT.out.qza )

    emit:
    classifier      = QIIME2_TRAIN.out.qza
    versions        = QIIME2_TRAIN.out.versions
}
