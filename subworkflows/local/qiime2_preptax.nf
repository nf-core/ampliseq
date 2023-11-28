/*
 * Training of a classifier with QIIME2
 */

include { UNTAR                 } from '../../modules/nf-core/untar/main'
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
    if (params.qiime_ref_tax_custom) {
        ch_qiime_ref_taxonomy
            .branch {
                tar: it.isFile() && ( it.getName().endsWith(".tar.gz") || it.getName().endsWith (".tgz") )
                dir: it.isDirectory()
                failed: true
            }.set { ch_qiime_ref_taxonomy }
        ch_qiime_ref_taxonomy.failed.subscribe { error "$it is neither a directory nor a file that ends in '.tar.gz' or '.tgz'. Please review input." }

        UNTAR (
            ch_qiime_ref_taxonomy.tar
                .map {
                    db ->
                        def meta = [:]
                        meta.id = val_qiime_ref_taxonomy
                        [ meta, db ] } )
        ch_qiime_db_dir = UNTAR.out.untar.map{ it[1] }
        ch_qiime_db_dir = ch_qiime_db_dir.mix(ch_qiime_ref_taxonomy.dir)

        ch_ref_database = ch_qiime_db_dir.map{ Channel.fromPath(it + "/*.tax").combine(Channel.fromPath(it + "/*.fna")) }
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
