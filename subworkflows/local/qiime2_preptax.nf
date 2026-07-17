/*
 * Training of a classifier with QIIME2
 */

include { UNTAR                 } from '../../modules/nf-core/untar/main'
include { PIGZ_UNCOMPRESS       } from '../../modules/nf-core/pigz/uncompress/main'
include { FORMAT_TAXONOMY_QIIME } from '../../modules/local/format_taxonomy_qiime'
include { QIIME2_EXTRACT        } from '../../modules/local/qiime2_extract'
include { QIIME2_TRAIN          } from '../../modules/local/qiime2_train'

workflow QIIME2_PREPTAX {
    take:
    ch_qiime_ref_taxonomy //channel, list of files
    val_qiime_ref_taxonomy //val
    fw_primer //val
    rv_primer //val

    main:
    if (params.qiime_ref_tax_custom) {
        // Handle case where we have been provided a pair of filepaths.
        if ("${params.qiime_ref_tax_custom}".contains(",")) {
            ch_qiime_ref_tax_branched =
                ch_qiime_ref_taxonomy.flatten()
                    .branch { it ->
                        compressed: it.isFile() && it.getName().endsWith(".gz")
                        decompressed: it.isFile() && ( it.getName().endsWith(".fna") || it.getName().endsWith(".tax") )
                        failed: true }
            ch_qiime_ref_tax_branched.failed.subscribe { it -> error "$it is neither a compressed (ends with `.gz`) or decompressed sequence (ends with `.fna`) or taxonomy file (ends with `.tax`). Please review input." }

            PIGZ_UNCOMPRESS(ch_qiime_ref_tax_branched.compressed.map{ it -> [[:], it] })

            ch_qiime_db_files = PIGZ_UNCOMPRESS.out.file.map{ it -> it[1] }
            ch_qiime_db_files = ch_qiime_db_files.mix(ch_qiime_ref_tax_branched.decompressed)

            ch_ref_database_fna = ch_qiime_db_files.filter { it ->
                it.getName().endsWith(".fna")
            }
            ch_ref_database_tax = ch_qiime_db_files.filter { it ->
                it.getName().endsWith(".tax")
            }

            ch_ref_database = ch_ref_database_fna.combine(ch_ref_database_tax)
        // Handle case we have been provided a single filepath (tarball or directory).
        } else {
            ch_qiime_ref_tax_branched =
                ch_qiime_ref_taxonomy.flatten()
                    .branch { it ->
                        tar: it.isFile() && ( it.getName().endsWith(".tar.gz") || it.getName().endsWith (".tgz") )
                        dir: it.isDirectory()
                        failed: true }
            ch_qiime_ref_tax_branched.failed.subscribe { it -> error "$it is neither a directory nor a file that ends in '.tar.gz' or '.tgz'. Please review input." }

            UNTAR (
                ch_qiime_ref_tax_branched.tar
                    .map {
                        db ->
                            def meta = [:]
                            meta.id = val_qiime_ref_taxonomy
                            [ meta, db ] } )
            ch_qiime_db_dir = UNTAR.out.untar.map{ it -> it[1] }
            ch_qiime_db_dir = ch_qiime_db_dir.mix(ch_qiime_ref_tax_branched.dir)

            ch_ref_database_fna = ch_qiime_db_dir.map{ dir ->
                file(dir.resolve("*.fna"), checkIfExists: true)
            } | filter { it ->
                if (it.size() > 1) log.warn "Found multiple fasta files for QIIME2 reference database."
                it.size() == 1
            }
            ch_ref_database_tax = ch_qiime_db_dir.map{ dir ->
                file(dir.resolve("*.tax"), checkIfExists: true)
            } | filter { it ->
                if (it.size() > 1) log.warn "Found multiple tax files for QIIME2 reference database."
                it.size() == 1
            }

            ch_ref_database = ch_ref_database_fna.combine(ch_ref_database_tax)
        }
    } else {
        FORMAT_TAXONOMY_QIIME ( ch_qiime_ref_taxonomy )
        ch_ref_database = FORMAT_TAXONOMY_QIIME.out.fasta.combine(FORMAT_TAXONOMY_QIIME.out.tax)
    }

    ch_ref_database =
        ch_ref_database
            .map {
                db ->
                    def meta = [:]
                    meta.fw_primer = fw_primer
                    meta.rv_primer = rv_primer
                    [ meta, db ] }

    QIIME2_EXTRACT ( ch_ref_database )
    QIIME2_TRAIN ( QIIME2_EXTRACT.out.qza )

    emit:
    classifier = QIIME2_TRAIN.out.qza
}
