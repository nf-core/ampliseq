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
    FW_primer //val
    RV_primer //val

    main:
    ch_qiime2_preptax_versions = Channel.empty()

    if (params.qiime_ref_tax_custom) {
        // Handle case where we have been provided a pair of filepaths.
        if ("${params.qiime_ref_tax_custom}".contains(",")) {
            ch_qiime_ref_taxonomy.flatten()
                .branch {
                    compressed: it.isFile() && it.getName().endsWith(".gz")
                    decompressed: it.isFile() && ( it.getName().endsWith(".fna") || it.getName().endsWith(".tax") )
                    failed: true
                }.set { ch_qiime_ref_tax_branched }
            ch_qiime_ref_tax_branched.failed.subscribe { error "$it is neither a compressed (ends with `.gz`) or decompressed sequence (ends with `.fna`) or taxonomy file (ends with `.tax`). Please review input." }

            PIGZ_UNCOMPRESS(ch_qiime_ref_tax_branched.compressed)
            ch_qiime2_preptax_versions = ch_qiime2_preptax_versions.mix(PIGZ_UNCOMPRESS.out.versions)

            ch_qiime_db_files = PIGZ_UNCOMPRESS.out.file
            ch_qiime_db_files = ch_qiime_db_files.mix(ch_qiime_ref_tax_branched.decompressed)

            ch_ref_database_fna = ch_qiime_db_files.filter {
                it.getName().endsWith(".fna")
            }
            ch_ref_database_tax = ch_qiime_db_files.filter {
                it.getName().endsWith(".tax")
            }

            ch_ref_database = ch_ref_database_fna.combine(ch_ref_database_tax)
        // Handle case we have been provided a single filepath (tarball or directory).
        } else {
            ch_qiime_ref_taxonomy.flatten()
                .branch {
                    tar: it.isFile() && ( it.getName().endsWith(".tar.gz") || it.getName().endsWith (".tgz") )
                    dir: it.isDirectory()
                    failed: true
                }.set { ch_qiime_ref_tax_branched }
            ch_qiime_ref_tax_branched.failed.subscribe { error "$it is neither a directory nor a file that ends in '.tar.gz' or '.tgz'. Please review input." }

            UNTAR (
                ch_qiime_ref_tax_branched.tar
                    .map {
                        db ->
                            def meta = [:]
                            meta.id = val_qiime_ref_taxonomy
                            [ meta, db ] } )
            ch_qiime2_preptax_versions = ch_qiime2_preptax_versions.mix(UNTAR.out.versions)

            ch_qiime_db_dir = UNTAR.out.untar.map{ it[1] }
            ch_qiime_db_dir = ch_qiime_db_dir.mix(ch_qiime_ref_tax_branched.dir)

            ch_ref_database_fna = ch_qiime_db_dir.map{ dir ->
                files = file(dir.resolve("*.fna"), checkIfExists: true)
            } | filter {
                if (it.size() > 1) log.warn "Found multiple fasta files for QIIME2 reference database."
                it.size() == 1
            }
            ch_ref_database_tax = ch_qiime_db_dir.map{ dir ->
                files = file(dir.resolve("*.tax"), checkIfExists: true)
            } | filter {
                if (it.size() > 1) log.warn "Found multiple tax files for QIIME2 reference database."
                it.size() == 1
            }

            ch_ref_database = ch_ref_database_fna.combine(ch_ref_database_tax)
        }
    } else {
        FORMAT_TAXONOMY_QIIME ( ch_qiime_ref_taxonomy )
        ch_qiime2_preptax_versions = ch_qiime2_preptax_versions.mix(FORMAT_TAXONOMY_QIIME.out.versions)

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
    ch_qiime2_preptax_versions = ch_qiime2_preptax_versions.mix(QIIME2_EXTRACT.out.versions)

    QIIME2_TRAIN ( QIIME2_EXTRACT.out.qza )
    ch_qiime2_preptax_versions = ch_qiime2_preptax_versions.mix(QIIME2_TRAIN.out.versions)

    emit:
    classifier = QIIME2_TRAIN.out.qza
    versions   = ch_qiime2_preptax_versions
}

