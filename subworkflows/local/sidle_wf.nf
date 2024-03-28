/*
 * Training of a classifier with QIIME2
 */

include { FORMAT_TAXONOMY_SIDLE } from '../../modules/local/format_taxonomy_sidle'
include { SIDLE_INDB            } from '../../modules/local/sidle_indb'
include { SIDLE_INDBALIGNED     } from '../../modules/local/sidle_indbaligned'
include { SIDLE_DBFILT          } from '../../modules/local/sidle_dbfilt'
include { SIDLE_IN              } from '../../modules/local/sidle_in'
include { SIDLE_TRIM            } from '../../modules/local/sidle_trim'
include { SIDLE_DBEXTRACT       } from '../../modules/local/sidle_dbextract'
include { SIDLE_ALIGN           } from '../../modules/local/sidle_align'
include { SIDLE_DBRECON         } from '../../modules/local/sidle_dbrecon'
include { SIDLE_TABLERECON      } from '../../modules/local/sidle_tablerecon'
include { SIDLE_TAXRECON        } from '../../modules/local/sidle_taxrecon'
include { SIDLE_FILTTAX         } from '../../modules/local/sidle_filttax'
include { SIDLE_SEQRECON        } from '../../modules/local/sidle_seqrecon'
include { SIDLE_TREERECON       } from '../../modules/local/sidle_treerecon'

workflow SIDLE_WF {
    take:
    ch_asv_tables_sequences
    ch_sidle_ref_taxonomy
    val_sidle_ref_taxonomy
    ch_db_tree

    main:
    ch_sidle_versions = Channel.empty()

    // DB
    if (!params.sidle_ref_tax_custom) {
        //standard ref taxonomy input from conf/ref_databases.config, one tar.gz / tgz with all files
        FORMAT_TAXONOMY_SIDLE ( ch_sidle_ref_taxonomy, val_sidle_ref_taxonomy )
        ch_sidle_versions = ch_sidle_versions.mix(FORMAT_TAXONOMY_SIDLE.out.versions)
        ch_db_sequences = FORMAT_TAXONOMY_SIDLE.out.seq
        ch_db_alignedsequences = FORMAT_TAXONOMY_SIDLE.out.alnseq
        ch_db_taxonomy = FORMAT_TAXONOMY_SIDLE.out.tax
    } else {
        //input from params.sidle_ref_tax_custom: it[0] = fasta = ch_db_sequences, it[1] = aligned fasta = ch_db_alignedsequences, it[2] = taxonomy txt = ch_db_taxonomy
        ch_db_sequences = ch_sidle_ref_taxonomy.map{ it[0] }
        ch_db_alignedsequences = ch_sidle_ref_taxonomy.map{ it[1] }
        ch_db_taxonomy = ch_sidle_ref_taxonomy.map{ it[2] }
    }
    SIDLE_INDB ( ch_db_sequences, ch_db_taxonomy )
    ch_sidle_versions = ch_sidle_versions.mix(SIDLE_INDB.out.versions)
    SIDLE_INDBALIGNED ( ch_db_alignedsequences )
    ch_sidle_versions = ch_sidle_versions.mix(SIDLE_INDBALIGNED.out.versions)
    SIDLE_DBFILT ( SIDLE_INDB.out.seq, SIDLE_INDB.out.tax )
    ch_sidle_versions = ch_sidle_versions.mix(SIDLE_DBFILT.out.versions)

    // ASV
    SIDLE_IN ( ch_asv_tables_sequences )
    ch_sidle_versions = ch_sidle_versions.mix(SIDLE_IN.out.versions)
    SIDLE_TRIM ( SIDLE_IN.out.table_seq )
    ch_sidle_versions = ch_sidle_versions.mix(SIDLE_TRIM.out.versions)

    // Combine & reconstruct
    SIDLE_DBEXTRACT (
        SIDLE_IN.out.table_seq
            .combine( SIDLE_DBFILT.out.seq )
            .combine( SIDLE_DBFILT.out.tax ) )
    ch_sidle_versions = ch_sidle_versions.mix(SIDLE_DBEXTRACT.out.versions)

    SIDLE_ALIGN ( SIDLE_DBEXTRACT.out.kmers.join(SIDLE_TRIM.out.seq).dump(tag: 'into_SIDLE_ALIGN') )
    ch_sidle_versions = ch_sidle_versions.mix(SIDLE_ALIGN.out.versions)

    SIDLE_DBEXTRACT.out.map
        .join(SIDLE_ALIGN.out.aligned_map)
        .multiMap { meta, map, aligned_map ->
            sampleid: meta.id
            map: map
            aligned_map: aligned_map
        }
        .set { ch_db_reconstruction }

    SIDLE_DBRECON (
        ch_db_reconstruction.sampleid.collect(),
        ch_db_reconstruction.map.collect(),
        ch_db_reconstruction.aligned_map.collect() )
    ch_sidle_versions = ch_sidle_versions.mix(SIDLE_DBRECON.out.versions)

    SIDLE_TRIM.out.table
        .join(SIDLE_ALIGN.out.aligned_map)
        .multiMap { meta, table, aligned_map ->
            sampleid: meta.id
            table: table
            aligned_map: aligned_map
        }
        .set { ch_table_reconstruction }

    // Abundance table
    SIDLE_TABLERECON (
        ch_table_reconstruction.sampleid.collect(),
        ch_table_reconstruction.table.collect(),
        ch_table_reconstruction.aligned_map.collect(),
        SIDLE_DBRECON.out.reconstruction_map,
        SIDLE_DBRECON.out.reconstruction_summary )
    ch_sidle_versions = ch_sidle_versions.mix(SIDLE_TABLERECON.out.versions)

    // Taxonomic classification
    SIDLE_TAXRECON (
        SIDLE_DBRECON.out.reconstruction_map,
        SIDLE_INDB.out.tax )
    ch_sidle_versions = ch_sidle_versions.mix(SIDLE_TAXRECON.out.versions)
    SIDLE_FILTTAX ( SIDLE_TAXRECON.out.tsv, SIDLE_TABLERECON.out.tsv )
    ch_sidle_versions = ch_sidle_versions.mix(SIDLE_FILTTAX.out.versions)

    // Reconstruct sequences/fragments
    // required: aligned sequences file: https://forum.qiime2.org/t/finding-alignment-files-for-sidle/23773/2
    SIDLE_SEQRECON (
        SIDLE_DBRECON.out.reconstruction_map,
        SIDLE_DBRECON.out.reconstruction_summary,
        SIDLE_INDBALIGNED.out.seq )
    // "The output of reconstruct-fragment-rep-seqs provides consensus sequences only if a reference sequence can't be resolved (ids that have a | symbol in them.) It's designed specifically to integrate with the fragment insertion and makes some downstream assumptions, including that you have the same database and insertion tree version.", see https://forum.qiime2.org/t/how-to-merge-q2-sidle-output-with-other-results/22823/2
    ch_sidle_versions = ch_sidle_versions.mix(SIDLE_SEQRECON.out.versions)

    // Reconstruct phylogenetic tree
    SIDLE_TREERECON (
        SIDLE_SEQRECON.out.qza,
        ch_db_tree )
    ch_sidle_versions = ch_sidle_versions.mix(SIDLE_TREERECON.out.versions)

    emit:
    tax_qza        = SIDLE_TAXRECON.out.qza
    tax_tsv        = SIDLE_FILTTAX.out.filtered
    tax_tsv_merged = SIDLE_FILTTAX.out.merged
    table_biom     = SIDLE_TABLERECON.out.biom
    table_qza      = SIDLE_TABLERECON.out.qza
    table_tsv      = SIDLE_TABLERECON.out.tsv
    tree_nwk       = SIDLE_TREERECON.out.nwk
    tree_qza       = SIDLE_TREERECON.out.qza
    versions       = ch_sidle_versions
}
