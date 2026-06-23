include { VSEARCH_USEARCHGLOBAL as VSEARCH_USEARCHGLOBAL_BM } from '../../modules/nf-core/vsearch/usearchglobal/main'
include { BENCHMARKING_MATCHES       } from '../../modules/local/benchmarking_matches'

workflow BENCHMARKING_WF {
    take:
    val_md5sum_version     // md5sum of params appended by pipeline version
    query_or_target        // region to evaluate
    ch_detected_sequences  // detected sequences (fasta)
    ch_detected_abundances // detected sequences (abundance table)
    ch_expected_sequences  // expected sequences (fasta)
    ch_expected_abundances // expected sequences (abundance table)

    main:
    // Compare detected sequences to expected sequences (global alignment)
    // alternative:-> "minimap2 -x asm5 -c reference.fasta query.fasta > alignments.paf" where "-cx asm5" are critical params
    def similarity_threshold = "0.80" // similarity threshold for alignment
    VSEARCH_USEARCHGLOBAL_BM (
        ch_detected_sequences.map { it = [ [id: val_md5sum_version], file(it) ] },
        ch_expected_sequences,
        similarity_threshold,
        'userout',
        "query+target+id+alnlen+mism+opens+qilo+qihi+tilo+tihi+ids+gaps+ql+tl+qstrand" )

    // Investigate mismatches per sample, plus barplot (y = number of sequences, x = number of mismatches)
    BENCHMARKING_MATCHES ( VSEARCH_USEARCHGLOBAL_BM.out.tsv, ch_detected_abundances, ch_expected_abundances.ifEmpty([]), similarity_threshold, query_or_target )

    // (2) input: val_md5sum_version, "${val_md5sum_version}_nucleotide_differences.tsv", ch_detected_abundances, ch_expected_abundances
    //     -> compare exact matches per sample: "${val_md5sum_version}_nucleotide_differences.tsv" = 0 & prevalence ch_detected_abundances vs ch_expected_abundances
    //     -> per sample count FN/FP/... & calculate F1/...
    //BENCHMARKING_SEQUENCES ( val_md5sum_version, BENCHMARKING_MATCHES.out.matches, ch_detected_abundances, ch_expected_abundances )

    // (3) input: val_md5sum_version, "${val_md5sum_version}_matches.txt", ch_detected_abundances, ch_expected_abundances
    //     -> compare abundance of exact matches
    //     -> per sample % total deviation from expected, potentially normalized somehow
    //BENCHMARKING_ABUNDANCES ( val_md5sum_version, "${val_md5sum_version}_exact_matches.txt", ch_detected_abundances, ch_expected_abundances )

    emit:
    mismatch_barplot_png = BENCHMARKING_MATCHES.out.png.collect()
    mismatch_barplot_svg = BENCHMARKING_MATCHES.out.svg.collect()
}
