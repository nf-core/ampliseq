include { VSEARCH_USEARCHGLOBAL as VSEARCH_USEARCHGLOBAL_BM } from '../../modules/nf-core/vsearch/usearchglobal/main'
include { COMPARE_SEQUENCES       } from '../../modules/local/compare_sequences'
include { COMPARE_PERFORMANCE     } from '../../modules/local/compare_performance'

workflow COMPARISON_WF {
    take:
    val_md5sum_version     // md5sum of params appended by pipeline version
    trace_report_suffix    // params.trace_report_suffix that allows linkage to files in pipeline_info, such as pipeline_info/params_<trace_report_suffix>.json
    query_or_target        // region to evaluate
    ch_observed_sequences  // observed sequences (fasta)
    ch_observed_abundances // observed sequences (abundance table)
    ch_expected_sequences  // expected sequences (fasta)
    ch_expected_abundances // expected sequences (abundance table)

    main:
    // Compare observed sequences to expected sequences (global alignment)
    // alternative:-> "minimap2 -x asm5 -c reference.fasta query.fasta > alignments.paf" where "-cx asm5" are critical params
    def similarity_threshold = "0.80" // similarity threshold for alignment
    VSEARCH_USEARCHGLOBAL_BM (
        ch_observed_sequences.map { it = [ [id: val_md5sum_version], file(it) ] },
        ch_expected_sequences,
        similarity_threshold,
        'userout',
        "query+target+ql+tl+qilo+qihi+tilo+tihi+gaps+mism+qstrand" )

    // Investigate mismatches per sample, plus barplot (y = number of sequences, x = number of mismatches)
    COMPARE_SEQUENCES ( VSEARCH_USEARCHGLOBAL_BM.out.tsv, ch_observed_abundances, ch_expected_abundances.ifEmpty([]), similarity_threshold, query_or_target, trace_report_suffix )
    COMPARE_SEQUENCES.out.warnings.subscribe{ it ->
            if( it.countLines() > 0 ) { log.warn "about comparing sequences\n\n" + it.splitText().join("") }
        }

    // Calculate absence/presence performance metrics per sample such as precision, recall, F1 score
    // Calculate abundance performance metrics per sample such as spearman's rho, RMSE, Bray-Curtis distance
    COMPARE_PERFORMANCE ( COMPARE_SEQUENCES.out.matches.map { it = [ [id: val_md5sum_version], file(it) ] }, ch_observed_abundances, ch_expected_abundances )

    emit:
    mismatch_barplot_png = COMPARE_SEQUENCES.out.png.collect()
    mismatch_barplot_svg = COMPARE_SEQUENCES.out.svg.collect()
}
