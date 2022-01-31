/*
 * Check input samplesheet or folder and get read channels
 */

include { CUTADAPT as CUTADAPT_BASIC                        } from '../../modules/nf-core/modules/cutadapt/main'
include { CUTADAPT as CUTADAPT_READTHROUGH                  } from '../../modules/nf-core/modules/cutadapt/main'
include { CUTADAPT as CUTADAPT_DOUBLEPRIMER                 } from '../../modules/nf-core/modules/cutadapt/main'
include { CUTADAPT_SUMMARY                                  } from '../../modules/local/cutadapt_summary'
include { CUTADAPT_SUMMARY as CUTADAPT_SUMMARY_DOUBLEPRIMER } from '../../modules/local/cutadapt_summary'
include { CUTADAPT_SUMMARY_MERGE                            } from '../../modules/local/cutadapt_summary_merge'

workflow CUTADAPT_WORKFLOW {
    take:
    ch_file
    illumina_pe_its
    double_primer
    main:
    CUTADAPT_BASIC ( ch_file ).reads.set { ch_trimmed_reads }
    CUTADAPT_BASIC.out.log
        .map {
            info, log ->
                def meta = [:]
                meta.single_end = info.single_end
                [ meta, log ] }
        .groupTuple(by: 0 )
        .set { ch_cutadapt_logs }
    CUTADAPT_SUMMARY ( "cutadapt_standard", ch_cutadapt_logs )

    if (illumina_pe_its) {
        CUTADAPT_READTHROUGH ( ch_trimmed_reads ).reads.set { ch_trimmed_reads }
    }

    if (double_primer) {
        CUTADAPT_DOUBLEPRIMER ( ch_trimmed_reads ).reads.set { ch_trimmed_reads }
        CUTADAPT_DOUBLEPRIMER.out.log
            .map {
                info, log ->
                    def meta = [:]
                    meta.single_end = info.single_end
                    [ meta, log ] }
            .groupTuple(by: 0 )
            .set { ch_cutadapt_doubleprimer_logs }
        CUTADAPT_SUMMARY_DOUBLEPRIMER ( "cutadapt_doubleprimer", ch_cutadapt_doubleprimer_logs )
        ch_summaries = CUTADAPT_SUMMARY.out.tsv.combine( CUTADAPT_SUMMARY_DOUBLEPRIMER.out.tsv )
        CUTADAPT_SUMMARY_MERGE ( "merge", ch_summaries )
    } else {
        CUTADAPT_SUMMARY_MERGE ( "copy", CUTADAPT_SUMMARY.out.tsv )
    }

    //Filter empty files
    ch_trimmed_reads
        .branch {
            failed: it[0].single_end ? it[1].size() < 1.KB : it[1][0].size() < 1.KB || it[1][1].size() < 1.KB
            passed: it[0].single_end ? it[1].size() >= 1.KB : it[1][0].size() >= 1.KB && it[1][1].size() >= 1.KB
        }
        .set { ch_trimmed_reads_result }
    ch_trimmed_reads_result.passed.set { ch_trimmed_reads_passed }
    ch_trimmed_reads_result.failed
        .map { meta, reads -> [ meta.id ] }
        .collect()
        .subscribe {
            samples = it.join("\n")
            log.error "the following samples had too small file size (<1KB) after trimming with cutadapt:\n$samples\nPlease check whether the correct primer sequences for trimming were supplied. Ignore that issue and samples using `--ignore_failed_trimming`."
            params.ignore_failed_trimming ? { log.warn "Ignoring failed samples and continue!" } : System.exit(1)
        }

    emit:
    reads   = ch_trimmed_reads_passed
    logs    = CUTADAPT_BASIC.out.log
    summary = CUTADAPT_SUMMARY_MERGE.out.tsv
    versions= CUTADAPT_BASIC.out.versions
}
