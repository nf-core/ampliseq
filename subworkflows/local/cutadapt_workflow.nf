/*
 * Check input samplesheet or folder and get read channels
 */

params.standard_options = [:]
params.readthrough_options = [:]
params.doubleprimer_options = [:]
params.summary_options = [:]
params.summary_merge_options = [:]

include { CUTADAPT                          } from '../../modules/nf-core/software/cutadapt/main' addParams( options: params.standard_options     )
include { CUTADAPT as CUTADAPT_READTHROUGH  } from '../../modules/nf-core/software/cutadapt/main' addParams( options: params.readthrough_options  )
include { CUTADAPT as CUTADAPT_DOUBLEPRIMER } from '../../modules/nf-core/software/cutadapt/main' addParams( options: params.doubleprimer_options )
include { CUTADAPT_SUMMARY                                  } from '../../modules/local/cutadapt_summary' addParams( options: params.summary_options     )
include { CUTADAPT_SUMMARY as CUTADAPT_SUMMARY_DOUBLEPRIMER } from '../../modules/local/cutadapt_summary' addParams( options: params.summary_options     )
include { CUTADAPT_SUMMARY_MERGE            } from '../../modules/local/cutadapt_summary_merge' addParams( options: params.summary_merge_options       )

workflow CUTADAPT_WORKFLOW {
    take:
    ch_file
	illumina_pe_its
    double_primer
    main:
    CUTADAPT ( ch_file ).reads.set { ch_trimmed_reads }
	CUTADAPT.out.log
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

    emit:
    reads   = ch_trimmed_reads
	logs    = CUTADAPT.out.log
	summary = CUTADAPT_SUMMARY_MERGE.out.tsv
	version = CUTADAPT.out.version
}
