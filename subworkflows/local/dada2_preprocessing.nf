/*
 * Preprocessing with DADA2
 */

include { DADA2_QUALITY as DADA2_QUALITY1 } from '../../modules/local/dada2_quality'
include { TRUNCLEN                        } from '../../modules/local/trunclen'
include { DADA2_FILTNTRIM                 } from '../../modules/local/dada2_filtntrim'
include { DADA2_QUALITY as DADA2_QUALITY2 } from '../../modules/local/dada2_quality'

workflow DADA2_PREPROCESSING {
    take:
    ch_trimmed_reads
    single_end
    find_truncation_values
    trunclenf
    trunclenr

    main:
    ch_versions_dada2_preprocessing = Channel.empty()

    //plot unprocessed, aggregated quality profile for forward and reverse reads separately
    if (single_end) {
        ch_trimmed_reads
            .map { meta, reads -> [ reads ] }
            .collect()
            .map { reads -> [ "single_end", reads ] }
            .set { ch_all_trimmed_reads }
    } else {
        ch_trimmed_reads
            .map { meta, reads -> [ reads[0] ] }
            .collect()
            .map { reads -> [ "FW", reads ] }
            .set { ch_all_trimmed_fw }
        ch_trimmed_reads
            .map { meta, reads -> [ reads[1] ] }
            .collect()
            .map { reads -> [ "RV", reads ] }
            .set { ch_all_trimmed_rv }
        ch_all_trimmed_fw
            .mix ( ch_all_trimmed_rv )
            .set { ch_all_trimmed_reads }
    }

    ch_DADA2_QUALITY1_SVG = Channel.empty()
    if ( !params.skip_dada_quality ) {
        DADA2_QUALITY1 ( ch_all_trimmed_reads.dump(tag: 'into_dada2_quality') )
        ch_versions_dada2_preprocessing = ch_versions_dada2_preprocessing.mix(DADA2_QUALITY1.out.versions)
        DADA2_QUALITY1.out.warning.subscribe { if ( it.baseName.toString().startsWith("WARNING") ) log.warn it.baseName.toString().replace("WARNING ","DADA2_QUALITY1: ") }
        ch_DADA2_QUALITY1_SVG = DADA2_QUALITY1.out.svg
    }

    //find truncation values in case they are not supplied
    if ( find_truncation_values ) {
        TRUNCLEN ( DADA2_QUALITY1.out.tsv )
        TRUNCLEN.out.trunc
            .toSortedList()
            .set { ch_trunc }
        ch_versions_dada2_preprocessing = ch_versions_dada2_preprocessing.mix(TRUNCLEN.out.versions.first())
        //add one more warning or reminder that trunclenf and trunclenr were chosen automatically
        ch_trunc.subscribe {
            if ( "${it[0][1]}".toInteger() + "${it[1][1]}".toInteger() <= 10 ) { log.warn "`--trunclenf` was set to ${it[0][1]} and `--trunclenr` to ${it[1][1]}, this is too low! Please either change `--trunc_qmin` (and `--trunc_rmin`), or set `--trunclenf` and `--trunclenr`." }
            else if ( "${it[0][1]}".toInteger() <= 10 ) { log.warn "`--trunclenf` was set to ${it[0][1]}, this is too low! Please either change `--trunc_qmin` (and `--trunc_rmin`), or set `--trunclenf` and `--trunclenr`." }
            else if ( "${it[1][1]}".toInteger() <= 10 ) { log.warn "`--trunclenr` was set to ${it[1][1]}, this is too low! Please either change `--trunc_qmin` (and `--trunc_rmin`), or set `--trunclenf` and `--trunclenr`." }
            else log.warn "Probably everything is fine, but this is a reminder that `--trunclenf` was set automatically to ${it[0][1]} and `--trunclenr` to ${it[1][1]}. If this doesnt seem reasonable, then please change `--trunc_qmin` (and `--trunc_rmin`), or set `--trunclenf` and `--trunclenr` directly."
        }
    } else {
        Channel.fromList( [['FW', trunclenf], ['RV', trunclenr]] )
            .toSortedList()
            .set { ch_trunc }
    }
    ch_trimmed_reads.combine(ch_trunc).set { ch_trimmed_reads }

    //filter reads
    DADA2_FILTNTRIM ( ch_trimmed_reads.dump(tag: 'into_filtntrim')  )
    ch_versions_dada2_preprocessing = ch_versions_dada2_preprocessing.mix(DADA2_FILTNTRIM.out.versions.first())

    //Filter empty files
    DADA2_FILTNTRIM.out.reads_logs_args
        .branch {
            failed: it[0].single_end ? it[1].countFastq() < params.min_read_counts : it[1][0].countFastq() < params.min_read_counts || it[1][1].countFastq() < params.min_read_counts
            passed: true
        }
        .set { ch_dada2_filtntrim_results }
    ch_dada2_filtntrim_results.passed.set { ch_dada2_filtntrim_results_passed }
    ch_dada2_filtntrim_results.failed
        .map { meta, reads, logs, args -> [ meta.id ] }
        .collect()
        .subscribe {
            samples = it.join("\n")
            if (params.ignore_failed_filtering) {
                log.warn "The following samples had too few reads (<$params.min_read_counts) after quality filtering with DADA2:\n$samples\nIgnoring failed samples and continue!\n"
            } else {
                error("The following samples had too few reads (<$params.min_read_counts) after quality filtering with DADA2:\n$samples\nPlease check whether the correct primer sequences for trimming were supplied. Ignore that samples using `--ignore_failed_filtering` or adjust the threshold with `--min_read_counts`.")
            }
        }

    // Break apart the reads and logs so that only the samples
    // which pass filtering are retained
    ch_dada2_filtntrim_results_passed
        .map{ meta, reads, logs, args -> [meta, reads] }
        .set{ ch_dada2_filtntrim_reads_passed }
    ch_dada2_filtntrim_results_passed
        .map{ meta, reads, logs, args -> [meta, logs] }
        .set{ ch_dada2_filtntrim_logs_passed }
    ch_dada2_filtntrim_results_passed
        .map{ meta, reads, logs, args -> args }
        .set{ ch_dada2_filtntrim_args_passed }

    //plot post-processing, aggregated quality profile for forward and reverse reads separately
    if (single_end) {
        ch_dada2_filtntrim_reads_passed
            .map { meta, reads -> [ reads ] }
            .collect()
            .map { reads -> [ "single_end", reads ] }
            .set { ch_all_preprocessed_reads }
    } else {
        ch_dada2_filtntrim_reads_passed
            .map { meta, reads -> [ reads[0] ] }
            .collect()
            .map { reads -> [ "FW", reads ] }
            .set { ch_all_preprocessed_fw }
        ch_dada2_filtntrim_reads_passed
            .map { meta, reads -> [ reads[1] ] }
            .collect()
            .map { reads -> [ "RV", reads ] }
            .set { ch_all_preprocessed_rv }
        ch_all_preprocessed_fw
            .mix ( ch_all_preprocessed_rv )
            .set { ch_all_preprocessed_reads }
    }

    ch_DADA2_QUALITY2_SVG = Channel.empty()
    if ( !params.skip_dada_quality ) {
        DADA2_QUALITY2 ( ch_all_preprocessed_reads.dump(tag: 'into_dada2_quality2') )
        ch_versions_dada2_preprocessing = ch_versions_dada2_preprocessing.mix(DADA2_QUALITY2.out.versions)
        DADA2_QUALITY2.out.warning.subscribe { if ( it.baseName.toString().startsWith("WARNING") ) log.warn it.baseName.toString().replace("WARNING ","DADA2_QUALITY2: ") }
        ch_DADA2_QUALITY2_SVG = DADA2_QUALITY2.out.svg
    }

    // group reads by sequencing run and region
    // 'groupTuple', 'size' or 'groupKey' should be used but to produce it we need to know how many elements to group but some can be lost here, so no way knowing before
    ch_dada2_filtntrim_reads_passed
        .map {
            info, reads ->
                def meta = info.subMap( info.keySet() - 'id' - 'sample' )
                [ meta, reads, info.id, info.sample ] }
        .groupTuple(by: 0 )
        .map {
            info, reads, ids, samples ->
                def meta = info + [id: ids.flatten().sort(), sample: samples.flatten().sort()]
                [ meta, reads.flatten().sort() ] }
        .set { ch_filt_reads }

    //group logs by sequencing run and region
    //for 'groupTuple', 'size' or 'groupKey' should be used but to produce it we need to know how many elements to group but some can be lost here, so no way knowing before
    ch_dada2_filtntrim_logs_passed
        .map {
            info, reads ->
                def meta = info.subMap( info.keySet() - 'id' - 'sample' )
                [ meta, reads, info.id, info.sample ] }
        .groupTuple(by: 0 )
        .map {
            info, reads, ids, samples ->
                def meta = info + [id: ids.flatten().sort(), sample: samples.flatten().sort()]
                [ meta, reads.flatten().sort() ] }
        .set { ch_filt_logs }

    emit:
    reads               = ch_filt_reads
    logs                = ch_filt_logs
    args                = ch_dada2_filtntrim_args_passed
    qc_svg              = ch_DADA2_QUALITY1_SVG.collect()
    qc_svg_preprocessed = ch_DADA2_QUALITY2_SVG.collect()
    versions            = ch_versions_dada2_preprocessing
}
