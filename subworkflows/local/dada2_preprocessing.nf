/*
 * Preprocessing with DADA2
 */

include { DADA2_QUALITY                   } from '../../modules/local/dada2_quality'
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

    if ( !params.skip_dada_quality ) {
        DADA2_QUALITY ( ch_all_trimmed_reads.dump(tag: 'into_dada2_quality') )
        ch_versions_dada2_preprocessing = ch_versions_dada2_preprocessing.mix(DADA2_QUALITY.out.versions.first())
        DADA2_QUALITY.out.warning.subscribe { if ( it.baseName.toString().startsWith("WARNING") ) log.warn it.baseName.toString().replace("WARNING ","DADA2_QUALITY: ") }
    }

    //find truncation values in case they are not supplied
    if ( find_truncation_values ) {
        TRUNCLEN ( DADA2_QUALITY.out.tsv )
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
        Channel.from( [['FW', trunclenf], ['RV', trunclenr]] )
            .toSortedList()
            .set { ch_trunc }
    }
    ch_trimmed_reads.combine(ch_trunc).set { ch_trimmed_reads }

    //filter reads
    DADA2_FILTNTRIM ( ch_trimmed_reads.dump(tag: 'into_filtntrim')  )
    ch_versions_dada2_preprocessing = ch_versions_dada2_preprocessing.mix(DADA2_FILTNTRIM.out.versions.first())

    //plot post-processing, aggregated quality profile for forward and reverse reads separately
    if (single_end) {
        DADA2_FILTNTRIM.out.reads
            .map { meta, reads -> [ reads ] }
            .collect()
            .map { reads -> [ "single_end", reads ] }
            .set { ch_all_preprocessed_reads }
    } else {
        DADA2_FILTNTRIM.out.reads
            .map { meta, reads -> [ reads[0] ] }
            .collect()
            .map { reads -> [ "FW", reads ] }
            .set { ch_all_preprocessed_fw }
        DADA2_FILTNTRIM.out.reads
            .map { meta, reads -> [ reads[1] ] }
            .collect()
            .map { reads -> [ "RV", reads ] }
            .set { ch_all_preprocessed_rv }
        ch_all_preprocessed_fw
            .mix ( ch_all_preprocessed_rv )
            .set { ch_all_preprocessed_reads }
    }
    if ( !params.skip_dada_quality ) {
        DADA2_QUALITY2 ( ch_all_preprocessed_reads.dump(tag: 'into_dada2_quality2') )
        DADA2_QUALITY2.out.warning.subscribe { if ( it.baseName.toString().startsWith("WARNING") ) log.warn it.baseName.toString().replace("WARNING ","DADA2_QUALITY2: ") }
    }

    //group by sequencing run
    DADA2_FILTNTRIM.out.reads
        .map {
            info, reads ->
                def meta = [:]
                meta.run = info.run
                meta.single_end = info.single_end
                [ meta, reads, info.id ] }
        .groupTuple(by: 0 )
        .map {
            info, reads, ids ->
                def meta = [:]
                meta.run = info.run
                meta.single_end = info.single_end
                meta.id = ids.flatten().sort()
                [ meta, reads.flatten().sort() ] }
        .set { ch_filt_reads }

    emit:
    reads    = ch_filt_reads
    logs      = DADA2_FILTNTRIM.out.log
    versions = ch_versions_dada2_preprocessing
}
