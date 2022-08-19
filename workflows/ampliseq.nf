/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowAmpliseq.initialise(params, log)

// Check input path parameters to see if they exist
// params.input may be: folder, samplesheet, fasta file, and therefore should not appear here (because tests only for "file")
def checkPathParamList = [ params.multiqc_config, params.metadata, params.classifier ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = [
                            file("$projectDir/assets/multiqc_config.yml", checkIfExists: true),
                            file("$projectDir/assets/nf-core-ampliseq_logo_light.png", checkIfExists: true)
                            ]
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    INPUT AND VARIABLES
========================================================================================
*/

// Input

if (params.metadata) {
    ch_metadata = Channel.fromPath("${params.metadata}", checkIfExists: true)
} else { ch_metadata = Channel.empty() }

if (params.classifier) {
    ch_qiime_classifier = Channel.fromPath("${params.classifier}", checkIfExists: true)
} else { ch_qiime_classifier = Channel.empty() }

if (params.dada_ref_tax_custom) {
    //custom ref taxonomy input from params.dada_ref_tax_custom & params.dada_ref_tax_custom_sp
    ch_assigntax = Channel.fromPath("${params.dada_ref_tax_custom}", checkIfExists: true)
    if (params.dada_ref_tax_custom_sp) {
        ch_addspecies = Channel.fromPath("${params.dada_ref_tax_custom_sp}", checkIfExists: true)
    }
    ch_dada_ref_taxonomy = Channel.empty()
} else if (params.dada_ref_taxonomy && !params.skip_taxonomy) {
    //standard ref taxonomy input from params.dada_ref_taxonomy & conf/ref_databases.config
    ch_dada_ref_taxonomy = Channel.fromList(params.dada_ref_databases[params.dada_ref_taxonomy]["file"]).map { file(it) }
    if (params.addsh) {
        ch_shinfo = Channel.fromList(params.dada_ref_databases[params.dada_ref_taxonomy]["shfile"]).map { file(it) }
    }
} else { ch_dada_ref_taxonomy = Channel.empty() }

if (params.qiime_ref_taxonomy && !params.skip_taxonomy && !params.classifier) {
    ch_qiime_ref_taxonomy = Channel.fromList(params.qiime_ref_databases[params.qiime_ref_taxonomy]["file"]).map { file(it) }
} else { ch_qiime_ref_taxonomy = Channel.empty() }


// Set non-params Variables

String[] fasta_extensions = [".fasta", ".fna", ".fa"] // this is the alternative ASV fasta input
is_fasta_input = WorkflowAmpliseq.checkIfFileHasExtension( params.input.toString().toLowerCase(), fasta_extensions )

single_end = params.single_end
if (params.pacbio || params.iontorrent) {
    single_end = true
}

trunclenf = params.trunclenf ?: 0
trunclenr = params.trunclenr ?: 0
if ( !single_end && !params.illumina_pe_its && (params.trunclenf == null || params.trunclenr == null) && !is_fasta_input ) {
    find_truncation_values = true
    log.warn "No DADA2 cutoffs were specified (`--trunclenf` & `--trunclenr`), therefore reads will be truncated where median quality drops below ${params.trunc_qmin} (defined by `--trunc_qmin`) but at least a fraction of ${params.trunc_rmin} (defined by `--trunc_rmin`) of the reads will be retained.\nThe chosen cutoffs do not account for required overlap for merging, therefore DADA2 might have poor merging efficiency or even fail.\n"
} else { find_truncation_values = false }

if ( !is_fasta_input && (!params.FW_primer || !params.RV_primer) && !params.skip_cutadapt ) {
    log.error "Incompatible parameters: `--FW_primer` and `--RV_primer` are required for primer trimming. If primer trimming is not needed, use `--skip_cutadapt`."
    System.exit(1)
}

//use custom taxlevels from --dada_assign_taxlevels or database specific taxlevels if specified in conf/ref_databases.config
if ( params.dada_ref_taxonomy ) {
    taxlevels = params.dada_assign_taxlevels ? "${params.dada_assign_taxlevels}" :
        params.dada_ref_databases[params.dada_ref_taxonomy]["taxlevels"] ?: ""
} else { taxlevels = params.dada_assign_taxlevels ? "${params.dada_assign_taxlevels}" : "" }

//make sure that taxlevels adheres to requirements when mixed with addSpecies
if ( params.dada_ref_taxonomy && !params.skip_dada_addspecies && !params.skip_taxonomy && taxlevels ) {
    if ( !taxlevels.endsWith(",Genus,Species") && !taxlevels.endsWith(",Genus") ) {
        log.error "Incompatible settings: To use exact species annotations, taxonomic levels must end with `,Genus,Species` or `,Genus,Species` but are currently `${taxlevels}`. Taxonomic levels can be set with `--dada_assign_taxlevels`. Skip exact species annotations with `--skip_dada_addspecies`.\n"
        System.exit(1)
    }
}

//only run QIIME2 when taxonomy is actually calculated and all required data is available
if ( !params.enable_conda && !params.skip_taxonomy && !params.skip_qiime ) {
    run_qiime2 = true
} else { run_qiime2 = false }

// Set cutoff to use for SH assignment
if ( params.addsh ) {
    vsearch_cutoff = 0.985
}

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RENAME_RAW_DATA_FILES         } from '../modules/local/rename_raw_data_files'
include { DADA2_FILTNTRIM               } from '../modules/local/dada2_filtntrim'
include { DADA2_QUALITY                 } from '../modules/local/dada2_quality'
include { TRUNCLEN                      } from '../modules/local/trunclen'
include { DADA2_ERR                     } from '../modules/local/dada2_err'
include { DADA2_DENOISING               } from '../modules/local/dada2_denoising'
include { DADA2_RMCHIMERA               } from '../modules/local/dada2_rmchimera'
include { DADA2_STATS                   } from '../modules/local/dada2_stats'
include { DADA2_MERGE                   } from '../modules/local/dada2_merge'
include { BARRNAP                       } from '../modules/local/barrnap'
include { FILTER_SSU                    } from '../modules/local/filter_ssu'
include { FILTER_LEN_ASV                } from '../modules/local/filter_len_asv'
include { MERGE_STATS as MERGE_STATS_FILTERSSU } from '../modules/local/merge_stats'
include { MERGE_STATS as MERGE_STATS_FILTERLENASV } from '../modules/local/merge_stats'
include { FORMAT_TAXONOMY               } from '../modules/local/format_taxonomy'
include { ITSX_CUTASV                   } from '../modules/local/itsx_cutasv'
include { MERGE_STATS                   } from '../modules/local/merge_stats'
include { DADA2_TAXONOMY                } from '../modules/local/dada2_taxonomy'
include { DADA2_ADDSPECIES              } from '../modules/local/dada2_addspecies'
include { ASSIGNSH                      } from '../modules/local/assignsh'
include { FORMAT_TAXRESULTS             } from '../modules/local/format_taxresults'
include { FORMAT_TAXRESULTS as FORMAT_TAXRESULTS_ADDSP } from '../modules/local/format_taxresults'
include { QIIME2_INSEQ                  } from '../modules/local/qiime2_inseq'
include { QIIME2_FILTERTAXA             } from '../modules/local/qiime2_filtertaxa'
include { QIIME2_INASV                  } from '../modules/local/qiime2_inasv'
include { FILTER_STATS                  } from '../modules/local/filter_stats'
include { MERGE_STATS as MERGE_STATS_FILTERTAXA } from '../modules/local/merge_stats'
include { QIIME2_BARPLOT                } from '../modules/local/qiime2_barplot'
include { METADATA_ALL                  } from '../modules/local/metadata_all'
include { METADATA_PAIRWISE             } from '../modules/local/metadata_pairwise'
include { QIIME2_INTAX                  } from '../modules/local/qiime2_intax'
include { PICRUST                       } from '../modules/local/picrust'
include { SBDIEXPORT                    } from '../modules/local/sbdiexport'
include { SBDIEXPORTREANNOTATE          } from '../modules/local/sbdiexportreannotate'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { PARSE_INPUT                   } from '../subworkflows/local/parse_input'
include { QIIME2_PREPTAX                } from '../subworkflows/local/qiime2_preptax'
include { QIIME2_TAXONOMY               } from '../subworkflows/local/qiime2_taxonomy'
include { CUTADAPT_WORKFLOW             } from '../subworkflows/local/cutadapt_workflow'
include { QIIME2_EXPORT                 } from '../subworkflows/local/qiime2_export'
include { QIIME2_DIVERSITY              } from '../subworkflows/local/qiime2_diversity'
include { QIIME2_ANCOM                  } from '../subworkflows/local/qiime2_ancom'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { CUTADAPT as CUTADAPT_TAXONOMY     } from '../modules/nf-core/modules/cutadapt/main'
include { FASTQC                            } from '../modules/nf-core/modules/fastqc/main'
include { MULTIQC                           } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS       } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { VSEARCH_USEARCHGLOBAL             } from '../modules/nf-core/modules/vsearch/usearchglobal/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report      = []

workflow AMPLISEQ {

    ch_versions = Channel.empty()

    //
    // Create a channel for input read files
    //
    PARSE_INPUT ( params.input, is_fasta_input, single_end, params.multiple_sequencing_runs, params.extension )
    ch_reads = PARSE_INPUT.out.reads
    ch_fasta = PARSE_INPUT.out.fasta

    //
    // MODULE: Rename files
    //
    RENAME_RAW_DATA_FILES ( ch_reads )
    ch_versions = ch_versions.mix(RENAME_RAW_DATA_FILES.out.versions.first())

    //
    // MODULE: Run FastQC
    //
    if (!params.skip_fastqc) {
        FASTQC ( RENAME_RAW_DATA_FILES.out.fastq )
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    }

    //
    // MODULE: Cutadapt
    //
    if (!params.skip_cutadapt) {
        CUTADAPT_WORKFLOW (
            RENAME_RAW_DATA_FILES.out.fastq,
            params.illumina_pe_its,
            params.double_primer
        ).reads.set { ch_trimmed_reads }
        ch_versions = ch_versions.mix(CUTADAPT_WORKFLOW.out.versions.first())
    } else {
        ch_trimmed_reads = RENAME_RAW_DATA_FILES.out.fastq
    }

    //
    // SUBWORKFLOW / MODULES : ASV generation with DADA2
    //
    //plot aggregated quality profile for forward and reverse reads separately
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
        DADA2_QUALITY ( ch_all_trimmed_reads )
        DADA2_QUALITY.out.warning.subscribe { if ( it.baseName.toString().startsWith("WARNING") ) log.warn it.baseName.toString().replace("WARNING ","DADA2_QUALITY: ") }
    }

    //find truncation values in case they are not supplied
    if ( find_truncation_values ) {
        TRUNCLEN ( DADA2_QUALITY.out.tsv )
        TRUNCLEN.out.trunc
            .toSortedList()
            .set { ch_trunc }
        ch_versions = ch_versions.mix(TRUNCLEN.out.versions.first())
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
    //ch_versions = ch_versions.mix(DADA2_FILTNTRIM.out.versions.first())

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

    DADA2_ERR ( ch_filt_reads )

    //group by meta
    ch_filt_reads
        .join( DADA2_ERR.out.errormodel )
        .set { ch_derep_errormodel }
    DADA2_DENOISING ( ch_derep_errormodel.dump(tag: 'into_denoising')  )
    ch_versions = ch_versions.mix(DADA2_DENOISING.out.versions.first())

    DADA2_RMCHIMERA ( DADA2_DENOISING.out.seqtab )

    //group by sequencing run & group by meta
    DADA2_FILTNTRIM.out.log
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
        .join( DADA2_DENOISING.out.denoised )
        .join( DADA2_DENOISING.out.mergers )
        .join( DADA2_RMCHIMERA.out.rds )
        .set { ch_track_numbers }
    DADA2_STATS ( ch_track_numbers )

    //merge if several runs, otherwise just publish
    DADA2_MERGE (
        DADA2_STATS.out.stats.map { meta, stats -> stats }.collect(),
        DADA2_RMCHIMERA.out.rds.map { meta, rds -> rds }.collect() )

    //merge cutadapt_summary and dada_stats files
    if (!params.skip_cutadapt) {
        MERGE_STATS (CUTADAPT_WORKFLOW.out.summary, DADA2_MERGE.out.dada2stats)
        ch_stats = MERGE_STATS.out.tsv
    } else {
        ch_stats = DADA2_MERGE.out.dada2stats
    }

    //
    // Modules : Filter rRNA
    // TODO: FILTER_SSU.out.stats needs to be merged still into "overall_summary.tsv"
    //
    if (!params.skip_barrnap && params.filter_ssu) {
        BARRNAP ( DADA2_MERGE.out.fasta )
        ch_versions = ch_versions.mix(BARRNAP.out.versions.ifEmpty(null))
        FILTER_SSU ( DADA2_MERGE.out.fasta, DADA2_MERGE.out.asv, BARRNAP.out.matches )
        MERGE_STATS_FILTERSSU ( ch_stats, FILTER_SSU.out.stats )
        ch_stats = MERGE_STATS_FILTERSSU.out.tsv
        ch_dada2_fasta = FILTER_SSU.out.fasta
        ch_dada2_asv = FILTER_SSU.out.asv
    } else if (!params.skip_barrnap && !params.filter_ssu) {
        BARRNAP ( DADA2_MERGE.out.fasta )
        ch_versions = ch_versions.mix(BARRNAP.out.versions.ifEmpty(null))
        ch_dada2_fasta =  DADA2_MERGE.out.fasta
        ch_dada2_asv = DADA2_MERGE.out.asv
    } else {
        ch_dada2_fasta =  DADA2_MERGE.out.fasta
        ch_dada2_asv = DADA2_MERGE.out.asv
    }

    //
    // Modules : amplicon length filtering
    //
    if (params.min_len_asv || params.max_len_asv) {
        FILTER_LEN_ASV ( ch_dada2_fasta,ch_dada2_asv )
        ch_versions = ch_versions.mix(FILTER_LEN_ASV.out.versions.ifEmpty(null))
        MERGE_STATS_FILTERLENASV ( ch_stats, FILTER_LEN_ASV.out.stats )
        ch_stats = MERGE_STATS_FILTERLENASV.out.tsv
        ch_dada2_fasta = FILTER_LEN_ASV.out.fasta
        ch_dada2_asv = FILTER_LEN_ASV.out.asv
    }

    //
    // SUBWORKFLOW / MODULES : Taxonomic classification with DADA2 and/or QIIME2
    //
    //Alternative entry point for fasta that is being classified
    if ( !is_fasta_input ) {
        ch_fasta = ch_dada2_fasta
    }

    //DADA2
    if (!params.skip_taxonomy) {
        if (!params.dada_ref_tax_custom) {
            //standard ref taxonomy input from conf/ref_databases.config
            FORMAT_TAXONOMY ( ch_dada_ref_taxonomy.collect() )
            ch_assigntax = FORMAT_TAXONOMY.out.assigntax
            ch_addspecies = FORMAT_TAXONOMY.out.addspecies
        }
        //Cut taxonomy to expected amplicon
        if (params.cut_dada_ref_taxonomy) {
            ch_assigntax
                .map {
                    db ->
                        def meta = [:]
                        meta.single_end = true
                        meta.id = "assignTaxonomy"
                        [ meta, db ] }
                .set { ch_assigntax }
            CUTADAPT_TAXONOMY ( ch_assigntax ).reads
                .map { meta, db -> db }
                .set { ch_assigntax }
        }
        if (params.cut_its == "none") {
            DADA2_TAXONOMY ( ch_fasta, ch_assigntax, 'ASV_tax.tsv', taxlevels )
            if (!params.skip_dada_addspecies) {
                DADA2_ADDSPECIES ( DADA2_TAXONOMY.out.rds, ch_addspecies, 'ASV_tax_species.tsv', taxlevels )
                if ( params.addsh ) {
                    ch_fasta
                        .map {
                            fasta ->
                                def meta = [:]
                                meta.id = "ASV.vsearch"
                                [ meta, fasta ] }
                        .set { ch_fasta_map }
                    VSEARCH_USEARCHGLOBAL( ch_fasta_map, ch_assigntax, vsearch_cutoff, 'blast6out', "" )
                    ch_versions = ch_versions.mix(VSEARCH_USEARCHGLOBAL.out.versions.ifEmpty(null))
                    ASSIGNSH( DADA2_ADDSPECIES.out.tsv, ch_shinfo.collect(), VSEARCH_USEARCHGLOBAL.out.txt, 'ASV_tax_species_SH.tsv')
                    ch_versions = ch_versions.mix(ASSIGNSH.out.versions.ifEmpty(null))
                    ch_dada2_tax = ASSIGNSH.out.tsv
                } else {
                    ch_dada2_tax = DADA2_ADDSPECIES.out.tsv
                }
            } else {
                if ( params.addsh ) {
                    ch_fasta
                        .map {
                            fasta ->
                                def meta = [:]
                                meta.id = "ASV.vsearch"
                                [ meta, fasta ] }
                        .set { ch_fasta_map }
                    VSEARCH_USEARCHGLOBAL( ch_fasta_map, ch_assigntax, vsearch_cutoff, 'blast6out', "" )
                    ch_versions = ch_versions.mix(VSEARCH_USEARCHGLOBAL.out.versions.ifEmpty(null))
                    ASSIGNSH( DADA2_TAXONOMY.out.tsv, ch_shinfo.collect(), VSEARCH_USEARCHGLOBAL.out.txt, 'ASV_tax_SH.tsv')
                    ch_versions = ch_versions.mix(ASSIGNSH.out.versions.ifEmpty(null))
                    ch_dada2_tax = ASSIGNSH.out.tsv
                 } else {
                     ch_dada2_tax = DADA2_TAXONOMY.out.tsv
                 }
            }
        //Cut out ITS region if long ITS reads
        } else {
            if (params.cut_its == "full") {
                outfile = params.its_partial ? "ASV_ITS_seqs.full_and_partial.fasta" : "ASV_ITS_seqs.full.fasta"
            }
            else if (params.cut_its == "its1") {
                outfile =  params.its_partial ? "ASV_ITS_seqs.ITS1.full_and_partial.fasta" : "ASV_ITS_seqs.ITS1.fasta"
            }
            else if (params.cut_its == "its2") {
                outfile =  params.its_partial ? "ASV_ITS_seqs.ITS2.full_and_partial.fasta" : "ASV_ITS_seqs.ITS2.fasta"
            }
            ITSX_CUTASV ( ch_fasta, outfile )
            ch_versions = ch_versions.mix(ITSX_CUTASV.out.versions.ifEmpty(null))
            ch_cut_fasta = ITSX_CUTASV.out.fasta
            DADA2_TAXONOMY ( ch_cut_fasta, ch_assigntax, 'ASV_ITS_tax.tsv', taxlevels )
            FORMAT_TAXRESULTS ( DADA2_TAXONOMY.out.tsv, ch_fasta, 'ASV_tax.tsv' )
            ch_versions = ch_versions.mix( FORMAT_TAXRESULTS.out.versions.ifEmpty(null) )
            if (!params.skip_dada_addspecies) {
                DADA2_ADDSPECIES ( DADA2_TAXONOMY.out.rds, ch_addspecies, 'ASV_ITS_tax_species.tsv', taxlevels )
                FORMAT_TAXRESULTS_ADDSP ( DADA2_ADDSPECIES.out.tsv, ch_fasta, 'ASV_tax_species.tsv' )
                if ( params.addsh ) {
                    ch_cut_fasta
                        .map {
                            fasta ->
                                def meta = [:]
                                meta.id = "ASV_cut.vsearch"
                                [ meta, fasta ] }
                        .set { ch_cut_fasta }
                    VSEARCH_USEARCHGLOBAL( ch_cut_fasta, ch_assigntax, vsearch_cutoff, 'blast6out', "" )
                    ch_versions = ch_versions.mix(VSEARCH_USEARCHGLOBAL.out.versions.ifEmpty(null))
                    ASSIGNSH( FORMAT_TAXRESULTS_ADDSP.out.tsv, ch_shinfo.collect(), VSEARCH_USEARCHGLOBAL.out.txt, 'ASV_tax_species_SH.tsv')
                    ch_versions = ch_versions.mix(ASSIGNSH.out.versions.ifEmpty(null))
                    ch_dada2_tax = ASSIGNSH.out.tsv
                } else {
                    ch_dada2_tax = FORMAT_TAXRESULTS_ADDSP.out.tsv
                }
           } else {
                if ( params.addsh ) {
                    ch_cut_fasta
                        .map {
                            fasta ->
                                def meta = [:]
                                meta.id = "ASV_cut.vsearch"
                                [ meta, fasta ] }
                        .set { ch_cut_fasta }
                    VSEARCH_USEARCHGLOBAL( ch_cut_fasta, ch_assigntax, vsearch_cutoff, 'blast6out', "" )
                    ch_versions = ch_versions.mix(VSEARCH_USEARCHGLOBAL.out.versions.ifEmpty(null))
                    ASSIGNSH( FORMAT_TAXRESULTS.out.tsv, ch_shinfo.collect(), VSEARCH_USEARCHGLOBAL.out.txt, 'ASV_tax_SH.tsv')
                    ch_versions = ch_versions.mix(ASSIGNSH.out.versions.ifEmpty(null))
                    ch_dada2_tax = ASSIGNSH.out.tsv
                } else {
                    ch_dada2_tax = FORMAT_TAXRESULTS.out.tsv
                }
            }
        }
    }

    //QIIME2
    if ( run_qiime2 ) {
        if (params.qiime_ref_taxonomy && !params.classifier) {
            QIIME2_PREPTAX (
                ch_qiime_ref_taxonomy.collect(),
                params.FW_primer,
                params.RV_primer
            )
            ch_qiime_classifier = QIIME2_PREPTAX.out.classifier
        }
        QIIME2_TAXONOMY (
            ch_fasta,
            ch_qiime_classifier
        )
        ch_versions = ch_versions.mix( QIIME2_TAXONOMY.out.versions.ifEmpty(null) ) //usually a .first() is here, dont know why this leads here to a warning
    }
    //
    // SUBWORKFLOW / MODULES : Downstream analysis with QIIME2
    //
    if ( run_qiime2 ) {
        //Import ASV abundance table and sequences into QIIME2
        QIIME2_INASV ( ch_dada2_asv )
        QIIME2_INSEQ ( ch_fasta )

        //Import taxonomic classification into QIIME2, if available
        if ( params.skip_taxonomy ) {
            log.info "Skip taxonomy classification"
            ch_tax = Channel.empty()
            tax_agglom_min = 1
            tax_agglom_max = 2
        } else if ( params.dada_ref_taxonomy ) {
            log.info "Use DADA2 taxonomy classification"
            ch_tax = QIIME2_INTAX ( ch_dada2_tax ).qza
            tax_agglom_min = params.dada_tax_agglom_min
            tax_agglom_max = params.dada_tax_agglom_max
        } else if ( params.qiime_ref_taxonomy || params.classifier ) {
            log.info "Use QIIME2 taxonomy classification"
            ch_tax = QIIME2_TAXONOMY.out.qza
            tax_agglom_min = params.qiime_tax_agglom_min
            tax_agglom_max = params.qiime_tax_agglom_max
        } else {
            log.info "Use no taxonomy classification"
            ch_tax = Channel.empty()
            tax_agglom_min = 1
            tax_agglom_max = 2
        }

        //Filtering by taxonomy & prevalence & counts
        if (params.exclude_taxa != "none" || params.min_frequency != 1 || params.min_samples != 1) {
            QIIME2_FILTERTAXA (
                QIIME2_INASV.out.qza,
                QIIME2_INSEQ.out.qza,
                ch_tax,
                params.min_frequency,
                params.min_samples,
                params.exclude_taxa
            )
            FILTER_STATS ( ch_dada2_asv, QIIME2_FILTERTAXA.out.tsv )
            ch_versions = ch_versions.mix( FILTER_STATS.out.versions.ifEmpty(null) )
            MERGE_STATS_FILTERTAXA (ch_stats, FILTER_STATS.out.tsv)
            ch_asv = QIIME2_FILTERTAXA.out.asv
            ch_seq = QIIME2_FILTERTAXA.out.seq
            ch_tsv = QIIME2_FILTERTAXA.out.tsv
        } else {
            ch_asv = QIIME2_INASV.out.qza
            ch_seq = QIIME2_INSEQ.out.qza
            ch_tsv = ch_dada2_asv
        }
        //Export various ASV tables
        if (!params.skip_abundance_tables) {
            QIIME2_EXPORT ( ch_asv, ch_seq, ch_tax, QIIME2_TAXONOMY.out.tsv, ch_dada2_tax, tax_agglom_min, tax_agglom_max )
        }

        if (!params.skip_barplot) {
            QIIME2_BARPLOT ( ch_metadata, ch_asv, ch_tax )
        }

        //Select metadata categories for diversity analysis & ancom
        if (params.metadata_category) {
            ch_metacolumn_all = Channel.from(params.metadata_category.tokenize(','))
            METADATA_PAIRWISE ( ch_metadata ).set { ch_metacolumn_pairwise }
            ch_metacolumn_pairwise = ch_metacolumn_pairwise.splitCsv().flatten()
            ch_metacolumn_pairwise = ch_metacolumn_all.join(ch_metacolumn_pairwise)
        } else if (!params.skip_ancom || !params.skip_diversity_indices) {
            METADATA_ALL ( ch_metadata ).set { ch_metacolumn_all }
            //return empty channel if no appropriate column was found
            ch_metacolumn_all.branch { passed: it != "" }.set { result }
            ch_metacolumn_all = result.passed
            ch_metacolumn_all = ch_metacolumn_all.splitCsv().flatten()
            METADATA_PAIRWISE ( ch_metadata ).set { ch_metacolumn_pairwise }
            ch_metacolumn_pairwise = ch_metacolumn_pairwise.splitCsv().flatten()
        } else {
            ch_metacolumn_all = Channel.empty()
            ch_metacolumn_pairwise = Channel.empty()
        }

        //Diversity indices
        if ( params.metadata && (!params.skip_alpha_rarefaction || !params.skip_diversity_indices) ) {
            QIIME2_DIVERSITY (
                ch_metadata,
                ch_asv,
                ch_seq,
                ch_tsv,
                ch_metacolumn_pairwise,
                ch_metacolumn_all,
                params.skip_alpha_rarefaction,
                params.skip_diversity_indices
            )
        }

        //Perform ANCOM tests
        if ( !params.skip_ancom && params.metadata ) {
            QIIME2_ANCOM (
                ch_metadata,
                ch_asv,
                ch_metacolumn_all,
                ch_tax,
                tax_agglom_min,
                tax_agglom_max
            )
        }
    }

    //
    // MODULE: Predict functional potential of a bacterial community from marker genes with Picrust2
    //
    if ( params.picrust ) {
        if ( run_qiime2 && !params.skip_abundance_tables && ( params.dada_ref_taxonomy || params.qiime_ref_taxonomy || params.classifier ) && !params.skip_taxonomy ) {
            PICRUST ( QIIME2_EXPORT.out.abs_fasta, QIIME2_EXPORT.out.abs_tsv, "QIIME2", "This Picrust2 analysis is based on filtered reads from QIIME2" )
        } else {
            PICRUST ( ch_fasta, ch_dada2_asv, "DADA2", "This Picrust2 analysis is based on unfiltered reads from DADA2" )
        }
        ch_versions = ch_versions.mix(PICRUST.out.versions.ifEmpty(null))
    }

    //
    // MODULE: Export data in SBDI's (Swedish biodiversity infrastructure) format
    //
    if ( params.sbdiexport ) {
        SBDIEXPORT ( ch_dada2_asv, ch_dada2_tax, ch_metadata )
        ch_versions = ch_versions.mix(SBDIEXPORT.out.versions.first())
        SBDIEXPORTREANNOTATE ( ch_dada2_tax )
    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        workflow_summary    = WorkflowAmpliseq.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
        ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
        if (!params.skip_fastqc) {
            ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
        }
        if (!params.skip_cutadapt) {
            ch_multiqc_files = ch_multiqc_files.mix(CUTADAPT_WORKFLOW.out.logs.collect{it[1]}.ifEmpty([]))
        }

        ch_multiqc_configs = Channel.from(ch_multiqc_config).mix(ch_multiqc_custom_config).ifEmpty([])

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_configs.collect()
        )
        multiqc_report = MULTIQC.out.report.toList()
        ch_versions    = ch_versions.mix(MULTIQC.out.versions)
    }

    //Save input in results folder
    input = file(params.input)
    if ( is_fasta_input || input.toString().toLowerCase().endsWith("tsv") ) {
        file("${params.outdir}/input").mkdir()
        input.copyTo("${params.outdir}/input")
    }
    //Save metadata in results folder
    if ( params.metadata ) { 
        file("${params.outdir}/input").mkdir()
        file("${params.metadata}").copyTo("${params.outdir}/input")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
