/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
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
if (!params.FW_primer) { exit 1, "Option --FW_primer missing" }
if (!params.RV_primer) { exit 1, "Option --RV_primer missing" }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
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

if (params.dada_ref_taxonomy && !params.skip_taxonomy) {
    ch_dada_ref_taxonomy = Channel.fromList(params.dada_ref_databases[params.dada_ref_taxonomy]["file"]).map { file(it) }
} else { ch_dada_ref_taxonomy = Channel.empty() }

if (params.qiime_ref_taxonomy && !params.skip_taxonomy && !params.classifier) {
    ch_qiime_ref_taxonomy = Channel.fromList(params.qiime_ref_databases[params.qiime_ref_taxonomy]["file"]).map { file(it) }
} else { ch_qiime_ref_taxonomy = Channel.empty() }

// Set non-params Variables

single_end = params.single_end
if (  params.pacbio || params.iontorrent ) {
   single_end = true
}

trunclenf = params.trunclenf ? params.trunclenf : 0 
trunclenr = params.trunclenr ? params.trunclenr : 0
if ( !single_end && !params.illumina_pe_its && (params.trunclenf == false || params.trunclenr == false) ) {
    find_truncation_values = true
    log.warn "No DADA2 cutoffs were specified (`--trunclenf` & --`trunclenr`), therefore reads will be truncated where median quality drops below ${params.trunc_qmin} (defined by `--trunc_qmin`) but at least a fraction of ${params.trunc_rmin} (defined by `--trunc_rmin`) of the reads will be retained.\nThe chosen cutoffs do not account for required overlap for merging, therefore DADA2 might have poor merging efficiency or even fail.\n"
} else { find_truncation_values = false }

//only run QIIME2 when taxonomy is actually calculated and all required data is available
if ( !params.enable_conda && !params.skip_taxonomy && !params.skip_qiime ) {
    run_qiime2 = true
} else { run_qiime2 = false }

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

def dada2_filtntrim_options = modules['dada2_filtntrim']
dada2_filtntrim_options.args       += single_end ? ", maxEE = $params.max_ee" : ", maxEE = c($params.max_ee, $params.max_ee)"
if (params.pacbio) {
    //PacBio data
    dada2_filtntrim_options.args   +=", trimLeft = 0, minLen = $params.min_len, maxLen = $params.max_len, rm.phix = FALSE"
} else if (params.iontorrent) {
    //Ion-torrent data
    dada2_filtntrim_options.args   += ", trimLeft = 15, minLen = $params.min_len, maxLen = $params.max_len, rm.phix = TRUE"
} else if (params.illumina_pe_its) {
    //Illumina ITS data or other sequences with high length variability
    dada2_filtntrim_options.args   += ", trimLeft = 0, minLen = $params.min_len, maxLen = $params.max_len, rm.phix = TRUE"
} else {
    //Illumina 16S data
    dada2_filtntrim_options.args   += ", trimLeft = 0, minLen = $params.min_len, maxLen = $params.max_len, rm.phix = TRUE"
}

def dada2_quality_options = modules['dada2_quality'] 

def trunclen_options = [:]
trunclen_options.args       ="$params.trunc_qmin $params.trunc_rmin"

def dada2_err_options = modules['dada2_err']
dada2_err_options.args   += params.pacbio ? ", errorEstimationFunction = PacBioErrfun" : ", errorEstimationFunction = loessErrfun"

def dada2_denoising_options = modules['dada2_denoising']
if (params.iontorrent) {
    //Ion-torrent data
    dada2_denoising_options.args   += ", BAND_SIZE = 32, HOMOPOLYMER_GAP_PENALTY = -1"
} else {
    dada2_denoising_options.args   += ", BAND_SIZE = 16, HOMOPOLYMER_GAP_PENALTY = NULL"
}
dada2_denoising_options.args         += params.sample_inference == "pseudo" ? ", pool = \"pseudo\"" : params.sample_inference == "pooled" ? ", pool = TRUE" : ", pool = FALSE"
dada2_denoising_options.args2        += params.concatenate_reads ? ", justConcatenate = TRUE" : ", justConcatenate = FALSE"

def dada2_rmchimera_options = modules['dada2_rmchimera']

def dada2_taxonomy_options  = modules['dada2_taxonomy']
dada2_taxonomy_options.args += params.pacbio ? ", tryRC = TRUE" : ""
dada2_taxonomy_options.args += params.iontorrent ? ", tryRC = TRUE" : ""

def dada2_addspecies_options  = modules['dada2_addspecies']
dada2_addspecies_options.args += params.pacbio ? ", tryRC = TRUE" : ""
dada2_addspecies_options.args += params.iontorrent ? ", tryRC = TRUE" : ""

include { RENAME_RAW_DATA_FILES         } from '../modules/local/rename_raw_data_files'
include { DADA2_FILTNTRIM               } from '../modules/local/dada2_filtntrim'              addParams( options: dada2_filtntrim_options         )
include { DADA2_QUALITY                 } from '../modules/local/dada2_quality'                addParams( options: dada2_quality_options           )
include { TRUNCLEN                      } from '../modules/local/trunclen'                     addParams( options: trunclen_options                )
include { DADA2_ERR                     } from '../modules/local/dada2_err'                    addParams( options: dada2_err_options               )
include { DADA2_DEREPLICATE             } from '../modules/local/dada2_dereplicate'            addParams( options: modules['dada2_dereplicate']    )
include { DADA2_DENOISING               } from '../modules/local/dada2_denoising'              addParams( options: dada2_denoising_options         )
include { DADA2_RMCHIMERA               } from '../modules/local/dada2_rmchimera'              addParams( options: dada2_rmchimera_options         )
include { DADA2_STATS                   } from '../modules/local/dada2_stats'                  addParams( options: modules['dada2_stats']          )
include { DADA2_MERGE                   } from '../modules/local/dada2_merge'                  addParams( options: modules['dada2_merge']          )
include { FORMAT_TAXONOMY               } from '../modules/local/format_taxonomy'
include { ITSX_CUTASV                   } from '../modules/local/itsx_cutasv'                  addParams( options: modules['itsx_cutasv']          )         
include { MERGE_STATS                   } from '../modules/local/merge_stats'                  addParams( options: modules['merge_stats']          )
include { DADA2_TAXONOMY                } from '../modules/local/dada2_taxonomy'               addParams( options: dada2_taxonomy_options          )
include { DADA2_ADDSPECIES              } from '../modules/local/dada2_addspecies'             addParams( options: dada2_addspecies_options        )
include { FORMAT_TAXRESULTS             } from '../modules/local/format_taxresults'
include { QIIME2_INSEQ                  } from '../modules/local/qiime2_inseq'                 addParams( options: modules['qiime2_inseq']         )
include { QIIME2_FILTERTAXA             } from '../modules/local/qiime2_filtertaxa'            addParams( options: modules['qiime2_filtertaxa']    )
include { QIIME2_INASV                  } from '../modules/local/qiime2_inasv'                 addParams( options: modules['qiime2_inasv']         )
include { FILTER_STATS                  } from '../modules/local/filter_stats'                 addParams( options: modules['filter_stats']         )
include { MERGE_STATS as MERGE_STATS_FILTERTAXA } from '../modules/local/merge_stats'          addParams( options: modules['merge_stats']          )
include { QIIME2_BARPLOT                } from '../modules/local/qiime2_barplot'               addParams( options: modules['qiime2_barplot']       )
include { METADATA_ALL                  } from '../modules/local/metadata_all'
include { METADATA_PAIRWISE             } from '../modules/local/metadata_pairwise'
include { QIIME2_INTAX                  } from '../modules/local/qiime2_intax'                 addParams( options: modules['qiime2_intax']         )
include { GET_SOFTWARE_VERSIONS         } from '../modules/local/get_software_versions'        addParams( options: [publish_files : ['tsv':'']]    )

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

if (params.pacbio) {
    //PacBio data
    cutadapt_options_args       = " --rc -g ${params.FW_primer}...${params.RV_primer}"
} else if (params.iontorrent) {
    //IonTorrent data
    cutadapt_options_args       = " --rc -g ${params.FW_primer}...${params.RV_primer}"
} else if (params.single_end) {
    //Illumina SE
    cutadapt_options_args       = " -g ${params.FW_primer}"
} else {
    //Illumina PE
    cutadapt_options_args       = " -g ${params.FW_primer} -G ${params.RV_primer}"
}

def cutadapt_options 			= modules['cutadapt']
cutadapt_options.args          += cutadapt_options_args
cutadapt_options.args          += params.retain_untrimmed ? '' : " --discard-untrimmed"

//prepare reverse complement primers to remove those in read-throughs
// Get the complement of a DNA sequence
// Complement table taken from http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html
def make_complement(String seq) {
    def complements = [ A:'T', T:'A', U:'A', G:'C', C:'G', Y:'R', R:'Y', S:'S', W:'W', K:'M', M:'K', B:'V', D:'H', H:'D', V:'B', N:'N' ]
    comp = seq.toUpperCase().collect { base -> complements[ base ] ?: 'X' }.join()
    return comp
}
FW_primer_RevComp = make_complement ( "${params.FW_primer}".reverse() )
RV_primer_RevComp = make_complement ( "${params.RV_primer}".reverse() )
def cutadapt_readthrough_options      = modules['cutadapt_readthrough']
cutadapt_readthrough_options.args    += " -a ${RV_primer_RevComp} -A ${FW_primer_RevComp}"

def cutadapt_doubleprimer_options        = modules['cutadapt_doubleprimer']
cutadapt_doubleprimer_options.args      += cutadapt_options_args

include { PARSE_INPUT                   } from '../subworkflows/local/parse_input'
include { QIIME2_PREPTAX                } from '../subworkflows/local/qiime2_preptax'           addParams( options: modules['qiime2_preptax']       )
include { QIIME2_TAXONOMY               } from '../subworkflows/local/qiime2_taxonomy'          addParams( options: modules['qiime2_taxonomy']      )
include { CUTADAPT_WORKFLOW             } from '../subworkflows/local/cutadapt_workflow'        addParams( standard_options: cutadapt_options, readthrough_options: cutadapt_readthrough_options,doubleprimer_options: cutadapt_doubleprimer_options,summary_options: modules['cutadapt_summary'],summary_merge_options: modules['cutadapt_summary_merge'] )
include { QIIME2_EXPORT                 } from '../subworkflows/local/qiime2_export'            addParams( absolute_options: modules['qiime2_export_absolute'], relasv_options: modules['qiime2_export_relasv'],reltax_options: modules['qiime2_export_reltax'],combine_table_options: modules['combine_table'] )
include { QIIME2_DIVERSITY              } from '../subworkflows/local/qiime2_diversity'         addParams( tree_options: modules['qiime2_tree'], alphararefaction_options: modules['qiime2_alphararefaction'], diversity_core_options: modules['qiime2_diversity_core'], diversity_alpha_options: modules['qiime2_diversity_alpha'], diversity_beta_options: modules['qiime2_diversity_beta'], diversity_betaord_options: modules['qiime2_diversity_betaord'] )
include { QIIME2_ANCOM                  } from '../subworkflows/local/qiime2_ancom'             addParams( filterasv_options: modules['qiime2_filterasv'], ancom_tax_options: modules['qiime2_ancom_tax'], ancom_asv_options: modules['qiime2_ancom_asv'] )

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
def fastqc_options              = modules['fastqc']

def cutadapt_taxonomy_options   = modules['cutadapt_taxonomy']
cutadapt_taxonomy_options.args += " -g ${params.FW_primer}...${RV_primer_RevComp}"

def multiqc_options         = modules['multiqc']
multiqc_options.args       += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

include { CUTADAPT as CUTADAPT_TAXONOMY     } from '../modules/nf-core/modules/cutadapt/main' addParams( options: cutadapt_taxonomy_options     )
include { FASTQC                            } from '../modules/nf-core/modules/fastqc/main'   addParams( options: fastqc_options                )
include { MULTIQC                           } from '../modules/nf-core/modules/multiqc/main'  addParams( options: multiqc_options               )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report      = []

workflow AMPLISEQ {
    ch_software_versions = Channel.empty()

    //
    // Create a channel for input read files
    //
    PARSE_INPUT ( params.input, single_end, params.multiple_sequencing_runs, params.extension )
    ch_reads = PARSE_INPUT.out.reads
    ch_fasta = PARSE_INPUT.out.fasta

    //
    // MODULE: Rename files
    //
    RENAME_RAW_DATA_FILES ( ch_reads )

    //
    // MODULE: Run FastQC
    //
    if (!params.skip_fastqc) {
        FASTQC ( RENAME_RAW_DATA_FILES.out ).html.set { fastqc_html }
        ch_software_versions = ch_software_versions.mix(FASTQC.out.version.first().ifEmpty(null))
    }

    //
    // MODULE: Cutadapt
    //
    CUTADAPT_WORKFLOW ( 
        RENAME_RAW_DATA_FILES.out,
        params.illumina_pe_its,
        params.double_primer
    ).reads.set { ch_trimmed_reads }
    ch_software_versions = ch_software_versions.mix(CUTADAPT_WORKFLOW.out.version.first().ifEmpty(null))

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
    DADA2_QUALITY ( ch_all_trimmed_reads )
    
    //find truncation values in case they are not supplied
    if ( find_truncation_values ) {
        TRUNCLEN ( DADA2_QUALITY.out.tsv )
        TRUNCLEN.out
            .toSortedList()
            .set { ch_trunc }
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
    DADA2_FILTNTRIM ( ch_trimmed_reads )
    ch_software_versions = ch_software_versions.mix(DADA2_FILTNTRIM.out.version.first().ifEmpty(null))

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

    DADA2_DEREPLICATE ( ch_filt_reads )

    //group by meta
    DADA2_DEREPLICATE.out.dereplicated
        .join( DADA2_ERR.out.errormodel )
        .set { ch_derep_errormodel }
    DADA2_DENOISING ( ch_derep_errormodel  )

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
    MERGE_STATS (CUTADAPT_WORKFLOW.out.summary, DADA2_MERGE.out.dada2stats)

    //
    // SUBWORKFLOW / MODULES : Taxonomic classification with DADA2 and/or QIIME2
    //
    //Alternative entry point for fasta that is being classified - the if clause needs to be the opposite (i.e. with !) of that in subworkflow/local/parse.nf
    if ( !(params.input.toString().toLowerCase().endsWith(".fasta") || params.input.toString().toLowerCase().endsWith(".fna") || params.input.toString().toLowerCase().endsWith(".fa") )) {
        ch_fasta = DADA2_MERGE.out.fasta
    }

    //DADA2
    if (!params.skip_taxonomy) {
        FORMAT_TAXONOMY ( ch_dada_ref_taxonomy.collect() )
        ch_assigntax = FORMAT_TAXONOMY.out.assigntax
        ch_addspecies = FORMAT_TAXONOMY.out.addspecies
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
        if (!params.cut_its) {
            DADA2_TAXONOMY ( ch_fasta, ch_assigntax, 'ASV_tax.tsv' )
            DADA2_ADDSPECIES ( DADA2_TAXONOMY.out.rds, ch_addspecies, 'ASV_tax_species.tsv' )
            ch_dada2_tax = DADA2_ADDSPECIES.out.tsv
        //Cut out ITS region if long ITS reads
        } else {
            ITSX_CUTASV ( ch_fasta )
            ch_software_versions = ch_software_versions.mix(ITSX_CUTASV.out.version.ifEmpty(null))
            ch_cut_fasta = ITSX_CUTASV.out.fasta
            DADA2_TAXONOMY ( ch_cut_fasta, ch_assigntax, 'ASV_ITS_tax.tsv' )
            DADA2_ADDSPECIES ( DADA2_TAXONOMY.out.rds, ch_addspecies, 'ASV_ITS_tax_species.tsv' )
            FORMAT_TAXRESULTS ( DADA2_TAXONOMY.out.tsv, DADA2_ADDSPECIES.out.tsv, ch_fasta )
            ch_dada2_tax = FORMAT_TAXRESULTS.out.tsv
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
        ch_software_versions = ch_software_versions.mix( QIIME2_TAXONOMY.out.version.ifEmpty(null) ) //usually a .first() is here, dont know why this leads here to a warning
    }
    //
    // SUBWORKFLOW / MODULES : Downstream analysis with QIIME2
    //
    if ( run_qiime2 ) {
        //Import ASV abundance table and sequences into QIIME2
        QIIME2_INASV ( DADA2_MERGE.out.asv )
        QIIME2_INSEQ ( ch_fasta )

        //Import taxonomic classification into QIIME2, if available
        if ( params.skip_taxonomy ) {
            log.info "Skip taxonomy classification"
            ch_tax = Channel.empty()
        } else if ( params.dada_ref_taxonomy ) {
            log.info "Use DADA2 taxonomy classification"
            ch_tax = QIIME2_INTAX ( ch_dada2_tax ).qza
        } else if ( params.qiime_ref_taxonomy || params.classifier ) {
            log.info "Use QIIME2 taxonomy classification"
            ch_tax = QIIME2_TAXONOMY.out.qza
        } else { 
            log.info "Use no taxonomy classification"
            ch_tax = Channel.empty() 
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
            FILTER_STATS ( DADA2_MERGE.out.asv, QIIME2_FILTERTAXA.out.tsv )
            MERGE_STATS_FILTERTAXA (MERGE_STATS.out.tsv, FILTER_STATS.out.tsv)
            ch_asv = QIIME2_FILTERTAXA.out.asv
            ch_seq = QIIME2_FILTERTAXA.out.seq
        } else {
            ch_asv = QIIME2_INASV.out.qza
            ch_seq = QIIME2_INSEQ.out.qza
        }
        //Export various ASV tables
        if (!params.skip_abundance_tables) {
            QIIME2_EXPORT ( ch_asv, ch_seq, ch_tax, QIIME2_TAXONOMY.out.tsv )
        }

        if (!params.skip_barplot) {
            QIIME2_BARPLOT ( ch_metadata, ch_asv, ch_tax )
        }

        //Select metadata categories for diversity analysis & ancom
        if (!params.skip_ancom || !params.skip_diversity_indices) {
            METADATA_ALL ( ch_metadata, params.metadata_category ).set { ch_metacolumn_all }
            //return empty channel if no appropriate column was found
            ch_metacolumn_all.branch { passed: it != "" }.set { result }
            ch_metacolumn_all = result.passed
    
            METADATA_PAIRWISE ( ch_metadata ).set { ch_metacolumn_pairwise }
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
                QIIME2_FILTERTAXA.out.tsv,
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
                ch_tax
            )
        }
    }

    //
    // MODULE: Pipeline reporting
    //
    ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }

    GET_SOFTWARE_VERSIONS (
        ch_software_versions.map { it }.collect()
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
        ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

        MULTIQC (
            ch_multiqc_files.collect()
        )
        multiqc_report       = MULTIQC.out.report.toList()
        ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
    }

}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/