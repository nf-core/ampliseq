////////////////////////////////////////////////////
/* --         LOCAL PARAMETER VALUES           -- */
////////////////////////////////////////////////////

params.summary_params = [:]

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

/*
 * Define pipeline steps
 */
params.onlyDenoising = false
if (params.onlyDenoising) {
	params.skip_abundance_tables = true
	params.skip_barplot = true
	params.skip_taxonomy = true
	params.skip_alpha_rarefaction = true
	params.skip_diversity_indices = true
	params.skip_ancom = true
} else {
	params.skip_abundance_tables = false
	params.skip_barplot = false
	params.skip_taxonomy = false
	params.skip_alpha_rarefaction = false
	params.skip_diversity_indices = false
	params.skip_ancom = false
}

/*
 * Import input files
 */
if (params.metadata) {
	ch_metadata = Channel.fromPath("${params.metadata}", checkIfExists: true)
		//.into { ch_metadata_for_barplot; ch_metadata_for_alphararefaction; ch_metadata_for_diversity_core; ch_metadata_for_alpha_diversity; ch_metadata_for_metadata_category_all; ch_metadata_for_metadata_category_pairwise; ch_metadata_for_beta_diversity; ch_metadata_for_beta_diversity_ordination; ch_metadata_for_ancom; ch_metadata_for_ancom_tax; ch_metadata_for_ancom_asv }
} else {
	Channel.from()
		.into { ch_metadata_for_barplot; ch_metadata_for_alphararefaction; ch_metadata_for_diversity_core; ch_metadata_for_alpha_diversity; ch_metadata_for_metadata_category_all; ch_metadata_for_metadata_category_pairwise; ch_metadata_for_beta_diversity; ch_metadata_for_beta_diversity_ordination; ch_metadata_for_ancom; ch_metadata_for_ancom_tax; ch_metadata_for_ancom_asv }
}

if (params.classifier) {
	Channel.fromPath("${params.classifier}", checkIfExists: true)
		   .set { ch_qiime_classifier }
}

if (params.dada_ref_taxonomy && !params.onlyDenoising) {
	Channel.fromPath("${params.dada_ref_taxonomy}", checkIfExists: true)
		   .set { ch_dada_ref_taxonomy }
} else { ch_dada_ref_taxonomy = Channel.empty() }

/*
 * Sanity check input values
 */

if (!params.FW_primer) { exit 1, "Option --FW_primer missing" }
if (!params.RV_primer) { exit 1, "Option --RV_primer missing" }
if (!params.input) { exit 1, "Option --input missing" }

if ("${params.split}".indexOf("_") > -1 ) {
	exit 1, "Underscore is not allowed in --split, please review your input."
}

single_end = params.single_end
if (params.pacbio) {
   single_end = true
}

//TRUE, FALSE, pseudo allowed, see https://benjjneb.github.io/dada2/pseudo.html#Pseudo-pooling
if (!["pooled", "independent", "pseudo"].contains(params.sample_inference)) {
	exit 1, "Please set --sample_inference to one of the following:\n\t-\"independent\" (lowest sensitivity and lowest resources),\n\t-\"pseudo\" (balance between required resources and sensitivity),\n\t-\"pooled\" (highest sensitivity and resources)."
}

if (params.double_primer && params.retain_untrimmed) { 
	exit 1, "Incompatible parameters --double_primer and --retain_untrimmed cannot be set at the same time."
}

if (!params.classifier){
        if (!(params.taxon_reference == 'silva' || params.taxon_reference == 'unite')) exit 1, "--taxon_reference need to be set to either 'silva' or 'unite'"
}

////////////////////////////////////////////////////
/* --          CONFIG FILES                    -- */
////////////////////////////////////////////////////

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

def dada_filter_and_trim_options = modules['dada_filter_and_trim']
dada_filter_and_trim_options.args       += single_end ? ", maxEE = $params.maxEE" : ", maxEE = c($params.maxEE, $params.maxEE)"
if (params.pacbio) {
	//PacBio data
	//TODO: "truncLen = $params.trunclenf" is in ampliseq 1.2.0 "$trunc" which is pruduced for "single_end" here: https://github.com/nf-core/ampliseq/blob/1e76f1c0c6d270572573a6ad76aa40bfc185c3da/main.nf#L944-L962
	dada_filter_and_trim_options.args   +=", truncQ = 2, truncLen = $params.trunclenf, minLen = $params.minLen, maxLen = $params.maxLen, rm.phix = FALSE"
} else if (params.illumina_its) {
	//Illumina ITS data or other sequences with high length variability
	dada_filter_and_trim_options.args   += ", truncQ = 2, minLen = $params.minLen, rm.phix = TRUE"
} else {
	//Illumina 16S data
	dada_filter_and_trim_options.args   += single_end ? ", truncLen = $params.trunclenf" : ", truncLen = c($params.trunclenf, $params.trunclenr)"
	dada_filter_and_trim_options.args   += ", truncQ = 2, rm.phix = TRUE"
}

//TODO: adjust for single_end & PacBio
def dada_error_model_options = modules['dada_error_model']
dada_error_model_options.args       += ", nbases = 1e8" //default is 1e8
if (params.pacbio) {
	//PacBio data
	dada_error_model_options.args   +=", errorEstimationFunction = PacBioErrfun"
} else {
	//Illumina data
	dada_error_model_options.args   +=", errorEstimationFunction = loessErrfun"
}

def dada_denoising_options = modules['dada_denoising']
dada_denoising_options.args         += params.sample_inference == "pseudo" ? ", pool = \"pseudo\"" : params.sample_inference == "pooled" ? ", pool = TRUE" : ", pool = FALSE"

def dada_chimera_removal_options = modules['dada_chimera_removal']

def multiqc_options         = modules['multiqc']
multiqc_options.args       += params.multiqc_title ? " --title \"$params.multiqc_title\"" : ''

def dada_assign_taxonomy_options  = modules['dada_assign_taxonomy']
dada_assign_taxonomy_options.args += ", tryRC = TRUE"

include { RENAME_RAW_DATA_FILES              } from './modules/local/process/rename_raw_data_files'
include { DADA_FILTER_AND_TRIM               } from './modules/local/process/dada'                        addParams( options: dada_filter_and_trim_options            )
include { DADA_ERROR_MODEL                   } from './modules/local/process/dada'                        addParams( options: dada_error_model_options                )
include { DADA_DEREPLICATE                   } from './modules/local/process/dada'
include { DADA_DENOISING                     } from './modules/local/process/dada'                        addParams( options: dada_denoising_options                  )
include { DADA_CHIMERA_REMOVAL               } from './modules/local/process/dada'                        addParams( options: dada_chimera_removal_options            )
include { DADA_STATS                         } from './modules/local/process/dada'
include { DADA_MERGE_AND_PUBLISH             } from './modules/local/process/dada'
include { DADA_ASSIGN_TAXONOMY               } from './modules/local/process/dada'                        addParams( options: dada_assign_taxonomy_options            )
include { MULTIQC                            } from './modules/local/process/multiqc'                     addParams( options: multiqc_options                                             )
include { GET_SOFTWARE_VERSIONS              } from './modules/local/process/get_software_versions'       addParams( options: [publish_files : ['csv':'']]                                )

/*
 * SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
 */

 //TODO

 ////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////

/*
 * MODULE: Installed directly from nf-core/modules
 */
def fastqc_options              = modules['fastqc']

//TODO: add another option: params.illumina_its which should be "-g ${params.FW_primer} -a ${params.RV_primer_RevComp} -G ${params.RV_primer} -A ${params.FW_primer_RevComp} -n 2"
//TODO: probably when params.illumina_its than also params.retain_untrimmed, have to test
def cutadapt_options_args       = !single_end ? " -g ${params.FW_primer} -G ${params.RV_primer}" : params.pacbio ? " --rc -g ${params.FW_primer}...${params.RV_primer}" : " -g ${params.FW_primer}"
def cutadapt_options 			= modules['cutadapt']
cutadapt_options.args          += cutadapt_options_args
cutadapt_options.args          += params.retain_untrimmed ? '' : " --discard-untrimmed"
def cutadapt_2nd_options        = [:]
cutadapt_2nd_options.args       = cutadapt_options_args
cutadapt_2nd_options.args      += " --discard-trimmed"
cutadapt_2nd_options.suffix     = "_double-primer"

//include { MULTIQC } from './modules/nf-core/software/multiqc/main' addParams( options: multiqc_options    )
include { FASTQC } from './modules/nf-core/software/fastqc/main' addParams( options: fastqc_options    )
include { CUTADAPT } from './modules/nf-core/software/cutadapt/main' addParams( options: cutadapt_options    )
include { CUTADAPT as CUTADAPT_2ND } from './modules/nf-core/software/cutadapt/main' addParams( options: cutadapt_2nd_options    )

/*
 * SUBWORKFLOW: Consisting entirely of nf-core/modules
 */

 ////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

// Info required for completion email and summary
def multiqc_report      = []

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def parse_samplesheet(LinkedHashMap row) {
	//Check if manifest contains column sampleID  & forwardReads
	if (row.sampleID == null || row.forwardReads == null) { 
		exit 1, "ERROR: Please check input samplesheet -> Column 'sampleID' and 'forwardReads' are required but not detected."
	}
	//Check if manifest contains a column for reverse reads
	if (row.reverseReads == null && !single_end) { 
		exit 1, "ERROR: Please check input samplesheet -> Column 'reverseReads' is missing. In case you do have only single ended reads, please specify '--single_end' or '--pacbio'."
	}
	//read meta info
    def meta = [:]
    meta.id           = row.sampleID
    meta.single_end   = single_end.toBoolean()
	meta.run          = row.run == null ? "1" : row.run
	//read data info
    def array = []
    if (!file(row.forwardReads).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Forward read FastQ file does not exist!\n${row.forwardReads}"
    }
    if (meta.single_end) {
        array = [ meta, [ file(row.forwardReads) ] ]
    } else {
		if (!file(row.reverseReads).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Reverse read FastQ file does not exist!\n${row.reverseReads}"
        }
        array = [ meta, [ file(row.forwardReads), file(row.reverseReads) ] ]
    }
    return array
}

workflow AMPLISEQ {

	/*
	* Create a channel for input read files
	*/
	if ( params.input.toString().toLowerCase().endsWith("tsv") ) {
		// Sample sheet input

		tsvFile = file(params.input).getName()
		// extracts read files from TSV and distribute into channels
		Channel
			.fromPath(params.input)
			.ifEmpty {exit 1, log.info "Cannot find path file ${tsvFile}"}
			.splitCsv(header:true, sep:'\t')
			.map { parse_samplesheet(it) }
			.set { ch_reads }
	} else {
		// Folder input

		//Check folders in folder when multipleSequencingRuns
		folders = params.multipleSequencingRuns ? "/*" : ""
		if ( single_end ) {
			//Get files - single end
			Channel
				.fromPath( params.input + folders + params.extension )
				.ifEmpty { exit 1, "Cannot find any reads matching: \"${params.input}${params.extension}\"\nPlease revise the input folder (\"--input\"): \"${params.input}\"\nand the input file pattern (\"--extension\"): \"${params.extension}\"\nIf you have multiple sequencing runs, please add \"--multipleSequencingRuns\".\nNB: Path needs to be enclosed in quotes!" }
				.map { read ->
						def meta = [:]
						meta.id           = read.baseName.toString().indexOf("_") != -1 ? read.baseName.toString().take(read.baseName.toString().indexOf("_")) : read.baseName
						meta.single_end   = single_end.toBoolean()
						meta.run          = params.multipleSequencingRuns ? read.take(read.findLastIndexOf{"/"})[-1] : "1"
						[ meta, read ] }
				.set { ch_reads }
		} else {
			//Get files - paired end
			Channel
				.fromFilePairs( params.input + folders + params.extension, size: 2 )
				.ifEmpty { exit 1, "Cannot find any reads matching: \"${params.input}${params.extension}\"\nPlease revise the input folder (\"--input\"): \"${params.input}\"\nand the input file pattern (\"--extension\"): \"${params.extension}\"\nIf you have multiple sequencing runs, please add \"--multipleSequencingRuns\".\nNB: Path needs to be enclosed in quotes!" }
				.map { name, reads ->
						def meta = [:]
						meta.id           = name.toString().indexOf("_") != -1 ? name.toString().take(name.toString().indexOf("_")) : name
						meta.single_end   = single_end.toBoolean()
						meta.run          = params.multipleSequencingRuns ? reads[0].take(reads[0].findLastIndexOf{"/"})[-1] : "1"
						[ meta, reads ] }
				.set { ch_reads }
		}
		if (params.multipleSequencingRuns) {
			//Get folder information
			ch_reads
				.flatMap { meta, reads -> [ meta.run ] }
				.unique()
				.set { ch_folders }
			//Report folders with sequencing files
			ch_folders
				.collect()
				.subscribe {
					String folders = it.toString().replace("[", "").replace("]","") 
					log.info "\nFound the folder(s) \"$folders\" containing sequencing read files matching \"${params.extension}\" in \"${params.input}\".\n" }
			//Stop if folder count is 1 and params.multipleSequencingRuns
			ch_folders
				.count()
				.subscribe { if ( it == 1 ) exit 1, "Found only one folder with read data but \"--multipleSequencingRuns\" was specified. Please review data input." }
			//Stop if folder names contain "_" or "${params.split}"
			//TODO: this might be not neccessary if params.split is removed
			ch_folders
				.subscribe { 
					if ( it.toString().indexOf("${params.split}") > -1 ) exit 1, "Folder name \"$it\" contains \"${params.split}\", but may not. Please review data input or choose another string using \"--split [str]\" (no underscore allowed!)."
					if ( it.toString().indexOf("_") > -1 ) exit 1, "Folder name \"$it\" contains \"_\", but may not. Please review data input." 
				}
		}
	}

    /*
     * MODULE: Rename files
     */
	//Check whether all sampleID = meta.id are unique
	//TODO: if not all sampleID unique, rename files with meta.run? //if ( params.multipleSequencingRuns ) { "meta.id" = "$meta.run${params.split}$meta.id" }
	//TODO: params.split might not be neccessary any more, decide when integrating QIIME2
	ch_reads
		.map { meta, reads -> [ meta.id ] }
		.count()
		.set { ch_ids }
	ch_reads
		.map { meta, reads -> [ meta.id ] }
		.unique()
		.count()
		.mix( ch_ids )
		.collect()
		.subscribe { k = it[0]; n = it[1];
			if ( k != n ) exit 1, "Please review data input, sampleIDs ($k) are not unique ($n).";
			}

	//rename files as in FASTQC to get same sample names everywhere!
	RENAME_RAW_DATA_FILES ( ch_reads )

    /*
     * MODULE: FastQC
     */
    if (!params.skip_fastqc) {
        FASTQC ( RENAME_RAW_DATA_FILES.out ).html.set { fastqc_html }
    }
    ch_software_versions = Channel.empty()
    ch_software_versions = ch_software_versions.mix(FASTQC.out.version.first().ifEmpty(null))

    /*
     * MODULE: Cutadapt
     */
    CUTADAPT ( RENAME_RAW_DATA_FILES.out ).reads.set { ch_trimmed_reads }
    ch_software_versions = ch_software_versions.mix(CUTADAPT.out.version.first().ifEmpty(null))

	if (params.double_primer) {
		CUTADAPT_2ND ( ch_trimmed_reads ).reads.set { ch_trimmed_reads }
	}

    /*
     * SUBWORKFLOW / MODULES : ASV generation with DADA2
     */
	DADA_FILTER_AND_TRIM ( ch_trimmed_reads )
	ch_software_versions = ch_software_versions.mix(DADA_FILTER_AND_TRIM.out.version.first().ifEmpty(null))

	//group by sequencing run
	DADA_FILTER_AND_TRIM.out.reads
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
				meta.id = ids.flatten()
				[ meta, reads.flatten() ] }
		.set { ch_filt_reads }

	DADA_ERROR_MODEL ( ch_filt_reads )

	DADA_DEREPLICATE ( ch_filt_reads )

	//group by meta
	DADA_DEREPLICATE.out.dereplicated
		.join( DADA_ERROR_MODEL.out.errormodel )
		.set { ch_derep_errormodel }
	DADA_DENOISING ( ch_derep_errormodel  )

	DADA_CHIMERA_REMOVAL ( DADA_DENOISING.out.seqtab )

	//group by sequencing run & group by meta
	DADA_FILTER_AND_TRIM.out.log
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
				meta.id = ids.flatten()
				[ meta, reads.flatten() ] }
		.join( DADA_DENOISING.out.denoised )
		.join( DADA_DENOISING.out.mergers )
		.join( DADA_CHIMERA_REMOVAL.out.rds )
		.set { ch_track_numbers }
	DADA_STATS ( ch_track_numbers )

	//merge if several runs, otherwise just publish
	DADA_MERGE_AND_PUBLISH ( 
		DADA_STATS.out.stats.map { meta, stats -> stats }.collect(), 
		DADA_CHIMERA_REMOVAL.out.rds.map { meta, rds -> rds }.collect() )

	//taxonomic classification with dada2
	DADA_ASSIGN_TAXONOMY ( DADA_MERGE_AND_PUBLISH.out.rds, ch_dada_ref_taxonomy )


    /*
     * MODULE: Pipeline reporting
     */
    GET_SOFTWARE_VERSIONS ( 
        ch_software_versions.map { it }.collect()
    )

    /*
     * MultiQC
     */
    if (!params.skip_multiqc) {
        workflow_summary    = Schema.params_summary_multiqc(workflow, params.summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        MULTIQC (
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            GET_SOFTWARE_VERSIONS.out.yaml.collect(),
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
            FASTQC.out.zip.collect{it[1]}.ifEmpty([]),
            CUTADAPT.out.log.collect{it[1]}.ifEmpty([])
        )
        multiqc_report = MULTIQC.out.report.toList()
    }

}

////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

workflow.onComplete {
    Completion.email(workflow, params, params.summary_params, projectDir, log, multiqc_report)
    Completion.summary(workflow, params, log)
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////