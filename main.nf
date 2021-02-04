#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/ampliseq
========================================================================================
 nf-core/ampliseq Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/ampliseq
----------------------------------------------------------------------------------------
*/

def helpMessage() {
	log.info nfcoreHeader()
	log.info"""
	Usage:

	The minimal command for running the pipeline is as follows:
	nextflow run nf-core/ampliseq -profile singularity --input "data" --FW_primer GTGYCAGCMGCCGCGGTAA --RV_primer GGACTACNVGGGTWTCTAAT

	In case of a timezone error, please specify "--qiime_timezone", e.g. --qiime_timezone 'Europe/Berlin'!

	Main arguments:
	  -profile [strings]            Use this parameter to choose a configuration profile. If not specified, runs locally and expects all software
	                                to be installed and available on the `PATH`. Otherwise specify a container engine, "docker" or "singularity" 
	                                and a specialized profile such as "binac".
	  --input [path/to/folder]      Folder containing paired-end demultiplexed fastq files
	                                Note: All samples have to be sequenced in one run, otherwise also specifiy "--multipleSequencingRuns"
	  --FW_primer [str]             Forward primer sequence
	  --RV_primer [str]             Reverse primer sequence
	  --metadata [path/to/file]     Path to metadata sheet, when missing most downstream analysis are skipped (barplots, PCoA plots, ...). 
	                                File extension is not relevant. Must have a comma separated list of metadata column headers.
	  --manifest [path/to/file]     Path to manifest.tsv table with the following labels in this exact order: sampleID, forwardReads, reverseReads. In case of single end reads, the labels should be: sampleID, Reads.
	                                Tab ('\t') must be the table separator. Multiple sequencing runs not supported by manifest at this stage.
	                                Default is FALSE. 
	  --qiime_timezone [str]        Needs to be specified to resolve a timezone error (default: 'Europe/Berlin')

	Other input options:
	  --extension [str]             Naming of sequencing files (default: "/*_R{1,2}_001.fastq.gz"). 
	                                The prepended "/" is required, also one "*" is required for sample names and "{1,2}" indicates read orientation
	  --multipleSequencingRuns      If samples were sequenced in multiple sequencing runs. Expects one subfolder per sequencing run
	                                in the folder specified by "--input" containing sequencing data of the specific run. These folders 
	                                may not contain underscores.
	  --split [str]                 A string that will be used between the prepended run/folder name and the sample name. (default: "-")
	                                May not be present in run/folder names and no underscore(s) allowed. Only used with "--multipleSequencingRuns"
	  --pacbio                      If PacBio data. Use this option together with --manifest.
	  --phred64                     If the sequencing data has PHRED 64 encoded quality scores (default: PHRED 33)

	Filters:
	  --exclude_taxa [str]          Comma separated list of unwanted taxa (default: "mitochondria,chloroplast")
	                                To skip taxa filtering use "none"
	  --min_frequency [int]         Remove entries from the feature table below an absolute abundance threshold (default: 1)
	  --min_samples [int]           Filtering low prevalent features from the feature table (default: 1)                   

	Cutoffs:
	  --double_primer               Cutdapt will be run twice, first to remove reads without primers (default), then a second time to remove reads that erroneously contain a second set of primers, not to be used with "--retain_untrimmed"
	  --retain_untrimmed            Cutadapt will retain untrimmed reads
	  --maxEE [number]              DADA2 read filtering option. After truncation, reads with higher than ‘maxEE’ "expected errors" will be discarded. We recommend (to start with) a value corresponding to approximately 1 expected error per 100-200 bp (default: 2)
	  --maxLen [int]                DADA2 read filtering option [PacBio only], remove reads with length greater than maxLen after trimming and truncation (default: 2999)
	  --minLen [int]                DADA2 read filtering option [PacBio only], remove reads with length less than minLen after trimming and truncation (default: 50)
	  --trunclenf [int]             DADA2 read truncation value for forward strand and single end reads, set this to 0 for no truncation
	  --trunclenr [int]             DADA2 read truncation value for reverse strand, set this to 0 for no truncation
	  --trunc_qmin [int]            If --trunclenf and --trunclenr are not set, these values will be automatically determined using this mean quality score (not preferred) (default: 25)
	  --trunc_rmin [float]          Assures that values chosen with --trunc_qmin will retain a fraction of reads (default: 0.75)

	References:                     If you have trained a compatible classifier before, or want to use a custom database
	  --classifier [path/to/file]   Path to QIIME2 classifier file (typically *-classifier.qza)
	  --classifier_removeHash       Remove all hash signs from taxonomy strings, resolves a rare ValueError during classification (process classifier)
	  --reference_database          Path to file with reference database with taxonomies, currently either a qiime compatible file Silva_132_release.zip, or a UNITE fasta file (default: "https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_132_release.zip")
	  --taxon_reference             Specify which database to use for taxonomic assignment. Either 'silva' or 'unite' (default: 'silva')

	Statistics:
	  --metadata_category [str]     Comma separated list of metadata column headers for statistics (default: false)
	                                If not specified, all suitable columns in the metadata sheet will be used.
	                                Suitable are columns which are categorical (not numerical) and have multiple  
	                                different values that are not all unique.

	Other options:
	  --untilQ2import               Skip all steps after importing into QIIME2, used for visually choosing DADA2 parameter
	  --Q2imported [path/to/file]   Path to imported reads (e.g. "demux.qza"), used after visually choosing DADA2 parameter
	  --onlyDenoising               Skip all steps after denoising, produce only sequences and abundance tables on ASV level
	  --keepIntermediates           Keep additional intermediate files, such as trimmed reads or various QIIME2 archives
	  --outdir [file]               The output directory where the results will be saved
	  --publish_dir_mode [str]      Mode for publishing results in the output directory. Available: symlink, rellink, link, copy, copyNoFollow, move (Default: copy)
	  --email [email]               Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
	  --email_on_fail [email]       Same as --email, except only send mail if the workflow is not successful
	  --max_multiqc_email_size [str]  Threshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB)
	  -name [str]                   Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

	Skipping steps:
	  --skip_fastqc                 Skip FastQC
	  --skip_alpha_rarefaction      Skip alpha rarefaction
	  --skip_taxonomy               Skip taxonomic classification
	  --skip_barplot                Skip producing barplot
	  --skip_abundance_tables       Skip producing any relative abundance tables
	  --skip_diversity_indices      Skip alpha and beta diversity analysis
	  --skip_ancom                  Skip differential abundance testing

	AWSBatch options:
	  --awsqueue [str]                The AWSBatch JobQueue that needs to be set when running on AWSBatch
	  --awsregion [str]               The AWS Region for your AWS Batch job to run on
	  --awscli [str]                  Path to the AWS CLI tool
	""".stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
if (params.help){
	helpMessage()
	exit 0
}

// Configurable variables
params.name = false
params.email = false
params.plaintext_email = false

ch_output_docs = Channel.fromPath("$projectDir/docs/output.md")
Channel.fromPath("$projectDir/assets/matplotlibrc")
	.into { ch_mpl_for_classifier_extract_seq; ch_mpl_for_classifier_train; ch_mpl_for_qiime_import; ch_mpl_for_ancom_asv; ch_mpl_for_ancom_tax; ch_mpl_for_ancom; ch_mpl_for_beta_diversity_ord; ch_mpl_for_beta_diversity; ch_mpl_for_alpha_diversity; ch_mpl_for_metadata_pair; ch_mpl_for_metadata_cat; ch_mpl_for_diversity_core; ch_mpl_for_alpha_rare; ch_mpl_for_tree; ch_mpl_for_barcode; ch_mpl_for_relreducetaxa; ch_mpl_for_relasv; ch_mpl_for_export_dada_output; ch_mpl_filter_taxa; ch_mpl_classifier; ch_mpl_dada; ch_mpl_dada_merge; ch_mpl_for_demux_visualize; ch_mpl_for_classifier }


/*
 * Define pipeline steps
 */
params.untilQ2import = false

params.Q2imported = false
if (params.Q2imported) {
	params.skip_fastqc = true
	params.skip_multiqc = true
} else {
	params.skip_multiqc = false
}

params.onlyDenoising = false
if (params.onlyDenoising || params.untilQ2import) {
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

params.manifest = false

/*
 * Import input files
 */
if (params.metadata) {
	Channel.fromPath("${params.metadata}", checkIfExists: true)
		.into { ch_metadata_for_barplot; ch_metadata_for_alphararefaction; ch_metadata_for_diversity_core; ch_metadata_for_alpha_diversity; ch_metadata_for_metadata_category_all; ch_metadata_for_metadata_category_pairwise; ch_metadata_for_beta_diversity; ch_metadata_for_beta_diversity_ordination; ch_metadata_for_ancom; ch_metadata_for_ancom_tax; ch_metadata_for_ancom_asv }
} else {
	Channel.from()
		.into { ch_metadata_for_barplot; ch_metadata_for_alphararefaction; ch_metadata_for_diversity_core; ch_metadata_for_alpha_diversity; ch_metadata_for_metadata_category_all; ch_metadata_for_metadata_category_pairwise; ch_metadata_for_beta_diversity; ch_metadata_for_beta_diversity_ordination; ch_metadata_for_ancom; ch_metadata_for_ancom_tax; ch_metadata_for_ancom_asv }
}

if (params.Q2imported) {
	Channel.fromPath("${params.Q2imported}", checkIfExists: true)
		   .into { ch_qiime_demux_import; ch_qiime_demux_vis; ch_qiime_demux_dada }
}

if (params.classifier) {
	Channel.fromPath("${params.classifier}", checkIfExists: true)
		   .set { ch_qiime_classifier }
}

/*
 * Sanity check input values
 */
if (!params.Q2imported) { 
	if (!params.FW_primer) { exit 1, "Option --FW_primer missing" }
	if (!params.RV_primer) { exit 1, "Option --RV_primer missing" }
	if (!params.input) { exit 1, "Option --input missing" }
}

if (params.Q2imported && params.untilQ2import) {
	exit 1, "Choose either to import data into a QIIME2 artefact and quit with --untilQ2import or use an already existing QIIME2 data artefact with --Q2imported."
}

if ("${params.split}".indexOf("_") > -1 ) {
	exit 1, "Underscore is not allowed in --split, please review your input."
}

if (params.multipleSequencingRuns && params.manifest) { 
	exit 1, "The manifest file does not support multiple sequencing runs at this point."
}

single_end = false
if (params.pacbio) {
   single_end = true
}

if (single_end && !params.manifest) {
        exit 1, "A manifest file is needed for single end reads such as PacBio data."
}

if (params.double_primer && params.retain_untrimmed) { 
	exit 1, "Incompatible parameters --double_primer and --retain_untrimmed cannot be set at the same time."
}

if (!params.classifier){
        if (!(params.taxon_reference == 'silva' || params.taxon_reference == 'unite')) exit 1, "--taxon_reference need to be set to either 'silva' or 'unite'"
}

// AWSBatch sanity checking
if(workflow.profile == 'awsbatch'){
	if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
	if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
}

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

// Check AWS batch settings
if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Stage config files
ch_multiqc_config = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs = file("$projectDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$projectDir/docs/images/", checkIfExists: true)


// Header log info
log.info nfcoreHeader()
def summary = [:]
summary['Pipeline Name']  = 'nf-core/ampliseq'
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
summary['Input']            = params.manifest ?: params.input
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']   = params.awsregion
    summary['AWS Queue']    = params.awsqueue
    summary['AWS CLI']      = params.awscli
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Profile Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Profile Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config Profile URL']         = params.config_profile_url
summary['Config Files'] = workflow.configFiles.join(', ')
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
    summary['MultiQC maxsize']   = params.max_multiqc_email_size
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

if( params.trunclenf == false || params.trunclenr == false ){
	if ( !params.untilQ2import ) log.info "\n######## WARNING: No DADA2 cutoffs were specified, therefore reads will be truncated where median quality drops below ${params.trunc_qmin} but at least a fraction of ${params.trunc_rmin} of the reads will be retained.\nThe chosen cutoffs do not account for required overlap for merging, therefore DADA2 might have poor merging efficiency or even fail.\n"
}
// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-ampliseq-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/ampliseq Workflow Summary'
    section_href: 'https://github.com/nf-core/ampliseq'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.indexOf(".csv") > 0) filename
                      else null
                }

    output:
    file 'software_versions_mqc.yaml' into ch_software_versions_yaml
    file "software_versions.csv"

	script:
	"""
	echo $workflow.manifest.version > v_pipeline.txt
	echo $workflow.nextflow.version > v_nextflow.txt
	fastqc --version > v_fastqc.txt
	multiqc --version > v_multiqc.txt
	cutadapt --version > v_cutadapt.txt
	qiime --version > v_qiime.txt
	scrape_software_versions.py &> software_versions_mqc.yaml
	"""
}


if (!params.Q2imported){
	/*
	* Create a channel for optional input manifest file
	*/
	 if (params.manifest && !single_end) {
		tsvFile = file(params.manifest).getName()
		// extracts read files from TSV and distribute into channels
		Channel
			.fromPath(params.manifest)
			.ifEmpty {exit 1, log.info "Cannot find path file ${tsvFile}"}
			.splitCsv(header:true, sep:'\t')
			.map { row -> [ row.sampleID, [ file(row.forwardReads, checkIfExists: true), file(row.reverseReads, checkIfExists: true) ] ] }
			.into { ch_read_pairs; ch_read_pairs_fastqc; ch_read_pairs_name_check }
	} else if ( single_end ) {
	       	  // Manifest file is currently the only available input option for single_end
		  tsvFile = file(params.manifest).getName()
		  // extracts read files from TSV and distribute into channels
		  Channel
			.fromPath(params.manifest)
			.ifEmpty {exit 1, log.info "Cannot find path file ${tsvFile}"}
			.splitCsv(header:true, sep:'\t')
			.map { row -> [ row.sampleID, file(row.Reads, checkIfExists: true) ] }
			.into { ch_read_pairs; ch_read_pairs_fastqc; ch_read_pairs_name_check }

	/*
	* Create a channel for input read files
	*/
	} else if (params.readPaths && params.input == "data${params.extension}" && !params.multipleSequencingRuns){
		//Test input for single sequencing runs, profile = test

		Channel
			.from(params.readPaths)
			.map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
			.ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
			.map { name, reads -> [ name.toString().indexOf("_") != -1 ? name.toString().take(name.toString().indexOf("_")) : name, reads ] }
			.into { ch_read_pairs; ch_read_pairs_fastqc; ch_read_pairs_name_check }

	} else if ( !params.readPaths && params.multipleSequencingRuns ) {
		//Standard input for multiple sequencing runs

		//Get files
		Channel
			.fromFilePairs( params.input + "/*" + params.extension, size: 2 )
			.ifEmpty { exit 1, "Cannot find any reads matching: ${params.input}/*${params.extension}\nNB: Path needs to be enclosed in quotes!" }
			.into { ch_extract_folders; ch_rename_key }

		//Get folder information
		ch_extract_folders
			.flatMap { key, files -> [files[0]] }
			.map { it.take(it.findLastIndexOf{"/"})[-1] }
			.unique()
			.into { ch_count_folders; ch_check_folders; ch_report_folders }

		//Report folders with sequencing files
		ch_report_folders
			.collect()
			.subscribe {
				String folders = it.toString().replace("[", "").replace("]","") 
				log.info "\nFound the folder(s) \"$folders\" containing sequencing read files matching \"${params.extension}\" in \"${params.input}\".\n" }

		//Stop if folder count is 1
		ch_count_folders
			.count()
			.subscribe { if ( it == 1 ) exit 1, "Found only one folder with read data but \"--multipleSequencingRuns\" was specified. Please review data input." }
		
		//Stop if folder names contain "_" or "${params.split}"
		ch_check_folders
			.subscribe { 
				if ( it.toString().indexOf("${params.split}") > -1 ) exit 1, "Folder name \"$it\" contains \"${params.split}\", but may not. Please review data input or choose another string using \"--split [str]\" (no underscore allowed!)."
				if ( it.toString().indexOf("_") > -1 ) exit 1, "Folder name \"$it\" contains \"_\", but may not. Please review data input." 
			}

		//Add folder information to sequence files
		ch_rename_key
			.map { key, files -> [ key, files, (files[0].take(files[0].findLastIndexOf{"/"})[-1]) ] }
			.into { ch_read_pairs; ch_read_pairs_fastqc }

	} else if ( params.readPaths && params.multipleSequencingRuns ) {
		//Test input for multiple sequencing runs, profile = test_multi

		Channel
			.from(params.readPaths)
			.map { row -> [ row[0], [file(row[1][0]), file(row[1][1])], row[2] ] }
			.ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
			.into { ch_read_pairs; ch_read_pairs_fastqc }
			
	} else {
		//Standard input

		Channel
			.fromFilePairs( params.input + params.extension, size: 2 )
			.ifEmpty { exit 1, "Cannot find any reads matching: ${params.input}${params.extension}\nNB: Path needs to be enclosed in quotes!" }
			.map { name, reads -> [ name.toString().indexOf("_") != -1 ? name.toString().take(name.toString().indexOf("_")) : name, reads ] }
			.into { ch_read_pairs; ch_read_pairs_fastqc }
	}

	/*
	 * fastQC
	 */
	if (!params.multipleSequencingRuns){
		process fastqc {
			tag "${pair_id}"
			publishDir "${params.outdir}/fastQC", mode: params.publish_dir_mode,
			saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

			input:
			set val(pair_id), file(reads) from ch_read_pairs_fastqc

			output:
			file "*_fastqc.{zip,html}" into ch_fastqc_results

			when:
			!params.skip_fastqc

			script: 
			"""
			fastqc -q ${reads}
			"""
		}
	} else {
		process fastqc_multi {
			tag "${folder}${params.split}${pair_id}"
			publishDir "${params.outdir}/fastQC", mode: params.publish_dir_mode,
			saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

			input:
			set val(pair_id), file(reads), val(folder) from ch_read_pairs_fastqc

			output:
			file "*_fastqc.{zip,html}" into ch_fastqc_results

			when:
			!params.skip_fastqc

			script: 
			"""
			#Rename files so that there is no possible overlap
			ln -s "${reads[0]}" "$folder${params.split}${reads[0]}"
			ln -s "${reads[1]}" "$folder${params.split}${reads[1]}"
			fastqc -q "$folder${params.split}${reads[0]}" "$folder${params.split}${reads[1]}"
			"""
		}
	}

	/*
	 * Trim each read or read-pair with cutadapt
	 */
	if (!params.multipleSequencingRuns){
		process trimming {
			tag "${pair_id}"  
			publishDir "${params.outdir}/trimmed", mode: params.publish_dir_mode,
				saveAs: {filename -> 
				if (filename.indexOf(".gz") == -1) "logs/$filename"
				else if(params.keepIntermediates) filename 
				else null}
		
			input:
			set val(pair_id), file(reads) from ch_read_pairs
		
			output:
			set val(pair_id), file ("trimmed/*.*") into ch_fastq_trimmed_manifest 
			file "trimmed/*.*" into (ch_fastq_trimmed, ch_fastq_trimmed_qiime)
			file "cutadapt_log_*.txt" into ch_fastq_cutadapt_log

			script:
			discard_untrimmed = params.retain_untrimmed ? '' : '--discard-untrimmed'
			primers = single_end ? "--rc -g ${params.FW_primer}...${params.RV_primer}" : "-g ${params.FW_primer} -G ${params.RV_primer}"
			in_2_files = single_end ? "second-trimming_${reads}" : "second-trimming_${reads[0]} second-trimming_${reads[1]}" 
			out_1_files = single_end ? "-o second-trimming_${reads}" : "-o second-trimming_${reads[0]} -p second-trimming_${reads[1]}"
			in_1_files = single_end ? "first-trimming_${reads}" : "first-trimming_${reads[0]} first-trimming_${reads[1]}" //these have to be symlinked below
			out_files = single_end ? "-o trimmed/${reads}" : "-o trimmed/${reads[0]} -p trimmed/${reads[1]}"
			in_files = single_end ? "${reads}" : "${reads[0]} ${reads[1]}"
			"""
			mkdir -p trimmed
			if [[ \"${params.double_primer}\" = \"true\" && \"${params.retain_untrimmed}\" = \"false\" ]]; then
				#rename files to list results correctly in MultiQC
				if [ \"${single_end}\" = \"true\" ]; then
					ln -s "${reads}" "first-trimming_${reads}"
				else
					ln -s "${reads[0]}" "first-trimming_${reads[0]}"
					ln -s "${reads[1]}" "first-trimming_${reads[1]}"
				fi

				cutadapt ${primers} ${discard_untrimmed} \
					${out_1_files} \
					${in_1_files} >> cutadapt_log_${pair_id}.txt
				cutadapt ${primers} --discard-trimmed \
					${out_files} \
					${in_2_files} >> cutadapt_log_${pair_id}.txt
			else
				cutadapt ${primers} ${discard_untrimmed} \
					${out_files} ${in_files} \
					>> cutadapt_log_${pair_id}.txt
			fi
			"""
		}

	} else {
		process trimming_multi {
			tag "$folder${params.split}$pair_id"  
			publishDir "${params.outdir}/trimmed", mode: params.publish_dir_mode,
				saveAs: {filename -> 
				if (filename.indexOf(".gz") == -1) "logs/$filename"
				else if(params.keepIntermediates) filename 
				else null}
		
			input:
			set val(pair_id), file(reads), val(folder) from ch_read_pairs
		
			output:
			file "trimmed/*.*" into (ch_fastq_trimmed, ch_fastq_trimmed_manifest, ch_fastq_trimmed_qiime)
			file "cutadapt_log_$folder${params.split}${pair_id}.txt" into ch_fastq_cutadapt_log

			script:
			discard_untrimmed = params.retain_untrimmed ? '' : '--discard-untrimmed'
			"""
			mkdir -p trimmed
			if [[ \"${params.double_primer}\" = \"true\" && \"${params.retain_untrimmed}\" = \"false\" ]]; then
				mkdir -p firstcutadapt
				#Rename files so that MultiQC will pick them up correctely as first trimming
				ln -s "${reads[0]}" "first-trimming_$folder${params.split}${reads[0]}"
				ln -s "${reads[1]}" "first-trimming_$folder${params.split}${reads[1]}"
				cutadapt -g ${params.FW_primer} -G ${params.RV_primer} ${discard_untrimmed} \
					-o firstcutadapt/$folder${params.split}${reads[0]} -p firstcutadapt/$folder${params.split}${reads[1]} \
					"first-trimming_$folder${params.split}${reads[0]}" "first-trimming_$folder${params.split}${reads[1]}" >> cutadapt_log_$folder${params.split}${pair_id}.txt
				#Rename files so that MultiQC will pick them up correctely as second trimming
				ln -s "firstcutadapt/$folder${params.split}${reads[0]}" "second-trimming_$folder${params.split}${reads[0]}"
				ln -s "firstcutadapt/$folder${params.split}${reads[1]}" "second-trimming_$folder${params.split}${reads[1]}"
				cutadapt -g ${params.FW_primer} -G ${params.RV_primer} --discard-trimmed \
                                        -o trimmed/$folder${params.split}${reads[0]} -p trimmed/$folder${params.split}${reads[1]} \
                                        second-trimming_$folder${params.split}${reads[0]} second-trimming_$folder${params.split}${reads[1]} >> cutadapt_log_$folder${params.split}${pair_id}.txt
			else
				#first, rename files so that MultiQC will pick them up correctely
				ln -s "${reads[0]}" "$folder${params.split}${reads[0]}"
				ln -s "${reads[1]}" "$folder${params.split}${reads[1]}"
				cutadapt -g ${params.FW_primer} -G ${params.RV_primer} ${discard_untrimmed} \
                               		-o trimmed/$folder${params.split}${reads[0]} -p trimmed/$folder${params.split}${reads[1]} \
                                	$folder${params.split}${reads[0]} $folder${params.split}${reads[1]} > cutadapt_log_$folder${params.split}${pair_id}.txt
			fi
			"""
		}
	}
		
	
	/*
	 * multiQC
	 */
	process multiqc {
		publishDir "${params.outdir}/MultiQC", mode: params.publish_dir_mode

		input:
		file (multiqc_config) from ch_multiqc_config
		file (mqc_custom_config) from ch_multiqc_custom_config.collect().ifEmpty([])
		file ('cutadapt/logs/*') from ch_fastq_cutadapt_log.collect().ifEmpty([])
		file ('fastqc/*') from ch_fastqc_results.collect().ifEmpty([])
		file ('software_versions/*') from ch_software_versions_yaml.collect().ifEmpty([])
		file workflow_summary from ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")

		output:
		file "*multiqc_report.html" into ch_multiqc_report
		file "*_data"

		when:
		!params.skip_multiqc

		script:
		rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
		rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
		custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
		"""
		multiqc --interactive -f $rtitle $rfilename $custom_config_file .
		"""
	}


	/*
	* Produce manifest file for QIIME2
	*/
	if (!params.multipleSequencingRuns && single_end){
		ch_fastq_trimmed_manifest
			.map { name, reads ->
				def sampleID = name 
				def Reads = reads
				[ "${sampleID}" +","+ "${Reads}" + ",forward" ]
			}
			.flatten()
			.collectFile(name: 'manifest.txt', newLine: true, storeDir: "${params.outdir}/demux", seed: "sample-id,absolute-filepath,direction")
			.set { ch_manifest }
	
	} else if (!params.multipleSequencingRuns){
		ch_fastq_trimmed_manifest
			.map { name, reads ->
				def sampleID = name 
				def fwdReads = reads [0]
				def revReads = reads [1]
				[ "${sampleID}" +","+ "${fwdReads}" + ",forward\n" + "${sampleID}" +","+ "${revReads}" +",reverse" ]
			}
			.flatten()
			.collectFile(name: 'manifest.txt', newLine: true, storeDir: "${params.outdir}/demux", seed: "sample-id,absolute-filepath,direction")
			.set { ch_manifest }
	
	} else {
		ch_fastq_trimmed_manifest
			.map { forward, reverse -> [ forward.drop(forward.findLastIndexOf{"/"})[0], forward, reverse ] } //extract file name
			.map { name, forward, reverse -> [ name.toString().indexOf("_") != -1 ? name.toString().take(name.toString().indexOf("_")) : name, forward, reverse ] } //extract sample name
			.map { name, forward, reverse -> [ name +","+ forward + ",forward\n" + name +","+ reverse +",reverse" ] } //prepare basic synthax
			.flatten()
			.collectFile(storeDir: "${params.outdir}", seed: "sample-id,absolute-filepath,direction\n") { item ->
				def folder = item.take(item.indexOf("${params.split}")) //re-extract folder
				[ "${folder}${params.split}manifest.txt", item + '\n' ]
			}
			.set { ch_manifest_file }
			
		ch_manifest_file
			.combine( ch_mpl_for_qiime_import )
			.set { ch_manifest }
	}

	/*
	* Import trimmed files into QIIME2 artefact
	*/
	if (!params.multipleSequencingRuns && !params.pacbio) {
		process qiime_import {
			publishDir "${params.outdir}/demux", mode: params.publish_dir_mode, 
			saveAs: { filename -> 
				params.keepIntermediates ? filename : null
				params.untilQ2import ? filename : null }

			input:
			file(manifest) from ch_manifest
			env MATPLOTLIBRC from ch_mpl_for_qiime_import
			file('*') from ch_fastq_trimmed_qiime.collect()

			output:
			file "demux.qza" into (ch_qiime_demux_import, ch_qiime_demux_vis, ch_qiime_demux_dada)

			when:
			!params.Q2imported
		
			script:
			input_format = params.phred64 ? "PairedEndFastqManifestPhred64" : "PairedEndFastqManifestPhred33"
			"""
			head -n 1 ${manifest} > header.txt
			tail -n+2 ${manifest} | cut -d, -f1 > col1.txt
			tail -n+2 ${manifest} | cut -d, -f2 | sed 's:.*/::' > col2.txt
			while read f; do
				realpath \$f >> full_path.txt
			done <col2.txt
			tail -n+2 ${manifest} | cut -d, -f3 > col3.txt
			paste -d, col1.txt full_path.txt col3.txt > cols.txt
			cat cols.txt >> header.txt && mv header.txt ${manifest}

			qiime tools import \
				--type 'SampleData[PairedEndSequencesWithQuality]' \
				--input-path ${manifest} \
				--output-path demux.qza \
				--input-format $input_format
			"""
		}
	} else if (params.pacbio) {
		process qiime_import_pacbio{
			publishDir "${params.outdir}/demux", mode: 'copy', 
			saveAs: { filename -> 
				params.keepIntermediates ? filename : null
				params.untilQ2import ? filename : null }

			input:
			file(manifest) from ch_manifest
			env MATPLOTLIBRC from ch_mpl_for_qiime_import

			output:
			file manifest into ch_dada_import
			file "demux.qza" into (ch_qiime_demux_import, ch_qiime_demux_vis)

			when:
			!params.Q2imported
		
			script:
			input_format = params.phred64 ? "SingleEndFastqManifestPhred64" : "SingleEndFastqManifestPhred33"
			"""
			qiime tools import \
				--type 'SampleData[SequencesWithQuality]' \
				--input-path ${manifest} \
				--output-path demux.qza \
				--input-format $input_format
			"""
		}
	} else {
		process qiime_import_multi {
			tag "${manifest}"

			publishDir "${params.outdir}", mode: params.publish_dir_mode, 
			saveAs: { filename ->
				params.keepIntermediates ? filename : null}

			input:
			set file(manifest), env(MATPLOTLIBRC) from ch_manifest
			file('*') from ch_fastq_trimmed_qiime.collect()

			output:
			file "*demux.qza" into (ch_qiime_demux_import, ch_qiime_demux_vis, ch_qiime_demux_dada) mode flatten

			when:
			!params.Q2imported

			script:
			input_format = params.phred64 ? "PairedEndFastqManifestPhred64" : "PairedEndFastqManifestPhred33"
			def folder = "${manifest}".take("${manifest}".indexOf("${params.split}"))
			"""
			head -n 1 ${manifest} > header.txt
			tail -n+2 ${manifest} | cut -d, -f1 > col1.txt
			tail -n+2 ${manifest} | cut -d, -f2 | sed 's:.*/::' > col2.txt
			while read f; do
				realpath \$f >> full_path.txt
			done <col2.txt
			tail -n+2 ${manifest} | cut -d, -f3 > col3.txt
			paste -d, col1.txt full_path.txt col3.txt > cols.txt
			cat cols.txt >> header.txt && mv header.txt ${manifest}

			qiime tools import \
				--type 'SampleData[PairedEndSequencesWithQuality]' \
				--input-path ${manifest} \
				--output-path ${folder}-demux.qza \
				--input-format $input_format
			"""
		}
	}
	ch_qiime_demux_vis  
		.combine( ch_mpl_for_demux_visualize )
		.set{ ch_qiime_demux_visualisation }
} 

/*
 * Download, unpack, extract and train classifier
 * Download, unpack, and extract classifier in one process, train classifier in following process
 * Use "--dereplication 90" for testing and "--dereplication 99" for real datasets
 * Requirements with "--dereplication 99": 1 core (seems not to scale with more?), ~35 Gb mem, ~2:15:00 walltime
 */

if( !params.classifier ){
	if( !params.onlyDenoising && !params.skip_taxonomy ){
		Channel.fromPath("${params.reference_database}")
			.set { ch_ref_database }
	} else {
		Channel.empty()
			.set { ch_ref_database }		
	}

	process classifier_extract_seq {

		input:
		file database from ch_ref_database
		env MATPLOTLIBRC from ch_mpl_for_classifier_extract_seq

		output:
		file("*.qza") into ch_qiime_pretrain
		stdout ch_message_classifier_removeHash

		when:
		!params.onlyDenoising || !params.untilQ2import 

		script:
	  
		"""
		export HOME="\${PWD}/HOME"

		if [ ${params.taxon_reference} = \"unite\" ]; then
                        create_unite_taxfile.py $database db.fa db.tax
			fasta=\"db.fa\"
		        taxonomy=\"db.tax\"
		else
			unzip -qq $database
			fasta=\"SILVA_132_QIIME_release/rep_set/rep_set_16S_only/${params.dereplication}/silva_132_${params.dereplication}_16S.fna\"
			taxonomy=\"SILVA_132_QIIME_release/taxonomy/16S_only/${params.dereplication}/consensus_taxonomy_7_levels.txt\"
		fi


		if [ \"${params.classifier_removeHash}\" = \"true\" ]; then
			sed \'s/#//g\' \$taxonomy >taxonomy-${params.dereplication}_removeHash.txt
			taxonomy=\"taxonomy-${params.dereplication}_removeHash.txt\"
			echo \"\n######## WARNING! The taxonomy file was altered by removing all hash signs!\"
		fi

		### Import
		qiime tools import --type \'FeatureData[Sequence]\' \
			--input-path \$fasta \
			--output-path ref-seq-${params.dereplication}.qza
		qiime tools import --type \'FeatureData[Taxonomy]\' \
			--input-format HeaderlessTSVTaxonomyFormat \
			--input-path \$taxonomy \
			--output-path ref-taxonomy-${params.dereplication}.qza

		#Extract sequences based on primers
		qiime feature-classifier extract-reads \
			--i-sequences ref-seq-${params.dereplication}.qza \
			--p-f-primer ${params.FW_primer} \
			--p-r-primer ${params.RV_primer} \
			--o-reads ${params.FW_primer}-${params.RV_primer}-${params.dereplication}-ref-seq.qza \
			--quiet
		"""
	}

        ch_message_classifier_removeHash
                .subscribe { log.info it }

        process classifier_train {
                publishDir "${params.outdir}/DB/", mode: params.publish_dir_mode,
                saveAs: {filename ->
                        if (filename.indexOf("${params.FW_primer}-${params.RV_primer}-${params.dereplication}-classifier.qza") == 0) filename
                        else if(params.keepIntermediates) filename
                        else null}

                input:
                file '*' from ch_qiime_pretrain
                env MATPLOTLIBRC from ch_mpl_for_classifier_train

                output:
                file("${params.FW_primer}-${params.RV_primer}-${params.dereplication}-classifier.qza") into ch_qiime_classifier

                when:
                !params.onlyDenoising || !params.untilQ2import

                script:

                """
                export HOME="\${PWD}/HOME"

		#Train classifier
		qiime feature-classifier fit-classifier-naive-bayes \
			--i-reference-reads ${params.FW_primer}-${params.RV_primer}-${params.dereplication}-ref-seq.qza \
			--i-reference-taxonomy ref-taxonomy-${params.dereplication}.qza \
			--o-classifier ${params.FW_primer}-${params.RV_primer}-${params.dereplication}-classifier.qza \
			--quiet
		"""
	}

}

/*
 * Import trimmed files into QIIME2 artefact
 */
if( !params.Q2imported ){
	process qiime_demux_visualize {
		tag "${demux.baseName}"
		publishDir "${params.outdir}", mode: params.publish_dir_mode

		input:
		set file(demux), env(MATPLOTLIBRC) from ch_qiime_demux_visualisation

		output:
		file("${demux.baseName}/*-seven-number-summaries.csv") into ch_csv_demux
		file("${demux.baseName}/*")
	  
		"""
		export HOME="\${PWD}/HOME"
		
		qiime demux summarize \
		--i-data ${demux} \
		--o-visualization ${demux.baseName}.qzv

		qiime tools export --input-path ${demux.baseName}.qzv --output-path ${demux.baseName}
		"""
	}
} else {
	process qiime_importdemux_visualize { 
		publishDir "${params.outdir}", mode: params.publish_dir_mode

		input:
		env MATPLOTLIBRC from ch_mpl_for_demux_visualize

		output:
		file("demux/*-seven-number-summaries.csv") into ch_csv_demux
		file("demux/*")
	  
		"""
		export HOME="\${PWD}/HOME"

		qiime demux summarize \
			--i-data ${params.Q2imported} \
			--o-visualization demux.qzv

		qiime tools export --input-path demux.qzv --output-path demux
		"""
	}
}


/*
 * Determine params.trunclenf and params.trunclenr where the median quality value drops below params.trunc_qmin
 * But at least the fraction of params.trunc_rmin reads is retained
 * "Warning massage" is printed
 */
if ( ! single_end ) {
   process dada_trunc_parameter { 

	input:
	file summary_demux from ch_csv_demux 

	output:
	stdout ch_dada_trunc

	when:
	!params.untilQ2import

	script:
	if( params.trunclenf == false || params.trunclenr == false ){
		"""
		dada_trunc_parameter.py ${summary_demux[0]} ${summary_demux[1]} ${params.trunc_qmin} ${params.trunc_rmin}
		"""
	}
	else
		"""
		printf "${params.trunclenf},${params.trunclenr}"
		"""
  }
} else {
  process dada_trunc_se {

  	  output:
	  stdout ch_dada_trunc

	  when:
	  !params.untilQ2import

	  script:
	  if ( params.trunclenf == false ) {
	     """
	     printf "0"
	     """
	  } else {
	     """
	     printf "${params.trunclenf}"
	     """
	  }
  }
}

if (params.multipleSequencingRuns){
	//find minimum dada truncation values
	ch_dada_trunc
		.into { dada_trunc_forward; dada_trunc_reverse }
	dada_trunc_forward
		.map { trunc -> (trunc.split(',')[0]) }
		.min()
		.set { dada_forward }
	dada_trunc_reverse
		.map { trunc -> (trunc.split(',')[1]) }
		.min()
		.set { dada_reverse }
	dada_forward
		.combine( dada_reverse )
		.set { dada_trunc_multi }
	//combine channels for dada_multi
	ch_qiime_demux_dada
		.combine( dada_trunc_multi )
		.combine( ch_mpl_dada )
		.set { ch_dada_multi }
}

 
/*
 * Find ASVs with DADA2
 *    (i) for single sequencing run
 *    (ii) for PacBio reads
 *    (iii) for multiple sequencing runs
 */
if (!params.multipleSequencingRuns && !params.pacbio){
	process dada_single {
		tag "$trunc"
		publishDir "${params.outdir}", mode: params.publish_dir_mode,
			saveAs: {filename -> 
					 if (filename.indexOf("dada_stats/stats.tsv") == 0)         "abundance_table/unfiltered/dada_stats.tsv"
				else if (filename.indexOf("dada_report.txt") == 0)              "abundance_table/unfiltered/dada_report.txt"
				else if (filename.indexOf("table.qza") == 0)                    "abundance_table/unfiltered/$filename"
				else if (filename.indexOf("rel-table/feature-table.biom") == 0) "abundance_table/unfiltered/rel-feature-table.biom"
				else if (filename.indexOf("table/feature-table.biom") == 0)     "abundance_table/unfiltered/feature-table.biom"
				else if (filename.indexOf("rel-feature-table.tsv") > 0)         "abundance_table/unfiltered/rel-feature-table.tsv"
				else if (filename.indexOf("feature-table.tsv") > 0)             "abundance_table/unfiltered/feature-table.tsv"
				else if (filename.indexOf("rep-seqs.qza") == 0)                 "representative_sequences/unfiltered/rep-seqs.qza"
				else if (filename.indexOf("unfiltered/*"))                      "representative_sequences/$filename"
				else null}

		input:
		file demux from ch_qiime_demux_dada
		val trunc from ch_dada_trunc
		env MATPLOTLIBRC from ch_mpl_dada

		output:
		file("table.qza") into ch_qiime_table_raw
		file("rep-seqs.qza") into (ch_qiime_repseq_raw_for_classifier,ch_qiime_repseq_raw_for_filter)
		file("table/feature-table.tsv") into ch_tsv_table_raw
		file("dada_stats/stats.tsv")
		file("table/feature-table.biom")
		file("rel-table/feature-table.biom")
		file("table/rel-feature-table.tsv")
		file("unfiltered/*")
		file("dada_report.txt")

		when:
		!params.untilQ2import

		script:
		def values = trunc.split(',')
		if (values[0].toInteger() + values[1].toInteger() <= 10) { 
			log.info "\n######## ERROR: Total read pair length is below 10, this is definitely too low.\nForward ${values[0]} and reverse ${values[1]} are chosen.\nPlease provide appropriate values for --trunclenf and --trunclenr or lower --trunc_qmin\n" }
		"""
		export HOME="\${PWD}/HOME"
		IFS=',' read -r -a trunclen <<< \"$trunc\"

		#denoise samples with DADA2 and produce
		qiime dada2 denoise-paired  \
			--i-demultiplexed-seqs ${demux}  \
			--p-trunc-len-f \${trunclen[0]} \
			--p-trunc-len-r \${trunclen[1]} \
			--p-max-ee-f ${params.maxEE} \
			--p-max-ee-r ${params.maxEE} \
			--p-n-threads 0  \
			--o-table table.qza  \
			--o-representative-sequences rep-seqs.qza  \
			--o-denoising-stats stats.qza \
			--verbose \
		>dada_report.txt

		#produce dada2 stats "dada_stats/stats.tsv"
		qiime tools export --input-path stats.qza \
			--output-path dada_stats

		#produce raw count table in biom format "table/feature-table.biom"
		qiime tools export --input-path table.qza  \
			--output-path table

		#produce raw count table
		biom convert -i table/feature-table.biom \
			-o table/feature-table.tsv  \
			--to-tsv

		#produce representative sequence fasta file
		qiime feature-table tabulate-seqs  \
			--i-data rep-seqs.qza  \
			--o-visualization rep-seqs.qzv
		qiime tools export --input-path rep-seqs.qzv  \
			--output-path unfiltered

		#convert to relative abundances
		qiime feature-table relative-frequency \
			--i-table table.qza \
			--o-relative-frequency-table relative-table-ASV.qza

		#export to biom
		qiime tools export --input-path relative-table-ASV.qza \
			--output-path rel-table

		#convert to tab separated text file
		biom convert \
			-i rel-table/feature-table.biom \
			-o table/rel-feature-table.tsv --to-tsv
		"""
	}
} else if (params.pacbio){
	process dada_pacBio {
		
		publishDir "${params.outdir}", mode: 'copy',
			saveAs: {filename -> 
					 if (filename.indexOf("dada_stats.tsv") == 0)         "abundance_table/unfiltered/dada_stats.tsv"
				else if (filename.indexOf("dada_report.txt") == 0)              "abundance_table/unfiltered/dada_report.txt"
				else if (filename.indexOf("rel-feature-table.tsv") == 0)         "abundance_table/unfiltered/rel-feature-table.tsv"
				else if (filename.indexOf("feature-table.tsv") == 0)             "abundance_table/unfiltered/feature-table.tsv"
				else if (filename.indexOf("feature-table.biom") == 0)     "abundance_table/unfiltered/feature-table.biom"
				else if (filename.indexOf("sequences.fasta") == 0)                      "representative_sequences/unfiltered/sequences.fasta"
				else if (filename.indexOf("rep-seqs.qza") == 0)                 "representative_sequences/unfiltered/rep-seqs.qza"
				else null}

		input:
		file demux from ch_dada_import
		val trunc from ch_dada_trunc

		output:
		file("table.qza") into ch_qiime_table_raw
		file("rep-seqs.qza") into (ch_qiime_repseq_raw_for_classifier,ch_qiime_repseq_raw_for_filter)
		file("feature-table.tsv") into ch_tsv_table_raw
		file("dada_stats.tsv")
		file("feature-table.biom")
		file("rel-feature-table.tsv")
		file("sequences.fasta")
		file("dada_report.txt")

		when:
		!params.untilQ2import

		script:
		"""
		# Quality filtering with DADA2 filterAndTrim
		dada2_filter_pacbio.r --infile ${demux} --filterDir dada2_filtered --maxEE ${params.maxEE} --truncLen ${trunc} --minLen ${params.minLen} --maxLen ${params.maxLen} --stats filter_stats.tsv --verbose

		# Estimation of error models with DADA2 learnErrors
		dada2_errmodels_pacbio.r --filterDir dada2_filtered > err.out

		# Denoise samples with DADA2
		dada2_denoise_pacbio.r --filterDir dada2_filtered --errModel err.rds --pool TRUE  --verbose > dd.out

		# Chimera removal with DADA2, and produce
		# * raw count table "feature-table.tsv"
		# * relative abundancies "rel-feature-table.tsv"
                # * DADA2 stats to file "denoise_stats.tsv"
		# * representative sequences "sequences.fasta"
		dada2_chimrem.r --manifest ${demux} --dadaObj dd.rds --method "pooled" --allowOneOff TRUE --table feature-table.tsv --reltable rel-feature-table.tsv --repseqs sequences.fasta --stats denoise_stats.tsv

		# Create qiime2 object from representative sequences
		qiime tools import --type \'FeatureData[Sequence]\' \
			--input-path sequences.fasta \
			--output-path rep-seqs.qza

		# Create qiime2 object for feature table
		biom convert -i feature-table.tsv -o feature-table.biom --to-hdf5
#		make_biom_from_tsv feature-table.tsv feature-table.biom
		qiime tools import \
		        --input-path feature-table.biom \
			--type 'FeatureTable[Frequency]' \
			--input-format BIOMV210Format \
			--output-path table.qza


		# Produce dada2 report "dada_report.txt" from err.out and dd.out
		make_dada2_report_pacbio.py err.out dd.out pooled

		# Produce dada2 stats "dada_stats.tsv" from filter_stats.tsv and denoise_stats.tsv
		make_dada2_stats_pacbio.py filter_stats.tsv denoise_stats.tsv

		"""
	}
} else {
	process dada_multi {
		tag "${demux.baseName} ${trunclenf} ${trunclenr}"

		input:
		set file(demux), val(trunclenf), val(trunclenr), env(MATPLOTLIBRC) from ch_dada_multi

		output:
		file("${demux.baseName}-table.qza") into ch_qiime_table
		file("${demux.baseName}-rep-seqs.qza") into ch_qiime_repseq
		file("${demux.baseName}-stats.tsv") into ch_dada_stats
		file("${demux.baseName}-report.txt") into ch_dada_reports

		when:
		!params.untilQ2import

		script:
		if (trunclenf.toInteger() + trunclenr.toInteger() <= 10) { 
			log.info "\n######## ERROR: Total read pair length is below 10, this is definitely too low.\nForward ${trunclenf} and reverse ${trunclenr} are chosen.\nPlease provide appropriate values for --trunclenf and --trunclenr or lower --trunc_qmin\n" }
		"""
		export HOME="\${PWD}/HOME"

		#denoise samples with DADA2 and produce
		qiime dada2 denoise-paired  \
			--i-demultiplexed-seqs ${demux}  \
			--p-trunc-len-f ${trunclenf} \
			--p-trunc-len-r ${trunclenr} \
			--p-max-ee-f ${params.maxEE} \
			--p-max-ee-r ${params.maxEE} \
			--p-n-threads 0  \
			--o-table ${demux.baseName}-table.qza  \
			--o-representative-sequences ${demux.baseName}-rep-seqs.qza  \
			--o-denoising-stats ${demux.baseName}-stats.qza \
			--verbose \
			>${demux.baseName}-report.txt

		#produce dada2 stats "${demux.baseName}-dada_stats/stats.tsv"
		qiime tools export --input-path ${demux.baseName}-stats.qza \
			--output-path ${demux.baseName}-dada_stats
		cp ${demux.baseName}-dada_stats/stats.tsv ${demux.baseName}-stats.tsv
		"""
	}

	process dada_merge {
		tag "${tables}"
		publishDir "${params.outdir}", mode: params.publish_dir_mode,
			saveAs: {filename -> 
					 if (filename.indexOf("stats.tsv") == 0)                    "abundance_table/unfiltered/dada_stats.tsv"
				else if (filename.indexOf("dada_report.txt") == 0)              "abundance_table/unfiltered/dada_report.txt"
				else if (filename.indexOf("table.qza") == 0)                    "abundance_table/unfiltered/$filename"
				else if (filename.indexOf("rel-table/feature-table.biom") == 0) "abundance_table/unfiltered/rel-feature-table.biom"
				else if (filename.indexOf("table/feature-table.biom") == 0)     "abundance_table/unfiltered/feature-table.biom"
				else if (filename.indexOf("rel-feature-table.tsv") > 0)         "abundance_table/unfiltered/rel-feature-table.tsv"
				else if (filename.indexOf("feature-table.tsv") > 0)             "abundance_table/unfiltered/feature-table.tsv"
				else if (filename.indexOf("rep-seqs.qza") == 0)                 "representative_sequences/unfiltered/rep-seqs.qza"
				else if (filename.indexOf("unfiltered/*"))                      "representative_sequences/$filename"
				else null}

		input:
		file tables from ch_qiime_table.collect()
		file repseqs from ch_qiime_repseq.collect()
		file stats from ch_dada_stats.collect()
		file reports from ch_dada_reports.collect()
		env MATPLOTLIBRC from ch_mpl_dada_merge

		output:
		file("table.qza") into ch_qiime_table_raw
		file("rep-seqs.qza") into (ch_qiime_repseq_raw_for_classifier,ch_qiime_repseq_raw_for_filter)
		file("table/feature-table.tsv") into ch_tsv_table_raw
		file("stats.tsv")
		file("table/feature-table.biom")
		file("rel-table/feature-table.biom")
		file("table/rel-feature-table.tsv")
		file("unfiltered/*")
		file("dada_report.txt")

		when:
		!params.untilQ2import

		script:
		def TABLES = ''
		def REPSEQ = ''
		def STAT = ''
		def REPORT = ''
		for (table in tables) { TABLES+= " --i-tables ${table}" }
		for (repseq in repseqs) { REPSEQ+= " --i-data ${repseq}" }
		for (stat in stats) { STAT+= " $stat" }
		for (report in reports) { REPORT+= " $report" }
		"""
		export HOME="\${PWD}/HOME"

		#concatenate tables
		#merge files
		qiime feature-table merge \
			${TABLES} \
			--o-merged-table table.qza \
			--quiet
		qiime feature-table merge-seqs \
			${REPSEQ} \
			--o-merged-data rep-seqs.qza \
			--quiet
		cat ${STAT} >stats.tsv
		cat ${REPORT} >dada_report.txt

		#produce raw count table in biom format "table/feature-table.biom"
		qiime tools export --input-path table.qza  \
			--output-path table

		#produce raw count table
		biom convert -i table/feature-table.biom \
			-o table/feature-table.tsv  \
			--to-tsv

		#produce representative sequence fasta file
		qiime feature-table tabulate-seqs  \
			--i-data rep-seqs.qza  \
			--o-visualization rep-seqs.qzv
		qiime tools export --input-path rep-seqs.qzv  \
			--output-path unfiltered

		#convert to relative abundances
		qiime feature-table relative-frequency \
			--i-table table.qza \
			--o-relative-frequency-table relative-table-ASV.qza

		#export to biom
		qiime tools export --input-path relative-table-ASV.qza \
			--output-path rel-table

		#convert to tab separated text file
		biom convert \
			-i rel-table/feature-table.biom \
			-o table/rel-feature-table.tsv --to-tsv
		"""
	}
}

/*
 * Assign taxonomy to ASV sequences
 * Requirements: many cores, ~35 Gb mem, walltime scales with no. of ASV and ${params.classifier} = trained_classifier size (~15 min to several hours)
 * USE NXF feature of file size introduced in 0.32.0 here!!!
 */

process classifier { 
	publishDir "${params.outdir}", mode: params.publish_dir_mode,
		saveAs: {filename -> 
			if (filename == "taxonomy/taxonomy.tsv") filename
			else if (filename == "taxonomy.qza") "taxonomy/$filename"}

	input:
	file repseq from ch_qiime_repseq_raw_for_classifier
	file trained_classifier from ch_qiime_classifier
	env MATPLOTLIBRC from ch_mpl_classifier

	output:
	file("taxonomy.qza") into (ch_qiime_taxonomy_for_filter,ch_qiime_taxonomy_for_relative_abundance_reduced_taxa,ch_qiime_taxonomy_for_barplot,ch_qiime_taxonomy_for_ancom,ch_qiime_taxonomy_for_export_filtered_dada_output)
	file("taxonomy/taxonomy.tsv") into ch_tsv_taxonomy

  
	"""
	export HOME="\${PWD}/HOME"

	qiime feature-classifier classify-sklearn  \
		--i-classifier ${trained_classifier}  \
		--p-n-jobs ${task.cpus}  \
		--i-reads ${repseq}  \
		--o-classification taxonomy.qza  \
		--verbose

	qiime metadata tabulate  \
		--m-input-file taxonomy.qza  \
		--o-visualization taxonomy.qzv  \
		--verbose

	#produce "taxonomy/taxonomy.tsv"
	qiime tools export --input-path taxonomy.qza  \
		--output-path taxonomy

	qiime tools export --input-path taxonomy.qzv  \
		--output-path taxonomy
	"""
}

/*
 * Filter out unwanted/off-target taxa
 */
if (params.exclude_taxa == "none" && !params.min_frequency && !params.min_samples) {

	ch_qiime_repseq_raw_for_filter
		.into{ ch_qiime_repseq_for_dada_output; ch_qiime_repseq_for_tree }

	ch_qiime_table_raw
		.into{ ch_qiime_table_for_filtered_dada_output; ch_qiime_table_for_relative_abundance_asv; ch_qiime_table_for_relative_abundance_reduced_taxa; ch_qiime_table_for_ancom; ch_qiime_table_for_barplot; ch_qiime_table_for_alpha_rarefaction; ch_qiime_table_for_diversity_core }

} else {
	process filter_taxa {
		tag "taxa:${params.exclude_taxa};min-freq:${params.min_frequency};min-samples:${params.min_samples}"

		publishDir "${params.outdir}", mode: params.publish_dir_mode,
		saveAs: {filename -> 
				 if (filename.indexOf("filtered-table.qza") == 0)      "abundance_table/filtered/table.qza"
			else if (filename.indexOf("filtered-sequences.qza") == 0)  "representative_sequences/filtered/rep-seqs.qza"
			else null}  

		input:
		file table from ch_qiime_table_raw
		file repseq from  ch_qiime_repseq_raw_for_filter
		file taxonomy from ch_qiime_taxonomy_for_filter
		env MATPLOTLIBRC from ch_mpl_filter_taxa

		output:
		file("filtered-table.qza") into (ch_qiime_table_for_filtered_dada_output, ch_qiime_table_for_relative_abundance_asv,ch_qiime_table_for_relative_abundance_reduced_taxa,ch_qiime_table_for_ancom,ch_qiime_table_for_barplot,ch_qiime_table_for_alpha_rarefaction, ch_qiime_table_for_diversity_core)
		file("filtered-sequences.qza") into (ch_qiime_repseq_for_dada_output,ch_qiime_repseq_for_tree)

		script:
		if ( "${params.min_frequency}" == "false" ) { minfrequency = 1 } else { minfrequency = "${params.min_frequency}" }
		if ( "${params.min_samples}" == "false" ) { minsamples = 1 } else { minsamples = "${params.min_samples}" }
		//if ( "${params.exclude_taxa}" == "none" ) { exclude = "" } else { exclude = "--p-exclude ${params.exclude_taxa} --p-mode contains " }
		"""
		export HOME="\${PWD}/HOME"

		if ! [ \"${params.exclude_taxa}\" = \"none\" ]; then
			#filter sequences
			qiime taxa filter-seqs \
				--i-sequences ${repseq} \
				--i-taxonomy ${taxonomy} \
				--p-exclude ${params.exclude_taxa} --p-mode contains \
				--o-filtered-sequences tax_filtered-sequences.qza

			#filter abundance table
			qiime taxa filter-table \
				--i-table ${table} \
				--i-taxonomy ${taxonomy} \
				--p-exclude ${params.exclude_taxa} --p-mode contains \
				--o-filtered-table tax_filtered-table.qza

			filtered_table="tax_filtered-table.qza"
			filtered_sequences="tax_filtered-sequences.qza"
		else
			filtered_table=${table}
			filtered_sequences=${repseq}
		fi

		qiime feature-table filter-features \
			--i-table \$filtered_table \
			--p-min-frequency ${minfrequency} \
			--p-min-samples ${minsamples} \
			--o-filtered-table filtered-table.qza
		
		qiime feature-table filter-seqs \
			--i-data \$filtered_sequences \
			--i-table filtered-table.qza \
			--o-filtered-data filtered-sequences.qza
		"""
	}
}

/*
 * Export qiime artefacts from filtered dada output
 */
process export_filtered_dada_output { 
	publishDir "${params.outdir}", mode: params.publish_dir_mode,
		saveAs: {filename -> 
				 if (filename.indexOf("table/feature-table.biom") == 0)  "abundance_table/filtered/feature-table.biom"
			else if (filename.indexOf("table/feature-table.tsv") == 0)   "abundance_table/filtered/feature-table.tsv"
			else if (filename.indexOf("abs-abund-table-") == 0)          "abundance_table/filtered/$filename"
			else if (filename.indexOf("filtered/*"))                     "representative_sequences/$filename"
			else null}   

	input:
	file table from ch_qiime_table_for_filtered_dada_output
	file repseq from ch_qiime_repseq_for_dada_output
	file taxonomy from ch_qiime_taxonomy_for_export_filtered_dada_output
	env MATPLOTLIBRC from ch_mpl_for_export_dada_output

	output:
	file("filtered/sequences.fasta") into ch_fasta_repseq
	file("table/feature-table.tsv") into (ch_tsv_table_for_alpha_rarefaction,ch_tsv_table_for_report_filter_stats,ch_tsv_table_for_diversity_core)
	file("table/feature-table.biom")
	file("filtered/*")
	file("abs-abund-table-*.tsv")

	"""
	export HOME="\${PWD}/HOME"

	#produce raw count table in biom format "table/feature-table.biom"
	qiime tools export --input-path ${table}  \
		--output-path table

	#produce raw count table "table/feature-table.tsv"
	biom convert -i table/feature-table.biom \
		-o table/feature-table.tsv  \
		--to-tsv

	#produce representative sequence fasta file "${params.outdir}/representative_sequences/sequences.fasta"
	qiime feature-table tabulate-seqs  \
		--i-data ${repseq}  \
		--o-visualization rep-seqs.qzv
	qiime tools export --input-path rep-seqs.qzv  \
		--output-path filtered

	##on several taxa level
	array=( 2 3 4 5 6 7 )
	for i in \${array[@]}
	do
		#collapse taxa
		qiime taxa collapse \
			--i-table ${table} \
			--i-taxonomy ${taxonomy} \
			--p-level \$i \
			--o-collapsed-table table-\$i.qza
		#export to biom
		qiime tools export --input-path table-\$i.qza \
			--output-path table-\$i
		#convert to tab separated text file
		biom convert \
			-i table-\$i/feature-table.biom \
			-o abs-abund-table-\$i.tsv --to-tsv
	done
	"""
}

/*
 * Report stats after taxa filtering
 */
process report_filter_stats { 
	publishDir "${params.outdir}/abundance_table/filtered", mode: params.publish_dir_mode     

	input:
	file 'unfiltered_table' from ch_tsv_table_raw
	file 'filtered_table' from ch_tsv_table_for_report_filter_stats

	output:
	file("count_table_filter_stats.tsv")
	
	"""
	count_table_filter_stats.py unfiltered_table filtered_table
	"""
}

/*
 * Export relative abundance tables on ASV level
 */
process RelativeAbundanceASV { 
	publishDir "${params.outdir}/rel_abundance_tables", mode: params.publish_dir_mode    

	input:
	file table from ch_qiime_table_for_relative_abundance_asv
	env MATPLOTLIBRC from ch_mpl_for_relasv

	output:
	file("rel-table-ASV.tsv") into ch_tsv_relASV_table

	when:
	!params.skip_abundance_tables

	"""
	export HOME="\${PWD}/HOME"

	#convert to relative abundances
	qiime feature-table relative-frequency \
		--i-table ${table} \
		--o-relative-frequency-table relative-table-ASV.qza

	#export to biom
	qiime tools export --input-path relative-table-ASV.qza --output-path relative-table-ASV

	#convert to tab separated text file "${params.outdir}/rel-table-ASV.tsv"
	biom convert -i relative-table-ASV/feature-table.biom \
		-o rel-table-ASV.tsv --to-tsv
	"""
}


/*
 * Export relative abundance tables based on taxonomic levels
 */
process RelativeAbundanceReducedTaxa { 
	publishDir "${params.outdir}/rel_abundance_tables", mode: params.publish_dir_mode    

	input:
	file table from ch_qiime_table_for_relative_abundance_reduced_taxa
	file taxonomy from ch_qiime_taxonomy_for_relative_abundance_reduced_taxa
	env MATPLOTLIBRC from ch_mpl_for_relreducetaxa

	output:
	file("*.tsv")

	when:
	!params.skip_abundance_tables && !params.skip_taxonomy

	"""
	export HOME="\${PWD}/HOME"
	##on several taxa level

	array=( 2 3 4 5 6 7 )
	for i in \${array[@]}
	do
		#collapse taxa
		qiime taxa collapse \
			--i-table ${table} \
			--i-taxonomy ${taxonomy} \
			--p-level \$i \
			--o-collapsed-table table-\$i.qza
		#convert to relative abundances
		qiime feature-table relative-frequency \
			--i-table table-\$i.qza \
			--o-relative-frequency-table relative-table-\$i.qza
		#export to biom
		qiime tools export --input-path relative-table-\$i.qza \
			--output-path relative-table-\$i
		#convert to tab separated text file
		biom convert \
			-i relative-table-\$i/feature-table.biom \
			-o rel-table-\$i.tsv --to-tsv
	done

	"""
}


/*
 * Produce a bar plot
 */
process barplot { 
	publishDir "${params.outdir}", mode: params.publish_dir_mode    

	input:
	file metadata from ch_metadata_for_barplot
	file table from ch_qiime_table_for_barplot
	file taxonomy from ch_qiime_taxonomy_for_barplot
	env MATPLOTLIBRC from ch_mpl_for_barcode

	output:
	file("barplot/*")

	when:
	!params.skip_barplot && !params.skip_taxonomy
  
	"""
	export HOME="\${PWD}/HOME"

	qiime taxa barplot  \
		--i-table ${table}  \
		--i-taxonomy ${taxonomy}  \
		--m-metadata-file ${metadata}  \
		--o-visualization taxa-bar-plots.qzv  \
		--verbose

	qiime tools export --input-path taxa-bar-plots.qzv  \
		--output-path barplot
	"""
}

/*
 * Produce a rooted tree
 */
process tree { 
	publishDir "${params.outdir}", mode: params.publish_dir_mode,
		saveAs: {filename -> 
			if (filename.indexOf("rooted-tree.qza") == 0)  "phylogenetic_tree/$filename"
			else filename }   

	input:
	file repseq from ch_qiime_repseq_for_tree
	env MATPLOTLIBRC from ch_mpl_for_tree

	output:
	file("rooted-tree.qza") into (ch_qiime_tree_for_diversity_core, ch_qiime_tree_for_alpha_rarefaction)
	file("phylogenetic_tree/tree.nwk")

	when:
	!params.skip_diversity_indices || !params.skip_alpha_rarefaction

  
	"""
	export HOME="\${PWD}/HOME"

	qiime alignment mafft \
		--i-sequences ${repseq} \
		--o-alignment aligned-rep-seqs.qza \
		--p-n-threads ${task.cpus}

	qiime alignment mask \
		--i-alignment aligned-rep-seqs.qza \
		--o-masked-alignment masked-aligned-rep-seqs.qza

	qiime phylogeny fasttree \
		--i-alignment masked-aligned-rep-seqs.qza \
		--p-n-threads ${task.cpus} \
		--o-tree unrooted-tree.qza

	qiime phylogeny midpoint-root \
		--i-tree unrooted-tree.qza \
		--o-rooted-tree rooted-tree.qza

	qiime tools export --input-path rooted-tree.qza  \
		--output-path phylogenetic_tree
	"""
}


/*
 * Alpha-rarefaction
 */
process alpha_rarefaction { 
	publishDir "${params.outdir}", mode: params.publish_dir_mode    

	input:
	file metadata from ch_metadata_for_alphararefaction
	file table from ch_qiime_table_for_alpha_rarefaction
	file tree from ch_qiime_tree_for_alpha_rarefaction
	file stats from ch_tsv_table_for_alpha_rarefaction
	env MATPLOTLIBRC from ch_mpl_for_alpha_rare

	output:
	file("alpha-rarefaction/*")

	when:
	!params.skip_alpha_rarefaction

	"""
	export HOME="\${PWD}/HOME"

	maxdepth=\$(count_table_minmax_reads.py $stats maximum 2>&1)

	#check values
	if [ \"\$maxdepth\" -gt \"75000\" ]; then maxdepth=\"75000\"; fi
	if [ \"\$maxdepth\" -gt \"5000\" ]; then maxsteps=\"250\"; else maxsteps=\$((maxdepth/20)); fi

	qiime diversity alpha-rarefaction  \
		--i-table ${table}  \
		--i-phylogeny ${tree}  \
		--p-max-depth \$maxdepth  \
		--m-metadata-file ${metadata}  \
		--p-steps \$maxsteps  \
		--p-iterations 10  \
		--o-visualization alpha-rarefaction.qzv

	qiime tools export --input-path alpha-rarefaction.qzv  \
		--output-path alpha-rarefaction
	"""
}

/*
 * Combine abundances, sequences and taxonomic classification into one table with R
 */
process combinetable { 
	publishDir "${params.outdir}/rel_abundance_tables", mode: params.publish_dir_mode    

	input:
	file TABLE from ch_tsv_relASV_table
	file SEQ from ch_fasta_repseq
	file TAXONOMY from ch_tsv_taxonomy

	output:
	file("qiime2_ASV_table.tsv")

	when:
	!params.skip_abundance_tables && !params.skip_taxonomy

	"""
	combineTable.r ${TABLE} ${SEQ} ${TAXONOMY}
	"""
}

/*
 * Compute diversity matrices
 */
process diversity_core { 
	publishDir "${params.outdir}", mode: params.publish_dir_mode,
	saveAs: {filename ->
		params.keepIntermediates ? filename : null} 

	input:
	file metadata from ch_metadata_for_diversity_core
	file table from ch_qiime_table_for_diversity_core
	file tree from ch_qiime_tree_for_diversity_core
	file stats from ch_tsv_table_for_diversity_core
	env MATPLOTLIBRC from ch_mpl_for_diversity_core

	output:
	file("diversity_core/*_pcoa_results.qza") into (ch_qiime_diversity_core_for_beta_diversity_ordination) mode flatten
	file("diversity_core/*_vector.qza") into ch_qiime_diversity_core_for_alpha_diversity mode flatten
	file("diversity_core/*_distance_matrix.qza") into ch_qiime_diversity_core_for_beta_diversity mode flatten
	stdout rarefaction_depth

	when:
	!params.skip_diversity_indices

	"""
	export HOME="\${PWD}/HOME"
	mindepth=\$(count_table_minmax_reads.py $stats minimum 2>&1)

	if [ \"\$mindepth\" -gt \"10000\" ]; then echo \"\nUse the sampling depth of \$mindepth for rarefaction\" ; fi
	if [ \"\$mindepth\" -lt \"10000\" -a \"\$mindepth\" -gt \"5000\" ]; then echo \"\n######## WARNING! The sampling depth of \$mindepth is quite small for rarefaction!\" ; fi
	if [ \"\$mindepth\" -lt \"5000\" -a \"\$mindepth\" -gt \"1000\" ]; then echo \"\n######## WARNING! The sampling depth of \$mindepth is very small for rarefaction!\" ; fi
	if [ \"\$mindepth\" -lt \"1000\" ]; then echo \"\n######## ERROR! The sampling depth of \$mindepth seems too small for rarefaction!\" ; fi
	
	qiime diversity core-metrics-phylogenetic \
		--m-metadata-file ${metadata} \
		--i-phylogeny ${tree} \
		--i-table ${table} \
		--p-sampling-depth \$mindepth \
		--output-dir diversity_core \
		--p-n-jobs ${task.cpus} \
		--quiet
	"""
}

rarefaction_depth
	.subscribe { log.info it }

/*
 * Capture all possible metadata categories for statistics
 */
process metadata_category_all { 
	input:
	file metadata from ch_metadata_for_metadata_category_all
	env MATPLOTLIBRC from ch_mpl_for_metadata_cat

	output:
	stdout into (ch_meta_category_all_for_alphadiversity, ch_meta_category_all_for_ancom)

	when:
	(!params.skip_ancom || !params.skip_diversity_indices) &&
	(!params.untilQ2import && !params.onlyDenoising)

	script:
	if( !params.metadata_category )
		"""
		metadataCategory.r ${metadata}
		"""
	else
		"""
		printf ${params.metadata_category}
		"""
}

/*
 * Capture all pairwise metadata categories for statistics
 */
process metadata_category_pairwise { 

	input:
	file metadata from ch_metadata_for_metadata_category_pairwise
	env MATPLOTLIBRC from ch_mpl_for_metadata_pair

	output:
	stdout ch_meta_category_pairwise

	when:
	!params.skip_diversity_indices

	"""
	metadataCategoryPairwise.r ${metadata}
	"""
}

/*
 * Combine channels for diversity analysis
 */

ch_metadata_for_alpha_diversity
	.combine( ch_qiime_diversity_core_for_alpha_diversity )
	.combine( ch_mpl_for_alpha_diversity )
	.combine( ch_meta_category_all_for_alphadiversity )
	.set{ ch_for_alpha_diversity }
ch_metadata_for_beta_diversity
	.combine( ch_qiime_diversity_core_for_beta_diversity )
	.combine( ch_meta_category_pairwise )
	.combine( ch_mpl_for_beta_diversity )
	.set{ ch_for_beta_diversity }
ch_metadata_for_beta_diversity_ordination
	.combine( ch_qiime_diversity_core_for_beta_diversity_ordination )
	.combine( ch_mpl_for_beta_diversity_ord )
	.set{ ch_for_beta_diversity_ordination }
	


/*
 * Compute alpha diversity indices
 */
process alpha_diversity { 
	tag "${core.baseName}"
	publishDir "${params.outdir}", mode: params.publish_dir_mode    

	input:
	set file(metadata), file(core), env(MATPLOTLIBRC), val(meta) from ch_for_alpha_diversity

	output:
	file("alpha-diversity/*")

	when:
	meta.length() > 0

	"""
	export HOME="\${PWD}/HOME"

	qiime diversity alpha-group-significance \
		--i-alpha-diversity ${core} \
		--m-metadata-file ${metadata} \
		--o-visualization ${core.baseName}-vis.qzv
	qiime tools export --input-path ${core.baseName}-vis.qzv \
		--output-path "alpha-diversity/${core.baseName}"
	"""
}


/*
 * Compute beta diversity indices
 */
process beta_diversity { 
	tag "${core.baseName}"
	publishDir "${params.outdir}", mode: params.publish_dir_mode     

	input:
	set file(meta), file(core), val(category), env(MATPLOTLIBRC) from ch_for_beta_diversity

	output:
	file("beta-diversity/*")

	when:
	category.length() > 0

	"""
	export HOME="\${PWD}/HOME"
	IFS=',' read -r -a metacategory <<< \"$category\"

	for j in \"\${metacategory[@]}\"
	do
		qiime diversity beta-group-significance \
			--i-distance-matrix $core \
			--m-metadata-file ${meta} \
			--m-metadata-column \"\$j\" \
			--o-visualization ${core.baseName}-\$j.qzv \
			--p-pairwise
		qiime tools export --input-path ${core.baseName}-\$j.qzv \
			--output-path beta-diversity/${core.baseName}-\$j
	done
	"""
}

/*
 * Compute beta diversity ordination
 */
process beta_diversity_ordination { 
	tag "${core.baseName}"
	publishDir "${params.outdir}", mode: params.publish_dir_mode

	input:
	set file(metadata), file(core), env(MATPLOTLIBRC) from ch_for_beta_diversity_ordination

	output:
	file("beta-diversity/*")

	"""
	export HOME="\${PWD}/HOME"

	qiime emperor plot \
		--i-pcoa ${core} \
		--m-metadata-file ${metadata} \
		--o-visualization ${core.baseName}-vis.qzv
	qiime tools export --input-path ${core.baseName}-vis.qzv \
		--output-path beta-diversity/${core.baseName}-PCoA
	"""
}


/*
 * Differential abundance analysis with ANCOM
 */
process prepare_ancom { 
	tag "${meta}"

	publishDir "${params.outdir}/ancom", mode: params.publish_dir_mode, 
	saveAs: {filename ->
		params.keepIntermediates ? filename : null}   

	input:
	file metadata from ch_metadata_for_ancom
	file table from ch_qiime_table_for_ancom
	val meta from ch_meta_category_all_for_ancom
	env MATPLOTLIBRC from ch_mpl_for_ancom

	output:
	file("*.qza") into (ch_meta_tables_tax, ch_meta_tables_asv) mode flatten

	when:
	!params.skip_ancom && (meta.length() > 0)

	"""
	export HOME="\${PWD}/HOME"
	IFS=',' read -r -a metacategory <<< \"$meta\"

	#remove samples that do not have any value
	for j in \"\${metacategory[@]}\"
	do
		qiime feature-table filter-samples \
			--i-table ${table} \
			--m-metadata-file ${metadata} \
			--p-where \"\$j<>\'\'\" \
			--o-filtered-table \$j.qza
	done
	"""
}

/*
 * Combine channels for ancom
 */

ch_taxlevel_tax = Channel.from( 2, 3, 4, 5, 6 )

ch_meta_tables_tax
	.combine( ch_taxlevel_tax )
	.combine( ch_qiime_taxonomy_for_ancom )
	.combine( ch_metadata_for_ancom_tax )
	.combine( ch_mpl_for_ancom_tax )
	.set{ ch_for_ancom_tax }

ch_meta_tables_asv
	.combine( ch_metadata_for_ancom_asv )
	.combine ( ch_mpl_for_ancom_asv )
	.set{ ch_for_ancom_asv }


/*
 * Differential abundance analysis with ANCOM on various taxonomic levels
 */
process ancom_tax { 
	tag "${table.baseName}-level${taxlevel}"

	publishDir "${params.outdir}", mode: params.publish_dir_mode    

	input:
	set file(table), val(taxlevel), file(taxonomy), file(metadata), env(MATPLOTLIBRC) from ch_for_ancom_tax

	output:
	file("ancom/*")

	when:
	!params.skip_ancom

	"""
	export HOME="\${PWD}/HOME"

	qiime taxa collapse \
		--i-table ${table} \
		--i-taxonomy ${taxonomy} \
		--p-level ${taxlevel} \
		--o-collapsed-table lvl${taxlevel}-${table}
	qiime composition add-pseudocount \
		--i-table lvl${taxlevel}-${table} \
		--o-composition-table comp-lvl${taxlevel}-${table}
	qiime composition ancom \
		--i-table comp-lvl${taxlevel}-${table} \
		--m-metadata-file ${metadata} \
		--m-metadata-column ${table.baseName} \
		--o-visualization comp-lvl${taxlevel}-${table.baseName}.qzv
	qiime tools export --input-path comp-lvl${taxlevel}-${table.baseName}.qzv \
		--output-path ancom/Category-${table.baseName}-level-${taxlevel}
	"""
}

/*
 * Differential abundance analysis with ANCOM on ASV level
 */
process ancom_asv { 
	tag "${table.baseName}"

	publishDir "${params.outdir}", mode: params.publish_dir_mode 

	input:
	set file(table), file(metadata), env(MATPLOTLIBRC) from ch_for_ancom_asv 

	output:
	file("ancom/*") 

	"""
	export HOME="\${PWD}/HOME"
	
	qiime composition add-pseudocount \
		--i-table ${table} \
		--o-composition-table comp-${table}
	qiime composition ancom \
		--i-table comp-${table} \
		--m-metadata-file ${metadata} \
		--m-metadata-column ${table.baseName} \
		--o-visualization comp-${table.baseName}.qzv
	qiime tools export --input-path comp-${table.baseName}.qzv \
		--output-path ancom/Category-${table.baseName}-ASV
	"""
}

/*
 * STEP 3 - Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    file output_docs from ch_output_docs
    file images from ch_output_docs_images

	output:
	file "results_description.html"

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/ampliseq] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[nf-core/ampliseq] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = ch_multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList) {
                log.warn "[nf-core/ampliseq] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/ampliseq] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$projectDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$projectDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, projectDir: "$projectDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$projectDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/ampliseq] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            if ( mqc_report.size() <= params.max_multiqc_email_size.toBytes() ) {
              mail_cmd += [ '-A', mqc_report ]
            }
            mail_cmd.execute() << email_html
            log.info "[nf-core/ampliseq] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/ampliseq]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/ampliseq]${c_red} Pipeline completed with errors${c_reset}-"
    }

}


def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/ampliseq v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}
