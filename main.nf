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
    log.info"""
    =======================================================
                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\'
        |\\ | |__  __ /  ` /  \\ |__) |__         }  {
        | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                              `._,._,\'

     nf-core/ampliseq : v${workflow.manifest.version}
    =======================================================
    
    Usage:

    The minimal command for running the pipeline is as follows:
    nextflow run nf-core/ampliseq --reads "data" --FW_primer GTGYCAGCMGCCGCGGTAA --RV_primer GGACTACNVGGGTWTCTAAT --metadata "Metadata.tsv"


    Required arguments:
      --reads [Path to folder]      Folder containing Casava 1.8 paired-end demultiplexed fastq files: *_L001_R{1,2}_001.fastq.gz
                                    Note: All samples have to be sequenced in one run, otherwise also specifiy --multipleSequencingRuns
      --FW_primer [str]             Forward primer sequence
      --RV_primer [str]             Reverse primer sequence
      --metadata                    Path to metadata sheet

    Filters:
      --exclude_taxa [str]          Comma seperated list of unwanted taxa (default: "mitochondria,chloroplast")
                                    To skip filtering use "none"
      --min_frequency [int]         Remove entries from the feature table below an absolute abundance threshold (default: 1)
      --min_samples [int]           Filtering low prevalent features from the feature table (default: 1)                   

    Cutoffs:
      --retain_untrimmed            Cutadapt will retain untrimmed reads
      --trunclenf [int]             DADA2 read truncation value for forward strand
      --trunclenr [int]             DADA2 read truncation value for reverse strand
      --trunc_qmin [int]            If --trunclenf and --trunclenr are not set, 
                                    these values will be automatically determined using 
                                    this mean quality score (not preferred) (default: 25)

    References:                     If you have trained a compatible classifier before
      --classifier                  Path to QIIME2 classifier file (typically *-classifier.qza)
      --classifier_removeHash       Remove all hash signs from taxonomy strings, resolves a rare ValueError during classification (process classifier)

    Statistics:
      --metadata_category           Diversity indices will be calculated using these groupings in the metadata sheet,
                                    all suitable columns in the metadata sheet will be used if not specified.
                                    Suitable are columns which are categorical (not numerical) and have multiple  
                                    different values which are not all unique.

    Other options:
      --untilQ2import               Skip all steps after importing into QIIME2, used for visually choosing DADA2 parameter
      --Q2imported [Path]           Path to imported reads (e.g. "demux.qza"), used after visually choosing DADA2 parameter
      --onlyDenoising               Skip all steps after denoising, produce only sequences and abundance tables on ASV level
      --multipleSequencingRuns      If samples were sequenced in multiple sequencing runs. Expects one subfolder per sequencing run
                                    in the folder specified by --reads containing sequencing data of the specific run. Also, fastQC
                                    is skipped because multiple sequencing runs might create overlapping file names that crash MultiQC.

    Skipping steps:
      --skip_fastqc                 Skip FastQC
      --skip_alpha_rarefaction      Skip alpha rarefaction
      --skip_taxonomy               Skip taxonomic classification
      --skip_barplot                Skip producing barplot
      --skip_abundance_tables       Skip producing any relative abundance tables
      --skip_diversity_indices      Skip alpha and beta diversity analysis
      --skip_ancom                  Skip differential abundance testing     

    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

// Configurable variables
params.name = false
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.email = false
params.plaintext_email = false

ch_multiqc_config = Channel.fromPath(params.multiqc_config)
ch_output_docs = Channel.fromPath("$baseDir/docs/output.md")
Channel.fromPath("$baseDir/assets/matplotlibrc")
    .into { ch_mpl_for_make_classifier; ch_mpl_for_qiime_import; ch_mpl_for_ancom_asv; ch_mpl_for_ancom_tax; ch_mpl_for_ancom; ch_mpl_for_beta_diversity_ord; ch_mpl_for_beta_diversity; ch_mpl_for_alpha_diversity; ch_mpl_for_metadata_pair; ch_mpl_for_metadata_cat; ch_mpl_for_diversity_core; ch_mpl_for_alpha_rare; ch_mpl_for_tree; ch_mpl_for_barcode; ch_mpl_for_relreducetaxa; ch_mpl_for_relasv; ch_mpl_for_export_dada_output; ch_mpl_filter_taxa; ch_mpl_classifier; ch_mpl_dada; ch_mpl_dada_merge; ch_mpl_for_demux_visualize; ch_mpl_for_classifier }

// Defines all parameters that are independent of a test run
params.trunc_qmin = 25 //to calculate params.trunclenf and params.trunclenr automatically
params.trunclenf = false
params.trunclenr = false
params.metadata_category = false
params.tree_cores = 2
params.diversity_cores = 2
params.retain_untrimmed = false
params.exclude_taxa = "mitochondria,chloroplast"
params.keepIntermediates = false
params.classifier_removeHash = false
params.min_frequency = false
params.min_samples = false
params.multipleSequencingRuns = false

//Database specific parameters
//currently only this is compatible with process make_SILVA_132_16S_classifier
params.reference_database = "https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_132_release.zip"
params.dereplication = 99


/*
 * Defines pipeline steps
 */

Channel.fromPath("${params.metadata}")
        .into { ch_metadata_for_barplot; ch_metadata_for_alphararefaction; ch_metadata_for_diversity_core; ch_metadata_for_alpha_diversity; ch_metadata_for_metadata_category_all; ch_metadata_for_metadata_category_pairwise; ch_metadata_for_beta_diversity; ch_metadata_for_beta_diversity_ordination; ch_metadata_for_ancom; ch_metadata_for_ancom_tax; ch_metadata_for_ancom_asv }

params.untilQ2import = false

params.Q2imported = false
if (params.Q2imported) {
    params.skip_fastqc = true
    params.skip_multiqc = true
    //Set up channel
    Channel.fromFile("${params.Q2imported}")
           .into { ch_qiime_demux_import; ch_qiime_demux_vis; ch_qiime_demux_dada }
    params.keepIntermediates = true
} else {
    params.skip_multiqc = false
}

//Currently, fastqc doesnt work for multiple runs when sample names are identical. These names are encoded in the sequencing file itself.
if (params.multipleSequencingRuns) {
    params.skip_fastqc = true
} else {
    params.skip_fastqc = false
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

/*
 * Sanity check input values
 * need to be extended eventually
 */
if (!params.Q2imported && (!params.FW_primer || !params.RV_primer || !params.metadata || !params.reads)) {
    println "${params.Q2imported}"
    println "\nERROR: Missing required input --Q2imported OR --FW_primer / --RV_primer / --metadata\n"
    helpMessage()
    exit 1
}


// AWSBatch sanity checking
if(workflow.profile == 'awsbatch'){
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Check workDir/outdir paths to be S3 buckets if running on AWSBatch
// related: https://github.com/nextflow-io/nextflow/issues/813
if( workflow.profile == 'awsbatch') {
    if(!workflow.workDir.startsWith('s3:') || !params.outdir.startsWith('s3:')) exit 1, "Workdir or Outdir not on S3 - specify S3 Buckets for each to run on AWSBatch!"
}

// Header log info
// TODO lets test this too - need to add more stuff as well here
log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'

nf-core/ampliseq v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']  = 'nf-core/ampliseq'
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Reads']        = params.reads
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(workflow.profile == 'awsbatch'){
   summary['AWS Region'] = params.awsregion
   summary['AWS Queue'] = params.awsqueue
}
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

if( !params.trunclenf || !params.trunclenr ){
    if ( !params.untilQ2import ) log.info "\n######## WARNING: No DADA2 cutoffs were specified, therefore reads will be truncated where median quality drops below ${params.trunc_qmin}.\nThe chosen cutoffs do not account for required overlap for merging, therefore DADA2 might have poor merging efficiency or even fail.\n"
}

def create_workflow_summary(summary) {

    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-ampliseq-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/ampliseq Workflow Summary'
    section_href: 'https://github.com/nf-core/ampliseq'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}


/*
 * Parse software version numbers
 */
process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    multiqc --version > v_multiqc.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}


if (!params.Q2imported){

    /*
    * Create a channel for input read files
    */
    if(params.readPaths && params.reads == "data${params.extension}"){
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { ch_read_pairs; ch_read_pairs_fastqc; ch_read_pairs_name_check }

    } else if ( params.multipleSequencingRuns ) {
        //Get files
        Channel
            .fromFilePairs( params.reads + "/*" + params.extension, size: 2 )
            .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}/*${params.extension}\nNB: Path needs to be enclosed in quotes!" }
            .into { ch_extract_folders; ch_rename_key }

        //Get folder information
        ch_extract_folders
            .flatMap { key, files -> [files[0]] }
            .map { it.take(it.findLastIndexOf{"/"})[-1] }
            .unique()
            .into { ch_folders; ch_check_folders; ch_report_folders }

        //Report folders with sequencing files
        ch_report_folders
            .collect()
            .subscribe {
                String folders = it.toString().replace("[", "").replace("]","") 
                log.info "\nFound the folder(s) \"$folders\" containing sequencing read files matching \"${params.extension}\" in \"${params.reads}\".\n" }

        //Stop if folder count is 1
        ch_check_folders
            .count()
            .subscribe { if ( it == 1 ) exit 1, "Found only one folder with read data but \"--multipleSequencingRuns\" was specified. Please review data input." }

        //Add folder information to sequence files
        ch_rename_key
            .map { key, files -> [ key, files, (files[0].take(files[0].findLastIndexOf{"/"})[-1]) ] }
            .into { ch_read_pairs; ch_read_pairs_fastqc; ch_read_pairs_name_check }

        //Check if key follows regex "^[a-zA-Z0-9-]+_[a-zA-Z0-9-]+$"
        ch_read_pairs_name_check
            .map { key, files, folder -> [ key ] }
            .subscribe { 
                if ( !(it =~ /[a-zA-Z0-9-]+_[a-zA-Z0-9-]+/) ) exit 1, "files starting with $it dont match the QIIME2 input requirements \"[a-zA-Z0-9-]+_[a-zA-Z0-9-]+_L[0-9][0-9][0-9]_R{1,2}_001.fastq.gz\". There might be more, just stopped here. \nPlease follow input requirements outlined in the documentation." }
            
    } else {
        Channel
            .fromFilePairs( params.reads + params.extension, size: 2 )
            .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}${params.extension}\nNB: Path needs to be enclosed in quotes!" }
            .into { ch_read_pairs; ch_read_pairs_fastqc; ch_read_pairs_name_check }

        //Check if key follows regex "^[a-zA-Z0-9-]+_[a-zA-Z0-9-]+$"
        ch_read_pairs_name_check
            .map { key, files -> [ key ] }
            .subscribe { 
                if ( !(it =~ /[a-zA-Z0-9-]+_[a-zA-Z0-9-]+/) ) exit 1, "files starting with $it dont match the QIIME2 input requirements \"[a-zA-Z0-9-]+_[a-zA-Z0-9-]+_L[0-9][0-9][0-9]_R{1,2}_001.fastq.gz\". There might be more, just stopped here. \nPlease follow input requirements outlined in the documentation." }
    }

	/*
	 * fastQC
	 */
    if (!params.multipleSequencingRuns){
        process fastqc {
            tag "$pair_id"
            publishDir "${params.outdir}/fastQC", mode: 'copy',
            saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

            input:
            set pair_id, file(reads) from ch_read_pairs_fastqc

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
            tag "$folder-$pair_id"
            publishDir "${params.outdir}/fastQC", mode: 'copy',
            saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

            input:
            set pair_id, file(reads), folder from ch_read_pairs_fastqc

            output:
            file "*_fastqc.{zip,html}" into ch_fastqc_results

            when:
            !params.skip_fastqc

            script: 
            """
            fastqc -q ${reads}
            """
        }
    }

	/*
	 * Trim each read-pair with cutadapt
	 */
    if (!params.multipleSequencingRuns){
        process trimming {
            tag "$pair_id"  
            publishDir "${params.outdir}/trimmed", mode: 'copy',
                saveAs: {filename -> 
                if (filename.indexOf(".gz") == -1) "logs/$filename"
                else if(params.keepIntermediates) filename 
                else null}
        
            input:
            set pair_id, file(reads) from ch_read_pairs
        
            output:
            file "trimmed/*.*" into ch_fastq_trimmed
            file "cutadapt_log_*.txt" into ch_fastq_cutadapt_log

            script:
            if( params.retain_untrimmed == false ){ 
                discard_untrimmed = "--discard-untrimmed"
            } else {
                discard_untrimmed = ""
            }
        
            """
            mkdir -p trimmed
            cutadapt -g ${params.FW_primer} -G ${params.RV_primer} $discard_untrimmed \
                -o trimmed/${reads[0]} -p trimmed/${reads[1]} \
                ${reads[0]} ${reads[1]} > cutadapt_log_${pair_id}.txt
            """
        }
    } else {
        process trimming_multi {
            tag "$folder-$pair_id"  
            publishDir "${params.outdir}/trimmed", mode: 'copy',
                saveAs: {filename -> 
                if (filename.indexOf(".gz") == -1) "logs/$filename"
                else if(params.keepIntermediates) filename 
                else null}
        
            input:
            set pair_id, file(reads), folder from ch_read_pairs
        
            output:
            file "trimmed/*.*" into ch_fastq_trimmed
            file "cutadapt_log_*.txt" into ch_fastq_cutadapt_log

            script:
            if( params.retain_untrimmed == false ){ 
                discard_untrimmed = "--discard-untrimmed"
            } else {
                discard_untrimmed = ""
            }
        
            """
            mkdir -p trimmed
            cutadapt -g ${params.FW_primer} -G ${params.RV_primer} $discard_untrimmed \
                -o trimmed/$folder-${reads[0]} -p trimmed/$folder-${reads[1]} \
                ${reads[0]} ${reads[1]} > cutadapt_log_${pair_id}.txt
            """
        }
	}

	/*
	 * multiQC
	 */
	process multiqc {
	    publishDir "${params.outdir}/MultiQC", mode: 'copy'

	    input:
	    file ('fastqc/*') from ch_fastqc_results.collect()
        file ('cutadapt/logs/*') from ch_fastq_cutadapt_log.collect()

	    output:
	    file "*multiqc_report.html" into multiqc_report
	    file "*_data"

	    when:
	    !params.skip_multiqc

	    script:
	    """
	    multiqc --force --interactive .
	    """
	}

    /*
    * Import trimmed files into QIIME2 artefact
    */
    if (!params.multipleSequencingRuns){
        process qiime_import {
            publishDir "${params.outdir}/demux", mode: 'copy', 
            saveAs: { filename -> 
                params.keepIntermediates ? filename : null
                params.untilQ2import ? filename : null }

            input:
            file(trimmed) from ch_fastq_trimmed.collect()
            env MATPLOTLIBRC from ch_mpl_for_qiime_import

            output:
            file "demux.qza" into (ch_qiime_demux_import, ch_qiime_demux_vis, ch_qiime_demux_dada)

            when:
            !params.Q2imported
        
            """
            qiime tools import  \
            --type 'SampleData[PairedEndSequencesWithQuality]'  \
            --input-path .  \
            --source-format CasavaOneEightSingleLanePerSampleDirFmt  \
            --output-path demux.qza
            """
        }
    } else {
        process qiime_import_multi {
            tag "${folders}"

            publishDir "${params.outdir}", mode: 'copy', 
            saveAs: {params.keepIntermediates ? filename : null}

            input:
            file(trimmed) from ch_fastq_trimmed.collect()
            val(folders) from ch_folders.collect()
            env MATPLOTLIBRC from ch_mpl_for_qiime_import

            output:
            file "*demux.qza" into (ch_qiime_demux_import, ch_qiime_demux_vis, ch_qiime_demux_dada) mode flatten

            when:
            !params.Q2imported

            script:
            """
            for folder in $folders
            do
                folder=\"\${folder//[],[]}\"
                mkdir \$folder
                mv \$folder-* \$folder/
                qiime tools import \
                --type 'SampleData[PairedEndSequencesWithQuality]' \
                --input-path \$folder \
                --source-format CasavaOneEightSingleLanePerSampleDirFmt \
                --output-path \$folder-demux.qza
            done
            """
        }
	}
    ch_qiime_demux_vis
        .combine( ch_mpl_for_demux_visualize )
        .set{ ch_qiime_demux_visualisation }

}


/*
 * Download, unpack, extract and train classifier
 * Use "--dereplication 90" for testing and "--dereplication 99" for real datasets
 * Requirements with "--dereplication 99": 1 core (seems not to scale with more?), ~35 Gb mem, ~2:15:00 walltime
 */

if( !params.classifier ){
    Channel.fromPath("${params.reference_database}")
        .set { ch_ref_database }

	process make_SILVA_132_16S_classifier {
        publishDir "${params.outdir}/DB/", mode: 'copy', 
        saveAs: {filename -> 
            if (filename.indexOf("${params.FW_primer}-${params.RV_primer}-${params.dereplication}-classifier.qza") == 0) filename
            else if(params.keepIntermediates) filename 
            else null}

        input:
        file database from ch_ref_database
        env MATPLOTLIBRC from ch_mpl_for_make_classifier

	    output:
	    file("${params.FW_primer}-${params.RV_primer}-${params.dereplication}-classifier.qza") into ch_qiime_classifier
        file("*.qza")
        stdout message_classifier_removeHash

	    when:
	    !params.onlyDenoising

	    script:
	  
	    """
	    unzip -qq $database

        fasta=\"SILVA_132_QIIME_release/rep_set/rep_set_16S_only/${params.dereplication}/silva_132_${params.dereplication}_16S.fna\"
        taxonomy=\"SILVA_132_QIIME_release/taxonomy/16S_only/${params.dereplication}/consensus_taxonomy_7_levels.txt\"

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
		--source-format HeaderlessTSVTaxonomyFormat \
		--input-path \$taxonomy \
		--output-path ref-taxonomy-${params.dereplication}.qza

	    #Extract sequences based on primers
	    qiime feature-classifier extract-reads \
		--i-sequences ref-seq-${params.dereplication}.qza \
		--p-f-primer ${params.FW_primer} \
		--p-r-primer ${params.RV_primer} \
		--o-reads ${params.FW_primer}-${params.RV_primer}-${params.dereplication}-ref-seq.qza \
        --quiet

	    #Train classifier
	    qiime feature-classifier fit-classifier-naive-bayes \
		--i-reference-reads ${params.FW_primer}-${params.RV_primer}-${params.dereplication}-ref-seq.qza \
		--i-reference-taxonomy ref-taxonomy-${params.dereplication}.qza \
		--o-classifier ${params.FW_primer}-${params.RV_primer}-${params.dereplication}-classifier.qza \
        --quiet
	    """
	}
    message_classifier_removeHash
        .subscribe { log.info it }
} else {
    Channel.fromPath("${params.classifier}")
           .set { ch_qiime_classifier }
}

/*
 * Import trimmed files into QIIME2 artefact
 */
if( !params.Q2imported ){
	process qiime_demux_visualize {
        tag "${demux.baseName}"
        publishDir "${params.outdir}", mode: 'copy'

	    input:
	    set file(demux), env(MATPLOTLIBRC) from ch_qiime_demux_visualisation

	    output:
        file("${demux.baseName}/*-seven-number-summaries.csv") into csv_demux
        file("${demux.baseName}/*")
	  
	    """
	    qiime demux summarize \
		--i-data $demux \
		--o-visualization ${demux.baseName}.qzv

	    qiime tools export ${demux.baseName}.qzv --output-dir ${demux.baseName}
	    """
	}
} else {
	process qiime_importdemux_visualize { 
        publishDir "${params.outdir}", mode: 'copy'

	    input:
        env MATPLOTLIBRC from ch_mpl_for_demux_visualize

	    output:
	    file("demux/*-seven-number-summaries.csv") into csv_demux
        file("demux/*")
	  
	    """
	    qiime demux summarize \
		--i-data ${params.Q2imported} \
		--o-visualization demux.qzv

	    qiime tools export demux.qzv --output-dir demux
	    """
	}
}


/*
 * Determine params.trunclenf and params.trunclenr where the median quality value drops below params.trunc_qmin
 * "Warning massage" should be printed but interferes with output: stdout
 * "Error and exit if too short" could be done in the python script itself?
 */
process dada_trunc_parameter { 

    input:
    file summary_demux from csv_demux 

    output:
    stdout dada_trunc

    when:
    !params.untilQ2import

    script:
    if( !params.trunclenf || !params.trunclenr ){
	    """
        dada_trunc_parameter.py ${summary_demux[0]} ${summary_demux[1]} ${params.trunc_qmin}
	    """
    }
    else
	    """
	    printf "${params.trunclenf},${params.trunclenr}"
	    """
}

if (params.multipleSequencingRuns){
    //find minimum dada truncation values
    dada_trunc
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
 * Find ASVs with DADA2 for single sequencing run
 */
if (!params.multipleSequencingRuns){
    process dada_single {
        tag "$trunc"
        publishDir "${params.outdir}", mode: 'copy',
            saveAs: {filename -> 
                     if (filename.indexOf("stats.tsv") > 0)                     "abundance_table/unfiltered/dada_stats.tsv"
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
        val trunc from dada_trunc
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
        IFS=',' read -r -a trunclen <<< \"$trunc\"

        #denoise samples with DADA2 and produce
        qiime dada2 denoise-paired  \
        --i-demultiplexed-seqs $demux  \
        --p-trunc-len-f \${trunclen[0]} \
        --p-trunc-len-r \${trunclen[1]} \
        --p-n-threads 0  \
        --o-table table.qza  \
        --o-representative-sequences rep-seqs.qza  \
        --o-denoising-stats stats.qza \
        --verbose \
        >dada_report.txt

        #produce dada2 stats "dada_stats/stats.tsv"
        qiime tools export stats.qza \
        --output-dir dada_stats

        #produce raw count table in biom format "table/feature-table.biom"
        qiime tools export table.qza  \
        --output-dir table

        #produce raw count table
        biom convert -i table/feature-table.biom \
        -o table/feature-table.tsv  \
        --to-tsv

        #produce representative sequence fasta file
        qiime feature-table tabulate-seqs  \
        --i-data rep-seqs.qza  \
        --o-visualization rep-seqs.qzv
        qiime tools export rep-seqs.qzv  \
        --output-dir unfiltered

        #convert to relative abundances
        qiime feature-table relative-frequency \
        --i-table table.qza \
        --o-relative-frequency-table relative-table-ASV.qza

        #export to biom
        qiime tools export relative-table-ASV.qza \
        --output-dir rel-table

        #convert to tab seperated text file
        biom convert \
        -i rel-table/feature-table.biom \
        -o table/rel-feature-table.tsv --to-tsv
        """
    }
} else {
    process dada_multi {
        tag "${demux.baseName} $trunclenf $trunclenr"

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
        #denoise samples with DADA2 and produce
        qiime dada2 denoise-paired  \
        --i-demultiplexed-seqs $demux  \
        --p-trunc-len-f ${trunclenf} \
        --p-trunc-len-r ${trunclenr} \
        --p-n-threads 0  \
        --o-table ${demux.baseName}-table.qza  \
        --o-representative-sequences ${demux.baseName}-rep-seqs.qza  \
        --o-denoising-stats ${demux.baseName}-stats.qza \
        --verbose \
        >${demux.baseName}-report.txt

        #produce dada2 stats "${demux.baseName}-dada_stats/stats.tsv"
        qiime tools export ${demux.baseName}-stats.qza \
        --output-dir ${demux.baseName}-dada_stats
        cp ${demux.baseName}-dada_stats/stats.tsv ${demux.baseName}-stats.tsv
        """
    }

    process dada_merge {
        tag "$tables"
        publishDir "${params.outdir}", mode: 'copy',
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
        for (table in tables) { TABLES+= " --i-tables $table" }
        for (repseq in repseqs) { REPSEQ+= " --i-data $repseq" }
        for (stat in stats) { STAT+= " $stat" }
        for (report in reports) { REPORT+= " $report" }
        """
        #concatenate tables
        #merge files
        qiime feature-table merge \
	        $TABLES \
	        --o-merged-table table.qza \
	        --quiet
        qiime feature-table merge-seqs \
	        $REPSEQ \
	        --o-merged-data rep-seqs.qza \
	        --quiet
        cat $STAT >stats.tsv
        cat $REPORT >dada_report.txt

        #produce raw count table in biom format "table/feature-table.biom"
        qiime tools export table.qza  \
        --output-dir table

        #produce raw count table
        biom convert -i table/feature-table.biom \
        -o table/feature-table.tsv  \
        --to-tsv

        #produce representative sequence fasta file
        qiime feature-table tabulate-seqs  \
        --i-data rep-seqs.qza  \
        --o-visualization rep-seqs.qzv
        qiime tools export rep-seqs.qzv  \
        --output-dir unfiltered

        #convert to relative abundances
        qiime feature-table relative-frequency \
        --i-table table.qza \
        --o-relative-frequency-table relative-table-ASV.qza

        #export to biom
        qiime tools export relative-table-ASV.qza \
        --output-dir rel-table

        #convert to tab seperated text file
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
    publishDir "${params.outdir}", mode: 'copy',
        saveAs: {filename -> 
            if (filename == "taxonomy/taxonomy.tsv") filename
            else if (filename == "taxonomy.qza") "taxonomy/$filename"}

    input:
    file repseq from ch_qiime_repseq_raw_for_classifier
    file trained_classifier from ch_qiime_classifier
    env MATPLOTLIBRC from ch_mpl_classifier

    output:
    file("taxonomy.qza") into (ch_qiime_taxonomy_for_filter,ch_qiime_taxonomy_for_relative_abundance_reduced_taxa,ch_qiime_taxonomy_for_barplot,ch_qiime_taxonomy_for_ancom)
    file("taxonomy/taxonomy.tsv") into ch_tsv_taxonomy

  
    """
    qiime feature-classifier classify-sklearn  \
	--i-classifier $trained_classifier  \
	--p-n-jobs ${task.cpus}  \
	--i-reads $repseq  \
	--o-classification taxonomy.qza  \
	--verbose

    qiime metadata tabulate  \
	--m-input-file taxonomy.qza  \
	--o-visualization taxonomy.qzv  \
	--verbose

    #produce "taxonomy/taxonomy.tsv"
    qiime tools export taxonomy.qza  \
	--output-dir taxonomy

    qiime tools export taxonomy.qzv  \
	--output-dir taxonomy
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

    	publishDir "${params.outdir}", mode: 'copy',
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
        if ! [ \"${params.exclude_taxa}\" = \"none\" ]; then
            #filter sequences
            qiime taxa filter-seqs \
            --i-sequences $repseq \
            --i-taxonomy $taxonomy \
            --p-exclude ${params.exclude_taxa} --p-mode contains \
            --o-filtered-sequences tax_filtered-sequences.qza

            #filter abundance table
            qiime taxa filter-table \
            --i-table $table \
            --i-taxonomy $taxonomy \
            --p-exclude ${params.exclude_taxa} --p-mode contains \
            --o-filtered-table tax_filtered-table.qza

            filtered_table="tax_filtered-table.qza"
            filtered_sequences="tax_filtered-sequences.qza"
        else
            filtered_table=$table
            filtered_sequences=$repseq
        fi

        qiime feature-table filter-features \
        --i-table \$filtered_table \
        --p-min-frequency $minfrequency \
        --p-min-samples $minsamples \
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
	publishDir "${params.outdir}", mode: 'copy',
        saveAs: {filename -> 
                 if (filename.indexOf("table/feature-table.biom") == 0)  "abundance_table/filtered/feature-table.biom"
            else if (filename.indexOf("table/feature-table.tsv") == 0)   "abundance_table/filtered/feature-table.tsv"
            else if (filename.indexOf("filtered/*"))                     "representative_sequences/$filename"
            else null}   

    input:
    file table from ch_qiime_table_for_filtered_dada_output
    file repseq from ch_qiime_repseq_for_dada_output
    env MATPLOTLIBRC from ch_mpl_for_export_dada_output

    output:
    file("filtered/sequences.fasta") into ch_fasta_repseq
    file("table/feature-table.tsv") into (ch_tsv_table_for_alpha_rarefaction,ch_tsv_table_for_report_filter_stats,ch_tsv_table_for_diversity_core)
    file("table/feature-table.biom")
    file("filtered/*")

    """
    #produce raw count table in biom format "table/feature-table.biom"
    qiime tools export $table  \
	--output-dir table

    #produce raw count table "table/feature-table.tsv"
    biom convert -i table/feature-table.biom \
	-o table/feature-table.tsv  \
	--to-tsv

    #produce representative sequence fasta file "${params.outdir}/representative_sequences/sequences.fasta"
    qiime feature-table tabulate-seqs  \
	--i-data $repseq  \
	--o-visualization rep-seqs.qzv
    qiime tools export rep-seqs.qzv  \
	--output-dir filtered
    """
}

/*
 * Report stats after taxa filtering
 */
process report_filter_stats { 
    publishDir "${params.outdir}/abundance_table/filtered", mode: 'copy'     

    input:
    file 'unfiltered_table' from ch_tsv_table_raw
    file 'filtered_table' from ch_tsv_table_for_report_filter_stats

    output:
    file("count_table_filter_stats.csv") into ch_csv_filter_stats
    
    """
    count_table_filter_stats.py unfiltered_table filtered_table
    """
}

/*
 * Export relative abundance tables on ASV level
 */
process RelativeAbundanceASV { 
    publishDir "${params.outdir}/rel_abundance_tables", mode: 'copy'    

    input:
    file table from ch_qiime_table_for_relative_abundance_asv
    env MATPLOTLIBRC from ch_mpl_for_relasv

    output:
    file("rel-table-ASV.tsv") into ch_tsv_relASV_table

    when:
    !params.skip_abundance_tables

    """
    ##onASV level

    #convert to relative abundances
    qiime feature-table relative-frequency \
	--i-table $table \
	--o-relative-frequency-table relative-table-ASV.qza

    #export to biom
    qiime tools export relative-table-ASV.qza --output-dir relative-table-ASV

    #convert to tab seperated text file "${params.outdir}/rel-table-ASV.tsv"
    biom convert -i relative-table-ASV/feature-table.biom \
	-o rel-table-ASV.tsv --to-tsv
    """
}


/*
 * Export relative abundance tables based on taxonomic levels
 */
process RelativeAbundanceReducedTaxa { 
    publishDir "${params.outdir}/rel_abundance_tables", mode: 'copy'    

    input:
    file table from ch_qiime_table_for_relative_abundance_reduced_taxa
    file taxonomy from ch_qiime_taxonomy_for_relative_abundance_reduced_taxa
    env MATPLOTLIBRC from ch_mpl_for_relreducetaxa

    output:
    file("*.tsv")

    when:
    !params.skip_abundance_tables && !params.skip_taxonomy

    """
    ##on several taxa level

    array=( 2 3 4 5 6 7 )
    for i in \${array[@]}
    do
	    #collapse taxa
	    qiime taxa collapse \
		    --i-table $table \
		    --i-taxonomy $taxonomy \
		    --p-level \$i \
		    --o-collapsed-table table-\$i.qza
	    #convert to relative abundances
	    qiime feature-table relative-frequency \
		    --i-table table-\$i.qza \
		    --o-relative-frequency-table relative-table-\$i.qza
	    #export to biom
	    qiime tools export relative-table-\$i.qza \
		    --output-dir relative-table-\$i
	    #convert to tab seperated text file
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
    publishDir "${params.outdir}", mode: 'copy'    

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
    qiime taxa barplot  \
	--i-table $table  \
	--i-taxonomy $taxonomy  \
	--m-metadata-file $metadata  \
	--o-visualization taxa-bar-plots.qzv  \
	--verbose

    qiime tools export taxa-bar-plots.qzv  \
	--output-dir barplot
    """
}

/*
 * Produce a rooted tree
 */
process tree { 
    publishDir "${params.outdir}", mode: 'copy',
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
    qiime alignment mafft \
	--i-sequences $repseq \
	--o-alignment aligned-rep-seqs.qza \
	--p-n-threads ${params.tree_cores}

    qiime alignment mask \
	--i-alignment aligned-rep-seqs.qza \
	--o-masked-alignment masked-aligned-rep-seqs.qza

    qiime phylogeny fasttree \
	--i-alignment masked-aligned-rep-seqs.qza \
	--p-n-threads ${params.tree_cores} \
	--o-tree unrooted-tree.qza

    qiime phylogeny midpoint-root \
	--i-tree unrooted-tree.qza \
	--o-rooted-tree rooted-tree.qza

    qiime tools export rooted-tree.qza  \
	--output-dir phylogenetic_tree
    """
}


/*
 * Alpha-rarefaction
 */
process alpha_rarefaction { 
    publishDir "${params.outdir}", mode: 'copy'    

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
    maxdepth=\$(count_table_minmax_reads.py $stats maximum 2>&1)

    #check values
    if [ \"\$maxdepth\" -gt \"75000\" ]; then maxdepth=\"75000\"; fi
    if [ \"\$maxdepth\" -gt \"5000\" ]; then maxsteps=\"250\"; else maxsteps=\$((maxdepth/20)); fi

    qiime diversity alpha-rarefaction  \
	--i-table $table  \
	--i-phylogeny $tree  \
	--p-max-depth \$maxdepth  \
	--m-metadata-file $metadata  \
	--p-steps \$maxsteps  \
	--p-iterations 10  \
	--o-visualization alpha-rarefaction.qzv

    qiime tools export alpha-rarefaction.qzv  \
	--output-dir alpha-rarefaction
    """
}

/*
 * Combine abundances, sequences and taxonomic classification into one table with R
 */
process combinetable { 
    publishDir "${params.outdir}/rel_abundance_tables", mode: 'copy'    

    input:
    file TABLE from ch_tsv_relASV_table
    file SEQ from ch_fasta_repseq
    file TAXONOMY from ch_tsv_taxonomy

    output:
    file("qiime2_ASV_table.csv")

    when:
    !params.skip_abundance_tables && !params.skip_taxonomy

    """
    combineTable.r $TABLE $SEQ $TAXONOMY
    """
}

/*
 * Compute diversity matrices
 */
process diversity_core { 
    publishDir "${params.outdir}/diversity_core", mode: 'copy',
    saveAs: {params.keepIntermediates ? filename : null}

    input:
    file metadata from ch_metadata_for_diversity_core
    file table from ch_qiime_table_for_diversity_core
    file tree from ch_qiime_tree_for_diversity_core
    file stats from ch_tsv_table_for_diversity_core
    env MATPLOTLIBRC from ch_mpl_for_diversity_core

    output:
    file("core/*_pcoa_results.qza") into (qiime_diversity_core_for_beta_diversity_ordination) mode flatten
    file("core/*_vector.qza") into qiime_diversity_core_for_alpha_diversity mode flatten
    file("core/*_distance_matrix.qza") into qiime_diversity_core_for_beta_diversity mode flatten
    stdout rarefaction_depth

    when:
    !params.skip_diversity_indices

    """
    mindepth=\$(count_table_minmax_reads.py $stats minimum 2>&1)

    if [ \"\$mindepth\" -gt \"10000\" ]; then echo \"\nUse the sampling depth of \$mindepth for rarefaction\" ; fi
    if [ \"\$mindepth\" -lt \"10000\" -a \"\$mindepth\" -gt \"5000\" ]; then echo \"\n######## WARNING! The sampling depth of \$mindepth is quite small for rarefaction!\" ; fi
    if [ \"\$mindepth\" -lt \"5000\" -a \"\$mindepth\" -gt \"1000\" ]; then echo \"\n######## WARNING! The sampling depth of \$mindepth is very small for rarefaction!\" ; fi
    if [ \"\$mindepth\" -lt \"1000\" ]; then echo \"\n######## ERROR! The sampling depth of \$mindepth seems too small for rarefaction!\" ; fi
    
    qiime diversity core-metrics-phylogenetic \
	--m-metadata-file $metadata \
	--i-phylogeny $tree \
	--i-table $table \
	--p-sampling-depth \$mindepth \
	--output-dir core \
	--p-n-jobs ${params.diversity_cores} \
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
    stdout into (meta_category_all,meta_category_all_for_ancom)

    when:
    !params.skip_ancom || !params.skip_diversity_indices
    when:
    !params.untilQ2import && !params.onlyDenoising

    script:
    if( !params.metadata_category )
	    """
	    metadataCategory.r $metadata
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
    stdout meta_category_pairwise

    when:
    !params.skip_diversity_indices

    """
    metadataCategoryPairwise.r $metadata
    """
}

/*
 * Combine channels for diversity analysis
 */

ch_metadata_for_alpha_diversity
    .combine( qiime_diversity_core_for_alpha_diversity )
    .combine( ch_mpl_for_alpha_diversity )
    .set{ ch_for_alpha_diversity }
ch_metadata_for_beta_diversity
    .combine( qiime_diversity_core_for_beta_diversity )
    .combine( meta_category_pairwise )
    .combine( ch_mpl_for_beta_diversity )
    .set{ ch_for_beta_diversity }
ch_metadata_for_beta_diversity_ordination
    .combine( qiime_diversity_core_for_beta_diversity_ordination )
    .combine( ch_mpl_for_beta_diversity_ord )
    .set{ ch_for_beta_diversity_ordination }
    


/*
 * Compute alpha diversity indices
 */
process alpha_diversity { 
    tag "${core.baseName}"
    publishDir "${params.outdir}", mode: 'copy'    

    input:
    set file(metadata), file(core), env(MATPLOTLIBRC) from ch_for_alpha_diversity

    output:
    file("alpha-diversity/*") into qiime_alphadiversity

    """
	qiime diversity alpha-group-significance \
        --i-alpha-diversity ${core} \
        --m-metadata-file ${metadata} \
        --o-visualization ${core.baseName}-vis.qzv
	qiime tools export ${core.baseName}-vis.qzv \
        --output-dir "alpha-diversity/${core.baseName}"
    """
}


/*
 * Compute beta diversity indices
 */
process beta_diversity { 
    tag "${core.baseName}"
    publishDir "${params.outdir}", mode: 'copy'     

    input:
    set file(meta), file(core), val(category), env(MATPLOTLIBRC) from ch_for_beta_diversity

    output:
    file "beta-diversity/*"

    """
    IFS=',' read -r -a metacategory <<< \"$category\"

	for j in \"\${metacategory[@]}\"
	do
		qiime diversity beta-group-significance \
		    --i-distance-matrix $core \
		    --m-metadata-file $meta \
		    --m-metadata-column \"\$j\" \
		    --o-visualization ${core.baseName}-\$j.qzv \
		    --p-pairwise
		qiime tools export ${core.baseName}-\$j.qzv \
		    --output-dir beta-diversity/${core.baseName}-\$j
    done
    """
}

/*
 * Compute beta diversity ordination
 */
process beta_diversity_ordination { 
    tag "${core.baseName}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    set file(metadata), file(core), env(MATPLOTLIBRC) from ch_for_beta_diversity_ordination

    output:
    file("beta-diversity/*")

    """
	qiime emperor plot \
        --i-pcoa ${core} \
        --m-metadata-file ${metadata} \
        --o-visualization ${core.baseName}-vis.qzv
	qiime tools export ${core.baseName}-vis.qzv \
        --output-dir beta-diversity/${core.baseName}-PCoA
    """
}


/*
 * Differential abundance analysis with ANCOM
 */
process prepare_ancom { 
    tag "$meta"

    publishDir "${params.outdir}/ancom", mode: 'copy', 
    saveAs: {params.keepIntermediates ? filename : null}   

    input:
    file metadata from ch_metadata_for_ancom
    file table from ch_qiime_table_for_ancom
    val meta from meta_category_all_for_ancom
    env MATPLOTLIBRC from ch_mpl_for_ancom

    output:
    file("*.qza") into (ch_meta_tables_tax, ch_meta_tables_asv) mode flatten

    when:
    !params.skip_ancom

    """
    IFS=',' read -r -a metacategory <<< \"$meta\"

    #remove samples that do not have any value
    for j in \"\${metacategory[@]}\"
    do
	    qiime feature-table filter-samples \
		    --i-table $table \
		    --m-metadata-file $metadata \
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
    tag "${table.baseName}-level$taxlevel"

    publishDir "${params.outdir}", mode: 'copy'    

    input:
    set file(table), val(taxlevel), file(taxonomy), file(metadata), env(MATPLOTLIBRC) from ch_for_ancom_tax

    output:
    file("ancom/*")

    when:
    !params.skip_ancom

    """
	qiime taxa collapse \
	    --i-table $table \
		--i-taxonomy $taxonomy \
		--p-level $taxlevel \
		--o-collapsed-table lvl$taxlevel-$table
	qiime composition add-pseudocount \
		--i-table lvl$taxlevel-$table \
		--o-composition-table comp-lvl$taxlevel-$table
	qiime composition ancom \
		--i-table comp-lvl$taxlevel-$table \
		--m-metadata-file $metadata \
		--m-metadata-column ${table.baseName} \
		--o-visualization comp-lvl$taxlevel-${table.baseName}.qzv
	qiime tools export comp-lvl$taxlevel-${table.baseName}.qzv \
	    --output-dir ancom/Category-${table.baseName}-level-$taxlevel
    """
}

/*
 * Differential abundance analysis with ANCOM on ASV level
 */
process ancom_asv { 
    tag "${table.baseName}"

    publishDir "${params.outdir}", mode: 'copy' 

    input:
    set file(table), file(metadata), env(MATPLOTLIBRC) from ch_for_ancom_asv 

    output:
    file("ancom/*") 

    """
    qiime composition add-pseudocount \
		--i-table $table \
		--o-composition-table comp-$table
	qiime composition ancom \
	    --i-table comp-$table \
	    --m-metadata-file $metadata \
	    --m-metadata-column ${table.baseName} \
	    --o-visualization comp-${table.baseName}.qzv
	qiime tools export comp-${table.baseName}.qzv \
	    --output-dir ancom/Category-${table.baseName}-ASV
    """
}

/*
 * STEP 3 - Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/Documentation", mode: 'copy'

    input:
    file output_docs from ch_output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r $output_docs results_description.html
    """
}



/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/ampliseq] Successful: $workflow.runName"
    if(!workflow.success){
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
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[nf-core/ampliseq] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nf-core/ampliseq] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/Documentation/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[nf-core/ampliseq] Pipeline Complete"

}
