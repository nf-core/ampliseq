#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/rrna-ampliseq
========================================================================================
 nf-core/rrna-ampliseq Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/rrna-ampliseq
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info"""
    =========================================
     nf-core/rrna-ampliseq v${manifest.pipelineVersion}
    =========================================
    
    Usage:

    The minimal command for running the pipeline is as follows:
    nextflow run qiime2.nf --reads "data/*_L001_R{1,2}_001.fastq.gz" --FW_primer GTGYCAGCMGCCGCGGTAA --RV_primer GGACTACNVGGGTWTCTAAT --metadata "$PWD/data/Metadata.tsv"

    The test command for running the pipeline is as follows:
    nextflow run qiime2.nf --minimumTest

    Required arguments:
      --reads [Path to folder]      Folder containing Casava 1.8 paired-end demultiplexed fastq files: *_L001_R{1,2}_001.fastq.gz
      --FW_primer [str]             Forward primer sequence
      --RV_primer [str]             Reverse primer sequence
      --metadata                    Absolute path to metadata sheet

    Filters:
      --exclude_taxa [str]          Comma seperated list of unwanted taxa (default: "mitochondria,chloroplast")
                                    To skip filtering use "none"

    Cutoffs:
      --retain_untrimmed            Cutadapt will retain untrimmed reads
      --trunclenf [int]             DADA2 read truncation value for forward strand
      --trunclenr [int]             DADA2 read truncation value for reverse strand
      --trunc_qmin [int]            If --trunclenf and --trunclenr are not set, 
                                    these values will be automatically determined using 
                                    this mean quality score (not preferred) (default: 25)

    References:                     If you have trained a compatible classifier before
      --classifier                  Path to QIIME2 classifier file (typically *-classifier.qza)

    Statistics:
      --metadata_category           Diversity indices will be calculated using these groupings in the metadata sheet,
                                    all suitable columns in the metadata sheet will be used if not specified.
                                    Suitable are columns which are categorical (not numerical) and have multiple  
                                    different values which are not all unique.

    Other options:
      --untilQ2import               Skip all steps after importing into QIIME2, used for visually choosing DADA2 parameter
      --Q2imported [Path]           Path to imported reads (e.g. "demux.qza"), used after visually choosing DADA2 parameter
      --onlyDenoising               Skip all steps after denoising, produce only sequences and abundance tables on ASV level

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
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.email = false
params.plaintext_email = false

multiqc_config = file(params.multiqc_config)
output_docs = file("$baseDir/docs/output.md")

// Defines all parameters that are independent of a test run
params.trunc_qmin = 25 //to calculate params.trunclenf and params.trunclenr automatically
params.trunclenf = false
params.trunclenr = false
params.metadata_category = false
params.qiimeimage = "$baseDir/qiime2_2018.6.simg"
params.tree_cores = 2
params.diversity_cores = 2
params.retain_untrimmed = false
params.exclude_taxa = "mitochondria,chloroplast"

/*
 * Defines pipeline steps
 */
params.untilQ2import = false

params.Q2imported = false
if (params.Q2imported) {
    params.skip_fastqc = true
    params.skip_multiqc = true
} else {
    params.skip_fastqc = false
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

// Validate inputs
if ( params.fasta ){
    fasta = file(params.fasta)
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
}
// AWSBatch sanity checking
if(workflow.profile == 'awsbatch'){
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
}
//
// NOTE - THIS IS NOT USED IN THIS PIPELINE, EXAMPLE ONLY
// If you want to use the above in a process, define the following:
//   input:
//   file fasta from fasta
//


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

/*
 * Create a channel for input read files
 */
 if(params.readPaths){
     if(params.singleEnd){
         Channel
             .from(params.readPaths)
             .map { row -> [ row[0], [file(row[1][0])]] }
             .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
             .into { read_files_fastqc; read_files_trimming }
     } else {
         Channel
             .from(params.readPaths)
             .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
             .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
             .into { read_files_fastqc; read_files_trimming }
     }
 } else {
     Channel
         .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
         .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
         .into { read_files_fastqc; read_files_trimming }
 }


// Header log info
log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'

nf-core/rrna-ampliseq v${manifest.pipelineVersion}"
======================================================="""
def summary = [:]
summary['Pipeline Name']  = 'nf-core/rrna-ampliseq'
summary['Pipeline Version'] = manifest.pipelineVersion
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Reads']        = params.reads
summary['Fasta Ref']    = params.fasta
summary['Data Type']    = params.singleEnd ? 'Single-End' : 'Paired-End'
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


def create_workflow_summary(summary) {

    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-rrna-ampliseq-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/rrna-ampliseq Workflow Summary'
    section_href: 'https://github.com/nf-core/rrna-ampliseq'
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
    echo $manifest.pipelineVersion > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    multiqc --version > v_multiqc.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}



/*
 * STEP 1 - FastQC
 */
process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from read_files_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    fastqc -q $reads
    """
}



/*
 * STEP 2 - MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config
    file ('fastqc/*') from fastqc_results.collect()
    file ('software_versions/*') from software_versions_yaml
    file workflow_summary from create_workflow_summary(summary)

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc -f $rtitle $rfilename --config $multiqc_config .
    """
}



/*
 * STEP 3 - Output Description HTML
 */
process output_documentation {
    tag "$prefix"
    publishDir "${params.outdir}/Documentation", mode: 'copy'

    input:
    file output_docs

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
    def subject = "[nf-core/rrna-ampliseq] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[nf-core/rrna-ampliseq] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = manifest.pipelineVersion
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
          log.info "[nf-core/rrna-ampliseq] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nf-core/rrna-ampliseq] Sent summary e-mail to $params.email (mail)"
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

    log.info "[nf-core/rrna-ampliseq] Pipeline Complete"

}
