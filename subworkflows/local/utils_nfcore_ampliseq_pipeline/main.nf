//
// Subworkflow with functionality specific to the nf-core/ampliseq pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFVALIDATION_PLUGIN } from '../../nf-core/utils_nfvalidation_plugin'
include { paramsSummaryMap          } from 'plugin/nf-validation'
include { fromSamplesheet           } from 'plugin/nf-validation'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { dashedLine                } from '../../nf-core/utils_nfcore_pipeline'
include { nfCoreLogo                } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { workflowCitation          } from '../../nf-core/utils_nfcore_pipeline'

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    help              // boolean: Display help text
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    pre_help_text = nfCoreLogo(monochrome_logs)
    post_help_text = '\n' + workflowCitation() + '\n' + dashedLine(monochrome_logs)
    def String workflow_command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"
    UTILS_NFVALIDATION_PLUGIN (
        help,
        workflow_command,
        pre_help_text,
        post_help_text,
        validate_params,
        "nextflow_schema.json"
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()
    // Check that keys for reference databases are valid
    if ( params.dada_ref_taxonomy && !params.skip_taxonomy && !params.skip_dada_taxonomy ) {
        dadareftaxonomyExistsError()
    }
    if ( params.sintax_ref_taxonomy && !params.skip_taxonomy ) {
        sintaxreftaxonomyExistsError()
    }
    if ( (params.qiime_ref_taxonomy || params.qiime_ref_tax_custom) && !params.skip_taxonomy && !params.classifier ) {
        qiimereftaxonomyExistsError()
    }
    if ( params.kraken2_ref_taxonomy  && !params.skip_taxonomy ) {
        kraken2reftaxonomyExistsError()
    }
    if ( params.sidle_ref_taxonomy && !params.skip_taxonomy ) {
        sidlereftaxonomyExistsError()
    }

    emit:
    versions    = ch_versions
}

/*
========================================================================================
    SUBWORKFLOW FOR PIPELINE COMPLETION
========================================================================================
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:

    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(summary_params, email, email_on_fail, plaintext_email, outdir, monochrome_logs, multiqc_report.toList())
        }

        completionSummary(monochrome_logs)

        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    if ( !params.input && !params.input_fasta && !params.input_folder ) {
        error("Missing input declaration: One of `--input`, `--input_fasta`, `--input_folder` is required.")
    }

    if ( !params.multiregion && !params.input_fasta && (!params.FW_primer || !params.RV_primer) && !params.skip_cutadapt ) {
        error("Incompatible parameters: `--FW_primer` and `--RV_primer` are required for primer trimming. If primer trimming is not needed, use `--skip_cutadapt`.")
    }

    if ( params.pacbio || params.iontorrent || params.single_end ) {
        if (params.trunclenr) { log.warn "Unused parameter: `--trunclenr` is ignored because the data is single end." }
    } else if (params.trunclenf && !params.trunclenr) {
        error("Invalid command: `--trunclenf` is set, but `--trunclenr` is not. Either both parameters `--trunclenf` and `--trunclenr` must be set or none.")
    } else if (!params.trunclenf && params.trunclenr) {
        error("Invalid command: `--trunclenr` is set, but `--trunclenf` is not. Either both parameters `--trunclenf` and `--trunclenr` must be set or none.")
    }

    if (!["pooled", "independent", "pseudo"].contains(params.sample_inference)) {
        def error_string = "Please set `--sample_inference` to one of the following:\n" +
            "\t-\"independent\" (lowest sensitivity and lowest resources),\n" +
            "\t-\"pseudo\" (balance between required resources and sensitivity),\n" +
            "\t-\"pooled\" (highest sensitivity and resources)."
        error(error_string)
    }

    if (params.double_primer && params.retain_untrimmed) {
        error("Incompatible parameters `--double_primer` and `--retain_untrimmed` cannot be set at the same time.")
    }

    if ( params.min_len_asv && params.max_len_asv && (params.min_len_asv > params.max_len_asv) ) {
        error("Incompatible parameters: `--min_len_asv` may not be greater than `--max_len_asv`.")
    }

    if ( params.skip_dada_quality && (params.trunclenf == null || params.trunclenr == null) ) {
        error("Incompatible parameters: `--skip_dada_quality` may not be used without setting `--trunclenf` and `--trunclenr`.")
    }

    if (params.tax_agglom_min > params.tax_agglom_max) {
        error("Incompatible parameters: `--tax_agglom_min` may not be greater than `--tax_agglom_max`.")
    }

    if (!params.dada_ref_tax_custom && params.dada_ref_tax_custom_sp) {
        error("Incompatible parameters: `--dada_ref_tax_custom_sp` requires `--dada_ref_tax_custom`.")
    }

    if (params.dada_ref_tax_custom && !params.dada_ref_tax_custom_sp && !params.skip_dada_addspecies) {
        error("Incompatible parameters: Either `--skip_dada_addspecies` or `--dada_ref_tax_custom_sp` is additionally required to `--dada_ref_tax_custom`.")
    }

    if (params.pplace_tree) {
        if (!params.pplace_aln) {
            error("Missing parameter: Phylogenetic placement requires in addition to `--pplace_tree` also `--pplace_aln`.")
        }
        if (!params.pplace_model) {
            error("Missing parameter: Phylogenetic placement requires in addition to `--pplace_tree` also `--pplace_model`.")
        }
    }

    if (params.dada_assign_taxlevels && params.sbdiexport && !params.sintax_ref_taxonomy) {
        error("Incompatible parameters: `--sbdiexport` expects specific taxonomics ranks (default) and therefore excludes modifying those using `--dada_assign_taxlevels`.")
    }

    if (params.skip_taxonomy && params.sbdiexport) {
        error("Incompatible parameters: `--sbdiexport` expects taxa annotation and therefore excludes `--skip_taxonomy`.")
    }

    if (params.skip_dada_taxonomy && params.sbdiexport) {
        if (!params.sintax_ref_taxonomy && (params.skip_qiime || (!params.qiime_ref_taxonomy && !params.qiime_ref_tax_custom))) {
            error("Incompatible parameters: `--sbdiexport` expects taxa annotation and therefore annotation with either DADA2, SINTAX, or QIIME2 is needed.")
        }
    }

    if ( (!params.FW_primer || !params.RV_primer) && (params.qiime_ref_taxonomy || params.qiime_ref_tax_custom) && !params.skip_qiime && !params.skip_taxonomy ) {
        error("Incompatible parameters: `--FW_primer` and `--RV_primer` are required for cutting the QIIME2 reference database to the amplicon sequences. Please specify primers or do not use `--qiime_ref_taxonomy`.")
    }

    if ( (!params.FW_primer || !params.RV_primer) && params.cut_dada_ref_taxonomy && !params.skip_taxonomy ) {
        error("Incompatible parameters: `--FW_primer` and `--RV_primer` are required for cutting the DADA2 reference database to the amplicon sequences. Please specify primers or do not use `--cut_dada_ref_taxonomy`.")
    }

    if ((params.qiime_ref_taxonomy || params.qiime_ref_tax_custom) && params.classifier) {
        error("Incompatible parameters: `--qiime_ref_taxonomy` and `--qiime_ref_tax_custom` will produce a classifier but `--classifier` points to a precomputed classifier, therefore, only use one of those.")
    }

    if (params.kraken2_ref_tax_custom && !params.kraken2_assign_taxlevels ) {
        error("Missing parameter: Taxonomic classification with a user provided database via `--kraken2_ref_tax_custom` requires `--kraken2_assign_taxlevels`")
    }

    if (params.filter_ssu && params.skip_barrnap) {
        error("Incompatible parameters: `--filter_ssu` cannot be used with `--skip_barrnap` because filtering for SSU's depends on barrnap.")
    }

    String[] sbdi_compatible_databases = ["coidb","coidb=221216","gtdb","gtdb=R08-RS214","gtdb=R07-RS207","gtdb=R06-RS202","gtdb=R05-RS95","midori2-co1","midori2-co1=gb250","pr2","pr2=5.0.0","pr2=4.14.0","pr2=4.13.0","rdp","rdp=18","sbdi-gtdb","sbdi-gtdb=R07-RS207-1","silva","silva=138","silva=132","unite-fungi","unite-fungi=9.0","unite-fungi=8.3","unite-fungi=8.2","unite-alleuk","unite-alleuk=9.0","unite-alleuk=8.3","unite-alleuk=8.2"]
    if (params.sbdiexport){
        if (params.sintax_ref_taxonomy ) {
            if (!Arrays.stream(sbdi_compatible_databases).anyMatch(entry -> params.sintax_ref_taxonomy.toString().equals(entry)) ) {
                error("Incompatible parameters: `--sbdiexport` does not work with the chosen database of `--sintax_ref_taxonomy` because the expected taxonomic levels do not match.")
            }
        } else if (!Arrays.stream(sbdi_compatible_databases).anyMatch(entry -> params.dada_ref_taxonomy.toString().equals(entry)) ) {
            error("Incompatible parameters: `--sbdiexport` does not work with the chosen database of `--dada_ref_taxonomy` because the expected taxonomic levels do not match.")
        }
    }

    if (params.addsh && !params.dada_ref_databases[params.dada_ref_taxonomy]["shfile"]) {
        def validDBs = ""
        for (db in params.dada_ref_databases.keySet()) {
            if (params.dada_ref_databases[db]["shfile"]) {
                validDBs += " " + db
            }
        }
        error("UNITE species hypothesis information is not available for the selected reference database, please use the option `--dada_ref_taxonomy` to select an appropriate database. Currently, the option `--addsh` can only be used together with the following UNITE reference databases:\n" + validDBs + ".")
    }

    if (params.addsh && params.cut_its == "none") {
        log.warn "Adding UNITE species hypothesis (SH) assignments is only feasible for ITS sequences. Please use option `--cut_its` to find ITS regions in the ASV sequences, unless the given sequences are already cut to the ITS region.\n"
    }

    // Error message for incompatible combination of --orf_start and --orf_end
    if ( params.orf_end && ( ( ( params.orf_end + 1 ) - params.orf_start ) % 3 != 0 ) ) {
        error("Incompatible parameters: The difference of  `--orf_end` and `--orf_start` must be a multiple of 3.")
    }

    // When multi-region analysis is used, some parameter combinations are required or not allowed:
    if ( params.multiregion ) {
        if ( !params.sidle_ref_taxonomy && !params.sidle_ref_tree_custom ) {
            log.warn "Missing parameter: Either use `--sidle_ref_taxonomy` or `--sidle_ref_tree_custom` to get (unified) taxonomic classifications"
        }
        if ( (params.dada_ref_tax_custom || params.dada_ref_taxonomy) && !params.skip_dada_taxonomy ) {
            error("Incompatible parameters: Multiple region analysis with `--multiregion` does not work with `--dada_ref_tax_custom`, `--dada_ref_taxonomy`")
        }
        if ( params.cut_its != "none" ) {
            error("Incompatible parameters: Multiple region analysis with `--multiregion` does not work with `--cut_its`")
        }
        if ( params.vsearch_cluster || params.filter_ssu || params.min_len_asv || params.max_len_asv || params.filter_codons ) {
            log.warn "Incompatible parameters: Multiple region analysis with `--multiregion` ignores any of `--vsearch_cluster`, `--filter_ssu`, `--min_len_asv`, `--max_len_asv`, `--filter_codons`, `--cut_its`"
        }
    }
}

//
// Prepare complement sequence - ultimately to make reverse complement primers
// Complement table taken from http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html
def makeComplement(seq) {
    def complements = [ A:'T', T:'A', U:'A', G:'C', C:'G', Y:'R', R:'Y', S:'S', W:'W', K:'M', M:'K', B:'V', D:'H', H:'D', V:'B', N:'N' ]
    def comp = seq.toUpperCase().collect { base -> complements[ base ] ?: 'X' }.join()
    return comp
}

//
// Exit pipeline if incorrect --dada_ref_taxonomy key provided
//
def dadareftaxonomyExistsError() {
    if (params.dada_ref_databases && params.dada_ref_taxonomy && !params.dada_ref_databases.containsKey(params.dada_ref_taxonomy)) {
        def error_string = "=============================================================================\n" +
            "  DADA2 reference database '${params.dada_ref_taxonomy}' not found in any config file provided to the pipeline.\n" +
            "  Currently, the available reference taxonomy keys for `--dada_ref_taxonomy` are:\n" +
            "  ${params.dada_ref_databases.keySet().join(", ")}\n" +
            "==================================================================================="
        error(error_string)
    }
}

//
// Exit pipeline if incorrect --sintax_ref_taxonomy key provided
//
def sintaxreftaxonomyExistsError() {
    if (params.sintax_ref_databases && params.sintax_ref_taxonomy && !params.sintax_ref_databases.containsKey(params.sintax_ref_taxonomy)) {
        def error_string = "=============================================================================\n" +
            "  SINTAX reference database '${params.sintax_ref_taxonomy}' not found in any config file provided to the pipeline.\n" +
            "  Currently, the available reference taxonomy keys for `--sintax_ref_taxonomy` are:\n" +
            "  ${params.sintax_ref_databases.keySet().join(", ")}\n" +
            "==================================================================================="
        error(error_string)
    }
}

//
// Exit pipeline if incorrect --qiime_ref_taxonomy key provided
//
def qiimereftaxonomyExistsError() {
    if (params.qiime_ref_databases && params.qiime_ref_taxonomy && !params.qiime_ref_databases.containsKey(params.qiime_ref_taxonomy)) {
        def error_string = "=============================================================================\n" +
            "  QIIME2 reference database '${params.qiime_ref_taxonomy}' not found in any config file provided to the pipeline.\n" +
            "  Currently, the available reference taxonomy keys for `--qiime_ref_taxonomy` are:\n" +
            "  ${params.qiime_ref_databases.keySet().join(", ")}\n" +
            "==================================================================================="
        error(error_string)
    }
}

//
// Exit pipeline if incorrect --kraken2_ref_taxonomy key provided
//
def kraken2reftaxonomyExistsError() {
    if (params.kraken2_ref_databases && params.kraken2_ref_taxonomy && !params.kraken2_ref_databases.containsKey(params.kraken2_ref_taxonomy)) {
        def error_string = "=============================================================================\n" +
            "  Kraken2 reference database '${params.kraken2_ref_taxonomy}' not found in any config file provided to the pipeline.\n" +
            "  Currently, the available reference taxonomy keys for `--kraken2_ref_taxonomy` are:\n" +
            "  ${params.kraken2_ref_databases.keySet().join(", ")}\n" +
            "==================================================================================="
        error(error_string)
    }
}

//
// Exit pipeline if incorrect --sidle_ref_taxonomy key provided
//
def sidlereftaxonomyExistsError() {
    if (params.sidle_ref_databases && params.sidle_ref_taxonomy && !params.sidle_ref_databases.containsKey(params.sidle_ref_taxonomy)) {
        def error_string = "=============================================================================\n" +
            "  Sidle reference database '${params.sidle_ref_taxonomy}' not found in any config file provided to the pipeline.\n" +
            "  Currently, the available reference taxonomy keys for `--sidle_ref_taxonomy` are:\n" +
            "  ${params.sidle_ref_databases.keySet().join(", ")}\n" +
            "==================================================================================="
        error(error_string)
    }
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ it.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
