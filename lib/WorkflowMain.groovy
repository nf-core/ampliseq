//
// This file holds several functions specific to the main.nf workflow in the nf-core/ampliseq pipeline
//

class WorkflowMain {

    //
    // Citation string for pipeline
    //
    public static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
            "* The pipeline publication\n" +
            "  https://doi.org/10.3389/fmicb.2020.550420\n\n" +
            "* The pipeline\n" +
            "  https://doi.org/10.5281/zenodo.1493841\n\n" +
            "* The nf-core framework\n" +
            "  https://doi.org/10.1038/s41587-020-0439-x\n\n" +
            "* Software dependencies\n" +
            "  https://github.com/${workflow.manifest.name}/blob/master/CITATIONS.md"
    }

    //
    // Print help to screen if required
    //
    public static String help(workflow, params, log) {
        def command = "nextflow run ${workflow.manifest.name} --input samplesheet.csv -profile docker"
        def help_string = ''
        help_string += NfcoreTemplate.logo(workflow, params.monochrome_logs)
        help_string += NfcoreSchema.paramsHelp(workflow, params, command)
        help_string += '\n' + citation(workflow) + '\n'
        help_string += NfcoreTemplate.dashedLine(params.monochrome_logs)
        return help_string
    }

    //
    // Print parameter summary log to screen
    //
    public static String paramsSummaryLog(workflow, params, log) {
        def summary_log = ''
        summary_log += NfcoreTemplate.logo(workflow, params.monochrome_logs)
        summary_log += NfcoreSchema.paramsSummaryLog(workflow, params)
        summary_log += '\n' + citation(workflow) + '\n'
        summary_log += NfcoreTemplate.dashedLine(params.monochrome_logs)
        return summary_log
    }

    //
    // Validate parameters and print summary to screen
    //
    public static void initialise(workflow, params, log) {
        // Print help to screen if required
        if (params.help) {
            log.info help(workflow, params, log)
            System.exit(0)
        }

        // Check that keys for reference databases are valid
        if (params.dada_ref_taxonomy && !params.skip_taxonomy) {
            dadareftaxonomyExistsError(params, log)
        }
        if (params.qiime_ref_taxonomy && !params.skip_taxonomy && !params.classifier) {
            qiimereftaxonomyExistsError(params, log)
        }

        // Validate workflow parameters via the JSON schema
        if (params.validate_params) {
            NfcoreSchema.validateParameters(workflow, params, log)
        }

        // Print parameter summary log to screen
        log.info paramsSummaryLog(workflow, params, log)

        // Check that a -profile or Nextflow config has been provided to run the pipeline
        NfcoreTemplate.checkConfigProvided(workflow, log)

        // Check that conda channels are set-up correctly
        if (params.enable_conda) {
            Utils.checkCondaChannels(log)
        }

        // Check AWS batch settings
        NfcoreTemplate.awsBatch(workflow, params)
    }

    //
    // Exit pipeline if incorrect --dada_ref_taxonomy key provided
    //
    private static void dadareftaxonomyExistsError(params, log) {
        if (params.dada_ref_databases && params.dada_ref_taxonomy && !params.dada_ref_databases.containsKey(params.dada_ref_taxonomy)) {
            log.error "=============================================================================\n" +
                "  DADA2 reference database '${params.dada_ref_taxonomy}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available reference taxonomy keys for `--dada_ref_taxonomy` are:\n" +
                "  ${params.dada_ref_databases.keySet().join(", ")}\n" +
                "==================================================================================="
            System.exit(1)
        }
    }

    //
    // Exit pipeline if incorrect --qiime_ref_taxonomy key provided
    //
    private static void qiimereftaxonomyExistsError(params, log) {
        if (params.qiime_ref_databases && params.qiime_ref_taxonomy && !params.qiime_ref_databases.containsKey(params.qiime_ref_taxonomy)) {
            log.error "=============================================================================\n" +
                "  QIIME2 reference database '${params.qiime_ref_taxonomy}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available reference taxonomy keys for `--qiime_ref_taxonomy` are:\n" +
                "  ${params.qiime_ref_databases.keySet().join(", ")}\n" +
                "==================================================================================="
            System.exit(1)
        }
    }
}
