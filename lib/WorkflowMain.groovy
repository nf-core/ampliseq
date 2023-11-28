//
// This file holds several functions specific to the main.nf workflow in the nf-core/ampliseq pipeline
//

import nextflow.Nextflow

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
    // Validate parameters and print summary to screen
    //
    public static void initialise(workflow, params, log) {

        // Check that keys for reference databases are valid
        if (params.dada_ref_taxonomy && !params.skip_taxonomy && !params.skip_dada_taxonomy) {
            dadareftaxonomyExistsError(params, log)
        }
        if (params.sintax_ref_taxonomy && !params.skip_taxonomy) {
            sintaxreftaxonomyExistsError(params, log)
        }
        if ((params.qiime_ref_taxonomy || params.qiime_ref_tax_custom) && !params.skip_taxonomy && !params.classifier) {
            qiimereftaxonomyExistsError(params, log)
        }

        // Print workflow version and exit on --version
        if (params.version) {
            String workflow_version = NfcoreTemplate.version(workflow)
            log.info "${workflow.manifest.name} ${workflow_version}"
            System.exit(0)
        }

        // Check that a -profile or Nextflow config has been provided to run the pipeline
        NfcoreTemplate.checkConfigProvided(workflow, log)

        // Check that conda channels are set-up correctly
        if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
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
            def error_string = "=============================================================================\n" +
                "  DADA2 reference database '${params.dada_ref_taxonomy}' not found in any config file provided to the pipeline.\n" +
                "  Currently, the available reference taxonomy keys for `--dada_ref_taxonomy` are:\n" +
                "  ${params.dada_ref_databases.keySet().join(", ")}\n" +
                "==================================================================================="
            Nextflow.error(error_string)
        }
    }
    //
    // Exit pipeline if incorrect --sintax_ref_taxonomy key provided
    //
    private static void sintaxreftaxonomyExistsError(params, log) {
        if (params.sintax_ref_databases && params.sintax_ref_taxonomy && !params.sintax_ref_databases.containsKey(params.sintax_ref_taxonomy)) {
            def error_string = "=============================================================================\n" +
                "  SINTAX reference database '${params.sintax_ref_taxonomy}' not found in any config file provided to the pipeline.\n" +
                "  Currently, the available reference taxonomy keys for `--sintax_ref_taxonomy` are:\n" +
                "  ${params.sintax_ref_databases.keySet().join(", ")}\n" +
                "==================================================================================="
            Nextflow.error(error_string)
        }
    }
    //
    // Exit pipeline if incorrect --qiime_ref_taxonomy key provided
    //
    private static void qiimereftaxonomyExistsError(params, log) {
        if (params.qiime_ref_databases && params.qiime_ref_taxonomy && !params.qiime_ref_databases.containsKey(params.qiime_ref_taxonomy)) {
            def error_string = "=============================================================================\n" +
                "  QIIME2 reference database '${params.qiime_ref_taxonomy}' not found in any config file provided to the pipeline.\n" +
                "  Currently, the available reference taxonomy keys for `--qiime_ref_taxonomy` are:\n" +
                "  ${params.qiime_ref_databases.keySet().join(", ")}\n" +
                "==================================================================================="
            Nextflow.error(error_string)
        }
    }
}
