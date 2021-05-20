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

nextflow.enable.dsl = 2

log.info Headers.nf_core(workflow, params.monochrome_logs)

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////

def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run nf-core/ampliseq --input 'samplesheet.tsv' -profile docker"
    log.info NfcoreSchema.params_help(workflow, params, json_schema, command)
    exit 0
}

////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////

def summary_params = NfcoreSchema.params_summary_map(workflow, params, json_schema)
log.info NfcoreSchema.params_summary_log(workflow, params, json_schema)

////////////////////////////////////////////////////
/* --         VALIDATE PARAMETERS              -- */
////////////////////////////////////////////////////
if (params.validate_params) {
    NfcoreSchema.validateParameters(params, json_schema, log)
}

////////////////////////////////////////////////////
/* --          PARAMETER CHECKS                -- */
////////////////////////////////////////////////////

// Check that conda channels are set-up correctly
if (params.enable_conda) {
    Checks.check_conda_channels(log)
}

// Check AWS batch settings
Checks.aws_batch(workflow, params)

// Check the hostnames against configured profiles
Checks.hostname(workflow, params, log)

///////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow {
        /*
         * SUBWORKFLOW: Run main nf-core/ampliseq analysis pipeline
         */
        include { AMPLISEQ } from './workflows/ampliseq' addParams( summary_params: summary_params )
        AMPLISEQ ()
}

////////////////////////////////////////////////////
/* --          CHECK PARAMETER ON ERROR        -- */
////////////////////////////////////////////////////

workflow.onError {
    // Print unexpected parameters - easiest is to just rerun validation
    NfcoreSchema.validateParameters(params, json_schema, log)
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
