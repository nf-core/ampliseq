#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/ampliseq
========================================================================================
    Github : https://github.com/nf-core/ampliseq
    Website: https://nf-co.re/ampliseq
    Slack  : https://nfcore.slack.com/channels/ampliseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { AMPLISEQ } from './workflows/ampliseq'

//
// WORKFLOW: Run main nf-core/ampliseq analysis pipeline
//
workflow NFCORE_AMPLISEQ {
    AMPLISEQ ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_AMPLISEQ ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
