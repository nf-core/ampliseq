// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided as a string i.e. "options.args"
//               where "params.options" is a Groovy Map that MUST be provided via the addParams section of the including workflow.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.

params.options = [:]
options        = initOptions(params.options)

process SBDIEXPORT {
    tag '${asvtable},${taxtable}'
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"SBDI", meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::r-tidyverse=1.2.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/r-tidyverse:1.2.1"
    } else {
        container "quay.io/biocontainers/r-tidyverse:1.2.1"
    }

    input:
    path asvtable
    path taxonomytable

    output:
    path "*.tsv", emit: sbditables

    script:
    def software = getSoftwareName(task.process)
    
    """
    sbdiexport.R $options.args 
    """
}
