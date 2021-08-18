// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SBDIEXPORTREANNOTATE {
    tag "${taxonomytable},${metadata}"
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
    path taxonomytable

    output:
    path "*.tsv", emit: sbdiannottables

    script:
    def software = getSoftwareName(task.process)

    """
    sbdiexportreannotate.R ${params.dada_ref_databases[params.dada_ref_taxonomy]["dbversion"]} $taxonomytable
    """
}
