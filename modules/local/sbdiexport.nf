// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SBDIEXPORT {
    tag "${asvtable},${taxonomytable},${metadata}"
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
    path metadata

    output:
    path "*.tsv", emit: sbditables

    script:
    def software = getSoftwareName(task.process)
    
    """
    sbdiexport.R $options.args $asvtable $taxonomytable $metadata
    """
}
